#!/usr/bin/env python3
"""Tune Samtools on A4; use collate -n 66 and omit stats (Python >= 3.6)."""

import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess
import time


GIB = 1024 ** 3
CANDIDATES = [8, 24, 48, 96, 127]
COLLATE_BINS = 66
# Available RAM gates, not Samtools allocation budgets. D1/D2 observed
# sort RSS peaks were 64.7/148.9 GiB; retain headroom for new thread counts.
SORT_AVAILABLE_GIB = {"D1": 96, "D2": 176, "D3": 64}
RESOURCE_WAIT_SECONDS = 30 * 60
APPS = ["flagstat", "filter", "collate", "fixmate", "sort", "markdup",
        "bam2sam", "sam2bam"]
OLD_THREADS = dict(flagstat=8, bam2sam=96, sam2bam=24,
                   filter=96, collate=96, fixmate=96, sort=24, markdup=127)
DATASETS = {
    "D1": {"label": "D1-WES-9G", "bytes": 9891133304,
           "records": 177798571, "sort_mib": 188736},
    "D2": {"label": "D2-WGS-25G", "bytes": 35255807945,
           "records": 365079492, "sort_mib": 196608},
    "D3": {"label": "D3-WGS-85G", "bytes": 89994058306,
           "records": 1210633370, "sort_mib": 49152},
}
FIELDS = ["phase", "dataset", "app", "threads", "repeat", "wall_s",
          "output_bytes", "records_checked", "sort_memory_mib",
          "maxrss_kb", "user_s", "sys_s", "cpu", "fs_inputs",
          "fs_outputs", "exit", "log"]
TIME_FORMAT = ("wall_s=%e user_s=%U sys_s=%S cpu=%P maxrss_kb=%M "
               "fs_inputs=%I fs_outputs=%O exit=%x")


def sort_memory(dataset, threads):
    return DATASETS[dataset]["sort_mib"] // max(1, threads)


def required_memory(dataset, app):
    return (SORT_AVAILABLE_GIB[dataset] if app in ("sort", "workflow") else 16) * GIB


def candidates(dataset, app):
    if app == "flagstat":
        return [0, 4, 8, 16, 24, 48, 96, 127]
    # The previous sort search still applies; collate changed its bin count.
    if dataset == "D3" and app == "sort":
        return [24, 96]
    return CANDIDATES[:]


def inherited_thread(dataset, app):
    return 96 if dataset == "D3" and app == "sort" else OLD_THREADS[app]


def shortlist(rows, inherited):
    ranked = sorted(rows, key=lambda row: (float(row["wall_s"]),
                                           int(row["threads"])))
    return sorted(set([int(row["threads"]) for row in ranked[:2]] + [inherited]))


def parse_time(path):
    values = dict(re.findall(r"(\w+)=([^\s]+)", path.read_text()))
    if "wall_s" not in values or values.get("exit") != "0":
        raise RuntimeError("missing successful GNU time result: " + str(path))
    return values


def memory_available():
    with open("/proc/meminfo") as source:
        for line in source:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    raise RuntimeError("MemAvailable is unavailable")


class Benchmark:
    def __init__(self, args):
        self.args = args
        self.root = Path(args.root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock = open(str(self.root.parent / "retune.lock"), "a")
        fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (self.root / "results.tsv").exists():
            raise RuntimeError("refusing to overwrite an existing benchmark")
        for name in ("logs", "work", "tmp", "prepared"):
            (self.root / name).mkdir()
        self.env = os.environ.copy()
        self.env["LD_LIBRARY_PATH"] = args.libdeflate + ":" + self.env.get(
            "LD_LIBRARY_PATH", "")
        self.env["TMPDIR"] = str(self.root / "tmp")
        (self.root / "driver.pid").write_text(str(os.getpid()) + "\n")
        self.rows = []
        self.cached_rows = {}
        self.selected = {}
        self.filter_count = {}
        self.writer_file = open(str(self.root / "results.tsv"), "w", newline="")
        self.writer = csv.DictWriter(self.writer_file, FIELDS, delimiter="\t")
        self.writer.writeheader()
        self.writer_file.flush()
        self.record_environment()
        if getattr(args, "resume_from", None):
            self.import_completed(Path(args.resume_from).resolve())
        plan = {"datasets": DATASETS, "apps": APPS,
                "D1_scope": "collate tuning + three ordinary workflows only",
                "collate_bins": COLLATE_BINS,
                "coarse_candidates": CANDIDATES,
                "flagstat_candidates": candidates("D2", "flagstat"),
                "D3_sort_confirm_candidates": [24, 96],
                "confirmation": "top two coarse + inherited; 3 interleaved rounds",
                "sort_budget": "per-thread -m = floor(total MiB / max(1, -@))",
                "resource_guard": {
                    "sort_available_gib": SORT_AVAILABLE_GIB,
                    "other_available_gib": 16,
                    "wait_seconds": getattr(args, "resource_wait_seconds", RESOURCE_WAIT_SECONDS),
                    "basis": "observed sort RSS plus headroom; Samtools -m unchanged"},
                "timing": "GNU time wall; real bigssd I/O included; validation excluded",
                "cache": "not dropped; coarse first round reverses by application",
                "workflows": args.workflows,
                "resume_from": getattr(args, "resume_from", None),
                "samtools": args.samtools}
        (self.root / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")

    @staticmethod
    def sample_key(dataset, app, threads, phase, repeat):
        return (dataset, app, str(threads), phase, str(repeat))

    def import_completed(self, source):
        plan = json.loads((source / "plan.json").read_text())
        if plan["samtools"] != self.args.samtools:
            raise RuntimeError("resume Samtools executable differs")
        same_bins = plan.get("collate_bins") == COLLATE_BINS
        datasets = ("D1", "D2", "D3") if same_bins else ("D2", "D3")
        for dataset in datasets:
            if plan["datasets"][dataset] != DATASETS[dataset]:
                raise RuntimeError("resume dataset metadata differs: " + dataset)
        with open(str(source / "results.tsv")) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FIELDS:
                raise RuntimeError("resume result schema differs")
            for row in reader:
                # Reuse collate only when its bin count matches. Workflows and
                # prepared inputs are regenerated in the new independent run.
                if (row["dataset"] not in datasets or
                        row["app"] not in APPS or
                        (row["app"] == "collate" and not same_bins) or
                        (row["dataset"] == "D1" and row["app"] != "collate") or
                        row["phase"] not in ("coarse", "confirm")):
                    continue
                if (row["exit"] != "0" or float(row["wall_s"]) < 0 or
                        not Path(row["log"]).is_file()):
                    raise RuntimeError("incomplete resume sample: " + str(row))
                key = self.sample_key(row["dataset"], row["app"], row["threads"],
                                      row["phase"], row["repeat"])
                if key in self.cached_rows:
                    raise RuntimeError("duplicate resume sample: " + str(key))
                if row["app"] == "filter" and row["records_checked"]:
                    count = int(row["records_checked"])
                    previous = self.filter_count.setdefault(row["dataset"], count)
                    if previous != count:
                        raise RuntimeError("resume filter counts differ")
                self.cached_rows[key] = row
                self.rows.append(row)
                self.writer.writerow(row)
        self.writer_file.flush()
        (self.root / "resume.json").write_text(json.dumps(
            dict(source=str(source), imported_samples=len(self.cached_rows),
                 collate_bins_match=same_bins,
                 excluded_apps=["stats", "workflow"] + ([] if same_bins else ["collate"]),
                 prepared_inputs="not reused; regenerate with collate -n 66"),
            indent=2) + "\n")
        self.log("IMPORTED {} unchanged samples from {}".format(
            len(self.cached_rows), source))

    def log(self, message):
        print("[{}] {}".format(time.strftime("%Y-%m-%d %H:%M:%S%z"), message),
              flush=True)
        (self.root / "state").write_text(message + "\n")

    def record_environment(self):
        with open(str(self.root / "environment.txt"), "w") as output:
            for command in (["date", "-Is"], ["hostname"],
                            [self.args.samtools, "--version"], ["lscpu"],
                            ["free", "-h"], ["df", "-hT", str(self.root)],
                            ["sha256sum", self.args.samtools, __file__]):
                output.write("$ " + " ".join(map(shlex.quote, command)) + "\n")
                output.flush()
                subprocess.run(command, stdout=output, stderr=output,
                               env=self.env, check=True)

    def input(self, dataset):
        label = DATASETS[dataset]["label"]
        return Path(self.args.data) / label / "initial" / (label + ".coordinate.bam")

    def record_command(self, label, command):
        with open(str(self.root / "commands.jsonl"), "a") as output:
            output.write(json.dumps({"label": label, "argv": list(map(str, command)),
                                     "time": time.time()}) + "\n")
        with open(str(self.root / "commands.sh"), "a") as output:
            output.write("# " + label + "\n")
            output.write(" ".join(shlex.quote(str(arg)) for arg in command) + "\n")

    def snapshot(self, label):
        with open(str(self.root / "system.jsonl"), "a") as output:
            output.write(json.dumps({"label": label, "time": time.time(),
                                     "loadavg": os.getloadavg(),
                                     "available_bytes": memory_available(),
                                     "disk_free_bytes": shutil.disk_usage(str(self.root)).free})
                         + "\n")
        with open(str(self.root / "processes.log"), "a") as output:
            output.write("\n" + label + "\n")
            output.flush()
            subprocess.run(["ps", "-eo", "pid,user,pcpu,pmem,etime,args", "--sort=-pcpu"],
                           stdout=output, stderr=output, check=True)

    def wait_resources(self, disk_bytes=100 * GIB, memory_bytes=16 * GIB):
        deadline = time.monotonic() + getattr(
            self.args, "resource_wait_seconds", RESOURCE_WAIT_SECONDS)
        while True:
            free = shutil.disk_usage(str(self.root)).free
            available = memory_available()
            load = os.getloadavg()[0]
            if free < disk_bytes:
                raise RuntimeError("insufficient free disk: {} < {}".format(free, disk_bytes))
            if available >= memory_bytes and load < 16:
                return
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise RuntimeError("resource wait timed out: available={} required={} load={:.2f}".format(
                    available, memory_bytes, load))
            self.log("WAIT_RESOURCES free={} available={} required={} load={:.2f} remaining_s={:.0f}".format(
                free, available, memory_bytes, load, remaining))
            time.sleep(min(60, remaining))

    def command(self, dataset, app, threads, out, prefix, input_path=None):
        s = self.args.samtools
        src = str(input_path or self.input(dataset))
        prep = self.root / "prepared" / dataset
        if app == "flagstat":
            return [s, app, "-@", str(threads), src]
        if app == "bam2sam":
            return [s, "view", "--no-PG", "-@", str(threads), "-h", "-O", "SAM",
                    "-o", str(out), src]
        if app == "sam2bam":
            return [s, "view", "--no-PG", "-@", str(threads), "-O", "BAM,level=1",
                    "-o", str(out), str(input_path or prep / "input.sam")]
        if app == "filter":
            return [s, "view", "--no-PG", "-@", str(threads), "-q", "30",
                    "-O", "BAM,level=1", "-o", str(out), src]
        if app == "collate":
            return [s, "collate", "--no-PG", "-@", str(threads), "-n", str(COLLATE_BINS), "-l", "1",
                    "-T", str(prefix), "-o", str(out), src]
        if app == "fixmate":
            return [s, "fixmate", "--no-PG", "-@", str(threads), "-m",
                    "--output-fmt", "BAM,level=1", str(input_path or prep / "name.bam"), str(out)]
        if app == "sort":
            return [s, "sort", "--no-PG", "-@", str(threads), "-m",
                    str(sort_memory(dataset, threads)) + "M", "-l", "1", "-T", str(prefix),
                    "-o", str(out), str(input_path or prep / "fixmate.bam")]
        if app == "markdup":
            return [s, "markdup", "--no-PG", "-@", str(threads), "-l", "300",
                    "-T", str(prefix), "--output-fmt", "BAM,level=1",
                    str(input_path or prep / "fixmate.sorted.bam"), str(out)]
        raise ValueError(app)

    def measured(self, label, command, stdout_path):
        time_log = self.root / "logs" / (label + ".time")
        err_log = self.root / "logs" / (label + ".stderr")
        self.record_command(label, command)
        self.snapshot(label + ".start")
        self.log("RUNNING " + label)
        with open(str(stdout_path), "w") as out, open(str(err_log), "w") as err:
            subprocess.run(["/usr/bin/time", "-o", str(time_log), "-f", TIME_FORMAT]
                           + command, stdout=out, stderr=err, env=self.env, check=True)
        self.snapshot(label + ".end")
        return parse_time(time_log), time_log, err_log

    def validate(self, path, full_count=False, expected=None):
        command = [self.args.samtools, "quickcheck", "-v", str(path)]
        self.record_command("quickcheck", command)
        subprocess.run(command, env=self.env, check=True)
        if not full_count:
            return ""
        command = [self.args.samtools, "view", "-@", "48", "-c", str(path)]
        self.record_command("record-count (outside measured time)", command)
        count = int(subprocess.check_output(command, env=self.env, universal_newlines=True))
        if expected is not None and count != expected:
            raise RuntimeError("record count {} != {} in {}".format(count, expected, path))
        return count

    def check_sort(self, dataset, stderr):
        text = stderr.read_text()
        match = re.search(r"merging from (\d+) files", text)
        files = int(match[1]) if match else 0
        if (dataset in ("D1", "D2") and files != 0) or (dataset == "D3" and files == 0):
            raise RuntimeError("unexpected sort mode in " + str(stderr))

    def prepare(self, dataset, app):
        prep = self.root / "prepared" / dataset
        prep.mkdir(exist_ok=True)
        needed = {"fixmate": ("name.bam", "collate"),
                  "sort": ("fixmate.bam", "fixmate"),
                  "markdup": ("fixmate.sorted.bam", "sort"),
                  "sam2bam": ("input.sam", "bam2sam")}
        if app not in needed:
            return
        name, producer = needed[app]
        out = prep / name
        if out.exists():
            return
        self.prepare(dataset, producer)
        self.wait_resources(disk_bytes=7 * DATASETS[dataset]["bytes"] + 100 * GIB,
                            memory_bytes=required_memory(dataset, producer))
        label = "prepare.{}.{}".format(dataset, producer)
        _, _, err = self.measured(label, self.command(dataset, producer, 48, out,
                                                   self.root / "tmp" / label),
                                  self.root / "logs" / (label + ".stdout"))
        if producer == "sort":
            self.check_sort(dataset, err)
        if out.suffix == ".bam":
            self.validate(out)
        self.log("PREPARED " + str(out))

    def run_app(self, dataset, app, threads, phase, repeat):
        key = self.sample_key(dataset, app, threads, phase, repeat)
        if key in self.cached_rows:
            self.log("REUSE {}".format(".".join(key)))
            return self.cached_rows[key]
        self.prepare(dataset, app)
        label = "{}.{}.{}.t{}.r{}".format(phase, dataset, app, threads, repeat)
        suffix = ".sam" if app == "bam2sam" else ".bam"
        if app == "flagstat":
            out = self.root / "logs" / (label + ".stdout")
        else:
            out = self.root / "work" / (label + suffix)
        memory = required_memory(dataset, app)
        disk = (7 * DATASETS[dataset]["bytes"] if app == "bam2sam" else
                3 * DATASETS[dataset]["bytes"]) + 100 * GIB
        self.wait_resources(disk_bytes=disk, memory_bytes=memory)
        timing, log, err = self.measured(
            label, self.command(dataset, app, threads, out, self.root / "tmp" / label),
            out if app == "flagstat" else self.root / "logs" / (label + ".stdout"))
        checked = ""
        if app == "sort":
            self.check_sort(dataset, err)
        if app not in ("flagstat", "bam2sam"):
            expected = None if app == "filter" else DATASETS[dataset]["records"]
            checked = self.validate(out, phase == "confirm" and repeat == 1, expected)
            if app == "filter" and checked != "":
                if dataset in self.filter_count and self.filter_count[dataset] != checked:
                    raise RuntimeError("filtered record count differs between configurations")
                self.filter_count[dataset] = checked
        elif app == "bam2sam":
            subprocess.run([self.args.samtools, "view", "-H", str(out)],
                           env=self.env, stdout=subprocess.DEVNULL, check=True)
        elif app == "flagstat":
            match = re.search(r"^(\d+) \+ (\d+) in total", out.read_text(), re.M)
            if not match or sum(map(int, match.groups())) != DATASETS[dataset]["records"]:
                raise RuntimeError("flagstat total mismatch")
            checked = DATASETS[dataset]["records"]
        row = dict(timing, phase=phase, dataset=dataset, app=app, threads=threads,
                   repeat=repeat, output_bytes=out.stat().st_size, records_checked=checked,
                   sort_memory_mib=sort_memory(dataset, threads) if app == "sort" else "",
                   log=str(log))
        self.writer.writerow(row)
        self.writer_file.flush()
        self.rows.append(row)
        self.cached_rows[key] = row
        if app != "flagstat":
            if app == "bam2sam" and not (self.root / "prepared" / dataset / "input.sam").exists():
                out.rename(self.root / "prepared" / dataset / "input.sam")
            else:
                out.unlink()
        self.log("DONE {} wall_s={}".format(label, timing["wall_s"]))
        return row

    def tune(self, dataset, app, index):
        coarse = []
        search = candidates(dataset, app)
        if not (dataset == "D3" and app == "sort"):
            for threads in search[::(-1 if index % 2 else 1)]:
                coarse.append(self.run_app(dataset, app, threads, "coarse", 1))
            search = shortlist(coarse, inherited_thread(dataset, app))
        rows = []
        for repeat in range(1, 4):
            order = search if repeat != 2 else list(reversed(search))
            for threads in order:
                rows.append(self.run_app(dataset, app, threads, "confirm", repeat))
        summary = []
        for threads in search:
            times = [float(row["wall_s"]) for row in rows if int(row["threads"]) == threads]
            summary.append(dict(dataset=dataset, app=app, threads=threads,
                                times_s=times, median_s=statistics.median(times)))
        selected = min(summary, key=lambda row: (row["median_s"], row["threads"]))
        self.selected[(dataset, app)] = selected["threads"]
        with open(str(self.root / "summary.jsonl"), "a") as output:
            output.write(json.dumps(dict(configurations=summary, selected=selected)) + "\n")
        self.log("SELECTED {} {} threads={} median_s={}".format(
            dataset, app, selected["threads"], selected["median_s"]))

    def workflow(self, dataset, repeat):
        label = "workflow.{}.r{}".format(dataset, repeat)
        work = self.root / "work" / label
        work.mkdir()
        paths = [self.input(dataset)] + [work / name for name in
                 ("name.bam", "fixmate.bam", "sorted.bam", "markdup.bam")]
        self.wait_resources(disk_bytes=6 * DATASETS[dataset]["bytes"] + 100 * GIB,
                            memory_bytes=required_memory(dataset, "workflow"))
        self.snapshot(label + ".start")
        self.log("RUNNING " + label)
        start = time.monotonic()
        stage_times = []
        for index, app in enumerate(("collate", "fixmate", "sort", "markdup")):
            threads = self.selected[(dataset, app)]
            command = self.command(dataset, app, threads, paths[index + 1],
                                   self.root / "tmp" / (label + "." + app),
                                   input_path=paths[index])
            stage_label = label + "." + app
            timing, _, err = self.measured(stage_label, command,
                                          self.root / "logs" / (stage_label + ".stdout"))
            if app == "sort":
                self.check_sort(dataset, err)
            stage_times.append(float(timing["wall_s"]))
            if index > 0:
                paths[index].unlink()
        wall = time.monotonic() - start
        count = self.validate(paths[-1], True, DATASETS[dataset]["records"])
        row = dict(phase="workflow", dataset=dataset, app="workflow",
                   threads="/".join(str(self.selected[(dataset, app)]) for app in
                                    ("collate", "fixmate", "sort", "markdup")),
                   repeat=repeat, wall_s="{:.6f}".format(wall),
                   output_bytes=paths[-1].stat().st_size, records_checked=count,
                   log=str(self.root / "logs" / (label + ".json")))
        Path(row["log"]).write_text(json.dumps(dict(row, stage_times_s=stage_times), indent=2) + "\n")
        self.writer.writerow(row)
        self.writer_file.flush()
        paths[-1].unlink()
        work.rmdir()
        self.log("DONE {} wall_s={:.6f}".format(label, wall))

    def run(self):
        for dataset in ("D1", "D2", "D3"):
            src = self.input(dataset)
            if src.stat().st_size != DATASETS[dataset]["bytes"]:
                raise RuntimeError("unexpected input file size: " + str(src))
            self.validate(src)
            if dataset == "D1":
                self.tune(dataset, "collate", 0)
                for app in ("fixmate", "sort", "markdup"):
                    self.selected[(dataset, app)] = OLD_THREADS[app]
                if self.args.workflows:
                    for repeat in range(1, 4):
                        self.workflow(dataset, repeat)
                self.log("DATASET_DONE " + dataset)
                continue
            for index, app in enumerate(APPS):
                if app == "bam2sam":
                    if self.args.workflows:
                        for repeat in range(1, 4):
                            self.workflow(dataset, repeat)
                    # Conversion needs a large SAM; release only our own prepared BAMs.
                    for name in ("name.bam", "fixmate.bam", "fixmate.sorted.bam"):
                        path = self.root / "prepared" / dataset / name
                        if path.exists():
                            path.unlink()
                self.tune(dataset, app, index)
            sam = self.root / "prepared" / dataset / "input.sam"
            if sam.exists():
                sam.unlink()
            self.log("DATASET_DONE " + dataset)
        self.log("COMPLETE")
        (self.root / "ALL_DONE").touch()


def self_test():
    assert "stats" not in APPS
    assert COLLATE_BINS == 66
    assert [required_memory(dataset, "sort") // GIB for dataset in ("D1", "D2", "D3")] == [96, 176, 64]
    assert required_memory("D2", "workflow") == 176 * GIB
    assert required_memory("D2", "collate") == 16 * GIB
    assert sort_memory("D1", 24) == 7864
    assert sort_memory("D2", 24) == 8192
    assert sort_memory("D3", 24) == 2048
    assert sort_memory("D3", 48) == 1024
    assert sort_memory("D3", 96) == 512
    assert sort_memory("D3", 127) * 127 <= 49152
    assert shortlist([dict(threads=8, wall_s="4"), dict(threads=24, wall_s="3"),
                      dict(threads=96, wall_s="2")], 127) == [24, 96, 127]
    assert candidates("D3", "sort") == [24, 96]
    assert candidates("D3", "collate") == CANDIDATES
    assert inherited_thread("D3", "sort") == 96
    bench = Benchmark.__new__(Benchmark)
    bench.args = argparse.Namespace(samtools="samtools", data="/data")
    bench.root = Path("/bench")
    for dataset in ("D1", "D2", "D3"):
        command = bench.command(dataset, "collate", 96, "out.bam", "tmp")
        assert command[command.index("-n") + 1] == "66"
    command = bench.command("D1", "sort", 24, "out.bam", "tmp")
    assert command[command.index("-m") + 1] == "7864M"
    try:
        bench.command("D2", "stats", 24, "out.txt", "tmp")
    except ValueError:
        pass
    else:
        raise AssertionError("stats must not be a benchmark application")
    print("self-test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root")
    parser.add_argument("--data", default="/home/bigssd/wzs_data")
    parser.add_argument("--samtools", default="/home/user_home/wzs/software/samtools-1.20/bin/samtools")
    parser.add_argument("--libdeflate", default="/home/user_home/wzs/software/libdeflate-1.20/lib64")
    parser.add_argument("--workflows", action="store_true")
    parser.add_argument("--resume-from", help="reuse compatible completed samples in a new run")
    parser.add_argument("--resource-wait-seconds", type=int, default=RESOURCE_WAIT_SECONDS,
                        help="maximum wait for RAM/load before stopping (default: 1800)")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.resource_wait_seconds < 0:
        parser.error("--resource-wait-seconds must be nonnegative")
    if args.self_test:
        self_test()
        return
    if not args.root:
        parser.error("--root is required")
    bench = None
    try:
        bench = Benchmark(args)
        bench.run()
    except Exception as error:
        print("FAILED: " + str(error), flush=True)
        if bench:
            (bench.root / "FAILED").write_text(str(error) + "\n")
            bench.log("FAILED " + str(error))
        raise


if __name__ == "__main__":
    main()
