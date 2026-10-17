#!/usr/bin/env python3
"""D2 external sort, collate -n 66, and ordinary four-command workflows."""

import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import statistics
import subprocess
import time

from bench_samtools_a4_retune import Benchmark, GIB, TIME_FORMAT, parse_time


STAGES = ("collate", "fixmate", "sort", "markdup")
FIELDS = ["phase", "dataset", "app", "threads", "repeat", "wall_s",
          "output_bytes", "records_checked", "duplicates_checked",
          "sort_memory_mib", "sort_runs", "temp_files_peak", "temp_bytes_peak",
          "maxrss_kb", "user_s", "sys_s", "cpu", "fs_inputs", "fs_outputs",
          "exit", "log"]


def sort_memory(total_mib, threads):
    return total_mib // max(1, threads)


def sort_runs(text):
    match = re.search(r"merging from (\d+) files", text)
    if not match or int(match.group(1)) <= 0:
        raise RuntimeError("sort did not generate external runs")
    return int(match.group(1))


class ExternalBenchmark(Benchmark):
    def __init__(self, args):
        self.args = args
        self.root = Path(args.root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock = open(str(self.root.parent / "retune.lock"), "a")
        fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (self.root / "results.tsv").exists():
            raise RuntimeError("refusing to overwrite an existing experiment")
        for name in ("logs", "work", "tmp", "prepared"):
            (self.root / name).mkdir()
        self.env = os.environ.copy()
        self.env["LD_LIBRARY_PATH"] = args.libdeflate + ":" + self.env.get("LD_LIBRARY_PATH", "")
        self.env["TMPDIR"] = str(self.root / "tmp")
        self.rows = []
        self.selected = dict(collate=127, fixmate=127, markdup=24)
        (self.root / "driver.pid").write_text(str(os.getpid()) + "\n")
        self.writer_file = open(str(self.root / "results.tsv"), "w", newline="")
        self.writer = csv.DictWriter(self.writer_file, FIELDS, delimiter="\t")
        self.writer.writeheader()
        self.writer_file.flush()
        plan = dict(dataset="D2", input=str(args.input), records=args.records,
                    duplicates=args.duplicates, input_bytes=args.input_bytes,
                    collate_bins=66, compression_level=1,
                    sort_budget_mib=args.sort_mib, candidates=args.candidates,
                    sort_budget="-m=floor(total MiB/max(1,-@)); approximate sort arena, not RSS",
                    collate_budget="no -m option; standard temporary-file implementation",
                    repeats=args.repeats, inherited_threads=self.selected,
                    timing="GNU time wall includes NVMe input, runs, and final output; validation excluded",
                    workflow_timing="sum the four measured command walls per round; median across rounds",
                    cache="not dropped; preparation before formal samples",
                    temp_sampling_seconds=args.temp_sample_interval,
                    samtools=args.samtools)
        (self.root / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
        self.record_environment()
        with open(str(self.root / "environment.txt"), "a") as output:
            subprocess.run(["sha256sum", __file__, args.input],
                           stdout=output, stderr=output, check=True)

    def input(self, dataset="D2"):
        return Path(self.args.input)

    def command(self, dataset, app, threads, out, prefix, input_path=None):
        s = self.args.samtools
        src = str(input_path or self.input())
        common = [s, app, "--no-PG", "-@", str(threads)]
        if app == "collate":
            return common + ["-n", "66", "-l", "1", "-T", str(prefix), "-o", str(out), src]
        if app == "fixmate":
            return common + ["-m", "--output-fmt", "BAM,level=1", src, str(out)]
        if app == "sort":
            return common + ["-m", str(sort_memory(self.args.sort_mib, threads)) + "M",
                             "-l", "1", "-T", str(prefix), "-o", str(out), src]
        if app == "markdup":
            return common + ["-l", "300", "-T", str(prefix),
                             "--output-fmt", "BAM,level=1", src, str(out)]
        raise ValueError(app)

    def measured(self, label, command, stdout_path, temporary):
        time_log = self.root / "logs" / (label + ".time")
        err_log = self.root / "logs" / (label + ".stderr")
        self.record_command(label, command)
        self.snapshot(label + ".start")
        self.log("RUNNING " + label)
        peak_files = peak_bytes = 0
        with open(str(stdout_path), "w") as out, open(str(err_log), "w") as err:
            process = subprocess.Popen(["/usr/bin/time", "-o", str(time_log),
                                        "-f", TIME_FORMAT] + command,
                                       stdout=out, stderr=err, env=self.env)
            while True:
                # Metadata sampling confirms collate really creates temp files.
                sizes = []
                for path in temporary.iterdir():
                    try:
                        if path.is_file():
                            sizes.append(path.stat().st_size)
                    except FileNotFoundError:
                        pass
                peak_files = max(peak_files, len(sizes))
                peak_bytes = max(peak_bytes, sum(sizes))
                try:
                    status = process.wait(timeout=self.args.temp_sample_interval)
                    break
                except subprocess.TimeoutExpired:
                    pass
            if status:
                raise RuntimeError("command failed: {} exit={}".format(label, status))
        self.snapshot(label + ".end")
        return parse_time(time_log), time_log, err_log, peak_files, peak_bytes

    def validate_final(self, path, label):
        output = self.root / "logs" / (label + ".flagstat")
        command = [self.args.samtools, "flagstat", "-@", "24", str(path)]
        self.record_command(label + " validation (outside timing)", command)
        with open(str(output), "w") as handle:
            subprocess.run(command, env=self.env, stdout=handle, check=True)
        text = output.read_text()
        counts = {}
        for key, suffix in (("records", "in total"), ("duplicates", "duplicates")):
            match = re.search(r"^(\d+) \+ (\d+) " + suffix + r"(?:\s|$)", text, re.M)
            if not match:
                raise RuntimeError("missing flagstat " + key)
            counts[key] = sum(map(int, match.groups()))
        if counts != dict(records=self.args.records, duplicates=self.args.duplicates):
            raise RuntimeError("final counts differ: " + repr(counts))
        stats = self.root / "logs" / (label + ".stats")
        command = [self.args.samtools, "stats", "-@", "24", str(path)]
        self.record_command(label + " sorted validation (outside timing)", command)
        with open(str(stats), "w") as handle:
            subprocess.run(command, env=self.env, stdout=handle, check=True)
        if not re.search(r"^SN\tis sorted:\t1(?:\s|$)", stats.read_text(), re.M):
            raise RuntimeError("final output is not coordinate sorted")
        return counts

    def sample(self, app, threads, phase, repeat, src, out, full=False):
        label = "{}.D2.{}.t{}.r{}".format(phase, app, threads, repeat)
        temporary = self.root / "tmp" / label
        temporary.mkdir()
        self.wait_resources(disk_bytes=400 * GIB,
                            memory_bytes=(64 if app == "sort" else 16) * GIB)
        command = self.command("D2", app, threads, out, temporary / "run", src)
        timing, log, err, files, size = self.measured(
            label, command, self.root / "logs" / (label + ".stdout"), temporary)
        runs = sort_runs(err.read_text()) if app == "sort" else ""
        if app == "collate" and not files:
            raise RuntimeError("collate did not create observed temporary files")
        checked = self.validate(out, full, self.args.records)
        duplicated = ""
        if app == "markdup":
            counts = self.validate_final(out, label)
            checked, duplicated = counts["records"], counts["duplicates"]
        row = dict(timing, phase=phase, dataset="D2", app=app, threads=threads,
                   repeat=repeat, output_bytes=out.stat().st_size,
                   records_checked=checked, duplicates_checked=duplicated,
                   sort_memory_mib=sort_memory(self.args.sort_mib, threads) if app == "sort" else "",
                   sort_runs=runs, temp_files_peak=files, temp_bytes_peak=size, log=str(log))
        self.writer.writerow(row)
        self.writer_file.flush()
        self.rows.append(row)
        self.log("DONE {} wall_s={} runs={}".format(label, timing["wall_s"], runs))
        if list(temporary.iterdir()):
            raise RuntimeError("unexpected remaining temporary files: " + str(temporary))
        temporary.rmdir()
        return row

    def tune_sort(self, src):
        coarse = []
        for threads in self.args.candidates:
            out = self.root / "work" / ("coarse.sort.t{}.bam".format(threads))
            coarse.append(self.sample("sort", threads, "coarse", 0, src, out))
            out.unlink()
        ranked = sorted(coarse, key=lambda row: float(row["wall_s"]))
        choices = sorted(set([int(row["threads"]) for row in ranked[:2]] + [24]))
        confirmed = []
        for repeat in range(1, self.args.repeats + 1):
            for threads in choices[::(-1 if repeat == 2 else 1)]:
                out = self.root / "work" / ("confirm.sort.t{}.r{}.bam".format(threads, repeat))
                confirmed.append(self.sample("sort", threads, "confirm", repeat, src, out,
                                             full=repeat == 1))
                out.unlink()
        summaries = []
        for threads in choices:
            values = [float(row["wall_s"]) for row in confirmed if row["threads"] == threads]
            summaries.append(dict(threads=threads, times_s=values,
                                  median_s=statistics.median(values)))
        selected = min(summaries, key=lambda row: row["median_s"])
        self.selected["sort"] = selected["threads"]
        (self.root / "sort_selection.json").write_text(json.dumps(
            dict(configurations=summaries, selected=selected), indent=2) + "\n")
        self.log("SELECTED sort threads={} median_s={}".format(
            selected["threads"], selected["median_s"]))

    def workflow(self, repeat):
        directory = self.root / "work" / ("workflow.r{}".format(repeat))
        directory.mkdir()
        paths = [self.input()] + [directory / (app + ".bam") for app in STAGES]
        times = []
        for index, app in enumerate(STAGES):
            row = self.sample(app, self.selected[app], "workflow", repeat,
                              paths[index], paths[index + 1])
            times.append(float(row["wall_s"]))
            if index:
                paths[index].unlink()
        with open(str(self.root / "workflows.tsv"), "a") as handle:
            handle.write("{}\t{}\t{:.6f}\n".format(repeat,
                         "\t".join("{:.6f}".format(value) for value in times), sum(times)))
        paths[-1].unlink()
        directory.rmdir()

    def run(self):
        if self.input().stat().st_size != self.args.input_bytes:
            raise RuntimeError("D2 input file size differs")
        self.validate(self.input())
        name = self.root / "prepared" / "name.bam"
        fixmate = self.root / "prepared" / "fixmate.bam"
        self.sample("collate", 127, "prepare", 0, self.input(), name, full=True)
        self.sample("fixmate", 127, "prepare", 0, name, fixmate, full=True)
        name.unlink()
        for repeat in range(1, self.args.repeats + 1):
            out = self.root / "work" / ("collate.r{}.bam".format(repeat))
            self.sample("collate", 127, "confirm", repeat, self.input(), out,
                        full=repeat == 1)
            out.unlink()
        self.tune_sort(fixmate)
        fixmate.unlink()
        (self.root / "workflows.tsv").write_text("repeat\tcollate_s\tfixmate_s\tsort_s\tmarkdup_s\ttotal_s\n")
        for repeat in range(1, self.args.repeats + 1):
            self.workflow(repeat)
        with open(str(self.root / "workflows.tsv")) as handle:
            workflows = list(csv.DictReader(handle, delimiter="\t"))
        collate = [float(row["wall_s"]) for row in self.rows
                   if row["app"] == "collate" and row["phase"] == "confirm"]
        summary = dict(collate_median_s=statistics.median(collate),
                       sort=self.selected["sort"],
                       workflow_times_s=[float(row["total_s"]) for row in workflows],
                       workflow_median_s=statistics.median(float(row["total_s"]) for row in workflows))
        (self.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        self.log("COMPLETE")
        (self.root / "ALL_DONE").touch()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--libdeflate", required=True)
    parser.add_argument("--sort-mib", type=int, default=49152)
    parser.add_argument("--candidates", type=int, nargs="+", default=[24, 48, 96, 127])
    parser.add_argument("--records", type=int, default=365079492)
    parser.add_argument("--duplicates", type=int, default=8894488)
    parser.add_argument("--input-bytes", type=int, default=35255807945)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--resource-wait-seconds", type=int, default=1800)
    parser.add_argument("--temp-sample-interval", type=float, default=1)
    args = parser.parse_args()
    if args.repeats <= 0 or args.temp_sample_interval <= 0 or any(value <= 0 for value in args.candidates) or \
            args.sort_mib < max(args.candidates + [24]):
        parser.error("positive repeats/threads and at least 1 MiB per sort thread required")
    bench = ExternalBenchmark(args)
    try:
        bench.run()
    except Exception as error:
        bench.log("FAILED " + str(error))
        raise
    finally:
        bench.writer_file.close()
        bench.lock.close()


if __name__ == "__main__":
    main()
