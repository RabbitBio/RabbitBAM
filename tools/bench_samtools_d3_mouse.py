#!/usr/bin/env python3
"""Independent mouse D3 library curves, eight-app tuning, and three workflows."""

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import time

import bench_samtools_a4_retune as base


DATASET = "D3-mouse-SRR7768905"
INPUT_BYTES = 14613289720
RECORDS = 165373770
KEPT = 134417603
DUPLICATES = 33171753
LIBRARY_THREADS = [0, 2, 4, 8, 16, 24, 32, 47, 48, 64, 96, 127]


def cpu_counters():
    fields = Path("/proc/stat").read_text().splitlines()[0].split()[1:9]
    values = list(map(int, fields))
    return sum(values), values[3] + values[4]


def current_cpu_busy():
    total, idle = cpu_counters()
    time.sleep(0.5)
    next_total, next_idle = cpu_counters()
    return 1.0 - float(next_idle - idle) / max(1, next_total - total)


class MouseBenchmark(base.Benchmark):
    def __init__(self, args):
        self.expected = dict(records=args.records, kept=args.kept, duplicates=args.duplicates)
        self.smoke = args.smoke
        base.DATASETS[DATASET] = dict(label=DATASET, bytes=Path(args.input).stat().st_size,
                                     records=args.records, sort_mib=128 if args.smoke else 98304)
        base.SORT_AVAILABLE_GIB[DATASET] = 0 if args.smoke else 112
        super().__init__(args)
        plan = dict(dataset=DATASET, input=args.input,
                    input_bytes=base.DATASETS[DATASET]["bytes"], expected=self.expected,
                    apps=base.APPS, library=["count", "bam2bam"],
                    library_threads=[0, 2] if args.smoke else LIBRARY_THREADS,
                    library_order="ascending/descending/ascending; both operations per point",
                    apps_search="coarse once; top two plus prior configuration confirmed three times",
                    app_candidates=[2, 4] if args.smoke else base.CANDIDATES,
                    flagstat_candidates=[0, 2] if args.smoke else base.candidates(DATASET, "flagstat"),
                    collate_bins=66, compression_level=1, sort_total_mib=base.DATASETS[DATASET]["sort_mib"],
                    workflows=3, workflow_timing="sum four measured command walls per round; then median",
                    timing="GNU time process wall, real NVMe I/O included; validation excluded",
                    cache="not dropped; warm preparation and coarse samples; no O_DIRECT or fsync promise",
                    resource_wait_seconds=args.resource_wait_seconds,
                    resource_guard="112 GiB available for sort; 16 GiB otherwise; current CPU busy <=12.5%",
                    smoke=args.smoke, samtools=args.samtools, prepared_reference=args.reference_root)
        (self.root / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
        with (self.root / "environment.txt").open("a") as out:
            subprocess.run(["sha256sum", __file__, args.input], stdout=out, stderr=out, check=True)
        self.library_file = (self.root / "library_results.tsv").open("w", newline="")
        self.library_writer = csv.DictWriter(self.library_file, [
            "app", "threads", "repeat", "wall_s", "output_bytes", "records_checked",
            "maxrss_kb", "user_s", "sys_s", "cpu", "fs_inputs", "fs_outputs", "exit", "log"], delimiter="\t")
        self.library_writer.writeheader()
        self.library_file.flush()
        self.library_rows = []

    def input(self, dataset):
        if dataset != DATASET:
            raise RuntimeError("only the new mouse D3 is allowed")
        return Path(self.args.input)

    def wait_resources(self, disk_bytes=100 * base.GIB, memory_bytes=16 * base.GIB):
        if self.smoke:
            return
        deadline = time.monotonic() + self.args.resource_wait_seconds
        while True:
            free = shutil.disk_usage(str(self.root)).free
            if free < disk_bytes:
                raise RuntimeError("insufficient disk space: {} < {}".format(free, disk_bytes))
            available = base.memory_available()
            busy = current_cpu_busy()
            if available >= memory_bytes and busy <= .125:
                return
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise RuntimeError("resource wait deadline: available={} required={} CPU_busy={:.1%}".format(
                    available, memory_bytes, busy))
            self.log("WAIT_RESOURCES available={} required={} CPU_busy={:.1%} remaining_s={:.0f}".format(
                available, memory_bytes, busy, remaining))
            time.sleep(min(30, remaining))

    def check_sort(self, dataset, stderr):
        match = re.search(r"merging from (\d+) files", stderr.read_text())
        if not match or int(match.group(1)) != 0:
            raise RuntimeError("D3 mouse sort did not demonstrate in-memory mode: " + str(stderr))

    def validate(self, path, full_count=False, expected=None):
        name = path.name
        if full_count and ".filter." in name:
            expected = self.expected["kept"]
        count = super().validate(path, full_count, expected)
        if full_count and (".markdup." in name or name == "markdup.bam"):
            label = name + ".validation"
            flags = self.root / "logs" / (label + ".flagstat.json")
            command = [self.args.samtools, "flagstat", "-@", "8", "-O", "json", str(path)]
            self.record_command("final duplicates validation (outside timing)", command)
            with flags.open("w") as out:
                subprocess.run(command, env=self.env, stdout=out, check=True)
            data = json.loads(flags.read_text())
            marked = sum(group["duplicates"] for group in data.values())
            if marked != self.expected["duplicates"]:
                raise RuntimeError("duplicate total {} != {}".format(marked, self.expected["duplicates"]))
        if full_count and (".sort." in name or name in ("sorted.bam", "markdup.bam") or ".markdup." in name):
            stats = self.root / "logs" / (name + ".validation.stats")
            command = [self.args.samtools, "stats", "-@", "8", str(path)]
            self.record_command("coordinate validation, not performance", command)
            with stats.open("w") as out:
                subprocess.run(command, env=self.env, stdout=out, check=True)
            if not re.search(r"^SN\tis sorted:\t1(?:\s|$)", stats.read_text(), re.M):
                raise RuntimeError("coordinate order validation failed")
        return count

    def measured(self, label, command, stdout_path):
        result = super().measured(label, command, stdout_path)
        if label.startswith("confirm.") and ".bam2sam." in label and label.endswith(".r1"):
            out = Path(command[command.index("-o") + 1])
            count_command = [self.args.samtools, "view", "-@", "8", "-c", str(out)]
            self.record_command("complete SAM count (outside timing)", count_command)
            count = int(subprocess.check_output(count_command, env=self.env, universal_newlines=True))
            if count != self.expected["records"]:
                raise RuntimeError("full SAM record count differs")
            with (self.root / "sam_validations.jsonl").open("a") as output:
                output.write(json.dumps(dict(sample=label, records=count, path=str(out))) + "\n")
        return result

    def library(self):
        candidates = [0, 2] if self.smoke else LIBRARY_THREADS
        for repeat in range(1, 4):
            for threads in candidates[::(-1 if repeat == 2 else 1)]:
                for app in ("count", "bam2bam"):
                    label = "library.{}.t{}.r{}".format(app, threads, repeat)
                    out = self.root / "work" / (label + ".bam")
                    stdout = self.root / "logs" / (label + ".stdout")
                    self.wait_resources(disk_bytes=100 * base.GIB, memory_bytes=16 * base.GIB)
                    command = [self.args.samtools, "view", "--no-PG", "-@", str(threads)]
                    if app == "count":
                        command += ["-c", str(self.input(DATASET))]
                    else:
                        command += ["-O", "BAM,level=1", "-o", str(out), str(self.input(DATASET))]
                    timing, log, _ = self.measured(label, command, stdout)
                    if app == "count":
                        checked = int(stdout.read_text().strip())
                        if checked != self.expected["records"]:
                            raise RuntimeError("library count differs")
                    else:
                        checked = self.validate(out, repeat == 1, self.expected["records"])
                    row = dict(timing, app=app, threads=threads, repeat=repeat,
                               output_bytes=out.stat().st_size if app == "bam2bam" else 0,
                               records_checked=checked, log=str(log))
                    self.library_writer.writerow(row)
                    self.library_file.flush()
                    self.library_rows.append(row)
                    if app == "bam2bam":
                        out.unlink()
                    self.log("DONE {} wall_s={}".format(label, timing["wall_s"]))
        summary = []
        for app in ("count", "bam2bam"):
            for threads in candidates:
                values = [float(row["wall_s"]) for row in self.library_rows
                          if row["app"] == app and row["threads"] == threads]
                if len(values) != 3:
                    raise RuntimeError("missing library repetitions")
                summary.append(dict(app=app, threads=threads, samples_s=values,
                                    median_s=statistics.median(values), min_s=min(values), max_s=max(values)))
        (self.root / "library_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    def tune(self, dataset, app, index):
        search = [2, 4] if self.smoke else base.candidates(dataset, app)
        if self.smoke and app == "flagstat":
            search = [0, 2]
        coarse = [self.run_app(dataset, app, threads, "coarse", 1)
                  for threads in search[::(-1 if index % 2 else 1)]]
        search = base.shortlist(coarse, 2 if self.smoke else base.inherited_thread(dataset, app))
        rows = []
        for repeat in range(1, 4):
            for threads in search[::(-1 if repeat == 2 else 1)]:
                rows.append(self.run_app(dataset, app, threads, "confirm", repeat))
        choices = []
        for threads in search:
            times = [float(row["wall_s"]) for row in rows if int(row["threads"]) == threads]
            if len(times) != 3:
                raise RuntimeError("missing confirmation repetitions")
            choices.append(dict(app=app, threads=threads, samples_s=times, median_s=statistics.median(times)))
        selected = min(choices, key=lambda row: (row["median_s"], row["threads"]))
        self.selected[(dataset, app)] = selected["threads"]
        with (self.root / "summary.jsonl").open("a") as output:
            output.write(json.dumps(dict(configurations=choices, selected=selected)) + "\n")
        self.log("SELECTED {} threads={} median_s={}".format(app, selected["threads"], selected["median_s"]))

    def workflow(self, dataset, repeat):
        super().workflow(dataset, repeat)
        label = "workflow.{}.r{}".format(dataset, repeat)
        path = self.root / "logs" / (label + ".json")
        row = json.loads(path.read_text())
        total = sum(row["stage_times_s"])
        with (self.root / "workflows.jsonl").open("a") as output:
            output.write(json.dumps(dict(repeat=repeat, stage_times_s=row["stage_times_s"],
                                        command_wall_sum_s=total, runner_wall_s=float(row["wall_s"]))) + "\n")

    def reuse_reference(self):
        if not self.args.reference_root:
            return
        source = Path(self.args.reference_root)
        metadata = json.loads((source / "full.json").read_text())
        if any(metadata[key] != self.expected[key] for key in ("records", "kept")):
            raise RuntimeError("prepared reference metadata differs")
        if metadata.get("marked") != self.expected["duplicates"]:
            raise RuntimeError("reference duplicate count differs")
        prep = self.root / "prepared" / DATASET
        prep.mkdir()
        for target, name in (("name.bam", "collate.bam"), ("fixmate.bam", "fixmate.bam"),
                             ("fixmate.sorted.bam", "sort.bam")):
            path = source / "full-reference" / name
            self.validate(path)
            (prep / target).symlink_to(path)
        (self.root / "prepared_reference.json").write_text(json.dumps(dict(
            source=str(source), purpose="prepared inputs only; no previous timing samples imported"), indent=2) + "\n")

    def run(self):
        src = self.input(DATASET)
        if not self.smoke and (src.stat().st_size != INPUT_BYTES or self.args.records != RECORDS):
            raise RuntimeError("unexpected mouse D3 input")
        self.validate(src)
        self.reuse_reference()
        self.library()
        for index, app in enumerate(base.APPS):
            if app == "bam2sam":
                for repeat in range(1, 4):
                    self.workflow(DATASET, repeat)
                for name in ("name.bam", "fixmate.bam", "fixmate.sorted.bam"):
                    path = self.root / "prepared" / DATASET / name
                    if path.exists():
                        path.unlink()
            self.tune(DATASET, app, index)
        workflows = [json.loads(line) for line in (self.root / "workflows.jsonl").read_text().splitlines()]
        if len(workflows) != 3:
            raise RuntimeError("missing workflow rounds")
        (self.root / "workflow_summary.json").write_text(json.dumps(dict(
            command_wall_sums=[row["command_wall_sum_s"] for row in workflows],
            median_s=statistics.median(row["command_wall_sum_s"] for row in workflows),
            timing="sum four measured command walls in each round, then median"), indent=2) + "\n")
        sam = self.root / "prepared" / DATASET / "input.sam"
        if sam.exists():
            sam.unlink()
        self.log("COMPLETE")
        (self.root / "ALL_DONE").touch()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--samtools", default="/home/user_home/wzs/software/samtools-1.20/bin/samtools")
    parser.add_argument("--libdeflate", default="/home/user_home/wzs/software/libdeflate-1.20/lib64")
    parser.add_argument("--reference-root")
    parser.add_argument("--records", type=int, default=RECORDS)
    parser.add_argument("--kept", type=int, default=KEPT)
    parser.add_argument("--duplicates", type=int, default=DUPLICATES)
    parser.add_argument("--resource-wait-seconds", type=int, default=1800)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    args.workflows = True
    args.resume_from = None
    if not args.smoke and (args.records, args.kept, args.duplicates) != (RECORDS, KEPT, DUPLICATES):
        parser.error("formal mouse D3 metadata must match the verified input")
    bench = None
    try:
        bench = MouseBenchmark(args)
        bench.run()
    except Exception as error:
        if bench:
            (bench.root / "FAILED").write_text(str(error) + "\n")
            bench.log("FAILED " + str(error))
        raise
    finally:
        if bench:
            bench.writer_file.close()
            bench.library_file.close()
            bench.lock.close()


if __name__ == "__main__":
    main()
