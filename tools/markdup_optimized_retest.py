#!/usr/bin/env python3
"""Retest only markdup and actual four-command workflows with the frozen optimized binary."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shlex
import statistics
import time

import probe_nonhuman_d3 as probe
from markdup_halo_metrics import summarize as halo_metrics


BENCH = probe.SW_DATA + "/swbam_bench"
NEW_APP = BENCH + "/markdup_batch_tuning_20261108/final/bin/RabbitBAM-MPI"
NEW_SHA = "d7c8db1c1d9d41d4e816dcab4adb089314dd8d68f5792d060714b81db5efce1e"
BASE_APP = BENCH + "/ram_D2_processing_matrix_20261013_114500/bin/RabbitBAM-MPI"
BASE_SHA = "81b68788a9d897b80f494fd1d84945f15b2fe8c7791a28f3559779c77c683cb1"
MOUSE_FM = BENCH + "/ram_D3_mouse_probe_20261016_023500/bin/RabbitBAM-MPI-fixmate-singleton"
MOUSE_FM_SHA = "fc7f37cbc87a6fc1a65b01cb87037c012da3df7851e2e28bf082f7692ea64602"
SCRIPTS = ("markdup_optimized_retest.py", "probe_nonhuman_d3.py", "markdup_halo_metrics.py")
SPECS = {
    "D1": dict(records=177798571, marked=50084605,
        input=probe.SW_DATA + "/D1-WES-9G/initial/D1-WES-9G.coordinate.bam",
        sorted=BENCH + "/D1_apps_current/apps_D1_n6_20261001_031800/out/fixmate.sorted.bam",
        pools=("8G", "12G", "8G", "8G"), verify_pool="4G", memories=("24G", "32G")),
    "D2": dict(records=365079492, marked=8894488,
        input=probe.SW_DATA + "/D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam",
        sorted=BENCH + "/D2_apps_current/apps_D2_D3_workflows_n6_20261010_003937_D2/out/fixmate.sorted.bam",
        pools=("16G", "24G", "16G", "24G"), verify_pool="8G", memories=("48G", "48G")),
    "D3": dict(records=165373770, marked=33171753,
        input=probe.SW_DATA + "/D3-mouse-SRR7768905/initial/D3-mouse-SRR7768905.coordinate.bam",
        sorted=BENCH + "/ram_D3_mouse_probe_20261016_023500/out/full-fm_fc7f37cbc87a/sort.bam",
        pools=("16G", "24G", "16G", "24G"), verify_pool="8G", memories=("48G", "48G")),
    "D2_external": dict(records=365079492, marked=8894488,
        input=probe.SW_DATA + "/D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam",
        pools=("64G", "24G", "64G", "24G"), verify_pool="8G", memories=("8G", "8G")),
}


def binaries(case):
    return {"RabbitBAM-MPI": (BASE_APP, BASE_SHA),
            "RabbitBAM-MPI-fixmate": (MOUSE_FM, MOUSE_FM_SHA) if case == "D3" else (BASE_APP, BASE_SHA),
            "RabbitBAM-MPI-markdup": (NEW_APP, NEW_SHA)}


def build_plan(root, case):
    spec = SPECS[case]
    app, fm, md = [root + "/bin/" + name for name in binaries(case)]
    external = case == "D2_external"
    ordering = "external" if external else "memory"
    plan = []

    def add(kind, name, repeat, capacity, check, argv, cleanup=()):
        plan.append(dict(kind=kind, app=name, repeat=repeat, ranks=6,
            capacity=capacity, check=check, argv=argv, cleanup=list(cleanup),
            label="%s_%s_n6_r%d" % (kind, name, repeat)))

    def validate(name, repeat, output, duplicates=False, cleanup=()):
        add("validation", name, repeat, spec["verify_pool"],
            "duplicates" if duplicates else "sorted", [app, "stats", "--basic", "-i", output], cleanup)

    for repeat in (1, 2, 3):
        if not external:
            output = root + "/out/app_markdup_r%d.bam" % repeat
            add("app", "markdup", repeat, spec["pools"][3], "markdup",
                [md, "markdup", "-i", spec["sorted"], "-o", output,
                 "--stream", "-l", "300", "-m", "8G", "--compress-level", "1"] +
                ([] if repeat == 1 else ["--no-dump"]))
            if repeat == 1:
                validate("app_markdup", repeat, output, duplicates=True)

        folder = root + "/out/workflow_r%d/" % repeat
        name, fixed, sorted_bam, marked = [folder + name for name in
            ("name.bam", "fixmate.bam", "sorted.bam", "markdup.bam")]
        add("workflow", "collate", repeat, spec["pools"][0], ordering,
            [app, "collate", "-i", spec["input"], "-o", name, "-n", "66",
             "-m", spec["memories"][0], "--compress-level", "1"])
        add("workflow", "fixmate", repeat, spec["pools"][1], "write",
            [fm, "fixmate", "-i", name, "-o", fixed, "-m", "--compress-level", "1"])
        add("workflow", "sort", repeat, spec["pools"][2], ordering,
            [app, "sort", "-i", fixed, "-o", sorted_bam,
             "-m", spec["memories"][1], "--compress-level", "1"])
        validate("workflow_sort", repeat, sorted_bam)
        add("workflow", "markdup", repeat, spec["pools"][3], "markdup",
            [md, "markdup", "-i", sorted_bam, "-o", marked,
             "--stream", "-l", "300", "-m", "8G", "--compress-level", "1"])
        validate("workflow_markdup", repeat, marked, duplicates=True,
            cleanup=[name, fixed, sorted_bam] + ([] if repeat == 1 else [marked]))
    return plan


def sample_metrics(log, sample, spec):
    check = sample["check"]
    result = probe.metrics(log, 6, spec["records"],
        "sorted" if check == "duplicates" else check, marked=spec["marked"])
    storage = probe.rank_lines(log, "ram-storage", 6)
    capacity = int(sample["capacity"][:-1]) * 1024 ** 3
    if any(int(row["capacity"]) != capacity or int(row["peak"]) > capacity for row in storage):
        raise RuntimeError("RAM pool differs from frozen capacity")
    if sample["kind"] != "validation" and (result["output_bytes"] <= 0 or
            len({row["output_bytes"] for row in storage}) != 1):
        raise RuntimeError("missing or inconsistent complete RAM output")
    if check == "duplicates" and not re.search(
            r"^SN\treads duplicated:\t" + str(spec["marked"]) + r"(?:\s|$)", log, re.M):
        raise RuntimeError("full reread duplicate total differs")
    core = re.findall(r"^Complete the 4\.3 .* cost ([0-9.]+)", log, re.M)
    if len(core) != 1:
        raise RuntimeError("missing or ambiguous core timing")
    result["core_s"] = float(core[0])
    if check == "markdup":
        result["halo"] = halo_metrics(log, 6, spec["records"], spec["marked"])
    if check == "external":
        runs = re.findall(r"^\s*external(?:_actual)?_sum .*?\bruns=(\d+)", log, re.M)
        peaks = re.findall(r"^\s*external(?:_memory)?_sum .*?\btracked_peak_max=(\d+)", log, re.M)
        if ("selected mode=external" not in log or "selected mode=memory" in log or
                len(runs) != 1 or int(runs[0]) <= 6 or not peaks or
                not 0 < max(map(int, peaks)) <= 8 * 1024 ** 3 or
                result["temp_read_bytes"] <= 0 or result["temp_write_bytes"] <= 0):
            raise RuntimeError("missing actual bounded external runs or temporary I/O")
        result["runs"] = int(runs[0])
        result["work_peak_bytes"] = max(map(int, peaks))
    return result


def summaries(rows):
    labels = [row["label"] for row in rows]
    if len(set(labels)) != len(labels):
        raise RuntimeError("duplicate accepted sample")
    series = []
    samples = sorted([row for row in rows if row["kind"] == "app"], key=lambda row: row["repeat"])
    if len(samples) == 3 and [row["repeat"] for row in samples] == [1, 2, 3]:
        series.append(dict(app="markdup", samples=samples,
            processing_median_s=statistics.median(row["processing_max_s"] for row in samples),
            core_median_s=statistics.median(row["core_s"] for row in samples)))
    workflows = []
    for repeat in (1, 2, 3):
        stages = [row for row in rows if row["kind"] == "workflow" and row["repeat"] == repeat]
        checks = {row["app"] for row in rows if row["kind"] == "validation" and row["repeat"] == repeat}
        if len(stages) != 4 or not {"workflow_sort", "workflow_markdup"}.issubset(checks):
            continue
        if {row["app"] for row in stages} != {"collate", "fixmate", "sort", "markdup"}:
            raise RuntimeError("invalid four-command round")
        workflows.append(dict(repeat=repeat, stages=stages,
            processing_total_s=sum(row["processing_max_s"] for row in stages),
            mpi_interval_total_s=sum(row["e2e_max_s"] for row in stages)))
    result = dict(accepted=len(rows), series=series, workflows=workflows)
    if len(workflows) == 3:
        result["workflow_processing_median_s"] = statistics.median(row["processing_total_s"] for row in workflows)
    return result


class Retest(probe.Probe):
    def state(self, phase, **values):
        value = dict(phase=phase, case=self.args.case, swls_root=self.sw,
            local_time=time.strftime("%Y-%m-%d %H:%M:%S %z"), **values)
        temporary = self.root / "state.tmp"
        temporary.write_text(json.dumps(value, indent=2) + "\n")
        temporary.replace(self.root / "state.json")
        print(json.dumps(value), flush=True)

    def initialize(self):
        spec = SPECS[self.args.case]
        plan = build_plan(self.sw, self.args.case)
        frozen = binaries(self.args.case)
        scripts = {name: hashlib.sha256((Path(__file__).parent / name).read_bytes()).hexdigest()
                   for name in SCRIPTS}
        config = dict(case=self.args.case, spec=spec, nodes=self.args.nodes,
            binaries={name: dict(source=path, sha256=digest) for name, (path, digest) in frozen.items()},
            commands=plan, scripts=scripts, repetitions=3,
            scope="markdup only + actual four-command rounds; no SDK/other app series",
            timing="rank-max processing; excludes physical preload/dump; sum four stages then median",
            retention="standalone r1 and workflow r1 final BAM retained; other workflow files removed only after readback passes")
        config = json.loads(json.dumps(config))
        path = self.root / "plan.json"
        if self.args.resume:
            if json.loads(path.read_text()) != config:
                raise RuntimeError("frozen plan, scripts, nodes or binaries changed")
        else:
            if path.exists() or (self.root / "jobs.jsonl").exists():
                raise RuntimeError("use a new root or explicit --resume")
            lines = ["set -eu", "test ! -e " + shlex.quote(self.sw),
                "mkdir -p " + " ".join(shlex.quote(self.sw + "/" + name)
                    for name in ("bin", "logs", "out", "sessions", "scripts"))]
            for name, (source, digest) in frozen.items():
                lines += ["test $(sha256sum " + shlex.quote(source) + " | awk '{print $1}') = " + digest,
                          "cp -p " + shlex.quote(source) + " " + shlex.quote(self.sw + "/bin/" + name)]
            lines += ["mkdir -p " + " ".join(shlex.quote(self.sw + "/out/workflow_r%d" % repeat)
                                              for repeat in (1, 2, 3))]
            probe.remote(self.args.swls, ["bash", "-s"], script="\n".join(lines))
            path.write_text(json.dumps(config, indent=2) + "\n")
        actual = probe.remote(self.args.swls, ["sha256sum"] + [self.sw + "/bin/" + name for name in frozen])
        if {Path(p).name: digest for digest, p in (line.split() for line in actual.splitlines())} != {
                name: digest for name, (_, digest) in frozen.items()}:
            raise RuntimeError("frozen binary mismatch")
        (self.root / "binaries.sha256").write_text(actual)
        inputs = [spec["input"]] + ([spec["sorted"]] if "sorted" in spec else [])
        info = probe.remote(self.args.swls, ["stat", "-c", "%n\t%s\t%Y"] + inputs)
        manifest = self.root / "inputs.tsv"
        if self.args.resume and manifest.read_text() != info:
            raise RuntimeError("input size/mtime changed")
        manifest.write_text(info)
        if not self.args.resume:
            for name in SCRIPTS:
                probe.copy(str(Path(__file__).parent / name), self.args.swls + ":" + self.sw + "/scripts/")
            for name in ("plan.json", "binaries.sha256", "inputs.tsv"):
                probe.copy(str(self.root / name), self.args.swls + ":" + self.sw + "/")
        self.state("CONFIGURATION_FROZEN", planned=len(plan))
        return plan

    def execute(self):
        plan = self.initialize()
        spec = SPECS[self.args.case]
        rows = []
        for sample in plan:
            check = "sorted" if sample["check"] == "duplicates" else sample["check"]
            self.run(sample["label"], sample["capacity"], sample["argv"], spec, check)
            text = (self.root / "logs" / (sample["label"] + ".log")).read_text()
            result = sample_metrics(text, sample, spec)
            with (self.root / "results.tsv").open() as stream:
                raw = [row for row in csv.DictReader(stream, delimiter="\t") if row["label"] == sample["label"]]
            if len(raw) != 1:
                raise RuntimeError("missing/repeated raw sample")
            for field in ("processing_max_s", "e2e_max_s", "output_bytes"):
                if abs(float(raw[0][field]) - result[field]) > .000002:
                    raise RuntimeError("raw metrics differ from log")
            row = dict(kind=sample["kind"], app=sample["app"], repeat=sample["repeat"],
                label=sample["label"], job=raw[0]["job"], capacity=sample["capacity"], **result)
            rows.append(row)
            (self.root / "accepted_samples.json").write_text(json.dumps(rows, indent=2) + "\n")
            scalar = [{key: value for key, value in row.items() if key != "halo"} for row in rows]
            keys = list(dict.fromkeys(key for row in scalar for key in row))
            with (self.root / "accepted_samples.tsv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, keys, delimiter="\t")
                writer.writeheader(); writer.writerows(scalar)
            summary = summaries(rows)
            (self.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            # Only this new run's generated files can be removed, after readback.
            for output in sample["cleanup"]:
                if not output.startswith(self.sw + "/out/workflow_r") or "/../" in output:
                    raise RuntimeError("unsafe cleanup path")
            if sample["cleanup"]:
                probe.remote(self.args.swls, ["rm", "-f", "--"] + sample["cleanup"])
            self.state("SAMPLE_ACCEPTED", label=sample["label"], job=row["job"],
                accepted=len(rows), planned=len(plan))
        if len(rows) != len(plan) or len(summary["workflows"]) != 3 or (
                self.args.case != "D2_external" and len(summary["series"]) != 1):
            raise RuntimeError("incomplete formal series")
        self.state("COMPLETE", accepted=len(rows), summary=str(self.root / "summary.json"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--case", choices=SPECS, required=True)
    parser.add_argument("--nodes", help="Optional six-node allocation, recorded in every command")
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.a4, args.retry_no_start, args.fixmate_app = "unused", False, None
    args.poll_seconds, args.job_prefix = 30, "mdopt_" + args.case + "_"
    retest = Retest(args)
    try:
        retest.execute()
    except Exception as error:
        retest.state("FAILED", job=retest.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
