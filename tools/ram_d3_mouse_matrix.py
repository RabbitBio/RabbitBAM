#!/usr/bin/env python3
"""Frozen mouse D3 RAM matrix: SDK curves, eight apps, real four-step rounds."""

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
from resume_markdup_markers import StartupWait


RECORDS, KEPT, MARKED = 165373770, 134417603, 33171753
EXPECTED = dict(records=RECORDS, kept=KEPT, marked=MARKED)
DATA_ROOT = Path("/home/wzs/work/data/D3-mouse-SRR7768905")
BASELINE = DATA_ROOT / "ram_D3_mouse_probe_20261016_023500"
MARKER_CHECK = DATA_ROOT / "markdup_markers_20261017_135600"
CONFIRMATION = DATA_ROOT / "markdup_marker_confirm_20261017_140500"
BASE_SW = probe.SW_DATA + "/swbam_bench/" + BASELINE.name
MARKER_SW = probe.SW_DATA + "/swbam_bench/" + MARKER_CHECK.name
INPUTS = {
    "bam": probe.SW_DATA + "/" + probe.DATASET + "/initial/" + probe.INPUT,
    "sam": BASE_SW + "/out/full/records.sam",
    "collated": BASE_SW + "/out/full/collate.bam",
    "fixmate": BASE_SW + "/out/full-fm_fc7f37cbc87a/fixmate.bam",
    "sorted": BASE_SW + "/out/full-fm_fc7f37cbc87a/sort.bam",
}
BINARIES = {
    "RabbitBAM-MPI": (BASE_SW + "/bin/RabbitBAM-MPI", probe.BINARIES["RabbitBAM-MPI"]),
    "RabbitBAM-MPI-fixmate": (BASE_SW + "/bin/RabbitBAM-MPI-fixmate-singleton",
                              "fc7f37cbc87a6fc1a65b01cb87037c012da3df7851e2e28bf082f7692ea64602"),
    "RabbitBAM-MPI-markdup": (MARKER_SW + "/bin/markdup-new",
                              "c7fcfe9e5e06aa88550dac2463d3b73b2308ea19fc2c63f58808ab85d5acb297"),
    "swbam-sdk-optimized-count": (BASE_SW + "/bin/swbam-sdk-optimized-count",
                                  probe.BINARIES["swbam-sdk-optimized-count"]),
    "swbam-sdk-optimized-bam2bam": (BASE_SW + "/bin/swbam-sdk-optimized-bam2bam",
                                    probe.BINARIES["swbam-sdk-optimized-bam2bam"]),
}
SDK_POOLS = {1: ("16G", "48G"), 2: ("12G", "24G"), 4: ("8G", "16G"),
             6: ("8G", "16G"), 12: ("4G", "8G")}
DRIVERS = ("ram_d3_mouse_matrix.py", "test_ram_d3_mouse_matrix.py", "probe_nonhuman_d3.py",
           "resume_markdup_markers.py", "markdup_marker_recheck.py",
           "markdup_incremental_recheck.py", "resume_markdup_lookup.py",
           "markdup_lookup_recheck.py", "verify_markdup_ready_outputs.py")


def build_plan(root):
    app, fm, md, count, copy = [root + "/bin/" + name for name in BINARIES]
    plan = []

    def add(kind, name, repeat, ranks, capacity, check, argv):
        plan.append(dict(kind=kind, app=name, repeat=repeat, ranks=ranks,
                         capacity=capacity, check=check, argv=argv,
                         label="%s_%s_n%d_r%d" % (kind, name, ranks, repeat)))

    def validate(name, repeat, path, check="sorted"):
        add("validation", name, repeat, 6, "8G", check, [app, "stats", "--basic", "-i", path])

    for repeat in (1, 2, 3):
        ranks_order = (12, 6, 4, 2, 1) if repeat == 2 else (1, 2, 4, 6, 12)
        for ranks in ranks_order:
            read_pool, write_pool = SDK_POOLS[ranks]
            add("sdk", "count", repeat, ranks, read_pool, "any", [count, INPUTS["bam"]])
            out = root + "/out/sdk_n%d_r%d.bam" % (ranks, repeat)
            save = [] if repeat == ranks == 1 else ["--no-dump"]
            add("sdk", "bam2bam", repeat, ranks, write_pool, "write",
                [copy, INPUTS["bam"], "-o", out] + save)
            if not save:
                validate("sdk", repeat, out)

        common = ["--compress-level", "1", "--no-dump"]
        out = root + "/out/"
        add("app", "flagstat", repeat, 6, "8G", "any", [app, "flagstat", "-i", INPUTS["bam"]])
        add("app", "filter", repeat, 6, "16G", "filter", [app, "run_all", "-i", INPUTS["bam"],
            "-o", out + "filter_r%d.bam" % repeat, "--min-mapq", "30"] + common)
        add("app", "collate", repeat, 6, "16G", "memory", [app, "collate", "-i", INPUTS["bam"],
            "-o", out + "collate_r%d.bam" % repeat, "-n", "66", "-m", "48G"] + common)
        add("app", "fixmate", repeat, 6, "24G", "write", [fm, "fixmate", "-i", INPUTS["collated"],
            "-o", out + "fixmate_r%d.bam" % repeat, "-m"] + common)
        add("app", "sort", repeat, 6, "16G", "memory", [app, "sort", "-i", INPUTS["fixmate"],
            "-o", out + "sort_r%d.bam" % repeat, "-m", "48G"] + common)
        add("app", "markdup", repeat, 6, "24G", "markdup", [md, "markdup", "-i", INPUTS["sorted"],
            "-o", out + "markdup_r%d.bam" % repeat, "--stream", "-l", "300", "-m", "8G"] + common)
        add("app", "bam2sam", repeat, 6, "40G", "write", [app, "run_all", "-i", INPUTS["bam"],
            "-o", out + "bam2sam_r%d.sam" % repeat] + common)
        converted = out + "sam2bam_r%d.bam" % repeat
        add("app", "sam2bam", repeat, 6, "40G", "write", [app, "run_all", "-i", INPUTS["sam"],
            "-o", converted, "--compress-level", "1"] + ([] if repeat == 1 else ["--no-dump"]))
        if repeat == 1:
            validate("sam2bam", repeat, converted)

        folder = out + "workflow_r%d/" % repeat
        name, fixed, sorted_bam, marked = [folder + name for name in
                                          ("name.bam", "fixmate.bam", "sorted.bam", "markdup.bam")]
        add("workflow", "collate", repeat, 6, "16G", "memory", [app, "collate", "-i", INPUTS["bam"],
            "-o", name, "-n", "66", "-m", "48G", "--compress-level", "1"])
        add("workflow", "fixmate", repeat, 6, "24G", "write", [fm, "fixmate", "-i", name,
            "-o", fixed, "-m", "--compress-level", "1"])
        add("workflow", "sort", repeat, 6, "16G", "memory", [app, "sort", "-i", fixed,
            "-o", sorted_bam, "-m", "48G", "--compress-level", "1"])
        validate("workflow_sort", repeat, sorted_bam)
        add("workflow", "markdup", repeat, 6, "24G", "markdup", [md, "markdup", "-i", sorted_bam,
            "-o", marked, "--stream", "-l", "300", "-m", "8G", "--compress-level", "1"])
        validate("workflow_markdup", repeat, marked, "duplicates")
    return plan


def sample_metrics(log, sample):
    check = sample["check"]
    result = probe.metrics(log, sample["ranks"], RECORDS,
                           "sorted" if check == "duplicates" else check, KEPT, MARKED)
    storage = probe.rank_lines(log, "ram-storage", sample["ranks"])
    capacity = int(sample["capacity"][:-1]) * 1024 ** 3
    if any(int(row["capacity"]) != capacity or int(row["peak"]) > capacity for row in storage):
        raise RuntimeError("RAM pool differs from frozen capacity")
    if check in ("write", "filter", "memory", "markdup") and (
            result["output_bytes"] <= 0 or len({row["output_bytes"] for row in storage}) != 1):
        raise RuntimeError("missing or inconsistent complete RAM output")
    if check == "duplicates" and not re.search(r"^SN\treads duplicated:\t33171753(?:\s|$)", log, re.M):
        raise RuntimeError("workflow output reread duplicate total differs")
    if check == "markdup" and (
            "boundary halo=fallback" not in log or
            sorted(map(int, re.findall(r"\[md-window-reuse rank=(\d+)\]", log))) != list(range(6))):
        raise RuntimeError("D3 markdup path/counters differ from frozen probe")
    core = re.findall(r"^sdk_(?:count|bam2bam) .*\bcore=([0-9.]+)|"
                      r"^Complete the 4\.3 .* cost ([0-9.]+)", log, re.M)
    if len(core) != 1:
        raise RuntimeError("missing or ambiguous diagnostic core timing")
    result["core_s"] = float(next(value for value in core[0] if value))
    return result


def summaries(rows):
    groups = {}
    for row in rows:
        if row["kind"] in ("sdk", "app"):
            groups.setdefault((row["kind"], row["app"], row["ranks"]), []).append(row)
    series = []
    for (kind, app, ranks), values in groups.items():
        if len(values) != 3:
            continue
        if sorted(value["repeat"] for value in values) != [1, 2, 3]:
            raise RuntimeError("duplicate repetitions in formal series")
        item = dict(kind=kind, app=app, ranks=ranks)
        for field in ("processing_max_s", "e2e_max_s", "core_s"):
            samples = [float(value[field]) for value in sorted(values, key=lambda row: row["repeat"])]
            item[field] = dict(samples=samples, median=statistics.median(samples), min=min(samples), max=max(samples))
        series.append(item)
    workflows = []
    for repeat in (1, 2, 3):
        stages = [row for row in rows if row["kind"] == "workflow" and row["repeat"] == repeat]
        verified = {row["app"] for row in rows if row["kind"] == "validation" and row["repeat"] == repeat}
        if len(stages) != 4:
            continue
        if not {"workflow_sort", "workflow_markdup"}.issubset(verified):
            continue
        if {row["app"] for row in stages} != {"collate", "fixmate", "sort", "markdup"}:
            raise RuntimeError("invalid four-step round")
        workflows.append(dict(repeat=repeat, stages=[dict(app=row["app"], job=row["job"],
            processing_s=row["processing_max_s"], mpi_interval_s=row["e2e_max_s"]) for row in stages],
            processing_total_s=sum(row["processing_max_s"] for row in stages),
            mpi_interval_total_s=sum(row["e2e_max_s"] for row in stages)))
    summary = dict(accepted_commands=len(rows), series=series, workflows=workflows)
    if len(workflows) == 3:
        summary["workflow_processing_median_s"] = statistics.median(row["processing_total_s"] for row in workflows)
        summary["workflow_mpi_interval_median_s"] = statistics.median(row["mpi_interval_total_s"] for row in workflows)
    return summary


class Matrix(probe.Probe):
    def state(self, phase, **values):
        values = dict(phase=phase, local_time=time.strftime("%Y-%m-%d %H:%M:%S %z"),
                      swls_root=self.sw, scope="mouse D3 formal RAM matrix only", **values)
        temporary = self.root / "state.tmp"
        temporary.write_text(json.dumps(values, indent=2) + "\n")
        temporary.replace(self.root / "state.json")
        print(json.dumps(values), flush=True)

    def gates(self):
        for root in (BASELINE, MARKER_CHECK, CONFIRMATION):
            if json.loads((root / "state.json").read_text())["phase"] != "COMPLETE":
                raise RuntimeError("prior capacity/correctness/control chain incomplete: " + str(root))
        if not json.loads((BASELINE / "small-verified.json").read_text())["passed"]:
            raise RuntimeError("prior field/AUX comparison did not pass")
        metadata = json.loads((BASELINE / "full.json").read_text())
        if any(metadata[key] != value for key, value in EXPECTED.items()):
            raise RuntimeError("reference record/filter/duplicate totals changed")
        confirmation = json.loads((CONFIRMATION / "confirmation_summary.json").read_text())
        if confirmation["binary_sha256"] != BINARIES["RabbitBAM-MPI-markdup"][1]:
            raise RuntimeError("markdup confirmation used a different version")
        self.state("PRIOR_GATES_VERIFIED")

    def initialize(self):
        self.gates()
        plan = build_plan(self.sw)
        source_dir = Path(__file__).resolve().parent
        scripts = "".join(hashlib.sha256((source_dir / name).read_bytes()).hexdigest() + "  " + name + "\n"
                          for name in DRIVERS)
        configuration = dict(dataset=probe.DATASET, expected=EXPECTED, input_paths=INPUTS,
            prior_gates=[str(path) for path in (BASELINE, MARKER_CHECK, CONFIRMATION)],
            binaries={name: dict(source=source, sha256=digest) for name, (source, digest) in BINARIES.items()},
            sdk_pools={str(key): value for key, value in SDK_POOLS.items()}, commands=plan,
            timing="rank-max processing and MPI interval measured separately; core diagnostic only",
            storage="full RAM input/temp/output; physical preload and post-Finalize save outside intervals",
            workflow="four actual commands each round; intermediate files saved and reloaded; sum then median",
            validation="eight full rereads; stats not a performance item; prior full small-data field/AUX gate",
            scope="30 SDK + 24 app + 12 workflow + 8 validation; no D3 external sort; no D1/D2 rerun")
        config_text = json.dumps(configuration, indent=2) + "\n"
        if not self.args.resume:
            lines = ["set -eu", "test ! -e " + shlex.quote(self.sw),
                     "test $(df -Pk " + shlex.quote(probe.SW_DATA) + " | awk 'NR==2{print $4}') -ge 419430400"]
            for path in INPUTS.values():
                lines.append("test -s " + shlex.quote(path))
            lines += ["mkdir -p " + " ".join(shlex.quote(self.sw + "/" + name) for name in
                      ("bin", "logs", "out", "sessions", "scripts"))]
            for name, (source, digest) in BINARIES.items():
                lines += ["test $(sha256sum " + shlex.quote(source) + " | awk '{print $1}') = " + digest,
                          "cp -p " + shlex.quote(source) + " " + shlex.quote(self.sw + "/bin/" + name)]
            lines += ["mkdir -p " + " ".join(shlex.quote(self.sw + "/out/workflow_r%d" % repeat)
                                             for repeat in (1, 2, 3))]
            probe.remote(self.args.swls, ["bash", "-s"], script="\n".join(lines))
        for name, value in (("plan.json", config_text), ("scripts.sha256", scripts)):
            path = self.root / name
            if path.exists() and path.read_text() != value:
                raise RuntimeError("frozen matrix configuration changed: " + name)
            if self.args.resume and not path.exists():
                raise RuntimeError("missing frozen manifest: " + name)
            path.write_text(value)
        actual = probe.remote(self.args.swls, ["sha256sum"] + [self.sw + "/bin/" + name for name in BINARIES])
        if {Path(path).name: digest for digest, path in (line.split() for line in actual.splitlines())} != {
                name: digest for name, (_, digest) in BINARIES.items()}:
            raise RuntimeError("frozen binary mismatch")
        manifest = self.root / "binaries.sha256"
        if manifest.exists() and manifest.read_text() != actual:
            raise RuntimeError("binary manifest changed")
        manifest.write_text(actual)
        info = probe.remote(self.args.swls, ["stat", "-c", "%n\t%s\t%Y"] + list(INPUTS.values()))
        inputs = self.root / "inputs.tsv"
        if inputs.exists() and inputs.read_text() != info:
            raise RuntimeError("input sizes or modification times changed")
        inputs.write_text(info)
        if not self.args.resume:
            for name in DRIVERS:
                probe.copy(str(source_dir / name), self.args.swls + ":" + self.sw + "/scripts/")
            for name in ("plan.json", "scripts.sha256", "binaries.sha256", "inputs.tsv"):
                probe.copy(str(self.root / name), self.args.swls + ":" + self.sw + "/")
        probe.remote(self.args.swls, ["bash", "-s"], script="set -eu\ncd " +
                     shlex.quote(self.sw + "/scripts") + "\nsha256sum -c ../scripts.sha256\n")
        self.state("CONFIGURATION_FROZEN", planned_commands=len(plan))
        return plan

    def execute(self):
        plan = self.initialize()
        rows = []
        original = probe.remote

        def observed(host, argv, timeout=120, script=None):
            output = original(host, argv, timeout=timeout, script=script)
            query = re.fullmatch(r"bjobs (\d+) 2>/dev/null \|\| true\n", script or "")
            if query and re.search(r"^" + query[1] + r"\s+RUN\s", output, re.M):
                path = self.root / ("allocation-" + query[1] + ".txt")
                if not path.exists():
                    path.write_text(output)
            return output

        probe.remote = StartupWait(observed, self.root / "scheduler-startup.jsonl")
        try:
            for sample in plan:
                self.run(sample["label"], sample["capacity"], sample["argv"], EXPECTED,
                         "sorted" if sample["check"] == "duplicates" else sample["check"], sample["ranks"])
                log = (self.root / "logs" / (sample["label"] + ".log")).read_text()
                result = sample_metrics(log, sample)
                with (self.root / "results.tsv").open() as stream:
                    raw = [row for row in csv.DictReader(stream, delimiter="\t") if row["label"] == sample["label"]]
                if len(raw) != 1 or any(abs(float(raw[0][key]) - value) > .000001
                                        for key, value in result.items() if key != "core_s"):
                    raise RuntimeError("sample log and raw metrics differ")
                rows.append(dict(kind=sample["kind"], app=sample["app"], repeat=sample["repeat"],
                                 ranks=sample["ranks"], capacity=sample["capacity"], label=sample["label"],
                                 job=raw[0]["job"], **result))
                with (self.root / "accepted_samples.tsv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, list(rows[0]), delimiter="\t")
                    writer.writeheader()
                    writer.writerows(rows)
                summary = summaries(rows)
                (self.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
                self.state("SAMPLE_ACCEPTED", accepted=len(rows), planned=len(plan),
                           label=sample["label"], job=raw[0]["job"])
        finally:
            probe.remote = original
        if len(rows) != 74 or len(summary["series"]) != 18 or len(summary["workflows"]) != 3:
            raise RuntimeError("incomplete formal matrix")
        self.state("COMPLETE", accepted=len(rows), summary=str(self.root / "summary.json"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.a4, args.retry_no_start, args.fixmate_app = "unused", False, None
    matrix = Matrix(args)
    try:
        matrix.execute()
    except Exception as error:
        matrix.state("FAILED", job=matrix.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
