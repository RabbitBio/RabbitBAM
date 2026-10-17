#!/usr/bin/env python3
"""Isolated small/D3 verification and D1/D2 old-new checks; never a formal matrix."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shlex
import time

from probe_nonhuman_d3 import Probe, SW_DATA, FROZEN, copy, remote
from verify_markdup_ready_outputs import verify


PROBE_ROOT = SW_DATA + "/swbam_bench/ram_D3_mouse_probe_20261016_023500"
INPUTS = {
    "small": PROBE_ROOT + "/out/small-fm_fc7f37cbc87a/sort.bam",
    "D3": PROBE_ROOT + "/out/full-fm_fc7f37cbc87a/sort.bam",
    "D1": SW_DATA + "/swbam_bench/D1_apps_current/apps_D1_n6_20261001_031800/out/fixmate.sorted.bam",
    "D2": SW_DATA + "/swbam_bench/D2_apps_current/apps_D2_D3_workflows_n6_20261010_003937_D2/out/fixmate.sorted.bam",
}
EXPECTED = {
    "small": dict(records=200000, marked=22259),
    "D3": dict(records=165373770, marked=33171753),
    "D1": dict(records=177798571, marked=50084605),
    "D2": dict(records=365079492, marked=8894488),
}
OLD_SHA = "3dc30d0a9ac3ba2a384fec617b66f32d7801a53d389288db0f51f110c7726a51"


class Recheck(Probe):
    def state(self, phase, **values):
        value = dict(phase=phase, local_time=time.strftime("%Y-%m-%d %H:%M:%S %z"),
                     swls_root=self.sw, scope="markdup checks only; no formal matrix", **values)
        temporary = self.root / "state.tmp"
        temporary.write_text(json.dumps(value, indent=2) + "\n")
        temporary.replace(self.root / "state.json")
        print(json.dumps(value), flush=True)

    def initialize(self):
        if not self.args.resume:
            remote(self.args.swls, ["bash", "-s"], script="\n".join([
                "set -eu", "test ! -e " + shlex.quote(self.sw),
                "mkdir -p " + " ".join(shlex.quote(self.sw + "/" + name)
                                          for name in ("bin", "logs", "out", "sessions")),
                "cp -p " + shlex.quote(FROZEN + "/RabbitBAM-MPI-markdup") + " " +
                shlex.quote(self.sw + "/bin/markdup-old"),
            ]))
            copy(str(self.args.new_app.resolve()), self.args.swls + ":" + self.sw + "/bin/markdup-new")
        actual = remote(self.args.swls, ["sha256sum", self.sw + "/bin/markdup-old",
                                        self.sw + "/bin/markdup-new"])
        fingerprints = {Path(path).name: checksum for checksum, path in
                        (line.split() for line in actual.splitlines())}
        new_sha = hashlib.sha256(self.args.new_app.read_bytes()).hexdigest()
        if fingerprints != {"markdup-old": OLD_SHA, "markdup-new": new_sha}:
            raise RuntimeError("binary fingerprint mismatch")
        manifest = self.root / "binaries.sha256"
        if manifest.exists() and manifest.read_text() != actual:
            raise RuntimeError("frozen programs changed")
        manifest.write_text(actual)
        files = [Path(__file__), Path(__file__).with_name("probe_nonhuman_d3.py")]
        revision = "".join(hashlib.sha256(path.read_bytes()).hexdigest() + "  " + path.name + "\n"
                           for path in files)
        scripts = self.root / "scripts.sha256"
        if scripts.exists() and scripts.read_text() != revision:
            raise RuntimeError("frozen check scripts changed")
        scripts.write_text(revision)
        self.state("CONFIGURATION_FROZEN", fingerprints=fingerprints)

    def markdup(self, dataset, variant="new", dump=False):
        label = dataset + "_" + variant
        output = self.sw + "/out/" + label + ".bam"
        argv = [self.sw + "/bin/markdup-" + variant, "markdup", "-i", INPUTS[dataset],
                "-o", output, "--stream", "-l", "300", "-m", "8G", "--compress-level", "1"]
        if not dump:
            argv.append("--no-dump")
        capacity = {"small": "2G", "D1": "8G", "D2": "24G", "D3": "24G"}[dataset]
        self.run(label, capacity, argv, EXPECTED[dataset], "markdup")
        log = (self.root / "logs" / (label + ".log")).read_text()
        if dataset in ("D1", "D2"):
            if "boundary halo=enabled" not in log:
                raise RuntimeError("historical halo layout changed: " + label)
            if variant == "new" and len(re.findall(r"lookup=linear index_peak=0", log)) != 6:
                raise RuntimeError("normal path unexpectedly uses fallback index: " + label)
        if dataset == "D3" and len(re.findall(r"lookup=indexed index_peak=\d+", log)) != 6:
            raise RuntimeError("D3 fallback index was not exercised")
        return output

    def checks(self):
        output = self.markdup("small", dump=True)
        local_output = self.root / "small.new.bam"
        copy(self.args.swls + ":" + output, str(local_output))
        verify([str(self.args.small_reference.resolve()), str(local_output)], 200000, 22259)
        self.state("SMALL_RECORD_BYTES_VERIFIED")

        output = self.markdup("D3", dump=True)
        self.run("D3_validate", "8G", [self.sw + "/bin/markdup-new", "stats", "--basic",
                 "-i", output], EXPECTED["D3"], "sorted")
        log = (self.root / "logs/D3_validate.log").read_text()
        if not re.search(r"^SN\treads duplicated:\t33171753(?:\s|$)", log, re.M):
            raise RuntimeError("full D3 reread duplicate total differs")

        for dataset, variants in (("D1", ("old", "new")), ("D2", ("new", "old"))):
            for variant in variants:
                self.markdup(dataset, variant)
        with (self.root / "results.tsv").open() as stream:
            samples = {row["label"]: float(row["processing_max_s"])
                       for row in csv.DictReader(stream, delimiter="\t")}
        ratios = {dataset: samples[dataset + "_new"] / samples[dataset + "_old"]
                  for dataset in ("D1", "D2")}
        summary = dict(processing_s=samples, new_over_old=ratios,
                       D3_previous_probe_s=388.807467,
                       D3_speedup=388.807467 / samples["D3_new"],
                       note="single paired probes, not formal medians; no samples replaced")
        (self.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        needs_review = any(ratio > 1.05 for ratio in ratios.values())
        self.state("CHECKS_NEED_PERFORMANCE_REVIEW" if needs_review else "COMPLETE", **summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--new-app", type=Path, required=True)
    parser.add_argument("--small-reference", type=Path, required=True)
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.a4 = "unused"
    args.retry_no_start = False
    args.fixmate_app = None
    check = Recheck(args)
    try:
        check.initialize()
        check.checks()
    except Exception as error:
        check.state("FAILED", job=check.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
