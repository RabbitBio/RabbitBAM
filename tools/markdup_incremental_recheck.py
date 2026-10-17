#!/usr/bin/env python3
"""Frozen incremental-boundary probe, reread, and normal-path paired checks."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

from markdup_lookup_recheck import EXPECTED, copy, verify
from resume_markdup_lookup import Recovery


class IncrementalCheck(Recovery):
    def initialize(self):
        super().initialize()
        repo = Path(__file__).resolve().parent.parent
        paths = [Path(__file__), repo / "tools/resume_markdup_lookup.py",
                 repo / "swbam-app/algorithms/markdup/swbam_markdup_mpi.cpp",
                 repo / "swbam-app/algorithms/markdup/swbam_markdup_stream_mpi.cpp",
                 repo / "swbam-app/algorithms/markdup/markdup-boundary-progress.h",
                 repo / "swbam-app/algorithms/markdup/markdup-pending-index.h",
                 repo / "swbam-app/include/swbam_mpi.h", repo / "tools/markdup_window_test.cpp"]
        revision = "".join(hashlib.sha256(path.read_bytes()).hexdigest() + "  " +
                           str(path.relative_to(repo)) + "\n" for path in paths)
        manifest = self.root / "incremental-sources.sha256"
        if manifest.exists() and manifest.read_text() != revision:
            raise RuntimeError("incremental probe sources changed")
        manifest.write_text(revision)

    def markdup(self, dataset, variant="new", dump=False):
        output = super().markdup(dataset, variant, dump)
        log = (self.root / "logs" / (dataset + "_" + variant + ".log")).read_text()
        matches = re.findall(r"\[md-boundary-state rank=(\d+)\] mode=incremental", log)
        if variant == "new" and dataset == "D3" and sorted(map(int, matches)) != list(range(6)):
            raise RuntimeError("incremental boundaries not exercised on all ranks")
        if dataset in ("D1", "D2") and matches:
            raise RuntimeError("normal halo path unexpectedly used incremental boundaries")
        return output

    def checks(self):
        output = self.markdup("small", dump=True)
        local = self.root / "small.new.bam"
        copy(self.args.swls + ":" + output, str(local))
        verify([str(self.args.small_reference.resolve()), str(local)], 200000, 22259)
        self.state("SMALL_RECORD_BYTES_VERIFIED")
        output = self.markdup("D3", dump=True)
        self.run("D3_validate", "8G", [self.sw + "/bin/markdup-new", "stats", "--basic",
                 "-i", output], EXPECTED["D3"], "sorted")
        log = (self.root / "logs/D3_validate.log").read_text()
        if not re.search(r"^SN\treads duplicated:\t33171753(?:\s|$)", log, re.M):
            raise RuntimeError("full incremental D3 duplicate total differs")
        for dataset, variants in (("D1", ("old", "new")), ("D2", ("new", "old"))):
            for variant in variants:
                self.markdup(dataset, variant)
        with (self.root / "results.tsv").open() as stream:
            samples = {row["label"]: float(row["processing_max_s"])
                       for row in csv.DictReader(stream, delimiter="\t")}
        ratios = {dataset: samples[dataset + "_new"] / samples[dataset + "_old"]
                  for dataset in ("D1", "D2")}
        summary = dict(processing_s=samples, new_over_old=ratios,
                       D3_previous_probe_s=292.167561,
                       D3_speedup=292.167561 / samples["D3_new"],
                       note="single probes; D3 previous sample is historical, not a paired ablation; no formal samples replaced")
        (self.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        self.state("CHECKS_NEED_PERFORMANCE_REVIEW" if any(ratio > 1.05 for ratio in ratios.values())
                   else "COMPLETE", **summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--new-app", type=Path, required=True)
    parser.add_argument("--small-reference", type=Path, required=True)
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.retry_no_start, args.a4, args.fixmate_app = False, "unused", None
    check = IncrementalCheck(args)
    try:
        check.initialize()
        check.checks()
    except Exception as error:
        check.state("FAILED", job=check.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
