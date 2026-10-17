#!/usr/bin/env python3
"""Confirm the frozen D3 markdup timing twice after the first full reread."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import statistics
import time

from markdup_lookup_recheck import EXPECTED, INPUTS
from markdup_marker_recheck import MarkerCheck
from probe_nonhuman_d3 import metrics


def accepted_first(root, checksum):
    binary = root / "binaries.sha256"
    results = root / "results.tsv"
    if not binary.exists() or not results.exists():
        return None
    fingerprints = {Path(path).name: digest for digest, path in
                    (line.split() for line in binary.read_text().splitlines())}
    if fingerprints.get("markdup-new") != checksum:
        raise RuntimeError("first sample used a different binary")
    with results.open() as stream:
        rows = {row["label"]: row for row in csv.DictReader(stream, delimiter="\t")}
    if "D3_new" not in rows or "D3_validate_retry1" not in rows:
        return None
    for label, check in (("D3_new", "markdup"), ("D3_validate_retry1", "sorted")):
        log = (root / "logs" / (label + ".log")).read_text()
        sample = metrics(log, 6, EXPECTED["D3"]["records"], check,
                         marked=EXPECTED["D3"]["marked"])
        if abs(sample["processing_max_s"] - float(rows[label]["processing_max_s"])) > 1e-6:
            raise RuntimeError("first sample summary differs from log")
        if check == "sorted" and not re.search(r"^SN\treads duplicated:\t33171753(?:\s|$)", log, re.M):
            raise RuntimeError("first output reread duplicate total differs")
    return rows["D3_new"]


def summarize(samples):
    if len(samples) != 3:
        raise RuntimeError("expected exactly three independent measurements")
    processing = [float(row["processing_max_s"]) for row in samples]
    e2e = [float(row["e2e_max_s"]) for row in samples]
    return dict(processing_samples_s=processing, processing_median_s=statistics.median(processing),
                processing_range_s=[min(processing), max(processing)],
                mpi_interval_samples_s=e2e, mpi_interval_median_s=statistics.median(e2e),
                samtools_wall_median_s=132.50,
                target_met=statistics.median(processing) <= 132.50,
                note="SWBAM RAM processing and MPI interval; samtools real-NVMe wall; preparation and post-Finalize save excluded only from SWBAM intervals")


class Confirmation(MarkerCheck):
    def initialize(self):
        super().initialize()
        revision = hashlib.sha256(Path(__file__).read_bytes()).hexdigest() + "  " + Path(__file__).name + "\n"
        manifest = self.root / "confirmation-driver.sha256"
        if manifest.exists() and manifest.read_text() != revision:
            raise RuntimeError("confirmation driver changed")
        manifest.write_text(revision)

    def checks(self):
        checksum = hashlib.sha256(self.args.new_app.read_bytes()).hexdigest()
        self.state("WAITING_FIRST_FULL_REREAD", first_root=str(self.args.first_root.resolve()))
        deadline = time.monotonic() + 2 * 3600
        first = None
        while time.monotonic() < deadline:
            first = accepted_first(self.args.first_root, checksum)
            if first:
                break
            status = self.args.first_root / "state.json"
            if status.exists() and json.loads(status.read_text()).get("phase") == "FAILED":
                raise RuntimeError("first sample verification chain failed")
            time.sleep(15)
        if not first:
            raise RuntimeError("first full reread not accepted within two hours")
        for repeat in (2, 3):
            self.run("D3_confirm_r%d" % repeat, "24G", [self.sw + "/bin/markdup-new",
                "markdup", "-i", INPUTS["D3"], "-o", self.sw + "/out/D3_confirm_r%d.bam" % repeat,
                "--stream", "-l", "300", "-m", "8G", "--compress-level", "1", "--no-dump"],
                EXPECTED["D3"], "markdup")
            log = (self.root / ("logs/D3_confirm_r%d.log" % repeat)).read_text()
            if "boundary halo=fallback" not in log or len(re.findall(r"lookup=indexed index_peak=", log)) != 6:
                raise RuntimeError("confirmation unexpectedly changed algorithm path")
            if len(re.findall(r"\[md-window-reuse rank=\d+\]", log)) != 6:
                raise RuntimeError("missing marker counters in confirmation")
        with (self.root / "results.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        result = summarize([first] + rows)
        result["first_root"] = str(self.args.first_root.resolve())
        result["binary_sha256"] = checksum
        (self.root / "confirmation_summary.json").write_text(json.dumps(result, indent=2) + "\n")
        self.state("COMPLETE", **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--first-root", type=Path, required=True)
    parser.add_argument("--new-app", type=Path, required=True)
    parser.add_argument("--small-reference", type=Path, required=True)
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.retry_no_start, args.a4, args.fixmate_app = False, "unused", None
    check = Confirmation(args)
    try:
        check.initialize()
        check.checks()
    except Exception as error:
        check.state("FAILED", job=check.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
