#!/usr/bin/env python3
"""D3 fallback marker compaction: full output, reread, D1/D2 controls."""

import argparse
import hashlib
import json
from pathlib import Path
import re

from markdup_incremental_recheck import IncrementalCheck


class MarkerCheck(IncrementalCheck):
    def initialize(self):
        super().initialize()
        revision = hashlib.sha256(Path(__file__).read_bytes()).hexdigest() + "  " + Path(__file__).name + "\n"
        manifest = self.root / "marker-driver.sha256"
        if manifest.exists() and manifest.read_text() != revision:
            raise RuntimeError("marker driver changed")
        manifest.write_text(revision)

    def markdup(self, dataset, variant="new", dump=False):
        output = super().markdup(dataset, variant, dump)
        log = (self.root / "logs" / (dataset + "_" + variant + ".log")).read_text()
        rows = re.findall(r"\[md-window-reuse rank=(\d+)\] marker_visits=(\d+) "
                          r"marker_retained=(\d+) marker_peak=(\d+)", log)
        if variant == "new" and dataset == "D3":
            if sorted(int(row[0]) for row in rows) != list(range(6)):
                raise RuntimeError("missing fallback marker counters")
            if not any(int(visits) > int(retained) for _, visits, retained, _ in rows):
                raise RuntimeError("D3 did not exercise marker compaction")
        if dataset in ("D1", "D2") and rows:
            raise RuntimeError("normal halo unexpectedly compacts fallback markers")
        return output

    def state(self, phase, **values):
        if "D3_previous_probe_s" in values:
            values["D3_previous_probe_s"] = 241.256594
            values["D3_speedup"] = 241.256594 / values["processing_s"]["D3_new"]
            values["samtools_D3_median_wall_s"] = 132.50
            values["note"] = "single probes; historical previous D3; D1/D2 paired controls; formal samples unchanged"
            (self.root / "summary.json").write_text(json.dumps(values, indent=2) + "\n")
        super().state(phase, **values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--new-app", type=Path, required=True)
    parser.add_argument("--small-reference", type=Path, required=True)
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    args.retry_no_start, args.a4, args.fixmate_app = False, "unused", None
    check = MarkerCheck(args)
    try:
        check.initialize()
        check.checks()
    except Exception as error:
        check.state("FAILED", job=check.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
