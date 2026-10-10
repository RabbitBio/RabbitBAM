#!/usr/bin/env python3
"""Validate ordinary four-step rounds and report the median of core sums."""

import argparse
import csv
import json
from pathlib import Path
import re
import statistics


STAGES = ("collate", "fixmate", "sort", "markdup")


def collect_rounds(rows, marked, dataset, repeats, records, expected_marked):
    samples = {}
    for row in rows:
        if row["kind"] != "bench":
            continue
        key = (int(row["repeat"]), row["app"])
        if (key in samples or key[1] not in STAGES or
                not 1 <= key[0] <= repeats or int(row["ranks"]) != 6 or
                int(row["records"]) != records or int(row["kept"]) != records or
                float(row["core_s"]) <= 0 or int(row["output_bytes"]) <= 0):
            raise ValueError("invalid or duplicate workflow sample: " + str(key))
        samples[key] = row
    if len(samples) != repeats * len(STAGES):
        raise ValueError("incomplete four-step workflow")
    result = []
    for repeat in range(1, repeats + 1):
        if marked[repeat] != (expected_marked, 0):
            raise ValueError("duplicate count differs from validated baseline")
        stage_rows = [samples[(repeat, app)] for app in STAGES]
        costs = [float(row["core_s"]) for row in stage_rows]
        walls = [r["submit_wall_s"] for r in stage_rows]
        row = dict(dataset=dataset, repeat=repeat,
                   core_sum_s=sum(costs),
                   stage_submit_wall_sum_s=(sum(int(w) for w in walls)
                                            if all(str(w).strip() for w in walls) else None),
                   marked=marked[repeat][0])
        row.update(dict((app + "_s", cost) for app, cost in zip(STAGES, costs)))
        result.append(row)
    return result


def summarize(root, dataset, repeats, records, expected_marked):
    with (root / "results.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    marked = {}
    for repeat in range(1, repeats + 1):
        text = (root / "logs" / ("bench_markdup_r%d.log" % repeat)).read_text()
        matches = re.findall(r"^\s+duplicates .*marked=(\d+).*removed=(\d+)", text, re.M)
        if len(matches) != 1:
            raise ValueError("missing unique markdup duplicate summary")
        marked[repeat] = tuple(map(int, matches[0]))
    result = collect_rounds(rows, marked, dataset, repeats, records, expected_marked)
    duplicate_counts = []
    for repeat in sorted(set((1, repeats))):
        checks = [r for r in rows if r["kind"] == "validation" and
                  r["app"] == "verify_markdup" and int(r["repeat"]) == repeat]
        if len(checks) != 1 or int(checks[0]["records"]) != records:
            raise ValueError("missing first/final complete output validation")
        text = (root / "logs" / ("validation_verify_markdup_r%d.log" % repeat)).read_text()
        sorted_values = re.findall(r"^SN\tis sorted:\t(\d+)", text, re.M)
        duplicates = re.findall(r"^SN\treads duplicated:\t(\d+)", text, re.M)
        if sorted_values != ["1"] or len(duplicates) != 1:
            raise ValueError("invalid final output statistics")
        duplicate_counts.append(int(duplicates[0]))
    if len(set(duplicate_counts)) != 1:
        raise ValueError("first/final duplicate flag counts differ")
    fields = ["dataset", "repeat"] + [s + "_s" for s in STAGES] + [
        "core_sum_s", "stage_submit_wall_sum_s", "marked"]
    with (root / "workflow_results.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(result)
    summary = dict(dataset=dataset, repeats=repeats, rounds=result,
                   median_core_sum_s=statistics.median(r["core_sum_s"] for r in result),
                   records=records, marked=expected_marked,
                   validated_duplicate_flags=duplicate_counts[0],
                   timing="median of per-round four-stage core sums; not real SSD wall time",
                   wall_diagnostic="sum of four submission walls; excludes validation, not whole workflow wall; null if a submission wall was lost")
    (root / "workflow_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


def self_test():
    rows = []
    for repeat, costs in enumerate(((1, 100, 2, 1), (100, 1, 2, 1), (1, 1, 100, 1)), 1):
        for app, cost in zip(STAGES, costs):
            rows.append(dict(kind="bench", app=app, repeat=repeat, ranks=6,
                             records=10, kept=10, core_s=cost, output_bytes=100,
                             submit_wall_s=200))
    marked = {1: (2, 0), 2: (2, 0), 3: (2, 0)}
    result = collect_rounds(rows, marked, "test", 3, 10, 2)
    assert statistics.median(r["core_sum_s"] for r in result) == 104
    assert sum(statistics.median(r[s + "_s"] for r in result) for s in STAGES) == 5
    missing_wall = [dict(row) for row in rows]
    missing_wall[0]["submit_wall_s"] = ""
    recovered = collect_rounds(missing_wall, marked, "test", 3, 10, 2)
    assert recovered[0]["core_sum_s"] == result[0]["core_sum_s"]
    assert recovered[0]["stage_submit_wall_sum_s"] is None
    assert recovered[1]["stage_submit_wall_sum_s"] == 800
    for invalid in (rows[:-1], rows + [rows[0]]):
        try:
            collect_rounds(invalid, marked, "test", 3, 10, 2)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid workflow samples accepted")
    try:
        collect_rounds(rows, marked, "test", 3, 10, 3)
    except ValueError:
        pass
    else:
        raise AssertionError("incorrect duplicate count accepted")
    print("workflow summary self-test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", type=Path)
    parser.add_argument("--dataset")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--records", type=int)
    parser.add_argument("--marked", type=int)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.root is None or args.dataset is None or args.records is None or args.marked is None:
        parser.error("root, --dataset, --records and --marked are required")
    summarize(args.root, args.dataset, args.repeats, args.records, args.marked)


if __name__ == "__main__":
    main()
