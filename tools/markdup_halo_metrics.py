#!/usr/bin/env python3
"""Summarize successful markdup halo probes; count record candidates separately from markers."""

import argparse
import json
import re
from pathlib import Path


def summarize(text, ranks, expected_records, expected_duplicates):
    if re.search(r"ERROR:|MPI_Abort|MPI_ABORT|Segmentation fault", text):
        raise ValueError("Log contains an execution error")

    def rows(prefix):
        result = {}
        for line in text.splitlines():
            if not line.startswith(prefix):
                continue
            fields = dict(re.findall(r"([a-z_]+)=([^\s\]]+)", line))
            rank = int(fields.pop("rank"))
            if rank not in range(ranks) or rank in result:
                raise ValueError("Invalid or repeated rank: " + prefix)
            result[rank] = fields
        return result

    def complete(values, label):
        if set(values) != set(range(ranks)):
            raise ValueError("Missing ranks: " + label)

    def total(values, key):
        numbers = [int(row[key]) for row in values.values()]
        if any(value < 0 for value in numbers):
            raise ValueError("Negative metric: " + key)
        return sum(numbers)

    progress = rows("[swbam-processing rank=")
    comm = rows("[md-comm rank=")
    halo = rows("[md-halo rank=")
    stats = {}
    for line in text.splitlines():
        match = re.match(r"\[rank (\d+)\] markdup_stats (.*)", line)
        if match:
            rank = int(match[1])
            if rank in stats or rank not in range(ranks):
                raise ValueError("Invalid or repeated markdup_stats rank")
            stats[rank] = dict(re.findall(r"([a-z_]+)=(\d+)", match[2]))
    for values, label in ((progress, "processing"), (comm, "communication"), (stats, "stats")):
        complete(values, label)
    if any(int(row["status"]) != 0 or float(row["seconds"]) <= 0 for row in progress.values()):
        raise ValueError("Processing did not succeed on all ranks")

    finish = re.findall(r"^OptimizedMarkdupMPI finished\. ranks=(\d+) blocks=\d+ records=(\d+)", text, re.M)
    duplicates = re.findall(r"^  duplicates .*?marked=(\d+)", text, re.M)
    if finish != [(str(ranks), str(expected_records))] or duplicates != [str(expected_duplicates)]:
        raise ValueError("Record or duplicate totals differ from the expected baseline")
    enabled = {int(row["halo_enabled"]) for row in comm.values()}
    if enabled == {1}:
        complete(halo, "halo")
        if any(int(row["enabled"]) != 1 for row in halo.values()):
            raise ValueError("Inconsistent halo status")
        partial = {int(row.get("partial", "0")) for row in halo.values()}
        if partial not in ({0}, {1}):
            raise ValueError("Inconsistent partial halo status")
        mode = "partial" if partial == {1} else "enabled"
        if mode == "partial" and (total(halo, "final_keep") != 0 or
                                  any(int(row["conservative"]) != 1 for row in comm.values())):
            raise ValueError("Partial halo must retain conservative winners without final KEEP")
    elif enabled == {0} and not halo:
        mode = "fallback"
    else:
        raise ValueError("Inconsistent halo status")

    main_candidates = total(comm, "main_candidates")
    main_markers = total(comm, "main_markers")
    if any(int(row["main_markers"]) > int(row["main_candidates"]) for row in comm.values()):
        raise ValueError("Marker count exceeds candidate count")
    pre_bytes = total(halo, "candidate_bytes")
    pre_results = total(halo, "result_bytes")
    main_bytes = total(comm, "main_candidate_bytes")
    main_results = total(comm, "main_result_bytes")
    if pre_bytes + main_bytes != total(stats, "mpi_candidate_bytes") or \
            pre_results + main_results != total(stats, "result_bytes"):
        raise ValueError("Phase communication bytes do not add up to cumulative stats")
    core = re.findall(r"^Complete the 4\.3 OptimizedMarkdupMPI cost ([0-9.]+)", text, re.M)
    if len(core) != 1:
        raise ValueError("Missing or repeated core time")

    return dict(
        ranks=ranks, records=expected_records, duplicates=expected_duplicates,
        halo=mode, conservative=any(int(row["conservative"]) for row in comm.values()),
        core_seconds=float(core[0]),
        processing_seconds=max(float(row["seconds"]) for row in progress.values()),
        preprocessing=dict(candidate_payload_bytes=pre_bytes, decision_payload_bytes=pre_results,
                           final_keep=total(halo, "final_keep"), final_dup=total(halo, "final_dup"),
                           pending=total(halo, "pending"), marker_pending=total(halo, "marker_pending")),
        main_loop=dict(record_candidate_entries=main_candidates - main_markers,
                       marker_entries=main_markers, all_candidate_entries=main_candidates,
                       candidate_payload_bytes=main_bytes, decision_payload_bytes=main_results),
        per_rank={str(rank): comm[rank] for rank in range(ranks)},
        byte_scope="sent candidate structs + QNAME; includes markers; excludes MPI count/status traffic",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--ranks", type=int, default=6)
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--duplicates", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.log.read_text(), args.ranks, args.records, args.duplicates), indent=2))


if __name__ == "__main__":
    main()
