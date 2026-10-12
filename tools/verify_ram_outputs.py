#!/usr/bin/env python3
"""Independent HTSlib/Pysam verification of a small RAM regression output set."""
import hashlib
import itertools
import pathlib
import sys

import pysam


def summary(path, coordinate=False):
    count = duplicates = scores = 0
    fingerprint = 0
    previous = (-1, -1)
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for record in bam:
            count += 1
            duplicates += record.is_duplicate
            scores += record.has_tag("ms")
            fingerprint = (fingerprint + int.from_bytes(
                hashlib.blake2b(record.to_string().encode(), digest_size=16).digest(), "little"
            )) % (1 << 128)
            key = (record.reference_id if record.reference_id >= 0 else 1 << 31,
                   record.reference_start)
            if coordinate:
                assert key >= previous, (path, count, previous, key)
            previous = key
    result = (count, duplicates, scores, fingerprint)
    print(path.name, result)
    return result


def check_collate(path):
    seen = set()
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for name, _ in itertools.groupby(bam, key=lambda record: record.query_name):
            digest = hashlib.blake2b(name.encode(), digest_size=16).digest()
            assert digest not in seen, (path, "noncontiguous QNAME", name)
            seen.add(digest)
    print("QNAME continuity PASS", len(seen))


def check_mates(path):
    checked = 0
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for name, group in itertools.groupby(bam, key=lambda record: record.query_name):
            primary = [r for r in group if r.is_paired and not r.is_secondary and not r.is_supplementary]
            if len(primary) != 2 or primary[0].is_unmapped or primary[1].is_unmapped:
                continue
            a, b = primary
            assert a.next_reference_id == b.reference_id and b.next_reference_id == a.reference_id, name
            assert a.next_reference_start == b.reference_start and b.next_reference_start == a.reference_start, name
            assert a.has_tag("MC") and b.has_tag("MC") and a.has_tag("ms") and b.has_tag("ms"), name
            assert a.get_tag("MC") == b.cigarstring and b.get_tag("MC") == a.cigarstring, name
            assert a.get_tag("ms") == sum(q for q in b.query_qualities if q >= 15), name
            assert b.get_tag("ms") == sum(q for q in a.query_qualities if q >= 15), name
            checked += 1
    assert checked > 0
    print("mate coordinates, MC and ms PASS", checked)


def main():
    source, root = map(pathlib.Path, sys.argv[1:])
    baseline = summary(source)
    assert summary(root / "sdk.bam") == baseline
    assert summary(root / "roundtrip.bam") == baseline
    with pysam.AlignmentFile(str(source), "rb") as bam:
        kept = sum(record.mapping_quality >= 30 for record in bam)
    assert summary(root / "filter.bam")[0] == kept
    assert summary(root / "collate.bam") == baseline
    check_collate(root / "collate.bam")
    check_mates(root / "fixmate.bam")
    fixed = summary(root / "fixmate.bam")
    assert summary(root / "sort.bam", coordinate=True) == fixed
    final = summary(root / "markdup.bam", coordinate=True)
    assert final[0] == baseline[0]
    assert final[1] > 0
    assert summary(root / "workflow.bam", coordinate=True) == final
    assert summary(root / "pipeline.bam", coordinate=True) == final
    print("Independent RAM output verification PASS; pysam", pysam.__version__)


if __name__ == "__main__":
    main()
