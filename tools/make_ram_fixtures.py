#!/usr/bin/env python3
"""Build block-aligned BAM and line-boundary SAM fixtures with Pysam."""
import pathlib
import struct
import sys

import pysam

root = pathlib.Path(sys.argv[1])
root.mkdir(parents=True, exist_ok=True)
text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
header = b"BAM\1" + struct.pack("<I", len(text)) + text
header += struct.pack("<II", 1, 5) + b"chr1\0" + struct.pack("<I", 1000)
# One mapped base, one 1M CIGAR, MAPQ 60, no mate or auxiliary data.
core = struct.pack("<iiIIiiii", 0, 0, (4681 << 16) | (60 << 8) | 2, 1, 1, -1, -1, 0)
payload = core + b"r\0" + struct.pack("<I", 16) + bytes([0x10, 30])
record = struct.pack("<I", len(payload)) + payload
for name, records in (("empty", b""), ("one", record)):
    with pysam.BGZFile(str(root / (name + ".bam")), "wb") as bgzf:
        bgzf.write(header)
        bgzf.flush()
        if records:
            bgzf.write(records)
    with pysam.AlignmentFile(str(root / (name + ".bam")), "rb") as bam:
        with pysam.AlignmentFile(str(root / (name + ".sam")), "w", header=bam.header) as sam:
            for alignment in bam:
                sam.write(alignment)
print("Block-aligned empty/single-record fixtures created")
