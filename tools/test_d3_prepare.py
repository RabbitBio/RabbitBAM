#!/usr/bin/env python3
"""Small regression for the private D3 preparation validator (no performance run)."""

import argparse
import gzip
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--samtools", required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="d3-prepare-test-") as temporary:
        root = Path(temporary)
        sequence, quality = "ACGT" * 25, "I" * 100
        fastq = []
        for mate in (1, 2):
            path = root / ("r%d.fastq.gz" % mate)
            with gzip.open(str(path), "wt") as output:
                output.write("@pair/{}\n{}\n+\n{}\n".format(mate, sequence, quality))
            fastq.append(str(path))
        result = subprocess.check_output([args.audit, "fastq"] + fastq + ["1"], text=True)
        stats = json.loads(result)
        assert stats["pairs"] == 1 and stats["bases"] == 200 and stats["paired_qnames"]

        sam = root / "input.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n"
            "pair\t99\t1\t1\t60\t100M\t=\t11\t110\t{}\t{}\tNM:i:0\tMD:Z:100\n"
            "pair\t147\t1\t11\t60\t100M\t=\t1\t-110\t{}\t{}\tNM:i:0\tMD:Z:100\n".format(
                sequence, quality, sequence, quality))
        bam = root / "output.bam"
        subprocess.run([args.samtools, "view", "-b", "-o", str(bam), str(sam)], check=True)
        result = subprocess.check_output([args.audit, "bam", str(bam), "2"], text=True)
        stats = json.loads(result)
        assert stats["records"] == stats["primary"] == stats["paired"] == 2
        assert stats["duplicates"] == 0 and stats["coordinate_sorted"]
        assert stats["record_aligned_bgzf"] and stats["crc"] and stats["eof"]
        assert len(stats["record_md5"]) == 32
        mismatch = subprocess.run([args.audit, "bam", str(bam), "3"], capture_output=True)
        assert mismatch.returncode != 0 and b"primary count" in mismatch.stderr
    print("PASS: FASTQ pair validation, complete BAM audit, and count gate")


if __name__ == "__main__":
    main()
