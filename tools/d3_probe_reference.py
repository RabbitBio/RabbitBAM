#!/usr/bin/env python3
"""One-off D3 samtools references and independent small-output checks on a4."""

import argparse
from collections import Counter
import hashlib
import itertools
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def command(args, root, label, argv, output=None):
    started = time.monotonic()
    with (root / "commands.jsonl").open("a") as log:
        log.write(json.dumps(dict(label=label, argv=argv)) + "\n")
    with (root / "logs" / (label + ".stderr")).open("w") as err:
        if output is None:
            subprocess.run(argv, check=True, stderr=err)
        else:
            with output.open("wb") as out:
                subprocess.run(argv, check=True, stdout=out, stderr=err)
    with (root / "samples.jsonl").open("a") as log:
        log.write(json.dumps(dict(label=label, wall_s=time.monotonic() - started)) + "\n")


def rows(samtools, path):
    process = subprocess.Popen([samtools, "view", "--no-PG", str(path)],
                               stdout=subprocess.PIPE, universal_newlines=True)
    try:
        result = []
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            result.append(tuple(fields[:11] + sorted(fields[11:])))
        if process.wait() != 0:
            raise RuntimeError("cannot fully parse " + str(path))
        return result
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.terminate()
            process.wait()


def digest(records):
    result = hashlib.sha256()
    for record in records:
        result.update(("\t".join(record) + "\n").encode())
    return result.hexdigest()


def four_steps(args, root, prefix, src, threads):
    s = args.samtools
    folder = root / prefix
    folder.mkdir()
    name, fixed, ordered, marked = [folder / (label + ".bam")
                                  for label in ("collate", "fixmate", "sort", "markdup")]
    command(args, root, prefix + ".collate", [s, "collate", "--no-PG", "-@", str(threads),
            "-n", "66", "-l", "1", "-T", str(folder / "temp.collate"), "-o", str(name), str(src)])
    command(args, root, prefix + ".fixmate", [s, "fixmate", "--no-PG", "-@", str(threads),
            "-m", "--output-fmt", "BAM,level=1", str(name), str(fixed)])
    sort_memory = "4096M" if prefix == "full-reference" else "1024M"
    command(args, root, prefix + ".sort", [s, "sort", "--no-PG", "-@", str(threads),
            "-m", sort_memory, "-l", "1", "-T", str(folder / "temp.sort"),
            "-o", str(ordered), str(fixed)])
    command(args, root, prefix + ".markdup", [s, "markdup", "--no-PG", "-@", str(threads),
            "-l", "300", "--output-fmt", "BAM,level=1", str(ordered), str(marked)])
    command(args, root, prefix + ".quickcheck", [s, "quickcheck", "-v"] +
            [str(path) for path in (name, fixed, ordered, marked)])
    return marked


def prepare_small(args, root):
    if (root / "small.json").exists():
        raise RuntimeError("refusing to overwrite small reference")
    source = root / "small.sam"
    with (root / "logs" / "extract.stderr").open("w") as err:
        process = subprocess.Popen([args.samtools, "view", "--no-PG", "-h", args.input],
                                   stdout=subprocess.PIPE, stderr=err)
        count = 0
        try:
            with source.open("wb") as out:
                for line in process.stdout:
                    if line.startswith(b"@"):
                        out.write(line)
                    elif count < args.small_records:
                        out.write(line)
                        count += 1
                    else:
                        break
            # A deliberately truncated, whole-record prefix, not a failed full download.
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
            status = process.wait(timeout=30)
            if count != args.small_records or status not in (0, 1, -signal.SIGPIPE, -signal.SIGTERM):
                raise RuntimeError("cannot extract requested complete-record prefix")
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
    bam = root / "small.bam"
    command(args, root, "small.encode", [args.samtools, "view", "--no-PG", "-b",
            "-o", str(bam), str(source)])
    baseline = rows(args.samtools, bam)
    final = four_steps(args, root, "small-reference", bam, 4)
    final_rows = rows(args.samtools, final)
    metadata = dict(records=count, kept=sum(int(row[4]) >= 30 for row in baseline),
                    marked=sum(bool(int(row[1]) & 1024) for row in final_rows),
                    ordered_sam_sha256=digest(baseline),
                    input=str(bam), sam=str(source), scope="small correctness only")
    write_json(root / "small.json", metadata)
    print(json.dumps(metadata), flush=True)


def prepare_full(args, root):
    deadline = time.monotonic() + 600
    while True:
        fields = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        available = int(fields["MemAvailable"].split()[0]) * 1024
        space = os.statvfs(str(root))
        if available >= 112 * 1024 ** 3 and space.f_bavail * space.f_frsize >= 200 * 1024 ** 3:
            break
        if time.monotonic() >= deadline:
            raise RuntimeError("resource gate not met within 10 minutes; no indefinite wait")
        time.sleep(15)
    kept = int(subprocess.check_output([args.samtools, "view", "-@", "8", "-c", "-q", "30",
                                      args.input]).strip())
    write_json(root / "full-filter.json", dict(records=args.records, kept=kept))
    final = four_steps(args, root, "full-reference", Path(args.input), 24)
    flagstat = root / "full.flagstat.json"
    command(args, root, "full.flagstat", [args.samtools, "flagstat", "-@", "8", "-O", "json", str(final)], flagstat)
    flags = json.loads(flagstat.read_text())
    total = sum(group["total"] for group in flags.values())
    marked = sum(group["duplicates"] for group in flags.values())
    if total != args.records:
        raise RuntimeError("full reference lost records")
    stderr = (root / "logs" / "full-reference.sort.stderr").read_text()
    if "merging from 0 files" not in stderr:
        raise RuntimeError("reference sort was not in-memory")
    write_json(root / "full.json", dict(records=total, kept=kept, marked=marked,
               threads=24, sort_memory="4096M per thread", full_reference=str(final),
               scope="one-off feasibility/correctness; not tuned or formal performance"))
    print("FULL_REFERENCE_COMPLETE", flush=True)


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def verify_small(args, root):
    base = rows(args.samtools, root / "small.bam")
    expected = json.loads((root / "small.json").read_text())
    out = root / "swbam-small"
    paths = [out / (name + suffix) for name, suffix in
             (("sdk", ".bam"), ("filter", ".bam"), ("collate", ".bam"), ("fixmate", ".bam"),
              ("sort", ".bam"), ("markdup", ".bam"), ("roundtrip", ".bam"), ("records", ".sam"))]
    command(args, root, "verify.small.quickcheck", [args.samtools, "quickcheck", "-v"] +
            [str(path) for path in paths[:-1]])
    results = {path.stem: rows(args.samtools, path) for path in paths}
    for name in ("sdk", "roundtrip", "records"):
        check(results[name] == base, name + " ordered fields/AUX mismatch")
    check(results["filter"] == [row for row in base if int(row[4]) >= 30], "filter mismatch")
    check(Counter(results["collate"]) == Counter(base), "collate record-content mismatch")
    seen = set()
    for name, _ in itertools.groupby(results["collate"], key=lambda row: row[0]):
        check(name not in seen, "non-contiguous collate QNAME " + name)
        seen.add(name)
    fixed_reference = rows(args.samtools, root / "small-reference" / "fixmate.bam")
    check(Counter(results["fixmate"]) == Counter(fixed_reference), "fixmate fields/AUX differ from samtools")
    check(Counter(results["sort"]) == Counter(results["fixmate"]), "sort changed records")
    header = subprocess.check_output([args.samtools, "view", "-H", str(root / "small.bam")],
                                     universal_newlines=True)
    references = [line.split("\t")[1][3:] for line in header.splitlines() if line.startswith("@SQ\t")]
    reference_ids = {name: i for i, name in enumerate(references)}
    for name in ("sort", "markdup"):
        keys = [(reference_ids.get(row[2], len(references)), int(row[3])) for row in results[name]]
        check(all(a <= b for a, b in zip(keys, keys[1:])), name + " coordinate order mismatch")
        check(len(results[name]) == expected["records"], name + " record count mismatch")
    marked = sum(bool(int(row[1]) & 1024) for row in results["markdup"])
    check(marked == expected["marked"], "markdup duplicate count differs from samtools")
    write_json(root / "small-verified.json", dict(records=len(base), marked=marked, passed=True,
               checks=["complete parse/header/EOF", "SDK and conversions ordered fields/AUX",
                       "filter", "collate content and QNAME continuity", "fixmate fields/AUX vs samtools",
                       "sort content/order", "markdup count/order/duplicate total vs samtools"]))
    print("SMALL_VERIFICATION_COMPLETE", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("small", "full", "verify"))
    parser.add_argument("root", type=Path)
    parser.add_argument("--input")
    parser.add_argument("--records", type=int, default=165373770)
    parser.add_argument("--small-records", type=int, default=200000)
    parser.add_argument("--samtools", default="/home/user_home/wzs/software/samtools-1.20/bin/samtools")
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    (args.root / "logs").mkdir(exist_ok=True)
    os.environ["LD_LIBRARY_PATH"] = "/home/user_home/wzs/software/libdeflate-1.20/lib64:" + os.environ.get("LD_LIBRARY_PATH", "")
    try:
        if args.mode == "small":
            prepare_small(args, args.root)
        elif args.mode == "full":
            prepare_full(args, args.root)
        else:
            verify_small(args, args.root)
    except Exception as error:
        write_json(args.root / (args.mode + ".failed.json"), dict(error=str(error)))
        raise


if __name__ == "__main__":
    main()
