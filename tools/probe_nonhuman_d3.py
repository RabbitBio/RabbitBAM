#!/usr/bin/env python3
"""Submit only the new mouse D3 capacity/correctness gates, not formal medians."""

import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import time


SW_DATA = "/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data"
A4_DATA = "/home/bigssd/wzs_data"
DATASET = "D3-mouse-SRR7768905"
INPUT = DATASET + ".coordinate.bam"
FROZEN = SW_DATA + "/swbam_bench/ram_D2_processing_matrix_20261013_114500/bin"
BINARIES = {
    "RabbitBAM-MPI": "81b68788a9d897b80f494fd1d84945f15b2fe8c7791a28f3559779c77c683cb1",
    "RabbitBAM-MPI-markdup": "3dc30d0a9ac3ba2a384fec617b66f32d7801a53d389288db0f51f110c7726a51",
    "swbam-sdk-optimized-count": "abdca8692b6b3fa6cb888da4da8ade0c7a268403628807e095a5a5c83a67a7f7",
    "swbam-sdk-optimized-bam2bam": "121d764161e8f6737aad9ada8f3ce5c85195f46456fa3daaba8f682b404df9e4",
}
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
       "-o", "ServerAliveInterval=20", "-o", "ServerAliveCountMax=3"]


def remote(host, argv, timeout=120, script=None):
    result = subprocess.run(SSH + [host, " ".join(shlex.quote(str(arg)) for arg in argv)],
                            input=script, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError("remote command failed: {}\n{}\n{}".format(host, result.stdout, result.stderr))
    return result.stdout


def copy(source, destination):
    subprocess.run(["rsync", "-a", "--partial", "-e", " ".join(SSH), source, destination],
                   check=True, timeout=1800)


def rank_lines(log, label, ranks):
    rows = re.findall(r"\[" + re.escape(label) + r" rank=(\d+)\]([^\n]*)", log)
    if len(rows) != ranks or {int(rank) for rank, _ in rows} != set(range(ranks)):
        raise RuntimeError("missing or duplicate " + label + " rank logs")
    return [dict(re.findall(r"(\w+)=([^\s]+)", rest)) for _, rest in rows]


def metrics(log, ranks, records, check="any", kept=None, marked=None):
    if re.search(r"ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT|exit_code=\s*[1-9]", log):
        raise RuntimeError("application error in log")
    e2e, processing, storage = [rank_lines(log, name, ranks) for name in
                               ("swbam-e2e", "swbam-processing", "ram-storage")]
    if any(row.get("status") != "0" for row in e2e + processing):
        raise RuntimeError("nonzero rank status")
    if not re.search(r"(?:total_records|records)=" + str(records) + r"(?:\s|$)", log):
        raise RuntimeError("wrong global record count")
    if check == "memory" and "selected mode=memory" not in log:
        raise RuntimeError("sort/collate not in-memory")
    if check == "filter" and not re.search(r"kept_records=" + str(kept) + r"(?:\s|$)", log):
        raise RuntimeError("filter differs from samtools")
    if check == "markdup" and not re.search(r"\bmarked[=: ]+" + str(marked) + r"(?:\s|$)", log):
        raise RuntimeError("markdup differs from samtools duplicate total")
    if check == "sorted" and not re.search(r"^SN\tis sorted:\t1(?:\s|$)", log, re.M):
        raise RuntimeError("full reread reports unsorted output")
    if check != "markdup" and any(int(row.get(field, 0)) for row in storage
                                  for field in ("temp_read_bytes", "temp_write_bytes")):
        raise RuntimeError("unexpected RAM temp traffic in in-memory probe")
    by_rank = {}
    for label in ("swbam-e2e", "swbam-processing"):
        by_rank[label] = {int(rank): float(value) for rank, value in re.findall(
            r"\[" + label + r" rank=(\d+)\].*?seconds=([0-9.]+)", log)}
    if any(by_rank["swbam-processing"][rank] <= 0 or
           by_rank["swbam-processing"][rank] > by_rank["swbam-e2e"][rank] + .000002
           for rank in range(ranks)):
        raise RuntimeError("invalid processing/MPI interval")
    prepared = re.findall(r"\[ram-prepare rank=\d+\].*?seconds=([0-9.]+)", log)
    return dict(processing_max_s=max(float(row["seconds"]) for row in processing),
                e2e_max_s=max(float(row["seconds"]) for row in e2e),
                prepare_max_s=max(map(float, prepared), default=0),
                pool_peak_max_bytes=max(int(row["peak"]) for row in storage),
                output_bytes=max(int(row.get("output_bytes", 0)) for row in storage),
                temp_read_bytes=sum(int(row.get("temp_read_bytes", 0)) for row in storage),
                temp_write_bytes=sum(int(row.get("temp_write_bytes", 0)) for row in storage))


def completed_log(fetch, ranks, expected, check="any", attempts=40, pause=time.sleep):
    # Scheduler DONE can precede the shared filesystem's complete log visibility.
    for attempt in range(attempts):
        log = fetch()
        try:
            result = metrics(log, ranks, expected["records"], check,
                             expected.get("kept"), expected.get("marked"))
            return log, result
        except RuntimeError:
            if re.search(r"ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT|exit_code=\s*[1-9]", log):
                raise
            if attempt + 1 == attempts:
                raise
            pause(3)


class Probe:
    def __init__(self, args):
        self.args = args
        self.root = args.root.resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.lock = (self.root / "driver.lock").open("a")
        fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        self.sw = SW_DATA + "/swbam_bench/" + self.root.name
        self.a4 = A4_DATA + "/samtools_bench/" + self.root.name
        self.records = 165373770
        self.current_job = None
        self.fixmate_tag = ""
        self.state("INITIALIZING")

    def state(self, phase, **values):
        value = dict(phase=phase, local_time=time.strftime("%Y-%m-%d %H:%M:%S %z"),
                     swls_root=self.sw, a4_root=self.a4, **values)
        temporary = self.root / "state.tmp"
        temporary.write_text(json.dumps(value, indent=2) + "\n")
        temporary.replace(self.root / "state.json")
        print(json.dumps(value), flush=True)

    def reference(self, mode, timeout=300):
        return remote(self.args.a4, ["python3", self.a4 + "/d3_probe_reference.py", mode,
                      self.a4, "--input", A4_DATA + "/" + DATASET + "/initial/" + INPUT], timeout)

    def initialize(self):
        if self.args.fixmate_app:
            if not self.args.fixmate_app.startswith("/"):
                raise RuntimeError("fixmate app must be an absolute frozen path")
            actual = remote(self.args.swls, ["sha256sum", self.args.fixmate_app])
            manifest = self.root / "fixmate.sha256"
            if manifest.exists() and manifest.read_text() != actual:
                raise RuntimeError("fixmate override fingerprint changed")
            manifest.write_text(actual)
            self.fixmate_tag = "fm_" + actual.split()[0][:12]
        if self.args.resume:
            actual = remote(self.args.swls, ["sha256sum"] + [self.sw + "/bin/" + name for name in BINARIES])
            old = {Path(path).name: checksum for checksum, path in
                   (line.split() for line in (self.root / "binaries.sha256").read_text().splitlines())}
            if {Path(path).name: checksum for checksum, path in
                (line.split() for line in actual.splitlines())} != old:
                raise RuntimeError("frozen binaries changed since initial submission")
            tools = Path(__file__).resolve().parent
            revision = "".join(hashlib.sha256((tools / name).read_bytes()).hexdigest() + "  " + name + "\n"
                               for name in ("probe_nonhuman_d3.py", "d3_probe_reference.py"))
            with (self.root / "script-revisions.log").open("a") as log:
                log.write(time.strftime("%Y-%m-%d %H:%M:%S %z") + " explicit resume\n" + revision)
            return
        if (self.root / "jobs.jsonl").exists():
            raise RuntimeError("use a new root; automatic resubmission of old jobs is prohibited")
        remote(self.args.swls, ["bash", "-s"], script="\n".join([
            "set -eu", "test -f " + shlex.quote(SW_DATA + "/" + DATASET + "/initial/.transfer_verified"),
            "test ! -e " + shlex.quote(self.sw),
            "mkdir -p " + " ".join(shlex.quote(self.sw + "/" + d) for d in ("bin", "logs", "out", "sessions", "small")),
        ] + ["cp -p " + shlex.quote(FROZEN + "/" + name) + " " + shlex.quote(self.sw + "/bin/" + name)
             for name in BINARIES]))
        remote(self.args.a4, ["bash", "-s"], script="\n".join([
            "set -eu", "test -f " + shlex.quote(A4_DATA + "/" + DATASET + "/initial/.transfer_verified"),
            "test ! -e " + shlex.quote(self.a4), "mkdir -p " + shlex.quote(self.a4)]))
        actual = remote(self.args.swls, ["sha256sum"] + [self.sw + "/bin/" + name for name in BINARIES])
        (self.root / "binaries.sha256").write_text(actual)
        for line in actual.splitlines():
            checksum, path = line.split()
            if BINARIES[Path(path).name] != checksum:
                raise RuntimeError("frozen binary checksum mismatch: " + path)
        tools = Path(__file__).resolve().parent
        copy(str(tools / "d3_probe_reference.py"), self.args.a4 + ":" + self.a4 + "/")
        (self.root / "scripts.sha256").write_text("".join(
            hashlib.sha256((tools / name).read_bytes()).hexdigest() + "  " + name + "\n"
            for name in ("probe_nonhuman_d3.py", "d3_probe_reference.py")))
        self.state("PREPARING_SMALL_REFERENCE")
        print(self.reference("small"), flush=True)
        copy(self.args.a4 + ":" + self.a4 + "/small.*", str(self.root) + "/")
        for filename in ("small.bam", "small.sam"):
            copy(str(self.root / filename), self.args.swls + ":" + self.sw + "/small/")
        full_command = ["python3", self.a4 + "/d3_probe_reference.py", "full", self.a4,
                        "--input", A4_DATA + "/" + DATASET + "/initial/" + INPUT]
        remote(self.args.a4, ["bash", "-s"], script="\n".join([
            "set -eu", "nohup " + " ".join(map(shlex.quote, full_command)) +
            " > " + shlex.quote(self.a4 + "/full.driver.log") + " 2>&1 < /dev/null &", "echo $!"
        ]))

    def expected(self, full=False):
        name = "full.json" if full else "small.json"
        if not full:
            return json.loads((self.root / name).read_text())
        deadline = time.monotonic() + 3600
        self.state("WAITING_A4_FULL_REFERENCE")
        while time.monotonic() < deadline:
            value = remote(self.args.a4, ["bash", "-s"], script="\n".join([
                "if test -f " + shlex.quote(self.a4 + "/full.failed.json") + "; then",
                "cat " + shlex.quote(self.a4 + "/full.failed.json"), "elif test -f " + shlex.quote(self.a4 + "/full.json") + "; then",
                "cat " + shlex.quote(self.a4 + "/full.json"), "fi"]))
            if value.strip():
                result = json.loads(value)
                if "error" in result:
                    raise RuntimeError("a4 full reference failed: " + result["error"])
                (self.root / name).write_text(value)
                return result
            time.sleep(20)
        raise RuntimeError("a4 reference not ready after one hour")

    def run(self, label, capacity, argv, expected, check="any", ranks=6):
        session = self.sw + "/sessions/" + label
        log_path = self.sw + "/logs/" + label + ".log"
        command = ["bsub", "-N", str(ranks), "-np", "1", "-health", "1", "-mpecg", "6", "-cgsp", "64",
                   "-share_size", "100", "-cross_size", "88000", "-cache_size", "32", "-xmalloc",
                   "-a", "mpi", "-b", "-K", "-q", "q_share", "-timelimit", "03:00:00",
                   "-J", "d3mouse_" + label, "-o", log_path] + argv + [
                   "--storage", "ram", "--ram-capacity", capacity, "--ram-ranks", str(ranks),
                   "--ram-session", session, "--ram-ready-timeout", "5400"]
        result_path = self.root / "results.tsv"
        if result_path.exists():
            with result_path.open() as stream:
                previous = next((row for row in csv.DictReader(stream, delimiter="\t") if row["label"] == label), None)
            if previous:
                log = (self.root / "logs" / (label + ".log")).read_text()
                metrics(log, ranks, expected["records"], check, expected.get("kept"), expected.get("marked"))
                self.state("REUSING_VERIFIED_GATE", label=label, job=previous["job"])
                return
        previous = None
        jobs_path = self.root / "jobs.jsonl"
        if jobs_path.exists():
            previous = next((row for row in reversed(list(map(json.loads, jobs_path.read_text().splitlines())))
                             if row["label"] == label), None)
        if previous and self.args.retry_no_start:
            status = remote(self.args.swls, ["bash", "-s"], script="bjobs " + previous["job"] + " 2>/dev/null || true\n")
            if re.search(r"^" + previous["job"] + r"\s+EXIT\s", status, re.M):
                remote(self.args.swls, ["bash", "-s"], script="\n".join([
                    "set -eu", "test ! -e " + shlex.quote(session),
                    "if test -f " + shlex.quote(previous["log"]) + "; then",
                    "grep -q ': not found' " + shlex.quote(previous["log"]),
                    "! grep -q '\\[ram-prepare rank=' " + shlex.quote(previous["log"]),
                    "mkdir -p " + shlex.quote(self.sw + "/retries/" + previous["job"]),
                    "mv " + shlex.quote(previous["log"]) + " " + shlex.quote(self.sw + "/retries/" + previous["job"] + "/"),
                    "fi"]))
                details = remote(self.args.swls, ["bjobs", "-l", previous["job"]])
                (self.root / ("failed-launch-" + previous["job"] + ".txt")).write_text(details)
                with (self.root / "failed-launches.jsonl").open("a") as stream:
                    stream.write(json.dumps(previous) + "\n")
                previous = None
        self.state("SUBMITTING", label=label, ranks=ranks, capacity=capacity)
        if previous:
            if previous["argv"] != command:
                raise RuntimeError("refusing to change an existing job's configuration")
            job = previous["job"]
        else:
            submission = remote(self.args.swls, ["bash", "-s"], script="\n".join([
                "set -eu", "cd " + shlex.quote(self.sw),
                "env SWBAM_NODE_MEMORY=88000M " + " ".join(map(shlex.quote, command))]))
            match = re.search(r"Job <(\d+)>", submission)
            if not match:
                raise RuntimeError("no submitted job id: " + submission)
            job = match.group(1)
            with (self.root / "jobs.jsonl").open("a") as stream:
                stream.write(json.dumps(dict(label=label, job=job, argv=command, log=log_path)) + "\n")
        self.current_job = job
        self.state("WAITING_JOB", label=label, job=job, ranks=ranks, capacity=capacity)
        deadline = time.monotonic() + 6 * 3600
        failures = 0
        unknown = 0

        def fetch_log():
            return remote(self.args.swls, ["bash", "-s"], script="test ! -f " + shlex.quote(log_path) +
                          " || cat " + shlex.quote(log_path) + "\n")

        while time.monotonic() < deadline:
            try:
                status = remote(self.args.swls, ["bash", "-s"], script="bjobs " + job + " 2>/dev/null || true\n")
                match = re.search(r"^" + job + r"\s+(\S+)", status, re.M)
                state = match.group(1) if match else ""
                failures = 0
                if state not in ("PEND", "RUN", "EXITING"):
                    log = fetch_log()
                    if state == "DONE" or (not state and "Job " + job + " has been finished" in log):
                        break
                    if not state:
                        # On explicit resume, old jobs may have left scheduler history.
                        # Reuse only a complete, strictly validated original log.
                        if previous and self.args.resume:
                            try:
                                metrics(log, ranks, expected["records"], check,
                                        expected.get("kept"), expected.get("marked"))
                                break
                            except RuntimeError:
                                pass
                        unknown += 1
                        if unknown >= 40:
                            raise RuntimeError("job absent from scheduler with no validated completion: " + job)
                    if state:
                        (self.root / "logs").mkdir(exist_ok=True)
                        (self.root / "logs" / (label + ".failed.log")).write_text(log)
                        raise RuntimeError("scheduler state " + state + " for " + job)
            except (subprocess.TimeoutExpired, OSError):
                failures += 1
                if failures >= 5:
                    raise
            time.sleep(10)
        else:
            remote(self.args.swls, ["bkill", job])
            raise RuntimeError("probe's own job exceeded 6-hour submission deadline")
        log, result = completed_log(fetch_log, ranks, expected, check)
        (self.root / "logs").mkdir(exist_ok=True)
        (self.root / "logs" / (label + ".log")).write_text(log)
        result.update(label=label, job=job, ranks=ranks, capacity=capacity)
        path = self.root / "results.tsv"
        with path.open("a", newline="") as stream:
            writer = csv.DictWriter(stream, list(result), delimiter="\t")
            if stream.tell() == 0:
                writer.writeheader()
            writer.writerow(result)
        self.current_job = None
        self.state("GATE_PASSED", **result)

    def suite(self, small, expected):
        phase = "small" if small else "full"
        source = self.sw + "/small/small.bam" if small else SW_DATA + "/" + DATASET + "/initial/" + INPUT
        out = self.sw + "/out/" + phase
        remote(self.args.swls, ["mkdir", "-p", out])
        fixed_out = out + "-" + self.fixmate_tag if self.fixmate_tag else out
        remote(self.args.swls, ["mkdir", "-p", fixed_out])
        app = self.sw + "/bin/RabbitBAM-MPI"
        fixmate_app = self.args.fixmate_app or app
        markdup = self.sw + "/bin/RabbitBAM-MPI-markdup"
        count = self.sw + "/bin/swbam-sdk-optimized-count"
        copier = self.sw + "/bin/swbam-sdk-optimized-bam2bam"
        read_pool, write_pool, fixed_pool, sam_pool = ("1G", "2G", "2G", "2G") if small else ("8G", "16G", "24G", "40G")
        algorithm = "1G" if small else "48G"

        def run(name, pool, argv, check="any"):
            suffix = "_" + self.fixmate_tag if self.fixmate_tag and name in (
                "fixmate", "sort", "validate_sort", "markdup", "validate_markdup") else ""
            self.run(phase + "_" + name + suffix, pool, argv, expected, check)

        def validate(name, path, sorted_output=False):
            run("validate_" + name, read_pool, [app, "stats", "--basic", "-i", path],
                "sorted" if sorted_output else "any")

        run("count", read_pool, [count, source])
        run("sdk", write_pool, [copier, source, "-o", out + "/sdk.bam"])
        validate("sdk", out + "/sdk.bam", True)
        run("flagstat", read_pool, [app, "flagstat", "-i", source])
        run("filter", write_pool, [app, "run_all", "-i", source, "-o", out + "/filter.bam",
                                   "--min-mapq", "30", "--compress-level", "1"], "filter")
        run("collate", write_pool, [app, "collate", "-i", source, "-o", out + "/collate.bam",
                                    "-n", "66", "-m", algorithm, "--compress-level", "1"], "memory")
        run("fixmate", fixed_pool, [fixmate_app, "fixmate", "-i", out + "/collate.bam", "-o", fixed_out + "/fixmate.bam",
                                    "-m", "--compress-level", "1"])
        run("sort", write_pool, [app, "sort", "-i", fixed_out + "/fixmate.bam", "-o", fixed_out + "/sort.bam",
                                 "-m", algorithm, "--compress-level", "1"], "memory")
        validate("sort", fixed_out + "/sort.bam", True)
        run("markdup", fixed_pool, [markdup, "markdup", "-i", fixed_out + "/sort.bam", "-o", fixed_out + "/markdup.bam",
                                    "--stream", "-l", "300", "-m", "1G" if small else "8G",
                                    "--compress-level", "1"], "markdup")
        validate("markdup", fixed_out + "/markdup.bam", True)
        run("bam2sam", sam_pool, [app, "run_all", "-i", source, "-o", out + "/records.sam"])
        run("sam2bam", sam_pool, [app, "run_all", "-i", out + "/records.sam", "-o", out + "/roundtrip.bam",
                                  "--compress-level", "1"])
        validate("roundtrip", out + "/roundtrip.bam", True)
        return fixed_out

    def execute(self):
        self.initialize()
        fixed_out = self.suite(True, self.expected())
        destination = self.root / "small-output"
        destination.mkdir(exist_ok=True)
        copy(self.args.swls + ":" + self.sw + "/out/small/", str(destination) + "/")
        if self.fixmate_tag:
            copy(self.args.swls + ":" + fixed_out + "/", str(destination) + "/")
        copy(str(destination) + "/", self.args.a4 + ":" + self.a4 + "/swbam-small/")
        self.state("VERIFYING_SMALL_FIELDS_ON_A4")
        print(self.reference("verify"), flush=True)
        copy(self.args.a4 + ":" + self.a4 + "/small-verified.json", str(self.root) + "/")
        expected = self.expected(True)
        self.suite(False, expected)
        source = SW_DATA + "/" + DATASET + "/initial/" + INPUT
        self.run("full_count_n1", "16G", [self.sw + "/bin/swbam-sdk-optimized-count", source], expected, ranks=1)
        self.run("full_sdk_n1", "48G", [self.sw + "/bin/swbam-sdk-optimized-bam2bam", source,
                 "-o", self.sw + "/out/sdk.n1.bam", "--no-dump"], expected, ranks=1)
        self.state("COMPLETE", note="capacity/correctness only; no formal matrix submitted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--a4", default="wzs@a4.rabbitbio.org")
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--resume", action="store_true", help="reuse only frozen, already-verified gates; never resubmit failed jobs")
    parser.add_argument("--retry-no-start", action="store_true",
                        help="explicitly retry only terminal EXIT jobs with neither log nor RAM session")
    parser.add_argument("--fixmate-app", help="separately frozen fixmate-only binary; dependent gates use versioned labels and output paths")
    args = parser.parse_args()
    probe = Probe(args)
    try:
        probe.execute()
    except Exception as error:
        probe.state("FAILED", error=str(error), job=probe.current_job)
        raise


if __name__ == "__main__":
    main()
