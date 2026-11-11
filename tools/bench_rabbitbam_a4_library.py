#!/usr/bin/env python3
"""Full-record RabbitBAM library tuning; compatible with A4 Python 3.6."""

import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import statistics
import subprocess
import time

GIB = 1024 ** 3
DATASETS = {
    "D1": ("D1-WES-9G", 177798571,
           "b7bc9220bc5e7051b51b758e7018907644b8ba323d67dec82e63acd2c80bb046"),
    "D2": ("D2-WGS-25G", 365079492,
           "1cbe4e17224d89ef331aeaa655835a0ec467ac94cf89511ddb07e6e9c5d041a6"),
    "D3": ("D3-mouse-SRR7768905", 165373770,
           "472593083f7e98acb40a84ec917341adb9578a64a89bf5b0ecc3c1669a83dc8c"),
}
FIELDS = ["dataset", "mode", "phase", "repeat", "nr", "nw", "wall_s",
          "user_s", "sys_s", "maxrss_kb", "records", "batches", "output_bytes", "label"]


def screen_configs(mode):
    if mode == "count":
        return [(n, 0) for n in (1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 120)]
    result = set()
    for total in (16, 32, 64, 96, 120):
        for divisor in (5, 3, 2):
            nr = max(1, total // divisor)
            result.add((nr, total - nr))
    for total in (32, 64, 120):
        nr = total * 2 // 3
        result.add((nr, total - nr))
    return sorted(result, key=lambda value: (sum(value), value))


def refine_configs(mode, best):
    nr, nw = best
    if mode == "count":
        step = max(1, nr // 4)
        return sorted({(n, 0) for n in (nr - step, nr + step, nr - 2, nr + 2)
                       if 1 <= n <= 120} - {best})
    dr, dw = max(1, nr // 4), max(1, nw // 4)
    return sorted({(r, w) for r in (nr - dr, nr, nr + dr)
                   for w in (nw - dw, nw, nw + dw)
                   if r >= 1 and w >= 1 and r + w <= 120} - {best})


def rank_configs(rows, dataset, mode, phases):
    groups = {}
    for row in rows:
        if row["dataset"] == dataset and row["mode"] == mode and row["phase"] in phases:
            key = (int(row["nr"]), int(row["nw"]))
            groups.setdefault(key, []).append(float(row["wall_s"]))
    return sorted(groups, key=lambda key: (statistics.median(groups[key]), sum(key), key))


def parse_metrics(text, expected):
    matches = re.findall(r"^BENCH records=(\d+) batches=(\d+) nr=(\d+) nw=(\d+) ", text, re.M)
    if len(matches) != 1 or int(matches[0][0]) != expected:
        raise RuntimeError("missing or incorrect complete-record count")
    return int(matches[0][0]), int(matches[0][1])


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


class Benchmark:
    def __init__(self, args):
        self.args, self.root = args, Path(args.root).resolve()
        for name in ("logs", "out", "validation"):
            (self.root / name).mkdir(exist_ok=True)
        self.lock = (self.root.parent / "library.lock").open("a")
        fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        self.env = os.environ.copy()
        self.env.update(OMP_DYNAMIC="FALSE", OMP_WAIT_POLICY="PASSIVE", LC_ALL="C")
        self.env["LD_LIBRARY_PATH"] = args.libdeflate + ":" + self.env.get("LD_LIBRARY_PATH", "")
        self.binary = str(self.root / "bin/rabbitbam-library")
        self.digest_binary = str(self.root / "bin/sam-stream")
        plan = dict(datasets=args.datasets, inputs=DATASETS, samtools=args.samtools,
                    libdeflate=args.libdeflate, timing="GNU time process wall including final close",
                    level=1, cache="buffered NVMe; no cache dropping; input SHA checks before tuning",
                    environment={k: self.env[k] for k in ("OMP_DYNAMIC", "OMP_WAIT_POLICY", "LC_ALL")},
                    screen={m: screen_configs(m) for m in ("count", "copy")},
                    selection="screen, winner-neighborhood refinement, top-three confirmation medians",
                    confirmation_order="forward/reverse/forward, three runs per candidate",
                    formal="three independent samples of confirmed winner, first copy output fully validated",
                    upstream_commit="2c2bf3572fa76ae775787d324f6660f06d6f6519",
                    parser_threads=4, serializer_threads=4, max_codec_workers=120)
        plan_path = self.root / "plan.json"
        encoded = json.dumps(plan, indent=2, sort_keys=True) + "\n"
        if plan_path.exists() and plan_path.read_text() != encoded:
            raise RuntimeError("refusing to resume with a different benchmark plan")
        plan_path.write_text(encoded)
        self.rows = []
        self.results = self.root / "results.tsv"
        if self.results.exists():
            with self.results.open() as stream:
                self.rows = list(csv.DictReader(stream, delimiter="\t"))
        else:
            with self.results.open("w") as stream:
                csv.DictWriter(stream, FIELDS, delimiter="\t").writeheader()
        self.log("START pid={} root={}".format(os.getpid(), self.root))

    def log(self, message):
        print(time.strftime("%Y-%m-%dT%H:%M:%S%z") + " " + message, flush=True)
        (self.root / "state").write_text(message + "\n")

    def input(self, dataset):
        label = DATASETS[dataset][0]
        return Path(self.args.data) / label / "initial" / (label + ".coordinate.bam")

    def call(self, label, argv, timeout=3600):
        with (self.root / "commands.jsonl").open("a") as stream:
            stream.write(json.dumps(dict(label=label, argv=argv, time=time.time())) + "\n")
        with (self.root / "logs" / (label + ".stdout")).open("w") as out, \
                (self.root / "logs" / (label + ".stderr")).open("w") as err:
            process = subprocess.Popen(argv, stdout=out, stderr=err, env=self.env,
                                       start_new_session=True)
            try:
                code = process.wait(timeout=timeout)
            except BaseException:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                raise
        if code:
            raise RuntimeError("{} exited {}; see logs".format(label, code))
        return (self.root / "logs" / (label + ".stdout")).read_text()

    def wait_resources(self):
        end = time.monotonic() + 1800
        while True:
            memory = {}
            for line in Path("/proc/meminfo").read_text().splitlines():
                name, value = line.split(":", 1)
                memory[name] = int(value.strip().split()[0]) * 1024
            disk = shutil.disk_usage(str(self.root)).free
            if disk < 120 * GIB:
                raise RuntimeError("less than 120 GiB NVMe free space; refusing to continue")
            def cpu():
                values = list(map(int, Path("/proc/stat").read_text().splitlines()[0].split()[1:9]))
                return sum(values), values[3] + values[4]
            a, ai = cpu()
            time.sleep(0.5)
            b, bi = cpu()
            busy = 1 - (bi - ai) / max(1, b - a)
            if memory.get("MemAvailable", 0) >= 32 * GIB and busy < 0.20:
                return
            self.log("WAIT_RESOURCES available_gib={:.1f} cpu_busy={:.1%}".format(
                memory.get("MemAvailable", 0) / GIB, busy))
            if time.monotonic() > end:
                raise RuntimeError("resource wait exceeded 30 minutes")
            time.sleep(30)

    def preflight(self):
        self.call("freeze-check", ["sha256sum", "-c", str(self.root / "frozen.sha256")])
        for label, command in (("hostname", ["hostname"]), ("cpu", ["lscpu"]),
                               ("memory", ["free", "-b"]), ("disk", ["df", "-h", str(self.root)]),
                               ("samtools-version", [self.args.samtools, "--version"])):
            if not (self.root / "logs" / (label + ".stdout")).exists():
                self.call(label, command)
        for dataset in self.args.datasets:
            path = self.input(dataset)
            marker = self.root / "validation" / (dataset + ".input.json")
            info = dict(path=str(path), bytes=path.stat().st_size,
                        mtime_ns=path.stat().st_mtime_ns, expected_sha256=DATASETS[dataset][2])
            if marker.exists() and json.loads(marker.read_text()) == info:
                continue
            self.log("VERIFY_INPUT " + dataset)
            if sha256(path) != DATASETS[dataset][2]:
                raise RuntimeError(dataset + " input BAM does not match frozen paper input")
            self.call("quickcheck-input-" + dataset, [self.args.samtools, "quickcheck", "-v", str(path)])
            marker.write_text(json.dumps(info, indent=2) + "\n")
        self.small_regression()

    def small_regression(self):
        marker = self.root / "validation/small.json"
        if marker.exists():
            return
        results = []
        for dataset in self.args.datasets:
            sam = self.root / "validation" / (dataset + ".small.sam")
            bam = sam.with_suffix(".bam")
            # Stop at a complete record; do not use shell pipelines/SIGPIPE as success.
            process = subprocess.Popen([self.args.samtools, "view", "--no-PG", "-h", str(self.input(dataset))],
                                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, env=self.env)
            with sam.open("wb") as out:
                count = 0
                for line in process.stdout:
                    out.write(line)
                    if not line.startswith(b"@"):
                        count += 1
                        if count == 250003:
                            break
            process.stdout.close()
            process.terminate()
            process.wait()
            self.call("small-encode-" + dataset, [self.args.samtools, "view", "--no-PG", "-b", "-o", str(bam), str(sam)])
            for nr, nw in ((2, 0), (2, 3), (9, 7)):
                mode = "copy" if nw else "count"
                output = str(bam) + ".{}-{}.bam".format(nr, nw)
                label = "small-{}-{}-{}".format(dataset, nr, nw)
                text = self.call(label, [self.binary, mode, str(bam), output if nw else "-",
                                         str(nr), str(nw), str(count)])
                parse_metrics(text, count)
                if nw:
                    self.verify_pair(label, bam, Path(output), count)
                    Path(output).unlink()
                results.append(dict(dataset=dataset, nr=nr, nw=nw, records=count))
        marker.write_text(json.dumps(results, indent=2) + "\n")

    def verify_pair(self, label, source, output, expected):
        self.call(label + ".quickcheck", [self.args.samtools, "quickcheck", "-v", str(output)])
        for mode in ("digest-header", "digest-records"):
            baseline = self.root / "validation" / (source.stem + "." + mode + ".json")
            if baseline.exists():
                reference = json.loads(baseline.read_text())["digest"]
            else:
                reference = self.call(label + ".input." + mode, [self.digest_binary, mode, str(source)]).strip()
                baseline.write_text(json.dumps(dict(input=str(source), digest=reference)) + "\n")
            actual = self.call(label + ".output." + mode, [self.digest_binary, mode, str(output)]).strip()
            if reference != actual or not re.fullmatch(r"[0-9a-f]{32}", reference):
                raise RuntimeError(label + " failed " + mode)
            if mode == "digest-records":
                text = (self.root / "logs" / (label + ".output." + mode + ".stderr")).read_text()
                if "records={}\n".format(expected) not in text:
                    raise RuntimeError(label + " output record count mismatch")
        (self.root / "validation" / (label + ".json")).write_text(
            json.dumps(dict(input=str(source), output=str(output), records=expected,
                            ordered_record_digest_equal=True, header_without_PG_equal=True), indent=2) + "\n")

    def run_sample(self, dataset, mode, config, phase, repeat):
        nr, nw = config
        label = "{}.{}.{}.nr{}.nw{}.r{}".format(dataset, mode, phase, nr, nw, repeat)
        for row in self.rows:
            if row["label"] == label:
                self.log("REUSE " + label)
                return row
        self.wait_resources()
        output = self.root / "out" / (label + ".bam")
        if output.exists():
            raise RuntimeError("unaccepted output already exists: " + str(output))
        command = [self.binary, mode, str(self.input(dataset)), str(output) if nw else "-",
                   str(nr), str(nw), str(DATASETS[dataset][1])]
        timing = self.root / "logs" / (label + ".time")
        self.log("RUNNING " + label)
        text = self.call(label, ["/usr/bin/time", "-o", str(timing), "-f", "%e\t%U\t%S\t%M\t%x"] + command)
        records, batches = parse_metrics(text, DATASETS[dataset][1])
        values = timing.read_text().strip().splitlines()[-1].split("\t")
        if len(values) != 5 or int(values[-1]) != 0:
            raise RuntimeError("invalid GNU time result for " + label)
        if nw:
            self.call(label + ".quickcheck", [self.args.samtools, "quickcheck", "-v", str(output)])
        if nw and phase == "formal" and repeat == 1:
            self.log("FULL_OUTPUT_VALIDATION " + label)
            self.verify_pair(label, self.input(dataset), output, records)
        row = dict(zip(FIELDS, [dataset, mode, phase, repeat, nr, nw] + values[:4] +
                       [records, batches, output.stat().st_size if nw else 0, label]))
        with self.results.open("a") as stream:
            csv.DictWriter(stream, FIELDS, delimiter="\t").writerow(row)
            stream.flush()
            os.fsync(stream.fileno())
        self.rows.append(row)
        if nw and not (phase == "formal" and repeat == 1):
            output.unlink()
        self.log("ACCEPTED {} wall_s={}".format(label, values[0]))
        return row

    def tune(self, dataset, mode):
        for config in screen_configs(mode):
            self.run_sample(dataset, mode, config, "screen", 1)
        ranked = rank_configs(self.rows, dataset, mode, {"screen"})
        existing = set(screen_configs(mode))
        for config in refine_configs(mode, ranked[0]):
            if config not in existing:
                self.run_sample(dataset, mode, config, "refine", 1)
        finalists = rank_configs(self.rows, dataset, mode, {"screen", "refine"})[:3]
        for repeat in (1, 2, 3):
            order = finalists if repeat != 2 else list(reversed(finalists))
            for config in order:
                self.run_sample(dataset, mode, config, "confirm", repeat)
        selected = rank_configs(self.rows, dataset, mode, {"confirm"})[0]
        self.log("SELECTED {} {} nr={} nw={}".format(dataset, mode, *selected))
        for repeat in (1, 2, 3):
            self.run_sample(dataset, mode, selected, "formal", repeat)
        self.summary()

    def summary(self):
        result = []
        for dataset in self.args.datasets:
            for mode in ("count", "copy"):
                rows = [r for r in self.rows if r["dataset"] == dataset and r["mode"] == mode
                        and r["phase"] == "formal"]
                if not rows:
                    continue
                values = [float(r["wall_s"]) for r in rows]
                result.append(dict(dataset=dataset, mode=mode, nr=int(rows[0]["nr"]), nw=int(rows[0]["nw"]),
                                   samples_s=values, median_s=statistics.median(values), complete=len(rows) == 3,
                                   throughput_m_records_s=DATASETS[dataset][1] / statistics.median(values) / 1e6))
        (self.root / "summary.json").write_text(json.dumps(result, indent=2) + "\n")

    def run(self):
        self.preflight()
        if self.args.preflight_only:
            self.log("PREFLIGHT_COMPLETE")
            return
        for dataset in self.args.datasets:
            for mode in ("count", "copy"):
                self.tune(dataset, mode)
        self.call("freeze-check-final", ["sha256sum", "-c", str(self.root / "frozen.sha256")])
        self.log("COMPLETE")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--data", default="/home/bigssd/wzs_data")
    parser.add_argument("--datasets", nargs="+", choices=sorted(DATASETS), default=["D1", "D2", "D3"])
    parser.add_argument("--samtools", default="/home/user_home/wzs/software/samtools-1.20/bin/samtools")
    parser.add_argument("--libdeflate", default="/home/user_home/wzs/software/libdeflate-1.20/lib64")
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    benchmark = Benchmark(args)
    try:
        benchmark.run()
    except BaseException as error:
        benchmark.log("FAILED " + repr(error))
        raise


if __name__ == "__main__":
    main()
