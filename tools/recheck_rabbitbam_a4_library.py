#!/usr/bin/env python3
"""Targeted, non-destructive stability checks of frozen RabbitBAM benchmarks."""
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import signal
import statistics
import subprocess
import time

import bench_rabbitbam_a4_library as base

FIELDS = ["dataset", "mode", "phase", "repeat", "policy", "nr", "nw", "wall_s",
          "user_s", "sys_s", "maxrss_kb", "input_blocks", "output_blocks",
          "voluntary_switches", "involuntary_switches", "major_faults", "minor_faults",
          "records", "batches", "label"]
TIME_FORMAT = "%e\t%U\t%S\t%M\t%I\t%O\t%w\t%c\t%F\t%R\t%x"


def policy_command(policy, nr, nw):
    if policy in ("baseline", "active"):
        return []
    if policy == "interleave":
        return ["numactl", "--interleave=all"]
    if policy in ("node0", "node1"):
        if nr + nw + 7 > 64:
            return None
        node = policy[-1]
        return ["numactl", "--cpunodebind=" + node, "--membind=" + node]
    raise ValueError(policy)


def metrics(text):
    values = text.strip().splitlines()[-1].split("\t")
    if len(values) != 11 or int(values[-1]):
        raise RuntimeError("invalid GNU time result")
    return dict(zip(FIELDS[7:17], values[:10]))


def candidate_configs(rows, dataset, mode):
    configs = base.rank_configs(rows, dataset, mode, {"confirm"})[:3]
    if not configs:
        raise RuntimeError("missing original confirmation samples")
    if not any(sum(config) + 7 <= 64 for config in configs):
        eligible = [row for row in rows if int(row["nr"]) + int(row["nw"]) + 7 <= 64]
        local = base.rank_configs(eligible, dataset, mode, {"screen", "refine", "confirm"})
        if local and local[0] not in configs:
            configs.append(local[0])
    return configs


def summarize(rows):
    groups = {}
    for row in rows:
        key = tuple(row[name] for name in ("dataset", "mode", "phase", "policy")) + (
            int(row["nr"]), int(row["nw"]))
        groups.setdefault(key, []).append(row)
    result = []
    for key, values in sorted(groups.items()):
        times = [float(row["wall_s"]) for row in values]
        median = statistics.median(times)
        result.append(dict(zip(("dataset", "mode", "phase", "policy", "nr", "nw"), key),
                           samples_s=times, median_s=median, minimum_s=min(times), maximum_s=max(times),
                           range_over_median_pct=100 * (max(times) - min(times)) / median,
                           median_cpu_s=statistics.median(float(v["user_s"]) + float(v["sys_s"]) for v in values),
                           median_input_blocks=statistics.median(int(v["input_blocks"]) for v in values),
                           median_voluntary_switches=statistics.median(int(v["voluntary_switches"]) for v in values),
                           median_involuntary_switches=statistics.median(int(v["involuntary_switches"]) for v in values)))
    return result


class Recheck:
    def __init__(self, args):
        self.args = args
        self.root, self.source = Path(args.root).resolve(), Path(args.source_root).resolve()
        self.root.mkdir(exist_ok=True)
        for name in ("logs", "out", "validation", "snapshots"):
            (self.root / name).mkdir(exist_ok=True)
        self.lock = (self.source.parent / "library.lock").open("a")
        fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        self.binary = str(self.source / "bin/rabbitbam-library")
        self.env = os.environ.copy()
        self.env.update(OMP_DYNAMIC="FALSE", OMP_WAIT_POLICY="PASSIVE", LC_ALL="C",
                        LD_LIBRARY_PATH=args.libdeflate + ":" + self.env.get("LD_LIBRARY_PATH", ""))
        with (self.source / "results.tsv").open() as stream:
            self.original = list(csv.DictReader(stream, delimiter="\t"))
        self.rows = []
        path = self.root / "results.tsv"
        if path.exists():
            with path.open() as stream:
                self.rows = list(csv.DictReader(stream, delimiter="\t"))
        else:
            with path.open("w") as stream:
                csv.DictWriter(stream, FIELDS, delimiter="\t").writeheader()
        self.plan = dict(source_root=str(self.source), datasets=args.datasets, policies=args.policies,
                         rounds=args.rounds, formal_rounds=5, stage=args.stage,
                         frozen_binary_sha256=base.sha256(Path(self.binary)),
                         source_input_checks="original validated SHA plus unchanged size and mtime",
                         cache="one complete buffered input read per dataset/mode before samples",
                         timing="full GNU time process wall; no fsync; diagnostics kept separate")
        plan = self.root / (args.stage + "-plan.json")
        encoded = json.dumps(self.plan, sort_keys=True, indent=2) + "\n"
        if plan.exists() and plan.read_text() != encoded:
            raise RuntimeError("plan mismatch; use a new directory")
        plan.write_text(encoded)

    def log(self, text):
        print(time.strftime("%Y-%m-%dT%H:%M:%S%z") + " " + text, flush=True)
        (self.root / "state").write_text(text + "\n")

    def input(self, dataset):
        label = base.DATASETS[dataset][0]
        return Path(self.args.data) / label / "initial" / (label + ".coordinate.bam")

    def run(self, label, argv, env=None):
        with (self.root / "commands.jsonl").open("a") as stream:
            stream.write(json.dumps(dict(label=label, argv=argv,
                                        OMP_WAIT_POLICY=(env or self.env)["OMP_WAIT_POLICY"])) + "\n")
        with (self.root / "logs" / (label + ".stdout")).open("w") as out, \
                (self.root / "logs" / (label + ".stderr")).open("w") as err:
            process = subprocess.Popen(argv, stdout=out, stderr=err, env=env or self.env,
                                       start_new_session=True)
            try:
                code = process.wait(timeout=1800)
            except BaseException:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                raise
        if code:
            raise RuntimeError("{} exit={}".format(label, code))
        return (self.root / "logs" / (label + ".stdout")).read_text()

    def snapshot(self, label):
        result = dict(time=time.time(), load=os.getloadavg())
        for name, path in (("meminfo", "/proc/meminfo"), ("stat", "/proc/stat"),
                           ("cpuinfo", "/proc/cpuinfo"), ("vmstat", "/proc/vmstat")):
            result[name] = Path(path).read_text()
        (self.root / "snapshots" / (label + ".json")).write_text(json.dumps(result))

    def preflight(self):
        self.run("recheck-freeze", ["sha256sum", "-c", str(self.root / "scripts.sha256")])
        self.run("source-freeze", ["sha256sum", "-c", str(self.source / "frozen.sha256")])
        self.run("source-clean", ["git", "-C", str(self.source / "upstream"), "diff", "--exit-code"])
        for dataset in self.args.datasets:
            path = self.input(dataset)
            expected = json.loads((self.source / "validation" / (dataset + ".input.json")).read_text())
            info = path.stat()
            if expected["bytes"] != info.st_size or expected["mtime_ns"] != info.st_mtime_ns:
                raise RuntimeError("input changed: " + dataset)
        self.snapshot("start")

    def warm(self, dataset):
        self.log("BUFFERED_INPUT_PREPARE " + dataset)
        with self.input(dataset).open("rb") as stream:
            while stream.read(8 * 1024 * 1024):
                pass

    def sample(self, dataset, mode, config, policy, phase, repeat):
        nr, nw = config
        prefix = policy_command(policy, nr, nw)
        if prefix is None:
            return
        label = "{}.{}.{}.{}.nr{}.nw{}.r{}".format(dataset, mode, phase, policy, nr, nw, repeat)
        for row in self.rows:
            if row["label"] == label:
                return
        # Reuse the original guard without sharing its output or timing state.
        guard = base.Benchmark.__new__(base.Benchmark)
        guard.root, guard.log = self.root, self.log
        guard.wait_resources()
        output = self.root / "out" / (label + ".bam")
        if output.exists():
            raise RuntimeError("unaccepted output exists: " + str(output))
        env = self.env.copy()
        if policy == "active":
            env["OMP_WAIT_POLICY"] = "ACTIVE"
        timing = self.root / "logs" / (label + ".time")
        argv = ["/usr/bin/time", "-o", str(timing), "-f", TIME_FORMAT] + prefix + [
            self.binary, mode, str(self.input(dataset)), str(output) if nw else "-",
            str(nr), str(nw), str(base.DATASETS[dataset][1])]
        self.snapshot(label + ".before")
        self.log("RUNNING " + label)
        text = self.run(label, argv, env)
        self.snapshot(label + ".after")
        count, batches = base.parse_metrics(text, base.DATASETS[dataset][1])
        row = dict(dataset=dataset, mode=mode, phase=phase, repeat=repeat, policy=policy,
                   nr=nr, nw=nw, records=count, batches=batches, label=label,
                   **metrics(timing.read_text()))
        if nw:
            self.run(label + ".quickcheck", [self.args.samtools, "quickcheck", "-v", str(output)])
        with (self.root / "results.tsv").open("a") as stream:
            csv.DictWriter(stream, FIELDS, delimiter="\t").writerow(row)
            stream.flush()
            os.fsync(stream.fileno())
        self.rows.append(row)
        if nw and not (phase == "formal" and repeat == 1):
            output.unlink()
        self.write_summary()
        self.log("ACCEPTED {} wall_s={} input_blocks={} voluntary_switches={}".format(
            label, row["wall_s"], row["input_blocks"], row["voluntary_switches"]))

    def write_summary(self):
        (self.root / "summary.json").write_text(json.dumps(summarize(self.rows), indent=2) + "\n")

    def verify(self, row):
        dataset, label = row["dataset"], row["label"]
        output = self.root / "out" / (label + ".bam")
        self.log("FULL_OUTPUT_VALIDATION " + label)
        result = dict(dataset=dataset, output=str(output))
        for mode in ("digest-header", "digest-records"):
            expected = json.loads((self.source / "validation" /
                                   (self.input(dataset).stem + "." + mode + ".json")).read_text())
            if expected["input"] != str(self.input(dataset)):
                raise RuntimeError("unexpected reference input")
            actual = self.run(label + "." + mode,
                              [str(self.source / "bin/sam-stream"), mode, str(output)]).strip()
            if actual != expected["digest"]:
                raise RuntimeError("full output mismatch " + mode)
            result[mode] = actual
            if mode == "digest-records":
                text = (self.root / "logs" / (label + "." + mode + ".stderr")).read_text()
                if "records={}\n".format(base.DATASETS[dataset][1]) not in text:
                    raise RuntimeError("full output record count mismatch")
        (self.root / "validation" / (label + ".json")).write_text(json.dumps(result, indent=2) + "\n")

    def execute(self):
        self.preflight()
        if self.args.stage == "diagnostic":
            self.warm("D1")
            for repeat in range(1, self.args.rounds + 1):
                policies = self.args.policies if repeat % 2 else list(reversed(self.args.policies))
                for policy in policies:
                    self.sample("D1", "count", (8, 0), policy, "diagnostic", repeat)
        else:
            for dataset in self.args.datasets:
                for mode in ("count", "copy"):
                    self.warm(dataset)
                    configs = candidate_configs(self.original, dataset, mode)
                    choices = [(config, policy) for config in configs for policy in self.args.policies
                               if policy_command(policy, *config) is not None]
                    for repeat in range(1, self.args.rounds + 1):
                        order = choices if repeat % 2 else list(reversed(choices))
                        for config, policy in order:
                            self.sample(dataset, mode, config, policy, "confirm", repeat)
                    candidates = [v for v in summarize(self.rows) if v["dataset"] == dataset
                                  and v["mode"] == mode and v["phase"] == "confirm"]
                    winner = min(candidates, key=lambda v: v["median_s"])
                    config, policy = (int(winner["nr"]), int(winner["nw"])), winner["policy"]
                    self.log("SELECTED {} {} {} nr={} nw={}".format(dataset, mode, policy, *config))
                    for repeat in range(1, 6):
                        self.sample(dataset, mode, config, policy, "formal", repeat)
                    if mode == "copy":
                        row = next(v for v in self.rows if v["dataset"] == dataset and v["mode"] == mode
                                   and v["phase"] == "formal" and int(v["repeat"]) == 1)
                        self.verify(row)
        self.run("final-freeze", ["sha256sum", "-c", str(self.source / "frozen.sha256")])
        self.run("final-recheck-freeze", ["sha256sum", "-c", str(self.root / "scripts.sha256")])
        self.log(self.args.stage.upper() + "_COMPLETE")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--source-root", required=True)
    parser.add_argument("--stage", choices=("diagnostic", "recheck"), required=True)
    parser.add_argument("--policies", nargs="+", default=["baseline", "active"],
                        choices=("baseline", "active", "interleave", "node0", "node1"))
    parser.add_argument("--datasets", nargs="+", default=["D1", "D2", "D3"], choices=sorted(base.DATASETS))
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--data", default="/home/bigssd/wzs_data")
    parser.add_argument("--samtools", default="/home/user_home/wzs/software/samtools-1.20/bin/samtools")
    parser.add_argument("--libdeflate", default="/home/user_home/wzs/software/libdeflate-1.20/lib64")
    args = parser.parse_args()
    if args.rounds < 3:
        parser.error("at least three rounds required")
    runner = Recheck(args)
    try:
        runner.execute()
    except BaseException as error:
        runner.log("FAILED " + repr(error))
        raise


if __name__ == "__main__":
    main()
