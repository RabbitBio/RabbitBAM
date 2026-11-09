#!/usr/bin/env python3
"""Compare decoded-payload cache limits without changing the production binary."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import time

import probe_nonhuman_d3 as probe
from markdup_halo_metrics import summarize
from markdup_optimized_retest import NEW_APP, NEW_SHA, SPECS


NODES = {("D2", "default"): "43334-43339", ("D2", "large"): "43340-43345",
         ("D3", "default"): "43346-43351", ("D3", "large"): "43316-43321"}


def details(log):
    def rows(pattern):
        values = {}
        for rank, line in re.findall(pattern, log, re.M):
            rank = int(rank)
            if rank in values:
                raise RuntimeError("duplicate diagnostic rank")
            values[rank] = dict(re.findall(r"(\w+)=([^\s]+)", line))
        if set(values) != set(range(6)):
            raise RuntimeError("incomplete diagnostic ranks")
        return values

    return dict(pending=rows(r"^\[md-pending-summary rank=(\d+)\] (.*)$"),
                scheduler=rows(r"^\[md-scheduler rank=(\d+)\] (.*)$"),
                stats=rows(r"^\[rank (\d+)\] markdup_stats (.*)$"),
                timing=rows(r"^\[rank (\d+)\] markdup_timing (.*)$"))


def submit(args, root, remote_root, state):
    state("INITIALIZING")
    probe.remote(args.swls, ["bash", "-s"], script="\n".join([
        "set -eu", "test ! -e " + shlex.quote(remote_root),
        "mkdir -p " + " ".join(shlex.quote(remote_root + "/" + name)
                                  for name in ("bin", "logs", "out", "sessions")),
        "test $(sha256sum " + shlex.quote(NEW_APP) + " | awk '{print $1}') = " + NEW_SHA,
        "cp -p " + shlex.quote(NEW_APP) + " " + shlex.quote(remote_root + "/bin/default"),
    ]))
    large_sha = hashlib.sha256(args.large_app.read_bytes()).hexdigest()
    probe.copy(str(args.large_app), args.swls + ":" + remote_root + "/bin/large")
    digests = probe.remote(args.swls, ["sha256sum", remote_root + "/bin/default", remote_root + "/bin/large"])
    if dict((Path(p).name, h) for h, p in map(str.split, digests.splitlines())) != {
            "default": NEW_SHA, "large": large_sha}:
        raise RuntimeError("binary fingerprint mismatch")
    (root / "binaries.sha256").write_text(digests)
    inputs = probe.remote(args.swls, ["stat", "-c", "%n\t%s\t%Y"] +
                          [SPECS[case]["sorted"] for case in ("D2", "D3")])
    (root / "inputs.tsv").write_text(inputs)
    jobs = []
    for (case, variant), nodes in NODES.items():
        label = case + "_" + variant
        logfile = remote_root + "/logs/" + label + ".log"
        command = ["bsub", "-node", nodes, "-N", "6", "-np", "1", "-health", "1",
                   "-mpecg", "6", "-cgsp", "64", "-share_size", "100", "-cross_size", "88000",
                   "-cache_size", "32", "-xmalloc", "-a", "mpi", "-b", "-q", "q_share",
                   "-timelimit", "01:00:00", "-J", "mdbuf_" + label, "-o", logfile,
                   "-x", "SWBAM_MARKDUP_PENDING_DIAG=1",
                   remote_root + "/bin/" + variant, "markdup", "-i", SPECS[case]["sorted"],
                   "-o", remote_root + "/out/" + label + ".bam", "--stream", "-l", "300",
                   "-m", "8G", "--compress-level", "1", "--storage", "ram",
                   "--ram-capacity", "24G", "--ram-ranks", "6", "--ram-session",
                   remote_root + "/sessions/" + label, "--ram-ready-timeout", "5400", "--no-dump"]
        command = ["env", "SWBAM_NODE_MEMORY=88000M"] + command
        submitted = probe.remote(args.swls, ["bash", "-s"], script=
            "set -eu\ncd " + shlex.quote(remote_root) + "\nexec " +
            " ".join(map(shlex.quote, command)) + "\n")
        match = re.search(r"Job <(\d+)>", submitted)
        if not match:
            raise RuntimeError("submission returned no job ID: " + submitted)
        jobs.append(dict(case=case, variant=variant, label=label, job=match[1], argv=command, log=logfile))
        (root / "jobs.json").write_text(json.dumps(jobs, indent=2) + "\n")
        state("SUBMITTED", **{key: jobs[-1][key] for key in ("label", "job")})
    return jobs


def collect(args, root, jobs, state):
    (root / "logs").mkdir(exist_ok=True)
    results, finished = [], set()
    deadline = time.monotonic() + 2 * 3600
    try:
        while time.monotonic() < deadline:
            for job in jobs:
                if job["job"] in finished:
                    continue
                status = probe.remote(args.swls, ["bash", "-s"], script=
                                      "bjobs " + job["job"] + " 2>/dev/null || true\n")
                match = re.search(r"^" + job["job"] + r"\s+(\S+)", status, re.M)
                phase = match[1] if match else ""
                if phase in ("PEND", "STARTING", "RUN", "EXITING"):
                    continue
                fetch = lambda: probe.remote(args.swls, ["cat", job["log"]])
                if phase not in ("", "DONE"):
                    (root / "logs" / (job["label"] + ".failed.log")).write_text(fetch())
                    raise RuntimeError("job " + job["job"] + " finished as " + phase)
                expected = dict(records=SPECS[job["case"]]["records"], marked=SPECS[job["case"]]["marked"])
                log, values = probe.completed_log(fetch, 6, expected, "markdup")
                (root / "logs" / (job["label"] + ".log")).write_text(log)
                values.update(job=job["job"], case=job["case"], variant=job["variant"],
                              halo=summarize(log, 6, expected["records"], expected["marked"]),
                              diagnostics=details(log))
                results.append(values)
                finished.add(job["job"])
                (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
                state("SAMPLE_PASSED", label=job["label"], seconds=values["processing_max_s"])
            if len(finished) == len(jobs):
                for case in ("D2", "D3"):
                    pair = [row for row in results if row["case"] == case]
                    if len({row["output_bytes"] for row in pair}) != 1:
                        raise RuntimeError("complete output sizes differ")
                state("COMPLETE", samples=len(results), full_record_comparison=False)
                return
            time.sleep(60)
        raise RuntimeError("diagnostic experiment exceeded two hours")
    except Exception as error:
        state("FAILED", error=str(error), jobs=[job["job"] for job in jobs],
              note="Inspect surviving jobs; no automatic resubmission or cancellation.")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--large-app", type=Path)
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--collect-only", action="store_true", help="Collect existing job IDs; never submit.")
    args = parser.parse_args()
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    remote_root = probe.SW_DATA + "/swbam_bench/" + root.name

    def state(phase, **extra):
        value = dict(phase=phase, local_time=time.strftime("%F %T %z"), **extra)
        (root / "state.json").write_text(json.dumps(value, indent=2) + "\n")
        print(json.dumps(value), flush=True)

    if args.collect_only:
        jobs = json.loads((root / "jobs.json").read_text())
        actual = probe.remote(args.swls, ["sha256sum", remote_root + "/bin/default", remote_root + "/bin/large"])
        if actual != (root / "binaries.sha256").read_text():
            raise RuntimeError("frozen binaries changed")
        state("COLLECTING", jobs=[job["job"] for job in jobs])
    else:
        if not args.large_app or (root / "jobs.json").exists():
            raise RuntimeError("require --large-app and a new experiment directory")
        jobs = submit(args, root, remote_root, state)
    collect(args, root, jobs, state)


if __name__ == "__main__":
    main()
