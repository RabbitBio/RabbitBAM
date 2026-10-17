#!/usr/bin/env python3
"""Resume frozen marker checks without treating scheduler STARTING as failure."""

import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import re
import time

from markdup_marker_recheck import MarkerCheck
import probe_nonhuman_d3 as probe


class StartupWait:
    def __init__(self, remote, audit, pause=time.sleep, attempts=30):
        self.remote, self.audit = remote, audit
        self.pause, self.attempts = pause, attempts

    def __call__(self, host, argv, timeout=120, script=None):
        query = re.fullmatch(r"bjobs (\d+) 2>/dev/null \|\| true\n", script or "")
        waiting = False
        for attempt in range(self.attempts):
            status = self.remote(host, argv, timeout=timeout, script=script)
            match = re.search(r"^" + query[1] + r"\s+(\S+)", status, re.M) if query else None
            starting = match is not None and match[1] == "STARTING"
            if starting or waiting:
                with self.audit.open("a") as stream:
                    stream.write(json.dumps(dict(time=time.strftime("%Y-%m-%d %H:%M:%S %z"),
                                                 job=query[1], raw_status=status)) + "\n")
            if not starting:
                return status
            waiting = True
            if attempt + 1 < self.attempts:
                self.pause(10)
        raise RuntimeError("scheduler STARTING persisted for five minutes: " + query[1])


@contextmanager
def wait_for_startup(audit):
    # Keep the frozen poller unchanged and return actual, unmodified status text.
    original = probe.remote
    probe.remote = StartupWait(original, audit)
    try:
        yield
    finally:
        probe.remote = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--new-app", type=Path, required=True)
    parser.add_argument("--small-reference", type=Path, required=True)
    parser.add_argument("--swls", default="swls")
    args = parser.parse_args()
    args.resume, args.retry_no_start, args.a4, args.fixmate_app = True, False, "unused", None
    if not (args.root / "jobs.jsonl").is_file():
        raise RuntimeError("recovery requires an existing frozen experiment")
    previous = (args.root / "state.json").read_text()
    recovery = MarkerCheck(args)
    archive = args.root / "state.before-startup-recovery.json"
    if not archive.exists():
        archive.write_text(previous)
    fingerprint = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest = args.root / "startup-recovery.sha256"
    revision = fingerprint + "  " + Path(__file__).name + "\n"
    try:
        if manifest.exists() and manifest.read_text() != revision:
            raise RuntimeError("startup recovery driver changed")
        manifest.write_text(revision)
        recovery.initialize()
        with wait_for_startup(args.root / "scheduler-startup.jsonl"):
            recovery.checks()
    except Exception as error:
        recovery.state("FAILED", job=recovery.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
