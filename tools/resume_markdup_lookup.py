#!/usr/bin/env python3
"""Resume frozen markdup checks with a separate, bounded reread retry."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import threading
import time

from markdup_lookup_recheck import Recheck, remote


class Recovery(Recheck):
    def run(self, label, capacity, argv, expected, check="any", ranks=6):
        if label != "D3_validate":
            return super().run(label, capacity, argv, expected, check, ranks)
        retry = label + "_retry1"
        stop = threading.Event()

        def watch():
            ready_since = None
            while not stop.wait(20):
                job = self.current_job
                if not job:
                    continue
                try:
                    log = remote(self.args.swls, ["cat", self.sw + "/logs/" + retry + ".log"])
                    prepared = set(re.findall(r"\[ram-prepare rank=(\d+)\]", log))
                    if len(prepared) != ranks:
                        continue
                    if "[swbam-processing rank=" in log:
                        return
                    if ready_since is None:
                        ready_since = time.monotonic()
                    if time.monotonic() - ready_since >= 900:
                        remote(self.args.swls, ["bkill", job])
                        (self.root / "reread-watchdog.json").write_text(json.dumps(
                            dict(job=job, reason="no processing completion within 15 minutes after all ranks prepared"),
                            indent=2) + "\n")
                        return
                except (RuntimeError, OSError):
                    continue

        thread = threading.Thread(target=watch, daemon=True)
        thread.start()
        try:
            super().run(retry, capacity, argv, expected, check, ranks)
            shutil.copyfile(self.root / "logs" / (retry + ".log"),
                            self.root / "logs" / (label + ".log"))
        finally:
            stop.set()
            thread.join(timeout=130)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--new-app", type=Path, required=True)
    parser.add_argument("--small-reference", type=Path, required=True)
    parser.add_argument("--swls", default="swls")
    args = parser.parse_args()
    args.resume, args.retry_no_start = True, False
    args.a4, args.fixmate_app = "unused", None
    previous = args.root / "state.json"
    archive = args.root / "state.before-reread-retry.json"
    if previous.exists() and not archive.exists():
        shutil.copyfile(previous, archive)
    fingerprint = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest = args.root / "recovery.sha256"
    revision = fingerprint + "  " + Path(__file__).name + "\n"
    if manifest.exists() and manifest.read_text() != revision:
        raise RuntimeError("recovery driver fingerprint changed")
    manifest.write_text(revision)
    recovery = Recovery(args)
    try:
        recovery.initialize()
        recovery.checks()
    except Exception as error:
        recovery.state("FAILED", job=recovery.current_job, error=str(error))
        raise


if __name__ == "__main__":
    main()
