#!/usr/bin/env python3
"""Wait for fat's validated D3, then relay identical files to a4 and swls."""

import argparse
import datetime
import json
from pathlib import Path
import shlex
import subprocess
import time


SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
       "-o", "ServerAliveInterval=20", "-o", "ServerAliveCountMax=3"]
RSYNC_SSH = shlex.join(SSH)


def command(args, timeout=3600):
    print("COMMAND", shlex.join(args), flush=True)
    subprocess.run(args, check=True, timeout=timeout)


def verify(directory, host=None):
    # The manifest uses relative paths so all three copies use the same checksums.
    code = (
        "import os,subprocess,pathlib; "
        "os.chdir(" + repr(str(directory)) + "); "
        "subprocess.check_call(['sha256sum','-c','SHA256SUMS']); "
        "pathlib.Path('.transfer_verified').touch()"
    )
    args = ["python3", "-c", code]
    if host:
        args = SSH + [host, shlex.join(args)]
    command(args, timeout=7200)


def sync(source, destination):
    args = ["rsync", "-a", "--partial", "--partial-dir=.rsync-partial",
            "--exclude=.transfer_verified",
            "--delay-updates", "--timeout=300", "--info=progress2", "-e", RSYNC_SSH,
            source, destination]
    for attempt in range(1, 6):
        try:
            command(args, timeout=8 * 3600)
            return
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            if attempt == 5:
                raise
            time.sleep(30)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fat", default="wzs@fat.rabbitbio.org")
    parser.add_argument("--fat-root", default="/home/user_home/wzs/wzs_data/D3-mouse-SRR7768905")
    parser.add_argument("--local", type=Path, required=True)
    parser.add_argument("--a4", default="wzs@a4.rabbitbio.org")
    parser.add_argument("--a4-root", default="/home/bigssd/wzs_data/D3-mouse-SRR7768905/initial")
    parser.add_argument("--swls", default="swls")
    parser.add_argument("--swls-root", default="/home/export/base/shisuan/swls-CFD/online/guoshi/wzs/wzs_data/D3-mouse-SRR7768905/initial")
    parser.add_argument("--wait-hours", type=float, default=24)
    args = parser.parse_args()
    args.local.mkdir(parents=True, exist_ok=True)
    status_path = args.local / "relay_status.json"

    def status(stage, **details):
        value = {"stage": stage, "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 **details}
        print(json.dumps(value), flush=True)
        temporary = status_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(value, indent=2) + "\n")
        temporary.replace(status_path)

    try:
        deadline = time.monotonic() + args.wait_hours * 3600
        while time.monotonic() < deadline:
            try:
                result = subprocess.run(
                    SSH + [args.fat, "cat " + shlex.quote(args.fat_root + "/stage")],
                    text=True, capture_output=True, timeout=45)
            except subprocess.TimeoutExpired:
                status("WAIT_FAT", error="SSH status query timed out")
                time.sleep(60)
                continue
            stage = result.stdout.strip()
            if stage.startswith("FAILED"):
                raise RuntimeError("fat preparation failed: " + stage)
            if result.returncode == 0 and stage == "READY":
                break
            status("WAIT_FAT", fat_stage=stage, ssh_returncode=result.returncode)
            time.sleep(60)
        else:
            raise TimeoutError("fat did not finish before the preparation deadline")

        initial = args.local / "initial"
        initial.mkdir(exist_ok=True)
        status("PULL_FAT")
        sync(args.fat + ":" + args.fat_root + "/initial/", str(initial) + "/")
        verify(initial)
        for name, host, directory in [("A4", args.a4, args.a4_root),
                                      ("SWLS", args.swls, args.swls_root)]:
            status("TRANSFER_" + name, host=host, destination=directory)
            command(SSH + [host, "mkdir -p -- " + shlex.quote(directory)], timeout=60)
            sync(str(initial) + "/", host + ":" + directory + "/")
            verify(directory, host)
        status("COMPLETE", a4=args.a4_root, swls=args.swls_root)
    except Exception as error:
        status("FAILED", error=str(error))
        raise


if __name__ == "__main__":
    main()
