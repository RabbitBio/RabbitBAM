#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: launch_rabbitbam_a4_library.sh ROOT}
cd "$ROOT"
if [[ -f driver.pid ]] && kill -0 "$(<driver.pid)" 2>/dev/null; then
    printf 'Driver is already alive: %s\n' "$(<driver.pid)" >&2
    exit 1
fi
nohup python3 -u scripts/bench_rabbitbam_a4_library.py --root "$ROOT" \
    > driver.log 2>&1 < /dev/null &
pid=$!
printf '%s\n' "$pid" > driver.pid
sleep 2
kill -0 "$pid"
printf 'driver_pid=%s\nroot=%s\n' "$pid" "$ROOT"
tail -n 8 driver.log
