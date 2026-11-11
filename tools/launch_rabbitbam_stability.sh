#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: launch_rabbitbam_stability.sh ROOT SOURCE_ROOT}
SOURCE=${2:?Frozen original experiment root required}
cd "$ROOT"
if [[ -f driver.pid ]] && kill -0 "$(<driver.pid)" 2>/dev/null; then
    printf 'Recheck driver is already running\n' >&2
    exit 1
fi
nohup python3 -u scripts/recheck_rabbitbam_a4_library.py \
    --root "$ROOT" --source-root "$SOURCE" --stage recheck \
    --policies baseline node0 interleave --rounds 3 \
    > driver.log 2>&1 < /dev/null &
pid=$!
printf '%s\n' "$pid" > driver.pid
sleep 2
kill -0 "$pid"
printf 'driver_pid=%s\nroot=%s\n' "$pid" "$ROOT"
tail -n 8 driver.log
