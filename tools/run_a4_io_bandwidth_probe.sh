#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:?Usage: run_a4_io_bandwidth_probe.sh NEW_ABSOLUTE_DIRECTORY [SIZE_MIB] [REPEATS] [full|buffered-write]}
SIZE_MIB=${2:-32768}
REPEATS=${3:-3}
MODE=${4:-full}
[[ "$ROOT" == /* && "$SIZE_MIB" =~ ^[0-9]+$ && "$REPEATS" =~ ^[0-9]+$ ]]
[[ "$MODE" == full || "$MODE" == buffered-write ]]
((SIZE_MIB >= 64 && SIZE_MIB <= 65536 && REPEATS >= 1 && REPEATS <= 10))
command -v fio >/dev/null
command -v python3 >/dev/null
mkdir "$ROOT"
cd "$ROOT"
mkdir logs
FILE="$ROOT/probe.bin"
MONITOR_PID=
cleanup() {
    if [[ -n "$MONITOR_PID" ]]; then
        kill "$MONITOR_PID" 2>/dev/null || true
        wait "$MONITOR_PID" 2>/dev/null || true
    fi
    rm -f -- "$FILE"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

{
    hostname
    date -Is
    fio --version
    free -h
    df -h "$ROOT"
    findmnt -no SOURCE,FSTYPE,OPTIONS --target "$ROOT"
    lsblk -d -o NAME,MODEL,SIZE,ROTA
    sysctl vm.dirty_ratio vm.dirty_background_ratio vm.dirty_bytes vm.dirty_background_bytes
    uptime
    ps -eo user,pid,pcpu,comm --sort=-pcpu | head -20 || true
} > environment.txt
sha256sum "$0" > script.sha256
if command -v iostat >/dev/null; then
    iostat -dx 2 > iostat.log &
    MONITOR_PID=$!
fi
printf 'name\trepeat\tchunk_bytes\tbytes\tseconds\tGBps\tjob_runtime_ms\n' > results.tsv
common=(--filename="$FILE" --size="$((SIZE_MIB * 1024 * 1024))"
    --ioengine=psync --iodepth=1 --numjobs=1 --thread=1
    --fallocate=posix --zero_buffers=0 --refill_buffers=0 --randrepeat=1
    --output-format=json --group_reporting=1)

run_case() {
    local name=$1 repeat=$2 chunk=$3 rw=$4 direct=$5 invalidate=$6
    local log="$ROOT/logs/${name}_r${repeat}.json"
    local sync=${7:-0}
    if [[ "$rw" == write && $# -lt 7 ]]; then sync=1; fi
    local command=(fio "${common[@]}" --name="$name" --rw="$rw" --bs="$chunk"
        --direct="$direct" --invalidate="$invalidate" --end_fsync="$sync" --output="$log")
    printf '%q ' "${command[@]}" >> commands.sh
    printf '\n' >> commands.sh
    printf '%s START %s repeat=%s\n' "$(date -Is)" "$name" "$repeat"
    timeout --signal=TERM --kill-after=10s 600s "${command[@]}"
    python3 - "$log" "$name" "$repeat" "$chunk" "$rw" "$((SIZE_MIB * 1024 * 1024))" <<'PY' | tee -a results.tsv
import json
import sys
path, name, repeat, chunk, direction, expected = sys.argv[1:]
with open(path) as stream:
    data = json.load(stream)
assert len(data['jobs']) == 1, 'expected one worker'
job = data['jobs'][0]
assert job['error'] == 0, 'fio reported an I/O error'
io = job[direction]
assert io['io_bytes'] == int(expected), 'incomplete file transfer'
assert io['runtime'] > 0, 'invalid I/O runtime'
seconds = io['runtime'] / 1000.0
print('{}\t{}\t{}\t{}\t{:.6f}\t{:.6f}\t{}'.format(
    name, repeat, chunk, io['io_bytes'], seconds,
    io['io_bytes'] / seconds / 1e9, job['job_runtime']))
PY
}

for ((repeat=1; repeat<=REPEATS; ++repeat)); do
    for chunk in 65536 8388608; do
        if [[ "$MODE" == buffered-write ]]; then
            rm -f -- "$FILE"
            run_case "buffered_write_no_sync_${chunk}" "$repeat" "$chunk" write 0 1 0
            # Complete persistence outside the sample before recycling this file.
            python3 - "$FILE" <<'PY' | tee -a flush-after-sample.log
import os
import sys
import time
fd = os.open(sys.argv[1], os.O_RDWR)
try:
    start = time.monotonic()
    os.fsync(fd)
    print('POST_SAMPLE_FSYNC seconds={:.6f} outside_bandwidth_timer=1'.format(
        time.monotonic() - start))
finally:
    os.close(fd)
PY
            continue
        fi
        # Only this file's cache is invalidated; no system-wide drop_caches.
        run_case "buffered_write_sync_${chunk}" "$repeat" "$chunk" write 0 1
        run_case "buffered_read_after_invalidate_${chunk}" "$repeat" "$chunk" read 0 1
        run_case "buffered_read_repeat_${chunk}" "$repeat" "$chunk" read 0 0
        run_case "direct_write_sync_${chunk}" "$repeat" "$chunk" write 1 1
        run_case "direct_read_${chunk}" "$repeat" "$chunk" read 1 1
    done
done

# Full-file CRC verification is separate from the bandwidth samples.
verify=(fio "${common[@]}" --name=content_verify --rw=write --bs=8388608
    --direct=1 --end_fsync=1 --verify=crc32c --do_verify=1 --verify_fatal=1
    --verify_state_save=0 --output="$ROOT/logs/content_verify.json")
printf '%q ' "${verify[@]}" >> commands.sh
printf '\n' >> commands.sh
printf '%s START full-file CRC verification\n' "$(date -Is)"
timeout --signal=TERM --kill-after=10s 600s "${verify[@]}"
stat -c 'file=%n size=%s allocated_512B_blocks=%b' "$FILE" > file-allocation.txt
python3 - "$ROOT/logs/content_verify.json" "$((SIZE_MIB * 1024 * 1024))" <<'PY'
import json
import sys
with open(sys.argv[1]) as stream:
    job = json.load(stream)['jobs'][0]
assert job['error'] == 0, 'CRC verification failed'
assert job['write']['io_bytes'] == int(sys.argv[2]), 'verification write incomplete'
assert job['read']['io_bytes'] == int(sys.argv[2]), 'verification read incomplete'
print('IO_CONTENT_VERIFIED bytes={}'.format(sys.argv[2]))
PY
python3 - results.tsv "$REPEATS" "$MODE" <<'PY' | tee summary.tsv
import csv
import statistics
import sys
from collections import defaultdict
samples = defaultdict(list)
with open(sys.argv[1]) as stream:
    for row in csv.DictReader(stream, delimiter='\t'):
        samples[row['name']].append(float(row['GBps']))
assert len(samples) == (10 if sys.argv[3] == 'full' else 2), 'missing test modes'
print('name\tsamples\tmedian_GBps\tmin_GBps\tmax_GBps')
for name, values in sorted(samples.items()):
    assert len(values) == int(sys.argv[2]), 'missing repeats'
    print('{}\t{}\t{:.6f}\t{:.6f}\t{:.6f}'.format(
        name, len(values), statistics.median(values), min(values), max(values)))
PY
cleanup
MONITOR_PID=
trap - EXIT
printf '%s IO_PROBE_COMPLETE status=0 data_file_removed=1\n' "$(date -Is)"
