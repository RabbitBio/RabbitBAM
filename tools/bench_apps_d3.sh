#!/usr/bin/env bash
# Frozen app binary; external smoke gates only the D3 external phase.
set -euo pipefail

ROOT=${1:?Usage: bench_apps_d3.sh absolute-run-directory}
[[ "$ROOT" = /* && "${ROOT##*/}" =~ ^[A-Za-z0-9_-]+$ ]]
[[ -x "$ROOT/bin/RabbitBAM-MPI" && -f "$ROOT/driver.sh" ]]
exec 8>"$ROOT/driver.lock"
flock -n 8 || exit 1
trap 'rc=$?; if ((rc != 0)); then printf "FAILED rc=%s\n" "$rc" > "$ROOT/state"; fi' EXIT
export SWBAM_APPS_DATASET=D3
export SWBAM_APPS_RANKS=6
export SWBAM_APPS_NODE_RANGE=auto
export SWBAM_APPS_IO_BACKEND=rank-memory
export SWBAM_APPS_COLLATE_DIAG=0
export SWBAM_APPS_SORT_MEMORY_DIAG=0
export SWBAM_APPS_MARKDUP_DIAG=0
export SWBAM_APPS_MARKDUP_MEMORY=8G
export SWBAM_APPS_BINS=66
export SWBAM_APPS_NO_CONVERSIONS=1
export SWBAM_APPS_RESUME=0
export SWBAM_APPS_PILOT_ONLY=0
export SWBAM_APPS_SIMPLE_ONLY=0
export SWBAM_APPS_PIPELINE_ONLY=0
export SWBAM_APPS_TAIL_ONLY=0
export SWBAM_APPS_COLLATE_FIXMATE_ONLY=0

phase_dir() {
    local dir="$ROOT/${ROOT##*/}_$1"
    mkdir -p "$dir/bin"
    cp "$ROOT/bin/RabbitBAM-MPI" "$dir/bin/RabbitBAM-MPI"
    printf '%s' "$dir"
}

printf 'EXTERNAL_SMOKE\n' > "$ROOT/state"
smoke=$(phase_dir external_smoke)
smoke_ok=0
if env SWBAM_APPS_INPUT=/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/code/RabbitBAM/performance_test/WES_0.25G.bam \
    SWBAM_APPS_EXPECTED_RECORDS=2668351 SWBAM_APPS_REPEATS=1 \
    SWBAM_APPS_PIPELINE_ONLY=1 SWBAM_APPS_REORDER_MODE=external \
    SWBAM_APPS_COLLATE_MEMORY=128M SWBAM_APPS_SORT_MEMORY=128M \
    SWBAM_APPS_TIME_LIMIT=00:15:00 /bin/bash "$ROOT/driver.sh" "$smoke"; then
    smoke_ok=1
else
    printf 'External smoke failed; D3 external jobs will not be submitted.\n' >&2
fi

printf 'D3_SIMPLE_MATRIX\n' > "$ROOT/state"
simple=$(phase_dir simple)
env SWBAM_APPS_EXPECTED_RECORDS=1210633370 SWBAM_APPS_REPEATS=3 \
    SWBAM_APPS_SIMPLE_ONLY=1 SWBAM_APPS_REORDER_MODE=memory \
    SWBAM_APPS_TIME_LIMIT=02:00:00 /bin/bash "$ROOT/driver.sh" "$simple"

if ((smoke_ok == 0)); then
    printf 'COMPLETE_SIMPLE_EXTERNAL_BLOCKED\n' > "$ROOT/state"
    exit 0
fi
printf 'D3_EXTERNAL_MATRIX\n' > "$ROOT/state"
external=$(phase_dir external)
env SWBAM_APPS_EXPECTED_RECORDS=1210633370 SWBAM_APPS_REPEATS=3 \
    SWBAM_APPS_PIPELINE_ONLY=1 SWBAM_APPS_REORDER_MODE=external \
    SWBAM_APPS_COLLATE_MEMORY=8G SWBAM_APPS_SORT_MEMORY=8G \
    SWBAM_APPS_TIME_LIMIT=06:00:00 /bin/bash "$ROOT/driver.sh" "$external"
printf 'COMPLETE_SELECTED_D3\n' > "$ROOT/state"
printf 'D3 selected app matrix complete %s\n' "$(date -Is)"
