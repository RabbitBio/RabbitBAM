#!/usr/bin/env bash
# Gate the full D3 matrix on real small-data roundtrip and paired D1 timing.
set -euo pipefail

ROOT=${1:?Usage: bench_apps_d3_conversions.sh absolute-run-directory}
[[ "$ROOT" = /* && "${ROOT##*/}" =~ ^[A-Za-z0-9_-]+$ ]]
[[ -x "$ROOT/bin/RabbitBAM-MPI" && -f "$ROOT/driver.sh" ]]
exec 8>"$ROOT/orchestrator.lock"
flock -n 8 || exit 1
trap 'rc=$?; if ((rc != 0)); then printf "FAILED rc=%s\n" "$rc" > "$ROOT/state"; fi' EXIT

export SWBAM_APPS_DATASET=D3
export SWBAM_APPS_RANKS=6
export SWBAM_APPS_NODE_RANGE=${SWBAM_APPS_NODE_RANGE:-auto}
export SWBAM_APPS_HEALTH_LEVEL=1
export SWBAM_APPS_IO_BACKEND=rank-memory
export SWBAM_APPS_CONVERSIONS_ONLY=1
export SWBAM_APPS_CONVERSION_RETAIN_ALL=0
export SWBAM_APPS_CONVERSION_CHUNK_BYTES=
export SWBAM_APPS_REPEATS=1
export SWBAM_APPS_EXPECTED_SORTED=0
export SWBAM_APPS_TIME_LIMIT=02:00:00
export SWBAM_APPS_RESUME=0
export SWBAM_APPS_NO_CONVERSIONS=0
export SWBAM_APPS_PILOT_ONLY=0
export SWBAM_APPS_SIMPLE_ONLY=0
export SWBAM_APPS_PIPELINE_ONLY=0
export SWBAM_APPS_TAIL_ONLY=0
export SWBAM_APPS_COLLATE_FIXMATE_ONLY=0
export SWBAM_APPS_BATCH=0
export SWBAM_APPS_COLLATE_DIAG=0
export SWBAM_APPS_VALIDATION_DIAG=0

SSH=(ssh -o BatchMode=yes -o ConnectTimeout=20
     -o ServerAliveInterval=30 -o ServerAliveCountMax=12 swls)
REMOTE_BASE=online/guoshi/wzs/wzs_data/swbam_bench/D3_apps_current
INPUT_BASE=/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs
SMALL="$INPUT_BASE/code/RabbitBAM/performance_test/WES_0.25G.bam"
D1="$INPUT_BASE/wzs_data/D1-WES-9G/initial/D1-WES-9G.coordinate.bam"
D3="$INPUT_BASE/wzs_data/D3-WGS-85G/initial/D3-WGS-85G.coordinate.bam"

phase_dir() {
    local dir="$ROOT/${ROOT##*/}_$1"
    mkdir -p "$dir/bin"
    cp "$ROOT/bin/RabbitBAM-MPI" "$dir/bin/RabbitBAM-MPI"
    printf '%s' "$dir"
}

run_phase() {
    local name=$1 input=$2 records=$3 chunk=$4
    local dir resume=0
    dir=$(phase_dir "$name")
    if [[ -s "$dir/results.tsv" ]]; then resume=1; fi
    printf 'RUNNING %s\n' "$name" > "$ROOT/state"
    env SWBAM_APPS_INPUT="$input" SWBAM_APPS_EXPECTED_RECORDS="$records" \
        SWBAM_APPS_RESUME="$resume" \
        SWBAM_APPS_CONVERSION_CHUNK_BYTES="$chunk" \
        /bin/bash "$ROOT/driver.sh" "$dir"
}

run_phase small_baseline "$SMALL" 2668351 ""
run_phase small_chunked "$SMALL" 2668351 8M
run_phase small_roundtrip "$INPUT_BASE/wzs_data/swbam_bench/D3_apps_current/${ROOT##*/}_small_chunked/out/sam2bam.bam" \
    2668351 ""

printf 'CHECKING_SMALL_SAM_CONTENT\n' > "$ROOT/state"
"${SSH[@]}" "cd '$REMOTE_BASE'
sha256sum '${ROOT##*/}_small_baseline/out/initial.sam' '${ROOT##*/}_small_chunked/out/initial.sam' '${ROOT##*/}_small_roundtrip/out/initial.sam'" \
    > "$ROOT/small_sam_sha256.txt"
awk 'NR == 1 { hash = $1 } $1 != hash { failed = 1 } END { exit failed || NR != 3 }' \
    "$ROOT/small_sam_sha256.txt"

# Verify the no-dump path on real data before depending on it for D3 repeats.
SWBAM_APPS_REPEATS=2 run_phase small_discard "$SMALL" 2668351 8M

export SWBAM_APPS_EXPECTED_SORTED=1
run_phase d1_baseline "$D1" 177798571 ""
run_phase d1_chunked "$D1" 177798571 1G
printf 'CHECKING_D1_TIMING\n' > "$ROOT/state"
awk -F '\t' '
    FNR == NR { if ($1 == "bench") baseline[$2] = $8; next }
    $1 == "bench" {
        if (!($2 in baseline) || baseline[$2] <= 0) { failed = 1; next }
        ratio = $8 / baseline[$2]
        printf "%s baseline_s=%.6f chunked_s=%.6f ratio=%.4f\n", $2, baseline[$2], $8, ratio
        if (ratio > 1.20) failed = 1
        ++count
    }
    END { exit failed || count != 2 }
' "$ROOT/${ROOT##*/}_d1_baseline/results.tsv" \
  "$ROOT/${ROOT##*/}_d1_chunked/results.tsv" > "$ROOT/d1_timing_gate.txt"

export SWBAM_APPS_REPEATS=3
export SWBAM_APPS_TIME_LIMIT=06:00:00
run_phase formal "$D3" 1210633370 1G
printf 'COMPLETE\n' > "$ROOT/state"
printf 'D3 conversion matrix complete %s\n' "$(date -Is)"
