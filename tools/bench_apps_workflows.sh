#!/usr/bin/env bash
# Run real, independent D2/D3 four-step rounds with the established memory model.
set -euo pipefail

ROOT=${1:?Usage: bench_apps_workflows.sh absolute-run-directory}
[[ "$ROOT" = /* && "${ROOT##*/}" =~ ^[A-Za-z0-9_-]+$ ]]
[[ -x "$ROOT/bin/RabbitBAM-MPI" && -f "$ROOT/driver.sh" && -f "$ROOT/summary.py" ]]
exec 8>"$ROOT/orchestrator.lock"
flock -n 8 || exit 1
trap 'rc=$?; if ((rc != 0)); then printf "FAILED rc=%s\n" "$rc" > "$ROOT/state"; fi' EXIT

export SWBAM_APPS_RANKS=6 SWBAM_APPS_BINS=66 SWBAM_APPS_HEALTH_LEVEL=1
export SWBAM_APPS_NODE_RANGE=${SWBAM_APPS_NODE_RANGE:-auto}
export SWBAM_APPS_BATCH=0 SWBAM_APPS_RESUME=0
export SWBAM_APPS_IO_BACKEND=rank-memory SWBAM_APPS_PIPELINE_ONLY=1
export SWBAM_APPS_SORT_INPUT=out/fixmate.bam SWBAM_APPS_EXPECTED_SORTED=1
export SWBAM_APPS_NO_CONVERSIONS=1 SWBAM_APPS_CONVERSIONS_ONLY=0
export SWBAM_APPS_CONVERSION_CHUNK_BYTES=
export SWBAM_APPS_PILOT_ONLY=0 SWBAM_APPS_SIMPLE_ONLY=0
export SWBAM_APPS_TAIL_ONLY=0 SWBAM_APPS_COLLATE_FIXMATE_ONLY=0
export SWBAM_APPS_COLLATE_DIAG=0 SWBAM_APPS_SORT_MEMORY_DIAG=0
export SWBAM_APPS_MARKDUP_DIAG=0 SWBAM_APPS_VALIDATION_DIAG=0

INPUT_BASE=/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs
SMALL="$INPUT_BASE/code/RabbitBAM/performance_test/WES_0.25G.bam"

run_phase() {
    local name=$1 dataset=$2 input=$3 records=$4 repeats=$5 mode=$6 memory=$7 marked=$8
    local markdup_memory=${9:-$memory} resume=0
    local dir="$ROOT/${ROOT##*/}_$name"
    mkdir -p "$dir/bin"
    cp "$ROOT/bin/RabbitBAM-MPI" "$dir/bin/RabbitBAM-MPI"
    if [[ -s "$dir/results.tsv" ]]; then resume=1; fi
    printf 'RUNNING %s\n' "$name" > "$ROOT/state"
    env SWBAM_APPS_DATASET="$dataset" SWBAM_APPS_INPUT="$input" \
        SWBAM_APPS_RESUME="$resume" \
        SWBAM_APPS_EXPECTED_RECORDS="$records" SWBAM_APPS_REPEATS="$repeats" \
        SWBAM_APPS_REORDER_MODE="$mode" SWBAM_APPS_COLLATE_MEMORY="$memory" \
        SWBAM_APPS_SORT_MEMORY="$memory" SWBAM_APPS_MARKDUP_MEMORY="$markdup_memory" \
        SWBAM_APPS_TIME_LIMIT=06:00:00 \
        /bin/bash "$ROOT/driver.sh" "$dir"
    PYTHONDONTWRITEBYTECODE=1 python3 "$ROOT/summary.py" "$dir" \
        --dataset "$dataset" --repeats "$repeats" --records "$records" --marked "$marked"
}

run_phase smoke_memory D2 "$SMALL" 2668351 1 memory 1G 15695
run_phase smoke_external D3 "$SMALL" 2668351 1 external 128M 15695 8G
run_phase D2 D2 "$INPUT_BASE/wzs_data/D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam" \
    365079492 3 memory 64G 8894488
run_phase D3 D3 "$INPUT_BASE/wzs_data/D3-WGS-85G/initial/D3-WGS-85G.coordinate.bam" \
    1210633370 3 external 8G 42383102
printf 'COMPLETE\n' > "$ROOT/state"
printf 'D2/D3 ordinary workflows complete %s\n' "$(date -Is)"
