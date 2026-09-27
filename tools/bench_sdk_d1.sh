#!/usr/bin/env bash
# Run fixed SDK snapshots serially; keep submission wall separate from core.
set -euo pipefail

RUN_DIR=${1:?Usage: bench_sdk_d1.sh absolute-run-directory}
RUN_ID=${RUN_DIR##*/}
REMOTE_DIR="online/guoshi/wzs/wzs_data/swbam_bench/D1_sdk_current/$RUN_ID"
# Compute nodes resolve the exported online tree relative to the job directory.
INPUT=../../../D1-WES-9G/initial/D1-WES-9G.coordinate.bam

[[ "$RUN_DIR" = /* && "$RUN_ID" =~ ^[A-Za-z0-9_-]+$ ]]
[[ -x "$RUN_DIR/bin/swbam-sdk-optimized-count" ]]
[[ -x "$RUN_DIR/bin/swbam-sdk-optimized-bam2bam" ]]
mkdir -p "$RUN_DIR/logs"
exec 9>"$RUN_DIR/driver.lock"
flock -n 9 || { echo 'Another driver owns this run.' >&2; exit 1; }
[[ ! -e "$RUN_DIR/results.tsv" ]] || {
    echo 'Existing results: use a new run directory, do not overwrite.' >&2
    exit 1
}
printf '%s\n' "$$" > "$RUN_DIR/driver.pid"
printf 'STARTING\n' > "$RUN_DIR/state"
trap 'rc=$?; if ((rc != 0)); then printf "FAILED rc=%d\n" "$rc" > "$RUN_DIR/state"; fi' EXIT
ssh -o BatchMode=yes -o ConnectTimeout=15 swls "mkdir -p '$REMOTE_DIR/logs'"
printf 'mode\tnodes\trepeat\tjob_id\tranks\trecords\tblocks\tinput_bytes\tbody_bytes\tcore_s\tsubmit_wall_s\n' > "$RUN_DIR/results.tsv"

run_case() {
    local mode=$1 n=$2 repeat=$3
    local exe=swbam-sdk-optimized-count prefix=sdk_count
    if [[ "$mode" = rw ]]; then
        exe=swbam-sdk-optimized-bam2bam
        prefix=sdk_bam2bam
    fi
    local label="${mode}_n${n}_r${repeat}"
    local log="$RUN_DIR/logs/${label}.submit.log"
    local start end job values
    printf 'RUNNING %s\n' "$label" > "$RUN_DIR/state"
    printf 'START %s %s\n' "$label" "$(date -Is)"
    start=$(date +%s)
    bsub.py swls -q q_share -I -b -K -J "sdk_d1_${label}" \
        -N "$n" -np 1 -mpecg 6 -cgsp 64 \
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
        -timelimit 00:30:00 -o "logs/${label}.log" \
        -x SWBAM_DIAGNOSTICS=0 \
        "$RUN_DIR/bin/${exe}@${REMOTE_DIR}" "$INPUT" > "$log" 2>&1
    end=$(date +%s)
    job=$(sed -n 's/.*Job <\([0-9][0-9]*\)> has been submitted.*/\1/p' "$log")
    [[ "$job" =~ ^[0-9]+$ ]]
    values=$(awk -v prefix="$prefix" -v n="$n" '
        $1 == prefix {
            ++seen
            for (i = 2; i <= NF; ++i) {
                split($i, field, "="); value[field[1]] = field[2]
            }
        }
        END {
            if (seen != 1 || value["ranks"] != n ||
                value["records"] != 177798571 || value["blocks"] != 814215 ||
                value["input_bytes"] != 9891133304 || value["core"] + 0 <= 0 ||
                (prefix == "sdk_bam2bam" && value["body_bytes"] + 0 <= 0)) exit 1
            printf "%s\t%s\t%s\t%s\t%s\t%s", value["ranks"],
                value["records"], value["blocks"], value["input_bytes"],
                prefix == "sdk_bam2bam" ? value["body_bytes"] : "0", value["core"]
        }' "$log")
    printf '%s\t%d\t%d\t%s\t%s\t%d\n' \
        "$mode" "$n" "$repeat" "$job" "$values" "$((end-start))" >> "$RUN_DIR/results.tsv"
    printf 'DONE %s job=%s wall=%ss metrics=%s %s\n' \
        "$label" "$job" "$((end-start))" "$values" "$(date -Is)"
}

# The gate runs count as formal repeats, so the full matrix is still 36 jobs.
run_case rw 6 1
run_case rw 6 2
if ! awk -F '\t' '
    $1 == "rw" && $2 == 6 {
        ++seen
        if (seen == 1) { first = $10; bytes = $9 }
        else {
            low = first < $10 ? first : $10
            high = first > $10 ? first : $10
            if ($9 != bytes || low <= 0 || high / low > 1.20) exit 1
        }
    }
    END { if (seen != 2) exit 1 }
' "$RUN_DIR/results.tsv"; then
    echo 'Six-group read-write stability gate failed; scaling was not started.' >&2
    exit 1
fi
printf 'SCALING_READY %s\n' "$(date -Is)" > "$RUN_DIR/scaling_ready"
printf 'STABILITY_OK %s\n' "$(date -Is)"

for n in 6 1 2 4 8 12; do
    run_case read "$n" 1
    if ((n != 6)); then run_case rw "$n" 1; fi
done
for n in 12 8 6 4 2 1; do
    if ((n != 6)); then run_case rw "$n" 2; fi
    run_case read "$n" 2
done
for n in 1 2 4 6 8 12; do
    run_case read "$n" 3
    run_case rw "$n" 3
done
printf 'COMPLETE\n' > "$RUN_DIR/state"
printf 'COMPLETE %s\n' "$(date -Is)"
