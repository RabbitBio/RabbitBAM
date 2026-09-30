#!/usr/bin/env bash
# Run fixed SDK snapshots serially; keep submission wall separate from core.
set -euo pipefail

RUN_DIR=${1:?Usage: bench_sdk_d1.sh absolute-run-directory}
RUN_ID=${RUN_DIR##*/}
DATASET=${SWBAM_BENCH_DATASET:-D1}
case "$DATASET" in
    D1)
        INPUT=../../../D1-WES-9G/initial/D1-WES-9G.coordinate.bam
        EXPECT_RECORDS=177798571
        EXPECT_INPUT_BYTES=9891133304
        EXPECT_BLOCKS=814215
        ;;
    D2)
        INPUT=../../../D2-WGS-25G/initial/D2-WGS-25G.coordinate.bam
        EXPECT_RECORDS=365079492
        EXPECT_INPUT_BYTES=35255807945
        EXPECT_BLOCKS=1927256
        ;;
    D3)
        INPUT=../../../D3-WGS-85G/initial/D3-WGS-85G.coordinate.bam
        EXPECT_RECORDS=1210633370
        EXPECT_INPUT_BYTES=89994058306
        EXPECT_BLOCKS=${SWBAM_EXPECT_BLOCKS:?Set SWBAM_EXPECT_BLOCKS from the D3 pilot}
        ;;
    *) echo "Unsupported dataset: $DATASET" >&2; exit 2 ;;
esac
REMOTE_DIR="online/guoshi/wzs/wzs_data/swbam_bench/${DATASET}_sdk_current/$RUN_ID"
# Compute nodes resolve the exported online tree relative to the job directory.
FIXED_NODE_BASE=${SWBAM_FIXED_NODE_BASE:-}
NODE_HEALTH=${SWBAM_NODE_HEALTH:-}
JOB_TIMEOUT_SECONDS=${SWBAM_JOB_TIMEOUT_SECONDS:-600}
JOB_TIMELIMIT=${SWBAM_JOB_TIMELIMIT:-00:30:00}
RESUME=${SWBAM_BENCH_RESUME:-0}
MAX_ATTEMPTS=${SWBAM_MAX_ATTEMPTS:-3}
D3_ALLOW_SLOW_IO=${SWBAM_D3_ALLOW_SLOW_IO:-0}
SSH_KEEPALIVE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/ssh_keepalive" && pwd)
SDK_OPTIONS=()
JOB_ENV=(-x SWBAM_DIAGNOSTICS=0)
if [[ "$DATASET" = D3 ]]; then
    SDK_OPTIONS=(--chunked-memory)
    JOB_ENV=(-x SWBAM_SDK_INDEX_CACHE=../D3-WGS-85G.bgzfi)
fi

[[ "$RUN_DIR" = /* && "$RUN_ID" =~ ^[A-Za-z0-9_-]+$ ]]
[[ -z "$FIXED_NODE_BASE" || "$FIXED_NODE_BASE" =~ ^[0-9]+$ ]]
[[ -z "$NODE_HEALTH" || "$NODE_HEALTH" =~ ^[012]$ ]]
[[ "$JOB_TIMEOUT_SECONDS" =~ ^[0-9]+$ && "$JOB_TIMEOUT_SECONDS" -gt 0 ]]
[[ "$RESUME" = 0 || "$RESUME" = 1 ]]
[[ "$D3_ALLOW_SLOW_IO" = 0 || "$D3_ALLOW_SLOW_IO" = 1 ]]
[[ "$MAX_ATTEMPTS" =~ ^[0-9]+$ && "$MAX_ATTEMPTS" -gt 0 ]]
[[ -x "$RUN_DIR/bin/swbam-sdk-optimized-count" ]]
[[ -x "$RUN_DIR/bin/swbam-sdk-optimized-bam2bam" ]]
mkdir -p "$RUN_DIR/logs"
exec 9>"$RUN_DIR/driver.lock"
flock -n 9 || { echo 'Another driver owns this run.' >&2; exit 1; }
if [[ -e "$RUN_DIR/results.tsv" && "$RESUME" != 1 ]]; then
    echo 'Existing results: set SWBAM_BENCH_RESUME=1 or use a new run directory.' >&2
    exit 1
fi
if [[ "$RESUME" = 1 ]]; then
    [[ -s "$RUN_DIR/results.tsv" ]] || {
        echo 'Cannot resume: results.tsv is missing or empty.' >&2; exit 1;
    }
    awk -F '\t' -v records="$EXPECT_RECORDS" -v blocks="$EXPECT_BLOCKS" \
        -v bytes="$EXPECT_INPUT_BYTES" '
        NR == 1 { if ($1 != "mode" || $2 != "nodes" || $10 != "core_s") exit 1; next }
        {
            key = $1 SUBSEP $2 SUBSEP $3
            if (seen[key]++ || ($1 != "read" && $1 != "rw") ||
                $2 !~ /^(1|2|4|6|8|12)$/ || $3 !~ /^[123]$/ ||
                $6 != records || $7 != blocks || $8 != bytes ||
                $10 + 0 <= 0) exit 1
        }
    ' "$RUN_DIR/results.tsv" || {
        echo 'Cannot resume: existing results failed validation.' >&2; exit 1;
    }
fi
printf '%s\n' "$$" > "$RUN_DIR/driver.pid"
printf 'STARTING\n' > "$RUN_DIR/state"
trap 'rc=$?; if ((rc != 0)); then printf "FAILED rc=%d\n" "$rc" > "$RUN_DIR/state"; fi' EXIT
ssh -o BatchMode=yes -o ConnectTimeout=15 swls "mkdir -p '$REMOTE_DIR/logs'"
if [[ "$RESUME" != 1 ]]; then
    printf 'mode\tnodes\trepeat\tjob_id\tranks\trecords\tblocks\tinput_bytes\tbody_bytes\tcore_s\tsimulated_total_s\tscan_plan_s\tindex_prepare_s\tstage_read_s\tstage_prepare_s\tbody_setup_s\tpeak_chunk_bytes\tmax_chunks\tsubmit_wall_s\n' > "$RUN_DIR/results.tsv"
fi

run_once() {
    local mode=$1 n=$2 repeat=$3 attempt=$4
    local exe=swbam-sdk-optimized-count prefix=sdk_count
    if [[ "$mode" = rw ]]; then
        exe=swbam-sdk-optimized-bam2bam
        prefix=sdk_bam2bam
    fi
    local label="${mode}_n${n}_r${repeat}"
    local attempt_label="${label}_a${attempt}"
    local job_name="sdk_${DATASET,,}_${attempt_label}"
    local log="$RUN_DIR/logs/${attempt_label}.submit.log"
    local start end job values submit_rc=0 submit_pid job_log
    local case_timeout=$JOB_TIMEOUT_SECONDS cap
    if [[ "$DATASET" = D3 && "$D3_ALLOW_SLOW_IO" != 1 ]]; then
        if ((n == 1)); then
            cap=$([[ "$mode" = rw ]] && printf 2400 || printf 1800)
        elif ((n == 2)); then
            cap=$([[ "$mode" = rw ]] && printf 1800 || printf 1200)
        elif ((n == 4)); then
            cap=$([[ "$mode" = rw ]] && printf 1500 || printf 1200)
        else
            cap=900
        fi
        if ((case_timeout > cap)); then case_timeout=$cap; fi
    fi
    local -a node_args=() health_args=()
    if [[ -n "$FIXED_NODE_BASE" ]]; then
        node_args=(-node "${FIXED_NODE_BASE}-$((FIXED_NODE_BASE + n - 1))")
    fi
    if [[ -n "$NODE_HEALTH" ]]; then health_args=(-health "$NODE_HEALTH"); fi
    printf 'RUNNING %s\n' "$attempt_label" > "$RUN_DIR/state"
    printf 'START %s timeout=%ss %s\n' "$attempt_label" "$case_timeout" "$(date -Is)"
    start=$(date +%s)
    timeout --signal=TERM --kill-after=10s "$case_timeout" \
    env PATH="$SSH_KEEPALIVE_DIR:$PATH" \
    bsub.py swls -q q_share -I -b -K -J "$job_name" \
        -N "$n" -np 1 -mpecg 6 -cgsp 64 \
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
        "${node_args[@]}" "${health_args[@]}" \
        -timelimit "$JOB_TIMELIMIT" -o "logs/${attempt_label}.log" \
        "${JOB_ENV[@]}" \
        "$RUN_DIR/bin/${exe}@${REMOTE_DIR}" "$INPUT" \
        "${SDK_OPTIONS[@]}" > "$log" 2>&1 &
    submit_pid=$!
    job=
    while kill -0 "$submit_pid" 2>/dev/null; do
        sleep 30
        if [[ -z "$job" ]]; then
            job=$(sed -n 's/.*Job <\([0-9][0-9]*\)> has been submitted.*/\1/p' "$log")
            if [[ -z "$job" ]]; then
                job=$(ssh -o BatchMode=yes swls \
                    "bjobs -w -J '$job_name' 2>/dev/null | awk 'NR==3 {print \$1}'" 2>/dev/null || true)
            fi
        fi
        if ssh -o BatchMode=yes swls \
            "grep -q '^${prefix} ranks=' '${REMOTE_DIR}/logs/${attempt_label}.log'" \
            2>/dev/null; then
            # Some LSF clients keep waiting after the result is written.
            sleep 30
            if kill -0 "$submit_pid" 2>/dev/null; then
                kill "$submit_pid" 2>/dev/null || true
            fi
            break
        fi
    done
    wait "$submit_pid" || submit_rc=$?
    end=$(date +%s)
    job_log="$RUN_DIR/logs/${attempt_label}.job.log"
    scp -q "swls:${REMOTE_DIR}/logs/${attempt_label}.log" "$job_log" || true
    if [[ ! -s "$job_log" ]]; then
        job_log=$log
    fi
    if [[ -z "$job" ]]; then
        job=$(sed -n 's/.*Job <\([0-9][0-9]*\)> has been submitted.*/\1/p' "$log")
    fi
    if [[ ! "$job" =~ ^[0-9]+$ ]]; then
        echo "Cannot identify job for $attempt_label (submit rc=$submit_rc)." >&2
        return 1
    fi
    values=$(awk -v prefix="$prefix" -v n="$n" -v dataset="$DATASET" \
        -v records="$EXPECT_RECORDS" -v blocks="$EXPECT_BLOCKS" \
        -v input_bytes="$EXPECT_INPUT_BYTES" '
        $1 == prefix {
            ++seen
            for (i = 2; i <= NF; ++i) {
                split($i, field, "="); value[field[1]] = field[2]
            }
        }
        END {
            if (seen != 1 || value["ranks"] != n ||
                value["records"] != records || value["blocks"] != blocks ||
                value["input_bytes"] != input_bytes || value["core"] + 0 <= 0 ||
                (prefix == "sdk_bam2bam" && value["body_bytes"] + 0 <= 0) ||
                (dataset == "D3" &&
                 (value["peak_chunk_bytes"] + 0 <= 0 ||
                  value["peak_chunk_bytes"] + 0 > 2147483648 ||
                  value["max_chunks"] + 0 < 1 ||
                  value["simulated_total"] + 0 < value["core"] + 0))) exit 1
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", value["ranks"],
                value["records"], value["blocks"], value["input_bytes"],
                prefix == "sdk_bam2bam" ? value["body_bytes"] : "0", value["core"],
                value["simulated_total"], value["scan_plan_wall"],
                value["index_prepare"], value["stage_read"],
                value["stage_prepare"], value["body_setup"],
                value["peak_chunk_bytes"], value["max_chunks"]
        }' "$job_log") || {
            ssh -o BatchMode=yes swls "bkill $job" || true
            echo "Invalid or missing SDK metrics: $attempt_label (submit rc=$submit_rc)." >&2
            return 1
        }
    printf '%s\t%d\t%d\t%s\t%s\t%d\n' \
        "$mode" "$n" "$repeat" "$job" "$values" "$((end-start))" >> "$RUN_DIR/results.tsv"
    printf 'DONE %s job=%s wall=%ss metrics=%s %s\n' \
        "$attempt_label" "$job" "$((end-start))" "$values" "$(date -Is)"
}

run_case() {
    local mode=$1 n=$2 repeat=$3 attempt
    if awk -F '\t' -v mode="$mode" -v n="$n" -v repeat="$repeat" '
        NR > 1 && $1 == mode && $2 == n && $3 == repeat { found = 1 }
        END { exit !found }
    ' "$RUN_DIR/results.tsv"; then
        printf 'SKIP %s_n%s_r%s (already complete)\n' "$mode" "$n" "$repeat"
        return 0
    fi
    for ((attempt = 1; attempt <= MAX_ATTEMPTS; ++attempt)); do
        if run_once "$mode" "$n" "$repeat" "$attempt"; then return 0; fi
        printf 'RETRY %s_n%s_r%s attempt=%d/%d\n' \
            "$mode" "$n" "$repeat" "$attempt" "$MAX_ATTEMPTS" >&2
    done
    return 1
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
