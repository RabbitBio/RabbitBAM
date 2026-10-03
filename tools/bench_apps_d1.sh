#!/usr/bin/env bash
# Six active core groups; physical preload/output stays outside app core timing.
set -euo pipefail

RUN_DIR=${1:?Usage: bench_apps_d1.sh absolute-run-directory}
RUN_ID=${RUN_DIR##*/}
[[ "$RUN_DIR" = /* && "$RUN_ID" =~ ^[A-Za-z0-9_-]+$ ]]
REMOTE_DIR="online/guoshi/wzs/wzs_data/swbam_bench/D1_apps_current/$RUN_ID"
INPUT=../../../D1-WES-9G/initial/D1-WES-9G.coordinate.bam
EXPECTED_RECORDS=177798571
RESUME=${SWBAM_APPS_RESUME:-0}
[[ "$RESUME" = 0 || "$RESUME" = 1 ]]
SSH=(ssh -o BatchMode=yes -o ConnectTimeout=20
     -o ServerAliveInterval=30 -o ServerAliveCountMax=12 swls)
OUTPUT_OPTIONS=(--io-backend memory --rank-body-backend memory
                --io-output-backend mpiio --compress-level 1)
CURRENT_JOB=
CURRENT_LOG=

mkdir -p "$RUN_DIR/logs"
exec 9>"$RUN_DIR/driver.lock"
flock -n 9 || exit 1
[[ -x "$RUN_DIR/bin/RabbitBAM-MPI" ]]
if [[ "$RESUME" = 1 ]]; then
    [[ -s "$RUN_DIR/results.tsv" ]]
else
    [[ ! -e "$RUN_DIR/results.tsv" ]]
fi
finish_driver() {
    local rc=$?
    if ((rc != 0)); then
        if [[ -z "$CURRENT_JOB" && -f "$CURRENT_LOG" ]]; then
            CURRENT_JOB=$(sed -n 's/.*Job <\([0-9][0-9]*\)> has been submitted.*/\1/p' "$CURRENT_LOG")
        fi
        printf 'FAILED rc=%s current_job=%s\n' "$rc" "$CURRENT_JOB" > "$RUN_DIR/state"
        if [[ "$CURRENT_JOB" =~ ^[0-9]+$ ]]; then
            "${SSH[@]}" "bkill $CURRENT_JOB" || true
        fi
    fi
}
trap finish_driver EXIT
printf 'STARTING\n' > "$RUN_DIR/state"
"${SSH[@]}" "mkdir -p '$REMOTE_DIR/logs' '$REMOTE_DIR/out' '$REMOTE_DIR/tmp'"
rsync -a -e 'ssh -o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=30' \
    "$RUN_DIR/bin/RabbitBAM-MPI" "swls:$REMOTE_DIR/RabbitBAM-MPI"
EXE_REMOTE=$("${SSH[@]}" "readlink -f '$REMOTE_DIR/RabbitBAM-MPI'")
if [[ "$RESUME" != 1 ]]; then
    printf 'kind\tapp\trepeat\tjob_id\tranks\trecords\tkept\tcore_s\tstage43_s\tload_s\tdump_s\toutput_bytes\tsubmit_wall_s\n' > "$RUN_DIR/results.tsv"
    printf 'app\trepeat\treason\n' > "$RUN_DIR/skipped.tsv"
fi

run_case() {
    local kind=$1 app=$2 repeat=$3 expected=$4 output=$5
    shift 5
    if awk -F '\t' -v kind="$kind" -v app="$app" -v repeat="$repeat" '
        NR > 1 && $1 == kind && $2 == app && $3 == repeat { found = 1 }
        END { exit !found }
    ' "$RUN_DIR/results.tsv"; then
        printf 'SKIP %s_%s_r%s (already complete)\n' "$kind" "$app" "$repeat"
        return 0
    fi
    local label="${kind}_${app}_r${repeat}"
    local log="$RUN_DIR/logs/$label.log"
    local job_name="d1app_${RUN_ID##*_}_${label}"
    local start end rc=0 command values output_bytes=0 job
    local -a submit=(bsub -q q_share -I -b -K -J "$job_name"
        -N 6 -np 1 -node 33449-33454 -mpecg 6 -cgsp 64
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
        -timelimit 01:30:00 -o "logs/$label.log"
        -x SWBAM_DIAGNOSTICS=0 "$EXE_REMOTE" "$@")
    printf -v command '%q ' "${submit[@]}"
    printf 'cd %q && %s\n' "$REMOTE_DIR" "$command" >> "$RUN_DIR/commands.sh"
    printf 'RUNNING %s\n' "$label" > "$RUN_DIR/state"
    CURRENT_LOG=$log
    printf 'START %s %s\n' "$label" "$(date -Is)"
    start=$(date +%s)
    timeout --signal=TERM --kill-after=30s 5700 \
        "${SSH[@]}" "cd '$REMOTE_DIR' && $command" > "$log" 2>&1 || rc=$?
    end=$(date +%s)
    job=$(sed -n 's/.*Job <\([0-9][0-9]*\)> has been submitted.*/\1/p' "$log")
    CURRENT_JOB=$job
    if ((rc != 0)) || [[ ! "$job" =~ ^[0-9]+$ ]]; then
        printf 'Submission failed: %s rc=%s job=%s\n' "$label" "$rc" "$job" >&2
        return 1
    fi
    values=$(awk -v expected="$expected" -v app="$app" '
        { sub(/\r$/, "") }
        /ERROR:|exit_code=[1-9]|exit with code: [1-9]/ { failed = 1 }
        /^Complete the 4[.]3 / { stage43 = $NF + 0 }
        /^Complete the total [(]4[.]1~4[.]6[)] cost / { core = $NF + 0 }
        /^444Complete the total body cost / { body = $NF + 0 }
        /^222Complete / { preload = $NF + 0 }
        /^555.*(MPI-IO|Dump memory)/ { dump += $NF + 0 }
        /^SN\tis sorted:/ { split($0, sn, "\t"); sorted = sn[3] + 0; have_sorted = 1 }
        /^Optimized[^ ]* finished[.]/ {
            ++summaries
            for (i = 1; i <= NF; ++i) {
                split($i, kv, "=")
                v[kv[1]] = kv[2]
            }
        }
        /^TOT TIME / { completed = 1 }
        END {
            records = v["total_records"]
            if (records == "") records = v["records"]
            if (records == "") records = v["local_records"]
            kept = v["kept_records"]
            if (kept == "") kept = records
            if (app == "flagstat" || app == "stats" || app ~ /^verify_/) core = body
            if (failed || summaries != 1 || !completed || v["ranks"] != 6 ||
                records != expected || core <= 0 || stage43 <= 0) exit 1
            if ((app == "sort" || app == "collate") && v["mode"] != "memory") exit 1
            if (app == "bam2bam" && kept != records) exit 1
            if (app == "filter" && (kept <= 0 || kept > records)) exit 1
            if ((app == "stats" || app ~ /^verify_/) && (!have_sorted || sorted != 1)) exit 1
            printf "%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f", v["ranks"],
                records, kept, core, stage43, preload, dump
        }
    ' "$log") || {
        printf 'Invalid summary: %s\n' "$label" >&2
        return 1
    }
    if [[ "$output" != - ]]; then
        output_bytes=$("${SSH[@]}" "stat -c %s '$REMOTE_DIR/$output'")
        [[ "$output_bytes" =~ ^[0-9]+$ && "$output_bytes" -gt 0 ]]
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$kind" "$app" "$repeat" "$job" "$values" "$output_bytes" "$((end-start))" \
        >> "$RUN_DIR/results.tsv"
    CURRENT_JOB=
    CURRENT_LOG=
    printf 'DONE %s job=%s wall=%ss metrics=%s %s\n' \
        "$label" "$job" "$((end-start))" "$values" "$(date -Is)"
}

for repeat in 1 2 3; do
    run_case bench flagstat "$repeat" "$EXPECTED_RECORDS" - \
        flagstat -i "$INPUT" --io-backend memory
    run_case bench stats "$repeat" "$EXPECTED_RECORDS" - \
        stats --basic -i "$INPUT" --io-backend memory
    run_case bench bam2bam "$repeat" "$EXPECTED_RECORDS" out/bam2bam.bam \
        run_all -i "$INPUT" -o out/bam2bam.bam "${OUTPUT_OPTIONS[@]}"
    run_case bench filter "$repeat" "$EXPECTED_RECORDS" out/filter.bam \
        run_all -i "$INPUT" -o out/filter.bam --min-mapq 30 "${OUTPUT_OPTIONS[@]}"
    run_case bench collate "$repeat" "$EXPECTED_RECORDS" out/name.bam \
        collate -i "$INPUT" -o out/name.bam -n 66 -m 32G \
        -T tmp/collate "${OUTPUT_OPTIONS[@]}"
    run_case bench fixmate "$repeat" "$EXPECTED_RECORDS" out/fixmate.bam \
        fixmate -m -i out/name.bam -o out/fixmate.bam "${OUTPUT_OPTIONS[@]}"
    run_case bench sort "$repeat" "$EXPECTED_RECORDS" out/fixmate.sorted.bam \
        sort -i out/fixmate.bam -o out/fixmate.sorted.bam -m 32G \
        -T tmp/sort "${OUTPUT_OPTIONS[@]}"
    run_case bench markdup "$repeat" "$EXPECTED_RECORDS" out/markdup.bam \
        markdup --stream -l 300 -m 8G -i out/fixmate.sorted.bam \
        -o out/markdup.bam "${OUTPUT_OPTIONS[@]}"
    run_case bench bam2sam "$repeat" "$EXPECTED_RECORDS" out/initial.sam \
        run_all -i "$INPUT" -o out/initial.sam "${OUTPUT_OPTIONS[@]}"
    sam_bytes=$("${SSH[@]}" "stat -c %s '$REMOTE_DIR/out/initial.sam'")
    [[ "$sam_bytes" =~ ^[0-9]+$ ]]
    if ((sam_bytes <= 72 * 1024 * 1024 * 1024)); then
        run_case bench sam2bam "$repeat" "$EXPECTED_RECORDS" out/sam2bam.bam \
            run_all -i out/initial.sam -o out/sam2bam.bam "${OUTPUT_OPTIONS[@]}"
    else
        printf 'sam2bam\t%s\tSAM size %s exceeds 72 GiB full-memory safety budget\n' \
            "$repeat" "$sam_bytes" >> "$RUN_DIR/skipped.tsv"
    fi
    if ((repeat == 1)); then
        for app in bam2bam markdup; do
            run_case validation "verify_$app" 1 "$EXPECTED_RECORDS" - \
                stats --basic -i "out/$app.bam" --io-backend memory
        done
        kept=$(awk -F '\t' '$1 == "bench" && $2 == "filter" && $3 == 1 {print $7}' "$RUN_DIR/results.tsv")
        run_case validation verify_filter 1 "$kept" - \
            stats --basic -i out/filter.bam --io-backend memory
        if ((sam_bytes <= 72 * 1024 * 1024 * 1024)); then
            run_case validation verify_sam2bam 1 "$EXPECTED_RECORDS" - \
                stats --basic -i out/sam2bam.bam --io-backend memory
        fi
        printf 'PILOT_COMPLETE %s\n' "$(date -Is)" > "$RUN_DIR/pilot_complete"
    fi
done
if [[ $(wc -l < "$RUN_DIR/skipped.tsv") -gt 1 ]]; then
    printf 'COMPLETE_WITH_SKIPS\n' > "$RUN_DIR/state"
else
    printf 'COMPLETE\n' > "$RUN_DIR/state"
fi
printf 'COMPLETE %s\n' "$(date -Is)"
