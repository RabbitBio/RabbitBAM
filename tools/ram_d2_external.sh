#!/usr/bin/env bash
set -euo pipefail
export SWBAM_NODE_MEMORY=${SWBAM_NODE_MEMORY:-88000M}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${1:?Usage: ram_d2_external.sh ABS_ROOT INPUT_BAM FIXMATE_BAM}
INPUT=${2:?Input BAM required}
FIXMATE=${3:?Fixmate BAM required}
APP="$ROOT/bin/RabbitBAM-MPI"
MARKDUP_APP="$ROOT/bin/RabbitBAM-MPI-markdup"
RECORDS=365079492
MARKED=8894488
POOL=${BENCH_EXTERNAL_POOL:-64G}
MEMORY=${BENCH_EXTERNAL_MEMORY:-8G}
REPEATS=${BENCH_REPEATS:-3}
CURRENT=initialization
CURRENT_JOB=
[[ "$ROOT" = /* && "$INPUT" = /* && "$FIXMATE" = /* &&
   -s "$INPUT" && -s "$FIXMATE" && -x "$APP" && -x "$MARKDUP_APP" ]]
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions"
mkdir "$ROOT/driver.lock"
trap 'rc=$?; if ((rc)); then printf "FAILED %s job=%s\n" "$CURRENT" "$CURRENT_JOB" > "$ROOT/state"; fi; rmdir "$ROOT/driver.lock"' EXIT
cd "$ROOT"
if [[ ${BENCH_RESUME:-0} == 1 ]]; then
    sha256sum -c binaries.sha256
    grep -Fx "dataset=D2 ranks=6 repeats=$REPEATS bins=66 compression=1 external_memory=$MEMORY external_pool=$POOL timing=mpi-interval/processing/core" configuration.txt >/dev/null
else
    [[ ! -e results.tsv ]]
    if [[ -e binaries.sha256 ]]; then sha256sum -c binaries.sha256; else
        sha256sum bin/RabbitBAM-MPI bin/RabbitBAM-MPI-markdup > binaries.sha256
    fi
    printf 'kind\tapp\trepeat\tranks\tjob_id\tcapacity\te2e_max_s\tcore_s\tprepare_max_s\tpool_peak_max_bytes\toutput_bytes\tprocessing_max_s\ttemp_read_bytes\ttemp_write_bytes\twork_peak_bytes\truns\n' > results.tsv
    printf 'repeat\tcollate_s\tfixmate_s\tsort_s\tmarkdup_s\ttotal_s\tcollate_processing_s\tfixmate_processing_s\tsort_processing_s\tmarkdup_processing_s\tprocessing_total_s\n' > workflows.tsv
    stat -c '%n\t%s bytes\tmtime=%Y' "$INPUT" "$FIXMATE" > inputs.txt
    printf 'dataset=D2 ranks=6 repeats=%s bins=66 compression=1 external_memory=%s external_pool=%s timing=mpi-interval/processing/core\n' "$REPEATS" "$MEMORY" "$POOL" > configuration.txt
    if [[ -n ${BENCH_FIRST_PROBE_JOB:-} ]]; then
        printf 'probe_collate\t%s\n' "$BENCH_FIRST_PROBE_JOB" >> jobs.tsv
    fi
fi

metrics() {
    local check=$1 log=$2
    [[ -s "$log" ]] || return 1
    if grep -Eq '(ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT)' "$log"; then
        tail -60 "$log" >&2
        return 1
    fi
    local filters=(-f "$SCRIPT_DIR/ram_metrics.awk")
    if [[ "$check" == external ]]; then
        filters+=(-f "$SCRIPT_DIR/ram_external_metrics.awk")
    fi
    awk -v ranks=6 -v records="$RECORDS" -v marked="$MARKED" \
        -v check="$check" "${filters[@]}" "$log"
}

run() {
    local kind=$1 app=$2 repeat=$3 capacity=$4 check=$5
    shift 5
    local label="${kind}_${app}" job= status= log data
    if [[ "$kind" != probe && "$kind" != probe_validation ]]; then label+="_r${repeat}"; fi
    CURRENT=$label
    CURRENT_JOB=
    log="$ROOT/logs/$label.log"
    if awk -F '\t' -v k="$kind" -v a="$app" -v r="$repeat" \
        '$1==k && $2==a && $3==r{found=1}END{exit !found}' results.tsv; then
        metrics "$check" "$log" >/dev/null
        printf 'REUSE verified %s\n' "$label"
        return
    fi
    job=$(awk -F '\t' -v label="$label" '$1==label{print $2}' jobs.tsv 2>/dev/null || true)
    printf 'RUNNING %s\n' "$label" > state
    local command=(bsub -N 6 -np 1 -health 1 -mpecg 6 -cgsp 64
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
        -a mpi -b -K -q q_share -timelimit 03:00:00
        -J "d2ram_external_$label" -o "$log"
        "$@" --storage ram --ram-capacity "$capacity" --ram-ranks 6
        --ram-session "$ROOT/sessions/$label" --ram-ready-timeout 5400)
    if [[ -z "$job" ]]; then
        [[ ! -e "$log" && ! -e "$ROOT/sessions/$label" ]]
        printf '%q ' "${command[@]}" >> commands.sh
        printf '\n' >> commands.sh
        printf '%s START %s capacity=%s\n' "$(date -Is)" "$label" "$capacity"
        local submission
        submission=$("${command[@]}")
        printf '%s\n' "$submission"
        job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
        [[ "$job" =~ ^[0-9]+$ ]]
        printf '%s\t%s\n' "$label" "$job" >> jobs.tsv
    elif [[ "$label" == probe_collate && ! -e commands.sh ]]; then
        printf '%q ' "${command[@]}" >> commands.sh
        printf '\n' >> commands.sh
    fi
    [[ "$job" =~ ^[0-9]+$ ]]
    CURRENT_JOB=$job
    local unknown=0
    while true; do
        status=$(bjobs "$job" 2>/dev/null | awk -v j="$job" '$1==j{print $2}' || true)
        case "$status" in
            PEND|RUN|EXITING) unknown=0; sleep 5 ;;
            DONE) break ;;
            "")
                if [[ -s "$log" ]] && grep -q "Job $job has been finished" "$log"; then break; fi
                unknown=$((unknown+1)); ((unknown<40)); sleep 5 ;;
            *) printf 'Scheduler failure %s job=%s\n' "$status" "$job" >&2; return 1 ;;
        esac
    done
    for ((tries=0; tries<40; ++tries)); do
        if [[ -s "$log" ]] && [[ $(awk '/\[ram-storage rank=/{n++}END{print n+0}' "$log") == 6 ]]; then break; fi
        sleep 3
    done
    data=$(metrics "$check" "$log")
    if [[ "$check" != external ]]; then data+=$'\t0\t0\t0\t0'; fi
    printf '%s\t%s\t%s\t6\t%s\t%s\t%s\n' "$kind" "$app" "$repeat" "$job" "$capacity" "$data" >> results.tsv
    printf '%s DONE %s job=%s metrics=%s\n' "$(date -Is)" "$label" "$job" "$data"
    CURRENT_JOB=
}

run probe collate 0 "$POOL" external "$APP" collate -i "$INPUT" \
    -o "$ROOT/out/probe.name.bam" -n 66 -m "$MEMORY" --compress-level 1
[[ -s "$ROOT/out/probe.name.bam" ]]
run probe fixmate 0 24G write "$APP" fixmate -i "$ROOT/out/probe.name.bam" \
    -o "$ROOT/out/probe.fixmate.bam" -m --compress-level 1
[[ -s "$ROOT/out/probe.fixmate.bam" ]]
run probe sort 0 "$POOL" external "$APP" sort -i "$ROOT/out/probe.fixmate.bam" \
    -o "$ROOT/out/probe.sorted.bam" -m "$MEMORY" --compress-level 1
[[ -s "$ROOT/out/probe.sorted.bam" ]]
run probe_validation sort 0 8G sorted "$APP" stats --basic -i "$ROOT/out/probe.sorted.bam"
run probe markdup 0 24G markdup "$MARKDUP_APP" markdup -i "$ROOT/out/probe.sorted.bam" \
    -o "$ROOT/out/probe.markdup.bam" --stream -l 300 -m 8G --compress-level 1
[[ -s "$ROOT/out/probe.markdup.bam" ]]
run probe_validation markdup 0 8G duplicates "$APP" stats --basic -i "$ROOT/out/probe.markdup.bam"
printf '%s CAPACITY AND CORRECTNESS GATES PASSED\n' "$(date -Is)" > capacity_gate.txt

for ((repeat=1; repeat<=REPEATS; ++repeat)); do
    for app in collate sort; do
        input="$INPUT"
        if [[ "$app" == sort ]]; then input="$FIXMATE"; fi
        output="$ROOT/out/${app}_r${repeat}.bam"
        # Bash 4.2 treats an empty array as unset under nounset.
        args=("$APP" "$app" -i "$input" -o "$output")
        if [[ "$app" == collate ]]; then args+=(-n 66); fi
        args+=(-m "$MEMORY" --compress-level 1)
        if ((repeat!=1)); then args+=(--no-dump); fi
        run app "$app" "$repeat" "$POOL" external "${args[@]}"
        if ((repeat==1)); then
            [[ -s "$output" ]]
            check=sorted
            validation_capacity=8G
            if [[ "$app" == collate ]]; then
                check=any
                # The 49 GiB collate output exceeds 8 GiB of input per rank.
                validation_capacity=12G
            fi
            run validation "$app" "$repeat" "$validation_capacity" "$check" \
                "$APP" stats --basic -i "$output"
        fi
    done
    round="$ROOT/out/workflow_r${repeat}"
    mkdir -p "$round"
    run workflow collate "$repeat" "$POOL" external "$APP" collate -i "$INPUT" \
        -o "$round/name.bam" -n 66 -m "$MEMORY" --compress-level 1
    [[ -s "$round/name.bam" ]]
    run workflow fixmate "$repeat" 24G write "$APP" fixmate -i "$round/name.bam" \
        -o "$round/fixmate.bam" -m --compress-level 1
    [[ -s "$round/fixmate.bam" ]]
    run workflow sort "$repeat" "$POOL" external "$APP" sort -i "$round/fixmate.bam" \
        -o "$round/sorted.bam" -m "$MEMORY" --compress-level 1
    [[ -s "$round/sorted.bam" ]]
    run validation workflow_sort "$repeat" 8G sorted "$APP" stats --basic -i "$round/sorted.bam"
    run workflow markdup "$repeat" 24G markdup "$MARKDUP_APP" markdup -i "$round/sorted.bam" \
        -o "$round/markdup.bam" --stream -l 300 -m 8G --compress-level 1
    [[ -s "$round/markdup.bam" ]]
    run validation workflow_markdup "$repeat" 8G duplicates "$APP" stats --basic -i "$round/markdup.bam"
    if ! awk -F '\t' -v r="$repeat" '$1==r{found=1}END{exit !found}' workflows.tsv; then
        awk -F '\t' -v r="$repeat" '
          $1=="workflow" && $3==r{t[$2]=$7+0;p[$2]=$12+0;n++}
          END{if(n!=4)exit 1;printf "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",r,t["collate"],t["fixmate"],t["sort"],t["markdup"],t["collate"]+t["fixmate"]+t["sort"]+t["markdup"],p["collate"],p["fixmate"],p["sort"],p["markdup"],p["collate"]+p["fixmate"]+p["sort"]+p["markdup"]}' results.tsv >> workflows.tsv
    fi
done
awk -F '\t' -v repeats="$REPEATS" 'NR>1{t[++n]=$6+0;p[n]=$11+0}END{if(n!=repeats)exit 1;for(i=1;i<=n;i++)for(j=i+1;j<=n;j++){if(t[i]>t[j]){v=t[i];t[i]=t[j];t[j]=v};if(p[i]>p[j]){v=p[i];p[i]=p[j];p[j]=v}};m=int((n+1)/2);printf "workflow_e2e_median_s=%.6f\nworkflow_processing_median_s=%.6f\n",t[m],p[m]}' workflows.tsv > workflow_summary.txt
printf 'COMPLETE\n' > state
printf '%s D2 EXTERNAL MATRIX COMPLETED\n' "$(date -Is)"
