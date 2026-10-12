#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Frozen RAM binaries, independent performance samples, and real four-step rounds.
ROOT=${1:?Usage: ram_matrix.sh ABS_ROOT BAM SAM COLLATED_BAM FIXMATE_BAM SORTED_BAM}
INPUT=${2:?BAM required}
SAM=${3:?Full SAM required}
COLLATED=${4:?Name-collated BAM required}
FIXMATE=${5:?Fixmate BAM required}
SORTED=${6:?Coordinate-sorted fixmate BAM required}
APP="$ROOT/bin/RabbitBAM-MPI"
COUNT="$ROOT/bin/swbam-sdk-optimized-count"
COPY="$ROOT/bin/swbam-sdk-optimized-bam2bam"
DATASET=${BENCH_DATASET:-D1}
case "$DATASET" in
    D1)
        RECORDS=177798571; KEPT=170359591; MARKED=50084605
        READ_POOL=2G; FILTER_POOL=8G; COLLATE_POOL=8G; FIXMATE_POOL=12G
        SORT_POOL=8G; MARKDUP_POOL=8G; SAM_POOL=24G; VERIFY_POOL=4G
        COLLATE_MEMORY=24G; SORT_MEMORY=32G ;;
    D2)
        RECORDS=365079492; KEPT=328059484; MARKED=8894488
        READ_POOL=8G; FILTER_POOL=16G; COLLATE_POOL=16G; FIXMATE_POOL=24G
        SORT_POOL=16G; MARKDUP_POOL=24G; SAM_POOL=40G; VERIFY_POOL=8G
        COLLATE_MEMORY=48G; SORT_MEMORY=48G ;;
    *) printf 'Unsupported dataset: %s\n' "$DATASET" >&2; exit 1 ;;
esac
STAGE=${BENCH_STAGE:-all}
[[ "$STAGE" == all || ("$DATASET" == D2 && "$STAGE" == capacity) ]]
CURRENT=initialization
CURRENT_JOB=
[[ "$ROOT" = /* ]]
for input in "$INPUT" "$SAM" "$COLLATED" "$FIXMATE" "$SORTED"; do
    [[ "$input" = /* && -s "$input" ]] || exit 1
done
for binary in "$APP" "$COUNT" "$COPY"; do [[ -x "$binary" ]] || exit 1; done
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions"
mkdir "$ROOT/driver.lock"
trap 'rc=$?; if ((rc)); then printf "FAILED %s job=%s\n" "$CURRENT" "$CURRENT_JOB" > "$ROOT/state"; fi; rmdir "$ROOT/driver.lock"' EXIT
cd "$ROOT"
if [[ ${BENCH_RESUME:-0} == 1 ]]; then
    [[ -s results.tsv && -s jobs.tsv && -s commands.sh ]]
    sha256sum -c binaries.sha256
else
    [[ ! -e results.tsv && ! -e commands.sh ]] || exit 1
    printf 'kind\tapp\trepeat\tranks\tjob_id\tcapacity\te2e_max_s\tcore_s\tprepare_max_s\tpool_peak_max_bytes\toutput_bytes\n' > results.tsv
    printf 'repeat\tcollate_s\tfixmate_s\tsort_s\tmarkdup_s\ttotal_s\n' > workflows.tsv
    sha256sum bin/RabbitBAM-MPI bin/swbam-sdk-optimized-count bin/swbam-sdk-optimized-bam2bam > binaries.sha256
    stat -c '%n\t%s bytes\tmtime=%Y' "$INPUT" "$SAM" "$COLLATED" "$FIXMATE" "$SORTED" > inputs.txt
    printf 'dataset=%s sdk_ranks=1/2/4/6/12 app_ranks=6 repeats=3 bins=66 level=1 min_mapq=30\n' "$DATASET" > configuration.txt
    printf 'sort_pool=%s sort_memory=%s collate_pool=%s collate_memory=%s baseline=%s wait_for=%s\n' \
        "$SORT_POOL" "$SORT_MEMORY" "$COLLATE_POOL" "$COLLATE_MEMORY" \
        "${BENCH_D2_BASELINE_ROOT:-none}" "${BENCH_WAIT_FOR:-none}" >> configuration.txt
fi

metrics_from_log() {
    local ranks=$1 check=$2 log=$3
    [[ -s "$log" ]] || return 1
    if grep -Eq '(ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT|exit_code=[[:space:]]*[1-9])' "$log"; then
        tail -60 "$log" >&2
        return 1
    fi
    awk -v ranks="$ranks" -v records="$RECORDS" -v kept="$KEPT" -v marked="$MARKED" \
        -v check="$check" -f "$SCRIPT_DIR/ram_metrics.awk" "$log"
}

run() {
    local kind=$1 app=$2 repeat=$3 ranks=$4 capacity=$5 check=$6
    shift 6
    local label="${kind}_${app}_n${ranks}_r${repeat}" job= status= metrics log
    CURRENT=$label
    CURRENT_JOB=
    log="$ROOT/logs/$label.log"
    if [[ ${BENCH_RESUME:-0} == 1 ]]; then
        if awk -F '\t' -v k="$kind" -v a="$app" -v r="$repeat" -v n="$ranks" \
            '$1==k && $2==a && $3==r && $4==n{found=1}END{exit !found}' results.tsv; then
            metrics_from_log "$ranks" "$check" "$log" >/dev/null
            printf 'REUSE verified %s\n' "$label"
            return 0
        fi
        job=$(awk -F '\t' -v name="$label" '$1==name{print $2}' jobs.tsv)
    fi
    printf 'RUNNING %s\n' "$label" > state
    if [[ -z "$job" ]]; then
        [[ ! -e "$log" && ! -e "$ROOT/sessions/$label" ]] || return 1
        local diagnostic_args=()
        if [[ "$DATASET:$kind:$app" == D2:probe:sort ]]; then
            diagnostic_args=(-x SWBAM_SORT_MEMORY_DIAG=1)
        fi
        local command=(bsub -N "$ranks" -np 1 -health 1 -mpecg 6 -cgsp 64
            -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
            -a mpi -b -K -q q_share -timelimit 03:00:00 -J "${DATASET,,}ram_$label" -o "$log"
            ${diagnostic_args[@]+"${diagnostic_args[@]}"}
            "$@" --storage ram --ram-capacity "$capacity" --ram-ranks "$ranks"
            --ram-session "$ROOT/sessions/$label" --ram-ready-timeout 5400)
        printf '%q ' "${command[@]}" >> commands.sh
        printf '\n' >> commands.sh
        printf '%s START %s capacity=%s\n' "$(date -Is)" "$label" "$capacity"
        local submission
        submission=$("${command[@]}")
        printf '%s\n' "$submission"
        job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
        [[ "$job" =~ ^[0-9]+$ ]] || return 1
        printf '%s\t%s\n' "$label" "$job" >> jobs.tsv
    fi
    [[ "$job" =~ ^[0-9]+$ ]] || return 1
    CURRENT_JOB=$job
    local unknown=0
    while true; do
        status=$(bjobs "$job" 2>/dev/null | awk -v job="$job" '$1==job{print $2}' || true)
        case "$status" in
            PEND|RUN|EXITING) unknown=0; sleep 5 ;;
            DONE) break ;;
            "")
                if [[ -s "$log" ]] && grep -q "Job $job has been finished" "$log"; then break; fi
                unknown=$((unknown + 1))
                ((unknown < 40)) || return 1
                sleep 5 ;;
            *) printf 'Scheduler failure: %s job=%s\n' "$status" "$job" >&2; return 1 ;;
        esac
    done
    for ((tries=0; tries<40; ++tries)); do
        if [[ -s "$log" ]]; then
            grep -Eq '(ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT)' "$log" && break
            [[ $(awk '/\[ram-storage rank=/{n++}END{print n+0}' "$log") == "$ranks" ]] && break
        fi
        sleep 3
    done
    metrics=$(metrics_from_log "$ranks" "$check" "$log") || {
        printf 'Invalid metrics: %s\n' "$label" >&2; return 1;
    }
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$kind" "$app" "$repeat" "$ranks" "$job" "$capacity" "$metrics" >> results.tsv
    printf '%s DONE %s job=%s metrics=%s\n' "$(date -Is)" "$label" "$job" "$metrics"
    CURRENT_JOB=
}

if [[ "$DATASET" == D2 ]]; then
    CURRENT=baseline_gate
    baseline=${BENCH_D2_BASELINE_ROOT:?D2 requires the prior successful six-rank probe directory}
    [[ "$baseline" = /* && "$baseline" != "$ROOT" ]]
    cmp binaries.sha256 "$baseline/binaries.sha256"
    diff <(tail -n +2 "$baseline/configuration.txt") \
         <(stat -c '%n\t%s bytes' "$INPUT" "$SAM" "$COLLATED" "$FIXMATE" "$SORTED")
    printf 'app\tjob_id\tverified_log\te2e_max_s\tcore_s\tprepare_max_s\tpool_peak_max_bytes\toutput_bytes\n' > prior_probes.tsv
    for name in count sdk_copy sdk_validate flagstat filter collate fixmate; do
        check=any
        case "$name" in
            sdk_copy|fixmate) check=write ;;
            sdk_validate) check=sorted ;;
            filter) check=filter ;;
            collate) check=memory ;;
        esac
        metrics=$(metrics_from_log 6 "$check" "$baseline/logs/$name.log")
        job=$(awk -F '\t' -v name="$name" '$1==name{print $2}' "$baseline/results.tsv")
        [[ "$job" =~ ^[0-9]+$ ]]
        printf '%s\t%s\t%s\t%s\n' "$name" "$job" "$baseline/logs/$name.log" "$metrics" >> prior_probes.tsv
    done

    run probe sort 1 6 "$SORT_POOL" memory "$APP" sort \
        -i "$FIXMATE" -o "$ROOT/out/probe.sort.bam" -m "$SORT_MEMORY" --compress-level 1
    [[ -s "$ROOT/out/probe.sort.bam" ]]
    run probe_validation sort 1 6 "$VERIFY_POOL" sorted "$APP" stats --basic -i "$ROOT/out/probe.sort.bam"
    run probe markdup 1 6 "$MARKDUP_POOL" markdup "$APP" markdup -i "$SORTED" \
        -o "$ROOT/out/probe.markdup.bam" --stream -l 300 -m 8G --compress-level 1
    [[ -s "$ROOT/out/probe.markdup.bam" ]]
    run probe_validation markdup 1 6 "$VERIFY_POOL" duplicates "$APP" stats --basic -i "$ROOT/out/probe.markdup.bam"
    run probe bam2sam 1 6 "$SAM_POOL" write "$APP" run_all -i "$INPUT" \
        -o "$ROOT/out/probe.sam" --compress-level 1 --no-dump
    run probe sam2bam 1 6 "$SAM_POOL" write "$APP" run_all -i "$SAM" \
        -o "$ROOT/out/probe.sam2bam.bam" --compress-level 1
    [[ -s "$ROOT/out/probe.sam2bam.bam" ]]
    run probe_validation sam2bam 1 6 "$VERIFY_POOL" sorted "$APP" stats --basic -i "$ROOT/out/probe.sam2bam.bam"
    run probe count 1 1 36G any "$COUNT" "$INPUT"
    run probe bam2bam 1 1 76G write "$COPY" "$INPUT" -o "$ROOT/out/probe.sdk.bam"
    [[ -s "$ROOT/out/probe.sdk.bam" ]]
    run probe_validation sdk 1 6 "$VERIFY_POOL" sorted "$APP" stats --basic -i "$ROOT/out/probe.sdk.bam"
    printf '%s CAPACITY AND CORRECTNESS GATES PASSED\n' "$(date -Is)" > capacity_gate.txt
    if [[ "$STAGE" == capacity ]]; then printf 'CAPACITY_COMPLETE\n' > state; exit 0; fi
fi

if [[ -n ${BENCH_WAIT_FOR:-} ]]; then
    CURRENT=waiting_for_prior_matrix
    [[ "$BENCH_WAIT_FOR" = /* && "$BENCH_WAIT_FOR" != "$ROOT" && -s "$BENCH_WAIT_FOR/state" ]]
    printf 'WAITING_FOR %s\n' "$BENCH_WAIT_FOR" > state
    wait_start=$SECONDS
    while [[ $(<"$BENCH_WAIT_FOR/state") != COMPLETE ]]; do
        if grep -q '^FAILED ' "$BENCH_WAIT_FOR/state"; then
            printf 'Prior matrix failed; D2 formal matrix not submitted.\n' >&2; exit 1
        fi
        ((SECONDS - wait_start < 86400)) || exit 1
        sleep 30
    done
fi

for repeat in 1 2 3; do
    points=(1 2 4 6 12)
    if ((repeat == 2)); then points=(12 6 4 2 1); fi
    for ranks in "${points[@]}"; do
        case "$DATASET:$ranks" in
            D1:1) read_capacity=12G; write_capacity=24G ;;
            D1:2) read_capacity=6G; write_capacity=16G ;;
            D1:4) read_capacity=4G; write_capacity=8G ;;
            D1:6) read_capacity=2G; write_capacity=8G ;;
            D1:12) read_capacity=2G; write_capacity=4G ;;
            D2:1) read_capacity=36G; write_capacity=76G ;;
            D2:2) read_capacity=20G; write_capacity=40G ;;
            D2:4) read_capacity=12G; write_capacity=24G ;;
            D2:6) read_capacity=8G; write_capacity=16G ;;
            D2:12) read_capacity=4G; write_capacity=8G ;;
        esac
        run sdk count "$repeat" "$ranks" "$read_capacity" any "$COUNT" "$INPUT"
        save=(--no-dump)
        if ((repeat == 1 && ranks == 1)); then save=(); fi
        run sdk bam2bam "$repeat" "$ranks" "$write_capacity" write "$COPY" "$INPUT" \
            -o "$ROOT/out/sdk_n${ranks}_r${repeat}.bam" ${save[@]+"${save[@]}"}
        if ((repeat == 1 && ranks == 1)); then
            [[ -s "$ROOT/out/sdk_n1_r1.bam" ]]
            run validation sdk "$repeat" 6 "$VERIFY_POOL" sorted "$APP" stats --basic -i "$ROOT/out/sdk_n1_r1.bam"
        fi
    done

    common=(--compress-level 1 --no-dump)
    run app flagstat "$repeat" 6 "$READ_POOL" any "$APP" flagstat -i "$INPUT"
    run app filter "$repeat" 6 "$FILTER_POOL" filter "$APP" run_all -i "$INPUT" \
        -o "$ROOT/out/filter_r${repeat}.bam" --min-mapq 30 "${common[@]}"
    run app collate "$repeat" 6 "$COLLATE_POOL" memory "$APP" collate -i "$INPUT" \
        -o "$ROOT/out/collate_r${repeat}.bam" -n 66 -m "$COLLATE_MEMORY" "${common[@]}"
    run app fixmate "$repeat" 6 "$FIXMATE_POOL" write "$APP" fixmate -i "$COLLATED" \
        -o "$ROOT/out/fixmate_r${repeat}.bam" -m "${common[@]}"
    run app sort "$repeat" 6 "$SORT_POOL" memory "$APP" sort -i "$FIXMATE" \
        -o "$ROOT/out/sort_r${repeat}.bam" -m "$SORT_MEMORY" "${common[@]}"
    run app markdup "$repeat" 6 "$MARKDUP_POOL" markdup "$APP" markdup -i "$SORTED" \
        -o "$ROOT/out/markdup_r${repeat}.bam" --stream -l 300 -m 8G "${common[@]}"
    run app bam2sam "$repeat" 6 "$SAM_POOL" write "$APP" run_all -i "$INPUT" \
        -o "$ROOT/out/bam2sam_r${repeat}.sam" "${common[@]}"
    save=(--no-dump)
    if ((repeat == 1)); then save=(); fi
    run app sam2bam "$repeat" 6 "$SAM_POOL" write "$APP" run_all -i "$SAM" \
        -o "$ROOT/out/sam2bam_r${repeat}.bam" --compress-level 1 ${save[@]+"${save[@]}"}
    if ((repeat == 1)); then
        [[ -s "$ROOT/out/sam2bam_r1.bam" ]]
        run validation sam2bam "$repeat" 6 "$VERIFY_POOL" sorted "$APP" stats --basic -i "$ROOT/out/sam2bam_r1.bam"
    fi

    round="$ROOT/out/workflow_r${repeat}"
    mkdir -p "$round"
    run workflow collate "$repeat" 6 "$COLLATE_POOL" memory "$APP" collate -i "$INPUT" \
        -o "$round/name.bam" -n 66 -m "$COLLATE_MEMORY" --compress-level 1
    [[ -s "$round/name.bam" ]]
    run workflow fixmate "$repeat" 6 "$FIXMATE_POOL" write "$APP" fixmate -i "$round/name.bam" \
        -o "$round/fixmate.bam" -m --compress-level 1
    [[ -s "$round/fixmate.bam" ]]
    run workflow sort "$repeat" 6 "$SORT_POOL" memory "$APP" sort -i "$round/fixmate.bam" \
        -o "$round/fixmate.sorted.bam" -m "$SORT_MEMORY" --compress-level 1
    [[ -s "$round/fixmate.sorted.bam" ]]
    run validation workflow_sort "$repeat" 6 "$VERIFY_POOL" sorted "$APP" stats --basic -i "$round/fixmate.sorted.bam"
    run workflow markdup "$repeat" 6 "$MARKDUP_POOL" markdup "$APP" markdup -i "$round/fixmate.sorted.bam" \
        -o "$round/markdup.bam" --stream -l 300 -m 8G --compress-level 1
    [[ -s "$round/markdup.bam" ]]
    run validation workflow_markdup "$repeat" 6 "$VERIFY_POOL" duplicates "$APP" stats --basic -i "$round/markdup.bam"
    if ! awk -F '\t' -v r="$repeat" '$1==r{found=1}END{exit !found}' workflows.tsv; then
        awk -F '\t' -v r="$repeat" '
          $1=="workflow" && $3==r{t[$2]=$7+0;n++}
          END{if(n!=4)exit 1;printf "%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",r,t["collate"],t["fixmate"],t["sort"],t["markdup"],t["collate"]+t["fixmate"]+t["sort"]+t["markdup"]}' \
          results.tsv >> workflows.tsv
    fi
done
awk -F '\t' 'NR>1{t[++n]=$6+0}END{if(n!=3)exit 1;for(i=1;i<=n;i++)for(j=i+1;j<=n;j++)if(t[i]>t[j]){v=t[i];t[i]=t[j];t[j]=v};printf "workflow_e2e_median_s=%.6f\n",t[2]}' workflows.tsv > workflow_summary.txt
printf 'COMPLETE\n' > state
printf '%s %s RAM MATRIX COMPLETED\n' "$(date -Is)" "$DATASET"
