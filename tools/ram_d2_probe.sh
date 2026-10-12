#!/usr/bin/env bash
set -euo pipefail

# One capacity/functional probe per command, not a performance matrix.
# Run on SWLS with frozen binaries in ROOT/bin and existing intermediate inputs.
ROOT=${1:?Usage: ram_d2_probe.sh ABS_ROOT BAM SAM COLLATED_BAM FIXMATE_BAM SORTED_BAM}
INPUT=${2:?D2 BAM required}
SAM=${3:?Existing full D2 SAM required}
COLLATED=${4:?Existing D2 name-collated BAM required}
FIXMATE=${5:?Existing D2 fixmate BAM required}
SORTED=${6:?Existing coordinate-sorted fixmate BAM required}
APP="$ROOT/bin/RabbitBAM-MPI"
RANKS=6
RECORDS=365079492
NODE_RANGE=${NODE_RANGE:-auto}
READY_TIMEOUT=${READY_TIMEOUT:-5400}
CURRENT=initialization
CURRENT_JOB=
[[ "$ROOT" = /* && "$READY_TIMEOUT" =~ ^[1-9][0-9]*$ ]]
for input in "$INPUT" "$SAM" "$COLLATED" "$FIXMATE" "$SORTED"; do
    [[ "$input" = /* && -s "$input" ]] || { printf 'Invalid input: %s\n' "$input" >&2; exit 1; }
done
for binary in "$APP" "$ROOT/bin/swbam-sdk-optimized-count" "$ROOT/bin/swbam-sdk-optimized-bam2bam"; do
    [[ -x "$binary" ]] || exit 1
done
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions"
mkdir "$ROOT/probe.lock"
trap 'rc=$?; if ((rc)); then printf "FAILED %s job=%s\n" "$CURRENT" "$CURRENT_JOB" > "$ROOT/state"; fi; rmdir "$ROOT/probe.lock"' EXIT
cd "$ROOT"
if [[ ${PROBE_RESUME:-0} == 1 ]]; then
    [[ -s results.tsv && -s commands.sh && -s jobs.tsv ]]
    sha256sum -c binaries.sha256
else
    [[ ! -e results.tsv && ! -e commands.sh ]] || exit 1
    printf 'name\tjob_id\tcapacity\te2e_max_s\tprepare_max_s\tpool_peak_max_bytes\toutput_bytes\n' > results.tsv
    sha256sum bin/RabbitBAM-MPI bin/swbam-sdk-optimized-count bin/swbam-sdk-optimized-bam2bam > binaries.sha256
    printf 'ranks=%s node_range=%s ready_timeout=%s\n' "$RANKS" "$NODE_RANGE" "$READY_TIMEOUT" > configuration.txt
    stat -c '%n\t%s bytes' "$INPUT" "$SAM" "$COLLATED" "$FIXMATE" "$SORTED" >> configuration.txt
fi
NODE_ARGS=(-health 1)
if [[ "$NODE_RANGE" != auto ]]; then
    [[ "$NODE_RANGE" =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]] || exit 1
    NODE_ARGS=(-node "$NODE_RANGE" -health 1)
fi

run() {
    local name=$1 capacity=$2 expected=$3 mode=$4
    shift 4
    local log="$ROOT/logs/$name.log" submission job metrics scheduler_status=
    if [[ ${PROBE_RESUME:-0} == 1 ]]; then
        if awk -F '\t' -v name="$name" '$1==name{found=1}END{exit !found}' results.tsv; then
            printf 'REUSE verified %s\n' "$name"
            return 0
        fi
        job=$(awk -F '\t' -v name="$name" '$1==name{print $2}' jobs.tsv)
    else
        job=
    fi
    CURRENT=$name
    printf 'RUNNING %s\n' "$name" > "$ROOT/state"
    if [[ -z "$job" ]]; then
        [[ ! -e "$log" && ! -e "$ROOT/sessions/$name" ]] || return 1
        printf '%s START %s capacity=%s\n' "$(date -Is)" "$name" "$capacity"
        local command=(bsub -N "$RANKS" -np 1 "${NODE_ARGS[@]}" -mpecg 6 -cgsp 64
            -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
            -a mpi -b -K -q q_share -timelimit 03:00:00 -J "ram_D2_$name" -o "$log"
            "$@" --storage ram --ram-capacity "$capacity" --ram-ranks "$RANKS"
            --ram-session "$ROOT/sessions/$name" --ram-ready-timeout "$READY_TIMEOUT")
        printf '%q ' "${command[@]}" >> commands.sh
        printf '\n' >> commands.sh
        submission=$("${command[@]}")
        printf '%s\n' "$submission"
        job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
        [[ "$job" =~ ^[0-9]+$ ]] || return 1
        printf '%s\t%s\n' "$name" "$job" >> jobs.tsv
    fi
    [[ "$job" =~ ^[0-9]+$ ]] || return 1
    CURRENT_JOB=$job
    while true; do
        scheduler_status=$(bjobs "$job" 2>/dev/null | awk -v job="$job" '$1==job{print $2}' || true)
        case "$scheduler_status" in
            PEND|RUN|EXITING) sleep 5 ;;
            DONE|"") break ;;
            *) printf 'Unexpected scheduler state: %s\n' "$scheduler_status" >&2; return 1 ;;
        esac
    done
    for ((tries=0; tries<40; ++tries)); do
        if [[ -s "$log" ]]; then
            grep -Eq '(ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT)' "$log" && break
            if [[ "$scheduler_status" == DONE ]] || grep -q "Job $job has been finished" "$log"; then
                [[ $(awk '/\[ram-storage rank=/{n++}END{print n+0}' "$log") == "$RANKS" ]] && break
            fi
        fi
        sleep 3
    done
    [[ -s "$log" ]] || return 1
    [[ "$scheduler_status" == DONE ]] || grep -q "Job $job has been finished" "$log" || return 1
    if grep -Eq '(ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT)' "$log"; then
        tail -60 "$log"
        return 1
    fi
    [[ $(awk '/\[swbam-e2e rank=[0-9]+\].* status=0/{n++}END{print n+0}' "$log") == "$RANKS" ]] || return 1
    [[ $(awk '/\[ram-storage rank=/{n++}END{print n+0}' "$log") == "$RANKS" ]] || return 1
    grep -Eq "(total_records=|records=)$expected([[:space:]]|$)" "$log" || return 1
    if [[ "$mode" == memory ]]; then grep -q 'selected mode=memory' "$log" || return 1; fi
    if [[ "$name" == filter ]]; then grep -q 'kept_records=328059484' "$log" || return 1; fi
    if [[ "$name" == markdup ]]; then grep -Eq 'marked[=: ]+8894488([[:space:]]|$)' "$log" || return 1; fi
    if [[ "$name" == sdk_validate ]]; then
        awk -F '\t' '$1=="SN" && $2=="is sorted:" && $3==1{ok=1}END{exit !ok}' "$log" || return 1
    fi
    metrics=$(awk '
      /\[swbam-e2e rank=/{for(i=1;i<=NF;i++)if($i~/^seconds=/){split($i,a,"=");if(a[2]>e)e=a[2]}}
      /\[ram-prepare rank=/{for(i=1;i<=NF;i++)if($i~/^seconds=/){split($i,a,"=");if(a[2]>p)p=a[2]}}
      /\[ram-storage rank=/{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="peak" && a[2]>k)k=a[2];if(a[1]=="output_bytes" && a[2]>o)o=a[2]}}
      END{printf "%.6f\t%.6f\t%.0f\t%.0f",e,p,k,o}' "$log")
    printf '%s\t%s\t%s\t%s\n' "$name" "$job" "$capacity" "$metrics" >> results.tsv
    printf '%s DONE %s metrics=%s\n' "$(date -Is)" "$name" "$metrics"
    CURRENT_JOB=
}

run count 8G "$RECORDS" any "$ROOT/bin/swbam-sdk-optimized-count" "$INPUT"
run sdk_copy 16G "$RECORDS" any "$ROOT/bin/swbam-sdk-optimized-bam2bam" "$INPUT" -o "$ROOT/out/D2.sdk.bam"
[[ -s "$ROOT/out/D2.sdk.bam" ]]
run sdk_validate 8G "$RECORDS" any "$APP" stats -i "$ROOT/out/D2.sdk.bam" --basic
run flagstat 8G "$RECORDS" any "$APP" flagstat -i "$INPUT"
run filter 16G "$RECORDS" any "$APP" run_all -i "$INPUT" -o "$ROOT/out/D2.filter.bam" --min-mapq 30 --no-dump
run collate 16G "$RECORDS" memory "$APP" collate -i "$INPUT" -o "$ROOT/out/D2.collate.bam" -n 66 -m 48G --no-dump
run fixmate 24G "$RECORDS" any "$APP" fixmate -i "$COLLATED" -o "$ROOT/out/D2.fixmate.bam" -m --no-dump
run sort 24G "$RECORDS" memory "$APP" sort -i "$FIXMATE" -o "$ROOT/out/D2.sort.bam" -m 48G --no-dump
run markdup 24G "$RECORDS" any "$APP" markdup -i "$SORTED" -o "$ROOT/out/D2.markdup.bam" --stream -m 8G --no-dump
run bam2sam 40G "$RECORDS" any "$APP" run_all -i "$INPUT" -o "$ROOT/out/D2.records.sam" --no-dump
run sam2bam 40G "$RECORDS" any "$APP" run_all -i "$SAM" -o "$ROOT/out/D2.roundtrip.bam" --no-dump
for name in filter collate fixmate sort markdup records.sam roundtrip; do
    [[ ! -e "$ROOT/out/D2.$name" && ! -e "$ROOT/out/D2.$name.bam" ]] || exit 1
done
printf 'COMPLETE\n' > "$ROOT/state"
printf '%s D2 CAPACITY PROBE COMPLETED (one-off, not formal medians)\n' "$(date -Is)"
