#!/usr/bin/env bash
set -euo pipefail

# Run on SWLS after copying an independent binary and this script into ROOT.
ROOT=${1:?Usage: markdup_ready_probe.sh ABS_ROOT small|D1|D2|D2_validate INPUT_BAM}
PHASE=${2:?Phase required}
INPUT=${3:?Coordinate-sorted fixmate BAM required}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
APP="$ROOT/bin/RabbitBAM-MPI"
export SWBAM_NODE_MEMORY=${SWBAM_NODE_MEMORY:-88000M}
[[ "$ROOT" = /* && "$INPUT" = /* && -s "$INPUT" && -x "$APP" ]]
case "$PHASE" in
    small) RANKS=2; CAPACITY=8G; RECORDS=2668351; MARKED=15695 ;;
    D1) RANKS=6; CAPACITY=8G; RECORDS=177798571; MARKED=50084605 ;;
    D2) RANKS=6; CAPACITY=24G; RECORDS=365079492; MARKED=8894488 ;;
    D2_validate) RANKS=6; CAPACITY=8G; RECORDS=365079492; MARKED=8894488 ;;
    *) printf 'Unsupported phase: %s\n' "$PHASE" >&2; exit 1 ;;
esac
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions"
cd "$ROOT"
mkdir driver.lock
CURRENT=initialization
trap 'rc=$?; if ((rc)); then printf "FAILED %s\n" "$CURRENT" > "$PHASE.state"; fi; rmdir driver.lock' EXIT
if [[ -e binaries.sha256 ]]; then
    sha256sum -c binaries.sha256
else
    sha256sum bin/RabbitBAM-MPI > binaries.sha256
fi
if [[ "$PHASE" != small ]]; then
    # The caller writes this only after independent small-output verification.
    [[ -s small.verified.sha256 ]]
    cmp binaries.sha256 small.verified.sha256
fi
[[ ! -e "$PHASE.results.tsv" ]]
stat -c '%n\t%s bytes\tmtime=%Y' "$INPUT" > "$PHASE.input.txt"
printf 'label\tjob_id\tranks\tstorage\te2e_s\tcore_s\tprepare_s\tpool_peak_bytes\toutput_bytes\tprocessing_s\n' > "$PHASE.results.tsv"

run() {
    local label=$1 storage=$2
    shift 2
    CURRENT=$label
    local log="$ROOT/logs/$label.log" submission job status unknown=0 metrics
    [[ ! -e "$log" && ! -e "$ROOT/sessions/$label" ]]
    local command=(bsub -N "$RANKS" -np 1 -health 1 -mpecg 6 -cgsp 64
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
        -a mpi -b -K -q q_share -timelimit 03:00:00
        -J "mdready_$label" -o "$log")
    if [[ -n ${PROBE_NODES:-} ]]; then command+=(-node "$PROBE_NODES"); fi
    if [[ ${PROBE_DIAGNOSTICS:-0} == 1 ]]; then
        command+=(-x SWBAM_MARKDUP_PENDING_DIAG=1)
    fi
    if [[ "$PHASE" == D2_validate ]]; then
        command+=("$APP" stats --basic -i "$INPUT" --storage "$storage" "$@")
    else
        command+=("$APP" markdup -i "$INPUT" -o "$ROOT/out/$label.bam"
            --stream -l 300 -m 8G --compress-level 1 --storage "$storage" "$@")
    fi
    if [[ "$storage" == ram ]]; then
        command+=(--ram-capacity "$CAPACITY" --ram-ranks "$RANKS"
            --ram-session "$ROOT/sessions/$label" --ram-ready-timeout 5400)
    fi
    printf '%q ' "${command[@]}" >> commands.sh
    printf '\n' >> commands.sh
    printf 'RUNNING %s\n' "$label" > "$PHASE.state"
    printf '%s START %s\n' "$(date -Is)" "$label"
    submission=$("${command[@]}")
    printf '%s\n' "$submission"
    job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
    [[ "$job" =~ ^[0-9]+$ ]]
    printf '%s\t%s\n' "$label" "$job" >> jobs.tsv
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
            *) printf 'Scheduler failure: job=%s status=%s\n' "$job" "$status" >&2; return 1 ;;
        esac
    done
    for ((i=0; i<40; ++i)); do
        if [[ -s "$log" ]] &&
            [[ $(awk '/\[swbam-processing rank=/{n++} END{print n+0}' "$log") == "$RANKS" ]] &&
            [[ $(awk '/\[swbam-e2e rank=/{n++} END{print n+0}' "$log") == "$RANKS" ]]; then break; fi
        sleep 3
    done
    ! grep -Eq '(ERROR:|RAM preparation failed|Segmentation fault|MPI_ABORT)' "$log"
    if [[ "$storage" == ram ]]; then
        local check=markdup
        if [[ "$PHASE" == D2_validate ]]; then check=duplicates; fi
        metrics=$(awk -v ranks="$RANKS" -v records="$RECORDS" -v marked="$MARKED" \
            -v check="$check" -f "$SCRIPT_DIR/ram_metrics.awk" "$log")
    else
        metrics=$(awk -v ranks="$RANKS" -v records="$RECORDS" -v marked="$MARKED" '
            /\[swbam-(e2e|processing) rank=/ {
                kind=$1; rank=$2; sub(/^rank=/,"",rank); sub(/\]$/,"",rank)
                if(rank !~ /^[0-9]+$/ || rank+0>=ranks || seen[kind,rank]++) bad=1
                ok=0; seconds=0
                for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="status"&&a[2]==0)ok=1;if(a[1]=="seconds")seconds=a[2]+0}
                if(!ok||seconds<=0)bad=1
                count[kind]++; if(seconds>max[kind])max[kind]=seconds
            }
            /^Complete the 4[.]3 / {core=$NF+0}
            /^OptimizedMarkdupMPI finished[.]/ {
                for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="records"&&a[2]==records)got_records=1;if(a[1]=="ranks"&&a[2]==ranks)got_ranks=1}
            }
            /^[[:space:]]*duplicates / {for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="marked"&&a[2]==marked)got_marked=1}}
            END {
                if(bad||count["[swbam-e2e"]!=ranks||count["[swbam-processing"]!=ranks||!got_records||!got_ranks||!got_marked||core<=0)exit 1
                printf "%.6f\t%.6f\t0\t0\t0\t%.6f",max["[swbam-e2e"],core,max["[swbam-processing"]
            }' "$log")
    fi
    if [[ "$PHASE" == small ]]; then [[ -s "$ROOT/out/$label.bam" ]]; fi
    if [[ "$PHASE" == D2 && ${PROBE_DUMP:-0} == 1 ]]; then
        [[ -s "$ROOT/out/$label.bam" ]]
        grep -q 'RAM output saved after MPI_Finalize:' "$log"
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$job" "$RANKS" "$storage" "$metrics" >> "$PHASE.results.tsv"
    printf '%s DONE %s metrics=%s\n' "$(date -Is)" "$label" "$metrics"
}

if [[ "$PHASE" == small ]]; then
    run small_ram ram
    run small_posix posix
    printf 'OUTPUTS_READY_FOR_INDEPENDENT_VERIFICATION\n' > "$PHASE.state"
elif [[ "$PHASE" == D2_validate ]]; then
    run D2_validate ram
    printf 'COMPLETE\n' > "$PHASE.state"
else
    # Complete output remains in the RAM pool; skip only the post-Finalize dump.
    if [[ ${PROBE_DUMP:-0} == 1 ]]; then
        run "${PHASE}_ram" ram
    else
        run "${PHASE}_ram" ram --no-dump
    fi
    printf 'COMPLETE\n' > "$PHASE.state"
fi
