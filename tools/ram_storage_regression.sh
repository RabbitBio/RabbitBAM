#!/usr/bin/env bash
set -euo pipefail
export SWBAM_NODE_MEMORY=${SWBAM_NODE_MEMORY:-88000M}

# Run on the SWLS login host, with binaries already copied into ROOT/bin.
ROOT=${1:?Usage: ram_storage_regression.sh ABS_ROOT INPUT_BAM [RANKS]}
INPUT=${2:?Input BAM required}
RANKS=${3:-2}
CAPACITY=${RAM_CAPACITY:-8G}
START_AT=${START_AT:-count}
SESSION_TAG=${SESSION_TAG:-}
APP="$ROOT/bin/RabbitBAM-MPI"
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions"
cd "$ROOT"
run() {
    local name=$1
    shift
    if [[ -n $START_AT ]]; then
        [[ $START_AT == "$name" ]] || return 0
        START_AT=
    fi
    printf '%s START %s\n' "$(date -Is)" "$name"
    local submission job log="$ROOT/logs/$name$SESSION_TAG.log"
    [[ ! -e $log ]] || { printf 'Refusing to reuse log: %s\n' "$log" >&2; return 1; }
    submission=$(bsub -N "$RANKS" -np 1 -health 1 -mpecg 6 -cgsp 64 \
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
        -a mpi -b -K -q q_share -J "ram_${name}" -o "$log" \
        "$@")
    printf '%s\n' "$submission"
    job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
    [[ -n $job ]] || return 1
    # SWLS -K controls launching, but does not block the submitting shell.
    while [[ $(bjobs "$job" 2>/dev/null) =~ [[:space:]](PEND|RUN|EXITING)[[:space:]] ]]; do
        sleep 5
    done
    for ((tries=0; tries<40; ++tries)); do
        [[ -s $log ]] && break
        sleep 3
    done
    [[ -s $log ]] || return 1
    [[ $(awk '/\[swbam-e2e rank=[0-9]+\].* status=0/{n++} END{print n+0}' "$log") == "$RANKS" ]] || {
        tail -60 "$log"
        return 1
    }
    if grep -Eq '(ERROR:|RAM preparation failed)' "$log"; then return 1; fi
    awk '/\[swbam-e2e rank=/{for(i=1;i<=NF;i++) if($i~/^seconds=/){split($i,a,"="); if(a[2]>max)max=a[2]}} END{printf "MPI interval max=%.6f s\n",max}' "$log"
    printf '%s DONE %s\n' "$(date -Is)" "$name"
}
ram() {
    local name=$1
    shift
    run "$name" "$@" --storage ram --ram-capacity "$CAPACITY" \
        --ram-ranks "$RANKS" --ram-session "$ROOT/sessions/$name$SESSION_TAG" \
        --ram-ready-timeout 300
}

ram count "$ROOT/bin/swbam-sdk-optimized-count" "$INPUT"
ram sdk_copy "$ROOT/bin/swbam-sdk-optimized-bam2bam" "$INPUT" -o "$ROOT/out/sdk.bam"
ram sdk_validate "$APP" stats -i "$ROOT/out/sdk.bam" --basic
ram io_check "$APP" io-check -i "$INPUT"
ram filter "$APP" run_all -i "$INPUT" -o "$ROOT/out/filter.bam" --min-mapq 30
ram filter_validate "$APP" stats -i "$ROOT/out/filter.bam" --basic
ram bam2sam "$APP" run_all -i "$INPUT" -o "$ROOT/out/records.sam"
ram sam2bam "$APP" run_all -i "$ROOT/out/records.sam" -o "$ROOT/out/roundtrip.bam"
ram roundtrip_validate "$APP" stats -i "$ROOT/out/roundtrip.bam" --basic
ram flagstat "$APP" flagstat -i "$INPUT"

# Force multiple runs and actual RAM temporary file traffic.
ram collate "$APP" collate -i "$INPUT" -o "$ROOT/out/collate.bam" -n 66 -m 128M
ram fixmate "$APP" fixmate -i "$ROOT/out/collate.bam" -o "$ROOT/out/fixmate.bam" -m
ram sort "$APP" sort -i "$ROOT/out/fixmate.bam" -o "$ROOT/out/sort.bam" -m 128M
ram markdup "$APP" markdup -i "$ROOT/out/sort.bam" -o "$ROOT/out/markdup.bam" --stream
ram final_validate "$APP" stats -i "$ROOT/out/markdup.bam" --basic
ram final_flagstat "$APP" flagstat -i "$ROOT/out/markdup.bam"

run native_count "$ROOT/bin/swbam-sdk-optimized-count" "$INPUT"
run native_copy "$ROOT/bin/swbam-sdk-optimized-bam2bam" "$INPUT" -o "$ROOT/out/native.bam"
run native_workflow "$APP" dedup-workflow -i "$INPUT" -o "$ROOT/out/workflow.bam" -n 66 -m 2G --temp-prefix "$ROOT/out/native-workflow"
run native_pipeline "$APP" dedup-pipeline -i "$INPUT" -o "$ROOT/out/pipeline.bam" -n 66 -m 2G
run native_workflow_validate "$APP" stats -i "$ROOT/out/workflow.bam" --basic
run native_pipeline_validate "$APP" stats -i "$ROOT/out/pipeline.bam" --basic

printf '%s REGRESSION COMPLETED\n' "$(date -Is)"
