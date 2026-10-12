#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: ram_storage_edges.sh ABS_ROOT SMALL_BAM}
INPUT=${2:?Small BAM required}
APP="$ROOT/bin/RabbitBAM-MPI"
COUNT="$ROOT/bin/swbam-sdk-optimized-count"
COPY="$ROOT/bin/swbam-sdk-optimized-bam2bam"
START_AT=${START_AT:-empty_count}
SESSION_TAG=${SESSION_TAG:-}
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions" "$ROOT/traces"
cd "$ROOT"
run() {
    local name=$1 expected=$2 job submission log="$ROOT/logs/edge_$1$SESSION_TAG.log"
    shift 2
    if [[ -n $START_AT ]]; then [[ $START_AT == "$name" ]] || return 0; START_AT=; fi
    [[ ! -e $log ]] || return 1
    submission=$(bsub -N 2 -np 1 -health 1 -mpecg 6 -cgsp 64 \
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
        -a mpi -b -K -q q_share -J "ram_edge_$name" -o "$log" "$@")
    printf '%s %s\n' "$name" "$submission"
    job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
    [[ -n $job ]] || return 1
    while [[ $(bjobs "$job" 2>/dev/null) =~ [[:space:]](PEND|RUN|EXITING)[[:space:]] ]]; do sleep 3; done
    for ((i=0;i<40;++i)); do [[ -s $log ]] && break; sleep 3; done
    [[ -s $log ]] || return 1
    if [[ $expected == pass ]]; then
        [[ $(awk '/\[swbam-e2e rank=[0-9]+\].* status=0/{n++}END{print n+0}' "$log") == 2 ]] || { tail -40 "$log"; return 1; }
        ! grep -Eq '(ERROR:|RAM preparation failed)' "$log"
    else
        grep -Eq '(ERROR:|RAM preparation failed|\[swbam-e2e.*status=[1-9])' "$log"
        ! grep -q 'RAM output saved' "$log"
    fi
    printf 'EDGE %s %s\n' "$name" "$expected"
}
ram() {
    local name=$1 expected=$2 capacity=$3
    shift 3
    run "$name" "$expected" "$@" --storage ram --ram-capacity "$capacity" \
        --ram-ranks 2 --ram-session "$ROOT/sessions/edge_$name$SESSION_TAG" --ram-ready-timeout 10
}
ram empty_count pass 64M "$COUNT" "$ROOT/fixtures/empty.bam"
ram empty_copy pass 64M "$COPY" "$ROOT/fixtures/empty.bam" -o "$ROOT/out/empty.bam"
ram one_count pass 64M "$COUNT" "$ROOT/fixtures/one.bam"
ram one_copy pass 64M "$COPY" "$ROOT/fixtures/one.bam" -o "$ROOT/out/one.bam"
ram one_sam pass 64M "$APP" run_all -i "$ROOT/fixtures/one.sam" -o "$ROOT/out/one_sam.bam"
ram no_dump pass 1G "$COPY" "$INPUT" -o "$ROOT/out/no_dump.bam" --no-dump
[[ ! -e $ROOT/out/no_dump.bam ]]
ram preload_capacity fail 16M "$COUNT" "$INPUT"
ram output_capacity fail 192M "$COPY" "$INPUT" -o "$ROOT/out/capacity.bam"
run missing_rank fail "$COUNT" "$ROOT/fixtures/one.bam" --storage ram --ram-capacity 64M \
    --ram-ranks 3 --ram-session "$ROOT/sessions/edge_missing_rank" --ram-ready-timeout 2
ram invalid_cli fail 64M "$APP" flagstat -i "$ROOT/fixtures/one.bam" --unknown-option
ram temporary_capacity fail 512M "$APP" collate -i "$INPUT" -o "$ROOT/out/temp_capacity.bam" -n 66 -m 128M
run reject_dedup fail "$APP" dedup-workflow -i "$INPUT" -o "$ROOT/out/rejected.bam" --storage ram \
    --ram-capacity 1G --ram-ranks 2 --ram-session "$ROOT/sessions/edge_rejected"
printf 'EDGE REGRESSION COMPLETED\n'
