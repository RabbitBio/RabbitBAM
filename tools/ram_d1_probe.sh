#!/usr/bin/env bash
set -euo pipefail
export SWBAM_NODE_MEMORY=${SWBAM_NODE_MEMORY:-88000M}
ROOT=${1:?Usage: ram_d1_probe.sh ABS_ROOT D1_BAM}
INPUT=${2:?D1 BAM required}
APP="$ROOT/bin/RabbitBAM-MPI"
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions"
cd "$ROOT"
run() {
    local name=$1 job submission log="$ROOT/logs/$1.log"
    shift
    [[ ! -e $log ]] || return 1
    printf '%s START %s\n' "$(date -Is)" "$name"
    submission=$(bsub -N 6 -np 1 -health 1 -mpecg 6 -cgsp 64 \
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
        -a mpi -b -K -q q_share -J "ram_D1_$name" -o "$log" "$@")
    printf '%s\n' "$submission"
    job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
    [[ -n $job ]] || return 1
    while [[ $(bjobs "$job" 2>/dev/null) =~ [[:space:]](PEND|RUN|EXITING)[[:space:]] ]]; do sleep 5; done
    for ((i=0;i<40;++i)); do [[ -s $log ]] && break; sleep 3; done
    [[ -s $log ]] && ! grep -Eq '(ERROR:|RAM preparation failed)' "$log" || { tail -40 "$log"; return 1; }
    if [[ $name == old_* ]]; then
        grep -q 'records=177798571' "$log"
    else
        [[ $(awk '/\[swbam-e2e rank=[0-9]+\].* status=0/{n++} END{print n+0}' "$log") == 6 ]] || { tail -40 "$log"; return 1; }
        awk '/\[swbam-e2e rank=/{for(i=1;i<=NF;i++)if($i~/^seconds=/){split($i,a,"=");if(a[2]>m)m=a[2]}}END{printf "MPI interval max=%.6f s\n",m}' "$log"
    fi
    printf '%s DONE %s\n' "$(date -Is)" "$name"
}
ram() {
    local name=$1 capacity=$2
    shift 2
    run "$name" "$@" --storage ram --ram-capacity "$capacity" --ram-ranks 6 \
        --ram-session "$ROOT/sessions/D1_$name" --ram-ready-timeout 900
}
ram count 2G "$ROOT/bin/swbam-sdk-optimized-count" "$INPUT"
run old_count "$ROOT/bin/sdk-count-before-ram-e2e" "$INPUT"
ram copy 8G "$ROOT/bin/swbam-sdk-optimized-bam2bam" "$INPUT" -o "$ROOT/out/D1.sdk.bam"
run old_copy "$ROOT/bin/sdk-bam2bam-before-ram-e2e" "$INPUT"
ram validate 8G "$APP" stats -i "$ROOT/out/D1.sdk.bam" --basic
ram flagstat 2G "$APP" flagstat -i "$INPUT"
ram filter 8G "$APP" run_all -i "$INPUT" -o "$ROOT/out/D1.filter.bam" --min-mapq 30 --no-dump
ram collate 8G "$APP" collate -i "$INPUT" -o "$ROOT/out/D1.collate.bam" -n 66 -m 24G --no-dump
ram sort 8G "$APP" sort -i "$INPUT" -o "$ROOT/out/D1.sort.bam" -m 24G --no-dump
printf '%s D1 PROBE COMPLETED\n' "$(date -Is)"
