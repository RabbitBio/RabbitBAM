#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:?Usage: submit_d3_markdup_check.sh ABS_ROOT diag|small|full|validate INPUT_BAM}
PHASE=${2:?Phase required}
INPUT=${3:?Input BAM required}
APP="$ROOT/bin/RabbitBAM-MPI"
[[ "$ROOT" = /* && "$INPUT" = /* && -x "$APP" && -s "$INPUT" ]]
case "$PHASE" in
    diag) CAPACITY=24G ;;
    small) CAPACITY=2G ;;
    full) CAPACITY=24G ;;
    validate) CAPACITY=8G ;;
    *) printf 'Unsupported phase: %s\n' "$PHASE" >&2; exit 1 ;;
esac
cd "$ROOT"
[[ ! -e job_id && ! -e "sessions/$PHASE" && ! -e "logs/$PHASE.log" ]]
mkdir -p logs sessions out
export SWBAM_NODE_MEMORY=88000M
sha256sum bin/RabbitBAM-MPI > binaries.sha256
sha256sum bin/submit_d3_markdup_check.sh > scripts.sha256
stat -c '%n %s bytes mtime=%Y' "$INPUT" > input.txt
command=(bsub -N 6 -np 1 -health 1 -mpecg 6 -cgsp 64
    -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
    -a mpi -b -K -q q_share -timelimit 02:00:00
    -J "d3_markdup_$PHASE" -o "$ROOT/logs/$PHASE.log")
if [[ "$PHASE" == diag ]]; then
    command+=(-x SWBAM_MARKDUP_PENDING_DIAG=1)
fi
if [[ "$PHASE" == validate ]]; then
    command+=("$APP" stats --basic -i "$INPUT")
else
    command+=("$APP" markdup -i "$INPUT" -o "$ROOT/out/markdup.bam"
        --stream -l 300 -m 8G --compress-level 1)
    if [[ "$PHASE" == diag ]]; then command+=(--no-dump); fi
fi
command+=(--storage ram --ram-capacity "$CAPACITY" --ram-ranks 6
    --ram-session "$ROOT/sessions/$PHASE" --ram-ready-timeout 5400)
printf '%q ' "${command[@]}" > command.sh
printf '\n' >> command.sh
submission=$("${command[@]}")
printf '%s\n' "$submission"
job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
[[ "$job" =~ ^[0-9]+$ ]]
printf '%s\n' "$job" > job_id
date -Is > started.txt
