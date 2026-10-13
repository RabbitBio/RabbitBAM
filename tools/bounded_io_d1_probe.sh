#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: bounded_io_d1_probe.sh ROOT OLD_ROOT BAM COLLATED SORTED}
OLD=${2:?Frozen previous benchmark root required}
INPUT=${3:?D1 BAM required}
COLLATED=${4:?D1 name-collated BAM required}
SORTED=${5:?D1 coordinate-sorted fixmate BAM required}
export SWBAM_NODE_MEMORY=${SWBAM_NODE_MEMORY:-88000M}
mkdir -p "$ROOT/logs" "$ROOT/out" "$ROOT/sessions"
cd "$ROOT"
run() {
    local label=$1 capacity=$2
    shift 2
    local log="$ROOT/logs/$label.log" submission job
    [[ ! -e $log ]]
    printf '%s START %s\n' "$(date -Is)" "$label"
    local nodes=(-health 1)
    if [[ -n ${PROBE_NODES:-} ]]; then nodes+=(-node "$PROBE_NODES"); fi
    submission=$(bsub -N 6 -np 1 -mpecg 6 -cgsp 64 \
        -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc \
        -a mpi -b -K -q q_share "${nodes[@]}" -J "bio_$label" -o "$log" \
        "$@" --storage ram --ram-capacity "$capacity" --ram-ranks 6 \
        --ram-session "$ROOT/sessions/$label" --no-dump)
    printf '%s\n' "$submission"
    job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
    [[ -n $job ]]
    printf '%s\t%s\n' "$label" "$job" >> "$ROOT/jobs.tsv"
    while [[ $(bjobs "$job" 2>/dev/null) =~ [[:space:]](PEND|RUN|EXITING)[[:space:]] ]]; do sleep 5; done
    for ((i=0; i<40; ++i)); do [[ -s $log ]] && break; sleep 3; done
    [[ $(awk '/\[swbam-e2e rank=[0-9]+\].*status=0/{n++} END{print n+0}' "$log") == 6 ]]
    ! grep -Eq '(ERROR:|Segmentation fault|MPI_ABORT)' "$log"
    awk -v label="$label" '/\[swbam-e2e rank=/{for(i=1;i<=NF;i++)if($i~/^seconds=/){split($i,a,"=");if(a[2]>m)m=a[2]}} END{printf "%s\t%.6f\n",label,m}' "$log" >> "$ROOT/e2e.tsv"
    printf '%s DONE %s\n' "$(date -Is)" "$label"
}
run new_count 2G "$ROOT/bin/swbam-sdk-optimized-count" "$INPUT"
run old_copy 8G "$OLD/bin/swbam-sdk-optimized-bam2bam" "$INPUT" -o "$ROOT/out/old-copy.bam"
run new_copy 8G "$ROOT/bin/swbam-sdk-optimized-bam2bam" "$INPUT" -o "$ROOT/out/new-copy.bam"
run old_fixmate 12G "$OLD/bin/RabbitBAM-MPI" fixmate -i "$COLLATED" -o "$ROOT/out/old-fixmate.bam" -m
run new_fixmate 12G "$ROOT/bin/RabbitBAM-MPI" fixmate -i "$COLLATED" -o "$ROOT/out/new-fixmate.bam" -m
run old_markdup 8G "$OLD/bin/RabbitBAM-MPI" markdup -i "$SORTED" -o "$ROOT/out/old-markdup.bam" --stream -m 8G
run new_markdup 8G "$ROOT/bin/RabbitBAM-MPI" markdup -i "$SORTED" -o "$ROOT/out/new-markdup.bam" --stream -m 8G
printf '%s PROBE COMPLETED\n' "$(date -Is)"
