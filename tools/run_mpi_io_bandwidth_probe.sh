#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: run_mpi_io_bandwidth_probe.sh ABS_ROOT SIX_NODES}
NODES=${2:?Six idle nodes required}
[[ "$ROOT" == /* && "$NODES" =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]]
BIN="$ROOT/bin/mpi-io-bandwidth-probe"
LOG="$ROOT/logs/probe.log"
[[ -x "$BIN" && ! -e "$LOG" ]]
mkdir -p "$ROOT/logs" "$ROOT/data"
cd "$ROOT"
sha256sum "$BIN" | tee "$ROOT/binary.sha256"
command=(bsub -N 6 -np 1 -node "$NODES" -health 1 -mpecg 6 -cgsp 64
    -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
    -a mpi -b -K -q q_share -timelimit 00:20:00 -J io_n6_probe -o "$LOG"
    "$BIN" "$ROOT/data" 512 512 3 16 6)
printf '%q ' "${command[@]}" > "$ROOT/commands.sh"
printf '\n' >> "$ROOT/commands.sh"
if ! submission=$("${command[@]}" 2>&1); then
    printf '%s\n' "$submission" | tee "$ROOT/submit.log"
    exit 1
fi
printf '%s\n' "$submission" | tee "$ROOT/submit.log"
job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
[[ "$job" =~ ^[0-9]+$ ]]
printf '%s\n' "$job" > "$ROOT/job.id"
start=$SECONDS
running_since=0
last_report=-60
while true; do
    state=$(bjobs "$job" 2>/dev/null | awk -v id="$job" '$1 == id {print $2}')
    if ((SECONDS - last_report >= 30)); then
        printf '%s job=%s state=%s\n' "$(date -Is)" "$job" "${state:-retired}"
        last_report=$SECONDS
    fi
    case "$state" in
        RUN|EXITING)
            if ((running_since == 0)); then
                running_since=$SECONDS
                cnload -j "$job" -l > "$ROOT/nodes.txt" || true
            fi
            if ((SECONDS - running_since > 1200)); then
                printf 'IO_PROBE_TIMEOUT job=%s; cancelling only this probe\n' "$job"
                bkill "$job" || true
                exit 1
            fi
            ;;
        PEND)
            if ((SECONDS - start > 300)); then
                printf 'IO_PROBE_QUEUE_TIMEOUT job=%s; cancelling only this probe\n' "$job"
                bkill "$job" || true
                exit 1
            fi
            ;;
        DONE|EXIT|KILLING|'') break ;;
        *) printf 'Unexpected probe state: %s\n' "$state"; exit 1 ;;
    esac
    sleep 5
done
for ((i=0; i<20; ++i)); do
    if [[ -s "$LOG" ]] && grep -q 'IO_PROBE_COMPLETE ranks=6 status=0' "$LOG"; then
        break
    fi
    sleep 2
done
grep -q 'IO_PROBE_COMPLETE ranks=6 status=0' "$LOG"
! grep -q 'IO_ERROR' "$LOG"
grep '^IO_SUMMARY ' "$LOG" | tee "$ROOT/summary.txt"
find "$ROOT/data" -maxdepth 1 -type f -printf '%f %s\n' > "$ROOT/remaining_files.txt"
[[ ! -s "$ROOT/remaining_files.txt" ]]
printf '%s IO_PROBE_VERIFIED job=%s\n' "$(date -Is)" "$job"
