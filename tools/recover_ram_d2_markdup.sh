#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:?Usage: recover_ram_d2_markdup.sh ABS_ROOT STUCK_JOB DRIVER_PID RETRY_NODES_OR_AUTO [RETRY_INDEX]}
STUCK_JOB=${2:?Stuck validation job required}
DRIVER_PID=${3:?Current driver PID required}
RETRY_NODES=${4:?Replacement nodes required}
RETRY_INDEX=${5:-1}
LABEL=validation_sdk_n6_r1
MARKDUP_APP="$ROOT/bin/RabbitBAM-MPI-markdup"
[[ "$ROOT" = /* && "$STUCK_JOB" =~ ^[1-9][0-9]*$ &&
   "$DRIVER_PID" =~ ^[1-9][0-9]*$ && "$RETRY_INDEX" =~ ^[1-9][0-9]*$ &&
   ("$RETRY_NODES" == auto || "$RETRY_NODES" =~ ^[0-9,-]+$) ]]
cd "$ROOT"
[[ $(<state) == "RUNNING $LABEL" && $(<driver.pid) == "$DRIVER_PID" ]]
[[ $(awk -F '\t' -v label="$LABEL" '$1==label{print $2}' jobs.tsv) == "$STUCK_JOB" ]]
[[ $(ps -p "$DRIVER_PID" -o args=) == *"$ROOT/bin/ram_matrix.sh"* ]]
old_status=$(bjobs "$STUCK_JOB" 2>/dev/null | awk -v id="$STUCK_JOB" '$1==id{print $2}')
[[ "$old_status" == RUN || "$old_status" == PEND ]]
[[ -x "$MARKDUP_APP" && -s out/sdk_n1_r1.bam && -s bin/ram_matrix.next.sh ]]
sha256sum -c binaries.sha256
sha256sum -c markdup.next.sha256

ARCHIVE="$ROOT/retries/$LABEL.stopped_job$STUCK_JOB"
[[ ! -e "$ARCHIVE" && ! -e "$ROOT/logs/$LABEL.retry$RETRY_INDEX.log" ]]
mkdir -p "$ARCHIVE"
driver_stopped=0
trap 'rc=$?; if ((rc && driver_stopped)); then kill -CONT "$DRIVER_PID" 2>/dev/null || true; fi' EXIT
kill -STOP "$DRIVER_PID"
driver_stopped=1
cp -a state results.tsv workflows.tsv jobs.tsv commands.sh binaries.sha256 \
    driver.log driver.pid scripts.sha256 bin/ram_matrix.sh "$ARCHIVE/"
if [[ -d "sessions/$LABEL" ]]; then cp -a "sessions/$LABEL" "$ARCHIVE/session"; fi
bjobs -l "$STUCK_JOB" > "$ARCHIVE/job.txt"
qload q_share -w -n 6 > "$ARCHIVE/nodes.txt"
if [[ -e "logs/$LABEL.log" || -L "logs/$LABEL.log" ]]; then
    mv "logs/$LABEL.log" "$ARCHIVE/original.log"
fi
printf '%s stopping job=%s driver=%s retry_nodes=%s\n' \
    "$(date -Is)" "$STUCK_JOB" "$DRIVER_PID" "$RETRY_NODES" > "$ARCHIVE/recovery.log"
bkill "$STUCK_JOB"
for ((attempt=0; attempt<60; ++attempt)); do
    status=$(bjobs "$STUCK_JOB" 2>/dev/null | awk -v id="$STUCK_JOB" '$1==id{print $2}' || true)
    case "$status" in EXIT|DONE|"") break ;; esac
    sleep 1
done
[[ "$status" == EXIT || "$status" == DONE || -z "$status" ]]
kill -TERM "$DRIVER_PID"
kill -CONT "$DRIVER_PID" 2>/dev/null || true
for ((attempt=0; attempt<20; ++attempt)); do
    process_state=$(ps -p "$DRIVER_PID" -o stat= || true)
    [[ -z "$process_state" || "$process_state" == Z* ]] && break
    sleep 1
done
[[ -z "$process_state" || "$process_state" == Z* ]]
driver_stopped=0
if [[ -d driver.lock ]]; then rmdir driver.lock; fi

cp bin/ram_matrix.next.sh bin/ram_matrix.sh
cp markdup.next.sha256 markdup.sha256
sha256sum -c binaries.sha256
sha256sum -c markdup.sha256
printf 'markdup_binary=%s migration=%s affected=app-markdup/workflow-markdup\n' \
    "$MARKDUP_APP" "$(date -Is)" >> configuration.txt
sha256sum launch.sh bin/ram_d2_matrix.sh bin/ram_matrix.sh bin/ram_metrics.awk > scripts.sha256

RETRY_LABEL="$LABEL.retry$RETRY_INDEX"
command=(bsub -N 6 -np 1 -health 1 -mpecg 6 -cgsp 64
    -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
    -a mpi -b -K -q q_share -timelimit 03:00:00)
if [[ "$RETRY_NODES" != auto ]]; then command+=(-node "$RETRY_NODES"); fi
command+=(-J "d2ram_$RETRY_LABEL" -o "$ROOT/logs/$RETRY_LABEL.log"
    "$ROOT/bin/RabbitBAM-MPI" stats --basic -i "$ROOT/out/sdk_n1_r1.bam"
    --storage ram --ram-capacity 8G --ram-ranks 6
    --ram-session "$ROOT/sessions/$RETRY_LABEL" --ram-ready-timeout 5400)
export SWBAM_NODE_MEMORY=88000M
printf '%q ' "${command[@]}" >> commands.sh
printf '\n' >> commands.sh
submission=$("${command[@]}")
printf '%s\n' "$submission"
job=$(sed -n 's/.*Job <\([0-9]*\)>.*/\1/p' <<< "$submission")
[[ "$job" =~ ^[1-9][0-9]*$ ]]
awk -F '\t' -v label="$LABEL" '$1!=label' jobs.tsv > jobs.next.tsv
printf '%s\t%s\n' "$LABEL" "$job" >> jobs.next.tsv
mv jobs.next.tsv jobs.tsv
ln -s "$RETRY_LABEL.log" "logs/$LABEL.log"
printf '%s\t%s\t%s\t%s\n' "$LABEL" "$STUCK_JOB" "$job" "$ARCHIVE" >> retries.tsv
printf 'RUNNING %s retry_job=%s\n' "$LABEL" "$job" > state
wait_limit=3600
if [[ "$RETRY_NODES" == auto ]]; then wait_limit=21600; fi
nohup env RETRY_WAIT_LIMIT="$wait_limit" bash "$ROOT/bin/retry_watchdog.sh" "$job" \
    > "$ROOT/retry_watchdog.$job.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$ROOT/retry_watchdog.$job.pid"
nohup env BENCH_RESUME=1 BENCH_MARKDUP_APP="$MARKDUP_APP" \
    bash "$ROOT/launch.sh" "$ROOT" > driver.log 2>&1 < /dev/null &
pid=$!
printf '%s\n' "$pid" > driver.pid
printf '%s retry_job=%s driver_pid=%s new_markdup=%s\n' \
    "$(date -Is)" "$job" "$pid" "$MARKDUP_APP" >> "$ARCHIVE/recovery.log"
printf 'retry_job=%s driver_pid=%s archive=%s\n' "$job" "$pid" "$ARCHIVE"
