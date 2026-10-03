#!/usr/bin/env bash
set -euo pipefail

RUN_DIR=${1:?Usage: bench_apps_d2_continue.sh absolute-run-directory pilot-unit}
PILOT_UNIT=${2:?Usage: bench_apps_d2_continue.sh absolute-run-directory pilot-unit}

while systemctl --user is-active --quiet "$PILOT_UNIT.service"; do
    sleep 30
done

if [[ $(<"$RUN_DIR/state") != PILOT_COMPLETE* ]]; then
    printf 'D2 capacity pilot did not complete; matrix was not started.\n' >&2
    exit 1
fi

printf 'D2 capacity pilot passed at %s; starting the matrix.\n' "$(date -Is)"
SWBAM_APPS_RESUME=1 /bin/bash "$RUN_DIR/driver.sh" "$RUN_DIR"
