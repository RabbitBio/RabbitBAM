#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:?Usage: markdup_ready_chain.sh ROOT SMALL_REFERENCE D1_SORTED D2_SORTED}
REFERENCE=${2:?Validated small markdup output required}
D1=${3:?D1 sorted fixmate input required}
D2=${4:?D2 sorted fixmate input required}
[[ "$ROOT" = /* && -s "$REFERENCE" && -s "$D1" && -s "$D2" ]]
mkdir "$ROOT/chain.lock"
STAGE=small
trap 'rc=$?; if ((rc)); then printf "FAILED %s\n" "$STAGE" > "$ROOT/chain.state"; fi; rmdir "$ROOT/chain.lock"' EXIT
printf 'WAITING_FOR_SMALL\n' > "$ROOT/chain.state"
deadline=$((SECONDS + 1800))
while true; do
    state=$(cat "$ROOT/small.state" 2>/dev/null || true)
    [[ "$state" != FAILED* ]] || exit 1
    if [[ "$state" == OUTPUTS_READY_FOR_INDEPENDENT_VERIFICATION && ! -d "$ROOT/driver.lock" ]]; then break; fi
    ((SECONDS < deadline)) || { printf 'Small regression wait timed out\n' >&2; exit 1; }
    sleep 5
done
STAGE=small_verification
printf 'RUNNING %s\n' "$STAGE" > "$ROOT/chain.state"
python3 "$ROOT/verify_markdup_ready_outputs.py" "$REFERENCE" \
    "$ROOT/out/small_ram.bam" "$ROOT/out/small_posix.bam" > "$ROOT/small.verification.log" 2>&1
cat "$ROOT/small.verification.log"
cp "$ROOT/binaries.sha256" "$ROOT/small.verified.sha256"

STAGE=D1
printf 'RUNNING %s\n' "$STAGE" > "$ROOT/chain.state"
bash "$ROOT/markdup_ready_probe.sh" "$ROOT" D1 "$D1" > "$ROOT/D1.driver.log" 2>&1
cat "$ROOT/D1.results.tsv"

STAGE=D2
printf 'RUNNING %s\n' "$STAGE" > "$ROOT/chain.state"
PROBE_DUMP=1 bash "$ROOT/markdup_ready_probe.sh" "$ROOT" D2 "$D2" > "$ROOT/D2.driver.log" 2>&1
cat "$ROOT/D2.results.tsv"

STAGE=D2_validate
printf 'RUNNING %s\n' "$STAGE" > "$ROOT/chain.state"
bash "$ROOT/markdup_ready_probe.sh" "$ROOT" D2_validate "$ROOT/out/D2_ram.bam" > "$ROOT/D2_validate.driver.log" 2>&1
cat "$ROOT/D2_validate.results.tsv"
printf 'COMPLETE\n' > "$ROOT/chain.state"
printf '%s MARKDUP READY PROBE COMPLETE\n' "$(date -Is)"
