#!/usr/bin/env bash
set -euo pipefail

BENCH=${1:?Usage: launch_ram_markdup_retests.sh ABS_BENCH D1_RETEST_ROOT}
D1_ROOT=${2:?New D1 retest directory required}
[[ "$BENCH" = /* && "$D1_ROOT" = "$BENCH"/* ]]
D1_PRIOR="$BENCH/ram_D1_processing_matrix_20261013_043000"
D2_ROOT="$BENCH/ram_D2_processing_matrix_20261013_114500"
VERIFIED="$BENCH/markdup_ready_20261013_181000"
EXPECTED=3dc30d0a9ac3ba2a384fec617b66f32d7801a53d389288db0f51f110c7726a51
[[ ! -e "$D1_ROOT/state" && $(<"$VERIFIED/chain.state") == COMPLETE ]]
[[ $(sha256sum "$VERIFIED/bin/RabbitBAM-MPI" | awk '{print $1}') == "$EXPECTED" ]]
cmp "$D1_PRIOR/binaries.sha256" "$D2_ROOT/binaries.sha256"
for binary in RabbitBAM-MPI swbam-sdk-optimized-count swbam-sdk-optimized-bam2bam; do
    cp -p "$D1_PRIOR/bin/$binary" "$D1_ROOT/bin/$binary"
done
cp -p "$VERIFIED/bin/RabbitBAM-MPI" "$D1_ROOT/bin/RabbitBAM-MPI-markdup"
cp -p "$VERIFIED/bin/RabbitBAM-MPI" "$D2_ROOT/bin/RabbitBAM-MPI-markdup"
cp -p "$D1_PRIOR/launch.sh" "$D1_ROOT/launch.sh"
cp -p "$BENCH/ram_D1_matrix_20261012_010800/retry_watchdog.sh" "$D2_ROOT/bin/retry_watchdog.sh"
sha256sum "$D2_ROOT/bin/RabbitBAM-MPI-markdup" > "$D2_ROOT/markdup.next.sha256"

bash -n "$D1_ROOT/bin/ram_matrix.sh" "$D2_ROOT/bin/ram_matrix.next.sh" \
    "$D2_ROOT/recover_ram_d2_markdup.sh"
(
    cd "$D1_ROOT"
    sha256sum launch.sh bin/ram_d1_matrix.sh bin/ram_matrix.sh bin/ram_metrics.awk > scripts.sha256
    date -Is > started.txt
    printf 'scope=markdup-three-samples/three-full-workflows commands=21 prior=%s\n' \
        "$D1_PRIOR" > retest_scope.txt
    nohup env BENCH_STAGE=markdup-workflow BENCH_MARKDUP_APP="$D1_ROOT/bin/RabbitBAM-MPI-markdup" \
        bash "$D1_ROOT/launch.sh" D1 "$D1_ROOT" > driver.log 2>&1 < /dev/null &
    printf '%s\n' "$!" > driver.pid
    printf 'D1 driver_pid=%s root=%s\n' "$(<driver.pid)" "$D1_ROOT"
)
bash "$D2_ROOT/recover_ram_d2_markdup.sh" "$D2_ROOT" 8543367 13138 43958-43963
