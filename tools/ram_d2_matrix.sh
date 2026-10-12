#!/usr/bin/env bash
set -euo pipefail
export BENCH_DATASET=D2
exec bash "$(dirname "${BASH_SOURCE[0]}")/ram_matrix.sh" "$@"
