#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: launch_a4_d3_mouse.sh ABS_ROOT REFERENCE_ROOT}
REFERENCE=${2:?Verified D3 reference root required}
INPUT=/home/bigssd/wzs_data/D3-mouse-SRR7768905/initial/D3-mouse-SRR7768905.coordinate.bam
[[ "$ROOT" = /* && "$REFERENCE" = /* && -s "$INPUT" ]]
[[ -f "$(dirname "$INPUT")/.transfer_verified" ]]
cd "$ROOT"
sha256sum -c binaries.sha256
sha256sum -c scripts.sha256
python3 -m unittest discover -s "$ROOT/bin" -p test_bench_samtools_d3_mouse.py -v
python3 -u "$ROOT/bin/bench_samtools_d3_mouse.py" \
    --root "$ROOT.smoke" --input "$REFERENCE/small.bam" \
    --samtools "$ROOT/bin/samtools" --smoke --records 200000 --kept 176298 --duplicates 22259
[[ -f "$ROOT.smoke/ALL_DONE" ]]
python3 -u "$ROOT/bin/bench_samtools_d3_mouse.py" \
    --root "$ROOT" --input "$INPUT" --samtools "$ROOT/bin/samtools" \
    --reference-root "$REFERENCE" --resource-wait-seconds 1800
