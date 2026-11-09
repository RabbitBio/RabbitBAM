#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:?Usage: submit_markdup_pos_halo.sh ABS_ROOT regression|small|D1|D2_final|D3_final|D3_partial}
PHASE=${2:?Phase required}
APP=${APP:-$ROOT/bin/RabbitBAM-MPI}
DATA=/home/export/online1/mdt00/shisuan/swls-CFD/guoshi/wzs/wzs_data
[[ "$ROOT" = /* && -d "$ROOT/logs" && -d "$ROOT/out" ]]
cd "$ROOT"
[[ ! -e "logs/$PHASE.log" ]]
command=(bsub -a mpi -b -q q_share -timelimit 01:00:00
    -J "md_pos_halo_$PHASE" -o "$ROOT/logs/$PHASE.log")
if [[ -n "${BENCH_NODES:-}" ]]; then
    command+=(-node "$BENCH_NODES")
fi
case "$PHASE" in
    regression|regression_final|regression_retry1|regression_paired|regression_partial)
        command+=(-N 1 -np 2 -cgsp 64 -share_size 14000 -cache_size 32
            "${HALO_PROBE:-$ROOT/bin/swbam-markdup-halo-regression}" "$ROOT/out")
        ;;
    small|small_retry1|small_final|small_partial)
        command+=(-N 1 -np 6 -cgsp 64 -share_size 14000 -cache_size 32
            -x SWBAM_NODE_MEMORY=88000M
            "$APP" markdup
            -i "$DATA/swbam_bench/ram_D3_mouse_probe_20261016_023500/out/small-fm_fc7f37cbc87a/sort.bam"
            -o "$ROOT/out/$PHASE.bam" --stream -l 300 -m 8G --compress-level 1 --storage posix)
        ;;
    D1|D1_final|D2_final|D3_final|D3_partial)
        case "$PHASE" in
            D1|D1_final)
                input="$DATA/swbam_bench/D1_apps_current/apps_D1_n6_20261001_031800/out/fixmate.sorted.bam"
                capacity=8G
                ;;
            D2_final)
                input="$DATA/swbam_bench/D2_apps_current/apps_D2_D3_workflows_n6_20261010_003937_D2/out/fixmate.sorted.bam"
                capacity=24G
                ;;
            D3_final|D3_partial)
                input="$DATA/swbam_bench/ram_D3_mouse_probe_20261016_023500/out/full-fm_fc7f37cbc87a/sort.bam"
                capacity=24G
                ;;
        esac
        [[ -s "$input" && ! -e "$ROOT/sessions/$PHASE" ]]
        command+=(-N 6 -np 1 -health 1 -mpecg 6 -cgsp 64
            -share_size 100 -cross_size 88000 -cache_size 32 -xmalloc
            -x SWBAM_NODE_MEMORY=88000M
            "$APP" markdup
            -i "$input"
            -o "$ROOT/out/$PHASE.bam" --stream -l 300 -m 8G --compress-level 1
            --storage ram --ram-capacity "$capacity" --ram-ranks 6
            --ram-session "$ROOT/sessions/$PHASE" --ram-ready-timeout 1800 --no-dump)
        ;;
    *) printf 'Unsupported phase: %s\n' "$PHASE" >&2; exit 1 ;;
esac
printf '%q ' "${command[@]}" > "$PHASE.command.sh"
printf '\n' >> "$PHASE.command.sh"
"${command[@]}"
