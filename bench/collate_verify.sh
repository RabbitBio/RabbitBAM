#!/bin/bash

set -euo pipefail
export LC_ALL=C

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 NAME.sam [TEMP_DIR]" >&2
    exit 2
fi

sam_path=$1
temp_root=${2:-${TMPDIR:-.}}
sort_memory=${SORT_MEMORY:-4G}

if [ ! -r "$sam_path" ]; then
    echo "ERROR: cannot read SAM file: $sam_path" >&2
    exit 1
fi

mkdir -p "$temp_root"
work_dir=$(mktemp -d "$temp_root/collate-verify.XXXXXX")
cleanup() {
    rm -rf "$work_dir" 2>/dev/null || true
}
trap cleanup EXIT

stats_file="$work_dir/stats.txt"
header_file="$work_dir/header.no_pg.sam"
unique_file="$work_dir/unique_qnames.txt"
ordered_md5_file="$work_dir/ordered_body.md5"
multiset_md5_file="$work_dir/record_multiset.md5"
ordered_fifo="$work_dir/ordered_body.fifo"
multiset_fifo="$work_dir/record_multiset.fifo"

mkfifo "$ordered_fifo" "$multiset_fifo"
: > "$header_file"

echo "collate_verify: scanning once and sorting once: $sam_path" >&2
echo "collate_verify: sort temp=$work_dir sort_memory=$sort_memory" >&2

md5sum < "$ordered_fifo" > "$ordered_md5_file" &
ordered_md5_pid=$!
md5sum < "$multiset_fifo" > "$multiset_md5_file" &
multiset_md5_pid=$!

awk -F '\t' -v stats_file="$stats_file" -v header_file="$header_file" '
    /^@/ {
        if ($0 !~ /^@PG/) print $0 > header_file
        next
    }
    {
        records++
        order = int($2 / 64) % 4
        if ($1 != qname) {
            groups++
            qname = $1
            previous_order = -1
        }
        if (order < previous_order) read_order_errors++
        previous_order = order
        print
    }
    END {
        printf("records %d\n", records + 0) > stats_file
        printf("groups %d\n", groups + 0) >> stats_file
        printf("read_order_errors %d\n", read_order_errors + 0) >> stats_file
    }
' "$sam_path" |
    tee "$ordered_fifo" |
    sort -S "$sort_memory" -T "$work_dir" |
    tee "$multiset_fifo" |
    awk -F '\t' '
        $1 != previous {
            unique_qnames++
            previous = $1
        }
        END { print unique_qnames + 0 }
    ' > "$unique_file"

wait "$ordered_md5_pid"
wait "$multiset_md5_pid"

read_stat() {
    awk -v key="$1" '$1 == key { print $2; found = 1 } END { if (!found) exit 1 }' "$stats_file"
}

records=$(read_stat records)
groups=$(read_stat groups)
read_order_errors=$(read_stat read_order_errors)
unique_qnames=$(cat "$unique_file")
continuity_errors=$((groups - unique_qnames))
header_without_pg_md5=$(md5sum "$header_file" | awk '{print $1}')
ordered_body_md5=$(awk '{print $1}' "$ordered_md5_file")
record_multiset_md5=$(awk '{print $1}' "$multiset_md5_file")

printf 'records %s\n' "$records"
printf 'groups %s\n' "$groups"
printf 'unique_qnames %s\n' "$unique_qnames"
printf 'continuity_errors %s\n' "$continuity_errors"
printf 'read_order_errors %s\n' "$read_order_errors"
printf 'header_without_pg_md5 %s\n' "$header_without_pg_md5"
printf 'ordered_body_md5 %s\n' "$ordered_body_md5"
printf 'record_multiset_md5 %s\n' "$record_multiset_md5"

echo "collate_verify: done" >&2
