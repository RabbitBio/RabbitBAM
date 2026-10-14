# Validate unique per-rank successes before accepting a performance sample.
{ sub(/\r$/, "") }
/\[swbam-e2e rank=/ {
    rank = $2; sub(/^rank=/, "", rank); sub(/\]$/, "", rank)
    if (rank !~ /^[0-9]+$/ || rank+0 >= ranks+0 || eids[rank]++) bad = 1
    ok = 0
    for (i = 1; i <= NF; i++) {
        split($i, a, "=")
        if (a[1] == "seconds") e_by_rank[rank] = a[2]+0
        if (a[1] == "seconds" && a[2] > e) e = a[2]
        if (a[1] == "status" && a[2] == 0) ok = 1
    }
    if (!ok) bad = 1
    ne++
}
/\[swbam-processing rank=/ {
    rank = $2; sub(/^rank=/, "", rank); sub(/\]$/, "", rank)
    if (rank !~ /^[0-9]+$/ || rank+0 >= ranks+0 || pids[rank]++) bad = 1
    ok = 0; seconds = ""
    for (i = 1; i <= NF; i++) {
        split($i, a, "=")
        if (a[1] == "seconds") seconds = a[2]
        if (a[1] == "status" && a[2] == 0) ok = 1
    }
    if (!ok || seconds !~ /^[0-9]+([.][0-9]+)?$/ || seconds+0 <= 0) bad = 1
    processing_by_rank[rank] = seconds+0
    if (seconds+0 > processing) processing = seconds+0
    np++
}
/\[ram-storage rank=/ {
    rank = $2; sub(/^rank=/, "", rank); sub(/\]$/, "", rank)
    if (rank !~ /^[0-9]+$/ || rank+0 >= ranks+0 || sids[rank]++) bad = 1
    ns++
    for (i = 1; i <= NF; i++) {
        split($i, a, "=")
        if (a[1] == "peak" && a[2] > peak) peak = a[2]
        if (a[1] == "output_bytes" && a[2] > bytes) bytes = a[2]
        # Markdup's bounded deferred payload may spill into actual RAM files.
        # Sort/collate still must remain in memory mode for this matrix.
        if ((a[1] == "temp_read_bytes" || a[1] == "temp_write_bytes") &&
            a[2] != 0 && check != "markdup") bad = 1
    }
}
/\[ram-prepare rank=/ {
    for (i = 1; i <= NF; i++) {
        split($i, a, "=")
        if (a[1] == "seconds" && a[2] > prepare) prepare = a[2]
    }
}
/^sdk_(count|bam2bam) |^Optimized[^ ]* finished[.]/ {
    for (i = 1; i <= NF; i++) {
        split($i, a, "=")
        if ((a[1] == "records" || a[1] == "total_records" || a[1] == "local_records") && a[2] == records) have_records = 1
        if (a[1] == "ranks" && a[2] == ranks) have_ranks = 1
        if (a[1] == "kept_records" && a[2] == kept) have_kept = 1
        if (a[1] == "core") core = a[2]+0
    }
}
/^Complete the 4[.]3 / { core = $NF+0 }
/selected mode=memory/ { memory = 1 }
/^[[:space:]]*duplicates / {
    for (i = 1; i <= NF; i++) {
        split($i, a, "=")
        if (a[1] == "marked" && a[2] == marked) have_marked = 1
    }
}
/^SN\tis sorted:\t1/ { sorted = 1 }
/^SN\treads duplicated:\t/ { split($0, a, "\t"); if (a[3] == marked) have_dup = 1 }
END {
    for (rank in processing_by_rank)
        if (processing_by_rank[rank] > e_by_rank[rank] + 0.000002) bad = 1
    if (bad || ne != ranks || np != ranks || ns != ranks || !have_records || !have_ranks || e <= 0 || core <= 0) exit 1
    if (check == "filter" && (!have_kept || bytes <= 0)) exit 1
    if (check == "memory" && (!memory || bytes <= 0)) exit 1
    if (check == "markdup" && (!have_marked || bytes <= 0)) exit 1
    if (check == "write" && bytes <= 0) exit 1
    if ((check == "sorted" || check == "duplicates") && !sorted) exit 1
    if (check == "duplicates" && !have_dup) exit 1
    printf "%.6f\t%.6f\t%.6f\t%.0f\t%.0f\t%.6f", e, core, prepare, peak, bytes, processing
}
