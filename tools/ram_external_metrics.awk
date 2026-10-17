# Append external-memory diagnostics after the common rank/timing checks.
/selected mode=external/ { external = 1 }
/selected mode=memory/ { wrong_mode = 1 }
/\[ram-storage rank=/ {
    for (i = 1; i <= NF; ++i) {
        split($i, a, "=")
        if (a[1] == "temp_read_bytes") temp_read += a[2]+0
        if (a[1] == "temp_write_bytes") temp_write += a[2]+0
    }
}
/external_sum |external_actual_sum |external_memory_sum / {
    for (i = 1; i <= NF; ++i) {
        split($i, a, "=")
        if (a[1] == "runs") runs = a[2]+0
        if (a[1] == "tracked_peak_max") work_peak = a[2]+0
    }
}
END {
    if (!external || wrong_mode || temp_read <= 0 || temp_write <= 0 ||
        runs <= ranks || work_peak <= 0) exit 1
    printf "\t%.0f\t%.0f\t%.0f\t%.0f", temp_read, temp_write, work_peak, runs
}
