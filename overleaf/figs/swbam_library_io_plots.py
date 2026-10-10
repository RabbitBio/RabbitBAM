import os
from PIL import Image
import matplotlib.pyplot as plt

records = {
    "D1": 177_798_571,
    "D2": 365_079_492,
    "D3": 1_210_633_370,
}
sam_cfg = [0, 2, 4, 8, 16, 24, 32, 47, 48, 64, 96, 127]
sw_cfg = [1, 2, 4, 6, 8, 12]

sam_times = {
    ("D1", "read"):  [48.06, 32.61, 16.59, 11.36, 12.21, 12.44, 11.86, 11.60, 12.19, 12.00, 11.88, 11.79],
    ("D2", "read"):  [153.97, 98.56, 52.55, 30.27, 29.76, 30.81, 27.72, 29.35, 28.86, 29.31, 28.40, 29.79],
    ("D3", "read"):  [434.73, 310.65, 162.47, 94.69, 95.16, 94.26, 96.15, 95.68, 97.36, 96.56, 99.99, 96.46],
    ("D1", "write"): [228.29, 121.95, 66.10, 36.31, 38.33, 37.07, 31.54, 31.49, 31.29, 31.35, 32.60, 32.10],
    ("D2", "write"): [691.11, 360.32, 203.85, 110.03, 91.54, 90.08, 80.42, 81.08, 79.71, 78.93, 78.46, 75.79],
    ("D3", "write"): [1926.24, 1027.26, 569.96, 307.68, 305.39, 275.33, 251.37, 252.33, 243.53, 244.50, 243.13, 246.19],
}
sw_times = {
    ("D1", "read"):  [14.762773, 7.393394, 3.709085, 2.469507, 1.860612, 1.242483],
    ("D2", "read"):  [41.171311, 20.701225, 10.428894, 6.967701, 5.240005, 3.505512],
    ("D3", "read"):  [125.631, 63.160, 31.816, 21.238, 15.972, 10.621],
    ("D1", "write"): [51.271696, 25.721544, 12.875070, 8.578293, 6.448508, 4.314136],
    ("D2", "write"): [144.242197, 72.248063, 36.279933, 24.249616, 18.181066, 12.142043],
    ("D3", "write"): [435.470, 218.574, 109.784, 73.298, 55.050, 36.651],
}

def throughput_m_records(n_records, times):
    return [n_records / t / 1e6 for t in times]

def build_panel(dataset, mode, out_path):
    n = records[dataset]
    y_sam = throughput_m_records(n, sam_times[(dataset, mode)])
    y_sw = throughput_m_records(n, sw_times[(dataset, mode)])

    x_sam = list(range(len(sam_cfg)))
    gap = 2
    x_sw = list(range(len(sam_cfg) + gap, len(sam_cfg) + gap + len(sw_cfg)))

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    ax.plot(x_sam, y_sam, marker="o", linewidth=2.5, markersize=7, label="SAMtools")
    ax.plot(x_sw, y_sw, marker="s", linewidth=2.5, markersize=7, label="SWBAM")

    sep_x = len(sam_cfg) - 0.5 + gap / 2
    ax.axvline(sep_x, linestyle="--", linewidth=1.5)

    all_x = x_sam + x_sw
    all_labels = [f"S{v}" for v in sam_cfg] + [f"W{v}" for v in sw_cfg]
    ax.set_xticks(all_x)
    ax.set_xticklabels(all_labels, rotation=45, ha="right")

    ax.set_ylabel("Throughput (M records/s)", fontsize=16)
    ax.set_xlabel("Parallel configuration (S: SAMtools -@, W: SWBAM active CGs)", fontsize=15)
    ax.set_title(f"{dataset} - {'Read-only' if mode == 'read' else 'Read-and-write'}", fontsize=18)
    ax.tick_params(axis="both", labelsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=14, loc="upper left")

    ymax = max(max(y_sam), max(y_sw)) * 1.18
    ax.set_ylim(0, ymax)

    ax.text((x_sam[0] + x_sam[-1]) / 2, ymax * 0.96, "SAMtools", ha="center", va="top", fontsize=14)
    ax.text((x_sw[0] + x_sw[-1]) / 2, ymax * 0.96, "SWBAM", ha="center", va="top", fontsize=14)

    best_sam_idx = max(range(len(y_sam)), key=lambda i: y_sam[i])
    best_sw_idx = max(range(len(y_sw)), key=lambda i: y_sw[i])
    ax.annotate(f"Peak {y_sam[best_sam_idx]:.2f}",
                (x_sam[best_sam_idx], y_sam[best_sam_idx]),
                xytext=(8, 8), textcoords="offset points", fontsize=12)
    ax.annotate(f"Peak {y_sw[best_sw_idx]:.2f}",
                (x_sw[best_sw_idx], y_sw[best_sw_idx]),
                xytext=(8, 8), textcoords="offset points", fontsize=12)

    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    out_dir = "/mnt/data/swbam_library_io_plots"
    os.makedirs(out_dir, exist_ok=True)
    for dataset in ["D1", "D2", "D3"]:
        for mode in ["read", "write"]:
            build_panel(dataset, mode, os.path.join(out_dir, f"{dataset}_{mode}.png"))
