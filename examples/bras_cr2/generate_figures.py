#!/usr/bin/env python3
"""Generate figures for the BRAS CR2 benchmark README.

Reads simulation results from results.csv (produced by the Rust benchmark)
and the room geometry from the AC3D file.

Produces:
  - results.png: Simulated vs measured room acoustic metrics (EDT, T20, C80, D50)
  - geometry.png: 3D wireframe of the CR2 seminar room with source/receiver positions

Usage:
  # First, run the benchmark to produce results.csv:
  cargo run --example bras_cr2 --release

  # Then generate figures:
  python3 examples/bras_cr2/generate_figures.py
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

SCRIPT_DIR = Path(__file__).parent

# ── Source and receiver positions (building3d z-up coordinates) ──

SOURCES = {
    "LS1": (0.931, -2.547, 0.723),
    "LS2": (0.119, 2.880, 0.723),
}
RECEIVERS = {
    "MP1": (-0.993, 1.426, 1.230),
    "MP2": (0.439, -0.147, 1.230),
    "MP3": (1.361, -0.603, 1.230),
    "MP4": (-1.110, -0.256, 1.230),
    "MP5": (-0.998, -1.409, 1.230),
}


def load_results_csv(path):
    """Load results.csv and return {metric: {"freqs": [...], "sim": [...], "meas": [...]}}."""
    data = defaultdict(lambda: {"freqs": [], "sim": [], "meas": []})
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            metric = row["metric"]
            data[metric]["freqs"].append(float(row["freq_hz"]))
            data[metric]["sim"].append(float(row["simulated"]))
            data[metric]["meas"].append(float(row["measured"]))
    return dict(data)


def freq_label(f):
    """Convert frequency to compact label: 1000 -> '1k', 125 -> '125'."""
    if f >= 1000:
        return f"{f / 1000:.0f}k"
    return f"{f:.0f}"


def generate_results_chart(data):
    """Generate a 2x2 comparison chart of simulated vs measured metrics."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle("BRAS CR2: Simulated vs Measured Room Acoustic Parameters",
                 fontsize=13, fontweight="bold", y=0.98)

    plots = [
        ("EDT", "EDT (s)", "Early Decay Time"),
        ("T20", "T20 (s)", "Reverberation Time (T20)"),
        ("C80", "C80 (dB)", "Clarity"),
        ("D50", "D50 (%)", "Definition"),
    ]

    width = 0.32

    for ax, (metric, ylabel, title) in zip(axes.flat, plots):
        d = data[metric]
        labels = [freq_label(f) for f in d["freqs"]]
        x = np.arange(len(labels))

        ax.bar(x - width / 2, d["sim"], width, label="Simulated",
               color="#4C72B0", edgecolor="white", linewidth=0.5)
        ax.bar(x + width / 2, d["meas"], width, label="Measured",
               color="#DD8452", edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Frequency (Hz)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.tick_params(axis="y", labelsize=8)
        ax.legend(fontsize=8, loc="best")
        ax.grid(axis="y", alpha=0.3, linewidth=0.5)
        ax.set_axisbelow(True)

        if "C80" in ylabel:
            ax.axhline(y=0, color="gray", linewidth=0.5, linestyle="--")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = SCRIPT_DIR / "results.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {out}")


def parse_ac3d(path):
    """Parse an AC3D file and return lists of (vertices, faces) per object.
    Each face is (vertex_indices, material_index). Coordinates are RAVEN y-up."""
    text = Path(path).read_text()
    objects = []

    parts = re.split(r"OBJECT poly", text)
    for part in parts[1:]:
        m = re.search(r"numvert\s+(\d+)", part)
        if not m:
            continue
        nv = int(m.group(1))
        vert_block = part[m.end():]
        verts = []
        for line in vert_block.strip().split("\n")[:nv]:
            coords = line.strip().split()
            if len(coords) >= 3:
                verts.append([float(c) for c in coords[:3]])
        verts = np.array(verts)

        faces = []
        for sm in re.finditer(r"SURF\s+\S+\s*\nmat\s+(\d+)\s*\nrefs\s+(\d+)", part):
            mat_idx = int(sm.group(1))
            nrefs = int(sm.group(2))
            ref_start = sm.end()
            ref_lines = part[ref_start:].strip().split("\n")[:nrefs]
            face = []
            for rl in ref_lines:
                idx = int(rl.strip().split()[0])
                face.append(idx)
            if len(face) >= 3:
                faces.append((face, mat_idx))

        objects.append((verts, faces))

    return objects


def generate_geometry_image():
    """Generate a 3D view of the CR2 room geometry with source/receiver positions."""
    ac3d_path = Path("validation/bras/informed_sim/RavenModels/scene9.ac")
    if not ac3d_path.exists():
        print(f"  WARNING: {ac3d_path} not found, skipping geometry figure")
        return

    objects = parse_ac3d(ac3d_path)

    # Material colors: concrete, windows, ceiling, plaster, floor
    mat_colors = [
        (0.65, 0.64, 0.62, 0.25),
        (0.30, 0.65, 0.60, 0.30),
        (0.75, 0.75, 0.75, 0.15),
        (0.70, 0.68, 0.62, 0.25),
        (0.60, 0.55, 0.42, 0.35),
    ]

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection="3d")

    def convert(verts):
        """RAVEN (x, y_up, z_raven) -> building3d (x, -z_raven, y_up)"""
        out = np.zeros_like(verts)
        out[:, 0] = verts[:, 0]
        out[:, 1] = -verts[:, 2]
        out[:, 2] = verts[:, 1]
        return out

    for verts_raw, faces in objects:
        verts = convert(verts_raw)
        for face_indices, mat_idx in faces:
            if len(face_indices) < 3:
                continue
            polygon = [verts[i] for i in face_indices]
            color = mat_colors[mat_idx] if mat_idx < len(mat_colors) else (0.5, 0.5, 0.5, 0.2)
            poly = Poly3DCollection([polygon], alpha=color[3],
                                    facecolor=color[:3],
                                    edgecolor=(0.3, 0.3, 0.3, 0.4),
                                    linewidth=0.3)
            ax.add_collection3d(poly)

    for name, (sx, sy, sz) in SOURCES.items():
        ax.scatter(sx, sy, sz, c="red", s=80, marker="^", zorder=10,
                   edgecolors="darkred", linewidths=0.8, depthshade=False)
        ax.text(sx + 0.15, sy, sz + 0.15, name, fontsize=8, fontweight="bold",
                color="red")

    for name, (rx, ry, rz) in RECEIVERS.items():
        ax.scatter(rx, ry, rz, c="blue", s=50, marker="o", zorder=10,
                   edgecolors="darkblue", linewidths=0.8, depthshade=False)
        ax.text(rx + 0.15, ry, rz + 0.1, name, fontsize=7, color="blue")

    ax.set_xlabel("x (m)", fontsize=9)
    ax.set_ylabel("y (m)", fontsize=9)
    ax.set_zlabel("z (m)", fontsize=9)
    ax.set_title("BRAS CR2: Seminar Room Geometry", fontsize=12, fontweight="bold")
    ax.tick_params(labelsize=7)

    all_verts = np.vstack([convert(v) for v, _ in objects])
    mid = (all_verts.max(axis=0) + all_verts.min(axis=0)) / 2
    span = (all_verts.max(axis=0) - all_verts.min(axis=0)).max() / 2 * 1.1
    ax.set_xlim(mid[0] - span, mid[0] + span)
    ax.set_ylim(mid[1] - span, mid[1] + span)
    ax.set_zlim(mid[2] - span, mid[2] + span)

    ax.view_init(elev=25, azim=-55)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="^", color="w", markerfacecolor="red",
               markersize=10, label="Sources (LS1, LS2)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="blue",
               markersize=8, label="Receivers (MP1--MP5)"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc="upper left")

    out = SCRIPT_DIR / "geometry.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {out}")


if __name__ == "__main__":
    results_csv = SCRIPT_DIR / "results.csv"
    if not results_csv.exists():
        print(f"ERROR: {results_csv} not found.")
        print("Run the benchmark first:  cargo run --example bras_cr2 --release")
        sys.exit(1)

    print("Generating BRAS CR2 figures...")
    data = load_results_csv(results_csv)
    generate_results_chart(data)
    generate_geometry_image()
    print("Done.")
