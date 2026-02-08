#!/usr/bin/env python3
"""Generate figures for the BRAS CR2 benchmark README.

Produces:
  - results.png: Simulated vs measured room acoustic metrics (EDT, RT60, C80, D50)
  - geometry.png: 3D wireframe of the CR2 seminar room with source/receiver positions

Usage:
  python3 examples/bras_cr2/generate_figures.py
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUTPUT_DIR = Path(__file__).parent

# ── Simulation results (5,000 rays, averaged over 10 source-receiver pairs) ──

FREQS = [125, 250, 500, 1000, 2000, 4000]
FREQ_LABELS = ["125", "250", "500", "1k", "2k", "4k"]

MEASURED_EDT = [1.501, 1.306, 2.018, 1.919, 1.758, 1.649]
SIM_EDT = [1.55, 1.48, 2.13, 2.03, 1.81, 1.76]

MEASURED_RT60 = [1.453, 1.345, 2.018, 1.934, 1.715, 1.604]
SIM_RT60 = [1.58, 1.55, 2.16, 2.05, 1.81, 1.76]

MEASURED_C80 = [0.551, 0.902, -1.862, -0.606, -0.483, 0.164]
SIM_C80 = [-0.08, 0.27, -1.93, -1.67, -1.02, -0.84]

MEASURED_D50 = [29.5, 39.4, 27.8, 32.1, 35.6, 36.6]
SIM_D50 = [34.0, 35.7, 26.0, 27.1, 29.8, 30.6]

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


def generate_results_chart():
    """Generate a 2x2 comparison chart of simulated vs measured metrics."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle("BRAS CR2: Simulated vs Measured Room Acoustic Parameters",
                 fontsize=13, fontweight="bold", y=0.98)

    x = np.arange(len(FREQS))
    width = 0.32

    datasets = [
        ("EDT (s)", SIM_EDT, MEASURED_EDT, "Early Decay Time"),
        ("RT60 (s)", SIM_RT60, MEASURED_RT60, "Reverberation Time"),
        ("C80 (dB)", SIM_C80, MEASURED_C80, "Clarity"),
        ("D50 (%)", SIM_D50, MEASURED_D50, "Definition"),
    ]

    for ax, (ylabel, sim, meas, title) in zip(axes.flat, datasets):
        bars_sim = ax.bar(x - width / 2, sim, width, label="Simulated",
                          color="#4C72B0", edgecolor="white", linewidth=0.5)
        bars_meas = ax.bar(x + width / 2, meas, width, label="Measured",
                           color="#DD8452", edgecolor="white", linewidth=0.5)
        ax.set_xlabel("Frequency (Hz)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(FREQ_LABELS, fontsize=8)
        ax.tick_params(axis="y", labelsize=8)
        ax.legend(fontsize=8, loc="best")
        ax.grid(axis="y", alpha=0.3, linewidth=0.5)
        ax.set_axisbelow(True)

        # Add a zero line for C80
        if "C80" in ylabel:
            ax.axhline(y=0, color="gray", linewidth=0.5, linestyle="--")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = OUTPUT_DIR / "results.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {out}")


def parse_ac3d(path):
    """Parse an AC3D file and return lists of (vertices, faces, material_index)
    per object. Coordinates are in RAVEN y-up format."""
    text = Path(path).read_text()
    objects = []

    # Split by OBJECT poly
    parts = re.split(r"OBJECT poly", text)
    for part in parts[1:]:  # skip header
        # Parse vertices
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

        # Parse surfaces
        faces = []
        mat_idx = 0
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

    # Material colors (matching AC3D material order: concrete, windows, ceiling, plaster, floor)
    mat_colors = [
        (0.65, 0.64, 0.62, 0.25),   # concrete - gray
        (0.30, 0.65, 0.60, 0.30),   # windows - teal/transparent
        (0.75, 0.75, 0.75, 0.15),   # ceiling - light gray
        (0.70, 0.68, 0.62, 0.25),   # plaster - warm gray
        (0.60, 0.55, 0.42, 0.35),   # floor - brown
    ]

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Convert from RAVEN y-up (x, y_up, z) to building3d z-up (x, -z, y_up)
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
            edge_color = (0.3, 0.3, 0.3, 0.4)
            poly = Poly3DCollection([polygon], alpha=color[3],
                                    facecolor=color[:3], edgecolor=edge_color,
                                    linewidth=0.3)
            ax.add_collection3d(poly)

    # Plot sources
    for name, (sx, sy, sz) in SOURCES.items():
        ax.scatter(sx, sy, sz, c="red", s=80, marker="^", zorder=10,
                   edgecolors="darkred", linewidths=0.8, depthshade=False)
        ax.text(sx + 0.15, sy, sz + 0.15, name, fontsize=8, fontweight="bold",
                color="red")

    # Plot receivers
    for name, (rx, ry, rz) in RECEIVERS.items():
        ax.scatter(rx, ry, rz, c="blue", s=50, marker="o", zorder=10,
                   edgecolors="darkblue", linewidths=0.8, depthshade=False)
        ax.text(rx + 0.15, ry, rz + 0.1, name, fontsize=7, color="blue")

    ax.set_xlabel("x (m)", fontsize=9)
    ax.set_ylabel("y (m)", fontsize=9)
    ax.set_zlabel("z (m)", fontsize=9)
    ax.set_title("BRAS CR2: Seminar Room Geometry", fontsize=12, fontweight="bold")
    ax.tick_params(labelsize=7)

    # Set equal aspect ratio
    all_verts = []
    for verts_raw, _ in objects:
        all_verts.append(convert(verts_raw))
    all_verts = np.vstack(all_verts)
    mid = (all_verts.max(axis=0) + all_verts.min(axis=0)) / 2
    span = (all_verts.max(axis=0) - all_verts.min(axis=0)).max() / 2 * 1.1
    ax.set_xlim(mid[0] - span, mid[0] + span)
    ax.set_ylim(mid[1] - span, mid[1] + span)
    ax.set_zlim(mid[2] - span, mid[2] + span)

    ax.view_init(elev=25, azim=-55)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="^", color="w", markerfacecolor="red",
               markersize=10, label="Sources (LS1, LS2)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="blue",
               markersize=8, label="Receivers (MP1--MP5)"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc="upper left")

    out = OUTPUT_DIR / "geometry.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {out}")


if __name__ == "__main__":
    print("Generating BRAS CR2 figures...")
    generate_results_chart()
    generate_geometry_image()
    print("Done.")
