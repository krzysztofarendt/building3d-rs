#!/usr/bin/env python3
"""Generate a comparison chart for the BESTEST 600 energy benchmark.

Usage:
  cargo run --example bestest_600_energy --release
  python3 examples/bestest_600_energy/generate_figures.py
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).parent
RESULTS = SCRIPT_DIR / "results.csv"


def load_rows(path: Path):
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows


def main():
    if not RESULTS.exists():
        raise SystemExit(f"Missing {RESULTS}. Run the Rust example first.")

    rows = load_rows(RESULTS)
    months = [r for r in rows if r["month"] != "annual"]

    def series(metric: str):
        ms = [r for r in months if r["metric"] == metric]
        ms.sort(key=lambda r: int(r["month"]))
        sim = np.array([float(r["simulated_kwh"]) for r in ms])
        ref = np.array([float(r["reference_kwh"]) for r in ms])
        return sim, ref

    sim_h, ref_h = series("heating")
    sim_c, ref_c = series("cooling")

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle("BESTEST 600: building3d vs OpenStudio/EnergyPlus (monthly loads)",
                 fontsize=12, fontweight="bold", y=0.98)

    x = np.arange(12)
    width = 0.35

    axes[0].bar(x - width / 2, sim_h, width, label="building3d", color="#4C72B0")
    axes[0].bar(x + width / 2, ref_h, width, label="OpenStudio/E+", color="#DD8452")
    axes[0].set_ylabel("Heating (kWh)")
    axes[0].grid(axis="y", alpha=0.3)
    axes[0].legend(loc="best", fontsize=9)

    axes[1].bar(x - width / 2, sim_c, width, label="building3d", color="#4C72B0")
    axes[1].bar(x + width / 2, ref_c, width, label="OpenStudio/E+", color="#DD8452")
    axes[1].set_ylabel("Cooling (kWh)")
    axes[1].grid(axis="y", alpha=0.3)

    axes[1].set_xticks(x)
    axes[1].set_xticklabels(
        ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    )

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = SCRIPT_DIR / "results.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()

