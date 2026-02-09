#!/usr/bin/env python3
"""Generate figures for the BESTEST energy suite.

Usage:
  bash examples/bestest_energy_suite/download_data.sh
  cargo run --example bestest_energy_suite --release
  python3 examples/bestest_energy_suite/generate_figures.py
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).parent
RESULTS = SCRIPT_DIR / "results.csv"

MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def load_rows(path: Path):
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)
    return rows


def series(rows, case: str, metric: str):
    rs = [r for r in rows if r["case"] == case and r["metric"] == metric and r["month"] != "annual"]
    rs.sort(key=lambda r: int(r["month"]))
    sim = np.array([float(r["simulated_kwh"]) for r in rs])
    ref = np.array([float(r["reference_kwh"]) for r in rs])
    return sim, ref


def annual_row(rows, case: str, metric: str):
    for r in rows:
        if r["case"] == case and r["metric"] == metric and r["month"] == "annual":
            return r
    return None


def fmt_pct(x: float) -> str:
    if not np.isfinite(x):
        return "n/a"
    return f"{x:+.1f}%"


def main():
    if not RESULTS.exists():
        raise SystemExit(f"Missing {RESULTS}. Run the Rust example first.")

    rows = load_rows(RESULTS)

    cases = ["600", "900"]
    metrics = ["heating", "cooling"]

    fig, axes = plt.subplots(len(cases), len(metrics), figsize=(12, 7), sharex=True)
    fig.suptitle("BESTEST Energy Suite: building3d vs OpenStudio/EnergyPlus (monthly loads)",
                 fontsize=12, fontweight="bold", y=0.98)

    x = np.arange(12)
    width = 0.36

    for i, case in enumerate(cases):
        for j, metric in enumerate(metrics):
            ax = axes[i, j] if len(cases) > 1 else axes[j]
            sim, ref = series(rows, case, metric)

            ax.bar(x - width / 2, sim, width, label="building3d", color="#4C72B0")
            ax.bar(x + width / 2, ref, width, label="OpenStudio/E+", color="#DD8452")

            ar = annual_row(rows, case, metric)
            pct = float(ar["error_pct"]) if ar and ar["error_pct"] else float("nan")
            ax.set_title(f"Case {case} — {metric} ({fmt_pct(pct)} annual)", fontsize=10, fontweight="bold")

            ax.set_ylabel("kWh")
            ax.grid(axis="y", alpha=0.3)
            ax.set_axisbelow(True)
            ax.set_xticks(x)
            ax.set_xticklabels(MONTH_LABELS, fontsize=8)

            if i == 0 and j == 0:
                ax.legend(loc="best", fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = SCRIPT_DIR / "results.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()

