#!/usr/bin/env python3
"""Plot BESTEST ASHRAE 140 suite results vs EnergyPlus reference.

Usage:
    python plot_results.py                        # reads results.csv in same dir
    python plot_results.py path/to/results.csv    # explicit path

Produces:
    bestest_600_series.png    – 600-series annual heating & cooling (2-panel)
    bestest_900_series.png    – 900-series annual heating & cooling (2-panel)
    bestest_peaks.png         – peak loads comparison (2-panel)
    bestest_freefloat.png     – free-float min/max zone temperatures
    bestest_error_summary.png – percentage error overview (all HVAC cases)
    bestest_600_monthly.png   – monthly profile for case 600
    bestest_900_monthly.png   – monthly profile for case 900
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── Style ────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.facecolor": "white",
    "axes.facecolor": "#fafafa",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})

C_REF = "#888888"       # grey – reference
C_HEAT = "#d63031"      # red – heating
C_COOL = "#0984e3"      # blue – cooling
C_HEAT_LT = "#fab1a0"   # light red
C_COOL_LT = "#74b9ff"   # light blue
C_GOOD = "#00b894"      # green – within tolerance
C_BAD = "#e17055"       # orange – outside tolerance


# ── Helpers ──────────────────────────────────────────────────────────────────

def read_csv(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f))


def safe_float(s: str) -> float | None:
    if s is None or s.strip() == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def extract_annual(rows: list[dict], cases: list[str], metric: str):
    """Return (sim, ref) arrays for the given metric across cases."""
    sim, ref = [], []
    for case in cases:
        found = False
        for r in rows:
            if r["case"] == case and r["month"] == "annual" and r["metric"] == metric:
                sim.append(safe_float(r["simulated_kwh"]))
                ref.append(safe_float(r["reference_kwh"]))
                found = True
                break
        if not found:
            sim.append(None)
            ref.append(None)
    return sim, ref


def extract_monthly(rows: list[dict], case: str, metric: str):
    """Return (months, sim, ref) for a single case's monthly data."""
    months, sim, ref = [], [], []
    for r in rows:
        if r["case"] == case and r["metric"] == metric:
            m = r["month"]
            if m == "annual":
                continue
            months.append(int(m))
            sim.append(safe_float(r["simulated_kwh"]))
            ref.append(safe_float(r["reference_kwh"]))
    return months, sim, ref


def pct_err(sim, ref):
    if ref is None or sim is None or ref == 0:
        return None
    return 100.0 * (sim - ref) / ref


def bar_label(ax, bars, fmt="{:.0f}", fontsize=8):
    """Add value labels above/below bars."""
    for bar in bars:
        h = bar.get_height()
        if h == 0:
            continue
        ax.text(bar.get_x() + bar.get_width() / 2, h,
                fmt.format(h), ha="center", va="bottom", fontsize=fontsize)


MONTH_LABELS = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]


# ── Plot: two-panel annual loads (heating on top, cooling on bottom) ─────────

def plot_series_annual(rows, cases, series_label, out_path):
    sim_h, ref_h = extract_annual(rows, cases, "heating")
    sim_c, ref_c = extract_annual(rows, cases, "cooling")

    x = np.arange(len(cases))
    w = 0.35

    fig, (ax_h, ax_c) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    fig.suptitle(f"BESTEST {series_label}: Annual Loads vs E+ Reference",
                 fontsize=15, fontweight="bold", y=0.98)

    # Heating panel
    b1 = ax_h.bar(x - w / 2, [v or 0 for v in ref_h], w,
                  label="E+ Reference", color=C_REF, alpha=0.6)
    b2 = ax_h.bar(x + w / 2, [v or 0 for v in sim_h], w,
                  label="Simulated", color=C_HEAT)
    ax_h.set_ylabel("Heating (kWh)")
    ax_h.legend(loc="upper left")
    ax_h.set_axisbelow(True)

    # Error annotations
    for i in range(len(cases)):
        e = pct_err(sim_h[i], ref_h[i])
        if e is not None and ref_h[i] and ref_h[i] > 0:
            y_pos = max(sim_h[i] or 0, ref_h[i] or 0) + 50
            color = C_GOOD if abs(e) <= 10 else C_BAD
            ax_h.text(x[i], y_pos, f"{e:+.0f}%",
                      ha="center", va="bottom", fontsize=10,
                      fontweight="bold", color=color)

    # Cooling panel
    b3 = ax_c.bar(x - w / 2, [v or 0 for v in ref_c], w,
                  label="E+ Reference", color=C_REF, alpha=0.6)
    b4 = ax_c.bar(x + w / 2, [v or 0 for v in sim_c], w,
                  label="Simulated", color=C_COOL)
    ax_c.set_ylabel("Cooling (kWh)")
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([f"Case {c}" for c in cases], fontsize=11)
    ax_c.legend(loc="upper left")
    ax_c.set_axisbelow(True)

    for i in range(len(cases)):
        e = pct_err(sim_c[i], ref_c[i])
        if e is not None and ref_c[i] and ref_c[i] > 0:
            y_pos = max(sim_c[i] or 0, ref_c[i] or 0) + 50
            color = C_GOOD if abs(e) <= 10 else C_BAD
            ax_c.text(x[i], y_pos, f"{e:+.0f}%",
                      ha="center", va="bottom", fontsize=10,
                      fontweight="bold", color=color)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  {out_path.name}")


# ── Plot: peak loads (two-panel) ─────────────────────────────────────────────

def plot_peaks(rows, cases_600, cases_900, out_path):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("BESTEST: Peak Loads vs E+ Reference",
                 fontsize=15, fontweight="bold", y=0.99)

    panels = [
        (axes[0, 0], cases_600, "peak_heating_w", "600-Series Peak Heating (W)", C_HEAT),
        (axes[0, 1], cases_600, "peak_cooling_w", "600-Series Peak Cooling (W)", C_COOL),
        (axes[1, 0], cases_900, "peak_heating_w", "900-Series Peak Heating (W)", C_HEAT),
        (axes[1, 1], cases_900, "peak_cooling_w", "900-Series Peak Cooling (W)", C_COOL),
    ]

    w = 0.35
    for ax, cases, metric, title, color in panels:
        sim, ref = extract_annual(rows, cases, metric)
        x = np.arange(len(cases))

        ax.bar(x - w / 2, [v or 0 for v in ref], w,
               label="E+ Ref", color=C_REF, alpha=0.6)
        ax.bar(x + w / 2, [v or 0 for v in sim], w,
               label="Sim", color=color)

        ax.set_xticks(x)
        ax.set_xticklabels(cases, fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.legend(loc="upper left", fontsize=9)
        ax.set_axisbelow(True)

        for i in range(len(cases)):
            e = pct_err(sim[i], ref[i])
            if e is not None and ref[i] and ref[i] > 0:
                y_pos = max(sim[i] or 0, ref[i] or 0) + 50
                clr = C_GOOD if abs(e) <= 15 else C_BAD
                ax.text(x[i], y_pos, f"{e:+.0f}%",
                        ha="center", va="bottom", fontsize=9,
                        fontweight="bold", color=clr)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  {out_path.name}")


# ── Plot: free-float temperatures ────────────────────────────────────────────

def plot_freefloat(rows, cases, out_path):
    sim_min, ref_min = extract_annual(rows, cases, "min_zone_temp_c")
    sim_max, ref_max = extract_annual(rows, cases, "max_zone_temp_c")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_title("Free-Float Cases: Min/Max Zone Temperature",
                 fontsize=14, fontweight="bold")

    x = np.arange(len(cases))
    w = 0.18

    # Draw paired bars for min temps (left half) and max temps (right half)
    # Min temps
    ax.bar(x - w / 2 - 0.12, [v or 0 for v in ref_min], w,
           label="Ref Tmin", color=C_COOL_LT, edgecolor=C_COOL, linewidth=0.8)
    ax.bar(x + w / 2 - 0.12, [v or 0 for v in sim_min], w,
           label="Sim Tmin", color=C_COOL, alpha=0.8)
    # Max temps
    ax.bar(x - w / 2 + 0.12, [v or 0 for v in ref_max], w,
           label="Ref Tmax", color=C_HEAT_LT, edgecolor=C_HEAT, linewidth=0.8)
    ax.bar(x + w / 2 + 0.12, [v or 0 for v in sim_max], w,
           label="Sim Tmax", color=C_HEAT, alpha=0.8)

    # Delta annotations
    for i in range(len(cases)):
        if sim_min[i] is not None and ref_min[i] is not None:
            d = sim_min[i] - ref_min[i]
            y = min(sim_min[i], ref_min[i]) - 2.5
            ax.text(x[i] - 0.12, y, f"{d:+.1f}\u00b0",
                    ha="center", va="top", fontsize=9, color=C_COOL, fontweight="bold")
        if sim_max[i] is not None and ref_max[i] is not None:
            d = sim_max[i] - ref_max[i]
            y = max(sim_max[i], ref_max[i]) + 0.5
            ax.text(x[i] + 0.12, y, f"{d:+.1f}\u00b0",
                    ha="center", va="bottom", fontsize=9, color=C_HEAT, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"Case {c}" for c in cases], fontsize=11)
    ax.set_ylabel("Zone Temperature (\u00b0C)")
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.legend(loc="upper left", ncol=2, fontsize=9)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  {out_path.name}")


# ── Plot: error summary ─────────────────────────────────────────────────────

def plot_error_summary(rows, cases_600, cases_900, out_path):
    all_cases = cases_600 + cases_900

    fig, (ax_h, ax_c) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    fig.suptitle("Annual Load Error vs E+ Reference",
                 fontsize=15, fontweight="bold", y=0.99)

    for ax, metric, title, color in [
        (ax_h, "heating", "Heating Error (%)", C_HEAT),
        (ax_c, "cooling", "Cooling Error (%)", C_COOL),
    ]:
        sim, ref = extract_annual(rows, all_cases, metric)
        errors = []
        colors = []
        for i in range(len(all_cases)):
            e = pct_err(sim[i], ref[i])
            if e is None or (ref[i] is not None and ref[i] == 0):
                errors.append(0)
                colors.append("#cccccc")
            else:
                errors.append(e)
                colors.append(C_GOOD if abs(e) <= 10 else color)

        y = np.arange(len(all_cases))
        bars = ax.barh(y, errors, color=colors, edgecolor="white", linewidth=0.5)

        # Value labels
        for i, (bar, e) in enumerate(zip(bars, errors)):
            if e == 0 and ref[i] is not None and ref[i] == 0:
                continue  # skip zero-ref cases
            ha = "left" if e >= 0 else "right"
            offset = 1.0 if e >= 0 else -1.0
            ax.text(e + offset, y[i], f"{e:+.1f}%", ha=ha, va="center", fontsize=9)

        # Tolerance bands
        ax.axvspan(-10, 10, alpha=0.08, color="green", zorder=0)
        ax.axvline(0, color="black", linewidth=0.8)

        ax.set_yticks(y)
        ax.set_yticklabels([f"Case {c}" for c in all_cases], fontsize=10)
        ax.set_xlabel("Error (%)")
        ax.set_title(title, fontsize=12)
        ax.set_axisbelow(True)

        # Add separator between 600 and 900 series
        sep = len(cases_600) - 0.5
        ax.axhline(sep, color="grey", linewidth=1, linestyle="--", alpha=0.5)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  {out_path.name}")


# ── Plot: monthly heating & cooling for a single case ────────────────────────

def plot_monthly(rows, case, out_path):
    months_h, sim_h, ref_h = extract_monthly(rows, case, "heating")
    months_c, sim_c, ref_c = extract_monthly(rows, case, "cooling")

    if not months_h:
        return

    has_ref = any(v is not None for v in ref_h)

    fig, (ax_h, ax_c) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle(f"Case {case} \u2013 Monthly Heating & Cooling",
                 fontsize=14, fontweight="bold", y=0.98)

    x = np.arange(12)
    w = 0.35

    if has_ref:
        ax_h.bar(x - w / 2, [v or 0 for v in ref_h], w,
                 label="E+ Reference", color=C_REF, alpha=0.6)
    ax_h.bar(x + (w / 2 if has_ref else 0), [v or 0 for v in sim_h], w,
             label="Simulated", color=C_HEAT)
    ax_h.set_ylabel("Heating (kWh)")
    ax_h.legend(loc="upper right", fontsize=9)
    ax_h.set_axisbelow(True)

    if has_ref:
        ax_c.bar(x - w / 2, [v or 0 for v in ref_c], w,
                 label="E+ Reference", color=C_REF, alpha=0.6)
    ax_c.bar(x + (w / 2 if has_ref else 0), [v or 0 for v in sim_c], w,
             label="Simulated", color=C_COOL)
    ax_c.set_ylabel("Cooling (kWh)")
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(MONTH_LABELS, fontsize=11)
    ax_c.legend(loc="upper right", fontsize=9)
    ax_c.set_axisbelow(True)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  {out_path.name}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results.csv"
    if not csv_path.exists():
        print(f"Error: {csv_path} not found. Run bestest_energy_suite first.")
        sys.exit(1)

    out_dir = csv_path.parent
    rows = read_csv(csv_path)
    print(f"Read {len(rows)} rows from {csv_path}")

    cases_600 = ["600", "610", "620", "630", "640", "650"]
    cases_900 = ["900", "910", "920", "930", "940", "950", "960"]
    cases_ff = ["600FF", "650FF", "900FF", "950FF"]

    # 1-2. Annual loads by series (2-panel each)
    plot_series_annual(rows, cases_600, "600-Series (Lightweight)",
                       out_dir / "bestest_600_series.png")
    plot_series_annual(rows, cases_900, "900-Series (Heavyweight)",
                       out_dir / "bestest_900_series.png")

    # 3. Peak loads (2x2 grid)
    plot_peaks(rows, cases_600, cases_900, out_dir / "bestest_peaks.png")

    # 4. Free-float temperatures
    plot_freefloat(rows, cases_ff, out_dir / "bestest_freefloat.png")

    # 5. Error summary (side-by-side heating/cooling)
    plot_error_summary(rows, cases_600, cases_900,
                       out_dir / "bestest_error_summary.png")

    # 6-7. Monthly profiles (cases with and without reference)
    for case in ["600", "900", "610", "620", "630", "640", "650",
                 "910", "920", "930", "940", "950", "960"]:
        plot_monthly(rows, case, out_dir / f"bestest_{case}_monthly.png")

    print("Done.")


if __name__ == "__main__":
    main()
