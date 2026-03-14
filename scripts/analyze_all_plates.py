"""
Analyze OD600 growth data from all 4 experiment plates (2 per team).

Plates:
  - CellAI Tutorial Plate: seeding density experiment (Cyclone media, 100-190 uL, 3 reps)
  - CellAI Experiment Plate: media composition optimization (various supplements)
  - ViNatX Tutorial Plate: seeding density experiment (CycloneX media, 100-190 uL, 3 reps)
  - ViNatX Experiment Plate: media composition optimization (various supplements)

Computes specific growth rate (mu) via ln(OD) linear regression and generates figures.
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ── helpers ──────────────────────────────────────────────────────────────────

def load_plate(csv_name):
    df = pd.read_csv(DATA_DIR / csv_name, parse_dates=["timestamp"])
    df = df.sort_values("timestamp").reset_index(drop=True)
    wells = [c for c in df.columns if c != "timestamp"]
    # compute elapsed hours from first reading
    t0 = df["timestamp"].iloc[0]
    df["hours"] = (df["timestamp"] - t0).dt.total_seconds() / 3600
    return df, wells


def fit_growth_rate(times, od_values):
    """Fit mu (1/hr) from ln(OD) vs time via linear regression."""
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan, np.nan
    ln_od = np.log(od)
    slope, intercept, r, p, se = linregress(t, ln_od)
    return slope, r**2, np.log(2) / slope * 60 if slope > 0 else np.nan  # doubling time in min


# ── Plate 1: CellAI Tutorial (seeding density) ──────────────────────────────

def analyze_cellai_tutorial():
    df, wells = load_plate("cellai_tutorial_plate.csv")
    # Conditions: wells in rows A/B/C, columns 1-10
    # Volumes: 100, 110, 120, 130, 140, 150, 160, 170, 180, 190 uL media
    # => seeding %: (200 - vol) / 200 * 100
    vol_map = {1: 100, 2: 110, 3: 120, 4: 130, 5: 140, 6: 150, 7: 160, 8: 170, 9: 180, 10: 190}
    conditions = {}
    for col_num, vol in vol_map.items():
        seed_pct = (200 - vol) / 200 * 100
        label = f"{seed_pct:.0f}%"
        reps = [f"{row}{col_num}" for row in "ABC"]
        conditions[label] = reps

    results = []
    for label, reps in conditions.items():
        mus = []
        for w in reps:
            if w in df.columns:
                mu, r2, td = fit_growth_rate(df["hours"], df[w])
                mus.append(mu)
                results.append({"plate": "CellAI Tutorial", "condition": label, "well": w,
                                "mu_per_hr": mu, "R2": r2, "doubling_min": td})

    return pd.DataFrame(results), df, conditions


def analyze_vinatx_tutorial():
    df, wells = load_plate("vinatx_tutorial_plate.csv")
    # Rows A-F, columns 2-4 (3 replicates)
    # Volumes: A=100, B=140, C=160, D=170, E=180, F=190 uL
    vol_map = {"A": 100, "B": 140, "C": 160, "D": 170, "E": 180, "F": 190}
    conditions = {}
    for row, vol in vol_map.items():
        seed_pct = (200 - vol) / 200 * 100
        label = f"{seed_pct:.0f}%"
        reps = [f"{row}{c}" for c in [2, 3, 4]]
        conditions[label] = reps

    results = []
    for label, reps in conditions.items():
        for w in reps:
            if w in df.columns:
                mu, r2, td = fit_growth_rate(df["hours"], df[w])
                results.append({"plate": "ViNatX Tutorial", "condition": label, "well": w,
                                "mu_per_hr": mu, "R2": r2, "doubling_min": td})

    return pd.DataFrame(results), df, conditions


# ── Plate 2: CellAI Experiment (media composition) ──────────────────────────

def analyze_cellai_experiment():
    df, wells = load_plate("cellai_experiment_plate.csv")
    # Each well is a unique media condition
    conditions = {
        "Semi-Def+Tryp+YE+MOPS (B2)": ["B2"],   # Semi-Defined 129 + Tryptone 26 + YE 14 + MOPS 11
        "Semi-Def+Tryp+YE+MOPS (B3)": ["B3"],   # Semi-Defined 128 + Tryptone 14 + YE 20 + MOPS 7 + Water 11
        "Semi-Def+Tryp+YE+MOPS (B4)": ["B4"],   # Semi-Defined 125 + Tryptone 20 + YE 26 + MOPS 9
        "HBDef+Tryp+YE+Glut (B5)": ["B5"],      # High Buffer Defined 110 + Tryptone 26 + YE 14 + Glutamate 10 + Water 20
        "HBDef+Tryp+YE+Glut (B6)": ["B6"],      # High Buffer Defined 110 + Tryptone 14 + YE 20 + Glutamate 6 + Water 30
        "HBDef+Tryp+YE+Glut (B7)": ["B7"],      # High Buffer Defined 110 + Tryptone 20 + YE 26 + Glutamate 8 + Water 16
        "LBv2+Gluc+KH2PO4+MOPS (B8)": ["B8"],   # LBv2 130 + Glucose 18 + KH2PO4 10 + MOPS 12 + Water 10
        "LBv2+Gluc+KH2PO4+MOPS (B9)": ["B9"],   # LBv2 130 + Glucose 12 + KH2PO4 6 + MOPS 16 + Water 16
        "LBv2+Gluc+KH2PO4+MOPS (B10)": ["B10"], # LBv2 130 + Glucose 8 + KH2PO4 8 + MOPS 18 + Water 16
        "Novel Bio Cyclone (B11)": ["B11"],       # Control: Cyclone 180
    }

    results = []
    for label, reps in conditions.items():
        for w in reps:
            if w in df.columns:
                mu, r2, td = fit_growth_rate(df["hours"], df[w])
                results.append({"plate": "CellAI Experiment", "condition": label, "well": w,
                                "mu_per_hr": mu, "R2": r2, "doubling_min": td})

    return pd.DataFrame(results), df, conditions


def analyze_vinatx_experiment():
    df, wells = load_plate("vinatx_experiment_plate.csv")
    # Each well is a unique media condition
    conditions = {
        "Def-Min (A1)": ["A1"],                    # Defined-Minimal 180
        "Def-Min+Gluc (A2)": ["A2"],               # Defined-Minimal 160 + Glucose 20
        "Def-Min+NaCl (A3)": ["A3"],               # Defined-Minimal 160 + NaCl 20
        "Def-Min+Gluc+NaCl (A4)": ["A4"],          # Defined-Minimal 140 + Glucose 20 + NaCl 20
        "Def-Min+MOPS (B1)": ["B1"],               # Defined-Minimal 160 + MOPS 20
        "Def-Min+Gluc+MOPS (B2)": ["B2"],          # Defined-Minimal 140 + Glucose 20 + MOPS 20
        "Def-Min+MOPS+NaCl (B3)": ["B3"],          # Defined-Minimal 140 + MOPS 20 + NaCl 20
        "Def-Min+Gluc+MOPS+NaCl (B4)": ["B4"],     # Defined-Minimal 120 + Glucose 20 + MOPS 20 + NaCl 20
        "NBxCyclone (C1)": ["C1"],                  # NBxCyclone 180
        "LBv2 (C2)": ["C2"],                        # LBv2 180
        "Semi-Defined (C3)": ["C3"],                 # Semi-Defined 180
        "Def-Glycerol (C4)": ["C4"],                 # Defined-Glycerol 180
    }

    results = []
    for label, reps in conditions.items():
        for w in reps:
            if w in df.columns:
                mu, r2, td = fit_growth_rate(df["hours"], df[w])
                results.append({"plate": "ViNatX Experiment", "condition": label, "well": w,
                                "mu_per_hr": mu, "R2": r2, "doubling_min": td})

    return pd.DataFrame(results), df, conditions


# ── Plotting ─────────────────────────────────────────────────────────────────

def plot_growth_curves_tutorial(df_cellai, cond_cellai, df_vinatx, cond_vinatx):
    """Growth curves for tutorial (seeding density) plates."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)

    colors = plt.cm.viridis(np.linspace(0, 0.9, max(len(cond_cellai), len(cond_vinatx))))

    for ax, df, conditions, title in [
        (axes[0], df_cellai, cond_cellai, "CellAI Tutorial\n(Seeding Density)"),
        (axes[1], df_vinatx, cond_vinatx, "ViNatX Tutorial\n(Seeding Density)"),
    ]:
        for i, (label, reps) in enumerate(conditions.items()):
            for j, w in enumerate(reps):
                if w in df.columns:
                    style = "-" if j == 0 else ("--" if j == 1 else ":")
                    ax.plot(df["hours"], df[w], style, color=colors[i],
                            label=label if j == 0 else None, alpha=0.8)
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("OD600")
        ax.set_title(title)
        ax.legend(fontsize=7, title="Seeding %")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "tutorial_growth_curves.png", dpi=150)
    plt.close()


def plot_growth_curves_experiment(df_cellai, cond_cellai, df_vinatx, cond_vinatx):
    """Growth curves for experiment (media composition) plates."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # CellAI experiment
    ax = axes[0]
    colors_c = plt.cm.tab10(np.linspace(0, 1, len(cond_cellai)))
    for i, (label, reps) in enumerate(cond_cellai.items()):
        for w in reps:
            if w in df_cellai.columns:
                ax.plot(df_cellai["hours"], df_cellai[w], "o-", color=colors_c[i],
                        label=label, markersize=3, alpha=0.8)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("OD600")
    ax.set_title("CellAI Experiment\n(Media Composition)")
    ax.legend(fontsize=6, loc="upper left")
    ax.grid(True, alpha=0.3)

    # ViNatX experiment
    ax = axes[1]
    colors_v = plt.cm.tab20(np.linspace(0, 1, len(cond_vinatx)))
    for i, (label, reps) in enumerate(cond_vinatx.items()):
        for w in reps:
            if w in df_vinatx.columns:
                ax.plot(df_vinatx["hours"], df_vinatx[w], "o-", color=colors_v[i],
                        label=label, markersize=3, alpha=0.8)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("OD600")
    ax.set_title("ViNatX Experiment\n(Media Composition)")
    ax.legend(fontsize=6, loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "experiment_growth_curves.png", dpi=150)
    plt.close()


def plot_growth_rate_comparison(all_results):
    """Bar chart comparing growth rates across all plates and conditions."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    plates = ["CellAI Tutorial", "ViNatX Tutorial", "CellAI Experiment", "ViNatX Experiment"]

    for ax, plate_name in zip(axes.flat, plates):
        sub = all_results[all_results["plate"] == plate_name].copy()
        if sub.empty:
            ax.set_title(f"{plate_name}\n(no data)")
            continue

        # group by condition
        grouped = sub.groupby("condition")["mu_per_hr"].agg(["mean", "std"]).reset_index()
        grouped = grouped.sort_values("mean", ascending=True)

        colors = plt.cm.RdYlGn(np.linspace(0.2, 0.9, len(grouped)))
        bars = ax.barh(grouped["condition"], grouped["mean"], xerr=grouped["std"],
                       color=colors, edgecolor="gray", capsize=3)
        ax.set_xlabel("Growth Rate mu (1/hr)")
        ax.set_title(plate_name)
        ax.grid(True, axis="x", alpha=0.3)

        # annotate
        for bar, val in zip(bars, grouped["mean"]):
            if np.isfinite(val):
                ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                        f"{val:.3f}", va="center", fontsize=8)

    plt.suptitle("Growth Rate Comparison Across All Plates", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "growth_rate_comparison.png", dpi=150)
    plt.close()


def plot_ln_od_fits(df, conditions, plate_name, ax):
    """Plot ln(OD) with linear fits on a given axis."""
    colors = plt.cm.tab10(np.linspace(0, 1, min(len(conditions), 10)))
    for i, (label, reps) in enumerate(conditions.items()):
        c = colors[i % len(colors)]
        for w in reps:
            if w not in df.columns:
                continue
            od = df[w].values
            t = df["hours"].values
            mask = np.isfinite(od) & (od > 0)
            if mask.sum() < 3:
                continue
            ln_od = np.log(od[mask])
            t_m = t[mask]
            ax.scatter(t_m, ln_od, color=c, s=15, alpha=0.6)
            slope, intercept, _, _, _ = linregress(t_m, ln_od)
            ax.plot(t_m, slope * t_m + intercept, "-", color=c, alpha=0.8,
                    label=f"{label}: mu={slope:.3f}")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("ln(OD600)")
    ax.set_title(plate_name)
    ax.legend(fontsize=5, loc="lower right")
    ax.grid(True, alpha=0.3)


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Analyze all plates
    res_ct, df_ct, cond_ct = analyze_cellai_tutorial()
    res_vt, df_vt, cond_vt = analyze_vinatx_tutorial()
    res_ce, df_ce, cond_ce = analyze_cellai_experiment()
    res_ve, df_ve, cond_ve = analyze_vinatx_experiment()

    all_results = pd.concat([res_ct, res_vt, res_ce, res_ve], ignore_index=True)
    all_results.to_csv(DATA_DIR / "all_growth_rates.csv", index=False)

    # Print summary
    print("=" * 80)
    print("GROWTH RATE ANALYSIS — ALL PLATES")
    print("=" * 80)

    for plate_name in ["CellAI Tutorial", "ViNatX Tutorial", "CellAI Experiment", "ViNatX Experiment"]:
        sub = all_results[all_results["plate"] == plate_name]
        print(f"\n{'─' * 60}")
        print(f"  {plate_name}")
        print(f"{'─' * 60}")
        grouped = sub.groupby("condition")["mu_per_hr"].agg(["mean", "std"]).reset_index()
        grouped = grouped.sort_values("mean", ascending=False)
        for _, row in grouped.iterrows():
            std_str = f" +/- {row['std']:.4f}" if pd.notna(row['std']) and row['std'] > 0 else ""
            print(f"  {row['condition']:>35s}:  mu = {row['mean']:.4f}{std_str} (1/hr)")

    # Find best overall
    print(f"\n{'=' * 80}")
    print("TOP 5 GROWTH RATES (by condition mean)")
    print(f"{'=' * 80}")
    top = all_results.groupby(["plate", "condition"])["mu_per_hr"].mean().reset_index()
    top = top.sort_values("mu_per_hr", ascending=False).head(5)
    for _, row in top.iterrows():
        print(f"  {row['plate']:>20s} | {row['condition']:>35s} | mu = {row['mu_per_hr']:.4f} 1/hr")

    # Generate all plots
    plot_growth_curves_tutorial(df_ct, cond_ct, df_vt, cond_vt)
    plot_growth_curves_experiment(df_ce, cond_ce, df_ve, cond_ve)
    plot_growth_rate_comparison(all_results)

    # ln(OD) fit plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    plot_ln_od_fits(df_ct, cond_ct, "CellAI Tutorial (Seeding Density)", axes[0, 0])
    plot_ln_od_fits(df_vt, cond_vt, "ViNatX Tutorial (Seeding Density)", axes[0, 1])
    plot_ln_od_fits(df_ce, cond_ce, "CellAI Experiment (Media Comp.)", axes[1, 0])
    plot_ln_od_fits(df_ve, cond_ve, "ViNatX Experiment (Media Comp.)", axes[1, 1])
    plt.suptitle("ln(OD600) Fits — Growth Rate Determination", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "all_ln_od_fits.png", dpi=150)
    plt.close()

    print(f"\nFigures saved to {FIG_DIR}/")
    print("  - tutorial_growth_curves.png")
    print("  - experiment_growth_curves.png")
    print("  - growth_rate_comparison.png")
    print("  - all_ln_od_fits.png")
    print(f"Data saved to {DATA_DIR / 'all_growth_rates.csv'}")
