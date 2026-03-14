"""
Fair growth rate comparison between CellAI and ViNatX experiment plates.

Computes mu using equal time windows so growth rates are directly comparable.
ViNatX has ~20 min of data (3 timepoints), so we truncate CellAI to its
first 3 timepoints (~20 min) as well. We also show mu computed over the
first ~30 min (4 timepoints) for CellAI since that's still close.
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)


def load_plate(csv_name):
    df = pd.read_csv(DATA_DIR / csv_name, parse_dates=["timestamp"])
    df = df.sort_values("timestamp").reset_index(drop=True)
    wells = [c for c in df.columns if c != "timestamp"]
    t0 = df["timestamp"].iloc[0]
    df["hours"] = (df["timestamp"] - t0).dt.total_seconds() / 3600
    return df, wells


def fit_growth_rate(times, od_values):
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    ln_od = np.log(od)
    slope, intercept, r, p, se = linregress(t, ln_od)
    return slope, r**2


# ── Load data ────────────────────────────────────────────────────────────────

df_ce, _ = load_plate("cellai_experiment_plate.csv")
df_ve, _ = load_plate("vinatx_experiment_plate.csv")

print(f"CellAI Experiment: {len(df_ce)} timepoints, span = {df_ce['hours'].iloc[-1]:.2f} hrs ({df_ce['hours'].iloc[-1]*60:.0f} min)")
print(f"ViNatX Experiment: {len(df_ve)} timepoints, span = {df_ve['hours'].iloc[-1]:.2f} hrs ({df_ve['hours'].iloc[-1]*60:.0f} min)")

# ── Define conditions ────────────────────────────────────────────────────────

cellai_conditions = {
    "Semi-Def+Tryp+YE+MOPS (B2)": "B2",
    "Semi-Def+Tryp+YE+MOPS (B3)": "B3",
    "Semi-Def+Tryp+YE+MOPS (B4)": "B4",
    "HBDef+Tryp+YE+Glut (B5)": "B5",
    "HBDef+Tryp+YE+Glut (B6)": "B6",
    "HBDef+Tryp+YE+Glut (B7)": "B7",
    "LBv2+Gluc+KH2PO4+MOPS (B8)": "B8",
    "LBv2+Gluc+KH2PO4+MOPS (B9)": "B9",
    "LBv2+Gluc+KH2PO4+MOPS (B10)": "B10",
    "Novel Bio Cyclone (B11)": "B11",
}

vinatx_conditions = {
    "Def-Min (A1)": "A1",
    "Def-Min+Gluc (A2)": "A2",
    "Def-Min+NaCl (A3)": "A3",
    "Def-Min+Gluc+NaCl (A4)": "A4",
    "Def-Min+MOPS (B1)": "B1",
    "Def-Min+Gluc+MOPS (B2)": "B2",
    "Def-Min+MOPS+NaCl (B3)": "B3",
    "Def-Min+Gluc+MOPS+NaCl (B4)": "B4",
    "NBxCyclone (C1)": "C1",
    "LBv2 (C2)": "C2",
    "Semi-Defined (C3)": "C3",
    "Def-Glycerol (C4)": "C4",
}

# ── Compute mu at different time windows ─────────────────────────────────────

# Match to ViNatX's span (use all CellAI data within that window, gap or not)
import math
vinatx_span_hrs = df_ve["hours"].iloc[-1]
match_hrs = math.ceil(vinatx_span_hrs * 60) / 60
print(f"\nViNatX span: {vinatx_span_hrs*60:.1f} min")
print(f"Matching time window: {match_hrs*60:.1f} min ({match_hrs:.3f} hrs)")

# Truncate both to the matching window
cellai_mask = df_ce["hours"] <= match_hrs + 0.001
n_matched = cellai_mask.sum()
print(f"CellAI timepoints within window: {n_matched} (of {len(df_ce)})")

df_ce_trunc = df_ce[cellai_mask].copy()
df_ve_trunc = df_ve[df_ve["hours"] <= match_hrs + 0.001].copy()

print(f"ViNatX timepoints within window: {len(df_ve_trunc)} (of {len(df_ve)})")

results = []

# ViNatX — truncated to matching window
for label, well in vinatx_conditions.items():
    mu_matched, r2_matched = fit_growth_rate(df_ve_trunc["hours"], df_ve_trunc[well])
    mu_full, r2_full = fit_growth_rate(df_ve["hours"], df_ve[well])
    results.append({
        "team": "ViNatX",
        "condition": label,
        "well": well,
        "mu_matched": mu_matched,
        "R2_matched": r2_matched,
        "mu_full": mu_full,
        "R2_full": r2_full,
        "window_min": match_hrs * 60,
    })

# CellAI — both truncated (matched) and full
for label, well in cellai_conditions.items():
    mu_matched, r2_matched = fit_growth_rate(df_ce_trunc["hours"], df_ce_trunc[well])
    mu_full, r2_full = fit_growth_rate(df_ce["hours"], df_ce[well])
    results.append({
        "team": "CellAI",
        "condition": label,
        "well": well,
        "mu_matched": mu_matched,
        "R2_matched": r2_matched,
        "mu_full": mu_full,
        "R2_full": r2_full,
        "window_min": match_hrs * 60,
    })

res_df = pd.DataFrame(results)
res_df.to_csv(DATA_DIR / "equal_time_growth_rates.csv", index=False)

# ── Print comparison ─────────────────────────────────────────────────────────

print(f"\n{'='*85}")
print(f"FAIR COMPARISON — mu computed over first {match_hrs*60:.0f} min for both teams")
print(f"{'='*85}")

for team in ["CellAI", "ViNatX"]:
    sub = res_df[res_df["team"] == team].sort_values("mu_matched", ascending=False)
    print(f"\n  {team}")
    print(f"  {'Condition':>40s}  {'mu(matched)':>12s}  {'R2':>6s}  {'mu(full)':>10s}")
    print(f"  {'─'*75}")
    for _, row in sub.iterrows():
        mu_m = f"{row['mu_matched']:.4f}" if np.isfinite(row['mu_matched']) else "N/A"
        r2_m = f"{row['R2_matched']:.3f}" if np.isfinite(row['R2_matched']) else "N/A"
        mu_f = f"{row['mu_full']:.4f}" if np.isfinite(row['mu_full']) else "N/A"
        print(f"  {row['condition']:>40s}  {mu_m:>12s}  {r2_m:>6s}  {mu_f:>10s}")

# ── Rank all conditions together ─────────────────────────────────────────────

print(f"\n{'='*85}")
print(f"OVERALL RANKING (matched window, {match_hrs*60:.0f} min)")
print(f"{'='*85}")

ranked = res_df.sort_values("mu_matched", ascending=False)
for i, (_, row) in enumerate(ranked.iterrows(), 1):
    mu_m = f"{row['mu_matched']:.4f}" if np.isfinite(row['mu_matched']) else "N/A"
    print(f"  {i:2d}. [{row['team']:>6s}] {row['condition']:>40s}  mu = {mu_m} 1/hr")

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(12, 8))

ranked = res_df.dropna(subset=["mu_matched"]).sort_values("mu_matched", ascending=True)

colors = []
for team in ranked["team"]:
    colors.append("#2196F3" if team == "CellAI" else "#FF9800")

bars = ax.barh(
    [f"[{r['team']}] {r['condition']}" for _, r in ranked.iterrows()],
    ranked["mu_matched"],
    color=colors,
    edgecolor="gray",
    alpha=0.85,
)

# Annotate bars
for bar, (_, row) in zip(bars, ranked.iterrows()):
    if np.isfinite(row["mu_matched"]):
        ax.text(
            bar.get_width() + 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{row['mu_matched']:.3f}",
            va="center",
            fontsize=8,
        )

# Also overlay full-time mu as markers for CellAI
cellai_ranked = ranked[ranked["team"] == "CellAI"]
if not cellai_ranked.empty:
    y_positions = []
    for _, r in cellai_ranked.iterrows():
        label = f"[{r['team']}] {r['condition']}"
        idx = list(ranked.index)
        pos = list(ranked.index).index(r.name)
        y_positions.append(pos)
    ax.scatter(
        cellai_ranked["mu_full"],
        [f"[{r['team']}] {r['condition']}" for _, r in cellai_ranked.iterrows()],
        color="red",
        marker="D",
        s=30,
        zorder=5,
        label=f"CellAI mu (full {df_ce['hours'].iloc[-1]*60:.0f} min)",
    )

ax.set_xlabel("Growth Rate mu (1/hr)")
ax.set_title(
    f"Fair Growth Rate Comparison — Equal Time Window ({match_hrs*60:.0f} min)\n"
    f"Blue = CellAI | Orange = ViNatX | Red diamonds = CellAI full-time mu",
    fontsize=11,
)
ax.legend(loc="lower right", fontsize=9)
ax.grid(True, axis="x", alpha=0.3)

plt.tight_layout()
plt.savefig(FIG_DIR / "equal_time_comparison.png", dpi=150)
plt.close()
print(f"\nFigure saved: {FIG_DIR / 'equal_time_comparison.png'}")
