"""
Segmented growth rate analysis: early vs late phase.

Splits each well's growth curve at a midpoint and computes mu separately
for early and late phases to detect deceleration or crash.
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
    t0 = df["timestamp"].iloc[0]
    df["hours"] = (df["timestamp"] - t0).dt.total_seconds() / 3600
    df["minutes"] = df["hours"] * 60
    return df


def fit_mu(t, od):
    mask = np.isfinite(od) & (od > 0)
    t_m, od_m = np.array(t)[mask], np.array(od)[mask]
    if len(t_m) < 3:
        return np.nan, np.nan
    ln_od = np.log(od_m)
    slope, intercept, r, _, _ = linregress(t_m, ln_od)
    return slope, r**2


df_ce = load_plate("cellai_experiment_plate.csv")
df_ve = load_plate("vinatx_experiment_plate.csv")

cellai_wells = {
    "HBDef+Tryp+YE+Glut (B5)": "B5",
    "HBDef+Tryp+YE+Glut (B6)": "B6",
    "HBDef+Tryp+YE+Glut (B7)": "B7",
    "Semi-Def+Tryp+YE+MOPS (B2)": "B2",
    "Semi-Def+Tryp+YE+MOPS (B3)": "B3",
    "Semi-Def+Tryp+YE+MOPS (B4)": "B4",
    "LBv2+Gluc+KH2PO4+MOPS (B8)": "B8",
    "LBv2+Gluc+KH2PO4+MOPS (B9)": "B9",
    "LBv2+Gluc+KH2PO4+MOPS (B10)": "B10",
    "Cyclone ctrl (B11)": "B11",
}

vinatx_wells = {
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

# ── Compute phase-specific growth rates ──────────────────────────────────────

def compute_phases(df, wells_dict, team_name):
    """Split data into early (first half) and late (second half) and compute mu for each."""
    total_time = df["hours"].iloc[-1]
    mid = total_time / 2

    results = []
    for label, well in wells_dict.items():
        if well not in df.columns:
            continue

        t = df["hours"].values
        od = df[well].values

        # Full
        mu_full, r2_full = fit_mu(t, od)

        # Early: first half
        early_mask = t <= mid
        mu_early, r2_early = fit_mu(t[early_mask], od[early_mask])

        # Late: second half
        late_mask = t >= mid
        mu_late, r2_late = fit_mu(t[late_mask], od[late_mask])

        # Instantaneous: last 3 points
        if len(t) >= 3:
            mu_last3, r2_last3 = fit_mu(t[-3:], od[-3:])
        else:
            mu_last3, r2_last3 = np.nan, np.nan

        # Latest OD
        latest_od = od[-1]

        # Deceleration ratio
        decel = mu_late / mu_early if mu_early > 0 and np.isfinite(mu_late) else np.nan

        results.append({
            "team": team_name,
            "condition": label,
            "well": well,
            "mu_full": mu_full,
            "mu_early": mu_early,
            "mu_late": mu_late,
            "mu_last3": mu_last3,
            "decel_ratio": decel,
            "latest_od": latest_od,
            "mid_time_min": mid * 60,
        })

    return pd.DataFrame(results)


res_ce = compute_phases(df_ce, cellai_wells, "CellAI")
res_ve = compute_phases(df_ve, vinatx_wells, "ViNatX")
all_res = pd.concat([res_ce, res_ve], ignore_index=True)

# ── Print results ────────────────────────────────────────────────────────────

print(f"CellAI: {len(df_ce)} timepoints, {df_ce['minutes'].iloc[-1]:.0f} min total, split at {res_ce['mid_time_min'].iloc[0]:.0f} min")
print(f"ViNatX: {len(df_ve)} timepoints, {df_ve['minutes'].iloc[-1]:.0f} min total, split at {res_ve['mid_time_min'].iloc[0]:.0f} min")

for team in ["CellAI", "ViNatX"]:
    sub = all_res[all_res["team"] == team].sort_values("decel_ratio", ascending=True)
    print(f"\n{'='*90}")
    print(f"  {team} — Phase Analysis")
    print(f"{'='*90}")
    print(f"  {'Condition':>35s}  {'mu_early':>9s}  {'mu_late':>8s}  {'mu_last3':>9s}  {'decel':>7s}  {'OD_now':>7s}  {'Status'}")
    print(f"  {'─'*85}")
    for _, row in sub.iterrows():
        decel = row['decel_ratio']
        if decel < 0:
            status = "CRASHING"
        elif decel < 0.5:
            status = "DECELERATING"
        elif decel < 0.8:
            status = "Slowing"
        elif decel < 1.2:
            status = "Steady"
        else:
            status = "Accelerating"

        mu_e = f"{row['mu_early']:.3f}" if np.isfinite(row['mu_early']) else "N/A"
        mu_l = f"{row['mu_late']:.3f}" if np.isfinite(row['mu_late']) else "N/A"
        mu_3 = f"{row['mu_last3']:.3f}" if np.isfinite(row['mu_last3']) else "N/A"
        dec = f"{decel:.2f}" if np.isfinite(decel) else "N/A"
        print(f"  {row['condition']:>35s}  {mu_e:>9s}  {mu_l:>8s}  {mu_3:>9s}  {dec:>7s}  {row['latest_od']:>7.3f}  {status}")

# ── Figure: Phase comparison ─────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(18, 14))

# Panel 1: Early vs Late mu scatter
ax = axes[0, 0]
for team, color, marker in [("CellAI", "#2196F3", "o"), ("ViNatX", "#FF9800", "s")]:
    sub = all_res[all_res["team"] == team]
    ax.scatter(sub["mu_early"], sub["mu_late"], c=color, marker=marker, s=80,
               edgecolors="k", label=team, zorder=3)
    for _, row in sub.iterrows():
        ax.annotate(row["well"], (row["mu_early"], row["mu_late"]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
lims = [-0.2, 1.2]
ax.plot(lims, lims, "k--", alpha=0.3, label="No change")
ax.plot(lims, [0, 0], "r-", alpha=0.3)
ax.axhline(0, color="red", alpha=0.3, linestyle="-")
ax.set_xlabel("Early-phase mu (1/hr)")
ax.set_ylabel("Late-phase mu (1/hr)")
ax.set_title("Early vs Late Growth Rate\n(below diagonal = decelerating)")
ax.legend()
ax.grid(True, alpha=0.3)

# Panel 2: Deceleration ratio bar chart
ax = axes[0, 1]
all_sorted = all_res.sort_values("decel_ratio")
colors = []
for _, row in all_sorted.iterrows():
    d = row["decel_ratio"]
    if d < 0:
        colors.append("#d32f2f")
    elif d < 0.5:
        colors.append("#ff9800")
    elif d < 0.8:
        colors.append("#fdd835")
    elif d < 1.2:
        colors.append("#66bb6a")
    else:
        colors.append("#2196F3")

bars = ax.barh(
    [f"[{r['team'][:2]}] {r['condition']}" for _, r in all_sorted.iterrows()],
    all_sorted["decel_ratio"],
    color=colors, edgecolor="gray",
)
ax.axvline(1.0, color="k", linestyle="--", alpha=0.5, label="No change")
ax.axvline(0, color="red", linestyle="-", alpha=0.5, label="Zero growth")
ax.set_xlabel("Deceleration Ratio (late mu / early mu)")
ax.set_title("Growth Phase Ratio\n(red=crashing, orange=decelerating, green=steady, blue=accelerating)")
ax.legend(fontsize=8)
ax.grid(True, axis="x", alpha=0.3)

# Panel 3: CellAI growth curves with phase split
ax = axes[1, 0]
mid_ce = df_ce["hours"].iloc[-1] / 2
colors_ce = plt.cm.Blues(np.linspace(0.3, 0.9, len(cellai_wells)))
for i, (label, well) in enumerate(cellai_wells.items()):
    if well in df_ce.columns:
        ax.plot(df_ce["minutes"], df_ce[well], "o-", color=colors_ce[i],
                markersize=4, label=f"{label}", alpha=0.8)
ax.axvline(mid_ce * 60, color="red", linestyle="--", alpha=0.7, label=f"Split: {mid_ce*60:.0f} min")
ax.set_xlabel("Elapsed Time (min)")
ax.set_ylabel("OD600")
ax.set_title("CellAI — OD Growth Curves")
ax.legend(fontsize=6, loc="upper left")
ax.grid(True, alpha=0.3)

# Panel 4: ViNatX growth curves with phase split
ax = axes[1, 1]
mid_ve = df_ve["hours"].iloc[-1] / 2
colors_ve = plt.cm.Oranges(np.linspace(0.3, 0.9, len(vinatx_wells)))
for i, (label, well) in enumerate(vinatx_wells.items()):
    if well in df_ve.columns:
        ax.plot(df_ve["minutes"], df_ve[well], "s-", color=colors_ve[i],
                markersize=4, label=f"{label}", alpha=0.8)
ax.axvline(mid_ve * 60, color="red", linestyle="--", alpha=0.7, label=f"Split: {mid_ve*60:.0f} min")
ax.set_xlabel("Elapsed Time (min)")
ax.set_ylabel("OD600")
ax.set_title("ViNatX — OD Growth Curves")
ax.legend(fontsize=6, loc="upper left")
ax.grid(True, alpha=0.3)

plt.suptitle("Phase Analysis — Early vs Late Growth Rates\nDetecting deceleration and crash",
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(FIG_DIR / "phase_analysis.png", dpi=150, bbox_inches="tight")
plt.close()

all_res.to_csv(DATA_DIR / "phase_analysis.csv", index=False)
print(f"\nFigure saved: {FIG_DIR / 'phase_analysis.png'}")
