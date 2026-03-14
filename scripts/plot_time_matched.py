"""
Time-matched growth rate plots: both CellAI and ViNatX truncated to the
same time window (ViNatX's span) so slopes are visually comparable.
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
    return df


def fit_mu(t, od):
    mask = np.isfinite(od) & (od > 0)
    t_m, od_m = np.array(t)[mask], np.array(od)[mask]
    if len(t_m) < 3:
        return np.nan, np.nan, np.nan
    ln_od = np.log(od_m)
    slope, intercept, r, _, _ = linregress(t_m, ln_od)
    return slope, intercept, r**2


df_ce = load_plate("cellai_experiment_plate.csv")
df_ve = load_plate("vinatx_experiment_plate.csv")

# Use ViNatX's full span, and for CellAI find the next available timepoint
# after ViNatX's span (to bridge CellAI's ~40 min gap)
import math
vinatx_span = df_ve["hours"].iloc[-1]
ve_match_hrs = math.ceil(vinatx_span * 60) / 60

# For CellAI: find the first timepoint >= ViNatX's span
ce_candidates = df_ce[df_ce["hours"] >= vinatx_span - 0.01]["hours"]
if len(ce_candidates) > 0:
    ce_match_hrs = ce_candidates.iloc[0] + 0.001  # include that point
else:
    ce_match_hrs = ve_match_hrs
match_hrs = ve_match_hrs  # for axis labels
match_min = match_hrs * 60

df_ce_m = df_ce[df_ce["hours"] <= ce_match_hrs].copy()
df_ve_m = df_ve[df_ve["hours"] <= ve_match_hrs + 0.001].copy()

print(f"Matched window: ViNatX={ve_match_hrs*60:.0f} min ({len(df_ve_m)} pts) | CellAI={df_ce_m['hours'].iloc[-1]*60:.0f} min ({len(df_ce_m)} pts)")
match_min = ve_match_hrs * 60

cellai_conds = [
    ("Semi-Def+Tryp+YE+MOPS", ["B2", "B3", "B4"]),
    ("HBDef+Tryp+YE+Glut", ["B5", "B6", "B7"]),
    ("LBv2+Gluc+KH2PO4+MOPS", ["B8", "B9", "B10"]),
    ("Cyclone (ctrl)", ["B11"]),
]

vinatx_conds = [
    ("Def-Min", ["A1"]),
    ("Def-Min+Gluc", ["A2"]),
    ("Def-Min+NaCl", ["A3"]),
    ("Def-Min+Gluc+NaCl", ["A4"]),
    ("Def-Min+MOPS", ["B1"]),
    ("Def-Min+Gluc+MOPS", ["B2"]),
    ("Def-Min+MOPS+NaCl", ["B3"]),
    ("Def-Min+Gluc+MOPS+NaCl", ["B4"]),
    ("NBxCyclone", ["C1"]),
    ("LBv2", ["C2"]),
    ("Semi-Defined", ["C3"]),
    ("Def-Glycerol", ["C4"]),
]

# ── Compute matched mu for all conditions ────────────────────────────────────

all_mus = []
for label, wells in cellai_conds:
    for w in wells:
        if w in df_ce_m.columns:
            mu, inter, r2 = fit_mu(df_ce_m["hours"], df_ce_m[w])
            all_mus.append(("CellAI", f"{label} ({w})", mu, r2, w, inter))

for label, wells in vinatx_conds:
    for w in wells:
        if w in df_ve_m.columns:
            mu, inter, r2 = fit_mu(df_ve_m["hours"], df_ve_m[w])
            all_mus.append(("ViNatX", f"{label} ({w})", mu, r2, w, inter))

all_mus.sort(key=lambda x: -x[2] if np.isfinite(x[2]) else -999)

# ── Figure: 2x2 layout ──────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

ce_colors = plt.cm.Blues(np.linspace(0.3, 0.9, 10))
ve_colors = plt.cm.Oranges(np.linspace(0.3, 0.9, 12))

# Top left: CellAI OD (matched window)
ax = axes[0, 0]
ci = 0
for label, wells in cellai_conds:
    for w in wells:
        if w in df_ce_m.columns:
            ax.plot(df_ce_m["hours"] * 60, df_ce_m[w], "o-", color=ce_colors[ci],
                    markersize=5, label=f"{label} ({w})")
            ci += 1
ax.set_xlabel("Elapsed Time (min)")
ax.set_ylabel("OD600")
ax.set_title(f"CellAI — OD600 (first {match_min:.0f} min)")
ax.legend(fontsize=6.5, loc="upper left")
ax.grid(True, alpha=0.3)

# Top right: ViNatX OD (matched window)
ax = axes[0, 1]
vi = 0
for label, wells in vinatx_conds:
    for w in wells:
        if w in df_ve_m.columns:
            ax.plot(df_ve_m["hours"] * 60, df_ve_m[w], "s-", color=ve_colors[vi],
                    markersize=5, label=f"{label} ({w})")
            vi += 1
ax.set_xlabel("Elapsed Time (min)")
ax.set_ylabel("OD600")
ax.set_title(f"ViNatX — OD600 ({match_min:.0f} min)")
ax.legend(fontsize=6.5, loc="upper left")
ax.grid(True, alpha=0.3)

# Bottom left: CellAI ln(OD) with fits
ax = axes[1, 0]
ci = 0
for label, wells in cellai_conds:
    for w in wells:
        if w not in df_ce_m.columns:
            continue
        od = df_ce_m[w].values
        t = df_ce_m["hours"].values
        mask = (od > 0) & np.isfinite(od)
        if mask.sum() < 3:
            ci += 1
            continue
        mu, inter, r2 = fit_mu(t, od)
        ax.scatter(t[mask] * 60, np.log(od[mask]), color=ce_colors[ci], s=30, zorder=3)
        t_fit = np.linspace(t[mask].min(), t[mask].max(), 50)
        ax.plot(t_fit * 60, mu * t_fit + inter, "-", color=ce_colors[ci], linewidth=2,
                label=f"{label} ({w}): mu={mu:.3f}")
        ci += 1
ax.set_xlabel("Elapsed Time (min)")
ax.set_ylabel("ln(OD600)")
ax.set_title(f"CellAI — ln(OD) Fits (first {match_min:.0f} min)")
ax.legend(fontsize=6.5, loc="lower right")
ax.grid(True, alpha=0.3)

# Bottom right: ViNatX ln(OD) with fits
ax = axes[1, 1]
vi = 0
for label, wells in vinatx_conds:
    for w in wells:
        if w not in df_ve_m.columns:
            continue
        od = df_ve_m[w].values
        t = df_ve_m["hours"].values
        mask = (od > 0) & np.isfinite(od)
        if mask.sum() < 3:
            vi += 1
            continue
        mu, inter, r2 = fit_mu(t, od)
        ax.scatter(t[mask] * 60, np.log(od[mask]), color=ve_colors[vi], s=30, zorder=3)
        t_fit = np.linspace(t[mask].min(), t[mask].max(), 50)
        ax.plot(t_fit * 60, mu * t_fit + inter, "-", color=ve_colors[vi], linewidth=2,
                label=f"{label} ({w}): mu={mu:.3f}")
        vi += 1
ax.set_xlabel("Elapsed Time (min)")
ax.set_ylabel("ln(OD600)")
ax.set_title(f"ViNatX — ln(OD) Fits ({match_min:.0f} min)")
ax.legend(fontsize=6.5, loc="lower right")
ax.grid(True, alpha=0.3)

# Match y-axes across teams for visual comparison
for row in [0, 1]:
    ymin = min(axes[row, 0].get_ylim()[0], axes[row, 1].get_ylim()[0])
    ymax = max(axes[row, 0].get_ylim()[1], axes[row, 1].get_ylim()[1])
    axes[row, 0].set_ylim(ymin, ymax)
    axes[row, 1].set_ylim(ymin, ymax)
    # Match x-axes to the larger of the two windows
    xmax = max(ce_match_hrs, ve_match_hrs) * 60 + 2
    axes[row, 0].set_xlim(-1, xmax)
    axes[row, 1].set_xlim(-1, xmax)

plt.suptitle(
    f"Time-Matched Growth Comparison — ViNatX {ve_match_hrs*60:.0f} min, CellAI {df_ce_m['hours'].iloc[-1]*60:.0f} min\n"
    f"CellAI (blue, left) vs ViNatX (orange, right) — same axes scales",
    fontsize=14, fontweight="bold",
)
plt.tight_layout()
plt.savefig(FIG_DIR / "time_matched_comparison.png", dpi=150, bbox_inches="tight")
plt.close()

print(f"\nFigure saved: {FIG_DIR / 'time_matched_comparison.png'}")
