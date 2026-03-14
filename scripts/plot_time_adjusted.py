"""
Time-adjusted growth rate plots for CellAI and ViNatX experiment plates.

Both plates are plotted on a common elapsed-time axis (hours from first reading)
so growth kinetics are directly comparable.
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
        return np.nan, np.nan, np.nan, np.array([]), np.array([])
    ln_od = np.log(od_m)
    slope, intercept, r, _, _ = linregress(t_m, ln_od)
    return slope, intercept, r**2, t_m, ln_od


# ── Load ─────────────────────────────────────────────────────────────────────

df_ce = load_plate("cellai_experiment_plate.csv")
df_ve = load_plate("vinatx_experiment_plate.csv")

cellai_conditions = {
    "Semi-Def+Tryp+YE+MOPS": {"wells": ["B2", "B3", "B4"], "group": "Semi-Defined"},
    "HBDef+Tryp+YE+Glut": {"wells": ["B5", "B6", "B7"], "group": "High Buffer Defined"},
    "LBv2+Gluc+KH2PO4+MOPS": {"wells": ["B8", "B9", "B10"], "group": "LBv2-based"},
    "Novel Bio Cyclone": {"wells": ["B11"], "group": "Control"},
}

vinatx_conditions = {
    "Def-Min": {"wells": ["A1"], "group": "Defined Minimal"},
    "Def-Min+Gluc": {"wells": ["A2"], "group": "Defined Minimal"},
    "Def-Min+NaCl": {"wells": ["A3"], "group": "Defined Minimal"},
    "Def-Min+Gluc+NaCl": {"wells": ["A4"], "group": "Defined Minimal"},
    "Def-Min+MOPS": {"wells": ["B1"], "group": "Defined+Buffer"},
    "Def-Min+Gluc+MOPS": {"wells": ["B2"], "group": "Defined+Buffer"},
    "Def-Min+MOPS+NaCl": {"wells": ["B3"], "group": "Defined+Buffer"},
    "Def-Min+Gluc+MOPS+NaCl": {"wells": ["B4"], "group": "Defined+Buffer"},
    "NBxCyclone": {"wells": ["C1"], "group": "Rich Media"},
    "LBv2": {"wells": ["C2"], "group": "Rich Media"},
    "Semi-Defined": {"wells": ["C3"], "group": "Rich Media"},
    "Def-Glycerol": {"wells": ["C4"], "group": "Rich Media"},
}

# ── Figure 1: OD600 growth curves, time-adjusted, side by side ───────────────

fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

# Color palettes per group
group_colors_ce = {
    "Semi-Defined": "#1f77b4",
    "High Buffer Defined": "#d62728",
    "LBv2-based": "#9467bd",
    "Control": "#2ca02c",
}
group_colors_ve = {
    "Defined Minimal": "#ff7f0e",
    "Defined+Buffer": "#e377c2",
    "Rich Media": "#8c564b",
}

# CellAI
ax = axes[0]
for label, info in cellai_conditions.items():
    c = group_colors_ce[info["group"]]
    for i, w in enumerate(info["wells"]):
        if w in df_ce.columns:
            ax.plot(df_ce["hours"], df_ce[w], "o-", color=c, markersize=4, alpha=0.8,
                    label=f"{label} ({w})" if True else None)
ax.set_xlabel("Elapsed Time (hours)", fontsize=12)
ax.set_ylabel("OD600", fontsize=12)
ax.set_title("CellAI Experiment", fontsize=14, fontweight="bold")
ax.legend(fontsize=7, loc="upper left")
ax.grid(True, alpha=0.3)

# ViNatX
ax = axes[1]
for label, info in vinatx_conditions.items():
    c = group_colors_ve[info["group"]]
    for w in info["wells"]:
        if w in df_ve.columns:
            ax.plot(df_ve["hours"], df_ve[w], "o-", color=c, markersize=4, alpha=0.8,
                    label=f"{label} ({w})")
ax.set_xlabel("Elapsed Time (hours)", fontsize=12)
ax.set_title("ViNatX Experiment", fontsize=14, fontweight="bold")
ax.legend(fontsize=7, loc="upper left")
ax.grid(True, alpha=0.3)

plt.suptitle("OD600 Growth Curves — Time-Adjusted", fontsize=15, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / "time_adjusted_od.png", dpi=150, bbox_inches="tight")
plt.close()

# ── Figure 2: ln(OD) with regression fits, time-adjusted ────────────────────

fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

# CellAI ln(OD)
ax = axes[0]
mu_list_ce = []
for label, info in cellai_conditions.items():
    c = group_colors_ce[info["group"]]
    for w in info["wells"]:
        if w not in df_ce.columns:
            continue
        mu, intercept, r2, t_m, ln_od = fit_mu(df_ce["hours"], df_ce[w])
        if np.isnan(mu):
            continue
        mu_list_ce.append((label, w, mu, r2))
        ax.scatter(t_m, ln_od, color=c, s=25, alpha=0.7, zorder=3)
        t_fit = np.linspace(t_m.min(), t_m.max(), 50)
        ax.plot(t_fit, mu * t_fit + intercept, "-", color=c, alpha=0.7, linewidth=1.5,
                label=f"{label} ({w}): mu={mu:.3f}")

ax.set_xlabel("Elapsed Time (hours)", fontsize=12)
ax.set_ylabel("ln(OD600)", fontsize=12)
ax.set_title("CellAI — ln(OD) Fits", fontsize=14, fontweight="bold")
ax.legend(fontsize=6.5, loc="lower right")
ax.grid(True, alpha=0.3)

# ViNatX ln(OD)
ax = axes[1]
mu_list_ve = []
for label, info in vinatx_conditions.items():
    c = group_colors_ve[info["group"]]
    for w in info["wells"]:
        if w not in df_ve.columns:
            continue
        mu, intercept, r2, t_m, ln_od = fit_mu(df_ve["hours"], df_ve[w])
        if np.isnan(mu):
            continue
        mu_list_ve.append((label, w, mu, r2))
        ax.scatter(t_m, ln_od, color=c, s=25, alpha=0.7, zorder=3)
        t_fit = np.linspace(t_m.min(), t_m.max(), 50)
        ax.plot(t_fit, mu * t_fit + intercept, "-", color=c, alpha=0.7, linewidth=1.5,
                label=f"{label} ({w}): mu={mu:.3f}")

ax.set_xlabel("Elapsed Time (hours)", fontsize=12)
ax.set_title("ViNatX — ln(OD) Fits", fontsize=14, fontweight="bold")
ax.legend(fontsize=6.5, loc="lower right")
ax.grid(True, alpha=0.3)

plt.suptitle("Growth Rate Fits (ln(OD) vs Time) — Time-Adjusted", fontsize=15, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / "time_adjusted_lnod_fits.png", dpi=150, bbox_inches="tight")
plt.close()

# ── Figure 3: Combined overlay on same axes ──────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Combined OD
ax = axes[0]
for label, info in cellai_conditions.items():
    c = group_colors_ce[info["group"]]
    for w in info["wells"]:
        if w in df_ce.columns:
            ax.plot(df_ce["hours"], df_ce[w], "o-", color=c, markersize=4, alpha=0.7,
                    label=f"[CE] {label} ({w})")
for label, info in vinatx_conditions.items():
    c = group_colors_ve[info["group"]]
    for w in info["wells"]:
        if w in df_ve.columns:
            ax.plot(df_ve["hours"], df_ve[w], "s--", color=c, markersize=4, alpha=0.7,
                    label=f"[VN] {label} ({w})")
ax.set_xlabel("Elapsed Time (hours)", fontsize=12)
ax.set_ylabel("OD600", fontsize=12)
ax.set_title("All Conditions — OD600 Overlay", fontsize=13, fontweight="bold")
ax.legend(fontsize=5, loc="upper left", ncol=2)
ax.grid(True, alpha=0.3)

# Combined ln(OD)
ax = axes[1]
for label, info in cellai_conditions.items():
    c = group_colors_ce[info["group"]]
    for w in info["wells"]:
        if w not in df_ce.columns:
            continue
        od = df_ce[w].values
        t = df_ce["hours"].values
        mask = (od > 0) & np.isfinite(od)
        if mask.sum() < 2:
            continue
        ax.plot(t[mask], np.log(od[mask]), "o-", color=c, markersize=4, alpha=0.7,
                label=f"[CE] {label} ({w})")
for label, info in vinatx_conditions.items():
    c = group_colors_ve[info["group"]]
    for w in info["wells"]:
        if w not in df_ve.columns:
            continue
        od = df_ve[w].values
        t = df_ve["hours"].values
        mask = (od > 0) & np.isfinite(od)
        if mask.sum() < 2:
            continue
        ax.plot(t[mask], np.log(od[mask]), "s--", color=c, markersize=4, alpha=0.7,
                label=f"[VN] {label} ({w})")
ax.set_xlabel("Elapsed Time (hours)", fontsize=12)
ax.set_ylabel("ln(OD600)", fontsize=12)
ax.set_title("All Conditions — ln(OD) Overlay", fontsize=13, fontweight="bold")
ax.legend(fontsize=5, loc="lower right", ncol=2)
ax.grid(True, alpha=0.3)

plt.suptitle("CellAI (circles, solid) vs ViNatX (squares, dashed) — Time-Adjusted",
             fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / "time_adjusted_combined.png", dpi=150, bbox_inches="tight")
plt.close()

# ── Print summary ────────────────────────────────────────────────────────────
print("CellAI growth rates (full time):")
for label, w, mu, r2 in sorted(mu_list_ce, key=lambda x: -x[2]):
    print(f"  {label:>30s} ({w}): mu = {mu:.4f} 1/hr  R2 = {r2:.3f}")

print("\nViNatX growth rates (full time):")
for label, w, mu, r2 in sorted(mu_list_ve, key=lambda x: -x[2]):
    print(f"  {label:>30s} ({w}): mu = {mu:.4f} 1/hr  R2 = {r2:.3f}")

print(f"\nFigures saved:")
print(f"  {FIG_DIR / 'time_adjusted_od.png'}")
print(f"  {FIG_DIR / 'time_adjusted_lnod_fits.png'}")
print(f"  {FIG_DIR / 'time_adjusted_combined.png'}")
