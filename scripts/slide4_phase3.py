"""
Phase 3 presentation slide: R1 vs R2 comparison + R3 current experiment.

Layout (2 rows):
  Top:    R1+R2 OD vs time (full data, 500+ min) | Growth rate bars | AUC bars
  Bottom: R3 OD vs time (all 14 datapoints)       | R3 growth rate bars
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures" / "presentation"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Well definitions ─────────────────────────────────────────────────────────

R1_CONDITIONS = {
    "Def-Min": "A1",
    "DM+Gluc": "A2",
    "DM+NaCl": "A3",
    "DM+Gluc+NaCl": "A4",
    "DM+MOPS": "B1",
    "DM+Gluc+MOPS": "B2",
    "DM+MOPS+NaCl": "B3",
    "DM+Gluc+MOPS+NaCl": "B4",
    "NBxCyclone": "C1",
    "LBv2": "C2",
    "Semi-Def": "C3",
    "Def-Glycerol": "C4",
}

R2_CONDITIONS = {
    "DM+MOPS+Glut(100mM)": "D1",
    "DM+MOPS+Tryp+YE": "D2",
    "DM+MOPS+Glut+Tryp+YE": "D3",
    "LBv2+MOPS+Glut(100mM)": "D4",
    "LBv2+Glut(100mM)": "E1",
    "HBDef+Glut+Tryp+YE": "E2",
    "LBv2+MOPS": "E3",
    "LBv2(ctrl)": "E4",
    "DM+MOPS+Tryp+Glut(DOE)": "F1",
}

R3_CONDITIONS = {
    "LBv2 ctrl": "F2",
    "LBv2+Glut100": "F3",
    "LBv2+Glut50+Gluc": "F4",
    "LBv2+Glut10": "G1",
    "LBv2+Glut25": "G2",
    "LBv2+Glut50": "G3",
    "LBv2+Glut75": "G4",
    "LBv2+Tryp+Gluc": "G5",
    "LBv2+Glut100+Gluc": "H1",
    "LBv2+Glut50+Tryp": "H2",
    "LBv2+Glut50+YE": "H3",
    "LBv2+Gluc": "H4",
}


def fit_growth_rate(times_hrs, od_values):
    """Fit ln(OD) vs time (hours) and return (mu, R2)."""
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    ln_od = np.log(od)
    slope, intercept, r, p, se = linregress(t, ln_od)
    return slope, r**2


def fit_exponential_phase(times_hrs, od_values, window_hrs=2.0):
    """Fit growth rate using only the exponential phase (first window_hrs from start of data)."""
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan

    # Use data within the first window_hrs
    t_start = t[0]
    exp_mask = t <= (t_start + window_hrs)
    t_exp = t[exp_mask]
    od_exp = od[exp_mask]

    if len(t_exp) < 3:
        return np.nan, np.nan

    ln_od = np.log(od_exp)
    slope, intercept, r, p, se = linregress(t_exp, ln_od)
    return slope, r**2


def compute_auc(times_hrs, od_values, max_hrs=None):
    """Compute area under OD curve (OD*hr) using trapezoidal rule."""
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if max_hrs is not None:
        time_mask = t <= max_hrs
        t = t[time_mask]
        od = od[time_mask]
    if len(t) < 2:
        return np.nan
    return trapezoid(od, t)


# ── Load data ────────────────────────────────────────────────────────────────

df = pd.read_csv(DATA_DIR / "vinatx_experiment_plate.csv", parse_dates=["timestamp"])
df = df.sort_values("timestamp").reset_index(drop=True)
t0_plate = df["timestamp"].iloc[0]
df["hours"] = (df["timestamp"] - t0_plate).dt.total_seconds() / 3600
df["minutes"] = df["hours"] * 60

# Determine round start times based on first valid data
r1_wells = list(R1_CONDITIONS.values())
r2_wells = list(R2_CONDITIONS.values())
r3_wells = list(R3_CONDITIONS.values())

# R1 starts at first timepoint (all R1 wells have data from start)
r1_start_hr = df["hours"].iloc[0]

# R2 starts when D1 first has data
r2_mask = df["D1"].notna() & (df["D1"] > 0)
r2_start_hr = df.loc[r2_mask, "hours"].iloc[0] if r2_mask.any() else np.nan

# R3 starts when F2 first has data
r3_mask = df["F2"].notna() & (df["F2"] > 0)
r3_start_hr = df.loc[r3_mask, "hours"].iloc[0] if r3_mask.any() else np.nan

print(f"Plate start: {t0_plate}")
print(f"R1 start: {r1_start_hr:.2f} hrs (minute 0)")
print(f"R2 start: {r2_start_hr:.2f} hrs ({(r2_start_hr - r1_start_hr)*60:.0f} min after R1)")
print(f"R3 start: {r3_start_hr:.2f} hrs ({(r3_start_hr - r1_start_hr)*60:.0f} min after R1)")

# Count datapoints per round
r1_pts = df[r1_wells[0]].notna().sum()
r2_pts = df.loc[r2_mask, "D1"].shape[0] if r2_mask.any() else 0
r3_pts = df.loc[r3_mask, "F2"].shape[0] if r3_mask.any() else 0
print(f"\nDatapoints: R1={r1_pts}, R2={r2_pts}, R3={r3_pts}")

# ── Compute growth rates and AUC ─────────────────────────────────────────────

# For R1: exponential phase = first 2 hrs from R1 start
# For R2: exponential phase = first 2 hrs from R2 start
# For R3: use all available data (mostly exponential still)

# Use R2's total duration as the matched AUC window for fair comparison
r2_last_hr = df.loc[r2_mask, "hours"].iloc[-1] if r2_mask.any() else np.nan
auc_window_hrs = r2_last_hr - r2_start_hr  # ~11 hrs
print(f"AUC matched window: {auc_window_hrs:.1f} hrs ({auc_window_hrs*60:.0f} min) from each round's start")

r1_results = {}
for label, well in R1_CONDITIONS.items():
    mask = df[well].notna() & (df[well] > 0)
    t_rel = (df.loc[mask, "hours"] - r1_start_hr).values
    od = df.loc[mask, well].values
    mu, r2 = fit_exponential_phase(t_rel, od, window_hrs=2.0)
    auc = compute_auc(t_rel, od, max_hrs=auc_window_hrs)
    r1_results[label] = {"well": well, "mu": mu, "R2": r2, "auc": auc}

r2_results = {}
for label, well in R2_CONDITIONS.items():
    mask = df[well].notna() & (df[well] > 0)
    t_rel = (df.loc[mask, "hours"] - r2_start_hr).values
    od = df.loc[mask, well].values
    mu, r2 = fit_exponential_phase(t_rel, od, window_hrs=2.0)
    auc = compute_auc(t_rel, od, max_hrs=auc_window_hrs)
    r2_results[label] = {"well": well, "mu": mu, "R2": r2, "auc": auc}

r3_results = {}
for label, well in R3_CONDITIONS.items():
    mask = df[well].notna() & (df[well] > 0)
    if not mask.any():
        continue
    t_rel = (df.loc[mask, "hours"] - r3_start_hr).values
    od = df.loc[mask, well].values
    # For R3 use all data (it's still in exponential phase)
    mu, r2 = fit_growth_rate(t_rel, od)
    r3_results[label] = {"well": well, "mu": mu, "R2": r2, "n_pts": len(od)}

# Print R3 results
print(f"\n{'='*60}")
print("ROUND 3 GROWTH RATES (updated with all datapoints)")
print(f"{'='*60}")
for label, data in sorted(r3_results.items(), key=lambda x: -x[1]["mu"] if np.isfinite(x[1]["mu"]) else -999):
    print(f"  {label:>25s} ({data['well']}): mu = {data['mu']:.4f} 1/hr, R² = {data['R2']:.3f}, n = {data['n_pts']} pts")

# ── Generate Figure ──────────────────────────────────────────────────────────

fig = plt.figure(figsize=(20, 14))
gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 1], hspace=0.35, wspace=0.3)

# Color schemes
R1_COLOR = "#E74C3C"  # Red
R2_COLOR = "#3498DB"  # Blue
R3_CMAP = plt.cm.viridis

# ── Top Left: R1 + R2 OD vs Time (500+ min) ─────────────────────────────────

ax1 = fig.add_subplot(gs[0, 0])

# Plot R1 wells (relative to R1 start, in minutes)
for label, well in R1_CONDITIONS.items():
    mask = df[well].notna() & (df[well] > 0)
    t_min = (df.loc[mask, "hours"] - r1_start_hr).values * 60
    od = df.loc[mask, well].values
    ax1.plot(t_min, od, "o-", color=R1_COLOR, markersize=2.5, alpha=0.5, linewidth=1)

# Plot R2 wells (relative to R2 start, in minutes)
for label, well in R2_CONDITIONS.items():
    mask = df[well].notna() & (df[well] > 0)
    t_min = (df.loc[mask, "hours"] - r2_start_hr).values * 60
    od = df.loc[mask, well].values
    ax1.plot(t_min, od, "s-", color=R2_COLOR, markersize=2.5, alpha=0.5, linewidth=1)

# Add legend handles
ax1.plot([], [], "o-", color=R1_COLOR, label=f"R1 ({r1_pts} pts, 12 wells)", markersize=4)
ax1.plot([], [], "s-", color=R2_COLOR, label=f"R2 ({r2_pts} pts, 9 wells)", markersize=4)
ax1.set_xlabel("Time from round start (min)", fontsize=11)
ax1.set_ylabel("OD600", fontsize=11)
ax1.set_title("R1 + R2 OD vs Time (full data)", fontsize=12, fontweight="bold")
ax1.legend(fontsize=9, loc="upper left")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(-10, max(700, (df["hours"].iloc[-1] - r1_start_hr) * 60 + 20))

# ── Top Middle: Growth Rate comparison (R1 vs R2) ───────────────────────────

ax2 = fig.add_subplot(gs[0, 1])

# Combine R1 and R2 and sort by mu
all_mu = []
for label, data in r1_results.items():
    if np.isfinite(data["mu"]):
        all_mu.append({"label": f"[R1] {label}", "mu": data["mu"], "round": "R1"})
for label, data in r2_results.items():
    if np.isfinite(data["mu"]):
        all_mu.append({"label": f"[R2] {label}", "mu": data["mu"], "round": "R2"})

all_mu.sort(key=lambda x: x["mu"])
bar_colors = [R1_COLOR if x["round"] == "R1" else R2_COLOR for x in all_mu]

bars = ax2.barh(
    [x["label"] for x in all_mu],
    [x["mu"] for x in all_mu],
    color=bar_colors,
    edgecolor="gray",
    alpha=0.85,
)
for bar, entry in zip(bars, all_mu):
    ax2.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height() / 2,
             f"{entry['mu']:.3f}", va="center", fontsize=7)

ax2.set_xlabel("mu (1/hr)", fontsize=11)
ax2.set_title("Growth Rate (exp. phase, 2 hr)", fontsize=12, fontweight="bold")
ax2.tick_params(axis="y", labelsize=7)
ax2.grid(True, axis="x", alpha=0.3)

# ── Top Right: AUC comparison (R1 vs R2) ────────────────────────────────────

ax3 = fig.add_subplot(gs[0, 2])

all_auc = []
for label, data in r1_results.items():
    if np.isfinite(data["auc"]):
        all_auc.append({"label": f"[R1] {label}", "auc": data["auc"], "round": "R1"})
for label, data in r2_results.items():
    if np.isfinite(data["auc"]):
        all_auc.append({"label": f"[R2] {label}", "auc": data["auc"], "round": "R2"})

all_auc.sort(key=lambda x: x["auc"])
bar_colors_auc = [R1_COLOR if x["round"] == "R1" else R2_COLOR for x in all_auc]

bars = ax3.barh(
    [x["label"] for x in all_auc],
    [x["auc"] for x in all_auc],
    color=bar_colors_auc,
    edgecolor="gray",
    alpha=0.85,
)
for bar, entry in zip(bars, all_auc):
    ax3.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2,
             f"{entry['auc']:.2f}", va="center", fontsize=7)

ax3.set_xlabel("AUC (OD*hr)", fontsize=11)
ax3.set_title(f"Total Biomass AUC (matched {auc_window_hrs:.0f} hr)", fontsize=12, fontweight="bold")
ax3.tick_params(axis="y", labelsize=7)
ax3.grid(True, axis="x", alpha=0.3)

# ── Bottom Left: R3 OD vs Time (all datapoints) ─────────────────────────────

ax4 = fig.add_subplot(gs[1, 0])

n_r3 = len(R3_CONDITIONS)
colors_r3 = [R3_CMAP(i / max(n_r3 - 1, 1)) for i in range(n_r3)]

for (label, well), color in zip(R3_CONDITIONS.items(), colors_r3):
    mask = df[well].notna() & (df[well] > 0)
    if not mask.any():
        continue
    t_min = (df.loc[mask, "hours"] - r3_start_hr).values * 60
    od = df.loc[mask, well].values
    ax4.plot(t_min, od, "o-", color=color, markersize=4, linewidth=1.5,
             label=f"{label} ({well})")

ax4.set_xlabel("Time from R3 start (min)", fontsize=11)
ax4.set_ylabel("OD600", fontsize=11)
ax4.set_title(f"Round 3: OD vs Time ({r3_pts} datapoints per well)", fontsize=12, fontweight="bold")
ax4.legend(fontsize=6.5, loc="upper left", ncol=2)
ax4.grid(True, alpha=0.3)

# ── Bottom Middle+Right: R3 Growth Rate bar chart ───────────────────────────

ax5 = fig.add_subplot(gs[1, 1:])

r3_sorted = sorted(r3_results.items(), key=lambda x: x[1]["mu"] if np.isfinite(x[1]["mu"]) else -999)
r3_labels = [f"{label} ({data['well']})" for label, data in r3_sorted]
r3_mus = [data["mu"] for _, data in r3_sorted]
r3_r2s = [data["R2"] for _, data in r3_sorted]

# Color by glutamate dose for visual grouping
glut_doses = {"LBv2 ctrl": 0, "LBv2+Glut10": 10, "LBv2+Glut25": 25,
              "LBv2+Glut50": 50, "LBv2+Glut75": 75, "LBv2+Glut100": 100,
              "LBv2+Glut50+Gluc": 50, "LBv2+Glut100+Gluc": 100,
              "LBv2+Glut50+Tryp": 50, "LBv2+Glut50+YE": 50,
              "LBv2+Tryp+Gluc": 0, "LBv2+Gluc": 0}
max_dose = 100
bar_colors_r3 = []
for label, _ in r3_sorted:
    dose = glut_doses.get(label, 0)
    bar_colors_r3.append(plt.cm.YlOrRd(0.2 + 0.7 * dose / max_dose))

bars = ax5.barh(r3_labels, r3_mus, color=bar_colors_r3, edgecolor="gray", alpha=0.85)
for bar, mu, r2 in zip(bars, r3_mus, r3_r2s):
    if np.isfinite(mu):
        ax5.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height() / 2,
                 f"{mu:.3f} (R²={r2:.2f})", va="center", fontsize=8)

ax5.set_xlabel("mu (1/hr)", fontsize=11)
ax5.set_title(f"Round 3: Growth Rates (all {r3_pts} datapoints, {(df.loc[r3_mask, 'hours'].iloc[-1] - r3_start_hr)*60:.0f} min span)",
              fontsize=12, fontweight="bold")
ax5.tick_params(axis="y", labelsize=9)
ax5.grid(True, axis="x", alpha=0.3)

# Add color bar legend for glutamate dose
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=plt.cm.YlOrRd(0.2), label="0 mM Glut"),
    Patch(facecolor=plt.cm.YlOrRd(0.55), label="50 mM Glut"),
    Patch(facecolor=plt.cm.YlOrRd(0.9), label="100 mM Glut"),
]
ax5.legend(handles=legend_elements, loc="lower right", fontsize=8, title="[Glutamate]")

fig.suptitle("Phase 3: R1 vs R2 Discussion + Round 3 Current Experiment",
             fontsize=16, fontweight="bold", y=0.98)

plt.savefig(FIG_DIR / "slide4_phase3.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {FIG_DIR / 'slide4_phase3.png'}")
