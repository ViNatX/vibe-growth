"""
Build the 5-slide presentation deck as a PDF.

Slide 1: Intro (existing PNG — diagram, not data-driven)
Slide 2: Phase 1 — R1 basal media screen (GENERATED)
Slide 3: Phase 2 — R2 supplement optimization (GENERATED)
Slide 4: Phase 3 — R1 vs R2 extended comparison (GENERATED)
Slide 5: R3 results + conclusions + lessons learned (GENERATED)
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures" / "presentation"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Well definitions ─────────────────────────────────────────────────────────

R1_CONDITIONS = {
    "Def-Min": "A1", "DM+Gluc": "A2", "DM+NaCl": "A3", "DM+Gluc+NaCl": "A4",
    "DM+MOPS": "B1", "DM+Gluc+MOPS": "B2", "DM+MOPS+NaCl": "B3", "DM+Gluc+MOPS+NaCl": "B4",
    "NBxCyclone": "C1", "LBv2": "C2", "Semi-Def": "C3", "Def-Glycerol": "C4",
}

R2_CONDITIONS = {
    "DM+MOPS+Glut(100mM)": "D1", "DM+MOPS+Tryp+YE": "D2",
    "DM+MOPS+Glut+Tryp+YE": "D3", "LBv2+MOPS+Glut(100mM)": "D4",
    "LBv2+Glut(100mM)": "E1", "HBDef+Glut+Tryp+YE": "E2",
    "LBv2+MOPS": "E3", "LBv2(ctrl)": "E4", "DM+MOPS+Tryp+Glut(DOE)": "F1",
}

R3_CONDITIONS = {
    "LBv2 ctrl": "F2", "LBv2+Glut100": "F3", "LBv2+Glut50+Gluc": "F4",
    "LBv2+Glut10": "G1", "LBv2+Glut25": "G2", "LBv2+Glut50": "G3",
    "LBv2+Glut75": "G4", "LBv2+Tryp+Gluc": "G5", "LBv2+Glut100+Gluc": "H1",
    "LBv2+Glut50+Tryp": "H2", "LBv2+Glut50+YE": "H3", "LBv2+Gluc": "H4",
}

# R1 color gradient (light to dark based on growth rate rank)
R1_BASE_COLORS = {
    "Def-Min": "#FFCC80", "DM+Gluc": "#FFCC80", "DM+NaCl": "#FFCC80", "DM+Gluc+NaCl": "#FFCC80",
    "DM+MOPS": "#EF5350", "DM+Gluc+MOPS": "#EF5350", "DM+MOPS+NaCl": "#EF5350", "DM+Gluc+MOPS+NaCl": "#EF5350",
    "NBxCyclone": "#C62828", "LBv2": "#C62828", "Semi-Def": "#C62828", "Def-Glycerol": "#C62828",
}

# R2 color gradient (light to dark based on base media)
R2_BASE_COLORS = {
    "DM+MOPS+Glut(100mM)": "#BBDEFB", "DM+MOPS+Tryp+YE": "#BBDEFB",
    "DM+MOPS+Glut+Tryp+YE": "#BBDEFB", "DM+MOPS+Tryp+Glut(DOE)": "#1565C0",
    "LBv2+MOPS+Glut(100mM)": "#42A5F5", "LBv2+Glut(100mM)": "#42A5F5",
    "LBv2+MOPS": "#42A5F5", "LBv2(ctrl)": "#1E88E5",
    "HBDef+Glut+Tryp+YE": "#64B5F6",
}


def fit_exponential_phase(times_hrs, od_values, window_hrs=2.0):
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    t_start = t[0]
    exp_mask = t <= (t_start + window_hrs)
    t_exp, od_exp = t[exp_mask], od[exp_mask]
    if len(t_exp) < 3:
        return np.nan, np.nan
    ln_od = np.log(od_exp)
    slope, intercept, r, p, se = linregress(t_exp, ln_od)
    return slope, r**2


def fit_growth_rate(times_hrs, od_values):
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    ln_od = np.log(od)
    slope, intercept, r, p, se = linregress(t, ln_od)
    return slope, r**2


def interpolate_od(times_hrs, od_values, target_hrs):
    """Interpolate OD at target time points, filling gaps linearly."""
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 2:
        return target_hrs, np.full_like(target_hrs, np.nan)
    od_interp = np.interp(target_hrs, t, od, left=np.nan, right=np.nan)
    return target_hrs, od_interp


AUC_WINDOW_MIN = 600  # matched AUC window in minutes
AUC_WINDOW_HRS = AUC_WINDOW_MIN / 60


def compute_auc_matched(times_hrs, od_values):
    """Compute AUC over a fixed 600-min window, interpolating gaps."""
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 2:
        return np.nan
    # Create a fine grid from 0 to AUC_WINDOW_HRS and interpolate
    t_grid = np.linspace(0, AUC_WINDOW_HRS, 300)
    od_grid = np.interp(t_grid, t, od, left=np.nan, right=np.nan)
    # Only integrate where we have valid interpolated values
    valid = np.isfinite(od_grid)
    if valid.sum() < 2:
        return np.nan
    return trapezoid(od_grid[valid], t_grid[valid])


# ── Load data ────────────────────────────────────────────────────────────────

df = pd.read_csv(DATA_DIR / "vinatx_experiment_plate.csv", parse_dates=["timestamp"])
df = df.sort_values("timestamp").reset_index(drop=True)
t0_plate = df["timestamp"].iloc[0]
df["hours"] = (df["timestamp"] - t0_plate).dt.total_seconds() / 3600

r1_start_hr = df["hours"].iloc[0]
r2_mask = df["D1"].notna() & (df["D1"] > 0)
r2_start_hr = df.loc[r2_mask, "hours"].iloc[0]
r3_mask = df["F2"].notna() & (df["F2"] > 0)
r3_start_hr = df.loc[r3_mask, "hours"].iloc[0]

r2_last_hr = df.loc[r2_mask, "hours"].iloc[-1]
auc_window_hrs = r2_last_hr - r2_start_hr

r1_pts = df["A1"].notna().sum()
r2_pts = r2_mask.sum()
r3_pts = r3_mask.sum()

# ── Compute metrics ──────────────────────────────────────────────────────────

R1_COLOR = "#E74C3C"
R2_COLOR = "#3498DB"

r1_results = {}
for label, well in R1_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0)
    t_rel = (df.loc[m, "hours"] - r1_start_hr).values
    od = df.loc[m, well].values
    mu, r2 = fit_exponential_phase(t_rel, od, window_hrs=2.0)
    auc = compute_auc_matched(t_rel, od)
    r1_results[label] = {"well": well, "mu": mu, "R2": r2, "auc": auc}

r2_results = {}
for label, well in R2_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0)
    t_rel = (df.loc[m, "hours"] - r2_start_hr).values
    od = df.loc[m, well].values
    mu, r2 = fit_exponential_phase(t_rel, od, window_hrs=2.0)
    auc = compute_auc_matched(t_rel, od)
    r2_results[label] = {"well": well, "mu": mu, "R2": r2, "auc": auc}

r3_results = {}
for label, well in R3_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0)
    if not m.any():
        continue
    t_rel = (df.loc[m, "hours"] - r3_start_hr).values
    od = df.loc[m, well].values
    mu, r2 = fit_growth_rate(t_rel, od)
    r3_results[label] = {"well": well, "mu": mu, "R2": r2, "n_pts": len(od)}


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 2: Phase 1 — R1 Basal Media Screen
# ══════════════════════════════════════════════════════════════════════════════

fig2 = plt.figure(figsize=(22, 10))
gs2 = GridSpec(1, 3, figure=fig2, width_ratios=[1.0, 1.2, 0.35], wspace=0.35)

fig2.suptitle("Phase 1: Which existing media formulations perform best?",
              fontsize=18, fontweight="bold", y=0.97)

# Panel 1: R1 growth curves (first 108 min)
ax_s2a = fig2.add_subplot(gs2[0, 0])

# Use R1 data only up to ~108 min from R1 start
r1_108_mask = (df["hours"] - r1_start_hr) <= 1.85  # ~111 min
group_colors = {
    "Defined Minimal": "#FFCC80", "Defined+Buffer": "#EF5350", "Rich Media": "#C62828",
}
group_map = {
    "Def-Min": "Defined Minimal", "DM+Gluc": "Defined Minimal",
    "DM+NaCl": "Defined Minimal", "DM+Gluc+NaCl": "Defined Minimal",
    "DM+MOPS": "Defined+Buffer", "DM+Gluc+MOPS": "Defined+Buffer",
    "DM+MOPS+NaCl": "Defined+Buffer", "DM+Gluc+MOPS+NaCl": "Defined+Buffer",
    "NBxCyclone": "Rich Media", "LBv2": "Rich Media",
    "Semi-Def": "Rich Media", "Def-Glycerol": "Rich Media",
}

for label, well in R1_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0) & r1_108_mask
    t_min = (df.loc[m, "hours"] - r1_start_hr).values * 60
    od = df.loc[m, well].values
    color = group_colors[group_map[label]]
    ax_s2a.plot(t_min, od, "o-", color=color, markersize=4, linewidth=1.5, alpha=0.8,
                label=f"{R1_CONDITIONS[label]}: {label}")

ax_s2a.set_xlabel("Time (min)", fontsize=12)
ax_s2a.set_ylabel("OD600", fontsize=12)
ax_s2a.set_title("Round 1: Basal Media Screen\n12 conditions, 108 min", fontsize=13, fontweight="bold")
ax_s2a.legend(fontsize=7, loc="upper left", ncol=2)
ax_s2a.grid(True, alpha=0.3)

# Panel 2: Growth Rate ranking bar chart
ax_s2b = fig2.add_subplot(gs2[0, 1])

r1_sorted = sorted(r1_results.items(), key=lambda x: x[1]["mu"] if np.isfinite(x[1]["mu"]) else -999)
r1_bar_labels = [f"{data['well']}: {label}" for label, data in r1_sorted]
r1_bar_mus = [data["mu"] for _, data in r1_sorted]
r1_bar_colors = [group_colors[group_map[label]] for label, _ in r1_sorted]

bars = ax_s2b.barh(r1_bar_labels, r1_bar_mus, color=r1_bar_colors, edgecolor="gray", alpha=0.85)
for bar, mu in zip(bars, r1_bar_mus):
    if np.isfinite(mu):
        ax_s2b.text(bar.get_width() - 0.01, bar.get_y() + bar.get_height() / 2,
                    f"{mu:.3f}", va="center", ha="right", fontsize=9,
                    color="white", fontweight="bold")

ax_s2b.set_xlabel("Growth Rate mu (1/hr)", fontsize=12)
ax_s2b.set_title("Growth Rate Ranking", fontsize=13, fontweight="bold")
ax_s2b.set_xlim(0, 0.85)
ax_s2b.tick_params(axis="y", labelsize=9)
ax_s2b.grid(True, axis="x", alpha=0.3)

# Panel 3: Key insights text
ax_s2c = fig2.add_subplot(gs2[0, 2])
ax_s2c.axis("off")

bbox_green = dict(boxstyle="round,pad=0.4", facecolor="#E8F5E9", edgecolor="#2E7D32", linewidth=2)
ax_s2c.text(0.05, 0.92, "LBv2 >> Defined Media\n(pre-made peptides win)",
            transform=ax_s2c.transAxes, fontsize=10, va="top", fontweight="bold", bbox=bbox_green)

bbox_blue = dict(boxstyle="round,pad=0.4", facecolor="#E3F2FD", edgecolor="#1565C0", linewidth=2)
ax_s2c.text(0.05, 0.68, "MOPS buffer helps\n(+0.05 on DM base)",
            transform=ax_s2c.transAxes, fontsize=10, va="top", fontweight="bold", bbox=bbox_blue)

bbox_red = dict(boxstyle="round,pad=0.4", facecolor="#FFEBEE", edgecolor="#C62828", linewidth=2)
ax_s2c.text(0.05, 0.44, "Glucose hurts\nNaCl hurts",
            transform=ax_s2c.transAxes, fontsize=10, va="top", fontweight="bold", bbox=bbox_red)

fig2.tight_layout(rect=[0, 0, 1, 0.94])
fig2.savefig(FIG_DIR / "slide2_phase1.png", dpi=150, bbox_inches="tight")
plt.close(fig2)
print("Saved slide2_phase1.png")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 3: Phase 2 — R2 Supplement Optimization
# ══════════════════════════════════════════════════════════════════════════════

fig3 = plt.figure(figsize=(22, 10))
gs3 = GridSpec(1, 3, figure=fig3, width_ratios=[1.0, 1.0, 0.55], wspace=0.35)

fig3.suptitle("Phase 2: Supplement optimization — what we learned at 2 hours + R3 design",
              fontsize=17, fontweight="bold", y=0.97)

# Panel 1: R2 growth curves (first 2 hours)
ax_s3a = fig3.add_subplot(gs3[0, 0])

r2_2hr_mask = (df["hours"] - r2_start_hr).between(0, 2.1) & r2_mask
for label, well in R2_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0) & r2_2hr_mask
    if not m.any():
        continue
    t_min = (df.loc[m, "hours"] - r2_start_hr).values * 60
    od = df.loc[m, well].values
    color = R2_BASE_COLORS.get(label, "#90CAF9")
    ax_s3a.plot(t_min, od, "o-", color=color, markersize=5, linewidth=2, alpha=0.85,
                label=f"{well}: {label}")

ax_s3a.set_xlabel("Time (min)", fontsize=12)
ax_s3a.set_ylabel("OD600", fontsize=12)
ax_s3a.set_title("R2 Growth Curves\n(first 2 hours)", fontsize=13, fontweight="bold")
ax_s3a.legend(fontsize=7, loc="upper left")
ax_s3a.grid(True, alpha=0.3)

# Panel 2: R2 Growth Rate bar chart
ax_s3b = fig3.add_subplot(gs3[0, 1])

r2_sorted = sorted(r2_results.items(), key=lambda x: x[1]["mu"] if np.isfinite(x[1]["mu"]) else -999)
r2_bar_labels = [label for label, _ in r2_sorted]
r2_bar_mus = [data["mu"] for _, data in r2_sorted]
r2_bar_colors = [R2_BASE_COLORS.get(label, "#90CAF9") for label, _ in r2_sorted]

bars = ax_s3b.barh(r2_bar_labels, r2_bar_mus, color=r2_bar_colors, edgecolor="gray", alpha=0.85)
for bar, mu in zip(bars, r2_bar_mus):
    if np.isfinite(mu):
        ax_s3b.text(bar.get_width() - 0.01, bar.get_y() + bar.get_height() / 2,
                    f"{mu:.2f}", va="center", ha="right", fontsize=9,
                    color="white", fontweight="bold")

ax_s3b.set_xlabel("mu (1/hr)", fontsize=12)
ax_s3b.set_title("R2 Growth Rate\n(2-hr snapshot)", fontsize=13, fontweight="bold")
ax_s3b.set_xlim(0, 0.85)
ax_s3b.tick_params(axis="y", labelsize=9)
ax_s3b.grid(True, axis="x", alpha=0.3)

# Panel 3: Conclusions + R3 design
ax_s3c = fig3.add_subplot(gs3[0, 2])
ax_s3c.axis("off")

conclusions_r2 = """CONCLUSIONS FROM EARLY DATA:

1. LBv2 + Glut 100mM = best mu (0.709)

2. MOPS didn't help on LBv2
   (D4=0.708 vs E1=0.709)

3. 100mM Glut on Def-Min killed growth
   (D1=0.366, likely osmotic stress)

4. HBDef+supps (E2=0.661) didn't
   beat bare LBv2"""

r3_design = """ROUND 3 DESIGN (12 wells, 40 transfers):

Glutamate dose-DOWN titration:
  G1: 10mM   G2: 25mM
  G3: 50mM   G4: 75mM

Glucose "supercharger" tests:
  H1: Glut 100mM + Glucose
  H4: Glucose only

Combo tests:
  H2: Glut 50mM + extra Tryp
  H3: Glut 50mM + extra YE

Controls:
  F2: Bare LBv2 (baseline drift)
  F3: Glut 100mM replicate"""

bbox_orange = dict(boxstyle="round,pad=0.4", facecolor="#FFF3E0", edgecolor="#E65100", linewidth=2)
ax_s3c.text(0.02, 0.98, conclusions_r2, transform=ax_s3c.transAxes,
            fontsize=8, va="top", fontfamily="monospace", bbox=bbox_orange)

bbox_green = dict(boxstyle="round,pad=0.4", facecolor="#E8F5E9", edgecolor="#2E7D32", linewidth=2)
ax_s3c.text(0.02, 0.40, r3_design, transform=ax_s3c.transAxes,
            fontsize=7.5, va="top", fontfamily="monospace", bbox=bbox_green)

fig3.tight_layout(rect=[0, 0, 1, 0.94])
fig3.savefig(FIG_DIR / "slide3_phase2.png", dpi=150, bbox_inches="tight")
plt.close(fig3)
print("Saved slide3_phase2.png")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 4: R1 vs R2 Extended Comparison
# ══════════════════════════════════════════════════════════════════════════════

fig4 = plt.figure(figsize=(24, 10))
gs4 = GridSpec(1, 3, figure=fig4, width_ratios=[0.8, 1.1, 1.1], wspace=0.45)

fig4.suptitle("Phase 3: Extended data — R1 vs R2 comparison",
              fontsize=18, fontweight="bold", y=0.97)

# Panel 1: OD vs Time
ax1 = fig4.add_subplot(gs4[0, 0])

for label, well in R1_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0)
    t_min = (df.loc[m, "hours"] - r1_start_hr).values * 60
    od = df.loc[m, well].values
    ax1.plot(t_min, od, "o-", color=R1_COLOR, markersize=2.5, alpha=0.45, linewidth=1)

for label, well in R2_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0)
    t_min = (df.loc[m, "hours"] - r2_start_hr).values * 60
    od = df.loc[m, well].values
    ax1.plot(t_min, od, "s-", color=R2_COLOR, markersize=2.5, alpha=0.45, linewidth=1)

ax1.plot([], [], "o-", color=R1_COLOR, label="R1 (12 wells)", markersize=4)
ax1.plot([], [], "s-", color=R2_COLOR, label="R2 (9 wells)", markersize=4)
ax1.set_xlabel("Time from round start (min)", fontsize=11)
ax1.set_ylabel("OD600", fontsize=11)
ax1.set_title("R1 + R2 OD vs Time (full data)", fontsize=13, fontweight="bold")
ax1.legend(fontsize=10, loc="upper left")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(-10, 1000)

# Panel 2: Growth Rate bars
ax2 = fig4.add_subplot(gs4[0, 1])

all_mu = []
for label, data in r1_results.items():
    if np.isfinite(data["mu"]):
        all_mu.append({"label": f"[R1] {label}", "mu": data["mu"], "round": "R1"})
for label, data in r2_results.items():
    if np.isfinite(data["mu"]):
        all_mu.append({"label": f"[R2] {label}", "mu": data["mu"], "round": "R2"})
all_mu.sort(key=lambda x: x["mu"])

bars = ax2.barh(
    [x["label"] for x in all_mu], [x["mu"] for x in all_mu],
    color=[R1_COLOR if x["round"] == "R1" else R2_COLOR for x in all_mu],
    edgecolor="gray", alpha=0.85,
)
for bar, entry in zip(bars, all_mu):
    ax2.text(bar.get_width() - 0.01, bar.get_y() + bar.get_height() / 2,
             f"{entry['mu']:.3f}", va="center", ha="right", fontsize=7, color="white", fontweight="bold")
ax2.set_xlabel("mu (1/hr)", fontsize=11)
ax2.set_title("Growth Rate (exp. phase, first 2 hr)", fontsize=13, fontweight="bold")
ax2.tick_params(axis="y", labelsize=7.5)
ax2.set_xlim(0, 0.85)
ax2.grid(True, axis="x", alpha=0.3)

# Panel 3: AUC bars
ax3 = fig4.add_subplot(gs4[0, 2])

all_auc = []
for label, data in r1_results.items():
    if np.isfinite(data["auc"]):
        all_auc.append({"label": f"[R1] {label}", "auc": data["auc"], "round": "R1"})
for label, data in r2_results.items():
    if np.isfinite(data["auc"]):
        all_auc.append({"label": f"[R2] {label}", "auc": data["auc"], "round": "R2"})
all_auc.sort(key=lambda x: x["auc"])

bars = ax3.barh(
    [x["label"] for x in all_auc], [x["auc"] for x in all_auc],
    color=[R1_COLOR if x["round"] == "R1" else R2_COLOR for x in all_auc],
    edgecolor="gray", alpha=0.85,
)
for bar, entry in zip(bars, all_auc):
    ax3.text(bar.get_width() - 0.1, bar.get_y() + bar.get_height() / 2,
             f"{entry['auc']:.1f}", va="center", ha="right", fontsize=7, color="white", fontweight="bold")
ax3.set_xlabel("AUC (OD*hr)", fontsize=11)
ax3.set_title(f"Total Biomass AUC (matched {AUC_WINDOW_MIN:.0f} min, interpolated)", fontsize=13, fontweight="bold")
ax3.tick_params(axis="y", labelsize=7.5)
ax3.grid(True, axis="x", alpha=0.3)

fig4.tight_layout(rect=[0, 0, 1, 0.94])
fig4.savefig(FIG_DIR / "slide4_phase3.png", dpi=150, bbox_inches="tight")
plt.close(fig4)
print("Saved slide4_phase3.png")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 5: R3 Results + Conclusions + Next Round
# ══════════════════════════════════════════════════════════════════════════════

fig5 = plt.figure(figsize=(24, 10))
gs5 = GridSpec(1, 3, figure=fig5, width_ratios=[1.0, 1.1, 0.8], wspace=0.35)

fig5.suptitle("Round 3 Results, Conclusions & Lessons Learned",
              fontsize=18, fontweight="bold", y=0.97)

# Panel 1: R3 OD vs Time
ax4 = fig5.add_subplot(gs5[0, 0])

R3_CMAP = plt.cm.viridis
n_r3 = len(R3_CONDITIONS)
colors_r3 = [R3_CMAP(i / max(n_r3 - 1, 1)) for i in range(n_r3)]

for (label, well), color in zip(R3_CONDITIONS.items(), colors_r3):
    m = df[well].notna() & (df[well] > 0)
    if not m.any():
        continue
    t_min = (df.loc[m, "hours"] - r3_start_hr).values * 60
    od = df.loc[m, well].values
    ax4.plot(t_min, od, "o-", color=color, markersize=4, linewidth=1.5,
             label=f"{label} ({well})")

ax4.set_xlabel("Time from R3 start (min)", fontsize=11)
ax4.set_ylabel("OD600", fontsize=11)
r3_span_min = (df.loc[r3_mask, "hours"].iloc[-1] - r3_start_hr) * 60
ax4.set_title(f"R3 Growth Curves ({r3_pts} pts, {r3_span_min:.0f} min)", fontsize=13, fontweight="bold")
ax4.legend(fontsize=6, loc="upper left", ncol=2)
ax4.grid(True, alpha=0.3)

# Panel 2: R3 Growth Rate bars (values inside bars)
ax5 = fig5.add_subplot(gs5[0, 1])

r3_sorted = sorted(r3_results.items(), key=lambda x: x[1]["mu"] if np.isfinite(x[1]["mu"]) else -999)
r3_labels = [f"{label} ({data['well']})" for label, data in r3_sorted]
r3_mus = [data["mu"] for _, data in r3_sorted]
r3_r2s = [data["R2"] for _, data in r3_sorted]

glut_doses = {"LBv2 ctrl": 0, "LBv2+Glut10": 10, "LBv2+Glut25": 25,
              "LBv2+Glut50": 50, "LBv2+Glut75": 75, "LBv2+Glut100": 100,
              "LBv2+Glut50+Gluc": 50, "LBv2+Glut100+Gluc": 100,
              "LBv2+Glut50+Tryp": 50, "LBv2+Glut50+YE": 50,
              "LBv2+Tryp+Gluc": 0, "LBv2+Gluc": 0}
bar_colors_r3 = [plt.cm.YlOrRd(0.2 + 0.7 * glut_doses.get(label, 0) / 100) for label, _ in r3_sorted]

bars = ax5.barh(r3_labels, r3_mus, color=bar_colors_r3, edgecolor="gray", alpha=0.85)
for bar, mu, r2 in zip(bars, r3_mus, r3_r2s):
    if np.isfinite(mu):
        ax5.text(bar.get_width() - 0.008, bar.get_y() + bar.get_height() / 2,
                 f"{mu:.3f} (R\u00b2={r2:.2f})", va="center", ha="right",
                 fontsize=7.5, color="white", fontweight="bold")

ax5.set_xlabel("mu (1/hr)", fontsize=11)
ax5.set_title(f"R3 Growth Rates (all {r3_pts} datapoints)", fontsize=13, fontweight="bold")
ax5.tick_params(axis="y", labelsize=8.5)
ax5.set_xlim(0, 0.80)
ax5.grid(True, axis="x", alpha=0.3)

legend_elements = [
    Patch(facecolor=plt.cm.YlOrRd(0.2), label="0 mM Glut"),
    Patch(facecolor=plt.cm.YlOrRd(0.55), label="50 mM Glut"),
    Patch(facecolor=plt.cm.YlOrRd(0.9), label="100 mM Glut"),
]
ax5.legend(handles=legend_elements, loc="lower right", fontsize=8, title="[Glutamate]")

# Panel 3: Conclusions + Lessons Learned
ax6 = fig5.add_subplot(gs5[0, 2])
ax6.axis("off")

conclusions_text = """KEY FINDINGS

1. Glutamate dose-response is clear
   100 > 75 > 50 > 25 > 10 mM
   Best: LBv2+Glut 100mM (mu=0.688)

2. Glutamate = best N-source
   Feeds directly into TCA cycle
   Outperforms tryptone and YE

3. Glucose HURTS growth rate
   Likely overflow metabolism

4. F3 replicates E1 (R2) well:
   mu 0.688 vs 0.709"""

lessons_text = """WHAT WE LEARNED ABOUT AI + BIO

AI got us to hypotheses fast —
initial DOE was wrong (glucose/NaCl
on Def-Min), but rapid iteration
let us pivot to glutamate + LBv2.

Integration was the bottleneck:
MCP, Notion, Monomer took effort.
Once connected, R2 and R3 moved
much faster than R1.

AI made mistakes a grad student
would also make, but lacked biases
against "unusual" experiments."""

bottom_text = """16% IMPROVEMENT IN < 24 HOURS
LBv2+Glut 100mM (0.688) vs
bare LBv2 (0.592) — not plateaued
Next: titrate Glut UP + add MOPS"""

bbox_orange = dict(boxstyle="round,pad=0.4", facecolor="#FFF3E0", edgecolor="#E65100", linewidth=2)
ax6.text(0.02, 0.98, conclusions_text, transform=ax6.transAxes,
         fontsize=8, va="top", fontfamily="monospace", bbox=bbox_orange)

bbox_blue = dict(boxstyle="round,pad=0.4", facecolor="#E3F2FD", edgecolor="#1565C0", linewidth=2)
ax6.text(0.02, 0.53, lessons_text, transform=ax6.transAxes,
         fontsize=7.5, va="top", fontfamily="monospace", bbox=bbox_blue)

bbox_green = dict(boxstyle="round,pad=0.4", facecolor="#E8F5E9", edgecolor="#2E7D32", linewidth=2)
ax6.text(0.02, 0.13, bottom_text, transform=ax6.transAxes,
         fontsize=8, va="top", fontfamily="monospace", fontweight="bold", bbox=bbox_green)

fig5.tight_layout(rect=[0, 0, 1, 0.94])
fig5.savefig(FIG_DIR / "slide5_conclusions.png", dpi=150, bbox_inches="tight")
plt.close(fig5)
print("Saved slide5_conclusions.png")


# ══════════════════════════════════════════════════════════════════════════════
# BUILD PDF DECK
# ══════════════════════════════════════════════════════════════════════════════

slide_files = [
    FIG_DIR / "slide1_intro.png",
    FIG_DIR / "slide2_phase1.png",
    FIG_DIR / "slide3_phase2.png",
    FIG_DIR / "slide4_phase3.png",
    FIG_DIR / "slide5_conclusions.png",
]

pdf_path = FIG_DIR / "ViNatX_Presentation.pdf"

with PdfPages(pdf_path) as pdf:
    for slide_path in slide_files:
        img = Image.open(slide_path)
        w, h = img.size
        fig_w = 16
        fig_h = fig_w * h / w
        fig = plt.figure(figsize=(fig_w, fig_h))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(img)
        ax.axis("off")
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        print(f"  Added to PDF: {slide_path.name}")

print(f"\nFull deck saved: {pdf_path}")
