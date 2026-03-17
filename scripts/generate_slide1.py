"""
Generate slide1_seeding_density.png — ViNatX Tutorial Plate analysis.

Left panel:  OD600 growth curves by seeding density (raw data, all reps)
Right panel: Growth rate (mu) bar chart, mean ± SD, n=3

Green annotation is placed OUTSIDE the bar chart axes (below x-axis label).
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures" / "presentation"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────

df = pd.read_csv(DATA_DIR / "vinatx_tutorial_plate.csv", parse_dates=["timestamp"])
df = df.sort_values("timestamp").reset_index(drop=True)
t0 = df["timestamp"].iloc[0]
df["minutes"] = (df["timestamp"] - t0).dt.total_seconds() / 60

# ── Seeding density map ───────────────────────────────────────────────────────
# Rows A-F, columns 2-4 (triplicates)
# Media volume → seeding density
vol_map = {"A": 100, "B": 140, "C": 160, "D": 170, "E": 180, "F": 190}
conditions = {}
for row, vol in vol_map.items():
    seed_pct = int((200 - vol) / 200 * 100)
    label = f"{seed_pct}%"
    conditions[label] = [f"{row}{c}" for c in [2, 3, 4]]

# ── Compute growth rates ──────────────────────────────────────────────────────

def fit_mu(t_min, od_vals):
    mask = np.isfinite(od_vals) & (od_vals > 0)
    t = t_min[mask]
    od = od_vals[mask]
    if len(t) < 3:
        return np.nan
    slope, *_ = linregress(t / 60, np.log(od))  # convert min → hr
    return slope

results = {}
for label, wells in conditions.items():
    mus = []
    for w in wells:
        if w in df.columns:
            mu = fit_mu(df["minutes"].values, df[w].values)
            mus.append(mu)
    results[label] = np.array([m for m in mus if np.isfinite(m)])

labels = list(conditions.keys())
means = np.array([results[l].mean() if len(results[l]) > 0 else np.nan for l in labels])
stds  = np.array([results[l].std()  if len(results[l]) > 1 else 0.0 for l in labels])

best_label = labels[np.nanargmax(means)]
best_mu    = np.nanmax(means)

# ── Figure ────────────────────────────────────────────────────────────────────

# Extra bottom margin for the green annotation outside the bar chart
fig = plt.figure(figsize=(14, 6.5))
fig.suptitle("Slide 1: Seeding Density Optimization", fontsize=13, fontweight="bold", y=0.98)

gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35,
                       left=0.07, right=0.97, top=0.91, bottom=0.22)

ax_curves = fig.add_subplot(gs[0])
ax_bars   = fig.add_subplot(gs[1])

# ── Left: OD vs Time ─────────────────────────────────────────────────────────
colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(conditions)))

for i, (label, wells) in enumerate(conditions.items()):
    for j, w in enumerate(wells):
        if w in df.columns:
            ax_curves.plot(df["minutes"], df[w],
                           color=colors[i],
                           linestyle=["-", "--", ":"][j],
                           alpha=0.75,
                           label=label if j == 0 else None)

ax_curves.set_xlabel("Time (min)")
ax_curves.set_ylabel("OD600")
ax_curves.set_title("Seeding Density Screen — OD vs Time")
ax_curves.legend(fontsize=7, title="Seeding Density", title_fontsize=7, loc="upper left")
ax_curves.grid(True, alpha=0.3)

# ── Right: Growth Rate bar chart ──────────────────────────────────────────────
bar_colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(labels)))
bars = ax_bars.bar(labels, means, yerr=stds, capsize=4,
                   color=bar_colors, edgecolor="gray", linewidth=0.8,
                   error_kw=dict(ecolor="black", elinewidth=1.2))

ax_bars.set_xlabel("Seeding Density")
ax_bars.set_ylabel("Growth Rate µ (1/hr)")
ax_bars.set_title("Growth Rate by Seeding Density\n(mean ± SD, n=3)")
ax_bars.grid(True, axis="y", alpha=0.3)

# Highlight best bar
best_idx = labels.index(best_label)
bars[best_idx].set_edgecolor("#2E7D32")
bars[best_idx].set_linewidth(2.5)

# Green annotation OUTSIDE the axes — placed in figure coordinates below the bar chart
ax_bars_bbox = ax_bars.get_position()

# Two lines below the bar chart axes
annotation = (
    f"Best: {best_mu:.3f} 1/hr  ({best_label} seeding)\n"
    "# Used 10% for all experiments"
)
fig.text(
    ax_bars_bbox.x0 + ax_bars_bbox.width / 2,
    0.13,
    annotation,
    ha="center", va="top",
    fontsize=9, color="#2E7D32", fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.4", facecolor="#E8F5E9",
              edgecolor="#2E7D32", linewidth=1.5),
)

fig.savefig(FIG_DIR / "slide1_seeding_density.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved {FIG_DIR / 'slide1_seeding_density.png'}")
