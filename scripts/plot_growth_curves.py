"""
Plot OD600 growth curves — ln(OD) vs time, one line per replicate.
Reads live data from cellai_data.csv and vinatx_data.csv (updated by fetch_data.py).
"""

import io
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

ROOT     = Path(__file__).parent.parent
DATA_DIR = ROOT / "data"
FIG_DIR  = ROOT / "figures"

CELLAI_CONDITIONS = {
    "50%": ["A1",  "B1",  "C1"],
    "45%": ["A2",  "B2",  "C2"],
    "40%": ["A3",  "B3",  "C3"],
    "35%": ["A4",  "B4",  "C4"],
    "30%": ["A5",  "B5",  "C5"],
    "25%": ["A6",  "B6",  "C6"],
    "20%": ["A7",  "B7",  "C7"],
    "15%": ["A8",  "B8",  "C8"],
    "10%": ["A9",  "B9",  "C9"],
     "5%": ["A10", "B10", "C10"],
}

VINATX_CONDITIONS = {
    "50%": ["A2", "A3", "A4"],
    "30%": ["B2", "B3", "B4"],
    "20%": ["C2", "C3", "C4"],
    "15%": ["D2", "D3", "D4"],
    "10%": ["E2", "E3", "E4"],
     "5%": ["F2", "F3", "F4"],
}


def load_plate(csv_path, conditions):
    df = pd.read_csv(csv_path)
    df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True)
    df = df.sort_values("timestamp").reset_index(drop=True)
    t0 = df["timestamp"].iloc[0]
    df["minutes"] = (df["timestamp"] - t0).dt.total_seconds() / 60

    result = {}
    for label, wells in conditions.items():
        available = [w for w in wells if w in df.columns]
        if available:
            result[label] = {
                "minutes": df["minutes"].values,
                "reps": [df[w].values.astype(float) for w in available],
            }
    return result, df["timestamp"].iloc[-1]


def plot_plate(ax, stats, title, cmap_name="viridis"):
    labels = list(stats.keys())
    colors = plt.colormaps[cmap_name](np.linspace(0.1, 0.9, len(labels)))
    for (label, s), color in zip(stats.items(), colors):
        t = s["minutes"]
        for i, rep in enumerate(s["reps"]):
            ax.plot(t, rep, marker="o", markersize=4, color=color,
                    label=label if i == 0 else "_nolegend_", lw=1.5, alpha=0.8)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Time (min from first read)", fontsize=11)
    ax.set_ylabel("OD600", fontsize=11)
    ax.legend(title="Cell %", fontsize=9, title_fontsize=9,
              loc="upper left", framealpha=0.8)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.grid(True, alpha=0.3)


cellai_path = DATA_DIR / "cellai_data.csv"
vinatx_path = DATA_DIR / "vinatx_data.csv"

if not cellai_path.exists() or not vinatx_path.exists():
    print("Data files not found. Run fetch_data.py first.")
    sys.exit(1)

cellai_stats, cellai_last = load_plate(cellai_path, CELLAI_CONDITIONS)
vinatx_stats, vinatx_last = load_plate(vinatx_path, VINATX_CONDITIONS)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(
    f"Seeding Density Growth Curves — March 14, 2026\n"
    f"Last read: CellAI {cellai_last.strftime('%H:%M PT')} · ViNatX {vinatx_last.strftime('%H:%M PT')}",
    fontsize=12,
)

n_cellai = len(next(iter(cellai_stats.values()))["minutes"])
n_vinatx = len(next(iter(vinatx_stats.values()))["minutes"])
plot_plate(ax1, cellai_stats, f"CellAI (10 conditions, {n_cellai} timepoints)", "plasma")
plot_plate(ax2, vinatx_stats, f"ViNatX (6 conditions, {n_vinatx} timepoints)", "viridis")

plt.tight_layout()
plt.savefig(FIG_DIR / "growth_curves.png", dpi=150, bbox_inches="tight")
print(f"Saved growth_curves.png  (CellAI: {n_cellai} pts, ViNatX: {n_vinatx} pts)")
