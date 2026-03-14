"""
Plot OD600 vs time (linear scale) with linear regression fits.
Two rows: top = raw OD replicates, bottom = linear fits with slope annotations.
"""

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


def load_plate(csv_path):
    df = pd.read_csv(csv_path)
    df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True)
    df = df.sort_values("timestamp").reset_index(drop=True)
    t0 = df["timestamp"].iloc[0]
    df["minutes"] = (df["timestamp"] - t0).dt.total_seconds() / 60
    return df


def plot_od(ax, df, conditions, title, cmap_name):
    """Raw OD600 replicates, one line per well."""
    labels = list(conditions.keys())
    colors = plt.colormaps[cmap_name](np.linspace(0.1, 0.9, len(labels)))
    t = df["minutes"].values
    for (label, wells), color in zip(conditions.items(), colors):
        available = [w for w in wells if w in df.columns]
        for i, w in enumerate(available):
            ax.plot(t, df[w].values.astype(float), "o-", markersize=4,
                    color=color, lw=1.5, alpha=0.8,
                    label=label if i == 0 else "_nolegend_")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Time (min from first read)", fontsize=10)
    ax.set_ylabel("OD600", fontsize=10)
    ax.legend(title="Cell %", fontsize=8, title_fontsize=8,
              loc="upper left", framealpha=0.8, ncol=2)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.grid(True, alpha=0.3)


def plot_od_fits(ax, df, conditions, title, cmap_name):
    """OD600 with linear regression fits and slope labels."""
    labels = list(conditions.keys())
    colors = plt.colormaps[cmap_name](np.linspace(0.1, 0.9, len(labels)))
    t = df["minutes"].values
    t_fit = np.linspace(t.min(), t.max(), 100)

    for (label, wells), color in zip(conditions.items(), colors):
        available = [w for w in wells if w in df.columns]
        # plot dots
        for i, w in enumerate(available):
            od = df[w].values.astype(float)
            ax.plot(t, od, "o", markersize=4, color=color, alpha=0.4,
                    label=label if i == 0 else "_nolegend_")
        # pooled linear fit
        t_all = np.concatenate([t for _ in available])
        od_all = np.concatenate([df[w].values.astype(float) for w in available])
        mask = od_all > 0
        if mask.sum() >= 2:
            coeffs = np.polyfit(t_all[mask], od_all[mask], 1)
            slope_per_hr = coeffs[0] * 60  # OD/hr
            ax.plot(t_fit, np.polyval(coeffs, t_fit), "-", color=color,
                    lw=2, alpha=0.9,
                    label=f"  {slope_per_hr:+.3f} OD/hr")

    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Time (min from first read)", fontsize=10)
    ax.set_ylabel("OD600", fontsize=10)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.85,
              ncol=2, title="Condition  /  slope", title_fontsize=8)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.grid(True, alpha=0.3)


cellai_df = load_plate(DATA_DIR / "cellai_data.csv")
vinatx_df = load_plate(DATA_DIR / "vinatx_data.csv")

fig, axes = plt.subplots(2, 2, figsize=(16, 10))
fig.suptitle("OD600 Growth — Linear Scale\nTop: raw replicates  ·  Bottom: linear regression fits",
             fontsize=13)

plot_od      (axes[0, 0], cellai_df, CELLAI_CONDITIONS,
              f"CellAI — raw OD ({len(cellai_df)} timepoints)", "plasma")
plot_od      (axes[0, 1], vinatx_df, VINATX_CONDITIONS,
              f"ViNatX — raw OD ({len(vinatx_df)} timepoints)", "viridis")
plot_od_fits (axes[1, 0], cellai_df, CELLAI_CONDITIONS,
              "CellAI — linear fits (OD/hr)", "plasma")
plot_od_fits (axes[1, 1], vinatx_df, VINATX_CONDITIONS,
              "ViNatX — linear fits (OD/hr)", "viridis")

plt.tight_layout()
plt.savefig(FIG_DIR / "growth_linear.png", dpi=150, bbox_inches="tight")
print("Saved growth_linear.png")
