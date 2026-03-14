"""
Plot ln(OD) vs time with linear regression fits used to calculate growth rates.
One subplot per plate, one line+fit per condition.
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


def plot_fits(ax, df, conditions, title, cmap_name="viridis"):
    labels = list(conditions.keys())
    colors = plt.colormaps[cmap_name](np.linspace(0.1, 0.9, len(labels)))
    t = df["minutes"].values
    t_fit = np.linspace(t.min(), t.max(), 100)

    for (label, wells), color in zip(conditions.items(), colors):
        available = [w for w in wells if w in df.columns]
        all_ln_od = []

        for i, w in enumerate(available):
            od = df[w].values.astype(float)
            mask = od > 0
            ln_od = np.where(mask, np.log(od), np.nan)
            ax.plot(t, ln_od, "o", markersize=4, color=color, alpha=0.5,
                    label=label if i == 0 else "_nolegend_")
            all_ln_od.append(ln_od[mask])

        # Fit on all valid points across replicates
        t_valid = np.concatenate([t[df[w].values.astype(float) > 0] for w in available])
        od_valid = np.concatenate([df[w].values.astype(float)[df[w].values.astype(float) > 0]
                                   for w in available])
        if len(t_valid) >= 2:
            ln_valid = np.log(od_valid)
            coeffs = np.polyfit(t_valid, ln_valid, 1)
            mu_hr = coeffs[0] * 60
            ax.plot(t_fit, np.polyval(coeffs, t_fit), "-", color=color,
                    lw=2, alpha=0.9, label=f"  fit μ={mu_hr:.2f}/hr")

    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Time (min from first read)", fontsize=11)
    ax.set_ylabel("ln(OD600)", fontsize=11)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.85,
              ncol=2, title="Condition / fit", title_fontsize=8)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.grid(True, alpha=0.3)


cellai_df = load_plate(DATA_DIR / "cellai_data.csv")
vinatx_df = load_plate(DATA_DIR / "vinatx_data.csv")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle("Growth Rate Fits — ln(OD600) vs Time\nDots = replicates, lines = linear regression fit",
             fontsize=12)

plot_fits(ax1, cellai_df, CELLAI_CONDITIONS,
          f"CellAI ({len(cellai_df)} timepoints)", "plasma")
plot_fits(ax2, vinatx_df, VINATX_CONDITIONS,
          f"ViNatX ({len(vinatx_df)} timepoints)", "viridis")

plt.tight_layout()
plt.savefig(FIG_DIR / "growth_rate_fits.png", dpi=150, bbox_inches="tight")
print("Saved growth_rate_fits.png")
