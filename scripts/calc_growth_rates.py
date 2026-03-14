"""
Calculate specific growth rates (μ) from OD600 time-series data.

Fits ln(OD) ~ μ·t + b for each replicate via linear regression.
Reports μ (per hour), doubling time (min), and mean ± std across replicates.
"""

import numpy as np
import pandas as pd
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


def growth_rate(t, od):
    """Fit ln(OD) vs time; return μ (per min), r², doubling time (min)."""
    od = np.asarray(od, dtype=float)
    mask = od > 0
    if mask.sum() < 2:
        return np.nan, np.nan, np.nan
    t_fit = t[mask]
    ln_od = np.log(od[mask])
    coeffs = np.polyfit(t_fit, ln_od, 1)
    mu = coeffs[0]  # per minute
    # r²
    ln_pred = np.polyval(coeffs, t_fit)
    ss_res = np.sum((ln_od - ln_pred) ** 2)
    ss_tot = np.sum((ln_od - ln_od.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    td = np.log(2) / mu if mu > 0 else np.nan
    return mu, r2, td


def calc_rates(df, conditions, plate_name):
    t = df["minutes"].values
    rows = []
    for label, wells in conditions.items():
        available = [w for w in wells if w in df.columns]
        mus, r2s, tds = [], [], []
        for w in available:
            mu, r2, td = growth_rate(t, df[w].values)
            mus.append(mu)
            r2s.append(r2)
            tds.append(td)
        mu_mean = np.nanmean(mus) * 60   # convert to per hour
        mu_std  = np.nanstd(mus)  * 60
        td_mean = np.nanmean(tds)        # already in minutes
        td_std  = np.nanstd(tds)
        r2_mean = np.nanmean(r2s)
        rows.append({
            "Plate":         plate_name,
            "Condition":     label,
            "μ (1/hr)":      round(mu_mean, 4),
            "μ std":         round(mu_std,  4),
            "Doubling (min)": round(td_mean, 1) if not np.isnan(td_mean) else "—",
            "Doubling std":  round(td_std,  1)  if not np.isnan(td_std)  else "—",
            "R²":            round(r2_mean, 3),
            "n reps":        len(available),
        })
    return rows


cellai_df = load_plate(DATA_DIR / "cellai_data.csv")
vinatx_df = load_plate(DATA_DIR / "vinatx_data.csv")

rows = calc_rates(cellai_df, CELLAI_CONDITIONS, "CellAI")
rows += calc_rates(vinatx_df, VINATX_CONDITIONS, "ViNatX")

results = pd.DataFrame(rows)

print("\n=== Growth Rates ===\n")
print(results.to_string(index=False))

print(f"\nCellAI timespan:  {cellai_df['minutes'].iloc[-1]:.0f} min  ({len(cellai_df)} timepoints)")
print(f"ViNatX timespan:  {vinatx_df['minutes'].iloc[-1]:.0f} min  ({len(vinatx_df)} timepoints)")

results.to_csv(DATA_DIR / "growth_rates.csv", index=False)
print("\nSaved growth_rates.csv")
