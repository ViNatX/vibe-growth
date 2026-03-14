"""
Round 2 + Round 3 analysis pipeline.

1. Fetches latest OD data for ViNatX experiment plate (now including D1-F1)
2. Calculates growth rates for all wells (Round 1 + Round 2)
3. Appends Round 2 growth rates to all_growth_rates.csv
4. Builds updated DOE response surface with base media as a factor
5. Finds DOE optimum and proposes Round 3 conditions

IMPORTANT: Round 3 proposals ONLY use reagents available on the ViNatX Reagent Plate:
  A1: Defined-Minimal Media
  A2: MOPS (400 mM, pH 7.0)
  A3: NaCl (650 mM)
  A4: Glucose (3% w/v)
  B1: NBxCyclone
  B2: LBv2
  B3: Semi-Defined
  B4: Defined-Glycerol
  C1: Glutamate (1 M Na L-Glutamate)
  C2: Tryptone (100 mg/mL)
  C3: Yeast Extract (100 mg/mL)
  C4: High Buffer Defined Media

Do NOT use CellAI-only reagents (KH2PO4, MgSO4, Trace Metals, FeSO4, Na Citrate).
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from scipy.optimize import minimize, LinearConstraint
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)


def fit_growth_rate(times, od_values):
    """Fit ln(OD) vs time (hours) and return (mu, R2)."""
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    ln_od = np.log(od)
    slope, intercept, r, p, se = linregress(t, ln_od)
    return slope, r**2


# ── Round 2 well formulations ────────────────────────────────────────────────
# Maps well -> (base_media, supplement volumes in uL)
# base_media: 0=Def-Min, 1=LBv2, 2=HBDef

ROUND2_WELLS = {
    "D1": {"base": "Def-Min",  "base_vol": 140, "MOPS": 20, "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 20, "NaCl": 0},
    "D2": {"base": "Def-Min",  "base_vol": 120, "MOPS": 20, "Glucose": 0,  "Tryptone": 20, "YE": 20, "Glutamate": 0,  "NaCl": 0},
    "D3": {"base": "Def-Min",  "base_vol": 120, "MOPS": 20, "Glucose": 0,  "Tryptone": 20, "YE": 10, "Glutamate": 10, "NaCl": 0},
    "D4": {"base": "LBv2",     "base_vol": 140, "MOPS": 20, "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 20, "NaCl": 0},
    "E1": {"base": "LBv2",     "base_vol": 160, "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 20, "NaCl": 0},
    "E2": {"base": "HBDef",    "base_vol": 130, "MOPS": 0,  "Glucose": 0,  "Tryptone": 20, "YE": 20, "Glutamate": 10, "NaCl": 0},
    "E3": {"base": "LBv2",     "base_vol": 160, "MOPS": 20, "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    "E4": {"base": "LBv2",     "base_vol": 180, "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    "F1": {"base": "Def-Min",  "base_vol": 134, "MOPS": 10, "Glucose": 0,  "Tryptone": 26, "YE": 0,  "Glutamate": 10, "NaCl": 0},
}

# Round 1 well formulations (for completeness in the expanded DOE)
# CellAI used different bases: B2-B4=Semi-Def, B5-B7=HBDef, B8-B10=LBv2, B11=NBxCyclone
# ViNatX Round 1: all used Def-Min base except C1=NBxCyclone, C2=LBv2, C3=Semi-Def, C4=Def-Glycerol

ROUND1_CELLAI = [
    {"well": "B2",  "base": "Semi-Def", "MOPS": 11, "Glucose": 0,  "Tryptone": 26, "YE": 14, "Glutamate": 0,  "NaCl": 0},
    {"well": "B3",  "base": "Semi-Def", "MOPS": 7,  "Glucose": 0,  "Tryptone": 14, "YE": 20, "Glutamate": 0,  "NaCl": 0},
    {"well": "B4",  "base": "Semi-Def", "MOPS": 9,  "Glucose": 0,  "Tryptone": 20, "YE": 26, "Glutamate": 0,  "NaCl": 0},
    {"well": "B5",  "base": "HBDef",    "MOPS": 0,  "Glucose": 0,  "Tryptone": 26, "YE": 14, "Glutamate": 10, "NaCl": 0},
    {"well": "B6",  "base": "HBDef",    "MOPS": 0,  "Glucose": 0,  "Tryptone": 14, "YE": 20, "Glutamate": 6,  "NaCl": 0},
    {"well": "B7",  "base": "HBDef",    "MOPS": 0,  "Glucose": 0,  "Tryptone": 20, "YE": 26, "Glutamate": 8,  "NaCl": 0},
    {"well": "B8",  "base": "LBv2",     "MOPS": 12, "Glucose": 18, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "B9",  "base": "LBv2",     "MOPS": 16, "Glucose": 12, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "B10", "base": "LBv2",     "MOPS": 18, "Glucose": 8,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "B11", "base": "NBxCyclone","MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
]

ROUND1_VINATX = [
    {"well": "A1",  "base": "Def-Min",    "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "A2",  "base": "Def-Min",    "MOPS": 0,  "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "A3",  "base": "Def-Min",    "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20},
    {"well": "A4",  "base": "Def-Min",    "MOPS": 0,  "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20},
    {"well": "B1",  "base": "Def-Min",    "MOPS": 20, "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "B2v", "base": "Def-Min",    "MOPS": 20, "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "B3v", "base": "Def-Min",    "MOPS": 20, "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20},
    {"well": "B4v", "base": "Def-Min",    "MOPS": 20, "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20},
    {"well": "C1",  "base": "NBxCyclone", "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "C2",  "base": "LBv2",       "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "C3",  "base": "Semi-Def",   "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
    {"well": "C4",  "base": "Def-Glycerol","MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0},
]

# Known Round 1 mu values (109-min matched window)
ROUND1_MU = {
    # CellAI
    "B2": 0.708, "B3": 0.716, "B4": 0.734, "B5": 0.762, "B6": 0.756,
    "B7": 0.750, "B8": 0.610, "B9": 0.740, "B10": 0.691, "B11": 0.641,
    # ViNatX
    "A1": 0.674, "A2": 0.613, "A3": 0.534, "A4": 0.549,
    "B1": 0.684, "B2v": 0.659, "B3v": 0.703, "B4v": 0.660,
    "C1": 0.683, "C2": 0.725, "C3": 0.625, "C4": 0.241,
}


def load_and_compute_round2_mu():
    """Load ViNatX experiment plate and compute mu for Round 2 wells."""
    csv_path = DATA_DIR / "vinatx_experiment_plate.csv"
    df = pd.read_csv(csv_path, parse_dates=["timestamp"])
    df = df.sort_values("timestamp").reset_index(drop=True)

    t0 = df["timestamp"].iloc[0]
    df["hours"] = (df["timestamp"] - t0).dt.total_seconds() / 3600

    round2_wells = [w for w in ROUND2_WELLS if w in df.columns]
    if not round2_wells:
        print("No Round 2 wells (D1-F1) found in CSV yet.")
        return {}

    results = {}
    for well in round2_wells:
        od = df[well].values
        valid = np.isfinite(od) & (od > 0)
        if valid.sum() < 3:
            print(f"  {well}: insufficient data ({valid.sum()} valid points)")
            continue
        mu, r2 = fit_growth_rate(df["hours"].values, od)
        results[well] = {"mu": mu, "R2": r2}
        print(f"  {well}: mu = {mu:.4f} 1/hr, R² = {r2:.3f}")

    return results


def build_expanded_design_matrix(round2_mu, cellai_r2_mu=None):
    """Build design matrix with base media as a factor.

    Inclusion policy:
      - CellAI Round 1: EXCLUDED (used KH2PO4, Trace Metals, etc. we don't have)
      - CellAI Round 2: INCLUDED (if data provided via cellai_r2_mu)
      - ViNatX Round 1: INCLUDED
      - ViNatX Round 2: INCLUDED
    """
    rows = []

    # CellAI Round 1 — EXCLUDED
    # These conditions used reagents not on ViNatX plate (KH2PO4, Trace Metals)
    # Including them would bias the model toward formulations we can't reproduce
    print("  CellAI Round 1: EXCLUDED (uses KH2PO4, Trace Metals not on our plate)")

    # CellAI Round 2 — include if available
    if cellai_r2_mu:
        print(f"  CellAI Round 2: INCLUDED ({len(cellai_r2_mu)} conditions)")
        for well, data in cellai_r2_mu.items():
            rows.append(data)
    else:
        print("  CellAI Round 2: no data yet")

    # ViNatX Round 1 — INCLUDED
    for entry in ROUND1_VINATX:
        well = entry["well"]
        if well in ROUND1_MU:
            row = {**entry, "mu": ROUND1_MU[well], "team": "ViNatX", "round": 1}
            rows.append(row)
    print(f"  ViNatX Round 1: INCLUDED ({len(ROUND1_VINATX)} conditions)")

    # ViNatX Round 2 — INCLUDED
    r2_count = 0
    for well, formulation in ROUND2_WELLS.items():
        if well in round2_mu:
            row = {
                "well": well,
                "base": formulation["base"],
                "MOPS": formulation["MOPS"],
                "Glucose": formulation["Glucose"],
                "Tryptone": formulation["Tryptone"],
                "YE": formulation["YE"],
                "Glutamate": formulation["Glutamate"],
                "NaCl": formulation["NaCl"],
                "mu": round2_mu[well]["mu"],
                "team": "ViNatX",
                "round": 2,
            }
            rows.append(row)
            r2_count += 1
    print(f"  ViNatX Round 2: INCLUDED ({r2_count} conditions)")

    df = pd.DataFrame(rows)
    print(f"\n  Total datapoints in model: {len(df)}")
    return df


def fit_response_surface(df):
    """Fit response surface with base media as one-hot encoded factor."""
    supplement_factors = ["MOPS", "Glucose", "Tryptone", "YE", "Glutamate", "NaCl"]

    # One-hot encode base media (drop Def-Min as reference)
    base_dummies = pd.get_dummies(df["base"], prefix="base", drop_first=False)
    # Keep all bases for interpretability, RidgeCV handles multicollinearity
    base_cols = list(base_dummies.columns)

    X_supps = df[supplement_factors].values
    X_bases = base_dummies.values
    X = np.hstack([X_supps, X_bases])
    y = df["mu"].values

    all_factors = supplement_factors + base_cols

    print(f"\nExpanded design matrix: {X.shape[0]} experiments x {len(all_factors)} factors")
    print(f"Response range: mu = {y.min():.3f} to {y.max():.3f} 1/hr")
    print(f"Base media types: {df['base'].value_counts().to_dict()}")

    # Fit model
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    poly = PolynomialFeatures(degree=2, include_bias=True)
    X_poly = poly.fit_transform(X_scaled)
    feature_names = poly.get_feature_names_out(all_factors)

    model = RidgeCV(alphas=np.logspace(-2, 3, 50), fit_intercept=False)
    model.fit(X_poly, y)
    y_pred = model.predict(X_poly)
    r2 = r2_score(y, y_pred)

    print(f"\nModel R² = {r2:.4f}  (Ridge alpha = {model.alpha_:.4f})")

    # Top coefficients
    coef_df = pd.DataFrame({"term": feature_names, "coef": model.coef_})
    coef_df["abs_coef"] = coef_df["coef"].abs()
    coef_df = coef_df.sort_values("abs_coef", ascending=False)
    print(f"\nTop 15 coefficients:")
    for _, row in coef_df.head(15).iterrows():
        print(f"  {row['term']:>35s}: {row['coef']:+.4f}")

    df["mu_pred"] = y_pred
    df["residual"] = y - y_pred

    return model, scaler, poly, all_factors, supplement_factors, base_cols, r2


def find_optimum(model, scaler, poly, all_factors, supplement_factors, base_cols):
    """Find optimal supplement volumes for each base media type.

    Only proposes reagents available on ViNatX Reagent Plate:
      MOPS (A2), Glucose (A4), Tryptone (C2), YE (C3), Glutamate (C1), NaCl (A3)
    Does NOT use KH2PO4, MgSO4, Trace Metals, FeSO4, Na Citrate (CellAI only).
    """
    # Bounds: MOPS, Glucose, Tryptone, YE, Glutamate, NaCl
    # Extended Tryp/YE bounds to explore beyond Round 1
    supp_bounds = [(0, 20), (0, 20), (0, 40), (0, 40), (0, 20), (0, 20)]
    n_supps = len(supplement_factors)

    # Encode each base media type
    base_types = {}
    for col in base_cols:
        base_name = col.replace("base_", "")
        vec = np.zeros(len(base_cols))
        vec[base_cols.index(col)] = 1.0
        base_types[base_name] = vec

    def predict_mu(supp_volumes, base_vec):
        x = np.concatenate([supp_volumes, base_vec]).reshape(1, -1)
        x_scaled = scaler.transform(x)
        x_poly = poly.transform(x_scaled)
        return model.predict(x_poly)[0]

    results = {}
    total_constraint = LinearConstraint(
        np.concatenate([np.ones(n_supps), np.zeros(len(base_cols))]),
        lb=0, ub=80
    )

    for base_name, base_vec in base_types.items():
        if base_name in ("Def-Glycerol", "NBxCyclone"):
            continue  # Skip non-viable bases

        def neg_mu(x):
            return -predict_mu(x, base_vec)

        best_mu = -np.inf
        best_x = None
        np.random.seed(42)
        for _ in range(500):
            x0 = np.array([np.random.uniform(lo, hi) for lo, hi in supp_bounds])
            if sum(x0) > 80:
                x0 = x0 * 80 / sum(x0)
            result = minimize(neg_mu, x0, bounds=supp_bounds, method="SLSQP")
            if result.success and -result.fun > best_mu:
                best_mu = -result.fun
                best_x = result.x

        results[base_name] = {
            "mu_pred": best_mu,
            "supplements": dict(zip(supplement_factors, best_x)),
        }
        print(f"\n  Optimum on {base_name}: mu = {best_mu:.3f} 1/hr")
        for name, val in zip(supplement_factors, best_x):
            print(f"    {name:>12s}: {val:5.1f} uL")

    return results


def propose_round3(optima, used_wells):
    """Propose Round 3 conditions based on DOE optima and expert judgment.

    CONSTRAINT: Only uses reagents available on ViNatX Reagent Plate.
    Available bases: Def-Min (A1), LBv2 (B2), HBDef (C4), Semi-Def (B3)
    Available supplements: MOPS (A2), Glucose (A4), Tryptone (C2),
                          YE (C3), Glutamate (C1), NaCl (A3)
    """
    VIABLE_BASES = {"Def-Min", "LBv2", "HBDef", "Semi-Def"}

    print(f"\n{'='*60}")
    print("ROUND 3 PROPOSED CONDITIONS")
    print(f"(Only ViNatX Reagent Plate reagents)")
    print(f"{'='*60}")

    # Filter optima to viable bases only
    viable_optima = {k: v for k, v in optima.items() if k in VIABLE_BASES}
    sorted_optima = sorted(viable_optima.items(), key=lambda x: x[1]["mu_pred"], reverse=True)
    best_base, best_opt = sorted_optima[0]

    proposals = []

    # 1. DOE optimum (best base)
    supps = best_opt["supplements"]
    proposals.append({
        "label": f"DOE optimum ({best_base})",
        "base": best_base,
        **{k: round(v, 1) for k, v in supps.items()},
    })

    # 2. Replicate of DOE optimum
    proposals.append({
        "label": f"DOE optimum replicate ({best_base})",
        "base": best_base,
        **{k: round(v, 1) for k, v in supps.items()},
    })

    # 3. DOE optimum on second-best base (deconfound)
    if len(sorted_optima) > 1:
        second_base, second_opt = sorted_optima[1]
        supps2 = second_opt["supplements"]
        proposals.append({
            "label": f"DOE optimum ({second_base})",
            "base": second_base,
            **{k: round(v, 1) for k, v in supps2.items()},
        })

    # 4. Best base + push Tryptone higher (if at bound)
    if supps.get("Tryptone", 0) >= 35:
        proposals.append({
            "label": f"Push Tryptone ({best_base})",
            "base": best_base,
            "MOPS": round(supps.get("MOPS", 0), 1),
            "Glucose": 0,
            "Tryptone": 50,
            "YE": round(supps.get("YE", 0), 1),
            "Glutamate": round(supps.get("Glutamate", 0), 1),
            "NaCl": 0,
        })

    # 5. Best base + push Glutamate higher (if at bound)
    if supps.get("Glutamate", 0) >= 15:
        proposals.append({
            "label": f"Push Glutamate ({best_base})",
            "base": best_base,
            "MOPS": round(supps.get("MOPS", 0), 1),
            "Glucose": 0,
            "Tryptone": round(supps.get("Tryptone", 0), 1),
            "YE": round(supps.get("YE", 0), 1),
            "Glutamate": 30,
            "NaCl": 0,
        })

    # 6. Cross: best supplements on remaining viable bases
    for base_name, opt in sorted_optima[2:]:
        supps_cross = opt["supplements"]
        proposals.append({
            "label": f"DOE optimum ({base_name})",
            "base": base_name,
            **{k: round(v, 1) for k, v in supps_cross.items()},
        })

    # Print proposals
    for i, p in enumerate(proposals, 1):
        label = p.pop("label")
        base = p.pop("base")
        supp_total = sum(v for k, v in p.items() if k in ["MOPS", "Glucose", "Tryptone", "YE", "Glutamate", "NaCl"])
        base_vol = 180 - supp_total
        print(f"\n  {i}. {label}")
        print(f"     Base: {base} ({base_vol:.0f} uL)")
        for k, v in p.items():
            if v > 0:
                print(f"     {k}: {v:.0f} uL")
        p["label"] = label
        p["base"] = base

    return proposals


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("ROUND 2 GROWTH RATE ANALYSIS")
    print("=" * 60)

    # Step 1: Compute Round 2 growth rates
    round2_mu = load_and_compute_round2_mu()

    if not round2_mu:
        print("\nNo Round 2 data available yet. Exiting.")
        print("Re-run this script after OD data for wells D1-F1 appears.")
        exit(0)

    # Step 1b: Check for CellAI Round 2 data
    cellai_r2_mu = None
    cellai_r2_csv = DATA_DIR / "cellai_experiment_plate2.csv"
    if cellai_r2_csv.exists():
        print(f"\nFound CellAI Round 2 data: {cellai_r2_csv}")
        # TODO: parse CellAI R2 plate, compute mu, map formulations
        # cellai_r2_mu = load_and_compute_cellai_r2_mu()
    else:
        print("\nNo CellAI Round 2 data file found (cellai_experiment_plate2.csv)")

    # Step 2: Build expanded design matrix
    # Policy: exclude CellAI R1 (uses reagents we don't have),
    #         include CellAI R2 + ViNatX R1 + ViNatX R2
    print(f"\n{'='*60}")
    print("EXPANDED RESPONSE SURFACE MODEL (with base media)")
    print(f"{'='*60}")
    design_df = build_expanded_design_matrix(round2_mu, cellai_r2_mu)
    design_df.to_csv(DATA_DIR / "doe_design_matrix_v2.csv", index=False)

    # Step 3: Fit response surface
    model, scaler, poly, all_factors, supp_factors, base_cols, r2 = fit_response_surface(design_df)

    # Step 4: Find optima per base media
    print(f"\n{'='*60}")
    print("OPTIMAL CONDITIONS BY BASE MEDIA")
    print(f"{'='*60}")
    optima = find_optimum(model, scaler, poly, all_factors, supp_factors, base_cols)

    # Step 5: Propose Round 3
    used_wells = list("A1 A2 A3 A4 B1 B2 B3 B4 C1 C2 C3 C4 D1 D2 D3 D4 E1 E2 E3 E4 F1".split())
    proposals = propose_round3(optima, used_wells)

    # Step 6: Append Round 2 to all_growth_rates.csv
    print(f"\n{'='*60}")
    print("APPENDING ROUND 2 TO all_growth_rates.csv")
    print(f"{'='*60}")
    existing = pd.read_csv(DATA_DIR / "all_growth_rates.csv")
    new_rows = []
    for well, data in round2_mu.items():
        formulation = ROUND2_WELLS[well]
        base = formulation["base"]
        supps = []
        for s in ["MOPS", "Glucose", "Tryptone", "YE", "Glutamate", "NaCl"]:
            if formulation[s] > 0:
                supps.append(f"{s} {formulation[s]}")
        condition = f"{base}+{'+'.join(supps)}" if supps else base
        new_rows.append({
            "plate": "ViNatX Experiment R2",
            "condition": f"{condition} ({well})",
            "well": well,
            "mu_per_hr": data["mu"],
            "R2": data["R2"],
            "doubling_min": np.log(2) / data["mu"] * 60 if data["mu"] > 0 else np.inf,
        })
    new_df = pd.DataFrame(new_rows)
    combined = pd.concat([existing, new_df], ignore_index=True)
    combined.to_csv(DATA_DIR / "all_growth_rates.csv", index=False)
    print(f"  Added {len(new_rows)} Round 2 entries. Total: {len(combined)} rows.")

    # Step 7: Generate figures
    print(f"\n{'='*60}")
    print("GENERATING FIGURES")
    print(f"{'='*60}")

    # Figure: Model diagnostics
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: Predicted vs Actual
    ax = axes[0]
    colors = {"CellAI": "#2196F3", "ViNatX": "#FF9800"}
    for team in design_df["team"].unique():
        sub = design_df[design_df["team"] == team]
        ax.scatter(sub["mu"], sub["mu_pred"], c=colors.get(team, "gray"),
                   s=80, edgecolors="k", zorder=3, label=team)
        for _, row in sub.iterrows():
            ax.annotate(row["well"], (row["mu"], row["mu_pred"]), fontsize=7,
                        xytext=(4, 4), textcoords="offset points")
    lims = [0.2, 0.9]
    ax.plot(lims, lims, "r--", alpha=0.5)
    ax.set_xlabel("Actual mu (1/hr)")
    ax.set_ylabel("Predicted mu (1/hr)")
    ax.set_title(f"Model Fit (R² = {r2:.3f})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 2: Round comparison
    ax = axes[1]
    for rnd in sorted(design_df["round"].unique()):
        sub = design_df[design_df["round"] == rnd]
        ax.scatter(sub.index, sub["mu"], s=80, edgecolors="k", zorder=3,
                   label=f"Round {rnd}", alpha=0.8)
        for _, row in sub.iterrows():
            ax.annotate(f"{row['well']}\n{row['base'][:4]}", (row.name, row["mu"]),
                        fontsize=6, ha="center", va="bottom")
    ax.set_ylabel("mu (1/hr)")
    ax.set_title("Growth Rates by Round")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle("Round 2 Response Surface Analysis", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "round2_response_surface.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: round2_response_surface.png")

    print(f"\nDone. Design matrix saved to doe_design_matrix_v2.csv")
