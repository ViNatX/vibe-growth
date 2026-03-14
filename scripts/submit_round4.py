"""
Round 4 protocol builder.

Design rationale:
- R3 showed clear glutamate dose-response (100 > 75 > 50 > 25 > 10 mM)
  with no plateau at 100 uL stock volume
- R1 showed MOPS-buffered wells had highest AUC over extended time
- R4 hypothesis: push glutamate higher + cross with MOPS for
  fast growth rate (glutamate) + sustained yield (MOPS buffering)

Layout: 12 wells, exactly 40 transfers (at the limit)

Usage:
  pixi run python scripts/submit_round4.py
"""

import json
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"

# Plate barcodes (Monomer Bio workcell)
EXPERIMENT_PLATE = "ViNatX Experiment Plate"
CELL_STOCK_PLATE = "ViNatX Cell Stock Plate"
REAGENT_NAME = "CellAI Reagent Plate"  # has LBv2, MOPS, and Glutamate

# CellAI Reagent Plate well map (verified from Monomer Cloud)
# LBv2 at A3 (2610 uL), MOPS pH7 at B3 (2860 uL), Na L-Glutamate at C3 (2756 uL)
BASE_TO_REAGENT = {
    "LBv2": "A3",
}

SUPP_TO_REAGENT = {
    "MOPS": "B3",
    "Glutamate": "C3",
}

# ── R4 Conditions ────────────────────────────────────────────────────────────
# All volumes in uL. Total media = 180 uL, + 20 uL inoculum = 200 uL/well.
# Glutamate stock at C1, MOPS (400 mM) at A2, LBv2 pre-mix at B2.

R4_CONDITIONS = [
    # --- Pure glutamate titration (no MOPS) ---
    {"label": "LBv2 ctrl",            "well": "A5",  "base": "LBv2", "Glutamate": 0,   "MOPS": 0},
    {"label": "LBv2+Glut100",         "well": "A6",  "base": "LBv2", "Glutamate": 100, "MOPS": 0},
    {"label": "LBv2+Glut125",         "well": "A7",  "base": "LBv2", "Glutamate": 125, "MOPS": 0},
    {"label": "LBv2+Glut150",         "well": "A8",  "base": "LBv2", "Glutamate": 150, "MOPS": 0},

    # --- Glutamate + MOPS 20 uL cross ---
    {"label": "LBv2+MOPS20",          "well": "B5",  "base": "LBv2", "Glutamate": 0,   "MOPS": 20},
    {"label": "LBv2+Glut75+MOPS",     "well": "B6",  "base": "LBv2", "Glutamate": 75,  "MOPS": 20},
    {"label": "LBv2+Glut100+MOPS",    "well": "B7",  "base": "LBv2", "Glutamate": 100, "MOPS": 20},
    {"label": "LBv2+Glut125+MOPS",    "well": "B8",  "base": "LBv2", "Glutamate": 125, "MOPS": 20},
    {"label": "LBv2+Glut150+MOPS",    "well": "C5",  "base": "LBv2", "Glutamate": 150, "MOPS": 20},

    # --- MOPS dose comparison ---
    {"label": "LBv2+MOPS40",          "well": "C6",  "base": "LBv2", "Glutamate": 0,   "MOPS": 40},
    {"label": "LBv2+Glut100+MOPS40",  "well": "C7",  "base": "LBv2", "Glutamate": 100, "MOPS": 40},

    # --- Replicate control ---
    {"label": "LBv2+Glut100 rep",     "well": "C8",  "base": "LBv2", "Glutamate": 100, "MOPS": 0},
]


def build_transfers(conditions):
    """Build transfer array for Monomer Bio workflow."""
    transfers = []
    supplement_names = ["Glutamate", "MOPS"]

    for cond in conditions:
        dst_well = cond["well"]
        base_reagent = BASE_TO_REAGENT[cond["base"]]

        supp_total = sum(cond.get(s, 0) for s in supplement_names)
        base_vol = 180 - supp_total

        if base_vol < 0:
            print(f"  WARNING: supplements exceed 180 uL for {dst_well}!")
            continue

        # Base media transfer
        transfers.append({
            "src_plate": "reagent", "src_well": base_reagent,
            "dst_plate": "experiment", "dst_well": dst_well,
            "volume": round(base_vol, 1), "new_tip": "once"
        })

        # Supplement transfers
        for supp in supplement_names:
            vol = cond.get(supp, 0)
            if vol > 0:
                transfers.append({
                    "src_plate": "reagent", "src_well": SUPP_TO_REAGENT[supp],
                    "dst_plate": "experiment", "dst_well": dst_well,
                    "volume": round(vol, 1), "new_tip": "once"
                })

    # Inoculum (always fresh tip, with mixing)
    for cond in conditions:
        transfers.append({
            "src_plate": "cell_culture_stock", "src_well": "A1",
            "dst_plate": "experiment", "dst_well": cond["well"],
            "volume": 20, "new_tip": "always",
            "pre_mix_volume": 50, "pre_mix_reps": 3,
            "post_mix_volume": 50, "post_mix_reps": 3
        })

    return transfers


def print_protocol_summary(conditions, transfers):
    """Print human-readable protocol summary."""
    print(f"\n{'='*70}")
    print("ROUND 4 PROTOCOL SUMMARY")
    print(f"{'='*70}")
    print(f"\nDesign: Glutamate titration UP + MOPS cross on LBv2")
    print(f"Hypothesis: Higher [glutamate] + MOPS buffering = optimal formulation")
    print(f"\n{'─'*70}")
    print(f"{'Well':<6} {'Label':<25} {'LBv2':>6} {'Glut':>6} {'MOPS':>6} {'Total':>6}")
    print(f"{'─'*70}")

    for cond in conditions:
        supp_total = cond.get("Glutamate", 0) + cond.get("MOPS", 0)
        base_vol = 180 - supp_total
        total = base_vol + supp_total + 20  # +20 inoculum
        print(f"{cond['well']:<6} {cond['label']:<25} {base_vol:>5.0f}  {cond.get('Glutamate', 0):>5.0f}  {cond.get('MOPS', 0):>5.0f}  {total:>5.0f}")

    media_transfers = [t for t in transfers if t["src_plate"] == "reagent"]
    inoc_transfers = [t for t in transfers if t["src_plate"] == "cell_culture_stock"]

    print(f"\n{'─'*70}")
    print(f"  Media/supplement transfers: {len(media_transfers)}")
    print(f"  Inoculum transfers:         {len(inoc_transfers)}")
    print(f"  TOTAL TRANSFERS:            {len(transfers)}")

    if len(transfers) > 40:
        print(f"  *** WARNING: exceeds 40 transfer limit by {len(transfers) - 40}! ***")
    else:
        print(f"  Within 40 transfer limit ({40 - len(transfers)} remaining)")

    print(f"\n{'='*70}")
    print("EXPERIMENTAL QUESTIONS:")
    print("  1. Does growth rate continue to increase past 100 uL glutamate?")
    print("  2. Does MOPS buffer sustain growth longer (higher AUC)?")
    print("  3. Is there a glutamate × MOPS interaction?")
    print("  4. What is the optimal glutamate:MOPS ratio?")
    print(f"{'='*70}")


if __name__ == "__main__":
    transfers = build_transfers(R4_CONDITIONS)
    print_protocol_summary(R4_CONDITIONS, transfers)

    # Save transfer JSON for Monomer Bio submission
    out_path = DATA_DIR / "round4_transfers.json"
    with open(out_path, "w") as f:
        json.dump(transfers, f, indent=2)
    print(f"\nTransfer JSON saved: {out_path}")
