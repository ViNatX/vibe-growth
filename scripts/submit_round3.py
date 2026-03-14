"""
Round 3 protocol builder.

Takes the output of analyze_round2.py (DOE optima + proposals) and builds
a Monomer Bio workflow definition for submission.

Usage:
  pixi run python scripts/submit_round3.py

This script:
1. Reads the Round 3 proposals from doe_design_matrix_v2.csv and DOE results
2. Maps proposals to reagent plate wells
3. Calculates transfer volumes (basal media fills gaps, no water)
4. Verifies total transfers <= 40
5. Prints the workflow definition code ready for submission
"""

import json
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"

# ViNatX Reagent Plate well map
REAGENT_MAP = {
    "Defined-Minimal Media": "A1",  # Def-Min base
    "MOPS (400 mM, pH 7.0)": "A2",
    "LBv2": "B2",
    "Glutamate": "C1",
    "Tryptone": "C2",
    "Yeast Extract": "C3",
    "High Buffer Defined Media": "C4",
}

# Base media -> reagent well (ViNatX Reagent Plate ONLY)
BASE_TO_REAGENT = {
    "Def-Min": "A1",
    "LBv2": "B2",
    "HBDef": "C4",
    "Semi-Def": "B3",
}

# Supplement -> reagent well (ViNatX Reagent Plate ONLY)
# Do NOT add CellAI-only reagents (KH2PO4, MgSO4, Trace Metals, FeSO4, Na Citrate)
SUPP_TO_REAGENT = {
    "MOPS": "A2",
    "Glucose": "A4",
    "Tryptone": "C2",
    "YE": "C3",
    "Glutamate": "C1",
    "NaCl": "A3",
}


def build_transfer_array(conditions, target_wells):
    """
    Build transfer array for a set of conditions.

    Each condition: {base, MOPS, Glucose, Tryptone, YE, Glutamate, NaCl}
    Volumes are adjusted so base media fills to 180 uL total (no water).
    """
    transfers = []
    supplement_names = ["MOPS", "Glucose", "Tryptone", "YE", "Glutamate", "NaCl"]

    # Group transfers by source for tip reuse
    base_transfers = {}  # reagent_well -> [(dst_well, volume)]
    supp_transfers = {}  # reagent_well -> [(dst_well, volume)]

    for cond, dst_well in zip(conditions, target_wells):
        base = cond["base"]
        base_reagent = BASE_TO_REAGENT.get(base)
        if not base_reagent:
            print(f"  WARNING: base {base} not available on reagent plate!")
            continue

        # Calculate supplement total
        supp_total = sum(cond.get(s, 0) for s in supplement_names)
        base_vol = 180 - supp_total  # Fill remainder with base

        if base_vol < 0:
            print(f"  WARNING: supplements exceed 180 uL for {dst_well}!")
            continue

        # Base media transfer
        if base_reagent not in base_transfers:
            base_transfers[base_reagent] = []
        base_transfers[base_reagent].append((dst_well, base_vol))

        # Supplement transfers
        for supp in supplement_names:
            vol = cond.get(supp, 0)
            if vol > 0:
                reagent_well = SUPP_TO_REAGENT.get(supp)
                if not reagent_well:
                    print(f"  WARNING: {supp} not available on reagent plate!")
                    continue
                if reagent_well not in supp_transfers:
                    supp_transfers[reagent_well] = []
                supp_transfers[reagent_well].append((dst_well, vol))

    # Build ordered transfer list: bases first, then supplements, then inoculum
    # Bases
    for src_well, dsts in base_transfers.items():
        for dst_well, vol in dsts:
            transfers.append({
                "src_plate": "reagent", "src_well": src_well,
                "dst_plate": "experiment", "dst_well": dst_well,
                "volume": round(vol, 1), "new_tip": "once"
            })

    # Supplements
    for src_well, dsts in supp_transfers.items():
        for dst_well, vol in dsts:
            transfers.append({
                "src_plate": "reagent", "src_well": src_well,
                "dst_plate": "experiment", "dst_well": dst_well,
                "volume": round(vol, 1), "new_tip": "once"
            })

    # Inoculum (always fresh tip, with mixing)
    for dst_well in target_wells:
        transfers.append({
            "src_plate": "cell_culture_stock", "src_well": "A1",
            "dst_plate": "experiment", "dst_well": dst_well,
            "volume": 20, "new_tip": "always",
            "pre_mix_volume": 50, "pre_mix_reps": 3,
            "post_mix_volume": 50, "post_mix_reps": 3
        })

    return transfers


def print_protocol_summary(conditions, target_wells, transfers):
    """Print human-readable protocol summary."""
    supplement_names = ["MOPS", "Glucose", "Tryptone", "YE", "Glutamate", "NaCl"]

    print(f"\n{'='*70}")
    print("ROUND 3 PROTOCOL SUMMARY")
    print(f"{'='*70}")

    for cond, well in zip(conditions, target_wells):
        base = cond["base"]
        supp_total = sum(cond.get(s, 0) for s in supplement_names)
        base_vol = 180 - supp_total
        supps = [f"{s} {cond[s]}" for s in supplement_names if cond.get(s, 0) > 0]
        print(f"  {well}: {base} {base_vol:.0f} uL + {' + '.join(supps) if supps else '(no supplements)'}")

    media_transfers = [t for t in transfers if t["src_plate"] == "reagent"]
    inoc_transfers = [t for t in transfers if t["src_plate"] == "cell_culture_stock"]

    print(f"\n  Media/supplement transfers: {len(media_transfers)}")
    print(f"  Inoculum transfers: {len(inoc_transfers)}")
    print(f"  TOTAL TRANSFERS: {len(transfers)}")

    if len(transfers) > 40:
        print(f"  *** WARNING: exceeds 40 transfer limit by {len(transfers) - 40}! ***")
        print(f"  Consider reducing wells or consolidating supplements.")
    else:
        print(f"  Within 40 transfer limit ({40 - len(transfers)} remaining)")


if __name__ == "__main__":
    # Example: will be called with actual DOE results
    print("This script provides helper functions for Round 3 protocol building.")
    print("Run analyze_round2.py first, then use these functions to build the protocol.")
