"""
Fetch latest OD data from Monomer Bio and re-run the equal-time comparison.

This script is designed to be called by Claude's /loop skill to auto-refresh.
It prints a compact summary so each loop iteration is easy to scan.
"""
# This script is meant to be run via: pixi run python scripts/refresh_and_analyze.py
# But the actual data fetching needs MCP access, which only Claude has.
# So this just re-runs the analysis on whatever CSVs are on disk.

import subprocess
import sys

# Run both analysis scripts
for script in ["analyze_all_plates.py", "compare_equal_time.py"]:
    print(f"\n{'='*60}")
    print(f"  Running {script}")
    print(f"{'='*60}")
    result = subprocess.run(
        [sys.executable, f"scripts/{script}"],
        capture_output=False,
    )
    if result.returncode != 0:
        print(f"  ERROR: {script} failed with code {result.returncode}")
