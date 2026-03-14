"""
Pitch deck: 5 slides telling the narrative arc of the hackathon.

Slide 1: The Challenge — what we started with
Slide 2: The Loop — how the autonomous cycle works
Slide 3: The Pivot — wrong fast, then right fast
Slide 4: The Result — glutamate dose-response, the 16% number
Slide 5: The Takeaway — what we learned about AI + bio
"""

import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures" / "pitch"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Shared style ─────────────────────────────────────────────────────────────

BG_COLOR = "#FAFAFA"
TITLE_COLOR = "#1A1A2E"
ACCENT_RED = "#E74C3C"
ACCENT_BLUE = "#2980B9"
ACCENT_GREEN = "#27AE60"
ACCENT_PURPLE = "#8E44AD"
ACCENT_ORANGE = "#E67E22"
MUTED_GRAY = "#7F8C8D"
DARK_TEXT = "#2C3E50"


def slide_figure():
    """Create a 16:9 figure with consistent styling."""
    fig = plt.figure(figsize=(20, 11.25), facecolor=BG_COLOR)
    return fig


def slide_title(fig, text, subtitle=None):
    """Add a large title and optional subtitle."""
    fig.text(0.5, 0.94, text, ha="center", va="top",
             fontsize=28, fontweight="bold", color=TITLE_COLOR,
             fontfamily="sans-serif")
    if subtitle:
        fig.text(0.5, 0.89, subtitle, ha="center", va="top",
                 fontsize=15, color=MUTED_GRAY, fontfamily="sans-serif",
                 style="italic")


def draw_rounded_box(ax, x, y, w, h, text, facecolor, edgecolor,
                     fontsize=12, text_color=DARK_TEXT, linewidth=2,
                     ha="center", fontweight="normal"):
    """Draw a rounded rectangle with centered text."""
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle="round,pad=0.02",
                         facecolor=facecolor, edgecolor=edgecolor,
                         linewidth=linewidth, transform=ax.transAxes,
                         clip_on=False)
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text,
            ha=ha if ha != "center" else "center", va="center",
            fontsize=fontsize, color=text_color, fontweight=fontweight,
            transform=ax.transAxes, wrap=True)


# ── Load data ────────────────────────────────────────────────────────────────

df = pd.read_csv(DATA_DIR / "vinatx_experiment_plate.csv", parse_dates=["timestamp"])
df = df.sort_values("timestamp").reset_index(drop=True)
t0_plate = df["timestamp"].iloc[0]
df["hours"] = (df["timestamp"] - t0_plate).dt.total_seconds() / 3600

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

r1_start_hr = df["hours"].iloc[0]
r2_mask = df["D1"].notna() & (df["D1"] > 0)
r2_start_hr = df.loc[r2_mask, "hours"].iloc[0]
r3_mask = df["F2"].notna() & (df["F2"] > 0)
r3_start_hr = df.loc[r3_mask, "hours"].iloc[0]


def fit_exponential_phase(times_hrs, od_values, window_hrs=2.0):
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    exp_mask = t <= (t[0] + window_hrs)
    t_exp, od_exp = t[exp_mask], od[exp_mask]
    if len(t_exp) < 3:
        return np.nan, np.nan
    slope, intercept, r, p, se = linregress(t_exp, np.log(od_exp))
    return slope, r**2


def fit_growth_rate(times_hrs, od_values):
    mask = np.isfinite(od_values) & (od_values > 0)
    t = np.array(times_hrs)[mask]
    od = np.array(od_values)[mask]
    if len(t) < 3:
        return np.nan, np.nan
    slope, intercept, r, p, se = linregress(t, np.log(od))
    return slope, r**2


# Compute all growth rates
r1_results = {}
for label, well in R1_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0)
    t_rel = (df.loc[m, "hours"] - r1_start_hr).values
    od = df.loc[m, well].values
    mu, r2 = fit_exponential_phase(t_rel, od)
    r1_results[label] = {"mu": mu, "R2": r2}

r3_results = {}
for label, well in R3_CONDITIONS.items():
    m = df[well].notna() & (df[well] > 0)
    if not m.any():
        continue
    t_rel = (df.loc[m, "hours"] - r3_start_hr).values
    od = df.loc[m, well].values
    mu, r2 = fit_growth_rate(t_rel, od)
    r3_results[label] = {"mu": mu, "R2": r2}


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 1: The Challenge
# ══════════════════════════════════════════════════════════════════════════════

fig1 = slide_figure()
slide_title(fig1, "Can an AI optimize a growth medium in 24 hours?",
            "Team ViNatX  —  Monomer Bio / Elnora AI Science Hackathon")

ax1 = fig1.add_axes([0.05, 0.05, 0.9, 0.78])
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 6)
ax1.axis("off")
ax1.set_facecolor(BG_COLOR)

# Three big elements: Organism, Robot, AI
# Organism
draw_rounded_box(ax1, 0.03, 0.55, 0.27, 0.35,
                 "Vibrio natriegens\n\nFastest-growing\nnon-pathogenic\nbacterium known\n\n~10 min doubling time",
                 facecolor="#E8F5E9", edgecolor=ACCENT_GREEN,
                 fontsize=14, fontweight="bold")

# Robot
draw_rounded_box(ax1, 0.36, 0.55, 0.27, 0.35,
                 "Monomer Bio Workcell\n\nOT Flex liquid handler\nTecan plate reader\nLiconic incubator\n\nAutomated OD600 reads",
                 facecolor="#E3F2FD", edgecolor=ACCENT_BLUE,
                 fontsize=14, fontweight="bold")

# AI
draw_rounded_box(ax1, 0.69, 0.55, 0.27, 0.35,
                 "Claude Code + MCP\n\nDesigns experiments\nSubmits protocols\nAnalyzes data\nProposes next round",
                 facecolor="#F3E5F5", edgecolor=ACCENT_PURPLE,
                 fontsize=14, fontweight="bold")

# The question at bottom
draw_rounded_box(ax1, 0.15, 0.08, 0.70, 0.22,
                 "THE GOAL\n\n"
                 "Starting from 12 candidate media formulations,\n"
                 "iterate autonomously to find the fastest-growing condition.\n\n"
                 "3 rounds  ·  33 conditions  ·  1 plate  ·  < 24 hours",
                 facecolor="#FFF8E1", edgecolor=ACCENT_ORANGE,
                 fontsize=14, fontweight="bold", linewidth=3)

# Arrows connecting the three
for x_start, x_end in [(0.31, 0.35), (0.64, 0.68)]:
    ax1.annotate("", xy=(x_end, 0.725), xytext=(x_start, 0.725),
                 xycoords="axes fraction", textcoords="axes fraction",
                 arrowprops=dict(arrowstyle="-|>", color=MUTED_GRAY,
                                 lw=3, mutation_scale=20))

fig1.savefig(FIG_DIR / "pitch_slide1.png", dpi=150, bbox_inches="tight",
             facecolor=BG_COLOR)
plt.close(fig1)
print("Saved pitch_slide1.png")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 2: The Loop — how the autonomous cycle works
# ══════════════════════════════════════════════════════════════════════════════

fig2 = slide_figure()
slide_title(fig2, "The Autonomous Loop",
            "Design → Execute → Analyze → Repeat — each round in ~2 hours")

ax2 = fig2.add_axes([0.05, 0.03, 0.9, 0.80])
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 6)
ax2.axis("off")
ax2.set_facecolor(BG_COLOR)

# Four nodes in a cycle — with gap between rows for arrows + timeline below
nodes = [
    (0.12, 0.55, 0.20, 0.28, "DESIGN\n\nClaude Code\n+ Elnora AI\npropose conditions\n+ plate layout",
     "#E3F2FD", ACCENT_BLUE),
    (0.55, 0.55, 0.20, 0.28, "EXECUTE\n\nRobot mixes media,\ninoculates,\nincubates",
     "#E8F5E9", ACCENT_GREEN),
    (0.55, 0.18, 0.20, 0.28, "ANALYZE\n\nAI fits growth\ncurves, ranks\nconditions",
     "#FFF3E0", ACCENT_ORANGE),
    (0.12, 0.18, 0.20, 0.28, "DECIDE\n\nHuman reviews,\nAI redesigns\nnext round",
     "#F3E5F5", ACCENT_PURPLE),
]

for x, y, w, h, text, fc, ec in nodes:
    draw_rounded_box(ax2, x, y, w, h, text, facecolor=fc, edgecolor=ec,
                     fontsize=14, fontweight="bold", linewidth=2.5)

# Arrows forming the cycle (clockwise)
arrow_style = dict(arrowstyle="-|>", color="#555555", lw=3, mutation_scale=22,
                   connectionstyle="arc3,rad=0.0")
# Design → Execute (top, right)
ax2.annotate("", xy=(0.54, 0.69), xytext=(0.33, 0.69),
             xycoords="axes fraction", textcoords="axes fraction",
             arrowprops={**arrow_style})
ax2.text(0.435, 0.73, "Protocol\nsubmission", ha="center", va="bottom",
         fontsize=11, color=MUTED_GRAY, transform=ax2.transAxes, style="italic")

# Execute → Analyze (right, down)
ax2.annotate("", xy=(0.65, 0.47), xytext=(0.65, 0.54),
             xycoords="axes fraction", textcoords="axes fraction",
             arrowprops={**arrow_style})
ax2.text(0.72, 0.505, "OD600\ndata", ha="left", va="center",
         fontsize=11, color=MUTED_GRAY, transform=ax2.transAxes, style="italic")

# Analyze → Decide (bottom, left)
ax2.annotate("", xy=(0.33, 0.32), xytext=(0.54, 0.32),
             xycoords="axes fraction", textcoords="axes fraction",
             arrowprops={**arrow_style})
ax2.text(0.435, 0.28, "Growth rates\n+ rankings", ha="center", va="top",
         fontsize=11, color=MUTED_GRAY, transform=ax2.transAxes, style="italic")

# Decide → Design (left, up)
ax2.annotate("", xy=(0.22, 0.54), xytext=(0.22, 0.47),
             xycoords="axes fraction", textcoords="axes fraction",
             arrowprops={**arrow_style})
ax2.text(0.09, 0.505, "New\nhypothesis", ha="right", va="center",
         fontsize=11, color=MUTED_GRAY, transform=ax2.transAxes, style="italic")

# Timeline at bottom — well below the boxes
timeline_y = 0.04
rounds_info = [
    (0.15, "Round 1", "12 conditions — Basal media screen", ACCENT_BLUE),
    (0.43, "Round 2", "9 conditions — Supplement optimization", ACCENT_GREEN),
    (0.71, "Round 3", "12 conditions — Glutamate titration", ACCENT_ORANGE),
]
for x, title, desc, color in rounds_info:
    ax2.text(x + 0.07, timeline_y + 0.06, title, ha="center", va="bottom",
             fontsize=14, fontweight="bold", color=color, transform=ax2.transAxes)
    ax2.text(x + 0.07, timeline_y - 0.01, desc, ha="center", va="top",
             fontsize=10, color=MUTED_GRAY, transform=ax2.transAxes)

# Connecting line for timeline
ax2.plot([0.15, 0.85], [timeline_y + 0.055, timeline_y + 0.055],
         color="#CCCCCC", linewidth=2, transform=ax2.transAxes, zorder=0)
for x, _, _, color in rounds_info:
    ax2.plot(x + 0.07, timeline_y + 0.055, "o", color=color, markersize=12,
             transform=ax2.transAxes, zorder=1)

fig2.savefig(FIG_DIR / "pitch_slide2.png", dpi=150, bbox_inches="tight",
             facecolor=BG_COLOR)
plt.close(fig2)
print("Saved pitch_slide2.png")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 3: The Pivot — wrong fast, then right fast
# ══════════════════════════════════════════════════════════════════════════════

fig3 = slide_figure()
slide_title(fig3, "Wrong Fast, Then Right Fast",
            "The AI's first hypothesis was wrong — but rapid iteration let us pivot")

gs3 = GridSpec(1, 2, figure=fig3, left=0.06, right=0.94, bottom=0.08, top=0.82,
               wspace=0.12)

# Left panel: R1 simplified bar chart — just top 5 and bottom 3
ax3a = fig3.add_subplot(gs3[0, 0])
ax3a.set_facecolor(BG_COLOR)

# Sort R1 by mu
r1_sorted = sorted(r1_results.items(), key=lambda x: x[1]["mu"] if np.isfinite(x[1]["mu"]) else -999)
# Show all 12 but highlight the story: defined media at bottom, LBv2 at top
r1_labels = [label for label, _ in r1_sorted]
r1_mus = [data["mu"] for _, data in r1_sorted]

# Color by group
group_map = {
    "Def-Min": "defined", "DM+Gluc": "defined", "DM+NaCl": "defined", "DM+Gluc+NaCl": "defined",
    "DM+MOPS": "buffered", "DM+Gluc+MOPS": "buffered", "DM+MOPS+NaCl": "buffered",
    "DM+Gluc+MOPS+NaCl": "buffered",
    "NBxCyclone": "rich", "LBv2": "rich", "Semi-Def": "rich", "Def-Glycerol": "rich",
}
group_colors = {"defined": "#FFCC80", "buffered": "#EF5350", "rich": "#C62828"}
bar_colors = [group_colors[group_map[label]] for label in r1_labels]

bars = ax3a.barh(r1_labels, r1_mus, color=bar_colors, edgecolor="gray", alpha=0.85)
for bar, mu in zip(bars, r1_mus):
    if np.isfinite(mu):
        ax3a.text(bar.get_width() - 0.01, bar.get_y() + bar.get_height() / 2,
                  f"{mu:.2f}", va="center", ha="right", fontsize=9,
                  color="white", fontweight="bold")

ax3a.set_xlabel("Growth Rate μ (1/hr)", fontsize=13)
ax3a.set_title("Round 1: The AI tested 12 media\nRich media (LBv2) crushed everything",
               fontsize=14, fontweight="bold", color=DARK_TEXT)
ax3a.set_xlim(0, 0.85)
ax3a.tick_params(axis="y", labelsize=10)
ax3a.grid(True, axis="x", alpha=0.3)

# Legend
legend_elements = [
    Patch(facecolor="#FFCC80", edgecolor="gray", label="Defined Minimal"),
    Patch(facecolor="#EF5350", edgecolor="gray", label="Defined + Buffer"),
    Patch(facecolor="#C62828", edgecolor="gray", label="Rich Media"),
]
ax3a.legend(handles=legend_elements, loc="lower right", fontsize=9)

# Right panel: The pivot narrative
ax3b = fig3.add_subplot(gs3[0, 1])
ax3b.axis("off")
ax3b.set_facecolor(BG_COLOR)

# "What AI designed" box (crossed out / wrong)
draw_rounded_box(ax3b, 0.02, 0.72, 0.96, 0.25,
                 "AI'S INITIAL HYPOTHESIS\n\n"
                 "\"Optimize glucose, NaCl, and MOPS buffering\n"
                 "on Defined Minimal base medium\"\n\n"
                 "✗  WRONG — defined media were the worst performers",
                 facecolor="#FFEBEE", edgecolor=ACCENT_RED,
                 fontsize=13, fontweight="bold", linewidth=3)

# Arrow pointing DOWN from initial hypothesis to pivot
ax3b.annotate("", xy=(0.5, 0.63), xytext=(0.5, 0.70),
              xycoords="axes fraction", textcoords="axes fraction",
              arrowprops=dict(arrowstyle="-|>", color=ACCENT_GREEN, lw=4,
                              mutation_scale=25))
ax3b.text(0.58, 0.665, "pivot in\nminutes", ha="left", va="center",
          fontsize=12, color=ACCENT_GREEN, fontweight="bold",
          transform=ax3b.transAxes, style="italic")

# "What data told us" box
draw_rounded_box(ax3b, 0.02, 0.35, 0.96, 0.25,
                 "DATA-DRIVEN PIVOT\n\n"
                 "\"LBv2 is the winner. Now optimize\n"
                 "nitrogen supplementation on LBv2.\"\n\n"
                 "✓  Round 2: tested glutamate, tryptone, YE on LBv2",
                 facecolor="#E8F5E9", edgecolor=ACCENT_GREEN,
                 fontsize=13, fontweight="bold", linewidth=3)

# Arrow pointing DOWN from pivot to confirmation
ax3b.annotate("", xy=(0.5, 0.26), xytext=(0.5, 0.33),
              xycoords="axes fraction", textcoords="axes fraction",
              arrowprops=dict(arrowstyle="-|>", color=ACCENT_GREEN, lw=4,
                              mutation_scale=25))

# "What we confirmed" box
draw_rounded_box(ax3b, 0.02, 0.01, 0.96, 0.22,
                 "CONFIRMED IN ROUND 3\n\n"
                 "\"Glutamate is the key driver.\n"
                 "Clear dose-response: more glutamate = faster growth.\"\n\n"
                 "Glucose, tryptone, YE → no help or actively harmful",
                 facecolor="#E3F2FD", edgecolor=ACCENT_BLUE,
                 fontsize=13, fontweight="bold", linewidth=3)

fig3.savefig(FIG_DIR / "pitch_slide3.png", dpi=150, bbox_inches="tight",
             facecolor=BG_COLOR)
plt.close(fig3)
print("Saved pitch_slide3.png")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 4: The Result — clean dose-response
# ══════════════════════════════════════════════════════════════════════════════

fig4 = slide_figure()
slide_title(fig4, "The Result: A Clear Dose-Response",
            "More glutamate = faster growth — and we haven't hit the ceiling")

gs4 = GridSpec(1, 2, figure=fig4, left=0.06, right=0.94, bottom=0.08, top=0.82,
               wspace=0.15, width_ratios=[1.2, 1.0])

# Left panel: Only the glutamate titration wells (the clean story)
ax4a = fig4.add_subplot(gs4[0, 0])
ax4a.set_facecolor(BG_COLOR)

# Select just the dose-response wells + control
dose_wells = ["LBv2 ctrl", "LBv2+Glut10", "LBv2+Glut25",
              "LBv2+Glut50", "LBv2+Glut75", "LBv2+Glut100"]
dose_mM = [0, 10, 25, 50, 75, 100]
dose_mus = [r3_results[label]["mu"] for label in dose_wells]
dose_r2s = [r3_results[label]["R2"] for label in dose_wells]

# Bar colors gradient
bar_colors_dose = [plt.cm.YlOrRd(0.15 + 0.75 * mM / 100) for mM in dose_mM]

bars = ax4a.bar(
    [f"{mM} mM" for mM in dose_mM], dose_mus,
    color=bar_colors_dose, edgecolor="gray", alpha=0.9, width=0.65,
)

# Value labels on top of bars
for bar, mu, r2 in zip(bars, dose_mus, dose_r2s):
    ax4a.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
              f"μ = {mu:.3f}", ha="center", va="bottom",
              fontsize=11, fontweight="bold", color=DARK_TEXT)

ax4a.set_xlabel("Glutamate Concentration (mM)", fontsize=14)
ax4a.set_ylabel("Growth Rate μ (1/hr)", fontsize=14)
ax4a.set_title("Glutamate Titration on LBv2",
               fontsize=15, fontweight="bold", color=DARK_TEXT)
ax4a.set_ylim(0, 0.78)
ax4a.grid(True, axis="y", alpha=0.3)
ax4a.tick_params(labelsize=11)

# Add a trend arrow
ax4a.annotate("", xy=(4.8, 0.70), xytext=(0.5, 0.60),
              arrowprops=dict(arrowstyle="-|>", color=ACCENT_GREEN, lw=3,
                              connectionstyle="arc3,rad=0.15",
                              mutation_scale=20))
ax4a.text(2.5, 0.72, "Not plateaued →", fontsize=13, color=ACCENT_GREEN,
          fontweight="bold", style="italic")

# Right panel: The key numbers
ax4b = fig4.add_subplot(gs4[0, 1])
ax4b.axis("off")
ax4b.set_facecolor(BG_COLOR)

# Big number: 16% improvement
ax4b.text(0.5, 0.90, "+16%", ha="center", va="top",
          fontsize=72, fontweight="bold", color=ACCENT_GREEN,
          transform=ax4b.transAxes)
ax4b.text(0.5, 0.70, "growth rate improvement", ha="center", va="top",
          fontsize=18, color=DARK_TEXT, transform=ax4b.transAxes)

# Before/after comparison
draw_rounded_box(ax4b, 0.05, 0.42, 0.40, 0.18,
                 "BEFORE\n\nBare LBv2\nμ = 0.592 /hr",
                 facecolor="#FFEBEE", edgecolor=ACCENT_RED,
                 fontsize=13, fontweight="bold")

ax4b.annotate("", xy=(0.54, 0.51), xytext=(0.46, 0.51),
              xycoords="axes fraction", textcoords="axes fraction",
              arrowprops=dict(arrowstyle="-|>", color=DARK_TEXT, lw=3,
                              mutation_scale=20))

draw_rounded_box(ax4b, 0.55, 0.42, 0.40, 0.18,
                 "AFTER\n\nLBv2 + Glut 100mM\nμ = 0.688 /hr",
                 facecolor="#E8F5E9", edgecolor=ACCENT_GREEN,
                 fontsize=13, fontweight="bold")

# Additional findings
findings = [
    "Glucose actively hurts growth (overflow metabolism)",
    "Extra tryptone / yeast extract don't help",
    "Good reproducibility: R2 vs R3 replicate ≈ 3% diff",
    "Dose-response hasn't plateaued — room to push higher",
]
for i, finding in enumerate(findings):
    marker_color = ACCENT_RED if i == 0 else (ACCENT_ORANGE if i == 1 else ACCENT_GREEN)
    ax4b.text(0.08, 0.32 - i * 0.08, f"▸  {finding}",
              fontsize=11, color=DARK_TEXT, transform=ax4b.transAxes,
              va="top")

fig4.savefig(FIG_DIR / "pitch_slide4.png", dpi=150, bbox_inches="tight",
             facecolor=BG_COLOR)
plt.close(fig4)
print("Saved pitch_slide4.png")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 5: The Takeaway
# ══════════════════════════════════════════════════════════════════════════════

fig5 = slide_figure()
slide_title(fig5, "What We Learned About AI + Biology",
            "The value isn't in getting it right the first time — it's in iterating fast")

ax5 = fig5.add_axes([0.05, 0.03, 0.9, 0.80])
ax5.set_xlim(0, 10)
ax5.set_ylim(0, 6)
ax5.axis("off")
ax5.set_facecolor(BG_COLOR)

# Three columns of insights
col_data = [
    ("AI Got the First\nHypothesis Wrong", ACCENT_BLUE, "#E3F2FD",
     "The initial DOE focused on glucose,\n"
     "NaCl, and buffering on Defined Minimal.\n\n"
     "The real drivers turned out to be\n"
     "nitrogen source (glutamate) and\n"
     "base media choice (LBv2).\n\n"
     "But AI helped us figure that out\n"
     "midstream and pivot in minutes,\n"
     "not days."),
    ("Integration Was\nthe Bottleneck", ACCENT_GREEN, "#E8F5E9",
     "Connecting MCP servers, Notion,\n"
     "Monomer Cloud, and the workcell\n"
     "took real effort up front.\n\n"
     "But once the pipeline was built,\n"
     "iteration accelerated dramatically.\n"
     "Rounds 2 and 3 moved much faster\n"
     "than Round 1.\n\n"
     "Context consolidation was key."),
    ("AI Makes Novel\nMistakes (and Bets)", ACCENT_PURPLE, "#F3E5F5",
     "AI made mistakes a grad student\n"
     "would also make: wrong time\n"
     "windows, mismatched comparisons.\n\n"
     "But it also lacked the biases that\n"
     "would discourage novel experiments.\n"
     "It proposed conditions we might\n"
     "have vetoed as 'too unusual.'\n\n"
     "Some of those could have worked."),
]

for i, (title, color, bgcolor, body) in enumerate(col_data):
    x = 0.03 + i * 0.33
    # Title
    draw_rounded_box(ax5, x, 0.78, 0.29, 0.18, title,
                     facecolor=bgcolor, edgecolor=color,
                     fontsize=16, fontweight="bold", linewidth=3)
    # Body box (draw empty)
    body_box = FancyBboxPatch((x, 0.30), 0.29, 0.44,
                              boxstyle="round,pad=0.02",
                              facecolor="white", edgecolor=color,
                              linewidth=1.5, transform=ax5.transAxes,
                              clip_on=False)
    ax5.add_patch(body_box)
    # Text aligned to top of box
    ax5.text(x + 0.145, 0.70, body,
             ha="center", va="top",
             fontsize=12, color=DARK_TEXT,
             transform=ax5.transAxes, linespacing=1.5)

# Bottom banner
draw_rounded_box(ax5, 0.08, 0.05, 0.84, 0.13,
                 "3 rounds  ·  33 conditions  ·  < 24 hours  ·  "
                 "16% growth rate improvement  ·  dose-response not plateaued\n"
                 "The right question isn't \"Was AI right?\" — "
                 "it's \"Did AI help us learn faster?\"  The answer is yes.",
                 facecolor="#FFF8E1", edgecolor=ACCENT_ORANGE,
                 fontsize=15, fontweight="bold", linewidth=3)

fig5.savefig(FIG_DIR / "pitch_slide5.png", dpi=150, bbox_inches="tight",
             facecolor=BG_COLOR)
plt.close(fig5)
print("Saved pitch_slide5.png")


# ══════════════════════════════════════════════════════════════════════════════
# BUILD PDF
# ══════════════════════════════════════════════════════════════════════════════

from PIL import Image

slide_files = [
    FIG_DIR / "pitch_slide1.png",
    FIG_DIR / "pitch_slide2.png",
    FIG_DIR / "pitch_slide3.png",
    FIG_DIR / "pitch_slide4.png",
    FIG_DIR / "pitch_slide5.png",
]

pdf_path = FIG_DIR / "ViNatX_Pitch.pdf"

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

print(f"\nPitch deck saved: {pdf_path}")
