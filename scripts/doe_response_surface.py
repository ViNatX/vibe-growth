"""
Design of Experiments — Response Surface Model

Maps reagent volumes (uL) to growth rate (mu) across both teams' experiment plates.
Uses StandardScaler + Ridge regression with cross-validated alpha.
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import Pipeline
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from pathlib import Path
from itertools import combinations
from scipy.optimize import minimize, LinearConstraint

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ── Build design matrix ──────────────────────────────────────────────────────

# Updated with 109-min matched-window mu values
cellai_data = [
    {"well": "B2",  "MOPS": 11, "Glucose": 0,  "Tryptone": 26, "YE": 14, "Glutamate": 0,  "NaCl": 0,  "mu": 0.708},
    {"well": "B3",  "MOPS": 7,  "Glucose": 0,  "Tryptone": 14, "YE": 20, "Glutamate": 0,  "NaCl": 0,  "mu": 0.716},
    {"well": "B4",  "MOPS": 9,  "Glucose": 0,  "Tryptone": 20, "YE": 26, "Glutamate": 0,  "NaCl": 0,  "mu": 0.734},
    {"well": "B5",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 26, "YE": 14, "Glutamate": 10, "NaCl": 0,  "mu": 0.762},
    {"well": "B6",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 14, "YE": 20, "Glutamate": 6,  "NaCl": 0,  "mu": 0.756},
    {"well": "B7",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 20, "YE": 26, "Glutamate": 8,  "NaCl": 0,  "mu": 0.750},
    {"well": "B8",  "MOPS": 12, "Glucose": 18, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.610},
    {"well": "B9",  "MOPS": 16, "Glucose": 12, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.740},
    {"well": "B10", "MOPS": 18, "Glucose": 8,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.691},
    {"well": "B11", "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.641},
]

vinatx_data = [
    {"well": "A1",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.674},
    {"well": "A2",  "MOPS": 0,  "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.613},
    {"well": "A3",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20, "mu": 0.534},
    {"well": "A4",  "MOPS": 0,  "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20, "mu": 0.549},
    {"well": "B1",  "MOPS": 20, "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.684},
    {"well": "B2v", "MOPS": 20, "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.659},
    {"well": "B3v", "MOPS": 20, "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20, "mu": 0.703},
    {"well": "B4v", "MOPS": 20, "Glucose": 20, "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 20, "mu": 0.660},
    {"well": "C1",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.683},
    {"well": "C2",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.725},
    {"well": "C3",  "MOPS": 0,  "Glucose": 0,  "Tryptone": 0,  "YE": 0,  "Glutamate": 0,  "NaCl": 0,  "mu": 0.625},
]

df = pd.DataFrame(cellai_data + vinatx_data)
factors = ["MOPS", "Glucose", "Tryptone", "YE", "Glutamate", "NaCl"]
X = df[factors].values
y = df["mu"].values

print(f"Design matrix: {X.shape[0]} experiments x {len(factors)} factors")
print(f"Response range: mu = {y.min():.3f} to {y.max():.3f} 1/hr\n")

# ── Fit model: scale -> polynomial -> ridge with CV ──────────────────────────

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

poly = PolynomialFeatures(degree=2, include_bias=True)
X_poly = poly.fit_transform(X_scaled)
feature_names = poly.get_feature_names_out(factors)

model = RidgeCV(alphas=np.logspace(-2, 3, 50), fit_intercept=False)
model.fit(X_poly, y)
y_pred = model.predict(X_poly)
r2 = r2_score(y, y_pred)

print(f"Model R² = {r2:.4f}  (Ridge alpha = {model.alpha_:.4f})")
print(f"\nTop coefficients:")
coef_df = pd.DataFrame({"term": feature_names, "coef": model.coef_})
coef_df["abs_coef"] = coef_df["coef"].abs()
coef_df = coef_df.sort_values("abs_coef", ascending=False)
for _, row in coef_df.head(12).iterrows():
    print(f"  {row['term']:>30s}: {row['coef']:+.4f}")

df["mu_pred"] = y_pred
df["residual"] = y - y_pred

print(f"\nPredicted vs Actual:")
for _, row in df.iterrows():
    print(f"  [{row['well']:>4s}] actual={row['mu']:.3f}  pred={row['mu_pred']:.3f}  resid={row['residual']:+.3f}")


# ── Helper: predict mu from raw volumes ──────────────────────────────────────

def predict_mu(volumes):
    """Predict mu from raw supplement volumes (uL)."""
    x_scaled = scaler.transform(volumes.reshape(1, -1))
    x_poly = poly.transform(x_scaled)
    return model.predict(x_poly)[0]


# ── Optimization ─────────────────────────────────────────────────────────────

def neg_mu(x):
    return -predict_mu(x)

# Bounds within tested ranges
bounds = [(0, 20), (0, 20), (0, 26), (0, 26), (0, 10), (0, 20)]
total_constraint = LinearConstraint(np.ones(len(factors)), lb=0, ub=80)

best_mu = -np.inf
best_x = None
np.random.seed(42)
for _ in range(2000):
    x0 = np.array([np.random.uniform(lo, hi) for lo, hi in bounds])
    if sum(x0) > 80:
        x0 = x0 * 80 / sum(x0)
    result = minimize(neg_mu, x0, bounds=bounds, constraints=[total_constraint], method="SLSQP")
    if result.success and -result.fun > best_mu:
        best_mu = -result.fun
        best_x = result.x

print(f"\n{'='*60}")
print(f"OPTIMAL REAGENT VOLUMES (within tested design space)")
print(f"{'='*60}")
total_supp = 0
for name, val, (lo, hi) in zip(factors, best_x, bounds):
    at_bound = " (at max)" if abs(val - hi) < 0.1 else (" (at min)" if val < 0.1 else "")
    print(f"  {name:>12s}: {val:5.1f} uL{at_bound}")
    total_supp += val
print(f"  {'Base media':>12s}: {180 - total_supp:.1f} uL")
print(f"  {'Inoculum':>12s}:  20.0 uL")
print(f"  {'Total':>12s}: 200.0 uL")
print(f"\n  Predicted mu: {best_mu:.3f} 1/hr")
print(f"  Best observed: {y.max():.3f} 1/hr (B5)")
print(f"  Improvement:   {(best_mu/y.max()-1)*100:+.1f}%")

# ── Figure 1: Model diagnostics + response surfaces ─────────────────────────

fig = plt.figure(figsize=(20, 14))

# Panel 1: Predicted vs Actual
ax = fig.add_subplot(2, 3, 1)
colors = ["#2196F3" if "B" in w and "v" not in w and w != "B1" else "#FF9800" for w in df["well"]]
ax.scatter(df["mu"], df["mu_pred"], c=colors, s=80, edgecolors="k", zorder=3)
for _, row in df.iterrows():
    ax.annotate(row["well"], (row["mu"], row["mu_pred"]), fontsize=7,
                xytext=(4, 4), textcoords="offset points")
lims = [0.5, 0.9]
ax.plot(lims, lims, "r--", alpha=0.5)
ax.set_xlabel("Actual mu (1/hr)")
ax.set_ylabel("Predicted mu (1/hr)")
ax.set_title(f"Model Fit (R² = {r2:.3f})\nBlue=CellAI, Orange=ViNatX")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.grid(True, alpha=0.3)

# Panels 2-6: Response surfaces for key factor pairs
key_pairs = [
    ("Tryptone", "Glutamate"),
    ("MOPS", "Glucose"),
    ("Tryptone", "YE"),
    ("MOPS", "NaCl"),
    ("Glutamate", "MOPS"),
]

for k, (fi, fj) in enumerate(key_pairs):
    ax = fig.add_subplot(2, 3, k + 2)
    i, j = factors.index(fi), factors.index(fj)

    grid_size = 60
    f1_range = np.linspace(bounds[i][0], bounds[i][1], grid_size)
    f2_range = np.linspace(bounds[j][0], bounds[j][1], grid_size)
    F1, F2 = np.meshgrid(f1_range, f2_range)

    Z = np.zeros_like(F1)
    for ii in range(grid_size):
        for jj in range(grid_size):
            x_test = best_x.copy()
            x_test[i] = F1[ii, jj]
            x_test[j] = F2[ii, jj]
            Z[ii, jj] = predict_mu(x_test)

    # Clip predictions to reasonable range
    Z = np.clip(Z, 0.4, 1.0)

    contour = ax.contourf(F1, F2, Z, levels=20, cmap="RdYlGn", vmin=0.5, vmax=0.9)
    plt.colorbar(contour, ax=ax, label="mu (1/hr)")
    ax.contour(F1, F2, Z, levels=8, colors="k", alpha=0.3, linewidths=0.5)

    ax.scatter(df[fi], df[fj], c=df["mu"], cmap="RdYlGn", edgecolors="k",
               s=60, zorder=3, vmin=0.5, vmax=0.9)
    ax.scatter(best_x[i], best_x[j], marker="*", c="red", s=250, zorder=4,
               edgecolors="k", linewidths=1.5)

    ax.set_xlabel(f"{fi} (uL)")
    ax.set_ylabel(f"{fj} (uL)")
    ax.set_title(f"{fi} vs {fj}")

plt.suptitle(
    f"Response Surface Model — Growth Rate Optimization\n"
    f"R² = {r2:.3f} | Predicted optimum mu = {best_mu:.3f} 1/hr | "
    f"Red star = predicted optimum",
    fontsize=14, fontweight="bold",
)
plt.tight_layout()
plt.savefig(FIG_DIR / "doe_response_surface.png", dpi=150, bbox_inches="tight")
plt.close()

# ── Figure 2: Main effects ──────────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
for k, factor in enumerate(factors):
    ax = axes.flat[k]
    idx = factors.index(factor)
    f_range = np.linspace(bounds[idx][0], bounds[idx][1], 100)
    mu_pred_line = []
    for val in f_range:
        x_test = best_x.copy()
        x_test[idx] = val
        mu_pred_line.append(predict_mu(x_test))

    mu_pred_line = np.clip(mu_pred_line, 0.4, 1.0)
    ax.plot(f_range, mu_pred_line, "b-", linewidth=2.5)
    ax.scatter(df[factor], df["mu"], c="orange", edgecolors="k", s=60, zorder=3)
    ax.axvline(best_x[idx], color="red", linestyle="--", alpha=0.7,
               label=f"Optimum = {best_x[idx]:.1f} uL")
    ax.set_xlabel(f"{factor} (uL)")
    ax.set_ylabel("mu (1/hr)")
    ax.set_title(f"Main Effect: {factor}")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.5, 0.9)

plt.suptitle("Main Effects — Growth Rate vs Individual Reagent Volumes\n"
             "(other factors held at optimum, orange dots = actual data)",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(FIG_DIR / "doe_main_effects.png", dpi=150, bbox_inches="tight")
plt.close()

df.to_csv(DATA_DIR / "doe_design_matrix.csv", index=False)
print(f"\nFigures saved: doe_response_surface.png, doe_main_effects.png")
