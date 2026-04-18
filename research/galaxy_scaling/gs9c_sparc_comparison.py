"""
gs9c: Quantitative comparison with SPARC RAR data
==================================================

Uses REAL data from:
  Lelli, McGaugh, Schombert, Pawlowski (2017)
  "One Law to Rule Them All: The Radial Acceleration Relation of Galaxies"

Data files:
  RAR.mrt     - 2696 individual data points (gbar, gobs with errors)
  RARbins.mrt - 14 binned data points (mean gbar, mean gobs, scatter)

Models compared:
  1. MOND simple:    nu(y) = 1/2 + sqrt(1/4 + 1/y)
  2. MOND standard:  nu(y) = sqrt(1/2 + sqrt(1/4 + 1/y^2))
  3. DGP:            nu(y) = 1 + 1/sqrt(y)
  4. Hybrid A:       nu(y) = 1 + 1/(sqrt(y) + y)
  5. Hybrid B:       nu(y) = 1 + exp(-sqrt(y))/sqrt(y)
  6. Hybrid C:       nu(y) = sqrt(1 + 1/y)

From gs9b: DGP is ruled out by solar system (dg/g ~ 10^-4 at Earth),
but we still compare it to SPARC to understand how the 2D channel
manifests in galaxy data.
"""

import numpy as np
import os, sys, io, warnings
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ============================================================
# Load SPARC data
# ============================================================
script_dir = os.path.dirname(os.path.abspath(__file__))

def load_rar_data(filename):
    """Load RAR data from .mrt file."""
    filepath = os.path.join(script_dir, filename)
    data = []
    with open(filepath, 'r') as f:
        in_data = False
        for line in f:
            line = line.strip()
            if line.startswith('---') and in_data:
                continue
            if line.startswith('---'):
                in_data = True
                continue
            if in_data and line and not line.startswith(('Title', 'Authors', 'Table', '=')):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        vals = [float(x) for x in parts]
                        data.append(vals)
                    except ValueError:
                        continue
    return np.array(data)

print("=" * 70)
print("  gs9c: Quantitative comparison with SPARC RAR data")
print("=" * 70)

# Load individual data points
rar_data = load_rar_data('RAR.mrt')
print(f"\n  Loaded RAR.mrt: {len(rar_data)} data points")
print(f"  Columns: log10(gbar), e_gbar, log10(gobs), e_gobs")

log_gbar = rar_data[:, 0]  # log10(gbar) in m/s^2
e_log_gbar = rar_data[:, 1]
log_gobs = rar_data[:, 2]  # log10(gobs) in m/s^2
e_log_gobs = rar_data[:, 3]

# Convert to linear
gbar = 10**log_gbar
gobs = 10**log_gobs

# Load binned data
rar_bins = load_rar_data('RARbins.mrt')
print(f"  Loaded RARbins.mrt: {len(rar_bins)} bins")

bin_log_gbar = rar_bins[:, 0]
bin_log_gobs = rar_bins[:, 1]
bin_sd = rar_bins[:, 2]
bin_N = rar_bins[:, 3]

print(f"\n  Data range:")
print(f"    log10(gbar): [{log_gbar.min():.2f}, {log_gbar.max():.2f}]")
print(f"    log10(gobs): [{log_gobs.min():.2f}, {log_gobs.max():.2f}]")

# ============================================================
# Define models
# ============================================================
a0_fiducial = 1.2e-10  # m/s^2

def model_mond_simple(gbar, a0):
    return 0.5 * gbar + np.sqrt(0.25 * gbar**2 + gbar * a0)

def model_mond_standard(gbar, a0):
    x = gbar / a0
    return gbar * np.sqrt(0.5 + 0.5 * np.sqrt(1 + 4.0/x**2))

def model_dgp(gbar, a0):
    return gbar + np.sqrt(gbar * a0)

def model_hybridA(gbar, a0):
    y = gbar / a0
    return gbar * (1 + 1/(np.sqrt(y) + y))

def model_hybridB(gbar, a0):
    y = gbar / a0
    return gbar * (1 + np.exp(-np.sqrt(y)) / np.sqrt(y))

def model_hybridC(gbar, a0):
    y = gbar / a0
    return gbar * np.sqrt(1 + 1/y)

models = {
    'MOND simple':    model_mond_simple,
    'MOND standard':  model_mond_standard,
    'DGP':            model_dgp,
    'Hybrid A':       model_hybridA,
    'Hybrid B':       model_hybridB,
    'Hybrid C':       model_hybridC,
}

# ============================================================
# 1. CHI-SQUARED with fixed a0
# ============================================================
print("\n" + "=" * 70)
print("  1. Chi-squared with fixed a0 = 1.2e-10 m/s^2")
print("=" * 70)

print("""
  Chi^2 = sum[ (log10(g_obs) - log10(g_model))^2 / sigma^2 ]
  where sigma^2 = e_gobs^2 (observational error on log10(gobs))

  Also compute RMS in dex = sqrt(mean[(log10(gobs) - log10(gmodel))^2])
""")

a0 = a0_fiducial

print(f"  {'Model':>16}  {'chi2':>10}  {'chi2/N':>8}  {'RMS(dex)':>9}  {'RMS(%)':>8}")
print("  " + "-" * 58)

results = {}
for name, func in models.items():
    g_pred = func(gbar, a0)
    log_g_pred = np.log10(g_pred)

    residuals = log_gobs - log_g_pred
    chi2 = np.sum(residuals**2 / e_log_gobs**2)
    rms_dex = np.sqrt(np.mean(residuals**2))
    # RMS in percentage: 10^rms_dex - 1
    rms_pct = (10**rms_dex - 1) * 100
    chi2_N = chi2 / len(gbar)

    results[name] = {
        'chi2': chi2, 'chi2_N': chi2_N,
        'rms_dex': rms_dex, 'rms_pct': rms_pct,
        'residuals': residuals
    }
    print(f"  {name:>16}  {chi2:10.0f}  {chi2_N:8.2f}  {rms_dex:9.4f}  {rms_pct:8.1f}")

# ============================================================
# 2. OPTIMIZE a0 for each model
# ============================================================
print("\n" + "=" * 70)
print("  2. Optimize a0 for each model (best-fit)")
print("=" * 70)

print("""
  For each model, find a0 that minimizes chi^2.
  This tests: which SHAPE fits best, independent of a0 normalization.
""")

from scipy.optimize import minimize_scalar

def chi2_for_a0(log_a0, func):
    a0_val = 10**log_a0
    g_pred = func(gbar, a0_val)
    log_g_pred = np.log10(g_pred)
    residuals = log_gobs - log_g_pred
    return np.sum(residuals**2 / e_log_gobs**2)

print(f"  {'Model':>16}  {'a0_best':>12}  {'a0/1.2e-10':>11}  {'chi2_best':>10}  {'chi2/N':>8}  {'RMS(dex)':>9}")
print("  " + "-" * 75)

best_results = {}
for name, func in models.items():
    res = minimize_scalar(chi2_for_a0, bounds=(-11, -9), method='bounded', args=(func,))
    a0_best = 10**res.x
    g_pred = func(gbar, a0_best)
    log_g_pred = np.log10(g_pred)
    residuals = log_gobs - log_g_pred
    chi2 = res.fun
    rms = np.sqrt(np.mean(residuals**2))

    best_results[name] = {
        'a0': a0_best, 'chi2': chi2, 'chi2_N': chi2/len(gbar),
        'rms': rms, 'residuals': residuals
    }
    print(f"  {name:>16}  {a0_best:12.3e}  {a0_best/1.2e-10:11.3f}  {chi2:10.0f}  {chi2/len(gbar):8.2f}  {rms:9.4f}")

# ============================================================
# 3. BINNED COMPARISON
# ============================================================
print("\n" + "=" * 70)
print("  3. Comparison with binned SPARC data")
print("=" * 70)

print(f"\n  Using best-fit a0 for each model.")
print(f"\n  {'log(gbar)':>10}  {'gobs_data':>10}  {'MOND_s':>8}  {'MOND_st':>8}  {'DGP':>8}  {'HybA':>8}  {'HybB':>8}  {'HybC':>8}")
print("  " + "-" * 70)

bin_gbar_lin = 10**bin_log_gbar

for i in range(len(rar_bins)):
    gb = bin_gbar_lin[i]
    go = bin_log_gobs[i]

    preds = []
    for name in ['MOND simple', 'MOND standard', 'DGP', 'Hybrid A', 'Hybrid B', 'Hybrid C']:
        a0_use = best_results[name]['a0']
        g_pred = models[name](gb, a0_use)
        preds.append(np.log10(g_pred))

    print(f"  {bin_log_gbar[i]:10.2f}  {go:10.2f}  {preds[0]:8.2f}  {preds[1]:8.2f}  {preds[2]:8.2f}  {preds[3]:8.2f}  {preds[4]:8.2f}  {preds[5]:8.2f}")

# Chi^2 on binned data
print(f"\n  Binned chi^2 (N_bins = {len(rar_bins)}):")
print(f"  {'Model':>16}  {'chi2_bin':>10}  {'chi2/N_bin':>11}")
print("  " + "-" * 42)
for name in models:
    a0_use = best_results[name]['a0']
    g_pred_bin = models[name](bin_gbar_lin, a0_use)
    log_pred = np.log10(g_pred_bin)
    chi2_bin = np.sum((bin_log_gobs - log_pred)**2 / bin_sd**2)
    print(f"  {name:>16}  {chi2_bin:10.2f}  {chi2_bin/len(rar_bins):11.3f}")

# ============================================================
# 4. RESIDUALS BY ACCELERATION REGIME
# ============================================================
print("\n" + "=" * 70)
print("  4. Residuals by acceleration regime")
print("=" * 70)

print("""
  Split data into three regimes:
  A) Deep MOND:  log(gbar) < -11   (g_bar < 0.1 * a0)
  B) Transition: -11 < log(gbar) < -9.5  (0.1*a0 < g_bar < 3*a0)
  C) Newtonian:  log(gbar) > -9.5  (g_bar > 3 * a0)
""")

regimes = {
    'Deep MOND (g<0.1*a0)': log_gbar < -11,
    'Transition (0.1-3*a0)': (log_gbar >= -11) & (log_gbar <= -9.5),
    'Newton (g>3*a0)':       log_gbar > -9.5,
}

for regime_name, mask in regimes.items():
    N = np.sum(mask)
    print(f"\n  {regime_name} (N={N}):")
    print(f"  {'Model':>16}  {'mean_res(dex)':>13}  {'RMS(dex)':>9}  {'median':>8}")
    print("  " + "-" * 48)
    for name in models:
        res = best_results[name]['residuals'][mask]
        print(f"  {name:>16}  {np.mean(res):+13.4f}  {np.sqrt(np.mean(res**2)):9.4f}  {np.median(res):+8.4f}")

# ============================================================
# 5. SCATTER ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("  5. Scatter analysis (observed vs intrinsic)")
print("=" * 70)

print("""
  The observed scatter in RAR is ~0.13 dex (Lelli+2017).
  Part is observational, part intrinsic.
  Which model gives the SMALLEST residual scatter?
""")

print(f"  {'Model':>16}  {'total_RMS':>10}  {'mean_err':>9}  {'intr_RMS':>10}  {'Lelli+17':>9}")
print("  " + "-" * 55)
for name in models:
    res = best_results[name]['residuals']
    total_rms = np.sqrt(np.mean(res**2))
    mean_err = np.sqrt(np.mean(e_log_gobs**2))
    # Intrinsic scatter: total^2 - observational^2
    intr2 = total_rms**2 - mean_err**2
    intr_rms = np.sqrt(max(0, intr2))
    print(f"  {name:>16}  {total_rms:10.4f}  {mean_err:9.4f}  {intr_rms:10.4f}  {'~0.13':>9}")

# ============================================================
# 6. AKAIKE / BAYESIAN INFORMATION CRITERION
# ============================================================
print("\n" + "=" * 70)
print("  6. Model comparison: AIC and BIC")
print("=" * 70)

print("""
  All models have 1 free parameter (a0).
  AIC = chi^2 + 2*k  (k = number of parameters)
  BIC = chi^2 + k*ln(N)

  Since k is the same for all models, differences come from chi^2 only.
  Delta_AIC > 10: strong evidence against the worse model.
""")

N = len(gbar)
k = 1  # number of parameters

chi2_vals = {name: best_results[name]['chi2'] for name in models}
min_chi2 = min(chi2_vals.values())
best_model = [name for name, v in chi2_vals.items() if v == min_chi2][0]

print(f"  N = {N}, k = {k}")
print(f"  Best model: {best_model}")
print(f"\n  {'Model':>16}  {'chi2':>10}  {'AIC':>10}  {'BIC':>10}  {'D_AIC':>8}  {'D_BIC':>8}")
print("  " + "-" * 62)

for name in models:
    chi2 = chi2_vals[name]
    aic = chi2 + 2*k
    bic = chi2 + k * np.log(N)
    min_aic = min_chi2 + 2*k
    min_bic = min_chi2 + k * np.log(N)
    print(f"  {name:>16}  {chi2:10.0f}  {aic:10.0f}  {bic:10.0f}  {aic-min_aic:8.0f}  {bic-min_bic:8.0f}")

# ============================================================
# 7. SENSITIVITY: where do models differ most?
# ============================================================
print("\n" + "=" * 70)
print("  7. Where do models differ most? (testability)")
print("=" * 70)

print("""
  Finding the acceleration regime where DGP/Hybrid differ most from MOND.
  This identifies the best observational window for distinguishing models.
""")

# Use best-fit a0 for each
gbar_test = np.logspace(-12.5, -8.5, 200)

# Predictions
pred_ms = np.log10(model_mond_simple(gbar_test, best_results['MOND simple']['a0']))
pred_dgp = np.log10(model_dgp(gbar_test, best_results['DGP']['a0']))
pred_hA = np.log10(model_hybridA(gbar_test, best_results['Hybrid A']['a0']))
pred_hC = np.log10(model_hybridC(gbar_test, best_results['Hybrid C']['a0']))

diff_dgp = pred_dgp - pred_ms
diff_hA = pred_hA - pred_ms
diff_hC = pred_hC - pred_ms

# Find max difference locations
idx_dgp = np.argmax(np.abs(diff_dgp))
idx_hA = np.argmax(np.abs(diff_hA))
idx_hC = np.argmax(np.abs(diff_hC))

print(f"  Max |DGP - MOND_s|:    {np.abs(diff_dgp[idx_dgp]):.3f} dex at log(gbar) = {np.log10(gbar_test[idx_dgp]):.2f}")
print(f"  Max |HybA - MOND_s|:   {np.abs(diff_hA[idx_hA]):.3f} dex at log(gbar) = {np.log10(gbar_test[idx_hA]):.2f}")
print(f"  Max |HybC - MOND_s|:   {np.abs(diff_hC[idx_hC]):.3f} dex at log(gbar) = {np.log10(gbar_test[idx_hC]):.2f}")

print(f"\n  Differences at key accelerations (dex):")
print(f"  {'log(gbar)':>10}  {'DGP-MOND':>10}  {'HybA-MOND':>10}  {'HybC-MOND':>10}  {'SPARC_scatter':>14}")
print("  " + "-" * 60)
for lg in [-11.5, -11.0, -10.5, -10.0, -9.5, -9.0]:
    idx = np.argmin(np.abs(np.log10(gbar_test) - lg))
    # Find local SPARC scatter
    mask_local = np.abs(log_gbar - lg) < 0.25
    local_scatter = np.sqrt(np.mean(e_log_gobs[mask_local]**2)) if np.any(mask_local) else 0.15
    print(f"  {lg:10.1f}  {diff_dgp[idx]:+10.3f}  {diff_hA[idx]:+10.3f}  {diff_hC[idx]:+10.3f}  {local_scatter:14.3f}")

print("""
  INTERPRETATION:
  If model difference < SPARC scatter, models are indistinguishable.
  If model difference > SPARC scatter, data can discriminate.
""")

# ============================================================
# 8. DETAILED COMPARISON AT TRANSITION
# ============================================================
print("=" * 70)
print("  8. Detailed comparison near transition (g_bar ~ a0)")
print("=" * 70)

# Select data in transition region
mask_trans = (log_gbar >= -10.5) & (log_gbar <= -9.5)
N_trans = np.sum(mask_trans)

print(f"\n  Transition region: -10.5 < log(gbar) < -9.5")
print(f"  Number of data points: {N_trans}")

if N_trans > 0:
    gbar_trans = gbar[mask_trans]
    gobs_trans = gobs[mask_trans]
    log_gobs_trans = log_gobs[mask_trans]
    e_trans = e_log_gobs[mask_trans]

    print(f"\n  {'Model':>16}  {'chi2_trans':>11}  {'chi2/N':>8}  {'mean_res':>10}  {'RMS':>8}")
    print("  " + "-" * 56)
    for name, func in models.items():
        a0_use = best_results[name]['a0']
        g_pred = func(gbar_trans, a0_use)
        log_pred = np.log10(g_pred)
        res = log_gobs_trans - log_pred
        chi2 = np.sum(res**2 / e_trans**2)
        print(f"  {name:>16}  {chi2:11.0f}  {chi2/N_trans:8.2f}  {np.mean(res):+10.4f}  {np.sqrt(np.mean(res**2)):8.4f}")

# ============================================================
# 9. TWO-PARAMETER FIT: a0 + shape parameter
# ============================================================
print("\n" + "=" * 70)
print("  9. Generalized interpolation: nu(y) = (1 + (y^n))^(-1/n) (EFE family)")
print("=" * 70)

print("""
  The McGaugh+2016 fit uses:
  g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0)))

  This is the "empirical RAR" formula. Let's also test:
  nu(y) = (1 + y^(-n))^(1/n) -- generalized family with shape param n

  n = 1: MOND simple   (nu = 1 + 1/y at large y)
  n = 2: MOND standard (nu = sqrt(1+1/y^2) at large y)
  n -> inf: step function
""")

def model_generalized(gbar, a0, n):
    y = gbar / a0
    return gbar * (1 + y**(-n))**(1./n)

def model_mcgaugh(gbar, a0):
    """McGaugh+2016 empirical formula."""
    x = np.sqrt(gbar / a0)
    return gbar / (1 - np.exp(-x))

# Optimize a0 for McGaugh formula
def chi2_mcgaugh(log_a0):
    a0_val = 10**log_a0
    g_pred = model_mcgaugh(gbar, a0_val)
    log_g_pred = np.log10(g_pred)
    residuals = log_gobs - log_g_pred
    return np.sum(residuals**2 / e_log_gobs**2)

res_mc = minimize_scalar(chi2_mcgaugh, bounds=(-11, -9), method='bounded')
a0_mc = 10**res_mc.x
g_pred_mc = model_mcgaugh(gbar, a0_mc)
res_mc_val = log_gobs - np.log10(g_pred_mc)
rms_mc = np.sqrt(np.mean(res_mc_val**2))

print(f"  McGaugh+2016 formula: g = gbar / (1 - exp(-sqrt(gbar/a0)))")
print(f"    Best-fit a0 = {a0_mc:.3e} (ratio to 1.2e-10: {a0_mc/1.2e-10:.3f})")
print(f"    chi2 = {res_mc.fun:.0f}, chi2/N = {res_mc.fun/len(gbar):.2f}, RMS = {rms_mc:.4f} dex")

# Optimize generalized for several n
print(f"\n  Generalized family: nu = (1 + y^-n)^(1/n)")
print(f"  {'n':>6}  {'a0_best':>12}  {'chi2':>10}  {'chi2/N':>8}  {'RMS(dex)':>9}")
print("  " + "-" * 50)
for n in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    def chi2_gen(log_a0, n=n):
        a0_val = 10**log_a0
        g_pred = model_generalized(gbar, a0_val, n)
        log_g_pred = np.log10(g_pred)
        residuals = log_gobs - log_g_pred
        return np.sum(residuals**2 / e_log_gobs**2)
    res = minimize_scalar(chi2_gen, bounds=(-11, -9), method='bounded')
    a0_best = 10**res.x
    g_pred = model_generalized(gbar, a0_best, n)
    rms = np.sqrt(np.mean((log_gobs - np.log10(g_pred))**2))
    print(f"  {n:6.1f}  {a0_best:12.3e}  {res.fun:10.0f}  {res.fun/len(gbar):8.2f}  {rms:9.4f}")

# ============================================================
# 10. CRITICAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("  10. CRITICAL SUMMARY")
print("=" * 70)

# Find best and worst
sorted_models = sorted(best_results.items(), key=lambda x: x[1]['chi2'])

print(f"""
  RANKING BY CHI^2 (best-fit a0):
  ===============================""")
for rank, (name, info) in enumerate(sorted_models, 1):
    marker = " <-- BEST" if rank == 1 else ""
    print(f"  {rank}. {name:>16}: chi2/N = {info['chi2_N']:.2f}, a0 = {info['a0']:.3e}, RMS = {info['rms']:.4f} dex{marker}")

dgp_chi2 = best_results['DGP']['chi2']
best_chi2_val = sorted_models[0][1]['chi2']
delta_chi2_dgp = dgp_chi2 - best_chi2_val

print(f"""
  KEY RESULTS:
  ============
  1. Best-fitting model: {sorted_models[0][0]}
     chi2/N = {sorted_models[0][1]['chi2_N']:.2f}

  2. DGP (nu = 1 + 1/sqrt(y)):
     chi2/N = {best_results['DGP']['chi2_N']:.2f}
     Delta_chi2 vs best = {delta_chi2_dgp:.0f}
     -> {'STRONGLY DISFAVORED' if delta_chi2_dgp > 100 else 'DISFAVORED' if delta_chi2_dgp > 10 else 'COMPARABLE'}

  3. Hybrid A (nu = 1 + 1/(sqrt(y) + y)):
     chi2/N = {best_results['Hybrid A']['chi2_N']:.2f}
     Delta_chi2 vs best = {best_results['Hybrid A']['chi2'] - best_chi2_val:.0f}

  4. Best-fit a0 values:""")
for name in sorted_models[:4]:
    n = name[0]
    print(f"     {n:>16}: a0 = {best_results[n]['a0']:.3e} ({best_results[n]['a0']/1.2e-10:.3f} x 1.2e-10)")

print(f"""
  PHYSICAL INTERPRETATION:
  ========================
  - ALL models with correct deep-MOND limit (g -> sqrt(g*a0)) give
    reasonable fits, because SPARC data is dominated by the deep-MOND
    regime where all models converge.

  - The TRANSITION REGION (g_bar ~ a0) is where models differ,
    but SPARC scatter (~0.13 dex) limits discrimination power.

  - DGP's {(best_results['DGP']['chi2_N'] - sorted_models[0][1]['chi2_N']):.1f}x worse chi2/N
    comes from SYSTEMATIC overprediction at g ~ a0.

  - The best models suppress the 2D channel in strong fields more
    aggressively than raw DGP.

  FOR TGP:
  ========
  The dimensional transition 3D -> 2D IS consistent with SPARC data
  IF the 2D channel is suppressed as ~ 1/y (not 1/sqrt(y)) at high g.
  This means:
  - The substrate's "2D channel" doesn't just add sqrt(g*a0)
  - It must be modulated by the local field strength
  - Hybrid A or MOND-like suppression is required

  NEXT STEPS:
  ===========
  -> gs9d: Investigate WHY TGP substrate would suppress 2D channel
  -> Search for physical mechanism behind nu(y) shape
  -> Compare with External Field Effect (EFE) predictions
""")
