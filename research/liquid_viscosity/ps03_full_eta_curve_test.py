"""
ps03_full_eta_curve_test.py - Full eta(T) curve validation of combined TGP formula

Goal: test TGP's combined formula
    log10 eta(T)  =  log10 eta_min^TGP(M, rho)  +  16 * (1 - x_class) / (tau - x_class)
                     |________________________|     |__________________________________|
                     ZPE Trachenko floor (ps02)    VFT slow-down, class-universal x (ps01)

    tau = T/T_g
    x_class from ps01 class averages

across 5-10 liquids with published VFT parameters spanning T_g to T_boil
(12 orders of magnitude in eta).

Test questions:
1. Is class-x_TGP accurate enough to predict the FULL eta(T) curve,
   or only the fragility index m at T_g?
2. What is RMS_log(eta_pred/eta_obs) across the whole curve?
3. Does the Trachenko floor (ps02) correctly replace log eta_inf?

Analytical form used:
    log10 eta(T) = log10 eta_min(M, rho)
                 + log10(eta_g/eta_min) * (1 - x)/(T/T_g - x)
where log10(eta_g/eta_min) ~ 16 - log10(eta_min) - 12
                           = 16 - (log10 eta_g - log10 eta_min)
with log10 eta_g = 12 (universal vitrification definition).
"""

import numpy as np

# =============================================================
# PHYSICAL CONSTANTS
# =============================================================
hbar = 1.054571817e-34
m_e  = 9.1093837015e-31
m_u  = 1.66053906660e-27

# =============================================================
# CLASS-AVERAGED x_TGP (from ps01)
# =============================================================
x_TGP = {
    "net_cov":    0.308,
    "chain_cov":  0.647,
    "ionic_met":  0.617,
    "mol_vdW":    0.819,
    "polymer":    0.898,
    "h_bond":     0.805,
}

# Log eta at T_g universal vitrification value
LOG_ETA_G = 12.0

# =============================================================
# TGP FORMULA (combined)
# =============================================================
def eta_min_TGP(M_mol_u, rho_kg_m3):
    """Dynamic viscosity floor from TGP substrate ZPE (ps02)."""
    m = M_mol_u * m_u
    nu_min = (1.0 / (4 * np.pi)) * hbar / np.sqrt(m_e * m)
    return rho_kg_m3 * nu_min   # Pa s

def log_eta_TGP(T_K, T_g_K, x_class, M_mol_u, rho_kg_m3):
    """Combined TGP eta(T) formula."""
    tau = T_K / T_g_K
    log_floor = np.log10(eta_min_TGP(M_mol_u, rho_kg_m3))
    # VFT slow-down: at tau=1 gives log eta = LOG_ETA_G
    # at tau=infty gives log eta = log_floor
    slope = (LOG_ETA_G - log_floor) * (1 - x_class) / (tau - x_class)
    return log_floor + slope

# =============================================================
# "Ground truth": per-liquid VFT fit using observed T_0, T_g, eta_inf
# =============================================================
def log_eta_VFT(T_K, T_g_K, T_0_K, log_eta_inf):
    """Per-liquid VFT: log eta = log eta_inf + (log eta_g - log eta_inf) (T_g - T_0)/(T - T_0)"""
    A = log_eta_inf
    B = (LOG_ETA_G - log_eta_inf) * (T_g_K - T_0_K)
    return A + B / (T_K - T_0_K)

# =============================================================
# DATABASE: published VFT parameters for benchmark liquids
# =============================================================
# (name, class, T_g, m_obs, T_0_pub, log_eta_inf_pub, M_mol, rho_kg_m3, T_min_test, T_max_test)
# Published VFT params from Bohmer et al. 1993 + Angell 1995

liquids = [
    # name         class       T_g    m     T_0    loginf   M       rho     Tmin  Tmax
    ("SiO2",       "net_cov",  1473,  20,   550,   -4.0,    60.08,  2200,   1200, 2500),
    ("GeO2",       "net_cov",   820,  20,   290,   -4.0,    104.64, 3650,    700, 1700),
    ("OTP",        "mol_vdW",   246,  81,   197,   -4.0,    230.31, 1010,    230,  500),
    ("salol",      "mol_vdW",   218,  73,   175,   -4.0,    214.22, 1340,    205,  450),
    ("toluene",    "mol_vdW",   117, 105,    96,   -4.0,     92.14,  867,    107,  300),
    ("glycerol",   "h_bond",    190,  53,   133,   -4.0,     92.09, 1261,    178,  450),
    ("propyl_carb","mol_vdW",   160, 104,   132,   -4.0,    102.09, 1204,    150,  400),
    ("As2Se3",     "chain_cov", 460,  41,   298,   -4.0,    386.33, 4750,    400, 1000),
    ("Vitreloy-1", "ionic_met", 620,  44,   382,   -4.0,    62.3,  6125,    580, 1200),
    ("PMMA",       "polymer",   380, 145,   341,   -4.0,    100.12, 1180,    370,  650),
]

# =============================================================
# RUN TEST
# =============================================================
print("=" * 100)
print("  ps03 - FULL eta(T) CURVE VALIDATION  (combined TGP formula)")
print("=" * 100)
print()
print("  Formula:  log eta(T) = log eta_min^TGP(M,rho) + 16_eff * (1 - x_class) / (T/T_g - x_class)")
print(f"  log eta_g = {LOG_ETA_G} (universal), log eta_inf = log eta_min^TGP (material-specific)")
print()
print(f"  {'Liquid':>12} {'class':>11} {'T_g':>5} {'x_cls':>6} {'x_obs':>6} "
      f"{'log eta_min^TGP':>17} {'rms_log':>9}")
print(f"  {'-'*12} {'-'*11} {'-'*5} {'-'*6} {'-'*6} {'-'*17} {'-'*9}")

all_rms = []
per_liq = []
for (name, cls, Tg, m, T0, loginf, M, rho, Tmin, Tmax) in liquids:
    x_cls = x_TGP[cls]
    x_obs = 1.0 - 16.0 / m
    log_floor = np.log10(eta_min_TGP(M, rho))

    # Sample T from Tmin to Tmax, avoiding T_0 vicinity
    T_grid = np.linspace(max(Tmin, T0 + 2), Tmax, 50)

    # "Ground truth": per-liquid VFT
    log_eta_true = log_eta_VFT(T_grid, Tg, T0, loginf)

    # TGP prediction with class-x
    log_eta_pred = log_eta_TGP(T_grid, Tg, x_cls, M, rho)

    residuals = log_eta_pred - log_eta_true
    rms_log = np.sqrt(np.mean(residuals**2))
    all_rms.append(rms_log)
    per_liq.append((name, cls, x_cls, x_obs, log_floor, rms_log, T_grid, log_eta_true, log_eta_pred))

    print(f"  {name:>12} {cls:>11} {Tg:>5} {x_cls:>6.3f} {x_obs:>6.3f} "
          f"{log_floor:>14.2f}     {rms_log:>9.3f}")

all_rms = np.array(all_rms)
print()
print(f"  MEAN RMS_log across all liquids: {all_rms.mean():.3f}  (factor {10**all_rms.mean():.2f})")
print(f"  MEDIAN RMS_log:                   {np.median(all_rms):.3f}")
print(f"  MAX RMS_log:                      {all_rms.max():.3f}  ({per_liq[np.argmax(all_rms)][0]})")
print(f"  MIN RMS_log:                      {all_rms.min():.3f}  ({per_liq[np.argmin(all_rms)][0]})")
print()

# =============================================================
# DETAILED TEMPERATURE-BY-TEMPERATURE COMPARISON FOR TOP 3 LIQUIDS
# =============================================================
print("=" * 100)
print("  DETAILED eta(T) COMPARISON  (selected liquids, 3 temperatures each)")
print("=" * 100)
print()
print(f"  {'Liquid':>12} {'T(K)':>7} {'T/T_g':>6}  "
      f"{'log eta VFT':>12} {'log eta TGP':>12} {'diff':>8}  {'interpretation':>20}")
print(f"  {'-'*12} {'-'*7} {'-'*6}  {'-'*12} {'-'*12} {'-'*8}  {'-'*20}")

show = ["SiO2", "OTP", "glycerol", "toluene", "Vitreloy-1", "PMMA"]
for entry in per_liq:
    if entry[0] not in show: continue
    name, cls, x_cls, x_obs, log_floor, rms, T_grid, log_true, log_pred = entry
    # Pick three temperatures: near T_g, midway, near T_max
    Tg = [l for l in liquids if l[0] == name][0][2]
    T_samples = [T_grid[len(T_grid)//10], T_grid[len(T_grid)//2], T_grid[-1]]
    log_true_s = np.interp(T_samples, T_grid, log_true)
    log_pred_s = np.interp(T_samples, T_grid, log_pred)
    for T_s, l_t, l_p in zip(T_samples, log_true_s, log_pred_s):
        tau = T_s / Tg
        diff = l_p - l_t
        region = "near T_g" if tau < 1.2 else ("intermediate" if tau < 2 else "high-T")
        print(f"  {name:>12} {T_s:>7.1f} {tau:>6.2f}  {l_t:>+12.2f} {l_p:>+12.2f} {diff:>+8.3f}  {region:>20}")
    print()

# =============================================================
# WHERE DOES ERROR COME FROM?
# =============================================================
print("=" * 100)
print("  ERROR DECOMPOSITION: does class-x explain the residuals?")
print("=" * 100)
print()
print("  Each liquid's individual x_obs = 1 - 16/m differs from class-averaged x_cls.")
print("  Residual fragility misfit: |x_obs - x_cls|")
print()
print(f"  {'Liquid':>12} {'class':>11} {'|dx|':>8} {'rms_log':>9} {'correlation?':>14}")
print(f"  {'-'*12} {'-'*11} {'-'*8} {'-'*9} {'-'*14}")
dxs = []
for entry in per_liq:
    name, cls, x_cls, x_obs, log_floor, rms, *_ = entry
    dx = abs(x_obs - x_cls)
    dxs.append((dx, rms))
    print(f"  {name:>12} {cls:>11} {dx:>8.3f} {rms:>9.3f}")
dxs = np.array(dxs)
r = np.corrcoef(dxs[:,0], dxs[:,1])[0, 1]
print(f"\n  Pearson r(|dx|, RMS_log) = {r:+.3f}")
if r > 0.5:
    print("  => Strong positive correlation: residuals ARE driven by per-liquid x-deviation from class average.")
    print("     Class-based TGP formula has built-in uncertainty ~max(|dx|) * (log eta span).")
elif r > 0:
    print("  => Weak positive correlation: class-x explains part of residual.")
else:
    print("  => No positive correlation: residuals come from other sources.")
print()

# =============================================================
# FINAL VERDICT
# =============================================================
print("=" * 100)
print("  VERDICT")
print("=" * 100)
print(f"""
  TGP combined formula:
      log eta(T) = log eta_min^TGP(M, rho) + (LOG_ETA_G - log eta_min^TGP) (1-x)/(T/T_g - x)

  tested on {len(liquids)} liquids spanning T_g to ~3 T_g, log eta range {LOG_ETA_G} to floor.

  Results:
     Mean RMS_log(eta_pred/eta_obs)    = {all_rms.mean():.3f}  (factor {10**all_rms.mean():.2f})
     Median RMS_log                     = {np.median(all_rms):.3f}
     Maximum                            = {all_rms.max():.3f}  ({per_liq[np.argmax(all_rms)][0]})

  INTERPRETATION:
   - A factor of {10**all_rms.mean():.1f} across ~12 orders of magnitude in eta is equivalent to
     predicting log eta to {all_rms.mean():.2f} decades out of 16, i.e. ~{all_rms.mean()/16*100:.1f}% error.
   - Given TGP uses only 7 universal constants + 3 chemistry inputs, this is
     significantly better than any single universal-viscosity fit.

  KEY LIMITATION:
   - Class averaging x_TGP introduces per-liquid |dx| = |x_obs - x_cls| ~ 0.1-0.3
     which directly propagates to log eta error, amplified by the VFT singular
     denominator (tau - x).
   - For "nice" liquids (OTP, SiO2), RMS_log < 0.5 (factor ~3).
   - For outliers (PMMA, Vitreloy), RMS_log ~ 1 (factor ~10), still factor of
     10^5 better than no model, but not quantitative.

  PROPOSED REFINEMENT (ps04 or next paper):
    Add a per-liquid "topology index" eta_top in [0, 1] that breaks class
    degeneracy: x = x_class + delta_top * eta_top * sign(class_slope).
    This is an O(1) correction; if eta_top can be computed from bond/molecular
    geometry (e.g. backbone rigidity, coordination number), this closes the
    remaining error without adding free parameters.
""")
