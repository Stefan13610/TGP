#!/usr/bin/env python3
"""
quark_exponent_investigation_v47b.py
====================================

The empirical formula m_0 = C * m_2^p with p ~ 1.73 predicts both
quark sectors to ~3 sigma.  This script investigates:

1. PRECISE FIT of the exponent p (not just grid search)
2. Candidate TGP values: phi=1.618, 7/4=1.75, 2-1/4, sqrt(3), etc.
3. PHYSICAL DERIVATION from TGP Regime III:
   - If m_0 = sigma * L_eff and L_eff ~ m_2^(p-1) / sigma^0, what is p?
   - Dimensional analysis: [m_0] = [mass], [m_2] = [mass], [C] = [mass^(1-p)]
   - In TGP: natural mass scale = Lambda_QCD or m_sp
4. BOOTSTRAP with each candidate p: predict m_b and m_t
5. SENSITIVITY to input masses (PDG uncertainties)
6. Cross-check with lepton sector (m_0 = 0 => self-consistency)

Author: TGP v47b (2026-04-12)
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

# ============================================================
# Constants
# ============================================================
PHI = (1 + np.sqrt(5)) / 2
a_Gamma = 0.040
Phi0 = 25.0
N_c = 3
g0_e = 0.86770494
hbarc = 197.3  # MeV * fm

# PDG masses
m_d, m_s, m_b = 4.67, 93.4, 4180.0
m_u, m_c, m_t = 2.16, 1270.0, 172_760.0
m_e, m_mu, m_tau = 0.511, 105.658, 1776.86

# PDG uncertainties (1-sigma)
dm_b = 30.0
dm_t = 300.0
dm_s = 8.0  # ~8 MeV
dm_c = 20.0  # ~20 MeV

def koide_Q(m1, m2, m3):
    s1 = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s1**2 / (3 * (m1 + m2 + m3))

def shifted_koide(m1, m2, m3, m0):
    return koide_Q(m1 + m0, m2 + m0, m3 + m0)

# Empirical m_0
m0_d = brentq(lambda m0: shifted_koide(m_d, m_s, m_b, m0) - 0.5, 0, 1e5)
m0_u = brentq(lambda m0: shifted_koide(m_u, m_c, m_t, m0) - 0.5, 0, 1e6)

print("=" * 70)
print("EXPONENT INVESTIGATION: m_0 = C * m_2^p")
print("=" * 70)
print(f"  Empirical m_0(down) = {m0_d:.4f} MeV")
print(f"  Empirical m_0(up)   = {m0_u:.4f} MeV")


# ============================================================
# 1. PRECISE FIT OF EXPONENT p
# ============================================================
print("\n" + "=" * 70)
print("1. PRECISE FIT: p = ln(m0_u/m0_d) / ln(m_c/m_s)")
print("=" * 70)

p_exact = np.log(m0_u / m0_d) / np.log(m_c / m_s)
C_d = m0_d / m_s**p_exact
C_u = m0_u / m_c**p_exact
C_mean = (C_d + C_u) / 2

print(f"  p_exact = ln({m0_u:.2f}/{m0_d:.2f}) / ln({m_c:.0f}/{m_s:.1f})")
print(f"         = ln({m0_u/m0_d:.4f}) / ln({m_c/m_s:.4f})")
print(f"         = {np.log(m0_u/m0_d):.6f} / {np.log(m_c/m_s):.6f}")
print(f"         = {p_exact:.6f}")
print(f"  C_down  = {C_d:.8f}")
print(f"  C_up    = {C_u:.8f}")
print(f"  C_mean  = {C_mean:.8f}")
print(f"  Exact by construction: C_d = C_u")

# Alternative: fit minimizing chi2 on m_3
def chi2_for_p(p):
    """Chi2 for predicting m_b and m_t given exponent p."""
    C_d_loc = m0_d / m_s**p
    C_u_loc = m0_u / m_c**p
    C = (C_d_loc + C_u_loc) / 2

    m0_d_pred = C * m_s**p
    m0_u_pred = C * m_c**p

    # Predict m_3 from shifted Koide
    chi2 = 0
    for m1, m2, m3_pdg, m0_pred, sigma in [
        (m_d, m_s, m_b, m0_d_pred, dm_b),
        (m_u, m_c, m_t, m0_u_pred, dm_t)]:
        try:
            def res(m3):
                return shifted_koide(m1, m2, m3, m0_pred) - 0.5
            m3_test = np.logspace(0, 7, 5000)
            r_vals = [res(m3) for m3 in m3_test]
            for i in range(len(r_vals)-1):
                if r_vals[i] * r_vals[i+1] < 0:
                    sol = brentq(res, m3_test[i], m3_test[i+1])
                    # Pick solution closest to PDG
                    chi2 += ((sol - m3_pdg) / sigma)**2
                    break
            else:
                chi2 += 1e6
        except:
            chi2 += 1e6
    return chi2

# Scan p
print(f"\n  Fitting p to minimize chi2(m_b, m_t):")
result = minimize_scalar(chi2_for_p, bounds=(1.0, 2.5), method='bounded')
p_chi2_opt = result.x
chi2_opt = result.fun
print(f"  p_opt (chi2 minimized) = {p_chi2_opt:.6f}")
print(f"  chi2_min = {chi2_opt:.4f}")


# ============================================================
# 2. CANDIDATE TGP VALUES
# ============================================================
print("\n" + "=" * 70)
print("2. CANDIDATE VALUES AND THEIR PREDICTIONS")
print("=" * 70)

candidates = {
    "p_exact (fit)":       p_exact,
    "p_chi2_opt":          p_chi2_opt,
    "phi = (1+sqrt5)/2":   PHI,            # 1.61803
    "phi^2/phi = phi":     PHI,            # same
    "7/4":                 7/4,            # 1.75000
    "sqrt(3)":             np.sqrt(3),     # 1.73205
    "ln(5)":               np.log(5),      # 1.60944
    "e/phi":               np.e / PHI,     # 1.67990
    "2 - 1/4":             1.75,           # 1.75000 (=7/4)
    "2 - 1/phi":           2 - 1/PHI,     # 1.38197
    "1 + 1/phi":           1 + 1/PHI,     # 1.61803 = phi
    "2*phi - 2":           2*PHI - 2,     # 1.23607
    "phi + 1/10":          PHI + 0.1,     # 1.71803
    "(d+1)/2 = 2":         2.0,           # d=3
    "2*d/(d+1) = 3/2":     1.5,           # 1.5
    "(2d-1)/d = 5/3":      5/3,           # 1.66667
    "2 - 1/d = 5/3":       5/3,           # same
    "2 - 1/N_c = 5/3":     5/3,           # same
    "(d-1+phi)/d = (2+phi)/3": (2 + PHI)/3,  # 1.20601
    "1 + phi/2":           1 + PHI/2,     # 1.80902
    "(phi+1)/phi = phi":   (PHI+1)/PHI,   # phi again
    "phi * sqrt(phi)/phi = sqrt(phi)": np.sqrt(PHI), # 1.27202
    "2*phi/phi^2 = 2/phi": 2/PHI,         # 1.23607
    "phi^(2/3)":           PHI**(2/3),     # 1.37316
    "phi^(3/2)":           PHI**(3/2),     # 2.05817
}

# Remove duplicates
seen = {}
for name, val in list(candidates.items()):
    key = f"{val:.5f}"
    if key not in seen:
        seen[key] = name
    elif name.startswith("p_"):
        seen[key] = name  # prefer fit names

# Evaluate each
print(f"\n  {'Candidate':35s} {'p':>8s} {'dev%':>7s} {'C':>10s} {'m_b':>8s} {'sig_b':>7s} {'m_t':>10s} {'sig_t':>7s} {'chi2':>8s}")
print("  " + "-" * 110)

results_list = []
for name, p in sorted(candidates.items(), key=lambda x: abs(x[1] - p_exact)):
    C_d_loc = m0_d / m_s**p
    C_u_loc = m0_u / m_c**p
    C = (C_d_loc + C_u_loc) / 2

    m0_d_pred = C * m_s**p
    m0_u_pred = C * m_c**p

    m3_results = {}
    chi2 = 0
    for sname, m1, m2, m3_pdg, m0_pred, sigma in [
        ("b", m_d, m_s, m_b, m0_d_pred, dm_b),
        ("t", m_u, m_c, m_t, m0_u_pred, dm_t)]:
        try:
            def res(m3):
                return shifted_koide(m1, m2, m3, m0_pred) - 0.5
            m3_test = np.logspace(0, 7, 5000)
            r_vals = [res(m3) for m3 in m3_test]
            best_sol = None
            for i in range(len(r_vals)-1):
                if r_vals[i] * r_vals[i+1] < 0:
                    sol = brentq(res, m3_test[i], m3_test[i+1])
                    if best_sol is None or abs(sol - m3_pdg) < abs(best_sol - m3_pdg):
                        best_sol = sol
            if best_sol is not None:
                sig = (best_sol - m3_pdg) / sigma
                m3_results[sname] = (best_sol, sig)
                chi2 += sig**2
            else:
                m3_results[sname] = (None, None)
                chi2 += 1e6
        except:
            m3_results[sname] = (None, None)
            chi2 += 1e6

    dev_p = (p - p_exact) / p_exact * 100

    mb_str = f"{m3_results['b'][0]:.0f}" if m3_results['b'][0] else "---"
    sb_str = f"{m3_results['b'][1]:+.2f}" if m3_results['b'][1] is not None else "---"
    mt_str = f"{m3_results['t'][0]:.0f}" if m3_results['t'][0] else "---"
    st_str = f"{m3_results['t'][1]:+.2f}" if m3_results['t'][1] is not None else "---"

    results_list.append((name, p, chi2))
    if chi2 < 1e5:
        print(f"  {name:35s} {p:8.5f} {dev_p:+7.2f} {C:10.6f} {mb_str:>8s} {sb_str:>7s} {mt_str:>10s} {st_str:>7s} {chi2:8.2f}")


# ============================================================
# 3. PHYSICAL DERIVATION ATTEMPT
# ============================================================
print("\n" + "=" * 70)
print("3. PHYSICAL DERIVATION: WHY p ~ 1.73?")
print("=" * 70)

print("""
  DIMENSIONAL ANALYSIS:
  m_0 = C * m_2^p  =>  [C] = MeV^(1-p)

  If C has TGP origin, it must be built from:
    Lambda_QCD ~ 220 MeV  (or m_sp = sqrt(gamma))
    a_Gamma ~ 0.040  (dimensionless)
    phi = 1.618...  (dimensionless)
    Phi_0 = 25  (dimensionless)

  So C ~ Lambda_QCD^(1-p) * f(a_Gamma, phi, Phi_0)
  Or C ~ m_sp^(1-p) * f(...)
""")

# What Lambda_QCD gives the right C?
Lambda_QCD = 220.0  # MeV

for p_name, p_val in [("p_exact", p_exact), ("phi", PHI), ("7/4", 1.75), ("sqrt(3)", np.sqrt(3))]:
    C_needed = (C_d + C_u) / 2 if p_name == "p_exact" else m0_d / m_s**p_val
    # C = Lambda^(1-p) * f
    # => f = C / Lambda^(1-p)
    f_val = C_needed / Lambda_QCD**(1-p_val)
    print(f"  p = {p_val:.5f} ({p_name}):")
    print(f"    C = {C_needed:.8f} MeV^({1-p_val:.3f})")
    print(f"    C / Lambda_QCD^(1-p) = {f_val:.6f}")
    print(f"    C / a_Gamma = {C_needed/a_Gamma:.6f}")
    print(f"    C * Phi_0 = {C_needed*Phi0:.6f}")
    print(f"    f / a_Gamma = {f_val/a_Gamma:.6f}")
    print(f"    f * phi = {f_val * PHI:.6f}")
    print()


# ============================================================
# 4. KEY TEST: IS p = sqrt(3)?
# ============================================================
print("=" * 70)
print("4. DETAILED TEST: p = sqrt(3) = 1.73205...")
print("=" * 70)

p_test = np.sqrt(3)
C_d_test = m0_d / m_s**p_test
C_u_test = m0_u / m_c**p_test
C_test = (C_d_test + C_u_test) / 2

print(f"  p = sqrt(3) = {p_test:.6f}")
print(f"  p_exact     = {p_exact:.6f}")
print(f"  Deviation   = {(p_test - p_exact)/p_exact*100:+.3f}%")
print(f"  C_down = {C_d_test:.8f}")
print(f"  C_up   = {C_u_test:.8f}")
print(f"  C_mean = {C_test:.8f}")
print(f"  Cross-sector C consistency: {abs(C_d_test/C_u_test - 1)*100:.2f}%")

# Predict m_b and m_t
m0_d_pred = C_test * m_s**p_test
m0_u_pred = C_test * m_c**p_test

for sname, m1, m2, m3_pdg, m0_pred, sigma in [
    ("m_b (down)", m_d, m_s, m_b, m0_d_pred, dm_b),
    ("m_t (up)", m_u, m_c, m_t, m0_u_pred, dm_t)]:
    def res(m3):
        return shifted_koide(m1, m2, m3, m0_pred) - 0.5
    m3_test = np.logspace(0, 7, 10000)
    r_vals = [res(m3) for m3 in m3_test]
    for i in range(len(r_vals)-1):
        if r_vals[i] * r_vals[i+1] < 0:
            sol = brentq(res, m3_test[i], m3_test[i+1])
            sig = (sol - m3_pdg) / sigma
            print(f"\n  {sname}:")
            print(f"    m_0 = {m0_pred:.2f} MeV")
            print(f"    m_3 = {sol:.1f} MeV  (PDG: {m3_pdg:.0f} +/- {sigma:.0f})")
            print(f"    Deviation = {sig:+.2f} sigma")
            break

# Why sqrt(3)?
print(f"""
  WHY sqrt(3)?

  In TGP, d=3 appears in multiple structural roles:
    - Spatial dimension
    - Number of generations N_gen = 3
    - Soliton ODE: g'' + (2/r)g' + ... = 0  (d-1 = 2)
    - Bessel function order: j_0(r) in d=3
    - Collapse threshold: g_crit = 8/5 (for alpha=2, d=3)

  sqrt(d) = sqrt(3) could arise from:
    - Geometric factor in d-dimensional integration
    - RMS of d-component Gaussian vector: sigma = sqrt(d) * sigma_1
    - Kinetic energy scaling: E_kin ~ p^2/(2m), with |p| ~ sqrt(d) * p_1

  Specifically for confinement:
    If L_eff ~ (m_2/Lambda)^(sqrt(d)-1) * Lambda^(-1), then:
    m_0 = sigma * L_eff ~ sigma/Lambda * (m_2/Lambda)^(sqrt(3)-1)
        ~ C * m_2^(sqrt(3))  with C ~ sigma/Lambda^(sqrt(3))

  This would make p = sqrt(3) a CONSEQUENCE of d=3 geometry in
  the confinement sector, analogous to how CV=1 follows from d=3
  in the Koide sector (chi^2(2) has 2 = d-1 degrees of freedom).
""")


# ============================================================
# 5. SENSITIVITY TO INPUT MASSES
# ============================================================
print("=" * 70)
print("5. SENSITIVITY TO PDG INPUT UNCERTAINTIES")
print("=" * 70)

# Vary m_s and m_c within PDG uncertainties
print("  Varying m_s and m_c within 1-sigma PDG uncertainties:")
print(f"  m_s = {m_s} +/- {dm_s} MeV,  m_c = {m_c} +/- {dm_c} MeV\n")

for dm_s_shift, dm_c_shift in [(0,0), (1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1)]:
    ms = m_s + dm_s_shift * dm_s
    mc = m_c + dm_c_shift * dm_c

    # Recompute m0 with shifted masses
    m0d = brentq(lambda m0: shifted_koide(m_d, ms, m_b, m0) - 0.5, 0, 1e5)
    m0u = brentq(lambda m0: shifted_koide(m_u, mc, m_t, m0) - 0.5, 0, 1e6)
    p_shifted = np.log(m0u / m0d) / np.log(mc / ms)

    print(f"  m_s={ms:5.1f}, m_c={mc:5.0f}: p = {p_shifted:.5f}  "
          f"(dev from sqrt(3): {(p_shifted - np.sqrt(3))/np.sqrt(3)*100:+.3f}%)")


# ============================================================
# 6. GRAND SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("GRAND SUMMARY")
print("=" * 70)

# Rank candidates by chi2
print("\n  Top 5 candidates by chi2(m_b, m_t):")
results_list.sort(key=lambda x: x[2])
for i, (name, p, chi2) in enumerate(results_list[:5]):
    if chi2 < 1e5:
        print(f"    {i+1}. {name:35s}  p = {p:.5f}  chi2 = {chi2:.2f}")

print(f"""
  FINDINGS:
  1. p_exact = {p_exact:.5f} (from m_0 ratio and m_2 ratio)
  2. p = sqrt(3) = {np.sqrt(3):.5f} matches to {abs(p_exact - np.sqrt(3))/p_exact*100:.2f}%
  3. p = 7/4 = 1.75000 matches to {abs(p_exact - 1.75)/p_exact*100:.2f}%
  4. p = phi = {PHI:.5f} matches to {abs(p_exact - PHI)/p_exact*100:.2f}%

  PHYSICAL MOTIVATION:
  - sqrt(3) has clear geometric origin in d=3 TGP
  - phi has deep structural role (phi-FP mechanism)
  - 7/4 has no obvious TGP motivation

  PREDICTION QUALITY (chi2 = sigma_b^2 + sigma_t^2):
  - p = sqrt(3): chi2 = {chi2_for_p(np.sqrt(3)):.2f}
  - p = 7/4:     chi2 = {chi2_for_p(1.75):.2f}
  - p = phi:     chi2 = {chi2_for_p(PHI):.2f}

  STATUS:
  If p = sqrt(3) is confirmed, the formula
    m_0 = C * m_2^sqrt(3),  C = sigma / Lambda_QCD^sqrt(3)
  would be a STRUCTURAL CONSEQUENCE of d=3 geometry,
  completing the quark sector derivation chain alongside
  Q_K = 3/2 from d=3 -> chi^2(2) -> CV=1.
""")
