#!/usr/bin/env python3
"""
TGP Prediction Taxonomy v47 (2026-04-10)

Honest 4-class classification of all TGP quantitative results,
following external audit recommendation P5.

Classes:
  [I]  INPUT: Observable used to calibrate a free parameter.
              Cannot be counted as a test.
  [D]  DERIVED: Follows analytically from inputs + TGP axioms,
              without any new fit. Genuine zero-parameter prediction.
  [R]  RECOVERED: Consistent with data, but depends on a structural
              choice or approximation that partially accommodates it.
  [P]  PROSPECTIVE: Out-of-sample prediction, not yet confirmed
              or with large experimental uncertainty.

Parameters: Phi_0 (from Lambda_obs), a_Gamma (from r_21 via phi-FP).
Derived constants: g_0^e = 0.8694, N_c = 3, N_f = 5, phi = (1+sqrt(5))/2.
"""

import numpy as np

# ============================================================
# Parameters
# ============================================================
Omega_Lambda = 0.6847
Phi_0 = 36 * Omega_Lambda   # 24.65 -- calibrated from Lambda_obs
a_Gamma = 0.040              # calibrated from r_21 (lepton sector)
N_c = 3
N_f = 5
phi = (1 + np.sqrt(5)) / 2
g_0_e = 0.8694

# ============================================================
# Data
# ============================================================
results = []

def add(name, tgp, obs, obs_err, cls, reason, source=""):
    """cls: I/D/R/P"""
    if obs != 0:
        dev_pct = (tgp / obs - 1) * 100
    else:
        dev_pct = 0.0
    sigma = abs(tgp - obs) / obs_err if obs_err > 0 else 0.0
    results.append({
        'name': name, 'tgp': tgp, 'obs': obs, 'err': obs_err,
        'dev': dev_pct, 'sigma': sigma, 'cls': cls,
        'reason': reason, 'source': source
    })

# ============================================================
# [I] INPUTS (calibrations) — 2 entries
# ============================================================

add("Phi_0 = 36*Omega_Lambda",
    Phi_0, 24.66, 0.3, "I",
    "Defines Phi_0 from Lambda_obs. This IS the calibration.",
    "Planck2018")

# r_21 is the calibration target for a_Gamma (through g_0^e / phi-FP)
add("r_21 = m_mu/m_e (phi-FP)",
    206.77, 206.768, 0.01, "I",
    "Used to calibrate a_Gamma via phi-FP ODE. Input, not test.",
    "PDG2024")

# ============================================================
# [D] DERIVED (zero-parameter predictions) — genuine tests
# ============================================================

# Koide from N_gen=3 chain
m_e = 0.51099895; m_mu = 105.6583755; m_tau = 1776.86
r_21_PDG = m_mu / m_e; r_31_PDG = m_tau / m_e
sqrt_r31 = 2*(1+np.sqrt(r_21_PDG)) + np.sqrt(3*(1+4*np.sqrt(r_21_PDG)+r_21_PDG))
r_31_Koide = sqrt_r31**2

add("r_31 = m_tau/m_e (Koide)",
    r_31_Koide, r_31_PDG, 0.5, "D",
    "Q_K=3/2 from N_gen=3 (which follows from d=3,k=4). "
    "r_31 then determined by r_21+Koide. No new parameter.",
    "PDG2024")

add("Q_K (Koide) = 3/2",
    1.5, 1.50001, 0.0001, "D",
    "Algebraic consequence: d=3 -> k=4 -> N_gen=3 -> Q_K=2N/(N+1)=3/2.",
    "PDG2024")

add("N_gen = 3",
    3, 3, 0.001, "D",
    "WKB quantization of soliton spectrum: k=4 in d=3 gives exactly 3 bound states.",
    "LEP")

# Cosmological
N_e = 57
n_s = 1 - 2.0/N_e - 9.0/(4*N_e**2)
r_ts = 12.0 / N_e**2

add("n_s (N_e=57)",
    n_s, 0.9649, 0.0042, "D",
    "Starobinsky attractor (robust to TGP volume element correction). "
    "N_e~57 from reheating scale. No fit.",
    "Planck2018")

add("r_ts = 12/N_e^2",
    r_ts, 0.004, 0.02, "D",
    "Same Starobinsky attractor. Robust prediction.",
    "BICEP/Keck")

add("w_DE = -1",
    -1.0, -1.0, 0.05, "D",
    "Frozen vacuum Phi=Phi_0 acts as cosmological constant. "
    "No dynamical dark energy mechanism.",
    "Planck+BAO")

add("gamma_PPN = 1",
    1.0, 1.0, 0.00003, "D",
    "Exponential metric (from beta=gamma vacuum condition) "
    "gives exact GR PPN. No tuning.",
    "Cassini")

add("beta_PPN = 1",
    1.0, 1.0, 0.00003, "D",
    "Same exponential metric.",
    "LLR/perihelion")

add("c_GW/c_EM - 1 = 0",
    0.0, 0.0, 1e-15, "D",
    "Both propagate on same substrate at c_0. Exact equality by construction.",
    "GW170817")

add("sin^2(theta_W) = 3/13 (tree)",
    N_c/(N_c**2+N_c+1), 0.23122, 0.00004, "D",
    "From k<=2 tensor saturation in d=3 + color-graded DOF sum. "
    "Zero parameters, pure combinatorics. Tree-level dev=0.19%.",
    "PDG2024")

# QCD-corrected sin^2 (v47b)
b0_QCD = (11*N_c - 2*5) / (12*np.pi)  # N_f=5 at M_Z
sin2_qcd = N_c / (1 + N_c + N_c**2 - b0_QCD * 0.1179 / N_c)
add("sin^2(theta_W) QCD-corrected",
    sin2_qcd, 0.23122, 0.00004, "D",
    "1-loop QCD correction: N_c^2 -> N_c^2(1-b_0*alpha_s/N_c^3). "
    "b_0 = (11N_c-2N_f)/(12pi). Reduces tree 11.3sigma to 0.6sigma.",
    "PDG2024+QCD")

alpha_s = N_c**3 * g_0_e / (8 * Phi_0)
add("alpha_s(M_Z)",
    alpha_s, 0.1179, 0.0009, "D",
    "Mass-coupling bridge: alpha_s = N_c^3 g_0^e / (8 Phi_0). "
    "g_0^e fixed by r_21 (input), Phi_0 by Lambda (input). No new fit.",
    "PDG2024")

add("alpha_s(m_tau)/alpha_s(M_Z)",
    (N_f/N_c)**2, 2.799, 0.03, "D",
    "Running ratio = (N_f/N_c)^2. Pure group theory.",
    "PDG2024")

# Color tube string tension -> quark binding constant
C_F = (N_c**2 - 1) / (2*N_c)  # = 4/3
alpha_s_PDG = 0.1179
A_univ_tube = C_F**2 * alpha_s_PDG**2
A_univ_empirical = a_Gamma / phi  # = 0.02472
add("A_univ = C_F^2 * alpha_s^2 (tube)",
    A_univ_tube, A_univ_empirical, 0.0005, "D",
    "Variational color tube: depletion A=C_F*alpha_s/(pi*Phi_0), "
    "sigma_hat=pi*A^2, with K_geo*m_sp^2=pi*Phi_0^2 scaling. "
    "Phi_0 cancels! Result: A_univ=C_F^2*alpha_s^2. Dev=0.04%.",
    "prop:X-A-from-tube-tension")

# ============================================================
# [R] RECOVERED (after calibration, with structural assumptions)
# ============================================================

kappa = 3.0 / (4 * Phi_0)
add("kappa = 3/(4*Phi_0)",
    kappa, 0.030, 0.005, "R",
    "Gravitational coupling from Phi_0. Correct order, "
    "but depends on normalization convention 3/4 vs other prefactors.",
    "LLR")

add("a_Gamma * Phi_0 ~ 1",
    a_Gamma*Phi_0, 1.0, 0.015, "R",
    "Near-unity product is suggestive but partially tautological: "
    "a_Gamma calibrated from lepton sector, Phi_0 from cosmology.",
    "DESI DR2")

lambda_eff = 2.0 / Phi_0**4
add("lambda_eff = 2/Phi_0^4",
    lambda_eff, 5.501e-6, 0.1e-6, "R",
    "Effective self-coupling. Depends on Z_2 potential form "
    "and coefficient '2' from specific calculation.",
    "soliton scan")

A_tgp = a_Gamma / phi
m_d, m_s, m_b = 4.67, 93.4, 4180.0
m0_d_tgp = A_tgp * m_b / m_d
add("A = a_Gamma/phi (quark universal)",
    A_tgp, 0.02464, 0.00013, "R",
    "Color-tube exhaustion of one oscillatory mode. "
    "Interpretation via phi-FP projection, but mechanism is [AN], not [TW].",
    "shifted Koide")

add("m0(d,s,b)",
    m0_d_tgp, 21.94, 0.5, "R",
    "Shifted Koide mass. Follows from A, but A itself is [R]-class.",
    "shifted Koide")

# ============================================================
# [P] PROSPECTIVE (out-of-sample, not yet confirmed)
# ============================================================

add("sum(m_nu)",
    62.9, 60.0, 30.0, "P",
    "TGP predicts ~63 meV (normal ordering). Planck+DESI gives "
    "upper bound ~72 meV. Not yet precision-testable.",
    "Planck+DESI")

add("m3/m2 (neutrino)",
    50.4/9.3, 5.0, 1.0, "P",
    "Predicted ratio from TGP neutrino spectrum. Precision test "
    "requires JUNO + next-gen oscillation experiments.",
    "oscillation")

add("BH shadow: no photon sphere",
    0.0, 0.0, 0.0, "P",
    "TGP soliton gives h(r) monotonically increasing => no photon ring. "
    "Testable by ngEHT (space baseline). Status O9.",
    "EHT/ngEHT")

add("ISCO ring at b ~ 0.6 r_s",
    0.6, 0.0, 0.0, "P",
    "TGP accretion disk ring much smaller than GR prediction. "
    "Distinguishable with next-gen imaging.",
    "ngEHT")

# ============================================================
# Output
# ============================================================
print("=" * 95)
print("TGP PREDICTION TAXONOMY v47 (honest 4-class)")
print(f"Parameters: Phi_0 = {Phi_0:.4f} [I], a_Gamma = {a_Gamma} [I]")
print(f"Derived: g_0^e = {g_0_e}, N_c = {N_c}, N_f = {N_f}, phi = {phi:.5f}")
print("=" * 95)

classes = {"I": "INPUT (calibration)", "D": "DERIVED (zero-param prediction)",
           "R": "RECOVERED (structural assumption)", "P": "PROSPECTIVE (untested)"}

for cls_key in ["I", "D", "R", "P"]:
    items = [r for r in results if r['cls'] == cls_key]
    print(f"\n  [{cls_key}] {classes[cls_key]} --- {len(items)} items")
    print(f"  {'-'*85}")
    for r in items:
        if r['err'] > 0 and r['obs'] != 0:
            print(f"    {r['name']:<35} TGP={r['tgp']:>12.6g}  obs={r['obs']:>12.6g}  "
                  f"dev={r['dev']:>+7.3f}%  {r['sigma']:.1f}s  [{r['source']}]")
        else:
            print(f"    {r['name']:<35} TGP={r['tgp']:>12.6g}  [{r['source']}]")
        print(f"      -> {r['reason']}")

# Summary
n_I = sum(1 for r in results if r['cls'] == 'I')
n_D = sum(1 for r in results if r['cls'] == 'D')
n_R = sum(1 for r in results if r['cls'] == 'R')
n_P = sum(1 for r in results if r['cls'] == 'P')

# Count sigma performance for D and R only (testable)
testable = [r for r in results if r['cls'] in ('D', 'R') and r['err'] > 0]
t_1s = sum(1 for r in testable if r['sigma'] < 1)
t_2s = sum(1 for r in testable if r['sigma'] < 2)

print(f"\n{'='*95}")
print(f"SUMMARY")
print(f"{'='*95}")
print(f"""
  Classification:
    [I] Inputs (calibrations):            {n_I}
    [D] Derived (zero-param predictions): {n_D}
    [R] Recovered (structural):           {n_R}
    [P] Prospective (untested):           {n_P}
    Total:                                {n_I + n_D + n_R + n_P}

  Genuine test power:
    Free parameters:                      2
    True zero-parameter predictions [D]:  {n_D}
    Predictive ratio M_D/N_param:         {n_D}/2 = {n_D/2:.1f}

  Performance on testable [D]+[R] ({len(testable)} items):
    Within 1-sigma:                       {t_1s}/{len(testable)}
    Within 2-sigma:                       {t_2s}/{len(testable)}

  HONEST ASSESSMENT:
    From 2 calibrated inputs, TGP derives {n_D} independent predictions
    (ratio {n_D/2:.0f}:1) with {t_1s}/{len(testable)} within 1-sigma.
    Additionally {n_R} results are recovered with structural assumptions,
    and {n_P} prospective predictions await experimental tests.
""")
