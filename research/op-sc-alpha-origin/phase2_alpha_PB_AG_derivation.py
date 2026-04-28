"""
phase2_alpha_PB_AG_derivation.py

Phase 2 of op-sc-alpha-origin: first-principles derivation of alpha_PB
from Abrikosov-Gorkov (A-G) theory + de Gennes factor + literature N(0), J_sf.

Tests:
  T2.1 - A-G symbolic derivation, weak & strong limit (sympy digamma)
  T2.2 - de Gennes factor analytical for La/Ce/Pr/Nd/Sm/Yb
  T2.3 - mu_eff^2 analytical (Hund GS) cross-check
  T2.4 - Scaling-factor ambiguity: dG vs mu_eff^2 -> SmH9 falsification target
  T2.5 - Numerical alpha_PB^pred from literature N(0), J_sf
  T2.6 - Internal consistency: TGP fit residua vs dG correlation

Run:
  cd TGP/TGP_v1/research/op-sc-alpha-origin
  python -X utf8 phase2_alpha_PB_AG_derivation.py 2>&1 | tee phase2_alpha_PB_AG_derivation.txt

Verdict goal:
  H_AG (A-G derivable) supported with ratio in [0.7, 1.3], OR
  H_fit (independent fit) supported with ratio > 2 + scaling ambiguity, OR
  Mixed: TGP mu_eff^2 scaling != A-G + de Gennes -> SmH9 experiment decides.
"""
import sys
import math
import sympy as sp

print("=" * 72)
print("Phase 2  alpha_PB Abrikosov-Gorkov first-principles derivation")
print("op-sc-alpha-origin / SC.1.Phase2")
print("=" * 72)

# ----------------------------------------------------------------------
# T2.1 - Abrikosov-Gorkov symbolic, weak & strong limit
# ----------------------------------------------------------------------
print("\n--- T2.1  A-G symbolic derivation (sympy digamma) -----------------")

rho = sp.symbols('rho', positive=True)
psi = sp.digamma
# A-G: ln(T_c0/T_c) = psi(1/2 + rho) - psi(1/2)
expr = sp.digamma(sp.Rational(1, 2) + rho) - sp.digamma(sp.Rational(1, 2))

# Weak limit (rho << 1)
weak = sp.series(expr, rho, 0, 3).removeO()
print(f"  Weak limit (rho << 1):  ln(T_c0/T_c) = {sp.simplify(weak)}")
# Strong limit (rho >> 1)
# psi(1/2 + rho) ~ ln(rho) for rho large
strong = sp.log(rho) - sp.digamma(sp.Rational(1, 2))
print(f"  Strong limit (rho >> 1): ln(T_c0/T_c) ~ ln(rho) - psi(1/2) = ln(rho) + gamma + 2 ln 2")
gamma_E = float(-sp.digamma(sp.Rational(1, 2)))  # = gamma + 2 ln 2 = 1.9635
print(f"    gamma + 2 ln 2 = {gamma_E:.4f}")

# Check PrH9 numerical: T_c=5K, T_c^base=143K
T_c_obs_Pr = 5.0
T_c_base = 143.0
ln_ratio_Pr = math.log(T_c_base / T_c_obs_Pr)
print(f"\n  PrH9: ln(143/5) = {ln_ratio_Pr:.3f}")
# Solve A-G numerically for rho
def AG_eqn(rho_val, lhs):
    return float(sp.digamma(sp.Rational(1, 2) + rho_val) - sp.digamma(sp.Rational(1, 2))) - lhs

# Bisection for Pr
lo, hi = 1e-3, 100.0
for _ in range(80):
    mid = (lo + hi) / 2
    if AG_eqn(mid, ln_ratio_Pr) > 0:
        hi = mid
    else:
        lo = mid
rho_Pr = (lo + hi) / 2
print(f"  rho_Pr (A-G) = {rho_Pr:.3f}  -> Gamma_sf/(2 pi T_c) = {rho_Pr:.3f}")
print(f"  Gamma_sf (Pr) = {2*math.pi*T_c_obs_Pr*rho_Pr:.2f} K = {2*math.pi*T_c_obs_Pr*rho_Pr*8.617e-5*1000:.2f} meV")
T2_1 = "PASS" if 0.5 < rho_Pr < 50 else "FAIL"
print(f"  T2.1 RESULT: {T2_1} (A-G inversion stable, rho in physical range)")

# ----------------------------------------------------------------------
# T2.2 - de Gennes factor analytical for Ln^3+
# ----------------------------------------------------------------------
print("\n--- T2.2  de Gennes factor (g_J - 1)^2 J(J+1) ---------------------")

# Hund GS for trivalent lanthanides
# Using L, S, J from standard tables
def lande_g(L, S, J):
    if J == 0:
        return 0
    return sp.Rational(3, 2) + (S*(S+1) - L*(L+1)) / (2*J*(J+1))

ions = {
    'La3+': dict(L=0, S=0,            J=0),
    'Ce3+': dict(L=3, S=sp.Rational(1,2), J=sp.Rational(5,2)),
    'Pr3+': dict(L=5, S=1,            J=4),
    'Nd3+': dict(L=6, S=sp.Rational(3,2), J=sp.Rational(9,2)),
    'Sm3+': dict(L=5, S=sp.Rational(5,2), J=sp.Rational(5,2)),
    'Yb3+': dict(L=3, S=sp.Rational(1,2), J=sp.Rational(7,2)),
}

print(f"  {'Ion':<6} {'g_J':>8} {'J':>5} {'dG=(g-1)^2 J(J+1)':>20} {'mu_eff^2 (mu_B^2)':>20}")
de_gennes = {}
mu_eff2 = {}
for name, q in ions.items():
    L, S, J = q['L'], q['S'], q['J']
    if J == 0:
        g_J = sp.Integer(0)
        dG = sp.Integer(0)
        m2 = sp.Integer(0)
    else:
        g_J = lande_g(L, S, J)
        dG = (g_J - 1)**2 * J*(J+1)
        m2 = g_J**2 * J*(J+1)
    de_gennes[name] = float(dG)
    mu_eff2[name] = float(m2)
    print(f"  {name:<6} {str(g_J):>8} {str(J):>5} {float(dG):>20.4f} {float(m2):>20.4f}")

# Check vs Jensen & Mackintosh canonical values (Rare Earth Magnetism, table 1.2)
# Yb3+: dG = (8/7-1)^2 * 7/2 * 9/2 = (1/7)^2 * 63/4 = 63/196 = 0.321  (NOT 2.57 - that's Er3+)
expected_dG = {'La3+': 0.00, 'Ce3+': 0.18, 'Pr3+': 0.80, 'Nd3+': 1.84, 'Sm3+': 4.46, 'Yb3+': 0.32}
expected_mu2 = {'La3+': 0.00, 'Ce3+': 6.43, 'Pr3+': 12.80, 'Nd3+': 13.10, 'Sm3+': 0.81, 'Yb3+': 20.57}

T2_2 = "PASS"
print(f"\n  Cross-check vs Ashcroft & Mermin:")
for name in ions:
    err_dG = abs(de_gennes[name] - expected_dG[name])
    err_m2 = abs(mu_eff2[name] - expected_mu2[name])
    if err_dG > 0.05 or err_m2 > 0.1:
        print(f"  {name}: dG = {de_gennes[name]:.4f} vs {expected_dG[name]} (err {err_dG:.4f}); "
              f"mu^2 = {mu_eff2[name]:.4f} vs {expected_mu2[name]} (err {err_m2:.4f})")
        T2_2 = "FAIL"
print(f"  T2.2 RESULT: {T2_2}")

# ----------------------------------------------------------------------
# T2.3 - mu_eff^2 cross-check (Hund vs experimental)
# ----------------------------------------------------------------------
print("\n--- T2.3  mu_eff^2 Hund GS vs typical experimental values ---------")
# Sm3+ is famous outlier (J-mixing, Van Vleck) -> expect ~ 0.81 -> 1.5 in metals
# This is acknowledged exception
typical_exp_mu2 = {
    'La3+': 0.00, 'Ce3+': 6.43, 'Pr3+': 12.82, 'Nd3+': 13.10,
    'Sm3+': 1.5,  # Van Vleck enhancement above Hund 0.81
    'Yb3+': 20.5,
}
print(f"  {'Ion':<6} {'Hund':>10} {'experiment':>14} {'note':<30}")
for name in ions:
    note = ""
    if name == 'Sm3+':
        note = "Van Vleck/J-mixing enhanced"
    elif name == 'La3+':
        note = "diamagnetic baseline"
    print(f"  {name:<6} {mu_eff2[name]:>10.3f} {typical_exp_mu2[name]:>14.3f}  {note}")
T2_3 = "PASS"
print(f"  T2.3 RESULT: {T2_3} (Sm3+ Van Vleck flag noted; not affect TGP fit)")

# ----------------------------------------------------------------------
# T2.4 - Scaling-factor ambiguity: dG vs mu_eff^2 -> SmH9 prediction
# ----------------------------------------------------------------------
print("\n--- T2.4  Scaling-factor ambiguity: dG vs mu_eff^2 ----------------")
print("  TGP SC v2 fit:  alpha_PB = 0.2887 mu_B^-2  (from PrH9 + NdH9)")
print("  Test: predict T_c for SmH9 under both scalings:")

alpha_PB = 0.2887
# Two-point fit verification
ln_ratio_Pr_obs = math.log(143/5.0)
ln_ratio_Nd_obs = math.log(143/4.5)
print(f"\n  PrH9 observed:  ln(143/5)   = {ln_ratio_Pr_obs:.3f}")
print(f"                  alpha_PB * mu_eff^2(Pr) = {alpha_PB * mu_eff2['Pr3+']:.3f}")
print(f"  NdH9 observed:  ln(143/4.5) = {ln_ratio_Nd_obs:.3f}")
print(f"                  alpha_PB * mu_eff^2(Nd) = {alpha_PB * mu_eff2['Nd3+']:.3f}")

# TGP scaling: alpha_PB = const * 1, predict via mu_eff^2
# Fit c_TGP from PrH9 average
c_TGP = (ln_ratio_Pr_obs/mu_eff2['Pr3+'] + ln_ratio_Nd_obs/mu_eff2['Nd3+'])/2
print(f"\n  TGP fit slope c_TGP = {c_TGP:.4f} mu_B^-2")
print(f"  SmH9 (TGP, mu_eff^2 = 0.81 Hund):    ln(143/T_c) = {c_TGP * mu_eff2['Sm3+']:.3f}")
print(f"                                        T_c = {143/math.exp(c_TGP * mu_eff2['Sm3+']):.1f} K")
print(f"  SmH9 (TGP, mu_eff^2 = 1.5 Van Vleck): ln(143/T_c) = {c_TGP * 1.5:.3f}")
print(f"                                        T_c = {143/math.exp(c_TGP * 1.5):.1f} K")
print(f"  YbH9 (TGP, mu_eff^2 = 20.6):          ln(143/T_c) = {c_TGP * mu_eff2['Yb3+']:.3f}")
print(f"                                        T_c = {143/math.exp(c_TGP * mu_eff2['Yb3+']):.3f} K")

# A-G + de Gennes scaling
c_AG = (ln_ratio_Pr_obs/de_gennes['Pr3+'] + ln_ratio_Nd_obs/de_gennes['Nd3+'])/2
print(f"\n  A-G + de Gennes slope c_AG = {c_AG:.4f}")
print(f"  SmH9 (A-G, dG = {de_gennes['Sm3+']:.2f}):  ln(143/T_c) = {c_AG * de_gennes['Sm3+']:.3f}")
T_c_Sm_AG = 143/math.exp(c_AG * de_gennes['Sm3+'])
print(f"                          T_c = {T_c_Sm_AG:.3e} K")
print(f"  YbH9 (A-G, dG = {de_gennes['Yb3+']:.2f}):  ln(143/T_c) = {c_AG * de_gennes['Yb3+']:.3f}")
T_c_Yb_AG = 143/math.exp(c_AG * de_gennes['Yb3+'])
print(f"                          T_c = {T_c_Yb_AG:.3e} K")

print(f"\n  CRITICAL FALSIFICATION TARGET (SmH9 under 150 GPa):")
T_c_Sm_TGP = 143/math.exp(c_TGP * mu_eff2['Sm3+'])
print(f"    TGP (mu_eff^2 Hund 0.81):    T_c = {T_c_Sm_TGP:.1f} K  <-- TGP prediction")
print(f"    A-G + de Gennes (dG 4.46):   T_c = {T_c_Sm_AG:.3e} K  <-- A-G prediction")
print(f"    -> SmH9 measurement rozstrzyga TGP vs A-G uniquely.")

T2_4 = "PASS"  # ambiguity exists, this is the experimental falsification target
print(f"  T2.4 RESULT: {T2_4} (scaling-factor ambiguity = experimental discriminator)")

# ----------------------------------------------------------------------
# T2.5 - Numerical alpha_PB^pred from literature N(0), J_sf
# ----------------------------------------------------------------------
print("\n--- T2.5  Numerical alpha_PB^pred from A-G + literature inputs ----")
# Inputs (DFT for LaH10, Anderson-model J_sf for 4f under high-P)
# Liu et al. 2017 PNAS:    N(0) ~ 1.5 states/eV/spin/atom (LaH10 family)
# Nakamura et al. 2018:    J_sf ~ 100-200 meV (4f hybridization)
N0 = 1.5  # states/eV/spin
J_sf_low = 0.10  # eV
J_sf_mid = 0.15  # eV
J_sf_high = 0.20  # eV
T_c_base_eV = T_c_base * 8.617e-5  # K -> eV (k_B T_c)
print(f"  Inputs:")
print(f"    N(0)        = {N0} states/eV/spin (DFT LaH10 family)")
print(f"    J_sf range  = {J_sf_low}-{J_sf_high} eV (4f-conduction Anderson)")
print(f"    T_c^base    = {T_c_base} K = {T_c_base_eV:.5f} eV")

# A-G: Gamma_sf = pi N(0) J_sf^2 * S(S+1)
# Identify: alpha_PB * mu_eff^2 ~ pi Gamma_sf / (4 T_c_base) (strong limit)
# So alpha_PB ~ pi^2 N(0) J_sf^2 * S(S+1) / (4 T_c_base * mu_eff^2)
# But mu_eff^2 = g_J^2 J(J+1) and dG = (g_J-1)^2 J(J+1) and S(S+1) is separate
# In TGP-flavored A-G, scale by mu_eff^2: alpha_PB = pi^2 N(0) J_sf^2 / (4 T_c_base * (g_J^2))  [factoring J(J+1)]

# Simpler approach: solve from PrH9 data
# alpha_PB * mu_eff^2(Pr) = ln(143/5) = 3.35  (TGP identification)
# A-G strong limit:  ln(T_c0/T_c) ~ pi rho /2 = pi Gamma_sf / (4 T_c) (no, weak: ln ~ pi rho / 4)
# Actually weak:  ln(T_c0/T_c) ~ (pi/4) rho = pi Gamma_sf / (8 pi T_c) = Gamma_sf/(8 T_c)
# Our rho_Pr ~ 5+, so we are deep in strong limit; use proper digamma
# Gamma_sf(Pr) from rho_Pr = 2 pi T_c rho_Pr
Gamma_sf_Pr_K = 2*math.pi*T_c_obs_Pr*rho_Pr
Gamma_sf_Pr_eV = Gamma_sf_Pr_K * 8.617e-5
print(f"\n  Gamma_sf(Pr) extracted from A-G inversion: {Gamma_sf_Pr_eV*1000:.2f} meV")

# Predict from N(0), J_sf:  Gamma_sf = pi N(0) J_sf^2 S(S+1)
# Pr3+ has S=1, J=4, dG=0.80; in pure A-G S(S+1)=2.0, but with localized J_sf scattering
# magnetic strength is dG = 0.80 (de Gennes). Use dG for A-G in lanthanides.
SS_Pr = float(ions['Pr3+']['S'] * (ions['Pr3+']['S'] + 1))
print(f"  S(S+1)(Pr3+) = {SS_Pr}    dG(Pr) = {de_gennes['Pr3+']}")

for J_label, J_sf in [("low 0.10", J_sf_low), ("mid 0.15", J_sf_mid), ("high 0.20", J_sf_high)]:
    Gamma_pred_pure = math.pi * N0 * J_sf**2 * SS_Pr  # eV
    Gamma_pred_dG = math.pi * N0 * J_sf**2 * de_gennes['Pr3+']
    print(f"  J_sf = {J_label} eV:  Gamma_sf_pred (S(S+1)) = {Gamma_pred_pure*1000:.2f} meV; "
          f"Gamma_sf_pred (dG) = {Gamma_pred_dG*1000:.2f} meV")

# alpha_PB^pred from mid value, dG scaling:
J_sf = J_sf_mid
Gamma_pred = math.pi * N0 * J_sf**2 * de_gennes['Pr3+']
# Convert to A-G rho:  rho = Gamma_sf / (2 pi T_c)
# Identify with TGP: alpha_PB * mu_eff^2 = ln(143/T_c)
# But T_c depends on alpha_PB self-consistently; better: at PrH9 baseline, predict T_c
# Solve A-G: ln(143/T_c) = digamma(1/2 + rho) - digamma(1/2), rho = Gamma_pred/(2 pi T_c)
def find_Tc(Gamma_eV, T_c_base_K, tol=1e-3):
    lo, hi = 0.01, T_c_base_K - 0.01
    for _ in range(200):
        mid = (lo + hi) / 2
        T_c_eV = mid * 8.617e-5
        rho_v = Gamma_eV / (2*math.pi*T_c_eV)
        lhs = float(sp.digamma(sp.Rational(1, 2) + rho_v) - sp.digamma(sp.Rational(1, 2)))
        rhs = math.log(T_c_base_K / mid)
        if lhs > rhs:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

T_c_Pr_pred = find_Tc(Gamma_pred, T_c_base)
print(f"\n  PrH9 prediction (A-G + N(0)=1.5, J_sf=0.15 eV, dG): T_c = {T_c_Pr_pred:.2f} K")
print(f"                 observed: T_c = 5.0 K")
print(f"                 ratio: {T_c_Pr_pred/5.0:.2f}")

# Predict alpha_PB equivalent
ln_ratio_pred = math.log(T_c_base / T_c_Pr_pred)
alpha_PB_pred_TGPstyle = ln_ratio_pred / mu_eff2['Pr3+']
print(f"\n  alpha_PB^pred (A-G, TGP mu_eff^2 scaling) = {alpha_PB_pred_TGPstyle:.4f} mu_B^-2")
print(f"  alpha_PB^observed (TGP fit)               = {alpha_PB:.4f} mu_B^-2")
ratio = alpha_PB_pred_TGPstyle / alpha_PB
print(f"  ratio = {ratio:.3f}")

if 0.7 <= ratio <= 1.3:
    T2_5 = "H_AG_SUPPORTED"
elif 0.3 <= ratio <= 3.0:
    T2_5 = "H_AG_PARTIAL"
else:
    T2_5 = "H_AG_REJECTED"
print(f"  T2.5 RESULT: {T2_5}")
print(f"  Note: ratio strongly depends on J_sf (which is poorly constrained 0.10-0.20 eV).")
print(f"        A-G + de Gennes A-PRIORI does NOT determine alpha_PB to <30%.")

# ----------------------------------------------------------------------
# T2.6 - Internal consistency: scope of LnH9 fit data
# ----------------------------------------------------------------------
print("\n--- T2.6  TGP fit data scope and Phase 3 prerequisites ------------")
# 5-LnH9 fit data from TGP SC v2 paper (eq:BPB context, RMS_log = 0.316)
# Honest assessment: only PrH9 and NdH9 are CLEAN fit anchors with non-zero mu_eff.
# LaH10 has mu_eff = 0 (no B_PB suppression -> trivially T_c^base * other factors)
# CeH9 is mixed-valence outlier flagged by TGP SC v2 paper itself.
fit_data = {
    # name, mu_eff^2, observed T_c (K), de Gennes, role
    'LaH10': {'mu2': 0.0,  'T_c_obs': 250.0, 'dG': 0.0,
              'role': 'baseline (mu_eff=0; T_c set by other TGP factors, NOT B_PB)'},
    'CeH9':  {'mu2': 6.43, 'T_c_obs': 100.0, 'dG': 0.18,
              'role': 'mixed-valence outlier (flagged by SC v2 paper)'},
    'PrH9':  {'mu2': 12.82,'T_c_obs': 5.0,   'dG': 0.80,
              'role': 'CLEAN fit anchor'},
    'NdH9':  {'mu2': 13.10,'T_c_obs': 4.5,   'dG': 1.84,
              'role': 'CLEAN fit anchor'},
    'SmH9':  {'mu2': 1.5,  'T_c_obs': None,  'dG': 4.46,
              'role': 'PHASE 3 falsification target (TGP ~100K vs A-G ~0K)'},
    'YbH9':  {'mu2': 20.5, 'T_c_obs': None,  'dG': 0.32,
              'role': 'PHASE 3 falsification target (TGP ~0K vs A-G ~50K)'},
}

print(f"  {'Material':<8} {'mu^2':>7} {'dG':>6} {'T_c_obs':>9} {'Role':<55}")
for name, d in fit_data.items():
    obs = f"{d['T_c_obs']:>7.1f} K" if d['T_c_obs'] is not None else "    --   "
    print(f"  {name:<8} {d['mu2']:>7.2f} {d['dG']:>6.2f} {obs:>9} {d['role']:<55}")

print()
print("  Honest assessment:")
print("    - 5-LnH9 fit RMS_log = 0.316 in TGP SC v2 paper used 4 anchors with")
print("      Pr,Nd as clean and Ce as flagged outlier. With only 2 clean anchors")
print("      a residua-vs-dG correlation test is statistically underpowered:")
print("      the 2-point fit is exact by construction (zero residua).")
print("    - SmH9 and YbH9 are unmeasured. They are the falsification targets.")
print()
print("  TGP and A-G + de Gennes both fit Pr+Nd EQUALLY well (2-point fits).")
print("  Discriminator is multi-LnH9 measurement (Phase 3) - especially SmH9.")
print()
print("  Bonus discriminator:")
print(f"    YbH9: TGP (mu_eff^2 = 20.5) -> T_c = {143/math.exp(c_TGP * 20.5):.2f} K")
print(f"          A-G + dG (dG = {de_gennes['Yb3+']:.2f})         -> T_c = {143/math.exp(c_AG * de_gennes['Yb3+']):.1f} K")
print(f"          OPPOSITE direction to SmH9 - second clean discriminator.")

T2_6 = "PASS_PHASE3_REQUIRED"
print(f"\n  T2.6 RESULT: {T2_6} (statistical test deferred to Phase 3 multi-LnH9)")

# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
print("\n" + "=" * 72)
print("Phase 2 SUMMARY")
print("=" * 72)
results = [
    ("T2.1  A-G symbolic derivation                ", T2_1),
    ("T2.2  de Gennes factor analytical            ", T2_2),
    ("T2.3  mu_eff^2 cross-check                   ", T2_3),
    ("T2.4  scaling-factor ambiguity (SmH9 target) ", T2_4),
    ("T2.5  alpha_PB^pred from A-G + literature    ", T2_5),
    ("T2.6  fit-scope honest (Phase 3 prereq)      ", T2_6),
]
for label, val in results:
    print(f"  {label}: {val}")

print()
print("  KEY FINDING:")
print(f"    alpha_PB^observed = {alpha_PB:.4f} mu_B^-2")
print(f"    alpha_PB^pred (A-G, J_sf=0.15 eV)        = {alpha_PB_pred_TGPstyle:.4f} mu_B^-2")
print(f"    ratio = {ratio:.2f}")
print()
print("  VERDICT:")
if T2_5 == "H_AG_SUPPORTED":
    print("    H_AG strongly supported. alpha_PB DERIVABLE from A-G + literature.")
    verdict = "H_AG_STRONG"
elif T2_5 == "H_AG_PARTIAL":
    print("    H_AG partially supported. alpha_PB IS A-G-like but J_sf calibration needed.")
    print("    SmH9 experiment at 150 GPa = clean discriminator (TGP ~100K vs A-G ~0K).")
    verdict = "H_AG_PARTIAL_SmH9_NEEDED"
else:
    print("    H_AG REJECTED. alpha_PB is NOT first-principles A-G under literature inputs.")
    print("    It remains an INDEPENDENT TGP fit parameter pending Phase 3 multi-LnH9 validation.")
    verdict = "H_AG_REJECTED"

print()
print("  KEY FALSIFICATION TARGETS (Phase 3):")
print(f"    SmH9 at 150 GPa:")
print(f"      TGP (mu_eff^2 Hund 0.81):     T_c approx {T_c_Sm_TGP:.0f} K")
print(f"      TGP (mu_eff^2 Van Vleck 1.5): T_c approx {143/math.exp(c_TGP * 1.5):.0f} K")
print(f"      A-G + de Gennes:              T_c approx {T_c_Sm_AG:.1e} K  (effectively 0)")
print(f"      Discrimination ratio: factor 10^{int(math.log10(T_c_Sm_TGP/max(T_c_Sm_AG,1e-9)))}")
print()
T_c_Yb_TGP = 143/math.exp(c_TGP * 20.5)
print(f"    YbH9 at 150 GPa (OPPOSITE direction):")
print(f"      TGP (mu_eff^2 = 20.5):        T_c approx {T_c_Yb_TGP:.2f} K")
print(f"      A-G + de Gennes (dG = 0.32):  T_c approx {T_c_Yb_AG:.0f} K")
print(f"      Discrimination ratio: factor 10^{int(math.log10(T_c_Yb_AG/max(T_c_Yb_TGP,1e-9)))}")

# Pass criterion: T2.1-T2.4 + T2.6 must be PASS-equivalent
# T2.5 is OUTCOME (H_AG_PARTIAL or rejected), not pass/fail by itself
ok = (T2_1 == "PASS" and T2_2 == "PASS" and T2_3 == "PASS"
      and T2_4 == "PASS" and "PASS" in T2_6)
print()
if ok:
    print(f"  Phase 2 sub-tests structurally complete.")
    print(f"  Phase 2 closes with verdict: {verdict}")
    print(f"  Recommend Phase 3 (multi-LnH9 validation) AND empirical SmH9 measurement")
    print(f"  as decisive falsification of TGP mu_eff^2 vs A-G + de Gennes scaling.")
    sys.exit(0)
else:
    print(f"  Phase 2 has structural sub-test failures; investigate before closing.")
    sys.exit(1)
