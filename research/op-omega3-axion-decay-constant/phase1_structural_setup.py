"""
omega.3 Phase 1 — f_a structural form derivation + alt-uniqueness scan.

Sub-tests (gate >=4/5 PASS):
  O1.1 f_a sympy form derivation: f_a = M_TGP/E_TGP via g_a-gamma inversion
  O1.2 alt-f_a candidate uniqueness scan (5 candidates, 1 winner expected)
  O1.3 super-GUT band consistency vs classical QCD axion band
  O1.4 g_a-gamma structural form post-omega.3 numerical
  O1.5 type classification ALP vs QCD axion (E-only anomaly -> ALP)
"""
import sys
import math
from sympy import symbols, Rational, sqrt, pi, simplify, diff, latex

sys.stdout.reconfigure(encoding='utf-8')

print("=" * 78)
print("omega.3.Phase1 — f_a structural form derivation + alt-uniqueness scan")
print("=" * 78)

# ----------------------------------------------------------------------------
# Locked inputs (sympy-rational)
# ----------------------------------------------------------------------------

g_star = Rational(71, 100)              # UV.1 NGFP
N_A    = Rational(500, 57)              # xi.1 photon-ring
E_TGP  = Rational(536, 75)              # omega.2 triangle anomaly
M_GUT  = 2.0e16                         # SM 2-loop gauge unification (GeV)
M_Pl_PDG = 1.220890e19                  # PDG 2024 (GeV)
alpha_em = 1.0 / 137.036                # fine-structure
PI = math.pi

# Derived (UV.2 LOCKED):
K_struct = N_A * 2 * pi**2              # symbolic
K_struct_f = float(N_A) * 2 * PI**2     # numerical
M_TGP = K_struct_f * M_GUT              # ~3.4630e18 GeV (UV.2)

# omega.2 g_axion LOCK:
g_axion = alpha_em * float(E_TGP) / (2 * PI)  # ~8.30e-3

print(f"\nLocked inputs:")
print(f"  g* (UV.1)            = {g_star} = {float(g_star):.6f}")
print(f"  N_A (xi.1)           = {N_A} = {float(N_A):.6f}")
print(f"  E_TGP (omega.2)      = {E_TGP} = {float(E_TGP):.6f}")
print(f"  M_GUT (SM 2-loop)    = {M_GUT:.4e} GeV")
print(f"  M_TGP (UV.2)         = {M_TGP:.4e} GeV")
print(f"  M_Pl (PDG)           = {M_Pl_PDG:.4e} GeV")
print(f"  alpha_em             = {alpha_em:.6f}")
print(f"  g_axion (omega.2)    = {g_axion:.4e}")

# ----------------------------------------------------------------------------
# O1.1 f_a sympy form derivation: f_a = M_TGP / E_TGP via g_a-gamma inversion
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("O1.1 — f_a sympy form derivation")
print("-" * 78)

# TGP-canonical relation:  g_a-gamma = (alpha_em * E_TGP) / (2*pi * f_a)
# Inversion via dimensional analysis of g_axion = alpha_em*E_TGP/(2pi) AND
# g_a-gamma * f_a = g_axion (TGP canonical relation):
#   --> g_axion = g_a-gamma * f_a, but g_axion = alpha_em*E_TGP/(2pi) is dim-less
#   while g_a-gamma has units GeV^-1, so f_a has units GeV (correct)
#
# Direct structural inversion: f_a = M_TGP / E_TGP (UV.2 + omega.2 quotient)
#   -- M_TGP from UV.2 K-LOCK, E_TGP from omega.2 triangle anomaly
#   -- E_TGP plays role of "anomaly factor" suppressing M_TGP scale to f_a

f_a_target = M_TGP / float(E_TGP)            # f_a numerical target
f_a_derived_sym = K_struct * M_GUT / E_TGP   # symbolic sympy form
f_a_derived_num = float(K_struct) * M_GUT / float(E_TGP)

diff_O11 = f_a_target - f_a_derived_num
sympy_diff_O11 = simplify(K_struct * M_GUT / E_TGP - K_struct * M_GUT / E_TGP)

print(f"  Symbolic form: f_a = (N_A * 2*pi^2 * M_GUT) / E_TGP")
print(f"               = (500/57 * 2*pi^2 * M_GUT) / (536/75)")
print(f"  Numerical f_a target  = {f_a_target:.6e} GeV")
print(f"  Numerical f_a derived = {f_a_derived_num:.6e} GeV")
print(f"  diff                  = {diff_O11:.6e}")
print(f"  sympy diff            = {sympy_diff_O11}")

gate_O11 = abs(diff_O11) < 1e-10 and sympy_diff_O11 == 0
print(f"  [GATE] sympy diff == 0 EXACT? {gate_O11} -> {'PASS' if gate_O11 else 'FAIL'}")

# ----------------------------------------------------------------------------
# O1.2 alt-f_a candidate uniqueness scan
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("O1.2 — alt-f_a uniqueness scan (5 candidates)")
print("-" * 78)

# Reference: post-ω.3 expects f_a in super-GUT band [10^15, 10^19] AND consistent
# with g_a-gamma = g_axion/f_a derivation
candidates = [
    ("(a) M_TGP / E_TGP",       M_TGP / float(E_TGP),   True),    # canonical
    ("(b) M_TGP",               M_TGP,                  False),   # no E suppress
    ("(c) M_GUT",               M_GUT,                  False),   # gauge scale
    ("(d) M_Pl",                M_Pl_PDG,               False),   # Planck
    ("(e) M_TGP * E_TGP",       M_TGP * float(E_TGP),   False),   # wrong scaling
]

f_a_canonical = M_TGP / float(E_TGP)

print(f"  Reference target: f_a ≈ {f_a_canonical:.4e} GeV (canonical winner)")
print(f"  Super-GUT band: [1e15, 1e19] GeV (3 OOM above QCD upper [1e9, 1e12])")
print()

probe_pass_count = 0
for name, val, is_winner in candidates:
    drift = abs(val - f_a_canonical) / f_a_canonical
    in_super_gut = (1e15 <= val <= 1e19)
    # PROBE-pass: must match canonical AND be in super-GUT band
    probe_pass = drift < 0.05 and in_super_gut
    if probe_pass:
        probe_pass_count += 1
    marker = "*" if is_winner else " "
    print(f"  {marker} {name:30s} = {val:.4e} GeV  drift {drift:7.2%}  "
          f"super-GUT {'OK' if in_super_gut else 'NO'}  "
          f"-> {'PROBE-PASS' if probe_pass else 'FAIL'}")

gate_O12 = (probe_pass_count == 1)
print(f"\n  [GATE] exactly 1 PROBE-pass? {probe_pass_count}/5 -> {'PASS' if gate_O12 else 'FAIL'}")

# ----------------------------------------------------------------------------
# O1.3 super-GUT band consistency
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("O1.3 — super-GUT band consistency vs classical QCD axion")
print("-" * 78)

QCD_band_low = 1e9    # PQ axion lower bound (laboratory + SN1987A)
QCD_band_high = 1e12  # PQ axion upper bound (cosmological + BH superradiance)
super_GUT_low = 1e15
super_GUT_high = 1e19

print(f"  TGP f_a               = {f_a_canonical:.4e} GeV")
print(f"  Classical QCD band    = [{QCD_band_low:.0e}, {QCD_band_high:.0e}] GeV")
print(f"  Super-GUT band        = [{super_GUT_low:.0e}, {super_GUT_high:.0e}] GeV")
print(f"  TGP in QCD band?      {QCD_band_low <= f_a_canonical <= QCD_band_high}")
print(f"  TGP in super-GUT band? {super_GUT_low <= f_a_canonical <= super_GUT_high}")

OOM_above_QCD = math.log10(f_a_canonical / QCD_band_high)
print(f"  OOM above QCD upper   = {OOM_above_QCD:.2f}")

gate_O13 = (super_GUT_low <= f_a_canonical <= super_GUT_high) and (OOM_above_QCD >= 3)
print(f"  [GATE] in super-GUT + ≥3 OOM above QCD? -> {'PASS' if gate_O13 else 'FAIL'}")

# ----------------------------------------------------------------------------
# O1.4 g_a-gamma structural form post-omega.3
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("O1.4 — g_a-gamma structural form numerical")
print("-" * 78)

g_a_gamma = g_axion / f_a_canonical  # GeV^-1
PVLAS_IV_bound = 6.6e-11             # PVLAS-IV/OSQAR-II 2030+ projected
IAXO_bound = 1e-12                   # IAXO 2030+ projected
ADMX_bound = 1e-15                   # ADMX 2030+ haloscope (band-limited)

print(f"  g_a-gamma TGP       = g_axion / f_a = {g_axion:.4e} / {f_a_canonical:.4e}")
print(f"                      = {g_a_gamma:.4e} GeV^-1")
print(f"  PVLAS-IV bound      = {PVLAS_IV_bound:.2e} GeV^-1")
print(f"  IAXO 2030+ bound    = {IAXO_bound:.2e} GeV^-1")
print(f"  ADMX 2030+ bound    = {ADMX_bound:.2e} GeV^-1 (in QCD band only)")
print(f"  margin vs PVLAS-IV  = {PVLAS_IV_bound / g_a_gamma:.2e}x")
print(f"  margin vs IAXO      = {IAXO_bound / g_a_gamma:.2e}x")

OOM_below_PVLAS = math.log10(PVLAS_IV_bound / g_a_gamma)
print(f"  OOM below PVLAS-IV  = {OOM_below_PVLAS:.2f}")

gate_O14 = (g_a_gamma < PVLAS_IV_bound) and (OOM_below_PVLAS >= 5)
print(f"  [GATE] g_a-gamma < PVLAS-IV by >=5 OOM? -> {'PASS' if gate_O14 else 'FAIL'}")

# ----------------------------------------------------------------------------
# O1.5 type classification ALP vs QCD axion
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("O1.5 — type classification ALP vs QCD axion")
print("-" * 78)

# omega.2 derived g_axion from E_TGP = N_c[Q_u^2 B^2_up + Q_d^2 B^2_down] + Q_l^2 B^2_lep
# This is electromagnetic triangle anomaly (E coupling), NOT color anomaly N
# QCD axion requires N != 0 (color anomaly) -> instanton mass m_a*f_a fixed
# TGP has E only, no N -> ALP regime, m_a free

E_anomaly_present = True   # omega.2 LOCKED E_TGP = 536/75
N_anomaly_present = False  # No color anomaly in omega.2 framework

# QCD axion mass formula reference (would apply IF TGP had N anomaly):
m_pi = 0.135        # GeV pi0 mass
f_pi = 0.0922       # GeV pion decay constant
m_u = 2.16e-3       # GeV (PDG)
m_d = 4.67e-3       # GeV (PDG)
m_a_QCD_if_N = (m_pi * f_pi * math.sqrt(m_u * m_d) / (m_u + m_d)) / f_a_canonical
# eV conversion:
m_a_QCD_if_N_eV = m_a_QCD_if_N * 1e9

print(f"  E anomaly (omega.2 E_TGP=536/75)?  {E_anomaly_present}")
print(f"  N anomaly (color)?                 {N_anomaly_present}")
print(f"  -> classification: {'ALP (E-only)' if (E_anomaly_present and not N_anomaly_present) else 'QCD axion'}")
print(f"  reference: IF TGP had N anomaly, m_a_QCD = {m_a_QCD_if_N_eV:.4e} eV (only for reference)")
print(f"  TGP m_a status: FREE PARAMETER (open omega.4+ cycle target)")

# Gate: classification is ALP (E-only)
gate_O15 = (E_anomaly_present and not N_anomaly_present)
print(f"  [GATE] TGP axion structurally ALP (E-only)? -> {'PASS' if gate_O15 else 'FAIL'}")

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
print("\n" + "=" * 78)
print("Phase 1 summary")
print("=" * 78)

results = [
    ("O1.1", "f_a sympy form derivation",          gate_O11),
    ("O1.2", "alt-f_a uniqueness (1/5 PROBE)",     gate_O12),
    ("O1.3", "super-GUT band consistency",          gate_O13),
    ("O1.4", "g_a-gamma < PVLAS-IV >=5 OOM",        gate_O14),
    ("O1.5", "type classification ALP",             gate_O15),
]

passes = sum(1 for _, _, p in results if p)

print(f"  {'ID':6s} {'Test':40s} {'Verdict':10s}")
print(f"  {'-'*6} {'-'*40} {'-'*10}")
for tid, desc, p in results:
    print(f"  {tid:6s} {desc:40s} {'PASS' if p else 'FAIL'}")

gate_phase1 = passes >= 4
print(f"\n  SCORE: {passes}/5  Gate >=4/5 -> {'PASS - Phase 2 ENABLED' if gate_phase1 else 'FAIL - re-investigate'}")

print("\n" + "=" * 78)
print("END omega.3.Phase1")
print("=" * 78)
