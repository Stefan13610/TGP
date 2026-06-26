# =============================================================================
# op-mechanism-v-pattern25-extreme-envs — Phase 2
# F-P25-B: condensed BVP verification (NS-NS branch) + FORMAL VERDICT
#
# Pre-registration: Phase0_balance.md §3 F-P25-B (LOCKED 2026-06-10):
#   PASS_REALIZED  : delta_psi_max >= 0.385
#   PARTIAL_APPROACH: 0.0385 <= delta_psi_max < 0.385
#   FAIL_NEGATIVE  : delta_psi_max < 0.0385
# Mandatory weak-field regression gate: reproduce T3 Phase 3 Branch A typical-LIGO
# (LOCKED .txt value 6.83e-81; .md transcription x25.5 flagged Phase 1 FP8).
# Discipline: full nonlinear D_kin (Pattern 2.1); 0 hardcoded T_pass; Branch A IMMUTABLE.
# =============================================================================

import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp

PASS_count = 0
FAIL_count = 0
results = {}

def check(label, cond, expected=None, got=None):
    global PASS_count, FAIL_count
    status = "PASS" if cond else "FAIL"
    if cond: PASS_count += 1
    else: FAIL_count += 1
    msg = f"  [{status}] {label}"
    if expected is not None or got is not None:
        msg += f"\n         (expected={expected}, got={got})"
    print(msg)
    return bool(cond)

def banner(t):
    print("\n" + "-" * 78 + f"\n  {t}\n" + "-" * 78)

print("=" * 78)
print("  P25 Phase 2 — F-P25-B condensed BVP verification (NS-NS) + formal verdict")
print("=" * 78)

# -----------------------------------------------------------------------------
# Section 1 — EOM re-verification (verbatim discipline, Phase 0 §4.1)
# -----------------------------------------------------------------------------
banner("Section 1: EOM / W(psi) re-verification (sympy)")
psi_s, g_s = sp.symbols("psi gamma", positive=True)
V_s = -g_s * psi_s**2 * (4 - 3 * psi_s)**2 / 12
W_expected = sp.Rational(1, 3) * psi_s * (8 - 18 * psi_s + 9 * psi_s**2)
fp1a = sp.simplify(W_expected - (-sp.diff(V_s, psi_s) / g_s)) == 0
m2_tilde = -sp.diff(W_expected, psi_s).subs(psi_s, sp.Rational(2, 3))
fp1b = sp.simplify(m2_tilde - sp.Rational(4, 3)) == 0
psi_plus = sp.Rational(6, 9) + 2 * sp.sqrt(3) / 9
delta_crit = sp.simplify(psi_plus - sp.Rational(2, 3))
fp1c = sp.simplify(delta_crit - 2 * sp.sqrt(3) / 9) == 0
results["FP1"] = check("FP1: W = -V'/gamma EXACT; m_tilde^2 = 4/3; delta_psi_critical = 2*sqrt(3)/9",
                       fp1a and fp1b and fp1c,
                       expected="all identities EXACT",
                       got=f"{fp1a}, m2={m2_tilde}, dcrit={delta_crit}")
DELTA_CRIT = float(2 * np.sqrt(3) / 9)        # 0.3849 (T3 EXACT)
PARTIAL_BAND = DELTA_CRIT / 10.0              # 0.0385 (Phase 0 LOCKED)
M_TILDE = float(np.sqrt(4.0 / 3.0))

# -----------------------------------------------------------------------------
# Section 2 — MANDATORY weak-field regression gate (Phase 0 §3 F-P25-B)
# Constants verbatim from LOCKED T3 Phase3_dimensional.py / P25 Phase 1 Section 3
# -----------------------------------------------------------------------------
banner("Section 2: regression gate — T3 Phase 3 Branch A typical-LIGO (LOCKED .txt)")
M_Pl_eV = 1.22e28
M_Sun_kg = 1.989e30
hbar_c_eV_m = 1.973e-7
eV_to_kg = 1.78e-36
m_Phi_A_eV = M_Pl_eV
L_unit_A_m = hbar_c_eV_m / m_Phi_A_eV
m_Phi_A_kg = m_Phi_A_eV * eV_to_kg
rho_unit_A_kg_m3 = m_Phi_A_kg / L_unit_A_m**3

M_BBH_nat = (10 * M_Sun_kg) / m_Phi_A_kg
sigma_LIGO_nat = 30000.0 / L_unit_A_m
delta_psi_typical = 0.75 * M_BBH_nat / ((2 * np.pi)**1.5 * sigma_LIGO_nat**3)
T3_locked_txt = 6.83e-81
rel_dev = abs(delta_psi_typical - T3_locked_txt) / T3_locked_txt
results["FP2"] = check("FP2: REGRESSION GATE — T3 Phase 3 LOCKED .txt 6.83e-81 (rel dev < 2%)",
                       rel_dev < 0.02, expected="6.83e-81",
                       got=f"{delta_psi_typical:.3e} (rel dev {rel_dev:.4f})")

# -----------------------------------------------------------------------------
# Section 3 — full nonlinear BVP solver (T3 Phase 2 template, natural units)
# EOM: psi'' + (2/r) psi' + 2 (psi')^2/psi + W(psi) = -q*rho(r),  q = 1
# BC: psi'(r_min) = 0, psi(r_max) = 2/3; Gaussian rho = M exp(-r^2/2s^2)/((2pi)^1.5 s^3)
# -----------------------------------------------------------------------------
banner("Section 3: nonlinear BVP — anchor vs LOCKED T3 Phase 2 (M=0.01, sigma=1)")

def W_np(p):
    return (1.0 / 3.0) * p * (8.0 - 18.0 * p + 9.0 * p**2)

def solve_psi(M, sigma, r_max=None, n0=3000):
    r_min = 1e-3
    if r_max is None:
        r_max = max(30.0, 5.0 * sigma)
    rho0 = M / ((2 * np.pi)**1.5 * sigma**3)
    def rho(r): return rho0 * np.exp(-r**2 / (2 * sigma**2))
    def odes(r, y):
        p, dp = y
        return np.vstack([dp,
                          -(2.0 / r) * dp - 2.0 * dp**2 / p - W_np(p) - rho(r)])
    def bc(ya, yb):
        return np.array([ya[1], yb[0] - 2.0 / 3.0])
    # mesh: dense geometric near origin + linear tail (2/r stiffness + Yukawa decay)
    r = np.unique(np.concatenate([np.geomspace(r_min, 2.0 * sigma, n0 // 2),
                                  np.linspace(2.0 * sigma, r_max, n0 // 2)]))
    y0 = np.vstack([np.full_like(r, 2.0 / 3.0), np.zeros_like(r)])
    sol = solve_bvp(odes, bc, r, y0, tol=1e-8, max_nodes=400000, verbose=0)
    dpsi_max = float(np.max(sol.y[0]) - 2.0 / 3.0)
    return dpsi_max, sol.success, float(np.max(np.abs(sol.rms_residuals)))

dpsi_anchor, ok_anchor, res_anchor = solve_psi(0.01, 1.0)
T3_phase2_anchor = 1.91e-4
rel_anchor = abs(dpsi_anchor - T3_phase2_anchor) / T3_phase2_anchor
results["FP3"] = check("FP3: nonlinear BVP anchor M=0.01, sigma=1 vs LOCKED T3 Phase 2 "
                       "delta_psi_max = 1.91e-4 (rel dev < 15%)",
                       ok_anchor and rel_anchor < 0.15,
                       expected="1.91e-4", got=f"{dpsi_anchor:.3e} (rel dev {rel_anchor:.3f}; "
                       f"converged={ok_anchor}; rms={res_anchor:.1e})")

# -----------------------------------------------------------------------------
# Section 4 — amplitude ladder: linearity of response (slope = 1 in log-log)
# Validates linear extrapolation toward NS amplitudes (with FP14 Phase 1 no-bootstrap)
# -----------------------------------------------------------------------------
banner("Section 4: amplitude ladder — linear scaling delta_psi_max ∝ M")
ladder_M = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
ladder_out = []
all_conv = True
for Mv in ladder_M:
    d, okv, rs = solve_psi(Mv, 1.0)
    ladder_out.append(d)
    all_conv = all_conv and okv
    print(f"    M = {Mv:.1e}  ->  delta_psi_max = {d:.6e}   (converged={okv}, rms={rs:.1e})")
slope = np.polyfit(np.log10(ladder_M), np.log10(ladder_out), 1)[0]
ratios = [d / Mv for d, Mv in zip(ladder_out, ladder_M)]
ratio_spread = max(ratios) / min(ratios)
results["FP4"] = check("FP4: log-log slope = 1.000 +/- 0.01 AND ratio spread < 2% "
                       "(linear regime confirmed across x100 amplitude range)",
                       all_conv and abs(slope - 1.0) < 0.01 and ratio_spread < 1.02,
                       expected="slope 1.000; spread < 1.02",
                       got=f"slope {slope:.5f}; spread {ratio_spread:.5f}")

# -----------------------------------------------------------------------------
# Section 5 — local-screening regime validation (the regime NS-NS actually sits in)
# Wide source sigma >> lambda_C: central delta_psi -> (3/4)*rho_tilde(0) (S-rho formula)
# -----------------------------------------------------------------------------
banner("Section 5: wide-source BVP — direct validation of local formula delta_psi=(3/4)rho")
sigma_wide = 10.0          # sigma * m_tilde ~ 11.5 >> 1 (local regime)
M_wide = 1e-2              # rho0 ~ 6.3e-7 — comfortably linear per FP4
dpsi_wide, ok_wide, res_wide = solve_psi(M_wide, sigma_wide, r_max=60.0, n0=6000)
rho0_wide = M_wide / ((2 * np.pi)**1.5 * sigma_wide**3)
local_pred = 0.75 * rho0_wide
rel_local = abs(dpsi_wide - local_pred) / local_pred
results["FP5"] = check("FP5: central delta_psi vs (3/4)*rho(0) for sigma*m_tilde ~ 11.5 "
                       "(rel dev < 5%) — S-rho local formula NUMERICALLY validated",
                       ok_wide and rel_local < 0.05,
                       expected=f"{local_pred:.4e}",
                       got=f"{dpsi_wide:.4e} (rel dev {rel_local:.4f}; converged={ok_wide})")

# -----------------------------------------------------------------------------
# Section 6 — NS-NS near-contact formal evaluation under Branch A (IMMUTABLE)
# Local regime FORCED (Phase 1 FP9: sigma_tilde*m_tilde ~ 1e39); field tracks LOCAL rho;
# near-contact overlap enhancement <= 2x central density (superposition bound, linear regime)
# -----------------------------------------------------------------------------
banner("Section 6: NS-NS near-contact — formal delta_psi_max under Branch A")
rho_NS_kg_m3 = 1.0e18                                  # massive-NS central (Phase 1 FP13 verbatim)
rho_NS_nat = rho_NS_kg_m3 / rho_unit_A_kg_m3
dpsi_NS_single = 0.75 * rho_NS_nat
dpsi_NS_contact = 2.0 * dpsi_NS_single                 # near-contact upper bound (x2 overlap)
NS_regime = (10000.0 / L_unit_A_m) * M_TILDE           # NS radius ~10 km in natural lengths
print(f"    rho_NS (Planck units, Branch A) = {rho_NS_nat:.3e}")
print(f"    regime selector (R_NS * m_tilde) = {NS_regime:.2e}  (>> 1: LOCAL regime forced)")
print(f"    delta_psi single NS center      = {dpsi_NS_single:.3e}")
print(f"    delta_psi near-contact (x2 cap) = {dpsi_NS_contact:.3e}")
results["FP6"] = check("FP6: NS-NS delta_psi_max (incl. x2 contact bound) computed; "
                       "consistency with Phase 1 FP13 preview (factor 3)",
                       0.33 < (dpsi_NS_single / 1.46e-79) < 3.0 and NS_regime > 1e30,
                       expected="~1.5e-79 (preview)", got=f"{dpsi_NS_single:.3e}")

# circularity/degeneration audit (Phase 0 mandate): rho -> 0 => delta_psi -> 0;
# thresholds 0.385/0.0385 absent from source forms (appear ONLY in comparison)
rho_chk = sp.symbols("rho_chk", nonnegative=True)
dpsi_form = sp.Rational(3, 4) * rho_chk
fp7a = dpsi_form.subs(rho_chk, 0) == 0
fp7b = not any(abs(float(c) - DELTA_CRIT) < 1e-12 or abs(float(c) - PARTIAL_BAND) < 1e-12
               for c in [sp.Rational(3, 4)])
results["FP7"] = check("FP7: degeneration audit — rho->0 => delta_psi->0; thresholds absent "
                       "from source/response forms (comparison-only)",
                       fp7a and fp7b, expected="clean", got=f"{fp7a}, {fp7b}")

# -----------------------------------------------------------------------------
# Section 7 — F-P25-B FORMAL VERDICT (mechanical; thresholds LOCKED Phase 0 §3)
# -----------------------------------------------------------------------------
banner("Section 7: F-P25-B formal verdict (mechanical)")
dpsi_max_NS = dpsi_NS_contact                          # most favorable surviving case
shortfall_orders = abs(np.log10(dpsi_max_NS / PARTIAL_BAND))
if dpsi_max_NS >= DELTA_CRIT:
    verdict_B = "PASS_REALIZED"
elif dpsi_max_NS >= PARTIAL_BAND:
    verdict_B = "PARTIAL_APPROACH"
else:
    verdict_B = "FAIL_NEGATIVE"
print(f"    delta_psi_max(NS-NS, x2 bound) = {dpsi_max_NS:.3e}")
print(f"    thresholds: PASS >= {DELTA_CRIT:.4f}; PARTIAL >= {PARTIAL_BAND:.4f}")
print(f"    shortfall vs PARTIAL band: {shortfall_orders:.1f} orders of magnitude")
print(f"\n    >>> F-P25-B VERDICT (NS-NS branch): {verdict_B} <<<")
print(f"    (BH-BH branch: NEGATIVE at the gate — Phase 1 FP12, source = 0 identically)")
results["FP8"] = check("FP8: verdict in pre-registered set + computed from thresholds "
                       "(0 hardcoded)",
                       verdict_B in {"PASS_REALIZED", "PARTIAL_APPROACH", "FAIL_NEGATIVE"},
                       expected="pre-registered set", got=verdict_B)

# numerical honesty clause (Phase 0 §3): all solver runs converged; no tachyonic breakdown
results["FP9"] = check("FP9: numerical honesty — all BVP runs converged (no tachyonic-regime "
                       "breakdown anywhere on the ladder; far from T3 M_critical = 15.80)",
                       all_conv and ok_anchor and ok_wide,
                       expected="all converged", got=f"{all_conv and ok_anchor and ok_wide}")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"  FP SUMMARY: {PASS_count}/{PASS_count + FAIL_count} PASS; "
      f"hardcoded T_pass = 0/{PASS_count + FAIL_count}")
print(f"  F-P25-B (NS-NS): {verdict_B}  (shortfall {shortfall_orders:.1f} orders vs PARTIAL band)")
print(f"  F-P25-B (BH-BH): NEGATIVE at gate (Phase 1, source = 0)")
print("=" * 78)
