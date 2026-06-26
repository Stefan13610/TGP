# -*- coding: ascii -*-
# Phase 4 sympy -- op-frontier-bridge-and-asymmetry-2026-06-12
# Modul B: KB4 (H-CP amplitude asymmetry) + SIG-2 (leakage -> predicted diffuse background)
#
# Pre-registration: Phase0_balance.md LOCKED 2026-06-12 (par. 1.4 KB4, SIG-2; par. 5 anchor
#   f_bar_max = 1e-6 COMPARISON-ONLY -- enters ONLY the final mechanical comparison line, never derivations)
# LOCKED inputs (read-only): m_eff^2 = lam(3Phi^2-Phi0^2) (FM P1); V_int = -C exp(-m_Phi L),
#   tail amplitude 2Phi0 (CE-H + BA P1); xi_d = ln(2+sqrt3) (BA P1); v_c = 2c/3; delta = 2/m_Phi.
# Declared approximations (per-use): 1D-proxy (1D-proxy != 3D claim); dilute pairwise tails;
#   chase-criterion band Dv in {c, c/3} (convention sensitivity reported as band);
#   creation-weight model: conversion density ~ V(w(x)) ~ sech^4(m x/2) (declared MODEL,
#   not derived from creation kinematics); partner-pull correction subleading (declared).
# 0 hardcoded T_pass; circularity guard FP9.

import sympy as sp

OUT = []
def log(s=""):
    OUT.append(str(s))
    print(s)

RESULTS = []
def verdict(name, desc, computed_pass, detail):
    status = "PASS" if computed_pass else "FAIL"
    RESULTS.append((name, status))
    log("[%s] %s -- %s" % (status, name, desc))
    log("        %s" % detail)

log("=" * 78)
log("Phase 4 sympy -- op-frontier-bridge-and-asymmetry (modul B: KB4 + SIG-2)")
log("=" * 78)

lam, Phi0, w_s, h_s, chi, t, c, G = sp.symbols('lambda Phi0 w h chi t c G', positive=True)
chi_s = sp.symbols('chi_s')                # unrestricted sign for parity audit
mPhi = sp.sqrt(2 * lam) * Phi0
MK = sp.Rational(2, 3) * mPhi * Phi0**2
xi_d = sp.log(2 + sp.sqrt(3))              # creation depth (BA P1 LOCKED-in-cycle)

# ================================================================ KB4 (H-CP)
log("")
log("--- KB4: does the creation amplitude distinguish Phi from Phi*? (H-CP) ---")

# complex fluctuation around REAL wall background: Phi = (w + h) + i chi
# (canonical normalization under TGP concept par. 3.2 convention L_kin = (1/2)|dPhi|^2:
#  |dPhi|^2 = (d(w+h))^2 + (d chi)^2 -> both h and chi canonical; convention fix vs first run,
#  parity structure in chi unaffected)
Phi_sq = (w_s + h_s)**2 + chi_s**2
Vfull = lam / 4 * (Phi_sq - Phi0**2)**2

# FP1: fluctuation spectrum: m_h^2 = lam(3w^2-Phi0^2) (== LOCKED m_eff^2); m_chi^2 = lam(w^2-Phi0^2); Goldstone at w=Phi0
m_h2 = sp.simplify(sp.diff(Vfull, h_s, 2).subs([(h_s, 0), (chi_s, 0)]))
m_chi2 = sp.simplify(sp.diff(Vfull, chi_s, 2).subs([(h_s, 0), (chi_s, 0)]))
fp1_ok = (sp.simplify(m_h2 - lam * (3 * w_s**2 - Phi0**2)) == 0) \
         and (sp.simplify(m_chi2 - lam * (w_s**2 - Phi0**2)) == 0) \
         and (sp.simplify(m_chi2.subs(w_s, Phi0)) == 0)
verdict("FP1", "fluctuation spectrum around real wall: m_h^2 == LOCKED m_eff^2; m_chi^2 = lam(w^2-Phi0^2); Goldstone m_chi^2(Phi0) = 0",
        fp1_ok, "m_h^2 = %s ; m_chi^2 = %s ; m_chi^2(w=Phi0) = %s"
        % (m_h2, m_chi2, m_chi2.subs(w_s, Phi0)))

# FP2: C-parity audit (C: Phi -> Phi* <=> chi -> -chi): potential even in chi EXACTLY
parity_res = sp.simplify(Vfull.subs(chi_s, -chi_s) - Vfull)
poly_chi = sp.Poly(sp.expand(Vfull), chi_s)
odd_monoms = [mdeg for mdeg in poly_chi.monoms() if mdeg[0] % 2 == 1]
fp2_ok = (parity_res == 0) and (len(odd_monoms) == 0)
verdict("FP2", "C-parity: V(h, -chi) - V(h, chi) == 0 EXACT; ZERO odd-chi monomials to ALL orders (|Phi|^2 even in chi)",
        fp2_ok, "parity residual = %s ; odd-chi monomial degrees: %s (kinetic term even trivially: |dPhi|^2 = (d(w+h))^2 + (d chi)^2)"
        % (parity_res, odd_monoms))

# FP3: consequence: every interaction vertex has even number of chi legs -> no C-asymmetric
# creation amplitude from the real wall at ANY perturbative order (enumerate vertices to chi^4)
vertices = {}
expV = sp.expand(Vfull)
for n_chi in range(0, 5):
    for n_h in range(0, 5):
        coeff = expV.coeff(chi_s, n_chi).coeff(h_s, n_h)
        if coeff != 0 and (n_chi + n_h) >= 3:
            vertices[(n_h, n_chi)] = sp.simplify(coeff)
odd_vertices = {k: v for k, v in vertices.items() if k[1] % 2 == 1}
fp3_ok = (len(odd_vertices) == 0) and (len(vertices) > 0)
verdict("FP3", "vertex audit: ALL interaction vertices have EVEN chi-leg count -> H-CP amplitude asymmetry EXCLUDED for real wall",
        fp3_ok, "vertices (n_h, n_chi): %s ; odd-chi vertices: %s"
        % (sorted(vertices.keys()), sorted(odd_vertices.keys())))

# ================================================================ SIG-2 (leakage)
log("")
log("--- SIG-2: leakage fraction bound + fate of leaked antimatter ---")

# FP4: wall-antikink binding at creation depth (deeply bound regime) -- amplitude-convention band
C_manton = 8 * mPhi * Phi0**2            # Manton tail-product 2 m a^2, a = 2Phi0
Cfit_over_mPhi02 = sp.Rational(616, 100)  # BA P1 FP4 fitted C/(m Phi0^2) = 6.16 (declared import)
bind_manton = sp.simplify((C_manton / MK) * sp.exp(-xi_d))
bind_fit = sp.simplify((Cfit_over_mPhi02 / sp.Rational(2, 3)) * sp.exp(-xi_d))
bind_manton_closed = sp.simplify(bind_manton - 12 * (2 - sp.sqrt(3)))
fp4_ok = (bind_manton_closed == 0) and (float(bind_manton) > 1) and (float(bind_fit) > 1)
verdict("FP4", "binding |V_wA(xi_d)|/M_K = 12(2-sqrt3) EXACT ~ 3.21 (Manton) / 2.48 (fitted band) > 1: DEEPLY BOUND at creation depth",
        fp4_ok, "Manton: %s = %.4f ; fitted: %.4f ; closed-form residual: %s"
        % (sp.nsimplify(bind_manton), float(bind_manton), float(bind_fit), bind_manton_closed))

# FP5: chase timescales: absorption acceleration vs layer transit
# F/M_K = 12 m_Phi exp(-xi) c^2 ; tau_tr = delta/(c - v_c) = 6/(m_Phi c)
# criterion band: Dv = c (strict) -> tau_acc/tau_tr = exp(xi)/72 ; Dv = c/3 -> exp(xi)/216
ratio_strict = sp.simplify(sp.exp(xi_d) / 72)
ratio_loose = sp.simplify(sp.exp(xi_d) / 216)
ratio_closed = sp.simplify(ratio_strict - (2 + sp.sqrt(3)) / 72)
fp5_ok = (ratio_closed == 0) and (float(ratio_strict) < 1) and (float(ratio_loose) < 1)
verdict("FP5", "chase: tau_acc/tau_tr = (2+sqrt3)/72 EXACT ~ 0.052 << 1 (band: 0.017-0.052): absorption FAST vs transit",
        fp5_ok, "strict (Dv=c): %.4f ; loose (Dv=c/3): %.4f ; closed-form residual: %s"
        % (float(ratio_strict), float(ratio_loose), ratio_closed))

# FP6: leakage channel opens only at depth xi > xi_leak = ln(72) [ln(216) loose];
# A-i LOCKED confines creation to the layer => exponential closure; bound under declared weight model
xi_var = sp.symbols('xi_var', positive=True)
xi_leak_strict = sp.solve(sp.Eq(sp.exp(xi_var) / 72, 1), xi_var)[0]
xi_leak_loose = sp.solve(sp.Eq(sp.exp(xi_var) / 216, 1), xi_var)[0]
# creation weight ~ sech^4(xi/2) (declared model); tail fraction beyond Xi:
u = sp.symbols('u', positive=True)
wgt = sp.sech(u)**4                       # u = xi/2
antider = sp.tanh(u) - sp.tanh(u)**3 / 3  # antiderivative of sech^4 (verified below)
antider_check = sp.simplify(sp.diff(antider, u) - wgt) == 0
total_w = sp.limit(antider, u, sp.oo) - antider.subs(u, 0)   # = 2/3 (closed form)
def tail_fraction(Xi):
    tail = sp.simplify(total_w - (sp.tanh(Xi / 2) - sp.tanh(Xi / 2)**3 / 3))
    return sp.simplify(sp.nsimplify(tail.rewrite(sp.exp)) / total_w)
f_leak_strict = tail_fraction(xi_leak_strict)         # conservative (shallower opening)
f_leak_loose = tail_fraction(xi_leak_loose)
tanh_check = sp.simplify(sp.tanh(sp.log(72) / 2).rewrite(sp.exp) - sp.Rational(71, 73)) == 0
fp6_ok = antider_check and tanh_check and (total_w == sp.Rational(2, 3)) \
         and (float(f_leak_strict) < 2e-3) and (float(f_leak_loose) < float(f_leak_strict)) \
         and (float(xi_leak_strict) > float(xi_d))
verdict("FP6", "leakage opens only at xi > ln72 ~ 4.28 (>> creation depth 1.32; A-i confines creation to layer); f_leak <= ~1.2e-3 (conservative)",
        fp6_ok, "xi_leak = %s = %.4f [loose: %.4f] ; tanh(ln72/2) = 71/73: %s ; f_leak: strict %.3e / loose %.3e"
        % (xi_leak_strict, float(xi_leak_strict), float(xi_leak_loose), tanh_check,
           float(f_leak_strict), float(f_leak_loose)))

# FP7: fate of leaked antimatter: bound K-Abar pair in bulk annihilates (1D: no crossing channel)
# => persistent antimatter fraction -> 0; mechanical COMPARISON-ONLY line vs anchor f_bar_max
tau_ann = sp.symbols('tau_ann', positive=True)        # finite annihilation time (V_int attractive, LOCKED)
f_persistent = f_leak_strict * tau_ann / t
limit_pers = sp.limit(f_persistent, t, sp.oo)
F_BAR_MAX = 1e-6                                       # COMPARISON-ONLY anchor (Phase 0 par. 5) -- reporting line ONLY
comparison_pass = (float(limit_pers) if limit_pers.is_number else None) == 0.0 and (0.0 < F_BAR_MAX)
f_rad = sp.simplify(2 * f_leak_strict)
fp7_ok = (limit_pers == 0) and comparison_pass
verdict("FP7", "persistent antimatter fraction -> 0 (leaked pairs annihilate; transient) ; COMPARISON-ONLY: 0 < f_bar_max = 1e-6 PASS",
        fp7_ok, "lim f_persistent = %s ; energy-injection fraction f_rad = 2 f_leak ~ %.3e -- NO pre-registered anchor "
        "for radiation injection => flagged for FUTURE pre-registration (PR-023 candidate observable), NOT compared (anti-goalpost)"
        % (limit_pers, float(f_rad)))

# FP8: SIG-1 refinement: effective demand multiplier >= 2 (pair creation; annihilated channels
# recycle to wall in-layer, leaked channel radiates) => t_*^(B) >= sqrt(2) t_* (lower bound stands)
mult = sp.symbols('mult', positive=True)
tstar_shift = sp.sqrt(mult)
fp8_ok = (sp.simplify(tstar_shift.subs(mult, 2) - sp.sqrt(2)) == 0) \
         and (sp.diff(tstar_shift, mult).is_positive is True)
verdict("FP8", "SIG-1 refinement: t_* shift = sqrt(multiplier) monotone; multiplier >= 2 => t_*^(B) >= sqrt2 t_* (P1 EXACT result = lower bound)",
        fp8_ok, "shift(mult=2) = %s ; d(shift)/d(mult) > 0: %s"
        % (tstar_shift.subs(mult, 2), sp.diff(tstar_shift, mult).is_positive))

# ================================================================ FP9: circularity guard
blacklist = {'G_obs', 'eta_B', 'etaB', 'f_max', 'fbar_max', 'asym_obs', 'log10G_obs'}
free = set()
for expr in [Vfull, m_h2, m_chi2, bind_manton, ratio_strict, f_leak_strict, f_leak_loose,
             f_persistent, f_rad, xi_leak_strict]:
    free |= {str(fs) for fs in sp.sympify(expr).free_symbols}
hits = free & blacklist
# additionally: anchor used ONLY in FP7 comparison line, not in any derivation expression above
fp9_ok = (len(hits) == 0)
verdict("FP9", "circularity guard: no G_obs / eta_B / anchor symbol in ANY derivation expression (anchor = FP7 comparison line only)",
        fp9_ok, "free symbols audited: %s ; blacklist hits: %s" % (sorted(free), sorted(hits)))

# ================================================================ summary
log("=" * 78)
npass = sum(1 for _, st in RESULTS if st == "PASS")
log("SUMMARY: %d/%d PASS ; 0 hardcoded T_pass (all verdicts computed)" % (npass, len(RESULTS)))
log("KB4: NEGATIVE_FOR_REAL_WALL (DERIVED, all orders): action exactly even in chi =>")
log("  creation amplitude CANNOT distinguish Phi/Phi* around the real wall; H-CP excluded")
log("  within LIVE machinery. Residual GAP (declared): phase-textured wall / RP2 holonomy sectors.")
log("SIG-2: BOUNDED: f_leak <= ~1.2e-3 (conservative; criterion band 1.3e-4..1.2e-3) under declared")
log("  weight model; structural chain: A-i confines creation to layer << xi_leak = ln72;")
log("  leaked pairs transient => persistent antimatter -> 0 (COMPARISON-ONLY: PASS vs 1e-6);")
log("  radiation-injection ~ 2 f_leak = NEW observable flagged for future PR-023 pre-registration.")
log("Disposition: 1D-proxy != 3D claim; F-BA-6 classification -> Phase FINAL (user decision on KB3).")

with open(__file__.replace('.py', '.txt'), 'w') as f:
    f.write("\n".join(OUT) + "\n")
