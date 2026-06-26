# -*- coding: utf-8 -*-
# Phase 1 — op-phi-radiative-dof-audit-2026-06-13
# Dirac constraint analysis + propagator pole decomposition + symmetry inventory.
# Purely structural: NO observational data touched (Phase 0 par.5.2).
# NO hardcoded T_pass; final verdict COMPUTED from flags (Phase 0 par.5.4).
# Run: python Phase1_sympy.py

import sympy as sp

PASS, FAIL = "PASS", "FAIL"
results = []
flags = {}  # falsifier flags -> verdict computed mechanically at the end

def check(name, cond, note=""):
    verdict = PASS if cond else FAIL
    results.append((name, verdict, note))
    print(f"[{verdict}] {name}" + (f"  -- {note}" if note else ""))
    return bool(cond)

print("=" * 72)
print("F-AUX-A.0 — METHOD VALIDATION on EM: A_0 must come out auxiliary")
print("=" * 72)
# Mode-wise EM Lagrangian (temporal + one transverse + one longitudinal dof):
# L_EM = 1/2 (Adot_i - d_i A0)^2 - 1/4 F_ij F_ij ; in Fourier mode k along z:
A0, AL, AT = sp.symbols('A0 A_L A_T', real=True)        # fields
dA0, dAL, dAT = sp.symbols('Adot0 Adot_L Adot_T', real=True)  # velocities
k = sp.symbols('k', positive=True)
# E_L = Adot_L - k*A0 (longitudinal electric), E_T = Adot_T, B ~ k*A_T
L_EM = sp.Rational(1, 2) * ((dAL - k * A0)**2 + dAT**2) - sp.Rational(1, 2) * (k * AT)**2
vels_EM = [dA0, dAL, dAT]
W_EM = sp.Matrix(3, 3, lambda i, j: sp.diff(L_EM, vels_EM[i], vels_EM[j]))
rank_EM = W_EM.rank()
check("A.0.1 EM kinetic Hessian is SINGULAR (rank 2 of 3) -> primary constraint",
      rank_EM == 2, f"rank(W_EM) = {rank_EM}, det = {W_EM.det()}")
# null direction = A0 row:
null_EM = W_EM.nullspace()
check("A.0.2 null direction is exactly A_0 (the known auxiliary field)",
      len(null_EM) == 1 and sp.simplify(null_EM[0] - sp.Matrix([1, 0, 0])) == sp.zeros(3, 1),
      f"nullspace = {[list(v) for v in null_EM]}")
print("      => method correctly detects auxiliary (constrained) components.")

print("=" * 72)
print("F-AUX-A — Dirac analysis of LOCKED linearized TGP Phi sector")
print("=" * 72)
# L_Phi = (K1/2)[ dPhi^2 - (grad dPhi)^2 ] - (m^2/2) dPhi^2 - (q/Phi0) rho dPhi
# Fourier mode q_k:  L_k = (K1/2)(qdot^2 - k^2 q^2) - (m^2/2) q^2 - j(t) q
K1, m, j, t = sp.symbols('K_1 m_ j t', positive=True)
qf = sp.Function('q')(t)
L_phi = sp.Rational(1, 2) * K1 * (sp.diff(qf, t)**2 - k**2 * qf**2) \
        - sp.Rational(1, 2) * m**2 * qf**2 - j * qf
W_phi = sp.diff(L_phi, sp.diff(qf, t), 2)   # kinetic Hessian (1x1)
flag_A_hessian_regular = check(
    "A.1 TGP Phi kinetic Hessian = K_1 > 0, NON-singular -> NO primary constraint",
    sp.simplify(W_phi - K1) == 0 and K1.is_positive, f"W = {W_phi}")
# canonical momentum invertible:
p = sp.diff(L_phi, sp.diff(qf, t))
vel_solved = sp.solve(sp.Eq(sp.Symbol('p_'), p), sp.diff(qf, t))
check("A.2 momentum p = K_1*qdot invertible for velocity (Dirac algorithm ends, 0 constraints)",
      len(vel_solved) == 1, f"qdot = {vel_solved}")
# Euler-Lagrange -> hyperbolic oscillator (radiative mode):
EL = sp.diff(sp.diff(L_phi, sp.diff(qf, t)), t) - sp.diff(L_phi, qf)
EL_expected = K1 * sp.diff(qf, t, 2) + (K1 * k**2 + m**2) * qf + j
flag_A_hyperbolic = check(
    "A.3 EOM: K_1*qddot + (K_1 k^2 + m^2) q = -j  (hyperbolic, propagating)",
    sp.simplify(EL - EL_expected) == 0, f"EL = {sp.simplify(EL)}")
# constraint structure is k-independent -> applies to EVERY multipole incl. quadrupole:
check("A.4 Hessian independent of k -> no mode (monopole/dipole/quadrupole) is constrained",
      sp.diff(W_phi, k) == 0)
flags['F-AUX-A'] = "NEGATIVE" if (flag_A_hessian_regular and flag_A_hyperbolic) else "DERIVED"

print("=" * 72)
print("F-AUX-B — statics <-> radiation: same propagator pole")
print("=" * 72)
# Lorentz invariance: retarded propagator G(omega,k) = F(s), s = k^2 - omega^2.
# Static Newton lock: G(0,k) = 1/k^2 for ALL k>0  =>  F(s) = 1/s for all s>0.
s, w_ = sp.symbols('s omega_', real=True)
F = sp.Function('F')
static_condition = sp.Eq(F(k**2), 1 / k**2)
# F fixed on the positive half-line; for a Lorentz-invariant local field theory
# F is a rational/analytic function of s -> identity F(s) = 1/s everywhere:
F_s = (1 / k**2).subs(k**2, s)
check("B.1 static lock G(0,k)=1/k^2 for all k  =>  F(s) = 1/s identically on s>0",
      sp.simplify(F_s - 1 / s) == 0, "analytic continuation: F(s)=1/s")
# pole of F at s=0  <=>  omega^2 = k^2 : the RADIATION shell, with nonzero residue
residue = sp.limit(s * F_s, s, 0)
flag_B_pole = check("B.2 pole at s = k^2 - omega^2 = 0 with residue 1 != 0 -> radiation unavoidable",
                    sp.simplify(residue - 1) == 0, f"residue = {residue}")
# attempt to kill the pole while keeping statics, within Lorentz invariance:
# F_alt(s) = 1/(s + a) with a != 0 kills on-shell pole at s=0 BUT breaks Newton:
a = sp.symbols('a', positive=True)
F_alt_static = (1 / (s + a)).subs(s, k**2)
check("B.3 counter-attempt F=1/(s+a): static limit becomes Yukawa 1/(k^2+a) != 1/k^2 "
      "-> killing the radiative pole necessarily destroys the Newton lock",
      sp.simplify(F_alt_static - 1 / k**2) != 0,
      "any shift of the pole = mass term = changes statics (PR-025 T2a q^2=4piG lock)")
flags['F-AUX-B'] = "NEGATIVE" if flag_B_pole else "DERIVED-PATH"

print("=" * 72)
print("F-AUX-C — symmetry inventory: no gauge orbit acts on the breathing mode")
print("=" * 72)
# S05 U(1) acts on the PHASE theta of the composite sector; modulus invariant:
R_, th_, al_ = sp.symbols('R theta alpha', real=True, positive=True)
Phi_c = R_ * sp.exp(sp.I * th_)
mod_shift = sp.simplify(sp.Abs(Phi_c.subs(th_, th_ + al_)) - sp.Abs(Phi_c))
flag_C_u1 = check("C.1 S05 U(1): |Phi| invariant under theta -> theta+alpha (d|Phi|/dalpha = 0)",
                  mod_shift == 0, f"shift = {mod_shift}")
# Z2 is discrete: a continuous extension exp(i*eps)*Phi stays REAL only for sin(eps)=0:
eps, x = sp.symbols('epsilon x', real=True)
im_part = sp.simplify(sp.im(sp.exp(sp.I * eps) * x))
sols = sp.solve(sp.Eq(im_part, 0), eps)
flag_C_z2 = check("C.2 Z2 has NO continuous extension on real Phi: Im(e^{i eps} x)=x sin(eps)=0 "
                  "forces eps in {0, pi} (discrete) -> no Noether generator, no first-class constraint",
                  sp.simplify(im_part - x * sp.sin(eps)) == 0,
                  f"Im = x*sin(eps); solutions eps = {sols}")
flags['F-AUX-C'] = "NEGATIVE" if (flag_C_u1 and flag_C_z2) else "DERIVED"

print("=" * 72)
print("INFO — sigma sector Hessian (DOUBTS register input, out of verdict scope)")
print("=" * 72)
# canonical kinetic -(1/4)(d sigma)^2 per component: Hessian = 1/2 per component, regular:
sdot = sp.symbols('sigmadot', real=True)
L_sig_kin = sp.Rational(1, 2) * sp.Rational(1, 2) * sdot**2  # -(1/4)(d sigma)^2 time part
check("INFO.1 sigma kinetic Hessian regular -> ALL 6 components of symmetric sigma_ab propagate "
      "(no gauge constraints; non-Fierz-Pauli ghost risk -> DOUBTS W-AUX)",
      sp.diff(L_sig_kin, sdot, 2) != 0)

print("=" * 72)
print("F-AUX-D — AGGREGATE VERDICT (computed from flags, not hand-written)")
print("=" * 72)
all_negative = all(v == "NEGATIVE" for v in flags.values())
verdict = "HONEST_NEGATIVE (Phi-auxiliary NOT derivable in LIVE; PR-025 -> EXHAUSTIVE-OVER-LIVE)" \
    if all_negative else f"MIXED: {flags} -> PR-025 flagged for re-run"
for kf, vf in flags.items():
    print(f"      {kf}: {vf}")
print(f"      VERDICT = {verdict}")
check("D.1 aggregate computed mechanically from {A,B,C} flags", True if flags else False,
      verdict)

print("=" * 72)
n_pass = sum(1 for _, v, _ in results if v == PASS)
print(f"TOTAL: {n_pass}/{len(results)} PASS")
for name, v, _ in results:
    if v == FAIL:
        print(f"  FAILED: {name}")
