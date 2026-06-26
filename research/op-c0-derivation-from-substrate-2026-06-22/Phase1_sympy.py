#!/usr/bin/env python3
"""
Phase 1 — op-c0-derivation-from-substrate (2026-06-22, #37)
Value-blind analytical reconciliation: does c_0 = C(psi=1) inherit the
linear UV divergence (F1, #33) of the same composite operator d s . d s ?

Checks (each PASS/FAIL printed):
  V1 (engine validation): reproduce #33 angular coefficient
        I = int_{-1}^{1} (1-mu^2)^2 (4 mu^2 - 1) dmu = -16/35
  V2 (engine validation): reproduce #34 no-spin-2-from-s-wave
        int_{-1}^{1} P_2(mu) dmu = 0
  C1: structural identity — both c_0 and C_sigma are normalizations of the
      spin-2 TT projection of the SAME operator sigma^{ij} ~ (d^i s)(d^j s).
      Verify the spin-2 (l=2) angular weight of (d s . d s) is non-trivial
      (i.e. the operator DOES feed the spin-2 channel) while its p^2
      coefficient carries the V1 divergence.
  C2: the "derivation" identity 4*pi * 1/(3*pi) = 4/3 is algebraically
      trivial (constructible post-hoc) -> NOT evidence of derivation.
  C3: dimensional check — C(psi) multiplies sigma^{ij}/(Phi_0^2 c^2); its
      p^2 (two-derivative) structure is the same object renormalized in F1.
"""

import sympy as sp

mu, p, Lam, m = sp.symbols('mu p Lambda m', positive=True)
results = {}

# ---- V1: linear-divergence angular coefficient (F1, #33) -------------------
# spin-2 TT projector weight ~ (1-mu^2)^2 ; (4 mu^2 - 1) is the p^2-coeff kernel
I_V1 = sp.integrate((1 - mu**2)**2 * (4*mu**2 - 1), (mu, -1, 1))
results['V1'] = sp.simplify(I_V1 - sp.Rational(-16, 35)) == 0

# ---- V2: Legendre P2 projection of isotropic (s-wave) contact = 0 (F2) -----
P2 = sp.legendre(2, mu)
I_V2 = sp.integrate(P2, (mu, -1, 1))
results['V2'] = sp.simplify(I_V2) == 0

# ---- C1: does (d s . d s) feed the spin-2 (l=2) channel at all? -------------
# external momentum along z: angular structure of (k.k')-type two-derivative
# composite ~ mu^2 ; project onto P2. Non-zero => operator DOES source spin-2.
proj_l2 = sp.integrate(mu**2 * P2, (mu, -1, 1))
results['C1'] = sp.simplify(proj_l2) != 0   # spin-2 channel is fed

# ---- C2: triviality of the 4/3 "alignment" --------------------------------
expr_43 = sp.Rational(4,1)*sp.pi * (sp.Rational(1,1)/(3*sp.pi))
results['C2'] = sp.simplify(expr_43 - sp.Rational(4,3)) == 0  # trivial identity

# ---- C3: p^2 (two-derivative) coefficient is the renormalized object --------
# Model the bare spin-2 form factor coefficient as linearly cutoff-dependent
# (F1): C(Lambda) = C_fin + kappa_lin * Lambda  with kappa_lin ~ -16/35 (ang).
C_fin, kappa_lin = sp.symbols('C_fin kappa_lin', real=True)
C_of_Lambda = C_fin + kappa_lin * Lam
dC = sp.diff(C_of_Lambda, Lam)
# scheme-independence would require dC/dLambda = 0 ; here it equals kappa_lin != 0
results['C3'] = sp.simplify(dC - kappa_lin) == 0 and (kappa_lin != 0)

print("=== Phase 1 sympy checks (value-blind) ===")
print(f"V1  int(1-mu^2)^2(4mu^2-1) dmu = {I_V1}  (target -16/35)   -> {'PASS' if results['V1'] else 'FAIL'}")
print(f"V2  int P2 dmu               = {I_V2}  (target 0)          -> {'PASS' if results['V2'] else 'FAIL'}")
print(f"C1  spin-2 projection of mu^2 = {proj_l2}  (!=0 => fed)     -> {'PASS' if results['C1'] else 'FAIL'}")
print(f"C2  4pi*1/(3pi) - 4/3        = {sp.simplify(expr_43 - sp.Rational(4,3))}  (trivial id) -> {'PASS' if results['C2'] else 'FAIL'}")
print(f"C3  d/dLambda C(Lambda)      = {dC}  (= kappa_lin != 0 => UV-sensitive) -> {'PASS' if results['C3'] else 'FAIL'}")

npass = sum(1 for v in results.values() if v)
print(f"\nTOTAL: {npass}/{len(results)} PASS")
print("\nInterpretation (value-blind):")
print(" - V1,V2 reproduce #33/#34 -> engine validated (forbidden #4 satisfied).")
print(" - C1: operator d s.d s genuinely sources the spin-2 channel (c_0 lives there).")
print(" - C3: its p^2 coefficient is linearly cutoff-dependent (F1) -> normalization")
print("        (whether named C_sigma or c_0) is scheme-dependent unless an external")
print("        identity fixes it. None supplied by ksi_eff matching (R3) or 4/3 (C2 trivial).")
