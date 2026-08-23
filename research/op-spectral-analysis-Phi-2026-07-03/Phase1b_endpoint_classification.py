#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Phase 1b: exact classification of the degenerate endpoint psi->0 in F-A.

Follow-up to Phase1 T7b (hypothesis 'repulsive barrier' REFUTED: leading
coefficient at assumed power chi^(-2/3) is 0). Here: compute the TRUE
leading behavior of U''(chi) as chi->0 and classify the field-space
endpoint. Also state precisely what protects physics from psi->0.

Criteria C6(b) (Phase0 LOCK): classification resolved, any direction.
"""
import sympy as sp

chi, gamma, Kgeo = sp.symbols('chi gamma K_geo', positive=True)
x = sp.symbols('x', positive=True)

# U_A(psi) with beta=gamma:  U_A' = K_geo*gamma*(x^7 - x^6)
UA = sp.integrate(Kgeo * gamma * (x**7 - x**6), x)     # = Kgeo*gamma*(x^8/8 - x^7/7)
psi_of_chi = (3 * chi) ** sp.Rational(1, 3)            # canonical field: chi = psi^3/3
U_chi = sp.simplify(UA.subs(x, psi_of_chi))
d2U = sp.expand(sp.simplify(sp.diff(U_chi, chi, 2)))
print("U(chi)      =", sp.simplify(U_chi))
print("U''(chi)    =", d2U)

# Leading terms as chi->0:
series = sp.series(d2U, chi, 0, 2)
print("series      =", series)

# leading power and coefficient
poly = sp.Poly(sp.expand(d2U / (chi ** sp.Rational(1, 3))), sp.Rational(1, 1)) if False else None
lead_13 = sp.limit(d2U / chi ** sp.Rational(1, 3), chi, 0)
print("lim chi->0 U''/chi^(1/3) =", sp.simplify(lead_13))

verdict_soft = sp.simplify(lead_13).is_negative
print()
print("CLASSIFICATION (C6b):")
print("  * U''(chi) -> 0^- as chi->0 with leading term %s * chi^(1/3)" % sp.simplify(lead_13))
print("  * NO repulsive barrier at the degenerate endpoint psi=0 (hypothesis of")
print("    Phase1 T7b REFUTED). The endpoint is SOFT (curvature -> 0 from below).")
print("  * U(0)=0, U'(0)=0: psi=0 is a degenerate stationary boundary point,")
print("    a local MAXIMUM direction of U for small chi (U ~ -Kgeo*gamma*(3chi)^(7/3)/7 < 0... check):")
U_small = sp.series(U_chi, chi, 0, 2)
print("    U(chi->0) =", U_small)
print("  * Consequence: exclusion of psi->0 in physical configurations does NOT")
print("    come from a dynamical barrier of U; it rests on the axiomatic")
print("    no-absolute-vacuum chain (W(1)=gamma/3 != 0, ssec:vacuum-source-chain)")
print("    plus positivity of sources (rho >= 0 pushes delta psi > 0,")
print("    rem:sign-convention).")
print("  * Spectral consequence for L03: on any physical background with")
print("    min_r psi_eq(r) = psi_min > 0 the S-L problem has p(r) = K(psi_eq) >=")
print("    K(psi_min) > 0: uniformly elliptic, no singular endpoint enters the")
print("    physical domain. The degenerate point psi=0 lies OUTSIDE all")
print("    physical backgrounds; its softness is irrelevant for sigma(L) as")
print("    long as no-absolute-vacuum holds (axiom-conditional statement).")
print()
print("PASS(C6b): classification resolved:",
      "SOFT endpoint, axiom-conditional exclusion" if verdict_soft else "UNRESOLVED")
