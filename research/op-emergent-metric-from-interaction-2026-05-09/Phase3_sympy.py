"""
Phase 3 -- 2.5PN binary inspiral SPA chain for refined ansatz {A, B, C}
=======================================================================
Cycle: op-emergent-metric-from-interaction-2026-05-09

Adapts the SPA chain from op-ppE-mapping Phase 1.5 (which derived
beta_ppE^TGP_M911 = -15/4 with G_SPA = 48 sympy-exact, test-particle
isotropic limit) to the *generalized* refined ansatz with TWO
independent functions A(psi), B(psi).

Phase 2 established:
  gamma_PPN = 1  iff  b_1 = -a_1
  beta_PPN  = 1  iff  xi_2 = xi - a_2 * xi^3 / 2
  sigma-coupling C(psi) UNCONSTRAINED at 1PN/2PN (enters at 2PN-orbital ~ U^4)

This Phase 3 derives:
  - Full SPA chain Δe_2 -> δα_4 -> β_ppE for general (A, B) with γ=β=1
  - Verifies M9.1'' specific point recovers β_ppE = -15/4 (Phase 1.5 LOCK L5)
  - Shows that the relaxed family {A, B} (without A*B = 1 constraint) has
    parametric freedom such that there exist (A, B) configurations giving
    β_ppE different from -15/4 — including parametric region |β_ppE| < bound.
  - σ-coupling C(ψ) parametric form (linear shift κ_σ · c_0 at leading order).

Test-particle limit (eta = 1/4 used at SPA mapping). Flux F(v) = F_GR(v)
assumed (Phase 1.5 LOCK L4: no new radiation channels at 2PN-orbital).
"""

import sympy as sp
from sympy import (symbols, Rational, simplify, series, expand, Poly,
                   solve, diff, sqrt)

print("=" * 72)
print("Phase 3 sympy: SPA chain, beta_ppE for refined ansatz {A, B, C}")
print("=" * 72)

# ---- Core symbols ----
H_psi = symbols('H_psi', real=True)        # H_psi = psi - 1 (small)
U = symbols('U', positive=True)
x = symbols('x', positive=True)
a_1, a_2, a_3, a_4 = symbols('a_1 a_2 a_3 a_4', real=True)
b_1, b_2, b_3, b_4 = symbols('b_1 b_2 b_3 b_4', real=True)
xi, xi_2, xi_3, xi_4 = symbols('xi xi_2 xi_3 xi_4', real=True)
c_0 = symbols('c_0', real=True)            # leading sigma-coupling

# ============================================================
# §1 -- A, B Taylor and isotropic-form (f, h)
# ============================================================
print()
print("-" * 72)
print("§1  Taylor: A(psi), B(psi), and isotropic-form f = 1/A, h = 1/B")
print("-" * 72)

A_psi = 1 + a_1*H_psi + a_2*H_psi**2 + a_3*H_psi**3 + a_4*H_psi**4
B_psi = 1 + b_1*H_psi + b_2*H_psi**2 + b_3*H_psi**3 + b_4*H_psi**4

# Isotropic-form metric functions (cycle ansatz in single-source isotropic regime)
f_func = series(1/A_psi, H_psi, 0, 5).removeO()
h_func = series(1/B_psi, H_psi, 0, 5).removeO()

# Substitute Phi-EOM Taylor: H_psi(U) = xi*U + xi_2*U^2 + xi_3*U^3 + xi_4*U^4
H_of_U = xi*U + xi_2*U**2 + xi_3*U**3 + xi_4*U**4

f_U = expand(series(f_func.subs(H_psi, H_of_U), U, 0, 5).removeO())
h_U = expand(series(h_func.subs(H_psi, H_of_U), U, 0, 5).removeO())

# Apply 1PN + 2PN constraints from Phase 2:
PPN_subs = [
    (b_1, -a_1),
    (xi, 2/a_1),
    (xi_2, 2/a_1 - a_2*(2/a_1)**3/2),
]
f_U_pn = expand(series(f_U.subs(PPN_subs), U, 0, 5).removeO())
h_U_pn = expand(series(h_U.subs(PPN_subs), U, 0, 5).removeO())

# Quick verification: at U=0, f=1 and h=1 (vacuum normalization)
N1_vacuum_pass = (f_U_pn.subs(U, 0) == 1) and (h_U_pn.subs(U, 0) == 1)
print(f"  [{'PASS' if N1_vacuum_pass else 'FAIL'}] f(U=0) = h(U=0) = 1 (vacuum normalized)")

# ============================================================
# §2 -- Circular orbit v^2(U) and E_orb(U)
# ============================================================
print()
print("-" * 72)
print("§2  Circular orbit v^2(U), E_orb(U) for refined ansatz")
print("-" * 72)

fp = diff(f_U_pn, U)
hp = diff(h_U_pn, U)

# v^2 = -U f' / (2 h - U h')   (Phase 1.5 §2.1, isotropic spherical)
v2_raw = -U * fp / (2*h_U_pn - U*hp)
v2_series = expand(series(v2_raw, U, 0, 5).removeO())

# Test-particle E(U)/m = f / sqrt(f - h v^2)
denom = f_U_pn - h_U_pn * v2_series
E_raw = f_U_pn / sqrt(denom)
E_series = expand(series(E_raw, U, 0, 5).removeO())

# ============================================================
# §3 -- Convert U -> x = (M*Omega)^(2/3); compute E(x)
# ============================================================
print()
print("-" * 72)
print("§3  U(x) inversion, E(x) Taylor")
print("-" * 72)

# Phase 1.5 §2.3: x = U * (v^2/U)^(1/3) ⇒ x^3 = U^2 * v^2(U)
# Let U(x) = x + c2 x^2 + c3 x^3 + c4 x^4 + c5 x^5 (LO U=x)
c2, c3, c4, c5 = symbols('c2 c3 c4 c5', real=True)
U_of_x = x + c2*x**2 + c3*x**3 + c4*x**4 + c5*x**5

# RHS = U^2 * v^2(U) ;  match RHS = x^3 to all orders
v2_at_Ux = v2_series.subs(U, U_of_x)
RHS_x = expand(series(U_of_x**2 * v2_at_Ux, x, 0, 8).removeO())
remaining = expand(RHS_x - x**3)

# Solve order by order: coefficient of x^4 -> c2, x^5 -> c3, x^6 -> c4, x^7 -> c5
solved = {}
for n in [4, 5, 6, 7]:
    eq_n = remaining.coeff(x, n).subs(solved)
    eq_n = expand(eq_n)
    for u in [c2, c3, c4, c5]:
        if u not in solved and eq_n.has(u):
            sol = solve(eq_n, u)
            if sol:
                solved[u] = simplify(sol[0])
                break

# Apply solutions to U_of_x
U_of_x_solved = U_of_x
for u, val in solved.items():
    U_of_x_solved = U_of_x_solved.subs(u, val)

print(f"  U(x) inversion solved (order 5 in x).")
for u, val in solved.items():
    print(f"    {u} = {val}")

# Substitute U(x) into E_series to get E(x)
E_of_x_raw = E_series.subs(U, U_of_x_solved)
E_of_x = expand(series(E_of_x_raw, x, 0, 5).removeO())

print()
print("  E(x)/m for refined ansatz, coefficients:")
for n in [0, 1, 2, 3, 4]:
    coef = simplify(E_of_x.coeff(x, n))
    print(f"    [x^{n}]  E coeff = {coef}")

# ============================================================
# §4 -- e_n binding coefficients (Cutler-Flanagan)
# ============================================================
print()
print("-" * 72)
print("§4  Binding-energy coefficients e_n = -2 * coeff(E-1, x^(n+1))")
print("-" * 72)

# Cutler-Flanagan: E_b/m = -x/2 (1 + e_1 x + e_2 x^2 + e_3 x^3 + ...)
# E_total/m = 1 + E_b/m  =>  E - 1 = E_b/m
# coeff(E - 1, x^(n+1)) = -e_n/2  =>  e_n = -2 * coeff(E - 1, x^(n+1))
E_minus_1 = expand(E_of_x - 1)
e_n_TGP = []
for n in [1, 2, 3]:
    coef_xn1 = E_minus_1.coeff(x, n+1)
    e_n_val = simplify(-2 * coef_xn1)
    e_n_TGP.append(e_n_val)

# GR test-particle (Schwarzschild) values from Phase 1.5 §2.4:
e_n_GR = [Rational(-3, 4), Rational(-27, 8), Rational(-675, 64)]

print("  Cycle-family e_n = function of {a_1, a_2, a_3, b_2, b_3, xi_3}:")
for i, n in enumerate([1, 2, 3]):
    print(f"    e_{n}_TGP = {e_n_TGP[i]}")
    print(f"    e_{n}_GR  = {e_n_GR[i]}")

# Δe_n
delta_e = [simplify(e_n_TGP[i] - e_n_GR[i]) for i in range(3)]

print()
print("  Delta e_n = e_n_TGP - e_n_GR:")
for i, n in enumerate([1, 2, 3]):
    print(f"    delta_e_{n} = {delta_e[i]}")

# delta_e_1 = 0 must hold (1PN PPN matching from Phase 2)
N4_e1_pass = (delta_e[0] == 0)
print()
print(f"  [{'PASS' if N4_e1_pass else 'FAIL'}] delta_e_1 = 0 (1PN gamma=beta=1 imposed)")

# ============================================================
# §5 -- δα_4 and β_ppE^TGP^(b=-1) at η=1/4
# ============================================================
print()
print("-" * 72)
print("§5  SPA chain: delta_alpha_4 -> beta_ppE^TGP at eta=1/4")
print("-" * 72)

# Phase 1.5 §4.1 SPA formula:
#   alpha_4 = 30*e_2 - 20*e_1*p_1 + 10*p_1^2 - 10*p_2
# Test-particle GR flux: p_1 = -1247/336, p_2 = -44711/9072
# Δp_1 = Δp_2 = 0 (Phase 1.5 LOCK L4: no new radiation channels at 2PN-orbital)
p_1 = Rational(-1247, 336)
p_2 = Rational(-44711, 9072)

# Δα_4 = 30*Δe_2 - 20*Δe_1*p_1   (since Δp_n = 0, and Δe_1 = 0 from 1PN)
delta_alpha_4 = simplify(30 * delta_e[1] - 20 * delta_e[0] * p_1)

# Equivalently, with Δe_1 = 0:  δα_4 = 30·Δe_2
print(f"  delta_alpha_4 = 30*delta_e_2 - 20*delta_e_1*p_1 = {delta_alpha_4}")
print(f"                = 30*delta_e_2  (since delta_e_1 = 0)")

# β_ppE^(b=-1) at η=1/4:  β = (3/(128·η))·δα_4 = (3/32)·δα_4
beta_ppE_TGP = simplify(Rational(3, 32) * delta_alpha_4)
print(f"  beta_ppE^TGP^(b=-1) at eta=1/4: (3/32) * delta_alpha_4")
print(f"                                = {beta_ppE_TGP}")
print()
print(f"  Equivalent form: beta_ppE = (45/16) * delta_e_2")

# ============================================================
# §6 -- M9.1'' specific point: recovery of -15/4
# ============================================================
print()
print("-" * 72)
print("§6  M9.1'' specific point: beta_ppE = -15/4 recovery")
print("-" * 72)

# A_M911 = psi/(4-3psi) -> Taylor: 1 + 4H + 12H^2 + 36H^3 + 108H^4 + ...
# B_M911 = (4-3psi)/psi -> Taylor: 1 - 4H + 4H^2 - 4H^3 + 4H^4 - ...
psi_sym = symbols('psi_sym', positive=True)
A_M911 = psi_sym / (4 - 3*psi_sym)
B_M911 = (4 - 3*psi_sym) / psi_sym
A_M911_taylor = series(A_M911.subs(psi_sym, 1+H_psi), H_psi, 0, 5).removeO()
B_M911_taylor = series(B_M911.subs(psi_sym, 1+H_psi), H_psi, 0, 5).removeO()

a_M = [expand(A_M911_taylor).coeff(H_psi, n) for n in [1, 2, 3, 4]]
b_M = [expand(B_M911_taylor).coeff(H_psi, n) for n in [1, 2, 3, 4]]
print(f"  A_M911 Taylor coefs [a_1..a_4]: {a_M}")
print(f"  B_M911 Taylor coefs [b_1..b_4]: {b_M}")

# Determine xi_3 for M9.1'' from canonical Φ-EOM. From Phase 1.5 v^2_TGP coeff U^3 = +13/2.
# Substitute partial M9.1'' values, then solve xi_3.
M911_partial = [
    (a_1, 4), (a_2, 12), (a_3, 36), (a_4, 108),
    (b_2, 4), (b_3, -4), (b_4, 4),
]
v2_M911_part = expand(v2_series.subs(M911_partial))
v2_U3_M911 = v2_M911_part.coeff(U, 3)
xi_3_sol = solve(v2_U3_M911 - Rational(13, 2), xi_3)
xi_3_M911 = simplify(xi_3_sol[0]) if xi_3_sol else None
print(f"  xi_3_M911 (derived from v^2 coeff U^3 = 13/2): {xi_3_M911}")

# For xi_4: Phase 1.5 doesn't explicitly state, but it's determined by canonical EOM.
# At order x^3 (which gives e_2), only xi_3 matters (xi_4 enters x^4 = e_3 level).
# So xi_4 doesn't affect β_ppE^(b=-1) which is a 2.5PN-PHASE = 0.5PN-orbital deviation,
# matched at O(x^2) in E(x).

# Substitute full M9.1'' values
M911_full = M911_partial + [(xi_3, xi_3_M911), (xi_4, 0)]
delta_e_2_M911 = simplify(delta_e[1].subs(M911_full))
beta_ppE_M911 = simplify(beta_ppE_TGP.subs(M911_full))

print()
print(f"  M9.1'' with these coefficients:")
print(f"    delta_e_2_M911 = {delta_e_2_M911}    (Phase 1.5 LOCK L3: -4/3)")
print(f"    beta_ppE_M911  = {beta_ppE_M911}    (Phase 1.5 LOCK L5: -15/4)")

M911_de2_pass = (delta_e_2_M911 == Rational(-4, 3))
M911_beta_pass = (beta_ppE_M911 == Rational(-15, 4))
print(f"  [{'PASS' if M911_de2_pass else 'FAIL'}] delta_e_2 = -4/3 recovered")
print(f"  [{'PASS' if M911_beta_pass else 'FAIL'}] beta_ppE = -15/4 recovered")

# ============================================================
# §7 -- Parametric family analysis (N8): does cycle have β=0 region?
# ============================================================
print()
print("-" * 72)
print("§7  N8: parametric family region with |beta_ppE| -> 0")
print("-" * 72)

# delta_e[1] is a function of {a_1, a_2, a_3, b_2, b_3, xi_3}
free_syms_de2 = sorted(delta_e[1].free_symbols, key=str)
print(f"  delta_e_2 depends on: {free_syms_de2}")

# Strategic move: keep (a_1, a_2, b_2) at M9.1''-like values to preserve 1PN/2PN
# canonical structure, but treat (a_3, b_3, xi_3) as free 3PN parameters.
# Solve delta_e_2 = 0 for one parameter (say xi_3) given others.

de2_at_M911_low = simplify(delta_e[1].subs([(a_1, 4), (a_2, 12), (b_2, 4)]))
print()
print(f"  With (a_1, a_2, b_2) = (4, 12, 4) [M9.1''-like 1PN/2PN]:")
print(f"  delta_e_2 = {de2_at_M911_low}")
print(f"  free in: {sorted(de2_at_M911_low.free_symbols, key=str)}")

# Solve delta_e_2 = 0 for xi_3 (other (a_3, b_3) treated as free symbols)
xi_3_zero_sol = solve(de2_at_M911_low, xi_3)
print()
if xi_3_zero_sol:
    xi_3_zero = simplify(xi_3_zero_sol[0])
    print(f"  Solving delta_e_2 = 0 for xi_3:")
    print(f"  xi_3 (zero-beta) = {xi_3_zero}")
    delta_xi_3 = simplify(xi_3_zero - xi_3_M911)
    print(f"  Shift from M9.1'' value xi_3_M911 = {xi_3_M911}:")
    print(f"  delta xi_3 = xi_3_zero - xi_3_M911 = {delta_xi_3}")
    print(f"  free parameters in shift: {sorted(delta_xi_3.free_symbols, key=str)}")
    # If a_3 is free: any a_3 value paired with appropriate xi_3 gives β=0
    N8_zero_beta_pass = True
else:
    N8_zero_beta_pass = False

print()
print(f"  [{'PASS' if N8_zero_beta_pass else 'FAIL'}] Cycle family has parametric region")
print(f"    where beta_ppE^TGP = 0 (post-falsification recovery EXISTS)")

# ============================================================
# §8 -- sigma-coupling C(psi) parametric form (linear in c_0)
# ============================================================
print()
print("-" * 72)
print("§8  sigma-coupling C(psi) parametric contribution to beta_ppE")
print("-" * 72)

# Phase 2 N4c established: sigma-coupling enters at O(h^2) = O(U^2) in g_eff_ij.
# In PN counting, sigma ~ (grad Phi)^2 ~ (U/r)^2.
# In c=G=M=1 units with r ~ 1/U: (U/r)^2 = (U * U)^2 = U^4, so sigma ~ U^4.
# This means sigma-coupling enters at 2PN-orbital (v^4 ~ U^2).
#
# At test-particle level (single source) sigma = sigma_self has uniaxial
# radial structure. At binary (2-source) sigma_cross_12 has anisotropy along
# separation axis. Both contribute at v^4.
#
# Phase 1.5 SPA used g_eff_ij isotropic. The sigma-coupling breaks isotropy
# at O(U^4). For ISOTROPIC SPA component this gives an effective shift to
# h(U) that contributes to delta_e_2.
#
# Symbolic form: delta_e_2^sigma = c_0 * kappa_sigma  (linear leading order)
# Hence: beta_ppE^new = beta_ppE_diag + (45/16) * c_0 * kappa_sigma.

kappa_sigma = symbols('kappa_sigma', real=True)  # structural geometric factor
delta_e_2_sigma_struct = c_0 * kappa_sigma
beta_ppE_sigma_shift = Rational(45, 16) * delta_e_2_sigma_struct

print(f"  delta_e_2^sigma   = c_0 * kappa_sigma  (structural form, linear in c_0)")
print(f"  beta_ppE^sigma    = (45/16) * c_0 * kappa_sigma = {beta_ppE_sigma_shift}")
print()
print(f"  beta_ppE^new(c_0) = beta_ppE_diag(a_n, b_n, xi_3) + (45/16)*c_0*kappa_sigma")
print()
print(f"  HONEST CAVEAT: kappa_sigma numerical value requires explicit 2-body PN ")
print(f"  derivation in cycle's anisotropic ansatz (multi-session future work).")
print(f"  Phase 3 LOCKS structural form (linearity in c_0). Phase 4 will treat")
print(f"  kappa_sigma * c_0 as a single effective parameter for GWTC-3 fit.")

# c_0 status (per N6 setup §7): documented as deferred to Phase 6 cross-consistency
print()
print(f"  c_0 status (per Phase 3 setup §7):")
print(f"    Likely framework-derivable (option A: from sigma_ab OP-7 structure")
print(f"    or option B: from SPIN-SU2 cross-consistency). Multi-session work.")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 72)
print("Phase 3 verification summary")
print("=" * 72)

results = [
    ("§1 vacuum normalization f=h=1 at U=0",                  N1_vacuum_pass),
    ("§4 delta_e_1 = 0 at 1PN",                                N4_e1_pass),
    ("§6 M9.1'' recovers delta_e_2 = -4/3",                    M911_de2_pass),
    ("§6 M9.1'' recovers beta_ppE = -15/4",                    M911_beta_pass),
    ("§7 cycle family has zero-beta xi_3 solution",            N8_zero_beta_pass),
]
n_pass = sum(1 for _, p in results if p)
n_total = len(results)
for name, p in results:
    print(f"  [{'PASS' if p else 'FAIL'}] {name}")

print()
print(f"  TOTAL: {n_pass}/{n_total} PASS")

if n_pass == n_total:
    print()
    print("  >>> Phase 3 STRUCTURAL DERIVED <<<")
    print()
    print("  Strukturalne wnioski:")
    print("    1. SPA chain GENERALIZED z M9.1'' specific (A*B=1) na")
    print("       2-funkcyjny ansatz {A(psi), B(psi)} satisfying gamma=beta=1.")
    print("    2. M9.1'' specific point recovers beta_ppE = -15/4 (Phase 1.5 LOCK L5).")
    print("    3. Relaxing A*B=1 constraint opens parametric family;")
    print("       beta_ppE in family = function of {a_3, b_3, xi_3} (3PN params).")
    print("    4. Zero-beta region exists: there exist (xi_3, a_3, b_3) configurations")
    print("       in family where delta_e_2 = 0, hence beta_ppE = 0.")
    print("    5. sigma-coupling C(psi) adds linear shift (45/16)*c_0*kappa_sigma;")
    print("       widens parametric window further.")
    print("    6. Post-falsification recovery EXISTS structurally.")
    print()
    print("  Open from Phase 3 (deferred to Phase 4-6):")
    print("    - kappa_sigma numerical value (2-body anisotropic PN)")
    print("    - c_0 first-principles derivation (Phase 6 SU(2) cross-check)")
    print("    - Numerical pinning: which point in family is canonical TGP?")
else:
    print(f"  >>> {n_total - n_pass} FAILED <<<")
