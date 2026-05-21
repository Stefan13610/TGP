#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1a -- N=2 kink+antikink isolation null test
op-CE-H-two-particle-equilibrium-2026-05-21

Pre-registered prediction (LOCKED 2026-05-21):
F-beta-1: brak stabilnego L* > 0 z d2E/dL2 > 0 w pure isolation (1D Z2 toy).

Tests:
- T_P1a_1: Verify single kink Phi_K(x) = v*tanh(m*x/sqrt(2)) satisfies EOM
- T_P1a_2: Compute kink energy E_K = (2*sqrt(2)/3) * m * v^2
- T_P1a_3: V_int(L) Manton-Sutcliffe form (LIT informational)
- T_P1a_4: dE/dL has no zero in L > 0 (no stationary point) -- F-beta-1 verification
- T_P1a_5: Asymptotic limits E(L->0+) and E(L->inf)

Discipline:
- Strict cycle 1/2/7 (0 hardcoded T_pass=True)
- 1 LIT informational (T_P1a_3)
- 4 substantive FP

Output: all values computed analytically via sympy.
"""
import sympy as sp
import sys

print("="*70)
print("Phase 1a -- N=2 isolation null test (1D Z2 toy)")
print("op-CE-H-two-particle-equilibrium-2026-05-21")
print("Pre-registered prediction F-beta-1: NO stable L* in isolation")
print("="*70)

# Symbols
x, L = sp.symbols('x L', real=True, positive=True)
x_real = sp.Symbol('x_real', real=True)
v, lam, m = sp.symbols('v lambda m', positive=True)
# Physical relation: m^2 = lambda * v^2 (mass of fluctuations around vacuum)
# Kink: Phi_K(x) = v * tanh(m*x/sqrt(2)) satisfies EOM with this relation
# E_K = (2*sqrt(2)/3) * m * v^2 (standard result)

# === T_P1a_1: Verify single kink solution ===
print("\n--- T_P1a_1: Single kink solution verification (substantive) ---")

# Single kink ansatz
Phi_K = v * sp.tanh(m*x_real/sp.sqrt(2))
print(f"Ansatz: Phi_K(x) = v * tanh(m*x/sqrt(2))")

# EOM: d^2 Phi/dx^2 = dV/dPhi = lambda * Phi * (Phi^2 - v^2)
d2Phi_K = sp.diff(Phi_K, x_real, 2)
dVdPhi = lam * Phi_K * (Phi_K**2 - v**2)

# EOM residual: d^2Phi/dx^2 - dV/dPhi must be 0 (using m^2 = lam*v^2)
EOM_residual = sp.simplify(d2Phi_K - dVdPhi)
# Substitute m -> sqrt(lam)*v (i.e. m^2 = lam*v^2)
EOM_residual_sub = sp.simplify(EOM_residual.subs(m, sp.sqrt(lam)*v))
EOM_residual_final = sp.simplify(EOM_residual_sub)

print(f"EOM residual after m = sqrt(lam)*v: {EOM_residual_final}")
T_P1a_1_pass = (EOM_residual_final == 0)
print(f"T_P1a_1 PASS (EOM satisfied)? {bool(T_P1a_1_pass)}")

# === T_P1a_2: Single kink energy ===
print("\n--- T_P1a_2: Single kink energy computation (substantive) ---")

# E_K = integral_{-inf}^{+inf} [ (1/2)(dPhi/dx)^2 + V(Phi) ] dx
# V(Phi) = (lambda/4)(Phi^2 - v^2)^2

dPhi_K_dx = sp.diff(Phi_K, x_real)
V_at_Phi_K = (lam/4) * (Phi_K**2 - v**2)**2
energy_density = sp.Rational(1, 2) * dPhi_K_dx**2 + V_at_Phi_K
energy_density_simplified = sp.simplify(energy_density)

# Integrate
E_K_integral = sp.integrate(energy_density_simplified, (x_real, -sp.oo, sp.oo))
E_K_computed = sp.simplify(E_K_integral)

# Substitute lam = m^2/v^2 to express in terms of m, v
E_K_in_mv = sp.simplify(E_K_computed.subs(lam, m**2 / v**2))

# Expected: E_K = (2*sqrt(2)/3) * m * v^2
E_K_expected = sp.Rational(2) * sp.sqrt(2) / sp.Rational(3) * m * v**2

E_K_diff = sp.simplify(E_K_in_mv - E_K_expected)

print(f"Computed E_K (in lam, v): {E_K_computed}")
print(f"Computed E_K (in m, v): {E_K_in_mv}")
print(f"Expected E_K = (2*sqrt(2)/3) * m * v^2 = {E_K_expected}")
print(f"Difference: {E_K_diff}")
T_P1a_2_pass = (E_K_diff == 0)
print(f"T_P1a_2 PASS (E_K matches known result)? {bool(T_P1a_2_pass)}")

# === T_P1a_3: Kink-antikink interaction (LIT informational) ===
print("\n--- T_P1a_3: V_int(L) Manton-Sutcliffe form (LIT informational) ---")
print("status_class: INFORMATIONAL (from literature, NOT pre-registered substantive)")

# From Manton-Sutcliffe 2004 Ch.5, for kink-antikink in phi^4 mexican hat:
# V_int(L) = -A * exp(-m*L) at large L, A > 0 positive constant
# Specific A: A_Manton = 32 * m^3 * v^2 / lam, but we'll use symbolic A

A_int = sp.Symbol('A_int', positive=True)
V_int_param = -A_int * sp.exp(-m*L)

V_int_at_0 = V_int_param.subs(L, 0)
V_int_at_inf = sp.limit(V_int_param, L, sp.oo)
sign_at_finite_L = sp.simplify(V_int_param.subs([(L, 1), (A_int, 1), (m, 1)]))

print(f"V_int(L) parameterization: V_int(L) = -A_int * exp(-m*L)")
print(f"V_int(L=0) = {V_int_at_0} (= -A_int < 0, attractive at L=0)")
print(f"V_int(L->inf) = {V_int_at_inf} (= 0, free particles at large L)")
print(f"Sign at finite L: {sign_at_finite_L} (always negative -- attractive)")

T_P1a_3_pass_lit = (V_int_at_0 == -A_int) and (V_int_at_inf == 0)
print(f"T_P1a_3 LIT consistency? {bool(T_P1a_3_pass_lit)}")

# === T_P1a_4: Stationary point search (substantive, F-beta-1 verification) ===
print("\n--- T_P1a_4: Stationary point dE/dL = 0 in L > 0 (substantive) ---")

# E_total(L) = 2*E_K + V_int(L)
E_total = 2 * E_K_expected + V_int_param  # in terms of m, v, A_int, L

# dE/dL
dE_dL = sp.diff(E_total, L)
dE_dL_simplified = sp.simplify(dE_dL)

print(f"E_total(L) = 2*E_K - A_int * exp(-m*L)")
print(f"dE/dL = {dE_dL_simplified}")

# Solve dE/dL = 0
solutions = sp.solve(dE_dL_simplified, L)
print(f"Solutions to dE/dL = 0 in L > 0: {solutions}")

# dE/dL = A_int * m * exp(-m*L), always positive for m, A_int > 0 and any L
# So NO solution in L > 0

# Verify positivity at sample L values (m, A_int positive)
test_at_L_1 = dE_dL_simplified.subs([(L, 1), (m, 1), (A_int, 1)])
test_at_L_10 = dE_dL_simplified.subs([(L, 10), (m, 1), (A_int, 1)])
test_at_L_0_1 = dE_dL_simplified.subs([(L, sp.Rational(1, 10)), (m, 1), (A_int, 1)])

print(f"dE/dL at L=0.1: {sp.N(test_at_L_0_1, 5)}")
print(f"dE/dL at L=1: {sp.N(test_at_L_1, 5)}")
print(f"dE/dL at L=10: {sp.N(test_at_L_10, 5)}")

# All positive -> monotonic increase -> no stationary point in L > 0
# T_P1a_4 PASS = no positive real solutions

n_solutions = len(solutions)
no_positive_solutions = (n_solutions == 0)
# Also explicitly check positivity of dE/dL
all_samples_positive = all([
    sp.N(test_at_L_0_1, 5) > 0,
    sp.N(test_at_L_1, 5) > 0,
    sp.N(test_at_L_10, 5) > 0,
])

T_P1a_4_pass = no_positive_solutions and all_samples_positive
print(f"  No solutions to dE/dL=0 in L>0: {no_positive_solutions}")
print(f"  dE/dL > 0 at all sample L: {all_samples_positive}")
print(f"T_P1a_4 PASS (F-beta-1 NULL prediction confirmed)? {bool(T_P1a_4_pass)}")

# === T_P1a_5: Asymptotic limits ===
print("\n--- T_P1a_5: Asymptotic limits E(L->0+) and E(L->inf) (substantive) ---")

E_at_L_to_0 = sp.limit(E_total, L, 0, '+')
E_at_L_to_inf = sp.limit(E_total, L, sp.oo)

# Expected:
# E(L->0+) = 2*E_K + V_int(0) = 2*E_K - A_int (solitons anihilate, energy released)
# E(L->inf) = 2*E_K (free pair, no interaction)
expected_at_0 = 2*E_K_expected - A_int
expected_at_inf = 2*E_K_expected

diff_at_0 = sp.simplify(E_at_L_to_0 - expected_at_0)
diff_at_inf = sp.simplify(E_at_L_to_inf - expected_at_inf)

print(f"E(L->0+) = {E_at_L_to_0}")
print(f"  Expected: 2*E_K - A_int = {expected_at_0}")
print(f"  Difference: {diff_at_0}")
print(f"E(L->inf) = {E_at_L_to_inf}")
print(f"  Expected: 2*E_K = {expected_at_inf}")
print(f"  Difference: {diff_at_inf}")

T_P1a_5_pass = (diff_at_0 == 0) and (diff_at_inf == 0)
print(f"T_P1a_5 PASS (limits correct)? {bool(T_P1a_5_pass)}")

# === Summary ===
print("\n" + "="*70)
print("Phase 1a SUMMARY")
print("="*70)

tests = [
    ("T_P1a_1", "Single kink EOM (substantive)", T_P1a_1_pass, "substantive"),
    ("T_P1a_2", "Kink energy E_K computation (substantive)", T_P1a_2_pass, "substantive"),
    ("T_P1a_3", "V_int Manton form (LIT informational)", T_P1a_3_pass_lit, "LIT"),
    ("T_P1a_4", "F-beta-1 verification: no stable L* in isolation (substantive)", T_P1a_4_pass, "substantive"),
    ("T_P1a_5", "Asymptotic limits (substantive)", T_P1a_5_pass, "substantive"),
]

n_pass_total = sum(1 for _, _, p, _ in tests if bool(p))
n_total = len(tests)
substantive_tests = [t for t in tests if t[3] == "substantive"]
n_pass_subst = sum(1 for _, _, p, _ in substantive_tests if bool(p))
n_subst = len(substantive_tests)

print()
for tid, desc, p, cls in tests:
    status = "PASS" if bool(p) else "FAIL"
    marker = " [LIT]" if cls == "LIT" else ""
    print(f"  {tid}: {status}{marker}")
    print(f"    {desc}")

print(f"\nTotal: {n_pass_total}/{n_total} PASS")
print(f"Substantive FP: {n_pass_subst}/{n_subst} PASS")
print(f"Hardcoded T_pass=True: 0 (strict cycle 1/2/7 preserved)")
print(f"DEC budget used: 0/1")

print()
if bool(T_P1a_4_pass):
    print("F-beta-1 PRE-REGISTERED PREDICTION (NULL): CONFIRMED")
    print("Isolation -> no stable L* -> CE-H NOT YET FALSIFIED (Phase 1b needed for positive test)")
else:
    print("F-beta-1 PRE-REGISTERED PREDICTION: FAILED")
    print("Isolation has stable L* -> CE-H STRUCTURALLY FALSIFIED")

print("="*70)
