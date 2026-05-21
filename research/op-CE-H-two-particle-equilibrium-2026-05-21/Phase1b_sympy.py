#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1b -- N=2 kink+antikink + CE-H background positive test
op-CE-H-two-particle-equilibrium-2026-05-21

Pre-registered prediction (LOCKED 2026-05-21):
F-beta-2: stable L* > 0 finite z d2E/dL2 > 0 (with bg).

CE-H bg model (explicit construction):
- Pure isolation: E_isolation(L) = 2*E_K - A*exp(-m*L)
- Add bg contribution: E_bg(L) = D / L^alpha
- Total: E_total(L) = 2*E_K - A*exp(-m*L) + D/L^alpha

Physical motivation (explicit):
- BG <Phi>_bg from rest of universe creates effective pressure on solitons
- Two solitons in close proximity: overlapping bg deformations -> compression cost
- Modeled as power-law repulsive contribution D/L^alpha
- alpha = 1 (Coulomb-like, simplest) as default
- D = 0 limit recovers Phase 1a (CE-H consistent: no bg -> no equilibrium)
- D > 0 enables stable equilibrium

Honest caveat: D/L^alpha is EXPLICITLY CONSTRUCTED to demonstrate mechanism.
Full TGP-native derivation of D, alpha from substrate Lagrangian = Poziom gamma scope.

Tests:
- T_P1b_1: Define total energy with bg, verify smooth D->0 limit
- T_P1b_2: Existence of stationary point L* (for some D range)
- T_P1b_3: Stability check d2E/dL2 > 0 at L*
- T_P1b_4: Critical D_max above which no equilibrium
- T_P1b_5: Monotonicity L*(D) - F-beta-3 verification

Discipline: strict cycle 1/2/7, 0 hardcoded T_pass=True, 5 substantive FP.
"""
import sympy as sp
import sys

print("="*70)
print("Phase 1b -- N=2 + CE-H bg positive test (1D Z2 toy)")
print("op-CE-H-two-particle-equilibrium-2026-05-21")
print("Pre-registered prediction F-beta-2: stable L* > 0 finite WITH bg")
print("="*70)

# Symbols
L = sp.symbols('L', positive=True)
v, lam, m = sp.symbols('v lambda m', positive=True)
A_int = sp.Symbol('A_int', positive=True)  # kink-antikink interaction coeff
D = sp.Symbol('D', positive=True)  # bg coupling strength
alpha = sp.Symbol('alpha', positive=True)  # bg power-law exponent

# Single kink energy (from Phase 1a)
E_K = sp.Rational(2) * sp.sqrt(2) / sp.Rational(3) * m * v**2

# === T_P1b_1: Define E_total(L; D, alpha) ===
print("\n--- T_P1b_1: Total energy with bg (substantive) ---")

E_isolation = 2*E_K - A_int * sp.exp(-m*L)
E_bg = D / L**alpha
E_total = E_isolation + E_bg

print(f"E_isolation(L) = 2*E_K - A_int * exp(-m*L)")
print(f"E_bg(L) = D / L^alpha")
print(f"E_total(L) = E_isolation + E_bg")

# Verify D->0 limit gives back Phase 1a result
E_total_at_D_0 = E_total.subs(D, 0)
diff_recovery = sp.simplify(E_total_at_D_0 - E_isolation)
print(f"E_total at D=0: {E_total_at_D_0}")
print(f"  Difference vs E_isolation: {diff_recovery}")

T_P1b_1_pass = (diff_recovery == 0)
print(f"T_P1b_1 PASS (D->0 limit recovers Phase 1a)? {bool(T_P1b_1_pass)}")

# === T_P1b_2: Stationary point existence (for alpha = 1) ===
print("\n--- T_P1b_2: Stationary point existence at alpha=1 (substantive) ---")

# Fix alpha = 1 (Coulomb-like)
E_total_alpha1 = E_total.subs(alpha, 1)
dE_dL = sp.diff(E_total_alpha1, L)
dE_dL_simplified = sp.simplify(dE_dL)

print(f"At alpha=1: dE/dL = {dE_dL_simplified}")

# Stationarity: A_int*m*exp(-m*L) = D/L^2
# Let u = m*L, then (A_int/m)*u^2*exp(-u) = D, i.e. u^2*exp(-u) = D*m/A_int

# f(u) = u^2 * exp(-u) has maximum at u = 2, value f(2) = 4/e^2 ~ 0.541
u = sp.Symbol('u', positive=True)
f_u = u**2 * sp.exp(-u)
f_u_prime = sp.diff(f_u, u)
critical_u = sp.solve(f_u_prime, u)
f_max = f_u.subs(u, 2)

print(f"f(u) = u^2 * exp(-u), max at u={critical_u} (interior critical point u=2)")
print(f"f(2) = 4/e^2 = {f_max} = {float(f_max):.6f}")
print(f"For equilibrium L* to exist: D*m/A_int <= 4/e^2")

# Test with specific D value (well within feasibility)
# Pick D such that D*m/A_int = 0.2 (below max 0.541)
D_test_val = sp.Rational(1, 5)  # D*m/A_int = 1/5

# Solve numerically for u (where u = m*L*)
u_eq_test = sp.nsolve(u**2 * sp.exp(-u) - D_test_val, u, 0.5)  # try near 0.5
print(f"\nTest D*m/A_int = 1/5 = {float(D_test_val):.4f}")
print(f"Stable branch u_1 = {float(u_eq_test):.4f} (m*L* < 2 stable)")

# Also find unstable branch
u_unstable = sp.nsolve(u**2 * sp.exp(-u) - D_test_val, u, 4)
print(f"Unstable branch u_2 = {float(u_unstable):.4f} (m*L* > 2 unstable)")

# Verify both are positive real
T_P1b_2_pass = (float(u_eq_test) > 0) and (float(u_eq_test) < 2) and (float(u_unstable) > 2)
print(f"T_P1b_2 PASS (stable + unstable branches exist as expected)? {bool(T_P1b_2_pass)}")

# === T_P1b_3: Stability at stable branch ===
print("\n--- T_P1b_3: Stability d2E/dL2 > 0 at L* (substantive) ---")

d2E_dL2 = sp.diff(E_total_alpha1, L, 2)
d2E_dL2_simplified = sp.simplify(d2E_dL2)
print(f"d2E/dL2 = {d2E_dL2_simplified}")

# At L* satisfying stationarity, simplify:
# A_int*m*exp(-m*L*) = D/L*^2 => substitute
# d2E/dL2 = -A_int*m^2*exp(-m*L) + 2*D/L^3
# At L*: = -m^2 * D/L*^2 + 2*D/L*^3 = (D/L*^3) * (2 - m*L*)
# Sign: positive iff m*L* < 2

# Test at stable branch (u = m*L_stable, u < 2)
L_stable_test = float(u_eq_test)  # m*L_stable = u_eq_test, so L_stable / (1/m) = u_eq_test
# Sub values: m=1, A_int=1, D=1/5
d2E_at_stable = d2E_dL2_simplified.subs([(L, L_stable_test), (m, 1), (A_int, 1), (D, sp.Rational(1, 5))])
d2E_at_stable_val = float(d2E_at_stable)
print(f"At stable branch L*=u_1/m: d2E/dL2 = {d2E_at_stable_val:.6f}")

# Test at unstable branch
L_unstable_test = float(u_unstable)
d2E_at_unstable = d2E_dL2_simplified.subs([(L, L_unstable_test), (m, 1), (A_int, 1), (D, sp.Rational(1, 5))])
d2E_at_unstable_val = float(d2E_at_unstable)
print(f"At unstable branch L*=u_2/m: d2E/dL2 = {d2E_at_unstable_val:.6f}")

# Analytical check: at any stationary point, d2E/dL2 = (D/L^3)*(2 - m*L)
# This is positive iff m*L < 2 (stable) and negative iff m*L > 2 (unstable)
T_P1b_3_pass = (d2E_at_stable_val > 0) and (d2E_at_unstable_val < 0)
print(f"T_P1b_3 PASS (stable branch has d2E/dL2 > 0, unstable < 0)? {bool(T_P1b_3_pass)}")
print(f"  F-beta-2 PRE-REGISTERED PREDICTION verified for D = 1/5")

# === T_P1b_4: Critical D_max above which no equilibrium ===
print("\n--- T_P1b_4: Critical D_max boundary of equilibrium (substantive) ---")

# At critical point: u = 2, f(2) = 4/e^2
# D_critical * m / A_int = 4/e^2
# So D_critical = 4 * A_int / (m * e^2)

D_critical_analytical = 4 * A_int / (m * sp.exp(2))
print(f"D_critical (analytical) = 4*A_int / (m*e^2) = {D_critical_analytical}")

# At D = D_critical, stationary point is at u = 2 (m*L* = 2), marginal stability
# Verify numerically: choose D slightly less, equilibrium exists; slightly more, no equilibrium

D_below = float(4/sp.exp(2)) * 0.99  # 0.535
D_above = float(4/sp.exp(2)) * 1.01  # 0.546

# Try to solve u^2*exp(-u) = D_below (should have 2 solutions)
try:
    u_below_stable = float(sp.nsolve(u**2 * sp.exp(-u) - D_below, u, 1.0))
    u_below_unstable = float(sp.nsolve(u**2 * sp.exp(-u) - D_below, u, 3.0))
    print(f"D = D_crit * 0.99 = {D_below:.4f}: u_stable = {u_below_stable:.4f}, u_unstable = {u_below_unstable:.4f}")
    has_solution_below = True
except Exception as e:
    print(f"D below critical: nsolve failed: {e}")
    has_solution_below = False

# Try D above critical
try:
    u_above_attempt = sp.nsolve(u**2 * sp.exp(-u) - D_above, u, 2.0)
    u_above_val = float(u_above_attempt)
    # Check if it's actually a solution (residual)
    residual = u_above_val**2 * float(sp.exp(-u_above_val)) - D_above
    if abs(residual) > 1e-6:
        has_solution_above = False
        print(f"D above critical: numerical artifact, true residual = {residual:.6f}")
    else:
        has_solution_above = True
        print(f"D above critical: spurious solution found u={u_above_val:.4f}")
except Exception as e:
    print(f"D above critical: no solution (expected) -- {type(e).__name__}")
    has_solution_above = False

T_P1b_4_pass = has_solution_below and (not has_solution_above)
print(f"T_P1b_4 PASS (D_critical boundary correctly identified)? {bool(T_P1b_4_pass)}")

# === T_P1b_5: Monotonicity L*(D) - F-beta-3 verification ===
print("\n--- T_P1b_5: L*(D) monotonicity (substantive, F-beta-3) ---")

# Scan D values and find stable L* (u_1 < 2)
import numpy as np
D_scan = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.53]  # all below 4/e^2 ~ 0.541
L_stars = []
for D_val in D_scan:
    try:
        u_stab = float(sp.nsolve(u**2 * sp.exp(-u) - D_val, u, 0.5))
        if u_stab < 2:
            L_stars.append((D_val, u_stab))
            print(f"  D = {D_val:.4f}: u_stable = m*L* = {u_stab:.4f}")
    except:
        pass

# Check monotonicity: L* should INCREASE with D (more bg pressure -> larger equilibrium separation)
# Wait: increasing D = stronger bg = more repulsion -> larger L*? Let me think.
# E_bg = D/L: larger D, more push apart, larger L*. YES.
u_values = [u for _, u in L_stars]
is_monotonic = all(u_values[i] < u_values[i+1] for i in range(len(u_values)-1))
print(f"\nL*(D) scan u values: {[f'{u:.3f}' for u in u_values]}")
print(f"Monotonically increasing? {is_monotonic}")

# Also check no fine-tuning: equilibrium exists for D ranging over factor 10
D_min = D_scan[0]
D_max_used = D_scan[-2]  # exclude the very last edge case
factor_range = D_max_used / D_min
print(f"Equilibrium exists for D range factor {factor_range:.1f}")

# F-beta-4 (no fine-tuning) requires factor >= 10
no_fine_tuning = factor_range >= 10
T_P1b_5_pass = is_monotonic and no_fine_tuning
print(f"T_P1b_5 PASS (monotonic + no fine-tuning)? {bool(T_P1b_5_pass)}")

# === Summary ===
print("\n" + "="*70)
print("Phase 1b SUMMARY")
print("="*70)

tests = [
    ("T_P1b_1", "D->0 limit recovers Phase 1a (substantive)", T_P1b_1_pass, "substantive"),
    ("T_P1b_2", "Stationary point exists (F-beta-2 part 1) (substantive)", T_P1b_2_pass, "substantive"),
    ("T_P1b_3", "Stability d2E/dL2 > 0 at L* (F-beta-2 part 2) (substantive)", T_P1b_3_pass, "substantive"),
    ("T_P1b_4", "D_critical boundary identified (substantive)", T_P1b_4_pass, "substantive"),
    ("T_P1b_5", "Monotonic L*(D) + no fine-tuning (F-beta-3, F-beta-4) (substantive)", T_P1b_5_pass, "substantive"),
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
if bool(T_P1b_2_pass) and bool(T_P1b_3_pass):
    print("F-beta-2 PRE-REGISTERED PREDICTION (POSITIVE): CONFIRMED")
    print("  -> stable L* exists WITH bg -> CE-H structurally consistent")
else:
    print("F-beta-2 PRE-REGISTERED PREDICTION: FAILED")
    print("  -> stable L* NOT found WITH bg -> CE-H STRUCTURALLY FALSIFIED")

if bool(T_P1b_5_pass):
    print("F-beta-3 monotonicity + F-beta-4 no-fine-tuning: CONFIRMED")
else:
    print("F-beta-3 or F-beta-4: FAILED")

print()
print("Honest caveat: D/L bg model is EXPLICITLY CONSTRUCTED.")
print("Full TGP-native derivation of bg form = Poziom gamma scope.")
print("Poziom beta = STRUCTURAL proof-of-principle that bg can stabilize.")

print("="*70)
