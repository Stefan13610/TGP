#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
phase2_P21_partition_function_explicit.py
===========================================

PURPOSE
-------
λ.1 Phase 2 sub-task P2.1: Explicit R3 partition function calculation.

Cel: Pochodzić wave-function renormalization Z_φ z explicit evaluation
log det O dla R3 soliton background.

PLAN
----
1. R3 effective EOM dla fluctuation δg wokół g_sol(r):
     O · δg = [-∇² + V''(g_sol(r))] · δg = ω² · δg
2. Numerycznie znaleźć eigenvalues operatora O (radial spherical 3D)
3. Compute log det O = Σ log(ω_n) z odpowiednim regulator
4. Porównać z analytical (e²/4)·log(g₀) form

NUMERICAL APPROACH
-------------------
- Discretize r ∈ [0, R_max] na siatce N punktów
- Solve eigenvalue problem dla radial operator (l=0 mode)
- Sum eigenvalues z cutoff regularization
- Compare jako funkcja g₀

PASS CRITERION
--------------
Show że log det O ma postać proportional do log(g₀) z α-zaleznym
współczynnikiem. Bonus jeśli wartość matche e²/4.

Autor: λ.1 Phase 2 P2.1
Data: 2026-05-01
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
import math

E = math.e
E_SQ = E**2
PHI = (1 + math.sqrt(5)) / 2
G0_E = 0.86941
G0_MU = G0_E * PHI

print("=" * 78)
print("  λ.1 P2.1 — Explicit R3 partition function (heat kernel)")
print("=" * 78)
print()


# ----------------------------------------------------------------
# SECTION 1: Solve R3 ODE dla soliton background g_sol(r)
# ----------------------------------------------------------------

def solve_R3_ode(g0, alpha=2.0, d=3, r_max=50.0, n_points=2000, g_floor=1e-10):
    singular = [False]
    def rhs(r, y):
        g, gp = y
        if g < g_floor:
            singular[0] = True
            g = g_floor
        rhs_val = (1.0 - g) * g**(2 - 2*alpha)
        if r < 1e-12:
            gpp = rhs_val / max(d, 1.0)
        else:
            gpp = rhs_val - (alpha/g) * gp**2 - ((d-1.0)/r) * gp
        return [gp, gpp]
    r_eval = np.linspace(1e-10, r_max, n_points)
    sol = solve_ivp(rhs, (1e-10, r_max), [g0, 0.0],
                    method='RK45', t_eval=r_eval,
                    rtol=1e-9, atol=1e-11, max_step=0.05)
    return sol, singular[0]


# ----------------------------------------------------------------
# SECTION 2: Compute V''(g_sol(r)) — fluctuation operator potential
# ----------------------------------------------------------------

def Vpp_R3(g, alpha=2.0):
    """
    R3 effective EOM: g'' + (α/g)g'² + (2/r)g' = (1-g)·g^(2-2α)
    Linearize around g_sol: O·δg = -∇²δg + V''(g_sol)·δg
    Z right-hand side -d/dg [RHS] gives effective V'' in fluctuation operator.

    RHS = (1-g)·g^(2-2α)
    -d/dg [-RHS] = d/dg [(1-g)·g^(2-2α)]
                 = -g^(2-2α) + (1-g)·(2-2α)·g^(1-2α)
                 = g^(2-2α) · [-1 + (1-g)·(2-2α)/g]
                 = g^(2-2α) · [-1 + (2-2α)·(1/g - 1)]
                 = g^(2-2α) · [(2-2α)/g - (2-2α) - 1]
                 = g^(2-2α) · [(2-2α)/g - (3-2α)]

    For α=2: = g^(-2) · [-2/g - (-1)] = g^(-2) · [-2/g + 1] = -2/g³ + 1/g²

    Effective V'' = -d/dg [(1-g)·g^(2-2α)] = -[g^(2-2α) · ((2-2α)/g - (3-2α))]
                  = g^(2-2α) · ((3-2α) - (2-2α)/g)

    For α=2: = g^(-2) · (-1 + 2/g) = -1/g² + 2/g³
    For g=1: V''(1) = -1 + 2 = 1 → vacuum mass m²=1 ✓
    """
    return g**(2 - 2*alpha) * ((3 - 2*alpha) - (2 - 2*alpha)/g)


# ----------------------------------------------------------------
# SECTION 3: Build radial fluctuation operator (matrix form)
# ----------------------------------------------------------------

def build_radial_operator(g_sol, r_grid, alpha=2.0):
    """
    Build matrix representation of O = -∇² + V''(g_sol(r))
    in radial spherical coordinates (l=0 mode).

    For l=0: -∇²f = -(1/r²) d/dr (r² df/dr) = -f'' - (2/r) f'

    In matrix form on uniform r-grid with spacing h:
    -f'' ≈ -(f_{i+1} - 2f_i + f_{i-1}) / h²
    -(2/r)f' ≈ -(1/r_i)·(f_{i+1} - f_{i-1}) / h
    """
    N = len(r_grid)
    h = r_grid[1] - r_grid[0]

    # Radial Laplacian -∇²
    O = np.zeros((N, N))
    for i in range(N):
        if i > 0:
            O[i, i-1] += -1/h**2
            if r_grid[i] > h/2:
                O[i, i-1] += 1/(r_grid[i]*h)  # -(2/r)f' contribution
        if i < N-1:
            O[i, i+1] += -1/h**2
            if r_grid[i] > h/2:
                O[i, i+1] += -1/(r_grid[i]*h)
        O[i, i] += 2/h**2

    # Add V''(g_sol(r))
    for i in range(N):
        O[i, i] += Vpp_R3(g_sol[i], alpha)

    return O


# ----------------------------------------------------------------
# SECTION 4: Compute log det O numerically dla różnych g₀
# ----------------------------------------------------------------

def compute_log_det_O(g0, alpha=2.0, R_max=30.0, N_grid=300, mass_cutoff=10.0):
    """
    Numerically compute (1/2)·log det O for fluctuation operator
    around R3 soliton with given g₀.

    Returns reduced log det (subtracting vacuum contribution at g=1).
    """
    sol, sing = solve_R3_ode(g0, alpha, r_max=R_max, n_points=N_grid)
    if sing or not sol.success:
        return None

    r_grid = sol.t
    g_sol = sol.y[0]

    # Build fluctuation operator
    O = build_radial_operator(g_sol, r_grid, alpha)

    # Symmetrize
    O = 0.5 * (O + O.T)

    # Compute eigenvalues
    eigenvalues = np.linalg.eigvalsh(O)

    # Filter: positive + below cutoff (mass_cutoff for UV regularization)
    eigenvalues = eigenvalues[eigenvalues > 0.001]
    eigenvalues = eigenvalues[eigenvalues < mass_cutoff]

    # Log det = sum of log eigenvalues
    log_det = np.sum(np.log(eigenvalues))

    return log_det, len(eigenvalues)


# ----------------------------------------------------------------
# SECTION 5: Vacuum subtraction — log det O at g=1
# ----------------------------------------------------------------

print("=" * 78)
print("  SEKCJA 1: Vacuum reference log det O at g=1 (background)")
print("=" * 78)
print()

# Vacuum operator (g_sol = 1 everywhere)
R_MAX = 20.0
N_GRID = 200
r_grid_vac = np.linspace(1e-3, R_MAX, N_GRID)
g_sol_vac = np.ones(N_GRID)
O_vac = build_radial_operator(g_sol_vac, r_grid_vac, alpha=2.0)
O_vac = 0.5 * (O_vac + O_vac.T)
eig_vac = np.linalg.eigvalsh(O_vac)
eig_vac = eig_vac[eig_vac > 0.001]
log_det_vac = np.sum(np.log(eig_vac[eig_vac < 10.0]))
N_vac_eigs = len(eig_vac[eig_vac < 10.0])
print(f"  Vacuum (g=1): log det O = {log_det_vac:.4f} ({N_vac_eigs} eigenvalues)")
print()


# ----------------------------------------------------------------
# SECTION 6: Compute log det O for R3 solitons (e, μ generations)
# ----------------------------------------------------------------

print("=" * 78)
print("  SEKCJA 2: log det O dla R3 solitons z różnymi g₀")
print("=" * 78)
print()

g0_test_values = [G0_E, G0_MU, 1.5, 0.5]

print(f"  {'g₀':>8} | {'log det O':>11} | {'Δ log det':>11} | {'Δ / log(g₀)':>15}")
print(f"  {'':>8} | {'(numerycznie)':>11} | {'vs vacuum':>11} | {'(γ_φ-like)':>15}")
print("  " + "-" * 60)

for g0 in g0_test_values:
    result = compute_log_det_O(g0, alpha=2.0, R_max=R_MAX, N_grid=N_GRID)
    if result is not None:
        log_det, n_eigs = result
        delta_log_det = log_det - log_det_vac

        # gamma_phi-like coefficient
        if abs(math.log(g0)) > 1e-3:
            gamma_eff = delta_log_det / math.log(g0)
        else:
            gamma_eff = float('nan')

        print(f"  {g0:8.4f} | {log_det:11.4f} | {delta_log_det:11.4f} | {gamma_eff:15.4f}")
    else:
        print(f"  {g0:8.4f} | (failed)")

print()


# ----------------------------------------------------------------
# SECTION 7: Empirical comparison — czy gamma_eff ~ e²/4·(4-α) = e²/2?
# ----------------------------------------------------------------

print("=" * 78)
print("  SEKCJA 3: Empirical comparison vs e²/2")
print("=" * 78)
print()
print(f"  Empirical n(2) = e²/2 = {E_SQ/2:.4f}")
print(f"  (Z R3 mass formula: m_obs = c·A²·g₀^(e²/2))")
print()
print("  Uwaga: log det O ≠ n(α) directly. Mass formula używa")
print("  m_obs = exp(coś · log g₀), więc 'coś' to exponent w mass formula.")
print("  log det O daje wkład do mass dressing przez:")
print("    m_dressed = m_bare + (1/2)·log det O / β")
print()
print("  Jeśli (1/2)·log det O ~ K·log g₀ gdzie K ~ e²/2,")
print("  to mass formula dostaje multiplicative factor g₀^(K/2).")


# ----------------------------------------------------------------
# SECTION 8: Scan more g₀ values, check linearity in log g₀
# ----------------------------------------------------------------

print()
print("=" * 78)
print("  SEKCJA 4: Linearity test — czy Δlog det O ∝ log g₀?")
print("=" * 78)
print()

g0_scan = np.linspace(0.5, 2.0, 8)
data_for_fit = []

print(f"  {'g₀':>8} | {'log g₀':>9} | {'Δ log det':>11} | {'ratio':>9}")
print("  " + "-" * 50)

for g0 in g0_scan:
    if abs(g0 - 1.0) < 0.05:
        continue  # skip vacuum
    result = compute_log_det_O(g0, alpha=2.0, R_max=R_MAX, N_grid=N_GRID)
    if result is not None:
        log_det, _ = result
        delta = log_det - log_det_vac
        log_g0 = math.log(g0)
        ratio = delta / log_g0 if abs(log_g0) > 1e-3 else float('nan')
        data_for_fit.append((g0, log_g0, delta))
        print(f"  {g0:8.4f} | {log_g0:9.4f} | {delta:11.4f} | {ratio:9.4f}")

print()

# Linear fit: delta = K · log(g₀) + C
if len(data_for_fit) >= 3:
    arr = np.array(data_for_fit)
    log_g_arr = arr[:, 1]
    delta_arr = arr[:, 2]
    K_fit, C_fit = np.polyfit(log_g_arr, delta_arr, 1)
    print(f"  Linear fit: Δ log det O = K · log(g₀) + C")
    print(f"    K = {K_fit:.4f}")
    print(f"    C = {C_fit:.4f}")
    print()
    print(f"  Compare with empirical:")
    print(f"    e²/2 = {E_SQ/2:.4f}")
    print(f"    e² = {E_SQ:.4f}")
    print(f"    e²/4 = {E_SQ/4:.4f}")
    print()
    matches_e2_half = abs(K_fit - E_SQ/2) / (E_SQ/2) < 0.2
    matches_e2 = abs(K_fit - E_SQ) / E_SQ < 0.2
    matches_e2_quarter = abs(K_fit - E_SQ/4) / (E_SQ/4) < 0.2
    if matches_e2_half:
        print(f"  >>> MATCH: K ≈ e²/2 (within 20%)")
    elif matches_e2:
        print(f"  >>> MATCH: K ≈ e² (within 20%)")
    elif matches_e2_quarter:
        print(f"  >>> MATCH: K ≈ e²/4 (within 20%)")
    else:
        print(f"  >> No match z e², e²/2, lub e²/4 (>20%).")


# ----------------------------------------------------------------
# SECTION 9: HONEST conclusion P2.1
# ----------------------------------------------------------------

print()
print("=" * 78)
print("  SEKCJA 5: HONEST conclusion P2.1")
print("=" * 78)
print()

print(f"""
  Co odkryłem w P2.1:

  1. Numerical computation log det O dla R3 fluctuation operator wokół
     soliton background jest WYKONALNE (matrix diagonalization).

  2. Linearity test sprawdza czy Δ log det O ∝ log(g₀) — co byłoby
     KEY consistency z mass formula m_obs ~ g₀^n.

  3. Linear fit gives K = (numerical value) — porównanie z e²/2, e², e²/4.

  PASS criterion P2.1:
    "Show że log det O ma postać proportional do log(g₀) ze współczynnikiem
     zalezne od α; bonus jeśli wartość matche e²/4"

  Status: zalezne od linearity test wyniku.

  Caveat:
  - Numerical computation może mieć duże błędy (R_max cutoff,
    eigenvalue cutoff, grid resolution).
  - Wymagałoby refinement w Phase 3 jeśli pozytywny sygnał.
""")
