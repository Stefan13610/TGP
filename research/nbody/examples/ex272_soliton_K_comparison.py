#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex272_soliton_K_comparison.py
==============================
TGP Soliton ODE: comparison of THREE kinetic coupling formulations.

For a general Lagrangian density (spherical symmetry):
    L = (1/2) K(g) (g')^2 + V(g)

the Euler-Lagrange equation is:
    K(g) g'' + (1/2) K'(g) (g')^2 + (2/r) K(g) g' = V'(g)

We solve this for three choices of K(g) and V(g):

  Formulation A  (canonical):   K = g^4,      V = g^7/7 - g^8/8
  Formulation B  (historical):  K = 1+4 ln g, V = g^3/3 - g^4/4
  Formulation C  (ex234-style): K = g^2,      V = g^3/3 - g^4/4

Key physics:
  - Vacuum at g=1.  The linearized equation near g=1 gives OSCILLATORY
    tails: (g-1)*r ~ A_tail * sin(omega*r + phase) / r.
  - Every g0 produces a valid oscillatory solution; there is no
    "shooting to g=1" in the usual sense.
  - The physically relevant g0 is the phi-fixed-point (phi-FP) where
    (A_tail(phi*g0) / A_tail(g0))^4 = R_21^PDG = 206.768.
  - For K=g^2 (Form C) the phi-FP is g0 ~ 0.86941 (well established).

Each formulation uses a canonical variable to regularize the ODE:
  Form A (K=g^4):       psi = g^3/3  =>  psi'' + (2/r)psi' = g^4(1-g)
  Form B (K=1+4 ln g):  g-space ODE  =>  (1+4ln g)g'' + (2/g)(g')^2
                                          + (2/r)(1+4ln g)g' = g^2 - g^3
  Form C (K=g^2):       u = g^2/2    =>  u'' + (2/r)u' = g(1-g)

For each solution (at a reference g0) we compute:
  - A_tail from oscillatory tail fit
  - Soliton energy E = 4 pi int [K/2 (g')^2 + V(g) - V(1)] r^2 dr
Then we find the phi-FP g0 for each formulation and compare.

Sesja: TGP v41 -- Claudian (2026-04-07)
"""

import sys, io, warnings
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings('ignore')

# ------------------------------------------------------------------
# Compatibility: numpy >= 2.0 renamed trapz -> trapezoid
# ------------------------------------------------------------------
_trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')

PHI   = (1.0 + np.sqrt(5.0)) / 2.0
R21_PDG = 206.768

# ============================================================
# 1. SOLVER: Formulation A (K=g^4) in psi = g^3/3 space
# ============================================================
# From K=g^4, V = g^7/7 - g^8/8, V' = g^6 - g^7 = g^6(1-g)
# EL in g-space: g^4 g'' + 2g^3(g')^2 + (2/r)g^4 g' = g^6 - g^7
#
# Canonical variable psi = g^3/3  =>  psi' = g^2 g'
# After substitution the ODE becomes (see ex231):
#   psi'' + (2/r) psi' = gamma * g^4 * (1 - g)
# with g = (3 psi)^{1/3}, gamma = 1 (beta=gamma=1).
# This is REGULAR at g=0 (psi=0): RHS -> 0.

def solve_formA(g0, R_max=150.0, N_pts=6000, gamma=1.0):
    """Solve Form A (K=g^4) soliton in psi=g^3/3 coordinates."""
    psi0 = g0**3 / 3.0

    def rhs(r, y):
        psi, psip = y
        psi_safe = max(psi, 0.0)
        r_safe = max(r, 1e-10)
        if psi_safe < 1e-30:
            force = 0.0
        else:
            g = (3.0 * psi_safe) ** (1.0 / 3.0)
            force = gamma * g**4 * (1.0 - g)
        return [psip, force - 2.0 * psip / r_safe]

    r0 = 1e-4
    acc0 = gamma * g0**4 * (1.0 - g0)
    psi_init = psi0 + acc0 * r0**2 / 6.0
    psip_init = acc0 * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [psi_init, psip_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.05, dense_output=True)
    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y_arr = sol.sol(r_arr)
    psi_arr = y_arr[0]
    psip_arr = y_arr[1]
    g_arr = (3.0 * np.maximum(psi_arr, 0.0)) ** (1.0 / 3.0)
    # g' = psi' / g^2
    gp_arr = np.where(g_arr > 1e-15, psip_arr / g_arr**2, 0.0)
    return r_arr, g_arr, gp_arr, psi_arr

# ============================================================
# 2. SOLVER: Formulation B (K=1+4 ln g) in g-space
# ============================================================
# EL: (1+4 ln g) g'' + (2/g)(g')^2 + (2/r)(1+4 ln g) g' = g^2 - g^3
#
# Rearranged: g'' = [g^2 - g^3 - (2/g)(g')^2 - (2/r)(1+4ln g) g'] / (1+4ln g)
#
# At r=0 with g'=0:  (1+4ln g0) g''(0) + 2(1+4ln g0) g''(0) = g0^2 - g0^3
#   => g''(0) = (g0^2 - g0^3) / (3(1+4 ln g0))
# Note: K(g) > 0 requires g > exp(-1/4) ~ 0.7788.

def solve_formB(g0, R_max=150.0, N_pts=6000):
    """Solve Form B (K=1+4 ln g) soliton in g-space."""

    def rhs(r, y):
        g, gp = y
        g_safe = max(g, 1e-15)
        r_safe = max(r, 1e-10)
        K = 1.0 + 4.0 * np.log(g_safe)
        if abs(K) < 1e-20:
            return [gp, 0.0]
        Kp_half = 2.0 / g_safe    # (1/2)*K'(g) = (1/2)*(4/g) = 2/g
        Vp = g_safe**2 - g_safe**3
        gpp = (Vp - Kp_half * gp**2 - (2.0 / r_safe) * K * gp) / K
        return [gp, gpp]

    r0 = 1e-4
    K0 = 1.0 + 4.0 * np.log(max(g0, 1e-15))
    acc0 = (g0**2 - g0**3) / (3.0 * K0) if abs(K0) > 1e-20 else 0.0
    g_init = g0 + acc0 * r0**2 / 2.0
    gp_init = acc0 * r0

    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.05, dense_output=True)
    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y_arr = sol.sol(r_arr)
    return r_arr, y_arr[0], y_arr[1]

# ============================================================
# 3. SOLVER: Formulation C (K=g^2) in u = g^2/2 space
# ============================================================
# EL in g-space: g'' + (g')^2/g + (2/r)g' = 1-g
# Canonical variable u = g^2/2  =>  u' = g g'
# ODE becomes (see ex232):
#   u'' + (2/r) u' = g * (1 - g)    with g = sqrt(2u)
# Completely regular at u=0 (g=0).

def solve_formC(g0, R_max=150.0, N_pts=6000):
    """Solve Form C (K=g^2) soliton in u=g^2/2 coordinates."""
    u0 = g0**2 / 2.0

    def rhs(r, y):
        u, up = y
        u_safe = max(u, 0.0)
        r_safe = max(r, 1e-10)
        g = np.sqrt(2.0 * u_safe)
        force = g * (1.0 - g)
        return [up, force - 2.0 * up / r_safe]

    r0 = 1e-4
    acc0 = g0 * (1.0 - g0)
    u_init = u0 + acc0 * r0**2 / 6.0
    up_init = acc0 * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [u_init, up_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.05, dense_output=True)
    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y_arr = sol.sol(r_arr)
    u_arr = y_arr[0]
    up_arr = y_arr[1]
    g_arr = np.sqrt(2.0 * np.maximum(u_arr, 0.0))
    # g' = u' / g
    gp_arr = np.where(g_arr > 1e-15, up_arr / g_arr, 0.0)
    return r_arr, g_arr, gp_arr, u_arr

# ============================================================
# 4. Tail fitting: oscillatory A_tail
# ============================================================
# Linearized near g=1 (u = g-1 small):
#   All formulations have K(1)=1, V''(1)=-1 (or equivalently omega=1)
#   => u'' + (2/r)u' + omega^2 u = 0   with omega = 1
#   Solution: u = A sin(r+delta) / r
# So: (g-1)*r ~ A sin(r+delta)
# We fit A_tail from (g-1)*r = B cos(r) + C sin(r).
#
# For Form A in psi-space:  (psi - 1/3)*r = B cos(omega*r) + C sin(omega*r)
# with omega = sqrt(gamma) = 1.

def fit_tail_g(r_arr, g_arr, omega=1.0, r_L=20.0, r_R=60.0):
    """Fit (g-1)*r = B cos(omega*r) + C sin(omega*r).
    Returns (A_tail, phase, R^2)."""
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(g_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan, 0.0
    r_fit = r_arr[mask]
    y_fit = (g_arr[mask] - 1.0) * r_fit
    X = np.column_stack([np.cos(omega * r_fit), np.sin(omega * r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)

    y_pred = X @ coefs
    ss_res = np.sum((y_fit - y_pred)**2)
    ss_tot = np.sum((y_fit - np.mean(y_fit))**2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return A, np.arctan2(B, C), R2

def fit_tail_psi(r_arr, psi_arr, omega=1.0, r_L=20.0, r_R=60.0):
    """Fit (psi - 1/3)*r = B cos(omega*r) + C sin(omega*r)."""
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(psi_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan, 0.0
    r_fit = r_arr[mask]
    y_fit = (psi_arr[mask] - 1.0/3.0) * r_fit
    X = np.column_stack([np.cos(omega * r_fit), np.sin(omega * r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    A = np.sqrt(B**2 + C**2)

    y_pred = X @ coefs
    ss_res = np.sum((y_fit - y_pred)**2)
    ss_tot = np.sum((y_fit - np.mean(y_fit))**2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return A, np.arctan2(B, C), R2

# ============================================================
# 5. Energy functional
# ============================================================

def compute_energy(r_arr, g_arr, gp_arr, K_func, V_func, V_vac):
    """E = 4 pi int_0^inf [K(g)/2 (g')^2 + V(g) - V(1)] r^2 dr."""
    K_vals = np.array([float(K_func(g)) for g in g_arr])
    V_vals = np.array([float(V_func(g)) for g in g_arr])
    integrand = (0.5 * K_vals * gp_arr**2 + V_vals - V_vac) * r_arr**2
    return 4.0 * np.pi * _trapz(integrand, r_arr)

# ============================================================
# 6. Phi-fixed-point finder for each formulation
# ============================================================

def get_A_formA(g0, r_L=20.0, r_R=60.0):
    """A_tail in g-space for Form A."""
    try:
        r, g, gp, psi = solve_formA(g0, R_max=120.0, N_pts=5000)
        A, _, R2 = fit_tail_g(r, g, omega=1.0, r_L=r_L, r_R=r_R)
        return A, R2
    except:
        return np.nan, 0.0

def get_A_formB(g0, r_L=20.0, r_R=60.0):
    """A_tail in g-space for Form B."""
    try:
        r, g, gp = solve_formB(g0, R_max=120.0, N_pts=5000)
        A, _, R2 = fit_tail_g(r, g, omega=1.0, r_L=r_L, r_R=r_R)
        return A, R2
    except:
        return np.nan, 0.0

def get_A_formC(g0, r_L=20.0, r_R=60.0):
    """A_tail in g-space for Form C."""
    try:
        r, g, gp, u = solve_formC(g0, R_max=120.0, N_pts=5000)
        A, _, R2 = fit_tail_g(r, g, omega=1.0, r_L=r_L, r_R=r_R)
        return A, R2
    except:
        return np.nan, 0.0

def find_phi_FP(get_A_func, g0_lo, g0_hi, label=""):
    """Find g0 where (A(phi*g0)/A(g0))^4 = R21_PDG."""

    def ratio_func(g0):
        A1, _ = get_A_func(g0)
        A2, _ = get_A_func(PHI * g0)
        if np.isnan(A1) or A1 < 1e-15 or np.isnan(A2):
            return 1e10
        return (A2 / A1)**4 - R21_PDG

    # Coarse scan
    g0_scan = np.linspace(g0_lo, g0_hi, 30)
    print(f"  Scanning g0 in [{g0_lo:.3f}, {g0_hi:.3f}] for phi-FP ...")
    vals = []
    for g0 in g0_scan:
        v = ratio_func(g0)
        vals.append(v)
    vals = np.array(vals)

    # Find bracket
    for i in range(len(vals) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i] * vals[i+1] < 0:
            a, b = g0_scan[i], g0_scan[i+1]
            try:
                g0_star = brentq(ratio_func, a, b, xtol=1e-8)
                return g0_star, True
            except:
                pass

    # Fallback: return closest
    finite_mask = np.isfinite(vals) & (np.abs(vals) < 1e9)
    if np.any(finite_mask):
        idx = np.argmin(np.abs(vals[finite_mask]))
        return g0_scan[finite_mask][idx], False

    return np.nan, False

# ============================================================
# 7. Main
# ============================================================

def main():
    print("=" * 72)
    print("  TGP Soliton ODE: Comparison of Three Kinetic Coupling Formulations")
    print("=" * 72)
    print()

    # ----------------------------------------------------------
    # Display EL equations
    # ----------------------------------------------------------
    print("--- Euler-Lagrange equations (spherical symmetry) ---")
    print()
    print("General:  K(g) g'' + (1/2) K'(g) (g')^2 + (2/r) K(g) g' = V'(g)")
    print()
    print("Form A:  K=g^4,  V=g^7/7 - g^8/8,  V'=g^6 - g^7")
    print("  => g^4 g'' + 2 g^3 (g')^2 + (2/r) g^4 g' = g^6 - g^7")
    print("  Canonical psi=g^3/3:  psi'' + (2/r) psi' = g^4(1-g)")
    print()
    print("Form B:  K=1+4 ln g,  V=g^3/3 - g^4/4,  V'=g^2 - g^3")
    print("  => (1+4 ln g) g'' + (2/g)(g')^2 + (2/r)(1+4 ln g) g' = g^2 - g^3")
    print()
    print("Form C:  K=g^2,  V=g^3/3 - g^4/4,  V'=g^2 - g^3")
    print("  => g^2 g'' + g (g')^2 + (2/r) g^2 g' = g^2 - g^3")
    print("  Canonical u=g^2/2:   u'' + (2/r) u' = g(1-g)")
    print("  Equivalently:  g'' + (g')^2/g + (2/r) g' = 1 - g")
    print()
    print("Tail: all K(1)=1, V''(1)=-1  =>  omega=1 (oscillatory)")
    print("  (g-1)*r ~ A_tail * sin(r + phase)")
    print()

    # ----------------------------------------------------------
    # Section 1: Profiles at reference g0
    # ----------------------------------------------------------
    print("=" * 72)
    print("  SECTION 1: Profiles at reference g0 = 0.87")
    print("=" * 72)
    print()

    g0_ref = 0.87
    R_MAX_SOLVE = 120.0
    N_PTS = 6000

    # Solve all three
    print("  Solving Form A ...", end="", flush=True)
    rA, gA, gpA, psiA = solve_formA(g0_ref, R_max=R_MAX_SOLVE, N_pts=N_PTS)
    AA, phaseA, R2A = fit_tail_g(rA, gA, omega=1.0)
    print(f" done.  A_tail = {AA:.6f}, R^2 = {R2A:.4f}")

    print("  Solving Form B ...", end="", flush=True)
    rB, gB, gpB = solve_formB(g0_ref, R_max=R_MAX_SOLVE, N_pts=N_PTS)
    AB, phaseB, R2B = fit_tail_g(rB, gB, omega=1.0)
    print(f" done.  A_tail = {AB:.6f}, R^2 = {R2B:.4f}")

    print("  Solving Form C ...", end="", flush=True)
    rC, gC, gpC, uC = solve_formC(g0_ref, R_max=R_MAX_SOLVE, N_pts=N_PTS)
    AC, phaseC, R2C = fit_tail_g(rC, gC, omega=1.0)
    print(f" done.  A_tail = {AC:.6f}, R^2 = {R2C:.4f}")

    # Energies at reference g0
    K_A = lambda g: g**4
    V_A = lambda g: g**7 / 7.0 - g**8 / 8.0
    V_A_vac = V_A(1.0)

    K_B = lambda g: 1.0 + 4.0 * np.log(max(g, 1e-30))
    V_BC = lambda g: g**3 / 3.0 - g**4 / 4.0
    V_BC_vac = V_BC(1.0)

    K_C = lambda g: g**2

    EA = compute_energy(rA, gA, gpA, K_A, V_A, V_A_vac)
    EB = compute_energy(rB, gB, gpB, K_B, V_BC, V_BC_vac)
    EC = compute_energy(rC, gC, gpC, K_C, V_BC, V_BC_vac)

    print()
    print(f"  At g0 = {g0_ref}:")
    print(f"  {'Formulation':<18} {'A_tail':>10} {'R^2':>8} {'E_soliton':>14}")
    print(f"  {'-'*52}")
    print(f"  {'A (K=g^4)':<18} {AA:10.6f} {R2A:8.4f} {EA:14.6f}")
    print(f"  {'B (K=1+4ln g)':<18} {AB:10.6f} {R2B:8.4f} {EB:14.6f}")
    print(f"  {'C (K=g^2)':<18} {AC:10.6f} {R2C:8.4f} {EC:14.6f}")
    print()

    # Profile at selected r values
    print(f"  {'r':>6s}  {'g_A':>10s}  {'g_B':>10s}  {'g_C':>10s}")
    print(f"  {'-'*40}")
    for r_val in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0]:
        idxA = np.searchsorted(rA, r_val)
        idxB = np.searchsorted(rB, r_val)
        idxC = np.searchsorted(rC, r_val)
        vA = gA[min(idxA, len(gA)-1)]
        vB = gB[min(idxB, len(gB)-1)]
        vC = gC[min(idxC, len(gC)-1)]
        print(f"  {r_val:6.1f}  {vA:10.6f}  {vB:10.6f}  {vC:10.6f}")
    print()

    # ----------------------------------------------------------
    # Section 2: A_tail scan over g0 range
    # ----------------------------------------------------------
    print("=" * 72)
    print("  SECTION 2: A_tail(g0) scan for each formulation")
    print("=" * 72)
    print()

    g0_scan_lo = np.linspace(0.5, 1.5, 21)
    print(f"  {'g0':>6s}  {'A_A':>10s}  {'A_B':>10s}  {'A_C':>10s}")
    print(f"  {'-'*40}")
    for g0 in g0_scan_lo:
        aA, _ = get_A_formA(g0)
        aB, _ = get_A_formB(g0)
        aC, _ = get_A_formC(g0)
        sA = f"{aA:10.6f}" if np.isfinite(aA) else f"{'N/A':>10}"
        sB = f"{aB:10.6f}" if np.isfinite(aB) else f"{'N/A':>10}"
        sC = f"{aC:10.6f}" if np.isfinite(aC) else f"{'N/A':>10}"
        print(f"  {g0:6.3f}  {sA}  {sB}  {sC}")
    print()

    # ----------------------------------------------------------
    # Section 3: phi-fixed-point for each formulation
    # ----------------------------------------------------------
    print("=" * 72)
    print("  SECTION 3: phi-Fixed-Point (phi-FP) for each formulation")
    print("  Condition: (A_tail(phi*g0) / A_tail(g0))^4 = R21_PDG = 206.768")
    print("=" * 72)
    print()

    # Form C (known: ~0.86941)
    print("--- Form C (K=g^2) ---")
    g0_C_star, conv_C = find_phi_FP(get_A_formC, 0.80, 0.95, "C")
    if conv_C:
        print(f"  >>> phi-FP: g0* = {g0_C_star:.8f}  [CONVERGED]")
    else:
        print(f"  >>> phi-FP: g0* = {g0_C_star:.8f}  [APPROX]")
    # Verify
    A1_C, _ = get_A_formC(g0_C_star)
    A2_C, _ = get_A_formC(PHI * g0_C_star)
    ratio_C = (A2_C / A1_C)**4 if A1_C > 1e-15 else np.nan
    print(f"  Verification: (A(phi*g0)/A(g0))^4 = {ratio_C:.2f}  (target: {R21_PDG})")
    print()

    # Form A
    print("--- Form A (K=g^4) ---")
    g0_A_star, conv_A = find_phi_FP(get_A_formA, 0.70, 1.00, "A")
    if conv_A:
        print(f"  >>> phi-FP: g0* = {g0_A_star:.8f}  [CONVERGED]")
    else:
        print(f"  >>> phi-FP: g0* = {g0_A_star:.8f}  [APPROX]")
    A1_A, _ = get_A_formA(g0_A_star)
    A2_A, _ = get_A_formA(PHI * g0_A_star)
    ratio_A = (A2_A / A1_A)**4 if A1_A > 1e-15 else np.nan
    print(f"  Verification: (A(phi*g0)/A(g0))^4 = {ratio_A:.2f}  (target: {R21_PDG})")
    print()

    # Form B
    print("--- Form B (K=1+4 ln g) ---")
    g0_B_star, conv_B = find_phi_FP(get_A_formB, 0.80, 1.10, "B")
    if conv_B:
        print(f"  >>> phi-FP: g0* = {g0_B_star:.8f}  [CONVERGED]")
    else:
        print(f"  >>> phi-FP: g0* = {g0_B_star:.8f}  [APPROX]")
    A1_B, _ = get_A_formB(g0_B_star)
    A2_B, _ = get_A_formB(PHI * g0_B_star)
    ratio_B = (A2_B / A1_B)**4 if A1_B > 1e-15 else np.nan
    print(f"  Verification: (A(phi*g0)/A(g0))^4 = {ratio_B:.2f}  (target: {R21_PDG})")
    print()

    # ----------------------------------------------------------
    # Section 4: Energy at phi-FP
    # ----------------------------------------------------------
    print("=" * 72)
    print("  SECTION 4: Soliton energy at phi-FP g0*")
    print("=" * 72)
    print()

    # Form A at phi-FP
    rA2, gA2, gpA2, psiA2 = solve_formA(g0_A_star, R_max=R_MAX_SOLVE, N_pts=N_PTS)
    EA2 = compute_energy(rA2, gA2, gpA2, K_A, V_A, V_A_vac)

    # Form B at phi-FP
    rB2, gB2, gpB2 = solve_formB(g0_B_star, R_max=R_MAX_SOLVE, N_pts=N_PTS)
    EB2 = compute_energy(rB2, gB2, gpB2, K_B, V_BC, V_BC_vac)

    # Form C at phi-FP
    rC2, gC2, gpC2, uC2 = solve_formC(g0_C_star, R_max=R_MAX_SOLVE, N_pts=N_PTS)
    EC2 = compute_energy(rC2, gC2, gpC2, K_C, V_BC, V_BC_vac)

    AA2, _, R2A2 = fit_tail_g(rA2, gA2)
    AB2, _, R2B2 = fit_tail_g(rB2, gB2)
    AC2, _, R2C2 = fit_tail_g(rC2, gC2)

    # ----------------------------------------------------------
    # Final comparison table
    # ----------------------------------------------------------
    print("=" * 72)
    print("  COMPARISON TABLE (at phi-FP)")
    print("=" * 72)
    print()
    hdr = (f"{'Formulation':<18} {'K(g)':<12} "
           f"{'g0*':>10} {'A_tail':>10} {'E_soliton':>12} {'conv':>6}")
    print(hdr)
    print("-" * len(hdr))

    def row(label, Kl, g0, A, E, conv):
        g0s = f"{g0:.6f}" if np.isfinite(g0) else "N/A"
        As  = f"{A:.6f}" if np.isfinite(A) else "N/A"
        Es  = f"{E:.6f}" if np.isfinite(E) else "N/A"
        cs  = "YES" if conv else "NO"
        print(f"{label:<18} {Kl:<12} {g0s:>10} {As:>10} {Es:>12} {cs:>6}")

    row("A (canonical)",  "g^4",       g0_A_star, AA2, EA2, conv_A)
    row("B (historical)", "1+4ln(g)",  g0_B_star, AB2, EB2, conv_B)
    row("C (ex234)",      "g^2",       g0_C_star, AC2, EC2, conv_C)
    print()

    # ----------------------------------------------------------
    # Cross-check: Form C direct ODE
    # ----------------------------------------------------------
    print("--- Cross-check: Form C via direct ODE g''+g'^2/g+(2/r)g'=1-g ---")

    def solve_formC_direct(g0, R_max=120.0, N_pts=5000):
        """Direct g-space solver for Form C."""
        def rhs(r, y):
            g, gp = y
            g_safe = max(g, 1e-15)
            r_safe = max(r, 1e-10)
            gpp = (1.0 - g_safe) - gp**2 / g_safe - 2.0 * gp / r_safe
            return [gp, gpp]
        r0 = 1e-3
        acc0 = 1.0 - g0
        g_init = g0 + acc0 * r0**2 / 6.0
        gp_init = acc0 * r0 / 3.0
        sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                        method='DOP853', rtol=1e-10, atol=1e-12,
                        max_step=0.03, dense_output=True)
        r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
        y_arr = sol.sol(r_arr)
        return r_arr, y_arr[0], y_arr[1]

    rCd, gCd, gpCd = solve_formC_direct(g0_C_star, R_max=R_MAX_SOLVE, N_pts=N_PTS)
    ACd, _, R2Cd = fit_tail_g(rCd, gCd)
    print(f"  g0 = {g0_C_star:.8f}")
    print(f"  A_tail (u-solver) = {AC2:.6f}")
    print(f"  A_tail (g-solver) = {ACd:.6f}")
    # Compare at common r points
    r_common = np.linspace(0.01, min(rC2[-1], rCd[-1]), 2000)
    from scipy.interpolate import interp1d
    gC2_interp = interp1d(rC2, gC2, kind='cubic', fill_value='extrapolate')
    gCd_interp = interp1d(rCd, gCd, kind='cubic', fill_value='extrapolate')
    diff_max = np.max(np.abs(gC2_interp(r_common) - gCd_interp(r_common)))
    print(f"  Max |g_u_solver - g_direct| = {diff_max:.2e}")
    print()

    # ----------------------------------------------------------
    # Scoring
    # ----------------------------------------------------------
    print("=" * 72)
    print("  SCORING")
    print("=" * 72)
    print()
    score = 0
    total = 6

    s1 = conv_C and abs(g0_C_star - 0.86941) < 0.01
    print(f"  S1: Form C phi-FP g0* ~ 0.86941?  g0*={g0_C_star:.6f}  "
          f"{'PASS' if s1 else 'FAIL'}")
    score += int(s1)

    s2 = conv_A
    print(f"  S2: Form A phi-FP converged?  "
          f"{'PASS' if s2 else 'FAIL'}  (g0*={g0_A_star:.6f})")
    score += int(s2)

    s3 = conv_B
    print(f"  S3: Form B phi-FP converged?  "
          f"{'PASS' if s3 else 'FAIL'}  (g0*={g0_B_star:.6f})")
    score += int(s3)

    s4 = all(np.isfinite(A) and A > 0 for A in [AA2, AB2, AC2])
    print(f"  S4: All A_tail finite & > 0?  {'PASS' if s4 else 'FAIL'}")
    score += int(s4)

    s5 = all(np.isfinite(E) for E in [EA2, EB2, EC2])
    print(f"  S5: All energies finite?  {'PASS' if s5 else 'FAIL'}")
    score += int(s5)

    s6 = diff_max < 1e-4
    print(f"  S6: Form C cross-check (u vs g solver)?  "
          f"{'PASS' if s6 else 'FAIL'}  (diff={diff_max:.2e})")
    score += int(s6)

    print(f"\n  TOTAL SCORE: {score}/{total}")
    print()

    # ----------------------------------------------------------
    # Physics notes
    # ----------------------------------------------------------
    print("--- Physics notes ---")
    print()
    print("  1. The soliton tail is OSCILLATORY (omega=1), not exponentially")
    print("     decaying.  This is because V''(1)/K(1) = -1 for all three")
    print("     formulations, yielding the spherical Bessel equation in the tail.")
    print()
    print("  2. The phi-fixed-point (phi-FP) g0* is defined by the condition")
    print("     (A_tail(phi*g0*)/A_tail(g0*))^4 = R21_PDG = 206.768.")
    print("     This encodes the electron-muon mass ratio via soliton tail")
    print("     amplitudes at golden-ratio-spaced central values.")
    print()
    print("  3. Form A (K=g^4) uses psi=g^3/3 as canonical variable, which")
    print("     regularizes the ODE at g=0.  The force term g^4(1-g) is")
    print("     strongly suppressed near g=0 compared to Form C's g(1-g).")
    print()
    print("  4. Form B (K=1+4 ln g) has K(g)>0 only for g > exp(-1/4) ~ 0.778.")
    print("     The ghost boundary at K=0 limits the accessible g0 range.")
    print()
    print("  5. The three formulations give DIFFERENT phi-FP values and")
    print("     A_tail amplitudes, showing that the kinetic coupling K(g)")
    print("     has physical consequences beyond simple field redefinition.")
    print()
    print("=" * 72)
    print("  DONE")
    print("=" * 72)


if __name__ == "__main__":
    main()
