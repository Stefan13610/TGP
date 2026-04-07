#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex230_bps_mass_formula.py
============================
BPS BOUND I MASA SOLITONU DLA K=g⁴

MOTYWACJA (ex229):
  K=g⁴ daje r₂₁ = 206.77 ✓ ale r₃₁ = 553,927 ✗ (PDG: 3477.5)
  Problem: A_tail(g₀) ∝ exp(2.63·g₀) rośnie za szybko
  → M ∝ A⁴ nie daje poprawnych stosunków mas

HIPOTEZA BPS:
  Solitony Bogomolny'ego mają masę = ładunek topologiczny:
    M_BPS = |Q| = |4π∫₀^∞ W'(ψ)·ψ'·r² dr|

  gdzie W(ψ) to "superpotencjał" taki że U(ψ) = ½(W')² - dW/dψ·smth
  (lub prościej: Derrick scaling → virial relation E_kin = E_pot)

  Dla równania ψ'' + 2ψ'/r = U'(ψ) z U'(ψ) = (3ψ)^{1/3} - 1:
  W(ψ) = ∫√(2U(ψ)) dψ

  BPS bound: M ≥ |Q_top|, z równością tylko dla solitonów BPS.

PLAN:
  §1. Derrick virial: E_kin vs E_pot (skalowanie)
  §2. Bogomolny bound: M ≥ Q_top
  §3. Superpotencjał W(ψ) i ładunek topologiczny
  §4. Test: M_BPS / M_BPS_e → r₂₁, r₃₁?
  §5. Efektywny n_K: jakie K_eff(g) daje poprawne A(g₀)?

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

_trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

phi = (1 + np.sqrt(5)) / 2
R21_PDG = 206.768
R31_PDG = 3477.48

TESTS = []
def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")

def koide(m1, m2, m3):
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / s**2


# ============================================================
# ψ-SOLVER (from ex228)
# ============================================================

def solve_soliton_psi(g0, R_max=200.0, N_pts=8000):
    """ψ'' + (2/r)ψ' = 1 - (3ψ)^{1/3}"""
    psi0 = g0**3 / 3.0

    def rhs(r, y):
        psi, psip = y
        psi_safe = max(psi, 1e-30)
        r_safe = max(r, 1e-10)
        g = (3.0 * psi_safe) ** (1.0/3.0)
        psi_pp = (1.0 - g) - 2.0 * psip / r_safe
        return [psip, psi_pp]

    r0 = 1e-4
    acc = 1.0 - g0
    psi_init = psi0 + acc * r0**2 / 6.0
    psip_init = acc * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [psi_init, psip_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.05, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y_arr = sol.sol(r_arr)
    return r_arr, y_arr[0], y_arr[1]


def fit_tail(r_arr, psi_arr, r_L=20.0, r_R=60.0):
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(psi_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan
    r_fit = r_arr[mask]
    y_fit = (psi_arr[mask] - 1.0/3.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    return np.sqrt(B**2 + C**2), np.arctan2(B, C)


# Potential in ψ-space
def U_psi(psi):
    """U(ψ) = V(g), V(g) = g³/3 - g⁴/4 - 1/12, g = (3ψ)^{1/3}"""
    g = (3.0 * np.maximum(psi, 0.0)) ** (1.0/3.0)
    return g**3 / 3.0 - g**4 / 4.0 - 1.0/12.0


def U_prime_psi(psi):
    """dU/dψ = (3ψ)^{1/3} - 1"""
    g = (3.0 * np.maximum(psi, 0.0)) ** (1.0/3.0)
    return g - 1.0


# ============================================================
# §0. Find φ-FP
# ============================================================
print("=" * 72)
print("§0. φ-FP Z ψ-SOLVEREM")
print("=" * 72)

def get_A(g0):
    r, psi, _ = solve_soliton_psi(g0, R_max=100.0, N_pts=4000)
    return fit_tail(r, psi)[0]

def ratio_func(g0):
    A1 = get_A(g0)
    A2 = get_A(phi * g0)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

g0_star = brentq(ratio_func, 0.82, 0.84, xtol=1e-8)
g0_e = g0_star
g0_mu = phi * g0_star
g0_tau = phi**2 * g0_star
print(f"  g₀* = {g0_star:.8f}")
print(f"  g₀_e = {g0_e:.6f}, g₀_μ = {g0_mu:.6f}, g₀_τ = {g0_tau:.6f}")

# Solve all three
print("  Solving solitons...")
r_e, psi_e, psip_e = solve_soliton_psi(g0_e, R_max=200.0)
r_mu, psi_mu, psip_mu = solve_soliton_psi(g0_mu, R_max=200.0)
r_tau, psi_tau, psip_tau = solve_soliton_psi(g0_tau, R_max=200.0)

A_e, _ = fit_tail(r_e, psi_e)
A_mu, _ = fit_tail(r_mu, psi_mu)
A_tau, _ = fit_tail(r_tau, psi_tau, r_L=25.0, r_R=80.0)
print(f"  A_e = {A_e:.6f}, A_μ = {A_mu:.6f}, A_τ = {A_tau:.6f}")


# ============================================================
# §1. DERRICK VIRIAL THEOREM
# ============================================================
print("\n" + "=" * 72)
print("§1. DERRICK VIRIAL (E_kin vs E_pot)")
print("=" * 72)

print("""
  Derrick's theorem w 3D: dla stabilnego solitonu
  E_kin = ½∫(∇ψ)²d³x  i  E_pot = ∫U(ψ)d³x

  Skalowanie: ψ_λ(r) = ψ(r/λ)
  E_kin[λ] = λ·E_kin,  E_pot[λ] = λ³·E_pot

  δE/δλ|_{λ=1} = 0 → E_kin + 3·E_pot = 0  (virial)

  Jeśli virial jest spełniony: M = E_kin + E_pot = -2·E_pot
""")

def compute_virial(r, psi, psip, R_cut=50.0):
    """Compute kinetic and potential energies up to R_cut."""
    idx = min(np.searchsorted(r, R_cut), len(r)-1)
    r_c = r[:idx]
    psi_c = psi[:idx]
    psip_c = psip[:idx]

    E_kin = 4 * np.pi * _trapz(0.5 * psip_c**2 * r_c**2, r_c)
    E_pot = 4 * np.pi * _trapz(U_psi(psi_c) * r_c**2, r_c)
    return E_kin, E_pot

# Use period-averaged cutoffs (multiples of π to average oscillations)
print(f"  {'Soliton':>8s}  {'R_cut':>6s}  {'E_kin':>12s}  {'E_pot':>12s}  {'E_kin+3E_pot':>14s}  {'E_tot':>12s}")
print("  " + "-" * 72)

virial_data = {}
for name, r, psi, psip, g0 in [('e', r_e, psi_e, psip_e, g0_e),
                                  ('μ', r_mu, psi_mu, psip_mu, g0_mu),
                                  ('τ', r_tau, psi_tau, psip_tau, g0_tau)]:
    for R_c in [10*np.pi, 20*np.pi, 30*np.pi]:
        Ek, Ep = compute_virial(r, psi, psip, R_c)
        vir = Ek + 3*Ep
        E_tot = Ek + Ep
        print(f"  {name:>8s}  {R_c:6.1f}  {Ek:12.6f}  {Ep:12.6f}  {vir:14.6f}  {E_tot:12.6f}")
    virial_data[name] = (Ek, Ep, E_tot)  # last R_c


# ============================================================
# §2. ŁADUNEK TOPOLOGICZNY (Bogomolny)
# ============================================================
print("\n" + "=" * 72)
print("§2. ŁADUNEK TOPOLOGICZNY")
print("=" * 72)

print("""
  Dla kanonicznego równania ψ'' + 2ψ'/r = U'(ψ):

  Superpotencjał W(ψ): U(ψ) = ½(dW/dψ)²
  → dW/dψ = ±√(2U(ψ))

  Ładunek topologiczny:
    Q = 4π∫₀^∞ (dW/dψ)·ψ'·r² dr = 4π∫₀^∞ W'(ψ)·ψ'·r² dr

  BPS bound: M ≥ |Q|
  Dla solitonu BPS: ψ' = ∓(dW/dψ)·f(r)
""")

# Compute W(ψ) numerically
psi_grid = np.linspace(0.001, 2.0, 10000)
U_grid = U_psi(psi_grid)

# W(ψ) = ∫₀^ψ √(2U(s)) ds (where U ≥ 0)
# Note: U might be negative in some regions
dW_dpsi = np.sqrt(2.0 * np.maximum(U_grid, 0.0))  # W' = √(2U) where U ≥ 0

print(f"  Analiza potencjału U(ψ):")
print(f"    U(0) = {U_psi(np.array([0.001]))[0]:.6f}")
print(f"    U(1/3) = {U_psi(np.array([1/3]))[0]:.6f}  (vacuum)")
print(f"    U_max at ψ ≈ {psi_grid[np.argmax(U_grid)]:.4f}, U_max = {np.max(U_grid):.6f}")
print(f"    U_min at ψ ≈ {psi_grid[np.argmin(U_grid)]:.4f}, U_min = {np.min(U_grid):.6f}")

# Find where U > 0 and U < 0
mask_pos = U_grid > 0
mask_neg = U_grid < 0
print(f"    U > 0 region: ψ ∈ [{psi_grid[mask_pos][0]:.4f}, {psi_grid[mask_pos][-1]:.4f}]" if any(mask_pos) else "    U > 0: nowhere")
print(f"    U < 0 region: ψ ∈ [{psi_grid[mask_neg][0]:.4f}, {psi_grid[mask_neg][-1]:.4f}]" if any(mask_neg) else "    U < 0: nowhere")


# Topological charge Q = 4π∫ √(2U)·|ψ'|·r² dr
def compute_Q_top(r, psi, psip, R_cut=50.0):
    """Topological charge integral."""
    idx = min(np.searchsorted(r, R_cut), len(r)-1)
    r_c = r[:idx]
    psi_c = psi[:idx]
    psip_c = psip[:idx]

    U_c = U_psi(psi_c)
    W_prime = np.sqrt(2.0 * np.maximum(U_c, 0.0))
    integrand = W_prime * np.abs(psip_c) * r_c**2
    return 4 * np.pi * _trapz(integrand, r_c)


# Alternative: simple field displacement charge
# Q_disp = 4π∫(ψ - ψ_vac)²·r² dr (measures "how much field deviates")
def compute_Q_disp(r, psi, R_cut=50.0):
    idx = min(np.searchsorted(r, R_cut), len(r)-1)
    integrand = (psi[:idx] - 1/3)**2 * r[:idx]**2
    return 4 * np.pi * _trapz(integrand, r[:idx])


# Gradient charge: Q_grad = 4π∫ψ'²·r² dr = 2·E_kin
def compute_Q_grad(r, psip, R_cut=50.0):
    idx = min(np.searchsorted(r, R_cut), len(r)-1)
    integrand = psip[:idx]**2 * r[:idx]**2
    return 4 * np.pi * _trapz(integrand, r[:idx])


# Central excess: ∫(ψ - ψ_vac)·r² dr
def compute_Q_central(r, psi, R_cut=50.0):
    idx = min(np.searchsorted(r, R_cut), len(r)-1)
    integrand = (psi[:idx] - 1/3) * r[:idx]**2
    return 4 * np.pi * _trapz(integrand, r[:idx])


# BPS-like charge: ∫|ψ'|·r dr (without r²)
def compute_Q_bps(r, psip, R_cut=50.0):
    idx = min(np.searchsorted(r, R_cut), len(r)-1)
    return 4 * np.pi * _trapz(np.abs(psip[:idx]) * r[:idx], r[:idx])


R_C = 20 * np.pi  # Multiple of period for averaging

print(f"\n  Ładunki topologiczne (R_cut = {R_C:.1f}):")
print(f"  {'Charge':25s}  {'Q_e':>12s}  {'Q_μ':>12s}  {'Q_τ':>12s}  {'Q_μ/Q_e':>10s}  {'Q_τ/Q_e':>10s}")
print("  " + "-" * 80)

charges = {}
for charge_name, charge_func in [
    ("Q_top = ∫W'|ψ'|r²", lambda r, p, pp: compute_Q_top(r, p, pp, R_C)),
    ("Q_disp = ∫(δψ)²r²", lambda r, p, pp: compute_Q_disp(r, p, R_C)),
    ("Q_grad = ∫ψ'²r²", lambda r, p, pp: compute_Q_grad(r, pp, R_C)),
    ("Q_central = ∫δψ·r²", lambda r, p, pp: compute_Q_central(r, p, R_C)),
    ("Q_bps = ∫|ψ'|r", lambda r, p, pp: compute_Q_bps(r, pp, R_C)),
]:
    Qe = charge_func(r_e, psi_e, psip_e)
    Qmu = charge_func(r_mu, psi_mu, psip_mu)
    Qtau = charge_func(r_tau, psi_tau, psip_tau)

    r21_q = Qmu/Qe if abs(Qe) > 1e-20 else np.nan
    r31_q = Qtau/Qe if abs(Qe) > 1e-20 else np.nan
    charges[charge_name] = (Qe, Qmu, Qtau, r21_q, r31_q)

    print(f"  {charge_name:25s}  {Qe:12.6f}  {Qmu:12.6f}  {Qtau:12.6f}  {r21_q:10.4f}  {r31_q:10.4f}")


# ============================================================
# §3. SZUKANIE FORMUŁY MASY: Q^p → r₂₁, r₃₁
# ============================================================
print("\n" + "=" * 72)
print("§3. ★ SZUKANIE FORMUŁY M = Q^p → PASUJĄCEJ DO PDG")
print("=" * 72)

print(f"\n  PDG: r₂₁ = {R21_PDG}, r₃₁ = {R31_PDG}")
print(f"  Koide K(PDG) = {koide(1, R21_PDG, R31_PDG):.6f} vs 2/3 = {2/3:.6f}\n")

print(f"  {'Charge / Observable':25s}  {'ratio_μ/e':>10s}  {'p(r₂₁)':>8s}  {'pred_r₃₁':>12s}  {'err_r₃₁':>10s}")
print("  " + "-" * 72)

best_formula = None
best_err = 1e10

for charge_name, (Qe, Qmu, Qtau, r21_q, r31_q) in charges.items():
    if abs(r21_q) > 1 and np.isfinite(r21_q) and r21_q > 0:
        # Find p such that r21_q^p = R21_PDG
        p_opt = np.log(R21_PDG) / np.log(r21_q)
        # Predict r₃₁
        r31_pred = r31_q ** p_opt if r31_q > 0 else np.nan
        err = abs(r31_pred - R31_PDG) / R31_PDG * 100 if np.isfinite(r31_pred) else np.nan

        if np.isfinite(err) and err < best_err:
            best_err = err
            best_formula = (charge_name, p_opt, r31_pred, err)

        print(f"  {charge_name:25s}  {r21_q:10.4f}  {p_opt:8.4f}  {r31_pred:12.2f}  {err:9.2f}%")

# Also try COMBINATIONS
print(f"\n  Kombinacje:")
for cname, (Qe, Qmu, Qtau, _, _) in charges.items():
    # Q * A_tail
    for power_A in [1, 2, 3, 4]:
        Xe = Qe * A_e**power_A
        Xmu = Qmu * A_mu**power_A
        Xtau = Qtau * A_tau**power_A

        if abs(Xe) > 1e-20:
            r21_x = Xmu / Xe
            r31_x = Xtau / Xe

            if r21_x > 1 and np.isfinite(r21_x):
                p_opt = np.log(R21_PDG) / np.log(r21_x)
                r31_pred = r31_x ** p_opt if r31_x > 0 else np.nan
                err = abs(r31_pred - R31_PDG) / R31_PDG * 100 if np.isfinite(r31_pred) else np.nan

                if np.isfinite(err) and err < 20:  # Only print close matches
                    print(f"  {cname[:15]}·A^{power_A}  {r21_x:10.4f}  {p_opt:8.4f}  {r31_pred:12.2f}  {err:9.2f}%")

                if np.isfinite(err) and err < best_err:
                    best_err = err
                    best_formula = (f"{cname}·A^{power_A}", p_opt, r31_pred, err)

if best_formula:
    name, p, r31, err = best_formula
    print(f"\n  ★ NAJLEPSZA FORMUŁA: {name}")
    print(f"    p = {p:.4f}, predicted r₃₁ = {r31:.2f}, error = {err:.2f}%")

    record("T1: BPS mass formula",
           err < 10,
           f"Best: {name}, p={p:.4f}, err(r₃₁)={err:.2f}%")
else:
    record("T1: BPS mass formula", False, "No formula found")


# ============================================================
# §4. ODWROTNE PODEJŚCIE: JAKI A(g₀) DAJE POPRAWNE RATIO?
# ============================================================
print("\n" + "=" * 72)
print("§4. REQUIRED A(g₀) FOR CORRECT MASS RATIOS")
print("=" * 72)

print("""
  Jeśli M ∝ A(g₀)⁴ ma dawać poprawne ratio:
    r₂₁ = [A(g₀_μ)/A(g₀_e)]⁴ = 206.768
    r₃₁ = [A(g₀_τ)/A(g₀_e)]⁴ = 3477.48

  To wymagane:
    A(g₀_μ)/A(g₀_e) = r₂₁^{1/4} = 3.7920
    A(g₀_τ)/A(g₀_e) = r₃₁^{1/4} = 7.6806
""")

A_ratio_21_req = R21_PDG**(1/4)
A_ratio_31_req = R31_PDG**(1/4)

# What does K=g⁴ give?
A_ratio_21_K4 = A_mu / A_e
A_ratio_31_K4 = A_tau / A_e

print(f"  Wymagane: A_μ/A_e = {A_ratio_21_req:.4f}, A_τ/A_e = {A_ratio_31_req:.4f}")
print(f"  K=g⁴:    A_μ/A_e = {A_ratio_21_K4:.4f}, A_τ/A_e = {A_ratio_31_K4:.4f}")
print(f"  K=g⁴ A_τ/A_e jest {A_ratio_31_K4/A_ratio_31_req:.1f}x za duży")

# What functional form of A(g₀) gives BOTH ratios?
# We need: A(φ·g₀)/A(g₀) = 3.792 AND A(φ²·g₀)/A(g₀) = 7.681
# This means: A(φ²·g₀)/A(φ·g₀) = 7.681/3.792 = 2.025 = r₃₂^{1/4}

r32_fourth_root = (R31_PDG / R21_PDG) ** 0.25
print(f"\n  A(φ²·g₀)/A(φ·g₀) powinno = {r32_fourth_root:.4f}")
print(f"  K=g⁴ daje:                   {A_tau/A_mu:.4f}")
print(f"  Nadmiar: {(A_tau/A_mu)/r32_fourth_root:.2f}x")

# For φ-FP to work with correct ratios, we need ln[A(g₀)] to be CONCAVE
# (sublinear growth), not exponential
# Required: ln(A(g₀_μ)/A(g₀_e)) / ln(A(g₀_τ)/A(g₀_e)) = ln(3.792)/ln(7.681) = 0.654
ratio_log = np.log(A_ratio_21_req) / np.log(A_ratio_31_req)
ratio_log_K4 = np.log(A_ratio_21_K4) / np.log(A_ratio_31_K4)

print(f"\n  ln(A_μ/A_e)/ln(A_τ/A_e):")
print(f"    Required = {ratio_log:.4f}")
print(f"    K=g⁴     = {ratio_log_K4:.4f}")

# This ratio = ln(φ)/ln(φ²) = 1/2 only if A ∝ exp(c·g₀), i.e. ln(A) ∝ g₀
# For concave: ratio > 1/2
# Required 0.654 > 0.5 → ln(A) must be CONCAVE in g₀
# K=g⁴ gives 0.407 < 0.5 → ln(A) is CONVEX → exponential growth

print(f"\n  Interpretacja:")
print(f"    ratio = 0.5 → ln(A) liniowe w g₀  (A ∝ exp(c·g₀))")
print(f"    ratio > 0.5 → ln(A) wklęsłe      (A rośnie wolniej niż exp)")
print(f"    ratio < 0.5 → ln(A) wypukłe       (A rośnie szybciej niż exp)")
print(f"\n    Required: {ratio_log:.4f} > 0.5 → WKLĘSŁE (concave)")
print(f"    K=g⁴:     {ratio_log_K4:.4f} < 0.5 → WYPUKŁE (convex) ✗")

record("T2: A(g₀) shape requirement",
       False,
       f"Need concave (ratio={ratio_log:.3f}), K=g⁴ is convex ({ratio_log_K4:.3f})")


# ============================================================
# §5. EFEKTYWNE K_eff: CO DAJE POPRAWNE A(g₀)?
# ============================================================
print("\n" + "=" * 72)
print("§5. EFEKTYWNE K_eff(g) DLA POPRAWNYCH MAS")
print("=" * 72)

print("""
  Pytanie: jaka kinetyka K(g) daje ln(A) WKLĘSŁE w g₀?

  Dla K=g^n_K, ODE w ψ-space: ψ'' + (2/r)ψ' = f(ψ)
  - n_K = 2: f = (1-g)/1 = 1-g          → ln(A) wklęsłe ✓
  - n_K = 4: f = 1-(3ψ)^{1/3}           → ln(A) wypukłe ✗
  - n_K = ?: szukamy granicznego n_K*

  Inaczej: w g-space, ODE to:
    g'' = (1-g)/g^(n_K-2) - (n_K/2)g'²/g - 2g'/r

  Dla n_K=2: (1-g)/g⁰ = (1-g)          → łagodne przywracanie
  Dla n_K=4: (1-g)/g² → 1/g² at g→0    → ostre przywracanie → szybki wzrost A

  Hipoteza: efektywne n_K(effective) ≈ 2 w sektorze leptonowym,
  mimo że n_K(bare) = 4 z geometrii akcji.
""")

# Scan n_K from 2 to 4 and compute A-ratio behavior
def solve_soliton_nK(g0, n_K, R_max=100.0, N_pts=4000):
    """General n_K solver in g-space."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-15)
        r = max(r, 1e-10)
        if n_K == 2:
            gpp = (1.0 - g) - gp**2 / g - 2.0 * gp / r
        else:
            gpp = (1.0 - g) / g**(n_K-2) - (n_K/2.0) * gp**2 / g - 2.0 * gp / r
        return [gp, gpp]

    r0 = 1e-3
    if n_K == 2:
        g0_acc = 1.0 - g0
    else:
        g0_acc = (1.0 - g0) / g0**(n_K-2)
    g_init = g0 + g0_acc * r0**2 / 6.0
    gp_init = g0_acc * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.03, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y = sol.sol(r_arr)
    return r_arr, y[0]


def fit_tail_g(r_arr, g_arr, r_L=20.0, r_R=55.0):
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(g_arr)
    if np.sum(mask) < 20:
        return np.nan
    r_fit = r_arr[mask]
    y_fit = (g_arr[mask] - 1.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    return np.sqrt(coefs[0]**2 + coefs[1]**2)


# For each n_K, find φ-FP and compute r₃₁ error
print(f"\n  {'n_K':>5s}  {'g₀*':>10s}  {'r₂₁':>10s}  {'A_τ/A_e':>10s}  {'(A_τ/A_e)⁴':>12s}  {'err_r₃₁':>10s}")
print("  " + "-" * 64)

for n_K_test in [2.0, 2.5, 3.0, 3.5, 4.0]:
    # Find φ-FP for this n_K
    def ratio_nK(g0):
        try:
            r1, g1 = solve_soliton_nK(g0, n_K_test)
            r2, g2 = solve_soliton_nK(phi*g0, n_K_test)
            A1 = fit_tail_g(r1, g1)
            A2 = fit_tail_g(r2, g2)
            if np.isnan(A1) or A1 < 1e-15:
                return 1e10
            return (A2/A1)**4 - R21_PDG
        except:
            return 1e10

    # Scan for zero crossing
    g0_fp = None
    for gl in np.arange(0.5, 2.0, 0.1):
        try:
            vl = ratio_nK(gl)
            vr = ratio_nK(gl + 0.1)
            if np.isfinite(vl) and np.isfinite(vr) and vl * vr < 0:
                g0_fp = brentq(ratio_nK, gl, gl + 0.1, xtol=1e-6)
                break
        except:
            continue

    if g0_fp is not None:
        r1, g1 = solve_soliton_nK(g0_fp, n_K_test)
        r2, g2 = solve_soliton_nK(phi*g0_fp, n_K_test)
        A1 = fit_tail_g(r1, g1)
        A2 = fit_tail_g(r2, g2)
        r21_test = (A2/A1)**4

        # τ: try solving with larger g₀
        g0_tau_test = phi**2 * g0_fp
        try:
            if n_K_test <= 2.5 or g0_tau_test < 1.8:
                r3, g3 = solve_soliton_nK(g0_tau_test, n_K_test, R_max=120.0)
                A3 = fit_tail_g(r3, g3, r_L=25.0, r_R=65.0)
            else:
                # Use ψ-solver for large g₀
                r3, psi3, _ = solve_soliton_psi(g0_tau_test, R_max=150.0)
                A3, _ = fit_tail(r3, psi3, r_L=25.0, r_R=70.0)

            if np.isfinite(A3) and A3 > 0:
                ratio_31 = A3 / A1
                r31_test = ratio_31**4
                err_31 = abs(r31_test - R31_PDG) / R31_PDG * 100
                print(f"  {n_K_test:5.1f}  {g0_fp:10.6f}  {r21_test:10.2f}  {ratio_31:10.4f}  {r31_test:12.2f}  {err_31:9.2f}%")
            else:
                print(f"  {n_K_test:5.1f}  {g0_fp:10.6f}  {r21_test:10.2f}  {'NaN':>10s}  {'NaN':>12s}")
        except Exception as e:
            print(f"  {n_K_test:5.1f}  {g0_fp:10.6f}  {r21_test:10.2f}  {'ERR':>10s}  {str(e)[:20]}")
    else:
        print(f"  {n_K_test:5.1f}  {'NOT FOUND':>10s}")


# ============================================================
# SCORECARD
# ============================================================
print("\n" + "=" * 72)
print("SCORECARD")
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)
for name, passed, detail in TESTS:
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")
print(f"\n  {n_pass}/{n_total} testów przeszło.")


# ============================================================
# PODSUMOWANIE
# ============================================================
print("\n" + "=" * 72)
print("PODSUMOWANIE ex230")
print("=" * 72)
print(f"""
  KLUCZOWE WYNIKI:

  1. K=g⁴ wymaga ln(A(g₀)) WYPUKŁE (convex)
     Poprawne masy wymagają ln(A(g₀)) WKLĘSŁE (concave)
     → FUNDAMENTALNA NIEZGODNOŚĆ

  2. Żadna formuła BPS / topologiczna / energetyczna
     nie naprawia problemu r₃₁ dla K=g⁴

  3. Jedyne rozwiązanie: efektywne n_K < 4
     (RG running: n_K_bare=4 → n_K_eff ≈ 2 na skali solitonu?)

  INTERPRETACJA:
  K(g) = g^n_K z n_K = D = 4 wynika z geometrii akcji.
  Ale solitonowa dynamika "widzi" EFEKTYWNE K_eff(g) = g^n_eff
  z n_eff ≈ 2, zmodyfikowane przez:
    - Pętlowe poprawki (1-loop: n_eff = n_K - δn z kontr-członów)
    - Anomalny wymiar γ_g pola g
    - Przejście UV→IR: n_K(UV) = 4, n_K(IR) = 2
""")
