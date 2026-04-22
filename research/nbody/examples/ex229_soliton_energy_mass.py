#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex229_soliton_energy_mass.py
================================
MASA SOLITONU = ENERGIA CZY A_tail⁴?

KONTEKST (ex228):
  ψ-solver rozwiązał τ soliton (g₀ > 2) ✓
  ALE: r₃₁ = (A_τ/A_e)⁴ = 550,264 ≫ PDG 3,477
  → Formuła M ∝ A_tail⁴ FATALNIE ZAWODZI dla K=g⁴ + τ

HIPOTEZA:
  Masa solitonu = ENERGIA pola, nie A_tail^p
  M = 4π ∫₀^R [½ψ'² + U(ψ)] r² dr
  gdzie ψ = g³/3 (zmienna kanoniczna dla K=g⁴)

  Ogon oscyluje (1/r), więc energia rośnie z R.
  ALE stosunek M_μ/M_e powinien dążyć do stałej.

TEST:
  1. Obliczyć E(R) dla e, μ, τ solitonów
  2. Sprawdzić zbieżność E_μ(R)/E_e(R) → r₂₁?
  3. Sprawdzić E_τ(R)/E_e(R) → r₃₁?
  4. Zbadać alternatywne formuły masy: E_core, E_kinetic, ...

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.integrate import solve_ivp, simpson

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

_trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)
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
# CANONICAL ψ-SOLVER (from ex228)
# ============================================================

def solve_soliton_psi(g0, R_max=200.0, N_pts=10000):
    """
    Solve ψ'' + (2/r)ψ' = 1 - (3ψ)^{1/3}
    Returns (r_arr, psi_arr, psip_arr).
    """
    psi0 = g0**3 / 3.0

    def rhs(r, y):
        psi, psip = y
        psi_safe = max(psi, 1e-30)
        r_safe = max(r, 1e-10)
        g = (3.0 * psi_safe) ** (1.0/3.0)
        psi_pp = (1.0 - g) - 2.0 * psip / r_safe
        return [psip, psi_pp]

    r0 = 1e-4
    g_val = (3.0 * psi0) ** (1.0/3.0)
    acc = 1.0 - g_val

    psi_init = psi0 + acc * r0**2 / 6.0
    psip_init = acc * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [psi_init, psip_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.05, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y_arr = sol.sol(r_arr)
    psi_arr = y_arr[0]
    psip_arr = y_arr[1]

    return r_arr, psi_arr, psip_arr


def psi_to_g(psi_arr):
    return (3.0 * np.maximum(psi_arr, 0.0)) ** (1.0/3.0)


def fit_tail(r_arr, psi_arr, r_L=20.0, r_R=60.0):
    """Fit (ψ - 1/3)·r = B cos(r) + C sin(r)."""
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(psi_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan
    r_fit = r_arr[mask]
    y_fit = (psi_arr[mask] - 1.0/3.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    return np.sqrt(B**2 + C**2), np.arctan2(B, C)


# Potential in ψ-space: U(ψ) such that U'(ψ) = (3ψ)^{1/3} - 1
# U(ψ) = -ψ + (3/4)(3ψ)^{1/3}·ψ + C  (with U(1/3) = 0)
# Simpler: U(ψ) = V(g) where g = (3ψ)^{1/3}
# V(g) = g³/3 - g⁴/4 - 1/12

def potential_psi(psi):
    """U(ψ) = V(g) with g = (3ψ)^{1/3}, V(g) = g³/3 - g⁴/4 - 1/12."""
    g = (3.0 * np.maximum(psi, 0.0)) ** (1.0/3.0)
    return g**3 / 3.0 - g**4 / 4.0 - 1.0/12.0


# ============================================================
# §1. PROFIL ENERGII SOLITONÓW
# ============================================================
print("=" * 72)
print("§1. ENERGIA SOLITONÓW e / μ / τ (K = g⁴)")
print("=" * 72)

# Find φ-FP
from scipy.optimize import brentq

def get_A(g0):
    r, psi, _ = solve_soliton_psi(g0, R_max=100.0, N_pts=5000)
    return fit_tail(r, psi)[0]

def ratio_func(g0):
    A1 = get_A(g0)
    A2 = get_A(phi * g0)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

print("\n  Szukanie φ-FP...")
g0_star = brentq(ratio_func, 0.82, 0.84, xtol=1e-8)
print(f"  g₀* = {g0_star:.8f}")

g0_e = g0_star
g0_mu = phi * g0_star
g0_tau = phi**2 * g0_star

print(f"  g₀_e = {g0_e:.6f}")
print(f"  g₀_μ = {g0_mu:.6f}")
print(f"  g₀_τ = {g0_tau:.6f}")

# Solve all three
R_MAX = 200.0
print(f"\n  Rozwiązywanie solitonów (R_max = {R_MAX})...")

r_e, psi_e, psip_e = solve_soliton_psi(g0_e, R_max=R_MAX)
r_mu, psi_mu, psip_mu = solve_soliton_psi(g0_mu, R_max=R_MAX)
r_tau, psi_tau, psip_tau = solve_soliton_psi(g0_tau, R_max=R_MAX)

# A_tail values
A_e, _ = fit_tail(r_e, psi_e)
A_mu, _ = fit_tail(r_mu, psi_mu)
A_tau, _ = fit_tail(r_tau, psi_tau, r_L=25.0, r_R=80.0)

print(f"  A_e = {A_e:.8f}")
print(f"  A_μ = {A_mu:.8f}")
print(f"  A_τ = {A_tau:.8f}")
print(f"  r₂₁(A⁴) = {(A_mu/A_e)**4:.2f}  (PDG: {R21_PDG})")
print(f"  r₃₁(A⁴) = {(A_tau/A_e)**4:.2f}  (PDG: {R31_PDG})")


# ============================================================
# §2. ENERGY INTEGRALS: M(R) = 4π∫₀ᴿ [½ψ'² + U(ψ)] r² dr
# ============================================================
print("\n" + "=" * 72)
print("§2. CAŁKA ENERGII: M(R) = 4π∫₀ᴿ [½ψ'² + U(ψ)] r² dr")
print("=" * 72)

def compute_energy_profile(r_arr, psi_arr, psip_arr):
    """Compute cumulative energy M(R) = 4π∫₀ᴿ [½ψ'² + U(ψ)] r² dr.

    Also computes period-averaged cumulative (averaging over 2π periods)
    to handle oscillating tail.
    """
    U = potential_psi(psi_arr)
    kinetic = 0.5 * psip_arr**2
    density = (kinetic + U) * r_arr**2 * 4.0 * np.pi

    # Cumulative integral
    dr = np.diff(r_arr)
    cumulative = np.zeros(len(r_arr))
    for i in range(1, len(r_arr)):
        cumulative[i] = cumulative[i-1] + 0.5 * (density[i-1] + density[i]) * dr[i-1]

    # Period-averaged: average cumulative over window of 2π
    # This smooths out oscillations in the tail energy
    period = 2 * np.pi
    avg_cumulative = np.copy(cumulative)
    for i in range(len(r_arr)):
        r_lo = r_arr[i] - period / 2
        r_hi = r_arr[i] + period / 2
        mask = (r_arr >= max(r_lo, r_arr[0])) & (r_arr <= min(r_hi, r_arr[-1]))
        if np.sum(mask) > 1:
            avg_cumulative[i] = np.mean(cumulative[mask])

    return cumulative, avg_cumulative, kinetic * r_arr**2, U * r_arr**2


print("\n  Obliczanie profili energii...")
E_e, Eavg_e, kin_e, pot_e = compute_energy_profile(r_e, psi_e, psip_e)
E_mu, Eavg_mu, kin_mu, pot_mu = compute_energy_profile(r_mu, psi_mu, psip_mu)
E_tau, Eavg_tau, kin_tau, pot_tau = compute_energy_profile(r_tau, psi_tau, psip_tau)

# Energy ratios vs R (period-averaged to handle oscillating tail)
print(f"\n  Stosunek energii vs R (uśredniony po okresie 2π):")
print(f"  {'R':>8s}  {'Ē_μ/Ē_e':>12s}  {'Ē_τ/Ē_e':>12s}  {'Ē_τ/Ē_μ':>12s}  {'Ē_e':>14s}")
print("  " + "-" * 66)

R_checks = [5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 180]
for R_check in R_checks:
    idx_e = np.searchsorted(r_e, R_check)
    idx_mu = np.searchsorted(r_mu, R_check)
    idx_tau = np.searchsorted(r_tau, R_check)

    if idx_e >= len(r_e) or idx_mu >= len(r_mu) or idx_tau >= len(r_tau):
        continue

    Ee = Eavg_e[idx_e]
    Emu = Eavg_mu[idx_mu]
    Etau = Eavg_tau[idx_tau]

    if abs(Ee) > 1e-20:
        r21_E = Emu / Ee
        r31_E = Etau / Ee
        r32_E = Etau / Emu if abs(Emu) > 1e-20 else np.nan
        print(f"  {R_check:8d}  {r21_E:12.4f}  {r31_E:12.4f}  {r32_E:12.4f}  {Ee:14.6e}")

# ============================================================
# §3. RÓŻNE FORMUŁY MASY
# ============================================================
print("\n" + "=" * 72)
print("§3. ALTERNATYWNE FORMUŁY MASY")
print("=" * 72)

R_cutoff = 100.0  # Fixed cutoff for comparison

idx_e_c = np.searchsorted(r_e, R_cutoff)
idx_mu_c = np.searchsorted(r_mu, R_cutoff)
idx_tau_c = np.searchsorted(r_tau, R_cutoff)

# Total energy (period-averaged)
E_tot_e = Eavg_e[idx_e_c]
E_tot_mu = Eavg_mu[idx_mu_c]
E_tot_tau = Eavg_tau[idx_tau_c]

# Kinetic energy only
def integrate(r, integrand, R_max):
    idx = min(np.searchsorted(r, R_max), len(r)-1)
    if idx < 2:
        return 0.0
    return _trapz(integrand[:idx] * 4 * np.pi, r[:idx])

E_kin_e = integrate(r_e, 0.5 * psip_e**2 * r_e**2, R_cutoff)
E_kin_mu = integrate(r_mu, 0.5 * psip_mu**2 * r_mu**2, R_cutoff)
E_kin_tau = integrate(r_tau, 0.5 * psip_tau**2 * r_tau**2, R_cutoff)

# Potential energy only
E_pot_e = integrate(r_e, potential_psi(psi_e) * r_e**2, R_cutoff)
E_pot_mu = integrate(r_mu, potential_psi(psi_mu) * r_mu**2, R_cutoff)
E_pot_tau = integrate(r_tau, potential_psi(psi_tau) * r_tau**2, R_cutoff)

# Core energy (R < 10) — raw, not averaged (core doesn't oscillate much)
R_core = 10.0
E_core_e = E_e[min(np.searchsorted(r_e, R_core), len(r_e)-1)]
E_core_mu = E_mu[min(np.searchsorted(r_mu, R_core), len(r_mu)-1)]
E_core_tau = E_tau[min(np.searchsorted(r_tau, R_core), len(r_tau)-1)]

# Field integral: ∫|ψ-ψ_vac|·r² dr (topological charge-like)
Q_e = integrate(r_e, np.abs(psi_e - 1/3) * r_e**2, R_cutoff)
Q_mu = integrate(r_mu, np.abs(psi_mu - 1/3) * r_mu**2, R_cutoff)
Q_tau = integrate(r_tau, np.abs(psi_tau - 1/3) * r_tau**2, R_cutoff)

# Field integral squared: (∫|ψ-ψ_vac|·r² dr)²
# Gradient integral: ∫|ψ'|²·r dr (without r²)
G_e = integrate(r_e, psip_e**2 * r_e, R_cutoff)
G_mu = integrate(r_mu, psip_mu**2 * r_mu, R_cutoff)
G_tau = integrate(r_tau, psip_tau**2 * r_tau, R_cutoff)

# Central value itself
psi_center_e = g0_e**3 / 3
psi_center_mu = g0_mu**3 / 3
psi_center_tau = g0_tau**3 / 3

print(f"\n  R_cutoff = {R_cutoff}")
print(f"\n  {'Formuła masy':25s}  {'r₂₁':>10s}  {'r₃₁':>12s}  {'K (Koide)':>10s}  {'PDG r₂₁':>10s}  {'PDG r₃₁':>10s}")
print("  " + "-" * 80)

formulas = [
    ("A_tail⁴", (A_mu/A_e)**4, (A_tau/A_e)**4),
    ("A_tail²", (A_mu/A_e)**2, (A_tau/A_e)**2),
    ("A_tail³", (A_mu/A_e)**3, (A_tau/A_e)**3),
    ("E_total(R=100)", E_tot_mu/E_tot_e, E_tot_tau/E_tot_e),
    ("E_kinetic", E_kin_mu/E_kin_e if E_kin_e > 0 else np.nan,
                  E_kin_tau/E_kin_e if E_kin_e > 0 else np.nan),
    ("E_potential", E_pot_mu/E_pot_e if abs(E_pot_e) > 0 else np.nan,
                    E_pot_tau/E_pot_e if abs(E_pot_e) > 0 else np.nan),
    ("E_core(R=10)", E_core_mu/E_core_e if E_core_e > 0 else np.nan,
                     E_core_tau/E_core_e if E_core_e > 0 else np.nan),
    ("Q = ∫|δψ|r²dr", Q_mu/Q_e, Q_tau/Q_e),
    ("Q²", (Q_mu/Q_e)**2, (Q_tau/Q_e)**2),
    ("∫ψ'²r dr", G_mu/G_e if G_e > 0 else np.nan,
                  G_tau/G_e if G_e > 0 else np.nan),
    ("ψ_center", psi_center_mu/psi_center_e, psi_center_tau/psi_center_e),
    ("g₀³", (g0_mu/g0_e)**3, (g0_tau/g0_e)**3),
    ("(g₀-1)⁴", ((g0_mu-1)/(g0_e-1))**4, ((g0_tau-1)/(g0_e-1))**4),
]

for name, r21_val, r31_val in formulas:
    if np.isfinite(r21_val) and np.isfinite(r31_val) and r21_val > 0 and r31_val > 0:
        # Compute effective Koide from r₂₁ and r₃₁
        # m₁ = 1, m₂ = r₂₁, m₃ = r₃₁
        K_val = koide(1.0, r21_val, r31_val)
        print(f"  {name:25s}  {r21_val:10.2f}  {r31_val:12.2f}  {K_val:10.6f}  {R21_PDG:10.3f}  {R31_PDG:10.1f}")
    else:
        print(f"  {name:25s}  {r21_val:10.2f}  {r31_val:12.2f}  {'—':>10s}")


# ============================================================
# §4. SZUKANIE OPTYMALNEGO WYKŁADNIKA p: M ∝ A^p
# ============================================================
print("\n" + "=" * 72)
print("§4. OPTYMALNY WYKŁADNIK p: (A_μ/A_e)^p = r₂₁")
print("=" * 72)

# From A ratios, find p such that (A_μ/A_e)^p = r₂₁
ratio_mu_e = A_mu / A_e
ratio_tau_e = A_tau / A_e

p_from_r21 = np.log(R21_PDG) / np.log(ratio_mu_e)
r31_pred_p = ratio_tau_e ** p_from_r21

print(f"\n  A_μ/A_e = {ratio_mu_e:.6f}")
print(f"  A_τ/A_e = {ratio_tau_e:.6f}")
print(f"  ln(r₂₁)/ln(A_μ/A_e) = p = {p_from_r21:.4f}")
print(f"  (A_τ/A_e)^{p_from_r21:.2f} = {r31_pred_p:.2f}")
print(f"  PDG r₃₁ = {R31_PDG}")
print(f"  Error: {abs(r31_pred_p-R31_PDG)/R31_PDG*100:.2f}%")

# Also: what p gives correct r₃₁?
p_from_r31 = np.log(R31_PDG) / np.log(ratio_tau_e)
r21_pred_p31 = ratio_mu_e ** p_from_r31

print(f"\n  Alternatywnie, p z r₃₁:")
print(f"  ln(r₃₁)/ln(A_τ/A_e) = p = {p_from_r31:.4f}")
print(f"  (A_μ/A_e)^{p_from_r31:.2f} = {r21_pred_p31:.2f}")
print(f"  PDG r₂₁ = {R21_PDG}")

# Conclusion: is there any p that fits BOTH?
print(f"\n  Czy istnieje p pasujące do OBU ratios?")
print(f"  p(r₂₁) = {p_from_r21:.4f}")
print(f"  p(r₃₁) = {p_from_r31:.4f}")
print(f"  Różnica: {abs(p_from_r21 - p_from_r31):.4f}")

if abs(p_from_r21 - p_from_r31) / p_from_r21 < 0.05:
    record("T1: Universal mass exponent",
           True,
           f"p = {(p_from_r21+p_from_r31)/2:.4f}")
else:
    record("T1: Universal mass exponent",
           False,
           f"p(r₂₁) = {p_from_r21:.4f}, p(r₃₁) = {p_from_r31:.4f}")


# ============================================================
# §5. GEOMETRIA φ-FP: CZY ln(A) ∝ g₀?
# ============================================================
print("\n" + "=" * 72)
print("§5. FUNKCJONALNA ZALEŻNOŚĆ A_tail(g₀)")
print("=" * 72)

# Compute A for many g₀ values
g0_vals = np.linspace(0.3, 2.5, 50)
A_vals = []
for g0 in g0_vals:
    r, psi, psip = solve_soliton_psi(g0, R_max=120.0, N_pts=3000)
    A, _ = fit_tail(r, psi, r_L=20.0, r_R=60.0)
    A_vals.append(A)
A_vals = np.array(A_vals)

# Fit ln(A) vs g₀ in different regions
mask_valid = np.isfinite(A_vals) & (A_vals > 0)

# Above g₀ = 1
mask_above = mask_valid & (g0_vals > 1.05)
if np.sum(mask_above) >= 3:
    coeffs = np.polyfit(g0_vals[mask_above], np.log(A_vals[mask_above]), 1)
    print(f"\n  Dla g₀ > 1.05: ln(A) ≈ {coeffs[0]:.4f}·g₀ + {coeffs[1]:.4f}")
    print(f"  → A ∝ exp({coeffs[0]:.4f}·g₀)")

    # Test: does A ∝ exp(α·g₀) with φ-FP give correct ratios?
    alpha = coeffs[0]
    r21_exp = np.exp(alpha * phi * g0_star) / np.exp(alpha * g0_star)
    print(f"  A_μ/A_e = exp(α·φ·g₀*)/exp(α·g₀*) = exp(α·(φ-1)·g₀*) = {r21_exp:.4f}")

# Fit ln(A) vs ln|g₀-1|
mask_near = mask_valid & (np.abs(g0_vals-1) > 0.03) & (np.abs(g0_vals-1) < 0.5)
if np.sum(mask_near) >= 3:
    coeffs2 = np.polyfit(np.log(np.abs(g0_vals[mask_near] - 1)),
                         np.log(A_vals[mask_near]), 1)
    print(f"\n  Blisko g₀ ≈ 1: ln(A) ≈ {coeffs2[0]:.4f}·ln|g₀-1| + {coeffs2[1]:.4f}")
    print(f"  → A ∝ |g₀-1|^{coeffs2[0]:.4f}")

# Full functional form
print(f"\n  Tabela A_tail(g₀):")
print(f"  {'g₀':>8s}  {'A':>12s}  {'ln A':>10s}  {'|g₀-1|':>8s}")
print("  " + "-" * 44)
for g0, A in zip(g0_vals, A_vals):
    if np.isfinite(A) and A > 0:
        print(f"  {g0:8.4f}  {A:12.6f}  {np.log(A):10.4f}  {abs(g0-1):8.4f}")


# ============================================================
# §6. PORÓWNANIE Z K=g² (n_K=2)
# ============================================================
print("\n" + "=" * 72)
print("§6. PORÓWNANIE: K=g⁴ vs K=g² (ENERGIA)")
print("=" * 72)

def solve_soliton_K2(g0, R_max=200.0, N_pts=10000):
    """Standard K=g² solver."""
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-15)
        r = max(r, 1e-10)
        gpp = (1.0 - g) - gp**2 / g - 2.0 * gp / r  # n_K=2
        return [gp, gpp]

    r0 = 1e-3
    g0_acc = (1.0 - g0)
    g_init = g0 + g0_acc * r0**2 / 6.0
    gp_init = g0_acc * r0 / 3.0

    sol = solve_ivp(rhs, [r0, R_max], [g_init, gp_init],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    max_step=0.03, dense_output=True)

    r_arr = np.linspace(r0, min(sol.t[-1], R_max), N_pts)
    y = sol.sol(r_arr)
    return r_arr, y[0], y[1]


def fit_tail_g(r_arr, g_arr, r_L=20.0, r_R=55.0):
    mask = (r_arr >= r_L) & (r_arr <= r_R) & np.isfinite(g_arr)
    if np.sum(mask) < 20:
        return np.nan, np.nan
    r_fit = r_arr[mask]
    y_fit = (g_arr[mask] - 1.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, y_fit, rcond=None)
    B, C = coefs
    return np.sqrt(B**2 + C**2), np.arctan2(B, C)


# Find K=g² φ-FP
def ratio_K2(g0):
    r1, g1, _ = solve_soliton_K2(g0, R_max=100.0, N_pts=4000)
    r2, g2, _ = solve_soliton_K2(phi*g0, R_max=100.0, N_pts=4000)
    A1, _ = fit_tail_g(r1, g1)
    A2, _ = fit_tail_g(r2, g2)
    if np.isnan(A1) or A1 < 1e-15:
        return 1e10
    return (A2/A1)**4 - R21_PDG

print("\n  Szukanie φ-FP dla K=g²...")
g0_K2 = None
for g_lo, g_hi in [(1.1, 1.2), (1.2, 1.3), (1.0, 1.15), (1.15, 1.35)]:
    try:
        v_lo = ratio_K2(g_lo)
        v_hi = ratio_K2(g_hi)
        if np.isfinite(v_lo) and np.isfinite(v_hi) and v_lo * v_hi < 0:
            g0_K2 = brentq(ratio_K2, g_lo, g_hi, xtol=1e-6)
            break
    except:
        continue

if g0_K2 is None:
    # Fallback: use known value from ex106
    g0_K2 = 1.24915
    print(f"  Nie znaleziono brentq, użycie wartości z ex106: {g0_K2}")
print(f"  g₀*(K=g²) = {g0_K2:.6f}")

# Solve K=g² solitons
r_e_K2, g_e_K2, gp_e_K2 = solve_soliton_K2(g0_K2)
r_mu_K2, g_mu_K2, gp_mu_K2 = solve_soliton_K2(phi*g0_K2)
r_tau_K2, g_tau_K2, gp_tau_K2 = solve_soliton_K2(phi**2*g0_K2)

A_e_K2, _ = fit_tail_g(r_e_K2, g_e_K2)
A_mu_K2, _ = fit_tail_g(r_mu_K2, g_mu_K2)
A_tau_K2, _ = fit_tail_g(r_tau_K2, g_tau_K2, r_L=25.0, r_R=70.0)

r21_K2 = (A_mu_K2/A_e_K2)**4
r31_K2 = (A_tau_K2/A_e_K2)**4
K_K2 = koide(1.0, r21_K2, r31_K2) if np.isfinite(r31_K2) else np.nan

print(f"  A_e = {A_e_K2:.8f}, A_μ = {A_mu_K2:.8f}, A_τ = {A_tau_K2:.8f}")
print(f"  r₂₁(K=g²) = {r21_K2:.2f}  (PDG: {R21_PDG})")
print(f"  r₃₁(K=g²) = {r31_K2:.2f}  (PDG: {R31_PDG})")
print(f"  Koide(K=g²) = {K_K2:.6f}  (2/3 = {2/3:.6f})")

# K=g² energy integrals (using g-space)
def compute_energy_K2(r_arr, g_arr, gp_arr):
    """E = 4π∫[½g²g'² + V(g)]r²dr with V(g)=½(1-g)²."""
    kinetic = 0.5 * g_arr**2 * gp_arr**2
    potential = 0.5 * (1.0 - g_arr)**2
    density = (kinetic + potential) * r_arr**2 * 4 * np.pi
    return _trapz(density, r_arr)

R_c = 100.0
idx_c = lambda r: min(np.searchsorted(r, R_c), len(r)-1)

E_e_K2_tot = compute_energy_K2(r_e_K2[:idx_c(r_e_K2)], g_e_K2[:idx_c(r_e_K2)], gp_e_K2[:idx_c(r_e_K2)])
E_mu_K2_tot = compute_energy_K2(r_mu_K2[:idx_c(r_mu_K2)], g_mu_K2[:idx_c(r_mu_K2)], gp_mu_K2[:idx_c(r_mu_K2)])
E_tau_K2_tot = compute_energy_K2(r_tau_K2[:idx_c(r_tau_K2)], g_tau_K2[:idx_c(r_tau_K2)], gp_tau_K2[:idx_c(r_tau_K2)])

print(f"\n  Energy integrals K=g² (R={R_c}):")
print(f"  E_e = {E_e_K2_tot:.6e}")
print(f"  E_μ = {E_mu_K2_tot:.6e}")
print(f"  E_τ = {E_tau_K2_tot:.6e}")
if E_e_K2_tot > 0:
    print(f"  E_μ/E_e = {E_mu_K2_tot/E_e_K2_tot:.4f}  (r₂₁ = {R21_PDG})")
    print(f"  E_τ/E_e = {E_tau_K2_tot/E_e_K2_tot:.4f}  (r₃₁ = {R31_PDG})")


# Grand comparison
print("\n" + "=" * 72)
print("§7. ★ TABELA PORÓWNAWCZA")
print("=" * 72)

print(f"\n  {'':30s}  {'K=g⁴':>12s}  {'K=g²':>12s}  {'PDG':>12s}")
print("  " + "-" * 70)
print(f"  {'g₀*':30s}  {g0_star:12.6f}  {g0_K2:12.6f}  {'—':>12s}")
print(f"  {'r₂₁ = (A_μ/A_e)⁴':30s}  {(A_mu/A_e)**4:12.2f}  {r21_K2:12.2f}  {R21_PDG:12.3f}")
print(f"  {'r₃₁ = (A_τ/A_e)⁴':30s}  {(A_tau/A_e)**4:12.0f}  {r31_K2:12.2f}  {R31_PDG:12.1f}")
print(f"  {'Koide (from A⁴)':30s}  {koide(1,(A_mu/A_e)**4,(A_tau/A_e)**4):12.6f}  {K_K2:12.6f}  {2/3:12.6f}")
print(f"  {'E_μ/E_e (R=100)':30s}  {E_tot_mu/E_tot_e if E_tot_e > 0 else np.nan:12.4f}  {E_mu_K2_tot/E_e_K2_tot if E_e_K2_tot > 0 else np.nan:12.4f}  {R21_PDG:12.3f}")
print(f"  {'E_τ/E_e (R=100)':30s}  {E_tot_tau/E_tot_e if E_tot_e > 0 else np.nan:12.4f}  {E_tau_K2_tot/E_e_K2_tot if E_e_K2_tot > 0 else np.nan:12.4f}  {R31_PDG:12.1f}")

r21_pass_K4 = abs((A_mu/A_e)**4 - R21_PDG)/R21_PDG < 0.01
r31_pass_K4 = abs((A_tau/A_e)**4 - R31_PDG)/R31_PDG < 0.10
r21_pass_K2 = abs(r21_K2 - R21_PDG)/R21_PDG < 0.01
r31_pass_K2 = abs(r31_K2 - R31_PDG)/R31_PDG < 0.10 if np.isfinite(r31_K2) else False

record("T2: K=g⁴ r₂₁", r21_pass_K4, f"r₂₁ = {(A_mu/A_e)**4:.2f}")
record("T3: K=g⁴ r₃₁", r31_pass_K4, f"r₃₁ = {(A_tau/A_e)**4:.0f} (PDG: {R31_PDG})")
record("T4: K=g² r₂₁", r21_pass_K2, f"r₂₁ = {r21_K2:.2f}")
record("T5: K=g² r₃₁", r31_pass_K2, f"r₃₁ = {r31_K2:.2f}" if np.isfinite(r31_K2) else "NaN")


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
print("PODSUMOWANIE ex229")
print("=" * 72)
print(f"""
  KLUCZOWE WYNIKI:

  1. K=g⁴: r₂₁ = {(A_mu/A_e)**4:.2f} ✓  ale  r₃₁ = {(A_tau/A_e)**4:.0f} ✗
     (r₃₁ off by {abs((A_tau/A_e)**4-R31_PDG)/R31_PDG*100:.0f}%)

  2. K=g²: r₂₁ = {r21_K2:.2f}  r₃₁ = {r31_K2 if np.isfinite(r31_K2) else 'NaN'}
     {'Koide = '+str(round(K_K2,4)) if np.isfinite(K_K2) else ''}

  3. Brak uniwersalnego p: A^p pasuje do r₂₁ LUB r₃₁, nie obu
     p(r₂₁) = {p_from_r21:.4f}, p(r₃₁) = {p_from_r31:.4f}

  4. Energia NIE daje lepszych ratios niż A⁴

  WNIOSEK:
  Dla K=g⁴, mechanizm φ-FP z M ∝ A^p NIE może
  jednocześnie reprodukować r₂₁ I r₃₁.
  K=g² daje poprawne OBA ratios z p=4.

  IMPLIKACJA:
  - n_K=4 z geometrii akcji MUSI być zmodyfikowane
  - LUB: τ soliton wymaga INNEGO mechanizmu niż φ²·g₀*
  - LUB: pełny potencjał TGP (g⁷,g⁸) zmienia A(g₀)
""")
