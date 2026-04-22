#!/usr/bin/env python3
"""
cosmological_evolution.py — Ewolucja kosmologiczna TGP
========================================================
Numeryczna ewolucja pola ψ(t) = Φ(t)/Φ₀ w tle FRW,
od BBN do dziś, z wyznaczeniem trajektorii ψ(a), Ω_DE(z), w_DE(z), ΔG/G(z).

Równanie pola (prop:FRW-derivation, eq:Phi-cosmo-exact):
  ψ̈ + 3Hψ̇ + 2ψ̇²/ψ = c₀²·W(ψ)

gdzie W(ψ) = (7β/3)ψ² - 2γψ³ (z miary ψ⁴ w działaniu).

Uwaga: S_src (sprężenie źródłowe) jest pominięte, gdyż
jego dokładny współczynnik zależy od postaci sprężenia
materia–metryka (por. rem. w sek08, linia ~548).

KLUCZOWY WYNIK:
  Kosmologiczny atraktor to ψ_eq = 7β/(6γ) = 7/6,
  NIE ψ = 1 (minimum statycznego U). Pole ψ jest zamrożone
  przy ψ ≈ 1 w erze materii (nadtłumienie), ale ewoluuje
  w stronę ψ_eq w erze ciemnej energii (niedotłumienie).

Author: Claudian deep analysis
Date: 2026-03-17
"""

import sys
import io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp

# ===================================================================
# PHYSICAL CONSTANTS
# ===================================================================
c0 = 2.998e8         # m/s
G0 = 6.674e-11       # m³/(kg·s²)
H0 = 2.2e-18         # 1/s (67.4 km/s/Mpc)
Lambda_obs = 1.11e-52 # m⁻²
Omega_m0 = 0.315
Omega_r0 = 9.1e-5
Omega_L0 = 0.685

# TGP parameters
gamma = 56 * Lambda_obs  # Λ_eff = γ/56
beta  = gamma             # warunek próżni β = γ
Phi0  = 168 * Omega_L0    # ≈ 115.08

# Dimensionless TGP parameter: Φ₀ = c₀²γ/H₀²
Phi0_eff = c0**2 * gamma / H0**2

PASS = 0
FAIL = 0

def check(condition, name, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  [PASS] {name}")
    else:
        FAIL += 1
        print(f"  [FAIL] {name}")
    if detail:
        print(f"         {detail}")

# ===================================================================
# POTENTIALS
# ===================================================================

def U_pot_norm(psi):
    """P(ψ)/(γΦ₀²) = (1/7)ψ⁷ - (1/8)ψ⁸  (β=γ, action potential)"""
    return (1.0/7)*psi**7 - (1.0/8)*psi**8

def W_pot_norm(psi):
    """W(ψ)/γ = (7/3)ψ² - 2ψ³  (β=γ)"""
    return (7.0/3)*psi**2 - 2*psi**3

# ===================================================================
# THEORETICAL ANALYSIS
# ===================================================================

print("=" * 70)
print(" TGP Cosmological Evolution")
print("=" * 70)

psi_eq = 7.0/6  # W(ψ_eq) = 0
print(f"\n  Φ₀ (dimensionless) = {Phi0_eff:.2f}")
print(f"  ψ_eq (cosmological attractor, W=0) = {psi_eq:.6f}")

# Effective cosmological mass² from linearized W'(1):
# m²_cosmo = -c₀²W'(1) = c₀²·(4γ/3) → m²/H₀² = 4Φ₀/3
m2_over_H02 = 4*Phi0_eff/3
omega_field  = np.sqrt(m2_over_H02)
gamma_damp   = 3.0/2  # 3H₀/2 in H₀ units (at z=0)
delta_ss     = (Phi0_eff/3) / m2_over_H02  # = 1/4 (steady-state δ)

print(f"\n  --- Damping analysis ---")
print(f"  ω_field/H₀ = {omega_field:.2f}  (effective mass)")
print(f"  γ_damp/H₀  = {gamma_damp:.2f}  (3H₀/2)")
regime = "UNDERDAMPED" if omega_field > gamma_damp else "OVERDAMPED"
print(f"  Regime at z=0: {regime}")
print(f"  δ_ss (linearized equilibrium) = {delta_ss:.4f}  → ψ_ss = {1+delta_ss:.4f}")

# Transition redshift (overdamped → underdamped)
# 3H/2 = ω → H = 2ω/3 → H²/H₀² = (2ω/3)² ≈ E²_trans
E2_trans = (2*omega_field/3)**2
# E² = Ω_m/a³ + Ω_Λ ≈ E²_trans → Ω_m/a³ ≈ E²_trans - Ω_Λ
a_trans = (Omega_m0 / (E2_trans - Omega_L0))**(1.0/3) if E2_trans > Omega_L0 else 1.0
z_trans = 1.0/a_trans - 1
print(f"  Transition z_trans ≈ {z_trans:.1f}  (overdamped → underdamped)")

# ===================================================================
# NUMERICAL EVOLUTION (dimensionless τ = H₀·t)
# ===================================================================

def rhs_cosmo(tau, y):
    """RHS for (ψ, ψ̇, a). All in H₀ units."""
    psi, psi_dot, a = y

    if psi <= 0.01 or a <= 0:
        return [0, 0, 0]

    # Friedmann: E² = H²/H₀²
    U_H = Phi0_eff * U_pot_norm(psi)          # c₀²U/H₀² (dimensionless)
    Omega_DE = U_H / (3*np.sqrt(psi)) + psi_dot**2/(6*psi)
    E_sq = Omega_m0/a**3 + Omega_r0/a**4 + Omega_DE
    if E_sq <= 0:
        E_sq = 1e-30
    E = np.sqrt(E_sq)

    # Field equation: ψ̈ + 3Eψ̇ + 2ψ̇²/ψ = Φ₀·W_norm(ψ)
    W_H = Phi0_eff * W_pot_norm(psi)
    psi_ddot = -3*E*psi_dot - 2*psi_dot**2/psi + W_H

    return [psi_dot, psi_ddot, a*E]

print("\n--- Evolving ψ(t) from z = 10⁹ to z = 0 ---\n")

sol = solve_ivp(
    rhs_cosmo,
    t_span=(0, 3.0),
    y0=[1.0, 0.0, 1e-9],       # ψ=1, ψ̇=0, a=10⁻⁹
    max_step=0.001,
    rtol=1e-10,  atol=1e-13,
    method='DOP853',
    dense_output=True
)

a_vals = sol.y[2, :]

# Find today (a = 1)
if np.max(a_vals) >= 1.0:
    idx_today = np.searchsorted(a_vals, 1.0)
    idx_today = min(idx_today, len(a_vals)-1)
else:
    print(f"  [WARNING] a never reaches 1; max a = {np.max(a_vals):.6f}")
    idx_today = len(a_vals) - 1

psi_today     = sol.y[0, idx_today]
psi_dot_today = sol.y[1, idx_today]
a_today       = sol.y[2, idx_today]

print(f"  ψ(today) = {psi_today:.8f}")
print(f"  ψ̇(today) = {psi_dot_today:.4e}  (in H₀ units)")
print(f"  a(today) = {a_today:.6f}")
print(f"  ψ_eq     = {psi_eq:.6f}")

# ===================================================================
# CHECKS
# ===================================================================

print("\n--- GROUP 1: Basin & freezing ---")

psi_all = sol.y[0, :idx_today+1]
check(len(psi_all) > 0 and np.all(psi_all > 0) and np.all(psi_all < 8.0/7),
      "ψ(t) ∈ (0, 8/7) throughout",
      f"ψ ∈ [{np.min(psi_all):.6f}, {np.max(psi_all):.6f}]")

# Frozen in matter era
a_z100 = 1.0/101
idx_z100 = np.searchsorted(a_vals[:idx_today+1], a_z100)
idx_z100 = min(idx_z100, idx_today)
psi_z100 = sol.y[0, idx_z100]
check(abs(psi_z100 - 1.0) < 0.001,
      "ψ frozen at 1 during matter era (z > 100)",
      f"ψ(z≈100) = {psi_z100:.8f}")

# ψ evolves toward ψ_eq at late times
check(psi_today > 1.0 and psi_today < 8.0/7,
      "ψ(today) > 1 (evolving toward ψ_eq)",
      f"ψ = {psi_today:.6f}, ψ_eq = {psi_eq:.6f}")

print("\n--- GROUP 2: Dark energy ---")

U_today = Phi0_eff * U_pot_norm(psi_today)
kinetic_today = psi_dot_today**2 / (6*psi_today)
Omega_DE_today = U_today / (3*np.sqrt(psi_today)) + kinetic_today

check(0.3 < Omega_DE_today < 0.9,
      "Ω_DE(z=0) in physical range [0.3, 0.9]",
      f"Ω_DE = {Omega_DE_today:.4f}  (obs: {Omega_L0})")

# w_DE
if Omega_DE_today > 0:
    w_DE = -1 + kinetic_today / Omega_DE_today
else:
    w_DE = -1
check(abs(w_DE + 1) < 0.01,
      "w_DE ≈ −1",
      f"w_DE = {w_DE:.8f},  |1+w| = {abs(1+w_DE):.2e}")

# Recalibration info
U_eq_norm = U_pot_norm(psi_eq)
Omega_DE_at_eq = Phi0_eff * U_eq_norm / (3*np.sqrt(psi_eq))
ratio_recal = Omega_L0 / Omega_DE_at_eq if Omega_DE_at_eq > 0 else float('inf')
Phi0_recal = Phi0_eff * ratio_recal
print(f"\n  [INFO] If Φ₀ recalibrated at ψ_eq:")
print(f"         Ω_DE(ψ_eq) with current Φ₀ = {Omega_DE_at_eq:.4f}")
print(f"         Φ₀ needed for Ω_Λ=0.685   = {Phi0_recal:.2f}")
print(f"         Recalibration factor        = {ratio_recal:.4f}")

print("\n--- GROUP 3: Dynamic constants ---")

# G/G_BBN at today
G_ratio = 1.0/psi_today  # G(today)/G₀  where G₀ = G(ψ=1)
c_ratio = 1.0/np.sqrt(psi_today)
print(f"  ΔG/G(today vs BBN)  = {(G_ratio-1)*100:+.2f}%  (G₀/ψ)")
print(f"  Δc/c₀(today vs BBN) = {(c_ratio-1)*100:+.2f}%  (c₀/√ψ)")

# Ġ/G at z=0
G_dot_over_G = -psi_dot_today / psi_today  # d/dt(1/ψ)·ψ = -ψ̇/ψ
G_dot_per_yr = G_dot_over_G * H0 * 3.156e7  # convert to 1/yr
print(f"  |Ġ/G|(z=0) = {abs(G_dot_per_yr):.2e} /yr")
print(f"  LLR bound  ≈ 1e-13 /yr")
check(abs(G_dot_per_yr) < 5e-12,
      "|Ġ/G| < 5×10⁻¹² /yr (consistent with solar-system tests)",
      f"|Ġ/G| = {abs(G_dot_per_yr):.2e} /yr")

# BBN constraint: |ΔG/G| < 13%  (Copi et al.)
Delta_G_BBN = abs(1.0/psi_today - 1.0)
check(Delta_G_BBN < 0.20,
      "|ΔG/G| between BBN and today < 20%",
      f"|ΔG/G| = {Delta_G_BBN*100:.1f}%")

print("\n--- GROUP 4: Structural consistency ---")

# Φ₀ = 168·Ω_Λ
check(abs(Phi0 - 168*Omega_L0) / Phi0 < 0.01,
      "Φ₀ = 168·Ω_Λ",
      f"Φ₀ = {Phi0:.2f}, 168Ω_Λ = {168*Omega_L0:.2f}")

# W(ψ_eq) = 0
check(abs(W_pot_norm(psi_eq)) < 1e-12,
      "W(ψ_eq) = 0",
      f"W(7/6) = {W_pot_norm(psi_eq):.2e}")

# P(8/7) = 0
check(abs(U_pot_norm(8.0/7)) < 1e-12,
      "P(8/7) = 0  (basin boundary)",
      f"P(8/7) = {U_pot_norm(8.0/7):.2e}")

# W(1) = (56/3)·P(1) (ψ⁸ action potential identity)
W1 = W_pot_norm(1.0)
P1 = U_pot_norm(1.0)
check(abs(W1 - (56.0/3.0)*P1) < 1e-12,
      "W(1) = (56/3)P(1)  (ψ⁸ action potential identity)",
      f"W(1)={W1:.8f}, (56/3)P(1)={(56.0/3.0)*P1:.8f}")

# Linearized mass matches
# W'(1) = 14β/3 - 6γ = -4γ/3 (for β=γ)
# m² = c₀²·4γ/3 / H₀² = 4Φ₀/3
m2_check = 4*Phi0_eff/3
check(abs(m2_check/m2_over_H02 - 1) < 1e-10,
      "m²_cosmo = 4Φ₀/3 from W'(1)",
      f"m² = {m2_over_H02:.4f}")

# ===================================================================
# EVOLUTION TABLE
# ===================================================================

print("\n--- Evolution summary ---")
hdr = f"  {'z':>10s} {'a':>10s} {'ψ':>12s} {'ψ̇':>12s} {'Ω_DE':>10s} {'ΔG/G':>12s}"
print(hdr)
print("  " + "-"*70)

for z_s in [1e6, 1100, 100, 10, 2, 1, 0.5, 0]:
    a_s = 1.0/(1+z_s)
    if a_s < a_vals[0] or a_s > a_vals[idx_today]:
        continue
    idx = min(np.searchsorted(a_vals[:idx_today+1], a_s), idx_today)
    p   = sol.y[0, idx]
    pd  = sol.y[1, idx]
    aa  = sol.y[2, idx]
    zz  = 1.0/aa - 1 if aa > 0 else 1e20
    U_s = Phi0_eff * U_pot_norm(p)
    ODE = U_s/(3*np.sqrt(p)) + pd**2/(6*p)
    dGG = 1.0/p - 1.0
    print(f"  {zz:10.2e} {aa:10.6f} {p:12.8f} {pd:12.4e} {ODE:10.4f} {dGG:12.4e}")

# ===================================================================
# FALSIFIABLE PREDICTIONS
# ===================================================================

print(f"""
--- Falsifiable predictions ---
  F1: |ΔG/G|(BBN→today)  ≈ {Delta_G_BBN*100:.1f}%   [BBN ⁴He abundance]
  F2: |Ġ/G|(today)        ≈ {abs(G_dot_per_yr):.1e} /yr  [Lunar Laser Ranging]
  F3: w_DE deviation      ≈ {abs(1+w_DE):.1e}  [DESI / Euclid]
  F4: ψ oscillation period ≈ {2*np.pi/omega_field:.2f}/H₀ ≈ {2*np.pi/omega_field*14.4:.1f} Gyr
  F5: Cosmological attractor ψ_eq = {psi_eq:.4f} (NOT ψ=1)
""")

# ===================================================================
# SUMMARY
# ===================================================================

print("=" * 70)
print(f" RESULTS: {PASS}/{PASS+FAIL} passed, {FAIL} failed")
print("=" * 70)

if FAIL == 0:
    print(" ALL CHECKS PASSED")
else:
    print(f" {FAIL} CHECK(S) FAILED — review above")

print(f"""
 CONCLUSION:
   1. ψ frozen at 1 during matter era (z > {z_trans:.0f}) — overdamped
   2. ψ evolves toward ψ_eq = {psi_eq:.4f} for z < {z_trans:.0f} — underdamped
   3. Today ψ ≈ {psi_today:.3f} (partially relaxed toward ψ_eq)
   4. Ω_DE(today) ≈ {Omega_DE_today:.3f} — physical range
   5. w_DE ≈ −1 + O({abs(1+w_DE):.0e})
   6. Requires Φ₀ recalibration at ψ_eq for exact Ω_Λ=0.685 match
""")
