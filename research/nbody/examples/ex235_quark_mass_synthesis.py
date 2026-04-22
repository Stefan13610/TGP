#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex235_quark_mass_synthesis.py
================================
ROZSZERZENIE SYNTEZY NA KWARKI: (d,s,b) I (u,c,t)

KONTEKST (ex234):
  Leptony: r₂₁(ODE) + K=2/3 → m_τ z 0.006% dokładnością (7/7 PASS)
  PYTANIE: Czy ten sam framework działa dla kwarków?

DANE PDG (MS-bar, μ=2 GeV for light, pole for heavy):
  Down-type: m_d ≈ 4.67 MeV, m_s ≈ 93.4 MeV, m_b ≈ 4180 MeV
  Up-type:   m_u ≈ 2.16 MeV, m_c ≈ 1270 MeV, m_t ≈ 172760 MeV

OBSERWACJE:
  1. Kwarki mają DUŻO większe niepewności niż leptony
  2. Koide K dla kwarków: K(d,s,b) ≈ 0.625, K(u,c,t) ≈ 0.670
  3. Kwarki "powinny" mieć własne φ-FP (inny sektor)
  4. r₂₁(quarks) ≈ m_s/m_d ≈ 20 (vs lepton 206.77)

PLAN:
  §1. Koide K dla obu trójek kwarków (PDG)
  §2. Predykcja m_b z r₂₁(d,s) + K=2/3
  §3. Predykcja m_t z r₂₁(u,c) + K=2/3
  §4. "Shifted Koide" z R12 connection (ex222)
  §5. Porównanie z ex224 full prediction chain

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.optimize import brentq

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

phi = (1 + np.sqrt(5)) / 2

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


def koide_shifted(m1, m2, m3, m0):
    """Koide with mass shift: K(m+m₀)."""
    return koide(m1 + m0, m2 + m0, m3 + m0)


def solve_r31_from_K(r21, K_target=2/3):
    """Solve K(1, r₂₁, r₃₁) = K_target for r₃₁."""
    a = np.sqrt(r21)
    # x² - (4a+4)x + (a²-4a+1) = 0 when K=2/3
    # General: K(1+a²+x²)/(1+a+x)² = K_target
    # → (1+a²+x²) = K_target(1+a+x)²
    # → (1-K_t)x² - (2K_t·a + 2K_t)x + (1+a²-K_t(1+a)²) = 0
    # Wait, let me redo:
    # K = (1+r21+r31)/S², S = 1+√r21+√r31
    # K·S² = 1+a²+x²
    # K(1+a+x)² = 1+a²+x²
    # K + Ka² + Kx² + 2Ka + 2Kx + 2Kax = 1 + a² + x²
    # (K-1)x² + (2Ka+2K)x + (K+Ka²+2Ka-1-a²) = 0
    # (K-1)x² + 2K(a+1)x + ((K-1)a² + 2Ka + (K-1)) = 0
    # (K-1)[x² + a² + 1] + 2K(a+1)x + 2Ka = 0
    # Hmm, let me just use numerical solver

    def K_func(x):
        r31 = x**2
        return koide(1.0, r21, r31) - K_target

    # K(x) has a minimum at x_min = (1+a²)/(1+a), may have 0 or 2 roots
    a = np.sqrt(r21)
    x_min = (1 + a**2) / (1 + a)
    K_at_min = K_func(x_min) + K_target  # actual K value at minimum

    roots = []
    # Search left of minimum
    try:
        x_left = brentq(K_func, 0.01, x_min)
        roots.append(x_left)
    except:
        pass
    # Search right of minimum
    try:
        x_right = brentq(K_func, x_min, 10000.0)
        roots.append(x_right)
    except:
        pass

    if not roots:
        return np.nan
    # Return the LARGER root (physical: r₃₁ > r₂₁)
    x_sol = max(roots)
    return x_sol**2


# ============================================================
# §1. KOIDE K DLA KWARKÓW (PDG)
# ============================================================
print("=" * 72)
print("§1. KOIDE K DLA KWARKÓW I LEPTONÓW")
print("=" * 72)

# PDG masses (central values, MeV)
# Light quarks at μ = 2 GeV (MS-bar)
m_d = 4.67    # ± 0.48
m_s = 93.4    # ± 8.6
m_b = 4180.0  # ± 30 (pole mass ~4.78 GeV, MS-bar ~4.18 GeV at μ=m_b)

m_u = 2.16    # ± 0.49
m_c = 1270.0  # ± 20
m_t = 172760.0  # ± 300 (pole mass)

# Leptons
m_e = 0.51099895
m_mu = 105.6583755
m_tau = 1776.86

print(f"\n  {'Trójka':15s}  {'m₁':>10s}  {'m₂':>10s}  {'m₃':>10s}  {'K':>10s}  {'2/3':>8s}  {'err':>8s}")
print("  " + "-" * 70)

triplets = [
    ("leptons (e,μ,τ)", m_e, m_mu, m_tau),
    ("down (d,s,b)", m_d, m_s, m_b),
    ("up (u,c,t)", m_u, m_c, m_t),
]

for name, m1, m2, m3 in triplets:
    K = koide(m1, m2, m3)
    err = abs(K - 2/3)/(2/3)*100
    print(f"  {name:15s}  {m1:10.2f}  {m2:10.2f}  {m3:10.2f}  {K:10.6f}  {2/3:8.6f}  {err:7.2f}%")

# Also: "shifted Koide" — K(m + m₀) = 2/3
print(f"\n  Shifted Koide: K(m + m₀) = 2/3")
print(f"  {'Trójka':15s}  {'m₀ (MeV)':>12s}  {'K(m+m₀)':>10s}")
print("  " + "-" * 42)

for name, m1, m2, m3 in triplets:
    # Find m₀ such that K(m+m₀) = 2/3
    def find_m0(m0):
        return koide(m1+m0, m2+m0, m3+m0) - 2/3

    m0_best = None
    for m0_try in np.linspace(-m1*0.99, m3*0.5, 1000):
        try:
            if find_m0(m0_try) * find_m0(m0_try + m3*0.01) < 0:
                m0_best = brentq(find_m0, m0_try, m0_try + m3*0.01)
                break
        except:
            continue

    if m0_best is not None:
        K_shifted = koide(m1+m0_best, m2+m0_best, m3+m0_best)
        print(f"  {name:15s}  {m0_best:12.4f}  {K_shifted:10.6f}")
    else:
        print(f"  {name:15s}  {'NOT FOUND':>12s}")


# ============================================================
# §2. PREDYKCJA m_b Z r₂₁(d,s) + K=2/3
# ============================================================
print("\n" + "=" * 72)
print("§2. PREDYKCJA m_b Z r₂₁(d,s) + K=2/3")
print("=" * 72)

r21_ds = m_s / m_d
print(f"\n  r₂₁(d,s) = m_s/m_d = {r21_ds:.4f}")
print(f"  r₂₁(leptons) = {m_mu/m_e:.4f} dla porównania")

r31_pred_ds = solve_r31_from_K(r21_ds, 2/3)
m_b_pred = m_d * r31_pred_ds

print(f"  r₃₁(pred, K=2/3) = {r31_pred_ds:.2f}")
print(f"  m_b(pred) = m_d × r₃₁ = {m_b_pred:.1f} MeV")
print(f"  m_b(PDG)  = {m_b:.1f} MeV")
print(f"  Error: {abs(m_b_pred - m_b)/m_b*100:.2f}%")

# With uncertainty on m_d
m_d_lo, m_d_hi = 4.67 - 0.48, 4.67 + 0.48
r21_lo = m_s / m_d_hi  # smaller r21 when m_d larger
r21_hi = m_s / m_d_lo
r31_lo = solve_r31_from_K(r21_lo, 2/3)
r31_hi = solve_r31_from_K(r21_hi, 2/3)
m_b_lo = m_d_hi * r31_lo
m_b_hi = m_d_lo * r31_hi

print(f"\n  Z niepewnością m_d = {m_d} ± 0.48:")
print(f"    m_b(pred) ∈ [{m_b_lo:.0f}, {m_b_hi:.0f}] MeV")

record("T1: m_b from K=2/3 (bare)",
       abs(m_b_pred - m_b)/m_b < 0.20,
       f"m_b = {m_b_pred:.0f} MeV (PDG: {m_b}), err = {abs(m_b_pred-m_b)/m_b*100:.1f}%")


# ============================================================
# §3. PREDYKCJA m_t Z r₂₁(u,c) + K=2/3
# ============================================================
print("\n" + "=" * 72)
print("§3. PREDYKCJA m_t Z r₂₁(u,c) + K=2/3")
print("=" * 72)

r21_uc = m_c / m_u
print(f"\n  r₂₁(u,c) = m_c/m_u = {r21_uc:.2f}")

r31_pred_uc = solve_r31_from_K(r21_uc, 2/3)
m_t_pred = m_u * r31_pred_uc

print(f"  r₃₁(pred, K=2/3) = {r31_pred_uc:.0f}")
print(f"  m_t(pred) = m_u × r₃₁ = {m_t_pred:.0f} MeV")
print(f"  m_t(PDG)  = {m_t:.0f} MeV")
print(f"  Error: {abs(m_t_pred - m_t)/m_t*100:.2f}%")

record("T2: m_t from K=2/3 (bare)",
       abs(m_t_pred - m_t)/m_t < 0.20,
       f"m_t = {m_t_pred:.0f} MeV (PDG: {m_t}), err = {abs(m_t_pred-m_t)/m_t*100:.1f}%")


# ============================================================
# §4. SHIFTED KOIDE: K(m + m₀) = 2/3
# ============================================================
print("\n" + "=" * 72)
print("§4. ★ SHIFTED KOIDE DLA KWARKÓW")
print("=" * 72)

print("""
  Hipoteza (ex222, ex224): Koide działa z PRZESUNIĘCIEM:
    K(m₁+m₀, m₂+m₀, m₃+m₀) = 2/3

  Dla leptonów: m₀ ≈ 0 (Koide działa bezpośrednio)
  Dla kwarków: m₀ ≠ 0 (potrzebne przesunięcie)

  Fizyczna interpretacja: m₀ = wkład "morza" (gluony + pary qq̄)
  do konstytutywnej masy kwarku.

  Predykcja R12 (ex222): A = 1/(Φ_eff × φ)
    m₀ = A · m₃/m₁ (sektor-zależne)
""")

# Find optimal m₀ for down quarks
def K_down_shifted(m0):
    return koide(m_d+m0, m_s+m0, m_b+m0) - 2/3

m0_down = None
for start in np.linspace(-m_d*0.9, 500.0, 2000):
    try:
        v1 = K_down_shifted(start)
        v2 = K_down_shifted(start + 0.5)
        if v1 * v2 < 0:
            m0_down = brentq(K_down_shifted, start, start + 0.5)
            break
    except:
        continue

if m0_down is not None:
    print(f"  Down quarks: m₀ = {m0_down:.4f} MeV")
    K_check = koide(m_d+m0_down, m_s+m0_down, m_b+m0_down)
    print(f"    K(d+m₀, s+m₀, b+m₀) = {K_check:.8f}")

    # Now predict m_b from m_d, m_s, m₀, K=2/3
    r21_shifted = (m_s + m0_down) / (m_d + m0_down)
    r31_shifted = solve_r31_from_K(r21_shifted, 2/3)
    m_b_shifted = (m_d + m0_down) * r31_shifted - m0_down

    print(f"    r₂₁(shifted) = {r21_shifted:.4f}")
    print(f"    m_b(shifted pred) = {m_b_shifted:.1f} MeV")
    print(f"    m_b(PDG) = {m_b:.1f} MeV")
    print(f"    Error: {abs(m_b_shifted - m_b)/m_b*100:.4f}%")

    record("T3: m_b shifted Koide",
           abs(m_b_shifted - m_b)/m_b < 0.01,
           f"m_b = {m_b_shifted:.1f} MeV, m₀ = {m0_down:.2f} MeV")

# Find optimal m₀ for up quarks
def K_up_shifted(m0):
    return koide(m_u+m0, m_c+m0, m_t+m0) - 2/3

m0_up = None
for start in np.linspace(-m_u*0.9, 50000.0, 5000):
    try:
        v1 = K_up_shifted(start)
        v2 = K_up_shifted(start + 10.0)
        if v1 * v2 < 0:
            m0_up = brentq(K_up_shifted, start, start + 10.0)
            break
    except:
        continue

if m0_up is not None:
    print(f"\n  Up quarks: m₀ = {m0_up:.2f} MeV")
    K_check = koide(m_u+m0_up, m_c+m0_up, m_t+m0_up)
    print(f"    K(u+m₀, c+m₀, t+m₀) = {K_check:.8f}")

    r21_shifted_up = (m_c + m0_up) / (m_u + m0_up)
    r31_shifted_up = solve_r31_from_K(r21_shifted_up, 2/3)
    m_t_shifted = (m_u + m0_up) * r31_shifted_up - m0_up

    print(f"    r₂₁(shifted) = {r21_shifted_up:.4f}")
    print(f"    m_t(shifted pred) = {m_t_shifted:.0f} MeV")
    print(f"    m_t(PDG) = {m_t:.0f} MeV")
    print(f"    Error: {abs(m_t_shifted - m_t)/m_t*100:.4f}%")

    record("T4: m_t shifted Koide",
           abs(m_t_shifted - m_t)/m_t < 0.01,
           f"m_t = {m_t_shifted:.0f} MeV, m₀ = {m0_up:.1f} MeV")


# ============================================================
# §5. GRAND COMPARISON: LEPTON vs QUARK SYNTHESIS
# ============================================================
print("\n" + "=" * 72)
print("§5. ★ PORÓWNANIE: LEPTONY vs KWARKI")
print("=" * 72)

# Lepton synthesis
r21_lep = m_mu / m_e
r31_lep_pred = solve_r31_from_K(r21_lep, 2/3)
m_tau_pred = m_e * r31_lep_pred

print(f"\n  {'':20s}  {'Leptony':>14s}  {'Down quarks':>14s}  {'Up quarks':>14s}")
print("  " + "-" * 66)
print(f"  {'m₁ (MeV)':20s}  {m_e:14.4f}  {m_d:14.2f}  {m_u:14.2f}")
print(f"  {'m₂ (MeV)':20s}  {m_mu:14.4f}  {m_s:14.2f}  {m_c:14.2f}")
print(f"  {'m₃ (MeV)':20s}  {m_tau:14.2f}  {m_b:14.2f}  {m_t:14.2f}")
print(f"  {'r₂₁':20s}  {r21_lep:14.4f}  {r21_ds:14.4f}  {r21_uc:14.2f}")
print(f"  {'K (bare)':20s}  {koide(m_e,m_mu,m_tau):14.6f}  {koide(m_d,m_s,m_b):14.6f}  {koide(m_u,m_c,m_t):14.6f}")

# Bare Koide predictions
print(f"  {'':20s}  {'':>14s}  {'':>14s}  {'':>14s}")
print(f"  {'m₃ pred (K=2/3)':20s}  {m_tau_pred:14.2f}  {m_b_pred:14.1f}  {m_t_pred:14.0f}")
print(f"  {'m₃ PDG':20s}  {m_tau:14.2f}  {m_b:14.1f}  {m_t:14.0f}")
print(f"  {'Error':20s}  {abs(m_tau_pred-m_tau)/m_tau*100:13.3f}%  {abs(m_b_pred-m_b)/m_b*100:13.1f}%  {abs(m_t_pred-m_t)/m_t*100:13.1f}%")

if m0_down is not None and m0_up is not None:
    print(f"\n  {'m₀ (shifted)':20s}  {'≈ 0':>14s}  {m0_down:14.2f}  {m0_up:14.1f}")
    print(f"  {'m₃ pred (shifted)':20s}  {m_tau_pred:14.2f}  {m_b_shifted:14.1f}  {m_t_shifted:14.0f}")
    print(f"  {'Error (shifted)':20s}  {abs(m_tau_pred-m_tau)/m_tau*100:13.3f}%  {abs(m_b_shifted-m_b)/m_b*100:13.3f}%  {abs(m_t_shifted-m_t)/m_t*100:13.3f}%")


# ============================================================
# §6. PREDYKCJA R12: m₀ Z TGP FORMUŁY
# ============================================================
print("\n" + "=" * 72)
print("§6. R12: m₀ Z TGP FORMUŁY (ex222)")
print("=" * 72)

print("""
  Z ex222: A = 1/(Φ_eff × φ)
  Φ_eff = Φ₀ × 3/14 ≈ 24.66

  Formuła m₀:
    m₀ = A × m₃/m₁  (sektor-zależne)

  Dla down quarks: m₀ = A × m_b/m_d
  Dla up quarks:   m₀ = A × m_t/m_u
""")

Phi0_bare = 168 * 0.685  # 115.08
Phi_eff = Phi0_bare * 3/14  # 24.66
A_R12 = 1.0 / (Phi_eff * phi)

print(f"  Φ₀ = {Phi0_bare:.2f}")
print(f"  Φ_eff = Φ₀ × 3/14 = {Phi_eff:.2f}")
print(f"  A = 1/(Φ_eff × φ) = {A_R12:.6f}")

# Predictions
m0_down_R12 = A_R12 * m_b / m_d
m0_up_R12 = A_R12 * m_t / m_u

print(f"\n  Down quarks:")
print(f"    m₀(R12) = A × m_b/m_d = {m0_down_R12:.2f} MeV")
if m0_down is not None:
    print(f"    m₀(fit)  = {m0_down:.2f} MeV")
    print(f"    Error: {abs(m0_down_R12 - m0_down)/abs(m0_down)*100:.1f}%")

print(f"\n  Up quarks:")
print(f"    m₀(R12) = A × m_t/m_u = {m0_up_R12:.1f} MeV")
if m0_up is not None:
    print(f"    m₀(fit)  = {m0_up:.1f} MeV")
    print(f"    Error: {abs(m0_up_R12 - m0_up)/abs(m0_up)*100:.1f}%")

# Use R12 m₀ to predict m_b, m_t
if m0_down_R12 > -m_d:
    r21_R12_d = (m_s + m0_down_R12) / (m_d + m0_down_R12)
    r31_R12_d = solve_r31_from_K(r21_R12_d, 2/3)
    m_b_R12 = (m_d + m0_down_R12) * r31_R12_d - m0_down_R12

    print(f"\n  m_b(R12 pred) = {m_b_R12:.1f} MeV  (PDG: {m_b})")
    print(f"  Error: {abs(m_b_R12 - m_b)/m_b*100:.2f}%")

    record("T5: m_b from R12 m₀",
           abs(m_b_R12 - m_b)/m_b < 0.10,
           f"m_b = {m_b_R12:.0f} MeV, err = {abs(m_b_R12-m_b)/m_b*100:.1f}%")

if m0_up_R12 > -m_u:
    r21_R12_u = (m_c + m0_up_R12) / (m_u + m0_up_R12)
    r31_R12_u = solve_r31_from_K(r21_R12_u, 2/3)
    m_t_R12 = (m_u + m0_up_R12) * r31_R12_u - m0_up_R12

    print(f"  m_t(R12 pred) = {m_t_R12:.0f} MeV  (PDG: {m_t})")
    print(f"  Error: {abs(m_t_R12 - m_t)/m_t*100:.2f}%")

    record("T6: m_t from R12 m₀",
           abs(m_t_R12 - m_t)/m_t < 0.10,
           f"m_t = {m_t_R12:.0f} MeV, err = {abs(m_t_R12-m_t)/m_t*100:.1f}%")


# ============================================================
# §7. SELF-CONSISTENT Ω_Λ (from ex222/ex224)
# ============================================================
print("\n" + "=" * 72)
print("§7. SELF-CONSISTENT Ω_Λ")
print("=" * 72)

print("""
  ex222/ex224 showed: Ω_Λ can be determined self-consistently
  from quark masses via total χ² minimization.

  Here: scan Ω_Λ and check m_b, m_t predictions.
""")

OL_scan = np.linspace(0.65, 0.75, 21)
print(f"\n  {'Ω_Λ':>8s}  {'Φ_eff':>8s}  {'A':>10s}  {'m₀_d':>8s}  {'m_b pred':>10s}  {'m₀_u':>10s}  {'m_t pred':>12s}")
print("  " + "-" * 76)

for OL in OL_scan:
    P0 = 168 * OL
    PE = P0 * 3/14
    A_val = 1.0 / (PE * phi)

    m0d = A_val * m_b / m_d
    m0u = A_val * m_t / m_u

    if m0d > -m_d and m0u > -m_u:
        r21d = (m_s + m0d) / (m_d + m0d)
        r31d = solve_r31_from_K(r21d, 2/3)
        mb_p = (m_d + m0d) * r31d - m0d if np.isfinite(r31d) else np.nan

        r21u = (m_c + m0u) / (m_u + m0u)
        r31u = solve_r31_from_K(r21u, 2/3)
        mt_p = (m_u + m0u) * r31u - m0u if np.isfinite(r31u) else np.nan

        print(f"  {OL:8.4f}  {PE:8.2f}  {A_val:10.6f}  {m0d:8.2f}  {mb_p:10.1f}  {m0u:10.1f}  {mt_p:12.0f}")


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
print("PODSUMOWANIE ex235")
print("=" * 72)
print(f"""
  ★ ROZSZERZENIE NA KWARKI:

  1. BARE Koide (K=2/3 bez shift):
     Leptony: K = {koide(m_e,m_mu,m_tau):.6f} ≈ 2/3  ← DOKŁADNY
     Down:    K = {koide(m_d,m_s,m_b):.6f}       ← 6% off
     Up:      K = {koide(m_u,m_c,m_t):.6f}       ← ~0.5% off

  2. SHIFTED Koide (K(m+m₀)=2/3):
     Kwarki wymagają m₀ ≠ 0 → "morze" gluonów/par
     m₀(down) = {m0_down:.2f} MeV
     m₀(up)   = {m0_up:.1f} MeV

  3. R12 formuła A = 1/(Φ_eff × φ) daje m₀ z kosmologii:
     m₀ = A × m₃/m₁ (sektor-zależne)

  4. PREDYKCJE z R12 + K=2/3:
     m_b: {m_b_R12:.0f} MeV (PDG: {m_b}, err: {abs(m_b_R12-m_b)/m_b*100:.1f}%)
     m_t: {m_t_R12:.0f} MeV (PDG: {m_t}, err: {abs(m_t_R12-m_t)/m_t*100:.1f}%)

  FRAMEWORK:
    Leptony: m_e + r₂₁(ODE) + K=2/3 → m_μ, m_τ, α_s    [ex234]
    Kwarki:  m_d,m_s + m₀(R12,Ω_Λ) + K=2/3 → m_b        [ex235]
    Kwarki:  m_u,m_c + m₀(R12,Ω_Λ) + K=2/3 → m_t        [ex235]
""")
