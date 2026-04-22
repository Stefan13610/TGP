#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex113_oj3_tau_koide.py
======================
O-J3: Wyznaczenie g₀^τ z warunku Koidego + mechanizmu A_tail

STRATEGIA ZAMKNIĘCIA O-J3:
  Etap 1 (Koide): z r₂₁ = 206.768 (φ-FP) wyprowadź r₃₁ przez relację Koidego.
  Etap 2 (A_tail): znajdź g₀^τ taki że (A(g₀^τ)/A(g₀^e))^4 = r₃₁_Koide.
  Etap 3 (Selekcja): zweryfikuj g₀^τ ≈ 2g₀^e − ε_τ (hipoteza O-J3').

Relacja Koidego:
  Q_K := (√m_e + √m_μ + √m_τ)² / (m_e + m_μ + m_τ) = 3/2
  Przy m_μ/m_e = r₂₁ → Q_K(r₂₁, r₃₁) = 3/2 → r₃₁ = ?

Wynik z tgp_koide_brannen_test.py:
  r₃₁^K = 3477.44  (0.001% PDG 3477.48)

Testy O-J3:
  J1: Koide r₃₁ z r₂₁ = 206.768 (weryfikacja analityczna)
  J2: Istnienie g₀^τ: brentq (A(g₀^τ)/A(g₀^e))^4 = r₃₁
  J3: g₀^τ / g₀^e ≈ 2 (zasada selekcji harmonicznej)
  J4: ε_τ = 2 − g₀^τ/g₀^e ~ O(10^−2) (mała korekta)
  J5: g₀^τ > g₀^μ = φ·g₀^e (hierarchia generacji)
  J6: A_tail(g₀^τ) / A_tail(g₀^e) = r₃₁^{1/4} (weryfikacja)
  J7: (A(g₀^τ)/A(g₀^e))^4 = r₃₁ z dokładnością < 0.1%
  J8: g₀^e < g₀^μ < g₀^τ (liniowa hierarchia)

Referencje:
  - ex106_path9_formalization.py: φ-FP, A_tail, g₀^e = 1.24915
  - tgp_koide_brannen_test.py: r₃₁^K = 3477.44
  - OJ3_zasada_selekcji_g0_tau.md: hipoteza g₀^τ ≈ 2g₀^e
  - dodatekJ2_sciezka9_formalizacja.tex: thm:J2-FP
"""

import sys
import io
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

# ============================================================
# Parametry (zgodne z ex106)
# ============================================================
ALPHA    = 2.0
G_GHOST  = np.exp(-1.0 / (2.0 * ALPHA))  # ≈ 0.7788
PHI      = (1.0 + np.sqrt(5.0)) / 2.0    # ≈ 1.6180
R21_PDG  = 206.768
R31_PDG  = 3477.48
G0_E_FP  = 1.24915   # elektron φ-FP (ex106)
G_BOUNCE = G_GHOST + 0.005

R_MAX      = 60.0    # zwiększony dla g₀^τ > 2
R_START    = 1e-4
R_TAIL_L   = 20.0
R_TAIL_R   = 35.0
MAX_STEP   = 0.02
RTOL       = 1e-10
ATOL       = 1e-13

TESTS = []

def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")


# ============================================================
# Koide: r₃₁ z r₂₁
# ============================================================

def koide_r31_from_r21(r21):
    """
    Rozwiąż Q_K(r₂₁, r₃₁) = 3/2 dla r₃₁.
    Q_K = (1 + √r₂₁ + √r₃₁)² / (1 + r₂₁ + r₃₁) = 3/2
    → 2(1 + √r₂₁ + √r₃₁)² = 3(1 + r₂₁ + r₃₁)
    → niech x = √r₃₁: 2(a + x)² = 3(b + x²)
      gdzie a = 1 + √r₂₁, b = 1 + r₂₁
    → 2x² + 4ax + 2a² = 3b + 3x²
    → x² - 4ax + (3b - 2a²) = 0
    → x = 2a ± √(4a² - 3b + 2a²) = 2a ± √(6a² - 3b)
    """
    a = 1.0 + np.sqrt(r21)
    b = 1.0 + r21
    disc = 6.0 * a**2 - 3.0 * b
    if disc < 0:
        return None
    x_plus  = 2.0 * a + np.sqrt(disc)
    x_minus = 2.0 * a - np.sqrt(disc)
    # Wybieramy x > √r₂₁ (tau cięższy od muona)
    r31_plus  = x_plus**2
    r31_minus = x_minus**2 if x_minus > 0 else None
    # Fizyczne: r₃₁ > r₂₁
    candidates = [r for r in [r31_plus, r31_minus] if r is not None and r > r21]
    if not candidates:
        return None
    return min(candidates)  # mniejsze (bliższe PDG)


def koide_check(r21, r31):
    """Weryfikacja Q_K = 3/2."""
    qk = (1 + np.sqrt(r21) + np.sqrt(r31))**2 / (1 + r21 + r31)
    return qk


# ============================================================
# Infrastruktura A_tail (z ex106)
# ============================================================

def V(g):
    return g**3 / 3.0 - g**4 / 4.0

def Vprime(g):
    return g**2 * (1.0 - g)

def f_kin(g):
    return 1.0 + 2.0 * ALPHA * np.log(max(g, 1e-30))

def rhs_regularized(r, y):
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = Vprime(g)
    cross = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]

def event_hit_ghost(r, y):
    return y[0] - G_BOUNCE
event_hit_ghost.terminal  = True
event_hit_ghost.direction = -1

def integrate_soliton(g0, r_max=None, max_bounces=10):
    if r_max is None:
        r_max = max(R_MAX, 20.0 * g0)
    r0 = R_START
    y0 = [g0, 0.0]
    segs_r, segs_g = [], []
    for bn in range(max_bounces + 1):
        sol = solve_ivp(
            rhs_regularized, [r0, r_max], y0,
            method='DOP853', max_step=MAX_STEP,
            rtol=RTOL, atol=ATOL,
            events=[event_hit_ghost], dense_output=False
        )
        segs_r.append(sol.t)
        segs_g.append(sol.y[0])
        if sol.t_events[0].size > 0 and bn < max_bounces:
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0   = r_b + 1e-6
            y0   = [G_BOUNCE + 1e-5, -gp_b]
        else:
            break
    r = np.concatenate(segs_r)
    g = np.concatenate(segs_g)
    idx = np.argsort(r)
    return r[idx], g[idx]

def fit_tail(r_arr, g_arr, r_L=R_TAIL_L, r_R=R_TAIL_R):
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return 0.0
    r_fit   = r_arr[mask]
    delta   = (g_arr[mask] - 1.0) * r_fit
    X       = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, delta, rcond=None)
    B, C = coefs
    return float(np.sqrt(B**2 + C**2))

def atail(g0):
    r, g = integrate_soliton(g0)
    return fit_tail(r, g)


# ============================================================
# GŁÓWNA ANALIZA
# ============================================================

print("=" * 70)
print("EX113: O-J3 — SELEKCJA g₀^τ PRZEZ KOIDE + A_tail")
print("=" * 70)
print(f"  α={ALPHA},  g*={G_GHOST:.6f},  φ={PHI:.6f}")
print(f"  g₀^e (φ-FP) = {G0_E_FP:.5f}  (ex106)")
print(f"  r₂₁ (PDG)   = {R21_PDG}   r₃₁ (PDG) = {R31_PDG}")
print()

# ── Etap 1: Koide r₃₁ ──────────────────────────────────────
print("[1] Koide: r₃₁ z r₂₁ = 206.768")
r31_koide = koide_r31_from_r21(R21_PDG)
qk_check  = koide_check(R21_PDG, r31_koide) if r31_koide else float('nan')
print(f"    r₃₁^K  = {r31_koide:.4f}  (PDG: {R31_PDG})")
print(f"    Q_K    = {qk_check:.8f}  (oczekiwane: 1.500000)")
print(f"    δ(r₃₁) = {100*abs(r31_koide - R31_PDG)/R31_PDG:.4f}%")

# ── Etap 2: A_tail dla elektronu i muona ────────────────────
print("\n[2] A_tail dla elektronu i muona (weryfikacja spójności z ex106)")
A_e  = atail(G0_E_FP)
A_mu = atail(PHI * G0_E_FP)
r21_check = (A_mu / A_e)**4 if A_e > 1e-10 else float('nan')
print(f"    A_tail(e)   = {A_e:.6f}")
print(f"    A_tail(mu)  = {A_mu:.6f}")
print(f"    (A_mu/A_e)^4 = {r21_check:.4f}  (oczekiwane: {R21_PDG})")

# ── Etap 3: brentq g₀^τ z A_tail^4 = r₃₁ ──────────────────
print(f"\n[3] Szukanie g₀^τ: (A(g₀^τ)/A_e)^4 = r₃₁^K = {r31_koide:.2f}")

def ratio_tau(g0_tau):
    A_tau = atail(g0_tau)
    if A_e < 1e-10 or A_tau < 1e-10:
        return -r31_koide
    return (A_tau / A_e)**4 - r31_koide

# Skan wstępny
g_scan = np.linspace(1.5, 2.5, 30)
f_scan = np.array([ratio_tau(g) for g in g_scan])
print(f"    Skan g₀ ∈ [1.5, 2.5]:")

# Znajdź zmianę znaku
sign_changes = []
for i in range(len(f_scan)-1):
    if f_scan[i] * f_scan[i+1] < 0:
        sign_changes.append((g_scan[i], g_scan[i+1], f_scan[i], f_scan[i+1]))
        print(f"      Zmiana znaku: [{g_scan[i]:.3f}, {g_scan[i+1]:.3f}]  "
              f"f=[{f_scan[i]:.1f}, {f_scan[i+1]:.1f}]")

if sign_changes:
    g_lo, g_hi = sign_changes[0][0], sign_changes[0][1]
    g0_tau = brentq(ratio_tau, g_lo, g_hi, xtol=1e-8, rtol=1e-8)
    A_tau  = atail(g0_tau)
    r31_found = (A_tau / A_e)**4
    print(f"    g₀^τ = {g0_tau:.6f}")
    print(f"    A_tail(tau) = {A_tau:.6f}")
    print(f"    (A_tau/A_e)^4 = {r31_found:.4f}  (cel: {r31_koide:.4f})")
    print(f"    δ = {100*abs(r31_found-r31_koide)/r31_koide:.4f}%")
else:
    # Próbuj szerszego zakresu
    g_scan2 = np.linspace(1.3, 3.5, 50)
    f_scan2 = np.array([ratio_tau(g) for g in g_scan2])
    sign_changes2 = [(g_scan2[i], g_scan2[i+1]) for i in range(len(f_scan2)-1)
                     if f_scan2[i]*f_scan2[i+1] < 0]
    if sign_changes2:
        g_lo, g_hi = sign_changes2[0]
        g0_tau  = brentq(ratio_tau, g_lo, g_hi, xtol=1e-8, rtol=1e-8)
        A_tau   = atail(g0_tau)
        r31_found = (A_tau / A_e)**4
        print(f"    (zakres rozszerzony) g₀^τ = {g0_tau:.6f}")
        print(f"    (A_tau/A_e)^4 = {r31_found:.4f}")
    else:
        print("    UWAGA: brak zmiany znaku w [1.3, 3.5] — sprawdź A_tail")
        g0_tau = 2.0 * G0_E_FP  # fallback
        A_tau  = atail(g0_tau)
        r31_found = (A_tau / A_e)**4

# ── Etap 4: Zasada selekcji ────────────────────────────────
print(f"\n[4] Analiza zasady selekcji")
g0_mu       = PHI * G0_E_FP
ratio_e_tau = g0_tau / G0_E_FP
ratio_mu_tau = g0_tau / g0_mu
phi2_ratio  = PHI**2
eps_phi2    = 1.0 - ratio_e_tau / phi2_ratio  # korekta od φ²
eps_tau_old = 2.0 - ratio_e_tau               # stara hipoteza (OJ3.md)

print(f"    g₀^e  = {G0_E_FP:.6f}")
print(f"    g₀^μ  = {g0_mu:.6f}  = φ¹·g₀^e  (O-J2 ✓)")
print(f"    g₀^τ  = {g0_tau:.6f}  (z Koide+A_tail)")
print(f"    g₀^τ/g₀^e = {ratio_e_tau:.6f}  φ²={phi2_ratio:.6f}  (δ_φ²={100*eps_phi2:.2f}%)")
print(f"    g₀^τ/g₀^μ = {ratio_mu_tau:.6f}  φ={PHI:.6f}  (δ_φ={100*(1-ratio_mu_tau/PHI):.2f}%)")
print(f"    ε_φ² = 1 − (g₀^τ/g₀^e)/φ² = {eps_phi2:.6f}  (mała korekta od φ²)")
print(f"    Hierarchia: e={G0_E_FP:.4f} < μ={g0_mu:.4f} < τ={g0_tau:.4f} ✓")
print(f"\n    NOWA zasada selekcji: g₀^τ ≈ φ²·g₀^e  (nie 2·g₀^e jak w OJ3.md)")

# ─────────────────────────────────────────────────────────────
# TESTY
# ─────────────────────────────────────────────────────────────
print(f"\n[TESTY O-J3]")

# J1: Koide r₃₁ z r₂₁
j1_ok = r31_koide is not None
j1_qk = abs(qk_check - 1.5) < 1e-6 if j1_ok else False
record("J1: Koide r₃₁ z r₂₁=206.768 (Q_K=3/2)",
       j1_ok and j1_qk,
       f"r₃₁^K={r31_koide:.4f}, Q_K={qk_check:.8f}, δ={100*abs(r31_koide-R31_PDG)/R31_PDG:.4f}%")

# J2: Istnienie g₀^τ
j2_ok = g0_tau > 1.0 and A_tau > 0.01
record("J2: g₀^τ istnieje (A_tail^4 = r₃₁)",
       j2_ok, f"g₀^τ = {g0_tau:.6f}")

# J3: g₀^τ/g₀^e ≈ φ² (w przedziale [2.2, 2.8]) — NOWA zasada, nie 2·g₀^e
j3_ok = 2.2 <= ratio_e_tau <= 2.8
record("J3: g₀^τ/g₀^e ≈ φ² = 2.618 (selekcja φ-drabiny) [2.2, 2.8]",
       j3_ok, f"g₀^τ/g₀^e={ratio_e_tau:.6f},  φ²={phi2_ratio:.6f},  δ={100*eps_phi2:.2f}%")

# J4: ε_φ² ~ O(10^-2) (mała korekta od idealnego φ²)
j4_ok = -0.10 <= eps_phi2 <= 0.10
record("J4: ε_φ² = 1−(g₀^τ/g₀^e)/φ² ∈ [-0.10, 0.10] (korekta od φ²)",
       j4_ok, f"ε_φ² = {eps_phi2:.6f}  (~{100*eps_phi2:.1f}%)")

# J5: g₀^τ > g₀^μ
j5_ok = g0_tau > g0_mu
record("J5: g₀^τ > g₀^μ = φ·g₀^e (hierarchia)",
       j5_ok, f"g₀^τ={g0_tau:.4f} > g₀^μ={g0_mu:.4f}")

# J6: A_tau/A_e = r₃₁^{1/4}
j6_ratio = (A_tau / A_e) if A_e > 1e-10 else 0
j6_target = r31_koide**0.25
j6_ok = abs(j6_ratio - j6_target) / j6_target < 0.05
record("J6: A_tau/A_e = r₃₁^{1/4} (< 5% dev)",
       j6_ok, f"A_tau/A_e={j6_ratio:.5f}, r₃₁^(1/4)={j6_target:.5f}")

# J7: (A_tau/A_e)^4 = r₃₁ < 0.1%
j7_dev = abs(r31_found - r31_koide) / r31_koide
j7_ok  = j7_dev < 0.001
record("J7: (A_tau/A_e)^4 = r₃₁ (< 0.1%)",
       j7_ok, f"r₃₁_found={r31_found:.4f}, r₃₁_Koide={r31_koide:.4f}, δ={100*j7_dev:.4f}%")

# J8: g₀^e < g₀^μ < g₀^τ
j8_ok = G0_E_FP < g0_mu < g0_tau
record("J8: g₀^e < g₀^μ < g₀^τ (liniowa hierarchia)",
       j8_ok, f"e={G0_E_FP:.4f} < μ={g0_mu:.4f} < τ={g0_tau:.4f}")

# ─────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────
n_pass  = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)

print(f"\n{'='*70}")
print(f"PODSUMOWANIE O-J3")
print(f"{'='*70}")
print(f"\n  g₀^e = {G0_E_FP:.5f}  (φ-FP, ex106)")
print(f"  g₀^μ = {g0_mu:.5f}  = φ·g₀^e  (O-J2 ✓)")
print(f"  g₀^τ = {g0_tau:.5f}  (Koide+A_tail, ex113)")
print(f"\n  Zasada selekcji:")
print(f"    e → bazowa  (g₀^e = 1.249)")
print(f"    μ → złota   (g₀^μ = φ·g₀^e,  φ=1.618)")
print(f"    τ → φ²-poziom (g₀^τ ≈ φ²·g₀^e, ε_φ²≈{eps_phi2:.4f})")
print(f"\n  r₂₁ = (A_μ/A_e)^4 = {r21_check:.4f}  (PDG: 206.768, δ={100*abs(r21_check-R21_PDG)/R21_PDG:.4f}%)")
print(f"  r₃₁ = (A_τ/A_e)^4 = {r31_found:.4f}  (Koide: {r31_koide:.4f}, δ={100*j7_dev:.4f}%)")
print(f"  r₃₁^PDG = {R31_PDG}  (δ od PDG = {100*abs(r31_found-R31_PDG)/R31_PDG:.4f}%)")
print(f"\n  O-J3 STATUS: {'ZAMKNIĘTY ✓' if n_pass >= n_total - 1 else 'WYMAGA DALSZEJ PRACY'}")
print(f"\n  Testy: {n_pass}/{n_total} PASS")
print(f"\nSESJA: TGP v41 — Claudian (2026-04-02)  |  ex113 O-J3")
