#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex125_qk32_phase_derivation.py
================================
T-OP1: Próba wyprowadzenia Q_K = 3/2 z dynamiki solitonu TGP.

PYTANIE (T-OP1):
  Dlaczego Q_K(m_e, m_μ, m_τ) = 3/2 dokładnie?
  φ-drabina daje Q_K ≈ 1.472 (1.9% za mało). Skąd korekta ξ*=2.553?

HIPOTEZA FAZOWA (Phase-Shift Hypothesis, PSH):
  Ogon solitonu: (g_k(r) - 1)·r ≈ A_k·cos(r + δ_k)   dla r ∈ [20,35]
  gdzie:  B_k = A_k·cos(δ_k),  C_k = -A_k·sin(δ_k),  δ_k = atan2(-C_k, B_k)

  Jeśli fazy δ_k tworzą postęp arytmetyczny z krokiem 2π/3:
    δ_μ - δ_e = δ_τ - δ_μ = 2π/3  (mod 2π)
  to wektory (A_k·e^{i·δ_k}) obracają się o 2π/3 → warunek Z₃ → Q_K = 3/2.

PLAN:
  1. Solitony dla g₀^e=1.24915, g₀^μ=φ·g₀^e=2.02064, g₀^τ=3.18912
  2. Fit ogona: (g-1)·r = B·cos(r) + C·sin(r) → A, δ dla każdej generacji
  3. Test PSH: δ_τ - δ_μ ≈? δ_μ - δ_e (mod 2π), jak blisko 2π/3?
  4. Test Brannena: A_k ∝ √m_k? → czy forma Brannena wyłania się naturalnie?
  5. Bliskie trafienie: r₂₁/r* = 206.77/22.956 ≈ 9.007 ≈ 9 = 3²?
  6. Skan g₀ → δ(g₀): jak zmienia się faza z g₀? Liniowo? Z jakimś krokiem?
  7. Wniosek o statusie T-OP1

TESTY N1..N14:
  N1:  Solitony dla e, μ, τ — A_tail > 0 dla wszystkich
  N2:  RMSE/A_tail < 10% (jakość fitu ogona)
  N3:  A_e : A_μ : A_τ — zgodność z √m (tj. A_k ∝ √√m_k = m_k^{1/4}?)
  N4:  Q_K z A_tail (przez Brannena): Q_K(A-ratios) = ?
  N5:  Fazy δ_e, δ_μ, δ_τ — obliczone
  N6:  PSH: Δδ₁₂ = δ_μ - δ_e i Δδ₂₃ = δ_τ - δ_μ (mod 2π) — różnice
  N7:  PSH precyzja: |Δδ₁₂ - 2π/3| i |Δδ₂₃ - 2π/3| w radianach
  N8:  PSH precyzja: |Δδ₁₂ - Δδ₂₃| < 0.5 rad (arytmetyczność)?
  N9:  r₂₁/r*: czy blisko całkowitej liczby N²?
  N10: Skan δ(g₀) dla g₀ ∈ [1.1, 3.5] — monotoniczny?
  N11: Krok δ na φ-drabinie: δ(φ·g₀) - δ(g₀) = ?
  N12: A_tail wyznaczone przez fit liniowy B/C vs g₀ (gęstość)
  N13: Bezpośrednie Q_K z (A_e, A_μ, A_τ) vs pomiarowe Q_K=1.5
  N14: Wniosek: status T-OP1 (OPEN/PARTIAL/RESOLVED)

Referencje:
  - ex114: mapa fazowa B, C, A_tail vs g₀
  - ex115: E_core ∝ A_tail^n, n≈2.15
  - ex119: Q_K = 3/(1+r²/2) → r=√2 ↔ Q_K=3/2
  - ex117: r* = 22.956, wieża leptonów
  - ex118: faktoryzacja Z₃: P(u) = (u²-5u+1)(u²+u+1)
"""

import sys
import io
import warnings
import math
import cmath
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

# ============================================================
# Stałe globalne
# ============================================================
ALPHA    = 2.0
PHI      = (1.0 + math.sqrt(5.0)) / 2.0     # ≈ 1.61803
G_GHOST  = math.exp(-1.0 / (2.0 * ALPHA))    # ≈ 0.77880
G_BOUNCE = G_GHOST + 0.005

# g₀ dla trzech generacji (z ex113/ex115)
G0_E   = 1.24915
G0_MU  = PHI * G0_E    # ≈ 2.02064
G0_TAU = 3.18912

# Masy PDG (MeV)
M_E_MEV   = 0.510999
M_MU_MEV  = 105.6584
M_TAU_MEV = 1776.86

# Okno ogona (standard z ex114)
R_TAIL_L = 20.0
R_TAIL_R = 35.0
R_MAX    = 60.0
R_START  = 1e-4
MAX_STEP = 0.02
RTOL     = 1e-10
ATOL     = 1e-13

# Znany r* (stały punkt Koide, z ex117)
R_STAR = (23.0 + 5.0 * math.sqrt(21.0)) / 2.0   # ≈ 22.9564

TESTS = []

def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")


# ============================================================
# ODE solitonu (ten sam co ex114/ex115)
# ============================================================

def Vprime(g):
    return g**2 * (1.0 - g)

def f_kin(g):
    return 1.0 + 2.0 * ALPHA * math.log(max(g, 1e-30))

def rhs(r, y):
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = Vprime(g)
    cross   = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]

def event_ghost(r, y):
    return y[0] - G_BOUNCE
event_ghost.terminal  = True
event_ghost.direction = -1

def integrate_soliton(g0, r_max=None, max_bounces=12):
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    r0 = R_START
    y0 = [g0, 0.0]
    segs_r, segs_g = [], []

    for bn in range(max_bounces + 1):
        sol = solve_ivp(
            rhs, [r0, r_max], y0,
            method='DOP853',
            max_step=MAX_STEP, rtol=RTOL, atol=ATOL,
            events=[event_ghost], dense_output=False,
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

    r_all = np.concatenate(segs_r)
    g_all = np.concatenate(segs_g)
    idx   = np.argsort(r_all)
    return r_all[idx], g_all[idx]


def fit_tail(r_arr, g_arr, r_L=R_TAIL_L, r_R=R_TAIL_R):
    """
    Fit: (g-1)·r = B·cos(r) + C·sin(r)
    Zwraca: A=√(B²+C²), delta=atan2(-C,B), B, C, rmse
    (konwencja: (g-1)·r ≈ A·cos(r + delta))
    """
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return 0.0, 0.0, 0.0, 0.0, float('nan')
    r_f  = r_arr[mask]
    y_f  = (g_arr[mask] - 1.0) * r_f
    X    = np.column_stack([np.cos(r_f), np.sin(r_f)])
    coef, _, _, _ = np.linalg.lstsq(X, y_f, rcond=None)
    B, C  = float(coef[0]), float(coef[1])
    y_hat = B * np.cos(r_f) + C * np.sin(r_f)
    rmse  = float(np.sqrt(np.mean((y_f - y_hat)**2)))
    A     = float(math.sqrt(B**2 + C**2))
    delta = float(math.atan2(-C, B))   # konwencja A·cos(r+delta)
    return A, delta, B, C, rmse


def koide_qk(m1, m2, m3):
    s = math.sqrt(m1) + math.sqrt(m2) + math.sqrt(m3)
    return s**2 / (m1 + m2 + m3)


def wrap_diff(a, b):
    """Różnica b-a zawijana do (-pi, pi]."""
    d = b - a
    return (d + math.pi) % (2 * math.pi) - math.pi


# ============================================================
# Sekcja 1: Obliczenie solitonów i ogonów dla e, μ, τ
# ============================================================
print("=" * 72)
print("EX125: PRÓBA WYPROWADZENIA Q_K=3/2 Z FAZY OGONA SOLITONU TGP")
print("       (T-OP1 Investigation)")
print("=" * 72)
print()

GENS = [
    ("e",   G0_E,   M_E_MEV),
    ("mu",  G0_MU,  M_MU_MEV),
    ("tau", G0_TAU, M_TAU_MEV),
]

print("[1] OBLICZANIE SOLITONÓW DLA (e, mu, tau)")
print("-" * 50)

tail_data = {}
for (lep, g0, mass) in GENS:
    r_arr, g_arr = integrate_soliton(g0)
    A, delta, B, C, rmse = fit_tail(r_arr, g_arr)
    tail_data[lep] = {
        'g0': g0, 'mass': mass,
        'A': A, 'delta': delta, 'B': B, 'C': C, 'rmse': rmse
    }
    print(f"  {lep:>3}: g₀={g0:.5f}  A={A:.6f}  δ={delta:.6f} rad "
          f"  B={B:.6f}  C={C:.6f}  RMSE={rmse:.6e}")

A_e   = tail_data['e']['A']
A_mu  = tail_data['mu']['A']
A_tau = tail_data['tau']['A']
d_e   = tail_data['e']['delta']
d_mu  = tail_data['mu']['delta']
d_tau = tail_data['tau']['delta']

# N1: A_tail > 0 dla wszystkich
record("N1: A_tail > 0 dla e, mu, tau",
       A_e > 0 and A_mu > 0 and A_tau > 0,
       f"A_e={A_e:.5f}, A_mu={A_mu:.5f}, A_tau={A_tau:.5f}")

# N2: RMSE/A < 10%
rmse_rel = [tail_data[k]['rmse'] / max(tail_data[k]['A'], 1e-10)
            for k in ('e', 'mu', 'tau')]
record("N2: RMSE/A_tail < 10% dla wszystkich generacji",
       all(rr < 0.10 for rr in rmse_rel),
       f"RMSE/A = {rmse_rel[0]:.4f}, {rmse_rel[1]:.4f}, {rmse_rel[2]:.4f}")


# ============================================================
# Sekcja 2: Stosunek A_tail vs √m_k
# ============================================================
print()
print("[2] STOSUNEK A_TAIL vs √m_k  (test: A_k ∝ m_k^p)")
print("-" * 50)

sqm_e   = math.sqrt(M_E_MEV)
sqm_mu  = math.sqrt(M_MU_MEV)
sqm_tau = math.sqrt(M_TAU_MEV)

# Jeśli A_k ∝ √m_k, to A_mu/A_e = sqrt(m_mu/m_e)
ratio_A_emu  = A_mu  / A_e
ratio_A_etau = A_tau / A_e
ratio_sqm_emu  = math.sqrt(M_MU_MEV / M_E_MEV)
ratio_sqm_etau = math.sqrt(M_TAU_MEV / M_E_MEV)

# Faktyczny wykładnik z dopasowania log-log
A_vals = np.array([A_e, A_mu, A_tau])
m_vals = np.array([M_E_MEV, M_MU_MEV, M_TAU_MEV])
log_A  = np.log(A_vals)
log_m  = np.log(m_vals)
p_fit  = float(np.polyfit(log_m, log_A, 1)[0])

print(f"  A_mu/A_e   = {ratio_A_emu:.4f}   (√(m_mu/m_e) = {ratio_sqm_emu:.4f})")
print(f"  A_tau/A_e  = {ratio_A_etau:.4f}  (√(m_tau/m_e) = {ratio_sqm_etau:.4f})")
print(f"  Wykładnik p (A∝m^p): p = {p_fit:.4f}  [√m → p=0.25]")

# N3: p blisko 0.25 (czyli A ∝ m^{1/4} = √√m)?
record("N3: wykładnik p (A_tail vs m) w [0.15, 0.40]",
       0.15 < p_fit < 0.40,
       f"p={p_fit:.4f}  (oczekiwane ~0.25 jeśli A ∝ m^{{1/4}})")


# ============================================================
# Sekcja 3: Q_K z A_tail przez formułę Brannena
# ============================================================
print()
print("[3] Q_K Z A_TAIL PRZEZ FORMUŁĘ BRANNENA")
print("-" * 50)

# Jeśli m_k ∝ A_k^4 (z p=0.25), to √m_k ∝ A_k².
# Q_K = (Σ√m_k)²/(Σm_k) = (ΣA_k²)²/(ΣA_k⁴) = koide_qk(A_e^4, ...)
# Parametry Brannena dla A_tail:
# √m_k ∝ A_k² → fit Brannena dla A_k²:
#   A_k² = M_B·(1 + r_B·cos(θ_B + 2πk/3))
# Poprawna formuła DFT:
#   Σ_{k=0}^{2} ε_k·e^{-i2πk/3} = (3/2)·r_B·e^{iθ_B}
# → r_B = |DFT|·2/3,  θ_B = arg(DFT)

A_vals3 = np.array([A_e, A_mu, A_tau])
sqm_B   = A_vals3**2     # ∝ √m_k (Brannena: √m_k = M_B(1+r·cos))
M_B     = float(np.mean(sqm_B))
eps     = sqm_B / M_B - 1.0  # ε_k = r_B·cos(θ_B + 2πk/3)

# Q_K z bezpośrednich mas PDG
qk_pdg = koide_qk(M_E_MEV, M_MU_MEV, M_TAU_MEV)
# Q_K z A_tail (m_k ∝ A_k^4, sprawdza czy A_k ∝ m_k^{1/4} jest spójne z Q_K)
qk_from_A4 = koide_qk(A_e**4, A_mu**4, A_tau**4)

print(f"  Q_K (PDG masy)                   = {qk_pdg:.6f}")
print(f"  Q_K z A_tail^4 (m∝A⁴, p=0.25)   = {qk_from_A4:.6f}")

fourier_eps = sum(eps[k] * cmath.exp(-2j * math.pi * k / 3.0) for k in range(3))
r_B_atail   = abs(fourier_eps) * 2.0 / 3.0
theta_B_atail = math.atan2(fourier_eps.imag, fourier_eps.real)

# Q_K z r_B: Q_K = 3/(1 + r_B²/2)
qk_from_rB = 3.0 / (1.0 + r_B_atail**2 / 2.0)
r_B_sqrt2  = math.sqrt(2.0)

print(f"  r_Brannen(A_tail)               = {r_B_atail:.6f}  [√2 = {r_B_sqrt2:.6f}]")
print(f"  θ_Brannen(A_tail)               = {math.degrees(theta_B_atail):.2f}°")
print(f"  Q_K(r_Brannen(A_tail))          = {qk_from_rB:.6f}  [cel: 1.5]")

# N4: Q_K z A_tail
record("N4: Q_K z A_tail bliskie 3/2 (±10%)?",
       abs(qk_from_rB - 1.5) / 1.5 < 0.10,
       f"Q_K={qk_from_rB:.4f}, δ={100*(qk_from_rB-1.5)/1.5:+.2f}%")


# ============================================================
# Sekcja 4: Analiza faz δ_k
# ============================================================
print()
print("[4] ANALIZA FAZ SOLITONU δ_e, δ_mu, δ_tau")
print("-" * 50)

print(f"  δ_e   = {d_e:.6f} rad  = {math.degrees(d_e):.3f}°")
print(f"  δ_mu  = {d_mu:.6f} rad  = {math.degrees(d_mu):.3f}°")
print(f"  δ_tau = {d_tau:.6f} rad  = {math.degrees(d_tau):.3f}°")

# Różnice faz (zawijane do (-π,π])
diff_12 = wrap_diff(d_e, d_mu)      # δ_mu - δ_e (mod 2π)
diff_23 = wrap_diff(d_mu, d_tau)    # δ_tau - δ_mu (mod 2π)

two_pi_3 = 2.0 * math.pi / 3.0     # ≈ 2.0944 rad = 120°

print(f"\n  Δδ₁₂ = δ_μ - δ_e    = {diff_12:.6f} rad = {math.degrees(diff_12):.3f}°")
print(f"  Δδ₂₃ = δ_τ - δ_μ    = {diff_23:.6f} rad = {math.degrees(diff_23):.3f}°")
print(f"  2π/3                 = {two_pi_3:.6f} rad = 120.000°")
print(f"  |Δδ₁₂ - 2π/3|       = {abs(diff_12 - two_pi_3):.6f} rad")
print(f"  |Δδ₂₃ - 2π/3|       = {abs(diff_23 - two_pi_3):.6f} rad")
print(f"  |Δδ₁₂ - Δδ₂₃|       = {abs(diff_12 - diff_23):.6f} rad")

# N5: fazy obliczone (numeryczne)
record("N5: fazy δ obliczone dla e, mu, tau",
       all(math.isfinite(d) for d in [d_e, d_mu, d_tau]),
       f"δ_e={d_e:.4f}, δ_μ={d_mu:.4f}, δ_τ={d_tau:.4f}")

# N6: PSH — różnice faz
record("N6: PSH — Δδ₁₂ i Δδ₂₃ obliczone",
       math.isfinite(diff_12) and math.isfinite(diff_23),
       f"Δδ₁₂={math.degrees(diff_12):.2f}°, Δδ₂₃={math.degrees(diff_23):.2f}°")

# N7: Precyzja PSH: jak blisko 2π/3?
PSH_tol_rad = 0.5  # 0.5 rad ≈ 28.6°
psh_err_12  = abs(diff_12 - two_pi_3)
psh_err_23  = abs(diff_23 - two_pi_3)
# Sprawdzamy też -2π/3 (mogą iść w odwrotną stronę)
psh_err_12m = abs(diff_12 + two_pi_3)
psh_err_23m = abs(diff_23 + two_pi_3)
psh_err_12  = min(psh_err_12, psh_err_12m)
psh_err_23  = min(psh_err_23, psh_err_23m)

psh_close = (psh_err_12 < PSH_tol_rad) and (psh_err_23 < PSH_tol_rad)
record("N7: PSH — |Δδ - 2π/3| < 0.5 rad dla obu różnic",
       psh_close,
       f"|Δδ₁₂-2π/3|={psh_err_12:.4f} rad, |Δδ₂₃-2π/3|={psh_err_23:.4f} rad")

# N8: Arytmetyczność: |Δδ₁₂ - Δδ₂₃| < 0.5 rad
arith_diff = abs(diff_12 - diff_23)
arith_ok   = arith_diff < 0.5
record("N8: Arytmetyczność faz — |Δδ₁₂ - Δδ₂₃| < 0.5 rad",
       arith_ok,
       f"|Δδ₁₂ - Δδ₂₃| = {arith_diff:.4f} rad = {math.degrees(arith_diff):.2f}°")


# ============================================================
# Sekcja 5: Bliskie trafienie r₂₁ / r*
# ============================================================
print()
print("[5] BLISKIE TRAFIENIE: r₂₁ / r*")
print("-" * 50)

R21_PDG  = M_MU_MEV / M_E_MEV   # masa ratio ≈ 206.768 (używamy masy, nie A_tail)
r21_over_rstar = R21_PDG / R_STAR
N_sq     = 3.0**2   # = 9

print(f"  r₂₁ (m_μ/m_e) = {R21_PDG:.4f}")
print(f"  r* (Koide FP) = {R_STAR:.6f}")
print(f"  r₂₁ / r*      = {r21_over_rstar:.6f}")
print(f"  N² = 3²       = {N_sq:.1f}")
print(f"  |r₂₁/r* - 9|  = {abs(r21_over_rstar - 9.0):.6f}")

# Czy r₂₁ / r* jest bliskie liczbie całkowitej?
nearest_int = round(r21_over_rstar)
frac_part   = abs(r21_over_rstar - nearest_int)
record("N9: r₂₁/r* bliskie całkowitej liczbie (±2%)",
       frac_part < 0.02 * nearest_int,
       f"r₂₁/r* = {r21_over_rstar:.4f} ≈ {nearest_int} (δ={100*frac_part/nearest_int:.3f}%)")


# ============================================================
# Sekcja 6: Skan δ(g₀) — jak zmienia się faza z g₀?
# ============================================================
print()
print("[6] SKAN δ(g₀) DLA g₀ ∈ [1.1, 3.5]")
print("-" * 50)

G0_SCAN = np.linspace(1.1, 3.5, 30)
scan_A   = []
scan_d   = []

for g0_s in G0_SCAN:
    r_s, g_s = integrate_soliton(g0_s)
    A_s, d_s, B_s, C_s, rmse_s = fit_tail(r_s, g_s)
    scan_A.append(A_s)
    scan_d.append(d_s)

scan_A = np.array(scan_A)
scan_d = np.array(scan_d)

# N10: Monotoniczny δ(g₀)?
# Sprawdzamy czy δ jest monotonicznie rosnący lub malejący
d_diffs = np.diff(scan_d)
n_pos   = int(np.sum(d_diffs > 0))
n_neg   = int(np.sum(d_diffs < 0))
monotone = (n_pos < 3 or n_neg < 3)  # w większości jeden kierunek
record("N10: δ(g₀) prawie monotoniczna (≤2 odwrócenia kierunku)",
       monotone,
       f"n_rosnący={n_pos}, n_malejący={n_neg} (z {len(d_diffs)} kroków)")

# Krok δ na φ-drabinie: δ(φ·g₀) - δ(g₀)
# Szukamy par (g₀, φ·g₀) w skanie
phi_steps = []
for i, g0_s in enumerate(G0_SCAN):
    g0_phi = PHI * g0_s
    if g0_phi <= G0_SCAN[-1]:
        # Interpolacja
        j = int(np.searchsorted(G0_SCAN, g0_phi))
        if 0 < j < len(G0_SCAN):
            # Liniowa interpolacja
            t = (g0_phi - G0_SCAN[j-1]) / (G0_SCAN[j] - G0_SCAN[j-1])
            d_phi = scan_d[j-1] + t * (scan_d[j] - scan_d[j-1])
            phi_steps.append(wrap_diff(scan_d[i], d_phi))

if phi_steps:
    phi_step_mean = float(np.mean(phi_steps))
    phi_step_std  = float(np.std(phi_steps))
    print(f"  Krok fazowy Δδ(φ·g₀): mean={math.degrees(phi_step_mean):.2f}°  "
          f"std={math.degrees(phi_step_std):.2f}°")
    print(f"  2π/3 = 120.00°")
    print(f"  |Δδ(φ·g₀) - 2π/3| = {math.degrees(abs(phi_step_mean - two_pi_3)):.2f}°")

    # N11: Czy krok δ na φ-drabinie ≈ 2π/3?
    record("N11: Krok fazowy Δδ(φ-drabina) ≈ 2π/3 (±30°)?",
           abs(phi_step_mean - two_pi_3) < math.radians(30),
           f"Δδ={math.degrees(phi_step_mean):.2f}°, "
           f"|Δδ-120°|={math.degrees(abs(phi_step_mean-two_pi_3)):.2f}°")
else:
    record("N11: Krok fazowy Δδ(φ-drabina) — brak danych", False, "skan za wąski")


# ============================================================
# Sekcja 7: Bezpośrednie Q_K z A_tail jako √m
# ============================================================
print()
print("[7] BEZPOŚREDNIE Q_K Z A_TAIL")
print("-" * 50)

# Hipoteza 1: m_k ∝ A_k^4 (z A ∝ m^{1/4}) → Q_K = (ΣA_k^2)^2 / (ΣA_k^4)
# (bo √m_k ∝ A_k^2, Σm_k ∝ ΣA_k^4, Q_K=(Σ√m)^2/(Σm) bez czynnika N)
qk_h1 = koide_qk(A_e**4, A_mu**4, A_tau**4)

# Hipoteza 2: m_k ∝ A_k^2 → Q_K = (ΣA_k)^2/(ΣA_k^2)
qk_h2 = koide_qk(A_e**2, A_mu**2, A_tau**2)

# Hipoteza 3: m_k ∝ A_k^(1/p_fit) (z wykładnika p_fit)
exp_h3 = 1.0 / p_fit
m_h3   = A_vals3**exp_h3
qk_h3  = koide_qk(m_h3[0], m_h3[1], m_h3[2])

print(f"  Q_K (PDG)                   = {qk_pdg:.6f}")
print(f"  Q_K H1 (m∝A⁴, √m∝A²)       = {qk_h1:.6f}  [δ={100*(qk_h1-1.5)/1.5:+.2f}%]")
print(f"  Q_K H2 (m∝A², √m∝A)         = {qk_h2:.6f}  [δ={100*(qk_h2-1.5)/1.5:+.2f}%]")
print(f"  Q_K H3 (m∝A^(1/p), p={p_fit:.3f}) = {qk_h3:.6f}  [δ={100*(qk_h3-1.5)/1.5:+.2f}%]")

# N12: Najlepsza hipoteza daje Q_K bliskie 1.5?
best_h = min([(abs(qk_h1-1.5), 'H1', qk_h1),
              (abs(qk_h2-1.5), 'H2', qk_h2),
              (abs(qk_h3-1.5), 'H3', qk_h3)])
record("N12: Najlepsza hipoteza A→m daje Q_K w [1.4, 1.6]",
       1.4 < best_h[2] < 1.6,
       f"Najlepsza: {best_h[1]}, Q_K={best_h[2]:.4f}")

# N13: Bezpośrednie Q_K PDG potwierdzone
record("N13: Q_K(PDG) ∈ [1.4998, 1.5002] (reprodukcja standardowa)",
       abs(qk_pdg - 1.5) < 0.002,
       f"Q_K(PDG)={qk_pdg:.6f}")


# ============================================================
# Sekcja 8: Podsumowanie statusu T-OP1
# ============================================================
print()
print("[8] PODSUMOWANIE STATUSU T-OP1")
print("=" * 72)

print("""
HIERARCHIA WYJAŚNIEŃ Q_K = 3/2:

  POZIOM 0 — FAKT OBSERWACYJNY:
    Q_K(m_e, m_μ, m_τ) = 3/2  [do 0.03%]

  POZIOM 1 — ALGEBRAICZNIE KONIECZNE (ex117):
    Warunek FP wieży Q_K(1, u², u⁴) = 3/2 → P(u) = u⁴-4u³-3u²-4u+1
    → r* = (23+5√21)/2 ≈ 22.956
    [ZALICZONY — ex117, ex118]

  POZIOM 2 — RÓWNOWAŻNE WARUNKI (ex119):
    Q_K=3/2 ↔ Brannen(r=√2) ↔ CV(√m)=1 ↔ r=√(N-1) dla N=3
    [ALGEBRAICZNIE POKAZANE — ex119]

  POZIOM 3 — SKĄD φ-DRABINA DAJE Q_K≈1.472? (ex114):
    g₀^μ = φ·g₀^e → A_tail(μ)/A_tail(e) = φ_A ≠ √2
    Korekta: ξ* = 2.553 zamiast φ² = 2.618 (2.5% różnica)
    [CZĘŚCIOWO — φ-drabina znana, ξ* niewyjaśnione]
""")

# PSH status
print("  HIPOTEZA FAZOWA (PSH) — wynik ex125:")
print(f"    Δδ₁₂ (δ_μ-δ_e)   = {math.degrees(diff_12):.2f}°  [cel: ±120°]")
print(f"    Δδ₂₃ (δ_τ-δ_μ)   = {math.degrees(diff_23):.2f}°  [cel: ±120°]")
print(f"    |Δδ₁₂-2π/3|      = {math.degrees(psh_err_12):.2f}°")
print(f"    |Δδ₂₃-2π/3|      = {math.degrees(psh_err_23):.2f}°")

psh_total_err = psh_err_12 + psh_err_23
if psh_total_err < 0.2:
    psh_verdict = "PSH: POTWIERDZONE (±6°) — Q_K=3/2 wynika z Z₃ fazy solitonu!"
    psh_status  = "PARTIAL RESOLUTION of T-OP1"
elif psh_total_err < 0.5:
    psh_verdict = "PSH: CZĘŚCIOWO — fazy niezbyt blisko Z₃ (±14°)"
    psh_status  = "WEAK EVIDENCE — T-OP1 nadal OPEN"
elif psh_total_err < 1.0:
    psh_verdict = "PSH: NEGATYWNY — fazy zbyt daleko od Z₃ (±30°)"
    psh_status  = "NEGATIVE — T-OP1 OPEN, hipoteza fazowa nie działa"
else:
    psh_verdict = "PSH: OBALONY — fazy niezgodne z Z₃"
    psh_status  = "REFUTED — T-OP1 OPEN, hipoteza fazowa błędna"

print(f"\n  {psh_verdict}")
print(f"  → STATUS: {psh_status}")

print("""
BLISKIE TRAFIENIE:
""")
print(f"  r₂₁/r* = {r21_over_rstar:.6f} ≈ {nearest_int}  "
      f"(|r₂₁/r* - 9| = {abs(r21_over_rstar-9):.4f})")

r21_rstar_pct = 100 * abs(r21_over_rstar - 9.0) / 9.0
print(f"  Dokładność: {r21_rstar_pct:.3f}% od 9 = 3²")
if r21_rstar_pct < 0.1:
    print("  → BLISKIE TRAFIENIE < 0.1% : może być nieprzypadkowe")
elif r21_rstar_pct < 1.0:
    print("  → BLISKIE TRAFIENIE 0.1-1% : prawdopodobnie numeryczne")
else:
    print("  → PRZYPADKOWA ZBIEŻNOŚĆ > 1%")

# N14: Wniosek T-OP1
psh_resolved = psh_total_err < 0.2
record("N14: T-OP1 — PSH pokazuje Q_K=3/2 z fazy solitonu (wynik główny)",
       psh_resolved,
       psh_verdict)


# ============================================================
# Końcowe podsumowanie
# ============================================================
print()
print("=" * 72)
print("WYNIKI TESTÓW")
print("=" * 72)
passed = sum(1 for _, p, _ in TESTS if p)
total  = len(TESTS)
for name, passed_t, detail in TESTS:
    mark = "PASS" if passed_t else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        print(f"         {detail}")

print()
print(f"WYNIK: {passed}/{total} testów PASS")
print()

# Szczegółowe podsumowanie ilościowe
print("SZCZEGÓŁOWE DANE ILOŚCIOWE:")
print(f"  A_tail:   A_e={A_e:.5f}, A_μ={A_mu:.5f}, A_τ={A_tau:.5f}")
print(f"  Fazy:     δ_e={math.degrees(d_e):.2f}°,  δ_μ={math.degrees(d_mu):.2f}°,  "
      f"δ_τ={math.degrees(d_tau):.2f}°")
print(f"  Δδ:       Δδ₁₂={math.degrees(diff_12):.2f}°,  Δδ₂₃={math.degrees(diff_23):.2f}°")
print(f"  PSH err:  |Δδ₁₂-120°|={math.degrees(psh_err_12):.2f}°,  "
      f"|Δδ₂₃-120°|={math.degrees(psh_err_23):.2f}°")
print(f"  r_Brannen(A_tail): {r_B_atail:.5f}  [√2={math.sqrt(2):.5f}]")
print(f"  Q_K:      PDG={qk_pdg:.4f}, H1={qk_h1:.4f}, H2={qk_h2:.4f}")
print(f"  r₂₁/r*:  {r21_over_rstar:.5f}  ≈ {nearest_int}  (δ={r21_rstar_pct:.3f}%)")
print()
print(f"STATUS T-OP1: {psh_status}")
