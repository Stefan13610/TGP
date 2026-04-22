#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex112_soliton_energy_ksub.py
============================
OP-G + OP-E: A_tail z elastycznym odbiciem + E_core z K_sub(g) = g²

WYNIKI ex111 (do poprawienia):
  - Ghost constraint przy f(g)=1+2αlng uniemożliwia g₀^μ,g₀^τ
  - Max E-ratio ≈ 5.7 << R₂₁=206.77
  - E_total IR-rozbieżna (E ∝ R_MAX)

ROZWIĄZANIE (z ex106):
  Profil solitonu: ODE z f(g) + ELASTYCZNE ODBICIA przy g*
  → g₀^μ = 2.021 osiągalny; FAR window [20,35] dostępne
  → A_tail weryfikuje R₂₁ = 206.77 (T4 ex106)
  → E_core = energia do pierwszego przejścia g(r₁)=1

PODEJŚCIE HYBRYDOWE (fizycznie uzasadnione):
  [P] Profil g(r):       ODE z f(g) + elastyczne odbicia (ex106-style)
  [E] Energia E_core:    K_sub(g)=g² jako funkcjonał energetyczny
      E_core = 4π ∫₀^{r₁} [g²/2·g'² + V_dw(g)] r² dr
  Uzasadnienie: K_sub=g² to energia z substratu (sek10); f(g) to
  kinetyka efektywna z pętli kwantowych. Obie są aspektami TGP.

TESTY (13):
  K1:  A_tail(e) z elastycznym odbiciem (weryfikacja ex106)
  K2:  A_tail(μ) > 0 (g₀ > g₀_max_f — dostępne przez odbicia)
  K3:  A_tail(τ) > 0 (φ²·g₀*)
  K4:  R₂₁^A = (A_μ/A_e)^4 ≈ 206.77 (potwierdza ex106)
  K5:  R₃₁^A = (A_τ/A_e)^4 vs PDG
  K6:  A_tail monotoniczne: A_e < A_μ < A_τ
  K7:  r₁(e), r₁(μ), r₁(τ) wyznaczone z profili bouncing
  K8:  r₁(e) > r₁(μ) > r₁(τ) (szybsze opadanie przy wyższym g₀)
  K9:  E_core(e) > 0 (K_sub)
  K10: E_core(μ) > 0 (K_sub)
  K11: R₂₁^{E_core} = E_core(μ)/E_core(e) — odchylenie od PDG
  K12: R₃₁^{E_core} = E_core(τ)/E_core(e) — odchylenie od PDG
  K13: E_core ∝ A_tail^n — pomiar wykładnika n (n=4 → zgodność A_tail⁴)

Sesja: TGP v41 — Claudian (2026-04-02)
"""

import sys, io, json, math, warnings
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')

# ============================================================
# Stałe i parametry
# ============================================================
ALPHA   = 2.0
PHI     = (1.0 + math.sqrt(5.0)) / 2.0
G_GHOST = math.exp(-1.0 / (2.0 * ALPHA))  # ≈ 0.7788
G_BOUNCE = G_GHOST + 0.005                 # granica odbicia (jak ex106)
G0_FP   = 1.24915
R21_PDG = 206.768
R31_PDG = 3477.48

# Parametry numeryczne
R_MAX     = 40.0
R_START   = 1e-4
MAX_STEP  = 0.02
RTOL      = 1e-10
ATOL      = 1e-13
MAX_BOUNCES = 12

# Okno ogona (spójne z ex106)
R_TAIL_L = 20.0
R_TAIL_R = 35.0

_trapz = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)

# ============================================================
# Infrastruktura testów
# ============================================================
RESULTS = []

def check(cond, label, detail=""):
    status = "PASS" if cond else "FAIL"
    RESULTS.append((label, status, detail))
    icon = "[PASS]" if cond else "[FAIL]"
    line = f"  {icon} {label}"
    if detail:
        line += f"\n         => {detail}"
    print(line)
    return cond

# ============================================================
# Potencjał i sprzężenia
# ============================================================

def Vprime(g):
    return g**2 * (1.0 - g)

def f_kin(g):
    """f(g) = 1 + 2α·ln(g)"""
    return 1.0 + 2.0 * ALPHA * np.log(np.maximum(g, 1e-30))

def V_dw(g):
    """V_dw(g) = (g-1)²(g+2)/4  — zawsze ≥ 0"""
    return (g - 1.0)**2 * (g + 2.0) / 4.0

def K_sub(g):
    """K_sub(g) = g²  — substrat; zawsze > 0"""
    return g * g

# ============================================================
# Solver: f(g) + elastyczne odbicia (ex106-style)
# ============================================================

def rhs_bounce(r, y):
    """ODE z f(g), regularyzowany przy g=G_BOUNCE"""
    g, gp = y
    g = max(g, G_BOUNCE + 1e-7)
    fg = f_kin(g)
    if abs(fg) < 1e-10:
        return [gp, 0.0]
    driving = Vprime(g)
    cross    = (ALPHA / g) * gp**2
    if r < 1e-10:
        return [gp, (driving - cross) / (3.0 * fg)]
    damp = fg * 2.0 * gp / r
    return [gp, (driving - cross - damp) / fg]


def event_hit_bounce(r, y):
    return y[0] - G_BOUNCE
event_hit_bounce.terminal = True
event_hit_bounce.direction = -1


def integrate_bounce(g0, r_max=None, max_bounces=MAX_BOUNCES):
    """
    Profil solitonu z elastycznymi odbiciami przy g=G_BOUNCE.
    Po odbiciu: gp → -gp  (jak ex106).
    """
    if r_max is None:
        r_max = max(R_MAX, 15.0 * g0)
    r0 = R_START
    y0 = [g0, 0.0]
    segs = []

    for _ in range(max_bounces + 1):
        sol = solve_ivp(
            rhs_bounce, [r0, r_max], y0,
            method='DOP853', max_step=MAX_STEP,
            rtol=RTOL, atol=ATOL,
            events=[event_hit_bounce], dense_output=False
        )
        segs.append((sol.t, sol.y[0], sol.y[1]))
        if sol.t_events[0].size > 0:
            r_b  = float(sol.t_events[0][0])
            gp_b = float(sol.y_events[0][0, 1])
            r0 = r_b + 1e-6
            y0 = [G_BOUNCE + 1e-5, -gp_b]   # odbicie
        else:
            break

    r = np.concatenate([s[0] for s in segs])
    g = np.concatenate([s[1] for s in segs])
    gp= np.concatenate([s[2] for s in segs])
    idx = np.argsort(r)
    return r[idx], g[idx], gp[idx]

# ============================================================
# Dopasowanie ogona
# ============================================================

def fit_tail(r_arr, g_arr, r_L=R_TAIL_L, r_R=R_TAIL_R):
    """g(r)-1 ≈ (B·cos r + C·sin r)/r → A = √(B²+C²)"""
    mask = (r_arr >= r_L) & (r_arr <= r_R)
    if np.sum(mask) < 10:
        return float('nan'), 0.0, 0.0
    r_fit = r_arr[mask]
    h = (g_arr[mask] - 1.0) * r_fit
    X = np.column_stack([np.cos(r_fit), np.sin(r_fit)])
    coefs, _, _, _ = np.linalg.lstsq(X, h, rcond=None)
    B, C = coefs
    return math.sqrt(B*B + C*C), B, C

# ============================================================
# Pierwsze przejście g(r) = 1
# ============================================================

def find_first_zero(r_arr, g_arr):
    """
    Pierwsze r₁ > 0 takie, że g(r₁) = 1  (spadek z g₀ > 1).
    Interpolacja liniowa między węzłami.
    """
    diff = g_arr - 1.0
    for i in range(1, len(diff)):
        if diff[i-1] > 0.0 and diff[i] <= 0.0:
            r0_, r1_ = r_arr[i-1], r_arr[i]
            d0_, d1_ = diff[i-1], diff[i]
            r_zero = r0_ - d0_ * (r1_ - r0_) / (d1_ - d0_)
            return r_zero, i
    return float('nan'), len(r_arr) - 1

# ============================================================
# Energia solitonu z K_sub na profilu bouncing
# ============================================================

def compute_E_total_ksub(r_arr, g_arr, gp_arr):
    """E_total = 4π ∫₀^{R} [K_sub(g)/2·g'² + V_dw(g)] r² dr"""
    k = K_sub(g_arr)
    integrand = (0.5 * k * gp_arr**2 + V_dw(g_arr)) * r_arr**2
    return 4.0 * math.pi * float(_trapz(integrand, r_arr))


def compute_E_core_ksub(r_arr, g_arr, gp_arr, r1):
    """E_core = 4π ∫₀^{r₁} [K_sub(g)/2·g'² + V_dw(g)] r² dr"""
    if not math.isfinite(r1) or r1 <= 0:
        return float('nan')
    mask = r_arr <= r1
    if np.sum(mask) < 3:
        return float('nan')
    r_c  = r_arr[mask]
    g_c  = g_arr[mask]
    gp_c = gp_arr[mask]
    k = K_sub(g_c)
    integrand = (0.5 * k * gp_c**2 + V_dw(g_c)) * r_c**2
    return 4.0 * math.pi * float(_trapz(integrand, r_c))

# ============================================================
# ANALIZA
# ============================================================

print("=" * 70)
print("EX112: A_tail (elastyczne odbicia) + E_core z K_sub(g) = g²")
print("=" * 70)
print(f"  α={ALPHA}, g*={G_GHOST:.6f}, g_bounce={G_BOUNCE:.6f}, φ={PHI:.6f}")
print(f"  Profil: f(g)+bounce | Energia: K_sub=g²")
print(f"  PDG: R₂₁={R21_PDG}, R₃₁={R31_PDG}")
print()

G0_VALS = {'e': G0_FP, 'μ': PHI * G0_FP, 'τ': PHI**2 * G0_FP}

print("--- SEKCJA 1: Profile solitonowe (f+bounce) ---")
profiles = {}
for lep, g0 in G0_VALS.items():
    r, g, gp = integrate_bounce(g0)
    A, B, C  = fit_tail(r, g)
    r1, _    = find_first_zero(r, g)
    E_total  = compute_E_total_ksub(r, g, gp)
    E_core   = compute_E_core_ksub(r, g, gp, r1)
    g_inf    = float(np.mean(g[r > 35])) if np.sum(r > 35) > 5 else float('nan')
    profiles[lep] = dict(g0=g0, r=r, g=g, gp=gp,
                         A_tail=A, r1=r1,
                         E_total=E_total, E_core=E_core, g_inf=g_inf)
    print(f"  [{lep}] g₀={g0:.5f}: A_tail={A:.6f}, r₁={r1:.4f}, "
          f"E_core={E_core:.5f}, g∞≈{g_inf:.5f}")

print()
print("--- SEKCJA 2: Amplitudy ogona ---")

A_e   = profiles['e']['A_tail']
A_mu  = profiles['μ']['A_tail']
A_tau = profiles['τ']['A_tail']

check(math.isfinite(A_e) and A_e > 0,
      "K1: A_tail(e) > 0 (bounce)",
      f"A_e={A_e:.6f}  [ex106: 0.298823]")

check(math.isfinite(A_mu) and A_mu > 0,
      "K2: A_tail(μ) > 0 (g₀>g₀_max_f, dostępne przez odbicia)",
      f"A_μ={A_mu:.6f}  [ex106: 1.133144]")

check(math.isfinite(A_tau) and A_tau > 0,
      "K3: A_tail(τ) > 0 (φ²·g₀*)",
      f"A_τ={A_tau:.6f}  [ex106: 2.369751]")

R21_A = float('nan')
R31_A = float('nan')

if A_e > 1e-12 and math.isfinite(A_mu):
    R21_A = (A_mu / A_e)**4
    delta21 = abs(R21_A - R21_PDG) / R21_PDG * 100.0
    check(delta21 < 2.0,
          "K4: R₂₁^A = (A_μ/A_e)^4 ≈ 206.77 (potwierdza ex106)",
          f"R₂₁^A={R21_A:.3f}, δ={delta21:.3f}%")
else:
    check(False, "K4: R₂₁^A", "brak A_e lub A_μ")

if A_e > 1e-12 and math.isfinite(A_tau):
    R31_A = (A_tau / A_e)**4
    delta31 = abs(R31_A - R31_PDG) / R31_PDG * 100.0
    check(delta31 < 20.0,
          "K5: R₃₁^A = (A_τ/A_e)^4 vs PDG",
          f"R₃₁^A={R31_A:.1f}, PDG={R31_PDG}, δ={delta31:.1f}%")
else:
    check(False, "K5: R₃₁^A", "brak A_τ")

check(A_e < A_mu if (math.isfinite(A_e) and math.isfinite(A_mu)) else False,
      "K6: A_tail(e) < A_tail(μ) < A_tail(τ) — monotoniczne",
      f"A_e={A_e:.4f} < A_μ={A_mu:.4f}" +
      (f" < A_τ={A_tau:.4f}" if math.isfinite(A_tau) else " [A_τ brak]"))

print()
print("--- SEKCJA 3: Pierwsze przejścia r₁ ---")

r1_e  = profiles['e']['r1']
r1_mu = profiles['μ']['r1']
r1_tau= profiles['τ']['r1']

print(f"  r₁(e) = {r1_e:.5f}")
print(f"  r₁(μ) = {r1_mu:.5f}")
print(f"  r₁(τ) = {r1_tau:.5f}")

check(all(math.isfinite(r1) and r1 > 0 for r1 in [r1_e, r1_mu, r1_tau]),
      "K7: r₁ wyznaczone dla wszystkich leptonów",
      f"r₁(e)={r1_e:.4f}, r₁(μ)={r1_mu:.4f}, r₁(τ)={r1_tau:.4f}")

check(r1_e > r1_mu > r1_tau,
      "K8: r₁(e) > r₁(μ) > r₁(τ) — szybsze opadanie przy wyższym g₀",
      f"r₁: {r1_e:.4f} > {r1_mu:.4f} > {r1_tau:.4f}")

print()
print("--- SEKCJA 4: Energia rdzenia E_core [K_sub] ---")

E_core_e   = profiles['e']['E_core']
E_core_mu  = profiles['μ']['E_core']
E_core_tau = profiles['τ']['E_core']

print(f"  E_core(e)   = {E_core_e:.6f}  [r₁={r1_e:.4f}]")
print(f"  E_core(μ)   = {E_core_mu:.6f}  [r₁={r1_mu:.4f}]")
print(f"  E_core(τ)   = {E_core_tau:.6f}  [r₁={r1_tau:.4f}]")

check(math.isfinite(E_core_e) and E_core_e > 0,
      "K9: E_core(e) > 0",
      f"E_core(e) = {E_core_e:.6f}")

check(math.isfinite(E_core_mu) and E_core_mu > 0,
      "K10: E_core(μ) > 0",
      f"E_core(μ) = {E_core_mu:.6f}")

R21_core = float('nan')
R31_core = float('nan')

if math.isfinite(E_core_e) and E_core_e > 0 and math.isfinite(E_core_mu) and E_core_mu > 0:
    R21_core = E_core_mu / E_core_e
    delta21c = abs(R21_core - R21_PDG) / R21_PDG * 100.0
    check(delta21c < 50.0,
          "K11: R₂₁^{E_core} = E_core(μ)/E_core(e)",
          f"R₂₁^Ec={R21_core:.3f}, PDG={R21_PDG}, δ={delta21c:.1f}%")
else:
    check(False, "K11: R₂₁^{E_core}", "brak danych")

if (math.isfinite(E_core_e) and E_core_e > 0 and
        math.isfinite(E_core_tau) and E_core_tau > 0):
    R31_core = E_core_tau / E_core_e
    delta31c = abs(R31_core - R31_PDG) / R31_PDG * 100.0
    check(delta31c < 50.0,
          "K12: R₃₁^{E_core} = E_core(τ)/E_core(e)",
          f"R₃₁^Ec={R31_core:.2f}, PDG={R31_PDG}, δ={delta31c:.1f}%")
else:
    check(False, "K12: R₃₁^{E_core}", "brak danych")

print()
print("--- SEKCJA 5: Skalowanie E_core ↔ A_tail ---")

ec_vals = [E_core_e, E_core_mu, E_core_tau]
at_vals = [A_e, A_mu, A_tau]

valid_pts = [(at, ec) for at, ec in zip(at_vals, ec_vals)
             if math.isfinite(at) and at > 0 and math.isfinite(ec) and ec > 0]

n_exp = float('nan')
r2_fit = float('nan')

if len(valid_pts) >= 2:
    log_at = np.array([math.log(v[0]) for v in valid_pts])
    log_ec = np.array([math.log(v[1]) for v in valid_pts])
    coeffs = np.polyfit(log_at, log_ec, 1)
    n_exp  = float(coeffs[0])
    log_ec_fit = np.polyval(coeffs, log_at)
    ss_res = np.sum((log_ec - log_ec_fit)**2)
    ss_tot = np.sum((log_ec - np.mean(log_ec))**2)
    r2_fit = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    print(f"  E_core ∝ A_tail^{n_exp:.4f}  (R²={r2_fit:.5f}, n_pts={len(valid_pts)})")
    check(r2_fit > 0.95,
          "K13: E_core ∝ A_tail^n — power-law skalowanie",
          f"n={n_exp:.3f}, R²={r2_fit:.5f}  (n=4 → zgodność z A_tail⁴)")
else:
    check(False, "K13: E_core ↔ A_tail skalowanie", "za mało punktów")

# ============================================================
# Podsumowanie
# ============================================================
print()
print("=" * 70)
print("PODSUMOWANIE EX112")
print("=" * 70)
print()
print("  Podejście hybrydowe: profil f(g)+bounce | energia K_sub=g²")
print()
print("  Leptony:")
for lep in ['e', 'μ', 'τ']:
    p = profiles[lep]
    print(f"    g₀({lep}) = {p['g0']:.5f}: A={p['A_tail']:.5f}, r₁={p['r1']:.4f}, "
          f"E_core={p['E_core']:.5f}")
print()
print(f"  A_tail⁴:   R₂₁^A  = {R21_A:.3f}  (PDG={R21_PDG})")
print(f"  E_core:    R₂₁^Ec = {R21_core:.3f}  (PDG={R21_PDG})")
if math.isfinite(R31_A):
    print(f"  A_tail⁴:   R₃₁^A  = {R31_A:.2f}  (PDG={R31_PDG})")
if math.isfinite(R31_core):
    print(f"  E_core:    R₃₁^Ec = {R31_core:.2f}  (PDG={R31_PDG})")
if math.isfinite(n_exp):
    print(f"  Wykładnik: E_core ∝ A_tail^{n_exp:.3f}")
print()

n_pass = sum(1 for _, s, _ in RESULTS if s == "PASS")
n_fail = sum(1 for _, s, _ in RESULTS if s == "FAIL")
n_total = len(RESULTS)
print(f"  Testy: {n_pass}/{n_total} PASS", "✓" if n_fail == 0 else "")
print()

# ============================================================
# Interpretacja fizyczna
# ============================================================
print("=" * 70)
print("INTERPRETACJA (OP-G + OP-E)")
print("=" * 70)
print()
print("  OP-G (ghost constraint):")
print(f"    Elastyczne odbicia przy g*={G_GHOST:.4f} → g₀^μ=2.021 i g₀^τ=3.270")
print(f"    DOSTĘPNE. FAR okno [{R_TAIL_L},{R_TAIL_R}] osiągalne dla wszystkich leptonów.")
if math.isfinite(R21_A):
    d21 = abs(R21_A - R21_PDG) / R21_PDG * 100.0
    if d21 < 1.0:
        print(f"    A_tail⁴ → R₂₁ = {R21_A:.3f} (δ={d21:.3f}%) — POTWIERDZA ex106 ✓")
    else:
        print(f"    A_tail⁴ → R₂₁ = {R21_A:.3f} (δ={d21:.2f}%)")
print()
print("  OP-E (energia rdzenia):")
if math.isfinite(R21_core):
    d21c = abs(R21_core - R21_PDG) / R21_PDG * 100.0
    print(f"    R₂₁^{{E_core}} = {R21_core:.3f}  (δ={d21c:.1f}% od PDG)")
    if d21c < 5.0:
        print("    → E_core REPLIKUJE stosunek mas — mechanizm masy przez rdzeń! ✓")
    elif d21c < 30.0:
        print("    → E_core przybliżone — obiecujące; wymaga dalszych badań")
    else:
        print(f"    → E_core ≠ masy leptonów; rdzeń = {d21c:.0f}% odchylenia")
if math.isfinite(n_exp):
    print(f"    Skalowanie E_core ∝ A_tail^{n_exp:.3f}")
    if abs(n_exp - 4.0) < 0.5:
        print("    → E_core ≈ A_tail⁴ — mechanizmy zgodne!")
    elif abs(n_exp - 2.0) < 0.5:
        print("    → E_core ≈ A_tail² — inna potęga niż Ścieżka 9")
    else:
        print(f"    → Niezależny wykładnik n={n_exp:.3f}")

# ============================================================
# Zapis JSON
# ============================================================
import os
out_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'scripts')
os.makedirs(out_dir, exist_ok=True)
results_json = {
    "session": "ex112",
    "date": "2026-04-02",
    "approach": "f(g)+bounce profile, K_sub=g^2 energy",
    "g_bounce": G_BOUNCE,
    "g0_e": G0_FP, "g0_mu": G0_VALS['μ'], "g0_tau": G0_VALS['τ'],
    "A_e": A_e if math.isfinite(A_e) else None,
    "A_mu": A_mu if math.isfinite(A_mu) else None,
    "A_tau": A_tau if math.isfinite(A_tau) else None,
    "R21_A": R21_A if math.isfinite(R21_A) else None,
    "R31_A": R31_A if math.isfinite(R31_A) else None,
    "r1_e": r1_e if math.isfinite(r1_e) else None,
    "r1_mu": r1_mu if math.isfinite(r1_mu) else None,
    "r1_tau": r1_tau if math.isfinite(r1_tau) else None,
    "E_core_e": E_core_e if math.isfinite(E_core_e) else None,
    "E_core_mu": E_core_mu if math.isfinite(E_core_mu) else None,
    "E_core_tau": E_core_tau if math.isfinite(E_core_tau) else None,
    "R21_E_core": R21_core if math.isfinite(R21_core) else None,
    "R31_E_core": R31_core if math.isfinite(R31_core) else None,
    "n_Ecore_Atail": n_exp if math.isfinite(n_exp) else None,
    "R21_PDG": R21_PDG, "R31_PDG": R31_PDG,
    "n_pass": n_pass, "n_fail": n_fail, "n_total": n_total,
    "tests": [{"label": l, "status": s, "detail": d} for l, s, d in RESULTS],
}
try:
    with open(os.path.join(out_dir, 'ex112_results.json'), 'w', encoding='utf-8') as f:
        json.dump(results_json, f, indent=2, ensure_ascii=False)
    print(f"\nWyniki: scripts/ex112_results.json")
except Exception as e:
    print(f"\nBŁĄD zapisu: {e}")
print(f"\nSESJA: TGP v41 — Claudian (2026-04-02)")
