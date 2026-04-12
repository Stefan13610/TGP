#!/usr/bin/env python3
"""
ex97_amplitude_convergence.py — Test zbieżności A_tail(φ²·z₀) z R_MAX
========================================================================
STATUS: LEGACY-TRANSLATIONAL

This file is part of the older tau/amplitude diagnostic chain around `G3`,
`alpha*`, and legacy matching rules. Treat it as historical exploratory
material, not as a canonical current `nbody` example.

Pytanie kluczowe po ex95/ex96:
  Z ex95: per-okno A_tau plateau ≈ 2.994 dla r∈[60,136] (R_MAX=150)
  Z ex96: G₃_far = (2.994/0.347)^4 ≈ 5540 >> R31=3477 WSZĘDZIE na [2.5,5.0]

  Ale czy plateau r∈[60,136] to prawdziwa A_∞?
  Czy soliton τ z g₀=φ²·z₀≈3.31 nadal ewoluuje dla r>136?

Plan:
  1. α=2.929 (ex93 α*₃): pełny profil A per-okno od r=20 do r=280 (R_MAX=300)
  2. α=3.500 (wyższy): to samo — kontrola zachowania dla dużych α
  3. α=2.500 (minimum G₃_far): to samo
  4. Porównanie Ae, Amu, Atau vs r dla wszystkich trzech α
  5. Dopasowanie A(r) = A_∞·(1 + a/r + b/r²) do całego zakresu
     → ekstrapolacja A_∞ z dużych r

Autor: Claudian (sesja v38, O-L5 diagnostyka)
Data: 2026-03-28
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import time
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE
# ================================================================
PHI  = (1.0 + np.sqrt(5.0)) / 2.0
PHI2 = PHI * PHI

M_TAU_MEV = 1776.86
M_MU_MEV  = 105.658
M_E_MEV   = 0.51100

R21 = M_MU_MEV / M_E_MEV
R31 = M_TAU_MEV / M_E_MEV

G_OFF   = 0.005
B_WIN_L = 28.0
B_WIN_R = 42.0
R_MAX   = 300   # Duże R_MAX dla testu zbieżności

# Okna od r=20 do r=280 (co 10, szerokie na 16)
PROBE_WINDOWS = [(r, r+16) for r in range(20, 271, 10)]

# ================================================================
# ODE — identyczna jak ex95/ex96
# ================================================================
def integrate_soliton(g0, alpha, r_max):
    gb = np.exp(-1.0 / (2.0 * alpha)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-8)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        Vprime = g * g * (1.0 - g)
        curl   = (alpha / g) * gp * gp
        if r < 1e-10: return [gp, (Vprime - curl) / (3.0 * fg)]
        return [gp, (Vprime - curl - fg * 2.0 * gp / r) / fg]
    def ev_bounce(r, y): return y[0] - gb
    ev_bounce.terminal = True; ev_bounce.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(30):  # więcej iteracji dla dużego R_MAX
        sol = solve_ivp(rhs, (r0, r_max), y0, events=ev_bounce,
                        dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.04)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else:
            break
    idx = np.argsort(np.concatenate(ra))
    return np.concatenate(ra)[idx], np.concatenate(ga)[idx]

def amplitude_per_window(r, g, wins=None):
    """Oblicz A per-okno; domyślnie PROBE_WINDOWS."""
    if wins is None:
        wins = PROBE_WINDOWS
    results = []
    for rL, rR in wins:
        if rR > r[-1] - 1: continue
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        A = float(np.sqrt(bc[0]**2 + bc[1]**2))
        results.append((float(rL), float(rR), A))
    return results

def B_coeff_func(g0, alpha, r_max):
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max=100, lo=1.05, hi=2.0, n=60):
    gs = np.linspace(lo, hi, n)
    Bs = []
    for g in gs:
        try: Bs.append(B_coeff_func(g, alpha, r_max))
        except: Bs.append(np.nan)
    Bs = np.array(Bs)
    for i in range(len(Bs) - 1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i] * Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff_func(x, alpha, r_max),
                              gs[i], gs[i+1], xtol=1e-9, maxiter=100)
            except: pass
    return None

def fit_Ainf(rvals, Avals):
    """Ekstrapoluj A_∞ z modelu A_∞(1+a/r+b/r²), zwróć A_∞ i opis."""
    rvals = np.array(rvals); Avals = np.array(Avals)
    if len(rvals) < 3:
        return float(np.median(Avals)), 0.0, 'median (n<3)'
    # Cv test
    cv = np.std(Avals) / (np.mean(Avals) + 1e-12)
    if cv < 0.003:
        return float(np.median(Avals)), float(np.std(Avals)), 'plateau'
    # Model fit
    try:
        p, pcov = curve_fit(lambda x, ai, a, b: ai * (1 + a/x + b/x**2),
                            rvals, Avals, p0=[Avals[-1], 0., 0.], maxfev=3000)
        err = float(np.sqrt(abs(np.diag(pcov)[0])))
        return float(p[0]), err, f'fit 1+a/r+b/r² (a={p[1]:.3f}, b={p[2]:.1f})'
    except:
        return float(np.median(Avals)), float(np.std(Avals)), 'median (fit failed)'

def analyze_alpha(alpha, label=''):
    """Kompletna analiza zbieżności A_tail dla danego α."""
    print(f"\n{'='*68}")
    print(f"  α = {alpha:.6f}  {label}")
    print(f"{'='*68}")

    z0 = find_z0(alpha)
    if z0 is None:
        print(f"  ⚠ z₀ nie znaleziona!")
        return None
    print(f"  z₀={z0:.6f}, φ·z₀={PHI*z0:.6f}, φ²·z₀={PHI2*z0:.6f}")

    # Integracje
    print(f"  Integruję solitony do R_MAX={R_MAX}...")
    t0 = time.time()
    r0, g0_sol = integrate_soliton(z0,      alpha, R_MAX)
    r1, g1_sol = integrate_soliton(PHI*z0,  alpha, R_MAX)
    r2, g2_sol = integrate_soliton(PHI2*z0, alpha, R_MAX)
    print(f"  Czas integracji: {time.time()-t0:.1f}s")
    print(f"  r zakres: [0, {r0[-1]:.1f}], [0, {r1[-1]:.1f}], [0, {r2[-1]:.1f}]")

    # Profil per okno
    wins_e  = amplitude_per_window(r0, g0_sol)
    wins_mu = amplitude_per_window(r1, g1_sol)
    wins_tau= amplitude_per_window(r2, g2_sol)

    print(f"\n  Profil A per-okno:")
    print(f"  {'rL':>6} {'rR':>6} | {'Ae':>10} {'Amu':>10} {'Atau':>10} | {'F':>8} {'G3':>10}")
    print(f"  {'-'*70}")

    r_profile = []; Ae_p=[]; Amu_p=[]; Atau_p=[]
    for (rLe, rRe, Ae_w), (rLm, _, Amu_w), (rLt, _, Atau_w) in zip(wins_e, wins_mu, wins_tau):
        F_w  = (Amu_w / Ae_w)**4  if Ae_w > 1e-8 else float('nan')
        G3_w = (Atau_w / Ae_w)**4 if Ae_w > 1e-8 else float('nan')
        print(f"  [{rLe:3.0f},{rRe:3.0f}] | {Ae_w:10.6f} {Amu_w:10.6f} {Atau_w:10.6f} | "
              f"{F_w:8.4f} {G3_w:10.4f}")
        r_profile.append(rLe); Ae_p.append(Ae_w); Amu_p.append(Amu_w); Atau_p.append(Atau_w)

    # Ekstrapolacja A_∞ dla trzech podziałów zakresu
    print(f"\n  Ekstrapolacja A_∞ z różnych zakresów r:")
    ranges = [
        ('all [20,286]', None),
        ('std [20,66]',  [(r,A) for r,A in zip(r_profile, Atau_p) if r <= 66]),
        ('mid [60,136]', [(r,A) for r,A in zip(r_profile, Atau_p) if 60 <= r <= 136]),
        ('far [120,220]',[(r,A) for r,A in zip(r_profile, Atau_p) if 120 <= r <= 220]),
        ('very far [200,286]', [(r,A) for r,A in zip(r_profile, Atau_p) if r >= 200]),
    ]

    results_Ainf = {}
    for name, subset in ranges:
        if subset is None:
            rv = r_profile; Av = Atau_p
        else:
            if not subset: continue
            rv, Av = zip(*subset)
            rv, Av = list(rv), list(Av)
        if len(rv) < 2: continue
        Ainf, err, desc = fit_Ainf(rv, Av)
        G3 = (Ainf / np.median(Ae_p))**4 if np.median(Ae_p) > 1e-8 else float('nan')
        print(f"    {name:25s}: A_tau_∞={Ainf:.6f} ±{err:.6f}  "
              f"→ G₃={G3:.2f} (G₃−R₃₁={G3-R31:+.2f})  [{desc}]")
        results_Ainf[name] = {'Ainf': Ainf, 'G3': G3}

    # Trend A_tau (ostatnie 50% okien)
    n = len(Atau_p)
    if n >= 4:
        half = n // 2
        rr = np.array(r_profile[half:])
        Arr = np.array(Atau_p[half:])
        slope = np.polyfit(rr, Arr, 1)[0]
        print(f"\n  Trend A_tau w drugiej połowie okien: dA/dr = {slope:+.6f} / unit_r")
        if abs(slope) < 1e-4:
            print(f"  → Plateau potwierdzone (|slope| < 1e-4)")
        elif slope < 0:
            print(f"  → A_tau maleje! (nadal transient)")
        else:
            print(f"  → A_tau rośnie (niestandardowe)")

    return results_Ainf

# ================================================================
# MAIN
# ================================================================
def main():
    t0_total = time.time()
    print("=" * 68)
    print("ex97_amplitude_convergence.py — Test zbieżności A_tau")
    print("=" * 68)
    print(f"  R_MAX = {R_MAX}")
    print(f"  R21={R21:.4f}, R31={R31:.4f}")
    print(f"  PROBE_WINDOWS: od [{PROBE_WINDOWS[0][0]},{PROBE_WINDOWS[0][1]}]"
          f" do [{PROBE_WINDOWS[-1][0]},{PROBE_WINDOWS[-1][1]}]")
    print(f"  {len(PROBE_WINDOWS)} okien, szerokość 16 każde")
    print()

    # Trzy testowe wartości α
    test_alphas = [
        (2.929524, "← ex93 α*₃ (artefakt STD)"),
        (2.500,    "← minimum G₃_far w ex96"),
        (3.500,    "← wysoka α, G₃_far duże"),
    ]

    all_results = {}
    for alpha, label in test_alphas:
        res = analyze_alpha(alpha, label)
        if res:
            all_results[alpha] = res

    # Podsumowanie
    print()
    print("=" * 68)
    print("PODSUMOWANIE ex97 — Zbieżność A_tail")
    print("=" * 68)
    print()
    print(f"  {'α':>8} | {'A_tau_far':>12} | {'G₃_far':>10} | {'G₃−R₃₁':>12} | trend")
    print(f"  {'-'*65}")
    for alpha, label in test_alphas:
        if alpha in all_results:
            r = all_results[alpha]
            # Preferuj 'far [120,220]' lub 'mid [60,136]'
            key = 'far [120,220]' if 'far [120,220]' in r else 'mid [60,136]'
            if key in r:
                Ainf = r[key]['Ainf']
                G3   = r[key]['G3']
                print(f"  {alpha:8.3f} | {Ainf:12.6f} | {G3:10.2f} | {G3-R31:+12.2f} | {key}")

    print()
    print(f"  R₃₁ = {R31:.4f}  (m_τ/m_e, nie czwarta potęga)")
    print(f"\nCzas całkowity: {time.time()-t0_total:.1f} s")
    print("=" * 68)

if __name__ == '__main__':
    main()
