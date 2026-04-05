#!/usr/bin/env python3
"""
ex93_tau_alpha_star.py — Precyzyjne α*₃ dla sektora tau (O-L5, v37)
====================================================================
Cel:
  Wyznaczenie α*₃ gdzie G₃(α*₃) = R31 = m_τ/m_e ≈ 3477.22

  G₃(α) ≡ (A_tail(φ²·z₀(α), α) / A_tail(z₀(α), α))^4

Z ex92: kandydat α*₃ ≈ 2.928 (G₃: 4100→3285 między α=2.815 a α=2.983).

Strategia:
  1. Gęsty skan G₃(α) na α∈[2.70, 3.30] z krokiem 0.01 przy R_MAX=120
  2. Identyfikacja przedziałów dla brentq
  3. Precyzyjne zera brentq (xtol=1e-8)
  4. Weryfikacja stabilności: R_MAX ∈ {120, 150, 200} dla α*₃
  5. Dodatkowa diagnostyka: G₃/F ratio, R32=m_τ/m_μ check

Parametry stałe:
  B_WIN = [28.0, 42.0]  (identyczne z ex88, ex92)
  FIT_WINS = [(20,36),(30,46),(40,56),(50,66)]
  PHI  = 1.6180...
  PHI2 = PHI² ≈ 2.6180

Autor: Claudian (sesja v37, O-L5)
Data: 2026-03-28
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
import multiprocessing as mp
import time
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE FIZYCZNE
# ================================================================
PHI  = (1.0 + np.sqrt(5.0)) / 2.0
PHI2 = PHI * PHI

M_TAU_MEV = 1776.86
M_MU_MEV  = 105.658
M_E_MEV   = 0.51100

R21  = M_MU_MEV / M_E_MEV    # ≈ 206.768
R31  = M_TAU_MEV / M_E_MEV   # ≈ 3477.221
R32  = M_TAU_MEV / M_MU_MEV  # ≈ 16.817

G_OFF    = 0.005
B_WIN_L  = 28.0
B_WIN_R  = 42.0
FIT_WINS = [(20, 36), (30, 46), (40, 56), (50, 66)]

# Zakres skanu
ALPHA_SCAN_LO = 2.70
ALPHA_SCAN_HI = 3.30
ALPHA_STEP    = 0.010   # gęsty skan

# R_MAX dla weryfikacji stabilności
RMAX_SCAN    = 120
RMAX_VERIFY  = [120, 150, 200]

print(f"STAŁE: PHI={PHI:.6f}, R21={R21:.4f}, R31={R31:.4f}, R32={R32:.4f}")
print()

# ================================================================
# ODE + AMPLITUDY
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
    for _ in range(20):   # więcej odbić dla większych g0
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

def fit_amplitude(r, g):
    Av, rv = [], []
    for rL, rR in FIT_WINS:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2))); rv.append(float(rL))
    if not Av: return 0.0, 0.0
    if len(Av) < 2: return Av[-1], 0.0
    rv = np.array(rv); Av_arr = np.array(Av)
    try:
        p, pcov = curve_fit(lambda x, ai, a, b: ai * (1 + a/x + b/x**2),
                            rv, Av_arr, p0=[Av_arr[-1], 0., 0.], maxfev=2000)
        return float(p[0]), float(np.sqrt(abs(np.diag(pcov)[0])))
    except:
        try:
            p, _ = curve_fit(lambda x, ai, a: ai * (1 + a/x),
                             rv, Av_arr, p0=[Av_arr[-1], 0.], maxfev=1000)
            return float(p[0]), 0.0
        except:
            return float(Av_arr[-1]), 0.0

def B_coeff(g0, alpha, r_max):
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max, lo=1.05, hi=2.0, n=60):
    """z₀(α) na gałęzi fizycznej (z₀≈1.1–1.3)."""
    gs = np.linspace(lo, hi, n)
    Bs = []
    for g in gs:
        try: Bs.append(B_coeff(g, alpha, r_max))
        except: Bs.append(np.nan)
    Bs = np.array(Bs)
    for i in range(len(Bs) - 1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i] * Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff(x, alpha, r_max),
                              gs[i], gs[i+1], xtol=1e-9, maxiter=100)
            except: pass
    return None

def compute_G3_F(alpha, r_max):
    """
    Oblicz G₃(α) i F(α) dla danego α i r_max.
    Zwraca słownik z pełnymi wynikami.
    """
    z0 = find_z0(alpha, r_max)
    if z0 is None:
        return {'alpha': alpha, 'r_max': r_max, 'z0': None,
                'F': None, 'G3': None, 'Ae': None, 'Amu': None, 'Atau': None}

    # Elektron
    r, g = integrate_soliton(z0, alpha, r_max)
    Ae, dAe = fit_amplitude(r, g)

    # Muon (φ·z₀)
    z0_mu = PHI * z0
    r2, g2 = integrate_soliton(z0_mu, alpha, r_max)
    Amu, dAmu = fit_amplitude(r2, g2)

    # Tau (φ²·z₀)
    z0_tau = PHI2 * z0
    r3, g3 = integrate_soliton(z0_tau, alpha, r_max)
    Atau, dAtau = fit_amplitude(r3, g3)

    if Ae < 1e-8:
        return {'alpha': alpha, 'r_max': r_max, 'z0': z0,
                'F': None, 'G3': None, 'Ae': Ae, 'Amu': Amu, 'Atau': Atau}

    F  = (Amu  / Ae) ** 4
    G3 = (Atau / Ae) ** 4

    return {
        'alpha': alpha, 'r_max': r_max, 'z0': z0,
        'F': F, 'G3': G3,
        'Ae': Ae, 'Amu': Amu, 'Atau': Atau,
        'dAe': dAe, 'dAmu': dAmu, 'dAtau': dAtau,
        'z0_mu': z0_mu, 'z0_tau': z0_tau
    }

def worker_scan(args):
    alpha, r_max = args
    return compute_G3_F(alpha, r_max)

# ================================================================
# PRECYZYJNE ZERO brentq
# ================================================================
def G3_minus_R31(alpha, r_max):
    """G₃(α) - R31 dla brentq."""
    res = compute_G3_F(alpha, r_max)
    if res['G3'] is None:
        return float('nan')
    return res['G3'] - R31

def precise_zero(a_lo, a_hi, r_max, label=""):
    """brentq dla G₃=R31 na [a_lo, a_hi]."""
    try:
        f_lo = G3_minus_R31(a_lo, r_max)
        f_hi = G3_minus_R31(a_hi, r_max)
        if np.isnan(f_lo) or np.isnan(f_hi):
            print(f"  [{label}] NaN w przedziale [{a_lo:.4f},{a_hi:.4f}]")
            return None
        if f_lo * f_hi > 0:
            print(f"  [{label}] Brak zmiany znaku: G3-R31({a_lo:.4f})={f_lo:.2f}, ({a_hi:.4f})={f_hi:.2f}")
            return None
        alpha_star = brentq(lambda a: G3_minus_R31(a, r_max),
                            a_lo, a_hi, xtol=1e-8, maxiter=200)
        return alpha_star
    except Exception as e:
        print(f"  [{label}] brentq błąd: {e}")
        return None

# ================================================================
# MAIN
# ================================================================
def main():
    t0 = time.time()
    print("=" * 68)
    print("ex93_tau_alpha_star.py — Precyzyjne α*₃ dla tauonu (O-L5)")
    print("=" * 68)
    print(f"  Skan: α∈[{ALPHA_SCAN_LO},{ALPHA_SCAN_HI}], krok={ALPHA_STEP}")
    print(f"  R_MAX scan = {RMAX_SCAN}, weryfikacja = {RMAX_VERIFY}")
    print(f"  R31 = {R31:.6f}  (m_τ/m_e)")
    print(f"  R21 = {R21:.6f}  (m_μ/m_e)")
    print(f"  R32 = {R32:.6f}  (m_τ/m_μ)")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Gęsty skan G₃(α) przy R_MAX=120
    # ----------------------------------------------------------------
    print("=" * 68)
    print(f"FAZA 1: Gęsty skan G₃(α) przy R_MAX={RMAX_SCAN}")
    print("=" * 68)

    alpha_arr = np.arange(ALPHA_SCAN_LO, ALPHA_SCAN_HI + ALPHA_STEP/2, ALPHA_STEP)
    args_list = [(a, RMAX_SCAN) for a in alpha_arr]

    ncpu = min(mp.cpu_count(), 6)
    print(f"  {len(alpha_arr)} punktów, {ncpu} rdzeni...")

    with mp.Pool(ncpu) as pool:
        scan_results = pool.map(worker_scan, args_list)

    good = [r for r in scan_results if r['G3'] is not None]
    bad  = [r for r in scan_results if r['G3'] is None]
    print(f"  Wyniki: {len(good)}/{len(scan_results)} poprawnych")
    if bad:
        bad_alphas = [f"{r['alpha']:.3f}" for r in bad[:8]]
        print(f"  Brak z0: alpha w {bad_alphas}")
    print()

    # Tabela wynikow
    alphas_s = np.array([r['alpha'] for r in good])
    G3s      = np.array([r['G3']    for r in good])
    Fs       = np.array([r['F']     for r in good])
    z0s      = np.array([r['z0']    for r in good])

    print(f"  {'α':>7} | {'z₀':>8} | {'F(α)':>9} | {'G₃(α)':>10} | {'G₃-R31':>10} | {'G₃/F':>7}")
    print(f"  {'-'*62}")
    for r in good:
        g3r31 = r['G3'] - R31
        g3f   = r['G3'] / r['F'] if r['F'] and r['F'] > 0 else float('nan')
        marker = " ←" if abs(g3r31) < 200 else ""
        print(f"  {r['alpha']:>7.3f} | {r['z0']:>8.5f} | {r['F']:>9.3f} | {r['G3']:>10.3f} | "
              f"{g3r31:>+10.3f} | {g3f:>7.3f}{marker}")
    print()

    # Znajdź zmiany znaku G₃-R31
    print("=" * 68)
    print("ZMIANY ZNAKU G₃-R31 (kandydaci brentq):")
    print("=" * 68)
    brackets = []
    for i in range(len(good) - 1):
        d1 = good[i]['G3']   - R31
        d2 = good[i+1]['G3'] - R31
        if d1 * d2 < 0:
            a1, a2 = good[i]['alpha'], good[i+1]['alpha']
            brackets.append((a1, a2, 'malejąco' if d1 > 0 else 'rosnąco'))
            print(f"  Bracket [{a1:.3f}, {a2:.3f}]  G₃={good[i]['G3']:.2f}→{good[i+1]['G3']:.2f}  ({brackets[-1][2]})")

    if not brackets:
        print("  Brak zmian znaku w przebadanym zakresie!")
        print(f"  G₃ zakres: [{G3s.min():.1f}, {G3s.max():.1f}]")

    print()

    # ----------------------------------------------------------------
    # FAZA 2: Precyzyjny brentq dla każdego bracketa
    # ----------------------------------------------------------------
    print("=" * 68)
    print(f"FAZA 2: Precyzyjny brentq (xtol=1e-8) przy R_MAX={RMAX_SCAN}")
    print("=" * 68)

    alpha_stars = []
    for a_lo, a_hi, direction in brackets:
        print(f"\n  Szukam w [{a_lo:.3f}, {a_hi:.3f}] ({direction})...")
        astar = precise_zero(a_lo, a_hi, RMAX_SCAN, label=f"R{RMAX_SCAN}")
        if astar is not None:
            # Weryfikacja
            res_check = compute_G3_F(astar, RMAX_SCAN)
            g3_check  = res_check['G3'] if res_check['G3'] else float('nan')
            f_check   = res_check['F']  if res_check['F']  else float('nan')
            print(f"  → α*₃ = {astar:.8f}")
            print(f"    G₃(α*₃) = {g3_check:.6f}  (R31={R31:.4f}, diff={g3_check-R31:+.6f})")
            print(f"    F(α*₃)  = {f_check:.6f}  (R21={R21:.4f}, diff={f_check-R21:+.6f})")
            print(f"    z₀      = {res_check['z0']:.6f}")
            alpha_stars.append({'alpha': astar, 'direction': direction,
                                'G3': g3_check, 'F': f_check,
                                'z0': res_check['z0']})
        else:
            print(f"  → brentq nie znalazł zera w tym przedziale")

    print()

    # ----------------------------------------------------------------
    # FAZA 3: Weryfikacja stabilności α*₃ względem R_MAX
    # ----------------------------------------------------------------
    if alpha_stars:
        print("=" * 68)
        print("FAZA 3: Stabilność α*₃ vs R_MAX")
        print("=" * 68)

        # Bierzemy pierwszy (główny) kandydat
        main_star = alpha_stars[0]
        a_lo_main = main_star['alpha'] - 0.05
        a_hi_main = main_star['alpha'] + 0.05

        print(f"  Kandydat główny α*₃ ≈ {main_star['alpha']:.6f}")
        print(f"  Bracket dla weryfikacji: [{a_lo_main:.4f}, {a_hi_main:.4f}]")
        print()

        print(f"  {'R_MAX':>6} | {'α*₃':>12} | {'G₃(α*₃)':>12} | {'G₃-R31':>10} | {'F(α*₃)':>10}")
        print(f"  {'-'*58}")

        rmax_results = []
        for rmax in RMAX_VERIFY:
            astar_r = precise_zero(a_lo_main, a_hi_main, rmax, label=f"R{rmax}")
            if astar_r is not None:
                res_r = compute_G3_F(astar_r, rmax)
                g3_r  = res_r['G3'] if res_r['G3'] else float('nan')
                f_r   = res_r['F']  if res_r['F']  else float('nan')
                print(f"  {rmax:>6} | {astar_r:>12.8f} | {g3_r:>12.4f} | {g3_r-R31:>+10.4f} | {f_r:>10.4f}")
                rmax_results.append({'rmax': rmax, 'alpha': astar_r, 'G3': g3_r})
            else:
                print(f"  {rmax:>6} | {'BRAK':>12} | {'—':>12} | {'—':>10} | {'—':>10}")

        # Ocena stabilności
        if len(rmax_results) >= 2:
            alphas_rmax = [r['alpha'] for r in rmax_results]
            sigma_alpha = np.std(alphas_rmax)
            spread = max(alphas_rmax) - min(alphas_rmax)
            print()
            print(f"  Spread α*₃ vs R_MAX: {spread:.2e}  (σ={sigma_alpha:.2e})")
            if spread < 0.01:
                print(f"  ✅ STABILNA: α*₃ niezależna od R_MAX w granicach {spread:.2e}")
            elif spread < 0.05:
                print(f"  ⚠️  UMIARKOWANA: α*₃ zmienia się o {spread:.4f} — potrzebna wyższa precyzja")
            else:
                print(f"  ❌ NIESTABILNA: α*₃ silnie zależy od R_MAX (Δ={spread:.4f})")
                print(f"     G₃(φ²·z₀) jest numerycznie niestabilna w tym zakresie")

    print()

    # ----------------------------------------------------------------
    # FAZA 4: Diagnostyka G₃/F — relacja tau/muon
    # ----------------------------------------------------------------
    print("=" * 68)
    print("FAZA 4: Diagnostyka G₃/F = (A_tau/A_mu)^4 vs R32=m_τ/m_μ")
    print("=" * 68)
    print(f"  R32 = m_τ/m_μ = {R32:.6f}")
    print(f"  φ^4 = {PHI**4:.6f},  φ^6 = {PHI**6:.6f},  φ^8 = {PHI**8:.6f}")
    print()
    print(f"  Szukam α gdzie G₃/F = R32 = {R32:.4f}:")
    print()
    print(f"  {'α':>7} | {'G₃(α)':>10} | {'F(α)':>9} | {'G₃/F':>8} | {'vs R32':>9}")
    print(f"  {'-'*55}")
    for r in good:
        if r['F'] and r['F'] > 0:
            g3f = r['G3'] / r['F']
            vs_r32 = (g3f - R32) / R32 * 100
            marker = " ←" if abs(vs_r32) < 10 else ""
            print(f"  {r['alpha']:>7.3f} | {r['G3']:>10.3f} | {r['F']:>9.3f} | {g3f:>8.4f} | {vs_r32:>+8.2f}%{marker}")

    # Zero G₃/F = R32?
    g3f_arr = np.array([r['G3']/r['F'] for r in good if r['F'] and r['F'] > 0])
    a_arr   = np.array([r['alpha']      for r in good if r['F'] and r['F'] > 0])
    print()
    print(f"  G₃/F zakres: [{g3f_arr.min():.4f}, {g3f_arr.max():.4f}]")
    if g3f_arr.min() <= R32 <= g3f_arr.max():
        print(f"  R32={R32:.4f} leży w zakresie G₃/F — zero możliwe!")
    else:
        print(f"  R32={R32:.4f} POZA zakresem G₃/F w α∈[{ALPHA_SCAN_LO},{ALPHA_SCAN_HI}]")

    print()

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    print("=" * 68)
    print("PODSUMOWANIE ex93 — O-L5")
    print("=" * 68)

    if alpha_stars:
        print(f"\n  Znalezione zera G₃=R31 (brentq, R_MAX={RMAX_SCAN}):")
        for i, s in enumerate(alpha_stars):
            print(f"  [{i+1}] α*₃ = {s['alpha']:.8f}  ({s['direction']})")
            print(f"       G₃  = {s['G3']:.6f}  (R31={R31:.4f})")
            print(f"       F   = {s['F']:.6f}  (R21={R21:.4f})")
            print(f"       z₀  = {s['z0']:.6f}")

        main = alpha_stars[0]
        print(f"\n  Kandydat główny α*₃ = {main['alpha']:.8f}")
        print(f"  Stosunek α*₃/α*₁ = {main['alpha']/2.43183767:.6f}  (α*₁=2.43183767)")
        print(f"  Stosunek α*₃/α*₂ = {main['alpha']/2.63557742:.6f}  (α*₂=2.63557742)")
        print(f"  α*₃ - α_TGP = {main['alpha'] - 2.0:+.6f}")
        print(f"  α*₁+α*₂+α*₃ = {2.43183767+2.63557742+main['alpha']:.8f}  (suma trzech)")
        print(f"  φ·α*₁ = {PHI*2.43183767:.6f}  (czy α*₃ = φ·α*₁?)")
    else:
        print("  Brak znalezionych zer G₃=R31 — diagnostyka konieczna")

    print(f"\nCzas całkowity: {time.time()-t0:.1f} s")
    print("=" * 68)


if __name__ == '__main__':
    mp.freeze_support()
    main()
