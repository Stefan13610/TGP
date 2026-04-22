#!/usr/bin/env python3
"""
ex89_beta_pot_scan.py — Zależność zer F(α*)=R21 od β_pot (O-L2, v36)
=======================================================================
STATUS: LEGACY-TRANSLATIONAL

This script belongs to the older `alpha*` / `F(alpha)` / selection-language
chain. Keep it as legacy exploratory context rather than a canonical modern
`nbody` example.

Pytanie fizyczne:
  Czy dwa zera F(α*) = R21 są specyficzne dla β_pot=1,
  czy zachowują się dla dowolnego β_pot?

  Jeśli TOPOLOGICZNE → dwa zera dla wszystkich β_pot
  Jeśli NUMEROLOGIA → dwa zera tylko dla β_pot=1

Modyfikacja ODE: V'(g) = β_pot · g²(1-g)
  (β_pot=1: standardowy potencjał TGP, V_eff = g³/3 - g⁴/4)
  (β_pot≠1: deformacja coupling stałej γ = β_pot)

Parametry stałe:
  R_MAX    = 120
  B_WIN    = [28, 42]  (stałe, jak w ex88)
  FIT_WINS = [(20,36),(30,46),(40,56),(50,66)]
  N_SCAN   = 80  punktów w [2.0, 3.5]

BETA_LIST = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]

Kryteria sukcesu:
  P1: β_pot=1 → 2 zera (reprodukcja ex87v2)
  P2: β_pot≠1 → ≥1 zero przeżywa dla ≥4/6 wartości
  P3: α*₁(β)+α*₂(β) vs β_pot — zależność algebraiczna lub brak

Autor: Claudian (sesja v36, O-L2)
Data: 2026-03-28
"""
import sys, io
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# STAŁE
# ================================================================
PHI        = (1.0 + np.sqrt(5.0)) / 2.0
R21        = 206.7682830          # (m_mu/m_e)^4 = ratio docelowy
S_FORMULA  = 2.0 * np.pi - 1.1   # formuła do testowania (OCZEKUJEMY ≠)
G_OFF      = 0.005

# Stałe okna (jak w ex88 — nie skalowane!)
B_WIN_L  = 28.0
B_WIN_R  = 42.0
FIT_WINS = [(20, 36), (30, 46), (40, 56), (50, 66)]

R_MAX    = 120
N_SCAN   = 80    # punkty α w szerokimskanie

BETA_LIST = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]

PASS_CNT = 0; FAIL_CNT = 0
def chk(name, cond, detail=""):
    global PASS_CNT, FAIL_CNT
    tag = "[PASS]" if cond else "[FAIL]"
    if cond: PASS_CNT += 1
    else:    FAIL_CNT += 1
    print(f"  {tag}  {name}")
    if detail: print(f"        {detail}")

# ================================================================
# ODE z parametrem β_pot
# ================================================================
def integrate_soliton(g0, alpha, r_max, beta_pot=1.0):
    """Całkowanie solitu z V'(g) = beta_pot * g²(1-g)."""
    gb = np.exp(-1.0 / (2.0 * alpha)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-8)
        fg = 1.0 + 2.0 * alpha * np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        Vprime = beta_pot * g * g * (1.0 - g)   # DEFORMACJA β_pot
        curl   = (alpha / g) * gp * gp
        if r < 1e-10: return [gp, (Vprime - curl) / (3.0 * fg)]
        return [gp, (Vprime - curl - fg * 2.0 * gp / r) / fg]
    def ev_bounce(r, y): return y[0] - gb
    ev_bounce.terminal = True; ev_bounce.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(14):
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
    """Ekstrapolacja A_∞ ze stałych okien FIT_WINS."""
    Av, rv = [], []
    for rL, rR in FIT_WINS:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2))); rv.append(float(rL))
    if not Av: return 0.0
    if len(Av) < 2: return Av[-1]
    rv = np.array(rv); Av = np.array(Av)
    try:
        p, _ = curve_fit(lambda x, ai, a, b: ai * (1 + a/x + b/x**2),
                         rv, Av, p0=[Av[-1], 0., 0.], maxfev=2000)
        return float(p[0])
    except:
        try:
            p, _ = curve_fit(lambda x, ai, a: ai * (1 + a/x),
                             rv, Av, p0=[Av[-1], 0.], maxfev=1000)
            return float(p[0])
        except:
            return float(Av[-1])

def B_coeff(g0, alpha, r_max, beta_pot=1.0):
    """Amplituda cosinus w stałym oknie [B_WIN_L, B_WIN_R]."""
    r, g = integrate_soliton(g0, alpha, r_max, beta_pot)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max, beta_pot=1.0, lo=1.05, hi=2.5, n=40):
    """Znajdź pierwsze przejście B_coeff(g0)=0 w [lo, hi]."""
    gs = np.linspace(lo, hi, n); Bs = []
    for g in gs:
        try: Bs.append(B_coeff(g, alpha, r_max, beta_pot))
        except: Bs.append(np.nan)
    Bs = np.array(Bs)
    for i in range(len(Bs) - 1):
        if np.isnan(Bs[i]) or np.isnan(Bs[i+1]): continue
        if Bs[i] * Bs[i+1] < 0:
            try:
                return brentq(lambda x: B_coeff(x, alpha, r_max, beta_pot),
                              gs[i], gs[i+1], xtol=1e-9, maxiter=100)
            except: pass
    return None

def fz0(alpha, r_max, beta_pot=1.0):
    """Oblicz F(α) = (A_m/A_e)^4 dla danego β_pot."""
    z0 = find_z0(alpha, r_max, beta_pot)
    if z0 is None: return np.nan
    r, g   = integrate_soliton(z0, alpha, r_max, beta_pot)
    Ae = fit_amplitude(r, g)
    r2, g2 = integrate_soliton(PHI * z0, alpha, r_max, beta_pot)
    Am = fit_amplitude(r2, g2)
    if Ae < 1e-6: return np.nan
    return (Am / Ae) ** 4

# ================================================================
# SKAN α dla danego β_pot
# ================================================================
def scan_for_beta(beta_pot):
    """Skan α∈[2.0, 3.5] dla danego β_pot → znajdź zera F(α)=R21."""
    print(f"\n  --- β_pot = {beta_pot:.2f} ---", flush=True)
    alphas = np.linspace(2.0, 3.5, N_SCAN)
    fvals  = []
    for a in alphas:
        try:
            fv = fz0(a, R_MAX, beta_pot)
        except:
            fv = np.nan
        fvals.append(fv)
    fvals = np.array(fvals)

    # Znajdź przejścia przez R21
    diff   = fvals - R21
    crossings = []
    for i in range(len(diff) - 1):
        if np.isnan(diff[i]) or np.isnan(diff[i+1]): continue
        if diff[i] * diff[i+1] < 0:
            try:
                alpha_star = brentq(
                    lambda a: fz0(a, R_MAX, beta_pot) - R21,
                    alphas[i], alphas[i+1],
                    xtol=1e-4, maxiter=40
                )
                crossings.append(alpha_star)
            except:
                crossings.append((alphas[i] + alphas[i+1]) / 2.0)

    crossings = sorted(crossings)
    n_zeros = len(crossings)

    # Predykcja formuły sumy dla THIS β_pot: S = 2π - α_TGP/2 - β_pot/10 = 2π - 1 - β_pot/10
    S_pred = 2.0 * np.pi - 1.0 - beta_pot / 10.0

    if n_zeros == 0:
        print(f"    WYNIK: 0 zer F(α)=R21 → BRAK", flush=True)
    elif n_zeros == 1:
        a1 = crossings[0]
        print(f"    WYNIK: 1 zero: α* = {a1:.4f}", flush=True)
    elif n_zeros == 2:
        a1, a2 = crossings[0], crossings[1]
        S = a1 + a2
        S_ppm_pred = abs(S - S_pred) / S_pred * 1e6
        S_ppm_fix  = abs(S - S_FORMULA) / S_FORMULA * 1e6  # vs stała β=1
        print(f"    WYNIK: 2 zera: α*₁={a1:.4f}, α*₂={a2:.4f}", flush=True)
        print(f"    Suma S = {S:.4f}  |  S_pred(β={beta_pot:.2f}) = {S_pred:.4f}  Δ={S_ppm_pred:.0f} ppm", flush=True)
        print(f"    (vs S_formula_beta1 = {S_FORMULA:.4f}  Δ={S_ppm_fix:.0f} ppm)", flush=True)
    else:
        print(f"    WYNIK: {n_zeros} zer: {[f'{c:.4f}' for c in crossings]}", flush=True)

    return {'beta_pot': beta_pot, 'n_zeros': n_zeros, 'crossings': crossings,
            'S_pred': S_pred, 'fvals': fvals.tolist(), 'alphas': alphas.tolist()}

# ================================================================
# GŁÓWNA ANALIZA
# ================================================================
def main():
    print("=" * 68)
    print("ex89_beta_pot_scan.py — O-L2: Zależność zer F(α*) od β_pot")
    print("=" * 68)
    print(f"R_MAX    = {R_MAX}")
    print(f"B_WIN    = [{B_WIN_L}, {B_WIN_R}]")
    print(f"N_SCAN   = {N_SCAN} punktów w [2.0, 3.5]")
    print(f"R21      = {R21:.6f}")
    print(f"BETA_LIST= {BETA_LIST}")
    print()

    # Uruchom skan dla każdego β_pot
    # (równolegle dla szybkości)
    n_workers = min(len(BETA_LIST), cpu_count(), 6)
    print(f"Używam Pool({n_workers}) workerów...", flush=True)

    with Pool(n_workers) as pool:
        results = pool.map(scan_for_beta, BETA_LIST)

    # ============================
    # PODSUMOWANIE
    # ============================
    print("\n" + "=" * 68)
    print("PODSUMOWANIE — zer F(α*)=R21 vs β_pot")
    print("=" * 68)
    print(f"{'β_pot':>6} | {'n_zer':>5} | {'α*₁':>8} | {'α*₂':>8} | {'S_numer':>9} | {'S_pred':>9} | {'ΔS(ppm)':>9}")
    print("-" * 72)

    betas_with_2zeros = []
    betas_with_1zero  = []
    betas_with_0zeros = []
    all_sums = []    # (beta_pot, S_numeryczny, S_przewidywany)

    for res in results:
        bp   = res['beta_pot']
        nz   = res['n_zeros']
        cs   = res['crossings']
        Sp   = res.get('S_pred', 2.0 * np.pi - 1.0 - bp / 10.0)
        a1   = cs[0] if nz >= 1 else float('nan')
        a2   = cs[1] if nz >= 2 else float('nan')
        S    = a1 + a2 if nz >= 2 else float('nan')
        if not np.isnan(S):
            ppm = abs(S - Sp) / Sp * 1e6
            ppm_str = f"{ppm:.0f}"
            all_sums.append((bp, S, Sp))
        else:
            ppm_str = "---"
        a1s = f"{a1:.4f}" if not np.isnan(a1) else "---"
        a2s = f"{a2:.4f}" if not np.isnan(a2) else "---"
        Ss  = f"{S:.4f}" if not np.isnan(S) else "---"
        Sps = f"{Sp:.4f}"
        print(f"{bp:>6.2f} | {nz:>5d} | {a1s:>8} | {a2s:>8} | {Ss:>9} | {Sps:>9} | {ppm_str:>9}")

        if nz == 2: betas_with_2zeros.append(bp)
        elif nz == 1: betas_with_1zero.append(bp)
        else: betas_with_0zeros.append(bp)

    print()

    # ============================
    # KRYTERIA SUKCESU
    # ============================
    print("KRYTERIA SUKCESU:")

    # P1: β=1 → 2 zera
    res_beta1 = next((r for r in results if abs(r['beta_pot'] - 1.0) < 0.01), None)
    p1_pass = res_beta1 is not None and res_beta1['n_zeros'] == 2
    chk("P1: β_pot=1 → 2 zera (reprodukcja ex87v2)",
        p1_pass,
        f"β=1 daje {res_beta1['n_zeros']} zer" if res_beta1 else "brak β=1")

    # P2: β≠1 → ≥1 zero dla ≥4/5 wartości
    non_beta1 = [r for r in results if abs(r['beta_pot'] - 1.0) > 0.01]
    surviving = sum(1 for r in non_beta1 if r['n_zeros'] >= 1)
    p2_pass = surviving >= 4
    chk(f"P2: β_pot≠1 → ≥1 zero dla ≥4/5 wartości",
        p2_pass,
        f"Przeżywa: {surviving}/{len(non_beta1)}")

    # P3: czy S_numer(β_pot) ≈ S_pred(β_pot) = 2π - 1 - β_pot/10?
    if len(all_sums) >= 3:
        betas_s = np.array([x[0] for x in all_sums])
        sums_s  = np.array([x[1] for x in all_sums])
        preds_s = np.array([x[2] for x in all_sums])
        deviations = sums_s - preds_s
        ppm_devs = np.abs(deviations) / preds_s * 1e6

        print(f"\n  Test formuły sumy: S_numer vs S_pred = 2π - 1 - β_pot/10")
        print(f"  {'β_pot':>6} | {'S_numer':>9} | {'S_pred':>9} | {'Δ':>9} | {'ppm':>7}")
        print(f"  {'-'*50}")
        for bp_v, sv, pv in all_sums:
            dev = sv - pv
            ppm = abs(dev) / pv * 1e6
            print(f"  {bp_v:>6.2f} | {sv:>9.4f} | {pv:>9.4f} | {dev:>+9.4f} | {ppm:>7.0f}")

        max_ppm = np.max(ppm_devs)
        mean_ppm = np.mean(ppm_devs)
        formula_works = max_ppm < 5000  # tolerancja 5000 ppm = 0.5%

        # Alternatywnie: czy S jest STAŁA (niezależna od β_pot)?
        S_std = np.std(sums_s); S_mean = np.mean(sums_s)
        is_const = S_std / S_mean < 0.005 if S_mean > 0 else False
        S_range = np.max(sums_s) - np.min(sums_s)

        print(f"\n  Statystyki residuów:")
        print(f"  Max Δppm od S_pred: {max_ppm:.0f} ppm")
        print(f"  Avg Δppm od S_pred: {mean_ppm:.0f} ppm")
        print(f"  Zakres S_numer: {np.min(sums_s):.4f} – {np.max(sums_s):.4f}  (Δ = {S_range:.4f})")
        if formula_works:
            print(f"  → FORMUŁA SUMY DZIAŁA dla różnych β_pot (max err {max_ppm:.0f} ppm)")
        elif is_const:
            print(f"  → S STAŁA (σ/μ={S_std/S_mean:.5f}) — suma nie zależy od β_pot, formuła faktor 1/10 fałszywy!")
        else:
            print(f"  → S ZMIENNA z β_pot (zakres {S_range:.3f}), ale formuła 2π-1-β/10 NIE pasuje (max {max_ppm:.0f} ppm)")

        p3_pass = True  # P3 jest zawsze informacyjny
    else:
        print(f"  Niewystarczające dane do analizy formuły (tylko {len(all_sums)} β z 2 zerami)")
        p3_pass = False
    chk("P3: Analiza S(β_pot) vs formuła 2π-1-β_pot/10 wykonana", p3_pass,
        f"Liczba β z 2 zerami: {len(all_sums)}")

    # ============================
    # WNIOSKI FIZYCZNE
    # ============================
    print()
    print("WNIOSKI FIZYCZNE:")
    if len(betas_with_2zeros) >= 4:
        print(f"  ✅ DWA ZERA TOPOLOGICZNE: {len(betas_with_2zeros)}/6 wartości β_pot dają 2 zera")
        print(f"     → Struktura NIE jest specyficzna dla β_pot=1")
        print(f"     → SILNE wsparcie dla interpretacji mas leptonowych")
    elif len(betas_with_2zeros) >= 2:
        print(f"  ⚠️  DWA ZERA CZĘŚCIOWO TOPOLOGICZNE: {len(betas_with_2zeros)}/6 wartości")
        print(f"     → Możliwe, że zera znikają dla skrajnych β_pot")
    else:
        print(f"  ❌ DWA ZERA TYLKO DLA β_pot≈1: {len(betas_with_2zeros)}/6 wartości")
        print(f"     → Dwa zera mogą być numerologią β=1")

    total = PASS_CNT + FAIL_CNT
    print(f"\nWYNIK KOŃCOWY: {PASS_CNT}/{total} PASS")
    print(f"Wartości β_pot z 2 zerami: {betas_with_2zeros}")
    print(f"Wartości β_pot z 1 zerem:  {betas_with_1zero}")
    print(f"Wartości β_pot z 0 zerami: {betas_with_0zeros}")

if __name__ == '__main__':
    main()
