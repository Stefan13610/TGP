#!/usr/bin/env python3
"""
ex98_tau_multiplier.py — Jaki mnożnik m* daje (A(m·z₀)/A(z₀))^4 = R₃₁?
==========================================================================
Wyniki ex96/ex97:
  G₃(φ²·z₀)_far > R₃₁ wszędzie na α∈[2.5,5.0]
  Plateau A_tau potwierdzone do r=286

Pytanie: jaki mnożnik m* ≠ φ² daje G₃=(A(m·z₀)/A(z₀))^4 = R₃₁?
  - Czy m*(α) jest stały (= uniwersalny mnożnik)?
  - Czy m* jest algebraicznie specjalny?
  - Czy warunek G₃/F=(A_tau/A_mu)^4=R₃₂ z m=φ² ma rozwiązanie?

Plan:
  FAZA 1: Precyzyjne α*₁_far, α*₂_far — brentq F_far=R₂₁
  FAZA 2: Profil A(g₀) vs g₀ przy α∈{α*₁_far, α*₂_far, 2.750, 3.000}
  FAZA 3: brentq m*(α) dla (A(m·z₀)/A(z₀))^4 = R₃₁
  FAZA 4: Test algebraicznych kandydatów na m*
  FAZA 5: Sprawdzenie warunku G₃/F = R₃₂ z FAR oknami

Autor: Claudian (sesja v38, O-L5 nowe)
Data: 2026-03-29
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
# STAŁE
# ================================================================
PHI  = (1.0 + np.sqrt(5.0)) / 2.0
PHI2 = PHI * PHI
PHI3 = PHI2 * PHI

M_TAU_MEV = 1776.86
M_MU_MEV  = 105.658
M_E_MEV   = 0.51100

R21 = M_MU_MEV / M_E_MEV    # ≈ 206.767
R31 = M_TAU_MEV / M_E_MEV   # ≈ 3477.221
R32 = M_TAU_MEV / M_MU_MEV  # ≈ 16.817

G_OFF   = 0.005
B_WIN_L = 28.0
B_WIN_R = 42.0
R_MAX   = 150

FIT_WINS_FAR = [(60, 76), (70, 86), (80, 96), (90, 106),
                (100, 116), (110, 126), (120, 136)]

# ================================================================
# ODE
# ================================================================
def integrate_soliton(g0, alpha, r_max=R_MAX):
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
    for _ in range(25):
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

def amplitude_far(r, g):
    """Amplituda asymptotyczna z okien FAR — z testem plateau."""
    Av, rv = [], []
    for rL, rR in FIT_WINS_FAR:
        if rR > r[-1]: break
        mask = (r >= rL) & (r <= rR)
        if mask.sum() < 15: continue
        rf = r[mask]; df = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
        Av.append(float(np.sqrt(bc[0]**2 + bc[1]**2)))
        rv.append(float(rL))
    if not Av: return 0.0
    Av = np.array(Av)
    cv = np.std(Av) / (np.mean(Av) + 1e-12)
    if cv < 0.005 or len(Av) < 3:
        return float(np.median(Av))
    try:
        p, _ = curve_fit(lambda x, ai, a, b: ai * (1+a/x+b/x**2),
                         np.array(rv), Av, p0=[Av[-1],0.,0.], maxfev=2000)
        return float(p[0])
    except:
        return float(np.median(Av))

def B_coeff_func(g0, alpha, r_max=R_MAX):
    r, g = integrate_soliton(g0, alpha, r_max)
    mask = (r >= B_WIN_L) & (r <= B_WIN_R)
    if mask.sum() < 15: return 0.0
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc, _, _, _ = np.linalg.lstsq(M, df, rcond=None)
    return float(bc[0])

def find_z0(alpha, r_max=R_MAX, lo=1.05, hi=2.0, n=60):
    gs = np.linspace(lo, hi, n)
    Bs = []
    for gg in gs:
        try: Bs.append(B_coeff_func(gg, alpha, r_max))
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

def compute_F_G3(alpha, r_max=R_MAX):
    """Oblicz F i G₃ z FAR. Zwraca dict lub None."""
    z0 = find_z0(alpha, r_max)
    if z0 is None: return None
    r0, g0_s = integrate_soliton(z0,      alpha, r_max)
    r1, g1_s = integrate_soliton(PHI*z0,  alpha, r_max)
    r2, g2_s = integrate_soliton(PHI2*z0, alpha, r_max)
    Ae  = amplitude_far(r0, g0_s)
    Amu = amplitude_far(r1, g1_s)
    Atau = amplitude_far(r2, g2_s)
    if Ae < 1e-8: return None
    return {'alpha': alpha, 'z0': z0,
            'Ae': Ae, 'Amu': Amu, 'Atau': Atau,
            'F': (Amu/Ae)**4, 'G3': (Atau/Ae)**4}

def worker_FG(args):
    alpha, r_max = args
    try: return compute_F_G3(alpha, r_max)
    except: return None

def amplitude_at_g0(g0, alpha, r_max=R_MAX):
    """Amplituda asymptotyczna solitonu startującego z g₀."""
    r, g = integrate_soliton(g0, alpha, r_max)
    return amplitude_far(r, g)

# ================================================================
# MAIN
# ================================================================
def main():
    t0_total = time.time()
    ncpu = 6

    print("=" * 70)
    print("ex98_tau_multiplier.py — Jaki mnożnik m* daje G₃=R₃₁? (O-L5 nowe)")
    print("=" * 70)
    print(f"  R21={R21:.6f}, R31={R31:.6f}, R32={R32:.6f}")
    print(f"  R21^(1/4)={R21**(0.25):.6f}, R31^(1/4)={R31**(0.25):.6f}")
    print(f"  PHI={PHI:.8f}, PHI²={PHI2:.8f}")
    print(f"  R_MAX={R_MAX}")
    print()

    # ----------------------------------------------------------------
    # FAZA 1: Precyzyjne α*₁_far, α*₂_far
    # ----------------------------------------------------------------
    print("=" * 70)
    print("FAZA 1: Precyzyjne α*₁_far, α*₂_far (F_far = R₂₁)")
    print("=" * 70)

    def F_minus_R21(alpha):
        res = compute_F_G3(alpha)
        if res is None: return float('nan')
        return res['F'] - R21

    # α*₁_far między 2.43 i 2.44 (z ex96 dense scan)
    alpha_star1_far, alpha_star2_far = None, None
    try:
        alpha_star1_far = brentq(F_minus_R21, 2.430, 2.445, xtol=1e-9, maxiter=200)
        res1 = compute_F_G3(alpha_star1_far)
        print(f"  α*₁_far = {alpha_star1_far:.9f}")
        if res1:
            print(f"    z₀={res1['z0']:.7f}, φ·z₀={PHI*res1['z0']:.7f}")
            print(f"    F_far={res1['F']:.6f} (R₂₁={R21:.6f}), G₃_far={res1['G3']:.4f}")
    except Exception as e:
        print(f"  α*₁_far: brentq failed: {e}")

    try:
        alpha_star2_far = brentq(F_minus_R21, 2.750, 2.760, xtol=1e-9, maxiter=200)
        res2 = compute_F_G3(alpha_star2_far)
        print(f"  α*₂_far = {alpha_star2_far:.9f}")
        if res2:
            print(f"    z₀={res2['z0']:.7f}, φ·z₀={PHI*res2['z0']:.7f}")
            print(f"    F_far={res2['F']:.6f} (R₂₁={R21:.6f}), G₃_far={res2['G3']:.4f}")
    except Exception as e:
        print(f"  α*₂_far: brentq failed: {e}")

    if alpha_star1_far and alpha_star2_far:
        S2_far = alpha_star1_far + alpha_star2_far
        print(f"\n  S₂_far = α*₁+α*₂ = {S2_far:.9f}")
        print(f"  S₂_std = {2.43183767 + 2.63557742:.9f}  (ex93)")
        print(f"  Δ(S₂)  = {S2_far - (2.43183767+2.63557742):+.6f}")

    # ----------------------------------------------------------------
    # FAZA 2: Profil A(g₀) vs g₀ przy kilku α
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 2: Profil A(g₀) vs g₀ przy wybranych α")
    print("=" * 70)

    test_alphas_2 = []
    if alpha_star1_far: test_alphas_2.append((alpha_star1_far, "α*₁_far"))
    if alpha_star2_far: test_alphas_2.append((alpha_star2_far, "α*₂_far"))
    test_alphas_2 += [(2.750, "2.750"), (3.000, "3.000")]

    # g₀ grid: z₀ do 5·z₀, krok 0.1·z₀
    profile_data = {}
    for alpha, label in test_alphas_2:
        z0 = find_z0(alpha)
        if z0 is None: continue
        print(f"\n  α={alpha:.6f} ({label}), z₀={z0:.6f}")
        print(f"  {'g₀':>8} {'g₀/z₀':>8} | {'A(g₀)':>10} | {'(A/Ae)^4':>12} | '(A/Ae)^4-R31':>12")
        print(f"  {'-'*68}")

        # Ae najpierw
        r0, g0_s = integrate_soliton(z0, alpha)
        Ae = amplitude_far(r0, g0_s)

        g0_list = [z0 * m for m in np.arange(1.0, 5.05, 0.1)]
        A_list  = []
        for g0_val in g0_list:
            A = amplitude_at_g0(g0_val, alpha)
            A_list.append(A)
            G3_val = (A/Ae)**4 if Ae > 1e-8 else float('nan')
            print(f"  {g0_val:8.5f} {g0_val/z0:8.4f} | {A:10.6f} | {G3_val:12.4f} | {G3_val-R31:+12.4f}")

        profile_data[alpha] = {'z0': z0, 'Ae': Ae, 'g0': g0_list, 'A': A_list}

    # ----------------------------------------------------------------
    # FAZA 3: brentq m*(α) — szukamy (A(m·z₀)/A(z₀))^4 = R₃₁
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 3: brentq m*(α) dla G₃=(A(m·z₀)/A(z₀))^4 = R₃₁")
    print("=" * 70)

    scan_alphas_3 = [2.43, 2.50, 2.60, 2.70, 2.75, 2.80, 2.90, 3.00, 3.10, 3.20, 3.50]
    if alpha_star1_far: scan_alphas_3.insert(0, alpha_star1_far)
    if alpha_star2_far:
        if alpha_star2_far not in scan_alphas_3:
            scan_alphas_3.append(alpha_star2_far)

    print(f"\n  {'α':>9} | {'z₀':>8} | {'m*':>10} | {'m*/φ':>8} | {'m*/φ²':>8} | {'G₃(m*)':>10} | Uwagi")
    print(f"  {'-'*78}")

    mstar_results = {}
    for alpha in sorted(set(scan_alphas_3)):
        z0 = find_z0(alpha)
        if z0 is None:
            print(f"  {alpha:9.5f} | z₀=None")
            continue
        Ae = amplitude_far(*integrate_soliton(z0, alpha))
        if Ae < 1e-8: continue

        target_A = R31**(0.25) * Ae   # A(m·z₀) potrzebne

        def G3_minus_R31(m):
            g0_val = m * z0
            A = amplitude_at_g0(g0_val, alpha)
            return (A/Ae)**4 - R31

        # Sprawdź granice: m∈[1.2, 5.0]
        try:
            v_lo = G3_minus_R31(1.2)
            v_hi = G3_minus_R31(5.0)
        except:
            print(f"  {alpha:9.5f} | {z0:.6f} | error computing bounds")
            continue

        if v_lo * v_hi > 0:
            # Szukaj dokładniej
            ms_test = np.linspace(1.2, 5.0, 40)
            vs_test = []
            for m in ms_test:
                try: vs_test.append(G3_minus_R31(m))
                except: vs_test.append(float('nan'))
            bracket = None
            for i in range(len(vs_test)-1):
                if np.isnan(vs_test[i]) or np.isnan(vs_test[i+1]): continue
                if vs_test[i] * vs_test[i+1] < 0:
                    bracket = (ms_test[i], ms_test[i+1])
                    break
            if bracket is None:
                # Może G₃ jest zawsze > R31 lub zawsze < R31
                min_G3 = min([(A/Ae)**4 for m in ms_test
                              for A in [amplitude_at_g0(m*z0, alpha)]
                              if A > 0], default=float('nan'))
                note = f"Brak — min(G₃)_scan={min_G3:.1f}" if not np.isnan(min_G3) else "Brak bracketa"
                print(f"  {alpha:9.5f} | {z0:.6f} | {'N/A':>10} | {'—':>8} | {'—':>8} | {'—':>10} | {note}")
                continue
        else:
            bracket = (1.2, 5.0)

        try:
            m_star = brentq(G3_minus_R31, bracket[0], bracket[1], xtol=1e-8, maxiter=200)
            G3_check = compute_F_G3_at_m(alpha, z0, Ae, m_star)
            note = f"m*/φ²={m_star/PHI2:.5f}, m*/φ={m_star/PHI:.5f}"
            print(f"  {alpha:9.5f} | {z0:.6f} | {m_star:10.7f} | {m_star/PHI:8.6f} | {m_star/PHI2:8.6f} | {G3_check:10.4f} | ←")
            mstar_results[alpha] = m_star
        except Exception as e:
            print(f"  {alpha:9.5f} | {z0:.6f} | brentq failed: {e}")

    # ----------------------------------------------------------------
    # FAZA 4: Kandydaci algebraiczni na m*
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 4: Analiza algebraiczna m*(α)")
    print("=" * 70)

    if mstar_results:
        mvals = np.array(list(mstar_results.values()))
        avals = np.array(list(mstar_results.keys()))
        print(f"\n  m* ∈ [{mvals.min():.6f}, {mvals.max():.6f}]")
        print(f"  m* mean = {mvals.mean():.6f}, std = {mvals.std():.6f}")
        print(f"  m* stały? cv = {mvals.std()/mvals.mean()*100:.2f}%")
        print()

        # Kandydaci
        candidates = [
            ("φ²",          PHI2),
            ("φ·√φ = φ^(3/2)", PHI**1.5),
            ("φ+1 = φ²",    PHI+1),
            ("e",           np.e),
            ("√6",          np.sqrt(6)),
            ("√5",          np.sqrt(5)),
            ("√7",          np.sqrt(7)),
            ("3φ-2",        3*PHI-2),
            ("2φ-1 = √5",   2*PHI-1),
            ("2φ",          2*PHI),
            ("1+φ+1/φ",     1+PHI+1/PHI),
            ("R31^(1/4)/R21^(1/4)", R31**(0.25)/R21**(0.25)),
            ("R32^(1/4)",   R32**(0.25)),
            ("R32^(1/4)+φ", R32**(0.25)+PHI),
            ("R31^(1/3)/π", R31**(1/3)/np.pi),
        ]

        print(f"  Porównanie m* (mean) z kandydatami algebraicznymi:")
        print(f"  {'Kandydat':30s} | {'Wartosc':10s} | {'D(mean)':>12s} | {'ppm':>10s}")
        print(f"  {'-'*68}")
        for name, val in sorted(candidates, key=lambda x: abs(x[1]-mvals.mean())):
            delta = mvals.mean() - val
            ppm = delta / mvals.mean() * 1e6
            print(f"  {name:30s} | {val:10.6f} | {delta:+12.6f} | {ppm:+10.1f}")

        # Czy m*(α) zależy od α?
        if len(avals) > 2:
            p = np.polyfit(avals, mvals, 1)
            print(f"\n  Liniowy fit m*(α) = {p[0]:.6f}·α + {p[1]:.6f}")
            print(f"  dm*/dα = {p[0]:+.6f} (0 = stały mnożnik)")

    # ----------------------------------------------------------------
    # FAZA 5: G₃/F = R₃₂ z FAR — szukamy α*_τμ
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("FAZA 5: Warunek G₃_far/F_far = R₃₂ — szukamy α*_τμ_far")
    print("=" * 70)

    def GF_minus_R32(alpha):
        res = compute_F_G3(alpha)
        if res is None or res['F'] < 1e-4: return float('nan')
        return res['G3'] / res['F'] - R32

    # Z ex96: G₃/F ≈ 23-25 wszędzie — ale sprawdźmy dokładnie szeroki zakres
    print(f"  Skan G₃_far/F_far na α∈[2.4, 5.0], krok=0.05...")
    args_gf = [(a, R_MAX) for a in np.arange(2.40, 5.01, 0.05)]
    with mp.Pool(ncpu) as pool:
        gf_results = pool.map(worker_FG, args_gf)

    sign_changes = []
    prev_sign = None
    prev_alpha = None
    print(f"\n  {'α':>7} | {'F_far':>9} | {'G3_far':>10} | {'G3/F':>9} | {'G3/F-R32':>10}")
    print(f"  {'-'*60}")
    for r in [x for x in gf_results if x]:
        a = r['alpha']
        gf = r['G3'] / r['F'] if r['F'] > 1e-4 else float('nan')
        diff = gf - R32
        print(f"  {a:7.3f} | {r['F']:9.4f} | {r['G3']:10.4f} | {gf:9.4f} | {diff:+10.4f}")
        cur_sign = np.sign(diff) if not np.isnan(diff) else None
        if prev_sign and cur_sign and cur_sign != prev_sign:
            sign_changes.append((prev_alpha, a))
        prev_sign = cur_sign; prev_alpha = a

    if sign_changes:
        print(f"\n  Zmiany znaku G₃/F−R₃₂: {sign_changes}")
        for a_lo, a_hi in sign_changes:
            try:
                a_star = brentq(GF_minus_R32, a_lo, a_hi, xtol=1e-8, maxiter=200)
                res = compute_F_G3(a_star)
                if res:
                    gf = res['G3']/res['F']
                    print(f"  ✓ α*_τμ_far = {a_star:.9f}")
                    print(f"    G₃/F={gf:.6f} (R₃₂={R32:.6f}), F={res['F']:.4f}, G₃={res['G3']:.4f}")
                    print(f"    z₀={res['z0']:.7f}")
                    S3 = 2.43183767 + 2.63557742 + a_star  # STD alpha*1,2
                    print(f"    S₃(ze STD α*₁,₂) = {S3:.9f}")
                    if alpha_star1_far and alpha_star2_far:
                        S3_far = alpha_star1_far + alpha_star2_far + a_star
                        print(f"    S₃(FAR α*₁,₂) = {S3_far:.9f}")
            except Exception as e:
                print(f"  brentq failed: {e}")
    else:
        print(f"\n  Brak zmiany znaku G₃/F−R₃₂ na [2.4,5.0].")
        gf_vals = [(r['alpha'], r['G3']/r['F']) for r in gf_results if r and r['F']>1e-4]
        if gf_vals:
            min_gf = min(gf_vals, key=lambda x: abs(x[1]-R32))
            print(f"  Min |G₃/F−R₃₂|: α={min_gf[0]:.3f}, G₃/F={min_gf[1]:.4f} (R₃₂={R32:.4f})")

    # ----------------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------------
    print()
    print("=" * 70)
    print("PODSUMOWANIE ex98 — Mnożnik m* dla sektora τ")
    print("=" * 70)
    print(f"Czas: {time.time()-t0_total:.1f} s")
    print()
    if alpha_star1_far and alpha_star2_far:
        print(f"  α*₁_far = {alpha_star1_far:.7f}")
        print(f"  α*₂_far = {alpha_star2_far:.7f}")
        print(f"  S₂_far  = {alpha_star1_far+alpha_star2_far:.7f}")
    if mstar_results:
        mvals = np.array(list(mstar_results.values()))
        print(f"  m*(mean)= {mvals.mean():.7f} ± {mvals.std():.7f}")
        print(f"  Wartość m* = {'stała' if mvals.std()/mvals.mean()<0.005 else 'zależy od α'}")
    print("=" * 70)


def compute_F_G3_at_m(alpha, z0, Ae, m):
    """Oblicz G₃ przy danym m (pomocnicza)."""
    A = amplitude_at_g0(m * z0, alpha)
    return (A/Ae)**4 if Ae > 1e-8 else float('nan')


if __name__ == '__main__':
    mp.freeze_support()
    main()
