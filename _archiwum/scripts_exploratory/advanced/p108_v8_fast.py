#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v8_fast.py  --  TGP v1
=====================================
Ścieżka 5: argmax K*₂/K*₁ — SZYBKA wersja

ZMIANA vs v7:
  v7 błąd: skan K na siatce n_K=30 punktów → spacing 0.001 > width 0.0005 zmian znaku → MISS
  v8 fix:  BRENTQ bezpośrednio w K (brak skanowania siatki K)
           wystarczy sprawdzić znaki g_U na dwóch końcach [K_lo, K_hi],
           potem brentq konwerguje niezależnie od szerokości zmiany znaku.

  Przyspieszenie: RTOL 1e-10→1e-7 (3-5× szybciej), brentq zamiast skanowania

Dane: 2026-03-26
"""
import sys, io, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Stałe ─────────────────────────────────────────────────────────────────────
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
AGAM   = 0.040
R_MAX  = 80.0
N_EVAL = 500
RTOL   = 1e-7    # Relaxed z 1e-10 → 3-5× szybciej
ATOL   = 1e-9

AK_REF  = 8.5616
K1_REF  = 0.010280;  PSI1_REF = 1.2419
K2_REF  = 0.09999348; PSI2_REF = 2.7624

DPSI_DK_B = 18.78   # dψ₀^B/dK w punkcie referencyjnym (z p107_v3)
PSI_B_MIN  = 1.80    # Branch B zawsze ψ₀ > 1.80 (Branch A < 1.70)
PSI_A_MAX  = 1.70    # Branch A zawsze ψ₀ < 1.70

# ── ODE ───────────────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def solve_one(K, psi0, ak):
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
        kfac = 1.0 + ak/phi
        return [dphi, dV_mod(phi)/kfac + ak*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL, t_eval=r_ev,
                        dense_output=False)
        if sol.t[-1] < R_MAX*0.99:
            return np.nan, np.nan
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        kfac = 1. + ak/phi
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Vv   = np.array([V_mod(float(q)) for q in phi])
        Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
        return sol.y[0, -1], Ek+Ep
    except Exception:
        return np.nan, np.nan


def find_psi0_brentq(K, ak, psi_c, dpsi, n_pts, psi_floor=0.80):
    """
    Znajdź ψ₀ w oknie [max(psi_floor, psi_c-dpsi), psi_c+dpsi]
    tak że φ(R_MAX) = 1. Wybiera rozwiązanie NAJBLIŻSZE psi_c.
    """
    lo = max(psi_floor, psi_c - dpsi)
    hi = psi_c + dpsi
    if lo >= hi: return np.nan, np.nan

    pts = np.linspace(lo, hi, n_pts)
    fv  = []
    for p in pts:
        pe, _ = solve_one(K, p, ak)
        fv.append(pe - 1.0 if np.isfinite(pe) else np.nan)

    candidates = []
    for i in range(len(fv)-1):
        fi, fi1 = fv[i], fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)):
            continue
        if fi * fi1 < 0:
            try:
                fn = lambda p: (solve_one(K, p, ak)[0] - 1.0)
                p0 = brentq(fn, pts[i], pts[i+1], xtol=1e-9, maxiter=40)
                candidates.append(p0)
            except Exception:
                t = -fi/(fi1-fi)
                candidates.append(pts[i] + t*(pts[i+1]-pts[i]))

    if not candidates:
        return np.nan, np.nan
    best = min(candidates, key=lambda p: abs(p - psi_c))
    _, E = solve_one(K, best, ak)
    return best, E


def gU_at_K(K, ak, psi_c, dpsi, n_pts, psi_floor=0.80):
    """g_U = E/(4πK) − 1 z trackingiem ψ₀."""
    p0, E = find_psi0_brentq(K, ak, psi_c, dpsi, n_pts, psi_floor=psi_floor)
    if not np.isfinite(p0): return np.nan, psi_c
    if not (np.isfinite(E) and K > 0): return np.nan, p0
    return E/(4*np.pi*K) - 1.0, p0


def find_K_star_brentq(ak, K_seed, psi_seed,
                        dpsi, n_pts, psi_floor,
                        K_frac=0.14, label=""):
    """
    Znajdź K* metodą brentq w K.
    Zakres: [K_seed*(1-K_frac), K_seed*(1+K_frac)]
    ψ₀ śledzony liniowo z DPSI_DK_B (dla Branch B) lub stały (dla Branch A).
    """
    K_lo = K_seed * (1 - K_frac)
    K_hi = K_seed * (1 + K_frac)

    # Inicjalne centra ψ₀ na końcach
    is_B = (psi_floor >= 1.5)
    if is_B:
        psi_lo = psi_seed + DPSI_DK_B * (K_lo - K_seed)
        psi_hi = psi_seed + DPSI_DK_B * (K_hi - K_seed)
    else:
        psi_lo = psi_seed
        psi_hi = psi_seed

    g_lo, p_lo = gU_at_K(K_lo, ak, psi_lo, dpsi, n_pts, psi_floor)
    g_hi, p_hi = gU_at_K(K_hi, ak, psi_hi, dpsi, n_pts, psi_floor)

    # Fallback: rozszerz zakres jeśli brak zmiany znaku
    if not (np.isfinite(g_lo) and np.isfinite(g_hi) and g_lo*g_hi < 0):
        K_lo2 = K_seed * (1 - 2*K_frac)
        K_hi2 = K_seed * (1 + 2*K_frac)
        if is_B:
            psi_lo2 = psi_seed + DPSI_DK_B*(K_lo2 - K_seed)
            psi_hi2 = psi_seed + DPSI_DK_B*(K_hi2 - K_seed)
        else:
            psi_lo2 = psi_seed; psi_hi2 = psi_seed
        g_lo2, p_lo2 = gU_at_K(K_lo2, ak, psi_lo2, dpsi*1.2, n_pts, psi_floor)
        g_hi2, p_hi2 = gU_at_K(K_hi2, ak, psi_hi2, dpsi*1.2, n_pts, psi_floor)
        if np.isfinite(g_lo2) and np.isfinite(g_hi2) and g_lo2*g_hi2 < 0:
            K_lo, K_hi = K_lo2, K_hi2
            g_lo, p_lo = g_lo2, p_lo2
            g_hi, p_hi = g_hi2, p_hi2
        else:
            return np.nan, psi_seed, np.nan

    # Brentq w K
    psi_ref = [p_lo if np.isfinite(p_lo) else psi_lo]
    dpsi_inner = dpsi * 0.85

    def fn(K):
        if is_B:
            psi_c = psi_ref[0] + DPSI_DK_B*(K - K_lo)
        else:
            psi_c = psi_ref[0]
        g, p = gU_at_K(K, ak, psi_c, dpsi_inner, n_pts, psi_floor)
        if np.isfinite(p): psi_ref[0] = p
        return g if np.isfinite(g) else g_lo

    try:
        K_star = brentq(fn, K_lo, K_hi, xtol=1e-8, maxiter=60)
        if is_B:
            psi_c_f = psi_ref[0] + DPSI_DK_B*(K_star - K_lo)
        else:
            psi_c_f = psi_ref[0]
        g_f, p_f = gU_at_K(K_star, ak, psi_c_f, dpsi*0.7, n_pts, psi_floor)
        if not np.isfinite(g_f): g_f = fn(K_star)
        p_out = p_f if np.isfinite(p_f) else psi_c_f
        return K_star, p_out, g_f
    except Exception as e:
        return np.nan, psi_seed, np.nan


# ─── Main ────────────────────────────────────────────────────────────────────
if __name__ == '__main__':

    print("="*72)
    print("p108_v8: Branch tracking z BRENTQ w K (bez skanowania siatki K)")
    print(f"  RTOL={RTOL:.0e}, N_EVAL={N_EVAL}")
    print(f"  Δα_K=0.05, K_frac=±14%")
    print(f"  Referencja: α_K={AK_REF}, K*₂={K2_REF:.8f}, ψ₀^B={PSI2_REF:.4f}")
    print("="*72)
    sys.stdout.flush()

    # ── Timing test ───────────────────────────────────────────────────────────
    t0 = time.time()
    _ = solve_one(K2_REF, PSI2_REF, AK_REF)
    t_call = time.time() - t0
    print(f"\n[TIMING] 1 wywołanie ODE: {t_call:.3f}s")
    sys.stdout.flush()

    # ── Weryfikacja Branch B przy AK_REF ─────────────────────────────────────
    print(f"\n[0] Weryfikacja Branch B przy α_K={AK_REF}:")
    sys.stdout.flush()

    for K_test, label_expected in [(0.090, "ujemne"), (0.110, "dodatnie")]:
        psi_exp = PSI2_REF + DPSI_DK_B*(K_test - K2_REF)
        g_t, p_t = gU_at_K(K_test, AK_REF, psi_exp, dpsi=0.32, n_pts=14,
                            psi_floor=PSI_B_MIN)
        print(f"  K={K_test:.3f}: psi_c={psi_exp:.4f}, ψ₀={p_t:.4f}, "
              f"g_U={g_t:.4f}  (oczekiwane {label_expected})")
        sys.stdout.flush()

    # K*₂ z brentq
    K2_v, psi2_v, gU2_v = find_K_star_brentq(
        AK_REF, K2_REF, PSI2_REF,
        dpsi=0.32, n_pts=14, psi_floor=PSI_B_MIN, K_frac=0.14, label="B_ref")
    if np.isfinite(K2_v):
        d2 = abs(K2_v - K2_REF)/K2_REF*100
        print(f"\n  K*₂ = {K2_v:.8f}  (ref={K2_REF:.8f})  Δ={d2:.3f}%  "
              f"g_U={gU2_v:.2e}  ψ₀={psi2_v:.4f}  "
              f"{'✓' if d2 < 2.0 else '✗ BAD'}")
        K2_start, PSI2_start = K2_v, psi2_v
    else:
        print(f"  ✗ Branch B nie znaleziono — używam ref")
        K2_start, PSI2_start = K2_REF, PSI2_REF

    # K*₁ Branch A
    K1_v, psi1_v, gU1_v = find_K_star_brentq(
        AK_REF, K1_REF, PSI1_REF,
        dpsi=0.15, n_pts=10, psi_floor=0.80, K_frac=0.18, label="A_ref")
    if np.isfinite(K1_v):
        d1 = abs(K1_v - K1_REF)/K1_REF*100
        print(f"  K*₁ = {K1_v:.8f}  Δ={d1:.3f}%  g_U={gU1_v:.2e}  ψ₀={psi1_v:.4f}")
        K1_start, PSI1_start = K1_v, psi1_v
    else:
        print(f"  ✗ Branch A nie znaleziono — używam ref")
        K1_start, PSI1_start = K1_REF, PSI1_REF

    print()
    sys.stdout.flush()

    # ── Scan ──────────────────────────────────────────────────────────────────
    DAK = 0.05

    def build_ak_list(start, stop, step):
        vals = [start]
        v = start
        while True:
            v = round(v + step, 8)
            if step > 0 and v > stop + 1e-9: break
            if step < 0 and v < stop - 1e-9: break
            vals.append(v)
        return vals

    AK_desc = build_ak_list(AK_REF, 6.0,  -DAK)
    AK_asc  = build_ak_list(AK_REF, 9.00, +DAK)

    print(f"Descending: {len(AK_desc)} kroków → α_K=6.0")
    print(f"Ascending:  {len(AK_asc)} kroków  → α_K=9.0")
    sys.stdout.flush()

    def run_scan(AK_values, k1_s, psi1_s, k2_s, psi2_s, label):
        results = []
        k1_c, psi1_c = k1_s, psi1_s
        k2_c, psi2_c = k2_s, psi2_s
        b2_fails = 0

        print(f"\n{'='*68}")
        print(f"Skan {label}:")
        print(f"  {'α_K':>7}  {'K*₁':>11}  {'K*₂':>12}  {'ratio':>9}  "
              f"{'ψ₀^B':>7}  {'dt':>5}")
        print(f"  {'-'*68}")
        sys.stdout.flush()

        for ak in AK_values:
            t0_ = time.time()

            # Branch A
            K1n, psi1n, gU1n = find_K_star_brentq(
                ak, k1_c, psi1_c,
                dpsi=0.15, n_pts=10, psi_floor=0.80, K_frac=0.20)
            if np.isfinite(K1n) and np.isfinite(gU1n) and abs(gU1n) < 0.1:
                k1_c, psi1_c = K1n, psi1n
            else:
                K1n = np.nan

            # Branch B
            K2n = np.nan
            if b2_fails < 6:
                K2n, psi2n, gU2n = find_K_star_brentq(
                    ak, k2_c, psi2_c,
                    dpsi=0.35, n_pts=14, psi_floor=PSI_B_MIN, K_frac=0.16)
                if np.isfinite(K2n) and np.isfinite(gU2n) and abs(gU2n) < 0.1:
                    k2_c, psi2_c = K2n, psi2n
                    b2_fails = 0
                else:
                    K2n = np.nan
                    b2_fails += 1

            ratio = (K2n/K1n
                     if (np.isfinite(K1n) and np.isfinite(K2n) and K1n > 0)
                     else np.nan)
            dt = time.time() - t0_

            k1s = f"{K1n:.6f}" if np.isfinite(K1n) else "---"
            k2s = f"{K2n:.8f}" if np.isfinite(K2n) else "---"
            rs  = f"{ratio:.4f}" if np.isfinite(ratio) else "---"
            ref = " <REF" if abs(ak - AK_REF) < 1e-3 else ""
            print(f"  {ak:>7.4f}  {k1s:>11}  {k2s:>12}  {rs:>9}  "
                  f"{psi2_c:>7.4f}  {dt:>4.0f}s{ref}")
            sys.stdout.flush()

            results.append({'ak': ak, 'K1': K1n, 'K2': K2n,
                            'psi2': psi2_c, 'ratio': ratio})

        return results

    t_total = time.time()
    res_d = run_scan(AK_desc, K1_start, PSI1_start, K2_start, PSI2_start, "MALEJĄCY")
    res_a = run_scan(AK_asc,  K1_start, PSI1_start, K2_start, PSI2_start, "ROSNĄCY")
    print(f"\nŁączny czas: {time.time()-t_total:.0f}s")
    sys.stdout.flush()

    # ── Tabela wynikowa ────────────────────────────────────────────────────────
    all_res = {}
    for r in res_d + res_a:
        ak = round(r['ak'], 8)
        if ak not in all_res:
            all_res[ak] = r
        elif np.isfinite(r['ratio']) and not np.isfinite(all_res[ak]['ratio']):
            all_res[ak] = r

    all_sorted = sorted(all_res.values(), key=lambda x: x['ak'])
    valid = [(r['ak'], r['K1'], r['K2'], r['ratio'])
             for r in all_sorted if np.isfinite(r['ratio'])]

    print("\n" + "="*72)
    print("TABELA WYNIKOWA:")
    print(f"  {'α_K':>7}  {'K*₁':>11}  {'K*₂':>12}  {'ratio K*₂/K*₁':>14}")
    for ak, K1, K2, rat in valid:
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        print(f"  {ak:>7.4f}  {K1:>11.6f}  {K2:>12.8f}  {rat:>14.4f}{ref}")

    if len(valid) >= 3:
        aks  = np.array([r[0] for r in valid])
        rats = np.array([r[3] for r in valid])
        i_mx = np.argmax(rats)
        print(f"\n  argmax α_K = {aks[i_mx]:.4f}  (ratio = {rats[i_mx]:.4f})")
        idx_r = np.argmin(np.abs(aks - AK_REF))
        print(f"  α_K_ref    = {AK_REF:.4f}  (ratio = {rats[idx_r]:.4f})")
        delta = aks[i_mx] - AK_REF
        print(f"  Δ          = {delta:.4f}")
        if abs(delta) < 0.1:
            print("  *** HIPOTEZA POTWIERDZONA (|Δ| < 0.1) ***")
        elif abs(delta) < 0.3:
            print(f"  *** CZĘŚCIOWE: |Δ|={abs(delta):.3f} ***")
        else:
            print(f"  *** HIPOTEZA OBALONA: |Δ|={abs(delta):.3f} ***")
    else:
        print(f"\n  Za mało danych ({len(valid)} punktów) do wyznaczenia argmax")

    print("="*72)
    print("KONIEC p108_v8")
