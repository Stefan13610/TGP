#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v7_narrow_tracking.py  --  TGP v1
========================================
Ścieżka 5: argmax K*₂/K*₁ — wąskie okno ψ₀ z trackingiem Branch B

DIAGNOZA BŁĘDÓW POPRZEDNICH WERSJI:
  - p108_v6 (dpsi=1.5): okno [2.46, 3.06] → ZAWIERA Branch A (ψ₀≈1.56)
    który jest bliżej k=0.090 i jest zwracany ZAMIAST Branch B (ψ₀≈2.58)!
    WYNIK: g_U = -104962 (zupełnie zły)

  - Poprawna diagnoza:
    Przy K=0.090 i α_K=8.5616:
      Branch B: ψ₀ = 2.763 - 18.78*(0.100-0.090) = 2.763 - 0.188 = 2.575
      Okno dpsi=0.30 wokół 2.763: [2.46, 3.06] → 2.575 JEST W ŚRODKU ✓
      Ale szerokie okno [1.26, 4.26] → PIERWSZE zero to Branch A (ψ₀≈1.56) ✗

  - FIX: dpsi=0.30 (wąskie) → wyklucza Branch A (ψ₀≈1.24-1.45)
    + wybór rozwiązania NAJBLIŻSZEGO centrum psi_c
    + ciągłe tracking psi_c między K-punktami

  - Zakres K w skan: [K*₂_prev * 0.85, K*₂_prev * 1.15] (±15%)
    Δψ₀ w ±15%: 18.78 × 0.15 × 0.10 = 0.28 ≈ 0.30 (graniczny)
    → Δα_K = 0.01 step → K*₂ zmienia się ≤5% → Δψ₀ ≤0.09 << 0.30 ✓

Data: 2026-03-26
"""
import sys
import io
import warnings
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
N_EVAL = 700
RTOL   = 1e-10
ATOL   = 1e-12

AK_REF  = 8.5616
K1_REF  = 0.010280;  PSI1_REF = 1.2419
K2_REF  = 0.09999348; PSI2_REF = 2.7624

# dψ₀^B/dK ≈ 18.78 w punkcie referencyjnym (z p107_v3)
DPSI_DK_B = 18.78

# ── ODE ───────────────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

_rhs_cache = {}
def get_rhs(ak):
    if ak not in _rhs_cache:
        def rhs(r, y):
            phi, dphi = y
            phi  = max(phi, 1e-10)
            kfac = 1.0 + ak/phi
            return [dphi, dV_mod(phi)/kfac + ak*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
        _rhs_cache[ak] = rhs
    return _rhs_cache[ak]


def solve_one(K, psi0, ak):
    rhs   = get_rhs(ak)
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL, t_eval=r_ev)
        if sol.t[-1] < R_MAX*0.99:
            return np.nan, np.nan
        r    = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        kfac = 1.+ak/phi
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Vv   = np.array([V_mod(float(q)) for q in phi])
        Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
        return sol.y[0,-1], Ek+Ep
    except Exception:
        return np.nan, np.nan


def find_psi0_near_center(K, ak, psi_c, dpsi=0.30, n_pts=20):
    """
    Brentq ψ₀ w [psi_c-dpsi, psi_c+dpsi] na φ(R_MAX)=1.
    Jeśli wiele zer → zwraca NAJBLIŻSZE psi_c (branch tracking).
    Zwraca (psi0, E) lub (nan, nan).
    """
    lo = max(psi_c - dpsi, 0.8)
    hi = psi_c + dpsi
    psis = np.linspace(lo, hi, n_pts)
    fv = []
    for p in psis:
        pe, _ = solve_one(K, p, ak)
        fv.append(pe - 1.0 if np.isfinite(pe) else np.nan)

    candidates = []
    for i in range(len(fv)-1):
        fi, fi1 = fv[i], fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            try:
                fn = lambda p: solve_one(K, p, ak)[0] - 1.0
                p0 = brentq(fn, psis[i], psis[i+1], xtol=1e-9, maxiter=50)
                candidates.append(p0)
            except Exception:
                t = -fi/(fi1-fi)
                candidates.append(psis[i]+t*(psis[i+1]-psis[i]))

    if not candidates:
        return np.nan, np.nan

    # Wybierz NAJBLIŻSZE psi_c
    best = min(candidates, key=lambda p: abs(p-psi_c))
    _, E = solve_one(K, best, ak)
    return best, E


def gU_tracked(K, ak, psi_c, dpsi=0.30):
    """g_U z trackingiem Branch B (wąskie okno, closest-to-center)."""
    p0, E = find_psi0_near_center(K, ak, psi_c, dpsi=dpsi)
    if not np.isfinite(p0):
        return np.nan, psi_c
    if not (np.isfinite(E) and K > 0):
        return np.nan, p0
    return E/(4*np.pi*K) - 1.0, p0


def find_K_star_tracking(ak, K_seed, psi_seed, dpsi=0.30, n_K=30, label=""):
    """
    Skan K w [K_seed*0.82, K_seed*1.18] z n_K punktami.
    ψ₀ śledzony z wąskim oknem dpsi.
    Zwraca (K_star, psi_star, gU_star).
    """
    # Oczekiwane Δψ₀ dla ΔK=±18%: 18.78 × 0.18 × K_seed
    dK_range = 0.18 * K_seed
    expected_dpsi = DPSI_DK_B * dK_range
    # Użyj dpsi = max(0.25, expected_dpsi + 0.10) ale nie mniej niż przekazany
    dpsi_use = max(dpsi, expected_dpsi + 0.08)

    K_lo = K_seed * 0.82
    K_hi = K_seed * 1.18
    K_arr = np.exp(np.linspace(np.log(K_lo), np.log(K_hi), n_K))

    gU_v   = np.full(n_K, np.nan)
    psi_v  = np.full(n_K, np.nan)
    psi_c  = psi_seed + DPSI_DK_B * (K_lo - K_seed)  # oczekiwane ψ₀ przy K_lo

    for i, K in enumerate(K_arr):
        g, p = gU_tracked(K, ak, psi_c, dpsi=dpsi_use)
        if np.isfinite(g):
            gU_v[i]  = g
            psi_v[i] = p
            psi_c    = p   # continuation
        # Jeśli nie znaleziono: aktualizuj psi_c liniowo (ekstrapolacja)
        else:
            if i > 0:
                psi_c = psi_c + DPSI_DK_B*(K_arr[i]-K_arr[i-1])

    # Szukaj zmiany znaku
    for i in range(len(gU_v)-1):
        gi, gi1 = gU_v[i], gU_v[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)): continue
        if gi*gi1 < 0:
            Ki, Ki1  = K_arr[i], K_arr[i+1]
            pm = [psi_v[i] if np.isfinite(psi_v[i]) else psi_seed]

            def fn(K_t):
                # psi_c: ekstrapoluj z Ki
                psi_est = pm[0] + DPSI_DK_B*(K_t - Ki)
                g_t, p_t = gU_tracked(K_t, ak, psi_est, dpsi=dpsi_use*0.8)
                if np.isfinite(p_t): pm[0] = p_t
                return g_t if np.isfinite(g_t) else gi

            try:
                K_star = brentq(fn, Ki, Ki1, xtol=1e-7, maxiter=50)
                psi_est = pm[0] + DPSI_DK_B*(K_star - Ki)
                g_fin, p_fin = gU_tracked(K_star, ak, psi_est, dpsi=0.25)
                return K_star, p_fin, g_fin
            except Exception:
                t = -gi/(gi1-gi)
                return Ki+t*(Ki1-Ki), psi_v[i], gi+t*(gi1-gi)

    return np.nan, psi_seed, np.nan


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    print("="*72)
    print("p108_v7: Branch B tracking z wąskim oknem dpsi=0.30")
    print(f"N_eval={N_EVAL}, RTOL={RTOL:.0e}")
    print(f"Referencja: α_K={AK_REF}, K*₂={K2_REF:.8f}, ψ₀^B={PSI2_REF:.4f}")
    print("="*72)

    # ── Weryfikacja startu ────────────────────────────────────────────────────
    print("\n[0] Diagnostyka przy α_K_ref=8.5616:")
    print(f"  Test K=0.090 (ψ₀^B oczekiwane ≈ {PSI2_REF + DPSI_DK_B*(0.090-K2_REF):.4f}):")
    K_test = 0.090
    psi_exp = PSI2_REF + DPSI_DK_B*(K_test - K2_REF)
    g_test, psi_test = gU_tracked(K_test, AK_REF, psi_exp, dpsi=0.30)
    print(f"  psi_c={psi_exp:.4f}, found ψ₀={psi_test:.4f}, g_U={g_test:.4f}")
    print(f"  (oczekiwane g_U ≈ -0.2 dla Branch B)")

    print(f"\n  Test K=0.110 (ψ₀^B oczekiwane ≈ {PSI2_REF + DPSI_DK_B*(0.110-K2_REF):.4f}):")
    K_test = 0.110
    psi_exp2 = PSI2_REF + DPSI_DK_B*(K_test - K2_REF)
    g_test2, psi_test2 = gU_tracked(K_test, AK_REF, psi_exp2, dpsi=0.30)
    print(f"  psi_c={psi_exp2:.4f}, found ψ₀={psi_test2:.4f}, g_U={g_test2:.4f}")
    print(f"  (oczekiwane g_U > 0 dla K > K*₂)")

    # Weryfikacja K*₂
    print(f"\n  Szukanie K*₂ z n_K=30 w [0.082, 0.118]:")
    K2_v, psi2_v, gU2_v = find_K_star_tracking(AK_REF, K2_REF, PSI2_REF)
    if np.isfinite(K2_v):
        delta = abs(K2_v - K2_REF)/K2_REF*100
        print(f"  K*₂ = {K2_v:.8f}  (ref={K2_REF:.8f})  g_U={gU2_v:.4e}  ψ₀={psi2_v:.4f}")
        print(f"  delta = {delta:.3f}%", "✓" if delta < 2.0 else "✗ BAD!")
        K2_start = K2_v; PSI2_start = psi2_v
    else:
        print(f"  ✗ K*₂ nie znaleziono — używam ref")
        K2_start = K2_REF; PSI2_start = PSI2_REF

    # Branch A verification
    print(f"\n  Weryfikacja Branch A:")
    K1_v, psi1_v, gU1_v = find_K_star_tracking(AK_REF, K1_REF, PSI1_REF, dpsi=0.15)
    if np.isfinite(K1_v):
        print(f"  K*₁ = {K1_v:.8f}  g_U={gU1_v:.4e}  ψ₀={psi1_v:.4f}")
        K1_start = K1_v; PSI1_start = psi1_v
    else:
        print(f"  ✗ Branch A nie znaleziono — używam ref")
        K1_start = K1_REF; PSI1_start = PSI1_REF

    # ── Descending scan ───────────────────────────────────────────────────────
    AK_desc = [AK_REF]
    ak = AK_REF
    while ak > 6.0:
        ak = round(ak - 0.01, 6)
        AK_desc.append(ak)

    # ── Ascending scan ────────────────────────────────────────────────────────
    AK_asc = [AK_REF]
    ak = AK_REF
    while ak < 8.98:
        ak = round(ak + 0.01, 6)
        AK_asc.append(ak)

    print(f"\nDescending: {len(AK_desc)} kroków (Δα_K=-0.01 do α_K=6.0)")
    print(f"Ascending:  {len(AK_asc)} kroków (Δα_K=+0.01 do α_K=8.98)")

    def run_scan(AK_values, k1_s, psi1_s, k2_s, psi2_s, label):
        results = []
        k1_c, psi1_c = k1_s, psi1_s
        k2_c, psi2_c = k2_s, psi2_s
        b2_consec_fails = 0

        print(f"\n{'='*60}")
        print(f"Skan {label}:")
        print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}  {'ψ₀^B':>8}")
        print(f"  {'-'*60}")

        for ak in AK_values:
            t0 = time.time()

            # Branch A (stable)
            K1n, psi1n, gU1n = find_K_star_tracking(ak, k1_c, psi1_c, dpsi=0.15)
            if np.isfinite(K1n) and abs(gU1n) < 0.05:
                k1_c, psi1_c = K1n, psi1n
            else:
                K1n = np.nan

            # Branch B (with tracking)
            K2n, psi2n, gU2n = np.nan, psi2_c, np.nan
            if b2_consec_fails < 5:
                K2n, psi2n, gU2n = find_K_star_tracking(ak, k2_c, psi2_c, dpsi=0.30)
                if np.isfinite(K2n) and abs(gU2n) < 0.05:
                    k2_c, psi2_c = K2n, psi2n
                    b2_consec_fails = 0
                else:
                    K2n = np.nan
                    b2_consec_fails += 1

            ratio = K2n/K1n if (np.isfinite(K1n) and np.isfinite(K2n) and K1n>0) else np.nan
            dt = time.time()-t0

            k1s = f"{K1n:.6f}" if np.isfinite(K1n) else "---"
            k2s = f"{K2n:.8f}" if np.isfinite(K2n) else "---"
            rs  = f"{ratio:.4f}" if np.isfinite(ratio) else "---"
            ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
            print(f"  {ak:>8.4f}  {k1s:>12}  {k2s:>12}  {rs:>10}  "
                  f"{psi2_c:>8.4f}  {dt:>2.0f}s{ref}")

            results.append({'ak': ak, 'K1': K1n, 'K2': K2n,
                            'psi2': psi2_c, 'ratio': ratio})

        return results

    t_total = time.time()
    res_d = run_scan(AK_desc, K1_start, PSI1_start, K2_start, PSI2_start, "MALEJĄCY")
    res_a = run_scan(AK_asc,  K1_start, PSI1_start, K2_start, PSI2_start, "ROSNĄCY")
    print(f"\nŁączny czas: {time.time()-t_total:.0f}s")

    # ── Analiza ───────────────────────────────────────────────────────────────
    all_res = {}
    for r in res_d + res_a:
        ak = round(r['ak'], 6)
        if ak not in all_res or (np.isfinite(r['ratio']) and not np.isfinite(all_res[ak]['ratio'])):
            all_res[ak] = r
    all_res = sorted(all_res.values(), key=lambda x: x['ak'])

    valid = [(r['ak'], r['K1'], r['K2'], r['ratio']) for r in all_res if np.isfinite(r['ratio'])]

    print("\n" + "="*72)
    print("TABELA:")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}")
    for ak, K1, K2, rat in valid:
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        print(f"  {ak:>8.4f}  {K1:>12.6f}  {K2:>12.6f}  {rat:>10.4f}{ref}")

    if len(valid) >= 3:
        aks  = np.array([r[0] for r in valid])
        rats = np.array([r[3] for r in valid])
        i_mx = np.argmax(rats)
        print(f"\n  argmax α_K = {aks[i_mx]:.4f}  (ratio = {rats[i_mx]:.4f})")
        idx_r = np.argmin(np.abs(aks-AK_REF))
        print(f"  α_K_ref    = {AK_REF:.4f}  (ratio = {rats[idx_r]:.4f})")
        delta = aks[i_mx] - AK_REF
        print(f"  Δ = {delta:.4f}")
        if abs(delta) < 0.05:
            print("  *** POTWIERDZONE ***")
        elif abs(delta) < 0.25:
            print(f"  *** CZĘŚCIOWE: |Δ|={abs(delta):.3f} ***")
        else:
            print(f"  *** OBALONA ***")
