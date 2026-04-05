#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v6_adaptive.py  --  TGP v1
=====================================
Ścieżka 5: argmax K*₂/K*₁ — adaptacyjna kontynuacja z szerokim oknem ψ₀

ROOT CAUSE poprzednich błędów:
  dψ₀^B/dK ≈ 18.5 przy α_K=8.5616 (z p107_v3 danych).
  Przy K odległym o 10% od K*₂, ψ₀ zmienia się o ±0.185.
  Poprzednie wersje używały okna ±0.2 → wystarczające dla 10% ale nie dalej.
  PRZY K_lo=K*₂*0.5 (50% poniżej): Δψ₀ = 18.5×0.5×K*₂ ≈ 0.93 → OUT OF ±0.2!

FIX: okno ψ₀ = [ψ₀_prev ± 1.5] z 30 punktami (spacing 0.10) na ψ₀.
Dla K w ±50% od K*₂: Δψ₀ ≤ 0.93 << 1.5 ✓

ALGORYTM KONTYNUACJI:
  1. Dla nowego α_K: oblicz g_U w K = K*₂_prev × {0.90, 1.10}
     Ψ₀ z szerokiego brentq [ψ₀_prev - 1.5, ψ₀_prev + 1.5]
  2. Jeśli zmiana znaku → brentq K w tym przedziale
  3. Jeśli brak zmiany → rozszerz do × {0.75, 1.30}, potem {0.60, 1.50}
  4. Jeśli dalej brak → Branch B martwa dla tego α_K

WERYFIKACJA STARTOWA: α_K=8.5616 musi dać K*₂=0.09999 (z p107_v3).

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


def find_psi0_wide(K, ak, psi_c, dpsi=1.5, n_pts=30):
    """
    Brentq ψ₀ w [psi_c-dpsi, psi_c+dpsi] na φ(R_MAX)=1.
    SZEROKI zakres ±1.5 by uchwycić Branch B nawet przy K daleko od K*₂.
    Zwraca (psi0, E) lub (nan, nan).
    """
    lo = max(psi_c - dpsi, 0.8)
    hi = psi_c + dpsi
    psis = np.linspace(lo, hi, n_pts)
    # Skan
    fv = []
    for p in psis:
        pe, _ = solve_one(K, p, ak)
        fv.append(pe - 1.0 if np.isfinite(pe) else np.nan)
    # Szukaj zmiany znaku
    for i in range(len(fv)-1):
        fi, fi1 = fv[i], fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            try:
                fn = lambda p: solve_one(K, p, ak)[0] - 1.0
                p0 = brentq(fn, psis[i], psis[i+1], xtol=1e-9, maxiter=50)
                _, E = solve_one(K, p0, ak)
                return p0, E
            except Exception:
                pass
    # Fallback: szukaj najbliższego zera (może wiele zmian znaku?)
    return np.nan, np.nan


def gU_at_K(K, ak, psi_c):
    """Oblicza g_U(K) z brentq ψ₀ w [psi_c±1.5]. Zwraca (g_U, psi0)."""
    p0, E = find_psi0_wide(K, ak, psi_c)
    if not np.isfinite(p0):
        return np.nan, psi_c
    if not (np.isfinite(E) and K > 0):
        return np.nan, p0
    return E/(4*np.pi*K) - 1.0, p0


def find_K_star_adaptive(ak, K_seed, psi_seed, label=""):
    """
    Znajdź K*₂ dla nowego α_K używając K_seed jako centrum.
    Adaptacyjne poszerzanie nawias K: [K_seed×f, K_seed/f] dla f=1.1,1.3,1.5,2.0
    Zwraca (K_star, psi_star, gU_star) lub (nan, nan, nan).
    """
    for factor in [1.10, 1.20, 1.35, 1.55, 2.00]:
        K_lo = K_seed / factor
        K_hi = K_seed * factor

        glo, plo = gU_at_K(K_lo, ak, psi_seed)
        ghi, phi = gU_at_K(K_hi, ak, psi_seed)

        if not (np.isfinite(glo) and np.isfinite(ghi)):
            continue

        if glo * ghi >= 0:
            # Brak zmiany znaku w tym przedziale — sprawdź czy oba po tej samej stronie
            # Może K*₂ jest poza przedziałem, spróbuj szerszy
            continue

        # Zmiana znaku znaleziona
        psi_c = [0.5*(plo+phi)]
        def fn(K_t):
            g_t, p_t = gU_at_K(K_t, ak, psi_c[0])
            if np.isfinite(p_t): psi_c[0] = p_t
            return g_t if np.isfinite(g_t) else glo

        try:
            K_star = brentq(fn, K_lo, K_hi, xtol=1e-7, maxiter=50)
            g_fin, p_fin = gU_at_K(K_star, ak, psi_c[0])
            return K_star, p_fin, g_fin
        except Exception:
            t = -glo/(ghi-glo)
            return K_lo+t*(K_hi-K_lo), 0.5*(plo+phi), glo+t*(ghi-glo)

    return np.nan, psi_seed, np.nan


# ── Weryfikacja w punkcie referencyjnym ──────────────────────────────────────

def verify_reference():
    print("[0] Weryfikacja w α_K_ref=8.5616 (K_seed=K*₂=0.09999, ψ₀=2.7624):")
    # Oblicz g_U w K_lo=0.090, K_hi=0.110 bezpośrednio
    g_lo, psi_lo = gU_at_K(0.090, AK_REF, PSI2_REF)
    g_hi, psi_hi = gU_at_K(0.110, AK_REF, PSI2_REF)
    print(f"  K=0.090: g_U={g_lo:.4f}  ψ₀={psi_lo:.4f}")
    print(f"  K=0.110: g_U={g_hi:.4f}  ψ₀={psi_hi:.4f}")

    K_v, psi_v, gU_v = find_K_star_adaptive(AK_REF, K2_REF, PSI2_REF, label="verify")
    if np.isfinite(K_v):
        print(f"  K*₂ found = {K_v:.8f}  (ref={K2_REF:.8f})  g_U={gU_v:.4e}  ψ₀={psi_v:.4f}")
        delta_K = abs(K_v - K2_REF)/K2_REF*100
        if delta_K < 1.0:
            print(f"  ✓ POPRAWNA WERYFIKACJA (delta = {delta_K:.3f}%)")
            return K_v, psi_v
        else:
            print(f"  ✗ BŁĘDNA WERYFIKACJA (delta = {delta_K:.1f}%) — używam ref")
    else:
        print(f"  ✗ Nie znaleziono K*₂ — używam ref")
    return K2_REF, PSI2_REF


# ── Skan kontynuacji ──────────────────────────────────────────────────────────

def run_scan(AK_values, k1_start, psi1_start, k2_start, psi2_start,
             also_branch_a=True, label=""):
    import time
    results = []
    k1_c, psi1_c = k1_start, psi1_start
    k2_c, psi2_c = k2_start, psi2_start
    b2_dead = False

    print(f"\n{'='*60}")
    print(f"Skan {label}: {len(AK_values)} α_K-wartości")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'g^A':>9}  {'K*₂':>12}  {'g^B':>9}  {'ratio':>8}  [s]")
    print(f"  {'-'*70}")

    for ak in AK_values:
        t0 = time.time()

        # Branch A
        K1n, psi1n, gU1n = np.nan, psi1_c, np.nan
        if also_branch_a:
            K1n, psi1n, gU1n = find_K_star_adaptive(ak, k1_c, psi1_c, label="A")
            if np.isfinite(K1n) and abs(gU1n) < 0.1:
                k1_c, psi1_c = K1n, psi1n
            else:
                K1n, psi1n, gU1n = np.nan, psi1_c, np.nan

        # Branch B
        K2n, psi2n, gU2n = np.nan, psi2_c, np.nan
        if not b2_dead:
            K2n, psi2n, gU2n = find_K_star_adaptive(ak, k2_c, psi2_c, label="B")
            if np.isfinite(K2n) and abs(gU2n) < 0.1:
                k2_c, psi2_c = K2n, psi2n
            else:
                K2n, psi2n, gU2n = np.nan, psi2_c, np.nan
                # Don't mark as dead immediately — might be temporary failure

        ratio = K2n/K1n if (np.isfinite(K1n) and np.isfinite(K2n) and K1n>0) else np.nan
        dt = time.time()-t0

        k1s = f"{K1n:.6f}" if np.isfinite(K1n) else "   ---  "
        g1s = f"{gU1n:.2e}" if np.isfinite(gU1n) else "   ---  "
        k2s = f"{K2n:.8f}" if np.isfinite(K2n) else "    ---     "
        g2s = f"{gU2n:.2e}" if np.isfinite(gU2n) else "   ---  "
        rs  = f"{ratio:.4f}" if np.isfinite(ratio) else "  ---  "
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        print(f"  {ak:>8.4f}  {k1s:>12}  {g1s:>9}  {k2s:>12}  {g2s:>9}  {rs:>8}  {dt:>3.0f}s{ref}")

        results.append({'ak': ak, 'K1': K1n, 'K2': K2n,
                        'psi1': psi1n, 'psi2': psi2n,
                        'gU1': gU1n, 'gU2': gU2n, 'ratio': ratio})

    return results


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    print("="*72)
    print("p108_v6: argmax K*₂/K*₁ — szeroke okno ψ₀ ±1.5, adaptacyjny nawias K")
    print(f"N_eval={N_EVAL}, RTOL={RTOL:.0e}")
    print(f"FIX: dψ₀/dK ≈ 18.5 → okno ±1.5 pokrywa K w ±50% od K*₂")
    print("="*72)

    # Weryfikacja
    K2_start, PSI2_start = verify_reference()
    print()

    # Branch A reference (reliable)
    print("[A_ref] Weryfikacja Branch A:")
    K1_start, PSI1_start, gU1_v = find_K_star_adaptive(AK_REF, K1_REF, PSI1_REF)
    if np.isfinite(K1_start):
        print(f"  K*₁ = {K1_start:.8f}  (ref={K1_REF:.6f})  g_U={gU1_v:.4e}")
    else:
        K1_start = K1_REF; PSI1_start = PSI1_REF
        print(f"  FALLBACK: używam K*₁_ref={K1_REF}")

    # Descending scan: od AK_REF w dół, krok 0.02
    # Nie za małe kroki bo wolno
    AK_desc = sorted([AK_REF,
                      8.54, 8.52, 8.50, 8.47, 8.44, 8.41, 8.38, 8.35,
                      8.30, 8.25, 8.20, 8.15, 8.10, 8.00, 7.90, 7.75,
                      7.50, 7.25, 7.00, 6.50, 6.00],
                     reverse=True)

    # Ascending scan: od AK_REF w górę, krok 0.02 (K*₂ szybko maleje)
    AK_asc = sorted([AK_REF,
                     8.58, 8.61, 8.64, 8.67, 8.70, 8.73, 8.76, 8.79,
                     8.82, 8.85, 8.88, 8.91, 8.94])

    print()
    t_total = time.time()

    res_desc = run_scan(AK_desc, K1_start, PSI1_start, K2_start, PSI2_start,
                        label="MALEJĄCY")
    res_asc  = run_scan(AK_asc,  K1_start, PSI1_start, K2_start, PSI2_start,
                        label="ROSNĄCY")

    print(f"\nŁączny czas: {time.time()-t_total:.0f}s = {(time.time()-t_total)/60:.1f}min")

    # ── Konsolidacja i tabela ─────────────────────────────────────────────────
    all_res_d = {}
    for r in res_desc + res_asc:
        ak = round(r['ak'], 6)
        if ak not in all_res_d or (np.isfinite(r['ratio']) and
                                    not np.isfinite(all_res_d[ak]['ratio'])):
            all_res_d[ak] = r

    all_res = sorted(all_res_d.values(), key=lambda x: x['ak'])

    print("\n" + "="*72)
    print("TABELA KOŃCOWA:")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}  {'ψ₀^B':>8}  {'g_U^B':>9}")
    print("  " + "-"*72)
    valid = []
    for r in all_res:
        ak = r['ak']
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        k1s = f"{r['K1']:.6f}" if np.isfinite(r['K1']) else "---"
        k2s = f"{r['K2']:.6f}" if np.isfinite(r['K2']) else "---"
        rs  = f"{r['ratio']:.4f}" if np.isfinite(r['ratio']) else "---"
        g2s = f"{r['gU2']:.2e}" if np.isfinite(r.get('gU2', np.nan)) else "---"
        ps2 = f"{r['psi2']:.4f}" if np.isfinite(r['psi2']) else "---"
        print(f"  {ak:>8.4f}  {k1s:>12}  {k2s:>12}  {rs:>10}  {ps2:>8}  {g2s:>9}{ref}")
        if np.isfinite(r['ratio']):
            valid.append(r)

    # ── Analiza ───────────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("ANALIZA MAKSIMUM K*₂/K*₁:")
    if len(valid) >= 3:
        aks  = np.array([r['ak']    for r in valid])
        rats = np.array([r['ratio'] for r in valid])
        i_mx = np.argmax(rats)
        ak_mx = aks[i_mx]; r_mx = rats[i_mx]
        idx_r = np.argmin(np.abs(aks-AK_REF))
        print(f"  argmax α_K = {ak_mx:.4f}  (ratio = {r_mx:.4f})")
        print(f"  α_K_ref    = {AK_REF:.4f}  (ratio = {rats[idx_r]:.4f})")
        print(f"  Δ = {ak_mx-AK_REF:.4f}")

        print("\n  Gradient (trend):")
        for i in range(len(valid)-1):
            dak = aks[i+1]-aks[i]
            dr  = rats[i+1]-rats[i]
            if abs(dak)>0.001:
                arrow = "↑" if dr>0 else "↓"
                print(f"    [{aks[i]:.3f}→{aks[i+1]:.3f}]: {arrow} Δratio={dr:+.3f}  grad={dr/dak:+.2f}")

        print("\n  Top-5:")
        for i in np.argsort(rats)[::-1][:5]:
            ref = " <REF" if abs(valid[i]['ak']-AK_REF)<1e-3 else ""
            print(f"    α_K={valid[i]['ak']:.4f}  ratio={rats[i]:.4f}  "
                  f"K*₁={valid[i]['K1']:.6f}  K*₂={valid[i]['K2']:.6f}{ref}")

        delta = ak_mx - AK_REF
        print(f"\n  WNIOSEK ŚCIEŻKA 5:")
        if abs(delta) < 0.05:
            print(f"  *** POTWIERDZONE: argmax={ak_mx:.4f} ≈ α_K_ref *** ")
        elif abs(delta) < 0.25:
            print(f"  *** CZĘŚCIOWE: |Δ|={abs(delta):.3f} ***")
        else:
            print(f"  *** OBALONA: argmax={ak_mx:.4f}, Δ={delta:.3f} ***")
    else:
        print("  Za mało punktów do analizy.")
