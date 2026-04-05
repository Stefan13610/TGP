#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v3_continuation.py  --  TGP v1
=====================================
Ścieżka 5 OP-3: argmax K*₂/K*₁ — continuation z małymi krokami Δα_K

Problem poprzednich wersji: zbyt szeroki skan K i ψ₀ → fałszywe minima g_U.

Rozwiązanie: dokładnie to samo co p107_v3 (skan K z brentq ψ₀), ale
powtarzane dla różnych α_K z WĄSKIM zakresem wokół poprzedniej wartości.

Metoda continuation:
  - Start: α_K=8.5616, K*₁=0.010280, ψ₀^A=1.2419, K*₂=0.09999, ψ₀^B=2.7624
  - Descending: α_K=8.5616→6.0 krokiem Δα_K=-0.10 (K*₂ wolno rośnie/maleje)
  - Ascending: α_K=8.5616→8.99 krokiem Δα_K=+0.02 (K*₂ szybko maleje!)
    → małe kroki bo blisko bifurkacji K*₂ → K*₁

Dla każdego α_K:
  1) Skan K^A w [K1_prev*0.5, K1_prev*3], ψ₀^A w [ψ1_prev-0.15, ψ1_prev+0.15]
  2) Skan K^B w [K2_prev*0.5, K2_prev*2.5], ψ₀^B w [ψ2_prev-0.20, ψ2_prev+0.20]
  3) Weryfikacja: g_U < threshold przy znalezionym K*

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
N_EVAL = 1200
RTOL   = 1e-10
ATOL   = 1e-12

AK_REF   = 8.5616
K1_REF   = 0.010280;  PSI1_REF = 1.2419
K2_REF   = 0.09999348; PSI2_REF = 2.7624
ALPHA_SN = 8.9964
GU_TOL   = 0.10   # g_U < GU_TOL do akceptacji K*

# ── ODE ───────────────────────────────────────────────────────────────────────

def V_mod(phi): return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def make_rhs(ak):
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
        kfac = 1.0 + ak/phi
        return [dphi,
                dV_mod(phi)/kfac + ak*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
    return rhs

_rhs_cache = {}
def get_rhs(ak):
    key = round(ak, 8)
    if key not in _rhs_cache:
        _rhs_cache[key] = make_rhs(ak)
    return _rhs_cache[key]


def solve_one(K, psi0, ak):
    """Integruje, zwraca (phi_end, E) lub (nan, nan)."""
    rhs   = get_rhs(ak)
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL, t_eval=r_ev)
        if sol.t[-1] < R_MAX*0.99:
            return np.nan, np.nan
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        kfac = 1.+ak/phi
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Vv   = np.array([V_mod(float(q)) for q in phi])
        Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
        return sol.y[0,-1], Ek+Ep
    except Exception:
        return np.nan, np.nan


def find_psi0_narrow(K, ak, psi_c, dpsi=0.18, n_pts=18):
    """Brentq ψ₀ w [psi_c-dpsi, psi_c+dpsi] na φ(r_max)=1."""
    lo = max(psi_c - dpsi, 0.5)
    hi = psi_c + dpsi
    psis = np.linspace(lo, hi, n_pts)
    Fv   = []
    for p in psis:
        pe, _ = solve_one(K, p, ak)
        Fv.append(pe - 1.0 if np.isfinite(pe) else np.nan)
    for i in range(len(Fv)-1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)):
            continue
        if fi*fi1 < 0:
            try:
                fn = lambda p: (solve_one(K, p, ak)[0] - 1.0)
                p0 = brentq(fn, psis[i], psis[i+1], xtol=1e-9, maxiter=60)
                return p0
            except Exception:
                pass
    return np.nan


def gU_single(K, ak, psi_c, dpsi=0.18):
    """g_U dla jednego K, z refinement ψ₀."""
    p0 = find_psi0_narrow(K, ak, psi_c, dpsi=dpsi)
    if not np.isfinite(p0):
        return np.nan, np.nan
    _, E = solve_one(K, p0, ak)
    if not np.isfinite(E):
        return np.nan, p0
    return E/(4*np.pi*K) - 1.0, p0


def find_K_star(ak, K_lo, K_hi, psi_c, dpsi=0.18, n_K=20, label=""):
    """
    Skan g_U na [K_lo,K_hi] z n_K punktami log-spaced.
    ψ₀ śledzone metodą continuation (wąskie okno wokół psi_c).
    Zwraca (K_star, psi_star, gU_star) lub (nan,nan,nan).
    """
    K_arr = np.exp(np.linspace(np.log(K_lo), np.log(K_hi), n_K))
    gU_v  = np.full(n_K, np.nan)
    psi_v = np.full(n_K, np.nan)
    psi_curr = psi_c

    for i, K in enumerate(K_arr):
        g, p = gU_single(K, ak, psi_curr, dpsi=dpsi)
        if np.isfinite(g):
            gU_v[i]  = g
            psi_v[i] = p
            psi_curr = p  # continuation
        # Else: psi_curr remains unchanged → narrow bracket might find next K

    # Znajdź zmianę znaku g_U
    for i in range(len(gU_v)-1):
        gi, gi1 = gU_v[i], gU_v[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)):
            continue
        if gi*gi1 < 0:
            Ki, Ki1  = K_arr[i], K_arr[i+1]
            pm = psi_v[i] if np.isfinite(psi_v[i]) else psi_c
            pm = [pm]

            def fn(K_t):
                g_t, p_t = gU_single(K_t, ak, pm[0], dpsi=dpsi*0.8)
                if np.isfinite(p_t): pm[0] = p_t
                return g_t if np.isfinite(g_t) else gi

            try:
                K_star = brentq(fn, Ki, Ki1, xtol=1e-7, maxiter=50)
                g_fin, p_fin = gU_single(K_star, ak, pm[0], dpsi=dpsi*0.6)
                if np.isfinite(g_fin) and np.isfinite(p_fin):
                    return K_star, p_fin, g_fin
                return K_star, pm[0], g_fin
            except Exception:
                t    = -gi/(gi1-gi)
                Klin = Ki + t*(Ki1-Ki)
                return Klin, psi_v[i], gi+t*(gi1-gi)

    return np.nan, psi_curr, np.nan


# ── Continuation scan ─────────────────────────────────────────────────────────

def run_continuation(AK_values, k1_start, psi1_start, k2_start, psi2_start,
                     label=""):
    """
    Skan listy AK_values w podanej kolejności.
    Używa continuation dla K*₁ i K*₂.
    """
    import time
    results = []
    k1_c, psi1_c = k1_start, psi1_start
    k2_c, psi2_c = k2_start, psi2_start
    b2_alive = True

    print(f"\n{'='*60}")
    print(f"Skan {label}: {len(AK_values)} punktów")
    print(f"  {'α_K':>8}  {'K*₁':>10}  {'g^A':>9}  {'K*₂':>12}  {'g^B':>9}  {'ratio':>8}  [s]")
    print(f"  {'-'*70}")

    for ak in AK_values:
        t0 = time.time()

        # ── Branch A ──────────────────────────────────────────────────────────
        K1_lo = max(k1_c*0.4, 5e-4)
        K1_hi = k1_c * 3.5
        K1n, psi1n, gU1n = find_K_star(ak, K1_lo, K1_hi, psi1_c,
                                        dpsi=0.18, n_K=22, label="A")
        if np.isfinite(K1n) and abs(gU1n) < GU_TOL:
            k1_c, psi1_c = K1n, psi1n
        else:
            K1n, psi1n, gU1n = np.nan, psi1_c, np.nan

        # ── Branch B ──────────────────────────────────────────────────────────
        K2n, psi2n, gU2n = np.nan, psi2_c, np.nan
        if b2_alive:
            K2_lo = max(k2_c*0.5, 0.002)
            K2_hi = min(k2_c*2.2, 0.50)
            K2n, psi2n, gU2n = find_K_star(ak, K2_lo, K2_hi, psi2_c,
                                            dpsi=0.22, n_K=25, label="B")
            if np.isfinite(K2n) and abs(gU2n) < GU_TOL:
                k2_c, psi2_c = K2n, psi2n
            else:
                K2n, psi2n, gU2n = np.nan, psi2_c, np.nan
                # Branch B died?
                # Mark as dead only if it was already marginal
                if not np.isfinite(K2n):
                    pass  # Try again next step with same seed

        ratio = K2n/K1n if (np.isfinite(K1n) and np.isfinite(K2n) and K1n>0) else np.nan
        dt = time.time() - t0

        # Display
        k1s = f"{K1n:.6f}" if np.isfinite(K1n) else "   ---  "
        g1s = f"{gU1n:.2e}" if np.isfinite(gU1n) else "   --- "
        k2s = f"{K2n:.8f}" if np.isfinite(K2n) else "     ---    "
        g2s = f"{gU2n:.2e}" if np.isfinite(gU2n) else "   --- "
        rs  = f"{ratio:.4f}" if np.isfinite(ratio) else "  ---  "
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        print(f"  {ak:>8.4f}  {k1s:>10}  {g1s:>9}  {k2s:>12}  {g2s:>9}  {rs:>8}  {dt:>4.0f}s{ref}")

        results.append({
            'ak': ak, 'K1': K1n, 'psi1': psi1n, 'gU1': gU1n,
            'K2': K2n, 'psi2': psi2n, 'gU2': gU2n, 'ratio': ratio
        })

    return results


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    print("="*72)
    print("p108_v3: argmax K*₂/K*₁ via continuation z wąskim oknem")
    print(f"N_eval={N_EVAL}, RTOL={RTOL:.0e}, g_U_tol={GU_TOL}")
    print(f"Referencja: α_K={AK_REF}, K*₁={K1_REF:.6f}, K*₂={K2_REF:.6f}")
    print("="*72)

    # ── Weryfikacja w punkcie referencyjnym ───────────────────────────────────
    print("\n[0] Weryfikacja w α_K_ref:")
    K1v, psi1v, gU1v = find_K_star(AK_REF, K1_REF*0.4, K1_REF*3, PSI1_REF,
                                    dpsi=0.18, n_K=22, label="ref_A")
    K2v, psi2v, gU2v = find_K_star(AK_REF, K2_REF*0.5, K2_REF*2.2, PSI2_REF,
                                    dpsi=0.22, n_K=25, label="ref_B")
    print(f"  Branch A: K*₁={K1v:.8f} (ref={K1_REF:.6f}), g_U={gU1v:.4e}")
    print(f"  Branch B: K*₂={K2v:.8f} (ref={K2_REF:.8f}), g_U={gU2v:.4e}")
    if np.isfinite(K1v) and abs(gU1v) < GU_TOL:
        K1_REF2, PSI1_REF2 = K1v, psi1v
    else:
        print("  UWAGA: Branch A nie potwierdzona! Używam ref.")
        K1_REF2, PSI1_REF2 = K1_REF, PSI1_REF
    if np.isfinite(K2v) and abs(gU2v) < GU_TOL:
        K2_REF2, PSI2_REF2 = K2v, psi2v
    else:
        print("  UWAGA: Branch B nie potwierdzona! Używam ref.")
        K2_REF2, PSI2_REF2 = K2_REF, PSI2_REF

    # ── Descending scan: α_K < AK_REF ────────────────────────────────────────
    # Krok 0.10 bo K*₂ zmienia się wolno (hipoteza)
    AK_desc = sorted([
        AK_REF,
        8.50, 8.45, 8.40, 8.35, 8.30, 8.25, 8.20, 8.15,
        8.10, 8.00, 7.90, 7.75, 7.50, 7.25, 7.00,
        6.75, 6.50, 6.00
    ], reverse=True)

    res_desc = run_continuation(AK_desc, K1_REF2, PSI1_REF2, K2_REF2, PSI2_REF2,
                                label="MALEJĄCY (α_K < ref)")

    # ── Ascending scan: α_K > AK_REF ─────────────────────────────────────────
    # Małe kroki 0.02 bo K*₂ szybko maleje blisko α_SN
    AK_asc = sorted([
        AK_REF,
        8.58, 8.60, 8.62, 8.65, 8.68, 8.70, 8.73, 8.75,
        8.78, 8.80, 8.83, 8.85, 8.88, 8.90, 8.92, 8.95
    ])

    res_asc = run_continuation(AK_asc, K1_REF2, PSI1_REF2, K2_REF2, PSI2_REF2,
                               label="ROSNĄCY (α_K > ref)")

    # ── Konsolidacja ──────────────────────────────────────────────────────────
    all_res = {}
    for r in res_desc + res_asc:
        ak = round(r['ak'], 6)
        if ak not in all_res or (np.isfinite(r['ratio']) and not np.isfinite(all_res[ak]['ratio'])):
            all_res[ak] = r
    all_res = sorted(all_res.values(), key=lambda x: x['ak'])

    # ── Tabela końcowa ────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("TABELA KOŃCOWA:")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}  {'g_U^A':>9}  {'g_U^B':>9}")
    print("  " + "-"*72)
    valid = []
    for r in all_res:
        ak = r['ak']
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        k1s = f"{r['K1']:.6f}" if np.isfinite(r['K1']) else "---"
        k2s = f"{r['K2']:.6f}" if np.isfinite(r['K2']) else "---"
        g1s = f"{r['gU1']:.2e}" if np.isfinite(r['gU1']) else "---"
        g2s = f"{r['gU2']:.2e}" if np.isfinite(r['gU2']) else "---"
        rs  = f"{r['ratio']:.4f}" if np.isfinite(r['ratio']) else "---"
        print(f"  {ak:>8.4f}  {k1s:>12}  {k2s:>12}  {rs:>10}  {g1s:>9}  {g2s:>9}{ref}")
        if np.isfinite(r['ratio']):
            valid.append(r)

    # ── Analiza ───────────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("ANALIZA:")
    if len(valid) >= 3:
        aks  = np.array([r['ak']   for r in valid])
        rats = np.array([r['ratio'] for r in valid])
        i_mx = np.argmax(rats)
        print(f"\n  argmax α_K = {aks[i_mx]:.4f}  (ratio = {rats[i_mx]:.4f})")
        print(f"  α_K_ref    = {AK_REF:.4f}  ", end="")
        idx_ref = np.argmin(np.abs(aks - AK_REF))
        print(f"(ratio = {rats[idx_ref]:.4f})")
        print(f"  Δ(argmax − ref) = {aks[i_mx] - AK_REF:.4f}")

        print("\n  Gradient (Δratio/Δα_K):")
        for i in range(len(valid)-1):
            dak = aks[i+1]-aks[i]
            dr  = rats[i+1]-rats[i]
            if abs(dak) > 1e-4:
                sign = "↑" if dr > 0 else "↓"
                print(f"    [{aks[i]:.3f}→{aks[i+1]:.3f}]: {sign}  Δratio={dr:+.3f}  grad={dr/dak:+.2f}")

        print("\n  Top-5 punktów wg ratio:")
        for i in np.argsort(rats)[::-1][:5]:
            ref = " <-- REF" if abs(valid[i]['ak']-AK_REF)<1e-3 else ""
            print(f"    α_K={valid[i]['ak']:.4f}  ratio={rats[i]:.4f}{ref}")

        delta = aks[i_mx] - AK_REF
        print(f"\n  WNIOSEK ŚCIEŻKA 5:")
        if abs(delta) < 0.08:
            print(f"  *** POTWIERDZONE: argmax ≈ α_K_ref (delta={delta:.3f}) ***")
        elif abs(delta) < 0.25:
            print(f"  *** CZĘŚCIOWE: argmax = {aks[i_mx]:.4f}, |delta|={abs(delta):.3f} ***")
        else:
            print(f"  *** OBALONA: argmax = {aks[i_mx]:.4f}, delta = {delta:.3f} ***")
    else:
        print("  Za mało punktów do analizy.")
