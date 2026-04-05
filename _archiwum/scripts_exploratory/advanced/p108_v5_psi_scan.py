#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v5_psi_scan.py  --  TGP v1
=====================================
Ścieżka 5 OP-3: argmax K*₂/K*₁ — skan ψ₀ zamiast K

KLUCZ: zamiast skanować K i szukać ψ₀ dla każdego K (problem: wiele
ψ₀-rozwiązań dla jednego K), skanujemy ψ₀ i szukamy K_bc(ψ₀) takie że
φ(R_MAX)=1. Wtedy g_U(ψ₀) = E/(4π·K_bc) − 1 to gładka krzywa 1D.

Dlaczego to działa:
- Dla danego ψ₀ > 1 i dużego α_K, ISTNIEJE dokładnie jeden K_bc (φ maleje z
  K, więc jest monotone i brentq zawsze znajdzie go w [K_lo, K_hi])
- Branch A: zero g_U(ψ₀) przy ψ₀ ≈ 1.24 (K_bc ≈ 0.01)
- Branch B (muon): zero g_U(ψ₀) przy ψ₀ ≈ 2.76 (K_bc ≈ 0.10) — ale
  dla α_K < 8.5616 to ψ₀ przesuwa się prawdopodobnie > 2.76

Zakres ψ₀: [1.0, 8.0] z 90 punktami → spacing ≈ 0.08
Dla każdego ψ₀: brentq K_bc w [0.001, 0.50], xtol=1e-7

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
N_EVAL = 900
RTOL   = 1e-10
ATOL   = 1e-12

AK_REF = 8.5616
K1_REF  = 0.010280
K2_REF  = 0.09999348
PSI1_REF = 1.2419
PSI2_REF = 2.7624

# ── ODE ───────────────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

_rhs_cache = {}
def get_rhs(ak):
    key = round(ak, 8)
    if key not in _rhs_cache:
        def rhs(r, y):
            phi, dphi = y
            phi  = max(phi, 1e-10)
            kfac = 1.0 + ak/phi
            return [dphi, dV_mod(phi)/kfac + ak*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
        _rhs_cache[key] = rhs
    return _rhs_cache[key]


def integrate_full(K, psi0, ak):
    """Zwraca (phi_end, E) lub (nan, nan)."""
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


def find_Kbc(psi0, ak, K_lo=5e-4, K_hi=0.50, n_pts=14):
    """
    Znajdź K_bc taki że φ(R_MAX; K_bc, ψ₀, α_K) = 1.
    K_lo → φ_end > 1 (małe K, pole "przelatuje nad 1")
    K_hi → φ_end < 1 (duże K, strome opadanie) — MUSI być spełnione
    Zwraca (K_bc, E) lub (nan, nan).
    """
    # Sprawdź granice
    phi_lo, _ = integrate_full(K_lo, psi0, ak)
    phi_hi, _ = integrate_full(K_hi, psi0, ak)

    if not (np.isfinite(phi_lo) and np.isfinite(phi_hi)):
        # Spróbuj innego zakresu
        for K_hi_try in [0.30, 0.20, 0.15, 0.10]:
            phi_hi, _ = integrate_full(K_hi_try, psi0, ak)
            if np.isfinite(phi_hi):
                K_hi = K_hi_try
                break
        if not np.isfinite(phi_hi):
            return np.nan, np.nan

    f_lo = phi_lo - 1.0
    f_hi = phi_hi - 1.0

    if not (np.isfinite(f_lo) and np.isfinite(f_hi)):
        return np.nan, np.nan

    # Dla ψ₀ > 1: f_lo > 0 typowo (mało K → duże φ_end) i f_hi < 0
    # Ale dla ψ₀ bardzo dużych może być inaczej
    if f_lo * f_hi >= 0:
        # Poszukaj zmiany znaku w gęstszej siatce
        Kv = np.exp(np.linspace(np.log(K_lo), np.log(K_hi), n_pts))
        fv = [integrate_full(K, psi0, ak)[0]-1.0 for K in Kv]
        for i in range(len(fv)-1):
            fi, fi1 = fv[i], fv[i+1]
            if np.isfinite(fi) and np.isfinite(fi1) and fi*fi1 < 0:
                f_lo, f_hi = fi, fi1
                K_lo, K_hi = Kv[i], Kv[i+1]
                break
        else:
            return np.nan, np.nan

    try:
        fn = lambda K: (integrate_full(K, psi0, ak)[0] - 1.0)
        K_bc = brentq(fn, K_lo, K_hi, xtol=1e-8, maxiter=60)
        _, E = integrate_full(K_bc, psi0, ak)
        return K_bc, E
    except Exception:
        return np.nan, np.nan


def gU_vs_psi(ak, psi_lo=1.05, psi_hi=8.0, n_psi=90):
    """
    Skanuje ψ₀ ∈ [psi_lo, psi_hi] i oblicza:
      K_bc(ψ₀), g_U(ψ₀) = E/(4π·K_bc) - 1
    Zwraca tablice (psi_arr, Kbc_arr, gU_arr).
    """
    psi_arr = np.linspace(psi_lo, psi_hi, n_psi)
    Kbc_arr = np.full(n_psi, np.nan)
    gU_arr  = np.full(n_psi, np.nan)
    E_arr   = np.full(n_psi, np.nan)

    for i, psi0 in enumerate(psi_arr):
        K_bc, E = find_Kbc(psi0, ak)
        if np.isfinite(K_bc) and np.isfinite(E) and K_bc > 0:
            Kbc_arr[i] = K_bc
            E_arr[i]   = E
            gU_arr[i]  = E/(4*np.pi*K_bc) - 1.0

    return psi_arr, Kbc_arr, gU_arr


def find_zeros(psi_arr, Kbc_arr, gU_arr, ak):
    """
    Znajdź wszystkie zera g_U(ψ₀) poprzez zmianę znaku.
    Zwraca listę (psi_star, K_star, gU_star).
    """
    zeros = []
    for i in range(len(gU_arr)-1):
        gi, gi1 = gU_arr[i], gU_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)): continue
        if gi*gi1 < 0:
            plo, phi = psi_arr[i], psi_arr[i+1]
            try:
                def fn(psi):
                    Kbc, E = find_Kbc(psi, ak)
                    if not (np.isfinite(Kbc) and np.isfinite(E)): return gi
                    return E/(4*np.pi*Kbc) - 1.0
                psi_star = brentq(fn, plo, phi, xtol=1e-7, maxiter=50)
                K_star, E_star = find_Kbc(psi_star, ak)
                gU_star = E_star/(4*np.pi*K_star) - 1.0 if np.isfinite(E_star) else np.nan
                zeros.append((psi_star, K_star, gU_star))
            except Exception:
                t = -gi/(gi1-gi)
                psi_lin = plo + t*(phi-plo)
                Klin = Kbc_arr[i] + t*(Kbc_arr[i+1]-Kbc_arr[i]) if np.isfinite(Kbc_arr[i+1]) else Kbc_arr[i]
                zeros.append((psi_lin, Klin, gi+t*(gi1-gi)))
    return zeros


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    print("="*72)
    print("p108_v5: g_U vs ψ₀ skan — K_bc(ψ₀) z brentq")
    print(f"R_MAX={R_MAX}, N_eval={N_EVAL}, RTOL={RTOL:.0e}")
    print(f"Skan ψ₀ ∈ [1.05, 8.0] z 90 punktami per α_K")
    print(f"Referencja: α_K={AK_REF}, K*₁={K1_REF:.6f}, K*₂={K2_REF:.8f}")
    print("="*72)

    # α_K do przeskanowania
    AK_LIST = [6.0, 7.0, 7.5, 8.0, 8.2, 8.3, 8.4, 8.45, 8.5,
               AK_REF, 8.6, 8.65, 8.68, 8.80, 8.90, 8.95]

    all_results = {}
    total_t0 = time.time()

    for ak in AK_LIST:
        t0 = time.time()
        print(f"\n{'─'*60}")
        print(f"α_K = {ak:.4f} — skan ψ₀...", end="", flush=True)

        psi_arr, Kbc_arr, gU_arr = gU_vs_psi(ak)
        dt = time.time() - t0
        n_valid = np.sum(np.isfinite(gU_arr))
        print(f" {dt:.0f}s, {n_valid} ważnych ψ₀")

        # Pokaż krzywą g_U vs ψ₀ (wybrane punkty)
        valid = np.isfinite(gU_arr) & np.isfinite(Kbc_arr)
        if np.any(valid):
            pv = psi_arr[valid]; Kv = Kbc_arr[valid]; gv = gU_arr[valid]
            print(f"  K_bc range: [{Kv.min():.5f}, {Kv.max():.5f}]")
            print(f"  g_U range:  [{gv.min():.3f}, {gv.max():.3f}]")
            # Pokaż kilka punktów co 10
            step = max(1, len(pv)//8)
            for i in range(0, len(pv), step):
                print(f"    ψ₀={pv[i]:.3f}  K_bc={Kv[i]:.6f}  g_U={gv[i]:.4f}")

        # Znajdź zera
        zeros = find_zeros(psi_arr, Kbc_arr, gU_arr, ak)
        print(f"  Znaleziono {len(zeros)} zer g_U:")
        for z in zeros:
            psi_z, K_z, gU_z = z
            print(f"    ψ₀*={psi_z:.6f}  K*={K_z:.8f}  g_U={gU_z:.4e}")

        all_results[ak] = {'zeros': zeros, 'psi_arr': psi_arr,
                           'Kbc_arr': Kbc_arr, 'gU_arr': gU_arr}

    total_dt = time.time() - total_t0
    print(f"\n{'='*72}")
    print(f"Łączny czas: {total_dt:.0f}s = {total_dt/60:.1f} min")

    # ── Tabela zbiorowa ───────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("TABELA: K*₁, K*₂, ratio vs α_K")
    print(f"  {'α_K':>8}  {'K*₁ (Branch A)':>16}  {'K*₂ (Branch B)':>16}  {'ratio':>10}")
    print("  " + "-"*60)
    valid_results = []
    for ak in AK_LIST:
        res = all_results[ak]
        zeros = res['zeros']
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""

        if not zeros:
            print(f"  {ak:>8.4f}  {'---':>16}  {'---':>16}  {'---':>10}  [brak zer]{ref}")
            continue

        # Sortuj zera po K_bc (rosnąco)
        zeros_sorted = sorted(zeros, key=lambda z: z[1] if np.isfinite(z[1]) else 1e9)

        # Branch A = zero z najmniejszym K (ψ₀ ≈ 1.2)
        # Branch B = zero z następnym K (ψ₀ ≈ 2.76 przy ref)
        K1 = zeros_sorted[0][1] if zeros_sorted else np.nan
        K2 = zeros_sorted[1][1] if len(zeros_sorted) >= 2 else np.nan
        psi1 = zeros_sorted[0][0] if zeros_sorted else np.nan
        psi2 = zeros_sorted[1][0] if len(zeros_sorted) >= 2 else np.nan
        gU1  = zeros_sorted[0][2] if zeros_sorted else np.nan
        gU2  = zeros_sorted[1][2] if len(zeros_sorted) >= 2 else np.nan
        ratio = K2/K1 if (np.isfinite(K1) and np.isfinite(K2) and K1>0) else np.nan

        # Szczegóły każdego zera
        for j, (pz, Kz, gUz) in enumerate(zeros_sorted):
            lab = f"zer_{j}"
            print(f"    [{lab}] ψ₀*={pz:.5f} K*={Kz:.6f} g_U={gUz:.3e}")

        rat_s = f"{ratio:.4f}" if np.isfinite(ratio) else "---"
        k1s = f"{K1:.6f}(ψ₀={psi1:.4f})" if np.isfinite(K1) else "---"
        k2s = f"{K2:.6f}(ψ₀={psi2:.4f})" if np.isfinite(K2) else "---"
        print(f"  {ak:>8.4f}  {k1s:>16}  {k2s:>16}  {rat_s:>10}{ref}")

        if np.isfinite(ratio):
            valid_results.append({'ak': ak, 'K1': K1, 'K2': K2,
                                   'psi1': psi1, 'psi2': psi2, 'ratio': ratio})

    # ── Finalna tabela ────────────────────────────────────────────────────────
    print("\nFINALNA TABELA (tylko kompletne wpisy):")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'ψ₀^A':>8}  {'K*₂':>12}  {'ψ₀^B':>8}  {'ratio':>10}")
    print("  " + "-"*72)
    for r in valid_results:
        ref = " <REF" if abs(r['ak']-AK_REF)<1e-3 else ""
        print(f"  {r['ak']:>8.4f}  {r['K1']:>12.6f}  {r['psi1']:>8.4f}  "
              f"{r['K2']:>12.6f}  {r['psi2']:>8.4f}  {r['ratio']:>10.4f}{ref}")

    # ── Analiza maksimum ──────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("ANALIZA MAKSIMUM K*₂/K*₁:")
    if len(valid_results) >= 3:
        aks  = np.array([r['ak']    for r in valid_results])
        rats = np.array([r['ratio'] for r in valid_results])
        i_mx = np.argmax(rats)
        ak_mx, r_mx = aks[i_mx], rats[i_mx]
        print(f"  argmax α_K  = {ak_mx:.4f}  (ratio = {r_mx:.4f})")
        idx_ref = np.argmin(np.abs(aks-AK_REF))
        print(f"  α_K_ref     = {AK_REF:.4f}  (ratio = {rats[idx_ref]:.4f})")
        delta = ak_mx - AK_REF
        print(f"  Δ = {delta:.4f}  ({abs(delta)/AK_REF*100:.2f}%)")

        print("\n  Pełna lista (ratio rosnąco):")
        for i in np.argsort(rats):
            ref = " <REF" if abs(valid_results[i]['ak']-AK_REF)<1e-3 else ""
            print(f"    α_K={valid_results[i]['ak']:.4f}  ratio={rats[i]:.4f}  "
                  f"K*₁={valid_results[i]['K1']:.6f}  K*₂={valid_results[i]['K2']:.6f}{ref}")

        print(f"\n  WNIOSEK ŚCIEŻKA 5:")
        if abs(delta) < 0.05:
            print(f"  *** POTWIERDZONE: argmax α_K={ak_mx:.4f} ≈ α_K_ref={AK_REF} ***")
        elif abs(delta) < 0.20:
            print(f"  *** CZĘŚCIOWE: argmax={ak_mx:.4f}, |Δ|={abs(delta):.3f} ***")
        else:
            print(f"  *** OBALONA: argmax={ak_mx:.4f}, Δ={delta:.4f} ***")
    else:
        print("  Za mało punktów do analizy!")
