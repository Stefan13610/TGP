#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v2_argmax_clean.py  --  TGP v1
=====================================
Ścieżka 5 OP-3: argmax K*₂/K*₁ — czysty skan

Problem p108_v1: continuation K*₂ zawodzi bo K*₂ spada 10× dla Δα_K=0.4,
a skan Branch A ma błędne zakresy.

Poprawka: dla każdego α_K szukaj K*₁ i K*₂ NIEZALEŻNIE przez:
  1) Skan K na gęstej siatce log-spacedowej [K_lo, K_hi]
  2) Dla każdego K: znajdź ψ₀ spełniające φ(R_MAX)=1 OSOBNO dla Branch A
     (ψ₀ ∈ [0.9, 2.0]) i Branch B (ψ₀ ∈ [1.9, 4.5])
  3) Oblicz g_U(K) dla obu gałęzi, znajdź zmianę znaku → K*

Kluczowe przyspieszenie: N_EVAL=700, RTOL=1e-9 (wystarczy dla scan),
                         brentq z małym n_scan=12 per ψ₀-range

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
GAM      = 1.0
V1       = GAM/3.0 - GAM/4.0
LAM      = 5.501357e-06
AGAM     = 0.040
R_MAX    = 80.0
N_EVAL   = 700        # szybki scan
RTOL     = 1e-9
ATOL     = 1e-11

# Referencja (p107_v3)
AK_REF   = 8.5616
K1_REF   = 0.010280
K2_REF   = 0.09999348
PSI1_REF = 1.2419
PSI2_REF = 2.7624
ALPHA_SN = 8.9964

# ── ODE ───────────────────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6

def dV_mod(phi):
    return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def make_rhs(ak):
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
        kfac = 1.0 + ak/phi
        return [dphi,
                dV_mod(phi)/kfac + ak*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
    return rhs


def integrate(K, psi0, ak):
    """Integruje ODE, zwraca (phi_end, E) lub (nan, nan)."""
    rhs  = make_rhs(ak)
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    for method in ('DOP853',):
        try:
            sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                            method=method, rtol=RTOL, atol=ATOL, t_eval=r_ev)
            if sol.t[-1] < R_MAX*0.99:
                continue
            r    = sol.t
            phi  = np.maximum(sol.y[0], 1e-10)
            dphi = sol.y[1]
            kfac = 1.+ak/phi
            Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
            Vv   = np.array([V_mod(float(q)) for q in phi])
            Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
            return sol.y[0,-1], Ek+Ep
        except Exception:
            pass
    return np.nan, np.nan


def find_psi0(K, ak, psi_lo, psi_hi, n_pts=16):
    """
    Brentq ψ₀ w zakresie [psi_lo, psi_hi] tak że φ(R_MAX)=1.
    Zwraca (psi0, phi_end) lub (nan, nan).
    """
    psis = np.linspace(psi_lo, psi_hi, n_pts)
    fvals = []
    for p in psis:
        pe, _ = integrate(K, p, ak)
        fvals.append(pe - 1.0 if np.isfinite(pe) else np.nan)

    # Szukaj zmiany znaku
    for i in range(len(fvals)-1):
        fi, fi1 = fvals[i], fvals[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)):
            continue
        if fi*fi1 < 0:
            try:
                def fn(psi):
                    pe, _ = integrate(K, psi, ak)
                    return pe - 1.0 if np.isfinite(pe) else fi
                psi_star = brentq(fn, psis[i], psis[i+1], xtol=1e-8, maxiter=50)
                pe_final, E = integrate(K, psi_star, ak)
                return psi_star, pe_final, E
            except Exception:
                pass
    return np.nan, np.nan, np.nan


def gU_from_E(K, E):
    """g_U = E/(4πK) − 1"""
    if np.isfinite(E) and K > 0:
        return E/(4*np.pi*K) - 1.0
    return np.nan


def find_K_star_branch(ak, K_lo, K_hi, psi_lo, psi_hi, n_K=25, label=""):
    """
    Szuka K* dla gałęzi A lub B:
    - Skanuje K w [K_lo, K_hi] z n_K punktami
    - Dla każdego K: find_psi0 w [psi_lo, psi_hi]
    - Oblicza g_U(K), szuka zmiany znaku → brentq
    Zwraca (K_star, psi0_star, gU_star) lub (nan, nan, nan).
    """
    K_arr  = np.exp(np.linspace(np.log(K_lo), np.log(K_hi), n_K))
    gU_arr = np.full(n_K, np.nan)
    psi_arr = np.full(n_K, np.nan)

    # Ciągłe ψ₀ — śledź między K-punktami
    psi_c = 0.5*(psi_lo + psi_hi)

    for i, K in enumerate(K_arr):
        # Szukaj ψ₀ adaptacyjnie wokół poprzedniej wartości
        lo_ = max(psi_c - 0.25, psi_lo)
        hi_ = min(psi_c + 0.25, psi_hi)
        p0, pe, E = find_psi0(K, ak, lo_, hi_, n_pts=14)
        if not np.isfinite(p0):
            # Fallback: szeroki zakres
            p0, pe, E = find_psi0(K, ak, psi_lo, psi_hi, n_pts=20)
        if np.isfinite(p0):
            psi_arr[i] = p0
            gU_arr[i]  = gU_from_E(K, E)
            psi_c = p0

    # Znajdź zmianę znaku g_U
    for i in range(len(gU_arr)-1):
        gi, gi1 = gU_arr[i], gU_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)):
            continue
        if gi*gi1 < 0:
            Ki, Ki1 = K_arr[i], K_arr[i+1]
            psi_mid = 0.5*(psi_arr[i] + (psi_arr[i+1] if np.isfinite(psi_arr[i+1]) else psi_arr[i]))

            psi_c2 = [psi_mid]

            def ratio_fn(K_t):
                lo_ = max(psi_c2[0] - 0.20, psi_lo)
                hi_ = min(psi_c2[0] + 0.20, psi_hi)
                p0, pe, E = find_psi0(K_t, ak, lo_, hi_, n_pts=14)
                if not np.isfinite(p0):
                    p0, pe, E = find_psi0(K_t, ak, psi_lo, psi_hi, n_pts=18)
                if np.isfinite(p0):
                    psi_c2[0] = p0
                    return gU_from_E(K_t, E)
                return gi  # fallback

            try:
                K_star = brentq(ratio_fn, Ki, Ki1, xtol=1e-7, maxiter=40)
                lo_ = max(psi_c2[0] - 0.15, psi_lo)
                hi_ = min(psi_c2[0] + 0.15, psi_hi)
                p_fin, pe_fin, E_fin = find_psi0(K_star, ak, lo_, hi_, n_pts=16)
                gU_fin = gU_from_E(K_star, E_fin)
                return K_star, p_fin, gU_fin
            except Exception as e:
                t = -gi/(gi1-gi)
                K_lin = Ki + t*(Ki1-Ki)
                return K_lin, psi_arr[i], gi + t*(gi1-gi)

    return np.nan, np.nan, np.nan


# ── Skan K*₂/K*₁ vs α_K ─────────────────────────────────────────────────────

def scan_one(ak, verbose=True):
    """
    Dla danego α_K: znajdź K*₁ i K*₂, oblicz stosunek.
    Branch A: K ∈ [0.003, 0.025], ψ₀ ∈ [0.9, 2.0]
    Branch B: K ∈ [0.005, 0.35], ψ₀ ∈ [1.8, 4.5]
             (szeroki zakres bo K*₂ mocno zmienia się z α_K!)
    """
    if verbose:
        print(f"\n  α_K = {ak:.4f}:", end="", flush=True)

    # Branch A
    K1, psi1, gU1 = find_K_star_branch(ak,
        K_lo=0.003, K_hi=0.025, psi_lo=0.9, psi_hi=2.0,
        n_K=25, label="A")

    if verbose:
        if np.isfinite(K1):
            print(f"  K*₁={K1:.6f}(g={gU1:.2e})", end="", flush=True)
        else:
            print(f"  K*₁=??? ", end="", flush=True)

    # Branch B: szeroki zakres K (od małego do dużego)
    K2, psi2, gU2 = find_K_star_branch(ak,
        K_lo=0.005, K_hi=0.35, psi_lo=1.8, psi_hi=4.5,
        n_K=30, label="B")

    if verbose:
        if np.isfinite(K2):
            print(f"  K*₂={K2:.6f}(g={gU2:.2e})", end="", flush=True)
        else:
            print(f"  K*₂=??? ", end="", flush=True)

    ratio = K2/K1 if (np.isfinite(K1) and np.isfinite(K2) and K1>0) else np.nan
    if verbose:
        if np.isfinite(ratio):
            print(f"  ratio={ratio:.4f}", flush=True)
        else:
            print(f"  ratio=???", flush=True)

    return K1, psi1, K2, psi2, ratio


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    print("="*72)
    print("p108_v2: argmax K*₂/K*₁ — niezależny skan dla każdego α_K")
    print(f"R_MAX={R_MAX}, N_eval={N_EVAL}, RTOL={RTOL:.0e}")
    print(f"Referencja: α_K={AK_REF}, K*₁={K1_REF:.6f}, K*₂={K2_REF:.6f}, ratio={K2_REF/K1_REF:.4f}")
    print("="*72)

    # Lista α_K: poniżej, przy, powyżej α_K_ref
    # (powyżej α_K_ref wiemy że ratio spada bo K*₂ zbliża się do K*₁)
    AK_VALUES = sorted([
        # Blisko α_SN (wiemy: ratio → 1)
        8.9500, 8.9000, 8.8000, 8.7000,
        # Powyżej referencji
        8.6500, 8.6000,
        # Referencja
        AK_REF,
        # Poniżej referencji — kluczowy zakres
        8.5000, 8.4500, 8.4000, 8.3500, 8.3000, 8.2000, 8.1000,
        8.0000, 7.8000, 7.5000, 7.0000, 6.5000, 6.0000,
    ])

    print(f"\nSkan {len(AK_VALUES)} wartości α_K od {min(AK_VALUES):.2f} do {max(AK_VALUES):.2f}")
    print("Skan niezależny: dla każdego α_K pełny skan K i ψ₀.")
    print("Proszę czekać...\n")

    results = []
    t0 = time.time()
    for ak in AK_VALUES:
        t1 = time.time()
        K1, psi1, K2, psi2, ratio = scan_one(ak, verbose=True)
        dt = time.time() - t1
        print(f"  [czas: {dt:.0f}s]")
        results.append((ak, K1, psi1, K2, psi2, ratio))

    total = time.time() - t0
    print(f"\nŁączny czas: {total:.0f}s = {total/60:.1f} min")

    # ── Tabela zbiorowa ───────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("TABELA: K*₁, K*₂, ratio vs α_K")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}  {'ψ₀^A':>8}  {'ψ₀^B':>8}")
    print("  " + "-"*70)
    valid = []
    for ak, K1, psi1, K2, psi2, rat in results:
        ref = " <-- REF" if abs(ak-AK_REF)<1e-3 else ""
        if np.isfinite(rat):
            print(f"  {ak:>8.4f}  {K1:>12.6f}  {K2:>12.6f}  {rat:>10.4f}  {psi1:>8.5f}  {psi2:>8.5f}{ref}")
            valid.append((ak, K1, K2, rat))
        else:
            k1s = f"{K1:.6f}" if np.isfinite(K1) else "---"
            k2s = f"{K2:.6f}" if np.isfinite(K2) else "---"
            print(f"  {ak:>8.4f}  {k1s:>12}  {k2s:>12}  {'---':>10}{ref}")

    # ── Analiza maksimum ──────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("ANALIZA MAKSIMUM:")
    if len(valid) < 3:
        print("  Za mało punktów do analizy.")
    else:
        aks   = np.array([r[0] for r in valid])
        rats  = np.array([r[3] for r in valid])
        i_max = np.argmax(rats)
        ak_max, _, _, r_max = valid[i_max]

        print(f"  argmax α_K = {ak_max:.4f}  (ratio = {r_max:.4f})")
        print(f"  α_K_ref    = {AK_REF:.4f}  (ratio = {rats[np.argmin(np.abs(aks-AK_REF))]::.4f})")
        print(f"  Δ(argmax − ref) = {ak_max - AK_REF:.4f}")

        # Monotonicity analysis
        print("\n  Gradient ratio (Δratio/Δα_K):")
        for i in range(len(valid)-1):
            dak = aks[i+1] - aks[i]
            dr  = rats[i+1] - rats[i]
            if abs(dak) > 0.001:
                print(f"    [{aks[i]:.3f}→{aks[i+1]:.3f}]:  Δratio={dr:+.3f}  grad={dr/dak:+.3f}")

        print("\n  Top 5 punktów wg ratio:")
        idx_s = np.argsort(rats)[::-1]
        for i in idx_s[:5]:
            ref = " <-- REF" if abs(valid[i][0]-AK_REF)<1e-3 else ""
            print(f"    α_K={valid[i][0]:.4f}  ratio={valid[i][3]:.4f}  K*₁={valid[i][1]:.6f}  K*₂={valid[i][2]:.6f}{ref}")

        print("\n  WNIOSEK:")
        delta = ak_max - AK_REF
        if abs(delta) < 0.05:
            print(f"  *** HIPOTEZA POTWIERDZONA: argmax α_K={ak_max:.4f} ≈ α_K_ref={AK_REF} ***")
        elif abs(delta) < 0.20:
            print(f"  *** CZĘŚCIOWE: argmax α_K={ak_max:.4f}, odchylenie {delta:.4f} ({abs(delta)/AK_REF*100:.2f}%) ***")
        else:
            print(f"  *** HIPOTEZA OBALONA: argmax α_K={ak_max:.4f}, Δ={delta:.4f} ***")
