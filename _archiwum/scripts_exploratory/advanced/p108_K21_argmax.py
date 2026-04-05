#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_K21_argmax.py  --  TGP v1
=================================
Ścieżka 5 OP-3: argmax K*₂/K*₁ jako zasada selekcji α_K

Hipoteza: α_K = argmax_{α_K}[K*₂(α_K)/K*₁(α_K)]
Fizyczna interpretacja: natura wybiera α_K przy którym dwa
solitonowe stany są NAJBARDZIEJ rozróżnialne.

Metoda:
- Skan α_K ∈ [5.0, α_SN≈8.996] krokiem adaptacyjnym
- Dla każdego α_K: brentq na g_U^A(K)=0 → K*₁
                   brentq na g_U^B(K)=0 → K*₂
- Continuation od referencji (α_K=8.5616, K*₁=0.010280, K*₂=0.09999)
- Szukanie maksimum ratio K*₂/K*₁

Parametry referencyjne: wyniki p107_v3 (K*₂=0.09999348 przy α_K=8.5616)
α_SN = 8.9964 ± 0.0004

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
N_EVAL   = 2000
RTOL     = 1e-11
ATOL     = 1e-13

ALPHA_SN = 8.9964    # górna granica istnienia Branch B

# Referencja przy α_K = 8.5616 (z p107_v3)
AK_REF   = 8.5616
K1_REF   = 0.010280;  PSI1_REF = 1.2417
K2_REF   = 0.09999348; PSI2_REF = 2.7624

# ── ODE ───────────────────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6

def dV_mod(phi):
    return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def make_rhs(alpha_k):
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
        kfac = 1.0 + alpha_k/phi
        return [dphi,
                dV_mod(phi)/kfac + alpha_k*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
    return rhs


def solve_sol(K, psi0, alpha_k, n_eval=None, rtol=None, atol=None):
    rhs = make_rhs(alpha_k)
    if n_eval is None: n_eval = N_EVAL
    if rtol   is None: rtol   = RTOL
    if atol   is None: atol   = ATOL
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, n_eval)
    for method in ('DOP853', 'Radau'):
        try:
            sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                            method=method, rtol=rtol, atol=atol, t_eval=r_ev)
            if sol.t[-1] >= R_MAX*0.99:
                return sol
        except Exception:
            pass
    return None


def phi_end(psi, K, alpha_k):
    sol = solve_sol(K, psi, alpha_k)
    return sol.y[0,-1] - 1.0 if sol else np.nan


def refine_psi0(K, alpha_k, psi_ref, dpsi=0.15, n_pts=30, label=""):
    """Brentq ψ₀ na φ(r_max)=1. Zakres: [psi_ref-dpsi, psi_ref+dpsi]."""
    lo = max(psi_ref - dpsi, 0.5)
    hi = psi_ref + dpsi
    try:
        flo = phi_end(lo, K, alpha_k)
        fhi = phi_end(hi, K, alpha_k)
        if np.isfinite(flo) and np.isfinite(fhi) and flo*fhi < 0:
            return brentq(lambda p: phi_end(p, K, alpha_k), lo, hi, xtol=1e-9, maxiter=80)
        psis = np.linspace(lo, hi, n_pts)
        Fv   = [phi_end(p, K, alpha_k) for p in psis]
        for i in range(len(Fv)-1):
            if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
                return brentq(lambda p: phi_end(p, K, alpha_k), psis[i], psis[i+1],
                              xtol=1e-9, maxiter=80)
    except Exception:
        pass
    return psi_ref  # fallback


def gU_for_K(K, alpha_k, psi_ref, dpsi=0.15):
    """Oblicza g_U(K) = E/(4πK)−1 dla danego α_K, K (z refinementem ψ₀)."""
    psi0 = refine_psi0(K, alpha_k, psi_ref, dpsi=dpsi)
    sol  = solve_sol(K, psi0, alpha_k)
    if sol is None:
        return np.nan, psi0
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.+alpha_k/phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vv   = np.array([V_mod(float(q)) for q in phi])
    Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
    E    = Ek + Ep
    return E/(4*np.pi*K) - 1.0, psi0


def find_K_star(alpha_k, K_lo, K_hi, psi_ref, dpsi=0.15, branch_label=""):
    """
    Szuka K* gdzie g_U(K)=0 w [K_lo, K_hi].
    Zwraca (K_star, psi0_star) lub (nan, nan).
    """
    # Skan grubej siatki
    n_scan = 12
    Kv = np.linspace(K_lo, K_hi, n_scan)
    gv = []
    psi_v = []
    psi_c = psi_ref
    for K in Kv:
        g, psi_c = gU_for_K(K, alpha_k, psi_c, dpsi=dpsi)
        gv.append(g)
        psi_v.append(psi_c)

    # Znajdź zmianę znaku
    for i in range(len(gv)-1):
        gi, gi1 = gv[i], gv[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)):
            continue
        if gi*gi1 < 0:
            Ki, Ki1 = Kv[i], Kv[i+1]
            psi_mid = 0.5*(psi_v[i]+psi_v[i+1])
            psi_cache = [psi_mid]
            try:
                def fn(K_t):
                    g_t, psi_t = gU_for_K(K_t, alpha_k, psi_cache[0], dpsi=dpsi*0.7)
                    psi_cache[0] = psi_t
                    return g_t if np.isfinite(g_t) else 0.0
                K_star = brentq(fn, Ki, Ki1, xtol=1e-7, maxiter=50)
                g_fin, psi_fin = gU_for_K(K_star, alpha_k, psi_cache[0], dpsi=dpsi*0.5)
                return K_star, psi_fin
            except Exception as e:
                # Interpolacja liniowa fallback
                t = -gv[i]/(gv[i+1]-gv[i])
                K_star = Ki + t*(Ki1-Ki)
                return K_star, psi_mid

    return np.nan, psi_ref


# ── Skan α_K ─────────────────────────────────────────────────────────────────

def scan_alpha_k(AK_list):
    """
    Skanuje listę α_K, dla każdego oblicza K*₁, K*₂ i ratio.
    Używa continuation (wartości z poprzedniego α_K jako start).
    """
    results = []

    # Skan malejący: od AK_REF w dół
    print("\n[A] Skan malejący α_K (od referencji w dół):")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}  {'ψ₀^A':>10}  {'ψ₀^B':>10}")
    print("  " + "-"*72)

    # Startujemy od referencji
    k1_c = K1_REF; psi1_c = PSI1_REF
    k2_c = K2_REF; psi2_c = PSI2_REF

    for ak in AK_list:
        if ak > AK_REF + 1e-6:
            continue  # pominięcie — to skan malejący

        # Branch A: szukaj w [K*₁/3, K*₁*3]
        K1_lo = max(k1_c * 0.3, 1e-4)
        K1_hi = k1_c * 4.0
        k1_new, psi1_new = find_K_star(ak, K1_lo, K1_hi, psi1_c, dpsi=0.20, branch_label="A")

        # Branch B: szukaj w [k2_c*0.4, k2_c*2.5] (ale nie poniżej K*₁)
        if np.isfinite(k2_c):
            K2_lo = max(k2_c * 0.4, (k1_new if np.isfinite(k1_new) else k1_c)*2)
            K2_hi = min(k2_c * 2.5, 0.5)
            k2_new, psi2_new = find_K_star(ak, K2_lo, K2_hi, psi2_c, dpsi=0.20, branch_label="B")
        else:
            k2_new, psi2_new = np.nan, psi2_c

        ratio = k2_new/k1_new if (np.isfinite(k1_new) and np.isfinite(k2_new) and k1_new>0) else np.nan

        flag = ""
        if np.isfinite(ratio):
            flag = " ***" if ratio > 9.5 else ""

        print(f"  {ak:>8.4f}  {k1_new:>12.6f}  {k2_new:>12.6f}  {ratio:>10.4f}  {psi1_new:>10.6f}  {psi2_new:>10.6f}{flag}")
        results.append((ak, k1_new, k2_new, ratio))

        # Aktualizuj continuation tylko jeśli wynik sensowny
        if np.isfinite(k1_new) and k1_new > 0:
            k1_c = k1_new; psi1_c = psi1_new
        if np.isfinite(k2_new) and k2_new > 0:
            k2_c = k2_new; psi2_c = psi2_new

    return results


def scan_alpha_k_ascending(AK_list, k1_start, psi1_start, k2_start, psi2_start):
    """Skan rosnący α_K (od referencji w górę, do α_SN)."""
    results = []
    print("\n[B] Skan rosnący α_K (od referencji w górę, do α_SN):")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}  {'ψ₀^A':>10}  {'ψ₀^B':>10}")
    print("  " + "-"*72)

    k1_c = k1_start; psi1_c = psi1_start
    k2_c = k2_start; psi2_c = psi2_start

    for ak in AK_list:
        if ak <= AK_REF + 1e-6:
            continue

        K1_lo = max(k1_c * 0.4, 1e-4)
        K1_hi = k1_c * 4.0
        k1_new, psi1_new = find_K_star(ak, K1_lo, K1_hi, psi1_c, dpsi=0.20, branch_label="A")

        if np.isfinite(k2_c):
            K2_lo = max(k2_c * 0.4, (k1_new if np.isfinite(k1_new) else k1_c)*1.5)
            K2_hi = min(k2_c * 2.0, 0.5)
            k2_new, psi2_new = find_K_star(ak, K2_lo, K2_hi, psi2_c, dpsi=0.20, branch_label="B")
        else:
            k2_new, psi2_new = np.nan, psi2_c

        ratio = k2_new/k1_new if (np.isfinite(k1_new) and np.isfinite(k2_new) and k1_new>0) else np.nan

        alive = "żyje" if (np.isfinite(k2_new) and np.isfinite(ratio)) else "DEAD"
        print(f"  {ak:>8.4f}  {k1_new:>12.6f}  {k2_new:>12.6f}  {ratio:>10.4f}  {psi1_new:>10.6f}  {psi2_new:>10.6f}  [{alive}]")
        results.append((ak, k1_new, k2_new, ratio))

        if np.isfinite(k1_new) and k1_new > 0:
            k1_c = k1_new; psi1_c = psi1_new
        if np.isfinite(k2_new) and k2_new > 0 and alive == "żyje":
            k2_c = k2_new; psi2_c = psi2_new
        else:
            k2_c = np.nan  # Branch B martwa

    return results


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("="*72)
    print("p108: argmax_{α_K}[K*₂/K*₁] jako zasada selekcji α_K")
    print(f"Referencja: α_K={AK_REF}, K*₁={K1_REF:.6f}, K*₂={K2_REF:.6f}")
    print(f"R_MAX={R_MAX}, N_eval={N_EVAL}, RTOL={RTOL:.0e}")
    print(f"α_SN ≈ {ALPHA_SN}")
    print("="*72)

    # ─── Lista α_K do przeskanowania ─────────────────────────────────────────
    # Poniżej referencji: gruba siatka (szybko maleją K*)
    AK_below = list(np.arange(4.0, 8.5616+0.01, 0.25))
    # Wokół referencji: gęsta siatka
    AK_near  = list(np.arange(8.10, 8.5616+0.01, 0.05))
    # Powyżej: do α_SN
    AK_above = list(np.arange(8.60, 8.99, 0.05))

    # Łączona lista (posortowana malejąco dla skan_A, rosnąco dla skan_B)
    AK_all_desc = sorted(set(AK_below + AK_near + [AK_REF]), reverse=True)
    AK_all_asc  = sorted(set(AK_above + [AK_REF]))

    print(f"\nPunkty skanowania (malejąco): {len(AK_all_desc)} punktów od α_K={max(AK_all_desc):.2f} do {min(AK_all_desc):.2f}")
    print(f"Punkty skanowania (rosnąco):  {len(AK_all_asc)} punktów od α_K={min(AK_all_asc):.2f} do {max(AK_all_asc):.2f}")

    # ── Skan malejący (od referencji w dół) ──────────────────────────────────
    results_desc = scan_alpha_k(AK_all_desc)

    # ── Skan rosnący (od referencji w górę) ──────────────────────────────────
    results_asc  = scan_alpha_k_ascending(AK_all_asc, K1_REF, PSI1_REF, K2_REF, PSI2_REF)

    # ── Połącz wyniki ─────────────────────────────────────────────────────────
    all_results = sorted(results_desc + results_asc, key=lambda x: x[0])
    # Usuń duplikat referencji
    seen = set()
    all_results_uniq = []
    for row in all_results:
        ak_r = round(row[0], 6)
        if ak_r not in seen:
            seen.add(ak_r)
            all_results_uniq.append(row)

    # ── Pełna tabela ─────────────────────────────────────────────────────────
    print("\n[C] Pełna tabela K*₂/K*₁ vs α_K:")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}")
    print("  " + "-"*50)
    valid = []
    for ak, k1, k2, rat in all_results_uniq:
        flag = " <-- REF" if abs(ak-AK_REF)<1e-3 else ""
        if np.isfinite(rat):
            print(f"  {ak:>8.4f}  {k1:>12.6f}  {k2:>12.6f}  {rat:>10.4f}{flag}")
            valid.append((ak, k1, k2, rat))
        else:
            print(f"  {ak:>8.4f}  {k1:>12.6f}  {'---':>12}  {'---':>10}  [Branch B dead]")

    # ── Analiza maksimum ──────────────────────────────────────────────────────
    print("\n[D] Analiza maksimum K*₂/K*₁:")
    if len(valid) < 3:
        print("  Za mało punktów do analizy!")
    else:
        ratios = np.array([r for _, _, _, r in valid])
        aks    = np.array([a for a, _, _, _ in valid])

        i_max = np.argmax(ratios)
        print(f"  Maksimum: ratio={ratios[i_max]:.4f} przy α_K={aks[i_max]:.4f}")
        print(f"  Referencja: α_K={AK_REF} → ratio={ratios[np.argmin(np.abs(aks-AK_REF))]:.4f}")

        # Sprawdź czy maksimum jest przy α_K = AK_REF
        delta = aks[i_max] - AK_REF
        print(f"  Odchylenie: α_K_max − α_K_ref = {delta:.4f}")
        if abs(delta) < 0.05:
            print("  *** HIPOTEZA POTWIERDZONA: maksimum blisko α_K_ref! ***")
        elif abs(delta) < 0.15:
            print("  *** CZĘŚCIOWE potwierdzenie: maksimum w pobliżu α_K_ref ***")
        else:
            print("  *** HIPOTEZA OBALONA: maksimum daleko od α_K_ref ***")

        # Pokaż top 5
        idx_sorted = np.argsort(ratios)[::-1]
        print("\n  Top 5 α_K wg ratio K*₂/K*₁:")
        for i in idx_sorted[:5]:
            mark = " <-- REF" if abs(aks[i]-AK_REF)<1e-3 else ""
            print(f"    α_K={aks[i]:.4f}  ratio={ratios[i]:.4f}  K*₁={valid[i][1]:.6f}  K*₂={valid[i][2]:.6f}{mark}")

        # Pochodna (czy funkcja jest monotoniczna?)
        print("\n  Gradient ratio (delta_ratio/delta_α_K):")
        for i in range(1, len(valid)):
            dak = aks[i] - aks[i-1]
            dr  = ratios[i] - ratios[i-1]
            if abs(dak) > 0.001:
                grad = dr/dak
                print(f"    [{aks[i-1]:.3f}→{aks[i]:.3f}]: Δratio={dr:+.3f}  grad={grad:+.3f}")

    # ── Podsumowanie ─────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("PODSUMOWANIE Ścieżka 5 (p108):")
    if valid:
        i_max = np.argmax([r for _, _, _, r in valid])
        ak_max, k1_max, k2_max, r_max = valid[i_max]
        print(f"  argmax α_K = {ak_max:.4f}")
        print(f"  max ratio  = {r_max:.4f}")
        print(f"  α_K_ref    = {AK_REF}")
        print(f"  Δ          = {ak_max - AK_REF:.4f}  ({abs(ak_max-AK_REF)/AK_REF*100:.2f}%)")
        if abs(ak_max - AK_REF) < 0.05:
            print("  WYNIK: α_K=8.5616 jest maksimum K*₂/K*₁ → HIPOTEZA POTWIERDZONA!")
        else:
            print(f"  WYNIK: maksimum przy α_K={ak_max:.4f} ≠ {AK_REF} → hipoteza wymaga rewizji")
