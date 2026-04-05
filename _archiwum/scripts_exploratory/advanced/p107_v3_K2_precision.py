#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p107_v3_K2_precision.py  --  TGP v1
=====================================
Krok 1: Precyzyjne wyznaczenie K*₂ dla Branch B.

Problem: p107_v2 daje g_U^B = 1.652e-02 (1.6%) z powodu katastrofalnego
kasowania (E_kin≈233, E_pot≈−232, E_tot≈1.28). Małe zmiany K lub ψ₀
bardzo silnie wpływają na E_tot.

Rozwiązanie: Skanuj K w zakresie [0.098, 0.104] przy ultra-wysokiej
precyzji (N_EVAL=3000, RTOL=1e-12), szukaj K*₂ gdzie g_U(K)=0 (brentq).

Krok 2: Oblicz M_conf(p) z poprawionymi K*₂, ψ₀.

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
ALPHA_K  = 8.5616
GAM      = 1.0
V1       = GAM/3.0 - GAM/4.0
LAM      = 5.501357e-06
AGAM     = 0.040
R_MAX    = 80.0
N_EVAL   = 3000       # ultra-wysokie zagęszczenie
RTOL     = 1e-12
ATOL     = 1e-14

TARGET   = 206.77   # m_μ/m_e

# Referencja Branch A (stabilna, bez problemów z precyzją)
K1_REF  = 0.010280;  PSI1_REF = 1.2417;  PSI1_LO = 1.22;  PSI1_HI = 1.27

# Referencja Branch B (do inicjalizacji szukania ψ₀)
K2_REF  = 0.100276;  PSI2_REF = 2.7676;  PSI2_LO = 2.65;  PSI2_HI = 2.90

# ── ODE ───────────────────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6

def dV_mod(phi):
    return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def rhs(r, y):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K/phi
    return [dphi,
            dV_mod(phi)/kfac + ALPHA_K*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]


def solve_sol(K, psi0, n_eval=None, rtol=None, atol=None):
    """Integruje ODE, zwraca sol lub None."""
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


def phi_end(psi, K):
    sol = solve_sol(K, psi)
    return sol.y[0,-1] - 1.0 if sol else np.nan


def refine_psi0(K, lo, hi, psi_ref, n_pts=60, label=""):
    """Brentq ψ₀ na φ(r_max)=1. Skanuje n_pts punktów jeśli brak znaku."""
    try:
        flo = phi_end(lo, K)
        fhi = phi_end(hi, K)
        if np.isfinite(flo) and np.isfinite(fhi) and flo*fhi < 0:
            return brentq(lambda p: phi_end(p, K), lo, hi, xtol=1e-10, maxiter=100)
        psis = np.linspace(lo, hi, n_pts)
        Fv   = [phi_end(p, K) for p in psis]
        for i in range(len(Fv)-1):
            if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
                return brentq(lambda p: phi_end(p, K), psis[i], psis[i+1],
                              xtol=1e-10, maxiter=100)
        print(f"    {label}: UWAGA — brak zmiany znaku w [{lo:.4f},{hi:.4f}], fallback ref")
        return psi_ref
    except Exception as e:
        print(f"    {label}: UWAGA — brentq błąd ({e}), fallback ref")
        return psi_ref


def compute_E_gU(K, psi0):
    """Oblicza E i g_U dla danego (K, ψ₀)."""
    sol = solve_sol(K, psi0)
    if sol is None:
        return np.nan, np.nan, np.nan
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.+ALPHA_K/phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vv   = np.array([V_mod(float(q)) for q in phi])
    Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
    E    = Ek + Ep
    gU   = E/(4*np.pi*K) - 1.0
    return E, gU, sol.y[0,-1]


def M_conf(sol, p):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.+ALPHA_K/phi
    Vv   = np.array([V_mod(float(q)) for q in phi])
    eps  = 0.5*dphi**2*kfac + (Vv-V1)
    return 4*np.pi*np.trapezoid(eps*phi**(3.*p)*r**2, r)


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("="*72)
    print("p107_v3: Precyzyjne K*₂ + M_conf(p)")
    print(f"alpha_K={ALPHA_K}, a_Gamma={AGAM}, R_max={R_MAX}")
    print(f"N_eval={N_EVAL}, RTOL={RTOL:.0e}, ATOL={ATOL:.0e}")
    print("="*72)

    # ── ETAP 1: Skan g_U(K) dla Branch B ─────────────────────────────────────
    print("\n[1] Skan g_U(K) dla Branch B (gruba siatka K):")
    print(f"  {'K':>10}  {'ψ₀':>12}  {'E':>12}  {'4πK':>12}  {'g_U':>12}  {'φ_end':>10}")

    # Gruba siatka, żeby zobaczyć kształt funkcji g_U(K)
    K_scan = [0.098, 0.099, 0.0995, 0.1000, 0.1003, 0.1005, 0.1010, 0.1020, 0.103, 0.104]

    gU_vals = []
    psi_prev = PSI2_REF
    lo_prev  = PSI2_LO
    hi_prev  = PSI2_HI

    for K in K_scan:
        # Adaptuj zakres szukania ψ₀ w zależności od K (K↑ → ψ₀↓ empirycznie)
        psi0 = refine_psi0(K, lo_prev, hi_prev, psi_prev, n_pts=50, label=f"K={K:.4f}")
        E, gU, phi_e = compute_E_gU(K, psi0)
        print(f"  {K:>10.5f}  {psi0:>12.8f}  {E:>12.6f}  {4*np.pi*K:>12.6f}  {gU:>12.4e}  {phi_e:>10.6f}")
        gU_vals.append((K, psi0, E, gU))
        if np.isfinite(psi0) and abs(psi0 - psi_prev) < 0.3:
            psi_prev = psi0

    # ── ETAP 2: Brentq na g_U(K) ─────────────────────────────────────────────
    print("\n[2] Szukanie K*₂ gdzie g_U(K)=0 (brentq):")

    # Znajdź przedział ze zmianą znaku g_U
    K_bracket = None
    for i in range(len(gU_vals)-1):
        Ki, _, _, gUi   = gU_vals[i]
        Ki1, _, _, gUi1 = gU_vals[i+1]
        if np.isfinite(gUi) and np.isfinite(gUi1) and gUi*gUi1 < 0:
            K_bracket = (Ki, Ki1)
            print(f"  Zmiana znaku g_U: K ∈ ({Ki:.5f}, {Ki1:.5f})")
            break

    if K_bracket is None:
        print("  BRAK zmiany znaku! Używam K*₂=K2_REF.")
        K2_star = K2_REF
        psi2_star = PSI2_REF
    else:
        # Funkcja g_U(K) — dla każdego K refinuje ψ₀
        psi_cache = {K: psi for K, psi, _, _ in gU_vals}

        def gU_fn(K_test):
            # Interpoluj ψ₀ startowy z cachowanej siatki
            Ks = sorted(psi_cache.keys())
            psi_start = psi_cache.get(K_test)
            if psi_start is None:
                # Interpoluj z sąsiadów
                import bisect
                idx = bisect.bisect_left(Ks, K_test)
                if idx == 0:
                    psi_start = psi_cache[Ks[0]]
                elif idx >= len(Ks):
                    psi_start = psi_cache[Ks[-1]]
                else:
                    K_lo, K_hi = Ks[idx-1], Ks[idx]
                    t = (K_test-K_lo)/(K_hi-K_lo)
                    psi_start = (1-t)*psi_cache[K_lo] + t*psi_cache[K_hi]
            # Refinement ψ₀ z wąskim zakresem wokół psi_start
            dp = 0.1
            psi0 = refine_psi0(K_test, max(psi_start-dp, 2.4), min(psi_start+dp, 3.0),
                               psi_start, n_pts=40, label=f"  brentq K={K_test:.6f}")
            E, gU, phi_e = compute_E_gU(K_test, psi0)
            psi_cache[K_test] = psi0
            return gU if np.isfinite(gU) else 0.0

        try:
            K2_star = brentq(gU_fn, K_bracket[0], K_bracket[1], xtol=1e-6, maxiter=40)
            print(f"  K*₂ = {K2_star:.8f}")
        except Exception as e:
            print(f"  Brentq K*₂ błąd: {e} — fallback K2_REF")
            K2_star = K2_REF

        # Precyzyjna ψ₀ dla K*₂
        psi2_star = refine_psi0(K2_star, PSI2_LO, PSI2_HI, PSI2_REF, n_pts=60, label="K*2_final")
        E2, gU2, phi_e2 = compute_E_gU(K2_star, psi2_star)
        print(f"  Weryfikacja: K*₂={K2_star:.8f}  ψ₀={psi2_star:.8f}  "
              f"E={E2:.6f}  4πK={4*np.pi*K2_star:.6f}  g_U={gU2:.4e}  φ_end={phi_e2:.6f}")

    # ── ETAP 3: Branch A ──────────────────────────────────────────────────────
    print("\n[3] Branch A (precyzyjna):")
    psi1_star = refine_psi0(K1_REF, PSI1_LO, PSI1_HI, PSI1_REF, n_pts=60, label="Branch A")
    E1, gU1, phi_e1 = compute_E_gU(K1_REF, psi1_star)
    print(f"  K*₁={K1_REF:.6f}  ψ₀={psi1_star:.8f}  "
          f"E={E1:.6f}  4πK={4*np.pi*K1_REF:.6f}  g_U={gU1:.4e}  φ_end={phi_e1:.6f}")

    # ── ETAP 4: ODE solucje ───────────────────────────────────────────────────
    print("\n[4] ODE (solucje finalne):")
    solA = solve_sol(K1_REF, psi1_star)
    solB = solve_sol(K2_star, psi2_star)
    if solA is None or solB is None:
        print("BŁĄD: ODE failed!"); sys.exit(1)

    r_A = solA.t; r_B = solB.t
    phi_A = solA.y[0]; phi_B = solB.y[0]
    print(f"  Branch A: φ₀={phi_A[0]:.6f}  φ_end={phi_A[-1]:.8f}  φ_max={phi_A.max():.4f}")
    print(f"  Branch B: φ₀={phi_B[0]:.6f}  φ_end={phi_B[-1]:.8f}  φ_max={phi_B.max():.4f}")

    # ── ETAP 5: Tabela M_conf ─────────────────────────────────────────────────
    P_VALS = [0, 0.25, 0.50, 0.75, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.25, 1.50, 2.00, 3.00]

    print("\n[5] Tabela M_conf^B(p) / M_conf^A(p):")
    print(f"  {'p':>6}  {'M_A':>14}  {'M_B':>14}  {'ratio':>10}  {'vs 206.77':>10}")
    print("  " + "-"*64)
    ratios_list = []
    for p in P_VALS:
        mA = M_conf(solA, p)
        mB = M_conf(solB, p)
        r  = mB/mA if mA != 0 else np.nan
        if np.isfinite(r):
            delta = abs(r-TARGET)/TARGET*100
            mk = "  <-- BLISKIE!" if delta < 3 else ("  <--" if delta < 10 else "")
            print(f"  {p:>6.2f}  {mA:>14.8f}  {mB:>14.6f}  {r:>10.4f}  {delta:>9.2f}%{mk}")
        else:
            print(f"  {p:>6.2f}  {mA:>14.8f}  {mB:>14.6f}  {'NaN':>10}")
        ratios_list.append((p, mA, mB, r))

    # ── ETAP 6: Szukanie p* ───────────────────────────────────────────────────
    print("\n[6] Szukanie p* gdzie ratio = 206.77:")
    p_arr = np.array([x[0] for x in ratios_list])
    r_arr = np.array([x[3] for x in ratios_list])
    p_star = None

    for i in range(len(r_arr)-1):
        ri, ri1 = r_arr[i], r_arr[i+1]
        if not (np.isfinite(ri) and np.isfinite(ri1)):
            continue
        if (ri-TARGET)*(ri1-TARGET) < 0:
            t = (TARGET-ri)/(ri1-ri)
            p_lin = p_arr[i] + t*(p_arr[i+1]-p_arr[i])
            print(f"  Bracket: p ∈ ({p_arr[i]:.3f}, {p_arr[i+1]:.3f}), lin.interp.: {p_lin:.4f}")

            def ratio_fn(p_t):
                mA_ = M_conf(solA, p_t)
                mB_ = M_conf(solB, p_t)
                return mB_/mA_ - TARGET if mA_ != 0 else np.nan

            try:
                p_star = brentq(ratio_fn, p_arr[i], p_arr[i+1], xtol=1e-6, maxiter=60)
                print(f"  Brentq: p* = {p_star:.7f}")
            except Exception as e:
                print(f"  Brentq błąd: {e}  →  używam lin. interp.")
                p_star = p_lin
            break

    if p_star is not None:
        nice = [(n/d, f"{n}/{d}") for n in range(1,13) for d in range(1,9)
                if 0 < n/d <= 3.0]
        nice += [(float(n), str(n)) for n in range(1,4)]
        nice = list({round(v,7):(v,s) for v,s in sorted(nice)}.values())
        nice.sort(key=lambda x: abs(x[0]-p_star))
        print(f"\n  *** p* = {p_star:.7f} ***")
        print("  Najbliższe 'ładne' frakcje:")
        for v, s in nice[:8]:
            check_r = M_conf(solB, v)/M_conf(solA, v)
            print(f"    p = {v:.6f} ({s:>6}):  odch = {abs(v-p_star):.6f}  "
                  f"ratio = {check_r:.4f}  (vs 206.77: {abs(check_r-TARGET)/TARGET*100:.2f}%)")
    else:
        rmin = np.nanmin(r_arr); rmax = np.nanmax(r_arr)
        print(f"  Ratio w zakresie [{rmin:.2f}, {rmax:.2f}] — cel 206.77 nie osiągnięty")

    # ── Podsumowanie ─────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("PODSUMOWANIE:")
    print(f"  K*₁ = {K1_REF:.6f}  ψ₀^A = {psi1_star:.8f}  g_U^A = {gU1:.4e}")
    print(f"  K*₂ = {K2_star:.8f}  ψ₀^B = {psi2_star:.8f}  g_U^B = {gU2 if K_bracket else float('nan'):.4e}")
    print(f"  K*₂/K*₁ = {K2_star/K1_REF:.6f}")

    idx0 = [i for i,x in enumerate(ratios_list) if abs(x[0]-0.0)<1e-6]
    idx1 = [i for i,x in enumerate(ratios_list) if abs(x[0]-1.0)<1e-6]
    if idx0: print(f"  M_conf ratio (p=0): {ratios_list[idx0[0]][3]:.6f}")
    if idx1: print(f"  M_conf ratio (p=1): {ratios_list[idx1[0]][3]:.4f}  (vs 206.77: {abs(ratios_list[idx1[0]][3]-TARGET)/TARGET*100:.2f}%)")
    if p_star is not None:
        print(f"  p* = {p_star:.7f}")
        print(f"  Odchylenie od p=1: {abs(p_star-1.0):.4f}  ({abs(p_star-1.0)/1.0*100:.2f}%)")
