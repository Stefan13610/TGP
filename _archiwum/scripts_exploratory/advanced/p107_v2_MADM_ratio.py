#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p107_v2_MADM_ratio.py  --  TGP v1
====================================
Ścieżka 4 OP-3: Masa konformalna z TGP-metryki

M_conf(p) = 4π ∫ ε(r)·φ(r)^(3p)·r²dr

gdzie ε(r) = ½(φ')²·kfac + (V_mod(φ)−V₁) to gęstość energii,
a (φ/φ₀)^(3p) = φ^(3p) (z φ₀=1) to element objętości metryki
g_ij = φ^(2p)·δ_ij.

Poprawki względem p107:
  - R_MAX=80 (zgodne z referencją p103_v4)
  - ψ₀ refinement w [lo,hi] — osobne zakresy dla A i B
  - Weryfikacja g_U ≈ 0 PRZED obliczeniami
  - Jeśli refinement zawiedzie → użyj wartości referencyjnej

Pytanie: czy istnieje "ładne" p takie że M_conf^B/M_conf^A = 206.77?

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
R_MAX    = 80.0     # ← kluczowe: zgodne z K* z p103_v4
N_EVAL   = 1500
RTOL     = 1e-10
ATOL     = 1e-12

# Wartości referencyjne (p103_v4, r_max=80)
K1_REF  = 0.010280;  PSI1_REF = 1.2417;  PSI1_LO = 1.22;  PSI1_HI = 1.27
K2_REF  = 0.100276;  PSI2_REF = 2.7676;  PSI2_LO = 2.70;  PSI2_HI = 2.83

TARGET  = 206.77   # m_μ/m_e

# Zakres p do skanowania
P_VALS  = [0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75,
           2.00, 2.50, 3.00, 4.00, 5.00]


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


def solve_sol(K, psi0):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    for method in ('DOP853', 'Radau'):
        try:
            sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                            method=method, rtol=RTOL, atol=ATOL, t_eval=r_ev)
            if sol.t[-1] >= R_MAX*0.99:
                return sol
        except Exception:
            pass
    return None


def phi_end_fn(psi, K):
    sol = solve_sol(K, psi)
    return sol.y[0,-1] - 1.0 if sol else np.nan


def refine_psi0(K, lo, hi, psi_ref):
    """Brentq ψ₀ na φ(r_max)=1. Fallback: psi_ref."""
    try:
        flo = phi_end_fn(lo, K)
        fhi = phi_end_fn(hi, K)
        if np.isfinite(flo) and np.isfinite(fhi) and flo*fhi < 0:
            return brentq(lambda p: phi_end_fn(p, K), lo, hi,
                          xtol=1e-9, maxiter=80)
        else:
            # Szukaj znaku zmiany w gęstszej siatce
            psis = np.linspace(lo, hi, 40)
            Fv   = [phi_end_fn(p, K) for p in psis]
            for i in range(len(Fv)-1):
                if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
                    return brentq(lambda p: phi_end_fn(p, K), psis[i], psis[i+1],
                                  xtol=1e-9, maxiter=60)
            print(f"    UWAGA: brak zmiany znaku w [{lo},{hi}] — fallback ref")
            return psi_ref
    except Exception as e:
        print(f"    UWAGA: brentq błąd ({e}) — fallback ref")
        return psi_ref


def check_gU(sol, K, label):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.+ALPHA_K/phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vv   = np.array([V_mod(float(q)) for q in phi])
    Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
    E    = Ek+Ep
    gU   = E/(4*np.pi*K)-1.0
    print(f"  {label}: psi0={sol.y[0,0]:.6f}  phi_end={sol.y[0,-1]:.6f}  "
          f"E={E:.6f}  4πK={4*np.pi*K:.6f}  g_U={gU:.3e}")
    return E, gU


# ── Integrały konforme ─────────────────────────────────────────────────────────

def M_conf(sol, p):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.+ALPHA_K/phi
    Vv   = np.array([V_mod(float(q)) for q in phi])
    eps  = 0.5*dphi**2*kfac + (Vv-V1)
    return 4*np.pi*np.trapezoid(eps*phi**(3.*p)*r**2, r)

def M_conf_kin(sol, p):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.+ALPHA_K/phi
    return 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*phi**(3.*p)*r**2, r)

def M_conf_pot(sol, p):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    Vv   = np.array([V_mod(float(q)) for q in phi])
    return 4*np.pi*np.trapezoid((Vv-V1)*phi**(3.*p)*r**2, r)


# ── Alternatywna masa: tylko z gęstości energii kinetycznej kfac·(φ')²·r² ─────

def M_kin_only(sol):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.+ALPHA_K/phi
    return 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)

def M_pot_only(sol):
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    Vv   = np.array([V_mod(float(q)) for q in phi])
    return 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("="*72)
    print("p107_v2: M_conf(p) = ∫ε·φ^(3p)·4πr²dr — poprawiona wersja")
    print(f"alpha_K={ALPHA_K}, a_Gamma={AGAM}, R_max={R_MAX}, N_eval={N_EVAL}")
    print("="*72)

    # 1. Refinement ψ₀
    print("\n[1] Refinement ψ₀ (brentq, r_max=80):")
    psi1 = refine_psi0(K1_REF, PSI1_LO, PSI1_HI, PSI1_REF)
    psi2 = refine_psi0(K2_REF, PSI2_LO, PSI2_HI, PSI2_REF)
    print(f"  Branch A: ψ₀ = {psi1:.8f}  (ref={PSI1_REF})")
    print(f"  Branch B: ψ₀ = {psi2:.8f}  (ref={PSI2_REF})")

    # 2. ODE
    print("\n[2] ODE:")
    solA = solve_sol(K1_REF, psi1)
    solB = solve_sol(K2_REF, psi2)
    if solA is None or solB is None:
        print("BŁĄD: ODE failed!"); sys.exit(1)
    print(f"  Branch A: phi_end={solA.y[0,-1]:.8f}  phi_max={solA.y[0].max():.4f}")
    print(f"  Branch B: phi_end={solB.y[0,-1]:.8f}  phi_max={solB.y[0].max():.4f}")

    # 3. Weryfikacja g_U
    print("\n[3] g_U = E/(4πK)−1 (MUSI być ≈0):")
    EA, gA = check_gU(solA, K1_REF, "Branch A")
    EB, gB = check_gU(solB, K2_REF, "Branch B")
    print(f"  E_B/E_A = {EB/EA:.6f}  (oczekiwane ≈ K*₂/K*₁ = {K2_REF/K1_REF:.6f})")
    if abs(gA) > 0.05 or abs(gB) > 0.05:
        print("  OSTRZEŻENIE: g_U jest duże — wyniki M_conf będą błędne!")

    # 4. Struktura energetyczna obu gałęzi
    print("\n[4] Struktura energii (p=0):")
    EkA = M_kin_only(solA); EpA = M_pot_only(solA)
    EkB = M_kin_only(solB); EpB = M_pot_only(solB)
    print(f"  Branch A: E_kin={EkA:.6f}  E_pot={EpA:.6f}  E_tot={EkA+EpA:.6f}")
    print(f"  Branch B: E_kin={EkB:.6f}  E_pot={EpB:.6f}  E_tot={EkB+EpB:.6f}")
    print(f"  Ratio kin: {EkB/EkA:.4f}  pot: {EpB/EpA:.4f}  tot: {(EkB+EpB)/(EkA+EpA):.4f}")

    # 5. Główna tabela M_conf
    print("\n[5] Tabela M_conf^B(p) / M_conf^A(p):")
    print(f"  {'p':>6}  {'M_A':>12}  {'M_B':>12}  {'ratio':>10}  {'vs 206.77':>10}")
    print("  " + "-"*60)
    ratios_list = []
    for p in P_VALS:
        mA = M_conf(solA, p)
        mB = M_conf(solB, p)
        r  = mB/mA if mA != 0 else np.nan
        delta = abs(r-TARGET)/TARGET*100
        mk = " <-- BLISKIE!" if delta < 5 else (" <--" if delta < 15 else "")
        print(f"  {p:>6.2f}  {mA:>12.6f}  {mB:>12.6f}  {r:>10.4f}  {delta:>9.2f}%{mk}")
        ratios_list.append((p, mA, mB, r))

    # 6. Interpolacja p*
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
            p_star = p_arr[i] + t*(p_arr[i+1]-p_arr[i])
            print(f"  Przejście w bracket p ∈ ({p_arr[i]:.3f}, {p_arr[i+1]:.3f})")
            print(f"  Interpolacja liniowa: p* ≈ {p_star:.4f}")

            # Bisect w brackecie dla precyzji
            def ratio_fn(p_test):
                mA_ = M_conf(solA, p_test)
                mB_ = M_conf(solB, p_test)
                return mB_/mA_ - TARGET if mA_ != 0 else np.nan
            try:
                p_star_exact = brentq(ratio_fn, p_arr[i], p_arr[i+1],
                                      xtol=1e-5, maxiter=60)
                print(f"  Brentq (precyzyjne): p* = {p_star_exact:.6f}")
                p_star = p_star_exact
            except Exception as e:
                print(f"  Brentq błąd: {e}, używam interpolacji")
            break

    if p_star is None:
        max_r = np.nanmax(r_arr)
        min_r = np.nanmin(r_arr)
        print(f"  Ratio w zakresie p=[{P_VALS[0]},{P_VALS[-1]}]: [{min_r:.2f}, {max_r:.2f}]")
        if max_r < TARGET:
            print(f"  Cel 206.77 NIE osiągnięty w tym zakresie p.")
            print(f"  Potrzebne p > {P_VALS[-1]} lub inna interpretacja.")
        else:
            print(f"  Ratio przekracza 206.77, ale przejście nie wykryte (błąd?)")
    else:
        print(f"\n  *** p* = {p_star:.4f} ***")
        # Lista "ładnych" frakcji
        nice = [(n/d, f"{n}/{d}") for n in range(1,9) for d in range(1,7)
                if 0 < n/d <= P_VALS[-1]]
        nice += [(float(n), str(n)) for n in range(1,6)]
        nice = list({round(v,6):(v,s) for v,s in sorted(nice)}.values())
        nice.sort(key=lambda x: abs(x[0]-p_star))
        print("  Najbliższe 'ładne' frakcje:")
        for v, s in nice[:6]:
            print(f"    p = {v:.4f} ({s:>5}):  odch = {abs(v-p_star):.5f}  "
                  f"ratio_check = {M_conf(solB,v)/M_conf(solA,v):.4f}")

    # 7. Szczegóły dla p=3/4
    print("\n[7] Analiza szczegółowa przy p=0.75 (najbliższe p* wg p107_v1):")
    for p_t in [0.75, 0.80, 0.85, 0.90]:
        mA_ = M_conf(solA, p_t);  mB_ = M_conf(solB, p_t)
        mkA = M_conf_kin(solA, p_t); mkB = M_conf_kin(solB, p_t)
        mpA = M_conf_pot(solA, p_t); mpB = M_conf_pot(solB, p_t)
        r_  = mB_/mA_
        print(f"  p={p_t:.2f}: M_A={mA_:.5f}  M_B={mB_:.4f}  "
              f"ratio={r_:.4f}  (kin_ratio={mkB/mkA:.4f})")

    # 8. Profil φ i gęstość energii
    print("\n[8] Porównanie profili [r/a_Γ, φ_A, φ_B, ε_A, ε_B, φ_B^(2.4)/φ_A^(2.4)]:")
    eps_A = 0.5*solA.y[1]**2*(1.+ALPHA_K/np.maximum(solA.y[0],1e-10)) + \
            np.array([V_mod(float(q)) for q in solA.y[0]]) - V1
    eps_B = 0.5*solB.y[1]**2*(1.+ALPHA_K/np.maximum(solB.y[0],1e-10)) + \
            np.array([V_mod(float(q)) for q in solB.y[0]]) - V1

    print(f"  {'r/a_Γ':>8}  {'φ_A':>8}  {'φ_B':>8}  {'ε_A':>10}  "
          f"{'ε_B':>10}  {'(φB/φA)^2.4':>12}")
    checkpts = [0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
    for rc in checkpts:
        r_target = rc * AGAM
        if r_target > R_MAX:
            break
        iA = np.argmin(np.abs(solA.t - r_target))
        iB = np.argmin(np.abs(solB.t - r_target))
        phiA = solA.y[0,iA]; phiB = solB.y[0,iB]
        eA   = eps_A[iA];    eB   = eps_B[iB]
        ratio_w = (phiB/phiA)**2.4 if phiA > 0 else np.nan
        print(f"  {rc:>8.1f}  {phiA:>8.4f}  {phiB:>8.4f}  {eA:>10.4e}  "
              f"{eB:>10.4e}  {ratio_w:>12.4f}")

    # 9. Podsumowanie
    print("\n" + "="*72)
    print("PODSUMOWANIE:")
    r0 = r_arr[0]
    print(f"  p=0 (standardowe): ratio = {r0:.6f}  "
          f"(oczekiwane K*₂/K*₁ = {K2_REF/K1_REF:.6f})")

    for p_n, desc in [(1.0,"g_ij=φ²"), (2.0,"g_ij=φ⁴"), (4.0,"K=φ⁴")]:
        idx = [i for i,x in enumerate(ratios_list) if abs(x[0]-p_n)<1e-6]
        if idx:
            r_n = ratios_list[idx[0]][3]
            print(f"  p={p_n:.1f} ({desc}): ratio = {r_n:.4f}  "
                  f"(vs 206.77: {abs(r_n-TARGET)/TARGET*100:.1f}%)")

    if p_star is not None:
        print(f"\n  p* = {p_star:.4f}  → ratio = {TARGET:.2f}")
        # Sprawdź czy jest to coś sensownego fizycznie
        if abs(p_star - 1.0) < 0.05:
            print("  INTERPRETACJA: p≈1 → metryka g_ij=φ²·δ_ij")
            print("  To jest STANDARDOWA TGP metryka! Wynik byłby strukturalnie naturalny.")
        elif abs(p_star - 2./3.) < 0.05:
            print("  INTERPRETACJA: p≈2/3 → metryka g_ij=φ^(4/3)·δ_ij")
        elif abs(p_star - 0.5) < 0.05:
            print("  INTERPRETACJA: p≈1/2 → metryka g_ij=φ·δ_ij")
        else:
            print(f"  INTERPRETACJA: p={p_star:.3f} nie jest oczywistą wartością")
    else:
        print("\n  Cel 206.77 nieosiągalny w standardowej metryce TGP.")
        print("  Wniosek: M_conf z metryki φ^(2p) nie odtwarza m_μ/m_e.")
