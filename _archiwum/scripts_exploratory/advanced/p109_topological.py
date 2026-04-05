#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p109_topological.py  --  TGP v1
=====================================
Ścieżka 6: Warunek topologiczny — diagnostyka energii

MOTYWACJA:
  Ścieżki 4 i 5 dały wyniki NEGATYWNE. Szukamy innego warunku który
  wyjaśniłby stosunek mas r₂₁=206.77.

  Z poprzednich wyników wiadomo że Branch B ma ekstremalną bliskość-do-zera:
  E_kin^B ≈ 233, E_pot^B ≈ -232, E_tot^B ≈ 1.26.

  HIPOTEZA do sprawdzenia:
  Może stosunek składowych energii E_kin lub E_pot daje 206.77?

  A) E_kin^B / E_kin^A = 206.77 ?
  B) |E_pot^B| / |E_pot^A| = 206.77 ?
  C) E_kin^B / E_tot^A = 206.77 ?
  D) "Topological condition": g_top(K) = E / [4π(ψ₀-1)] - 1 = 0
     Czy ta nowa warunek self-consistency daje K*_top_A, K*_top_B
     z innym stosunkiem K*_top_B/K*_top_A?

PLAN OBLICZEŃ:
  1. Oblicz E_kin^A, E_pot^A, E_tot^A dla Branch A przy K*₁=0.010280, ψ₀^A=1.2419
  2. Oblicz E_kin^B, E_pot^B, E_tot^B dla Branch B przy K*₂=0.09999348, ψ₀^B=2.7624
  3. Raportuj wszystkie stosunki i sprawdź które zbliżają się do 206.77
  4. Skan g_top(K) = E/[4π(ψ₀-1)] - 1 dla obu gałęzi
  5. Oblicz dE/dK wzdłuż krzywej (totalna pochodna)

Data: 2026-03-27
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
N_EVAL = 800    # Wysoka dokładność dla E_kin/E_pot
RTOL   = 1e-10
ATOL   = 1e-12

AK_REF   = 8.5616
K1_REF   = 0.010280;  PSI1_REF = 1.2419
K2_REF   = 0.09999348; PSI2_REF = 2.7624

# ── ODE + energia ─────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def solve_full(K, psi0, ak):
    """Rozwiązuje ODE i zwraca (phi_end, E_kin, E_pot, E_tot, r_arr, phi_arr, dphi_arr)."""
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
        kfac = 1.0 + ak/phi
        return [dphi, dV_mod(phi)/kfac + ak*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL, t_eval=r_ev)
        if sol.t[-1] < R_MAX*0.99:
            return np.nan, np.nan, np.nan, np.nan, None, None, None
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        kfac = 1. + ak/phi
        Vv   = np.array([V_mod(float(q)) for q in phi])
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
        return sol.y[0,-1], Ek, Ep, Ek+Ep, r, phi, dphi
    except:
        return np.nan, np.nan, np.nan, np.nan, None, None, None


# ── Pomocnicze ─────────────────────────────────────────────────────────────────
def find_psi0(K, ak, psi_c, dpsi=0.30, n_pts=20, psi_min=0.8):
    """Znajdź ψ₀ s.t. φ(R_MAX)=1."""
    lo = max(psi_min, psi_c-dpsi); hi = psi_c+dpsi
    pts = np.linspace(lo, hi, n_pts)
    fv  = [(solve_full(K,p,ak)[0]-1.0) for p in pts]
    fv  = [v if np.isfinite(v) else np.nan for v in fv]
    best = None
    for i in range(len(fv)-1):
        fi, fi1 = fv[i], fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            fn = lambda p: solve_full(K,p,ak)[0]-1.0
            try:
                p0 = brentq(fn, pts[i], pts[i+1], xtol=1e-10, maxiter=60)
                if best is None or abs(p0-psi_c)<abs(best-psi_c): best=p0
            except: pass
    return best


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("="*72)
    print("p109: Diagnostyka energii — warunek topologiczny")
    print(f"  α_K={AK_REF}")
    print("="*72)

    # ── [1] Energia Branch A i Branch B przy referencji ──────────────────────
    print(f"\n[1] Energia przy K*, ψ₀* referencyjnych (α_K={AK_REF}):")
    print(f"  RTOL={RTOL:.0e}, N_EVAL={N_EVAL}")

    t0 = time.time()
    phi_end_A, Ek_A, Ep_A, Et_A, rA, phiA, dphiA = solve_full(K1_REF, PSI1_REF, AK_REF)
    phi_end_B, Ek_B, Ep_B, Et_B, rB, phiB, dphiB = solve_full(K2_REF, PSI2_REF, AK_REF)
    dt = time.time()-t0

    print(f"\n  Branch A (K*₁={K1_REF}, ψ₀^A={PSI1_REF}):")
    print(f"    φ(R_MAX) = {phi_end_A:.6f}  (oczekiwane ≈1.000)")
    print(f"    E_kin^A  = {Ek_A:.6f}")
    print(f"    E_pot^A  = {Ep_A:.6f}")
    print(f"    E_tot^A  = {Et_A:.6f}  (powinno ≈ 4π×{K1_REF:.5f} = {4*np.pi*K1_REF:.6f})")
    print(f"    g_U^A    = {Et_A/(4*np.pi*K1_REF)-1:.4e}")
    print(f"    Q_top^A  = ψ₀^A-1 = {PSI1_REF-1:.4f}")
    print(f"    E/Q_top^A = {Et_A/(PSI1_REF-1):.4f}")

    print(f"\n  Branch B (K*₂={K2_REF}, ψ₀^B={PSI2_REF}):")
    print(f"    φ(R_MAX) = {phi_end_B:.6f}  (oczekiwane ≈1.000)")
    print(f"    E_kin^B  = {Ek_B:.6f}")
    print(f"    E_pot^B  = {Ep_B:.6f}")
    print(f"    E_tot^B  = {Et_B:.6f}  (powinno ≈ 4π×{K2_REF:.5f} = {4*np.pi*K2_REF:.6f})")
    print(f"    g_U^B    = {Et_B/(4*np.pi*K2_REF)-1:.4e}")
    print(f"    Q_top^B  = ψ₀^B-1 = {PSI2_REF-1:.4f}")
    print(f"    E/Q_top^B = {Et_B/(PSI2_REF-1):.4f}")
    print(f"  (czas: {dt:.1f}s)")

    # ── [2] Stosunki energii ──────────────────────────────────────────────────
    print(f"\n[2] Stosunki energii (B/A):")
    TARGET = 206.7683

    ratios = {
        "E_kin^B / E_kin^A"         : Ek_B / Ek_A,
        "|E_pot^B| / |E_pot^A|"     : abs(Ep_B) / abs(Ep_A),
        "E_tot^B / E_tot^A"         : Et_B / Et_A,
        "E_kin^B / E_tot^A"         : Ek_B / Et_A,
        "E_kin^B / |E_pot^A|"       : Ek_B / abs(Ep_A),
        "|E_pot^B| / E_tot^A"       : abs(Ep_B) / Et_A,
        "E_kin^B / E_tot^B"         : Ek_B / Et_B,  # near-cancellation factor B
        "E_kin^A / E_tot^A"         : Ek_A / Et_A,  # near-cancellation factor A
        "(E_kin^B/E_tot^B)/(E_kin^A/E_tot^A)" : (Ek_B/Et_B)/(Ek_A/Et_A),  # ratio of nc factors
        "(K*₂/K*₁)"                : K2_REF/K1_REF,
    }

    print(f"  Target: r₂₁ = {TARGET}")
    print(f"  {'Stosunek':<45}  {'Wartość':>12}  {'%od206.77':>10}")
    for name, val in ratios.items():
        if np.isfinite(val):
            pct = abs(val-TARGET)/TARGET*100
            flag = " *** BLISKIE 206.77!" if pct < 5.0 else ""
            print(f"  {name:<45}  {val:>12.4f}  {pct:>9.3f}%{flag}")
        else:
            print(f"  {name:<45}  {'---':>12}")

    # ── [3] Q_top i nowy warunek ──────────────────────────────────────────────
    print(f"\n[3] Nowy warunek topologiczny g_top = E/[4π·Q_top] - 1:")
    Q_top_A = PSI1_REF - 1.0
    Q_top_B = PSI2_REF - 1.0
    g_top_A = Et_A / (4*np.pi*Q_top_A) - 1.0
    g_top_B = Et_B / (4*np.pi*Q_top_B) - 1.0
    print(f"  Branch A: g_top^A = {Et_A:.4f}/[4π×{Q_top_A:.4f}] - 1 = {g_top_A:.4f}")
    print(f"  Branch B: g_top^B = {Et_B:.4f}/[4π×{Q_top_B:.4f}] - 1 = {g_top_B:.4f}")
    print(f"  (g_top=0 jeśli E=4π×Q_top; Branch A potrzebuje E={4*np.pi*Q_top_A:.4f}, ma {Et_A:.4f})")
    print(f"  (g_top=0 jeśli E=4π×Q_top; Branch B potrzebuje E={4*np.pi*Q_top_B:.4f}, ma {Et_B:.4f})")

    # ── [4] Skan g_top(K) wzdłuż krzywej (z śledzeniem ψ₀) ──────────────────
    print(f"\n[4] Skan g_top(K) = E/[4π·(ψ₀(K)-1)] - 1 wzdłuż krzywej φ(R_MAX)=1:")
    print(f"    (ψ₀(K) wyznaczane przez brentq dla każdego K)")

    def compute_g_top(K, ak, psi_c, dpsi=0.40, branch="B"):
        psi_min = 1.80 if branch=="B" else 0.80
        dpsi_b = 0.38 if branch=="B" else 0.15
        p0 = find_psi0(K, ak, psi_c, dpsi=dpsi_b, n_pts=16, psi_min=psi_min)
        if p0 is None or not np.isfinite(p0): return np.nan, np.nan, np.nan, np.nan, np.nan
        _, Ek, Ep, Et, _, _, _ = solve_full(K, p0, ak)
        if not np.isfinite(Et): return np.nan, np.nan, np.nan, np.nan, p0
        Q_top = abs(p0 - 1.0)
        g_U   = Et/(4*np.pi*K) - 1.0
        g_top = Et/(4*np.pi*Q_top) - 1.0 if Q_top > 1e-6 else np.nan
        return g_U, g_top, Ek, Ep, p0

    DPSI_DK_B = 18.78
    DPSI_DK_A = 2.0   # rough estimate for Branch A

    # Branch B skan
    print(f"\n  Branch B:")
    print(f"  {'K':>10}  {'ψ₀^B':>8}  {'E_kin':>10}  {'E_pot':>10}  {'E_tot':>10}  "
          f"{'g_U':>10}  {'g_top':>10}")

    K_list_B = [0.020, 0.040, 0.060, 0.080, 0.100, 0.120, 0.140, 0.180]
    psi_c_B = PSI2_REF
    for K in K_list_B:
        psi_c_B_K = PSI2_REF + DPSI_DK_B*(K - K2_REF)
        g_U, g_top, Ek, Ep, p0 = compute_g_top(K, AK_REF, psi_c_B_K, branch="B")
        p0s = f"{p0:.4f}" if np.isfinite(p0) else "---"
        Eks = f"{Ek:.3f}" if np.isfinite(Ek) else "---"
        Eps = f"{Ep:.3f}" if np.isfinite(Ep) else "---"
        Ets = f"{(Ek+Ep) if (np.isfinite(Ek) and np.isfinite(Ep)) else np.nan:.3f}" if np.isfinite(g_U) else "---"
        gUs = f"{g_U:.4f}" if np.isfinite(g_U) else "---"
        gTs = f"{g_top:.4f}" if np.isfinite(g_top) else "---"
        print(f"  {K:>10.4f}  {p0s:>8}  {Eks:>10}  {Eps:>10}  {Ets:>10}  "
              f"{gUs:>10}  {gTs:>10}")

    # Branch A skan
    print(f"\n  Branch A:")
    print(f"  {'K':>10}  {'ψ₀^A':>8}  {'E_kin':>10}  {'E_pot':>10}  {'E_tot':>10}  "
          f"{'g_U':>10}  {'g_top':>10}")

    K_list_A = [0.005, 0.008, 0.010, 0.012, 0.015, 0.020, 0.025, 0.030]
    for K in K_list_A:
        psi_c_A = PSI1_REF + DPSI_DK_A*(K - K1_REF)
        g_U, g_top, Ek, Ep, p0 = compute_g_top(K, AK_REF, psi_c_A, branch="A")
        p0s = f"{p0:.4f}" if np.isfinite(p0) else "---"
        Eks = f"{Ek:.4f}" if np.isfinite(Ek) else "---"
        Eps = f"{Ep:.4f}" if np.isfinite(Ep) else "---"
        Ets = f"{(Ek+Ep) if (np.isfinite(Ek) and np.isfinite(Ep)) else np.nan:.4f}" if np.isfinite(g_U) else "---"
        gUs = f"{g_U:.4f}" if np.isfinite(g_U) else "---"
        gTs = f"{g_top:.4f}" if np.isfinite(g_top) else "---"
        print(f"  {K:>10.4f}  {p0s:>8}  {Eks:>10}  {Eps:>10}  {Ets:>10}  "
              f"{gUs:>10}  {gTs:>10}")

    # ── [5] Profil solitonu — gdzie leży centrum ──────────────────────────────
    if rA is not None and rB is not None:
        print(f"\n[5] Profil solitonu:")

        # R_50 — promień gdzie φ osiąga połowę drogi od ψ₀ do 1
        def r_half(r_arr, phi_arr, psi0):
            phi_mid = (psi0 + 1.0)/2.0
            for i in range(len(phi_arr)-1):
                if phi_arr[i] >= phi_mid >= phi_arr[i+1]:
                    t = (phi_mid - phi_arr[i])/(phi_arr[i+1]-phi_arr[i])
                    return r_arr[i] + t*(r_arr[i+1]-r_arr[i])
            return np.nan

        rh_A = r_half(rA, phiA, PSI1_REF)
        rh_B = r_half(rB, phiB, PSI2_REF)

        # Centroid of kinetic energy density r_kin = ∫ρ_kin·r·r²dr / ∫ρ_kin·r²dr
        def centroid(r_arr, dphi_arr, phi_arr, ak):
            kfac = 1. + ak/np.maximum(phi_arr, 1e-10)
            rho  = 0.5*dphi_arr**2 * kfac
            num  = np.trapezoid(rho*r_arr**3, r_arr)
            den  = np.trapezoid(rho*r_arr**2, r_arr)
            return num/den if den != 0 else np.nan

        r_cent_A = centroid(rA, dphiA, phiA, AK_REF)
        r_cent_B = centroid(rB, dphiB, phiB, AK_REF)

        print(f"  Branch A: R_½ = {rh_A:.3f},  r_kin_centroid = {r_cent_A:.3f}")
        print(f"  Branch B: R_½ = {rh_B:.3f},  r_kin_centroid = {r_cent_B:.3f}")
        print(f"  R_½^B/R_½^A = {rh_B/rh_A:.4f}")
        print(f"  r_kin^B/r_kin^A = {r_cent_B/r_cent_A:.4f}")

    # ── [6] Podsumowanie i wnioski ─────────────────────────────────────────────
    print(f"\n[6] PODSUMOWANIE — co daje wynik bliski 206.77?")
    print(f"  Target: 206.7683")
    print()

    checks = {
        "E_kin^B / E_kin^A" : Ek_B/Ek_A,
        "E_kin^B / E_tot^A" : Ek_B/Et_A,
        "|E_pot^B| / E_tot^A" : abs(Ep_B)/Et_A,
        "E_kin^B / |E_pot^A|" : Ek_B/abs(Ep_A),
        "E_tot^B / E_tot^A (= K*₂/K*₁)" : Et_B/Et_A,
    }
    found_any = False
    for name, val in checks.items():
        pct = abs(val-TARGET)/TARGET*100
        flag = " *** BLISKIE!!!" if pct<2.0 else ("  (ok)" if pct<10 else "")
        print(f"  {name}: {val:.4f}  (Δ={pct:.2f}%){flag}")
        if pct < 2.0: found_any = True

    if found_any:
        print("\n  *** ZNALEZIONO stosunek bliski 206.77 — hipoteza topologiczna MA SENS ***")
    else:
        print(f"\n  Żaden stosunek nie jest bliski 206.77 przy α_K={AK_REF}.")
        print(f"  WNIOSEK: Prosta wersja warunku topologicznego E_kin^B/E_kin^A ≠ 206.77.")

    print()
    print("="*72)
    print("KONIEC p109")
