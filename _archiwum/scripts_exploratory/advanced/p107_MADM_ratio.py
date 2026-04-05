#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p107_MADM_ratio.py  --  TGP v1
================================
Ścieżka 4 OP-3: Masa konformalna M_conf(p) = ∫ε(r)·(φ/φ₀)^(3p)·4πr²dr

Motywacja:
  W TGP metryka przestrzenna g_ij = (φ/φ₀)^(2p)·δ_ij, więc element
  objętości √γ = (φ/φ₀)^(3p). Fizyczna energia (covariant w tej metryce) to
      M_conf(p) = ∫T₀₀·√γ·d³x = ∫ε(r)·φ(r)^(3p)·4πr²dr
  (z φ₀=1). Dla p=0: M_conf = E (standardowa energia). Pytanie: czy
  M_conf^B(p) / M_conf^A(p) = 206.77 dla jakiegoś "naturalnego" p?

  Naturalne wartości p w TGP:
    p=1  : g_ij = φ²·δ_ij  → weight φ³
    p=2  : g_ij = φ⁴·δ_ij  → weight φ⁶  (spójne z K(φ)=φ⁴?)
    p=1/2: g_ij = φ·δ_ij   → weight φ^(3/2)

Dane wejściowe (p103_v4, α_K=8.5616, a_Γ=0.040):
  Branch A: K*₁=0.010280, ψ₁=1.2417
  Branch B: K*₂=0.100276, ψ₂=2.7676

Wynik: tabela M_conf^B/M_conf^A(p) dla p=0..5; znalezienie p* gdzie ratio=206.77.

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
ALPHA_K = 8.5616
GAM     = 1.0
V1      = GAM/3.0 - GAM/4.0
LAM     = 5.501357e-06
AGAM    = 0.040
R_MAX   = 200.0     # duże r_max dla zbieżności całek
N_EVAL  = 2000      # gęsta siatka dla precyzji całkowania
RTOL    = 1e-10
ATOL    = 1e-12

# Znane rozwiązania (p103_v4)
K1_REF  = 0.010280;  PSI1_REF = 1.2417
K2_REF  = 0.100276;  PSI2_REF = 2.7676

# Zakres p (wykładnik konformalny)
P_VALS  = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0]

TARGET_RATIO = 206.77   # m_μ/m_e


# ── Potencjał i ODE ───────────────────────────────────────────────────────────

def V_mod(phi):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + LAM/6.0*(phi - 1.0)**6

def dV_mod(phi):
    return GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1.0)**5

def rhs(r, y):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K / phi
    return [dphi,
            dV_mod(phi)/kfac + ALPHA_K*dphi**2/(2.0*phi**2*kfac) - 2.0/r*dphi]


def solve_soliton(K, psi0, verbose=False):
    """Rozwiązuje ODE solitonu dla danych K i psi0. Zwraca sol."""
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    for method in ('DOP853', 'Radau'):
        try:
            sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                            method=method, rtol=RTOL, atol=ATOL, t_eval=r_ev)
            if sol.t[-1] >= R_MAX * 0.99:
                if verbose:
                    print(f"    Solver: {method}, sukces, r_end={sol.t[-1]:.1f}, "
                          f"phi_end={sol.y[0,-1]:.6f}")
                return sol
        except Exception as e:
            if verbose:
                print(f"    Solver: {method}, błąd: {e}")
    return None


def energy_density(sol):
    """Gęstość energii ε(r) = ½(φ')²·kfac + (V_mod(φ)-V1)."""
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.0 + ALPHA_K / phi
    Vv   = np.array([V_mod(float(p)) for p in phi])
    return 0.5 * dphi**2 * kfac + (Vv - V1)


def M_conf(sol, p):
    """M_conf(p) = ∫ε(r)·φ(r)^(3p)·4πr²dr = 4π·∫eps·phi^(3p)·r²dr."""
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    eps  = energy_density(sol)
    # Waga konformalna: (φ/φ₀)^(3p) z φ₀=1
    w    = phi**(3.0*p)
    integrand = eps * w * r**2
    return 4.0 * np.pi * np.trapezoid(integrand, r)


def M_kinetic(sol, p):
    """Sama część kinetyczna: M_kin(p) = 4π·∫½(φ')²·kfac·φ^(3p)·r²dr."""
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.0 + ALPHA_K / phi
    w    = phi**(3.0*p)
    return 4.0*np.pi*np.trapezoid(0.5*dphi**2*kfac*w*r**2, r)


def M_potential(sol, p):
    """Sama część potencjalna: M_pot(p) = 4π·∫(V-V1)·φ^(3p)·r²dr."""
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    Vv   = np.array([V_mod(float(q)) for q in phi])
    w    = phi**(3.0*p)
    return 4.0*np.pi*np.trapezoid((Vv - V1)*w*r**2, r)


# ── Weryfikacja samospójności ──────────────────────────────────────────────────

def check_self_consistency(sol, K, label):
    """Sprawdza warunek g_U = E/(4πK)-1 = 0."""
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    kfac = 1.0 + ALPHA_K / phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
    Vv   = np.array([V_mod(float(q)) for q in phi])
    Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
    E    = Ek + Ep
    g_U  = E/(4*np.pi*K) - 1.0
    print(f"  {label}: K={K:.6f}  psi0={sol.y[0,0]:.4f}  "
          f"phi_end={sol.y[0,-1]:.6f}  E={E:.6f}  g_U={g_U:.2e}")
    return E, g_U


# ── Refinement K* i psi0 ──────────────────────────────────────────────────────

def refine_psi0(K, psi_lo, psi_hi, tol=1e-8):
    """Dokładne ψ₀ przez brentq na warunku φ(r_max)=1."""
    def f(psi):
        sol = solve_soliton(K, psi)
        return sol.y[0,-1] - 1.0 if sol else np.nan
    try:
        return brentq(f, psi_lo, psi_hi, xtol=tol, maxiter=80)
    except Exception:
        return (psi_lo+psi_hi)/2.0


# ── Główna analiza ────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("="*72)
    print("p107: Masa konformalna M_conf(p) = ∫ε·φ^(3p)·4πr²dr")
    print(f"Pytanie: czy M_conf^B(p) / M_conf^A(p) = {TARGET_RATIO} dla ładnego p?")
    print(f"alpha_K={ALPHA_K}, a_Gamma={AGAM}, R_max={R_MAX}, N_eval={N_EVAL}")
    print("="*72)

    # ── Krok 1: Precyzyjne psi0 dla obu gałęzi ───────────────────────────────
    print("\n[1] Refinement ψ₀ (brentq, tol=1e-8):")
    psi1 = refine_psi0(K1_REF, 1.20, 1.30)
    psi2 = refine_psi0(K2_REF, 2.60, 2.90)
    print(f"  Branch A: K={K1_REF:.6f}  ψ₀={psi1:.8f}  (ref: {PSI1_REF})")
    print(f"  Branch B: K={K2_REF:.6f}  ψ₀={psi2:.8f}  (ref: {PSI2_REF})")

    # ── Krok 2: Rozwiązanie ODE ───────────────────────────────────────────────
    print("\n[2] Rozwiązanie ODE:")
    solA = solve_soliton(K1_REF, psi1, verbose=True)
    solB = solve_soliton(K2_REF, psi2, verbose=True)

    if solA is None or solB is None:
        print("BŁĄD: nie udało się rozwiązać ODE!")
        sys.exit(1)

    # ── Krok 3: Weryfikacja samospójności ─────────────────────────────────────
    print("\n[3] Weryfikacja samospójności g_U = E/(4πK)-1:")
    EA, gA = check_self_consistency(solA, K1_REF, "Branch A")
    EB, gB = check_self_consistency(solB, K2_REF, "Branch B")
    print(f"  Stosunek energii E_B/E_A = {EB/EA:.6f}  "
          f"(= K*₂/K*₁ = {K2_REF/K1_REF:.6f} — sprawdzenie)")

    # ── Krok 4: Profile φ(r) ──────────────────────────────────────────────────
    print("\n[4] Profile solitonów:")
    print(f"  Branch A: φ_max={solA.y[0].max():.4f}  "
          f"φ_min={solA.y[0].min():.4f}  r_max_used={solA.t[-1]:.1f}")
    print(f"  Branch B: φ_max={solB.y[0].max():.4f}  "
          f"φ_min={solB.y[0].min():.4f}  r_max_used={solB.t[-1]:.1f}")

    # Promienie połowiczne (r gdzie φ = 0.5*(psi0-1)+1 = (psi0+1)/2)
    for label, sol, psi_c in [("A", solA, psi1), ("B", solB, psi2)]:
        mid = (psi_c + 1.0) / 2.0
        try:
            idx = np.where(sol.y[0] <= mid)[0]
            r_half = sol.t[idx[0]] if len(idx) > 0 else np.nan
        except Exception:
            r_half = np.nan
        print(f"  Branch {label}: ψ₀={psi_c:.4f}  "
              f"φ_mid={mid:.4f}  r_half={r_half:.4f}  ({r_half/AGAM:.1f} a_Γ)")

    # ── Krok 5: Tabela M_conf(p) ─────────────────────────────────────────────
    print("\n[5] Tabela M_conf^B(p) / M_conf^A(p):")
    print(f"  {'p':>6}  {'M_A':>12}  {'M_B':>12}  {'ratio':>10}  "
          f"{'vs 206.77':>10}  {'waga_ctr':>10}")
    print("  " + "-"*70)

    ratios = []
    for p in P_VALS:
        mA = M_conf(solA, p)
        mB = M_conf(solB, p)
        r  = mB / mA if mA > 0 else np.nan
        # Przybliżony czynnik wagi w centrum: φ_A^(3p) vs φ_B^(3p)
        w_ctr = (psi2/psi1)**(3.0*p)
        delta = abs(r - TARGET_RATIO)/TARGET_RATIO*100 if np.isfinite(r) else np.nan
        marker = " <--" if delta < 5.0 else ""
        print(f"  {p:>6.2f}  {mA:>12.6f}  {mB:>12.6f}  {r:>10.4f}  "
              f"{delta:>9.2f}%{marker}  {w_ctr:>10.3f}")
        ratios.append((p, mA, mB, r))

    # ── Krok 6: Interpolacja p* gdzie ratio=206.77 ──────────────────────────
    print("\n[6] Wyznaczenie p* gdzie M_conf^B/M_conf^A = 206.77:")
    p_arr = np.array([x[0] for x in ratios])
    r_arr = np.array([x[3] for x in ratios])

    # Znajdź bracket gdzie ratio przechodzi przez 206.77
    p_star = None
    for i in range(len(r_arr)-1):
        ri, ri1 = r_arr[i], r_arr[i+1]
        if not (np.isfinite(ri) and np.isfinite(ri1)):
            continue
        if (ri - TARGET_RATIO) * (ri1 - TARGET_RATIO) < 0:
            # Interpolacja liniowa
            t = (TARGET_RATIO - ri) / (ri1 - ri)
            p_star = p_arr[i] + t * (p_arr[i+1] - p_arr[i])
            print(f"  Bracket: p ∈ ({p_arr[i]:.4f}, {p_arr[i+1]:.4f})")
            print(f"  Interpolacja liniowa: p* ≈ {p_star:.4f}")
            break

    if p_star is None:
        # Sprawdź czy ratio jest powyżej lub poniżej targetelu w całym zakresie
        if np.all(r_arr < TARGET_RATIO):
            print(f"  Ratio zawsze < {TARGET_RATIO} w zakresie p=[0,{P_VALS[-1]}]")
            print(f"  Max ratio = {np.nanmax(r_arr):.2f} przy p = {p_arr[np.nanargmax(r_arr)]:.2f}")
            print(f"  Potrzebne większe p.")
        elif np.all(r_arr > TARGET_RATIO):
            print(f"  Ratio zawsze > {TARGET_RATIO} w zakresie p=[0,{P_VALS[-1]}]")
            print(f"  p* < {P_VALS[0]}")
        else:
            print("  Nie znaleziono przejścia (NaN?)")
    else:
        # Ocena czy p* jest "ładną" liczbą
        nice = {0.5: "1/2", 1.0: "1", 1.5: "3/2", 2.0: "2", 2.5: "5/2",
                3.0: "3", 4.0: "4 (jak K=φ⁴)"}
        print(f"\n  p* = {p_star:.4f}")
        best_nice = min(nice.keys(), key=lambda x: abs(x - p_star))
        dist = abs(best_nice - p_star)
        print(f"  Najbliższa 'ładna' wartość: p = {best_nice} ({nice[best_nice]}), "
              f"odległość = {dist:.4f}")
        if dist < 0.1:
            print(f"  -> POTENCJALNIE NATURALNA wartość! (odch < 0.1)")
        elif dist < 0.3:
            print(f"  -> Blisko naturalnej (odch < 0.3), wymaga sprawdzenia")
        else:
            print(f"  -> Nie jest oczywistą naturalną wartością")

    # ── Krok 7: Szczegółowa analiza dla "naturalnych" p ─────────────────────
    print("\n[7] Analiza dla naturalnych p (1, 2, 4):")
    for p_nat, desc in [(1.0, "g_ij=φ²·δ_ij"), (2.0, "g_ij=φ⁴·δ_ij"),
                        (4.0, "K(φ)=φ⁴ analogy")]:
        mA = M_conf(solA, p_nat)
        mB = M_conf(solB, p_nat)
        mAk = M_kinetic(solA, p_nat)
        mBk = M_kinetic(solB, p_nat)
        mAp = M_potential(solA, p_nat)
        mBp = M_potential(solB, p_nat)
        r_tot  = mB/mA
        r_kin  = mBk/mAk
        r_pot  = mBp/mAp if mAp > 0 else np.nan
        print(f"\n  p = {p_nat} ({desc}):")
        print(f"    M_conf^A = {mA:.6f}  (kin={mAk:.6f}, pot={mAp:.6f})")
        print(f"    M_conf^B = {mB:.6f}  (kin={mBk:.6f}, pot={mBp:.6f})")
        print(f"    Ratio total = {r_tot:.4f}  kin = {r_kin:.4f}  pot = {r_pot:.4f}")
        print(f"    vs target {TARGET_RATIO}: delta = {abs(r_tot-TARGET_RATIO)/TARGET_RATIO*100:.2f}%")

    # ── Krok 8: Profil φ(r) — wartości centralne i zakres ───────────────────
    print("\n[8] Porównanie profili φ(r):")
    print(f"  {'r/a_Γ':>8}  {'φ_A(r)':>10}  {'φ_B(r)':>10}  {'φ_B/φ_A':>8}  {'(φ_B/φ_A)^3':>12}")
    print("  " + "-"*60)
    r_checkpoints = [0, 1, 2, 5, 10, 20, 50, 100]
    r_vals = [AGAM * rc for rc in r_checkpoints]
    for i, rc in enumerate(r_checkpoints):
        r_target = r_vals[i]
        idxA = np.argmin(np.abs(solA.t - r_target))
        idxB = np.argmin(np.abs(solB.t - r_target))
        phiA = solA.y[0, idxA]
        phiB = solB.y[0, idxB]
        ratio_phi = phiB/phiA if phiA > 0 else np.nan
        ratio_phi3 = ratio_phi**3 if np.isfinite(ratio_phi) else np.nan
        print(f"  {rc:>8.0f}  {phiA:>10.4f}  {phiB:>10.4f}  "
              f"{ratio_phi:>8.4f}  {ratio_phi3:>12.4f}")

    # ── Wnioski ───────────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("WNIOSKI:")
    r_p0 = [x[3] for x in ratios if x[0] == 0][0]
    print(f"  p=0 (standardowe E): ratio = {r_p0:.4f}  (= K*₂/K*₁, sprawdzenie)")

    found_near_natural = False
    for p_nat, desc in [(1.0, "g_ij=φ²"), (2.0, "g_ij=φ⁴")]:
        mA = M_conf(solA, p_nat)
        mB = M_conf(solB, p_nat)
        r_nat = mB/mA
        delta = abs(r_nat - TARGET_RATIO)/TARGET_RATIO*100
        print(f"  p={p_nat} ({desc}): ratio = {r_nat:.4f}  (vs {TARGET_RATIO}: delta={delta:.1f}%)")
        if delta < 10.0:
            found_near_natural = True
            print(f"    *** BLISKIE CELOWI — wymaga dalszej analizy! ***")

    if p_star is not None:
        print(f"\n  p* (ratio=206.77) ≈ {p_star:.4f}")
        nice = {0.5: "1/2", 1.0: "1", 4/3: "4/3", 1.5: "3/2", 2.0: "2",
                8/3: "8/3", 3.0: "3", 4.0: "4"}
        for pn, desc in sorted(nice.items(), key=lambda x: abs(x[0]-p_star)):
            print(f"    p={pn:.4f} ({desc}): odch = {abs(pn-p_star):.4f}")
            if abs(pn-p_star) > 0.5:
                break
    else:
        max_r = np.nanmax(r_arr)
        max_p = p_arr[np.nanargmax(r_arr)]
        print(f"\n  Maksymalne ratio w zakresie p=[0,{P_VALS[-1]}]: "
              f"{max_r:.2f} przy p={max_p:.2f}")
        print(f"  Cel 206.77 NIE osiągnięty — potrzebne p > {P_VALS[-1]} lub "
              f"inna interpretacja.")
