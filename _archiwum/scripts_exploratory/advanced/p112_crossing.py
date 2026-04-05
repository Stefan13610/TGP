#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p112_crossing.py  --  TGP v1
=====================================
Ścieżka 7 (finale): Precyzyjne wyznaczenie α_K* s.t. E_kin^B/E_tot^B = 206.77

Z p111: ratio=205.60 przy α_K=7.70 (Δ=-0.56%).
Trend: ~37.4 per unit α_K → crossing w α_K ≈ 7.669.

PLAN:
  1. Fine scan α_K ∈ [7.70, 7.60], krok=0.010
  2. Interpolacja liniowa na crossing
  3. Weryfikacja precyzyjna przy α_K*

Data: 2026-03-26
"""
import sys, io, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)

# ── Stałe ─────────────────────────────────────────────────────────────────────
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
AGAM   = 0.040
R_MAX  = 80.0
N_EVAL = 400
RTOL   = 1e-8
ATOL   = 1e-10
TARGET = 206.7683

def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def solve_full(K, psi0, ak):
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
        if sol.t[-1] < R_MAX*0.99: return np.nan, np.nan, np.nan, np.nan
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        kfac = 1. + ak/phi
        Vv   = np.array([V_mod(float(q)) for q in phi])
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Ep   = 4*np.pi*np.trapezoid((Vv - V1)*r**2, r)
        return sol.y[0,-1], Ek, Ep, Ek+Ep
    except:
        return np.nan, np.nan, np.nan, np.nan


def find_psi0(K, ak, psi_c, dpsi=1.0, n_pts=12, psi_lo_abs=1.25):
    lo  = max(psi_lo_abs, psi_c - dpsi)
    hi  = psi_c + dpsi
    pts = np.linspace(lo, hi, n_pts)
    fv  = [(solve_full(K, p, ak)[0] - 1.0) for p in pts]
    fv  = [v if np.isfinite(v) else np.nan for v in fv]
    best = None
    for i in range(len(fv)-1):
        fi, fi1 = fv[i], fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            try:
                p0 = brentq(lambda p: solve_full(K, p, ak)[0]-1.0,
                            pts[i], pts[i+1], xtol=1e-9, maxiter=50)
                if best is None or abs(p0-psi_c) < abs(best-psi_c):
                    best = p0
            except: pass
    return best


def find_K2_energy(ak, K_seed, psi_seed, DPSI_DK=18.78,
                   dpsi=1.5, n_pts=12, psi_lo_abs=1.30,
                   f_lo=0.85, f_hi=1.18, n_bracket=6):
    """Stateless brentq-in-K."""
    def fn(K):
        psi_c = max(psi_lo_abs, psi_seed + DPSI_DK*(K - K_seed))
        p0 = find_psi0(K, ak, psi_c, dpsi, n_pts, psi_lo_abs)
        if p0 is None: return np.nan
        _, Ek, Ep, Et = solve_full(K, p0, ak)
        if not np.isfinite(Et): return np.nan
        return Et/(4*np.pi*K) - 1.0

    K_lo = K_seed * f_lo; K_hi = K_seed * f_hi
    K_pts = np.linspace(K_lo, K_hi, n_bracket)
    gv = [fn(K_) for K_ in K_pts]

    found_lo, found_hi = None, None
    for i in range(len(gv)-1):
        gi, gi1 = gv[i], gv[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)): continue
        if gi*gi1 < 0:
            found_lo = K_pts[i]; found_hi = K_pts[i+1]; break

    if found_hi is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    try:
        K2 = brentq(fn, found_lo, found_hi, xtol=1e-8, rtol=1e-8, maxiter=50)
    except:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    psi_c2 = max(psi_lo_abs, psi_seed + DPSI_DK*(K2 - K_seed))
    psi2 = find_psi0(K2, ak, psi_c2, dpsi, n_pts+4, psi_lo_abs)
    if psi2 is None: return np.nan, np.nan, np.nan, np.nan, np.nan

    _, Ek, Ep, Et = solve_full(K2, psi2, ak)
    if not np.isfinite(Et): return np.nan, np.nan, np.nan, np.nan, np.nan

    gU = Et/(4*np.pi*K2) - 1.0
    return K2, psi2, Ek, Et, gU


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("="*72)
    print("p112: Fine scan α_K ∈ [7.70, 7.60] — crossing E_kin^B/E_tot^B=206.77")
    print(f"  RTOL={RTOL:.0e}, N_EVAL={N_EVAL}, Target={TARGET}")
    print("="*72)

    # Seed z p111 α_K=7.70
    K_cur   = 0.550812
    psi_cur = 11.4013

    print(f"\n[1] Fine scan (seed: K*₂={K_cur:.6f}, ψ₀={psi_cur:.4f} @ α_K=7.70):")
    print(f"  {'α_K':>8}  {'K*₂':>10}  {'ψ₀^B':>8}  {'E_kin':>12}  "
          f"{'E_tot':>9}  {'ratio':>8}  {'Δ%':>7}  {'g_U':>10}  {'dt':>5}")
    print("  " + "-"*90)

    results = []

    for ak in np.arange(7.70, 7.59, -0.01):
        # Adaptacyjny dpsi: duże przy dużym ψ₀
        dpsi = max(1.2, 0.15*(psi_cur-1.5)+0.8)
        dpsi = min(dpsi, 3.5)

        t0 = time.time()
        K2, psi2, Ek, Et, gU = find_K2_energy(
            ak, K_cur, psi_cur,
            dpsi=dpsi, n_pts=12,
            psi_lo_abs=1.30,
            f_lo=0.90, f_hi=1.12,
            n_bracket=6)
        dt = time.time()-t0

        if np.isfinite(K2) and Et > 0:
            ratio = Ek/Et
            delta = (ratio/TARGET - 1)*100
            results.append((ak, K2, psi2, Ek, Et, ratio))
            print(f"  {ak:8.4f}  {K2:10.6f}  {psi2:8.4f}  {Ek:12.4f}  "
                  f"{Et:9.4f}  {ratio:8.4f}  {delta:+7.4f}%  {gU:10.3e}  {dt:.1f}s")
            K_cur = K2; psi_cur = psi2
        else:
            print(f"  {ak:8.4f}  {'FAIL':>10}  {'---':>8}  {'---':>12}  "
                  f"{'---':>9}  {'---':>8}  {'---':>7}  {'---':>10}  {dt:.1f}s")

    # Interpolacja
    print(f"\n[2] Interpolacja crossing:")
    valid = [(ak, K2, psi2, Ek, Et, r) for (ak, K2, psi2, Ek, Et, r) in results
             if np.isfinite(r)]

    crossing = False
    for i in range(len(valid)-1):
        r_i  = valid[i][5]; r_i1 = valid[i+1][5]
        if (r_i - TARGET)*(r_i1 - TARGET) < 0:
            ak_lo = valid[i][0]; ak_hi = valid[i+1][0]
            frac  = (TARGET - r_i)/(r_i1 - r_i)
            ak_star  = ak_lo + frac*(ak_hi - ak_lo)
            K2_star  = valid[i][1] + frac*(valid[i+1][1]-valid[i][1])
            psi_star = valid[i][2] + frac*(valid[i+1][2]-valid[i][2])
            Ek_star  = valid[i][3] + frac*(valid[i+1][3]-valid[i][3])
            Et_star  = valid[i][4] + frac*(valid[i+1][4]-valid[i][4])
            print(f"  Crossing [{ak_lo:.3f}, {ak_hi:.3f}]: "
                  f"ratio=[{r_i:.4f}, {r_i1:.4f}]")
            print(f"  α_K* ≈ {ak_star:.5f}")
            print(f"  K*₂* ≈ {K2_star:.6f}   (E_tot ≈ 4π×K*₂ = {4*np.pi*K2_star:.4f})")
            print(f"  ψ₀*  ≈ {psi_star:.4f}")
            print(f"  E_kin* ≈ {Ek_star:.4f}")
            print(f"  E_tot* ≈ {Et_star:.4f}")
            print(f"  Check: E_kin*/E_tot* ≈ {Ek_star/Et_star:.4f}  (target={TARGET})")
            crossing = True

    if not crossing:
        print("  Brak crossing w zakresie.")
        if valid:
            r_vals = [r for (_,_,_,_,_,r) in valid]
            ak_vals = [a for (a,_,_,_,_,_) in valid]
            print(f"  min ratio={min(r_vals):.4f} @ α_K={ak_vals[np.argmin(r_vals)]:.4f}")
            print(f"  max ratio={max(r_vals):.4f} @ α_K={ak_vals[np.argmax(r_vals)]:.4f}")

    # [3] Kontekst: Branch A przy α_K*
    if crossing and valid:
        print(f"\n[3] Kontekst fizyczny przy α_K*={ak_star:.4f}:")
        print(f"  K*₁ (Branch A) ≈ ??? (wymaga osobnego obliczenia)")
        print(f"  K*₂ (Branch B) ≈ {K2_star:.6f}")
        # Estymacja K*₁ — Branch A przy tym α_K
        # K*₁ zmienia się wolno z α_K; estymacja: ~0.012 przy α_K=7.67
        K1_est = 0.010280 + (ak_star - 8.5616) * (-0.0005)  # gruba estymacja
        print(f"  K*₁ ≈ {K1_est:.6f}  (gruba estymacja)")
        print(f"  K*₂/K*₁ ≈ {K2_star/K1_est:.2f}  (nie 206.77!)")
        print(f"\n  INTERPRETACJA:")
        print(f"  - E_kin^B/E_tot^B = {TARGET} przy α_K*={ak_star:.4f}")
        print(f"  - ALE K*₂/K*₁ ≈ {K2_star/K1_est:.1f} ≠ 206.77")
        print(f"  - Warunek energetyczny ≠ warunek masowy")
        print(f"  - Hipoteza wymaga fizycznej reinterpretacji")

    print("\n" + "="*72)
    print("KONIEC p112")
    print()
