#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p111_ratio_wide.py  --  TGP v1
=====================================
Ścieżka 7 (cont.): Szeroki skan E_kin^B/E_tot^B vs α_K

MOTYWACJA:
  p110 wykazał: ratio=188.1 przy α_K=8.50 i MALEJE z rosnącym α_K.
  Ale p110 stracił Branch B przy α_K=8.595 → przeskoczył na Branch A.

  Dwa obszary do zbadania:
  A) ROSNĄCY: α_K = 8.60 → 8.92 (blisko SN bifurkacji ~8.93)
     Seeds z p108_v8c: K*₂={0.08514, 0.07130, 0.05887, ..., 0.02252}
     Czy ratio rośnie przy małym K*₂ (bliżej SN)?

  B) MALEJĄCY: α_K = 8.50 → 7.70 (Branch B ma duże K*₂, duże ψ₀)
     Seed z p110 α_K=8.50: K*₂=0.121176, ψ₀=3.1291
     Czy ratio rośnie przy dużym K*₂ (niższe α_K)?

EFEKTYWNOŚĆ:
  - Stateless fn(K) z seedem z poprzedniego kroku (jak p108_v8b)
  - Tight K bracket: [K_prev×f_lo, K_prev×f_hi]
  - RTOL=1e-8, N_EVAL=350 (kompromis szybkość/dokładność)
  - n_pts=12 (mało punktów do psi₀-scan, bo dpsi małe)

Data: 2026-03-26
"""
import sys, io, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

warnings.filterwarnings('ignore')
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)  # Flush każda linia

# ── Stałe ─────────────────────────────────────────────────────────────────────
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
AGAM   = 0.040
R_MAX  = 80.0
N_EVAL = 350
RTOL   = 1e-8
ATOL   = 1e-10
TARGET = 206.7683

# ── ODE + energia ─────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def solve_full(K, psi0, ak):
    """ODE → (phi_end, E_kin, E_pot, E_tot). Szybka wersja."""
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
            return np.nan, np.nan, np.nan, np.nan
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


def find_psi0(K, ak, psi_c, dpsi=0.45, n_pts=12, psi_lo_abs=1.25):
    """Szukaj ψ₀ s.t. φ(R_MAX)=1, blisko psi_c. Zwraca best lub None."""
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


def find_K2_and_energy(ak, K_seed, psi_seed, dpsi=0.45, n_pts=12,
                       f_lo=0.80, f_hi=1.22, n_bracket=6,
                       psi_lo_abs=1.25, DPSI_DK=18.78):
    """
    Dla α_K=ak: znajdź K*₂ (g_U=0) i oblicz E_kin, E_tot.
    STATELESS: fn(K) używa stałego (K_seed, psi_seed) do ekstrapolacji ψ₀.
    """
    def fn(K):
        psi_c = psi_seed + DPSI_DK*(K - K_seed)
        psi_c = max(psi_lo_abs, psi_c)
        p0 = find_psi0(K, ak, psi_c, dpsi, n_pts, psi_lo_abs)
        if p0 is None: return np.nan
        _, Ek, Ep, Et = solve_full(K, p0, ak)
        if not np.isfinite(Et): return np.nan
        return Et/(4*np.pi*K) - 1.0

    K_lo = K_seed * f_lo
    K_hi = K_seed * f_hi

    # Skanuj n_bracket punktów K dla wykrycia znaku
    K_pts = np.linspace(K_lo, K_hi, n_bracket)
    gv = [fn(K_) for K_ in K_pts]

    found_lo, found_hi = None, None
    for i in range(len(gv)-1):
        gi, gi1 = gv[i], gv[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)): continue
        if gi*gi1 < 0:
            found_lo = K_pts[i]; found_hi = K_pts[i+1]
            break

    if found_hi is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    try:
        K2 = brentq(fn, found_lo, found_hi, xtol=1e-8, rtol=1e-8, maxiter=50)
    except:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # Oblicz finalną ψ₀ i energię przy K2
    psi_c = psi_seed + DPSI_DK*(K2 - K_seed)
    psi_c = max(psi_lo_abs, psi_c)
    psi2  = find_psi0(K2, ak, psi_c, dpsi, n_pts+4, psi_lo_abs)  # więcej punktów dla ψ₀
    if psi2 is None: return np.nan, np.nan, np.nan, np.nan, np.nan

    _, Ek, Ep, Et = solve_full(K2, psi2, ak)
    if not np.isfinite(Et): return np.nan, np.nan, np.nan, np.nan, np.nan

    gU = Et/(4*np.pi*K2) - 1.0
    return K2, psi2, Ek, Et, gU


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("="*72)
    print("p111: Szeroki skan E_kin^B/E_tot^B vs α_K")
    print(f"  RTOL={RTOL:.0e}, N_EVAL={N_EVAL}, Target={TARGET}")
    print("="*72)

    header = (f"  {'α_K':>8}  {'K*₂':>10}  {'ψ₀^B':>7}  {'E_kin':>10}  "
              f"{'E_tot':>9}  {'ratio':>8}  {'g_U':>10}  {'dt':>5}")
    sep    = "  " + "-"*85

    all_results = []   # (ak, K2, psi2, Ek, Et, ratio, scan_label)

    # ──────────────────────────────────────────────────────────────────────────
    # SKAN A: ROSNĄCY — α_K ∈ [8.60, 8.93), seeds z p108_v8c
    # ──────────────────────────────────────────────────────────────────────────
    # Seeds z p108_v8c (K*₂ at various α_K, ψ₀ estimated via DPSI_DK=18.78 from ref)
    v8c_seeds = {
        # α_K  :  (K*₂,      ψ₀_est)
        8.6116 : (0.08514,   2.763 + 18.78*(0.08514-0.10040)),  # ≈ 2.481
        8.6616 : (0.07130,   2.763 + 18.78*(0.07130-0.10040)),  # ≈ 2.217
        8.7116 : (0.05887,   2.0501),   # z v8c wprost
        8.7616 : (0.04785,   1.874),    # interpolacja
        8.8116 : (0.03811,   1.687),    # z v8c
        8.8616 : (0.02970,   1.542),    # interpolacja
        8.9116 : (0.02252,   1.4137),   # z v8c
    }
    # Oblicz brakujące ψ₀ dla v8c_seeds
    PSI2_REF = 2.7624; K2_REF = 0.10040
    for ak_key in v8c_seeds:
        K2_s, psi_s = v8c_seeds[ak_key]
        if psi_s < 1.0:  # zły estymator — zastąp
            v8c_seeds[ak_key] = (K2_s, max(1.30, PSI2_REF + 18.78*(K2_s - K2_REF)))

    ak_asc = [8.6116, 8.6616, 8.7116, 8.7616, 8.8116, 8.8616, 8.9116]

    print(f"\n[A] SKAN ROSNĄCY (α_K ∈ [8.61, 8.91], seedy z v8c):")
    print(header); print(sep)

    for ak in ak_asc:
        K_seed, psi_seed = v8c_seeds[ak]
        psi_lo = max(1.25, psi_seed - 0.60)

        # Adaptacyjny dpsi: większy gdy ψ₀ duże
        dpsi = max(0.35, 0.07*(psi_seed - 1.5) + 0.35)

        t0 = time.time()
        K2, psi2, Ek, Et, gU = find_K2_and_energy(
            ak, K_seed, psi_seed,
            dpsi=dpsi, n_pts=12,
            f_lo=0.75, f_hi=1.30,
            n_bracket=8, psi_lo_abs=psi_lo)
        dt = time.time()-t0

        if np.isfinite(K2) and Et > 0:
            ratio = Ek/Et
            all_results.append((ak, K2, psi2, Ek, Et, ratio, 'A'))
            print(f"  {ak:8.4f}  {K2:10.6f}  {psi2:7.4f}  {Ek:10.4f}  "
                  f"{Et:9.4f}  {ratio:8.4f}  {gU:10.3e}  {dt:.1f}s")
        else:
            print(f"  {ak:8.4f}  {'FAIL':>10}  {'---':>7}  {'---':>10}  "
                  f"{'---':>9}  {'---':>8}  {'---':>10}  {dt:.1f}s")

    # ──────────────────────────────────────────────────────────────────────────
    # SKAN B: MALEJĄCY — α_K = 8.50 → 7.70, seed z p110
    # ──────────────────────────────────────────────────────────────────────────
    print(f"\n[B] SKAN MALEJĄCY (α_K ∈ [8.50, 7.70], Δ=0.05, seed z p110):")
    print(header); print(sep)

    # Seed startowy z p110 α_K=8.50
    K_cur   = 0.121176
    psi_cur = 3.1291

    for ak in np.arange(8.50, 7.69, -0.05):
        dpsi = max(0.50, 0.10*(psi_cur - 1.5) + 0.40)
        dpsi = min(dpsi, 2.0)  # cap

        t0 = time.time()
        K2, psi2, Ek, Et, gU = find_K2_and_energy(
            ak, K_cur, psi_cur,
            dpsi=dpsi, n_pts=12,
            f_lo=0.82, f_hi=1.20,
            n_bracket=7, psi_lo_abs=1.30)
        dt = time.time()-t0

        if np.isfinite(K2) and Et > 0:
            ratio = Ek/Et
            all_results.append((ak, K2, psi2, Ek, Et, ratio, 'B'))
            print(f"  {ak:8.4f}  {K2:10.6f}  {psi2:7.4f}  {Ek:10.4f}  "
                  f"{Et:9.4f}  {ratio:8.4f}  {gU:10.3e}  {dt:.1f}s")
            K_cur   = K2
            psi_cur = psi2
        else:
            print(f"  {ak:8.4f}  {'FAIL':>10}  {'---':>7}  {'---':>10}  "
                  f"{'---':>9}  {'---':>8}  {'---':>10}  {dt:.1f}s")
            # Nie aktualizuj seed — spróbuj dalej z poprzednim

    # ──────────────────────────────────────────────────────────────────────────
    # ANALIZA: szukaj przejścia przez TARGET=206.77
    # ──────────────────────────────────────────────────────────────────────────
    print(f"\n[C] ANALIZA: szukanie przejścia przez ratio = {TARGET}:")

    valid = [(ak, K2, psi2, Ek, Et, r, lab)
             for (ak, K2, psi2, Ek, Et, r, lab) in all_results
             if np.isfinite(r) and Et > 0]

    # Sortuj po α_K
    valid.sort(key=lambda x: x[0])

    if not valid:
        print("  Brak danych!")
    else:
        r_vals  = [r for (_,_,_,_,_,r,_) in valid]
        ak_vals = [a for (a,_,_,_,_,_,_) in valid]
        print(f"  Zakres ratio: [{min(r_vals):.4f}, {max(r_vals):.4f}]")
        print(f"  Zakres α_K:   [{min(ak_vals):.4f}, {max(ak_vals):.4f}]")

        crossing = False
        for i in range(len(valid)-1):
            r_i  = valid[i][5]
            r_i1 = valid[i+1][5]
            if (r_i - TARGET) * (r_i1 - TARGET) < 0:
                ak_lo = valid[i][0]; ak_hi = valid[i+1][0]
                print(f"\n  *** CROSSING w [{ak_lo:.4f}, {ak_hi:.4f}]: "
                      f"ratio=[{r_i:.4f}, {r_i1:.4f}] ***")
                frac = (TARGET - r_i)/(r_i1 - r_i)
                ak_star = ak_lo + frac*(ak_hi - ak_lo)
                K2_star = valid[i][1] + frac*(valid[i+1][1]-valid[i][1])
                print(f"  Interpolacja: α_K* ≈ {ak_star:.4f}, K*₂* ≈ {K2_star:.6f}")
                crossing = True

        if not crossing:
            print(f"\n  Brak przejścia przez {TARGET} w zakresie.")
            # Pokaż najbliższe punkty
            diffs = [(abs(r-TARGET), ak, r) for ak,_,_,_,_,r,_ in valid]
            diffs.sort()
            print(f"  Najbliższe punkty do target:")
            for diff, ak_, r_ in diffs[:5]:
                print(f"    α_K={ak_:.4f}  ratio={r_:.4f}  Δ={diff:.4f} ({(r_/TARGET-1)*100:+.2f}%)")

    # ──────────────────────────────────────────────────────────────────────────
    # TABELA KOMPLETNA
    # ──────────────────────────────────────────────────────────────────────────
    print(f"\n[D] TABELA KOMPLETNA (posortowana po α_K):")
    print(f"  {'α_K':>8}  {'K*₂':>10}  {'E_kin':>10}  {'E_tot':>9}  "
          f"{'ratio':>8}  {'Δ_target%':>10}  {'scan':>5}")
    print("  " + "-"*70)
    for (ak, K2, psi2, Ek, Et, ratio, lab) in valid:
        delta = (ratio/TARGET - 1)*100
        print(f"  {ak:8.4f}  {K2:10.6f}  {Ek:10.4f}  {Et:9.4f}  "
              f"{ratio:8.4f}  {delta:+10.4f}%  {lab:>5}")

    print("\n" + "="*72)
    print("KONIEC p111")
    print()
