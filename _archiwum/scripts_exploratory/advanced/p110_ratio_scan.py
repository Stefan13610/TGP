#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p110_ratio_scan.py  --  TGP v1
=====================================
Ścieżka 7: Szukanie α_K s.t. E_kin^B / E_tot^B = 206.77

MOTYWACJA:
  p109 wykazał: przy α_K=8.5616 E_kin^B/E_tot^B = 187.4 (Δ=9.35% od 206.77).
  To jest najbliższy stosunek energii znaleziony do tej pory.

  HIPOTEZA: Istnieje α_K przy którym E_kin^B(α_K) / E_tot^B(α_K) = 206.77 DOKŁADNIE.

  Obserwacja:
    E_tot^B ≈ 4π·K*₂  (warunek self-consistency)
    E_kin^B maleje wolno gdy K*₂ rośnie (bo K*₂ maleje gdy α_K rośnie → E_kin rośnie)
    Więc E_kin^B/E_tot^B rośnie gdy α_K rośnie (ponad 8.5616)

  Szacunek: Target ratio=206.77 → E_tot=E_kin/206.77≈233/206.77≈1.127
    → K*₂≈1.127/(4π)≈0.0897 → między α_K=8.5616 (K*₂=0.1004)
      a α_K=8.6116 (K*₂=0.0851)

PLAN:
  1. Skan α_K ∈ [8.50, 8.70], krok 0.005
  2. Dla każdego α_K: znajdź K*₂ (brentq stateless), oblicz E_kin^B, E_tot^B
  3. Wyznacz α_K* s.t. ratio = 206.77 (brentq po α_K)
  4. Weryfikacja numeryczna przy α_K*

Data: 2026-03-26
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

# Dokładność: wyższa niż p108 (1e-7), niższa niż p109 (1e-10) dla równowagi
N_EVAL = 500
RTOL   = 1e-9
ATOL   = 1e-11

# Punkt referencyjny (α_K=8.5616)
AK_REF   = 8.5616
K2_REF   = 0.10040394
PSI2_REF = 2.7624       # ψ₀^B przy referencji

# Pochodna dψ₀^B/dK (z p108_v8b)
DPSI_DK_B = 18.78

# Cel
TARGET_RATIO = 206.7683

# ── ODE + energia ─────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def solve_full(K, psi0, ak):
    """Rozwiązuje ODE, zwraca (phi_end, E_kin, E_pot, E_tot)."""
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
        Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
        return sol.y[0,-1], Ek, Ep, Ek+Ep
    except:
        return np.nan, np.nan, np.nan, np.nan


def find_psi0_B(K, ak, psi_c, dpsi=0.45, n_pts=20):
    """Znajdź ψ₀^B s.t. φ(R_MAX)=1, blisko psi_c."""
    lo = max(1.05, psi_c - dpsi)
    hi = psi_c + dpsi
    pts = np.linspace(lo, hi, n_pts)
    fv  = [(solve_full(K, p, ak)[0] - 1.0) for p in pts]
    fv  = [v if np.isfinite(v) else np.nan for v in fv]
    best = None
    for i in range(len(fv)-1):
        fi, fi1 = fv[i], fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            fn_ = lambda p: (solve_full(K, p, ak)[0] - 1.0)
            try:
                p0 = brentq(fn_, pts[i], pts[i+1], xtol=1e-10, maxiter=60)
                if best is None or abs(p0-psi_c) < abs(best-psi_c):
                    best = p0
            except: pass
    return best


def compute_B_at_alpha(ak, K_seed, psi_seed,
                       K_lo=0.01, K_hi=0.25,
                       dpsi=0.45, n_pts=20):
    """
    Dla danego α_K znajdź K*₂ (g_U=0) i oblicz E_kin^B, E_tot^B.
    Zwraca (K2, psi2, Ek, Et, gU).
    Używa STATELESS extrapolacji: psi_c = psi_seed + DPSI_DK_B*(K-K_seed).
    """
    # Funkcja g_U(K) — STATELESS (brentq może ją wołać w dowolnej kolejności)
    def fn(K):
        psi_c = psi_seed + DPSI_DK_B*(K - K_seed)
        psi_c = max(1.10, psi_c)  # nie poniżej min
        p0 = find_psi0_B(K, ak, psi_c, dpsi, n_pts)
        if p0 is None:
            return np.nan
        _, Ek, Ep, Et = solve_full(K, p0, ak)
        if not np.isfinite(Et):
            return np.nan
        return Et/(4*np.pi*K) - 1.0

    # Sprawdź endpoints
    g_lo = fn(K_lo)
    g_hi = fn(K_hi)

    if not (np.isfinite(g_lo) and np.isfinite(g_hi)):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # Szukaj przedziału ze zmianą znaku
    found_lo, found_hi = None, None
    for K_test in np.linspace(K_lo, K_hi, 40):
        g = fn(K_test)
        if not np.isfinite(g): continue
        if found_lo is None:
            found_lo = K_test; g_found_lo = g
        else:
            if g_found_lo * g < 0:
                found_hi = K_test; g_found_hi = g
                break
            # Aktualizuj tylko jeśli ten sam znak
            if g * g_found_lo >= 0:
                found_lo = K_test; g_found_lo = g

    if found_hi is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # brentq na znalezionym przedziale
    try:
        K2 = brentq(fn, found_lo, found_hi, xtol=1e-8, rtol=1e-8, maxiter=50)
    except Exception as e:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # Oblicz pełną energię przy K2
    psi_c = psi_seed + DPSI_DK_B*(K2 - K_seed)
    psi_c = max(1.10, psi_c)
    psi2 = find_psi0_B(K2, ak, psi_c, dpsi, n_pts)
    if psi2 is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    _, Ek, Ep, Et = solve_full(K2, psi2, ak)
    if not np.isfinite(Et):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    gU = Et/(4*np.pi*K2) - 1.0
    return K2, psi2, Ek, Et, gU


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("="*72)
    print("p110: Szukanie α_K* s.t. E_kin^B / E_tot^B = 206.77")
    print(f"  RTOL={RTOL:.0e}, N_EVAL={N_EVAL}")
    print(f"  Target ratio = {TARGET_RATIO}")
    print("="*72)

    # ── [1] Weryfikacja przy α_K_REF ─────────────────────────────────────────
    print(f"\n[1] Weryfikacja przy α_K={AK_REF} (seed K*₂={K2_REF}, ψ₀={PSI2_REF}):")
    t0 = time.time()
    K2v, psi2v, Ekv, Etv, gUv = compute_B_at_alpha(
        AK_REF, K2_REF, PSI2_REF, K_lo=0.05, K_hi=0.20)
    dt = time.time()-t0
    if np.isfinite(K2v):
        print(f"  K*₂={K2v:.6f}, ψ₀={psi2v:.4f}, E_kin={Ekv:.4f}, "
              f"E_tot={Etv:.4f}, g_U={gUv:.3e}, ratio={Ekv/Etv:.4f}  ({dt:.1f}s)")
    else:
        print(f"  FAIL ({dt:.1f}s)")

    # ── [2] Skan α_K ─────────────────────────────────────────────────────────
    ak_vals  = np.arange(8.50, 8.71, 0.005)
    results  = []  # (ak, K2, psi2, Ek, Et, ratio)

    print(f"\n[2] Skan α_K ∈ [{ak_vals[0]:.3f}, {ak_vals[-1]:.3f}], krok=0.005:")
    print(f"  {'α_K':>8}  {'K*₂':>10}  {'ψ₀^B':>7}  {'E_kin':>10}  "
          f"{'E_tot':>9}  {'ratio':>8}  {'g_U':>10}  {'dt':>5}")
    print("  " + "-"*85)

    # Seed startowy (z referencji)
    k_seed_cur  = K2_REF
    psi_seed_cur = PSI2_REF

    for ak in ak_vals:
        t0 = time.time()
        # Dynamiczny bracket: K*₂ maleje gdy α_K rośnie
        K_lo = max(0.01, k_seed_cur - 0.08)
        K_hi = min(0.30, k_seed_cur + 0.08)

        K2, psi2, Ek, Et, gU = compute_B_at_alpha(
            ak, k_seed_cur, psi_seed_cur, K_lo=K_lo, K_hi=K_hi)
        dt = time.time()-t0

        if np.isfinite(K2) and Et > 0:
            ratio = Ek/Et
            results.append((ak, K2, psi2, Ek, Et, ratio))
            print(f"  {ak:8.4f}  {K2:10.6f}  {psi2:7.4f}  {Ek:10.4f}  "
                  f"{Et:9.4f}  {ratio:8.4f}  {gU:10.3e}  {dt:.1f}s")
            # Aktualizuj seed dla następnego kroku (tracking Branch B)
            k_seed_cur   = K2
            psi_seed_cur = psi2
        else:
            print(f"  {ak:8.4f}  {'---':>10}  {'---':>7}  {'---':>10}  "
                  f"{'---':>9}  {'---':>8}  {'---':>10}  {dt:.1f}s")

    # ── [3] Interpolacja: gdzie ratio = 206.77 ────────────────────────────────
    print(f"\n[3] Interpolacja: gdzie ratio = {TARGET_RATIO}:")
    valid = [(ak, K2, psi2, Ek, Et, r) for (ak, K2, psi2, Ek, Et, r) in results if np.isfinite(r)]

    # Znajdź przedział gdzie ratio przechodzi przez TARGET_RATIO
    crossing_found = False
    for i in range(len(valid)-1):
        r_i  = valid[i][5]
        r_i1 = valid[i+1][5]
        if (r_i - TARGET_RATIO) * (r_i1 - TARGET_RATIO) < 0:
            ak_lo = valid[i][0];  ak_hi = valid[i+1][0]
            r_lo  = r_i;          r_hi  = r_i1
            print(f"  Crossing w [{ak_lo:.4f}, {ak_hi:.4f}]: "
                  f"ratio=[{r_lo:.4f}, {r_hi:.4f}]")
            # Interpolacja liniowa
            frac = (TARGET_RATIO - r_lo) / (r_hi - r_lo)
            ak_star = ak_lo + frac*(ak_hi - ak_lo)
            K2_star = valid[i][1] + frac*(valid[i+1][1]-valid[i][1])
            print(f"  Interpolacja liniowa: α_K* ≈ {ak_star:.5f}, K*₂* ≈ {K2_star:.6f}")
            crossing_found = True

    if not crossing_found:
        print("  Brak przejścia przez 206.77 w przeskanowanym zakresie.")
        print("  Ekstrema ratio:")
        if valid:
            r_vals = [r for (_,_,_,_,_,r) in valid]
            ak_vals_ = [a for (a,_,_,_,_,_) in valid]
            i_min = np.argmin(r_vals); i_max = np.argmax(r_vals)
            print(f"    min ratio = {r_vals[i_min]:.4f} przy α_K={ak_vals_[i_min]:.4f}")
            print(f"    max ratio = {r_vals[i_max]:.4f} przy α_K={ak_vals_[i_max]:.4f}")

    # ── [4] Weryfikacja przy α_K* (precyzyjna) ────────────────────────────────
    if crossing_found:
        print(f"\n[4] Weryfikacja precyzyjna przy α_K*={ak_star:.5f}:")

        # Seed z sąsiednich punktów
        ak_idx = next(j for j,(a,_,_,_,_,_) in enumerate(valid) if a >= ak_lo)
        k_seed_v  = valid[ak_idx-1][1]  if ak_idx > 0 else K2_REF
        psi_seed_v = valid[ak_idx-1][2] if ak_idx > 0 else PSI2_REF

        t0 = time.time()
        K2f, psi2f, Ekf, Etf, gUf = compute_B_at_alpha(
            ak_star, k_seed_v, psi_seed_v,
            K_lo=max(0.01, k_seed_v-0.05),
            K_hi=min(0.25, k_seed_v+0.05))
        dt = time.time()-t0

        if np.isfinite(K2f) and Etf > 0:
            ratio_f = Ekf/Etf
            delta_pct = (ratio_f/TARGET_RATIO - 1)*100
            print(f"  α_K*  = {ak_star:.5f}")
            print(f"  K*₂   = {K2f:.8f}")
            print(f"  ψ₀^B  = {psi2f:.6f}")
            print(f"  E_kin = {Ekf:.6f}")
            print(f"  E_tot = {Etf:.6f}")
            print(f"  g_U   = {gUf:.3e}")
            print(f"  ratio = E_kin/E_tot = {ratio_f:.6f}")
            print(f"  Δ od {TARGET_RATIO} = {delta_pct:+.4f}%")
            print(f"  4π·K*₂ = {4*np.pi*K2f:.6f}  (powinno ≈ E_tot)")
            print(f"  ({dt:.1f}s)")
        else:
            print(f"  FAIL ({dt:.1f}s)")

    # ── [5] Tabela ratio vs α_K (podsumowanie) ───────────────────────────────
    print(f"\n[5] TABELA: E_kin^B/E_tot^B vs α_K:")
    print(f"  {'α_K':>8}  {'K*₂':>10}  {'E_kin':>10}  {'E_tot':>9}  "
          f"{'ratio':>8}  {'Δ_target%':>10}")
    print("  " + "-"*65)
    for (ak, K2, psi2, Ek, Et, ratio) in results:
        if Et > 0:
            delta = (ratio/TARGET_RATIO - 1)*100
            print(f"  {ak:8.4f}  {K2:10.6f}  {Ek:10.4f}  {Et:9.4f}  "
                  f"{ratio:8.4f}  {delta:+10.4f}%")

    print("\n" + "="*72)
    print("KONIEC p110")
    print()
