#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v8b_fix.py  --  TGP v1
=====================================
Ścieżka 5: argmax K*₂/K*₁ — BRENTQ w K, POPRAWIONY tracking

BUG w v8:
  fn(K) aktualizował psi_ref[0] → brentq wywołuje fn() nietonotonicznie
  → okno ψ₀ dryfuje na złą gałąź → spurious K*₂=0.10185, g_U=-1.08, ψ₀=9.47

FIX:
  fn(K) używa CZYSTEJ ekstrapolacji z psi_seed (bez aktualizacji stanu):
    psi_c = psi_seed + DPSI_DK_B*(K - K_seed)
  + Weryfikacja: |g_U(K*)| < 0.05 i 1.80 < ψ₀ < 5.0

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
N_EVAL = 500
RTOL   = 1e-7
ATOL   = 1e-9

AK_REF   = 8.5616
K1_REF   = 0.010280;  PSI1_REF = 1.2419
K2_REF   = 0.09999348; PSI2_REF = 2.7624
DPSI_DK_B = 18.78
PSI_B_MIN  = 1.80
PSI_B_MAX  = 5.50

# ── ODE ───────────────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def solve_one(K, psi0, ak):
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
        if sol.t[-1] < R_MAX*0.99: return np.nan, np.nan
        r    = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        kfac = 1. + ak/phi
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Vv   = np.array([V_mod(float(q)) for q in phi])
        Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
        return sol.y[0,-1], Ek+Ep
    except: return np.nan, np.nan


def find_psi0_near(K, ak, psi_c, dpsi, n_pts, psi_min=0.80):
    """Znajdź ψ₀ w [max(psi_min, psi_c-dpsi), psi_c+dpsi] s.t. φ(R_MAX)=1."""
    lo = max(psi_min, psi_c - dpsi)
    hi = psi_c + dpsi
    if lo >= hi: return np.nan, np.nan
    pts = np.linspace(lo, hi, n_pts)
    fv  = []
    for p in pts:
        pe, _ = solve_one(K, p, ak)
        fv.append(pe - 1.0 if np.isfinite(pe) else np.nan)
    candidates = []
    for i in range(len(fv)-1):
        fi, fi1 = fv[i], fv[i+1]
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            try:
                fn = lambda p: solve_one(K, p, ak)[0] - 1.0
                p0 = brentq(fn, pts[i], pts[i+1], xtol=1e-9, maxiter=40)
                candidates.append(p0)
            except:
                candidates.append(pts[i] - fv[i]*(pts[i+1]-pts[i])/(fv[i+1]-fv[i]))
    if not candidates: return np.nan, np.nan
    best = min(candidates, key=lambda p: abs(p-psi_c))
    _, E = solve_one(K, best, ak)
    return best, E


def gU_at_K(K, ak, psi_c, dpsi, n_pts, psi_min=0.80):
    p0, E = find_psi0_near(K, ak, psi_c, dpsi, n_pts, psi_min=psi_min)
    if not (np.isfinite(p0) and np.isfinite(E) and K > 0): return np.nan, psi_c
    return E/(4*np.pi*K) - 1.0, p0


def find_K_star(ak, K_seed, psi_seed,
                dpsi, n_pts, psi_min,
                K_frac=0.14):
    """
    Brentq w K na [K_seed*(1-K_frac), K_seed*(1+K_frac)].
    KLUCZOWA POPRAWKA: fn(K) używa CZYSTEJ ekstrapolacji z (K_seed, psi_seed).
    Bez aktualizacji stanu → bezpieczne dla brentq.
    """
    K_lo = K_seed * (1 - K_frac)
    K_hi = K_seed * (1 + K_frac)
    is_B = (psi_min >= 1.5)

    def psi_center(K):
        """Stała ekstrapolacja z punktu referencyjnego (K_seed, psi_seed)."""
        if is_B:
            return psi_seed + DPSI_DK_B*(K - K_seed)
        return psi_seed  # Branch A: ψ₀ stabilne

    # Sprawdź znaki na końcach
    g_lo, _ = gU_at_K(K_lo, ak, psi_center(K_lo), dpsi, n_pts, psi_min)
    g_hi, _ = gU_at_K(K_hi, ak, psi_center(K_hi), dpsi, n_pts, psi_min)

    # Fallback: rozszerz 2×
    if not (np.isfinite(g_lo) and np.isfinite(g_hi) and g_lo*g_hi < 0):
        K_lo2 = K_seed*(1-2*K_frac)
        K_hi2 = K_seed*(1+2*K_frac)
        dpsi2 = dpsi * 1.25
        g_lo2, _ = gU_at_K(K_lo2, ak, psi_center(K_lo2), dpsi2, n_pts, psi_min)
        g_hi2, _ = gU_at_K(K_hi2, ak, psi_center(K_hi2), dpsi2, n_pts, psi_min)
        if np.isfinite(g_lo2) and np.isfinite(g_hi2) and g_lo2*g_hi2 < 0:
            K_lo, K_hi, g_lo, dpsi = K_lo2, K_hi2, g_lo2, dpsi2
        else:
            return np.nan, psi_seed, np.nan

    # Brentq — CZYSTA ekstrapolacja, bez aktualizacji stanu
    def fn(K):
        g, _ = gU_at_K(K, ak, psi_center(K), dpsi*0.85, n_pts, psi_min)
        return g if np.isfinite(g) else g_lo  # fallback zachowuje znak z K_lo

    try:
        K_star = brentq(fn, K_lo, K_hi, xtol=1e-8, maxiter=60)
    except Exception:
        return np.nan, psi_seed, np.nan

    # Weryfikacja końcowa
    psi_f = psi_center(K_star)
    g_f, p_f = gU_at_K(K_star, ak, psi_f, dpsi*0.75, n_pts, psi_min)

    # SANITY CHECK: musi być bliskie 0 i ψ₀ w rozsądnym zakresie
    if not np.isfinite(g_f):
        return np.nan, psi_seed, np.nan
    if abs(g_f) > 0.08:
        return np.nan, psi_seed, np.nan
    if is_B and not (PSI_B_MIN < p_f < PSI_B_MAX):
        return np.nan, psi_seed, np.nan

    return K_star, p_f, g_f


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':

    print("="*72)
    print("p108_v8b: Branch B tracking — POPRAWIONY brentq (czysta ekstrapolacja)")
    print(f"  RTOL={RTOL:.0e}, N_EVAL={N_EVAL}")
    print(f"  Referencja: α_K={AK_REF}, K*₂={K2_REF:.8f}")
    print("="*72)
    sys.stdout.flush()

    # Timing
    t0 = time.time()
    _ = solve_one(K2_REF, PSI2_REF, AK_REF)
    print(f"\n[TIMING] 1 ODE: {time.time()-t0:.3f}s")
    sys.stdout.flush()

    # ── Weryfikacja przy AK_REF ───────────────────────────────────────────────
    print(f"\n[0] Weryfikacja Branch B przy α_K={AK_REF}:")
    sys.stdout.flush()

    for K_test, lbl in [(0.090, "ujemne"), (0.110, "dodatnie")]:
        pc = PSI2_REF + DPSI_DK_B*(K_test - K2_REF)
        g_, p_ = gU_at_K(K_test, AK_REF, pc, dpsi=0.35, n_pts=14, psi_min=PSI_B_MIN)
        print(f"  K={K_test:.3f}: psi_c={pc:.4f} → ψ₀={p_:.4f}, g_U={g_:+.4f}  (oczekiwane {lbl})")
        sys.stdout.flush()

    # K*₂ referencyjne
    K2_v, psi2_v, gU2_v = find_K_star(
        AK_REF, K2_REF, PSI2_REF,
        dpsi=0.35, n_pts=14, psi_min=PSI_B_MIN, K_frac=0.14)

    if np.isfinite(K2_v):
        d2 = abs(K2_v - K2_REF)/K2_REF*100
        ok2 = d2 < 2.0 and abs(gU2_v) < 0.05
        print(f"\n  K*₂ = {K2_v:.8f}  (ref={K2_REF:.8f})")
        print(f"       Δ={d2:.3f}%  |g_U|={abs(gU2_v):.2e}  ψ₀={psi2_v:.4f}  {'✓' if ok2 else '✗ BAD'}")
        K2_s, PSI2_s = (K2_v, psi2_v) if ok2 else (K2_REF, PSI2_REF)
    else:
        print(f"\n  ✗ K*₂ nie znaleziono — używam ref")
        K2_s, PSI2_s = K2_REF, PSI2_REF

    # K*₁ referencyjne
    K1_v, psi1_v, gU1_v = find_K_star(
        AK_REF, K1_REF, PSI1_REF,
        dpsi=0.15, n_pts=10, psi_min=0.80, K_frac=0.20)
    if np.isfinite(K1_v) and abs(gU1_v) < 0.05:
        print(f"  K*₁ = {K1_v:.8f}  Δ={abs(K1_v-K1_REF)/K1_REF*100:.3f}%  "
              f"|g_U|={abs(gU1_v):.2e}  ψ₀={psi1_v:.4f}  ✓")
        K1_s, PSI1_s = K1_v, psi1_v
    else:
        print(f"  ✗ K*₁ nie znaleziono — używam ref")
        K1_s, PSI1_s = K1_REF, PSI1_REF
    sys.stdout.flush()

    # ── Skan ──────────────────────────────────────────────────────────────────
    DAK = 0.05

    def build_ak(start, stop, step):
        vals = [start]
        v = start
        while True:
            v = round(v + step, 8)
            if step>0 and v > stop+1e-9: break
            if step<0 and v < stop-1e-9: break
            vals.append(v)
        return vals

    AK_desc = build_ak(AK_REF, 6.0,  -DAK)
    AK_asc  = build_ak(AK_REF, 9.00, +DAK)

    print(f"\nDescending: {len(AK_desc)} kroków")
    print(f"Ascending:  {len(AK_asc)} kroków")
    sys.stdout.flush()

    def run_scan(AK_list, k1s, psi1s, k2s, psi2s, label):
        res = []
        k1c, psi1c = k1s, psi1s
        k2c, psi2c = k2s, psi2s
        fails2 = 0

        print(f"\n{'='*68}")
        print(f"Skan {label}:")
        print(f"  {'α_K':>7}  {'K*₁':>11}  {'K*₂':>12}  {'ratio':>9}  "
              f"{'ψ₀^B':>7}  {'dt':>5}")
        print(f"  {'-'*68}")
        sys.stdout.flush()

        for ak in AK_list:
            t0_ = time.time()

            # Branch A
            K1n, psi1n, gU1n = find_K_star(ak, k1c, psi1c,
                                             dpsi=0.15, n_pts=10,
                                             psi_min=0.80, K_frac=0.22)
            if np.isfinite(K1n) and abs(gU1n) < 0.08:
                k1c, psi1c = K1n, psi1n
            else:
                K1n = np.nan

            # Branch B
            K2n = np.nan
            if fails2 < 6:
                K2n, psi2n, gU2n = find_K_star(ak, k2c, psi2c,
                                                 dpsi=0.38, n_pts=14,
                                                 psi_min=PSI_B_MIN, K_frac=0.16)
                if np.isfinite(K2n) and abs(gU2n) < 0.08:
                    k2c, psi2c = K2n, psi2n
                    fails2 = 0
                else:
                    K2n = np.nan
                    fails2 += 1

            ratio = K2n/K1n if (np.isfinite(K1n) and np.isfinite(K2n) and K1n>0) else np.nan
            dt = time.time()-t0_

            k1s_ = f"{K1n:.6f}" if np.isfinite(K1n) else "---"
            k2s_ = f"{K2n:.8f}" if np.isfinite(K2n) else "---"
            rs_  = f"{ratio:.4f}" if np.isfinite(ratio) else "---"
            ref  = " <REF" if abs(ak-AK_REF)<1e-3 else ""
            print(f"  {ak:>7.4f}  {k1s_:>11}  {k2s_:>12}  {rs_:>9}  "
                  f"{psi2c:>7.4f}  {dt:>4.0f}s{ref}")
            sys.stdout.flush()

            res.append({'ak': ak, 'K1': K1n, 'K2': K2n,
                        'psi2': psi2c, 'ratio': ratio})
        return res

    t_total = time.time()
    rd = run_scan(AK_desc, K1_s, PSI1_s, K2_s, PSI2_s, "MALEJĄCY")
    ra = run_scan(AK_asc,  K1_s, PSI1_s, K2_s, PSI2_s, "ROSNĄCY")
    print(f"\nŁączny czas: {time.time()-t_total:.0f}s")
    sys.stdout.flush()

    # ── Analiza ───────────────────────────────────────────────────────────────
    seen = {}
    for r in rd + ra:
        ak = round(r['ak'], 8)
        if ak not in seen:
            seen[ak] = r
        elif np.isfinite(r['ratio']) and not np.isfinite(seen[ak]['ratio']):
            seen[ak] = r

    valid = [(r['ak'], r['K1'], r['K2'], r['ratio'])
             for r in sorted(seen.values(), key=lambda x: x['ak'])
             if np.isfinite(r['ratio'])]

    print("\n" + "="*72)
    print("TABELA:")
    print(f"  {'α_K':>7}  {'K*₁':>11}  {'K*₂':>12}  {'ratio':>10}")
    for ak, K1, K2, rat in valid:
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        print(f"  {ak:>7.4f}  {K1:>11.6f}  {K2:>12.8f}  {rat:>10.4f}{ref}")

    if len(valid) >= 3:
        aks  = np.array([r[0] for r in valid])
        rats = np.array([r[3] for r in valid])
        i_mx = np.argmax(rats)
        print(f"\n  argmax α_K = {aks[i_mx]:.4f}  (ratio = {rats[i_mx]:.4f})")
        idx_r = np.argmin(np.abs(aks - AK_REF))
        print(f"  α_K_ref    = {AK_REF:.4f}  (ratio = {rats[idx_r]:.4f})")
        delta = aks[i_mx] - AK_REF
        print(f"  Δ = {delta:+.4f}")
        if   abs(delta) < 0.1:   print("  *** HIPOTEZA POTWIERDZONA ***")
        elif abs(delta) < 0.3:   print(f"  *** CZĘŚCIOWE: |Δ|={abs(delta):.3f} ***")
        else:                     print(f"  *** OBALONA ***")
    else:
        print(f"\n  Za mało punktów ({len(valid)}) — Branch B mógł nie istnieć w tym zakresie α_K")
    print("="*72)
    print("KONIEC p108_v8b")
