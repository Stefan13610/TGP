#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v8c_extend.py  --  TGP v1
=====================================
Ścieżka 5: rozszerzony skan Branch B po v8b

v8b ograniczenia:
  - Descending zatrzymał się na α_K=8.2116 (ψ₀^B=5.32 > PSI_B_MAX=5.5)
  - Ascending zatrzymał się na α_K=8.7116 (ψ₀^B=2.05 < PSI_B_MIN=1.80)

v8c: brak sztywnych ograniczeń PSI_B_MIN/MAX, adaptacyjne dpsi.

Skan descending: 8.2116 → 3.0 w Δα_K=0.10
Skan ascending:  8.7116 → 8.98 w Δα_K=0.02 (bliżej bifurkacji α_SN≈8.9964)

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
N_EVAL = 500
RTOL   = 1e-7
ATOL   = 1e-9

AK_REF    = 8.5616
DPSI_DK_B = 18.78

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


def find_psi0_near(K, ak, psi_c, dpsi, n_pts, psi_min=0.8):
    lo = max(psi_min, psi_c - dpsi)
    hi = psi_c + dpsi
    if lo >= hi or lo <= 0: return np.nan, np.nan
    pts = np.linspace(lo, hi, n_pts)
    fv  = [((solve_one(K, p, ak)[0] - 1.0) if np.isfinite(solve_one(K, p, ak)[0]) else np.nan)
           for p in pts]
    # Memoize
    memo = {}
    for p in pts:
        pe, E = solve_one(K, p, ak)
        memo[p] = (pe - 1.0 if np.isfinite(pe) else np.nan, E)
    fv = [memo[p][0] for p in pts]

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
                t = -fi/(fi1-fi)
                candidates.append(pts[i]+t*(pts[i+1]-pts[i]))
    if not candidates: return np.nan, np.nan
    best = min(candidates, key=lambda p: abs(p-psi_c))
    _, E = solve_one(K, best, ak)
    return best, E


def gU_at_K(K, ak, psi_c, dpsi, n_pts, psi_min=0.8):
    p0, E = find_psi0_near(K, ak, psi_c, dpsi, n_pts, psi_min=psi_min)
    if not (np.isfinite(p0) and np.isfinite(E) and K > 0): return np.nan, psi_c
    return E/(4*np.pi*K) - 1.0, p0


def find_K_star(ak, K_seed, psi_seed, dpsi, n_pts, psi_min, K_frac):
    """Brentq w K — czysta ekstrapolacja z (K_seed, psi_seed)."""
    K_lo = K_seed * (1 - K_frac)
    K_hi = K_seed * (1 + K_frac)

    def pc(K): return psi_seed + DPSI_DK_B*(K - K_seed)

    g_lo, _ = gU_at_K(K_lo, ak, pc(K_lo), dpsi, n_pts, psi_min)
    g_hi, _ = gU_at_K(K_hi, ak, pc(K_hi), dpsi, n_pts, psi_min)

    if not (np.isfinite(g_lo) and np.isfinite(g_hi) and g_lo*g_hi < 0):
        for fac in [2.0, 3.0]:
            Kl = K_seed*(1-fac*K_frac); Kh = K_seed*(1+fac*K_frac)
            d2 = dpsi*(1+0.3*(fac-1))
            gl, _ = gU_at_K(Kl, ak, pc(Kl), d2, n_pts, psi_min)
            gh, _ = gU_at_K(Kh, ak, pc(Kh), d2, n_pts, psi_min)
            if np.isfinite(gl) and np.isfinite(gh) and gl*gh < 0:
                K_lo, K_hi, g_lo, dpsi = Kl, Kh, gl, d2
                break
        else:
            return np.nan, psi_seed, np.nan

    def fn(K):
        g, _ = gU_at_K(K, ak, pc(K), dpsi*0.85, n_pts, psi_min)
        return g if np.isfinite(g) else g_lo

    try:
        K_star = brentq(fn, K_lo, K_hi, xtol=1e-8, maxiter=60)
    except: return np.nan, psi_seed, np.nan

    g_f, p_f = gU_at_K(K_star, ak, pc(K_star), dpsi*0.75, n_pts, psi_min)
    if not np.isfinite(g_f) or abs(g_f) > 0.12: return np.nan, psi_seed, np.nan
    return K_star, (p_f if np.isfinite(p_f) else pc(K_star)), g_f


# ── Wyniki v8b (punkty referencyjne) ─────────────────────────────────────────
V8B = [
    (8.2116, 0.010632, 0.24616948, 23.1539, 5.3245),
    (8.2616, 0.010576, 0.22141456, 20.9347, 4.8819),
    (8.3116, 0.010522, 0.19788817, 18.8065, 4.4656),
    (8.3616, 0.010470, 0.17564842, 16.7769, 4.0753),
    (8.4116, 0.010419, 0.15474525, 14.8526, 3.7108),
    (8.4616, 0.010370, 0.13521872, 13.0397, 3.3720),
    (8.5116, 0.010323, 0.11709856, 11.3434, 3.0586),
    (8.5616, 0.010279, 0.10040394,  9.7676, 2.7701),
    (8.6116, 0.010239, 0.08513977,  8.3152, 2.5062),
    (8.6616, 0.010203, 0.07130227,  6.9880, 2.2663),
    (8.7116, 0.010174, 0.05887320,  5.7866, 2.0501),
]

# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("="*72)
    print("p108_v8c: Rozszerzony skan — bez ograniczeń PSI_B_MIN/MAX")
    print(f"  RTOL={RTOL:.0e}, N_EVAL={N_EVAL}")
    print("="*72)
    sys.stdout.flush()

    t0 = time.time()
    solve_one(0.246, 5.32, 8.2116)
    print(f"\n[TIMING] 1 ODE: {time.time()-t0:.3f}s")
    sys.stdout.flush()

    def build_ak(start, stop, step):
        vals = [start]
        v = start
        while True:
            v = round(v + step, 8)
            if step>0 and v > stop+1e-9: break
            if step<0 and v < stop-1e-9: break
            vals.append(v)
        return vals

    # ── DESCENDING ─────────────────────────────────────────────────────────────
    # Startujemy od najniższego dobrego punktu v8b
    D_AK0   = 8.2116
    D_K2_0  = 0.24616948
    D_PSI_0 = 5.3245

    AK_desc = build_ak(D_AK0, 3.0, -0.10)   # Δα_K=0.10 (szybszy)
    AK_asc  = build_ak(8.7116, 8.99, +0.02)  # Δα_K=0.02 (gęsty, blisko α_SN)

    print(f"\nDescending: {len(AK_desc)} kroków (Δα_K=0.10, od {D_AK0} do 3.0)")
    print(f"Ascending:  {len(AK_asc)} kroków (Δα_K=0.02, od 8.7116 do 8.99)")
    sys.stdout.flush()

    def run_scan(AK_list, K2_s, psi2_s, K1_s, psi1_s, label):
        res = []
        k2c, psi2c = K2_s, psi2_s
        k1c, psi1c = K1_s, psi1_s
        fails2 = 0

        print(f"\n{'='*68}")
        print(f"Skan {label}:")
        print(f"  {'α_K':>7}  {'K*₁':>11}  {'K*₂':>12}  {'ratio':>9}  "
              f"{'ψ₀^B':>8}  {'dpsi':>6}  {'dt':>5}")
        print(f"  {'-'*72}")
        sys.stdout.flush()

        for ak in AK_list:
            t0_ = time.time()

            # Branch A
            K1n, psi1n, gU1n = find_K_star(ak, k1c, psi1c,
                                             dpsi=0.15, n_pts=10,
                                             psi_min=0.80, K_frac=0.25)
            if np.isfinite(K1n) and abs(gU1n) < 0.12:
                k1c, psi1c = K1n, psi1n
            else:
                K1n = np.nan

            # Branch B — bez PSI_B_MAX, adaptacyjne dpsi
            K2n = np.nan
            dpsi_b = 0.45
            if fails2 < 7:
                # Adaptacyjne dpsi: większe dla dużego ψ₀^B
                dpsi_b = max(0.45, 0.08*(psi2c - 1.5) + 0.38)
                dpsi_b = min(dpsi_b, 4.0)

                K2n, psi2n, gU2n = find_K_star(
                    ak, k2c, psi2c,
                    dpsi=dpsi_b, n_pts=14, psi_min=1.05, K_frac=0.20)
                if np.isfinite(K2n) and abs(gU2n) < 0.12:
                    k2c, psi2c = K2n, psi2n
                    fails2 = 0
                else:
                    K2n = np.nan
                    fails2 += 1

            ratio = (K2n/K1n if (np.isfinite(K1n) and np.isfinite(K2n) and K1n>0)
                     else np.nan)
            dt = time.time()-t0_

            k1s_ = f"{K1n:.6f}" if np.isfinite(K1n) else "---"
            k2s_ = f"{K2n:.8f}" if np.isfinite(K2n) else "---"
            rs_  = f"{ratio:.2f}" if np.isfinite(ratio) else "---"
            print(f"  {ak:>7.4f}  {k1s_:>11}  {k2s_:>12}  {rs_:>9}  "
                  f"{psi2c:>8.4f}  {dpsi_b:>6.2f}  {dt:>4.0f}s")
            sys.stdout.flush()

            res.append({'ak': ak, 'K1': K1n, 'K2': K2n,
                        'psi2': psi2c, 'ratio': ratio})
        return res

    # Seed dla Branch A przy D_AK0 i A_START_AK
    K1_at_D, psi1_at_D, g = find_K_star(D_AK0, 0.010632, 1.2418,
                                          dpsi=0.15, n_pts=10, psi_min=0.80, K_frac=0.25)
    if not (np.isfinite(K1_at_D) and abs(g) < 0.12): K1_at_D, psi1_at_D = 0.010632, 1.2418

    K1_at_A, psi1_at_A, g = find_K_star(8.7116, 0.010174, 1.2418,
                                          dpsi=0.15, n_pts=10, psi_min=0.80, K_frac=0.25)
    if not (np.isfinite(K1_at_A) and abs(g) < 0.12): K1_at_A, psi1_at_A = 0.010174, 1.2418

    t_total = time.time()

    rd = run_scan(AK_desc, D_K2_0, D_PSI_0, K1_at_D, psi1_at_D, "MALEJĄCY (od 8.2116 do 3.0)")
    ra = run_scan(AK_asc,  0.05887320, 2.0501, K1_at_A, psi1_at_A, "ROSNĄCY (od 8.7116 do 8.99)")

    print(f"\nŁączny czas: {time.time()-t_total:.0f}s")
    sys.stdout.flush()

    # ── Kompletna tabela ───────────────────────────────────────────────────────
    seen = {}
    for (ak, K1, K2, rat, psi2) in V8B:
        seen[round(ak,4)] = {'ak':ak, 'K1':K1, 'K2':K2, 'ratio':rat, 'psi2':psi2}
    for r in rd + ra:
        ak = round(r['ak'],4)
        if ak not in seen:
            seen[ak] = r
        elif np.isfinite(r.get('ratio',np.nan)) and not np.isfinite(seen[ak].get('ratio',np.nan)):
            seen[ak] = r

    valid = [(r['ak'], r['K1'], r['K2'], r['ratio'])
             for r in sorted(seen.values(), key=lambda x: x['ak'])
             if np.isfinite(r.get('ratio', np.nan))]

    print("\n" + "="*72)
    print("KOMPLETNA TABELA (v8b + v8c):")
    print(f"  {'α_K':>7}  {'K*₁':>11}  {'K*₂':>13}  {'ratio':>10}")
    for ak, K1, K2, rat in valid:
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        print(f"  {ak:>7.4f}  {K1:>11.6f}  {K2:>13.8f}  {rat:>10.4f}{ref}")

    if len(valid) >= 3:
        aks  = np.array([r[0] for r in valid])
        rats = np.array([r[3] for r in valid])
        i_mx = np.argmax(rats)
        print(f"\n  argmax α_K = {aks[i_mx]:.4f}  (ratio = {rats[i_mx]:.4f})")
        ref_idx = np.argmin(np.abs(aks - AK_REF))
        print(f"  α_K_ref    = {AK_REF:.4f}  (ratio = {rats[ref_idx]:.4f})")
        delta = aks[i_mx] - AK_REF
        print(f"  Δ = {delta:+.4f}")

        if np.max(rats) >= 206.77:
            for i in range(len(rats)-1):
                if (rats[i]-206.77)*(rats[i+1]-206.77) <= 0:
                    t = (206.77-rats[i])/(rats[i+1]-rats[i])
                    ak_i = aks[i]+t*(aks[i+1]-aks[i])
                    print(f"\n  *** ratio=206.77 osiągnięte przy α_K≈{ak_i:.4f} ***")
        else:
            print(f"\n  Maksymalne ratio={np.max(rats):.2f} < 206.77")

        if abs(delta) < 0.15:   print("  *** HIPOTEZA POTWIERDZONA (|Δ|<0.15) ***")
        elif abs(delta) < 0.5:  print(f"  *** CZĘŚCIOWE: |Δ|={abs(delta):.3f} ***")
        else:                    print(f"  *** OBALONA ***")

    print("="*72)
    print("KONIEC p108_v8c")
