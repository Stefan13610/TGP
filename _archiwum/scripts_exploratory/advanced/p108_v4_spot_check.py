#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p108_v4_spot_check.py  --  TGP v1
=====================================
Ścieżka 5 OP-3: argmax K*₂/K*₁ — spot check dla wybranych α_K

Diagnoza problemów p108_v1/v2/v3:
 - p107_v3 znalazł K*₂=0.09999 przy α_K=8.5616 dzięki WĄSKIEMU skanowi K
   w oknie [0.098, 0.104] (10 punktów, spacing~0.0007) + ψ₀ brentq w [2.65, 2.90]
 - Poprzednie p108 używały za rzadkiego K-scanu (spacing~0.006) i trafiały
   w FAŁSZYWE zera g_U przy K≈0.068

Podejście p108_v4: dla każdego α_K:
  1) SZEROKI prescan g_U(K) na 60 log-spaced punktach [0.005, 0.30]
     → identificuje PRZYBLIŻONE K*₁ i K*₂ przez wiele potencjalnych zer
  2) Dla każdego kandydata (K_lo, K_hi) z zerem: GĘSTY lokalny skan (30 pts)
     → brentq z xtol=1e-7 by uchwycić sharp sign change
  3) Weryfikacja: |g_U| < 0.02 po brentq

KLUCZ: ψ₀ z brentq osobno dla każdego K, zakres [1.0, 3.5] z 40 punktami.
Drogie ale wiarygodne.

Dane referencyjne (p107_v3):
  α_K=8.5616: K*₁=0.010280 (g_U=-7e-8), K*₂=0.09999 (g_U=1.6e-6)
  ratio = 9.727

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
GAM    = 1.0
V1     = GAM/3.0 - GAM/4.0
LAM    = 5.501357e-06
AGAM   = 0.040
R_MAX  = 80.0
N_EVAL = 1000
RTOL   = 1e-10
ATOL   = 1e-12

AK_REF  = 8.5616
K1_REF  = 0.010280;  PSI1_REF = 1.2419
K2_REF  = 0.09999348; PSI2_REF = 2.7624

# ── ODE ───────────────────────────────────────────────────────────────────────
def V_mod(phi):  return GAM/3.*phi**3 - GAM/4.*phi**4 + LAM/6.*(phi-1.)**6
def dV_mod(phi): return GAM*phi**2 - GAM*phi**3 + LAM*(phi-1.)**5

def make_rhs(ak):
    def rhs(r, y):
        phi, dphi = y
        phi  = max(phi, 1e-10)
        kfac = 1.0 + ak/phi
        return [dphi, dV_mod(phi)/kfac + ak*dphi**2/(2.*phi**2*kfac) - 2./r*dphi]
    return rhs

_rhs_cache = {}
def get_rhs(ak):
    key = round(ak, 8)
    if key not in _rhs_cache:
        _rhs_cache[key] = make_rhs(ak)
    return _rhs_cache[key]


def integrate_full(K, psi0, ak):
    """Zwraca (phi_end, E) lub (nan, nan)."""
    rhs   = get_rhs(ak)
    dphi0 = -K / AGAM**2
    r_ev  = AGAM*(R_MAX/AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(rhs, [AGAM, R_MAX], [psi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL, t_eval=r_ev)
        if sol.t[-1] < R_MAX*0.99:
            return np.nan, np.nan
        r    = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        kfac = 1.+ak/phi
        Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*kfac*r**2, r)
        Vv   = np.array([V_mod(float(q)) for q in phi])
        Ep   = 4*np.pi*np.trapezoid((Vv-V1)*r**2, r)
        return sol.y[0,-1], Ek+Ep
    except Exception:
        return np.nan, np.nan


def find_psi0_wide(K, ak, psi_lo=1.0, psi_hi=3.5, n_pts=40):
    """
    Scana ψ₀ w [psi_lo, psi_hi] z n_pts punktami.
    Zwraca WSZYSTKIE (psi, g_U) gdzie g_U zmienia znak.
    """
    psis = np.linspace(psi_lo, psi_hi, n_pts)
    phi_ends = []
    Es = []
    for p in psis:
        pe, E = integrate_full(K, p, ak)
        phi_ends.append(pe)
        Es.append(E)

    # Szukaj wszystkich zer phi_end - 1
    solutions = []
    for i in range(len(psis)-1):
        fi  = phi_ends[i]  - 1.0 if np.isfinite(phi_ends[i])  else np.nan
        fi1 = phi_ends[i+1]- 1.0 if np.isfinite(phi_ends[i+1]) else np.nan
        if not (np.isfinite(fi) and np.isfinite(fi1)): continue
        if fi*fi1 < 0:
            try:
                fn = lambda p: integrate_full(K, p, ak)[0] - 1.0
                p0 = brentq(fn, psis[i], psis[i+1], xtol=1e-9, maxiter=60)
                pe_f, E_f = integrate_full(K, p0, ak)
                gU = E_f/(4*np.pi*K) - 1.0 if np.isfinite(E_f) and K>0 else np.nan
                solutions.append((p0, pe_f, gU))
            except Exception:
                pass
    return solutions  # lista (psi0, phi_end, g_U)


def scan_gU_vs_K(ak, K_lo=0.003, K_hi=0.30, n_K=60, psi_lo=1.0, psi_hi=3.5, n_psi=40):
    """
    Pełny prescan: dla każdego K z n_K punktami (log-spaced):
      - Szuka wszystkich ψ₀ rozwiązań φ(R_MAX)=1 w [psi_lo, psi_hi]
      - Dla każdego rozwiązania oblicza g_U(K)
    Zwraca listę (K, psi0, g_U) dla wszystkich znalezionych rozwiązań.
    """
    K_arr = np.exp(np.linspace(np.log(K_lo), np.log(K_hi), n_K))
    all_sols = []
    for K in K_arr:
        sols = find_psi0_wide(K, ak, psi_lo=psi_lo, psi_hi=psi_hi, n_pts=n_psi)
        for (psi0, pe, gU) in sols:
            if np.isfinite(gU):
                all_sols.append((K, psi0, gU))
    return all_sols


def find_K_star_precise(ak, K_lo_hint, K_hi_hint, psi_hint, label=""):
    """
    Precyzyjne K*: dense scan w [K_lo_hint, K_hi_hint] + brentq.
    ψ₀ śledzone z wąskim oknem wokół psi_hint.
    """
    # Dense scan w oknie
    n_dense = 40
    K_dense = np.linspace(K_lo_hint, K_hi_hint, n_dense)
    gU_d    = np.full(n_dense, np.nan)
    psi_d   = np.full(n_dense, np.nan)
    pc = psi_hint

    for i, K in enumerate(K_dense):
        # Wąski skan ψ₀ wokół pc
        dpsi = 0.25
        psis = np.linspace(max(pc-dpsi, 0.8), pc+dpsi, 25)
        fv   = [integrate_full(K, p, ak)[0]-1.0 for p in psis]
        for j in range(len(fv)-1):
            if np.isfinite(fv[j]) and np.isfinite(fv[j+1]) and fv[j]*fv[j+1]<0:
                try:
                    fn = lambda p: integrate_full(K, p, ak)[0]-1.0
                    p0 = brentq(fn, psis[j], psis[j+1], xtol=1e-9, maxiter=50)
                    _, E = integrate_full(K, p0, ak)
                    gU_d[i] = E/(4*np.pi*K)-1.0
                    psi_d[i] = p0
                    pc = p0  # update
                    break
                except Exception:
                    pass

    # Szukaj zmiany znaku
    for i in range(len(gU_d)-1):
        gi, gi1 = gU_d[i], gU_d[i+1]
        if not (np.isfinite(gi) and np.isfinite(gi1)): continue
        if gi*gi1 < 0:
            Ki, Ki1 = K_dense[i], K_dense[i+1]
            pm = [psi_d[i] if np.isfinite(psi_d[i]) else psi_hint]
            def fn(K_t):
                psis2 = np.linspace(max(pm[0]-0.15, 0.8), pm[0]+0.15, 20)
                fv2   = [integrate_full(K_t, p, ak)[0]-1.0 for p in psis2]
                for jj in range(len(fv2)-1):
                    if np.isfinite(fv2[jj]) and np.isfinite(fv2[jj+1]) and fv2[jj]*fv2[jj+1]<0:
                        try:
                            fn2 = lambda p: integrate_full(K_t, p, ak)[0]-1.0
                            p0  = brentq(fn2, psis2[jj], psis2[jj+1], xtol=1e-9, maxiter=40)
                            pm[0] = p0
                            _, E  = integrate_full(K_t, p0, ak)
                            return E/(4*np.pi*K_t)-1.0
                        except Exception:
                            pass
                return gi
            try:
                K_star = brentq(fn, Ki, Ki1, xtol=1e-7, maxiter=40)
                psis3 = np.linspace(max(pm[0]-0.15, 0.8), pm[0]+0.15, 20)
                fv3   = [integrate_full(K_star, p, ak)[0]-1.0 for p in psis3]
                for jj in range(len(fv3)-1):
                    if np.isfinite(fv3[jj]) and np.isfinite(fv3[jj+1]) and fv3[jj]*fv3[jj+1]<0:
                        try:
                            fn3 = lambda p: integrate_full(K_star, p, ak)[0]-1.0
                            p0  = brentq(fn3, psis3[jj], psis3[jj+1], xtol=1e-9, maxiter=40)
                            _, E = integrate_full(K_star, p0, ak)
                            return K_star, p0, E/(4*np.pi*K_star)-1.0
                        except Exception:
                            pass
                return K_star, pm[0], np.nan
            except Exception:
                t = -gi/(gi1-gi)
                return K_dense[i]+t*(K_dense[i+1]-K_dense[i]), psi_d[i], gi+t*(gi1-gi)

    return np.nan, psi_hint, np.nan


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    print("="*72)
    print("p108_v4: spot-check K*₁, K*₂ dla wybranych α_K")
    print(f"N_eval={N_EVAL}, RTOL={RTOL:.0e}")
    print(f"Referencja: α_K={AK_REF}, K*₁={K1_REF:.6f}, K*₂={K2_REF:.8f}")
    print("="*72)

    # α_K do sprawdzenia: referencja + kilka poniżej i powyżej
    AK_CHECK = [6.0, 7.0, 7.5, 8.0, 8.2, 8.3, 8.4, 8.5, AK_REF,
                8.62, 8.68, 8.75, 8.80, 8.90, 8.95]

    results = []
    for ak in AK_CHECK:
        t0 = time.time()
        print(f"\n{'─'*60}")
        print(f"α_K = {ak:.4f}")

        # ── PRESCAN K ─────────────────────────────────────────────────────
        print(f"  Prescan K ∈ [0.003, 0.30], ψ₀ ∈ [1.0, 3.5], 60 K × 40 ψ₀...", end="", flush=True)
        t1 = time.time()
        all_sols = scan_gU_vs_K(ak, K_lo=0.003, K_hi=0.30,
                                  n_K=60, psi_lo=1.0, psi_hi=3.5, n_psi=40)
        print(f" {time.time()-t1:.0f}s, {len(all_sols)} rozwiązań", flush=True)

        if not all_sols:
            print("  BRAK rozwiązań w prescanie!")
            results.append({'ak': ak, 'K1': np.nan, 'K2': np.nan, 'ratio': np.nan})
            continue

        # ── Identyfikacja kandydatów K* ────────────────────────────────────
        # Szukaj zmian znaku g_U w posortowanej liście (K, gU)
        # Sortuj wg K, a dla każdego K weź NAJNIŻSZY ψ₀ (Branch A) i następny (Branch B)
        from collections import defaultdict
        K_sols = defaultdict(list)
        for (K, psi0, gU) in all_sols:
            K_sols[round(K, 8)].append((psi0, gU))

        # Stwórz osobne krzywe dla Branch A (niskie ψ₀) i Branch B (wyższe ψ₀)
        K_sorted = sorted(K_sols.keys())
        branchA = []  # (K, psi0, gU) dla najniższego ψ₀ przy każdym K
        branchB = []  # (K, psi0, gU) dla drugiego ψ₀ przy każdym K

        for K in K_sorted:
            slist = sorted(K_sols[K], key=lambda x: x[0])  # sort by psi0
            if slist:
                branchA.append((K, slist[0][0], slist[0][1]))
            if len(slist) >= 2:
                branchB.append((K, slist[1][0], slist[1][1]))
            elif len(slist) == 1 and slist[0][0] > 2.0:
                # Tylko jedno rozwiązanie ale wysoki ψ₀ → może Branch B
                branchB.append((K, slist[0][0], slist[0][1]))

        print(f"  Branch A: {len(branchA)} punktów K")
        print(f"  Branch B: {len(branchB)} punktów K")

        # Debug: pokaż gU dla Branch B
        if branchB:
            Kb = np.array([x[0] for x in branchB])
            gUb = np.array([x[2] for x in branchB])
            psi_b = np.array([x[1] for x in branchB])
            print(f"  Branch B gU range: [{np.nanmin(gUb):.3f}, {np.nanmax(gUb):.3f}]")
            print(f"  Branch B K range:  [{np.nanmin(Kb):.5f}, {np.nanmax(Kb):.5f}]")
            # Pokaż kilka punktów Branch B
            step = max(1, len(branchB)//6)
            for i in range(0, len(branchB), step):
                print(f"    K={Kb[i]:.5f}  ψ₀={psi_b[i]:.4f}  g_U={gUb[i]:.4f}")

        # ── Znajdź K* przez zmianę znaku g_U ─────────────────────────────
        def find_zero_in_branch(branch, label=""):
            if len(branch) < 2:
                return np.nan, np.nan, np.nan
            K_arr  = np.array([x[0] for x in branch])
            gU_arr = np.array([x[2] for x in branch])
            psi_arr= np.array([x[1] for x in branch])

            for i in range(len(gU_arr)-1):
                gi, gi1 = gU_arr[i], gU_arr[i+1]
                if not (np.isfinite(gi) and np.isfinite(gi1)): continue
                if gi*gi1 < 0:
                    Ki, Ki1  = K_arr[i], K_arr[i+1]
                    psi_hint = psi_arr[i]
                    print(f"  {label}: zmiana znaku g_U w [{Ki:.5f}, {Ki1:.5f}], "
                          f"gU=[{gi:.4f},{gi1:.4f}], ψ₀≈{psi_hint:.4f}")
                    K_prec, psi_prec, gU_prec = find_K_star_precise(
                        ak, Ki*0.95, Ki1*1.05, psi_hint, label=label)
                    return K_prec, psi_prec, gU_prec
            return np.nan, np.nan, np.nan

        # Branch A
        K1_found, psi1_found, gU1_found = find_zero_in_branch(branchA, "Branch A")
        # Branch B
        K2_found, psi2_found, gU2_found = find_zero_in_branch(branchB, "Branch B")

        dt = time.time() - t0

        ratio = K2_found/K1_found if (np.isfinite(K1_found) and np.isfinite(K2_found) and K1_found>0) else np.nan
        print(f"  WYNIK: K*₁={K1_found:.6f}(g={gU1_found:.2e})  "
              f"K*₂={K2_found:.6f}(g={gU2_found:.2e})  ratio={ratio:.4f}  [{dt:.0f}s]")

        results.append({
            'ak': ak, 'K1': K1_found, 'K2': K2_found,
            'psi1': psi1_found, 'psi2': psi2_found,
            'gU1': gU1_found, 'gU2': gU2_found, 'ratio': ratio
        })

    # ── Tabela końcowa ────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("TABELA KOŃCOWA:")
    print(f"  {'α_K':>8}  {'K*₁':>12}  {'K*₂':>12}  {'ratio':>10}  {'g_U^A':>9}  {'g_U^B':>9}")
    print("  " + "-"*70)
    valid = []
    for r in results:
        ak = r['ak']
        ref = " <REF" if abs(ak-AK_REF)<1e-3 else ""
        k1s = f"{r['K1']:.6f}" if np.isfinite(r['K1']) else "---"
        k2s = f"{r['K2']:.6f}" if np.isfinite(r['K2']) else "---"
        g1s = f"{r['gU1']:.2e}" if np.isfinite(r.get('gU1', np.nan)) else "---"
        g2s = f"{r['gU2']:.2e}" if np.isfinite(r.get('gU2', np.nan)) else "---"
        rs  = f"{r['ratio']:.4f}" if np.isfinite(r['ratio']) else "---"
        print(f"  {ak:>8.4f}  {k1s:>12}  {k2s:>12}  {rs:>10}  {g1s:>9}  {g2s:>9}{ref}")
        if np.isfinite(r['ratio']):
            valid.append(r)

    # ── Analiza ───────────────────────────────────────────────────────────────
    print("\n" + "="*72)
    print("ANALIZA MAKSIMUM:")
    if len(valid) >= 3:
        aks  = np.array([r['ak']    for r in valid])
        rats = np.array([r['ratio'] for r in valid])
        i_mx = np.argmax(rats)
        print(f"  argmax α_K = {aks[i_mx]:.4f}  (ratio = {rats[i_mx]:.4f})")
        print(f"  α_K_ref    = {AK_REF:.4f}  (ratio = {rats[np.argmin(np.abs(aks-AK_REF))]:.4f})")
        delta = aks[i_mx] - AK_REF
        print(f"  Δ = {delta:.4f}")

        print("\n  Top-5:")
        for i in np.argsort(rats)[::-1][:5]:
            ref = " <REF" if abs(valid[i]['ak']-AK_REF)<1e-3 else ""
            print(f"    α_K={valid[i]['ak']:.4f}  ratio={rats[i]:.4f}{ref}")

        if abs(delta) < 0.05:
            print("\n  *** HIPOTEZA POTWIERDZONA ***")
        elif abs(delta) < 0.20:
            print(f"\n  *** CZĘŚCIOWE: |delta|={abs(delta):.3f} ***")
        else:
            print(f"\n  *** HIPOTEZA OBALONA: argmax≠α_K_ref ***")
