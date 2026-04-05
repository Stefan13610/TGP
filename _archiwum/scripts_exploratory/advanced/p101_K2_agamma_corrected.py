#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p101_K2_agamma_corrected.py  --  TGP v1
=========================================
Cel: OP-3 Sciezka 2 (korekta p100) — czy K*_2 fizyczne (psi_0~5) istnieje
     dla a_Gam < 0.040? Waska metoda psi zamiast "ostatniego zera".

Diagnoza p100: "ostatnie zero psi" w [1.001,8.0] wybieralo niefizyczne galezi.
Korekta: osobne okna psi dla K*_1 i K*_2.

Metoda:
  K*_1 branch:  psi scan w [1.001, 1.8]  (fizyczny K*_1 ma psi_0~1.24)
  K*_2 branch:  psi scan w [3.5, 7.0]    (fizyczny K*_2 ma psi_0~5.0 z p99)

  Dla K*_1: skan K w [K1_lo, K1_hi] (proporcjonalny do a_Gam)
  Dla K*_2: skan K w [K2_lo, K2_hi] (szerszy zakres)

  a_Gam w {0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.050, 0.060}
  r_max = 80 (jak p99)

Testy:
  P1: K*_1 znalezione przy wszystkich a_Gam (sanity check)
  P2: K*_1 ~ K*_1_ref * (a_Gam/0.040) (skalowanie Schwarzschilda, tol 8%)
  P3: K*_2 znalezione przy a_Gam = 0.040 (walidacja)
  P4: K*_2/K*_1 przy a_Gam = 0.040 zblizony do p99 (~21-25)
  P5: K*_2 NIE istnieje przy a_Gam <= a_c_estim dla jakiegos a_c

Data: 2026-03-26
"""

import sys, io, warnings, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# ─── parametry fizyczne ──────────────────────────────────────────────────────
ALPHA_K   = 8.5616
LAM       = 5.501357e-06
GAM       = 1.0
V1        = GAM/3 - GAM/4
K1_INF    = 0.010029        # K*_1(inf) przy a_Gam=0.040 (p97)
K2_P99    = 0.219088        # K*_2 przy a_Gam=0.040, r_max=80 (p99)
AGAM_PHYS = 0.040
R_MAX     = 80.0

# ─── skan a_Gam ──────────────────────────────────────────────────────────────
AGAM_VALUES = [0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.050, 0.060]

# ─── parametry numeryczne ─────────────────────────────────────────────────────
# Osobne okna psi dla K*_1 i K*_2 !
PSI1_LO = 1.001;  PSI1_HI = 1.80   # dla K*_1 (psi_0 ~ 1.24)
PSI2_LO = 3.50;   PSI2_HI = 7.00   # dla K*_2 (psi_0 ~ 5.0)
N_PSI   = 60
N_K1    = 30
N_K2    = 45
N_EVAL  = 700
RTOL    = 1e-8
ATOL    = 1e-10

# ─── ODE ─────────────────────────────────────────────────────────────────────
def dV_mod(phi):
    return GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1)**5

def V_mod(phi):
    return GAM/3*phi**3 - GAM/4*phi**4 + LAM/6*(phi - 1)**6

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K / phi
    return [dphi, dV_mod(phi)/kfac + ALPHA_K*dphi**2/(2*phi**2*kfac) - 2.0/r*dphi]

# ─── strzelanie ──────────────────────────────────────────────────────────────
def phi_at_rmax(psi, K, a_Gam, r_max=R_MAX):
    phi0  = psi
    dphi0 = -K / a_Gam**2
    r_ev  = a_Gam * (r_max / a_Gam)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(ode_rhs, [a_Gam, r_max], [phi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] < r_max * 0.99:
            return np.nan
        return float(sol.y[0, -1])
    except Exception:
        return np.nan

def energy_soliton(psi, K, a_Gam, r_max=R_MAX):
    phi0  = psi
    dphi0 = -K / a_Gam**2
    r_ev  = a_Gam * (r_max / a_Gam)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(ode_rhs, [a_Gam, r_max], [phi0, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] < r_max * 0.99:
            return np.nan
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        if not np.all(np.isfinite(phi)) or not np.all(np.isfinite(dphi)):
            return np.nan
        kfac = 1.0 + ALPHA_K / phi
        Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2 * kfac * r**2, r)
        Ep   = 4*np.pi * np.trapezoid((V_mod(phi) - V1) * r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan

def g_at_psi(psi, K, a_Gam):
    E = energy_soliton(psi, K, a_Gam)
    return E / (4*np.pi*K) - 1.0 if np.isfinite(E) else np.nan

# ─── znajdz zero psi w zadanym oknie ─────────────────────────────────────────
def find_psi_zero_in_window(K, a_Gam, psi_lo, psi_hi, n_psi=N_PSI):
    """Znajdz PIERWSZE zero phi(r_max)=1 w [psi_lo, psi_hi].
    Zwraca psi_zero lub nan."""
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv   = [phi_at_rmax(p, K, a_Gam) - 1.0 for p in psis]
    for i in range(len(Fv) - 1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not np.isfinite(fi) or not np.isfinite(fi1):
            continue
        if fi * fi1 < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p, K, a_Gam) - 1.0,
                            psis[i], psis[i+1], xtol=1e-7, maxiter=50)
                return pz
            except Exception:
                pass
    return np.nan

def g_U_narrow(K, a_Gam, psi_lo, psi_hi):
    """g_U uzywajac pierwszego zera psi w [psi_lo, psi_hi]."""
    pz = find_psi_zero_in_window(K, a_Gam, psi_lo, psi_hi)
    if np.isnan(pz):
        return np.nan, np.nan
    return g_at_psi(pz, K, a_Gam), pz

# ─── skan K i znajdz K* ───────────────────────────────────────────────────────
def find_Kstar_narrow(a_Gam, K_lo, K_hi, psi_lo, psi_hi, n_K, verbose=False):
    """Znajdz K* w [K_lo, K_hi] uzywajac okna psi [psi_lo, psi_hi]."""
    Ks  = np.linspace(K_lo, K_hi, n_K)
    gv  = []
    psv = []
    for K in Ks:
        g, pz = g_U_narrow(K, a_Gam, psi_lo, psi_hi)
        gv.append(g)
        psv.append(pz)

    # Szukaj zmian znaku
    sign_changes = []
    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if np.isfinite(gi) and np.isfinite(gi1) and gi * gi1 < 0:
            sign_changes.append((i, Ks[i], Ks[i+1], gv[i], gv[i+1], psv[i]))

    if verbose and sign_changes:
        for idx, Kl, Kr, gl, gr, pz in sign_changes:
            print(f"      zmiana znaku K=[{Kl:.5f},{Kr:.5f}] "
                  f"g=[{gl:.4f},{gr:.4f}] psi~{pz:.4f}")

    if not sign_changes:
        # Sprawdz czy g_U jest zawsze nan (psi zero nie istnieje w oknie)
        n_nan = sum(np.isnan(g) for g in gv)
        if verbose:
            n_valid = len([g for g in gv if np.isfinite(g)])
            print(f"      Brak zmian znaku: {n_valid}/{len(gv)} valid g_U")
        return np.nan, np.nan

    # Bierz pierwsze zero
    idx, Kl, Kr, gl, gr, _ = sign_changes[0]
    try:
        Kstar = brentq(lambda K: g_U_narrow(K, a_Gam, psi_lo, psi_hi)[0],
                       Kl, Kr, xtol=1e-5, maxiter=40)
    except Exception:
        Kstar = Kl - gl * (Kr - Kl) / (gr - gl)

    _, psi_at_Kstar = g_U_narrow(Kstar, a_Gam, psi_lo, psi_hi)
    return Kstar, psi_at_Kstar

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 68)
    print("p101_K2_agamma_corrected.py  --  TGP v1")
    print("OP-3 Sciezka 2 (korekta p100): waska metoda psi")
    print("=" * 68)
    print(f"alpha_K={ALPHA_K}, r_max={R_MAX}")
    print(f"K*_1 psi okno: [{PSI1_LO},{PSI1_HI}]")
    print(f"K*_2 psi okno: [{PSI2_LO},{PSI2_HI}]")
    print(f"a_Gam values: {AGAM_VALUES}")
    print()

    results = {}

    for a_Gam in AGAM_VALUES:
        t0 = time.time()
        print(f"─── a_Gam = {a_Gam:.3f} ───────────────────────────────────")
        scale = a_Gam / AGAM_PHYS

        # K*_1 skan: [0.4 * K1_ref, 2.5 * K1_ref]
        K1_ref = K1_INF * scale
        K1_lo  = K1_ref * 0.40
        K1_hi  = K1_ref * 2.50
        print(f"  K*_1 skan: [{K1_lo:.5f},{K1_hi:.5f}], psi=[{PSI1_LO},{PSI1_HI}]")
        K1, psi1 = find_Kstar_narrow(a_Gam, K1_lo, K1_hi,
                                     PSI1_LO, PSI1_HI, N_K1, verbose=True)

        # K*_2 skan: szerszy zakres, zakotw. do K2_P99 * scale (jesli K*_2 ~ K*_1)
        # Uzywamy tez szerokiego zakresu na wypadek innego skalowania
        K2_lo = max(K1_hi * 2.0, K2_P99 * scale * 0.10)
        K2_hi = K2_P99 * scale * 5.0
        print(f"  K*_2 skan: [{K2_lo:.4f},{K2_hi:.4f}], psi=[{PSI2_LO},{PSI2_HI}]")
        K2, psi2 = find_Kstar_narrow(a_Gam, K2_lo, K2_hi,
                                     PSI2_LO, PSI2_HI, N_K2, verbose=True)

        dt = time.time() - t0
        results[a_Gam] = {'K1': K1, 'K2': K2, 'psi1': psi1, 'psi2': psi2}
        K1s = f"{K1:.6f}" if np.isfinite(K1) else "BRAK"
        K2s = f"{K2:.6f}" if np.isfinite(K2) else "BRAK"
        p1s = f"{psi1:.4f}" if np.isfinite(psi1) else "nan"
        p2s = f"{psi2:.4f}" if np.isfinite(psi2) else "nan"
        print(f"  => K*_1={K1s}  (psi={p1s})")
        print(f"     K*_2={K2s}  (psi={p2s})")
        print(f"  Czas: {dt:.0f}s")
        print()

    # ── tabela ──────────────────────────────────────────────────────────────
    print("=" * 68)
    print("TABELA: K*_1, K*_2, K*_2/K*_1 vs a_Gam")
    print("=" * 68)
    print(f"{'a_Gam':>8} {'K*_1':>12} {'K*_1/ref':>10} {'K*_2':>12} "
          f"{'K*_2/K*_1':>12} {'K*_2?':>8}")
    print("-" * 68)
    for a_Gam in AGAM_VALUES:
        res   = results[a_Gam]
        K1    = res['K1'];  K2 = res['K2']
        scale = a_Gam / AGAM_PHYS
        ref1  = K1_INF * scale
        r1    = K1/ref1 if np.isfinite(K1) and ref1>0 else np.nan
        r21   = K2/K1   if np.isfinite(K1) and np.isfinite(K2) and K1>0 else np.nan
        K1s   = f"{K1:.6f}" if np.isfinite(K1) else "  BRAK  "
        K2s   = f"{K2:.6f}" if np.isfinite(K2) else "  BRAK  "
        r1s   = f"{r1:.4f}" if np.isfinite(r1) else " BRAK"
        r21s  = f"{r21:.3f}" if np.isfinite(r21) else " BRAK"
        K2ex  = "TAK" if np.isfinite(K2) else "NIE"
        print(f"{a_Gam:>8.3f} {K1s:>12} {r1s:>10} {K2s:>12} {r21s:>12} {K2ex:>8}")

    # ── wyznaczenie a_c ─────────────────────────────────────────────────────
    print()
    K2_ex = {a: np.isfinite(results[a]['K2']) for a in AGAM_VALUES}
    sorted_a = sorted(AGAM_VALUES)
    a_c = None
    a_c_prev = None
    for i, a in enumerate(sorted_a):
        if K2_ex[a]:
            a_c = a
            a_c_prev = sorted_a[i-1] if i > 0 else None
            break

    # ── PASS/FAIL ───────────────────────────────────────────────────────────
    print("=" * 68)
    print("TESTY PASS/FAIL")
    print("=" * 68)
    passes = []

    # P1: K*_1 przy wszystkich a_Gam
    n1 = sum(np.isfinite(results[a]['K1']) for a in AGAM_VALUES)
    p1 = n1 == len(AGAM_VALUES)
    passes.append(p1)
    print(f"P1 [K*_1 przy wszystkich]:       {'PASS' if p1 else 'FAIL'}  ({n1}/{len(AGAM_VALUES)})")

    # P2: K*_1 ~ ref (Schwarzschild, tol 8%)
    devs = []
    for a in AGAM_VALUES:
        K1 = results[a]['K1']
        ref = K1_INF * (a / AGAM_PHYS)
        if np.isfinite(K1) and ref > 0:
            devs.append(abs(K1/ref - 1))
    if devs:
        mx = max(devs)
        p2 = mx < 0.08
        print(f"P2 [K*_1 skalowanie Schwarzschil.]: {'PASS' if p2 else 'FAIL'}  "
              f"(max odch={mx*100:.1f}%, prog=8%)")
    else:
        p2 = False
        print(f"P2 [K*_1 skalowanie]:            FAIL  (brak danych)")
    passes.append(p2)

    # P3: K*_2 przy a_Gam=0.040
    p3 = np.isfinite(results[AGAM_PHYS]['K2'])
    passes.append(p3)
    K2p = results[AGAM_PHYS]['K2']
    print(f"P3 [K*_2 przy a_Gam=0.040]:      {'PASS' if p3 else 'FAIL'}  "
          f"(K*_2={'%.5f' % K2p if p3 else 'BRAK'})")

    # P4: K*_2/K*_1 przy 0.040 ~ 21-25 (z p99)
    if np.isfinite(results[AGAM_PHYS]['K1']) and np.isfinite(results[AGAM_PHYS]['K2']):
        r21_phys = results[AGAM_PHYS]['K2'] / results[AGAM_PHYS]['K1']
        p4 = 15.0 < r21_phys < 35.0
        print(f"P4 [K*_2/K*_1 ~ 21-25 (p99)]:   {'PASS' if p4 else 'FAIL'}  "
              f"(ratio={r21_phys:.3f}, prog=[15,35])")
    else:
        p4 = False
        print(f"P4 [K*_2/K*_1 ~ 21-25]:          FAIL  (brak danych)")
    passes.append(p4)

    # P5: K*_2 NIE istnieje dla jakiegos a_Gam < 0.040
    a_below = [a for a in AGAM_VALUES if a < AGAM_PHYS]
    n_no_K2 = sum(not K2_ex[a] for a in a_below)
    p5 = n_no_K2 > 0
    passes.append(p5)
    print(f"P5 [K*_2 znika dla mniejszych]:  {'PASS' if p5 else 'FAIL'}  "
          f"(brak K*_2 przy {n_no_K2}/{len(a_below)} pkt < 0.040)")

    print(f"\nWYNIK: {sum(passes)}/{len(passes)} PASS")

    # ── podsumowanie ─────────────────────────────────────────────────────────
    print("\n" + "=" * 68)
    print("PODSUMOWANIE FIZYCZNE  --  OP-3 Sciezka 2")
    print("=" * 68)
    if a_c is not None and a_c_prev is not None:
        print(f"  K*_2 pojawia sie przy a_Gam >= a_c, a_c w ({a_c_prev:.3f},{a_c:.3f}]")
        if abs(a_c - AGAM_PHYS) < 0.006:
            print(f"  ** a_c ~ a_Gam_fizyczne = 0.040  =>  OP-3 Sciezka 2 POTWIERDZONA!")
        else:
            print(f"  ** a_c != 0.040  =>  OP-3 Sciezka 2 NIE potwierdzona")
    elif a_c is None:
        print(f"  K*_2 NIE istnieje przy zadnym a_Gam w [{min(AGAM_VALUES)},{max(AGAM_VALUES)}]!")
        print(f"  OP-3 Sciezka 2 nie dotyczy bifurkacji K*_2 (brak 2. solitonu?)")
    else:
        print(f"  K*_2 istnieje juz przy a_Gam={a_c:.3f} <= min (wszedzie)")
        print(f"  OP-3 Sciezka 2 NIE potwierdzona (K*_2 istnieje dla wszystkich a_Gam)")
    print("=" * 68)
