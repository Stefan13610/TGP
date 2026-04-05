#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p100_K2_vs_agamma.py  --  TGP v1
====================================
Cel: OP-3 Sciezka 2 (relacyjna) — czy K*_2 istnieje dla a_Gam < 0.040?

Kontekst:
  p99: K*_2 istnieje przy a_Gam=0.040 (K*_2/K*_1(inf)~21-25)
  Hipoteza Sciezki 2: a_c(alpha_K) = a_Gam = 0.040 (krytyczna wart. substratu)
  Mechanizm (skalowanie): ODE w wsp. bezwym. rho=r/a_Gam zawiera czlon
    potencjalu ~a_Gam^2 * V_mod -- dla malego a_Gam V_mod zanikajacy;
    K*_2 (o duzym psi_0~5) moze nie istniec ponizej progu a_c.

Pytanie: Dla jakich a_Gam istnieje K*_2?

Metoda:
  1. Ustalony alpha_K = 8.5616, r_max = 80 (kompromis dokladnosc/czas)
  2. Skan a_Gam w {0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.050, 0.060}
  3. Dla kazdego a_Gam:
     a. K*_1 (zakres K proporcjonalny do a_Gam, ref: K*_1(inf)=0.010029 przy a=0.040)
     b. K*_2 (zakres K proporcjonalny do a_Gam, ref: K*_2~0.24 przy a=0.040)
     c. czy K*_2 istnieje? (zmiana znaku g_U)
  4. Wyznaczenie a_c: najmniejsze a_Gam przy ktorym K*_2 istnieje
  5. Sprawdzenie: czy a_c ~ 0.040 (= a_Gam_fizyczne)?

Testy:
  P1: K*_1 znalezione przy wszystkich a_Gam (sanity check)
  P2: K*_1 ~ K*_1_ref * (a_Gam/0.040) (skalowanie Schwarzschilda)
  P3: K*_2 istnieje przy a_Gam = 0.040 (powt. p99)
  P4: K*_2 NIE istnieje dla a_Gam <= a_c < 0.040 (test hipotezy)
  P5: a_c(alpha_K) ~ 0.040 (Sciezka 2 spelniona)

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
V1        = GAM/3 - GAM/4   # = 1/12
K1_INF    = 0.010029        # K*_1(inf) przy a_Gam=0.040 (p97)
K2_REF    = 0.237           # K*_2 szacunek przy a_Gam=0.040 (p99, r=120)
AGAM_PHYS = 0.040
R_MAX     = 80.0            # kompromis dokladnosc/czas

# ─── skan a_Gam ──────────────────────────────────────────────────────────────
AGAM_VALUES = [0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.050, 0.060]

# ─── parametry numeryczne ─────────────────────────────────────────────────────
N_K    = 40     # punktow w skanie K
N_PSI  = 70     # punktow w skanie psi
PSI_HI = 8.0    # rozszerzone (p99 pokazalo psi~5; probujemy do 8)
PSI_LO = 1.001
N_EVAL = 600
RTOL   = 1e-8
ATOL   = 1e-10

# ─── ODE ─────────────────────────────────────────────────────────────────────
def dV_mod(phi):
    return GAM*phi**2 - GAM*phi**3 + LAM*(phi - 1)**5

def V_mod(phi):
    return GAM/3*phi**3 - GAM/4*phi**4 + LAM/6*(phi - 1)**6

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K / phi
    d2phi = dV_mod(phi)/kfac + ALPHA_K*dphi**2/(2*phi**2*kfac) - 2.0/r*dphi
    return [dphi, d2phi]

# ─── strzelanie ──────────────────────────────────────────────────────────────
def phi_at_rmax(psi, K, a_Gam, r_max=R_MAX):
    """phi(r_max; psi, K, a_Gam)."""
    phi0  = psi   # bezwymiarowe psi = phi(a_Gam) [nie mnozenie przez a_Gam!]
    # UWAGA: psi tu jest juz phi(a_Gam), bez skalowania przez a_Gam
    # Zachowujemy konwencje p98: psi to wartosc poczatkowa phi(a_Gam)
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
    """E[phi] = 4pi int [kfac/2 (dphi)^2 + (V-V1)] r^2 dr."""
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
        kfac     = 1.0 + ALPHA_K / phi
        Ek       = 4*np.pi * np.trapezoid(0.5*dphi**2 * kfac * r**2, r)
        Ep       = 4*np.pi * np.trapezoid((V_mod(phi) - V1) * r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan

def g_U(psi, K, a_Gam):
    """g = E/(4piK) - 1."""
    E = energy_soliton(psi, K, a_Gam)
    return E / (4*np.pi*K) - 1.0 if np.isfinite(E) else np.nan

def find_all_psi_zeros(K, a_Gam, psi_lo=PSI_LO, psi_hi=PSI_HI, n_psi=N_PSI):
    """Wszystkie psi: phi(r_max; psi,K,a_Gam) = 1."""
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv   = [phi_at_rmax(p, K, a_Gam) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv) - 1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not np.isfinite(fi) or not np.isfinite(fi1):
            continue
        if fi * fi1 < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p, K, a_Gam) - 1.0,
                            psis[i], psis[i+1], xtol=1e-7, maxiter=50)
                roots.append(pz)
            except Exception:
                pass
    return roots

def g_U_last(K, a_Gam):
    """g_U uzywajac ostatniego zera psi (Galaz U)."""
    roots = find_all_psi_zeros(K, a_Gam)
    if not roots:
        return np.nan, np.nan
    psi_last = roots[-1]
    return g_U(psi_last, K, a_Gam), psi_last

def find_Kstar(a_Gam, K_lo, K_hi, n_K=N_K, verbose=False):
    """
    Znajdz WSZYSTKIE zera g_U(K) w [K_lo, K_hi].
    Zwraca liste (K*, psi0_at_Kstar).
    """
    Ks  = np.linspace(K_lo, K_hi, n_K)
    gv  = []
    psv = []
    for K in Ks:
        g, ps = g_U_last(K, a_Gam)
        gv.append(g)
        psv.append(ps)

    zeros = []
    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if not np.isfinite(gi) or not np.isfinite(gi1):
            continue
        if gi * gi1 < 0:
            try:
                Kz = brentq(lambda K: g_U_last(K, a_Gam)[0],
                            Ks[i], Ks[i+1], xtol=1e-5, maxiter=40)
                _, pz = g_U_last(Kz, a_Gam)
                zeros.append((Kz, pz))
            except Exception:
                Kz_lin = Ks[i] - gv[i]*(Ks[i+1]-Ks[i])/(gv[i+1]-gv[i])
                _, pz = g_U_last(Kz_lin, a_Gam)
                zeros.append((Kz_lin, pz))

    if verbose:
        print(f"    Skan K [{K_lo:.4f},{K_hi:.4f}]: {len(zeros)} zer g_U")
        for Kz, pz in zeros:
            print(f"      K*={Kz:.6f}, psi0={pz:.4f}")

    return zeros

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 68)
    print("p100_K2_vs_agamma.py  --  TGP v1")
    print("OP-3 Sciezka 2: czy K*_2 istnieje dla a_Gam < 0.040?")
    print("=" * 68)
    print(f"alpha_K={ALPHA_K}, r_max={R_MAX}, psi in [{PSI_LO},{PSI_HI}]")
    print(f"a_Gam values: {AGAM_VALUES}")
    print()

    results = {}  # a_Gam -> {'K1': K1, 'K2': K2, 'psi1': ..., 'psi2': ...}

    for a_Gam in AGAM_VALUES:
        t0 = time.time()
        print(f"─── a_Gam = {a_Gam:.3f} ───────────────────────────────────")

        # Skalowane zakresy K
        scale = a_Gam / AGAM_PHYS
        K1_lo = K1_INF * scale * 0.5
        K1_hi = K1_INF * scale * 3.0
        K2_lo = K2_REF * scale * 0.3
        K2_hi = K2_REF * scale * 3.5

        print(f"  K*_1 skan: [{K1_lo:.5f}, {K1_hi:.5f}]")
        zeros1 = find_Kstar(a_Gam, K1_lo, K1_hi, n_K=30, verbose=True)

        print(f"  K*_2 skan: [{K2_lo:.5f}, {K2_hi:.5f}]")
        zeros2 = find_Kstar(a_Gam, K2_lo, K2_hi, n_K=N_K, verbose=True)

        K1  = zeros1[0][0]  if zeros1 else np.nan
        psi1= zeros1[0][1]  if zeros1 else np.nan
        K2  = zeros2[0][0]  if zeros2 else np.nan
        psi2= zeros2[0][1]  if zeros2 else np.nan

        dt = time.time() - t0
        results[a_Gam] = {'K1': K1, 'K2': K2, 'psi1': psi1, 'psi2': psi2}
        K1_s2 = f"{K1:.6f}" if np.isfinite(K1) else "BRAK"
        K2_s2 = f"{K2:.6f}" if np.isfinite(K2) else "BRAK"
        print(f"  => K*_1={K1_s2}, K*_2={K2_s2}  ({dt:.0f}s)")
        print()

    # ── tabela ──────────────────────────────────────────────────────────────
    print("=" * 68)
    print("TABELA: K*_1, K*_2 vs a_Gam")
    print("=" * 68)
    print(f"{'a_Gam':>8} {'K*_1':>12} {'K*_1/ref':>10} {'K*_2':>12} {'K*_2/K*_1':>12} {'K*_2 exist':>12}")
    print("-" * 70)
    for a_Gam in AGAM_VALUES:
        res   = results[a_Gam]
        K1    = res['K1']
        K2    = res['K2']
        scale = a_Gam / AGAM_PHYS
        K1_ref = K1_INF * scale
        K1_ratio = K1 / K1_ref if np.isfinite(K1) and K1_ref > 0 else np.nan
        K2K1  = K2 / K1 if np.isfinite(K1) and np.isfinite(K2) and K1 > 0 else np.nan
        K2_ex = "TAK" if np.isfinite(K2) else "NIE"
        K1_s  = f"{K1:.6f}" if np.isfinite(K1) else "  BRAK  "
        K2_s  = f"{K2:.6f}" if np.isfinite(K2) else "  BRAK  "
        K1r_s = f"{K1_ratio:.4f}" if np.isfinite(K1_ratio) else "  BRAK"
        K2K1_s= f"{K2K1:.3f}" if np.isfinite(K2K1) else "  BRAK"
        print(f"{a_Gam:>8.3f} {K1_s:>12} {K1r_s:>10} {K2_s:>12} {K2K1_s:>12} {K2_ex:>12}")

    # ── wyznaczenie a_c ─────────────────────────────────────────────────────
    print()
    K2_exists = {a: np.isfinite(results[a]['K2']) for a in AGAM_VALUES}
    print("K*_2 istnieje przy a_Gam:")
    for a in AGAM_VALUES:
        print(f"  a_Gam={a:.3f}: {'TAK' if K2_exists[a] else 'NIE'}")

    a_c = None
    sorted_a = sorted(AGAM_VALUES)
    for i, a in enumerate(sorted_a):
        if K2_exists[a]:
            a_c = a
            if i > 0:
                a_c_prev = sorted_a[i-1]
                print(f"\na_c(alpha_K) w przedziale ({a_c_prev:.3f}, {a:.3f}]")
            else:
                print(f"\na_c(alpha_K) <= {a:.3f} (K*_2 istnieje juz przy najmniejsym a_Gam)")
            break

    if a_c is None:
        print("\nK*_2 NIE istnieje przy zadnym a_Gam w badanym zakresie!")
        a_c = float('inf')

    # ══════════════════════════════════════════════════════════════════════
    # PASS / FAIL
    # ══════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 68)
    print("TESTY PASS/FAIL")
    print("=" * 68)
    passes = []

    # P1: K*_1 przy wszystkich a_Gam
    p1 = all(np.isfinite(results[a]['K1']) for a in AGAM_VALUES)
    passes.append(p1)
    n1 = sum(np.isfinite(results[a]['K1']) for a in AGAM_VALUES)
    print(f"P1 [K*_1 przy wszystkich a_Gam]: {'PASS' if p1 else 'FAIL'}  "
          f"({n1}/{len(AGAM_VALUES)})")

    # P2: K*_1 ~ ref * scale (Schwarzschild)
    deviations = []
    for a in AGAM_VALUES:
        K1 = results[a]['K1']
        ref = K1_INF * (a / AGAM_PHYS)
        if np.isfinite(K1) and ref > 0:
            deviations.append(abs(K1/ref - 1))
    if deviations:
        max_dev = max(deviations)
        p2 = max_dev < 0.10  # 10% tolerancja (bias numeryczny z r_max=80)
        print(f"P2 [K*_1 ~ K_ref * scale]:       {'PASS' if p2 else 'FAIL'}  "
              f"(max odch={max_dev*100:.1f}%, prog=10%)")
    else:
        p2 = False
        print(f"P2 [K*_1 ~ K_ref * scale]:       FAIL  (brak danych)")
    passes.append(p2)

    # P3: K*_2 istnieje przy a_Gam = 0.040
    p3 = np.isfinite(results[AGAM_PHYS]['K2'])
    passes.append(p3)
    K2_phys = results[AGAM_PHYS]['K2']
    print(f"P3 [K*_2 przy a_Gam=0.040]:      {'PASS' if p3 else 'FAIL'}  "
          f"(K*_2={'%.5f' % K2_phys if p3 else 'BRAK'})")

    # P4: K*_2 nie istnieje dla a_Gam < a_Gam_phys
    a_below_phys = [a for a in AGAM_VALUES if a < AGAM_PHYS]
    p4 = all(not K2_exists[a] for a in a_below_phys) if a_below_phys else False
    n_below = len(a_below_phys)
    n_no_K2 = sum(not K2_exists[a] for a in a_below_phys)
    passes.append(p4)
    print(f"P4 [K*_2 brak dla a_Gam<0.040]:  {'PASS' if p4 else 'FAIL'}  "
          f"(brak K*_2 przy {n_no_K2}/{n_below} punktach ponizej 0.040)")

    # P5: a_c ~ 0.040
    if p3 and p4 and a_below_phys:
        a_just_below = max(a for a in a_below_phys)
        # a_c jest w (a_just_below, AGAM_PHYS]
        p5 = abs(AGAM_PHYS - a_just_below) / AGAM_PHYS < 0.25
        print(f"P5 [a_c ~ 0.040]:                {'PASS' if p5 else 'FAIL'}  "
              f"(a_c w ({a_just_below:.3f},{AGAM_PHYS:.3f}])")
    elif p3 and not a_below_phys:
        p5 = False
        print(f"P5 [a_c ~ 0.040]:                FAIL  (brak punkt. ponizej 0.040)")
    else:
        p5 = False
        print(f"P5 [a_c ~ 0.040]:                FAIL  (K*_2 brak przy 0.040 lub K*_2 wsz.)")
    passes.append(p5)

    n_pass = sum(passes)
    print(f"\nWYNIK: {n_pass}/{len(passes)} PASS")

    # ── podsumowanie ─────────────────────────────────────────────────────────
    print("\n" + "=" * 68)
    print("PODSUMOWANIE FIZYCZNE  --  OP-3 Sciezka 2")
    print("=" * 68)
    if p3 and p4:
        print(f"  ** K*_2 pojawia sie przy a_Gam >= a_c, gdzie a_c w ({a_just_below:.3f},{AGAM_PHYS:.3f}]")
        if p5:
            print(f"  ** a_c ~ a_Gam_fizyczne = 0.040  =>  OP-3 Sciezka 2 POTWIERDZONA!")
        else:
            print(f"  ** a_c daleko od 0.040  =>  OP-3 Sciezka 2 NIE potwierdzona")
    elif p3 and not p4:
        print(f"  K*_2 istnieje rowniez dla a_Gam < 0.040  =>  brak krytycznej a_c przy 0.040")
        print(f"  OP-3 Sciezka 2 NIE potwierdzona (K*_2 istnieje wszedzie)")
    else:
        print(f"  Blad numeryczny lub K*_2 nie znaleziono nawet przy 0.040")
    print("=" * 68)
