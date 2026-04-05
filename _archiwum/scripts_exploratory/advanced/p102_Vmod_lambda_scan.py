#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p102_Vmod_lambda_scan.py  --  TGP v1
======================================
Cel: Sciezka 1 Koide — skan parametru lambda potencjalu V_mod.
     Pytanie: jaka wartosc lambda daje K*_2/K*_1 = 206.77 (stosunek m_mu/m_e)?

Motywacja:
  Standardowe V_mod (lambda=5.50e-6) daje K*_2/K*_1 ~16.7 (p101, gal. psi_0~4.1)
  lub ~21-25 (p99, gal. psi_0~5.0).
  Cel: K*_2/K*_1 = r_21 = 206.77 (Koide / PDG: m_mu/m_e).
  Hipoteza: zwiekszenie lambda zwieksza energie solitonu K*_2 (duze psi_0 => duze
  (phi-1)^6 wkl.) przy prawie stalym K*_1 (psi_0~1.24 => (0.24)^6 ~ 3e-4, male).

V_mod(phi; lambda) = gamma/3 * phi^3 - gamma/4 * phi^4 + lambda/6 * (phi-1)^6
dV_mod(phi; lambda) = gamma*phi^2 - gamma*phi^3 + lambda*(phi-1)^5

Metoda:
  - Skan lambda w [1e-7, 1e-2] (log scale, 20 punktow)
  - a_Gam = 0.040, alpha_K = 8.5616, r_max = 80 (jak p101)
  - Dla kazdego lambda: K*_1 (psi in [1.001, 1.80]) i K*_2 (psi in [2.0, 20.0])
  - Oblicz K*_2/K*_1(lambda) i interpoluj lambda* dla K*_2/K*_1 = 206.77

Testy:
  P1: K*_1 stabilne (max odch. < 5%) przez caly zakres lambda
  P2: K*_2 istnieje przy lambda_phys = 5.5e-6 (walidacja z p101)
  P3: K*_2/K*_1 przy lambda_phys ~16.7 (konsystencja z p101)
  P4: K*_2/K*_1 jest monotonicznie rosnaca funkcja lambda
  P5: lambda* interpolowane dla K*_2/K*_1 = 206.77 (lub ograniczenie dolne)

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
ALPHA_K    = 8.5616
GAM        = 1.0
V1         = GAM/3.0 - GAM/4.0      # = 1/12 (vacuum; niezalezne od lambda)
LAM_PHYS   = 5.501357e-06
K1_INF     = 0.010029               # K*_1(inf) przy a_Gam=0.040 (p97)
K2_P101    = 0.171272               # K*_2 przy lambda_phys, a_Gam=0.040 (p101)
AGAM       = 0.040
R_MAX      = 80.0
R21_TARGET = 206.77                 # m_mu / m_e (Koide / PDG)

# ─── skan lambda ─────────────────────────────────────────────────────────────
# Zakres 1e-7 do 1e-2 (5 dekad, 20 punktow)
LAM_VALUES = np.logspace(-7, -2, 20)

# ─── parametry numeryczne ─────────────────────────────────────────────────────
PSI1_LO = 1.10;   PSI1_HI = 1.80    # dla K*_1 (psi_0 ~ 1.24; psi<1.10 to spurialne rozw. przy prozni)
PSI2_LO = 3.50;   PSI2_HI = 7.00    # dla K*_2 (psi_0 ~ 4.1; szerszy = spurialne zera w [2,3.5])
N_PSI1  = 50
N_PSI2  = 60
N_K1    = 25
N_K2    = 45
N_EVAL  = 600
RTOL    = 1e-8
ATOL    = 1e-10

# K*_2 skan: od K2_P101*0.25 do K2_P101*4 (pokrywa K*_2 +/- 4x od wartosci p101)
# Uzasadnienie: zmiana lambda < x1000 przesuwa K*_2 o < 5% (szacunki analityczne)
K2_LO  = K2_P101 * 0.25    # ~ 0.043
K2_HI  = K2_P101 * 4.0     # ~ 0.685  (r21_max ~ 66; jesli K*_2/K*_1 > 66, rozszerz)

# ─── ODE i V_mod (parametryzowane przez lambda) ───────────────────────────────
def dV_mod_lam(phi, lam):
    return GAM*phi**2 - GAM*phi**3 + lam*(phi - 1.0)**5

def V_mod_lam(phi, lam):
    return GAM/3.0*phi**3 - GAM/4.0*phi**4 + lam/6.0*(phi - 1.0)**6

def ode_rhs_lam(r, y, lam):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K / phi
    phi2 = max(phi, 1e-10)**2
    return [dphi,
            dV_mod_lam(phi, lam)/kfac
            + ALPHA_K*dphi**2/(2.0*phi2*kfac)
            - 2.0/r*dphi]

# ─── strzelanie: phi na r_max ─────────────────────────────────────────────────
def phi_at_rmax(psi, K, lam, r_max=R_MAX):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (r_max / AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs_lam(r, y, lam),
                        [AGAM, r_max], [psi, dphi0],
                        method='DOP853', rtol=RTOL, atol=ATOL,
                        t_eval=r_ev, dense_output=False)
        if sol.t[-1] < r_max * 0.99:
            return np.nan
        return float(sol.y[0, -1])
    except Exception:
        return np.nan

# ─── energia solitonu ──────────────────────────────────────────────────────────
def energy_soliton(psi, K, lam, r_max=R_MAX):
    dphi0 = -K / AGAM**2
    r_ev  = AGAM * (r_max / AGAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs_lam(r, y, lam),
                        [AGAM, r_max], [psi, dphi0],
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
        Ep   = 4*np.pi * np.trapezoid((V_mod_lam(phi, lam) - V1) * r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan

def g_at_psi(psi, K, lam):
    E = energy_soliton(psi, K, lam)
    return E / (4.0*np.pi*K) - 1.0 if (np.isfinite(E) and K > 0) else np.nan

# ─── znajdz pierwsze zero psi w oknie ────────────────────────────────────────
def find_psi_zero_in_window(K, lam, psi_lo, psi_hi, n_psi):
    """Pierwsze zero phi(r_max)=1 w [psi_lo, psi_hi]. Zwraca psi_zero lub nan."""
    psis = np.linspace(psi_lo, psi_hi, n_psi)
    Fv   = [phi_at_rmax(p, K, lam) - 1.0 for p in psis]
    for i in range(len(Fv) - 1):
        fi, fi1 = Fv[i], Fv[i+1]
        if not np.isfinite(fi) or not np.isfinite(fi1):
            continue
        if fi * fi1 < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p, K, lam) - 1.0,
                            psis[i], psis[i+1], xtol=1e-7, maxiter=50)
                return pz
            except Exception:
                pass
    return np.nan

def g_U_narrow(K, lam, psi_lo, psi_hi, n_psi):
    pz = find_psi_zero_in_window(K, lam, psi_lo, psi_hi, n_psi)
    if np.isnan(pz):
        return np.nan, np.nan
    return g_at_psi(pz, K, lam), pz

# ─── skan K i znajdz K* ───────────────────────────────────────────────────────
def find_Kstar(lam, K_lo, K_hi, psi_lo, psi_hi, n_K, n_psi, verbose=False):
    """Znajdz K* w [K_lo, K_hi] uzywajac okna psi [psi_lo, psi_hi]."""
    Ks  = np.linspace(K_lo, K_hi, n_K)
    gv  = []
    psv = []
    for K in Ks:
        g, pz = g_U_narrow(K, lam, psi_lo, psi_hi, n_psi)
        gv.append(g)
        psv.append(pz)

    # Szukaj zmian znaku g_U
    sign_changes = []
    for i in range(len(gv) - 1):
        gi, gi1 = gv[i], gv[i+1]
        if np.isfinite(gi) and np.isfinite(gi1) and gi * gi1 < 0:
            sign_changes.append((Ks[i], Ks[i+1], gv[i], gv[i+1], psv[i]))

    if verbose:
        n_valid = sum(np.isfinite(g) for g in gv)
        print(f"        g_U valid: {n_valid}/{len(gv)}, "
              f"zmian znaku: {len(sign_changes)}")
        # Pokaz kilka wartosci g_U dla diagnozy
        valid_pairs = [(Ks[i], gv[i], psv[i]) for i in range(len(gv))
                       if np.isfinite(gv[i])]
        if valid_pairs and len(valid_pairs) <= 8:
            for K_v, g_v, pz_v in valid_pairs:
                print(f"          K={K_v:.5f} g={g_v:.4f} psi={pz_v:.4f}")
        elif valid_pairs:
            for K_v, g_v, pz_v in valid_pairs[:4]:
                print(f"          K={K_v:.5f} g={g_v:.4f} psi={pz_v:.4f}")
            print(f"          ... ({len(valid_pairs)-4} wiecej) ...")
        if sign_changes:
            Kl, Kr, gl, gr, pz0 = sign_changes[0]
            print(f"        1. zmiana znaku: K=[{Kl:.5f},{Kr:.5f}], "
                  f"g=[{gl:.4f},{gr:.4f}], psi~{pz0:.4f}")

    if not sign_changes:
        return np.nan, np.nan

    # Bierz pierwsze zero
    Kl, Kr, gl, gr, _ = sign_changes[0]
    try:
        Kstar = brentq(
            lambda K: g_U_narrow(K, lam, psi_lo, psi_hi, n_psi)[0],
            Kl, Kr, xtol=1e-5, maxiter=40)
    except Exception:
        Kstar = Kl - gl * (Kr - Kl) / (gr - gl)

    _, psi_at_K = g_U_narrow(Kstar, lam, psi_lo, psi_hi, n_psi)
    return Kstar, psi_at_K

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=" * 70)
    print("p102_Vmod_lambda_scan.py  --  TGP v1")
    print("Sciezka 1 Koide: K*_2/K*_1 vs lambda")
    print("=" * 70)
    print(f"alpha_K={ALPHA_K}, a_Gam={AGAM}, r_max={R_MAX}")
    print(f"K*_1 psi okno: [{PSI1_LO},{PSI1_HI}], K skan: [{K1_INF*0.4:.4f},{K1_INF*2.5:.4f}]")
    print(f"K*_2 psi okno: [{PSI2_LO},{PSI2_HI}], K skan: [{K2_LO:.4f},{K2_HI:.4f}]")
    print(f"lambda skan: {LAM_VALUES[0]:.1e} do {LAM_VALUES[-1]:.1e} ({len(LAM_VALUES)} pkt)")
    print(f"Cel: K*_2/K*_1 = {R21_TARGET}")
    print()

    results = {}

    t_total = time.time()
    for lam in LAM_VALUES:
        t0 = time.time()
        print(f"--- lambda = {lam:.3e} ---")

        # K*_1 skan
        K1_lo = K1_INF * 0.40
        K1_hi = K1_INF * 2.50
        K1, psi1 = find_Kstar(lam, K1_lo, K1_hi, PSI1_LO, PSI1_HI,
                               N_K1, N_PSI1, verbose=True)

        # K*_2 skan (szeroki zakres K)
        K2, psi2 = find_Kstar(lam, K2_LO, K2_HI, PSI2_LO, PSI2_HI,
                               N_K2, N_PSI2, verbose=True)

        dt  = time.time() - t0
        r21 = (K2/K1) if (np.isfinite(K1) and np.isfinite(K2) and K1 > 0) else np.nan

        results[lam] = dict(K1=K1, K2=K2, psi1=psi1, psi2=psi2, r21=r21)

        K1s  = f"{K1:.6f}"  if np.isfinite(K1)  else "BRAK"
        K2s  = f"{K2:.6f}"  if np.isfinite(K2)  else "BRAK"
        p1s  = f"{psi1:.4f}" if np.isfinite(psi1) else "nan"
        p2s  = f"{psi2:.4f}" if np.isfinite(psi2) else "nan"
        r21s = f"{r21:.2f}" if np.isfinite(r21) else "BRAK"
        print(f"  => K*_1={K1s} (psi={p1s}), K*_2={K2s} (psi={p2s}), r21={r21s}  [{dt:.0f}s]")
        print()

    dt_total = time.time() - t_total
    print(f"Czas laczny: {dt_total:.0f}s ({dt_total/60:.1f} min)\n")

    # ── tabela wynikow ────────────────────────────────────────────────────────
    print("=" * 70)
    print("TABELA: K*_1, K*_2, K*_2/K*_1 vs lambda")
    print("=" * 70)
    print(f"{'lambda':>12} {'K*_1':>12} {'K*_1/ref':>10} {'K*_2':>12} "
          f"{'psi_0(K*2)':>12} {'r21':>10}")
    print("-" * 70)
    for lam in LAM_VALUES:
        res  = results[lam]
        K1   = res['K1']; K2 = res['K2']
        psi2 = res['psi2']; r21 = res['r21']
        ref1 = K1_INF   # referencyjna K*_1 przy lambda_phys
        r1   = K1/ref1  if np.isfinite(K1) else np.nan
        K1s  = f"{K1:.6f}"  if np.isfinite(K1)  else "  BRAK  "
        K2s  = f"{K2:.6f}"  if np.isfinite(K2)  else "  BRAK  "
        r1s  = f"{r1:.4f}"  if np.isfinite(r1)  else " BRAK"
        r21s = f"{r21:.3f}" if np.isfinite(r21) else " BRAK"
        p2s  = f"{psi2:.4f}" if np.isfinite(psi2) else "  nan  "
        print(f"{lam:>12.3e} {K1s:>12} {r1s:>10} {K2s:>12} {p2s:>12} {r21s:>10}")

    # ── analiza trendu K*_2/K*_1(lambda) ─────────────────────────────────────
    print()
    lam_valid = [lam for lam in LAM_VALUES if np.isfinite(results[lam]['r21'])]
    r21_valid = [results[lam]['r21'] for lam in lam_valid]

    print("Analiza trendu K*_2/K*_1(lambda):")
    if len(r21_valid) < 2:
        print("  Zbyt malo danych do analizy trendu.")
        trend_monotone = False
        r21_max        = max(r21_valid) if r21_valid else np.nan
        lam_at_max     = lam_valid[np.argmax(r21_valid)] if r21_valid else np.nan
    else:
        r21_arr  = np.array(r21_valid)
        lam_arr  = np.array(lam_valid)
        r21_max  = np.nanmax(r21_arr)
        lam_at_max = lam_arr[np.argmax(r21_arr)]
        # Sprawdz monotonnosc
        diffs = np.diff(r21_arr)
        n_pos = np.sum(diffs > 0)
        n_neg = np.sum(diffs < 0)
        trend_monotone = (n_neg == 0) or (n_pos == 0)
        print(f"  Zakres r21: [{r21_arr.min():.3f}, {r21_arr.max():.3f}]")
        print(f"  Max r21={r21_max:.3f} przy lambda={lam_at_max:.3e}")
        print(f"  Monotonnosc: {'TAK' if trend_monotone else 'NIE'} "
              f"(rosnace={n_pos}, malejace={n_neg})")

    # ── interpolacja lambda* dla r21 = R21_TARGET ─────────────────────────────
    print()
    lam_star = np.nan
    if len(r21_valid) >= 2:
        r21_arr = np.array(r21_valid)
        lam_arr = np.array(lam_valid)
        # Szukaj przejscia przez R21_TARGET
        for i in range(len(r21_arr) - 1):
            a, b = r21_arr[i], r21_arr[i+1]
            if (a - R21_TARGET) * (b - R21_TARGET) < 0:
                # interpolacja liniowa w log(lambda)
                la, lb = np.log10(lam_arr[i]), np.log10(lam_arr[i+1])
                lam_star = 10**(la + (lb - la) * (R21_TARGET - a) / (b - a))
                print(f"  K*_2/K*_1 = {R21_TARGET} przy lambda* ~ {lam_star:.4e}")
                print(f"  (log10(lambda*) = {np.log10(lam_star):.3f})")
                break

        if np.isnan(lam_star):
            if r21_max < R21_TARGET:
                print(f"  BRAK przeciecia: max r21={r21_max:.2f} < {R21_TARGET}")
                print(f"  K*_2/K*_1=206.77 NIE osiagalne w zakresie lambda=[1e-7,1e-2]")
                print(f"  Wymagana lambda > {lam_at_max:.3e} lub inne modyfikacje V_mod")
            else:
                print(f"  BRAK przeciecia mimo max r21={r21_max:.2f} > {R21_TARGET}")
                print(f"  (numeryczne problemy lub niemonotonosc)")

    # ── PASS/FAIL ─────────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("TESTY PASS/FAIL")
    print("=" * 70)
    passes = []

    # P1: K*_1 stabilne (max odch < 5%)
    K1_vals = [results[lam]['K1'] for lam in LAM_VALUES if np.isfinite(results[lam]['K1'])]
    if K1_vals:
        devs = [abs(k/K1_INF - 1) for k in K1_vals]
        mx   = max(devs)
        p1   = (mx < 0.05) and (len(K1_vals) == len(LAM_VALUES))
        passes.append(p1)
        print(f"P1 [K*_1 stabilne, max dev<5%]:  {'PASS' if p1 else 'FAIL'}  "
              f"(max odch={mx*100:.1f}%, znaleziono {len(K1_vals)}/{len(LAM_VALUES)})")
    else:
        passes.append(False)
        print(f"P1 [K*_1 stabilne]:              FAIL  (brak K*_1)")

    # P2: K*_2 przy lambda_phys
    # znajdz najblizszy lambda do LAM_PHYS
    idx_phys = np.argmin(np.abs(np.array(LAM_VALUES) - LAM_PHYS))
    lam_near = LAM_VALUES[idx_phys]
    K2_near  = results[lam_near]['K2']
    p2 = np.isfinite(K2_near)
    passes.append(p2)
    K2ns = f"{K2_near:.5f}" if p2 else "BRAK"
    print(f"P2 [K*_2 przy lambda~phys]:      {'PASS' if p2 else 'FAIL'}  "
          f"(lambda={lam_near:.2e}, K*_2={K2ns})")

    # P3: K*_2/K*_1 przy lambda~phys ~16.7
    r21_near = results[lam_near]['r21']
    p3 = np.isfinite(r21_near) and (10.0 < r21_near < 30.0)
    passes.append(p3)
    r21ns = f"{r21_near:.3f}" if np.isfinite(r21_near) else "BRAK"
    print(f"P3 [r21 przy lambda_phys~16.7]:  {'PASS' if p3 else 'FAIL'}  "
          f"(r21={r21ns}, prog=[10,30])")

    # P4: K*_2/K*_1 monotonicznie rosnace z lambda
    p4 = trend_monotone and (len(r21_valid) >= 3) and (n_pos > n_neg if len(r21_valid)>=2 else False)
    passes.append(p4)
    print(f"P4 [r21 monotonicznie rosnace]:  {'PASS' if p4 else 'FAIL'}")

    # P5: lambda* dla r21=206.77 albo jasne ograniczenie
    p5 = np.isfinite(lam_star) or (len(r21_valid) >= 5)
    passes.append(p5)
    lss = f"{lam_star:.4e}" if np.isfinite(lam_star) else f"brak (max r21={r21_max:.1f})"
    print(f"P5 [lambda* lub ograniczenie]:   {'PASS' if p5 else 'FAIL'}  "
          f"(lambda*={lss})")

    print(f"\nWYNIK: {sum(passes)}/{len(passes)} PASS")

    # ── podsumowanie fizyczne ────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("PODSUMOWANIE FIZYCZNE  --  Sciezka 1 Koide")
    print("=" * 70)
    r21_near_s = f"{r21_near:.1f}" if np.isfinite(r21_near) else "?"
    print(f"Standard lambda_phys={LAM_PHYS:.3e}: K*_2/K*_1 ~ {r21_near_s}")
    print(f"Cel Koide: K*_2/K*_1 = {R21_TARGET}")
    if np.isfinite(lam_star):
        print(f"=> lambda* = {lam_star:.4e}  ({lam_star/LAM_PHYS:.1f} * lambda_phys)")
        print(f"   Modyfikacja V_mod: zmiana lambda z {LAM_PHYS:.2e} -> {lam_star:.2e}")
        print(f"   Pytanie: czy ta lambda ma fizyczne uzasadnienie w TGP?")
    else:
        if len(r21_valid) > 0:
            print(f"Max K*_2/K*_1 = {r21_max:.2f} przy lambda = {lam_at_max:.3e}")
            print(f"=> V_mod klasy (phi-1)^6 nie odtwarza r21=206.77 w zakresie lambda=[1e-7,1e-2]")
            print(f"   Potrzebna zmiana formy V_mod (np. inna potega (phi-1)^n, n>6)")
        else:
            print("Brak danych — problemy numeryczne.")
    print("=" * 70)
