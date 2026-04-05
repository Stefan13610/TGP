#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p77_K2_scan.py  --  TGP v1  *  OP-3 Path 1 (Koide) — weryfikacja K*_2, K*_3
=============================================================================
Cel: Sprawdzic czy w aktualnym V_mod (bez modyfikacji) istnieja K*_2, K*_3
     jako samospojne solitony (g = E/(4*pi*K) - 1 = 0).

Motywacja: Sekcja D skryptu p76_alpha_K_bifurcation.py wykryla zmiany znaku
g(K) miedzy:
  - K=0.0570 (g=+2.42) i K=0.0889 (g=-0.81)  → potencjalny K*_2 ~ 0.07
  - K=0.0889 (g=-0.81) i K=0.1387 (g=+2.05)  → potencjalny K*_3 ~ 0.12

Pytanie kluczowe: Czy sa to PRAWDZIWE zera g(K), czy artefakty przełaczenia
galezi w algorytmie compute_g (kryterium min|g|)?

Metoda:
  1. Gesты skan K ∈ [0.015, 0.25] przy alpha=8.5616, a_Gam=0.040
     -- 60 punktow K (logarytmicznie rownomierne)
  2. Dla kazdego K: compute_g zwraca g_best (min|g|) i lista g_all (wszystkie
     galezi). Drukujemy obie.
  3. Jesli znaki g_best zmieniaja sie monotonicznie (nie przeskakuja) → zera
     sa prawdziwe.
  4. Jesli g_best skacze z powrotem → artefakt przełaczenia galezi.
  5. Jezeli K*_2, K*_3 istnieja: oblicz Q_TGP = Koide ratio.

Testy PASS/FAIL:
  P1: g(K) zmienia znak w zakresie [0.015, 0.25]                (K*_2 istnieje?)
  P2: Zmiana znaku jest MONOTONICZNA (nie artefakt)
  P3: K*_2 znaleziony dokladnie (brentq)
  P4: |K*_2 / K*_1 - R21_PDG| / R21_PDG < 30%                  (skala muon?)
  P5: K*_3 znaleziony; Q_TGP(K*_1, K*_2, K*_3) obliczone

Data: 2026-03-25
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# ─────────────────────────────────────────────────────────────────────────────
# Stale TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K_REF  = 8.5616
A_GAM        = 0.040
K1_REF       = 0.010414
R21_PDG      = 206.77        # m_mu / m_e
R31_PDG      = 1776.86 / 0.51100  # m_tau / m_e
LAM          = 5.501357e-06
GAMMA_V      = 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Potencjal V_mod  (STANDARDOWY — bez modyfikacji)
# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):
    return GAMMA_V/3 * p**3 - GAMMA_V/4 * p**4 + LAM/6 * (p - 1)**6

def dV_mod(p):
    return GAMMA_V * p**2 - GAMMA_V * p**3 + LAM * (p - 1)**5

V1 = GAMMA_V/3 - GAMMA_V/4   # = 1/12


def ode_rhs(r, y):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + ALPHA_K_REF / phi
    ddphi = (dV_mod(phi) / kfac
             + ALPHA_K_REF * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]


# ─────────────────────────────────────────────────────────────────────────────
# Funkcja pomocnicza: phi(r_max) dla danego (psi_core, K)
# ─────────────────────────────────────────────────────────────────────────────
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA_K_REF)
R_MAX  = max(40.0, A_GAM * 5, 6.0 / M_EFF)  # stale dla wszystkich K


def phi_at_rmax(psi_core, K, r_max=R_MAX, rtol=1e-7, atol=1e-9, n_eval=800):
    dphi0 = -K / A_GAM**2

    def ev_lo(r, y):
        return y[0] - 1e-5
    ev_lo.terminal  = True
    ev_lo.direction = -1

    r_eval = A_GAM * (r_max / A_GAM) ** np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(
            ode_rhs, [A_GAM, r_max], [psi_core, dphi0],
            method='DOP853', rtol=rtol, atol=atol,
            t_eval=r_eval, events=[ev_lo]
        )
        if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
            return np.nan
        return float(sol.y[0, -1])
    except Exception:
        return np.nan


def energy(psi_z, K, r_max=R_MAX):
    """Oblicza energie E = Ek + Ep dla profilu phi(psi_z, K)."""
    dphi0 = -K / A_GAM**2

    def ev_lo(r, y):
        return y[0] - 1e-5
    ev_lo.terminal  = True
    ev_lo.direction = -1

    r_ev = A_GAM * (r_max / A_GAM) ** np.linspace(0, 1, 4000)
    try:
        sol = solve_ivp(
            ode_rhs, [A_GAM, r_max], [psi_z, dphi0],
            method='DOP853', rtol=1e-9, atol=1e-11,
            t_eval=r_ev, events=[ev_lo]
        )
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        Ek = 4*np.pi * np.trapezoid(0.5*dphi**2*(1+ALPHA_K_REF/phi)*r**2, r)
        Ep = 4*np.pi * np.trapezoid((V_mod(phi) - V1)*r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan


def compute_g_all(K, n_psi=120, r_max=R_MAX):
    """
    Dla danego K:
      - Skan psi_core in [1.001, psi_max] (n_psi punktow)
      - Znajdz WSZYSTKIE galezi (zmiany znaku F)
      - Dla kazdej gałęzi oblicz g = E/(4piK) - 1
      - Zwroc (g_best, g_all, psi_all) gdzie g_best = min|g|, g_all = lista
    """
    psi_max  = max(3.0, 1.0 + 2.0 * K / A_GAM)
    psi_vals = np.linspace(1.001, min(psi_max, 300.0), n_psi)

    F_vals = np.array([phi_at_rmax(p, K, r_max) - 1.0 for p in psi_vals])

    branches = []
    for i in range(len(F_vals) - 1):
        Fi, Fj = F_vals[i], F_vals[i+1]
        if np.isfinite(Fi) and np.isfinite(Fj) and Fi * Fj < 0:
            try:
                psi_z = brentq(
                    lambda p: phi_at_rmax(p, K, r_max) - 1.0,
                    psi_vals[i], psi_vals[i+1],
                    xtol=1e-5, rtol=1e-5, maxiter=40
                )
                E = energy(psi_z, K, r_max)
                if np.isfinite(E):
                    g = E / (4*np.pi*K) - 1.0
                    branches.append((psi_z, g))
            except Exception:
                pass

    if not branches:
        return np.nan, [], []

    g_all   = [b[1] for b in branches]
    psi_all = [b[0] for b in branches]
    best    = min(branches, key=lambda x: abs(x[1]))
    return best[1], g_all, psi_all


# ─────────────────────────────────────────────────────────────────────────────
# NAGLOWEK
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 72)
print("TGP v1  *  p77_K2_scan.py  --  OP-3 Path 1: weryfikacja K*_2, K*_3")
print("alpha = 8.5616,  a_Gam = 0.040,  V_mod STANDARD (bez modyfikacji)")
print("=" * 72)
print(f"\n[PARAMETRY]")
print(f"  K*_1_ref = {K1_REF:.6f}")
print(f"  R21_PDG  = {R21_PDG:.2f}  (m_mu/m_e)")
print(f"  R31_PDG  = {R31_PDG:.2f}  (m_tau/m_e)")
print(f"  R_MAX    = {R_MAX:.2f}  (granica calkowania)")
print(f"  M_eff    = {M_EFF:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA A: Gesты skan g(K) w [0.015, 0.25] — 60 punktow
#           Drukuj g_best I wszystkie galezi
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SEKCJA A: Gesты skan g(K) — 60 punktow w [0.015, 0.25]")
print("          Sprawdzamy ciaglosc i galezi przy kazdym K")
print("=" * 72)

K_scan_A = np.logspace(np.log10(0.015), np.log10(0.25), 60)

print(f"\n  {'K':>10}  {'g_best':>9}  {'n_gal':>6}  {'g_all (wszystkie gałęzie)':}")
print("  " + "-" * 70)

g_best_arr  = np.full(len(K_scan_A), np.nan)
g_all_arr   = []
psi_all_arr = []

for idx, K in enumerate(K_scan_A):
    g_best, g_all, psi_all = compute_g_all(K, n_psi=120)
    g_best_arr[idx] = g_best
    g_all_arr.append(g_all)
    psi_all_arr.append(psi_all)

    n_g = len(g_all)
    g_all_str = "  ".join(f"{g:.4f}(ψ={p:.3f})"
                          for g, p in zip(g_all, psi_all)) if n_g else "---"
    g_best_str = f"{g_best:.4f}" if np.isfinite(g_best) else "---"
    print(f"  {K:>10.5f}  {g_best_str:>9}  {n_g:>6}  {g_all_str}")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA B: Analiza zmian znaku g_best(K)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SEKCJA B: Analiza zmian znaku g_best(K)")
print("=" * 72)

sign_changes = []
for i in range(len(g_best_arr) - 1):
    gi, gj = g_best_arr[i], g_best_arr[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
        sign_changes.append((i, K_scan_A[i], K_scan_A[i+1], gi, gj))

print(f"\n  Znaleziono {len(sign_changes)} zmian znaku g_best(K):")
for idx, (i, Ka, Kb, ga, gb) in enumerate(sign_changes):
    print(f"\n  [{idx+1}] K ∈ [{Ka:.5f}, {Kb:.5f}]")
    print(f"       g_best({Ka:.5f}) = {ga:.4f}")
    print(f"       g_best({Kb:.5f}) = {gb:.4f}")
    print(f"       --> sugerowany K*_{idx+2} ≈ {0.5*(Ka+Kb):.5f}")

    # Sprawdz ciaglosc galezi:
    n_a = len(g_all_arr[i])
    n_b = len(g_all_arr[i+1])
    print(f"       Liczba galezi: {n_a} przy Ka, {n_b} przy Kb")
    if n_a == 1 and n_b == 1:
        print(f"       -> 1 galaz po obu stronach — zmiana MONOTONCZNA (brak przeskok)")
    elif n_a != n_b:
        print(f"       -> UWAGA: rozna liczba galezi ({n_a} vs {n_b}) — mozliwy przeskok!")
    else:
        print(f"       -> Ta sama liczba galezi — sprawdz recznic psi_z")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA C: Dokladne wyznaczenie K*_2 (i K*_3) przez brentq
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SEKCJA C: Dokladne wyznaczenie K*_2, K*_3 (brentq na g_best)")
print("=" * 72)

K_stars = [K1_REF]  # K*_1 z referencji

def g_best_func(K):
    g, _, _ = compute_g_all(K, n_psi=120)
    return g if np.isfinite(g) else np.nan

for idx, (i, Ka, Kb, ga, gb) in enumerate(sign_changes):
    print(f"\n  Szukam K*_{idx+2} w [{Ka:.5f}, {Kb:.5f}] ...")
    try:
        K_star = brentq(g_best_func, Ka, Kb,
                        xtol=Ka * 0.005, rtol=0.005, maxiter=25)
        g_at_star, g_all_s, psi_all_s = compute_g_all(K_star, n_psi=120)
        print(f"  K*_{idx+2} = {K_star:.6f}  (g = {g_at_star:.2e})")
        if psi_all_s:
            print(f"  psi_core* = {psi_all_s[0]:.4f}")
        print(f"  K*_{idx+2} / K*_1_ref = {K_star / K1_REF:.2f}")
        K_stars.append(K_star)
    except Exception as e:
        print(f"  Blad brentq: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA D: Koide ratio Q_TGP
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("SEKCJA D: Koide ratio Q_TGP(K*_1, K*_2, K*_3)")
print("=" * 72)

def koide_Q(E1, E2, E3):
    """Q = (E1+E2+E3) / (sqrt(E1)+sqrt(E2)+sqrt(E3))^2"""
    s_sum = np.sqrt(E1) + np.sqrt(E2) + np.sqrt(E3)
    return (E1 + E2 + E3) / s_sum**2

print(f"\n  Znalezione K* (masy solitonow): {K_stars}")
print()

if len(K_stars) >= 3:
    K1, K2, K3 = K_stars[0], K_stars[1], K_stars[2]

    Q = koide_Q(K1, K2, K3)
    R21_num = K2 / K1
    R31_num = K3 / K1

    print(f"  K*_1 = {K1:.6f}  (referencja)")
    print(f"  K*_2 = {K2:.6f}")
    print(f"  K*_3 = {K3:.6f}")
    print()
    print(f"  K*_2 / K*_1 = {R21_num:.2f}  (PDG: {R21_PDG:.2f})")
    print(f"  K*_3 / K*_1 = {R31_num:.2f}  (PDG: {R31_PDG:.2f})")
    print()
    print(f"  Q_TGP = {Q:.6f}  (Koide: 2/3 = {2/3:.6f})")
    print(f"  odch(Q_TGP - 2/3) = {abs(Q - 2/3)*100 / (2/3):.2f}%")
    print()

    if abs(Q - 2/3) < 0.01:
        print("  [!!!] Q_TGP blisko 2/3 -- solitony PASUJA do formuly Koidego!")
    else:
        print(f"  Q_TGP = {Q:.4f} -- nie odpowiada Koide 2/3.")
        print(f"  Masy solitonow nie odpowiadaja masam leptonow:")
        print(f"  R21_num = {R21_num:.1f} vs R21_PDG = {R21_PDG:.1f}")
        print(f"  Solitony V_mod STANDARD nie daja hierarchii leptonow.")
        print(f"  (Wymaga modyfikacji V_mod -- p15 / Path 1 dalej otwarta)")

elif len(K_stars) == 2:
    print(f"  Znaleziono tylko K*_1, K*_2 -- brak K*_3 dla testu Koidego.")
    print(f"  K*_1 = {K_stars[0]:.6f}")
    print(f"  K*_2 = {K_stars[1]:.6f}")
    print(f"  K*_2/K*_1 = {K_stars[1]/K_stars[0]:.2f}  (PDG: {R21_PDG:.2f})")
else:
    print(f"  Za malo K* ({len(K_stars)}) -- test Koidego niemozliwy.")
    print(f"  Zmian znaku g brak lub tylko K*_1 znane.")


# ─────────────────────────────────────────────────────────────────────────────
# PASS / FAIL
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PASS / FAIL")
print("=" * 72)

PASS_COUNT = 0
FAIL_COUNT = 0

def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok:
        PASS_COUNT += 1
    else:
        FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")


ok_P1 = len(sign_changes) >= 1
record("P1: g(K) zmienia znak w [0.015, 0.25]", ok_P1,
       f"{len(sign_changes)} zmian znaku" if ok_P1 else "brak zmian znaku")

# P2: sprawdz ciaglosc (brak przeskokow galezi w okolicach zmian znaku)
ok_P2 = False
if sign_changes:
    monotone = True
    for i, Ka, Kb, ga, gb in sign_changes:
        na = len(g_all_arr[i])
        nb = len(g_all_arr[i+1])
        if na != nb:
            monotone = False
    ok_P2 = monotone
record("P2: Zmiana znaku MONOTONICZNA (nie artefakt galezi)", ok_P2,
       "ta sama l. galezi po obu stronach" if ok_P2
       else "rozna l. galezi -- mozliwy artefakt")

ok_P3 = len(K_stars) >= 2
K_star2_found = K_stars[1] if len(K_stars) >= 2 else np.nan
record("P3: K*_2 wyznaczony przez brentq", ok_P3,
       f"K*_2 = {K_star2_found:.6f}" if ok_P3 else "K*_2 nieznaleziony")

ok_P4 = (ok_P3 and
         abs(K_star2_found / K1_REF - R21_PDG) / R21_PDG < 0.30)
if ok_P3:
    record("P4: K*_2/K*_1 blizej R21_PDG (<30%)", ok_P4,
           f"K*_2/K*_1 = {K_star2_found/K1_REF:.2f}  (PDG: {R21_PDG:.2f})")
else:
    record("P4: K*_2/K*_1 blizej R21_PDG (<30%)", False, "K*_2 nieznaleziony")

ok_P5 = len(K_stars) >= 3
if ok_P5:
    K1, K2, K3 = K_stars[:3]
    Q_val = koide_Q(K1, K2, K3)
    record("P5: K*_3 znaleziony + Q_TGP obliczone", ok_P5,
           f"Q_TGP = {Q_val:.5f}  (Koide = {2/3:.5f},"
           f" odch = {abs(Q_val-2/3)*100/(2/3):.1f}%)")
else:
    record("P5: K*_3 znaleziony + Q_TGP obliczone", False,
           f"tylko {len(K_stars)} K* znalezionych")

print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("PODSUMOWANIE  p77_K2_scan.py")
print("=" * 72)
print()
print(f"  alpha = {ALPHA_K_REF}, a_Gam = {A_GAM}, V_mod STANDARD")
print(f"  K*_1_ref = {K1_REF:.6f}")
print()
if ok_P1 and ok_P2:
    print("  K*_2 znaleziony (zmiana znaku monotoniczna) -- prawdziwy soliton")
    if not ok_P4:
        print(f"  JEDNAK K*_2/K*_1 = {K_star2_found/K1_REF:.2f} << R21_PDG = {R21_PDG:.2f}")
        print("  => Solitony V_mod standard NIE odtwarzaja hierarchii leptonow!")
        print("  => Path 1 wymaga modyfikacji V_mod (skala K*_2 jest zla).")
    else:
        print(f"  K*_2/K*_1 ≈ R21_PDG = {R21_PDG:.2f} -- masa muona odtworzona!")
elif ok_P1 and not ok_P2:
    print("  Zmiana znaku g jest artefaktem przełaczenia galezi, nie prawdziwym K*_2.")
    print("  Wniosek: V_mod STANDARD ma tylko 1 samospojny soliton (K*_1).")
    print("  Path 1 (Koide) wymaga modyfikacji V_mod -- p15.")
else:
    print("  Brak zmian znaku g w [0.015, 0.25] -- K*_2 nie istnieje w V_mod standard.")
    print("  Wyniki Sekcji D w p76 byly artefaktem (zbyt mala gestosc skanu).")
    print("  Path 1 (Koide) wymaga modyfikacji V_mod -- p15.")
print()
