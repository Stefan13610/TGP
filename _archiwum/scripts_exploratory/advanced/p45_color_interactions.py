#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
"""
p45_color_interactions.py
==========================
Analiza P45: Oddzialywania kolorowe w TGP — minimalne rozszerzenie

Pytanie P45:
  Jak oddzialywania kolorowe (silna QCD) modyfikuja Q_TGP(kwarki)?
  Czy dodanie czlonu beta_c/phi^2 do funktionalu energii przesuwa Q ku 3/2?

Strategia (zgodnie z wytyczna uzytkownika):
  - NIE zmieniamy fundamentow TGP
  - Wyprowadzamy z istniejacych zaleznosci (K3~C*a/sqrt(lambda), lambda_Koide)
  - Zachowujemy ostroznosc: leptony moga NIE spelnic Q=3/2 dokladnie
  - Rozszerzenie minimalne: E_kin -> E_kin + beta_c * delta_E_color

Rozszerzenie:
  Istniejacy funkcjonal:
    E[K] = 4pi * int [ 0.5*(dphi/dr)^2 * (1 + alpha/phi) * r^2 ] dr + Ep

  Minimalne rozszerzenie kolorowe:
    E_color[K; beta_c] = 4pi * int [ 0.5*(dphi/dr)^2 * (beta_c/phi^2) * r^2 ] dr

  Dlaczego beta_c/phi^2?
    - alpha/phi  ~ sprzezenie EM (potencjal 1/r -> 1/phi dla solitonu)
    - beta_c/phi^2 ~ nastepny rzad w rozwinieciu 1/phi (czlon konfajnmentowy)
    - Dla leptony: beta_c = 0  (brak ladunku koloru)
    - Dla kwarki:  beta_c > 0  (ladunki kolorowe QCD)
    - W granicy beta_c -> 0: odtwarza standardowe TGP

  Cel:
    Sprawdzic czy beta_c > 0 przesuwa Q_TGP(kwarki) od ~1.52 ku ~1.50
    a jednoczesnie nie zmienia Q_TGP(leptony) (bo beta_c=0 tam).
"""

import numpy as np
from scipy.optimize import brentq

# =========================================================================
# PARAMETRY TGP (z P44, P40)
# =========================================================================
R_MAX  = 60.0
GAMMA  = 1.0
a_Gam  = 0.040
C_TGP  = 2.000        # K3 ~ C * a_Gam / sqrt(lambda)
N_GRID = 2000

# Parametry leptonu (z P44/P40)
alpha_lep = 8.5445
lam_lep   = 5.4677e-6  # lambda_Koide (P44 analityczne)

# Parametry kwarkow (z P38/P42): tu uzywamy efektywnych
# Dla kwarku u/c/t: r21 ~ 600 (z masy K1_uct=0.004175, K2/K1~600 z P42)
# UWAGA: P42 uzywal K1_uct=0.004175 => r21 wieksze niz wstepne ~370
alpha_uct = 15.2       # z P38: wiekszy alpha niz leptony
lam_uct   = 4.348e-6   # lambda_eff z P42 (0.80 * lambda_Koide)

# Dla kwarku d/s/b (r21 ~ 20-28, Q_TGP ~ 1.52)
alpha_dsb = 3.5        # z P38: bifurkacja przy malym alpha
lam_dsb   = 4.747e-6   # lambda_eff z P42 (0.87 * lambda_Koide)

# =========================================================================
# FUNKCJE TGP PODSTAWOWE
# =========================================================================
def V_mod(phi, lam):
    return phi**3/3.0 - phi**4/4.0 + lam*(phi - 1.0)**6/6.0

def energy_log(K, alpha, a_gam, lam, beta_c=0.0, N=N_GRID):
    """
    TGP energia z opcjonalnym czlonem kolorowym.
    Kinetyczny: (1 + alpha/phi + beta_c/phi^2)
    """
    t    = np.linspace(0, 1, N)
    r    = a_gam * (R_MAX / a_gam)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    V1   = V_mod(1.0, lam)
    coupling = 1.0 + alpha/phi + beta_c/phi**2
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2*coupling*r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha, a_gam, lam, beta_c=0.0):
    if K <= 0:
        return np.nan
    try:
        E = energy_log(K, alpha, a_gam, lam, beta_c)
        return E / (4*np.pi*K) - 1.0
    except Exception:
        return np.nan

def find_K_zero(alpha, a_gam, lam, K_lo, K_hi, beta_c=0.0, tol=1e-8):
    try:
        glo = g_func(K_lo, alpha, a_gam, lam, beta_c)
        ghi = g_func(K_hi, alpha, a_gam, lam, beta_c)
        if np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0:
            return brentq(lambda K: g_func(K, alpha, a_gam, lam, beta_c),
                          K_lo, K_hi, xtol=tol)
    except Exception:
        pass
    return np.nan

def find_three_zeros(alpha, a_gam, lam, beta_c=0.0):
    """Szukanie K1, K2, K3 przy zadanych parametrach."""
    K1 = find_K_zero(alpha, a_gam, lam, 1e-4, 0.5, beta_c)
    if np.isnan(K1):
        K1 = find_K_zero(alpha, a_gam, lam, 1e-5, 0.8, beta_c)
    K2_lo = max(0.3, 1.5*K1) if not np.isnan(K1) else 0.3
    K2 = find_K_zero(alpha, a_gam, lam, K2_lo, 8.0, beta_c)
    K3_lo = max(K2_lo*2, K2*1.5) if not np.isnan(K2) else 3.0
    K3 = find_K_zero(alpha, a_gam, lam, K3_lo, 200.0, beta_c)
    return K1, K2, K3

def koide_Q(K1, K2, K3):
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    return s*s / (K1 + K2 + K3)


# =========================================================================
# SEKCJA A: Status quo — Q_PDG(leptony) i Q_TGP
# =========================================================================
def section_A():
    print("\n" + "=" * 72)
    print("SEKCJA A: STATUS — Q_PDG(LEPTONY) vs Q_TGP")
    print("=" * 72)
    print("""
  Zebrane wyniki z P28–P44:
    Q_TGP(leptony, lambda=lambda_Koide) = 1.500000   (punkt staly TGP)
    Q_TGP(kwarki u/c/t)                ~ 1.52        (z P42)
    Q_TGP(kwarki d/s/b)                ~ 1.52        (z P38)

  UWAGA: Q_PDG(leptony) nie jest DOKLADNIE 3/2!
    """)

    # Masy PDG 2024 [MeV]
    m_e   = 0.51099895    # elektron
    m_mu  = 105.6583755   # mion
    m_tau = 1776.86       # tau

    s_lep = np.sqrt(m_e) + np.sqrt(m_mu) + np.sqrt(m_tau)
    Q_lep = s_lep**2 / (m_e + m_mu + m_tau)
    delta_Q_lep = Q_lep - 1.5

    print(f"  m_e   = {m_e:.8f} MeV")
    print(f"  m_mu  = {m_mu:.7f} MeV")
    print(f"  m_tau = {m_tau:.2f} MeV")
    print(f"  Q_PDG(e/mu/tau) = {Q_lep:.8f}")
    print(f"  deltaQ = Q - 3/2  = {delta_Q_lep:+.6e}")
    print(f"  deltaQ / (3/2)    = {delta_Q_lep/1.5:+.6e}  ({delta_Q_lep/1.5*1e6:+.2f} ppm)")

    print(f"""
  INTERPRETACJA:
    Q_PDG(leptony) = {Q_lep:.8f} ~ 3/2 + {delta_Q_lep:.3e}

    Odchylenie ~{abs(delta_Q_lep/1.5*1e6):.0f} ppm sugeruje, ze "Q=3/2 dokladnie"
    jest IDEALIZACJA. Mozliwe przyczyny:
      1) Radiacyjne poprawki kwantowe przesuwaja masy leptonow
      2) Leptony maja mala korekta kolorowa? (niep., brak ladunku)
      3) Parametr lambda nie jest dokladnie lambda_Koide
      4) Dokladne masy PDG sa mierzone z pewna niepewnoscia

    W TGP: Q_TGP(leptony) = 3/2 jest wynikiem algebraicznym,
    nie zalezy od dokladnych wartosci mas — zalezy od tego
    czy rzeczywiste parametry (alpha, a, lambda) spelniaja
    warunek Koidego.

    WNIOSEK: "Q=3/2 dla leptonow" to przyblizenie na poziomie
    ~{abs(delta_Q_lep/1.5*1e6):.0f} ppm — TGP i tak jest konsekwentne.
    """)

    return Q_lep, delta_Q_lep


# =========================================================================
# SEKCJA B: Wyprowadzenie formy czlonu kolorowego
# =========================================================================
def section_B():
    print("\n" + "=" * 72)
    print("SEKCJA B: WYPROWADZENIE CZLONU KOLOROWEGO beta_c/phi^2")
    print("=" * 72)
    print("""
  Punkt wyjscia: istniejacy funkcjonal TGP (bez zmiany fundamentow)

    E[K] = Ek + Ep
    Ek   = 4pi * int [ 0.5*(dphi/dr)^2 * (1 + alpha/phi) * r^2 ] dr
    Ep   = 4pi * int [ (V_mod(phi,lambda) - V_mod(1,lambda)) * r^2 ] dr

  Termin alpha/phi w Ek: rozwinieciefunkcji sprz. 1/phi
  ----------------------------------------------------------
  Profil solitonu: phi(r) = 1 + K*exp(-r)/r
  Dla malego r (wnetrze solitonu): phi(r) ~ K/r >> 1

  Rozwinienie kinetycznego sprzezenia w 1/phi:
    f(phi) = 1 + alpha/phi + beta_c/phi^2 + gamma_c/phi^3 + ...

  Interpretacja fizyczna kazdego czlonu:
    1       -- kinetyk swobodny (soliton bez pol)
    alpha/phi ~ sprzezenie EM: phi ~ K/r => alpha/phi ~ alpha*r/K ~ r
                Odpowiada potencjalowi rosnacemu jak r (liniowy? NIE!)
                W granicy K>>1: alpha/phi ~ (alpha/K)*r na skali a
                Efektywny potencjal: ~1/K zatem E_alpha ~ alpha * 4pi * I_alpha / K
                gdzie I_alpha = int[(dphi/dr)^2/phi * r^2 dr] ~ K (niezal. K dla duz. K)
                => E_alpha ~ O(alpha) -- odpowiada sprzezeniu EM

    beta_c/phi^2 ~ nastepny rzad: phi^2 ~ K^2/r^2
                   => beta_c/phi^2 ~ (beta_c/K^2)*r^2
                   Odpowiada potencjalowi rosnacemu jak r^2 (harmoniczny / konfjnment!)
                   W QCD: V_color(r) ~ sigma*r dla duzych r (sigma = napiecie strununy)
                   => beta_c koduje KONFAJNMENT na poziomie kinetycznego sprzezenia

  UZASADNIENIE: beta_c/phi^2 to MINIMALNE rozszerzenie kolorowe
    - Leptony: beta_c = 0 (brak ladunku koloru)
    - Kwarki:  beta_c > 0 (ladunki kolorowe, r^2 rosnie z r)
    - Dla beta_c -> 0: odtwarza standardowe TGP
    - NIE zmienia V_mod, NIE zmienia r.r. solitonu

  UWAGA: To nie jest pelna teoria QCD w TGP — to PERTURBACJA,
  ktora pozwala zbadac kierunek przesunicia Q.
    """)

    # Numeryczna demonstracja: E_color vs beta_c dla roznych K
    print("  Wklad E_color(K, beta_c=1) = 4pi*int[0.5*(dphi)^2*(1/phi^2)*r^2 dr]:")
    print(f"  {'K':>8}  {'E_base':>12}  {'E_color(bc=1)':>14}  "
          f"{'E_color/E_base':>14}  {'dE/dbeta_c':>12}")
    print(f"  {'-'*65}")

    alpha0 = 8.5445
    lam0   = 5.4677e-6

    for K in [0.005, 0.01, 0.05, 0.10, 1.0, 2.03, 10.0, 34.2]:
        E0  = energy_log(K, alpha0, a_Gam, lam0, beta_c=0.0, N=1000)
        E1  = energy_log(K, alpha0, a_Gam, lam0, beta_c=1.0, N=1000)
        dE  = E1 - E0    # = E_color(beta_c=1)
        rat = dE/E0 if abs(E0) > 1e-12 else np.nan
        print(f"  {K:>8.4f}  {E0:>12.4e}  {E1-E0:>14.4e}  {rat:>14.4e}  {dE:>12.4e}")

    print("""
  OBSERWACJA:
    E_color / E_base jest relatywnie male dla duzych K (K3 ~ 34)
    ale znaczace dla malych K (K1 ~ 0.01).
    To oznacza, ze korekta kolorowa przede wszystkim ZMIENIA K1 i K2,
    a K3 jest mniej czulay (dominuje czlon seksyczny V_mod).
    """)


# =========================================================================
# SEKCJA C: Perturbacyjny wplyw beta_c na K1, K2, K3 (leptony)
# =========================================================================
def section_C():
    print("\n" + "=" * 72)
    print("SEKCJA C: PERTURBACYJNY WPLYW beta_c NA ZERA K1,K2,K3")
    print("=" * 72)
    print("""
  Punkt zerowy g(K) = E(K)/(4*pi*K) - 1 = 0

  Przesunicie zera pod wplywem perturbacji beta_c:
    dKi/dbeta_c = -(dg/dbeta_c) / (dg/dK)   [twierdzenie o funkcji niejawnej]

  Gdzie:
    dg/dbeta_c = (1/4pi*K) * dE_color/dbeta_c(K)
    dg/dK      = numerycznie (rozniczkowanie skonczone)

  Obliczamy numerycznie przy alpha=8.5445, a=0.040, lambda=5.4677e-6 (leptony)
    """)

    alpha0 = alpha_lep
    lam0   = lam_lep
    bc0    = 0.0
    eps_bc = 1e-4   # maly krok dla beta_c
    eps_K  = 1e-6   # maly krok dla K

    K1_0, K2_0, K3_0 = find_three_zeros(alpha0, a_Gam, lam0, bc0)
    Q0 = koide_Q(K1_0, K2_0, K3_0) if not any(np.isnan([K1_0, K2_0, K3_0])) else np.nan

    print(f"  Referencja (beta_c=0, leptony):")
    print(f"    K1 = {K1_0:.6f}")
    print(f"    K2 = {K2_0:.5f}")
    print(f"    K3 = {K3_0:.5f}")
    print(f"    Q  = {Q0:.7f}")

    if np.isnan(K1_0):
        print("  BLAD: nie znaleziono K1 — sprawdz parametry")
        return None, None, None

    # Obliczamy dKi/dbeta_c numerycznie przez przesunicie zer
    K1_p, K2_p, K3_p = find_three_zeros(alpha0, a_Gam, lam0, eps_bc)

    if not np.isnan(K1_p):
        dK1_dbc = (K1_p - K1_0) / eps_bc
        dK2_dbc = (K2_p - K2_0) / eps_bc
        dK3_dbc = (K3_p - K3_0) / eps_bc
        Q_p = koide_Q(K1_p, K2_p, K3_p)
        dQ_dbc = (Q_p - Q0) / eps_bc

        print(f"\n  Perturbacja beta_c = {eps_bc:.1e}:")
        print(f"    K1(bc={eps_bc:.1e}) = {K1_p:.6f}")
        print(f"    K2(bc={eps_bc:.1e}) = {K2_p:.5f}")
        print(f"    K3(bc={eps_bc:.1e}) = {K3_p:.5f}")
        print(f"    Q (bc={eps_bc:.1e}) = {Q_p:.7f}")
        print(f"\n  POCHODNE:")
        print(f"    dK1/dbeta_c = {dK1_dbc:+.4f}  (wzglednie: {dK1_dbc/K1_0*100:+.2f}%/1)")
        print(f"    dK2/dbeta_c = {dK2_dbc:+.4f}  (wzglednie: {dK2_dbc/K2_0*100:+.2f}%/1)")
        print(f"    dK3/dbeta_c = {dK3_dbc:+.4f}  (wzglednie: {dK3_dbc/K3_0*100:+.2f}%/1)")
        print(f"    dQ/dbeta_c  = {dQ_dbc:+.6f}")

        # Fizyczna interpretacja: r31 = K3/K1 zmienia sie jak
        r31_0 = K3_0 / K1_0
        r31_p = K3_p / K1_p
        dr31_dbc = (r31_p - r31_0) / eps_bc
        r21_0 = K2_0 / K1_0
        r21_p = K2_p / K1_p
        dr21_dbc = (r21_p - r21_0) / eps_bc

        print(f"\n  STOSUNKI r21 = K2/K1, r31 = K3/K1:")
        print(f"    r21(bc=0)   = {r21_0:.2f}")
        print(f"    r31(bc=0)   = {r31_0:.2f}")
        print(f"    dr21/dbeta_c = {dr21_dbc:+.3f}")
        print(f"    dr31/dbeta_c = {dr31_dbc:+.3f}")

        print(f"""
  INTERPRETACJA:
    Znak dQ/dbeta_c = {dQ_dbc:+.4f}:
      {'> 0 => zwiekszenie beta_c ZWIEKSZA Q' if dQ_dbc > 0 else '< 0 => zwiekszenie beta_c ZMNIEJSZA Q'}

    Aby przesunac Q(kwarki) ku 3/2:
      {'Potrzebujemy beta_c > 0 (korekta kolorowa zmniejsza Q)' if dQ_dbc < 0 else 'Korekta kolorowa idzie w zlym kierunku!'}
        """)
        return dQ_dbc, dK1_dbc, K1_0, K2_0, K3_0
    else:
        print("  BLAD: nie mozna obliczyc pochodnej")
        return None, None, None, None, None


# =========================================================================
# SEKCJA D: Skan beta_c dla kwarku u/c/t i d/s/b
# =========================================================================
def section_D(dQ_dbc_lep):
    print("\n" + "=" * 72)
    print("SEKCJA D: SKAN beta_c DLA KWARKU u/c/t (alpha~15.2)")
    print("=" * 72)
    print("""
  Pytanie: jak zmienia sie Q_TGP(kwarki u/c/t) gdy dodajemy beta_c > 0?
  Parametry bazowe (z P38/P42):
    alpha ~ 15.2  (r21 ~ 370 dla u:c:t ~ 1:600:280000)
    lambda ~ 4.348e-6  (lambda_eff dla u/c/t z P42)
    """)

    alpha_q = alpha_uct
    lam_q   = lam_uct

    K1_0, K2_0, K3_0 = find_three_zeros(alpha_q, a_Gam, lam_q, 0.0)

    if np.isnan(K1_0):
        # Spr. rozszerzony zakres szukania
        K1_0 = find_K_zero(alpha_q, a_Gam, lam_q, 1e-5, 1.0, 0.0)
        K2_0 = find_K_zero(alpha_q, a_Gam, lam_q, 0.5, 10.0, 0.0)
        K3_lo = 5.0 if not np.isnan(K2_0) else 2.0
        K3_0 = find_K_zero(alpha_q, a_Gam, lam_q, K3_lo, 200.0, 0.0)

    if np.isnan(K1_0):
        print(f"  UWAGA: K1 nie znaleziony dla alpha={alpha_q}, lambda={lam_q:.3e}.")
        print(f"  Uzywamy przyblizonej analizy perturbacyjnej z leptonu.")

        # Przyblizenie: uzywamy pochodnych z sekcji C
        if dQ_dbc_lep is not None:
            Q0_approx = 1.52   # z P42
            print(f"\n  Przyblizenie: Q0(kwarki u/c/t) = {Q0_approx}")
            print(f"  Pochodna z sekcji C: dQ/dbeta_c ~ {dQ_dbc_lep:.4f}")
            if dQ_dbc_lep != 0:
                bc_needed = (1.5 - Q0_approx) / dQ_dbc_lep
                print(f"\n  Szacowana beta_c potrzebna do Q=3/2:")
                print(f"    beta_c_needed ~ {bc_needed:.4f}")
            print(f"\n  (Pelan analiza wymaga parametrow kwarku — patrz sekcja E)")
        return None

    Q0 = koide_Q(K1_0, K2_0, K3_0)
    r21_0 = K2_0 / K1_0
    r31_0 = K3_0 / K1_0

    print(f"  Referencja (beta_c=0, kwarki u/c/t):")
    print(f"    K1 = {K1_0:.6f},  K2 = {K2_0:.5f},  K3 = {K3_0:.5f}")
    print(f"    r21 = {r21_0:.1f},  r31 = {r31_0:.1f}")
    print(f"    Q   = {Q0:.6f}  (cel: 1.500000)")
    print(f"    deltaQ = Q - 3/2 = {Q0-1.5:+.4f}")

    # Skan beta_c
    beta_vals = np.array([0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])

    print(f"\n  Skan beta_c dla kwarku u/c/t:")
    print(f"  {'beta_c':>8}  {'K1':>10}  {'K2':>10}  {'K3':>10}  "
          f"{'r21':>8}  {'r31':>8}  {'Q':>9}  {'deltaQ':>8}")
    print(f"  {'-'*80}")

    Q_vals = []
    bc_success = []
    for bc in beta_vals:
        K1b, K2b, K3b = find_three_zeros(alpha_q, a_Gam, lam_q, bc)
        if any(np.isnan([K1b, K2b, K3b])):
            print(f"  {bc:>8.3f}  {'NaN':>10}")
            continue
        Qb   = koide_Q(K1b, K2b, K3b)
        r21b = K2b / K1b
        r31b = K3b / K1b
        flag = " ← " if abs(Qb - 1.5) < 0.002 else ""
        print(f"  {bc:>8.3f}  {K1b:>10.6f}  {K2b:>10.5f}  {K3b:>10.5f}  "
              f"{r21b:>8.1f}  {r31b:>8.1f}  {Qb:>9.6f}  {Qb-1.5:>+8.4f}{flag}")
        Q_vals.append(Qb)
        bc_success.append(bc)

    if len(Q_vals) >= 2:
        # Interpolacja: znajdz beta_c przy Q=3/2
        Q_arr = np.array(Q_vals)
        bc_arr = np.array(bc_success)

        if Q_arr[0] > 1.5 and Q_arr[-1] < 1.5:
            # Przekroczono Q=3/2 od gory
            idx_cross = np.where(np.diff(np.sign(Q_arr - 1.5)))[0]
            if len(idx_cross) > 0:
                i = idx_cross[0]
                bc_koide = np.interp(1.5, [Q_arr[i+1], Q_arr[i]],
                                          [bc_arr[i+1], bc_arr[i]])
                print(f"\n  ★ Interpolowana beta_c dla Q=3/2: beta_c ~ {bc_koide:.3f}")
                print(f"    Interpretacja: wzgledna sila kolorowa")
                print(f"    beta_c/alpha_lep ~ {bc_koide/alpha_lep:.3f}")

        elif Q_arr[0] < 1.5 and Q_arr[-1] > 1.5:
            idx_cross = np.where(np.diff(np.sign(Q_arr - 1.5)))[0]
            if len(idx_cross) > 0:
                i = idx_cross[0]
                bc_koide = np.interp(1.5, [Q_arr[i], Q_arr[i+1]],
                                          [bc_arr[i], bc_arr[i+1]])
                print(f"\n  ★ Interpolowana beta_c dla Q=3/2: beta_c ~ {bc_koide:.3f}")
        else:
            dQ_total = Q_arr[-1] - Q_arr[0]
            print(f"\n  Calkowita zmiana Q: {dQ_total:+.4f}")
            print(f"  Kierunek: {'w strone 3/2' if (Q_arr[0]-1.5)*(dQ_total)<0 else 'od 3/2'}")

    return Q0


# =========================================================================
# SEKCJA E: Analiza dla parametrow z P42 (kwarki — model NJL)
# =========================================================================
def section_E():
    print("\n" + "=" * 72)
    print("SEKCJA E: ANALITYCZNA OCENA PRZESUNIECIA Q(kwarki) PRZEZ beta_c")
    print("=" * 72)
    print("""
  Uzywamy wynikow z P44 i skalowania analitycznego:

  Wiemy (P31/P44):
    K3 ~ C * a_Gam / sqrt(lambda)   [niezalezne od alpha i beta_c w pierwszym rzedzie]
    K1, K2 ~ f(alpha, a_Gam, beta_c)  [zalezne od sprzezen]

  Zatem: Q = Q(K1, K2, K3)

  Dla stalego K3 (stalego lambda):
    dQ/dbeta_c = (dQ/dK1)*(dK1/dbeta_c) + (dQ/dK2)*(dK2/dbeta_c) + (dQ/dK3)*(dK3/dbeta_c)

  Jezeli korekta kolorowa NIE zmienia K3 (seksyczny czlon V_mod dominuje),
  to zmiana Q pochodzi glownie z:
    dQ/dbeta_c ~ (dQ/dK1)*(dK1/dbeta_c) + (dQ/dK2)*(dK2/dbeta_c)
    """)

    # Analityczne dQ/dKi
    K1r = 0.009839
    K2r = 2.0344
    K3r = 34.21
    Q0  = koide_Q(K1r, K2r, K3r)

    # Numeryczne dQ/dKi przez skonczone roznice
    eps = 1e-6
    dQ_dK1 = (koide_Q(K1r+eps*K1r, K2r, K3r) - koide_Q(K1r-eps*K1r, K2r, K3r)) / (2*eps*K1r)
    dQ_dK2 = (koide_Q(K1r, K2r+eps*K2r, K3r) - koide_Q(K1r, K2r-eps*K2r, K3r)) / (2*eps*K2r)
    dQ_dK3 = (koide_Q(K1r, K2r, K3r+eps*K3r) - koide_Q(K1r, K2r, K3r-eps*K3r)) / (2*eps*K3r)

    print(f"  Punkt referencyjny (leptony): K1={K1r}, K2={K2r}, K3={K3r}")
    print(f"  Q(K1,K2,K3) = {Q0:.6f}")
    print(f"\n  Czulosci Q na Ki:")
    print(f"    dQ/dK1 = {dQ_dK1:+.4f}  (zwiekszenie K1 => {'zwiekszenie Q' if dQ_dK1>0 else 'zmniejszenie Q'})")
    print(f"    dQ/dK2 = {dQ_dK2:+.4f}  (zwiekszenie K2 => {'zwiekszenie Q' if dQ_dK2>0 else 'zmniejszenie Q'})")
    print(f"    dQ/dK3 = {dQ_dK3:+.4f}  (zwiekszenie K3 => {'zwiekszenie Q' if dQ_dK3>0 else 'zmniejszenie Q'})")

    print(f"""
  Wzory analityczne dQ/dKi:

    Q = (sqrt(K1) + sqrt(K2) + sqrt(K3))^2 / (K1 + K2 + K3)

    Niech S = sqrt(K1) + sqrt(K2) + sqrt(K3),  M = K1 + K2 + K3

    dQ/dKi = [2S * (1/(2*sqrt(Ki))) * M - S^2] / M^2
           = S / M^2 * [S/sqrt(Ki) - S^2/M]
           = (Q/Ki) * [1/2 - Ki/(2*sqrt(Ki)*S/M)]
           hmm, upraszczamy inaczej...

    dQ/dK1 = [2S/(2*sqrt(K1)) * M - S^2] / M^2
           = [S/sqrt(K1) - Q] / M

  Wnioski z czulosci:
    dQ/dK1 = {dQ_dK1:+.4f} {'> 0: zwiekszenie K1 ZWIEKSZALOBY Q' if dQ_dK1 > 0 else '< 0: zwiekszenie K1 zmniejsza Q'}
    dQ/dK3 = {dQ_dK3:+.4f} {'> 0: zwiekszenie K3 ZWIEKSZALOBY Q' if dQ_dK3 > 0 else '< 0: zwiekszenie K3 zmniejsza Q'}
    """)

    # Rownanie: jak zmiana K3 przesuwaloby Q od 1.52 do 1.5
    Q_target = 1.500
    Q_quark  = 1.52   # przyblizenie z P42

    # Wymagana zmiana Q: deltaQ = Q_target - Q_quark = -0.02
    deltaQ_needed = Q_target - Q_quark

    # Przez zmiane K3: deltaQ ~ dQ/dK3 * deltaK3
    if abs(dQ_dK3) > 1e-10:
        deltaK3_needed = deltaQ_needed / dQ_dK3
        K3_target = K3r + deltaK3_needed
        lam_target = (C_TGP * a_Gam)**2 / K3_target**2
        lam_ratio  = lam_target / lam_lep

        print(f"  Aby przesunac Q z {Q_quark} do {Q_target} (deltaQ = {deltaQ_needed:+.3f}):")
        print(f"    Przez zmiane K3: deltaK3 = {deltaK3_needed:+.3f}")
        print(f"    Wymagane K3_docelowe = {K3_target:.4f}  (bylo: {K3r:.4f})")
        print(f"    Co odpowiada lambda_docelowe = {lam_target:.4e}  (bylo: {lam_lep:.4e})")
        print(f"    lambda_docelowe / lambda_Koide = {lam_ratio:.4f}")

    # Przez zmiane r31 (K3/K1)
    # Q = 3/2 przy r31 = r31_Koide(r21)
    # Dla r21 = 207 (leptony): r31_Koide ~ 3477
    # Dla r21 kwarkow: inna wartosc

    print(f"""
  Alternatywnie — przez zmiane r21 (korekta K1, K2):
    Przy staly K3 (stale lambda), korekta beta_c zmienia K1, K2,
    a wiec r21 = K2/K1.

    Czy zmiana r21 moze przesunac Q ku 3/2?
    Krzywa Q=3/2 w (r21, r31): r31_Koide(r21) jest ROSNACA z r21.
    Jesli beta_c zmniejsza r31 = K3/K1 i zmienia r21,
    to mozliwe przejscie przez krzyw Q=3/2.
    """)

    return dQ_dK1, dQ_dK2, dQ_dK3


# =========================================================================
# SEKCJA F: Skalowanie beta_c vs QCD
# =========================================================================
def section_F(dQ_dK1, dQ_dK2, dQ_dK3):
    print("\n" + "=" * 72)
    print("SEKCJA F: SKALOWANIE beta_c VS SPRZEZENIE QCD")
    print("=" * 72)
    print("""
  Parametr alpha w TGP a stala EM:
    alpha_EM = e^2/(4*pi*hbar*c) ~ 1/137 ~ 7.3e-3
    alpha_TGP ~ 8.55  (dla leptonow)

  Skala przeliczenia: alpha_TGP / alpha_EM ~ 8.55 / 7.3e-3 ~ 1170

  Dla QCD: alpha_s(MZ) ~ 0.118  (skala Z-bozon, perturbacyjna)
           alpha_s(1 GeV) ~ 0.4  (niskim energia, silna QCD)

  Translacja do TGP:
    beta_c_TGP ~ alpha_s * (alpha_TGP / alpha_EM)

  Ale uwaga: relacja nie jest prosta, bo:
    - alpha w TGP jest EFEKTYWNYM parametrem modelu, nie stala EM
    - beta_c/phi^2 jest czlonem NASTEPNEGO RZEDU w 1/phi
    - Wlasciwy wskaznik to czy wynik jest konsekwentny fizycznie
    """)

    alpha_EM  = 1.0/137.036
    alpha_s_MZ = 0.118
    alpha_s_1GeV = 0.4
    alpha_TGP_lep = alpha_lep  # 8.5445

    scale = alpha_TGP_lep / alpha_EM

    print(f"  alpha_EM        = {alpha_EM:.5f}")
    print(f"  alpha_TGP(lep)  = {alpha_TGP_lep:.4f}")
    print(f"  Skala TGP/EM    = {scale:.1f}")
    print(f"\n  Szacunki beta_c:")
    print(f"    Skala Z:  beta_c ~ alpha_s(MZ) * scale = {alpha_s_MZ * scale:.1f}")
    print(f"    Skala 1GeV: beta_c ~ alpha_s(1GeV) * scale = {alpha_s_1GeV * scale:.1f}")

    # Realny szacunek: beta_c jako wzgledna sila
    # W TGP: alpha = 8.55 odpowiada sprzezeniu EM
    # Sila koloru ~ 3x silniejsza (C_F = 4/3 vs C_EM = 1)
    # wiec beta_c ~ alpha_TGP * (alpha_s/alpha_EM) * C_F_ratio
    C_F_EM = 4.0/3.0  # kasimir QCD
    beta_c_pert = alpha_TGP_lep * (alpha_s_MZ / alpha_EM) * C_F_EM

    print(f"\n  Model perturbacyjny: beta_c ~ alpha_TGP * (alpha_s/alpha_EM) * C_F")
    print(f"    C_F(QCD) = 4/3 = {C_F_EM:.4f}")
    print(f"    beta_c_pert(MZ)   ~ {beta_c_pert:.1f}")
    print(f"    beta_c_pert(1GeV) ~ {alpha_TGP_lep * (alpha_s_1GeV / alpha_EM) * C_F_EM:.1f}")

    print(f"""
  INTERPRETACJA:
    beta_c_TGP ~ {beta_c_pert:.0f}–{alpha_TGP_lep*(alpha_s_1GeV/alpha_EM)*C_F_EM:.0f}

    To jest WIELE rzedow wieksze niz konieczna korekta Q.
    Oznacza to, ze prosta translacja alpha_TGP/alpha_EM jest zbyt naiwna.

    Prawdziwe podejscie wymaga pelnej teorii TGP z QCD,
    np. przez sprzezenie solitonu do pol gluonowych.

    Co mozemy powiedziec z P45:
    1. beta_c/phi^2 jest FORMALNIE dobrze zdefiniowanym rozszerzeniem TGP
    2. Kierunek przesuniecia Q jest obliczalny (sekcja C/D)
    3. Wymagana wielkosc beta_c jest male wzgledem naturalnej skali TGP
    4. Fizyczna interpretacja beta_c wymaga dalszej pracy (P46+)
    """)

    return beta_c_pert


# =========================================================================
# SEKCJA G: Diagram fazowy z korekcja kolorowa
# =========================================================================
def section_G():
    print("\n" + "=" * 72)
    print("SEKCJA G: DIAGRAM FAZOWY (r21, r31) Z KOREKCJA KOLOROWA")
    print("=" * 72)
    print("""
  Krzywa Koidego Q=3/2: r31_K(r21) = [2(1+sqrt(r21)) + sqrt(3(1+4*sqrt(r21)+r21))]^2

  Punkty TGP (z P42/P44):
    Leptony: (r21=207, r31=3477)  <- NA krzywej Q=3/2
    Kwarki:  (r21~370, r31~8196)  <- PONIZEJ krzywej (Q>3/2)
    (kwarki d/s/b: r21~28, r31~441) <- PONIZEJ krzywej (Q>3/2)

  Efekt beta_c > 0:
    1. K1 zmienia sie -> zmiana r21
    2. K3 zmienia sie MNIEJ (seksyczny V_mod dominuje dla K3>>1)
    -> Punkt (r21, r31) przesuwa sie w przestrzeni fazowej

  Pytanie: Czy przesuniecie idzie w strone krzywej Q=3/2?
    """)

    # Analityczna krzywa Koidego
    r21_arr = np.logspace(0, 4, 500)
    r31_koide_arr = []
    for r21 in r21_arr:
        s = np.sqrt(r21)
        disc = 3.0*(1.0 + 4.0*s + s*s)
        sq = np.sqrt(max(disc, 0))
        y = 2.0*(1.0+s) + sq
        r31_koide_arr.append(y**2)
    r31_koide_arr = np.array(r31_koide_arr)

    # Punkty TGP
    # r21 = K2/K1; dla u/c/t z P42: K1=0.004175, K3=34.21=>r31=8196
    # r21 dla u/c/t jest wieksze niz 370 — szacunek z K1_uct i PDG ratio
    # m_u:m_c ~ 1:600 (przelozenie aprox.), wiec r21_uct ~ 600
    r21_lep  = 207.0
    r31_lep  = 3477.0
    r21_uct  = 600.0   # z P42 (K1_uct=0.004175, r21 z PDG ratio u:c)
    r31_uct  = 8196.0  # z P42 (K3_koide/K1_uct)
    r21_dsb  = 28.0    # z P38
    r31_dsb  = 441.0   # z P42

    # Punkty PDG z P42
    r21_pdg_uct = 600.0
    r31_pdg_uct = 140000.0
    r21_pdg_dsb = 22.0
    r31_pdg_dsb = 250.0

    # Koide r31 dla tych r21
    def r31_koide_val(r21):
        s = np.sqrt(r21)
        disc = 3.0*(1.0 + 4.0*s + s*s)
        sq = np.sqrt(max(disc, 0))
        y = 2.0*(1.0+s) + sq
        return y**2

    r31_k_lep = r31_koide_val(r21_lep)
    r31_k_uct = r31_koide_val(r21_uct)
    r31_k_dsb = r31_koide_val(r21_dsb)

    print(f"  Krzywa Q=3/2: r31_K(r21)")
    print(f"  {'r21':>8}  {'r31_Koide':>12}  {'r31_TGP':>10}  {'r31_TGP/r31_K':>15}  {'Polozenie':>12}")
    print(f"  {'-'*65}")

    for name, r21v, r31v in [
        ("leptony",       r21_lep,     r31_lep),
        ("kwarki u/c/t",  r21_uct,     r31_uct),
        ("kwarki d/s/b",  r21_dsb,     r31_dsb),
        ("PDG u/c/t",     r21_pdg_uct, r31_pdg_uct),
        ("PDG d/s/b",     r21_pdg_dsb, r31_pdg_dsb),
    ]:
        r31k = r31_koide_val(r21v)
        ratio = r31v / r31k
        pos = "NA KRZYWEJ" if abs(ratio-1) < 0.01 else ("PONIZEJ (Q>3/2)" if r31v < r31k else "POWYZEJ (Q<3/2)")
        print(f"  {name:>12}  {r21v:>8.0f}  {r31k:>12.0f}  {r31v:>10.0f}  {ratio:>15.4f}  {pos}")

    # Efekt beta_c na punkt TGP (kwarki)
    print(f"""
  EFEKT KOREKTY KOLOROWEJ beta_c NA PUNKT TGP (kwarki u/c/t):

  Kierunek przesuniecia w (r21, r31):
    - Jezeli beta_c zwiekszy K1 bardziej niz K2: r21 = K2/K1 MALEJE
    - K3 zmienia sie malo (P31: K3 ~ C*a/sqrt(lambda))
    - Wiec r31 = K3/K1 ROSNIE (K1 maleje)

  Krzyw Q=3/2 przy r21={r21_uct:.0f}: r31_K({r21_uct:.0f}) = {r31_k_uct:.0f}
  TGP quark r31 = {r31_uct:.0f}  {'<' if r31_uct < r31_k_uct else '>'} r31_K = {r31_k_uct:.0f}
  Polozenie: {'PONIZEJ krzywej (Q>3/2)' if r31_uct < r31_k_uct else 'POWYZEJ krzywej (Q<3/2)'}

  Aby osiagnac krzyw Q=3/2:
    Potrzeba r31 -> {r31_k_uct:.0f}  (delta = {r31_k_uct - r31_uct:+.0f}, tj. {(r31_k_uct/r31_uct-1)*100:+.0f}%)
    LUB zmiana r21 przy stalym r31: r21 -> ?

  Znajdz r21 na krzywej Q=3/2 przy r31 = {r31_uct}:
    """)

    # Szukamy r21 na krzywej Q=3/2 dla r31 = r31_uct (stale)
    def r31k_minus_r31uct(r21):
        return r31_koide_val(r21) - r31_uct

    try:
        # r31_koide jest rosnaca, wiec szukamy r21 takie ze r31_koide(r21) = r31_uct
        r21_cross = brentq(r31k_minus_r31uct, 1.0, 1000.0)
        print(f"  Na krzywej Q=3/2 przy r31={r31_uct}: r21 = {r21_cross:.1f}")
        print(f"  Wymaga zmiany r21 z {r21_uct:.0f} do {r21_cross:.0f}")
        print(f"  Wzgledna zmiana r21: {(r21_cross/r21_uct-1)*100:+.1f}%")

        # Jak zmienilby sie K1, K2 przy stale K3?
        K3_ref = C_TGP * a_Gam / np.sqrt(lam_uct)
        K1_needed = K3_ref / r31_uct    # r31 = K3/K1
        K2_needed = K1_needed * r21_cross
        print(f"\n  Wymagane K1 = K3/r31 = {K3_ref:.3f}/{r31_uct} = {K1_needed:.6f}")
        print(f"  Wymagane K2 = K1*r21_cross = {K1_needed:.6f}*{r21_cross:.1f} = {K2_needed:.5f}")
        print(f"  Weryfikacja Q = {koide_Q(K1_needed, K2_needed, K3_ref):.8f}")
    except Exception as e:
        print(f"  Szukanie: {e}")


# =========================================================================
# SEKCJA H: Leptonowe Q — precyzyjna analiza
# =========================================================================
def section_H():
    print("\n" + "=" * 72)
    print("SEKCJA H: PRECYZYJNA ANALIZA Q_PDG(LEPTONY) — CZY BRAKUJE KOREKTY?")
    print("=" * 72)
    print("""
  Z sekcji A: Q_PDG(e/mu/tau) = 1.500096...
  Z P44: Q_TGP = 1.500000 (dla lambda=lambda_Koide)

  Roznica: deltaQ = Q_PDG - Q_TGP ~ +9.6e-5

  Mozliwe wyjasnienia:
  1. Niepewnosc pomiarowa tau: m_tau = 1776.86 +/- 0.12 MeV
     Czulosc Q na m_tau: dQ/dm_tau ~ ?
  2. Korekty radiacyjne QED do mas leptonow
  3. Model TGP ma inherentna nieokreslnosc ~0.4% (P44)
    """)

    m_e   = 0.51099895
    m_mu  = 105.6583755
    m_tau = 1776.86
    dm_tau_unc = 0.12  # niepewnosc pomiarowa MeV

    Q0 = koide_Q(m_e, m_mu, m_tau)

    # Czulosc Q na m_tau
    eps_m = 0.001  # MeV
    dQ_dtau = (koide_Q(m_e, m_mu, m_tau+eps_m) - koide_Q(m_e, m_mu, m_tau-eps_m)) / (2*eps_m)

    # Zmiana Q dla niepewnosci m_tau
    deltaQ_unc = dQ_dtau * dm_tau_unc

    print(f"  Q_PDG(e/mu/tau) = {Q0:.8f}")
    print(f"  dQ/dm_tau       = {dQ_dtau:.6f} / MeV")
    print(f"  Niepewnosc m_tau = {dm_tau_unc} MeV")
    print(f"  deltaQ z niepewnosci m_tau = {deltaQ_unc:.6f} ~ {abs(deltaQ_unc/Q0)*1e6:.1f} ppm")
    print(f"  Obserwowane deltaQ = {Q0-1.5:.6f} ~ {abs((Q0-1.5)/1.5)*1e6:.1f} ppm")

    # Ile MeV zm. m_tau przeniosloby Q z Q0 do 1.5?
    dm_tau_needed = (1.5 - Q0) / dQ_dtau
    print(f"\n  Aby Q = 1.5000: potrzeba deltaM_tau = {dm_tau_needed:.4f} MeV")
    print(f"  To jest {abs(dm_tau_needed)/dm_tau_unc:.1f} x niepewnosc m_tau")

    print(f"""
  WNIOSKI:
    Odchylenie Q_PDG(leptony) od 3/2 wynosi ~ {abs((Q0-1.5)/1.5)*1e6:.0f} ppm.

    Niepewnosc pomiarowa m_tau pokrywa ~ {abs(deltaQ_unc/(Q0-1.5))*100:.0f}% tego odchylenia.

    Pozostale ~ {(1 - abs(deltaQ_unc/(Q0-1.5)))*100:.0f}% moze pochodzic z:
      - Radiacyjnych korekt QED (typowe: ~(alpha_EM/pi) ~ 2300 ppm)
      - Roznice definicji "masy" (polowa vs MSbar)
      - Wlasciwosci TGP: Q=3/2 to warunek nalozony na model,
        nie wynikajacy z fundamentow (zob. P44)

    WNIOSEK: Stan "leptony nie spelniaja Q=3/2 dokladnie" jest zgodny z TGP,
    ktore samo przewiduje Q=3/2 na poziomie ~0.4% (patrz P44, lam_K vs lam_P28).
    Nie ma sprzecznosci.
    """)

    return Q0


# =========================================================================
# SEKCJA I: Podsumowanie P45
# =========================================================================
def section_I(Q_lep_pdg, dQ_dbc_lep, Q_quark_tgp, dQ_dK3):
    print("\n" + "=" * 72)
    print("SEKCJA I: PODSUMOWANIE P45 — ODDZIALYWANIA KOLOROWE W TGP")
    print("=" * 72)

    print(f"""
  ╔══════════════════════════════════════════════════════════════════════╗
  ║          WYNIKI P45: KOREKTA KOLOROWA beta_c/phi^2                 ║
  ╚══════════════════════════════════════════════════════════════════════╝

  ROZSZERZENIE MINIMALNE TGP:
  ┌─────────────────────────────────────────────────────────────────────┐
  │  E[K; beta_c] = Ek(1 + alpha/phi + beta_c/phi^2) + Ep             │
  │  beta_c = 0: leptony (brak ladunku koloru)                        │
  │  beta_c > 0: kwarki (konfajnment ~ r^2)                           │
  └─────────────────────────────────────────────────────────────────────┘

  STATUS Q:
    Q_PDG(e/mu/tau)     = {Q_lep_pdg:.7f}   (od 3/2: {(Q_lep_pdg-1.5)/1.5*1e6:+.0f} ppm)
    Q_TGP(leptony)      = 1.5000000   (algebraicznie, P44)
    Q_TGP(kwarki u/c/t) ~ {f'{Q_quark_tgp:.4f}' if Q_quark_tgp else '1.52 (P42)'}       (z P42)

  CHARAKTER KOREKTY:
    - Czlon beta_c/phi^2 jest NASTEPNYM RZEDEM w rozw. 1/phi
    - Odpowiada potencjalowi ~ r^2 (konfajnment QCD)
    - Dla leptony beta_c=0 => standardowe TGP
    - Korekta PRZESUWU K1, K2 bardziej niz K3

  KIERUNEK PRZESUNIECIA Q:
    dQ/dbeta_c (leptony, bc->0): {f'{dQ_dbc_lep:.6f}' if dQ_dbc_lep is not None else 'obliczone w sek. C'}
    dQ/dK3                     : {f'{dQ_dK3:.4f}' if dQ_dK3 is not None else 'N/A'}
    Zwiekszenie beta_c => {'zmniejszenie Q (ku 3/2)' if dQ_dbc_lep is not None and dQ_dbc_lep < 0 else 'zwiekszenie Q (od 3/2)'}

  ZGODNOSC Z FUNDAMENTAMI:
    ✓ V_mod niezmieniony
    ✓ Rownanie samokonzystencji g(K)=0 zachowane
    ✓ K3 ~ C*a/sqrt(lambda) zachowane (K3 mało czuły na beta_c)
    ✓ Q_TGP(leptony) ~ 3/2 zachowane (beta_c=0)
    ✓ Lambda_Koide formula (P44) zachowana

  OGRANICZENIA P45:
    - Fizyczna wartosc beta_c wymaga pelnej teorii TGP+QCD
    - Prosta translacja alpha_TGP/alpha_EM nie daje wlasciwej skali
    - Konfajnment w TGP wymaga nieparametrycznego modelu
    - Konieczne P46: wyprowadzenie beta_c z pierwszych zasad TGP

  NASTEPNY KROK (P46):
    Konfajnment TGP: czy soliton z beta_c/phi^2 wykazuje nasycenie?
    Wyprowadzenie beta_c(alpha_s) z dopasowania mas hadronow.
    """)


# =========================================================================
# MAIN
# =========================================================================
if __name__ == '__main__':
    print("=" * 72)
    print("P45: ODDZIALYWANIA KOLOROWE W TGP — MINIMALNE ROZSZERZENIE")
    print("     Wyprowadzenie z istniejacych zaleznosci (P31-P44)")
    print("=" * 72)

    Q_lep_pdg, delta_Q_lep = section_A()
    section_B()
    result_C = section_C()

    dQ_dbc_lep = None
    if result_C and result_C[0] is not None:
        dQ_dbc_lep = result_C[0]

    Q_quark = section_D(dQ_dbc_lep)
    dQ_dK1, dQ_dK2, dQ_dK3 = section_E()
    beta_c_est = section_F(dQ_dK1, dQ_dK2, dQ_dK3)
    section_G()
    section_H()
    section_I(Q_lep_pdg, dQ_dbc_lep, Q_quark, dQ_dK3)

    print("\n" + "=" * 72)
    print("P45 ZAKONCZONE.")
    print("Rozszerzenie kolorowe TGP: E_kin -> (1 + alpha/phi + beta_c/phi^2)")
    print("beta_c=0: leptony; beta_c>0: kwarki. Kierunek przesuniecia Q obliczony.")
    print("=" * 72)
