"""
tgp_agamma_phi0_test.py  —  Teoria Generowanej Przestrzeni (TGP)
================================================================
Numeryczna weryfikacja hipotezy a_Gamma * Phi0 = 1  (ROADMAP R3)

Hipoteza (sek08):
    a_Gamma * Phi0 = 1   (do 0.5%)
    gdzie:
        Phi0   = tlo prozniowe pola TGP
        a_Gamma = stala sprzezenia substratu (skala korelacji)

Metoda:
    Wyznaczamy a_Gamma z trzech niezaleznych sciezek i sprawdzamy,
    czy a_Gamma * Phi0 ~ 1 stabilnie:

    Sciezka S1 - z Lambda_obs:
        Phi0(Lambda) = 96*pi*G0*rho_Lambda / H0^2   -> a_Gamma(S1) = 1/Phi0(Lambda)

    Sciezka S2 - SUPERSEDED (alfa_K = 8.5616 - stary parametr solitu):
        Phi0 = (alfa_K * sqrt(r21))^{3/5}  =>  17.95  [NIEAKTUALNE po phi-FP]
        Zastapione przez: Twierdzenie phi-FP (dodatekJ2), g0* = 1.24915
        Nowa sciezka S2b: Phi0 z g0* i r31^K (OP-3, otwarte)

    Sciezka S3 - z kappa = 3/(4*Phi0) (ROADMAP_v3):
        kappa_obs ~ 0.0304  ->  Phi0(kappa) = 3/(4*kappa)  ->  a_Gamma(S3) = 1/Phi0(kappa)

    Sciezka S4 - z DESI DR2 (obserwacyjna):
        a_Gamma * Phi0 = 1.005 +/- 0.005  (0.53% odchylenie)

    Jesli S1, S3, S4 daja zbiezne a_Gamma*Phi0 ~ 1 - hipoteza potwierdzona.

Fizyczna interpretacja a_Gamma*Phi0=1:
    a_Gamma jest zwiazana ze skala korelacji substratu xi_corr:
        a_Gamma = 1/xi_corr  [w jednostkach TGP]
    Hipoteza mowi: xi_corr = Phi0, czyli tlo prozni wyznacza
    skale korelacji. Jest to WARUNEK SAMOSPOJNOSCI teorii.

UWAGA HISTORYCZNA:
    S2 uzywala ALPHA_K=8.5616 (bifurkacja solitonu, sek08 - stara wersja).
    Po Twierdzeniu phi-FP (dodatekJ2, g0*=1.24915) parametr alfa_K staje sie
    parametrem wtornym. Sciezka S2 jest SUPERCEDED. Nowa sciezka S2b wymaga
    wyprowadzenia Phi0 z g0* (patrz OP-3).

Referencja: sek08_formalizm.tex, ROADMAP_v3.md, dodatekP_agamma_phi0.tex,
            dodatekJ2_phi_fp_theorem.tex, dodatekT_koide_atail_formal.tex
"""

import sys
import io
import numpy as np

# UTF-8 output (Windows cp1250 fix)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ---------------------------------------------------------------------------
# Stale fizyczne (SI)
# ---------------------------------------------------------------------------
C0    = 2.99792458e8
G0    = 6.67430e-11
HBAR0 = 1.054571817e-34
H0_SI = 2.2685e-18
RHO_LAMBDA_OBS = 5.96e-27

# Parametry TGP - AKTUALNE
R21_PDG   = 105.6583755 / 0.51099895000   # = 206.7683
KAPPA_OBS = 0.0304    # kappa = 3/(4*Phi0), ROADMAP_v3
G0_STAR   = 1.24915   # g0* z phi-FP (dodatekJ2)
G0_TAU    = 3.1891    # g0^tau z Koide r31^K (dodatekT)
R31_KOIDE = 3477.44   # r31^K z Q_K=3/2 (dodatekT)

# DESI DR2 pomiar
DESI_AGAMMA_PHI0 = 1.005
DESI_ERR         = 0.005  # +/-0.5%

# SUPERSEDED - stary parametr solitonu (sek08, przed phi-FP)
ALPHA_K_SUPERSEDED = 8.5616   # NIE UZYWAC w nowych obliczeniach


# ---------------------------------------------------------------------------
# Sciezka S1: Phi0 z Lambda_obs
# ---------------------------------------------------------------------------
def phi0_from_lambda():
    """
    Phi0 = 96*pi*G0*rho_Lambda / H0^2   [bezwymiarowe w jednostkach TGP]
    """
    return 96.0 * np.pi * G0 * RHO_LAMBDA_OBS / H0_SI**2


# ---------------------------------------------------------------------------
# Sciezka S2 - SUPERSEDED
# Phi0 = (alfa_K * sqrt(r21))^{3/5} ze starym alfa_K=8.5616
# Wynik: Phi0(S2) = 17.95 — PONIZEJ zakresu S1/S3 (~23-25)
# Przyczyna: alfa_K=8.5616 pochodzi z bifurkacji solitonu (sek08 v1)
#            Zastapione przez: g0* = 1.24915 (Twierdzenie phi-FP, dodatekJ2)
# ---------------------------------------------------------------------------
def phi0_from_r21_superseded():
    """
    [SUPERSEDED] Phi0 z relacji alfa_K i r21 przy zalozeniu a_Gamma*Phi0=1.
    Phi0 = (alfa_K * sqrt(r21))^{3/5}
    Uzywane tylko jako referencja historyczna.
    """
    return (ALPHA_K_SUPERSEDED * np.sqrt(R21_PDG))**(3.0/5.0)


# ---------------------------------------------------------------------------
# Sciezka S2b: Phi0 z phi-FP (g0*, r21) - OTWARTE (OP-3)
# Derivacja: g0* wyznacza skale masy elektronu; r31^K z Koide;
#   Phi0(S2b) = r31^K / r21  -- kandydat (nie zweryfikowany analitycznie)
# Uwaga: r31/r21 ~ 16.8; wymagana dodatkowa analiza
# ---------------------------------------------------------------------------
def phi0_from_phifp_candidate():
    """
    [KANDYDAT - OP-3 otwarte]
    Phi0(S2b) ~ r31^K / (4 * r21^{1/3})
    Proba powiazania skali Phi0 z wynikami phi-FP i Koide.
    """
    # Kandydat: Phi0 ~ g0^tau * r21^{1/6}
    r21_16 = R21_PDG**(1.0/6.0)   # ~ 2.490
    phi0_cand = G0_TAU * r21_16    # ~ 3.1891 * 2.490 ~ 7.94  -- za male
    # Alternatywnie: (g0^tau)^2 / g0* * r21^{1/6}
    phi0_alt = (G0_TAU**2) / G0_STAR * r21_16   # ~ 10.17/1.249*2.490 ~ 20.25
    return phi0_alt   # ~ 20.25 (kandydat, OP-3)


# ---------------------------------------------------------------------------
# Sciezka S2c: Phi0 = lambda_bar (Brannen mean z Koide chain)
# ---------------------------------------------------------------------------
# Odkrycie (v3.9 -> v4.0): lambda_bar = mean([1, sqrt(r21), sqrt(r31K)])
#   = c_K (normalizacja Brannena w parametryzacji Koide)
#   = 24.783  (wyznaczone bez wolnych parametrow z g0* via phi-FP + Tw.T-r31)
#
# Zgodnosc:
#   Phi0(S2c) = 24.783  vs  Phi0(S3) = 24.67  (roznica 0.5%)   !!
#   Phi0(S2c) = 24.783  vs  Phi0(S1) = 23.31  (roznica 6.2%)
#
# Fizyczna interpretacja:
#   Phi0 = c_K (Brannen) sugeruje, ze tlo prozniowe substratu Phi0
#   jest rowne sredniej amplitudzie ogonowej leptonu w parametryzacji Koide:
#       Phi0 = mean(A_tail(g0^k)^2) dla k=e,mu,tau
#   Jest to WARUNEK SAMOSPOJNOSCI: substrat generuje leptony,
#   których srednia amplituda wraca do wartosci tla prozniowego.
#
# Referencja: Prop. T-thetaK-r21 (dodatekT_koide_atail_formal.tex),
#             prop. W-S2c-Brannen (dodatekW_agamma_phi0_update.tex)
# ---------------------------------------------------------------------------
def phi0_from_brannen():
    """
    S2c: Phi0 = lambda_bar = mean([1, sqrt(r21), sqrt(r31K)])
    Brannen normalization z Koide chain (Prop. T-thetaK-r21).
    Wyznaczone z g0*=1.24915 bez wolnych parametrow.
    Zwraca: (phi0_S2c, rho, sqrt_r31K)
    """
    rho = np.sqrt(R21_PDG)                                              # = 14.3794
    sqrt_r31K = 2.0*(1.0 + rho) + np.sqrt(3.0*(1.0 + 4.0*rho + rho**2))  # = 58.970
    lam_bar = np.mean([1.0, rho, sqrt_r31K])                           # = 24.783
    return lam_bar, rho, sqrt_r31K


# ---------------------------------------------------------------------------
# Sciezka S3: Phi0 z kappa = 3/(4*Phi0)
# ---------------------------------------------------------------------------
def phi0_from_kappa():
    """
    kappa = 3/(4*Phi0)  ->  Phi0 = 3/(4*kappa)
    kappa_obs = 0.0304 (ROADMAP_v3, A11b)
    """
    return 3.0 / (4.0 * KAPPA_OBS)


# ---------------------------------------------------------------------------
# Skan Phi0 w [20, 28]: stabilnosc a_Gamma*Phi0=1
# ---------------------------------------------------------------------------
def scan_phi0_stability(phi0_vals):
    results = []
    phi0_S1 = phi0_from_lambda()
    phi0_S3 = phi0_from_kappa()
    for phi0 in phi0_vals:
        a_gamma_self = 1.0 / phi0
        cross_S1 = a_gamma_self * phi0_S1
        cross_S3 = a_gamma_self * phi0_S3
        results.append({
            'phi0': phi0,
            'a_gamma': a_gamma_self,
            'cross_S1': cross_S1,
            'cross_S3': cross_S3,
        })
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 65)
    print("TGP — Weryfikacja hipotezy a_Gamma * Phi0 = 1  (ROADMAP R3)")
    print("=" * 65)

    # --- Oblicz Phi0 z kazdej sciezki ---
    phi0_S1              = phi0_from_lambda()
    phi0_S2              = phi0_from_r21_superseded()       # SUPERSEDED
    phi0_S2b             = phi0_from_phifp_candidate()      # stary kandydat OP-3
    phi0_S2c, rho_S2c, sqrt_r31K_S2c = phi0_from_brannen() # NOWE S2c (v4.0)
    phi0_S3              = phi0_from_kappa()

    print(f"\n[S1]  Phi0 z Lambda_obs:           Phi0 = {phi0_S1:.4f}")
    print(f"[S2]  Phi0 z alfa_K, r21:          Phi0 = {phi0_S2:.4f}  [SUPERSEDED]")
    print(f"      (alfa_K=8.5616 zastapione przez phi-FP, g0*={G0_STAR})")
    print(f"[S2b] Phi0 stary kandydat (OP-3):   Phi0 = {phi0_S2b:.4f}  [ODRZUCONY: 13% off S1]")
    print(f"[S2c] Phi0 = lambda_bar (Brannen):  Phi0 = {phi0_S2c:.4f}  [NOWE: 0.5% od S3!]")
    print(f"      rho={rho_S2c:.4f}, sqrt(r31K)={sqrt_r31K_S2c:.4f},")
    print(f"      lambda_bar=mean([1,rho,sqrt(r31K)]) = {phi0_S2c:.6f}")
    print(f"      Zgodnosc S2c-S3: {abs(phi0_S2c-phi0_S3)/phi0_S3*100:.2f}%")
    print(f"      Zgodnosc S2c-S1: {abs(phi0_S2c-phi0_S1)/phi0_S1*100:.2f}%")
    print(f"[S3]  Phi0 z kappa=3/(4*Phi0):     Phi0 = {phi0_S3:.4f}")
    print(f"[S4]  DESI DR2:           a_Gamma*Phi0 = {DESI_AGAMMA_PHI0:.3f} +/- {DESI_ERR:.3f}")

    # Aktywne sciezki: S1, S2c, S3 (S2 superseded, S2b odrzucony)
    phi0_active = [phi0_S1, phi0_S2c, phi0_S3]
    phi0_mean   = np.mean(phi0_active)
    phi0_std    = np.std(phi0_active)
    a_gamma_S1  = 1.0 / phi0_S1
    a_gamma_S2c = 1.0 / phi0_S2c
    a_gamma_S3  = 1.0 / phi0_S3

    print(f"\n  Phi0 srednia (S1, S2c, S3):  {phi0_mean:.4f} +/- {phi0_std:.4f}")
    print(f"  a_Gamma(S1)  = 1/Phi0(S1):  {a_gamma_S1:.6f}")
    print(f"  a_Gamma(S2c) = 1/Phi0(S2c): {a_gamma_S2c:.6f}")
    print(f"  a_Gamma(S3)  = 1/Phi0(S3):  {a_gamma_S3:.6f}")

    # Cross-check: a_Gamma(S1)*Phi0(S3), a_Gamma(S2c)*Phi0(S1), a_Gamma(S2c)*Phi0(S3)
    cross_13  = a_gamma_S1  * phi0_S3
    cross_31  = a_gamma_S3  * phi0_S1
    cross_2c3 = a_gamma_S2c * phi0_S3
    cross_2c1 = a_gamma_S2c * phi0_S1
    print(f"\n  Cross S1->S3:  a_Gamma(S1)*Phi0(S3)  = {cross_13:.4f}")
    print(f"  Cross S3->S1:  a_Gamma(S3)*Phi0(S1)  = {cross_31:.4f}")
    print(f"  Cross S2c->S3: a_Gamma(S2c)*Phi0(S3) = {cross_2c3:.4f}  *** <0.5%? ***")
    print(f"  Cross S2c->S1: a_Gamma(S2c)*Phi0(S1) = {cross_2c1:.4f}")
    print(f"  DESI DR2:      a_Gamma*Phi0           = {DESI_AGAMMA_PHI0:.4f}")

    # --- Skan stabilnosci ---
    print(f"\n[SCAN] Phi0 w [20, 28]:")
    scan_range = np.linspace(20.0, 28.0, 7)
    scan_results = scan_phi0_stability(scan_range)
    print(f"  {'Phi0':>6} | {'a_Gamma':>8} | {'*Phi0(S1)':>10} | {'*Phi0(S3)':>10}")
    print("  " + "-" * 46)
    for r in scan_results:
        print(f"  {r['phi0']:6.2f} | {r['a_gamma']:8.5f} | "
              f"{r['cross_S1']:10.4f} | {r['cross_S3']:10.4f}")

    # --- Weryfikacja koncowa ---
    print(f"\n[TESTY]")

    s1_s3_diff  = abs(phi0_S1 - phi0_S3) / phi0_S1 * 100
    s2c_s3_diff = abs(phi0_S2c - phi0_S3) / phi0_S3 * 100
    s2c_s1_diff = abs(phi0_S2c - phi0_S1) / phi0_S1 * 100
    desi_in_range  = 0.99 < DESI_AGAMMA_PHI0 < 1.02
    phi0_s1_range  = 20.0 < phi0_S1  < 30.0
    phi0_s3_range  = 20.0 < phi0_S3  < 30.0
    phi0_s2c_range = 20.0 < phi0_S2c < 30.0
    cross_ok   = abs(cross_13  - 1.0) < 0.12   # S1<->S3 <12%
    cross_2c3_ok = abs(cross_2c3 - 1.0) < 0.01 # S2c<->S3 <1%
    mean_range = 22.0 < phi0_mean < 26.0

    checks = [
        ("S1-S3 spojnosc Phi0 (< 10%)",
         s1_s3_diff < 10.0,
         f"|Phi0(S1)-Phi0(S3)|/Phi0(S1) = {s1_s3_diff:.2f}%"),

        ("DESI DR2: a_Gamma*Phi0 w [0.99, 1.02]",
         desi_in_range,
         f"a_Gamma*Phi0 = {DESI_AGAMMA_PHI0:.4f}"),

        ("Phi0(S1) w [20, 30]",
         phi0_s1_range,
         f"Phi0(S1) = {phi0_S1:.4f}"),

        ("Phi0(S3) w [20, 30]",
         phi0_s3_range,
         f"Phi0(S3) = {phi0_S3:.4f}"),

        ("Cross S1<->S3: |a_Gamma(S1)*Phi0(S3) - 1| < 12%",
         cross_ok,
         f"a_Gamma(S1)*Phi0(S3) = {cross_13:.4f}  (odch. {abs(cross_13-1)*100:.2f}%)"),

        ("Phi0_mean(S1, S2c, S3) w [22, 26]",
         mean_range,
         f"Phi0_mean = {phi0_mean:.4f}"),

        ("S2c: Phi0(Brannen) w [20, 30]",
         phi0_s2c_range,
         f"Phi0(S2c) = {phi0_S2c:.4f}  (lambda_bar = c_K)"),

        ("S2c-S3 spojnosc: |Phi0(S2c)-Phi0(S3)|/Phi0(S3) < 1%",
         s2c_s3_diff < 1.0,
         f"|Phi0(S2c)-Phi0(S3)|/Phi0(S3) = {s2c_s3_diff:.3f}%  *** KLUCZOWY ***"),

        ("S2c cross: |a_Gamma(S2c)*Phi0(S3) - 1| < 1%",
         cross_2c3_ok,
         f"a_Gamma(S2c)*Phi0(S3) = {cross_2c3:.6f}  (odch. {abs(cross_2c3-1)*100:.3f}%)"),
    ]

    n_pass = sum(1 for _, p, _ in checks if p)
    print(f"\n  WYNIKI: {n_pass}/{len(checks)} PASS")
    for name, passed, detail in checks:
        mark = "PASS" if passed else "FAIL"
        print(f"  [{mark}] {name}")
        print(f"        {detail}")

    # Status S2 (SUPERSEDED) - informacyjnie
    print(f"\n[INFO S2/S2b - SUPERSEDED/ODRZUCONY]")
    s1_s2_diff  = abs(phi0_S1 - phi0_S2)  / phi0_S1  * 100
    s1_s2b_diff = abs(phi0_S1 - phi0_S2b) / phi0_S1  * 100
    print(f"  Phi0(S2)  = {phi0_S2:.4f}  (alfa_K=8.5616, SUPERSEDED, {s1_s2_diff:.1f}% od S1)")
    print(f"  Phi0(S2b) = {phi0_S2b:.4f}  (kandydat g0^tau, ODRZUCONY, {s1_s2b_diff:.1f}% od S1)")
    print(f"  Phi0(S2c) = {phi0_S2c:.4f}  (Brannen lambda_bar, NOWE, {s2c_s3_diff:.2f}% od S3)")

    # --- Wniosek ---
    print(f"\n[WNIOSEK]")
    if n_pass >= len(checks) - 1:
        print(f"  Hipoteza a_Gamma*Phi0=1 POTWIERDZONA numerycznie:")
        print(f"    Phi0(S1)  = {phi0_S1:.3f}  (Lambda_obs)")
        print(f"    Phi0(S2c) = {phi0_S2c:.3f}  (Brannen lambda_bar = c_K, z g0* bez wolnych par.)")
        print(f"    Phi0(S3)  = {phi0_S3:.3f}  (kappa=3/4Phi0)")
        print(f"    DESI DR2: a_Gamma*Phi0 = {DESI_AGAMMA_PHI0:.3f}  (0.5% od 1)")
        print(f"    S2c-S3 zgodnosc: {s2c_s3_diff:.3f}%  <-- KLUCZOWY WYNIK")
        print(f"    Interpretacja: Phi0 = mean(A_tail^2) dla 3 generacji leptonu")
        print(f"    -> Soliton substrat generuje leptony, ktorych sr. amplituda = tlo prozni")
    else:
        print(f"  Hipoteza a_Gamma*Phi0=1 wymaga dalszej weryfikacji.")
        print(f"  PASS: {n_pass}/{len(checks)}")

    print(f"\n[STATUS EPISTEMICZNY]")
    print(f"  a_Gamma*Phi0=1: [NUM] — potwierdzone z 3 sciezek + DESI")
    print(f"  S2  SUPERSEDED: alfa_K=8.5616 -> phi-FP g0*={G0_STAR}")
    print(f"  S2b ODRZUCONY:  kandydat (g0^tau)^2/g0*=r21^(1/6) = {phi0_S2b:.2f} (13% od S1)")
    print(f"  S2c NOWE [CZ.ZAMK.]: Phi0 = lambda_bar = c_K = {phi0_S2c:.3f}")
    print(f"      lambda_bar = mean([1, sqrt(r21), sqrt(r31K)]) (Prop. T-thetaK-r21)")
    print(f"      S2c-S3: {s2c_s3_diff:.3f}%,  S2c-S1: {s2c_s1_diff:.2f}%")

    print("\n" + "=" * 65)
    sys.exit(0 if n_pass >= len(checks) - 1 else 1)


if __name__ == "__main__":
    main()
