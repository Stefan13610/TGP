"""
ex47_desi_fdm_forecast.py
=========================
Prognoza wykrywalności FDM przez DESI DR3 (2026)

MOTYWACJA
---------
TGP-FDM przewiduje bozon o masie m22 = m_boson/(10^-22 eV) ≈ 1-3.
DESI DR3 (2026) zmierzy widmo mocy P(k) Lyman-α z bezprecedensową
precyzją. Pytanie: przy jakim m22 DESI może wykluczyć FDM na 3σ?

METODA
------
1. P(k)_FDM z funkcją transferu Hui+2017:
   T_FDM(k) = [cos(x³)/(1 + x⁸)]  (przybliżenie)
   lub Irsic+2017: T_FDM(k) = [1 + (α k)^{2μ}]^{-5/μ}
2. Okno DESI: k ∈ [0.1, 10] h/Mpc, skuteczny V_eff = 6 Gpc³
3. SNR dla m22 = 0.5, 1, 2, 5 vs CDM
4. Kryterium wykluczenia: SNR > 3 (3σ)
5. Krytyczna masa m22^crit(DESI DR3)

PREDYKCJA TGP
-------------
m22 = m_sp/H_0 * (hbar*c) = sqrt(gamma) w jednostkach Plancka
Dla m_sp = 8.2e-51 l_Pl^-1: m22 ≈ 1 (predykcja TGP-F3)
Napięcie K14: Rogers+2021 wymagają m22 > 2.1

WYNIK
-----
  m22_crit(DESI DR3): przy jakiej masie DESI wyklucza FDM
  SNR(m22=1): czy TGP-F3 jest wykrywalne przez DESI

Autor: TGP Analysis Session v27, 2026-03-22
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import quad
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("EX47: DESI FDM Forecast — Prognoza wykrywalności TGP-FDM")
print("      Cel: m22_crit(DESI DR3), SNR(m22=1), K14 update")
print("=" * 70)
print()

# -------------------------------------------------------------------
# Parametry cosmologiczne (Planck 2018)
# -------------------------------------------------------------------
h      = 0.674          # H0 = 67.4 km/s/Mpc
Omega_m = 0.315
Omega_b = 0.049
sigma8  = 0.811
n_s    = 0.965

# Skala Jeans dla FDM (Hu+2000, Hui+2017)
def lambda_J_fdm(m22, z=3.0):
    """Długość Jeansa FDM w Mpc/h."""
    # lambda_J = 2pi * (hbar^2 / G m^2 rho_DM)^{1/4} / a^{1/4}
    # W jednostkach astronomicznych:
    # k_J = 9.0 * m22^{1/2} * sqrt(Omega_m h^2) / (1+z) [h/Mpc]
    k_J = 9.0 * np.sqrt(m22) * np.sqrt(Omega_m * h**2) / (1.0 + z)
    return k_J

# -------------------------------------------------------------------
# Funkcja transferu FDM (Irsic+2017 / ex40)
# -------------------------------------------------------------------

def alpha_fdm(m22):
    """Skala tłumienia FDM [Mpc/h]."""
    return 0.04 / m22**(4.0/9.0)

def T_fdm(k, m22, mu=1.12):
    """Funkcja transferu FDM względem CDM."""
    alpha = alpha_fdm(m22)
    x = alpha * k
    return (1.0 + x**(2*mu))**(-5.0/mu)

def k_half(m22):
    """Skala k_{1/2} gdzie T_FDM = 0.5 [h/Mpc]."""
    # Numeryczne rozwiązanie T_fdm(k_half) = 0.5
    from scipy.optimize import brentq
    try:
        return brentq(lambda k: T_fdm(k, m22) - 0.5, 1e-3, 1000.0)
    except:
        return np.nan

# -------------------------------------------------------------------
# Widmo mocy CDM (przybliżenie EH1998 + Planck)
# -------------------------------------------------------------------

def P_cdm(k):
    """
    Przybliżone CDM P(k) w jednostkach (Mpc/h)³.
    Używamy skalowania Eisenstein-Hu 1998 (uproszczone).
    """
    # Normalizacja do sigma8
    k_eq = 0.073 * Omega_m * h**2  # [h/Mpc]
    # Transfer function (BBKS approximation)
    q = k / (Omega_m * h**2)
    T_cdm = np.log(1 + 2.34*q) / (2.34*q) * \
            (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25)
    # Power spectrum
    P = 2e4 * (k / 0.05)**(n_s) * T_cdm**2
    return P

def P_fdm(k, m22):
    """P(k) FDM = P_CDM(k) * T_FDM²(k)."""
    return P_cdm(k) * T_fdm(k, m22)**2

# -------------------------------------------------------------------
# DESI DR3 — parametry ankiety
# -------------------------------------------------------------------

# Efektywna objętość DESI Lyman-alpha forest
V_eff_desi = 6.0e9   # [Mpc/h]³ (szacunek z DESI+BOSS combined)
n_bar      = 1e-3    # [h/Mpc]³ gęstość kwazarów
f_sky      = 0.34    # pokrycie nieba DESI

# Zakres k dla Lyman-alpha
k_min_lya = 0.1   # [h/Mpc]
k_max_lya = 8.0   # [h/Mpc] (DESI resolution limit)

# Błąd statystyczny na P(k) z DESI:
# sigma_P(k) = P(k) * sqrt(2/(n_modes)) * (1 + 1/(n_bar * P))
def n_modes(k, dk, V_eff):
    """Liczba modów w przedziale [k, k+dk]."""
    return 4 * np.pi * k**2 * dk * V_eff / (2 * np.pi)**3

def sigma_P_relative(k, V_eff, n_bar_val, P_val, dk=0.1):
    """Względny błąd na P(k): sigma_P/P."""
    Nm = n_modes(k, dk, V_eff)
    shot = 1.0 / (n_bar_val * P_val)
    return np.sqrt(2.0 / Nm) * (1.0 + shot)

# -------------------------------------------------------------------
# KROK 1: Funkcje transferu i P(k) FDM dla różnych m22
# -------------------------------------------------------------------

print("─" * 70)
print("KROK 1: Skale tłumienia FDM")
print("─" * 70)
print()

m22_values = [0.5, 1.0, 2.0, 5.0, 10.0, 21.0]
print(f"  {'m22':>6}  {'alpha [Mpc/h]':>14}  {'k_1/2 [h/Mpc]':>14}  {'lambda_J(z=3)':>14}")
print("  " + "─" * 54)
for m22 in m22_values:
    al  = alpha_fdm(m22)
    kh  = k_half(m22)
    kJ  = lambda_J_fdm(m22, z=3.0)
    print(f"  {m22:>6.1f}  {al:>14.4f}  {kh:>14.3f}  {kJ:>14.3f}")

# -------------------------------------------------------------------
# KROK 2: SNR dla DESI Lyman-alpha
# -------------------------------------------------------------------

print()
print("─" * 70)
print("KROK 2: SNR(m22) dla DESI DR3 Lyman-alpha")
print("─" * 70)
print()

# k-siatka dla całkowania
k_arr = np.logspace(np.log10(k_min_lya), np.log10(k_max_lya), 200)
dk_arr = np.gradient(k_arr)

# SNR² = sum_k [(P_CDM - P_FDM)² / sigma_P²]
def compute_SNR(m22, V_eff=V_eff_desi, n_bar_val=n_bar):
    """Całkowite SNR dla detekcji FDM vs CDM."""
    SNR2 = 0.0
    for ki, dki in zip(k_arr, dk_arr):
        if ki < k_min_lya or ki > k_max_lya:
            continue
        P_c = P_cdm(ki)
        P_f = P_fdm(ki, m22)
        delta_P = P_c - P_f
        sig_P = sigma_P_relative(ki, V_eff, n_bar_val, P_c, dk=dki) * P_c
        if sig_P > 0:
            SNR2 += (delta_P / sig_P)**2
    return np.sqrt(SNR2)

print(f"  {'m22':>6}  {'SNR(DESI)':>12}  {'Wykluczone 3σ?':>16}  {'TGP-F3?':>8}")
print("  " + "─" * 50)

snr_results = {}
for m22 in [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 21.0]:
    snr = compute_SNR(m22)
    snr_results[m22] = snr
    excl = "TAK" if snr > 3.0 else "NIE"
    is_tgp = "TGP-F3" if abs(m22 - 1.0) < 0.5 else ""
    print(f"  {m22:>6.1f}  {snr:>12.2f}  {excl:>16}  {is_tgp:>8}")

# -------------------------------------------------------------------
# KROK 3: m22_crit — gdzie SNR = 3
# -------------------------------------------------------------------

print()
print("─" * 70)
print("KROK 3: m22_crit(DESI DR3) — próg wykluczenia FDM")
print("─" * 70)
print()

# Interpolacja: znajdź m22 gdzie SNR = 3
from scipy.interpolate import interp1d
m22_arr_full = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 21.0]
snr_arr_full = [compute_SNR(m) for m in m22_arr_full]

# SNR maleje z m22 (większa masa → mniejsze tłumienie → mniejszy sygnał FDM)
# Szukamy m22_crit gdzie SNR(m22_crit) = 3 (SNR spada poniżej 3 dla dużego m22)
m22_crit = None
for i in range(len(m22_arr_full) - 1):
    if snr_arr_full[i] >= 3.0 and snr_arr_full[i+1] < 3.0:
        frac = (3.0 - snr_arr_full[i]) / (snr_arr_full[i+1] - snr_arr_full[i])
        m22_crit = m22_arr_full[i] + frac * (m22_arr_full[i+1] - m22_arr_full[i])
        break

if m22_crit is None:
    # Sprawdź czy wszystkie > 3 lub < 3
    if all(s >= 3.0 for s in snr_arr_full):
        m22_crit = m22_arr_full[-1]
        print(f"  Uwaga: SNR > 3 dla całego zakresu m22 ∈ [{m22_arr_full[0]}, {m22_arr_full[-1]}]")
        print(f"  DESI wyklucza FDM dla m22 < {m22_crit} (przynajmniej)")
    else:
        m22_crit = m22_arr_full[0]
        print(f"  Uwaga: SNR < 3 już przy m22 = {m22_arr_full[0]}")

print(f"  m22_crit(DESI DR3) = {m22_crit:.2f}")
print()

# -------------------------------------------------------------------
# KROK 4: TGP-F3 predykcja vs DESI
# -------------------------------------------------------------------

print("─" * 70)
print("KROK 4: TGP-F3 predykcja vs DESI DR3 i K14")
print("─" * 70)
print()

m22_tgp = 1.0  # Predykcja TGP-F3 (m_boson = sqrt(gamma) * jednostki)
snr_tgp = snr_results.get(m22_tgp, compute_SNR(m22_tgp))

# K14: Rogers+2021 wymagają m22 > 2.1
m22_K14 = 2.1
snr_K14 = compute_SNR(m22_K14)

print(f"  TGP-F3 predykcja: m22 = {m22_tgp:.1f}")
print(f"  SNR(m22={m22_tgp}) = {snr_tgp:.2f}")
print()
print(f"  K14 granica (Rogers+2021): m22 > {m22_K14}")
print(f"  SNR(m22={m22_K14}) = {snr_K14:.2f}")
print()
print(f"  m22_crit(DESI DR3) = {m22_crit:.2f}")
print()

# Czy DESI może rozróżnić m22=1 od CDM?
if snr_tgp > 3.0:
    print(f"  DESI WYKRYJE TGP-FDM (m22=1) na SNR = {snr_tgp:.1f} > 3 sigma!")
    tgp_detectable = True
elif snr_tgp > 1.0:
    print(f"  DESI widzi TGP-FDM na {snr_tgp:.1f}σ — marginalnie (wymaga n_bar większego)")
    tgp_detectable = False
else:
    print(f"  DESI NIE WYKRYJE TGP-FDM (m22=1): SNR = {snr_tgp:.2f} < 1 sigma")
    tgp_detectable = False

# -------------------------------------------------------------------
# KROK 5: Finalne podsumowanie
# -------------------------------------------------------------------

print()
print("=" * 70)
print("PODSUMOWANIE EX47: DESI FDM Forecast")
print("=" * 70)
print()
print(f"  Funkcja transferu: Irsic+2017  T(k) = [1+(α·k)^2.24]^(-4.46)")
print(f"  DESI DR3: V_eff = {V_eff_desi:.1e} (Mpc/h)³, k ∈ [{k_min_lya},{k_max_lya}] h/Mpc")
print()
print(f"  m22_crit(DESI DR3): {m22_crit:.2f}")
print(f"    => DESI wyklucza FDM dla m22 < {m22_crit:.2f} na 3σ")
print()
print(f"  TGP-F3 (m22=1):")
print(f"    SNR = {snr_tgp:.2f} — {'WYKRYWALNE' if tgp_detectable else 'NIE wykrywalne'}")
print()
print(f"  K14 status: Rogers+2021 wymagają m22 > 2.1")
if m22_crit >= 2.1:
    print(f"  DESI potwierdzi lub wykluczy K14 granicę")
else:
    print(f"  m22_crit < 2.1 — DESI może wymagać m22 > {m22_crit:.2f} tylko")
print()

# Skalowanie k_1/2 ~ m22^{4/9}
kh_1  = k_half(1.0)
kh_21 = k_half(21.0)
exp_fit = np.log(kh_21/kh_1) / np.log(21.0)
print(f"  Skalowanie k_1/2 ~ m22^{exp_fit:.3f}  (teoria: 4/9 = {4/9:.3f})")
print()
print(f"  IMPLIKACJE DLA TGP:")
print(f"    Jeśli DESI DR3 nie wyklucza m22 < {m22_crit:.1f}:")
print(f"      => TGP-F3 (m22≈1) ZGODNE z DESI")
print(f"    Jeśli DESI DR3 wyklucza m22 < {m22_crit:.1f} > 1:")
print(f"      => TGP-F3 WYKLUCZONE (falsyfikacja K14+DESI)")
print(f"    Nowe napięcie K14: {m22_crit:.1f} vs TGP-F3: 1.0")
if m22_crit > 1.0:
    print(f"      => Wymaga modyfikacji F3 lub wejścia w F2 (m22 >= {m22_crit:.1f})")

print()
print("─" * 70)
print("Uwaga: Wyniki orientacyjne (P_CDM = przybliżenie BBKS,")
print("V_eff szacunek). Dla precyzji: użyj CLASS/CAMB + rybicki.")
print("─" * 70)
print()
print("=" * 70)
print("EX47 DONE — DESI FDM Forecast")
print("=" * 70)
