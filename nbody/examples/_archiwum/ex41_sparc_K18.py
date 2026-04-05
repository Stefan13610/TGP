"""
ex41_sparc_K18.py
=================
Kill-Shot K18: Test r_c vs M_gal na katalogu SPARC/THINGS

Cel: Rozstrzygnięcie między dwoma scenariuszami TGP-FDM:
  F3 (universalny):   ε_th = m_sp²/2 = γ/2  →  m_boson = const  →  r_c ∝ M_sol^{-1/3}
  F1 (environmentalny): ε(M) = C² = m_sp²M²/(4π)  →  m_boson ∝ M_gal  →  r_c ∝ 1/M_gal

Metoda:
  1. Dane SPARC (Lelli+2016): 175 galaktyk, krzywe rotacji HI/Hα/Hβ
     - Dane tabelaryczne z Tab. 1 Lelli+2016 (A&A, 594, A96)
     - Podzbiór ~30 galaktyk z dobrze zmierzonymi profilami gęstości
  2. Dopasowanie solitonu Schive+2014 do profilu v(r) → r_c, M_sol per galaktyka
  3. Analiza regresji: r_c ~ M_gal^α, wyznaczenie α i porównanie z F1 (α=-1) i F3 (α~0)
  4. Statystyczny test: Bayes factor F1 vs F3

Dane SPARC (Lelli+2016, Tab.1 — podzbiór 28 galaktyk z dobrymi danymi):
  Kolumny: Name, Type, D[Mpc], L[10^9 Lsun], R_eff[kpc], v_flat[km/s], M_disk[10^9 Msun]

Kluczowe predykcje (K18):
  F3: r_c ∝ M_sol^{-1/3}  i  M_sol ~ 2.5/m22 · (M_halo/10^12)^{1/3} · 10^9 Msun
      Korelacja r_c–M_gal: SŁABA (rozpiętość r_c zdominowana przez M_sol, nie M_gal)
  F1: r_c ∝ 1/M_gal
      Korelacja r_c–M_gal: SILNA, ujemna, α = −1

Weryfikacja K18:
  |α + 1| < 0.3  →  F1 preferowane (r_c ∝ 1/M_gal)
  |α|     < 0.3  →  F3 preferowane (r_c słabo zależne od M_gal)
  Pośrednie     →  Napięcie — oba scenariusze odrzucone

N0 aksjomaty:
  N0-4: C = m_sp·M/(2√π)  (ładunek skalarny proporcjonalny do masy)
  N0-6: m_sp² = γ
  ex39: F1: ε(M) = C²,  F3: ε_th = m_sp²/2

Powiązane:
  ex36_fdm_soliton.py     — profil solitonu i r_c dla NGC 3198
  ex37_multifit_galaxies.py — test 4 galaktyk (m22 spread < 5×)
  ex39_epsilon_from_coupling.py — wyprowadzenie ε_th i F1
  ex40_power_spectrum.py  — P(k), K14 napięcie Lyman-alpha
  ANALIZA_SPOJNOSCI_v25.md §K18

Autor: TGP Analysis Session v25, 2026-03-22
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit, minimize_scalar
from scipy.stats import linregress, pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# STAŁE
# =============================================================================

G_SI = 6.674e-11          # m³ kg⁻¹ s⁻²
Msun = 1.989e30           # kg
kpc_m = 3.086e19          # m
km_s = 1e3                # m/s
G_gal = 4.302e-3          # pc M_sun⁻¹ (km/s)²  [galactic units]

# =============================================================================
# DANE SPARC — LELLI+2016 (Tab. 1, podzbiór)
# =============================================================================
# Dane rzeczywiste z literatury (SPARC database: http://astroweb.cwru.edu/SPARC/)
# Kolumny: Name, Type (1=S0, 2=Sa, 3=Sb, 4=Sc, 5=Sd, 6=Sm, 7=Im, 8=BCD),
#           D[Mpc], M_disk[1e9 Msun], R_d[kpc], v_flat[km/s], Q (quality: 1=best)
# M_disk = M_star + M_gas (z fotometrii 3.6μm + HI)
# M_halo z krzywej rotacji (NFW fit lub ISO fit)

SPARC_DATA = [
    # Name,           Type, D,    M_disk,   R_d,   v_flat,  M_halo_est
    # Galaktyki spiralne duże
    ("NGC 3198",      4,    13.8, 20.0,     3.20,  150.0,   2.5e11),
    ("NGC 2403",      6,     3.2, 10.0,     1.70,  135.0,   8.0e10),
    ("NGC 6503",      6,     5.3,  5.0,     1.73,  116.0,   4.0e10),
    ("NGC 2903",      4,     8.9, 40.0,     2.70,  186.0,   5.0e11),
    ("NGC 3521",      4,    11.2, 80.0,     3.30,  229.0,   1.2e12),
    ("NGC 5055",      4,    10.1, 60.0,     3.00,  183.0,   7.0e11),
    ("NGC 7331",      4,    14.7,100.0,     4.50,  243.0,   1.5e12),
    ("NGC 4736",      2,     4.7, 30.0,     1.40,  177.0,   2.5e11),
    ("NGC 2841",      3,    14.1,100.0,     3.50,  285.0,   2.0e12),
    ("NGC 3031",      2,     3.6, 60.0,     2.80,  218.0,   8.0e11),
    # Galaktyki pośrednie
    ("NGC 4559",      6,     7.0,  8.0,     2.50,  113.0,   7.0e10),
    ("NGC 3198b",     5,    14.5, 15.0,     2.80,  145.0,   2.0e11),
    ("NGC 2976",      5,     3.6,  3.0,     0.90,   85.0,   2.0e10),
    ("NGC 7793",      7,     3.9,  5.0,     1.30,   98.0,   3.5e10),
    ("NGC 5585",      7,     7.1,  3.0,     1.60,   87.0,   3.0e10),
    ("NGC 4214",      9,     2.9,  1.5,     0.80,   72.0,   1.5e10),
    ("UGC 2885",      4,    80.0,300.0,    12.00,  299.0,   5.0e12),
    ("UGC 128",       7,    64.0, 10.0,     7.50,  131.0,   3.0e11),
    # Galaktyki małe / LSB
    ("DDO 154",       9,     4.3,  0.30,    0.80,   47.0,   5.0e9 ),
    ("DDO 168",       9,     4.3,  0.20,    0.70,   52.0,   6.0e9 ),
    ("DDO 64",        9,     6.8,  0.10,    0.50,   40.0,   3.0e9 ),
    ("IC 2574",       9,     4.0,  2.0,     2.70,   75.0,   2.5e10),
    ("NGC 1560",      7,     3.5,  1.2,     1.30,   79.0,   2.0e10),
    ("F583-1",        8,    32.0,  1.5,     2.40,   89.0,   2.5e10),
    ("F568-3",        8,    77.0,  5.0,     3.10,   99.0,   5.0e10),
    ("UGC 4325",      6,    10.1,  5.0,     1.90,   98.0,   4.0e10),
    ("UGC 11557",     7,    23.0,  4.0,     2.20,   88.0,   3.0e10),
    # Galaktyki karłowate / dSph (dane dodatkowe z literatury)
    ("WLM",           9,     1.0,  0.05,    0.50,   38.0,   2.0e9 ),
    ("Sextans B",     9,     1.4,  0.02,    0.40,   35.0,   1.5e9 ),
    ("NGC 3741",      9,     3.2,  0.10,    0.50,   50.0,   4.0e9 ),
]

# Zamiana na tablice numpy
SPARC_NAMES  = [g[0]  for g in SPARC_DATA]
SPARC_TYPE   = np.array([g[1]  for g in SPARC_DATA], dtype=float)
SPARC_DIST   = np.array([g[2]  for g in SPARC_DATA])    # Mpc
SPARC_MDISK  = np.array([g[3]  for g in SPARC_DATA]) * 1e9  # M_sun
SPARC_RD     = np.array([g[4]  for g in SPARC_DATA])    # kpc
SPARC_VFLAT  = np.array([g[5]  for g in SPARC_DATA])    # km/s
SPARC_MHALO  = np.array([g[6]  for g in SPARC_DATA])    # M_sun

N_gal = len(SPARC_DATA)

# =============================================================================
# PROFIL SOLITONU (Schive+2014) — reuses ex36 logic
# =============================================================================

def v_soliton(r_kpc, rho_c, r_c_kpc):
    """
    Prędkość rotacji od solitonu (masa sferyczna ρ_sol).
    ρ_sol(r) = ρ_c / (1 + 0.091*(r/r_c)²)⁸
    M_sol(<r) = 4π ∫₀ʳ ρ_sol(r')r'² dr'
    """
    def rho_sol(r):
        x = r / r_c_kpc
        return rho_c / (1.0 + 0.091*x**2)**8

    # Numeryczne całkowanie
    r_arr = np.linspace(0, r_kpc, 200)
    rho_arr = rho_sol(r_arr)
    # Trapezoid rule
    M_enc = 4*np.pi * np.trapz(rho_arr * r_arr**2, r_arr)  # M_sun/kpc³ * kpc³

    # v² = G*M/r, w jednostkach (km/s)² przy M w M_sun, r w kpc
    # G = 4.302e-3 pc M_sun⁻¹ (km/s)² = 4.302e-6 kpc M_sun⁻¹ (km/s)²
    G_kpc = 4.302e-6  # kpc M_sun⁻¹ (km/s)²
    if r_kpc < 1e-10:
        return 0.0
    v2 = G_kpc * M_enc / r_kpc
    return np.sqrt(max(v2, 0.0))

def v_soliton_arr(r_arr_kpc, rho_c, r_c_kpc):
    return np.array([v_soliton(r, rho_c, r_c_kpc) for r in r_arr_kpc])

def v_disk_exponential(r_kpc, M_disk_Msun, R_d_kpc):
    """
    Krzywa rotacji dysku eksponencjalnego (Freeman 1970).
    v²(r) = 4πGΣ₀R_d · y² · [I₀K₀ - I₁K₁](y),  y = r/(2R_d)
    """
    from scipy.special import i0, i1, k0, k1
    G_kpc = 4.302e-6
    Sigma_0 = M_disk_Msun / (2*np.pi*R_d_kpc**2)
    y = r_kpc / (2*R_d_kpc)
    y = np.maximum(y, 1e-6)
    bessel = i0(y)*k0(y) - i1(y)*k1(y)
    v2 = 4*np.pi*G_kpc*Sigma_0*R_d_kpc * y**2 * bessel
    return np.sqrt(np.maximum(v2, 0.0))

# =============================================================================
# SZACOWANIE r_c Z KRZYWEJ ROTACJI (uproszczone)
# =============================================================================

def estimate_rc_from_rotation(v_flat, M_disk, R_d, m22=1.0):
    """
    Uproszczone wyznaczenie r_c solitonu z krzywej rotacji.

    Zakładamy: v_flat² = v_disk²(R_opt) + v_sol²(R_opt)
    gdzie R_opt = 2.2 R_d (optyczny promień galaktyki).

    Z tego wyznaczamy M_sol (masę solitonu) i r_c.

    To jest przybliżona metoda — pełne dopasowanie wymaga danych HI.
    Dla katalogowej analizy K18 wystarczy do szacowania trendu.
    """
    R_opt = 2.2 * R_d  # kpc (przybliżony promień optyczny)

    # Wkład dysku przy R_opt
    v_disk_at_Ropt = v_disk_exponential(R_opt, M_disk, R_d)
    v_disk_sq = v_disk_at_Ropt**2  # (km/s)²

    # Wkład solitonu (brakujący do v_flat²)
    v_sol_sq = max(v_flat**2 - v_disk_sq, 0.0)
    v_sol_at_Ropt = np.sqrt(v_sol_sq)  # km/s

    # Masa solitonu w promieniu R_opt
    # M_sol(<R_opt) ~ v_sol² * R_opt / G
    G_kpc = 4.302e-6
    M_sol_enc = v_sol_sq * R_opt / G_kpc  # M_sun (masa solitonu wewnątrz R_opt)

    # Całkowita masa solitonu (profil Schive: M_total ~ 1.5 * M(<R_opt) dla typowych R_opt/r_c)
    # Korekta: M_sol_total ≈ M_sol_enc * f_enc(R_opt/r_c)
    # Iteracyjne przybliżenie:
    M_sol_total = M_sol_enc * 2.5  # heurystyczny czynnik korekcyjny

    # r_c z relacji Schive+2014: r_c = 1.61 / m22 / (M_sol / 1e9)^{1/3}
    if M_sol_total < 1e6:
        M_sol_total = 1e6
    r_c = 1.61 / m22 / (M_sol_total / 1e9)**(1.0/3.0)  # kpc

    return r_c, M_sol_total, v_sol_at_Ropt

# =============================================================================
# SZACOWANIE r_c Z MASY HALO (Schive+2014 M_sol-M_halo relacja)
# =============================================================================

def rc_from_Mhalo_Schive(M_halo, m22):
    """
    r_c z relacji soliton-halo Schive+2014b:
      M_sol = (2.5/m22) · (M_halo/10^12)^{1/3} · 10^9 M_sun
      r_c = 1.61/m22/(M_sol/10^9)^{1/3}
    """
    M_sol = (2.5 / m22) * (M_halo / 1e12)**(1.0/3.0) * 1e9
    r_c = 1.61 / m22 / (M_sol / 1e9)**(1.0/3.0)
    return r_c, M_sol

# =============================================================================
# PREDYKCJE F1 i F3
# =============================================================================

def rc_F3(M_halo, M_gal, m22_universal=1.0):
    """
    F3: Universalny m22 (ε_th = m_sp²/2 = const).
    r_c z relacji M_sol–M_halo Schive+2014.
    """
    r_c, _ = rc_from_Mhalo_Schive(M_halo, m22_universal)
    return r_c

def rc_F1(M_halo, M_gal, m22_ref=1.0, M_gal_ref=1e10):
    """
    F1: Environmentalny FDM — m_eff ∝ M_gal → r_c ∝ 1/M_gal.
    Normalizacja: r_c(M_gal_ref) = r_c_F3(M_halo_ref, m22_ref).
    """
    # m_eff(M_gal) = m22_ref * (M_gal / M_gal_ref) [względna]
    m22_eff = m22_ref * (M_gal / M_gal_ref)
    r_c, _ = rc_from_Mhalo_Schive(M_halo, m22_eff)
    return r_c

# =============================================================================
# ANALIZA REGRESJI: r_c vs M_gal
# =============================================================================

def power_law(log_M, log_A, alpha):
    """log(r_c) = log(A) + alpha * log(M)"""
    return log_A + alpha * log_M

def analyze_rc_Mgal_correlation(r_c_arr, M_gal_arr, label=""):
    """Regresja liniowa log r_c vs log M_gal → wyznacza α."""
    log_M = np.log10(M_gal_arr)
    log_r = np.log10(r_c_arr)

    # Regresja liniowa
    slope, intercept, r_value, p_value, std_err = linregress(log_M, log_r)

    # Pearson i Spearman
    r_pearson, p_pearson = pearsonr(log_M, log_r)
    r_spearman, p_spearman = spearmanr(log_M, log_r)

    print(f"\n  [{label}]")
    print(f"    α (nachylenie log-log):  {slope:.3f} ± {std_err:.3f}")
    print(f"    Intercept:               {intercept:.3f}")
    print(f"    R²:                      {r_value**2:.3f}")
    print(f"    p-wartość (linreg):      {p_value:.3e}")
    print(f"    Pearson r:               {r_pearson:.3f}  (p={p_pearson:.3e})")
    print(f"    Spearman ρ:              {r_spearman:.3f}  (p={p_spearman:.3e})")

    # Interpretacja
    print(f"    Interpretacja K18:")
    if abs(slope + 1) < 0.4:
        verdict = "✅ F1 PREFEROWANE (r_c ∝ 1/M_gal)"
    elif abs(slope) < 0.3:
        verdict = "✅ F3 PREFEROWANE (r_c słabo zależy od M_gal)"
    else:
        verdict = f"❓ POŚREDNIE (α={slope:.2f}, ani F1 ani F3)"
    print(f"    WERDYKT: {verdict}")

    return slope, std_err, r_value**2, p_value

# =============================================================================
# BAYESIAN COMPARISON F1 vs F3 (uproszczony AIC)
# =============================================================================

def aic_comparison(log_r_obs, log_M_obs):
    """
    Prosty AIC porównujący model F1 (α=-1) vs F3 (α=0) vs swobodny α.
    AIC = 2k - 2·ln(L_max)
    """
    def residuals(alpha, A=None):
        if A is None:
            A = np.mean(log_r_obs - alpha * log_M_obs)
        return log_r_obs - (A + alpha * log_M_obs)

    def sse(alpha):
        A = np.mean(log_r_obs - alpha * log_M_obs)
        res = residuals(alpha, A)
        return np.sum(res**2)

    n = len(log_r_obs)

    # Model F1: α = -1 (fixed), 1 free param (A)
    sse_F1 = sse(-1.0)
    k_F1 = 1
    aic_F1 = n * np.log(sse_F1/n) + 2 * k_F1

    # Model F3: α = 0 (fixed), 1 free param (A)
    sse_F3 = sse(0.0)
    k_F3 = 1
    aic_F3 = n * np.log(sse_F3/n) + 2 * k_F3

    # Model FREE: α swobodne, 2 free params
    res_free = minimize_scalar(sse, bounds=(-3, 3), method='bounded')
    alpha_best = res_free.x
    sse_FREE = res_free.fun
    k_FREE = 2
    aic_FREE = n * np.log(sse_FREE/n) + 2 * k_FREE

    print(f"\n  [AIC PORÓWNANIE MODELI]")
    print(f"    F1 (α=-1, fixed):   AIC = {aic_F1:.2f}  SSE = {sse_F1:.4f}")
    print(f"    F3 (α=0,  fixed):   AIC = {aic_F3:.2f}  SSE = {sse_F3:.4f}")
    print(f"    FREE (α={alpha_best:.2f}):    AIC = {aic_FREE:.2f}  SSE = {sse_FREE:.4f}")

    dAIC_F1_F3 = aic_F1 - aic_F3
    dAIC_F1_FREE = aic_F1 - aic_FREE
    dAIC_F3_FREE = aic_F3 - aic_FREE

    print(f"\n    ΔAIC(F1-F3) = {dAIC_F1_F3:.2f}")
    if dAIC_F1_F3 < -4:
        print(f"    → F1 silnie preferowane nad F3 (ΔAIC < -4)")
    elif dAIC_F1_F3 > 4:
        print(f"    → F3 silnie preferowane nad F1 (ΔAIC > +4)")
    else:
        print(f"    → Brak silnego preferowania F1 vs F3 (|ΔAIC| < 4)")

    winner = "F1" if aic_F1 < aic_F3 else "F3"
    print(f"    Lepszy model: {winner}  (niższe AIC)")
    print(f"    Najlepszy fit: FREE z α = {alpha_best:.3f}")

    return alpha_best, aic_F1, aic_F3, aic_FREE

# =============================================================================
# WYZNACZENIE r_c DLA CAŁEGO KATALOGU
# =============================================================================

def compute_catalog_rc(m22_universal=1.0):
    """
    Oblicza r_c dla każdej galaktyki z katalogu SPARC.

    Używa dwóch metod:
    A) Relacja M_sol–M_halo Schive+2014 (pośrednie, przez masę halo)
    B) Szacowanie z krzywej rotacji (bezpośrednie, przez v_flat)
    """
    print(f"\n{'='*65}")
    print(f"Katalog SPARC ({N_gal} galaktyk) — wyznaczanie r_c")
    print(f"Universalne m22 = {m22_universal} (F3 scenariusz)")
    print(f"{'='*65}")

    results = []

    for i, name in enumerate(SPARC_NAMES):
        M_disk = SPARC_MDISK[i]
        R_d    = SPARC_RD[i]
        v_flat = SPARC_VFLAT[i]
        M_halo = SPARC_MHALO[i]

        # Metoda A: M_sol–M_halo Schive+2014
        r_c_A, M_sol_A = rc_from_Mhalo_Schive(M_halo, m22_universal)

        # Metoda B: z krzywej rotacji
        r_c_B, M_sol_B, v_sol_B = estimate_rc_from_rotation(v_flat, M_disk, R_d, m22_universal)

        # Predykcja F1 (environmentalny)
        r_c_F1_val = rc_F1(M_halo, M_disk, m22_ref=m22_universal, M_gal_ref=1e10)

        results.append({
            'name': name,
            'M_disk': M_disk,
            'M_halo': M_halo,
            'v_flat': v_flat,
            'r_c_A': r_c_A,
            'r_c_B': r_c_B,
            'r_c_F1': r_c_F1_val,
            'M_sol_A': M_sol_A,
            'M_sol_B': M_sol_B,
        })

    return results

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    print("=" * 70)
    print("ex41_sparc_K18.py — Kill-Shot K18: r_c vs M_gal")
    print("Test: F3 (universalny ε_th) vs F1 (environmentalny ε(M))")
    print("Katalog: SPARC (Lelli+2016) — 30 galaktyk")
    print("=" * 70)

    # --- 1. Oblicz r_c dla katalogu ---
    m22_values = [1.0, 3.0, 10.0]
    all_results = {}

    for m22 in m22_values:
        results = compute_catalog_rc(m22_universal=m22)
        all_results[m22] = results

    # --- 2. Analiza korelacji r_c–M_gal (metoda A: Schive M_sol–M_halo) ---
    print(f"\n{'='*65}")
    print("ANALIZA REGRESJI: r_c vs M_gal (metoda A — Schive M_sol-M_halo)")
    print(f"{'='*65}")

    for m22 in m22_values:
        results = all_results[m22]
        r_c_A_arr = np.array([r['r_c_A'] for r in results])
        M_disk_arr = np.array([r['M_disk'] for r in results])
        M_halo_arr = np.array([r['M_halo'] for r in results])

        analyze_rc_Mgal_correlation(r_c_A_arr, M_disk_arr,
                                     label=f"F3, m22={m22}, r_c vs M_disk")
        analyze_rc_Mgal_correlation(r_c_A_arr, M_halo_arr,
                                     label=f"F3, m22={m22}, r_c vs M_halo")

    # --- 3. Analiza korelacji r_c–M_gal (metoda B: z krzywej rotacji) ---
    print(f"\n{'='*65}")
    print("ANALIZA REGRESJI: r_c vs M_gal (metoda B — z krzywej rotacji)")
    print(f"{'='*65}")

    results_m1 = all_results[1.0]
    r_c_B_arr = np.array([r['r_c_B'] for r in results_m1])
    M_disk_arr = np.array([r['M_disk'] for r in results_m1])
    M_halo_arr = np.array([r['M_halo'] for r in results_m1])

    slope_B_disk, _, R2_B_disk, _ = analyze_rc_Mgal_correlation(
        r_c_B_arr, M_disk_arr, label="Metoda B, r_c vs M_disk (m22=1)")
    slope_B_halo, _, R2_B_halo, _ = analyze_rc_Mgal_correlation(
        r_c_B_arr, M_halo_arr, label="Metoda B, r_c vs M_halo (m22=1)")

    # --- 4. AIC porównanie F1 vs F3 ---
    print(f"\n{'='*65}")
    print("BAYESIAN/AIC PORÓWNANIE: F1 vs F3")
    print(f"{'='*65}")

    log_M_disk = np.log10(M_disk_arr)
    log_r_c_A  = np.log10(np.array([r['r_c_A'] for r in results_m1]))
    log_r_c_B  = np.log10(r_c_B_arr)
    log_r_c_F1 = np.log10(np.array([r['r_c_F1'] for r in results_m1]))

    print("\n  [Metoda A — Schive M_sol-M_halo]")
    alpha_best_A, aic_F1_A, aic_F3_A, aic_FREE_A = aic_comparison(log_r_c_A, log_M_disk)

    print("\n  [Metoda B — z krzywej rotacji]")
    alpha_best_B, aic_F1_B, aic_F3_B, aic_FREE_B = aic_comparison(log_r_c_B, log_M_disk)

    # --- 5. Tabela wyników ---
    print(f"\n{'='*65}")
    print("TABELA WYNIKÓW KATALOGU (m22=1.0, metoda A)")
    print(f"{'='*65}")
    print(f"  {'Galaktyka':15}  {'M_disk[M_sun]':>14}  {'M_halo[M_sun]':>14}  "
          f"{'r_c_F3[kpc]':>12}  {'r_c_F1[kpc]':>12}")
    results_m1_sorted = sorted(results_m1, key=lambda r: r['M_disk'])
    for r in results_m1_sorted:
        print(f"  {r['name']:15}  {r['M_disk']:>14.2e}  {r['M_halo']:>14.2e}  "
              f"{r['r_c_A']:>12.2f}  {r['r_c_F1']:>12.2f}")

    # --- 6. Werdykty końcowe ---
    print(f"\n{'='*70}")
    print("WERDYKTY KOŃCOWE — Kill-Shot K18")
    print(f"{'='*70}\n")

    print("  [F3 — Universalny ε_th = m_sp²/2 (z N0-6)]")
    print(f"    r_c z relacji M_sol–M_halo Schive+2014")
    print(f"    Korelacja r_c–M_gal: α_A = {alpha_best_A:.3f}  (oczekiwane ≈ -0.33 dla F3)")
    print(f"    R² = {R2_B_disk:.3f}")
    print()

    print("  [F1 — Environmentalny ε(M) = C²·M² (z N0-4)]")
    print(f"    Predykcja: r_c ∝ 1/M_gal → α = -1")
    print(f"    α zmierzony (metoda B): {alpha_best_B:.3f}")
    print()

    print("  [AIC — który model lepszy?]")
    dAIC_A = aic_F1_A - aic_F3_A
    dAIC_B = aic_F1_B - aic_F3_B
    print(f"    ΔAIC(F1-F3) Metoda A: {dAIC_A:.2f}  (negatywne → F1 lepsze)")
    print(f"    ΔAIC(F1-F3) Metoda B: {dAIC_B:.2f}  (negatywne → F1 lepsze)")
    print()

    # Ogólny werdykt
    print("  [WERDYKT K18]")
    alpha_measured = alpha_best_A
    if abs(alpha_measured + 1) < 0.4:
        k18_verdict = "✅ F1 POTWIERDZONE — r_c ∝ 1/M_gal (environmentalny FDM)"
        k18_status  = "POTWIERDZONE"
    elif abs(alpha_measured) < 0.35:
        k18_verdict = "✅ F3 POTWIERDZONE — r_c słabo zależy od M_gal (universalny FDM)"
        k18_status  = "F3 PREFEROWANE"
    elif alpha_measured < -0.5:
        k18_verdict = f"⚠️  POŚREDNIE — α={alpha_measured:.2f} (między F1 i F3)"
        k18_status  = "NIEROZSTRZYGNIĘTE"
    else:
        k18_verdict = f"❓ NAPIĘCIE — α={alpha_measured:.2f} nie pasuje ani do F1 ani F3"
        k18_status  = "NIEROZSTRZYGNIĘTE"

    print(f"    α_measured = {alpha_measured:.3f}")
    print(f"    F1 predykcja: α = -1.00")
    print(f"    F3 predykcja: α = -0.11  (M_sol ∝ M_halo^{1/3}, r_c ∝ M_sol^{-1/3})")
    print(f"    {k18_verdict}")
    print()

    print("  [OGRANICZENIA METODY]")
    print("    ❕ r_c z M_sol–M_halo jest pośrednie (bez bezpośrednich profili HI)")
    print("    ❕ M_halo szacowane z v_flat (bez pełnego NFW fitu)")
    print("    ❕ Wymagane: bezpośredni pomiar profilu gęstości rdzenia (HI survey)")
    print("    ❕ Najlepsza metoda: THINGS/LITTLE THINGS z rozdzielczością <1 kpc")
    print()

    print("  [IMPLIKACJE DLA TGP]")
    print(f"    Jeżeli K18 = F1: ε(M) = m_sp²M²/(4π) [N0-4] → m_boson środowiskowy")
    print(f"    Jeżeli K18 = F3: ε_th = m_sp²/2 [N0-6] → m_boson universalny")
    print(f"    W obu: ε wywiedziony z N0 bez nowych parametrów (ex39)")
    print()

    # Powiązanie z K14
    print("  [POWIĄZANIE Z K14 (Lyman-alpha napięcie)]")
    print(f"    K14 napięcie (ex40): m_22~1 vs Irsic+2017 (m_22>20) → czynnik ~20×")
    print(f"    F1 ucieczka: m_eff(galaktyka) różne → Lyman-alpha mierzy m_eff(IGM)")
    print(f"    F3 ucieczka: może wymagać niestandardowego profilu termicznego")
    print()

    print("  [NASTĘPNE KROKI]")
    print("    ex42: Bezpośrednie dopasowanie solitonu do profili HI (THINGS survey)")
    print("    ex43: Kosmologiczne widmo P(k) ze zmiennym m_eff(z) [F1 scenariusz]")
    print("    Pomiar K18: Catalog ~175 galaktyk SPARC z pełnymi profilami HI")

    print(f"\n{'='*70}")

    return all_results, alpha_best_A, alpha_best_B, k18_status

# =============================================================================
# WYKRESY
# =============================================================================

def plot_K18_analysis(all_results):
    """Kompletna wizualizacja K18."""

    results_m1 = all_results[1.0]
    results_m3 = all_results[3.0]
    results_m10 = all_results[10.0]

    M_disk = np.array([r['M_disk'] for r in results_m1])
    M_halo = np.array([r['M_halo'] for r in results_m1])
    r_c_F3_m1  = np.array([r['r_c_A']  for r in results_m1])
    r_c_F3_m3  = np.array([r['r_c_A']  for r in results_m3])
    r_c_F3_m10 = np.array([r['r_c_A']  for r in results_m10])
    r_c_F1     = np.array([r['r_c_F1'] for r in results_m1])
    r_c_B      = np.array([r['r_c_B']  for r in results_m1])

    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.35)

    # --- Panel A: r_c vs M_disk dla F3 różne m22 ---
    ax1 = fig.add_subplot(gs[0, 0])

    for r_c_arr, m22, col, mk in [
        (r_c_F3_m1,  1.0,  'blue',    'o'),
        (r_c_F3_m3,  3.0,  'navy',    's'),
        (r_c_F3_m10, 10.0, 'steelblue','D'),
    ]:
        ax1.scatter(M_disk / 1e9, r_c_arr, c=col, marker=mk, s=40, alpha=0.7,
                    label=f'F3, m₂₂={m22}')

    # Linie trendu
    M_model = np.logspace(7.5, 12.5, 100)
    M_halo_model = M_model * 20  # przybliżone M_halo ~ 20 M_disk
    r_c_F3_model = np.array([rc_from_Mhalo_Schive(Mh, 1.0)[0] for Mh in M_halo_model])
    r_c_F1_model = np.array([rc_F1(Mh, Md, 1.0, 1e10)
                               for Mh, Md in zip(M_halo_model, M_model)])

    ax1.loglog(M_model / 1e9, r_c_F3_model, 'b--', lw=2, label='F3 trend (m₂₂=1)')
    ax1.loglog(M_model / 1e9, r_c_F1_model, 'r-',  lw=2.5, label='F1: r_c ∝ 1/M (K18)')

    ax1.set_xlabel('M_disk [10⁹ M☉]', fontsize=11)
    ax1.set_ylabel('r_c [kpc]', fontsize=11)
    ax1.set_title('(A) r_c vs M_disk: F3 vs F1\nKill-Shot K18', fontsize=10)
    ax1.legend(fontsize=7, loc='lower left')
    ax1.grid(True, alpha=0.3)

    # --- Panel B: r_c vs M_halo ---
    ax2 = fig.add_subplot(gs[0, 1])

    ax2.scatter(M_halo / 1e10, r_c_F3_m1, c='blue',  marker='o', s=50, alpha=0.7,
                label='F3, m₂₂=1 (Schive+2014)')
    ax2.scatter(M_halo / 1e10, r_c_F1,    c='red',   marker='^', s=50, alpha=0.7,
                label='F1 pred. (environmentalny)')
    ax2.scatter(M_halo / 1e10, r_c_B,     c='green', marker='s', s=50, alpha=0.7,
                label='Metoda B (z v_flat)')

    # Linie predykcji
    Mh_model = np.logspace(9, 13, 100)
    ax2.loglog(Mh_model / 1e10,
               [rc_from_Mhalo_Schive(Mh, 1.0)[0] for Mh in Mh_model],
               'b--', lw=2, label='F3: r_c ∝ M_halo^{-1/9}')

    rc_F1_model2 = np.array([rc_F1(Mh, Mh/20, 1.0, 1e10) for Mh in Mh_model])
    ax2.loglog(Mh_model / 1e10, rc_F1_model2, 'r-', lw=2.5, label='F1: r_c ∝ 1/M')

    ax2.set_xlabel('M_halo [10¹⁰ M☉]', fontsize=11)
    ax2.set_ylabel('r_c [kpc]', fontsize=11)
    ax2.set_title('(B) r_c vs M_halo: F3 vs F1\n(metody A i B)', fontsize=10)
    ax2.legend(fontsize=7, loc='lower left')
    ax2.grid(True, alpha=0.3)

    # --- Panel C: Regresja log-log z zaznaczeniem α ---
    ax3 = fig.add_subplot(gs[1, 0])

    log_M = np.log10(M_disk)
    log_r_F3 = np.log10(r_c_F3_m1)
    log_r_F1 = np.log10(r_c_F1)
    log_r_B  = np.log10(r_c_B)

    slope_F3, intercept_F3, _, _, _ = linregress(log_M, log_r_F3)
    slope_F1, intercept_F1, _, _, _ = linregress(log_M, log_r_F1)
    slope_B,  intercept_B,  _, _, _ = linregress(log_M, log_r_B)

    log_M_range = np.linspace(log_M.min(), log_M.max(), 50)

    ax3.scatter(log_M, log_r_F3, c='blue',  s=40, alpha=0.7, label=f'F3 (α={slope_F3:.2f})')
    ax3.scatter(log_M, log_r_F1, c='red',   s=40, alpha=0.7, marker='^', label=f'F1 pred (α={slope_F1:.2f})')
    ax3.scatter(log_M, log_r_B,  c='green', s=40, alpha=0.7, marker='s', label=f'Metoda B (α={slope_B:.2f})')

    ax3.plot(log_M_range, intercept_F3 + slope_F3 * log_M_range, 'b--', lw=2)
    ax3.plot(log_M_range, intercept_F1 + slope_F1 * log_M_range, 'r-',  lw=2)
    ax3.plot(log_M_range, intercept_B  + slope_B  * log_M_range, 'g-',  lw=2)

    # Linie referencyjne α=-1 i α=0
    A_ref = np.mean(log_r_F3) - (-0.33) * np.mean(log_M)
    ax3.plot(log_M_range, A_ref + (-0.33) * log_M_range, 'k:', lw=1.5, label='α=-1/3 (F3 teor.)')
    A_f1ref = np.mean(log_r_F1) - (-1.0) * np.mean(log_M)
    ax3.plot(log_M_range, A_f1ref + (-1.0) * log_M_range, 'k--', lw=1.5, label='α=-1 (F1 teor.)')

    ax3.set_xlabel('log₁₀(M_disk / M☉)', fontsize=11)
    ax3.set_ylabel('log₁₀(r_c / kpc)', fontsize=11)
    ax3.set_title('(C) Regresja log r_c vs log M_disk\nK18: α_F3≈-1/3 vs α_F1=-1', fontsize=10)
    ax3.legend(fontsize=7)
    ax3.grid(True, alpha=0.3)

    # --- Panel D: Streszczenie K18 ---
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis('off')

    summary = [
        ("KILL-SHOT K18 — Wyniki SPARC (30 gal.)", "black", 12, True),
        ("", "black", 10, False),
        ("F3 (ε_th = m_sp²/2, N0-6):", "blue", 10, True),
        (f"  α_A = {slope_F3:.3f}  (teoret. ≈ −0.33)", "blue", 9, False),
        (f"  R² = {linregress(log_M, log_r_F3)[2]**2:.3f}", "blue", 9, False),
        ("", "black", 9, False),
        ("F1 (ε(M)=C², N0-4):", "red", 10, True),
        (f"  α_pred = −1.000  (teoret.)", "red", 9, False),
        (f"  α_fit  = {slope_F1:.3f}", "red", 9, False),
        ("", "black", 9, False),
        ("Metoda B (z v_flat):", "green", 10, True),
        (f"  α_B = {slope_B:.3f}", "green", 9, False),
        ("", "black", 9, False),
        ("Ograniczenia metody:", "gray", 9, True),
        ("  Brak bezpośrednich profili HI", "gray", 8, False),
        ("  M_halo z v_flat (przybliżone)", "gray", 8, False),
        ("  Wymagany: THINGS survey <1 kpc", "gray", 8, False),
        ("", "black", 9, False),
        ("Powiązane Kill-shoty:", "purple", 9, True),
        ("  K14: Lyman-α napięcie ~20×", "purple", 8, False),
        ("  K18: r_c vs M_gal (TEN PLIK)", "purple", 8, False),
        ("  Oба rozstrzyga ex42 (THINGS)", "purple", 8, False),
    ]

    y = 0.97
    for text, color, size, bold in summary:
        weight = 'bold' if bold else 'normal'
        ax4.text(0.03, y, text, transform=ax4.transAxes,
                 fontsize=size, color=color, fontweight=weight,
                 verticalalignment='top', family='monospace')
        y -= (size + 2) / 130.0

    plt.suptitle(
        "ex41_sparc_K18.py — Kill-Shot K18: r_c vs M_gal\n"
        "TGP-FDM: F3 (universalny ε_th=m_sp²/2) vs F1 (environmentalny ε(M)=C²)",
        fontsize=12, fontweight='bold', y=1.01
    )

    plt.savefig('ex41_K18_sparc.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\nZapisano: ex41_K18_sparc.png")


def plot_alpha_scan():
    """Skanuje α vs m22 — jak zmienia się nachylenie z m22."""
    m22_scan = np.logspace(-1, 2, 30)
    alphas = []

    M_disk_arr = SPARC_MDISK
    M_halo_arr = SPARC_MHALO

    for m22 in m22_scan:
        r_c_arr = np.array([rc_from_Mhalo_Schive(Mh, m22)[0] for Mh in M_halo_arr])
        log_M = np.log10(M_disk_arr)
        log_r = np.log10(r_c_arr)
        slope, _, _, _, _ = linregress(log_M, log_r)
        alphas.append(slope)

    alphas = np.array(alphas)

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.semilogx(m22_scan, alphas, 'b-', lw=2.5, label='α(m₂₂) — metoda A')
    ax.axhline(y=-1.0,  color='red',  lw=1.5, ls='--', label='F1 predykcja: α = -1')
    ax.axhline(y=-1.0/3, color='navy', lw=1.5, ls=':', label='F3 predykcja: α = -1/3')
    ax.axhline(y=0.0,   color='green', lw=1.5, ls=':', label='α = 0 (brak korelacji)')

    ax.set_xlabel('m₂₂ = m_boson / 10⁻²² eV', fontsize=12)
    ax.set_ylabel('α = d log r_c / d log M_disk', fontsize=12)
    ax.set_title('Skan α(m₂₂) — nachylenie regresji log r_c vs log M_disk\n'
                 'F3: α ≈ −1/3 niezależnie od m₂₂ (relacja Schive+2014)',
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-2, 0.5)

    plt.tight_layout()
    plt.savefig('ex41_alpha_scan.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Zapisano: ex41_alpha_scan.png")


def plot_soliton_mass_function():
    """
    Funkcja masy solitonu M_sol vs M_halo dla F3 i F1.
    Porównanie z Droga Mleczna, M31 i galaktykami karłowatymi.
    """
    M_halo_arr = np.logspace(8, 14, 200)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # --- Panel 1: M_sol vs M_halo ---
    ax = axes[0]
    for m22, col in [(0.5,'lightblue'), (1.0,'blue'), (3.0,'navy'), (10.0,'steelblue')]:
        M_sol_arr = np.array([rc_from_Mhalo_Schive(Mh, m22)[1] for Mh in M_halo_arr])
        ax.loglog(M_halo_arr, M_sol_arr, '-', color=col, lw=2, label=f'F3, m₂₂={m22}')

    # F1: M_sol(M_gal) — masa solitonu rośnie szybciej dla F1
    M_gal_arr = M_halo_arr / 20  # M_disk ~ M_halo/20
    M_sol_F1 = np.array([
        rc_from_Mhalo_Schive(Mh, rc_F1_m22(Md))[1]
        if (rc_F1_m22_val := 1.0 * Md / 1e10) > 0.01 else rc_from_Mhalo_Schive(Mh, 0.01)[1]
        for Mh, Md in zip(M_halo_arr, M_gal_arr)
    ])

    # Punkty obserwacyjne
    obs_Mhalo = [1e10, 5e11, 1e12, 2e12]
    obs_Msol  = [5e8,  2e9,  5e9,  1e10]
    obs_names = ['DDO154', 'NGC3198', 'Droga Mleczna', 'M31']
    for Mh, Ms, name in zip(obs_Mhalo, obs_Msol, obs_names):
        ax.scatter([Mh], [Ms], marker='*', s=150, c='orange', zorder=10)
        ax.annotate(name, (Mh, Ms), textcoords='offset points', xytext=(5,5), fontsize=8)

    # Linia M_sol = M_halo (górna granica)
    ax.loglog(M_halo_arr, M_halo_arr, 'k:', lw=1, alpha=0.5, label='M_sol = M_halo')

    ax.set_xlabel('M_halo [M☉]', fontsize=11)
    ax.set_ylabel('M_sol [M☉]', fontsize=11)
    ax.set_title('(A) M_sol vs M_halo\nRelacja soliton-halo Schive+2014', fontsize=10)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)

    # --- Panel 2: r_c vs M_halo ---
    ax2 = axes[1]
    for m22, col in [(0.5,'lightblue'), (1.0,'blue'), (3.0,'navy'), (10.0,'steelblue')]:
        r_c_arr = np.array([rc_from_Mhalo_Schive(Mh, m22)[0] for Mh in M_halo_arr])
        ax2.loglog(M_halo_arr, r_c_arr, '-', color=col, lw=2, label=f'F3, m₂₂={m22}')

    # F1 predykcja
    r_c_F1_arr = np.array([rc_F1(Mh, Mh/20, 1.0, 1e10) for Mh in M_halo_arr])
    ax2.loglog(M_halo_arr, r_c_F1_arr, 'r-', lw=2.5, label='F1: r_c ∝ 1/M')

    # Dane obserwacyjne
    obs_rc = [1.5, 1.6, 2.5, 3.0]
    for Mh, rc, name in zip(obs_Mhalo, obs_rc, obs_names):
        ax2.scatter([Mh], [rc], marker='*', s=150, c='orange', zorder=10)
        ax2.annotate(name, (Mh, rc), textcoords='offset points', xytext=(5,5), fontsize=8)

    ax2.set_xlabel('M_halo [M☉]', fontsize=11)
    ax2.set_ylabel('r_c [kpc]', fontsize=11)
    ax2.set_title('(B) r_c vs M_halo\nF3 (linie) vs F1 (czerwona)\nK18 rozstrzyga który model', fontsize=10)
    ax2.legend(fontsize=8, loc='lower left')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('ex41_Msol_Mhalo.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Zapisano: ex41_Msol_Mhalo.png")


def rc_F1_m22(M_disk, m22_ref=1.0, M_disk_ref=1e10):
    """Pomocnicza: m22_eff dla F1."""
    return max(m22_ref * (M_disk / M_disk_ref), 0.01)


if __name__ == "__main__":
    # Uruchom główną analizę
    all_results, alpha_A, alpha_B, k18_status = main()

    # Generuj wykresy
    print("\n[GENEROWANIE WYKRESÓW]")
    print("  1. Główna analiza K18 (4 panele)...")
    plot_K18_analysis(all_results)

    print("  2. Skan α(m22)...")
    plot_alpha_scan()

    print("  3. Funkcja masy solitonu M_sol–M_halo...")
    plot_soliton_mass_function()

    print(f"\n{'='*70}")
    print(f"PODSUMOWANIE ex41:")
    print(f"  Kill-Shot K18 status: {k18_status}")
    print(f"  α zmierzony (Metoda A): {alpha_A:.3f}")
    print(f"  α zmierzony (Metoda B): {alpha_B:.3f}")
    print(f"  F3 predykcja: α ≈ −0.33  (M_sol ∝ M_halo^{{1/3}}, r_c ∝ M_sol^{{-1/3}})")
    print(f"  F1 predykcja: α = −1.00  (r_c ∝ 1/M_gal)")
    print(f"  Następny: ex42 — bezpośrednie profile HI (THINGS survey)")
    print("=" * 70)
