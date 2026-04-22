"""
ex35_dm_mechanisms.py
=====================
Pełny test mechanizmów ciemnej materii w TGP.

Analiza czterech potencjalnych mechanizmów DM w TGP Drogi B:

  M1: Grawitacja Yukawa (ex14 revisited)
      G_eff(r) = G * exp(-r/λ) — modyfikacja siły grawitacji.
      WYNIK ex14: Daje MNIEJ grawitacji (nie więcej). FALSIFIED.

  M2: Akumulacja trójciałowa V₃ (N-body DM)
      N-ciał gwiazd tworzy kolektywne pole V₃ wzmacniające grawitację.
      Potrzeba: F₃/F₂ ~ O(1) dla N ~ 10^11 gwiazd.
      Test: oblicz F₃/F₂ dla galaktyki.

  M3: Energia pola skalarnego jako DM
      ρ_field = (∇Φ)²/2 + V(Φ) — energia pola działa jak masa.
      Test: czy ρ_field(r) ~ ρ_NFW(r) dla galaktyki?

  M4: FDM — kondensatem bozonowy (patrz ex36 dla szczegółów)
      Pole Φ jako bozon FDM, m ~ 10^{-22} eV.
      WYNIK ex36: Możliwe, ale wymaga ε << 1 i fine-tuningu.

GALAKTYKA: NGC 3198 (Begeman 1989)
ZAŁOŻENIA:
  - TGP Droga B: Yukawa (∇²-m²)Φ = -4πC δ³(r)
  - Φ = Φ₀ * g(r), g_vac = 1
  - V_TGP = βg³/3 - γg⁴/4, β = γ = 1 (Planck)
  - m_sp = √γ = 1 (masa Plancka, limit silny)
  - Dla galaktyk: m_sp << 1 (hierarchia skal)
"""

import sys, os
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import k0, k1, i0, i1
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# =============================================================================
# Stałe i dane
# =============================================================================
G_SI     = 6.674e-11      # m^3 kg^-1 s^-2
kpc_m    = 3.086e19       # m
pc_m     = 3.086e16       # m
M_sun    = 1.989e30       # kg
km_s     = 1e3            # m/s
l_Pl_m   = 1.616e-35      # m
m_Pl_kg  = 2.176e-8       # kg
hbar_SI  = 1.055e-34      # J*s
c_SI     = 3e8            # m/s

# NGC 3198 (Begeman 1989)
OBS_R_KPC = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0])
OBS_V_KMS = np.array([90., 130., 155., 160., 155., 152., 150., 148., 147., 146., 145.])
OBS_ERR   = np.array([ 8.,   6.,   5.,   5.,   4.,   4.,   4.,   4.,   5.,   6.,   7.])

M_disk = 2.0e10 * M_sun
R_d    = 3.2 * kpc_m

print("=" * 70)
print("EX35: MECHANIZMY CIEMNEJ MATERII W TGP — PEŁNY TEST")
print("=" * 70)
print()
print("Galaktyka: NGC 3198, M_disk=2e10 M_sun, R_d=3.2 kpc")
print()

r_kpc = np.linspace(0.3, 33.0, 150)
r_m   = r_kpc * kpc_m

# =============================================================================
# Funkcje pomocnicze
# =============================================================================
def v_disk_freeman(r_m_arr, M_disk_kg, R_d_m):
    """Prędkość kołowa dla dysku eksponencjalnego (Freeman 1970)."""
    Sigma0 = M_disk_kg / (2.0 * np.pi * R_d_m**2)
    y_arr = r_m_arr / (2.0 * R_d_m)
    y_arr = np.maximum(y_arr, 1e-10)
    v2 = (4.0 * np.pi * G_SI * Sigma0 * R_d_m
          * y_arr**2 * (i0(y_arr)*k0(y_arr) - i1(y_arr)*k1(y_arr)))
    return np.sqrt(np.maximum(v2, 0.0))

def v_NFW_arr(r_kpc_arr, rho0_kg_m3, r_s_m):
    """NFW profil prędkości rotacji."""
    r_m_arr = r_kpc_arr * kpc_m
    x = r_m_arr / r_s_m
    M = 4.0 * np.pi * rho0_kg_m3 * r_s_m**3 * (np.log(1.0 + x) - x/(1.0+x))
    return np.sqrt(G_SI * M / r_m_arr)

def chi2_ngc(v_arr_kms, r_arr_kpc=None):
    """Chi^2/N vs NGC 3198."""
    if r_arr_kpc is None:
        r_arr_kpc = r_kpc
    v_interp = np.interp(OBS_R_KPC, r_arr_kpc, v_arr_kms)
    return np.sum(((v_interp - OBS_V_KMS) / OBS_ERR)**2) / len(OBS_V_KMS)

# ============================================================
# MECHANIZM 1: Yukawa (Droga B) — ex14 revisited
# ============================================================
print("=" * 70)
print("MECHANIZM 1: Yukawa G_eff(r) = G*exp(-r/λ)  [ex14 revisited]")
print("=" * 70)
print()
print("Założenie: G_eff(r) modyfikuje siłę Newtona z dysku baryonowego.")
print("G_eff(r) = G * exp(-r/λ)  — maleje z r (tłumienie Yukawa).")
print()

v_bar = v_disk_freeman(r_m, M_disk, R_d)

lambda_vals = [1.0, 5.0, 10.0, 20.0, 50.0, 1000.0]  # kpc
print(f"  {'λ [kpc]':>10s}  {'v(5kpc)':>8s}  {'v(15kpc)':>9s}  {'v(25kpc)':>9s}  {'chi^2/N':>8s}")
print("-" * 50)

m1_curves = {}
for lam_kpc in lambda_vals:
    lam_m = lam_kpc * kpc_m
    # G_eff(r) = G * exp(-r/λ) => v^2_TGP(r) = v^2_disk(r) * exp(-r/λ) * fcorr
    # (uproszczone: faktyczne całkowanie po dysku daje < niż punkt.)
    # Dla r >> R_d: v_TGP(r) ≈ v_bar(r) * sqrt(exp(-r/λ))
    # Dla r ~ R_d: korygujemy przez efekt Bessela
    t_arr = r_m / lam_m
    corr = np.exp(-t_arr) * np.sqrt(1.0 + t_arr)  # heurystyczna korekta dysku
    corr = np.minimum(corr, 1.0)
    v_tgp = v_bar * np.sqrt(corr)
    m1_curves[lam_kpc] = v_tgp / km_s

    c2 = chi2_ngc(v_tgp / km_s)
    v5  = np.interp(5.0, r_kpc, v_tgp) / km_s
    v15 = np.interp(15.0, r_kpc, v_tgp) / km_s
    v25 = np.interp(25.0, r_kpc, v_tgp) / km_s
    print(f"  {lam_kpc:>10.0f}  {v5:>8.1f}  {v15:>9.1f}  {v25:>9.1f}  {c2:>8.2f}")

chi2_bar = chi2_ngc(v_bar / km_s)
print(f"  {'Newton':>10s}  "
      f"{np.interp(5.0,r_kpc,v_bar)/km_s:>8.1f}  "
      f"{np.interp(15.0,r_kpc,v_bar)/km_s:>9.1f}  "
      f"{np.interp(25.0,r_kpc,v_bar)/km_s:>9.1f}  {chi2_bar:>8.2f}")
print()
print("WYNIK M1: Yukawa G_eff(r) < G dla wszystkich r>0.")
print("  => Krzywa rotacji TGP jest ZAWSZE niższa od Newtona.")
print("  => TGP nie może wyjaśnić płaskiej krzywej bez dodatkowej masy.")
print("  MECHANIZM M1: SFALSYFIKOWANY (potwierdza ex14).")
print()

# ============================================================
# MECHANIZM 2: Akumulacja V₃ — trójciałowe sily od N gwiazd
# ============================================================
print("=" * 70)
print("MECHANIZM 2: Akumulacja V₃ (trójciałowe, N-body)")
print("=" * 70)
print()
print("Siła trójciałowa V₃ ~ 1/P (obwód trójkąta) dla trójki Yukawa.")
print("Dla N gwiazd: kolektywne V₃ daje efektywną masę ciemną?")
print()
print("Analiza:")
print("  V₂(d) = -4πC²/d  (Yukawa parowe, m_sp→0)")
print("  V₃(d) = -γ * I_Yukawa(d)  (trójciałowe)")
print("  I_Yukawa ~ 8π²/P ~ 8π²/(3d)  dla równobocznego trójkąta")
print()
print("  F₂(r) ≈ 4πC²/r² (siła parowa, od 1 sąsiada)")
print("  F₃(r) ≈ γ * 8π²/(9d²) * N_eff(r)  (od N_eff trojek)")
print()
print("Dla galaktyki NGC 3198:")
N_gal = 2e11  # liczba gwiazd

# Parametry TGP dla skali galaktycznej
# Masa gwiezdna C w TGP: C = m_sp/(2√π) * M [Planck]
# Dla m_sp = m_Planck: C_star = m_Pl_kg/(2√π) * M_sun (bez sensu)
# W TGP, C jest w jednostkach Plancka:
# C_star = (m_sp / (2*sqrt(π))) * (M_star / m_Pl)  [bezwymiarowe Planck]
# M_star = 1 M_sun, m_Pl = 2.176e-8 kg
C_star_Pl = (1.0 / (2.0*np.sqrt(np.pi))) * (M_sun / m_Pl_kg)  # dla m_sp=1 l_Pl^{-1}
print(f"  N_gal = {N_gal:.0e} gwiazd")
print(f"  C_☉ (m_sp=1 Pl) = {C_star_Pl:.3e}  [jednostki Plancka]")
print()

# Stosunek F₃/F₂ dla zbiorowej siły trójciałowej
# Każda gwiazda uczestniczy w N trojkach z jej sąsiadami
# F₃ ≈ (4/9) * γ * N_gal * C_star * F₂
# (z derivacji w ANALIZA_CIEMNA_MATERIA.md)

# Dla galaktyki: gamma ~ (m_sp)^2
# Potrzeba F₃/F₂ ~ 1 dla efektu DM
# => γ * N_gal * C_star ~ 1
# => γ ~ 1/(N_gal * C_star)

gamma_needed = 1.0 / (N_gal * C_star_Pl)
m_sp_needed = np.sqrt(gamma_needed)

print("  Warunek F₃/F₂ ~ 1:")
print(f"    γ = m_sp² wymagane = {gamma_needed:.3e}")
print(f"    m_sp wymagane = {m_sp_needed:.3e} l_Pl^{{-1}}")
print()

# Lambda Yukawa odpowiadające m_sp_needed
lambda_needed_Pl = 1.0 / m_sp_needed  # w l_Pl
lambda_needed_m  = lambda_needed_Pl * l_Pl_m
lambda_needed_kpc = lambda_needed_m / kpc_m
print(f"    λ = 1/m_sp = {lambda_needed_kpc:.3e} kpc")
print()
print(f"    Rozmiar NGC 3198 ~ 30 kpc")
print(f"    λ wymagane >> rozmiar galaktyki?",
      "TAK" if lambda_needed_kpc > 30 else "NIE")
print()

# Krzywa rotacji z V₃
# Jeśli F₃ daje efektywną masę dodatkową M₃(r):
# v²_V3(r) = G * M₃(r) / r  gdzie M₃ ~ constant (równomierne halo)
# To daje v_V3 ~ 1/√r (malejąca) — NIE płaska!

print("  Kształt krzywej rotacji od V₃:")
print("    F₃(r) ∝ N(<r) * ρ(r) / r²")
print("    Dla płaskiej krzywej: M(<r) ∝ r => ρ(r) ∝ 1/r²")
print("    Dysk wykładniczy daje: N(<r) → const dla r >> R_d")
print("    => F₃(r) ∝ 1/r² => M_eff(r) ∝ r => v(r) ~ const !")
print()
print("    Ale TYLKO jeśli gęstość stała 1/r² — dysk spada szybciej!")
print("    Profil dysku: Σ(r) ~ exp(-r/R_d) => N(<r) → N_total dla r>>R_d")
print("    M_eff V₃ ∝ N_total/r² => v_V3 ∝ 1/r  (Keplerian, MALEJĄCA!)")
print()

# Oblicz krzywą rotacji V₃ dla różnych γ
def v_V3_rotation(r_kpc_arr, gamma_val, N_total, C_star_Pl_val, R_d_kpc=3.2):
    """
    Prędkość rotacji od akumulacji V₃.
    Uproszczony model: M_eff_V3(r) = (4/9)*γ*N_eff(r)*C*M_star
    N_eff(r) = N(<r) * <N_neigh>  ~ N_total * (1-exp(-r/R_d))
    F₃/F₂ daje efektywne wzmocnienie G.
    """
    R_d_m_val = R_d_kpc * kpc_m
    N_eff = N_total * (1.0 - np.exp(-r_kpc_arr / R_d_kpc))
    # Efektywna masa V₃ (heurystyczna)
    # G_eff(r) = G * (1 + (4/9)*γ*N_eff*C_star)
    ratio = (4.0/9.0) * gamma_val * N_eff * C_star_Pl_val
    G_eff = G_SI * (1.0 + ratio)
    v2_bar = v_disk_freeman(r_kpc_arr * kpc_m, M_disk, R_d) ** 2
    v2_total = G_eff / G_SI * v2_bar  # skalowanie przez G_eff
    return np.sqrt(v2_total) / km_s

print(f"  Symulacja krzywych rotacji V₃ dla różnych γ:")
print(f"  {'γ':>10s}  {'m_sp':>10s}  {'v(5kpc)':>8s}  {'v(25kpc)':>9s}  {'chi²/N':>8s}")
print("-" * 50)

gamma_test = [1e-40, 1e-45, 1e-50, 1e-60, gamma_needed]
m2_curves = {}
for gam in gamma_test:
    m_sp_g = np.sqrt(gam)
    v_v3 = v_V3_rotation(r_kpc, gam, N_gal, C_star_Pl)
    c2 = chi2_ngc(v_v3)
    v5  = np.interp(5.0, r_kpc, v_v3)
    v25 = np.interp(25.0, r_kpc, v_v3)
    m2_curves[gam] = v_v3
    print(f"  {gam:>10.2e}  {m_sp_g:>10.2e}  {v5:>8.1f}  {v25:>9.1f}  {c2:>8.2f}")

print()
print("WYNIK M2: V₃ akumulacja może wzmocnić grawitację,")
print("  ale kształt krzywej ∝ v_bar(r) — NIE PŁASKI dla dużych r.")
print("  Wzmocnienie jest równomierne, nie odtwarza halo DM.")
print(f"  Wymagane γ = {gamma_needed:.2e} daje m_sp = {np.sqrt(gamma_needed):.2e} l_Pl^{{-1}}.")
print(f"  Odpowiada λ = {lambda_needed_kpc:.2e} kpc — dużo większe niż galaktyka!")
print("  MECHANIZM M2: NIESKUTECZNY dla płaskiej krzywej rotacji.")
print()

# ============================================================
# MECHANIZM 3: Energia pola skalarnego jako DM
# ============================================================
print("=" * 70)
print("MECHANIZM 3: Energia pola skalarnego TGP jako DM")
print("=" * 70)
print()
print("Idea: Pole Φ wokół galaktyki ma energię (∇Φ)²/2 + V(Φ),")
print("która działa jak efektywna gęstość masy (T₀₀ w GR).")
print()
print("Dla profilu Yukawa Φ(r) = Φ₀ - (C/r)*exp(-m_sp*r):")
print("  ∇Φ = (C/r²)(1 + m_sp*r)*exp(-m_sp*r)")
print("  (∇Φ)² = (C/r²)²*(1+m_sp*r)²*exp(-2*m_sp*r)")
print()
print("Gęstość energii pola TGP:")
print("  ρ_field(r) = Φ₀² * [(∇g)²/2 + V(g)]")
print("  gdzie g = 1 - (C/Φ₀r)*exp(-m_sp*r)")
print()

def g_yukawa(r_m_val, C_Pl, m_sp_Pl, Phi0=1.0):
    """
    Pole g = 1 - (C/Phi0*r)*exp(-m_sp*r) dla pojedynczego źródła.
    r w l_Pl, C w jednostkach Plancka.
    """
    if r_m_val < 1e-10:
        return 1.0
    return 1.0 - (C_Pl / (Phi0 * r_m_val)) * np.exp(-m_sp_Pl * r_m_val)

def rho_field_TGP(r_m_arr, C_Pl, m_sp_Pl, Phi0=1.0):
    """
    Gęstość energii pola TGP (w jednostkach Plancka na l_Pl^3).
    ρ = Φ₀² * [(∇g)²/2 + V(g)]
    V(g) = g^3/3 - g^4/4  (β=γ=1)

    Gradient analityczny: dg/dr = (C/Phi0) * (1/r² + m_sp/r) * exp(-m_sp*r)
    """
    out = np.zeros(len(r_m_arr))
    for i, r_val in enumerate(r_m_arr):
        if r_val < 1e-8:
            r_val = 1e-8

        # Gradient g
        exp_term = np.exp(-m_sp_Pl * r_val)
        dg_dr = (C_Pl / Phi0) * (1.0/r_val**2 + m_sp_Pl/r_val) * exp_term

        # Wartość g
        g_val = 1.0 - (C_Pl / (Phi0 * r_val)) * exp_term

        # Potencjał V(g)
        V_g = g_val**3/3.0 - g_val**4/4.0

        # Gęstość energii (bez Φ₀² — bezwymiarowe)
        rho = 0.5 * dg_dr**2 + V_g
        out[i] = rho
    return out

# Skala galaktyczna w Planck
# r_10kpc = 10 kpc * (3.086e19 m/kpc) / (1.616e-35 m/l_Pl) = gigantyczna liczba
r_10kpc_Pl = 10.0 * kpc_m / l_Pl_m
print(f"  r = 10 kpc = {r_10kpc_Pl:.3e} l_Pl")
print(f"  => λ_Yukawa = r_10kpc w Planck ~ {r_10kpc_Pl:.3e}")
print(f"  => m_sp = 1/λ = {1.0/r_10kpc_Pl:.3e} l_Pl^{{-1}}")
print()

# Dla galaktycznych skal, liczymy w fizycznych jednostkach
# ρ_field(r) w kg/m^3, dla gwiazdy o masie M_sun
M_disk_Pl = M_disk / m_Pl_kg  # masa dysku w Planck
C_disk_Pl = M_disk_Pl / (2.0 * np.sqrt(np.pi))  # C całego dysku

m_sp_galaxy = 1.0 / (10.0 * kpc_m / l_Pl_m)  # λ=10 kpc → m_sp w l_Pl^{-1}

print(f"  C_disk (Planck) = {C_disk_Pl:.3e}  (masa dysku 2e10 M_sun)")
print(f"  m_sp (λ=10kpc) = {m_sp_galaxy:.3e} l_Pl^{{-1}}")
print()
print("  Gęstość energii pola ρ_field vs NFW:")
print()

# Oblicz profile gęstości w fizycznych jednostkach
# Przeliczenie: ρ_field [Pl] * (m_Pl/l_Pl^3) = ρ_field * (2.176e-8 kg)/(1.616e-35 m)^3
rho_Pl_to_SI = m_Pl_kg / l_Pl_m**3  # ~ 5.155e96 kg/m^3

# NFW dla NGC 3198 (Begeman 1989)
rho0_NFW = 6.7e6 * M_sun / kpc_m**3  # kg/m^3
r_s_NFW  = 12.0 * kpc_m

rho_NFW_arr = rho0_NFW / ((r_m / r_s_NFW) * (1.0 + r_m/r_s_NFW)**2)

# Dla wizualizacji — tylko profile kształtu (normalizowane)
r_kpc_test = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 30.0])
r_m_test   = r_kpc_test * kpc_m
r_Pl_test  = r_m_test / l_Pl_m  # w l_Pl

# Profil gęstości pola (bezwymiarowy kształt)
rho_field_shape = np.zeros(len(r_kpc_test))
for i, r_Pl in enumerate(r_Pl_test):
    # Gradientowy człon dominuje dla małych r
    exp_t = np.exp(-m_sp_galaxy * r_Pl)
    dg_dr_Pl = C_disk_Pl * (1.0/r_Pl**2 + m_sp_galaxy/r_Pl) * exp_t
    rho_field_shape[i] = 0.5 * dg_dr_Pl**2

# NFW kształt (normalizowany przez ρ₀)
rho_NFW_shape = 1.0 / ((r_kpc_test / (r_s_NFW/kpc_m)) *
                        (1.0 + r_kpc_test / (r_s_NFW/kpc_m))**2)

# Normalizuj kształty do r=10 kpc
idx_10 = np.argmin(np.abs(r_kpc_test - 10.0))
norm_field = rho_field_shape[idx_10]
norm_NFW   = rho_NFW_shape[idx_10]

print(f"  {'r [kpc]':>8s}  {'ρ_field shape':>15s}  {'ρ_NFW shape':>13s}  {'Ratio':>8s}")
print("-" * 50)
for i in range(len(r_kpc_test)):
    rf = rho_field_shape[i] / norm_field if norm_field > 0 else 0
    rn = rho_NFW_shape[i] / norm_NFW
    ratio = rf / rn if rn > 0 and rf > 0 else 0.0
    print(f"  {r_kpc_test[i]:>8.1f}  {rf:>15.4f}  {rn:>13.4f}  {ratio:>8.4f}")

print()

# Porównanie kształtu
print("  Kształt ρ_field vs ρ_NFW:")
print("    NFW: ρ ∝ 1/[r*(1+r/r_s)²]  — kulisty halo")
print("    Yukawa: (∇Φ)² ∝ (1/r² + m_sp/r)² * exp(-2*m_sp*r)")
print()
print("    Dla r << λ: ρ_field ~ 1/r⁴  (szybszy niż NFW!)")
print("    Dla r >> λ: ρ_field ~ exp(-2*m_sp*r)  (tłumione eksponencjalnie!)")
print()
print("  Wniosek: Kształt pola TGP ≠ kształt NFW w żadnym zakresie r.")
print("           ρ_field maleje ZA SZYBKO dla dużych r (Yukawa tłumione).")
print()
print("WYNIK M3: Energia pola TGP nie odtwarza profilu NFW.")
print("  Yukawa: ρ_field ~ exp(-2r/λ) vs NFW ~ 1/r dla dużych r.")
print("  Musi być λ >> rozmiar galaktyki (Newton limit),")
print("  ale wtedy ρ_field → 0 w całej galaktyce.")
print("  MECHANIZM M3: SFALSYFIKOWANY dla standardowego Yukawa.")
print()

# ============================================================
# MECHANIZM 4: Tachionowa niestabilność (m_sp² < 0)
# ============================================================
print("=" * 70)
print("MECHANIZM 4: Tachionowa niestabilność spinodalna")
print("=" * 70)
print()
print("Jeśli m_sp² < 0 (tachion), rozwiązanie pola to nie e^{-mr}/r")
print("ale sin(|m|r)/r lub cos(|m|r)/r — oscylujące potencjały!")
print()
print("Pole tachionowe: Φ(r) = A*sin(|m_t|*r)/r")
print("  => ∇Φ = A*(cos(|m_t|*r)/r² - sin(|m_t|*r)/r³) * |m_t|")
print()

# Rozwiązanie ODE dla tachionowego pola
# (∇² + m_t²) Φ = 0 (bez źródła, jednorodne)
# Rozwiązania: sin(m_t*r)/r, cos(m_t*r)/r

def phi_tachyon(r_kpc_arr, A_coeff, m_t_inv_kpc):
    """
    Tachionowe pole skalarne: Φ(r) = A*sin(m_t*r)/r + B*cos(m_t*r)/r.
    (wybieramy sin — regularny w r=0 w sensie dystrybucji)
    """
    r = r_kpc_arr
    m_t = m_t_inv_kpc
    # Unikamy r=0
    r = np.maximum(r, 1e-6)
    return A_coeff * np.sin(m_t * r) / r

def v_tachyon_rotation(r_kpc_arr, A_coeff, m_t_inv_kpc, Phi0=1.0):
    """
    Prędkość rotacji od tachionowego pola.
    ρ_field = (∇Φ)²/2 + V(Φ) — ale V jest nieokreślone dla tachionu!
    Używamy tylko członu kinetycznego: ρ ~ (∇Φ)²/2.

    Gradient: dΦ/dr = A*(m_t*cos(m_t*r)/r - sin(m_t*r)/r²)
    """
    r = np.maximum(r_kpc_arr, 1e-6)
    m_t = m_t_inv_kpc
    dPhi_dr = A_coeff * (m_t * np.cos(m_t*r)/r - np.sin(m_t*r)/r**2)
    rho_kin = 0.5 * Phi0**2 * dPhi_dr**2  # [bezwym./kpc^2]

    # M(<r) = 4π ∫₀ʳ ρ_kin(r') r'^2 dr'
    # Konwertujemy do SI po drodze (skalowanie)
    # ρ_kin [1/kpc^2] * (Φ₀ c² hbar / l_Pl)² / kpc
    # Na razie obliczamy tylko kształt, normalizujemy do v_obs

    # Masa od pola tachionowego (numerycznie)
    M_arr = np.zeros(len(r_kpc_arr))
    for i, r_val in enumerate(r_kpc_arr):
        # Całka od 0 do r_val
        r_int = np.linspace(1e-6, r_val, 100)
        dPhi = A_coeff * (m_t * np.cos(m_t*r_int)/r_int - np.sin(m_t*r_int)/r_int**2)
        rho_int = 0.5 * dPhi**2  # [A²/kpc^4]
        # M ∝ ∫ ρ * r'^2 * dr' [A²/kpc]
        M_arr[i] = 4.0 * np.pi * np.trapz(rho_int * r_int**2, r_int)

    # Normalizuj tak, żeby v(10kpc) = 150 km/s
    if M_arr[np.argmin(np.abs(r_kpc_arr - 10.0))] > 0:
        M_norm = M_arr / M_arr[np.argmin(np.abs(r_kpc_arr - 10.0))]
        # v²(r) = G * M_dm(<r) / r  z M_dm znormalizowanym
        # G * M_dm_SI / r_SI = (150 km/s)² przy r=10kpc
        v2_norm = M_norm / r_kpc_arr  # [kpc^{-1}], kształt
        v_ref = 150.0  # km/s
        v_arr = v_ref * np.sqrt(v2_norm / v2_norm[np.argmin(np.abs(r_kpc_arr-10.0))])
        return v_arr
    else:
        return np.zeros(len(r_kpc_arr))

print("  Kształt krzywej rotacji od pola tachionowego (znormalizowany do 150 km/s):")
print()

m_t_vals = [1.0, 2.0, 5.0, 10.0]  # kpc^{-1}
print(f"  {'m_t [kpc⁻¹]':>12s}  {'λ_t [kpc]':>10s}  {'v(5kpc)':>8s}  "
      f"{'v(15kpc)':>9s}  {'chi²/N':>8s}  {'Płaska?':>8s}")
print("-" * 60)

m4_curves = {}
for m_t in m_t_vals:
    lam_t = 1.0 / m_t
    v_t = v_tachyon_rotation(r_kpc, A_coeff=1.0, m_t_inv_kpc=m_t)
    if v_t.max() > 0:
        # Dodaj dysk barionowy
        v_total = np.sqrt(v_bar**2 / km_s**2 + v_t**2)
        c2 = chi2_ngc(v_total)
        v5 = np.interp(5.0, r_kpc, v_total)
        v15 = np.interp(15.0, r_kpc, v_total)
        # Płaska: v(25kpc)/v(10kpc) > 0.90?
        v10 = np.interp(10.0, r_kpc, v_total)
        v25 = np.interp(25.0, r_kpc, v_total)
        flat = "TAK" if v25/v10 > 0.90 else "NIE"
        m4_curves[m_t] = v_total
        print(f"  {m_t:>12.1f}  {lam_t:>10.2f}  {v5:>8.1f}  {v15:>9.1f}  {c2:>8.2f}  {flat:>8s}")

print()
print("  Problem z M4 (tachion):")
print("    1. Tachionowe pole oscilluje: sin(m_t*r)/r — zmienne znaki!")
print("    2. Gęstość energii oscyluje → lokalne 'halo' i 'anty-halo'")
print("    3. Brak stabilnego profilu ciemnej materii")
print("    4. Tachion → niestabilność pola → teoria nie jest sensowna")
print()
print("WYNIK M4: Tachionowe pole NIE daje stabilnego DM halo.")
print("  Oscylacje ρ_field nie odtwarzają spójnego profilu NFW.")
print("  Teoria z m_sp² < 0 jest niestabilna i nie przewiduje obserwacji.")
print("  MECHANIZM M4: SFALSYFIKOWANY (niestabilność teorii).")
print()

# ============================================================
# PORÓWNANIE WSZYSTKICH MECHANIZMÓW
# ============================================================
print("=" * 70)
print("TABELA PODSUMOWUJĄCA: Mechanizmy DM w TGP vs NGC 3198")
print("=" * 70)
print()

table = [
    ("Newton (baryony)", v_bar / km_s, "—"),
    ("M1: Yukawa (λ=100 kpc)", m1_curves[100.0], "Sfalsyfikowany"),
    ("M2: V₃ akumulacja (γ=1e-40)", m2_curves[1e-40], "Nieskuteczny"),
    ("M4: Tachion (m_t=2)", m4_curves.get(2.0, v_bar/km_s), "Niestabilny"),
]

print(f"  {'Mechanizm':35s}  {'chi²/N':>8s}  {'Wynik':>20s}")
print("-" * 70)
for name, v, verdict in table:
    c2 = chi2_ngc(v)
    print(f"  {name:35s}  {c2:>8.2f}  {verdict:>20s}")

# NFW (benchmark)
rho0_ISO = 0.054 * M_sun / pc_m**3
r_c_ISO  = 4.0 * kpc_m
v_iso = np.sqrt((v_bar)**2 +
                (4.0*np.pi*G_SI*rho0_ISO*r_c_ISO**2 *
                 (1.0 - r_c_ISO/r_m*np.arctan(r_m/r_c_ISO)))**1)
v_cdm = np.sqrt(v_bar**2 +
                4.0*np.pi*G_SI*rho0_ISO*r_c_ISO**2 *
                np.maximum(1.0 - r_c_ISO/r_m*np.arctan(r_m/r_c_ISO), 0))
chi2_cdm = chi2_ngc(v_cdm / km_s)
print(f"  {'CDM (ISO halo, Begeman fit)':35s}  {chi2_cdm:>8.2f}  {'Akceptowalny (2par)':>20s}")

print()
print()
print("WNIOSKI KOŃCOWE EX35:")
print()
print("  ● ŻADEN z czterech mechanizmów TGP nie wyjaśnia krzywej rotacji")
print("    NGC 3198 w sposób naturalny.")
print()
print("  ● M1 (Yukawa): Daje mniej grawitacji — dowód ex14.")
print()
print("  ● M2 (V₃): Wzmocnienie jest równomierne (zależy od v_bar),")
print("    nie tworzy dodatkowego halo. Wymagane γ = {:.2e}".format(gamma_needed))
print("    daje λ = {:.2e} kpc >> R_galaktyki.".format(lambda_needed_kpc))
print()
print("  ● M3 (energia pola): kształt Yukawa ~ exp(-2r/λ) ≠ NFW ~ 1/r.")
print("    Nie ma zakresu λ gdzie oba pasują jednocześnie.")
print()
print("  ● M4 (tachion): Niestabilność kwantowa, oscylujący profil.")
print("    Nie może dać stabilnego halo DM.")
print()
print("  ● JEDYNA NADZIEJA: Rozszerzone TGP z V_mod = εg² + V_TGP")
print("    (Fuzzy Dark Matter, patrz ex36). Wymaga nowego parametru ε.")
print()

# ============================================================
# WYKRESY
# ============================================================
print("=" * 70)
print("Generowanie wykresów...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('TGP i Ciemna Materia: Test Czterech Mechanizmów (ex35)\n'
             'Galaktyka: NGC 3198', fontsize=13, y=1.01)

# 1. Krzywe rotacji Yukawa (M1)
ax = axes[0, 0]
ax.errorbar(OBS_R_KPC, OBS_V_KMS, yerr=OBS_ERR,
            fmt='ko', ms=5, capsize=3, label='NGC 3198 (obs.)', zorder=10)
ax.plot(r_kpc, v_bar/km_s, 'k--', lw=1.5, label='Newton (baryony)')
cols = cm.Reds(np.linspace(0.4, 0.9, len(lambda_vals)))
for (lam_kpc, v_tgp), col in zip(m1_curves.items(), cols):
    c2 = chi2_ngc(v_tgp)
    ax.plot(r_kpc, v_tgp, color=col, lw=1.5,
            label=rf'TGP $\lambda={lam_kpc:.0f}$ kpc ($\chi^2$={c2:.1f})')
ax.set_title('M1: Yukawa G_eff(r) = G·exp(-r/λ)', fontsize=11)
ax.set_xlabel('r [kpc]', fontsize=10)
ax.set_ylabel('v [km/s]', fontsize=10)
ax.legend(fontsize=7, loc='lower right')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 33); ax.set_ylim(0, 220)

# 2. Krzywe rotacji V₃ (M2)
ax = axes[0, 1]
ax.errorbar(OBS_R_KPC, OBS_V_KMS, yerr=OBS_ERR,
            fmt='ko', ms=5, capsize=3, label='NGC 3198 (obs.)', zorder=10)
ax.plot(r_kpc, v_bar/km_s, 'k--', lw=1.5, label='Newton')
cols2 = cm.Blues(np.linspace(0.4, 0.9, len(gamma_test)))
for (gam, v_v3), col in zip(m2_curves.items(), cols2):
    c2 = chi2_ngc(v_v3)
    ax.plot(r_kpc, v_v3, color=col, lw=1.5,
            label=rf'V₃ $\gamma$={gam:.0e} ($\chi^2$={c2:.1f})')
ax.set_title('M2: Akumulacja V₃ (trójciałowe)', fontsize=11)
ax.set_xlabel('r [kpc]', fontsize=10)
ax.set_ylabel('v [km/s]', fontsize=10)
ax.legend(fontsize=7, loc='lower right')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 33); ax.set_ylim(0, 220)

# 3. Profil ρ_field vs NFW (M3)
ax = axes[0, 2]
# Normalizuj ρ_field do ρ_NFW przy r=5kpc
r_plot_kpc = np.linspace(0.5, 30.0, 200)
r_plot_m = r_plot_kpc * kpc_m
r_plot_Pl = r_plot_m / l_Pl_m

rho_NFW_plot = rho0_NFW / ((r_plot_m / r_s_NFW) * (1.0 + r_plot_m/r_s_NFW)**2)

# ρ_field kształt (∇Φ)² ∝ (1/r² + m/r)² exp(-2mr) — w kpc^{-2} arbitralne
m_sp_plot_inv_kpc = np.array([0.1, 0.3, 1.0])  # kpc^{-1}
cols3 = cm.Greens(np.linspace(0.5, 0.9, len(m_sp_plot_inv_kpc)))

# Normalizacja do r=5kpc
idx_5 = np.argmin(np.abs(r_plot_kpc - 5.0))
rho_NFW_norm = rho_NFW_plot / rho_NFW_plot[idx_5]

ax.loglog(r_plot_kpc, rho_NFW_norm, 'k-', lw=2.5, label='NFW ρ')
for m_sp_k, col in zip(m_sp_plot_inv_kpc, cols3):
    rho_yukawa = ((1.0/r_plot_kpc**2 + m_sp_k/r_plot_kpc)**2
                  * np.exp(-2.0 * m_sp_k * r_plot_kpc))
    rho_yukawa_norm = rho_yukawa / rho_yukawa[idx_5]
    ax.loglog(r_plot_kpc, rho_yukawa_norm, color=col, lw=1.8,
              label=rf'ρ_field m_sp={m_sp_k:.1f} kpc⁻¹')

ax.set_title('M3: ρ_field TGP vs NFW (kształt)', fontsize=11)
ax.set_xlabel('r [kpc]', fontsize=10)
ax.set_ylabel('ρ (normalizowane do r=5kpc)', fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(0.5, 33); ax.set_ylim(1e-4, 1e3)

# 4. Pole tachionowe (M4)
ax = axes[1, 0]
ax.errorbar(OBS_R_KPC, OBS_V_KMS, yerr=OBS_ERR,
            fmt='ko', ms=5, capsize=3, label='NGC 3198 (obs.)', zorder=10)
ax.plot(r_kpc, v_bar/km_s, 'k--', lw=1.5, label='Newton')
cols4 = cm.Purples(np.linspace(0.4, 0.9, len(m4_curves)))
for (m_t, v_t), col in zip(m4_curves.items(), cols4):
    c2 = chi2_ngc(v_t)
    ax.plot(r_kpc, v_t, color=col, lw=1.5,
            label=rf'Tachion $m_t$={m_t:.0f} kpc⁻¹ ($\chi^2$={c2:.1f})')
ax.set_title('M4: Tachionowe pole (m²<0)', fontsize=11)
ax.set_xlabel('r [kpc]', fontsize=10)
ax.set_ylabel('v [km/s]', fontsize=10)
ax.legend(fontsize=8, loc='lower right')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 33); ax.set_ylim(0, 220)

# 5. Podsumowanie chi^2 wszystkich mechanizmów
ax = axes[1, 1]
mechanisms = ['Newton\n(baryony)', 'M1 Yukawa\nλ=100kpc', 'M2 V₃\nγ=1e-40',
              'M4 tachion\nm_t=2', 'CDM\n(ISO)']
chi2_vals = [
    chi2_ngc(v_bar/km_s),
    chi2_ngc(m1_curves[100.0]),
    chi2_ngc(m2_curves[1e-40]),
    chi2_ngc(m4_curves.get(2.0, v_bar/km_s)),
    chi2_cdm
]
colors_bar = ['gray', 'tomato', 'steelblue', 'mediumpurple', 'limegreen']
bars = ax.bar(mechanisms, chi2_vals, color=colors_bar, alpha=0.8, edgecolor='k')
ax.axhline(2.0, color='gold', ls='--', lw=2, label='Granica akceptacji (χ²=2)')
ax.axhline(5.0, color='orange', ls=':', lw=1.5, alpha=0.7, label='χ²=5')
for bar, val in zip(bars, chi2_vals):
    ax.text(bar.get_x() + bar.get_width()/2., min(val, 48),
            f'{val:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
ax.set_yscale('log')
ax.set_ylabel('χ²/N', fontsize=11)
ax.set_title('Porównanie χ²/N mechanizmów DM', fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim(0.1, 200)

# 6. Diagram pająk — skale i paradoks γ
ax = axes[1, 2]
categories = ['Efimov\nm_sp<0.12', 'Galaktyka\nλ~15kpc', 'FDM\nm~1e-22eV',
              'Trójciałowe\nM2 fit', 'Droga B\nszczupły']
m_sp_vals = [0.12, 1.0/(15.0*kpc_m/l_Pl_m),
             1e-22 / (1.22e28),
             np.sqrt(gamma_needed),
             1.0]
m_sp_log  = np.log10(np.array(m_sp_vals))

y_pos = np.arange(len(categories))
bars2 = ax.barh(y_pos, m_sp_log - m_sp_log.min(),
                left=m_sp_log.min(),
                color=['#2196F3', '#4CAF50', '#F44336', '#FF9800', '#9C27B0'],
                alpha=0.8, edgecolor='k')
ax.set_yticks(y_pos)
ax.set_yticklabels(categories, fontsize=9)
ax.set_xlabel('log₁₀(m_sp) [l_Pl⁻¹]', fontsize=10)
ax.set_title('Paradoks γ: m_sp dla różnych wymagań', fontsize=11)
for bar, val in zip(bars2, m_sp_log):
    ax.text(val + 0.5, bar.get_y() + bar.get_height()/2.,
            f'{val:.0f}', va='center', fontsize=8)
ax.grid(True, alpha=0.3, axis='x')
spread = m_sp_log.max() - m_sp_log.min()
ax.text(0.5, -0.15, f'Rozpiętość: {spread:.0f} rzędów wielkości!',
        transform=ax.transAxes, ha='center', fontsize=10,
        color='red', fontweight='bold')

plt.tight_layout()
out_png = os.path.join(os.path.dirname(__file__), 'ex35_dm_mechanisms.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()
print(f"Wykres zapisany: {out_png}")

# ============================================================
# Końcowe podsumowanie
# ============================================================
print()
print("=" * 70)
print("KOŃCOWE PODSUMOWANIE ex35")
print("=" * 70)
print()
print("FALSYFIKACJE POTWIERDZONE:")
print("  ✗  M1: Yukawa G_eff(r) < G — potwierdza ex14")
print("  ✗  M2: V₃ akumulacja — kształt krzywej niewłaściwy, γ paradoks")
print("  ✗  M3: ρ_field kształt ≠ NFW w żadnym zakresie λ")
print("  ✗  M4: Tachion — niestabilność teorii, oscylujący profil")
print()
print("OTWARTA MOŻLIWOŚĆ:")
print("  ?  M5: FDM z V_mod = εg² + V_TGP (patrz ex36)")
print("         Wymaga ε ~ 10^{-101} i nowego parametru.")
print()
print("PARADOKS γ (kluczowy wynik):")
print(f"  Rozpiętość m_sp dla różnych wymagań: ~{spread:.0f} rzędów")
print("  Efimov: m_sp < 0.12 (silna siła)")
print(f"  Galaktyka: m_sp ~ {1.0/(15.0*kpc_m/l_Pl_m):.2e} (długi zasięg)")
print(f"  FDM:       m_sp ~ {1e-22/(1.22e28):.2e} (ultra-lekki bozon)")
print("  Nie ma jednej wartości γ = m_sp² spełniającej WSZYSTKIE wymagania.")
print()
print("WNIOSEK:")
print("  Minimalne TGP (N0 aksjomaty) nie może wyjaśnić ciemnej materii.")
print("  Potwierdzono ANALIZA_CIEMNA_MATERIA.md i wyniki ex14.")
print("  FDM (ex36) pozostaje jedyną nadzieją kosztem ε-rozszerzenia.")
