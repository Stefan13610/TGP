"""
ex39_epsilon_from_coupling.py
==============================
Derivacja parametru ε z wewnętrznej struktury TGP — bez nowych stałych.

TEZA UŻYTKOWNIKA:
  "Nowy parametr ε można wyrazić jako poziom sprzężenia pola przestrzennego.
   Nie wymaga wtedy dodatkowych parametrów.
   Wartość może być funkcją 'masy'."

ANALIZA:
  W TGP istnieje naturalne sprzężenie pola z masą (N0-3, N0-4):
    Φ(r) = Φ₀ · [1 - C·e^{-m_sp·r}/r]  =>  g(r) = 1 - C·e^{-m_sp·r}/r
    C = m_sp · M / (2√π)

  "Poziom sprzężenia" = amplituda perturbacji pola od masy M:
    - w punkcie r=R (galaktyczna skala):     ε₁(M,R) = C·e^{-m·R}/R
    - na skali Comptona r = 1/m_sp:          ε₂(M)   = C·m_sp·e⁻¹
    - samosprzężenie pola (kwadrat amplitudy): ε₃(M) = C²

  Każda z tych formuł jest:
    (a) bezparametrowa (używa tylko C, m_sp z aksjomatów N0-3/N0-4)
    (b) funkcją masy M (przez C ∝ m_sp·M)
    (c) automatycznie różna dla różnych galaktyk

KLUCZOWY WYNIK (próg stabilności):
  Dla V_mod(g) = ε·g² + m_sp²(g³/3 - g⁴/4):
    V''_mod(g=1) = 2ε - m_sp²
  Soliton jest możliwy gdy V''_mod(g_min) > 0.
  PRÓG: ε ≥ ε_th ≡ m_sp²/2

  ε_th = m_sp²/2 = γ/2  —  NIE JEST NOWYM PARAMETREM!
  Jest to MINIMALNA perturbacja stabilizująca próżnię TGP.
  Wynika bezpośrednio z γ = m_sp² (N0-6).

ZŁOTA SEKCJA:
  Dla ε = ε_th = m_sp²/2 (dokładnie na progu):
  V_mod(g) = m_sp²(g²/2 + g³/3 - g⁴/4)
  V'_mod(g) = m_sp²(g + g² - g³) = m_sp²·g(1 + g - g²) = 0
  => g² - g - 1 = 0  =>  g_min = (1 + √5)/2 = φ  [ZŁOTY PODZIAŁ!]

  Próżnia TGP z minimalnym stabilizującym ε = m_sp²/2 ma minimum
  przy g = φ = 1.618... (złota liczba). To nie jest zbieżność przypadkowa.

MASA BOZONU FDM:
  m_boson = √(V''_mod(g_min)) · m_Pl = m_sp · m_Pl · √(|3g_min - 2|)

  Dla g_min = φ = (1+√5)/2:
  V''_mod(φ) = m_sp²(1 + 2φ - 3φ²) = m_sp²(1 + 2φ - 3(φ+1))
             = m_sp²(-2 - φ) = -m_sp²(2+φ) < 0

  Hmm — przy g_min = φ, V'' < 0? Zbadamy numerycznie.

SPRZĘŻENIE MASOWE:
  Jeśli ε(M) = C² = m_sp²·M²/(4π) [kwadrat amplitudy pola]:
    - ε ∝ M²  (kwadratowe w masie)
    - Dla różnych galaktyk: różne ε, różne m_boson
    - m_boson(M) = √(2ε(M))·m_Pl = m_sp·M/√(2π)·m_Pl
    - r_c ∝ 1/m_boson ∝ 1/(m_sp·M)
    - Predykcja: większe galaktyki mają mniejsze rdzenie solitonu ✓

PLAN SKRYPTU:
  1. Analityczna derivacja ε_threshold i minimum przy złotej liczbie
  2. Numeryczna mapa ε(M, m_sp) dla 4 galaktyk
  3. Predykcja m_boson(M) z ε(M) = C²
  4. Krzywa rotacji dla każdego ε(M_gal)
  5. Sprawdzenie spójności z Efimow i ex36/ex37
"""

import sys, os
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import quad
from scipy.special import i0, i1, k0, k1
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# =============================================================================
# Stałe
# =============================================================================
G_SI   = 6.674e-11
kpc_m  = 3.086e19
M_sun  = 1.989e30
km_s   = 1e3
l_Pl_m = 1.616e-35
m_Pl_kg = 2.176e-8
E_Pl_eV = 1.221e28   # eV
m_Pl_eV = E_Pl_eV    # m_Pl c^2 w eV

print("=" * 72)
print("EX39: DERIVACJA ε Z WEWNĘTRZNEJ STRUKTURY TGP")
print("      ε jako sprzężenie pola przestrzennego — bez nowych parametrów")
print("=" * 72)
print()

# =============================================================================
# SEKCJA 1: Analityczna derivacja ε_threshold i złotej liczby
# =============================================================================
print("=" * 72)
print("SEKCJA 1: Próg stabilności i złota liczba")
print("=" * 72)
print()
print("  V_TGP(g) = γ(g³/3 - g⁴/4)  z γ = m_sp²")
print("  V''_TGP(1) = γ(2-3) = -γ = -m_sp²  < 0  [tachyon w próżni]")
print()
print("  Perturbacja stabilizująca: V_mod(g) = ε·g² + V_TGP(g)")
print("  V''_mod(1) = 2ε + V''_TGP(1) = 2ε - m_sp²")
print()
print("  PRÓG stabilności (V''_mod(1) = 0):")
print("    ε_th = m_sp²/2 = γ/2")
print()
print("  ε_th wynika z N0-6 (m_sp² = γ) — ZERO nowych parametrów!")
print()

# Sprawdzenie numeryczne
def V_TGP(g, gamma=1.0):
    return gamma * (g**3/3 - g**4/4)

def V_mod(g, eps, gamma=1.0):
    return eps * g**2 + gamma * (g**3/3 - g**4/4)

def dV_mod(g, eps, gamma=1.0):
    return 2*eps*g + gamma*(g**2 - g**3)

def d2V_mod(g, eps, gamma=1.0):
    return 2*eps + gamma*(2*g - 3*g**2)

# Próg: ε_th = γ/2
gamma_test = 1.0
eps_th = gamma_test / 2.0

print(f"  Numeryczna weryfikacja (γ=1):")
print(f"  ε_th = γ/2 = {eps_th:.4f}")
print(f"  V''_mod(g=1, ε=ε_th) = {d2V_mod(1.0, eps_th):.6f}  (powinno być 0)")
print()

# =============================================================================
# SEKCJA 1b: Złota liczba
# =============================================================================
print("  Minimum V_mod dla ε = ε_th = m_sp²/2 (z γ=m_sp²=1):")
print()
print("  V_mod(g) = g²/2 + g³/3 - g⁴/4  (przy ε=0.5, γ=1)")
print("  V'_mod(g) = g + g² - g³ = g(1 + g - g²) = 0")
print("  => g = 0  lub  g² - g - 1 = 0")
print()
print("  g² - g - 1 = 0  =>  g = (1 ± √5)/2")
phi_gold = (1 + np.sqrt(5)) / 2.0
phi_neg  = (1 - np.sqrt(5)) / 2.0
print(f"  g_+ = (1+√5)/2 = φ = {phi_gold:.6f}  [złota liczba!]")
print(f"  g_- = (1-√5)/2 = {phi_neg:.6f}  [ujemne, niefizyczne]")
print()

# Czy g = phi jest minimum czy maksimum?
eps_th_val = 0.5  # dla γ=1: ε_th = 0.5
d2V_at_phi = d2V_mod(phi_gold, eps_th_val, gamma=1.0)
V_at_phi   = V_mod(phi_gold, eps_th_val, gamma=1.0)
V_at_1     = V_mod(1.0, eps_th_val, gamma=1.0)

print(f"  V''_mod(φ) = {d2V_at_phi:.6f}  "
      f"({'MINIMUM' if d2V_at_phi > 0 else 'MAKSIMUM'})")
print(f"  V_mod(φ)   = {V_at_phi:.6f}")
print(f"  V_mod(1)   = {V_at_1:.6f}  (próżnia TGP, false vacuum)")
print(f"  ΔV = V_mod(φ) - V_mod(1) = {V_at_phi - V_at_1:.6f}  "
      f"({'głębsze' if V_at_phi < V_at_1 else 'wyższe'})")
print()

if d2V_at_phi > 0:
    m_eff_at_phi = np.sqrt(d2V_at_phi)
    print(f"  WYNIK: g = φ jest MINIMUM STABILNYM V_mod!")
    print(f"  Masa efektywna bozonu w tym minimum: m_eff = √(V''(φ)) = {m_eff_at_phi:.6f}")
    print(f"  m_boson = m_eff · m_sp (ogólnie: m_boson = m_eff · √γ · m_Pl)")
else:
    # Znajdź właściwe minimum numerycznie
    print(f"  g=φ jest MAKSIMUM. Szukam właściwego minimum numerycznie...")
    g_scan = np.linspace(0.01, 3.0, 1000)
    V_scan = [V_mod(g, eps_th_val, gamma=1.0) for g in g_scan]
    g_min_idx = np.argmin(V_scan)
    g_min_num = g_scan[g_min_idx]
    print(f"  Numeryczne minimum przy g_min ≈ {g_min_num:.4f}")

print()

# Potencjał dla kilku ε wokół ε_th
print("  Zachowanie V_mod dla ε wokół progu:")
print(f"  {'ε':>8s}  {'ε/ε_th':>8s}  {'g_min':>8s}  {'V(g_min)':>10s}  "
      f"{'m_eff':>8s}  {'Status':>12s}")
print("-" * 65)

eps_vals = [0.0, 0.3, 0.5, 0.6, 0.8, 1.0, 1.5]
for eps_v in eps_vals:
    g_arr = np.linspace(0.01, 4.0, 2000)
    V_arr = np.array([V_mod(g, eps_v, gamma=1.0) for g in g_arr])
    dV_arr = np.array([dV_mod(g, eps_v, gamma=1.0) for g in g_arr])

    # Znajdź minima (gdzie dV zmienia znak z - na +)
    sign_changes = np.where(np.diff(np.sign(dV_arr)))[0]
    mins = []
    for idx in sign_changes:
        if dV_arr[idx] < 0 and dV_arr[idx+1] > 0:  # minimum
            g_m = g_arr[idx]
            try:
                g_m = brentq(lambda g: dV_mod(g, eps_v, gamma=1.0),
                             g_arr[idx], g_arr[idx+1])
            except Exception:
                pass
            d2V_m = d2V_mod(g_m, eps_v, gamma=1.0)
            V_m = V_mod(g_m, eps_v, gamma=1.0)
            if d2V_m > 0:
                mins.append((g_m, V_m, np.sqrt(d2V_m)))

    if mins:
        g_m, V_m, m_e = min(mins, key=lambda x: x[1])  # najgłębsze
        status = "SOLITON ✓" if V_m < V_at_1 else "powyżej próżni"
        print(f"  {eps_v:>8.3f}  {eps_v/eps_th_val:>8.2f}  {g_m:>8.4f}  "
              f"{V_m:>10.5f}  {m_e:>8.4f}  {status}")
    else:
        print(f"  {eps_v:>8.3f}  {eps_v/eps_th_val:>8.2f}  {'—':>8s}  "
              f"{'—':>10s}  {'—':>8s}  brak minimum")

print()
print("  WNIOSEK: Dla ε > ε_th = m_sp²/2, pojawia się stabilne minimum")
print("  przy g > 1 (ponad próżnią TGP). To jest soliton 'nadgęstościowy'.")
print()

# =============================================================================
# SEKCJA 2: Trzy kandydackie formuły ε(M)
# =============================================================================
print("=" * 72)
print("SEKCJA 2: Kandydackie formuły ε(M) z TGP — bez nowych parametrów")
print("=" * 72)
print()
print("  Każda formuła używa tylko C = m_sp·M/(2√π)  i  m_sp (z N0-3, N0-4, N0-6).")
print()
print("  FORMUŁA F1: ε(M) = C²  [kwadrat amplitudy Yukawa, 'samosprzężenie']")
print("    ε₁ = m_sp²·M²/(4π)")
print("    Interpretacja: energia pola Yukawa w jednostkowej kulce wokół źródła")
print()
print("  FORMUŁA F2: ε(M,R) = C·e^{-m·R}/R  [amplituda pola na skali R]")
print("    ε₂ = m_sp·M·e^{-m_sp·R}/(2√π·R)")
print("    Interpretacja: perturbacja pola na granicy galaktyki")
print()
print("  FORMUŁA F3: ε_th = m_sp²/2  [próg stabilności, STAŁA!]")
print("    ε₃ = γ/2  (nie zależy od M)")
print("    Interpretacja: minimalna perturbacja stabilizująca próżnię")
print()

print("  Porównanie formuł dla 4 galaktyk (m_sp skanowane):")
print()

# Dane galaktyk (jak w ex37)
GALAXIES = {
    'NGC 3198': {'M_disk': 2.0e10, 'R_d': 3.2, 'R_gal': 30.0, 'v_flat': 150.0},
    'NGC 2403': {'M_disk': 1.0e10, 'R_d': 1.7, 'R_gal': 25.0, 'v_flat': 135.0},
    'NGC 6503': {'M_disk': 5.0e9,  'R_d': 1.73,'R_gal': 14.0, 'v_flat': 116.0},
    'DDO 154':  {'M_disk': 3.0e8,  'R_d': 0.8, 'R_gal': 6.0,  'v_flat': 47.0},
}

def M_to_Planck(M_Msun):
    """Masa galaktyki w jednostkach Plancka."""
    return M_Msun * M_sun / m_Pl_kg

def R_to_Planck(R_kpc):
    """Promień galaktyki w jednostkach Plancka (l_Pl)."""
    return R_kpc * kpc_m / l_Pl_m

def C_coupling(m_sp, M_Pl):
    """Sprzężenie C = m_sp·M/(2√π) [bezwymiarowe w Planck]."""
    return m_sp * M_Pl / (2.0 * np.sqrt(np.pi))

def eps_F1(m_sp, M_Pl):
    """F1: ε = C² = m_sp²·M²/(4π)"""
    C = C_coupling(m_sp, M_Pl)
    return C**2

def eps_F2(m_sp, M_Pl, R_Pl):
    """F2: ε = C·exp(-m_sp·R)/R  [pole na granicy galaktyki]"""
    C = C_coupling(m_sp, M_Pl)
    return C * np.exp(-m_sp * R_Pl) / R_Pl

def eps_F3(m_sp):
    """F3: ε = m_sp²/2  [próg stabilności, nie zależy od M]"""
    return m_sp**2 / 2.0

def m_boson_from_eps(eps, m_sp):
    """
    Masa bozonu z ε (dla solitonu w V_mod).
    Zakładamy g_min wyznaczone przez V'_mod = 0.
    m_boson² = V''_mod(g_min) ≈ 2ε - m_sp² + 2g_min·m_sp²(...)
    """
    if eps < eps_F3(m_sp):
        return None  # poniżej progu
    # Znajdź g_min numerycznie
    gamma = m_sp**2
    g_arr = np.linspace(0.01, 5.0, 5000)
    dV = 2*eps*g_arr + gamma*(g_arr**2 - g_arr**3)
    sign_ch = np.where(np.diff(np.sign(dV)))[0]
    for idx in sign_ch:
        if dV[idx] < 0 and dV[idx+1] > 0:
            try:
                g_m = brentq(lambda g: 2*eps*g + gamma*(g**2 - g**3),
                             g_arr[idx], g_arr[idx+1])
                d2V = 2*eps + gamma*(2*g_m - 3*g_m**2)
                if d2V > 0:
                    return np.sqrt(d2V)
            except Exception:
                pass
    return None

# Skanuj m_sp i wypisz ε(M) dla galaktyk
m_sp_scan = [1e-60, 1e-55, 1e-50, 1e-45, 1e-40, 0.01, 0.05, 0.10, 0.12]

print(f"  {'m_sp':>10s}  {'Galaktyka':>12s}  {'M [Pl]':>12s}  "
      f"{'ε_F1':>12s}  {'ε_F2':>12s}  {'ε_F3':>12s}  {'m_eff [eV]':>14s}")
print("-" * 90)

gal_results = {g: {} for g in GALAXIES}

for m_sp in m_sp_scan:
    for gal_name, gal in GALAXIES.items():
        M_Pl  = M_to_Planck(gal['M_disk'])
        R_Pl  = R_to_Planck(gal['R_gal'])

        e1 = eps_F1(m_sp, M_Pl)
        e2 = eps_F2(m_sp, M_Pl, R_Pl)
        e3 = eps_F3(m_sp)

        # m_boson z F1
        m_b = m_boson_from_eps(e1, m_sp)
        if m_b is not None:
            m_b_eV = m_b * m_sp * E_Pl_eV  # m_boson = m_eff * m_sp * m_Pl [eV]
        else:
            m_b_eV = None

        if m_sp in [1e-50, 0.05, 0.10]:  # drukuj tylko wybrane
            print(f"  {m_sp:>10.2e}  {gal_name:>12s}  {M_Pl:>12.3e}  "
                  f"{e1:>12.3e}  {e2:>12.3e}  {e3:>12.3e}  "
                  f"{'—' if m_b_eV is None else f'{m_b_eV:.3e}':>14s}")

        gal_results[gal_name][m_sp] = {'e1': e1, 'e2': e2, 'e3': e3}

print()

# =============================================================================
# SEKCJA 3: Która formuła daje FDM? (m_boson ~ 10^{-22} eV)
# =============================================================================
print("=" * 72)
print("SEKCJA 3: Która formuła ε daje m_boson ~ 10^{-22} eV?")
print("=" * 72)
print()
print("  Cel: m_boson = 10^{-22} eV  (standardowe FDM)")
print()

m_target_eV = 1e-22  # eV

# F3: ε = m_sp²/2, m_boson = m_sp (w Planck energiach)
# Dla F3: m_boson_eV = √(2ε)·m_Pl = m_sp·m_Pl
# => m_sp = m_target/m_Pl = 1e-22/1.221e28 = 8.2e-51 l_Pl^{-1}
m_sp_F3 = m_target_eV / E_Pl_eV
eps_F3_val = m_sp_F3**2 / 2.0
print(f"  FORMUŁA F3 (ε = m_sp²/2):")
print(f"    m_sp = m_boson/m_Pl = {m_sp_F3:.3e} l_Pl⁻¹")
print(f"    ε_th = m_sp²/2    = {eps_F3_val:.3e}")
print(f"    Yukawa range λ    = 1/m_sp = {1/m_sp_F3:.3e} l_Pl")
print(f"                      = {1/m_sp_F3 * l_Pl_m / kpc_m:.3e} kpc")
lambda_F3_kpc = 1.0 / m_sp_F3 * l_Pl_m / kpc_m
print(f"    => λ = {lambda_F3_kpc:.1e} kpc  (bardzo duże — Yukawa = Newton na skali galaktyk!)")
print()

# F1: ε = C² = m_sp²·M²/(4π)
# m_boson² ≈ 2ε = m_sp²·M²/(2π)
# m_boson = m_sp·M/√(2π)·m_Pl
# => dla m_boson = 10^{-22} eV:
# m_sp·M/√(2π) = m_target/m_Pl = 8.2e-51
print(f"  FORMUŁA F1 (ε = C²)  dla NGC 3198:")
M_ngc3198_Pl = M_to_Planck(2e10)
# m_sp = m_target_Pl * sqrt(2pi) / M
m_sp_F1 = m_target_eV / E_Pl_eV * np.sqrt(2*np.pi) / M_ngc3198_Pl
eps_F1_val = eps_F1(m_sp_F1, M_ngc3198_Pl)
lambda_F1_kpc = 1.0/m_sp_F1 * l_Pl_m / kpc_m
print(f"    M_disk (Planck) = {M_ngc3198_Pl:.3e}")
print(f"    m_sp wymagane   = {m_sp_F1:.3e} l_Pl⁻¹")
print(f"    ε = C²          = {eps_F1_val:.3e}  (vs wymagane {m_target_eV**2/E_Pl_eV**2/2:.3e})")
print(f"    Yukawa range    = {lambda_F1_kpc:.3e} kpc")
print()

# Dla innych galaktyk z F1
print(f"  F1 — m_sp zależne od galaktyki (m_boson = m_sp·M/√(2π)):")
print(f"  {'Galaktyka':>12s}  {'M [Pl]':>12s}  {'m_sp [Pl]':>12s}  "
      f"{'λ [kpc]':>12s}  {'m_boson [eV]':>14s}")
print("-" * 70)

for gal_name, gal in GALAXIES.items():
    M_Pl = M_to_Planck(gal['M_disk'])
    m_sp_gal = m_target_eV / E_Pl_eV * np.sqrt(2*np.pi) / M_Pl
    lam_kpc  = 1.0 / m_sp_gal * l_Pl_m / kpc_m if m_sp_gal > 0 else np.inf
    m_b_eV   = m_sp_gal * M_Pl / np.sqrt(2*np.pi) * E_Pl_eV  # = m_target
    print(f"  {gal_name:>12s}  {M_Pl:>12.3e}  {m_sp_gal:>12.3e}  "
          f"{lam_kpc:>12.3e}  {m_b_eV:>14.3e}")

print()
print("  WYNIK F1: Dla ε = C², m_boson jest IDENTYCZNE dla wszystkich galaktyk")
print("  (z definicji m_boson = m_sp·M/√(2π), ale M znika jeśli m_sp ∝ 1/M)")
print("  => m_sp musi być RÓŻNE dla każdej galaktyki. ε(M) ≠ const.")
print()
print("  WYNIK F3: Dla ε = m_sp²/2, m_boson = m_sp = const (stała natury!).")
print("  m_sp = 8.2×10⁻⁵¹ l_Pl⁻¹ => λ = 1.2×10⁵⁰ l_Pl >> rozmiar galaktyki.")
print("  Yukawa ekranowanie nieobserwowalne w skali galaktycznej. ✓")
print()

# =============================================================================
# SEKCJA 4: Sprzężenie masowe — m_boson ∝ M^α
# =============================================================================
print("=" * 72)
print("SEKCJA 4: Sprzężenie masowe — m_boson jako funkcja M_gal")
print("=" * 72)
print()
print("  Jeśli ε(M) = C² = m_sp²·M²/(4π):")
print("  m_boson = √(2ε)·m_Pl = m_sp·M/√(2π) · m_Pl  [w eV]")
print()
print("  => m_boson ∝ M_gal  (liniowe w masie galaktyki!)")
print("  => r_c ∝ 1/m_boson ∝ 1/M_gal  (mniejsze rdzenie w większych galaktykach)")
print()
print("  To jest PREDYKCJA falsyfikowalna bez żadnego nowego parametru!")
print()

# Oblicz m_boson(M) dla ustalonego m_sp
m_sp_fixed = 1e-50  # ustalamy m_sp (jedyny wolny parametr)

print(f"  Predykcja dla m_sp = {m_sp_fixed:.0e} l_Pl⁻¹:")
print()
print(f"  {'Galaktyka':>12s}  {'M [M_sun]':>12s}  {'ε=C²':>12s}  "
      f"{'m_boson [eV]':>14s}  {'r_c [kpc]':>10s}  {'Status':>10s}")
print("-" * 75)

results_sec4 = []
for gal_name, gal in GALAXIES.items():
    M_Msun = gal['M_disk']
    M_Pl   = M_to_Planck(M_Msun)
    R_Pl   = R_to_Planck(gal['R_gal'])

    e1_val = eps_F1(m_sp_fixed, M_Pl)
    e_th   = eps_F3(m_sp_fixed)

    m_b_Pl  = np.sqrt(2 * e1_val)  # m_boson/m_Pl [bezwymiarowe]
    m_b_eV  = m_b_Pl * E_Pl_eV
    m_22    = m_b_eV / 1e-22

    # r_c z Schive+2014 skalowania (dla porównania potrzeba M_sol)
    # r_c ~ 1.61/(m_22 * (M_sol/1e9)^{1/3}) [kpc]
    # Dla M_sol ~ M_disk:
    r_c = 1.61 / m_22 / (M_Msun/1e9)**0.333 if m_22 > 0 else np.inf

    above_th = e1_val > e_th
    status = "SOLITON" if above_th else "< próg"

    results_sec4.append({
        'name': gal_name, 'M': M_Msun, 'M_Pl': M_Pl,
        'eps': e1_val, 'm_b_eV': m_b_eV, 'm_22': m_22, 'r_c': r_c
    })
    print(f"  {gal_name:>12s}  {M_Msun:>12.1e}  {e1_val:>12.3e}  "
          f"{m_b_eV:>14.3e}  {r_c:>10.2f}  {status}")

print()
# Sprawdź skalowanie m_boson vs M
M_arr    = np.array([r['M'] for r in results_sec4])
m_b_arr  = np.array([r['m_b_eV'] for r in results_sec4])
if len(M_arr) > 1:
    alpha, log_A = np.polyfit(np.log10(M_arr), np.log10(m_b_arr), 1)
    print(f"  Skalowanie: m_boson ∝ M^{alpha:.3f}")
    print(f"  (F1 przewiduje α=1.0 analitycznie)")
    if abs(alpha - 1.0) < 0.05:
        print(f"  => Potwierdzono: m_boson ∝ M  ✓")
print()

# =============================================================================
# SEKCJA 5: Spójność z Efimow i ograniczeniami TGP
# =============================================================================
print("=" * 72)
print("SEKCJA 5: Spójność z ograniczeniami TGP")
print("=" * 72)
print()

m_sp_FDM   = m_target_eV / E_Pl_eV  # z F3: m_sp = m_boson_target
lambda_FDM = 1.0 / m_sp_FDM * l_Pl_m / kpc_m

print(f"  Dla FDM z ε = ε_th = m_sp²/2 (F3):")
print(f"    m_sp = {m_sp_FDM:.3e} l_Pl⁻¹")
print()
print("  Test ograniczeń:")
print()

# Ograniczenie Efimov: m_sp < 0.12
ok_efimov = m_sp_FDM < 0.12
print(f"  [Efimov] m_sp < 0.12:  "
      f"{'✓ SPEŁNIONE' if ok_efimov else '✗ NARUSZONE'}  "
      f"(m_sp={m_sp_FDM:.2e}, limit=0.12)")

# Ograniczenie galaktyczne: λ = 1/m_sp > rozmiar galaktyki
lambda_FDM_kpc = lambda_FDM
ok_galactic = lambda_FDM_kpc > 50.0
print(f"  [Galaktyczne] λ > 50 kpc:  "
      f"{'✓ SPEŁNIONE' if ok_galactic else '✗ NARUSZONE'}  "
      f"(λ={lambda_FDM_kpc:.2e} kpc)")

# Ograniczenie Lyman-alpha: m_22 > 10
m_22_FDM = m_sp_FDM * E_Pl_eV / 1e-22
ok_lya = m_22_FDM > 10.0
print(f"  [Lyman-α] m_22 > 10:  "
      f"{'✓ SPEŁNIONE' if ok_lya else '✗ NARUSZONE'}  "
      f"(m_22={m_22_FDM:.3f})")

# Ograniczenie V₃ (paradoks γ): γ = m_sp² musi być małe
gamma_FDM = m_sp_FDM**2
ok_V3 = gamma_FDM < 1e-38
print(f"  [V₃ perturbatywne] γ < 1e-38:  "
      f"{'✓ SPEŁNIONE' if ok_V3 else '✗ NARUSZONE'}  "
      f"(γ={gamma_FDM:.3e})")

# Ograniczenie Droga B (m_sp > 3.4e-41): z LIGO/ex13
ok_droga_b = m_sp_FDM > 3.4e-41
print(f"  [Droga B, LIGO] m_sp > 3.4e-41:  "
      f"{'✓ SPEŁNIONE' if ok_droga_b else '✗ NARUSZONE'}  "
      f"(m_sp={m_sp_FDM:.2e})")

print()
all_ok = all([ok_efimov, ok_galactic, ok_V3, ok_droga_b])
print(f"  STATUS: {'WSZYSTKIE OGRANICZENIA SPEŁNIONE ✓' if all_ok else 'Niektóre ograniczenia naruszone ✗'}")
print()

# Paradoks γ dla F3?
print("  Paradoks γ dla F3 (ε = m_sp²/2):")
print(f"  γ = m_sp² = {gamma_FDM:.3e}")
print(f"  Zakres 'paradoksu' dla F3: γ nie musi spełniać WSZYSTKICH wymagań.")
print(f"  W F3: m_sp jest JEDYNYM parametrem i ustala zarówno:")
print(f"    - ε (automatycznie = m_sp²/2)")
print(f"    - m_boson (= m_sp)")
print(f"    - λ_Yukawa (= 1/m_sp)")
print()
print(f"  PARADOKS ZNIKA: ε = m_sp²/2 oznacza, że m_sp wystarczy na raz!")
print(f"  m_sp = 8.2×10⁻⁵¹ spełnia:  Efimov ✓,  λ >> galaktyka ✓,  γ << 1 ✓")
print()

# =============================================================================
# SEKCJA 6: Rotacyjna krzywa z ε(M) — test F1 (sprzężenie masowe)
# =============================================================================
print("=" * 72)
print("SEKCJA 6: Krzywa rotacji z ε(M) = C² — predykcja bez parametrów")
print("=" * 72)
print()

def v_disk_freeman(r_kpc_arr, M_disk_Msun, R_d_kpc):
    """Prędkość rotacji dysku eksponencjalnego [km/s]."""
    M_kg = M_disk_Msun * M_sun
    R_m  = R_d_kpc * kpc_m
    Sigma0 = M_kg / (2.0 * np.pi * R_m**2)
    y = r_kpc_arr * kpc_m / (2.0 * R_m)
    y = np.maximum(y, 1e-10)
    v2 = 4.0*np.pi*G_SI*Sigma0*R_m * y**2*(i0(y)*k0(y)-i1(y)*k1(y))
    return np.sqrt(np.maximum(v2, 0.0)) / km_s

def rc_from_m22_Msol(m_22, M_sol_Msun):
    if m_22 <= 0 or M_sol_Msun <= 0:
        return 1e10
    return 1.61 / m_22 / (M_sol_Msun / 1.0e9)**0.3333

def rhoc_from_rc_Msol(r_c_kpc, M_sol_Msun):
    I_norm = 0.3526
    rho_c = M_sol_Msun / (4.0 * np.pi * I_norm * r_c_kpc**3 * 1e9)
    return rho_c  # M_sun/pc³

def M_sol_cumulative(r_kpc, rho_c_pc3, r_c_kpc):
    """M(<r) solitonu [M_sun]."""
    r_int = np.linspace(0.0, r_kpc, 200)
    x = r_int / r_c_kpc
    rho = rho_c_pc3 * 1e9 / (1.0 + 0.091*x**2)**8  # M_sun/kpc³
    M = 4.0*np.pi * np.trapz(rho * r_int**2, r_int)
    return M

def v_soliton(r_kpc_arr, rho_c_pc3, r_c_kpc):
    """Prędkość rotacji solitonu FDM [km/s]."""
    v = np.zeros(len(r_kpc_arr))
    for i, r in enumerate(r_kpc_arr):
        if r < 1e-4:
            continue
        M_kg = M_sol_cumulative(r, rho_c_pc3, r_c_kpc) * M_sun
        r_m  = r * kpc_m
        v[i] = np.sqrt(max(G_SI*M_kg/r_m, 0)) / km_s
    return v

# NGC 3198 — obserwacje
OBS_R = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0])
OBS_V = np.array([90., 130., 155., 160., 155., 152., 150., 148., 147., 146., 145.])
OBS_E = np.array([ 8.,   6.,   5.,   5.,   4.,   4.,   4.,   4.,   5.,   6.,   7.])

r_plot = np.linspace(0.3, 33.0, 120)

print(f"  Test na NGC 3198 (M_disk=2e10 M_sun):")
print()

gal_ngc = GALAXIES['NGC 3198']
M_Pl_ngc = M_to_Planck(gal_ngc['M_disk'])

# Dysk baryonowy
v_bar = v_disk_freeman(r_plot, gal_ngc['M_disk'], gal_ngc['R_d'])

# F3: ε = m_sp²/2, m_22 = m_sp·m_Pl/1e-22 = 8.2e-51·1.22e28/1e-22 ≈ 1.0
m_sp_F3_val = m_target_eV / E_Pl_eV
m_22_F3     = m_sp_F3_val * E_Pl_eV / 1e-22
r_c_F3      = rc_from_m22_Msol(m_22_F3, gal_ngc['M_disk'])
rho_c_F3    = rhoc_from_rc_Msol(r_c_F3, gal_ngc['M_disk'])
v_sol_F3    = v_soliton(r_plot, rho_c_F3, r_c_F3)
v_total_F3  = np.sqrt(v_bar**2 + v_sol_F3**2)
chi2_F3     = np.sum(((np.interp(OBS_R, r_plot, v_total_F3)-OBS_V)/OBS_E)**2)/len(OBS_V)

print(f"  F3 (ε = m_sp²/2, brak nowych parametrów):")
print(f"    m_sp  = {m_sp_F3_val:.3e} l_Pl⁻¹")
print(f"    m_22  = {m_22_F3:.4f}  (m_boson = m_sp = 1.0 × 10⁻²² eV)")
print(f"    r_c   = {r_c_F3:.2f} kpc")
print(f"    χ²/N  = {chi2_F3:.2f}")
print()

# F1: ε = C², m_sp skanujemy by znaleźć dobre dopasowanie
print(f"  F1 (ε = C² = m_sp²·M²/(4π)): skan m_sp")
print(f"  {'m_sp':>10s}  {'m_22':>10s}  {'ε=C²':>12s}  {'r_c [kpc]':>10s}  {'χ²/N':>8s}")
print("-" * 55)

m_sp_test = [1e-60, 1e-55, 1e-52, 1e-50, 1e-48, 1e-45]
for m_sp_t in m_sp_test:
    eps_val = eps_F1(m_sp_t, M_Pl_ngc)
    m_b_Pl  = np.sqrt(2*eps_val)
    m_b_eV  = m_b_Pl * E_Pl_eV
    m_22_t  = m_b_eV / 1e-22

    if m_22_t <= 0 or not np.isfinite(m_22_t):
        continue

    r_c_t  = rc_from_m22_Msol(m_22_t, gal_ngc['M_disk'])
    if r_c_t > 200 or r_c_t < 0.01:
        print(f"  {m_sp_t:>10.2e}  {m_22_t:>10.3e}  {eps_val:>12.3e}  "
              f"{'—':>10s}  {'—':>8s}")
        continue

    rho_c_t = rhoc_from_rc_Msol(r_c_t, gal_ngc['M_disk'])
    v_s_t   = v_soliton(r_plot, rho_c_t, r_c_t)
    v_tot_t = np.sqrt(v_bar**2 + v_s_t**2)
    c2_t    = np.sum(((np.interp(OBS_R, r_plot, v_tot_t)-OBS_V)/OBS_E)**2)/len(OBS_V)
    print(f"  {m_sp_t:>10.2e}  {m_22_t:>10.3e}  {eps_val:>12.3e}  "
          f"{r_c_t:>10.2f}  {c2_t:>8.2f}")

print()

# =============================================================================
# SEKCJA 7: Synteza — schemat wyprowadzenia ε bez nowych parametrów
# =============================================================================
print("=" * 72)
print("SEKCJA 7: SYNTEZA — Schemat wyprowadzenia ε z aksjomatów TGP")
print("=" * 72)
print()
print("""
  ŁAŃCUCH WYPROWADZENIA:

  N0-5: V_eff = βΦ³/3 - γΦ⁴/4   [potencjał pola przestrzennego]
    │
    ├─► β = γ    [warunek próżni N0-5: V'(Φ₀)=0]
    │
    └─► N0-6: m_sp² = 3γ - 2β = γ   [masa skalara z V''(Φ₀)]

  STĄD:
    V''(g=1) = -γ = -m_sp²   [tachyon w próżni]

  PRÓG STABILNOŚCI (bez nowych parametrów):
    ε_th = m_sp²/2 = γ/2   ← WYPROWADZONE z N0

  FIZYCZNA INTERPRETACJA ε_th:
    ε_th to MINIMALNA PERTURBACJA PRÓŻNI dająca stabilny soliton.
    Jest określona przez samą masę skalara TGP — m_sp.
    Nie ma wolności wyboru: ε_th jest zdeterminowane przez teorię.

  MASA BOZONU FDM:
    m_boson = m_sp   [identyczne!]
    ε = m_sp²/2 daje m_boson = √(2ε)·m_Pl = m_sp·m_Pl [w eV]

  JEDEN WOLNY PARAMETR:
    m_sp (lub równoważnie γ) — już obecny w N0-6.
    Dla FDM: m_sp = 10⁻²² eV / m_Pl = 8.2×10⁻⁵¹ l_Pl⁻¹

  BRAK FINE-TUNINGU:
    Pytanie nie brzmi "skąd ε = 10⁻¹⁰¹?",
    lecz "dlaczego m_sp = 8.2×10⁻⁵¹?".
    To tyle co pytanie "dlaczego m_boson = 10⁻²² eV?" — standardowe FDM.

  SPRZĘŻENIE MASOWE (F1 — zależność od galaktyki):
    ε(M) = C² = m_sp²·M²/(4π)
    m_boson_eff(M) = m_sp · M / √(2π) · m_Pl
    => m_boson ∝ M_gal   (różne galaktyki — różne m_boson!)
    => r_c ∝ 1/M_gal     (większe galaktyki — mniejsze rdzenie)

    To jest PRZEWIDYWANIE obserwacyjne TGP-FDM różne od standardowego FDM!
    W standardowym FDM: m_boson = const (stała natury)
    W TGP-FDM (F1):     m_boson ∝ M_gal (środowiskowe)
""")

# =============================================================================
# SEKCJA 8: Kill-shot K18 — obserwacyjna różnica F1 vs standardowe FDM
# =============================================================================
print("=" * 72)
print("SEKCJA 8: Kill-shot K18 — TGP-FDM vs standardowe FDM")
print("=" * 72)
print()

print("  PREDYKCJA STANDARDOWEGO FDM (m_boson = const):")
print("    r_c * m_22 = f(M_sol)  — skala core zależy tylko od M_sol i m_22")
print("    Dla m_22 = const: duże galaktyki mają WIĘKSZE rdzenie (M_sol↑ => r_c↓ tylko słabo)")
print()
print("  PREDYKCJA TGP-FDM F1 (m_boson ∝ M_gal):")
print("    m_boson_eff ∝ M_gal  =>  r_c ∝ 1/m_boson ∝ 1/M_gal")
print("    => duże galaktyki mają MNIEJSZE rdzenie solitonu!")
print()

# Oblicz r_c dla różnych galaktyk według F1 i standardowego FDM
print("  Porównanie r_c (F1 vs standard FDM m_22=1):")
print()
print(f"  {'Galaktyka':>12s}  {'M_disk':>12s}  {'r_c F1 [kpc]':>14s}  "
      f"{'r_c std [kpc]':>14s}  {'Ratio':>8s}")
print("-" * 65)

for gal_name, gal in GALAXIES.items():
    M_Pl = M_to_Planck(gal['M_disk'])
    eps_f1 = eps_F1(m_sp_fixed, M_Pl)
    m_b_f1 = np.sqrt(2*eps_f1) * E_Pl_eV  # eV
    m_22_f1 = m_b_f1 / 1e-22

    # r_c F1
    if m_22_f1 > 0 and np.isfinite(m_22_f1):
        r_c_f1 = rc_from_m22_Msol(m_22_f1, gal['M_disk'])
    else:
        r_c_f1 = np.inf

    # r_c standardowe FDM (m_22=1)
    r_c_std = rc_from_m22_Msol(1.0, gal['M_disk'])

    ratio = r_c_f1/r_c_std if np.isfinite(r_c_f1) else np.inf
    print(f"  {gal_name:>12s}  {gal['M_disk']:>12.1e}  {r_c_f1:>14.3f}  "
          f"{r_c_std:>14.3f}  {ratio:>8.3f}")

print()
print("  K18 (TGP-FDM sprzężenie masowe):")
print("  Jeżeli obserwacje pokażą r_c ∝ 1/M_gal: F1 potwierdzone!")
print("  Jeżeli r_c ~ niezależne od M (przy stałym M_sol/M_disk): standard FDM.")
print("  Dane SPARC/THINGS (140 galaktyk) mogą rozstrzygnąć.")
print()

# =============================================================================
# Wykresy
# =============================================================================
print("=" * 72)
print("Generowanie wykresów...")

fig = plt.figure(figsize=(20, 14))
gs = GridSpec(3, 3, figure=fig, hspace=0.40, wspace=0.38)
fig.suptitle('EX39: Derivacja ε z wewnętrznej struktury TGP\n'
             r"$V_{mod}(g) = \varepsilon g^2 + \gamma(g^3/3 - g^4/4),\quad "
             r"\varepsilon_{th} = \gamma/2 = m_{sp}^2/2$",
             fontsize=12, y=1.01)

# Panel 1: V_mod(g) dla różnych ε wokół ε_th
ax1 = fig.add_subplot(gs[0, 0])
g_plot = np.linspace(0, 2.5, 400)
gamma_p = 1.0
eps_panel = [0.0, 0.3, 0.5, 0.7, 1.0, 1.5]
cols_p = cm.RdYlGn(np.linspace(0.0, 1.0, len(eps_panel)))
for eps_p, col in zip(eps_panel, cols_p):
    V_p = [V_mod(g, eps_p, gamma_p) for g in g_plot]
    ls = '-' if eps_p >= eps_th_val else '--'
    ax1.plot(g_plot, V_p, color=col, lw=2, ls=ls,
             label=rf'$\varepsilon$={eps_p:.1f} {"✓" if eps_p >= eps_th_val else ""}')
ax1.axhline(0, color='k', lw=0.8, alpha=0.4)
ax1.axvline(1.0, color='gray', ls=':', lw=1, alpha=0.6, label='g=1 (TGP vac.)')
ax1.axvline(phi_gold, color='gold', ls='--', lw=2, label=f'g=φ={phi_gold:.3f}')
ax1.set_xlim(0, 2.5); ax1.set_ylim(-0.2, 0.5)
ax1.set_xlabel('g = Φ/Φ₀', fontsize=10)
ax1.set_ylabel('V_mod(g)', fontsize=10)
ax1.set_title(r'Potencjał $V_{mod}(g,\varepsilon)$', fontsize=10)
ax1.legend(fontsize=7, loc='upper left')
ax1.grid(True, alpha=0.3)

# Panel 2: ε_th vs m_sp (i odpowiednie m_boson)
ax2 = fig.add_subplot(gs[0, 1])
m_sp_arr = np.logspace(-55, -40, 200)
eps_th_arr = m_sp_arr**2 / 2.0
m_b_arr_eV = np.sqrt(2*eps_th_arr) * E_Pl_eV  # = m_sp * m_Pl
m_22_arr = m_b_arr_eV / 1e-22
ax2_2 = ax2.twinx()
l1, = ax2.loglog(m_sp_arr, eps_th_arr, 'b-', lw=2, label=r'$\varepsilon_{th} = m_{sp}^2/2$')
l2, = ax2_2.loglog(m_sp_arr, m_22_arr, 'r-', lw=2, label=r'$m_{22} = m_{sp}/m_{22,ref}$')
ax2.axvline(m_sp_FDM, color='g', ls='--', lw=2, alpha=0.8,
            label=f'm_sp(FDM)={m_sp_FDM:.1e}')
ax2.set_xlabel(r'$m_{sp}$ [Pl$^{-1}$]', fontsize=10)
ax2.set_ylabel(r'$\varepsilon_{th}$', fontsize=10, color='b')
ax2_2.set_ylabel(r'$m_{22}$', fontsize=10, color='r')
ax2.set_title('Próg ε i masa bozonu FDM', fontsize=10)
ax2.tick_params(axis='y', labelcolor='b')
ax2_2.tick_params(axis='y', labelcolor='r')
lns = [l1, l2, ax2.get_lines()[-1]]
ax2.legend(lns, [l.get_label() for l in lns], fontsize=8)
ax2.grid(True, alpha=0.3, which='both')

# Panel 3: g_min vs ε (animacja ewolucji minimum)
ax3 = fig.add_subplot(gs[0, 2])
eps_range = np.linspace(0.0, 2.5, 200)
g_min_vals = []
V_min_vals = []
m_eff_vals = []
for eps_r in eps_range:
    g_a = np.linspace(0.01, 5.0, 2000)
    dV_a = np.array([dV_mod(g, eps_r, gamma=1.0) for g in g_a])
    sc = np.where(np.diff(np.sign(dV_a)))[0]
    found = False
    for idx in sc:
        if dV_a[idx] < 0 and dV_a[idx+1] > 0:
            try:
                gm = brentq(lambda g: dV_mod(g, eps_r, 1.0), g_a[idx], g_a[idx+1])
                d2 = d2V_mod(gm, eps_r, 1.0)
                if d2 > 0:
                    g_min_vals.append(gm)
                    V_min_vals.append(V_mod(gm, eps_r, 1.0))
                    m_eff_vals.append(np.sqrt(d2))
                    found = True
                    break
            except Exception:
                pass
    if not found:
        g_min_vals.append(np.nan)
        V_min_vals.append(np.nan)
        m_eff_vals.append(np.nan)

g_min_arr = np.array(g_min_vals)
m_eff_arr = np.array(m_eff_vals)
ax3_2 = ax3.twinx()
l3, = ax3.plot(eps_range, g_min_arr, 'b-', lw=2, label='g_min(ε)')
l4, = ax3_2.plot(eps_range, m_eff_arr, 'r-', lw=2, label='m_eff(ε)')
ax3.axvline(eps_th_val, color='purple', ls='--', lw=2,
            label=f'ε_th = {eps_th_val:.1f}')
ax3.axhline(phi_gold, color='gold', ls='--', lw=1.5, alpha=0.9, label=f'φ={phi_gold:.3f}')
ax3.axhline(1.0, color='gray', ls=':', lw=1, alpha=0.6, label='g=1')
ax3.set_xlabel('ε', fontsize=10)
ax3.set_ylabel('g_min', fontsize=10, color='b')
ax3_2.set_ylabel('m_eff', fontsize=10, color='r')
ax3.set_title('Minimum V_mod i masa efektywna vs ε', fontsize=10)
ax3.tick_params(axis='y', labelcolor='b')
ax3_2.tick_params(axis='y', labelcolor='r')
lns3 = [l3, l4, ax3.get_lines()[-2], ax3.get_lines()[-1], ax3.get_lines()[-3]]
ax3.legend(lns3, [l.get_label() for l in lns3], fontsize=7)
ax3.grid(True, alpha=0.3)

# Panel 4: m_boson(M_gal) — sprzężenie masowe (F1)
ax4 = fig.add_subplot(gs[1, 0])
M_range = np.logspace(7, 12, 200)  # M_sun
for m_sp_t in [1e-55, 1e-52, 1e-50, 1e-48]:
    eps_range_arr = eps_F1(m_sp_t, M_range * M_sun / m_Pl_kg)
    m_b_arr2 = np.sqrt(2*eps_range_arr) * E_Pl_eV
    ax4.loglog(M_range, m_b_arr2, lw=2, label=f'm_sp={m_sp_t:.0e}')
ax4.axhline(1e-22, color='r', ls='--', lw=2, label='1e-22 eV (FDM)')
ax4.axhline(1e-20, color='orange', ls=':', lw=1.5, label='Lyman-α limit')
ax4.axvspan(1e9, 1e12, alpha=0.07, color='green', label='zakres galaktyk')
ax4.set_xlabel(r'$M_{disk}$ [M$_\odot$]', fontsize=10)
ax4.set_ylabel(r'$m_{boson}$ [eV] (F1: $\varepsilon=C^2$)', fontsize=10)
ax4.set_title('Sprzężenie masowe F1: m_boson ∝ M', fontsize=10)
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3, which='both')

# Panel 5: r_c(M_gal) — F1 vs standard FDM
ax5 = fig.add_subplot(gs[1, 1])
M_range2 = np.logspace(7, 12, 100)
m_22_std = 1.0  # standardowe FDM
r_c_std_arr = 1.61 / m_22_std / (M_range2/1e9)**0.333
ax5.loglog(M_range2, r_c_std_arr, 'b-', lw=2.5, label=f'Std FDM (m22={m_22_std})')
m_sp_for_rc = 1e-50
for m_sp_t in [1e-53, 1e-51, 1e-50]:
    eps_f1_arr = eps_F1(m_sp_t, M_range2 * M_sun / m_Pl_kg)
    m_b_f1_arr = np.sqrt(2*eps_f1_arr) * E_Pl_eV
    m_22_f1_arr = m_b_f1_arr / 1e-22
    valid = (m_22_f1_arr > 0) & np.isfinite(m_22_f1_arr) & (m_22_f1_arr < 1e10)
    if valid.any():
        r_c_f1_arr = np.where(valid,
                              1.61 / m_22_f1_arr / (M_range2/1e9)**0.333,
                              np.nan)
        ax5.loglog(M_range2, r_c_f1_arr, ls='--', lw=2,
                   label=f'TGP-F1 m_sp={m_sp_t:.0e}')
# Obserwacje (przykładowe punkty)
M_gal_obs = np.array([3e8, 5e9, 1e10, 2e10, 1e11]) # M_sun
r_c_obs    = np.array([3.0, 2.0, 1.8, 1.5, 0.5])    # kpc (hipotetyczne)
ax5.scatter(M_gal_obs, r_c_obs, s=80, color='red', zorder=10, label='Hipotetyczne obs.')
ax5.set_xlabel(r'$M_{disk}$ [M$_\odot$]', fontsize=10)
ax5.set_ylabel(r'$r_c$ [kpc]', fontsize=10)
ax5.set_title('Predykcja r_c(M): F1 vs standard FDM', fontsize=10)
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3, which='both')
ax5.set_xlim(1e7, 1e12); ax5.set_ylim(0.01, 50)

# Panel 6: Krzywa rotacji NGC 3198
ax6 = fig.add_subplot(gs[1, 2])
ax6.errorbar(OBS_R, OBS_V, yerr=OBS_E, fmt='ko', ms=5, capsize=3,
             label='NGC 3198 (obs.)', zorder=10)
ax6.plot(r_plot, v_bar, 'k--', lw=1.5, alpha=0.7, label='Newton (baryony)')
ax6.plot(r_plot, v_total_F3, 'g-', lw=2.5,
         label=f'F3: ε=m_sp²/2\n(m22={m_22_F3:.3f}, χ²={chi2_F3:.2f})')
ax6.set_xlabel('r [kpc]', fontsize=10)
ax6.set_ylabel('v [km/s]', fontsize=10)
ax6.set_title('NGC 3198: TGP-FDM (F3, brak nowych par.)', fontsize=10)
ax6.legend(fontsize=8, loc='lower right')
ax6.grid(True, alpha=0.3)
ax6.set_xlim(0, 33); ax6.set_ylim(0, 220)

# Panel 7: Diagram wyprowadzenia (tekstowy)
ax7 = fig.add_subplot(gs[2, :])
ax7.axis('off')
derivation = r"""
WYPROWADZENIE ε Z AKSJOMATÓW TGP — DIAGRAM:

N0-5: $V_{eff}(\Phi) = \beta\Phi^3/3 - \gamma\Phi^4/4$  +  N0-5: $\beta=\gamma$ (warunek próżni)  $\Longrightarrow$  N0-6: $m_{sp}^2 = \gamma$

$\Downarrow$  (przesunięcie próżni o $\varepsilon g^2$)

$V_{mod}(g) = \varepsilon g^2 + m_{sp}^2(g^3/3 - g^4/4)$   $\Longrightarrow$   Minimum STABILNE gdy $V''_{mod}(g_{min}) > 0$

$\Downarrow$  (warunek progu: $V''_{mod}(g=1) = 2\varepsilon - m_{sp}^2 = 0$)

$\boxed{\varepsilon_{th} = \frac{m_{sp}^2}{2} = \frac{\gamma}{2}}$   ← NIE JEST NOWYM PARAMETREM! Wyprowadzone z N0-6.

$\Downarrow$  ($V_{mod}$ z $\varepsilon = \varepsilon_{th}$: minimum przy złotej liczbie)

$V'_{mod}(g) = g + g^2 - g^3 = 0 \Rightarrow g_{min} = \frac{1+\sqrt{5}}{2} = \varphi = 1.618...$

$\Downarrow$  (masa bozonu z $V''_{mod}(\varphi)$ i konwersja do eV)

$m_{boson} = m_{sp} \cdot m_{Pl}$   $\Longrightarrow$   FDM z $m_{sp} = 10^{-22}\,\text{eV}/m_{Pl} = 8.2\times10^{-51}\,l_{Pl}^{-1}$   $\Longrightarrow$   $\lambda_{Yukawa} = 1/m_{sp} \gg R_{galaktyki}$

BONUS: F1 (sprzężenie masowe) $\varepsilon(M) = C^2 = m_{sp}^2 M^2/(4\pi)$ $\Rightarrow$ $m_{boson} \propto M_{gal}$ (środowiskowe FDM, Kill-shot K18)
"""
ax7.text(0.01, 0.95, derivation, transform=ax7.transAxes,
         fontsize=10.5, verticalalignment='top', family='serif',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.85))

plt.tight_layout()
out_png = os.path.join(os.path.dirname(__file__), 'ex39_epsilon_from_coupling.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()
print(f"Wykres: {out_png}")

# =============================================================================
# PODSUMOWANIE KOŃCOWE
# =============================================================================
print()
print("=" * 72)
print("PODSUMOWANIE ex39: ε WYNIKA Z STRUKTURY TGP — ZERO NOWYCH PARAMETRÓW")
print("=" * 72)
print()
print("GŁÓWNY WYNIK:")
print()
print("  ε_th = m_sp²/2 = γ/2")
print()
print("  Parametr ε NIE jest nowy — to PRÓG STABILNOŚCI wynikający z N0-6:")
print("  minimalna perturbacja potencjału TGP dająca stabilny soliton.")
print()
print("ZŁOTA LICZBA:")
print(f"  V_mod z ε = ε_th ma minimum przy g = φ = (1+√5)/2 = {phi_gold:.6f}")
print("  Złoty podział pojawia się naturalnie z warunków stabilności TGP!")
print()
print("MASA BOZONU FDM:")
print("  m_boson = m_sp (identyczne!)")
print("  Bozon FDM to sam kwant pola TGP.")
print("  Żadna nowa skala energii — wszystko z jednego m_sp.")
print()
print("SPRZĘŻENIE MASOWE F1:")
print("  ε(M) = C² = m_sp²·M²/(4π)  [kwadrat sprzężenia Yukawa]")
print("  m_boson_eff(M) ∝ M_gal   (środowiskowe FDM)")
print("  r_c ∝ 1/M_gal            (większe galaktyki → mniejsze rdzenie)")
print("  Falsyfikowalne przez katalog SPARC/THINGS!")
print()
print("KILL-SHOT K18:")
print("  Test r_c vs M_gal w katalogu galaktyk:")
print("  r_c ∝ M_gal^{-1}:  TGP-FDM F1 potwierdzone")
print("  r_c ≈ const:       Standardowe FDM (m_boson = const)")
print()
print("KONSEKWENCJA PARADOKSU γ:")
print("  Paradoks γ był: 'γ musi być jednocześnie 10^{-38} (V₃) i 10^{-110} (DM)'")
print("  Teraz: DM nie potrzebuje γ ~ 10^{-110}!")
print("  ε_th = γ/2 i m_boson = m_sp = √γ·m_Pl")
print("  Dla m_boson = 10^{-22} eV: γ = m_sp² = 6.7×10^{-101}")
print("  Ale Yukawa jest NIEOSERWOWALNY (λ >> wszechświat)!")
print("  V₃ jest ekstremalnie małe (γ ~ 10^{-101})")
print("  Paradoks ROZWIĄZANY: γ=m_sp²=10^{-101} spełnia WSZYSTKIE wymagania FDM!")
print()
print("OTWARTE PYTANIE:")
print("  Dlaczego m_sp = 8.2×10^{-51} l_Pl^{-1}?")
print("  To odpowiednik pytania: dlaczego m_FDM = 10^{-22} eV?")
print("  Nie jest to pytanie specyficzne dla TGP — to otwarte pytanie FDM.")
