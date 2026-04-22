"""
ex37_multifit_galaxies.py
=========================
Test predykcyjności TGP-FDM: jeden parametr ε (→ m_22) dla wielu galaktyk.

PYTANIE KLUCZOWE:
  Jeżeli TGP-FDM (ex36) wyjaśnia NGC 3198 z m_22 ~ 1, to czy TA SAMA
  masa bozonu m_22 pasuje do INNYCH galaktyk?

  FDM jest predykcyjne: m_boson jest STAŁA natury, niezależna od galaktyki.
  Test: skanuj m_22 i szukaj wartości minimalizującej chi^2 dla WSZYSTKICH galaktyk jednocześnie.

GALAKTYKI TESTOWE (4):
  1. NGC 3198 — Begeman 1989        spiralna, v_flat ~ 150 km/s, R_d=3.2 kpc
  2. NGC 2403 — Begeman 1991        mniejsza spiralna, v_flat ~ 135 km/s, R_d=1.7 kpc
  3. NGC 6503 — de Blok et al 1996  mała spiralna, v_flat ~ 116 km/s, R_d=1.73 kpc
  4. DDO 154  — Carignan & Purton   karłowata, DM-zdominowana, v_flat ~ 47 km/s

HIPOTEZY:
  H1 (FDM universalne): m_22 = const, M_sol ∝ M_disk → jednoparametrowe skalowanie
  H2 (FDM galaktyczne): m_22 = const, M_sol wolna dla każdej galaktyki (więcej wolności)
  H3 (brak FDM): m_22 różna dla każdej galaktyki → nie jest stałą natury

SKALOWANIE FDM (Schive+2014, Chan+2018):
  r_c = 1.61 / m_22 / (M_sol/1e9)^{1/3}  [kpc]
  Relacja masa-promień solitonu dla galaktyki spiralnej:
    M_sol ~ 1.4e9 * m_22^{-2} * (sigma_v/10km/s)^{-4}  M_sun
  gdzie sigma_v = dyspersja prędkości gwiazd (∝ v_flat dla spiralnych).

OCZEKIWANE WYNIKI:
  Jeśli FDM prawdziwe: m_22 ~ 0.5–3 (zgodne z obserwacjami Lyman-alpha i struktur kpc)
  Jeśli TGP-FDM specyficzne: m_22 inna niż standardowa FDM
  Jeśli H3: TGP-FDM wymaga różnych ε dla każdej galaktyki → nie jest teorią

STATUS: ex37 = kluczowy test predykcyjności TGP-FDM.
"""

import sys, os
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.integrate import quad
from scipy.special import k0, k1, i0, i1
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# =============================================================================
# Stałe fizyczne
# =============================================================================
G_SI  = 6.674e-11
kpc_m = 3.086e19
M_sun = 1.989e30
km_s  = 1e3
pc_m  = 3.086e16

print("=" * 70)
print("EX37: TEST PREDYKCYJNOŚCI TGP-FDM NA CZTERECH GALAKTYKACH")
print("=" * 70)
print()

# =============================================================================
# Dane obserwacyjne galaktyk
# =============================================================================
GALAXIES = {

    'NGC 3198': {
        'ref': 'Begeman 1989',
        'type': 'Sb (spiralna)',
        'M_disk': 2.0e10,   # M_sun (maximum disk fit)
        'R_d':    3.2,      # kpc (skala dysku)
        'r_obs':  np.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0]),
        'v_obs':  np.array([90., 130., 155., 160., 155., 152., 150., 148., 147., 146., 145.]),
        'v_err':  np.array([ 8.,   6.,   5.,   5.,   4.,   4.,   4.,   4.,   5.,   6.,   7.]),
        'v_flat': 150.0,    # km/s
    },

    'NGC 2403': {
        'ref': 'Begeman 1991',
        'type': 'Scd (spiralna)',
        'M_disk': 1.0e10,   # M_sun
        'R_d':    1.7,      # kpc
        'r_obs':  np.array([0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0]),
        'v_obs':  np.array([55., 80., 110., 125., 130., 132., 133., 135., 135., 135., 134.]),
        'v_err':  np.array([ 6.,  5.,   5.,   5.,   5.,   4.,   4.,   4.,   5.,   6.,   7.]),
        'v_flat': 135.0,
    },

    'NGC 6503': {
        'ref': 'de Blok et al 1996',
        'type': 'Sc (mała spiralna)',
        'M_disk': 5.0e9,    # M_sun
        'R_d':    1.73,     # kpc
        'r_obs':  np.array([0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 14.0]),
        'v_obs':  np.array([30., 60., 85., 100., 108., 112., 114., 116., 116., 116.]),
        'v_err':  np.array([ 5.,  4.,  4.,   4.,   4.,   3.,   3.,   3.,   4.,   5.]),
        'v_flat': 116.0,
    },

    'DDO 154': {
        'ref': 'Carignan & Purton 1998',
        'type': 'karłowata (DM-zdominowana)',
        'M_disk': 3.0e8,    # M_sun
        'R_d':    0.8,      # kpc
        'r_obs':  np.array([0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        'v_obs':  np.array([10., 22., 33., 42., 45., 47., 47., 47.]),
        'v_err':  np.array([ 4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.]),
        'v_flat': 47.0,
    },
}

# =============================================================================
# Profil FDM solitonu (Schive+2014) + dysk baryonowy
# =============================================================================

def v_disk_freeman(r_kpc_arr, M_disk_Msun, R_d_kpc):
    """Prędkość rotacji dysku eksponencjalnego (Freeman 1970)."""
    M_disk_kg = M_disk_Msun * M_sun
    R_d_m     = R_d_kpc * kpc_m
    Sigma0 = M_disk_kg / (2.0 * np.pi * R_d_m**2)
    y = r_kpc_arr * kpc_m / (2.0 * R_d_m)
    y = np.maximum(y, 1e-10)
    v2 = (4.0 * np.pi * G_SI * Sigma0 * R_d_m
          * y**2 * (i0(y)*k0(y) - i1(y)*k1(y)))
    return np.sqrt(np.maximum(v2, 0.0)) / km_s

def rc_from_m22_Msol(m_22, M_sol_Msun):
    """r_c = 1.61 / m_22 / (M_sol/1e9)^(1/3) [kpc] — Schive+2014."""
    return 1.61 / m_22 / (M_sol_Msun / 1.0e9)**(1.0/3.0)

def rhoc_from_rc_Msol(r_c_kpc, M_sol_Msun):
    """Gęstość centralna solitonu z r_c i M_sol (normalizacja Schive+2014)."""
    # I = 4π ∫₀^∞ r²/(1+0.091r²)^8 dr ≈ 0.3526 (w j. r_c)
    I_norm = 0.3526
    rho_c_Msun_kpc3 = M_sol_Msun / (4.0 * np.pi * I_norm * r_c_kpc**3)
    return rho_c_Msun_kpc3 / 1.0e9  # → M_sun/pc^3

def M_soliton_cumulative(r_kpc_arr, rho_c_Msun_pc3, r_c_kpc):
    """Kumulatywna masa solitonu M(<r) dla tablicy r [M_sun]."""
    M_arr = np.zeros(len(r_kpc_arr))
    for i, r_val in enumerate(r_kpc_arr):
        if r_val < 1e-6:
            continue
        r_int = np.linspace(0.0, r_val, 300)
        x = r_int / r_c_kpc
        rho_int = rho_c_Msun_pc3 * 1.0e9 / (1.0 + 0.091*x**2)**8  # M_sun/kpc^3
        M_arr[i] = 4.0 * np.pi * np.trapz(rho_int * r_int**2, r_int)
    return M_arr

def v_soliton_arr(r_kpc_arr, rho_c_Msun_pc3, r_c_kpc):
    """Prędkość rotacji od solitonu FDM [km/s]."""
    M_kg_arr = M_soliton_cumulative(r_kpc_arr, rho_c_Msun_pc3, r_c_kpc) * M_sun
    r_m_arr  = r_kpc_arr * kpc_m
    v2 = G_SI * M_kg_arr / r_m_arr
    return np.sqrt(np.maximum(v2, 0.0)) / km_s

def compute_rotation_curve_fdm(r_kpc_arr, M_disk_Msun, R_d_kpc,
                                m_22, M_sol_Msun):
    """Pełna krzywa rotacji: dysk Freeman + soliton Schive+2014 [km/s]."""
    v_bar = v_disk_freeman(r_kpc_arr, M_disk_Msun, R_d_kpc)
    r_c   = rc_from_m22_Msol(m_22, M_sol_Msun)
    rho_c = rhoc_from_rc_Msol(r_c, M_sol_Msun)

    if r_c < 0.01 or r_c > 500 or rho_c <= 0 or not np.isfinite(r_c):
        return v_bar  # fallback do czystego baryonu

    v_sol = v_soliton_arr(r_kpc_arr, rho_c, r_c)
    return np.sqrt(v_bar**2 + v_sol**2)

def chi2_galaxy(m_22, M_sol_Msun, gal_data, r_model_kpc=None):
    """Chi^2/N dopasowania do obserwacji galaktyki."""
    if r_model_kpc is None:
        r_model_kpc = np.linspace(0.1, gal_data['r_obs'].max() * 1.1, 120)

    v_model = compute_rotation_curve_fdm(
        r_model_kpc, gal_data['M_disk'], gal_data['R_d'], m_22, M_sol_Msun)
    v_interp = np.interp(gal_data['r_obs'], r_model_kpc, v_model)
    return np.sum(((v_interp - gal_data['v_obs']) / gal_data['v_err'])**2) / len(gal_data['v_obs'])

# =============================================================================
# Sekcja 1: Optymalne (m_22, M_sol) dla każdej galaktyki oddzielnie
# =============================================================================
print("=" * 70)
print("SEKCJA 1: Optymalne parametry FDM dla każdej galaktyki oddzielnie")
print("=" * 70)
print()

m22_grid   = np.logspace(-1, 1.5, 25)    # 0.1 → ~32
Msol_grid  = np.logspace(8, 12, 25)      # 1e8 → 1e12 M_sun

best_per_galaxy = {}

for gal_name, gal in GALAXIES.items():
    r_model = np.linspace(0.1, gal['r_obs'].max() * 1.2, 100)
    best_c2 = np.inf
    best_m22 = 1.0
    best_Msol = 1e10

    for m22 in m22_grid:
        for Msol in Msol_grid:
            try:
                c2 = chi2_galaxy(m22, Msol, gal, r_model)
                if np.isfinite(c2) and c2 < best_c2:
                    best_c2 = c2
                    best_m22 = m22
                    best_Msol = Msol
            except Exception:
                pass

    best_rc = rc_from_m22_Msol(best_m22, best_Msol)
    best_per_galaxy[gal_name] = {
        'm_22': best_m22, 'M_sol': best_Msol,
        'r_c': best_rc, 'chi2': best_c2
    }

    print(f"  {gal_name}  ({gal['type']}):")
    print(f"    m_22   = {best_m22:.3f}   (m_boson = {best_m22*1e-22:.2e} eV)")
    print(f"    M_sol  = {best_Msol:.2e} M_sun")
    print(f"    r_c    = {best_rc:.2f} kpc")
    print(f"    chi2/N = {best_c2:.2f}")
    print()

# Sprawdź rozpiętość m_22 między galaktykami
m22_vals = [v['m_22'] for v in best_per_galaxy.values()]
m22_spread = max(m22_vals) / min(m22_vals)
print(f"  Rozpiętość m_22: {min(m22_vals):.3f} – {max(m22_vals):.3f}")
print(f"  Stosunek max/min: {m22_spread:.1f}x")
print(f"  Rozpiętość w log: {np.log10(m22_spread):.2f} dekad")
print()
is_universal = m22_spread < 5.0
print(f"  Czy m_22 jest 'universal' (< 5x rozpiętość)?",
      "TAK ✓" if is_universal else "NIE ✗")
print()

# =============================================================================
# Sekcja 2: Hipoteza H1 — jeden m_22 dla wszystkich galaktyk
# =============================================================================
print("=" * 70)
print("SEKCJA 2: Hipoteza H1 — wspólne m_22 (stała natury)")
print("=" * 70)
print()

def chi2_all_galaxies_fixed_m22(m_22, M_sol_free=True):
    """
    Chi^2 sumowane po wszystkich galaktykach dla ustalonego m_22.
    Jeśli M_sol_free=True: dla każdej galaktyki optymalizujemy M_sol.
    """
    total_chi2 = 0.0
    n_pts = 0

    for gal_name, gal in GALAXIES.items():
        r_model = np.linspace(0.1, gal['r_obs'].max() * 1.2, 100)

        if M_sol_free:
            # Znajdź najlepsze M_sol dla tej galaktyki przy danym m_22
            def chi2_Msol(log_Msol):
                M_sol = 10**log_Msol
                try:
                    return chi2_galaxy(m_22, M_sol, gal, r_model)
                except Exception:
                    return 1e10

            # Skanuj M_sol
            log_Msol_arr = np.linspace(8, 12, 20)
            c2_arr = [chi2_Msol(lm) for lm in log_Msol_arr]
            best_idx = np.argmin(c2_arr)
            best_chi2_gal = c2_arr[best_idx]
        else:
            # Użyj skalowania M_sol ∝ M_disk
            M_sol = gal['M_disk'] * 0.5  # heurystyczne skalowanie
            best_chi2_gal = chi2_galaxy(m_22, M_sol, gal, r_model)

        total_chi2 += best_chi2_gal * len(gal['v_obs'])
        n_pts += len(gal['v_obs'])

    return total_chi2 / n_pts

print("  Skan m_22 → łączne chi²/N (M_sol optymalizowane per galaktyka):")
print()

m22_scan = np.logspace(-1.5, 2.0, 40)
chi2_scan_total = []
for m22 in m22_scan:
    c2 = chi2_all_galaxies_fixed_m22(m22, M_sol_free=True)
    chi2_scan_total.append(c2)
    print(f"  m_22 = {m22:.3f}  chi2_total/N = {c2:.2f}")

chi2_scan_arr = np.array(chi2_scan_total)
idx_best = np.argmin(chi2_scan_arr)
m22_universal = m22_scan[idx_best]
chi2_universal = chi2_scan_arr[idx_best]

print()
print(f"  NAJLEPSZE WSPÓLNE m_22 = {m22_universal:.3f}")
print(f"  m_boson = {m22_universal*1e-22:.3e} eV")
print(f"  Łączne chi2/N = {chi2_universal:.2f}")
print()

# Sprawdź per-galaktyka dla m22_universal
print(f"  Dopasowanie per galaktyka dla m_22 = {m22_universal:.3f}:")
print(f"  {'Galaktyka':15s}  {'M_sol opty. [M_sun]':22s}  {'r_c [kpc]':10s}  {'chi2/N':8s}")
print("-" * 62)

best_Msol_for_universal = {}
for gal_name, gal in GALAXIES.items():
    r_model = np.linspace(0.1, gal['r_obs'].max() * 1.2, 100)
    log_Msol_arr = np.linspace(8, 12, 30)
    c2_arr = []
    for lm in log_Msol_arr:
        M_sol = 10**lm
        try:
            c2 = chi2_galaxy(m22_universal, M_sol, gal, r_model)
        except Exception:
            c2 = 1e10
        c2_arr.append(c2)
    idx = np.argmin(c2_arr)
    M_sol_opt = 10**log_Msol_arr[idx]
    chi2_opt  = c2_arr[idx]
    r_c_opt   = rc_from_m22_Msol(m22_universal, M_sol_opt)
    best_Msol_for_universal[gal_name] = M_sol_opt
    print(f"  {gal_name:15s}  {M_sol_opt:22.2e}  {r_c_opt:10.2f}  {chi2_opt:8.2f}")

# =============================================================================
# Sekcja 3: Hipoteza H1b — M_sol ∝ M_disk (jedno wolne skalowanie)
# =============================================================================
print()
print("=" * 70)
print("SEKCJA 3: Hipoteza H1b — M_sol = f * M_disk (2 parametry globalne)")
print("=" * 70)
print()
print("  FDM przewiduje: M_sol ~ f * M_disk (masa solitonu proporcjonalna")
print("  do masy dysku baryonowego). Testujemy ten ansatz.")
print()

def chi2_all_h1b(params):
    """Chi^2 dla modelu H1b: m_22 globalne, M_sol = f * M_disk."""
    log_m22, log_f = params
    m22 = 10**log_m22
    f   = 10**log_f

    total_c2 = 0.0
    n_pts = 0
    for gal_name, gal in GALAXIES.items():
        M_sol = f * gal['M_disk']
        if M_sol < 1e7 or M_sol > 1e13:
            return 1e10
        r_model = np.linspace(0.1, gal['r_obs'].max() * 1.2, 100)
        try:
            c2 = chi2_galaxy(m22, M_sol, gal, r_model)
        except Exception:
            return 1e10
        total_c2 += c2 * len(gal['v_obs'])
        n_pts += len(gal['v_obs'])
    return total_c2 / n_pts

# Skanuj siatkę 2D
log_m22_grid = np.linspace(-1, 1.5, 15)
log_f_grid   = np.linspace(-2, 2, 15)
chi2_h1b = np.zeros((len(log_m22_grid), len(log_f_grid)))

for i, lm22 in enumerate(log_m22_grid):
    for j, lf in enumerate(log_f_grid):
        chi2_h1b[i, j] = chi2_all_h1b([lm22, lf])

best_ij = np.unravel_index(np.nanargmin(chi2_h1b), chi2_h1b.shape)
m22_h1b = 10**log_m22_grid[best_ij[0]]
f_h1b   = 10**log_f_grid[best_ij[1]]
chi2_h1b_best = chi2_h1b[best_ij]

print(f"  Wynik H1b (m_22 globalne + M_sol = f * M_disk):")
print(f"    m_22 = {m22_h1b:.3f}")
print(f"    f    = {f_h1b:.3f}  (M_sol / M_disk)")
print(f"    chi2_total/N = {chi2_h1b_best:.2f}")
print()

# Per-galaktyka dla H1b
print(f"  Dopasowanie per galaktyka dla H1b (m_22={m22_h1b:.3f}, f={f_h1b:.3f}):")
print(f"  {'Galaktyka':15s}  {'M_sol [M_sun]':16s}  {'r_c [kpc]':10s}  {'chi2/N':8s}  {'v(v_flat)':10s}")
print("-" * 65)
for gal_name, gal in GALAXIES.items():
    M_sol = f_h1b * gal['M_disk']
    r_c   = rc_from_m22_Msol(m22_h1b, M_sol)
    r_model = np.linspace(0.1, gal['r_obs'].max() * 1.2, 100)
    try:
        c2 = chi2_galaxy(m22_h1b, M_sol, gal, r_model)
    except Exception:
        c2 = 999.0
    v_at_flat = np.interp(gal['r_obs'].max(), r_model,
                          compute_rotation_curve_fdm(r_model, gal['M_disk'],
                                                      gal['R_d'], m22_h1b, M_sol))
    print(f"  {gal_name:15s}  {M_sol:16.2e}  {r_c:10.2f}  {c2:8.2f}  "
          f"{v_at_flat:8.1f}/{gal['v_flat']:.0f}")

# =============================================================================
# Sekcja 4: Masa solitonu vs masa dysku — relacja skalowania
# =============================================================================
print()
print("=" * 70)
print("SEKCJA 4: Relacja skalowania M_sol vs M_disk")
print("=" * 70)
print()

M_disk_vals = [GALAXIES[g]['M_disk'] for g in GALAXIES]
M_sol_vals  = [best_Msol_for_universal[g] for g in GALAXIES]
names = list(GALAXIES.keys())

print(f"  {'Galaktyka':15s}  {'M_disk [M_sun]':16s}  {'M_sol opt. [M_sun]':20s}  {'M_sol/M_disk':12s}")
print("-" * 68)
for name, Md, Ms in zip(names, M_disk_vals, M_sol_vals):
    ratio = Ms / Md
    print(f"  {name:15s}  {Md:16.2e}  {Ms:20.2e}  {ratio:12.3f}")

# Fit relacji M_sol ~ M_disk^alpha
log_Md = np.log10(M_disk_vals)
log_Ms = np.log10(M_sol_vals)
if len(log_Md) > 1:
    alpha, log_A = np.polyfit(log_Md, log_Ms, 1)
    A = 10**log_A
    print()
    print(f"  Relacja skalowania: M_sol ~ {A:.3e} * M_disk^{alpha:.3f}")
    print(f"  FDM teoria: M_sol ~ M_disk (alpha=1) lub M_sol ∝ M_disk^(1/3) (jeśli r_c=const)")
    print()
    # Sprawdź czy alpha ~ 1 (proporcjonalność)
    if 0.5 < alpha < 1.5:
        print("  => Relacja quasi-liniowa (alpha ~ 1): SPÓJNE z FDM")
    elif 0.1 < alpha < 0.5:
        print("  => Relacja sublinearna: M_sol rośnie wolniej niż M_disk")
    else:
        print(f"  => Relacja nielinearna (alpha={alpha:.2f}): ANOMALNY")

# =============================================================================
# Sekcja 5: Predykcja Drogi Mlecznej
# =============================================================================
print()
print("=" * 70)
print("SEKCJA 5: Predykcja dla Drogi Mlecznej (Milky Way)")
print("=" * 70)
print()

# Parametry Drogi Mlecznej (Bovy 2015)
MW_data = {
    'M_disk': 6.0e10,   # M_sun (baryony + gaz)
    'R_d':    2.5,      # kpc
    'v_LSR':  220.0,    # km/s (Local Standard of Rest)
    'R_sun':  8.2,      # kpc
    'r_obs':  np.array([3.0, 5.0, 8.2, 12.0, 15.0, 20.0, 25.0]),
    'v_obs':  np.array([200., 215., 220., 222., 225., 230., 232.]),
    'v_err':  np.array([10.,   8.,   6.,   8.,   10.,  15.,  20.]),
}

print(f"  Droga Mleczna (Bovy 2015):")
print(f"    M_disk = {MW_data['M_disk']:.1e} M_sun, R_d = {MW_data['R_d']} kpc")
print(f"    v(R_sun={MW_data['R_sun']} kpc) = {MW_data['v_LSR']} km/s")
print()

# Predykcja FDM dla MW używając m_22_universal
M_sol_MW_H1b = f_h1b * MW_data['M_disk']
M_sol_MW_best = best_Msol_for_universal.get('NGC 3198', 1e10)  # fallback

r_MW = np.linspace(0.5, 30.0, 120)

v_MW_bar = v_disk_freeman(r_MW, MW_data['M_disk'], MW_data['R_d'])

v_MW_fdm_H1b = compute_rotation_curve_fdm(
    r_MW, MW_data['M_disk'], MW_data['R_d'],
    m22_h1b, M_sol_MW_H1b)

v_MW_fdm_free = compute_rotation_curve_fdm(
    r_MW, MW_data['M_disk'], MW_data['R_d'],
    m22_universal, 2.0 * MW_data['M_disk'])  # M_sol ~ 2*M_disk jako przykład

v_RLSR_bar = np.interp(MW_data['R_sun'], r_MW, v_MW_bar)
v_RLSR_H1b = np.interp(MW_data['R_sun'], r_MW, v_MW_fdm_H1b)
v_RLSR_free = np.interp(MW_data['R_sun'], r_MW, v_MW_fdm_free)

print(f"  v(R_sun = {MW_data['R_sun']} kpc):")
print(f"    Obserwacja:     {MW_data['v_LSR']:.0f} km/s")
print(f"    Newton only:    {v_RLSR_bar:.0f} km/s")
print(f"    FDM H1b:        {v_RLSR_H1b:.0f} km/s  (M_sol={M_sol_MW_H1b:.2e} M_sun, r_c={rc_from_m22_Msol(m22_h1b, M_sol_MW_H1b):.1f} kpc)")
print(f"    FDM wolne:      {v_RLSR_free:.0f} km/s  (M_sol=2*M_disk)")
print()

# Czy MW jest spójna z m22_universal?
c2_MW_H1b = np.sum(((np.interp(MW_data['r_obs'], r_MW, v_MW_fdm_H1b) - MW_data['v_obs']) / MW_data['v_err'])**2) / len(MW_data['v_obs'])
print(f"  Chi²/N dla MW (H1b, m_22={m22_h1b:.3f}): {c2_MW_H1b:.2f}")
print()

# =============================================================================
# Sekcja 6: Test Universalności — werdykt
# =============================================================================
print("=" * 70)
print("SEKCJA 6: WERDYKT — Universalność FDM w TGP")
print("=" * 70)
print()

print("  TABELA PORÓWNAWCZA (m_22 optymalne indywidualnie vs globalnie):")
print()
print(f"  {'Galaktyka':15s}  {'m_22 indyw.':12s}  {'m_22 global':12s}  "
      f"{'chi2 indyw.':12s}  {'chi2 global':12s}  {'Koszt H1':8s}")
print("-" * 78)

chi2_global_per_gal = {}
for gal_name, gal in GALAXIES.items():
    r_model = np.linspace(0.1, gal['r_obs'].max() * 1.2, 100)
    M_sol_opt = best_Msol_for_universal[gal_name]
    c2_g = chi2_galaxy(m22_universal, M_sol_opt, gal, r_model)
    chi2_global_per_gal[gal_name] = c2_g

    m22_indiv = best_per_galaxy[gal_name]['m_22']
    c2_indiv  = best_per_galaxy[gal_name]['chi2']
    cost = c2_g / max(c2_indiv, 0.01) - 1.0

    print(f"  {gal_name:15s}  {m22_indiv:12.3f}  {m22_universal:12.3f}  "
          f"{c2_indiv:12.2f}  {c2_g:12.2f}  {cost:+8.1%}")

print()

# Ocena spójności
chi2_vals_indiv = [best_per_galaxy[g]['chi2'] for g in GALAXIES]
chi2_vals_global = list(chi2_global_per_gal.values())

all_ok_indiv  = all(c < 5.0 for c in chi2_vals_indiv)
all_ok_global = all(c < 5.0 for c in chi2_vals_global)
some_bad_global = any(c > 10.0 for c in chi2_vals_global)

print("  DIAGNOZA:")
print()
print(f"  m_22 indywidualnie < 5× różnica: {is_universal}")
print(f"  Wszystkie galaktyki OK (chi²<5) indywidualnie: {all_ok_indiv}")
print(f"  Wszystkie galaktyki OK (chi²<5) dla m_22={m22_universal:.3f}: {all_ok_global}")
print(f"  Jakaś galaktyka bardzo zła (chi²>10) globalnie: {some_bad_global}")
print()

if all_ok_global and is_universal:
    verdict_h1 = "TAK ✓ — TGP-FDM jest PREDYKCYJNE"
    detail = (f"m_22 ~ {m22_universal:.2f} pasuje do wszystkich 4 galaktyk\n"
              f"  przy optymalnym M_sol per galaktyka.")
elif all_ok_indiv and not all_ok_global:
    verdict_h1 = "CZĘŚCIOWO — m_22 różni się między galaktykami"
    detail = ("m_22 musi być dopasowane per galaktyka.\n"
              "  FDM TGP może nie być prawdziwą stałą natury.")
else:
    verdict_h1 = "NIE ✗ — TGP-FDM nie jest predykcyjne globalnie"
    detail = "Brak wspólnego m_22 spójnego ze wszystkimi galaktykami."

print(f"  HIPOTEZA H1 (m_22 = stała natury): {verdict_h1}")
print(f"  {detail}")
print()

# =============================================================================
# Sekcja 7: Ograniczenie Lyman-alpha na m_22
# =============================================================================
print("=" * 70)
print("SEKCJA 7: Ograniczenia zewnętrzne na m_22")
print("=" * 70)
print()
print("  Niezależne obserwacyjne ograniczenia na m_boson FDM:")
print()
print("  [1] Lyman-alpha forest (Irsic+2017, Armengaud+2017):")
print("      m_22 > 10–20  (wyklucza zbyt lekkie bozony)")
print("      => m_boson > 1.5×10^{-21} eV")
print()
print("  [2] Galaktyki karłowate MW + podstruktura halo (Schutz 2020):")
print("      m_22 > 2.9   (Droga Mleczna satellites)")
print()
print("  [3] Substruktura w halo MW (Safarzadeh & Spergel 2020):")
print("      m_22 > 0.4–2.0")
print()
print("  [4] Profil solitonu w galaktykach karłowatych Fornax, Sculptor:")
print("      m_22 ~ 0.5–5.0 (DOPASOWANE do profilu)")
print()
print("  [5] Semi-analityczne modele galaktyk (Lovell+2021):")
print("      m_22 ~ 1–3 (najlepsze dopasowanie do SDSS)")
print()

print(f"  Wynik ex37: m_22_universal = {m22_universal:.3f}")
print()

lya_lower = 10.0  # Lyman-alpha bound (konserwatywny)
satellite_lower = 2.9
soliton_range = (0.5, 5.0)

in_lya = m22_universal > lya_lower
in_satellite = m22_universal > satellite_lower
in_soliton = soliton_range[0] <= m22_universal <= soliton_range[1]

print(f"  Zgodność z ograniczeniami:")
print(f"    Lyman-alpha (m_22 > {lya_lower:.0f}):   {'✓' if in_lya else '✗'}  "
      f"(m_22={m22_universal:.2f} {'>' if in_lya else '<'} {lya_lower:.0f})")
print(f"    MW satellites (m_22 > {satellite_lower:.1f}): {'✓' if in_satellite else '✗'}  "
      f"(m_22={m22_universal:.2f} {'>' if in_satellite else '<'} {satellite_lower:.1f})")
print(f"    Soliton profil ({soliton_range[0]:.1f}–{soliton_range[1]:.1f}): "
      f"{'✓' if in_soliton else '✗'}")
print()

if in_lya:
    lya_msg = "SPÓJNE z Lyman-alpha"
else:
    lya_msg = f"NIEZGODNE z Lyman-alpha (m_22 zbyt małe o {np.log10(lya_lower/m22_universal):.1f} dekad)"
print(f"  WNIOSEK: {lya_msg}")
print()

# =============================================================================
# Wykresy
# =============================================================================
print("=" * 70)
print("Generowanie wykresów...")

fig = plt.figure(figsize=(20, 14))
gs = GridSpec(3, 3, figure=fig, hspace=0.40, wspace=0.35)
fig.suptitle('TGP-FDM: Test Universalności na 4 Galaktykach (ex37)\n'
             'Profil solitonu Schive+2014, dysk Freeman 1970',
             fontsize=13, y=1.01)

# --- 4 krzywe rotacji ---
gal_names = list(GALAXIES.keys())
gal_list  = list(GALAXIES.values())
axes_rcs  = [fig.add_subplot(gs[i//2, i%2]) for i in range(4)]

gal_colors = ['#1976D2', '#388E3C', '#F57C00', '#7B1FA2']

for idx, (ax, gal_name, gal, col) in enumerate(zip(axes_rcs, gal_names, gal_list, gal_colors)):
    r_plot = np.linspace(0.1, gal['r_obs'].max() * 1.2, 150)

    # Dane obs.
    ax.errorbar(gal['r_obs'], gal['v_obs'], yerr=gal['v_err'],
                fmt='o', color=col, ms=5, capsize=3,
                label=f"{gal_name}\n({gal['ref']})", zorder=10)

    # Newton
    v_n = v_disk_freeman(r_plot, gal['M_disk'], gal['R_d'])
    ax.plot(r_plot, v_n, 'k--', lw=1.5, alpha=0.6, label='Newton')

    # FDM indywidualne (najlepsze)
    bp = best_per_galaxy[gal_name]
    v_fdm_ind = compute_rotation_curve_fdm(r_plot, gal['M_disk'], gal['R_d'],
                                            bp['m_22'], bp['M_sol'])
    ax.plot(r_plot, v_fdm_ind, '-', color=col, lw=2.5,
            label=f"FDM indyw. m22={bp['m_22']:.2f}\n(χ²={bp['chi2']:.1f})")

    # FDM globalne
    M_sol_g = best_Msol_for_universal[gal_name]
    v_fdm_g = compute_rotation_curve_fdm(r_plot, gal['M_disk'], gal['R_d'],
                                          m22_universal, M_sol_g)
    c2_g = chi2_global_per_gal[gal_name]
    ax.plot(r_plot, v_fdm_g, '--', color=col, lw=1.8, alpha=0.8,
            label=f"FDM glob. m22={m22_universal:.2f}\n(χ²={c2_g:.1f})")

    ax.set_xlabel('r [kpc]', fontsize=10)
    ax.set_ylabel('v [km/s]', fontsize=10)
    ax.set_title(gal_name, fontsize=11, color=col, fontweight='bold')
    ax.legend(fontsize=7, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

# --- chi²(m_22) scan ---
ax5 = fig.add_subplot(gs[2, 0])
ax5.semilogx(m22_scan, chi2_scan_arr, 'k-', lw=2.5, label='χ² globalny (M_sol wolne)')
ax5.axvline(m22_universal, color='r', ls='--', lw=2,
            label=f'opt. m_22={m22_universal:.2f}')
ax5.axhline(2.0, color='gold', ls=':', lw=1.5, alpha=0.8, label='χ²=2')
ax5.axvline(10, color='gray', ls='-.', lw=1.2, alpha=0.7, label='Lyman-α limit')
ax5.axvline(2.9, color='purple', ls=':', lw=1.2, alpha=0.7, label='MW sat. limit')
ax5.set_xlabel(r'$m_{22}$ (globalne)', fontsize=11)
ax5.set_ylabel(r'$\chi^2/N$ (suma 4 galaktyk)', fontsize=11)
ax5.set_title(r'Skan globalny $m_{22}$', fontsize=11)
ax5.legend(fontsize=7)
ax5.grid(True, alpha=0.3, which='both')
ax5.set_ylim(0, 30)

# --- M_sol vs M_disk ---
ax6 = fig.add_subplot(gs[2, 1])
M_disk_arr = np.array([gal['M_disk'] for gal in gal_list])
M_sol_arr  = np.array([best_Msol_for_universal[gn] for gn in gal_names])

for gn, Md, Ms, col in zip(gal_names, M_disk_arr, M_sol_arr, gal_colors):
    ax6.scatter(Md, Ms, color=col, s=150, zorder=10, label=gn, edgecolors='k')

# Linia fit
if len(M_disk_arr) > 1:
    x_fit = np.logspace(np.log10(M_disk_arr.min()*0.5),
                         np.log10(M_disk_arr.max()*2), 50)
    y_fit = A * x_fit**alpha
    ax6.loglog(x_fit, y_fit, 'r--', lw=2, label=rf'Fit: $M_{{sol}} \propto M_{{disk}}^{{{alpha:.2f}}}$')

ax6.loglog([1e8, 1e12], [1e8, 1e12], 'k:', lw=1, alpha=0.5, label=r'$M_{sol}=M_{disk}$')
ax6.set_xlabel(r'$M_{disk}$ [M$_\odot$]', fontsize=11)
ax6.set_ylabel(r'$M_{sol}$ [M$_\odot$] (optymalne)', fontsize=11)
ax6.set_title('Skalowanie masy solitonu', fontsize=11)
ax6.legend(fontsize=7)
ax6.grid(True, alpha=0.3, which='both')

# --- Pasek ograniczeń m_22 ---
ax7 = fig.add_subplot(gs[2, 2])
constraints = [
    ('Lyman-α (Irsic+17)', 10.0, 1000.0, '#FF5252', 0.3),
    ('MW satellites (Schutz+20)', 2.9, 1000.0, '#FF9800', 0.3),
    ('Soliton profil (Fornax)', 0.5, 5.0, '#4CAF50', 0.4),
    ('LSB galaktyki (Chen+17)', 0.2, 5.0, '#2196F3', 0.3),
]
y_positions = np.arange(len(constraints))
for i, (label, low, high, col, alpha_val) in enumerate(constraints):
    ax7.barh(i, min(high, 100) - low, left=low, height=0.6,
             color=col, alpha=alpha_val, edgecolor='k')
    ax7.text(max(low, 0.08), i, f'  {label}', va='center', fontsize=8)

ax7.axvline(m22_universal, color='r', lw=3, ls='--', zorder=10,
            label=f'TGP-FDM m_22={m22_universal:.2f}')
ax7.set_xscale('log')
ax7.set_xlim(0.07, 100)
ax7.set_yticks([])
ax7.set_xlabel(r'$m_{22}$', fontsize=11)
ax7.set_title('Ograniczenia obserwacyjne FDM', fontsize=11)
ax7.legend(fontsize=9)
ax7.grid(True, alpha=0.3, axis='x', which='both')

out_png = os.path.join(os.path.dirname(__file__), 'ex37_multifit_galaxies.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()
print(f"Wykres zapisany: {out_png}")

# =============================================================================
# Podsumowanie finalne
# =============================================================================
print()
print("=" * 70)
print("PODSUMOWANIE ex37: TEST UNIVERSALNOŚCI TGP-FDM")
print("=" * 70)
print()
print("WYNIK 1: Optymalne m_22 per galaktyka:")
for gn in gal_names:
    bp = best_per_galaxy[gn]
    print(f"  {gn:15s}: m_22 = {bp['m_22']:.3f}  (chi2={bp['chi2']:.2f})")
print()
print(f"WYNIK 2: Globalne m_22 = {m22_universal:.3f}  (łączne chi2/N = {chi2_universal:.2f})")
print()
print(f"WYNIK 3: H1b (M_sol ∝ M_disk): m_22={m22_h1b:.3f}, f={f_h1b:.3f}, chi2={chi2_h1b_best:.2f}")
print()
print("WYNIK 4: Relacja skalowania M_sol vs M_disk:")
print(f"  M_sol ~ {A:.2e} * M_disk^{alpha:.2f}")
print(f"  (FDM teoria: alpha=1 dla proportional scaling)")
print()
print("WYNIK 5: Zgodność z obserwacyjnymi ograniczeniami:")
print(f"  Lyman-alpha (m_22>10): {'✓' if in_lya else '✗'}")
print(f"  MW satellites (m_22>2.9): {'✓' if in_satellite else '✗'}")
print(f"  Soliton profil (0.5<m_22<5): {'✓' if in_soliton else '✗'}")
print()
print("WNIOSEK KOŃCOWY:")
print()
if is_universal and (not some_bad_global):
    print("  TGP-FDM jest PREDYKCYJNE — jeden m_22 pasuje do wszystkich galaktyk.")
    print("  Ale m_22 optimal może być niezgodne z Lyman-alpha.")
    print("  Jeśli m_22 < 10: TGP-FDM wykluczone przez Lyman-alpha.")
    print("  Potrzeba pełnej analizy kosmologicznej (power spectrum).")
elif not is_universal:
    print("  m_22 różni się >5× między galaktykami.")
    print("  TGP-FDM wymaga różnego ε dla każdej galaktyki.")
    print("  => NIE jest prawdziwą stałą natury w standardowym sensie.")
    print("  => Alternatywnie: m_22 zależy od środowiska (mechanizm chameleonowy).")
else:
    print("  TGP-FDM CZĘŚCIOWO predykcyjne — m_22 rozpiętość akceptowalna,")
    print("  ale niektóre galaktyki słabo dopasowane globalnie.")
    print("  Wymagana dalsza analiza.")
print()
print("NASTĘPNE KROKI:")
print("  ex38: Mechanizm spinodalny — czy V''<0 halo daje struktury galaktyczne?")
print("  ex39: Kosmologiczne widmo mocy TGP-FDM (vs ΛCDM, Lyman-alpha)")
print("  Teoria: skąd pochodzi m_22 z aksjomatów N0? (V_mod = εg²+V_TGP)")
