#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ps5_global_fit.py
=================

Program P5 - problem #5: globalny fit formuly T_c do 19 znanych SC,
z trzema wariantami wolnosci parametrow.

Motywacja:
  ps3 dal tylko r(log-log) = 0.306 bo uzywalo g_0 = g_0^e dla wszystkich
  materialow + Lambda_E skalibrowane tylko z Al. To zbyt upraszczajace.

  ps4 uzylo R_orb = (A_orb/A_s)^2 z mapowaniem orbitali na solitony P4:
     s <-> e-soliton (A_e  = 0.125)
     d <-> mu-soliton (A_mu = 0.472)   <- ZALOZENIE
     f <-> tau-soliton (A_tau = 0.956) <- ZALOZENIE

  Pytanie: czy to mapowanie pasuje ilosciowo do obserwacji T_c?

Strategia ps5 (3 warianty):
  5a (PURE-TGP):  fix A_s, A_d, A_f z P4, fit tylko {Lambda_E, sigma_a}
                   Sprawdza czy P4-P5 mapowanie daje sensowne T_c.
  5b (SP-FREE):   fix A_s, A_d, A_f z P4, fit {Lambda_E, A_sp, sigma_a}
                   Sprawdza czy sp-hybryd wymaga osobnej amplitudy.
  5c (ALL-FREE):  fit {Lambda_E, A_s, A_sp, A_d, A_f, sigma_a}
                   Porownanie z P4: czy fitowane A_d, A_f pasuja do mu/tau-solitonow?

Kryterium sukcesu:
  r(log-log) > 0.6   sensowny fit
  r(log-log) > 0.8   silny fit -> pewne predykcje dla kandydatow

Wyjscie: ps5_results.txt + predykcje T_c dla Yb, Ce, Y, FeTe, CaB6 (pod cisnieniem)
"""

import numpy as np
from scipy.optimize import minimize

# ---- stale z P4/P5 ----
PHI = (1.0 + np.sqrt(5.0)) / 2.0
G0_E = 0.869470
A_E_P4 = 0.124587    # z P4 (lepton e-soliton, ps1)
A_MU_P4 = 0.472198   # z P4 (mu-soliton, g_0^mu = phi * g_0^e)
A_TAU_P4 = 0.956027  # z P4 (tau-soliton)

C_0 = 48.8222        # J* = C_0 * A^2 (ps2)
A_STAR = 7.725       # pref. stala sieci (ps1/ps2)

A_BOHR = 0.52917721067  # Bohr w Angstromach
K_B = 8.617333e-5       # Boltzmann w eV/K


# ==============================================================================
# Baza 19 znanych SC + 7 kandydatow non-SC (dla predykcji)
# ==============================================================================
# Pola: (symbol, T_c[K], a[A], orbital, z, komentarz)
# z = liczba koordynacyjna (SC=6, BCC=8, FCC/HCP=12, A15=14, clathrate=20)
# T_c = 0 dla non-SC przy standardowym cisnieniu

KNOWN_SC = [
    ('Al',    1.175, 4.046, 'sp', 12, 'BCS klasyczny'),
    ('Zn',    0.850, 2.665, 'sp', 12, 'HCP'),
    ('Sn',    3.722, 5.831, 'sp',  8, 'white Sn, tetragon'),
    ('In',    3.408, 4.599, 'sp', 12, 'FCC'),
    ('Hg',    4.154, 2.992, 'sp',  6, 'pierwszy SC'),
    ('Pb',    7.196, 4.950, 'sp', 12, 'FCC, silnie sprz.'),
    ('Nb',    9.26,  3.301, 'd',   8, 'BCC, pierwiastek rekord'),
    ('V',     5.30,  3.027, 'd',   8, 'BCC'),
    ('Ta',    4.47,  3.303, 'd',   8, 'BCC'),
    ('Tc',    7.80,  2.739, 'd',  12, 'HCP radioaktywny'),
    ('NbN',  16.1,   4.388, 'd',  12, 'NaCl binarny'),
    ('V3Si', 17.1,   4.725, 'd',  14, 'A15'),
    ('Nb3Sn',18.3,   5.290, 'd',  14, 'A15'),
    ('Nb3Ge',23.2,   5.156, 'd',  14, 'A15'),
    ('MgB2', 39.0,   3.086, 'sp', 12, 'dwu-przerwowy'),
    ('YBCO', 92.0,   3.820, 'd',   6, 'cuprate, 2D-like'),
    ('BiSCCO',110.0, 3.815, 'd',   6, 'Bi2Sr2CaCu2O8'),
    ('H3S',  203.0,  3.089, 'sp',  8, '155 GPa'),
    ('LaH10',250.0,  3.320, 'df', 20, 'clathrate 170 GPa'),
]

CANDIDATES_NONSC = [
    ('Yb',    3.883, 'f',  12, 'f-electron, ambient para'),
    ('Ce',    3.650, 'f',  12, 'f-electron'),
    ('Y',     3.647, 'd',  12, 'SC pod cisnieniem ~12 GPa T_c=2K'),
    ('FeTe',  3.813, 'd',   8, 'parent compound'),
    ('CaB6',  4.149, 'sp',  8, 'borek wapnia'),
    ('BeB2',  3.058, 'sp', 12, ''),
    ('Ba',    5.018, 's',   8, 'SC pod 5.5 GPa T_c=5K'),
]

# ==============================================================================
# Helpers
# ==============================================================================

OUT = []
def P(s=''):
    OUT.append(str(s)); print(s)


def nearest_tgp_harmonic(a_A):
    """Zwraca (n_best, a_star_n, Delta_a = a/a_0 - a_star_n) w jedn. substr."""
    ratio = a_A / A_BOHR
    best_n, best_resid, best_a = 0, 1e9, A_STAR
    for n in range(0, 8):
        a_star_n = A_STAR + 2 * np.pi * n
        resid = abs(ratio - a_star_n)
        if resid < best_resid:
            best_resid = resid
            best_n = n
            best_a = a_star_n
    return best_n, best_a, ratio - best_a, ratio


def k_d_from_z(z):
    """3D XY Monte Carlo T_c/J, interpolacja liniowa po z."""
    # Referencje MC: SC(z=6)=2.202, BCC(z=8)=2.936, FCC(z=12)=4.403
    # Liniowo: k_d(z) = 2.20 * z / 6 (przybliz.)
    return 2.20168 * (z / 6.0)


def compute_Tc_Kelvin(a_A, orb, z, Lambda_E_meV, A_orb_map, sigma_a):
    """T_c pred = k_d(z) * C_0 * A(orb)^2 * suppress(a) * Lambda_E / k_B"""
    _, _, delta_a, _ = nearest_tgp_harmonic(a_A)
    suppress = np.exp(-delta_a**2 / sigma_a**2)
    A = A_orb_map[orb]
    J_star = C_0 * A**2
    Tc_substr = k_d_from_z(z) * J_star * suppress
    # Lambda_E_meV * 1e-3 eV / K_B eV/K = T [K]
    Tc_K = Tc_substr * (Lambda_E_meV * 1e-3) / K_B
    return Tc_K


def build_A_map(variant, params):
    """Zwraca slownik {orb -> A_orb} z parametrami fitowanymi."""
    if variant == '5a':
        Lambda_E, sigma_a = params
        A = {
            's':  A_E_P4,
            'sp': A_E_P4,            # sp = s-dominated, uzywa e-soliton
            'd':  A_MU_P4,
            'f':  A_TAU_P4,
            'df': np.sqrt(A_MU_P4 * A_TAU_P4),  # geo. mean
        }
    elif variant == '5b':
        Lambda_E, A_sp, sigma_a = params
        A = {
            's':  A_E_P4,
            'sp': A_sp,              # swobodne
            'd':  A_MU_P4,
            'f':  A_TAU_P4,
            'df': np.sqrt(A_MU_P4 * A_TAU_P4),
        }
    elif variant == '5c':
        Lambda_E, A_s, A_sp, A_d, A_f, sigma_a = params
        A = {
            's':  A_s,
            'sp': A_sp,
            'd':  A_d,
            'f':  A_f,
            'df': np.sqrt(A_d * A_f),
        }
    else:
        raise ValueError(f'unknown variant {variant}')
    return A, Lambda_E, sigma_a


def loss(params, variant, materials):
    """Suma (log10 T_c_pred - log10 T_c_obs)^2."""
    A_map, Lambda_E_meV, sigma_a = build_A_map(variant, params)
    if Lambda_E_meV <= 0 or sigma_a <= 0:
        return 1e10
    total = 0.0
    for mat in materials:
        sym, Tc_obs, a_A, orb, z, _ = mat
        if Tc_obs <= 0:
            continue
        Tc_p = compute_Tc_Kelvin(a_A, orb, z, Lambda_E_meV, A_map, sigma_a)
        if Tc_p <= 1e-6:
            return 1e10
        total += (np.log10(Tc_p) - np.log10(Tc_obs)) ** 2
    return total


def fit_variant(variant, x0, bounds=None):
    """Przeprowadz minimalizacje i zwroc (params, loss_final, status)."""
    res = minimize(
        loss, x0, args=(variant, KNOWN_SC),
        method='Nelder-Mead',
        options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 50000, 'maxfev': 50000},
    )
    return res.x, res.fun, res.success


def stats_logspace(materials, params, variant):
    """Oblicz Pearson r, RMS, MAE w log10-space."""
    A_map, Lambda_E, sigma_a = build_A_map(variant, params)
    xs, ys = [], []
    for mat in materials:
        sym, Tc_obs, a_A, orb, z, _ = mat
        if Tc_obs <= 0:
            continue
        Tc_p = compute_Tc_Kelvin(a_A, orb, z, Lambda_E, A_map, sigma_a)
        if Tc_p <= 1e-6:
            continue
        xs.append(np.log10(Tc_obs))
        ys.append(np.log10(Tc_p))
    xs = np.array(xs)
    ys = np.array(ys)
    n = len(xs)
    r = np.corrcoef(xs, ys)[0, 1]
    rms = np.sqrt(np.mean((xs - ys) ** 2))
    mae = np.mean(np.abs(xs - ys))
    return r, rms, mae, n


# ==============================================================================
# MAIN
# ==============================================================================

P("=" * 78)
P("  ps5_global_fit.py")
P("=" * 78)
P()
P("  Program P5 #5: globalny fit T_c do 19 znanych SC")
P("  Warianty: 5a (pure-TGP), 5b (sp-free), 5c (all-free)")
P()

# -----------------------------------------------------------------------------
# Part A. Baseline: ps3 one-param fit (Lambda_E calibrated na Al)
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part A.  Baseline (ps3): Lambda_E tylko z Al, A_orb = A_e (bez orbitali)")
P("=" * 78)
P()

# Skalibruj Lambda_E jak w ps3: uzywaja sie e-soliton-A dla wszystkich
# T_c(Al) = 1.175 K; k_d(12) = 4.40; A_e^2 * C_0 = 0.758
# Tc_substr = 4.40 * 0.758 = 3.335 -> Lambda_E = 1.175 * k_B / 3.335
k_d_Al = k_d_from_z(12)
Tc_substr_Al = k_d_Al * C_0 * A_E_P4**2
Lambda_E_baseline_meV = 1.175 * K_B / Tc_substr_Al * 1e3  # -> meV
P(f"  Lambda_E (baseline, Al FCC z=12, A=A_e):  {Lambda_E_baseline_meV:.4f} meV")
P()

# Baseline statystyki: orbital-naiwny model (A=A_e dla wszystkich)
def baseline_Tc(a_A, z):
    _, _, delta_a, _ = nearest_tgp_harmonic(a_A)
    suppress = np.exp(-delta_a**2 / (np.pi/2)**2)
    Tc_substr = k_d_from_z(z) * C_0 * A_E_P4**2 * suppress
    return Tc_substr * (Lambda_E_baseline_meV * 1e-3) / K_B

xs_bl, ys_bl = [], []
for mat in KNOWN_SC:
    sym, Tc_obs, a_A, orb, z, _ = mat
    Tc_p = baseline_Tc(a_A, z)
    if Tc_obs > 0 and Tc_p > 1e-6:
        xs_bl.append(np.log10(Tc_obs))
        ys_bl.append(np.log10(Tc_p))

xs_bl = np.array(xs_bl)
ys_bl = np.array(ys_bl)
r_bl = np.corrcoef(xs_bl, ys_bl)[0, 1]
rms_bl = np.sqrt(np.mean((xs_bl - ys_bl) ** 2))
mae_bl = np.mean(np.abs(xs_bl - ys_bl))
P(f"  Statystyki baseline (ps3 styl):")
P(f"    r(log-log) = {r_bl:.4f}")
P(f"    RMS_log    = {rms_bl:.4f}  (czynnik 10^RMS = {10**rms_bl:.2f}x)")
P(f"    MAE_log    = {mae_bl:.4f}")
P()

# -----------------------------------------------------------------------------
# Part B. Wariant 5a: pure-TGP, fix A z P4, fit {Lambda_E, sigma_a}
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part B.  Wariant 5a (PURE-TGP): A_s=A_e, A_d=A_mu, A_f=A_tau fix z P4")
P("=" * 78)
P()

x0 = [0.06, np.pi/2]  # Lambda_E [meV], sigma_a
params_5a, loss_5a, _ = fit_variant('5a', x0)
r_5a, rms_5a, mae_5a, n_5a = stats_logspace(KNOWN_SC, params_5a, '5a')

P(f"  Parametry fit:")
P(f"    Lambda_E = {params_5a[0]:.4f} meV  ({params_5a[0]*1e3:.2f} ueV)")
P(f"    sigma_a  = {params_5a[1]:.4f}     ({params_5a[1]/np.pi:.3f} * pi)")
P()
P(f"  Fit ilosciowy:")
P(f"    r(log-log) = {r_5a:.4f}  (baseline: {r_bl:.4f})")
P(f"    RMS_log    = {rms_5a:.4f}  (baseline: {rms_bl:.4f})")
P(f"    MAE_log    = {mae_5a:.4f}  (baseline: {mae_bl:.4f})")
P(f"    loss       = {loss_5a:.4f}")
P()

# -----------------------------------------------------------------------------
# Part C. Wariant 5b: fit {Lambda_E, A_sp, sigma_a}
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part C.  Wariant 5b (SP-FREE): A_sp swobodne, reszta fix z P4")
P("=" * 78)
P()

x0 = [0.06, A_E_P4, np.pi/2]
params_5b, loss_5b, _ = fit_variant('5b', x0)
r_5b, rms_5b, mae_5b, n_5b = stats_logspace(KNOWN_SC, params_5b, '5b')

P(f"  Parametry fit:")
P(f"    Lambda_E = {params_5b[0]:.4f} meV")
P(f"    A_sp     = {params_5b[1]:.4f}   (A_e_P4 = {A_E_P4:.4f})")
P(f"              stosunek A_sp / A_e = {params_5b[1]/A_E_P4:.3f}")
P(f"    sigma_a  = {params_5b[2]:.4f}")
P()
P(f"  Fit ilosciowy:")
P(f"    r(log-log) = {r_5b:.4f}")
P(f"    RMS_log    = {rms_5b:.4f}")
P(f"    MAE_log    = {mae_5b:.4f}")
P(f"    loss       = {loss_5b:.4f}")
P()

# -----------------------------------------------------------------------------
# Part D. Wariant 5c: wszystkie amplitudy swobodne
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part D.  Wariant 5c (ALL-FREE): {Lambda_E, A_s, A_sp, A_d, A_f, sigma_a}")
P("=" * 78)
P()

x0 = [0.06, A_E_P4, A_E_P4, A_MU_P4, A_TAU_P4, np.pi/2]
params_5c, loss_5c, _ = fit_variant('5c', x0)
r_5c, rms_5c, mae_5c, n_5c = stats_logspace(KNOWN_SC, params_5c, '5c')

Lambda_E_c, A_s_c, A_sp_c, A_d_c, A_f_c, sigma_c = params_5c

P(f"  Parametry fit:")
P(f"    Lambda_E = {Lambda_E_c:.4f} meV")
P(f"    A_s      = {A_s_c:.4f}   (P4: A_e   = {A_E_P4:.4f},  stosunek {A_s_c/A_E_P4:.3f})")
P(f"    A_sp     = {A_sp_c:.4f}   (ref A_e  = {A_E_P4:.4f},  stosunek {A_sp_c/A_E_P4:.3f})")
P(f"    A_d      = {A_d_c:.4f}   (P4: A_mu  = {A_MU_P4:.4f},  stosunek {A_d_c/A_MU_P4:.3f})")
P(f"    A_f      = {A_f_c:.4f}   (P4: A_tau = {A_TAU_P4:.4f},  stosunek {A_f_c/A_TAU_P4:.3f})")
P(f"    sigma_a  = {sigma_c:.4f}")
P()
P(f"  Fit ilosciowy:")
P(f"    r(log-log) = {r_5c:.4f}")
P(f"    RMS_log    = {rms_5c:.4f}")
P(f"    MAE_log    = {mae_5c:.4f}")
P(f"    loss       = {loss_5c:.4f}")
P()

# -----------------------------------------------------------------------------
# Part E. Podsumowanie porownawcze
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part E.  Podsumowanie 4 wariantow")
P("=" * 78)
P()
P(f"  {'Wariant':12s} {'params':>6s}  {'r':>7s}  {'RMS_log':>8s}  {'loss':>9s}")
P(f"  {'-'*12:12s} {'-'*6:>6s}  {'-'*7:>7s}  {'-'*8:>8s}  {'-'*9:>9s}")
P(f"  {'Baseline':12s} {'1':>6s}  {r_bl:7.4f}  {rms_bl:8.4f}  {'-':>9s}")
P(f"  {'5a (pure)':12s} {'2':>6s}  {r_5a:7.4f}  {rms_5a:8.4f}  {loss_5a:9.4f}")
P(f"  {'5b (sp-free)':12s} {'3':>6s}  {r_5b:7.4f}  {rms_5b:8.4f}  {loss_5b:9.4f}")
P(f"  {'5c (all-free)':12s} {'6':>6s}  {r_5c:7.4f}  {rms_5c:8.4f}  {loss_5c:9.4f}")
P()

# -----------------------------------------------------------------------------
# Part F. Tabela predykcji dla 19 SC (3 warianty)
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part F.  Predykcja T_c dla 19 znanych SC - wszystkie warianty")
P("=" * 78)
P()
P(f"  {'SC':>8s}  {'obs':>7s}  {'5a':>8s}  {'5b':>8s}  {'5c':>8s}  {'r5c':>6s}")
P(f"  {'-'*8:>8s}  {'-'*7:>7s}  {'-'*8:>8s}  {'-'*8:>8s}  {'-'*8:>8s}  {'-'*6:>6s}")

A_map_a, Lam_a, sig_a = build_A_map('5a', params_5a)
A_map_b, Lam_b, sig_b = build_A_map('5b', params_5b)
A_map_c, Lam_c, sig_c = build_A_map('5c', params_5c)

for mat in KNOWN_SC:
    sym, Tc_obs, a_A, orb, z, _ = mat
    Tc_a = compute_Tc_Kelvin(a_A, orb, z, Lam_a, A_map_a, sig_a)
    Tc_b = compute_Tc_Kelvin(a_A, orb, z, Lam_b, A_map_b, sig_b)
    Tc_c = compute_Tc_Kelvin(a_A, orb, z, Lam_c, A_map_c, sig_c)
    ratio_c = Tc_c / Tc_obs if Tc_obs > 0 else 0
    P(f"  {sym:>8s}  {Tc_obs:7.2f}  {Tc_a:8.2f}  {Tc_b:8.2f}  {Tc_c:8.2f}  {ratio_c:6.2f}")
P()

# -----------------------------------------------------------------------------
# Part G. Predykcja T_c dla kandydatow non-SC (z 5c)
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part G.  Predykcja T_c dla kandydatow non-SC (model 5c)")
P("=" * 78)
P()
P(f"  Uwaga: Tc_pred dotyczy sytuacji gdy material jest fazowo-koherentny.")
P(f"        Dla non-SC przy ambient presumujemy: inny mechanizm blokuje faze")
P(f"        (magnetyzm, lokalizacja, Mott). Pod cisnieniem moze sie odblokowac.")
P()
P(f"  {'Mat':>6s}  {'a[A]':>6s}  {'orb':>4s}  {'z':>3s}  {'Tc_5c[K]':>9s}  {'Komentarz':<35s}")
P(f"  {'-'*6:>6s}  {'-'*6:>6s}  {'-'*4:>4s}  {'-'*3:>3s}  {'-'*9:>9s}  {'-'*35:<35s}")
for cand in CANDIDATES_NONSC:
    sym, a_A, orb, z, com = cand
    Tc_c = compute_Tc_Kelvin(a_A, orb, z, Lam_c, A_map_c, sig_c)
    P(f"  {sym:>6s}  {a_A:6.3f}  {orb:>4s}  {z:3d}  {Tc_c:9.3f}  {com:<35s}")
P()

# -----------------------------------------------------------------------------
# Part H. Reszty per-materia (5c) - gdzie model pada?
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part H.  Reszty per-material (5c)")
P("=" * 78)
P()
P(f"  {'SC':>8s}  {'obs':>7s}  {'5c':>8s}  {'log10(p/o)':>11s}  {'komentarz':<25s}")
P(f"  {'-'*8:>8s}  {'-'*7:>7s}  {'-'*8:>8s}  {'-'*11:>11s}  {'-'*25:<25s}")
for mat in KNOWN_SC:
    sym, Tc_obs, a_A, orb, z, com = mat
    Tc_c = compute_Tc_Kelvin(a_A, orb, z, Lam_c, A_map_c, sig_c)
    if Tc_obs > 0 and Tc_c > 1e-6:
        r_log = np.log10(Tc_c / Tc_obs)
    else:
        r_log = float('nan')
    P(f"  {sym:>8s}  {Tc_obs:7.2f}  {Tc_c:8.2f}  {r_log:+11.4f}  {com:<25s}")
P()

# -----------------------------------------------------------------------------
# Part I. Werdykt
# -----------------------------------------------------------------------------
P("=" * 78)
P("  Part I.  Werdykt ps5")
P("=" * 78)
P()

P(f"  Poprawa r(log-log) od baseline: {r_bl:.3f}")
for name, r in [('5a', r_5a), ('5b', r_5b), ('5c', r_5c)]:
    delta = r - r_bl
    P(f"    {name}: r = {r:.4f}  (poprawa Delta = {delta:+.4f})")
P()

best_variant = max([('5a', r_5a, loss_5a), ('5b', r_5b, loss_5b), ('5c', r_5c, loss_5c)],
                    key=lambda t: t[1])
P(f"  Najwyzsze r: wariant {best_variant[0]} (r = {best_variant[1]:.4f})")
P()

if best_variant[1] > 0.8:
    status = 'SILNY fit -> pewne predykcje dla kandydatow'
elif best_variant[1] > 0.6:
    status = 'SENSOWNY fit -> predykcje kwalitatywne OK'
elif best_variant[1] > 0.4:
    status = 'SLABY fit -> potrzebne dalsze poprawki (ps6+)'
else:
    status = 'NIEWYSTARCZAJACY fit -> baza niedostatecznie opisowa'
P(f"  Status ilosciowy: {status}")
P()

# Sprawdzenie P4-P5 mapowania (5c vs P4)
rat_Ad = A_d_c / A_MU_P4
rat_Af = A_f_c / A_TAU_P4
P(f"  Test P4-P5 mapowania (wariant 5c vs particle sector P4):")
P(f"    A_d (fit) / A_mu (P4) = {rat_Ad:.3f}")
P(f"    A_f (fit) / A_tau (P4) = {rat_Af:.3f}")
P()

if 0.7 < rat_Ad < 1.4 and 0.7 < rat_Af < 1.4:
    p45_status = 'P4-P5 mapowanie POTWIERDZONE numerycznie'
elif 0.4 < rat_Ad < 2.5 and 0.4 < rat_Af < 2.5:
    p45_status = 'P4-P5 mapowanie JAKOSCIOWE (rzad wielkosci OK)'
else:
    p45_status = 'P4-P5 mapowanie ZAKWESTIONOWANE - solitony leptonowe != orbital amplitudy'
P(f"  {p45_status}")
P()

P(f"  STATUS ps5: Hipoteza -> {'HIPOTEZA (silna)' if best_variant[1] > 0.6 else 'HIPOTEZA (slaba)'}")
P(f"              Formula T_c pasuje z r = {best_variant[1]:.3f} na 19 SC.")
P()

# -----------------------------------------------------------------------------
# Zapis
# -----------------------------------------------------------------------------
P("=" * 78)
P("  ps5 complete.")
P("=" * 78)

with open('ps5_results.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(OUT) + '\n')
