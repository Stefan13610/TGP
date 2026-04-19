"""
ps17_full_p6_validation.py  -  P6 COMPLETE validation

Cel:
  Spina wszystkie mechanizmy P6 (A+B+C+D) w jednej formule
  i testuje na zjednoczonym datasecie 25+ materialow.

Architektura:
  - Cuprates: P6.A (d-wave + Zhang-Rice + van Hove + layers)
  - Wszystko inne: P6.B (phonon) * P6.D (blocking) * A_eff(eta) z P6.C

Walidacja:
  1. Per-class r, RMS_log
  2. Globalny r na wszystkich 25+ materialach
  3. Lista outlierow (|dlog| > 0.4) -> kandydaci do P7
  4. Prediction vs observed scatter plot (ASCII)
"""

import numpy as np

# =============================================================
# Stale (z poprzednich ps)
# =============================================================

K_B = 8.617333e-5
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088
C_0 = 48.8222
sigma_a = 2.5856

A_map = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}

# P6.A (cuprates, z ps10)
K_dw = 3.498
Lambda_E_cup = 0.0513  # meV
A_ZR_sq = 0.181

# P6.B (phonon, z ps12)
alpha_P6B = 1.04
Lambda_0_P6B = 0.0962   # meV
omega_0 = 15.0          # meV

# P6.D (blocking, z ps15)
beta_P6D = 2.527

# P6.C (orbital switching, z ps16)
# eta dla Ce: P_scale = 5.8 GPa

def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma_a**2)


def Tc_cuprate(a_A, n_layers, z_planar=8):
    """P6.A formula."""
    layer_factor = n_layers ** 0.5
    M = M_gauss(a_A)
    Tc_substr = K_dw * k_d(z_planar) * C_0 * A_ZR_sq * M * layer_factor
    return Tc_substr * (Lambda_E_cup * 1e-3) / K_B


def Tc_phonon(a_A, orb_or_eta, z, omega_phonon, lambda_sf):
    """P6.B * P6.D * P6.C uniwersalna formula.

    orb_or_eta:
      - string 's', 'sp', 'd', 'f' -> odpowiednie A_map
      - float w [0,1] -> A_eff = eta * A_d (P6.C)
    """
    if isinstance(orb_or_eta, str):
        A = A_map[orb_or_eta]
    else:
        A = orb_or_eta * A_map["d"]

    J = C_0 * A**2
    M = M_gauss(a_A)
    boost = (omega_phonon / omega_0) ** alpha_P6B
    Lambda_eff = Lambda_0_P6B * boost
    B_mag = 1.0 / (1.0 + beta_P6D * lambda_sf)

    Tc_substr = k_d(z) * J * M
    return Tc_substr * (Lambda_eff * 1e-3) / K_B * B_mag


# =============================================================
# ZJEDNOCZONY DATASET
# =============================================================

# (nazwa, klasa, parametry, T_obs)
# klasa: 'cuprate' (uzywa Tc_cuprate), 'phonon' (uzywa Tc_phonon)

# Cuprates (P6.A, 8 material)
cuprates = [
    # (nazwa, a, n_layers, T_obs)
    ("La2CuO4",    3.780, 1, 38.0),
    ("YBCO",       3.820, 2, 92.0),
    ("BiSCCO2212", 3.820, 2, 85.0),
    ("Tl2212",     3.850, 2, 108.0),
    ("Hg1223",     3.855, 3, 138.0),
    ("Tl2223",     3.820, 3, 125.0),
    ("Nd2CuO4",    3.945, 1, 24.0),
    ("Bi2201",     3.810, 1, 34.0),
]

# Phonon-mediated + other (P6.B+C+D, ~17 material)
phonon = [
    # (nazwa, a, orb_or_eta, z, omega, lambda_sf, T_obs)
    # --- BCS klasyczne ---
    ("Al",          4.046, "s",  12,  15.0, 0.00,  1.18),
    ("Pb",          4.950, "sp", 12,   8.3, 0.00,  7.20),
    ("Nb",          3.301, "d",   8,  22.0, 0.20,  9.26),
    ("V",           3.024, "d",   8,  31.0, 0.60,  5.30),
    ("Hg",          2.989, "sp",  6,   9.0, 0.00,  4.15),

    # --- MgB2 ---
    ("MgB2",        3.086, "sp",  5,  75.0, 0.00, 39.00),

    # --- FeSe ---
    ("FeSe_bulk",   3.770, "d",   8,  25.0, 0.90,  8.00),
    ("FeSe/STO",    3.770, "d",   8,  80.0, 0.20, 65.00),

    # --- Pnictidy Fe ---
    ("Ba122-Co",    3.960, "d",   8,  30.0, 0.30, 22.00),
    ("LaFeAsO",     4.040, "d",   8,  35.0, 0.50, 26.00),
    ("NdFeAsO-F",   3.970, "d",   8,  40.0, 0.30, 55.00),

    # --- Intermetalliki ---
    ("Nb3Sn",       5.290, "d",  12,  22.0, 0.10, 18.30),
    ("NbTi",        3.300, "d",   8,  25.0, 0.15, 10.00),

    # --- Hydrydy high-P (z P6.C eta=1) ---
    ("H3S",         3.100, "sp",  8, 175.0, 0.00, 203.0),
    ("LaH10",       5.100, 1.0,  12, 250.0, 0.00, 250.0),
    ("CeH9",        3.500, 1.0,   8, 135.0, 0.00, 100.0),
    ("CeH10",       3.500, 1.0,   8, 175.0, 0.00, 115.0),

    # --- f-metale ambient (z P6.C eta<1) ---
    ("La_amb",      3.770, 1.0,  12,  12.0, 0.10,  6.00),
    ("Y_amb",       3.648, 1.0,  12,  17.0, 0.15,  1.30),
    ("Th_amb",      5.080, 1.0,  12,  10.0, 0.10,  1.38),
    ("Ce_amb",      5.160, 0.5,  12,  11.0, 0.20,  0.01),
    ("Yb_amb",      5.485, 0.0,  12,  16.0, 0.10,  0.01),
    ("Ce_5GPa",     4.900, 0.8,  12,  12.0, 0.30,  1.70),
]


# =============================================================
# Walidacja per-class
# =============================================================

print("=" * 78)
print("  ps17_full_p6_validation.py  -  PELNE P6 (A+B+C+D)")
print("=" * 78)
print()

# --- Cuprates ---
print("=" * 78)
print("  CLASS 1: Cuprates (P6.A)")
print("=" * 78)
print()
print(f"  {'Material':>12} {'a':>5} {'n':>2} {'T_obs':>7} {'T_pred':>7} {'dlog':>7}")
log_obs_cup, log_pred_cup = [], []
outliers_cup = []
for name, a, n, Tobs in cuprates:
    Tp = Tc_cuprate(a, n)
    lo_ = np.log10(Tobs)
    lp = np.log10(max(Tp, 1e-6))
    log_obs_cup.append(lo_)
    log_pred_cup.append(lp)
    dl = lp - lo_
    mark = " *" if abs(dl) > 0.4 else ""
    if abs(dl) > 0.4:
        outliers_cup.append((name, dl))
    print(f"  {name:>12} {a:>5.3f} {n:>2d} {Tobs:>7.2f} {Tp:>7.2f} {dl:>+6.3f}{mark}")
print()

log_obs_cup = np.array(log_obs_cup)
log_pred_cup = np.array(log_pred_cup)
r_cup = np.corrcoef(log_obs_cup, log_pred_cup)[0, 1]
rms_cup = np.sqrt(np.mean((log_pred_cup - log_obs_cup)**2))
print(f"  Cuprates (N=8): r={r_cup:.4f}, RMS_log={rms_cup:.4f}")
print()

# --- Phonon-mediated + others ---
print("=" * 78)
print("  CLASS 2: Phonon-mediated + others (P6.B + P6.C + P6.D)")
print("=" * 78)
print()
print(f"  {'Material':>12} {'a':>5} {'orb/eta':>7} {'z':>2} {'omega':>5} {'lam_sf':>6} "
      f"{'T_obs':>7} {'T_pred':>7} {'dlog':>7}")

log_obs_ph, log_pred_ph = [], []
outliers_ph = []
for name, a, oe, z, om, lam, Tobs in phonon:
    Tp = Tc_phonon(a, oe, z, om, lam)
    if Tobs < 0.05:
        continue  # skip non-SC dla log analizy
    lo_ = np.log10(Tobs)
    lp = np.log10(max(Tp, 1e-6))
    log_obs_ph.append(lo_)
    log_pred_ph.append(lp)
    dl = lp - lo_
    mark = " *" if abs(dl) > 0.4 else ""
    if abs(dl) > 0.4:
        outliers_ph.append((name, dl))
    oe_str = str(oe) if isinstance(oe, str) else f"{oe:.2f}"
    print(f"  {name:>12} {a:>5.3f} {oe_str:>7} {z:>2d} {om:>5.0f} {lam:>6.2f} "
          f"{Tobs:>7.2f} {Tp:>7.2f} {dl:>+6.3f}{mark}")
print()

log_obs_ph = np.array(log_obs_ph)
log_pred_ph = np.array(log_pred_ph)
r_ph = np.corrcoef(log_obs_ph, log_pred_ph)[0, 1]
rms_ph = np.sqrt(np.mean((log_pred_ph - log_obs_ph)**2))
print(f"  Phonon-others (N={len(log_obs_ph)}): r={r_ph:.4f}, RMS_log={rms_ph:.4f}")
print()

# --- Globalny ---
print("=" * 78)
print("  GLOBAL r, RMS_log (cuprates + phonon-others)")
print("=" * 78)
print()
log_obs_all = np.concatenate([log_obs_cup, log_obs_ph])
log_pred_all = np.concatenate([log_pred_cup, log_pred_ph])
r_all = np.corrcoef(log_obs_all, log_pred_all)[0, 1]
rms_all = np.sqrt(np.mean((log_pred_all - log_obs_all)**2))
print(f"  N_total = {len(log_obs_all)}")
print(f"  r(log)   = {r_all:.4f}")
print(f"  RMS_log  = {rms_all:.4f}")
print()

# --- Outliers lista ---
print("=" * 78)
print("  Outliers (|dlog| > 0.4) -> kandydaci P7")
print("=" * 78)
print()
all_outliers = [(n, dl, "cup") for n, dl in outliers_cup] + [(n, dl, "phon") for n, dl in outliers_ph]
if all_outliers:
    all_outliers.sort(key=lambda x: -abs(x[1]))
    for name, dl, cls in all_outliers:
        sign = "OVER" if dl > 0 else "UNDER"
        print(f"  {name:>14} ({cls:>4}) dlog = {dl:+.3f}  ({sign})")
else:
    print("  Brak outlierow - pelny P6 domknal zbior.")
print()

# --- Porownanie z P5 i etapami P6 ---
print("=" * 78)
print("  Progres domkniecia: P5 -> P6.B -> P6.B+D -> pelne P6")
print("=" * 78)
print()
print(f"  {'Etap':>20} {'N':>3} {'r(log)':>8} {'RMS_log':>8}")
print(f"  {'P5 ps5 5c':>20} {'19':>3} {'0.480':>8} {'0.62':>8}")
print(f"  {'P6.A cuprates':>20} {'8':>3} {'0.957':>8} {'0.19':>8}")
print(f"  {'P6.B uniwersalny':>20} {'16':>3} {'0.877':>8} {'0.32':>8}")
print(f"  {'P6.B + P6.D':>20} {'16':>3} {'0.930':>8} {'0.23':>8}")
print(f"  {'Pelne P6 (A+B+C+D)':>20} {len(log_obs_all):>3d} {r_all:>8.3f} {rms_all:>8.3f}")
print()

print("=" * 78)
print("  ps17 complete.")
print("=" * 78)
