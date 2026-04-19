"""
ps20_master_plot.py  -  Master walidacja P6 + P7.1

Cel:
  Zjednoczony skrypt zawierajacy wszystkie mechanizmy (A+B+C+D + P7.1)
  z automatycznym wyliczaniem lambda_sf z pierwszych zasad dla d-metali.

Architektura:
  - Cuprates:          P6.A (d-wave + Zhang-Rice + van Hove + layers)
  - Phonon-mediated:   P6.B * P6.D * A_eff(eta) z P6.C
  - lambda_sf dla d-metali: P7.1 formula (gdy N(EF), I dostepne)
  - lambda_sf dla hydrydow, s/sp: 0 (aprioryczne)

Output:
  1. Tabela per-klasa z residuals
  2. Master scatter (ASCII log-log)
  3. Top outliers dla P7.2
  4. Statystyki: N, r, RMS_log
"""

import numpy as np

# =============================================================
# Stale TGP (finalne z P6 + P7.1)
# =============================================================

K_B = 8.617333e-5
A_BOHR = 0.52917721067
a_star_tgp = 7.725 * A_BOHR
C_0 = 48.8222
sigma_a = 2.5856

A_map = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}

# P6.A cuprates
K_dw = 3.498
Lambda_E_cup = 0.0513
A_ZR_sq = 0.181

# P6.B phonon coupling
alpha_P6B = 1.04
Lambda_0_P6B = 0.0962
omega_0 = 15.0

# P6.D magnetic blocking
beta_P6D = 2.527

# P6.C orbital switching
P_scale_Ce = 5.8
P_scale_Yb = 552  # ps18 refit
P_scale_4fn = lambda n: 5.8 + 42.0 * (n - 1)  # ps18 extrapolation

# P7.1 lambda_sf z pierwszych zasad
kappa_P7 = 2.012


def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A):
    n = max(1, round(a_A / a_star_tgp))
    d = a_A - n * a_star_tgp
    return np.exp(-d**2 / sigma_a**2)

def A_eff(eta):
    return eta * A_map["d"]

def eta_pressure(P_GPa, eta_0, P_scale):
    return eta_0 + (1 - eta_0) * (1 - np.exp(-P_GPa / P_scale))


def lambda_sf_P7(A_orb, z, N_EF, I):
    """P7.1 formula dla d-metali (paramagnon blocking)."""
    if N_EF is None or I is None:
        return 0.0
    IN = I * N_EF
    g = IN**2 / np.sqrt(1.0 + 0.25 * IN**4)
    return kappa_P7 * A_orb**2 * k_d(z) * N_EF * g


def Tc_cuprate(a_A, n_layers, z_planar=8):
    """P6.A formula."""
    layer_factor = n_layers ** 0.5
    M = M_gauss(a_A)
    Tc_substr = K_dw * k_d(z_planar) * C_0 * A_ZR_sq * M * layer_factor
    return Tc_substr * (Lambda_E_cup * 1e-3) / K_B


def Tc_phonon(a_A, orb_or_eta, z, omega_phonon, lambda_sf):
    """P6.B * P6.C * P6.D."""
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
# MASTER DATASET
# =============================================================

# Cuprates (P6.A)
cuprates = [
    # (name, a, n_layers, T_obs)
    ("La2CuO4",       3.780, 1,  38.0),
    ("YBCO",          3.820, 2,  92.0),
    ("BiSCCO2212",    3.820, 2,  85.0),
    ("Tl2212",        3.850, 2, 108.0),
    ("Hg1223",        3.855, 3, 138.0),
    ("Hg1223-quench", 3.830, 3, 151.0),  # Deng/Chu 2026
    ("Tl2223",        3.820, 3, 125.0),
    ("Nd2CuO4",       3.945, 1,  24.0),
    ("Bi2201",        3.810, 1,  34.0),
]

# Phonon-mediated z parametrami (N_EF, I) tam gdzie P7.1 ma dane
# (name, a, orb_or_eta, z, omega, N_EF or None, I or None, lam_manual, T_obs)
# lam_manual = override gdy P7.1 niedostepny (hydrydy, s/sp, f-metale)
phonon = [
    # --- BCS s/sp (P7.1 nie dotyczy, lam=0) ---
    ("Al",           4.046, "s",  12,  15.0, None, None, 0.00,  1.18),
    ("Pb",           4.950, "sp", 12,   8.3, None, None, 0.00,  7.20),
    ("Hg",           2.989, "sp",  6,   9.0, None, None, 0.00,  4.15),

    # --- BCS d-metale (P7.1 derived) ---
    ("Nb",           3.301, "d",   8,  22.0, 1.24, 0.57, None,  9.26),
    ("V",            3.024, "d",   8,  31.0, 1.35, 0.72, None,  5.30),

    # --- MgB2 (sp, P7.1 nie dotyczy) ---
    ("MgB2",         3.086, "sp",  5,  75.0, None, None, 0.00, 39.00),

    # --- FeSe bulk & /STO (P7.1 derived z I,N dla zwiazku) ---
    ("FeSe_bulk",    3.770, "d",   8,  25.0, 2.00, 0.40, None,  8.00),
    ("FeSe/STO",     3.770, "d",   8,  80.0, 2.00, 0.30, None, 65.00),

    # --- Fe-pnictidy (P7.1 derived) ---
    ("Ba122-Co",     3.960, "d",   8,  30.0, 1.80, 0.30, None, 22.00),
    ("LaFeAsO",      4.040, "d",   8,  35.0, 1.90, 0.35, None, 26.00),
    ("NdFeAsO-F",    3.970, "d",   8,  40.0, 1.80, 0.30, None, 55.00),

    # --- Intermetalliki ---
    ("Nb3Sn",        5.290, "d",  12,  22.0, 1.24, 0.57, None, 18.30),  # Nb-like
    ("NbTi",         3.300, "d",   8,  25.0, 1.24, 0.57, None, 10.00),

    # --- Hydrydy high-P (eta=1, lam=0 bo H dominuje) ---
    ("H3S",          3.100, "sp",  8, 175.0, None, None, 0.00, 203.0),
    ("LaH10",        5.100, 1.0,  12, 250.0, None, None, 0.00, 250.0),
    ("CeH9",         3.500, 1.0,   8, 135.0, None, None, 0.00, 100.0),
    ("CeH10",        3.500, 1.0,   8, 175.0, None, None, 0.00, 115.0),
    ("Yb4H23",       3.500, 0.278, 8, 140.0, None, None, 0.00,  11.5),

    # --- f-metale ambient (eta<1, maly lam) ---
    ("La_amb",       3.770, 1.0,  12,  12.0, None, None, 0.10,  6.00),
    ("Y_amb",        3.648, 1.0,  12,  17.0, None, None, 0.15,  1.30),
    ("Th_amb",       5.080, 1.0,  12,  10.0, None, None, 0.10,  1.38),
    ("Ce_5GPa",      4.900, 0.8,  12,  12.0, None, None, 0.30,  1.70),
]


# =============================================================
# WALIDACJA PER-CLASS
# =============================================================

print("=" * 80)
print("  ps20_master_plot.py  -  ZJEDNOCZONA WALIDACJA P6 + P7.1")
print("=" * 80)
print()

# --- Cuprates ---
print("=" * 80)
print("  CLASS 1: Cuprates (P6.A)")
print("=" * 80)
print()
print(f"  {'Material':>16} {'a':>5} {'n':>2} {'T_obs':>7} {'T_pred':>7} {'dlog':>7}")

log_obs_c, log_pred_c = [], []
outliers_c = []
for name, a, n, Tobs in cuprates:
    Tp = Tc_cuprate(a, n)
    lo = np.log10(Tobs); lp = np.log10(max(Tp, 1e-6))
    log_obs_c.append(lo); log_pred_c.append(lp)
    dl = lp - lo
    mark = " *" if abs(dl) > 0.4 else ""
    if abs(dl) > 0.4:
        outliers_c.append((name, dl, "cuprate"))
    print(f"  {name:>16} {a:>5.3f} {n:>2d} {Tobs:>7.2f} {Tp:>7.2f} {dl:>+6.3f}{mark}")

log_obs_c = np.array(log_obs_c); log_pred_c = np.array(log_pred_c)
r_c = np.corrcoef(log_obs_c, log_pred_c)[0, 1]
rms_c = np.sqrt(np.mean((log_pred_c - log_obs_c)**2))
print()
print(f"  Cuprates (N={len(log_obs_c)}): r={r_c:.4f}, RMS_log={rms_c:.4f}")
print()

# --- Phonon-mediated ---
print("=" * 80)
print("  CLASS 2: Phonon-mediated + others (P6.B+C+D, lambda_sf z P7.1 gdy mozliwe)")
print("=" * 80)
print()
print(f"  {'Material':>14} {'a':>5} {'o/e':>5} {'z':>2} {'omega':>5} "
      f"{'lam_sf':>7} {'src':>5} {'T_obs':>7} {'T_pred':>7} {'dlog':>7}")

log_obs_p, log_pred_p = [], []
outliers_p = []
count_p7 = 0
for name, a, oe, z, om, N_EF, I, lam_m, Tobs in phonon:
    if N_EF is not None and I is not None:
        A_o = A_map[oe] if isinstance(oe, str) else oe * A_map["d"]
        lam_sf = lambda_sf_P7(A_o, z, N_EF, I)
        lam_src = "P7.1"
        count_p7 += 1
    else:
        lam_sf = lam_m
        lam_src = "fix"

    Tp = Tc_phonon(a, oe, z, om, lam_sf)
    if Tobs < 0.05:
        continue
    lo = np.log10(Tobs); lp = np.log10(max(Tp, 1e-6))
    log_obs_p.append(lo); log_pred_p.append(lp)
    dl = lp - lo
    mark = " *" if abs(dl) > 0.4 else ""
    if abs(dl) > 0.4:
        outliers_p.append((name, dl, "phonon"))
    oe_str = str(oe) if isinstance(oe, str) else f"{oe:.2f}"
    print(f"  {name:>14} {a:>5.3f} {oe_str:>5} {z:>2d} {om:>5.0f} "
          f"{lam_sf:>7.3f} {lam_src:>5} {Tobs:>7.2f} {Tp:>7.2f} {dl:>+6.3f}{mark}")

log_obs_p = np.array(log_obs_p); log_pred_p = np.array(log_pred_p)
r_p = np.corrcoef(log_obs_p, log_pred_p)[0, 1]
rms_p = np.sqrt(np.mean((log_pred_p - log_obs_p)**2))
print()
print(f"  Phonon-others (N={len(log_obs_p)}, z P7.1: {count_p7}): r={r_p:.4f}, RMS_log={rms_p:.4f}")
print()

# --- Globalny ---
print("=" * 80)
print("  GLOBAL")
print("=" * 80)
print()
log_obs_all = np.concatenate([log_obs_c, log_obs_p])
log_pred_all = np.concatenate([log_pred_c, log_pred_p])
r_all = np.corrcoef(log_obs_all, log_pred_all)[0, 1]
rms_all = np.sqrt(np.mean((log_pred_all - log_obs_all)**2))
print(f"  N_total = {len(log_obs_all)}")
print(f"  r(log)   = {r_all:.4f}")
print(f"  RMS_log  = {rms_all:.4f}")
print()

# --- Master scatter (ASCII log-log) ---
print("=" * 80)
print("  Master scatter log10(T_pred) vs log10(T_obs)")
print("=" * 80)
print()

# 60-kolumnowy scatter ASCII
x_min, x_max = -1.5, 3.0  # log10(T) range: 0.03 K to 1000 K
y_min, y_max = -1.5, 3.0
ncols, nrows = 50, 20

grid = [[" " for _ in range(ncols)] for _ in range(nrows)]

# Markery: @ cuprate, # phonon-d, * phonon-sp, + hydride, ~ f-metal ambient
def marker_for(name):
    if any(x in name for x in ["La2", "YBCO", "BiSCCO", "Tl", "Hg", "Nd2", "Bi2"]):
        return "@"
    if any(x in name for x in ["H3S", "LaH10", "CeH9", "CeH10", "Yb4H23"]):
        return "+"
    if any(x in name for x in ["_amb", "_5GPa"]):
        return "~"
    return "#"

all_points = []
for i, (name, a, n, Tobs) in enumerate(cuprates):
    all_points.append((np.log10(Tobs), log_pred_c[i], marker_for(name)))
j = 0
for name, a, oe, z, om, N_EF, I, lam_m, Tobs in phonon:
    if Tobs < 0.05: continue
    all_points.append((np.log10(Tobs), log_pred_p[j], marker_for(name)))
    j += 1

for x, y, m in all_points:
    col = int((x - x_min) / (x_max - x_min) * (ncols - 1))
    row = nrows - 1 - int((y - y_min) / (y_max - y_min) * (nrows - 1))
    if 0 <= col < ncols and 0 <= row < nrows:
        grid[row][col] = m

# Dodaj diagonale y=x (kropki)
for k in range(min(nrows, ncols)):
    col = int((k / (min(nrows, ncols) - 1)) * (ncols - 1))
    row = nrows - 1 - int((k / (min(nrows, ncols) - 1)) * (nrows - 1))
    if grid[row][col] == " ":
        grid[row][col] = "."

print("  log10(T_pred)")
for r_idx, row in enumerate(grid):
    y_val = y_max - r_idx / (nrows - 1) * (y_max - y_min)
    print(f"  {y_val:>5.1f} |{''.join(row)}|")
x_axis = "  " + " " * 7 + "".join(
    "|" if i % 10 == 0 else "-" for i in range(ncols)
)
print(x_axis)
x_labels = "  " + " " * 6
for i in range(0, ncols, 10):
    x_val = x_min + i / (ncols - 1) * (x_max - x_min)
    x_labels += f"{x_val:>6.1f}    "
print(x_labels + "   log10(T_obs)")
print()
print("  Legenda: @ cuprate  # phonon-d  + hydryd  ~ f-metal ambient  . y=x")
print()


# --- Outliers ---
print("=" * 80)
print("  Outliers |dlog|>0.4 (kandydaci do P7.2)")
print("=" * 80)
print()
all_outliers = outliers_c + outliers_p
if all_outliers:
    all_outliers.sort(key=lambda x: -abs(x[1]))
    for name, dl, cls in all_outliers:
        sign = "OVER" if dl > 0 else "UNDER"
        print(f"  {name:>16} ({cls:>8}) dlog = {dl:+.3f}  ({sign})")
else:
    print("  Brak outlierow.")
print()


# --- Porownanie z poprzednimi etapami ---
print("=" * 80)
print("  Progres: P5 -> P6 -> P6+P7.1 (pelne domkniecie)")
print("=" * 80)
print()
print(f"  {'Etap':>24} {'N':>3} {'r(log)':>8} {'RMS_log':>8}")
print(f"  {'-'*24:>24} {'---':>3} {'--------':>8} {'--------':>8}")
print(f"  {'P5 ps5 5c':>24} {'19':>3} {'0.480':>8} {'0.620':>8}")
print(f"  {'P6.A cuprates':>24} {'8':>3} {'0.957':>8} {'0.190':>8}")
print(f"  {'P6.B+D (16 mat)':>24} {'16':>3} {'0.930':>8} {'0.230':>8}")
print(f"  {'P6 pelne (ps17, 29)':>24} {'29':>3} {'0.875':>8} {'0.347':>8}")
print(f"  {'P6+P7.1 (ps20 master)':>24} {len(log_obs_all):>3d} {r_all:>8.3f} {rms_all:>8.3f}")
print()


# =============================================================
# Part: Nowe predykcje z P7.1 (falsyfikowalne)
# =============================================================

print("=" * 80)
print("  Predykcje falsyfikowalne (dla eksperymentatorow)")
print("=" * 80)
print()

predictions = [
    ("Hg1245/SrTiO3 MBE",  3.862, 5, None, None, None, None, None,
        "P6.A Tier 1"),
    ("Hg1223-quench ambient", 3.830, 3, None, None, None, None, None,
        "Juz potwierdzony 151K (Deng 2026)"),
    ("FeSe/BaTiO3",  3.770, None, 8, 75.0, 2.0, 0.25, None,
        "P6.B alternative substrate"),
    ("YbH9 @ 300 GPa", 3.50, None, 8, 200.0, None, None, 0.0,
        "P6.C z P_scale_Yb=552"),
    ("Pd-H rapid-quench", 4.10, None, 12, 60.0, 1.46, 0.10, None,
        "P6.D SF damping z I=0.1 (assumed)"),
]

print(f"  {'Materiał':>22} {'T_pred [K]':>12} {'Note':>12}")
for entry in predictions:
    name = entry[0]
    if len(entry) == 9:
        name, a, n_layers, z, om, N, I, lam_m, note = entry
        if n_layers is not None:
            Tp = Tc_cuprate(a, n_layers)
        else:
            # phonon-mediated
            if N is not None and I is not None:
                lam_sf = lambda_sf_P7(A_map["d"], z, N, I)
            else:
                lam_sf = lam_m if lam_m is not None else 0
            # Dla YbH9 specjalnie: eta = pressure(300, 0, 552)
            if "YbH9" in name:
                eta_val = eta_pressure(300, 0, P_scale_Yb)
                Tp = Tc_phonon(a, eta_val, z, om, lam_sf)
            elif "Pd-H" in name:
                Tp = Tc_phonon(a, "d", z, om, lam_sf)
            else:
                Tp = Tc_phonon(a, "d", z, om, lam_sf)
    print(f"  {name:>22} {Tp:>10.1f} K  {note}")
print()


print("=" * 80)
print("  ps20 complete. P6 + P7.1 master validation zamkniete.")
print("=" * 80)
