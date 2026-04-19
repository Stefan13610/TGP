"""
ps6_cosine_modulation.py  -  Program P5 #6

Cel:
  Sprawdzic czy cosinusowa kara modulacji
      M_cos(a) = 0.5 * (1 + cos(2*pi*(a - a_near_harm)/a*_tgp))
  daje lepszy fit niz gaussowska
      M_gauss(a) = exp(-delta_a^2 / sigma^2)
  uzyta w ps5.

Fizycznie uzasadnienie:
  W siatce periodycznej z a* strukture Blocha opisuje funkcja cosinusowa,
  nie gaussowska. Gauss to lokalne rozwiniecie cos wokol maksimum.
  Dla duzych delta_a cos przewiduje PERIODYCZNE odrodzenie (bliskosc
  kolejnej harmoniki), gauss monotonicznie tlumi.

Plan:
  1. Wariant 6a: (pure-TGP + cos) - A_s=A_e, A_d=A_mu, A_f=A_tau, fit Lambda_E
  2. Wariant 6c: (all-free + cos) - fit {Lambda_E, A_s, A_sp, A_d, A_f}
  3. Porownac r, RMS_log z ps5 (gauss)
  4. Sprawdzic czy outliery (H3S, YBCO) maja lepszy residual
"""

import numpy as np
from scipy.optimize import minimize

# =============================================================
# Stale i P4 fix
# =============================================================

K_B = 8.617333262e-2  # meV/K
A_BOHR_ANGSTROM = 0.52917721067

# Soliton zero-order
g0_e = 0.869470
g0_mu = 1.406833
g0_tau = 1.729615

# Amplitudy P4
A_e_P4 = 0.124587
A_mu_P4 = 0.472198
A_tau_P4 = 0.956027

# TGP harmonika
a_star_tgp_substr = 7.725  # jedn. Bohra
a_star_tgp_A = a_star_tgp_substr * A_BOHR_ANGSTROM  # ~4.088 A
C_0 = 48.8222

# =============================================================
# k_d(z)
# =============================================================

def k_d_from_z(z):
    table = {8: 2.936, 12: 4.403, 6: 2.202, 4: 0.893}
    return table.get(z, 2.936)

# =============================================================
# Modulacje
# =============================================================

def nearest_tgp_harmonic_delta(a_A):
    """Zwroc delta_a = (a - najblizsza n*a_star_tgp_A)."""
    n = round(a_A / a_star_tgp_A)
    n = max(1, n)
    a_harm = n * a_star_tgp_A
    return a_A - a_harm

def M_gauss(a_A, sigma_a):
    d = nearest_tgp_harmonic_delta(a_A)
    return np.exp(-d**2 / sigma_a**2)

def M_cos(a_A):
    """Czysto periodyczny bez wolnego parametru sigma."""
    d = nearest_tgp_harmonic_delta(a_A)
    return 0.5 * (1.0 + np.cos(2.0 * np.pi * d / a_star_tgp_A))

# =============================================================
# Model T_c
# =============================================================

def compute_Tc_cos(a_A, orb, z, Lambda_E_meV, A_map):
    M = M_cos(a_A)
    A = A_map[orb]
    J_star = C_0 * A**2
    Tc_substr = k_d_from_z(z) * J_star * M
    Tc_K = Tc_substr * (Lambda_E_meV * 1e-3) / K_B
    return max(Tc_K, 1e-9)

# =============================================================
# Materialy (19 SC)
# =============================================================

materials = [
    # (nazwa, a[A], orb, z, Tc_obs [K])
    ("Al",    4.046, "s",  12,   1.18),
    ("Zn",    2.665, "s",  12,   0.85),
    ("Sn",    5.831, "sp",  6,   3.72),
    ("In",    3.252, "sp", 12,   3.41),
    ("Hg",    2.989, "sp",  6,   4.15),
    ("Pb",    4.950, "sp", 12,   7.20),
    ("Nb",    3.301, "d",   8,   9.26),
    ("V",     3.024, "d",   8,   5.30),
    ("Ta",    3.303, "d",   8,   4.47),
    ("Tc",    2.735, "d",  12,   7.80),
    ("NbN",   4.394, "d",   6,  16.10),
    ("V3Si",  4.725, "d",  12,  17.10),  # A15
    ("Nb3Sn", 5.290, "d",  12,  18.30),  # A15
    ("Nb3Ge", 5.140, "d",  12,  23.20),  # A15
    ("MgB2",  3.086, "sp", 12,  39.00),
    ("YBCO",  3.820, "d",   6,  92.00),  # in-plane
    ("BiSCCO",3.820, "d",   6, 110.00),
    ("H3S",   3.100, "sp",  6, 203.00),  # 155 GPa
    ("LaH10", 5.100, "f",  12, 250.00),  # 170 GPa
]

def build_A_map(variant, params):
    if variant == "6a":
        Lambda_E = params[0]
        A_map = {"s": A_e_P4, "sp": A_e_P4, "d": A_mu_P4, "f": A_tau_P4}
    elif variant == "6c":
        Lambda_E, A_s, A_sp, A_d, A_f = params
        A_map = {"s": A_s, "sp": A_sp, "d": A_d, "f": A_f}
    else:
        raise ValueError(variant)
    return A_map, Lambda_E

def loss(params, variant):
    A_map, Lambda_E = build_A_map(variant, params)
    if Lambda_E <= 0:
        return 1e12
    total = 0.0
    for name, a_A, orb, z, Tc_obs in materials:
        Tc_p = compute_Tc_cos(a_A, orb, z, Lambda_E, A_map)
        total += (np.log10(Tc_p) - np.log10(Tc_obs)) ** 2
    return total

def stats(variant, params):
    A_map, Lambda_E = build_A_map(variant, params)
    log_obs = []
    log_pred = []
    rows = []
    for name, a_A, orb, z, Tc_obs in materials:
        Tc_p = compute_Tc_cos(a_A, orb, z, Lambda_E, A_map)
        log_obs.append(np.log10(Tc_obs))
        log_pred.append(np.log10(Tc_p))
        rows.append((name, Tc_obs, Tc_p, np.log10(Tc_p) - np.log10(Tc_obs)))
    lo = np.array(log_obs); lp = np.array(log_pred)
    r = np.corrcoef(lo, lp)[0, 1]
    resid = lp - lo
    return r, float(np.sqrt(np.mean(resid**2))), float(np.mean(np.abs(resid))), rows

# =============================================================
# Main
# =============================================================

print("=" * 78)
print("  ps6_cosine_modulation.py")
print("=" * 78)
print()
print("  Program P5 #6: cosinusowa kara modulacji (zamiast gauss w ps5)")
print("  M_cos(a) = 0.5 * (1 + cos(2pi * delta_a / a*_tgp))")
print()

# Fit 6a
print("=" * 78)
print("  Part A. Wariant 6a (pure-TGP + cos): fit tylko Lambda_E")
print("=" * 78)
res_a = minimize(loss, x0=[0.1], args=("6a",), method="Nelder-Mead",
                 options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 2000})
Lambda_E_6a = res_a.x[0]
r6a, rms6a, mae6a, rows6a = stats("6a", res_a.x)
print(f"\n  Lambda_E = {Lambda_E_6a:.4f} meV")
print(f"  r(log-log) = {r6a:.4f}")
print(f"  RMS_log    = {rms6a:.4f}")
print(f"  MAE_log    = {mae6a:.4f}")
print(f"  loss       = {res_a.fun:.4f}")

# Fit 6c
print()
print("=" * 78)
print("  Part B. Wariant 6c (all-free + cos): fit {Lambda_E, A_s, A_sp, A_d, A_f}")
print("=" * 78)
x0 = [0.1, A_e_P4, A_e_P4, A_mu_P4, A_tau_P4]
res_c = minimize(loss, x0=x0, args=("6c",), method="Nelder-Mead",
                 options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 20000})
Lambda_E_6c, A_s_6c, A_sp_6c, A_d_6c, A_f_6c = res_c.x
r6c, rms6c, mae6c, rows6c = stats("6c", res_c.x)
print(f"\n  Lambda_E = {Lambda_E_6c:.4f} meV")
print(f"  A_s      = {A_s_6c:.4f}   (P4 A_e   = {A_e_P4:.4f})")
print(f"  A_sp     = {A_sp_6c:.4f}   (P4 A_e   = {A_e_P4:.4f})")
print(f"  A_d      = {A_d_6c:.4f}   (P4 A_mu  = {A_mu_P4:.4f})")
print(f"  A_f      = {A_f_6c:.4f}   (P4 A_tau = {A_tau_P4:.4f})")
print(f"  r(log-log) = {r6c:.4f}")
print(f"  RMS_log    = {rms6c:.4f}")
print(f"  MAE_log    = {mae6c:.4f}")
print(f"  loss       = {res_c.fun:.4f}")

# Porownanie
print()
print("=" * 78)
print("  Part C. Porownanie ps5 (gauss) vs ps6 (cos)")
print("=" * 78)
print()
print("  Wariant     r(gauss ps5)   r(cos ps6)   RMS(ps5)  RMS(ps6)")
print("  ---------   ------------   ----------   --------  --------")
print(f"  pure-TGP       0.4190       {r6a:.4f}     0.7146    {rms6a:.4f}")
print(f"  all-free       0.4799       {r6c:.4f}     0.6050    {rms6c:.4f}")

# Reszty 6c
print()
print("=" * 78)
print("  Part D. Reszty per-material (6c cos)")
print("=" * 78)
print(f"\n  {'SC':>8} {'obs':>8} {'pred':>8}  {'log10(p/o)':>10}")
print(f"  {'--------':>8} {'--------':>8} {'--------':>8}  {'----------':>10}")
for name, obs, pred, dlog in rows6c:
    print(f"  {name:>8} {obs:>8.2f} {pred:>8.2f}  {dlog:>+10.4f}")

# Werdykt
print()
print("=" * 78)
print("  Part E. Werdykt ps6")
print("=" * 78)
print()
if r6c > 0.55:
    print(f"  Wariant 6c: r = {r6c:.4f}  >> ps5 5c (0.480)")
    print("  Cosinus DAJE istotna poprawe -> ps6 nalezy uwzglednic w dodatku P5.")
elif r6c > 0.490:
    print(f"  Wariant 6c: r = {r6c:.4f}  marginalnie lepiej niz ps5 5c (0.480)")
    print("  Cosinus nie jest kluczowy, ale nieco lepsze dopasowanie.")
else:
    print(f"  Wariant 6c: r = {r6c:.4f}  nie lepiej niz ps5 5c (0.480)")
    print("  Cosinus nie naprawia problemu. Brak efektu.")

# Sprawdz czy outliery
h3s = [row for row in rows6c if row[0] == "H3S"][0]
ybco = [row for row in rows6c if row[0] == "YBCO"][0]
mgb2 = [row for row in rows6c if row[0] == "MgB2"][0]
print()
print(f"  Outliery ps5:                        ps6 (cos):")
print(f"    H3S:   pred=5.46,  off=-1.57   ->  pred={h3s[2]:.2f}, off={h3s[3]:+.2f}")
print(f"    YBCO:  pred=15.06, off=-0.79   ->  pred={ybco[2]:.2f}, off={ybco[3]:+.2f}")
print(f"    MgB2:  pred=8.16,  off=-0.68   ->  pred={mgb2[2]:.2f}, off={mgb2[3]:+.2f}")

print()
print("=" * 78)
print("  ps6 complete.")
print("=" * 78)
