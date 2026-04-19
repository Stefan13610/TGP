"""
ps10_cuprates_2d_dwave.py  -  Program P6.A #2

Cel:
  Cuprates (YBCO, BiSCCO, Hg-1223, LSCO, Tl-2212) maja T_c = 35-138 K
  przy ambient pressure. Model ps5 daje tylko 15 K. Fit zawodzi o ~6x.

Hipoteza fizyczna:
  Cuprates to quasi-2D SC w plaszczyznach CuO2. Mechanizm:
    1. Zhang-Rice singlet: Cu 3d(x^2-y^2) + O 2p (in-plane) -> hybryda
       A_ZR^2 = A_d^2 + 2 * A_sp^2  (2 tlenki na Cu)
    2. d-wave parowanie: gap Delta(k) = Delta_0 * (cos(kx a) - cos(ky a))
       Wzmacnia sprzezenie w kierunkach (1,0) i (0,1) Fermi surface.
    3. Van Hove singularity: DOS diverges przy optymalnym domieszkaniu
       -> effective lambda_eff * (1 + ln(DOS)) - ekwiwalent boost ~ 2-3
    4. Quasi-2D: inter-planar weak coupling utrzymuje 3D LRO,
       ale k_d jest enhancowane vs pure 3D BKT.

Model TGP P6.A:
  T_c^cuprate = K_dw * k_d^eff(z) * C_0 * A_ZR^2 * M(a_inplane) * Lambda_E

  gdzie K_dw to d-wave+vH booster (fit).
  A_ZR^2 = A_d^2 + N_O * A_p^2  (N_O = 2 in CuO2 plane).
  A_p ~ A_sp (no separate sp/p distinction in P4).

Plan:
  1. Fit K_dw * Lambda_E_eff razem na 5-6 cupratach przy ambient.
  2. Porownaj predykcje z exp.
  3. Sprawdz uniwersalnosc K_dw - ten sam wspolczynnik dla wszystkich cuprat'ow?
"""

import numpy as np
from scipy.optimize import minimize

# Stale
K_B = 8.617333e-5  # eV/K
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088 A

A_e_P4  = 0.124587
A_mu_P4 = 0.472198
A_tau_P4 = 0.956027

Lambda_E_5c = 0.1309  # meV
A_map_5c = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}
sigma_a_5c = 2.5856
C_0 = 48.8222

def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A, sigma=sigma_a_5c):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma**2)

# Zhang-Rice amplitude
A_d = A_map_5c["d"]
A_p = A_map_5c["sp"]
A_ZR_sq = A_d**2 + 2.0 * A_p**2
print(f"  Zhang-Rice: A_ZR^2 = A_d^2 + 2*A_p^2 = {A_d**2:.4f} + 2*{A_p**2:.4f} = {A_ZR_sq:.4f}")
print(f"              A_ZR   = {np.sqrt(A_ZR_sq):.4f}  (vs A_d = {A_d:.4f})")

# =============================================================
# Baza cuprates (ambient pressure, optymalne domieszkanie)
# =============================================================

# (nazwa, a_inplane [A], n_CuO2_layers, T_c_obs [K], komentarz)
cuprates = [
    ("La214",    3.780, 1,  38.0, "LSCO opt x=0.15, 1 CuO2/cell"),
    ("YBCO",     3.820, 2,  92.0, "YBa2Cu3O7, 2 CuO2/cell"),
    ("BiSCCO",   3.820, 2, 110.0, "Bi2Sr2CaCu2O8+d, 2 CuO2"),
    ("Tl2212",   3.855, 2, 105.0, "Tl2Ba2CaCu2O8, 2 CuO2"),
    ("Hg1223",   3.855, 3, 138.0, "HgBa2Ca2Cu3O8, 3 CuO2 - record ambient!"),
    ("Tl2223",   3.855, 3, 127.0, "Tl2Ba2Ca2Cu3O10, 3 CuO2"),
    ("Nd214",    3.940, 1,  24.0, "Nd2CuO4 (electron-dop.)"),
    ("Bi2201",   3.820, 1,  34.0, "Bi2Sr2CuO6+d, 1 CuO2"),
]

# =============================================================
# Model Tc
# =============================================================

def Tc_cuprate(a_A, n_layers, K_dw, Lambda_E_meV, z_planar=8):
    """
    T_c^cuprate = K_dw * k_d(z_planar) * C_0 * A_ZR^2 * M(a) * layer_factor * Lambda_E
    layer_factor = (n_layers)^p   p~0.5 (saturation)
    """
    layer_factor = n_layers ** 0.5   # empiryczne: T_c(3-layer)/T_c(1-layer) ~ sqrt(3)
    M = M_gauss(a_A)
    Tc_substr = K_dw * k_d(z_planar) * C_0 * A_ZR_sq * M * layer_factor
    return Tc_substr * (Lambda_E_meV * 1e-3) / K_B

def loss(params):
    K_dw, Lambda_E = params
    if K_dw <= 0 or Lambda_E <= 0:
        return 1e10
    total = 0.0
    for name, a, n, Tobs, _ in cuprates:
        Tp = Tc_cuprate(a, n, K_dw, Lambda_E)
        total += (np.log10(max(Tp, 1e-6)) - np.log10(Tobs))**2
    return total

# =============================================================
# Fit
# =============================================================

print()
print("=" * 78)
print("  ps10: Fit modelu cuprates 2D-ZR-dwave")
print("=" * 78)
print()

res = minimize(loss, x0=[3.0, 0.2], method="Nelder-Mead",
               options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 5000})
K_dw_fit, Lam_fit = res.x

# Statystyki
log_obs, log_pred = [], []
print(f"  Fit: K_dw = {K_dw_fit:.3f},  Lambda_E = {Lam_fit:.4f} meV")
print()
print(f"  {'Cuprate':>9} {'a[A]':>5} {'n':>3} {'T_obs':>7} {'T_pred':>7}  {'log10(p/o)':>11}")
print(f"  {'---------':>9} {'-----':>5} {'---':>3} {'-------':>7} {'-------':>7}  {'-----------':>11}")

for name, a, n, Tobs, comment in cuprates:
    Tp = Tc_cuprate(a, n, K_dw_fit, Lam_fit)
    log_obs.append(np.log10(Tobs))
    log_pred.append(np.log10(max(Tp, 1e-6)))
    print(f"  {name:>9} {a:>5.3f} {n:>3d} {Tobs:>7.1f} {Tp:>7.2f}  {np.log10(max(Tp,1e-6))-np.log10(Tobs):>+11.4f}")

lo, lp = np.array(log_obs), np.array(log_pred)
r = np.corrcoef(lo, lp)[0, 1]
rms = float(np.sqrt(np.mean((lp - lo)**2)))
mae = float(np.mean(np.abs(lp - lo)))
print()
print(f"  r(log-log) = {r:.4f}")
print(f"  RMS_log    = {rms:.4f}")
print(f"  MAE_log    = {mae:.4f}")

# =============================================================
# Interpretacja K_dw_fit
# =============================================================
print()
print("=" * 78)
print("  Part B. Interpretacja K_dw")
print("=" * 78)
print()
print(f"  K_dw = {K_dw_fit:.3f}")
print()
print("  Sklad K_dw ~ (d-wave boost) * (vH enhancement) * (ZR correction):")
print("  - d-wave boost: ~1.5 (kalkulacja BCS na Fermi surface z D_0(cos kx - cos ky))")
print("  - vH enhancement: ln(D(E_F)/D_0) ~ 2-3 dla domieszkan optymalnych")
print("  - ZR correction: inkl. oxygen-bridge mediated hopping ~ 1.2")
print(f"  Oczekiwane: 1.5 * 2.5 * 1.2 = 4.5. Fit: {K_dw_fit:.2f}.")

# =============================================================
# Part C. Predykcja wysokotemp. przy ambient - skan
# =============================================================
print()
print("=" * 78)
print("  Part C. Predykcja: kandydaci ambient-pressure high-Tc")
print("=" * 78)
print()
print("  Szukamy materialow z:")
print("    a_inplane bliskie a*_TGP = 4.088 A")
print("    CuO2-like plaszczyzny (d-orbital przewodnictwa)")
print("    n_layers = 3-5 (wiecej -> wyzsze T_c)")
print()
print("  Propozycje literaturowe (eksperymentalne targets):")

# Kandydaci: nickelates, inne layered d-SC
candidates = [
    # (nazwa, a, n, komentarz, obs_if_any)
    ("NdNiO2-inf", 3.916, 1, "infinite-layer nickelate", "T_c ~10K (thin film)"),
    ("La3Ni2O7",   3.830, 2, "bilayer Nd nickelate P>14 GPa", "T_c=80K @ 14GPa"),
    ("Sr2CuO3",    3.910, 1, "simple Cu-chain-oxide", "non-SC"),
    ("HgCuO2-max", 3.855, 5, "hypothetic 5-layer Hg cuprate", "nie syntetyzowane"),
    ("HgCuO2-inf", 3.855,10, "HYPOTHESIS infinite-layer Hg", "przy akademicka spekulacja"),
    ("(Ba,Sr)CuO2",4.010, 1, "apical-free parent", "non-SC ambient"),
    ("Na-Tl2212",  3.800, 2, "alkali-doped Tl", "mozliwe"),
    ("La114-ovd",  3.750, 1, "LSCO overdoped x=0.25", "T_c=25K"),
]

print(f"  {'Material':>12} {'a[A]':>5} {'n':>3} {'T_pred[K]':>9}  {'komentarz':>40}")
print(f"  {'------------':>12} {'-----':>5} {'---':>3} {'---------':>9}  {'----------------------------------------':>40}")
for name, a, n, desc, exp_info in candidates:
    Tp = Tc_cuprate(a, n, K_dw_fit, Lam_fit)
    info = f"{desc}; exp: {exp_info}"
    print(f"  {name:>12} {a:>5.3f} {n:>3d} {Tp:>9.2f}  {info[:40]:>40}")

# =============================================================
# Part D. Sweet spot - idealizowany cuprate
# =============================================================
print()
print("=" * 78)
print("  Part D. Idealny kandydat: maksymalizacja T_c przy ambient")
print("=" * 78)
print()

# Scan a_inplane and n_layers
a_scan = np.linspace(3.60, 4.30, 15)
n_scan = [1, 2, 3, 4, 5, 6, 8, 10]

print(f"  {'a[A]':>6}", end="")
for n in n_scan:
    print(f" {'n='+str(n):>6}", end="")
print()

for a in a_scan:
    print(f"  {a:>6.3f}", end="")
    for n in n_scan:
        Tp = Tc_cuprate(a, n, K_dw_fit, Lam_fit)
        print(f" {Tp:>6.1f}", end="")
    print()

# Maksimum
print()
best = (0, 0, 0)  # (Tc, a, n)
for a in np.linspace(3.50, 4.50, 100):
    for n in range(1, 11):
        Tp = Tc_cuprate(a, n, K_dw_fit, Lam_fit)
        if Tp > best[0]:
            best = (Tp, a, n)
print(f"  Maksimum mapy: T_c = {best[0]:.1f} K przy a = {best[1]:.3f} A, n = {best[2]}")
print(f"  Dla porownania: Hg1223 (obs rekord ambient) = 138 K przy a=3.855, n=3.")

# =============================================================
# Part E. Wniosek koncowy
# =============================================================
print()
print("=" * 78)
print("  Part E. Werdykt ps10")
print("=" * 78)
print()
if r > 0.85:
    print(f"  r(log-log) = {r:.3f} >> ps5 5c 0.48 - BARDZO DOBRY fit cuprates.")
    print("  Model P6.A dziala dla high-Tc cupratow przy ambient.")
elif r > 0.70:
    print(f"  r(log-log) = {r:.3f} > ps5 5c 0.48 - DOBRY fit cuprates.")
elif r > 0.50:
    print(f"  r(log-log) = {r:.3f} > ps5 5c 0.48 - umiarkowany fit.")
else:
    print(f"  r(log-log) = {r:.3f} - slaby fit, model wymaga dalszej pracy.")

print()
print(f"  K_dw universal = {K_dw_fit:.2f}")
print(f"  Lambda_E cuprate = {Lam_fit:.3f} meV  (vs ps5 5c: {Lambda_E_5c} meV)")
print()
if abs(Lam_fit - Lambda_E_5c) / Lambda_E_5c < 0.3:
    print("  Lambda_E konsystentne z ps5 - uniwersalna stala TGP.")
else:
    print(f"  Lambda_E cuprate != Lambda_E BCS ({Lam_fit/Lambda_E_5c:.2f}x).")
    print("  Sygnal ze Lambda_E zalezy od klasy SC -> P6.B (coupling-dependent).")

print()
print("=" * 78)
print("  ps10 complete.")
print("=" * 78)
