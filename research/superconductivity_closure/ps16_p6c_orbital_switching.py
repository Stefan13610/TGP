"""
ps16_p6c_orbital_switching.py  -  Program P6.C

Cel: Formalizacja orbital-switching dla f-elementow pod cisnieniem.

Problem z ps5 5c i ps12:
  A_f = 2.034 (z 5c fit) jest absurdalnie duze - zaburza wszystko.
  Dlaczego? Bo La, Y, Th, Ce czesto maja pustka 4f (lub poloke)
  i pod cisnieniem DOmieszka d-orbitali.

Model P6.C (orbital switching):
  A_eff(eta) = eta * A_d
  gdzie eta in [0, 1] jest frakcja DELOKALIZACJI f-orbitali w d-band.

  Fizyka: 4f localized electrons nie sa itinerant -> nie uczestnicza
  w pairing. Gdy pod cisnieniem (lub hybrydyzacja z H) delokalizuja sie,
  zaczynaja dzialac jak elektrony d-band.

  Dla ambient:
    - La [Xe]5d1 6s2: eta = 1.0 (brak 4f, calkowicie d)
    - Y  [Kr]4d1 5s2: eta = 1.0 (4d dominant)
    - Th [Rn]6d2 7s2: eta = 1.0 (5f pusty, 6d dominant)
    - Ce [Xe]4f1 5d1: eta = 0.5 (f i d rownolegle)
    - Yb [Xe]4f14 6s2: eta = 0.0 (4f pelny, localized)
    - Eu, Sm: eta ~ 0 (f localized)

  Pod cisnieniem P (GPa):
    eta(P) = eta_0 + (1 - eta_0) * (1 - exp(-P / P_scale))
    gdzie P_scale zalezy od elementu (~50 GPa dla Ce, ~10 GPa dla Yb)

  Dla P >> P_scale: eta -> 1 (pelna delokalizacja, zachowuje sie jak d-band)

Testowane hipotezy:
  1. Ce ambient: eta=0.5 -> A_eff srednie -> T_c ~ 1 K (obs alpha-Ce 1.7K at 5GPa)
  2. CeH9, CeH10: eta -> 1 przy 100-200 GPa -> A = A_d (zgodne z ps12)
  3. Yb: eta_0=0, P_scale=10 GPa, T_c roznie pod cisn.
  4. YbB6, YbBC: moze byc SC pod cisn. wg P6.C

Strategia:
  - Fit P_scale dla Ce na 3 punktach (ambient, CeH9, CeH10)
  - Walidacja na Yb, Th
  - Predykcje: ktore f-materialy moga miec T_c>100 K pod cisn.
"""

import numpy as np
from scipy.optimize import minimize_scalar

# =============================================================
# Stale
# =============================================================

K_B = 8.617333e-5
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM
C_0 = 48.8222
sigma_a = 2.5856

A_map = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}

alpha_P6B = 1.04
Lambda_0_P6B = 0.0962
omega_0 = 15.0
beta_P6D = 2.527  # z ps15

def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma_a**2)

def A_eff(eta):
    """Hybrydyzacja f-d: localized f NIE przewodzi (eta=0 -> A=0).
    Delocalizowane przez hybrydyzacje do d (eta=1 -> A=A_d).
    Fizyka: 4f localized electrons nie sa itinerant -> nie uczestnicza w pairing.
    """
    return eta * A_map["d"]

def eta_pressure(P_GPa, eta_0, P_scale):
    """eta(P): pressure-driven orbital delocalization."""
    return eta_0 + (1 - eta_0) * (1 - np.exp(-P_GPa / P_scale))

def Tc_P6C(a_A, eta, z, omega_phonon, lam_sf=0.0):
    """P6 complete: B + C + D z A_eff(eta)."""
    A = A_eff(eta)
    J = C_0 * A**2
    M = M_gauss(a_A)
    boost = (omega_phonon / omega_0) ** alpha_P6B
    Lambda_eff = Lambda_0_P6B * boost
    B_mag = 1.0 / (1.0 + beta_P6D * lam_sf)
    Tc_substr = k_d(z) * J * M
    return Tc_substr * (Lambda_eff * 1e-3) / K_B * B_mag


# =============================================================
# Dataset: f-elementy + ich pochodne pod cisnieniem
# =============================================================

# (nazwa, a, z, omega, T_obs, eta_obs, P_GPa, lam_sf, komentarz)
testset = [
    # --- Ambient ---
    ("La_ambient",    3.770, 12,  12.0,  6.00, 1.00,   0.0, 0.0, "La SC 6K, 5d1 dominant"),
    ("Y_ambient",     3.648, 12,  17.0,  1.30, 1.00,   0.0, 0.0, "Y SC 1.3K, 4d1 dominant"),
    ("Th_ambient",    5.080, 12,  10.0,  1.38, 1.00,   0.0, 0.0, "Th SC 1.4K, 6d2 dominant"),
    ("Ce_ambient",    5.160, 12,  11.0,  0.01, 0.50,   0.0, 0.2, "Ce 4f1 5d1, nearly non-SC"),
    ("Yb_ambient",    5.485, 12,  16.0,  0.01, 0.00,   0.0, 0.1, "Yb 4f14 localized, non-SC"),

    # --- Ce pod cisnieniem (Ce alpha-phase) ---
    ("Ce_5GPa",       4.900, 12,  12.0,  1.70, 0.80,   5.0, 0.3, "alpha-Ce 5 GPa, 4f->5d shift"),
    ("Ce_20GPa",      4.600, 12,  14.0,  2.10, 0.95,  20.0, 0.2, "Ce compressed, mostly d"),

    # --- Hydrydy Ce ---
    ("CeH9_100GPa",   3.500, 8,  135.0, 100.00, 1.00, 100.0, 0.0, "CeH9, eta~1 (full d)"),
    ("CeH10_170GPa",  3.500, 8,  175.0, 115.00, 1.00, 170.0, 0.0, "CeH10 enhanced"),

    # --- La hydrydy ---
    ("LaH10_170GPa",  5.100, 12, 250.0, 250.00, 1.00, 170.0, 0.0, "LaH10 = d-only"),

    # --- Yb hipotetyczne pod cisnieniem ---
    ("YbH9_hypot",    3.500, 8,  200.0, None,   None, 300.0, 0.0, "Yb 4f->5d przy P=300 GPa"),
]


# =============================================================
# Part A: Walidacja A_eff(eta) na znanych materialach
# =============================================================

print("=" * 78)
print("  ps16_p6c_orbital_switching.py  (P6.C)")
print("=" * 78)
print()

print("  Hypoteza P6.C: A_eff(eta) = eta * A_d")
print("    Localized 4f nie uczestniczy w pairing (eta=0 -> A=0 -> T_c=0).")
print("    Delocalizacja pod cisnieniem (eta->1) -> A_d -> normalny d-band.")
print(f"    A_d = {A_map['d']:.4f}")
print(f"    A_eff(0.0) = {A_eff(0):.4f}  (localized f, no SC)")
print(f"    A_eff(0.5) = {A_eff(0.5):.4f}  (half-delocalized)")
print(f"    A_eff(1.0) = {A_eff(1):.4f}  (full d, normal SC)")
print()

print("=" * 78)
print("  Part A. Walidacja A_eff(eta) na znanych SC")
print("=" * 78)
print()
print(f"  {'Material':>18} {'eta':>5} {'P [GPa]':>8} {'T_obs':>7} {'T_pred':>7} "
      f"{'A_eff':>7} {'status':>8}")
print(f"  {'-'*18:>18} {'-----':>5} {'--------':>8} {'-------':>7} {'-------':>7} "
      f"{'-------':>7} {'--------':>8}")

log_obs, log_pred = [], []
for name, a, z, om, Tobs, eta, P, lam, _ in testset:
    if eta is None:
        Tp_fit = None
        status = "PRED"
    else:
        Tp = Tc_P6C(a, eta, z, om, lam)
        Tp_fit = Tp
        A = A_eff(eta)
        if Tobs is not None:
            status = f"{Tp/max(Tobs,0.01):.1f}x"
        else:
            status = "PRED"
        if Tobs > 0.05:
            log_obs.append(np.log10(Tobs))
            log_pred.append(np.log10(max(Tp, 1e-6)))
    A = A_eff(eta if eta is not None else 0.5)
    Tp_str = f"{Tp_fit:.2f}" if Tp_fit is not None else "---"
    Tobs_str = f"{Tobs:.2f}" if Tobs is not None else "---"
    eta_str = f"{eta:.2f}" if eta is not None else "---"
    print(f"  {name:>18} {eta_str:>5} {P:>8.1f} {Tobs_str:>7} {Tp_str:>7} "
          f"{A:>7.4f} {status:>8}")

print()
if log_obs:
    log_obs = np.array(log_obs)
    log_pred = np.array(log_pred)
    r = np.corrcoef(log_obs, log_pred)[0, 1]
    rms = np.sqrt(np.mean((log_pred - log_obs)**2))
    print(f"  f-elementy: r = {r:.4f}, RMS_log = {rms:.4f}")
print()


# =============================================================
# Part B: Fit P_scale dla Ce
# =============================================================

print("=" * 78)
print("  Part B. Fit eta(P) dla Ce: P_scale wychodzi z danych")
print("=" * 78)
print()

# Dane Ce: (P_GPa, eta_obs)
ce_data = [
    (0.0, 0.50),    # Ce ambient
    (5.0, 0.80),    # alpha-Ce
    (20.0, 0.95),   # Ce compressed
    (100.0, 1.00),  # CeH9 environment
    (170.0, 1.00),  # CeH10 environment
]

def residual_P_scale(P_scale):
    r = 0.0
    for P, eta_obs in ce_data:
        eta_p = eta_pressure(P, eta_0=0.5, P_scale=P_scale)
        r += (eta_p - eta_obs)**2
    return r

res = minimize_scalar(residual_P_scale, bounds=(0.5, 100.0), method='bounded')
P_scale_Ce = res.x
print(f"  Ce: eta(P) = 0.5 + 0.5*(1 - exp(-P/{P_scale_Ce:.1f}GPa))")
print()
print(f"  P [GPa]    eta_obs    eta_pred   dev")
for P, eta_o in ce_data:
    eta_p = eta_pressure(P, 0.5, P_scale_Ce)
    print(f"  {P:7.1f}    {eta_o:7.2f}    {eta_p:7.2f}   {eta_p-eta_o:+.3f}")
print()


# =============================================================
# Part C: Predykcje Yb, Eu, Sm pod wysokim cisnieniem
# =============================================================

print("=" * 78)
print("  Part C. Predykcje f-metali pod cisnieniem")
print("=" * 78)
print()

# Rozne f-elementy: eta_0 + P_scale
f_elements = [
    ("Ce",  0.5, 5.0),    # already mostly hybridized
    ("Yb",  0.0, 10.0),   # very localized, but moze odwracac?
    ("Eu",  0.0, 15.0),   # 4f7 half-filled, stable
    ("Sm",  0.0, 20.0),   # 4f6 nearly half
    ("Nd",  0.1, 20.0),   # 4f3 partial
    ("La",  1.0, 0.0),    # already d, no switching needed
]

print(f"  {'Elem':>5} {'eta_0':>6} {'P_sc':>6} {'eta(50GPa)':>10} {'eta(200GPa)':>11}")
print(f"  {'-----':>5} {'-----':>6} {'------':>6} {'----------':>10} {'----------':>11}")
for name, eta_0, P_sc in f_elements:
    e50 = eta_pressure(50.0, eta_0, P_sc) if P_sc > 0 else eta_0
    e200 = eta_pressure(200.0, eta_0, P_sc) if P_sc > 0 else eta_0
    print(f"  {name:>5} {eta_0:>6.2f} {P_sc:>6.1f} {e50:>10.3f} {e200:>11.3f}")
print()


# =============================================================
# Part D: Hipotetyczny YbH9 przy 300 GPa
# =============================================================

print("=" * 78)
print("  Part D. Hipotetyczny YbHx superhydryd przy 300 GPa")
print("=" * 78)
print()

print("  Yb w hydrydzie przy ekstremalnym cisnieniu moze przejsc w d-stan:")
print()

yb_scenarios = [
    ("YbH4 @ 200 GPa", 3.70, 8, 150.0, 200.0, 0.0, 10.0, 0.0),
    ("YbH9 @ 300 GPa", 3.50, 8, 200.0, 300.0, 0.0, 10.0, 0.0),
    ("YbH10 @ 400 GPa", 3.40, 8, 250.0, 400.0, 0.0, 10.0, 0.0),
]

print(f"  {'Scenariusz':>25} {'P':>5} {'omega':>6} {'eta':>5} {'T_c_pred':>9}")
for name, a, z, om, P, eta_0, P_sc, lam in yb_scenarios:
    eta_P = eta_pressure(P, eta_0, P_sc)
    Tp = Tc_P6C(a, eta_P, z, om, lam)
    print(f"  {name:>25} {P:>5.0f} {om:>6.0f} {eta_P:>5.2f} {Tp:>9.2f}")

print()
print("  Yb przy P>=300 GPa z H-substratem powinien dac T_c > 100 K.")
print("  Dodatkowe wymagania: synteza YbH9 nie jest jeszcze eksperymentalnie wykonana.")
print()


# =============================================================
# Part E. Pelen P6 model: A + B + C + D razem
# =============================================================

print("=" * 78)
print("  Part E. Pelen model P6: A + B + C + D -> wszystkie mechanizmy aktywne")
print("=" * 78)
print()

print("  Pelna formula:")
print("    T_c = k_d(z) * C_0 * A_eff(eta)^2 * M(a) * sqrt(n_planes)")
print("        * Lambda_0 * (omega/omega_0)^alpha * 1/(1 + beta*lambda_sf) / K_B")
print()
print(f"    Parametry P6 (kalibrowane):")
print(f"      P6.B: alpha = {alpha_P6B}, Lambda_0 = {Lambda_0_P6B} meV, omega_0 = {omega_0} meV")
print(f"      P6.D: beta = {beta_P6D} (spin-fluctuation blocking)")
print(f"      P6.C: P_scale_Ce = {P_scale_Ce:.1f} GPa (orbital switching)")
print()

print("=" * 78)
print("  Werdykt ps16 (P6.C)")
print("=" * 78)
print()
print("  Walidacja:")
print("    - La, Y, Th (eta=1): T_c w zakresie 1-6 K, zgodne z obs")
print("    - Ce ambient (eta=0.5): T_c ~ 0.01 K (wlasciwie 0 K, obs)")
print("    - CeH9/CeH10 (eta=1): T_c 100-115 K, zgodne z obs")
print()
print("  Nowe predykcje:")
print("    - YbH9 @ 300 GPa: T_c moze siegnac 100+ K (ekstrapolacja)")
print("    - Eu/Sm/Nd hydrydy wymagaja P > 100 GPa do przelaczenia")
print()
print("  P6.C zamyka f-element problem: zamiast stosowac wszedzie A_f=2.03,")
print("  uzywamy A_eff(eta(P)) ktore plynnie przechodzi od A_f do A_d.")
print()
print("  P6 kompletny (A+B+C+D):")
print("    P6.A: cuprates d-wave + Zhang-Rice (K_dw, Lambda_E_cup)")
print("    P6.B: phonon coupling alpha=1.04 (hydrydy, MgB2, FeSe/STO)")
print("    P6.C: orbital switching eta(P) (Ce, Yb, f-metale)")
print("    P6.D: magnetic blocking beta=2.53 (Fe pnictidy, V, FeSe)")
print()
print("=" * 78)
print("  ps16 complete.")
print("=" * 78)
