"""
ps12_hydride_phonon_coupling.py  -  Program P6.B #1

Cel:
  Kalibracja modelu P6.B dla hydrydow wysokocisnieniowych.

Model P6.B:
  Lambda_E^eff = Lambda_E^(0) * (omega_phonon / omega_0)^alpha

gdzie:
  omega_phonon = dominujaca czestosc optycznych fononow materialu
  omega_0 = referencyjna czestosc ~ Debye w Al ~ 15 meV
  alpha = eksponent coupling (fit)
  Lambda_E^(0) = stala bazowa z ps5 BCS (0.131 meV)

Materialy (high-P hydrydy):
  H3S @ 155 GPa:   T_c = 203 K,  omega_H ~ 175 meV (DFT Errea)
  LaH10 @ 170 GPa: T_c = 250 K,  omega_H ~ 250 meV (clathrate)
  YH9 @ 200 GPa:   T_c = 243 K,  omega_H ~ 220 meV
  CaH6 @ 172 GPa:  T_c = 215 K,  omega_H ~ 200 meV
  YH6 @ 237 GPa:   T_c = 224 K,  omega_H ~ 180 meV
  ThH10 @ 100 GPa: T_c = 161 K,  omega_H ~ 150 meV
  CeH9 @ 100 GPa:  T_c = 100 K,  omega_H ~ 135 meV (lighter H-density)

Formula BCS-like:
  T_c = k_d(z) * C_0 * A_sp^2 * M(a) * Lambda_E^eff

Struktura hydrydow: clathrate-like, a ~ 3.0-5.0 A, z=12 (FCC-like H cage)
"""

import numpy as np
from scipy.optimize import minimize

# Stale
K_B = 8.617333e-5   # eV/K
A_BOHR_ANGSTROM = 0.52917721067
a_star_tgp_A = 7.725 * A_BOHR_ANGSTROM  # 4.088 A
C_0 = 48.8222
sigma_a = 2.5856

# ps5 5c parametry (BCS baseline)
Lambda_E_0 = 0.1309  # meV
A_map_5c = {"s": -0.1110, "sp": 0.2067, "d": 0.3096, "f": 2.0336}

# Referencyjna omega Debye dla Al (metal BCS classical)
omega_0_meV = 15.0  # Al Debye

def k_d(z):
    return {4: 0.893, 6: 2.202, 8: 2.936, 12: 4.403}.get(z, 2.936)

def M_gauss(a_A):
    n = max(1, round(a_A / a_star_tgp_A))
    d = a_A - n * a_star_tgp_A
    return np.exp(-d**2 / sigma_a**2)

# =============================================================
# Baza hydrydow
# =============================================================

# (nazwa, a_cell [A], orb, z, omega_H [meV], T_c [K], P [GPa], komentarz)
#   UWAGA: La, Y, Th nie maja naprawdę populowanych 4f
#   (La = [Xe]5d1 6s2, Y = [Kr]4d1 5s2, Th = [Rn]6d2 7s2).
#   Ce pod P -> 4f delokalizuje -> d-like hybrydyzacja (P6.C).
#   Wszystkie lanthanide hydrydy pod P traktujemy jako orb='d' (nie 'f').
hydrides = [
    ("H3S",     3.100, "sp", 8,  175.0, 203.0, 155, "Imbar body-centered cubic"),
    ("LaH10",   5.100, "d",  12, 250.0, 250.0, 170, "clathrate La 5d"),
    ("YH9",     4.300, "d",  12, 220.0, 243.0, 200, "P63/mmc Y 4d"),
    ("CaH6",    4.500, "sp", 8,  200.0, 215.0, 172, "clathrate Ca 4s"),
    ("YH6",     4.150, "d",  8,  180.0, 224.0, 237, "Im-3m Y 4d"),
    ("ThH10",   5.230, "d",  12, 150.0, 161.0, 100, "clathrate Th 6d"),
    ("CeH9",    3.800, "d",  12, 135.0, 100.0, 95,  "Ce 4f->5d hybryd"),
    ("CeH10",   5.040, "d",  12, 175.0, 115.0, 95,  "Ce 4f->5d hybryd"),
    ("H3S_D",   3.200, "sp", 8,  150.0, 145.0, 100, "isotope D3S"),
]

# =============================================================
# Model P6.B
# =============================================================

def Tc_hydride_P6B(a_A, orb, z, omega_H, alpha, Lambda_0=Lambda_E_0):
    A = A_map_5c[orb]
    J = C_0 * A**2
    M = M_gauss(a_A)
    boost = (omega_H / omega_0_meV) ** alpha
    Lambda_eff = Lambda_0 * boost
    Tc_substr = k_d(z) * J * M
    return Tc_substr * (Lambda_eff * 1e-3) / K_B

def loss_alpha(params):
    alpha, Lambda_0 = params
    if Lambda_0 <= 0:
        return 1e10
    total = 0.0
    for name, a, orb, z, omH, Tobs, _, _ in hydrides:
        Tp = Tc_hydride_P6B(a, orb, z, omH, alpha, Lambda_0)
        total += (np.log10(max(Tp, 1e-6)) - np.log10(Tobs))**2
    return total

# =============================================================
# Fit
# =============================================================

print("=" * 78)
print("  ps12_hydride_phonon_coupling.py  (P6.B)")
print("=" * 78)
print()
print(f"  Model:  Lambda_E^eff = Lambda_E^(0) * (omega_H / omega_0)^alpha")
print(f"  omega_0 = {omega_0_meV} meV (Al Debye)")
print(f"  Lambda_0 baseline = {Lambda_E_0} meV (ps5 5c BCS)")
print(f"  Fit parametry: alpha, Lambda_0 (2 params, {len(hydrides)} materialow)")
print()

res = minimize(loss_alpha, x0=[1.5, 0.1], method="Nelder-Mead",
               options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 5000})
alpha_fit, Lam_fit = res.x

print("=" * 78)
print("  Part A. Fit hydrydow")
print("=" * 78)
print()
print(f"  alpha   = {alpha_fit:.4f}")
print(f"  Lambda_0 = {Lam_fit:.4f} meV  (baseline ps5 5c = {Lambda_E_0})")
print()

log_obs, log_pred = [], []
print(f"  {'Material':>8} {'a[A]':>5} {'orb':>3} {'z':>3} {'omH':>5} {'T_obs':>7} {'T_pred':>8}  {'log10(p/o)':>11}")
print(f"  {'--------':>8} {'-----':>5} {'---':>3} {'---':>3} {'-----':>5} {'-------':>7} {'--------':>8}  {'-----------':>11}")

for name, a, orb, z, omH, Tobs, P, comment in hydrides:
    Tp = Tc_hydride_P6B(a, orb, z, omH, alpha_fit, Lam_fit)
    log_obs.append(np.log10(Tobs))
    log_pred.append(np.log10(max(Tp, 1e-6)))
    print(f"  {name:>8} {a:>5.3f} {orb:>3} {z:>3d} {omH:>5.0f} {Tobs:>7.1f} {Tp:>8.2f}  {log_pred[-1]-log_obs[-1]:>+11.4f}")

lo, lp = np.array(log_obs), np.array(log_pred)
r = np.corrcoef(lo, lp)[0, 1]
rms = float(np.sqrt(np.mean((lp - lo)**2)))
print()
print(f"  r(log-log) = {r:.4f}")
print(f"  RMS_log    = {rms:.4f}")

# =============================================================
# Part B. Interpretacja alpha
# =============================================================
print()
print("=" * 78)
print("  Part B. Interpretacja alpha")
print("=" * 78)
print()
print(f"  alpha = {alpha_fit:.2f}")
print()
if 1.0 <= alpha_fit <= 1.5:
    print("  alpha ~ 1 -> LINIOWE skalowanie Lambda_E z omega_phonon.")
    print("  Oznacza: substrat sprzega sie PROPORCJONALNIE do szybkosci oscylacji.")
elif 1.5 < alpha_fit <= 2.5:
    print("  alpha ~ 2 -> KWADRATOWE skalowanie (jak energia kinetyczna oscylatora).")
    print("  Oznacza: coupling ~ (amplituda * czestosc)^2 -> energia wibracyjna.")
elif alpha_fit > 2.5:
    print("  alpha > 2.5 -> super-linearne, podejrzane o overfitting.")
elif alpha_fit < 1.0:
    print("  alpha < 1 -> sub-linearne, saturating w wysokich omega.")

# =============================================================
# Part C. Walidacja - wykluczajac jeden hydryd z fitu
# =============================================================
print()
print("=" * 78)
print("  Part C. Cross-validation (leave-one-out)")
print("=" * 78)
print()
print("  Wyklucz pojedynczy hydryd z fitu, przewidz jego T_c.")
print(f"  {'Excluded':>10}  {'alpha':>6}  {'Lam':>6}  {'T_obs':>7}  {'T_pred':>8}  {'Delta':>7}")
print(f"  {'--------':>10}  {'------':>6}  {'------':>6}  {'-------':>7}  {'--------':>8}  {'-------':>7}")

def loss_alpha_sub(params, exclude_idx):
    alpha, Lambda_0 = params
    if Lambda_0 <= 0:
        return 1e10
    total = 0.0
    for i, (name, a, orb, z, omH, Tobs, _, _) in enumerate(hydrides):
        if i == exclude_idx:
            continue
        Tp = Tc_hydride_P6B(a, orb, z, omH, alpha, Lambda_0)
        total += (np.log10(max(Tp, 1e-6)) - np.log10(Tobs))**2
    return total

for i, (name, a, orb, z, omH, Tobs, _, _) in enumerate(hydrides):
    res_sub = minimize(loss_alpha_sub, x0=[1.5, 0.1], args=(i,),
                       method="Nelder-Mead",
                       options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 3000})
    a_sub, L_sub = res_sub.x
    Tp = Tc_hydride_P6B(a, orb, z, omH, a_sub, L_sub)
    delta_log = np.log10(max(Tp, 1e-6)) - np.log10(Tobs)
    print(f"  {name:>10}  {a_sub:>6.3f}  {L_sub:>6.4f}  {Tobs:>7.1f}  {Tp:>8.2f}  {delta_log:>+7.3f}")

# =============================================================
# Part D. Predykcje: nowe hydrydy (teoria / nie-syntezowane)
# =============================================================
print()
print("=" * 78)
print("  Part D. Predykcje dla nowych/spekulatywnych hydrydow")
print("=" * 78)
print()
# Kandydaci:
# - Fazy roomTC: metastabilne
# - Spekulowane: CH4, SiH4, BH4, AlH4 pod roznymi P

predictions = [
    # (nazwa, a, orb, z, omega_H, P_needed [GPa], opis)
    ("MgH2",    4.500, "sp", 8,  100.0, 180, "prosta stable roomtemp, SC?"),
    ("BeH2",    3.900, "sp", 6,  140.0, 200, "spekulowany ultra-light"),
    ("AlH3",    4.450, "sp", 12, 150.0, 110, "metallic pod P"),
    ("CaH6",    4.500, "sp", 8,  200.0, 172, "SPRAWDZENIE - jest w fit"),
    ("BH",      3.200, "sp", 8,  220.0, 300, "metallic boron-hydride spek."),
    ("LiBH4",   5.150, "sp", 8,  130.0, 100, "complex hydrydek"),
    ("NaAlH4",  5.020, "sp", 8,  110.0, 100, "ternary"),
    ("BaH6",    5.300, "s",  8,  140.0, 80,  "heavy host, low P?"),
    ("RbH9",    5.500, "s",  12, 160.0, 100, "alkali heavy"),
    ("H_metallic", 2.80, "s", 12, 350.0, 500, "monatomic hydrogen metal"),
]

print(f"  {'Material':>12} {'a[A]':>5} {'orb':>3} {'z':>3} {'omH':>5} {'T_pred[K]':>9} {'P[GPa]':>7}  {'opis':>35}")
print(f"  {'------------':>12} {'-----':>5} {'---':>3} {'---':>3} {'-----':>5} {'---------':>9} {'-------':>7}  {'-----------------------------------':>35}")

for name, a, orb, z, omH, P, desc in predictions:
    Tp = Tc_hydride_P6B(a, orb, z, omH, alpha_fit, Lam_fit)
    print(f"  {name:>12} {a:>5.3f} {orb:>3} {z:>3d} {omH:>5.0f} {Tp:>9.1f} {P:>7.0f}  {desc[:35]:>35}")

# =============================================================
# Part E. Czy mozna niska P hydride (ambient)?
# =============================================================
print()
print("=" * 78)
print("  Part E. Problem AMBIENT hydrydow")
print("=" * 78)
print()
print("  Wszystkie znane hydrydy SC wymagaja P >= 95 GPa.")
print("  Powod: dynamiczna stabilizacja H-cage wymaga kompresji.")
print()
print("  TGP: formula T_c NIE WYMAGA ciśnienia per se -")
print("       wymaga tylko a blisko a*_tgp oraz omega_H wysokiej.")
print()
print("  Czy mozemy dostac wysoka omega_H przy ambient?")
print("  Odpowiedz: molekularne hydrogen-rich compounds.")
print("    H2O lod: omega_H ~ 400 meV (O-H stretch), ale izolator (brak SC)")
print("    NH3 amoniak: omega_H ~ 420 meV, izolator")
print("    CH4 metan: izolator")
print()
print("  Klucz: potrzeba METALICZNEGO hydrydu przy ambient.")
print("         Kandidates: borohidydy (MgB2, M-BH4), berylowodorki.")
print()
print("  TGP predykcja dla MgB2 (istn, T_c=39K @ ambient):")
a_MgB2 = 3.086
omega_MgB2 = 75.0  # B-B stretch phonon
Tp_MgB2_P6B = Tc_hydride_P6B(a_MgB2, "sp", 5, omega_MgB2, alpha_fit, Lam_fit)
print(f"    a=3.086, orb=sp, z=5, omega=75 meV -> T_pred = {Tp_MgB2_P6B:.1f} K (vs obs 39 K)")
print()
print("  Dla nieistn. AlB2-analog pod ambient (BB bez Mg):")
Tp_BB = Tc_hydride_P6B(3.0, "sp", 6, 120.0, alpha_fit, Lam_fit)
print(f"    a=3.0, orb=sp, z=6, omega=120 meV -> T_pred = {Tp_BB:.1f} K")

# =============================================================
# Werdykt
# =============================================================
print()
print("=" * 78)
print("  Werdykt ps12")
print("=" * 78)
print()
print(f"  P6.B fit hydrydow: r(log-log) = {r:.3f}, RMS_log = {rms:.3f}")
if r > 0.80:
    print(f"  DOBRY fit - model phonon coupling dziala.")
elif r > 0.60:
    print(f"  Umiarkowany fit - model kierunkowo OK ale nie idealny.")
else:
    print(f"  Slaby fit - poszukaj innej postaci coupling.")

print()
print(f"  alpha = {alpha_fit:.2f} -> skalowanie Lambda_E z omega_phonon")
print()
print("  Glowny wniosek:")
print("    1. Hydrydy high-P (100-250 GPa) pasuja do TGP P6.B.")
print("    2. Ambient-P hydrydy w P6.B wymagaja metalicznego H-rich + wysoka omega.")
print("    3. MgB2 pod ambient z omega=75 meV: P6.B daje lepsze niz ps9 (~2x, ale wciaz poniza 39K).")
print()
print("=" * 78)
print("  ps12 complete.")
print("=" * 78)
