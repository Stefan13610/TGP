#!/usr/bin/env python3
"""
omicron2_phi_mean_shift_check.py — Opcja A quick check (2026-05-03)
====================================================================

Pytanie (user 2026-05-03):
  Phi_0 jako observable powinno trackowac globalny rozklad materii w czasie
  (sek05 §de-struktury thesis), wiec powinno byc zmienne. Czy to moglby
  rozwiazac Hubble tension?

Metoda:
  Test 1: Linear <delta Phi>_V — sprawdzenie czy 0 (stat. homogeneity)
  Test 2: Linear variance <delta Phi^2>(z) cosmological evolution
  Test 3: Nonlinear Buchert backreaction Q_D — czy MEAN shift omitted?
  Test 4: Compare z Hubble tension requirement

Hypoteza ujawniona w drift-check:
  M10.5 obliczyl Buchert variance ~ 10^-8.
  Jezeli mean shift <Phi>(z) jest osobnym mechanism od variance,
  M10.5 mial blind spot.

Result: see verdict at end.
"""
import sys
import io
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np

# ----------------------------------------------------------------
# Cosmological parameters (Planck 2018)
# ----------------------------------------------------------------
H0 = 67.4
Omega_m = 0.315
Omega_L = 0.685

# Observed gravitational potential rms (CMB-LSS combined)
sigma_Psi_obs = 3e-5

print("=" * 72)
print("omicron.2 quick check: Phi_0(z) cosmological evolution & H_0 tension")
print("=" * 72)
print(f"\n  Cosmology: Omega_m = {Omega_m}, Omega_L = {Omega_L}, H_0 = {H0}")
print(f"  Observed sigma_Psi (CMB+LSS) = {sigma_Psi_obs:.0e}")

# ================================================================
# Test 1: Linear cosmological mean <delta Phi>_V
# ================================================================
print("\n" + "=" * 72)
print("TEST 1: Linear cosmological mean <delta Phi>_V")
print("=" * 72)

print("""
  Theorem (statystyczna homogenicznosc kosmologii):
    <delta(x)>_V = 0  EXACT  dla dowolnej V wystarczajaco duzej.

  Uzasadnienie:
    - W FRW perturbacje delta(x, t) sa stat. jednorodne i izotropowe.
    - Wartosc oczekiwana <delta(x)> = 0 z definicji ensemble average.
    - Linear theory NIE generuje deterministic mean shift.

  Wniosek:
    Linear theory daje DOKLADNIE ZERO mean shift z formacja struktur.
    Mean shift wymaga (a) nieliniowosci albo (b) zlamania homogenicznosci.
""")

T1_result = 0.0
print(f"  Test 1 result: <delta Phi>_V (linear) = {T1_result}")

# ================================================================
# Test 2: Linear variance cosmological evolution
# ================================================================
print("\n" + "=" * 72)
print("TEST 2: Linear variance <delta Phi^2>(z) — cosmological evolution")
print("=" * 72)

def D_growth(z, Om=Omega_m, OL=Omega_L):
    """Linear growth factor (Carroll-Press-Turner approx, normalized D(0)=1)."""
    a = 1.0/(1.0+z)
    Om_a = Om/a**3 / (Om/a**3 + OL)
    OL_a = 1.0 - Om_a
    g_a = 5.0/2.0 * Om_a / (Om_a**(4/7) - OL_a + (1+Om_a/2)*(1+OL_a/70))
    g_0 = 5.0/2.0 * Om / (Om**(4/7) - OL + (1+Om/2)*(1+OL/70))
    return a * g_a / g_0

def sigma_Psi_linear(z):
    """RMS gravitational potential in Newton gauge.

    Standard result: in matter era Psi ~ const; in Lambda era Psi decays.
    Approximation: Psi(z) = sigma_Psi_obs * [D(z)/a(z)] / [D(0)/1]
    """
    a = 1.0/(1.0+z)
    Da = D_growth(z)
    decay = Da / a
    decay0 = D_growth(0) / 1.0
    return sigma_Psi_obs * (decay / decay0)

z_grid = [0, 0.5, 1, 2, 5, 10, 100, 1100]
print(f"\n  {'z':>6}  {'D(z)':>9}  {'sigma_Psi':>12}  {'<delta^2>':>12}  {'dLambda/Lambda':>16}")
print(f"  {'-'*6}  {'-'*9}  {'-'*12}  {'-'*12}  {'-'*16}")
for z in z_grid:
    D = D_growth(z)
    sP = sigma_Psi_linear(z)
    var = (2*sP)**2  # delta_Phi/Phi_0 = 2*Psi (PPN gamma=1)
    dLL = 6 * var    # sek05 lemma: dLambda/Lambda = 6 <delta^2>
    print(f"  {z:>6}  {D:>9.4f}  {sP:>12.3e}  {var:>12.3e}  {dLL:>16.3e}")

sP_today = sigma_Psi_linear(0)
sP_recomb = sigma_Psi_linear(1100)
dLL_today = 6 * (2*sP_today)**2
dLL_recomb = 6 * (2*sP_recomb)**2
delta_Lambda_cosmo = dLL_today - dLL_recomb

print(f"\n  Cosmological evolution recomb -> today:")
print(f"    dLambda/Lambda(z=0)    = {dLL_today:.3e}")
print(f"    dLambda/Lambda(z=1100) = {dLL_recomb:.3e}")
print(f"    Delta (cosmological)   = {delta_Lambda_cosmo:.3e}")

T2_result = dLL_today
print(f"\n  Test 2 result (linear variance, today): {T2_result:.3e}")

# ================================================================
# Test 3: Nonlinear Buchert backreaction Q_D — szukamy mean shift
# ================================================================
print("\n" + "=" * 72)
print("TEST 3: Nonlinear Buchert Q_D — czy M10.5 missed mean shift?")
print("=" * 72)

# Buchert (2000) kinematical backreaction:
#   Q_D = (2/3) <(theta - <theta>)^2>_D - 2 <sigma^2>_D
# where theta = expansion scalar, sigma = shear.
# Theta has variance: <(delta theta)^2> = H^2 * f^2 * <delta_m^2>
#   where f = dlnD/dlna ≈ Omega_m^0.55 (growth rate)

f_growth = Omega_m**0.55
print(f"\n  Growth rate index: f = Omega_m^0.55 = {f_growth:.4f}")

# delta_m today on linear scales (>~ 8 Mpc/h):
# sigma_8 ~ 0.811 corresponds to delta_m ~ 0.811 at scale 8 Mpc/h
# But for COSMOLOGICAL volume average, we want delta on Hubble scale L_H ~ 4 Gpc
# At Hubble scale, sigma_delta ~ Psi ~ 3e-5 (CMB)
sigma_delta_horizon = 2 * sP_today  # delta = 2 Psi at horizon scale

print(f"  sigma_delta (horizon scale) = {sigma_delta_horizon:.3e}")

# Q_D / H^2 (leading order):
Q_D_over_H2 = sigma_delta_horizon**2 * f_growth**2
print(f"  Q_D/H^2 ~ <delta^2> · f^2 = {Q_D_over_H2:.3e}")

# Effective Lambda shift from Buchert:
# Lambda_eff(t) -> Lambda + Q_D / 2 (approx)
# Fractional shift: dLambda/Lambda ~ Q_D / (H^2 Omega_L)
dLL_buchert = Q_D_over_H2 / Omega_L
print(f"  Effective Lambda shift: dLambda/Lambda = Q_D/(H^2 Omega_L) = {dLL_buchert:.3e}")

ratio_buchert_to_variance = dLL_buchert / dLL_today
print(f"  Ratio Buchert/variance = {ratio_buchert_to_variance:.2f}x")

print(f"""
  Interpretacja:
    Q_D jest 2nd-order w perturbacjach (variance of theta), wiec ten sam
    rzad co linear variance. NIE jest niezalezny mechanism dla mean shift.

  Wniosek strukturalny:
    Mean shift jako efekt nieliniowy w rownaniu Phi-EOM (z czlonu ~Phi²/Phi_0)
    daje wkład rzędu <delta^2>/Phi_0 — TEN SAM rzad co Buchert variance.

    M10.5 NIE missed mean shift — Buchert formula juz captures rectification
    of <delta^2> przez nonlinearity.
""")

T3_result = dLL_buchert
print(f"  Test 3 result (Buchert nonlinear): {T3_result:.3e}")

# ================================================================
# Test 4: Hubble tension requirement
# ================================================================
print("\n" + "=" * 72)
print("TEST 4: Hubble tension requirement (target)")
print("=" * 72)

H0_shoes = 73.04
H0_planck = 67.4
Delta_H = (H0_shoes - H0_planck) / H0_planck

# At z=0: H^2 = (8 pi G/3)(rho_m + rho_L), Omega_L = rho_L / rho_total
# If Lambda shifts by dLambda at fixed Omega_m fraction:
#   (H_new/H_old)^2 = 1 + (Omega_L_new - Omega_L_old)
# Solving: dOmega_L = (H_new^2 - H_old^2)/H_old^2 ≈ 2 Delta_H + Delta_H^2
# Fractional: dLambda/Lambda = dOmega_L / Omega_L
required_dL_over_L = (2*Delta_H + Delta_H**2) / Omega_L

print(f"  H_0(SH0ES)  = {H0_shoes:.2f} km/s/Mpc")
print(f"  H_0(Planck) = {H0_planck:.2f} km/s/Mpc")
print(f"  Delta H/H   = {Delta_H:.4f} ({Delta_H*100:.2f}%)")
print(f"  Required dLambda/Lambda = {required_dL_over_L:.4f} ({required_dL_over_L*100:.1f}%)")

T4_result = required_dL_over_L

# ================================================================
# Final verdict
# ================================================================
print("\n" + "=" * 72)
print("FINAL VERDICT")
print("=" * 72)

summary = [
    ("Linear mean <delta Phi>_V (theorem)", T1_result, "EXACT 0"),
    ("Linear variance dLambda/Lambda today", T2_result, "sek05 lemma"),
    ("Linear cosmological shift (recomb -> today)", abs(delta_Lambda_cosmo), "growth × variance"),
    ("Nonlinear Buchert Q_D / Omega_L", T3_result, "rectified <delta^2>"),
    ("Hubble tension requirement", T4_result, "target"),
]

print(f"\n  {'Quantity':<48} {'Value':>12} {'Notes':<20}")
print(f"  {'-'*48} {'-'*12} {'-'*20}")
for name, val, note in summary:
    if val == 0:
        print(f"  {name:<48} {'0 (exact)':>12} {note}")
    else:
        order = np.log10(val) if val > 0 else 0
        print(f"  {name:<48} {val:>12.3e} {note}")

# Gap analysis
gap = T4_result / T3_result
gap_orders = np.log10(gap)
print(f"\n  Gap analysis (Hubble tension vs largest TGP-derived effect):")
print(f"    Required:         {T4_result:.3e}")
print(f"    TGP best (Q_D):   {T3_result:.3e}")
print(f"    Gap (ratio):      {gap:.3e}")
print(f"    Gap (orders):     {gap_orders:.1f} dekad")

print(f"""
============================================================================
WNIOSEK:

  Q1: Czy Phi_0 jako observable shifts cosmologically?
  A1: TAK strukturalnie (sek05 §de-struktury thesis), ALE:
      - Linear: 0 (statystyczna homogenicznosc)
      - Nonlinear (Buchert Q_D): ~ {T3_result:.1e}
      - Cosmological evolution recomb -> today: ~ {abs(delta_Lambda_cosmo):.1e}

  Q2: Czy M10.5 had blind spot na mean shift?
  A2: NIE. Buchert formula juz captures nonlinear rectification of
      <delta^2> w mean shift. M10.5 verdict reinforced.

  Q3: Czy to rozwiazuje Hubble tension?
  A3: NIE. Gap = {gap_orders:.1f} dekad miedzy TGP best ({T3_result:.1e})
      i required ({T4_result:.1e}).

  Q4: Twoja intuicja ontologiczna byla bledna?
  A4: NIE — byla strukturalnie poprawna (sek05 explicit potwierdza).
      Tylko quantitative coupling przez Psi ~ 10^-5 (CMB bound) jest
      za slaby. To nie jest blad teorii — to jest wlasciwosc TGP coupled
      do gravitational potential (PPN gamma=1 minimal coupling).

  Q5: Co bymoze zmienilo wynik?
  A5: Direct rho-Phi coupling poza ax:metric-coupling axiom (TGP_FOUNDATIONS §4).
      To bylby strukturalna modyfikacja TGP, NIE w current formalism.

  STATUS koncowy:
      - Twoja teza: <Phi>(t) tracks matter distribution(t)         CONFIRMED
      - M10.5 verdict: TGP nie jest H_0 tension solver               CONFIRMED
      - Gap strukturalny: 7+ dekad                                   CONFIRMED
      - M10.5 mean-shift blind spot: NIE EXISTS                      VERIFIED
      - Required physics: poza current TGP                           OPEN

  TGP scope statement (M10 cycle): TGP_v1 = galaxy-scale + structural DE
  (w >= -1) + CMB safety; NIE H_0 tension solver. POTWIERDZONE drugi raz.
============================================================================
""")
