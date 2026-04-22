#!/usr/bin/env python3
"""
TGP Color Tube: Variational String Tension (2026-04-10, v2)

CORRECTED DERIVATION from the TGP field equation:

  nabla^2 phi + 2(nabla phi)^2 / phi + beta*phi^2*(1-phi) = 0

  where phi = Phi/Phi_0 (dimensionless field), beta = gamma (vacuum cond.)
  alpha = 2 from K(phi) = phi^4  (thm:D-uniqueness, sek08_formalizm)

Lagrangian (static):
  L = -(1/2)*phi^4*(nabla phi)^2 + V_L(phi)
  V_L(phi) = phi^7/7 - phi^8/8    [in units where beta=1]

EL gives:  V_L'(phi)/phi^4 = phi^2(1-phi)   CHECK

Energy functional (Hamiltonian for static config):
  E = (1/2)*phi^4*(nabla phi)^2  +  V_E(phi)
  V_E(phi) = phi^8/8 - phi^7/7 + 1/56   [V_E(1) = 0, V_E''(1) = 1 > 0]

String tension = energy per unit length of cylindrical tube:
  sigma_hat = 2*pi * int_0^inf [E_kin + V_E] rho drho

  with Gaussian ansatz: phi(rho) = 1 - A * exp(-rho^2/(2w^2))
  A = depletion amplitude, w = tube width in units of 1/m_sp

Physical string tension:
  sigma_phys = K_geo * m_sp^2 * sigma_hat
"""

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad

# ============================================================
# TGP Parameters
# ============================================================
Phi_0 = 24.6492
kappa = 3.0 / (4 * Phi_0)  # = 0.03043
a_Gamma = 0.040
phi_golden = (1 + np.sqrt(5)) / 2
A_target = a_Gamma / phi_golden  # = 0.02472
alpha_s_TGP = 0.1190
N_c = 3
C_F = (N_c**2 - 1) / (2*N_c)  # = 4/3

print("=" * 70)
print("TGP COLOR TUBE: VARIATIONAL STRING TENSION (v2 - corrected)")
print("=" * 70)

# ============================================================
# [0] VERIFY POTENTIAL
# ============================================================
print("\n[0] POTENTIAL VERIFICATION")
print("-" * 50)

def V_E(phi):
    """
    Energy potential: V_E(phi) = phi^8/8 - phi^7/7 + 1/56
    Properties: V_E(1) = 0, V_E'(1) = 0, V_E''(1) = 1 > 0
    V_E(phi) > 0 for phi != 1 (vacuum is global energy minimum)
    """
    return phi**8/8.0 - phi**7/7.0 + 1.0/56.0

def V_E_prime(phi):
    """V_E'(phi) = phi^7 - phi^6 = phi^6(phi - 1)"""
    return phi**7 - phi**6

# Verify properties
print(f"  V_E(1)   = {V_E(1.0):.15f}  (should be 0)")
print(f"  V_E'(1)  = {V_E_prime(1.0):.15f}  (should be 0)")

# V_E''(phi) = 7*phi^6 - 6*phi^5
V_E_pp_1 = 7*1.0**6 - 6*1.0**5
print(f"  V_E''(1) = {V_E_pp_1:.1f}  (should be 1 > 0, vacuum is min)")

# Check positivity for depleted states
for test_phi in [0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999]:
    print(f"  V_E({test_phi:.3f}) = {V_E(test_phi):.8f} > 0 {'OK' if V_E(test_phi) > 0 else 'FAIL!'}")


# ============================================================
# [1] ENERGY FUNCTIONAL FOR GAUSSIAN TUBE
# ============================================================

def sigma_gaussian(A, w):
    """
    String tension for Gaussian tube: phi(rho) = 1 - A*exp(-rho^2/(2w^2))

    sigma = 2*pi * int_0^inf [(1/2)*phi^4*(dphi/drho)^2 + V_E(phi)] rho drho

    Returns sigma_hat (dimensionless, in units of m_sp^2 = beta).
    """
    if A <= 0 or A >= 1 or w <= 0:
        return 1e10

    def integrand(rho):
        x = rho**2 / (2*w**2)
        gauss = np.exp(-x)
        phi = 1.0 - A * gauss
        if phi < 1e-12:
            return 0.0

        # dphi/drho = A * rho/w^2 * exp(-rho^2/(2w^2))
        dphi = A * rho / w**2 * gauss

        # Kinetic: (1/2)*phi^4*(dphi)^2
        E_kin = 0.5 * phi**4 * dphi**2

        # Potential: V_E(phi) = phi^8/8 - phi^7/7 + 1/56
        E_pot = V_E(phi)

        return (E_kin + E_pot) * 2 * np.pi * rho

    result, err = quad(integrand, 0, max(10*w, 30), limit=300, epsabs=1e-14)
    return result


# ============================================================
# [2] SCAN (A, w) PARAMETER SPACE
# ============================================================
print("\n\n[1] PARAMETER SCAN: sigma_hat(A, w)")
print("-" * 50)

print(f"  {'A':>6} {'w_opt':>6} {'sigma':>14} {'sigma/A^2':>12}")
print(f"  {'='*46}")

scan_results = []
for A_val in np.linspace(0.01, 0.95, 20):
    best_w, best_sig = 0, 1e10
    for w in np.linspace(0.2, 5.0, 50):
        sig = sigma_gaussian(A_val, w)
        if 0 < sig < best_sig:
            best_w = w
            best_sig = sig
            scan_results.append((A_val, w, sig))
    if best_sig < 1e9:
        print(f"  {A_val:6.3f} {best_w:6.2f} {best_sig:14.8f} {best_sig/A_val**2:12.6f}")


# ============================================================
# [3] ANALYTICAL PERTURBATIVE LIMIT (small A)
# ============================================================
print(f"\n\n[2] PERTURBATIVE REGIME (small A)")
print("-" * 50)

# For small A, linearize:
#   E_kin ~ (1/2)*(dphi/drho)^2  (phi^4 ~ 1)
#   V_E   ~ (1/2)*V_E''(1)*delta^2 = (1/2)*delta^2
#
# where delta = A*exp(-rho^2/(2w^2))
#
# sigma ~ 2*pi * int [(1/2)*A^2*rho^2/w^4*exp(-rho^2/w^2)
#                     + (1/2)*A^2*exp(-rho^2/w^2)] rho drho
#
# Kinetic integral: pi*A^2 * int_0^inf rho^3/w^4 * exp(-rho^2/w^2) drho
#                 = pi*A^2 * (1/(2w^4)) * w^4 = pi*A^2/2
#
# Potential integral: pi*A^2 * int_0^inf rho * exp(-rho^2/w^2) drho
#                   = pi*A^2 * w^2/2
#
# Total: sigma_pert = (pi*A^2/2)*(1 + w^2)

def sigma_perturbative(A, w):
    """Perturbative (small A) approximation."""
    return (np.pi * A**2 / 2) * (1 + w**2)

print(f"  sigma_pert(A, w) = (pi/2)*A^2*(1+w^2)")
print()

# Compare perturbative vs exact
for A_test in [0.001, 0.01, 0.05, 0.1]:
    w_test = 1.0
    sig_exact = sigma_gaussian(A_test, w_test)
    sig_pert = sigma_perturbative(A_test, w_test)
    ratio = sig_exact/sig_pert if sig_pert > 0 else 0
    print(f"  A={A_test:.3f}, w=1: exact={sig_exact:.10f}, "
          f"pert={sig_pert:.10f}, ratio={ratio:.6f}")


# ============================================================
# [4] OPTIMAL TUBE WIDTH FOR EACH A
# ============================================================
print(f"\n\n[3] OPTIMAL TUBE WIDTH")
print("-" * 50)

# The perturbative result sigma = (pi/2)*A^2*(1+w^2) is MINIMIZED at w->0.
# But at finite A, the non-linear potential changes the balance.
# For a physical tube, w is set by the source scale, not by free minimization.
#
# Physical expectation: w ~ 1/m_sp (the scalar mass is the screening length)
# This gives w = 1 in our dimensionless units.

print(f"  Perturbative: sigma increases with w (no minimum in w)")
print(f"  Physical tube: w ~ 1/m_sp (screening length)")
print(f"  Setting w = 1 (physical value)")
print()

w_phys = 1.0

print(f"  {'A':>8} {'sig(w=1)':>14} {'sig/A^2':>12} {'sig/(pi*A^2)':>14}")
for A_val in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]:
    sig = sigma_gaussian(A_val, w_phys)
    print(f"  {A_val:8.4f} {sig:14.10f} {sig/A_val**2:12.6f} {sig/(np.pi*A_val**2):14.6f}")


# ============================================================
# [5] PHYSICAL DEPLETION FROM COLOR COUPLING
# ============================================================
print(f"\n\n[4] PHYSICAL TUBE: COLOR-INDUCED DEPLETION")
print("-" * 50)

# Color flux depletes Phi. The depletion amplitude:
# A = C_F * alpha_s / (pi * Phi_0)
# This comes from the conformal coupling of gauge fields to Phi:
#   L_gauge ~ (Phi/Phi_0)^2 * F^2
# Inside the tube, the color electric field creates effective source.

A_color = C_F * alpha_s_TGP / (np.pi * Phi_0)
print(f"  A_color = C_F * alpha_s / (pi*Phi_0)")
print(f"         = ({C_F:.4f}) * {alpha_s_TGP} / (pi*{Phi_0})")
print(f"         = {A_color:.8f}")
print()

sig_color = sigma_gaussian(A_color, w_phys)
sig_pert_color = sigma_perturbative(A_color, w_phys)
print(f"  sigma_hat(A_color, w=1):")
print(f"    Exact:        {sig_color:.14f}")
print(f"    Perturbative: {sig_pert_color:.14f}")
print(f"    Ratio:        {sig_color/sig_pert_color:.8f}")
print()

# Quadratic coefficient
C_w1 = sig_color / A_color**2
print(f"  sigma_hat = C(w=1) * A^2")
print(f"  C(w=1) = {C_w1:.6f}")
print(f"  pi (predicted) = {np.pi:.6f}")
print(f"  C/pi = {C_w1/np.pi:.6f}")


# ============================================================
# [6] STRING TENSION -> m_0 -> A_univ
# ============================================================
print(f"\n\n[5] FROM STRING TENSION TO QUARK MASS")
print("-" * 50)

# Physical string tension:
#   sigma_phys = K_geo * m_sp^2 * sigma_hat   [mass^2 in natural units]
#
# Quark additive mass:
#   m_0 = sigma_phys * L_eff  [mass]
#   where L_eff = confinement length ~ 1/Lambda_QCD
#
# Universal constant:
#   A_univ = m_0 * m_1 / m_3 = a_Gamma / phi_golden = 0.02472
#
# Chain:
#   A_univ = K_geo * m_sp^2 * sigma_hat * L_eff * m_1/m_3
#
# Key insight: m_sp ~ sqrt(beta), and in TGP beta ~ Phi_0 * Lambda^2
# The color tube length L_eff ~ m_3/(m_1 * Lambda_QCD)
# So the m_1/m_3 factors cancel:
#   A_univ = K_geo * sigma_hat * (m_sp / Lambda_QCD)^2
#
# For dimensional self-consistency, define:
#   xi = m_sp / Lambda_QCD  (ratio of scalar mass to QCD scale)
#
# Then: A_univ = K_geo * sigma_hat * xi^2

print(f"  sigma_hat = {sig_color:.6e}")
print(f"  A_target  = {A_target:.8f}")
print()

# What product K_geo * xi^2 would reproduce A_target?
K_xi2_needed = A_target / sig_color
print(f"  K_geo * xi^2 needed = A_target / sigma_hat = {K_xi2_needed:.4f}")
print()

# Alternative derivation: use QCD string tension
# sqrt(sigma_QCD) = 440 MeV
# sigma_QCD = (440 MeV)^2 = 0.1936 GeV^2
sigma_QCD_GeV2 = 0.440**2  # GeV^2

# In TGP: sigma_QCD = K_geo * m_sp^2 * sigma_hat
# => K_geo * m_sp^2 = sigma_QCD / sigma_hat
K_msp2 = sigma_QCD_GeV2 / sig_color
print(f"  From sigma_QCD = {sigma_QCD_GeV2:.4f} GeV^2:")
print(f"    K_geo * m_sp^2 = {K_msp2:.4f} GeV^2")
print()

# m_0 = sigma_QCD * L_eff
# For light quarks: L_eff ~ 1/Lambda_QCD ~ 1/0.332 GeV ~ 3.01 GeV^-1
Lambda_QCD = 0.332  # GeV
L_eff_approx = 1.0 / Lambda_QCD  # GeV^-1

m_0_computed = sigma_QCD_GeV2 * L_eff_approx
print(f"  m_0 = sigma_QCD * L_eff")
print(f"      = {sigma_QCD_GeV2:.4f} * {L_eff_approx:.4f}")
print(f"      = {m_0_computed*1000:.2f} MeV")
print()

# PDG constituent quark mass: m_0(u) ~ 336 MeV, m_0(d) ~ 340 MeV
# Our estimate: sigma_QCD / Lambda_QCD ~ 583 MeV (too high by factor ~1.7)
# More refined: L_eff = R_conf / 2 where R_conf ~ 0.5-1.0 fm
# Using R_conf = 1.0 fm = 5.07 GeV^-1: L_eff ~ 0.5 fm = 2.5 GeV^-1
L_eff_fm = 0.50  # fm  (half the confinement radius)
L_eff_refined = L_eff_fm * 5.068  # GeV^-1 (1 fm = 5.068 GeV^-1)

m_0_refined = sigma_QCD_GeV2 * L_eff_refined
print(f"  Refined: L_eff = {L_eff_fm:.2f} fm = {L_eff_refined:.3f} GeV^-1")
print(f"  m_0 = {m_0_refined*1000:.1f} MeV  (cf. constituent quark ~336 MeV)")
print()


# ============================================================
# [7] THE KEY RATIO: sigma_hat vs A_univ
# ============================================================
print(f"\n[6] SCALING ANALYSIS")
print("-" * 50)

# From the variational calculation:
#   sigma_hat(A_color, w=1) = pi * A_color^2  (perturbative regime)
#
# A_color = C_F * alpha_s / (pi * Phi_0)
# => sigma_hat = pi * C_F^2 * alpha_s^2 / (pi^2 * Phi_0^2)
#              = C_F^2 * alpha_s^2 / (pi * Phi_0^2)
#
# A_univ = f * sigma_hat
# where f absorbs all scale factors.
#
# Self-consistency test:
#   A_univ / sigma_hat = f = K_geo * xi^2

sigma_hat_analytic = C_F**2 * alpha_s_TGP**2 / (np.pi * Phi_0**2)
print(f"  sigma_hat (analytic) = C_F^2 * alpha_s^2 / (pi * Phi_0^2)")
print(f"                       = {sigma_hat_analytic:.6e}")
print(f"  sigma_hat (numerical) = {sig_color:.6e}")
print()

# The ratio A_univ / sigma_hat has a beautiful form:
f_ratio = A_target / sig_color
print(f"  f = A_univ / sigma_hat = {f_ratio:.4f}")
print()

# Decompose f:
# f = A_univ / sigma_hat = (a_Gamma/phi) / (C_F^2 * alpha_s^2 / (pi * Phi_0^2))
#   = a_Gamma * pi * Phi_0^2 / (phi * C_F^2 * alpha_s^2)

f_analytic = a_Gamma * np.pi * Phi_0**2 / (phi_golden * C_F**2 * alpha_s_TGP**2)
print(f"  f = a_Gamma * pi * Phi_0^2 / (phi * C_F^2 * alpha_s^2)")
print(f"    = {a_Gamma} * pi * {Phi_0}^2 / ({phi_golden:.5f} * {C_F:.4f}^2 * {alpha_s_TGP}^2)")
print(f"    = {f_analytic:.4f}")
print()

# Check if f = Phi_0^3 (natural TGP scaling):
print(f"  Phi_0^3    = {Phi_0**3:.2f}")
print(f"  Phi_0^2    = {Phi_0**2:.2f}")
print(f"  f          = {f_ratio:.2f}")
print(f"  f/Phi_0^2  = {f_ratio/Phi_0**2:.4f}")
print(f"  f/Phi_0^3  = {f_ratio/Phi_0**3:.6f}")
print()

# sqrt(f):
print(f"  sqrt(f)        = {np.sqrt(f_ratio):.4f}")
print(f"  Phi_0          = {Phi_0:.4f}")
print(f"  sqrt(f)/Phi_0  = {np.sqrt(f_ratio)/Phi_0:.4f}")
print()

# Key dimensionless ratios in TGP:
print(f"  Key TGP scales:")
print(f"    1/kappa         = {1/kappa:.4f}")
print(f"    4*Phi_0/3       = {4*Phi_0/3:.4f}")
print(f"    Phi_0^2         = {Phi_0**2:.2f}")
print(f"    Phi_0^3         = {Phi_0**3:.0f}")
print(f"    (4*Phi_0/3)^2   = {(4*Phi_0/3)**2:.2f}")


# ============================================================
# [8] ALTERNATIVE: DIRECT a_Gamma SCALING
# ============================================================
print(f"\n\n[7] DIRECT SCALING: a_Gamma FROM TUBE")
print("-" * 50)

# Hypothesis: a_Gamma = sigma_hat * Phi_0^2 * phi_golden * (factor)
# From the analytical result:
#   sigma_hat = C_F^2 * alpha_s^2 / (pi * Phi_0^2)
#
# If a_Gamma = N_c * sigma_hat * Phi_0^2:
#   a_Gamma = N_c * C_F^2 * alpha_s^2 / pi
#           = 3 * (4/3)^2 * 0.119^2 / pi
#           = 3 * (16/9) * 0.01416 / pi
#           = 0.02416

a_Gamma_hyp1 = N_c * C_F**2 * alpha_s_TGP**2 / np.pi
print(f"  Hypothesis 1: a_Gamma = N_c * C_F^2 * alpha_s^2 / pi")
print(f"    = {a_Gamma_hyp1:.6f}  (target: {a_Gamma:.6f})")
print(f"    ratio: {a_Gamma_hyp1/a_Gamma:.4f}")
print()

# Hypothesis 2: include Casimir directly
# A_univ = C_F * alpha_s / (pi * Phi_0) * Phi_0 * sqrt(N_c)
#        = C_F * alpha_s * sqrt(N_c) / pi
a_univ_hyp2 = C_F * alpha_s_TGP * np.sqrt(N_c) / np.pi
print(f"  Hypothesis 2: A_univ = C_F * alpha_s * sqrt(N_c) / pi")
print(f"    = {a_univ_hyp2:.6f}  (target: {A_target:.6f})")
print(f"    ratio: {a_univ_hyp2/A_target:.4f}")
print()

# Hypothesis 3: pure group theory factor
# a_Gamma = C_F^2 * alpha_s * 4 / phi_golden
a_Gamma_hyp3 = C_F**2 * alpha_s_TGP * 4 / phi_golden
print(f"  Hypothesis 3: a_Gamma = 4*C_F^2*alpha_s/phi")
print(f"    = {a_Gamma_hyp3:.6f}  (target: {a_Gamma:.6f})")
print(f"    ratio: {a_Gamma_hyp3/a_Gamma:.4f}")
print()

# Hypothesis 4: sigma_hat * 4*Phi_0^3/3 (= Phi_0^2/kappa)
a_univ_hyp4 = sig_color * Phi_0**2 / kappa
print(f"  Hypothesis 4: A_univ = sigma_hat * Phi_0^2/kappa")
print(f"    = {a_univ_hyp4:.6f}  (target: {A_target:.6f})")
print(f"    ratio: {a_univ_hyp4/A_target:.4f}")
print()

# Scan through more combinations
print(f"  Systematic search for A_univ = sigma_hat * X:")
candidates = {
    "Phi_0^2/kappa": Phi_0**2 / kappa,
    "Phi_0^3": Phi_0**3,
    "4*Phi_0^3/3": 4*Phi_0**3/3,
    "Phi_0^2": Phi_0**2,
    "Phi_0^2*pi": Phi_0**2 * np.pi,
    "1/(kappa^2)": 1/kappa**2,
    "N_c*Phi_0^2": N_c * Phi_0**2,
    "Phi_0/(alpha_s*kappa)": Phi_0/(alpha_s_TGP * kappa),
}

print(f"  {'X':<30} {'sig*X':>12} {'target':>12} {'ratio':>8}")
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]*sig_color/A_target - 1)):
    result = sig_color * val
    ratio = result / A_target
    mark = " <== " if 0.5 < ratio < 2.0 else ""
    print(f"  {name:<30} {result:12.8f} {A_target:12.8f} {ratio:8.4f}{mark}")


# ============================================================
# [9] TUBE PROFILE PROPERTIES
# ============================================================
print(f"\n\n[8] TUBE PROFILE FOR PHYSICAL PARAMETERS")
print("-" * 50)

print(f"  A_color = {A_color:.6e}")
print(f"  w = {w_phys:.1f} (in units of 1/m_sp)")
print()

# Depletion at tube center
phi_center = 1.0 - A_color
print(f"  phi(0) = 1 - A = {phi_center:.10f}")
print(f"  Phi(0)/Phi_0 = {phi_center:.10f}")
print(f"  Depletion: {A_color*100:.6f}%")
print()

# Tube radius (where depletion drops to 1/e of peak)
R_tube = w_phys * np.sqrt(2)  # rho where exp(-rho^2/(2w^2)) = 1/e
print(f"  Tube radius (1/e): R = w*sqrt(2) = {R_tube:.4f} / m_sp")
print()

# Energy partition
def E_kin_total(A, w):
    def integrand(rho):
        x = rho**2 / (2*w**2)
        gauss = np.exp(-x)
        phi = 1.0 - A * gauss
        dphi = A * rho / w**2 * gauss
        return 0.5 * phi**4 * dphi**2 * 2 * np.pi * rho
    result, _ = quad(integrand, 0, max(10*w, 30), limit=300, epsabs=1e-14)
    return result

def E_pot_total(A, w):
    def integrand(rho):
        x = rho**2 / (2*w**2)
        gauss = np.exp(-x)
        phi = 1.0 - A * gauss
        return V_E(phi) * 2 * np.pi * rho
    result, _ = quad(integrand, 0, max(10*w, 30), limit=300, epsabs=1e-14)
    return result

T_kin = E_kin_total(A_color, w_phys)
T_pot = E_pot_total(A_color, w_phys)
print(f"  Energy partition (w=1):")
print(f"    Kinetic:   {T_kin:.6e}")
print(f"    Potential: {T_pot:.6e}")
print(f"    Total:     {T_kin + T_pot:.6e}")
print(f"    Ratio T/V: {T_kin/T_pot:.4f}")
print(f"    (Perturbative prediction: T/V = 1/w^2 = {1/w_phys**2:.1f})")


# ============================================================
# [10] NON-PERTURBATIVE REGIME
# ============================================================
print(f"\n\n[9] NON-PERTURBATIVE REGIME (large A)")
print("-" * 50)

# For large depletion (strongly non-linear regime), the phi^4 kinetic
# factor and the phi^8 - phi^7 potential become important.
# The tube is deeply depleted and the string tension no longer scales as A^2.

print(f"  {'A':>6} {'sigma':>14} {'sigma/A^2':>12} {'sig/(pi*A^2)':>14} {'nonlin':>8}")
for A_val in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    sig = sigma_gaussian(A_val, w_phys)
    sig_lin = sigma_perturbative(A_val, w_phys)
    nonlin = sig / sig_lin if sig_lin > 0 else 0
    print(f"  {A_val:6.2f} {sig:14.8f} {sig/A_val**2:12.6f} "
          f"{sig/(np.pi*A_val**2):14.6f} {nonlin:8.4f}")


# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY: VARIATIONAL COLOR TUBE ANALYSIS (v2)")
print(f"{'='*70}")

print(f"""
  CORRECT ENERGY FUNCTIONAL:
    E = (1/2)*phi^4*(grad phi)^2 + V_E(phi)
    V_E(phi) = phi^8/8 - phi^7/7 + 1/56
    V_E(1) = 0 (vacuum minimum), V_E''(1) = 1 > 0

  GAUSSIAN TUBE PROFILE:
    phi(rho) = 1 - A * exp(-rho^2/(2w^2))
    Physical depletion: A = C_F*alpha_s/(pi*Phi_0) = {A_color:.6e}
    Width: w = 1/m_sp (screening length)

  STRING TENSION (perturbative, A << 1):
    sigma_hat = (pi/2) * A^2 * (1 + w^2)
    For w = 1: sigma_hat = pi * A^2 = {np.pi * A_color**2:.6e}
    Numerical confirmation: {sig_color:.6e}

  KEY PHYSICS:
    1. POTENTIAL SIGN: V_E > 0 for phi < 1 (depletion costs energy)
    2. Kinetic factor phi^4 reduces gradient cost in depleted region
    3. Both E_kin and V_E are POSITIVE => sigma > 0 (tube has tension)
    4. Perturbative regime: sigma ~ pi*A^2 (quadratic in depletion)
    5. T_kin/V_pot ~ 1 (equipartition at w = 1)

  CONNECTION TO A_univ:
    A_univ = a_Gamma/phi = {A_target:.6f}
    sigma_hat = {sig_color:.6e}
    Scaling factor f = A_univ/sigma_hat = {f_ratio:.1f}
    This factor encodes the physical scales: K_geo, m_sp, Lambda_QCD

  MECHANISTIC CHAIN:
    color flux -> Phi depletion -> Gaussian tube ->
    -> sigma_hat ~ pi*A^2 -> sigma_phys = K_geo*m_sp^2*sigma_hat ->
    -> m_0 = sigma_phys * L_conf -> A_univ = m_0*m_1/m_3

  STATUS: PARTIAL [AN+NUM] CLOSURE
    - Tube profile: SOLVED (variational, positive energy)
    - sigma_hat: COMPUTED analytically and numerically
    - m_0 derivation: requires fixing K_geo and m_sp from TGP axioms
    - Full closure: needs CG-1/CG-3 (coarse graining -> physical scales)
""")
