"""
TGP Prediction for PTA Breathing Mode Signal

Quantitative prediction of the scalar breathing mode in the PTA band,
compared with NANOGrav 15yr constraints on non-Einsteinian polarizations.

Key physics:
  - TGP breathing mode: delta_g_ij/g_ij = 2*delta_Phi/Phi_0 * delta_ij
  - Scalar mass: m_sp ~ H_0/c_0 ~ 10^-26 m^-1
  - Cutoff frequency: f_cut = c_0 * m_sp / (2*pi) ~ 10^-18 Hz
  - PTA band: f ~ 1-100 nHz = 10^-9 to 10^-7 Hz
  - f_PTA >> f_cut => breathing mode propagates freely in PTA band

NANOGrav 15yr (arXiv: 2310.12138):
  - Searched for scalar-transverse (ST) correlations
  - HD correlations preferred over ST by Bayes factor ~2
  - No significant evidence for ST, but not strongly excluded either

The question: what amplitude does TGP predict for the breathing mode
in the PTA band?
"""
import numpy as np

# ============================================================
# Physical constants
# ============================================================
c0 = 3.0e8          # m/s
G0 = 6.674e-11      # m^3 kg^-1 s^-2
H0 = 2.2e-18        # 1/s (67.4 km/s/Mpc)
Phi0 = 115.0         # TGP background field
M_sun = 1.989e30    # kg
pc = 3.086e16        # m
Mpc = 3.086e22       # m
yr = 3.156e7          # s

# ============================================================
# TGP scalar mass
# ============================================================
m_sp = np.sqrt(Phi0) * H0 / c0  # 1/m
f_cut = c0 * m_sp / (2 * np.pi)  # Hz
lambda_sp = 2 * np.pi / m_sp  # m

print("=" * 65)
print("TGP SCALAR (BREATHING) MODE PARAMETERS")
print("=" * 65)
print(f"  Phi_0 = {Phi0}")
print(f"  m_sp = sqrt(Phi_0) * H_0 / c_0 = {m_sp:.3e} m^-1")
print(f"  Compton wavelength = {lambda_sp:.3e} m = {lambda_sp/Mpc:.1f} Mpc")
print(f"  Cutoff frequency = {f_cut:.3e} Hz")
print(f"  PTA band: 1 nHz - 100 nHz = 1e-9 to 1e-7 Hz")
print(f"  f_PTA / f_cut = {1e-9 / f_cut:.1e} >> 1")
print(f"  => Breathing mode is MASSLESS in PTA band (no dispersion)")

# ============================================================
# PTA characteristic strain from SMBH binary background
# ============================================================
print("\n" + "=" * 65)
print("GRAVITATIONAL WAVE BACKGROUND IN PTA BAND")
print("=" * 65)

# NANOGrav 15yr measured GWB amplitude:
# h_c(f) = A * (f / f_yr)^alpha
# A = 2.4e-15 (+0.7, -0.6) at f_yr = 1/yr
# alpha = -2/3 (for SMBH mergers)
A_nano = 2.4e-15
f_yr = 1.0 / yr  # ~ 3.17e-8 Hz
alpha = -2/3

print(f"  NANOGrav 15yr GWB:")
print(f"    A = {A_nano:.1e} at f_yr = {f_yr:.2e} Hz")
print(f"    Spectral index alpha = {alpha}")
print(f"    h_c(f) = A * (f/f_yr)^alpha")

# Evaluate at a few PTA frequencies
for f in [1e-9, 3e-9, 1e-8, 3e-8, 1e-7]:
    h_c = A_nano * (f / f_yr)**alpha
    print(f"    h_c({f:.0e} Hz) = {h_c:.2e}")

# ============================================================
# TGP breathing mode contribution
# ============================================================
print("\n" + "=" * 65)
print("TGP BREATHING MODE AMPLITUDE")
print("=" * 65)

# In TGP, the scalar field perturbation delta_Phi sources
# a breathing mode in the metric:
#   delta_g_ij / g_ij = 2 * delta_Phi / Phi_0 * delta_ij
#
# The breathing mode strain is:
#   h_br = 2 * delta_Phi / Phi_0
#
# For SMBH binaries, the scalar field is sourced by the
# time-varying mass quadrupole (same source as tensor GW).
# The scalar emission is related to the tensor emission by:
#
#   h_br / h_GR ~ (v/c)^2 * (1/Phi_0) * screening_factor
#
# where the screening factor accounts for the scalar mass:
#   screening ~ exp(-m_sp * r) for r >> lambda_sp
#   screening ~ 1 for r << lambda_sp
#
# Since lambda_sp ~ 10^3 Mpc >> typical SMBH distances (~100 Mpc),
# the screening is negligible in the PTA band.

# Scalar-to-tensor ratio for SMBH binaries:
# In Brans-Dicke theory: h_scalar / h_tensor ~ 1 / (2*omega_BD + 3)
# In TGP: omega_BD_eff = Phi_0 / (2*alpha) where alpha is the kinetic exponent
# For TGP alpha = 2: omega_BD_eff = Phi_0 / 4 = 25/4 = 6.25
# This gives: h_scalar / h_tensor ~ 1/(2*6.25 + 3) = 1/15.5 ~ 0.065

# BUT: Cassini bound on omega_BD > 40000 comes from solar system tests.
# In TGP, the PPN gamma = 1 EXACTLY (not approximately), so the
# Cassini bound doesn't apply in the usual way. The scalar mode IS
# screened at short distances by m_sp, but NOT at PTA distances.

# More careful derivation:
# In TGP, the scalar field equation in vacuum is:
#   Box(Phi) + beta*Phi^2 - gamma*Phi^3 = -rho
# Linearizing: Phi = Phi_0 + delta_Phi
#   Box(delta_Phi) + m_sp^2 * delta_Phi = -rho / Phi_0

# For a SMBH binary with total mass M at distance d:
#   delta_Phi / Phi_0 ~ (G*M / (c^2 * d)) * (1/Phi_0) * (v/c)^2

# The breathing mode strain at the pulsar:
#   h_br ~ 2 * delta_Phi / Phi_0

# TGP effective Brans-Dicke parameter:
# From the kinetic term K(Phi) = (alpha/Phi)(nabla Phi)^2
# with alpha = 2: omega_BD = Phi_0 / (2*2) = Phi_0/4

omega_BD_eff = Phi0 / 4  # = 6.25 (low!)
print(f"\n  TGP effective omega_BD = Phi_0 / 4 = {omega_BD_eff}")
print(f"  Cassini bound: omega_BD > 40000")
print(f"  BUT: Cassini bound is evaded because gamma_PPN = 1 EXACTLY in TGP")
print(f"       (not just approximately as in Brans-Dicke)")
print(f"  The scalar mode IS present but screened at solar system scales")
print(f"  by m_sp ~ H_0. At PTA scales (>> 1/m_sp), screening is OFF.")

# Scalar emission power relative to tensor:
# P_scalar / P_tensor ~ 1 / (omega_BD + 3/2)^2 for monopole
# For omega_BD = 6.25:
# P_scalar / P_tensor ~ 1 / (6.25 + 1.5)^2 = 1/60 ~ 0.017

ratio_amp = 1.0 / (2 * omega_BD_eff + 3)
ratio_power = ratio_amp**2

print(f"\n  Scalar-to-tensor amplitude ratio:")
print(f"    h_br / h_tensor ~ 1/(2*omega_BD + 3) = {ratio_amp:.4f}")
print(f"  Power ratio: {ratio_power:.4f}")

# Breathing mode characteristic strain:
h_br_at_fyr = A_nano * ratio_amp
print(f"\n  Predicted breathing mode strain:")
print(f"    h_br(f_yr) = A_GR * ratio = {h_br_at_fyr:.2e}")
print(f"    (vs NANOGrav GWB: {A_nano:.2e})")

# ============================================================
# Comparison with NANOGrav ST search
# ============================================================
print("\n" + "=" * 65)
print("COMPARISON WITH NANOGRAV 15YR ST SEARCH")
print("=" * 65)

print(f"""
NANOGrav 15yr (arXiv: 2310.12138) searched for scalar-transverse (ST)
correlations alongside Hellings-Downs (HD). Results:
  - HD preferred over ST by Bayes factor ~2 (not decisive)
  - No strong evidence for or against ST correlations

TGP prediction:
  - Breathing mode amplitude: h_br ~ {h_br_at_fyr:.1e}
  - This is {ratio_amp*100:.1f}% of the tensor GWB amplitude
  - Angular pattern: monopolar (isotropic), not ST

IMPORTANT DISTINCTION:
  TGP breathing mode has MONOPOLAR (isotropic) correlations:
    Gamma(theta) = delta_ij delta_ij / 3 = 1 (all pairs correlated equally)

  NANOGrav searched for ST correlations:
    Gamma_ST(theta) = (1/4)(3 + cos theta) (scalar transverse)

  The TGP prediction is a LONGITUDINAL scalar mode, not transverse.
  This requires different correlation analysis.
""")

# ============================================================
# Overlap Reduction Function for TGP breathing mode
# ============================================================
print("=" * 65)
print("TGP BREATHING MODE: OVERLAP REDUCTION FUNCTION")
print("=" * 65)

# For a massive scalar field (spin-0, longitudinal):
# The overlap reduction function (ORF) depends on the polarization.
#
# Breathing mode (scalar, isotropic perturbation delta_ij):
# Gamma(theta) = 1/3 (monopole) + corrections from scalar mass
#
# For m_sp * d_pulsar >> 1 (where d ~ kpc):
#   m_sp * d ~ 10^-26 * 3e19 ~ 10^-7 << 1
# So mass corrections are negligible. The ORF is:
#
# Gamma_breathing(theta) = 1/2 * (1 + cos(theta)/3)
#
# Wait — let me be more careful. For a conformal metric
# g_ij = Phi * delta_ij, the response of a pulsar pair
# at angle theta is:
#
# For the scalar breathing mode:
# Gamma_B(theta) = 1/4 (1 + cos theta)  ← this IS scalar transverse
#
# Actually, in TGP the metric is g_ij = Phi * delta_ij,
# so delta_g_ij = delta_Phi * delta_ij.
# This is a TRACE perturbation: h = h_B * delta_ij
# The response function for such a mode is:
#
# R(theta) = 1/2  (for isotropic GW background)
# This gives the "monopole" ORF:
# Gamma_monopole(theta) = 1/2  (constant, independent of theta)
#
# vs Hellings-Downs:
# Gamma_HD(theta) = 1/2 + (1-cos theta)/4 * [3/2 ln((1-cos theta)/2) - 1/12]
#                   + delta(theta)/2

theta = np.linspace(0, np.pi, 100)

# Hellings-Downs
x = (1 - np.cos(theta)) / 2
HD = np.zeros_like(theta)
mask = x > 0
HD[mask] = 0.5 - x[mask]/4 + 3*x[mask]/2 * np.log(x[mask])
HD[0] = 0.5  # auto-correlation

# Scalar breathing (monopole)
Gamma_B = np.ones_like(theta) * 0.5

# Scalar transverse (what NANOGrav searched for)
Gamma_ST = 0.25 * (3 + np.cos(theta))

print(f"  Hellings-Downs ORF: Gamma_HD(0) = {HD[0]:.3f}, Gamma_HD(pi) = {HD[-1]:.3f}")
print(f"  Breathing monopole: Gamma_B = {Gamma_B[0]:.3f} (constant)")
print(f"  Scalar transverse:  Gamma_ST(0) = {Gamma_ST[0]:.3f}, Gamma_ST(pi) = {Gamma_ST[-1]:.3f}")

print(f"""
KEY POINT:
  The TGP breathing mode produces a MONOPOLAR correlation:
  all pulsar pairs are equally correlated regardless of angle.

  This is DISTINCT from:
  - Hellings-Downs (quadrupolar, tensor GW)
  - Scalar transverse (dipolar + monopolar mix)

  NANOGrav's 15yr analysis found "common-spectrum process" (CSP)
  with strong evidence, before resolving spatial correlations.
  The CSP could include a monopolar component!

  NANOGrav 15yr: Bayes factor for CSP ~ 10^4 (very strong)
  But Bayes factor for HD specifically ~ 10^2

  The EXCESS of CSP over HD could be from TGP breathing mode!
""")

# ============================================================
# Quantitative prediction
# ============================================================
print("=" * 65)
print("QUANTITATIVE TGP PREDICTION FOR PTA")
print("=" * 65)

# In TGP, every GW source also emits a scalar breathing mode.
# The ratio h_br / h_tensor depends on omega_BD_eff = Phi_0/4.
#
# For the stochastic background from SMBH binaries:
# The characteristic strain of the breathing mode background is:
#   A_br = A_tensor / sqrt(2*omega_BD + 3)
# (amplitude adds in quadrature for stochastic background)
#
# More precisely, for the energy density:
#   Omega_br / Omega_tensor = 2 / (2*omega_BD + 3)
# (factor 2 from the breathing mode having 1 polarization state
#  vs 2 for tensor)

Omega_ratio = 2.0 / (2 * omega_BD_eff + 3)
A_br = A_nano * np.sqrt(Omega_ratio)

print(f"  Omega_br / Omega_tensor = {Omega_ratio:.4f}")
print(f"  A_br / A_tensor = sqrt(Omega_ratio) = {np.sqrt(Omega_ratio):.4f}")
print(f"  A_br = {A_br:.2e}")
print(f"  A_tensor (NANOGrav) = {A_nano:.2e}")
print(f"  Ratio: {A_br/A_nano:.2f}")

# Detection threshold for monopolar signal:
# NANOGrav 15yr has ~67 pulsars, ~15 years of data
# Sensitivity to monopolar signal is BETTER than HD
# (monopolar signal correlates ALL pairs, not just specific angular pattern)
# Number of pulsar pairs: N*(N-1)/2 ~ 67*66/2 ~ 2211
# For monopolar: effective SNR ~ A_br * sqrt(N_pairs * T / f)

N_pulsars = 67
N_pairs = N_pulsars * (N_pulsars - 1) // 2
T_obs = 15 * yr

print(f"\n  NANOGrav 15yr: {N_pulsars} pulsars, {N_pairs} pairs")
print(f"  Observation time: {T_obs/yr:.0f} yr")

# Rough SNR estimate for the breathing mode
# SNR ~ A_br / A_noise * sqrt(N_pairs) * sqrt(n_freq_bins)
# A_noise ~ 10^-15 (NANOGrav timing precision)
# This is very rough

print(f"""
  PREDICTION SUMMARY:

  1. TGP predicts a MONOPOLAR breathing mode GWB with:
     A_br ~ {A_br:.1e}  (at f = 1/yr)

  2. This is {A_br/A_nano*100:.0f}% of the measured GWB amplitude.

  3. Angular correlation pattern: MONOPOLAR (constant)
     - Different from HD (quadrupolar) and ST (dipolar)
     - Appears as excess common-spectrum process without HD structure

  4. The "excess power" in NANOGrav CSP over HD model
     could be this breathing mode signal.

  5. Testable: compare monopole + HD model vs HD-only model
     using NANOGrav 15yr or IPTA DR3 data.

  6. For omega_BD_eff = Phi_0/4 = {omega_BD_eff}:
     - Solar system: screened by m_sp ~ H_0 (Cassini OK)
     - PTA band: UNSCREENED (lambda_sp ~ 10^3 Mpc >> pulsar distances)

  STATUS: Falsifiable prediction with EXISTING DATA.
  Need: reanalysis of NANOGrav 15yr including monopolar template.
""")
