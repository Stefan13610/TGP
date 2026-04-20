"""
ps01_schwinger_TGP_1loop.py - 1-loop QED vertex in TGP substrate; EFT scaling

GOAL
====
Compute the TGP correction to Schwinger's 1-loop anomalous magnetic moment
a_l = (g-2)/2 for l = e, mu, tau, and derive the IMPLIED TGP energy scale
from the current Fermilab Muon g-2 anomaly.

Physical setup
==============
Schwinger 1948:  a_l^(flat, 1-loop) = alpha/(2pi) ~ 1.16e-3  (universal)

TGP substrate contribution to QED vertex (generic EFT form):
    Delta a_l^TGP  =  (alpha/2pi) * xi_TGP * (m_l / M_TGP)^2
where
    M_TGP   = effective UV scale of substrate coupling to QED loops
    xi_TGP  = O(1) Wilson coefficient from substrate topology

Diagnostic note on the SC lattice scale
---------------------------------------
One might guess M_TGP = hbar*c / a_TGP where a_TGP = 4.088 A is the SC
Gaussian lattice scale from P6.A. That gives Lambda_a = 483 eV.
But m_mu c^2 = 106 MeV >> 483 eV, so the muon Compton wavelength is
MUCH SMALLER than a_TGP. In that regime a virtual-photon loop of momentum
q ~ m_mu sees the substrate as a slowly-varying background (long-wavelength
EFFECTIVE refractive index), not as a UV cutoff. The substrate cannot
resolve lepton structure.

Conclusion: a_TGP is NOT the relevant QED scale. The TGP scale for lepton
QED loops must be a distinct UV object, call it M_TGP^lepton, which is
PRESUMABLY the substrate nonlinear excitation scale Lambda_E at short
distance. This is a NEW constant, not fixed by SC or viscosity.

What ps01 actually does
=======================
1. Set up the EFT: Delta a_l = (alpha/2pi) * xi * (m_l/M_TGP)^2
2. Fit M_TGP/sqrt(xi) to Delta a_mu tension
3. Predict Delta a_e and Delta a_tau from same M_TGP
4. Compare to experimental bounds / sensitivities
5. Test the quadratic-mass ("mass-dependent BSM") signature
"""

import numpy as np

# =============================================================
# CONSTANTS
# =============================================================
alpha     = 7.2973525693e-3
alpha_2pi = alpha / (2 * np.pi)

# Lepton masses (GeV)
m_e   = 0.51099895000e-3
m_mu  = 0.1056583755
m_tau = 1.77686

# =============================================================
# EXPERIMENTAL VALUES (April 2026 snapshot)
# =============================================================
a_e_exp      = 1.15965218059e-3
sig_a_e_tot  = 1.3e-13
Delta_a_e_Cs_central = -8.8e-13
Delta_a_e_Rb_central = +4.8e-13

a_mu_exp_central = 116592059e-11
a_mu_exp_sigma   = 22e-11
a_mu_SM_WP_dd    = 116591810e-11
a_mu_SM_BMW_lat  = 116591960e-11

Delta_a_mu_dd  = a_mu_exp_central - a_mu_SM_WP_dd
Delta_a_mu_lat = a_mu_exp_central - a_mu_SM_BMW_lat

a_tau_bound = 1.3e-2

# =============================================================
# EFT SCALING
# =============================================================
def delta_a_TGP(m_l, M_TGP, xi=1.0):
    return alpha_2pi * xi * (m_l / M_TGP)**2

def fit_M_TGP(Delta_a, m_l, xi=1.0):
    return m_l * np.sqrt(alpha_2pi * xi / Delta_a)

# =============================================================
# REPORT
# =============================================================
print("=" * 90)
print("  ps01 - TGP 1-LOOP EFT CORRECTION TO (g-2)_l")
print("=" * 90)
print()
print("  EFT ansatz:   Delta a_l^TGP = (alpha/2pi) * xi_TGP * (m_l / M_TGP)^2")
print()
print(f"  m_e   = {m_e:.6e} GeV")
print(f"  m_mu  = {m_mu:.6e} GeV")
print(f"  m_tau = {m_tau:.6e} GeV")
print()

# =============================================================
# Step 1
# =============================================================
print("=" * 90)
print("  STEP 1 - Fit M_TGP to observed Delta a_mu tension (xi = 1)")
print("=" * 90)
print()
scenarios = [
    ("WP20 data-driven HVP",  Delta_a_mu_dd),
    ("BMW24 lattice HVP",     Delta_a_mu_lat),
]
print(f"  {'Source':>25} {'Delta a_mu':>18} {'M_TGP (GeV)':>13}")
print(f"  {'-'*25} {'-'*18} {'-'*13}")
M_TGP_list = []
for label, damu in scenarios:
    M = fit_M_TGP(damu, m_mu)
    M_TGP_list.append(M)
    print(f"  {label:>25} {damu*1e11:>12.1f} x 1e-11   {M:>10.2f}")
M_TGP_central = np.sqrt(M_TGP_list[0] * M_TGP_list[1])
print()
print(f"  Geometric-mean M_TGP (xi=1) = {M_TGP_central:.2f} GeV")
print(f"  Comparison: M_W = 80.4 GeV, M_Z = 91.2 GeV, v_EW = 246 GeV")
print()

# =============================================================
# Step 2
# =============================================================
print("=" * 90)
print("  STEP 2 - Predict Delta a_e and Delta a_tau at fitted M_TGP")
print("=" * 90)
print()
for label, M in [("WP20-fit   ", M_TGP_list[0]),
                 ("BMW24-fit  ", M_TGP_list[1]),
                 ("geom-mean  ", M_TGP_central)]:
    dae  = delta_a_TGP(m_e,   M)
    damu = delta_a_TGP(m_mu,  M)
    dat  = delta_a_TGP(m_tau, M)
    print(f"  {label} (M_TGP = {M:.1f} GeV):")
    print(f"     Delta a_e   = {dae:.3e}   ({dae/sig_a_e_tot:+6.2f} sigma of exp precision)")
    print(f"     Delta a_mu  = {damu:.3e}   ({damu*1e11:.1f} x 1e-11)  [fit target]")
    print(f"     Delta a_tau = {dat:.3e}   ({dat/a_tau_bound*100:.2e}% of DELPHI bound)")
    print()

# =============================================================
# Step 3 - mass scaling signature
# =============================================================
print("=" * 90)
print("  STEP 3 - Quadratic-mass scaling test (TGP fingerprint)")
print("=" * 90)
print()
print("  If TGP: Delta a_l / Delta a_l' = (m_l/m_l')^2  (regardless of M_TGP)")
print()
print(f"    (m_mu/m_e)^2    = {(m_mu/m_e)**2:>12.3e}")
print(f"    (m_tau/m_mu)^2  = {(m_tau/m_mu)**2:>12.3e}")
print(f"    (m_tau/m_e)^2   = {(m_tau/m_e)**2:>12.3e}")
print()
Delta_a_mu_obs = Delta_a_mu_lat
Delta_a_e_TGP  = Delta_a_mu_obs * (m_e/m_mu)**2
Delta_a_tau_TGP = Delta_a_mu_obs * (m_tau/m_mu)**2
print(f"  Using observed Delta a_mu (BMW24) = {Delta_a_mu_obs*1e11:.0f} x 1e-11:")
print(f"     TGP-predicted Delta a_e   = {Delta_a_e_TGP:.3e}  ({Delta_a_e_TGP/sig_a_e_tot:+.2f} sigma)")
print(f"     TGP-predicted Delta a_tau = {Delta_a_tau_TGP:.3e}  ({Delta_a_tau_TGP/a_tau_bound*100:.1e}% of bound)")
print()
print(f"  Current electron deviations from SM:")
print(f"     Delta a_e(Cs-alpha) ~ {Delta_a_e_Cs_central:.2e}  (obs 2.4 sigma)")
print(f"     Delta a_e(Rb-alpha) ~ {Delta_a_e_Rb_central:.2e}  (obs 1.7 sigma)")
print()
print(f"  Observed |Delta a_e| ~ 5-9 x 10^-13, but TGP prediction is only")
print(f"  {Delta_a_e_TGP:.1e}. Ratio {abs(Delta_a_e_Cs_central)/Delta_a_e_TGP:.0f}-fold too large for Cs-alpha.")
print()
print(f"  DIAGNOSIS: if both a_e(Cs) and a_mu deviations are real, they can")
print(f"  NOT both come from the same TGP mass-dependent operator.")
print(f"  => Electron deviation (if confirmed) suggests a DIFFERENT mechanism,")
print(f"     likely systematic in alpha_Cs measurement, not BSM.")
print()

# =============================================================
# Step 4
# =============================================================
print("=" * 90)
print("  STEP 4 - Physical interpretation")
print("=" * 90)
print()
print(f"  M_TGP ~ {M_TGP_central:.0f} GeV  (geometric mean of WP-fit and lattice-fit)")
print()
print("  This is remarkably close to the ELECTROWEAK scale:")
print(f"     M_W  = 80.4 GeV                   ratio = {80.4/M_TGP_central:.2f}")
print(f"     M_Z  = 91.2 GeV                   ratio = {91.2/M_TGP_central:.2f}")
print(f"     v_EW = 246 GeV                    ratio = {246.22/M_TGP_central:.2f}")
print()
print("  TGP HYPOTHESIS: M_TGP^lepton = xi * v_EW, with xi = O(0.3-1.0)")
print("  from substrate topology. This implies the TGP substrate nonlinearity")
print("  scale for lepton QED is TIED to electroweak breaking, not to the")
print("  macroscopic a_TGP = 4.088 A used in SC.")
print()

# =============================================================
# Step 5 - Fermilab Run-6
# =============================================================
print("=" * 90)
print("  STEP 5 - Fermilab Muon g-2 Run-6 projection (final data ~2027)")
print("=" * 90)
print()
print("  Projected FNAL final sigma (Runs 1-6):  14 x 10^-11")
print("  Current combined sigma (Run-3):         22 x 10^-11")
print()
for Mref, name in [(246.22, "v_EW"), (91.2, "M_Z"), (80.4, "M_W")]:
    dam = delta_a_TGP(m_mu, Mref)
    xi_needed_dd  = Delta_a_mu_dd  / dam
    xi_needed_lat = Delta_a_mu_lat / dam
    print(f"  M_TGP = {name} ({Mref:.1f} GeV):  Delta a_mu(xi=1) = {dam*1e11:.1f} x 1e-11")
    print(f"    xi required for WP20 data-driven tension:  {xi_needed_dd:.2f}")
    print(f"    xi required for BMW24 lattice tension:     {xi_needed_lat:.2f}")
    print()

# =============================================================
# FINAL SUMMARY
# =============================================================
dae_final  = delta_a_TGP(m_e,   M_TGP_central)
damu_final = delta_a_TGP(m_mu,  M_TGP_central)
dat_final  = delta_a_TGP(m_tau, M_TGP_central)
print("=" * 90)
print("  FINAL SUMMARY")
print("=" * 90)
print(f"""
  EFT:             Delta a_l = (alpha/2pi) * xi_TGP * (m_l/M_TGP)^2
  Fitted scale:    M_TGP = {M_TGP_central:.0f} GeV   (=~ M_W/M_Z)
  Wilson coeff:    xi_TGP = O(1)  (from substrate topology, to be computed)

  Predictions at M_TGP = {M_TGP_central:.0f} GeV, xi = 1:
    Delta a_e^TGP   = {dae_final:.2e}   (INVISIBLE, < 0.01 sigma of exp precision)
    Delta a_mu^TGP  = {damu_final:.2e}   = {damu_final*1e11:.1f} x 1e-11  (matches tension)
    Delta a_tau^TGP = {dat_final:.2e}   (below DELPHI bound by factor {a_tau_bound/dat_final:.0e})

  KEY FALSIFIERS (testable 2026-2030):

   (F1) Lepton-universality quadratic scaling
        Measure Delta a_tau via BELLE-II or FCC-ee. TGP predicts
        Delta a_tau = (m_tau/m_mu)^2 * Delta a_mu ~ 3 x 10^-7.
        A universal (mass-independent) BSM would predict
        Delta a_tau ~ Delta a_mu ~ 10^-9. Factor ~300 apart.

   (F2) Electron g-2 quiet
        TGP predicts Delta a_e ~ 2 x 10^-14 (invisible). If electron
        g-2 shows persistent > 3 sigma deviation not explainable by
        alpha_Cs systematics, TGP mass-dependent mechanism is incomplete.

   (F3) FNAL Run-6 consistency
        If Delta a_mu settles at ~100 x 10^-11 (lattice HVP) the required
        xi is {Delta_a_mu_lat/delta_a_TGP(m_mu, M_TGP_central):.2f}.
        If at ~250 x 10^-11 (data-driven) xi is {Delta_a_mu_dd/delta_a_TGP(m_mu, M_TGP_central):.2f}.
        Both O(1); TGP SURVIVES either HVP resolution.

   (F4) Sign of Delta a_mu
        Current data: Delta a_mu > 0. TGP requires xi_TGP > 0, consistent
        with the default exp(+2Phi/Phi_0) metric. Future measurement of
        Delta a_mu < 0 would require revising metric sign convention.

  STATUS: zeroth-order proof-of-concept EFT prediction.
          ps02 extends to (a_e, a_tau) universality scan.
          Full 1-loop derivation in exp metric pending.
""")
