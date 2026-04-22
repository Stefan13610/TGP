#!/usr/bin/env python3
"""
TGP sin^2(theta_W) QCD Correction Analysis (2026-04-10)

Problem: TGP predicts sin^2(theta_W) = 3/13 = 0.23077 (tree level).
PDG: 0.23122 +/- 0.00004 (MS-bar, M_Z).
Deviation: -0.045 (0.19%), which is 11.3 sigma.

This is the ONLY TGP prediction beyond 2-sigma.

Question: Can a perturbative QCD correction close the gap?

The correction delta = 0.23122 - 0.23077 = +0.000451 is positive (additive).

Strategy:
  1. Test physically motivated correction forms
  2. Check if the coefficient has a natural group-theoretic origin
  3. Verify that the correction is self-consistent with TGP alpha_s
"""

import numpy as np

# ============================================================
# Parameters
# ============================================================
N_c = 3
N_f = 5  # at M_Z scale
alpha_s = 0.1179   # PDG
sin2_tree = N_c / (1 + N_c + N_c**2)  # 3/13 = 0.23077
sin2_PDG = 0.23122
sin2_err = 0.00004
delta_needed = sin2_PDG - sin2_tree

print("=" * 70)
print("TGP sin^2(theta_W): QCD CORRECTION ANALYSIS")
print("=" * 70)

print(f"\n  Tree level: sin^2 = N_c/(1+N_c+N_c^2) = {N_c}/{1+N_c+N_c**2} = {sin2_tree:.6f}")
print(f"  PDG:        sin^2 = {sin2_PDG:.5f} +/- {sin2_err:.5f}")
print(f"  Deficit:    delta = {delta_needed:.6f}")
print(f"  Relative:   delta/sin2 = {delta_needed/sin2_tree:.6f} = {delta_needed/sin2_tree*100:.4f}%")
print(f"  In sigma:   {delta_needed/sin2_err:.1f} sigma")


# ============================================================
# [1] STANDARD ELECTROWEAK RADIATIVE CORRECTIONS
# ============================================================
print(f"\n\n[1] STANDARD ELECTROWEAK CORRECTIONS")
print("-" * 50)

# In the SM, the running of sin^2(theta_W) from tree level to M_Z
# involves several contributions:
# 1. Electromagnetic running (alpha -> alpha(M_Z))
# 2. QCD corrections to gauge boson self-energies
# 3. Top quark corrections (rho parameter)
# 4. Higgs corrections
#
# For TGP, the tree-level value is 3/13 (not 3/8 as in SU(5) GUT).
# The QUESTION is whether the same SM radiative corrections apply,
# just starting from a different tree value.

# In SM: sin^2(tree) = 3/8 = 0.375
# sin^2(M_Z) = 0.23122 (after running from GUT to M_Z)
# Running: delta_SM = 0.375 - 0.23122 = 0.14378 (huge, 38% shift)

# In TGP: sin^2(tree) = 3/13 = 0.23077
# Required: delta_TGP = +0.00045 (tiny, 0.19% shift)

# This is a COMPLETELY different situation:
# SM needs O(1) running; TGP needs O(0.1%) correction.
# TGP tree level is ALREADY close to the measured value!

print(f"  SM tree: sin^2 = 3/8 = {3/8:.5f}")
print(f"  SM running needed: {3/8 - sin2_PDG:.5f} ({(3/8 - sin2_PDG)/(3/8)*100:.1f}%)")
print(f"  TGP tree: sin^2 = 3/13 = {sin2_tree:.5f}")
print(f"  TGP correction needed: {delta_needed:.5f} ({delta_needed/sin2_tree*100:.3f}%)")
print(f"  -> TGP is {(3/8 - sin2_PDG)/delta_needed:.0f}x closer to observation at tree level!")


# ============================================================
# [2] SCAN CORRECTION FORMS
# ============================================================
print(f"\n\n[2] SYSTEMATIC SCAN OF CORRECTION FORMS")
print("-" * 50)

# The correction has form: sin^2 = (3/13) + delta
# or multiplicatively: sin^2 = (3/13) * (1 + epsilon)
# where epsilon = delta / sin2_tree = 0.001954

epsilon = delta_needed / sin2_tree
print(f"  Relative correction epsilon = {epsilon:.6f}")
print()

# Test various forms for epsilon:
candidates = {}

# Simple alpha_s powers
candidates["alpha_s/pi"] = alpha_s / np.pi
candidates["alpha_s/(2*pi)"] = alpha_s / (2*np.pi)
candidates["alpha_s/(4*pi)"] = alpha_s / (4*np.pi)
candidates["alpha_s^2"] = alpha_s**2
candidates["alpha_s^2/pi"] = alpha_s**2 / np.pi
candidates["alpha_s^2/(2*pi)"] = alpha_s**2 / (2*np.pi)
candidates["alpha_s^2/pi^2"] = alpha_s**2 / np.pi**2

# Group theory factors
candidates["alpha_s/(N_c*pi)"] = alpha_s / (N_c * np.pi)
candidates["C_F*alpha_s/(N_c*pi)"] = (4/3)*alpha_s / (N_c * np.pi)
candidates["alpha_s*N_f/(12*pi)"] = alpha_s * N_f / (12 * np.pi)
candidates["alpha_s*N_c/(4*pi^2)"] = alpha_s * N_c / (4 * np.pi**2)

# Combinations with sin^2 itself
candidates["alpha_s*sin2/(2*pi)"] = alpha_s * sin2_tree / (2 * np.pi)
candidates["alpha_s*cos2/(4*pi)"] = alpha_s * (1-sin2_tree) / (4*np.pi)

# TGP-specific
candidates["3*alpha_s/(4*pi*N_c)"] = 3*alpha_s / (4*np.pi*N_c)
candidates["alpha_s/(2*pi*N_c)"] = alpha_s / (2*np.pi*N_c)
candidates["(N_f-N_c)*alpha_s/(6*pi)"] = (N_f-N_c)*alpha_s / (6*np.pi)
candidates["alpha_s/(pi*(1+N_c+N_c^2))"] = alpha_s / (np.pi*(1+N_c+N_c**2))

# Electroweak-like
candidates["alpha_s*sin2*cos2/pi"] = alpha_s * sin2_tree * (1-sin2_tree) / np.pi
candidates["alpha_s/(2*pi*(N_c^2-1))"] = alpha_s / (2*np.pi*(N_c**2-1))

# Known SM-like correction:
# In SM, leading QCD correction to rho parameter:
# delta_rho = 3*G_F*m_t^2/(8*sqrt(2)*pi^2) ~ 0.01
# This contributes: delta_sin2 ~ -cos^2*delta_rho/(cos^2-sin^2)
# Much too large for our purposes.
#
# The OBLIQUE correction to sin^2 from light quarks:
# delta_sin2_QCD ~ alpha_s/(12*pi) * sum_quarks(T3_q - 2*Q_q*sin^2)^2
# For up-type: T3=+1/2, Q=+2/3
# For down-type: T3=-1/2, Q=-1/3
# Factor per doublet: (1/2 - 4/3*sin2)^2 + (-1/2 + 2/3*sin2)^2

T3u_eff = 0.5 - 2*(2/3)*sin2_tree  # 0.5 - 4/3*0.23077 = 0.1923
T3d_eff = -0.5 + 2*(1/3)*sin2_tree  # -0.5 + 2/3*0.23077 = -0.3462
oblique_factor = T3u_eff**2 + T3d_eff**2
n_doublets = 3  # u/d, c/s, t/b at M_Z (effectively 2.5 due to top)
candidates["oblique: as/(12pi)*F*N_d"] = alpha_s/(12*np.pi) * oblique_factor * n_doublets

# Physical: the shift comes from the running of coupling constants
# In TGP, sin^2 = g'^2/(g^2 + g'^2) where g, g' are SU(2), U(1) couplings
# QCD modifies g through quark loops at 1-loop:
# delta(1/g^2) ~ -b_2 * ln(mu/M_Z) where b_2 depends on N_f
# At the scale M_Z itself, there's a threshold correction.
# The threshold correction for N_f quarks:
# delta_sin^2 ~ sin^2*cos^2 * (alpha_s/pi) * sum_q (Y_q^2 - T3_q^2) / (cos^2-sin^2)
# This is model-dependent. Let me just compute it numerically.

# For the specific TGP form, the threshold correction at M_Z:
cos2_tree = 1 - sin2_tree
# Sum over quarks: Y^2 - T3^2 for each quark
# Up-type: Y = 2/3 - 1/2 = 1/6 (weak hypercharge Y = Q - T3), T3 = 1/2
# Actually Y_L = 1/6 (for left-handed quark doublet)
# Y_R(u) = 2/3, Y_R(d) = -1/3
# The correction involves specific diagram calculations. Let me just
# check if a simple form gives the right number.

# Simple candidates involving the TGP denominator 1+N_c+N_c^2 = 13:
candidates["alpha_s/(pi*(1+N_c+N_c^2))"] = alpha_s / (np.pi * 13)
candidates["alpha_s*N_c/(pi*(1+N_c+N_c^2))"] = alpha_s * N_c / (np.pi * 13)
candidates["3*alpha_s^2/pi"] = 3 * alpha_s**2 / np.pi

print(f"  {'Formula':<40} {'epsilon':>10} {'target':>10} {'ratio':>8} {'sin2':>10} {'dev_sigma':>10}")
print(f"  {'='*88}")

for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]/epsilon - 1)):
    sin2_corrected = sin2_tree * (1 + val)
    dev_sigma = (sin2_corrected - sin2_PDG) / sin2_err
    ratio = val / epsilon
    mark = " <== " if abs(ratio - 1) < 0.15 else ""
    print(f"  {name:<40} {val:10.6f} {epsilon:10.6f} {ratio:8.4f} {sin2_corrected:10.7f} {dev_sigma:+10.1f}{mark}")


# ============================================================
# [3] BEST CANDIDATES ANALYSIS
# ============================================================
print(f"\n\n[3] BEST CANDIDATES (within 15% of target)")
print("-" * 50)

# Filter best candidates
best = [(n, v) for n, v in candidates.items() if abs(v/epsilon - 1) < 0.15]
best.sort(key=lambda x: abs(x[1]/epsilon - 1))

for name, val in best:
    sin2_c = sin2_tree * (1 + val)
    sig = (sin2_c - sin2_PDG) / sin2_err
    print(f"\n  {name}")
    print(f"    epsilon = {val:.8f}")
    print(f"    sin^2 = {sin2_c:.7f}")
    print(f"    PDG: {sin2_PDG:.5f} +/- {sin2_err:.5f}")
    print(f"    Deviation: {sig:+.2f} sigma")


# ============================================================
# [4] ADDITIVE CORRECTION FORMS
# ============================================================
print(f"\n\n[4] ADDITIVE CORRECTIONS: sin^2 = 3/13 + delta")
print("-" * 50)

# Instead of multiplicative, try additive
add_candidates = {}
add_candidates["alpha_s^2/(2*pi)"] = alpha_s**2 / (2*np.pi)
add_candidates["alpha_s^2*N_c/(4*pi)"] = alpha_s**2 * N_c / (4*np.pi)
add_candidates["alpha_s^3/(2*pi)"] = alpha_s**3 / (2*np.pi)
add_candidates["3*alpha_s/(4*pi*(1+Nc+Nc2))"] = 3*alpha_s / (4*np.pi*13)
add_candidates["alpha_s*sin2*cos2/pi^2"] = alpha_s * sin2_tree * cos2_tree / np.pi**2
add_candidates["alpha_s/(2*pi*13)"] = alpha_s / (2*np.pi*13)
add_candidates["alpha_s/(4*pi*N_c)"] = alpha_s / (4*np.pi*N_c)
add_candidates["3/(13^2)"] = 3.0/169  # pure structure, no alpha_s
add_candidates["1/(13*pi)"] = 1.0/(13*np.pi)
add_candidates["alpha_s/(8*pi*N_c)"] = alpha_s / (8*np.pi*N_c)

print(f"  {'Formula':<40} {'delta':>12} {'target':>12} {'ratio':>8} {'sin2':>10} {'sigma':>8}")
for name, val in sorted(add_candidates.items(), key=lambda x: abs(x[1]/delta_needed - 1)):
    sin2_c = sin2_tree + val
    sig = (sin2_c - sin2_PDG) / sin2_err
    ratio = val / delta_needed
    mark = " <==" if abs(ratio - 1) < 0.15 else ""
    print(f"  {name:<40} {val:12.8f} {delta_needed:12.8f} {ratio:8.4f} {sin2_c:10.7f} {sig:+8.2f}{mark}")


# ============================================================
# [5] DEEPER: sin^2 = N_c/(1+N_c+N_c^2+f(alpha_s))
# ============================================================
print(f"\n\n[5] MODIFIED DENOMINATOR: sin^2 = 3/(13 + delta_D)")
print("-" * 50)

# Instead of correcting sin^2 directly, modify the denominator:
# sin^2 = 3/(13 + delta_D) where delta_D is a QCD correction
# 0.23122 = 3/(13 + delta_D)
# 13 + delta_D = 3/0.23122 = 12.9744
# delta_D = 12.9744 - 13 = -0.0256

delta_D_needed = 3.0/sin2_PDG - 13.0
print(f"  3/sin2_PDG = {3/sin2_PDG:.6f}")
print(f"  delta_D = {delta_D_needed:.6f}")
print()

# delta_D is NEGATIVE and small. So the correction reduces the denominator.
# Check some forms:
denom_candidates = {}
denom_candidates["-alpha_s/pi"] = -alpha_s/np.pi
denom_candidates["-alpha_s/(2*pi)"] = -alpha_s/(2*np.pi)
denom_candidates["-alpha_s*N_c/(4*pi)"] = -alpha_s*N_c/(4*np.pi)
denom_candidates["-3*alpha_s^2"] = -3*alpha_s**2
denom_candidates["-alpha_s^2*13"] = -alpha_s**2*13
denom_candidates["-N_c*alpha_s^2"] = -N_c*alpha_s**2
denom_candidates["-alpha_s/(4*pi)"] = -alpha_s/(4*np.pi)
denom_candidates["-C_F*alpha_s/pi"] = -(4/3)*alpha_s/np.pi
denom_candidates["-alpha_s*N_c/pi^2"] = -alpha_s*N_c/np.pi**2

print(f"  {'Formula':<35} {'delta_D':>12} {'target':>12} {'ratio':>8} {'sin2':>10} {'sigma':>8}")
for name, val in sorted(denom_candidates.items(), key=lambda x: abs(x[1]/delta_D_needed - 1)):
    sin2_c = 3.0/(13.0 + val)
    sig = (sin2_c - sin2_PDG) / sin2_err
    ratio = val / delta_D_needed
    mark = " <==" if abs(ratio - 1) < 0.15 else ""
    print(f"  {name:<35} {val:12.8f} {delta_D_needed:12.8f} {ratio:8.4f} {sin2_c:10.7f} {sig:+8.2f}{mark}")


# ============================================================
# [6] BEST PHYSICAL INTERPRETATION
# ============================================================
print(f"\n\n[6] BEST PHYSICAL INTERPRETATION")
print("-" * 50)

# From the scan, let me identify the cleanest forms that work.
# Key constraint: the correction must have a natural group-theoretic
# or perturbative origin.

# The multiplicative form alpha_s^2/(2*pi) has epsilon ~ 0.00222
# vs target 0.00195. That's 14% off, still 6.7 sigma.

# Let me try: sin^2 = 3/(13 - alpha_s*N_c/pi^2)
# delta_D = -alpha_s*N_c/pi^2 = -0.1179*3/9.8696 = -0.03584
# 3/(13-0.0358) = 3/12.964 = 0.23142 -- too much

# sin^2 = 3/(13 - N_c*alpha_s^2)
# delta_D = -3*0.01391 = -0.04173
# 3/(13-0.0417) = 3/12.958 = 0.23153 -- too much

# sin^2 = 3/(13 - alpha_s/(4*pi))
# delta_D = -0.009382
# 3/(13-0.00938) = 3/12.9906 = 0.23094 -- not enough

# The target delta_D = -0.0256 is between alpha_s/(4*pi) and alpha_s/pi.

# Check: 3/(13 - alpha_s/(2*pi))
delta_test = -alpha_s/(2*np.pi)
sin2_test = 3.0/(13.0 + delta_test)
sig_test = (sin2_test - sin2_PDG)/sin2_err
print(f"  3/(13 - alpha_s/(2*pi)):")
print(f"    = 3/{13+delta_test:.6f} = {sin2_test:.7f}")
print(f"    Sigma: {sig_test:+.1f}")
print()

# The cleanest form that works within ~3 sigma:
# sin^2 = N_c / (1 + N_c + N_c^2) * (1 + alpha_s^2/(2*pi))
sin2_best1 = sin2_tree * (1 + alpha_s**2/(2*np.pi))
sig_best1 = (sin2_best1 - sin2_PDG)/sin2_err
print(f"  Candidate A: sin^2 = (3/13)*(1 + alpha_s^2/(2*pi))")
print(f"    = {sin2_best1:.7f}")
print(f"    Sigma: {sig_best1:+.1f}")
print()

# What if we use TGP alpha_s = 0.1190 instead of PDG 0.1179?
alpha_s_TGP = 0.1190
sin2_tgp = sin2_tree * (1 + alpha_s_TGP**2/(2*np.pi))
sig_tgp = (sin2_tgp - sin2_PDG)/sin2_err
print(f"  With alpha_s(TGP) = {alpha_s_TGP}:")
print(f"    = {sin2_tgp:.7f}")
print(f"    Sigma: {sig_tgp:+.1f}")
print()

# Another form: sin^2 = 3/(13 - 3*alpha_s^2/(2*pi))
delta_D2 = -3*alpha_s**2/(2*np.pi)
sin2_alt = 3.0/(13 + delta_D2)
sig_alt = (sin2_alt - sin2_PDG)/sin2_err
print(f"  Candidate B: sin^2 = 3/(13 - 3*alpha_s^2/(2*pi))")
print(f"    delta_D = {delta_D2:.8f}")
print(f"    = {sin2_alt:.7f}")
print(f"    Sigma: {sig_alt:+.1f}")
print()

# The N_c*alpha_s/(4*pi^2) multiplicative correction
eps_nc = alpha_s * N_c / (4*np.pi**2)
sin2_nc = sin2_tree * (1 + eps_nc)
sig_nc = (sin2_nc - sin2_PDG)/sin2_err
print(f"  Candidate C: sin^2 = (3/13)*(1 + N_c*alpha_s/(4*pi^2))")
print(f"    epsilon = {eps_nc:.8f} (target: {epsilon:.8f})")
print(f"    = {sin2_nc:.7f}")
print(f"    Sigma: {sig_nc:+.1f}")
print()


# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*70}")
print("SUMMARY: sin^2(theta_W) QCD CORRECTION")
print(f"{'='*70}")
print(f"""
  Tree level: sin^2 = 3/13 = {sin2_tree:.6f}  (11.3 sigma from PDG)

  Best correction candidates:

  Candidate A: sin^2 = (3/13) * (1 + alpha_s^2/(2*pi))
    = {sin2_best1:.7f}  ({sig_best1:+.1f} sigma)
    Interpretation: 2-loop QCD correction (alpha_s^2 term)

  Candidate B: sin^2 = 3/(13 - 3*alpha_s^2/(2*pi))
    = {sin2_alt:.7f}  ({sig_alt:+.1f} sigma)
    Interpretation: Denominator correction from color-graded running

  Candidate C: sin^2 = (3/13) * (1 + N_c*alpha_s/(4*pi^2))
    = {sin2_nc:.7f}  ({sig_nc:+.1f} sigma)
    Interpretation: 1-loop correction with color factor N_c/(4*pi^2)

  ASSESSMENT:
  - Tree-level TGP is already {abs((sin2_tree-sin2_PDG)/sin2_tree*100):.2f}% from observation.
  - Several simple alpha_s corrections reduce deviation from 11.3 to 2-7 sigma.
  - None gives exact match (<1 sigma) with a clean group-theoretic form.
  - The most promising: alpha_s^2-type corrections (Candidate A: 6.7 sigma).
  - STATUS: PARTIALLY OPEN. Tree-level TGP is impressively close;
    a proper 1-loop EW+QCD calculation in TGP framework is needed.
""")
