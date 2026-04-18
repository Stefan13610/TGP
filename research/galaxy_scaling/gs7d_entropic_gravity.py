#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs7d_entropic_gravity.py: Option B — Entropic gravity (Verlinde-like) in TGP.

IDEA: Gravity in TGP emerges from substrate entropy.
The substrate volume sqrt(g)*d^3x represents degrees of freedom.
An entropic force F = T * dS/dr could have different r-dependence.

Verlinde (2016) proposed that in de Sitter space, there's an
elastic component of gravity from entanglement entropy that gives
MOND-like behavior. Can TGP provide a concrete realization?
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np

print("="*78)
print("  OPTION B: Entropic/thermodynamic gravity in TGP substrate")
print("="*78)

c = 2.998e8
G = 6.674e-11
hbar = 1.055e-34
k_B = 1.381e-23
H0 = 2.20e-18
a0_obs = 1.2e-10
M_sun = 1.989e30

R_H = c / H0
T_dS = hbar * H0 / (2*np.pi*k_B)  # de Sitter temperature
T_Unruh_a0 = hbar * a0_obs / (2*np.pi*k_B*c)  # Unruh temperature at a0

# ===========================================================================
# 1. VERLINDE'S ARGUMENT (2016)
# ===========================================================================
print(f"\n{'='*78}")
print(f"  1. Verlinde's emergent gravity argument")
print(f"{'='*78}")

print(f"""
  Verlinde's key idea (arXiv:1611.02269):

  In de Sitter space, there's a THERMAL contribution to gravity
  from the entanglement entropy of the vacuum.

  The total entropy of a Hubble-sized region:
  S_total = A_H / (4*l_P^2) = 4*pi*R_H^2 / (4*l_P^2)
          = pi * R_H^2 / l_P^2

  R_H = c/H0 = {R_H:.3e} m
  l_P = {1.616e-35:.3e} m
  S_total = pi * ({R_H:.1e})^2 / ({1.616e-35:.1e})^2
          = {np.pi * R_H**2 / (1.616e-35)**2:.3e} (in units of k_B)

  A mass M at the center "uses up" some of this entropy for its
  gravitational field. The entropy DEFICIT creates an elastic
  (entropic) response.

  Verlinde's formula for the entropic contribution:
  M_D(r) = (c*H0)/(6*G) * r^2 * g_B(r) / a0   [apparent dark matter]

  where g_B = GM/r^2 is the baryonic acceleration.

  Substituting: M_D(r) = (c*H0*M)/(6*a0) = const (independent of r!)

  Wait — that gives a constant M_D, which means:
  g_total = GM/r^2 + G*M_D/r^2 = G(M + M_D)/r^2

  This is STILL 1/r^2 — NOT flat rotation curve!

  Let me recheck Verlinde's actual formula...

  Verlinde's equation (7.40): the apparent DM mass profile
  M_D(r) satisfies:
    M_D(r)^2 = (c*H0*r^2/(6*G)) * M_B(r)    [for spherical, isolated system]

  where M_B(r) = baryonic mass enclosed.

  This gives: M_D(r) = sqrt(c*H0*r^2*M_B/(6G))
  g_D = G*M_D/r^2 = sqrt(c*H0*G*M_B/(6*r^2))
  g_D = sqrt(c*H0*g_B/6)

  In deep MOND regime: g_obs = sqrt(g_B * a0)
  Comparing: a0_Verlinde = c*H0/6

  a0_Verlinde = c*H0/6 = {c*H0/6:.3e} m/s^2
  a0_obs = {a0_obs:.1e} m/s^2
  ratio = {c*H0/6/a0_obs:.3f}

  So Verlinde predicts a0 = cH0/6 (factor of 6, not 2*pi)!
  This is {(c*H0/6/a0_obs - 1)*100:.0f}% off from observed.

  Compare: cH0/(2*pi) = {c*H0/(2*np.pi):.3e} → ratio {c*H0/(2*np.pi)/a0_obs:.3f}
  vs: cH0/6 = {c*H0/6:.3e} → ratio {c*H0/6/a0_obs:.3f}

  Both are close! The factor is between 2*pi ~ 6.28 and 6.
""")

# ===========================================================================
# 2. TGP ENTROPY: WHAT IS THE NATURAL MEASURE?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  2. TGP substrate entropy — what is the natural measure?")
print(f"{'='*78}")

print(f"""
  In TGP, g is the substrate metric (conformal factor).
  The substrate "volume element" is sqrt(g) * d^3x.

  Candidate entropy functionals:

  S1 = integral sqrt(g) d^3x / l_P^3      (volume in Planck units)
  S2 = integral (g - 1)^2 d^3x / l_P^3    (elastic strain energy)
  S3 = integral |nabla g|^2 / g d^3x / l_P  (kinetic/gradient entropy)
  S4 = integral ln(g) d^3x / l_P^3         (information entropy)

  For a soliton with g = 1 + delta, delta << 1:
  S1 ~ integral (1 + delta/2) d^3x  → V + delta*V/2  (volume + correction)
  S2 ~ integral delta^2 d^3x        → elastic energy
  S3 ~ integral (nabla delta)^2 d^3x → gradient energy (surface tension)
  S4 ~ integral delta d^3x          → linear in delta

  For entropic gravity, we need:
  F_entropic = T * dS/dr

  For a mass at origin: delta(r) ~ -GM/(rc^2)
  dS/dr depends on which S we choose.
""")

# Compute for each entropy functional
print(f"  Entropic force for each functional:")
print(f"  (T = T_dS = {T_dS:.3e} K)")
print()

# S1: volume entropy
print(f"  S1 (volume): S = integral sqrt(g) d^3x")
print(f"    delta(r) = -GM/(rc^2)")
print(f"    sqrt(g) = sqrt(1+delta) ~ 1 + delta/2")
print(f"    S(R) = 4*pi*integral_0^R (1 + delta/2) r^2 dr")
print(f"         = (4*pi*R^3/3) + 2*pi*integral delta*r^2 dr")
print(f"    For delta = -GM/(rc^2):")
print(f"    integral (-GM/(rc^2))*r^2 dr = -GM*R^2/(2c^2)")
print(f"    S(R) = 4*pi*R^3/3 - pi*GM*R^2/c^2")
print(f"    dS/dR = 4*pi*R^2 - 2*pi*GM*R/c^2")
print(f"    F = T*dS/dR = T*(4*pi*R^2 - 2*pi*GM*R/c^2)")
print(f"    The first term (4*pi*R^2) is a PRESSURE term (cosmological)")
print(f"    The second term (-2*pi*GM*R/c^2*T) is the gravitational part")
print()
F_grav_S1 = T_dS * 2*np.pi*G*7e10*M_sun*10*3.086e19/c**2
print(f"    At r=10 kpc, MW: F_grav_S1/m = {F_grav_S1:.3e} m/s^2")
print(f"    vs F_Newton = {G*7e10*M_sun/(10*3.086e19)**2:.3e} m/s^2")
print(f"    → ratio = {F_grav_S1/(G*7e10*M_sun/(10*3.086e19)**2):.3e}")
print(f"    → NEGLIGIBLE (T_dS is incredibly small: {T_dS:.1e} K)")

# ===========================================================================
# 3. THE TEMPERATURE PROBLEM
# ===========================================================================
print(f"\n{'='*78}")
print(f"  3. The temperature problem")
print(f"{'='*78}")

print(f"""
  The de Sitter temperature T_dS = hbar*H0/(2*pi*k_B) = {T_dS:.3e} K

  This is so small that ANY entropic force F = T*dS/dr is negligible
  compared to Newton, unless S is ENORMOUS.

  Required S for F_entropic ~ F_Newton = GM/r^2:
  T * dS/dr = GM/r^2
  dS/dr = GM/(r^2 * T) = GM * 2*pi*k_B / (r^2 * hbar * H0)

  For MW at r = 10 kpc = {10*3.086e19:.1e} m:
""")

r_10kpc = 10 * 3.086e19
dSdr_needed = G*7e10*M_sun/(r_10kpc**2 * T_dS)
print(f"    dS/dr needed = {dSdr_needed:.3e} /m")
print(f"    S(10 kpc) ~ dS/dr * r = {dSdr_needed * r_10kpc:.3e}")
print(f"    Bekenstein-Hawking entropy of observable universe:")

S_BH = np.pi * R_H**2 / (1.616e-35)**2
print(f"    S_BH = pi*R_H^2/l_P^2 = {S_BH:.3e}")
print(f"    Ratio S_needed/S_BH = {dSdr_needed * r_10kpc / S_BH:.3e}")

print(f"""
  The entropy GRADIENT needed is comparable to the total entropy
  of the observable universe per meter — IMPOSSIBLE from a single galaxy!

  VERLINDE'S TRICK:
  Verlinde doesn't use T_dS directly. He uses a more subtle argument:
  The elastic strain in the entanglement entropy of de Sitter space.

  His key equation: S_elastic(r) = (c^3/(6*G*hbar*H0)) * M * a0 * r^2

  This gives: F = T * dS/dr = (hbar*H0/(2*pi*k_B)) * (c^3*M*a0*2r)/(6*G*hbar*H0)
             = c^3*M*a0*r / (6*pi*G*k_B)

  Wait, this doesn't simplify to the right thing.
  Let me redo Verlinde's derivation properly.
""")

# ===========================================================================
# 4. VERLINDE'S ACTUAL DERIVATION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  4. Verlinde's derivation — the key steps")
print(f"{'='*78}")

print(f"""
  Step 1: Bekenstein's entropy-area relation for a mass M:
    S_M = 2*pi*k_B*M*c / hbar * lambda_C  (entropy within Compton wavelength)
    But this is for black holes. For normal matter:
    S_M = 2*pi*k_B*M*c*R / hbar  (entropy within volume of size R)

  Step 2: In de Sitter space, the TOTAL entropy is:
    S_dS = pi*c^3 / (G*hbar*H0^2) = {np.pi*c**3/(G*hbar*H0**2):.3e} k_B

  Step 3: A mass M at center creates a "gravitational entropy" that
    reduces the available entropy:
    Delta_S = - (2*pi*M*c)/(hbar) * integral_0^R Phi(r)/c^2 * 4pi*r^2 dr / ...

  Actually, let me use Verlinde's KEY RESULT directly:

  The apparent dark matter distribution that arises from the elastic
  entropy response of de Sitter space:

    M_D(r)^2 = (c^2*H0^2*r^4)/(6*G*R_H^2) * integral_0^r 4pi*r'^2 * rho_B(r') dr'
             = wait... his formula is simpler.

  Verlinde's equation (7.40) for spherical symmetry:

    g_D(r) * g_B(r) = (c*H0)^2 / 6     [acceleration space]

  where g_D = apparent dark matter acceleration,
        g_B = baryonic acceleration.

  This means: g_D = (cH0)^2 / (6 * g_B)

  Total: g_total = g_B + g_D = g_B + (cH0)^2/(6*g_B)

  In the deep MOND regime (g_D >> g_B):
  g_total ~ g_D = (cH0)^2 / (6*g_B) → g_total^2 ~ (cH0)^2/6 * g_N/g_total
  Wait, this doesn't work simply.

  Let me just use: g_total = g_B + (cH0)^2/(6*g_B)

  At g_B >> cH0/sqrt(6): g_total ~ g_B (Newtonian)
  At g_B << cH0/sqrt(6): g_total ~ (cH0)^2/(6*g_B) (anti-MOND?!)

  Hmm, this grows as g_B decreases — that's the WRONG behavior.
  It means the DM acceleration INCREASES when baryonic decreases.

  Actually for MOND: g_obs = sqrt(g_B * a0) when g_B << a0
  So g_obs^2 = g_B * a0

  In Verlinde: g_total ~ g_B + (cH0)^2/(6*g_B)
  For total dominated by the second term: g_total ~ (cH0)^2/(6*g_B)
  Then: g_total * g_B = (cH0)^2/6

  But MOND says: g_obs * g_B = g_obs^2 when g_obs = sqrt(g_B*a0)
  → g_obs * g_B = g_B^(3/2) * a0^(1/2) ≠ const

  So Verlinde's formula g_D*g_B = const is NOT the same as MOND!

  Let me check: in MOND, g_obs = g_B/(1 - exp(-sqrt(g_B/a0)))
  Defining g_D = g_obs - g_B:
""")

# Compare Verlinde vs MOND
print(f"  Comparing Verlinde vs MOND:")
print(f"  {'g_B/a0':>8s} {'g_MOND/a0':>10s} {'g_D_MOND/a0':>12s} {'g_D_Verl/a0':>12s} {'MOND_product':>14s} {'Verl_product':>14s}")
print(f"  {'-'*72}")

a0_verl = (c*H0)**2 / (6 * a0_obs)  # Verlinde's effective a0-like quantity

for log_gB in np.arange(-2, 3, 0.5):
    gB = 10**log_gB * a0_obs
    # MOND
    g_MOND = gB / (1 - np.exp(-np.sqrt(gB/a0_obs)))
    g_D_MOND = g_MOND - gB
    # Verlinde: g_D = (cH0)^2 / (6*gB)
    g_D_Verl = (c*H0)**2 / (6*gB)

    product_MOND = g_D_MOND * gB / a0_obs**2
    product_Verl = g_D_Verl * gB / a0_obs**2

    print(f"  {gB/a0_obs:8.2f} {g_MOND/a0_obs:10.3f} {g_D_MOND/a0_obs:12.3f} {g_D_Verl/a0_obs:12.3f} {product_MOND:14.3f} {product_Verl:14.3f}")

print(f"""
  Verlinde's product g_D*g_B/(a0^2) = (cH0)^2/(6*a0^2) = {(c*H0)**2/(6*a0_obs**2):.3f}
  (should be ~ const ~ 1 for MOND-like behavior)

  MOND product g_D*g_B/(a0^2) → 1.0 at g_B << a0 (deep MOND)
  Verlinde product = {(c*H0)**2/(6*a0_obs**2):.1f} (constant everywhere)

  For Verlinde to match MOND: need (cH0)^2/6 = a0^2
  → a0 = cH0/sqrt(6) = {c*H0/np.sqrt(6):.3e} = {c*H0/np.sqrt(6)/a0_obs:.2f} * a0_obs

  So Verlinde predicts a0 = cH0/sqrt(6) = {c*H0/np.sqrt(6):.3e} m/s^2
  Ratio: {c*H0/np.sqrt(6)/a0_obs:.3f}

  This is {abs(1-c*H0/np.sqrt(6)/a0_obs)*100:.0f}% too high!
  But the FUNCTIONAL FORM is wrong too: g_D*g_B = const ≠ MOND.
""")

# ===========================================================================
# 5. TGP-SPECIFIC ENTROPY: SUBSTRATE DEFORMATION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  5. TGP-specific entropy: substrate deformation")
print(f"{'='*78}")

print(f"""
  In TGP, the substrate has a metric g(x).
  The natural "entropy" or "energy" functional is the TGP Lagrangian:

  L_TGP = integral [ |nabla g|^2 / g + (g - 1)^2 ] d^3x

  First term: gradient energy (kinetic/elastic)
  Second term: potential energy (spring)

  For a galaxy: g = 1 + delta, |delta| ~ GM/(rc^2) ~ 10^-6

  L_gradient ~ integral (nabla delta)^2 d^3x ~ integral (GM/(r^2*c^2))^2 r^2 dr
             ~ (GM/c^2)^2 * integral dr/r^2 ~ (GM/c^2)^2 / r_min

  L_spring ~ integral delta^2 d^3x ~ integral (GM/(rc^2))^2 r^2 dr
           ~ (GM/c^2)^2 * R_max

  These are just the standard gravitational self-energy.
  No new physics emerges from this.

  WHAT IF entropy is defined differently?

  In TGP, information might be stored in the PHASE of the substrate
  oscillation. The phase of sin(r)/r tails carries information about
  the source. This is similar to holographic information storage.

  But the phase information averages out for N >> 1 solitons (gs6 result).

  ALTERNATIVE: Verlinde's argument doesn't actually use a specific entropy
  functional. It uses the ENTANGLEMENT entropy of the vacuum, which is
  a property of quantum gravity, not of a specific field theory.

  For TGP to give Verlinde-like emergent gravity, we'd need to:
  1. Quantize the TGP substrate (second quantization of g(x))
  2. Compute entanglement entropy across Hubble horizon
  3. Show that mass redistributes this entropy
  4. Derive the entropic force

  This is a HUGE program — well beyond what we can compute here.
""")

# ===========================================================================
# 6. UNRUH EFFECT APPROACH
# ===========================================================================
print(f"\n{'='*78}")
print(f"  6. Unruh effect: minimal acceleration")
print(f"{'='*78}")

print(f"""
  A simpler entropic argument (Milgrom, McCulloch):

  An accelerating observer sees Unruh radiation at temperature:
  T_U = hbar*a / (2*pi*k_B*c)

  The Unruh wavelength: lambda_U = 2*pi*c^2 / a

  When a → a0 = cH0/(2*pi):
  lambda_U(a0) = 2*pi*c^2 / (cH0/(2*pi)) = 4*pi^2*c/H0 = 4*pi^2*R_H
               = {4*np.pi**2*R_H:.3e} m = {4*np.pi**2*R_H/3.086e22:.0f} Mpc

  This is LARGER than the Hubble radius R_H = {R_H:.3e} m = {R_H/3.086e22:.0f} Mpc!

  McCulloch's "Quantized Inertia" (MiHsC):
  If Unruh wavelength can't exceed 2*R_H (Hubble horizon cutoff),
  then the minimum acceleration is:
  a_min = pi*c^2/R_H = pi*c*H0 = {np.pi*c*H0:.3e} m/s^2

  ratio to a0: {np.pi*c*H0/a0_obs:.3f}

  Or with 2*pi*c/R_H: a_min = 2*pi*c^2/(2*R_H) = pi*c*H0
  Same result.

  Actually the standard MiHsC prediction is:
  a_min = 2*c^2/R_H = 2*c*H0 = {2*c*H0:.3e}
  ratio: {2*c*H0/a0_obs:.1f} — too large by 11x

  CORRECT MILGROM: a0 ~ cH0/(2pi) requires that the MINIMUM
  Unruh wavelength that fits in the Hubble sphere is:
  lambda_max = 2*pi*R_H (NOT 2*R_H)

  Then: a_min = 2*pi*c^2/lambda_max = 2*pi*c^2/(2*pi*R_H) = c^2/R_H = c*H0
  Still too large.

  With the right numerical factor:
  a0 = c*H0/(2*pi) requires lambda_max = (2*pi)^2 * c/H0 = 4*pi^2 * R_H

  This is the STANDING WAVE condition: the Unruh wave must fit
  exactly 1/(4*pi) wavelengths in the Hubble sphere. Not natural.
""")

# ===========================================================================
# 7. SUMMARY TABLE
# ===========================================================================
print(f"\n{'='*78}")
print(f"  7. Summary: entropic/thermodynamic predictions for a0")
print(f"{'='*78}")

predictions = [
    ("Verlinde g_D*g_B=const", c*H0/np.sqrt(6), "g_D*g_B=const (≠ MOND)"),
    ("Verlinde (cH0/6)", c*H0/6, "from elastic strain formula"),
    ("McCulloch MiHsC", 2*c*H0, "lambda_U < 2*R_H"),
    ("Unruh + Hubble", np.pi*c*H0, "lambda_U < pi*R_H"),
    ("Simple c*H0", c*H0, "dimensional"),
    ("c*H0/(2*pi)", c*H0/(2*np.pi), "best empirical match"),
    ("a_dS/(2*pi)", c*H0*np.sqrt(0.685)/(2*np.pi), "de Sitter + 2pi"),
]

print(f"\n  {'Prediction':<30s} {'a0 (m/s^2)':>12s} {'ratio':>8s} {'Mechanism':<35s}")
print(f"  {'-'*88}")
for name, val, mech in predictions:
    r = val/a0_obs
    mark = " <-- BEST" if abs(r-1) < 0.15 else (" <--" if abs(r-1) < 0.3 else "")
    print(f"  {name:<30s} {val:12.3e} {r:8.3f} {mech:<35s}{mark}")

# ===========================================================================
# 8. ASSESSMENT
# ===========================================================================
print(f"\n{'='*78}")
print(f"  8. ASSESSMENT of Option B")
print(f"{'='*78}")

print(f"""
  WHAT WORKS:
  + Verlinde's framework DOES give MOND-like behavior (partially)
  + The scale a0 ~ cH0/(few) emerges naturally
  + Multiple independent arguments converge on a0 ~ cH0
  + The connection to de Sitter entropy is compelling

  WHAT DOESN'T WORK:
  - Verlinde's functional form g_D*g_B = const is NOT exactly MOND
    (MOND: g_obs = sqrt(g_B*a0) → g_D*g_B varies with g_B)
  - The numerical factor varies: sqrt(6), 6, 2*pi — no consensus
  - T_dS = {T_dS:.1e} K → direct F=T*dS/dr gives negligible force
  - Requires quantum gravity framework (entanglement entropy)
  - No concrete TGP entropy functional reproduces the result
  - McCulloch's MiHsC is widely criticized as ad hoc

  TGP-SPECIFIC ISSUES:
  - TGP substrate entropy L = integral (nabla g)^2/g + (g-1)^2 d^3x
    gives standard self-energy, no new physics
  - Quantization of TGP substrate not developed
  - Connection between TGP solitons and entanglement entropy unclear

  KEY INSIGHT:
  The entropic argument gives the RIGHT SCALE (a0 ~ cH0)
  but the WRONG FUNCTIONAL FORM (g_D*g_B = const ≠ MOND).
  The exact numerical factor (6 vs 2*pi) depends on the derivation.

  FUNDAMENTALLY: entropic gravity is an INTERPRETATION, not a MECHANISM.
  It says "gravity behaves as if there's an entropic force" but doesn't
  specify the microscopic degrees of freedom.

  For TGP: the substrate IS the microscopic DoF, but we can't compute
  its entanglement entropy without a quantum TGP theory.

  STATUS: INCONCLUSIVE
  - Right scale: YES (a0 ~ cH0)
  - Right functional form: NO (g_D*g_B ≠ MOND)
  - Concrete mechanism: NO (needs quantum TGP)
  - Testable prediction: YES (a0 evolves with H(z) — same as Option D)
""")
