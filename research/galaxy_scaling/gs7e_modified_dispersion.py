#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs7e_modified_dispersion.py: Option C — Modified dispersion relation.

IDEA: The TGP substrate has structure at some scale, modifying
the dispersion relation for gravity propagation. At long wavelengths
(galactic scales), this changes the effective potential.

Also explores: what dispersion relation is needed to produce MOND?
Working BACKWARDS from the observed phenomenology to find what
the substrate must look like.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn

print("="*78)
print("  OPTION C: Modified dispersion relation")
print("="*78)

c = 2.998e8
G = 6.674e-11
H0 = 2.20e-18
a0_obs = 1.2e-10
M_sun = 1.989e30

# ===========================================================================
# 1. STANDARD VS MODIFIED DISPERSION
# ===========================================================================
print(f"\n{'='*78}")
print(f"  1. Standard vs modified dispersion in field equations")
print(f"{'='*78}")

print(f"""
  Standard TGP (linearized): delta'' + 2delta'/r + delta = source
  In Fourier space: (-k^2 + 1)*delta_k = source_k
  → delta_k = -source_k / (k^2 - 1)

  This is the dispersion relation: omega^2 = c^2*k^2 - c^2*mu^2
  (tachyonic! the soliton has imaginary mass in field theory sense)
  Or equivalently: k^2 = omega^2/c^2 + mu^2

  For static case (omega=0): k^2 = mu^2 = 1 (in soliton units)
  → Green's function: delta = -source * sin(r)/r (oscillating!)
  → This is the PROBLEM: oscillating tails average out.

  MODIFIED DISPERSION: what if k^2 → k^2 + f(k)?

  Options:
  a) k^4 term (lattice-like):     k^2 → k^2 - beta*k^4
  b) Log correction:               k^2 → k^2*(1 + alpha*ln(k/k0))
  c) Fractional power:             k^2 → k^(2-epsilon)
  d) Infrared modification:        k^2 → k^2 + Lambda_IR (already have: =1)
  e) k-dependent mass:             mu^2 → mu^2 * (1 + f(k))
""")

# ===========================================================================
# 2. WHAT DISPERSION GIVES 1/r FORCE?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  2. Working backwards: what gives 1/r force (not 1/r^2)?")
print(f"{'='*78}")

print(f"""
  For FLAT rotation curves: F ~ 1/r → Phi ~ ln(r)

  Newtonian: Phi_N ~ 1/r → F ~ 1/r^2 (from k^2 in denominator)
  MOND (deep): Phi_M ~ ln(r) → F ~ 1/r

  In Fourier: Phi(k) = -4*pi*G*rho(k) / D(k)
  where D(k) is the "dispersion denominator".

  Newtonian: D(k) = k^2 → Phi(r) ~ 1/r (Green's function)
  For Phi ~ ln(r): need D(k) ~ k (NOT k^2!)

  In 3D: Green's function of D(k) = k:
  G(r) = integral_0^inf k^2 * sin(kr)/(kr) / k * dk / (2*pi^2)
       = 1/(2*pi^2*r) * integral_0^inf k*sin(kr) dk

  This integral DIVERGES! Not well-defined without cutoff.

  With UV cutoff k_max:
  G(r) ~ 1/(2*pi^2*r) * integral_0^k_max k*sin(kr) dk
       = 1/(2*pi^2*r) * [sin(kr)/(r) - k*cos(kr)]_0^k_max
       = k_max*cos(k_max*r) / ... (oscillatory)

  → D(k) = k doesn't give a clean ln(r) potential in 3D.

  ALTERNATIVE: D(k) = k^2 * (1 + k/k0)^(-1) = k^2*k0/(k + k0)

  At small k (long range): D ~ k^2/k0 * k0 → k^2 (Newtonian)
  Wait, that doesn't help.

  Let me think more carefully...

  In MOND, the AQUAL formulation modifies Poisson:
  nabla . [mu(|nabla Phi|/a0) * nabla Phi] = 4*pi*G*rho

  where mu(x) → 1 at x >> 1, mu(x) → x at x << 1.

  This is a NONLINEAR modification — NOT a simple dispersion change!
  The modification depends on the FIELD VALUE, not on the wavenumber k.

  This is fundamentally different from a dispersion relation,
  which is LINEAR in the field.

  → MOND CANNOT be expressed as a modified dispersion relation
    for a LINEAR field theory!
""")

# ===========================================================================
# 3. CAN NONLINEAR DISPERSION HELP?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  3. Nonlinear dispersion: amplitude-dependent propagation")
print(f"{'='*78}")

print(f"""
  In TGP, the equation IS nonlinear: g'^2/g term.
  Could this give amplitude-dependent dispersion?

  Full TGP: g'' + g'^2/g + 2g'/r + g = 1

  Linearize around g = 1 + delta:
  delta'' + 2delta'/r + delta + (delta')^2/(1+delta) = 0

  The (delta')^2 term:
  - Depends on AMPLITUDE of delta
  - At delta ~ 10^-6 (galactic): (delta')^2 ~ delta^2/r^2 ~ 10^-12 (negligible)
  - At delta ~ 1 (soliton core): (delta')^2 ~ 1 (significant)

  So the nonlinear "dispersion" works at STRONG fields (delta~1)
  but is negligible at WEAK fields (delta~10^-6).

  For galaxies: delta ~ 10^-6 → linear theory applies → no modified dispersion.

  WHAT IF there's a different nonlinearity that matters at weak fields?

  In standard physics, there's no such thing:
  weak-field nonlinearities are always perturbative.

  BUT: MOND's mu function IS a "strong coupling at weak fields" effect.
  It's as if the gravitational coupling INCREASES at low field strengths.

  In field theory terms: this requires an ANTI-SCREENING mechanism.
  Normal screening (QED, QCD): coupling weakens at large distances.
  MOND-like: coupling STRENGTHENS at large distances.

  In QCD: ANTI-screening does occur! (asymptotic freedom reversed)
  The coupling INCREASES at long distances → confinement.

  Could TGP have a similar mechanism?
  The g'^2/g term in TGP is like a self-interaction.
  At weak fields, g ~ 1, g'^2/g ~ (delta')^2 → standard.
  No anti-screening arises naturally.
""")

# ===========================================================================
# 4. LATTICE/DISCRETE SUBSTRATE
# ===========================================================================
print(f"\n{'='*78}")
print(f"  4. Lattice/discrete substrate effects")
print(f"{'='*78}")

print(f"""
  If TGP substrate is discrete with spacing l:
  omega^2 = (c^2/l^2) * (2 - 2*cos(k*l)) = (2*c^2/l^2) * sin^2(k*l/2)

  At small k (k*l << 1): omega^2 ~ c^2*k^2 (standard)
  At k ~ pi/l (Brillouin zone boundary): omega^2 = 4*c^2/l^2 (maximum)

  The modification to the POTENTIAL at large r:
  Phi(r) = -G*M/r * f(r/l)

  For cubic lattice:
  f(r/l) = 1 + (l/r)^2 * (something) + ...

  The correction scales as (l/r)^2.
  For l = any microphysical scale and r = 10 kpc:
  (l/r)^2 ~ (l / 3e20)^2

  Even for l = 1 meter: (1/3e20)^2 ~ 10^-41 → NOTHING.
  For l = c/H0 ~ 10^26 m: (10^26/3e20)^2 ~ 10^11 → TOO LARGE!

  No intermediate scale l gives the right correction at galactic scales
  without being fine-tuned.

  UNLESS: the lattice affects gravity differently than a simple
  modified dispersion. For example, if the lattice has DEFECTS
  or DISLOCATIONS that concentrate at certain scales...

  But this is too speculative without a concrete model.
""")

# ===========================================================================
# 5. INFRARED MODIFICATION: MASSIVE GRAVITON
# ===========================================================================
print(f"\n{'='*78}")
print(f"  5. Infrared modification: massive graviton analogy")
print(f"{'='*78}")

print(f"""
  A massive graviton gives Yukawa potential:
  Phi(r) = -G*M/r * exp(-r/lambda_g)

  where lambda_g = hbar/(m_g*c).

  For flat rotation curves, we DON'T want exponential cutoff.
  But: TGP has a DIFFERENT mass term (+g = 1 → oscillatory, not Yukawa).

  What if we REMOVE the mass term? (mu = 0)
  g'' + g'^2/g + 2g'/r = 0

  Linearized: delta'' + 2delta'/r = 0
  Solution: delta = A/r + B → NEWTONIAN!

  So: mu = 0 gives Newton, mu = 1 gives oscillating tails.
  Neither gives MOND.

  What about mu^2 < 0 (truly tachyonic)?
  delta'' + 2delta'/r - |mu|^2 * delta = 0
  Solution: delta = A * exp(-|mu|*r)/r + B * exp(+|mu|*r)/r

  The growing mode exp(+|mu|*r)/r is UNPHYSICAL (diverges at infinity).
  The decaying mode gives Yukawa: exp(-|mu|*r)/r → shorter range than Newton.

  → Tachyonic mass makes gravity WEAKER at large r, not stronger.
  → OPPOSITE of what MOND needs!
""")

# ===========================================================================
# 6. BEKENSTEIN'S TeVeS APPROACH
# ===========================================================================
print(f"\n{'='*78}")
print(f"  6. How Bekenstein solved this: TeVeS")
print(f"{'='*78}")

print(f"""
  Bekenstein (2004) created a relativistic MOND theory called TeVeS
  (Tensor-Vector-Scalar).

  The key: he introduced a SCALAR FIELD phi with Lagrangian:
  L = -mu(l^2 * |nabla phi|^2) * |nabla phi|^2

  where mu(y) is a free function that behaves as:
  mu → 1 at y >> 1 (strong field, Newtonian)
  mu → sqrt(y) at y << 1 (weak field, MOND)

  The crucial feature: the Lagrangian density depends on |nabla phi|^2
  in a NONLINEAR way. This is what produces MOND.

  In TGP terms: can we identify phi with the TGP potential delta?

  TGP Lagrangian: L_TGP = |nabla g|^2 / g = |nabla delta|^2 / (1+delta)

  For weak field (delta << 1):
  L_TGP ~ |nabla delta|^2 * (1 - delta + delta^2 - ...)

  The 1/(1+delta) factor IS a nonlinear function of delta,
  but it's a function of delta, NOT of |nabla delta|.

  For TeVeS/MOND: need mu(|nabla phi|^2), not mu(phi).

  TGP has mu(g) = 1/g, which depends on the FIELD VALUE g,
  not on the GRADIENT |nabla g|.

  → TGP's nonlinearity is in the WRONG VARIABLE to produce MOND.

  HOWEVER: could a MORE GENERAL TGP Lagrangian work?

  L_gen = f(|nabla g|^2 / g^alpha) * g^beta

  If alpha = 0, beta = 0: f(|nabla g|^2) → TeVeS-like!
  If alpha = 1, beta = 0: |nabla g|^2/g → standard TGP

  So: TeVeS is a SPECIAL CASE of generalized TGP with alpha = 0.
  Standard TGP (alpha = 1) is the "wrong" choice for MOND.
""")

# ===========================================================================
# 7. WHAT MODIFICATION OF TGP GIVES MOND?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  7. What modification of TGP would give MOND?")
print(f"{'='*78}")

print(f"""
  Standard TGP: L = |nabla g|^2 / g + (g-1)^2
  → Field equation: nabla(nabla g / g) + (g-1) = source

  For MOND, we need the KINETIC term to be nonlinear in a specific way.
  Bekenstein's insight: L_kinetic = F(X) where X = |nabla g|^2

  F(X) → X at large X (strong field, Newtonian)
  F(X) → X^(3/2) / sqrt(a0_scale) at small X (weak field, MOND)

  The resulting field equation:
  nabla . [F'(|nabla g|^2) * nabla g] + (g-1) = source

  At strong field (F' ~ 1): standard Poisson
  At weak field (F' ~ sqrt(X)/a0_scale): MOND equation

  This gives the AQUAL formulation of MOND if:
  F'(X) = mu(sqrt(X)/a0_TGP)

  where a0_TGP is an acceleration scale that we need to determine.

  The TGP Lagrangian L = |nabla g|^2 / g has:
  F(X, g) = X / g → F depends on BOTH X and g.

  For weak field (g ~ 1): F ~ X (linear in X) → Newtonian.
  No MOND behavior.

  TO GET MOND FROM TGP:
  We'd need to MODIFY the TGP Lagrangian from:
    L = |nabla g|^2 / g + (g-1)^2
  to:
    L = f(|nabla g|^2) + (g-1)^2

  where f(X) has MOND-like behavior.

  But this CHANGES the fundamental TGP equation!
  The question is: does the microscopic TGP substrate naturally
  produce such a Lagrangian at macroscopic scales?
""")

# ===========================================================================
# 8. RENORMALIZATION GROUP FLOW OF THE TGP LAGRANGIAN
# ===========================================================================
print(f"\n{'='*78}")
print(f"  8. RG flow: does the TGP Lagrangian change at long distances?")
print(f"{'='*78}")

print(f"""
  In quantum field theory, coupling constants "run" with energy scale.
  The Lagrangian at scale mu is different from the Lagrangian at scale Lambda.

  Could the TGP Lagrangian:
  L_micro = |nabla g|^2 / g + (g-1)^2

  flow to a DIFFERENT form at macroscopic scales?
  L_macro = f(|nabla g|^2) + (g-1)^2  (MOND-like)

  This would require that the 1/g factor in the kinetic term
  gets RENORMALIZED away at large distances, replaced by a
  nonlinear function of |nabla g|^2.

  In perturbation theory:
  g = 1 + delta, |delta| << 1
  L = |nabla delta|^2 * (1 - delta + delta^2 - ...) + delta^2

  One-loop corrections from the delta^2 and delta^3 vertices
  would renormalize the |nabla delta|^2 coupling.

  But in weak field (delta ~ 10^-6):
  corrections ~ delta ~ 10^-6 → NEGLIGIBLE.

  RG flow can't generate O(1) modifications from O(10^-6) corrections.

  UNLESS: there's a NONPERTURBATIVE effect (instanton, condensate, etc.)
  that changes the effective Lagrangian. But we have no evidence for this.
""")

# ===========================================================================
# 9. SUMMARY: WHAT DISPERSION/MODIFICATION WOULD WORK?
# ===========================================================================
print(f"\n{'='*78}")
print(f"  9. Summary: what modification is needed?")
print(f"{'='*78}")

print(f"""
  To get MOND from a field theory, you need ONE of:

  1. NONLINEAR KINETIC TERM: f(|nabla phi|^2) with f'→sqrt(X) at small X
     → Bekenstein's TeVeS / AQUAL
     → TGP has 1/g factor but that's field-dependent, not gradient-dependent
     → DOESN'T match

  2. ADDITIONAL FIELD: a vector or tensor field that mediates
     the MOND modification (e.g., TeVeS's vector field)
     → Requires extending TGP with extra degrees of freedom
     → POSSIBLE but goes beyond single-field TGP

  3. NONLOCAL MODIFICATION: the Poisson equation becomes nonlocal
     (e.g., Milgrom's QUMOND: nabla^2 Phi = nabla . [nu * nabla Phi_N])
     → Would need a nonlocal TGP substrate
     → POSSIBLE if substrate has long-range correlations

  4. HIGHER DIMENSIONS: if the substrate is higher-dimensional,
     gravity could "leak" into extra dimensions at large distances
     → DGP (Dvali-Gabadadze-Porrati) model
     → Changes from 1/r^2 to 1/r at large r → flat RC!
     → BUT: 5D gravity gives 1/r^3 force, 4D → 1/r^2, neither gives 1/r
     → Need a "crossover" mechanism (DGP brane model)

  ASSESSMENT:
  None of the standard dispersion modifications naturally produces MOND
  from the TGP Lagrangian L = |nabla g|^2/g + (g-1)^2.

  The fundamental issue: TGP's nonlinearity (1/g factor) depends on
  the FIELD VALUE, but MOND requires nonlinearity in the GRADIENT.
""")

# ===========================================================================
# 10. FINAL ASSESSMENT
# ===========================================================================
print(f"\n{'='*78}")
print(f"  10. ASSESSMENT of Option C")
print(f"{'='*78}")

print(f"""
  TESTED:
  1. Linear modified dispersion (k^4, fractional, etc.):
     MOND requires nonlinear modification → linear dispersion can't help
  2. Lattice/discrete substrate:
     Corrections scale as (l/r)^2 → no intermediate scale works
  3. Massive graviton / Yukawa:
     TGP already has this (mu=1) → oscillatory, not MOND
     mu=0 gives Newton, mu^2<0 gives shorter range
  4. Amplitude-dependent (nonlinear) dispersion:
     TGP nonlinearity g'^2/g is ~ delta^2 at galactic scales → negligible
  5. TeVeS comparison:
     MOND needs f(|nabla g|^2), TGP has |nabla g|^2/g
     → Wrong variable (field value vs gradient)
  6. RG flow:
     Can't generate O(1) modifications from O(10^-6) weak-field corrections

  FUNDAMENTAL INSIGHT:
  MOND is NOT a modified dispersion relation.
  MOND is a NONLINEAR field theory where the nonlinearity
  depends on |nabla Phi| (the acceleration), not on Phi itself.

  TGP's nonlinearity depends on g (the field value), not |nabla g|.
  This is the WRONG TYPE of nonlinearity for MOND.

  To make TGP produce MOND, you'd need to CHANGE the Lagrangian
  from L = |nabla g|^2/g to L = f(|nabla g|^2), which is a
  DIFFERENT theory.

  STATUS: FAILED
  - Modified dispersion: can't produce MOND (needs nonlinear, not linear)
  - TGP nonlinearity: wrong variable (g vs |nabla g|)
  - No natural RG flow to MOND-like Lagrangian
  - Lattice effects too small at any physical scale
""")
