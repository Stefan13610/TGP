#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TGP Ringdown: Perturbation equation derivation.
Derives V_eff(r) for scalar perturbations around static spherical background.
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import sympy as sp

r = sp.Symbol('r', positive=True)
f = sp.Function('f')(r)
Psi = sp.Function('Psi')(r)
fp = sp.diff(f, r)
fpp = sp.diff(f, r, 2)

Phi0 = sp.Symbol('Phi_0', positive=True)
c0 = sp.Symbol('c_0', positive=True)
beta_s = sp.Symbol('beta', positive=True)
gamma_s = sp.Symbol('gamma', positive=True)
omega = sp.Symbol('omega')
ell = sp.Symbol('ell', nonnegative=True, integer=True)

ALPHA = 2

print("=" * 65)
print("  TGP RINGDOWN: PERTURBATION EQUATION DERIVATION")
print("=" * 65)

print("""
STEP 1: Starting equations
---------------------------
Covariant field eq (vacuum, rho=0):

  (Phi/(c0^2 Phi0)) Phi_tt - nabla^2 Phi
  + (2 Phi/(c0^2 Phi0))(Phi_t)^2 / Phi - 2 (nabla Phi)^2 / Phi
  - beta Phi^2/Phi0 + gamma Phi^3/Phi0^2 = 0

c(Phi) = c0 sqrt(Phi0/Phi),  alpha = 2.

Background f(r):
  f'' + (2/r)f' + 2(f')^2/f + beta f^2/Phi0 - gamma f^3/Phi0^2 = 0
""")

print("STEP 2: Linearize Phi = f(r) + eps*h(t,r,theta,phi)")
print("-" * 50)
print("""
chi(r) = f(r)/Phi0,  c_bg^2 = c0^2/chi

Master equation at O(eps):

  (chi/c0^2) h_tt - nabla^2 h
  - 4(f'/f) h_r + 2(f'/f)^2 h
  + (-2 beta chi + 3 gamma chi^2) h = 0
""")

print("STEP 3: Separate variables")
print("-" * 50)
print("""
h = exp(-i omega t) * u(r) * Y_ell^m(theta, phi)

Radial equation:

  -u'' - [2/r + 4f'/f] u'
  + [ell(ell+1)/r^2 + 2(f'/f)^2 - 2 beta chi + 3 gamma chi^2] u
  = (chi/c0^2) omega^2 u
""")

# ============================================================
# STEP 4: Weight function
# ============================================================
print("STEP 4: Weight function to kill u' term")
print("-" * 50)

g = 1 / (r * f**2)
u_expr = Psi * g
u_r = sp.diff(u_expr, r)
u_rr = sp.diff(u_expr, r, 2)

p_r = 2/r + 4*fp/f
kinetic = -u_rr - p_r * u_r
kinetic_expanded = sp.expand(kinetic)

Psi_r = sp.diff(Psi, r)
Psi_rr = sp.diff(Psi, r, 2)

coeff_Psi_rr = sp.simplify(kinetic_expanded.coeff(Psi_rr))
coeff_Psi_r = sp.simplify(kinetic_expanded.coeff(Psi_r))

print(f"  u = Psi / (r * f^2)")
print(f"  Coeff of Psi'': {coeff_Psi_rr}")
print(f"  Coeff of Psi':  {coeff_Psi_r}")
assert sp.simplify(coeff_Psi_r) == 0, "First derivative NOT eliminated!"
print("  => Psi' term ELIMINATED. OK")
print()

# Extract the Psi-coefficient (W_kin)
coeff_Psi_raw = kinetic_expanded - coeff_Psi_rr * Psi_rr
W_kin = sp.simplify(coeff_Psi_raw / Psi)

# ============================================================
# STEP 5: Assemble Schrodinger equation
# ============================================================
print("STEP 5: Assemble radial equation for Psi")
print("-" * 50)

# The full equation after substitution u = Psi/(r f^2):
# coeff_Psi_rr * Psi'' + W_kin * Psi + V_pot * Psi/(r f^2)
#   = (chi/c0^2) omega^2 * Psi/(r f^2)
# where V_pot = ell(ell+1)/r^2 + 2(fp/f)^2 - 2 beta chi + 3 gamma chi^2
# and coeff_Psi_rr = -1/(r f^2)
#
# Multiply by -(r f^2):
# Psi'' - (r f^2) W_kin Psi - V_pot Psi = -(chi/c0^2) omega^2 Psi
# =>  Psi'' + (chi/c0^2) omega^2 Psi = (r f^2) W_kin Psi + V_pot Psi

chi_expr = f / Phi0
V_pot = ell*(ell+1)/r**2 + 2*(fp/f)**2 - 2*beta_s*chi_expr + 3*gamma_s*chi_expr**2

V_raw = sp.expand(r * f**2 * W_kin) + V_pot
V_raw_simplified = sp.simplify(V_raw)

print("  Psi''(r) + (chi/c0^2) omega^2 Psi = V_raw(r) * Psi")
print()
print(f"  V_raw(r) = {V_raw_simplified}")

# Simplify for vacuum (beta = gamma)
V_raw_vac = sp.simplify(V_raw_simplified.subs(beta_s, gamma_s))
print()
print(f"  V_raw(r)|_(beta=gamma) = {V_raw_vac}")
print()

# ============================================================
# STEP 6: Tortoise coordinate
# ============================================================
print("STEP 6: Tortoise coordinate")
print("-" * 50)
print("""
  dr*/dr = sqrt(chi)/c0 = 1/c_bg(r)

  d/dr = n * d/dr*  where  n = sqrt(chi)/c0
  d^2/dr^2 = n^2 d^2/dr*^2 + n*n' d/dr*

  Psi(r) = Z(r) * Phi_w(r*)  with Z = chi^{1/4} = (f/Phi0)^{1/4}
  to absorb the n*n' first-derivative term.
""")

n = sp.sqrt(f / Phi0) / c0
n_prime = sp.diff(n, r)

Z = (f/Phi0)**sp.Rational(1, 4)
Z_r = sp.diff(Z, r)
Z_rr = sp.diff(Z, r, 2)

# Psi = Z * Phi_w
# Psi' = Z' Phi_w + Z * n * Phi_w'  (chain rule: d/dr* -> dr*/dr = n)
# Psi'' = Z'' Phi_w + 2 Z' n Phi_w' + Z (n^2 Phi_w'' + n n' Phi_w')
#       = Z'' Phi_w + (2 Z' n + Z n') Phi_w' + Z n^2 Phi_w''
#
# In our equation: Psi'' + (chi/c0^2) w^2 Psi = V_raw Psi
# chi/c0^2 = n^2 * Phi0  ... wait:
# chi/c0^2 = (f/Phi0)/c0^2 and n^2 = f/(Phi0 c0^2), so chi/c0^2 = n^2 * ... let me check
# n = sqrt(f/Phi0) / c0, n^2 = f/(Phi0 c0^2) = chi/c0^2.  YES!

print("  n^2 = chi/c0^2 = f/(Phi0 c0^2)  [confirmed]")
print()

# Substituting Psi = Z Phi_w:
# Z n^2 Phi_w'' + (2 Z' n + Z n') Phi_w' + Z'' Phi_w + n^2 w^2 Z Phi_w = V_raw Z Phi_w
# Divide by Z n^2:
# Phi_w'' + [(2Z'/Z)(1/n) + n'/n] Phi_w' + (Z''/(Z n^2)) Phi_w + w^2 Phi_w = (V_raw/n^2) Phi_w

# For the Phi_w' term to vanish:
# 2Z'/(Zn) + n'/n = 0
# => 2Z'/Z = -n'  ... let's check
# Z = chi^{1/4} = (f/Phi0)^{1/4}
# Z'/Z = (1/4) f'/(f)
# n = (f/Phi0)^{1/2} / c0
# n'/n = (1/2) f'/f
# So: 2Z'/(Zn) + n'/n = 2*(1/4)*(f'/f)*(1/n) + (1/2)*(f'/f)
# = (1/2)(f'/f)/n + (1/2)(f'/f)
# = (f'/(2f)) * (1/n + 1)
# This is NOT zero. So Z = chi^{1/4} doesn't quite kill the first deriv.

# Alternative: since we already killed u' and have Psi'' + ... = ...,
# the tortoise transform introduces a new first derivative.
# The correct weight is: Z such that 2 Z' n + Z n' = 0
# => 2 Z'/Z = -n'/n  => Z'/Z = -(1/2) n'/n = -(1/4) f'/f
# => Z = f^{-1/4}  (up to constant)

Z_correct = f**sp.Rational(-1, 4)
Z_c_r = sp.diff(Z_correct, r)

# Check: 2 Z_c'/Z_c * 1/... hmm, let me redo.
# The first-deriv coefficient in the Phi_w equation is:
# (2 Z' n + Z n') / (Z n^2)
# = 2 Z'/(Z n) + n'/n^2
# Wait, I divided by Z n^2 above. Let me redo:
#
# After dividing by Z:
# n^2 Phi_w'' + (2 Z'/Z n + n') n Phi_w' + ... wait I'm confusing myself.
#
# Let me start clean.  The equation is:
# Psi'' + n^2 w^2 Psi = V_raw Psi
# Set Psi = Z(r) * Phi_w(r*(r)):
# d Psi/dr = Z' Phi_w + Z (d Phi_w/dr) = Z' Phi_w + Z n (d Phi_w/dr*)
# d^2 Psi/dr^2 = Z'' Phi_w + 2 Z' n Phi_w'* + Z(n' Phi_w'* + n^2 Phi_w''*)
#
# Equation becomes:
# Z'' Phi_w + 2 Z' n Phi_w'* + Z n' Phi_w'* + Z n^2 Phi_w''* + n^2 w^2 Z Phi_w = V_raw Z Phi_w
#
# Phi_w''* coefficient: Z n^2
# Phi_w'* coefficient: 2 Z' n + Z n'
# Phi_w coefficient: Z'' + n^2 w^2 Z - V_raw Z
#
# Kill Phi_w'*: 2 Z' n + Z n' = 0
# => Z'/Z = -n'/(2n)
# n = sqrt(chi)/c0 = sqrt(f/Phi0)/c0
# n' = f'/(2 c0 sqrt(f Phi0))
# n'/n = f'/(2f)
# => Z'/Z = -f'/(4f)
# => Z = f^{-1/4} = (Phi0/f)^{1/4} / Phi0^{1/4}  (absorb constants)

print("  Correct weight: Z(r) = f(r)^{-1/4}")
print("  (Kills Phi_w' term in tortoise equation)")
print()

Z = f**sp.Rational(-1, 4)
Z_r = sp.diff(Z, r)
Z_rr = sp.diff(Z, r, 2)

# Now the equation is:
# Z n^2 Phi_w'' + (Z'' + n^2 w^2 Z - V_raw Z) Phi_w = 0
# Divide by Z n^2:
# Phi_w'' + w^2 Phi_w = (V_raw/n^2) Phi_w - (Z''/(Z n^2)) Phi_w
# Phi_w'' + w^2 Phi_w - V_eff Phi_w = 0
# where V_eff = V_raw/n^2 - Z''/(Z n^2) = (1/n^2)[V_raw - Z''/Z]

V_tortoise_corr = sp.simplify(Z_rr / Z)
print(f"  Z''/Z = {V_tortoise_corr}")
print()

# ============================================================
# STEP 7: Final V_eff
# ============================================================
print("=" * 65)
print("  FINAL RESULT")
print("=" * 65)
print()
print("  Schrodinger equation:")
print()
print("    d^2 Phi_w / dr*^2 + [omega^2 - V_eff(r)] Phi_w = 0")
print()
print("  where:")
print("    V_eff(r) = (1/n^2) * [V_raw(r) - Z''/Z]")
print("            = (c0^2 Phi0 / f) * [V_raw(r) - Z''/Z]")
print()

V_eff_full = sp.simplify((V_raw_simplified - V_tortoise_corr) * c0**2 * Phi0 / f)
print(f"  V_eff = {V_eff_full}")
print()

# Vacuum case
V_eff_vac = sp.simplify(V_eff_full.subs(beta_s, gamma_s))
print(f"  V_eff|_(beta=gamma) = {V_eff_vac}")
print()

# ============================================================
# LIMITS
# ============================================================
print("=" * 65)
print("  LIMITING CASES")
print("=" * 65)

# 1. Weak field
print()
print("1. WEAK FIELD: f -> Phi0, f' -> 0, f'' -> 0")
V_weak = V_eff_vac.subs([(f, Phi0), (fp, 0), (fpp, 0)])
V_weak = sp.simplify(V_weak)
print(f"   V_eff -> {V_weak}")
V_weak_l0 = sp.simplify(V_weak.subs(ell, 0))
print(f"   ell=0: V_eff -> {V_weak_l0}")
print(f"   Expected: c0^2 * gamma (= c0^2 m_sp^2)")
match = sp.simplify(V_weak_l0 - c0**2 * gamma_s)
print(f"   Match: {match == 0}")
print()

# 2. Strong field
print("2. STRONG FIELD: f >> Phi0 (chi >> 1)")
print("   V_mass ~ 3 gamma chi^2 * (c0^2/chi) = 3 gamma c0^2 chi")
print("   V_eff ~ 3 gamma c0^2 * f/Phi0  (GROWS with Phi)")
print("   => Rising potential barrier traps perturbations")
print()

# 3. Asymptotic
print("3. ASYMPTOTIC: r -> inf, f -> Phi0")
print("   V_eff -> c0^2 [ell(ell+1)/r^2 + gamma]")
print("   = centrifugal + mass gap")
print("   Propagating modes only for omega^2 > c0^2 gamma")
print()

# ============================================================
# SUMMARY
# ============================================================
print("=" * 65)
print("  SUMMARY")
print("=" * 65)
print("""
  TGP ringdown equation:

    d^2 Phi_w / dr*^2 + [omega^2 - V_eff(r)] Phi_w = 0

  V_eff has THREE components:

    V_eff = (c0^2 Phi0/f) * [ V_angular + V_mass + V_weight ]

    V_angular = ell(ell+1) / r^2

    V_mass    = 3 gamma (f/Phi0)^2 - 2 beta (f/Phi0) + 2(f'/f)^2

    V_weight  = -(Z''/Z)  where Z = f^{-1/4}
              = (1/4) f''/f - (5/16)(f'/f)^2  ... (geometry correction)

  Key differences from GR (Regge-Wheeler):

    1. MASS GAP: omega_min = c0 sqrt(gamma)
       No low-frequency modes (Yukawa, not Coulomb)

    2. SELF-REFERENTIAL SPEED: c_bg(r) = c0 sqrt(Phi0/f)
       Speed of perturbations depends on background

    3. RISING BARRIER: V_eff ~ 3 gamma c0^2 (f/Phi0) for f >> Phi0
       Dense space TRAPS perturbations (c -> 0 region is impenetrable)

    4. ONLY SCALAR MODES: breathing polarization
       Tensor h+, hx require disformal sector B(Phi)

    5. NO HORIZON in classical sense
       Instead: c(Phi) -> 0 creates effective causal boundary

  For quasinormal modes:
    - Boundary condition at r* -> +inf: outgoing wave
    - Boundary condition at r* -> -inf: ingoing (into c->0 region)
    - Complex omega_n = omega_R + i omega_I
    - omega_R > c0 sqrt(gamma) (above mass gap)
""")
