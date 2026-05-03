#!/usr/bin/env python3
"""
Z_precise_derivation.py — Opcja Z: precyzyjna derywacja s_natural z sek08a action

Cel: zamknac niepewnosc s_natural ~ 2.5e-3 z handwave dimensional analysis.
Zapuszczam pelna derywacje:
  1. sek08a S[Phi, rho] symbolically (sympy)
  2. EL equation -> Phi-EOM with matter source
  3. FRW reduction -> dimensionless s_natural in de2 normalization
  4. Plug in G.0 + UV.3 + gamma.1 + delta.1+2 closures
  5. Re-run FRW with precise s_natural
  6. Final verdict: czy TGP rozwiazuje Hubble tension?

Source folder cross-refs:
  research/op-uv3-phi0-renormalization/         (Z_Phi = 14/3)
  research/op-gamma1-phi-eff-anchor-resolution/ (Omega_L^pure = 2pi/9)
  research/op-delta1-g-tilde-derivation/        (g_tilde = N_f e^2/(12 pi))
  research/op-delta2-Nf-derivation/             (N_f = 5)
  closure_2026-04-26/Lambda_from_Phi0/          (rho_vac = M_Pl^2 H_0^2/12)

Convention tracking:
  q*Phi_0 = ?  (need to fix from Newton limit)
  K(phi) = K_geo * phi^4, K_geo = ?
  V(psi) = (beta/3) psi^3 - (gamma/4) psi^4 with beta=gamma
  rho_vac,TGP = gamma/12 = M_Pl^2 H_0^2 g_tilde / 12 (T-Lambda)
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


# ============================================================================
# STEP 1: Symbolic derivation of Phi-EOM in FRW with matter source
# ============================================================================
print("=" * 78)
print("STEP 1: Symbolic Phi-EOM derivation z sek08a action")
print("=" * 78)

# Define symbols
t, a = sp.symbols('t a', positive=True, real=True)
phi, phi_dot, phi_ddot = sp.symbols('phi phi_dot phi_ddot', real=True)
rho, H = sp.symbols('rho H', positive=True, real=True)
Phi_0, K_geo, beta, gamma_param, q, c0, G = sp.symbols(
    'Phi_0 K_geo beta gamma q c_0 G', positive=True, real=True)

# Field psi = Phi/Phi_0 (dimensionless)
# Action: S = integral d^4x sqrt(-g) [(1/2) K(phi) g^mu_nu d_mu phi d_nu phi - V(phi) - (q/Phi_0) phi rho]
# For FRW (homogeneous): only time derivative
# K(phi) = K_geo * phi^4
# V(phi) = (beta/3) phi^3 - (gamma/4) phi^4

K_func = K_geo * phi**4
K_prime = sp.diff(K_func, phi)
V_func = (beta/3) * phi**3 - (gamma_param/4) * phi**4
V_prime = sp.diff(V_func, phi)

# Vacuum condition: V'(phi=1) = 0 -> beta = gamma
# So substitute beta = gamma in V_prime
V_prime_vacuum = V_prime.subs(beta, gamma_param).simplify()
print(f"\n  V(phi)         = {sp.simplify(V_func.subs(beta, gamma_param))}")
print(f"  V'(phi)        = {sp.simplify(V_prime_vacuum)}")
print(f"  V'(phi=1)      = {V_prime_vacuum.subs(phi, 1)}")
print(f"  V''(phi=1)     = {sp.diff(V_prime_vacuum, phi).subs(phi, 1)}")

# EOM: K * phi_ddot + 3 H K * phi_dot + (1/2) K' * phi_dot^2 + V'(phi) + (q/Phi_0) * rho = 0
EOM = (K_func * phi_ddot
       + 3 * H * K_func * phi_dot
       + sp.Rational(1, 2) * K_prime * phi_dot**2
       + V_prime_vacuum
       + (q / Phi_0) * rho)

print(f"\n  Full EOM (Phi-coupled):")
print(f"    {sp.simplify(EOM)}")

# Solve for phi_ddot
phi_ddot_solved = sp.solve(EOM, phi_ddot)[0]
print(f"\n  phi_ddot = {sp.simplify(phi_ddot_solved)}")


# ============================================================================
# STEP 2: Convert to de2 dimensionless form
# ============================================================================
print("\n" + "=" * 78)
print("STEP 2: De2 dimensionless reduction")
print("=" * 78)

# In de2: V(psi) = U(psi)/U(1) = (4 psi^3 - 3 psi^4) (V_norm)
# Original V(psi) = (gamma/12) * (4 psi^3 - 3 psi^4)
# So V_orig = (gamma/12) * V_norm
# In de2: rho_psi = (1/2) psi_dot^2 + V_0 * V_norm(psi), where V_0 ~ Omega_DE
# That means V_0 absorbs (gamma/12) factor: V_0 = gamma/12 * (Phi_0^2 c^2 ... etc)

# Let me be careful. de2 uses:
#   rho_psi_dim = (1/2) psi_dot^2 + V_0 V_norm(psi)
# where psi_dot is in units (H_0 t = tau)
# and rho_psi_dim is in units rho_crit_today
# This implies absorption of factors to make it dimensionless

# Original action: rho_psi_phys = (1/2) K(psi) (Phi_0 dot{psi}/c)^2 + V_orig(psi)
# wait, careful: phi_dot in action has dimensions [phi]/[time]; phi is dimensionless,
# so phi_dot has [time]^-1
# and gradient terms in action use g^mu_nu so timelike component gets c^-2

# Actually let me use a simplification. Just work in units c=1, H_0=1 (de2 convention).
# Then everything is dimensionless and we can match coefficients.

print("""
  Konwencja de2:
    - Units: H_0 = 1, c = 1, 8 pi G/3 = 1, rho_crit_0 = 1
    - psi = Phi / Phi_0
    - tau = H_0 * t
    - Friedmann: H^2 = Omega_r/a^4 + Omega_m/a^3 + rho_psi
    - rho_psi = (1/2) psi_dot^2 + V_0 V_norm(psi)
    - V_0 = Omega_DE0 (shooting parameter; absorbs gamma/12 normalization)
    - V_norm(psi) = 4 psi^3 - 3 psi^4 (V_norm(1)=1, V_norm'(1)=0)
""")

# In de2 units, sek08a EOM transforms to:
#   K(psi) * [psi_ddot + 3 H psi_dot] + (1/2) K'(psi) psi_dot^2 + V'(psi)/(c^2 ... factor) + source = 0
# After absorbing V_0 normalization:
#   psi_ddot + 3 H psi_dot + (2/psi) psi_dot^2 (from K'/K = 4/psi)
#   + V_0 * V_norm'(psi) / [psi^4 K_geo * normalization]
#   + s * Omega_m / [psi^4 K_geo * normalization]
#   = 0

# At vacuum psi ~ 1, K = K_geo * 1 = K_geo. So leading-order in slow-roll:
#   psi_ddot + 3 H psi_dot + V_0 V_norm'(psi) / K_geo + s_eff Omega_m / K_geo = 0

# The matter source coefficient in dimensionless EOM:
#   s_eff = (q/Phi_0) rho_phys / [V_0 * factor for V']
# For de2 normalization V_0 V_norm'(psi) ~ V_orig'(psi)/(rho_crit), so:
#   s_eff = (q/Phi_0) rho_crit Omega_m / rho_crit = (q/Phi_0) for Omega_m factor in rho_m_dim

# Wait I need to be more careful with dimensional analysis.

# In PHYSICAL units, EOM (with K_geo absorbed):
#   psi_ddot + 3 H psi_dot + V_orig'(psi)/Phi_0^2 + (q/Phi_0) rho_phys / Phi_0^2 + ... = 0
# Wait no. Let me redo.

# Physical action: S = integral d^4x sqrt(-g) [(1/2) K(phi)(d phi)^2 - V(phi) - (q/Phi_0) phi rho]
# Where phi = Phi/Phi_0 dimensionless, V has units of energy density.
# V(phi) = (gamma/12) * (4 phi^3 - 3 phi^4) where gamma ~ rho_vac (energy density).

# EOM divided through: phi_ddot + ... + V'(phi)/K + (q/Phi_0) rho/K = 0
# Multiply by 1/H^2 to make dimensionless:
#   (phi_ddot/H^2) + (3 H phi_dot)/H^2 + V'(phi)/(K H^2) + (q/Phi_0) rho/(K H^2) = 0

# In de2 dimensionless time tau = H_0 t:
#   psi_ddot_tau = psi_ddot / H_0^2 (where psi_ddot now in tau-units)
# So everything scaled by H_0^2.

# At vacuum K(1) = K_geo: V'(psi)/(K_geo H_0^2) = V'(psi) / (K_geo H_0^2)
# V_norm'(psi) (de2) absorbs gamma/12 factor:
#   V_orig(psi) = (gamma/12) V_norm(psi)
#   V_orig'(psi) = (gamma/12) V_norm'(psi)
# de2 EOM: psi_ddot + 3H psi_dot + V_0 V_norm'(psi) = 0
# Comparing: V_0 = (gamma/12) / (K_geo H_0^2)

# T-Lambda closure: gamma/12 = M_Pl^2 H_0^2 g_tilde / 12 (with g_tilde ~ 1)
# so gamma = M_Pl^2 H_0^2 g_tilde

# Physical V_0:
#   V_0 = (gamma/12) / (K_geo H_0^2) = M_Pl^2 g_tilde / (12 K_geo)

# But de2 sets V_0 = Omega_DE0 by shooting. So:
#   Omega_DE0 = M_Pl^2 g_tilde / (12 K_geo)
# This implies: K_geo = M_Pl^2 g_tilde / (12 Omega_DE0)
# In Planck units (M_Pl = 1): K_geo = g_tilde / (12 Omega_DE0) ~ 1/8.2 ~ 0.12

# Source term: (q/Phi_0) rho / (K_geo H_0^2)
#   q*Phi_0 = 4 pi G / c^2 (Newton limit)
#   q/Phi_0 = 4 pi G / (c^2 Phi_0^2)
# Source/H_0^2 = (4 pi G / c^2) rho / (Phi_0^2 K_geo H_0^2)
# Using rho/rho_crit = Omega_m and rho_crit = 3 H_0^2 c^2/(8 pi G c^4):
# Hmm let me carefully use units.

# In dimensionful units:
#   rho_crit = 3 H_0^2 / (8 pi G) [c=1 units, mass density]
#   rho_m = Omega_m rho_crit = Omega_m * 3 H_0^2 / (8 pi G)
#   (q/Phi_0) rho_m = (4 pi G / Phi_0^2) * Omega_m * 3 H_0^2 / (8 pi G) = 3 Omega_m H_0^2 / (2 Phi_0^2)
# So source/H_0^2 = 3 Omega_m / (2 Phi_0^2)
# Compare to V_0 V_norm'(psi):
#   V_0 = Omega_DE0 ~ 0.685
#   V_norm'(psi) ~ O(0.3) at slow-roll psi ~ 0.97
# V_0 * V_norm' ~ 0.2

# Source/restoring = 3 Omega_m / (2 Phi_0^2) / [V_0 V_norm'(psi)]
#                  = (3/2)(0.315)/(24.65^2)/(0.2) ~ 0.012

# So source is ~1% of restoring force. But wait, my numerical run gave dpsi ~ 3% with
# this s. Discrepancy?

# Let me re-examine. The "source coupling s" in de2-normalized EOM is:
#   psi_ddot + 3 H psi_dot + V_0 V_norm'(psi) + s_de2 Omega_m_dim(a) = 0
# where Omega_m_dim(a) = Omega_m_0/a^3 in dimensionless units.

# From my derivation:
#   s_de2 (in K_geo) = source/H_0^2 / Omega_m(a) = (3/2)/Phi_0^2 / K_geo

# Wait! I forgot K_geo! Let me redo:
# Source term in physical EOM: (q/Phi_0) rho [Pa K_geo factor]
# Dividing through by K(psi=1) = K_geo:
#   (q/Phi_0) rho / K_geo
# In dimensionless de2:
#   (q/Phi_0) rho / (K_geo H_0^2) = 3 Omega_m / (2 Phi_0^2 K_geo)
# So s_de2 = 3 / (2 Phi_0^2 K_geo)

# With K_geo ~ 1/8.2 ~ 0.12 (from V_0 = Omega_DE matching above):
# s_de2 = 3 / (2 * 24.65^2 * 0.12) ~ 0.020

# That's 8x BIGGER than my earlier estimate (2.5e-3)!
# K_geo factor was missed in initial dimensional analysis.

print("""
  KORRECTA dimensional analysis (z K_geo non-canonical kinetic):
    Source term in physical EOM: (q/Phi_0) rho_m [from L_mat]
    Divided through by K(psi=1) = K_geo:
       Effective source = (q/Phi_0) rho_m / K_geo

    From Newton limit: q*Phi_0 = 4 pi G / c^2
    So: q/Phi_0 = 4 pi G / (c^2 Phi_0^2)

    Source in dimensionless de2 form:
       s_de2 * Omega_m = (4 pi G rho_m) / (c^2 Phi_0^2 K_geo H_0^2)
       Using rho_m = Omega_m * rho_crit, rho_crit = 3 H_0^2 c^2 / (8 pi G c^4):
       s_de2 = 3 / (2 Phi_0^2 K_geo)
""")

# K_geo derivation from T-Lambda + de2 V_0:
# V_0 (de2 shooting) = gamma/(12 K_geo H_0^2) = Omega_DE0
# T-Lambda: gamma = M_Pl^2 H_0^2 g_tilde (with g_tilde ~ 1)
# Therefore: V_0 = M_Pl^2 g_tilde / (12 K_geo)
# In dimensionless units M_Pl = 1: K_geo = g_tilde / (12 Omega_DE0)
print(f"  K_geo derivation:")
print(f"    From T-Lambda: gamma = M_Pl^2 H_0^2 g_tilde")
print(f"    V_0 (de2) = gamma/(12 K_geo H_0^2) = Omega_DE0")
print(f"    => K_geo = g_tilde * M_Pl^2 / (12 Omega_DE0)")

# In Planck units (M_Pl = 1):
g_tilde = 5.0 * np.exp(2)/(12 * np.pi)  # delta.1: g_tilde = 5 e^2 / (12 pi)
# Actually e^2 in delta.1 is Euler's e squared
# Wait, let me re-read. delta.1: g_tilde = N_f * e^2 / (12 pi) with N_f = 5
# Phi_eff^corr = (10/3) e^2 = 8 pi g_tilde, so g_tilde = (10/3) e^2 / (8 pi) = 5 e^2 / (12 pi)
# e here is Euler's number
e_euler = np.exp(1)
g_tilde_val = 5 * e_euler**2 / (12 * np.pi)
print(f"    g_tilde = 5 e^2/(12 pi) = {g_tilde_val:.6f}  (delta.1+delta.2 closure)")

OMEGA_DE0 = 0.6847
M_PL = 1.0  # Planck units

K_geo_val = g_tilde_val * M_PL**2 / (12 * OMEGA_DE0)
print(f"    K_geo (Planck units, M_Pl=1) = {K_geo_val:.6f}")

# But wait — Phi_0 is in same Planck units?
# Phi_0 has dimensions [length]^-2 = [mass]^2 in c=hbar=1
# In Planck units, Phi_0 ~ M_Pl^2 ? Or Phi_0 = 24.65 dimensionless?

# From sek00: Phi_eff = 24.65 (dimensionless? or in some unit?)
# Looking at sek00 formula: Omega_L = Phi_0/36 (linia 401)
# That implies Phi_0 has units such that Omega_L is dimensionless.
# Omega_L = rho_vac / rho_crit = (gamma/12)/(3 H_0^2/(8 pi G c^4))
#        = (gamma * 8 pi G c^4) / (36 H_0^2)
# Setting this = Phi_0/36:
#   Phi_0 = (gamma * 8 pi G c^4) / H_0^2

# In Planck units c=hbar=G=1, H_0 = ~1e-60 (in M_Pl units), so:
#   Phi_0 = (gamma * 8 pi) / H_0^2

# With gamma = M_Pl^2 H_0^2 g_tilde = H_0^2 g_tilde (M_Pl = 1):
#   Phi_0 = 8 pi g_tilde

# Numerically: 8 pi * g_tilde ≈ 8 pi * 0.308 ≈ 7.7
# But sek00 says Phi_eff ≈ 24.65. Discrepancy by factor 3.2.

# Looking at sek00 more carefully: "Phi_eff^pure = 8 pi" and "Phi_eff^corr = (10/3) e^2 ≈ 24.6302"
# So in TGP normalization, Phi_eff = 8 pi g_tilde * (some_factor).
# Actually: 10/3 e^2 = 8 pi * (e^2)/(8 pi 3/10) = ... let me just verify numerically:
# 8 pi g_tilde with g_tilde = 5 e^2/(12 pi) = 8 pi * 5 e^2/(12 pi) = 40 e^2 / 12 = (10/3) e^2 ✓

# So Phi_eff = 8 pi * g_tilde = (10/3) e^2 ≈ 24.6302. Confirms!

# Therefore in Planck units (M_Pl=1, c=hbar=G=1, H_0 ~ 1e-60):
# Phi_0 (= Phi_eff) ≈ 24.65 (dimensionless natural number!)

print(f"\n  Cross-check: Phi_eff = 8 pi * g_tilde = {8 * np.pi * g_tilde_val:.4f}")
print(f"  Sek00 value: Phi_eff = 24.6302 (10/3 * e^2)")

# Note: Phi_0 is dimensionless when expressed as Phi_0 = Phi_eff (the IR scale).
# It's NOT in M_Pl^2 units.

# So in s_de2 = 3 / (2 Phi_0^2 K_geo):
PHI_0_val = 24.6302
s_de2_naive = 3.0 / (2 * PHI_0_val**2 * K_geo_val)
print(f"\n  s_de2 = 3 / (2 Phi_0^2 K_geo) = 3 / (2 * {PHI_0_val:.2f}^2 * {K_geo_val:.4f}) = {s_de2_naive:.4e}")

# But wait - this says K_geo is ~0.04 in our units, making s ~0.02.
# Let me check this differently.


# ============================================================================
# STEP 3: Cross-check using sek08a Lagrangian directly
# ============================================================================
print("\n" + "=" * 78)
print("STEP 3: Direct cross-check from sek08a Lagrangian densities")
print("=" * 78)

# Energy density of Phi field: rho_phi = (1/2) K(phi) phi_dot^2 + V_orig(phi)
# At vacuum (phi=1, phi_dot=0): rho_phi_vac = V_orig(1) = gamma/12
# Observation: rho_phi_vac = Omega_DE0 * rho_crit

# In Planck units (M_Pl=1, H_0 small):
# rho_crit = 3 H_0^2 / (8 pi G) = 3 H_0^2 / (8 pi)
# So rho_phi_vac = Omega_DE0 * 3 H_0^2 / (8 pi) = 0.685 * 3 H_0^2 / (8 pi)
# = gamma / 12
# => gamma = 12 * 0.685 * 3 H_0^2 / (8 pi) = (12*0.685*3)/(8 pi) * H_0^2 = 0.981 H_0^2
# Coincidentally g_tilde ≈ 0.98 ✓

# So gamma = 0.981 H_0^2 in Planck units (M_Pl=1).

# Now K(phi) factor: at vacuum K(1) = K_geo. The kinetic energy at vacuum is 0 (phi_dot=0),
# but the K_geo factor enters the EOM through (1/2) K phi_dot^2 + V.
# When we go to dimensionless de2 form: rho_phi_dim = rho_phi / rho_crit
# (1/2) K(1) (Phi_0 dot{psi})^2 / rho_crit + V(psi) / rho_crit
# The first term: (1/2) K_geo Phi_0^2 psi_dot^2 / rho_crit = (1/2) psi_dot^2 if K_geo Phi_0^2 = rho_crit
# i.e., K_geo = rho_crit / Phi_0^2 = 3 H_0^2 / (8 pi Phi_0^2)
# In Planck units (M_Pl=1, H_0=H_0): K_geo = 3 H_0^2 / (8 pi Phi_0^2)
# Numerically (in units H_0^2): K_geo = 3 / (8 pi * 24.63^2) = 1.96e-4
# Hmm tiny.

# But wait: de2 uses (1/2) psi_dot^2 with psi_dot in units 1/(H_0 t) so psi_dot is order 1.
# Physical units phi_dot has [1/time]. dimensionless psi_dot_de2 = phi_dot / H_0.
# So (1/2) psi_dot_de2^2 (de2 dimensionless) = (1/2) phi_dot^2 / H_0^2
# Multiplied by H_0^2 to get physical kinetic energy: (1/2) phi_dot^2 / H_0^2 * H_0^2 = (1/2) phi_dot^2 (no units?)

# OK this is getting tangled. Let me just write the Lagrangian explicitly with all units and match.

print("""
  Direct method: match physical and dimensionless EOMs term by term.

  Physical EOM at vacuum (K = K_geo, V_0 = gamma/12, V_norm = de2):
    K_geo phi_ddot + 3H K_geo phi_dot + V_0' V_norm'(phi) + (q/Phi_0) rho_m = 0

  De2 form: psi_ddot + 3 H psi_dot + V_0_de2 V_norm'(psi) + s_de2 Omega_m = 0

  Identifying: psi = phi (dimensionless), tau = H_0 t (so psi_ddot in tau is phi_ddot / H_0^2)

  For physical gamma/12 V_norm' to match V_0_de2 V_norm':
    V_0_de2 = gamma/(12 K_geo H_0^2)

  Setting V_0_de2 = Omega_DE0 (de2 shooting):
    Omega_DE0 = gamma/(12 K_geo H_0^2)

  T-Lambda closure: gamma = g_tilde M_Pl^2 H_0^2:
    Omega_DE0 = g_tilde M_Pl^2 / (12 K_geo)

  Solving: K_geo = g_tilde M_Pl^2 / (12 Omega_DE0)
""")

# In Planck units M_Pl = 1, so K_geo = g_tilde / (12 Omega_DE0):
print(f"  K_geo (Planck units M_Pl=1) = g_tilde / (12 Omega_DE0) = {g_tilde_val / (12 * OMEGA_DE0):.6f}")

# Source term in dimensionless EOM:
#   (q/Phi_0) rho_m / (K_geo H_0^2) = s_de2 * Omega_m_a(a)
# At today (rho_m = Omega_m rho_crit):
#   s_de2 = (q/Phi_0)(rho_crit) / (K_geo H_0^2 Omega_m) * Omega_m = (q/Phi_0)(rho_crit)/(K_geo H_0^2)

# q*Phi_0 = 4 pi G / c^2 (Newton limit, c=1 units): q = 4 pi G / Phi_0 (in units c=1)
# q/Phi_0 = 4 pi G / Phi_0^2

# rho_crit = 3 H_0^2 / (8 pi G) (units c=1)
# So (q/Phi_0)(rho_crit) = (4 pi G / Phi_0^2) * 3 H_0^2 / (8 pi G) = 3 H_0^2 / (2 Phi_0^2)

# Therefore: s_de2 = 3 H_0^2 / (2 Phi_0^2 K_geo H_0^2) = 3 / (2 Phi_0^2 K_geo)

# With K_geo = g_tilde / (12 Omega_DE0):
# s_de2 = 3 * 12 * Omega_DE0 / (2 Phi_0^2 g_tilde) = 18 Omega_DE0 / (Phi_0^2 g_tilde)

s_de2_precise = 18 * OMEGA_DE0 / (PHI_0_val**2 * g_tilde_val)
print(f"\n  s_de2 (precise) = 18 Omega_DE0 / (Phi_0^2 g_tilde)")
print(f"                  = 18 * {OMEGA_DE0} / ({PHI_0_val:.2f}^2 * {g_tilde_val:.4f})")
print(f"                  = {s_de2_precise:.4e}")

# Compare to my earlier estimate s_natural = 3/(2 Phi_0^2) = 2.47e-3:
s_naive = 3.0 / (2 * PHI_0_val**2)
print(f"\n  Compare to naive (no K_geo): s_naive = 3/(2 Phi_0^2) = {s_naive:.4e}")
print(f"  Amplification factor: K_geo correction = s_precise/s_naive = {s_de2_precise/s_naive:.2f}")


# ============================================================================
# STEP 4: Algebraic simplification using closure cascade
# ============================================================================
print("\n" + "=" * 78)
print("STEP 4: Closed-form expression using closure cascade")
print("=" * 78)

# We have:
#   s_de2 = 18 Omega_DE0 / (Phi_0^2 g_tilde)
#   Phi_0 = 8 pi g_tilde (gamma.1 corrected; Phi_eff = 8 pi g_tilde)
# Therefore:
#   s_de2 = 18 Omega_DE0 / ((8 pi g_tilde)^2 g_tilde) = 18 Omega_DE0 / (64 pi^2 g_tilde^3)

s_de2_closed = 18 * OMEGA_DE0 / (64 * np.pi**2 * g_tilde_val**3)
print(f"\n  Closed form: s_de2 = 18 Omega_DE0 / (64 pi^2 g_tilde^3)")
print(f"                    = 18 * {OMEGA_DE0} / (64 * pi^2 * {g_tilde_val:.4f}^3)")
print(f"                    = {s_de2_closed:.4e}")

# Cross-check: using Omega_DE = Phi_0/36 (sek05 K3 formula):
# Phi_0 = 36 * Omega_DE
# s_de2 = 18 * Phi_0/36 / (Phi_0^2 g_tilde) = (1/2) / (Phi_0 g_tilde)
s_de2_omega_form = 0.5 / (PHI_0_val * g_tilde_val)
print(f"\n  Alt closed form: s_de2 = 1/(2 Phi_0 g_tilde) = {s_de2_omega_form:.4e}")
print(f"  (Cross-check should match above)")

# Final value
S_PRECISE = s_de2_precise
print(f"\n  >>> S_PRECISE (using all closures) = {S_PRECISE:.6e}")
print(f"  >>> S_NAIVE (initial estimate)     = {s_naive:.6e}")
print(f"  >>> Ratio = {S_PRECISE/s_naive:.2f}x")


# ============================================================================
# STEP 5: Re-run FRW solver with precise s_natural
# ============================================================================
print("\n" + "=" * 78)
print("STEP 5: Re-run FRW solver with precise s coupling")
print("=" * 78)

OMEGA_M0  = 0.315
OMEGA_R0  = 9.1e-5

V_0_global = OMEGA_DE0

def V(psi):
    return 4.0*psi**3 - 3.0*psi**4

def dVdpsi(psi):
    return 12.0 * psi**2 * (1.0 - psi)

def Omega_m_a(a):
    return OMEGA_M0 / a**3

def Omega_r_a(a):
    return OMEGA_R0 / a**4

def eom_N(N, y, s_c):
    psi, u = y
    a = np.exp(N)
    rho_total = OMEGA_R0/a**4 + OMEGA_M0/a**3 + 0.5*u**2 + V_0_global*V(psi)
    H = np.sqrt(max(rho_total, 1e-30))
    if H < 1e-30:
        return [0, 0]
    dpsi_dN = u / H
    du_dN = (-3.0*H*u - V_0_global*dVdpsi(psi) - s_c*Omega_m_a(a)) / H
    return [dpsi_dN, du_dN]

a_init = 1e-5
N_init = np.log(a_init)
N_eval = np.linspace(N_init, 0.0, 5000)

DELTA_H_REQUIRED = 0.0837

scenarios = [
    ('s = 0 (no source, baseline)',           0.0),
    ('s = naive (no K_geo, ~2.5e-3)',        s_naive),
    ('s = PRECISE (with K_geo)',              S_PRECISE),
    ('s = 2 x precise (uncertainty band)',    2 * S_PRECISE),
    ('s = 0.5 x precise (uncertainty band)',  0.5 * S_PRECISE),
]

print(f"\n  Required Hubble shift: |dH/H| = {DELTA_H_REQUIRED:.4f} ({DELTA_H_REQUIRED*100:.2f}%)\n")
print(f"  {'Scenario':<40} {'psi_today':>10} {'psi_recomb':>11} {'dpsi':>11} {'dL/L':>11} {'|dH/H|':>10}")
print(f"  {'-'*40} {'-'*10} {'-'*11} {'-'*11} {'-'*11} {'-'*10}")

results = []
for label, s_value in scenarios:
    sol = solve_ivp(eom_N, [N_init, 0.0], [1.0 - 1e-3, 0.0],
                    args=(s_value,), t_eval=N_eval, method='RK45',
                    rtol=1e-9, atol=1e-12)
    if not sol.success:
        print(f"  {label:<40} FAILED")
        continue

    psi = sol.y[0]
    u = sol.y[1]
    z = 1.0/np.exp(N_eval) - 1.0

    psi_today = psi[-1]
    idx_recomb = np.argmin(np.abs(z - 1100.0))
    psi_recomb = psi[idx_recomb]
    Vt = V(psi_today)
    Vr = V(psi_recomb)
    dLL = (Vt - Vr) / Vt if Vt != 0 else 0
    dH = 0.5 * abs(dLL) * OMEGA_DE0

    print(f"  {label:<40} {psi_today:>10.6f} {psi_recomb:>11.6f} {psi_today-psi_recomb:>+11.3e} {dLL:>+11.3e} {dH:>10.3e}")
    results.append({'label': label, 's': s_value, 'dLL': dLL, 'dH': dH, 'psi_today': psi_today, 'psi_recomb': psi_recomb})


# ============================================================================
# STEP 6: Final verdict
# ============================================================================
print("\n" + "=" * 78)
print("STEP 6: FINAL VERDICT")
print("=" * 78)

# Find precise scenario
precise_result = next((r for r in results if abs(r['s'] - S_PRECISE) < 1e-10), None)
if precise_result:
    ratio = precise_result['dH'] / DELTA_H_REQUIRED
    print(f"""
  Precise s_natural derivation (z UV.3 + gamma.1 + delta.1+2 + K_geo non-canonical):
    s_de2 = 18 Omega_DE0 / (Phi_eff^2 g_tilde)
          = 18 * {OMEGA_DE0} / ({PHI_0_val:.4f}^2 * {g_tilde_val:.4f})
          = {S_PRECISE:.4e}

  Hubble tension shift:
    Required:   |dH/H| = {DELTA_H_REQUIRED:.4f}
    TGP gives:  |dH/H| = {precise_result['dH']:.4e}
    Ratio:      TGP / required = {ratio:.4f} ({ratio*100:.1f}%)
""")

    if ratio > 0.5:
        print("  >>> TGP rozwiazuje H_0 tension w wiekszosci - **DUZA DROGE PRZESZLISMY**")
    elif ratio > 0.1:
        print("  >>> TGP daje istotny czesciowy wklad - z innymi efektami moglo by zamknac tension")
    elif ratio > 0.01:
        print("  >>> TGP daje partial wklad rzadu kilku procent - undershoot")
    else:
        print("  >>> TGP wklad nadal niewystarczajacy")

print("\n" + "=" * 78)
print("Done.")
