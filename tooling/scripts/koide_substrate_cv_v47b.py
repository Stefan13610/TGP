#!/usr/bin/env python3
"""
koide_substrate_cv_v47b.py -- PROPOSAL #10: SUBSTRATE -> CV=1 -> KOIDE

QUESTION:
  Does the discrete substrate (3D Ising at criticality) produce
  field moments that enforce CV(sqrt(m)) = 1 for generation masses?

APPROACH:
  The substrate is a 3D lattice with Z2 symmetry.
  Coarse-graining produces a continuous field Phi(x).
  Solitons in this field have masses m_k ~ A_tail(g0_k)^4.

  The KEY LINK missing: how does the substrate determine g0_k?

  CHAIN:
    Substrate (Ising 3D) -> coarse-grained field Phi(x) -> soliton ODE
    -> g0_k (central values) -> A_tail(g0_k) -> m_k = A_tail^4

  For CV(sqrt(m)) = 1, we need CV(A_tail^2) = 1, i.e.,
  std(A_k^2) = mean(A_k^2) for k=1,2,3.

  POSSIBLE MECHANISMS:
  1. CRITICAL FLUCTUATIONS: At the Wilson-Fisher fixed point,
     the field phi has specific moment ratios.
     The Binder cumulant U = 1 - <phi^4>/(3<phi^2>^2) has a
     universal value U* at criticality.
     Could this relate to CV = 1?

  2. COARSE-GRAINING LEVELS: Different block sizes L produce
     different Phi values. Three "generations" could correspond
     to three coarse-graining scales L_1 < L_2 < L_3.
     If Phi(L_k) = g0_k in some mapping, then the scale hierarchy
     determines mass hierarchy.

  3. SPIN CLUSTERS: Near T_c, the spin system has clusters of
     all sizes (scale-free). Three dominant cluster scales could
     map to three generations.

  4. FIELD DISTRIBUTION AT CRITICALITY: The probability P(phi)
     at the WF point has specific moments. The kurtosis, skewness,
     and higher moments are universal. Could the moment ratios
     map to CV = 1?

  We test each mechanism quantitatively.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI_GOLD = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0

# ==== 3D Ising critical exponents (universal) ====
NU = 0.6300          # correlation length exponent
ETA = 0.0363         # anomalous dimension
BETA_C = 0.3265      # order parameter exponent
GAMMA_C = 1.2372     # susceptibility exponent
DELTA_C = 4.789      # equation of state
ALPHA_C = 0.110      # specific heat (2 - 3*nu)

# Wilson-Fisher fixed point (normalized)
R_STAR = -2.251
U_STAR = 3.917

# Universal amplitude ratios (3D Ising)
# Binder cumulant at T_c:
U4_STAR = 0.6233     # U4* = 1 - <m^4>/(3<m^2>^2) at T_c (Ising, sc)
# This means: <m^4>/<m^2>^2 = 3(1-U4*) = 3*0.3767 = 1.1301


def solver_A(g0, r_max=400):
    def rhs(r, y):
        g, gp = y
        g = max(g, 1e-10)
        source = g**2 * (1.0 - g)
        cross = (2.0 / g) * gp**2
        if r < 1e-10:
            return [gp, (source - cross) / 3.0]
        return [gp, source - cross - 2.0 * gp / r]
    sol = solve_ivp(rhs, (1e-6, r_max), [g0, 0.0],
                    rtol=1e-12, atol=1e-14, max_step=0.02)
    return sol.t, sol.y[0]


def A_tail(g0, r_min=50, r_max=300):
    if g0 >= G0_CRIT - 0.001 or g0 <= 0.01:
        return 0.0
    try:
        r, g = solver_A(g0, r_max=r_max+50)
        mask = (r > r_min) & (r < r_max)
        if np.sum(mask) < 100:
            return 0.0
        rf = r[mask]
        h = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, h, rcond=None)[0]
        return np.sqrt(bc[0]**2 + bc[1]**2)
    except:
        return 0.0


def QK_masses(A_vals):
    A = np.array(A_vals)
    A2, A4 = A**2, A**4
    if np.sum(A4) < 1e-30:
        return 0.0
    return np.sum(A2)**2 / np.sum(A4)


def CV(vals):
    v = np.array(vals)
    return np.std(v) / np.mean(v) if np.mean(v) > 0 else 0


# ================================================================
print("=" * 70)
print("PROPOSAL #10: SUBSTRATE DISCRETE MOMENTS -> CV=1 -> KOIDE")
print("=" * 70)

# ================================================================
print("\n" + "=" * 70)
print("SECTION 1: UNIVERSAL MOMENT RATIOS AT CRITICALITY")
print("=" * 70)

# At the WF fixed point, the field distribution P(phi) has universal
# moment ratios. Key quantity:
# <phi^4> / <phi^2>^2 = 3(1 - U4*) where U4* is the Binder cumulant.

ratio_42 = 3 * (1 - U4_STAR)  # = <phi^4>/<phi^2>^2
print(f"\n  3D Ising at criticality:")
print(f"  Binder cumulant U4* = {U4_STAR:.4f}")
print(f"  <phi^4>/<phi^2>^2 = {ratio_42:.4f}")
print(f"  For Gaussian: this ratio = 3.0000")
print(f"  Deviation from Gaussian: {(ratio_42/3 - 1)*100:+.2f}%")

# The field distribution at WF is NON-Gaussian with kurtosis:
kurtosis_excess = ratio_42 / 3.0 - 1.0  # excess kurtosis (normalized)
print(f"  Excess kurtosis: {kurtosis_excess:.4f}")

# For Koide: CV(A^2) = 1 for three "effective masses."
# A^2 is the sqrt(mass). CV = 1 means std = mean.
# Is there a connection between kurtosis and CV?

# For a distribution with moments mu_2 = <x^2>, mu_4 = <x^4>:
# If we draw 3 samples from this distribution, the expected CV is:
# E[CV^2] = (sigma^2/mu^2) where sigma^2 = mu_2 - mu_1^2, mu = mu_1
# But at T_c, <phi> = 0 (symmetric phase), so mu_1 = 0.
# We should consider |phi| or phi^2 as the relevant variable.

print(f"\n  Connection to CV = 1:")
print(f"  At T_c: <phi> = 0 (Z2 symmetry). Relevant variable: phi^2.")
print(f"  For X = phi^2: <X> = <phi^2>, <X^2> = <phi^4>")
print(f"  CV(X) = sqrt(<X^2>/<X>^2 - 1) = sqrt({ratio_42:.4f} - 1) = {np.sqrt(ratio_42 - 1):.4f}")
print(f"  TARGET: CV = 1.0000")
print(f"  Deviation: {(np.sqrt(ratio_42 - 1) - 1.0)*100:+.2f}%")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 2: WHAT CV=1 REQUIRES FROM MOMENTS")
print("=" * 70)

# CV(X) = 1 for X = phi^2 means:
# sqrt(<phi^4>/<phi^2>^2 - 1) = 1
# <phi^4>/<phi^2>^2 = 2
# Binder cumulant: U4 = 1 - <phi^4>/(3<phi^2>^2) = 1 - 2/3 = 1/3

U4_needed = 1.0 / 3.0
print(f"\n  CV(phi^2) = 1 requires:")
print(f"  <phi^4>/<phi^2>^2 = 2.0000")
print(f"  Binder cumulant U4 = {U4_needed:.4f}")
print(f"  Actual Ising U4* = {U4_STAR:.4f}")
print(f"  Difference: {(U4_STAR - U4_needed)/U4_needed*100:+.2f}%")

# RESULT: U4* = 0.6233, needed = 0.3333
# These are VERY different (87% off).
# So the DIRECT moment ratio of phi^2 at criticality does NOT give CV=1.

# But maybe it's not phi^2 directly. Maybe it's through the ODE map.
# The ODE map g0 -> A_tail is highly nonlinear.
# Could the nonlinearity compensate?


# ================================================================
print("\n" + "=" * 70)
print("SECTION 3: COARSE-GRAINING SCALES AS GENERATIONS")
print("=" * 70)

# Idea: three generations correspond to THREE coarse-graining levels
# L_1 < L_2 < L_3. At each level, the block-averaged field is:
# Phi(L) ~ L^{-(d-2+eta)/2} * phi_0  (finite-size scaling)
# In d=3: Phi(L) ~ L^{-(1+eta)/2} ~ L^{-0.518}

# If g0_k ~ Phi(L_k) and L_k = phi^k * L_0 (golden ratio scaling):
# g0_k ~ L_k^{-0.518} ~ (phi^k)^{-0.518} * g0_0

# This gives: g0_mu/g0_e = phi^{-0.518} (NOT phi!)
# And g0_tau/g0_e = phi^{-1.036}

exponent = (1 + ETA) / 2
print(f"\n  Finite-size scaling: Phi(L) ~ L^{{-(d-2+eta)/2}} = L^{{-{exponent:.4f}}}")
print(f"  If L_k = phi^k * L_0:")
print(f"  g0_k/g0_0 ~ phi^{{-k*{exponent:.4f}}}")
print(f"  g0_mu/g0_e = phi^{{-{exponent:.4f}}} = {PHI_GOLD**(-exponent):.6f}")
print(f"  Actual: g0_mu/g0_e = phi = {PHI_GOLD:.6f}")
print(f"  Does NOT match! Exponent is -0.518, not +1.")

# Alternative: maybe L_k = phi^{k/exponent} * L_0?
# Then Phi(L_k) ~ phi^{-k} * phi_0 -> g0_k ~ phi^{-k} * g0_0
# But g0_mu = phi * g0_e, not phi^{-1} * g0_e.

# The sign is wrong: larger L gives SMALLER Phi, but we need
# g0_mu > g0_e (heavier = larger g0 = closer to collapse).

# If smaller L (more UV, higher resolution) -> larger g0:
# g0_k ~ L_k^{+exponent} * const -> g0 grows with smaller L

# The hierarchy would be: tau = smallest block (most UV),
# electron = largest block (most IR).

print(f"\n  Reversed: UV (small L) = heavy generation, IR (large L) = light")
print(f"  g0_k ~ L_k^{{+{exponent:.4f}}}")
print(f"  For phi-ratio: L_tau/L_e = (g0_tau/g0_e)^{{1/{exponent:.4f}}}")

ratio_g0 = 1.569627 / 0.86770494
print(f"  g0_tau/g0_e = {ratio_g0:.4f}")
print(f"  L_tau/L_e = {ratio_g0**(1/exponent):.4f}")
print(f"  phi^2 = {PHI_GOLD**2:.4f}")
print(f"  phi^{1/exponent:.2f} = {PHI_GOLD**(1/exponent):.4f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 4: CRITICAL EXPONENTS AND CV")
print("=" * 70)

# Another approach: the mass spectrum of excitations at the WF point.
# In the ordered phase (T < T_c), the field has:
# - Goldstone mode (massless if continuous symmetry, but Ising is discrete)
# - Higgs mode (massive, m_H ~ xi^{-1})
# - Higher modes from the lattice

# But Ising has only Z2, no Goldstone. So the single critical mode has:
# m ~ (T_c - T)^nu ~ xi^{-1}

# The RATIO of moments at different scales gives:
# <phi^{2n}>/<phi^2>^n = universal numbers (Ising universality class)

# Key universal ratios for 3D Ising (from conformal bootstrap):
# <phi^6>/<phi^2>^3 ~ 15 * A_6 (with A_6 universal)
# Let me use the known cumulant ratios.

print(f"\n  Universal moment ratios at WF (3D Ising):")
print(f"  <phi^2>:  defines the scale")
print(f"  <phi^4>/<phi^2>^2 = {ratio_42:.4f} (kurtosis ratio)")

# From conformal bootstrap / MC:
# The ratio <phi^6>/<phi^2>^3 is related to the sextic coupling
# At WF: g_6*/g_4*^2 is known but depends on normalization.
# Let's use the connected correlators instead.

# A different approach: the MAGNETIZATION distribution at T_c
# P(M) for a finite system of size L^d at T_c has a known universal shape.
# <|M|^k> / <|M|>^k are universal.
# For 3D Ising: <|M|^2>/<|M|>^2 ~ 1.13 (Binder paper)

# The ratio R = <|M|^2>/<|M|>^2 = pi/(3*sqrt(3)) + O(1/L) ???
# Actually, for the magnetization distribution at T_c:
# P(m) ~ |m|^{(beta/nu)} at small m

# Let me compute from known exponents:
# <|M|^k> ~ L^{-k*beta/nu}  (finite-size scaling)
# beta/nu = (d - 2 + eta)/2 = (1 + eta)/2 = 0.5182

bnu = BETA_C / NU  # = 0.3265/0.6300 = 0.5183
print(f"\n  beta/nu = {bnu:.4f}")

# For THREE generations at scales L_1, L_2, L_3:
# m_k ~ A_k^4 ~ L_k^{-4*beta/nu}  (if A_tail scales like magnetization)
# Then sqrt(m_k) ~ L_k^{-2*beta/nu}
# CV(sqrt(m)) = std/mean of {L_k^{-2*bnu}}

# If L_k = phi^{k-1} * L_0 (for k=1,2,3):
# sqrt(m_k) ~ phi^{-2*bnu*(k-1)} = phi^{-1.037*(k-1)}

# Then x_k = phi^{-1.037*(k-1)} for k=0,1,2
x_vals = [PHI_GOLD**(-2*bnu*k) for k in range(3)]
cv_x = CV(x_vals)
print(f"\n  Test: L_k = phi^k, x_k = L_k^{{-2*beta/nu}}:")
print(f"  x = [{x_vals[0]:.6f}, {x_vals[1]:.6f}, {x_vals[2]:.6f}]")
print(f"  CV(x) = {cv_x:.6f}")
print(f"  Target CV = 1.0000")
print(f"  Deviation: {(cv_x - 1.0)*100:+.2f}%")

# What exponent alpha gives CV = 1?
# x_k = phi^{-alpha*k} for k=0,1,2
# CV = 1 means std(x) = mean(x)
# This is transcendental in alpha. Solve numerically.
def cv_from_alpha(alpha):
    x = [PHI_GOLD**(-alpha*k) for k in range(3)]
    return CV(x) - 1.0

from scipy.optimize import brentq as brentq_opt
alpha_cv1 = brentq_opt(cv_from_alpha, 0.1, 10.0)
print(f"\n  CV = 1 requires alpha = {alpha_cv1:.6f}")
print(f"  2*beta/nu = {2*bnu:.6f}")
print(f"  alpha/2*beta/nu = {alpha_cv1 / (2*bnu):.6f}")
print(f"  Expected from TGP: m ~ A^4, A ~ phi^k -> sqrt(m) ~ phi^{2*1}*k")
print(f"  In TGP: alpha should = 2 (because sqrt(m) ~ A^2 ~ g0^2 ~ phi^2k)")

# With alpha = 2:
x_alpha2 = [PHI_GOLD**(-2*k) for k in range(3)]
cv_alpha2 = CV(x_alpha2)
print(f"\n  With alpha = 2 (TGP mass law):")
print(f"  x = [{x_alpha2[0]:.6f}, {x_alpha2[1]:.6f}, {x_alpha2[2]:.6f}]")
print(f"  CV = {cv_alpha2:.6f}")
print(f"  Deviation: {(cv_alpha2 - 1.0)*100:+.2f}%")

# The TGP mass law gives sqrt(m_k) ~ A_k^2, and A_k = A_tail(g0_k).
# From PATH 18b: g0_tau != phi^2 * g0_e.
# So pure geometric scaling x_k = phi^{-alpha*k} is WRONG for tau.
# The alpha that gives CV=1 is not achievable from pure phi-scaling.

# Let's check: what alpha gives CV = 1?
print(f"\n  The alpha for CV=1 is {alpha_cv1:.6f}.")
print(f"  Geometric hierarchy ratio: phi^alpha = {PHI_GOLD**alpha_cv1:.6f}")
print(f"  Actual mass ratios:")
g0_e_val = 0.86770494
g0_mu_val = PHI_GOLD * g0_e_val
Ae = A_tail(g0_e_val)
Am = A_tail(g0_mu_val)

def find_g0_tau_koide():
    def resid(g0t):
        At = A_tail(g0t)
        if At <= 0: return 10.0
        xe, xm, xt = Ae**2, Am**2, At**2
        return (xe+xm+xt)**2/(xe**2+xm**2+xt**2) - 1.5
    return brentq(resid, 1.50, 1.59)

g0_tau_val = find_g0_tau_koide()
At = A_tail(g0_tau_val)
print(f"  r_21 = m_mu/m_e = {(Am/Ae)**4:.2f}")
print(f"  r_31 = m_tau/m_e = {(At/Ae)**4:.2f}")
print(f"  phi^(4*alpha) = {PHI_GOLD**(4*alpha_cv1):.2f}")
print(f"  phi^8 = {PHI_GOLD**8:.2f}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 5: ISING MOMENT CONSTRAINT ON SOLITON MASSES")
print("=" * 70)

# The most direct connection would be:
# If soliton masses come from the PARTITION FUNCTION of the substrate,
# then the mass ratios are determined by the free energy landscape.

# In mean-field Ising near T_c:
# F(M) = (T-T_c)/T_c * M^2/2 + u* * M^4/24
# The minima are at M_0 = sqrt(6|t|/u*) where t = (T_c-T)/T_c

# The FLUCTUATIONS around the minimum have:
# <(dM)^2> = T / (|t| + 3*u*M_0^2/2) = T/(3|t|) near T_c

# For THREE generations, we'd need THREE minima (different magnetizations).
# But Ising has only TWO minima (+M_0, -M_0) related by Z2.
# Three generations from Z2 Ising requires SYMMETRY BREAKING.

# Possible resolution:
# 1. Three solitons sit at DIFFERENT positions in the soliton spectrum
#    (not three Ising minima).
# 2. The mass spectrum of solitons is determined by the ODE,
#    not directly by Ising moments.

# Let me test a concrete model:
# Suppose the "effective Binder cumulant" of the soliton mass distribution
# matches the Ising value.

# For 3 values {m_1, m_2, m_3}, define "soliton Binder":
# U_sol = 1 - <m^2>/(3*<m>^2) where <> = 1/3 sum
m_e = Ae**4
m_mu = Am**4
m_tau = At**4
m_arr = np.array([m_e, m_mu, m_tau])

U_sol = 1 - np.mean(m_arr**2) / (3 * np.mean(m_arr)**2)
print(f"\n  Soliton 'Binder cumulant':  U_sol = {U_sol:.6f}")
print(f"  Ising U4* = {U4_STAR:.6f}")

# For sqrt(m):
sq_m = np.sqrt(m_arr)
U_sol_sq = 1 - np.mean(sq_m**2) / (3 * np.mean(sq_m)**2)
print(f"\n  For sqrt(m) values:")
print(f"  U_sol(sqrt(m)) = {U_sol_sq:.6f}")
print(f"  Note: CV=1 <=> U_sol(sqrt(m)) = 1 - <x^2>/(3*<x>^2)")
print(f"       = 1 - (mean(x^2))/(3*mean(x)^2)")
R_sq = np.mean(sq_m**2) / np.mean(sq_m)**2
print(f"  <x^2>/<x>^2 = {R_sq:.6f}")
print(f"  CV(x)^2 + 1 = {R_sq:.6f}  -> CV = {np.sqrt(R_sq - 1):.6f}")
print(f"  (Should be 1.0000 for Koide)")

# For CV=1: <x^2>/<x>^2 = 2 -> U_sol = 1 - 2/3 = 1/3
print(f"\n  CV = 1 requires: U(sqrt(m)) = 1/3 = {1/3:.6f}")
print(f"  Actual: U(sqrt(m)) = {U_sol_sq:.6f}")
print(f"  Match: {abs(U_sol_sq - 1/3) < 0.001}")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 6: CLT ARGUMENT (from ex126 P9/P13)")
print("=" * 70)

# From ex126: r = sqrt(N-1) is interpreted as:
# "variance of N-1 independent Bernoulli-like variables"
# sigma^2 = r^2 = N-1

# This suggests: if sqrt(m_k) = c*(1 + X_k) where X_k is a sum of
# (N-1) = 2 independent contributions, then var(X_k) = N-1 = 2,
# and CV(1+X) = sqrt(var(X)) / (1 + mean(X)) = sqrt(2) / 1 = sqrt(2)
# Hmm, that gives CV = sqrt(2), not 1.

# Actually in Brannen: sqrt(m_k) = c*(1 + b*cos(theta + 2pi*k/3))
# The "mean" of sqrt(m) = c (since sum of cos terms = 0 for N=3)
# The "std" of sqrt(m) = c*b/sqrt(2) (from cos^2 average = 1/2, N=3 terms)
# So CV = (b/sqrt(2)) / 1 = b/sqrt(2)
# CV = 1 <=> b = sqrt(2)

# CLT interpretation: b = sqrt(2) = sqrt(N-1) for N=3
# As if each cos component carries unit variance, and there are N-1 = 2
# independent "directions" contributing.

# In the substrate: could these "directions" be the TWO independent
# fluctuation modes near T_c?
# In 3D Ising with Z2: there's ONE order parameter mode and
# the transverse (Goldstone-like, but discrete Z2 has none) modes.

# Actually: N-1 = 2 "directions" could be:
# - The two independent components of the cos/sin decomposition of the tail
# - The two polarization states of the oscillatory tail
# - Or simply a mathematical coincidence with N-1 degrees of freedom

print(f"\n  CLT argument: b = sqrt(N-1) = sqrt(2) for N=3")
print(f"  Interpretation: 'N-1 independent fluctuation modes'")
print(f"  In 3D Ising at WF:")
print(f"    - Z2 symmetry -> 1 order parameter mode")
print(f"    - In d=3: 3 spatial directions")
print(f"    - Soliton tail has 2 independent components (cos, sin)")
print(f"    - N-1 = 2 matches 'two tail components'")

# Test: for the soliton tail h(r) = A*sin(r-delta)/r
# = (A*cos(delta))*sin(r)/r - (A*sin(delta))*cos(r)/r
# Two independent components: a = A*cos(delta), b_comp = -A*sin(delta)
# A^2 = a^2 + b_comp^2 (two degrees of freedom)

# If these two components have equal independent fluctuations:
# <a^2> = <b^2> = sigma^2 (from substrate randomness)
# Then A^2 = a^2 + b^2 follows chi^2(2) distribution
# Chi^2(2) = Exponential distribution!
# And Exp distribution has CV = 1!

print(f"\n  KEY CONNECTION:")
print(f"  Soliton tail: h(r) = [a*sin(r) + b*cos(r)] / r")
print(f"  A^2 = a^2 + b^2 (Pythagoras)")
print(f"  If a, b are independent Gaussians: A^2 ~ chi^2(2) = Exp(1/sigma^2)")
print(f"  Exponential distribution has CV = 1!")
print(f"  This is EXACTLY the Koide condition!")
print(f"\n  But: a_k, b_k for the THREE solitons are NOT independent draws.")
print(f"  They are DETERMINISTIC functions of g0_k.")
print(f"  The 'random Gaussian' interpretation requires justification.")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 7: CHI-SQUARED HYPOTHESIS TEST")
print("=" * 70)

# If A_k^2 follows chi^2(2), then for k=1,2,3:
# The sample of 3 values should be consistent with Exp distribution.
# CV(Exp) = 1, which is what we see.
# But N=3 is too small for a real test.

# Instead: are the A^2 values consistent with Exp spacing?
# For Exp(lambda): the ordered statistics have known ratios.
# X_(1):X_(2):X_(3) for Exp ~ 1:1:1 in the differences.
# E[X_(1)] = 1/(3*lambda), E[X_(2)] = 1/(3*lambda) + 1/(2*lambda),
# E[X_(3)] = 1/(3*lambda) + 1/(2*lambda) + 1/lambda = 11/(6*lambda)

A2_sorted = np.sort([Ae**2, Am**2, At**2])  # sorted sqrt(m) values
diffs = np.diff(A2_sorted)
print(f"\n  Sorted A^2 values (= sqrt(m)):")
print(f"  A_e^2  = {A2_sorted[0]:.8f}")
print(f"  A_mu^2 = {A2_sorted[1]:.8f}")
print(f"  A_tau^2 = {A2_sorted[2]:.8f}")
print(f"\n  Spacings:")
print(f"  A_mu^2 - A_e^2  = {diffs[0]:.8f}")
print(f"  A_tau^2 - A_mu^2 = {diffs[1]:.8f}")
print(f"  Ratio: {diffs[1]/diffs[0]:.4f}")
print(f"  For Exp order stats: expected ratio of spacings = {(1/2)/(1/3):.4f} (i.e., 1.5)")
print(f"  For uniform spacings: ratio = 1.0")
print(f"  Actual ratio: {diffs[1]/diffs[0]:.4f}")

# Exp spacing test: the differences D_k = X_(k) - X_(k-1) should have
# D_k ~ Exp((N-k+1)*lambda). So D_1/D_2 ~ (3/2) * Exp/Exp
# This is too noisy with N=3.

# More meaningful: does the chi^2(2) interpretation work mechanistically?
# a_k = A_k * cos(delta_k), b_k = A_k * sin(delta_k)
# If a_k, b_k are "effectively random" (from substrate noise):
# A_k^2 = a_k^2 + b_k^2 ~ chi^2(2)

# The phases delta_k are DIFFERENT for each soliton (from ODE).
# So the a,b decomposition DOES have two independent components.
# But the MAGNITUDES A_k are set by g0_k through the ODE.

# The question reduces to: is the ODE map g0 -> A_tail such that
# the (a, b) components at the three g0 values "look like" Gaussian draws?

# Compute a_k, b_k for each soliton
from scipy.integrate import solve_ivp as sivp

def extract_ab(g0, r_min=50, r_max=300):
    r, g = solver_A(g0, r_max=r_max+50)
    mask = (r > r_min) & (r < r_max)
    if np.sum(mask) < 100:
        return 0.0, 0.0
    rf = r[mask]
    h = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, h, rcond=None)[0]
    return bc[0], bc[1]  # a = coeff of cos, b = coeff of sin

a_e, b_e = extract_ab(g0_e_val)
a_mu, b_mu = extract_ab(g0_mu_val)
a_tau, b_tau = extract_ab(g0_tau_val)

print(f"\n  Tail components (a, b) where h = a*cos(r)/r + b*sin(r)/r:")
print(f"  {'':>6s} {'a':>12s} {'b':>12s} {'A=sqrt(a^2+b^2)':>16s} {'A^2':>12s}")
print(f"  {'e':>6s} {a_e:12.8f} {b_e:12.8f} {np.sqrt(a_e**2+b_e**2):16.8f} {a_e**2+b_e**2:12.8f}")
print(f"  {'mu':>6s} {a_mu:12.8f} {b_mu:12.8f} {np.sqrt(a_mu**2+b_mu**2):16.8f} {a_mu**2+b_mu**2:12.8f}")
print(f"  {'tau':>6s} {a_tau:12.8f} {b_tau:12.8f} {np.sqrt(a_tau**2+b_tau**2):16.8f} {a_tau**2+b_tau**2:12.8f}")

# Are the a, b values consistent with Gaussian(0, sigma)?
a_all = np.array([a_e, a_mu, a_tau])
b_all = np.array([b_e, b_mu, b_tau])
print(f"\n  Statistics of a-components: mean = {np.mean(a_all):.6f}, std = {np.std(a_all):.6f}")
print(f"  Statistics of b-components: mean = {np.mean(b_all):.6f}, std = {np.std(b_all):.6f}")
print(f"  Combined sigma from a: {np.std(a_all):.6f}")
print(f"  Combined sigma from b: {np.std(b_all):.6f}")
print(f"  (For Gaussian: sigma_a ~ sigma_b)")


# ================================================================
print("\n" + "=" * 70)
print("SECTION 8: THE EXPONENTIAL DISTRIBUTION MECHANISM")
print("=" * 70)

# SYNTHESIS: The chain of reasoning would be:
# 1. Substrate (3D Ising) at criticality has universal fluctuations
# 2. Soliton tails have 2 independent components (sin, cos)
# 3. Each component is a functional of the substrate configuration
# 4. At criticality, these functionals become "effectively Gaussian"
#    (by universality / CLT on the block-averaged field)
# 5. A^2 = a^2 + b^2 follows chi^2(2) = Exponential
# 6. Exponential distribution has CV = 1
# 7. CV = 1 <=> Q_K = 3/2 (Koide)

# The MISSING STEPS:
# - Step 4 is the hardest: WHY are a, b "effectively Gaussian"?
# - The soliton is a DETERMINISTIC solution, not a random draw.
# - The three generations have FIXED g0 values, not random.
# - But: if the g0 values come from phi-FP + threshold,
#   the resulting a, b HAPPEN to look like Gaussian draws.

# Test: is the (a, b) pattern consistent with i.i.d. Gaussian?
print(f"\n  Pattern in (a, b) space:")
for name, a, b in [("e", a_e, b_e), ("mu", a_mu, b_mu), ("tau", a_tau, b_tau)]:
    angle = np.degrees(np.arctan2(b, a))
    radius = np.sqrt(a**2 + b**2)
    print(f"  {name:>4s}: r = {radius:.6f}, theta = {angle:.1f} deg")

print(f"\n  For i.i.d. Gaussian (a,b):")
print(f"  - Angles should be uniformly distributed (not clustered)")
print(f"  - Radii should follow Rayleigh distribution")
print(f"  - With only N=3, this is not statistically testable")


# ================================================================
print("\n" + "=" * 70)
print("FINAL SYNTHESIS: PROPOSAL #10")
print("=" * 70)
print(f"""
  RESULTS:

  1. DIRECT MOMENT MATCH: The Ising Binder cumulant U4* = {U4_STAR}
     does NOT match the CV=1 requirement U4 = 1/3.
     The substrate phi^2 distribution at criticality has CV = {np.sqrt(ratio_42-1):.3f},
     not 1. Deviation: {(np.sqrt(ratio_42-1)-1)*100:+.1f}%.

  2. COARSE-GRAINING SCALES: Finite-size scaling gives Phi ~ L^{{-{exponent:.3f}}}.
     This does NOT match the phi-FP scaling (exponent should be ~1, not 0.518).

  3. CRITICAL EXPONENT MATCH: 2*beta/nu = {2*bnu:.4f}, while CV=1
     requires alpha = {alpha_cv1:.4f}. These differ significantly.

  4. **KEY INSIGHT -- CHI-SQUARED CONNECTION**:
     Soliton tails have TWO independent components (a*cos + b*sin).
     A^2 = a^2 + b^2.
     IF a, b are "effectively Gaussian" (from substrate randomness):
       A^2 ~ chi^2(2) = Exponential
       Exp distribution has CV = 1 = KOIDE!

     This would provide a STATISTICAL derivation of Q_K = 3/2:
     - 2 components (sin/cos of tail) -> chi^2(2) -> CV = 1
     - The "2" comes from d=3 solitons having oscillatory tails
       (for d=2 no oscillation, for d>3 more components?)

  5. PROBLEM: The three solitons are DETERMINISTIC, not random.
     The (a, b) values are fixed by g0_k through the ODE.
     The "effective Gaussian" interpretation requires:
     - Either a STATISTICAL argument from substrate averaging
     - Or a STRUCTURAL argument that the ODE naturally produces
       chi^2(2)-like spread in A^2 values

  STATUS: Proposal #10 does NOT provide a direct substrate->CV=1 derivation.
  But it reveals a potentially deep connection:
    TWO tail components -> chi^2(2) -> Exponential -> CV=1 -> Q_K=3/2
  This is the most promising structural link found so far, but requires
  justifying WHY deterministic soliton amplitudes behave as if they
  were draws from an exponential distribution.
""")
