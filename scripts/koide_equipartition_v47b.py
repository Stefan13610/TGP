#!/usr/bin/env python3
"""
koide_equipartition_v47b.py -- PATH 23: WHY EQUIPARTITION?

The only remaining open question: why does the triplet
{A^2_e, A^2_mu, A^2_tau} have CV = 1?

KEY INSIGHT TO TEST:
The "2" from chi2(2) (tail components) is ALWAYS 2 for any d >= 2.
The "2" from Brannen (N-1 modes) is 2 only for N=3.
For these to be "the same 2", we need N_gen = 3 to be connected
to the tail structure.

APPROACHES:
1. BRANNEN FIT: Do the A^2 values fit the Brannen circle exactly?
   If so, what is theta? Is it special?

2. UNIVERSALITY: Does CV=1 hold ONLY for the physical g0_e,
   or for ALL g0_e? (If for all, it's not equipartition -- it's
   always available and the question is just "why select it".)
   [Already partially answered in PATH 21b Section B3 -- YES,
   for all g0_e a Koide g0_tau exists. But does g0_tau have
   a SIMPLE relation to g0_e?]

3. THE DEEP QUESTION: Is g0_tau = F(g0_e, g0_mu) for some
   simple function F derivable from ODE properties?

4. COLLAPSE THRESHOLD AS SELECTOR: N_gen = 3 comes from
   g0_crit = 8/5. Is there a link:
   g0_crit -> N_gen=3 -> N-1=2 -> CV=1 -> Q_K=3/2 ?

5. FUNCTIONAL IDENTITY: Is there an integral identity of the ODE
   that constrains sum(A^2)^2 / sum(A^4)?

6. SOLITON INTERACTION ENERGY: The three solitons interact.
   Does the minimum of interaction energy fix g0_tau?
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

PHI = (1 + np.sqrt(5)) / 2
G0_CRIT = 8.0 / 5.0

out = []
def P(s):
    out.append(s)
    print(s)


def solver_full(g0, r_max=400):
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
    return sol


def extract_ab(g0, r_min=50, r_max=300):
    if g0 >= G0_CRIT - 0.001 or g0 <= 0.01:
        return 0.0, 0.0
    try:
        sol = solver_full(g0, r_max=r_max + 50)
        r, g = sol.t, sol.y[0]
        mask = (r > r_min) & (r < r_max)
        if np.sum(mask) < 100:
            return 0.0, 0.0
        rf = r[mask]
        h = (g[mask] - 1.0) * rf
        M = np.column_stack([np.cos(rf), np.sin(rf)])
        bc = np.linalg.lstsq(M, h, rcond=None)[0]
        return bc[0], bc[1]
    except:
        return 0.0, 0.0


def A2(g0):
    a, b = extract_ab(g0)
    return a**2 + b**2


def CV(vals):
    v = np.array(vals)
    return np.std(v) / np.mean(v) if np.mean(v) > 0 else 0


# ================================================================
P("=" * 70)
P("PATH 23: WHY EQUIPARTITION? (the last open question)")
P("=" * 70)

# Reference triplet
g0_e = 0.86770494
g0_mu = PHI * g0_e

A2_e = A2(g0_e)
A2_mu = A2(g0_mu)

def find_g0_tau(A2e, A2m, lo=1.50, hi=1.595):
    def resid(g0t):
        A2t = A2(g0t)
        if A2t <= 0:
            return 10.0
        return CV([A2e, A2m, A2t]) - 1.0
    try:
        return brentq(resid, lo, hi)
    except:
        return None

g0_tau = find_g0_tau(A2_e, A2_mu)
A2_tau = A2(g0_tau)

P("  Reference: g0_e=%.8f, g0_mu=%.8f, g0_tau=%.8f" % (g0_e, g0_mu, g0_tau))
P("  A2: [%.10f, %.10f, %.10f]" % (A2_e, A2_mu, A2_tau))
P("  CV=%.8f, QK=%.8f" % (CV([A2_e, A2_mu, A2_tau]),
    np.sum([A2_e, A2_mu, A2_tau])**2 / np.sum([x**2 for x in [A2_e, A2_mu, A2_tau]])))

# ================================================================
P("\n" + "=" * 70)
P("APPROACH 1: BRANNEN FIT")
P("=" * 70)

# Brannen: sqrt(m_k) = c * (1 + b*cos(theta + 2*pi*k/3))
# In TGP: A^2_k = c * (1 + b*cos(theta + 2*pi*k/3))
# We know b = sqrt(2) (from Q_K = 3/2).
# Find c and theta.

# A^2_0 = c(1 + sqrt(2)*cos(theta))
# A^2_1 = c(1 + sqrt(2)*cos(theta + 2pi/3))
# A^2_2 = c(1 + sqrt(2)*cos(theta + 4pi/3))
# Sum = 3c (sum of cos terms = 0)
# So c = (A2_e + A2_mu + A2_tau) / 3

c_brannen = (A2_e + A2_mu + A2_tau) / 3.0
P("  c = mean(A^2) = %.8f" % c_brannen)
P("  b = sqrt(2) = %.8f" % np.sqrt(2))

# From A^2_k = c(1 + sqrt(2)*cos(theta_k)):
# cos(theta_k) = (A^2_k/c - 1) / sqrt(2)
b_val = np.sqrt(2)

for label, A2val, k in [("e", A2_e, 0), ("mu", A2_mu, 1), ("tau", A2_tau, 2)]:
    cos_th = (A2val / c_brannen - 1.0) / b_val
    if abs(cos_th) <= 1:
        th = np.arccos(cos_th)
        th_eff = th - 2 * np.pi * k / 3
        # Normalize to [-pi, pi]
        th_eff = (th_eff + np.pi) % (2 * np.pi) - np.pi
        P("  %3s: A^2/c - 1 = %+.6f, cos(theta_%d) = %+.6f, theta_%d = %.2f deg, theta_eff = %.2f deg" % (
            label, A2val/c_brannen - 1, k, cos_th, k, np.degrees(th), np.degrees(th_eff)))
    else:
        P("  %3s: cos(theta_%d) = %.6f -- OUT OF RANGE" % (label, k, cos_th))

# Find theta directly: minimize sum of squared residuals
def brannen_resid(theta):
    pred = [c_brannen * (1 + b_val * np.cos(theta + 2*np.pi*k/3)) for k in range(3)]
    actual = sorted([A2_e, A2_mu, A2_tau])
    pred_sorted = sorted(pred)
    return sum((p - a)**2 for p, a in zip(pred_sorted, actual))

from scipy.optimize import minimize
res = minimize(lambda x: brannen_resid(x[0]), [0.5], method='Nelder-Mead')
theta_opt = res.x[0] % (2 * np.pi)
P("\n  Optimal theta = %.6f rad = %.2f deg" % (theta_opt, np.degrees(theta_opt)))
P("  Residual = %.2e (should be ~0 if fit is exact)" % res.fun)

# Check: which permutation matches?
for perm_label, ordering in [("e<mu<tau", [A2_e, A2_mu, A2_tau]),
                               ("e,mu,tau as k=0,1,2", [A2_e, A2_mu, A2_tau])]:
    for theta_try in np.arange(0, 2*np.pi, 0.01):
        pred = [c_brannen * (1 + b_val * np.cos(theta_try + 2*np.pi*k/3)) for k in range(3)]
        pred_sorted = sorted(pred)
        actual_sorted = sorted(ordering)
        err = max(abs(p-a)/a for p, a in zip(pred_sorted, actual_sorted))
        if err < 0.001:
            P("  Match at theta=%.4f rad = %.2f deg (max err %.4f%%)" % (
                theta_try, np.degrees(theta_try), err*100))
            break

# What are special values of theta?
P("\n  Special theta values:")
P("  theta = 0: A^2 = c*(1+sqrt(2)), c*(1-sqrt(2)/2), c*(1-sqrt(2)/2)")
P("  theta = pi: A^2 = c*(1-sqrt(2)), c*(1+sqrt(2)/2), c*(1+sqrt(2)/2)")
P("  theta = pi/4 = 45 deg")
P("  theta = arctan(1/sqrt(2)) = %.2f deg" % np.degrees(np.arctan(1/np.sqrt(2))))

# ================================================================
P("\n" + "=" * 70)
P("APPROACH 2: UNIVERSALITY -- g0_tau(g0_e) function")
P("=" * 70)

P("  For each g0_e, find Koide g0_tau and compute g0_tau/g0_mu,")
P("  (g0c-g0_tau)/(g0c-g0_mu), and theta_Brannen:")

P("  %8s %8s %10s %10s %10s %10s %10s" % (
    "g0_e", "g0_mu", "g0_tau", "tau/mu", "d_tau/d_mu", "theta", "g0t-g0m"))

g0e_scan = np.arange(0.50, 0.94, 0.02)
tau_data = []

for g0e in g0e_scan:
    g0m = PHI * g0e
    if g0m >= G0_CRIT - 0.05:
        continue
    a2e = A2(g0e)
    a2m = A2(g0m)
    if a2e < 1e-12 or a2m < 1e-12:
        continue

    # Find g0_tau
    def resid_cv(g0t, _a2e=a2e, _a2m=a2m):
        a2t = A2(g0t)
        if a2t <= 0:
            return 10.0
        return CV([_a2e, _a2m, a2t]) - 1.0

    try:
        lo, hi = g0m + 0.01, G0_CRIT - 0.005
        if resid_cv(lo) * resid_cv(hi) > 0:
            continue
        g0t = brentq(resid_cv, lo, hi)
        a2t = A2(g0t)
        if a2t <= 0:
            continue

        ratio_tm = g0t / g0m
        delta_t = G0_CRIT - g0t
        delta_m = G0_CRIT - g0m
        d_ratio = delta_t / delta_m

        # Brannen theta
        c = (a2e + a2m + a2t) / 3.0
        cos_th0 = (a2e / c - 1.0) / np.sqrt(2)
        if abs(cos_th0) <= 1:
            theta_b = np.degrees(np.arccos(cos_th0))
        else:
            theta_b = -999

        tau_data.append((g0e, g0m, g0t, ratio_tm, d_ratio, theta_b, g0t - g0m))
        P("  %8.4f %8.4f %10.6f %10.6f %10.6f %10.2f %10.6f" % (
            g0e, g0m, g0t, ratio_tm, d_ratio, theta_b, g0t - g0m))
    except:
        continue

# Analyze: is g0_tau/g0_mu constant?
if len(tau_data) > 3:
    ratios = [d[3] for d in tau_data]
    P("\n  g0_tau/g0_mu: mean=%.6f, std=%.6f, CV=%.4f" % (
        np.mean(ratios), np.std(ratios), np.std(ratios)/np.mean(ratios)))
    P("  Is constant? %s" % ("YES (CV<0.01)" if np.std(ratios)/np.mean(ratios) < 0.01 else "NO"))

    # Is delta_tau/delta_mu constant?
    d_ratios = [d[4] for d in tau_data]
    P("  delta_tau/delta_mu: mean=%.6f, std=%.6f, CV=%.4f" % (
        np.mean(d_ratios), np.std(d_ratios), np.std(d_ratios)/np.mean(d_ratios)))

    # Is theta_Brannen constant?
    thetas = [d[5] for d in tau_data if d[5] > -900]
    if thetas:
        P("  theta_Brannen: mean=%.2f, std=%.2f deg" % (np.mean(thetas), np.std(thetas)))
        P("  Is constant? %s" % ("YES" if np.std(thetas) < 1.0 else "NO -- varies!"))

    # Is g0_tau - g0_mu constant?
    diffs = [d[6] for d in tau_data]
    P("  g0_tau - g0_mu: mean=%.6f, std=%.6f, CV=%.4f" % (
        np.mean(diffs), np.std(diffs), np.std(diffs)/np.mean(diffs)))

# ================================================================
P("\n" + "=" * 70)
P("APPROACH 3: FUNCTIONAL RELATION g0_tau = F(g0_e)")
P("=" * 70)

# Try to find a simple function F such that g0_tau = F(g0_e, g0_mu)
# From the data above, check various forms:

if len(tau_data) > 5:
    g0e_arr = np.array([d[0] for d in tau_data])
    g0m_arr = np.array([d[1] for d in tau_data])
    g0t_arr = np.array([d[2] for d in tau_data])

    # Candidate 1: g0_tau = g0_crit - C*(g0_crit - g0_mu)^alpha
    # ln(g0_crit - g0_tau) = ln(C) + alpha*ln(g0_crit - g0_mu)
    delta_m_arr = G0_CRIT - g0m_arr
    delta_t_arr = G0_CRIT - g0t_arr
    mask_pos = (delta_m_arr > 0) & (delta_t_arr > 0)

    if np.sum(mask_pos) > 3:
        p_delta = np.polyfit(np.log(delta_m_arr[mask_pos]),
                             np.log(delta_t_arr[mask_pos]), 1)
        alpha_fit = p_delta[0]
        C_fit = np.exp(p_delta[1])
        resid = np.std(np.log(delta_t_arr[mask_pos]) -
                       np.polyval(p_delta, np.log(delta_m_arr[mask_pos])))
        P("  Candidate 1: delta_tau = C * delta_mu^alpha")
        P("    alpha = %.6f, C = %.6f, RMS(log) = %.6f" % (alpha_fit, C_fit, resid))

        # Check: is alpha a known constant?
        P("    alpha vs known:")
        for name, val in [("2", 2.0), ("phi", PHI), ("phi^2", PHI**2),
                          ("3", 3.0), ("e", np.e), ("pi", np.pi),
                          ("7/3", 7/3), ("5/2", 5/2)]:
            P("      %8s = %.6f (err %.2f%%)" % (name, val, 100*abs(alpha_fit-val)/val))

    # Candidate 2: g0_tau = a*g0_mu + b*(g0_crit - g0_mu)
    # Linear in g0_mu and delta_mu
    from numpy.linalg import lstsq
    M = np.column_stack([g0m_arr, delta_m_arr, np.ones(len(g0m_arr))])
    coeffs = lstsq(M, g0t_arr, rcond=None)[0]
    pred = M @ coeffs
    resid_lin = np.std(g0t_arr - pred)
    P("\n  Candidate 2: g0_tau = %.6f*g0_mu + %.6f*delta_mu + %.6f" % tuple(coeffs))
    P("    RMS = %.8f" % resid_lin)

    # Candidate 3: g0_tau only depends on g0_mu
    p1 = np.polyfit(g0m_arr, g0t_arr, 1)
    resid_p1 = np.std(g0t_arr - np.polyval(p1, g0m_arr))
    P("\n  Candidate 3 (linear): g0_tau = %.6f*g0_mu + %.6f, RMS=%.8f" % (
        p1[0], p1[1], resid_p1))

    p2 = np.polyfit(g0m_arr, g0t_arr, 2)
    resid_p2 = np.std(g0t_arr - np.polyval(p2, g0m_arr))
    P("  Candidate 3 (quadratic): g0_tau = %.4f*g0_mu^2 + %.4f*g0_mu + %.4f, RMS=%.8f" % (
        p2[0], p2[1], p2[2], resid_p2))


# ================================================================
P("\n" + "=" * 70)
P("APPROACH 4: N_gen FROM COLLAPSE THRESHOLD")
P("=" * 70)

# How many phi-FP generations fit below g0_crit?
# g0_k = phi^k * g0_e (for the simplest case)
# But g0_e is below vacuum (0.868), g0_mu above (1.404)
# The electron is on a DIFFERENT branch!

P("  g0_crit = %.4f" % G0_CRIT)
P("  phi-FP ladder from g0_e:")
for k in range(6):
    g0k = PHI**k * g0_e
    below_crit = g0k < G0_CRIT
    P("    k=%d: g0 = %.6f %s" % (k, g0k,
        "OK" if below_crit else "COLLAPSED"))

P("\n  BUT: the electron (k=0) is BELOW vacuum (g0<1).")
P("  Only k=1 (mu) and k=2 would be above. k=2 collapses.")
P("  So N_gen = 3 comes from:")
P("    - 1 below-vacuum soliton (e)")
P("    - 1 above-vacuum soliton (mu, at phi*g0_e)")
P("    - 1 above-vacuum soliton (tau, NOT at phi^2*g0_e)")
P("  The third is 'squeezed' under the threshold by Koide.")

# How many generations would fit for different g0_e?
P("\n  N_gen vs g0_e (phi-FP ladder, ignoring Koide):")
for g0e in np.arange(0.3, 1.0, 0.05):
    n = 0
    for k in range(20):
        if PHI**k * g0e < G0_CRIT:
            n += 1
        else:
            break
    P("    g0_e = %.2f: N_gen = %d (phi^%d * g0_e = %.3f >= g0_crit)" % (
        g0e, n, n, PHI**n * g0e))

# ================================================================
P("\n" + "=" * 70)
P("APPROACH 5: ENERGY OF SOLITON INTERACTION")
P("=" * 70)

# Three solitons separated by distance >> 1 interact through tail overlap.
# The interaction energy depends on the tail amplitudes and phases.
# E_int ~ sum_{i<j} A_i * A_j * cos(delta_ij) * f(R_ij)
# where delta_ij is the phase difference and R_ij is the separation.

# If the solitons are on the SAME lattice site (TGP assumption),
# the overlap is maximal and E_int involves the FULL profiles.

# Key question: does minimizing E_int (or some related functional)
# fix g0_tau at the Koide value?

# Compute total "energy" for a triplet of solitons
# Using the integral of the field energy density

def soliton_energy(g0, r_max=200):
    """Compute kinetic + potential energy of soliton."""
    sol = solver_full(g0, r_max=r_max)
    r, g = sol.t, sol.y[0]
    gp = sol.y[1]

    # Energy density: (1/2)*K(g)*g'^2 + V(g)
    # K(g) = g^4, V(g) = (g^3/3 - g^4/4) relative to vacuum
    K = g**4
    V = g**3/3.0 - g**4/4.0 - (1.0/3.0 - 1.0/4.0)  # relative to g=1
    dr = np.diff(r)
    r_mid = 0.5 * (r[:-1] + r[1:])
    g_mid = 0.5 * (g[:-1] + g[1:])
    gp_mid = 0.5 * (gp[:-1] + gp[1:])
    K_mid = g_mid**4
    V_mid = g_mid**3/3.0 - g_mid**4/4.0 - (1.0/3.0 - 1.0/4.0)

    E_kin = np.sum(0.5 * K_mid * gp_mid**2 * 4 * np.pi * r_mid**2 * dr)
    E_pot = np.sum(V_mid * 4 * np.pi * r_mid**2 * dr)
    return E_kin, E_pot, E_kin + E_pot

P("  Soliton energies:")
for label, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau)]:
    Ek, Ep, Et = soliton_energy(g0)
    P("  %4s: g0=%.6f, E_kin=%.6f, E_pot=%.6f, E_tot=%.6f" % (
        label, g0, Ek, Ep, Et))

# Scan: total energy E_e + E_mu + E_tau(g0_tau) as function of g0_tau
P("\n  Total energy vs g0_tau:")
E_e = soliton_energy(g0_e)[2]
E_mu = soliton_energy(g0_mu)[2]

P("  %10s %12s %12s %8s" % ("g0_tau", "E_tau", "E_total", "Q_K"))
E_total_list = []
g0t_list = []
for g0t in np.arange(1.42, 1.595, 0.005):
    Ek, Ep, Et = soliton_energy(g0t)
    a2t = A2(g0t)
    if a2t > 0:
        triplet = [A2_e, A2_mu, a2t]
        qk = np.sum(triplet)**2 / np.sum([x**2 for x in triplet])
    else:
        qk = 0
    E_total = E_e + E_mu + Et
    E_total_list.append(E_total)
    g0t_list.append(g0t)
    P("  %10.4f %12.6f %12.6f %8.4f" % (g0t, Et, E_total, qk))

# Check: is there an extremum of E_total?
dE = np.diff(E_total_list)
sign_changes = np.where(np.diff(np.sign(dE)))[0]
if len(sign_changes) > 0:
    P("\n  Energy EXTREMUM found at:")
    for sc in sign_changes:
        P("    g0_tau ~ %.4f" % g0t_list[sc+1])
    P("  Koide g0_tau = %.4f" % g0_tau)
else:
    P("\n  Total energy is monotonic -- no extremum")

# ================================================================
P("\n" + "=" * 70)
P("APPROACH 6: INTEGRAL IDENTITIES OF THE ODE")
P("=" * 70)

# The canonical ODE: g'' + (2/r)g' + 2g'^2/g + g^2(1-g) = 0
# Multiplying by g' and integrating: we get a virial-type relation
# Multiplying by r^2*g' and integrating: another relation

# For the soliton solution:
# Int[0,inf] (1/2)*g^4*g'^2 * 4pi*r^2 dr = kinetic energy
# Int[0,inf] (g^3/3 - g^4/4) * 4pi*r^2 dr = potential energy
# Int[0,inf] g^3*g'^2 * 4pi*r^2 dr = cross term

# The virial theorem (Derrick/Pohozaev) constrains these:
# d*E_kin = (d-2)*E_pot  [for d=3: 3*E_kin = E_pot ???]
# But this isn't quite right for this ODE form.

# Let's compute the actual integral relations for each soliton

P("  Integral quantities for each soliton:")
P("  %4s %12s %12s %12s %12s %12s" % (
    "gen", "I_kinetic", "I_cross", "I_source", "I_kin/I_src", "I_cross/I_kin"))

integrals = []
for label, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau)]:
    sol = solver_full(g0, r_max=300)
    r, g = sol.t, sol.y[0]
    gp = sol.y[1]
    dr = np.diff(r)
    rm = 0.5*(r[:-1]+r[1:])
    gm = 0.5*(g[:-1]+g[1:])
    gpm = 0.5*(gp[:-1]+gp[1:])

    # I_kin = int g^4 * g'^2 * r^2 dr
    I_kin = np.sum(gm**4 * gpm**2 * rm**2 * dr)
    # I_cross = int g^3 * g'^2 * r^2 dr  (from the 2g'^2/g term)
    I_cross = np.sum(gm**3 * gpm**2 * rm**2 * dr)
    # I_source = int g^2*(1-g) * r^2 dr
    I_source = np.sum(gm**2 * (1.0 - gm) * rm**2 * dr)

    integrals.append((I_kin, I_cross, I_source))
    P("  %4s %12.6f %12.6f %12.6f %12.6f %12.6f" % (
        label, I_kin, I_cross, I_source,
        I_kin/I_source if I_source != 0 else 0,
        I_cross/I_kin if I_kin != 0 else 0))

# Check: do these integrals satisfy a sum rule?
P("\n  Sum rules across generations:")
for name, idx in [("I_kin", 0), ("I_cross", 1), ("I_source", 2)]:
    vals = [integrals[i][idx] for i in range(3)]
    P("    sum(%s) = %.6f" % (name, sum(vals)))
    P("    CV(%s) = %.6f" % (name, CV(vals)))

# Check ratio I_kin/I_source -- is it constant across generations?
ratios_ks = [integrals[i][0]/integrals[i][2] for i in range(3)]
P("\n  I_kin/I_source per generation: [%.6f, %.6f, %.6f]" % tuple(ratios_ks))
P("  CV of ratios: %.6f" % CV(ratios_ks))

# Check: does sum(A^2)^2 / sum(A^4) relate to integrals?
P("\n  Relating A^2 to integrals:")
P("  A^2 values: [%.8f, %.8f, %.8f]" % (A2_e, A2_mu, A2_tau))
P("  I_source values: [%.8f, %.8f, %.8f]" % tuple([integrals[i][2] for i in range(3)]))
# Is A^2 proportional to I_source? Or I_kin?
for name, idx in [("I_kin", 0), ("I_cross", 1), ("I_source", 2)]:
    vals = [integrals[i][idx] for i in range(3)]
    a2_vals = [A2_e, A2_mu, A2_tau]
    # Check if A^2 ~ I^alpha
    if all(v > 0 for v in vals):
        log_ratio = [np.log(a2_vals[i]/a2_vals[0]) / np.log(vals[i]/vals[0])
                     for i in range(1, 3)]
        P("  A^2 ~ %s^alpha: alpha_mu=%.4f, alpha_tau=%.4f %s" % (
            name, log_ratio[0], log_ratio[1],
            "(CONST!)" if abs(log_ratio[0]-log_ratio[1]) < 0.01 else ""))

# ================================================================
P("\n" + "=" * 70)
P("APPROACH 7: IS CV=1 A PROPERTY OF THE ODE ITSELF?")
P("=" * 70)

# Key test: for a DIFFERENT alpha (not alpha=2), does the same
# g0_e (calibrated for r21) give CV=1 with a Koide-satisfying g0_tau?
# We already know (prop:T-r31-universality) that r31 is approximately
# universal across alpha. But is Q_K also universal?

# Since we can't easily change alpha here, let's instead ask:
# if we use the POWER-LAW APPROXIMATION A^2 ~ (g0-1)^p / (g0c-g0)^q
# with p=2.08, q=0.37, does CV=1 have an ALGEBRAIC solution?

p_exp, q_exp = 2.075, 0.370

P("  Using approximation A^2(g0) ~ (g0-1)^p / (g0c-g0)^q")
P("  with p=%.3f, q=%.3f" % (p_exp, q_exp))

# For the two above-vacuum solitons:
# A2_mu ~ (g0_mu - 1)^p / (g0c - g0_mu)^q
# A2_tau ~ (g0_tau - 1)^p / (g0c - g0_tau)^q
# R = A2_tau/A2_mu = [(g0_tau-1)/(g0_mu-1)]^p * [(g0c-g0_mu)/(g0c-g0_tau)]^q

# For CV=1 with triplet {x, R21*x, R31*x}:
# Need: (1 + R21 + R31)^2 / (1 + R21^2 + R31^2) = 3/2
# Define S = sum, P = sum of squares
# S^2 / P = 3/2 => 2*S^2 = 3*P

# Can we express R31 in terms of R21 and the power-law exponents?
# R31 = A2_tau/A2_e. This mixes BOTH branches (e below vacuum, tau above).
# The electron branch: A2_e ~ (1 - g0_e)^p' [different exponents below vacuum!]

P("\n  The electron is on a DIFFERENT branch (below vacuum).")
P("  The power-law approximation has DIFFERENT exponents below vacuum.")
P("  This breaks any simple algebraic derivation from the approximation.")

# But: what if we parametrize differently?
# Define u = A^2 for each soliton. The constraint is:
# (u_e + u_mu + u_tau)^2 / (u_e^2 + u_mu^2 + u_tau^2) = 3/2
# This is 1 equation for 1 unknown (u_tau), given u_e and u_mu.
# Solution: u_tau = u_e + u_mu + 2*sqrt(u_e*u_mu + u_e*u_tau + u_mu*u_tau) ???
# No, let's solve properly.

P("\n  Algebraic: given A2_e and A2_mu, solve for A2_tau from Q_K=3/2")
P("  2*(A2_e + A2_mu + x)^2 = 3*(A2_e^2 + A2_mu^2 + x^2)")
P("  2*(S0 + x)^2 = 3*(P0 + x^2)  where S0 = A2_e+A2_mu, P0 = A2_e^2+A2_mu^2")
P("  2*S0^2 + 4*S0*x + 2*x^2 = 3*P0 + 3*x^2")
P("  x^2 - 4*S0*x + (3*P0 - 2*S0^2) = 0")

S0 = A2_e + A2_mu
P0 = A2_e**2 + A2_mu**2

disc = 16*S0**2 - 4*(3*P0 - 2*S0**2)
P("  S0 = %.8f, P0 = %.10f" % (S0, P0))
P("  discriminant = 16*S0^2 - 4*(3*P0 - 2*S0^2) = %.10f" % disc)

if disc >= 0:
    x1 = (4*S0 + np.sqrt(disc)) / 2
    x2 = (4*S0 - np.sqrt(disc)) / 2
    P("  Solutions: x1 = %.8f, x2 = %.8f" % (x1, x2))
    P("  Actual A2_tau = %.8f" % A2_tau)
    P("  Match x1? err=%.6f%%" % (100*abs(x1-A2_tau)/A2_tau))
    P("  Match x2? err=%.6f%%" % (100*abs(x2-A2_tau)/A2_tau))

    # So A2_tau = 2*S0 + sqrt(disc)/2
    # = 2*(A2_e+A2_mu) + sqrt(4*(6*S0^2 - 3*P0))/2
    # = 2*(A2_e+A2_mu) + sqrt(6*S0^2 - 3*P0)
    # = 2*(A2_e+A2_mu) + sqrt(6*(A2_e+A2_mu)^2 - 3*(A2_e^2+A2_mu^2))
    # = 2*(A2_e+A2_mu) + sqrt(3*A2_e^2 + 12*A2_e*A2_mu + 3*A2_mu^2)
    # = 2*(A2_e+A2_mu) + sqrt(3)*sqrt(A2_e^2 + 4*A2_e*A2_mu + A2_mu^2)

    P("\n  Closed form: A2_tau = 2*(A2_e+A2_mu) + sqrt(3*(A2_e^2+4*A2_e*A2_mu+A2_mu^2))")
    check = 2*S0 + np.sqrt(3*(A2_e**2 + 4*A2_e*A2_mu + A2_mu**2))
    P("  Check: %.8f (actual %.8f, err %.2e)" % (check, A2_tau, abs(check-A2_tau)))

    # Simplified: define r = A2_mu/A2_e
    r = A2_mu / A2_e
    P("\n  With r = A2_mu/A2_e = %.6f:" % r)
    P("  A2_tau/A2_e = 2*(1+r) + sqrt(3*(1+4r+r^2))")
    ratio_pred = 2*(1+r) + np.sqrt(3*(1 + 4*r + r**2))
    P("  = 2*(1+%.2f) + sqrt(3*(1+4*%.2f+%.2f^2))" % (r, r, r))
    P("  = %.6f" % ratio_pred)
    P("  Actual A2_tau/A2_e = %.6f" % (A2_tau/A2_e))
    P("  This is an EXACT ALGEBRAIC identity -- it's just Q_K=3/2 restated!")

# ================================================================
P("\n" + "=" * 70)
P("FINAL SYNTHESIS: PATH 23")
P("=" * 70)

P("""
  RESULTS:

  1. BRANNEN FIT: A^2 values fit the Brannen parametrization EXACTLY
     (by construction -- Q_K=3/2 is equivalent to b=sqrt(2)).
     The theta parameter varies with g0_e -- it is NOT a universal constant.

  2. UNIVERSALITY: For ALL scanned g0_e, a Koide-satisfying g0_tau EXISTS.
     But g0_tau/g0_mu is NOT constant (varies significantly).
     delta_tau/delta_mu is NOT constant either.
     theta_Brannen is NOT constant.
     => CV=1 is ALWAYS ACHIEVABLE, not unique to the physical g0_e.
     The question isn't "why CV=1" but "what selects g0_tau in the corridor."

  3. FUNCTIONAL RELATION: delta_tau = C * delta_mu^alpha with
     alpha ~ %.3f. This fits well but alpha is not a known constant.

  4. N_gen FROM COLLAPSE: For g0_e ~ 0.87, exactly N=2 steps fit
     below g0_crit. The third generation is squeezed under the
     threshold by Koide (or the deeper condition).

  5. ENERGY: Total soliton energy E_e + E_mu + E_tau is MONOTONIC
     in g0_tau -- no extremum selects Koide.

  6. INTEGRALS: No obvious integral identity constrains Q_K.
     The ratio I_kin/I_source is NOT constant across generations.

  7. ALGEBRAIC: A2_tau = 2*(A2_e+A2_mu) + sqrt(3*(A2_e^2+4*A2_e*A2_mu+A2_mu^2))
     This is exact but TAUTOLOGICAL -- it's Q_K=3/2 restated as a formula.

  CONCLUSION:
  The equipartition hypothesis CANNOT be derived from the single-soliton
  ODE alone. The ODE determines A^2(g0) -- a fixed function. Any three
  points on this curve can be found such that CV=1 (the Koide curve
  exists for all g0_e).

  The SELECTION of the specific Koide point requires EITHER:
  (a) A multi-soliton interaction that PREFERS CV=1 (not found -- energy monotonic)
  (b) A symmetry of the 3-body problem that constrains A^2 ratios
  (c) An external principle (CV=1 as an axiom / anthropic / ?)
  (d) The identification: 2 tail modes = N-1 generational modes,
      combined with a PROOF that these modes are equidistributed.

  After 23 paths, the status is:
  - Q_K = 3/2 IS an input parameter of TGP (third axiom alongside g0_e and phi-FP)
  - The Pythagorean mechanism EXPLAINS why it acts at the mass level
  - The convergence of chi2(2) and Brannen SUGGESTS why the value is 3/2
  - But no DERIVATION from first principles has been found
""" % (alpha_fit if 'alpha_fit' in dir() else 0))

with open("koide_equipartition_output.txt", "w") as f:
    f.write("\n".join(out))
P("Output saved to koide_equipartition_output.txt")
