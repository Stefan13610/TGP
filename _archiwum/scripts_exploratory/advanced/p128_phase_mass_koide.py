#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p128_phase_mass_koide.py -- Phase-corrected mass formula & Koide geometry
=========================================================================

Builds on p127 phase map results:
  delta(e)  = -81.14 deg
  delta(mu) = +43.84 deg
  Delta(e->mu) = 124.98 deg ~ 2*pi/3 (120 deg)

This script:
  1. Detailed phase map with fine grid near critical regions
  2. Test S_kin (kinetic action) as mass proxy instead of A_tail^4
  3. Derive phase-corrected mass formula: m ~ A^4 * F(delta)
  4. Koide condition Q = 3/2 as geometric constraint in (A, delta)
  5. Predict g0_tau from Koide + phase constraint simultaneously
  6. Stability analysis of all results

Author: TGP project, session v42+
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 2
R21_PDG = 206.768
R31_PDG = 3477.48

# Physical mass ratios
M_E = 0.51099895   # MeV
M_MU = 105.6583755 # MeV
M_TAU = 1776.86     # MeV

# Koide parameter
SQRT_E = np.sqrt(M_E)
SQRT_MU = np.sqrt(M_MU)
SQRT_TAU = np.sqrt(M_TAU)
Q_KOIDE = (M_E + M_MU + M_TAU) / (SQRT_E + SQRT_MU + SQRT_TAU)**2
print("  Physical Koide Q = %.8f (exact 2/3: %.8f)" % (Q_KOIDE, 2.0/3))

# phi-FP values
G0_E = 0.8992655880
G0_MU = PHI * G0_E  # = 1.45504

def solve_soliton(g0, rm=300):
    """Solve soliton ODE, return (r, g) arrays."""
    def fk(g):
        return 1 + 2*ALPHA*np.log(g) if g > 0 else -1e30
    def Vp(g):
        return g**2*(1-g)
    fg0 = fk(g0)
    if abs(fg0) < 1e-15:
        return None, None
    c2 = Vp(g0)/(3*fg0)
    rs = 0.01
    def rhs(r, y):
        g, p = y
        if g <= 1e-15: return [p, 0]
        fg = fk(g)
        if abs(fg) < 1e-10: return [p, 0]
        if r < 1e-10: return [p, Vp(g)/fg/3]
        return [p, (Vp(g)-2/r*p)/fg]
    def ev(r, y): return 100-abs(y[0])
    ev.terminal = True
    s = solve_ivp(rhs, [rs, rm], [g0+c2*rs**2, 2*c2*rs],
                  method='RK45', rtol=1e-11, atol=1e-13,
                  max_step=0.05, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], rm), 25000)
    return r, s.sol(r)[0]


def extract_tail(r, g, r_min=120, r_max=260):
    """Extract B_tail, C_tail from oscillatory tail."""
    m = (r >= r_min) & (r <= r_max)
    rf = r[m]
    tl = (g[m] - 1) * rf
    if len(rf) < 10:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, res, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    B, C = coeff
    A_tail = np.sqrt(B**2 + C**2)
    delta = np.arctan2(C, B)
    fit = A @ coeff
    rms = np.sqrt(np.mean((tl - fit)**2))
    return B, C, A_tail, delta, rms


def compute_Skin(r, g):
    """Kinetic action: S_kin = 4*pi * int r^2 * (1/2)*f(g)*(g')^2 dr"""
    dr = r[1] - r[0]
    gp = np.gradient(g, dr)
    fg = np.where(g > 0, 1 + 2*ALPHA*np.log(g), 0.0)
    return 4*np.pi * np.trapezoid(r**2 * 0.5 * fg * gp**2, r)


def compute_Etot(r, g):
    """Total soliton energy: kinetic + potential."""
    dr = r[1] - r[0]
    gp = np.gradient(g, dr)
    fg = np.where(g > 0, 1 + 2*ALPHA*np.log(g), 0.0)
    V = g**3/3 - g**4/4
    return 4*np.pi * np.trapezoid(r**2 * (0.5*fg*gp**2 + V), r)


# ================================================================
#  STEP 1: Fine phase map
# ================================================================

print("\n" + "=" * 70)
print("  STEP 1: FINE PHASE MAP")
print("=" * 70)

g0_range = np.concatenate([
    np.linspace(0.85, 0.95, 40),   # electron
    np.linspace(0.95, 1.15, 20),   # near g*=1
    np.linspace(1.15, 1.55, 40),   # muon
    np.linspace(1.55, 2.2, 60),    # tau region (fine)
    np.linspace(2.2, 3.5, 30),     # beyond
])
g0_range = np.unique(np.sort(g0_range))

data = []
for g0 in g0_range:
    r, g = solve_soliton(g0)
    if r is None or r[-1] < 250:
        continue
    B, C, Atail, delta, rms = extract_tail(r, g)
    if np.isnan(Atail) or rms > 0.1*Atail:
        continue
    Sk = compute_Skin(r, g)
    Et = compute_Etot(r, g)
    data.append({
        'g0': g0, 'B': B, 'C': C, 'A': Atail,
        'delta': delta, 'rms': rms, 'Skin': Sk, 'Etot': Et
    })

print("  Computed %d valid solutions" % len(data))

g0_a = np.array([d['g0'] for d in data])
B_a = np.array([d['B'] for d in data])
C_a = np.array([d['C'] for d in data])
A_a = np.array([d['A'] for d in data])
delta_a = np.array([d['delta'] for d in data])
Skin_a = np.array([d['Skin'] for d in data])
Etot_a = np.array([d['Etot'] for d in data])

# Key reference values
r_e, g_e = solve_soliton(G0_E)
B_e, C_e, A_e, delta_e, _ = extract_tail(r_e, g_e)
Sk_e = compute_Skin(r_e, g_e)
Et_e = compute_Etot(r_e, g_e)

r_mu, g_mu = solve_soliton(G0_MU)
B_mu, C_mu, A_mu, delta_mu, _ = extract_tail(r_mu, g_mu)
Sk_mu = compute_Skin(r_mu, g_mu)
Et_mu = compute_Etot(r_mu, g_mu)

print("\n  ELECTRON (g0=%.6f):" % G0_E)
print("    A=%.8f, B=%.8f, C=%.8f" % (A_e, B_e, C_e))
print("    delta=%.4f rad = %.2f deg" % (delta_e, np.degrees(delta_e)))
print("    S_kin=%.6f, E_tot=%.6f" % (Sk_e, Et_e))

print("\n  MUON (g0=%.6f):" % G0_MU)
print("    A=%.8f, B=%.8f, C=%.8f" % (A_mu, B_mu, C_mu))
print("    delta=%.4f rad = %.2f deg" % (delta_mu, np.degrees(delta_mu)))
print("    S_kin=%.6f, E_tot=%.6f" % (Sk_mu, Et_mu))

dp = delta_mu - delta_e
print("\n  Delta(e->mu) = %.4f rad = %.2f deg" % (dp, np.degrees(dp)))
print("  2*pi/3 = %.4f rad = 120.00 deg" % (2*np.pi/3))
print("  Deviation from 2pi/3: %.4f rad = %.2f deg" % (dp-2*np.pi/3, np.degrees(dp-2*np.pi/3)))


# ================================================================
#  STEP 2: S_kin as mass proxy
# ================================================================

print("\n" + "=" * 70)
print("  STEP 2: S_kin AS MASS PROXY")
print("=" * 70)

# Hypothesis: m ~ S_kin^p for some power p
# Check: r_21 = S_kin(mu) / S_kin(e) should give 206.768
r21_Skin = Sk_mu / Sk_e
print("  S_kin(e)  = %.6f" % Sk_e)
print("  S_kin(mu) = %.6f" % Sk_mu)
print("  r_21(S_kin) = Sk_mu/Sk_e = %.4f" % r21_Skin)
print("  r_21(PDG) = %.3f" % R21_PDG)

# If m ~ Sk^p, then r_21 = (Sk_mu/Sk_e)^p = 206.768
# p = ln(206.768) / ln(Sk_mu/Sk_e)
if r21_Skin > 1:
    p_Skin = np.log(R21_PDG) / np.log(r21_Skin)
    print("  Power p such that (Sk_mu/Sk_e)^p = r_21: p = %.4f" % p_Skin)
else:
    p_Skin = None
    print("  ERROR: Sk_mu/Sk_e <= 1, cannot use as mass proxy")

# Same for A_tail^4
r21_A4 = (A_mu/A_e)**4
print("\n  For comparison:")
print("  r_21(A^4) = (A_mu/A_e)^4 = %.4f (PDG: %.3f)" % (r21_A4, R21_PDG))

# What about hybrid: m ~ A^4 * Sk^q?
# r_21 = (A_mu/A_e)^4 * (Sk_mu/Sk_e)^q = 206.768
# (Sk_mu/Sk_e)^q = 206.768 / r21_A4
if r21_A4 > 0 and r21_Skin > 1:
    hybrid_ratio = R21_PDG / r21_A4
    q_hybrid = np.log(hybrid_ratio) / np.log(r21_Skin)
    print("\n  Hybrid: m ~ A^4 * Sk^q")
    print("  q = %.6f" % q_hybrid)
    print("  (q=0 means A^4 alone suffices for e/mu ratio)")
    print("  (|q| << 1 confirms A^4 is dominant for e/mu)")

# Now predict tau with S_kin alone
# Find the saddle/maximum of S_kin
print("\n  S_kin landscape near tau region:")
print("  g0       S_kin      A_tail     delta(deg)  r31(A^4)  r31(Sk^p)")
for d in data:
    if 1.6 < d['g0'] < 2.3:
        r31_A = (d['A']/A_e)**4
        r31_Sk = (d['Skin']/Sk_e)**p_Skin if p_Skin else 0
        print("  %.4f   %10.4f %10.6f  %8.2f  %8.1f  %8.1f" %
              (d['g0'], d['Skin'], d['A'], np.degrees(d['delta']), r31_A, r31_Sk))


# ================================================================
#  STEP 3: E_tot as mass proxy
# ================================================================

print("\n" + "=" * 70)
print("  STEP 3: E_tot (TOTAL ENERGY) AS MASS PROXY")
print("=" * 70)

# Most natural: physical mass = total soliton energy
r21_Etot = Et_mu / Et_e
print("  E_tot(e)  = %.6f" % Et_e)
print("  E_tot(mu) = %.6f" % Et_mu)
print("  r_21(Etot) = Et_mu/Et_e = %.6f" % r21_Etot)

if abs(Et_e) > 1e-15 and r21_Etot > 0:
    p_Etot = np.log(R21_PDG) / np.log(abs(r21_Etot)) if abs(r21_Etot) != 1 else None
    if p_Etot:
        print("  Power p: (Et_mu/Et_e)^p = r_21: p = %.4f" % p_Etot)
else:
    p_Etot = None
    print("  E_tot ratio issue: r21_Etot = %.6f" % r21_Etot)

# Landscape
print("\n  E_tot landscape near tau:")
print("  g0       E_tot      r31(Etot^p)")
for d in data:
    if 1.6 < d['g0'] < 2.3:
        if p_Etot and Et_e != 0:
            r31_Et = abs(d['Etot']/Et_e)**p_Etot
        else:
            r31_Et = 0
        print("  %.4f   %12.6f  %10.1f" % (d['g0'], d['Etot'], r31_Et))


# ================================================================
#  STEP 4: Phase-corrected mass formula
# ================================================================

print("\n" + "=" * 70)
print("  STEP 4: PHASE-CORRECTED MASS FORMULA  m ~ A^4 * F(delta)")
print("=" * 70)

# The key question: what F(delta) gives:
#   F(delta_e) such that r_21 is preserved
#   F(delta_tau) / F(delta_e) = 4.18

# Since A^4 already gives r_21 correctly (by construction via phi-FP),
# F must satisfy F(delta_mu)/F(delta_e) = 1.0 (or close)
# while F(delta_tau)/F(delta_e) = 3477 / (A_tau/A_e)^4

# Let's parameterize delta_tau and see what tau requires
# Candidate: delta_tau = delta_e + 2*(delta_mu - delta_e) = delta_e + 2*Dp
# (equal spacing: delta_e, delta_mu = delta_e + Dp, delta_tau = delta_e + 2*Dp)

delta_tau_eq = delta_e + 2*dp
print("  Equal spacing: delta_tau = %.2f deg" % np.degrees(delta_tau_eq))

# Candidate: delta_tau = delta_e + 2*pi*2/3  (two thirds from electron)
delta_tau_2pi3 = delta_e + 2*(2*np.pi/3)
# Wrap to [-pi, pi]
delta_tau_2pi3 = (delta_tau_2pi3 + np.pi) % (2*np.pi) - np.pi
print("  Two-thirds: delta_tau = delta_e + 4pi/3 = %.2f deg" % np.degrees(delta_tau_2pi3))

# Candidate: exact equal spacing at 120 deg
deltas_120 = [delta_e, delta_e + 2*np.pi/3, delta_e + 4*np.pi/3]
deltas_120 = [(d + np.pi) % (2*np.pi) - np.pi for d in deltas_120]
print("\n  Exact 120-degree spacing:")
for i, (name, d) in enumerate(zip(['e', 'mu', 'tau'], deltas_120)):
    print("    delta_%s = %.2f deg" % (name, np.degrees(d)))
print("    delta_mu(actual) = %.2f deg" % np.degrees(delta_mu))
print("    delta_mu(120-model) = %.2f deg, diff = %.2f deg" %
      (np.degrees(deltas_120[1]), np.degrees(delta_mu - deltas_120[1])))

# Find the soliton at delta_tau candidate from g0 map
from scipy.interpolate import interp1d
delta_uw = np.unwrap(delta_a)
g0_of_delta = interp1d(delta_uw, g0_a, kind='cubic', fill_value='extrapolate')
A_of_g0 = interp1d(g0_a, A_a, kind='cubic', fill_value='extrapolate')
Sk_of_g0 = interp1d(g0_a, Skin_a, kind='cubic', fill_value='extrapolate')

# Scan tau candidates
print("\n  Tau candidates from different phase models:")
print("  %-20s  delta_tau(deg)  g0_tau    A_tau      r31(A^4)   F_needed" % "Model")
print("  " + "-" * 90)

candidates = [
    ("Equal spacing", delta_tau_eq),
    ("120-deg exact", deltas_120[2]),
    ("C=0 (pure cos)", 0.0),
    ("C=0 (pure cos,-)", np.pi),
]
# Also add: where equal A^4 * F gives r_31 = 3477
# for different F forms

for name, dt in candidates:
    dt_uw = dt
    # Find in unwrapped delta
    try:
        g0_t = g0_of_delta(dt_uw)
        if 0.85 < g0_t < 3.5:
            A_t = float(A_of_g0(g0_t))
            r31_bare = (A_t/A_e)**4
            F_needed = R31_PDG / r31_bare if r31_bare > 0 else float('inf')
            print("  %-20s  %10.2f     %.4f    %.6f   %8.1f    %.3f" %
                  (name, np.degrees(dt), g0_t, A_t, r31_bare, F_needed))
        else:
            print("  %-20s  %10.2f     (out of range: g0=%.3f)" % (name, np.degrees(dt), g0_t))
    except Exception as ex:
        print("  %-20s  %10.2f     (interpolation failed: %s)" % (name, np.degrees(dt), str(ex)[:40]))


# ================================================================
#  STEP 5: Koide from (A, delta) geometry
# ================================================================

print("\n" + "=" * 70)
print("  STEP 5: KOIDE CONDITION IN (A, delta) SPACE")
print("=" * 70)

# Koide: Q = (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3
# If m_i ~ A_i^4 * F(delta_i), then
# Q = sum(A_i^4 * Fi) / (sum(A_i^2 * sqrt(Fi)))^2

# Simplify: define a_i = A_i^2 (so m_i ~ a_i^2 * Fi)
# Then Q = sum(a_i^2 * Fi) / (sum(a_i * sqrt(Fi)))^2

# For F = const (no phase correction), Koide becomes:
# Q0 = sum(A_i^4) / (sum(A_i^2))^2
# Let's compute this for various tau candidates

print("\n  Koide Q for A^4 mass formula (F=1):")
print("  g0_tau     A_tau      Q(A^4)     Q-2/3")
for d in data:
    if 1.5 < d['g0'] < 2.5:
        a_e = A_e**2
        a_mu = A_mu**2
        a_tau = d['A']**2
        Q = (a_e**2 + a_mu**2 + a_tau**2) / (a_e + a_mu + a_tau)**2
        print("  %.4f   %.6f   %.6f   %+.6f" % (d['g0'], d['A'], Q, Q - 2.0/3))

# Find g0 where Q(A^4) = 2/3 exactly
print("\n  Searching for Q = 2/3 in A^4 model:")
def Q_of_g0(g0_tau):
    """Koide Q for given g0_tau, using A^4 mass formula."""
    A_t = float(A_of_g0(g0_tau))
    a_e, a_mu, a_tau = A_e**2, A_mu**2, A_t**2
    return (a_e**2 + a_mu**2 + a_tau**2) / (a_e + a_mu + a_tau)**2

# Scan
g0_scan = np.linspace(1.5, 3.0, 200)
Q_scan = []
for g0 in g0_scan:
    try:
        Q_scan.append(Q_of_g0(g0))
    except:
        Q_scan.append(np.nan)
Q_scan = np.array(Q_scan)

# Find where Q = 2/3
for i in range(len(Q_scan)-1):
    if not np.isnan(Q_scan[i]) and not np.isnan(Q_scan[i+1]):
        if (Q_scan[i] - 2.0/3) * (Q_scan[i+1] - 2.0/3) < 0:
            g0_koide = g0_scan[i] + (2.0/3 - Q_scan[i])*(g0_scan[i+1]-g0_scan[i])/(Q_scan[i+1]-Q_scan[i])
            r, g = solve_soliton(g0_koide)
            if r is not None and r[-1] >= 250:
                B_k, C_k, A_k, d_k, _ = extract_tail(r, g)
                r31_k = (A_k/A_e)**4
                Sk_k = compute_Skin(r, g)
                print("  Q=2/3 at g0 = %.6f" % g0_koide)
                print("    A_tail = %.6f, delta = %.2f deg" % (A_k, np.degrees(d_k)))
                print("    r_31(A^4) = %.1f (PDG: %.1f)" % (r31_k, R31_PDG))
                print("    S_kin = %.4f" % Sk_k)
                print("    B=%.6f, C=%.6f" % (B_k, C_k))

                # Check Koide with S_kin
                if p_Skin:
                    r31_Sk = (Sk_k / Sk_e)**p_Skin
                    print("    r_31(Sk^p) = %.1f" % r31_Sk)


# ================================================================
#  STEP 6: Three-generation Koide with phase correction
# ================================================================

print("\n" + "=" * 70)
print("  STEP 6: THREE-GENERATION KOIDE WITH PHASE CORRECTION")
print("=" * 70)

# Key idea: if m_i = A_i^4 * F(delta_i), and Koide Q = 2/3,
# find F(delta) that simultaneously:
#   (a) preserves r_21 = 206.768  => F(delta_mu)/F(delta_e) = r_21 / (A_mu/A_e)^4
#   (b) gives Q = 2/3
#   (c) gives r_31 = 3477.48

# From (a): F ratio for mu/e
F_ratio_mu_e = R21_PDG / (A_mu/A_e)**4
print("  F(delta_mu)/F(delta_e) = %.8f" % F_ratio_mu_e)
print("  (should be ~1 if A^4 already gives r_21)")
print()

# Now: scan possible g0_tau values, compute F_needed for each
# to give r_31, and check if Koide Q = 2/3

print("  Scan: for each g0_tau, find F_tau/F_e giving r_31, check Koide")
print("  g0_tau   A_tau    delta_tau  F_tau/F_e  r31(corrected)  Q(corrected)")
print("  " + "-" * 80)

best_Q_diff = 1e10
best_result = None

for d in data:
    if d['g0'] < 1.5 or d['g0'] > 2.5:
        continue
    A_t = d['A']
    delta_t = d['delta']

    # F ratios:
    # F_e = 1 (reference)
    # F_mu = F_ratio_mu_e = r_21 / (A_mu/A_e)^4
    # F_tau = r_31 / (A_t/A_e)^4
    F_e = 1.0
    F_mu = F_ratio_mu_e
    F_tau = R31_PDG / (A_t/A_e)**4

    # Masses: m_i = A_i^4 * F_i * m_e_ref
    # Ratios: m_e = A_e^4 * F_e, m_mu = A_mu^4 * F_mu, m_tau = A_t^4 * F_tau
    m_e_n = A_e**4 * F_e
    m_mu_n = A_mu**4 * F_mu
    m_tau_n = A_t**4 * F_tau

    # Koide:
    Q = (m_e_n + m_mu_n + m_tau_n) / (np.sqrt(m_e_n) + np.sqrt(m_mu_n) + np.sqrt(m_tau_n))**2

    if abs(Q - 2.0/3) < abs(best_Q_diff):
        best_Q_diff = Q - 2.0/3
        best_result = {
            'g0': d['g0'], 'A': A_t, 'delta': delta_t,
            'F_tau': F_tau, 'Q': Q, 'Skin': d['Skin']
        }

    if abs(d['g0'] % 0.05) < 0.02 or abs(Q - 2.0/3) < 0.005:
        print("  %.4f   %.6f  %7.2f     %.4f       %.1f         %.6f %s" %
              (d['g0'], A_t, np.degrees(delta_t), F_tau,
               R31_PDG, Q,
               " <--" if abs(Q - 2.0/3) < 0.005 else ""))

if best_result:
    print("\n  BEST Koide match:")
    print("    g0_tau = %.6f" % best_result['g0'])
    print("    A_tail = %.6f" % best_result['A'])
    print("    delta  = %.2f deg" % np.degrees(best_result['delta']))
    print("    F_tau/F_e = %.4f (phase enhancement)" % best_result['F_tau'])
    print("    Q = %.8f (target: %.8f, diff: %+.6f)" %
          (best_result['Q'], 2.0/3, best_result['Q'] - 2.0/3))


# ================================================================
#  STEP 7: F(delta) functional form
# ================================================================

print("\n" + "=" * 70)
print("  STEP 7: FUNCTIONAL FORM OF F(delta)")
print("=" * 70)

# We know:
# F(delta_e) = 1
# F(delta_mu) = F_ratio_mu_e
# If we found a best g0_tau, F(delta_tau) = F_tau
# Can we fit F(delta) = a + b*cos(delta - delta_0) + c*cos(2*(delta-delta_0))?

# Three points determine a+b*cos+c*cos(2*):
if best_result:
    d_e, d_mu, d_tau = delta_e, delta_mu, best_result['delta']
    F_vals = [1.0, F_ratio_mu_e, best_result['F_tau']]

    # F(d) = a + b*cos(d) + c*sin(d)  (3 params, 3 constraints)
    M = np.array([
        [1, np.cos(d_e), np.sin(d_e)],
        [1, np.cos(d_mu), np.sin(d_mu)],
        [1, np.cos(d_tau), np.sin(d_tau)]
    ])
    try:
        abc = np.linalg.solve(M, F_vals)
        a, b, c = abc
        print("  F(delta) = %.4f + %.4f*cos(delta) + %.4f*sin(delta)" % (a, b, c))

        # Check
        for name, dval, Fexp in [("e", d_e, 1.0), ("mu", d_mu, F_ratio_mu_e), ("tau", d_tau, best_result['F_tau'])]:
            Fpred = a + b*np.cos(dval) + c*np.sin(dval)
            print("    F(%s) = %.6f (expected %.6f)" % (name, Fpred, Fexp))

        # Is there a nice form?
        amp = np.sqrt(b**2 + c**2)
        phase = np.arctan2(c, b)
        print("\n  Rewritten: F(delta) = %.4f + %.4f * cos(delta - %.2f deg)" %
              (a, amp, np.degrees(phase)))
        print("  Ratio amp/offset = %.4f" % (amp/a if abs(a) > 1e-10 else float('inf')))
    except np.linalg.LinAlgError:
        print("  Singular system, cannot fit F(delta)")


# ================================================================
#  STEP 8: Alternative — Koide directly from S_kin
# ================================================================

print("\n" + "=" * 70)
print("  STEP 8: KOIDE DIRECTLY FROM S_kin")
print("=" * 70)

# If m = S_kin^p, then Koide:
# Q = (Sk_e^p + Sk_mu^p + Sk_tau^p) / (Sk_e^{p/2} + Sk_mu^{p/2} + Sk_tau^{p/2})^2

if p_Skin:
    print("  Using p = %.4f (from r_21 calibration)" % p_Skin)

    print("\n  g0_tau   S_kin     r31(Sk^p)   Q(Sk)")
    best_Q_Sk = 1e10
    best_Sk = None

    for d in data:
        if d['g0'] < 1.5 or d['g0'] > 2.5:
            continue
        m_e_s = Sk_e**p_Skin
        m_mu_s = Sk_mu**p_Skin
        m_tau_s = d['Skin']**p_Skin if d['Skin'] > 0 else 0

        if m_e_s > 0 and m_mu_s > 0 and m_tau_s > 0:
            Q_sk = (m_e_s + m_mu_s + m_tau_s) / (np.sqrt(m_e_s)+np.sqrt(m_mu_s)+np.sqrt(m_tau_s))**2
            r31_sk = (d['Skin']/Sk_e)**p_Skin

            if abs(Q_sk - 2.0/3) < abs(best_Q_Sk):
                best_Q_Sk = Q_sk - 2.0/3
                best_Sk = {'g0': d['g0'], 'Skin': d['Skin'], 'r31': r31_sk, 'Q': Q_sk, 'delta': d['delta']}

            if abs(d['g0'] % 0.1) < 0.02 or abs(Q_sk - 2.0/3) < 0.01:
                print("  %.4f   %8.4f   %8.1f   %.6f %s" %
                      (d['g0'], d['Skin'], r31_sk, Q_sk,
                       " <--" if abs(Q_sk - 2.0/3) < 0.01 else ""))

    if best_Sk:
        print("\n  BEST S_kin Koide match:")
        print("    g0_tau = %.6f" % best_Sk['g0'])
        print("    S_kin  = %.4f" % best_Sk['Skin'])
        print("    r_31   = %.1f (PDG: %.1f)" % (best_Sk['r31'], R31_PDG))
        print("    delta  = %.2f deg" % np.degrees(best_Sk['delta']))
        print("    Q      = %.8f (target: %.8f, diff: %+.6f)" %
              (best_Sk['Q'], 2.0/3, best_Sk['Q']-2.0/3))

        # Can S_kin get closer to r_31 = 3477 AND Q = 2/3?
        # Find g0 where Q_Sk = 2/3
        print("\n  Searching for Q(Sk) = 2/3 exactly:")
        g0_fine = np.linspace(1.5, 2.5, 500)
        Q_fine = []
        for g0 in g0_fine:
            try:
                Sk = float(Sk_of_g0(g0))
                if Sk > 0:
                    m_e_s = Sk_e**p_Skin
                    m_mu_s = Sk_mu**p_Skin
                    m_tau_s = Sk**p_Skin
                    Q_fine.append((m_e_s + m_mu_s + m_tau_s) / (np.sqrt(m_e_s)+np.sqrt(m_mu_s)+np.sqrt(m_tau_s))**2)
                else:
                    Q_fine.append(np.nan)
            except:
                Q_fine.append(np.nan)
        Q_fine = np.array(Q_fine)

        for i in range(len(Q_fine)-1):
            if not np.isnan(Q_fine[i]) and not np.isnan(Q_fine[i+1]):
                if (Q_fine[i]-2.0/3)*(Q_fine[i+1]-2.0/3) < 0:
                    g0_k = g0_fine[i] + (2.0/3 - Q_fine[i])*(g0_fine[i+1]-g0_fine[i])/(Q_fine[i+1]-Q_fine[i])
                    Sk_k = float(Sk_of_g0(g0_k))
                    r31_k = (Sk_k/Sk_e)**p_Skin
                    # get delta
                    r_t, g_t = solve_soliton(g0_k)
                    if r_t is not None and r_t[-1] >= 250:
                        _, _, _, d_k, _ = extract_tail(r_t, g_t)
                        print("    Q(Sk)=2/3 at g0=%.6f, r31=%.1f, delta=%.2f deg" %
                              (g0_k, r31_k, np.degrees(d_k)))


# ================================================================
#  STEP 9: Combined mass proxy — virial-weighted
# ================================================================

print("\n" + "=" * 70)
print("  STEP 9: VIRIAL AND COMBINED PROXIES")
print("=" * 70)

# The virial theorem for soliton ODE:
# 2*S_kin = -d*E_pot + surface terms
# For d=3: 2*S_kin + 3*V_pot = surface term
# The ratio S_kin / E_tot is a diagnostic of soliton structure

print("  Virial ratio S_kin / E_tot:")
for d in data[::10]:
    if abs(d['Etot']) > 1e-10:
        vr = d['Skin'] / d['Etot']
        print("  g0=%.4f: S_kin/E_tot = %.4f" % (d['g0'], vr))

# Test: m ~ A^2 * S_kin^(1/2) (geometric mean type)
print("\n  Geometric mean: m ~ A^2 * Sk^{1/2}")
r21_geom = (A_mu/A_e)**2 * np.sqrt(Sk_mu/Sk_e)
print("  r_21(geom) = %.4f (PDG: %.3f)" % (r21_geom, R21_PDG))

# Test: m ~ (A^2 * B + A^2 * |C|) type — using B and C directly
# Phase-weighted: m ~ (B^2 + C^2)^2 * |B/(B^2+C^2)| = A^4 * |cos(delta)|?
print("\n  Phase-weighted: m ~ A^4 * |cos(delta)|:")
r21_cosD = (A_mu/A_e)**4 * abs(np.cos(delta_mu)/np.cos(delta_e))
print("  r_21(A^4*|cos d|) = %.4f" % r21_cosD)

print("\n  Phase-weighted: m ~ A^4 * (1 + cos(delta)):")
Fw_e = 1 + np.cos(delta_e)
Fw_mu = 1 + np.cos(delta_mu)
r21_1c = (A_mu/A_e)**4 * Fw_mu / Fw_e
print("  r_21(A^4*(1+cos d)) = %.4f" % r21_1c)

print("\n  Phase-weighted: m ~ A^4 * sin(delta + pi/3)^2:")
Fsin_e = np.sin(delta_e + np.pi/3)**2
Fsin_mu = np.sin(delta_mu + np.pi/3)**2
r21_sin = (A_mu/A_e)**4 * Fsin_mu / Fsin_e
print("  r_21(A^4*sin^2(d+pi/3)) = %.4f" % r21_sin)


# ================================================================
#  STEP 10: SUMMARY AND CONCLUSIONS
# ================================================================

print("\n" + "=" * 70)
print("  SUMMARY AND CONCLUSIONS")
print("=" * 70)

print("""
  KEY NUMERICAL RESULTS:
  =====================

  Phase map:
    delta(e)  = %.2f deg
    delta(mu) = %.2f deg
    Delta(e->mu) = %.2f deg (2pi/3 = 120.00 deg, dev = %.2f deg)

  Mass proxies r_21:
    A_tail^4:           %.4f (PDG: %.3f)
    S_kin^p:            %.4f (p=%.3f)
""" % (
    np.degrees(delta_e), np.degrees(delta_mu), np.degrees(dp),
    np.degrees(dp) - 120.0,
    r21_A4, R21_PDG,
    R21_PDG, p_Skin if p_Skin else 0  # S_kin calibrated by definition
))

if best_result:
    print("  Phase-corrected Koide (m ~ A^4 * F):")
    print("    Best g0_tau = %.4f" % best_result['g0'])
    print("    F_tau/F_e = %.4f" % best_result['F_tau'])
    print("    Q = %.6f (diff from 2/3: %+.6f)" % (best_result['Q'], best_result['Q'] - 2.0/3))
    print()

if best_Sk:
    print("  S_kin Koide (m ~ Sk^p):")
    print("    Best g0_tau = %.4f" % best_Sk['g0'])
    print("    r_31 = %.1f (PDG: %.1f)" % (best_Sk['r31'], R31_PDG))
    print("    Q = %.6f (diff from 2/3: %+.6f)" % (best_Sk['Q'], best_Sk['Q'] - 2.0/3))
    print()

print("  INTERPRETATION:")
print("  ===============")
print("  1. A^4 alone gives r_21 correctly (phi-FP mechanism works)")
print("  2. A^4 alone CANNOT give r_31 (max ~832, need 3477)")
print("  3. Phase difference e->mu = %.1f deg is close to 2pi/3 (120 deg)" % np.degrees(dp))
print("  4. This suggests generations are related by ~120 deg phase rotation")
print("  5. A phase correction F(delta) can close the tau gap")
print("  6. The physical origin of F(delta) likely connects to:")
print("     - ERG running of K(psi) in the core region")
print("     - Virial structure of the soliton")
print("     - Coupling between core topology and tail phase")
