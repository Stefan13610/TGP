#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p127_phase_selection_analysis.py -- Phase selection hypothesis for leptons
==========================================================================

Hipoteza robocza: selekcja generacji leptonów sterowana jest nie tylko
przez A_tail, ale przez fazę ogona (B_tail, C_tail).

Ten skrypt:
  1. Buduje pełną mapę (g0) -> (A_tail, B_tail, C_tail, delta)
  2. Sprawdza warunki fazowe: mu (B=C), tau (C~0), e (B~0?)
  3. Bada stabilność warunków fazowych
  4. Testuje poprawkę fazową do masy: m ~ A^4 * F(B,C)
  5. Szuka warunku Koidego w przestrzeni (A, delta)

Kluczowa uwaga: tau to studnia/siodło, nie zero.

Author: TGP project, session v42+
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

PHI = (1 + np.sqrt(5)) / 2
ALPHA = 2
GG = np.exp(-1/(2*ALPHA))
R21_PDG = 206.768
R31_PDG = 3477.48

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
    # Residual (fit quality)
    fit = A @ coeff
    rms = np.sqrt(np.mean((tl - fit)**2))
    return B, C, A_tail, delta, rms


def extract_tail_multiwindow(r, g):
    """Extract tail params at multiple windows for stability check."""
    windows = [(80,180), (100,200), (120,260), (140,260), (100,260)]
    results = []
    for rmin, rmax in windows:
        B, C, Atail, delta, rms = extract_tail(r, g, rmin, rmax)
        results.append({
            'window': (rmin, rmax), 'B': B, 'C': C,
            'A': Atail, 'delta': delta, 'rms': rms
        })
    return results


# ================================================================
#  STEP 1: Complete phase map
# ================================================================

print("=" * 60)
print("  PHASE MAP OF SOLITON TAIL: (g0) -> (B, C, A, delta)")
print("=" * 60)

g0_range = np.concatenate([
    np.linspace(0.80, 1.0, 50),   # electron region
    np.linspace(1.0, 1.6, 50),    # muon region
    np.linspace(1.6, 2.5, 50),    # tau region
    np.linspace(2.5, 4.0, 30),    # beyond tau
])
g0_range = np.unique(np.sort(g0_range))

phase_data = []
for g0 in g0_range:
    r, g = solve_soliton(g0)
    if r is None or r[-1] < 250:
        continue
    B, C, Atail, delta, rms = extract_tail(r, g)
    if np.isnan(Atail):
        continue
    phase_data.append({
        'g0': g0, 'B': B, 'C': C, 'A': Atail,
        'delta': delta, 'rms': rms
    })

print("  Computed %d valid soliton solutions" % len(phase_data))
print()

# Print table for key points
print("  g0       B_tail     C_tail     A_tail     delta(rad)  delta(deg)")
print("  " + "-" * 70)
for d in phase_data[::5]:  # every 5th point
    print("  %.4f  %10.6f %10.6f %10.6f %10.4f  %8.2f" %
          (d['g0'], d['B'], d['C'], d['A'], d['delta'], np.degrees(d['delta'])))

# ================================================================
#  STEP 2: Detailed analysis at e, mu, tau
# ================================================================

print("\n" + "=" * 60)
print("  PHASE AT LEPTON FIXED POINTS")
print("=" * 60)

for label, g0 in [("electron", G0_E), ("muon", G0_MU)]:
    r, g = solve_soliton(g0)
    B, C, Atail, delta, rms = extract_tail(r, g)
    print("\n  %s (g0 = %.6f):" % (label, g0))
    print("    B_tail = %.8f" % B)
    print("    C_tail = %.8f" % C)
    print("    A_tail = %.8f" % Atail)
    print("    delta  = %.6f rad = %.2f deg" % (delta, np.degrees(delta)))
    print("    B/C    = %.6f" % (B/C if abs(C) > 1e-15 else float('inf')))
    print("    C/B    = %.6f" % (C/B if abs(B) > 1e-15 else float('inf')))
    print("    rms    = %.2e (fit quality)" % rms)

    # Stability check
    results = extract_tail_multiwindow(r, g)
    deltas = [res['delta'] for res in results]
    print("    Phase stability: delta = %.4f +/- %.4f (over %d windows)" %
          (np.mean(deltas), np.std(deltas), len(deltas)))


# ================================================================
#  STEP 3: Find special phase conditions
# ================================================================

print("\n" + "=" * 60)
print("  SPECIAL PHASE CONDITIONS")
print("=" * 60)

g0_arr = np.array([d['g0'] for d in phase_data])
B_arr = np.array([d['B'] for d in phase_data])
C_arr = np.array([d['C'] for d in phase_data])
A_arr = np.array([d['A'] for d in phase_data])
delta_arr = np.array([d['delta'] for d in phase_data])

# Condition 1: B = 0 (pure sine)
print("\n  --- Condition B_tail = 0 (pure sine) ---")
for i in range(len(g0_arr)-1):
    if B_arr[i] * B_arr[i+1] < 0:
        # Linear interpolation
        g0_cross = g0_arr[i] - B_arr[i] * (g0_arr[i+1]-g0_arr[i]) / (B_arr[i+1]-B_arr[i])
        # Refine with binary search
        r, g = solve_soliton(g0_cross)
        B_c, C_c, A_c, d_c, _ = extract_tail(r, g) if r is not None else (np.nan,)*5
        print("  B=0 crossing at g0 = %.6f: B=%.2e, C=%.6f, A=%.6f, delta=%.2f deg" %
              (g0_cross, B_c, C_c, A_c, np.degrees(d_c)))

# Condition 2: C = 0 (pure cosine)
print("\n  --- Condition C_tail = 0 (pure cosine) ---")
for i in range(len(g0_arr)-1):
    if C_arr[i] * C_arr[i+1] < 0:
        g0_cross = g0_arr[i] - C_arr[i] * (g0_arr[i+1]-g0_arr[i]) / (C_arr[i+1]-C_arr[i])
        r, g = solve_soliton(g0_cross)
        B_c, C_c, A_c, d_c, _ = extract_tail(r, g) if r is not None else (np.nan,)*5
        print("  C=0 crossing at g0 = %.6f: B=%.6f, C=%.2e, A=%.6f, delta=%.2f deg" %
              (g0_cross, B_c, C_c, A_c, np.degrees(d_c)))

# Condition 3: B = C (quadrature, delta=pi/4)
print("\n  --- Condition B_tail = C_tail (quadrature, delta=pi/4) ---")
diff_BC = B_arr - C_arr
for i in range(len(g0_arr)-1):
    if diff_BC[i] * diff_BC[i+1] < 0:
        g0_cross = g0_arr[i] - diff_BC[i] * (g0_arr[i+1]-g0_arr[i]) / (diff_BC[i+1]-diff_BC[i])
        r, g = solve_soliton(g0_cross)
        B_c, C_c, A_c, d_c, _ = extract_tail(r, g) if r is not None else (np.nan,)*5
        print("  B=C crossing at g0 = %.6f: B=%.6f, C=%.6f, A=%.6f, delta=%.2f deg" %
              (g0_cross, B_c, C_c, A_c, np.degrees(d_c)))

# Condition 4: B = -C (delta = -pi/4 or 3pi/4)
print("\n  --- Condition B_tail = -C_tail (delta=-pi/4) ---")
sum_BC = B_arr + C_arr
for i in range(len(g0_arr)-1):
    if sum_BC[i] * sum_BC[i+1] < 0:
        g0_cross = g0_arr[i] - sum_BC[i] * (g0_arr[i+1]-g0_arr[i]) / (sum_BC[i+1]-sum_BC[i])
        r, g = solve_soliton(g0_cross)
        B_c, C_c, A_c, d_c, _ = extract_tail(r, g) if r is not None else (np.nan,)*5
        print("  B=-C crossing at g0 = %.6f: B=%.6f, C=%.6f, A=%.6f, delta=%.2f deg" %
              (g0_cross, B_c, C_c, A_c, np.degrees(d_c)))

# Condition 5: delta = n*pi/3 (Koide-like thirds)
print("\n  --- Condition delta = n*pi/3 (120-degree spacing) ---")
delta_uw = np.unwrap(delta_arr)
for n in range(-3, 4):
    target = n * np.pi / 3
    for i in range(len(g0_arr)-1):
        if (delta_uw[i] - target) * (delta_uw[i+1] - target) < 0:
            g0_cross = g0_arr[i] + (target - delta_uw[i]) * (g0_arr[i+1]-g0_arr[i]) / (delta_uw[i+1]-delta_uw[i])
            if 0.8 < g0_cross < 4.0:
                r, g = solve_soliton(g0_cross)
                if r is not None and r[-1] >= 250:
                    B_c, C_c, A_c, d_c, _ = extract_tail(r, g)
                    print("  delta=%d*pi/3 = %.2f deg at g0 = %.6f: A=%.6f" %
                          (n, np.degrees(target), g0_cross, A_c))


# ================================================================
#  STEP 4: Phase difference between e and mu
# ================================================================

print("\n" + "=" * 60)
print("  PHASE DIFFERENCES BETWEEN GENERATIONS")
print("=" * 60)

r_e, g_e = solve_soliton(G0_E)
B_e, C_e, A_e, delta_e, _ = extract_tail(r_e, g_e)

r_mu, g_mu = solve_soliton(G0_MU)
B_mu, C_mu, A_mu, delta_mu, _ = extract_tail(r_mu, g_mu)

dp = delta_mu - delta_e
print("  delta_e  = %.6f rad = %.2f deg" % (delta_e, np.degrees(delta_e)))
print("  delta_mu = %.6f rad = %.2f deg" % (delta_mu, np.degrees(delta_mu)))
print("  Delta(e->mu) = %.6f rad = %.2f deg" % (dp, np.degrees(dp)))
print()
print("  Comparison with special values:")
for name, val in [("pi/3", np.pi/3), ("pi/2", np.pi/2), ("2pi/3", 2*np.pi/3),
                   ("pi", np.pi), ("3pi/4", 3*np.pi/4), ("5pi/6", 5*np.pi/6)]:
    print("    %6s = %.4f rad: diff = %.4f rad (%.1f deg)" %
          (name, val, abs(dp - val), abs(np.degrees(dp - val))))

# What g0_tau would have delta_tau such that delta_mu - delta_e = delta_tau - delta_mu?
# (equal phase spacing)
delta_tau_target = 2 * delta_mu - delta_e
print("\n  If equal phase spacing:")
print("    delta_tau (target) = %.4f rad = %.2f deg" % (delta_tau_target, np.degrees(delta_tau_target)))

# Find g0 with this delta
delta_interp = interp1d(g0_arr, np.unwrap(delta_arr), kind='cubic', fill_value='extrapolate')
from scipy.optimize import brentq
try:
    g0_tau_phase = brentq(lambda g: delta_interp(g) - delta_tau_target, 1.5, 3.5)
    r_tau, g_tau = solve_soliton(g0_tau_phase)
    B_tau, C_tau, A_tau, delta_tau_actual, _ = extract_tail(r_tau, g_tau)
    r31_phase = (A_tau / A_e)**4

    print("    g0_tau (equal spacing) = %.6f" % g0_tau_phase)
    print("    A_tau = %.6f" % A_tau)
    print("    r_31 = %.1f (Koide: 3477)" % r31_phase)
    print("    delta_tau = %.4f rad = %.2f deg" % (delta_tau_actual, np.degrees(delta_tau_actual)))
except:
    print("    Could not find g0_tau with equal phase spacing")


# ================================================================
#  STEP 5: Saddle/well structure near tau
# ================================================================

print("\n" + "=" * 60)
print("  SADDLE/WELL STRUCTURE NEAR TAU")
print("=" * 60)
print("  (Tau is a saddle or well, not zero!)")

# Study d(A_tail)/d(g0) and d(delta)/d(g0) near the A_max
A_interp = interp1d(g0_arr, A_arr, kind='cubic')

# Find A_tail maximum
from scipy.optimize import minimize_scalar
res = minimize_scalar(lambda g: -A_interp(g), bounds=(1.5, 2.5), method='bounded')
g0_Amax = res.x
A_max = -res.fun

print("  A_tail maximum at g0 = %.6f, A_max = %.6f" % (g0_Amax, A_max))
print("  r_31(A_max) = %.1f" % ((A_max/A_e)**4))
print()

# Check the landscape near A_max
print("  Landscape near A_max:")
for dg in [-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3]:
    g0_test = g0_Amax + dg
    if 0.8 < g0_test < 4.0:
        r_t, g_t = solve_soliton(g0_test)
        if r_t is not None and r_t[-1] >= 250:
            B_t, C_t, A_t, d_t, _ = extract_tail(r_t, g_t)
            r31_t = (A_t/A_e)**4
            print("    g0=%.4f: A=%.6f, r31=%.1f, B=%.6f, C=%.6f, delta=%.2f deg" %
                  (g0_test, A_t, r31_t, B_t, C_t, np.degrees(d_t)))

# Is there a SADDLE in the (B,C) plane?
# Plot B(g0) and C(g0) to see the structure
print("\n  B(g0) and C(g0) behavior near saddle/maximum:")
g0_fine = np.linspace(1.5, 2.5, 100)
B_fine = []
C_fine = []
for g0 in g0_fine:
    r, g = solve_soliton(g0)
    if r is None or r[-1] < 250:
        B_fine.append(np.nan)
        C_fine.append(np.nan)
        continue
    B, C, _, _, _ = extract_tail(r, g)
    B_fine.append(B)
    C_fine.append(C)
B_fine = np.array(B_fine)
C_fine = np.array(C_fine)

# Find where dB/dg0 = 0 and dC/dg0 = 0 (critical points)
valid = ~np.isnan(B_fine) & ~np.isnan(C_fine)
g0_v = g0_fine[valid]
B_v = B_fine[valid]
C_v = C_fine[valid]

dB = np.gradient(B_v, g0_v)
dC = np.gradient(C_v, g0_v)

# Find B extrema
for i in range(len(g0_v)-1):
    if dB[i]*dB[i+1] < 0:
        g0_ext = g0_v[i] - dB[i]*(g0_v[i+1]-g0_v[i])/(dB[i+1]-dB[i])
        print("  B extremum at g0 ~ %.4f" % g0_ext)
for i in range(len(g0_v)-1):
    if dC[i]*dC[i+1] < 0:
        g0_ext = g0_v[i] - dC[i]*(g0_v[i+1]-g0_v[i])/(dC[i+1]-dC[i])
        print("  C extremum at g0 ~ %.4f" % g0_ext)


# ================================================================
#  STEP 6: Phase-corrected mass formula
# ================================================================

print("\n" + "=" * 60)
print("  PHASE-CORRECTED MASS FORMULA")
print("=" * 60)

# Hypothesis: m ~ A^4 * F(delta)
# If F(delta) can enhance the tau mass beyond A_max^4...
# Test F(delta) = 1 + k * cos(delta - delta_0) or similar

# First: what F(delta) would be needed?
# At tau saddle: A_tau ~ A_max = 0.547, but need m_tau/m_e = 3477
# m_tau/m_e = (A_tau/A_e)^4 * F(delta_tau)/F(delta_e)
# 3477 = 832 * F(delta_tau)/F(delta_e)
# F(delta_tau)/F(delta_e) = 3477/832 = 4.18

print("  Enhancement factor needed: F(tau)/F(e) = %.2f" % (R31_PDG / ((A_max/A_e)**4)))
print()

# Could the phase correction come from the INNER region?
# The mass formula m ~ A_tail^4 assumes the tail IS the coupling to external fields.
# But the CORE of the soliton (inner region near g*) also contributes.
# Perhaps: m ~ integral of f(g) * (g')^2 * r^2 dr (kinetic action)
# This integral depends on BOTH the core structure AND the tail.

print("  Testing kinetic action as mass proxy:")
for label, g0 in [("e", G0_E), ("mu", G0_MU), ("saddle", g0_Amax)]:
    r, g = solve_soliton(g0)
    if r is None:
        continue
    dr = r[1] - r[0]
    gp = np.gradient(g, dr)
    fg = np.array([1 + 2*ALPHA*np.log(gi) if gi > 0 else 0 for gi in g])
    # Kinetic action: S_kin = 4*pi * int r^2 * (1/2)*f(g)*(g')^2 dr
    S_kin = 4*np.pi * np.trapezoid(r**2 * 0.5 * fg * gp**2, r)
    # Total energy
    V_vals = np.array([gi**3/3 - gi**4/4 for gi in g])
    E_tot = 4*np.pi * np.trapezoid(r**2 * (0.5*fg*gp**2 + V_vals), r)

    B_t, C_t, A_t, d_t, _ = extract_tail(r, g)
    r_ratio = (A_t/A_e)**4 if label != "e" else 1.0
    print("    %s (g0=%.4f): S_kin=%.4f, E_tot=%.4f, A_tail=%.6f, r_mass=%.1f" %
          (label, g0, S_kin, E_tot, A_t, r_ratio))


# ================================================================
#  STEP 7: Node count and topology
# ================================================================

print("\n" + "=" * 60)
print("  SOLITON NODE COUNT AND TOPOLOGY")
print("=" * 60)

for label, g0 in [("e", G0_E), ("mu", G0_MU), ("saddle", g0_Amax),
                   ("1.2", 1.2), ("1.5", 1.5), ("2.5", 2.5)]:
    r, g = solve_soliton(g0)
    if r is None:
        continue
    # Count zero-crossings of (g-1) in inner region
    inner = r < 50
    gi = g[inner]
    crossings = sum(1 for i in range(len(gi)-1) if (gi[i]-1)*(gi[i+1]-1) < 0)
    # Count local extrema of g
    gp = np.gradient(g, r[1]-r[0])
    extrema_inner = sum(1 for i in range(1, sum(inner)-1) if gp[inner][i-1]*gp[inner][i+1] < 0)

    B_t, C_t, A_t, d_t, _ = extract_tail(r, g)
    print("  %s (g0=%.4f): nodes(r<50)=%d, extrema(r<50)=%d, delta=%.2f deg, A=%.5f" %
          (label, g0, crossings, extrema_inner, np.degrees(d_t), A_t))


# ================================================================
#  SUMMARY
# ================================================================

print("\n" + "=" * 60)
print("  PHASE SELECTION ANALYSIS SUMMARY")
print("=" * 60)
print()
print("  Key results:")
print("    - Phase map computed for %d soliton solutions" % len(phase_data))
print("    - Special phase crossings identified")
print()
print("  Phase at lepton points:")
print("    e:  delta = %.2f deg" % np.degrees(delta_e))
print("    mu: delta = %.2f deg" % np.degrees(delta_mu))
print("    Delta(e->mu) = %.2f deg" % np.degrees(dp))
print()
print("  Saddle structure:")
print("    A_max = %.6f at g0 = %.4f" % (A_max, g0_Amax))
print("    This is a saddle (not zero) in the (B,C) landscape")
print()
print("  Status: analyzing phase selection mechanism")
