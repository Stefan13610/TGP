#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p135_g0tau_constraint.py - What constrains g0_tau?
====================================================
Systematic investigation of candidate mechanisms:

A) Algebraic: g0_tau = phi^n * g0_e  or  g0_tau/g0_mu = phi^k
B) Energy: minimize soliton energy E[g0] = int [f(g)*g'^2/2 + V(g)] r^2 dr
C) Topological: winding number or quantization of action
D) Stability: 3D perturbative stability (negative eigenvalues)
E) Phase: combined phase condition Delta(e->tau) = special angle
F) Potential: V(g0_tau) = special value, or V_eff with running alpha

Key data from p134f:
  g0_e  = 0.905481, g0_mu = phi*g0_e = 1.465099
  At g0_tau=4.0: r_31=3477.8 (matches Koide)
  g0_tau/g0_e = 4.0/0.905481 = 4.417
  g0_tau/g0_mu = 4.0/1.465099 = 2.730

Author: TGP project
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
ETA_K = 12.067
G0_E = 0.905481
G0_MU = PHI * G0_E
R21_PDG = 206.768
R31_PDG = 3477.48

def solve_soliton(g0, eta_K=ETA_K, rm=300):
    def fk(g):
        a = 2.0 / (1 + eta_K * (g - 1)**2)
        return 1 + 2*a*np.log(g) if g > 0 else -1e30
    def Vp(g):
        return g**2*(1-g)
    fg0 = fk(g0)
    if abs(fg0) < 1e-15: return None, None, None
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
    r = np.linspace(rs, min(s.t[-1], rm), 15000)
    g_arr = s.sol(r)[0]
    gp_arr = s.sol(r)[1]
    return r, g_arr, gp_arr

def extract_tail(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan, np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    B, C = coeff
    return np.sqrt(B**2 + C**2), np.arctan2(C, B)

def V_func(g):
    return g**3/3.0 - g**4/4.0

def fk_func(g, eta_K=ETA_K):
    a = 2.0 / (1 + eta_K * (g - 1)**2)
    return 1 + 2*a*np.log(g) if g > 0 else -1e30

# Get reference A_e
r_e, g_e, gp_e = solve_soliton(G0_E)
A_e, d_e = extract_tail(r_e, g_e)

# =====================================================================
print("="*70)
print("  PART A: ALGEBRAIC RELATIONS")
print("="*70)

# Test various algebraic relations for g0_tau that gives r_31=3477
# First: what exact g0_tau gives r_31 = 3477.48?
print("\n  Finding exact g0_tau for r_31 = 3477.48...")

def r31_for_g0tau(g0_t):
    r_t, g_t, _ = solve_soliton(g0_t)
    if r_t is None: return 0.0
    A_t, _ = extract_tail(r_t, g_t)
    if np.isnan(A_t): return 0.0
    return (A_t / A_e)**4

# Bisect between 3.9 and 4.1
g_lo, g_hi = 3.9, 4.1
for i in range(40):
    g_mid = (g_lo + g_hi) / 2.0
    r31 = r31_for_g0tau(g_mid)
    if r31 > R31_PDG:
        g_hi = g_mid
    else:
        g_lo = g_mid
    if abs(g_hi - g_lo) < 1e-8:
        break

g0_tau_exact = (g_lo + g_hi) / 2.0
r31_check = r31_for_g0tau(g0_tau_exact)
print(f"  g0_tau(exact) = {g0_tau_exact:.10f}")
print(f"  r_31 = {r31_check:.2f}")

ratio_te = g0_tau_exact / G0_E
ratio_tm = g0_tau_exact / G0_MU
print(f"\n  g0_tau / g0_e  = {ratio_te:.6f}")
print(f"  g0_tau / g0_mu = {ratio_tm:.6f}")
print(f"  g0_mu / g0_e   = {PHI:.6f} = phi")

# Test algebraic candidates for the ratios
print(f"\n  Candidates for g0_tau/g0_e = {ratio_te:.6f}:")
candidates_te = {
    "phi^3": PHI**3,
    "phi^2 + 1": PHI**2 + 1,
    "2*phi + 1": 2*PHI + 1,
    "e + phi": np.e + PHI,
    "phi^3 + phi^(-3)": PHI**3 + PHI**(-3),
    "5 - 1/phi": 5 - 1/PHI,
    "4 + phi/3": 4 + PHI/3,
    "sqrt(5)*phi": np.sqrt(5)*PHI,
    "3 + phi": 3 + PHI,
    "pi*sqrt(2)": np.pi*np.sqrt(2),
    "4*phi - 2": 4*PHI - 2,
    "2*e": 2*np.e,
    "phi^2*phi": PHI**3,
    "F(8)/F(5)": 21.0/5.0,   # Fibonacci
    "F(9)/F(5)": 34.0/5.0,
    "L(4)": 7.0,   # Lucas
    "L(5)/F(3)": 11.0/2.0,
    "3*phi^2/2": 3*PHI**2/2,
}
results = [(abs(v - ratio_te)/ratio_te*100, k, v) for k, v in candidates_te.items()]
results.sort()
for dev, name, val in results[:10]:
    print(f"    {name:25s} = {val:.6f} (dev {dev:.3f}%)")

print(f"\n  Candidates for g0_tau/g0_mu = {ratio_tm:.6f}:")
candidates_tm = {
    "phi^2": PHI**2,
    "e": np.e,
    "phi + 1": PHI + 1,
    "phi^2 + 1/phi": PHI**2 + 1/PHI,
    "sqrt(2*pi)": np.sqrt(2*np.pi),
    "3 - 1/phi": 3 - 1/PHI,
    "pi - 1/3": np.pi - 1.0/3,
    "8/3": 8.0/3,
    "phi*phi": PHI**2,
    "phi + 1/phi + 1/phi^2": PHI + 1/PHI + 1/PHI**2,
    "F(7)/F(5)": 13.0/5.0,
    "1 + phi": 1 + PHI,
    "5/phi": 5.0/PHI,
    "sqrt(phi)*2": np.sqrt(PHI)*2,
    "phi*sqrt(phi)": PHI*np.sqrt(PHI),
}
results = [(abs(v - ratio_tm)/ratio_tm*100, k, v) for k, v in candidates_tm.items()]
results.sort()
for dev, name, val in results[:10]:
    print(f"    {name:25s} = {val:.6f} (dev {dev:.3f}%)")

# =====================================================================
print("\n" + "="*70)
print("  PART B: SOLITON ENERGY")
print("="*70)

# E[g0] = 4*pi * int_0^inf [f(g)*g'^2/2 + V(g) - V(1)] r^2 dr
# If there's a minimum in E(g0), that determines the preferred soliton.

def soliton_energy(g0, eta_K=ETA_K):
    r, g, gp = solve_soliton(g0, eta_K)
    if r is None: return np.nan, np.nan, np.nan
    # Kinetic energy density
    fk_arr = np.array([fk_func(gi, eta_K) for gi in g])
    T_kin = 0.5 * fk_arr * gp**2
    # Potential energy density (shifted so vacuum = 0)
    V_pot = np.array([V_func(gi) - V_func(1.0) for gi in g])
    # Total integrand * r^2
    integrand_T = T_kin * r**2
    integrand_V = V_pot * r**2
    integrand_tot = (T_kin + V_pot) * r**2
    E_kin = 4*np.pi * trapezoid(integrand_T, r)
    E_pot = 4*np.pi * trapezoid(integrand_V, r)
    E_tot = 4*np.pi * trapezoid(integrand_tot, r)
    return E_kin, E_pot, E_tot

print(f"\n  Soliton energy E(g0) = 4*pi * int [f*g'^2/2 + V(g)-V(1)] r^2 dr")
print(f"\n  {'g0':>6s} {'E_kin':>12s} {'E_pot':>12s} {'E_tot':>12s} {'E_tot/g0^2':>12s}")
print(f"  {'-'*58}")

g0_test = [0.5, 0.7, 0.85, G0_E, 1.2, G0_MU, 2.0, 2.5, 3.0, 3.5, g0_tau_exact, 4.5, 5.0, 6.0]
energies = {}
for g0 in g0_test:
    ek, ep, et = soliton_energy(g0)
    energies[g0] = (ek, ep, et)
    mark = ""
    if abs(g0 - G0_E) < 0.01: mark = " (e)"
    elif abs(g0 - G0_MU) < 0.01: mark = " (mu)"
    elif abs(g0 - g0_tau_exact) < 0.02: mark = " (tau?)"
    if not np.isnan(et):
        print(f"  {g0:6.3f} {ek:12.2f} {ep:12.2f} {et:12.2f} {et/g0**2:12.4f}{mark}")
    else:
        print(f"  {g0:6.3f}  {'NaN':>12s}")

# Is there a minimum in E(g0)?
print(f"\n  Looking for E_tot minimum in [1.5, 8.0]...")
g0_scan = np.linspace(1.5, 8.0, 40)
e_scan = []
for g0 in g0_scan:
    _, _, et = soliton_energy(g0)
    e_scan.append(et)
e_scan = np.array(e_scan)

if not np.all(np.isnan(e_scan)):
    valid = ~np.isnan(e_scan)
    if np.any(valid):
        min_idx = np.nanargmin(e_scan)
        print(f"  E_tot min at g0 = {g0_scan[min_idx]:.3f}, E = {e_scan[min_idx]:.2f}")
        # Is it a true local minimum?
        if 0 < min_idx < len(e_scan)-1:
            if e_scan[min_idx] < e_scan[min_idx-1] and e_scan[min_idx] < e_scan[min_idx+1]:
                print(f"  -> TRUE local minimum!")
            else:
                print(f"  -> Boundary minimum (monotonic)")
        # Check if E is monotonic
        diffs = np.diff(e_scan[valid])
        if np.all(diffs > 0):
            print(f"  E_tot is MONOTONICALLY INCREASING in [1.5, 8.0]")
        elif np.all(diffs < 0):
            print(f"  E_tot is MONOTONICALLY DECREASING in [1.5, 8.0]")
        else:
            # Find sign changes
            sign_changes = np.where(np.diff(np.sign(diffs)))[0]
            print(f"  E_tot has {len(sign_changes)} extrema in [1.5, 8.0]")

# =====================================================================
print("\n" + "="*70)
print("  PART C: SOLITON ACTION (Derrick-like virial)")
print("="*70)

# Derrick's theorem in d=3: for a static soliton to exist,
# the virial identity must hold:
#   E_kin + 3*E_pot = 0  (in d=3, for standard kinetic)
# With modified kinetic f(g), the virial becomes:
#   int f(g)*g'^2 r^2 dr + 3 * int V_shifted r^2 dr = 0
# i.e., E_kin = -3*E_pot for the stable soliton.
#
# Check if there's a special g0 where this virial is satisfied.

print(f"\n  Derrick virial: E_kin + 3*E_pot = 0 for stable soliton in d=3")
print(f"\n  {'g0':>6s} {'E_kin':>10s} {'E_pot':>10s} {'E_kin+3*E_pot':>14s} {'ratio E_k/E_p':>14s}")
print(f"  {'-'*58}")

virial_vals = []
for g0 in np.linspace(1.5, 8.0, 30):
    ek, ep, _ = soliton_energy(g0)
    if not np.isnan(ek) and abs(ep) > 1e-10:
        derrick = ek + 3*ep
        ratio = ek / ep
        virial_vals.append((g0, derrick, ratio))
        mark = ""
        if abs(g0 - g0_tau_exact) < 0.15: mark = " <-- tau"
        if abs(derrick) < abs(ek)*0.05: mark += " VIRIAL!"
        print(f"  {g0:6.3f} {ek:10.2f} {ep:10.2f} {derrick:14.2f} {ratio:14.4f}{mark}")

# Check if Derrick virial crosses zero
if len(virial_vals) > 1:
    prev = virial_vals[0]
    for curr in virial_vals[1:]:
        if prev[1] * curr[1] < 0:
            # Bisect
            g_lo, g_hi = prev[0], curr[0]
            for _ in range(30):
                g_mid = (g_lo + g_hi) / 2
                ek, ep, _ = soliton_energy(g_mid)
                if not np.isnan(ek):
                    d_mid = ek + 3*ep
                    if d_mid * prev[1] < 0:
                        g_hi = g_mid
                    else:
                        g_lo = g_mid
                else:
                    g_lo = g_mid
                if abs(g_hi - g_lo) < 1e-6:
                    break
            g_virial = (g_lo + g_hi) / 2
            ek, ep, et = soliton_energy(g_virial)
            r31 = r31_for_g0tau(g_virial)
            print(f"\n  DERRICK VIRIAL ZERO at g0 = {g_virial:.6f}")
            print(f"  E_kin = {ek:.4f}, E_pot = {ep:.4f}, E_tot = {et:.4f}")
            print(f"  r_31 = {r31:.1f} (Koide: {R31_PDG:.1f})")
            print(f"  m_tau = {0.51099895 * r31:.1f} MeV (obs: 1776.86)")
        prev = curr

# =====================================================================
print("\n" + "="*70)
print("  PART D: STABILITY / NEGATIVE MODES")
print("="*70)

# A soliton with n negative modes of the fluctuation operator
# H = -f(g)*d^2/dr^2 - (2/r)*f(g)*d/dr + V''(g) + f''(g)*g'^2/2
# is unstable if n > 1 (the single negative mode corresponds
# to the translation zero mode in d > 1).
#
# Count "zero crossings" of the soliton profile (g-1) to estimate modes.
# In 1D, the number of nodes = number of bound states below threshold.

print(f"\n  Counting zero crossings of (g-1) for each soliton:")
print(f"  {'g0':>6s} {'crossings':>10s} {'last_r':>10s} {'comment':>20s}")

for g0 in [G0_E, G0_MU, 2.0, 3.0, g0_tau_exact, 5.0, 6.0, 8.0]:
    r, g, _ = solve_soliton(g0)
    if r is None:
        print(f"  {g0:6.3f}  FAILED")
        continue
    # Count zero crossings of (g-1)
    diff = g - 1.0
    crossings = np.sum(np.abs(np.diff(np.sign(diff))) > 0)
    # Last crossing radius
    cross_idx = np.where(np.abs(np.diff(np.sign(diff))) > 0)[0]
    last_r = r[cross_idx[-1]] if len(cross_idx) > 0 else 0.0
    mark = ""
    if abs(g0 - G0_E) < 0.01: mark = "electron"
    elif abs(g0 - G0_MU) < 0.01: mark = "muon"
    elif abs(g0 - g0_tau_exact) < 0.02: mark = "tau(Koide)"
    print(f"  {g0:6.3f} {crossings:10d} {last_r:10.1f} {mark:>20s}")

# =====================================================================
print("\n" + "="*70)
print("  PART E: ACTION QUANTIZATION")
print("="*70)

# The Euclidean action S[g0] might be quantized: S = n * pi or S = n * pi/2
# Check if at g0_tau_exact, S takes a special value.

print(f"\n  Euclidean action S = 4*pi * int [f*g'^2/2 + V(g)-V(1)] r^2 dr")
print(f"  (Same as energy, but looking for quantized values)")
print(f"\n  {'g0':>6s} {'S/pi':>10s} {'S/(pi/2)':>10s} {'near_int?':>10s}")

for g0 in [G0_E, G0_MU, 2.0, 3.0, g0_tau_exact, 5.0, 6.0]:
    _, _, et = soliton_energy(g0)
    if np.isnan(et): continue
    s_over_pi = et / np.pi
    s_over_halfpi = et / (np.pi/2)
    near = abs(s_over_pi - round(s_over_pi)) < 0.1
    mark = ""
    if abs(g0 - G0_E) < 0.01: mark = " (e)"
    elif abs(g0 - G0_MU) < 0.01: mark = " (mu)"
    elif abs(g0 - g0_tau_exact) < 0.02: mark = " (tau)"
    print(f"  {g0:6.3f} {s_over_pi:10.4f} {s_over_halfpi:10.4f} {'YES' if near else 'no':>10s}{mark}")

# =====================================================================
print("\n" + "="*70)
print("  PART F: EFFECTIVE POTENTIAL WITH RUNNING ALPHA")
print("="*70)

# With running alpha, the effective kinetic function is:
# f_eff(g) = 1 + 2 * [alpha/(1+eta*(g-1)^2)] * ln(g)
# The "effective potential" in the action is:
# S = int [f_eff(g)*g'^2/2 + V(g)] r^2 dr
#
# One possibility: the soliton mass M(g0) = A(g0)^4 has a
# minimum at some g0 > 1, which would select the tau.
# (For g0 < 1, i.e. electron/muon, the phi-FP mechanism works.)

# Actually, the mass is proportional to A^4 where A is the tail amplitude.
# We already know A grows monotonically. But does A^4 / E_tot have structure?

print(f"\n  Mass-to-energy ratio: A^4 / E_tot")
print(f"  If M=A^4 is 'mass' and E_tot is 'classical energy',")
print(f"  quantum corrections might select M/E = special value.")
print(f"\n  {'g0':>6s} {'A':>10s} {'A^4':>10s} {'E_tot':>10s} {'A^4/E_tot':>10s}")
print(f"  {'-'*50}")

for g0 in np.linspace(1.5, 8.0, 25):
    r, g, gp = solve_soliton(g0)
    if r is None: continue
    A, _ = extract_tail(r, g)
    _, _, et = soliton_energy(g0)
    if np.isnan(A) or np.isnan(et) or abs(et) < 1e-10: continue
    a4 = A**4
    ratio = a4 / et
    mark = ""
    if abs(g0 - g0_tau_exact) < 0.15: mark = " <-- tau"
    print(f"  {g0:6.3f} {A:10.6f} {a4:10.2f} {et:10.2f} {ratio:10.6f}{mark}")

# =====================================================================
print("\n" + "="*70)
print("  SUMMARY")
print("="*70)

print(f"""
  Exact g0_tau for r_31 = {R31_PDG}: g0_tau = {g0_tau_exact:.8f}
  g0_tau / g0_e  = {ratio_te:.6f}
  g0_tau / g0_mu = {ratio_tm:.6f}

  The mechanism fixing g0_tau will determine whether TGP can
  predict m_tau from first principles or needs it as input.
""")
print("DONE")
