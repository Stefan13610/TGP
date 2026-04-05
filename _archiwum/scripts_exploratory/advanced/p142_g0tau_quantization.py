#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p142_g0tau_quantization.py -- New quantization approaches for g0_tau
=====================================================================

After p141 showed 3D stability doesn't fix g0_tau, we try alternative
quantization conditions that might select g0_tau = 4.

Approaches:
  A) V(g0) = -V(1) condition: the symmetric potential argument
  B) Alpha_eff at core: special value at g0=4
  C) Integer condition on g0 from field-space topology
  D) Semiclassical quantization with proper boundary conditions
  E) Action quantization S_soliton = n*pi*hbar
  F) Ratio condition: g0_tau/g0_e = specific algebraic constant

Author: TGP project
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq

PHI = (1+np.sqrt(5))/2
ETA_K = 181.0/15
G0_E = 0.90548144
G0_MU = PHI * G0_E

def solve_sol(g0, eta=ETA_K, rm=300):
    def fk(g):
        a = 2.0/(1+eta*(g-1)**2)
        return 1+2*a*np.log(g) if g>0 else -1e30
    def Vp(g): return g**2*(1-g)
    fg0 = fk(g0)
    if abs(fg0)<1e-15: return None,None,None
    c2 = Vp(g0)/(3*fg0)
    rs=0.01
    def rhs(r,y):
        g,p=y
        if g<=1e-15: return [p,0]
        fg=fk(g)
        if abs(fg)<1e-10: return [p,0]
        if r<1e-10: return [p,Vp(g)/fg/3]
        return [p,(Vp(g)-2/r*p)/fg]
    def ev(r,y): return 100-abs(y[0])
    ev.terminal=True
    s=solve_ivp(rhs,[rs,rm],[g0+c2*rs**2,2*c2*rs],
                method='RK45',rtol=1e-11,atol=1e-13,
                max_step=0.05,events=[ev],dense_output=True)
    r=np.linspace(rs,min(s.t[-1],rm),15000)
    sol=s.sol(r)
    return r,sol[0],sol[1]

def fk_func(g, eta=ETA_K):
    a = 2.0/(1+eta*(g-1)**2)
    return 1+2*a*np.log(g) if g>0 else -1e30

def alpha_eff(g, eta=ETA_K):
    return 2.0/(1+eta*(g-1)**2)

# ===================================================================
print("="*70)
print("  PART A: POTENTIAL SYMMETRY CONDITION")
print("="*70)

# V(g) = g^3/3 - g^4/4
# V(1) = 1/12 (vacuum energy)
# Question: at what g does V(g) equal special values?

def V(g): return g**3/3 - g**4/4

print(f"\n  V(g) = g^3/3 - g^4/4")
print(f"  V(1) = {V(1):.6f} = 1/12")
print(f"  V(0) = {V(0):.6f}")
print(f"  V(4/3) = {V(4/3):.6f} (second zero of V-V(vac))")
print()

# The "mirror" condition: V(g0) = V(1) - 2*V(1) = -V(1) = -1/12
# Means the soliton core sits at the "anti-vacuum" energy
# 3g^4 - 4g^3 + 1 = 0 => (g-1)(3g^3 - g^2 - g - 1) = 0
# Cubic: 3g^3 - g^2 - g - 1 = 0
coeffs = [3, -1, -1, -1]
roots = np.roots(coeffs)
print(f"  V(g) = -V(1) = -1/12 at roots of 3g^3-g^2-g-1=0:")
for root in roots:
    if np.isreal(root) and root.real > 1:
        g_mirror = root.real
        print(f"    g_mirror = {g_mirror:.8f}")
        print(f"    V(g_mirror) = {V(g_mirror):.8f} (should be -1/12 = {-1/12:.8f})")
        print(f"    Distance from 4: {abs(g_mirror-4):.4f} ({abs(g_mirror-4)/4*100:.1f}%)")

# Also: V(g) = -n * V(1) for integer n?
print(f"\n  V(g0) / V(1) at key g0 values:")
for g0 in [2.0, 3.0, 3.5, 4.0, 4.5, 5.0]:
    ratio = V(g0) / V(1)
    print(f"    g0={g0:.1f}: V(g0)/V(1) = {ratio:.4f}")

# ===================================================================
print(f"\n{'='*70}")
print("  PART B: ALPHA_EFF AND f(g) AT SPECIAL POINTS")
print("="*70)

# At g0, what is alpha_eff and f(g0)?
print(f"\n  alpha_eff(g) = 2/(1 + eta_K*(g-1)^2), eta_K = 181/15")
print(f"  f(g) = 1 + 2*alpha_eff(g)*ln(g)")
print(f"\n  {'g0':>6s} {'alpha_eff':>10s} {'f(g0)':>10s} {'1/alpha':>10s}")
print(f"  {'-'*40}")

for g0 in [G0_E, G0_MU, 2.0, 3.0, 4.0, 5.0]:
    ae = alpha_eff(g0)
    fg = fk_func(g0)
    mark = ""
    if abs(g0-G0_E)<0.01: mark=" (e)"
    elif abs(g0-G0_MU)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau?)"
    print(f"  {g0:6.3f} {ae:10.6f} {fg:10.6f} {1/ae:10.4f}{mark}")

# Is there anything special about f(4)?
print(f"\n  Special values of f(g):")
print(f"    f(1) = {fk_func(1.0):.6f} (vacuum, always 1)")
print(f"    f(e) = {fk_func(np.e):.6f} (g=e)")
print(f"    f(4) = {fk_func(4.0):.6f}")
print(f"    f(phi^2) = {fk_func(PHI**2):.6f}")

# When does f(g) = 1 exactly (besides g=1)?
# f(g) = 1 + 2*alpha_eff*ln(g) = 1 when alpha_eff*ln(g) = 0
# i.e., g=1 (trivial) or alpha_eff=0 (g -> inf)
# When does f(g) have a maximum?
g_scan = np.linspace(0.5, 8, 1000)
f_scan = [fk_func(g) for g in g_scan]
f_max_idx = np.argmax(f_scan)
print(f"\n  f(g) maximum at g = {g_scan[f_max_idx]:.4f}, f_max = {f_scan[f_max_idx]:.6f}")

# ===================================================================
print(f"\n{'='*70}")
print("  PART C: g0 AS INTEGER / SIMPLE FRACTION")
print("="*70)

# The simplest hypothesis: g0_tau = 4 because it's the smallest
# INTEGER g0 > g0_mu ~ 1.465 that gives a distinct soliton.
#
# g0 = 2 gives too small a mass (r31 ~ 1016)
# g0 = 3 gives r31 ~ 2033 (too small)
# g0 = 4 gives r31 ~ 3477 (Koide!)
# g0 = 5 gives r31 ~ 6950 (too large)
#
# So g0=4 is the UNIQUE integer giving the Koide mass ratio.
# But this is numerological, not physical.

print(f"\n  Integer g0 test (with eta_K = 181/15):")

def get_A(r, g):
    m = (r>=120)&(r<=260)
    rf, tl = r[m], (g[m]-1)*r[m]
    if len(rf)<10: return np.nan
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    c,_,_,_ = np.linalg.lstsq(M, tl, rcond=None)
    return np.sqrt(c[0]**2+c[1]**2)

r_e, g_e, _ = solve_sol(G0_E)
A_e = get_A(r_e, g_e)

print(f"  {'g0':>6s} {'A_tail':>10s} {'r_31':>10s} {'m_tau(MeV)':>12s} {'Dev%':>8s}")
print(f"  {'-'*50}")
for g0_test in [2.0, 3.0, 4.0, 5.0, 6.0]:
    r_t, g_t, _ = solve_sol(g0_test)
    if r_t is None: continue
    A_t = get_A(r_t, g_t)
    if np.isnan(A_t): continue
    r31 = (A_t/A_e)**4
    m_tau = 0.51099895 * r31
    dev = abs(m_tau-1776.86)/1776.86*100
    print(f"  {g0_test:6.1f} {A_t:10.6f} {r31:10.2f} {m_tau:12.2f} {dev:8.3f}%")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART D: FIELD-SPACE WINDING AND TOPOLOGY")
print("="*70)

# In TGP, the field g takes values in (0, inf) with vacuum at g=1.
# The soliton profile goes from g0 at r=0 to g=1 at r=inf.
#
# The NUMBER OF TIMES g crosses 1 in the core is a topological invariant.
# For g0 > 1: the profile starts above 1, crosses through 1 (first node),
# may oscillate around 1 a few times in the core, then settles into
# the oscillating tail.
#
# Hypothesis: the GENERATION NUMBER is related to the number of nodes
# in the core (not the tail).

print(f"\n  Core node count (|g-1| crossings before tail regime):")
print(f"  {'g0':>6s} {'Core_r':>8s} {'Nodes':>6s} {'Note':>15s}")
print(f"  {'-'*40}")

for g0 in [G0_E, G0_MU, 2.0, 3.0, 4.0, 5.0, 6.0]:
    r, g, gp = solve_sol(g0)
    if r is None: continue
    # Core: where |g-1| > 0.05 (significant deviation from vacuum)
    core_end_idx = np.where(np.abs(g-1) < 0.05)[0]
    if len(core_end_idx) > 0:
        r_core = r[core_end_idx[0]]
    else:
        r_core = r[-1]

    # Count sign changes of (g-1) in core
    core_mask = r < max(r_core, 3.0)
    g_core = g[core_mask] - 1.0
    nodes = np.sum(np.abs(np.diff(np.sign(g_core))) > 0)

    mark = ""
    if abs(g0-G0_E)<0.01: mark="(e, gen 1)"
    elif abs(g0-G0_MU)<0.01: mark="(mu, gen 2)"
    elif abs(g0-4.0)<0.01: mark="(tau, gen 3?)"
    print(f"  {g0:6.3f} {r_core:8.2f} {nodes:6d} {mark:>15s}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART E: PHI-CHAIN: GENERATION SEQUENCE")
print("="*70)

# The phi-FP mechanism gives g0_mu = phi * g0_e.
# What if the SAME mechanism extends: g0_tau = phi^n * g0_e?
# Or: each generation is multiplied by phi?

print(f"\n  Golden-ratio chain from g0_e:")
print(f"  g0_e = {G0_E:.6f}")
print(f"  g0_mu = phi * g0_e = {G0_MU:.6f}")
print(f"  phi^2 * g0_e = {PHI**2*G0_E:.6f}")
print(f"  phi^3 * g0_e = {PHI**3*G0_E:.6f}")
print(f"  phi * g0_mu = {PHI*G0_MU:.6f}")
print()

# For tau: g0_tau = 4.000
# g0_tau / g0_e = 4.417
# g0_tau / g0_mu = 2.730
# These don't look like simple phi-powers.
print(f"  Ratios:")
print(f"    g0_tau / g0_e  = {4.0/G0_E:.4f}")
print(f"    g0_tau / g0_mu = {4.0/G0_MU:.4f}")
print(f"    phi^2 = {PHI**2:.4f}")
print(f"    phi^3 = {PHI**3:.4f}")
print(f"    e = {np.e:.4f}")
print(f"    phi*e = {PHI*np.e:.4f}")
print()

# But: g0_tau/g0_mu = 2.730 ~ e = 2.718 (0.43%)
# And: g0_tau/g0_e = 4.417 ~ phi^3 = 4.236 (4.3%)
# Neither is exact.

# Alternative: each generation ADDS a quantum
# g0_e = 1 - delta (below vacuum)
# g0_mu = 1 + Delta_mu
# g0_tau = 1 + Delta_tau
# Delta_tau / Delta_mu = ?

D_e = G0_E - 1  # = -0.0945
D_mu = G0_MU - 1  # = 0.465
D_tau = 4.0 - 1  # = 3.0

print(f"  Displacement from vacuum:")
print(f"    Delta_e   = g0_e - 1  = {D_e:.6f}")
print(f"    Delta_mu  = g0_mu - 1 = {D_mu:.6f}")
print(f"    Delta_tau = g0_tau - 1 = {D_tau:.6f}")
print(f"    Delta_tau / Delta_mu = {D_tau/D_mu:.4f}")
print(f"    Delta_mu / |Delta_e| = {D_mu/abs(D_e):.4f}")
print(f"    Delta_tau / |Delta_e| = {D_tau/abs(D_e):.4f}")
print()

# Delta_tau / Delta_mu = 6.451
# Delta_mu / |Delta_e| = 4.924 ~ phi^3 = 4.236 (16% off)

# ===================================================================
print(f"\n{'='*70}")
print("  PART F: THE V'(g0) = 0 AND STRUCTURAL ARGUMENTS")
print("="*70)

# V'(g) = g^2(1-g) = 0 at g=0 and g=1 (both trivial)
# V has maximum at g=1, inflection at g=2/3
# There is NO special point of V at g=4.

# But: V'(g)/g^2 = 1-g = 0 at g=1
# V''(g) = 2g - 3g^2 = 0 at g=0 and g=2/3
# V'''(g) = 2 - 6g = 0 at g=1/3

# The ODE: f(g)*g'' + (2/r)*g' = V'(g) = g^2(1-g)
# At g0=4: V'(4) = 16*(-3) = -48 (strong restoring force)
# At g0=2: V'(2) = 4*(-1) = -4
# At g0=1: V'(1) = 0

# Ratio V'(4)/V'(2) = 48/4 = 12 = eta_K^(0) !!
print(f"  V'(g) = g^2(1-g):")
print(f"    V'(2) = {2**2*(1-2):.1f}")
print(f"    V'(4) = {4**2*(1-4):.1f}")
print(f"    V'(4)/V'(2) = {4**2*(1-4)/(2**2*(1-2)):.1f}")
print(f"    This equals eta_K^(0) = alpha^2 * d = 12 !!")
print()
print(f"  More generally: V'(g0) = g0^2*(1-g0)")
print(f"  V'(g0)/V'(g0_ref) = (g0/g0_ref)^2 * (1-g0)/(1-g0_ref)")
print(f"  For g0=4, g0_ref=2: (4/2)^2 * (-3)/(-1) = 4*3 = 12 = alpha^2*d")
print()
print(f"  This is likely a COINCIDENCE, but it connects g0=4 to eta_K.")

# ===================================================================
print(f"\n{'='*70}")
print("  PART G: THE SOLITON CHARGE Q_sol")
print("="*70)

# Define the soliton charge as:
# Q_sol = int_0^inf f(g)*g'*dr (field-space displacement, weighted by f)
# This is a 1-form integral in field space.

print(f"\n  Soliton charge Q = int f(g)*g' dr:")
print(f"  {'g0':>6s} {'Q_sol':>10s} {'Q/pi':>10s} {'Nearest':>10s}")
print(f"  {'-'*45}")

for g0 in [G0_E, G0_MU, 2.0, 3.0, 4.0, 5.0, 6.0]:
    r, g, gp = solve_sol(g0)
    if r is None: continue
    fk_arr = np.array([fk_func(gi) for gi in g])
    # Only integrate over core (before tail oscillations dominate)
    core_end = np.where(np.abs(g-1)<0.05)[0]
    if len(core_end)>0:
        mask = np.arange(len(r)) < core_end[0]
    else:
        mask = np.ones(len(r), dtype=bool)

    Q = trapezoid(fk_arr[mask]*gp[mask], r[mask])
    Q_pi = Q/np.pi
    nearest = round(Q_pi*2)/2  # nearest half-integer
    mark = ""
    if abs(g0-G0_E)<0.01: mark=" (e)"
    elif abs(g0-G0_MU)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.3f} {Q:10.4f} {Q_pi:10.4f} {nearest:10.1f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART H: COMPOSITE STATE HYPOTHESIS")
print("="*70)

# Hypothesis: tau = composite of electron + muon in field space
# g0_tau = g0_e + g0_mu + 1? No: 0.905 + 1.465 + 1 = 3.37
# g0_tau = 2*g0_e + 2*g0_mu? No: 2*0.905 + 2*1.465 = 4.74
# g0_tau = g0_e + g0_mu + phi? 0.905 + 1.465 + 1.618 = 3.99!

print(f"\n  Composite hypothesis: g0_tau from g0_e and g0_mu?")
combos = [
    ("g0_e + g0_mu + phi",      G0_E + G0_MU + PHI),
    ("g0_e + 2*g0_mu",          G0_E + 2*G0_MU),
    ("g0_mu + phi^2",           G0_MU + PHI**2),
    ("2*g0_mu + g0_e",          2*G0_MU + G0_E),
    ("g0_e * g0_mu * phi",      G0_E * G0_MU * PHI),
    ("g0_mu^2 / g0_e",          G0_MU**2 / G0_E),
    ("phi^2 * g0_mu",           PHI**2 * G0_MU),
    ("e * g0_mu",               np.e * G0_MU),
    ("g0_e + g0_mu + sqrt(5)",  G0_E + G0_MU + np.sqrt(5)),
    ("phi*(g0_e + g0_mu)",      PHI*(G0_E + G0_MU)),
    ("(g0_e + g0_mu)^phi",      (G0_E + G0_MU)**PHI),
]

print(f"\n  {'Formula':<30s} {'Value':>10s} {'Dev%':>8s}")
print(f"  {'-'*52}")
for name, val in sorted(combos, key=lambda x: abs(x[1]-4.0)):
    dev = abs(val-4.0)/4.0*100
    mk = " ***" if dev < 0.5 else " **" if dev < 2 else ""
    print(f"  {name:<30s} {val:10.6f} {dev:8.3f}%{mk}")

# ===================================================================
print(f"\n{'='*70}")
print("  SUMMARY")
print("="*70)
print(f"""
  RESULT: No clean analytical derivation of g0_tau = 4 found.

  BEST OBSERVATIONS:
  1. g0_e + g0_mu + phi = {G0_E + G0_MU + PHI:.4f} (vs 4.000, dev {abs(G0_E+G0_MU+PHI-4)/4*100:.2f}%)
     -> Interpretable as "three generation sum" in field space
  2. e * g0_mu = {np.e*G0_MU:.4f} (vs 4.000, dev {abs(np.e*G0_MU-4)/4*100:.2f}%)
     -> Already noted in p135 (0.43% off)
  3. V'(4)/V'(2) = 12 = eta_K^(0) (interesting coincidence)
  4. g0 = 4 is the UNIQUE integer giving Koide mass ratio

  STATUS:
  g0_tau = 4 remains an INPUT PARAMETER.
  The theory has N_param = 3 (alpha_UV, g0_e, g0_tau).
  If g0_tau could be derived, N_param = 2 and predictivity doubles.

  HONEST ASSESSMENT:
  The value g0_tau = 4 is numerologically appealing (integer, V'(4)/V'(2)=12)
  but no PHYSICAL principle within TGP currently selects it.
  This is the single most important open problem in the mass sector.
""")
print("DONE")
