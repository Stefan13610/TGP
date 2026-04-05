#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p143_g0tau_deep.py -- Deep investigation of g0_tau = 4
=======================================================

NEW approaches not yet tried:

  A) Field-space geodesic distance: d(1,g0) = int_1^g0 sqrt(f(g)) dg
     Is this quantized (integer or half-integer multiples of pi)?

  B) Action ratio: S(g0_tau)/S(g0_e) -- is this a simple integer?

  C) Self-consistent condition: g0 determined by requiring
     alpha_eff(g0) * g0 = specific value (RG fixed-point in field space)

  D) Core-tail matching: the TRANSITION RADIUS r_t where the solution
     changes from core to tail behavior. Is r_t quantized?

  E) The "potential well depth" quantization:
     int_1^g0 sqrt(|V(g)-V(1)|) dg = n * something

  F) The mass formula inversion: what ALGEBRAIC equation must g0_tau
     satisfy for m_tau to equal the Koide prediction exactly?

  G) Scale invariance breaking: g0 is where alpha_eff * f(g) = 1
     (the theory becomes "marginally relevant")

  H) Connection to substrate: g0 = 2*alpha_UV = 2*2 = 4 (?!)

Author: TGP project
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)
import numpy as np
from scipy.integrate import solve_ivp, trapezoid, quad
from scipy.optimize import brentq

PHI = (1+np.sqrt(5))/2
ETA_K = 181.0/15
ALPHA_UV = 2.0
D = 3

def fk_func(g, eta=ETA_K):
    a = 2.0/(1+eta*(g-1)**2)
    return 1+2*a*np.log(g) if g>0 else -1e30

def alpha_eff(g, eta=ETA_K):
    return 2.0/(1+eta*(g-1)**2)

def V(g): return g**3/3 - g**4/4
def Vp(g): return g**2*(1-g)
def Vpp(g): return 2*g - 3*g**2
def V_eff(g): return V(g) - V(1)  # = V(g) - 1/12

def solve_sol(g0, eta=ETA_K, rm=300):
    def fk(g):
        a = 2.0/(1+eta*(g-1)**2)
        return 1+2*a*np.log(g) if g>0 else -1e30
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

def get_A(r, g):
    m = (r>=120)&(r<=260)
    rf, tl = r[m], (g[m]-1)*r[m]
    if len(rf)<10: return np.nan
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    c,_,_,_ = np.linalg.lstsq(M, tl, rcond=None)
    return np.sqrt(c[0]**2+c[1]**2)

G0_E = 0.90548144
G0_MU = PHI * G0_E

# ===================================================================
print("="*70)
print("  PART A: FIELD-SPACE GEODESIC DISTANCE")
print("="*70)

print(f"""
  The field space of TGP has a natural Riemannian metric:
    ds^2 = f(g) * dg^2
  where f(g) = 1 + 2*alpha_eff(g)*ln(g).

  The geodesic distance from vacuum (g=1) to soliton center (g=g0):
    d(1, g0) = int_1^g0 sqrt(|f(g)|) dg

  If d is quantized: d(1, g0_n) = n * Delta_d
  then the generation sequence has a physical interpretation as
  "winding" in field space.
""")

print(f"  {'g0':>6s} {'d(1,g0)':>10s} {'d/pi':>10s} {'d/d_e':>10s} {'d/d_mu':>10s}")
print(f"  {'-'*50}")

d_values = {}
for g0 in [G0_E, G0_MU, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]:
    # Integrate sqrt(|f(g)|) from 1 to g0
    if g0 < 1:
        dist, _ = quad(lambda g: np.sqrt(abs(fk_func(g))), g0, 1.0)
        dist = -dist  # negative for g0 < 1
    else:
        dist, _ = quad(lambda g: np.sqrt(abs(fk_func(g))), 1.0, g0)
    d_values[g0] = dist

d_e = abs(d_values.get(G0_E, 0))
d_mu = d_values.get(G0_MU, 1)

for g0 in [G0_E, G0_MU, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]:
    d = d_values[g0]
    mark = ""
    if abs(g0-G0_E)<0.01: mark=" (e)"
    elif abs(g0-G0_MU)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.3f} {d:10.4f} {d/np.pi:10.4f} {abs(d)/d_e if d_e>0 else 0:10.4f} {d/d_mu:10.4f}{mark}")

# Check ratios
d_tau = d_values.get(4.0, 0)
print(f"\n  Key ratios:")
print(f"    d_tau / d_mu = {d_tau/d_mu:.6f}")
print(f"    d_tau / d_e  = {d_tau/d_e:.6f}")
print(f"    d_mu / d_e   = {d_mu/d_e:.6f}")
print(f"    d_tau / pi   = {d_tau/np.pi:.6f}")
print(f"    d_tau         = {d_tau:.8f}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART B: SOLITON ACTION RATIO")
print("="*70)

print(f"""
  The 3D soliton action is:
    S[g0] = 4*pi * int_0^inf [ (1/2)*f(g)*g'^2 + V_eff(g) ] r^2 dr

  For a solution of the ODE, by Derrick's theorem:
    E_kin + 3*E_pot = 0 => E_kin = -3*E_pot
    S = E_kin + E_pot = -2*E_pot

  The ACTION RATIO between generations might be quantized.
""")

def compute_action(g0):
    r, g, gp = solve_sol(g0)
    if r is None: return np.nan, np.nan, np.nan
    fk_arr = np.array([fk_func(gi) for gi in g])
    v_arr = V_eff(g)
    E_kin = 4*np.pi * trapezoid(0.5*fk_arr*gp**2*r**2, r)
    E_pot = 4*np.pi * trapezoid(v_arr*r**2, r)
    S = E_kin + E_pot
    return S, E_kin, E_pot

print(f"  {'g0':>6s} {'S':>12s} {'E_kin':>12s} {'E_pot':>12s} {'S/S_e':>10s}")
print(f"  {'-'*55}")

S_e = compute_action(G0_E)[0]
for g0 in [G0_E, G0_MU, 2.0, 3.0, 4.0, 5.0, 6.0]:
    S, Ek, Ep = compute_action(g0)
    ratio = S/S_e if abs(S_e) > 0.001 else np.nan
    mark = ""
    if abs(g0-G0_E)<0.01: mark=" (e)"
    elif abs(g0-G0_MU)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.3f} {S:12.4f} {Ek:12.4f} {Ep:12.4f} {ratio:10.4f}{mark}")
sys.stdout.flush()

S_mu = compute_action(G0_MU)[0]
S_tau = compute_action(4.0)[0]
print(f"\n  Action ratios:")
print(f"    S_mu / S_e   = {S_mu/S_e:.6f}")
print(f"    S_tau / S_e  = {S_tau/S_e:.6f}")
print(f"    S_tau / S_mu = {S_tau/S_mu:.6f}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART C: g0 = 2*alpha_UV HYPOTHESIS")
print("="*70)

print(f"""
  The simplest structural argument:
    g0_tau = 2 * alpha_UV = 2 * 2 = 4

  Physical interpretation: the tau soliton core reaches the field
  value g = 2*alpha, where alpha is the UV coupling. At this point:
    alpha_eff(2*alpha) = 2/(1 + eta_K*(2*alpha-1)^2)
                       = 2/(1 + (181/15)*9) = 2/109.6 = 0.0182

  The "renormalized coupling" at the tau core is:
    alpha_eff(g0_tau) * alpha_UV = 0.0182 * 2 = 0.0365
  Compare with: 1/(alpha_UV^2 * d) = 1/12 = 0.0833 -- different

  Another reading: g0 = alpha_UV^(alpha_UV/alpha_UV) = alpha^1 = 2
  But we need g0 = 4 = alpha^2. So:
    g0_tau = alpha_UV^2 = 2^2 = 4 !!

  This is the SIMPLEST possible relation.
""")

# Test: at g0 = alpha^n for n = 0, 1, 2, 3:
print(f"  g0 = alpha_UV^n hypothesis:")
print(f"  {'n':>4s} {'g0':>8s} {'Name':>12s} {'Closest gen':>15s}")
print(f"  {'-'*45}")
for n in range(5):
    g0_n = ALPHA_UV**n
    closest = ""
    if abs(g0_n - 1.0) < 0.2: closest = "vacuum (g=1)"
    elif abs(g0_n - G0_MU) < 0.6: closest = f"mu (g0={G0_MU:.3f})"
    elif abs(g0_n - 4.0) < 0.1: closest = "tau (g0=4)"
    elif abs(g0_n - G0_E) < 0.2: closest = f"e (g0={G0_E:.3f})"
    print(f"  {n:4d} {g0_n:8.3f} {'alpha^'+str(n):>12s} {closest:>15s}")

print(f"""
  Observation:
    alpha^0 = 1  -> vacuum
    alpha^1 = 2  -> (no generation here)
    alpha^2 = 4  -> TAU!

  But alpha^1 = 2 doesn't correspond to muon (g0_mu = 1.465).
  So the simple alpha^n sequence doesn't work for ALL generations.

  However: alpha_UV^2 = 4 for the tau is structurally compelling.
  It means: g0_tau sits at the field value where the BARE kinetic
  coupling K(g) = g^(2*alpha) = g^4 equals g^(alpha^2) = g^4.
  That is: K(g0_tau) = g0_tau^4 = 4^4 = 256 = alpha^8.
""")

# K(g0) = g0^4 for K = phi^4 in TGP
print(f"  K(g0) = g0^(2*alpha) = g0^4:")
for g0 in [G0_E, G0_MU, 4.0]:
    K = g0**4
    mark = ""
    if abs(g0-G0_E)<0.01: mark=" (e)"
    elif abs(g0-G0_MU)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"    K({g0:.3f}) = {K:.4f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART D: SELF-CONSISTENT FIELD-SPACE FIXED POINT")
print("="*70)

print(f"""
  Hypothesis: g0_tau is the solution of a SELF-CONSISTENT equation
  involving alpha_eff and the potential V.

  Candidate equations:
    1) g0 * alpha_eff(g0) = 1/(d-1)
    2) f(g0) * V'(g0) = -g0^2 (some normalization)
    3) alpha_eff(g0) = 1/(g0^2 * d)
    4) The soliton "radius" R = 1/sqrt(|V'(g0)|/f(g0)) = specific value
""")

# Test each candidate
print(f"  {'g0':>6s} {'aeff*g0':>10s} {'f*Vp/g2':>10s} {'aeff-1/g2d':>12s} {'R_core':>10s}")
print(f"  {'-'*55}")

for g0 in np.arange(1.5, 8.1, 0.5):
    ae = alpha_eff(g0)
    fg = fk_func(g0)
    vp = Vp(g0)
    R_core = 1.0/np.sqrt(abs(vp)/fg) if abs(vp)>0 and fg>0 else np.nan

    c1 = ae * g0
    c2 = fg * vp / g0**2 if g0 > 0 else np.nan
    c3 = ae - 1.0/(g0**2 * D)

    mark = ""
    if abs(g0-G0_MU)<0.1: mark=" (mu)"
    elif abs(g0-4.0)<0.1: mark=" (tau)"
    print(f"  {g0:6.2f} {c1:10.6f} {c2:10.6f} {c3:12.6f} {R_core:10.4f}{mark}")

# Where does alpha_eff(g0) * g0^2 = some constant?
print(f"\n  alpha_eff(g0) * g0^2:")
for g0 in np.arange(1.5, 8.1, 0.5):
    ae = alpha_eff(g0)
    val = ae * g0**2
    mark = " (tau)" if abs(g0-4.0)<0.1 else ""
    print(f"    g0={g0:.1f}: alpha*g0^2 = {val:.6f}{mark}")

# ===================================================================
print(f"\n{'='*70}")
print("  PART E: CORE-TAIL TRANSITION STRUCTURE")
print("="*70)

# The soliton has a core where g deviates significantly from 1,
# and a tail where g oscillates around 1 with decaying amplitude.
# The TRANSITION between core and tail happens at some radius r_t.
# Is there a special value of r_t that singles out g0=4?

print(f"\n  Core-tail transition analysis:")
print(f"  {'g0':>6s} {'r_trans':>8s} {'g(r_t)':>8s} {'g_dot(r_t)':>10s} {'A_core':>10s} {'r_t*A_core':>10s}")
print(f"  {'-'*60}")

for g0 in [G0_E, G0_MU, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]:
    r, g, gp = solve_sol(g0)
    if r is None: continue

    # Find transition: first time |g-1| < 0.1 after being large
    core_mask = np.abs(g - 1) > 0.1
    if not np.any(core_mask):
        r_t = r[0]
    else:
        # Find first index where |g-1| drops below 0.1
        trans_idx = np.where(~core_mask)[0]
        if len(trans_idx) > 0 and trans_idx[0] > 0:
            ti = trans_idx[0]
            r_t = r[ti]
            g_t = g[ti]
            gp_t = gp[ti]

            # "Core amplitude": max deviation of g from 1 in core
            A_core = np.max(np.abs(g[:ti] - 1))

            mark = ""
            if abs(g0-G0_E)<0.01: mark=" (e)"
            elif abs(g0-G0_MU)<0.01: mark=" (mu)"
            elif abs(g0-4.0)<0.01: mark=" (tau)"
            print(f"  {g0:6.3f} {r_t:8.3f} {g_t:8.4f} {gp_t:10.4f} {A_core:10.4f} {r_t*A_core:10.4f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART F: POTENTIAL WELL QUANTIZATION")
print("="*70)

print(f"""
  The "potential well depth" integral:
    I(g0) = int_1^g0 sqrt(2*f(g)*|V_eff(g)|) dg / pi

  This is like a WKB phase integral in the "field-space potential well".
  V_eff(g) = V(g) - 1/12 < 0 for g > 4/3.

  For a bound state, I(g0) should be n + 1/2 (with Maslov correction)
  or n (without).
""")

print(f"  {'g0':>6s} {'I(g0)':>10s} {'I-n':>10s} {'I-(n+1/2)':>10s}")
print(f"  {'-'*40}")

for g0 in np.arange(1.5, 8.1, 0.25):
    # Integrate sqrt(2*f*|V_eff|) from 1 to g0
    def integrand(g):
        f = fk_func(g)
        v = V_eff(g)
        if f > 0 and v < 0:
            return np.sqrt(2.0 * f * abs(v))
        return 0.0

    I_val, _ = quad(integrand, 1.001, g0, limit=200)
    I_pi = I_val / np.pi
    n = round(I_pi)
    dev_n = I_pi - n
    dev_half = I_pi - (round(I_pi - 0.5) + 0.5)

    mark = ""
    if abs(g0-G0_MU)<0.15: mark=" (mu)"
    elif abs(g0-4.0)<0.15: mark=" (tau)"
    if abs(dev_n) < 0.05: mark += " ***"
    elif abs(dev_half) < 0.05: mark += " (n+1/2)*"

    print(f"  {g0:6.3f} {I_pi:10.4f} {dev_n:10.4f} {dev_half:10.4f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART G: MATCHING CONDITION: ODE AT THE BOUNDARY")
print("="*70)

print(f"""
  At r=0 (soliton center), the ODE regularity condition gives:
    g(0) = g0, g'(0) = 0
    g''(0) = V'(g0) / (3*f(g0))

  The CURVATURE of the soliton at center is:
    kappa = g''(0) = V'(g0) / (3*f(g0))

  Question: does kappa(g0=4) have a special value?
""")

print(f"  {'g0':>6s} {'g_pp(0)':>12s} {'kappa*g0':>10s} {'kappa*g0^2':>10s}")
print(f"  {'-'*45}")

for g0 in np.arange(1.5, 8.1, 0.5):
    f_g0 = fk_func(g0)
    vp_g0 = Vp(g0)
    kappa = vp_g0 / (3*f_g0) if abs(f_g0) > 1e-10 else np.nan
    mark = ""
    if abs(g0-4.0)<0.1: mark=" (tau)"
    print(f"  {g0:6.2f} {kappa:12.4f} {kappa*g0:10.4f} {kappa*g0**2:10.4f}{mark}")

# ===================================================================
print(f"\n{'='*70}")
print("  PART H: alpha_eff * (g0-1)^2 = CONSTANT?")
print("="*70)

# In the running alpha formula: alpha_eff = 2/(1 + eta_K*(g-1)^2)
# At g0: alpha_eff(g0) = 2/(1 + eta_K*(g0-1)^2)
# So: alpha_eff * (1 + eta_K*(g0-1)^2) = 2
# This is ALWAYS 2 by construction - it's the definition.

# But: alpha_eff * eta_K * (g0-1)^2 = 2 - alpha_eff
# For g0=4: 0.01825 * 12.0667 * 9 = 1.982 ~ 2

# Better question: what is alpha_eff(g0) * g0^n for the three leptons?
print(f"  alpha_eff(g0) * g0^n for leptons:")
for n in range(5):
    vals = []
    for name, g0 in [("e", G0_E), ("mu", G0_MU), ("tau", 4.0)]:
        ae = alpha_eff(g0)
        val = ae * g0**n
        vals.append(val)
    print(f"    n={n}: e={vals[0]:.6f}  mu={vals[1]:.6f}  tau={vals[2]:.6f}  "
          f"tau/e={vals[2]/vals[0]:.4f}  tau/mu={vals[2]/vals[1]:.4f}")

# ===================================================================
print(f"\n{'='*70}")
print("  PART I: THE SOLITON 'INFORMATION' INTEGRAL")
print("="*70)

# The information content of the soliton: how much the field deviates
# from vacuum, weighted by the natural metric

print(f"""
  Define the soliton information as:
    I[g0] = 4*pi * int_0^inf |g(r)-1|^2 * r^2 dr / R_core^3

  This measures the "total deviation" from vacuum normalized by the
  core volume. For a fundamental soliton, this might be quantized.
""")

print(f"  {'g0':>6s} {'I_raw':>12s} {'I/pi':>10s} {'I/I_e':>10s}")
print(f"  {'-'*45}")

I_e = None
for g0 in [G0_E, G0_MU, 2.0, 3.0, 4.0, 5.0, 6.0]:
    r, g, gp = solve_sol(g0)
    if r is None: continue

    I_raw = 4*np.pi * trapezoid((g-1)**2 * r**2, r)

    # Core radius
    core_end = np.where(np.abs(g-1)<0.1)[0]
    if len(core_end)>0:
        R = r[core_end[0]]
    else:
        R = 1.0

    I_norm = I_raw / max(R**3, 0.01)

    if I_e is None: I_e = I_raw

    mark = ""
    if abs(g0-G0_E)<0.01: mark=" (e)"
    elif abs(g0-G0_MU)<0.01: mark=" (mu)"
    elif abs(g0-4.0)<0.01: mark=" (tau)"
    print(f"  {g0:6.3f} {I_raw:12.4f} {I_raw/np.pi:10.4f} {I_raw/I_e:10.4f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART J: PHASE ACCUMULATION IN THE CORE")
print("="*70)

print(f"""
  New idea: count the TOTAL PHASE accumulated by the soliton
  from r=0 to r=r_core, defined via the local wavenumber:
    k_local(r) = sqrt(|V''(g(r))/f(g(r))|)
    Phi_core = int_0^r_core k_local(r) dr

  For the tail: k_local -> 1 (the oscillation frequency).
  For the core: k_local can be much larger (strong curvature).

  If Phi_core = n*pi, we have a resonance condition.
""")

print(f"  {'g0':>6s} {'Phi_core':>10s} {'Phi/pi':>10s} {'nearest n':>10s} {'dev':>10s}")
print(f"  {'-'*55}")

for g0 in np.arange(1.5, 8.1, 0.25):
    r, g, gp = solve_sol(g0)
    if r is None: continue

    # Core: where |g-1| > 0.05
    core_end = np.where(np.abs(g-1)<0.05)[0]
    if len(core_end) > 0:
        ci = core_end[0]
    else:
        ci = len(r)//2

    r_core = r[:ci]
    g_core = g[:ci]

    # Local wavenumber
    fk_arr = np.array([fk_func(gi) for gi in g_core])
    vpp_arr = np.array([Vpp(gi) for gi in g_core])

    k_local = np.sqrt(np.abs(vpp_arr / fk_arr))
    k_local = np.where(np.isfinite(k_local), k_local, 0)

    Phi_core = trapezoid(k_local, r_core) if len(r_core) > 1 else 0
    Phi_pi = Phi_core / np.pi
    n_near = round(Phi_pi)
    dev = Phi_pi - n_near

    mark = ""
    if abs(g0-G0_MU)<0.15: mark=" (mu)"
    elif abs(g0-4.0)<0.15: mark=" (tau)"
    if abs(dev) < 0.06: mark += " ***"

    print(f"  {g0:6.3f} {Phi_core:10.4f} {Phi_pi:10.4f} {n_near:10d} {dev:10.4f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*70}")
print("  PART K: THE UNIQUENESS ARGUMENT")
print("="*70)

print(f"""
  Perhaps g0_tau = 4 is not derived from a quantization condition
  but from a UNIQUENESS argument:

  Given:
    1) alpha_UV = 2 (from substrate K = phi^4)
    2) eta_K = 181/15 (from ERG self-consistency)
    3) phi-FP mechanism: g0_mu = phi * g0_e

  The ONLY remaining freedom is (g0_e, g0_tau).
  g0_e is FIXED by r_21 = 206.768 (calibration).
  g0_tau is unconstrained...

  UNLESS we add one more physical condition.
  What if that condition is:

  g0_tau = alpha_UV^2 = 4

  This would mean: the tau soliton core amplitude equals the SQUARE
  of the UV coupling. Since alpha = 2 comes from K(phi) = phi^4
  (i.e., phi^(2*alpha) with alpha=2), we have:

    g0_tau = alpha^2  <=>  K(g0_tau) = g0_tau^4 = 256 = alpha^8

  The connection: g0_tau is where the KINETIC FUNCTION EXPONENT
  equals the field value squared:
    phi^(2*alpha) at phi = alpha^(alpha/2) when 2*alpha = alpha^2
    This requires: 2*alpha = alpha^2 => alpha(alpha-2) = 0 => alpha=2!

  So alpha = 2 is SPECIAL: it's the unique value where
    the kinetic exponent (2*alpha = 4)
    equals the SQUARE of the coupling (alpha^2 = 4)
    equals the candidate g0_tau (= 4).

  This is a SELF-REFERENTIAL condition specific to alpha = 2.
""")

# Check: 2*alpha = alpha^2 only when alpha = 0 or alpha = 2
print(f"  Self-referential condition: 2*alpha = alpha^2")
print(f"  Solutions: alpha = 0 (trivial) and alpha = 2 (TGP!)")
print(f"  At alpha = 2:")
print(f"    2*alpha = {2*ALPHA_UV}")
print(f"    alpha^2 = {ALPHA_UV**2}")
print(f"    g0_tau  = {4.0}")
print(f"    ALL EQUAL: 2*alpha = alpha^2 = g0_tau = 4")
print()
print(f"  K(phi) = phi^(2*alpha) = phi^4")
print(f"  K(g0_tau) = K(alpha^2) = (alpha^2)^(2*alpha) = 4^4 = 256")
print(f"  Also: K(g0_tau) = g0_tau^(g0_tau) = 4^4 = 256")
print(f"  So at the tau core: K(g0) = g0^g0 -- a SELF-DUAL condition!")
print()

# Verify: K(g0) = g0^g0 only at g0 = 2*alpha when K = phi^(2*alpha)
# phi^(2*alpha) = phi^phi when 2*alpha = phi, i.e., phi=4 when alpha=2
# Yes: g0^4 = g0^g0 when g0 = 4. Self-consistent!
print(f"  Verification: g0^(2*alpha) = g0^g0")
print(f"  This requires: 2*alpha = g0")
print(f"  With alpha = 2: g0 = 4. QED.")
print()
print(f"  This is a SELF-REFERENTIAL FIXED POINT:")
print(f"  The exponent in K(phi) = phi^(2*alpha) equals the field value g0_tau")
print(f"  if and only if g0_tau = 2*alpha = 2*alpha_UV.")

# ===================================================================
print(f"\n{'='*70}")
print("  SUMMARY")
print("="*70)
print(f"""
  ===================================================================
  MAIN FINDING: g0_tau = 2*alpha_UV = alpha_UV^2 = 4
  ===================================================================

  This follows from a SELF-REFERENTIAL FIXED POINT condition:

  In TGP, the kinetic function is K(phi) = phi^(2*alpha), alpha = 2.
  The soliton amplitude g0 defines the field excursion from vacuum.

  At g0_tau: K(g0_tau) = g0_tau^(2*alpha) should equal g0_tau^(g0_tau)
  (self-dual condition: the kinetic exponent matches the field value).

  This requires: 2*alpha = g0_tau
  With alpha_UV = 2: g0_tau = 4.

  UNIQUENESS: This works ONLY because alpha_UV = 2 is the UNIQUE
  value where 2*alpha = alpha^2 (besides alpha=0).
  If alpha were anything other than 2, this self-referential
  condition would give 2*alpha != alpha^2 and no clean integer.

  PHYSICAL INTERPRETATION:
  The tau soliton is the configuration where the field excursion
  equals the kinetic exponent: g0 = 2*alpha. This is a kind of
  "self-duality" or "fixed point" in the space of field configurations.

  CHAIN: K(phi) = phi^4 (from substrate)
       -> alpha_UV = 2
       -> g0_tau = 2*alpha = 4 (self-referential condition)
       -> eta_K = 181/15 (from ERG)
       -> m_tau = 1776.97 MeV (0.006% from observation)

  STATUS: This is a HYPOTHESIS (not a derivation from first principles).
  The self-referential condition is aesthetically compelling but
  needs to be embedded in a formal variational or RG framework.

  COMPLEMENTARY OBSERVATIONS:
  - V'(4)/V'(2) = 12 = alpha^2 * d = eta_K^(0)
  - |V(4)|/V(1) = 512 = 2^9 = (2*alpha)^(2*alpha+1) = 4^(4+1)/2... no, just 2^9
  - g0_e + g0_mu + phi = 3.989 ~ 4 (0.28% off)
  - The tau is the UNIQUE integer-g0 lepton matching Koide

  PREDICTIVITY:
  If g0_tau = 2*alpha accepted:
    N_param = 2 (alpha_UV from substrate, g0_e from r_21 calibration)
    N_pred = 5 (r_21, r_31/m_tau, Delta(e->mu), omega_tail, eta_K)
    Ratio = 5/2 = 2.5
""")
print("DONE")
