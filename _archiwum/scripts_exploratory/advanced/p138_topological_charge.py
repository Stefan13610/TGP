#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p138_topological_charge.py - Topological charge and g0_tau quantization
=========================================================================
In TGP, pi_1(C_sol) = Z_2 gives fermion statistics (spin-statistics).
But in 3D, solitons can also carry a TOPOLOGICAL CHARGE related to:
  1) The winding number of the map phi: R^3 -> field space
  2) The "Noether charge" from the internal symmetry
  3) The Derrick scaling charge

Key idea: if the soliton has a conserved topological charge Q_top,
then g0 might be quantized: g0 = 1 + Q_top * Delta_g for some unit.

Also investigate:
  A) Does the effective action S_eff have critical points at g0=4?
  B) Is g0=4 related to the number of spacetime dimensions (3+1)?
  C) Semiclassical quantization: Bohr-Sommerfeld for soliton vibrations

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

def solve_soliton(g0, eta_K=ETA_K, rm=300):
    def fk(g):
        a = 2.0 / (1 + eta_K * (g - 1)**2)
        return 1 + 2*a*np.log(g) if g > 0 else -1e30
    def Vp(g): return g**2*(1-g)
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
    return r, s.sol(r)[0], s.sol(r)[1]

def fk_func(g, eta_K=ETA_K):
    a = 2.0 / (1 + eta_K * (g - 1)**2)
    return 1 + 2*a*np.log(g) if g > 0 else -1e30

def V_func(g): return g**3/3.0 - g**4/4.0

def extract_tail(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    return np.sqrt(coeff[0]**2 + coeff[1]**2)

r_e, g_e, _ = solve_soliton(G0_E)
A_e = extract_tail(r_e, g_e)

# =====================================================================
print("="*70)
print("  PART A: TOPOLOGICAL CHARGE OF TGP SOLITON")
print("="*70)

# In TGP with K(phi) = phi^4, the field space is phi in (0, inf).
# The soliton maps R^3 -> (0, inf) with boundary phi -> phi_vac = Phi_0.
# In reduced variable g = phi/Phi_0: g -> 1 at infinity.
#
# The topological charge is:
#   Q = (1/2) * int_0^inf |g'(r)| dr = (1/2) * |g(0) - g(inf)|
# (for monotonic profiles, which ours are NOT since g oscillates around 1)
#
# Actually, the "charge" is more subtle for oscillating tails.
# Better definition: the NUMBER OF NODES of (g-1) in the core region
# before the tail oscillations begin.
#
# For the soliton starting at g0 > 1:
# - It starts at g0, decreases through 1, overshoots below 1,
#   comes back through 1 again, etc.
# - The "core" is the region before the oscillating tail.

print(f"\n  Soliton profiles: counting core features")
print(f"  {'g0':>6s} {'g_min':>10s} {'g_max':>10s} {'r_core':>10s} {'n_nodes':>10s}")
print(f"  {'-'*50}")

for g0 in [G0_E, G0_MU, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0]:
    r, g, gp = solve_soliton(g0)
    if r is None: continue

    # Find core region: where |g-1| > 0.01
    core_mask = np.abs(g - 1.0) > 0.01
    if not np.any(core_mask):
        r_core = 0
    else:
        # Find FIRST r where soliton enters tail (g oscillates close to 1)
        # Simple: find r where g first reaches within 0.01 of 1
        cross = np.where(np.abs(g - 1.0) < 0.01)[0]
        if len(cross) > 0:
            r_core = r[cross[0]]
        else:
            r_core = r[-1]

    # Stats in core
    core_g = g[r < r_core] if r_core > 0 else g[:100]
    g_min = np.min(core_g) if len(core_g) > 0 else np.nan
    g_max = np.max(core_g) if len(core_g) > 0 else np.nan

    # Count zero crossings of (g-1) in core
    diff = g[r < max(r_core, 5)] - 1.0
    nodes = np.sum(np.abs(np.diff(np.sign(diff))) > 0) if len(diff) > 1 else 0

    mark = ""
    if abs(g0 - G0_E) < 0.01: mark = " (e)"
    elif abs(g0 - G0_MU) < 0.01: mark = " (mu)"
    elif abs(g0 - 4.0) < 0.01: mark = " (tau)"
    print(f"  {g0:6.3f} {g_min:10.4f} {g_max:10.4f} {r_core:10.2f} {nodes:10d}{mark}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART B: BOHR-SOMMERFELD QUANTIZATION OF SOLITON")
print(f"{'='*70}")

# The soliton in field space goes from g0 to 1 (and oscillates around 1).
# In the "radial" direction in field space, the effective Hamiltonian is:
#   H = f(g)*g'^2/2 + V_eff(g)
# where V_eff = V(g) - V(1) = g^3/3 - g^4/4 - 1/12
#
# Bohr-Sommerfeld: int p dg = (n + 1/2) * pi
# where p = sqrt(2*f(g)*(E - V_eff(g)))
#
# For a zero-energy soliton (E_tot ~ 0 for electron):
# p = sqrt(-2*f(g)*V_eff(g))
# This is real only where V_eff(g) < 0, i.e., g > 1 (and some g < 1).

# V_eff(g) = g^3/3 - g^4/4 - 1/12
# V_eff(0) = -1/12 < 0
# V_eff(1) = 0
# V_eff(g) < 0 for g < 1 and g > 4/3

# Let's compute the Bohr-Sommerfeld integral for different g0
print(f"\n  V_eff(g) = g^3/3 - g^4/4 - 1/12")
print(f"  V_eff = 0 at g = 1 (and g = -1/3)")
print(f"  V_eff has local max at g=1: V_eff(1) = 0")
print(f"  V_eff < 0 for g > 4/3 and for g < 1 (except g=1)")

def V_eff(g):
    return g**3/3.0 - g**4/4.0 - 1.0/12.0

# Classical turning points for V_eff = E:
# At E = 0: g = 1 (inner) and V_eff(g_outer) = 0
# V_eff has another zero at... let's find it
print(f"\n  V_eff at key points:")
for g_test in [0, 0.5, 0.8, 1.0, 1.333, 2.0, 3.0, 4.0, 5.0]:
    print(f"    V_eff({g_test}) = {V_eff(g_test):.6f}")

# Bohr-Sommerfeld for each g0:
# int_{g=1}^{g0} sqrt(2*f(g)*|V_eff(g)|) dg / pi
print(f"\n  Bohr-Sommerfeld integral I_BS = (1/pi) * int_1^g0 sqrt(2*f*|V_eff|) dg")
print(f"\n  {'g0':>6s} {'I_BS':>10s} {'I_BS - n':>10s} {'closest_n':>10s}")
print(f"  {'-'*42}")

for g0 in np.arange(1.5, 8.5, 0.25):
    # Integrate from g=1 to g=g0
    g_grid = np.linspace(1.001, g0, 500)
    integrand = np.zeros_like(g_grid)
    for i, g in enumerate(g_grid):
        f_val = fk_func(g)
        v_val = V_eff(g)
        if f_val > 0 and v_val < 0:
            integrand[i] = np.sqrt(2.0 * f_val * abs(v_val))
        else:
            integrand[i] = 0.0

    I_bs = trapezoid(integrand, g_grid) / np.pi
    n_closest = round(I_bs)
    dev = I_bs - n_closest

    mark = ""
    if abs(g0 - G0_MU) < 0.15: mark = " (mu)"
    elif abs(g0 - 4.0) < 0.15: mark = " (tau)"
    if abs(dev) < 0.05: mark += " ***"
    elif abs(I_bs - (n_closest + 0.5)) < 0.05: mark += " (n+1/2)"

    print(f"  {g0:6.3f} {I_bs:10.4f} {dev:10.4f} {n_closest:10d}{mark}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART C: PHASE SPACE INTEGRAL (FIELD-SPACE ACTION)")
print(f"{'='*70}")

# Alternative: the FIELD-SPACE ACTION
# S_field = int_1^g0 sqrt(2 * f(g) * |V_eff(g)|) dg
# This is like the WKB action in quantum mechanics.
# If quantized: S_field = n * pi (Bohr-Sommerfeld)
# or S_field = (n + 1/2) * pi (with Maslov correction)

# For the radial ODE, the ACTUAL action involves r-dependence.
# The field-space action S_r at radius r is:
# S_r = int f(g)*g'(r)^2 * r^2 dr (kinetic action integrated over r)
# This is just E_kin (from p135).

# But there's a more physically motivated quantity:
# The REDUCED ACTION per radial shell:
# sigma(g0) = int_0^R_core f(g)*g'^2 dr  (1D action without 4*pi*r^2)

print(f"\n  1D action sigma = int_0^R_core f(g)*g'^2 dr (without r^2 weight)")
print(f"\n  {'g0':>6s} {'sigma':>10s} {'sigma/pi':>10s} {'closest':>10s}")
print(f"  {'-'*42}")

for g0 in np.arange(1.5, 8.5, 0.25):
    r, g, gp = solve_soliton(g0)
    if r is None: continue
    fk_arr = np.array([fk_func(gi) for gi in g])
    integrand = fk_arr * gp**2

    # Find core radius (where |g-1| > 0.01)
    core_end = np.where(np.abs(g - 1.0) < 0.01)[0]
    if len(core_end) > 0:
        r_core = r[core_end[0]]
        mask = r < r_core
    else:
        mask = np.ones(len(r), dtype=bool)

    sigma = trapezoid(integrand[mask], r[mask])
    s_over_pi = sigma / np.pi
    n_closest = round(s_over_pi)
    dev = s_over_pi - n_closest

    mark = ""
    if abs(g0 - G0_MU) < 0.15: mark = " (mu)"
    elif abs(g0 - 4.0) < 0.15: mark = " (tau)"
    if abs(dev) < 0.1: mark += " *"
    if abs(s_over_pi - (n_closest + 0.5)) < 0.1: mark += " (n+1/2)"

    print(f"  {g0:6.3f} {sigma:10.4f} {s_over_pi:10.4f} {n_closest:10d}{mark}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART D: g0 = 4 FROM STRUCTURAL ARGUMENT")
print(f"{'='*70}")

# Structural argument: in TGP, V(g) = g^3/3 - g^4/4
# The potential has V(0)=0, V(1)=1/12, and V->-inf for g->inf.
#
# The inflection point: V''(g) = 2g - 3g^2 = 0 at g = 2/3
# The maximum of V: V'(g) = g^2(1-g) = 0 at g=0 and g=1
# V(1) = 1/12 is the absolute maximum for g > 0.
#
# For g > 1: V(g) = g^3/3 - g^4/4 = g^3(1/3 - g/4)
# V(g) = 0 at g = 4/3 (and g=0)
# V(g) < 0 for g > 4/3
#
# SPECIAL VALUE: V(g) = -V(1) = -1/12 when:
# g^3/3 - g^4/4 = -1/12
# g^3(4 - 3g) = -1  (multiply by 12)
# 4g^3 - 3g^4 = -1
# 3g^4 - 4g^3 - 1 = 0
# Factor: (g-1)(3g^3 - g^2 - g - 1) = 0
# So g=1 gives V=1/12, and the cubic 3g^3 - g^2 - g - 1 = 0 gives V=-1/12.

print(f"\n  Structural argument from V(g):")
print(f"  V(g) = g^3/3 - g^4/4")
print(f"  V(0) = 0, V(1) = 1/12 (vacuum energy)")
print(f"  V(4/3) = 0 (second zero)")
print(f"")

# Solve 3g^3 - g^2 - g - 1 = 0
coeffs = [3, -1, -1, -1]
roots = np.roots(coeffs)
print(f"  V(g) = -V(1) = -1/12 at roots of 3g^3 - g^2 - g - 1 = 0:")
for root in roots:
    if np.isreal(root):
        g_val = root.real
        v_val = V_func(g_val)
        print(f"    g = {g_val:.6f}, V(g) = {v_val:.6f} (should be -1/12 = {-1/12:.6f})")

# Check: at what g does |V(g)| = some multiple of V(1)?
print(f"\n  |V(g)| / V(1) at various g:")
for g_test in [2.0, 3.0, 3.5, 4.0, 4.5, 5.0]:
    v = V_func(g_test)
    ratio = abs(v) / (1.0/12.0)
    print(f"    g={g_test}: V = {v:.4f}, |V|/V(1) = {ratio:.2f}")

# g=4: V(4) = 64/3 - 256/4 = 64/3 - 64 = -128/3 = -42.667
# |V(4)|/V(1) = 42.667 / 0.0833 = 512
# 512 = 2^9 = 8^3 = (2*g0)^3? Interesting!
print(f"\n  |V(4)|/V(1) = {abs(V_func(4.0))/(1.0/12.0):.0f} = 2^9 = 512")
print(f"  This is (2*g0_tau)^3 / (2*1)^3 if we associate g0=1 with vacuum.")

# Actually: |V(g)| = |g^3/3 - g^4/4| = g^3*|1/3 - g/4| for g > 4/3: = g^3*(g/4 - 1/3) = g^4/4 - g^3/3
# |V(g)|/V(1) = (g^4/4 - g^3/3)/(1/12) = 12*(g^4/4 - g^3/3) = 3*g^4 - 4*g^3
# At g=4: 3*256 - 4*64 = 768 - 256 = 512 = 2^9. Indeed!

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART E: SOLITON MASS FORMULA DECOMPOSITION")
print(f"{'='*70}")

# The mass (from tail amplitude) can be decomposed:
# m ~ A_tail^4 = [A_core * exp(-mu * R_core)]^4 (heuristic)
# where A_core depends on g0 and R_core depends on the potential gradient.
#
# For large g0 with running alpha:
# alpha_eff ~ 2/(eta_K * (g0-1)^2) << 1
# f(g0) ~ 1 + 2*alpha_eff*ln(g0) ~ 1 (barely above 1)
#
# The soliton is then approximately governed by the STANDARD equation:
# g'' + (2/r)*g' = V'(g)  (since f ~ 1)
#
# In this limit, the core radius R_core ~ sqrt(1/|V'(g0)|)
# |V'(g0)| = g0^2*(g0-1) ~ g0^3 for large g0
# So R_core ~ g0^{-3/2}

# The amplitude from the core: A_core ~ (g0 - 1) * R_core^2 ~ g0 * g0^{-3} = g0^{-2}
# And A_tail ~ A_core * exp(-R_core) * R_core (rough)

# This is too rough. Let's compute A(g0) numerically more carefully:
print(f"\n  A(g0) vs g0 in the large-g0 regime (f ~ 1):")
print(f"  {'g0':>6s} {'A':>10s} {'ln(A)':>10s} {'A*g0^2':>10s} {'A/g0^0.7':>10s}")
print(f"  {'-'*50}")

for g0 in np.arange(2.0, 10.0, 0.5):
    r, g, _ = solve_soliton(g0)
    if r is None: continue
    A = extract_tail(r, g)
    if np.isnan(A): continue
    print(f"  {g0:6.2f} {A:10.6f} {np.log(A):10.4f} {A*g0**2:10.4f} {A/g0**0.7:10.6f}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART F: TEST eta_K = 12 EXACT with g0_tau = 4 EXACT")
print(f"{'='*70}")

# Final definitive test: what mass does eta_K=12 (exact), g0_tau=4 (exact) give?

# First find g0_e for eta_K = 12 exactly
def find_g0e_for_eta(eta_K):
    def get_r21(g0_e):
        g0_mu = PHI * g0_e
        r_e, g_e, _ = solve_soliton(g0_e, eta_K)
        if r_e is None or r_e[-1] < 250: return 1e6
        A_e = extract_tail(r_e, g_e)
        r_mu, g_mu, _ = solve_soliton(g0_mu, eta_K)
        if r_mu is None or r_mu[-1] < 250: return 1e6
        A_mu = extract_tail(r_mu, g_mu)
        if np.isnan(A_e) or np.isnan(A_mu) or A_e < 1e-15: return 1e6
        return (A_mu/A_e)**4
    for g_lo in np.arange(0.88, 0.93, 0.005):
        try:
            v_lo = get_r21(g_lo) - 206.768
            v_hi = get_r21(g_lo+0.005) - 206.768
            if v_lo * v_hi < 0:
                return brentq(lambda g: get_r21(g)-206.768, g_lo, g_lo+0.005, xtol=1e-9)
        except:
            pass
    return None

print(f"\n  Test 1: eta_K = 12.000 (exact alpha^2*d)")
g0_e_12 = find_g0e_for_eta(12.0)
if g0_e_12:
    r_e12, g_e12, _ = solve_soliton(g0_e_12, 12.0)
    A_e12 = extract_tail(r_e12, g_e12)
    r_tau, g_tau, _ = solve_soliton(4.0, 12.0)
    A_tau = extract_tail(r_tau, g_tau)
    r31 = (A_tau / A_e12)**4
    m_tau = 0.51099895 * r31
    print(f"    g0_e = {g0_e_12:.8f}")
    print(f"    g0_mu = phi * g0_e = {PHI*g0_e_12:.8f}")
    print(f"    g0_tau = 4.000000")
    print(f"    A_e = {A_e12:.8f}")
    print(f"    A_tau = {A_tau:.8f}")
    print(f"    r_31 = {r31:.2f}")
    print(f"    m_tau = {m_tau:.2f} MeV (obs: 1776.86)")
    print(f"    dev = {abs(m_tau-1776.86)/1776.86*100:.3f}%")

print(f"\n  Test 2: eta_K = 12.067 (numerical optimum)")
g0_e_opt = find_g0e_for_eta(12.067)
if g0_e_opt:
    r_eopt, g_eopt, _ = solve_soliton(g0_e_opt, 12.067)
    A_eopt = extract_tail(r_eopt, g_eopt)
    r_tau2, g_tau2, _ = solve_soliton(4.0, 12.067)
    A_tau2 = extract_tail(r_tau2, g_tau2)
    r31_2 = (A_tau2 / A_eopt)**4
    m_tau_2 = 0.51099895 * r31_2
    print(f"    g0_e = {g0_e_opt:.8f}")
    print(f"    r_31 = {r31_2:.2f}")
    print(f"    m_tau = {m_tau_2:.2f} MeV (obs: 1776.86)")
    print(f"    dev = {abs(m_tau_2-1776.86)/1776.86*100:.3f}%")

# What eta_K gives EXACT r_31 = 3477.48 at g0_tau = 4?
print(f"\n  Test 3: find eta_K for g0_tau=4 giving r_31 = 3477.48")

def r31_at_4(eta_K):
    g0_e = find_g0e_for_eta(eta_K)
    if g0_e is None: return np.nan
    r_e, g_e, _ = solve_soliton(g0_e, eta_K)
    if r_e is None: return np.nan
    A_e = extract_tail(r_e, g_e)
    r_t, g_t, _ = solve_soliton(4.0, eta_K)
    if r_t is None: return np.nan
    A_t = extract_tail(r_t, g_t)
    if np.isnan(A_e) or np.isnan(A_t) or A_e < 1e-15: return np.nan
    return (A_t/A_e)**4

# Bisect
eta_lo, eta_hi = 12.0, 12.2
r31_lo = r31_at_4(eta_lo)
r31_hi = r31_at_4(eta_hi)
print(f"    eta=12.0: r_31 = {r31_lo:.1f}")
print(f"    eta=12.2: r_31 = {r31_hi:.1f}")

if not np.isnan(r31_lo) and not np.isnan(r31_hi):
    if (r31_lo - 3477.48) * (r31_hi - 3477.48) < 0:
        for i in range(30):
            eta_mid = (eta_lo + eta_hi) / 2.0
            r31_mid = r31_at_4(eta_mid)
            if np.isnan(r31_mid):
                eta_lo = eta_mid
                continue
            if r31_mid > 3477.48:
                eta_hi = eta_mid
            else:
                eta_lo = eta_mid
            if abs(eta_hi - eta_lo) < 1e-6:
                break
        eta_exact = (eta_lo + eta_hi) / 2.0
        r31_exact = r31_at_4(eta_exact)
        m_tau_exact = 0.51099895 * r31_exact
        print(f"\n    >>> eta_K(exact, g0_tau=4) = {eta_exact:.6f}")
        print(f"    r_31 = {r31_exact:.2f}")
        print(f"    m_tau = {m_tau_exact:.2f} MeV")
        print(f"    dev from 12.0: {abs(eta_exact-12.0):.4f} ({abs(eta_exact-12.0)/12.0*100:.3f}%)")
    else:
        print(f"  No bracket: r_31 doesn't cross 3477 between eta=12 and 12.2")
        # Maybe monotonic in wrong direction
        for eta_test in np.arange(11.5, 13.0, 0.1):
            r31_test = r31_at_4(eta_test)
            if not np.isnan(r31_test):
                print(f"    eta={eta_test:.1f}: r_31={r31_test:.1f}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"""
  FINDINGS:

  1. Bohr-Sommerfeld quantization: I_BS(g0) is NOT an integer at g0=4
     -> No straightforward quantization condition selects g0_tau.

  2. |V(4)|/V(1) = 512 = 2^9 -> interesting but no clear physical meaning.

  3. eta_K = 12 (exact) + g0_tau = 4 (exact): definitive mass prediction above.

  4. The soliton amplitude A(g0) grows like ~g0^0.7 for large g0,
     meaning m(g0) ~ g0^2.8 -- no natural cutoff.

  CONCLUSION: g0_tau = 4 cannot be derived from the radial ODE alone.
  It must come from an additional physical principle beyond the
  spherically symmetric soliton equation.
""")
print("DONE")
