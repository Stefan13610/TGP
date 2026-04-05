#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p145_FK_analytical_extremum.py -- Analytical derivation of F_K extremum at g0=2*alpha
=======================================================================================

GOAL: Prove that dF_K/dg0 = 0 at g0 = 2*alpha = 4, upgrading g0_tau from [HR] to [AN].

STRATEGY:
  The kinetic fraction F_K(g0) = E_kin / (E_kin + |E_pot|) has a maximum near g0=4.

  We use three complementary approaches:

  APPROACH A: VARIATIONAL ANSATZ
    Use a trial profile g(r) = 1 + (g0-1)*sech(r/R)^p with variational R(g0).
    Compute E_kin(g0,R) and E_pot(g0,R) semi-analytically.
    Show dF_K/dg0 = 0 at g0 = 2*alpha.

  APPROACH B: SCALING ANALYSIS
    Under rescaling g0 -> lambda*g0 (with appropriate R scaling),
    E_kin ~ f(g0)*(g0-1)^2 * R^(d-2)
    E_pot ~ V_eff(g0) * R^d
    The virial gives R^2 ~ f(g0)*(g0-1)^2 / V_eff(g0)
    Then F_K can be expressed as function of alpha_eff(g0) only.

  APPROACH C: EXACT ODE IDENTITY
    From the soliton ODE: f(g)*g'' + (2/r)*g' = V'(g)
    Multiply by g' and integrate -> energy identity.
    The ratio E_kin/E_pot involves integrals that can be
    evaluated at the extremum of F_K.

  APPROACH D: DIMENSIONLESS FORMULATION
    Define u = (g-1)/(g0-1), xi = r/R.
    The ODE becomes universal (g0-independent) up to corrections
    proportional to alpha_eff(g0). The g0-dependence enters only
    through alpha_eff(g0) and the potential shape.
    dF_K/dg0 = 0 when d(alpha_eff)/dg0 satisfies a specific condition.

Author: TGP project
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)
import numpy as np
from scipy.integrate import solve_ivp, trapezoid, quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gamma as Gamma

PHI = (1+np.sqrt(5))/2
ETA_K = 181.0/15
ALPHA = 2.0
D = 3

def alpha_eff(g, eta=ETA_K):
    return ALPHA/(1+eta*(g-1)**2)

def fk(g, eta=ETA_K):
    a = alpha_eff(g, eta)
    return 1+2*a*np.log(g) if g > 0 else -1e30

def V(g): return g**3/3 - g**4/4
def Vp(g): return g**2*(1-g)
def V_eff(g): return V(g) - V(1)

# ===================================================================
print("="*72)
print("  APPROACH A: VARIATIONAL ANSATZ FOR F_K EXTREMUM")
print("="*72)
# ===================================================================

print("""
  Trial profile: g(r) = 1 + (g0-1) * sech(r/R)^2

  In d=3 spherical:
    E_kin = 4*pi * int_0^inf 0.5*f(g)*g'^2 * r^2 dr
    E_pot = 4*pi * int_0^inf V_eff(g) * r^2 dr

  Key insight: f(g) varies with g(r), but in the SOLITON CORE
  where g ~ g0, f ~ f(g0). Near the tail where g ~ 1, f ~ 1.

  DOMINANT CONTRIBUTION: The core region dominates both E_kin and E_pot.
  So we approximate f(g(r)) ~ f_avg(g0) = weighted average.

  Better: split f(g) = 1 + 2*alpha_eff(g)*ln(g) into:
    f_1 = 1 (universal part)
    f_alpha = 2*alpha_eff(g)*ln(g) (coupling-dependent part)

  Then E_kin = E_kin^(1) + E_kin^(alpha), where E_kin^(alpha) depends on g0
  through alpha_eff.
""")

# First: compute exact F_K numerically for comparison
def solve_sol(g0, eta=ETA_K, rm=300):
    def fk_loc(g):
        a = ALPHA/(1+eta*(g-1)**2)
        return 1+2*a*np.log(g) if g>0 else -1e30
    fg0 = fk_loc(g0)
    if abs(fg0)<1e-15: return None,None,None
    c2 = Vp(g0)/(3*fg0)
    rs=0.01
    def rhs(r,y):
        g,p=y
        if g<=1e-15: return [p,0]
        fg=fk_loc(g)
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

def compute_energies(g0, eta=ETA_K):
    """Compute E_kin, E_pot, and their decomposition."""
    r, g, gp = solve_sol(g0, eta)
    if r is None: return None

    fk_arr = np.array([fk(gi, eta) for gi in g])
    ae_arr = np.array([alpha_eff(gi, eta) for gi in g])

    # Total kinetic
    E_kin = 4*np.pi * trapezoid(0.5*fk_arr*gp**2*r**2, r)
    # "Unit" kinetic (f=1 part)
    E_kin_1 = 4*np.pi * trapezoid(0.5*gp**2*r**2, r)
    # Alpha kinetic (f-1 part)
    f_alpha = 2*ae_arr*np.log(np.maximum(g, 1e-30))
    E_kin_a = 4*np.pi * trapezoid(0.5*f_alpha*gp**2*r**2, r)

    # Potential
    v_arr = V_eff(g)
    E_pot = 4*np.pi * trapezoid(v_arr*r**2, r)

    return {
        'E_kin': E_kin, 'E_pot': E_pot,
        'E_kin_1': E_kin_1, 'E_kin_a': E_kin_a,
        'F_K': E_kin/(E_kin+abs(E_pot)),
        'fk_core': fk(g0, eta),
        'alpha_core': alpha_eff(g0, eta),
    }

print(f"  {'g0':>6s} {'F_K':>10s} {'E_kin':>12s} {'E_kin_1':>12s} {'E_kin_a':>12s}"
      f" {'|E_pot|':>12s} {'fk(g0)':>8s} {'a_eff':>8s}")
print(f"  {'-'*90}")

FK_data = []
for g0 in np.arange(1.5, 8.01, 0.25):
    res = compute_energies(g0)
    if res is None: continue
    FK_data.append((g0, res['F_K'], res['E_kin'], res['E_kin_1'],
                     res['E_kin_a'], res['E_pot'], res['fk_core'], res['alpha_core']))
    mark = " <-- tau" if abs(g0-4.0)<0.13 else ""
    print(f"  {g0:6.2f} {res['F_K']:10.6f} {res['E_kin']:12.4f} {res['E_kin_1']:12.4f}"
          f" {res['E_kin_a']:12.4f} {abs(res['E_pot']):12.4f}"
          f" {res['fk_core']:8.4f} {res['alpha_core']:8.5f}{mark}")
sys.stdout.flush()

# Fine scan near g0=4
print(f"\n  Fine scan near g0=4:")
print(f"  {'g0':>8s} {'F_K':>12s} {'dF_K':>12s}")
print(f"  {'-'*35}")

FK_fine = []
prev_FK = None
for g0 in np.arange(3.0, 5.51, 0.1):
    res = compute_energies(g0)
    if res is None: continue
    dFK = res['F_K'] - prev_FK if prev_FK is not None else 0
    prev_FK = res['F_K']
    FK_fine.append((g0, res['F_K'], dFK))
    mark = ""
    if abs(g0-4.0)<0.05: mark = " <-- tau"
    if len(FK_fine)>1 and FK_fine[-2][2]>0 and dFK<0:
        mark += " ** sign change **"
    print(f"  {g0:8.2f} {res['F_K']:12.8f} {dFK:12.8f}{mark}")
sys.stdout.flush()

# Precise sign change via bisection
print(f"\n  Bisecting dF_K/dg0 = 0 ...")
def FK_at(g0):
    res = compute_energies(g0)
    return res['F_K'] if res else np.nan

# Use finite differences
def dFK_dg0(g0, h=0.02):
    return (FK_at(g0+h) - FK_at(g0-h))/(2*h)

# Find sign change
g_lo, g_hi = 3.5, 5.0
for _ in range(40):
    g_mid = (g_lo+g_hi)/2
    if dFK_dg0(g_mid) > 0:
        g_lo = g_mid
    else:
        g_hi = g_mid
g0_star = (g_lo+g_hi)/2
print(f"  F_K maximum at g0* = {g0_star:.6f}")
print(f"  Deviation from 2*alpha=4: {g0_star-4.0:.6f} ({(g0_star-4.0)/4.0*100:.4f}%)")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH B: SCALING / DIMENSIONAL ANALYSIS")
print("="*72)
# ===================================================================

print("""
  For the soliton ODE: f(g)*g'' + (2/r)*g' = V'(g)

  Define dimensionless variables:
    u = (g-1)/(g0-1),  xi = r/R
  where R is the effective soliton radius.

  The virial relation (from Derrick scaling) gives:
    (d-2)*E_kin + d*E_pot = surface terms (from oscillating tail)

  For TGP solitons with oscillating tail, the standard virial
  is violated. But we can use the LOCAL virial at the core.

  KEY OBSERVATION: In the core (r << R_tail), the profile is
  monotonically decreasing and the standard virial holds approximately:
    E_kin^core ~ -(d/(d-2)) * E_pot^core  (for d=3: E_kin ~ -3*E_pot)

  The TAIL contribution (r >> R_core) adds:
    E_kin^tail ~ A^2 * int r^2 * sin^2(r)/r^2 dr ~ A^2 * R_tail
    E_pot^tail ~ A^2 * int r^2 * (g-1)^2/r^2 dr  ~ A^2 * R_tail

  So F_K = E_kin/(E_kin+|E_pot|) depends on the core/tail PARTITION.

  The partition depends on alpha_eff(g0):
  - Large alpha_eff -> strong kinetic term -> extended core -> more E_kin
  - Small alpha_eff -> weak kinetic term -> compact core -> more E_pot

  alpha_eff(g0) = 2/(1 + eta_K*(g0-1)^2)

  So alpha_eff is a DECREASING function of g0 (for g0 > 1).
  F_K should INCREASE with alpha_eff (more kinetic coupling = more E_kin).
  But f(g0) = 1 + 2*alpha_eff*ln(g0) ALSO matters: for large g0,
  even though alpha_eff is small, ln(g0) grows, partially compensating.
""")

# Compute the "alpha_eff * ln(g0)" product
print(f"  {'g0':>6s} {'alpha_eff':>10s} {'ln(g0)':>8s} {'a*ln(g0)':>10s} {'f(g0)':>8s}"
      f" {'f(g0)-1':>10s}")
print(f"  {'-'*60}")

for g0 in np.arange(1.5, 8.1, 0.5):
    ae = alpha_eff(g0)
    lg = np.log(g0)
    fg = fk(g0)
    mark = " <-- tau" if abs(g0-4.0)<0.1 else ""
    print(f"  {g0:6.2f} {ae:10.5f} {lg:8.4f} {ae*lg:10.5f} {fg:8.4f} {fg-1:10.5f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH C: ANALYTICAL F_K FROM f(g0) STRUCTURE")
print("="*72)
# ===================================================================

print("""
  KEY INSIGHT: F_K depends on g0 primarily through the RATIO:

    rho(g0) = [f(g0) - 1] / f(g0) = 2*alpha_eff(g0)*ln(g0) / f(g0)

  This is the fraction of f(g0) coming from the coupling term.

  When rho = 0 (f = 1): pure "free field" kinetic term
  When rho -> 1 (f >> 1): kinetic term dominated by coupling

  HYPOTHESIS: F_K ~ F_0 + C * rho(g0) where F_0, C are g0-independent.
  Then dF_K/dg0 = 0 when d(rho)/dg0 = 0.

  Let's compute rho(g0) and check if its maximum is at g0 = 2*alpha.
""")

def rho(g0, eta=ETA_K):
    fg = fk(g0, eta)
    return (fg - 1)/fg if fg > 0 else 0

print(f"  {'g0':>6s} {'rho':>12s} {'d_rho':>12s}")
print(f"  {'-'*35}")

rho_prev = None
rho_data = []
for g0 in np.arange(1.2, 8.01, 0.1):
    r_val = rho(g0)
    dr = r_val - rho_prev if rho_prev is not None else 0
    rho_prev = r_val
    rho_data.append((g0, r_val, dr))
    if int(g0*10) % 5 == 0 or abs(g0-4.0)<0.05:
        mark = " <-- tau" if abs(g0-4.0)<0.05 else ""
        if len(rho_data)>1 and rho_data[-2][2]>0 and dr<0:
            mark += " ** max **"
        print(f"  {g0:6.2f} {r_val:12.8f} {dr:12.8f}{mark}")

# Find rho maximum precisely
def rho_neg(g0):
    return -rho(g0)
res = minimize_scalar(rho_neg, bounds=(1.5, 10.0), method='bounded')
g0_rho_max = res.x
print(f"\n  Maximum of rho at g0 = {g0_rho_max:.6f}")
print(f"  rho_max = {rho(g0_rho_max):.8f}")
print(f"  Deviation from 2*alpha=4: {g0_rho_max-4.0:.6f} ({(g0_rho_max-4.0)/4.0*100:.3f}%)")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH C2: ANALYTICAL DERIVATIVE d(rho)/dg0 = 0")
print("="*72)
# ===================================================================

print("""
  rho(g0) = 2*a(g0)*ln(g0) / (1 + 2*a(g0)*ln(g0))
  where a(g0) = alpha/(1 + eta*(g0-1)^2)

  Define u = 2*a(g0)*ln(g0) = kinetic coupling * field excursion.
  Then rho = u/(1+u), so d(rho)/dg0 = (du/dg0)/(1+u)^2.

  dF_K/dg0 = 0 <=> du/dg0 = 0  (since (1+u)^2 > 0)

  u(g0) = 2*alpha*ln(g0) / (1 + eta*(g0-1)^2)

  du/dg0 = 2*alpha * [ 1/g0 * (1+eta*(g0-1)^2) - ln(g0)*2*eta*(g0-1) ]
                     / (1+eta*(g0-1)^2)^2

  Setting du/dg0 = 0:
    1/g0 * (1+eta*(g0-1)^2) = 2*eta*(g0-1)*ln(g0)

  Define x = g0-1 (excursion from vacuum). Then g0 = 1+x, and:
    (1+eta*x^2) / (1+x) = 2*eta*x*ln(1+x)

  For alpha = 2, eta = 181/15:
""")

eta = ETA_K

def du_condition(g0):
    """Returns LHS - RHS of the extremum condition."""
    x = g0 - 1
    if x <= 0: return 1e10
    lhs = (1 + eta*x**2) / (1+x)
    rhs = 2*eta*x*np.log(1+x)
    return lhs - rhs

# Find zero
print(f"  {'g0':>8s} {'LHS':>12s} {'RHS':>12s} {'LHS-RHS':>12s}")
print(f"  {'-'*48}")
for g0 in np.arange(1.5, 8.1, 0.25):
    x = g0-1
    lhs = (1+eta*x**2)/(1+x)
    rhs = 2*eta*x*np.log(1+x)
    mark = " <-- tau" if abs(g0-4.0)<0.13 else ""
    print(f"  {g0:8.3f} {lhs:12.4f} {rhs:12.4f} {lhs-rhs:12.4f}{mark}")

g0_u_max = brentq(du_condition, 1.5, 10.0)
print(f"\n  Analytical extremum of u(g0) at g0 = {g0_u_max:.8f}")
print(f"  Deviation from 4: {g0_u_max-4.0:.6f}")
print(f"  Deviation from 2*alpha: {g0_u_max-2*ALPHA:.6f}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH C3: UNIVERSAL FORM OF THE EXTREMUM CONDITION")
print("="*72)
# ===================================================================

print("""
  The condition du/dg0 = 0 is:
    (1 + eta*x^2)/(1+x) = 2*eta*x*ln(1+x)    where x = g0-1

  Substitute eta = alpha^2*d + 1/((alpha^2+1)*d) = alpha^2*d*(1 + epsilon)
  where epsilon = 1/((alpha^2+1)*d*alpha^2*d) << 1.

  At LEADING ORDER (eta ~ alpha^2*d):
    (1 + alpha^2*d*x^2)/(1+x) = 2*alpha^2*d*x*ln(1+x)

  For large eta*x^2 >> 1 (which holds for x=3, eta~12):
    eta*x^2/(1+x) ~ 2*eta*x*ln(1+x)
    x/(1+x) ~ 2*ln(1+x)

  Define t = 1+x = g0. Then x = t-1:
    (t-1)/t = 2*ln(t)
    1 - 1/t = 2*ln(t)

  This is a UNIVERSAL equation (independent of alpha, eta, d)!
  Let's solve it.
""")

def universal_eq(t):
    """1 - 1/t - 2*ln(t) = 0"""
    return 1 - 1/t - 2*np.log(t)

# Scan
print(f"  {'t=g0':>8s} {'1-1/t':>10s} {'2*ln(t)':>10s} {'diff':>12s}")
print(f"  {'-'*44}")
for t in np.arange(1.1, 6.1, 0.2):
    d = universal_eq(t)
    mark = ""
    if abs(t-4.0)<0.1: mark = " <-- tau"
    print(f"  {t:8.2f} {1-1/t:10.6f} {2*np.log(t):10.6f} {d:12.6f}{mark}")

# The universal equation 1-1/t = 2*ln(t) has solution near t ~ 1.76
# NOT at t = 4. So the leading-order approximation misses it.

# Let's be more careful. The FULL condition is:
# (1 + eta*x^2)/(1+x) = 2*eta*x*ln(1+x)
# Don't drop the "1" term.

print(f"\n  The universal equation 1-1/t = 2*ln(t) has NO solution at t=4.")
print(f"  (Solution is near t ~ 1.76)")
print(f"  So the leading-order approximation is TOO crude.")
print(f"\n  The FULL condition gives g0 = {g0_u_max:.6f}, NOT universal.")
print(f"  This means the F_K extremum position depends on eta (and thus on alpha, d).")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH D: PARAMETRIC STUDY — eta DEPENDENCE OF g0*")
print("="*72)
# ===================================================================

print("""
  If g0*(eta) is the F_K extremum, does g0* = 2*alpha require
  eta = alpha^2*d specifically?

  The condition du/dg0 = 0 depends on eta. Let's compute g0*(eta)
  for a range of eta values and check if g0* = 4 occurs at eta = 181/15.
""")

print(f"  {'eta':>8s} {'g0*':>10s} {'g0*/2alpha':>12s} {'note':>20s}")
print(f"  {'-'*55}")

for eta_test in [1, 2, 4, 6, 8, 10, 11, 12, 12.0667, 13, 14, 16, 20, 30, 50]:
    def cond(g0, et=eta_test):
        x = g0-1
        if x <= 0: return 1e10
        return (1+et*x**2)/(1+x) - 2*et*x*np.log(1+x)
    try:
        g0s = brentq(cond, 1.01, 50.0)
    except:
        g0s = np.nan
    ratio = g0s/(2*ALPHA)
    note = ""
    if abs(eta_test-ETA_K)<0.01: note = "<-- TGP"
    if abs(g0s-4.0)<0.05: note += " g0*~4!"
    print(f"  {eta_test:8.4f} {g0s:10.6f} {ratio:12.6f}  {note}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH E: EXACT SELF-CONSISTENCY CONDITION")
print("="*72)
# ===================================================================

print("""
  From Approach D, we see g0*(eta) is a smooth function.
  The question is: at eta = 181/15, does g0* = 4 EXACTLY?

  The condition is:
    (1 + eta*x^2)/(1+x) = 2*eta*x*ln(1+x),  x = g0-1 = 3

  With eta = 181/15, x = 3:
    LHS = (1 + 181/15 * 9) / 4 = (1 + 1629/15) / 4 = (15/15 + 1629/15) / 4
        = 1644/(15*4) = 1644/60 = 27.4
    RHS = 2 * 181/15 * 3 * ln(4) = 1086/15 * ln(4) = 72.4 * 1.38629...
        = 72.4 * 1.38629 = 100.37...

  These are NOT equal! So g0=4 is NOT the rho-maximum for eta=181/15.

  Let me check what g0* actually is...
""")

x = 3.0  # g0 = 4
lhs = (1 + ETA_K*x**2)/(1+x)
rhs = 2*ETA_K*x*np.log(1+x)
print(f"  At g0=4 (x=3): LHS = {lhs:.6f}, RHS = {rhs:.6f}")
print(f"  Ratio LHS/RHS = {lhs/rhs:.6f}")
print(f"  g0*(rho) = {g0_u_max:.6f} (NOT 4.0)")

print(f"""
  CONCLUSION: The rho(g0) = [f(g0)-1]/f(g0) maximum is at g0 = {g0_u_max:.4f},
  NOT at g0 = 4.0. So rho is NOT the right proxy for F_K.

  But the NUMERICAL F_K has its maximum near g0 = {g0_star:.4f}.
  The discrepancy means that F_K depends on MORE than just rho(g0).

  The profile shape (how g(r) depends on g0) also matters.
  Let's investigate the FULL numerical F_K more carefully.
""")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH F: DECOMPOSE F_K INTO CORE AND TAIL")
print("="*72)
# ===================================================================

print("""
  Split the radial integral at r_cut = R_core (where g crosses 1+epsilon).
  F_K = F_K^core + F_K^tail.

  The CORE part should follow the rho analysis.
  The TAIL part adds a correction that shifts the maximum.
""")

def compute_FK_decomposed(g0, eta=ETA_K, r_cut=None):
    r, g, gp = solve_sol(g0, eta)
    if r is None: return None

    # Find r_cut: where g first crosses 1.01 (heading down from g0)
    if r_cut is None:
        cross_idx = np.where(np.diff(np.sign(g - 1.01)))[0]
        if len(cross_idx) > 0:
            r_cut = r[cross_idx[0]]
        else:
            r_cut = r[-1]*0.3

    m_core = r <= r_cut
    m_tail = r > r_cut

    fk_arr = np.array([fk(gi, eta) for gi in g])
    v_arr = V_eff(g)

    ek_core = 4*np.pi * trapezoid(0.5*fk_arr[m_core]*gp[m_core]**2*r[m_core]**2, r[m_core]) if m_core.sum()>2 else 0
    ek_tail = 4*np.pi * trapezoid(0.5*fk_arr[m_tail]*gp[m_tail]**2*r[m_tail]**2, r[m_tail]) if m_tail.sum()>2 else 0
    ep_core = 4*np.pi * trapezoid(v_arr[m_core]*r[m_core]**2, r[m_core]) if m_core.sum()>2 else 0
    ep_tail = 4*np.pi * trapezoid(v_arr[m_tail]*r[m_tail]**2, r[m_tail]) if m_tail.sum()>2 else 0

    E_kin = ek_core + ek_tail
    E_pot = ep_core + ep_tail

    return {
        'E_kin': E_kin, 'E_pot': E_pot, 'r_cut': r_cut,
        'ek_core': ek_core, 'ek_tail': ek_tail,
        'ep_core': ep_core, 'ep_tail': ep_tail,
        'F_K': E_kin/(E_kin+abs(E_pot)),
        'frac_core_kin': ek_core/E_kin if E_kin>0 else 0,
        'frac_core_pot': abs(ep_core)/abs(E_pot) if abs(E_pot)>0 else 0,
    }

print(f"  {'g0':>6s} {'F_K':>10s} {'EkCore%':>8s} {'EpCore%':>8s} "
      f"{'Rcut':>6s} {'F_K_core':>10s}")
print(f"  {'-'*55}")

for g0 in np.arange(2.0, 7.1, 0.25):
    res = compute_FK_decomposed(g0)
    if res is None: continue
    FK_core = res['ek_core']/(res['ek_core']+abs(res['ep_core'])) if (res['ek_core']+abs(res['ep_core']))>0 else 0
    mark = " <-- tau" if abs(g0-4.0)<0.13 else ""
    print(f"  {g0:6.2f} {res['F_K']:10.6f} {res['frac_core_kin']*100:8.1f} "
          f"{res['frac_core_pot']*100:8.1f} {res['r_cut']:6.2f} {FK_core:10.6f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH G: THE f*V' BALANCE CONDITION")
print("="*72)
# ===================================================================

print("""
  New idea: at the SOLITON CENTER (r=0), the ODE gives:
    g''(0) = V'(g0) / (3*f(g0))

  The "driving force per unit kinetic inertia" is:
    Sigma(g0) = V'(g0) / f(g0) = g0^2*(1-g0) / (1 + 2*alpha_eff(g0)*ln(g0))

  This is the EFFECTIVE CURVATURE at the soliton center.

  Define the NORMALIZED driving force:
    sigma_n(g0) = |Sigma(g0)| / (g0-1)^2

  This measures how efficiently the potential drives the soliton
  per unit field excursion squared.

  Also define the KINETIC EFFICIENCY:
    kappa(g0) = f(g0) * (g0-1) / (g0^2 * |g0-1|)
              = f(g0) / g0^2  (for g0 > 1)

  The soliton energy F_K should be related to kappa.
""")

print(f"  {'g0':>6s} {'|Sigma|':>12s} {'sigma_n':>12s} {'kappa':>12s} {'kappa/max':>10s}")
print(f"  {'-'*55}")

kappas = []
for g0 in np.arange(1.5, 8.1, 0.25):
    fg = fk(g0)
    vp = Vp(g0)
    sigma = abs(vp/fg)
    sigma_n = sigma/(g0-1)**2
    kappa = fg/g0**2
    kappas.append((g0, kappa))
    mark = " <-- tau" if abs(g0-4.0)<0.13 else ""
    print(f"  {g0:6.2f} {sigma:12.6f} {sigma_n:12.6f} {kappa:12.6f}{mark}")

kappas = np.array(kappas)
max_idx = np.argmax(kappas[:,1])
print(f"\n  Maximum of kappa = f(g0)/g0^2 at g0 = {kappas[max_idx,0]:.4f}")
print(f"  Deviation from 4: {kappas[max_idx,0]-4.0:.4f}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH H: d/dg0 [f(g0)/g0^2] = 0  ANALYTICALLY")
print("="*72)
# ===================================================================

print("""
  kappa(g0) = f(g0)/g0^2 = [1 + 2*a(g0)*ln(g0)] / g0^2

  d(kappa)/dg0 = 0 gives:
    g0 * f'(g0) = 2 * f(g0)

  where f'(g0) = d/dg0 [1 + 2*a(g0)*ln(g0)]
               = 2*[a'(g0)*ln(g0) + a(g0)/g0]
  and a'(g0) = -2*alpha*eta*(g0-1) / (1+eta*(g0-1)^2)^2

  So the condition is:
    g0 * 2 * [a'*ln(g0) + a/g0] = 2*(1 + 2*a*ln(g0))
    g0*a'*ln(g0) + a = 1 + 2*a*ln(g0)

  Substituting a = alpha/(1+eta*x^2), a' = -2*alpha*eta*x/(1+eta*x^2)^2
  where x = g0-1:
""")

def kappa_condition(g0, eta_val=ETA_K):
    """Returns the value of [g0*f'(g0) - 2*f(g0)], should be 0 at extremum."""
    x = g0 - 1
    s = 1 + eta_val*x**2
    a = ALPHA/s
    ap = -2*ALPHA*eta_val*x/s**2
    lg = np.log(g0)

    f_val = 1 + 2*a*lg
    fp_val = 2*(ap*lg + a/g0)

    return g0*fp_val - 2*f_val

# Find zero
g0_kappa = brentq(kappa_condition, 1.5, 10.0)
print(f"  kappa extremum at g0 = {g0_kappa:.8f}")
print(f"  Deviation from 4: {g0_kappa-4.0:.6f} ({(g0_kappa-4.0)/4.0*100:.4f}%)")

# Check at g0 = 4 exactly
val_at_4 = kappa_condition(4.0)
print(f"  kappa_condition(4.0) = {val_at_4:.8f} (should be 0 if exact)")

# What eta would make kappa maximum at exactly g0=4?
def kappa_at_4(eta_test):
    return kappa_condition(4.0, eta_test)

try:
    eta_exact = brentq(kappa_at_4, 1.0, 100.0)
    print(f"\n  eta that puts kappa max at g0=4: eta = {eta_exact:.8f}")
    print(f"  Compare: 181/15 = {181/15:.8f}")
    print(f"  Ratio: {eta_exact/(181/15):.6f}")
    # Express as fraction
    from fractions import Fraction
    frac = Fraction(eta_exact).limit_denominator(1000)
    print(f"  Nearest simple fraction: {frac} = {float(frac):.8f}")
except:
    print(f"  Could not find eta for kappa max at g0=4")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH I: SELF-CONSISTENT CONDITION g0 = 2*alpha AT KAPPA MAX")
print("="*72)
# ===================================================================

print("""
  KEY QUESTION: Is there a STRUCTURAL reason why the kappa maximum
  should be at g0 = 2*alpha, independent of the specific eta value?

  Let's substitute g0 = 2*alpha into the kappa condition:
    g0*f'(g0) = 2*f(g0)

  At g0 = 2*alpha, with a(g0) = alpha/(1+eta*(2*alpha-1)^2):

  For alpha = 2:
    g0 = 4, x = 3
    a(4) = 2/(1+9*eta)
    a'(4) = -2*2*eta*3/(1+9*eta)^2 = -12*eta/(1+9*eta)^2
    ln(4) = 2*ln(2)
    f(4) = 1 + 2*a*2*ln(2) = 1 + 4*ln(2)*2/(1+9*eta)
    f'(4) = 2*[a'*2*ln(2) + a/4]

  Condition: 4*f'(4) = 2*f(4)
    4*2*[-12*eta/(1+9*eta)^2 * 2*ln(2) + 2/(1+9*eta)/4] = 2*[1 + 8*ln(2)/(1+9*eta)]
    8*[-24*eta*ln(2)/(1+9*eta)^2 + 1/(2*(1+9*eta))] = 2 + 16*ln(2)/(1+9*eta)

  Let s = 1+9*eta:
    8*[-24*(s-1)/9*ln(2)/s^2 + 1/(2*s)] = 2 + 16*ln(2)/s
    8*[-24*(s-1)*ln(2)/(9*s^2) + 1/(2*s)] = 2 + 16*ln(2)/s

  This is getting messy. Let me just verify numerically across alpha values.
""")

# Generalize: for any alpha, what eta makes kappa max at g0 = 2*alpha?
print(f"  {'alpha':>6s} {'g0=2a':>6s} {'eta*':>12s} {'eta*/a^2d':>10s} {'note':>20s}")
print(f"  {'-'*60}")

for alpha_test in [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
    def kappa_cond_gen(g0, eta_val, al=alpha_test):
        x = g0-1
        s = 1 + eta_val*x**2
        a = al/s
        ap = -2*al*eta_val*x/s**2
        lg = np.log(g0)
        f_val = 1+2*a*lg
        fp_val = 2*(ap*lg + a/g0)
        return g0*fp_val - 2*f_val

    def find_eta(eta_val, al=alpha_test):
        g0_target = 2*al
        return kappa_cond_gen(g0_target, eta_val, al)

    try:
        eta_s = brentq(find_eta, 0.01, 1000.0)
        ratio = eta_s / (alpha_test**2 * D)
        note = ""
        if abs(alpha_test-2.0)<0.01: note = f"<-- TGP (181/15={181/15:.4f})"
        print(f"  {alpha_test:6.2f} {2*alpha_test:6.1f} {eta_s:12.6f} {ratio:10.6f}  {note}")
    except:
        print(f"  {alpha_test:6.2f} {2*alpha_test:6.1f} {'N/A':>12s}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  APPROACH J: DIRECT PROOF — F_K vs kappa CORRELATION")
print("="*72)
# ===================================================================

print("""
  Check: is the kappa = f(g0)/g0^2 maximum a good proxy for the
  ACTUAL F_K maximum (computed from full ODE solutions)?
""")

# Compare g0*(kappa) vs g0*(F_K) for various eta
print(f"  {'eta':>8s} {'g0*(kappa)':>12s} {'g0*(F_K)':>12s} {'diff':>10s}")
print(f"  {'-'*48}")

for eta_test in [4.0, 8.0, 10.0, 12.0, 12.0667, 14.0, 18.0, 25.0]:
    # kappa maximum
    try:
        g0_k = brentq(lambda g: kappa_condition(g, eta_test), 1.5, 15.0)
    except:
        g0_k = np.nan

    # F_K maximum (numerical)
    def neg_FK(g0, et=eta_test):
        r, g, gp = solve_sol(g0, et)
        if r is None: return 0
        fk_arr = np.array([1+2*(ALPHA/(1+et*(gi-1)**2))*np.log(gi) if gi>0 else 1 for gi in g])
        v_arr = V_eff(g)
        ek = 4*np.pi*trapezoid(0.5*fk_arr*gp**2*r**2, r)
        ep = 4*np.pi*trapezoid(v_arr*r**2, r)
        return -(ek/(ek+abs(ep)))

    try:
        res = minimize_scalar(neg_FK, bounds=(1.5, 10.0), method='bounded',
                             args=(eta_test,), options={'xatol': 0.01})
        g0_f = res.x
    except:
        g0_f = np.nan

    diff = g0_k - g0_f if not (np.isnan(g0_k) or np.isnan(g0_f)) else np.nan
    note = " <-- TGP" if abs(eta_test-ETA_K)<0.01 else ""
    print(f"  {eta_test:8.4f} {g0_k:12.6f} {g0_f:12.4f} {diff:10.4f}{note}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  SYNTHESIS AND CONCLUSIONS")
print("="*72)
# ===================================================================

print(f"""
  RESULTS SUMMARY:

  1. F_K (numerical, full ODE) maximum at g0 = {g0_star:.4f}
  2. rho = (f-1)/f maximum at g0 = {g0_u_max:.4f}
  3. kappa = f/g0^2 maximum at g0 = {g0_kappa:.4f}

  KEY FINDING: kappa = f(g0)/g0^2 is an ANALYTICAL proxy for F_K.
  Its extremum condition d(kappa)/dg0 = 0 gives:

    g0*f'(g0) = 2*f(g0)

  which is the condition that the kinetic function grows at the
  SAME RATE as g0^2 (the potential's natural scale).

  This condition depends on eta. For eta = 181/15 (TGP):
    g0*(kappa) = {g0_kappa:.6f}

  INTERPRETATION:
  - The kappa proxy gives g0 ~ {g0_kappa:.2f}, which is CLOSE to 4
    but NOT exactly 4.
  - The deviation ({g0_kappa-4.0:.4f}) may be within the approximation
    of using kappa as F_K proxy.
  - The K-self-duality condition g0 = 2*alpha = 4 remains the CLEANEST
    algebraic condition, with F_K/kappa extremum as SUPPORTING EVIDENCE.

  STATUS ASSESSMENT:
  - g0_tau = 2*alpha = 4 from K-self-duality: ALGEBRAIC [HR]
  - F_K extremum near g0 ~ 4: NUMERICAL SUPPORT
  - kappa = f/g0^2 extremum near g0 ~ {g0_kappa:.2f}: ANALYTICAL SUPPORT
  - Full analytical proof that g0* = 4 EXACTLY: NOT YET ACHIEVED

  The most promising path to [AN] closure:
  Show that the K-self-duality condition K(g0) = g0^g0 is equivalent
  to a VARIATIONAL PRINCIPLE for the soliton action.
""")

print("DONE")
