#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p145b_FK_kappa_analysis.py -- Precision analysis of F_K and kappa extrema
==========================================================================

KEY CORRECTION from p145:
  F_K maximum is at g0 = 3.845, NOT 4.0.
  The previous claim (p144) was based on coarse sampling.

GOAL: Find what analytical quantity (if any) peaks at EXACTLY g0 = 4.0
for eta_K = 181/15, and understand the relationship.

APPROACH:
  1. Precise F_K(g0) scan and its maximum
  2. kappa(g0) = f(g0)/g0^2 and its maximum
  3. Other analytical proxies
  4. eta_K dependence: at what eta does F_K peak at 4.0?
  5. Self-consistent condition analysis

Author: TGP project
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import brentq, minimize_scalar

PHI = (1+np.sqrt(5))/2
ETA_K = 181.0/15
ALPHA = 2.0
D = 3

def alpha_eff(g, eta=ETA_K):
    return ALPHA/(1+eta*(g-1)**2)

def fk(g, eta=ETA_K):
    a = alpha_eff(g, eta)
    return 1+2*a*np.log(g) if g > 0 else -1e30

def Vp(g): return g**2*(1-g)
def V(g): return g**3/3 - g**4/4
def V_eff(g): return V(g) - V(1)

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

def compute_FK(g0, eta=ETA_K):
    r, g, gp = solve_sol(g0, eta)
    if r is None: return np.nan
    fk_arr = np.array([fk(gi, eta) for gi in g])
    v_arr = V_eff(g)
    E_kin = 4*np.pi * trapezoid(0.5*fk_arr*gp**2*r**2, r)
    E_pot = 4*np.pi * trapezoid(v_arr*r**2, r)
    return E_kin / (E_kin + abs(E_pot)) if (E_kin + abs(E_pot)) > 0 else np.nan

# ===================================================================
print("="*72)
print("  PART 1: PRECISE F_K MAXIMUM")
print("="*72)
# ===================================================================

# Bisect dF_K/dg0 = 0 with fine steps
def dFK(g0, eta=ETA_K, h=0.005):
    return (compute_FK(g0+h, eta) - compute_FK(g0-h, eta))/(2*h)

g_lo, g_hi = 2.5, 5.5
for _ in range(50):
    g_mid = (g_lo+g_hi)/2
    if dFK(g_mid) > 0:
        g_lo = g_mid
    else:
        g_hi = g_mid

g0_FK_max = (g_lo+g_hi)/2
FK_max = compute_FK(g0_FK_max)
FK_at_4 = compute_FK(4.0)

print(f"\n  F_K maximum: g0* = {g0_FK_max:.6f}")
print(f"  F_K(g0*)    = {FK_max:.10f}")
print(f"  F_K(4.0)    = {FK_at_4:.10f}")
print(f"  F_K(4)/F_K* = {FK_at_4/FK_max:.8f}")
print(f"  Deviation from 4: {g0_FK_max-4.0:.6f} ({(g0_FK_max-4.0)/4.0*100:.3f}%)")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  PART 2: kappa = f(g0)/g0^n FOR VARIOUS n")
print("="*72)
# ===================================================================

print("""
  kappa_n(g0) = f(g0)/g0^n.
  Find n such that kappa_n maximum is at g0 = 4.0.
""")

def kappa_n_max(n, eta=ETA_K):
    """Find g0 that maximizes f(g0)/g0^n."""
    def neg_kn(g0):
        return -fk(g0, eta)/g0**n
    res = minimize_scalar(neg_kn, bounds=(1.1, 15.0), method='bounded')
    return res.x

# Find n that gives g0* = 4
print(f"  {'n':>6s} {'g0*(kappa_n)':>14s} {'dev from 4':>12s}")
print(f"  {'-'*35}")
for n in np.arange(0.5, 4.1, 0.25):
    g0s = kappa_n_max(n)
    print(f"  {n:6.2f} {g0s:14.6f} {g0s-4.0:12.4f}")

# Bisect for n that gives g0=4
def cond_n(n):
    return kappa_n_max(n) - 4.0

try:
    n_star = brentq(cond_n, 0.1, 4.0)
    print(f"\n  n* such that kappa_n max at g0=4: n* = {n_star:.8f}")
    # Simple fraction?
    from fractions import Fraction
    frac = Fraction(n_star).limit_denominator(100)
    print(f"  Nearest fraction: {frac} = {float(frac):.8f}")
    print(f"  Deviation: {n_star - float(frac):.6e}")
except Exception as e:
    print(f"  Error: {e}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  PART 3: ETA_K DEPENDENCE OF F_K MAXIMUM")
print("="*72)
# ===================================================================

print("""
  At what eta_K does the F_K maximum land on g0 = 4.0?
""")

print(f"  {'eta_K':>10s} {'g0*(F_K)':>12s} {'dev from 4':>12s}")
print(f"  {'-'*40}")

for eta_test in [0, 2, 4, 6, 8, 10, 11, 12, 12.0667, 13, 14, 16, 20, 30]:
    # Bisect
    try:
        g_lo, g_hi = 1.5, 10.0
        for _ in range(40):
            g_mid = (g_lo+g_hi)/2
            if dFK(g_mid, eta_test) > 0:
                g_lo = g_mid
            else:
                g_hi = g_mid
        g0s = (g_lo+g_hi)/2
    except:
        g0s = np.nan
    note = " <-- TGP" if abs(eta_test-ETA_K)<0.01 else ""
    if abs(g0s-4.0)<0.05: note += " g0~4!"
    print(f"  {eta_test:10.4f} {g0s:12.6f} {g0s-4.0:12.4f}{note}")
sys.stdout.flush()

# Find eta that gives F_K max at g0=4
def FK_max_at_4(eta_test):
    g_lo, g_hi = 2.0, 8.0
    for _ in range(40):
        g_mid = (g_lo+g_hi)/2
        if dFK(g_mid, eta_test) > 0:
            g_lo = g_mid
        else:
            g_hi = g_mid
    return (g_lo+g_hi)/2 - 4.0

# Identify the crossing
print(f"\n  Searching eta that puts F_K max at g0 = 4.0 ...")
try:
    eta_for_4 = brentq(FK_max_at_4, 0.1, 100.0)
    print(f"  eta*(F_K max at 4) = {eta_for_4:.6f}")
    print(f"  Compare eta_K(TGP) = {ETA_K:.6f}")
    print(f"  Ratio: {eta_for_4/ETA_K:.6f}")
except Exception as e:
    print(f"  No crossing found: {e}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  PART 4: ANALYTICAL PROXIES THAT MIGHT PEAK AT g0=4")
print("="*72)
# ===================================================================

print("""
  Test various analytical functions Q(g0) for maximum at g0=4:
""")

proxies = {}

# Proxy 1: f(g0) * (g0-1)^2 / g0^4  (kinetic weight * excursion / K)
def proxy1(g0):
    return fk(g0)*(g0-1)**2/g0**4

# Proxy 2: alpha_eff(g0) * |V'(g0)| * g0
def proxy2(g0):
    return alpha_eff(g0)*abs(Vp(g0))*g0

# Proxy 3: f(g0) * |V'(g0)| / g0^3  (dimensionless driving)
def proxy3(g0):
    return fk(g0)*abs(Vp(g0))/g0**3

# Proxy 4: (f(g0)-1) * (g0-1)  (coupling * excursion)
def proxy4(g0):
    return (fk(g0)-1)*(g0-1)

# Proxy 5: alpha_eff * ln(g0) * (g0-1)
def proxy5(g0):
    return alpha_eff(g0)*np.log(g0)*(g0-1)

# Proxy 6: f(g0) / (g0 * ln(g0))
def proxy6(g0):
    lg = np.log(g0)
    return fk(g0)/(g0*lg) if lg > 0.01 else np.nan

# Proxy 7: (f-1)/f * (g0-1)/g0 = rho * (1-1/g0)
def proxy7(g0):
    fg = fk(g0)
    return (fg-1)/fg * (1-1/g0)

# Proxy 8: E_kin_alpha / E_kin_total (coupling fraction of kinetic energy)
# This is numerical but might reveal structure
# Skip for now, use analytical approximation:
# E_kin_a/E_kin ~ (f-1)/f at core
def proxy8(g0):
    fg = fk(g0)
    # Weight by (g0-1)^2 to account for excursion amplitude
    return (fg-1)/fg * (g0-1)**2

# Proxy 9: d[ln f]/d[ln g] at g=g0 (logarithmic derivative)
def proxy9(g0):
    a = alpha_eff(g0)
    fg = fk(g0)
    # d(ln f)/d(ln g) = g * f'/f
    # f = 1+2*a*ln(g), f' = 2*(a'*ln(g)+a/g)
    x = g0-1
    s = 1+ETA_K*x**2
    ap = -2*ALPHA*ETA_K*x/s**2
    fp = 2*(ap*np.log(g0)+a/g0)
    return g0*fp/fg

# Proxy 10: Ratio V'(g0) / (g0 * V'(g0/2))
def proxy10(g0):
    v1 = abs(Vp(g0))
    v2 = abs(Vp(g0/2))
    if v2 < 1e-30: return np.nan
    return v1/(g0*v2)

all_proxies = [
    ("f*(g0-1)^2/g0^4", proxy1),
    ("a_eff*|V'|*g0", proxy2),
    ("f*|V'|/g0^3", proxy3),
    ("(f-1)*(g0-1)", proxy4),
    ("a*ln(g0)*(g0-1)", proxy5),
    ("f/(g0*ln(g0))", proxy6),
    ("rho*(1-1/g0)", proxy7),
    ("rho*(g0-1)^2", proxy8),
    ("d(ln f)/d(ln g)", proxy9),
    ("V'(g0)/(g0*V'(g0/2))", proxy10),
]

for name, func in all_proxies:
    try:
        res = minimize_scalar(lambda g: -func(g), bounds=(1.5, 10.0), method='bounded')
        g0s = res.x
        mark = " *** CLOSE ***" if abs(g0s-4.0)<0.15 else ""
        print(f"  {name:>25s}: max at g0 = {g0s:.4f}  (dev: {g0s-4:.3f}){mark}")
    except:
        print(f"  {name:>25s}: ERROR")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  PART 5: THE V'(g0)/V'(g0/alpha) RATIO")
print("="*72)
# ===================================================================

print("""
  From p143: V'(4)/V'(2) = 12 = alpha^2 * d = eta_K^(0).
  This is a REMARKABLE coincidence. Let's understand it.

  V'(g) = g^2*(1-g)
  V'(4) = 16*(-3) = -48
  V'(2) = 4*(-1) = -4
  Ratio: -48/-4 = 12 = alpha^2 * d ✓

  V'(2*alpha) / V'(alpha) = (2*alpha)^2*(1-2*alpha) / (alpha^2*(1-alpha))
  For alpha=2: = 16*(-3)/(4*(-1)) = 48/4 = 12 = 4*3 = alpha^2*d ✓

  General: V'(2a)/V'(a) = 4*(1-2a)/(1-a)
  For a=2: 4*(-3)/(-1) = 12 = a^2*d when a=2, d=3.

  Is 4*(1-2a)/(1-a) = a^2*d a coincidence or structural?
  4*(2a-1)/(a-1) = a^2*d
  For a=2, d=3: 4*3/1 = 12 = 4*3 ✓✓✓

  This is: 4*(2a-1)/(a-1) = a^2*d
  With a=2: 4*3 = 4*3 = 12. This is ALWAYS true for a=2, ANY d!
  Because: 4*(2*2-1)/(2-1) = 4*3/1 = 12 = 4*3 = 2^2*3 = a^2*3.

  Wait: 4*(2a-1)/(a-1) at a=2 gives 12, and a^2*d = 4d.
  So 12 = 4d => d=3. This FIXES d=3!

  Or equivalently: given a=2 (from K=phi^4), the ratio V'(2a)/V'(a)=12
  equals a^2*d only when d=3.

  This is interesting but circular (we input d=3).
""")

# Verify
a = ALPHA
ratio = 4*(2*a-1)/(a-1)
print(f"  V'(2a)/V'(a) = 4*(2a-1)/(a-1) = {ratio:.1f}")
print(f"  a^2*d = {a**2*D:.1f}")
print(f"  Match: {'YES' if abs(ratio-a**2*D)<0.01 else 'NO'}")
print(f"  This relation requires d = 4*(2a-1)/((a-1)*a^2) = {4*(2*a-1)/((a-1)*a**2):.1f}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  PART 6: INFORMATION-THEORETIC CONDITION")
print("="*72)
# ===================================================================

print("""
  Define the KINETIC INFORMATION at field value g:
    I_K(g) = ln K(g) / ln g = 2*alpha  (for K = g^(2*alpha))

  This is constant (= 2*alpha = 4) for the bare theory.

  With running alpha_eff:
    I_K^eff(g) = ln K_eff(g) / ln g
  where K_eff(g) = g^(2*alpha_eff(g)).

  I_K^eff(g) = 2*alpha_eff(g) = 2*alpha/(1+eta*(g-1)^2)

  The K-self-duality condition is:
    I_K^bare(g0) = g0  <=> 2*alpha = g0

  With running:
    I_K^eff(g0) = g0 has solution g0 = 2*alpha/(1+eta*(g0-1)^2) = g0
  which gives 2*alpha = g0*(1+eta*(g0-1)^2) -- a cubic in g0.

  But the K-self-duality uses the BARE coupling, not running.
  WHY? Because the soliton amplitude g0 is set by the UV structure.

  PHYSICAL ARGUMENT:
  The soliton core field excursion (g0-1) probes the UV regime
  of the kinetic function. The running alpha_eff is an IR correction.
  The MAXIMAL soliton is limited by the UV kinetic structure,
  not the IR running.

  The UV kinetic structure says: K(g) = g^(2*alpha_UV).
  The maximal soliton consistent with this structure is at g0 = 2*alpha_UV,
  where the kinetic exponent matches the field value.

  For g0 > 2*alpha: the kinetic function K = g^4 grows SLOWER than g^g0,
  meaning the kinetic term cannot keep up with the field excursion.
  The soliton becomes "kinetically under-supported".

  For g0 < 2*alpha: the kinetic function grows FASTER than g^g0,
  meaning the kinetic term over-supports the field excursion.

  The BALANCE is at g0 = 2*alpha.

  THIS IS THE PHYSICAL CONTENT OF K-SELF-DUALITY.
""")

# Verify: K(g0) vs g0^g0
print(f"  {'g0':>6s} {'K(g0)=g0^4':>14s} {'g0^g0':>14s} {'K/g0^g0':>10s} {'support':>12s}")
print(f"  {'-'*60}")
for g0 in np.arange(2.0, 7.1, 0.5):
    K = g0**(2*ALPHA)
    gg = g0**g0
    ratio = K/gg
    support = "over" if ratio > 1 else "under" if ratio < 1 else "balanced"
    mark = " <-- tau" if abs(g0-4.0)<0.1 else ""
    print(f"  {g0:6.1f} {K:14.1f} {gg:14.1f} {ratio:10.6f} {support:>12s}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  PART 7: COMPARISON OF g0_tau PREDICTIONS")
print("="*72)
# ===================================================================

def get_A(r, g):
    m = (r>=120)&(r<=260)
    rf, tl = r[m], (g[m]-1)*r[m]
    if len(rf)<10: return np.nan
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    c,_,_,_ = np.linalg.lstsq(M, tl, rcond=None)
    return np.sqrt(c[0]**2+c[1]**2)

G0_E = 0.90548144
G0_MU = PHI * G0_E

r_e, g_e, _ = solve_sol(G0_E)
A_e = get_A(r_e, g_e)

print(f"\n  Predictions for various g0_tau candidates:")
print(f"  {'g0_tau':>8s} {'source':>30s} {'m_tau(MeV)':>12s} {'dev(%)':>8s}")
print(f"  {'-'*65}")

candidates = [
    (2*ALPHA, "K-self-duality: 2*alpha"),
    (g0_FK_max, f"F_K maximum: {g0_FK_max:.4f}"),
    (3.989, "g0_e+g0_mu+phi = 3.989"),
    (np.e * G0_MU, f"e*g0_mu = {np.e*G0_MU:.4f}"),
    (np.pi*np.sqrt(2)*G0_E, f"pi*sqrt(2)*g0_e = {np.pi*np.sqrt(2)*G0_E:.4f}"),
]

for g0_tau, source in candidates:
    r_t, g_t, _ = solve_sol(g0_tau)
    if r_t is None:
        print(f"  {g0_tau:8.4f} {source:>30s} {'FAIL':>12s}")
        continue
    A_t = get_A(r_t, g_t)
    r31 = (A_t/A_e)**4
    m_tau = 0.51099895 * r31
    dev = abs(m_tau-1776.86)/1776.86*100
    print(f"  {g0_tau:8.4f} {source:>30s} {m_tau:12.2f} {dev:8.4f}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  PART 8: SELF-CONSISTENT EXTREMUM — COMBINED F_K + KOIDE")
print("="*72)
# ===================================================================

print("""
  Question: Is g0=4 the point where F_K = const * m_tau?
  i.e., is g0=4 special when we combine the MASS PREDICTION
  with the KINETIC EFFICIENCY?

  Define: Lambda(g0) = F_K(g0) * r31(g0) / r31(Koide)
  This is the "efficiency-weighted mass". Does it have an extremum at g0=4?
""")

print(f"  {'g0':>6s} {'F_K':>10s} {'r31':>10s} {'m_tau':>10s} {'Lambda':>12s}")
print(f"  {'-'*55}")

lambdas = []
for g0 in np.arange(3.0, 5.51, 0.1):
    FK = compute_FK(g0)
    r_t, g_t, _ = solve_sol(g0)
    if r_t is None: continue
    A_t = get_A(r_t, g_t)
    r31 = (A_t/A_e)**4
    m_tau = 0.51099895 * r31
    lam = FK * r31 / 3477.48
    lambdas.append((g0, lam))
    mark = " <-- tau" if abs(g0-4.0)<0.05 else ""
    print(f"  {g0:6.2f} {FK:10.6f} {r31:10.2f} {m_tau:10.2f} {lam:12.6f}{mark}")
sys.stdout.flush()

# ===================================================================
print(f"\n{'='*72}")
print("  SYNTHESIS: HONEST ASSESSMENT")
print("="*72)
# ===================================================================

print(f"""
  CORRECTION: F_K maximum is at g0 = {g0_FK_max:.4f}, NOT 4.0.
  The deviation is {abs(g0_FK_max-4.0):.3f} ({abs(g0_FK_max-4.0)/4*100:.1f}%).

  NO analytical proxy tested peaks at EXACTLY g0 = 4.0.

  The K-self-duality remains the SOLE exact algebraic condition:
    K(g0) = g0^g0  <=>  2*alpha = g0  <=>  g0_tau = 4.

  Physical interpretation (REFINED):
  - K(g0) = g0^(2*alpha) is the "kinetic weight" at field value g0.
  - g0^g0 is the field's "self-exponentiation" (maximal self-reference).
  - K = g0^g0 means: the kinetic structure matches the field's self-reference.
  - For g0 > 2*alpha: K(g0) < g0^g0 — kinetically under-supported.
  - For g0 < 2*alpha: K(g0) > g0^g0 — kinetically over-supported.
  - The BALANCE point is g0 = 2*alpha = 4.

  STATUS ASSESSMENT:
  ┌─────────────────────────────────────────────────────────────┐
  │  K-self-duality: g0 = 2*alpha = 4                   [HR]   │
  │  Algebraically exact, structurally motivated                │
  │  Physical: kinetic balance at soliton core                  │
  │  Unique for alpha = 2 (2*alpha = alpha^2)                   │
  │                                                             │
  │  F_K extremum at g0 = {g0_FK_max:.3f}                  [NUM]  │
  │  CLOSE but not EXACT support (3.9% off)                     │
  │  NOT sufficient for [AN] upgrade                            │
  │                                                             │
  │  V'(4)/V'(2) = 12 = alpha^2*d               [AN coincidence]│
  │  Connects g0=4 to eta_K leading order                       │
  │                                                             │
  │  CONCLUSION: g0_tau = 4 remains [HR].                       │
  │  The K-self-duality is the strongest argument,              │
  │  but it is not derivable from a variational principle.      │
  │  No tested proxy gives g0=4 EXACTLY.                        │
  │  This is an HONEST assessment.                              │
  └─────────────────────────────────────────────────────────────┘

  WHAT WOULD BE NEEDED FOR [AN]:
  1. A variational principle S[g0] with d(S)/d(g0)=0 at g0=4, OR
  2. An RG fixed-point equation with g0=4 as its only solution, OR
  3. A topological argument (quantized charge, winding number), OR
  4. Derivation of K-self-duality from the soliton EOM.

  None of these have been found. The problem remains OPEN.
""")

# Final: correct F_K claim for documents
print(f"  DOCUMENT CORRECTION NEEDED:")
print(f"    Old: 'F_K sign change at exactly g0=4.0'")
print(f"    New: 'F_K maximum at g0={g0_FK_max:.3f} (3.9% below 4.0)'")
print(f"    The F_K extremum provides APPROXIMATE, not exact, support.")

print("\nDONE")
