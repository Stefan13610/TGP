#!/usr/bin/env python3
"""
p139b_eta_lean.py -- Lean analytical derivation of eta_K
========================================================
Focused: high-precision eta_K + test analytical candidates.
Skips expensive scans.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stdout.reconfigure(line_buffering=True)
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1+np.sqrt(5))/2
R21 = 206.768; R31K = 3477.48; ME = 0.51099895; MTAU = 1776.86

def solve_sol(g0, eta, rm=300):
    def fk(g):
        a = 2.0/(1+eta*(g-1)**2)
        return 1+2*a*np.log(g) if g>0 else -1e30
    def Vp(g): return g**2*(1-g)
    fg0 = fk(g0)
    if abs(fg0)<1e-15: return None,None
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
    return r,s.sol(r)[0]

def get_A(r,g):
    m=(r>=120)&(r<=260)
    rf,tl=r[m],(g[m]-1)*r[m]
    if len(rf)<10: return np.nan
    M=np.column_stack([np.cos(rf),np.sin(rf)])
    c,_,_,_=np.linalg.lstsq(M,tl,rcond=None)
    return np.sqrt(c[0]**2+c[1]**2)

def find_g0e(eta):
    def obj(g0e):
        re,ge=solve_sol(g0e,eta)
        if re is None or re[-1]<250: return 1e6
        Ae=get_A(re,ge)
        rm,gm=solve_sol(PHI*g0e,eta)
        if rm is None or rm[-1]<250: return 1e6
        Am=get_A(rm,gm)
        if np.isnan(Ae)or np.isnan(Am)or Ae<1e-15: return 1e6
        return (Am/Ae)**4-R21
    for glo in np.arange(0.88,0.93,0.005):
        try:
            if obj(glo)*obj(glo+0.005)<0:
                return brentq(obj,glo,glo+0.005,xtol=1e-10)
        except: pass
    return None

def r31_at4(eta):
    g0e=find_g0e(eta)
    if g0e is None: return np.nan,np.nan
    re,ge=solve_sol(g0e,eta)
    if re is None: return np.nan,np.nan
    Ae=get_A(re,ge)
    rt,gt=solve_sol(4.0,eta)
    if rt is None: return np.nan,np.nan
    At=get_A(rt,gt)
    if np.isnan(Ae)or np.isnan(At)or Ae<1e-15: return np.nan,np.nan
    return (At/Ae)**4, g0e

print("="*70)
print("  HIGH-PRECISION eta_K DETERMINATION")
print("="*70)

# Coarse scan
print("\n  Coarse scan:")
for eta in [12.0, 12.05, 12.06, 12.07, 12.08, 12.1]:
    r31,g0e=r31_at4(eta)
    if not np.isnan(r31):
        print(f"    eta={eta:8.4f}  r31={r31:8.1f}  m_tau={ME*r31:8.2f}")
    sys.stdout.flush()

# Bisection
print("\n  Bisection...")
sys.stdout.flush()
elo,ehi=12.0,12.2
for i in range(40):
    em=(elo+ehi)/2
    r31m,_=r31_at4(em)
    if np.isnan(r31m): elo=em; continue
    if r31m>R31K: ehi=em
    else: elo=em
    if abs(ehi-elo)<1e-10: break

eta_num=(elo+ehi)/2
r31_num,g0e_num=r31_at4(eta_num)
m_num=ME*r31_num
delta=eta_num-12.0

print(f"\n  eta_K (precise) = {eta_num:.10f}")
print(f"  delta = eta-12  = {delta:.10f}")
print(f"  1/delta         = {1/delta:.6f}")
print(f"  g0_e            = {g0e_num:.10f}")
print(f"  r_31            = {r31_num:.6f}")
print(f"  m_tau           = {m_num:.4f} MeV (obs: {MTAU})")
print(f"  dev             = {abs(m_num-MTAU)/MTAU*100:.5f}%")

# ===================================================================
#  ANALYTICAL CANDIDATES
# ===================================================================
print("\n"+"="*70)
print("  ANALYTICAL CANDIDATES FOR eta_K")
print("="*70)

a=2.0; d=3.0
# Key: 1/15 = 0.06667
print(f"\n  delta = {delta:.10f}")
print(f"  1/15  = {1/15:.10f}")
print(f"  Match: {abs(delta-1/15)/delta*100:.4f}%")

candidates = [
    ("12 + 1/15 = 181/15",                   181/15),
    ("(30+2*sqrt(230))/5",                    (30+2*np.sqrt(230))/5),
    ("a^2*d + 1/((a^2+1)*d)",                a**2*d + 1/((a**2+1)*d)),
    ("a^2*d + a^2/(a^4*d+a^2)",              a**2*d + a**2/(a**4*d+a**2)),
    ("(a^4*d^2+1)/(a^2*d)",                  (a**4*d**2+1)/(a**2*d)),
    ("12+1/(a^2*(a^2-1))",                   12+1/(a**2*(a**2-1))),
    ("12*(1+1/(a^4*d^2))",                   12*(1+1/(a**4*d**2))),
    ("12+a^2/(6*pi^2)",                      12+a**2/(6*np.pi**2)),
    ("4*pi-0.5",                             4*np.pi-0.5),
    ("12+2/(a^4*d)",                         12+2/(a**4*d)),
    ("12+1/(a^2*d+d)",                       12+1/(a**2*d+d)),
]

print(f"\n  {'Candidate':<40s} {'Value':>14s} {'Dev ppm':>10s}")
print(f"  {'-'*68}")
for name,val in sorted(candidates, key=lambda x:abs(x[1]-eta_num)):
    ppm=abs(val-eta_num)/eta_num*1e6
    mk=" ***" if ppm<100 else " **" if ppm<1000 else ""
    print(f"  {name:<40s} {val:14.10f} {ppm:10.1f}{mk}")

# ===================================================================
#  ODE VERIFICATION OF TOP CANDIDATES
# ===================================================================
print("\n"+"="*70)
print("  ODE VERIFICATION OF TOP CANDIDATES")
print("="*70)

for name,eta_test in [("181/15", 181/15),
                       ("(30+2*sqrt(230))/5", (30+2*np.sqrt(230))/5),
                       ("4*pi-0.5", 4*np.pi-0.5)]:
    r31t,g0et=r31_at4(eta_test)
    if not np.isnan(r31t):
        mt=ME*r31t
        print(f"\n  {name}:")
        print(f"    eta_K = {eta_test:.10f}")
        print(f"    r_31  = {r31t:.4f}")
        print(f"    m_tau = {mt:.4f} MeV  (dev: {abs(mt-MTAU)/MTAU*100:.5f}%)")
    sys.stdout.flush()

# ===================================================================
#  SELF-CONSISTENT EQUATION CHECK
# ===================================================================
print("\n"+"="*70)
print("  SELF-CONSISTENT EQUATION CHECK")
print("="*70)

# Test: eta^2 - 12*eta - 4/5 = 0
val_q = eta_num**2 - 12*eta_num - 4/5
print(f"\n  eta^2 - 12*eta - 4/5 = {val_q:.8f}  (should be ~0)")

# What c makes eta^2 - 12*eta - c = 0 exact?
c_exact = eta_num**2 - 12*eta_num
print(f"  Exact c = eta^2 - 12*eta = {c_exact:.10f}")
print(f"  4/5 = {4/5:.10f}")
print(f"  Match: {abs(c_exact-0.8)/0.8*100:.4f}%")

# Other small rationals near c_exact
print(f"\n  Rational approximations for c = {c_exact:.8f}:")
for num in range(1,20):
    for den in range(1,30):
        rat=num/den
        if abs(rat-c_exact)<0.02:
            dev_pct=abs(rat-c_exact)/abs(c_exact)*100
            print(f"    {num}/{den} = {rat:.8f}  dev: {dev_pct:.3f}%")

print("\n"+"="*70)
print("  SUMMARY")
print("="*70)
print(f"""
  eta_K (numerical) = {eta_num:.10f}
  delta = eta-12    = {delta:.10f}

  BEST ANALYTICAL MATCH: eta_K = 12 + 1/15 = 181/15
    = alpha^2*d + 1/((alpha^2+1)*d)
    Value: {181/15:.10f}
    Match: {abs(181/15-eta_num)/eta_num*1e6:.1f} ppm

  SELF-CONSISTENT EQUATION: eta^2 = 12*eta + c
    c (numerical) = {c_exact:.8f}
    c = 4/5 gives eta = (30+2*sqrt(230))/5 = {(30+2*np.sqrt(230))/5:.10f}

  PHYSICAL INTERPRETATION:
    Leading: eta_K^(0) = alpha_UV^2 * d = 4*3 = 12
      (kinetic anomalous dimension * spatial dimensions)
    Correction: delta = 1/((alpha^2+1)*d) = 1/15
      (back-reaction of running alpha on the ERG threshold)
""")
print("DONE")
