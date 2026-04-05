#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p134d_lean.py - Lean eta_K derivation: precise value + algebraic candidates
============================================================================
Outputs to vault file. Optimized: fewer tau scan points, focused bisection.
"""
import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + np.sqrt(5)) / 2
R21_PDG = 206.768
R31_PDG = 3477.48

# Output file in vault
OUT = os.path.join(os.path.dirname(__file__), "p134d_output.txt")
lines = []
def log(s=""):
    print(s)
    lines.append(s)

def solve_soliton(g0, eta_K=0.0, rm=300):
    def fk(g):
        a = 2.0 / (1 + eta_K * (g - 1)**2)
        return 1 + 2*a*np.log(g) if g > 0 else -1e30
    def Vp(g):
        return g**2*(1-g)
    fg0 = fk(g0)
    if abs(fg0) < 1e-15: return None, None
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
    return r, s.sol(r)[0]

def extract_tail(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan, np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    B, C = coeff
    return np.sqrt(B**2 + C**2), np.arctan2(C, B)

def find_g0e(eta_K):
    def obj(g0_e):
        g0_mu = PHI * g0_e
        r_e, g_e = solve_soliton(g0_e, eta_K)
        if r_e is None or r_e[-1] < 250: return 1e6
        A_e, _ = extract_tail(r_e, g_e)
        r_mu, g_mu = solve_soliton(g0_mu, eta_K)
        if r_mu is None or r_mu[-1] < 250: return 1e6
        A_mu, _ = extract_tail(r_mu, g_mu)
        if np.isnan(A_e) or np.isnan(A_mu) or A_e < 1e-15: return 1e6
        return (A_mu/A_e)**4 - R21_PDG
    for g_lo in np.arange(0.88, 0.93, 0.005):
        try:
            v_lo, v_hi = obj(g_lo), obj(g_lo+0.005)
            if v_lo * v_hi < 0:
                return brentq(obj, g_lo, g_lo+0.005, xtol=1e-9)
        except:
            pass
    return None

def get_max_r31(g0_e, eta_K):
    r_e, g_e = solve_soliton(g0_e, eta_K)
    if r_e is None: return np.nan
    A_e, _ = extract_tail(r_e, g_e)
    if np.isnan(A_e) or A_e < 1e-15: return np.nan
    best = 0
    # Focused scan around expected tau region
    for g0_t in np.linspace(3.0, 5.0, 25):
        r_t, g_t = solve_soliton(g0_t, eta_K)
        if r_t is None or r_t[-1] < 250: continue
        A_t, _ = extract_tail(r_t, g_t)
        if np.isnan(A_t): continue
        r31 = (A_t/A_e)**4
        if r31 > best: best = r31
    return best

# =====================================================================
log("="*70)
log("  PART 1: VERIFY p131 result (eta_K=12.067)")
log("="*70)

g0_e = find_g0e(12.067)
log(f"  g0_e = {g0_e}")
if g0_e:
    r31 = get_max_r31(g0_e, 12.067)
    log(f"  r_31 = {r31:.1f} (target: {R31_PDG:.1f})")

# =====================================================================
log("\n" + "="*70)
log("  PART 2: FINE BISECTION (eta in [12.0, 12.2])")
log("="*70)

# Quick check endpoints
r31_12 = get_max_r31(find_g0e(12.0), 12.0) if find_g0e(12.0) else np.nan
r31_122 = get_max_r31(find_g0e(12.2), 12.2) if find_g0e(12.2) else np.nan
log(f"  eta=12.0: r_31={r31_12:.1f}")
log(f"  eta=12.2: r_31={r31_122:.1f}")

# Determine bracket
if not np.isnan(r31_12) and not np.isnan(r31_122):
    if r31_12 < R31_PDG and r31_122 > R31_PDG:
        eta_lo, eta_hi = 12.0, 12.2
    elif r31_12 > R31_PDG and r31_122 < R31_PDG:
        eta_lo, eta_hi = 12.2, 12.0  # swap
    else:
        # Widen search
        log("  No bracket in [12.0, 12.2], widening...")
        eta_lo, eta_hi = 11.0, 14.0

    log(f"\n  Bisecting in [{eta_lo}, {eta_hi}]:")
    for i in range(25):
        eta_mid = (eta_lo + eta_hi) / 2.0
        g0 = find_g0e(eta_mid)
        if g0 is None:
            eta_lo = eta_mid
            continue
        r31 = get_max_r31(g0, eta_mid)
        if np.isnan(r31):
            eta_lo = eta_mid
            continue
        if r31 > R31_PDG:
            eta_hi = eta_mid
        else:
            eta_lo = eta_mid
        if i % 3 == 0:
            log(f"    iter {i:2d}: eta={eta_mid:.8f}, r_31={r31:.1f}, g0_e={g0:.7f}")
        if abs(eta_hi - eta_lo) < 1e-5:
            break

    eta_opt = (eta_lo + eta_hi) / 2.0
    log(f"\n  >>> RESULT: eta_K = {eta_opt:.8f}")
    g0_final = find_g0e(eta_opt)
    r31_final = get_max_r31(g0_final, eta_opt) if g0_final else np.nan
    log(f"  g0_e = {g0_final}")
    log(f"  r_31 = {r31_final:.1f}")
else:
    eta_opt = 12.067
    log("  Could not bracket, using p131 value: 12.067")

# =====================================================================
log("\n" + "="*70)
log("  PART 3: ALGEBRAIC CANDIDATES")
log("="*70)

candidates = {
    "12":                    12.0,
    "4*pi - 1/2":            4*np.pi - 0.5,
    "12 + 1/15":             12.0 + 1.0/15.0,
    "181/15":                181.0/15.0,
    "12 + pi^2/150":         12.0 + np.pi**2/150.0,
    "2*(pi/2)^4":            2*(np.pi/2)**4,
    "pi^4/8":                np.pi**4/8,
    "12 + 1/phi^4":          12.0 + 1.0/PHI**4,
    "3*(4+1/pi^2)":          3*(4+1/np.pi**2),
    "4*(3+1/(6*pi))":        4*(3+1/(6*np.pi)),
    "12+1/(3*pi)":           12.0+1/(3*np.pi),
    "12+2/(3*pi^2)":         12.0+2/(3*np.pi**2),
    "12+ln(phi)/pi":         12.0+np.log(PHI)/np.pi,
    "48/(pi+1)":             48.0/(np.pi+1),
    "145/12":                145.0/12.0,
    "3*phi^3":               3*PHI**3,
    "36/e":                  36.0/np.e,
    "12+phi-1":              12.0+PHI-1,
    "12+1/e^2":              12.0+1/np.e**2,
    "4*pi/atan(2)":          4*np.pi/np.arctan(2),
    "4*(pi-1/sqrt(pi))":     4*(np.pi-1/np.sqrt(np.pi)),
}

log(f"\n  Reference: eta_K = {eta_opt:.8f}")
log(f"\n  {'Expression':25s} {'Value':>12s} {'Dev %':>10s}")
log(f"  {'-'*50}")
results = [(abs(v-eta_opt)/eta_opt*100, k, v) for k,v in candidates.items()]
results.sort()
for dev, name, val in results[:15]:
    marker = " ***" if dev < 0.005 else " <--" if dev < 0.05 else " ~" if dev < 0.5 else ""
    log(f"  {name:25s} {val:12.8f} {dev:10.4f}%{marker}")

# =====================================================================
log("\n" + "="*70)
log("  PART 4: PHYSICAL INTERPRETATION")
log("="*70)
log(f"""
  eta_K = {eta_opt:.4f} parametrizes field-space screening:
    alpha_eff(g) = 2 / (1 + {eta_opt:.2f} * (g-1)^2)

  At vacuum g=1: alpha_eff = 2.0 (full coupling)
  At e-core g=0.905: alpha_eff = {2.0/(1+eta_opt*0.905**2 - 2*eta_opt*0.905 + eta_opt):.3f}
  At tau g=4.0: alpha_eff = {2.0/(1+eta_opt*9):.4f}

  Field-space Debye radius: 1/sqrt(eta_K) = {1/np.sqrt(eta_opt):.4f}

  Leading order: eta_K = alpha_UV^2 * d = 4*3 = 12
  Subleading: delta = {eta_opt - 12:.4f} ({(eta_opt-12)/eta_opt*100:.2f}%)

  STATUS: eta_K is a numerical constant of TGP.
  Leading analytical form: alpha^2*d = 12 captures 99.4%.
""")

log("="*70)
log("  DONE")
log("="*70)

# Save to vault
with open(OUT, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
log(f"\n  Output saved to: {OUT}")
