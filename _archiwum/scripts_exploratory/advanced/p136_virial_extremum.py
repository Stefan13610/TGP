#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p136_virial_extremum.py - Virial ratio extremum as tau selector
================================================================
From p135: E_kin/E_pot has an extremum near g0~3.9-4.0 (tau region!)
This could be the missing constraint on g0_tau.

Also test: g0_tau = e * g0_mu (algebraic relation, 0.43% match)
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import minimize_scalar, brentq

PHI = (1 + np.sqrt(5)) / 2
ETA_K = 12.067
G0_E = 0.905481
G0_MU = PHI * G0_E
R31_PDG = 3477.48

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

def extract_tail(r, g):
    m = (r >= 120) & (r <= 260)
    rf, tl = r[m], (g[m] - 1) * r[m]
    if len(rf) < 10: return np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    coeff, _, _, _ = np.linalg.lstsq(A, tl, rcond=None)
    return np.sqrt(coeff[0]**2 + coeff[1]**2)

def V_func(g): return g**3/3.0 - g**4/4.0

def fk_func(g, eta_K=ETA_K):
    a = 2.0 / (1 + eta_K * (g - 1)**2)
    return 1 + 2*a*np.log(g) if g > 0 else -1e30

r_e, g_e, gp_e = solve_soliton(G0_E)
A_e = extract_tail(r_e, g_e)

def soliton_energies(g0):
    r, g, gp = solve_soliton(g0)
    if r is None: return np.nan, np.nan
    fk_arr = np.array([fk_func(gi) for gi in g])
    T = 0.5 * fk_arr * gp**2
    V = np.array([V_func(gi) - V_func(1.0) for gi in g])
    E_kin = 4*np.pi * trapezoid(T * r**2, r)
    E_pot = 4*np.pi * trapezoid(V * r**2, r)
    return E_kin, E_pot

def virial_ratio(g0):
    """E_kin / |E_pot| — look for extremum"""
    ek, ep = soliton_energies(g0)
    if np.isnan(ek) or abs(ep) < 1e-10: return np.nan
    return ek / abs(ep)

# =====================================================================
print("="*70)
print("  PART 1: FINE SCAN OF VIRIAL RATIO E_kin/|E_pot|")
print("="*70)

print(f"\n  {'g0':>6s} {'E_kin/|E_pot|':>14s} {'dR/dg0':>10s}")
print(f"  {'-'*35}")

g0_vals = np.linspace(2.0, 7.0, 60)
vr_vals = []
for g0 in g0_vals:
    vr = virial_ratio(g0)
    vr_vals.append(vr)
vr_vals = np.array(vr_vals)

# Compute numerical derivative
for i, g0 in enumerate(g0_vals):
    if np.isnan(vr_vals[i]): continue
    if i > 0 and i < len(g0_vals)-1 and not np.isnan(vr_vals[i-1]) and not np.isnan(vr_vals[i+1]):
        dvr = (vr_vals[i+1] - vr_vals[i-1]) / (g0_vals[i+1] - g0_vals[i-1])
    else:
        dvr = np.nan
    mark = ""
    if abs(g0 - 4.0) < 0.1: mark = " <-- tau region"
    if not np.isnan(dvr) and abs(dvr) < 0.0005: mark += " EXTREMUM?"
    # Only print every 3rd for readability
    if i % 3 == 0 or (not np.isnan(dvr) and abs(dvr) < 0.001):
        print(f"  {g0:6.3f} {vr_vals[i]:14.8f} {dvr if not np.isnan(dvr) else 0:10.6f}{mark}")

# Find exact extremum
print(f"\n  Finding exact extremum of E_kin/|E_pot|...")

# Use scipy minimize_scalar (we negate to find maximum)
def neg_virial(g0):
    vr = virial_ratio(g0)
    return -vr if not np.isnan(vr) else 0.0

result = minimize_scalar(neg_virial, bounds=(3.0, 6.0), method='bounded',
                         options={'xatol': 1e-6})
g0_extremum = result.x
vr_extremum = virial_ratio(g0_extremum)

print(f"  Virial ratio MAXIMUM at g0 = {g0_extremum:.6f}")
print(f"  E_kin/|E_pot| = {vr_extremum:.8f}")

# What r_31 does this give?
r_t, g_t, _ = solve_soliton(g0_extremum)
A_t = extract_tail(r_t, g_t)
r31_virial = (A_t / A_e)**4
m_tau_virial = 0.51099895 * r31_virial

print(f"  A_tau = {A_t:.8f}")
print(f"  r_31 = {r31_virial:.1f} (Koide: {R31_PDG:.1f})")
print(f"  m_tau = {m_tau_virial:.1f} MeV (obs: 1776.86)")
print(f"  Deviation: {abs(m_tau_virial - 1776.86)/1776.86*100:.2f}%")

# Ratio to g0_e and g0_mu
print(f"\n  g0(virial) / g0_e  = {g0_extremum/G0_E:.6f}")
print(f"  g0(virial) / g0_mu = {g0_extremum/G0_MU:.6f}")
print(f"  Compare: e = {np.e:.6f}")
print(f"  Compare: phi^2 = {PHI**2:.6f}")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART 2: OTHER VIRIAL-TYPE FUNCTIONALS")
print(f"{'='*70}")

# Try other functionals that might have extrema at tau:
# F1: E_kin / E_tot
# F2: E_kin * g0^2 / E_pot
# F3: (E_kin - E_pot) / (E_kin + E_pot)  [virial parameter]
# F4: d(ln E_tot)/d(ln g0)  [scaling exponent]

def functional_scan(g0):
    ek, ep = soliton_energies(g0)
    if np.isnan(ek) or np.isnan(ep): return {}
    et = ek + ep
    return {
        'E_kin': ek,
        'E_pot': ep,
        'E_tot': et,
        'F1': ek / et if abs(et) > 1e-10 else np.nan,
        'F2': ek * g0**2 / abs(ep) if abs(ep) > 1e-10 else np.nan,
        'F3': (ek - abs(ep)) / (ek + abs(ep)) if (ek + abs(ep)) > 0 else np.nan,
    }

print(f"\n  {'g0':>6s} {'F1:Ek/Et':>12s} {'F2:Ek*g^2/Ep':>14s} {'F3:virial_p':>12s}")
print(f"  {'-'*50}")

f1_vals, f2_vals, f3_vals = [], [], []
g0_fine = np.linspace(2.0, 7.0, 50)
for g0 in g0_fine:
    fs = functional_scan(g0)
    f1_vals.append(fs.get('F1', np.nan))
    f2_vals.append(fs.get('F2', np.nan))
    f3_vals.append(fs.get('F3', np.nan))

f1_vals = np.array(f1_vals)
f2_vals = np.array(f2_vals)
f3_vals = np.array(f3_vals)

# Find extrema of each
for name, vals in [('F1: E_kin/E_tot', f1_vals), ('F2: E_kin*g^2/|E_pot|', f2_vals), ('F3: virial param', f3_vals)]:
    if np.all(np.isnan(vals)): continue
    valid = ~np.isnan(vals)
    if not np.any(valid): continue

    # Check for sign changes in derivative
    dv = np.diff(vals[valid])
    sign_changes = np.where(np.diff(np.sign(dv)))[0]
    g_valid = g0_fine[valid]

    if len(sign_changes) > 0:
        for sc in sign_changes:
            g_ext = g_valid[sc+1]
            v_ext = vals[valid][sc+1]
            r_t2, g_t2, _ = solve_soliton(g_ext)
            if r_t2 is not None:
                A_ext = extract_tail(r_t2, g_t2)
                r31_ext = (A_ext/A_e)**4 if not np.isnan(A_ext) else np.nan
                print(f"  {name}: extremum at g0={g_ext:.3f}, val={v_ext:.6f}, r_31={r31_ext:.0f}")
    else:
        print(f"  {name}: monotonic in [2, 7] (no extremum)")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART 3: ALGEBRAIC TEST: g0_tau = e * g0_mu")
print(f"{'='*70}")

g0_tau_e = np.e * G0_MU
print(f"  g0_tau = e * g0_mu = {np.e:.6f} * {G0_MU:.6f} = {g0_tau_e:.6f}")
r_t3, g_t3, _ = solve_soliton(g0_tau_e)
A_t3 = extract_tail(r_t3, g_t3)
r31_e = (A_t3 / A_e)**4
m_tau_e = 0.51099895 * r31_e
print(f"  r_31 = {r31_e:.1f}")
print(f"  m_tau = {m_tau_e:.1f} MeV (obs: 1776.86, dev: {abs(m_tau_e-1776.86)/1776.86*100:.2f}%)")

# Also test g0_tau = phi^2 * g0_mu
g0_tau_phi2 = PHI**2 * G0_MU
print(f"\n  g0_tau = phi^2 * g0_mu = {PHI**2:.6f} * {G0_MU:.6f} = {g0_tau_phi2:.6f}")
r_t4, g_t4, _ = solve_soliton(g0_tau_phi2)
A_t4 = extract_tail(r_t4, g_t4)
r31_phi2 = (A_t4 / A_e)**4
m_tau_phi2 = 0.51099895 * r31_phi2
print(f"  r_31 = {r31_phi2:.1f}")
print(f"  m_tau = {m_tau_phi2:.1f} MeV (obs: 1776.86, dev: {abs(m_tau_phi2-1776.86)/1776.86*100:.2f}%)")

# g0_tau = pi*sqrt(2) * g0_e
g0_tau_pi = np.pi * np.sqrt(2) * G0_E
print(f"\n  g0_tau = pi*sqrt(2) * g0_e = {np.pi*np.sqrt(2):.6f} * {G0_E:.6f} = {g0_tau_pi:.6f}")
r_t5, g_t5, _ = solve_soliton(g0_tau_pi)
A_t5 = extract_tail(r_t5, g_t5)
r31_pi = (A_t5 / A_e)**4
m_tau_pi = 0.51099895 * r31_pi
print(f"  r_31 = {r31_pi:.1f}")
print(f"  m_tau = {m_tau_pi:.1f} MeV (obs: 1776.86, dev: {abs(m_tau_pi-1776.86)/1776.86*100:.2f}%)")

# =====================================================================
print(f"\n{'='*70}")
print(f"  PART 4: ENERGY PER UNIT MASS")
print(f"{'='*70}")

# If mass ~ A^4, then energy/mass = E_tot / A^4
# The most efficient soliton (minimum E/m) would be preferred
# (analogous to BPS bound in SUSY theories)

print(f"\n  E_tot / A^4 (energy per unit mass) — minimum = most efficient")
print(f"\n  {'g0':>6s} {'E_tot':>10s} {'A^4':>10s} {'E/m':>10s}")
print(f"  {'-'*40}")

em_vals = []
for g0 in np.linspace(2.0, 8.0, 50):
    r, g, gp = solve_soliton(g0)
    if r is None: continue
    A = extract_tail(r, g)
    fk_arr = np.array([fk_func(gi) for gi in g])
    T = 0.5 * fk_arr * gp**2
    V = np.array([V_func(gi) - V_func(1.0) for gi in g])
    et = 4*np.pi * trapezoid((T + V) * r**2, r)
    if np.isnan(A) or A < 1e-15: continue
    a4 = A**4
    em = et / a4
    em_vals.append((g0, et, a4, em))

# Find minimum
if em_vals:
    min_em = min(em_vals, key=lambda x: x[3])
    for g0, et, a4, em in em_vals:
        mark = ""
        if abs(g0 - min_em[0]) < 0.15: mark = " <-- MIN"
        if abs(g0 - 4.0) < 0.15: mark += " tau"
        if len(mark) > 0 or em_vals.index((g0, et, a4, em)) % 4 == 0:
            print(f"  {g0:6.3f} {et:10.2f} {a4:10.4f} {em:10.4f}{mark}")

    print(f"\n  MINIMUM E/m at g0 = {min_em[0]:.3f}")
    print(f"  E/m = {min_em[3]:.4f}")
    r_min, g_min, _ = solve_soliton(min_em[0])
    A_min = extract_tail(r_min, g_min)
    r31_min = (A_min / A_e)**4
    print(f"  r_31 = {r31_min:.1f} (Koide: {R31_PDG})")
    print(f"  m_tau = {0.51099895*r31_min:.1f} MeV (obs: 1776.86)")

# =====================================================================
print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"""
  CANDIDATE MECHANISMS FOR FIXING g0_tau:

  1. VIRIAL RATIO EXTREMUM: E_kin/|E_pot| max at g0 = {g0_extremum:.3f}
     r_31 = {r31_virial:.0f}, m_tau = {m_tau_virial:.0f} MeV

  2. ALGEBRAIC: g0_tau = e * g0_mu = {g0_tau_e:.4f}
     r_31 = {r31_e:.0f}, m_tau = {m_tau_e:.0f} MeV (dev {abs(m_tau_e-1776.86)/1776.86*100:.1f}%)

  3. ALGEBRAIC: g0_tau = phi^2 * g0_mu = {g0_tau_phi2:.4f}
     r_31 = {r31_phi2:.0f}, m_tau = {m_tau_phi2:.0f} MeV (dev {abs(m_tau_phi2-1776.86)/1776.86*100:.1f}%)

  4. ALGEBRAIC: g0_tau = pi*sqrt(2) * g0_e = {g0_tau_pi:.4f}
     r_31 = {r31_pi:.0f}, m_tau = {m_tau_pi:.0f} MeV (dev {abs(m_tau_pi-1776.86)/1776.86*100:.1f}%)

  Exact: g0_tau = 3.99987 for r_31 = Koide.
""")
print("DONE")
