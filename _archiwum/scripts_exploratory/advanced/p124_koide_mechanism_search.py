#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p124_koide_mechanism_search.py -- Search for Koide Q=3/2 mechanism in TGP
==========================================================================

The Koide formula: Q = (sum sqrt(m_i))^2 / (sum m_i) = 3/2
is satisfied to <0.02% by charged lepton masses.

In TGP, masses come from soliton tail amplitudes:
  m_n ~ A_tail(g_0^(n))^4
The phi-FP gives g0_mu = phi*g0_e with r_21 = 206.768 (exact PASS).
But r_31 = 3477 requires either:
  (a) a new selection principle for g0_tau, or
  (b) a structural property of the soliton ODE that enforces Q=3/2.

This script:
  1. Finds g0_tau by bisection to match Koide r_31
  2. Tests multiple selection principles for g0_tau
  3. Analyzes phase, symmetry, and geometric properties

Uses exact same ODE solver as phi_fp_bisection.py for consistency.

Author: TGP project, session v42+
Date: 2026-03-31
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.interpolate import interp1d

# === Constants ===
PHI    = (1 + np.sqrt(5)) / 2       # golden ratio = 1.6180339887
ALPHA  = 2
GG     = np.exp(-1 / (2 * ALPHA))   # ghost point g* = 0.77880
R21_PDG = 206.768
R31_PDG = 3477.48

# Koide-predicted r_31 from r_21
_a_K = 1 + np.sqrt(R21_PDG)
_disc_K = 6 * _a_K**2 - 3 - 3 * R21_PDG
R31_KOIDE = (2 * _a_K + np.sqrt(_disc_K))**2  # = 3477.44

# === Soliton ODE solver (identical to phi_fp_bisection.py) ===

def fk(g):
    return 1 + 2 * ALPHA * np.log(g) if g > 0 else -1e30

def Vp(g):
    return g**2 * (1 - g)

def V(g):
    return g**3/3 - g**4/4

def solve_soliton(g0, rm=300, np_=25000):
    """Solve soliton ODE: f(g)*g'' + (2/r)*g' = V'(g)."""
    fg0 = fk(g0)
    if abs(fg0) < 1e-15:
        return None, None
    c2 = Vp(g0) / (3 * fg0)
    rs = 0.01
    def rhs(r, y):
        g, p = y
        if g <= 1e-15:
            return [p, 0]
        fg = fk(g)
        if abs(fg) < 1e-10:
            return [p, 0]
        if r < 1e-10:
            return [p, Vp(g) / fg / 3]
        return [p, (Vp(g) - 2/r * p) / fg]
    def ev(r, y):
        return 100 - abs(y[0])
    ev.terminal = True
    s = solve_ivp(rhs, [rs, rm], [g0 + c2*rs**2, 2*c2*rs],
                  method='RK45', rtol=1e-11, atol=1e-13,
                  max_step=0.05, events=[ev], dense_output=True)
    r = np.linspace(rs, min(s.t[-1], rm), np_)
    return r, s.sol(r)[0]

def tail_amplitude(r, g, r_min=120, r_max=260):
    """Extract A_tail from oscillatory tail."""
    m = (r >= r_min) & (r <= r_max)
    rf = r[m]
    tl = (g[m] - 1) * rf
    if len(rf) < 10:
        return np.nan, np.nan, np.nan
    A = np.column_stack([np.cos(rf), np.sin(rf)])
    B, C = np.linalg.lstsq(A, tl, rcond=None)[0]
    amp = np.sqrt(B**2 + C**2)
    phase = np.arctan2(C, B)
    return amp, phase, np.sqrt(np.mean((tl - A @ [B, C])**2))

def mass_ratio(g0_1, g0_2):
    """Compute mass ratio m_2/m_1 = (A_2/A_1)^4."""
    r1, g1 = solve_soliton(g0_1)
    r2, g2 = solve_soliton(g0_2)
    if r1 is None or r2 is None:
        return np.nan
    if r1[-1] < 250 or r2[-1] < 250:
        return np.nan
    A1, _, _ = tail_amplitude(r1, g1)
    A2, _, _ = tail_amplitude(r2, g2)
    if A1 < 1e-15:
        return np.nan
    return (A2 / A1)**4


# ================================================================
#  STEP 0: Establish baseline (phi-FP)
# ================================================================

def step0_baseline():
    """Reproduce phi-FP result: g0_e, g0_mu, r_21."""
    print("="*60)
    print("  STEP 0: Baseline phi-FP")
    print("="*60)

    g0_e = 0.8992655880
    g0_mu = PHI * g0_e

    r1, g1 = solve_soliton(g0_e)
    r2, g2 = solve_soliton(g0_mu)
    A_e, phase_e, _ = tail_amplitude(r1, g1)
    A_mu, phase_mu, _ = tail_amplitude(r2, g2)
    r21 = (A_mu / A_e)**4

    print(f"  g0_e  = {g0_e:.10f}")
    print(f"  g0_mu = {g0_mu:.10f} (= phi*g0_e)")
    print(f"  A_e   = {A_e:.10f},  phase_e  = {phase_e:.6f} rad")
    print(f"  A_mu  = {A_mu:.10f},  phase_mu = {phase_mu:.6f} rad")
    print(f"  r_21  = {r21:.6f} (PDG: {R21_PDG})")
    print(f"  [PASS]" if abs(r21 - R21_PDG) < 0.01 else f"  [FAIL]")

    return g0_e, g0_mu, A_e, A_mu, phase_e, phase_mu


# ================================================================
#  STEP 1: Find g0_tau by bisection (target: r_31 = Koide)
# ================================================================

def step1_find_g0_tau(g0_e, A_e):
    """Find g0_tau such that r_31 = R31_KOIDE."""
    print("\n" + "="*60)
    print("  STEP 1: Find g0_tau for Koide r_31")
    print("="*60)

    print(f"  Target r_31 (Koide) = {R31_KOIDE:.2f}")
    print(f"  Target r_31 (PDG)   = {R31_PDG}")

    # First: scan r_31(g0_tau) landscape
    g0_scan = np.linspace(0.80, 4.0, 200)
    r31_scan = []
    for g0 in g0_scan:
        r, g = solve_soliton(g0)
        if r is None or r[-1] < 250:
            r31_scan.append(np.nan)
            continue
        A, _, _ = tail_amplitude(r, g)
        if A < 1e-15 or np.isnan(A):
            r31_scan.append(np.nan)
            continue
        r31_scan.append((A / A_e)**4)
    r31_scan = np.array(r31_scan)

    valid = ~np.isnan(r31_scan) & (r31_scan > 0)
    print(f"\n  Scan: {np.sum(valid)}/{len(g0_scan)} valid solutions")
    print(f"  r_31 range: [{np.nanmin(r31_scan[valid]):.1f}, {np.nanmax(r31_scan[valid]):.1f}]")

    # Find where r_31 crosses the Koide target
    # r_31 is monotonically increasing with g0 (away from g0_e)
    crossings = []
    for i in range(len(g0_scan)-1):
        if valid[i] and valid[i+1]:
            if (r31_scan[i] - R31_KOIDE) * (r31_scan[i+1] - R31_KOIDE) < 0:
                crossings.append((g0_scan[i], g0_scan[i+1]))

    if not crossings:
        # Try wider range
        print("  No crossing found in [0.80, 4.0]. Extending to [4.0, 10.0]...")
        g0_ext = np.linspace(4.0, 10.0, 200)
        for g0 in g0_ext:
            r, g = solve_soliton(g0, rm=300)
            if r is None or r[-1] < 250:
                continue
            A, _, _ = tail_amplitude(r, g)
            if A < 1e-15 or np.isnan(A):
                continue
            r31_val = (A / A_e)**4
            g0_scan = np.append(g0_scan, g0)
            r31_scan = np.append(r31_scan, r31_val)
            if len(r31_scan) >= 2 and not np.isnan(r31_scan[-2]):
                if (r31_scan[-2] - R31_KOIDE) * (r31_val - R31_KOIDE) < 0:
                    crossings.append((g0_ext[max(0, len(g0_ext)-2)], g0))

    if crossings:
        lo, hi = crossings[0]
        print(f"  Crossing bracket: [{lo:.4f}, {hi:.4f}]")

        # Bisection
        def r31_func(g0_tau):
            r, g = solve_soliton(g0_tau)
            if r is None or r[-1] < 250:
                return np.nan
            A, _, _ = tail_amplitude(r, g)
            return (A / A_e)**4 - R31_KOIDE

        try:
            g0_tau = brentq(r31_func, lo, hi, xtol=1e-10, maxiter=200)
            r_tau, g_tau = solve_soliton(g0_tau)
            A_tau, phase_tau, resid = tail_amplitude(r_tau, g_tau)
            r31_check = (A_tau / A_e)**4

            print(f"\n  g0_tau (Koide) = {g0_tau:.10f}")
            print(f"  A_tau          = {A_tau:.10f}")
            print(f"  phase_tau      = {phase_tau:.6f} rad")
            print(f"  r_31           = {r31_check:.2f} (target: {R31_KOIDE:.2f})")
            return g0_tau, A_tau, phase_tau
        except Exception as e:
            print(f"  Bisection failed: {e}")
    else:
        print("  No crossing found. r_31 may not reach Koide target.")
        # Report maximum r_31
        valid2 = ~np.isnan(r31_scan)
        if np.any(valid2):
            idx_max = np.nanargmax(r31_scan)
            print(f"  Max r_31 = {r31_scan[idx_max]:.1f} at g0 = {g0_scan[idx_max]:.4f}")

    return None, None, None


# ================================================================
#  STEP 2: Analyze g0_tau properties
# ================================================================

def step2_analyze(g0_e, g0_mu, g0_tau, phase_e, phase_mu, phase_tau):
    """Analyze what makes g0_tau special."""
    print("\n" + "="*60)
    print("  STEP 2: Analyze g0_tau Geometric Properties")
    print("="*60)

    if g0_tau is None:
        print("  [SKIPPED: g0_tau not found]")
        return {}

    results = {}

    # --- 2a: Ratios and algebraic properties ---
    print("\n  --- 2a: Algebraic ratios ---")
    ratio_tau_e = g0_tau / g0_e
    ratio_tau_mu = g0_tau / g0_mu
    ratio_mu_e = g0_mu / g0_e

    print(f"  g0_tau / g0_e  = {ratio_tau_e:.6f}")
    print(f"  g0_tau / g0_mu = {ratio_tau_mu:.6f}")
    print(f"  g0_mu / g0_e   = {ratio_mu_e:.6f} (= phi = {PHI:.6f})")

    # Check common ratios
    specials = {
        'phi^2': PHI**2,
        'phi^3': PHI**3,
        'e': np.e,
        'pi': np.pi,
        'phi*e/2': PHI*np.e/2,
        '2*phi': 2*PHI,
        '3': 3.0,
        'phi+1': PHI + 1,  # = phi^2
        'sqrt(2)*phi': np.sqrt(2)*PHI,
        'phi*(1+1/phi)': PHI*(1+1/PHI),  # = phi^2
    }
    for name, val in specials.items():
        if abs(ratio_tau_e - val) < 0.05 * val:
            print(f"  ** g0_tau/g0_e ~ {name} = {val:.6f} (diff {abs(ratio_tau_e-val)/val*100:.2f}%)")

    # --- 2b: Relation to ghost point g* ---
    print("\n  --- 2b: Relation to ghost point ---")
    delta_e = g0_e - GG
    delta_mu = g0_mu - GG
    delta_tau = g0_tau - GG

    print(f"  g0_e  - g* = {delta_e:.6f}")
    print(f"  g0_mu - g* = {delta_mu:.6f}")
    print(f"  g0_tau- g* = {delta_tau:.6f}")
    print(f"  (g0_mu-g*)/(g0_e-g*)  = {delta_mu/delta_e:.6f}")
    print(f"  (g0_tau-g*)/(g0_e-g*) = {delta_tau/delta_e:.6f}")
    print(f"  (g0_tau-g*)/(g0_mu-g*) = {delta_tau/delta_mu:.6f}")

    # --- 2c: Phase analysis ---
    print("\n  --- 2c: Phase spacing ---")
    dp_emu = phase_mu - phase_e
    dp_mutau = phase_tau - phase_mu
    dp_etau = phase_tau - phase_e

    print(f"  phase_e   = {phase_e:.6f} rad")
    print(f"  phase_mu  = {phase_mu:.6f} rad")
    print(f"  phase_tau = {phase_tau:.6f} rad")
    print(f"  Delta_phase(e->mu)  = {dp_emu:.6f} rad")
    print(f"  Delta_phase(mu->tau) = {dp_mutau:.6f} rad")
    print(f"  Delta_phase(e->tau) = {dp_etau:.6f} rad")
    print(f"  Ratio dp(mu->tau)/dp(e->mu) = {dp_mutau/dp_emu:.6f}" if abs(dp_emu) > 1e-10 else "  dp(e->mu) ~ 0")
    print(f"  dp(e->mu) / (2*pi/3) = {dp_emu/(2*np.pi/3):.6f}")
    print(f"  dp(e->tau) / (4*pi/3) = {dp_etau/(4*np.pi/3):.6f}")

    results['ratio_tau_e'] = ratio_tau_e
    results['delta_phase_emu'] = dp_emu
    results['delta_phase_mutau'] = dp_mutau

    # --- 2d: Soliton energy analysis ---
    print("\n  --- 2d: Soliton energy and action ---")
    for label, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau)]:
        r, g = solve_soliton(g0)
        if r is None:
            continue
        dr = r[1] - r[0]
        gp = np.gradient(g, dr)
        fg = np.array([fk(gi) for gi in g])

        # Energy density: E = (1/2)*f(g)*(g')^2 + V(g)
        E_dens = 0.5 * fg * gp**2 + np.array([V(gi) for gi in g])
        # Total energy (spherical): E_tot = 4*pi*int r^2 * E_dens dr
        E_tot = 4 * np.pi * np.trapz(r**2 * E_dens, r)

        # Action = int (1/2)*f(g)*(g')^2 * 4*pi*r^2 dr
        S_kin = 4 * np.pi * np.trapz(r**2 * 0.5 * fg * gp**2, r)

        print(f"    {label}: E_tot = {E_tot:.6f},  S_kin = {S_kin:.6f}")

    # --- 2e: Koide Q verification ---
    print("\n  --- 2e: Koide Q verification ---")
    r, g = solve_soliton(g0_e)
    A_e, _, _ = tail_amplitude(r, g)
    r, g = solve_soliton(g0_mu)
    A_mu, _, _ = tail_amplitude(r, g)
    r, g = solve_soliton(g0_tau)
    A_tau, _, _ = tail_amplitude(r, g)

    m_e = A_e**4
    m_mu = A_mu**4
    m_tau = A_tau**4

    Q = (np.sqrt(m_e) + np.sqrt(m_mu) + np.sqrt(m_tau))**2 / (m_e + m_mu + m_tau)
    print(f"  Q = (sum sqrt(m))^2 / (sum m) = {Q:.10f}")
    print(f"  Q exact (Koide) = 1.5")
    print(f"  Deviation = {abs(Q - 1.5)/1.5*100:.6f}%")
    print(f"  [PASS]" if abs(Q - 1.5) < 0.001 else f"  [FAIL]")

    results['Q'] = Q
    return results


# ================================================================
#  STEP 3: Selection principle candidates
# ================================================================

def step3_selection_principles(g0_e, g0_mu, g0_tau_koide, A_e):
    """Test candidate selection principles for g0_tau."""
    print("\n" + "="*60)
    print("  STEP 3: Selection Principle Candidates")
    print("="*60)

    if g0_tau_koide is None:
        print("  [SKIPPED: g0_tau not found]")
        return

    # --- 3a: Variational principle on total soliton energy ---
    print("\n  --- 3a: Minimum total soliton energy ---")
    # Is g0_tau at a minimum of E_tot(g0)?
    g0_range = np.linspace(g0_mu + 0.1, g0_tau_koide + 2.0, 100)
    E_vals = []
    for g0 in g0_range:
        r, g = solve_soliton(g0)
        if r is None or r[-1] < 250:
            E_vals.append(np.nan)
            continue
        dr = r[1] - r[0]
        gp = np.gradient(g, dr)
        fg = np.array([fk(gi) for gi in g])
        E_dens = 0.5 * fg * gp**2 + np.array([V(gi) for gi in g])
        E_tot = 4 * np.pi * np.trapz(r**2 * E_dens, r)
        E_vals.append(E_tot)
    E_vals = np.array(E_vals)

    valid = ~np.isnan(E_vals)
    if np.sum(valid) > 5:
        idx_min = np.nanargmin(E_vals[valid])
        g0_E_min = g0_range[valid][idx_min]
        print(f"  E_tot minimum at g0 = {g0_E_min:.4f}")
        print(f"  g0_tau (Koide)      = {g0_tau_koide:.4f}")
        print(f"  Match: {'YES' if abs(g0_E_min - g0_tau_koide) < 0.1 else 'NO'}")
        print(f"  (Delta = {abs(g0_E_min - g0_tau_koide):.4f})")

    # --- 3b: Action quantization ---
    print("\n  --- 3b: Action quantization S_n = n * S_1 ---")
    actions = {}
    for label, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau_koide)]:
        r, g = solve_soliton(g0)
        if r is None:
            continue
        dr = r[1] - r[0]
        gp = np.gradient(g, dr)
        fg = np.array([fk(gi) for gi in g])
        S = 4 * np.pi * np.trapz(r**2 * 0.5 * fg * gp**2, r)
        actions[label] = S
        print(f"    S({label}) = {S:.6f}")

    if 'e' in actions and 'mu' in actions and 'tau' in actions:
        Se, Smu, Stau = actions['e'], actions['mu'], actions['tau']
        print(f"    S(mu)/S(e) = {Smu/Se:.4f}")
        print(f"    S(tau)/S(e) = {Stau/Se:.4f}")
        print(f"    S(tau)/S(mu) = {Stau/Smu:.4f}")

    # --- 3c: Phase quantization (Bohr-Sommerfeld) ---
    print("\n  --- 3c: Phase spacing quantization ---")
    g0_range2 = np.linspace(0.82, 5.0, 300)
    phase_data = []
    A_data = []

    for g0 in g0_range2:
        r, g = solve_soliton(g0)
        if r is None or r[-1] < 250:
            phase_data.append(np.nan)
            A_data.append(np.nan)
            continue
        A, ph, _ = tail_amplitude(r, g)
        phase_data.append(ph)
        A_data.append(A)

    phase_data = np.array(phase_data)
    A_data = np.array(A_data)
    valid = ~np.isnan(phase_data) & ~np.isnan(A_data)

    if np.sum(valid) > 20:
        phases_uw = np.unwrap(phase_data[valid])
        g0_v = g0_range2[valid]
        ph_interp = interp1d(g0_v, phases_uw, kind='cubic', fill_value='extrapolate')

        ph_e = ph_interp(g0_e)
        ph_mu = ph_interp(g0_mu)
        ph_tau = ph_interp(g0_tau_koide)

        dp_12 = ph_mu - ph_e
        dp_23 = ph_tau - ph_mu
        dp_13 = ph_tau - ph_e

        print(f"    Phase(e) = {ph_e:.6f}")
        print(f"    Phase(mu) = {ph_mu:.6f}")
        print(f"    Phase(tau) = {ph_tau:.6f}")
        print(f"    Delta(e->mu) = {dp_12:.6f}")
        print(f"    Delta(mu->tau) = {dp_23:.6f}")
        print(f"    Delta(e->tau) = {dp_13:.6f}")
        print(f"    Delta(mu->tau)/Delta(e->mu) = {dp_23/dp_12:.4f}" if abs(dp_12) > 0.01 else "    dp_12 ~ 0")

        # Check if dp_12 ~ n * pi for some small n
        for n in range(1, 7):
            for base in [np.pi, 2*np.pi/3, np.pi/2, np.pi/3]:
                if abs(dp_12 - n * base) < 0.2:
                    print(f"    ** Delta(e->mu) ~ {n}*{base/np.pi:.3f}*pi (diff = {abs(dp_12 - n*base):.4f})")

    # --- 3d: Zero-crossing count (topological quantum number) ---
    print("\n  --- 3d: Soliton zero-crossing count (g=1 crossings) ---")
    for label, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau_koide)]:
        r, g = solve_soliton(g0)
        if r is None:
            continue
        # Count crossings of g=1
        crossings = 0
        for i in range(len(g)-1):
            if (g[i]-1)*(g[i+1]-1) < 0:
                crossings += 1
        # Count nodes of (g-1) in the "inner" region r < 50
        inner_mask = r < 50
        inner_crossings = 0
        gi = g[inner_mask]
        for i in range(len(gi)-1):
            if (gi[i]-1)*(gi[i+1]-1) < 0:
                inner_crossings += 1
        print(f"    {label}: total crossings = {crossings}, inner (r<50) = {inner_crossings}")

    # --- 3e: Soliton charge / topological invariant ---
    print("\n  --- 3e: Soliton topological charge ---")
    for label, g0 in [("e", g0_e), ("mu", g0_mu), ("tau", g0_tau_koide)]:
        r, g = solve_soliton(g0)
        if r is None:
            continue
        # Topological charge: Q_top = g(0) - g(inf) = g0 - 1 (for g -> 1)
        Q_top = g0 - 1
        # Winding: integral of d(ln g)
        dg = np.diff(g)
        g_mid = 0.5 * (g[:-1] + g[1:])
        g_mid[g_mid <= 0] = 1e-30
        winding = np.sum(dg / g_mid)
        print(f"    {label}: g0-1 = {Q_top:.6f}, winding = {winding:.6f}")

    # --- 3f: Log-spacing (geometric sequence) ---
    print("\n  --- 3f: Log-spacing ---")
    lng_e = np.log(g0_e)
    lng_mu = np.log(g0_mu)
    lng_tau = np.log(g0_tau_koide)
    d1 = lng_mu - lng_e
    d2 = lng_tau - lng_mu
    print(f"    ln(g0_e)  = {lng_e:.6f}")
    print(f"    ln(g0_mu) = {lng_mu:.6f}")
    print(f"    ln(g0_tau)= {lng_tau:.6f}")
    print(f"    d1 = ln(g0_mu/g0_e)  = {d1:.6f} (= ln(phi) = {np.log(PHI):.6f})")
    print(f"    d2 = ln(g0_tau/g0_mu) = {d2:.6f}")
    print(f"    d2/d1 = {d2/d1:.6f}")
    print(f"    d2 - d1 = {d2-d1:.6f}")

    # Check if g0_tau = g0_e * phi^k for some k
    k_tau = np.log(g0_tau_koide/g0_e) / np.log(PHI)
    print(f"    g0_tau/g0_e = phi^{k_tau:.6f}")
    print(f"    Nearest integer k: {round(k_tau)}, fractional part: {k_tau - round(k_tau):.6f}")

    # --- 3g: Alpha_K dependence ---
    print("\n  --- 3g: Dependence on kinetic coupling ---")
    print(f"    f(g) = 1 + 2*alpha*ln(g), alpha = {ALPHA}")
    print(f"    f(g0_e)   = {fk(g0_e):.6f}")
    print(f"    f(g0_mu)  = {fk(g0_mu):.6f}")
    print(f"    f(g0_tau) = {fk(g0_tau_koide):.6f}")

    # f(g0_tau)/f(g0_e)
    ratio_f = fk(g0_tau_koide) / fk(g0_e)
    print(f"    f(tau)/f(e) = {ratio_f:.6f}")
    print(f"    f(tau)/f(mu) = {fk(g0_tau_koide)/fk(g0_mu):.6f}")


# ================================================================
#  STEP 4: Comprehensive scan for g0_tau selection
# ================================================================

def step4_selection_scan(g0_e, A_e):
    """Scan for selection principles that give r_31 ~ 3477."""
    print("\n" + "="*60)
    print("  STEP 4: Systematic selection principle scan")
    print("="*60)

    candidates = {}

    # Candidate 1: phi^k * g0_e
    print("\n  --- Candidate 1: phi^k * g0_e ---")
    for k_name, k_val in [("phi^2", PHI**2), ("phi^3", PHI**3),
                           ("phi^(5/2)", PHI**2.5), ("phi^(3/2)", PHI**1.5),
                           ("phi^(ln3/ln_phi)", PHI**(np.log(3)/np.log(PHI)))]:
        g0_try = g0_e * k_val
        r31 = mass_ratio(g0_e, g0_try)
        print(f"    {k_name}: g0_tau = {g0_try:.4f}, r_31 = {r31:.1f}" if not np.isnan(r31) else f"    {k_name}: g0_tau = {g0_try:.4f}, r_31 = NaN")
        if not np.isnan(r31):
            candidates[k_name] = (g0_try, r31, abs(r31 - R31_KOIDE)/R31_KOIDE)

    # Candidate 2: arithmetic progressions in g0
    print("\n  --- Candidate 2: arithmetic progression ---")
    delta_12 = g0_e * (PHI - 1)  # = g0_mu - g0_e
    g0_arith = [g0_e + n * delta_12 for n in [2, 3, 4, 5]]
    for n, g0_try in zip([2,3,4,5], g0_arith):
        r31 = mass_ratio(g0_e, g0_try)
        print(f"    n={n}: g0 = {g0_try:.4f}, r_31 = {r31:.1f}" if not np.isnan(r31) else f"    n={n}: g0 = {g0_try:.4f}, ODE diverged")

    # Candidate 3: g0_tau such that ln(g0_tau) = a + b*n^2
    print("\n  --- Candidate 3: quadratic in ln(g0) ---")
    # ln(g0_e) = a, ln(g0_mu) = a + b, ln(g0_tau) = a + 4*b
    lng_e = np.log(g0_e)
    b_log = np.log(PHI)  # since g0_mu/g0_e = phi
    # ln(g0_tau) = ln(g0_e) + 4*b = lng_e + 4*ln(phi)
    g0_quad = np.exp(lng_e + 4 * b_log)
    r31_quad = mass_ratio(g0_e, g0_quad)
    print(f"    ln(g0_n) = a + b*n^2: g0_tau = {g0_quad:.4f}, r_31 = {r31_quad:.1f}" if not np.isnan(r31_quad) else f"    Diverged")

    # Candidate 4: g0_tau from Koide eigenvalue condition
    # Koide: (sum sqrt(m))^2 / (sum m) = 3/2
    # In terms of amplitudes: A^2 replaces sqrt(m)
    # (A_e^2 + A_mu^2 + A_tau^2)^2 / (A_e^4 + A_mu^4 + A_tau^4) = 3/2
    print("\n  --- Candidate 4: Direct Koide on A_tail ---")
    r_e, g_e = solve_soliton(g0_e)
    A_e_val, _, _ = tail_amplitude(r_e, g_e)
    r_mu, g_mu = solve_soliton(PHI * g0_e)
    A_mu_val, _, _ = tail_amplitude(r_mu, g_mu)

    # Find A_tau such that Koide holds
    # (A_e^2 + A_mu^2 + A_tau^2)^2 = 1.5 * (A_e^4 + A_mu^4 + A_tau^4)
    # Let x = A_tau^2
    # (A_e^2 + A_mu^2 + x)^2 = 1.5 * (A_e^4 + A_mu^4 + x^2)
    a2 = A_e_val**2
    b2 = A_mu_val**2
    # (a2 + b2 + x)^2 = 1.5*(a2^2 + b2^2 + x^2)
    # a2^2 + b2^2 + x^2 + 2*a2*b2 + 2*a2*x + 2*b2*x = 1.5*a2^2 + 1.5*b2^2 + 1.5*x^2
    # -0.5*x^2 + 2*(a2+b2)*x + (a2+b2)^2 - 1.5*(a2^2+b2^2) = 0
    # -0.5*x^2 + 2*(a2+b2)*x + a2^2 + 2*a2*b2 + b2^2 - 1.5*a2^2 - 1.5*b2^2 = 0
    # -0.5*x^2 + 2*(a2+b2)*x - 0.5*a2^2 + 2*a2*b2 - 0.5*b2^2 = 0
    # x^2 - 4*(a2+b2)*x + a2^2 - 4*a2*b2 + b2^2 = 0
    # x^2 - 4*(a2+b2)*x + (a2-b2)^2 - 2*a2*b2 = 0 -- hmm let me just solve numerically

    coeff_a = 1.0
    coeff_b = -4*(a2 + b2)
    coeff_c = (a2 + b2)**2 - 1.5*(a2**2 + b2**2)  # wait, let me redo
    # Going back: -0.5x^2 + 2(a2+b2)x + (a2+b2)^2 - 1.5(a2^2+b2^2) = 0
    # multiply by -2: x^2 - 4(a2+b2)x - 2(a2+b2)^2 + 3(a2^2+b2^2) = 0
    # = x^2 - 4(a2+b2)x + 3(a2^2+b2^2) - 2(a2^2+2*a2*b2+b2^2)
    # = x^2 - 4(a2+b2)x + a2^2 - 4*a2*b2 + b2^2
    # = x^2 - 4(a2+b2)x + (a2-b2)^2 - 2*a2*b2  -- hmm this isn't factoring nicely

    A_coef = 1.0
    B_coef = -4*(a2+b2)
    C_coef = a2**2 - 4*a2*b2 + b2**2
    disc = B_coef**2 - 4*A_coef*C_coef
    if disc >= 0:
        x1 = (-B_coef + np.sqrt(disc)) / (2*A_coef)
        x2 = (-B_coef - np.sqrt(disc)) / (2*A_coef)
        for x in [x1, x2]:
            if x > 0:
                A_tau_koide = np.sqrt(x)
                r31_from_A = (A_tau_koide / A_e_val)**4
                print(f"    A_tau (Koide) = {A_tau_koide:.8f}, r_31 = {r31_from_A:.2f}")

    # Best candidates summary
    print("\n  --- Best candidates summary ---")
    if candidates:
        best = min(candidates.items(), key=lambda kv: kv[1][2])
        print(f"  Best: {best[0]}, g0={best[1][0]:.4f}, r_31={best[1][1]:.1f}, " +
              f"dev={best[1][2]*100:.1f}%")
        for name, (g0, r31, dev) in sorted(candidates.items(), key=lambda kv: kv[1][2]):
            print(f"    {name}: dev = {dev*100:.1f}%")


# ================================================================
#  MAIN
# ================================================================

def main():
    print("="*60)
    print("  O-K1: Search for Koide Q=3/2 Mechanism in TGP")
    print("  (Corrected version with phi-FP parameters)")
    print("="*60)

    # Step 0
    g0_e, g0_mu, A_e, A_mu, phase_e, phase_mu = step0_baseline()

    # Step 1: find g0_tau
    g0_tau, A_tau, phase_tau = step1_find_g0_tau(g0_e, A_e)

    # Step 2: analyze
    results = step2_analyze(g0_e, g0_mu, g0_tau, phase_e, phase_mu, phase_tau)

    # Step 3: selection principles
    step3_selection_principles(g0_e, g0_mu, g0_tau, A_e)

    # Step 4: systematic scan
    step4_selection_scan(g0_e, A_e)

    # Summary
    print("\n" + "="*60)
    print("  O-K1 FINAL SUMMARY")
    print("="*60)

    if g0_tau is not None:
        ratio = g0_tau / g0_e
        k_phi = np.log(ratio) / np.log(PHI)
        print(f"  g0_e   = {g0_e:.10f}")
        print(f"  g0_mu  = {g0_mu:.10f} (= phi * g0_e)")
        print(f"  g0_tau = {g0_tau:.10f} (for Koide Q=3/2)")
        print(f"  g0_tau/g0_e = {ratio:.6f} = phi^{k_phi:.4f}")
        if results and 'Q' in results:
            print(f"  Q (Koide) = {results['Q']:.10f}")
        print(f"\n  Status: g0_tau FOUND numerically.")
        print(f"  Selection principle: OPEN (none of the simple principles match)")
    else:
        print(f"  g0_tau NOT FOUND in scanned range.")
        print(f"  The soliton ODE may not reach r_31 = {R31_KOIDE:.0f}")
        print(f"  This suggests Koide requires ADDITIONAL physics")
        print(f"  (multi-soliton, substrate corrections, or higher-order terms)")

    return results


if __name__ == '__main__':
    main()
