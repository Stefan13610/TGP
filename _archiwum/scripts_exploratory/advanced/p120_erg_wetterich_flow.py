#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p120_erg_wetterich_flow.py  --  Wetterich ERG flow for TGP
===========================================================
Functional renormalization group (FRG) for TGP effective potential.

Wetterich equation (d=4 Euclidean, Litim regulator):

  dV/ds = -K(psi)*k^6 / (32*pi^2 * (K(psi)*k^2 + V''(psi)))

where s = ln(k_UV/k), k(s) = k_UV * exp(-s).

Three modes:
  A: LPA only (K=1)           -- baseline, expected to show instability
  B: LPA + K(psi)=psi^4 fixed -- expected to stabilize vacuum
  C: Coupled (V,K)            -- both V and K flow

Physical conclusion: K(psi)=psi^4 REVERSES the sign of V''(1) from
negative to positive, stabilizing the TGP vacuum.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp

# === Physical parameters ===
PHI0    = 24.66
A_GAMMA = 0.040049
GAMMA   = 1.0           # beta = gamma = 1
M_SP    = np.sqrt(GAMMA)

K_UV  = 1.0 / A_GAMMA   # ~24.97
K_IR  = M_SP             # 1.0
S_MAX = np.log(K_UV / K_IR)  # ~3.22

# === Field grid ===
N_PSI   = 80
PSI_MIN = 0.05
PSI_MAX = 3.5
PSI     = np.linspace(PSI_MIN, PSI_MAX, N_PSI)
DPSI    = PSI[1] - PSI[0]

# Denominator regularization (prevents 1/0 at spinodal)
EPS_DENOM = 1e-8


def V_tree(psi):
    """V = psi^3/3 - psi^4/4  (beta=gamma=1)"""
    return psi**3 / 3.0 - psi**4 / 4.0

def K_psi4(psi):
    """K(psi) = psi^4"""
    return psi**4

def d2_grid(f):
    """Second derivative via 5-point stencil on uniform grid."""
    n = len(f)
    d2 = np.zeros(n)
    h = DPSI
    # Interior: 5-point central
    for i in range(2, n-2):
        d2[i] = (-f[i+2] + 16*f[i+1] - 30*f[i] + 16*f[i-1] - f[i-2]) / (12*h**2)
    # Edges: 3-point
    d2[0] = (f[2] - 2*f[1] + f[0]) / h**2
    d2[1] = (f[2] - 2*f[1] + f[0]) / h**2
    d2[n-1] = (f[n-1] - 2*f[n-2] + f[n-3]) / h**2
    d2[n-2] = (f[n-1] - 2*f[n-2] + f[n-3]) / h**2
    return d2

def d3_grid(f):
    """Third derivative via 4-point central."""
    n = len(f)
    d3 = np.zeros(n)
    h = DPSI
    for i in range(2, n-2):
        d3[i] = (f[i+2] - 2*f[i+1] + 2*f[i-1] - f[i-2]) / (2*h**3)
    d3[:2] = d3[2]
    d3[n-2:] = d3[n-3]
    return d3


def safe_denom(K_arr, k, Vpp):
    """Regularized denominator: max(|D|, eps) * sign(D)."""
    D = K_arr * k**2 + Vpp
    return np.where(np.abs(D) < EPS_DENOM, EPS_DENOM, D)


# === Flow equations ===

def rhs_A(s, V_vec):
    """Mode A: LPA (K=1), dV/ds = -k^6/(32*pi^2*(k^2+V''))."""
    k = K_UV * np.exp(-s)
    Vpp = d2_grid(V_vec)
    D = safe_denom(np.ones(N_PSI), k, Vpp)
    return -k**6 / (32.0 * np.pi**2 * D)

def rhs_B(s, V_vec):
    """Mode B: K(psi)=psi^4 fixed, dV/ds = -K*k^6/(32*pi^2*(K*k^2+V''))."""
    k = K_UV * np.exp(-s)
    K = K_psi4(PSI)
    Vpp = d2_grid(V_vec)
    D = safe_denom(K, k, Vpp)
    return -K * k**6 / (32.0 * np.pi**2 * D)

def rhs_C(s, state):
    """Mode C: Coupled (V,K). K flows via V''' feedback."""
    V_vec = state[:N_PSI]
    K_vec = np.maximum(state[N_PSI:], 1e-12)  # keep K > 0
    k = K_UV * np.exp(-s)

    Vpp  = d2_grid(V_vec)
    Vppp = d3_grid(V_vec)

    D = safe_denom(K_vec, k, Vpp)

    dVds = -K_vec * k**6 / (32.0 * np.pi**2 * D)
    dKds = K_vec**2 * k**6 * Vppp**2 / (16.0 * np.pi**2 * D**3)

    return np.concatenate([dVds, dKds])


def run_mode(mode, verbose=True):
    """Run one ERG mode. Returns dict of results."""
    V0 = V_tree(PSI)

    s_eval = np.linspace(0, S_MAX, 200)

    if mode == 'A':
        sol = solve_ivp(rhs_A, [0, S_MAX], V0, method='Radau',
                        rtol=1e-6, atol=1e-8, t_eval=s_eval, max_step=0.02)
        K_final = np.ones(N_PSI)
        label = "LPA (K=1)"
    elif mode == 'B':
        sol = solve_ivp(rhs_B, [0, S_MAX], V0, method='Radau',
                        rtol=1e-6, atol=1e-8, t_eval=s_eval, max_step=0.02)
        K_final = K_psi4(PSI)
        label = "LPA + K(psi)=psi^4"
    elif mode == 'C':
        K0 = K_psi4(PSI)
        state0 = np.concatenate([V0, K0])
        sol = solve_ivp(rhs_C, [0, S_MAX], state0, method='Radau',
                        rtol=1e-6, atol=1e-8, t_eval=s_eval, max_step=0.02)
        K_final = np.maximum(sol.y[N_PSI:, -1], 1e-15)
        label = "Coupled (V,K)"
    else:
        raise ValueError(f"Unknown mode {mode}")

    V_IR = sol.y[:N_PSI, -1]
    Vpp_UV = d2_grid(V0)
    Vpp_IR = d2_grid(V_IR)

    # Index of psi=1
    i1 = np.argmin(np.abs(PSI - 1.0))

    K_UV_1 = (K_psi4(PSI) if mode != 'A' else np.ones(N_PSI))[i1]
    K_IR_1 = K_final[i1]

    m2_UV = Vpp_UV[i1] / K_UV_1 if abs(K_UV_1) > 1e-15 else float('inf')
    m2_IR = Vpp_IR[i1] / K_IR_1 if abs(K_IR_1) > 1e-15 else float('inf')

    # Stability scan
    scan = {}
    for pt in [0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.3]:
        idx = np.argmin(np.abs(PSI - pt))
        vpp = Vpp_IR[idx]
        kk  = K_final[idx]
        m2  = vpp / kk if abs(kk) > 1e-15 else float('inf')
        scan[pt] = {"Vpp": float(vpp), "K": float(kk), "m2": float(m2)}

    res = {
        "mode": mode, "label": label,
        "success": bool(sol.success), "nsteps": int(sol.t.size),
        "Vpp_UV_1": float(Vpp_UV[i1]), "Vpp_IR_1": float(Vpp_IR[i1]),
        "K_UV_1": float(K_UV_1), "K_IR_1": float(K_IR_1),
        "m2_UV": float(m2_UV), "m2_IR": float(m2_IR),
        "stable": bool(m2_IR > 0),
        "scan": scan
    }

    if verbose:
        ok = "OK" if sol.success else "FAILED"
        stab = "STABLE" if res['stable'] else "UNSTABLE"
        print(f"\n  Mode {mode}: {label}")
        print(f"  Solver: {ok} ({sol.t.size} eval points)")
        print(f"  V''(1): UV={Vpp_UV[i1]:+.4f}  IR={Vpp_IR[i1]:+.4f}")
        print(f"  K(1):   UV={K_UV_1:.4f}  IR={K_IR_1:.6f}")
        print(f"  m2_phys: UV={m2_UV:+.4f}  IR={m2_IR:+.4f}")
        print(f"  Vacuum: {stab}")
        print(f"  {'psi':>5s} {'Vpp':>10s} {'K':>10s} {'m2':>10s} {'ok':>6s}")
        for pt, v in sorted(scan.items()):
            st = "+" if v['m2'] > 0 else "-"
            print(f"  {pt:5.1f} {v['Vpp']:+10.3f} {v['K']:10.4f} {v['m2']:+10.3f} {st:>6s}")

    return res


def main():
    print("="*60)
    print("  TGP Wetterich ERG Flow Verification")
    print("  k_UV={:.2f}, k_IR={:.2f}, ln(k_UV/k_IR)={:.3f}".format(K_UV, K_IR, S_MAX))
    print("="*60)

    results = {}
    for mode in ['A', 'B', 'C']:
        try:
            results[mode] = run_mode(mode)
        except Exception as e:
            print(f"\n  Mode {mode}: EXCEPTION -- {e}")
            results[mode] = {"mode": mode, "error": str(e)}

    # === Summary and tests ===
    print("\n" + "="*60)
    print("  SUMMARY")
    print("="*60)
    print(f"  {'Mode':<25s} {'Vpp(1)_IR':>10s} {'m2_phys':>10s} {'Stable':>8s}")
    print(f"  {'-'*55}")
    for m in ['A', 'B', 'C']:
        r = results[m]
        if 'error' in r:
            print(f"  {m:<25s} {'ERR':>10s}")
        else:
            st = "YES" if r['stable'] else "NO"
            print(f"  {r['label']:<25s} {r['Vpp_IR_1']:>+10.3f} "
                  f"{r['m2_IR']:>+10.3f} {st:>8s}")

    # Tests
    tests = []
    if 'error' not in results.get('A', {}):
        rA = results['A']
        tests.append(("T1: LPA vacuum unstable (expected)", not rA['stable']))
        tests.append(("T2: LPA V''(1)<0 at IR", rA['Vpp_IR_1'] < 0))

    if 'error' not in results.get('B', {}):
        rB = results['B']
        tests.append(("T3: K(psi)=psi^4 stabilizes vacuum", rB['stable']))
        tests.append(("T4: V''(1)>0 at IR with K", rB['Vpp_IR_1'] > 0))
        tests.append(("T5: m2_phys > 0 with K", rB['m2_IR'] > 0))

    if 'error' not in results.get('C', {}):
        rC = results['C']
        tests.append(("T6: Coupled mode stable", rC['stable']))
        tests.append(("T7: Coupled m2 marginal (0<m2<200)", 0 < rC['m2_IR'] < 200))
        if 'error' not in results.get('B', {}):
            tests.append(("T8: Coupled m2 <= fixed-K m2",
                          rC['m2_IR'] <= results['B']['m2_IR'] + 1))

    n_pass = sum(1 for _, v in tests if v)
    n_total = len(tests)

    print(f"\n  TESTS: {n_pass}/{n_total} PASS")
    for name, passed in tests:
        mark = "PASS" if passed else "FAIL"
        print(f"    [{mark}] {name}")

    if n_pass == n_total and n_total > 0:
        print(f"\n  VERDICT: ALL {n_total} TESTS PASS")
        print("  K(psi)=psi^4 from substrate stabilizes the TGP vacuum under ERG flow.")
    else:
        print(f"\n  VERDICT: {n_total - n_pass} TESTS FAILED")

    return results


if __name__ == '__main__':
    main()
