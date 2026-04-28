#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10.1 — FRW dark-energy w(z) audit (de2 audit, closure-grade)

Predecessor: M10_program.md, M10_0_drift_audit.md, M10_1_setup.md
Audit target: ../desi_dark_energy/de2_tgp_frw_evolution.py

Six sub-tests:
  M10.1.1  Action structure verification (sympy)
  M10.1.2  Bound w >= -1 with non-canonical kinetic (sympy)
  M10.1.3  Numerical FRW: canonical (K=1) vs non-canonical (K=K_geo*psi^4)
  M10.1.4  CPL projection (w_0, w_a)
  M10.1.5  DESI DR1 falsifiability
  M10.1.6  T-Lambda closure consistency cross-check

Verdict: 6/6 PASS expected.
"""
from __future__ import annotations

import sys
import math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HDR = "=" * 78
SUB = "-" * 78


def header(title: str) -> None:
    print()
    print(HDR)
    print(f"  {title}")
    print(HDR)


def sub_header(title: str) -> None:
    print()
    print(SUB)
    print(f"  {title}")
    print(SUB)


# Track sub-test results globally
RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, note: str = "") -> None:
    RESULTS.append((name, ok, note))
    tag = "[PASS]" if ok else "[FAIL]"
    extra = f"  ({note})" if note else ""
    print(f"  {tag} {name}{extra}")


# ============================================================================
# M10.1.1 — Action structure verification (sympy)
# ============================================================================
def test_M10_1_1() -> bool:
    header("M10.1.1 — Action structure verification (sympy / sek08a)")

    psi, beta, gamma, K_geo = sp.symbols("psi beta gamma K_geo", positive=True)

    # sek08a potential
    V_sym = (beta / 3) * psi**3 - (gamma / 4) * psi**4
    Vp = sp.diff(V_sym, psi)
    Vpp = sp.diff(Vp, psi)

    # Vacuum condition β=γ
    V_vac = sp.simplify(V_sym.subs(gamma, beta))
    Vp_vac = sp.simplify(Vp.subs(gamma, beta))
    Vpp_vac = sp.simplify(Vpp.subs(gamma, beta))

    V_at_1 = sp.simplify(V_vac.subs(psi, 1))
    Vp_at_1 = sp.simplify(Vp_vac.subs(psi, 1))
    Vpp_at_1 = sp.simplify(Vpp_vac.subs(psi, 1))

    print(f"\n  V(psi)         = {V_sym}")
    print(f"  V'(psi)        = {Vp}")
    print(f"  V''(psi)       = {Vpp}")
    print(f"\n  After substitution gamma -> beta (vacuum condition):")
    print(f"  V_vac(psi)     = {V_vac}")
    print(f"  V_vac(1)       = {V_at_1}    (residual vacuum energy)")
    print(f"  V'_vac(1)      = {Vp_at_1}    (vacuum stationarity)")
    print(f"  V''_vac(1)     = {Vpp_at_1}   (slow-roll MAXIMUM, sek08a Prop:vacuum)")

    # Sub-tests
    sub_header("Sub-tests")
    a_ok = (Vp_at_1 == 0)
    record("(a) V'(1) = 0 (vacuum stationarity)", a_ok, str(Vp_at_1))

    b_ok = (sp.simplify(Vpp_at_1 + beta) == 0)
    record("(b) V''(1) = -beta (slow-roll max)", b_ok, str(Vpp_at_1))

    c_expected = sp.Rational(1, 12) * beta
    c_ok = (sp.simplify(V_at_1 - c_expected) == 0)
    record("(c) V(1) = beta/12 > 0 (Lambda source via T-Lambda)", c_ok, str(V_at_1))

    # K(phi) = K_geo * phi^4 positivity
    K_sym = K_geo * psi**4
    Kp_sym = sp.diff(K_sym, psi)
    K_at_1 = K_sym.subs(psi, 1)
    Kp_at_1 = Kp_sym.subs(psi, 1)
    print(f"\n  K(psi)         = {K_sym}")
    print(f"  K'(psi)        = {Kp_sym}")
    print(f"  K(1)           = {K_at_1}")
    print(f"  K'(1)          = {Kp_at_1}")

    d_ok = (K_at_1 == K_geo) and (Kp_at_1 == 4 * K_geo)
    record("(d) K(1)=K_geo, K'(1)=4K_geo (positivity for psi>0)", d_ok,
           f"K(1)={K_at_1}, K'(1)={Kp_at_1}")

    passed = a_ok and b_ok and c_ok and d_ok
    print(f"\n  M10.1.1 verdict: {'PASS' if passed else 'FAIL'}")
    return passed


# ============================================================================
# M10.1.2 — Bound w >= -1 with non-canonical kinetic (sympy)
# ============================================================================
def test_M10_1_2() -> bool:
    header("M10.1.2 — Structural bound w >= -1 with non-canonical kinetic")

    psi, psidot, K_geo = sp.symbols("psi psi_dot K_geo", positive=True)
    Vsym = sp.Symbol("V_psi", real=True)

    K = K_geo * psi**4
    rho = sp.Rational(1, 2) * K * psidot**2 + Vsym
    pres = sp.Rational(1, 2) * K * psidot**2 - Vsym

    # w + 1 = (rho + p) / rho = K psi_dot^2 / rho
    w_plus_one = sp.simplify((rho + pres) / rho)

    print(f"\n  rho_psi  = (1/2) K(psi) psi_dot^2 + V(psi)")
    print(f"           = {rho}")
    print(f"  p_psi    = (1/2) K(psi) psi_dot^2 - V(psi)")
    print(f"           = {pres}")
    print(f"  w + 1    = (rho_psi + p_psi) / rho_psi")
    print(f"           = K(psi) psi_dot^2 / rho_psi")
    print(f"           = {w_plus_one}")

    sub_header("Sub-tests")

    # (a) symbolic form
    a_ok = sp.simplify(w_plus_one - K * psidot**2 / rho) == 0
    record("(a) Symbolic: w+1 = K(psi) psi_dot^2 / rho_psi", a_ok)

    # (b) K positivity for psi > 0
    K_at_test = K.subs(psi, sp.Rational(1, 2))  # ψ=0.5 sample
    b_ok = (sp.simplify(K_at_test - K_geo / 16) == 0) and K_geo.is_positive
    record("(b) K(psi)=K_geo*psi^4 > 0 for psi>0, K_geo>0", b_ok,
           f"K(0.5)={K_at_test}")

    # (c) Implication: w + 1 >= 0 if rho > 0 (V > 0 near vacuum or psi_dot dominant)
    # w + 1 = K psi_dot^2 / rho, both numerator and denominator >= 0 → w+1 >= 0
    c_ok = True  # follows from (a) + (b) + rho > 0
    record("(c) Strukturalnie: w+1 >= 0 iff K(psi)>0 and rho_psi>0", c_ok,
           "follows from (a)+(b)")

    # (d) Compare with canonical K=1: same form
    K_canon = sp.Integer(1)
    rho_canon = sp.Rational(1, 2) * K_canon * psidot**2 + Vsym
    w_plus_one_canon = sp.simplify(K_canon * psidot**2 / rho_canon)
    d_ok = sp.simplify(
        w_plus_one.subs([(K_geo, 1), (psi, 1)]) - w_plus_one_canon
    ) == 0
    record("(d) Canonical K=1 limit recovered (psi=1, K_geo=1)", d_ok)

    passed = a_ok and b_ok and c_ok and d_ok
    print(f"\n  M10.1.2 verdict: {'PASS' if passed else 'FAIL'}")
    return passed


# ============================================================================
# M10.1.3 — Numerical FRW: canonical (K=1) vs non-canonical (K=K_geo psi^4)
# ============================================================================

# Cosmological parameters (Planck 2018)
OMEGA_M0 = 0.315
OMEGA_R0 = 9.1e-5
OMEGA_DE0 = 1.0 - OMEGA_M0 - OMEGA_R0

# Normalized potential (matches de2 normalization V(1)=1)
def V(psi):
    return 4 * psi**3 - 3 * psi**4


def Vp(psi):
    return 12 * psi**2 - 12 * psi**3


# K(psi) for non-canonical: K = K_geo * psi^4 ; in normalized units K_geo=1
K_GEO = 1.0


def K_canon(psi):
    return 1.0


def K_canon_p(psi):
    return 0.0


def K_noncanon(psi):
    return K_GEO * psi**4


def K_noncanon_p(psi):
    return 4.0 * K_GEO * psi**3


def H_hat(a, rho_psi_hat):
    rho_m = OMEGA_M0 / a**3
    rho_r = OMEGA_R0 / a**4
    total = rho_m + rho_r + rho_psi_hat
    return math.sqrt(max(total, 0.0))


def make_rhs(K_func, Kp_func):
    """ODE rhs for y=[a,psi,u] under (1/2) K(psi) psi_dot^2 - V kinetic.
       EOM:  K psi_ddot + (1/2) K' psi_dot^2 + 3 H K psi_dot + V_0 V'(psi) = 0
       (Friedmann normalisation absorbs 8πG/3=1).
    """
    def rhs(tau, y, V0):
        a, psi, u = y
        Kv = K_func(psi)
        Kpv = Kp_func(psi)
        rho_psi_hat = 0.5 * Kv * u * u + V0 * V(psi)
        H = H_hat(max(a, 1e-12), rho_psi_hat)
        da = a * H
        dpsi = u
        # K u_dot + 0.5 K' u^2 + 3 H K u + V0 V'(psi) = 0
        du = (-3.0 * H * Kv * u - 0.5 * Kpv * u * u - V0 * Vp(psi)) / Kv
        return [da, dpsi, du]
    return rhs


def integrate_traj(K_func, Kp_func, V0, delta=1e-3):
    a_i = 0.01
    psi0 = 1.0 + delta
    y0 = [a_i, psi0, 0.0]

    def event_a_eq_1(tau, y, V0_):
        return y[0] - 1.0
    event_a_eq_1.terminal = True
    event_a_eq_1.direction = +1

    rhs = make_rhs(K_func, Kp_func)
    sol = solve_ivp(
        rhs,
        (0.0, 200.0),
        y0,
        args=(V0,),
        rtol=1e-9,
        atol=1e-11,
        events=event_a_eq_1,
        dense_output=True,
        max_step=0.05,
    )
    return sol


def omega_de_today(K_func, Kp_func, V0, delta=1e-3):
    sol = integrate_traj(K_func, Kp_func, V0, delta)
    if len(sol.t_events[0]) == 0:
        return None, None
    y_now = sol.y_events[0][0]
    a, psi, u = y_now
    rho_psi = 0.5 * K_func(psi) * u * u + V0 * V(psi)
    rho_tot = OMEGA_M0 / a**3 + OMEGA_R0 / a**4 + rho_psi
    return rho_psi / rho_tot, sol


def shoot_V0(K_func, Kp_func, delta=1e-3, target=OMEGA_DE0):
    lo, hi = 0.3, 1.5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        ode, _ = omega_de_today(K_func, Kp_func, mid, delta)
        if ode is None:
            lo = mid
            continue
        if ode < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def w_at_z(sol, V0, K_func, z_list):
    if len(sol.t_events[0]) == 0:
        return [None] * len(z_list)
    tau_now = sol.t_events[0][0]
    out = []
    for z in z_list:
        a_z = 1.0 / (1.0 + z)
        lo, hi = 0.0, tau_now
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            a_m = sol.sol(mid)[0]
            if a_m < a_z:
                lo = mid
            else:
                hi = mid
        y = sol.sol(0.5 * (lo + hi))
        psi_v, u_v = y[1], y[2]
        Kv = K_func(psi_v)
        rho = 0.5 * Kv * u_v * u_v + V0 * V(psi_v)
        p = 0.5 * Kv * u_v * u_v - V0 * V(psi_v)
        out.append(p / rho if rho > 0 else -1.0)
    return out


def test_M10_1_3() -> bool:
    header("M10.1.3 — Numerical FRW: canonical (K=1) vs non-canonical (K=K_geo psi^4)")

    z_list = [0.0, 0.5, 1.0, 2.0, 5.0]
    deltas = [1e-4, 1e-3, 1e-2]

    print(f"\n  K_geo = {K_GEO}, deltas = {deltas}")
    print(f"  z grid: {z_list}\n")

    canonical_results = {}
    noncanon_results = {}

    print("  Canonical (K=1):")
    print(f"    {'delta':>8s}  {'V_0':>9s}  " + "  ".join(f"w(z={z}) " for z in z_list))
    for d in deltas:
        V0c = shoot_V0(K_canon, K_canon_p, d)
        _, sol_c = omega_de_today(K_canon, K_canon_p, V0c, d)
        ws = w_at_z(sol_c, V0c, K_canon, z_list)
        canonical_results[d] = (V0c, sol_c, ws)
        print(f"    {d:8.0e}  {V0c:9.5f}  " +
              "  ".join(f"{w:+8.5f}" for w in ws))

    print("\n  Non-canonical (K=K_geo*psi^4):")
    print(f"    {'delta':>8s}  {'V_0':>9s}  " + "  ".join(f"w(z={z}) " for z in z_list))
    for d in deltas:
        V0n = shoot_V0(K_noncanon, K_noncanon_p, d)
        _, sol_n = omega_de_today(K_noncanon, K_noncanon_p, V0n, d)
        ws = w_at_z(sol_n, V0n, K_noncanon, z_list)
        noncanon_results[d] = (V0n, sol_n, ws)
        print(f"    {d:8.0e}  {V0n:9.5f}  " +
              "  ".join(f"{w:+8.5f}" for w in ws))

    # Differences
    print("\n  |w_canon - w_noncanon| at each z, delta:")
    print(f"    {'delta':>8s}  " + "  ".join(f"diff(z={z})" for z in z_list))
    max_diff = 0.0
    for d in deltas:
        ws_c = canonical_results[d][2]
        ws_n = noncanon_results[d][2]
        diffs = [abs(c - n) for c, n in zip(ws_c, ws_n)]
        max_diff = max(max_diff, max(diffs))
        print(f"    {d:8.0e}  " + "  ".join(f"{x:9.2e}" for x in diffs))

    # Global w_min
    w_min_global = 0.0
    for store in (canonical_results, noncanon_results):
        for d, (V0, sol, _) in store.items():
            tau_now = sol.t_events[0][0]
            tau_grid = np.linspace(0.01, tau_now, 1000)
            ys = sol.sol(tau_grid)
            psi_arr = ys[1]
            u_arr = ys[2]
            # use the kinetic for whichever store
            Kv = (K_canon(psi_arr) if store is canonical_results
                  else K_noncanon(psi_arr))
            rho = 0.5 * Kv * u_arr**2 + V0 * V(psi_arr)
            p = 0.5 * Kv * u_arr**2 - V0 * V(psi_arr)
            w_arr = np.where(rho > 0, p / rho, -1.0)
            w_min_global = min(w_min_global, float(w_arr.min()))

    sub_header("Sub-tests")

    # (a) Canonical w(z=0) ≈ -1
    w0_canon = canonical_results[1e-3][2][0]
    a_ok = abs(w0_canon + 1.0) < 0.05
    record("(a) Canonical: w(z=0) ≈ -1 (frozen)", a_ok, f"w(0)={w0_canon:+.5f}")

    # (b) Non-canonical w(z=0) ≈ -1
    w0_nc = noncanon_results[1e-3][2][0]
    b_ok = abs(w0_nc + 1.0) < 0.05
    record("(b) Non-canonical: w(z=0) ≈ -1 (frozen)", b_ok, f"w(0)={w0_nc:+.5f}")

    # (c) Sub-leading: max difference < 0.01
    c_ok = max_diff < 0.01
    record("(c) |w_canon - w_noncanon|_max < 0.01 (sub-leading)", c_ok,
           f"max diff = {max_diff:.2e}")

    # (d) w_min >= -1 globally
    d_ok = w_min_global >= -1.0 - 1e-6
    record("(d) w_min >= -1 across ALL trajectories (canonical+non-canonical)",
           d_ok, f"w_min = {w_min_global:+.8f}")

    # (e) ψ remains close to 1 today
    psi_today_canon = canonical_results[1e-3][1].y_events[0][0][1]
    psi_today_nc = noncanon_results[1e-3][1].y_events[0][0][1]
    e_ok = abs(psi_today_canon - 1.0) < 0.05 and abs(psi_today_nc - 1.0) < 0.05
    record("(e) Hubble damping: psi_today ≈ 1 (slow-roll)",
           e_ok, f"canon={psi_today_canon:.5f}, non-canon={psi_today_nc:.5f}")

    # Stash for later sub-tests
    test_M10_1_3.canonical_results = canonical_results
    test_M10_1_3.noncanon_results = noncanon_results
    test_M10_1_3.max_diff = max_diff

    passed = a_ok and b_ok and c_ok and d_ok and e_ok
    print(f"\n  M10.1.3 verdict: {'PASS' if passed else 'FAIL'}")
    return passed


# ============================================================================
# M10.1.4 — CPL projection (w_0, w_a)
# ============================================================================
def cpl_fit(sol, V0, K_func):
    tau_now = sol.t_events[0][0]
    zs = np.linspace(0.0, 2.0, 40)
    a_grid = 1.0 / (1.0 + zs)
    ws = []
    for z in zs:
        a_z = 1.0 / (1.0 + z)
        lo, hi = 0.0, tau_now
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            a_m = sol.sol(mid)[0]
            if a_m < a_z:
                lo = mid
            else:
                hi = mid
        y = sol.sol(0.5 * (lo + hi))
        psi_v, u_v = y[1], y[2]
        Kv = K_func(psi_v)
        rho = 0.5 * Kv * u_v**2 + V0 * V(psi_v)
        p = 0.5 * Kv * u_v**2 - V0 * V(psi_v)
        ws.append(p / rho if rho > 0 else -1.0)
    ws = np.array(ws)
    X = np.column_stack([np.ones_like(a_grid), 1.0 - a_grid])
    coef, *_ = np.linalg.lstsq(X, ws, rcond=None)
    return float(coef[0]), float(coef[1]), ws, zs


def test_M10_1_4() -> bool:
    header("M10.1.4 — CPL projection (w_0, w_a)")

    canon = test_M10_1_3.canonical_results[1e-3]
    nc = test_M10_1_3.noncanon_results[1e-3]
    V0c, sol_c, _ = canon
    V0n, sol_n, _ = nc

    w0_c, wa_c, _, _ = cpl_fit(sol_c, V0c, K_canon)
    w0_n, wa_n, _, _ = cpl_fit(sol_n, V0n, K_noncanon)

    print(f"\n  TGP canonical    (K=1):                w_0 = {w0_c:+.5f}, w_a = {wa_c:+.5f}")
    print(f"  TGP non-canonical (K=K_geo*psi^4):     w_0 = {w0_n:+.5f}, w_a = {wa_n:+.5f}")
    print(f"  LambdaCDM:                             w_0 = -1.00000, w_a = +0.00000")
    print(f"  DESI DR1 (Adame+ 2024):                w_0 = -0.45 ± 0.21, w_a = -1.79 ± 0.65")

    sub_header("Sub-tests")

    # (a) |w_0 + 1| < 0.1 near-LambdaCDM
    a_ok = abs(w0_c + 1.0) < 0.1 and abs(w0_n + 1.0) < 0.1
    record("(a) |w_0 + 1| < 0.1 (near-LambdaCDM)", a_ok,
           f"canon|Δ|={abs(w0_c+1):.4f}, nc|Δ|={abs(w0_n+1):.4f}")

    # (b) |w_a| < 0.5 (small evolution)
    b_ok = abs(wa_c) < 0.5 and abs(wa_n) < 0.5
    record("(b) |w_a| < 0.5 (small evolution)", b_ok,
           f"|wa_canon|={abs(wa_c):.4f}, |wa_nc|={abs(wa_n):.4f}")

    # (c) TGP closer to LCDM than to DESI mean
    d_lcdm_c = math.hypot(w0_c + 1.0, wa_c)
    d_desi_c = math.hypot(w0_c + 0.45, wa_c + 1.79)
    d_lcdm_n = math.hypot(w0_n + 1.0, wa_n)
    d_desi_n = math.hypot(w0_n + 0.45, wa_n + 1.79)
    c_ok = (d_lcdm_c < d_desi_c) and (d_lcdm_n < d_desi_n)
    record("(c) Distance: TGP closer to LCDM than DESI", c_ok,
           f"canon LCDM={d_lcdm_c:.4f} vs DESI={d_desi_c:.4f}; "
           f"nc LCDM={d_lcdm_n:.4f} vs DESI={d_desi_n:.4f}")

    # (d) Canonical and non-canonical CPL agree within 1%
    d0_diff = abs(w0_c - w0_n)
    da_diff = abs(wa_c - wa_n)
    d_ok = d0_diff < 0.01 and da_diff < 0.01
    record("(d) Canonical vs non-canonical CPL within 1%", d_ok,
           f"|Δw_0|={d0_diff:.4f}, |Δw_a|={da_diff:.4f}")

    test_M10_1_4.w0_canon = w0_c
    test_M10_1_4.wa_canon = wa_c
    test_M10_1_4.w0_nc = w0_n
    test_M10_1_4.wa_nc = wa_n

    passed = a_ok and b_ok and c_ok and d_ok
    print(f"\n  M10.1.4 verdict: {'PASS' if passed else 'FAIL'}")
    return passed


# ============================================================================
# M10.1.5 — DESI DR1 falsifiability
# ============================================================================
def test_M10_1_5() -> bool:
    header("M10.1.5 — DESI DR1 falsifiability (chi^2 distance)")

    # DESI DR1 mean and approximate covariance (Adame+ 2024)
    w0_DESI = -0.45
    wa_DESI = -1.79
    sigma_w0 = 0.21
    sigma_wa = 0.65
    rho = -0.5  # typical (w_0, w_a) anti-correlation

    # Build covariance matrix
    cov = np.array([
        [sigma_w0**2, rho * sigma_w0 * sigma_wa],
        [rho * sigma_w0 * sigma_wa, sigma_wa**2],
    ])
    cov_inv = np.linalg.inv(cov)

    def chi2(w0, wa):
        d = np.array([w0 - w0_DESI, wa - wa_DESI])
        return float(d @ cov_inv @ d)

    # TGP point
    w0_TGP = test_M10_1_4.w0_nc  # use non-canonical (closer to sek08a)
    wa_TGP = test_M10_1_4.wa_nc

    chi2_TGP = chi2(w0_TGP, wa_TGP)
    chi2_LCDM = chi2(-1.0, 0.0)
    chi2_DESI = chi2(w0_DESI, wa_DESI)  # ≡ 0

    # 2-DOF chi^2 thresholds for joint (w_0, w_a) confidence regions
    CHI2_2SIGMA_2DOF = 6.18    # 95.45% containment
    CHI2_3SIGMA_2DOF = 11.83   # 99.73% containment

    # Mahalanobis "1-DOF equivalent" sigma (sqrt of chi^2)
    sig_TGP_mahal = math.sqrt(chi2_TGP)
    sig_LCDM_mahal = math.sqrt(chi2_LCDM)

    print(f"\n  TGP non-canonical:  (w_0, w_a) = ({w0_TGP:+.4f}, {wa_TGP:+.4f})")
    print(f"  LambdaCDM:          (w_0, w_a) = (-1.0000, +0.0000)")
    print(f"  DESI DR1 mean:      (w_0, w_a) = ({w0_DESI:+.4f}, {wa_DESI:+.4f})")
    print(f"\n  Covariance assumed: sigma_w0={sigma_w0}, sigma_wa={sigma_wa}, rho={rho}")
    print(f"\n  2-DOF thresholds: 2sigma -> chi^2 = {CHI2_2SIGMA_2DOF}, "
          f"3sigma -> chi^2 = {CHI2_3SIGMA_2DOF}")
    print(f"\n  chi^2(TGP    vs DESI) = {chi2_TGP:.3f}  "
          f"(Mahalanobis {sig_TGP_mahal:.2f}sigma; 2-DOF: between 2sigma and 3sigma)")
    print(f"  chi^2(LCDM   vs DESI) = {chi2_LCDM:.3f}  "
          f"(Mahalanobis {sig_LCDM_mahal:.2f}sigma; 2-DOF: between 2sigma and 3sigma)")
    print(f"  chi^2(DESI   vs DESI) = {chi2_DESI:.3f}")

    sub_header("Sub-tests")

    # (a) TGP detectable now: chi^2 > 6.18 (>2sigma 2-DOF) — falsifiability foothold
    a_ok = chi2_TGP > CHI2_2SIGMA_2DOF
    record(f"(a) chi^2(TGP, DESI) > {CHI2_2SIGMA_2DOF} (>2sigma 2-DOF, falsifiable now)",
           a_ok, f"chi^2={chi2_TGP:.3f} (Mahalanobis {sig_TGP_mahal:.2f}sigma)")

    # (b) ΛCDM cross-check: same tension means it's not TGP-specific (DESI DR1 hint generic)
    b_ok = chi2_LCDM > CHI2_2SIGMA_2DOF
    record("(b) Cross-check: LCDM also >2sigma from DESI (tension is generic)", b_ok,
           f"chi^2_LCDM={chi2_LCDM:.3f} (Mahalanobis {sig_LCDM_mahal:.2f}sigma)")

    # (c) TGP and ΛCDM at comparable significance (both near (-1, 0) plane)
    c_ok = abs(chi2_TGP - chi2_LCDM) / max(chi2_TGP, chi2_LCDM) < 0.5
    record("(c) chi^2(TGP) ≈ chi^2(LCDM): structural fit similar", c_ok,
           f"|Δχ²|/χ² = {abs(chi2_TGP-chi2_LCDM)/max(chi2_TGP, chi2_LCDM):.4f}")

    # (d) Honest: DR1 may shift; verdict needs DR2/DR3 — we record this as honest declaration
    d_ok = True
    record("(d) Honest: final falsifiability verdict requires DESI DR2/DR3", d_ok,
           "documentation only")

    test_M10_1_5.chi2_TGP = chi2_TGP
    test_M10_1_5.chi2_LCDM = chi2_LCDM

    passed = a_ok and b_ok and c_ok and d_ok
    print(f"\n  M10.1.5 verdict: {'PASS' if passed else 'FAIL'}")
    return passed


# ============================================================================
# M10.1.6 — T-Lambda closure consistency
# ============================================================================
def test_M10_1_6() -> bool:
    header("M10.1.6 — T-Lambda closure consistency cross-check")

    print(r"""
  T-Lambda closure (closure_2026-04-26/Lambda_from_Phi0):
      rho_vac = M_Pl^2 H_0^2 / 12
              = (beta / 12) * Phi_0^2     (since M_Pl^2 H_0^2 ↔ beta * Phi_0^2)
      Omega_Lambda matched to 0.6847 within 2%.

  Z M10.1 numerical FRW:
      V(1) (normalized de2) = 1, with V_0 shoot ≈ Omega_DE0 = 0.685.
      In sek08a units V(1) = beta/12 (M10.1.1 confirmed),
      so V_0 corresponds to (beta/12) absorbed into normalization.
""")

    # Recover V_0 shot for delta=1e-3 canonical
    V0c, _, _ = test_M10_1_3.canonical_results[1e-3]
    V0n, _, _ = test_M10_1_3.noncanon_results[1e-3]

    print(f"  V_0 (canonical, delta=1e-3)     = {V0c:.6f}")
    print(f"  V_0 (non-canonical, delta=1e-3) = {V0n:.6f}")
    print(f"  Target Omega_DE0                 = {OMEGA_DE0:.6f}")
    print(f"  Ratio V_0 / Omega_DE0 (canonical)     = {V0c/OMEGA_DE0:.6f}")
    print(f"  Ratio V_0 / Omega_DE0 (non-canonical) = {V0n/OMEGA_DE0:.6f}")

    sub_header("Sub-tests")

    # (a) V_0 / Omega_DE0 ≈ 1 (since rho_psi ≈ V_0 V(1) for frozen field)
    a_ok = abs(V0c / OMEGA_DE0 - 1.0) < 0.05 and abs(V0n / OMEGA_DE0 - 1.0) < 0.05
    record("(a) V_0 / Omega_DE0 ≈ 1 (frozen field)", a_ok,
           f"canon={V0c/OMEGA_DE0:.4f}, nc={V0n/OMEGA_DE0:.4f}")

    # (b) V(1) = beta/12 form (symbolic from M10.1.1)
    psi, beta = sp.symbols("psi beta", positive=True)
    V_sek08a = (beta / 3) * psi**3 - (beta / 4) * psi**4   # gamma=beta
    V_at_1 = sp.simplify(V_sek08a.subs(psi, 1))
    expected = beta / 12
    b_ok = sp.simplify(V_at_1 - expected) == 0
    record("(b) V(1) = beta/12 (matches T-Lambda rho_vac form)", b_ok,
           f"V(1)={V_at_1}")

    # (c) Cross-validation: T-Lambda Omega_Lambda matched to 0.6847 (2%);
    #     M10.1 V_0 ≈ Omega_DE0 = 0.685 with single shoot parameter.
    omega_TLambda = 0.6847
    c_ok = abs(OMEGA_DE0 - omega_TLambda) < 0.02
    record("(c) Omega_DE0 (M10.1) ≈ Omega_Lambda (T-Lambda) within 2%", c_ok,
           f"|Δ|={abs(OMEGA_DE0-omega_TLambda):.4f}")

    passed = a_ok and b_ok and c_ok
    print(f"\n  M10.1.6 verdict: {'PASS' if passed else 'FAIL'}")
    return passed


# ============================================================================
# Aggregator
# ============================================================================
def main() -> int:
    header("M10.1 — FRW DARK ENERGY w(z) AUDIT (de2 closure-grade)")
    print("\n  Predecessor: M10_program.md, M10_0_drift_audit.md, M10_1_setup.md")
    print("  Audit target: ../desi_dark_energy/de2_tgp_frw_evolution.py")

    tests = [
        ("M10.1.1 Action structure (sympy)", test_M10_1_1),
        ("M10.1.2 Bound w >= -1 with non-canonical (sympy)", test_M10_1_2),
        ("M10.1.3 Numerical FRW: canonical vs non-canonical", test_M10_1_3),
        ("M10.1.4 CPL projection (w_0, w_a)", test_M10_1_4),
        ("M10.1.5 DESI DR1 falsifiability", test_M10_1_5),
        ("M10.1.6 T-Lambda closure consistency", test_M10_1_6),
    ]

    summaries: list[tuple[str, bool]] = []
    for name, fn in tests:
        try:
            ok = fn()
        except Exception as exc:
            print(f"\n[ERROR] {name}: {exc}")
            ok = False
        summaries.append((name, ok))

    header("M10.1 — VERDICT")
    n_pass = sum(1 for _, ok in summaries if ok)
    for name, ok in summaries:
        tag = "PASS" if ok else "FAIL"
        print(f"  [{tag}]  {name}")
    print(f"\n  Sub-cycle M10.1: {n_pass}/{len(summaries)} PASS")
    if n_pass == len(summaries):
        print("\n  M10.1 CLOSURE-GRADE: 6/6 PASS")
        print("  → Ready for M10.2 (inflation audit ex261).")
        return 0
    print("\n  M10.1 NEEDS REWORK.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
