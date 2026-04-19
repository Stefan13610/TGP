#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
de2: TGP FRW EVOLUTION FROM sek05 POTENTIAL — DERIVE w(z) FIRST PRINCIPLES

Extension of de1, which used a B_psi(a) ansatz. Here we derive w(z) by
direct numerical integration of the scalar field equation in FRW from the
sek05 self-coupling potential:

    U(psi) = (beta/3) * psi^3 - (gamma/4) * psi^4,   with beta = gamma

  -> U(1) = gamma/12       (residual vacuum energy = Lambda_eff source)
  -> U'(1) = 0             (vacuum condition)
  -> U''(1) = -gamma < 0   (slow-roll maximum)

Field equation in FRW (flat):
    psi_ddot + 3 H psi_dot + U'(psi) = 0                           (1)
    H^2 = (8 pi G / 3) (rho_m + rho_r + rho_psi)                    (2)
    rho_psi = (1/2) psi_dot^2 + U(psi)                              (3)
    p_psi   = (1/2) psi_dot^2 - U(psi)                              (4)
    w_psi(a) = p_psi / rho_psi                                     (5)

Key structural prediction (derived, not ansatz):
    p_psi - (-rho_psi) = psi_dot^2 >= 0
    =>  w_psi + 1 = psi_dot^2 / rho_psi >= 0
    =>  w_psi >= -1  ALWAYS

No ghost fields, no phantom crossing. This is a strict bound from
the canonical kinetic term. To produce w < -1 would require non-canonical
kinetic, which sek02/sek05 does NOT have.

Parts:
    A. Dimensionless formulation
    B. Initial conditions: slow-roll from unstable max at psi = 1
    C. Numerical integration of coupled ODEs from z = 100 to z = 0
    D. Compute w(z), compare to DESI DR1 CPL fit
    E. Verify w >= -1 structural bound numerically
    F. CPL projection: fit w_0 + w_a (1-a) to TGP w(z)
    G. Falsification summary (baseline for de3, de4)

Units: dimensionless with H0 = 1, Omega_m0 + Omega_DE0 = 1.
"""
from __future__ import annotations

import sys
import numpy as np
from scipy.integrate import solve_ivp, quad

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

hdr = "=" * 78
sub = "-" * 78


def header(title: str) -> None:
    print()
    print(hdr)
    print(f"  {title}")
    print(hdr)


def sub_header(title: str) -> None:
    print()
    print(sub)
    print(f"  {title}")
    print(sub)


# ==========================================================================
# Cosmological parameters (Planck 2018)
# ==========================================================================
OMEGA_M0  = 0.315     # matter today
OMEGA_R0  = 9.1e-5    # radiation today (negligible for DE)
OMEGA_DE0 = 1.0 - OMEGA_M0 - OMEGA_R0   # dark energy today

# ==========================================================================
# Part A. Dimensionless formulation
# ==========================================================================
header("Part A. Dimensionless sek05 scalar field in FRW")

print(r"""
Work in units where:
    8 pi G / 3 = 1,   H_0 = 1,   rho_crit_0 = 1.

Then the Friedmann equation becomes:
    H^2(a) = Omega_r0/a^4 + Omega_m0/a^3 + rho_psi / rho_crit_0

Normalize the scalar field such that at today (psi = psi_0 close to 1):
    rho_psi(today) = Omega_DE0

With sek05 potential U(psi) = (gamma/12) [4 psi^3 - 3 psi^4] (beta = gamma),
we can absorb gamma into the normalization by choosing:

    V(psi) = U(psi) / U(1) = (4 psi^3 - 3 psi^4)

Then V(1) = 1, and normalize the TGP field amplitude so that the
dimensionless scalar density:

    rho_psi_hat = (1/2) (dpsi/dt)^2 + V_0 * V(psi),
        with V_0 = Omega_DE0 if psi is EXACTLY frozen at psi = 1.

For slow-roll psi != 1, V(psi) differs slightly from 1 and V_0 is
re-calibrated by requiring rho_psi_hat(today) = Omega_DE0 (one shooting
parameter; provides initial conditions that match observed Omega_DE0).

Key equations (tau = H_0 t):
    dpsi/d(tau) = u                                                   (A1)
    du/d(tau)  = -3 H_hat u - V_0 V'(psi)                             (A2)
    H_hat(a, psi, u) = sqrt(Omega_r0/a^4 + Omega_m0/a^3 + rho_psi_hat) (A3)
    da/d(tau) = a H_hat                                               (A4)
""")


# Potential and derivatives (normalized so V(1) = 1)
def V(psi):
    return 4*psi**3 - 3*psi**4


def Vp(psi):
    return 12*psi**2 - 12*psi**3  # = 12 psi^2 (1 - psi)


def Vpp(psi):
    return 24*psi - 36*psi**2


print("Potential sanity check:")
for psi_val in [0.9, 0.95, 1.0, 1.05, 1.1]:
    print(f"  psi = {psi_val:5.2f}:  V = {V(psi_val):+.4f}, "
          f"V' = {Vp(psi_val):+.4f}, V'' = {Vpp(psi_val):+.4f}")


# ==========================================================================
# Part B. Initial conditions
# ==========================================================================
header("Part B. Initial conditions from slow-roll on the max at psi=1")

print(r"""
The sek05 potential has psi = 1 as an UNSTABLE MAXIMUM (V''(1) = -12 in
this normalization after rescaling).  Slow-roll dark-energy dynamics
requires the field to roll AWAY from the max on one side; V is symmetric
but we choose the psi < 1 side (arbitrary; same physics for psi > 1).

Slow-roll conditions:
    u = dpsi/dtau ~ -V'(psi) / (3 H)
    rho_psi ~ V_0 * V(psi)  (kinetic term negligible)

At early times (a << 1), H is radiation/matter-dominated and LARGE, so
u is extremely suppressed (frozen field). As a -> 1, H decreases, u grows.
This mimics quintessence slow-roll.

We set INITIAL conditions at a_i = 0.01 (z = 99):
    psi_i = 1 - delta,        with delta ~ 10^-3 (small displacement)
    u_i   = 0                  (frozen; Hubble drag locks it)

Then integrate forward to a = 1 and shoot on V_0 so that Omega_DE0 = 0.685.
""")


def friedmann_H_hat(a, rho_psi_hat):
    """Dimensionless Hubble parameter."""
    rho_m = OMEGA_M0 / a**3
    rho_r = OMEGA_R0 / a**4
    total = rho_m + rho_r + rho_psi_hat
    return np.sqrt(max(total, 0.0))


def rhs(tau, y, V0):
    """ODE system. y = [a, psi, u]."""
    a, psi, u = y
    V_psi = V(psi)
    Vp_psi = Vp(psi)
    rho_psi_hat = 0.5 * u**2 + V0 * V_psi
    H_hat = friedmann_H_hat(max(a, 1e-10), rho_psi_hat)
    da = a * H_hat
    dpsi = u
    du = -3 * H_hat * u - V0 * Vp_psi
    return [da, dpsi, du]


# ==========================================================================
# Part C. Numerical integration + shooting on V_0
# ==========================================================================
header("Part C. Coupled ODE integration from z=99 to z=0")

a_i = 0.01
z_i = 1.0 / a_i - 1.0
tau_max = 200.0    # units of 1/H0; integrates past a=1

def integrate_trajectory(V0, delta=1e-3, psi_branch="+"):
    sign = +1 if psi_branch == "+" else -1
    psi0 = 1.0 + sign * delta
    y0 = [a_i, psi0, 0.0]

    def a_eq_1(tau, y, V0_):
        return y[0] - 1.0
    a_eq_1.terminal = True
    a_eq_1.direction = +1

    sol = solve_ivp(
        rhs,
        (0.0, tau_max),
        y0,
        args=(V0,),
        rtol=1e-9,
        atol=1e-11,
        events=a_eq_1,
        dense_output=True,
        max_step=0.05,
    )
    return sol


def omega_de_at_today(V0, delta=1e-3, psi_branch="+"):
    sol = integrate_trajectory(V0, delta, psi_branch)
    if len(sol.t_events[0]) == 0:
        return None, None
    tau_now = sol.t_events[0][0]
    y_now = sol.y_events[0][0]
    a_now, psi_now, u_now = y_now
    rho_psi = 0.5 * u_now**2 + V0 * V(psi_now)
    rho_tot = OMEGA_M0 / a_now**3 + OMEGA_R0 / a_now**4 + rho_psi
    return rho_psi / rho_tot, sol


def shoot_V0(target_omega_de=OMEGA_DE0, psi_branch="+", delta=1e-3):
    """Bisection on V0 so that Omega_DE(today) = target."""
    lo, hi = 0.3, 1.5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        ode, _ = omega_de_at_today(mid, delta, psi_branch)
        if ode is None:
            lo = mid
            continue
        if ode < target_omega_de:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


deltas_to_try = [1e-4, 1e-3, 1e-2]
results = []
for delta in deltas_to_try:
    V0_shot = shoot_V0(OMEGA_DE0, "+", delta)
    omega_check, sol = omega_de_at_today(V0_shot, delta, "+")
    print(f"  delta = {delta:.1e}:  V_0 = {V0_shot:.6f},  "
          f"Omega_DE(today) = {omega_check:.6f} (target {OMEGA_DE0:.6f})")
    results.append((delta, V0_shot, sol))


# ==========================================================================
# Part D. Compute w(z), compare to DESI
# ==========================================================================
header("Part D. w(z) from TGP FRW trajectories")

print(f"\n  {'delta':>10s}  "
      f"{'z=0':>10s} {'z=0.5':>10s} {'z=1':>10s} "
      f"{'z=2':>10s} {'z=5':>10s}")

w_at_z = {}
for delta, V0_shot, sol in results:
    if sol is None or len(sol.t_events[0]) == 0:
        continue
    tau_now = sol.t_events[0][0]
    # Dense output back in z
    zs = [0.0, 0.5, 1.0, 2.0, 5.0]
    ws = []
    for z in zs:
        a_z = 1.0 / (1.0 + z)
        # Solve for tau where a = a_z
        # Use dense output and bisection
        lo_tau, hi_tau = 0.0, tau_now
        for _ in range(80):
            mid_tau = 0.5*(lo_tau + hi_tau)
            a_mid = sol.sol(mid_tau)[0]
            if a_mid < a_z:
                lo_tau = mid_tau
            else:
                hi_tau = mid_tau
        y_z = sol.sol(0.5*(lo_tau + hi_tau))
        a_v, psi_v, u_v = y_z
        V_psi = V(psi_v)
        rho_psi = 0.5 * u_v**2 + V0_shot * V_psi
        p_psi = 0.5 * u_v**2 - V0_shot * V_psi
        w_psi = p_psi / rho_psi if rho_psi > 0 else -1.0
        ws.append(w_psi)
    w_at_z[delta] = ws
    print(f"  {delta:10.1e}  "
          + "  ".join(f"{w:>8.5f}" for w in ws))

print(r"""
Observations:
  * All trajectories give w >= -1 at all redshifts.
  * Small delta: effectively frozen field, w ~ -1 (cosmological constant).
  * Larger delta: non-trivial evolution, w > -1 by small amount at low z.
""")


# ==========================================================================
# Part E. Structural bound w >= -1 — numerical verification
# ==========================================================================
header("Part E. Verify w >= -1 structural bound across redshift grid")

w_min_global = 0.0
for delta, V0_shot, sol in results:
    if sol is None:
        continue
    tau_now = sol.t_events[0][0]
    tau_grid = np.linspace(0.01, tau_now, 1000)
    ys = sol.sol(tau_grid)
    a_arr = ys[0]
    psi_arr = ys[1]
    u_arr = ys[2]
    rho_arr = 0.5 * u_arr**2 + V0_shot * V(psi_arr)
    p_arr   = 0.5 * u_arr**2 - V0_shot * V(psi_arr)
    w_arr = np.where(rho_arr > 0, p_arr / rho_arr, -1.0)
    w_min = w_arr.min()
    w_max = w_arr.max()
    print(f"  delta = {delta:.1e}:  w_min = {w_min:.8f},  w_max = {w_max:.8f}")
    if w_min < w_min_global:
        w_min_global = w_min

print(f"\n  Global w_min across all trajectories: {w_min_global:.8f}")
assert w_min_global >= -1.0 - 1e-6, "STRUCTURAL BOUND VIOLATED"
print("  [OK] Structural bound w >= -1 holds numerically across all trajectories.")


# ==========================================================================
# Part F. CPL projection of TGP w(z)
# ==========================================================================
header("Part F. Fit CPL (w0, wa) to TGP w(z) for direct DESI comparison")

# Pick the trajectory with intermediate delta = 1e-3
target_delta = 1e-3
for d, v0, s in results:
    if d == target_delta:
        sol_target = s
        V0_target = v0
        break

tau_now_t = sol_target.t_events[0][0]


def a_of_tau_wrapper(tau):
    return sol_target.sol(tau)[0]


def w_of_z(z):
    a_z = 1.0 / (1.0 + z)
    lo, hi = 0.0, tau_now_t
    for _ in range(80):
        mid = 0.5*(lo + hi)
        a_m = a_of_tau_wrapper(mid)
        if a_m < a_z:
            lo = mid
        else:
            hi = mid
    y_z = sol_target.sol(0.5*(lo + hi))
    psi_v, u_v = y_z[1], y_z[2]
    rho = 0.5*u_v**2 + V0_target * V(psi_v)
    p   = 0.5*u_v**2 - V0_target * V(psi_v)
    return p/rho if rho > 0 else -1.0


zs_grid = np.linspace(0.0, 2.0, 40)
a_grid = 1.0 / (1.0 + zs_grid)
w_TGP = np.array([w_of_z(z) for z in zs_grid])

# Fit CPL: w(a) = w0 + wa (1 - a)
from numpy.linalg import lstsq
X = np.column_stack([np.ones_like(a_grid), 1.0 - a_grid])
coef, *_ = lstsq(X, w_TGP, rcond=None)
w0_fit, wa_fit = coef[0], coef[1]

print(f"\n  TGP (delta=1e-3) fit to CPL:")
print(f"    w_0  = {w0_fit:+.5f}")
print(f"    w_a  = {wa_fit:+.5f}")

print(f"\n  DESI DR1 (2024) CPL fit: w_0 = -0.45 +/- 0.21, w_a = -1.79 +/- 0.65")
print(f"  LambdaCDM:              w_0 = -1.00,             w_a =  0.00")

print(f"""
Comparison:
    TGP w_0   = {w0_fit:+.4f}  (close to -1, slightly above)
    TGP w_a   = {wa_fit:+.4f}  (very small, near zero)
    DESI w_0  = -0.45 +/- 0.21  (far from -1)
    DESI w_a  = -1.79 +/- 0.65  (strongly negative)

  Distance in (w_0, w_a) plane:
    TGP to LambdaCDM:  {np.sqrt(w0_fit**2 + 2*w0_fit + 1 + wa_fit**2):.4f}
    TGP to DESI mean:  {np.sqrt((w0_fit+0.45)**2 + (wa_fit+1.79)**2):.4f}

  TGP is MUCH closer to LambdaCDM than to DESI DR1 CPL best-fit.
  If DESI DR1 values hold up to full-sky DR2/DR3 at high significance,
  TGP is FALSIFIED.  Quantitative falsification metric in de4.
""")


# ==========================================================================
# Part G. Summary & roadmap
# ==========================================================================
header("Part G. Summary & roadmap to de3/de4")

print(rf"""
TGP DARK-ENERGY PREDICTIONS (first-principles, from sek05):

  1. w(z) >= -1 for ALL z:
     PROVEN analytically:  w + 1 = (psi_dot)^2 / rho_psi >= 0
     VERIFIED numerically:  global w_min = {w_min_global:.8f}

  2. Near-LambdaCDM behavior today:
     w_0 TGP = {w0_fit:+.4f}   (vs LambdaCDM -1, vs DESI DR1 -0.45)
     w_a TGP = {wa_fit:+.4f}   (vs LambdaCDM 0, vs DESI DR1 -1.79)

  3. No phantom crossing:
     STRUCTURAL: follows from canonical kinetic term in sek05 action.
     Phantom (w < -1) would require wrong-sign kinetic (ghost).
     sek02+sek05 has NO ghost degrees of freedom (verified gs08b).

  4. Matches observed Omega_DE0 = 0.685 with single shoot parameter V_0:
     V_0 = {results[0][1]:.4f} (for delta=1e-4, effectively cosmological constant)

PATH OF DE EVOLUTION:
  * Early: Hubble drag freezes psi at psi = 1 + O(delta).
  * Late: H(a) decreases, psi slowly rolls.
  * Today: small kinetic contribution; w still very close to -1.
  * Future: unclear (depends on sign of delta and potential shape).

FALSIFICATION POSITION:

  TGP IS FALSIFIED IF:
    (a) DESI DR2/DR3 confirms phantom crossing w(z) < -1 at any z
        with > 3 sigma significance.
    (b) DESI DR3 combined with CMB/SN constrains (w_0, w_a) to a
        region more than 3 sigma from TGP prediction (w_0 ~ -1, w_a ~ 0).

  TGP IS SUPPORTED IF:
    (a) w(z) asymptotes to -1 (consistent with cosmological constant).
    (b) CPL fit tightens around (w_0 = -1, w_a = 0).

NEXT STEPS:
  de3: Fit TGP template, CPL, and LambdaCDM to DESI DR1 BAO + CMB + SN.
       Chi^2 comparison, model selection.
  de4: Forecast DESI DR2/DR3 precision.
       What phantom-crossing sigma kills TGP?
       What measurement converges with TGP?
""")

header("de2 complete. Derived TGP w(z) from sek05 FRW. w >= -1 verified.")
