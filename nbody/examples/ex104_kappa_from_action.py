#!/usr/bin/env python3
"""
ex104_kappa_from_action.py
==========================
Verification: corrected matter-coupling constant kappa from unified TGP action.

Solves the cosmological evolution equation for psi(t) = Phi(t)/Phi_0:

    psi'' + (3 - eps_H) psi' + 2 (psi')^2 / psi
        = [ -V_self(psi) - kappa * Omega_m / a^3 ] / H^2

where primes are d/dN  (N = ln a),  and

    V_self(psi) = gamma * psi^2 * (1 - psi)
    H^2         = (Omega_r/a^4 + Omega_m/a^3 + Omega_Lambda) / psi
    eps_H       = -H'/H   (slow-roll parameter)

Two coupling prescriptions:
    kappa_old = 3 / (2 Phi_0)   ~  0.0608   (current theory)
    kappa_new = 3 / (4 Phi_0)   ~  0.0304   (corrected from unified action)

Observable:  |Gdot/G| / H0 = |psi'(0)| / psi(0)
LLR bound:   |Gdot/G| / H0 < 0.02

Author : TGP collaboration (auto-generated verification)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ──────────────────────────── physical parameters ────────────────────────────
PHI0    = 24.66
GAMMA   = 1.0          # normalised V_self coefficient
OMEGA_R = 9.1e-5
OMEGA_M = 0.315
OMEGA_L = 0.685
PSI_INI = 7.0 / 6.0    # initial psi  (attractor value)

Z_START = 1.0e9         # redshift at BBN
Z_END   = 0.0

N_START = np.log(1.0 / (1.0 + Z_START))   # ~ -20.72
N_END   = 0.0

LLR_BOUND = 0.02       # |Gdot/G| / H0 upper limit


# ──────────────────────────── RHS of the ODE system ──────────────────────────
def make_rhs(kappa):
    """Return the RHS  f(N, y)  for  y = [psi, psi']."""

    def rhs(N, y):
        psi, dpsi = y

        # scale factor
        a = np.exp(N)
        a3 = a**3
        a4 = a3 * a

        # Friedmann: H^2 (in H0^2 units) divided by psi
        rho_tot = OMEGA_R / a4 + OMEGA_M / a3 + OMEGA_L
        H2 = rho_tot / psi          # H^2 in H0^2 units

        # slow-roll parameter  eps_H = -H'/H
        # H^2 = rho_tot / psi  =>  d(H^2)/dN = drho/dN / psi - rho_tot dpsi / psi^2
        # drho/dN = -4 Omega_r/a^4 - 3 Omega_m/a^3   (Omega_L const)
        drho_dN = -4.0 * OMEGA_R / a4 - 3.0 * OMEGA_M / a3
        dH2_dN  = (drho_dN * psi - rho_tot * dpsi) / psi**2
        # eps_H = -H'/H = - (1/(2H^2)) dH^2/dN
        eps_H = -0.5 * dH2_dN / H2

        # self-interaction potential
        V_self = GAMMA * psi**2 * (1.0 - psi)

        # source term
        source = (-V_self - kappa * OMEGA_M / a3) / H2

        # second-order equation  =>  dpsi' = source - (3 - eps_H) psi' - 2 (psi')^2 / psi
        ddpsi = source - (3.0 - eps_H) * dpsi - 2.0 * dpsi**2 / psi

        return [dpsi, ddpsi]

    return rhs


# ──────────────────────────── integrator wrapper ─────────────────────────────
def evolve_psi(kappa, psi_ini=PSI_INI, dpsi_ini=0.0,
               N_start=N_START, N_end=N_END):
    """Integrate psi(N) from N_start to N_end.  Returns sol object."""
    rhs = make_rhs(kappa)
    y0  = [psi_ini, dpsi_ini]

    sol = solve_ivp(rhs, [N_start, N_end], y0,
                    method='Radau',
                    rtol=1e-10, atol=1e-12,
                    dense_output=True,
                    max_step=0.5)
    if not sol.success:
        raise RuntimeError(f"Integration failed for kappa={kappa:.6f}: {sol.message}")
    return sol


def observables(sol):
    """Extract observables at N=0  (z=0, today)."""
    psi_today, dpsi_today = sol.sol(0.0)
    Gdot_over_G_H0 = abs(dpsi_today) / psi_today
    G_ratio = psi_today / PSI_INI       # G_BBN / G_today  (since G ~ 1/Phi ~ 1/psi)
    delta_psi = abs(psi_today - PSI_INI)
    return {
        'psi_0'          : psi_today,
        'dpsi_0'         : dpsi_today,
        'delta_psi'      : delta_psi,
        'Gdot_G_over_H0' : Gdot_over_G_H0,
        'G_BBN_over_G0'  : 1.0 / G_ratio,   # G_BBN/G_0 = psi_0 / psi_ini  (G~1/psi => ratio inverts)
        'pass_LLR'       : Gdot_over_G_H0 < LLR_BOUND,
    }


# ──────────────────────────── main comparison ────────────────────────────────
def main():
    kappa_old = 3.0 / (2.0 * PHI0)
    kappa_new = 3.0 / (4.0 * PHI0)

    print("=" * 72)
    print("  TGP cosmological kappa verification  (ex104)")
    print("=" * 72)
    print(f"  Phi_0       = {PHI0}")
    print(f"  gamma       = {GAMMA}")
    print(f"  psi_ini     = {PSI_INI:.6f}  (= 7/6)")
    print(f"  Omega_r     = {OMEGA_R}")
    print(f"  Omega_m     = {OMEGA_M}")
    print(f"  Omega_L     = {OMEGA_L}")
    print(f"  z_start     = {Z_START:.0e}  (N_start = {N_START:.4f})")
    print(f"  LLR bound   = {LLR_BOUND}")
    print(f"  kappa_old   = 3/(2 Phi_0) = {kappa_old:.6f}")
    print(f"  kappa_new   = 3/(4 Phi_0) = {kappa_new:.6f}")
    print("=" * 72)

    # ── solve for both kappa values ──────────────────────────────────────
    results = {}
    for label, kappa in [("kappa_old", kappa_old), ("kappa_new", kappa_new)]:
        print(f"\n  Integrating  {label} = {kappa:.6f} ...")
        sol = evolve_psi(kappa)
        obs = observables(sol)
        results[label] = obs
        print(f"    done  ({sol.t.size} steps)")

    # ── summary table ────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("  RESULTS SUMMARY")
    print("=" * 72)
    header = f"  {'Observable':<28s} {'kappa_old':>14s} {'kappa_new':>14s}  {'LLR':>8s}"
    print(header)
    print("  " + "-" * 68)

    r_old = results["kappa_old"]
    r_new = results["kappa_new"]

    rows = [
        ("kappa",              f"{kappa_old:.6f}",               f"{kappa_new:.6f}",            ""),
        ("psi(z=0)",           f"{r_old['psi_0']:.8f}",         f"{r_new['psi_0']:.8f}",       ""),
        ("delta_psi",          f"{r_old['delta_psi']:.2e}",     f"{r_new['delta_psi']:.2e}",   ""),
        ("|Gdot/G|/H0",       f"{r_old['Gdot_G_over_H0']:.6f}", f"{r_new['Gdot_G_over_H0']:.6f}", f"< {LLR_BOUND}"),
        ("G_BBN / G_today",   f"{r_old['G_BBN_over_G0']:.6f}", f"{r_new['G_BBN_over_G0']:.6f}", ""),
        ("pass LLR?",         "YES" if r_old['pass_LLR'] else "NO",
                              "YES" if r_new['pass_LLR'] else "NO",        ""),
    ]
    for name, v1, v2, bound in rows:
        print(f"  {name:<28s} {v1:>14s} {v2:>14s}  {bound:>8s}")

    # ── kappa scan to find critical value ────────────────────────────────
    print("\n" + "=" * 72)
    print("  KAPPA SCAN  (finding kappa_crit where |Gdot/G|/H0 = 0.02)")
    print("=" * 72)

    kappa_scan = np.linspace(0.01, 0.08, 15)
    gdot_vals  = []

    print(f"\n  {'kappa':>10s}  {'|Gdot/G|/H0':>14s}  {'pass LLR?':>10s}")
    print("  " + "-" * 38)

    for k in kappa_scan:
        try:
            sol = evolve_psi(k)
            obs = observables(sol)
            gv  = obs['Gdot_G_over_H0']
            gdot_vals.append(gv)
            tag = "YES" if gv < LLR_BOUND else "NO"
            print(f"  {k:10.5f}  {gv:14.6f}  {tag:>10s}")
        except RuntimeError as e:
            gdot_vals.append(np.nan)
            print(f"  {k:10.5f}  {'FAILED':>14s}  {'---':>10s}")

    # ── bisection for kappa_crit (both crossings) ────────────────────────
    def gdot_minus_bound(k):
        sol = evolve_psi(k)
        obs = observables(sol)
        return obs['Gdot_G_over_H0'] - LLR_BOUND

    # find ALL sign changes (non-monotonic profile!)
    kappa_arr = np.array(kappa_scan)
    gdot_arr  = np.array(gdot_vals)
    valid     = ~np.isnan(gdot_arr)
    kv = kappa_arr[valid]
    gv = gdot_arr[valid]
    diff_sign = gv - LLR_BOUND

    crossings = []
    for i in range(len(diff_sign) - 1):
        if diff_sign[i] * diff_sign[i + 1] < 0:
            k_lo, k_hi = kv[i], kv[i + 1]
            try:
                kc = brentq(gdot_minus_bound, k_lo, k_hi,
                            xtol=1e-8, rtol=1e-10)
                crossings.append(kc)
            except Exception:
                pass

    # find minimum of |Gdot/G|/H0
    i_min = np.argmin(gv)
    kappa_opt = kv[i_min]
    gdot_min  = gv[i_min]

    print("\n  " + "-" * 38)
    print(f"  Minimum |Gdot/G|/H0 = {gdot_min:.6f}  at kappa = {kappa_opt:.5f}")
    if len(crossings) >= 2:
        print(f"  LLR window:  kappa in [{crossings[0]:.6f}, {crossings[1]:.6f}]")
        print(f"  kappa_new  = {kappa_new:.6f}  (3/(4 Phi_0))")
        in_window = crossings[0] <= kappa_new <= crossings[1]
        if in_window:
            print(f"  => kappa_new is INSIDE the LLR window  =>  PASSES")
        else:
            print(f"  => kappa_new is OUTSIDE the LLR window  =>  FAILS")
    elif len(crossings) == 1:
        print(f"  Single crossing at kappa_crit = {crossings[0]:.6f}")
        print(f"  kappa_new = {kappa_new:.6f}: {'PASSES' if r_new['pass_LLR'] else 'FAILS'}")
    else:
        print("  No crossing found in scan range.")

    # ── final verdict ────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("  VERDICT")
    print("=" * 72)
    if r_new['pass_LLR'] and not r_old['pass_LLR']:
        print("  The corrected coupling  kappa = 3/(4 Phi_0)  from the unified action")
        print("  brings |Gdot/G|/H0 within the LLR bound, while the old value")
        print("  kappa = 3/(2 Phi_0)  violates it.  The new derivation is confirmed.")
    elif r_new['pass_LLR'] and r_old['pass_LLR']:
        print("  Both kappa values pass the LLR bound.")
    elif not r_new['pass_LLR'] and not r_old['pass_LLR']:
        print("  Neither kappa value passes the LLR bound.")
        print("  Check parameters or numerical convergence.")
    else:
        print("  Unexpected: old passes but new fails. Check derivation.")
    print("=" * 72)


if __name__ == "__main__":
    main()
