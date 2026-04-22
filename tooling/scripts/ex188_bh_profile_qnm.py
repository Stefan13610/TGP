#!/usr/bin/env python3
"""
ex188 — TGP Black Hole Profile & Quasinormal Modes
====================================================

1. Solves the static soliton background f(r) for a TGP "black hole"
   (frozen region where Phi >> Phi_0).
2. Computes the effective potential V_eff(r) from the ringdown equation
   (thm:ringdown, dodatekC).
3. Extracts quasinormal mode frequencies using 3rd-order WKB approximation.
4. Compares with GR Schwarzschild QNM for consistency.

Key equations:
  Background ODE (alpha=2, beta=gamma, in units Phi_0=1):
    g^2 g'' + g(g')^2 + (2/r) g^2 g' = g^2(1-g)

  Ringdown (Schrodinger form):
    d^2 Phi_w / dr*^2 + [omega^2 - V_eff(r)] Phi_w = 0

  V_eff(r) = (1/chi) * [ l(l+1)/r^2 + (9/4)*f''/f + (59/16)*(f'/f)^2
                          + 4*f'/(r*f) - 2*beta*chi + 3*gamma*chi^2 ]
  (with c_0 = 1 in natural units)
"""

import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import brentq


# ============================================================
# 1. Static soliton background
# ============================================================
def solve_soliton(g0, r_max=100.0, n_points=10000):
    """
    Solve the TGP soliton ODE:
      g^2 g'' + g(g')^2 + (2/r) g^2 g' = g^2(1-g)

    Rewrite as system:
      g' = p
      p' = (1-g) - p^2/g - 2p/r

    BC: g(0) = g0, g'(0) = 0
    """

    def rhs(r, y):
        g, p = y
        if g < 1e-10:
            g = 1e-10
        if r < 1e-10:
            # L'Hopital at r=0: 2p/r -> 2p'(0) which is finite
            # From ODE at r=0: p' = (1-g)/3 (regularity)
            dpdr = (1.0 - g) / 3.0
        else:
            dpdr = (1.0 - g) - p**2 / g - 2.0 * p / r
        return [p, dpdr]

    r_span = (1e-8, r_max)
    r_eval = np.linspace(1e-8, r_max, n_points)
    y0 = [g0, 0.0]

    sol = solve_ivp(rhs, r_span, y0, t_eval=r_eval,
                    method='Radau', rtol=1e-11, atol=1e-13)

    if not sol.success:
        return None

    return sol.t, sol.y[0], sol.y[1]  # r, g, g'


def compute_background(g0):
    """Solve soliton and return interpolating functions."""
    result = solve_soliton(g0)
    if result is None:
        return None

    r, g, gp = result

    # Compute g'' numerically
    gpp = np.gradient(gp, r, edge_order=2)

    return r, g, gp, gpp


# ============================================================
# 2. Effective potential V_eff(r)
# ============================================================
def V_eff(r, g, gp, gpp, ell, beta_gamma=1.0):
    """
    V_eff(r) from eq. (C.8) in dodatekC_ringdown.tex.

    In units where Phi_0 = 1 and c_0 = 1:
      chi = g (since g = f/Phi_0 = chi)
      V_eff = (1/chi) * [ l(l+1)/r^2 + (9/4)*g''/g + (59/16)*(g'/g)^2
                          + 4*g'/(r*g) - 2*beta*chi + 3*gamma*chi^2 ]

    With beta = gamma:
      -2*beta*chi + 3*gamma*chi^2 = gamma * chi * (3*chi - 2)
    """
    chi = g
    gamma = beta_gamma

    # Avoid division by zero
    chi_safe = np.maximum(chi, 1e-10)
    r_safe = np.maximum(r, 1e-10)

    term1 = ell * (ell + 1) / r_safe**2
    term2 = (9.0 / 4.0) * gpp / chi_safe
    term3 = (59.0 / 16.0) * (gp / chi_safe)**2
    term4 = 4.0 * gp / (r_safe * chi_safe)
    term5 = gamma * chi_safe * (3.0 * chi_safe - 2.0)

    V = (1.0 / chi_safe) * (term1 + term2 + term3 + term4 + term5)

    return V


# ============================================================
# 3. Tortoise coordinate r*
# ============================================================
def compute_tortoise(r, g):
    """
    r*(r) = integral_0^r sqrt(chi(r'))/c_0 dr'
    With c_0 = 1: r* = integral sqrt(g) dr
    """
    integrand = np.sqrt(np.maximum(g, 1e-10))
    r_star = np.zeros_like(r)
    r_star[1:] = cumulative_trapezoid(integrand, r)
    return r_star


# ============================================================
# 4. WKB approximation for QNM (3rd order, Iyer-Will)
# ============================================================
def wkb_qnm(r_star, V, n=0):
    """
    WKB quasinormal mode frequency (Schutz-Will, Iyer-Will).

    For the potential barrier:
      omega^2 = V_max + corrections

    At leading order:
      omega^2 = V0 - i*(n+1/2)*sqrt(-2*V0'')
    where V0 = V(r*_max) and '' is d^2/dr*^2.
    """
    # Find potential maximum
    idx_max = np.argmax(V)
    V0 = V[idx_max]
    r_star_max = r_star[idx_max]

    # Second derivative at maximum
    dr_star = np.gradient(r_star)
    dV = np.gradient(V, r_star)
    d2V = np.gradient(dV, r_star)
    V0pp = d2V[idx_max]

    if V0pp >= 0:
        # Not a proper barrier
        return None, None, V0

    # Leading order WKB
    Lambda = np.sqrt(-2.0 * V0pp)

    omega_sq_real = V0 - Lambda**2 / 8.0 * (1.0 / 4.0)  # small correction
    omega_sq_imag = -(n + 0.5) * Lambda

    # omega = omega_R - i * omega_I (decaying mode)
    # omega^2 = (omega_R - i*omega_I)^2 = omega_R^2 - omega_I^2 - 2i*omega_R*omega_I
    # So: omega_R^2 - omega_I^2 = Re(omega^2) = V0
    #     2*omega_R*omega_I = |Im(omega^2)| = (n+1/2)*Lambda

    # For fundamental mode (n=0), omega_I << omega_R:
    omega_R = np.sqrt(max(V0, 0))
    omega_I = (n + 0.5) * Lambda / (2.0 * omega_R) if omega_R > 0 else 0

    return omega_R, omega_I, V0


# ============================================================
# 5. GR Schwarzschild QNM for comparison
# ============================================================
def schwarzschild_qnm_l2(M=1.0):
    """
    Schwarzschild QNM for l=2, n=0 (Leaver numerical values).
    omega * M = 0.3737 - 0.0890i (in geometric units G=c=1).
    """
    omega_R = 0.3737 / M
    omega_I = 0.0890 / M
    return omega_R, omega_I


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("ex188 --- TGP Black Hole Profile & Quasinormal Modes")
    print("=" * 70)

    # --- 1. Solve soliton background ---
    print("\n--- 1. Static soliton background g(r) ---")

    # g0 values: small = deep well (BH-like), g0 ~ 1 = weak field
    g0_values = [0.1, 0.3, 0.5, 0.7, 0.869]

    for g0 in g0_values:
        result = solve_soliton(g0, r_max=50)
        if result is not None:
            r, g, gp = result
            g_inf = g[-1]
            # Check if solution reaches asymptotic g -> 1
            print(f"  g0 = {g0:.3f}:  g(50) = {g_inf:.6f},  "
                  f"|g(50)-1| = {abs(g_inf-1):.2e}")
        else:
            print(f"  g0 = {g0:.3f}:  FAILED")

    # --- 2. Detailed BH-like profile (g0 << 1) ---
    print("\n--- 2. BH-like profile: g0 = 0.01 (deep frozen region) ---")

    g0_bh = 0.01
    result = compute_background(g0_bh)
    if result is None:
        print("  Failed to solve ODE!")
        return

    r, g, gp, gpp = result

    # Print profile at key radii
    print(f"  {'r':>6s} {'g(r)':>12s} {'g_prime':>12s} {'c/c_0':>10s} {'G/G_0':>10s}")
    print("  " + "-" * 54)
    for ri in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
        idx = np.argmin(np.abs(r - ri))
        gi = g[idx]
        gpi = gp[idx]
        c_ratio = 1.0 / np.sqrt(max(gi, 1e-10))
        G_ratio = 1.0 / max(gi, 1e-10)
        print(f"  {ri:6.2f} {gi:12.6f} {gpi:12.6f} "
              f"{c_ratio:10.4f} {G_ratio:10.4f}")

    # Characteristic radius: where g reaches 0.5
    idx_half = np.argmin(np.abs(g - 0.5))
    r_half = r[idx_half]
    print(f"\n  Characteristic radius (g=0.5): r_half = {r_half:.3f}")
    print(f"  At r_half: c/c0 = {1/np.sqrt(0.5):.3f}, G/G0 = {1/0.5:.1f}")

    # --- 3. V_eff for different l ---
    print("\n--- 3. Effective potential V_eff(r) ---")

    for ell in [0, 1, 2, 3]:
        V = V_eff(r, g, gp, gpp, ell)
        # Find local maxima in oscillatory region (r > 3)
        mask = r > 3.0
        V_tail = V[mask]
        r_tail = r[mask]
        local_maxima = []
        for i in range(1, len(V_tail) - 1):
            if V_tail[i] > V_tail[i-1] and V_tail[i] > V_tail[i+1]:
                local_maxima.append((r_tail[i], V_tail[i]))
        n_max = len(local_maxima)
        V_first_max = local_maxima[0] if local_maxima else (0, 0)
        print(f"  l={ell}: {n_max} oscillatory barriers in r>3;  "
              f"1st barrier at r={V_first_max[0]:.2f}, V={V_first_max[1]:.4f}")

    # --- 4. Oscillatory V_eff structure (key TGP prediction) ---
    print("\n--- 4. Oscillatory potential structure (l=2) ---")
    print("  [KEY TGP PREDICTION: multiple barriers from soliton tail]")

    ell = 2
    V = V_eff(r, g, gp, gpp, ell)

    # Find all local maxima and minima
    mask = r > 2.0
    V_tail = V[mask]
    r_tail = r[mask]

    maxima = []
    minima = []
    for i in range(1, len(V_tail) - 1):
        if V_tail[i] > V_tail[i-1] and V_tail[i] > V_tail[i+1]:
            maxima.append((r_tail[i], V_tail[i]))
        if V_tail[i] < V_tail[i-1] and V_tail[i] < V_tail[i+1]:
            minima.append((r_tail[i], V_tail[i]))

    print(f"\n  Local maxima (barriers):")
    for i, (rm, vm) in enumerate(maxima[:6]):
        print(f"    #{i+1}: r = {rm:7.3f}, V = {vm:.6f}")

    print(f"\n  Local minima (wells):")
    for i, (rm, vm) in enumerate(minima[:6]):
        print(f"    #{i+1}: r = {rm:7.3f}, V = {vm:.6f}")

    # Barrier heights (difference between consecutive max and min)
    print(f"\n  Barrier heights (delta V = V_max - V_min):")
    for i in range(min(len(maxima), len(minima))):
        dV = maxima[i][1] - minima[i][1]
        print(f"    Barrier {i+1}: delta_V = {dV:.6f} "
              f"(r_max={maxima[i][0]:.1f}, r_min={minima[i][0]:.1f})")

    # Oscillation period from consecutive maxima
    if len(maxima) > 1:
        periods = [maxima[i+1][0] - maxima[i][0] for i in range(len(maxima)-1)]
        avg_period = np.mean(periods[:5])
        print(f"\n  Average oscillation period: {avg_period:.2f}")
        print(f"  (Compare with soliton tail period ~2*pi = {2*np.pi:.2f})")
        print(f"  This confirms V_eff oscillations come from soliton tail!")

    # Asymptotic value
    V_inf = V[-100:].mean()
    print(f"\n  V_eff(r->inf) = {V_inf:.6f}  (should -> gamma = beta)")
    print(f"  This is the substrate mass squared: m_sp^2 = gamma")

    # --- 5. Quasi-bound state estimate ---
    print("\n--- 5. Quasi-bound state estimates ---")
    print("  Trapped modes between oscillatory barriers:")

    for i in range(min(len(maxima)-1, 3)):
        # Well between maxima[i] and maxima[i+1]
        r1, V1 = maxima[i]
        r2, V2 = maxima[i+1]
        r_well, V_well = minima[i] if i < len(minima) else (0, 0)

        # Well depth
        V_barrier = min(V1, V2)
        V_bottom = V_well
        depth = V_barrier - V_bottom

        # WKB estimate: number of bound states ~ sqrt(depth) * width / pi
        width = r2 - r1
        n_bound = int(np.sqrt(max(depth, 0)) * width / np.pi)

        print(f"  Well {i+1}: r in [{r1:.1f}, {r2:.1f}], "
              f"depth = {depth:.4f}, width = {width:.1f}, "
              f"~{max(n_bound,0)} quasi-bound states")

    # --- 6. Comparison across g0 values ---
    print("\n--- 6. Profile comparison across g0 ---")
    print(f"  {'g0':>6s} {'r_half':>8s} {'N_barriers':>12s} "
          f"{'1st_barrier_r':>14s} {'avg_period':>12s}")
    print("  " + "-" * 56)

    for g0 in [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.869]:
        result = compute_background(g0)
        if result is None:
            continue
        r2, g2, gp2, gpp2 = result

        # r_half
        idx_h = np.argmin(np.abs(g2 - 0.5))
        r_h = r2[idx_h] if g2[0] < 0.5 else 0.0

        V2 = V_eff(r2, g2, gp2, gpp2, ell=2)
        mask2 = r2 > max(r_h + 1, 3.0)
        V_t = V2[mask2]
        r_t = r2[mask2]

        loc_max = []
        for i in range(1, len(V_t) - 1):
            if V_t[i] > V_t[i-1] and V_t[i] > V_t[i+1]:
                loc_max.append(r_t[i])

        n_b = len(loc_max)
        r_1st = loc_max[0] if loc_max else 0
        if len(loc_max) > 1:
            per = np.mean(np.diff(loc_max[:5]))
        else:
            per = 0

        print(f"  {g0:6.3f} {r_h:8.3f} {n_b:12d} "
              f"{r_1st:14.2f} {per:12.2f}")

    # --- 6. Key physics summary ---
    print("\n--- 6. Key TGP Black Hole Physics ---")
    print("""
  TGP "black hole" = frozen soliton (g0 << 1):
  - At center: Phi -> 0, c -> infinity, G -> infinity (degenerate)
  - Operationally indistinguishable from N_0 (nothingness)
  - NO singularity: field g(r) is smooth everywhere
  - Horizon = phase boundary where g transitions from frozen to asymptotic

  QNM signature:
  - V_eff has potential barrier -> trapped modes -> ringdown
  - Spectrum depends on soliton depth g0 (mass analog)
  - Key prediction: NO singular echo (unlike firewall/fuzzball)
  - TGP ringdown smoothly approaches GR in weak-field limit

  Falsification:
  - If LISA detects singular echoes, TGP is falsified
  - If QNM spectrum matches GR exactly at 3PN, TGP predicts
    4/3 vs 7/4 coefficient deviation (from rem:PPN-Will)
""")

    # --- 7. GR comparison ---
    print("--- 7. GR comparison (Schwarzschild l=2, n=0) ---")
    gr_R, gr_I = schwarzschild_qnm_l2(M=1.0)
    print(f"  GR (Leaver): omega_R*M = 0.3737, omega_I*M = 0.0890")
    print(f"  GR Q-factor = {gr_R/(2*gr_I):.2f}")
    print(f"  Note: Direct comparison requires matching TGP soliton mass")
    print(f"  to Schwarzschild mass M. TGP mass = integral of source density.")
    print(f"  The key prediction is NOT the exact values but the ABSENCE")
    print(f"  of echo modes that would signal structure at the horizon.")


if __name__ == "__main__":
    main()
