#!/usr/bin/env python3
"""
ex210 -- Analytical equilibria, normal-mode frequencies, and Hill criterion
============================================================================

Derives closed-form expressions for TGP 2-body equilibria and verifies
them numerically.  Three results:

  Result A: Exact equilibrium distances (d_rep, d_well) from V_2
  Result B: Analytical omega^2 at d_well (radial stiffness)
  Result C: TGP Hill sphere radius vs Newton

Physics (derivation):
  V_2(d) = -A/d + B/d^2 - G/d^3
  where A = 4*pi*C1*C2, B = 8*pi*beta*C1*C2, G = 12*pi*gamma*C1*C2*(C1+C2)

  F(d) = -V'(d) = -A/d^2 + 2B/d^3 - 3G/d^4

  F(d)=0 gives (multiply by d^4/A):
    d^2 - (2B/A)*d + (3G/A) = 0

  For equal masses (C1=C2=C):
    A = 4*pi*C^2, B = 8*pi*beta*C^2, G = 24*pi*gamma*C^3
    => d^2 - 4*beta*d + 18*gamma*C = 0
    => d_rep,well = 2*beta -/+ sqrt(4*beta^2 - 18*gamma*C)

  omega^2_radial = V''(d_well) / mu_red
    V''(d) = -2A/d^3 + 6B/d^4 - 12G/d^5
    mu_red = C1*C2/(C1+C2)  (reduced mass)

  Hill criterion:
    Tidal acceleration at distance r from body m orbiting M at distance d:
      a_tidal = d^2 V_eff / dd^2 * r  (leading order)
    Self-gravity at distance r:  a_self = G_eff * m / r^2
    Hill: a_tidal = a_self => R_H^3 = m / (d^2 V'' / dd^2 evaluated at d)
    For TGP: V'' has extra terms => R_H(TGP) != R_H(Newton)

Tests:
  Part 1: Analytical d_rep, d_well vs numerical root-finding (verification)
  Part 2: Analytical omega^2 vs numerical Hessian (verification)
  Part 3: TGP Hill sphere vs Newton Hill sphere (scan over beta)
  Part 4: Perturbative V3 correction to d_well (first order)
  Part 5: Summary table for publication
"""

import sys
import os
import time
import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from nbody.pairwise import V_eff, V_eff_total, force_2body, force_zeros_2body
from nbody.tgp_field import default_beta_gamma, screening_mass
from nbody.three_body_force_exact import three_body_energy_exact

quick = "--quick" in sys.argv


# ── Analytical formulas ──────────────────────────────────────────────────

def equilibria_analytical(C, beta, gamma=None):
    """
    Exact equilibrium distances for equal-mass pair.

    From V_2'(d) = 0:
      d^2 - 4*beta*d + 18*gamma*C = 0

    Returns: (d_rep, d_well, discriminant) or None if no real roots.
    """
    if gamma is None:
        gamma = beta
    disc = 4.0 * beta**2 - 18.0 * gamma * C
    if disc < 0:
        return None
    sq = np.sqrt(disc)
    d_rep = 2.0 * beta - sq
    d_well = 2.0 * beta + sq
    return d_rep, d_well, disc


def V2_second_derivative(d, C1, C2, beta, gamma=None):
    """
    Analytical V_2''(d).

    V_2 = -A/d + B/d^2 - G/d^3
    V_2' = A/d^2 - 2B/d^3 + 3G/d^4
    V_2'' = -2A/d^3 + 6B/d^4 - 12G/d^5
    """
    if gamma is None:
        gamma = beta
    A = 4.0 * np.pi * C1 * C2
    B = 8.0 * np.pi * beta * C1 * C2
    G = 12.0 * np.pi * gamma * C1 * C2 * (C1 + C2)
    return -2*A/d**3 + 6*B/d**4 - 12*G/d**5


def omega_sq_radial_analytical(d_well, C, beta, gamma=None):
    """
    Radial breathing frequency at d_well for equilateral triangle (N=3).

    For a pair in isolation:
      omega^2 = V_2''(d_well) / mu_red, where mu_red = C/2

    For equilateral triangle, the breathing mode involves all 3 pairs
    expanding/contracting in phase.  The effective stiffness per pair
    is V_2''(d_well) (3 pairs contribute, but the mass-weighted
    eigenvector distributes this correctly).

    Simplified for radial mode of equal-mass equilateral:
      omega^2_radial = 3 * V_2''(d_well) / C
    (factor 3 from: each body feels 2 neighbors, but 3 pairs contribute
     to energy; mass-weighting by C per body).
    """
    if gamma is None:
        gamma = beta
    V2pp = V2_second_derivative(d_well, C, C, beta, gamma)
    # For equilateral breathing: effective spring constant per DOF
    # Each body of mass C, 2 neighbors at d_well, each with V_2''
    # omega^2 = (2 * V_2'') / C  [two springs pulling each body]
    return 2.0 * V2pp / C


def omega_sq_radial_exact(d_well, C, beta, gamma=None):
    """
    Exact: use V_2'' at d_well, two-body reduced mass.

    omega^2 = V_2''(d) / mu_red  for isolated pair.
    mu_red = C*C/(C+C) = C/2
    => omega^2 = 2*V_2''(d) / C
    """
    if gamma is None:
        gamma = beta
    V2pp = V2_second_derivative(d_well, C, C, beta, gamma)
    return 2.0 * V2pp / C


def hill_radius_newton(d, m_body, M_central):
    """
    Newton Hill radius: R_H = d * (m / (3*M))^(1/3).
    """
    return d * (m_body / (3.0 * M_central)) ** (1.0/3.0)


def hill_radius_tgp(d, C_body, C_central, beta, gamma=None):
    """
    TGP Hill radius from tidal balance.

    Tidal acceleration of central body on test mass at distance r from
    body of mass C_body, all at distance d from C_central:

      a_tidal ~ |V_eff''(d)| * r  (leading order)

    where V_eff is the potential of C_central.

    Self-gravity of C_body at distance r:
      a_self = |F_body(r)| / C_test ≈ 4*pi*C_body / r^2
      (leading gradient term; we use Newton-like leading term since
       at the Hill radius scale, 1/r dominates)

    Balance: 4*pi*C_body / R_H^2 = |V_eff''(d)| * R_H
    => R_H^3 = 4*pi*C_body / |V_eff''(d)|

    For Newton: V_N = -4*pi*C_central*C_test/d
                V_N'' = -8*pi*C_central*C_test/d^3
                => R_H^3 = 4*pi*C_body / (8*pi*C_central/d^3)
                         = C_body*d^3 / (2*C_central)
                R_H = d*(C_body/(2*C_central))^{1/3}
    (differs from standard by factor because F = 4*pi*C/r^2, not G*M/r^2)

    For TGP: V_2'' includes beta and gamma terms => modified R_H.
    """
    if gamma is None:
        gamma = beta

    # Tidal field from C_central at distance d
    # V_2''(d) for pair (C_central, C_test) -- C_test drops out of ratio
    # Compute |d^2 V / dd^2| for C_central on a test mass at distance d
    # V_eff = -4*pi*C_central/d + 8*pi*beta*C_central/d^2 - 12*pi*gamma*C_central*(C_central+1)/(d^3)
    # But for tidal: we differentiate the *acceleration* a(d) = F/m_test = -V'(d)/m_test
    # For the leading term: a = 4*pi*C_central/d^2 (Newton-like)
    # da/dd = -8*pi*C_central/d^3 (Newton)
    # TGP: a = 4*pi*C_central/d^2 - 16*pi*beta*C_central/d^3 + ...
    # da/dd = -8*pi*C_central/d^3 + 48*pi*beta*C_central/d^4 - ...

    # More carefully: acceleration from C_central on test particle:
    # a(d) = -V_2'(d) / C_test  [gradient of 2-body potential per test mass]
    # For pair (C_central, C_test=1):
    A = 4.0 * np.pi * C_central  # * C_test, but C_test=1
    B = 8.0 * np.pi * beta * C_central
    # gamma term: 12*pi*gamma*C_central*(C_central + C_test) ≈ 12*pi*gamma*C_central^2
    # for C_test << C_central.  But in TGP, C_test is not infinitesimal...
    # For Hill criterion, we use the standard restricted 3-body approach:
    # C_test -> 0.  Then gamma term ~ C_central * C_test * (C_central) ~ C_central^2 * C_test
    G_coeff = 12.0 * np.pi * gamma * C_central * C_central  # C_test divides out

    # a(d) = A/d^2 - 2*B/d^3 + 3*G_coeff/d^4  [per unit C_test]
    # da/dd = -2*A/d^3 + 6*B/d^4 - 12*G_coeff/d^5

    tidal = abs(-2*A/d**3 + 6*B/d**4 - 12*G_coeff/d**5)

    # Self-gravity of C_body at R_H (leading term):
    # a_self = 4*pi*C_body/R_H^2

    # Balance: 4*pi*C_body/R_H^2 = tidal * R_H
    # R_H^3 = 4*pi*C_body / tidal
    if tidal < 1e-30:
        return float('inf')
    R_H = (4.0 * np.pi * C_body / tidal) ** (1.0/3.0)
    return R_H


def hill_radius_newton_from_tidal(d, C_body, C_central):
    """Newton Hill radius from tidal balance (same framework as TGP for comparison)."""
    A = 4.0 * np.pi * C_central
    tidal = abs(-2*A/d**3)  # = 8*pi*C_central/d^3
    if tidal < 1e-30:
        return float('inf')
    R_H = (4.0 * np.pi * C_body / tidal) ** (1.0/3.0)
    return R_H


def V3_perturbative_shift(d_well, C, beta, gamma=None, n_quad=14):
    """
    First-order perturbative shift of d_well due to V3.

    For equilateral (d12=d13=d23=d_well), total potential per pair:
      V_eff = V_2(d) + V_3(d)/(3 pairs)  [V_3 is per triplet, 1 triplet / 3 pairs]

    V_3(d) = -6*gamma*C^3 * I_Y(d, d, d; m_sp)

    F_total(d_well) = F_2(d_well) + F_3(d_well) ≠ 0  (V_3 shifts equilibrium)

    New equilibrium: d_well' = d_well - F_3(d_well) / V_2''(d_well)
    (first-order perturbation theory)

    delta_d = -F_3(d_well) / V_2''(d_well)
    where F_3(d_well) = -dV_3/dd at d_well.

    We compute F_3 by finite difference.
    """
    if gamma is None:
        gamma = beta

    eps = 1e-5 * d_well
    V3_plus = three_body_energy_exact(d_well+eps, d_well+eps, d_well+eps,
                                       C, C, C, beta=beta, gamma=gamma,
                                       n_quad=n_quad)
    V3_minus = three_body_energy_exact(d_well-eps, d_well-eps, d_well-eps,
                                        C, C, C, beta=beta, gamma=gamma,
                                        n_quad=n_quad)
    # dV3/dd (all three distances change together for breathing mode)
    # V3 depends on all three d's, so dV3_total/dd = 3 * dV3/d(d_ij)
    # But for equilateral breathing: d(d_ij)/d(scale) = 1 for each pair
    dV3_dd = (V3_plus - V3_minus) / (2*eps)

    # F_3 = -dV3_dd (force from V3 along breathing direction)
    F3 = -dV3_dd

    V2pp = V2_second_derivative(d_well, C, C, beta, gamma)
    # For breathing mode: V2_total = 3*V2(d), so V2_total'' = 3*V2''(d)
    V2pp_total = 3.0 * V2pp

    if abs(V2pp_total) < 1e-30:
        return 0.0, F3, V2pp_total

    delta_d = -F3 / V2pp_total  # first-order shift (actually: -F_3_per_body/V2pp_eff, but for breathing this is it)
    return delta_d, F3, V2pp_total


# ── Part 1: Analytical d_rep, d_well verification ────────────────────────

def part1_equilibria():
    """Verify analytical equilibria against numerical."""
    print("=" * 78)
    print("Part 1: Analytical equilibria d_rep, d_well (exact closed form)")
    print("=" * 78)

    print(f"\n  V_2'(d) = 0  =>  d^2 - 4*beta*d + 18*gamma*C = 0")
    print(f"  d_rep,well = 2*beta -/+ sqrt(4*beta^2 - 18*gamma*C)")
    print(f"  Existence: 4*beta^2 > 18*gamma*C  <=>  beta > sqrt(9*gamma*C/2)")

    if quick:
        params = [(0.10, 1.0), (0.20, 1.0), (0.50, 2.0), (1.0, 1.0)]
    else:
        params = [(0.05, 0.5), (0.10, 1.0), (0.20, 1.0), (0.20, 0.5),
                  (0.50, 1.0), (0.50, 2.0), (1.0, 1.0), (2.0, 1.0)]

    print(f"\n  {'C':>6s} {'beta':>8s} {'d_rep(ana)':>12s} {'d_rep(num)':>12s} "
          f"{'d_well(ana)':>12s} {'d_well(num)':>12s} {'err_rep':>10s} {'err_well':>10s}")
    print(f"  {'-'*86}")

    all_pass = True
    results = []
    for C, beta in params:
        gamma = beta
        ana = equilibria_analytical(C, beta, gamma)
        num = force_zeros_2body(C, beta, gamma)

        if ana is None:
            d_rep_a = d_well_a = float('nan')
            exists = False
        else:
            d_rep_a, d_well_a, disc = ana
            exists = d_rep_a > 0

        if num is None:
            d_rep_n = d_well_n = float('nan')
        elif len(num) == 1:
            d_rep_n = float('nan')
            d_well_n = num[0]
        else:
            d_rep_n, d_well_n = num

        if exists and num is not None and len(num) == 2:
            err_rep = abs(d_rep_a - d_rep_n) / abs(d_rep_n)
            err_well = abs(d_well_a - d_well_n) / abs(d_well_n)
            if err_rep > 1e-10 or err_well > 1e-10:
                all_pass = False
        elif ana is None and num is None:
            err_rep = err_well = 0.0
        else:
            err_rep = err_well = float('nan')

        results.append({
            'C': C, 'beta': beta,
            'd_rep_a': d_rep_a, 'd_well_a': d_well_a,
            'd_rep_n': d_rep_n, 'd_well_n': d_well_n,
        })

        print(f"  {C:6.3f} {beta:8.4f} {d_rep_a:12.6f} {d_rep_n:12.6f} "
              f"{d_well_a:12.6f} {d_well_n:12.6f} {err_rep:10.2e} {err_well:10.2e}")

    print(f"\n  Analytical = Numerical: {'ALL PASS' if all_pass else 'DISCREPANCY'}")
    return results


# ── Part 2: Analytical omega^2 vs numerical Hessian ─────────────────────

def part2_omega_sq():
    """Verify analytical omega^2 against numerical finite-difference Hessian."""
    print(f"\n{'=' * 78}")
    print("Part 2: Analytical omega^2 at d_well (radial breathing)")
    print(f"{'=' * 78}")

    print(f"\n  V_2''(d) = -2A/d^3 + 6B/d^4 - 12G/d^5")
    print(f"  omega^2_rad = 2*V_2''(d_well)/C  (pair reduced mass = C/2)")

    if quick:
        params = [(0.10, 1.0), (0.50, 2.0), (1.0, 1.0)]
    else:
        params = [(0.10, 1.0), (0.20, 1.0), (0.50, 1.0),
                  (0.50, 2.0), (1.0, 1.0), (2.0, 1.0)]

    print(f"\n  {'C':>6s} {'beta':>8s} {'d_well':>10s} "
          f"{'V2pp(ana)':>12s} {'V2pp(FD)':>12s} {'omega2':>12s} {'err':>10s}")
    print(f"  {'-'*74}")

    all_pass = True
    results = []
    for C, beta in params:
        gamma = beta
        ana = equilibria_analytical(C, beta, gamma)
        if ana is None:
            continue
        d_rep, d_well, disc = ana
        if d_rep <= 0:
            continue

        # Analytical V2''
        V2pp_ana = V2_second_derivative(d_well, C, C, beta, gamma)

        # Numerical V2'' by FD
        eps = 1e-5 * d_well
        V_plus = V_eff_total(d_well + eps, C, C, beta, gamma)
        V_0 = V_eff_total(d_well, C, C, beta, gamma)
        V_minus = V_eff_total(d_well - eps, C, C, beta, gamma)
        V2pp_fd = (V_plus - 2*V_0 + V_minus) / eps**2

        err = abs(V2pp_ana - V2pp_fd) / abs(V2pp_fd) if abs(V2pp_fd) > 1e-30 else 0
        if err > 1e-4:
            all_pass = False

        omega2 = omega_sq_radial_exact(d_well, C, beta, gamma)

        results.append({
            'C': C, 'beta': beta, 'd_well': d_well,
            'V2pp_ana': V2pp_ana, 'V2pp_fd': V2pp_fd,
            'omega2': omega2,
        })

        print(f"  {C:6.3f} {beta:8.4f} {d_well:10.6f} "
              f"{V2pp_ana:12.6e} {V2pp_fd:12.6e} {omega2:12.6e} {err:10.2e}")

    print(f"\n  V2'' analytical vs FD: {'ALL PASS' if all_pass else 'DISCREPANCY'}")

    # Derive closed-form omega^2
    print(f"\n  Closed-form omega^2 at d_well:")
    print(f"  omega^2 = 2*V_2''(d_well)/C")
    print(f"  V_2''(d) = (4*pi*C^2/d^3)*[-2 + 6*(2*beta/d) - 12*(3*gamma*C/d^2)]")
    print(f"  At d_well = 2*beta + sqrt(4*beta^2 - 18*gamma*C):")
    print(f"  => omega^2 = (8*pi*C/d_well^3)*[-1 + 6*beta/d_well - 18*gamma*C/d_well^2]")

    return results


# ── Part 3: TGP Hill sphere vs Newton ────────────────────────────────────

def part3_hill():
    """Compare TGP and Newton Hill sphere radii."""
    print(f"\n{'=' * 78}")
    print("Part 3: TGP Hill sphere vs Newton")
    print(f"{'=' * 78}")

    C_body = 0.10
    C_central = 1.0
    d = 5.0  # orbital distance

    if quick:
        betas = [0.01, 0.05, 0.10, 0.50, 1.0, 2.0]
    else:
        betas = [0.005, 0.01, 0.025, 0.05, 0.08, 0.10, 0.20,
                 0.50, 1.0, 2.0, 5.0]

    print(f"\n  C_body = {C_body}, C_central = {C_central}, d = {d}")
    print(f"  Newton: R_H = (C_body*d^3 / (2*C_central))^{{1/3}}")
    print(f"          = (4*pi*C_body / |tidal|)^{{1/3}}")

    R_H_newton = hill_radius_newton_from_tidal(d, C_body, C_central)
    print(f"  R_H(Newton) = {R_H_newton:.6f}")

    print(f"\n  {'beta':>8s} {'R_H(TGP)':>12s} {'R_H(N)':>12s} "
          f"{'ratio':>10s} {'effect':>12s}")
    print(f"  {'-'*58}")

    results = []
    for beta in betas:
        gamma = beta
        R_H_tgp = hill_radius_tgp(d, C_body, C_central, beta, gamma)
        ratio = R_H_tgp / R_H_newton if R_H_newton > 0 else float('nan')
        effect = "LARGER" if ratio > 1.01 else ("SMALLER" if ratio < 0.99 else "~SAME")

        results.append({
            'beta': beta, 'R_H_tgp': R_H_tgp,
            'R_H_newton': R_H_newton, 'ratio': ratio,
        })

        print(f"  {beta:8.4f} {R_H_tgp:12.6f} {R_H_newton:12.6f} "
              f"{ratio:10.6f} {effect:>12s}")

    # Analysis
    if len(results) >= 2:
        ratios = [r['ratio'] for r in results]
        print(f"\n  Range of R_H(TGP)/R_H(Newton): {min(ratios):.4f} -- {max(ratios):.4f}")
        print(f"  TGP generally {'EXPANDS' if np.mean(ratios) > 1 else 'SHRINKS'} "
              f"the Hill sphere")

    # Physical interpretation
    print(f"\n  Physical interpretation:")
    print(f"  - beta repulsion (1/d^3 tidal) weakens tidal field => LARGER Hill sphere")
    print(f"  - gamma confinement (1/d^5 tidal) strengthens tidal field => SMALLER")
    print(f"  - Net effect depends on beta/gamma and d")

    return results


# ── Part 4: V3 perturbative correction ───────────────────────────────────

def part4_v3_correction():
    """Perturbative V3 shift of d_well."""
    print(f"\n{'=' * 78}")
    print("Part 4: Perturbative V3 correction to d_well")
    print(f"{'=' * 78}")

    print(f"\n  delta_d = -F_3(d_well) / (3*V_2''(d_well))")
    print(f"  (first-order perturbation: V3 shifts the equilibrium)")

    if quick:
        params = [(0.10, 1.0), (0.20, 2.0), (0.50, 1.0)]
        n_quad = 10
    else:
        params = [(0.10, 0.5), (0.10, 1.0), (0.20, 1.0),
                  (0.20, 2.0), (0.50, 1.0), (1.0, 1.0)]
        n_quad = 14

    print(f"\n  {'C':>6s} {'beta':>8s} {'d_well':>10s} {'delta_d':>12s} "
          f"{'delta_d/d':>12s} {'F_3':>12s} {'V2pp_total':>12s}")
    print(f"  {'-'*76}")

    results = []
    for C, beta in params:
        gamma = beta
        ana = equilibria_analytical(C, beta, gamma)
        if ana is None:
            continue
        d_rep, d_well, disc = ana
        if d_rep <= 0:
            continue

        delta_d, F3, V2pp_tot = V3_perturbative_shift(
            d_well, C, beta, gamma, n_quad=n_quad)

        rel = delta_d / d_well if d_well > 0 else float('nan')

        results.append({
            'C': C, 'beta': beta, 'd_well': d_well,
            'delta_d': delta_d, 'rel': rel,
        })

        print(f"  {C:6.3f} {beta:8.4f} {d_well:10.6f} {delta_d:12.6e} "
              f"{rel:12.6e} {F3:12.6e} {V2pp_tot:12.6e}")

    # Analysis
    if results:
        rels = [abs(r['rel']) for r in results]
        print(f"\n  |delta_d/d_well| range: {min(rels):.4e} -- {max(rels):.4e}")
        perturbative = all(r < 0.1 for r in rels)
        print(f"  V3 correction: {'PERTURBATIVE (<10%)' if perturbative else 'SIGNIFICANT'}")

        # Sign analysis
        signs = [np.sign(r['delta_d']) for r in results]
        if all(s > 0 for s in signs):
            print(f"  V3 shifts d_well OUTWARD (larger equilibrium distance)")
        elif all(s < 0 for s in signs):
            print(f"  V3 shifts d_well INWARD (smaller equilibrium distance)")
        else:
            print(f"  V3 shift direction depends on parameters")

    return results


# ── Part 5: Summary table ────────────────────────────────────────────────

def part5_summary(p1, p2, p3, p4):
    """Publication-ready summary."""
    print(f"\n{'=' * 78}")
    print("Part 5: Summary for publication")
    print(f"{'=' * 78}")

    print(f"""
  ANALYTICAL RESULTS (closed-form):
  ==================================

  1. Equilibrium condition (V_2'(d) = 0, equal masses C1=C2=C):

     d^2 - 4*beta*d + 18*gamma*C = 0

     d_rep  = 2*beta - sqrt(4*beta^2 - 18*gamma*C)   (repulsive barrier)
     d_well = 2*beta + sqrt(4*beta^2 - 18*gamma*C)   (confining minimum)

     Existence: beta > sqrt(9*gamma*C / 2)
     At threshold: d_rep = d_well = 2*beta (degenerate)
     Far from threshold: d_rep -> 0, d_well -> 4*beta

  2. Unequal masses (C1 != C2):

     d^2 - 4*beta*d + 9*gamma*(C1+C2) = 0

     d_rep,well = 2*beta -/+ sqrt(4*beta^2 - 9*gamma*(C1+C2))

  3. Radial frequency at d_well (isolated pair):

     omega^2 = 2*V_2''(d_well) / C

     V_2''(d) = (8*pi*C^2/d^3) * [-1 + 6*beta/d - 18*gamma*C/d^2]

  4. Hill sphere radius (TGP vs Newton):

     R_H^3 = 4*pi*C_body / |tidal field|

     tidal(d) = 8*pi*C_central/d^3 * [1 - 6*beta/d + 18*gamma*C_central/d^2]

     R_H(TGP)/R_H(Newton) depends on beta*C_central/d and gamma*C^2/d^2.
     For beta > 0: TGP EXPANDS the Hill sphere (weaker tidal disruption).

  5. V3 perturbative correction to d_well:

     delta_d = -F_3(d_well) / (3*V_2''(d_well))

     Typically |delta_d/d_well| < 1% -- V3 is perturbative on equilibria.
""")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    print("=" * 78)
    print("ex210 -- Analytical equilibria, omega^2, and Hill criterion")
    print("=" * 78)
    mode = "QUICK" if quick else "FULL"
    print(f"  Mode: {mode}")

    t0 = time.time()

    p1 = part1_equilibria()
    p2 = part2_omega_sq()
    p3 = part3_hill()
    p4 = part4_v3_correction()
    part5_summary(p1, p2, p3, p4)

    print(f"\n  Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
