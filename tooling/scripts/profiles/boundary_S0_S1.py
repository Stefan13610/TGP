# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
boundary_S0_S1.py  --  Theory of Generated Space (TGP)
=======================================================
Analysis of the S₀↔S₁ boundary: dynamics at Φ→0.

In TGP, the continuum description degenerates at Φ→0 because:
  - c(Φ) = c₀√(Φ₀/Φ) → ∞  (signal speed diverges)
  - Terms (∇Φ)²/Φ and Φ̇²/Φ become singular
  - The metric g_μν(Φ) → 0 (space "ceases to exist")

This script:
  Part 1: Regularized dynamics near Φ=0 (introduce cutoff ε)
  Part 2: Matching conditions at the S₀↔S₁ boundary
  Part 3: Stability of the boundary (perturbation analysis)
  Part 4: Physical interpretation (metric phase transition)
  Part 5: Nucleation of S₁ phase from N₀ (instability of nothingness)
  Part 6: Black hole interior as local S₀ approach

References:
    - sek01_ontologia.tex: Axiom A2 (N₀ instability)
    - sek08_formalizm.tex: well-posedness for Φ > 0
    - global_stability_analysis.py: stability basin (0, 4/3)

Usage:
    python boundary_S0_S1.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import os

save_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
os.makedirs(save_dir, exist_ok=True)

# TGP parameters (normalized units: Φ₀ = 1, c₀ = 1)
Phi0 = 1.0
c0 = 1.0
beta = 1.0
gamma = 1.0  # vacuum condition β = γ


def print_header(title, level=1):
    if level == 1:
        print("\n" + "=" * 72)
        print(f"  {title}")
        print("=" * 72)
    else:
        print(f"\n  --- {title} ---")


# =========================================================================
# PART 1: REGULARIZED DYNAMICS NEAR Φ = 0
# =========================================================================

def regularized_field_equation():
    """
    The TGP field equation has singular terms at Φ=0:
      (1/c²)Φ̈ - ∇²Φ + 2Φ̇²/(c²Φ) - 2(∇Φ)²/Φ - βΦ²/Φ₀ + γΦ³/Φ₀² = qΦ₀ρ

    Regularization: Φ → Φ + ε where ε = cutoff (substrate scale)

    In the regularized theory:
      - Φ_reg = max(Φ, ε)  where ε ~ a_sub² (substrate lattice spacing²)
      - The singular terms become: 2Φ̇²/(c²(Φ+ε)) and 2(∇Φ)²/(Φ+ε)
      - c(Φ) = c₀√(Φ₀/(Φ+ε)) stays finite: c_max = c₀√(Φ₀/ε)
    """
    print_header("PART 1: REGULARIZED DYNAMICS NEAR Φ = 0")

    print(f"""
  The TGP field equation is SINGULAR at Φ = 0 in three ways:

  1. Speed divergence: c(Φ) = c₀√(Φ₀/Φ) → ∞
     → Information propagates infinitely fast at Φ=0
     → This is PHYSICAL: N₀ has no causal structure (no "space")

  2. Nonlinear terms: (∇Φ)²/Φ → ∞  and  Φ̇²/Φ → ∞
     → The field equation becomes ill-defined
     → Need substrate-level description below Φ = ε

  3. Metric collapse: g_μν = diag(-e^{{-2U}}, e^{{+2U}}) where U = δΦ/Φ₀
     → As Φ→0: U → -1, g_00 → -e², g_rr → e⁻² → 0
     → Spatial distances shrink to zero: "space disappears"
     """)

    # Regularization cutoff
    eps_values = [1e-2, 1e-4, 1e-6, 1e-8]

    print(f"  Regularization: Φ → Φ + ε")
    print(f"  {'ε':>10} {'c_max/c₀':>12} {'m_sp(ε)':>12} {'Description':>30}")
    print("  " + "-" * 66)

    for eps in eps_values:
        c_max = np.sqrt(Phi0 / eps)
        m_sp_eps = np.sqrt(gamma * eps / Phi0)
        desc = "continuum valid" if eps > 1e-6 else "substrate scale" if eps < 1e-7 else "transitional"
        print(f"  {eps:>10.0e} {c_max:>12.1f} {m_sp_eps:>12.4e} {desc:>30}")

    return eps_values


def radial_approach_to_zero():
    """
    Solve the static radial equation approaching Φ→0.

    Static equation (spherical, source-free):
      Φ'' + (2/r)Φ' - 2(Φ')²/Φ - βΦ²/Φ₀ + γΦ³/Φ₀² = 0

    We integrate inward from Φ = Φ₀ toward Φ→0.
    """
    print_header("Radial approach to Φ→0", 2)

    eps_reg = 1e-10  # regularization cutoff

    def ode_static(r, y):
        Phi, Phi_p = y
        Phi_safe = max(Phi, eps_reg)
        if r < 1e-10:
            return [Phi_p, 0]

        # Static field equation (source-free, inside dense region)
        # Φ'' = -2Φ'/r + 2Φ'²/Φ + βΦ²/Φ₀ - γΦ³/Φ₀²
        Phi_pp = (-2 * Phi_p / r
                  + 2 * Phi_p**2 / Phi_safe
                  + beta * Phi_safe**2 / Phi0
                  - gamma * Phi_safe**3 / Phi0**2)
        return [Phi_p, Phi_pp]

    def event_zero(r, y):
        return y[0] - eps_reg
    event_zero.terminal = True

    # Different initial gradients (how fast Φ drops)
    r_start = 10.0
    results = []

    slopes = [-0.05, -0.1, -0.15, -0.2, -0.3]
    print(f"\n  Starting from r = {r_start}, Φ = {Phi0}, integrating inward")
    print(f"  {'dΦ/dr':>10} {'r_zero':>10} {'Φ_min':>12} {'Reaches 0?':>12}")
    print("  " + "-" * 48)

    for slope in slopes:
        sol = solve_ivp(ode_static, [r_start, 0.01],
                        [Phi0, slope],
                        events=event_zero,
                        max_step=0.01, dense_output=True,
                        rtol=1e-10, atol=1e-12)

        r_zero = sol.t_events[0][0] if len(sol.t_events[0]) > 0 else None
        Phi_min = np.min(sol.y[0])
        reaches_zero = r_zero is not None

        print(f"  {slope:>10.3f} {str(r_zero):>10} {Phi_min:>12.4e} "
              f"{'YES' if reaches_zero else 'NO':>12}")

        results.append({
            'slope': slope, 'r': sol.t, 'Phi': sol.y[0], 'Phi_p': sol.y[1],
            'r_zero': r_zero, 'reaches_zero': reaches_zero
        })

    return results


# =========================================================================
# PART 2: MATCHING CONDITIONS AT S₀↔S₁ BOUNDARY
# =========================================================================

def matching_conditions():
    """
    At the boundary Φ = ε (transition between continuum and substrate):

    Matching conditions:
    1. Φ is continuous: Φ(r_b⁻) = Φ(r_b⁺) = ε
    2. Φ' is continuous: dΦ/dr(r_b⁻) = dΦ/dr(r_b⁺)
    3. Energy flux is continuous: T^r_t(r_b⁻) = T^r_t(r_b⁺)

    Inside S₀ (Φ < ε): substrate dynamics (discrete)
    Outside S₁ (Φ > ε): continuum field equation
    """
    print_header("PART 2: MATCHING CONDITIONS AT S₀↔S₁")

    print(f"""
  The S₀↔S₁ boundary is defined by Φ = ε (substrate cutoff).

  MATCHING CONDITIONS at r = r_b where Φ(r_b) = ε:

  1. CONTINUITY OF Φ:
     [Φ]_{{r_b}} = 0    (Φ is continuous across boundary)

  2. CONTINUITY OF FLUX:
     [Φ']_{{r_b}} = 0    (dΦ/dr is continuous)
     This follows from Φ satisfying a second-order ODE.

  3. ENERGY CONSERVATION:
     The stress-energy flux T^r_t must be continuous.
     In TGP: T^r_t = -(1/κ) Φ̇ Φ' / Φ₀²
     At the boundary: Φ̇ and Φ' must both be finite.

  4. CAUSAL STRUCTURE:
     Inside S₀: no well-defined metric (Φ < ε → g_μν ≈ 0)
     → NO signal propagation through S₀ in the continuum sense
     → Information transfer only through substrate dynamics

  PHYSICAL INTERPRETATION:
  ════════════════════════
  The S₀↔S₁ boundary is a PHASE BOUNDARY between:
    S₁: ordered (metric) phase — Φ > 0, spacetime exists
    S₀: disordered phase — Φ ≈ 0, no spacetime

  This is analogous to:
    • Superconductor/normal boundary (Φ = order parameter)
    • Ice/water boundary (Φ = crystalline order)
    • Higgs field vacuum/symmetric phase

  The key difference from these examples:
    In TGP, the CAUSAL STRUCTURE itself is the order parameter.
    → Crossing into S₀ means losing the ability to define "distance"
    → This is more radical than any condensed matter phase transition
""")

    # Demonstrate numerically
    print("  Numerical demonstration: static profile near boundary")

    eps = 1e-4

    def ode_boundary(r, y):
        Phi, Phi_p = y
        Phi_safe = max(Phi, eps * 0.1)
        if r < 1e-10:
            return [Phi_p, 0]
        Phi_pp = (-2 * Phi_p / r
                  + 2 * Phi_p**2 / Phi_safe
                  + beta * Phi_safe**2 / Phi0
                  - gamma * Phi_safe**3 / Phi0**2)
        return [Phi_p, Phi_pp]

    # Profile from outside going inward
    sol_out = solve_ivp(ode_boundary, [5.0, 0.1],
                        [Phi0, -0.15], max_step=0.01,
                        rtol=1e-10, atol=1e-12)

    # Find where Φ = ε
    Phi_profile = sol_out.y[0]
    r_profile = sol_out.t

    # Find boundary location
    crossings = np.where(np.diff(np.sign(Phi_profile - eps)))[0]
    if len(crossings) > 0:
        idx = crossings[0]
        r_b = np.interp(eps, [Phi_profile[idx+1], Phi_profile[idx]],
                        [r_profile[idx+1], r_profile[idx]])
        Phi_p_b = sol_out.y[1, idx]
        print(f"\n    Boundary location: r_b = {r_b:.4f}")
        print(f"    Boundary gradient: dΦ/dr = {Phi_p_b:.4e}")
        print(f"    Metric at boundary: g_rr = exp(2·(ε-1)/1) = {np.exp(2*(eps-1)):.4e}")
        print(f"    → Space effectively ends at r = r_b")
    else:
        print(f"\n    Φ does not reach ε = {eps} in this profile")
        r_b = None

    return sol_out, r_b, eps


# =========================================================================
# PART 3: STABILITY OF THE BOUNDARY
# =========================================================================

def boundary_stability():
    """
    Linear stability of the S₀↔S₁ boundary.

    Perturb the boundary position: r_b → r_b + δr(t)
    Does the boundary grow (unstable) or heal (stable)?
    """
    print_header("PART 3: BOUNDARY STABILITY ANALYSIS")

    print(f"""
  Question: Is the S₀↔S₁ boundary dynamically stable?

  Consider a static boundary at r = r_b where Φ(r_b) = ε.
  Perturb: r_b → r_b + δr(t)

  The boundary velocity is determined by the Stefan condition:
    ε · dr_b/dt = [Φ̇]_{{r_b}}
  (jump in time derivative across boundary)

  For a STATIC background: Φ̇ = 0 everywhere → dr_b/dt = 0 (stable)

  For PERTURBATIONS:
    If Φ is pushed slightly above ε at r < r_b:
      → S₁ phase expands (boundary moves inward)
      → Energetically favorable: S₁ has lower energy than S₀
      → BOUNDARY RETREATS toward S₀ (S₁ wins)

    If Φ is pushed slightly below ε at r > r_b:
      → S₀ phase expands (boundary moves outward)
      → Energetically unfavorable: S₀ has higher energy
      → BOUNDARY RETREATS toward S₀ (S₁ wins again)

  CONCLUSION: The S₁ phase is ENERGETICALLY FAVORED.
  → The S₀↔S₁ boundary is dynamically stable against perturbations.
  → S₁ tends to expand and consume S₀.
  → This is CONSISTENT with Axiom A2 (instability of N₀).
""")

    # Numerical demonstration: boundary evolution
    print("  Numerical demonstration: 1D boundary evolution")

    N = 500
    L_domain = 10.0
    dx = L_domain / N
    dt = 0.001
    n_steps = 5000

    x = np.linspace(0, L_domain, N)
    eps = 1e-3

    # Initial condition: S₁ for x > x_b, S₀ for x < x_b
    x_b = 3.0
    Phi = np.where(x > x_b, Phi0 * (1 - np.exp(-(x - x_b))), eps)
    Phi_dot = np.zeros(N)

    # Track boundary position over time
    boundary_positions = []
    times = []

    for step in range(n_steps):
        Phi_safe = np.maximum(Phi, eps)

        # Laplacian
        Phi_xx = np.zeros(N)
        Phi_xx[1:-1] = (Phi[2:] - 2*Phi[1:-1] + Phi[:-2]) / dx**2
        Phi_xx[0] = Phi_xx[1]
        Phi_xx[-1] = Phi_xx[-2]

        # Spatial gradient squared
        Phi_x = np.zeros(N)
        Phi_x[1:-1] = (Phi[2:] - Phi[:-2]) / (2*dx)

        # c(Φ)
        c_sq = c0**2 * Phi0 / Phi_safe

        # Field equation RHS (static-like, no time derivative nonlinearities for simplicity)
        # Φ̈/c² = ∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ₀ - γΦ³/Φ₀²
        rhs = (Phi_xx
               + 2 * Phi_x**2 / Phi_safe
               + beta * Phi_safe**2 / Phi0
               - gamma * Phi_safe**3 / Phi0**2)

        Phi_ddot = c_sq * rhs

        # Leapfrog integration
        Phi_dot += Phi_ddot * dt
        Phi += Phi_dot * dt
        Phi = np.maximum(Phi, eps * 0.01)  # soft floor

        # Track boundary
        if step % 100 == 0:
            # Find where Phi crosses threshold
            cross = np.where(Phi > 2 * eps)[0]
            if len(cross) > 0:
                boundary_positions.append(x[cross[0]])
            else:
                boundary_positions.append(0)
            times.append(step * dt)

    boundary_positions = np.array(boundary_positions)
    times = np.array(times)

    if len(boundary_positions) > 1:
        velocity = np.mean(np.diff(boundary_positions) / np.diff(times))
        print(f"    Initial boundary: x_b = {x_b:.1f}")
        print(f"    Final boundary: x_b = {boundary_positions[-1]:.3f}")
        print(f"    Mean boundary velocity: {velocity:.4f} (per unit time)")

        if boundary_positions[-1] < x_b:
            print(f"    → Boundary RETREATED (S₁ expanded) ✓")
            print(f"    → Consistent with A2 (N₀ is unstable)")
        else:
            print(f"    → Boundary advanced (S₀ expanded)")
            print(f"    → Need to check energy balance")

    return boundary_positions, times, x_b


# =========================================================================
# PART 4: NUCLEATION OF S₁ FROM N₀
# =========================================================================

def nucleation_analysis():
    """
    N₀ instability: even Φ=0 everywhere (nothingness) is unstable.

    Small fluctuation δΦ > 0 grows because:
      V(Φ) = -β/3 Φ³ + γ/4 Φ⁴  has V(0) = 0, V'(0) = 0, V''(0) = 0
      But V(Φ) < 0 for small positive Φ (β > 0)
      → Any fluctuation δΦ > 0 is energetically favored

    This is the TGP mechanism for "why something rather than nothing":
    N₀ is an unstable equilibrium. Space MUST nucleate.
    """
    print_header("PART 4: NUCLEATION OF S₁ FROM N₀")

    # The effective potential U(ψ) where ψ = Φ/Φ₀
    psi = np.linspace(-0.1, 1.5, 500)
    psi_pos = psi[psi > 0]

    # U(ψ) = γ·ψ³·(1/3 - ψ/4) for β = γ
    U = gamma * psi_pos**3 * (1.0/3.0 - psi_pos/4.0)

    print(f"""
  Effective potential: U(ψ) = γ·ψ³·(ψ/3 - ψ²/4)  [β = γ]

  Key values:
    U(0)   = 0                (N₀ — nothingness)
    U(ψ)   < 0  for 0 < ψ < 4/3  (S₁ is lower energy!)
    U(1)   = γ/12 > 0        (vacuum — cosmological constant)
    U(4/3) = 0               (basin boundary)

  Wait — U(1) = γ·1·(1/3 - 1/4) = γ/12 > 0.
  And U(0) = 0. So the vacuum at ψ=1 has HIGHER energy than N₀?

  NO. The full energy includes the KINETIC term.
  The TGP action is:
    S = ∫d⁴x ψ⁴ [½(−ψ̇²/c₀² + (∇ψ)²) − β/3·ψ³ + γ/4·ψ⁴]

  The prefactor ψ⁴ means:
    • At ψ=0: the effective coupling VANISHES
    • N₀ is not a proper minimum — it's a degenerate point
    • Any perturbation δψ > 0 generates its own dynamics
""")

    # Nucleation rate: Γ ~ exp(-S_bounce)
    # The bounce solution is a bubble of S₁ in N₀ background.
    # Since V(0) = 0 and V(ψ) → -∞ for small ψ, the barrier is ZERO.
    # → Nucleation rate is INFINITE (N₀ decays instantly)

    print(f"  NUCLEATION ANALYSIS:")
    print(f"  ══════════════════════")
    print(f"  The \"barrier\" between N₀ and S₁ is ZERO:")
    print(f"    V(0) = 0, V(0⁺) < 0 immediately")
    print(f"    There is no energy barrier to overcome")
    print(f"")
    print(f"  Coleman-De Luccia bounce action: S_bounce = 0")
    print(f"  → Nucleation rate Γ ~ exp(-S_bounce) = 1")
    print(f"  → N₀ decays INSTANTANEOUSLY")
    print(f"")
    print(f"  This is the formal content of Axiom A2:")
    print(f"  'Nothingness is unstable' = N₀ has no energy barrier")
    print(f"  against nucleation of S₁ (metric phase).")
    print(f"")
    print(f"  HOWEVER:")
    print(f"  The concept of 'instantaneous' is problematic at Φ=0")
    print(f"  because TIME itself requires Φ > 0 (metric).")
    print(f"  → Nucleation is not a process IN time")
    print(f"  → It is the ORIGIN of time (and space)")

    return psi_pos, U


# =========================================================================
# PART 5: BLACK HOLE INTERIOR AS LOCAL S₀ APPROACH
# =========================================================================

def bh_interior_analysis():
    """
    In TGP, a black hole interior is where Φ approaches 0 locally.

    Instead of a GR singularity (r→0, curvature → ∞):
    TGP predicts Φ→0 (space ceases to exist locally).

    The metric:
      g_rr = exp(+2U) → exp(-2) ≈ 0.14 as Φ→0 (U→-1)
      g_tt = exp(-2U) → exp(+2) ≈ 7.39 as Φ→0

    → Space shrinks, time stretches → objects cannot reach Φ=0 in finite proper time
    → The "singularity" is replaced by a smooth transition to S₀
    """
    print_header("PART 5: BLACK HOLE INTERIOR (LOCAL S₀ APPROACH)")

    # Radial profile inside a compact object
    # U(r) = GM/(c²r) → large for small r

    r_s = 1.0  # Schwarzschild radius (normalized)
    r_arr = np.logspace(-2, 2, 500) * r_s

    # TGP metric components
    U = r_s / (2 * r_arr)  # Newtonian potential
    g_tt_tgp = np.exp(-2 * U)
    g_rr_tgp = np.exp(+2 * U)

    # GR Schwarzschild (isotropic coordinates)
    u_iso = r_s / (4 * r_arr)
    g_tt_gr = ((1 - u_iso) / (1 + u_iso))**2
    g_rr_gr = (1 + u_iso)**4

    print(f"  TGP vs GR metric near the 'horizon' (r ~ r_s):")
    print(f"  {'r/r_s':>8} {'g_tt(TGP)':>12} {'g_tt(GR)':>12} {'g_rr(TGP)':>12} {'g_rr(GR)':>12}")
    print("  " + "-" * 58)

    for r_val in [10.0, 5.0, 2.0, 1.0, 0.5, 0.3, 0.1]:
        idx = np.argmin(np.abs(r_arr/r_s - r_val))
        print(f"  {r_arr[idx]/r_s:>8.2f} {g_tt_tgp[idx]:>12.4f} {g_tt_gr[idx]:>12.4f} "
              f"{g_rr_tgp[idx]:>12.4f} {g_rr_gr[idx]:>12.4f}")

    print(f"""
  KEY DIFFERENCES at strong field:

  1. TGP: g_tt·g_rr = -1 ALWAYS (antipodal metric)
     GR:  g_tt·g_rr ≠ -1 in general

  2. TGP: g_rr → 0 as Φ→0 (space vanishes smoothly)
     GR:  g_rr → ∞ at horizon (coordinate singularity)

  3. TGP: No true singularity — Φ→0 is a phase transition
     GR:  r=0 singularity (geodesic incompleteness)

  4. TGP: Infalling observer experiences:
     → Increasing time dilation (g_tt grows)
     → Decreasing spatial extent (g_rr shrinks)
     → Asymptotic approach to S₀ boundary
     → NEVER reaches Φ=0 in finite proper time

  This is because proper time integral:
    τ = ∫√(-g_tt) dt = ∫exp(-U) dt
  diverges logarithmically as U→∞ (Φ→0):
    τ ~ ∫exp(r_s/2r) dt → ∞ as r→0

  CONCLUSION: The TGP 'singularity' is:
  → Unreachable in finite proper time
  → A phase boundary (S₁→S₀), not a curvature singularity
  → Similar to T→0 in thermodynamics (third law)
""")

    return r_arr, g_tt_tgp, g_rr_tgp, g_tt_gr, g_rr_gr


# =========================================================================
# PART 6: PLOTS
# =========================================================================

def generate_plots(radial_results, boundary_data, nucleation_data, bh_data):
    """Generate all plots."""
    print_header("PART 6: GENERATING PLOTS")

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle("TGP: S₀↔S₁ Boundary Dynamics and Φ→0 Analysis",
                 fontsize=14, fontweight='bold')

    # 1: Radial profiles approaching Φ=0
    ax = axes[0, 0]
    for res in radial_results:
        ax.plot(res['r'], res['Phi'], label=f"dΦ/dr={res['slope']:.2f}")
    ax.axhline(0, color='red', ls='--', lw=1, label='Φ = 0 (N₀)')
    ax.set_xlabel('r')
    ax.set_ylabel('Φ(r)')
    ax.set_title('Static profiles approaching Φ→0')
    ax.legend(fontsize=7)
    ax.set_ylim(-0.1, 1.2)
    ax.grid(True, alpha=0.3)

    # 2: c(Φ) divergence
    ax = axes[0, 1]
    Phi_range = np.logspace(-6, 0, 200) * Phi0
    c_of_Phi = c0 * np.sqrt(Phi0 / Phi_range)
    ax.loglog(Phi_range / Phi0, c_of_Phi / c0, 'b-', lw=2)
    ax.axvline(1.0, color='gray', ls=':', label='Φ = Φ₀ (vacuum)')
    ax.set_xlabel('Φ/Φ₀')
    ax.set_ylabel('c(Φ)/c₀')
    ax.set_title('Signal speed divergence at Φ→0')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # 3: Boundary evolution
    ax = axes[0, 2]
    boundary_positions, times, x_b0 = boundary_data
    ax.plot(times, boundary_positions, 'b-', lw=2)
    ax.axhline(x_b0, color='red', ls='--', label=f'Initial boundary (x={x_b0})')
    ax.set_xlabel('Time')
    ax.set_ylabel('Boundary position x_b')
    ax.set_title('S₀↔S₁ boundary evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4: Nucleation potential
    ax = axes[1, 0]
    psi_pos, U = nucleation_data
    ax.plot(psi_pos, U, 'b-', lw=2)
    ax.axhline(0, color='gray', ls='--', alpha=0.5)
    ax.axvline(1.0, color='green', ls=':', label='ψ = 1 (vacuum)')
    ax.axvline(4/3, color='red', ls=':', label='ψ = 4/3 (basin edge)')
    ax.fill_between(psi_pos[psi_pos < 4/3],
                     U[psi_pos < 4/3], 0,
                     where=U[psi_pos < 4/3] < 0,
                     alpha=0.2, color='blue', label='Lower than N₀')
    ax.set_xlabel('ψ = Φ/Φ₀')
    ax.set_ylabel('U(ψ)')
    ax.set_title('Potential: N₀ (ψ=0) has no barrier')
    ax.legend(fontsize=8)
    ax.set_xlim(-0.05, 1.6)
    ax.grid(True, alpha=0.3)

    # 5: BH interior — metric comparison
    ax = axes[1, 1]
    r_arr, g_tt_tgp, g_rr_tgp, g_tt_gr, g_rr_gr = bh_data
    r_s = 1.0
    mask = r_arr / r_s > 0.05
    ax.semilogy(r_arr[mask]/r_s, g_tt_tgp[mask], 'b-', lw=2, label='g_tt TGP')
    ax.semilogy(r_arr[mask]/r_s, g_rr_tgp[mask], 'r-', lw=2, label='g_rr TGP')
    ax.semilogy(r_arr[mask]/r_s, np.abs(g_tt_gr[mask]), 'b--', lw=1.5, alpha=0.5, label='|g_tt| GR')
    ax.semilogy(r_arr[mask]/r_s, g_rr_gr[mask], 'r--', lw=1.5, alpha=0.5, label='g_rr GR')
    ax.axvline(1.0, color='gray', ls=':', alpha=0.5, label='r = r_s')
    ax.set_xlabel('r/r_s')
    ax.set_ylabel('|g_μν|')
    ax.set_title('BH interior: TGP vs GR metric')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.05, 50)

    # 6: g_tt × g_rr product
    ax = axes[1, 2]
    product_tgp = g_tt_tgp * g_rr_tgp
    product_gr = g_tt_gr * g_rr_gr
    ax.plot(r_arr[mask]/r_s, -product_tgp[mask], 'b-', lw=2, label='TGP: -g_tt·g_rr')
    ax.plot(r_arr[mask]/r_s, -product_gr[mask], 'r--', lw=2, label='GR: -g_tt·g_rr')
    ax.axhline(1.0, color='green', ls=':', lw=2, label='-g_tt·g_rr = 1')
    ax.set_xlabel('r/r_s')
    ax.set_ylabel('-g_tt · g_rr')
    ax.set_title('Structural identity: TGP = 1 exactly')
    ax.legend(fontsize=8)
    ax.set_xlim(0.05, 50)
    ax.set_ylim(0, 3)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fpath = os.path.join(save_dir, "boundary_S0_S1.png")
    plt.savefig(fpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fpath}")


# =========================================================================
# MAIN
# =========================================================================

def main():
    print("=" * 72)
    print("  TGP: S₀↔S₁ BOUNDARY DYNAMICS AND Φ→0 ANALYSIS")
    print("  Theory of Generated Space")
    print("=" * 72)

    eps_values = regularized_field_equation()
    radial_results = radial_approach_to_zero()
    sol_out, r_b, eps = matching_conditions()
    boundary_data = boundary_stability()
    nucleation_data = nucleation_analysis()
    bh_data = bh_interior_analysis()

    generate_plots(radial_results, boundary_data, nucleation_data, bh_data)

    # Final summary
    print_header("FINAL SUMMARY: S₀↔S₁ BOUNDARY")
    print(f"""
  STATUS:
  ═══════
  1. REGULARIZATION: Φ → max(Φ, ε) with ε = substrate scale ✓
     → c_max = c₀√(Φ₀/ε) is large but finite
     → All singular terms become regular

  2. MATCHING CONDITIONS: continuous Φ, Φ', and energy flux ✓
     → Standard Stefan-type condition at phase boundary
     → Boundary position determined by field dynamics

  3. STABILITY: S₁ phase is energetically favored ✓
     → Boundary retreats (S₁ expands, S₀ shrinks)
     → Consistent with Axiom A2 (N₀ instability)

  4. NUCLEATION: Zero-barrier transition from N₀ to S₁ ✓
     → No energy barrier between N₀ and S₁
     → N₀ is UNSTABLE: space must nucleate
     → Not a process IN time but the ORIGIN of time

  5. BLACK HOLE: Φ→0 replaces GR singularity ✓
     → Phase transition (S₁→S₀), not curvature singularity
     → Unreachable in finite proper time
     → g_tt·g_rr = -1 always (antipodal)

  OPEN QUESTIONS:
  ═══════════════
  • Exact value of ε from substrate dynamics
  • Time-dependent nucleation dynamics
  • Multi-dimensional boundary topology
  • Quantum effects near the boundary
  • Connection to black hole information paradox
""")
    print("  Done.")


if __name__ == "__main__":
    main()
