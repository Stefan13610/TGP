"""
M11.S — Mini-PoC: single soliton Φ_sol(r) feasibility check
================================================================

Goal: BEFORE committing to full M11.S audit, verify that a non-trivial
classical configuration Φ_sol(r) sourced by a point mass exists, that
its asymptotic structure matches M9.3.1 Yukawa linearization, and that
its linearization spectrum is well-defined (basic stability indicator).

Tests:
  PoC.1 — Linear Yukawa benchmark (M9.3.1 spatial linearization)
          (∇² - μ²)δφ = -(q/Φ_0)ρ/K_geo,  μ² = β/K_geo
          Verifies sign convention consistent with M9.3.1.

  PoC.2 — Nonlinear ODE solver Φ_sol(r) with regularized point source
          K(φ)[φ'' + (2/r)φ'] + ½K'(φ)(φ')² + V'(φ) + (q/Φ_0)ρ(r) = 0
          From sek08a action; β=γ vacuum; M9.1'' hyperbolic-effective sign.

  PoC.3 — Linear vs nonlinear comparison
          Verify nonlinear correction is small at large r (Yukawa tail OK),
          potentially significant near core.

  PoC.4 — Linearization spectrum (radial s-wave) around Φ_sol
          Sturm-Liouville problem: -d/dr[r² K(Φ_sol) du/dr] + r² M_eff²(r) u
                                   = ω² r² K(Φ_sol) u
          Verify: lowest eigenvalues finite, no obvious instability
                  in s-wave sector (full spectrum is M11.S audit job).

  PoC.5 — H1/H2/H3 verdict (existence, asymptotic, basic stability)

Pre-test hypotheses (from M11_branch_strategy.md §2.1):
  H1: Φ_sol(r) nontrivial classical solution exists numerically
  H2: Φ_sol(0+) regular, Φ_sol(∞) → Φ_0 = 1
  H3: Spectrum {ω²_n} ≥ 0 in s-wave (modulo zero modes from translations)

Units: β = 1, K_geo = 1, Φ_0 = 1, q/Φ_0 = 1
       so Compton range λ_C = √(K_geo/β) = 1.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
import sys


# ============================================================================
# CONSTANTS (TGP, in PoC units)
# ============================================================================

beta = 1.0          # V''(Φ_0) magnitude (λ_C^-2 in K_geo=1 units)
K_geo = 1.0         # kinetic coefficient
Phi_0 = 1.0         # cosmological vacuum
q_over_Phi0 = 1.0   # source coupling
M_source = 1.0      # point-source "mass" amplitude
a_source = 0.15     # Gaussian regularization (small fraction of λ_C = 1)

# Numerical grid
r_min = 1e-3
r_max = 12.0        # ~12 λ_C — well into asymptotic regime


# ============================================================================
# POTENTIAL & KINETIC (sek08a, β=γ vacuum)
# ============================================================================

def V(phi):
    """V(φ) = (β/3)φ³ - (γ/4)φ⁴, β=γ"""
    return (beta/3.0) * phi**3 - (beta/4.0) * phi**4

def Vp(phi):
    """V'(φ) = β·φ²(1-φ)"""
    return beta * phi**2 * (1.0 - phi)

def Vpp(phi):
    """V''(φ) = β·φ(2-3φ)"""
    return beta * phi * (2.0 - 3.0*phi)

def K(phi):
    """K(φ) = K_geo · φ⁴"""
    return K_geo * phi**4

def Kp(phi):
    """K'(φ) = 4·K_geo · φ³"""
    return 4.0 * K_geo * phi**3


# ============================================================================
# SOURCE (regularized point mass — Gaussian)
# ============================================================================

def rho_source(r):
    """Gaussian-regularized point source.  ∫ρ dV = M_source."""
    norm = M_source / ((2.0*np.pi*a_source**2)**1.5)
    return norm * np.exp(-0.5 * (r/a_source)**2)


# ============================================================================
# PoC.1 — LINEAR YUKAWA BENCHMARK
# ============================================================================
# Linearized EOM around Φ_0=1:  (∇² - μ²)δφ = -(q/Φ_0)ρ/K_geo
# μ² = β/K_geo,  μ = √(β/K_geo)  →  Yukawa mass / inverse Compton range
#
# Green's fn: (∇² - μ²)G = -δ³ ⇒ G(r) = exp(-μr)/(4πr)
# So δφ(r) = (q/(Φ_0·K_geo)) · ∫ G(|r-r'|) ρ(r') d³r'
# For Gaussian source, this convolution has known closed form (smeared Yukawa)

def yukawa_linear(r, mu):
    """
    Smeared Yukawa response to Gaussian source
    (closed form: convolution of G with Gaussian).
    For sharp δ-source: δφ(r) = (1/4πK)·(q/Φ_0)·M·exp(-μr)/r
    For Gaussian source of width a: smeared variant, but for r >> a → same form.
    """
    sharp = (q_over_Phi0 / K_geo) * M_source * np.exp(-mu * r) / (4.0 * np.pi * r)
    # Note: for r ≪ a_source the sharp formula diverges; smearing regularizes.
    # For PoC PoC.1 we use sharp form since we verify in r >> a regime.
    return sharp

mu_yukawa = np.sqrt(beta / K_geo)  # = 1.0 in PoC units


# ============================================================================
# PoC.2 — NONLINEAR EOM SOLVER
# ============================================================================
# Static, spherically symmetric, in flat-spatial (M9.1''-effective sign):
#
#   K(φ)·[φ'' + (2/r)φ'] + ½K'(φ)·(φ')² + V'(φ) + (q/Φ_0)ρ(r) = 0
#
# Solving for φ''(r):
#   φ'' = -[(2/r) K(φ) φ' + ½ K'(φ) (φ')² + V'(φ) + (q/Φ_0)ρ(r)] / K(φ)
#
# Strategy: shoot from r_max (where φ ≈ Φ_0 + Yukawa tail) inward to r_min.
# BC at r_max: φ(r_max) = Φ_0 + δφ_lin(r_max),
#              φ'(r_max) = -μ·δφ_lin(r_max) (matching exponential decay).

def eom_static(r, y):
    phi, phip = y
    if r < 1e-9:
        return [phip, 0.0]
    K_v = K(phi)
    if abs(K_v) < 1e-12:
        return [phip, 0.0]
    Kp_v = Kp(phi)
    Vp_v = Vp(phi)
    src = q_over_Phi0 * rho_source(r)

    phi_dd = -((2.0/r) * K_v * phip + 0.5 * Kp_v * phip**2 + Vp_v + src) / K_v
    return [phip, phi_dd]


def solve_phi_sol():
    """Shoot inward from r=r_max. Returns r_grid, phi_grid, phip_grid."""
    delta_max = yukawa_linear(r_max, mu_yukawa)  # tail value
    phi_init = Phi_0 + delta_max
    phip_init = -mu_yukawa * delta_max  # d/dr of Yukawa exp decay ≈ -μ·δφ for r ≫ 1/μ

    sol = solve_ivp(
        eom_static,
        t_span=[r_max, r_min],
        y0=[phi_init, phip_init],
        method='RK45',
        rtol=1e-9,
        atol=1e-11,
        dense_output=True,
        max_step=0.05,
    )
    if not sol.success:
        raise RuntimeError(f"ODE failed: {sol.message}")

    r_grid = np.linspace(r_min, r_max, 600)
    y_grid = sol.sol(r_grid)
    return r_grid, y_grid[0], y_grid[1]


# ============================================================================
# PoC.4 — LINEARIZATION SPECTRUM (radial s-wave around Φ_sol)
# ============================================================================
# Static fluctuation eq for s-wave u(r), δΦ(r,t) = u(r)·e^{-iωt}/r:
# (chosen so kinetic boundary terms vanish at r=0 with u(0)=0)
#
# Working from the action expanded to O(δΦ²):
#   S_quad = ½ ∫ [K(Φ_sol)·(δΦ̇)² - K(Φ_sol)·(∂_r δΦ)² - M²_eff(r)·(δΦ)²]·4πr² dr dt
#
# where M²_eff(r) = V''(Φ_sol(r)) - sign convention as in M9.3.1
#                 ≡ -Vpp(Φ_sol(r))   [so that around Φ_sol→Φ_0=1: M²_eff = +β]
#
# Eigenvalue equation:
#   -d/dr [r²·K(Φ_sol)·du/dr] + r²·M²_eff(r)·u = ω²·r²·K(Φ_sol)·u
#
# Sign choice: M²_eff(r) = -V''(Φ_sol(r)) is the M9.1''-hyperbolic
# effective curvature, consistent with M9.3.1 Yukawa stability (M_eff²=+β
# at Φ_sol → Φ_0 = 1).
#
# Discretization on radial grid → generalized eigenvalue problem A u = ω² B u
# where A = -d/dr[r²K d/dr] + r²M²_eff (Hermitian under boundary conds)
#       B = r²·K(Φ_sol) (positive-definite weight).

def build_spectrum(r_grid, phi_grid):
    N = len(r_grid)
    dr = np.gradient(r_grid)

    K_at = K(phi_grid)
    Vpp_at = Vpp(phi_grid)
    M2_eff = -Vpp_at   # M9.1'' hyperbolic-effective sign (gives +β at Φ_0)

    r2K = r_grid**2 * K_at

    # Build finite-difference operator
    # Approximate -d/dr [r²K · du/dr] by 3-point Laplacian
    # plus M²_eff term (diagonal)
    A = np.zeros((N, N))
    B = np.zeros((N, N))

    for i in range(1, N-1):
        h_minus = r_grid[i] - r_grid[i-1]
        h_plus = r_grid[i+1] - r_grid[i]
        h_avg = 0.5*(h_minus + h_plus)

        # r²K at midpoints
        rK_left = 0.5*(r2K[i] + r2K[i-1]) / h_minus
        rK_right = 0.5*(r2K[i+1] + r2K[i]) / h_plus

        # -d/dr[r²K du/dr] at point i:
        A[i, i-1] = -rK_left / h_avg
        A[i, i+1] = -rK_right / h_avg
        A[i, i]   = (rK_left + rK_right) / h_avg
        # + r²M²_eff·u term
        A[i, i] += r_grid[i]**2 * M2_eff[i]

        # Weight on RHS:
        B[i, i] = r_grid[i]**2 * K_at[i]

    # Boundary conds: u(0)=0 (regular at origin), u(r_max)=0 (decay)
    # Already enforced by NOT updating row 0 and row N-1 (Dirichlet, B=0 there)
    A[0, 0] = 1.0
    A[-1, -1] = 1.0
    B[0, 0] = 0.0
    B[-1, -1] = 0.0

    # Solve generalized eigenvalue problem on interior
    A_int = A[1:-1, 1:-1]
    B_int = B[1:-1, 1:-1]

    # Symmetrize numerically
    A_int = 0.5*(A_int + A_int.T)
    B_int = 0.5*(B_int + B_int.T)

    # Solve A·u = ω²·B·u
    try:
        eigenvalues, eigenvectors = eigh(A_int, B_int)
    except np.linalg.LinAlgError:
        eigenvalues = np.array([np.nan])
        eigenvectors = None

    return eigenvalues, eigenvectors


# ============================================================================
# RUN ALL TESTS
# ============================================================================

def run_tests():
    print("=" * 78)
    print("M11.S — Mini-PoC: single soliton Φ_sol(r) feasibility")
    print("=" * 78)
    print(f"Units: β=K_geo=Φ_0=q/Φ_0=1; M_source={M_source}, a_source={a_source}")
    print(f"Yukawa mass μ = √(β/K_geo) = {mu_yukawa:.4f}  (λ_C = 1.0)")
    print(f"Grid: r ∈ [{r_min}, {r_max}], 600 points")
    print()

    results = {}

    # --- PoC.1 ---
    print("-" * 78)
    print("PoC.1 — Linear Yukawa benchmark")
    print("-" * 78)
    r_test = np.array([0.5, 1.0, 2.0, 5.0, 10.0])
    yuk_vals = yukawa_linear(r_test, mu_yukawa)
    print(f"Linear δφ_Yuk(r) at r ∈ {list(r_test)}:")
    for r, y in zip(r_test, yuk_vals):
        print(f"   r={r:5.2f}:  δφ = {y:+.6e}   (Φ_0 + δφ = {1.0+y:.6f})")
    decay_check = abs(yuk_vals[1] / yuk_vals[0] - np.exp(-mu_yukawa*0.5) * 0.5/1.0) < 0.05
    print(f"  Yukawa decay rate verified (e^-μ·Δr ratio): {decay_check}")
    results['PoC.1'] = decay_check

    # --- PoC.2 ---
    print()
    print("-" * 78)
    print("PoC.2 — Nonlinear ODE solver Φ_sol(r)")
    print("-" * 78)
    try:
        r_grid, phi_grid, phip_grid = solve_phi_sol()
        # Verify: phi_grid(r_max) ≈ Φ_0
        phi_at_rmax = phi_grid[-1]
        phi_at_rmin = phi_grid[0]
        phi_minus_1 = phi_grid - Phi_0
        nontrivial = np.max(np.abs(phi_minus_1)) > 1e-4
        regular_origin = np.isfinite(phi_at_rmin)
        asymptotic_ok = abs(phi_at_rmax - Phi_0) < 1e-3

        print(f"  Φ_sol(r_min={r_min}) = {phi_at_rmin:.6f}")
        print(f"  Φ_sol(r=1.0)         = {phi_grid[np.argmin(np.abs(r_grid-1.0))]:.6f}")
        print(f"  Φ_sol(r=2.0)         = {phi_grid[np.argmin(np.abs(r_grid-2.0))]:.6f}")
        print(f"  Φ_sol(r_max={r_max}) = {phi_at_rmax:.6f}")
        print(f"  max|Φ_sol - Φ_0|     = {np.max(np.abs(phi_minus_1)):.4e}")
        print(f"  finite at origin: {regular_origin}")
        print(f"  asymptotic Φ_0:   {asymptotic_ok}")
        print(f"  nontrivial:       {nontrivial}")
        h1_pass = regular_origin and asymptotic_ok and nontrivial
        results['PoC.2'] = h1_pass
        print(f"  H1 (existence) + H2 (asymptotic) verdict: {'PASS' if h1_pass else 'FAIL'}")
    except Exception as e:
        print(f"  Solver failed: {e}")
        r_grid, phi_grid = None, None
        results['PoC.2'] = False

    # --- PoC.3 ---
    if phi_grid is not None:
        print()
        print("-" * 78)
        print("PoC.3 — Linear vs nonlinear comparison")
        print("-" * 78)
        delta_lin = yukawa_linear(r_grid, mu_yukawa)
        delta_nonlin = phi_grid - Phi_0
        # Restrict to r > 3·a_source to avoid source core
        mask = r_grid > 3.0 * a_source
        rel_diff = np.abs(delta_nonlin[mask] - delta_lin[mask]) / (np.abs(delta_lin[mask]) + 1e-12)
        max_rel = np.max(rel_diff)
        print(f"  max relative diff (lin vs nonlin) for r > {3*a_source}: {max_rel*100:.2f}%")
        print(f"  At r=2.0:   lin={delta_lin[np.argmin(np.abs(r_grid-2.0))]:+.4e},"
              f"  nonlin={delta_nonlin[np.argmin(np.abs(r_grid-2.0))]:+.4e}")
        print(f"  At r=5.0:   lin={delta_lin[np.argmin(np.abs(r_grid-5.0))]:+.4e},"
              f"  nonlin={delta_nonlin[np.argmin(np.abs(r_grid-5.0))]:+.4e}")
        # Nonlinear should agree with linear in tail (large r)
        tail_mask = r_grid > 4.0
        tail_diff = np.max(np.abs(delta_nonlin[tail_mask] - delta_lin[tail_mask]) / np.abs(delta_lin[tail_mask] + 1e-12))
        print(f"  tail (r > 4 λ_C) max relative diff: {tail_diff*100:.2f}%")
        tail_ok = tail_diff < 0.10  # within 10% in tail
        results['PoC.3'] = tail_ok
        print(f"  Tail consistency verdict: {'PASS' if tail_ok else 'FAIL'}")

    # --- PoC.4 ---
    if phi_grid is not None:
        print()
        print("-" * 78)
        print("PoC.4 — Linearization spectrum (radial s-wave)")
        print("-" * 78)
        try:
            eigenvalues, _ = build_spectrum(r_grid, phi_grid)
            eigenvalues = eigenvalues[np.isfinite(eigenvalues)]
            n_neg = np.sum(eigenvalues < -1e-6)
            n_zero = np.sum(np.abs(eigenvalues) < 1e-6)
            print(f"  Number of negative ω² modes:  {n_neg}")
            print(f"  Number of near-zero modes:    {n_zero}")
            print(f"  Lowest 5 ω² values:           {sorted(eigenvalues)[:5]}")
            print(f"  Asymptotic ω² (continuum):    expect ω² → β + k² ≈ 1 + k²")
            print(f"  Highest ω² in grid:           {eigenvalues.max():.3f}")
            # H3: stable means n_neg = 0 in s-wave (zero mode is L=1 vector, not s-wave)
            h3_pass = n_neg == 0
            results['PoC.4'] = h3_pass
            print(f"  H3 (s-wave stability) verdict: {'PASS' if h3_pass else 'FAIL'}")
            if n_neg > 0:
                print(f"  WARNING: {n_neg} unstable mode(s) found in s-wave —")
                print(f"  may indicate sign convention issue or genuine instability.")
        except Exception as e:
            print(f"  Spectrum solver failed: {e}")
            results['PoC.4'] = False

    # --- PoC.5: verdict ---
    print()
    print("=" * 78)
    print("PoC.5 — VERDICT (M11.S feasibility for full audit)")
    print("=" * 78)
    all_pass = all(results.values())
    n_pass = sum(results.values())
    n_total = len(results)
    print(f"  Tests passed: {n_pass}/{n_total}")
    for test, status in results.items():
        symbol = "✓" if status else "✗"
        print(f"    {symbol} {test}: {'PASS' if status else 'FAIL'}")
    print()
    if all_pass:
        print("  ✅ M11.S FEASIBILITY CONFIRMED")
        print("  Φ_sol(r) classical solution exists, asymptotically matches Yukawa,")
        print("  s-wave linearization spectrum is non-negative.")
        print("  → READY for full M11.S audit (Branch I bottom-up start)")
    else:
        print("  ⚠️  M11.S FEASIBILITY: PARTIAL")
        print("  Some PoC tests failed — investigate before committing to full audit.")
        print("  Diagnostic info above; possible issues:")
        print("    (a) sign convention for spatial Laplacian vs M9.1'' hyperbolic")
        print("    (b) source regularization width a_source vs λ_C")
        print("    (c) numerical method (RK45 may need adaptive, or BVP solver)")
    print()

    return results, all_pass


if __name__ == "__main__":
    results, all_pass = run_tests()
    sys.exit(0 if all_pass else 1)
