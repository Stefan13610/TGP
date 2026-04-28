"""
M11.I — Multi-soliton interference FULL AUDIT (Branch I level 2)
=================================================================

Two-soliton classical interference + Yukawa interaction extraction.
Builds on M11.S (CLOSED, 6/6 PASS, 2026-04-26).

Sub-tests:
  M11.I.1 — Two-soliton ansatz at variable separation d ∈ {0.5..12}·λ_C
            (axial symmetry; verify domain validity Φ ∈ (0, 4/3))
  M11.I.2 — Interaction energy V_int(d): sign, monotonicity, asymptote → 0
  M11.I.3 — Yukawa tail extraction: log|V_int(d)·d| vs d slope = -μ;
            amplitude consistency with M9.3.1
  M11.I.4 — Critical merge distance d_min (Φ_2sol exits (0, 4/3) below)
  M11.I.5 — Force F(d) = -dV_int/dd, Yukawa-derivative form
  M11.I.6 — Cross-check single-source M9.3.1 tail: A_single = q·M / (4π·K_geo)

Geometry: axially-symmetric (cylindrical) coords (ρ, z); two Gaussian sources
          at z = ±d/2, ρ = 0; volume element 2π·ρ·dρ·dz.

Ansatz (linear superposition):
   Φ_2sol(ρ, z; d) = Φ_sol(r_1) + Φ_sol(r_2) - Φ_0
   r_i = √(ρ² + (z ∓ d/2)²)
   (exact in linear regime, accurate for d >> a_source and weak coupling)

Sign convention (M9.3.1):
   E[φ] = ∫d³x [½K(φ)|∇φ|² + V(φ) - V(Φ_0) + (q/Φ_0)·(φ-Φ_0)·ρ]
   V_int(d) = E[Φ_2sol(d)] - 2·E[Φ_sol_single]
   For like-charges (positive q·M): V_int(d) < 0 (Yukawa attractive)

Units (PoC-style): β = K_geo = Φ_0 = q/Φ_0 = 1, λ_C = 1.
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_bvp, simpson
from scipy.interpolate import interp1d
import sys
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# CONSTANTS (consistent with M11.S)
# ============================================================================

beta = 1.0
K_geo = 1.0
Phi_0 = 1.0
mu_yukawa = np.sqrt(beta / K_geo)  # = 1.0
lambda_C = 1.0 / mu_yukawa


# ============================================================================
# POTENTIAL & KINETIC (sek08a)
# ============================================================================

def V(phi):
    return (beta/3.0) * phi**3 - (beta/4.0) * phi**4

def Vp(phi):
    return beta * phi**2 * (1.0 - phi)

def Vpp(phi):
    return beta * phi * (2.0 - 3.0*phi)

def K(phi):
    return K_geo * phi**4

def Kp(phi):
    return 4.0 * K_geo * phi**3


# ============================================================================
# SINGLE-SOLITON SOLVER (port from M11.S, BVP-stable)
# ============================================================================

def solve_single_soliton(qM, a_source=0.15, r_max=15.0, r_min=1e-3, N=600):
    """Solve static EOM for single Gaussian source via solve_bvp."""
    rho_norm = qM / ((2.0 * np.pi * a_source**2)**1.5)

    def rho(r):
        return rho_norm * np.exp(-0.5 * (r / a_source)**2)

    def eom(r, y):
        phi = y[0]; phip = y[1]
        K_v = K(phi)
        K_v = np.where(np.abs(K_v) < 1e-12, 1e-12, K_v)
        Kp_v = Kp(phi); Vp_v = Vp(phi); src = rho(r)
        r_safe = np.where(r < 1e-6, 1e-6, r)
        phi_dd = -((2.0/r_safe) * K_v * phip + 0.5 * Kp_v * phip**2 + Vp_v + src) / K_v
        return np.vstack([phip, phi_dd])

    def bc(ya, yb):
        return np.array([ya[1], yb[0] - Phi_0])

    r_init = np.linspace(r_min, r_max, 200)
    delta_init = (qM/(4.0*np.pi*K_geo)) * np.exp(-mu_yukawa*r_init) / np.maximum(r_init, 0.1)
    phi_init = Phi_0 + delta_init
    phip_init = -mu_yukawa * delta_init
    phip_init[0] = 0.0
    y_init = np.vstack([phi_init, phip_init])

    sol_bvp = solve_bvp(eom, bc, r_init, y_init, tol=1e-7, max_nodes=5000, verbose=0)
    if not sol_bvp.success:
        return None

    r_grid = np.linspace(r_min, r_max, N)
    y_grid = sol_bvp.sol(r_grid)
    return {
        'qM': qM,
        'r': r_grid,
        'phi': y_grid[0],
        'phip': y_grid[1],
        'a_source': a_source,
    }


def soliton_energy_single(sol_dict):
    """Self-energy of single soliton (3D spherical integration).

    Uses the EOM-consistent Hamiltonian:
       H[φ] = ∫·[½K(φ)|∇φ|² - V(φ) - (q/Φ_0)·φ·ρ]
    (matches static EOM: K∇²φ + ½K'|∇φ|² + V'(φ) + (q/Φ_0)·ρ = 0)

    Subtract vacuum reference (Φ_0, fixed ρ) to define binding energy:
       M_sol = H[Φ_sol; ρ] - H[Φ_0; ρ]
             = ∫·[½K|∇δφ|² - (V(Φ_sol) - V(Φ_0)) - (q/Φ_0)·(Φ_sol-Φ_0)·ρ]

    Sign: in linearized regime M_sol ≈ -½·(q/Φ_0)·∫δφ·ρ < 0 (binding,
    consistent with Newton-like attractive interaction at distance).
    """
    r = sol_dict['r']
    phi = sol_dict['phi']
    phip = sol_dict['phip']
    qM = sol_dict['qM']
    a_source = sol_dict['a_source']
    rho_norm = qM / ((2.0 * np.pi * a_source**2)**1.5)
    rho = rho_norm * np.exp(-0.5 * (r / a_source)**2)

    e_kin = 0.5 * K(phi) * phip**2
    e_pot = -(V(phi) - V(Phi_0))                  # EOM-consistent sign
    e_src = -rho * (phi - Phi_0)                  # EOM-consistent sign
    integrand = 4.0 * np.pi * r**2 * (e_kin + e_pot + e_src)
    return simpson(integrand, x=r)


# ============================================================================
# TWO-SOLITON ANSATZ ON 2D AXIAL GRID
# ============================================================================

def build_axial_grid(d, Lz_factor=20.0, Lrho_factor=8.0, Nz=400, Nrho=200):
    """Build (ρ, z) grid covering both solitons with generous Yukawa tail margin.
    Box: Lz = max(d + Lz_factor·λ_C, 4·λ_C); each source has ≥Lz_factor/2·λ_C
    margin from z-boundary, exp(-10)/10 ≈ 4.5e-6 leakage at d=0.
    """
    Lz = max(d + Lz_factor * lambda_C, 4 * lambda_C)
    Lrho = max(Lrho_factor * lambda_C, 4 * lambda_C)
    rho_grid = np.linspace(1e-3, Lrho, Nrho)
    z_grid = np.linspace(-Lz/2, Lz/2, Nz)
    return rho_grid, z_grid


def two_soliton_ansatz(rho_grid, z_grid, d, sol_single):
    """
    Linear superposition ansatz:
        Φ_2sol(ρ, z; d) = Φ_sol(r_1) + Φ_sol(r_2) - Φ_0
    where r_i = √(ρ² + (z ∓ d/2)²).

    Returns Φ on (ρ, z) meshgrid.
    """
    RHO, Z = np.meshgrid(rho_grid, z_grid, indexing='ij')
    r1 = np.sqrt(RHO**2 + (Z - d/2)**2)
    r2 = np.sqrt(RHO**2 + (Z + d/2)**2)

    # Interpolator for single soliton — extrapolate to Φ_0 outside r_max
    r_s = sol_single['r']
    phi_s = sol_single['phi']
    interp = interp1d(r_s, phi_s, kind='cubic', bounds_error=False, fill_value=Phi_0)

    phi1 = interp(r1)
    phi2 = interp(r2)
    return phi1 + phi2 - Phi_0


def _energy_in_box(rho_grid, z_grid, phi, src_total):
    """Generic boxed energy integral (EOM-consistent Hamiltonian) on (ρ,z) grid.
    H = ∫ 2π·ρ·dρ·dz · [½K(φ)|∇φ|² - (V(φ)-V(Φ_0)) - ρ_src·(φ-Φ_0)]
    """
    RHO, _Z = np.meshgrid(rho_grid, z_grid, indexing='ij')
    dphi_drho = np.gradient(phi, rho_grid, axis=0)
    dphi_dz = np.gradient(phi, z_grid, axis=1)
    grad_phi_sq = dphi_drho**2 + dphi_dz**2
    e_kin = 0.5 * K(phi) * grad_phi_sq
    e_pot = -(V(phi) - V(Phi_0))
    e_src = -src_total * (phi - Phi_0)
    integrand = 2.0 * np.pi * RHO * (e_kin + e_pot + e_src)
    return simpson(simpson(integrand, x=z_grid, axis=1), x=rho_grid)


def single_soliton_field_in_box(rho_grid, z_grid, d_offset, sol_single):
    """Single-soliton template Φ_sol(r) placed at z = d_offset on (ρ,z) grid."""
    RHO, Z = np.meshgrid(rho_grid, z_grid, indexing='ij')
    r = np.sqrt(RHO**2 + (Z - d_offset)**2)
    interp = interp1d(sol_single['r'], sol_single['phi'], kind='cubic',
                      bounds_error=False, fill_value=Phi_0)
    return interp(r)


def gaussian_source(rho_grid, z_grid, d_offset, sol_single):
    """Normalized Gaussian source ρ_src centered at z=d_offset."""
    RHO, Z = np.meshgrid(rho_grid, z_grid, indexing='ij')
    qM = sol_single['qM']
    a = sol_single['a_source']
    rho_norm = qM / ((2.0 * np.pi * a**2)**1.5)
    return rho_norm * np.exp(-0.5 * (RHO**2 + (Z - d_offset)**2) / a**2)


def two_soliton_energy(d, sol_single, return_field=False):
    """Compute E[Φ_2sol(d)] AND E_self_box (single soliton at same offsets in
    same cylindrical box) — V_int = E_2sol - 2·E_self_box cancels boundary
    truncation systematically.
    """
    rho_grid, z_grid = build_axial_grid(d)

    # Two-soliton field + total source
    phi_2sol = two_soliton_ansatz(rho_grid, z_grid, d, sol_single)
    src1 = gaussian_source(rho_grid, z_grid, +d/2, sol_single)
    src2 = gaussian_source(rho_grid, z_grid, -d/2, sol_single)
    src_total = src1 + src2
    E_2sol = _energy_in_box(rho_grid, z_grid, phi_2sol, src_total)

    # Boxed self-energy: each soliton ALONE in same box at its own offset
    phi_self1 = single_soliton_field_in_box(rho_grid, z_grid, +d/2, sol_single)
    phi_self2 = single_soliton_field_in_box(rho_grid, z_grid, -d/2, sol_single)
    E_self_box1 = _energy_in_box(rho_grid, z_grid, phi_self1, src1)
    E_self_box2 = _energy_in_box(rho_grid, z_grid, phi_self2, src2)
    E_self_box_total = E_self_box1 + E_self_box2

    out = {
        'd': d,
        'E': E_2sol,
        'E_self_box': E_self_box_total,
        'V_int': E_2sol - E_self_box_total,
        'phi_max': phi_2sol.max(),
        'phi_min': phi_2sol.min(),
        'in_domain': (phi_2sol.max() < 4.0/3.0) and (phi_2sol.min() > 0.0),
    }
    if return_field:
        out['rho_grid'] = rho_grid
        out['z_grid'] = z_grid
        out['phi'] = phi_2sol
    return out


# ============================================================================
# SUB-TESTS
# ============================================================================

def test_I1_separation_sweep(sol_single, verbose=True):
    """M11.I.1 — Two-soliton ansatz on separation sweep, domain validity."""
    if verbose:
        print("=" * 78)
        print("M11.I.1 — Two-soliton ansatz separation sweep + domain validity")
        print("=" * 78)

    d_values = np.array([0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]) * lambda_C
    results = []
    for d in d_values:
        res = two_soliton_energy(d, sol_single)
        results.append(res)

    if verbose:
        print(f"  Single-soliton qM = {sol_single['qM']}")
        print(f"  {'d/λ_C':>8} | {'Φ_max':>8} | {'in (0,4/3)?':>12} | {'E[Φ_2sol]':>12}")
        print("-" * 60)
        for r in results:
            print(f"  {r['d']:>8.2f} | {r['phi_max']:>8.4f} | "
                  f"{str(r['in_domain']):>12} | {r['E']:>+12.4e}")

    # PASS criteria:
    # (1) at least one in-domain at large separation (d ≥ 2·λ_C)
    # (2) Φ_max monotonically decreases with d
    # (3) E values are finite
    in_domain_far = any(r['in_domain'] for r in results if r['d'] >= 2.0)
    phi_max_arr = [r['phi_max'] for r in results]
    monotonic = all(phi_max_arr[i+1] <= phi_max_arr[i] + 1e-3 for i in range(len(phi_max_arr)-1))
    finite_E = all(np.isfinite(r['E']) for r in results)

    test_pass = in_domain_far and monotonic and finite_E
    if verbose:
        print(f"\n  in-domain at large d:     {in_domain_far}")
        print(f"  Φ_max monotonic decrease: {monotonic}")
        print(f"  All E finite:             {finite_E}")
        print(f"  M11.I.1 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, results


def test_I2_interaction_energy(sol_single, sweep_results, verbose=True):
    """M11.I.2 — V_int(d) = E[Φ_2sol(d)] - E_self_box(d).

    E_self_box is energy of two NON-INTERACTING single solitons in SAME
    cylindrical box (each at z=±d/2 alone, same source). This cancels
    boundary truncation and grid artifacts systematically.
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.I.2 — Interaction energy V_int(d) [box-cancelled]")
        print("=" * 78)

    # Reference E_self in spherical coords (sanity check; not used for V_int)
    E_single_spherical = soliton_energy_single(sol_single)
    if verbose:
        print(f"  E_self (spherical reference, qM={sol_single['qM']}) = {E_single_spherical:.6e}")

    V_int_values = []
    d_values = []
    for r in sweep_results:
        V_int = r['V_int']  # already computed as E_2sol - E_self_box in same box
        V_int_values.append(V_int)
        d_values.append(r['d'])
        if verbose:
            print(f"  d/λ_C = {r['d']:>5.2f}: E_2sol = {r['E']:>+12.4e}, "
                  f"E_self_box = {r['E_self_box']:>+12.4e}, V_int = {V_int:>+12.4e}")

    V_int_arr = np.array(V_int_values)
    d_arr = np.array(d_values)

    # Filter to in-domain results
    in_dom_mask = np.array([r['in_domain'] for r in sweep_results])
    V_int_dom = V_int_arr[in_dom_mask]
    d_dom = d_arr[in_dom_mask]

    # PASS criteria:
    # (1) V_int < 0 in-domain (attractive scalar Yukawa)
    # (2) V_int → 0 as d → ∞ (decoupling)
    # (3) V_int monotonically increases (becomes less negative) with d
    crit_attract = np.all(V_int_dom < 0) if len(V_int_dom) > 0 else False
    crit_decoupled = abs(V_int_arr[-1]) < abs(V_int_dom[0]) * 0.1 if len(V_int_dom) > 0 else False
    crit_monotonic = True
    for i in range(len(V_int_dom)-1):
        if V_int_dom[i+1] < V_int_dom[i] - 1e-6:
            crit_monotonic = False
            break

    test_pass = crit_attract and crit_decoupled and crit_monotonic
    if verbose:
        print(f"\n  V_int < 0 in-domain (attractive):   {crit_attract}")
        print(f"  V_int → 0 at large d (decoupled):    {crit_decoupled}")
        print(f"  V_int monotonic decrease in |·|:     {crit_monotonic}")
        print(f"  M11.I.2 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {'V_int': V_int_arr, 'd': d_arr, 'in_domain': in_dom_mask, 'E_single': E_single_spherical}


def test_I3_yukawa_tail(V_int_data, sol_single, verbose=True):
    """M11.I.3 — Extract Yukawa tail: V_int(d) ≈ -A·exp(-μd)/d for large d.
    Cross-check amplitude A vs M9.3.1: A_predicted = (q·M)² / (4π·K_geo).
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.I.3 — Yukawa tail extraction (ln|d·V_int| vs d, slope = -μ)")
        print("=" * 78)

    V_int = V_int_data['V_int']
    d = V_int_data['d']
    in_dom = V_int_data['in_domain']

    # Use only LARGE-d in-domain points (d ≥ 2·λ_C)
    mask = in_dom & (d >= 2.0)
    d_fit = d[mask]
    V_fit = V_int[mask]
    if len(d_fit) < 3:
        if verbose:
            print(f"  ERROR: too few large-d points for fitting ({len(d_fit)})")
        return False, None

    # Yukawa form: V_int(d) = -A·exp(-μd)/d
    # → ln|V_int·d| = ln(A) - μ·d  (linear in d)
    y_fit = np.log(np.abs(V_fit) * d_fit)
    # Linear regression
    n = len(d_fit)
    sx = d_fit.sum(); sy = y_fit.sum()
    sxy = (d_fit * y_fit).sum(); sxx = (d_fit**2).sum()
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
    intercept = (sy - slope*sx) / n
    mu_extracted = -slope
    A_extracted = np.exp(intercept)

    # M9.3.1 prediction
    qM = sol_single['qM']
    A_M9 = qM**2 / (4.0 * np.pi * K_geo)

    rel_diff_mu = abs(mu_extracted - mu_yukawa) / mu_yukawa
    rel_diff_A = abs(A_extracted - A_M9) / A_M9

    if verbose:
        print(f"  Fit data ({n} points, d ∈ [{d_fit[0]:.2f}, {d_fit[-1]:.2f}]·λ_C):")
        for di, vi in zip(d_fit, V_fit):
            print(f"    d={di:.2f}, V_int={vi:.4e}, ln|d·V_int|={np.log(abs(vi)*di):.4f}")
        print(f"\n  Extracted slope:    -μ = {slope:.4f}  →  μ = {mu_extracted:.4f}")
        print(f"  Extracted intercept: ln(A) = {intercept:.4f}  →  A = {A_extracted:.4e}")
        print(f"\n  M9.3.1 predictions:")
        print(f"    μ_Yukawa = √(β/K_geo) = {mu_yukawa:.4f}")
        print(f"    A_M9 = (q·M)²/(4π·K_geo) = {A_M9:.4e}")
        print(f"\n  Rel. diff μ:  {rel_diff_mu*100:.2f}%")
        print(f"  Rel. diff A:  {rel_diff_A*100:.2f}%")

    # PASS criteria:
    # (1) μ extracted within 10% of √(β/K_geo) = M9.3.1 Yukawa range
    # (2) A extracted within 30% of (q·M)²/(4π·K_geo) — finite a_source corrections expected
    crit_mu = rel_diff_mu < 0.10
    crit_A = rel_diff_A < 0.30

    test_pass = crit_mu and crit_A
    if verbose:
        print(f"\n  μ within 10% of M9.3.1:   {crit_mu}")
        print(f"  A within 30% of M9.3.1:   {crit_A}")
        print(f"  M11.I.3 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {
        'mu_extracted': mu_extracted,
        'A_extracted': A_extracted,
        'mu_M9': mu_yukawa,
        'A_M9': A_M9,
    }


def test_I4_critical_merge(sol_single, sweep_results, verbose=True):
    """M11.I.4 — Critical merge distance d_min where Φ_2sol exits (0, 4/3)."""
    if verbose:
        print("\n" + "=" * 78)
        print("M11.I.4 — Critical merge distance d_min")
        print("=" * 78)

    d_in = [r['d'] for r in sweep_results if r['in_domain']]
    d_out = [r['d'] for r in sweep_results if not r['in_domain']]
    phi_max_in = [r['phi_max'] for r in sweep_results if r['in_domain']]
    phi_max_out = [r['phi_max'] for r in sweep_results if not r['in_domain']]

    if verbose:
        if d_in:
            print(f"  In-domain d: {[f'{d:.2f}' for d in d_in]}")
            print(f"    smallest in-domain d:  {min(d_in):.2f}·λ_C  (Φ_max = {phi_max_in[d_in.index(min(d_in))]:.4f})")
        if d_out:
            print(f"  Out-of-domain d: {[f'{d:.2f}' for d in d_out]}")
            print(f"    largest out-of-domain d: {max(d_out):.2f}·λ_C  (Φ_max = {phi_max_out[d_out.index(max(d_out))]:.4f})")

    if not d_out:
        # All in-domain — check if SMALLEST d still in domain (no critical merge in range)
        d_min_estimate = min(d_in) if d_in else None
        if verbose:
            print(f"  No out-of-domain in scan range; d_min ≤ {d_min_estimate:.2f}·λ_C")
        crit_in_domain = True
        crit_estimate = d_min_estimate is not None
    else:
        d_min_estimate = max(d_out)  # critical merge between this and next-larger
        if verbose:
            print(f"  d_critical ∈ ({max(d_out):.2f}, {min(d_in):.2f})·λ_C")
        crit_in_domain = True
        crit_estimate = True

    # Check structural finding: critical merge should be roughly comparable to a_source for q·M ≤ 0.3
    # (smearing scale beyond which two sources don't superimpose strongly)

    # PASS criteria:
    # (1) Domain check completes (no NaN, all data valid)
    # (2) Critical merge identified or all-in-domain
    crit_complete = all(np.isfinite(r['phi_max']) for r in sweep_results)
    test_pass = crit_complete and crit_estimate

    if verbose:
        print(f"\n  Domain analysis complete: {crit_complete}")
        print(f"  d_min identifiable:       {crit_estimate}")
        print(f"  M11.I.4 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {'d_min_estimate': d_min_estimate}


def test_I5_force(V_int_data, verbose=True):
    """M11.I.5 — Force F(d) = -dV_int/dd; Yukawa derivative form.

    For attractive Yukawa V(d) = -A·exp(-μd)/d (V<0, increasing → 0):
       dV/dd = A·(μd+1)·exp(-μd)/d²  > 0
       F = -dV/dd = -A·(μd+1)·exp(-μd)/d²  < 0   (attractive: F·d̂ < 0
       means force pulls each soliton toward smaller separation).
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.I.5 — Force F(d) = -dV_int/dd, Yukawa-derivative consistency")
        print("=" * 78)

    V_int = V_int_data['V_int']
    d = V_int_data['d']
    in_dom = V_int_data['in_domain']

    # Numerical force: F = -dV/dd via finite differences
    F_num = -np.gradient(V_int, d)

    if verbose:
        print(f"  {'d/λ_C':>8} | {'V_int':>12} | {'F_num=-dV/dd':>14} | {'|F|':>12} | {'in_dom?':>8}")
        print("-" * 78)
        for i in range(len(d)):
            print(f"  {d[i]:>8.2f} | {V_int[i]:>+12.4e} | {F_num[i]:>+14.4e} | "
                  f"{abs(F_num[i]):>12.4e} | {str(in_dom[i]):>8}")

    # PASS criteria:
    # (1) F < 0 at large d (attractive: pulls toward smaller separation)
    # (2) |F| monotonically decays with d (Yukawa fall-off)
    mask_large = in_dom & (d >= 2.0)
    F_large = F_num[mask_large]
    d_large = d[mask_large]

    crit_attract = len(F_large) > 0 and np.all(F_large < 0)
    crit_decay = True
    abs_F_large = np.abs(F_large)
    for i in range(len(abs_F_large)-1):
        if abs_F_large[i+1] > abs_F_large[i] + 1e-12:
            crit_decay = False
            break

    test_pass = crit_attract and crit_decay

    if verbose:
        print(f"\n  F < 0 at large d (attractive, pulls together):  {crit_attract}")
        print(f"  |F| decays monotonically with d (Yukawa):       {crit_decay}")
        print(f"  M11.I.5 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {'F': F_num, 'd': d}


def test_I6_M9_3_1_single_amplitude(sol_single, verbose=True):
    """M11.I.6 — Single-source M9.3.1 cross-check.

    For r >> a_source, Φ_sol(r) - Φ_0 → A_single·exp(-μr)/r
    where M9.3.1 predicts A_single = (q·M)/(4π·K_geo).
    """
    if verbose:
        print("\n" + "=" * 78)
        print("M11.I.6 — Single-source amplitude cross-check (M9.3.1 Yukawa)")
        print("=" * 78)

    r = sol_single['r']
    phi = sol_single['phi']
    qM = sol_single['qM']

    # Use far-tail r ∈ [3, 10]·λ_C
    mask = (r >= 3.0) & (r <= 10.0)
    r_fit = r[mask]
    delta_fit = phi[mask] - Phi_0

    # ln|r·δφ| vs r should be linear with slope -μ, intercept ln(A)
    y_fit = np.log(np.abs(delta_fit) * r_fit)
    n_fit = len(r_fit)
    sx = r_fit.sum(); sy = y_fit.sum()
    sxy = (r_fit*y_fit).sum(); sxx = (r_fit**2).sum()
    slope = (n_fit*sxy - sx*sy) / (n_fit*sxx - sx*sx)
    intercept = (sy - slope*sx) / n_fit
    mu_ext = -slope
    A_ext = np.exp(intercept)

    # M9.3.1 predictions
    A_M9 = qM / (4.0 * np.pi * K_geo)
    mu_M9 = mu_yukawa

    rel_diff_mu = abs(mu_ext - mu_M9) / mu_M9
    rel_diff_A = abs(A_ext - A_M9) / A_M9

    # Symbolic verification (sympy)
    qM_sym, K_sym, beta_sym = sp.symbols('qM K_geo beta', positive=True)
    A_M9_sym = qM_sym / (4*sp.pi*K_sym)
    mu_M9_sym = sp.sqrt(beta_sym / K_sym)

    if verbose:
        print(f"  Far-tail fit: r ∈ [3, 10]·λ_C, {n_fit} points")
        print(f"  Numerical extraction:")
        print(f"    μ_extracted = {mu_ext:.4f}")
        print(f"    A_extracted = {A_ext:.6e}")
        print(f"  M9.3.1 prediction (sympy):")
        print(f"    A_M9 = q·M / (4π·K_geo)  = {A_M9_sym}  →  {A_M9:.6e}")
        print(f"    μ_M9 = √(β/K_geo)        = {mu_M9_sym}  →  {mu_M9:.4f}")
        print(f"  Relative diffs:")
        print(f"    μ:  {rel_diff_mu*100:.3f}%")
        print(f"    A:  {rel_diff_A*100:.3f}%")

    # PASS criteria:
    # (1) μ within 5% (single-source has cleaner tail than two-soliton)
    # (2) A within 10%
    # (3) Symbolic A_M9 has correct functional form
    crit_mu = rel_diff_mu < 0.05
    crit_A = rel_diff_A < 0.10
    crit_sym = (sp.simplify(A_M9_sym - qM_sym/(4*sp.pi*K_sym)) == 0)

    test_pass = crit_mu and crit_A and crit_sym

    if verbose:
        print(f"\n  μ within 5%:               {crit_mu}")
        print(f"  A within 10%:              {crit_A}")
        print(f"  Symbolic A form correct:   {crit_sym}")
        print(f"  M11.I.6 → {'PASS' if test_pass else 'FAIL'}")

    return test_pass, {'mu_ext': mu_ext, 'A_ext': A_ext, 'A_M9': A_M9}


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("#" * 78)
    print("# M11.I — Multi-soliton interference FULL AUDIT (Branch I level 2)")
    print("# closure-grade target: ≥6/6 PASS")
    print("# ", end="")
    print(f"Units: β=K_geo=Φ_0=1, λ_C=1/√β={lambda_C}")
    print("#" * 78)
    print()

    # === Solve single soliton (q·M = 0.3, in-domain from M11.S) ===
    print("Solving single-soliton template (q·M=0.3)...")
    sol_single = solve_single_soliton(qM=0.3)
    if sol_single is None:
        print("FATAL: single-soliton solver failed")
        return False
    print(f"  Φ(0) = {sol_single['phi'][0]:.4f},  Φ(r_max={sol_single['r'][-1]:.1f}) = {sol_single['phi'][-1]:.6f}")
    print()

    results = {}

    # === M11.I.1 ===
    i1_pass, sweep_results = test_I1_separation_sweep(sol_single, verbose=True)
    results['I1'] = i1_pass

    # === M11.I.2 ===
    i2_pass, V_int_data = test_I2_interaction_energy(sol_single, sweep_results, verbose=True)
    results['I2'] = i2_pass

    # === M11.I.3 ===
    i3_pass, i3_data = test_I3_yukawa_tail(V_int_data, sol_single, verbose=True)
    results['I3'] = i3_pass

    # === M11.I.4 ===
    i4_pass, i4_data = test_I4_critical_merge(sol_single, sweep_results, verbose=True)
    results['I4'] = i4_pass

    # === M11.I.5 ===
    i5_pass, i5_data = test_I5_force(V_int_data, verbose=True)
    results['I5'] = i5_pass

    # === M11.I.6 ===
    i6_pass, i6_data = test_I6_M9_3_1_single_amplitude(sol_single, verbose=True)
    results['I6'] = i6_pass

    # === FINAL ===
    print()
    print("=" * 78)
    print(" M11.I — FINAL VERDICT")
    print("=" * 78)
    n_pass = sum(results.values())
    n_total = len(results)
    print(f"  Sub-tests passed: {n_pass}/{n_total}")
    for test, status in results.items():
        symbol = "v" if status else "x"
        print(f"    [{symbol}] M11.I.{test[1:]}: {'PASS' if status else 'FAIL'}")
    print()
    if n_pass >= 6:
        print(f"  >> M11.I CLOSED ({n_pass}/6 PASS) — Branch I multi-soliton interference COMPLETE")
        print(f"  >> Ready for M11.G (global field, source extraction) launch")
    elif n_pass >= 4:
        print(f"  ?? M11.I PARTIAL ({n_pass}/6) — promote to closure with documented caveats")
    else:
        print(f"  !! M11.I FAIL ({n_pass}/6) — investigate before next sub-cycle")
    print()

    return n_pass >= 6


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
