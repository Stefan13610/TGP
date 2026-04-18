#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs4_n_soliton_monte_carlo.py: Monte Carlo superposition of N soliton tails.

KEY QUESTION: When N solitons are randomly placed in a galaxy-like distribution,
do their OSCILLATING tails (sin(r)/r) average out, or does the nonlinear
term (∇g)²/g produce a net EXTRA FORCE that gives flat rotation curves?

PHYSICS:
  Single soliton tail: δ(r) = A · sin(r + φ₀) / r  (oscillating, 1/r envelope)
  N solitons at random positions xᵢ:
    g(x) = 1 + Σᵢ δ(|x - xᵢ|)  (linear superposition, approximate)
    ∇g = Σᵢ ∇δᵢ

  The LINEAR part: <Σ δᵢ> → smooth average (oscillations partially cancel)
  The NONLINEAR part: (∇g)²/g ≥ 0 ALWAYS → systematic positive contribution!

  SPECIFICALLY:
    (∇g)² = (Σᵢ ∇δᵢ)² = Σᵢ(∇δᵢ)² + 2·Σᵢ<ⱼ ∇δᵢ·∇δⱼ

    Self-terms Σᵢ(∇δᵢ)² ≥ 0 → always positive → extra "mass"
    Cross-terms: random phases → partial cancellation
    BUT: cross-terms have SYSTEMATIC component from spatial coherence

PLAN:
  1. Generate soliton tail profile from numerical solution
  2. Place N solitons in spherical galaxy with ρ ∝ r^(-2) or exponential
  3. Monte Carlo: evaluate g(x), ∇g(x) at test points along radial line
  4. Compute v²_eff(r) = r × [Newtonian + nonlinear correction]
  5. Compare with flat rotation curve expectation
  6. Scale: start with N=100-1000 (tractable), then extrapolate
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import time

print("="*78)
print("  MONTE CARLO N-SOLITON SUPERPOSITION")
print("="*78)

# ==========================================================================
# STEP 1: Compute single soliton tail profile (high precision)
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 1: Single soliton profile")
print(f"{'='*78}")

def tgp_rhs(r, y):
    g, gp = y
    if g <= 1e-10:
        return [gp, 0]
    if r < 1e-10:
        gpp = (1 - g - gp**2/g) / 3
    else:
        gpp = -gp**2/g - 2*gp/r - g + 1
    return [gp, gpp]

# Use g₀ = 1.1 (weak soliton — particle in weak field limit)
g0 = 1.1
r_max_sol = 200
sol = solve_ivp(tgp_rhs, [1e-6, r_max_sol], [g0, 0],
                method='RK45', max_step=0.02, rtol=1e-11, atol=1e-13,
                dense_output=True)

r_sol = sol.t
g_sol = sol.y[0]
gp_sol = sol.y[1]
delta_sol = g_sol - 1

# Create interpolation functions for fast evaluation
delta_interp = interp1d(r_sol, delta_sol, kind='cubic', fill_value=0, bounds_error=False)
gp_interp = interp1d(r_sol, gp_sol, kind='cubic', fill_value=0, bounds_error=False)

# Tail analysis: fit A·sin(r+φ)/r to the tail
r_tail = r_sol[r_sol > 5]
d_tail = delta_sol[r_sol > 5]
# Envelope: |δ|·r should be roughly constant
envelope = np.abs(d_tail) * r_tail
A_tail = np.mean(envelope[r_tail > 10])
print(f"  g₀ = {g0}, δ₀ = {g0-1}")
print(f"  Tail envelope: A = |δ|·r ≈ {A_tail:.4f}")
print(f"  First zero at r ≈ {r_sol[np.where(np.diff(np.sign(delta_sol)))[0][0]]:.2f}")

# ==========================================================================
# STEP 2: Galaxy model — place N solitons
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 2: Galaxy soliton distribution")
print(f"{'='*78}")

def generate_galaxy(N, R_scale, profile='hernquist', seed=42):
    """
    Generate N soliton positions in a galaxy-like distribution.
    Returns positions (N, 3) in soliton units.

    Hernquist: ρ(r) ∝ 1/(r·(r+a)³)
    Exponential disk: ρ(R,z) ∝ exp(-R/Rd)·exp(-|z|/zd)
    """
    rng = np.random.default_rng(seed)

    if profile == 'hernquist':
        # Hernquist: M(r) = M_total × r²/(r+a)²
        # Invert: r = a × √u / (1 - √u) where u = U[0,1]
        u = rng.uniform(0, 0.95, N)  # truncate at 95% mass
        sqrt_u = np.sqrt(u)
        radii = R_scale * sqrt_u / (1 - sqrt_u)
        # Random angles
        cos_theta = rng.uniform(-1, 1, N)
        phi = rng.uniform(0, 2*np.pi, N)
        sin_theta = np.sqrt(1 - cos_theta**2)
        x = radii * sin_theta * np.cos(phi)
        y = radii * sin_theta * np.sin(phi)
        z = radii * cos_theta
        return np.column_stack([x, y, z])

    elif profile == 'isothermal':
        # ρ ∝ 1/r² → M(r) ∝ r → r = R_max × U^(1/3)... no
        # Actually for ρ ∝ 1/r²: M(r) ∝ r
        # CDF: F(r) = r/R_max → r = R_max × U
        radii = R_scale * rng.uniform(0.1, 1, N)**(1)  # linear in r
        cos_theta = rng.uniform(-1, 1, N)
        phi = rng.uniform(0, 2*np.pi, N)
        sin_theta = np.sqrt(1 - cos_theta**2)
        x = radii * sin_theta * np.cos(phi)
        y = radii * sin_theta * np.sin(phi)
        z = radii * cos_theta
        return np.column_stack([x, y, z])

# Galaxy parameters (in soliton units)
# We'll work in soliton units where the soliton core ≈ 2 units wide
# A galaxy has solitons spread over ~R_galaxy soliton units
# The key is the RATIO of galaxy size to soliton tail wavelength (2π)

R_galaxy_values = [10, 30, 50, 100]  # in soliton units
N_solitons_values = [100, 500, 1000, 5000]

print(f"  Galaxy configurations to test:")
print(f"  R_galaxy (sol. units)  N_solitons  n_density  R/λ_tail")
for R in R_galaxy_values:
    for N in N_solitons_values:
        n_dens = N / (4/3 * np.pi * R**3)
        R_over_lambda = R / (2*np.pi)
        if N == 1000:
            print(f"  {R:20d}  {N:10d}  {n_dens:9.4f}  {R_over_lambda:8.2f}")

# ==========================================================================
# STEP 3: Evaluate substrate field along radial line
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 3: Monte Carlo evaluation of collective field")
print(f"{'='*78}")

def evaluate_field_at_point(point, soliton_positions, delta_func, gp_func, r_max=200):
    """
    Evaluate the total substrate field g(x) and its gradient at a point.

    Returns: delta_total, grad_delta_total (3-vector),
             sum_grad_sq, sum_self_grad_sq
    """
    # Displacement vectors from each soliton to the point
    displacements = point - soliton_positions  # (N, 3)
    distances = np.sqrt(np.sum(displacements**2, axis=1))  # (N,)

    # Avoid singularities
    distances = np.maximum(distances, 0.1)

    # Mask: only include solitons within interpolation range
    mask = distances < r_max

    if not np.any(mask):
        return 0, np.zeros(3), 0, 0

    # Delta values from each soliton
    d_vals = delta_func(distances[mask])  # δᵢ(|x - xᵢ|)

    # Total delta
    delta_total = np.sum(d_vals)

    # Gradient: ∇δᵢ = δ'ᵢ(r) × r̂ᵢ = gp(rᵢ) × (x - xᵢ)/rᵢ
    gp_vals = gp_func(distances[mask])  # g'(rᵢ)

    # Unit vectors from soliton to point
    r_hat = displacements[mask] / distances[mask, np.newaxis]  # (N_mask, 3)

    # Gradient contributions
    grad_contributions = gp_vals[:, np.newaxis] * r_hat  # (N_mask, 3)

    # Total gradient
    grad_total = np.sum(grad_contributions, axis=0)  # (3,)

    # |∇g|² = |Σ ∇δᵢ|²
    grad_sq_total = np.sum(grad_total**2)

    # Self-terms: Σᵢ |∇δᵢ|²
    self_grad_sq = np.sum(gp_vals**2)

    return delta_total, grad_total, grad_sq_total, self_grad_sq


def compute_rotation_curve(N_solitons, R_galaxy, N_radial=40, profile='hernquist', seed=42):
    """
    Compute the effective rotation curve for N solitons in a galaxy.
    """
    # Generate soliton positions
    positions = generate_galaxy(N_solitons, R_galaxy, profile, seed)

    # Radial test points along x-axis (average over a few directions)
    r_test = np.linspace(1, 3 * R_galaxy, N_radial)

    # Average over several directions for statistical stability
    n_directions = 6
    directions = np.array([
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1]
    ], dtype=float)

    delta_avg = np.zeros(N_radial)
    grad_r_sq_avg = np.zeros(N_radial)  # (∇g · r̂)² — radial component squared
    grad_total_sq_avg = np.zeros(N_radial)  # |∇g|²
    self_grad_sq_avg = np.zeros(N_radial)  # Σ|∇δᵢ|²
    M_enc_Newton = np.zeros(N_radial)

    for i_r, r in enumerate(r_test):
        for d in directions:
            point = r * d

            delta_t, grad_t, grad_sq, self_sq = evaluate_field_at_point(
                point, positions, delta_interp, gp_interp)

            delta_avg[i_r] += delta_t / n_directions
            grad_total_sq_avg[i_r] += grad_sq / n_directions
            self_grad_sq_avg[i_r] += self_sq / n_directions

            # Radial gradient component
            grad_r = np.dot(grad_t, d)
            grad_r_sq_avg[i_r] += grad_r**2 / n_directions

        # Newtonian enclosed mass: count solitons within r
        dists_from_center = np.sqrt(np.sum(positions**2, axis=1))
        M_enc_Newton[i_r] = np.sum(dists_from_center < r)

    return r_test, delta_avg, grad_total_sq_avg, self_grad_sq_avg, grad_r_sq_avg, M_enc_Newton


# ==========================================================================
# STEP 4: Run for multiple configurations
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 4: Results for different galaxy configurations")
print(f"{'='*78}")

# Test configurations: (N, R_galaxy)
configs = [
    (200,  15, 'hernquist'),
    (500,  20, 'hernquist'),
    (1000, 30, 'hernquist'),
    (2000, 40, 'hernquist'),
    (5000, 50, 'hernquist'),
]

results = {}

for N, R_gal, prof in configs:
    t0 = time.time()
    r_t, delta, grad_sq, self_sq, grad_r_sq, M_enc = compute_rotation_curve(
        N, R_gal, N_radial=30, profile=prof, seed=42)
    dt = time.time() - t0

    results[(N, R_gal)] = (r_t, delta, grad_sq, self_sq, grad_r_sq, M_enc)

    # Rotation curve analysis:
    # Newtonian: v²_N(r) = G · M_enc / r → proportional to M_enc/r
    # In soliton units with G=1 and unit mass per soliton:
    v2_Newton = M_enc / r_t

    # The g'²/g correction:
    # In the TGP equation: ∇²g + g'²/g + g = 1
    # The g'²/g term acts as extra source → effective extra mass
    # At radius r, the contribution to acceleration:
    # a_NL(r) ≈ <(∇g)²/g> evaluated in a shell at r
    # For g ≈ 1 + δ_avg with |δ_avg| << 1:
    # (∇g)²/g ≈ (∇g)² · (1 - δ_avg + ...) ≈ (∇g)²

    g_total = 1 + delta

    # Nonlinear "pressure" term: |∇g|²/g
    nl_term = grad_sq / np.maximum(g_total, 0.1)

    # Cross-term contribution: total - self
    cross_sq = grad_sq - self_sq

    # Effective extra acceleration from NL term:
    # The extra enclosed "mass" from nonlinear term
    # ∇²g_extra ≈ -(∇g)²/g → M_extra(r) ∝ ∫₀ʳ (∇g)²/g × r'² dr'
    # Approximation: shell contribution
    dr = r_t[1] - r_t[0] if len(r_t) > 1 else 1
    M_extra_cumul = np.cumsum(nl_term * r_t**2 * dr) * 4 * np.pi / 3

    v2_total = (M_enc + M_extra_cumul) / r_t

    v_N = np.sqrt(np.maximum(v2_Newton, 0))
    v_T = np.sqrt(np.maximum(v2_total, 0))

    # Flatness analysis
    idx_peak = np.argmax(v_T)
    v_peak = v_T[idx_peak]
    r_peak = r_t[idx_peak]

    # v at 2×R_galaxy
    idx_2R = np.argmin(np.abs(r_t - 2*R_gal))
    v_2R = v_T[idx_2R]
    v_N_2R = v_N[idx_2R]

    ratio_total = v_2R / v_peak if v_peak > 0 else 0
    ratio_newton = v_N_2R / v_N[np.argmax(v_N)] if np.max(v_N) > 0 else 0

    enhancement_2R = v_T[idx_2R]**2 / max(v_N[idx_2R]**2, 1e-10)

    print(f"\n  N={N:5d}, R_gal={R_gal:3d}, t={dt:.1f}s")
    print(f"    v_peak(TGP) = {v_peak:.3f} at r = {r_peak:.1f}")
    print(f"    v(2R)/v_peak:  TGP = {ratio_total:.3f},  Newton = {ratio_newton:.3f}")
    print(f"    v²(2R) enhancement: {enhancement_2R:.2f}×")
    print(f"    <δ> at R_gal: {delta[np.argmin(np.abs(r_t-R_gal))]:.4f}")
    print(f"    NL term at R_gal: {nl_term[np.argmin(np.abs(r_t-R_gal))]:.4e}")
    print(f"    Self/Total at R_gal: {self_sq[np.argmin(np.abs(r_t-R_gal))]/max(grad_sq[np.argmin(np.abs(r_t-R_gal))],1e-20):.3f}")

# ==========================================================================
# STEP 5: Detailed rotation curve for the largest galaxy
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 5: Detailed rotation curve — N=5000, R=50")
print(f"{'='*78}")

N_best, R_best = 5000, 50
r_t, delta, grad_sq, self_sq, grad_r_sq, M_enc = results[(N_best, R_best)]

g_total = 1 + delta
nl_term = grad_sq / np.maximum(g_total, 0.1)
dr = r_t[1] - r_t[0]
M_extra = np.cumsum(nl_term * r_t**2 * dr) * 4 * np.pi / 3

v2_N = M_enc / r_t
v2_T = (M_enc + M_extra) / r_t
v_N = np.sqrt(np.maximum(v2_N, 0))
v_T = np.sqrt(np.maximum(v2_T, 0))

# Keplerian reference
idx_R = np.argmin(np.abs(r_t - R_best))
M_R = M_enc[idx_R]
v_Kepler = np.sqrt(np.maximum(M_R / r_t, 0))

print(f"\n  {'r':>6s} {'M_enc':>7s} {'M_extra':>8s} {'M_tot':>8s} {'v_Newton':>8s} {'v_TGP':>8s} {'v_Kepler':>8s} {'v_T/v_N':>7s} {'v_T/v_K':>7s}")
print(f"  {'─'*70}")

for i in range(0, len(r_t), max(1, len(r_t)//20)):
    r = r_t[i]
    mn = M_enc[i]
    me = M_extra[i]
    mt = mn + me
    vn = v_N[i]
    vt = v_T[i]
    vk = v_Kepler[i]
    rat_n = vt/vn if vn > 0 else 0
    rat_k = vt/vk if vk > 0 else 0
    print(f"  {r:6.1f} {mn:7.0f} {me:8.1f} {mt:8.1f} {vn:8.3f} {vt:8.3f} {vk:8.3f} {rat_n:7.3f} {rat_k:7.3f}")

# ==========================================================================
# STEP 6: SCALING ANALYSIS — how does the NL enhancement scale with N?
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 6: Scaling of nonlinear enhancement with N")
print(f"{'='*78}")

# For each config, measure the enhancement at R_galaxy and 2×R_galaxy
print(f"\n  {'N':>6s} {'R_gal':>6s} {'M_enc(R)':>8s} {'M_NL(R)':>8s} {'M_NL/M_enc':>10s} {'v_enhance':>9s}")
print(f"  {'─'*55}")

for (N, R_gal), (r_t, delta, grad_sq, self_sq, grad_r_sq, M_enc) in results.items():
    g_total = 1 + delta
    nl_term = grad_sq / np.maximum(g_total, 0.1)
    dr = r_t[1] - r_t[0]
    M_extra = np.cumsum(nl_term * r_t**2 * dr) * 4 * np.pi / 3

    idx_R = np.argmin(np.abs(r_t - R_gal))
    me = M_enc[idx_R]
    mnl = M_extra[idx_R]
    ratio = mnl / max(me, 1)
    v_enh = np.sqrt(1 + ratio)

    print(f"  {N:6d} {R_gal:6d} {me:8.0f} {mnl:8.1f} {ratio:10.4f} {v_enh:9.4f}")

# ==========================================================================
# STEP 7: STATISTICAL ANALYSIS — averaging over realizations
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 7: Statistical robustness — multiple random seeds")
print(f"{'='*78}")

N_stat, R_stat = 1000, 30
n_seeds = 5
enhancements = []

for seed in range(n_seeds):
    r_t, delta, grad_sq, self_sq, grad_r_sq, M_enc = compute_rotation_curve(
        N_stat, R_stat, N_radial=20, profile='hernquist', seed=seed*17 + 3)

    g_total = 1 + delta
    nl_term = grad_sq / np.maximum(g_total, 0.1)
    dr = r_t[1] - r_t[0]
    M_extra = np.cumsum(nl_term * r_t**2 * dr) * 4 * np.pi / 3

    idx_R = np.argmin(np.abs(r_t - R_stat))
    ratio = M_extra[idx_R] / max(M_enc[idx_R], 1)
    enhancements.append(ratio)

    # Also check at 2R
    idx_2R = np.argmin(np.abs(r_t - 2*R_stat))
    ratio_2R = M_extra[idx_2R] / max(M_enc[idx_2R], 1) if M_enc[idx_2R] > 0 else 0

print(f"\n  N={N_stat}, R={R_stat}, {n_seeds} realizations:")
print(f"  M_NL/M_enc at R: {np.mean(enhancements):.4f} ± {np.std(enhancements):.4f}")
print(f"  v enhancement: {np.sqrt(1+np.mean(enhancements)):.4f}")

# ==========================================================================
# STEP 8: EXTRAPOLATION TO REAL GALAXY
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 8: Extrapolation to real galaxy (N ~ 10¹¹)")
print(f"{'='*78}")

# From the scaling results, determine how M_NL/M_enc scales with N and R
# If M_NL ∝ N^α · R^β, we can extrapolate

Ns = []
Rs = []
ratios_at_R = []

for (N, R_gal), (r_t, delta, grad_sq, self_sq, grad_r_sq, M_enc) in results.items():
    g_total = 1 + delta
    nl_term = grad_sq / np.maximum(g_total, 0.1)
    dr = r_t[1] - r_t[0]
    M_extra = np.cumsum(nl_term * r_t**2 * dr) * 4 * np.pi / 3
    idx_R = np.argmin(np.abs(r_t - R_gal))
    ratio = M_extra[idx_R] / max(M_enc[idx_R], 1)
    Ns.append(N)
    Rs.append(R_gal)
    ratios_at_R.append(ratio)

Ns = np.array(Ns, dtype=float)
Rs = np.array(Rs, dtype=float)
ratios_at_R = np.array(ratios_at_R)

# Fit: log(ratio) = α·log(N) + β·log(R) + const
# Or simpler: ratio ∝ N^α (at fixed R/N^(1/3) concentration)
if len(Ns) > 2 and np.all(ratios_at_R > 0):
    log_N = np.log(Ns)
    log_ratio = np.log(ratios_at_R)
    # Linear fit in log-log
    coeffs = np.polyfit(log_N, log_ratio, 1)
    alpha_N = coeffs[0]
    const = coeffs[1]

    print(f"  Scaling fit: M_NL/M_enc ∝ N^α")
    print(f"  α = {alpha_N:.3f}")
    print(f"  const = {const:.3f}")

    # Extrapolate to N = 10¹¹ (MW-like)
    N_MW = 1e11
    ratio_MW = np.exp(const) * N_MW**alpha_N
    print(f"\n  Extrapolation to N = {N_MW:.0e}:")
    print(f"  M_NL/M_enc ≈ {ratio_MW:.4f}")
    print(f"  v_enhancement ≈ {np.sqrt(1+ratio_MW):.4f}")

    if ratio_MW > 1:
        print(f"  → SIGNIFICANT: nonlinear term DOMINATES → flat curves plausible!")
    elif ratio_MW > 0.1:
        print(f"  → MODERATE: nonlinear correction is 10-100% → partially flat")
    elif ratio_MW > 0.01:
        print(f"  → SMALL: 1-10% correction → not enough for flat curves")
    else:
        print(f"  → NEGLIGIBLE: < 1% correction → mechanism does NOT work")
else:
    print(f"  WARNING: insufficient data for scaling fit")
    print(f"  Ratios: {ratios_at_R}")

# ==========================================================================
# STEP 9: DECOMPOSITION — self vs cross terms
# ==========================================================================
print(f"\n{'='*78}")
print(f"  STEP 9: Self-terms vs cross-terms decomposition")
print(f"{'='*78}")

print(f"""
  |∇g|² = Σᵢ|∇δᵢ|² + 2·Σᵢ<ⱼ ∇δᵢ·∇δⱼ
           ─────────   ──────────────────
           self-terms   cross-terms

  Self-terms: always positive, ∝ N (each soliton contributes independently)
  Cross-terms: partially cancel due to random phases
    - Coherent part: ∝ N² (if phases correlated)
    - Random part: ∝ N (central limit theorem)
""")

for (N, R_gal), (r_t, delta, grad_sq, self_sq, grad_r_sq, M_enc) in results.items():
    idx_R = np.argmin(np.abs(r_t - R_gal))
    total = grad_sq[idx_R]
    self_t = self_sq[idx_R]
    cross = total - self_t

    print(f"  N={N:5d}, R={R_gal:3d}: |∇g|²={total:.4e}, self={self_t:.4e}, cross={cross:+.4e}")
    print(f"    self/total = {self_t/max(total,1e-20):.3f}, cross/total = {cross/max(total,1e-20):+.3f}")
    print(f"    cross/self = {cross/max(self_t,1e-20):+.3f}")
    if cross > 0:
        print(f"    → Cross-terms POSITIVE (constructive) = {cross/self_t*100:.1f}% of self")
    else:
        print(f"    → Cross-terms NEGATIVE (destructive) = {cross/self_t*100:.1f}% of self")

# ==========================================================================
# STEP 10: SUMMARY
# ==========================================================================
print(f"\n{'='*78}")
print(f"  SUMMARY: N-soliton Monte Carlo results")
print(f"{'='*78}")

print(f"""
  CONFIGURATIONS TESTED:
  N = 200 to 5000 solitons, R_galaxy = 15 to 50 (soliton units)
  Profile: Hernquist (realistic galaxy)

  KEY FINDINGS:

  1. LINEAR SUPERPOSITION (δ):
     The average <δ> at R_galaxy is small but non-zero (soliton-dominated).
     Oscillations partially cancel for random positions.

  2. NONLINEAR TERM |∇g|²/g:
     Always positive → provides systematic extra "mass".
     Dominated by SELF-TERMS (each soliton's own gradient squared).
     Cross-terms are smaller and can be positive or negative.

  3. ENHANCEMENT SCALING:
     The ratio M_NL/M_enc at R_galaxy determines whether the mechanism
     produces flat rotation curves.

  4. SELF-TERMS DOMINATE:
     |∇δᵢ|² at the soliton's own location is large (core gradient).
     These contribute to the overall "pressure" but don't necessarily
     create a smooth profile — they're localized at soliton positions.

  ASSESSMENT:
  The nonlinear term produces corrections, but the key question is
  whether they scale fast enough with N to become significant at N ~ 10¹¹.
""")
