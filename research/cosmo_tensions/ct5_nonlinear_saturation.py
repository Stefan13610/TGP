#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ct5_nonlinear_saturation.py: Nonlinear saturation of TGP tachyonic mode.

KEY INSIGHT from ct4: linear perturbation analysis EXPLODES (B/H₀² ~ 10⁹).
This means:
1. The tachyonic instability at ψ=1 is REAL and STRONG
2. The nonlinear term 3ψ̇²/ψ MUST provide saturation
3. The actual δψ is set by BALANCE: tachyon ↔ geometric friction

SATURATION ANALYSIS:
At equilibrium (homogeneous, no spatial gradients):
  3H·ψ̇ + 3ψ̇²/ψ ≈ γc₀²(ψ - 1)

This is a NONLINEAR ODE for ψ(t). The "field" oscillates around ψ=1
with amplitude set by the tachyonic energy vs friction dissipation.

PLAN:
1. Solve the FULL nonlinear background equation (not linearized)
   with various initial kicks δψ₀
2. Find the equilibrium oscillation amplitude
3. Compute σ²_ψ from the time-averaged ψ variance
4. Compute B_ψ = 3<ψ̇²/ψ> - 3<ψ̇>²/<ψ>
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.integrate import solve_ivp
import math

print("="*78)
print("  NONLINEAR SATURATION OF TGP TACHYONIC MODE")
print("="*78)

# ==========================================================================
# PARAMETERS
# ==========================================================================
H0 = 67.4
Om_m = 0.315
Om_L = 0.685
gamma_H0 = 36 * Om_L   # γc₀²/H₀² ≈ 24.7
m_tach = math.sqrt(gamma_H0)  # ≈ 4.97

print(f"\n  γc₀²/H₀² = {gamma_H0:.2f}")
print(f"  m_tach/H₀ = {m_tach:.3f}")

# ==========================================================================
# PART 1: Nonlinear oscillation around ψ = 1
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 1: Full nonlinear dynamics near ψ=1")
print("="*78)

def H_over_H0(a):
    return np.sqrt(Om_m * a**(-3) + Om_L)

# Full nonlinear field equation in ln(a):
# dψ/d(ln a) = χ
# dχ/d(ln a) = -(3 + d ln H/d ln a)χ - 3χ²/ψ + γ(ψ-1)/H²(a)
#
# But for the saturation analysis, we can work at FIXED H (quasi-static)
# since the saturation timescale (1/m_tach ~ 0.2/H₀) << Hubble time

# Simplification: at fixed a=1 (today), H = H₀√Ω_Λ
H_today = H_over_H0(1.0)

def nonlinear_rhs_t(t_H0, y):
    """
    Full nonlinear field equation in cosmic time (units of H₀⁻¹).
    y = [ψ, dψ/dt (in H₀ units)]
    """
    psi, psi_dot = y

    if psi < 1e-6: return [psi_dot, 0]  # prevent blowup

    H = H_today  # quasi-static approximation

    # ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = γ(ψ - 1)
    psi_ddot = -3 * H * psi_dot - 3 * psi_dot**2 / psi + gamma_H0 * (psi - 1)

    return [psi_dot, psi_ddot]

# Solve with various initial kicks
print(f"\n  Solving nonlinear equation with H = H₀·{H_today:.4f} (today)")
print(f"  ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = γ(ψ - 1)")

kicks = [0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
t_span = [0, 20]  # in H₀⁻¹ units (several Hubble times)

print(f"\n  {'δψ₀':>8s}  {'ψ_final':>10s}  {'ψ̇_final':>12s}  {'<(ψ-1)²>':>12s}  {'<ψ̇²/ψ>':>12s}  {'Stable?':>8s}")

for kick in kicks:
    psi0 = 1.0 + kick
    y0 = [psi0, 0.0]

    sol = solve_ivp(nonlinear_rhs_t, t_span, y0, method='DOP853',
                    rtol=1e-12, atol=1e-15, max_step=0.01, dense_output=True)

    if sol.success:
        # Compute time averages over last 10 H₀⁻¹
        t_avg = np.linspace(10, 20, 10000)
        psi_t = sol.sol(t_avg)[0]
        psidot_t = sol.sol(t_avg)[1]

        psi_final = psi_t[-1]
        psidot_final = psidot_t[-1]
        var_psi = np.mean((psi_t - 1)**2)
        mean_psidot2_over_psi = np.mean(psidot_t**2 / psi_t)

        stable = "YES" if abs(psi_final - 1) < 0.1 and abs(psidot_final) < 0.1 else "NO"
        if np.any(psi_t < 0): stable = "CRASH"

        print(f"  {kick:8.3f}  {psi_final:10.6f}  {psidot_final:12.6f}  {var_psi:12.6e}  {mean_psidot2_over_psi:12.6e}  {stable:>8s}")

# ==========================================================================
# PART 2: Attractor behavior — where does ψ end up?
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 2: Long-term evolution (attractor)")
print("="*78)

# Solve for 100 Hubble times to find attractor
t_long = [0, 100]

for kick in [0.01, 0.1, 0.5, 1.0]:
    psi0 = 1.0 + kick
    sol = solve_ivp(nonlinear_rhs_t, t_long, [psi0, 0.0], method='DOP853',
                    rtol=1e-12, atol=1e-15, max_step=0.01, dense_output=True)

    if sol.success:
        t_sample = np.linspace(80, 100, 5000)
        psi_s = sol.sol(t_sample)[0]
        psidot_s = sol.sol(t_sample)[1]

        psi_mean = np.mean(psi_s)
        psi_std = np.std(psi_s)
        psidot_rms = np.sqrt(np.mean(psidot_s**2))

        print(f"  kick={kick:.2f}: <ψ>={psi_mean:.8f}, σ_ψ={psi_std:.2e}, "
              f"ψ̇_rms={psidot_rms:.2e}")

# ==========================================================================
# PART 3: Include matter source continuously
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 3: Continuously sourced nonlinear equation")
print("="*78)

# In reality, matter perturbations continuously kick the field.
# Model: ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = γ(ψ-1) + S(t)
# where S(t) represents stochastic matter forcing
#
# For a simple model: S = S₀ (constant source, representing average forcing)
# Then the equilibrium is: γ(ψ_eq - 1) = -S₀ → ψ_eq = 1 - S₀/γ
#
# The source amplitude: S ~ 4πGρ_m/c₀² ~ (3/2)Ω_m H₀²
# In our units: S₀ = 1.5 × Om_m = 0.4725

S0 = 1.5 * Om_m  # matter source in H₀² units
print(f"\n  Matter source amplitude: S₀ = 1.5·Ω_m = {S0:.4f} (H₀² units)")

# New equilibrium: γ(ψ_eq - 1) + S₀ = 0
# → ψ_eq = 1 - S₀/γ = 1 - 0.4725/24.66 = 0.981
psi_eq = 1 - S0 / gamma_H0
print(f"  Equilibrium shift: ψ_eq = 1 - S₀/γ = {psi_eq:.6f}")
print(f"  δψ_eq = {psi_eq - 1:.6f}")

# This is tiny! The field barely moves from 1.
# BUT: the fluctuations around this equilibrium are what matter.

# Solve with continuous source
def sourced_nonlinear(t_H0, y):
    psi, psi_dot = y
    if psi < 1e-6: return [psi_dot, 0]
    H = H_today
    psi_ddot = -3*H*psi_dot - 3*psi_dot**2/psi + gamma_H0*(psi-1) + S0
    return [psi_dot, psi_ddot]

sol_sourced = solve_ivp(sourced_nonlinear, [0, 100], [1.0, 0.0],
                        method='DOP853', rtol=1e-12, atol=1e-15,
                        max_step=0.01, dense_output=True)

if sol_sourced.success:
    t_s = np.linspace(50, 100, 10000)
    psi_s = sol_sourced.sol(t_s)[0]
    psidot_s = sol_sourced.sol(t_s)[1]

    print(f"\n  With constant source S₀ = {S0:.4f}:")
    print(f"    <ψ> = {np.mean(psi_s):.8f}")
    print(f"    σ_ψ = {np.std(psi_s):.2e}")
    print(f"    ψ_eq (predicted) = {psi_eq:.6f}")
    print(f"    ψ̇_rms = {np.sqrt(np.mean(psidot_s**2)):.2e}")

# ==========================================================================
# PART 4: The REAL backreaction — inhomogeneous ψ
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 4: Backreaction from spatial inhomogeneity")
print("="*78)

print(f"""
  The previous parts show that the HOMOGENEOUS ψ is stable at ψ_eq ≈ 1.
  The backreaction comes from SPATIAL VARIATIONS of ψ.

  In regions with matter (galaxies, clusters): ψ pulled to ψ < 1 (deficit)
    or ψ > 1 (excess) depending on soliton type
  In voids: ψ → 1 (vacuum state)

  The variance σ²_ψ is NOT from temporal oscillations,
  but from SPATIAL inhomogeneity: ψ(x) varies across space.

  Key: the TGP field equation with SPATIAL gradients:
    ψ̈ + 3Hψ̇ + 3ψ̇²/ψ = c₀²[∇²ψ/(a²) + γ(ψ-1)]

  The ∇²ψ term acts as RESTORING FORCE for spatial fluctuations.
  At equilibrium (quasi-static):
    c₀²∇²ψ/a² ≈ -γc₀²(ψ-1) - matter source
    → ψ - 1 ≈ Φ_N/c₀²  (Newtonian potential)

  So δψ ~ Φ/c² ~ 10⁻⁵ (as expected from weak field limit).
  The variance: σ²_ψ ~ <(Φ/c²)²> ~ 10⁻¹⁰.

  BUT: this is for the STATIC (equilibrium) response.
  The DYNAMIC backreaction from 3ψ̇²/ψ involves TIME DERIVATIVES.

  In TGP, ψ̇ ≠ 0 even in quasi-equilibrium because:
  1. Structure formation: ψ changes as matter accretes
  2. Hubble flow: cosmological expansion changes ψ
  3. Tachyonic amplification: perturbations grow on m_tach ~ 5H₀ timescale

  The dynamic term ψ̇ at each point:
  ψ̇ = ψ̇_bg + δψ̇_growth + δψ̇_infall
  where:
  - ψ̇_bg ≈ 0 (background stable at ψ=1)
  - δψ̇_growth ~ m_tach·δψ ~ 5H₀·Φ/c² ~ 5H₀·10⁻⁵ = 5×10⁻⁵ H₀
  - δψ̇_infall ~ H·f·δψ ~ H₀·0.5·10⁻⁵ ~ 0.5×10⁻⁵ H₀

  So ψ̇ ~ 5×10⁻⁵ H₀ (dominated by tachyonic growth).
""")

# Backreaction: B_ψ = 3[<ψ̇²/ψ> - <ψ̇>²/<ψ>]
# With δψ ~ 10⁻⁵ and δψ̇ ~ m·δψ ~ 5·10⁻⁵ H₀:
# <ψ̇²/ψ> ≈ <(δψ̇)²>/1 ≈ (5·10⁻⁵)² H₀² = 25·10⁻¹⁰ H₀²
# B_ψ ≈ 3 × 25·10⁻¹⁰ H₀² = 7.5×10⁻⁹ H₀²

delta_psi_typ = 1e-5  # Φ/c² typical
m_tach_val = m_tach    # ~ 5 H₀
psi_dot_typ = m_tach_val * delta_psi_typ  # H₀ units

B_naive = 3 * psi_dot_typ**2
print(f"  Naive (linear) backreaction:")
print(f"    δψ ~ {delta_psi_typ:.1e} (from Φ/c²)")
print(f"    δψ̇ ~ m·δψ = {m_tach_val:.1f}·{delta_psi_typ:.1e} = {psi_dot_typ:.1e} (H₀ units)")
print(f"    B_ψ/H₀² ~ 3·(δψ̇)² = {B_naive:.2e}")
print(f"    Required: 0.174")
print(f"    → SMALL by factor {0.174/B_naive:.0e}")

# ==========================================================================
# PART 5: THE KEY — nonlinear amplification beyond weak field
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 5: Beyond weak field — nonlinear TGP substrate")
print("="*78)

print(f"""
  The weak field approximation δψ ≈ 2Φ/c² gives δψ ~ 10⁻⁵.
  This is TOO SMALL for significant backreaction.

  BUT: TGP is NOT just "GR with a scalar field."
  The key difference:

  In GR: Φ is determined by Poisson ∇²Φ = 4πGρ → Φ ~ GM/(Rc²)
  In TGP: ψ is determined by ∇²ψ + γ(ψ-1) = -source

  The γ(ψ-1) term is the TACHYONIC POTENTIAL.
  For the Poisson-like equation with tachyonic mass:
    (∇² - m²)ψ = -source    where m² = -γ < 0

  This is a SCREENED Poisson equation, but with NEGATIVE mass² (anti-screening).
  The solution is NOT the usual Yukawa exp(-mr)/r but OSCILLATORY/growing.

  For a point source in this equation:
    ψ(r) - 1 = (G·M / r) · exp(+μr) / c₀²    where μ = √γ/c₀

  The exp(+μr) factor means the field GROWS with distance!
  This is bounded only by the cosmological horizon (H⁻¹).

  At the Hubble radius R_H ~ c₀/H₀:
    enhancement = exp(μ·R_H) = exp(√γ · c₀/(c₀·H₀)) = exp(√(γ/H₀²))
                = exp({m_tach:.1f}) = {math.exp(m_tach):.1e}

  This gives: δψ at Hubble scale ~ (Φ_N/c²) × exp(m_tach/H₀)
            ~ 10⁻⁵ × {math.exp(m_tach):.0e} = {1e-5 * math.exp(m_tach):.2e}

  WAIT — this means the STATIC solution for ψ already has
  exponentially enhanced tails! δψ is NOT 10⁻⁵ but much larger!
""")

# Check: solve static TGP equation with point source
# ψ'' + (2/r)ψ' - γ·ψ = 0   (homogeneous, far from source)
# ψ(r) = A·sinh(μr)/(μr) + B·cosh(μr)/(μr)
# where μ = √γ / c₀ (in proper distance)
# In comoving coordinates with H₀ units:
mu_comov = math.sqrt(gamma_H0)  # = m_tach ~ 5 (in units of H₀/c₀)

print(f"\n  Static TGP field around a mass:")
print(f"    ψ - 1 ~ (GM/rc₀²) · exp(μr)/(μr)")
print(f"    μ = √γ/c₀ = {mu_comov:.3f} H₀/c₀")
print(f"    At r = 1 Mpc: μr = {mu_comov * 1/4451:.6f}")
print(f"    At r = 100 Mpc: μr = {mu_comov * 100/4451:.4f}")
print(f"    At r = 4451 Mpc (c/H₀): μr = {mu_comov:.2f}")
print(f"    exp(μr) at c/H₀ = {math.exp(mu_comov):.1e}")

# This exponential growth of the static solution is UNPHYSICAL
# for a single source — it must be cut off by:
# 1. Other sources (screening by large-scale structure)
# 2. Cosmological boundary conditions (ψ → 1 at infinity)
# 3. The nonlinear 3ψ̇²/ψ term (dynamic stabilization)

# The ACTUAL static field is determined by:
# The LINEAR operator (∇² - m²) has NEGATIVE eigenvalues
# → the static solution does NOT exist in the usual sense
# → the field MUST be time-dependent (oscillating/growing)

# This is the SAME tachyonic instability seen in Part 2 of ct4.
# It confirms: ψ cannot sit at ψ = 1 in an inhomogeneous universe.

# ==========================================================================
# PART 6: Self-consistent equilibrium
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 6: Self-consistent equilibrium state")
print("="*78)

print(f"""
  The tachyonic instability means ψ = 1 is UNSTABLE when perturbed.
  The field must settle into a DIFFERENT state, determined by:

  1. Tachyonic potential pushes ψ away from 1
  2. Geometric friction 3ψ̇²/ψ resists change
  3. Hubble damping 3Hψ̇ dissipates energy
  4. Spatial gradients ∇²ψ provide restoring force at small scales

  The self-consistent state has ψ ≠ 1, with amplitude
  set by the BALANCE of these forces.

  For a SPHERICAL void (no matter, radius R_void):
    Inside void: ψ = 1 (vacuum, stable)
    At void boundary: ψ jumps to ψ_wall = 1 + δψ_wall

  For a HALO (mass M, radius R_halo):
    ψ_halo = 1 + 2Φ/c₀² + corrections from tachyonic term
    Corrections enhance δψ by factor ~ 1 + γR²/(c₀²/a²)
    For R_halo ~ 1 Mpc: γR²·(a/c₀)² = {gamma_H0} × (1/{4451})² ≈ {gamma_H0/(4451**2):.2e}
    → correction negligible at galaxy/cluster scales!

  KEY REALIZATION:
  The tachyonic mass μ ~ 5H₀/c₀ corresponds to
  λ_tach = 2π/μ ~ {2*math.pi/mu_comov * 4451:.0f} Mpc.

  At scales R << λ_tach ({2*math.pi/mu_comov * 4451:.0f} Mpc):
    The tachyonic term is NEGLIGIBLE → standard Newtonian gravity
    δψ ≈ 2Φ/c₀² ~ 10⁻⁵ (unchanged)

  At scales R >> λ_tach:
    The tachyonic term DOMINATES → exponential enhancement
    But these are SUPER-HUBBLE scales → not observable!

  CONCLUSION: The tachyonic instability operates at COSMOLOGICAL scales
  (hundreds to thousands of Mpc), not at galaxy/cluster scales.
""")

# ==========================================================================
# PART 7: Cosmological backreaction — correct accounting
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 7: Correct backreaction computation")
print("="*78)

# The backreaction B_ψ depends on the variance of ψ over the HUBBLE VOLUME.
# At galaxy/cluster scales: δψ ~ Φ/c₀² ~ 10⁻⁵ (standard)
# At void scales (100-300 Mpc): δψ enhanced by tachyonic effect

# The TOTAL variance:
# σ²_ψ = ∫₀^∞ P_δψ(k) k² dk/(2π²)
#
# P_δψ(k) = |T_TGP(k)|² × P_Φ(k)
# where T_TGP(k) is the TGP transfer function:
# T_TGP(k) = k² / (k² - μ²)  (for static, μ² < 0 → enhancement)
#
# At k >> μ: T_TGP → 1 (standard)
# At k → μ: T_TGP → ∞ (resonance!)
# At k < μ: T_TGP → k²/(k²-μ²) (enhanced)

# The integral has a POLE at k = μ (tachyonic resonance).
# In practice, this pole is regularized by:
# - Time-dependence (growth rate finite)
# - Nonlinear saturation
# - Finite age of universe

# Regularized: T_TGP(k) ≈ k² / (k² - μ² + iΓμ)
# where Γ ~ 3H (Hubble damping) → effective width

Gamma = 3 * H_over_H0(1.0)  # damping rate in H₀ units
mu_sq = gamma_H0  # μ² in (H₀/c₀)²

# The enhanced variance near the pole:
# Δσ²_ψ ~ |T(k≈μ)|² × P_Φ(μ) × μ² × Δk/(2π²)
# |T(μ)| ≈ μ²/|iΓμ| = μ/Γ
T_at_pole = mu_comov / Gamma
print(f"\n  Transfer function at tachyonic pole:")
print(f"    μ = √γ = {mu_comov:.3f} H₀/c₀")
print(f"    Γ = 3H = {Gamma:.3f} H₀")
print(f"    |T(k≈μ)| ≈ μ/Γ = {T_at_pole:.3f}")
print(f"    Enhancement: |T|² = {T_at_pole**2:.2f}")

# The gravitational potential power spectrum at k ~ μ ~ 5 H₀/c₀:
# λ ~ 2π/μ · (c₀/H₀) ~ 5600 Mpc → SUPER-HUBBLE
# P_Φ(k) ≈ (3Ωm/2)² × (H₀/k)⁴ × P_matter(k)
# At super-Hubble scales: P_matter(k) ~ k^(n_s-1) · T²(k)
# T(k) → 1 for k << k_eq, so P_matter ~ k^(n_s-1) ~ k^(-0.033)

# The Φ power spectrum on large scales:
# P_Φ ~ As × (k/k_pivot)^(ns-1) × (3Ωm/2)² × (H₀/k)⁴
# where As ≈ 2.1×10⁻⁹ (scalar amplitude), k_pivot = 0.05/Mpc
As = 2.1e-9  # scalar amplitude
k_pivot_H0 = 0.05 * 4451  # k_pivot in H₀/c₀ units ≈ 223

# At k = μ ≈ 5 (H₀/c₀ units):
P_Phi_at_mu = As * (mu_comov / k_pivot_H0)**(0.967-1) * (1.5*Om_m)**2 * (1/mu_comov)**4
delta_psi_sq_at_mu = T_at_pole**2 * P_Phi_at_mu * mu_comov**2  # per log k

print(f"\n  Power spectrum at tachyonic scale:")
print(f"    k_pivot = {k_pivot_H0:.0f} (H₀/c₀)")
print(f"    As = {As:.1e}")
print(f"    P_Φ(k=μ) ≈ {P_Phi_at_mu:.2e}")
print(f"    P_δψ(k=μ) = |T|² × P_Φ = {T_at_pole**2 * P_Phi_at_mu:.2e}")
print(f"    Variance per Δln k ~ 1: {delta_psi_sq_at_mu:.2e}")

# The backreaction:
# B_ψ/H₀² ≈ 3 × m_tach² × σ²_ψ  (since ψ̇ ~ m·δψ)
B_from_pole = 3 * mu_comov**2 * delta_psi_sq_at_mu
print(f"\n  Backreaction from tachyonic pole:")
print(f"    B_ψ/H₀² ≈ 3 × m² × σ²_ψ = {B_from_pole:.2e}")
print(f"    Required: 0.174")

# ==========================================================================
# PART 8: Enhanced model — proper integral
# ==========================================================================
print(f"\n{'='*78}")
print("  PART 8: Full variance integral with TGP transfer function")
print("="*78)

# Integrate σ²_ψ = ∫ |T(k)|² × P_Φ(k) × k² dk/(2π²)
# with T(k) = k²/(k² - μ² + iΓμ)

from scipy.integrate import quad

def integrand_sigma2(lnk):
    """Integrand for σ²_ψ in ln k."""
    k = np.exp(lnk)

    # TGP transfer function (regularized)
    T_sq = (k**4) / ((k**2 - mu_sq)**2 + (Gamma*mu_comov)**2)

    # Matter power spectrum (simplified: scale-invariant + Eisenstein-Hu transfer)
    # P_Φ(k) = As × (k/k_pivot)^(ns-1) × (1.5Ωm)² × (1/k)⁴ × T²_EH(k)
    # For simplicity: T_EH ≈ 1 for k < k_eq, T_EH ~ (k_eq/k)² for k > k_eq
    k_eq = 0.01 * 4451  # matter-radiation equality in H₀/c₀ units ≈ 45

    if k < k_eq:
        T_matter = 1.0
    else:
        T_matter = (k_eq / k)**2

    P_Phi = As * (k/k_pivot_H0)**(0.967-1) * (1.5*Om_m)**2 * (1/k)**4 * T_matter**2

    # Integrand: P_δψ × k³/(2π²) = |T|² × P_Φ × k³/(2π²)
    return T_sq * P_Phi * k**3 / (2 * math.pi**2)

# Integrate from k = 0.001 to k = 100 (in H₀/c₀ units)
result, error = quad(integrand_sigma2, np.log(0.001), np.log(100),
                     limit=200, epsabs=1e-15, epsrel=1e-8)

sigma2_psi_total = result
print(f"\n  σ²_ψ = ∫ P_δψ(k) k² dk/(2π²) = {sigma2_psi_total:.6e}")
print(f"  √σ² = {math.sqrt(abs(sigma2_psi_total)):.6e}")

# With ψ̇ ~ m_tach × δψ:
# <ψ̇²> ~ m² × σ² → B_ψ = 3m²σ²
B_total = 3 * mu_sq * sigma2_psi_total
print(f"\n  Backreaction:")
print(f"    B_ψ/H₀² = 3 × m² × σ² = 3 × {mu_sq:.1f} × {sigma2_psi_total:.2e} = {B_total:.6e}")
print(f"    Required: 0.174")

# ΔH₀/H₀ from this backreaction:
if B_total > 0:
    Delta_H = math.sqrt(1 + B_total/3) - 1
    H0_pred = H0 * (1 + Delta_H)
    print(f"\n  H₀ prediction:")
    print(f"    H₀_eff = H₀ × √(1 + B/3) = {H0:.1f} × {math.sqrt(1 + B_total/3):.6f}")
    print(f"    H₀_eff = {H0_pred:.2f} km/s/Mpc")
    print(f"    ΔH₀/H₀ = {Delta_H*100:.4f}%")

# ==========================================================================
# SUMMARY
# ==========================================================================
print(f"\n{'='*78}")
print("  SUMMARY")
print("="*78)
print(f"""
  TACHYONIC AMPLIFICATION MECHANISM — QUANTITATIVE RESULTS:

  1. Tachyonic scale: λ_tach = {2*math.pi/mu_comov * 4451:.0f} Mpc (super-Hubble)
     All sub-Hubble modes are tachyonic (m²_eff < 0)

  2. Transfer function enhancement at pole: |T|² = {T_at_pole**2:.1f}×
     Regularized by Hubble damping Γ = 3H ≈ {Gamma:.1f}H₀

  3. Integrated σ²_ψ = {sigma2_psi_total:.2e}
     With tachyonic ψ̇ amplification: B_ψ/H₀² = {B_total:.2e}

  4. The effect is {B_total/0.174*100:.1f}% of what's needed for H₀ tension.

  ASSESSMENT:
  - If B_ψ/H₀² >> 0.17: TGP OVER-PREDICTS → needs additional damping
  - If B_ψ/H₀² ~ 0.17: PERFECT → TGP naturally explains H₀ tension
  - If B_ψ/H₀² << 0.17: TOO SMALL → need nonlinear amplification

  CAVEATS:
  - Used simplified matter power spectrum
  - Quasi-static approximation for transfer function
  - Regularization by Γ = 3H is approximate
  - Missing: nonlinear backreaction on the transfer function itself
  - Missing: scale-dependent growth from full TGP perturbation theory
""")
print(f"{'='*78}")
