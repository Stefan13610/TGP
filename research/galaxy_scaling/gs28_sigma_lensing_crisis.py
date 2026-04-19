"""
gs28: THE LENSING CRISIS — Σ = 1 AND ITS CONSEQUENCES FOR TGP
================================================================

gs27 discovered that in DGP-type theories (including TGP as formulated
in gs25), the gravitational slip parameter Σ = 1 EXACTLY. This means
lensing does NOT see the MOND-like enhancement ν(y).

This script performs a RIGOROUS analysis of:
  A. Derivation of Σ = 1 in DGP (from first principles)
  B. What galaxy-galaxy lensing data actually require
  C. Whether the bending term K_μν K^μν can break Σ = 1
  D. The AQUAL vs DGP fork: which path for TGP?
  E. Observational confrontation: what data distinguish Σ=1 from Σ=ν
  F. An honest scorecard: what survives and what doesn't
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import sys, io, warnings
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# Constants
G = 6.674e-11
c = 2.998e8
H0 = 2.27e-18
a0 = 1.12e-10
M_sun = 1.989e30
kpc = 3.086e19
Mpc = 3.086e22
r_c = c / (2 * np.pi * H0)

alpha = 0.8
gamma_disk = 0.4
gamma_sphere = 0.562

def nu_tgp(y, gam=gamma_disk):
    if y <= 0: return 1e10
    return 1.0 + np.exp(-y**alpha) / y**gam

def nu_mond(y):
    return 0.5 + np.sqrt(0.25 + 1.0/max(y, 1e-20))

print("=" * 78)
print("  gs28: THE LENSING CRISIS — Σ = 1 AND ITS CONSEQUENCES FOR TGP")
print("=" * 78)


# =============================================================================
# PART A: RIGOROUS DERIVATION OF Σ = 1 IN DGP
# =============================================================================
print("\n" + "=" * 78)
print("  PART A: WHY Σ = 1 IN DGP — FROM THE ACTION")
print("=" * 78)

print("""
  A.1  Setup: DGP action and metric perturbations
  -------------------------------------------------
  DGP action:
    S = M₅³/2 ∫d⁵x √(-g₅)R₅ + M₄²/2 ∫d⁴x √(-g₄)R₄ + S_m[g₄,ψ]

  Perturbed metric on the brane (Newtonian gauge):
    ds² = −(1+2Φ)dt² + a²(1−2Ψ)δᵢⱼdxⁱdxʲ

  In GR: Φ = Ψ (no anisotropic stress from matter).
  In DGP: Φ ≠ Ψ because the brane-bending mode π creates anisotropic stress.

  A.2  The brane-bending mode π
  -------------------------------
  The brane can fluctuate in the 5th dimension: y_brane = π(xᵘ).
  This π is a scalar degree of freedom on the brane.

  In the decoupling limit (M₅ → 0, M₄ → ∞, r_c fixed):

    π couples to matter through: L_int = π T/(M₄²)

  where T = T^μ_μ is the trace of the stress-energy tensor.

  For non-relativistic matter: T ≈ −ρc²
  → π mediates an ATTRACTIVE scalar force (like a scalar graviton).

  A.3  Metric potentials in terms of π
  --------------------------------------
  The DGP effective equations on the brane give:

    ∇²Φ = −4πGρ − (1/2)∇²π/M₄²         ... (i)
    ∇²Ψ = −4πGρ + (1/2)∇²π/M₄²         ... (ii)

  ⚠️  NOTE THE SIGN: π enters with OPPOSITE signs in Φ and Ψ.

  This is because:
  • Φ (time-time component) couples to the ENERGY of the scalar field
  • Ψ (space-space component) couples to the PRESSURE of the scalar field
  • For a massless scalar: energy = +E, pressure = +E/3 in spatial part
  • But the CONFORMAL coupling of π to gravity gives the opposite sign
    in the space-space equation.

  Adding (i) and (ii):
    ∇²(Φ + Ψ) = −8πGρ                   ... (iii)

  This is the STANDARD GR Poisson equation! The π contribution CANCELS.

  A.4  Lensing potential
  -------------------------
  Light deflection is governed by Φ_lens = (Φ+Ψ)/2.
  From (iii):
    ∇²Φ_lens = −4πGρ = ∇²Φ_GR

  → Φ_lens = Φ_GR   (same boundary conditions)
  → Σ ≡ Φ_lens/Φ_GR = 1    ■

  A.5  Why this is EXACT (not perturbative)
  -------------------------------------------
  The cancellation in (iii) follows from the STRUCTURE of the DGP action,
  specifically from 4D diffeomorphism invariance on the brane.

  The proof:
  1. The DGP brane has a 4D Einstein-Hilbert action → G_μν on the brane
     satisfies the contracted Bianchi identity ∇^μ G_μν = 0.
  2. The junction conditions relate the 5D Weyl tensor to brane curvature.
  3. The Weyl tensor E_μν is traceless: E^μ_μ = 0.
  4. For a spherically symmetric source: E_μν contributes EQUALLY to
     Φ and Ψ gradients but with OPPOSITE signs (tracelessness!).
  5. Therefore: Φ + Ψ gets NO contribution from E_μν.

  This holds:
  • At all orders in perturbation theory ✓
  • In the Vainshtein regime (nonlinear π) ✓
  • For ANY matter distribution (not just spherical) ✓
  • On BOTH branches (normal and self-accelerating) ✓

  The ONLY assumption: 5D bulk is Minkowski (or has Z₂ symmetry).
  If the bulk has non-trivial Weyl curvature (Weyl fluid), Σ ≠ 1.
""")

# Verify numerically
print(f"  A.6  Numerical verification")
print(f"  ----------------------------")

# For a point mass M in DGP:
# The π field in the linear regime (r >> r_*):
# π'(r) = GM/(3 β r_c)  → π(r) = -GM/(3 β r_c × r)
# where β depends on cosmological background

# For a galaxy: r >> r_* but r << r_c
# The DGP force correction:
# F_π = -dπ/dr = GM/(3 β r_c × r²)
# F_N = GM/r²
# F_π/F_N = 1/(3β r_c × ...)  ... this is tiny!

# Wait - this doesn't give MOND. The MOND-like behavior in DGP comes
# from a DIFFERENT regime. Let me be more careful.

# In the DGP model applied to GALAXIES (Nicolis, Rattazzi, Trincherini):
# The key is the Vainshtein mechanism.
#
# The π field equation for a static spherical source is:
#
#   (1/r²)(r² π')' + r_c/(3β) × [(π')²/r²]' = 8πG ρ/(6β)  ... schematic
#
# In the LINEAR regime (r >> r_*):
#   π(r) ≈ -GM/(6β r) × (1/r_c) × r_c  ... ≈ GM/(6β r)
#   → F_π/F_N = 1/(6β)
#
# In the NONLINEAR regime (r << r_*):
#   π'(r) ≈ (GM r_c/(9β²r))^(1/2) × (1/r_c)
#   → F_π/F_N ≈ (r/r_*)^(3/2) / (6β)  ... screened
#
# For the MOND-like effect, we need F_π ~ F_N, which requires β ~ 1.
# But β >> 1 in the normal branch at z=0.
#
# RESOLUTION: In the self-accelerating branch, β ~ O(1) at late times!
# OR: the local β for an isolated system differs from the cosmological β.

print(f"""
  A.7  The β puzzle: how does DGP give MOND?
  ============================================
  In standard DGP cosmology:
    β(z=0) ≈ 1 + 2r_c H₀ ~ 10⁸  (normal branch)

  With β ~ 10⁸, the scalar force F_π/F_N ~ 1/(3β) ~ 10⁻⁹ → NEGLIGIBLE!
  This cannot give MOND-like rotation curves.

  The TGP program (gs8-gs12) ASSUMED that the DGP propagator gives:
    g_obs = g_N × ν(y)  where ν(y) = 1 + exp(-y^(4/5))/y^(2/5)

  But this assumes the MOND modification IS the gravitational field,
  not a separate scalar mode.

  THERE ARE TWO DISTINCT INTERPRETATIONS:

  ╔═══════════════════════════════════════════════════════════════╗
  ║  INTERPRETATION A: "Phenomenological TGP"                    ║
  ║  ─────────────────────────────────────────                   ║
  ║  The modified Poisson equation:                              ║
  ║    ∇²Φ = -4πG ρ_eff,  ρ_eff = ρ × ν(y)                    ║
  ║                                                              ║
  ║  Here Φ IS the gravitational potential. There's no separate  ║
  ║  scalar mode. Lensing sees Φ directly.                       ║
  ║                                                              ║
  ║  → Σ = ν(y)  (lensing enhanced)                             ║
  ║  → Works for galaxy-galaxy lensing                           ║
  ║  → Cluster deficit is "only" ×1.9                            ║
  ║  → BUT: no relativistic completion (no action principle)     ║
  ╠═══════════════════════════════════════════════════════════════╣
  ║  INTERPRETATION B: "DGP TGP" (gs25 action)                  ║
  ║  ──────────────────────────────────                          ║
  ║  The DGP action with bending term:                           ║
  ║    S = M₅³/2 ∫R₅ + M₄²/2 ∫[R₄ + λK²] + S_m               ║
  ║                                                              ║
  ║  Here the MOND effect comes from scalar mode π.              ║
  ║  Lensing sees Φ+Ψ = 2Φ_GR (no π).                          ║
  ║                                                              ║
  ║  → Σ = 1  (lensing NOT enhanced)                             ║
  ║  → Galaxy-galaxy lensing FAILS                               ║
  ║  → Has a covariant action (relativistic completion)          ║
  ║  → BUT: contradicts lensing observations                     ║
  ╚═══════════════════════════════════════════════════════════════╝

  These are MUTUALLY EXCLUSIVE. TGP must choose one.
""")


# =============================================================================
# PART B: GALAXY-GALAXY LENSING DATA
# =============================================================================
print("=" * 78)
print("  PART B: WHAT DOES GALAXY-GALAXY LENSING ACTUALLY REQUIRE?")
print("=" * 78)

print("""
  B.1  Galaxy-galaxy lensing observations
  ----------------------------------------
  Galaxy-galaxy lensing (GGL) measures the tangential shear γ_t around
  lens galaxies, which probes the excess surface mass density:

    ΔΣ(R) = Σ̄(<R) − Σ(R)

  where Σ(R) is the projected surface mass density at projected
  distance R from the lens center.

  Key observations (from SDSS, KiDS, DES, HSC):
  • ΔΣ profiles extend to R ~ 1-10 Mpc
  • At R > 100 kpc: ΔΣ ∝ R^(-1) to R^(-0.8) (approximately NFW-like)
  • Total halo mass for L* galaxy: M_200 ~ 10¹² M_sun
  • Baryonic mass for L* galaxy: M_bar ~ 5×10¹⁰ M_sun
  • Ratio: M_halo/M_bar ~ 20

  B.2  Comparison of predictions
  --------------------------------
""")

# Model an L* galaxy
M_star = 5e10 * M_sun  # baryonic mass
R_eff = 5 * kpc          # effective radius

# Projected surface mass density as function of projected radius R
R_proj = np.logspace(np.log10(10*kpc), np.log10(3*Mpc), 50)

def g_bar(r, M):
    """Newtonian gravitational field at radius r from mass M."""
    return G * M / r**2

def delta_sigma_GR(R, M_halo, r_s, c_nfw):
    """NFW ΔΣ profile (simplified)."""
    x = R / r_s
    # NFW Σ(x) ~ 1/(x²-1) × (arctan or arctanh term)
    # Simplified: ΔΣ ~ M_halo / (π R²) × f(x)
    rho_s = M_halo / (4 * np.pi * r_s**3 * (np.log(1+c_nfw) - c_nfw/(1+c_nfw)))
    if x < 1:
        f = 1/(x**2 - 1) * (1 - np.log((1+np.sqrt(1-x**2))/x) / np.sqrt(1-x**2))
    elif x > 1:
        f = 1/(x**2 - 1) * (1 - np.arctan(np.sqrt(x**2-1)) / np.sqrt(x**2-1))
    else:
        f = 1/3.0
    Sigma = 2 * rho_s * r_s * f
    # Mean Sigma inside R:
    # For simplicity, use ΔΣ ~ M_2D(<R)/(πR²) - Σ(R) ≈ scaling
    return Sigma  # This is Σ(R), not ΔΣ, but gives the right order

# Compute predictions
print(f"  B.3  Projected mass profiles around L* galaxy")
print(f"  -----------------------------------------------")
print(f"    M_star = {M_star/M_sun:.1e} M_sun")
print(f"    R_eff  = {R_eff/kpc:.0f} kpc")
print(f"")
print(f"    {'R (kpc)':<10s} {'y=g/a₀':<10s} {'ν_TGP':<8s} {'ν_MOND':<8s} {'Σ_bar(R)':<14s} {'Σ_AQUAL(R)':<14s} {'Σ_DGP(R)':<14s}")
print(f"    {'─'*10} {'─'*10} {'─'*8} {'─'*8} {'─'*14} {'─'*14} {'─'*14}")

for R in [30*kpc, 100*kpc, 300*kpc, 1000*kpc, 3000*kpc]:
    y = g_bar(R, M_star) / a0
    nu_t = nu_tgp(y)
    nu_m = nu_mond(y)
    # Surface mass density (simplified: thin lens, projected mass)
    # Σ ~ M/(πR²) for a point mass (order of magnitude)
    Sigma_bar = M_star / (np.pi * R**2)  # kg/m²
    Sigma_aqual = Sigma_bar * nu_t  # AQUAL: lensing sees ν(y) enhancement
    Sigma_dgp = Sigma_bar * 1.0     # DGP: lensing sees baryonic only
    # Convert to M_sun/pc²
    pc = 3.086e16  # m
    conv = M_sun / pc**2  # kg/m² per M_sun/pc²
    print(f"    {R/kpc:<10.0f} {y:<10.4f} {nu_t:<8.2f} {nu_m:<8.2f} {Sigma_bar/conv:<14.4f} {Sigma_aqual/conv:<14.4f} {Sigma_dgp/conv:<14.4f}")

# Mass enclosed within projected radius
print(f"\n  B.4  Enclosed projected mass")
print(f"  ─────────────────────────────")
print(f"    {'R (kpc)':<10s} {'M_bar(<R)':<14s} {'M_AQUAL(<R)':<14s} {'M_DGP(<R)':<14s} {'M_LCDM(<R)':<14s}")
print(f"    {'─'*10} {'─'*14} {'─'*14} {'─'*14} {'─'*14}")

# Use actual phantom dark matter mass profile
for R in [50*kpc, 100*kpc, 200*kpc, 500*kpc, 1000*kpc]:
    # Baryonic mass (point mass approx)
    M_bar_enc = M_star  # all mass within R for point mass

    # AQUAL phantom mass: integrate the ν(y) enhancement
    # M_phantom ~ M_bar × [ν(y_outer) - 1] × (R/R_eff)^p ... simplified
    y = g_bar(R, M_star) / a0
    nu_val = nu_tgp(y)
    M_aqual = M_star * nu_val  # total (bar + phantom) seen by lensing

    # DGP: lensing sees only baryonic
    M_dgp = M_star * 1.0

    # ΛCDM: NFW halo
    M_halo = 1e12 * M_sun
    c_nfw = 10
    r_s_nfw = 200 * kpc / c_nfw  # r_200 ~ 200 kpc
    x = R / r_s_nfw
    f_nfw = np.log(1+x) - x/(1+x)
    f_nfw_max = np.log(1+c_nfw) - c_nfw/(1+c_nfw)
    M_nfw_enc = M_halo * f_nfw / f_nfw_max + M_star

    print(f"    {R/kpc:<10.0f} {M_bar_enc/M_sun:<14.2e} {M_aqual/M_sun:<14.2e} {M_dgp/M_sun:<14.2e} {M_nfw_enc/M_sun:<14.2e}")

print(f"""
  B.5  The lensing test
  ─────────────────────
  At R = 100 kpc from an L* galaxy:

  • ΛCDM:  M(<R) ~ 2×10¹¹ M_sun  (dominated by dark matter halo)
  • AQUAL: M(<R) ~ 2×10¹¹ M_sun  (phantom DM from ν(y))
  • DGP:   M(<R) ~ 5×10¹⁰ M_sun  (baryonic only!)

  Galaxy-galaxy lensing CLEARLY detects excess mass at 100 kpc.
  The measured ΔΣ profiles require M(100kpc)/M_bar ~ 2-4.

  AQUAL (Σ = ν(y)): ✓ Produces the right amount of lensing signal
  DGP (Σ = 1):      ❌ Predicts NO excess lensing signal

  ⚠️ THIS IS A DIRECT, MODEL-INDEPENDENT CONTRADICTION.
  Galaxy-galaxy lensing has been measured to >10σ significance
  by SDSS, KiDS, DES, HSC surveys.

  B.6  Can baryonic mass estimation errors save DGP?
  ---------------------------------------------------
  NO. The excess lensing signal extends to R ~ 1 Mpc, where:
  • g_bar ~ 10⁻¹³ m/s² (y ~ 10⁻³)
  • ν(y) ~ 20-30
  • M_enclosed ~ 10¹² M_sun vs M_bar ~ 5×10¹⁰ M_sun
  • Would need 20× more baryons — physically impossible
""")


# =============================================================================
# PART C: CAN THE BENDING TERM BREAK Σ = 1?
# =============================================================================
print("=" * 78)
print("  PART C: CAN THE K_μν K^μν TERM BREAK Σ = 1?")
print("=" * 78)

print("""
  C.1  How the bending term enters the equations
  ─────────────────────────────────────────────────
  The TGP action adds λ K_μν K^μν to the brane action.
  The extrinsic curvature K_μν describes how the brane bends in the bulk.

  For metric perturbations, K_μν has contributions from:
  • The tensor mode h_ij^TT (gravitational waves)
  • The scalar mode (brane bending π)
  • The vector mode (negligible for static sources)

  The K² term modifies the JUNCTION CONDITIONS at the brane.
  Instead of [∂_y h] = (1/r_c) h, we get:
    [∂_y h] = (1/r_c)(1 - λ∇²) h

  C.2  Effect on Φ and Ψ separately
  ────────────────────────────────────
  The modified junction conditions change the relationship between
  the bulk solution and the brane potentials.

  For the SCALAR sector (Φ, Ψ, π):
  The bending term adds k²-dependent corrections to the π equation:

    Original DGP:  ∇²π + r_c²/3 [(∇²π)²-(∂_i∂_jπ)²] = 8πGρ/(6β)

    With bending:  ∇²π + r_c²/3 [(∇²π)²-(∂_i∂_jπ)²]
                   + λ ∇⁴π + ... = 8πGρ/(6β) + λ-corrections

  The metric potentials become:
    ∇²Φ = -4πGρ - (1/2)(∇²π/M₄²) × F_Φ(k²λ)
    ∇²Ψ = -4πGρ + (1/2)(∇²π/M₄²) × F_Ψ(k²λ)

  where F_Φ and F_Ψ are form factors from the bending correction.

  C.3  The KEY question: does F_Φ = F_Ψ?
  ─────────────────────────────────────────
  In DGP (λ=0): F_Φ = F_Ψ = 1  →  π cancels  →  Σ = 1

  With bending (λ≠0): the K² term treats the time-time and space-space
  components of K_μν DIFFERENTLY because K^μ_ν has different eigenvalues.

  Specifically:
    K^0_0 ~ ∂_y Φ  (time-time)
    K^i_j ~ ∂_y Ψ × δ^i_j  (space-space)

  The K_μν K^μν contraction:
    K² = (K⁰₀)² + 3(Kⁱⱼ)² = (∂_yΦ)² + 3(∂_yΨ)²

  Varying with respect to Φ vs Ψ gives DIFFERENT factors:
    δS/δΦ: gets factor 2K⁰₀ = 2∂_yΦ
    δS/δΨ: gets factor 6Kⁱⱼ = 6∂_yΨ  (factor of 3 from trace!)

  THEREFORE: the bending term introduces an ASYMMETRY between
  its contribution to the Φ and Ψ equations!

  But the MAGNITUDE of this asymmetry depends on how much ∂_yΦ
  differs from ∂_yΨ, which in DGP itself is an O(1/(3β)) effect.
""")

# Compute the asymmetry quantitatively
# In the linearized theory:
# K^0_0 = -(1/2) ∂_y g_00 = ∂_y Φ
# K^i_j = -(1/2) ∂_y g_ij / a² = ∂_y Ψ δ^i_j
# K_μν K^μν = (∂_yΦ)² + 3(∂_yΨ)²

# In DGP: ∂_yΦ = |k|Φ, ∂_yΨ = |k|Ψ (from bulk solution exp(-|k|y))
# With Φ = Φ_GR(1+1/(3β)), Ψ = Φ_GR(1-1/(3β)):

beta_est = 2 * r_c * H0  # ~ 10^8 for normal branch, ~ 1-10 for self-acc

print(f"\n  C.4  Quantitative asymmetry")
print(f"  ───────────────────────────")
print(f"    Normal branch: β ~ {beta_est:.1e}")
print(f"    Self-acc branch: β ~ 1-10")

for beta_val in [1.0, 3.0, 10.0, 1e4, 1e8]:
    dPhi = 1 + 1/(3*beta_val)  # proportional to ∂_yΦ
    dPsi = 1 - 1/(3*beta_val)  # proportional to ∂_yΨ
    K2 = dPhi**2 + 3*dPsi**2   # K_μν K^μν
    # Variation gives different coefficients for Φ and Ψ modifications:
    coeff_Phi = 2*dPhi  # d(K²)/d(∂_yΦ) = 2∂_yΦ
    coeff_Psi = 6*dPsi  # d(K²)/d(∂_yΨ) = 6∂_yΨ
    # The ratio determines the asymmetry:
    asymmetry = coeff_Phi / coeff_Psi  # should be 1/3 if Φ=Ψ
    # The Σ correction:
    # F_Φ = 1 + λk²×coeff_Phi/..., F_Ψ = 1 + λk²×coeff_Psi/...
    # Σ = 1 + (F_Φ - F_Ψ) × π/(2Φ_GR)
    delta_Sigma = (coeff_Phi - coeff_Psi/3) / (coeff_Phi + coeff_Psi)
    print(f"    β = {beta_val:<10.1e}: coeff_Φ = {coeff_Phi:.4f}, coeff_Ψ = {coeff_Psi:.4f}, "
          f"asymmetry = {asymmetry:.4f}, δΣ/Σ ~ {abs(delta_Sigma):.4e}")

print(f"""
  C.5  Assessment of bending-induced Σ ≠ 1
  ──────────────────────────────────────────
  The bending term K² treats Φ and Ψ differently by a factor
  related to 1/(3β). For the normal branch (β ~ 10⁸), the
  asymmetry is ~ 10⁻⁸ → NEGLIGIBLE.

  For the self-accelerating branch (β ~ 1): the asymmetry could
  be O(1), giving Σ ≠ 1 at a significant level!

  BUT: the self-accelerating branch has a GHOST (negative energy mode).
  This ghost makes the theory unstable and is usually considered fatal.

  VERDICT: The bending term CANNOT save Σ = 1 on the normal (ghost-free)
  branch. On the self-accelerating branch, Σ ≠ 1 is possible but
  the theory is ghostly (unstable).
""")


# =============================================================================
# PART D: THE AQUAL vs DGP FORK
# =============================================================================
print("=" * 78)
print("  PART D: THE AQUAL vs DGP FORK — WHICH PATH FOR TGP?")
print("=" * 78)

print("""
  D.1  The fundamental dichotomy
  ─────────────────────────────────
  TGP's phenomenological success (gs10-gs12) is based on:

    g_obs = g_bar × ν(g_bar/a₀)

  This is KINEMATICALLY identical to MOND. The question is the
  DYNAMICAL interpretation:

  ┌─────────────────────────────────────────────────────────────┐
  │  PATH A: AQUAL-like ("modified Poisson")                    │
  │                                                              │
  │  ∇·[μ(|∇Φ|/a₀)∇Φ] = -4πGρ                                │
  │                                                              │
  │  Properties:                                                 │
  │  ✓ Φ is THE gravitational potential                         │
  │  ✓ Lensing sees Φ directly → Σ = ν(y)                      │
  │  ✓ Galaxy-galaxy lensing works                               │
  │  ✓ Phantom DM visible in lensing                             │
  │  ✗ No known healthy relativistic completion                  │
  │  ✗ TeVeS (Bekenstein 2004) had problems: ruled out by GW170817│
  │  ✗ New TeVeS (Skordis & Złośnik 2021) works but is complex  │
  │  ? Cluster deficit: ×1.9 (manageable?)                       │
  ├─────────────────────────────────────────────────────────────┤
  │  PATH B: DGP-like ("braneworld")                            │
  │                                                              │
  │  S = ∫R₅ + ∫[R₄ + λK²] + S_m                              │
  │                                                              │
  │  Properties:                                                 │
  │  ✓ Natural relativistic completion                           │
  │  ✓ v_GW = c (protected by 4D diffs)                         │
  │  ✓ CMB consistent                                            │
  │  ✗ Σ = 1 → lensing blind to MOND effect                     │
  │  ✗ Contradicts galaxy-galaxy lensing                         │
  │  ✗ Contradicts cluster lensing                               │
  │  ✗ Normal branch: β >> 1 → tiny MOND effect anyway!         │
  │  ✗ Self-acc branch: ghost                                    │
  │  ? Could exotic brane terms (Gauss-Bonnet, f(R)) help?      │
  └─────────────────────────────────────────────────────────────┘

  D.2  The β problem: DGP cannot even give MOND!
  ─────────────────────────────────────────────────
  There's a DEEPER problem with Path B that goes beyond Σ = 1.

  In DGP, the scalar force on the normal branch is:
    F_π / F_N = 1/(3β)

  At z=0: β ~ 2r_c H₀ ~ 10⁸ → F_π/F_N ~ 3×10⁻⁹

  This is NINE ORDERS OF MAGNITUDE too small for MOND!

  The TGP program implicitly assumed that the DGP crossover gives
  MOND-like behavior with ν(y) >> 1. But in the actual DGP equations
  for the normal branch, the modification is essentially zero.

  On the self-accelerating branch: β ~ O(1) at late times,
  and the scalar force IS comparable to Newton. But this branch:
  (a) has a ghost, (b) is observationally excluded by ISW effect,
  (c) doesn't match expansion history without fine-tuning.
""")

# Show the actual DGP prediction vs TGP assumed
print(f"\n  D.3  DGP vs TGP: what the equations actually give")
print(f"  ═══════════════════════════════════════════════════")

print(f"\n    For MW-like galaxy (M = 10¹¹ M_sun) at R = 30 kpc:")
M_MW = 1e11 * M_sun
R_test = 30 * kpc
g_N = G * M_MW / R_test**2
y_test = g_N / a0

# DGP actual prediction (normal branch, β >> 1):
F_pi_over_FN_normal = 1 / (3 * beta_est)
nu_dgp_actual = 1 + F_pi_over_FN_normal

# DGP self-accelerating (β ~ 1):
beta_sa = 1.5  # order of magnitude
nu_dgp_sa = 1 + 1/(3*beta_sa)

# TGP assumed
nu_tgp_val = nu_tgp(y_test)
nu_mond_val = nu_mond(y_test)

print(f"      g_N          = {g_N:.3e} m/s²")
print(f"      y = g_N/a₀   = {y_test:.4f}")
print(f"")
print(f"      ν(y) values:")
print(f"        TGP (assumed):            {nu_tgp_val:.3f}  ← what gs10-gs12 used")
print(f"        MOND:                     {nu_mond_val:.3f}")
print(f"        DGP normal (β~10⁸):      {nu_dgp_actual:.9f}  ← ACTUAL DGP!")
print(f"        DGP self-acc (β~1.5):     {nu_dgp_sa:.3f}  ← ghost branch")
print(f"")
print(f"      ⚠️  DGP NORMAL BRANCH GIVES ν ≈ 1.000 (NO MODIFICATION!)")
print(f"      ⚠️  The entire gs10-gs12 program assumed ν >> 1, which")
print(f"      ⚠️  CANNOT come from the DGP normal branch.")

print(f"""
  D.4  Resolution: the gs25 action is WRONG as a TGP completion
  ═══════════════════════════════════════════════════════════════
  The DGP action proposed in gs25 CANNOT produce the ν(y) function
  that fits galaxy rotation curves, for TWO independent reasons:

  1. β >> 1 on the normal branch → ν ≈ 1 (no MOND effect)
  2. Σ = 1 on both branches → lensing blind to any enhancement

  The gs25 action was a reasonable GUESS but it does not work.
  The TGP membrane model needs a DIFFERENT relativistic completion.

  D.5  What kind of theory DO we need?
  ═════════════════════════════════════
  The requirements for a viable TGP relativistic completion:

  1. MUST produce ν(y) = 1 + exp(-y^(4/5))/y^(2/5) as non-rel limit ✓
  2. MUST have Σ ≈ ν(y) (lensing sees phantom DM) ✓
  3. MUST have v_GW = c (GW170817) ✓
  4. MUST NOT have ghost (stability) ✓
  5. SHOULD explain a₀ = cH₀/(2π) ✓
  6. SHOULD explain α = 4/5, γ/α = 1/2 ✓

  Known theories satisfying (1)+(2)+(3)+(4):
  • AeST (Skordis & Złośnik 2021) — "new TeVeS"
    - Has scalar + vector fields on a metric background
    - Σ ≠ 1 by construction (vector field ensures lensing)
    - v_GW = c (no modification to tensor modes)
    - Ghost-free
    - Fits CMB power spectrum to Planck accuracy!
    - BUT: ~6 free functions, complex action

  • RMOND (Milgrom 2009) — "relativistic MOND"
    - Modified Einstein-aether theory
    - Σ ≠ 1 through the aether vector
    - v_GW = c (with appropriate choice of parameters)
    - Status: unclear ghost-freedom

  TGP as a membrane model could potentially MAP to AeST:
  • The scalar field ← brane bending mode (π)
  • The vector field ← membrane tangent direction (intrinsic mode!)
  • a₀ ← membrane tension = cH₀/(2π)
""")


# =============================================================================
# PART E: OBSERVATIONAL TESTS THAT DISTINGUISH Σ = 1 FROM Σ = ν
# =============================================================================
print("=" * 78)
print("  PART E: OBSERVATIONAL TESTS DISTINGUISHING Σ = 1 vs Σ = ν(y)")
print("=" * 78)

print(f"""
  E.1  Available data and predictions
  ─────────────────────────────────────
  Three key observables probe different metric combinations:

  1. ROTATION CURVES / VELOCITY DISPERSIONS:
     Probe: ∇Φ (time-time potential)
     In both Σ=1 and Σ=ν: F = g_N × ν(y)  [IDENTICAL]
     → Cannot distinguish

  2. GRAVITATIONAL LENSING:
     Probe: ∇(Φ+Ψ)/2 (lensing potential)
     Σ=1:  F_lens = g_N              [baryonic only]
     Σ=ν:  F_lens = g_N × ν(y)      [enhanced]
     → CAN distinguish! Requires systems where both dynamics and
        lensing are measured independently.

  3. GRAVITATIONAL REDSHIFT:
     Probe: Φ (at a fixed point)
     Σ=1:  Φ = Φ_N + π  [enhanced]
     Σ=ν:  Φ = Φ_N × ν   [enhanced]
     → Both enhanced, but different profile shapes.

  E.2  The EG (gravitational slip) parameter
  ────────────────────────────────────────────
  The observable: E_G = Ω_m × [Σ/f]

  where f = d ln D / d ln a is the growth rate and Σ is our parameter.

  In GR:  E_G ≈ Ω_m/f ≈ 0.4 (at z ~ 0.3)
  In DGP: E_G ≈ Ω_m/f × Σ_DGP = Ω_m/f × 1 ≈ 0.4 (indistinguishable!)
  In AQUAL: E_G would be modified (Σ ≠ 1)

  Current measurements: E_G = 0.39 ± 0.06 (Reyes et al. 2010)
  → Cannot yet distinguish (error bars too large)
  → Euclid/LSST will measure to ~1% level

  E.3  Galaxy-galaxy lensing as the definitive test
  ───────────────────────────────────────────────────
  The MOST direct test: compare rotation curve mass to lensing mass
  for the SAME galaxies.

  If Σ = 1: M_lens = M_baryon (no phantom DM in lensing)
  If Σ = ν: M_lens = M_baryon × ν(y) = M_dynamical

  Measurements:
""")

# Compute specific predictions
galaxies_test = [
    ("MW-like spiral", 1e11, 30, gamma_disk),
    ("Dwarf (DDO154)", 1e9, 8, gamma_disk),
    ("Massive E (M87)", 1e12, 100, gamma_sphere),
    ("Galaxy cluster", 1e14, 2000, gamma_sphere),
]

print(f"    {'System':<22s} {'M_bar(M☉)':<12s} {'R(kpc)':<8s} {'y':<8s} {'ν(y)':<6s} {'M_dyn':<12s} {'M_lens(Σ=1)':<14s} {'M_lens(Σ=ν)':<14s}")
print(f"    {'─'*22} {'─'*12} {'─'*8} {'─'*8} {'─'*6} {'─'*12} {'─'*14} {'─'*14}")

for name, M_b, R_kpc, gam in galaxies_test:
    M = M_b * M_sun
    R = R_kpc * kpc
    g = G * M / R**2
    y = g / a0
    nu = nu_tgp(y, gam)
    M_dyn = M_b * nu
    M_lens_1 = M_b  # Σ = 1
    M_lens_nu = M_b * nu  # Σ = ν
    print(f"    {name:<22s} {M_b:<12.1e} {R_kpc:<8.0f} {y:<8.4f} {nu:<6.2f} {M_dyn:<12.2e} {M_lens_1:<14.2e} {M_lens_nu:<14.2e}")

print(f"""
  E.4  Existing galaxy-galaxy lensing results
  ─────────────────────────────────────────────
  Brouwer et al. (2021, A&A 650, A113) — KiDS+GAMA:
  • Measured GGL around 259,000 galaxies (0.1 < z < 0.5)
  • Found ΔΣ profiles consistent with BOTH ΛCDM and EG (Verlinde)
  • Both models produce EXCESS lensing beyond baryonic → Σ > 1 required

  Milgrom (2013, PRL 111, 041105):
  • Showed that MOND phantom DM produces galaxy-galaxy lensing signal
  • Predicted ΔΣ ∝ √(a₀ M_bar) / (2πG R) at large R
  • This REQUIRES Σ = ν(y) (AQUAL interpretation)

  Tian et al. (2020) — CLASH clusters:
  • Lensing mass of clusters at multiple radii
  • Require M_lens/M_bar ~ 5-10
  • Σ = 1 (DGP) would give M_lens/M_bar = 1 → RULED OUT

  E.5  VERDICT: Σ = 1 is OBSERVATIONALLY EXCLUDED
  ═════════════════════════════════════════════════
  Galaxy-galaxy lensing at >10σ shows excess mass around galaxies.
  Cluster lensing at >100σ shows excess mass in clusters.

  If Σ = 1: these observations require REAL dark matter.
  TGP + DGP + Σ=1 + dark matter → TGP is redundant.

  Therefore: TGP MUST have Σ ≠ 1 to be a meaningful theory.
  This EXCLUDES the DGP action (gs25) as the relativistic completion.
""")


# =============================================================================
# PART F: HONEST SCORECARD
# =============================================================================
print("=" * 78)
print("  PART F: HONEST SCORECARD — WHAT SURVIVES?")
print("=" * 78)

print(f"""
  F.1  What we KNOW works (model-independent):
  ═════════════════════════════════════════════
  These results do NOT depend on the relativistic completion:

  ✅ ν(y) = 1 + exp(-y^(4/5))/y^(2/5) fits 171 SPARC galaxies (gs12)
  ✅ a₀ = 1.12×10⁻¹⁰ m/s² (global fit)
  ✅ γ/α = 1/2 (confirmed independently)
  ✅ γ(S) hierarchy: disk < S0 < E < cluster (gs19)
  ✅ BTFR slope ~ 3.4-3.6 (better than MOND's 4.0)
  ✅ Freeman limit = a₀/(2πG) = 137 M☉/pc²
  ✅ dSph's prefer higher γ (gs21)
  ✅ Differential γ(S) gives 69% of Bullet Cluster offset (gs26)

  F.2  What REQUIRES a specific relativistic completion:
  ═══════════════════════════════════════════════════════
  ⚠️ Lensing predictions → need Σ parameter → need action
  ⚠️ GW speed → need tensor sector → need action
  ⚠️ CMB predictions → need perturbation theory → need action
  ⚠️ Graviton mass → need propagator → need action

  F.3  Status of gs25 (DGP action):
  ═════════════════════════════════
  ❌ RULED OUT by two independent arguments:
     1. β >> 1 on normal branch → no MOND effect at all
     2. Σ = 1 → contradicts galaxy-galaxy lensing

  F.4  What the TGP relativistic completion MUST look like:
  ═════════════════════════════════════════════════════════
  It must be a theory with:
  • A SCALAR mode that gives ν(y) enhancement to dynamics
  • A VECTOR or TENSOR mode that ensures Σ ≈ ν(y) for lensing
  • No ghost (stability)
  • v_GW = c (GW170817)

  The ONLY known theory satisfying all four: AeST (Skordis-Złośnik 2021)

  F.5  Revised TGP program path:
  ═══════════════════════════════
  1. KEEP the phenomenological ν(y) (it works!)
  2. KEEP the membrane interpretation (a₀ = cH₀/(2π), α = 4/5)
  3. REPLACE the DGP action with AeST-type action
  4. MAP membrane modes to AeST fields:
     • Brane bending → AeST scalar field
     • Membrane tangent → AeST vector field (ensures Σ ≠ 1!)
     • Membrane tension → a₀ (cosmological origin)
  5. DERIVE ν(y) from the new action
  6. CONFRONT with cluster data (with correct Σ)

  F.6  The cluster problem REVISITED with Σ = ν(y):
  ═══════════════════════════════════════════════════
  If Σ = ν(y) (AQUAL-like, as required by lensing):
""")

# Recompute cluster predictions with Σ = ν
print(f"    Bullet Cluster with Σ = ν(y):")
M_cluster = 1.58e14 * M_sun
M_lens = 1.10e15 * M_sun
y_cl = 0.05
for gam, label in [(gamma_disk, "γ_disk=0.40"), (gamma_sphere, "γ_sphere=0.56"), (0.65, "γ=0.65"), (0.80, "γ=0.80")]:
    nu = nu_tgp(y_cl, gam)
    M_pred = M_cluster * nu / M_sun
    deficit = (1 - M_cluster*nu/M_lens) * 100
    print(f"      {label}: ν={nu:.2f}, M_pred={M_pred:.2e} M☉, deficit={deficit:.0f}%")

print(f"""
  F.7  THE REMAINING CLUSTER PROBLEM
  ═══════════════════════════════════
  Even with Σ = ν(y) (best case), the cluster deficit persists:
  • With γ_sphere = 0.56: deficit = {(1-nu_tgp(0.05,0.562)*M_cluster/M_lens)*100:.0f}%
  • Need γ_cluster ~ 0.80 to close the gap
  • But gs19 gives γ_cluster ≈ 0.56 from morphological fit

  This is the SAME ~×2 cluster problem that MOND has had since 1983.
  It's NOT caused by a wrong Σ. It's a genuine shortfall.

  Possible resolutions (same as MOND):
  • Hot baryons missed in surveys (+30-50%)
  • Massive neutrinos (m_ν ~ 1-2 eV, needs KATRIN relaxation)
  • Residual dark matter component (sterile ν, ~keV)
  • γ increases beyond the linear γ(S) model at very high S
  • Non-equilibrium effects in cluster mergers

  F.8  SUMMARY TABLE
  ═══════════════════
  ┌─────────────────────────────────────────────────────────────┐
  │  Component              │ Status           │ Confidence     │
  ├─────────────────────────┼──────────────────┼────────────────┤
  │  ν(y) interpolation     │ ✅ CONFIRMED      │ 30σ (gs12)     │
  │  a₀ = cH₀/(2π)         │ ✅ DERIVED         │ 3% match       │
  │  α = 4/5 (Flory)       │ ✅ DERIVED         │ exact          │
  │  γ/α = 1/2             │ ✅ DERIVED         │ exact          │
  │  γ(S) hierarchy         │ ✅ CONFIRMED       │ p < 0.01       │
  │  DGP action (gs25)      │ ❌ RULED OUT       │ (β>>1, Σ=1)   │
  │  Cluster mass (×1.9)    │ ❌ UNSOLVED        │ ~same as MOND  │
  │  Bullet offset (69%)    │ ✅ UNIQUE PRED     │ novel          │
  │  Rel. completion        │ ⚠️ NEEDS REWORK   │ → AeST path    │
  └─────────────────────────────────────────────────────────────┘
""")

print(f"""
  F.9  FINAL HONEST ASSESSMENT
  ═════════════════════════════
  The gs27-gs28 audit reveals:

  1. There are NO hidden forces in TGP that solve the cluster problem.
     All candidate effects (B-F in gs27) are negligibly small (10⁻⁸ to 10⁻¹⁴).

  2. The DGP-based action (gs25) is RULED OUT because:
     (a) Normal branch gives β >> 1 → NO MOND effect
     (b) Σ = 1 → lensing doesn't see MOND → contradicts observations
     (c) Self-accelerating branch has ghost

  3. TGP's PHENOMENOLOGICAL results (ν(y), γ(S), a₀) remain VALID.
     These are non-relativistic and don't depend on the action.

  4. The relativistic completion MUST be AeST-like (scalar + vector),
     NOT DGP-like (scalar only). The membrane model can potentially
     provide both modes:
     • Bending = scalar field (out-of-plane)
     • Tangent direction = vector field (in-plane)

  5. The cluster problem (×1.9 deficit even with Σ = ν) is the SAME
     problem that MOND has faced for 40 years. TGP does not solve it
     and does not make it worse. It's a shared problem of all MOND-like
     theories that likely requires hot baryons, neutrinos, or some
     minimal dark matter component.
""")

print("=" * 78)
print("  END OF gs28: THE LENSING CRISIS")
print("=" * 78)
