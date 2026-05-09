"""
Stage 1 — FRW tracker γ(z) simulation dla kandydata L

Test hipotezy: γ_eff(z) = γ_0 · [α + (1-α)·(H(z)/H_0)²]
  α=1.0  → canonical TGP_v1 (T-Λ as is)
  α<1.0  → mixed tracker (kandydat L)
  α=0.0  → full tracker quintessence

Cel: czy istnieje α ∈ [0.5, 0.99] dający:
  (a) ΔH/H_0 ∈ [5%, 10%]  (target H₀ tension solution)
  (b) ρ_Λ(z=1100)/ρ_total(z=1100) < 0.05  (CMB safety)
  (c) Δr_s/r_s acceptable
  (d) w_eff(z=0) ≥ -1

Implementuje ansatz (A) z Stage0_results §6.1.

Self-consistent FRW solver: H(z) ↔ γ(z) ↔ ρ_Λ(z) iteruje aż do konwergencji.

Output: stage1_tracker_gamma_FRW.txt z tabelą wyników per α.

Reference values (Planck 2018):
  H_0 = 67.36 km/s/Mpc
  Ω_m = 0.315, Ω_r = 9.18e-5, Ω_Λ = 0.6847
  z_recomb = 1090
  r_s_LCDM (z_recomb) = 144.43 Mpc

Reference SH0ES: H_0 = 73.04 ⇒ tension Δ = 8.4%
"""

import numpy as np
from scipy.integrate import quad

# ============================================================
# Cosmological parameters (Planck 2018)
# ============================================================
H0_PLANCK = 67.36       # km/s/Mpc
H0_SH0ES = 73.04        # km/s/Mpc
TENSION = (H0_SH0ES - H0_PLANCK) / H0_PLANCK  # ≈ 0.0843

OMEGA_M = 0.315
OMEGA_R = 9.18e-5       # radiation today (CMB+neutrinos)
OMEGA_L = 0.6847        # baseline LCDM
Z_RECOMB = 1090.0
Z_EQ = 3402.0           # matter-radiation equality

# Algebraic prediction TGP (closure_2026-04-26 + γ.1 + δ.1 + δ.2)
# Ω_Λ_TGP = 5e²/54 ≈ 0.6842 (independent of H_0 vs H(z) choice)
E_EULER = np.e
G_TILDE = 5 * E_EULER**2 / (12 * np.pi)  # ≈ 0.98003
OMEGA_L_TGP_ALGEBRAIC = (2 * np.pi / 9) * G_TILDE  # ≈ 0.6842
PHI_EFF = 8 * np.pi * G_TILDE  # ≈ 24.6302

# ============================================================
# H(z) functions
# ============================================================

def H_LCDM(z):
    """Standard LCDM H(z)/H_0 ratio (z OΩ_Λ = const)."""
    return np.sqrt(OMEGA_M * (1+z)**3 + OMEGA_R * (1+z)**4 + OMEGA_L)


def rho_total_over_rhocrit_today(z):
    """ρ_total(z) / ρ_crit(0) — total matter+radiation density z LCDM."""
    return OMEGA_M * (1+z)**3 + OMEGA_R * (1+z)**4


def H_TGP_self_consistent(z, alpha, n_iter=50, tol=1e-8):
    """
    Self-consistent solver dla H(z) z γ(z) tracker.

    Ansatz (A): γ(z) = γ_0 · [α + (1-α)·(H(z)/H_0)²]
    ⇒ ρ_Λ_TGP(z) / ρ_Λ_TGP(0) = α + (1-α)·(H(z)/H_0)²

    Friedmann (TGP cosmological background):
    (H(z)/H_0)² = Ω_m·(1+z)³ + Ω_r·(1+z)⁴ + Ω_Λ_TGP · [α + (1-α)·(H(z)/H_0)²]

    Define x = (H(z)/H_0)², matter+radiation = ρ_mr(z), then:
    x = ρ_mr(z) + Ω_Λ_TGP · [α + (1-α)·x]
    x · [1 - Ω_Λ_TGP·(1-α)] = ρ_mr(z) + Ω_Λ_TGP·α
    x = [ρ_mr(z) + Ω_Λ_TGP·α] / [1 - Ω_Λ_TGP·(1-α)]

    Solid analytic solution! No iteration needed.
    """
    Omega_L_eff = OMEGA_L_TGP_ALGEBRAIC  # use algebraic prediction

    rho_mr = rho_total_over_rhocrit_today(z)
    denom = 1.0 - Omega_L_eff * (1 - alpha)

    if denom <= 0:
        # α tak małe, że tracker dominuje i powstaje fixed point divergence
        return np.full_like(np.asarray(z, dtype=float), np.nan)

    x = (rho_mr + Omega_L_eff * alpha) / denom
    return np.sqrt(x)


def rho_Lambda_TGP_over_rhocrit_today(z, alpha):
    """ρ_Λ_TGP(z) / ρ_crit(0) z trackerem."""
    Omega_L_eff = OMEGA_L_TGP_ALGEBRAIC
    H_ratio_sq = H_TGP_self_consistent(z, alpha)**2
    return Omega_L_eff * (alpha + (1 - alpha) * H_ratio_sq)


# ============================================================
# Sound horizon r_s
# ============================================================

def c_s_over_c(z):
    """Adiabatic sound speed of photon-baryon fluid."""
    OMEGA_B_H2 = 0.02237  # Planck 2018
    OMEGA_G_H2 = 2.473e-5
    R = (3 * OMEGA_B_H2) / (4 * OMEGA_G_H2 * (1 + z))
    return 1.0 / np.sqrt(3 * (1 + R))


def r_s_integrand_LCDM(z):
    """Sound horizon integrand: c_s(z) / [(1+z) · H(z)]."""
    return c_s_over_c(z) / ((1 + z) * H_LCDM(z))


def r_s_integrand_TGP(z, alpha):
    H_ratio = H_TGP_self_consistent(z, alpha)
    if np.isnan(H_ratio):
        return np.nan
    return c_s_over_c(z) / ((1 + z) * H_ratio)


def sound_horizon_LCDM():
    """r_s(z_recomb) for LCDM, in units of c/H_0."""
    val, _ = quad(r_s_integrand_LCDM, Z_RECOMB, 1e8, limit=500)
    return val  # multiply by c/H_0 for Mpc; we keep dimensionless


def sound_horizon_TGP(alpha):
    val, _ = quad(lambda z: r_s_integrand_TGP(z, alpha),
                  Z_RECOMB, 1e8, limit=500)
    return val


# ============================================================
# Equation of state w(z) for tracker
# ============================================================

def w_DE_TGP(z, alpha, dz=0.001):
    """
    Effective equation of state w_DE(z) for tracker:
    w_DE(z) = -1 - (1/3) · d ln ρ_Λ(z) / d ln(1+z)

    For ρ_Λ(z) ∝ [α + (1-α)·H(z)²/H_0²]:
    d ln ρ_Λ / d ln(1+z) — numerical derivative
    """
    rho1 = rho_Lambda_TGP_over_rhocrit_today(z + dz, alpha)
    rho0 = rho_Lambda_TGP_over_rhocrit_today(max(z - dz, 0.001), alpha)
    dlnrho_dz = np.log(rho1/rho0) / (np.log(1 + z + dz) - np.log(1 + max(z - dz, 0.001)))
    return -1 - dlnrho_dz / 3.0


# ============================================================
# Main analysis per alpha
# ============================================================

def analyze_alpha(alpha):
    results = {'alpha': alpha}

    # CMB-era safety: ρ_Λ(z=1100) / ρ_total(z=1100)
    rho_L_recomb = rho_Lambda_TGP_over_rhocrit_today(Z_RECOMB, alpha)
    rho_total_recomb = rho_total_over_rhocrit_today(Z_RECOMB) + rho_L_recomb
    results['rho_L_over_rho_total_recomb'] = rho_L_recomb / rho_total_recomb

    # Equality-era: ρ_Λ(z=3402) / ρ_total(z=3402)
    rho_L_eq = rho_Lambda_TGP_over_rhocrit_today(Z_EQ, alpha)
    rho_total_eq = rho_total_over_rhocrit_today(Z_EQ) + rho_L_eq
    results['rho_L_over_rho_total_eq'] = rho_L_eq / rho_total_eq

    # Sound horizon shift Δr_s / r_s
    rs_LCDM = sound_horizon_LCDM()
    try:
        rs_TGP = sound_horizon_TGP(alpha)
        delta_rs = (rs_TGP - rs_LCDM) / rs_LCDM
    except Exception:
        delta_rs = np.nan
    results['delta_rs_over_rs'] = delta_rs

    # Inferred H_0 shift: ΔH_0_inf / H_0 ≈ -Δr_s / r_s
    # (CMB measures θ_s = r_s/D_A; if r_s smaller, inferred H_0 larger to match D_A)
    results['delta_H0_inferred_over_H0'] = -delta_rs

    # Coverage of H_0 tension (8.4% target)
    results['tension_coverage'] = -delta_rs / TENSION

    # w(z=0) — equation of state today
    results['w_today'] = w_DE_TGP(0.001, alpha)

    # Δw at z=0.5 (DESI relevant)
    results['w_z05'] = w_DE_TGP(0.5, alpha)

    # H(z=1100) shift vs LCDM
    H_recomb_LCDM = H_LCDM(Z_RECOMB)
    H_recomb_TGP = H_TGP_self_consistent(Z_RECOMB, alpha)
    results['H_recomb_TGP_over_LCDM'] = H_recomb_TGP / H_recomb_LCDM

    # H(z=0) check (must be ≈ 1 by construction)
    results['H_today_TGP_over_LCDM'] = H_TGP_self_consistent(0.001, alpha) / H_LCDM(0.001)

    return results


# ============================================================
# Run for all alphas
# ============================================================

ALPHAS = [1.0, 0.99, 0.95, 0.9, 0.75, 0.5, 0.25]

print("=" * 80)
print("Stage 1 — FRW tracker γ(z) simulation")
print("op-void-flat-modes-h0-2026-05-06")
print("=" * 80)
print()
print("Algebraic predictions (γ.1 + δ.1 + δ.2):")
print(f"  g̃ = 5e²/(12π) = {G_TILDE:.6f}")
print(f"  Φ_eff = 8π·g̃ = {PHI_EFF:.4f}")
print(f"  Ω_Λ_TGP = (2π/9)·g̃ = {OMEGA_L_TGP_ALGEBRAIC:.6f}  vs Planck 0.6847")
print()
print(f"H₀ tension target: ΔH/H = {TENSION*100:.2f}%")
print(f"CMB safety constraint: ρ_Λ(z=1100)/ρ_total(z=1100) < 0.05")
print()
print("Test 5+ values of α in mixing γ(z) = γ_0·[α + (1-α)·(H(z)/H_0)²]")
print("=" * 80)

print()
print(f"{'α':>5} | {'ρ_Λ(rec)/ρ_tot':>14} | {'ρ_Λ(eq)/ρ_tot':>14} | "
      f"{'Δr_s/r_s':>10} | {'ΔH₀_inf/H₀':>11} | {'Coverage':>9} | "
      f"{'w(z=0)':>8} | {'w(z=0.5)':>9} | {'H_rec_TGP/LCDM':>14}")
print("-" * 130)

all_results = []
for alpha in ALPHAS:
    r = analyze_alpha(alpha)
    all_results.append(r)
    print(f"{r['alpha']:>5.3f} | "
          f"{r['rho_L_over_rho_total_recomb']:>14.5e} | "
          f"{r['rho_L_over_rho_total_eq']:>14.5e} | "
          f"{r['delta_rs_over_rs']*100:>9.3f}% | "
          f"{r['delta_H0_inferred_over_H0']*100:>10.3f}% | "
          f"{r['tension_coverage']*100:>8.2f}% | "
          f"{r['w_today']:>8.4f} | "
          f"{r['w_z05']:>9.4f} | "
          f"{r['H_recomb_TGP_over_LCDM']:>14.6f}")

print()
print("=" * 80)
print("VERDICT per α:")
print("=" * 80)

for r in all_results:
    alpha = r['alpha']
    issues = []

    if r['rho_L_over_rho_total_recomb'] > 0.05:
        issues.append(f"CMB safety FAIL (ρ_Λ/ρ_tot = {r['rho_L_over_rho_total_recomb']*100:.2f}% > 5%)")
    if r['tension_coverage'] < 0.5:
        issues.append(f"Tension coverage too small ({r['tension_coverage']*100:.1f}% < 50%)")
    if r['tension_coverage'] > 1.5:
        issues.append(f"Tension coverage overshoot ({r['tension_coverage']*100:.1f}% > 150%)")
    if r['w_today'] < -1.001:
        issues.append(f"w_today < -1 (phantom, {r['w_today']:.4f})")
    if abs(r['H_today_TGP_over_LCDM'] - 1.0) > 0.001:
        issues.append(f"H_0(today) drift {(r['H_today_TGP_over_LCDM']-1)*100:.3f}%")

    status = "PASS" if not issues else "FAIL"
    print(f"\nα = {alpha}: {status}")
    if issues:
        for issue in issues:
            print(f"   ✗ {issue}")
    else:
        print(f"   ✓ All conditions met: ΔH₀/H₀ = {r['delta_H0_inferred_over_H0']*100:.2f}%, "
              f"w(0) = {r['w_today']:.4f}, ρ_Λ_recomb/ρ_tot = {r['rho_L_over_rho_total_recomb']*100:.2f}%")

print()
print("=" * 80)
print("STAGE 1 SUMMARY")
print("=" * 80)

# Find best alpha (max coverage subject to CMB safety + w >= -1)
best_alpha = None
best_coverage = -np.inf
for r in all_results:
    cmb_safe = r['rho_L_over_rho_total_recomb'] < 0.05
    w_ok = r['w_today'] >= -1.001
    in_target = 0.5 <= r['tension_coverage'] <= 1.5
    if cmb_safe and w_ok and in_target and r['tension_coverage'] > best_coverage:
        best_coverage = r['tension_coverage']
        best_alpha = r['alpha']

if best_alpha is not None:
    print(f"\n✅ BEST α = {best_alpha}: coverage = {best_coverage*100:.1f}% — passing all conditions")
else:
    print("\n❌ NO α satisfies all conditions simultaneously")
    print("\nDiagnostic — what blocks each α:")
    for r in all_results:
        issues = []
        if r['rho_L_over_rho_total_recomb'] > 0.05:
            issues.append("CMB-fail")
        if r['tension_coverage'] < 0.5:
            issues.append("under-coverage")
        if r['tension_coverage'] > 1.5:
            issues.append("over-coverage")
        if r['w_today'] < -1.001:
            issues.append("phantom")
        block = ",".join(issues) if issues else "OK"
        print(f"   α={r['alpha']:.3f}: {block} (coverage {r['tension_coverage']*100:.1f}%)")

print()
print("Note on self-consistency (B2 from M10.5/omicron2):")
print("This solver assumes γ(z) is FUNDAMENTAL coupling (no observer-frame rescaling).")
print("If sek04 self-consistency principle G(ψ) = G_0·(Φ_0/Φ) cancels boost,")
print("results above are upper bounds. Stage 1' addendum needed to verify.")
print()
print("=" * 80)
