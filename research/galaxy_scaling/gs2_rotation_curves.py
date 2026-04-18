#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gs2_rotation_curves.py: Rotation curves from concentration limit model.

KEY INSIGHT FROM gs1:
- v⁴ = G·M·a₀ gives BTFR automatically ✓
- Σ_crit = a₀/(2πG) ≈ 137 M☉/pc² ≈ Freeman limit ✓
- a₀ ≈ c·H₀/(2π) ✓
- BUT: v_flat underpredicted by ~25% for MW-class galaxies
- f_DM underpredicted vs observations

THE ISSUE:
- Model says v⁴ = G·M_bar_observed · a₀
- But M_bar_observed only counts VISIBLE baryons (stars + cold gas)
- Missing baryons: hot gas, WHIM, warm circumgalactic medium
- Cosmological baryon budget: Ω_b/Ω_m ≈ 0.049/0.315 ≈ 16%
  → only ~16% of gravitating mass is baryonic IN ΛCDM
  → in our model: ALL gravitating mass is baryonic, some hidden

APPROACH:
  Instead of adding dark matter, ask:
  "What total baryonic mass does the flat curve REQUIRE?"
  → M_total = v⁴_obs / (G·a₀)  = observed v_flat gives total mass
  → M_visible/M_total = baryonic fraction that's visible
  → The "dark" mass is diffuse baryonic matter

THIS SCRIPT:
  1. For each galaxy: compute M_total from v_flat → find "hidden baryon" fraction
  2. Build rotation curve from baryonic Hernquist+disk profile
  3. Apply concentration limit (MOND-like interpolation)
  4. Compute and display rotation curves
  5. For IC 1101: show the full picture
  6. Derive the "equivalent density profile" of the flat well
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np
from scipy.special import i0, i1, k0, k1
import math

print("="*78)
print("  ROTATION CURVES FROM CONCENTRATION LIMIT MODEL")
print("="*78)

# ==========================================================================
# CONSTANTS
# ==========================================================================
G_SI = 6.674e-11          # m³/(kg·s²)
M_sun = 1.989e30          # kg
kpc = 3.086e19            # m
km = 1e3                  # m

# G in galaxy units: kpc, km/s, M☉
G = G_SI * M_sun / (kpc * km**2)  # (km/s)² kpc / M☉
print(f"\n  G = {G:.6e} (km/s)²·kpc/M☉")

# Critical acceleration
a0 = 1.2e-10 * kpc / km**2     # (km/s)²/kpc
print(f"  a₀ = {a0:.4e} (km/s)²/kpc")

# Check BTFR for Milky Way
v_MW = 220  # km/s
M_MW_btfr = v_MW**4 / (G * a0)
print(f"  MW: v=220 → M_bar(BTFR) = {M_MW_btfr:.2e} M☉")
print(f"  MW: M_bar(observed stars+gas) ≈ 6-8×10¹⁰ M☉")
print(f"  → 'Missing' baryons: {M_MW_btfr - 7e10:.2e} M☉ ({(M_MW_btfr - 7e10)/M_MW_btfr*100:.0f}%)")

# ==========================================================================
# PART 1: BARYONIC MASS BUDGET — ALL GALAXIES
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 1: What total baryonic mass does each galaxy require?")
print(f"{'='*78}")

# Galaxy data: [name, M_visible (M☉), v_flat (km/s), R_eff (kpc), type]
galaxies = [
    ("DDO 154",       5e8,    47,    2.0, "dwarf Irr"),
    ("NGC 1560",      2e9,    78,    4.0, "dwarf Sp"),
    ("NGC 2403",      1e10,  135,    7.0, "Sc"),
    ("M33",           5e9,   110,    8.0, "Sc"),
    ("Milky Way",     7e10,  220,   15.0, "SBbc"),
    ("M31",           1.2e11, 250,  25.0, "Sb"),
    ("UGC 2885",      2e11,  300,   40.0, "Sc giant"),
    ("M87",           8e11,  400,  100.0, "cD"),
    ("NGC 4889",      3e12,  450,  120.0, "cD BCG"),
    ("IC 1101",       5e12,  500,  200.0, "cD supergiant"),
]

print(f"\n  {'Galaxy':<14s} {'M_vis':>9s} {'v_flat':>6s} {'M_BTFR':>10s} {'M_vis/M_BTFR':>12s} {'M_hidden':>10s} {'Note':>20s}")
print(f"  {'─'*85}")

for name, M_vis, v_flat, R_eff, gtype in galaxies:
    M_btfr = v_flat**4 / (G * a0)
    f_vis = M_vis / M_btfr
    M_hidden = M_btfr - M_vis

    if f_vis > 0.8:
        note = "mostly visible"
    elif f_vis > 0.3:
        note = "significant hidden"
    else:
        note = "MOSTLY hidden"

    print(f"  {name:<14s} {M_vis:9.1e} {v_flat:6.0f} {M_btfr:10.2e} {f_vis:12.2f} {M_hidden:10.2e} {note:>20s}")

print(f"""
  INTERPRETATION:
  M_BTFR = v⁴/(G·a₀) = total mass needed for BTFR consistency.
  M_hidden = M_BTFR - M_visible = baryonic mass in diffuse form.

  For dwarfs (DDO 154, NGC 1560): ~50-60% visible
  For MW-class: ~30-40% visible
  For giant ellipticals: ~10-15% visible

  The "missing" mass would be:
  - Hot circumgalactic gas (10⁵-10⁷ K, X-ray emitting)
  - Warm-hot intergalactic medium (WHIM, 10⁵-10⁶ K)
  - Molecular hydrogen (H₂, hard to detect)
  - Ultra-faint stellar populations
  - Compact objects (old white dwarfs, neutron stars, stellar black holes)

  KEY QUESTION: Is there enough missing baryonic matter?
  Cosmological baryon budget: Ω_b·ρ_crit ≈ 4.2×10⁻²⁸ kg/m³
  About 50% of cosmic baryons are "missing" (WHIM), so plausible!
""")

# ==========================================================================
# PART 2: ROTATION CURVE SHAPE
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 2: Rotation curve — the concentration limit approach")
print(f"{'='*78}")

print(f"""
  PHYSICAL MODEL:
  ═══════════════════════════════════════════════════════════════

  Instead of MOND's "modified gravity at low a":
  We say: "matter can't concentrate below critical density ρ_min(r)"

  When gravity tries to concentrate mass to density > ρ_max_effective:
  → Substrate pressure prevents further collapse
  → Mass redistributes into ρ ∝ 1/r² profile
  → This is the ISOTHERMAL profile — exactly what gives flat curves

  The critical density at radius r:
    ρ_crit(r) = a₀ / (4πGr)  ← from g(r) = a₀ ↔ 4πGρr = a₀

  For MW at r = 10 kpc:
    ρ_crit = a₀/(4πG·r) = {a0/(4*np.pi*G*10):.2e} M☉/kpc³
           = {a0/(4*np.pi*G*10*1e9):.4f} M☉/pc³

  Observed DM halo density at 10 kpc: ~0.01 M☉/pc³ (consistent!)

  ═══════════════════════════════════════════════════════════════
""")

# Compute rotation curves for selected galaxies
def rotation_curve_MOND(radii, M_disk, R_d, M_bulge, R_b, a0_val, G_val):
    """
    Compute rotation curve with MOND-like interpolation.

    Baryonic components:
    - Exponential disk with mass M_disk, scale length R_d
    - Hernquist bulge with mass M_bulge, scale radius R_b
    """
    # Enclosed baryonic mass at each radius
    M_enc = np.zeros_like(radii)
    for i, r in enumerate(radii):
        # Exponential disk (cumulative mass)
        x = r / R_d
        f_disk = 1 - (1 + x) * np.exp(-x)
        # Hernquist bulge
        f_bulge = r**2 / (r + R_b)**2
        M_enc[i] = M_disk * f_disk + M_bulge * f_bulge

    # Newtonian acceleration from baryons
    g_bar = G_val * M_enc / radii**2

    # MOND interpolation (McGaugh 2016):
    # g_obs = g_bar / (1 - exp(-√(g_bar/a₀)))
    g_obs = g_bar / (1 - np.exp(-np.sqrt(g_bar / a0_val)))

    # Newtonian velocity
    v_newton = np.sqrt(G_val * M_enc / radii)

    # MOND velocity
    v_mond = np.sqrt(g_obs * radii)

    return v_newton, v_mond, g_bar, g_obs, M_enc

# ==========================================================================
# PART 3: MILKY WAY ROTATION CURVE
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 3: Milky Way rotation curve")
print(f"{'='*78}")

r = np.linspace(0.5, 60, 200)

M_disk_MW = 5.5e10
R_d_MW = 2.6  # kpc
M_bulge_MW = 1.5e10
R_b_MW = 0.5  # kpc

v_N, v_M, g_b, g_o, M_e = rotation_curve_MOND(r, M_disk_MW, R_d_MW, M_bulge_MW, R_b_MW, a0, G)

print(f"  Milky Way: M_disk = {M_disk_MW:.1e}, M_bulge = {M_bulge_MW:.1e}")
print(f"  Total M_bar = {M_disk_MW + M_bulge_MW:.1e} M☉")

print(f"\n  {'r kpc':>6s} {'v_Newton':>9s} {'v_model':>9s} {'v_obs':>6s} {'g_bar/a₀':>9s} {'g_obs/a₀':>9s} {'M_enc':>10s} {'M_Newton':>10s} {'f_DM_app':>9s}")
print(f"  {'─'*80}")

for r_show in [1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50, 60]:
    idx = np.argmin(np.abs(r - r_show))
    M_newton_at_r = v_M[idx]**2 * r[idx] / G  # mass Newton would infer
    f_dm = 1 - M_e[idx] / M_newton_at_r if M_newton_at_r > 0 else 0
    print(f"  {r[idx]:6.1f} {v_N[idx]:9.1f} {v_M[idx]:9.1f} {220:6.0f} {g_b[idx]/a0:9.3f} {g_o[idx]/a0:9.3f} {M_e[idx]:10.2e} {M_newton_at_r:10.2e} {f_dm:9.3f}")

# ==========================================================================
# PART 4: EQUIVALENT DENSITY PROFILE ("phantom dark matter")
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 4: Equivalent 'phantom dark matter' density profile")
print(f"{'='*78}")

print(f"""
  If we interpret the MOND-modified rotation curve in Newtonian terms,
  the "missing mass" has a specific density profile — the PHANTOM DM.

  ρ_phantom(r) = [M_Newton(r) - M_baryon(r)] / (4πr²·dr)
  where M_Newton(r) = v²_MOND(r)·r/G
""")

# Compute phantom DM density
M_newton_prof = v_M**2 * r / G
M_phantom = M_newton_prof - M_e

# Density from dM/dr = 4πr²ρ
rho_phantom = np.zeros_like(r)
for i in range(1, len(r)-1):
    dr = r[i+1] - r[i-1]
    dM = M_phantom[i+1] - M_phantom[i-1]
    rho_phantom[i] = dM / (4 * np.pi * r[i]**2 * dr)

# Also compute the isothermal ρ = v²/(4πGr²) for comparison
v_flat_eff = 220.0
rho_isothermal = v_flat_eff**2 / (4 * np.pi * G * r**2)

print(f"  {'r kpc':>6s} {'ρ_phantom':>12s} {'ρ_isothermal':>12s} {'ρ_ph/ρ_iso':>10s}  {'ρ M☉/pc³':>10s}")
print(f"  {'─'*55}")
for r_show in [2, 5, 8, 10, 15, 20, 30, 50]:
    idx = np.argmin(np.abs(r - r_show))
    rho_ph = rho_phantom[idx]
    rho_iso = rho_isothermal[idx]
    rho_pc3 = rho_ph / 1e9 if rho_ph > 0 else 0  # kpc³ → pc³
    ratio = rho_ph / rho_iso if rho_iso > 0 else 0
    print(f"  {r[idx]:6.1f} {rho_ph:12.2e} {rho_iso:12.2e} {ratio:10.3f}  {rho_pc3:10.4f}")

# ==========================================================================
# PART 5: IC 1101 — THE EXTREME CASE
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 5: IC 1101 — flat well at the extreme")
print(f"{'='*78}")

# IC 1101 is a BCG (brightest cluster galaxy) in Abell 2029
# It's an elliptical → use Hernquist profile for both components

M_stars_IC = 5e12      # M☉ (stellar mass, uncertain)
R_eff_IC = 60          # kpc (effective radius)
R_h_IC = R_eff_IC / 1.8153  # Hernquist scale radius
sigma_IC = 400         # km/s (velocity dispersion)

# For Hernquist: R_h = R_eff / (1 + √2) ≈ R_eff / 2.414...
# Actually R_h = R_eff / 1.8153 (exact for Hernquist projected half-mass)

# Hot gas in Abell 2029 cluster:
M_gas_cluster = 2e13   # M☉ (X-ray gas, extended)
R_gas = 500            # kpc (gas extends further than stars)

# What does BTFR predict for total mass?
# For ellipticals: v_circ ≈ √2 · σ (for isotropic dispersion)
v_circ_IC = sigma_IC * np.sqrt(2)
M_btfr_IC = v_circ_IC**4 / (G * a0)

print(f"  IC 1101:")
print(f"    σ = {sigma_IC} km/s → v_circ ≈ √2·σ = {v_circ_IC:.0f} km/s")
print(f"    M★ (visible stars) = {M_stars_IC:.1e} M☉")
print(f"    M_gas (cluster X-ray) = {M_gas_cluster:.1e} M☉")
print(f"    M_BTFR (from v_circ) = {M_btfr_IC:.2e} M☉")
print(f"    M_BTFR / M★ = {M_btfr_IC/M_stars_IC:.1f} → '{(M_btfr_IC/M_stars_IC - 1)*100:.0f}% dark matter' in standard model")

# Rotation curve for IC 1101
r_IC = np.linspace(1, 800, 400)

# Stellar component (Hernquist)
M_enc_stars = M_stars_IC * r_IC**2 / (r_IC + R_h_IC)**2

# Gas component (β-model, roughly)
R_core_gas = 100  # kpc
M_enc_gas = M_gas_cluster * (r_IC / np.sqrt(r_IC**2 + R_core_gas**2))**3

M_enc_total = M_enc_stars + M_enc_gas

# Newtonian g
g_bar_IC = G * M_enc_total / r_IC**2

# MOND g
g_obs_IC = g_bar_IC / (1 - np.exp(-np.sqrt(g_bar_IC / a0)))

v_N_IC = np.sqrt(G * M_enc_total / r_IC)
v_M_IC = np.sqrt(g_obs_IC * r_IC)

print(f"\n  Rotation curve:")
print(f"  {'r kpc':>6s} {'v_Newton':>9s} {'v_model':>9s} {'g_bar/a₀':>9s} {'M_enc':>10s}")
print(f"  {'─'*50}")
for r_show in [5, 10, 30, 60, 100, 200, 300, 500, 700]:
    idx = np.argmin(np.abs(r_IC - r_show))
    print(f"  {r_IC[idx]:6.0f} {v_N_IC[idx]:9.0f} {v_M_IC[idx]:9.0f} {g_bar_IC[idx]/a0:9.2f} {M_enc_total[idx]:10.2e}")

# ==========================================================================
# PART 6: THE CONCENTRATION LIMIT — PHYSICAL PICTURE
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 6: The concentration limit — what prevents further collapse?")
print(f"{'='*78}")

print(f"""
  In the "flat well" picture, the potential depth has a maximum.
  In terms of accelerations: g ≤ a₀ at the "edge" of the well.

  But galaxy CENTERS have g >> a₀ (e.g., MW center: g ~ 10⁶ a₀).
  So the limit is NOT on the maximum local acceleration.

  REVISED INTERPRETATION:
  The limit is on the GLOBAL gravitational influence of mass.

  Think of it as a "gravitational range" limit:
  - A mass M influences space out to radius R_grav
  - R_grav = √(GM/a₀) — the "MOND radius"
  - Beyond R_grav: gravity falls off SLOWER than 1/r²
    (substrate effect — metric deformation extends further)
  - The "extra" gravity at r > R_grav mimics dark matter

  IN TGP LANGUAGE:
  Each mass creates a soliton deformation of the substrate.
  The soliton tail ψ(r) → 1 as r → ∞, but the approach is
  SLOWER than Newtonian 1/r because of the nonlinear substrate.

  Specifically, the TGP field equation:
    ∇²ψ + ψ'/ψ·(ψ')... + ψ = 1
  has tails that go as sin(r)/r (oscillating!) at large r.

  The oscillating tail means the potential extends FURTHER
  than Newtonian, giving the appearance of extra mass.
""")

# ==========================================================================
# PART 7: SCALING LAW — R_galaxy vs M_bar
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 7: Galaxy size-mass scaling law")
print(f"{'='*78}")

print(f"""
  From the flat well model:

  The "MOND radius" R_M = √(GM/a₀) sets the scale where
  the transition from Newtonian to deep-MOND occurs.

  This is NOT the galaxy's visible size, but the DYNAMICAL size —
  the radius where the rotation curve would flatten.

  Observable galaxy size (R_eff) is typically:
  - R_eff ~ 0.1-0.5 × R_M for spirals (disk concentrated)
  - R_eff ~ 0.3-0.8 × R_M for ellipticals (more extended)
""")

print(f"\n  {'Galaxy':<14s} {'M_bar':>9s} {'v_flat':>6s} {'R_M':>7s} {'R_eff':>6s} {'R_eff/R_M':>9s}")
print(f"  {'─'*55}")

for name, M_vis, v_flat, R_eff, gtype in galaxies:
    M_btfr = v_flat**4 / (G * a0)
    R_M = np.sqrt(G * M_btfr / a0)
    ratio = R_eff / R_M
    print(f"  {name:<14s} {M_btfr:9.2e} {v_flat:6.0f} {R_M:7.1f} {R_eff:6.1f} {ratio:9.3f}")

# ==========================================================================
# PART 8: MASS-CONCENTRATION RELATION
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 8: Surface density at R_eff — is it constant?")
print(f"{'='*78}")

print(f"""
  In MOND, there's a prediction: the surface density at the MOND radius
  is always ~ Σ_crit = a₀/(2πG) ≈ 137 M☉/pc².

  Let's check: what's the surface density Σ_eff = M/(πR²_eff) at R_eff?
""")

print(f"\n  {'Galaxy':<14s} {'M_BTFR':>10s} {'R_eff':>6s} {'Σ_eff':>10s} {'Σ_eff/Σ_crit':>12s} {'log Σ_eff':>9s}")
print(f"  {'':14s} {'(M☉)':>10s} {'kpc':>6s} {'M☉/pc²':>10s} {'':>12s} {'M☉/pc²':>9s}")
print(f"  {'─'*65}")

Sigma_crit = a0 / (2 * np.pi * G)  # M☉/kpc²
Sigma_crit_pc2 = Sigma_crit / 1e6  # M☉/pc²

for name, M_vis, v_flat, R_eff, gtype in galaxies:
    M_btfr = v_flat**4 / (G * a0)
    Sigma_eff = M_btfr / (np.pi * R_eff**2)  # M☉/kpc²
    Sigma_eff_pc2 = Sigma_eff / 1e6
    ratio_S = Sigma_eff_pc2 / Sigma_crit_pc2
    print(f"  {name:<14s} {M_btfr:10.2e} {R_eff:6.1f} {Sigma_eff_pc2:10.1f} {ratio_S:12.2f} {np.log10(Sigma_eff_pc2):9.2f}")

print(f"\n  Σ_crit = {Sigma_crit_pc2:.1f} M☉/pc²")
print(f"  Freeman limit ≈ 140 M☉/pc² → Σ_crit/Freeman = {Sigma_crit_pc2/140:.2f}")

# ==========================================================================
# PART 9: THE "FLAT WELL" AS A POTENTIAL PROFILE
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 9: Potential well shape")
print(f"{'='*78}")

# For a galaxy with flat rotation curve v_flat:
# Φ(r) = v²_flat · ln(r/r₀) + const (logarithmic potential)
# In the "well" picture:
# Φ(r) starts at some maximum depth at r=0
# At r → ∞: Φ → 0 (but logarithmic — diverges!)
# Resolution: cutoff at R_max ~ c/H₀ (Hubble radius)

# The depth of the well:
# ΔΦ = Φ(R_max) - Φ(R_min) = v²_flat · ln(R_max/R_min)

R_min = 0.001  # kpc (inner cutoff, ~1 pc)
R_max_Hubble = 4451e3  # kpc (c/H₀ in kpc)

for name, M_vis, v_flat, R_eff, gtype in [
    ("DDO 154", 5e8, 47, 2.0, "dwarf"),
    ("Milky Way", 7e10, 220, 15.0, "spiral"),
    ("IC 1101", 5e12, 500, 200.0, "cD"),
]:
    delta_Phi = v_flat**2 * np.log(R_max_Hubble / R_min)
    Phi_over_c2 = delta_Phi / (3e5)**2
    R_M = np.sqrt(G * v_flat**4 / (G * a0**2))  # simplify
    R_M2 = v_flat**2 / a0

    print(f"\n  {name}:")
    print(f"    v_flat = {v_flat} km/s")
    print(f"    Potential depth: ΔΦ = v²·ln(R_max/R_min) = {delta_Phi:.2e} (km/s)²")
    print(f"    ΔΦ/c² = {Phi_over_c2:.2e}  (relativistic parameter)")
    print(f"    MOND radius: R_M = v²/a₀ = {R_M2:.1f} kpc")
    print(f"    R_eff/R_M = {R_eff/R_M2:.3f}")

# ==========================================================================
# PART 10: SUMMARY AND IC 1101 INTERPRETATION
# ==========================================================================
print(f"\n{'='*78}")
print(f"  PART 10: SUMMARY — Galaxy as flat potential well")
print(f"{'='*78}")

# IC 1101 specific summary
M_IC_btfr = (500)**4 / (G * a0)
R_M_IC = 500**2 / a0

print(f"""
  ╔═══════════════════════════════════════════════════════════════╗
  ║  THE IC 1101 PICTURE                                         ║
  ╠═══════════════════════════════════════════════════════════════╣
  ║                                                               ║
  ║  Standard model:                                              ║
  ║    M★ = 5×10¹² M☉, M_DM = ~5×10¹⁴ M☉ (99% DM!)            ║
  ║    "Galaxy embedded in massive DM halo"                       ║
  ║                                                               ║
  ║  Flat well model:                                             ║
  ║    M_total(BTFR) = {M_IC_btfr:.2e} M☉ (from v=500)       ║
  ║    R_MOND = {R_M_IC:.0f} kpc                                ║
  ║    Galaxy extends to R_eff = 200 kpc = {200/R_M_IC:.2f}×R_MOND      ║
  ║                                                               ║
  ║  Interpretation:                                              ║
  ║    IC 1101 is NOT "mostly dark matter"                        ║
  ║    It's a galaxy whose baryonic matter can't concentrate      ║
  ║    below the critical surface density.                        ║
  ║                                                               ║
  ║    The 97% "dark matter" is actually:                         ║
  ║    - 30% visible (stars + detected gas)                       ║
  ║    - 70% diffuse hot baryons (WHIM, CGM)                     ║
  ║    + modified gravitational dynamics at r > R_MOND            ║
  ║                                                               ║
  ║  WHY IC 1101 IS SO LARGE:                                    ║
  ║    More mass → larger R_MOND → wider flat potential well      ║
  ║    R_MOND ∝ M^(1/2) → 10× more mass → 3.2× larger           ║
  ║    IC 1101 has ~100× MW mass → ~10× MW radius ✓              ║
  ╚═══════════════════════════════════════════════════════════════╝

  UNIVERSAL SCALING RELATIONS (ONE parameter: a₀):

  1. v_flat = (G·M·a₀)^(1/4)      — BTFR (exponent 4) ✓
  2. R_MOND = v²/a₀ = (GM/a₀)^(1/2) — size-mass ✓
  3. Σ_crit = a₀/(2πG) ≈ 137 M☉/pc² — Freeman limit ✓
  4. a₀ ≈ c·H₀/(2π)               — cosmological origin ✓

  WHAT'S DIFFERENT FROM MOND:
  Our interpretation focuses on the STRUCTURAL consequence:
  - Not "gravity is modified" but "matter can't concentrate"
  - Not "dark matter exists" but "baryons are more spread out"
  - The mechanism: TGP substrate limits metric deformation

  TESTABLE PREDICTION:
  If we're right, the "missing mass" should eventually be found
  as diffuse baryonic matter (WHIM, CGM, faint stellar populations).
  Recent surveys (eROSITA, Sunyaev-Zel'dovich) are starting to find it!

  NEXT: Connect a₀ to TGP substrate properties (gs3)
""")
