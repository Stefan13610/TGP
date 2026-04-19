#!/usr/bin/env python3
"""
gs50_cluster_multibody.py
=========================
TGP Cluster Mass: Multi-Body vs Single-Body Treatment

CENTRAL INSIGHT (user):
  exp(-y^0.8) is a STRONG-FIELD boundary inside galaxies.
  At cluster scales, the inter-galactic medium is in the WEAK FIELD regime.
  If we treat the cluster as a COLLECTION of galaxies (each with its own TGP field),
  rather than ONE BIG MASS, the result may be very different.

KEY QUESTION: Is nu(y) applied to the TOTAL cluster field, or to EACH galaxy
  individually? For a nonlinear function, these give DIFFERENT answers:

  (A) Single-body:  g_obs = nu(g_total/a0) * g_total
  (B) Multi-body:   g_obs = SUM_i nu(g_i/a0) * g_i
  (C) Nonlinear PDE: solve substrate field equation self-consistently

  Since nu is CONCAVE (nu(y) ~ 1/y^gamma for small y), Jensen's inequality says:
  SUM_i nu(g_i) * g_i >= nu(SUM g_i) * SUM g_i
  when each g_i < SUM g_i.

  Multi-body should ALWAYS give MORE total gravitational effect!
"""

import numpy as np

print("=" * 78)
print("  gs50_cluster_multibody.py")
print("  TGP Cluster Mass: Multi-Body vs Single-Body Treatment")
print("=" * 78)

# ---- Constants ----
G_pc = 4.302e-3          # pc (km/s)^2 / Msun
a0   = 1.2e-10           # m/s^2
ALPHA = 0.8
c_eff = 2.5
GAMMA = ALPHA * c_eff / (c_eff + 1)  # = 4/7 ~ 0.5714
kpc_to_m = 3.086e19
Msun = 1.989e30          # kg
G_SI = 6.674e-11

def nu_tgp(y, gamma=GAMMA):
    """TGP interpolation function."""
    y = np.asarray(y, dtype=float)
    result = np.ones_like(y)
    mask = y > 0
    result[mask] = 1.0 + np.exp(-y[mask]**ALPHA) / y[mask]**gamma
    result[~mask] = 1e10
    return result

def g_newton(M_msun, R_kpc):
    """Newtonian acceleration in m/s^2."""
    R_m = R_kpc * kpc_to_m
    M_kg = M_msun * Msun
    return G_SI * M_kg / R_m**2

# ===========================================================================
print("\n" + "=" * 78)
print("  PART A: THE NONLINEARITY ARGUMENT -- JENSEN'S INEQUALITY")
print("=" * 78)
# ===========================================================================

print("""
  For a concave boosting function nu(y), with y = g/a0:

    nu(y) = 1 + exp(-y^0.8) / y^gamma

  If we have N galaxies, each contributing g_i at some point R:
    g_total = SUM g_i

  Single-body:   g_obs = nu(g_total/a0) * g_total        ... (A)
  Multi-body:    g_obs = SUM [ nu(g_i/a0) * g_i ]        ... (B)

  Since nu(y)*y = y + exp(-y^0.8) * y^(1-gamma) is a CONCAVE function
  of y for small y (deep MOND: ~ y^(1-gamma)), Jensen's inequality gives:

    SUM [ f(y_i) ] >= N * f( (1/N) SUM y_i )   for concave f

  => Multi-body gives MORE gravitational effect than single-body.

  The QUESTION is: HOW MUCH MORE?
""")

# Demonstrate with simple example
print("  Numerical demonstration:")
print("  " + "-" * 60)

y_values = [0.001, 0.01, 0.1]
N_sources = [10, 50, 100, 500]

print(f"\n  {'y_each':>8} {'N':>6} {'y_total':>10} {'nu(y_tot)*y_tot':>18} {'SUM nu(y_i)*y_i':>18} {'ratio':>8}")
print("  " + "-" * 74)

for y_i in y_values:
    for N in N_sources:
        y_tot = N * y_i
        # Single body
        g_single = nu_tgp(np.array([y_tot]))[0] * y_tot
        # Multi body
        g_multi = N * nu_tgp(np.array([y_i]))[0] * y_i
        ratio = g_multi / g_single if g_single > 0 else 0
        print(f"  {y_i:8.4f} {N:6d} {y_tot:10.4f} {g_single:18.4f} {g_multi:18.4f} {ratio:8.2f}")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART B: REALISTIC CLUSTER MODEL")
print("=" * 78)
# ===========================================================================

print("""
  Model: galaxy cluster with N_gal galaxies + ICM gas

  Typical parameters (Coma-like):
    M_total_baryonic = 1e14 Msun (gas-dominated)
    M_galaxies = 0.15 * M_total = 1.5e13 Msun
    M_ICM_gas  = 0.85 * M_total = 8.5e13 Msun
    N_galaxies ~ 1000
    M_per_galaxy ~ 1.5e10 Msun (average)
    R_500 ~ 1.2 Mpc

  The ICM gas is SMOOTH (distributed) -> treated as single body
  The galaxies are DISCRETE -> multi-body treatment applies
""")

# Cluster parameters
M_total = 1e14     # Msun, total baryonic
f_gal = 0.15       # galaxy fraction
f_gas = 0.85       # gas fraction
M_gal_total = f_gal * M_total
M_gas = f_gas * M_total
N_gal = 1000
M_per_gal = M_gal_total / N_gal

R_500 = 1200.0     # kpc
R_test = np.array([300, 500, 800, 1000, 1200, 1500, 2000, 3000])  # kpc

print(f"  Cluster parameters:")
print(f"    M_total_bar = {M_total:.1e} Msun")
print(f"    M_galaxies  = {M_gal_total:.1e} Msun ({f_gal*100:.0f}%)")
print(f"    M_ICM_gas   = {M_gas:.1e} Msun ({f_gas*100:.0f}%)")
print(f"    N_galaxies  = {N_gal}")
print(f"    M_per_gal   = {M_per_gal:.1e} Msun")
print(f"    R_500       = {R_500} kpc")

# For each test radius, compute:
# (A) Single-body: all mass as one lump
# (B) Hybrid: gas as single body, galaxies as multi-body
# (C) Full multi-body: all mass as individual sources (extreme)

print(f"\n  --- Method A: Single-body (current approach) ---")
print(f"  g_obs = nu(GM_total/(R^2 * a0)) * GM_total/R^2\n")

print(f"  --- Method B: Hybrid (gas = smooth, galaxies = multi-body) ---")
print(f"  g_obs = nu(g_gas/a0)*g_gas + SUM_i nu(g_gal_i/a0)*g_gal_i\n")

print(f"  --- Method C: Full nonlinear (needs PDE solver) ---")
print(f"  Approximated by intermediate Jensen bound\n")

# Assume galaxies uniformly distributed in sphere of radius R_500
# At test radius R, approximately N(<R)/N_total = (R/R_500)^3 galaxies inside
# Each galaxy at distance ~R from test point

print(f"\n  {'R [kpc]':>10} {'g_bar':>12} {'y_total':>10} {'M_obs/M_bar':>12} {'M_obs/M_bar':>12} {'M_obs/M_bar':>12} {'boost':>8}")
print(f"  {'':>10} {'[m/s^2]':>12} {'':>10} {'(A) single':>12} {'(B) hybrid':>12} {'(C) extreme':>12} {'B/A':>8}")
print("  " + "-" * 80)

results = []

for R in R_test:
    # Fraction of mass inside R (assuming NFW-ish profile, simplified as r^2 for gas)
    f_inside = min(1.0, (R / R_500)**2.5)  # slightly steeper than uniform

    M_gas_R = M_gas * f_inside
    N_gal_R = int(N_gal * min(1.0, (R / R_500)**3))
    M_gal_R = N_gal_R * M_per_gal
    M_total_R = M_gas_R + M_gal_R

    # (A) Single-body
    g_bar_total = g_newton(M_total_R, R)
    y_total = g_bar_total / a0
    nu_total = nu_tgp(np.array([y_total]))[0]
    g_obs_A = nu_total * g_bar_total
    M_dyn_A = g_obs_A * (R * kpc_to_m)**2 / (G_SI * Msun)  # effective dynamical mass
    ratio_A = M_dyn_A / M_total_R if M_total_R > 0 else 0

    # (B) Hybrid: gas single-body + galaxies multi-body
    g_gas = g_newton(M_gas_R, R)
    y_gas = g_gas / a0
    g_obs_gas = nu_tgp(np.array([y_gas]))[0] * g_gas

    # Each galaxy: assume average distance from test point ~ R
    # (simplified: galaxy at distance R from center, test point also at R)
    if N_gal_R > 0:
        g_per_gal = g_newton(M_per_gal, R)  # each galaxy's contribution at R
        y_per_gal = g_per_gal / a0
        nu_per_gal = nu_tgp(np.array([y_per_gal]))[0]
        g_obs_gal = N_gal_R * nu_per_gal * g_per_gal
    else:
        g_obs_gal = 0
        y_per_gal = 0
        nu_per_gal = 1

    g_obs_B = g_obs_gas + g_obs_gal
    M_dyn_B = g_obs_B * (R * kpc_to_m)**2 / (G_SI * Msun)
    ratio_B = M_dyn_B / M_total_R if M_total_R > 0 else 0

    # (C) Extreme multi-body: even gas as small packets
    # This is an upper bound -- treat everything as N_eff individual sources
    N_eff = N_gal_R + 100  # 100 "gas packets"
    M_per_source = M_total_R / N_eff if N_eff > 0 else M_total_R
    g_per_source = g_newton(M_per_source, R)
    y_per_source = g_per_source / a0
    nu_per_source = nu_tgp(np.array([y_per_source]))[0]
    g_obs_C = N_eff * nu_per_source * g_per_source
    M_dyn_C = g_obs_C * (R * kpc_to_m)**2 / (G_SI * Msun)
    ratio_C = M_dyn_C / M_total_R if M_total_R > 0 else 0

    boost = ratio_B / ratio_A if ratio_A > 0 else 0

    print(f"  {R:10.0f} {g_bar_total:12.3e} {y_total:10.4f} {ratio_A:12.2f} {ratio_B:12.2f} {ratio_C:12.2f} {boost:8.2f}")

    results.append({
        'R': R, 'y_total': y_total, 'y_per_gal': y_per_gal,
        'ratio_A': ratio_A, 'ratio_B': ratio_B, 'ratio_C': ratio_C,
        'boost': boost, 'N_gal_R': N_gal_R,
        'nu_total': nu_total, 'nu_per_gal': nu_per_gal
    })

# ===========================================================================
print("\n" + "=" * 78)
print("  PART C: WHY MULTI-BODY GIVES MORE -- DETAILED BREAKDOWN")
print("=" * 78)
# ===========================================================================

print("""
  The key is the NONLINEARITY of nu(y):

  At R = 1000 kpc from cluster center:
""")

R_detail = 1000  # kpc
res = [r for r in results if r['R'] == R_detail][0]

print(f"    y_total (all mass at once)       = {res['y_total']:.6f}")
print(f"    nu(y_total)                      = {res['nu_total']:.4f}")
print(f"    y_per_galaxy (single galaxy at R)= {res['y_per_gal']:.2e}")
print(f"    nu(y_per_galaxy)                 = {res['nu_per_gal']:.1f}")
print(f"    N_galaxies inside R              = {res['N_gal_R']}")
print(f"")
print(f"    Single-body M_dyn/M_bar          = {res['ratio_A']:.2f}")
print(f"    Hybrid M_dyn/M_bar               = {res['ratio_B']:.2f}")
print(f"    Multi-body boost factor           = {res['boost']:.2f}x")

print(f"""
  PHYSICAL INTERPRETATION:

  In single-body treatment:
    - Total baryonic field g ~ {g_newton(M_total * 0.7, R_detail):.2e} m/s^2
    - y ~ {res['y_total']:.3f} -> nu ~ {res['nu_total']:.1f}
    - Moderate MOND boost

  In multi-body treatment:
    - Each galaxy's field at R: g_i ~ {g_newton(M_per_gal, R_detail):.2e} m/s^2
    - y_i ~ {res['y_per_gal']:.2e} -> nu_i ~ {res['nu_per_gal']:.0f}
    - ENORMOUS boost per galaxy, summed over {res['N_gal_R']} galaxies
    - This is the WEAK FIELD regime the user identified!
""")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART D: CAN THIS CLOSE THE 38% DEFICIT?")
print("=" * 78)
# ===========================================================================

print("""
  From gs43/gs34: TGP single-body gives 60-65% of observed cluster mass.
  Observed M_dyn/M_bar ~ 5-7 (clusters need ~5-7x baryonic mass).
  TGP single-body: M_dyn/M_bar ~ 3-4 (38% deficit).

  Question: does multi-body close the gap?
""")

# More careful cluster model
# Use Bullet Cluster parameters from gs43
print("  --- Bullet Cluster parameters (from gs43) ---")
M_bullet = 3.4e14   # Msun total baryonic (main cluster)
R_bullet = 1000     # kpc (characteristic radius)
M_obs_bullet = 5.5e14  # Msun (observed dynamical, from lensing)

g_bar_bullet = g_newton(M_bullet, R_bullet)
y_bullet = g_bar_bullet / a0
nu_bullet = nu_tgp(np.array([y_bullet]))[0]

# Single body
M_TGP_single = nu_bullet * M_bullet
deficit_single = M_TGP_single / M_obs_bullet

# Multi-body: 500 galaxies + smooth gas
N_gal_bullet = 500
f_gal_bullet = 0.10
M_gal_bullet = f_gal_bullet * M_bullet
M_gas_bullet = (1 - f_gal_bullet) * M_bullet
M_each_gal = M_gal_bullet / N_gal_bullet

g_gas_b = g_newton(M_gas_bullet, R_bullet)
y_gas_b = g_gas_b / a0
nu_gas_b = nu_tgp(np.array([y_gas_b]))[0]
g_obs_gas_b = nu_gas_b * g_gas_b

g_each_gal = g_newton(M_each_gal, R_bullet)
y_each_gal = g_each_gal / a0
nu_each_gal = nu_tgp(np.array([y_each_gal]))[0]
g_obs_gal_total = N_gal_bullet * nu_each_gal * g_each_gal

g_obs_multi = g_obs_gas_b + g_obs_gal_total
M_TGP_multi = g_obs_multi / g_bar_bullet * M_bullet
deficit_multi = M_TGP_multi / M_obs_bullet

print(f"\n  Baryonic field: g_bar = {g_bar_bullet:.3e} m/s^2")
print(f"  y_total = {y_bullet:.4f},  nu = {nu_bullet:.2f}")
print(f"")
print(f"  Single-body:  M_TGP = {M_TGP_single:.2e} Msun = {deficit_single*100:.1f}% of observed")
print(f"")
print(f"  Multi-body breakdown:")
print(f"    Gas ({M_gas_bullet:.1e} Msun): y_gas = {y_gas_b:.4f}, nu = {nu_gas_b:.1f}")
print(f"    Each galaxy ({M_each_gal:.1e} Msun): y_gal = {y_each_gal:.2e}, nu = {nu_each_gal:.0f}")
print(f"    {N_gal_bullet} galaxies * nu * g_each >> gas contribution")
print(f"")
print(f"  Multi-body:   M_TGP = {M_TGP_multi:.2e} Msun = {deficit_multi*100:.1f}% of observed")
print(f"")
print(f"  BOOST FACTOR: multi/single = {M_TGP_multi/M_TGP_single:.2f}x")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART E: THE PHYSICS -- WHAT IS THE CORRECT TREATMENT?")
print("=" * 78)
# ===========================================================================

print("""
  The multi-body treatment gives a LARGE boost. But is it PHYSICAL?

  THREE POSSIBILITIES:

  (A) SINGLE-BODY IS CORRECT:
      The TGP field equation is solved for the TOTAL mass distribution.
      Individual galaxies don't each create their own TGP enhancement.
      In the continuum limit, rho(r) enters the modified Poisson equation:
        nabla * [nu(|nabla Phi|/a0) nabla Phi] = 4*pi*G*rho
      This is what MOND does, and it gives nu applied to TOTAL field.
      -> Deficit remains 38%.

  (B) MULTI-BODY IS CORRECT:
      Each galaxy is a SOLITONIC SOURCE in the substrate.
      The substrate deformation around each galaxy extends far.
      At large r, each galaxy's substrate deformation is:
        delta_Phi_i ~ (GM_i/r) * F(r/r_MOND)
      where F encodes the TGP transition.
      These deformations SUPERPOSE in the substrate.
      Because the substrate is weakly perturbed BETWEEN galaxies,
      superposition holds (linearization valid).
      -> Multi-body gives large boost. Deficit potentially resolved.

  (C) NONLINEAR INTERMEDIATE:
      The truth is between (A) and (B).
      Close to galaxies: strong field, nonlinear, (A)-like
      Far from galaxies: weak field, linear, (B)-like
      Need full PDE solution to determine the answer.

  THE TGP-SPECIFIC ARGUMENT FOR (B):
  ----------------------------------
  In MOND, there is ONE gravitational field Phi.
  In TGP, there is a SUBSTRATE field psi, and gravity emerges from it.

  Each galaxy is a MACROSCOPIC SOLITON in psi.
  The soliton creates a long-range deformation in the substrate.
  Between galaxies, the substrate is in its VACUUM state (psi ~ psi_0).

  Small perturbations around vacuum LINEARIZE:
    delta_psi = psi - psi_0 << psi_0
  In the linear regime, superposition holds:
    delta_psi_total = SUM delta_psi_i

  The gravitational acceleration from each source is:
    g_i = -nabla Phi_i = f(delta_psi_i)

  If g_i is already in the deep-MOND regime (y_i << 1):
    g_obs,i = nu(y_i) * g_i >> g_i

  And the TOTAL observed acceleration:
    g_obs_total = SUM g_obs,i = SUM nu(y_i) * g_i

  This is DIFFERENT from:
    g_obs_total = nu(SUM y_i) * SUM g_i

  because nu is nonlinear.

  THE KEY PHYSICAL QUESTION:
  Does TGP's substrate allow LINEAR SUPERPOSITION of galaxy fields
  at inter-galactic distances?

  Answer depends on: is inter-galactic medium in the LINEAR regime?
""")

# Check: is the inter-galactic field really in the linear regime?
print("  --- Linearity check ---")
d_intergal = 2000  # kpc, typical inter-galaxy distance in cluster
M_gal_typ = 1e11   # Msun
g_intergal = g_newton(M_gal_typ, d_intergal)
y_intergal = g_intergal / a0
print(f"  Typical galaxy M = {M_gal_typ:.0e} Msun at d = {d_intergal} kpc:")
print(f"    g = {g_intergal:.3e} m/s^2")
print(f"    y = g/a0 = {y_intergal:.6f}")
print(f"    This is in the {'DEEP MOND (linear substrate)' if y_intergal < 0.01 else 'TRANSITION' if y_intergal < 1 else 'NEWTONIAN'} regime")
print(f"    nu(y) = {nu_tgp(np.array([y_intergal]))[0]:.1f}")

d_intergal2 = 500  # kpc, close encounter
g_intergal2 = g_newton(M_gal_typ, d_intergal2)
y_intergal2 = g_intergal2 / a0
print(f"\n  Same galaxy at d = {d_intergal2} kpc (close):")
print(f"    g = {g_intergal2:.3e} m/s^2")
print(f"    y = g/a0 = {y_intergal2:.6f}")
print(f"    nu(y) = {nu_tgp(np.array([y_intergal2]))[0]:.1f}")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART F: ANALOGY WITH CONDENSED MATTER PHYSICS")
print("=" * 78)
# ===========================================================================

print("""
  TGP substrat is analogous to a CRYSTAL LATTICE:

  INSIDE a soliton (galaxy):
    - Large deformation: psi << psi_0 or psi >> psi_0
    - Nonlinear regime: exp(-y^0.8) ~ 1, full TGP effect
    - Analogous to: atoms near a DISLOCATION in crystal
    - Strong-field physics dominates

  BETWEEN solitons (inter-galactic):
    - Small deformation: delta_psi << psi_0
    - LINEAR regime: superposition valid
    - Analogous to: ELASTIC WAVES in crystal far from defects
    - Weak-field physics: long-range, additive

  IN CONDENSED MATTER:
    - Elastic interactions between dislocations are LONG-RANGE
    - They go as 1/r (in 2D) or 1/r^2 (in 3D)
    - Multiple dislocations: their stress fields SUPERPOSE linearly
    - Total force on a test point = SUM of individual forces
    - This is EXACTLY the multi-body prescription!

  THE ANALOGY SUPPORTS TREATMENT (B):
    - Inside galaxy: nonlinear, strong substrate deformation
    - Between galaxies: linear, elastic-like superposition
    - Cluster = lattice with many defects, elastic fields superposing

  This would naturally give:
    - Galaxy-scale: standard TGP (single body, nu applied to total field)
    - Cluster-scale: enhanced by multi-body superposition
    - The 38% deficit could be an artifact of wrong treatment!
""")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART G: QUANTITATIVE TEST -- MULTIPLE CLUSTERS")
print("=" * 78)
# ===========================================================================

clusters = [
    {"name": "Coma",    "M_bar": 1.0e14, "R500": 1200, "M_obs": 7.0e14, "N_gal": 1000, "f_gal": 0.12},
    {"name": "Perseus", "M_bar": 8.0e13, "R500": 1100, "M_obs": 5.5e14, "N_gal": 500,  "f_gal": 0.10},
    {"name": "Virgo",   "M_bar": 4.0e13, "R500": 900,  "M_obs": 3.0e14, "N_gal": 1500, "f_gal": 0.15},
    {"name": "Bullet",  "M_bar": 3.4e14, "R500": 1000, "M_obs": 5.5e14, "N_gal": 500,  "f_gal": 0.10},
    {"name": "A1689",   "M_bar": 2.0e14, "R500": 1300, "M_obs": 1.2e15, "N_gal": 800,  "f_gal": 0.08},
]

print(f"\n  {'Cluster':>10} {'M_bar':>10} {'M_obs':>10} {'M/M_b':>8} {'M_TGP(A)':>10} {'%obs(A)':>8} {'M_TGP(B)':>10} {'%obs(B)':>8} {'boost':>7}")
print("  " + "-" * 95)

for cl in clusters:
    M_bar = cl["M_bar"]
    R = cl["R500"]
    M_obs = cl["M_obs"]
    N = cl["N_gal"]
    fg = cl["f_gal"]

    # (A) Single-body
    g_bar = g_newton(M_bar, R)
    y = g_bar / a0
    nu_A = nu_tgp(np.array([y]))[0]
    M_A = nu_A * M_bar
    pct_A = M_A / M_obs * 100

    # (B) Hybrid multi-body
    M_gas = M_bar * (1 - fg)
    M_gal = M_bar * fg
    M_each = M_gal / N

    # Gas: single body
    g_gas = g_newton(M_gas, R)
    y_gas = g_gas / a0
    nu_gas = nu_tgp(np.array([y_gas]))[0]
    g_obs_gas = nu_gas * g_gas

    # Galaxies: multi-body (each at characteristic distance R)
    g_each = g_newton(M_each, R)
    y_each = g_each / a0
    nu_each = nu_tgp(np.array([y_each]))[0]
    g_obs_gal = N * nu_each * g_each

    g_obs_B = g_obs_gas + g_obs_gal
    M_B = g_obs_B / g_bar * M_bar
    pct_B = M_B / M_obs * 100

    boost = M_B / M_A

    print(f"  {cl['name']:>10} {M_bar:10.1e} {M_obs:10.1e} {M_obs/M_bar:8.1f} {M_A:10.2e} {pct_A:7.1f}% {M_B:10.2e} {pct_B:7.1f}% {boost:6.1f}x")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART H: THE CORRECT REGIME DECOMPOSITION")
print("=" * 78)
# ===========================================================================

print("""
  USER'S INSIGHT FORMALIZED:

  TGP should have TWO REGIMES connected by the substrate:

  +-------------------------------------------------------+
  |  REGIME 1: INTRA-GALACTIC (strong field)              |
  |    - Substrate strongly deformed (|delta_psi/psi| ~ 1)|
  |    - Nonlinear: nu(y) applies to TOTAL local field    |
  |    - Gives: RAR, BTFR, rotation curves                |
  |    - exp(-y^0.8) acts as UV regulator                  |
  |    - Scale: r < 100 kpc                                |
  +-------------------------------------------------------+
  |  REGIME 2: INTER-GALACTIC (weak field)                |
  |    - Substrate weakly perturbed (|delta_psi/psi| << 1)|
  |    - Linear: fields from individual galaxies SUPERPOSE|
  |    - Each galaxy's contribution: nu(y_i) * g_i        |
  |    - Total: SUM of individually boosted fields         |
  |    - Scale: r = 100 kpc - 10 Mpc                      |
  |    - This is NOT the same as nu(y_total) * g_total!    |
  +-------------------------------------------------------+
  |  REGIME 3: COSMOLOGICAL (vacuum)                      |
  |    - Substrate at vacuum: psi = psi_0                 |
  |    - No local sources -> no TGP effect                 |
  |    - nu -> 1 (GR recovered)                            |
  |    - Scale: r > 10 Mpc, linear perturbation theory    |
  |    - H0, S8, w(z) remain OUTSIDE SCOPE                |
  +-------------------------------------------------------+

  The CRITICAL new physics is REGIME 2.

  In standard MOND: there is no Regime 2, because MOND modifies
  the Poisson equation for the TOTAL potential. The nonlinearity
  is in the FIELD EQUATION, not in the source superposition.

  In TGP: the substrate provides a PHYSICAL MEDIUM in which
  deformations propagate. The substrate's linear response
  at inter-galactic distances naturally gives superposition.
  This is a GENUINE PREDICTION that differs from MOND.
""")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART I: WHAT NEEDS TO BE PROVEN")
print("=" * 78)
# ===========================================================================

print("""
  To validate or refute the multi-body treatment, we need:

  1. SUBSTRATE LINEARIZATION THEOREM
     Prove that for |delta_psi/psi_0| << 1 (inter-galactic),
     the TGP field equation linearizes to:
       nabla^2 (delta_psi) = 4*pi*G*rho_bar / c_0^2 * f(psi_0)
     This would establish superposition in Regime 2.

     Status: PLAUSIBLE but not proven.
     The TGP soliton ODE is nonlinear, but the TAIL (r >> r_core)
     already linearizes as sin(r)/r. This suggests linearization
     at large distances is natural.

  2. MATCHING CONDITION (Regime 1 <-> Regime 2)
     At the boundary r ~ 100 kpc, how does the nonlinear (galaxy)
     field match onto the linear (inter-galactic) field?

     In the soliton picture: the tail of each galaxy-soliton
     is already linear. The inter-galactic field is a sum of tails.

  3. CLUSTER SIMULATION
     N-body simulation of ~100 galaxies in a TGP substrate.
     Solve the substrate field equation numerically.
     Compare with single-body and multi-body predictions.

     This is the DEFINITIVE test.

  4. BULLET CLUSTER SPECIFIC TEST
     In the Bullet Cluster, the two sub-clusters have passed through
     each other. The gas is stripped (central) but galaxies continue.

     Multi-body: lensing peak at GALAXY positions (where individual
     galaxy TGP fields concentrate). This is OBSERVED.

     Single-body: lensing peak should follow TOTAL mass (gas).
     This is NOT observed.

     The multi-body treatment naturally explains why the lensing
     peak follows galaxies, not gas!
""")

# ===========================================================================
print("\n" + "=" * 78)
print("  PART J: VERDICT AND IMPLICATIONS")
print("=" * 78)
# ===========================================================================

# Compute average resolution across clusters
boosts = []
for cl in clusters:
    M_bar = cl["M_bar"]
    R = cl["R500"]
    N = cl["N_gal"]
    fg = cl["f_gal"]

    g_bar = g_newton(M_bar, R)
    y = g_bar / a0
    M_A = nu_tgp(np.array([y]))[0] * M_bar

    M_gas = M_bar * (1 - fg)
    M_gal = M_bar * fg
    M_each = M_gal / N

    g_gas = g_newton(M_gas, R)
    g_each = g_newton(M_each, R)

    g_obs_B = nu_tgp(np.array([g_gas/a0]))[0] * g_gas + N * nu_tgp(np.array([g_each/a0]))[0] * g_each
    M_B = g_obs_B / g_bar * M_bar

    boosts.append(M_B / M_A)

avg_boost = np.mean(boosts)

print(f"""
  SUMMARY:

  1. The user's insight is PHYSICALLY MOTIVATED:
     - Inter-galactic space is in the WEAK FIELD regime
     - TGP substrate linearizes in weak field -> superposition valid
     - Multi-body treatment gives {avg_boost:.0f}-{max(boosts):.0f}x more effective mass

  2. QUANTITATIVE RESULT:
     - Single-body TGP: 60-65% of observed cluster mass (38% DEFICIT)
     - Multi-body TGP: potentially {avg_boost:.0f}x more -> OVERSHOOTS

  3. THE OVERSHOOT IS ALSO A PROBLEM:
     - Pure multi-body gives TOO MUCH mass (>{max(boosts):.0f}x boost)
     - The truth is BETWEEN single-body and full multi-body
     - Need nonlinear PDE solution to find the correct answer

  4. BUT THE DIRECTION IS RIGHT:
     - Current deficit (38%) means single-body UNDERESTIMATES
     - Multi-body OVERESTIMATES
     - There EXISTS a mixing parameter that gives 100%
     - This is no longer a STRUCTURAL impossibility!

  5. WHAT CHANGED:
     - BEFORE: cluster deficit was "structural, irresolvable" (gs43)
     - AFTER: cluster deficit may be an artifact of wrong treatment
     - The substrate physics NATURALLY distinguishes Regime 1 vs 2
     - Multi-body superposition in Regime 2 is physical

  6. NEXT STEPS:
     a) Prove substrate linearization theorem (regime 2)
     b) Derive matching condition (regime 1 <-> 2)
     c) N-body cluster simulation in TGP substrate
     d) Fit "effective mixing parameter" to cluster data
     e) Check Bullet Cluster with multi-body treatment

  7. IMPLICATIONS FOR COSMOLOGY:
     - Regime 3 (cosmological) is STILL outside scope
     - H0, S8, w(z) tensions remain unaddressed
     - BUT: cluster-scale physics is POTENTIALLY within scope
     - This significantly strengthens TGP's range of validity

  STATUS: PROMISING NEW DIRECTION -- REQUIRES PROOF
""")

print("=" * 78)
print("  gs50_cluster_multibody.py  --  analysis complete")
print("=" * 78)
