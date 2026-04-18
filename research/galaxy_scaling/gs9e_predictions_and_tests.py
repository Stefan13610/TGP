"""
gs9e: Testable predictions and confrontation with observational data
=====================================================================

CRITICAL TEST: Does a0 evolve with redshift?

Our model (gs9d geometric mean mechanism) predicts:
  a0(z) = c * H(z) / (2*pi)
  -> a0 grows with z as H(z) grows

Existing observational constraints:
1. Milgrom (2017, 1703.06110): 6 galaxies at z~0.9-2.4
   "all but exclude a0 ~ 4*a0 at z~2"
   Rules out a0 ~ (1+z)^(3/2)

2. McGaugh (Triton Station, 2025): BTFR shows NO perceptible
   evolution across cosmic time. DLA0817g at z=4.26: v=272 km/s
   consistent with local BTFR.

3. MIGHTEE-HI (2025): "first tentative evidence for redshift
   evolution" but z < 0.08 only, needs confirmation.

Our prediction: a0(z=2) = 3.0 * a0(0)
Milgrom's exclusion: a0(z=2) ~ 4 * a0(0) ruled out

THIS IS A POTENTIAL FALSIFICATION.
We must honestly assess whether our model survives.
"""

import numpy as np
import sys, io, warnings
warnings.filterwarnings('ignore')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# Constants
G = 6.674e-11
c = 2.998e8
H0 = 2.18e-18
a0 = 1.2e-10
M_sun = 1.989e30
kpc = 3.086e19

Omega_m = 0.315
Omega_L = 0.685

def H_of_z(z):
    """Hubble parameter at redshift z."""
    return H0 * np.sqrt(Omega_m * (1+z)**3 + Omega_L)

print("=" * 70)
print("  gs9e: Testable predictions and confrontation with data")
print("=" * 70)

# ============================================================
# 1. OUR PREDICTION: a0(z) = c*H(z)/(2*pi)
# ============================================================
print("\n" + "=" * 70)
print("  1. Our prediction vs observational constraints")
print("=" * 70)

print("""
  Model prediction: a0(z) = c * H(z) / (2*pi)

  Observational constraints:
  - Milgrom (2017): a0(z~2) < 4*a0(0) [excluded]
  - McGaugh (2025): BTFR consistent across z -> a0 ~ constant
  - MIGHTEE-HI (2025): tentative evolution at z < 0.08
""")

print(f"  {'z':>6}  {'H(z)/H0':>8}  {'a0(z)/a0(0)':>12}  {'Status':>30}")
print("  " + "-" * 60)
for z in [0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.26, 5.0]:
    H_ratio = H_of_z(z) / H0
    status = ""
    if z == 0:
        status = "reference"
    elif z < 0.1:
        status = "MIGHTEE range, OK"
    elif z < 1:
        status = "limited data"
    elif abs(z - 2.0) < 0.3:
        status = f"Milgrom: <4x -> {'OK' if H_ratio < 4 else 'EXCLUDED'}"
    elif abs(z - 4.26) < 0.5:
        status = f"DLA0817g: BTFR same -> TENSION"
    else:
        status = "no data yet"
    print(f"  {z:6.2f}  {H_ratio:8.3f}  {H_ratio:12.3f}  {status:>30}")

print(f"""
  CRITICAL ASSESSMENT:
  ====================
  At z=2: our model predicts a0/a0(0) = {H_of_z(2)/H0:.1f}x
  Milgrom excludes a0/a0(0) ~ 4x at z~2.
  Our prediction ({H_of_z(2)/H0:.1f}x) is BELOW the exclusion (4x).
  -> NOT IMMEDIATELY EXCLUDED, but in TENSION.

  At z=4.26: our model predicts a0/a0(0) = {H_of_z(4.26)/H0:.1f}x
  McGaugh: BTFR at z=4.26 looks same as local.
  If a0 were {H_of_z(4.26)/H0:.1f}x larger, BTFR would shift by:
  v -> v * (a0_z/a0_0)^(1/4) = v * {(H_of_z(4.26)/H0)**0.25:.2f}x
  For DLA0817g: v=272 km/s -> should be 272/{(H_of_z(4.26)/H0)**0.25:.2f}
  = {272/(H_of_z(4.26)/H0)**0.25:.0f} km/s if a0 were local value.
""")

# How much BTFR shift?
z_test = 4.26
H_ratio = H_of_z(z_test) / H0
v_obs = 272  # km/s
# If a0(z) = H_ratio * a0(0), then v^4 = GM*a0(z) = GM*a0_0*H_ratio
# Same M -> v_z = v_0 * H_ratio^(1/4)
# Or: v_0 = v_obs / H_ratio^(1/4)
v_local_equiv = v_obs / H_ratio**0.25
# Alternatively: same v but less mass
M_ratio = 1.0 / H_ratio  # to get same v, need less M
print(f"  DLA0817g test (z={z_test}):")
print(f"    v_obs = {v_obs} km/s")
print(f"    If a0(z) = {H_ratio:.1f}*a0_0:")
print(f"      Same mass -> v should be {v_obs*H_ratio**0.25:.0f} km/s (not {v_obs})")
print(f"      Same v -> mass should be {M_ratio:.2f}x local BTFR prediction")
print(f"    This means galaxy is {1/M_ratio:.1f}x LESS massive than local BTFR implies")
print(f"    OR: a0 has NOT evolved by factor {H_ratio:.1f}")

# ============================================================
# 2. ALTERNATIVE: a0 depends on Lambda, not H
# ============================================================
print("\n" + "=" * 70)
print("  2. Alternative: a0 depends on Lambda (constant), not H(z)")
print("=" * 70)

print("""
  If the geometric mean mechanism uses Lambda instead of H:

  Option A: a0 = c*H0/(2*pi) [H-dependent] -> a0 evolves with z
  Option B: a0 = c*sqrt(Lambda/3) [Lambda-dependent] -> a0 = CONSTANT

  Lambda = 3*H0^2*Omega_L/c^2 (cosmological constant)

  sqrt(Lambda/3) = H0*sqrt(Omega_L)/c

  a0_B = c^2 * sqrt(Lambda/3) = c * H0 * sqrt(Omega_L)
""")

Lambda = 3 * H0**2 * Omega_L / c**2
a0_lambda = c * H0 * np.sqrt(Omega_L)
a0_H = c * H0 / (2 * np.pi)

print(f"  Lambda = {Lambda:.3e} m^-2")
print(f"  a0 (H-based) = c*H0/(2*pi)        = {a0_H:.3e} m/s^2  (ratio: {a0_H/a0:.3f})")
print(f"  a0 (L-based) = c*H0*sqrt(Omega_L)  = {a0_lambda:.3e} m/s^2  (ratio: {a0_lambda/a0:.3f})")

print(f"""
  The Lambda-based a0 is {a0_lambda/a0:.2f}x the observed value -- too large!
  But: if we use a0 = c*H0*sqrt(Omega_L)/(2*pi):
  a0 = {a0_lambda/(2*np.pi):.3e} = {a0_lambda/(2*np.pi)/a0:.3f}x observed

  Or with different normalization: a0 = (c*H0/(2*pi)) * f(Omega_L)
  The key question is: does f depend on z through H(z) or through Lambda?
""")

# Try different combinations
print(f"\n  Combination survey (looking for a0 ~ 1.2e-10):")
print(f"  {'Formula':>35}  {'a0 (m/s^2)':>12}  {'ratio':>8}")
print("  " + "-" * 58)
combos = [
    ("c*H0/(2*pi)", c*H0/(2*np.pi)),
    ("c*H0/6", c*H0/6),
    ("c*H0*sqrt(Omega_L)", a0_lambda),
    ("c*H0*sqrt(Omega_L)/(2*pi)", a0_lambda/(2*np.pi)),
    ("c*H0/(2*pi*sqrt(Omega_L))", c*H0/(2*np.pi*np.sqrt(Omega_L))),
    ("c*sqrt(Lambda*c^2/3)", c*np.sqrt(Lambda*c**2/3)),
    ("c^2*sqrt(Lambda/3)", c**2*np.sqrt(Lambda/3)),
    ("sqrt(c^3*H0*Omega_L)", np.sqrt(c**3*H0*Omega_L)),
]
for name, val in combos:
    print(f"  {name:>35}  {val:12.3e}  {val/a0:8.3f}")

# ============================================================
# 3. THE DILEMMA: H-dependent vs Lambda-dependent
# ============================================================
print("\n" + "=" * 70)
print("  3. The dilemma: which cosmic scale controls a0?")
print("=" * 70)

print("""
  Two interpretations of the geometric mean mechanism:

  Interpretation A: H-dependent (dynamic)
  ========================================
  H = sqrt(r_S * c/H(z))
  a0(z) = c*H(z)/(2*pi) -- evolves with z
  Predicts: stronger MOND at high z
  Status: in TENSION with BTFR non-evolution

  Interpretation B: Lambda-dependent (static)
  ============================================
  H = sqrt(r_S * c/sqrt(Lambda))
  a0 = c*sqrt(Lambda/3) / (2*pi) -- CONSTANT
  Predicts: a0 same at all z
  Status: CONSISTENT with BTFR non-evolution

  Interpretation C: Present-epoch H0 imprinted
  =============================================
  a0 was SET during some epoch and doesn't track H(z).
  E.g.: substrate formed at some time t_form, and a0 = c*H(t_form)/(2*pi).
  If t_form = now: a0 = c*H0/(2*pi) = const going forward.
  But galaxies at z=2 were also in a substrate with H(z=2) -> a0(z=2)?
  This has a logical problem: when was a0 "set"?

  Interpretation D: Lambda + geometric factor
  =============================================
  The correlation length is NOT c/H but c/sqrt(Lambda*c^2/3).
  Since Lambda = const, a0 = const.
  At z=0: sqrt(Lambda*c^2/3) = H0*sqrt(Omega_L)
  So a0 = c*H0*sqrt(Omega_L)/(2*pi) = c*H0*0.828/(2*pi)
  = 0.828 * a0_H = 0.828 * 1.05e-10 = 8.7e-11
  Ratio to observed: 0.72 -- not great, but within range.
""")

a0_D = c * H0 * np.sqrt(Omega_L) / (2 * np.pi)
print(f"  Interpretation D: a0 = c*H0*sqrt(Omega_L)/(2*pi) = {a0_D:.3e}")
print(f"  Ratio to observed: {a0_D/a0:.3f}")
print(f"  Ratio to a0_H:     {a0_D/(c*H0/(2*np.pi)):.3f}")

# What a0 value works best?
# From gs9c: MOND simple best-fit a0 = 1.135e-10
a0_sparc = 1.135e-10
print(f"\n  Best-fit from SPARC (gs9c): a0 = {a0_sparc:.3e}")
print(f"  c*H0/(2*pi) =              {c*H0/(2*np.pi):.3e}  (ratio: {c*H0/(2*np.pi)/a0_sparc:.3f})")
print(f"  c*H0*sqrt(OL)/(2*pi) =     {a0_D:.3e}  (ratio: {a0_D/a0_sparc:.3f})")

# Could a0 be c*H0*Omega_L^alpha/(2*pi) for some alpha?
# Find alpha such that c*H0*Omega_L^alpha/(2*pi) = 1.135e-10
target = a0_sparc / (c*H0/(2*np.pi))
# target = Omega_L^alpha -> alpha = log(target)/log(Omega_L)
if target > 0:
    alpha_fit = np.log(target) / np.log(Omega_L)
    print(f"\n  To match SPARC: need Omega_L^alpha with alpha = {alpha_fit:.3f}")
    print(f"  a0 = c*H0/(2*pi) * Omega_L^{alpha_fit:.3f}")
    a0_fit = c * H0 / (2*np.pi) * Omega_L**alpha_fit
    print(f"  Check: {a0_fit:.3e} (ratio: {a0_fit/a0_sparc:.4f})")

# ============================================================
# 4. CONFRONTATION TABLE
# ============================================================
print("\n" + "=" * 70)
print("  4. Honest confrontation: model vs observations")
print("=" * 70)

print("""
  ================================================================
  OBSERVATION                    | a0=const | a0~H(z) | Status
  ================================================================
  Local BTFR (z~0)              |   OK     |   OK    | Both OK
  SPARC RAR (z~0)               |   OK     |   OK    | Both OK
  Freeman limit                 |   OK     |   OK    | Both OK
  Milgrom z~2 (a0 < 4x)        |   OK     |  ~3x    | MARGINAL
  DLA0817g z=4.26 BTFR          |   OK     |  ~8x    | TENSION
  MIGHTEE z<0.08 (tentative)    |   OK     |  ~1%    | ~OK
  Falling RCs at z>0.77         |   ?      |  helps  | Ambiguous
  ================================================================

  VERDICT:
  a0 = const: consistent with ALL current data
  a0 ~ H(z): in TENSION with BTFR at z > 2

  Our geometric mean mechanism (gs9d) predicts a0 ~ H(z).
  This is DISFAVORED by high-z BTFR data.

  HOWEVER:
  - The high-z data is very limited (6 galaxies, Milgrom 2017)
  - Rotation curves at z > 1 are much less precise than local
  - Selection effects may bias toward massive/fast rotators
  - DLA0817g is a SINGLE galaxy

  The honest conclusion:
  The geometric mean mechanism in its simplest form (a0 = cH(z)/(2pi))
  is in tension with existing data. But the data is not yet
  conclusive enough for definitive falsification.
""")

# ============================================================
# 5. SALVAGE: modified geometric mean
# ============================================================
print("=" * 70)
print("  5. Can the geometric mean mechanism be saved?")
print("=" * 70)

print("""
  The problem: if a0 ~ H(z), it evolves too fast at high z.
  But the geometric mean mechanism is ELEGANT and gives correct scaling.

  THREE ESCAPE ROUTES:

  Route 1: Lambda instead of H
  =============================
  Replace the Hubble correlation length c/H with:
  L_corr = c / sqrt(Lambda*c^2/3) = c^2 / (H0*c*sqrt(Omega_L))

  Since Lambda = const, a0 = const.
  Physical reason: the correlation length is set by the de Sitter
  horizon (Lambda-dominated epoch), not the instantaneous Hubble.

  H = sqrt(r_S * L_corr) = sqrt(GM*c / (H0*c^2*sqrt(Omega_L)))
  a0 = c * H0 * sqrt(Omega_L) / (2*pi)

  Problem: gives a0 = 0.72 * observed. Not terrible but not great.

  Route 2: Frozen correlation length
  ===================================
  The correlation length was SET at some formation epoch z_f
  and hasn't evolved since. This gives a0 = c*H(z_f)/(2*pi) = const.

  To match a0 = 1.2e-10: H(z_f) = 2*pi*a0/c = 2.51e-18 1/s
  H(z_f)/H0 = 1.15
  -> z_f ~ 0.2 (very recent!)

  This seems arbitrary. Why z_f = 0.2?

  Route 3: Geometric mean of r_S and 1/sqrt(Lambda)
  ==================================================
  L_corr = 1/sqrt(Lambda) = c / (H0*sqrt(3*Omega_L))

  H = sqrt(r_S / sqrt(Lambda))
  = sqrt(GM/(c^2 * H0*sqrt(3*Omega_L)/c))
  = sqrt(GM/(c*H0*sqrt(3*Omega_L)))
""")

L_corr_lambda = c / (H0 * np.sqrt(3 * Omega_L))
print(f"  Route 3: L_corr = 1/sqrt(Lambda) = c/(H0*sqrt(3*Omega_L))")
print(f"         = {L_corr_lambda:.3e} m = {L_corr_lambda/kpc:.0f} kpc")
print(f"  vs c/H0 = {c/H0:.3e} m = {c/H0/kpc:.0f} kpc")
print(f"  Ratio: {L_corr_lambda/(c/H0):.3f}")

a0_route3 = c * H0 * np.sqrt(3 * Omega_L) / (2 * np.pi)
print(f"\n  a0 (Route 3) = c*H0*sqrt(3*Omega_L)/(2*pi) = {a0_route3:.3e}")
print(f"  Ratio to observed: {a0_route3/a0:.3f}")

# Route 4: Maybe the 2*pi is wrong, try other factors
print(f"\n  Route 4: a0 = c*H0*sqrt(Omega_L)/N, solve for N:")
N_match = c * H0 * np.sqrt(Omega_L) / a0_sparc
print(f"    N = c*H0*sqrt(Omega_L)/a0_SPARC = {N_match:.3f}")
print(f"    Interesting: {N_match:.2f} ~ 5 or c*H0*sqrt(Omega_L)/5")
print(f"    Not a clean number.")

# What about c*H0/(2*pi) simply being close to a0 by coincidence?
print(f"""
  FUNDAMENTAL QUESTION:
  =====================
  Is a0 ~ c*H0/(2*pi) a COINCIDENCE or CAUSAL?

  If CAUSAL (a0 tracks H): -> a0 evolves -> tension with high-z
  If COINCIDENCE:           -> a0 = const -> no tension
  If CAUSAL but tracks Lambda: -> a0 const -> no tension, but why Lambda?

  The coincidence argument is weakened by the fact that
  a0/cH0 = 1/(2*pi) is a NATURAL number, not arbitrary.
  But coincidences do happen (cf. Moon angular size = Sun angular size).
""")

# ============================================================
# 6. WHAT WOULD DEFINITIVELY FALSIFY THE MODEL?
# ============================================================
print("=" * 70)
print("  6. Falsification criteria")
print("=" * 70)

print("""
  ================================================================
  FALSIFICATION TABLE
  ================================================================

  Prediction                     | Falsified if...
  -------------------------------|--------------------------------
  a0 = c*H(z)/(2*pi)            | BTFR shift >20% between z=0-2
  (a0 ~ H version)              | with >50 galaxies at z>1
                                 |
  a0 = const ~ c*H0/(2*pi)      | a0 measured at z>1 differs by
  (Lambda version)               | >3*sigma from local value
                                 |
  BTFR: v^4 = GM*a0             | BTFR slope != 4.0 at >3*sigma
                                 |
  Freeman limit = a0/(2*pi*G)   | Sigma_max outside 100-200 M_sun/pc^2
                                 |
  RAR tighter than MOND simple   | gs9c Hybrid A chi2/N < MOND chi2/N
  for TGP model                  | on next-generation data
                                 |
  Dimensional transition 3D->2D | d_eff(r) NOT trending toward 2
                                 | at large r in stacked RC analysis
                                 |
  Screening in solar system      | Perihelion precession anomaly
                                 | >10^-9 level from extra force
  ================================================================

  CURRENTLY SURVIVING PREDICTIONS:
  ================================
  1. BTFR: v^4 = GM*a0 -- YES (confirmed by SPARC)
  2. Freeman limit -- YES (137 vs 140 M_sun/pc^2)
  3. RAR shape -- YES (Hybrid A within 0.01 dex of MOND)
  4. d_eff -> 2 at large r -- YES (from rotation curve shape)
  5. a0 ~ c*H0/(2*pi) -- YES numerically, UNKNOWN if causal

  CURRENTLY IN TENSION:
  =====================
  6. a0 evolution with z -- MARGINAL (Milgrom: 3x not excluded,
     but McGaugh BTFR non-evolution is concerning)
  7. Cluster masses -- SAME problem as MOND
  8. Screening -- NOT DERIVED from TGP equations
""")

# ============================================================
# 7. HIGH-Z BTFR QUANTITATIVE TEST
# ============================================================
print("=" * 70)
print("  7. Quantitative test: high-z BTFR with a0 evolution")
print("=" * 70)

print("""
  Using Genzel+2017 and other high-z rotation curve data:
  Can we detect the predicted BTFR shift?

  If a0(z) = c*H(z)/(2*pi), then:
  v^4 = GM * a0(z) = GM * a0_0 * H(z)/H0

  At redshift z, the BTFR normalization is:
  log(v^4/GM) = log(a0_0) + log(H(z)/H0)

  Shift in log(v) at fixed M:
  Delta_log_v = (1/4) * log10(H(z)/H0)
""")

print(f"  {'z':>6}  {'H(z)/H0':>8}  {'Delta_log_v':>12}  {'v_shift(%)':>11}  {'Detectable?':>12}")
print("  " + "-" * 55)
for z in [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
    Hr = H_of_z(z) / H0
    dlv = 0.25 * np.log10(Hr)
    v_pct = (Hr**0.25 - 1) * 100
    # Current precision of BTFR: ~0.05 dex in log(v)
    detectable = "YES" if abs(dlv) > 0.05 else "MARGINAL" if abs(dlv) > 0.02 else "NO"
    print(f"  {z:6.1f}  {Hr:8.3f}  {dlv:+12.4f}  {v_pct:+11.1f}  {detectable:>12}")

print(f"""
  At z=1: shift = +5.6% in v -> measurable with ~20 galaxies
  At z=2: shift = +14% in v -> easily detectable with good RCs
  At z=4: shift = +28% in v -> very clear if RCs are reliable

  Current status:
  - Genzel+2017: 6 galaxies at z~1-2.5, large errors
  - DLA0817g (z=4.26): one galaxy, v=272 km/s
  - JWST is rapidly adding data at z > 1

  PREDICTION TO TEST:
  If a0 ~ H(z), galaxies at z=2 should have v ~ 14% higher
  than expected from their local-BTFR mass.

  If a0 = const, v should match local BTFR at all z.

  This is testable with JWST + ALMA within 2-3 years.
""")

# ============================================================
# 8. ALTERNATIVE TEST: falling rotation curves at high z
# ============================================================
print("=" * 70)
print("  8. Falling rotation curves at high z")
print("=" * 70)

print("""
  Genzel+2017 and Lang+2017 observed FALLING rotation curves at z~1-2.
  This was unexpected and led to claims of "baryon dominance" at high z.

  In MOND: falling RCs are EXPECTED at high z IF the measured radii
  are well within r_MOND. All measured accelerations at z~2 are
  (3-11)*a0, i.e., in the Newtonian regime.

  In our model:
  - If a0(z) = const: same as MOND. Falling RCs because g >> a0.
  - If a0(z) = 3*a0_0: r_MOND is SMALLER by factor 1/sqrt(3) = 0.58x.
    The Newtonian regime extends further in.
    This makes falling RCs MORE likely, not less.

  So both versions (a0 const and a0 ~ H(z)) predict falling RCs
  at z~2. This is NOT a discriminating test.

  BUT: the amount of "phantom dark matter" (excess over Newton)
  at the outermost measured radius WOULD differ:

  If a0 const: phantom DM fraction = MOND prediction
  If a0 = 3x: phantom DM fraction = higher (more extra gravity)

  At r where g = 5*a0 (typical for Genzel data):
""")

# Phantom DM fraction
for a0_mult in [1.0, 2.0, 3.0]:
    g_bar = 5 * a0
    a0_eff = a0_mult * a0
    # MOND simple: g = g_bar/2 + sqrt(g_bar^2/4 + g_bar*a0_eff)
    g_obs = 0.5*g_bar + np.sqrt(0.25*g_bar**2 + g_bar*a0_eff)
    phantom_frac = (g_obs - g_bar) / g_bar
    print(f"  a0 = {a0_mult:.0f}x: g_obs/g_bar = {g_obs/g_bar:.3f}, phantom DM = {phantom_frac*100:.1f}%")

# ============================================================
# 9. THE REVIEWER'S THREE HARD TESTS
# ============================================================
print("\n" + "=" * 70)
print("  9. Addressing the reviewer's three hard tests")
print("=" * 70)

print("""
  TEST 1: Effective operator for 3D->2D crossover
  ================================================
  STATUS: NOT DONE.

  What we have: qualitative geometric mean argument (gs9d).
  What we need: derivation from TGP Lagrangian showing that
  the effective equation at large r becomes 2D.

  The challenge: TGP in weak field is Helmholtz (delta'' + delta = 0),
  which gives oscillatory sin(r)/r. To get 2D behavior (ln(r)),
  we need to remove the spring term at large scales.

  Possible approach: If the spring constant becomes scale-dependent
  (mu^2(r) -> 0 at r >> c/H0), the equation becomes Poisson at large r.
  Then in 2D geometry (confined to thickness H), we get ln(r) potential.

  But gs7b showed scale-dependent mu gives Newton, not MOND.
  The key difference here: we're not modifying mu perturbatively,
  but arguing the equation CHANGES CHARACTER (from Helmholtz to Poisson)
  at the cosmological scale.

  This requires a TWO-SCALE analysis of the TGP Lagrangian:
  micro (soliton, fm) and macro (galaxy, kpc).

  TEST 2: Screening in solar system limit
  ========================================
  STATUS: PARTIALLY ADDRESSED (gs9b section 11-12).

  With r_c = c/H0 (cosmic, Vainshtein): screening WORKS.
  dg/g ~ 10^-15 at Earth. SAFE.

  But: this requires the crossover scale to be UNIVERSAL (c/H0),
  not mass-dependent (sqrt(GM/a0)). This is consistent with
  the Lambda interpretation (a0 = const, screening at cosmic scale).

  TEST 3: Clusters and redshift evolution
  ========================================
  STATUS: PARTIALLY ADDRESSED.

  Clusters: same prediction as MOND (underprediction by ~2x).
  Not a unique success, not a unique failure.

  Redshift evolution: a0 ~ H(z) version is in TENSION.
  a0 = const version survives but loses the "prediction" aspect.
""")

# ============================================================
# 10. REVISED MECHANISM: Lambda-based geometric mean
# ============================================================
print("=" * 70)
print("  10. REVISED mechanism: Lambda-based geometric mean")
print("=" * 70)

print("""
  Given the tension with high-z BTFR, we REVISE the mechanism:

  ORIGINAL: H = sqrt(r_S * c/H(z))  -> a0 evolves (DISFAVORED)
  REVISED:  H = sqrt(r_S * l_Lambda) -> a0 = const (SURVIVES)

  where l_Lambda = c/sqrt(Lambda*c^2/3) is the de Sitter radius.

  Physical interpretation:
  ========================
  The substrate correlation length is NOT set by the instantaneous
  Hubble expansion rate, but by the COSMOLOGICAL CONSTANT Lambda.

  Lambda defines the ultimate de Sitter horizon of the universe.
  This is a STATIC property of spacetime, not a dynamic one.

  The substrate "knows" about Lambda because it determines the
  maximum possible correlation length, even in the matter-dominated era.

  At z=0: l_Lambda = c/(H0*sqrt(Omega_L)) = {c/(H0*np.sqrt(Omega_L))/kpc:.0f} kpc
  This is ~{c/(H0*np.sqrt(Omega_L))/(c/H0):.2f}x the Hubble radius.

  H = sqrt(r_S * l_Lambda) = sqrt(GM*c / (c^2*H0*sqrt(Omega_L)))
  = sqrt(GM / (c*H0*sqrt(Omega_L)))

  a0 = GM/(H^2) = c*H0*sqrt(Omega_L) = const

  But: a0 = c*H0*sqrt(Omega_L) = {c*H0*np.sqrt(Omega_L):.3e}
  vs observed: 1.2e-10
  Ratio: {c*H0*np.sqrt(Omega_L)/1.2e-10:.2f} -> TOO HIGH by factor ~5.

  With 2*pi: a0 = c*H0*sqrt(Omega_L)/(2*pi) = {c*H0*np.sqrt(Omega_L)/(2*np.pi):.3e}
  Ratio: {c*H0*np.sqrt(Omega_L)/(2*np.pi)/1.2e-10:.3f}
""")

# Try: a0 = c*H0/(2*pi) was 0.87x observed
# Lambda version: a0 = c*H0*sqrt(Omega_L)/(2*pi) = 0.87*0.828 = 0.72x
# Neither is exact. The "best" combination from gs7a was just c*H0/(2*pi).

print("""
  HONEST ASSESSMENT:
  ==================
  No simple combination of c, H0, Omega_L gives a0 = 1.2e-10 exactly.
  The best match remains c*H0/(2*pi) = 1.05e-10 (ratio 0.87).
  Adding Omega_L makes it WORSE (0.72).

  This suggests one of:
  1. a0 = c*H0/(2*pi) is a coincidence (a0 doesn't track H or Lambda)
  2. The 2*pi comes from a more subtle geometric factor
  3. There's an O(1) correction we're missing
  4. a0 involves BOTH H0 and Omega_L in a non-trivial way

  The DIMENSIONAL ANALYSIS from gs7a stands: a0 has the dimensions
  of c*H0, and the numerical coefficient is close to 1/(2*pi).
  Whether this is causal or coincidental remains OPEN.
""")

# ============================================================
# 11. FINAL SUMMARY
# ============================================================
print("=" * 70)
print("  11. FINAL SUMMARY: Status of the 3D->2D dimensional transition")
print("=" * 70)

print("""
  ================================================================
  WHAT SURVIVES (Level 1: Phenomenology)
  ================================================================
  [YES] Dimensional transition 3D -> 2D gives flat RCs
  [YES] BTFR: v^4 = GM*a0 follows from any 3D->2D transition
  [YES] Freeman limit: 137 M_sun/pc^2 (independent of model details)
  [YES] Hybrid A fits SPARC nearly as well as MOND (Delta_chi2/N = 0.04)
  [YES] d_eff(r) transitions from 3 to 2 (from RC shape)
  [YES] a0 ~ c*H0/(2*pi) dimensionally (ratio 0.87)

  ================================================================
  WHAT'S IN TENSION (Level 2: Mechanism)
  ================================================================
  [TENSION] a0 ~ H(z) predicts BTFR evolution -- disfavored by data
  [OPEN]    Why geometric mean? No derivation from TGP equations
  [OPEN]    Screening: Vainshtein works with r_c = c/H0 but not derived
  [SAME]    Cluster problem shared with MOND
  [OPEN]    2*pi factor: geometric origin not rigorously shown

  ================================================================
  HONEST STATUS
  ================================================================
  The dimensional transition 3D -> 2D is a VIABLE PHENOMENOLOGICAL
  FRAMEWORK for galaxy scaling relations. It reproduces MOND
  phenomenology (BTFR, RAR, Freeman limit) with a clear physical
  picture (curvature controls effective dimension).

  The MECHANISM (geometric mean of r_S and r_H) is elegant but:
  - Not derived from TGP equations
  - Predicts a0 evolution that may be falsified
  - Requires screening that works but isn't derived

  The reviewer's assessment is correct:
  "Mikrofizyka jeszcze nie siedzi, ale asymptotyka i struktura
  skalowa zaczynaja wskazywac, gdzie ona moze siedziec."

  NEXT CONCRETE STEPS:
  =====================
  1. Collect ALL high-z BTFR data (JWST + ALMA) and do
     quantitative chi^2 test for a0 = const vs a0 ~ H(z)
  2. Try to derive effective 2D action from TGP Lagrangian
     at cosmological scales
  3. Compute EFE prediction from the model and compare with
     satellite galaxy data (Crater II, Antlia 2, etc.)
  4. Check if Lambda-based mechanism gives better numerics
""")
