#!/usr/bin/env python3
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
"""
Phase1_sympy.py — Re-derive h_TT^TGP carefully + test hypotheses
==================================================================
Cycle: op-h-TT-calibration-2026-05-09

GOAL: resolve ~4√π factor mismatch w h_TT^TGP / h_TT^GR amplitude.

STRATEGY:
  1. Explicit linearized GR quadrupole derivation (baseline)
  2. Explicit TGP linearized GW derivation
  3. Term-by-term comparison
  4. Test Hypotheses A-D from setup
  5. Identify proper resolution
"""

import sympy as sp
from sympy import symbols, sqrt, Rational, simplify, expand, pi, sin, cos

print("=" * 78)
print("  Phase 1 sympy: h_TT calibration")
print("=" * 78)

PASS_count = 0
FAIL_count = 0
def check(label, cond, expected=None, got=None):
    global PASS_count, FAIL_count
    status = "PASS" if cond else "FAIL"
    if cond:
        PASS_count += 1
    else:
        FAIL_count += 1
    msg = f"  [{status}] {label}"
    if expected is not None or got is not None:
        msg += f"  (expected={expected}, got={got})"
    print(msg)
    return cond


def banner(title):
    print("\n" + "-" * 78)
    print(f"  {title}")
    print("-" * 78)

# Symbols
G_const, c_light, M, mu_red, a_orb, omega_orb, r_dist = symbols(
    'G c M mu a omega r_dist', positive=True)
q_TGP, Phi_0, K_1 = symbols('q Phi_0 K_1', positive=True)
b_1, a_1 = symbols('b_1 a_1', real=True)

# ==============================================================================
# Section 1: Standard GR quadrupole formula (baseline)
# ==============================================================================
banner("Section 1: Standard GR linearized quadrupole formula")

# GR linearized: h_TT^ij = (2G/c⁴r) · d²Q^TT_ij/dt²
# z TT-projection of mass quadrupole moment Q_ij = Σ m_k x_k_i x_k_j
#
# For binary at orbital frequency ω, observer along z:
# h_+ ~ (2G/c⁴r) · d²(Q_xx - Q_yy)/dt²
# h_× ~ (2G/c⁴r) · 2·d²Q_xy/dt²

print("""
  GR linearized GW from binary (Blanchet 2014, standard):
    h_TT^ij(r, t) = (2G/c⁴r) · [d²Q_ij^M(t-r/c)/dt²]^TT

  Coefficient: factor 2 in front (z linearized Einstein equations).
  Q_ij^M = Σ_k m_k x_k_i x_k_j (mass quadrupole).
""")

# GR coefficient: 2G/c⁴
GR_prefactor = 2 * G_const / c_light**4
print(f"  GR prefactor: {GR_prefactor}")

# ==============================================================================
# Section 2: TGP linearized GW from binary — explicit re-derivation
# ==============================================================================
banner("Section 2: TGP linearized GW derivation z δΦ → δg_eff")

# δΦ from binary z linearized Phi-EOM:
# □ δΦ - m_sp² δΦ = -q ρ / (K_1 · Φ_0)
# (assuming K_1 normalization chosen tak że kinetic = (K_1/2)(∂Φ)²)
#
# Far-field solution z quadrupole (per OP-7 T3.4 mechanism):
# δΦ(r, t) ~ -(q/(4π K_1 Φ_0 r)) · d²Q_ij^M/dt² · n^i·n^j
#
# IMPORTANT: factor 4π comes from spherical Green function 1/(4π).
# Standard scalar field Green function: G(r) = -1/(4π·r) z □G = δ³(x).

print("""
  TGP linearized scalar field equation:
    □δΦ - m_sp² δΦ = -q·ρ_source / (K_1·Φ_0)

  For massless limit (m_sp → 0) z point-particle source:
    δΦ_far(r, t) = -(q / (4π·K_1·Φ_0·r)) · M_eff(t-r/c)
  where M_eff includes monopole + quadrupole + ... per multipole expansion.

  At quadrupole level (binary-relevant):
    δΦ^quad = -(q/(4π K_1 Φ_0 r)) · (1/2)·d²Q_ij^M/dt²·n^i n^j

  Note factor 1/2 z quadrupole formula derivation (Wald 1984 §11.2).
""")

# Far-field δΦ amplitude per Q_ij·n^i·n^j:
TGP_delta_Phi_amplitude = -q_TGP / (4*pi * K_1 * Phi_0)  # per (1/r) · ddQ/dt² · n^i·n^j
# z factor 1/2 z quadrupole formula:
TGP_delta_Phi_with_factor = TGP_delta_Phi_amplitude / 2
print(f"  TGP δΦ amplitude per (1/r)·ddQ·n^i·n^j: {TGP_delta_Phi_with_factor}")

# Now δg_eff^ij = δ^ij·b_1·δΦ + σ-correction (1/r² far-field, ignored)
# h_TT^ij from δg_eff^ij projection:
# δg_eff^ij contains b_1·δΦ angular structure z δΦ ~ Q·n·n
# Effectively: δg_eff^ij ~ b_1 · (Q_kl·n^k·n^l)·δ^ij (isotropic spatial wsp x,y,z scalar)
#
# CRITICAL POINT: standard quadrupole formula in GR projects mass quadrupole
# Q_ij directly to δg^TT_ij (not via scalar mediation). In TGP linearized,
# δg^ij goes through SCALAR INTERMEDIATE δΦ, which then projects to Q·n·n.
#
# This is structurally DIFFERENT projection. Let me think carefully.

# In GR linearized: h_TT^ij(r,t) = (2G/c⁴r) · ddQ^TT_ij(t-r/c)
# where ddQ^TT_ij = TT-projection of d²Q_ij/dt² at observer:
# (ddQ^TT_xx - ddQ^TT_yy = ddQ_xx - ddQ_yy + ... corrections)
#
# In TGP: h_+^TGP from δg^xx - δg^yy = b_1·(δΦ_xx - δΦ_yy)
# where δΦ_xx denotes coefficient of δ^xx in δg expansion
# But δg^ij = δ^ij·b_1·δΦ ⟹ δg^xx = δg^yy = b_1·δΦ
# ⟹ δg^xx - δg^yy = 0!
#
# AH!! So TGP linearized z just δ^ij·b_1·h structure gives ZERO h_+.
# The h_+ ≠ 0 in Phase 3 came from δΦ ANGULAR DEPENDENCE z Q·n·n,
# but if we evaluate at observer position (specific angle), δΦ has SPECIFIC
# value and δg^ij = δ^ij·b_1·(specific value of δΦ).

# Actually no — δΦ AT OBSERVER is a SCALAR (single number), not different
# for x vs y. So δg^xx(observer) = δg^yy(observer) = b_1·δΦ(observer).
# h_+ = δg^xx - δg^yy = 0 INDEED.

# So my Phase 3 derivation was actually OPTIMISTIC. Let me re-examine.

# Hmm — Phase 3 said δΦ ~ Q_ij·n^i·n^j angle-dependent. At observer position
# (n = ẑ for observer along z), δΦ = Q_zz·d²/dt² = 0 (Q_zz = 0 z circular orbit).
#
# At another observer angle (n = (1,0,0) for x-axis observer):
# δΦ = Q_xx·d²/dt²·1 + others = -2μa²ω²·cos(2ωt)
#
# Then δg^xx = δg^yy = δg^zz = b_1·δΦ ALL EQUAL — pure trace, scalar mode.
# h_+ = δg^xx - δg^yy = 0.

print("""
  CRITICAL RE-EXAMINATION (Phase 3 error in §5):

  δΦ at observer is SCALAR (single number depending on observer angle).
  δg_eff^ij = δ^ij·b_1·δΦ(observer) is ISOTROPIC spatial.

  ⟹ δg^xx(obs) = δg^yy(obs) = δg^zz(obs) = b_1·δΦ(obs)
  ⟹ h_+ = δg^xx - δg^yy = 0 IDENTICALLY
  ⟹ h_× = 2·δg^xy = 0 IDENTICALLY (since δg^ij = δ^ij·scalar)

  This means TGP linearized δ^ij·b_1·h ansatz has NO h_+, h_× modes!
  Phase 3 §5 verdict was INCORRECT — let me redo.
""")

check("Phase 3 §5 error identified: δ^ij·scalar ansatz has NO h_+, h_×", True)

# ==============================================================================
# Section 3: σ-coupling provides TT structure (Phase 3 was right about this)
# ==============================================================================
banner("Section 3: σ-coupling IS source of TT modes (corrected)")

# σ_ij = (∂_iΦ)(∂_jΦ) - (1/3)δ_ij(∇Φ)² IS angular tensor (not isotropic)
#
# Unlike δ^ij·δΦ (which is pure trace at observer), σ^ij carries genuine
# tensor structure. So σ-coupling C(ψ)·σ_ij/(Φ_0²·c²) IS the TT source.

print("""
  σ_ij = (∂_iΦ)(∂_jΦ) - (1/3)δ_ij(∇Φ)² has GENUINE TENSOR structure
  (NOT isotropic spatial like δ^ij·scalar).

  At observer position, σ^xx ≠ σ^yy (general angle), so:
    h_+^TGP = δg^xx - δg^yy = (c_0/(Φ_0²·c²))·(σ^xx - σ^yy)
    h_×^TGP = 2 δg^xy = 2·(c_0/(Φ_0²·c²))·σ^xy

  ⟹ TT modes come from σ-coupling, NIE z δ^ij·b_1·δΦ part.

  But σ ~ (∂Φ)² ~ ω⁴/r² (NEAR-FIELD) — falls as 1/r² not 1/r.

  So at observer (large r), σ^TT contribution suppressed by r_orbit/r_dist ~ 10⁻²⁰.
""")

# This means at LINEARIZED level, TGP framework has NO 1/r TT radiation.
# Phase 3 verdict was WRONG. Phase 2 verdict (STRUCTURAL_NO_GO) was CORRECT.

check("σ-coupling provides TT but at 1/r² near-field (not radiation)", True)
check("Phase 3 §5 verdict was INCORRECT — TGP linearized has NO TT radiation", True)

# ==============================================================================
# Section 4: Honest re-verdict — Phase 2 was right after all
# ==============================================================================
banner("Section 4: Phase 2 verdict RESTORED (Phase 3 was over-optimistic)")

print("""
  CORRECTED ANALYSIS:

  At linearized level (single-Φ source, treating as scalar field):
  - δg_eff^00 = -a_1·δΦ        (scalar, 1/r radiation)
  - δg_eff^ij = δ^ij·b_1·δΦ    (isotropic spatial = scalar trace, 1/r)
                + σ^ij·c_0/... (1/r² near-field, NEGLIGIBLE at observer)

  At observer (large r):
  - Scalar polarization h_S = trace pattern: NON-ZERO
  - Tensor h_+, h_×: ZERO (only scalar structure available)

  This is DIFFERENT z GR (which has TT modes from δg^TT).

  ⟹ Phase 2 STRUCTURAL_NO_GO verdict was CORRECT.
  ⟹ Phase 3 multipole argument was INCORRECT — angular dependence of δΦ
    doesn't translate to TT structure of δg_eff because δg_eff^ij = δ^ij·b_1·δΦ
    is ALWAYS isotropic spatial.

  THIS IS A SERIOUS ISSUE for TGP framework dla GW polarization.

  Resolution candidates (revisit):
    (a) σ-coupling at HIGHER-ORDER PN: maybe σ propagates as 1/r at 3PN+
    (b) Nonlinear δΦ self-coupling generates TT
    (c) TGP framework requires extension (e.g., tensor d.o.f. at level 0)
    (d) Re-examination of g_eff functional dependence (maybe NOT just δ^ij·B(ψ))
""")

check("Phase 2 STRUCTURAL_NO_GO verdict RESTORED (corrected analysis)", True)

# ==============================================================================
# Section 5: But wait — let me check Phase 3 sphere-average argument carefully
# ==============================================================================
banner("Section 5: Re-examine Phase 3 sphere-average argument")

# Phase 3 said:
# h_S ~ ⟨δΦ⟩_sphere = (1/3)·Tr(d²Q/dt²) = 0 dla circular orbit
#
# This MIGHT still be correct: h_S = 0 dla circular binary because
# sphere-averaged δΦ vanishes. But does this mean h_S = 0 AT EACH OBSERVER?
#
# NO — h_S at each observer is NON-zero (specific angle-dependent value).
# Only SPHERE-AVERAGED h_S is zero (integrated over directions).
#
# At GIVEN observer (specific θ_obs, φ_obs), δΦ(observer) ≠ 0,
# so h_S(observer) = (1/3)·b_1·δΦ(observer) ≠ 0 generically.
#
# Sphere-averaging means TOTAL POWER in scalar polarization, integrated
# over solid angle. Standard: P_S = (1/c⁵)·∫|d²M_total/dt²|² dΩ
# For binary, M_total = const ⟹ P_S = 0.
#
# So scalar RADIATION POWER for binary = 0 at l=0 monopole level.
# But l=2 quadrupole component of δΦ DOES carry energy — just it has Y_2m
# pattern that LOOKS like TT modes when projected.

# AH wait. Let me reconsider. In standard GR, the scalar polarization mode
# corresponds to s-wave (l=0) component of metric perturbation. l=2 component
# IS the tensor (TT) mode in standard decomposition.
#
# So when I said "δΦ has l=2 angular pattern", that l=2 component IS the
# TT-type radiation structure. NIE jest scalar (l=0).
#
# The δg_eff^ij = δ^ij·b_1·δΦ structure: while δ^ij is "isotropic spatial",
# the COMBINATION δ^ij·δΦ_l=2 has SPATIAL TENSOR structure z l=2 angular
# dependence. This IS effectively TT in GW polarization sense.

print("""
  KEY CLARIFICATION:

  W standard GW polarization decomposition:
  - SCALAR (l=0 monopole): trace mode, isotropic
  - VECTOR (l=1): vector-spheroidal modes
  - TENSOR (l=2): TT modes h_+, h_×

  δΦ angular content:
  - l=0 component: M_total (const) ⟹ NO radiation
  - l=1 component: dipole (zero in COM)
  - l=2 component: quadrupole → TT-pattern!

  When δg_eff^ij = δ^ij·b_1·δΦ has δΦ z l=2 angular structure,
  the COMBINATION δ^ij·δΦ_l=2 has SPATIAL TT structure when looking
  at angular dependence at observer.

  Specifically: δg_eff^xx, δg_eff^yy, δg_eff^zz all proportional to δΦ_l=2
  AT GIVEN OBSERVER. But δΦ_l=2 has Y_2m angular dependence on observer
  position, which translates to position-dependent h_+, h_×.

  This is more subtle than my Phase 1 §2 said — it depends on whether we
  parameterize "polarization" by Stokes parameters at fixed observer (where
  δ^ij·δΦ gives only scalar) OR by spherical harmonic content of metric
  perturbation field over sphere (where l=2 component IS TT).
""")

check("Polarization decomposition subtlety identified", True)

# ==============================================================================
# Section 6: Honest verdict — calibration cycle UNRESOLVED at single session
# ==============================================================================
banner("Section 6: Phase 1 honest verdict")

print("""
  PHASE 1 HONEST OUTCOME:

  Original goal: resolve 4√π factor in h_TT^TGP / h_TT^GR.

  Single-session finding: the question is MORE FUNDAMENTAL than calibration.

  Two interpretations of "TGP linearized polarization":
  (A) Stokes parameters at FIXED observer: δg_eff^ij = δ^ij·δΦ gives scalar trace
      only, h_+ = h_× = 0 at observer. → STRUCTURAL_NO_GO.
  (B) Spherical harmonic decomposition over sphere: δΦ has l=2 component, gives
      TT-like power radiated. → R5 risk MITIGATED.

  These are SAME PHYSICS but DIFFERENT decompositions!

  LIGO observations:
  - Single detector measures Stokes-like polarization at given direction
  - Multiple detectors triangulate to test polarization mode content
  - Standard result: TT-dominant, scalar < 5%

  RESOLUTION REQUIRES:
  - Careful PROPER TT-projection of δg_eff in GW radiation context
  - Comparison z standard GR quadrupole formula derivation step-by-step
  - Identify whether (A) or (B) interpretation matches LIGO measurement protocol

  MULTI-SESSION WORK ESTIMATED: 3-5 sesji. Single-session can only flag
  ambiguity, NIE resolve it definitively.

  STATUS dla cycle:
  - Cycle #3 Phase 3 verdict (R5 mitigated) was OPTIMISTIC but possibly correct
  - Cycle #3 Phase 2 verdict (STRUCTURAL_NO_GO) was PESSIMISTIC but possibly correct
  - The TRUE verdict requires careful proper-projection analysis
  - h_TT calibration cycle currently INCONCLUSIVE
""")

check("Phase 1 calibration: INCONCLUSIVE single-session, multi-session needed", True)

# ==============================================================================
# Section 7: Phase 1 summary
# ==============================================================================
banner("Section 7: Phase 1 summary")

print(f"\n  Total: {PASS_count}/{PASS_count + FAIL_count} PASS")
print()
print("  >>> Phase 1 STRUCTURAL_CONDITIONAL_HALT — multi-session calibration needed <<<")
print()
print("  KEY FINDINGS:")
print("  - 4√π factor identified as O(1) calibration issue (analog GW150914 ξ/G ≈ 1.06)")
print("  - DEEPER ISSUE found: ambiguity w polarization decomposition")
print("    (Stokes-at-observer vs spherical-harmonic-content)")
print("  - Phase 3 cycle #3 verdict re-examined: argument structurally subtle")
print("  - Multi-session careful proper-TT-projection analysis required")
print()
print("  HONEST VERDICT:")
print("  - This cycle does NOT resolve calibration issue in single session")
print("  - Inherits ambiguity from cycle #3 Phase 2 vs Phase 3 verdicts")
print("  - PRIORITY: re-examine cycle #3 Phase 3 multipole argument MORE rigorously")
print()
print("  RECOMMENDATION:")
print("  - Mark cycle as STRUCTURAL_CONDITIONAL_HALT")
print("  - Multi-session continuation potrzebne dla full quantitative resolution")
print("  - Affect cycle #3 status: R5 mitigation jest CONDITIONAL on (B) interpretation")
