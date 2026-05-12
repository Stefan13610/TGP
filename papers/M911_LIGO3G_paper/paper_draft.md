---
title: "M9.1'' canonical metric ansatz of TGP — observational falsification by GWTC-3 and recovery framework via emergent-metric substrate"
authors: ["Mateusz Serafin"]
date: 2026-05-07 (v1 draft) / 2026-05-10 (v2 revision: forecast → negative result)
journal_target: "Phys. Rev. D (or Class. Quantum Grav.)"
type: paper-draft
status: v2-draft-pending-figures
revision_history:
  - v1 (2026-05-07): "Strong-field test of M9.1'': 2PN-phase deviation predictions for ET/CE" — predictive forecast based on Phase 1 OOM heuristic β=-5/64 (G_SPA ≈ 1)
  - v1.1-update (2026-05-09): Phase 1.5 G_SPA=48 sympy-LOCK + GWTC-3 RE-RUN 5.02σ FALSIFIED → CRITICAL UPDATE banner added at top of v1; v1 historical content preserved as APPENDIX
  - v2 (2026-05-10): full revision per native-first methodology; reframed from "predictive forecast" to "negative result + factor-48 methodological finding + recovery via emergent-metric framework"
keywords: ["gravitational waves", "post-Newtonian expansion", "ppE", "GWTC-3", "modified gravity", "TGP", "scalar-tensor", "stationary phase approximation", "factor-48 G_SPA", "structural-modification theories", "negative result"]
related:
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
  - "[[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]]"
  - "[[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]"
  - "[[../../research/op-emergent-metric-from-interaction-2026-05-09]]"
  - "[[../../research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/ADDENDUM_2026-05-10_native_observables_first.md]]"
  - "[[../../meta/PPN_AS_PROJECTION.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
---

# M9.1'' canonical metric ansatz of TGP — observational falsification by GWTC-3 and recovery framework via emergent-metric substrate

> ## v2 REVISION 2026-05-10 — primary content (forecast → negative result)
>
> This is **draft v2** (2026-05-10, native-first revision per
> `meta/PPN_AS_PROJECTION.md` binding methodology). v1 framed
> M9.1'' as a *predictive forecast* for ET/CE 2035+; v2 reframes as
> *observational falsification* (already occurred in current GWTC-3
> data) plus *methodological finding* (factor-48 G_SPA correction)
> plus *recovery framework* (via emergent-metric substrate Phase 4
> zero-β region). v1 abstract + sections 1-6 + acknowledgments +
> references + submission notes are preserved as **APPENDIX A**
> (historical record of the predictive-forecast framing now superseded).
>
> **v2 sections (PRIMARY content below):**
>
> - §A2 Abstract v2 (negative result + recovery framing)
> - §1' Introduction (TGP program + M9.1'' canonical ansatz + falsification
>   via GWTC-3)
> - §2' The (5/6) U³ deviation + sympy-LOCKED Phase 1 + Phase 1.5
> - §3' Phase 1.5 G_SPA=48 sympy-LOCK + factor-48 methodological finding
> - §4' GWTC-3 RE-RUN 5.02σ falsification of M9.1'' specific f(ψ)=(4-3ψ)/ψ
> - §5' Recovery framework: emergent-metric Phase 4 zero-β region +
>   S07 alternative f(ψ) reset
> - §6' Discussion: native-first methodology + form-meaning case study +
>   meta-methodological lessons for SPA chain in structural-modification
>   theories
> - §7' Conclusion (negative result outcome + program-level path forward)

---

## §A2 Abstract v2 (native-first, post-Phase 1.5 + post-GWTC-3 RE-RUN)

We report the observational falsification of the canonical metric ansatz
M9.1'' of the Theory of Geometry of Pulsation (TGP) program, with
hyperbolic time-time component `g_tt = -c²·(4-3ψ)/ψ` (ψ being the
fundamental TGP scalar field). Through sympy-locked PN expansion of the
α=2 vacuum scalar EOM (op-newton-momentum/M9_1_pp_P1, sympy LOCK 5/5
PASS), M9.1'' produces native Taylor coefficients of `g_tt[Φ]`:
{c_1, c_2, c_3, c_4, c_5} = {+1, +1, **-5/6**, **-23/12**, **+337/72**},
matching GR exactly through 1PN (β_PPN = γ_PPN = 1) but deviating
structurally at U³ with magnitude (5/6). Mapping this native deviation to
the Yunes-Pretorius parametrized post-Einsteinian (ppE) chart at
2PN-phase (b_ppE = -1), the stationary phase approximation (SPA) chain
prefactor G_SPA derives exactly as **G_SPA = 48** (sympy-LOCK 5/5 PASS +
4-level verification: hand-calc + numerical sanity + alternative SPA
orthogonal route), giving β_ppE^TGP^(b=-1) = **-15/4 ≈ -3.75** for
equal-mass binaries. This is **factor ~48× larger** than the
small-perturbation approximation G_SPA ≈ 1 (Sampson-Yunes-Cornish 2013)
would suggest — a notable methodological finding: the SPA chain for
*structural-modification* theories (O(1) coupling, not perturbative)
amplifies through cross-terms in the e_n → α_n derivation. Re-running
the LIGO/Virgo TIGER-framework Bayes inference on GWTC-3 (~90 BBH events)
with this corrected β prior yields BF_TGP/GR = 3.5·10⁻⁶ (log10 BF =
-5.45, "OVERWHELMING GR preference"), corresponding to **5.02σ
observational falsification of M9.1'' specific f(ψ)=(4-3ψ)/ψ ansatz**
in current public data. Crucially, the falsification applies to the
specific M9.1'' functional form, NOT the TGP framework as a whole: the
native L1 substrate-emergent-metric machinery is preserved, and the
parallel cycle `op-emergent-metric-from-interaction-2026-05-09` Phase 4
provides a STRUCTURALLY DERIVED zero-β region in
(a_n, ξ_n, b_n, c_0, κ_σ) substrate parameter space (Path 1: c_0 = 0
substrate-scaling; Path 2: κ_σ = 4/3 canonical coupling) where
2PN-phase deviation vanishes, i.e. the framework can recover GR
agreement without sacrificing its native L1 emergent-metric structure.
We discuss methodological implications for ppE catalogs of structural-
modification theories and advocate for native-first three-layer
presentation (L1 native source / L2 projection chart / L3 falsifier)
in modified-gravity falsification reports.

## §1' Introduction (v2)

The advent of high-statistics gravitational-wave astronomy with the
LIGO-Virgo-KAGRA O3 run (GWTC-2 / GWTC-3 catalogs) has enabled
quantitative tests of post-Newtonian (PN) waveform predictions of
General Relativity (GR) at strong field [LVK 2021, 2023]. Within the
parametrized post-Einsteinian (ppE) framework of Yunes & Pretorius
[YP 2009], a generic phase deviation `δΨ(f) = β_ppE · u^b` (with
u = (πMf)^(1/3)) provides a unified language for testing modified-
gravity predictions. The b_ppE = -1 entry corresponds to 2PN-phase
modifications, sourced by U³ corrections in the static metric `g_tt`.

The Theory of Geometry of Pulsation (TGP) [Serafin 2025-2026] is a
covariant scalar-tensor program built on a single scalar field ψ over
a discrete substrate, with axion-like Z₂ symmetry and emergent
effective metric `g_eff[Φ]`. The TGP master metric ansatz, denoted
M9.1'' in the program ledger, takes the canonical hyperbolic form
[Eq. (1) in §A1 below]:

   ds² = -c²·(4-3ψ)/ψ · dt² + ψ/(4-3ψ)·δ_ij dx^i dx^j     (1)

Through sympy-locked PN expansion of the α=2 vacuum scalar EOM (∇²ε
+ 2(∇ε)²/(1+ε) = 0 with ε ≡ ψ-1, op-newton-momentum/M9_1_pp_P1
sympy LOCK 5/5 PASS), M9.1'' reproduces GR exactly at 0PN-1PN (β_PPN =
γ_PPN = 1) but predicts a structural deviation at U³ with magnitude
**(5/6)**:

   Δα_3_metric = α_3^TGP − α_3^GR = -7/3 − (-3/2) = **-5/6**     (2)

This deviation has been the core falsifiable prediction of the TGP
strong-field sector since 2026-05-04 (audyt T01_LIGO3G_falsifier).

**Purpose of this paper (v2 framing):** to report the **observational
falsification** of M9.1'' specific ansatz f(ψ) = (4-3ψ)/ψ by GWTC-3
(2019-2020 LIGO O3 catalog), made possible by **two cascading findings
discovered 2026-05-09:**

1. **Phase 1.5 G_SPA = 48 sympy-LOCK** (op-ppE-mapping Phase 1.5):
   the SPA chain prefactor mapping (5/6) U³ in `g_tt` to β_ppE^(b=-1)
   in the inspiral phase derives exactly as **G_SPA = 48**
   (sympy-exact, test-particle limit), NOT G_SPA ≈ 1 as Sampson-Yunes-
   Cornish 2013 [SYC 2013] small-perturbation framework would suggest.
   This is a **factor 48× upward correction** to the predicted ppE
   coefficient β_ppE^TGP^(b=-1) at η=1/4 from -5/64 ≈ -0.078 (Phase 1
   heuristic) to **-15/4 ≈ -3.75**.

2. **GWTC-3 RE-RUN 5.02σ falsification** (op-GWTC3-reanalysis Phase 2
   RERUN): re-running the LIGO TIGER-framework Bayes inference on
   GWTC-3 ~90 BBH events with the corrected β prior gives
   BF_TGP/GR = 3.5·10⁻⁶, log10 BF = -5.45 ("OVERWHELMING GR
   preference"), corresponding to **5.02σ observational falsification**
   of M9.1'' specific (4-3ψ)/ψ ansatz.

The paper is structured as follows:

- **§2'** establishes the (5/6) U³ deviation as native L1 prediction
  (sympy-LOCKED).
- **§3'** derives Phase 1.5 G_SPA = 48 and discusses why Phase 1's
  G_SPA ≈ 1 heuristic fails for structural-modification theories — a
  **notable methodological finding** independent of TGP fate.
- **§4'** reports the GWTC-3 RE-RUN 5.02σ falsification of M9.1''
  specific f(ψ).
- **§5'** outlines the recovery framework: emergent-metric Phase 4
  STRUCTURALLY DERIVED zero-β region in substrate parameter space, plus
  S07 alternative f(ψ) reset path.
- **§6'** discusses native-first methodology, form-meaning case study
  (Phase 1 heuristic as a "true formula structure, false numerical
  assumption" sibling to ψ.1.v1), and meta-methodological lessons.
- **§7'** concludes.

We adopt **phase-PN convention** throughout (Cutler-Flanagan 1994 /
Yunes-Pretorius 2009): the (5/6) U³ deviation in `g_tt` corresponds to
the b_ppE = -1 entry (2PN-phase) in the ppE parametrization. Sympy
verification scripts are reproducible; raw LOCK transcripts are linked
in op-ppE-mapping Phase 1.5 README.

## §2' The (5/6) U³ deviation — native L1 prediction (sympy-LOCKED)

Native L1 source: TGP scalar field equation in vacuum (single-Φ axiom):

   ∇²ε + 2(∇ε)²/(1+ε) = 0,    ε ≡ ψ − 1.        (3)

Asymptotic expansion ε(r) = a₁/r + a₂/r² + ... in dimensionless
PN parameter η_pn = a₁/r yields, via sympy LOCK 5/5 PASS
(op-newton-momentum/M9_1_pp_P1 §3.2):

   ε(η_pn) = η_pn − η_pn² + (5/3)η_pn³ − (10/3)η_pn⁴
            + (22/3)η_pn⁵ − (154/9)η_pn⁶ + ...   (4)

Substituting through f(1+ε(η_pn)) with M9.1'' f(ψ) = (4-3ψ)/ψ and
applying Newton-matching a₁ = U/2 yields the analytically-locked PN
expansion of `g_tt^TGP / (-c²)`:

| n | α_n^TGP (g_tt) | α_n^GR | Δα_n native |
|---|---|---|---|
| 0 | +1 | +1 | 0 |
| 1 | -2 | -2 | 0 (PPN γ=1) |
| 2 | +2 | +2 | 0 (PPN β=1) |
| **3** | **-7/3** | **-3/2** | **-5/6** ← LEADING TGP DEVIATION |
| 4 | +35/12 | +1 | +23/12 |
| 5 | -91/24 | -5/8 | -19/6 |
| 6 | +91/18 | +3/8 | +337/72 |

(GR Schwarzschild isotropic: g_tt^GR / (-c²) = (1 - U/2)²/(1 + U/2)² —
sympy verified.)

This is the **L1 native prediction** in the three-layer methodology
of `meta/PPN_AS_PROJECTION.md`: native Taylor coefs are *forced* by
α=2 + M9.1'' f(ψ) ansatz; **zero free parameters**. Falsification is
sharp: the native derivation makes a definite numerical prediction
which is either confirmed or ruled out at all orders {U³, U⁴, U⁵, ...}.

## §3' Phase 1.5 G_SPA = 48 sympy-LOCK + factor-48 methodological finding

### §3.1' SPA chain in M9.1''

The mapping from native L1 (metric coefficient Δα_3 = -5/6) to L2 chart
projection (β_ppE^TGP at b_ppE = -1) proceeds through the standard SPA
chain: (metric perturbation) → (modified orbital binding energy E_orb(v))
→ (modified energy flux F(v)) → (modified phase Ψ(f)).

In compact form (op-ppE-mapping Phase 1.5 §4):

   β_ppE^TGP^(b=-1) = -(3/(128 η)) · Δα_3 · G_SPA      (5)

where G_SPA encapsulates the SPA chain amplification from metric to
TaylorF2 phase coefficient α_4. **The crux of the v1 → v2 correction
lies in the value of G_SPA.**

### §3.2' Phase 1 heuristic (G_SPA ≈ 1, OBSOLETE)

The original v1 draft cited Sampson-Yunes-Cornish 2013 [SYC 2013] to
justify "G_SPA ≈ 1 for metric-only modifications", giving β_ppE^TGP =
-5/64 ≈ -7.81·10⁻² at η=1/4 and OOM window |β| ∈ [5.5·10⁻², 1.2·10⁻¹].
This placed M9.1'' at the *borderline* of LIGO-O3 sensitivity (bound
~10⁻¹) and *decisive* in the ET-D + CE 2035+ era — the predictive-
forecast framing of v1.

### §3.3' Phase 1.5 sympy-LOCK G_SPA = 48 (sympy-exact, test-particle)

A careful re-derivation of the SPA chain in M9.1'' isotropic test-
particle limit (op-ppE-mapping Phase 1.5 §2-§4) gives:

| LOCK | Quantity | TGP value | Verification |
|------|----------|-----------|--------------|
| L1 | f_TGP, h_TGP reproduce α_n^TGP (7/7 OK) | sympy LOCK 7/7 | Phase 1 |
| L2 | Δv²(U) at U³ = +5/4 (1PN matches GR exactly) | sympy LOCK | Phase 1.5 §2.1 |
| L3 | E_TGP(U) at U³ = 49/48; ΔE/m = 11/24; **Δe_2 = -4/3** | sympy LOCK + hand-calc 49/48 + numerical sanity at U=0.1 | Phase 1.5 §2.2-§2.4 |
| L4 | F_TGP(v) = F_GR(v) at leading 2PN-orbital (cross-channel via GW1+GW2) | by-construction PASS | Phase 1.5 §3 |
| L5 | **Δα_4 = -40 (sympy-exact); G_SPA = 48 (sympy-exact)** | sympy LOCK + alternative SPA orthogonal route | Phase 1.5 §4.3-§4.5 |

Specifically:

   Δe_2 / Δα_3 = (-4/3) / (-5/6) = 8/5     (metric → orbital amplification)
   Δα_4 = 30·Δe_2 + cross-terms = -40       (sympy-exact rational)
   G_SPA = Δα_4 / Δα_3 = -40 / (-5/6) = **48**     (sympy-exact)
   β_ppE^TGP^(b=-1) at η=1/4 = -(3/32)·40 = **-15/4 ≈ -3.75**

with η-correction ±25% from test-particle approximation, giving OOM
window |β_TGP| ∈ [2.81, 4.69].

### §3.4' Why Phase 1 G_SPA ≈ 1 fails: regime-of-validity error

The SYC 2013 framework [SYC 2013] derives "G_SPA ≈ 1 for metric-only
deviations" within a **small-perturbation regime**, where the metric
deviation is suppressed by a small dimensionless coupling parameter:

| Theory | Δg_tt scaling | Coupling bound |
|--------|---------------|----------------|
| Brans-Dicke | ∝ 1/ω_BD | ω_BD > 4·10⁴ (Cassini) |
| dCS gravity | ∝ ζ_dCS | ζ_dCS bounded |
| sGB gravity | ∝ ζ_sGB | similar |

In this regime, cross-terms in the SPA chain (`α_4 = 30·e_2 - 20·e_1·p_1
+ 10·p_1² - 10·p_2`) reduce to leading-order metric contribution, and
the prefactor reduces to ~1.

**TGP M9.1'' is fundamentally different:** f(ψ) = (4-3ψ)/ψ is a
**structural** modification with Δα_3 = -5/6 being O(1), NOT a small
perturbation. The cross-terms in the SPA chain do **not** reduce to
1 — instead they amplify through `α_4 = 30·e_2 + ...`, giving a
chain factor of **30·(8/5) = 48** between Δα_3 and Δα_4.

**Phase 1's heuristic G_SPA ≈ 1 was applied OUTSIDE its regime of
validity.** This is a notable methodological caveat for ppE catalog
work (Yunes-Yagi-Pretorius 2016 review etc.) — **structural-modification
theories require explicit SPA chain derivation** rather than the
small-perturbation approximation.

### §3.5' 4-level verification of G_SPA = 48

| Level | Method | Result |
|-------|--------|--------|
| L1 | Sympy LOCK 5/5 PASS (Phase 1.5 main derivation) | β = -15/4, G_SPA = 48 ✓ |
| L2 | Independent hand-calculation E_TGP(U³)/m = 49/48 from M9.1'' isotropic geodesic | Matches sympy LOCK ✓ |
| L3 | Numerical sanity: direct M9.1'' metric eval at U=0.1 vs series prediction | 6 decimal places match ✓ |
| L4 | Alternative SPA derivation orthogonal route (substitute concrete e_n, p_n FIRST, integrate dt/dv, extract β directly) | β = -15/4 EXACT match ✓ |

Two refutational possibilities ruled out (Phase 1.5 §10 sign-off):

- **(B) Methodological error in Phase 1.5:** REFUTED by L4 (alternative
  SPA gives identical β = -15/4).
- **(C) LIGO TIGER convention different than assumed:** REFUTED by
  WebSearch verification (TIGER δφ̂_n IS fractional deviation;
  conversion β_ppE = 4.336·δφ̂_4 at η=1/4 is correct).

**Possibility (A) — TGP M9.1'' truly observationally falsified at 5σ —
CONFIRMED.**

## §4' GWTC-3 RE-RUN: 5.02σ observational falsification

### §4.1' Re-running the TIGER framework Bayes with corrected β prior

The LIGO TIGER framework [LVK 2023, GWTC-3 ToGR] reports posteriors on
fractional PN-coefficient deviations δφ̂_n in the TaylorF2 inspiral
phase. The mapping to ppE convention at 2PN-phase (b_ppE = -1, η=1/4):

   β_ppE = (3/(128 η)) · Δα_4 = (3/32) · α_4_GR · δφ̂_4 ≈ 4.336 · δφ̂_4

GWTC-3 ~90 BBH events combined posterior at b=-1 (2PN-phase) gives
σ_β ≈ 0.78 (single-coef Bayes prior).

With Phase 1.5 LOCK β_ppE^TGP = -3.75 (η=1/4, ±25% test-p):

   Z-score: |β_TGP| / σ_β = 3.75 / 0.78 ≈ 4.81 (from posterior shift)
   Bayes factor (combined): BF_TGP/GR = **3.5·10⁻⁶**
   log10 BF: -5.45 ("OVERWHELMING GR preference" per Jeffreys scale)
   σ-level: **5.02σ FALSIFIED-OBSERVATIONAL**

(op-GWTC3-reanalysis Phase 2 RERUN, [[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]
for full posterior chains, Bayes factor calculation, and per-event
contribution table.)

### §4.2' Contrast with v1 Phase 1 GWTC-3 verdict (2026-05-07)

The v1 GWTC-3 reanalysis (op-GWTC3-reanalysis Phase 2 ORIGINAL,
2026-05-07) reported BF_TGP/GR ≈ 0.97 INCONCLUSIVE, deviation 0.37σ
from observed centroid — based on the Phase 1 heuristic β = -5/64.

The corrected RE-RUN gives BF = 3.5·10⁻⁶, σ = 5.02 — a **factor ~3·10⁵
shift in Bayes factor**, driven entirely by the factor-48 correction
to the predicted β. This dramatic change is itself a **methodological
finding**: heuristic SPA chain assumptions (G_SPA ≈ 1) outside their
regime of validity can shift falsification verdicts by ~5 OOM.

### §4.3' Scope of falsification

The 5.02σ verdict applies to **the specific M9.1'' ansatz f(ψ) =
(4-3ψ)/ψ**, not the TGP framework as a whole. Specifically:

| Element | Status post-2026-05-09 |
|---------|------------------------|
| Single-Φ axiom | PRESERVED (no new radiative DOF, GW1 c_T=c_s, GW2 3 DOF still hold) |
| α=2 vacuum EOM | PRESERVED (Φ-EOM is the L1 source in TGP formalism) |
| Substrate emergent metric machinery | PRESERVED (op-emergent-metric Phase 4 STRUCTURAL DERIVED) |
| **M9.1'' specific f(ψ) = (4-3ψ)/ψ** | **FALSIFIED-OBSERVATIONAL 5.02σ** |
| Multi-coefficient ratios M911-P2 | WITHDRAWN-needs-rederivation (Phase 1 ratios incorrect) |
| 4-channel orthogonal pattern M911-P3 | PARTIAL-FALSIFIED (2/4 channels invalid; BH5 + ε.1 remain LIVE) |

## §5' Recovery framework

### §5.1' Emergent-metric Phase 4 zero-β region (STRUCTURAL DERIVED)

The parallel research cycle `op-emergent-metric-from-interaction-2026-05-09`
Phase 4 derives the emergent effective metric `g_eff[Φ]` from substrate
2-point function in (a_n, ξ_n, b_n, c_0, κ_σ) parameter space. Phase 4
STRUCTURAL DERIVED [Phase4_falsifier_recovery.md] identifies a
**zero-β region** where the 2PN-phase ppE coefficient β_ppE^(b=-1)
identically vanishes, with two natural anchor paths:

- **Path 1: c_0 = 0** (substrate-scaling vanishes) → `g_eff[Φ]` reduces
  to GR-form at 2PN-phase order; β = 0 forced.
- **Path 2: κ_σ = 4/3 (canonical coupling)** + (a_n, ξ_n) tuned to
  preserve PPN γ=β=1 → β = 0 + native L1 c_n structure preserved.

Both paths are STRUCTURAL DERIVED (analytic), not phenomenological
fits. They demonstrate that the TGP substrate-emergent-metric machinery
**can recover GR agreement at 2PN-phase without sacrificing native L1
structure** — i.e., the framework is not falsified, only the M9.1''
specific functional form is.

### §5.2' S07 alternative f(ψ) reset path

Audit cycle `audyt/S07_M911_derivation` provides the parallel
mathematical path: derivation of `f(ψ)` from first-substrate principles
(rather than postulating M9.1'' (4-3ψ)/ψ form). Post-falsification of
(4-3ψ)/ψ specifically, S07 reset paths A (algebraic minimization) and
B (geometric-symmetry) become priority — alternative f(ψ) ansatze
preserving the conformal-budget identity f·h = 1 (forces PPN γ=β=1)
but with different higher-order coefficients (different c_3, c_4, …).

### §5.3' Cross-cycle convergence diagnostic

Native-first reframe identifies a **5-fold + 7-fold cross-cycle
convergence diagnostic** for the recovery framework
([[../../audyt/T01_LIGO3G_falsifier/ADDENDUM_2026-05-10_native_observables_first.md]]
§4):

| Cycle | Diagnosis layer | Convergence on separable sector structure |
|-------|------------------|-------------------------------------------|
| L01 ρ-bridge | numerical | Q3 ρ_EM_quantum/ρ_NS ~ 10⁻¹² |
| τ.3 substrate-clock | mechanism | TT10 tests L4, decoupled from ρ |
| ψ.1 substrate-light | operator class | dim-6 EFT disjoint from EM trace anomaly |
| Q2 vacuum budget | vacuum-level | matter sector NOT additive to bare Λ |
| Φ-vacuum scale | dual-V | substrate vs matter sector |
| **T01 (this paper)** | **chart-level** | β_ppE projection vs native c_n source |
| op-newton-momentum/B9 | regime-level | WEP composition test 6/6 PASS |

Seven independent diagnoses converge on the **separable sector structure**
of the TGP framework — a *structural* property, not post-hoc tuning.
This convergence is the framework-level argument that recovery (via S07
reset OR Phase 4 zero-β region) is mathematically natural, not ad hoc.

## §6' Discussion

### §6.1' Native-first three-layer methodology (binding 2026-05-10+)

Per `meta/PPN_AS_PROJECTION.md` §3.1, modified-gravity falsification
reports should adopt the three-layer presentation:

- **L1 native** = source-level observables of the theory (e.g., Taylor
  coefficients of `g_tt[Φ]`, scalar-field EOM).
- **L2 projection chart** = derived parametrizations via SPA / PPN /
  ppE / Bayesian framework (e.g., β_ppE, b_ppE, G_SPA).
- **L3 falsifier** = empirical bounds (e.g., GWTC-3 BF, σ-level,
  detector thresholds).

This paper's v2 framing makes the layering explicit. v1 conflated L1
(structural prediction) with L2 (preliminary chart projection) and
made forecasts on L3. The factor-48 correction is a **L2-level**
methodological event: native L1 c_n coefficients did not change between
Phase 1 and Phase 1.5; only the SPA chain projection from L1 to L2 was
corrected.

### §6.2' Form-meaning case study: Phase 1 heuristic = "true formula structure, false numerical assumption"

In the form-meaning lens of `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md`
§4, Phase 1's β_ppE^TGP = -5/64 (G_SPA ≈ 1 heuristic) was the same
**form** as Phase 1.5's β_ppE^TGP = -15/4 (G_SPA = 48 sympy-exact) —
the underlying formula `β_ppE = -(3/(128 η)) · Δα_3 · G_SPA` is
identical. The error was at the level of a **numerical assumption**
(specifically, the value of G_SPA) being applied outside its regime
of validity (SYC 2013 small-perturbation framework not applicable to
TGP structural O(1) modifications).

This is a **sibling case study to ψ.1.v1** ("BD-form Z(x)·F² varying-α
form, false 44-OOM hierarchy artifact in TGP context"): both involve
a **valid functional form being applied outside its regime of validity**
and producing an artifact (44 OOM in ψ.1.v1; factor-48 in T01 v1).
The form-meaning lens identifies these as **methodologically traceable
errors** (not random) — and resolution comes through native-first
re-derivation directly from L1 (here, sympy-LOCKED Phase 1.5 SPA chain
in M9.1''; in ψ.1.v3, dim-6 EFT operator class disjoint from EM
trace anomaly).

### §6.3' Meta-methodological lesson for ppE catalog work

ppE catalogs (Yunes-Yagi-Pretorius 2016, Tahura & Yagi 2018, ...) often
list theoretical predictions for β_ppE^(b=-1) from various modified-
gravity theories using approximate G_SPA mapping. **TGP M9.1'' is an
illustrative case where the SYC 2013 G_SPA ≈ 1 approximation fails by
~factor 50**, with major falsification implications. Catalog work for
**structural-modification theories** should explicitly derive G_SPA
in the relevant regime (e.g., test-particle exact + η-corrections)
rather than relying on small-perturbation heuristics.

### §6.4' Independent significance of factor-48 finding

Even setting aside TGP fate, the Phase 1.5 G_SPA = 48 derivation
constitutes a **notable methodological result**:

- It quantifies the regime-of-validity boundary of SYC 2013 framework.
- It demonstrates that structural O(1) modifications (Δα_3 ~ unity)
  amplify SPA chain prefactors by factors ~10²-10² depending on the
  specific cross-term coefficients.
- It advocates explicit sympy-LOCK 4-level verification (sympy +
  hand-calc + numerical + alternative orthogonal route) for high-stakes
  phenomenological predictions.

These methodological lessons are publishable independent of TGP fate
and should be incorporated into ppE catalog updates.

## §7' Conclusion (v2)

The canonical metric ansatz M9.1'' of TGP, with hyperbolic time-time
component f(ψ) = (4-3ψ)/ψ, is **observationally falsified at 5.02σ**
by GWTC-3 (LIGO-O3 ~90 BBH events) via TIGER-framework Bayes inference
at 2PN-phase (b_ppE = -1) with corrected Phase 1.5 G_SPA = 48 prior.
The falsification was made possible by two cascading findings:

1. **Phase 1.5 sympy-LOCK G_SPA = 48** (sympy-exact, test-particle,
   4-level verified) — factor 48× upward correction to the predicted
   β_ppE^TGP from -5/64 (Phase 1 heuristic) to -15/4 (Phase 1.5 LOCK).

2. **GWTC-3 RE-RUN 2026-05-09** with corrected β prior gives
   BF_TGP/GR = 3.5·10⁻⁶, log10 BF = -5.45, σ-level = 5.02σ
   (FALSIFIED-OBSERVATIONAL).

Crucially, the falsification applies to the **specific M9.1'' f(ψ)
form**, NOT the TGP framework as a whole. The native L1 substrate-
emergent-metric machinery is preserved (single-Φ axiom + α=2 vacuum
EOM + emergent g_eff[Φ] + GW1/GW2 cross-channel consistency), and
**recovery paths are STRUCTURALLY DERIVED**:

- **Path 1 (c_0 = 0):** substrate-scaling vanishes; `g_eff[Φ]` reduces
  to GR-form at 2PN-phase; β = 0 forced.
- **Path 2 (κ_σ = 4/3 canonical coupling):** β = 0 with native L1 c_n
  structure preserved.

Both paths are documented in op-emergent-metric Phase 4 STRUCTURAL DERIVED
([Phase4_falsifier_recovery.md]). Plus the **S07 alternative f(ψ)
reset** path provides a parallel mathematical route via first-principles
derivation of f(ψ) from substrate symmetry.

**Methodologically**, the factor-48 G_SPA finding is a **publishable
result independent of TGP fate**: the SPA chain prefactor for
*structural-modification* theories does NOT obey the small-perturbation
G_SPA ≈ 1 approximation (SYC 2013) — a notable caveat for ppE catalog
work that should be incorporated in future modified-gravity tests of
GR.

**Native-first three-layer methodology** is advocated for modified-
gravity falsification reports: L1 native (source-level), L2 chart
(projection), L3 falsifier (empirical). Conflating layers (e.g.,
"chart projection makes future predictions" as in v1 framing) leads
to brittle forecasts; layered framing makes the regime-of-validity
of each step explicit.

The TGP program continues with: (i) Phase 1.6 equal-mass DJS 2-body
Lagrangian (5% precision lock on G_SPA(η)); (ii) M911-P2 multi-coefficient
re-derivation; (iii) S07 alternative f(ψ) ansatz exploration;
(iv) emergent-metric Phase 4 zero-β region empirical validation in
post-Path-1/2 substrate setups; (v) BH5 + ε.1 cross-channel orthogonal
predictions remain LIVE.

## §A1 Appendix A — v1 historical content (predictive forecast framing, 2026-05-07)

> *The following sections are the original v1 paper draft from
> 2026-05-07, preserved for revision history. They describe M9.1''
> as a "predictive forecast for ET/CE 2035+" based on Phase 1
> heuristic β_ppE^TGP = -5/64 (G_SPA ≈ 1). This framing was
> SUPERSEDED on 2026-05-09 by Phase 1.5 + GWTC-3 RE-RUN findings
> (factor 48× correction → 5.02σ FALSIFIED-OBSERVATIONAL). The v1
> content is retained as a record of the predictive forecast framing
> that v2 reframes. Do not cite v1 sections without v2 context.*

---

## Abstract

We present a falsifiable prediction of the Theory of Geometry of
Pulsation (TGP) program: the canonical metric ansatz M9.1'' (with
hyperbolic time-time component `g_tt = -c²(4-3ψ)/ψ`, where ψ is the
fundamental TGP scalar field) reproduces General Relativity exactly
through 1st post-Newtonian order (β_PPN = γ_PPN = 1) but deviates
explicitly at the 3PN-energy level (or, equivalently, 2PN-phase in
the inspiral waveform): `|Δg_tt|/c² = -(5/6) U³`, where U = GM/(rc²)
is the dimensionless Newtonian potential. We map this deviation
into the parametrized post-Einsteinian (ppE) framework, obtaining
a 2PN-phase coefficient `β_ppE^TGP^(b=-1) = -(3/(128 η)) · (5/6) ·
G_SPA` with sympy-locked rational structure (no fitting). For
equal-mass binaries (η = 1/4) and central SPA chain prefactor
G_SPA = 1 the prediction is `β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 ·
10⁻²`, within an OOM window |β| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹]. Multi-
coefficient ratios at higher PN orders (3PN-, 4PN-, 5PN-phase)
are forced to {-23/10, -38/23, +337/228} respectively, providing
a TGP-distinguishing signature with zero free parameters that
discriminates against dCS, sGB, Einstein-Æther, and Brans-Dicke
modifications. Fisher-matrix forecasts indicate (i) Cosmic Explorer
will detect the M9.1'' deviation at >5σ in a single loud BBH event
(M_tot = 30 M_⊙, d_L = 1 Gpc), (ii) Einstein Telescope will reach
the same significance with a stack of ~10 events, and (iii) LIGO-O5
A+ design will provide the first decisive bound from a stack of
~250 BBH events accumulated over ~2.5 years. We discuss the
window of testability 2027–2035, the compatibility with LIGO O3
existing GWTC-3 ppE constraints, and propose a re-analysis of GWTC-3
with TGP-specific β_ppE prior as a near-term confirmation channel.

## 1. Introduction

The advent of third-generation gravitational-wave detectors —
Einstein Telescope (Maggiore et al. 2020) and Cosmic Explorer (Reitze
et al. 2019) — will dramatically improve the precision of
post-Newtonian (PN) waveform tests of General Relativity (GR) at
strong field. Sensitivity gains of 10× over LIGO O3 are expected,
enabling detection of phase deviations at the 10⁻⁴ level in stacked
analyses (Chamberlain & Yunes 2017; LIGO/Virgo/KAGRA 2023).

The Theory of Geometry of Pulsation (TGP) is a covariant scalar-tensor
program that produces an emergent effective metric from a single
scalar field ψ on a discrete substrate, with axion-like Z₂ symmetry
(Serafin 2025; Serafin 2026). The TGP master metric ansatz, denoted
M9.1'' in the program ledger, takes the canonical form

   ds² = -c² (4-3ψ)/ψ · dt² + ψ/(4-3ψ) · δ_ij dx^i dx^j     (1)

(in isotropic coordinates), with the conformal-budget identity
f(ψ) · h(ψ) = 1 enforcing γ_PPN = 1 exactly through 1PN. The
TGP α=2 vacuum scalar EOM produces a definite asymptotic
expansion of the effective potential ε(r), and substitution into
the metric (1) yields the analytically-locked PN expansion

   g_tt^TGP/(-c²) = 1 − 2U + 2U² − (7/3)U³ + (35/12)U⁴
                   − (91/24)U⁵ + (91/18)U⁶ + ...   (2)

while the standard Schwarzschild GR result in isotropic coordinates is

   g_tt^GR/(-c²)  = 1 − 2U + 2U² − (3/2)U³ + (1)U⁴
                   − (5/8)U⁵ + (3/8)U⁶ + ...        (3)

The first non-trivial structural deviation appears at the U³ level
with magnitude **(5/6)**:

   Δα₃ = α₃^TGP − α₃^GR = -7/3 − (-3/2) = **−5/6**  (4)

This is the central prediction of this paper.

## 2. The (5/6) U³ deviation: derivation and convention

### 2.1 Derivation from α=2 vacuum

The TGP scalar field equation in vacuum reduces to

   ∇²ε + 2(∇ε)²/(1+ε) = 0         (5)

where ε ≡ ψ - 1. The asymptotic expansion ε(r) = a₁/r +
a₂/r² + ... in the local PN parameter η_pn = a₁/r yields,
through sympy LOCK 5/5 PASS (Serafin 2026, op-newton-momentum):

   ε(η_pn) = η_pn − η_pn² + (5/3)η_pn³ − (10/3)η_pn⁴
             + (22/3)η_pn⁵ − (154/9)η_pn⁶ + ...     (6)

with c_n = a_n/a₁^n satisfying a deterministic recursion derivable
from (5). The Newton-matching condition a₁ = U/(2r) (equivalently
η_pn = U/2) propagates through f(1+ε(η_pn)) to (2). The verification
chain f(1+ε) → expansion in U is independently sympy-locked
(7/7 OK in op-ppE-mapping Phase 1).

### 2.2 PN counting convention

The U³ coefficient in g_tt is variously called "3PN" (in the
energy/PPN convention used by Will 1971 for static-metric expansion)
and "2PN-phase" (in the Cutler-Flanagan 1994 convention adopted
universally in gravitational-wave waveform literature). For
clarity and registry consistency we adopt **phase-PN** as
primary in this paper: the (5/6) U³ deviation in g_tt corresponds
to the b_ppE = -1 entry in the Yunes-Pretorius (2009) ppE
parametrization. See Serafin 2026 (audyt/T01_LIGO3G_falsifier/
CONVENTION_DECISION.md) for the full mapping.

## 3. Mapping to ppE: locked β_ppE^TGP^(b=-1)

In the stationary phase approximation (SPA) for binary inspiral
(Cutler & Flanagan 1994), a deviation in g_tt at order U^N
generates a phase modification δΨ(f) = β_ppE · u^b where
u = (πMf)^(1/3) and b = 2N - 7. For U³ this gives b = -1 (2PN-phase).

The mapping from metric coefficient to ppE phase coefficient
proceeds through the standard SPA chain: (metric perturbation)
→ (modified orbital energy E_orb(v)) → (modified energy flux dE/dt)
→ (modified phase Ψ(f)). For metric-only deviations (no new
radiative degrees of freedom — consistent with TGP's single-Φ
axiom and GW1, GW2 entries in the TGP predictions registry that
verify the c_T = c_s and 3-DOF polarization observed at GW170817),
the SPA chain prefactor G_SPA is O(1) for typical equal-mass
binaries (Sampson, Yunes, Cornish 2013). We obtain

   β_ppE^TGP^(b=-1) = -(3/(128 η)) · (5/6) · G_SPA       (7)

For equal-mass binaries (η = 1/4) and central G_SPA = 1:

   β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 · 10⁻²             (8)

The OOM window for G_SPA ∈ [0.7, 1.5] (uncertainty from the
modified quadrupole formula in M9.1''):

   |β_ppE^TGP^(b=-1)| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹]      (9)

### 3.1 Multi-coefficient TGP-distinguishing signature

Higher-order ppE coefficients follow the same SPA mapping with the
respective Δα_n from (2)-(3):

   Δα_4 = +23/12,  Δα_5 = -19/6,  Δα_6 = +337/72.

The corresponding 3PN-, 4PN-, 5PN-phase ppE coefficients yield
fixed ratios:

   β_(3PN-phase)/β_(2PN-phase) = -23/10                   (10a)
   β_(4PN-phase)/β_(3PN-phase) = -38/23                   (10b)
   β_(5PN-phase)/β_(4PN-phase) = +337/228                 (10c)

These ratios are deterministic functions of α=2 and the hyperbolic
form of f(ψ); they involve **zero free parameters**. Existing
modified-gravity competitors with b = -1 (dCS gravity, scalar-Gauss-
Bonnet, Einstein-Æther, Brans-Dicke higher-PN) all require ≥1
free coupling and cannot reproduce (10) without fitting. This
provides a sharp **multi-coefficient TGP signature** for ET+CE
era simultaneous Bayes inference on multiple ppE coefficients.

## 4. Detection forecasts

We perform Fisher-matrix forecasts using the TaylorF2 + ppE
single-coefficient deviation waveform model (Buonanno et al. 2009),
with analytical fits to the design ASD curves of LIGO A+, ET-D
(Hild 2010), and Cosmic Explorer (Reitze 2019). The Fisher 1σ
uncertainty σ_β includes a degeneracy factor of 5 capturing the
realistic post-marginalization covariance with M_chirp, η, and
χ_eff (Yagi & Yunes 2016 review).

### 4.1 Single-event 5σ thresholds (calibrated)

After calibration to literature SNRs (Maggiore et al. 2020 for
ET-D, Reitze et al. 2019 for CE):

| Detector / scenario | β_5σ | β_TGP / β_5σ | Verdict |
|---------------------|--------|----------------|---------|
| LIGO-O3 (GW150914-like) | ~10⁻¹ | ~0.78 | borderline |
| LIGO-O5 single (loud BBH M=30, 200 Mpc) | ~3·10⁻² | ~2.6 | **YES first decisive single-event 2027+** |
| ET-D single (loud BBH M=30, 1 Gpc) | ~10⁻² | ~7.8 | YES (>20σ) |
| CE single (loud BBH M=30, 1 Gpc) | ~3·10⁻³ | ~26 | **YES (>50σ) decisive** |
| ET+CE network single | ~2·10⁻³ | ~40 | YES (>>50σ) |

### 4.2 Stacked thresholds

Stacking N independent events scales σ_β as 1/√N. First decisive
detection (β_5σ_stack = β_TGP_central):

   N(LIGO-O5)   ≈ 16 events     (~2 months at A+ rate ~100 BBH/yr)
   N(ET-D)      ≈ 1 event       (single loud BBH)
   N(CE)        ≈ 1 event       (single loud BBH)
   N(ET+CE)     ≈ 1 event       (single loud BBH)

### 4.3 Borderline status of LIGO O3

The LIGO O3 single-event ToGR bound at 2PN-phase (Abbott et al.
2021, GWTC-2) is ~10⁻¹, comparable to our prediction
β_ppE^TGP ≈ 8 · 10⁻². The TGP prediction is therefore at the
**boundary of current LIGO sensitivity** — neither falsified nor
fully constrained. A re-analysis of GWTC-3 with the TGP-specific
β_ppE = -5/64 prior (instead of generic ppE prior) should
produce a Bayes factor of order unity between TGP and GR, with
modest discriminating power. We estimate ~1-3σ tentative signal
from public GWTC-3 data is achievable; this re-analysis is
recommended as a near-term confirmation channel.

## 5. Discussion

### 5.1 Status of M9.1''

The deviation (5/6)U³ is rigorously derived from M9.1'' under the
following assumptions, all internally consistent: (i) α = 2 ODE
kinetic exponent (TGP F2 LOCKED); (ii) hyperbolic f(ψ) form
(M9_1_pp_setup, structurally selected by anti-podal budget +
β_PPN = γ_PPN = 1 + ψ = 4/3 horizon coincidence with R3 mass
spectrum at 4-digit precision); (iii) single-Φ axiom (TGP F1
LOCKED — no new radiative DOF).

We acknowledge that M9.1'' itself is currently postulated rather
than derived from first-substrate principles; this is the subject
of a parallel derivation cycle (Serafin 2026, audyt/S07_M911_derivation/).
However, falsifiability does not require derivation: even as a
postulate, M9.1'' makes the sharp prediction (8) and (10), which
ET+CE will test or rule out within the 2027-2035 window.

### 5.2 Window of testability

The most aggressive timeline:
- **2027** (LIGO-O5 first observing run, A+ design): single-event
  detection borderline; stack ~16 events for first decisive.
- **2030** (LIGO-O5 + LISA via cross-channel BH5 QNM ringdown):
  cross-channel orthogonal verification possible.
- **2035** (ET-D + CE first observing run): single-event decisive
  detection in CE; ET-D single borderline-YES, ET-D stack ~10
  events confirms.
- **2036+** (1 yr of ET+CE operation): >700σ confidence on
  multi-coefficient TGP signature M911-P2.

### 5.3 Falsification clause

M9.1'' is sfalsified iff at least one of the following is verified:
- LIGO-O5 with N ≥ 250 BBH events does NOT detect β_ppE^(b=-1)
  in the window |β_ppE^TGP| ∈ [5.5·10⁻², 1.2·10⁻¹] at 5σ;
- ET-D with stack 10 BBH events does NOT detect the same window;
- CE with single loud BBH event (M_chirp ≥ 30 M_⊙, SNR ≥ 150)
  does NOT detect β_ppE^(b=-1) > 5·10⁻² at 5σ.

A positive detection of the predicted window with simultaneous
verification of the multi-coefficient ratios (10) at any of LIGO-O5,
ET-D, CE, or ET+CE will promote M9.1'' from `LIVE` to
`DERIVED-OBSERVATIONAL` status.

## 6. Conclusion

The M9.1'' canonical metric of TGP makes a sharp, falsifiable,
zero-free-parameter prediction at 2PN-phase: β_ppE^TGP^(b=-1) =
-5/64 ≈ -7.81·10⁻², with a multi-coefficient ratio signature
{-23/10, -38/23, +337/228} that discriminates TGP from all known
modified-gravity competitors. Cosmic Explorer will detect this
deviation in a single loud BBH event at >5σ; Einstein Telescope
within ~10 events; LIGO-O5 A+ within ~250 events. This prediction
provides a clean public falsification contract for the strong-field
sector of TGP, with first decisive test possible by 2027 and
decisive confirmation by 2035.

## Acknowledgments

Sympy LOCK verifications: op-newton-momentum (M9_1_pp_P1, c_n=2..7
α=2 vacuum) and op-ppE-mapping (Phase 1, β_ppE^TGP and multi-
coefficient pattern). Fisher-matrix forecasts: op-LIGO-3G-deviation
(Phase 2-3). Audit framework: T01_LIGO3G_falsifier (NEEDS, FALSIFIER
STATEMENT, CONVENTION DECISION, FINDINGS).

## References

[1] Yunes, N., Pretorius, F., 2009, Phys. Rev. D 80:122003.
[2] Cutler, C., Flanagan, É., 1994, Phys. Rev. D 49:2658.
[3] Buonanno, A., Iyer, B. R., Ochsner, E., Pan, Y., Sathyaprakash,
    B. S., 2009, Phys. Rev. D 80:084043.
[4] Blanchet, L., 2014, Living Rev. Relativ. 17:2.
[5] Yunes, N., Yagi, K., Pretorius, F., 2016, Phys. Rev. D 94:084002.
[6] Sampson, L., Yunes, N., Cornish, N., 2013, Phys. Rev. D 88:064056.
[7] Will, C. M., 2014, Living Rev. Relativ. 17:4.
[8] Maggiore, M., et al., 2020, JCAP 03:050 (Einstein Telescope
    Science Case).
[9] Reitze, D., et al., 2019, Bull. Am. Astron. Soc. 51:035 (Cosmic
    Explorer).
[10] Chamberlain, K., Yunes, N., 2017, Phys. Rev. D 96:084039.
[11] Yagi, K., Yunes, N., 2016, arXiv:1602.04674 (review).
[12] LIGO/Virgo Collaboration, 2021, Phys. Rev. D 103:122002 (GWTC-2
     ToGR).
[13] LIGO/Virgo/KAGRA Collaboration, 2023, Phys. Rev. X 13:041039
     (GWTC-3 ToGR).
[14] Vallisneri, M., 2008, Phys. Rev. D 77:042001 (Fisher matrix
     pitfalls).
[15] Hild, S., et al., 2010, Class. Quantum Grav. 27:015003 (ET-D
     design).
[16] Damour, T., Jaranowski, P., Schäfer, G., 2014, Phys. Rev. D
     89:064058 (4PN ADM).
[17] Mishra, C. K., Iyer, B. R., Sundararajan, P. A., 2016, Phys.
     Rev. D 93:084054 (3PN waveform).
[18] Serafin, M., 2025-2026, TGP_FOUNDATIONS.md, sek08c_metryka_z_
     substratu.tex, sek08a_akcja_zunifikowana.tex, op-newton-
     momentum/M9_1_pp_setup.md, op-newton-momentum/M9_1_pp_P1_
     results.md, op-ppE-mapping/Phase{0,1,2,3}_*.md, op-LIGO-3G-
     deviation/Phase{0,1,2,3}_*.md, audyt/T01_LIGO3G_falsifier/.

---

## Submission notes

**Target journals (in preference order):**
1. **Phys. Rev. D** — strong-field falsifier paper; expected length
   ~12-18 pages PRD format with appendices. Section structure
   above is paper-ready; needs LaTeX porting and figure work.
2. **Class. Quantum Grav.** — secondary target; comparable scope.
3. **Phys. Rev. Lett.** — possible if reduced to 4 pages with
   strongest result (multi-coefficient signature + ET/CE forecast)
   as headline.

**Pre-submission tasks:**
1. Produce 3 figures: (i) M9.1'' vs GR g_tt expansion residuals;
   (ii) β_ppE^TGP vs detector sensitivity curves (Bayes prediction
   bands); (iii) multi-coefficient ratio test in {β_2PN,β_3PN}-plane
   showing TGP point + competitor exclusion regions.
2. Tighter G_SPA lock (Phase 1.5 future cycle, ~5% precision):
   reduces OOM uncertainty on β_TGP from ~30% to ~5%.
3. Production-grade Fisher matrix re-run with bilby/pycbc and
   official noise curves: replaces analytical ASD fits.
4. GWTC-3 reanalysis (Tier 5 side-cycle) — provides independent
   pre-publication confirmation if positive.
5. (Optional) Bayes-factor calculation TGP vs GR for selected
   GWTC-3 events.

**Estimated time to PRD-ready submission:** ~1-2 months from this
draft, contingent on figure work and Phase 1.5 G_SPA lock.

**Author contributions:** Mateusz Serafin (theory + program design).

**Funding & data availability:** TBD per author preference.

**Co-author candidates:** TBD per author network.
