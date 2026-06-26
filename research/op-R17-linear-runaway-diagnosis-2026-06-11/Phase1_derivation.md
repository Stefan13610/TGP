---
title: "Phase 1 — F-R17-A/B/C/D: runaway diagnosis — ARTIFACT_PARTIAL"
type: phase_result
status: PHASE1_COMPLETE
phase: 1
cycle: op-R17-linear-runaway-diagnosis-2026-06-11
created_date: 2026-06-11
authorization: "User 2026-06-11: 'działaj z R1 #17, zobaczymy jakie będą wyniki'"
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
sympy: "13/13 PASS; 0 hardcoded T_pass=True"
falsifiers_resolved: "F-R17-A PASS (regression gate) / F-R17-B INCONSISTENT_O1 / F-R17-C {C1 FAIL_LOW, C2a PARTIAL_BAND, C2b FAIL_LOW, C2c PARTIAL_BAND} / F-R17-D ARTIFACT_PARTIAL"
anti_lakatos_lock: PRESERVED
---

# Phase 1 — R1 #17 runaway diagnosis

## §0 — Verdict at a glance

**F-R17-D = ARTIFACT_PARTIAL** (mechanical application of Phase 0 §1.3 decision tree):

| Falsifier | Verdict | Key number |
|---|---|---|
| F-R17-A (regression gate) | **PASS** | ε_G = 1.7056 (ref 1.71); runaway reproduced log₁₀G = 214.09 (ref 213.78, tol ±1) |
| F-R17-B (background audit) | **INCONSISTENT_O1** | Δ(τ) = ε_G/(3τ): **0.57 today, 2.07×10⁴ at recombination** (threshold 0.1) |
| F-R17-B.2 (lemma) | **PASS (exact)** | φ′(τ) = √(3Δ(τ))/τ — runaway generated EXACTLY by the unbounded residual |
| F-R17-C | C1 FAIL_LOW · **C2a PARTIAL (log₁₀G = 1.40)** · C2b FAIL_LOW · **C2c PARTIAL (log₁₀G = 4.10)** | observed log₁₀G = 3.0 **bracketed** by C2a/C2c |
| F-R17-D | **ARTIFACT_PARTIAL** | B inconsistent + best routes PARTIAL, none in PASS band [2,4] |

**Plain-language summary:** the 10²¹³ runaway is **not a property of TGP cosmology** — it is the
integrated effect of an internally inconsistent transcription (M_univ = const forces a matter
source that violates, by up to 4 orders of magnitude at recombination, the very background
dynamics the growth equation presupposes). Self-consistent TGP-native transcriptions built on
EQ-5 frontier creation (M ∝ t) give clean **power-law** growth with O(1) exponents that bracket
the observed growth factor 10³ from both sides (10¹·⁴ and 10⁴·¹) — within ~1.6 OOM instead of
210 OOM — but **no pre-declared route lands inside the factor-10 PASS band**.

## §1 — F-R17-A: regression gate (PASS)

Independent re-implementation of the γ-7 Phase 3 LOCKED transcription:
- ε_G(t₀) = 3G_eff·M_univ/(c³t₀) = 1.7056 (rel dev 0.26% vs LOCKED 1.71) — FP1 ✓
- log-form ODE integration (Radau, rtol 10⁻¹⁰): log₁₀(δ_present/δ_rec) = **214.09** vs LOCKED
  reference 213.78 (6×10²¹³), within ±1.0 tolerance — FP2 ✓
  (scipy emits benign numerical-jacobian overflow RuntimeWarnings near τ_init; integration
  converged, `success=True`, result inside tolerance — verdict-irrelevant)
- symbolic WKB: dominant balance 2p−2 = −3 ⇒ p = −1/2; λ = ±2√ε_G — exactly the γ-7 Phase 3
  §3.1 analytic structure — FP3 ✓

## §2 — F-R17-B: the audited transcription is internally inconsistent (INCONSISTENT_O1)

The growth equation δ̈ + 2Hδ̇ − 4πGρ̄δ = 0 presupposes a background satisfying the Newtonian
Friedmann pair. The audited transcription imposes a ∝ t (⇒ ä = 0) while carrying
ρ̄ = 3M/(4πc³t³) with M = const. The residual (FP4, symbolic):

```
Δ(τ) ≡ |ä/a + (4πG/3)ρ̄| / H²  =  G·M/(c³·t₀·τ)  =  ε_G(t₀)/(3τ)
```

- Δ(1) = 0.569 — already O(1) **today**
- Δ(τ_rec) = 2.07×10⁴ — the matter source the transcription feeds into the perturbation
  equation would, at recombination, demand background dynamics **20 000× stronger** than the
  kinematically imposed a ∝ t allows

**Verdict: INCONSISTENT_O1** (threshold 0.1; LOCKED Phase 0 §3).

### §2.2 — The lemma (FP5, exact symbolic identity)

The WKB phase derivative of the runaway mode satisfies **φ′(τ) = √(3Δ(τ))/τ exactly**.
Consequences:
- Δ(τ) unbounded (∝ 1/τ) ⟺ super-power-law (runaway) mode exists
- Δ(τ) bounded ⟹ φ grows at most logarithmically ⟹ **pure power-law growth, no runaway**

The runaway is therefore **identically the integrated inconsistency** — not an independent
dynamical prediction of TGP. This is the cleanest possible artifact diagnosis: remove the
inconsistency and the runaway mode is structurally incapable of appearing.

## §3 — F-R17-C: self-consistent TGP-native routes (pre-declared, CLOSED set)

All route ODEs derived symbolically (FP6-FP10) from the perturbed {continuity + Euler + Poisson}
system with creation source S = αHρ̄ (α = 1 ⇒ M ∝ t per EQ-5/§10.6 hyp-Q3); the general ODE
reduces at α → 0 to the audited runaway equation (mandatory limit check FP6 ✓). All routes are
scale-invariant (equidimensional) ⇒ exact indicial power-law solutions; numerical cross-check
rel dev < 10⁻¹³ (FP11 ✓). G_obs never enters any derivation (FP12: free symbols = {ε_G} only ✓).

| Route | Growth ODE (τ-units) | p₊ | log₁₀G (τ_rec→1) | Verdict |
|---|---|---|---|---|
| C1 zero-active-mass | δ″ + (2/τ)δ′ = 0 | 0 | 0.00 | **FAIL_LOW** |
| C2a unclustered, no drag | δ″ + (3/τ)δ′ + (1−ε_G)/τ² δ = 0 | −1+√ε_G = 0.306 | 1.40 | **PARTIAL_BAND** |
| C2b unclustered, drag | δ″ + (4/τ)δ′ + (2−ε_G)/τ² δ = 0 | −0.102 | 0.00 | **FAIL_LOW** |
| C2c comoving-clustered | δ″ + (2/τ)δ′ − (ε_G/τ²) δ = 0 | (−1+√(1+4ε_G))/2 = 0.898 | 4.10 | **PARTIAL_BAND** |

Observed: log₁₀G_obs = 3.0 (LOCKED comparison target). PASS band [2,4]; PARTIAL [1,5].

**Structural picture:** under any frontier-creation transcription the growth is a **power law
δ ∝ τ^p with p = O(1)** — the qualitative pathology is gone. The quantitative outcome is decided
by the (currently underived) momentum/clustering treatment of frontier-created matter:
unclustered creation dilutes δ (C2a low side, 10¹·⁴), comoving-clustered creation does not
(C2c high side, 10⁴·¹). The observed 10³ sits between the two limiting treatments.

### §3.2 — Borderline note (honest)

C2c misses the PASS band by 0.10 dex (10⁴·¹⁰ vs band edge 10⁴·⁰); C2a by 0.60 dex. Bands were
LOCKED Phase 0 §3 before execution; verdicts applied mechanically. No threshold adjustment
(forbidden move #7).

### §3.3 — Sensitivity (INFORMATIONAL; R-R17-6; NOT a re-fit)

ε_G inherits M_univ = 10⁵³ kg (γ-7 LOCKED, O(2) rough). M_univ × 0.5 puts C2c at 10²·⁵ (inside
PASS band); M_univ × 2 puts C2a at 10³·⁹ (inside PASS band). Inverse note: log₁₀G = 3 would
require M_univ = 6.4×10⁵² kg (C2c) or 1.6×10⁵³ kg (C2a) — both within factor 1.6 of the LOCKED
rough value. **Reported for transparency only; adopting a re-fit M_univ is forbidden move #3/#10.
The PARTIAL verdicts stand on the LOCKED value.**

## §4 — F-R17-D aggregate: **ARTIFACT_PARTIAL**

Mechanical per Phase 0 §1.3: gate PASS + F-R17-B INCONSISTENT_O1 + best route PARTIAL_BAND
(no PASS_BAND) ⇒ **ARTIFACT_PARTIAL**.

**Pre-declared consequence (README "Honest pre-disclosed outcomes"):**
- R1 #17 downgrade **CRITICAL → HIGH**, re-scoped: ~~"TGP linear theory predicts unphysical
  runaway"~~ → "TGP-native structure formation theory is OPEN; consistent transcriptions give
  power-law growth bracketing observation within ~1.6 OOM; the discriminating unknown is the
  momentum/clustering treatment of frontier-created matter + derivation of S_creation (hyp-Q3)"
- All quantitative C2 results **CONDITIONAL on hyp-Q3** (S_creation ∝ H — concept paper §10.6
  open question Q3, NOT derived) — per Phase 0 §3.4 this conditionality is mandatory in the verdict

## §5 — What this does and does NOT change (predecessor invariance)

| Element | Status post Phase 1 |
|---|---|
| γ-7 HALT-B + F-γ-7-A/B/C/D + F8 FAIL ×4 | **UNCHANGED** (γ-7 Phase 3 SCENARIO B used *observed* growth — insensitive to this diagnosis by construction) |
| γ-3 / γ-3' / γ-5 B+ | **UNCHANGED** (background kinematics untouched; this cycle audited a *perturbation-level transcription*, not the γ-3 background claims) |
| PR-017/018/019/020 | UNCHANGED; **no PR append** (diagnostic cycle, per Phase 0) |
| R1 #17 | **DOWNGRADE CRITICAL → HIGH + re-scope** (pending FINAL ceremony) |
| ζ-cycle viability | **UNBLOCKED in principle** — cosmological perturbation work no longer faces a structural runaway, but requires the C2-class consistent transcription + hyp-Q3 resolution |
| New open item | **candidate cycle: op-frontier-creation-rate-derivation** (derive S_creation from substrate dynamics; discriminate C2a vs C2c momentum treatment) — would convert PARTIAL bracket into a sharp test |

## §6 — Anti-Lakatos verification (Phase 1): 10/10 COMPLIANT ✓

- ✓ 0/13 hardcoded T_pass; all verdicts computed
- ✓ Routes = pre-declared CLOSED set (C1/C2a/C2b/C2c); none added/removed post-hoc
- ✓ Bands LOCKED before execution; borderline C2c (0.10 dex) NOT promoted (§3.2)
- ✓ G_obs used only as comparison target (FP12 circularity guard: coefficients depend on ε_G only)
- ✓ S_creation = Hρ̄ only (no α-scans, no exponent tuning; forbidden move #3)
- ✓ 0 new constants; M_univ sensitivity reported INFORMATIONAL, not adopted (§3.3)
- ✓ EQ-5 provenance pre-exists R1 #17 (concept paper + γ-7 Phase 3 §3.3 #2/#4)
- ✓ F8/Ω_DE untouched (no acceleration claims anywhere in cycle)
- ✓ All predecessor verdicts PRESERVED (§5)
- ✓ Honest negatives delivered as registered: 2/4 routes FAIL_LOW; aggregate NOT inflated to RESOLVED

## §7 — Recommended next step

Per Phase 0 §9: **Phase FINAL closure** with claim_status `CLOSED-RESOLVED ARTIFACT_PARTIAL`
(R1 #17 → HIGH re-scoped; candidate follow-up cycle op-frontier-creation-rate-derivation
registered as proposal, NOT activated). Awaiting user reaction to this report.
