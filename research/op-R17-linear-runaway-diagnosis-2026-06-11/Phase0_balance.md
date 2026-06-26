---
title: "Phase 0 — pre-registration balance sheet: op-R17-linear-runaway-diagnosis"
type: phase0_balance
status: LOCKED
locked_date: 2026-06-11
cycle: op-R17-linear-runaway-diagnosis-2026-06-11
authorization: "User 2026-06-11: 'działaj z R1 #17, zobaczymy jakie będą wyniki' (Phase 0 + Phase 1 jointly authorized; FINAL ceremony pending user reaction)"
methodology_binding: "CALIBRATION_PROTOCOL §3.6 BINDING; TGP_NATIVE_COMPUTATIONAL_PATTERNS BINDING; anti-Lakatos LOCK inherited from γ-3+γ-3'+γ-5+γ-7+B+A+D+S07-INT+Mech-v-enum+P25 sequence"
anti_lakatos_lock: PRESERVED
---

# Phase 0 — pre-registration: R1 #17 runaway diagnosis

## §1 — Scope and decision structure

### §1.1 Primary question

Is the R1 #17 runaway (δ_present/δ_rec ≈ 6×10²¹³, γ-7 Phase 3 §3.1, LOCKED) a property of
**TGP cosmology itself**, or of the **specific transcription** used to obtain it?

### §1.2 The transcription under audit (inherited verbatim, LOCKED γ-7 Phase 3)

| Ingredient | Source | Form |
|---|---|---|
| Background kinematics | γ-3 LOCKED | a ∝ t, H = 1/t |
| Gravitational coupling | γ-5 LOCKED | G_eff = 6.674×10⁻¹¹ SI |
| Matter conservation | γ-7 Phase 3 assumption | M_univ = const = 10⁵³ kg → ρ̄ = 3M/(4πc³t³) ∝ t⁻³ |
| Growth equation | standard sub-Hubble Newtonian | δ̈ + 2Hδ̇ − 4πG_eff·ρ̄·δ = 0 |
| ε_G(t₀) | γ-7 Phase 3 LOCKED | 3G_eff·M_univ/(c³t₀) ≈ 1.71 |
| τ_init (recombination) | γ-7 Phase 3 LOCKED | 2.75×10⁻⁵ |

**Critical observation defining the audit:** the standard growth equation is DERIVED for a
background satisfying the Newtonian Friedmann pair {continuity: ρ̄̇ + 3Hρ̄ = 0; acceleration:
ä/a = −(4πG/3)ρ̄}. The transcription imposes a ∝ t (⇒ ä = 0) **kinematically** while keeping
ρ̄ ≠ 0 — whether these are jointly consistent is exactly falsifier F-R17-B.

### §1.3 Decision structure (pre-declared)

```
F-R17-A (regression gate) ──FAIL──► HALT (pipeline invalid; no verdicts)
        │PASS
F-R17-B (background self-consistency audit)
        ├─ CONSISTENT (Δ < 0.1 at all τ ∈ [τ_init, 1]) ──► GENUINE_PATHOLOGY (F-R17-C still run, informational)
        └─ INCONSISTENT_O1 (Δ ≥ 0.1 anywhere) ──► runaway = ARTIFACT_CANDIDATE ──►
F-R17-C (consistent-transcription growth, routes C1/C2a/C2b/C2c)
        ├─ ≥1 route in PASS band   ──► ARTIFACT_RESOLVED (conditional flag §3.4)
        ├─ best route in PARTIAL    ──► ARTIFACT_PARTIAL
        └─ no route ≥ PARTIAL       ──► ARTIFACT_OPEN
F-R17-D = aggregate (mechanical application of the above; no discretion)
```

## §2 — Mandatory inputs (read/used before any verdict)

1. γ-7 Phase 3: [[../op-CE-H-gamma-7-clumping-acceleration-2026-05-24/Phase3_derivation.md]] §3.1-3.3, §7.1 (R1 #17 origin) + Phase3_sympy.py lines 46-145 (LOCKED numerics)
2. Concept paper cosmology note: [[../../meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md]] (EQ-5 + §10.6 Q3)
3. STATE.md sesja #8 (γ-7 closure + R1 #17 CRITICAL registration)
4. CALIBRATION_PROTOCOL §3.6 (verdict discipline, PARTIAL semantics, §3.6.13 constants)

## §3 — Falsifiers (LOCKED 2026-06-11; thresholds immutable)

### F-R17-A — regression gate (mandatory precondition)

Re-implement the γ-7 Phase 3 LOCKED transcription independently and reproduce:
- ε_G(t₀) = 3G_eff·M_univ/(c³t₀) within ±2% of 1.71
- numerical ODE growth factor: |log₁₀(G) − 213.8| ≤ 1.0 (LOCKED reference 6×10²¹³)
- analytic asymptotic exponent structure δ ~ exp(−2√ε_G/√τ) confirmed symbolically

**FAIL → HALT** (no downstream verdicts; pipeline invalid).

### F-R17-B — background self-consistency audit (PRIMARY, decides classification axis)

Define residual of the acceleration equation under the audited transcription:

```
Δ(τ) ≡ |ä/a − (−(4πG_eff/3)·ρ̄(τ))| / H²(τ)        [dimensionless]
```

With a ∝ t (ä = 0): Δ(τ) = (4πG_eff/3)·ρ̄(τ)·t² /H²·… computed exactly in sympy (no hand
simplification trusted). Verdict classes (pre-declared):
- **CONSISTENT**: max Δ(τ) < 0.1 over τ ∈ [2.75×10⁻⁵, 1]
- **INCONSISTENT_O1**: max Δ(τ) ≥ 0.1

Additional mandatory sub-test (B2, analytic link): show symbolically whether the runaway
exponent 2√ε_G/√τ is generated exactly by the residual term (i.e., whether setting the residual
to zero removes the super-power-law mode). PASS/FAIL recorded; informs interpretation only.

### F-R17-C — consistent-transcription growth (QUANTITATIVE)

**Routes pre-declared and CLOSED (post-hoc additions = forbidden move #6).** Every route must:
(a) use ONLY pre-existing TGP elements (γ-3 H = 1/t; γ-5 G_eff; concept paper EQ-5/frontier creation —
pre-flagged at R1 #17 creation, γ-7 Phase 3 §3.3 #2/#4); (b) have a self-consistent background
(Δ-residual of its OWN background equations = 0 by construction or verified); (c) introduce 0 new
constants and 0 free exponents.

| Route | Background | Physical assumption (declared) |
|---|---|---|
| **C1** zero-active-mass | a ∝ t exact ⇒ effective source (ρ+3p)_eff = 0 | gravity term in growth eq. vanishes (standard R_h=ct discipline); drag 2H only |
| **C2a** frontier creation, unclustered, no momentum drag | S_creation = H·ρ̄ ⇒ M(t) ∝ t, ρ̄ ∝ t⁻² | created matter uniform (δ_S = 0), carries Hubble-flow momentum |
| **C2b** frontier creation, unclustered, with momentum dilution | as C2a | as C2a + created matter dilutes peculiar momentum (θ-drag term) |
| **C2c** frontier creation, comoving-clustered | as C2a | created matter inherits local δ (δ_S = δ) and local velocity |

**Derivation rule:** the linear growth ODE for each route MUST be derived symbolically in sympy
from the perturbed {continuity + Euler + Poisson} system with the declared source terms — no
hand-inserted coefficients. Limit check mandatory: S_creation → 0 must recover the standard
equation (else route FAILS structurally).

**Observable target (LOCKED, inherited γ-7 Phase 3 §3.2):** G_obs ≡ δ_present/δ_rec = 10³
(δ: 10⁻⁵ → 10⁻²). Bands per project convention (factor-10 primary, factor-100 PARTIAL — same
convention as γ-7/P25; NOT tuned to any anticipated result):
- **PASS_BAND**: log₁₀ G ∈ [2, 4]
- **PARTIAL_BAND**: log₁₀ G ∈ [1, 5] (excluding PASS band)
- **FAIL_LOW / FAIL_HIGH**: below/above PARTIAL band (direction recorded)

### F-R17-D — aggregate diagnosis (mechanical)

Per §1.3 decision tree. Pre-registered verdict names: GENUINE_PATHOLOGY / ARTIFACT_RESOLVED /
ARTIFACT_PARTIAL / ARTIFACT_OPEN / INDETERMINATE.

### §3.4 — Conditionality flag (mandatory in any ARTIFACT_* verdict)

S_creation = H·ρ̄ is concept paper §10.6 **open question Q3**, NOT a derived result. Any C2 route
landing in PASS/PARTIAL band therefore yields a verdict **CONDITIONAL on hyp-Q3**; the full
closure of R1 #17 would additionally require a future derivation cycle for S_creation
(candidate: op-frontier-creation-rate-derivation). Treating a conditional band-hit as
"TGP predicts structure formation" is forbidden move #8.

## §4 — Predecessor invariance LOCK

Regardless of outcome, ALL of the following are PRESERVED unchanged:
γ-3 B+ / γ-3' B+ / γ-5 B+ / **γ-7 HALT-B** / F8 FAIL_LITERAL ×4 / PR-017 / PR-018 / PR-019 /
PR-020 / P25 Phase 1 PARTIAL_SOURCE_NS_ONLY (cycle WIP-paused, separate user decision) /
γ classification (γ) OBSERVATIONAL_ANCHOR / Branch A mapping.

In particular: an ARTIFACT_* verdict does NOT reopen F8 (the γ-7 FAIL was about Ω_DE magnitude
via ξ_clump; this cycle's scope is the internal consistency of the growth equation only).
γ-7 Phase 3 SCENARIO B (empirical ξ_clump fallback) used the OBSERVED growth — its verdicts are
insensitive to this cycle's outcome by construction.

## §5 — Forbidden moves register (10)

1. Modifying any LOCKED predecessor verdict (§4 list)
2. Citing this cycle as F8 progress/rescue; any Ω_DE/z_onset/V_eff computation
3. Tuning S_creation normalization or time-dependence to fit G_obs (ONLY S = H·ρ̄ from §10.6 Q3 allowed; no free exponents, no S = αHρ̄ scans)
4. Introducing new constants (§3.6.13 budget: 0)
5. Press-Schechter / ΛCDM growth-function borrowing (inherited forbidden move γ-7 #18)
6. Post-hoc route additions beyond C1/C2a/C2b/C2c
7. Threshold loosening beyond §3 bands
8. Interpreting conditional band-hit as unconditional prediction (§3.4)
9. Hardcoded T_pass = True (cycle 1/2/7 discipline)
10. Using G_obs = 10³ as INPUT to any derivation (comparison target only)

## §6 — Constants classification (§3.6.13)

| # | Constant | Class | Source |
|---|---|---|---|
| 1 | G_eff = 6.674×10⁻¹¹ | LOCKED inherited | γ-5 |
| 2 | M_univ = 10⁵³ kg | LOCKED inherited (audited, not re-fit) | γ-7 Phase 3 |
| 3 | c = 3×10⁸ m/s | (α) fundamental | — |
| 4 | H₀ = 2.3×10⁻¹⁸ s⁻¹ (t₀ = 1/H₀) | (γ) OBSERVATIONAL_ANCHOR | per PR-019 classification |
| 5 | τ_init = 2.75×10⁻⁵ | LOCKED inherited | γ-7 Phase 3 |
| 6 | δ_rec = 10⁻⁵, δ_present = 10⁻² | observational comparison targets (NOT inputs) | γ-7 Phase 3 §3.2 |

**NEW constants: 0.** ✓

## §7 — Risk register (6)

| # | Risk | Severity | Mitigation |
|---|---|---|---|
| R-R17-1 | Sign/term error in source-modified perturbed fluid equations | HIGH | symbolic sympy derivation from first equations; mandatory S→0 limit check per route |
| R-R17-2 | Momentum treatment of created matter ambiguous | MEDIUM | three sub-cases C2a/C2b/C2c ALL pre-declared and ALL computed; no selection |
| R-R17-3 | "Consistency" criterion contestable (TGP ≠ GR) | MEDIUM | criterion = the SAME Newtonian framework γ-7 used to derive the growth eq.; internal-consistency audit is framework-neutral |
| R-R17-4 | Border-band outcome (log G near 2 or 4) | MEDIUM | bands pre-LOCKED §3; verdict mechanical; borderline distance reported honestly |
| R-R17-5 | Yukawa screening neglected (sub-Hubble unscreened, as γ-7) | LOW | assumption inherited verbatim from audited transcription; flagged, not modified |
| R-R17-6 | M_univ = 10⁵³ kg rough (O(2) uncertainty → shifts ε_G → shifts p) | MEDIUM | inherited LOCKED; sensitivity reported as INFORMATIONAL (no re-fit) |

## §8 — Anti-Lakatos checklist (Phase 0)

- ✓ R1 #17 cited as DIAGNOSIS TARGET (a registered flag), not as motivation to rescue any FAIL
- ✓ F8 FAILs NOT cited as motivation; Ω_DE out of scope (§5 #2)
- ✓ EQ-5/frontier-creation provenance: concept paper (pre-exists γ-7) + γ-7 Phase 3 §3.3 possible-resolutions #2/#4 (flagged AT R1 #17 creation — not post-hoc)
- ✓ Routes CLOSED set; selection rule mechanical (§1.3)
- ✓ Thresholds = project-convention bands, declared before execution
- ✓ Honest negatives pre-registered (GENUINE_PATHOLOGY, ARTIFACT_OPEN, INDETERMINATE)
- ✓ Conditionality of hyp-Q3 pre-flagged (§3.4)
- ✓ 0 new constants; no PR append under any outcome
- ✓ Predecessor invariance LOCK explicit (§4)

## §9 — Phase plan

- **Phase 1** (THIS session, authorized): F-R17-A/B/C/D full execution in sympy + numerical ODE; Phase1_derivation.md
- **Phase FINAL** (pending user reaction to Phase 1 report): closure ceremony, R1 #17 register
  update per pre-declared downgrade rules (README "Honest pre-disclosed outcomes"), STATE.md entry

## §10 — Anticipated outcomes (INFORMATIONAL ONLY — not pre-registered as verdict)

- F-R17-B: INCONSISTENT_O1 expected (ä = 0 vs −4πGρ̄/3 ≠ 0 appears irreconcilable at face value;
  residual expected to grow ∝ 1/τ toward early times — would explain why runaway concentrates at small τ)
- F-R17-C: routes expected to BRACKET observed growth (C1 no growth → FAIL_LOW; C2 power-law
  growth with O(1) exponents → anywhere in [PARTIAL_LOW … PARTIAL_HIGH]); whether ANY route lands
  in PASS band is genuinely open — bimodal between ARTIFACT_RESOLVED and ARTIFACT_PARTIAL
- Aggregate most likely: ARTIFACT_* family with conditional hyp-Q3 flag; GENUINE_PATHOLOGY would
  require the residual to vanish, which would be a surprising structural result in its own right
