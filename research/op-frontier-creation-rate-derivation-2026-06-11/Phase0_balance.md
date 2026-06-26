---
title: "Phase 0 — pre-registration: op-frontier-creation-rate-derivation"
type: phase0_balance
status: LOCKED
locked_date: 2026-06-11
cycle: op-frontier-creation-rate-derivation-2026-06-11
authorization: "User 2026-06-11: 'op-frontier-creation-rate-derivation' → cycle activation (R17 Phase_FINAL §3 registered proposal)"
methodology_binding: "CALIBRATION_PROTOCOL §3.6 BINDING; anti-Lakatos LOCK inherited (sequence … + P25 + R17)"
anti_lakatos_lock: PRESERVED
---

# Phase 0 — frontier creation rate: derivation cycle

## §1 — Scope and decision structure

### §1.1 Primary question

Can the three inputs left open by R17 (CLOSED-RESOLVED ARTIFACT_PARTIAL) be derived
TGP-internally, converting the C2 growth bracket into a **parameter-free prediction**?

1. **S_creation** (matter-sector creation rate; concept paper §10.6 hyp-Q3 "S ∝ H?")
2. **M_univ(t) relation** (removes the rough M_univ = 10⁵³ kg input)
3. **Bulk perturbation-equation form** (selects among R17 routes C2a/C2b/C2c)

### §1.2 Pre-declared derivation routes (CLOSED set)

**F-FCR-B (M_univ relation) — decides ε_G:**
- **B1 (zero-energy/horizon condition):** R = c·t ∧ R = 2G·M_tot/c² ⇒ M_tot = c³t/(2G)
  ⇒ ρ_tot = 3H²/(8πG) (critical density identity) ⇒ **ε_G(tot) = 3/2 EXACT**.
  Provenance: concept paper open task "Derive Schwarzschild R_s z critical density Ω → 1"
  (TGP_GENERATED_SPACE_COSMOLOGY, final roadmap) — pre-exists this cycle. Status to record
  honestly: DERIVED vs STRUCTURAL_POSTULATE (if the zero-energy condition is consistent and
  minimal but not forced by LOCKED machinery).
- **B2 (E2 equilibrium fraction):** concept paper NATIVE claim Ω_m ≈ 0.31 = "equilibrium density
  (E2 stability condition)" (status: 'do verification'; F5 pre-registered there, factor 2)
  ⇒ ε_G(matter) = (3/2)·Ω_m ≈ 0.465.
- ε_G = (3/2)·Ω_m is the EXACT skeleton; Ω_m ∈ {1 (B1: matter = total), 0.31 (B2)} pre-declared.
  **No other Ω values; no ξ-tuning in M = ξc³t/G (ξ = 1/2 canonical horizon coefficient only;
  ξ = 1 reported as INFORMATIONAL sensitivity only).**

**F-FCR-A (creation rate):**
- A1: from B-relation: S_matter/ρ̄ = Ṁ/M = H EXACT (M ∝ t) — resolves hyp-Q3 POSITIVELY
  conditional on B
- A2 (substrate side): EQ-5 stationarity of E2 vacuum (d⟨Φ⟩/dt = 0) ⇒ S_Φ = 3H⟨Φ⟩ ∝ H —
  consistent scaling; the Φ-substrate → matter-sector bridge is NOT specified in the concept
  paper → if so, declare the bridge as a gap (PARTIAL_concept_mismatch candidate), NOT improvise

**F-FCR-C (bulk equation form):** frontier creation is boundary-localized (EQ-5 note: "kreacja
tylko na granicy"). Pre-declared possible verdicts:
- BULK_CLEAN — boundary-localized creation leaves sub-horizon bulk equations source-free
  (= R17 C2c-form: δ″ + (2/τ)δ′ − (ε_G/τ²)δ = 0) — but requires a transport/homogenization
  statement absent so far
- BULK_UNIFORM_NODRAG (= C2a-form) / BULK_UNIFORM_DRAG (= C2b-form) — effective uniform creation
- **PARTIAL_concept_mismatch** — concept paper lacks the bulk-transport specification needed to
  select; record honestly; all three forms evaluated and reported as a matrix

**Native τ_init (mandatory item):** under γ-3 a ∝ t: 1 + z = t₀/t ⇒ τ_init = 1/(1 + z_rec),
z_rec = 1090 (observational anchor, NEW γ-class constant — declared §5). NOTE: γ-7/R17 inherited
τ_init = 2.75×10⁻⁵ = ΛCDM age-at-recombination — a borrowed mapping; flag as R1 candidate;
**NO retroactive modification of γ-7/R17 LOCKED verdicts** (forbidden move #3).

### §1.3 Falsifier F-FCR-D (aggregate; mechanical)

Prediction matrix: 3 C-forms × 2 B-values; growth factor G = (1/τ_init)^p where p = positive
indicial root. **Comparison target (LOCKED): log₁₀G_obs = 3.0**; bands (project convention,
inherited R17): **PASS [2,4] / PARTIAL [1,5] / FAIL outside** — applied per cell.

Aggregate verdicts (pre-registered):
- **PREDICTION_REALIZED** — A+B+C all DERIVED (no postulate/mismatch) AND the selected cell in
  PASS band → PR-022 LOCK candidate (genuine pre-registered prediction)
- **STRUCTURAL_CONDITIONAL** — skeleton derived (ε_G = (3/2)Ω_m; S ∝ H; native τ_init) but ≥1
  input remains postulate/claim/mismatch → NO PR-022; prediction recorded as conditional with
  explicit missing-piece list
- **FAIL_NEGATIVE** — all derivable cells outside PARTIAL band → frontier-creation pathway
  cannot reproduce observed growth even conditionally
- **INDETERMINATE** — honest stall

## §2 — Mandatory inputs

R17 full cycle (esp. Phase1_derivation §3 route ODEs — machinery reused verbatim);
TGP_GENERATED_SPACE_COSMOLOGY (EQ-5, §10.6 Q3, Ω_m E2 claim, roadmap task);
γ-3 LOCKED (a ∝ t); γ-5 LOCKED (G_eff); γ-7 Phase 3 (τ_init provenance audit).

## §3 — Forbidden moves (9)

1. ξ-tuning in M = ξc³t/G (B1 = 1/2 canonical only; ξ = 1 INFORMATIONAL)
2. Ω_m values beyond pre-declared {1, 0.31}; no Ω fitting to G_obs
3. Retroactive modification of γ-7/R17/any LOCKED verdict (τ_init flag is forward-looking only)
4. C-form cherry-picking (all 3 forms reported; selection only by derivation, else mismatch declared)
5. PR-022 append unless PREDICTION_REALIZED (all-derived + PASS band)
6. G_obs as input anywhere (comparison only; circularity guard mandatory)
7. New fundamental constants (budget 0; z_rec = observational anchor declared §5)
8. Hardcoded T_pass
9. Improvised Φ→matter bridge for A2 (gap declared, not filled ad hoc)

## §4 — Predecessor invariance LOCK

ALL PRESERVED: R17 ARTIFACT_PARTIAL; P25 NEGATIVE; γ-3/γ-3'/γ-5 B+; γ-7 HALT-B; F8 ×4;
PR-017..020; γ OBSERVATIONAL_ANCHOR; Branch A. F8/Ω_DE out of scope (growth ≠ acceleration).

## §5 — Constants (§3.6.13)

G_eff (γ-5 LOCKED); c; H₀ (γ-anchor); **z_rec = 1090 — NEW (γ) OBSERVATIONAL_ANCHOR** (CMB epoch;
used for τ_init mapping, not for any dynamics); δ targets 10⁻⁵→10⁻² (comparison only).
**New fundamental constants: 0.**

## §6 — Risk register

R-FCR-1 (HIGH): B1 zero-energy condition may be unforced → STRUCTURAL_POSTULATE status honest;
R-FCR-2 (HIGH): C-form selection likely PARTIAL_concept_mismatch (transport unspecified);
R-FCR-3 (MED): Ω_m = 0.31 concept-claim unverified (its own F5 pending); R-FCR-4 (LOW): τ_init
flag optics — guarded by forbidden move #3 + forward-only scope; R-FCR-5 (MED): borderline band
cells — mechanical application, distances reported.

## §7 — Anticipated outcomes (INFORMATIONAL, not pre-registered as verdict)

Most likely **STRUCTURAL_CONDITIONAL**: skeleton ε_G = (3/2)Ω_m EXACT + S ∝ H + native τ_init
derived; B1 likely postulate-status; C likely mismatch → matrix reported, conditional prediction
recorded, PR-022 NOT appended. PREDICTION_REALIZED possible only if concept-paper texts force
both B1 and a unique C-form — genuinely open.
