---
title: "Phase 0 — Balance sheet + 6/6 gate + scope (L07 zero-sum Z₂ derivation)"
date: 2026-05-16
parent: "[[./README.md]]"
type: phase-balance
phase: 0
status: 🟢 BALANCE OK — proceed to Phase 1 (honest partial expectation B+ pre-registered)
---

# Phase 0 — Balance sheet

## §1 — Wprowadzane wielkości (this cycle)

| Symbol | Status | Definicja | Źródło |
|---|---|---|---|
| `Δ(x)` | **DEFINED** | local Z₂ order parameter (signed): Δ = ⟨φ⟩ in chiral phase, or domain density difference | sek01_ontologia remark:zero-precyzacja ZS1; Phase 1 T2 |
| `Φ(x)` | **INHERITED** | spatial-density field; relacja Φ = (φ/φ_ref)²·Φ₀ z substratem | sek01_ontologia eq.Phi-from-phi; Phase 1 T5 |
| `δφ(x)` | **DEFINED** | fluctuation around symmetric/domain vacuum: δφ = φ - v gdzie v = ±φ_ref (vacuum value) | Phase 1 T6 |
| `δΦ(x)` | **DERIVED** | fluctuation of Φ around vacuum Φ₀: δΦ = Φ - Φ₀ | Phase 1 T6 (z δφ expansion) |
| `ZS1_LHS` | **TARGET** | ∫_Σ ⟨Δ(x)⟩√h d³x | Phase 1 T4 (derived = 0 z Z₂) |
| `ZS2_LHS` | **TARGET** | ∫_Σ ⟨δΦ(x)⟩√h d³x = linear_part + quadratic_part | Phase 1 T7 + T8 |
| `ZS1_status` | **DERIVED-TARGET** | ∈ {Z₂-identity, axiom-stands, falsified} | Phase 1 T4 verdict |
| `ZS2_status` | **DERIVED-TARGET** | ∈ {Z₂-identity, partially-derived+boundary, axiom-stands, falsified} | Phase 1 T9-T10 verdict |

## §2 — Inherited LOCKs (z poprzednich cykli i core)

| Inheritance | Status | Źródło |
|---|---|---|
| TGP_FOUNDATIONS §1: single field Φ z Z₂ substrate axiom | LIVE | TGP_FOUNDATIONS.md |
| Substrate Z₂ symmetry of H_Γ (φ → -φ) | LIVE | sek01_ontologia §83, eq.H-Gamma |
| Relacja Φ = (φ/φ_ref)²·Φ₀ | LIVE | sek01_ontologia eq.Phi-from-phi (chain in stw:Phi-properties) |
| Φ > 0 everywhere (positivity) | LIVE | sek02_pole §116-117 |
| T-Λ closure: γ/12 = M_Pl²·H₀² scale fixed | LIVE | closure_2026-04-26/Lambda_from_Phi0 |
| prop:Lambda-positive uses ZS2 dla ⟨U(φ_min)⟩_Σ > 0 | LIVE (DEPENDS ON THIS CYCLE) | sek05_ciemna_energia §240-293 |
| Q2 F1 substrate-vacuum decoupling | LIVE | op-Q2-vacuum-budget-2026-05-10 |
| ρ-T^μ_μ bridge full SM coverage | LIVE | op-L01-rho-stress-energy-bridge-2026-05-04 (N1-N5 all CLOSED) |

## §3 — 6/6 P-requirements gate (target)

| P# | Requirement | Verification source (Phase 1) |
|---|---|---|
| P1 | H_Γ Z₂-invariance under φ → -φ explicit | T1 sympy |
| P2 | Z₂-invariant ground state ⇒ ⟨Δ(x)⟩ = 0 pointwise OR globally | T3 sympy |
| P3 | ZS1 ∫⟨Δ⟩√h = 0 derived AS Z₂-tożsamość (Path A audit closure) | T4 sympy |
| P4 | Φ expansion δΦ = (Φ₀/v²)·(2v·δφ + (δφ)²) explicit linear+quadratic split | T6 sympy |
| P5 | ZS2 linear part vanishes IF δφ averages to zero (parallel ZS1 dla δφ) | T7 sympy |
| P6 | ZS2 quadratic part status documented (boundary condition character; non-zero) | T8-T10 sympy + declarative |

## §4 — Risks (R-flags)

| R# | Risk | Mitigation |
|---|---|---|
| R1 | Spontaneous Z₂ breaking ±v domains — ZS1 wymaga domain balance globally | T3 + T4 handle Z₂-symmetric superposition over domains explicit |
| R2 | Φ = (φ/φ_ref)²·Φ₀ jest Z₂-EVEN — ZS2 NIE jest pure Z₂-identity | T5 explicit; T8 verifies quadratic remainder ≥ 0 |
| R3 | Compatibility Φ > 0 (sek02) z ZS2 wymagającym ⟨δΦ⟩ = 0 | T9 explicit: ZS2 jest na poziomie integralu, NIE pointwise; δΦ ∈ R; Φ = Φ₀ + δΦ pozostaje > 0 jeśli max(\|δΦ\|) < Φ₀ |
| R4 | Higher-order δφ⁴/v⁴ corrections | Outside scope (extension cycle); leading O(δφ²/v²) handled in T8 |
| R5 | Cosmological FRW boundary conditions (closed Σ topology) — ZS2 quadratic remainder | T10 + T11 declarative: this jest boundary condition w cyklicznym FRW Σ; otwarte do extension nonlocal cycle (audit L07 Path D) |
| R6 | prop:Lambda-positive coupling dla deficyt φ < 0 zakaz (sek02 §116-117) | T9 clarification: "deficyt" w sek05 odnosi się do δΦ < 0, NIE Φ < 0; consistent |

## §5 — Scope (in/out)

**IN:**
- Z₂-symmetry H_Γ pod φ → -φ symbolic verification
- ZS1 (chiralna) Z₂-tożsamość derivation (Path A audit L07)
- ZS2 (przestrzenna) linear-quadratic decomposition + status documentation
- Compatibility check z Φ > 0 + prop:Lambda-positive
- Comparison z QCD chiral analog ⟨q̄γ⁵q⟩=0

**OUT:**
- ZS2 quadratic remainder full structural origin (likely Path D nonlocal foundations) — deferred
- Higher-order δφ/v corrections (Wilson coefficients) — extension cycle
- Time-dependent cosmological evolution of ZS2 in FRW — separate cycle
- Lagrange multiplier (Path B) explicit construction — alternative path, NIE this cycle
- φ_eff = φ - ⟨φ⟩_Σ redefinition (Path C) — alternative, NIE this cycle

## §6 — Native vs Standard physics split

Per CYCLE_KICKOFF_TEMPLATE §0.3:

- **Native part:** Substrate Z₂ symmetry (TGP-specific); relacja Φ = (φ/φ_ref)²·Φ₀
  (TGP S05 single-Φ z Z₂ substrate); ZS1/ZS2 jako native warunki dla
  prop:Lambda-positive (sek05 dark energy mechanism)
- **Standard part:** SSB framework Goldstone-Nambu (1960-61); operator-level
  expectation value calculus; Wess-Zumino consistency; chiral analog QCD
  ⟨q̄γ⁵q⟩=0
- **Mode:** Native-with-mapping — standard SSB tools applied to TGP-specific
  Z₂ substrate; L2 reduction (comparison z QCD chiral structure)
  jest OPTIONAL last stage

## §7 — Validator gate compatibility

Validator `tooling/validate_kickoff.py` checks (mental verification):
- `contract::L1_native::output_observable` non-empty ✓ (ZS1/ZS2 LHS + status)
- `contract::L1_native::measurement_instrument` non-empty ✓ (symbolic Z₂ operator calculus)
- `contract::L1_native::falsification_rule` pre-registered ✓ (with date 2026-05-16)
- `pre_registration_date` matches cycle date ✓
- 6/6 P-requirements declared ✓

**Status:** Phase 0 balance OK; proceed to Phase 1 symbolic derivation.

## §8 — Pre-registered honest partial expectation

**Reasoning for B+ pre-registration:**

ZS1 (chiralna): derivation from Z₂-invariance H_Γ + Z₂-invariant ground state jest
**clean structural argument** — bezpośrednia tożsamość operatorowa, identyczna
strukturalnie z QCD ⟨q̄γ⁵q⟩=0. **A− achievable dla ZS1.**

ZS2 (przestrzenna): wynika z relacji Φ = (φ/φ_ref)²·Φ₀ która jest Z₂-EVEN. Czysta
Z₂-tożsamość NIE działa dla ZS2 quadratic part. **Linear part** (ZS2 ~ ∫δφ√h) może
być zerowany przez Z₂-symmetric averaging (analog do ZS1), ale **quadratic part**
(ZS2 ~ ∫(δφ)²√h) jest **positive-semi-definite** i wymaga additional argument.

Reasonable candidate dla quadratic part: **boundary condition / gauge fixing** —
"mean Φ = Φ₀ across hypersurface" jako gauge-like condition fixing global zero-mode
of Φ. To jest standard physics technique (analog do flat-direction gauge fixing
w supersymmetric models). **B+ verdict appropriate** jeśli we cleanly identify
ZS2 quadratic remainder jako "boundary condition" rather than "Z₂ identity."

**A− verdict only if:** unexpectedly find structural Z₂-identity dla quadratic part
(low probability, but worth checking explicit).

**HALT-B verdict only if:** discover Z₂-invariance H_Γ doesn't imply ZS1 (extremely
unlikely; would falsify entire L07 Path A).

## Cross-references

- [[./README.md]] — kickoff contract
- [[./Phase1_sympy.py]] — symbolic derivation (next deliverable)
- [[../../audyt/L07_zero_sum_axiom/README.md]] — audit issue cycle addresses
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] ax:zero + remark zero-precyzacja
- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]] §240-293 prop:Lambda-positive
