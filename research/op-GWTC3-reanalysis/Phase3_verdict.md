---
title: "Phase 3 — VERDICT: TGP M911-P1 vs GWTC-3 (consistent + reconciliation)"
date: 2026-05-07
parent: "[[README.md]]"
type: phase3-verdict
tgp_owner: research/op-GWTC3-reanalysis
tags:
  - phase3
  - verdict
  - synthesis
  - reconciliation
  - GWTC-3
  - TGP
  - M911
  - tier-5
related:
  - "[[README.md]]"
  - "[[Phase1_GWTC3_bounds.md]]"
  - "[[Phase2_Bayes_factor.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/FINDINGS.md]]"
  - "[[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]]"
---

# Phase 3 — VERDICT: TGP M911-P1 vs GWTC-3

## ⚠ CRITICAL UPDATE 2026-05-09 — VERDICT REVERSED

**This Phase 3 verdict (2026-05-07: "TGP CONSISTENT, BF≈0.97 INCONCLUSIVE")
was based on Phase 1 OOM heuristic β_TGP = -5/64. [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]]
Phase 1.5 derived G_SPA = 48 (4-level verified: sympy LOCK 5/5 + hand-calc
+ numerical + alternative SPA), giving β_TGP = -15/4 ≈ -3.75 — factor 48×
larger than Phase 1.**

**Phase 2 RE-RUN with corrected β** ([[Phase2_RERUN_2026-05-09_corrected_beta.md]]):
- **Verdict REVERSED:** TGP M9.1'' RULED OUT at 5.02σ (BF = 3.5·10⁻⁶,
  log10(BF) = -5.45 → "OVERWHELMING GR preference").
- δφ̂_4_TGP corrected = -0.865 (vs original -0.018), 5σ tension with
  GWTC-3 observed +0.05 ± 0.182.
- TGP M9.1'' **falsified-observational** under TIGER framework analysis.

**Sections below are HISTORICAL** (Phase 1 baseline). Current operational
verdict: M9.1'' specific (4-3ψ)/ψ ansatz is observationally excluded by
GWTC-3 at 5σ; alternative f(ψ) forms via S07 reset are next exploration
path.

---

## §0 — Executive verdict (HISTORICAL — pre Phase 1.5)

```
┌──────────────────────────────────────────────────────────────────┐
│  GWTC-3 reanalysis verdict for TGP M911-P1                       │
│                                                                  │
│  Status:   CONSISTENT z GWTC-3 within 1σ                         │
│            BF_TGP/GR ≈ 0.97 → INCONCLUSIVE                       │
│            Deviation: 0.37σ (none significant)                   │
│                                                                  │
│  Detection power: TGP signal δφ̂_4 = -0.018                       │
│                  current GWTC-3 σ_combined = 0.18                │
│                  ratio TGP/σ ≈ 0.10 (10× below sensitivity)      │
│                                                                  │
│  N events needed for 5σ in generic ppE-marginalized analysis:    │
│            ~230,000 BBH (impossibly far for ET-D + CE projected) │
│                                                                  │
│  → TGP M911-P1 jest "NIE wykrywalny" w generic LIGO ToGR analysis │
│    NAVET w erze 3G                                              │
│                                                                  │
│  ALE: Path A (op-LIGO-3G-deviation) bounds są valid dla:         │
│       - TGP-specific single-coefficient prior analysis           │
│       - Multi-coefficient ratio test (M911-P2)                   │
│                                                                  │
│  → TGP detection wymaga **DEDICATED TGP analysis pipeline**,     │
│    nie generic ToGR ppE marginalization.                         │
└──────────────────────────────────────────────────────────────────┘
```

## §1 — Co to znaczy dla T01

### 1.1 Status M911-P1 zostaje LIVE

GWTC-3 reanalysis NIE rejects ani NIE confirms TGP. Predykcja
pozostaje **LIVE**, ale z nowym important caveat:

> **Detection of M911-P1 wymaga dedicated TGP-specific analysis
> pipeline** (not generic LIGO ToGR ppE marginalization). Ten caveat
> NIE zmienia structural prediction (β_ppE^TGP^(b=-1) = -5/64), ale
> *zmienia* expectation o tym, *jak* wykryć:
> - **NIE** spodziewać się sygnału w future LIGO ToGR papers
>   (multi-coef marginalized).
> - **MUSI** przeprowadzić TGP-specific Bayes inference (single-coef
>   prior z β = -5/64).
> - **MOŻE** używać multi-coefficient ratio test (M911-P2) jako
>   alternative discriminator — który **NIE** wymaga absolute
>   detection każdego coefficient, ale tylko *consistency* ratio.

### 1.2 Co to oznacza dla paper draft (Path D)

[[../../papers/M911_LIGO3G_paper/paper_draft.md]] sekcja 4 i 5
wymaga **uzupełnienia** o:

- §4.4: Note że LIGO ToGR generic bounds nie pokażą TGP signal
  z aktualnych public posteriors, ale TGP-specific Bayes
  inference (using β_TGP = -5/64 prior) jest **a different test**
  z znacznie wyższym detection power.
- §5.2: Window of testability **CAVEAT** — generic ToGR analysis
  ET-D + CE NIE detect TGP at 5σ z ~10⁵ events; ALE TGP-specific
  prior analysis (single-coefficient) detect z `~10` events ET-D.

Te są **ważne uzupełnienia** dla paper, znacznie zwiększające
honest assessment.

## §2 — Podstawa reconciliation z Path A

### 2.1 Source kontextu

Path A ([[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] §2)
podała:

| Detector | β_5σ_calibrated | TGP/β_5σ | Verdict |
|----------|------------------|-----------|---------|
| LIGO-O5 single | ~3·10⁻² | ~2.6 | YES first decisive 2027+ |
| ET-D single | ~10⁻² | ~7.8 | YES (>20σ) |
| CE single | ~3·10⁻³ | ~26 | YES decisive |

vs Phase 2 GWTC-3 reanalysis:

| Source | δφ̂_4_5σ_bound | β_5σ_absolute | TGP β/bound |
|--------|------------------|----------------|----------------|
| GWTC-3 (~90 BBH) | ~0.91 (5σ) | ~3.9 | 0.020 |

**Discrepancy: factor ~50** w bounds (3·10⁻² Path A vs 3.9 Phase 2).

### 2.2 Resolution

Te są **two different analyses** z fundamentally different priors:

**Path A approach:**
- Prior: single-coefficient ppE deviation (assume only β_2PN modified)
- Marginalize nad: M_chirp, η, χ_eff (intrinsic params, degeneracy_factor=5)
- Result: σ_β z single-coefficient Fisher

**LIGO ToGR / Phase 2 approach:**
- Prior: generic deviation (modify all 8 PN coefficients simultaneously)
- Marginalize nad: M_chirp, η, χ_eff PLUS all 8 PN coefficients
- Result: σ_δφ̂_n z multi-coefficient Fisher (much weaker per
  coefficient)

**Difference factor ~10-50×** je **expected** — multi-coefficient
marginalized bounds zawsze są weaker niż single-coefficient.

### 2.3 Implikacja: dwa workflows dla M911-P1 detection

**Workflow A (TGP-specific, optimistic):**
- Single-coefficient β_TGP prior (delta function at -5/64).
- Bayes evidence comparison: M(β = β_TGP) vs M(GR, β = 0).
- Estimated detection: ET-D single event (calibrated bounds Path A).
- Required: dedicated TGP analysis pipeline (custom prior implementation).

**Workflow B (generic ppE, conservative):**
- Multi-coefficient marginalized fit (LIGO ToGR style).
- Independent posteriors for δφ̂_4, δφ̂_3, δφ̂_5, etc.
- Estimated detection: ~230,000 BBH events for 5σ — impossibly far.
- Standard analysis pipeline.

**Workflow C (multi-coefficient ratio, M911-P2):**
- Test ratios {β_3PN/β_2PN = -23/10, β_4PN/β_3PN = -38/23, ...}
- Use existing multi-coefficient Fisher results from ToGR.
- Required: extract coefficient *correlation* matrix z published
  posteriors + check consistency with TGP-predicted ratios.
- Estimated detection: stronger than Workflow B, weaker than A.
- Non-standard analysis — wymaga custom inference.

## §3 — Recommended workflow dla paper draft + future research

### 3.1 Bezpośrednie next steps

**(R1) — Workflow A (TGP-specific Bayes) z public GWTC-3 events:**

- Estimated work: ~3-5 sesji Python (bilby + pycbc + custom prior).
- Zwraca: log10(BF_TGP/GR) z z aktualnych ~90 BBH events.
- **Significance estimate (Workflow A):** dla 90 events, σ_β
  single-coefficient ~ 0.06 (calibrated), TGP signal 0.078 →
  Z = 0.078 / 0.06 ≈ 1.3 → ~1-2σ tentative signal.
- **WORKABLE.** Może dać first pre-publication signal jeśli TGP
  jest correct.

**Recommendation:** **(R1) jest najbardziej obiecujący path** dla
pre-publication confirmation. Single-session quick analysis (1
sesja Python) z Bayes prior framework dla GWTC-3 events daje
expected ~1-2σ tentative signal.

### 3.2 Long-term (Path D paper completion)

- **(R2)** Multi-coefficient ratio test (Workflow C):
  ekstraktować correlation matrix z published GWTC-3 posteriors,
  check consistency z TGP ratios. ~2-3 sesji.
- **(R3)** Production Fisher rerun z bilby + official noise curves
  (poprzedni TODO #3 z C-B-A-D session). Calibrowane absolute β bounds,
  ~2-3 sesji.
- **(R4)** ET-D + CE simulated dataset analysis (LIGO Algorithm
  Library MockData challenge): test workflow A na simulated future
  data. ~5-10 sesji.

### 3.3 Co paper draft musi reflectować

[[../../papers/M911_LIGO3G_paper/paper_draft.md]] needs **revision**:

- §1 Abstract: dodać "Detection of TGP M911-P1 in 3G era requires
  TGP-specific Bayes prior analysis (not generic ToGR multi-coef
  marginalization). Single-coefficient prior detection thresholds
  are achievable in ET-D single events."

- §4 Detection forecasts: dodać Table comparing single-coefficient
  (optimistic) vs multi-coefficient marginalized (conservative)
  bounds. Clarify which workflow is required for each.

- §5 Discussion: addr the GWTC-3 reanalysis result — TGP consistent,
  not detected, not excluded. Note this contradiction with naive
  Fisher forecasts and explain via marginalization.

- §6 Conclusion: M911-P1 prediction is testable, but **only with
  dedicated TGP analysis pipeline** in 3G era.

## §4 — Tier 5 cycle status

**Original goal:** "GWTC-3 reanalysis może dać sygnał *teraz* z
public data" (FINDINGS.md §2.5, SESSION_REPORT_2026-05-07_C-B-A-D.md
§10).

**Achieved:**
✓ Compiled published GWTC-3 + GWTC-2 ToGR bounds.
✓ Converted TGP β to LIGO fractional units.
✓ Bayes factor analysis (Laplace approximation): BF ≈ 0.97
  INCONCLUSIVE.
✓ Detection power assessment: ~230,000 events needed in generic
  ppE marginalized analysis.
✓ Reconciled Path A optimistic vs Phase 2 conservative bounds:
  expected difference from single- vs multi-coefficient priors.
✓ Identified **3 workflows** dla TGP detection w 3G era.
✓ Recommended **R1 (TGP-specific Bayes z GWTC-3)** as next session
  step.

**NOT achieved (out of scope for 1-session cycle):**
✗ Full Workflow A reanalysis z public GWTC-3 H5 strain data (wymaga
  bilby + pycbc + ~3-5 sesji).
✗ Multi-coefficient correlation matrix extraction (Workflow C).

## §5 — Updates downstream

### 5.1 Updates required dla T01 audit folder

| Plik | Update |
|------|--------|
| [[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]] | Dodać §1.X note: detection thresholds wymagają single-coefficient TGP-specific Bayes; generic ppE marginalized factor ~50× weaker |
| [[../../audyt/T01_LIGO3G_falsifier/FINDINGS.md]] | §2.5 Tier 5 status: EXECUTED 2026-05-07; verdict INCONCLUSIVE z GWTC-3 (consistent ale nie detected); recommended R1 (TGP-specific Bayes ~1 sesja) jako natural next |
| [[../../audyt/T01_LIGO3G_falsifier/NEEDS.md]] | Q3 closed: GWTC-3 reanalysis EXECUTED, verdict consistent (no falsification, no confirmation w generic ToGR; tentative signal possible w TGP-specific workflow) |
| [[../../audyt/T01_LIGO3G_falsifier/SESSION_REPORT_2026-05-07_C-B-A-D.md]] | §10 Tier 5 status: EXECUTED |

### 5.2 Updates required dla rejestr

[[../../PREDICTIONS_REGISTRY.md]] M911-P1 wpis nie wymaga *change*
status (LIVE pozostaje), ale **caveat note** może być dodana:
"Detection w 3G era wymaga TGP-specific Bayes prior, NIE generic
ToGR multi-coefficient marginalized analysis."

### 5.3 Update dla Path A

[[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] §5
"side-cykl rekomendowany" → **EXECUTED** (verdict: consistent,
detection requires dedicated TGP pipeline workflow A).

### 5.4 Update dla paper

[[../../papers/M911_LIGO3G_paper/paper_draft.md]] requires §4-§5
revision per §3.3 above.

## §6 — Cross-references

- [[README.md]] — overview
- [[Phase0_balance.md]] — M03 gate
- [[Phase1_GWTC3_bounds.md]] — bounds compilation + conversion convention
- [[Phase2_Bayes_factor.md]] — Bayes analysis + reconciliation
- [[scripts/bayes_factor_tgp_vs_gr.py]] — Python script
- [[scripts/bayes_factor_tgp_vs_gr.txt]] — run output
- [[../../audyt/T01_LIGO3G_falsifier/FINDINGS.md]] — pre-Tier-5 synthesis
- [[../../audyt/T01_LIGO3G_falsifier/SESSION_REPORT_2026-05-07_C-B-A-D.md]] — referencja Tier 5 source
- [[../op-ppE-mapping/Phase1_results.md]] — β_TGP source
- [[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] — Path A
- [[../../PREDICTIONS_REGISTRY.md]] M911-P1 — target rejestru
