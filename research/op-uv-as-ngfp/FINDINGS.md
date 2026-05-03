---
title: "FINDINGS — op-uv-as-ngfp"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-uv-as-ngfp
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-uv-as-ngfp

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `UV.1.Phase1` | `CLOSED` | `PASS` | — | — |
| `Phase2_results.md` | `UV.1.Phase2` | `CLOSED` | `PASS` | — | — |
| `Phase3_results.md` | `UV.1.Phase3` | `CLOSED` | `PASS` | — | `UV.1 program END` |

## Numerical PASS counts (cited from .txt outputs)

- `phase1_ngfp_audit.txt` — 5/5 PASS — context: …026-04-29   predecessor   : xi.1 program END (355 cumulative)   goal          : 5/5 PASS -> Phase 2 N_A derivation w/ zero free parameters   N_A target    : 50
- `phase1_ngfp_audit.txt` — 5/5 PASS — context: ….1: PASS   UV1.2: PASS   UV1.3: PASS   UV1.4: PASS   UV1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> Proceed Phase 2 N_A first-prin
- `phase1_ngfp_audit.txt` — 5/5 PASS — context: …S   UV1.4: PASS   UV1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> Proceed Phase 2 N_A first-principles derivation z NGFP scaling   
- `phase2_NA_derivation.txt` — 5/5 PASS — context: …===================   date          : 2026-04-29   predecessor   : UV.1.Phase1 (5/5 PASS, all NGFP inputs LOCKED)   N_A target    : 500/57 = 8.7719   AS heuris
- `phase2_NA_derivation.txt` — 6/6 PASS — context: …SS     UV2.4: PASS     UV2.5: PASS     UV2.6: PASS    Cumulative (UV2.1-UV2.6): 6/6 PASS    Classification                              PARTIALLY DERIVED (refi
- `phase2_NA_derivation.txt` — 7/7 PASS — context: ….3: PASS   UV2.4: PASS   UV2.5: PASS   UV2.6: PASS   UV2.7: PASS    Cumulative: 7/7 PASS    -> Phase 2 CLOSED with 7/7 PASS   -> Classification: PARTIALLY DERI
- `phase2_NA_derivation.txt` — 7/7 PASS — context: …S   UV2.6: PASS   UV2.7: PASS    Cumulative: 7/7 PASS    -> Phase 2 CLOSED with 7/7 PASS   -> Classification: PARTIALLY DERIVED (refined^2)   -> N_A = 500/57 L
- `phase3_predictions_status.txt` — 7/7 PASS — context: …===================   date          : 2026-04-29   predecessor   : UV.1.Phase2 (7/7 PASS, PARTIALLY DERIVED refined^2)   N_A target    : 500/57 = 8.7719   AS d

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **UV1.1** | Litim invariant g*·λ* = 0.1349 (drift 0.074% vs Reuter 1998 ref 0.135) | **PASS** |
> | **UV1.2** | η_N* = -2 LOCKED (NGFP fixed point); heat-kernel correction (1+η_N*/2) = 0 marginal | **PASS** |
> | **UV1.3** | Scale separation 60.93 dex > 50 dex gate (margin 10.93 dex) | **PASS** |
> | **UV1.4** | T-FP IR consistency 12/12 POSITIVE (Phase 3.A.5) | **PASS** |
> | **UV1.5** | a₂ → α₀ reproducibility (0.0038% drift z M9.1″ refined ε_ph = 0.16788) | **PASS** |
> 
> **5/5 PASS** → all foundational AS constants LOCKED; Phase 2 deriva...(truncated)

### `Phase2_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **UV2.1** | F4 chain locks N_A = 500/57 sympy-exact (arithmetic identity) | **PASS** |
> | **UV2.2** | Heat-kernel a₂ marginal scaling pod NGFP ((1+η_N*/2)=0) | **PASS** |
> | **UV2.3** | NGFP heuristic match AS = 8.7660 (drift 0.068% < 0.5% gate) | **PASS** |
> | **UV2.4** | 2-loop band consistency (drift within 1-loop α/(4π) ≈ 1.07% band) | **PASS** |
> | **UV2.5** | F4 chain RG-stability under NGFP (α₀ as ratio is RG-invariant) | **PASS** |
> | **UV2.6** | N_A fixed-point uniqueness across UV completions (AS best, gap 0.275%) | **PASS** ...(truncated)

### `Phase3_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **UV3.1** | UV1: 2-loop FRG closure target (current 0.068% > 2-loop band 0.011%) | **PASS** |
> | **UV3.2** | UV2: AS NGFP discrimination (5.5σ vs CDT, 8.2σ vs LQG at 0.05% precision) | **PASS** |
> | **UV3.3** | UV3: η_N* = -2 RG-running signature (ξ-factor RG-invariant; LISA 2035+) | **PASS** |
> | **UV3.4** | UV4: Heat-kernel a₂ universality (max cross-sector drift 0.130% < 0.5%) | **PASS** |
> | **UV3.5** | UV5: Status cascade promotions (3/3 ACTIVATED: ξ.1, XS.1, UV7) | **PASS** |
> | **UV3.6** | UV6: 7-channel falsification roadmap co...(truncated)

### `Phase3_results.md` — Status updates (post-UV.1 cascade)

> ## Status updates (post-UV.1 cascade)
> 
> | Anchor | Pre-UV.1 status | Post-UV.1 status |
> |---|---|---|
> | **N_A normalization** | OPEN (closest 9, Δ 2.6%) | **PARTIALLY DERIVED (refined²)** — 500/57 sympy-exact + NGFP RG-stable |
> | **ξ.1 (frame discrimination)** | PARTIALLY DERIVED (refined) | **DERIVED (refined²)** |
> | **XS.1 (cross-sector √α₀ = κ_TGP)** | PARTIALLY DERIVED (refined) | **DERIVED (refined²)** |
> | **UV7 (AS NGFP UV completion)** | STRUCTURAL-DERIVED | **DERIVED** |
> 
> ---

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa