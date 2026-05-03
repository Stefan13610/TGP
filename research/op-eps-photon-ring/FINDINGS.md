---
title: "FINDINGS — op-eps-photon-ring"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-eps-photon-ring
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-eps-photon-ring

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `ε.1.Phase1` | `CLOSED` | `PASS` | — | — |
| `Phase2_results.md` | `ε.1.Phase2` | `CLOSED` | `PASS` | — | — |
| `Phase3_results.md` | `ε.1.Phase3` | `CLOSED` | `PASS` | — | — |

## Numerical PASS counts (cited from .txt outputs)

- `phase1_eps_audit.txt` — 5/5 PASS — context: …026-04-29   predecessor   : UV.1 program END (373 cumulative)   goal          : 5/5 PASS -> Phase 2 structural decomposition  =================================
- `phase1_eps_audit.txt` — 5/5 PASS — context: …E1.1: PASS   E1.2: PASS   E1.3: PASS   E1.4: PASS   E1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> eps_ph numerical landscape fully
- `phase1_eps_audit.txt` — 5/5 PASS — context: …ASS   E1.4: PASS   E1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> eps_ph numerical landscape fully LOCKED   -> Proceed Phase 2 (str
- `phase2_eps_decomposition.txt` — 5/5 PASS — context: …==============   date          : 2026-04-29   predecessor   : epsilon.1.Phase1 (5/5 PASS, psi_ph=160/137)   hypothesis    : eps_ph = 23/137 (prime-137 denomina
- `phase2_eps_decomposition.txt` — 6/7 PASS — context: …hypothesis    : eps_ph = 23/137 (prime-137 denominator)   goal          : >= 6/7 PASS -> PARTIALLY DERIVED (refined)  =========================================
- `phase2_eps_decomposition.txt` — 6/6 PASS — context: …3: PASS     E2.4: PASS     E2.5: PASS     E2.6: PASS    Cumulative (E2.1-E2.6): 6/6 PASS    Classification                                  PARTIALLY DERIVED (
- `phase2_eps_decomposition.txt` — 7/7 PASS — context: …E2.3: PASS   E2.4: PASS   E2.5: PASS   E2.6: PASS   E2.7: PASS    Cumulative: 7/7 PASS    -> Phase 2 CLOSED with 7/7 PASS   -> Classification: PARTIALLY DERIVE
- `phase2_eps_decomposition.txt` — 7/7 PASS — context: …ASS   E2.6: PASS   E2.7: PASS    Cumulative: 7/7 PASS    -> Phase 2 CLOSED with 7/7 PASS   -> Classification: PARTIALLY DERIVED (refined)   -> eps_ph = 23/137 

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **E1.1** | ε_ph = 4197/25000 sympy-exact rational | **PASS** |
> | **E1.2** | ψ_ph = 4/(3+0.4250) = 160/137; ε_ph = ψ_ph − 1 = 23/137 (drift 0.0019%) | **PASS** |
> | **E1.3** | ε_ph² sympy exact; M9.2-D coarse = 441/15625 (clean) | **PASS** |
> | **E1.4** | F4 chain consistency α₀ = target_shift/ε_ph² = 4.04489 (drift 0.0038%) | **PASS** |
> | **E1.5** | Refinement gain M9.2-D → M9.1″ 36.3× (gate ≥30×) | **PASS** |
> 
> **5/5 PASS** → ε_ph numerical landscape fully LOCKED; Phase 2 proceeds.
> 
> ---

### `Phase2_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **E2.1** | ε_ph = ψ_ph − 1 = 23/137 structural identity (drift 0.0019%) | **PASS** |
> | **E2.2** | 5 alternative candidates falsification (drifts 10.7%–181.8%) | **PASS** |
> | **E2.3** | F4 chain implicit lock ε_ph² = 529/18769 (drift 0.0019%) | **PASS** |
> | **E2.4** | Heat-kernel a₂ frame quadratic consistency (ε_ph² in F4 chain) | **PASS** |
> | **E2.5** | M9.2-D vs M9.1″ refinement audit (drift 0.0100% << 0.527% a₂ band) | **PASS** |
> | **E2.6** | NGFP RG-stability of ε_ph² ratio (drift 0.0000%) | **PASS** |
> | **E2.7** | Classificat...(truncated)

### `Phase3_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **E3.1** | E1: ngEHT 2030+ r_ph = (160/137)·r_g, 0.1% precision (margin 5×) | **PASS** |
> | **E3.2** | E2: Cross-sector ε_ph² consistency, max drift 0.25% < 0.5% | **PASS** |
> | **E3.3** | E3: 5 identity candidates falsified (drifts 10.7%–181.8%) | **PASS** |
> | **E3.4** | E4: ε_ph² ratio RG-invariant (drift 0%), LISA 2035+ gate 0.5% | **PASS** |
> | **E3.5** | E5: F4 closure unique (single positive root, 0 alternatives <5%) | **PASS** |
> | **E3.6** | E6: 5-channel falsification convergence (margin 1 over ≥4 criterion) | **PASS** |
> 
> **6...(truncated)

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa