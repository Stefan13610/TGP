---
title: "FINDINGS — op-theta-quark-koide"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-theta-quark-koide
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-theta-quark-koide

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `θ.1.Phase1` | `CLOSED` | `PASS` | — | — |
| `Phase2_results.md` | `θ.1.Phase2` | `CLOSED` | `PASS` | — | — |
| `Phase3_results.md` | `θ.1.Phase3` | `CLOSED` | `PASS` | — | `END` |

## Numerical PASS counts (cited from .txt outputs)

- `phase1_kquark_audit.txt` — 5/5 PASS — context: …orana B^2=1)   lambda_C ref  : 0.22550  (zeta.1 cross-sector)   goal          : 5/5 PASS -> proceed Phase 2 K_quark first-principles  =========================
- `phase1_kquark_audit.txt` — 5/5 PASS — context: …T1.1: PASS   T1.2: PASS   T1.3: PASS   T1.4: PASS   T1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> K_up   = 0.874559  LOCKED (PDG M
- `phase1_kquark_audit.txt` — 5/5 PASS — context: …ASS   T1.4: PASS   T1.5: PASS    Cumulative: 5/5 PASS    -> Phase 1 CLOSED with 5/5 PASS   -> K_up   = 0.874559  LOCKED (PDG MS-bar @ M_Z)   -> K_down = 0.7399
- `phase2_kquark_derivation.txt` — 7/7 PASS — context: …7  (Dirac B^2=2)   K_nu     ref  : 0.500000  (Majorana B^2=1)   goal          : 7/7 PASS -> Phase 3 predictions Q1-Q6  ========================================
- `phase2_kquark_derivation.txt` — 7/7 PASS — context: …T2.3: PASS   T2.4: PASS   T2.5: PASS   T2.6: PASS   T2.7: PASS    Cumulative: 7/7 PASS    -> Phase 2 CLOSED with 7/7 PASS   -> K_up = 7/8 sympy LOCKED via B^2_
- `phase2_kquark_derivation.txt` — 7/7 PASS — context: …ASS   T2.6: PASS   T2.7: PASS    Cumulative: 7/7 PASS    -> Phase 2 CLOSED with 7/7 PASS   -> K_up = 7/8 sympy LOCKED via B^2_up = 13/4 chirality-counting   ->
- `phase3_theta_predictions.txt` — 6/6 PASS — context: ….046% vs PDG)   K_down = 37/50: best candidate (drift 0.014%)   goal          : 6/6 PASS -> theta.1 program END (ledger 427)  =================================
- `phase3_theta_predictions.txt` — 6/6 PASS — context: …(Q3): PASS   T3.4 (Q4): PASS   T3.5 (Q5): PASS   T3.6 (Q6): PASS    Cumulative: 6/6 PASS    -> Phase 3 CLOSED with 6/6 PASS   -> 6 predictions Q1-Q6 LIVE (Bell

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **T1.1** | K_up = 0.8746 sympy z PDG MS-bar @ M_Z (drift 0.0047%) | **PASS** |
> | **T1.2** | K_down = 0.7398 sympy z PDG MS-bar @ M_Z (drift 0.0136%) | **PASS** |
> | **T1.3** | K_quark RG-invariance under common β-rescaling (max 10⁻²⁹%) | **PASS** |
> | **T1.4** | 4-sector K-taxonomy distinct (min pair drift 9.9% > 1%) | **PASS** |
> | **T1.5** | Wolfenstein λ-cascade vs PDG (λ_C drift 0.44% < 1%) | **PASS** |
> 
> **5/5 PASS** → K_quark numerical landscape LOCKED; Phase 2 first-principles proceeds.
> 
> ---

### `Phase2_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **T2.1** | K_up = 7/8 sympy candidate (drift 0.046%) | **PASS** |
> | **T2.2** | K_down candidates ranking (best 37/50, drift 0.014%) | **PASS** |
> | **T2.3** | B²_up = 13/4, B²_down = 2.4388 chirality-counting decomp | **PASS** |
> | **T2.4** | Cross-sector V_us = λ_C (drift 0.22%) + ratio √2 (drift 6.77%) | **PASS** |
> | **T2.5** | 5 alternative K_up formulas FALSIFIED (drifts 1.6-42.8%) | **PASS** |
> | **T2.6** | NGFP RG-stability common β-rescaling (max drift 10⁻²⁹%) | **PASS** |
> | **T2.7** | Classification PARTIALLY DERIVED (refined...(truncated)

### `Phase3_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **T3.1 (Q1)** | Belle II 2027+ |V_ub| = 0.00348 within window [0.00340, 0.00400] | **PASS** |
> | **T3.2 (Q2)** | LHCb Run 4 Jarlskog J = 2.93·10⁻⁵ within window [2.85, 3.30]·10⁻⁵ | **PASS** |
> | **T3.3 (Q3)** | EIC 2030+ proton mass-radius cross-check (universal γ_m) | **PASS** |
> | **T3.4 (Q4)** | Lepton-quark cross-sector ratio sin θ_C/sin θ₁₃ vs √2 (drift 7.5%) | **PASS** |
> | **T3.5 (Q5)** | K-taxonomy 4-sector completion (universal pattern verified) | **PASS** |
> | **T3.6 (Q6)** | 4-channel falsification convergence 4/4 LIVE | **P...(truncated)

### `Phase3_results.md` — Status cascade post-θ.1

> ## Status cascade post-θ.1
> 
> | Anchor | Pre-θ.1 | Post-θ.1 |
> |---|---|---|
> | K_up = 7/8 | OPEN | **PARTIALLY DERIVED (refined)** |
> | K_down ≈ 0.74 | OPEN | **STRUCTURAL refined** (best 37/50) |
> | B² chirality decomposition | OPEN | **DERIVED** |
> | Cross-sector λ_C anchor | DERIVED (ζ.1) | **DERIVED refined** |
> | 4-sector K-taxonomy | partial | **completed** |
> 
> ---

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa