---
title: "FINDINGS — op-eta-wolfenstein"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-eta-wolfenstein
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-eta-wolfenstein

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `η.1.Phase1` | `CLOSED` | `PASS` | — | — |
| `Phase2_results.md` | `η.1.Phase2` | `CLOSED` | `PASS` | — | — |
| `Phase3_results.md` | `η.1.Phase3` | `CLOSED` | `PASS` | — | `END` |

## Numerical PASS counts (cited from .txt outputs)

- `phase1_wolfenstein_audit.txt` — 5/5 PASS — context: …====   T1.1: PASS   T1.2: PASS   T1.3: PASS   T1.4: PASS   T1.5: PASS    Total: 5/5 PASS   → η.1.Phase2 proceeds z best-candidate triple. =====================
- `phase2_wolfenstein_derivation.txt` — 7/7 PASS — context: …PASS   T2.3: PASS   T2.4: PASS   T2.5: PASS   T2.6: PASS   T2.7: PASS    Total: 7/7 PASS   → η.1.Phase3 proceeds; (A, ρ̄, η̄) triple LOCKED. ==================
- `phase3_eta_predictions.txt` — 6/6 PASS — context: …PASS   T3.2: PASS   T3.3: PASS   T3.4: PASS   T3.5: PASS   T3.6: PASS    Total: 6/6 PASS   → η.1 program END (6/6 PASS).   → Cumulative 18/18 PASS (Phase1 5 + 
- `phase3_eta_predictions.txt` — 6/6 PASS — context: …T3.4: PASS   T3.5: PASS   T3.6: PASS    Total: 6/6 PASS   → η.1 program END (6/6 PASS).   → Cumulative 18/18 PASS (Phase1 5 + Phase2 7 + Phase3 6).   → Ledger 
- `phase3_eta_predictions.txt` — 18/18 PASS — context: …T3.6: PASS    Total: 6/6 PASS   → η.1 program END (6/6 PASS).   → Cumulative 18/18 PASS (Phase1 5 + Phase2 7 + Phase3 6).   → Ledger 427 → 445 (+18).   → Cla

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **T1.1** | A = 79/100 best, top-5 ranked, drift 0.0% | **PASS** |
> | **T1.2** | ρ̄ = 11/78 best, top-5 ranked, drift 0.018% | **PASS** |
> | **T1.3** | η̄ = 5/14 #1, drift 0.040% < 0.1% | **PASS** |
> | **T1.4** | Unitarity triangle α+β+γ = 180° exact closure | **PASS** |
> | **T1.5** | Jarlskog J_TGP = 2.93·10⁻⁵, drift 4.57% < 5% | **PASS** |
> 
> **5/5 PASS** → η.1.Phase2 proceeds z best-candidate triple (A, ρ̄, η̄).
> 
> ---

### `Phase2_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **T2.1** | A = 64/81 = 8²/3⁴ sympy candidate (drift 0.0156%) | **PASS** |
> | **T2.2** | ρ̄ = 11/78 sympy candidate (drift 0.0182%) | **PASS** |
> | **T2.3** | η̄ = 5/14 sympy LOCKED (drift 0.0400%) | **PASS** |
> | **T2.4** | Cross-sector denom-family analysis (81, 78, 14 coprime) | **PASS** |
> | **T2.5** | 5/5 alternative triples FALSIFIED at 0.5% threshold | **PASS** |
> | **T2.6** | V_ub_TGP refined cascade (drift 8.93% < 9%) | **PASS** |
> | **T2.7** | Classification PARTIALLY DERIVED (refined) | **PASS** |
> 
> **7/7 PASS** → η.1.Phase3 pr...(truncated)

### `Phase3_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **T3.1 (H1)** | Belle II 2027+ |V_ub| = 0.00348 within window [0.00340, 0.00400] | **PASS** |
> | **T3.2 (H2)** | LHCb Run 4 J = 2.93·10⁻⁵ within window [2.85, 3.30]·10⁻⁵ | **PASS** |
> | **T3.3 (H3)** | sin(2β)_TGP = 0.7090, drift 1.43% within window [0.65, 0.75] | **PASS** |
> | **T3.4 (H4)** | Cross-sector denom-prime sharing (3 ↔ K_lepton, 7 ↔ K_up) | **PASS** |
> | **T3.5 (H5)** | |V_td/V_ts|_TGP = 0.2098, drift 2.33% within [0.195, 0.220] | **PASS** |
> | **T3.6 (H6)** | 4-channel η.1 convergence 4/4 LIVE | **PASS** |
> 
> **6/6 PASS** → ...(truncated)

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa