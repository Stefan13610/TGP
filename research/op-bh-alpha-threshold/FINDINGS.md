---
title: "FINDINGS — op-bh-alpha-threshold"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-bh-alpha-threshold
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-bh-alpha-threshold

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `BH.1.Phase1` | `CLOSED` | `H0_REJECTED` | — | — |
| `Phase2_results.md` | `BH.1.Phase2` | `CLOSED` | `T_ALPHA_UPGRADED` | — | — |
| `Phase3_results.md` | `BH.1.Phase3` | `CLOSED` | `PASS` | — | — |

## Numerical PASS counts (cited from .txt outputs)

- `phase1_multi_source_audit.txt` — 5/5 PASS — context: …a(psi) threshold) is UNIQUE viable minimal extension.            Phase 1 closes 5/5 PASS.  Phase 0+ multi-source ISSUE FORMALIZED.            Ready for Phase 2
- `phase3_multi_source_falsification.txt` — 7/7 PASS — context: …esign     : PASS   T3.7  PREDICTIONS_REGISTRY BH4-BH9        : PASS    VERDICT: 7/7 PASS.  Multi-source falsification map COMPLETE.            6 new prediction

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — TL;DR

> ## TL;DR
> 
> > Pełna analiza wymiarowa Candidate D action `S_int ∝ α · T^μν J_μ J_ν √-g d⁴x`
> > pokazuje że **dim[α] nie jest dimensionless w SI** pod żadną z trzech
> > naturalnych konwencji (mass-flux J, EM-like J, geometric-natural).
> > α dziedzicznie niesie wymiar `[m²]` (lub ekwiwalent), co bezpośrednio
> > wymusza scaling α_SI = α_geom · R_S² → **M² scaling** dla Schwarzschilda.
> >
> > Numeryczna demonstracja dla 4 źródeł potwierdza idealny slope = 2.000
> > (log10 α_SI vs log10 M_BH/M_⊙) i 19.33 orders of magnitude w α_SI dla
> > 9.67 orders w M_BH.
> >
> > Trzy proste resolution paths (L_TGP intrinsic len...(truncated)

### `Phase2_results.md` — TL;DR

> ## TL;DR
> 
> > Phase 2 upgraduje T-α closure (2026-04-26) o:
> > 1. **EFT/symmetry argument** dla n=2: Taylor expansion + Z₂ reflection symmetry around ψ_th + threshold-vanishing condition → tylko parzyste potęgi (ψ-ψ_th)^{2k} dopuszczalne; minimum k=1 → n=2 jako leading order.
> > 2. **WEP-MICROSCOPE-2 lower bound**: η_TGP = α₀(ψ_Earth-1)^n. Dla projected MICROSCOPE-2 sensitivity 10⁻¹⁷, n=1 daje η=2.8·10⁻⁹ (fails by 8 orders), n=2 daje η=1.97·10⁻¹⁸ (passes z margin 5.1×). **n ≥ 2 ściśle wymagane**.
> > 3. **Non-overkill upper bound**: n=3 daje η=10⁻²⁷ (margin 10⁹×, sygnał TGP undetectable forever); n=...(truncated)

### `Phase3_results.md` — Verdict

> ## Verdict
> 
> | Sub-test | Description | Result |
> |---|---|---|
> | **T3.1** | ngEHT photon ring multi-source map (10 SMBH) | **PASS** |
> | **T3.2** | LIGO O5 / LISA ringdown frequency shift | **PASS** |
> | **T3.3** | NICER pulsar M-R relation | **PASS** |
> | **T3.4** | MICROSCOPE-2 WEP test | **PASS** |
> | **T3.5** | Solar PPN Cassini-class precision | **PASS** |
> | **T3.6** | Cross-sector consistency falsification design | **PASS** |
> | **T3.7** | PREDICTIONS_REGISTRY entries BH4–BH9 generation | **PASS** |
> 
> **Cumulative BH.1:** 5 (Phase 1) + 7 (Phase 2) + 7 (Phase 3) = **19 sub-tests**.
> 
> ---

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa