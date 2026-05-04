---
title: "D01 — liczbowy dryft (α_s, m_H, Φ₀, Σm_ν, g₀^e, …)"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/D01_drifting_numbers
tags:
  - audit
  - numerical-drift
  - parameters
  - anchor-lock
  - propagation
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../L04_ODE_dualism_alpha]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
  - "[[../../meta/CALIBRATION_PROTOCOL.md]]"
tgp_status:
  folder_status: audit
  level: L2
  kind: audit
  core_compatibility: partial
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["B3-v2-alpha-s", "C10-v2-Phi0", "B4-sigma-m-nu"]
  depends_on:
    - "[[../L04_ODE_dualism_alpha]]"
  impacts: []
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §B.3, §B.4, §C.10, §C.9"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# D01 — liczbowy dryft parametrów TGP

## Klasa: LICZBOWE NIESPÓJNOŚCI (anchor-lock pending) • Priorytet: P3

## Diagnoza

**TGP twierdzi `40 predykcji z 3 inputów`** (g₀^e, Ω_Λ, N=3). Faktyczne
liczby kluczowych parametrów dryfują w rdzeniu o 0.4–6%. Każdy spread
osobno jest mały, razem tworzą „shopping list" dla autora wybierającego
najlepszą wartość pod konkretną predykcję.

## Pełna lista driftów

### α_s(M_Z) — 4 wartości w żywych plikach

| Wartość | Lokalizacja | Sigma vs PDG 0.1180 ± 0.0009 |
|---------|-------------|-------------------------------|
| **0.1190** | README pre-lock + AUDYT (oryginalny) | 1.1σ |
| **0.1184** | sek00 / ROADMAP_v3:883 / slownik (B3-locked 2026-05-01) | 0.4σ |
| **0.1174** | LP-5 13/13 PASS (`lp5_quark_m0_from_string_tension.py`) | 0.7σ |
| **0.1171** | LP-7 6/6 PASS (`lp7_defect_phase_emergence.py`) | 1.0σ |

Audit § B.3 — **CLOSED structurally** (anchor lock 0.1184); **B3-v2 pending**:
*„pełna propagacja przez tooling/scripts + sek00/companion (15+ plików dotkniętych)"*.
README.md już zaktualizowany 2026-05-01, ale skrypty LP-5, LP-7 wciąż
podają inne wartości.

### m_H — 4 wartości

| Wartość | Lokalizacja | Sigma vs PDG 2024 125.20 ± 0.11 |
|---------|-------------|----------------------------------|
| **125.31** | README, TGP claim (`v × 57/112`) | 1.0σ |
| **125.25** | older PDG anchor (was used in older sek09) | 0.3σ (ze starym PDG) |
| **125.20** | PDG 2024 latest | 0σ (target) |
| **125.1** | dodatekU, sek09 (jeszcze starsze) | 0.9σ |

Audit § C.9 — CLOSED structurally (PDG re-anchor done w README:97, 280).
sek09 i dodatekU pozostają nie-zsynchronizowane.

### Φ₀ — 3 wartości w ROADMAP_v3

| Wartość | Lokalizacja | Pochodzenie |
|---------|-------------|-------------|
| **24.66** | ROADMAP_v3:65 | kalibracja z Λ_obs |
| **24.783** | ROADMAP_v3:122, 896 | F11 mass formula consistency (m_H = v × 57/112) |
| **25.0** | inne miejsca | algebraiczne (`N_f² = (2N_c-1)²`) lub `a_Γ·Φ₀ = 1` |

Spread 0.4%. α_s drugiego rzędu zależny od Φ₀ → propaguje 0.4%
niejednoznaczność w wszystkich predykcjach kosmologicznych.

Audit § C.10 — **CLOSED annotation-only** (anchor proposed Φ₀ = 24.783);
**C10-v2 pending**: pełny propagacja przez 15+ plików nie wykonana.

### Σm_ν — 3 wartości

| Wartość (meV) | Lokalizacja | Status |
|---------------|-------------|--------|
| **59.01** | PREDICTIONS_REGISTRY Z1 anchor (B4-anchored 2026-05-01) | LOCKED post-B4 |
| **59.6** | README + zeroth-order target | pre-bisection |
| **62.9** | historical | obsolete |

Audit § B.4 — CLOSED annotation-only. README wciąż używa 59.6.

### g₀^e — 4 wartości (LP-6 N-1)

| Wartość | Formulacja ODE |
|---------|----------------|
| **0.834** | Pełne f(g) |
| **0.869** | Substratowa (K=g²) — LP-5/LP-6 preferowana |
| **0.870** | (zaokrąglona) — README, claim „TGP coupling g₀^e = 0.86941" |
| **0.899** | Uproszczone |

Spread 7% w g₀^e — *kluczowy wejściowy parametr teorii*. Wartość zależy
od wyboru formulacji ODE ([[../L04_ODE_dualism_alpha]]).

### α (kinetic) — dualizm

- **α=1** (substratowa K=g²): faktycznie używana w skryptach LP-4, LP-5, LP-6
- **α=2** (kanoniczna K=g⁴): deklarowana w README jako „algebraic theorem"

Patrz [[../L04_ODE_dualism_alpha]].

### k (mass exponent) — pivot

- **k=4** (LP-4 oryginalne, 9/9 PASS): argument konwergencyjny
- **p=5−α** (R3 2026-05-01): zależne od formulacji

Patrz [[../L05_mass_exponent_drift]].

## Co to znaczy dla predyktywności

README claim: *„From three inputs … the theory derives 40 quantitative
predictions"* + LS-8 audit: *„Predictivity ratio = 11/2 = 5.5"*.

Ta ratio jest podana **dla wybranej trajektorii wartości**:

- α_s = 0.1184 (lock z 4 wartości)
- Φ₀ = 24.783 (lock z 3 wartości)
- m_H = 125.31 / 125.20 (PDG 2024)
- g₀^e = 0.869 (substratowa)
- k = 4 lub p = 3 (zależnie od α)
- Σm_ν = 59.01 lub 59.6

Każda inna trajektoria daje inne predykcje. Bez globalnego anchor lock
i propagacji, *predictivity ratio jest księgowo zawyżony*.

## Wpływ na falsifiability

Każdy z parametrów ma niezerowy spread, więc każda predykcja zależna
ma wbudowaną niepewność wyboru anchora. Dla 5 falsyfikatorów TGP:

- **Inverted neutrino mass ordering** — depends on Σm_ν anchor
- **Proton decay observed** — niezależne, OK
- **w_DE < -1 (>5σ)** — depends on Φ₀ → DE prediction
- **c_GW ≠ c (>5σ)** — depends on metric form (S01, S02)
- **No breathing mode in 1000+ events** — depends on σ_ab decision (S05)

3 z 5 falsyfikatorów wisi na anchor choices.

## Status w audycie

| Parametr | Audit ref | Status |
|----------|-----------|--------|
| α_s | §B.3 | CLOSED structurally; B3-v2 pending |
| m_H | §C.9 | CLOSED structurally; sek09/dodatekU sync pending |
| Φ₀ | §C.10 | CLOSED annotation-only; C10-v2 pending |
| Σm_ν | §B.4 | CLOSED annotation-only |
| g₀^e (formulacje) | LP-6 N-1 | rekomendacja unifikacji nie wcielona |

## Rekomendacja

Otworzyć dedykowany cykl `op-anchor-lock-propagation/`:

### Phase 1 — globalna decyzja per parametr

Po decyzji L04 (formulacja ODE):

- **α_s = 0.1184** (anchor już wybrany 2026-05-01)
- **m_H** TGP-side = 125.31 (formal `v × 57/112`); reference = 125.20 (PDG 2024)
- **Φ₀ = 24.783** (centralna z F11)
- **Σm_ν = 59.01 meV** (Z1 anchor)
- **g₀^e = 0.86941** (z wybranej formulacji)
- **k vs p = 4 vs 3** (po L04 decyzji)

### Phase 2 — globalna propagacja przez tooling

Globalne grep + replace przez:

- `core/sek*.tex` (~30 plików)
- `companion.tex`, `letter.tex`
- `arxiv_submission/`, `paper_*/`
- `scripts/` (~40 skryptów)
- `nbody/` (~20 plików)
- `PREDICTIONS_REGISTRY.md`, `INDEX.md`, `README.md`

### Phase 3 — sanity check

Globalna re-evaluacja predyktywności:

- Re-run wszystkich verification scripts z lockowanymi wartościami
- Re-derive `predictivity ratio` dla *jednej* trajektorii
- Confirm że falsyfikatory są dobrze ustawione

**Estymata:** 2 tygodnie globalnej propagacji + 1 tydzień sanity check.

## Pliki dotknięte (estymata)

15+ plików zidentyfikowanych przez audit jako pending B3-v2 i C10-v2.
Pełna lista po wykonaniu Phase 1.

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §III
- [[../PRIORITY_MATRIX.md]] klaster C
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §B.3, §B.4, §C.9, §C.10
- [[../../meta/CALIBRATION_PROTOCOL.md]]
- [[../L04_ODE_dualism_alpha]] (g₀^e zależy od formulacji)
- [[../L05_mass_exponent_drift]] (k vs p)
