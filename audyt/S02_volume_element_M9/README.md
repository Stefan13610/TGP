---
title: "S02 — √(-g) niespójny z M9.1'' (efekt domina w M9.x)"
date: 2026-05-04
parent: "[[../README.md]]"
type: audit-issue
tgp_owner: audyt/S02_volume_element_M9
tags:
  - audit
  - volume-element
  - M9
  - critical
  - structural-contradiction
related:
  - "[[../SUMMARY_2026-05-04.md]]"
  - "[[../S01_metric_four_forms]]"
  - "[[../../research/op-newton-momentum]]"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tgp_status:
  folder_status: audit
  level: L2
  kind: audit
  core_compatibility: broken
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges: ["B6-M9.x-re-run"]
  depends_on:
    - "[[../S01_metric_four_forms]]"
  impacts:
    - "[[../S03_beta_PPN_convention]]"
  source_of_status:
    - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §A.2, §M.5"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# S02 — `√(-g)` niespójny z M9.1''

## Klasa: STRUKTURALNA SPRZECZNOŚĆ • Priorytet: P1

## Diagnoza

Dla M9.1'' poprawny element objętościowy:

```
det(g) = -c² · ψ²/(4-3ψ)²
√(-g) = c·ψ/(4-3ψ) = c · h(ψ)
```

**Wszystkie żywe wyprowadzenia w `core/` używają obsoletnego
`√(-g) = c·ψ`** skopiowanego z formy (I) `g_tt = -c²/ψ` (FALSIFIED
2026-04-25). Ta niespójność uderza w wszystkie predykcje grawitacyjne
TGP na poziomie M9.x.

## Dotknięte wyprowadzenia

| Element | Plik | Co używa |
|---------|------|----------|
| `κ = 3/(4Φ_0)` | `sek08a_akcja_zunifikowana.tex` lin. 451+ | `√(-g) = c·ψ` |
| Φ-EOM `prop:field-eq-from-action` | `sek08a` + `sek08_formalizm` | volume z (I) |
| Φ-EOM `eq:field-eq-reproduced` | `sek08a` lin. 80–86 | volume z (I) |
| `m_field` w M9.2 | `research/op-newton-momentum/M9_2_results.md` | volume z (I) |
| Kwadrupol Peters-Mathews w M9.3 | `research/op-newton-momentum/M9_3_results.md` | volume z (I) |
| c₂ koefitsent PPN | `research/op-newton-momentum/M9_1_pp_P3_results.md` | volume z (I) |

## Wpływ na status PASS

Audit § A.2 explicit:

> To **unieważnia strukturalnie 5/5 PASS w M9.2 i 5/5 PASS w M9.3**.

Aktualne 35/35 PASS w „M9.1''+M9.2+M9.3 phase-ledger" jest formalnie
**INTACT structurally**, ale **PRE-AUDIT NUMERICALLY** dla M9.2 i M9.3
(audit § M.4 tabela: „Numeryczne wartości M9.2 m_field, M9.3
Peters-Mathews → PRE-AUDIT").

## Wpływ na predykcje liczbowe

Skala efektu dla obserwowalnych:

- **M9.1 1PN PPN** (γ=β=1 w słabym polu): zmiana `√(-g)` o czynnik
  `1/(4-3ψ)` jest O(U) w słabym polu — γ_PPN nie zmienia się w wiodącym
  rzędzie, ale wkład do c₂ (i co za tym idzie β_PPN przez master formula)
  wymaga ponownego policzenia.
- **m_field**: wartość masy pola Φ wokół Φ_eq — zmiana `√(-g)` modyfikuje
  ekstrakcję `m²` z linearyzacji EOM. Skala efektu O(1).
- **Kwadrupol Peters-Mathews** (M9.3): formuła `dE/dt ∝ ⟨...⟩√(-g)` z
  poprawną metryką (IV) ma dodatkowy czynnik `1/(4-3ψ)` w volume. Dla
  pól słabych (binary inspiral) to korekta ~O(U) — może być istotna dla
  GW170817 sub-percent test `c_GW = c`.

## Status w `meta/AUDYT_TGP_2026-05-01.md`

§A.2 — **CRITICAL**, oznaczony „CLOSED structural" (§ M.6) z explicit
pending: **B6 dedicated repair cycle** (§ O.7).

> A2: B6 / dedicated repair cycle pending — pełny re-run M9.x z
> `√(-g) = c·ψ/(4-3ψ)`. Wymaga re-derivation κ, c_2, β_PPN, m_field,
> kwadrupol Peters-Mathews; structural M9 results 5/5 + 5/5 + 5/5 PASS
> intact, ale numeryczne wartości PRE-AUDIT.

Audit § M.4 explicit acknowledge: `√(-g) w M9.x derivations = OBSOLETE`.

## Rekomendacja: cykl B6 M9.x re-run

Pełny dedykowany cykl naprawczy `op-newton-momentum-M9-rerun-B6/`:

1. **B6.1** — re-derivacja κ z `√(-g)=c·ψ/(4-3ψ)` w sek08a.
2. **B6.2** — re-derivacja Φ-EOM z poprawnym volume element.
3. **B6.3** — re-extract `c₂` PPN coefficient z M9.1'' linearyzacji.
4. **B6.4** — re-run M9.2 m_field z poprawną metryką.
5. **B6.5** — re-run M9.3 Peters-Mathews z poprawną metryką;
   re-fit GW170817 ograniczenia `|c_GW − c|`.
6. **B6.6** — re-run M9.2 WEP test (MICROSCOPE 10⁻¹⁵).

**Estymata:** 2–3 tygodnie pracy (po S01 zamkniętym strukturalnie).

**Zależność:** S02 NIE może iść bez S01 — bez kanonicznej derywacji
M9.1'' w sek08c, re-run nie ma punktu odniesienia.

## Pliki dotknięte

| Plik | Zakres edycji |
|------|---------------|
| `sek08a_akcja_zunifikowana.tex` | re-derivacja κ z poprawnym √(-g) |
| `sek08_formalizm.tex` | sync Φ-EOM |
| `research/op-newton-momentum/M9_1_pp_P3_results.md` | re-extract c₂ |
| `research/op-newton-momentum/M9_2_results.md` | re-run m_field, WEP |
| `research/op-newton-momentum/M9_3_results.md` | re-run quadrupole, GW170817 |

## Open NEEDS

Patrz [[NEEDS.md]].

## Cross-references

- [[../SUMMARY_2026-05-04.md]] §I.S2
- [[../PRIORITY_MATRIX.md]] klaster A
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.2, §M.5, §O.7
- [[../S01_metric_four_forms]] (zależność wstępna)
- [[../S03_beta_PPN_convention]] (zależność wynikowa)
