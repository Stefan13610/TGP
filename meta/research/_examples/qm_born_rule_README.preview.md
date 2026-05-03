<!--
PODGLĄD (Sesja 2). Pokazuje, jak README.md folderu `research/qm_born_rule/`
WYGLĄDAŁBY po zaaplikowaniu szablonu z meta/research/templates/.

Oryginał: `research/qm_born_rule/README.md` (40 linii, plain text bez frontmattera).

Cel podglądu: zwalidować, czy schemat YAML pasuje do MINIMALNEGO topic-folderu
(brak `program.md`, brak `Phase*_*.md`, jeden README + 0 skryptów). To jest
dolna granica scaffoldingu — folder klasyfikuje się jako `needs-bridge`
albo `sandbox`.

Źródła użyte do wypełnienia:
  - research/qm_born_rule/README.md (jedyny plik w folderze)
  - meta/research/AUDIT_RESEARCH_S1.md (link_total=2, PASS=0, mtime
    2026-04-15, klasa heurystyczna candidate-needs-bridge-or-sandbox)
-->

---
title: "Born Rule from Soliton Tail Interference"
date: 2026-04-15
parent: "[[INDEX.md]]"
related:
  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"
  - "[[../qm_superposition/README.md]]"
  - "[[../qm_foundations/README.md]]"
tags:
  - TGP
  - quantum-mechanics
  - born-rule
  - soliton
tgp_status:
  folder_status: needs-bridge
  level: L1
  kind: derivation
  core_compatibility: unknown
  last_reviewed_against_core: unknown
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: false
  open_bridges:
    - "Soliton tail amplitude → exact |psi|^2 form (not just quadratic scaling)"
    - "Normalization derivation (probabilities sum to 1)"
    - "Cross-term emergence for tail superposition"
    - "Proportionality constant origin"
  depends_on:
    - "research/atom_from_soliton"             # Q1: soliton existence + tail profile
    - "research/qm_superposition"              # weak-field linearity
  impacts:
    - "research/qm_measurement"
    - "research/qm_foundations"
  source_of_status:
    - "README.md §Status: 'Depends on Q1 results. Preliminary analysis suggests quadratic dependence is robust, but full derivation requires confirmed soliton profiles from numerics.'"
    - "Audit S1: PASS=0, CLOSED=0, mtime 2026-04-15 (najstarszy w klastrze qm_*)"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-03
---

# Born Rule from Soliton Tail Interference

## Cel

Wyprowadzić regułę Borna `P = |ψ|²` z dynamiki solitonu w TGP, bez
postulowania jej jako osobnego aksjomatu. Mechanizm proponowany: amplituda
ogona solitonu A_tail gra rolę ψ; całka nakładania ogonów daje
oddziaływanie ~ |A_tail|².

## Stan

**Pre-derivation, L1.** Folder zawiera tylko opis hipotezy (README, 40 linii).
Brak skryptów weryfikujących, brak Phase{1,2,3} struktury. Heurystyka
quadratic scaling jest "robust" jakościowo, ale pełna derivacja wymaga
potwierdzonych profili solitonu z `research/atom_from_soliton/` (Q1).

Source: `README.md §Status` (2026-04-15).

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `README.md` | hipoteza + 4 pytania kluczowe + dependencies | DRAFT |

(folder nie ma jeszcze skryptów, `program.md`, ani `Phase*_*.md`.)

## Wejścia (depends_on)

- [[../atom_from_soliton]] — Q1: soliton existence + tail amplitude profile
- [[../qm_superposition]] — linearność reżimu słabego pola

## Wyjścia (impacts)

- [[../qm_measurement]] — pomiar konsumuje regułę Borna
- [[../qm_foundations]] — fundament QM zawiera Borna

## Otwarte luki

Pełna lista w [[NEEDS.md]]. Headlines (4 pytania kluczowe z README):

- **N1**: czy z ODE solitonu wyprowadzamy dokładną formę |ψ|², czy tylko skalowanie kwadratowe?
- **N2**: jak emerguje normalizacja (Σ P = 1)?
- **N3**: czy wzór interferencji (cross-terms) wychodzi poprawnie dla superpozycji ogonów?
- **N4**: co ustala stałą proporcjonalności?

## Cross-references

- [[meta/PLAN_RESEARCH_WORKFLOW_v1.md]] — workflow
- [[meta/research/CANDIDATE_BRIDGES.md]] — folder jest źródłem 4 otwartych mostów
