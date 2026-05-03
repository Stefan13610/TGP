<!--
README.template.md
Szablon README.md dla folderu w `research/<name>/`.

Reguły wypełniania:
  1. Frontmatter `tgp_status:` jest OBOWIĄZKOWY. Schemat:
     [[meta/research/templates/STATUS_BLOCK.yaml]].
  2. Jeśli folder MA już inne pola frontmattera (title/date/tags/...) —
     zachowujemy je, dodajemy `tgp_status:` jako kolejne pole. NIE
     nadpisujemy istniejących pól.
  3. Sekcja "## Cel" jest obowiązkowa (1–3 zdania, czego dotyczy folder).
  4. Sekcja "## Stan" cytuje weryfikowalne źródło (anty-overclaim).
  5. Pozostałe sekcje są SUGEROWANE. Jeśli folder ma już rozbudowany README,
     nie nadpisujemy treści — agent dopisuje tylko brakujące sekcje
     (Cel / Stan / Pliki / Dependencies).
  6. Wikilinks używamy z formą `[[research/<inny-folder>/...]]` lub
     bez ścieżki dla portable plików.
-->

---
title: "<TYTUŁ FOLDERU — krótki, zgodny z konwencją>"
date: <YYYY-MM-DD utworzenia / istotnej aktualizacji>
parent: "[[INDEX.md]]"
related:
  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"
tags:
  - TGP
  # - <topic-tags>
tgp_status:
  folder_status: <active | core-ready | core-promoted | needs-bridge | needs-migration | audit | review | program-doc | sandbox | archive>
  level: <L0 | L1 | L2 | L3 | L4 | mixed | unknown>
  kind: <derivation | numerical | phenomenology | bridge | audit | review | program-doc | closure-aggregator>
  core_compatibility: <current | partial | stale | broken | unknown>
  last_reviewed_against_core: <YYYY-MM-DD | unknown>
  may_edit_core: false
  exports_findings: <true | false>
  has_needs_file: <true | false>
  has_findings_file: <true | false>
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status: []
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: <YYYY-MM-DD>
---

# <TYTUŁ FOLDERU>

## Cel

<1–3 zdania: jakiego problemu dotyczy folder, jaką lukę domyka.>

## Stan

<Jednolinijkowy werdykt + cytat weryfikowalnego źródła.>

Przykłady poprawnych zapisów:
- "Phase 3 zamknięte 94/97 PASS, 2026-04-25 (zob. `Phase3_results.md`)."
- "Numeryka istnieje (3/3 PASS w `m2b_results.txt`), brak analitycznego mostu — patrz [[NEEDS.md]]."
- "Polluted incydentem 74394a8; status zamrożony do forward-patch
  (zob. [[meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]])."

Przykłady nadużyć (NIE wolno):
- "FULL CONVERGENCE 18/18" bez podanego pliku z PASS-ami.
- "DERIVED" w stanie, gdy `level: L2`.
- "Wszystko działa" bez liczby + pliku.

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `program.md` | <opis cyklu / planu fazowego>  | <ACTIVE / CLOSED / DRAFT> |
| `Phase1_setup.md` / `Phase1_results.md` | <faza 1>  | <—> |
| <skrypty>  | <co robią>  | <—> |

## Wejścia (depends_on)

<Foldery / pliki, których wyników ten folder potrzebuje. Każdy z linkiem.>

- [[research/<inny-folder>/FINDINGS.md]] — co dokładnie jest brane na wejściu

## Wyjścia (impacts)

<Foldery / sekcje core, które mogą skonsumować wyniki tego folderu.>

- [[research/<konsument-folder>]] — co eksportujemy
- (jeśli core-ready / core-promoted) `core/sekX_*` — referencja docelowa

## Otwarte luki

Pełna lista w [[NEEDS.md]]. Tu tylko jednolinijkowy headline:

- <luka 1>
- <luka 2>

## Cross-references

- [[meta/PLAN_RESEARCH_WORKFLOW_v1.md]] — workflow
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalny indeks
- [[meta/research/RESEARCH_BUS.md]] — broadcast
- (opcjonalne) [[meta/AUDYT_TGP_2026-05-01.md]] / inne audyty
