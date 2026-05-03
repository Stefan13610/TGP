---
title: "_archive — archiwum obsoletnych ścieżek badawczych"
date: 2026-05-03
type: archive-readme
parent: "[[../INDEX.md]]"
related:
  - "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/research/MIGRATION_LOG.md]]"
  - "[[../../_archiwum/README.md]]"
tags:
  - archive
  - meta
tgp_status:
  folder_status: archive
  level: unknown
  kind: program-doc
  core_compatibility: unknown
  last_reviewed_against_core: unknown
  may_edit_core: false
  exports_findings: false
  has_needs_file: false
  has_findings_file: false
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "Sesja 3 — utworzenie folderu zgodnie z PLAN_RESEARCH_WORKFLOW_v1.md §3.3"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-03
---

# `_archive/` — archiwum obsoletnych ścieżek badawczych

> **Czym jest:** miejsce na **stare ścieżki badawcze / próby
> derivacyjne**, które są oznaczone jako WITHDRAWN, FALSIFIED,
> OBSOLETE lub SUPERSEDED. Folder zachowany **jako historia**, nie jako
> żywy research.
>
> **Różnica vs `_archiwum/` (TGP_v1 level):**
>
> | Folder | Co tu trafia | Przykład |
> |--------|--------------|----------|
> | `research/_archive/` | **Stare ścieżki badawcze** (research-level) | obsoletny `op-X/` przed pivot |
> | `_archiwum/` (TGP_v1) | **Stare plany / sesje / analizy meta** (meta-level) | `PLAN_DOMKNIECIA_v1.md`, `SESSION_v37.md`, `research_program_docs_2026-04/` |

---

## Reguły wpisu (sekcja 3.4 PLAN)

Folder trafia do `_archive/` **wyłącznie** gdy spełnia jeden z warunków:

1. Treść jest oznaczona w samym folderze jako `WITHDRAWN`, `FALSIFIED`,
   `OBSOLETE`, `SUPERSEDED` w sposób potwierdzony **co najmniej dwoma
   źródłami** (np. notka w `README.md` + wpis w `KNOWN_ISSUES.md`
   lub w jakimś `AUDYT_*.md`).
2. **Lub** człowiek explicite zdecydował o archiwizacji.

**Nie trafiają tu:**
- foldery aktualnie używane jako historia/uzasadnienie (zostają in
  place z `folder_status: core-promoted`),
- foldery z otwartym researchem (zostają in place z
  `active`/`needs-bridge`),
- foldery typu `closure_2026-04-26/` (closure-aggregator, in place),
- legacy `op6/`, `op7/`, `op1-op2-op4/` (zbyt gęsto linkowane,
  in place).

---

## Procedura archiwizacji (Sesja 8.5)

Każde przeniesienie do `_archive/` jest **fizycznym `git mv`** + patchem
linków + wpisem do `MIGRATION_LOG.md` z akceptacją człowieka **per
folder**:

```
git mv research/<folder> research/_archive/<folder>
# + patch linków w *.md, *.tex, *.py
# + update tgp_status w README:
#     folder_status: archive
#     archived_date: YYYY-MM-DD
#     archived_reason: "<quote z 2 źródeł>"
#     archived_by_human_decision: YYYY-MM-DD
# + git commit (osobny per folder)
# + sanity check: pdflatex main.tex && tooling/build_deps_graph.py
```

Pełna procedura: [[meta/PLAN_RESEARCH_WORKFLOW_v1.md]] §6 Sesja 8.5.

---

## Anty-overclaim w `_archive/`

Foldery w `_archive/` zachowują swój oryginalny content. **NIE poprawiamy**
historycznych werdyktów. Jeśli stary folder zawiera "FULL CONVERGENCE
18/18 PASS" — zostaje, **z dopisanym** w README header'em:

```yaml
tgp_status:
  folder_status: archive
  archived_date: YYYY-MM-DD
  archived_reason: "<powód>"
  pre_archive_status_invalidated: true   # informacja, że stary status w treści NIE odzwierciedla aktualnego stanu
```

Reguła: "stary tekst zostaje, status YAML mówi, że jest stary".

---

## Aktualna zawartość

(folder utworzony 2026-05-03 — pusty)

Pierwsza fala kandydatów do archiwizacji **prawdopodobnie nie wystąpi**
przed Sesją 8.5. Z audytu Sesji 1: 6 folderów w klasie heurystycznej
"candidate-archive (verify)" — wszystkie z PASS > 0, więc to false
positive na heurystykę słów-kluczy `WITHDRAWN`/`SUPERSEDED` używanych
w wewnętrznym wersjonowaniu (np. "v1 WITHDRAWN, replaced by v2").
Sesja 4 zweryfikuje ręcznie.

---

## Cross-references

- [[meta/PLAN_RESEARCH_WORKFLOW_v1.md]] §3.4 — reguły wpisu
- [[meta/research/MIGRATION_LOG.md]] — log fizycznych ruchów
- [[meta/research/AGENT_PROTOCOL.md]] — globalne reguły
- [[../_sandbox/README.md]] — komplementarny folder dla luźnych eksploracji
- [[../../_archiwum/README.md]] — meta-archiwum (TGP_v1 level)
