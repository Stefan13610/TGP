---
title: "_sandbox — piaskownica dla agentów (NIE część teorii)"
date: 2026-05-03
type: sandbox-readme
parent: "[[../INDEX.md]]"
related:
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
tags:
  - sandbox
  - meta
tgp_status:
  folder_status: sandbox
  level: L0
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

# `_sandbox/` — piaskownica dla agentów

> **Czym jest:** bezpieczna przestrzeń dla luźnych eksploracji
> agentów. Pomysły, hipotezy, "a może by tak…", testy szablonów,
> draft-i nowych folderów badawczych przed promocją do `research/`.
>
> **Czym NIE jest:** częścią aktywnej teorii. Folder `_sandbox/` jest
> wyraźnie odrębny od reszty `research/` — jego zawartość **nie**
> wlicza się do statusu teorii, **nie** generuje predykcji, **nie**
> wpływa na `INDEX.md` / `DEPENDENCIES.md`.

---

## Reguły

1. **Co tu trzymamy:**
   - Eksploracyjne notatki agenta (np. „warto by sprawdzić X")
   - Draft-y nowych folderów badawczych przed promocją
   - Eksperymenty z szablonami `meta/research/templates/`
   - Tymczasowe scratch-pady multi-agent rozmów

2. **Czego tu NIE trzymamy:**
   - Wyników, które są już w pełnym `research/<folder>/`
     (jeśli wynik dojrzał — folder przenosimy do `research/`, nie
     duplikujemy)
   - Plików z otwartymi zależnościami (`depends_on:`) wskazujących
     na `_sandbox/` z innych folderów (zależność na sandbox =
     sygnał, że folder powinien wyjść z sandbox)
   - Czegokolwiek z `level: L1` lub wyższym — to oznacza, że folder
     jest gotowy na `research/<folder>/`

3. **Awansowanie z `_sandbox/` do `research/`:**
   - Folder dostaje `tgp_status` zgodny z [[meta/research/templates/STATUS_BLOCK.yaml]]
   - Wpisuje się do `MIGRATION_LOG.md`
   - Reorganizujemy linki, jeśli były (zwykle nie ma — sandbox jest
     izolowany)

4. **Anty-overclaim w sandbox:**
   - **NIE wolno** używać "DERIVED", "PASS", "CLOSED" bez
     weryfikowalnego pliku — nawet w sandboxie. Reguła protokołu jest
     globalna.
   - Sandbox nie jest "no rules zone".

5. **Maintenance:**
   - Folder `_sandbox/` jest scanowany w Sesji 8 (regression).
   - Pliki bez aktywności > 90 dni są oznaczane do przeglądu.
   - Człowiek decyduje, czy promować do `research/` czy skasować.

---

## Aktualna zawartość

(folder utworzony 2026-05-03 — pusty)

---

## Cross-references

- [[meta/PLAN_RESEARCH_WORKFLOW_v1.md]] §3.3 — definicja sandbox w wariancie minimalnym
- [[meta/research/AGENT_PROTOCOL.md]] — globalne reguły agentów (obowiązują też tu)
- [[../_archive/README.md]] — komplementarny folder dla obsoletnych podejść
