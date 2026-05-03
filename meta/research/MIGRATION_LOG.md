---
title: "MIGRATION_LOG — log fizycznych ruchów / promocji / kasowań w research/ + meta/"
date: 2026-05-03
type: log
status: LIVING (Sesja 3 initial seed)
parent: "[[meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
related:
  - "[[meta/research/AGENT_PROTOCOL.md]]"
  - "[[meta/research/CORE_CANDIDATES.md]]"
  - "[[meta/core/CORE_INTAKE.md]]"
tags:
  - log
  - migration
  - audit-trail
---

# MIGRATION_LOG — log fizycznych ruchów

> **Cel pliku:** chronologiczny ślad każdej fizycznej zmiany lokalizacji
> folderu / pliku w `research/` lub `meta/`. Punkt audytowy dla człowieka
> i dla regression review w Sesji 8.
>
> **Każdy wpis** ma: datę, typ akcji, ścieżki source/target, ID akceptacji
> człowieka (jeśli dotyczy), liczbę zaktualizowanych linków, hash commita
> (jeśli git).

---

## 0. Reguły wpisu

1. **Trigger**: każde `git mv` / `shutil.move` / `rmtree` / utworzenie
   nowego fizycznego folderu w `research/` lub `meta/`.
2. **NIE** loguje się: edycji plików (to jest w git history).
   Tylko ruchy strukturalne.
3. **Format**: ID rosnąco (`MIG-001`, `MIG-002`, ...), data, typ akcji.
4. **Akceptacja człowieka** (`human_accepted_at:`) jest WYMAGANA dla:
   - Fizycznego `mv` w Sesji 8.5
   - Kasowania folderu / pliku
   - Awansu folderu do `_archive/`
5. **Bez akceptacji człowieka** dopuszczalne (przez agenta):
   - Tworzenie nowych pomocniczych folderów (`_sandbox`, `_archive`,
     `_examples`, `templates`, `intake`)
   - Tworzenie nowych pomocniczych plików `meta/research/*` i
     `meta/core/*` (Sesja 3)
6. **Link patches**: każdy `mv` ma listę `links_updated` z liczbą i
   plikami źródłowymi.

---

## 1. Statystyka

| Typ akcji | Liczba |
|-----------|-------:|
| `delete` | 1 |
| `mv` (file) | 4 |
| `mv` (folder) | 0 |
| `mkdir` (new structural folder) | 4 |
| `mkdir` (new helper folder) | 5 |
| `git mv` (Sesja 8.5) | 0 |
| `archive` (do `_archive/`) | 0 |
| **Razem** | **14** |

Ostatnia aktualizacja: 2026-05-03.

---

## 2. Log wpisów (chronologicznie)

### MIG-001 — 2026-05-03 — `delete` ./research/ (vault root)

- **Action**: `shutil.rmtree`
- **Source**: `./research/` (vault root, OUTSIDE `TGP/TGP_v1/`)
- **Target**: (deleted)
- **Reason**: Decyzja człowieka Q1 z 2026-05-03 (zob.
  [[PLAN_RESEARCH_WORKFLOW_v1.md]] §9.1). Vault-root `./research/`
  zawierał 11 stub-folderów (14 KB), tylko stub `phase1_*_audit.txt`
  ze stderr nieudanych runów. Brak unikalnej treści — zaweryfikowane
  przez `_audit_s1.py vault_diff` (only_in_vault: 0 dla każdego folderu
  wspólnego).
- **Human accepted at**: 2026-05-03 (user message: "to jest do skasowania")
- **Links updated**: 0 (vault-root `./research/` nie był linkowany z
  `TGP/TGP_v1/`).
- **Side effects**: `_audit_s1_raw.json` po regeneracji
  (vault_research_exists: false).

### MIG-002 → MIG-005 — 2026-05-03 — `mv` 4 stale top-level docs

- **Action**: `shutil.move` (4 plików)
- **Source files**:
  - `research/TGP_STATUS_2026-04-18.md`
  - `research/TGP_STATUS_2026-04-19.md`
  - `research/REDIRECT_PROGRAM_2026-04-19.md`
  - `research/NEW_DIRECTIONS_2026-04-20.md`
- **Target**: `_archiwum/research_program_docs_2026-04/`
- **Reason**: Decyzja człowieka Q3 z 2026-05-03 (top-level docs
  outdated, no point updating, move to archive). Wszystkie 4 są
  point-in-time research-program snapshots z 2026-04, superseded
  przez `meta/AUDYT_TGP_2026-05-01.md`.
- **Human accepted at**: 2026-05-03
- **Links updated**: 1 plik
  - `meta/PLAN_RESEARCH_WORKFLOW_v1.md` — 2 eksplicytne `[[research/...]]`
    paths zmienione na `[[_archiwum/research_program_docs_2026-04/...]]`
- **Wikilinks not updated**: 7 plików w research/ używa wikilinków
  bez ścieżki (`[[FILENAME.md]]`); Obsidian rozwiązuje po nazwie pliku,
  więc działają nadal:
  - `research/atomic_shells_closure/ATOMIC_SHELLS_VERDICT.md`
  - `research/desi_dark_energy/README.md`
  - `research/liquid_viscosity/SUMMARY.md`
  - `research/galaxy_scaling/CLOSURE_2026-04-19.md`
  - `research/casimir_mof/SUMMARY.md`
  - `research/particle_sector_closure/README.md`
  - `research/muon_g_minus_2/SUMMARY.md`

### MIG-006 — 2026-05-03 — `mkdir` `meta/research/`

- **Action**: `mkdir -p`
- **Target**: `TGP/TGP_v1/meta/research/`
- **Reason**: Sesja 1 audytu — utworzenie warstwy meta dla badań.
- **Human accepted at**: 2026-05-03 (user "działaj")

### MIG-007 — 2026-05-03 — `mkdir` `meta/research/templates/`

- **Action**: `mkdir -p`
- **Target**: `TGP/TGP_v1/meta/research/templates/`
- **Reason**: Sesja 2 — szablony.
- **Human accepted at**: 2026-05-03

### MIG-008 — 2026-05-03 — `mkdir` `meta/research/_examples/`

- **Action**: `mkdir -p`
- **Target**: `TGP/TGP_v1/meta/research/_examples/`
- **Reason**: Sesja 2 — 3 case-study previews.
- **Human accepted at**: 2026-05-03

### MIG-009 — 2026-05-03 — `mkdir` `meta/core/`

- **Action**: `mkdir -p`
- **Target**: `TGP/TGP_v1/meta/core/`
- **Reason**: Sesja 3 — warstwa meta dla rdzenia (decyzja człowieka
  #4 z 2026-05-02: `meta/research/` + `meta/core/`).
- **Human accepted at**: 2026-05-03

### MIG-010 — 2026-05-03 — `mkdir` `meta/core/intake/`

- **Action**: `mkdir -p` + `.gitkeep`
- **Target**: `TGP/TGP_v1/meta/core/intake/`
- **Reason**: Sesja 3 — folder na wnioski INTAKE.
- **Human accepted at**: 2026-05-03

### MIG-011 — 2026-05-03 — `mkdir` `_archiwum/research_program_docs_2026-04/`

- **Action**: `mkdir -p`
- **Target**: `TGP/TGP_v1/_archiwum/research_program_docs_2026-04/`
- **Reason**: Cel ruchu MIG-002→MIG-005.
- **Human accepted at**: 2026-05-03

### MIG-012 — 2026-05-03 — `mkdir` `research/_sandbox/` (Sesja 3)

- **Action**: `mkdir -p` + `README.md`
- **Target**: `TGP/TGP_v1/research/_sandbox/`
- **Reason**: Sesja 3 — sandbox dla agentów (decyzja człowieka #2/#3
  z 2026-05-02).
- **Human accepted at**: 2026-05-03

### MIG-013 — 2026-05-03 — `mkdir` `research/_archive/` (Sesja 3)

- **Action**: `mkdir -p` + `README.md`
- **Target**: `TGP/TGP_v1/research/_archive/`
- **Reason**: Sesja 3 — research-level archiwum (decyzja człowieka #3
  z 2026-05-02). Różne od `_archiwum/` (TGP_v1 level — meta-archiwum).
- **Human accepted at**: 2026-05-03

---

## 3. Reguły wpisu (do Sesji 8.5+)

Każdy `git mv` w Sesji 8.5 (selektywna archiwizacja) musi mieć:

```
### MIG-NNN — YYYY-MM-DD — `git mv` <folder> → `_archive/`

- **Action**: `git mv research/<folder> research/_archive/<folder>`
- **Source folder**: research/<folder>
- **Target folder**: research/_archive/<folder>
- **Reason**: <2 źródła zgodnie z sekcją 3.4 PLAN — np.
  "WITHDRAWN per [[../KNOWN_ISSUES.md]] §X" + "OBSOLETE per
  [[meta/AUDYT_TGP_*.md]] §Y">
- **Human accepted at**: YYYY-MM-DD (per-folder akceptacja)
- **Links updated**: <N> plików (lista poniżej)
  - `path/to/file.md` — <count> referencji zmieniono z
    `research/<folder>/...` na `research/_archive/<folder>/...`
- **YAML changes**: research/_archive/<folder>/README.md:
  - `folder_status: archive`
  - `archived_date: YYYY-MM-DD`
  - `archived_reason: <quote>`
  - `archived_by_human_decision: YYYY-MM-DD`
- **Sanity checks passed**:
  - [ ] `pdflatex main.tex && pdflatex main.tex` → OK
  - [ ] `tooling/build_deps_graph.py` → orphans=0
  - [ ] grep szczątkowych referencji w *.txt / *.py → 0
- **Commit hash**: <git short SHA>
```
