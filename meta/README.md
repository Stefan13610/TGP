# Dokumenty meta TGP

Dokumenty meta-poziomu: audyty calosci, plany domkniecia, plany rozwoju, plany publikacji oraz propozycje programow badawczych. Nie zawieraja dowodow — sluza do zarzadzania rozprawa.

## Zawartosc — point-in-time audyty i master plany

- `AUDYT_TGP_2026-04-09.md` — audyt TGP (kwiecien 2026, v1)
- `AUDYT_TGP_2026-04-14_v2.md` — audyt TGP (kwiecien 2026, v2)
- `AUDYT_TGP_2026-05-01.md` — audyt TGP (maj 2026, 8 CRITICAL + 12 HIGH)
- `SUBAGENT_AUDIT_74394a8_2026-05-02.md` — audyt incydentu agentowego (commit 74394a8)
- `ANALIZA_ALPHAK_STATUS_I_DOMKNIECIE_2026-04-12.md` — analiza statusu alpha_K i plan domkniecia
- `ANALIZA_KRYTYCZNA_v6.md` — analiza krytyczna teorii (v6)
- `PLAN_DOMKNIECIA_MASTER.md` — master-plan domkniecia otwartych punktow
- `PLAN_NUMERYCZNY_CG3_CG4.md` — plan numeryczny CG3/CG4 (coarse-graining)
- `PLAN_PUBLIKACJI_MASTER.md` — master-plan publikacji
- `PLAN_ROZWOJU_v3.md` — plan rozwoju TGP (wersja 3)
- `PLAN_ROZWOJU_v4.md` — plan rozwoju TGP (wersja 4)
- `ROADMAP_v3.md` — roadmap (wersja 3)
- `KOIDE_OPEN_PROPOSALS.md` — otwarte propozycje dla sektora Koidego
- `PROPOZYCJA_MOST_GAMMA_DO_PHI.md` — propozycja mostu Gamma -> Phi

## Workflow `research/` + `core/` (od 2026-05-02)

Master plan: [`PLAN_RESEARCH_WORKFLOW_v1.md`](PLAN_RESEARCH_WORKFLOW_v1.md) — wieloфazowy
plan organizacji folderu `research/` z warstwą statusową, multi-agent regułami i bezpieczną
ścieżką promocji do core. Wybrany wariant: **minimalny** (status w YAMLu, bez fizycznych
nadfolderów statusowych). Akceptowane decyzje człowieka: zob. §9 i §9.1 planu.

### `meta/research/` — warstwa workflow dla badań

Pliki **żywe** (aktualizowane przez agentów + skrypty + człowieka):
- [`research/AGENT_PROTOCOL.md`](research/AGENT_PROTOCOL.md) — twarde reguły dla agentów (anti-overclaim, kwarantanna 74394a8)
- [`research/AUDIT_RESEARCH_S1.md`](research/AUDIT_RESEARCH_S1.md) — snapshot audytu Sesji 1 (86 folderów badawczych)
- [`research/RESEARCH_BUS.md`](research/RESEARCH_BUS.md) — tablica ogłoszeń międzyfolderowych
- [`research/CANDIDATE_BRIDGES.md`](research/CANDIDATE_BRIDGES.md) — match `NEEDS` × `FINDINGS`
- [`research/CORE_CANDIDATES.md`](research/CORE_CANDIDATES.md) — kandydaci do awansu do core
- [`research/IMPACT_MATRIX.md`](research/IMPACT_MATRIX.md) — macierz wpływu zmian
- [`research/FOLDER_STATUS_INDEX.md`](research/FOLDER_STATUS_INDEX.md) — globalna mapa statusu (placeholder, wypełnia Sesja 4)
- [`research/MIGRATION_LOG.md`](research/MIGRATION_LOG.md) — log fizycznych ruchów

Pliki **statyczne** (zmienia tylko człowiek):
- [`research/GLOSSARY.md`](research/GLOSSARY.md) — terminy organizacyjne workflow
- [`research/templates/`](research/templates/) — szablony README/NEEDS/FINDINGS/STATUS_BLOCK
- [`research/_examples/`](research/_examples/) — 3 case-study previews (op-uv-as-ngfp, qm_born_rule, closure_2026-04-26)

### `meta/core/` — warstwa workflow dla rdzenia

- [`core/CORE_INVENTORY.md`](core/CORE_INVENTORY.md) — mapa „co jest w core, skąd przyszło"
- [`core/CORE_HOTSPOTS.md`](core/CORE_HOTSPOTS.md) — living checklist 8 CRITICAL + 12 HIGH z `AUDYT_TGP_2026-05-01`
- [`core/CORE_INTAKE.md`](core/CORE_INTAKE.md) — protokół promocji `research/` → `core/` (twardy human gate)
- [`core/intake/`](core/intake/) — wnioski promocyjne `INTAKE_<folder>_<date>.md` (Sesja 7+)

### Najważniejsze decyzje człowieka (zob. PLAN §9, §9.1)

| # | Decyzja |
|---|---------|
| Q-2026-05-02 #1 | Kanoniczny `research/` = `TGP/TGP_v1/research/`. Vault-root `./research/` skasowany 2026-05-03. |
| Q-2026-05-02 #2 | Wariant **minimalny** (status w YAMLu, bez fizycznych nadfolderów statusowych). |
| Q-2026-05-02 #3 | Drugie archiwum `research/_archive/` z rygorystycznymi regułami. |
| Q-2026-05-02 #4 | Warstwy meta: `meta/research/` + `meta/core/`. |
| Q-2026-05-02 #5 | Sesja 4 batch = 10 folderów. |
| Q-2026-05-02 #6 | Agent **może** proponować awans do L4 przez `INTAKE_*.md` (gate ludzki). |
| Q-2026-05-03 #1 | Vault-root `./research/` skasowany. |
| Q-2026-05-03 #2 | Closure-aggregator: tylko rodzic dostaje pełną warstwę YAML/NEEDS/FINDINGS. |
| Q-2026-05-03 #3 | Top-level docs (`TGP_STATUS_*`, `REDIRECT_*`, `NEW_DIRECTIONS_*`) przeniesione do `_archiwum/research_program_docs_2026-04/`. |
| Q-2026-05-03 #4 | Pre-existing `FINDINGS.md`/`NEEDS.md` nie są nadpisywane (`pre_existing: true` flag). |
