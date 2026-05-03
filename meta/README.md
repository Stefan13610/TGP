# Dokumenty meta TGP

Dokumenty meta-poziomu: audyty calosci, plany domkniecia, plany rozwoju, plany publikacji oraz propozycje programow badawczych. Nie zawieraja dowodow — sluza do zarzadzania rozprawa.

## Zawartosc — aktywne audyty i master plany (post-2026-05-04 cleanup)

### Audyty (canonical)
- `AUDYT_TGP_2026-05-01.md` — canonical audit (maj 2026, 8 CRITICAL + 12 HIGH, all closed §J-§AB)
- `SUBAGENT_AUDIT_74394a8_2026-05-02.md` — audyt incydentu agentowego (CLOSED 2026-05-04 §8 forward-patch)
- `CALIBRATION_PROTOCOL.md` — anti-overclaim discipline (BINDING 2026-05-04+)

### Master plany (active)
- `PLAN_DOMKNIECIA_MASTER.md` — master-plan domkniecia otwartych punktow (10/10 closed)
- `PLAN_PUBLIKACJI_MASTER.md` — master-plan publikacji
- `PLAN_ROZWOJU_v4.md` — plan rozwoju TGP (current dev roadmap)
- `PLAN_RESEARCH_WORKFLOW_v1.md` — workflow `research/` (multi-agent + status layer)

### Zarchiwizowane 2026-05-04 (zob. [`../_archiwum/old_meta_2026-04/`](../_archiwum/old_meta_2026-04/))

Następujące plany/analizy/audyty z kwietnia 2026 zostały zarchiwizowane jako
**superseded** przez canonical AUDYT_TGP_2026-05-01 + PLAN_ROZWOJU_v4 +
PLAN_DOMKNIECIA_MASTER:
- `AUDYT_TGP_2026-04-09.md` (superseded by AUDYT_2026-05-01)
- `AUDYT_TGP_2026-04-14_v2.md` (superseded by AUDYT_2026-05-01)
- `ANALIZA_ALPHAK_STATUS_I_DOMKNIECIE_2026-04-12.md` (alpha_K closed via cycles ε.1+α.1)
- `ANALIZA_KRYTYCZNA_v6.md` (superseded by AUDYT_2026-05-01)
- `PLAN_NUMERYCZNY_CG3_CG4.md` (CG3/CG4 done w `research/continuum_limit/`)
- `PLAN_ROZWOJU_v3.md` (superseded by PLAN_ROZWOJU_v4)
- `ROADMAP_v3.md` (superseded by PLAN_ROZWOJU_v4 + AUDYT_2026-05-01 §B-cluster;
  cytowane historycznie w 6× AUDYT references + 5× tooling/scripts comments)
- `KOIDE_OPEN_PROPOSALS.md` (open issues addressed via cycles θ.1, ζ.1, η.1, η.2, κ.1, ι.1)
- `PROPOZYCJA_MOST_GAMMA_DO_PHI.md` (propozycja zrealizowana w `partial_proofs/most_gamma_phi/dodatekQ2_*.tex`)

> **Reading rule:** zarchiwizowane pliki są **read-only historical record**.
> Nie cytować jako authoritative — używać canonical replacements. Patrz
> [`../_archiwum/README.md`](../_archiwum/README.md) §"Mapping: superseded → current".

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
