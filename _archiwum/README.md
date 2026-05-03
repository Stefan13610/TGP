---
title: "_archiwum — frozen historical drafts and superseded plans"
date: 2026-05-04
type: archive-readme
status: FROZEN — historical record only
parent: "[[../INDEX.md]]"
related:
  - "[[../meta/PLAN_DOMKNIECIA_MASTER.md]] (current closure plan)"
  - "[[../meta/PLAN_ROZWOJU_v4.md]] (current development roadmap)"
  - "[[../meta/PLAN_PUBLIKACJI_MASTER.md]] (current publication plan)"
  - "[[../meta/AUDYT_TGP_2026-05-01.md]]"
  - "[[../meta/PLAN_RESEARCH_WORKFLOW_v1.md]]"
tags:
  - archive
  - frozen
  - historical
  - meta
---

# `_archiwum/` — frozen historical drafts

> **STATUS: FROZEN 2026-05-04. Read-only historical record.**
>
> **Do NOT edit, do NOT promote into current cycles.** All files in this
> folder are **superseded** by current canonical documents in
> [[../meta/]] root and [[../README.md]].

## Why this folder exists

Earlier development cycles (2026-03 → 2026-04-19) accumulated draft plans,
session summaries, and consistency analyses. After the closure_2026-04-26
+ M10 + audit_2026-05-01 + research workflow v1 (Sesja 1-8) cleanup, these
files were moved to `_archiwum/` to keep the working repo clean while
preserving audit trail.

## Mapping: superseded → current

### Top-level (`_archiwum/*.md`)

| Archived file/family | Status | Superseded by |
|---|---|---|
| `PLAN_ROZWOJU_v1.md` / `v2.md` | superseded | [[../meta/PLAN_ROZWOJU_v4.md]] |
| `PLAN_DOMKNIECIA.md` / `_v1.md` | superseded | [[../meta/PLAN_DOMKNIECIA_MASTER.md]] |
| `PLAN_AKTUALIZACJI_v26.md` | superseded | covered in [[../meta/PLAN_ROZWOJU_v4.md]] |
| `ANALIZA_BRAKUJACYCH_ELEMENTOW.md` | superseded | [[../meta/AUDYT_TGP_2026-05-01.md]] §A-D |
| `ANALIZA_SPRZEZONY.md` | superseded | covered in current sek08* |
| `KOIDE_OPEN_PROPOSALS_v47b_archive.md` | superseded | `old_meta_2026-04/KOIDE_OPEN_PROPOSALS.md` (also archived 2026-05-04) |
| `WYPROWADZENIE_r21.md` | historical | research/zero_mode + partial_proofs/zero_mode |
| `tgp_letter_draft.md` | superseded | [[../tgp_letter.tex]] (PRL format final) |
| `SESSION_v35.md` / `v37.md` | superseded | covered in canonical audits |
| `PLAN_ANALITYCZNY_KOIDE.md` | superseded | covered via cycles θ.1, ζ.1, κ.1 |

### `_archiwum/old_meta_2026-04/` (archived 2026-05-04)

Plans/audits/analyses from April 2026 that were superseded by AUDYT_TGP_2026-05-01
+ PLAN_ROZWOJU_v4 + PLAN_DOMKNIECIA_MASTER cluster:

| Archived file | Superseded by |
|---|---|
| `AUDYT_TGP_2026-04-09.md` | [[../meta/AUDYT_TGP_2026-05-01.md]] |
| `AUDYT_TGP_2026-04-14_v2.md` | [[../meta/AUDYT_TGP_2026-05-01.md]] |
| `ANALIZA_ALPHAK_STATUS_I_DOMKNIECIE_2026-04-12.md` | cycles ε.1 + α.1 (closed) |
| `ANALIZA_KRYTYCZNA_v6.md` | [[../meta/AUDYT_TGP_2026-05-01.md]] |
| `PLAN_NUMERYCZNY_CG3_CG4.md` | [[../research/continuum_limit/]] (closed) |
| `PLAN_ROZWOJU_v3.md` | [[../meta/PLAN_ROZWOJU_v4.md]] |
| `ROADMAP_v3.md` | [[../meta/PLAN_ROZWOJU_v4.md]] + [[../meta/AUDYT_TGP_2026-05-01.md]] §B-cluster (historical line-number citations preserved verbatim in audit) |
| `KOIDE_OPEN_PROPOSALS.md` | open issues addressed via cycles θ.1, ζ.1, η.1, η.2, κ.1, ι.1, ν.1 |
| `PROPOZYCJA_MOST_GAMMA_DO_PHI.md` | [[../partial_proofs/most_gamma_phi/dodatekQ2_most_gamma_phi_lematy.tex]] |

### `_archiwum/research_program_docs_2026-04/` (archived earlier)

| Archived file | Superseded by |
|---|---|
| `REDIRECT_PROGRAM_2026-04-19.md` | [[../meta/PLAN_RESEARCH_WORKFLOW_v1.md]] |
| `TGP_STATUS_2026-04-18.md` / `2026-04-19.md` | [[../README.md]] header + [[../INDEX.md]] |
| `NEW_DIRECTIONS_2026-04-20.md` | [[../meta/PLAN_ROZWOJU_v4.md]] |

### Subfolders

| Folder | Status |
|---|---|
| `stare_analizy/` (ANALIZA_SPOJNOSCI_v1-v5, SESSION_v36-v40, etc.) | superseded; covered in canonical audits |
| `scripts_exploratory/` | exploratory only; promoted scripts moved to research/ or scripts/ |
| `OJ3_zasada_selekcji_g0_tau.md` | historical | covered in `op-theta-quark-koide/` + `op-zeta-mass-spectrum/` |

## Reading rules for agents

1. **Do NOT edit any file in `_archiwum/`** — these are immutable historical
   record.
2. **Do NOT cite `_archiwum/` files as authoritative** — always cite the
   canonical replacement (per table above).
3. **Do NOT promote** any `_archiwum/` content back into research cycles
   without explicit operator approval.
4. **OK to read** for historical context (e.g. understanding why a
   particular formulation was abandoned).
5. **OK to reference** in audit cross-cites if explaining a transition
   (e.g. "previously X, now Y per …").

## Migration log reference

Folder contents were migrated per [[../meta/research/MIGRATION_LOG.md]]
during Sesja 2-3 (2026-05-02 → 2026-05-03) of the research workflow v1
deployment.

---

**Status:** FROZEN — frozen-archive section of TGP_v1 repository.
**Last sweep:** 2026-05-04 (CALIBRATION audit cleanup phase).
