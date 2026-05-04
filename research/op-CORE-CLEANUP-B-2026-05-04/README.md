---
title: "CORE-CLEANUP-B — hybrydowe wyczyszczenie rdzenia (literal cleanup body + dodatek_pivot_history)"
date: 2026-05-04
cycle: CORE-CLEANUP-B
type: cleanup-cycle
status: in_progress (multi-session)
parent: "[[../../audyt/CLEANUP_INVENTORY_2026-05-04.md]]"
predecessors:
  - "[[../op-g0-r3-from-canonical-projection/README.md]]" (G.0 closure 2026-05-02)
  - "[[../op-L01-rho-stress-energy-bridge-2026-05-04]]" (L01 EXECUTED)
  - "[[../op-L04-ODE-canonicalization-2026-05-04]]" (L04 RESOLVED)
related:
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]"
  - "[[../../main.tex]]"
tags:
  - TGP
  - core-cleanup
  - literal-cleanup
  - pivot-history
  - publication-ready
  - non-breaking
tgp_status:
  folder_status: active
  level: L4
  kind: derivation
  core_compatibility: in_progress
  last_reviewed_against_core: 2026-05-04
  may_edit_core: true
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: ["sek08a-cleanup", "sek08c-cleanup", "cross-ref-audit"]
  depends_on:
    - "[[../op-g0-r3-from-canonical-projection]]"
    - "[[../op-L04-ODE-canonicalization-2026-05-04]]"
  impacts:
    - "[[../../core/sek08a_akcja_zunifikowana]]"
    - "[[../../core/sek08c_metryka_z_substratu]]"
    - "[[../../core/formalizm/]]"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# CORE-CLEANUP-B — hybrydowe wyczyszczenie rdzenia

## Cel

Wyczyścić rdzeń LaTeX (sek08a, sek08c) z deprecated propositions/equations
przenosząc je do **dedicated appendix** `dodatek_pivot_history.tex` z explicit
`[HISTORYCZNE — FALSIFIED/OBSOLETE]` markerami, zachowując **wszystkie labels**
dla cross-reference compatibility.

## Zasada nadrzędna: NON-BREAKING

- ✓ **Wszystkie istniejące labels zachowane** w appendix (cross-refs działają)
- ✓ **Żadna derywacja merytorycznie usunięta** — tylko relokacja z body do appendix
- ✓ **pdflatex compile clean** po każdej fazie
- ✓ **Body czysty** — czytelnik widzi tylko kanoniczne G.0+L01+L04 formuły
- ✗ Brak modyfikacji propositions/equations (struktura zachowana, tylko miejsce)

## Plan wielofazowy

### Phase 0 — setup (TA SESJA)

- ✓ Utworzony cykl folder
- ✓ PRE_PATCH_SNAPSHOT (hash plików rdzenia)
- ✓ PLAN szczegółowy (per-file lista deprecated propositions)
- ✓ Inwentaryzacja cross-references (skopiowana z P33 audit G.0)

### Phase 1 — `dodatek_pivot_history.tex` skeleton (TA SESJA)

- ✓ Utworzyć `core/formalizm/dodatek_pivot_history.tex`
- ✓ Dodać include do `main.tex` po dodatekV (lin. 90)
- ✓ Skeleton z 3 sekcjami:
  - Section 1: Pre-G.0 sek08a v1.x (V_orig, √(-g)=c·φ, κ z M9.1)
  - Section 2: Pre-G.0 sek08c forms (I), (II), (III)
  - Section 3: Pre-Phase 2 mass formulas LP-4 / R5 K² (specjalne case α=1)
- ✓ pdflatex compile check (skeleton sam compile clean)

### Phase 2 — sek08c body cleanup (NEXT SESSION lub w tej, jeśli czas)

Mniejszy plik (592 linii), 4 formy metryki, cleaner refactor:

- 2.1: Skopiować formy (I), (II), (III) z labels do `dodatek_pivot_history.tex` Section 2
- 2.2: Usunąć z body `eq:metric-full-derived` (lin. ~131-200), formę eksponencjalną
  (lin. ~208-211), `thm:antipodal-uniqueness` (lin. ~345-371)
- 2.3: Body refactor: tylko M9.1'' kanoniczna derivation (z preamble notki +
  TGP_FOUNDATIONS:64-69 expanded)
- 2.4: pdflatex compile check po każdym kroku

### Phase 3 — sek08a body cleanup (NEXT SESSION)

Większy plik (1021 linii), v1.x → v2.0 reorganizacja:

- 3.1: Skopiować v1.x propositions/equations z labels do `dodatek_pivot_history.tex` Section 1
- 3.2: Usunąć z body deprecated:
  - `eq:V-selfinterference` (lin. 83-93) → V_orig (β/3)φ³−(γ/4)φ⁴
  - `eq:sqrt-g-eff` (lin. 141-145) → √(-g)=c·φ
  - `prop:vacuum-condition` (β=γ z M9.1)
  - `prop:kappa-corrected` (κ z falsified √(-g))
  - inne deprecated v1.x propositions
- 3.3: Body refactor: użycie wprost v2.0 propositions:
  - `prop:V-M911-canonical` jako główny
  - `eq:sqrt-g-eff = c·ψ/(4-3ψ)` jako kanoniczne
  - `prop:kappa-corrected-G0` jako kanoniczne
- 3.4: pdflatex compile check

### Phase 4 — cross-ref audit (NEXT SESSION)

~30 HIGH-impact files (P33 audit z op-g0-r3):

- 4.1: Sprawdzić wszystkie `\ref{}` do deprecated labels
- 4.2: Większość powinna nadal działać (labels w appendix)
- 4.3: Update gdzie potrzebne — referencja do new G.0 propositions
- 4.4: Weryfikacja podczas pdflatex compile

### Phase 5 — compile check + finalization (NEXT SESSION)

- 5.1: Pełny `pdflatex main` + `bibtex` + `pdflatex` × 2
- 5.2: Sprawdzić: 0 nowych undefined references, 0 nowych multiply-defined labels
- 5.3: Strony porównanie (G.0 = 537 stron clean; cleanup target = ≤537 stron)
- 5.4: POST_PATCH_VERIFICATION + final summary

## Estymata

| Faza | Czas | Status sesja 2026-05-04 |
|------|------|--------------------------|
| 0 — Setup | 30 min | **DONE** |
| 1 — Skeleton dodatek | 1h | **DOING** |
| 2 — sek08c cleanup | 2-3h | NEXT |
| 3 — sek08a cleanup | 4-6h | NEXT |
| 4 — Cross-ref audit | 2-3h | NEXT |
| 5 — Compile check | 1h | NEXT |

**Razem estymowane:** 10-15h pracy. W jednej sesji można zrobić Phase 0+1 +
zacząć Phase 2.

## Status sesji 2026-05-04

**Wykonane w tej sesji:**

1. Phase 0: setup folder + plan
2. Phase 1: skeleton `dodatek_pivot_history.tex`
3. Phase 2 partial: sek08c forms (I), (II), (III) → appendix (jeśli czas)

**Pozostałe sesje:**

- Phase 2 finish + Phase 3 (sek08a) + Phase 4 (cross-ref) + Phase 5 (compile)

## Pliki w cyklu

| Plik | Opis | Status |
|------|------|--------|
| [[README.md]] | (ten plik) | DONE |
| [[PLAN.md]] | per-file lista deprecated propositions | TODO |
| [[PRE_PATCH_SNAPSHOT.md]] | hashes + line counts pre-cleanup | TODO |
| [[POST_PATCH_VERIFICATION.md]] | checklist akceptacyjny + smoke tests | TODO |
| [[patch_log.md]] | log każdej edycji | DOING |
| [[FINDINGS.md]] | eksportowalne wyniki | TODO |
| [[NEEDS.md]] | open issues | TODO |

## Cross-references

- [[../../audyt/CLEANUP_INVENTORY_2026-05-04.md]] — inwentaryzacja źródłowa
- [[../op-g0-r3-from-canonical-projection/README.md]] §7.4 — strategia ADDENDUM (kontekst)
- [[../audyt_cosmology_drift_2026-05-03/README.md]] — wzorzec NON-BREAKING patches
