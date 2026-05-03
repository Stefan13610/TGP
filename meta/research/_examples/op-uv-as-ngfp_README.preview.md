<!--
PODGLĄD (Sesja 2). Pokazuje, jak README.md folderu `research/op-uv-as-ngfp/`
WYGLĄDAŁBY po zaaplikowaniu szablonu z meta/research/templates/.

NIE ZASTĘPUJE oryginalnego README. Oryginał jest w
`research/op-uv-as-ngfp/program.md` (ten folder nie ma jeszcze README.md).

Cel podglądu: zwalidować, czy schemat YAML pasuje do realnego folderu
3-fazowego op-* z Phase{1,2,3}_{setup,results}.md.

Źródła użyte do wypełnienia (anty-overclaim — wszystko cytowane):
  - research/op-uv-as-ngfp/program.md
  - research/op-uv-as-ngfp/Phase1_results.md
  - research/op-uv-as-ngfp/Phase2_results.md
  - research/op-uv-as-ngfp/Phase3_results.md
  - meta/research/AUDIT_RESEARCH_S1.md (link_total=9, PASS=164,
    LOCKED=77, mtime 2026-04-29)
-->

---
title: "UV.1 — AS NGFP first-principles closure for N_A normalization"
date: 2026-04-29
parent: "[[INDEX.md]]"
related:
  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"
  - "[[meta/research/RESEARCH_BUS.md]]"
  - "[[../op-xi-photon-ring/Phase3_results.md]]"
  - "[[../op-phase3-uv-completion/Phase3_R_final_results.md]]"
  - "[[../op-cross-sector-charge/Phase3_results.md]]"
tags:
  - TGP
  - uv-completion
  - asymptotic-safety
  - NGFP
  - Reuter
tgp_status:
  folder_status: active
  level: L3
  kind: derivation
  core_compatibility: current
  last_reviewed_against_core: 2026-05-03
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges:
    - "UV1 2-loop FRG residual closure (current 0.068% vs target 0.011%)"
  depends_on:
    - "research/op-xi-photon-ring"             # ξ.1 N_A target 500/57
    - "research/op-phase3-uv-completion"       # AS NGFP KEYSTONE 12/12
    - "research/op-cross-sector-charge"        # XS.1 cross-sector
    - "research/closure_2026-04-26"            # F4 rational + α₀ chain
  impacts:
    - "research/op-bh-alpha-threshold"         # AS RG signature LIGO O5
    - "research/op-eps-photon-ring"            # ngEHT 2030+ photon-ring
    - "research/op-uv-renormalizability-research"  # UV2-UV7 long-term
  source_of_status:
    - "Phase1_results.md (5/5 PASS, 2026-04-29)"
    - "Phase2_results.md (PARTIALLY DERIVED, 2-loop residual 0.068%, 2026-04-29)"
    - "Phase3_results.md (6/6 PASS, UV.1 program END, 2026-04-29)"
    - "program.md §Predecessors (cross-link ξ.1.Phase3 N_A target)"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-03
---

# UV.1 — AS NGFP first-principles closure for N_A normalization

## Cel

Domknąć single residual N_A = 8.7719 algebraic provenance pozostawiony przez
ξ.1.Phase3 z asymptotic-safety NGFP first-principles (Reuter 1998 / Eichhorn
2018). Próba pochodzenia N_A bezpośrednio z Litim invariant g*·λ* = 0.1349 +
η_N* = -2 + heat-kernel a₂ pod NGFP RG flow.

## Stan

**UV.1 program END** (2026-04-29) — Phase 3 zamknięte 6/6 PASS,
status cascade ACTIVATED:
- ξ.1: PARTIALLY DERIVED (refined) → **DERIVED (refined²)**
- XS.1: PARTIALLY DERIVED (refined) → **DERIVED (refined²)**
- UV7: STRUCTURAL-DERIVED → **DERIVED**

Phase 2 zakończona z residualem 2-loop 0.068% (target 0.011%) — pozostawia
single open bridge (zob. [[NEEDS.md]] N1).

Source: `Phase3_results.md` lin. 21–45.

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| `program.md` | UV.1 3-phase mini-cycle (5+7+6 sub-tests) | CLOSED |
| `Phase1_setup.md` / `Phase1_results.md` | AS NGFP foundational audit | CLOSED 5/5 |
| `Phase2_setup.md` / `Phase2_results.md` | N_A first-principles derivation | PARTIALLY DERIVED |
| `Phase3_setup.md` / `Phase3_results.md` | UV1–UV6 predictions + status cascade | CLOSED 6/6 |
| `phase1_ngfp_audit.py` + `.txt` | AS NGFP foundational sanity | OK |
| `phase2_NA_derivation.py` + `.txt` | N_A z g*, λ*, η_N* | OK |
| `phase3_predictions_status.py` + `.txt` | predictions UV1–UV6 + cascade | OK |

## Wejścia (depends_on)

- [[../op-xi-photon-ring/Phase3_results.md]] — ξ.1 N_A target 500/57 = 8.7719
- [[../op-phase3-uv-completion/Phase3_R_final_results.md]] — AS NGFP synthesis 4-of-4
- [[../op-cross-sector-charge/Phase3_results.md]] — XS.1 PARTIALLY DERIVED baseline
- [[../closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]] — F4 rational, α₀ chain

## Wyjścia (impacts)

- [[../op-bh-alpha-threshold]] — AS RG signature for LIGO O5 dispersion
- [[../op-eps-photon-ring]] — ngEHT 2030+ photon-ring sharpening band
- [[../op-uv-renormalizability-research]] — UV2–UV7 long-term track (UV1 closed)

## Otwarte luki

Pełna lista w [[NEEDS.md]]. Headlines:

- **N1 (bridge)**: 2-loop FRG residual closure (current 0.068% > 2-loop band 0.011%) —
  Sesja 8 regression target.

## Cross-references

- [[meta/PLAN_RESEARCH_WORKFLOW_v1.md]] — workflow
- [[meta/research/RESEARCH_BUS.md]] — broadcast UV.1 cascade
- [[meta/AUDYT_TGP_2026-05-01.md]] — audyt całości
