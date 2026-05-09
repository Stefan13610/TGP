---
title: "L01 ρ stress-energy bridge — formal kowariantna definicja gęstości materii"
date: 2026-05-04
cycle: L01
type: audit-resolution
status: ANALYTICAL DECISION DOCUMENTED — ρ ≡ -T^μ_μ/c_0² jako formal definition + mapping SM
parent: "[[../../audyt/L01_rho_operational/README.md]]"
predecessors:
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]" (§A.4, §C.3, §N — Option-2)
  - "[[../op-newton-momentum/B9_wep_microscope_composition_results.md]]" (B9 closure 2026-05-01)
related:
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]" (lin. 102-129, eq:L-mat-unified)
  - "[[../../core/sek08_formalizm/sek08_formalizm.tex]]" (lin. 11257-11337, ax:metric-coupling)
  - "[[../../TGP_FOUNDATIONS.md]]" §4 (warstwa 3a/3b/3c)
tags:
  - TGP
  - L01
  - rho-formal-definition
  - stress-energy
  - matter-coupling
  - SM-sector-mapping
  - photon-treatment
  - audit-resolution
tgp_status:
  folder_status: closed-resolved
  level: L1
  kind: derivation
  core_compatibility: current
  last_reviewed_against_core: 2026-05-04
  may_edit_core: true
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: ["quantum-trace-anomaly"]
  depends_on: []
  impacts:
    - "[[../../audyt/L01_rho_operational]]"
    - "[[../../audyt/S04_metric_coupling_axiom]]"
  source_of_status:
    - "ax:metric-coupling (sek08_formalizm.tex:11257)"
    - "eq:L-mat-unified (sek08a:103-111)"
    - "AUDYT 2026-05-01 §A.4 Option-2 decision"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L01 — formal kowariantna definicja `ρ`

## Cel

Operacyjnie domknąć status `ρ` w `L_mat = -(q/Φ_0)·φ·ρ` (sek08a:103)
poprzez:

1. Explicit kowariantna definicja `ρ ≡ -T^μ_μ/c_0²` w sek08a
2. Mapping na sektory SM: Dirac (m·ψ̄ψ), scalar (gradient terms), EM (T^μ_μ=0)
3. Treatment fotonów (`ρ_EM = 0` — co to znaczy fizycznie)
4. NON-BREAKING addytywna edycja sek08a

## Werdykt

**`ρ ≡ -T^μ_μ/c_0²`** jest formal kowariantną definicją gęstości
materii w `L_mat` sek08a. Dla:

- **Dirac fermion** (m·ψ̄ψ): `ρ_Dirac = m·ψ̄ψ/c_0²` ≥ 0 ✓ → standardowa
  gęstość masy bezwładnościowej
- **Massless EM field**: `T^μ_μ_EM = 0` → `ρ_EM = 0` → **fotony nie
  generują pola Φ przez `L_mat`**, ale propagują się geodezyjnie po
  `g_eff` (`prop:coupling-consequences` 2)
- **Massive scalar**: `T^μ_μ_scalar` zawiera `m²|φ|²` term + gradient
  contributions
- **Yang-Mills gauge**: `T^μ_μ_YM = 0` (klasycznie); quantum trace anomaly
  modyfikuje to przez β-funkcję

Ta definicja **ratuje strukturalną spójność** `ax:metric-coupling`
(audit Option-2): `L_mat = -(q/Φ_0)·φ·ρ` jest **derived consequence**
volume element `√(-g_eff) ∝ φ` plus stress-energy tensor materii w
kanonicznym sprzęganiu metrycznym.

## Pliki w cyklu

| Plik | Opis |
|------|------|
| [[README.md]] | (ten plik) |
| [[formal_definition.md]] | explicit `ρ ≡ -T^μ_μ/c_0²` derivation z `L_mat[ψ_m, g_eff]` |
| [[SM_sector_mapping.md]] | Dirac/scalar/EM/gauge — explicit ρ dla każdego sektora |
| [[photon_treatment.md]] | T^μ_μ_EM=0 i implikacje (GW170817 c_GW=c, no fifth force from photons) |
| [[FINDINGS.md]] | eksportowalne wyniki |
| [[NEEDS.md]] | open problems (quantum trace anomaly) |

## Co wykonano w sesji

1. **Synthesis 5 plików rdzenia + audit** — sek08a §74-129, sek08_formalizm
   §11257-11337, AUDYT 2026-05-01 §A.4+§C.3+§N, B9 closure 2026-05-01.
2. **Formal derivation** `ρ = -T^μ_μ/c_0²` z `L_mat[ψ_m, g_eff]` —
   demonstracja, że Option-2 audytu jest *wyprowadzeniem*, nie *postulatem*.
3. **Mapping SM sector-by-sector** — Dirac, scalar, EM, Yang-Mills.
4. **Treatment fotonów** — explicit `ρ_EM = 0` z T^μ_μ_EM = 0; brak
   piątej siły *od fotonów* (konsystentne z B9 MICROSCOPE 6/6 PASS).
5. **NON-BREAKING addytywne edits** w sek08a (planowane, opcjonalne).

## Status w audicie L01

Patrz [[../../audyt/L01_rho_operational/POST_ACTION_UPDATE_2026-05-04.md]].

## Open problems

- **N1: Quantum trace anomaly** — `T^μ_μ ≠ 0` w QFT na curved background
  pochodzi od *trace anomaly* (Birrell-Davies 1982); modyfikacja
  ρ przez `O(R²)` corrections wymaga dedicated cycle dla **renormalized
  ρ** w obecności `g_eff[Φ]`.
- **N2: SPARC fits** — N-body symulacje używają `ρ = ρ_baryon` (HI + stars
  + bulge); explicit weryfikacja że to jest non-relativistic limit
  `T^μ_μ = -ρ_rest·c²` jest pożądana (cosmetic, low priority).

Patrz [[NEEDS.md]] dla pełnej listy.

## Cross-references

- [[../../audyt/L01_rho_operational]]
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] §`eq:L-mat-unified`
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`ax:metric-coupling`
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.4, §C.3, §N
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]]
