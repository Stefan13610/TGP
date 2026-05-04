---
title: "L04 ODE-canonicalization — kanoniczna formulacja TGP α=2 (analytical decision)"
date: 2026-05-04
cycle: L04
type: audit-resolution
status: ANALYTICAL DECISION DOCUMENTED — α=2 (TGP-canonical) JEST jedyną kanoniczną formulacją
parent: "[[../../audyt/L04_ODE_dualism_alpha/README.md]]"
predecessors:
  - "[[../why_n3/PHASE2_n_alpha_derivation.md]]"
  - "[[../why_n3/CORRECTIONS_2026-05-01.md]]"
  - "[[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]]"
  - "[[../op-g0-r3-from-canonical-projection/README.md]]"
related:
  - "[[../../core/sek08_formalizm/sek08_formalizm.tex]]" (thm:D-uniqueness lin. 956–1048)
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]]"
  - "[[../../meta/PLAN_DOMKNIECIA_MASTER.md]]" (LP-4, LP-6, N-1)
tags:
  - TGP
  - L04
  - ODE-canonicalization
  - alpha-selection
  - m_obs-vs-M_full
  - mass-formula-universal
  - PHASE2-e-squared
  - audit-resolution
tgp_status:
  folder_status: audit
  level: L1
  kind: derivation
  core_compatibility: current
  last_reviewed_against_core: 2026-05-04
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: ["X = e²/4 RG derivation"]
  depends_on:
    - "[[../why_n3/PHASE2_n_alpha_derivation.md]]"
    - "[[../op-g0-r3-from-canonical-projection]]"
  impacts:
    - "[[../../audyt/L04_ODE_dualism_alpha]]"
    - "[[../../audyt/L05_mass_exponent_drift]]"
    - "[[../../audyt/D01_drifting_numbers]]"
  source_of_status:
    - "thm:D-uniqueness (sek08_formalizm.tex:958)"
    - "PHASE2 e² discovery (why_n3/PHASE2_n_alpha_derivation.md)"
    - "R5 PHASE2 bridge (mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md)"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L04 — kanoniczna formulacja TGP: **α=2 jest WYBRANE STRUKTURALNIE**

## Cel

Rozstrzygnąć rzekomy „dualizm" formulacji ODE (α=1 K=g² „substratowa"
vs α=2 K=g⁴ „kanoniczna") wskazany w [[../../audyt/L04_ODE_dualism_alpha]],
**wykorzystując distinction `m_obs` vs `M_full`** (insight użytkownika
2026-05-01) jako klucz interpretacyjny.

## Werdykt strukturalny

**TGP-canonical α=2 (K=K_geo·φ⁴) JEST jedyną poprawną kanoniczną
formulacją** z trzech niezależnych powodów:

1. **Twierdzenie strukturalne** (`thm:D-uniqueness`,
   `core/sek08_formalizm.tex:958`): α=2 jest *jednoznacznie wybrane*
   przez 3 warunki klasy operatorów (C1) stałość α + (C2) K(0)=0 +
   (C3) geometric substrate coupling z `prop:substrate-action`.
2. **Uniwersalna mass formula** (PHASE2 closure 2026-05-01):
   `m_obs = c · A² · g₀^[e²(1−α/4)]` daje <0.1% PDG dla α=2
   ([[../why_n3/PHASE2_n_alpha_derivation.md]]).
3. **R5 = Phase 2 IFF α=1** (analytical theorem,
   [[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]]):
   formuła `m=c·A⁴` jest **specjalnym przypadkiem** Phase 2 dla α=1, NIE
   uniwersalnym mechanizmem.

**„Dualizm" jest pozorny:** α=1 nie jest *równoważną* formulacją z α=2,
tylko *specjalnym przypadkiem* uniwersalnej formuły z α-zależnym
wykładnikiem.

## Centralna distinction: `m_obs` vs `M_full`

[[m_obs_vs_M_full.md]] — pełna analiza fizyczna.

**Krótko:**

| Wielkość | Definicja | Co charakteryzuje |
|----------|-----------|-------------------|
| `M_full` | Pełna energia struktury wewnętrznej solitonu: K + V_eff | Strukturalna własność ODE — **bariera g₀_crit operuje na M_full** |
| `m_obs` | „Waga z dystansu" — masa mierzona przez asymptotyczną projekcję A_tail | Obserwowalna w eksperymencie — formula zależna od α |

Analogia (z [[../why_n3/CORRECTIONS_2026-05-01.md]]):

- **GR**: masa ADM (asymptotyczna) ≠ masa Komara (wewnętrzna)
- **QFT**: bare mass (UV) ≠ renormalized mass (IR)
- **EM**: ładunek (asymptotyczny) ≠ energia pola (objętościowa)

W TGP:
```
m_obs(g₀, α) = c_M · A_tail²(g₀, α) · g₀^[e²(1−α/4)]   ← Phase 2 universal
M_full(g₀, α) = K + V_eff (z ODE energy integral)      ← strukturalne
```

## Pliki w cyklu

| Plik | Opis |
|------|------|
| [[README.md]] | (ten plik) — werdykt + indeks |
| [[m_obs_vs_M_full.md]] | Pełna analiza fizyczna distinction |
| [[canonical_form_evidence.md]] | 3 niezależne dowody α=2 strukturalnej |
| [[ODE_class_taxonomy.md]] | Taxonomia klas operatorów C1-C3 i selekcja α |
| [[mass_formula_unification.md]] | Phase 2 universal vs LP-4/LP-6 specjalne przypadki |
| [[FINDINGS.md]] | Eksportowalne wyniki (anti-overclaim) |
| [[NEEDS.md]] | Otwarte luki (głównie X=e²/4 RG derivation) |

## Co zostało zrobione w tej sesji (2026-05-04)

1. **Synteza przeprowadzona** — pięć niezależnych analiz pokazuje
   spójność TGP-canonical α=2:
   - thm:D-uniqueness (sek08_formalizm)
   - Phase 2 e² discovery (why_n3 2026-05-01)
   - R5 ↔ Phase 2 analytical bridge (mass_scaling_k4 2026-05-02)
   - G.0 closure (op-g0-r3-from-canonical-projection 2026-05-02)
   - r3_alpha2_full_closure 6/6 PASS

2. **Decyzja autorska sformulowana**: α=2 jest kanoniczna; α=1 to
   specjalny case LP-4 (k=4 = 2(d-1)/(d-2) jest argument *konwergencyjny
   tylko dla K=φ²*).

3. **Open problem zidentyfikowany**: derywacja X = e²/4 z RG flow lub
   Hobart-Derrick balance pozostaje OPEN (Phase 6 Q5 R⁵-bridge NEGATIVE).

4. **Brak modyfikacji core**: cykl jest **analytical-decision-doc**,
   nie modyfikuje sek08*.

## Status w audycie L04

Patrz [[../../audyt/L04_ODE_dualism_alpha/POST_ACTION_UPDATE_2026-05-04.md]].

## Cross-references

### Strukturalne (theorem-level)
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`thm:D-uniqueness` (lin. 956–1048)
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`rem:alpha2-pivot-status-pl` (lin. 1050–1068)
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] §`thm:D-uniqueness` ref (lin. 39, 78, 817)

### Numerical evidence
- [[../why_n3/r3_alpha2_full_closure.txt]] (6/6 PASS dla α=2 + p=3)
- [[../why_n3/r3_observable_vs_full_mass.txt]] (p(α) = 5−α empirical scan)
- [[../why_n3/r3_p_alpha_analytical.txt]] (PHASE2 anchor + n(α) linear fit)
- [[../mass_scaling_k4/r5_phase2_analytical_bridge.txt]] (R5 ≡ Phase 2 IFF α=1)

### Audit context
- [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-4, LP-6, N-1
- [[../../meta/AUDYT_TGP_2026-05-01.md]] (post-audit context)
- [[../../audyt/L04_ODE_dualism_alpha]] (audit luki)
- [[../../audyt/L05_mass_exponent_drift]] (powiązany — k=4 vs p=5−α)
