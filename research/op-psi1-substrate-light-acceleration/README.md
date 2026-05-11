---
title: "op-psi1-substrate-light-acceleration"
date: 2026-05-03
parent: "[[../INDEX.md]]"
related:
  - "[[meta/research/FOLDER_STATUS_INDEX.md]]"
tags:
  - TGP
tgp_status:
  folder_status: paused
  level: L2
  kind: derivation
  core_compatibility: "unknown"
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "Phase{1..7}_results.md PASS=169, CLOSED=23, FAIL=38"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

# op-psi1-substrate-light-acceleration

> **Sesja 4 auto-generated README** (2026-05-03).
> Folder nie miał wcześniej README. Treść poniżej jest minimalna —
> uzupełnij ręcznie zgodnie z [[meta/research/templates/README.template.md]].

## Cel

Substrate-light photon-substrate coupling: derivation EFT-grade dim-6 operatorów
sprzęgających pole substratu X (≡ Φ/Φ_0) z F_μν gauge field photon. Cykl odpowiada
na pytanie: *"czy istnieje native TGP mechanism dla anisotropic c_eff(θ) lub
vacuum birefringence z gradientów substratu?"* — TAK (v3 PASS-CLOSED, 2-element
canonical basis).

## Stan

| Faza | Cycle version | Status | Główny rezultat |
|---|---|---|---|
| Phase 1 | v1 | **INVALIDATED** | Z(x)·F² scalar prefactor — gauge-equivalent std EM (varying-α, NIE varying-c); A6 audit invalidation 2026-05-01 |
| Phase 2 | v1 | INVALIDATED | (consequence of v1 framework) |
| Phase 3 | v1 | INVALIDATED | (consequence of v1 framework) |
| Phase 4 | v2 | PASS | Tensor L₅'_a `(∂lnX)(∂lnX)·F·F` canonical — strukturalnie różny od v1, modyfikuje light cone genuinely |
| Phase 5 | v2 | PASS | Eikonal dispersion LOCK; anisotropic c_eff(θ) sympy-LOCK |
| Phase 6 | v2 | PASS | TT19-TT23 5/5 PASS; β_g sign FORCED < 0 (Adams DECISIVE Phase 6.T6.5) |
| Phase 7 | v3 | **PASS — CLOSED** | Hilbert-series-style enumeration → 2-element canonical basis {L₅'_a parity-even, L₅'_b parity-odd}; C8 audit closure |

🟡 **Cycle status: paused** (folder_status). v3 closed (Phase 7); pending Phase 8
(parity-odd Adams analysis dla β̃_g sign — non-blocking follow-up).

## Status (2026-05-10 update)

> **ADDENDUM 2026-05-10 (interpretive overlay + L01-Q1 resolution + form-meaning lens):**
> [[ADDENDUM_2026-05-10_native_observables_first.md]] — aplikuje
> [[../../meta/PPN_AS_PROJECTION.md]] methodology (binding 2026-05-10+) plus form-meaning
> framework z [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §4.
>
> Kluczowe statementy:
>
> 1. **v1 invalidation jako form-meaning case study:** v1 forma `Z(x)·F²` była
>    BD-form (scalar prefactor) z FALSE TGP-meaning (zinterpretowana jako varying-c, gdy
>    strukturalnie jest gauge-equivalent varying-α). v2 zastąpiła tensorowym L₅'_a —
>    *strukturalna re-derivation*, nie reparametryzacja.
> 2. **Three-layer specification dla TT19-TT23 + L₅'_b vacuum birefringence:** L1 native
>    (anisotropic c_eff lub L/R photon speed difference) / L2 projection (Sagnac, TOF,
>    CMB θ rotation, GRB pol-rotation) / L3 falsifikator (TT19 NULL lab; CMB θ pending C6).
> 3. **L01-Q1 RESOLVED:** ψ.1.v3 dim-6 EFT operator class jest *disjoint* od quantum
>    trace anomaly EM (`α·F²` pure-photon, NO ∂lnX leg). L01 N1 zrealizowany 2026-05-11
>    [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] (STRUCTURAL_DERIVED) —
>    Theorem 2.1 (Disjointness) potwierdzona *konstruktywnie*.
> 4. **Native parameter audit:** ~3 indep (β_g sign-LOCKED < 0 z Adams; β̃_g pending
>    Phase 8; Λ free). Sterile sectors (varying-α, axion-like) explicit excluded.
> 5. **Cross-cycle convergence (post-N4 update 2026-05-11):** ψ.1 ADDENDUM §3 +
>    τ.3 ADDENDUM §2 + L01 ADDENDUM §3.2 + Q2 cycle + **op-L01-N1 + N2 + N3 + N4
>    cycles (2026-05-11 *konstruktywne* derivations)** zbieżne na **separable
>    sector structure** (**8-fold diagnostic**, 4× SM sektory × 2 diagnostic
>    methods). Theorem 2.1 (Disjointness) z N1 inherits to N2 (QCD), N3 (SPARC),
>    N4 (Higgs) sektory — patrz ψ.1 ADDENDUM §3.5.

## Pliki

| Plik | Cycle version | Opis |
|---|---|---|
| [[Phase1_results.md]] | v1 | INVALIDATED — Z(x)·F² scalar (varying-α, nie varying-c) |
| [[Phase2_results.md]] | v1 | INVALIDATED |
| [[Phase3_results.md]] | v1 | INVALIDATED |
| [[Phase4_results.md]] | v2 | L₅'_a tensor canonical |
| [[Phase5_results.md]] | v2 | Eikonal + anisotropic c_eff(θ) LOCK |
| [[Phase6_results.md]] | v2 | TT19-TT23 + β_g sign Adams DECISIVE |
| [[Phase7_results.md]] | v3 | Hilbert-series-style enumeration → 2-elem basis; C8 closure |
| [[FINDINGS.md]] | — | eksportowalne wyniki |
| [[NEEDS.md]] | — | otwarte luki |
| [[ADDENDUM_2026-05-10_native_observables_first.md]] | — | **interpretive overlay (2026-05-10):** native-first methodology + L01-Q1 resolution + form-meaning lens |

## Cross-references

- [[ADDENDUM_2026-05-10_native_observables_first.md]] — native-first overlay (2026-05-10)
- [[../../meta/PPN_AS_PROJECTION.md]] — parent methodology (binding 2026-05-10+)
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] — form-meaning framework (v1 case study)
- [[../op-L01-rho-stress-energy-bridge-2026-05-04]] — Q1 source (resolved by §3 tego addendum)
- [[../op-tau3-substrate-clock-acceleration]] — siostrzany cykl, downstream consumer photon-substrate sector
- [[../op-omega1-substrate-em-coupling]] — siostrzany sector (ω.1 axion-like, sterile w ψ.1 v3 enumeration)
- [[meta/research/RESEARCH_BUS.md]] — broadcast wyników
