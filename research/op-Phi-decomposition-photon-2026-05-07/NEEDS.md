---
title: "NEEDS — op-Phi-decomposition-photon (Stage 2)"
date: 2026-05-07
parent: "[[README.md]]"
type: needs-list
cycle: Stage 2
tags:
  - needs
  - Stage2
  - photon-ontology
  - open-questions
---

# NEEDS — Stage 2 photon ontology

## Phase 1 needs (formal Φ̄+δΦ decomposition)

- **N1**: Definition operacyjna `Φ̄(t)`. Czy `<Φ>_volume` po jakim
  zakresie? Cosmological FRW box (Hubble volume H_0⁻¹)? Lokalnie wokół
  obserwatora? Czy jest jednoznaczne dla każdego obserwatora?
  - Status: OPEN
  - Decyzja: prawdopodobnie cosmological avg (Hubble volume) dla globalnego
    Φ̄(t); lokalnie δΦ obejmuje wszystko poza tym

- **N2**: Czy K(Φ) musi być rozwijane w szereg czy traktować jako stałe
  K(Φ̄) na poziomie liniowym?
  - Status: OPEN
  - Decyzja: zacząć od constant K(Φ̄), sprawdzić czy K'(Φ̄)·δΦ wprowadza
    isnatne nieliniowe efekty

- **N3**: Wybór tła ψ̄. Domyślnie ψ̄ = 1 (obecna epoka, zamknięte M9.1''
  weak-field). Ale w erze radiacyjnej ψ̄ < 1 (Phase 4). Czy linearyzacja
  wokół ψ̄ < 1 nadal działa?
  - Status: OPEN
  - Decyzja: Phase 1 robić wokół ψ̄ = 1; Phase 4 weryfikować inny ψ̄

- **N4**: Czy V''(ψ̄=1) jest dodatnie czy ujemne?
  V(ψ) = (γ/3)ψ³ - (γ/4)ψ⁴
  V'(ψ) = γψ² - γψ³ = γψ²(1-ψ)
  V''(ψ) = 2γψ - 3γψ² = γψ(2-3ψ)
  V''(1) = γ·(2-3) = -γ < 0
  - Status: KRYTYCZNY — V''(ψ̄=1) = -γ < 0 oznacza **niestabilność
    tachyonowa** dookoła ψ̄=1!
  - Decyzja: weryfikacja vs closure_2026-04-26 T-Λ; być może T-Λ
    dodaje stabilizujące kontrybucje, lub β=γ vacuum condition zmienia
    efektywne V; jeśli nie — Phase 1 wykryje krytyczną
    niezgodność (tachyon → instability → STRUCTURAL_NO_GO)

## Phase 2 needs (foton jako mod δΦ)

- **N5**: Kanoniczna kwantyzacja w curved spacetime z varying c(Φ̄).
  Standard QFT zakłada flat Minkowski; Phase 2 musi obsłużyć metrykę
  M9.1'' background.
  - Status: OPEN
  - Decyzja: zacząć od locally-flat tangent space approximation, ψ̄≈1

- **N6**: Czy stress-energy δΦ ma T^μ_μ = 0?
  Skalarne pole z masą m: T^μ_μ = -m²·φ² (nie zero!)
  Skalarne pole bez masy: T^μ_μ = 0 (Weyl-niezmiennicze 4D)
  Z kosmologiczną stałą V(ψ̄): T^μ_μ ≠ 0 ogólnie
  - Status: KRYTYCZNY — Może istnieć konflikt z L01 (ρ_EM = 0)
  - Decyzja: jeśli T^μ_μ_(δΦ) ≠ 0, może wymagać re-open L01 albo
    ścieżki β (gradient) lub γ (TT-mod), które mają inne stress-energy

- **N7**: Energia próżniowa modów δΦ. Suma 1/2·ℏω_k po wszystkich k
  daje rozbieżność UV. Standard renormalizacja.
  - Status: OPEN, manageable
  - Decyzja: zero-point subtraction standardowy; cosmological constant
    z T-Λ closure jest niezależnie ustalona (β·H_0² ~ Λ_today)

## Phase 3 needs (POLARIZATION — KRYTYCZNE)

- **N8**: Wybór mechanizmu polaryzacji (α/β/γ).
  - α: longitudinal-only — FAIL obserwacyjnie (foton transverse only)
  - β: gradient ∇δΦ daje 3 DOF; longitudinal mode się łatwo eliminuje
    bo gauge condition (∂_μ A^μ = 0)? — ale to wymaga że δΦ + ∇δΦ
    razem stanowią coś jak A_μ
  - γ: TT-mod ∂_i∂_j δΦ — TT projekcja daje 2 DOF, ale propagacja
    fali wymaga że ω² = c²k², co dla ∂² operatora wprowadza ω⁴
    zależność (problematyczne)
  - Status: KRYTYCZNY OPEN
  - Decyzja: Phase 3 musi zdecydować, **jest to make-or-break dla Stage 2**

- **N9**: Czy alternatywa β (gradient ∇δΦ) wymaga wprowadzenia osobnego
  pola wektorowego A_μ? Jeśli tak, narusza S05 (single-Φ axiom).
  - Status: OPEN, S05 dependency
  - Decyzja: jeśli ∇δΦ tworzy efektywny A_μ NIE niezależny field, to
    spójne z S05; w przeciwnym razie re-open S05

- **N10**: Spin fotonu = 1. Czy mod δΦ niesie spin = 1?
  - Skalarny δΦ → spin 0 (FAIL!)
  - ∇δΦ → spin 1 (możliwe, ale longitudinal+transverse)
  - ∂²δΦ TT → spin 2 (FAIL — to grawiton)
  - Status: KRYTYCZNY OPEN
  - Decyzja: Phase 3 musi dopasować spin = 1; być może wymaga modu
    pochodnej pierwszego rzędu

## Phase 4 needs (BBN return — CONDITIONAL)

- **N11**: Jeśli ρ_rad = density of δΦ-modes, jak to wpływa na Friedmann
  eq w erze radiacyjnej?
  - Status: OPEN, Phase 4 scope
  - Decyzja: H_TGP² = (8πG(Φ̄)/3)·(ρ_matter + ρ_rad_(δΦ-modes)) +
    Φ̄-kinetic terms

- **N12**: Czy w erze radiacyjnej ψ̄ << 1 dekompozycja Φ̄+δΦ nadal działa?
  - Status: OPEN, depends on Phase 1 N3
  - Decyzja: jeśli ψ̄ << 1, może wymagać innego rozwinięcia (nie wokół
    ψ̄ = 1); Phase 4 musi to obsłużyć

- **N13**: Czy reactivacja varying-c effects ratuje BBN drift?
  EXT-1 ścieżka A: drift 99%. Stage 2 z explicit ρ_rad: TBD.
  - Status: OPEN, Phase 4 critical decision
  - Decyzja: jeśli drift nadal > 5%, EXT-1 STRUCTURAL_NO_GO utrzymany;
    Stage 2 ma własną wartość niezależnie od EXT-1

## Cross-cycle integration needs

- **N14**: Aktualizacja sek04_stale.tex jeśli Stage 2 pokazuje że
  c jest funkcją Φ̄ (tła) a NIE Φ (lokalnego)
  - Status: OPEN, post-Phase-1
  - Decyzja: edit sek04 prop:c-from-metric annotation

- **N15**: Aktualizacja sek08a_akcja_zunifikowana.tex jeśli Φ-EOM
  rozdzielić explicitly na background+perturbation
  - Status: OPEN, post-Phase-2
  - Decyzja: dodać sekcję "Linearization around background"

- **N16**: L08 kink-fermion connection — explicit unifikacja:
  - statyczne δΦ → fermion (L08)
  - oscylujące δΦ → foton (this cycle)
  - obojga to perturbations of common background Φ̄
  - Status: OPEN, post-Phase-3
  - Decyzja: cross-link Stage 2 ↔ L08 audit

## Status summary

- Phase 1 needs: N1-N4 (N4 KRYTYCZNY tachyon question)
- Phase 2 needs: N5-N7 (N6 KRYTYCZNY T^μ_μ question)
- Phase 3 needs: N8-N10 (N8, N10 KRYTYCZNE polarization+spin)
- Phase 4 needs: N11-N13 (N13 KRYTYCZNY BBN return)
- Cross-cycle: N14-N16

**Total: 16 NEEDS, 5 KRYTYCZNYCH (N4, N6, N8, N10, N13).**

**Najwiekszą strukturalną przeszkodą jest N10 (spin fotonu = 1 z modu
δΦ).** Skalarny δΦ ma spin 0; trzeba pokazać że obserwowany spin = 1
wynika z modu δΦ, NIE z osobnego pola A_μ.

To jest największe ryzyko dla Stage 2 STRUCTURAL_NO_GO i wymaga
najwcześniejszej weryfikacji (przed pełną Phase 1 może warto
quick-check Phase 3 alternatyw).
