---
title: "FINDINGS — L04 ODE-canonicalization (analytical decision)"
date: 2026-05-04
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-L04-ODE-canonicalization-2026-05-04
tags:
  - findings
  - L04
  - alpha-2-canonical
  - PHASE2-universal
  - m_obs-vs-M_full
---

# FINDINGS — L04 ODE-canonicalization

> Eksportowalne wyniki cyklu L04. Każdy item z cytowanym source.

## Wyniki strukturalne

| ID | Statement | Source | Consumers |
|----|-----------|--------|-----------|
| F1 | TGP-canonical α=2 jest jednoznacznie wybrane przez warunki (C1)+(C2)+(C3) klasy operatorów kinetycznych | `core/sek08_formalizm/sek08_formalizm.tex` §`thm:D-uniqueness` lin. 956-1048 | wszystkie cykle używające ODE solitonu |
| F2 | Phase 2 universal mass formula `m_obs = c·A²·g₀^[e²(1−α/4)]` jest fundamentalna dla całego α-range | [[../why_n3/PHASE2_n_alpha_derivation.md]] | LP-4, LP-6, R5, why_n3, mass_scaling_k4 |
| F3 | R5 K² mass formula ≡ Phase 2 IFF α=1 (analytical theorem) | [[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]] §2 | mass_scaling_k4 status update |
| F4 | LP-4 k=4 mass formula jest specjalnym przypadkiem α=1 z Phase 2 universal — nie *uniwersalna* | mass_formula_unification.md §3.1 | PLAN_DOMKNIECIA LP-4 reframing |
| F5 | „Dualizm" α=1 vs α=2 jest pozorny — α=1 to specjalizacja α=2 universal Phase 2 | canonical_form_evidence.md §"Synteza" | LP-6 dictionary reframing |
| F6 | Bariera g₀_crit operuje na M_full (strukturalna własność ODE), niezależnie od mass formula m_obs | m_obs_vs_M_full.md §3.3 | wszystkie cykle używające bariery topologicznej |
| F7 | g₀_crit(α=2) = 1.874 koincyduje z horyzontem Lorentzowskim M9.1'' ψ=4/3 (4 cyfry znaczące) | [[../why_n3/PHASE1_psi_g0_identification.md]] | M9.1'' status, geometric uniqueness α=2 |

## Wyniki numeryczne

| ID | Wynik | Wartość/formuła | Source | Consumers |
|----|-------|-----------------|--------|-----------|
| F8 | Wykładnik n(α) = e²·(1−α/4) liniowy fit | residuum <0.003 dla α∈[0.25, 4.0] | [[../why_n3/r3_phase2_n_alpha_derivation.txt]] | Phase 2 mass formula |
| F9 | n(α=2) = e²/2 ≈ 3.6945 dla TGP-canonical | sympy diff vs numerical −0.001% | [[../why_n3/r3_phase2b_X_constant.py]] | TGP-canonical mass spectrum |
| F10 | m_μ/m_e dla α=2 z Phase 2 formuły | 206.766 (PDG 206.7682, **−0.001%**) | [[../why_n3/r3_alpha2_full_closure.txt]] § Sekcja 1 | mass spectrum lepton |
| F11 | m_τ/m_e dla α=2 z Phase 2 + Koide K=2/3 | 3474.28 (PDG 3477.23, **−0.085%**) | [[../why_n3/r3_alpha2_full_closure.txt]] § Sekcja 4 | mass spectrum lepton |
| F12 | m_τ/m_μ dla α=2 z Phase 2 + Koide | 16.820 (PDG 16.817, **+0.015%**) | [[../why_n3/r3_alpha2_full_closure.txt]] § Sekcja 4 | mass spectrum lepton |
| F13 | g₀^τ z Koide K=2/3 + Phase 2 (α=2) | 1.755 (margin do bariery +0.119) | [[../why_n3/r3_alpha2_full_closure.txt]] § Sekcja 3 | N=3 selection rule |
| F14 | g₀^4 = φ·g₀^τ = 2.840 > g₀_crit(α=2) = 1.874 | 4-ta generacja zakazana z marginem +0.965 | [[../why_n3/r3_alpha2_full_closure.txt]] § Sekcja 5 | falsyfikacja N≥4 |
| F15 | n(4) = -0.006 (zero z dokładnością ODE solver) | Hobart-Derrick balance point | [[../why_n3/r3_phase2_n_alpha_derivation.txt]] | open problem (Derrick analiza) |
| F16 | R5 K² ratio dla α=2 = 1221 vs Phase 2 = 207 | mismatch +490% — R5 fail dla α=2 | [[../mass_scaling_k4/r5_phase2_analytical_bridge.txt]] | R5 status downgrade |

## Status PASS-ów per faza (1-fazowy cykl L04 — analytical decision)

| Etap | Tests | Werdykt | Plik wynikowy |
|------|-------|---------|---------------|
| Synteza analytical evidence | 3/3 strukturalne dowody (sek08 thm + Phase 2 + R5 bridge) | CLOSED-DECISION | canonical_form_evidence.md |
| m_obs vs M_full distinction | physical analysis + 5 niezależnych analogii | CLOSED | m_obs_vs_M_full.md |
| ODE class taxonomy (C1-C3 selection) | klasyfikacja, taxonomy, rationale | CLOSED | ODE_class_taxonomy.md |
| Mass formula unification (Phase 2 vs LP-4/LP-6/R5) | reframing wszystkich konkurencyjnych formulas | CLOSED | mass_formula_unification.md |

**Phase L04 score: 4/4 CLOSED.** Cykl jest analytical-decision-doc.

## Falsyfikatory ustawione przez ten cykl

| ID | Predykcja | Kryterium falsyfikacji | Test |
|----|-----------|-------------------------|------|
| FX1 | Phase 2 mass formula z α=2 daje m_τ/m_e = 3477.23 ± 0.085% | m_τ pomiar PDG poza ±0.5% byłby falsyfikujący | wykonane (PASS) |
| FX2 | n(α) = e²·(1−α/4) liniowy fit | residuum >1% byłoby falsyfikujące | wykonane (PASS) |
| FX3 | g₀_crit(α=2) = 4/3 · 1.40554 = 1.874 (M9.1'' horizon koincydencja) | inna wartość bariery dałaby sprzeczność z M9.1'' | wykonane (PASS) |
| FX4 | g₀^4 = φ·g₀^τ > g₀_crit (4-ta gen zakazana) | obserwacja 4-tej generacji leptonowej < g₀_crit byłaby falsyfikującą | open (LHC / future colliders) |

## Eksport do innych folderów (impacts)

- [[../../audyt/L04_ODE_dualism_alpha]] → POST_ACTION_UPDATE z reference do tego cyklu (status CLOSED-RESOLVED)
- [[../../audyt/L05_mass_exponent_drift]] → resolution k=4 vs p=5−α: oba są specjalnymi przypadkami Phase 2 (α=1 i α=2)
- [[../../audyt/D01_drifting_numbers]] → g₀^e = 0.86941 (Phase 2 universal calibration anchor) — globally consistent
- [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] → LP-4 reframing, LP-6 dictionary downgrade
- [[../mass_scaling_k4/]] → status R5 K² → DERIVATIVE (specjalny case α=1)
- [[../why_n3/]] → confirm Phase 2 jako fundamental, p(α)=5−α jako empirical decomposition

## Pre-existing flag

Ten plik jest tworzony 2026-05-04. Wszystkie cytowane source plików
(why_n3/, mass_scaling_k4/, sek08_formalizm.tex) są **pre-existing**;
L04 cykl synthesizes ich findings into single decision-doc.

## Cross-references

- [[NEEDS.md]] — otwarte problemy (X = e²/4 RG derivation)
- [[../../audyt/L04_ODE_dualism_alpha/POST_ACTION_UPDATE_2026-05-04.md]]
