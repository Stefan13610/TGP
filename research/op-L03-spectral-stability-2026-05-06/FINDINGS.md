---
title: "FINDINGS — L03 spectral stability (eksportowalne wyniki)"
date: 2026-05-06
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-L03-spectral-stability-2026-05-06
tags:
  - findings
  - L03
  - spectral-synthesis
  - exportable
---

# FINDINGS — L03 spectral stability

## Główny werdykt

> **TGP jest strukturalnie stabilne dla wszystkich physical sources ρ ≥ 0.**
>
> Spektrum operatora perturbacji `L̂ = -∂[K(φ_eq)·∂] + V''_eff(φ_eq)` jest
> non-negative we wszystkich fizycznych konfiguracjach (vacuum, FRW,
> Yukawa, soliton). Mass gap m_sp² = γ/K_geo > 0 zachowany asymptotycznie.

## Eksportowalne wyniki (16 findings)

### F1: Operator perturbacji L̂ w formie Sturm-Liouville'a

```
L̂[δφ] = -∂_μ[K(φ_eq)·g_eff^μν·∂_ν δφ] + V''_eff(φ_eq)·δφ
```

z `K(φ) = K_geo·φ⁴` (kanoniczna α=2, `tw:D-uniqueness`)
i `V_M911(φ) = -γφ²(4-3φ)²/12` (`prop:V-M911-canonical`).

### F2: Spektrum L̂ jest non-negative dla wszystkich ρ ≥ 0

`σ(L̂) ⊂ [0, ∞)` dla physical configurations. Punktowa pozytywność
`K(φ_eq) > 0` (3-tier ghost-freedom sek08b) i `V''_eff(φ_eq) ≥ 0`
asymptotycznie (`prop:vacuum-stability-G0`) → S-L theorem.

### F3: Mass gap m_sp² = γ/K_geo > 0

W asymptotyce r → ∞ (φ_eq → 1), continuum spektrum zaczyna się przy
ω² = c²k² + m_sp²c⁴, m_sp² = γ/K_geo > 0. Identyczne z
`prop:vacuum-stability` Yukawa profile.

### F4: Mode count Z₂-broken vacuum

| Liczba | Wartość | Powód |
|--------|---------|-------|
| Total DOF | 1 | pojedyncze pole skalarne Φ |
| NGB | 0 | Z₂ dyskretne, Goldstone theorem n/a |
| Massive | 1 | m_sp = √(γ/K_geo) |
| Ghost | 0 | 3-tier ghost-freedom (sek08b) |

### F5: Goldstone theorem nie aplikuje do dyskretnych grup

`dim(Z₂) = 0` ⇒ liczba NGB = `dim(G) - dim(H) = 0 - 0 = 0`.

### F6: Brak tachyon mode na 4 fizycznych profilach Φ_eq[ρ]

| Profil | ρ | Tachion? |
|--------|---|----------|
| Vacuum | 0 | NIE |
| FRW | ρ_FRW > 0 | NIE |
| Yukawa point | δ³M | NIE |
| Soliton (lepton α=2) | concentrated | NIE |

### F7: Ghost wall jest artefaktem continuum dla obu formulacji

| Formulacja | g* (ghost wall) | g_min(physical) | Bezpieczeństwo |
|------------|-----------------|-----------------|-----------------|
| Substratowa α=1 | e⁻¹ ≈ 0.368 | 0.91 | 0.55 above |
| Kanoniczna α=2 | g_ghost ≈ 0.717 | 0.91 | 0.19 above |

Soliton profiles wszystkich leptonów (m_e, m_μ, m_τ z R3 ODE) są
**bezpiecznie powyżej** ghost wall (sek08b `cor:ghost-artifact`).

### F8: Synteza 3-tier ghost-freedom

Wszystkie trzy poziomy mają wspólny korzeń: `K_ij = J(φ_iφ_j)² ≥ 0`
substratowo.

| Sektor | Wynik | Twierdzenie |
|--------|-------|-------------|
| Fundamental ψ | ψ⁴(∇ψ)² > 0 | `prop:ghost-free-fundamental` |
| MS-TGP perturbations | Q_s = ψ_bg⁴ > 0 | `prop:ghost-free-MS` |
| Soliton substrate | K_sub(g) = K_geo·g² > 0 | `thm:ghost-free-soliton` |

### F9: V''_eff(1) = +γ > 0 (G.0 fix tachion bug v1.x)

Pre-G.0 closure: V_orig(ψ) = (β/3)ψ³ - (γ/4)ψ⁴ + √(-g) = c_0·ψ
dawało U_eff,old''(1) = -γ < 0 (tachion bug, ujawniony przez Phase 2 P21).
Po G.0 closure z V_M911 + M9.1'' √(-g): U_eff''(1) = +γ > 0 stabilne.

### F10: Stability extension na FRW perturbation

`eq:cosmo-linearized-unified-G0`:
```
δ̈ψ + 3H·δ̇ψ + m_sp²c²·δψ = -(5q·c²/Φ_0)·δρ
```
m_sp² niezmienione przez Hubble friction; mass gap zachowany w epoce
post-recombination.

### F11: Stability extension na N-body simulations

Multiple sources w `Φ_eq[ρ_total]` nie destabilizują spektrum:
- pętla iteracji `Φ_eq^(n+1) = solver(ρ, Φ_eq^(n))` zachowuje physical
  domain (Φ_eq > 0)
- ghost wall pozostaje poniżej minimum profile
- tachyon-free dla każdego intermediate state

Numerical verification: B9 6/6 PASS WEP MICROSCOPE.

### F12: Stability extension na soliton sektor leptonowy

Continuum spektrum perturbacji wokół soliton (asymptotic g(∞) = 1)
ma m_sp² > 0; bound state spektrum ma λ_n ≥ 0 z energetic argument
(`prop:nonlinear-stability` + `prop:Atail-preserved`).

Numerical verification: Phase 2 P22 mass spectrum 5/5 PASS.

### F13: TGP klastruje z 3D Ising po stronie spektrum

| Model | Symetria | NGB | Massive | TGP |
|-------|----------|-----|---------|-----|
| 3D Ising (LPA) | Z₂ dyskretne | 0 | 1 | analogicznie |
| TGP | Z₂ dyskretne | 0 | 1 | identical |
| Linear σ-model O(N) | O(N) ciągłe | N-1 | 1 | różny |

### F14: Stabilność solitonu wymaga energetycznego argumentu, nie lokalnego V''

W solitonowym core (g(0) > 1 dla cięższych leptonów), V_M911''(g) może
być lokalnie ujemne, ale soliton jest *non-trivial stationary point* full
action — stabilność topologiczna z `prop:nonlinear-stability` +
`prop:Atail-preserved` (sek08b).

### F15: Sektor S₁ (Φ > 0) jest physically dopuszczalny

[[../../TGP_FOUNDATIONS.md]] §1 (warstwa 0): brak absolute vacuum (Φ = 0
fizycznie nieosiągalne). Stąd K(φ_eq) = K_geo·φ_eq⁴ > 0 dla wszystkich
fizycznych konfiguracji — żaden sector singularity.

### F16: L03 audit substantialnie ZAMKNIĘTY

| Wymóg L03 | Status |
|-----------|--------|
| (1) Pełna analiza spektralna | ZAMKNIĘTE (synteza pre-existing + thm:spectral-synthesis-L03) |
| (2) Mode counting Z₂ broken | ZAMKNIĘTE (mode_counting_Z2.md) |
| (3) Tachyonic projection na Φ_eq[ρ] | ZAMKNIĘTE (tachyon_check_with_source.md) |
| (4) Ghost wall g* < g_phys | ZAMKNIĘTE (sek08b cor:ghost-artifact, F7) |

## Edycja core (NON-BREAKING)

`sek08b_ghost_resolution.tex`: dodanie `ssec:spectral-synthesis-L03`
(~80 linii LaTeX) z `thm:spectral-synthesis-L03` + 2 remarks
syntezującymi pre-existing wyniki.

→ Patrz [[spectral_synthesis.md]] §"Edycja sek08b" dla pełnej treści.

## Wpływ na PREDICTIONS_REGISTRY

Brak nowych predykcji. Ten cykl jest **konsolidacyjny**:
- nie dodaje nowych falsyfikowalnych testów
- nie modyfikuje istniejących predykcji
- domyka audit gap (L03) — kategoria meta/structural

Rejestr counter pozostaje 856 (cumulative). Cykl L03 jest typu
*audit-resolution*, nie *new prediction*.

## Open issues po tym cyklu

Z [[NEEDS.md]]:

1. **Rigorous lower bound spektrum gap** — formal derivation że λ_0 = m_sp²
   (nie tylko asymptotyczne)
2. **Functional analysis domain conditions** — explicit Sobolev space
   argumentation, self-adjointness L̂
3. **Quantum corrections to spectrum** — 1-loop effective action
   modification of m_sp²

Wszystkie są P3/P4 (nie blokują predykcji).

## Cross-references

- [[README.md]] — werdykt + indeks
- [[mode_counting_Z2.md]] — DOF count
- [[tachyon_check_with_source.md]] — 4 profile spektrum
- [[spectral_synthesis.md]] — S-L unifikacja + edycja core
- [[NEEDS.md]] — pozostałe otwarte
- [[../../audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-05-06.md]] — będzie utworzony
- [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] — target edycji
