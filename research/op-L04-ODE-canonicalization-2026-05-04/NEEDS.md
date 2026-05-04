---
title: "NEEDS — L04 ODE-canonicalization (open problems)"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: research/op-L04-ODE-canonicalization-2026-05-04
tags:
  - needs
  - L04
  - open-problems
  - X-e-squared
  - Hobart-Derrick
---

# NEEDS — L04 (otwarte problemy)

## Otwarte luki

### N1: Analytical derivation X = e²/4

**Status:** OPEN — najpoważniejsze.

**Co wiemy:** n(α) = e²·(1−α/4) jest **empirical discovery** z residuum
<0.1% w fit dla α∈[0.25, 4.0]. Numerycznie n(α=2) = e²/2 jest dokładne
do 0.001%.

**Co brakuje:** analitycznie wyprowadzić X = e²/4 z fundamentalnej fizyki.
Trzy hipotezy ([[../why_n3/PHASE2_n_alpha_derivation.md]] §5):

- **H1 RG fixed point**: e² jako anomalous dimension wave-function
  renormalization Z(g₀, α) na fixed point. Phase 6 Q5 R⁵-bridge:
  **NEGATIVE** ([[../why_n3/PHASE6_Q5_R5_bridge_first_attempt.md]]
  lin. 170-179: X = e²/4 downgraded do "leading candidate, not definitive")
- **H2 Gaussian/exponential integral**: e^x z classical Gaussian-weighted
  core integrals z odpowiednią szerokością
- **H3 Topological winding**: RP² Berry phase π z [[../why_n3/PHASE3_RP2_defect_quantization.md]]
  → e^(2π) lub similar

**Source:** [[../why_n3/r3_phase2b_X_constant.py]] (36 candidates scan,
e² leader),
[[../why_n3/PHASE6_Q5_R5_bridge_first_attempt.md]] (NEGATIVE attempt)

**Kandydat dostawcy:** dedicated cycle `op-L04b-X-e2-derivation`
(propozycja: 4-6 tygodni, formal physics)

**Typ:** derivation (analytical from RG flow / Gaussian integral / topology)

### N2: Hobart-Derrick balance point α=4

**Status:** OPEN — n(4) = 0 numerically exact, ale brak formalnej derywacji.

**Co wiemy:** w skanie ODE dla α=4, n(α) = 0 numerycznie. To znaczy że
mass formula degeneruje do `m_obs = c · A²` (czyste, no core dressing).

**Co brakuje:** formal Hobart-Derrick analysis z TGP `K(φ) = K_geo·φ^(2α)`.
Standard Hobart-Derrick scaling argument w d wymiarach daje warunek na
stable solitony z `K=const`; dla `K=φ^(2α)` Derrick scaling powinno dać
explicit α_Derrick(d).

**Hipoteza:** α_Derrick(d=3) = 4 — jeśli prawda, sugeruje że TGP-canonical
α=2 jest *podawalne* w sensie kinetycznego balance (między „too kinetic"
α<2 a Derrick-degenerate α≥4).

**Source:** [[../why_n3/r3_phase2_n_alpha_derivation.txt]] α=4 scan + audit
hipoteza [[../../audyt/L04_ODE_dualism_alpha/README.md]]

**Kandydat dostawcy:** dedicated formal physics cycle `op-Hobart-Derrick-TGP`

**Typ:** derivation (Derrick-scaling z `K(φ) = φ^(2α)`)

### N3: Convergence argument k(α, d) generalization

**Status:** OPEN — argument konwergencyjny z LP-4 nie jest universal.

**Co wiemy:** LP-4 (PLAN_DOMKNIECIA) zamknął się 9/9 PASS na **k=4** z
trzech argumentów:
- Zero-mode (E_kin = E_pot virial)
- Konwergencja wymiarowa: `k = 2(d-1)/(d-2) = 4` w d=3
- Dyskryminacja: tylko k=4 daje r₂₁∈[200,210]

ALE: LP-4 used substrate ODE K=g² (α=1). Dla α=2 (TGP-canonical) wzór
empirycznie p=3 (z r3_alpha2_full_closure), nie k=4.

**Co brakuje:** **uogólniona** relacja konwergencyjna `k(α, d)`. Dla TGP
w d=3 z `K=g^(2α)`, warunek skończoności energii powinien dać
explicit `k(α, 3)`.

**Hipoteza:** k(α, d) = 2(d−1)/(d−2α) (?) — daje k=4 dla α=1 i k=indef
dla α=d/2. Wymaga formal check.

**Source:** [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] LP-4 +
[[../why_n3/r3_observable_vs_full_mass.txt]] α-scan

**Kandydat dostawcy:** dedicated mass-formula derivation cycle
`op-L05-mass-exponent-formal`

**Typ:** derivation (asymptotyczna analiza energy convergence)

### N4: Single-exponent vs two-exponent reconciliation

**Status:** OPEN — empirical p(α)=5−α very close ale nie exact dla wszystkich α.

**Co wiemy:** Z [[../why_n3/r3_observable_vs_full_mass.txt]]:

| α | p_emp | 5−α | diff% |
|---|-------|-----|-------|
| **1.00** | **4.001** | **4.000** | **+0.02%** |
| **2.00** | **3.001** | **3.000** | **+0.02%** |
| 0.50 | 4.809 | 4.500 | +6.87% |
| 2.50 | 2.668 | 2.500 | −6.3% |

`p(α)=5−α` jest dokładne *tylko* dla α=1 i α=2; rozmywa do 6% dla
ekstremów.

Phase 2 dwuwykładnikowa (n(α)=e²·(1−α/4)) jest dokładna do <0.1% dla
całego α∈[0.25, 4.0].

**Co brakuje:** analytical bridge między single-exponent `p(α)=5−α` a
two-exponent `m=c·A²·g₀^n(α)`. Czy istnieje *perturbacyjna ekspansja*
która daje p(α)=5−α jako leading-order, a Phase 2 jako full?

**Source:** [[../why_n3/r3_p_alpha_analytical.txt]] §6-7

**Kandydat dostawcy:** część cyklu N3 (mass-exponent-formal)

**Typ:** analytical-bridge (single ↔ two-exponent decomposition)

### N5: Stability of TGP solitons (Derrick theorem)

**Status:** OPEN — niezaadresowane w R3 i mass_scaling_k4.

**Problem:** Twierdzenie Derricka mówi, że w d≥2 brak stabilnych
solitonów skalarnych z `K=const` i kanonicznym `V(φ)`. Soliton 3D
z czystym skalarem powinien być niestabilny.

W TGP: `K(φ) = K_geo·φ⁴` (α=2), więc Derrick scaling jest *modified*.
Standard Derrick: scaling x → λx daje `E(λ) = λ^(d−2)·E_K + λ^d·E_V`.
Dla TGP z K(φ) = φ^(2α): kinetic czynnik staje się `λ^(d−2)·E_K[φ(x)]
· (jak skaluje się sam K(φ))` co dodatkowo komplikuje balance.

**Pytanie:** czy R3 / TGP solitony są:
- stable
- metastable (ze skończonym lifetime przez vacuum decay)
- artefaktem ODE (fake stationary points)

Bez tego nie wiemy czy m_e, m_μ, m_τ są *fizycznymi cząstkami* czy
matematycznymi konstrukcjami.

**Source:** [[../why_n3/CORRECTIONS_2026-05-01.md]] § B.5 (SPECULATIVE/OPEN)

**Kandydat dostawcy:** dedicated `op-derrick-stability-TGP` cycle

**Typ:** derivation (full Derrick analysis dla `K=φ^(2α)`)

### N6: Geometric coincidence g₀_crit(α=2) ≡ ψ_horizon(M9.1'')

**Status:** OPEN — koincydencja na 4 cyfry, ale brak formalnej derywacji.

**Co wiemy:** g₀_crit(α=2) = 1.874 jest *strukturalną* własnością ODE
solitonu (ODE blow-up). ψ_horizon(M9.1'') = 4/3 jest *strukturalną*
własnością metryki (Lorentzian horizon `g_tt = 0`).

Z [[../why_n3/PHASE1_psi_g0_identification.md]]: liniowa identyfikacja
ψ = 0.3814·g + 0.6186 daje:
```
g₀_crit ≡ (4/3 − 0.6186) / 0.3814 = 1.874   ✓ 4 cyfry
```

**Co brakuje:** formal derivation tej linearnej identyfikacji ψ↔g.
Czy to jest *wymuszona* korespondencja (z aksjomatów), czy *koincydencja*?

Jeśli *wymuszona*: TGP-canonical α=2 jest *strukturalnie unique* przez
to że jego bariera ODE pokrywa się z horyzontem metryki M9.1''.

Jeśli *koincydencja*: pozostaje numerologiczna obserwacja.

**Source:** [[../why_n3/PHASE1_psi_g0_identification.md]] +
[[../op-g0-r3-from-canonical-projection/Phase1_results.md]]

**Kandydat dostawcy:** dedicated `op-psi-g-linear-identification` cycle

**Typ:** derivation (linear ψ↔g jako konsekwencja `K=φ⁴` + M9.1'')

## Pytania otwarte

- **Q1**: Czy α=4 (Hobart-Derrick balance) jest *fizycznym ograniczeniem*
  (TGP nie może mieć α≥4) czy tylko *degenerate point* (TGP może mieć α=4
  z m_obs = c·A² czyste)?

- **Q2**: Czy istnieje *fundamental* reason for X = e²/4 (RG / topology /
  Gaussian) lub jest to numerical coincidence?

- **Q3**: Czy LP-6 audit decyzja "obie formulacje równoważne w słabym polu"
  jest fenomenologicznie poprawna, biorąc pod uwagę że R5 K² fail dla α=2?
  (Odpowiedź per L04: NIE — równoważne są tylko predykcje *dwuwykładnikowe*
  Phase 2; specjalne przypadki α=1 vs α=2 dają różne mass formulas
  *dla tej samej metryki słabego pola*)

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | N1 (X = e²/4 RG derivation) | dedicated R⁵-bridge cycle (formal RG / Hobart-Derrick balance) | Phase 6 Q5 NEGATIVE |
| B2 | N2 (Hobart-Derrick α=4) | formal Derrick scaling z K(φ)=φ^(2α) | open |
| B3 | N3 (k(α,d) general) | formal asymptotic energy convergence analysis | LP-4 |
| B4 | N4 (single ↔ two-exponent) | część dedicated mass-formula cycle | r3_p_alpha_analytical |
| B5 | N5 (Derrick stability) | dedicated stability cycle TGP | open |
| B6 | N6 (psi-g linear identification) | dedicated `op-psi-g-linear` cycle | PHASE1 |

## Closed needs (po L04 cyklu)

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| α=1 vs α=2 dualism | jednolita decyzja kanoniczna α=2 | L04 cykl analytical decision | 2026-05-04 |
| LP-6 dictionary status | reframing: α=1 to specjalny case α=2 | mass_formula_unification.md | 2026-05-04 |
| LP-4 k=4 universality | reframing: k=4 to specjalny case α=1 | mass_formula_unification.md §3.1 | 2026-05-04 |
| R5 K² status | downgrade do DERIVATIVE specjalny case α=1 | analytical theorem (R5 ≡ Phase 2 IFF α=1) | 2026-05-02 (R5 bridge), L04 confirm 2026-05-04 |

## Cross-references

- [[README.md]] — cykl L04 indeks
- [[FINDINGS.md]] — eksportowalne wyniki
- [[../../audyt/L04_ODE_dualism_alpha/]]
- [[../why_n3/PHASE6_Q5_R5_bridge_first_attempt.md]] (X = e²/4 NEGATIVE)
