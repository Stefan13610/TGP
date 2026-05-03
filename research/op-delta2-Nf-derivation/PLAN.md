---
title: "δ.2 PLAN — Derivation N_f = 5 z TGP first principles + cosmological-gauge bridge"
date: 2026-05-02
cycle: δ.2
status: PLAN — propozycja cyklu, nieuruchomiona
parent: TGP-program portfolio
predecessor:
  - "[[../op-delta1-g-tilde-derivation/README.md]]"
  - "[[../op-gamma1-phi-eff-anchor-resolution/README.md]]"
trigger: δ.1 PARTIAL POSITIVE z H_NF (N_f=5 hypothesis dla "5" w g̃ = 5e²/(12π)). Open problem § 3.4 i § 4.2: pełna derivation wymaga (i) cosmological-gauge bridge, (ii) N_f=5 z TGP first principles
related:
  - "[[../op-N0/]]"
  - "[[../closure_2026-04-26/Lambda_from_Phi0/results.md]]"
tags:
  - TGP
  - delta2
  - N_f-derivation
  - cosmological-gauge-bridge
  - top-decoupling
  - QCD-flavors
  - plan-only
---

# δ.2 — Derivation `N_f = 5` z TGP first principles + cosmological-gauge bridge

> **Status:** PLAN. Nie uruchomiony. Decyzja go/no-go po user review.
>
> **Trigger:** δ.1 H_NF zidentyfikowało `g̃ = N_f·e²/(12π)` z N_f=5 (active flavors at M_Z scale).
> Open problems pozostawione:
> 1. **Czemu N_f = 5** (a nie 4, 6, 3)? Jakie TGP first-principles?
> 2. **Czemu cosmological Λ jest powiązane z M_Z scale** (UV gauge)?
>
> δ.2 ma rozwiązać oba lub document pełen scope problem.

---

## 1. PROBLEM (do rozwiązania)

### 1.1 Co δ.1 zostawiło OPEN

δ.1 H_NF identyfikuje:
$$\tilde{g} = \frac{N_f \cdot e^2}{12\pi}, \quad \Phi_{\text{eff}} = \frac{2 N_f e^2}{3}, \quad \Omega_\Lambda = \frac{N_f \cdot e^2}{2 N_c^3}$$

z **N_f = 5** (active flavors at M_Z scale).

**Niezamknięte:**
1. **N_f = 5 jest empirical** — pochodzi z standard model phenomenology (top mass m_t ≈ 173 GeV decoupled at M_Z)
2. **M_Z scale dla cosmological Λ** — TGP-substrate Λ jest infrared (~H₀), ale g̃ correction jest QCD-flavored (M_Z scale)
3. **Bridge mechanism** — co linkuje cosmological vacuum z gauge sector?

### 1.2 Trzy poziomy "derivation"

**Level A (najambitniejsze):** Pełna first-principles N_f=5 derivation z TGP-substrate.
- Wymaga: derivation top mass m_t z TGP, derivation M_Z, derivation decoupling argument
- ROI: TGP zyskuje closed structural prediction Ω_Λ
- Prawdopodobieństwo: niskie (5%)

**Level B (medium ambition):** N_f=5 z TGP partial structural argument.
- Wymaga: argument że top jest natural exclusion w TGP-substrate (np. instability)
- Argument że M_Z jest natural reference scale (np. EWSB)
- ROI: PARTIAL POSITIVE δ.2 closure
- Prawdopodobieństwo: średnie (30%)

**Level C (najmniej ambitne):** Acknowledge N_f=5 jako empirical input.
- δ.2 closes confirming że N_f=5 jest empirical, ale identifies dlaczego
  taki choice jest natural w TGP context
- ROI: structural justification bez pełnej derivacji
- Prawdopodobieństwo: wysokie (50%)

### 1.3 Dlaczego δ.2 jest valuable nawet at Level C

Even Level C closure jest valuable bo:
- Identyfikuje **specific structural argument** dla N_f=5 (nie czysto numerical)
- Distinguishes between "numerologia" (γ.1 prereframe) i "structural choice" (δ.2 post-closure)
- Daje TGP-program **falsifiable boundary conditions** dla future cycles

---

## 2. POTENCJALNE KIERUNKI DERYVACJI

### Hipoteza H_decouple (PRIMARY)

**Założenie:** N_f = 5 z natural decoupling argument:
- Top quark m_t ≈ 173 GeV >> M_Z = 91.2 GeV
- Top decouples z RG running below M_Z
- Active flavors at M_Z: u, d, s, c, b → N_f=5

**TGP first-principles question:**
- Czy TGP deriviuje hierarchy m_t >> M_Z?
- Czy TGP wybór "M_Z jako reference scale" jest natural?

**Sub-hipotezy:**
- **D1.a:** M_Z jest TGP natural bo electroweak symmetry breaking
- **D1.b:** M_Z jest TGP natural bo coupling unification scale
- **D1.c:** M_Z jest TGP natural bo Higgs VEV scale

### Hipoteza H_topfree (SECONDARY)

**Założenie:** Top quark jest **strukturalnie excluded** z TGP-substrate vacuum:
- Top decays at lifetime ~10⁻²⁵ s < QCD hadronization scale
- W TGP-substrate, vacuum coherence wymaga że soliton istnieje long enough do bound state
- Top "fails" tę test → excluded z N_f counting

**TGP first-principles question:**
- Czy TGP może derive top instability z R3 ODE?
- Czy istnieje natural mass cutoff w TGP-substrate vacuum?

### Hipoteza H_dim (ALGEBRAIC)

**Założenie:** N_f=5 z algebraic structure SU(3) lub TGP-substrate:
- 5 = dim(SU(3)) − dim(Cartan) = 8 − 3 (off-diagonal generators)
- 5 = N_c² − (N_c−1) = 9 − 4 (algebraic combination)
- 5 = N_c + 2 = N_c + dim(diagonal pair)

**TGP first-principles question:**
- Czy któraś z decompositions matche z R3 ODE algebra?
- Czy TGP-substrate ma natural 5-dimensional tensor structure?

### Hipoteza H_geom (BRIDGE)

**Założenie:** Cosmological-gauge bridge przez geometric ratio:
- Φ_eq = H₀ (substrate macro-scale, T-Λ)
- M_Z = ?? (UV anchor)
- Ratio H₀/M_Z ~ 10⁻⁴¹

**TGP first-principles question:**
- Czy istnieje natural bridge formula H₀ ↔ M_Z?
- Czy substrate ma "two-scale" structure?

---

## 3. KORZYŚCI (jeśli się powiedzie)

### 3.1 Direct (Level A or B)

1. **δ.1 H_NF zyskuje pełną derivation** — Ω_Λ = 5e²/54 staje się derived
2. **TGP-program zyskuje structural Λ-derivation** z N_f, N_c, e²
3. **Cosmological-gauge bridge** documentowany (open problem TGP od T-Λ)

### 3.2 Strukturalne

1. **TGP zyskuje multi-scale framework** (cosmological + electroweak coupling)
2. **Foundation dla future Λ-prediction tests** (DESI/Euclid)
3. **Closure dla TGP-program quark sector** (z connection do op-N0)

### 3.3 Predykcyjne

- Jeśli H_decouple wins: TGP musi predict m_t/M_Z ratio (testable)
- Jeśli H_topfree: TGP predict natural mass cutoff
- Jeśli H_dim: TGP predict alternative N_f scenarios (4 lub 6 below/above different scales)

---

## 4. RYZYKA (uczciwa lista)

### 4.1 Krytyczne

**R1. N_f=5 jest fundamentally empirical (60%).**
   Top mass m_t ≈ 173 GeV jest measured value, nie derivable z TGP-substrate.
   Bez derivation m_t, decoupling argument przy M_Z jest empirical.
   **Mitigacja:** Level C closure (acknowledge empirical) zaakceptować jako valuable.

**R2. Cosmological-gauge bridge jest open problem QFT (50%).**
   To jest hard problem — TGP może nie rozwiązać go bezpośrednio.
   **Mitigacja:** Document problem precisely, identify gdzie δ.2 zatrzymuje się.

**R3. δ.2 wymaga op-N0 closure (40%).**
   Jeśli N_f wymaga N_c=3 derived, a op-N0 (N_c=3) nie jest closed,
   δ.2 ma blocked dependency.
   **Mitigacja:** Sprawdź op-N0 status w P1.

### 4.2 Duże

**R4. δ.2 może powtórzyć δ.1 fate** (PARTIAL bez pełnej derivacji).
   **Konsekwencja:** Acceptable — TGP-program continues z PARTIAL chain
   (λ.1 → γ.1 → δ.1 → δ.2 z incremental progress).

**R5. Hipoteza H_topfree może wymagać top mass derivation z R3 ODE α=2.**
   To samo jak λ.1 X = e²/2 wymóg — może NEG.
   **Mitigacja:** Sprawdź czy R3 ODE może give mass cutoff naturally.

### 4.3 Małe

**R6. Sukces δ.2 wymaga update sek00/sek09 z derivation block.**
   Procedura już exists z γ.1 i δ.1 implementations.

---

## 5. CYKL δ.2 — proponowana struktura

### Phase 1 (Foundation) — Map context

**P1.1 op-N0 status check:** Czy op-N0 ma N_c=3 derivation? Czy N_f related?

**P1.2 sek09 SU(3) deep read:** explicit appearances 5 (off-diagonal), N_c+2, decoupling.

**P1.3 Cosmological-gauge bridge survey:** sek00, sek08a, sek08c — does TGP have natural M_Z derivation?

**P1.4 Top mass / instability literature scan:** Czy TGP ma derivation top mass lub instability?

**GATE:** Co najmniej 1/4 hipotez (H_decouple, H_topfree, H_dim, H_geom) ma plausible
TGP-internal foundation. Inaczej Level C closure.

### Phase 2 (Hypothesis testing)

**P2.1 H_decouple test:** czy M_Z jako reference scale ma TGP-natural argument?

**P2.2 H_topfree test:** czy top exclusion z TGP-vacuum coherence?

**P2.3 H_dim test:** czy 5 = SU(3) algebraic decomposition?

**P2.4 H_geom test:** czy ratio H₀/M_Z ma TGP-derivation?

**GATE:** ≥1 hipoteza zwraca structural argument (nie czysto numerical).
Inaczej Level C closure.

### Phase 3 (Verification)

**P3.1 Sympy verification** preferred hipothesis.

**P3.2 Cross-check:** czy implication (np. m_t/M_Z) dają correct values?

**P3.3 Cosmological-gauge bridge document:** explicit RG argument lub structural identification.

### Phase 4 (Consolidation)

**P4.1 Level A (POSITIVE):** pełna derivation, update sek00/sek09 z derivation.
**P4.2 Level B (PARTIAL):** structural argument identified, document boundaries.
**P4.3 Level C (acknowledge):** N_f=5 stays empirical, ale TGP-context for choice
documented.

---

## 6. KOSZT — estimated effort

| Phase | Estimated time | Resources |
|-------|---------------|-----------|
| P1 (Foundation) | 1 dzień | Comprehensive grep, literature scan |
| P2 (Hypothesis testing) | 1-2 dni | sympy, structural analysis |
| P3 (Verification) | 0.5-1 dnia | Cross-check, bridge construction |
| P4 (Consolidation) | 0.5 dnia | Updates lub Level C documentation |
| **Total** | **~3-4 dni focused** | |

---

## 7. DECYZJA — kryterium go/no-go

**GO** jeśli:
- Plan się spina structurally
- P1 znajduje co najmniej 1 plausible foundation
- User chce continue chain λ.1 → γ.1 → δ.1 → δ.2

**NO-GO** jeśli:
- P1 ujawnia że wszystkie hipothesy require op-N0 closure first
- Top mass derivation jest blokujące (out of scope)

**SOFT NO-GO** (revisit later):
- Po Level C closure z δ.2, future cycles (op-N0 success) mogłyby reopen δ.2

---

## 8. RELACJA do innych cykli

- **δ.1** — δ.2 directly follows; δ.1 H_NF interpretation pozostaje, δ.2 deepens
- **γ.1** — algebraic identification stays, δ.2 dodaje structural arg
- **λ.1** — P2.3 reframed via δ.1, δ.2 może give pełne POSITIVE P2.3
- **op-N0** — potential dependency dla N_c, also dla N_f

---

## 9. STATUS PLANU

**Niniejszy dokument:** PLAN. **Nieuruchomiony.** Czeka na user GO.

**Oczekiwane outcomes:**

| Level | Probability | Outcome |
|-------|-------------|---------|
| **A (full)** | 5% | Pełna derivation N_f=5; sek00/09 z derivation |
| **B (partial)** | 30% | Structural argument bez pełnej derivacji |
| **C (acknowledge)** | 50% | Empirical input + TGP-context documented |
| **NEG** | 15% | Wszystkie hipotezy fail, δ.2 closes NEG |

**W każdym przypadku:** δ.2 daje closure dla open problem δ.1 §4.2.

---

## 10. SAMOOCENA planu (uczciwa)

**Strengths:**
- Konkretny problem z δ.1 chain
- 4 hipothesis testable
- Realistic gate criteria (Level C acceptable)

**Weaknesses:**
- **R1+R2 razem (combined ~75%)** mogą zmusić Level C closure
- Top mass derivation może być blokujące (out of scope)
- Cosmological-gauge bridge jest hard QFT problem

**Honest verdict:** δ.2 jest **probably PARTIAL or Level C** ale daje
incremental progress w chain λ.1 → γ.1 → δ.1 → δ.2. Worth doing dla closure
δ.1 open problems, even jeśli result jest acknowledgment empirical.

**Powinien być zrobiony jeśli user chce:**
1. Document TGP scope/boundaries dla N_f=5 question
2. Identify gdzie pełna derivation by needed (foundation dla future cycles)
3. Continue progressive structural deepening
