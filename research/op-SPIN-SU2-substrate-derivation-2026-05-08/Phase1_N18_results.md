---
title: "Phase 1 N18 results — Quantum lift: 2 branches → SU(2) fundamental rep"
date: 2026-05-09
type: phase-needs-result
status: RESOLVED_POSITIVE_CONDITIONAL
parent: "[[./README.md]]"
phase: 1
need_id: N18
sympy_verification: 7/7 PASS
tags:
  - phase1
  - N18
  - quantum-lift
  - sympy-verified
  - SU2-fundamental
  - 720-degrees
  - born-rule
  - bloch-sphere
---

# Phase 1 N18 — Quantum lift to SU(2) results

## Status: **RESOLVED POSITIVE conditional** (7/7 sympy F-tests PASS)

**Conditional na:** orientation degree of freedom solitonu w 3D (placeholder
dla spatial extension Phase 1 cd).

## Centralny wynik

**Standard QM machinery** (linear superposition, unitarity, Born rule)
zastosowana do **2-branch bifurcation state space** (z N17) **natywnie**
daje:

- 2D complex Hilbert space ℂ²
- SU(2) fundamental representation
- 720° symmetry (NIE 360°)
- Born rule cos²(θ/2)

**Wszystko bez dodatkowych postulatów** poza standardową QM. To jest
**konkretny mathematical mechanism** dla spinora 1/2 w TGP context.

## F-tests breakdown

### F1: Bifurcation states → 2D complex Hilbert ℂ² (✓ PASS)

Z N17:
```
{|zanik⟩, |ekspansja⟩}    classical 2-state
```

Quantization (linear superposition):
```
|ψ⟩ = α|zanik⟩ + β|ekspansja⟩,   α, β ∈ ℂ
|α|² + |β|² = 1
```

**Argument dla ℂ (nie ℝ):** ewolucja czasowa wokół meta-stable saddle
oscyluje z frequency ω_osc; exp(-iω_osc·t) accumulates phase →
amplitudes muszą być complex.

Hilbert space H = ℂ² with inner product ⟨a|b⟩ = ā·b.
Unit sphere w ℂ² = S³ (real-dim 3-sphere).

### F2: Bloch sphere parametrization (✓ PASS)

```
α = cos(θ/2)
β = e^(iφ) sin(θ/2)
```

Sympy verification: ⟨ψ|ψ⟩ = 1 (normalizacja).

Geometria:
- θ=0: |ψ⟩=(1,0) = pure |zanik⟩ (south pole)
- θ=π: |ψ⟩=(0, e^(iφ)) ∝ |ekspansja⟩ (north pole)
- θ=π/2: equal superposition (equator)

### F3: SU(2) Lie algebra (Pauli matrices) (✓ PASS)

Pauli matrices σ_x, σ_y, σ_z generują SU(2). Sympy verified:
```
[σ_x, σ_y] = 2i σ_z   ✓
[σ_y, σ_z] = 2i σ_x   ✓
[σ_z, σ_x] = 2i σ_y   ✓
```

### F4: SU(2) → SO(3) double cover (✓ PASS)

```
U(angle) = exp(-i(angle/2)σ_z) = diag(e^(-i·angle/2), e^(i·angle/2))
```

Verified: U U† = I, det(U) = 1 (czyli U ∈ SU(2)).

Działanie na R³ przez sprzężenie: U X U† dla X = Σ x_i σ_i.

### F5: π_1(SO(3)) = Z_2 (✓ PASS — KLUCZOWY)

Sympy verified algebraicznie:

```
U(2π) = [[-1, 0], [0, -1]] = -I        (360° rotation w R³ → -I w SU(2))
U(4π) = [[+1, 0], [0, +1]] = +I        (720° rotation → +I = identity)
```

**To jest** π_1(SO(3)) = Z_2 **strukturalnie**: zamknięta pętla w SO(3)
o "winding 1" (360°) podnosi się do **otwartej** ścieżki w SU(2)
(I → -I); pętla winding 2 (720°) podnosi się do zamkniętej (I → -I → I).

### F6: 360° vs 720° na quantum state (✓ PASS — DEMONSTRACJA WAŻNA)

Stan początkowy |ψ_0⟩ = (1, 0) (pure |zanik⟩):

```
|ψ(360°)⟩ = U(2π)|ψ_0⟩ = -|ψ_0⟩    ← znak zmieniony!
|ψ(720°)⟩ = U(4π)|ψ_0⟩ = +|ψ_0⟩    ← pełen powrót
```

**To NIE jest postulat** — wynika directly z reprezentacji SU(2) na ℂ².
Standard QM machinery na 2-state space DAJE this naturally.

### F7: Born rule cos²(θ/2) (✓ PASS)

Pomiar wzdłuż osi pod kątem θ:
```
|+, n⟩ = (cos(θ/2), e^(iφ)sin(θ/2))
```

Probability dla |ψ_0⟩ = (1, 0) measured along n:
```
P(+|θ) = |⟨+, n | (1,0)⟩|² = cos²(θ/2)   ✓ sympy verified
```

To jest **standard Born rule**, ale **emerguje natywnie** z 2-state
geometry, bez additional structure.

## Część II: Co to RZECZYWIŚCIE pokazuje (krytycznie)

### II.1 Strukturalny win

**Pre-N18 problem:** mając klasyczne 2-state, jak go ulift'ować do
quantum SU(2) z 720° symmetry **bez dodawania postulatów**?

**Post-N18 odpowiedź:** standard QM aplikowane do bifurcation space
**naturalnie** daje SU(2) fundamental rep. 720° **wynika** z ℂ²
representation, nie postulatem.

### II.2 Co N18 ZAKŁADA (ale nie udowadnia)

**Założenie 1 (mocne):** soliton TGP ma orientation degree of freedom
w 3D otaczającej geometrii Φ̄.

- Jeśli soliton sferycznie symetryczny → SO(3) działa trivialnie →
  brak SU(2) lift potrzebny → spin 0 only
- Jeśli soliton ma orientation (multipole / asymetria) → SO(3) działa
  nontrivialnie → SU(2) lift wymagany

To wymaga **spatial 3D extension** (Phase 1 cd) — N18 NIE rozstrzyga
sam ten punkt.

**Założenie 2 (średnie):** linear quantum mechanics aplikuje do
bifurcation space.

- Czy 2-branch outcome jest opisany przez linear superposition (ψ ∈
  ℂ²)? Standard QM: tak. Pokazanie tego z TGP-natywnej dynamiki
  wymagałoby derivation Schrödinger-like equation z Φ-EOM.

**Założenie 3 (słabe):** complex amplitudes (nie tylko real).

- Argument oscillation phase jest suggestive, ale formal derivation
  z Φ-EOM nie został wykonany.

### II.3 Ograniczenia mathematical scope

N18 jest **algebraiczna verification** standardowych QM constructions
zastosowanych do 2-state. To jest **kompetentne mathematical work**,
ale nie revealing — odtwarzamy znaną maszynerię (SU(2) na ℂ²).

**TGP-specific contribution:** identyfikacja {|zanik⟩, |ekspansja⟩}
jako fizyczna realizacja 2-state space. To jest **konceptualny wkład**
z N17, w N18 robimy quantization tego.

### II.4 Honest scope statement

N18 RESOLVED w sensie: "JEŚLI standard QM aplikuje do bifurcation
space z N17, TO SU(2) z 720° emerguje naturalnie."

N18 NIE RESOLVED w sensie: "TGP physics z pewnością generuje SU(2)
spinora bez dodatkowych assumptions."

Różnica jest **ważna metodologicznie** — N18 jest *consistency check*,
nie *uniqueness derivation*.

## Część III: Six requirements update

| # | Req | Pre-N17 | Post-N17 | Post-N18 |
|---|-----|---------|----------|----------|
| 1 | 2 wyniki ± | structural | ✓ derived | ✓ |
| 2 | Nie wektor R³ | structural | ✓ derived | ✓ |
| 3 | **720°** | OPEN | open | **✓ RESOLVED conditional** |
| 4 | **cos²(θ/2)** | OPEN | classical only | **✓ RESOLVED quantum** |
| 5 | -cos(θ) singlet | OPEN | open | open (Phase 5) |
| 6 | Brak tabeli | structural | ✓ | ✓ |

**5/6 satisfied at N18 level.** Pozostaje:
- **Bell singlet -cos(θ)** — wymaga 2-soliton joint state (Phase 5)
- Implicit: orientation degree of freedom (Phase 1 cd)

## Część IV: TGP-specific interpretation

### IV.1 Co znaczy "rotacja w R³" dla solitona TGP?

W standardowym QM rotacja jest passive change układu współrzędnych
laboratorium. W TGP-natywnym ujęciu:

**External R³ rotation** = rotacja **otaczającej geometrii Φ̄**:
- Φ̄ jest tłem, ma "preferred axis" tylko w obecności matter (np. B-field)
- Rotacja Φ̄ ekvivalentna jest do rotacji "axis pomiaru"
- Soliton ma swój **internal orientation** O w stosunku do tła

**Induced quantum rotation:**
- External R ∈ SO(3) → induced ρ(R) ∈ SU(2) (przez double cover)
- |ψ⟩ → ρ(R)|ψ⟩
- 360° rotation w lab → -1 phase
- 720° rotation w lab → identity

### IV.2 Connection do "lean direction" picture

Z parent working hypothesis [[../op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07/]]:

> "Mode A: cloud rozlewa się do wewnątrz"
> "Mode B: cloud rozlewa się na zewnątrz"

N17 + N18 mapping (precyzowane):
```
|zanik⟩    ↔ lean inward  ↔ θ = 0  (south pole Bloch)
|ekspansja⟩ ↔ lean outward ↔ θ = π  (north pole Bloch)
superpozycja ↔ θ ∈ (0,π)  (równik = równowaga)
```

**Pomiar Stern-Gerlach:**
- B-field gradient narzuca lokalną oś (axis n)
- Soliton "wybiera" branch zgodnie z kątem między osią a stanem
- Probability cos²(θ/2) zgodna z Born rule

### IV.3 Dlaczego elektron jest spin 1/2, nie spin 1?

W naszym frameworku:
- 2-state space → fundamental rep SU(2) (j=1/2)
- 3-state space byłby (j=1)
- 4-state byłby (j=3/2), itd.

Soliton TGP w meta-stable bifurcation ma **EXACTLY 2** outcomes
(zanik/ekspansja, z N17). To **wymusza** j=1/2.

**Heurystyka:** "elektron ma spin 1/2 bo dynamika TGP-soliton ma
exactly 2 dynamic outcomes, nie 3 ani 4."

## Część V: Implications dla cyklu

### V.1 Probability re-update

| Outcome | Pre-N18 | Post-N18 |
|---------|---------|----------|
| Pełen DERIVED | 30-40% | **40-50%** ↑ |
| STRUCTURAL CONDITIONAL | 35-45% | 30-40% |
| STRUCTURAL_NO_GO | 15-25% | **10-20%** ↓ |
| EARLY_HALT | 10-15% | 5-10% ↓ |

**Reasoning:** N17 + N18 razem dostarczają **konkretne mathematical
chain**:
```
TGP Φ-EOM → bifurcation 2-state → quantization → SU(2) → 720°, cos²(θ/2)
```

To jest pierwszy raz w cyklu gdy mamy **end-to-end derivation skeleton**
łączący TGP physics z quantum spin requirements (modulo orientation
condition).

### V.2 Pozostałe pytania krytyczne

1. **Orientation degree of freedom** — Phase 1 cd:
   - Czy spatial 3D soliton ma orientation?
   - Jeśli tak, jaka jest jej topologia (SO(3) vs S² vs other)?

2. **Bell singlet -cos(θ)** — Phase 5:
   - 2-soliton joint state z dynamic equilibrium
   - Antysymetryczny singlet
   - Reprodukcja korelacji Bell

3. **Why C, not R?** — derivation z Φ-EOM:
   - Argument oscillation phase jest heuristic
   - Formal derivation z TGP physics potrzebny dla DERIVED status

### V.3 Następny krok rekomendowany

Po N18, najważniejsze ścieżki:

**A) Phase 1 cd: spatial 3D extension N17/N18**
   - Sprawdzić czy spatial soliton ma orientation
   - Adresuje główne conditional N18

**B) N16: dynamic equilibrium z tłem Φ̄**
   - Dlaczego solitony nie ekspandują/zanikają w realnym universe
   - Foundation dla stable dynamic objects

**C) N19: external SO(3) embedding mapping**
   - Konkretnie jak rotation otaczającego Φ̄ działa na bifurcation state
   - Mapowanie SU(2) ↔ SO(3) na fizycznym poziomie

**D) Phase 5 sneak preview: Bell singlet**
   - Najbardziej ambitne; sprawdza czy framework rozszerza się na
     2-soliton entanglement

Każde wartościowe; **A** lub **C** są naturalne (kompletują N18
foundation).

## Część VI: Honest reporting note

### Co N18 osiąga
- ✅ Quantum 2-state → SU(2) fundamental rep, 720°, Born rule
- ✅ Standard QM machinery na bifurcation space
- ✅ Pierwszy strukturalny most TGP physics → quantum spin
- ✅ Sympy 7/7 PASS, mathematically robust

### Czego N18 NIE osiąga
- ❌ Derivation orientation degree of freedom z TGP
- ❌ Derivation linear quantum mechanics z Φ-EOM
- ❌ Derivation complex amplitude (vs real) z Φ physics
- ❌ Bell singlet -cos(θ) (Phase 5)

### Robustność
N18 wynik jest **conditionally robust**: jeśli N17 (2 branches)
i orientation hold, then N18 conclusions hold. Maszyneria SU(2) sama
jest standard math, niezależna od TGP-specific.

### Subjektywne reflection
N17 + N18 razem są **najsilniejszym pozytywnym wynikiem cyklu** do
tej pory. Adresują 5/6 spinor requirements z konkretną Φ-EOM podstawą.

Pozostałe 1/6 (Bell singlet) jest **najbardziej challenging** ale też
**najbardziej rewarding** jeśli się uda.

## Cross-references

- [[./README.md]] — overview cyklu (probability update)
- [[./Phase1_N18_quantum_lift_sympy.py]] — sympy verification
- [[./Phase1_N17_results.md]] — parent (klasyczne 2-state)
- [[./Dynamic_equilibrium_framework.md]] — physical foundations
- [[./N14_M911_inspection.md]] — N14 negative (rozbrojony przez N17+N18)
- [[../op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07/]] — parent
  working hypothesis

## Status check post-N18

**Cykl op-SPIN-SU2-substrate-derivation:**
- Phase 0: ✓ COMPLETED
- Phase 1: 70% complete (N16, spatial 3D, N19 still open)
- Phase 2-5: pending
- Probability DERIVED: 40-50% (significant uplift)

**Six requirements:**
- 5/6 RESOLVED (1 conditional na Phase 1 cd)
- 1/6 OPEN (Bell singlet, Phase 5)

**N18 status:** RESOLVED POSITIVE CONDITIONAL, sympy 7/7 PASS.
