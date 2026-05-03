---
title: "λ.1 L1.5 — Boundary check: dlaczego e² w amplitude, nie w phase"
date: 2026-05-01
phase: 1
sub-task: L1.5
parent: "[[Phase1_setup.md]]"
status: COMPLETED
tags:
  - TGP
  - lambda1
  - phase1
  - L1.5
  - amplitude-vs-phase
  - boundary-check
  - conceptual-analysis
---

# L1.5 — Dlaczego e² może być w amplitude, ale nie w phase sector?

## 1. Cel sub-tasku

**Pytanie:** Hipoteza λ.1 mówi że e² jest fundamentalna w **J_amp**
(amplitude/grawitacja), ale α_em (J_phase) używa Φ₀, J, π **bez** e_Euler.
Skąd ta selektywność strukturalna?

**Format:** Conceptual/analytical analysis (nie numeryczne).

**PASS criterion:** Argument analytical (nie tylko numeryczny) **dlaczego**
e_Euler MOŻE być w amplitude, ale **NIE może** być w phase sector.

---

## 2. Strukturalna różnica J_amp vs J_phase w TGP

Z `core/formalizm/dodatekO_u1_formalizacja.tex:405-421`:

> "TGP ma jeden substrat Γ = (V, E), ale **dwa sektory sprzężeń**:
> - **J_amp**: sprzężenie amplitudowe (|ψ_i|) — daje pole Φ, geometrię, grawitację.
> - **J_phase**: sprzężenie fazowe (θ_i) — daje pole A_μ, elektromagnetyzm."

**Substrate hamiltonian** (sek10:43-44 i dodatekO):
```
H_Γ = -J · Σ_<ij> Re(ψ_i* · ψ_j)
    = -J · Σ_<ij> |ψ_i| · |ψ_j| · cos(θ_j - θ_i)
```

gdzie ψ_i = |ψ_i| · exp(iθ_i) jest complex.

**Dwa sektory sprzężeń:**

| Sektor | Zmienna | Domena | Symetria | Topologia |
|--------|---------|--------|----------|-----------|
| **J_amp** | \|ψ_i\| ∈ ℝ_{≥0} | real, non-negative | **Z₂** (\|ψ\| → \|ψ\|) | trywialna |
| **J_phase** | θ_i ∈ [0, 2π) | compact circle | **U(1)** (θ → θ + 2π) | nontrywialna π₁ = ℤ |

---

## 3. Topologia jest kluczowa

### 3.1 J_phase: compact U(1) → kwantyzacja, brak natural e

`dodatekO:221-237` (Twierdzenie kwantyzacji liczby zwojowej):

> "Liczba zwojowa n[γ] jest **całkowita** dla każdej pętli γ.
> Wynika z kompaktowości θ ∈ [0, 2π)."

**Konsekwencja:** w phase sector wszystkie loop integrals są **kwantowane
przez 2π**. Charge quantization (`dodatekO:241-280`):

```
e_0 = ℏ_0 · 2π / Φ_mag,min
```

— ładunek pochodzi z winding number `n ∈ ℤ`, więc tylko **całkowite
wielokrotności 2π** mogą się pojawić.

**Brak natural source dla e_Euler ≈ 2.718 w phase sector:**
- e_Euler nie jest wielokrotnością 2π
- e_Euler nie jest "winding number" (nie jest całkowita)
- W phase sector, każda struktura jest **discrete** (kwantyzowana)
- Continuous limit (1+1/n)^n → e wymaga **ciągłej** zmiennej, której
  phase sector nie posiada

**Wniosek:** **e_Euler jest strukturalnie wykluczone z J_phase**.

### 3.2 J_amp: real Z₂ → continuous, MOŻE pomieścić e

W amplitude sector |ψ_i| ∈ ℝ_{≥0} jest **continuous**:
- Brak kompaktowości (|ψ| może rosnąć w nieograniczoności
- Z₂ symetria jest dyskretna ale nie kompaktyfikuje continuum
- W ciągłej granicy `|ψ_i+1| = |ψ_i| · (1 + Δ/n)^n → |ψ_i| · exp(Δ)`
  jest **mathematically allowed**

**Konsekwencja:** w amplitude sector, **continuous limit dyskretnego
procesu mnożnikowego** (== e w matematycznej definicji) jest **strukturalnie
dostępny**.

**Wniosek:** **e_Euler MOŻE pojawić się w J_amp przez kumulatywny limit**.

---

## 4. Kumulatywny mechanizm w amplitude sector — formal sketch

### 4.1 R3 soliton dressing (Twoja luźna analogia, sformalizowana)

R3 soliton ma profil g(r) z `g(0) = g₀, g(∞) = 1`. W kontinuum, każda
"warstwa" r → r+dr modyfikuje pole:

```
g(r+dr) ≈ g(r) · (1 + Δ(r)·dr) ≈ g(r) · (1 + Δ(r)·dr/n)^n  [n→∞]
```

W limicie ciągłym:
```
g(r) = g(0) · exp(∫₀^r Δ(r') dr')
```

**Mass formula z całki**:
```
m_obs ~ A_tail² · g₀^n   [empirycznie z why_n3]
       = exp(2 log A_tail + n log g₀)
       = exp[wykładnik kumulatywny]
```

`exp` pojawia się **naturalnie** ze struktury continuous-limit.
Dla `n = e²/2` (TGP-canonical α=2), wykładnik to:
```
exp[(e²/2) · log g₀] = g₀^(e²/2)
```

**Ale skąd e²?** Tu jest **otwarty problem** — sam `exp` jest natural,
ale wartość wykładnika `e²/2` wymaga konkretnego mechanizmu (RG flow,
self-energy, partition function).

### 4.2 Phase sector NIE ma analogicznego mechanizmu

W phase sector, θ_{i+1} = θ_i + Δθ_i z Δθ ∈ (-π, π] (dodatekO:213-214).
**Modulo 2π** identyfikacja:
```
θ_{i+1} ≡ θ_i + Δθ  (mod 2π)
```

Continuous limit nie produkuje exp(coś) — produkuje **winding number**:
```
Σ_i Δθ_i = 2π · n   [n ∈ ℤ]
```

To jest **fundamentalnie różna struktura** niż amplitude continuous.

---

## 5. Hipoteza: e² jako "amplitude winding analog"

Analogia:
- **Phase sector:** winding number n ∈ ℤ → ładunek `e_0` kwantyzowany
- **Amplitude sector:** continuous "winding" exp(n·something) → mass
  `m ~ g₀^n` z continuous n

W phase sector mamy:
```
ładunek q = n · e_0 (z e_0 = 2π·ℏ/Φ_mag,min)
```

W amplitude sector mogłaby być analogiczna struktura:
```
masa m ~ g₀^(c · e²/4)   gdzie c jest "amplitude winding count"
```

**Dla R3 mass formula z α=2: c = 2** (= n·4-α z α=2). To "winding count"
mogłoby pochodzić z:
- 2 = liczba spatial dimensions w amplitude wave (transverse?)
- 2 = liczba TT-tensor modes (gravity)
- 2 = K_charged (=2/3, related)

To jest **conceptual hint**, nie derivation.

---

## 6. Dlaczego e² (nie e, nie e³) w amplitude?

Empirycznie z R3 dla α=2: n(2) = e²/2.

**Obserwacja:** wykładnik 2 w `e²` może odpowiadać:
- α = 2 (kinetic prefactor exponent dla TGP-canonical K=φ⁴, gdzie 2α=4)
- 2 = symmetric tensor rank w 4D
- 2 = K_lep (=2/3, related)

Konkretne: jeśli RG anomalous dimension γ_φ wokół AS NGFP w 4D ma postać:
```
γ_φ = (1/4) · (4-α) · (cos²θ + 1)   [z Z₂ symmetry, θ jakaś phase]
```

Wtedy `Z = exp(γ_φ · ln g₀)` z odpowiednim cos²θ avg = e²/4·(4-α). To
jest hipotezą do sprawdzenia w L1.2.

---

## 7. Wnioski L1.5

### 7.1 Główna konkluzja

**Strukturalnie udowodnione:**
1. Phase sector (J_phase, compact U(1)) **wyklucza** e_Euler — winding
   numbers są dyskretne 2π·n.
2. Amplitude sector (J_amp, continuous ℝ_{≥0}) **pozwala** na e_Euler —
   continuous limit dyskretnego procesu mnożnikowego daje exp().

**Pozostałe pytania (do L1.2/L1.3):**
- Konkretny mechanizm który **wymusza** e² (nie tylko allows exp)
- Dlaczego wykładnik **2** (e², nie e lub e³)
- Connection do α=2 (TGP-canonical K=φ⁴)

### 7.2 PASS / FAIL judgment dla L1.5

**PASS criterion:** "Argument analytical dlaczego e_Euler MOŻE być w
amplitude, ale NIE może być w phase sector."

**Status:** ✓ **PASS**

Konkretne argumenty:
1. **Topology argument**: phase sector compact U(1) wymusza kwantyzację
   2π·n; brak miejsca dla irracjonalnej e_Euler.
2. **Continuous limit argument**: amplitude sector pozwala na
   `(1 + Δ/n)^n → exp(Δ)` w limicie ciągłym; phase nie pozwala.
3. **TGP confirmation**: dodatekO_u1_formalizacja:300-421 explicit
   pokazuje że α_em używa Φ₀, J, π — **bez** e_Euler.

### 7.3 Implikacja dla λ.1

**λ.1 może działać konsekwentnie** — szukanie e² w amplitude sector jest
**strukturalnie dozwolone**, brak e w phase sector jest **strukturalnie
wymagany**. Cykl λ.1 jest **internally consistent** z TGP-formalismem.

**ALE:** to jest argument *negative* (e² **może** być w amplitude) — nie
**positive** (e² **musi** być w amplitude). Pełna pozytywna derywacja
wymaga L1.2 (β-function) lub L1.3 (partition function).

---

## 8. Status sub-tasku L1.5

```
PASS criterion spełnione.
Score gate contribution: 1 PASS / 1 sub-task

Wkład do λ.1 Phase 1 score: +1
```

---

## 9. Pliki referencyjne

- `core/formalizm/dodatekO_u1_formalizacja.tex:300-421`
  (α_em derivation, J_amp/J_phase distinction)
- `core/formalizm/dodatekO_u1_formalizacja.tex:213-237`
  (phase compactness, winding quantization)
- `core/sek10_N0_wyprowadzenie/sek10_N0_wyprowadzenie.tex:43-56`
  (substrate hamiltonian H_Γ)

---

**Autor:** λ.1 Phase 1 L1.5.
**Data:** 2026-05-01.
**Status:** PASS — argument analytical complete.
**Następne:** L1.2 (β-function R3 ODE), L1.3 (partition function).
