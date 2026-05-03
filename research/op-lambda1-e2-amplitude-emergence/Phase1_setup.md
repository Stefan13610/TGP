---
title: "λ.1 Phase 1 — Structural setup + falsification ranges"
date: 2026-05-01
phase: 1
parent: "[[README.md]]"
status: ACTIVE
tags:
  - TGP
  - lambda1
  - phase1
  - structural-setup
---

# λ.1 Phase 1 — Structural setup + falsification ranges

## Cel Phase 1

Wyłonić **najbardziej obiecujące mechanizmy** dla derywacji e² w R3
amplitude sector poprzez 5 niezależnych sub-tasks. Score gate ≥3/5 PASS
przed Phase 2 (sympy LOCK + analytical derivation).

---

## Sub-tasks (5)

### L1.1 — Algebraic search dla 10/3 w TGP integers

**Cel:** Sprawdzić czy `10/3` z hipotezy `Φ_eff = (10/3)·e²` ma natural
TGP-algebraic origin.

**Metoda:**
- Systematic enumeration TGP integers: {3, 4, 7, 8, 12, 14, 36, 56, 168, 184}
  + lepton/neutrino masy + Φ_0 wartości
- Test rationals z różnym mianownikiem ≤ 10
- Test relacji 10 ↔ stopnie swobody w 4D (T_4 = sum 1..4 = 10, lub
  symmetric 2-tensor 4×5/2 = 10)
- 3 ↔ N_gen lub d=3 spatial

**PASS criterion:** Znalezienie kombinacji TGP-integers z explicit
fizyczną motywacją która daje 10/3 z diff <0.5% Φ_eff/e².

**Plik:** `phase1L1_algebraic_search.py`

---

### L1.2 — Explicit β-function dla R3 ODE jako effective theory

**Cel:** Pochodzić β-function dla coupling w R3 ODE (potential
`V(g) = -1/g - ln(g)` z Faza 1 why_n3) i sprawdzić czy 1-loop daje
γ_φ = e²/4·(4-α).

**Metoda:**
- Linearyzacja R3 ODE wokół vacuum g=1
- Identyfikacja "running coupling" w R3 ODE (α-zalezne)
- Sympy 1-loop self-energy w 4D substrate
- Porównanie γ_φ z empirical n(α)/(4-α) = e²/4

**PASS criterion:** Pokazać że 1-loop β-function R3 ODE daje
γ_φ proporcjonalne do (4-α) (kierunek, nie dokładna stała).

**Plik:** `phase1L2_beta_function_R3.py`

---

### L1.3 — Partition function setup dla R3 amplitude fluctuations

**Cel:** Setup statystyczno-mechaniczny dla fluktuacji amplitudowych
wokół R3 solitonu; sprawdzić czy `e^(-βE)` factor pojawia się
w mass formula naturally.

**Metoda:**
- Hamilton substrate H = -J·Σ(φ_i·φ_j)² (sek10:43-44)
- Linearyzacja wokół vacuum φ = √Φ₀ + δφ
- Z = ∫Dδφ exp(-S[δφ; soliton bg])
- Sprawdzić leading term w limicie classical → quantum

**PASS criterion:** Strukturalne pojawienie się exp(coś) w mass
contribution z fluktuacji.

**Plik:** `phase1L3_partition_function.py`

---

### L1.4 — Cross-test ζ.1 neutrino spectrum dla e²-hint

**Cel:** Sprawdzić czy neutrino mass spectrum (Σm_ν = 59.6 meV)
wykazuje e²-trace analogiczny do charged lepton ratios.

**Metoda:**
- Pobrać ζ.1 neutrino masses (m_1 = 0.80, m_2 = 8.65, m_3 = 50.11 meV)
- Zaaplikować R3 mass formula z α=2 + n(2)=e²/2
- Sprawdzić czy g₀_νi spełnia analogiczne relacje (φ-drabinka? Koide K=1/2?)
- Test residuum vs PDG / NuFit

**PASS criterion:** Neutrino masses match z e²/2 mechanism do <5%
(neutrina są precyzyjnie znane do ~10%, więc próg luźniejszy).

**Plik:** `phase1L4_neutrino_cross_test.py`

---

### L1.5 — Boundary check: dlaczego e² w amplitude, nie w phase

**Cel:** Conceptual/analytical proof dla strukturalnej różnicy między
J_amp i J_phase która wymusza e² selektywnie w amplitude sector.

**Metoda:**
- Analiza dodatekO_u1_formalizacja:405-421 (J_amp vs J_phase)
- Strukturalne różnice: amplitude → real Φ z Z₂, phase → compact θ z U(1)
- U(1) compactness daje 2π winding (kwantyzacja); amplitude nie ma
  analogicznej topologii dla e (e nie jest kwantyzowane)
- Identyfikacja: dlaczego amplitude może mieć continuous limit
  exp(coś), a phase tylko dyskretne 2π·n

**PASS criterion:** Argument analytical (nie tylko numeryczny) dlaczego
e_Euler MOŻE być w amplitude, ale NIE może być w phase sector.

**Plik:** `phase1L5_amplitude_phase_separation.md`

---

## Score gate Phase 1

```
Score = sum (L1.x PASS)
≥ 3/5 PASS → Phase 2 forward
< 3/5 PASS → λ.1 paused, write postmortem
```

**Kill criterion:** Jeśli żaden z mechanizmów nie produkuje e²-trace,
zarzucamy hipotezę λ.1 i flagujemy X = e²/4 jako EMPIRICAL z explicit
"no fundamental derivation found" w R3 dokumentach.

---

## Timing

- L1.1: ~1 dzień (algebraic search, sympy)
- L1.2: ~3 dni (β-function R3 ODE)
- L1.3: ~3 dni (partition function setup)
- L1.4: ~1 dzień (neutrino cross-test, leverage existing ζ.1 data)
- L1.5: ~2 dni (conceptual analysis)

**Total Phase 1: ~10 dni** (lub 2 tygodnie kalendarzowe).

---

**Status:** ACTIVE — startujemy 2026-05-01.
**Order:** L1.1 → L1.4 → L1.5 → L1.2 → L1.3 (od najlżejszych do najcięższych).
