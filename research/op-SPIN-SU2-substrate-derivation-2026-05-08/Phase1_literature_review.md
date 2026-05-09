---
title: "Phase 1 Literature Review — solitony w real scalar field, Derrick threat, TGP rescue paths"
date: 2026-05-08
type: phase-literature-review
status: WIP
parent: "[[./README.md]]"
phase: 1
tags:
  - phase1
  - literature-review
  - solitons
  - derrick-theorem
  - skyrme
  - oscillons
  - tgp-rescue
---

# Phase 1 Literature Review

## Cel review

Przed numerical scan dla N2 (existence non-spherical solitons), trzeba
ustalić **co już wiadomo** o solitonach w real scalar field. Konkretnie:

1. Czy w 3D z czystym single real scalar Φ + standard kinetic + V(Φ)
   istnieją stable static localized solutions?
2. Jakie struktury dodatkowe pozwalają na takie solitony?
3. Które z tych struktur są **TGP-natywne** (nie wymagają nowych pól)?
4. Które ansätze konkretnie testować w Phase 1 numerical scan?

## ⚠️ KEY THREAT: Derrick's theorem (1964)

**Statement (Derrick 1964):** Niech L = ½(∂_μ φ)² - V(φ) z V(φ) ≥ 0,
real scalar field, D_s spatial dimensions. Wtedy **żadna** static
finite-energy localized solution NIE jest stable względem rescaling
x → λx.

### Dowód (krótki)

Energia statyczna:
```
E[φ] = T[φ] + U[φ]
T[φ] = ½ ∫ d^{D_s}x (∇φ)²
U[φ] = ∫ d^{D_s}x V(φ)
```

Pod rescalingiem φ_λ(x) ≡ φ(λx):
```
T[φ_λ] = λ^{2-D_s} T[φ]
U[φ_λ] = λ^{-D_s} U[φ]
E(λ) = λ^{2-D_s} T + λ^{-D_s} U
```

Warunek stacjonarności:
```
dE/dλ |_{λ=1} = (2-D_s)·T - D_s·U = 0
```

**Dla D_s = 3:**
```
(2-3)·T - 3·U = 0  ⇒  -T = 3U  ⇒  U = -T/3
```

Ponieważ T ≥ 0 strukturalnie i U ≥ 0 z założenia V ≥ 0, **JEDYNA**
możliwość to T = U = 0, czyli trywialny vacuum φ = const.

**Konkluzja:** W 3D z standardowym kinetic + V ≥ 0, **NIE istnieją**
stable static finite-energy solitony.

### Implication dla TGP Path A

To jest **bezpośrednie zagrożenie** dla H1:
- Path A wymaga localized δΦ-soliton z orientation degree of freedom
- Derrick mówi: w 3D ze standardowym Lagrangianem **takie solitony
  nie istnieją**
- Bez solitonów, nie ma moduli space, nie ma SU(2), całość **upada**

**Surowy honest assessment:** jeśli TGP framework używa standardowego
kinetic term i V(Φ) ≥ 0, **Path A jest STRUCTURAL_NO_GO już na poziomie
Derricka**.

## Workarounds Derricka — landscape

Literatura zna **kilka** sposobów ominięcia Derricka:

### W1: 1D solutions (ścieżka trywialna, niezdolna w 3D)
- φ⁴ kink: φ(x) = v·tanh(x/δ) — stable w 1D
- sine-Gordon: φ(x) = 4·arctan(exp(x/δ))
- **Problem:** strictly 1D, nie generalizują na 3D bez dodatkowej struktury

### W2: Multi-component scalar field (Skyrmion path)
**Skyrme (1961, 1962):** φ ∈ S³ (multi-component), Lagrangian z
**higher-derivative term** (Skyrme term ~ (∇φ × ∇φ)²) plus winding
number stabilizujący.

**Status w TGP:** ❌ wymaga multi-component Φ → narusza S05.

### W3: Q-balls (Coleman 1985)
**Mechanism:** complex scalar field z U(1) symmetry, time-dependent
phase ω, conserved charge Q. Q-ball: φ(x,t) = e^{iωt}·f(|x|).

**Status w TGP:** ❌ wymaga complex scalar → narusza single real Φ.

### W4: Boson stars (Kaup 1968, Ruffini-Bonazzola 1969)
**Mechanism:** complex scalar + gravity. Gravitational binding stabilizes.

**Status w TGP:** częściowo — TGP ma emergent gravity z M9.1''(Φ),
ale wciąż wymaga complex field.

### W5: Oscillons (quasi-stable, time-periodic)
**Bogolyubsky-Makhankov 1976, Gleiser 1994.** Real scalar field + V(Φ)
może mieć **time-periodic** quasi-stable lokalizowane konfiguracje:
```
Φ(x,t) ≈ Φ_0(x) cos(ωt) + harmonics
```

**Properties:**
- NIE są strict static — Derrick nie aplikuje
- Quasi-stable: very long lifetime (>>10^6 oscillation periods),
  ale ostatecznie radiate energy
- Wymagają specific potential (asymmetric V(Φ) z negative quartic
  region)
- Generic w realistic potentials (e.g. φ⁴ Symmetry-breaking)

**Status w TGP:** ✅ **POTENTIALLY VIABLE** — Φ-EOM z M9.1''(Φ) może
mieć asymmetric potential dopuszczający oscillons.

**Critical caveat:** oscillons nie są **strictly stable** — emit
radiation, decay. Ale lifetime może być **kosmologicznie długi** (longer
than universe age dla niektórych parametrów).

### W6: Higher-derivative kinetic (k-essence, Skyrme-like)
**Mechanism:** Lagrangian ma terms wyższego rzędu pochodnych:
```
L = K(X) + L₄((∂φ)⁴) - V(φ)
gdzie X = ½(∂φ)²
```

**Derrick scaling z k-essence K(X) ~ X^n:**
```
T_n[φ_λ] = λ^{2n-D_s} T_n[φ]
```
Dla D_s=3: stable solution wymaga balansu między różnymi terms o
różnych scaling dimensions. Możliwe gdy są **co najmniej dwa terms**
o przeciwnych znakach scaling exponent.

**Klasyczny przykład:** Skyrme model w 3D ma K₂ + K₄ + V, gdzie K₂
(standard kinetic, scaling -1) i K₄ (Skyrme term, scaling +1) balansują,
plus V (potential, scaling -3) wprowadza skala długości.

**Status w TGP:** ✅ **POTENTIALLY VIABLE** — jeśli M9.1''(Φ) wprowadza
**non-canonical kinetic structure**, derivation może omijać Derricka.

**Specific:** jeśli effective Lagrangian ma postać:
```
L_eff = ½ g^{μν}(Φ̄) ∂_μδΦ ∂_νδΦ + α (∂_μδΦ ∂^μδΦ)² + V_eff(δΦ)
```
gdzie g^{μν}(Φ̄) jest non-trivial Φ̄-dependent metric, i α ≠ 0 wynika
z M9.1''(Φ̄) expansion, to mamy **TGP-natywny analog Skyrme**.

### W7: Non-trivial background metric
**Mechanism:** scalar field na curved spacetime z fixed g_μν może mieć
solitony nawet w 3D jeśli geometry zapewnia odpowiednie boundary
conditions.

**Status w TGP:** Φ̄ generuje tło M9.1''(Φ̄). Effective EOM dla δΦ
**działa na tym tle**. To jest **forma W6** (effective non-canonical
kinetic z g^{μν}(Φ̄)), ale fundamentally może dać dodatkowe efekty.

### W8: Time-dependent rotating solutions (spinning Q-balls, vortons)
Wymagają complex field z winding — niedostępne w TGP single-real-Φ.

## Decision tree dla TGP Path A

```
Czy M9.1''(Φ̄) wprowadza non-canonical kinetic dla δΦ?
   │
   ├─ TAK → effective Lagrangian może mieć Skyrme-like balance
   │        → solitony 3D POTENTIALLY EXIST
   │        → Phase 1 numerical scan ma sens
   │
   └─ NIE → standard Lagrangian dla δΦ
            → Derrick FORBIDS static localized 3D solitons
            → Tylko opcje:
              ├─ Oscillons (W5): quasi-stable time-periodic
              │   → orientation może istnieć ALE state nie jest static
              │   → wymaga rethink moduli space
              │
              └─ STRUCTURAL_NO_GO Path A
```

**Najpierwszy konkretny test (preceding numerical scan):**

> Czy effective Lagrangian dla δΦ wokół tła Φ̄, wynikający z M9.1''(Φ),
> zawiera **higher-derivative terms** (skala (∂δΦ)⁴ lub wyższa)?

**Source dla odpowiedzi:**
- L08 kink-fermion closure (audyt OPEN) — formal Φ-Lagrangian
- M9.1'' canonical metric definition — explicit form
- S07 M9.1'' derivation (audyt OPEN) — derivation z substratu

**Status:** OPEN — wymagana inspection M9.1''(Φ̄) struktury w istniejących
audytach.

## Konkretne literaturowe wnioski dla Phase 1

### L1: Sferyczny hedgehog wariant
W literaturze (Skyrme model, monopole), hedgehog ansätze mają postać:
```
φ_a(x) = f(r) · x_a / r        (spatial vector indices)
```
Dla **single real scalar** brak dim a → trivial sferyczny przypadek
φ(r) jest jedynym dostępnym sferycznym ansätzem.

**Implikacja:** w czystym single-Φ, jedyny sferyczny soliton jest
"radial breathing", BEZ orientation degree of freedom.

### L2: Anizotropowe ansätze
Możliwe ansätze:
```
Dipole:    φ(x) = f₁(r) + f₂(r)·cos(θ)
Quadrupole: φ(x) = f₁(r) + f₂(r)·(3cos²θ-1)
General:   φ(x) = ∑_{l,m} f_{lm}(r)·Y_{lm}(θ,φ)
```

Pytanie: czy któreś z tych są stable w Φ-EOM?

**Wynik literatury (general):** dla **standardowego** Lagrangianu, Derrick
wyklucza wszystkie. Dla non-standard (W6), możliwe że NON-spherical mają
NIŻSZĄ energię niż spherical (standard Skyrme: hedgehog jest lokalnym
minimum, ale w pewnych potencjałach NON-hedgehog są lower-energy).

### L3: Oscillon family
Oscillons w realistic potentials (e.g. V = ½m²Φ² - λΦ⁴ + g²Φ⁶) wykazują:
- Spherical oscillons (standard)
- Anizotropowe oscillons (mniej stabilne, krótsze życie)
- "Multi-bumps" (klastery oscillonów)

**Critical for TGP:** czy oscillon ma **orientation degree of freedom**?
- Spherical oscillon: oczywiście NIE
- Multi-bump: tak (orientation cluster axis)
- Single anisotropic oscillon: TAK (jeśli stable)

### L4: Time-dependent moduli spaces
Jeśli soliton jest oscillonem (time-periodic), moduli space jest
**bardziej skomplikowany**:
- Static moduli (spatial position, orientation jeśli istnieje)
- Phase moduli (offset czasowy oscylacji)
- Frequency modulus (jeśli rodzina o różnych ω)

**Implikacja dla SU(2):** może być konieczne rozszerzenie z SO(3) na
SO(3) × U(1) (orientation × phase) → quantum lift na SU(2) × U(1) lub
podobne. Wymaga dokładnej analizy w Phase 2.

## Re-assessment Phase 0 risk after literature review

### Aktualizacja R1
**R1 (orig):** Localized δΦ-soliton sferycznie symetryczny → spin 0.
**R1 (post-review):** **Dwa scenariusze threat-level:**
- (a) Standardowy Lagrangian → **Derrick**: brak static solitons w ogóle
  (nie tylko sferycznych)
- (b) Non-canonical Lagrangian (M9.1''-derived) → solitony mogą istnieć,
  Path A żywy

**Probability assessment update (subiektywna):**
- Pełen DERIVED: 10-15% (obniżone z 15-20%)
- STRUCTURAL CONDITIONAL: 30-40%
- STRUCTURAL_NO_GO: 40-50% (podniesione z 30-40%)
- EARLY_HALT: 15-20%

**Powód:** Derrick to konkretne, znane ograniczenie. Workaround wymaga
non-canonical structure z M9.1'' którego istnienie jest currently OPEN.

### Nowe NEEDS

#### N14: M9.1''(Φ) effective Lagrangian inspection
**Priority:** CRITICAL (BLOCKER)
**Phase:** 1
**Status:** OPEN — IMMEDIATE ACTION REQUIRED

**Description:** Wymagana inspection effective Lagrangian dla δΦ-fluctuations
wokół background Φ̄, wyprowadzonej z M9.1''(Φ) canonical metric:
- Czy zawiera **higher-derivative terms** ((∂δΦ)⁴, etc.)?
- Czy effective kinetic g^{μν}(Φ̄) jest non-trivial?
- Czy effective potential V_eff może być negatywny lokalnie?

**Resolution path:**
- Inspect istniejących plików S07_M911_derivation (audyt OPEN)
- Inspect L08_kink_fermion_closure
- Inspect [[../op-Phi-decomposition-photon-2026-05-07/Phase1_results.md]]
  — może już zawiera relevant expansion

**Bez tego:** Phase 1 numerical scan jest premature — może testujemy
ansätze które Derrick i tak wyklucza.

#### N15: Oscillon path eksploracja
**Priority:** IMPORTANT
**Phase:** 1
**Status:** OPEN

**Description:** Jeśli static soliton wykluczony, czy oscillon jest
viable alternative? Wymagana:
- Existence oscillonów w Φ-EOM (numerical scan)
- Orientation degree of freedom dla oscillonów (anisotropic shape lub
  cluster)
- Quantum quantization rotating/oscillating field configuration

## Action plan post-review

### Immediate (przed numerical scan)
1. Inspection M9.1''(Φ) effective Lagrangian (N14)
2. Decision: standard kinetic vs non-canonical?

### Conditional na N14
- **N14 RESOLVED z non-canonical TAK:**
  → kontynuujemy z static soliton ansätze (Phase 1 numerical scan)
- **N14 RESOLVED z standardowym kinetic:**
  → pivot na oscillon path (N15)
- **N14 niemożliwe do rozstrzygnięcia z istniejącej literatury:**
  → wymagana derivation effective Lagrangian de novo (significant work)

### Decision branching
Per CALIBRATION_PROTOCOL honest reporting, jeśli N14 → standard kinetic
i N15 → oscillon nie daje orientation, **STRUCTURAL_NO_GO** Path A
**niezależnie od ilości włożonej pracy**.

To jest sygnał metodologiczny: **honest reporting > sunk cost**.

## Cytaty key references (do dodania w bibliography)

- **Derrick GH (1964)** "Comments on Nonlinear Wave Equations as Models
  for Elementary Particles" J Math Phys 5:1252.
  → fundamentalne ograniczenie 3D scalar solitons

- **Skyrme TH (1962)** "A unified field theory of mesons and baryons"
  Nucl Phys 31:556.
  → pierwszy stable 3D soliton (multi-component + higher-derivative)

- **Coleman S (1985)** "Q-Balls" Nucl Phys B 262:263.
  → Q-ball construction (complex field)

- **Bogolyubsky IL, Makhankov VG (1976)** "Lifetime of pulsating solitons"
  JETP Lett 24:12.
  → first oscillon discovery

- **Gleiser M (1994)** "Pseudostable bubbles" Phys Rev D 49:2978.
  → modern oscillon analysis

- **Manton NS, Sutcliffe P (2004)** "Topological Solitons" Cambridge
  Monographs.
  → comprehensive textbook reference

## Honest reporting note

Niniejszy review zmienił landscape probability assessment Phase 0:
**STRUCTURAL_NO_GO probability wzrosło** w świetle Derricka.

To jest **honest result** — odkrycie konkretnego, znanego ograniczenia
matematycznego które trzeba addresować przed kontynuacją.

Nie dyskwalifikuje cyklu, ale **rekonfiguruje strategię**:
- N14 staje się PRIMARY GATE (przed numerical scan)
- Oscillon path staje się **alternative primary** (jeśli N14 → standard)
- Probability DERIVED jest niższe niż przed review

**To jest TGP-natywny sposób pracy:** honest assessment > optimism, każdy
nowy fakt może rekonfigurować strategię.

## Cross-references

- [[./README.md]] — overview cyklu
- [[./Phase0_balance.md]] — balance sheet (Phase 0 zaliczony)
- [[./NEEDS.md]] — needs list (N14, N15 do dodania)
- [[../op-Phi-decomposition-photon-2026-05-07/Phase1_results.md]] —
  potential source dla N14 (effective δΦ Lagrangian)
- [[../../audyt/S07_M911_derivation/]] — M9.1'' canonical metric definition
- [[../../audyt/L08_kink_fermion_closure/]] — kink-fermion model
