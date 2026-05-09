---
title: "op-SPIN-SU2-substrate-derivation — formal derivation SU(2) jako moduli space sklejeń z Φ-EOM"
date: 2026-05-08
type: research-cycle
status: CLOSED_STRUCTURAL_DERIVED
phase: Phase6_close
classification: STRUCTURAL_DERIVED
sympy_total: 47/47 PASS
parent: "[[../op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07/README.md]]"
related_audit:
  - "[[../../audyt/L08_kink_fermion_closure/]]"
  - "[[../../audyt/S05_tensor_sector_singleField/]]"
related_research:
  - "[[../op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07/]]"
  - "[[../op-Phi-decomposition-photon-2026-05-07/]]"
tgp_owner: research/op-SPIN-SU2-substrate-derivation-2026-05-08
tags:
  - research-cycle
  - spin
  - SU2
  - moduli-space
  - bell-entanglement
  - substrate-derivation
  - Phase0
  - phi-soliton
  - kink-fermion
---

# op-SPIN-SU2-substrate-derivation-2026-05-08

## Status

**OPEN — Phase 1 STRUCTURAL_EXPLORATION** (re-scoped 2026-05-08 post-N14).

Cykl badawczy zainicjowany 2026-05-08 jako **formalizacja sklejenie picture**
z working hypothesis [[../op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07/README.md]]
(sesja 2026-05-08, sekcja "Sklejenie picture").

**Re-scoping 2026-05-08 (post-N14):** N14 inspection wykazała że canonical
M9.1'' action NIE rozwiązuje Derricka, oraz że M9.1'' sama jest postulatem
(S07 P2 OPEN). Plus autor cyklu zaproponował **internal/external geometric
duality** ([[./Internal_external_geometry_proposal.md]]) jako alternative
framework — soliton ma "internal" geometrię która z external view daje
SU(2), ale wewnętrznie może NIE być SU(2). To rozbraja problem na poziomie
założeń.

**Status efektywny:** STRUCTURAL_EXPLORATION (nie pure derivation cycle).
Cel: identyfikacja konkretnej struktury matematycznej (Path-IE-A/B/C/D)
przed jakimkolwiek formal derivation.

## Geneza

Post-Stage-2 (`op-Phi-decomposition-photon-2026-05-07`, STRUCTURAL_NO_GO
dla photon-as-δΦ) konceptualna eksploracja wyłoniła schemat:

> Spin w TGP nie jest stanem "A albo B", ale **strukturą kompatybilności
> sklejenia przestrzeni**. A i B są wynikami sklejenia względem konkretnego
> warunku brzegowego. Ta sama geometria może dawać inne wyniki przy różnych
> osiach pomiaru.
>
> (autor cyklu, 2026-05-08)

Niniejszy cykl formalizuje ten schemat na konkretną strukturę matematyczną.

## Centralna hipoteza H1

**H1 (substrate derivation):** Moduli space localized δΦ-soliton (kink-fermion
z L08) ma topologię SO(3); kwantowa wave function żyje na podwójnym nakryciu
SU(2). Ten obiekt geometryczny realizuje **wszystkie podpisy spinora 1/2**
QM, a w szczególności reprodukuje:

1. Dwa wyniki pomiaru ±1
2. Strukturę inną niż wektor R³
3. 720°, nie 360°, symetrię
4. Projekcję cos²(θ/2)
5. Singletową korelację -cos(θ)
6. Niemożność tabeli predefined values (Kochen-Specker)

## Sześć wymagań — match do SU(2)

| # | Wymaganie | SU(2) realizacja | Status weryfikacji |
|---|-----------|------------------|--------------------|
| 1 | Dwa wyniki + / − | Reprezentacja fundamentalna 2D (spinor); σ_n eigenvalues ±1 | Standard QM ✓ |
| 2 | Nie wektor R³ | Element SU(2) = kwaternion jednostkowy = macierz unitarna 2×2 | Standard math ✓ |
| 3 | 720°, nie 360° | π_1(SO(3))=Z_2; SU(2) → SO(3) podwójne nakrycie; rotacja 360° → -ψ | Standard math ✓ |
| 4 | cos²(θ/2) projekcja | \|⟨↑\|ψ_θ⟩\|² = cos²(θ/2) — half-angle structure SU(2) | Standard QM ✓ |
| 5 | -cos(θ) Bell singlet | Antisymmetric (\|↑↓⟩-\|↓↑⟩)/√2; E(α,β) = -cos(α-β) | Standard QM ✓ |
| 6 | Brak globalnej tabeli | Kochen-Specker (1967): no consistent {0,1}-assignment dla SU(2)/dim≥3 | Standard math ✓ |

**Konkluzja weryfikacyjna:** SU(2) **trywialnie** spełnia wszystkie 6
wymagań — to są **podpisy spinora**. Nietrywialna jest część TGP:

> Czy SU(2) **emerguje z Φ-EOM** jako moduli space lokalizowanego soliton,
> czy jest dodatkową strukturą postulowaną?

## Trzy ścieżki realizacji — analiza preliminarna

### Ścieżka A: External orientation localized soliton

**Mechanizm:**
- Localized δΦ-soliton (kink-fermion z L08) **łamie sferyczną symetrię** —
  ma niezerowy multipole (dipole / quadrupole / hedgehog asymmetry)
- Orientation = element SO(3) (3 kąty Eulera)
- Quantum lift na podwójne nakrycie = SU(2)
- Single-valuedness wave function NA SU(2) (nie SO(3)) → spinor structure
- 720° emerguje **naturalnie** z geometrii nakrycia

**Szanse a priori:** 🟢 **WYSOKIE** — to jest standardowy argument Pauli/Dirac;
TGP wnosi **derivation orientation degree of freedom z Φ-EOM** (zamiast postulatu).

**Najsłabsze ogniwo:** czy localized δΦ-soliton **rzeczywiście ma
orientation degree of freedom**? Sferycznie symetryczny hedgehog miałby
trywialną orientation → spin 0. Wymaga pokazania że ground state NIE jest
sferycznie symetryczny (lub że zlamana sphericalism jest natywna).

**Status:** **PRIMARY PATH** dla niniejszego cyklu.

### Ścieżka B: Skyrmion-like internal winding

**Mechanizm:** δΦ-konfiguracja owija się wokół S²/S³ w polu wartości.

**Problem fundamentalny:** Φ jest jednowymiarowe (real scalar), brak
target manifold S²/S³ w którym może winding się wystąpić.

**Szanse a priori:** 🔴 **NISKIE** — wymagałoby rozszerzenia S05
(multi-component Φ), co jest CLOSED axiom (S05 PATH B).

**Status:** **OUT OF SCOPE** (re-open S05 wymagałby osobnego cyklu).

### Ścieżka C: Spacetime SO(3,1) embedding

**Mechanizm:** Pełna grupa Lorentza ma SU(2)_L × SU(2)_R w (j_L, j_R)
reprezentacji; kink-fermion siedzi w (1/2,0)⊕(0,1/2).

**Szanse a priori:** 🟡 **ŚREDNIE** — wymaga że solitonowa konfiguracja
wybiera subreprezentację Lorentza. To jest dokładnie pytanie L08.

**Status:** **CONTINGENT NA L08 CLOSURE** — śledzimy jako secondary, ale
nie zaczynamy.

## Plan cyklu Phase 0-6

### Phase 0: Balance sheet
**Status:** WIP (niniejszy plik + Phase0_balance.md + NEEDS.md).

Wykonanie: 8/8 gate criteria check, claims C1-Cn, falsifiability,
deal-breakers identification.

### Phase 1: Φ-EOM solitonowa konfiguracja
- Explicit Φ-EOM dla static localized configurations
- Existence proof / numerical scan: czy istnieją non-spherical solitons?
- Identification orientation degree of freedom
- Sympy / numerics verification

**Najsłabsze ogniwo:** N1, N2 (NEEDS) — wymaga Φ-EOM solver capability.

### Phase 2: Moduli space topologia
- Family rozwiązań Φ-EOM differing tylko orientation
- Computation moduli space topology
- Cel: pokazanie M ≅ SO(3) (lub equivalently S³/Z_2)

### Phase 3: Quantum lift na SU(2)
- Wave functional ψ[Φ] na konfiguracji
- Single-valuedness condition
- Pokazanie że spinor structure (SU(2)) jest wymagana, nie tylko dopuszczalna

### Phase 4: Born rule projekcja → cos²(θ/2)
- Pomiar = Φ-EOM coupling do external axis (B-field analog)
- Projekcja stanu na eigen-konfigurację względem osi
- Cel: derivation cos²(θ/2) z dynamiki Φ-EOM

### Phase 5: Singlet → Bell -cos(θ)
- Konfiguracja dwóch solitonów w stanie antysymetrycznym
- Wspólna geometria sklejenia (nielokalne)
- Cel: derivation E(α,β) = -cos(α-β)
- Plus: explicit no-signaling check (marginalne stats niezależne od second axis)

### Phase 6: ABSOLUTE BINDING gate
Kryteria zaliczenia per CALIBRATION_PROTOCOL:
- Czy SU(2) emerguje z Φ-EOM bez additional postulates?
- Czy 720° emerguje z geometrii (nie postulatu)?
- Czy Bell cos(θ) wynika z derivation?
- Czy no-signaling jest preservowane?
- Czy S05 single-Φ axiom zachowane?

**Klasyfikacja końcowa:** DERIVED / STRUCTURAL / STRUCTURAL_NO_GO.

## Probability assessment (subiektywna)

### Pre-cycle (initial)
| Outcome | Prob |
|---------|------|
| Pełen DERIVED | 15-20% |
| STRUCTURAL CONDITIONAL | 35-45% |
| STRUCTURAL_NO_GO | 30-40% |
| EARLY_HALT | 10-15% |

### Post-N14 (Derrick wins canonical M9.1'')
| Outcome | Prob |
|---------|------|
| Pełen DERIVED | 10-15% ↓ |
| STRUCTURAL CONDITIONAL | 30-40% |
| STRUCTURAL_NO_GO | **40-50%** ↑ |
| EARLY_HALT | 15-20% ↑ |

### Post-IE-proposal (Internal/external duality)
**Niepewne** — zależy od konkretnej Path-IE-A/B/C/D realizacji.

### Post-dynamic-equilibrium framework (2026-05-09)
| Outcome | Prob |
|---------|------|
| Pełen DERIVED | 25-35% (z 10-15%) |
| STRUCTURAL CONDITIONAL | 35-45% |
| STRUCTURAL_NO_GO | 20-30% (z 40-50%) |
| EARLY_HALT | 10-15% |

### Post-N17 (bifurcation analysis sympy 7/7 PASS)
| Outcome | Prob |
|---------|------|
| Pełen DERIVED | 30-40% ↑ |
| STRUCTURAL CONDITIONAL | 35-45% |
| STRUCTURAL_NO_GO | 15-25% ↓ |
| EARLY_HALT | 10-15% |

### Post-N18 (quantum lift sympy 7/7 PASS)
| Outcome | Prob |
|---------|------|
| Pełen DERIVED | 40-50% ↑ |
| STRUCTURAL CONDITIONAL | 30-40% |
| STRUCTURAL_NO_GO | 10-20% ↓ |
| EARLY_HALT | 5-10% ↓ |

### Post-Phase-1cd (spatial 3D extension, orientation OPEN)
| Outcome | Prob |
|---------|------|
| Pełen DERIVED | **30-40%** ↓ |
| STRUCTURAL CONDITIONAL | **35-45%** ↑ |
| STRUCTURAL_NO_GO | **20-30%** ↑ |
| EARLY_HALT | 5-10% |

**Reasoning Phase 1 cd:** spatial 3D analiza wykazała że orientation
degree of freedom **NIE JEST automatic consequence** N17. Musi być
zapewniona przez jeden z trzech mechanisms (A: SSB, B: dynamic temporal,
C: M9.1'' horizon deformations). Mechanism C jest **najsilniejszym
TGP-natywnym kandydatem** ale wymaga osobnej analizy. Honest acknowledgment
podnosi NO_GO probability i wprowadza CONDITIONAL na N18.

Per [[./Phase1_spatial3D_results.md]].

**Subiektywnie:** to jest najtrudniejszy cykl. Obecnie ma **konkretną
mathematical chain** (N17 → N18 → 5/6) ale z **konkretnym otwartym
pytaniem** (orientation mechanism). Cykl proceedss z **conditional**
flagą.

## Oczekiwane wartości każdego outcome

- **DERIVED:** TGP zyskuje natywną derivation spinora, Bell, splątania.
  Możliwy game-changer dla foundations.
- **STRUCTURAL CONDITIONAL:** ważny postęp, zarys formalizmu plus
  konkretne brakujące ogniwa do dalszych cykli.
- **STRUCTURAL_NO_GO:** także wartościowy — concretizes że TGP wymaga
  rozszerzenia (np. S05 reopen) lub że spin jest fundamentalnie irreducible
  do Φ-EOM. Honest negative result.
- **EARLY_HALT:** szybki zwrot z minimalnym zużyciem zasobów; lesson learned.

## Falsifiability — eksperymentalne testy

Niezależnie od formalizacji matematycznej, każdy claim cyklu musi mieć
empiryczne falsifiery (per CALIBRATION_PROTOCOL):

1. **Bell inequality test (CHSH):** standardowe eksperymenty (Aspect 1982,
   Hensen 2015 loophole-free) — TGP musi reprodukować S = 2√2.
2. **g_e ≈ 2:** anomalous magnetic moment electron (g-2 measurement)
3. **Stern-Gerlach:** dwie wiązki ±, projekcja cos²(θ/2)
4. **Spinor 720° (neutron interferometry):** Werner 1975, Rauch 2015 —
   neutron daje phase shift -1 po 360°, +1 po 720°
5. **Pauli exclusion:** atomic spectra, periodic table structure
6. **Singlet entanglement:** GHZ states, quantum teleportation experiments

**Critical:** TGP NIE może odbiegać od żadnego z powyższych w precyzji
~10⁻⁶ lub lepiej (eksperymentalne accuracy QM jest brutalnie wysokie).

## Cross-references

- [[../op-SPIN-MAG-leakage-lean-hypothesis-2026-05-07/README.md]] —
  parent working hypothesis (sklejenie picture, sesja 2026-05-08)
- [[../op-Phi-decomposition-photon-2026-05-07/Phase3_results.md]] —
  motivation: photon-as-δΦ STRUCTURAL_NO_GO sparked spin exploration
- [[../../audyt/L08_kink_fermion_closure/]] — **kluczowa zależność**:
  niniejszy cykl wymaga formal kink-fermion modelu z L08
- [[../../audyt/S05_tensor_sector_singleField/]] — single-Φ axiom
  preservation requirement (Path A musi zachować)
- [[../op-tensor-modes-Phi-FUTURE/]] — pokrewny placeholder dla
  geometric perturbation modes
- [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 0-6 protocol

## Kolejność sub-pliki

- `README.md` (niniejszy) — overview
- `Phase0_balance.md` — balance sheet, claims, gate criteria
- `NEEDS.md` — lista konkretnych potrzeb dla Phase 1-6
- (Phase1_*) — gdy Phase 0 zaliczony
