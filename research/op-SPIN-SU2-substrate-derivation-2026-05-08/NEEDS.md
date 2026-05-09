---
title: "NEEDS — op-SPIN-SU2-substrate-derivation"
date: 2026-05-08
type: needs-list
status: WIP
parent: "[[./README.md]]"
tags:
  - needs
  - phase0
  - spin
  - SU2
---

# NEEDS — explicit wymagania dla cyklu

## Cel

Per CALIBRATION_PROTOCOL, Phase 0 wymaga explicit listy potrzeb. Każdy
NEED ma:
- **ID** (Nx)
- **Priority** (CRITICAL / IMPORTANT / NICE-TO-HAVE)
- **Phase** (kiedy jest potrzebny)
- **Status** (OPEN / RESOLVED / DEFERRED)
- **Description**
- **Resolution path** (jak adresować)

## CRITICAL needs (5)

### N1: Φ-EOM solver dla static localized solitonów
**Priority:** CRITICAL
**Phase:** 1
**Status:** OPEN

**Description:** Wymagany numerical solver dla Φ-EOM:
```
□Φ - q · F(Φ) = 0      (full nonlinear)
```
gdzie F(Φ) jest potential term (z M9.1''(Φ)). Dla static localized
configurations: ∇²Φ - q · F(Φ) = 0 z boundary conditions
Φ → Φ̄ at |x| → ∞.

**Resolution path:**
- (a) Implementacja relaxation method (Gauss-Seidel / multigrid) w Pythonie
- (b) Use of existing PDE solvers (FEniCS, deal.II)
- (c) Sympy dla simple symmetric ansätze (analytical)

**Bez tego:** Phase 1 niemożliwe — całość cyklu blokowana.

### N2: Existence non-spherical solitons w Φ-EOM
**Priority:** CRITICAL (deal-breaker risk)
**Phase:** 1
**Status:** OPEN

**Description:** Czy Φ-EOM dopuszcza **stable** lub **metastable**
localized solutions które są **NIEsferycznie symetryczne** (mają
niezerowy multipole — dipole, quadrupole, hedgehog asymmetry)?

Jeśli ground state jest sferycznie symetryczny → orientation degree of
freedom = 0 → cykl FAIL (R1).

**Resolution path:**
- N1 jako prerequisite
- Numerical scan po anizotropowych ansätze
- Stability analysis (linearize fluctuations wokół candidate solution)
- Comparison energy: spherical vs non-spherical

**Decision branch:**
- N2 RESOLVED z TAK → kontynuujemy Phase 1+
- N2 RESOLVED z NIE → cykl pivot lub STRUCTURAL_NO_GO

### N3: Topologiczna analiza moduli space
**Priority:** CRITICAL
**Phase:** 2
**Status:** OPEN

**Description:** Po znalezieniu rodziny rozwiązań (N2 RESOLVED), trzeba
matematycznie scharakteryzować moduli space M:
- Wymiar (oczekiwany: 3 dla spatial orientation)
- Topologia (oczekiwana: M ≅ SO(3) lub S³/Z_2)
- Czy jest gładką varietas

**Resolution path:**
- Identification group of symmetries acting on solutions
- Quotient: M = (configuration space) / (gauge equivalence)
- Computation π_1(M), π_0(M)

### N4: Quantum lift na podwójne nakrycie
**Priority:** CRITICAL
**Phase:** 3
**Status:** OPEN

**Description:** Wave functional ψ[Φ] na konfiguracji moduli space.
Czy single-valuedness wymaga lift na SU(2) (nie SO(3))?

**Resolution path:**
- Standard argument Pauli/Dirac
- Adaptacja do TGP context — czy działa bez additional assumptions?
- Explicit construction wave functional na covering space

### N5: Φ-EOM coupling do external axis (B-field analog)
**Priority:** CRITICAL
**Phase:** 4
**Status:** OPEN

**Description:** Pomiar = Φ-EOM coupling do external axis. Wymagana:
- Definicja "external axis" w czystym Φ language (gradient Φ̄? perturbation
  pattern?)
- Mechanizm projekcji: jak external axis selects eigen-konfiguracja?
- Derivation cos²(θ/2) z dynamiki

**Resolution path:**
- Stage 2 wynik: photon = standard A_μ; więc B-field = standard EM
- Coupling δΦ-soliton ↔ A_μ (jeśli soliton ma efektywny "ładunek
  geometryczny")
- LUB pure Φ analog: lokalny gradient Φ̄ jako "axis"

## IMPORTANT needs (4)

### N6: Singletu konstrukcja w sklejenie picture
**Priority:** IMPORTANT
**Phase:** 5
**Status:** ✅ **RESOLVED STRUCTURAL DERIVED** (2026-05-09, sympy 8/8 PASS)

**Description:** Jak skonstruować "wspólną geometrię sklejenia" dla 2
solitonów na odległość?

**Resolution path:**
- Tensor product moduli spaces M_1 × M_2
- Antysymetryzacja → singletu state
- Explicit derivation E(α,β) = -cos(α-β)

**Resolution:** [[./Phase5_singletu_Bell_results.md]] + [[./Phase5_singletu_Bell_sympy.py]]
- |Ψ⟩ = (1/√2)(|+−⟩ - |−+⟩) ✓
- Antisymmetric P₁₂|Ψ⟩ = -|Ψ⟩ ✓
- J²|Ψ⟩ = 0 (singletu = j=0) ✓
- E(α,β) = -cos(α-β) ✓
- CHSH = 2√2 (Tsirelson bound) ✓

**Caveat:** TGP-natywne fizyczne mapping pair production wymaga osobnej analizy.

### N7: No-signaling formal proof
**Priority:** IMPORTANT
**Phase:** 5
**Status:** ✅ **RESOLVED DERIVED** (2026-05-09, sympy w Phase 5)

**Description:** Marginalne stats P(+|α) niezależne od wyboru osi β
detektora 2 (i symetrycznie).

**Resolution path:**
- Trace out subsystem 2
- Verify P(+|α) wynik niezależnie od czego się dzieje na 2

**Resolution:** [[./Phase5_singletu_Bell_results.md]]
- P(+|a) = 1/2 (rotation invariant) ✓
- Reduced ρ₁ = Tr₂(|Ψ⟩⟨Ψ|) = (1/2)I (maximally mixed) ✓
- Trace formal property → no-signaling automatic z framework

### N8: g_e ≈ 2 derivation w niniejszym frameworku
**Priority:** IMPORTANT
**Phase:** 4
**Status:** OPEN (linked to working hypothesis blinking section)

**Description:** Czy SU(2) framework reprodukuje g_e ≈ 2 z working
hypothesis blinking mechanism (static + dynamic = 2)?

**Resolution path:**
- Magnetic moment μ = ⟨ψ|σ|ψ⟩ z Φ-EOM coupling
- Ratio do orbital μ → g factor
- Cel: g = 2 w lowest order, deviation g-2 z radiative corrections (Φ̄
  fluctuations)

### N9: Stern-Gerlach mechanism
**Priority:** IMPORTANT
**Phase:** 4-5
**Status:** OPEN

**Description:** Pełen mechanizm Stern-Gerlach: niejednorodne B-field
splita wiązkę na dwie. W TGP-natywnym ujęciu — czy to wynika z dynamiki
soliton-Φ̄?

**Resolution path:**
- Force on soliton z gradient B
- Force depends on lean (sklejenie) state względem osi B
- Two trajectories → SG split

## CRITICAL needs added post-literature-review (N14, N15)

### N14: M9.1''(Φ) effective Lagrangian — non-canonical structure?
**Priority:** CRITICAL — **BLOCKER**
**Phase:** 1 (immediate, przed numerical scan)
**Status:** OPEN
**Added:** 2026-05-08 post-literature-review

**Description:** Derrick's theorem (1964) wyklucza static localized
solitony w 3D dla standardowego Lagrangianu z V≥0. Path A wymaga
non-canonical kinetic structure (Skyrme-like) lub oscillon mechanism.

**Konkretne pytanie:** czy effective Lagrangian dla δΦ-fluctuations
wokół Φ̄ wyprowadzony z M9.1''(Φ) canonical metric zawiera:
- Higher-derivative terms ((∂δΦ)⁴ lub wyższe)?
- Non-trivial Φ̄-dependent metric g^{μν}(Φ̄) inducing non-canonical
  kinetic strength?
- Effective potential V_eff z negative regions lokalnie?

**Resolution path:**
- Inspection [[../op-Phi-decomposition-photon-2026-05-07/Phase1_results.md]]
  (zawiera Φ̄+δΦ expansion)
- Inspection [[../../audyt/S07_M911_derivation/]]
  (M9.1'' canonical metric)
- Możliwa derivation de novo effective Lagrangianu z pełnej Φ-action
  z M9.1''(Φ)

**Decision:**
- N14 → non-canonical (TAK) → Phase 1 z static soliton ansätze
- N14 → standard kinetic → pivot na N15 oscillon path
- N14 → niemożliwe rozstrzygnąć z istniejących źródeł → wymagana
  derivation de novo

### N15: Oscillon path — alternative primary jeśli N14 fails
**Priority:** CRITICAL (jeśli N14 → standard kinetic)
**Phase:** 1
**Status:** OPEN

**Description:** Jeśli static solitony wykluczone Derrickiem, oscillons
(time-periodic quasi-stable) są naturalną alternatywą:
- Existence oscillonów w Φ-EOM dla M9.1''-derived V(Φ)
- Orientation degree of freedom dla oscillonów (anisotropic shape /
  multi-bump cluster)
- Quantization framework dla time-dependent solitons
- Czy moduli space oscillonu zawiera SO(3)?

**Resolution path:**
- Numerical scan oscillonów dla generic V(Φ) (Gleiser-style)
- Specific dla TGP: V(Φ) z M9.1''(Φ)
- Stability analysis (lifetime estimate)
- Orientation analysis (anisotropic vs spherical)

**Caveat:** oscillons NIE są strictly stable — emit radiation. Question
czy lifetime jest **kosmologicznie wystarczający** dla representing
stable elementary particles (e.g. electron lifetime > 10²⁹ years).

## CRITICAL needs added post-Phase-1cd (N20, N21)

### N20: Orientation mechanism — numerical PDE solver
**Priority:** CRITICAL (deal-breaker risk)
**Phase:** deferred to future cycle
**Status:** OPEN

**Description:** Pełen numerical solver dla 1D radial Φ-EOM:
- Test sferyczny vs anizotropowy ground state
- Stability multipole modes (l=0, 1, 2, ...)
- Identify which orientation mechanism (A/B/C) active

**Scope:** osobny cykl numerical PDE
**Bez tego:** N18 SU(2) lift pozostaje conditional

### N21: Mechanizm C — M9.1'' horizon deformations
**Priority:** IMPORTANT
**Phase:** 2 (proceed pomimo N20 OPEN)
**Status:** ✅ **RESOLVED PARTIAL DERIVED** (2026-05-09, sympy 11/11 PASS)

**Description:** Analytical analysis horyzontu ψ=4/3 jako boundary
configuration:
- 2-sphere topology
- Multipole expansion deformations
- Connection do quantum SU(2) lift

**Plus:** najmocniejsza TGP-natywna ścieżka, niezależna od full PDE solver.

**Resolution:** [[./Phase2_N21_results.md]] + [[./Phase2_N21_horizon_deformations_sympy.py]]
- Horizon condition ψ=4/3 → g_tt=0 ✓
- 2-sphere topology with SO(3) isometry ✓
- l=1 dipole = 3-vector orientation degree of freedom ✓
- SO(3)/Z₂ → SU(2) double cover (720° symmetry) ✓
- Centrifugal stability l=1 mode (soft Goldstone-like) ✓

**Convergence z N17/N18:** dwie niezależne ścieżki do SU(2), wzmocniona confidence.

**Limitations (open):** preferencja jednego l=1 mode wymaga external trigger (consistent z QM measurement). Backreaction na M9.1'' background = future work.

## NICE-TO-HAVE needs (4)

### N10: Connection do L08 kink-fermion closure
**Priority:** NICE-TO-HAVE
**Phase:** 1-2
**Status:** OPEN (L08 audyt OPEN)

**Description:** Czy formal kink-fermion z L08 jest dokładnie obiektem
naszego niniejszego cyklu? Jeśli tak, L08 closure i niniejszy cykl
mogą się **wzajemnie zasilać**.

### N11: Connection do tensor modes (op-tensor-modes-Phi-FUTURE)
**Priority:** NICE-TO-HAVE
**Phase:** N/A
**Status:** placeholder cycle existing

**Description:** Tensor modes = geometric perturbations Φ̄. Czy soliton
spin-1/2 ma natywne sprzężenie z tensor modes (analog gravitational
spin-orbit)?

### N12: Kochen-Specker explicit demonstration w SU(2)
**Priority:** NICE-TO-HAVE
**Phase:** 5
**Status:** OPEN

**Description:** Pokazanie że SU(2) state w niniejszym frameworku NIE
ma global hidden variable assignment dla wszystkich osi (zgodnie z
KS theorem).

**Resolution path:** standardowa konstrukcja KS (33 vector, etc.) +
adaptacja do sklejenie language.

### N13: Phenomenological cross-check g-2 anomaly
**Priority:** NICE-TO-HAVE
**Phase:** 6
**Status:** OPEN (extension)

**Description:** Anomaly g-2 ≈ 0.00231930... — czy TGP framework daje
przewidywanie zgodne z QED + new physics tension (~5σ Fermilab/BNL)?

**Resolution path:** computation z Φ̄ fluctuation contributions; może
być **TGP-specific prediction** odrzucalny eksperymentalnie.

## CRITICAL needs added post-dynamic-equilibrium (N16-N19)

Dodane po insightcie autora że soliton intrinsically meta-stable
(zanik vs ekspansja bifurcation), stabilizowany przez interakcję z Φ̄.
Per [[./Dynamic_equilibrium_framework.md]].

### N16: Dynamic equilibrium soliton-background formalization
**Priority:** CRITICAL
**Phase:** 1
**Status:** ✅ **RESOLVED STRUCTURAL DERIVED** (2026-05-09, sympy 3/3 PASS)

**Description:** Energy budget E_sol + E_int + E_rad. Equilibrium
condition. Stability under perturbation. Scaling analysis of E_int
under soliton rescaling (klucz dla Derrick rozbrojenia).

**Resolution:** [[./Phase1_N16_results.md]] + [[./Phase1_N16_dynamic_equilibrium_sympy.py]]
- Standard Derrick (T/λ + U/λ³) confirmed FAILS dla static stable ✓
- Z background coupling E_int = S λ² (skin/boundary):
  - Equilibrium 2S λ⁵ - T λ² - 3U = 0 ma real positive root ✓
  - λ* ≈ 1.169 dla T=U=S=1 ✓
  - d²E/dλ² > 0 STABLE local minimum ✓
- Parameter scan S∈{0.1,...,5}: smooth λ*(S) ✓
- Connection do N17 bifurcation: S→0 limit recovers 2-outcome dynamics ✓

**Key insight:** Derrick rozbrojony **strukturalnie** (TGP nie szuka static
stable solutions, szuka dynamic equilibria z background). E_int z α>0 scaling
dostarcza brakujący term w Derrick scaling argument.

**Limitations:** specific coupling mechanism (α=2), numerical T,U,S dla TGP
wymagają solitonowego ansatzu (out of scope analytical), backreaction = future.

### N17: Bifurcation analysis (zanik vs ekspansja)
**Priority:** CRITICAL
**Phase:** 1-2
**Status:** OPEN

**Description:** Confirm że izolowany soliton ma exactly 2 dynamic
branches (decay vs growth). Phase plane analysis Φ-EOM. Show 2-state
structure jest fundamental dla niniejszej hipotezy.

### N18: Map bifurcation → SU(2) fundamental rep
**Priority:** CRITICAL
**Phase:** 2
**Status:** OPEN

**Description:** Pokazanie że 2-branch state space ma SU(2) strukturę:
continuous superposition |α|²+|β|²=1, SU(2) action jako bifurcation
operators, Bloch sphere geometry, half-angle structure.

### N19: External SO(3) embedding → induced SU(2)
**Priority:** IMPORTANT
**Phase:** 3
**Status:** ✅ **RESOLVED STRUCTURAL DERIVED** (2026-05-09, sympy 11/11 PASS)

**Description:** External R³ rotation indukuje SU(2) transformation
na bifurcation space. "Lean direction" w 3D = (θ,φ) na Bloch sphere.
Spinor structure jako naturalna konsekwencja embeddingu.

**Resolution:** [[./Phase3_N19_results.md]] + [[./Phase3_N19_external_embedding_sympy.py]]
- |+,n̂⟩ = (cos(θ/2), sin(θ/2)e^(iφ)) eigenstate σ·n̂ ✓
- U(α,m̂) = exp(-iα/2 σ·m̂) induced from external R(α,m̂) ✓
- U†U = I, det(U) = 1 (SU(2)) ✓
- Group homomorphism U(α₁+α₂) = U(α₁)U(α₂) ✓
- Double cover: U(2π)=-I, U(4π)=+I ✓

**Trzecia niezależna ścieżka do SU(2)** (z N17/N18 i N21).
Spinor structure NIE jest postulatem - emerguje z geometric embedding.

## Summary by phase

| Phase | Critical needs | Important | Nice |
|-------|---------------|-----------|------|
| 1 | **N16, N17 (dynamic equilibrium)**, N14 (RESOLVED neg), N1, N2, N15 | — | N10 |
| 2 | **N18** (bifurc → SU(2)), N3 | — | — |
| 3 | N4, **N19** (external embedding) | — | — |
| 4 | N5 | N8, N9 | — |
| 5 | — | N6, N7 | N12 |
| 6 | — | — | N13, N11 |

**Status update post-dynamic-equilibrium framework (2026-05-09):**
- N14 BLOCKER **RESOLVED z negative** dla canonical M9.1'' static path
- **PIVOT** na dynamic equilibrium framework (N16-N19) jako primary path
- Dynamic equilibrium **rozbraja Derricka** strukturalnie (nie technically),
  dopuszczając wartościowe progress przy mniejszych assumptions o M9.1''

## Decision criteria — Phase 0 → Phase 1

**Wymagane przed startem Phase 1:**
- ☑ Phase 0 balance sheet completed (gate criteria 8/8)
- ☑ NEEDS document complete (niniejszy plik)
- ☐ N1 setup plan (Φ-EOM solver) — może być rozwiązane przy starcie Phase 1
- ☐ User confirmation Path A primary

**Critical decision points w cyklu:**

1. **Po N2** (Phase 1): czy non-spherical solitons istnieją?
   - TAK → Phase 2
   - NIE → cycle pivot (rozważyć Path C lub STRUCTURAL_NO_GO)

2. **Po N3** (Phase 2): czy moduli space ≅ SO(3)?
   - TAK → Phase 3
   - INNE topology (np. SO(3) × inne) → analiza co to oznacza
   - Trywialna → STRUCTURAL_NO_GO

3. **Po N6, N7** (Phase 5): czy Bell -cos(θ) i no-signaling?
   - OBA → DERIVED (best case)
   - JEDEN → STRUCTURAL CONDITIONAL
   - ŻADEN → STRUCTURAL_NO_GO

## Honest reporting note

Lista NEEDS jest **wstępna** (Phase 0). Może rosnąć w Phase 1+ jak
ujawnia się więcej technicznych wymagań. Każdy nowy NEED dodawany
incrementalnie z timestamp i justification.

**Status NEEDS document:** WIP, początkowa wersja 2026-05-08.
