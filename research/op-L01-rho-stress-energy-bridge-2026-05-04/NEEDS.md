---
title: "NEEDS — L01 ρ stress-energy bridge (open problems)"
date: 2026-05-04
parent: "[[README.md]]"
type: needs
tgp_owner: research/op-L01-rho-stress-energy-bridge-2026-05-04
tags:
  - needs
  - L01
  - quantum-trace-anomaly
  - QCD-vacuum
  - SPARC-consistency
---

# NEEDS — L01 (otwarte problemy)

## Otwarte luki

### N1: Quantum trace anomaly EM (1-loop QED na curved background)

**Status:** OPEN.

**Problem:** Klasycznie `T^μ_μ_EM = 0`, ale w QFT 1-loop:

```
T^μ_μ_EM,1-loop = (β(α)/(2α)) · F_μν·F^μν + ... · R²
```

z `β(α) = α²/(3π)` (1-loop QED running of fine structure).

W typowych warunkach lab: efektywnie zero. W ekstremalnych (Schwinger
limit, magnetary): potencjalnie obserwowalne.

**Co potrzebne:**

1. Pełne 1-loop kowariantne renormalization w `g_eff = M9.1''`
2. Mapping na ψ.1.v2 dim-6 EFT basis (parity-even L₅'_a + parity-odd L₅'_b)
3. Verification, że quantum `ρ_EM` nie generuje 5-th force violation
   poza MICROSCOPE 10⁻¹⁵
4. Weryfikacja dla magnetar regimes (B ~ 10¹¹ T) — czy obserwowalne

**Source:** [[photon_treatment.md]] §5

**Kandydat dostawcy:** dedicated cycle `op-EM-trace-anomaly-TGP`
(estymata: ~3-4 tygodnie formal physics)

**Typ:** derivation (1-loop QFT na curved background)

### N2: QCD trace anomaly (gluon condensate)

**Status:** OPEN.

**Problem:** W QCD: `T^μ_μ_QCD ~ Λ_QCD⁴` (non-perturbative gluon
condensate from trace anomaly). To **DOES** sprzęga z Φ przez
`ρ_QCD ~ Λ_QCD⁴/c_0² ~ (300 MeV)⁴`.

**Konsekwencje:**

- W epoce QCD phase transition (T ~ 200 MeV w wczesnym Wszechświecie):
  Φ-evolution może być source-dominated by QCD vacuum
- W kosmologii: faza-przejście QCD wpływa na TGP w specyficznym epoce
- T-Λ closure: `ρ_vac,TGP = M_Pl²·H₀²/12` daje DE; dodatkowy QCD
  contribution może być nieuwzględniony

**Co potrzebne:**

1. Explicit formula `ρ_QCD(T)` jako funkcja temperatury (running of α_s)
2. Włączenie do FRW Friedmann equation
3. Sprawdzenie konsystencji z BBN (T ~ MeV) i CMB (T ~ eV)
4. Implikacje dla T-Λ closure

**Source:** [[SM_sector_mapping.md]] §4

**Kandydat dostawcy:** dedicated cycle `op-QCD-trace-anomaly-cosmology`
(estymata: ~4-6 tygodni, requires lattice QCD inputs)

**Typ:** derivation (non-perturbative QCD + cosmology)

### N3: SPARC fits consistency check

**Status:** OPEN — low priority cosmetic.

**Problem:** N-body symulacje TGP (SPARC galaxy rotation, `nbody/`)
używają `ρ = ρ_baryon` (HI + stars + bulge) bez explicit T^μ_μ derivation.

**Pytanie:** czy `ρ_baryon` w SPARC is the same as `ρ = -T^μ_μ/c_0²`
in non-relativistic limit?

**Analiza:**

Dla nierelatywistycznej materii (galactic stars + HI gas):

```
T^00 = ρ_rest · c² + (1/2)·ρ_rest·v²  ≈  ρ_rest·c²   (v² ≪ c²)
T^ii = ρ_rest·v²  ≈  0
T^μ_μ = T^00 - T^ii  ≈  -ρ_rest·c²
```

z minusem (signature):

```
ρ = -T^μ_μ/c² = ρ_rest    ✓
```

Czyli **klasycznie konsystentne**. SPARC fits powinny działać poprawnie.

**Co potrzebne:**

1. Explicit verification że `ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0²` z dokładnością <1%
2. Sanity check że nie ma double-counting `ρ_DM` (dark matter z TGP
   emergent vs separate component)
3. Update `nbody/` dokumentacji z explicit T^μ_μ → ρ mapping

**Source:** [[SM_sector_mapping.md]] §1, [[../../audyt/L01_rho_operational/NEEDS.md]] N4

**Kandydat dostawcy:** krótki re-derivation w SPARC manuscript (1-2 dni)

**Typ:** numerical (re-fit z explicit derivation)

### N4: Higgs sector explicit derivation

**Status:** OPEN — partial.

**Problem:** Dla Higgs scalar w stanie próżni, klasycznie `T^μ_μ = 0`
po SSB cancellation. Po kwantyzacji: trace anomaly dla skalarnego pola
+ fluktuacje Higgsa h(x).

**Co potrzebne:**

1. Explicit `T^μ_μ_Higgs` w 1-loop QFT
2. Mapping na `ρ_Higgs` w fluktuacjach h(x)
3. Konsekwencje dla TGP: czy Higgs vacuum (non-zero in Standard Model)
   sprzęga z Φ?

**Hipoteza:** Higgs vacuum TGP-decoupled w klasycznej granicy (T^μ_μ=0
exact), ale 1-loop daje `m_H²·h²` term sprzęgający fluktuacje.

**Source:** [[SM_sector_mapping.md]] §2

**Kandydat dostawcy:** część cyklu N1 (EM + Higgs trace anomaly together)

**Typ:** derivation (1-loop Higgs sektor)

### N5: Gauge field anomaly w SU(2) i U(1) (electroweak)

**Status:** OPEN — analogous do N1 i N2.

**Problem:** Yang-Mills SU(2)×U(1) ma analogiczne trace anomaly co
QCD (per non-Abelian sector). Czy electroweak sektor sprzęga z Φ
przez `ρ_EW`?

**Co potrzebne:**

1. β-function SU(2)×U(1) (slope: dla SU(2) b₀=−19/12, dla U(1) b₀=+41/12)
2. Trace anomaly contribution
3. Konsekwencje kosmologiczne (faza-przejście EW na T ~ 100 GeV)

**Source:** [[SM_sector_mapping.md]] §4 (analogiczne do QCD)

**Kandydat dostawcy:** część cyklu N2 (cosmological phase transitions
w TGP)

**Typ:** derivation (non-perturbative EW + cosmology)

## Pytania otwarte

- **Q1**: Czy quantum trace anomaly QED (β/(2α) F²) jest poprawnie
  zaadresowane przez ψ.1.v2 dim-6 EFT operatory, czy wymaga osobnego
  cyklu?

- **Q2**: Jak `ρ_QCD = Λ_QCD⁴` interaguje z T-Λ closure (ρ_vac,TGP =
  M_Pl²·H₀²/12)? Czy to są te same gęstości energii vacuum, czy oddzielne?

- **Q3**: Dla magnetarów (B ~ 10¹¹ T), oczekiwany observable Φ-shift
  od `ρ_EM_quantum` jest jaki rząd? — to jest *concrete falsifiability
  test* dla magnetar polar shift TT10 z τ.3 cyclu.

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | N1 (EM trace anomaly) | dedicated 1-loop QFT cycle | quantum field theory |
| B2 | N2 (QCD trace anomaly) | non-perturbative QCD inputs (lattice) + cosmology | T-Λ closure |
| B3 | N3 (SPARC consistency) | krótki re-derivation w nbody/ | low priority |
| B4 | N4 (Higgs sector) | część N1 cycle | depends on N1 |
| B5 | N5 (EW anomaly) | część N2 cycle | depends on N2 |

## Closed needs (po L01 cyklu 2026-05-04)

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| Formal kowariantna definicja ρ | ρ ≡ -T^μ_μ/c_0² (formal_definition.md) | L01 cykl | 2026-05-04 |
| Mapping SM sektorów | 5 sektorów (Dirac, scalar, EM, YM, fluid) | SM_sector_mapping.md | 2026-05-04 |
| Photon treatment (ρ_EM=0) | analiza 3 implikacji + 1 open (quantum) | photon_treatment.md | 2026-05-04 |
| Spójność z ax:metric-coupling | Option-2 audytu *DERIVED* nie *postulowane* | formal_definition.md §5 | 2026-05-04 |
| S04 N1 (formal kowariantna derywacja L_mat) | wyprowadzenie z ax:metric-coupling + perturbation | formal_definition.md §4 | 2026-05-04 |

## Cross-references

- [[README.md]] — cykl L01 indeks
- [[FINDINGS.md]] — eksportowalne wyniki
- [[../../audyt/L01_rho_operational/]]
- [[../../audyt/S04_metric_coupling_axiom/]] (powiązany — N1 formal derivation)
- [[../op-psi1-substrate-light-acceleration]] (ψ.1.v2 dim-6 EFT)
- [[../op-tau3-substrate-clock-acceleration]] (τ.3 magnetar polar)
- [[../op-omega1-substrate-em-coupling]] (ω.1 EM coupling)
- [[../closure_2026-04-26/Lambda_from_Phi0]] (T-Λ vacuum energy)
