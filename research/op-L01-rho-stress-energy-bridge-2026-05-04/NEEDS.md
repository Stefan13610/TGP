---
title: "NEEDS — L01 ρ stress-energy bridge (open problems)"
date: 2026-05-04
last_addendum: 2026-05-10
parent: "[[README.md]]"
type: needs
tgp_owner: research/op-L01-rho-stress-energy-bridge-2026-05-04
tags:
  - needs
  - L01
  - quantum-trace-anomaly
  - QCD-vacuum
  - SPARC-consistency
  - native-observables-first
---

# NEEDS — L01 (otwarte problemy)

## §0 — Status methodology (post-2026-05-10)

> **Native-observables-first overlay binding 2026-05-10+:**
> [[ADDENDUM_2026-05-10_native_observables_first.md]] aplikuje
> [[../../meta/PPN_AS_PROJECTION.md]] methodology do tego cyklu. ρ jest **natywnym
> source field** dla Φ-EOM; PPN/ppE/cosmology są L2 projekcjami. NEEDS N1–N5 są
> *deferred ρ-corrections*; trójwarstwowa specyfikacja (L1 native obserwable / L2
> projection chart / L3 falsifikator) per NEEDS jest w sekcji [§T — Three-layer
> specification per NEEDS](#T) niżej. Treść otwartych luk poniżej (N1–N5, Q1–Q3,
> blokery) jest **nieruszona** — addendum jedynie dodaje warstwę interpretacyjną.

## Otwarte luki

### N1: Quantum trace anomaly EM (1-loop QED na curved background)

**Status:** **CLOSED** (2026-05-11) przez dedicated cycle
[[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] (STRUCTURAL_DERIVED, 16/16 sympy PASS,
6/6 P-requirements RESOLVED).

**Problem:** Klasycznie `T^μ_μ_EM = 0`, ale w QFT 1-loop:

```
T^μ_μ_EM,1-loop = (β(α)/(2α)) · F_μν·F^μν + ... · R²
```

z `β(α) = (2/(3π))·α²` (1-loop QED, single Dirac fermion; Capper-Duff-Halpern 1974).

**Status post-2026-05-11 — closure summary:**

- **Trace anomaly prefactor LOCK** (Phase 1 sympy T2): `β(α)/(2α) = α/(3π) ≈ 7.74·10⁻⁴`
  (sympy LOCK; cytowana w §3.2 Q3 wartość `7.7·10⁻⁷` jest **typo** — correct = 7.74·10⁻⁴).
- **Renormalized ρ_EM_quantum[{Φ_i}]** explicit form (Phase 1 §1.3 + Phase 2 §1.4):
  `ρ_EM_quantum = -[(α/(3π))·F² + γ_i·(curvature × F²) + Riegert local σ_eff = function(ψ)] / c_0²`
- **g_eff = G[{Φ_i}]** (NIE M9.1''; M9.1'' FALSIFIED 5σ GWTC-3 2026-05-09).
- **Theorem 2.1 (Disjointness)** (Phase 2 §2 + sympy T4): trace anomaly TGP-reduced
  operator classes są strukturalnie DISJOINT od ψ.1.v3 canonical basis B = {L₅'_a, L₅'_b}.
  Q1 *konstruktywnie* potwierdzony.
- **Lab regime** (B=1 T): ρ_EM_quantum ~ 7·10⁻¹⁵ kg/m³ inside magnet; niewykrywalne.
- **Magnetar typical** (B=10¹⁰ T): ρ_EM_quantum/ρ_NS ~ 1.7·10⁻¹².
- **Magnetar extreme** (B=10¹¹ T): ratio ~1.7·10⁻¹⁰ z corrected α/(3π) (NIE 10⁻¹²).
- **MICROSCOPE η bound:** `η_TGP_EM_quantum (Pt vs Ti) = 0` strukturalnie z universal
  coupling structure (S05) + B9 baseline 1.32·10⁻²⁶ preserved.
- **GW170817:** `Δc/c ~ 10⁻⁸⁰` ≪ 9·10⁻²² (~58 OOM margin).
- **R6 (Asorey-2015 QEP violations) closed strukturalnie** — TGP universal coupling immune.
- **R5 (perturbative breakdown B ≳ B_QED ≈ 4.4·10⁹ T)** honestly documented; deferred.

Patrz [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase_FINAL_close.md]].

**Source:** [[photon_treatment.md]] §5

**Closure:** dedicated cycle `op-L01-N1-EM-trace-anomaly-TGP-2026-05-11` (STRUCTURAL_DERIVED).

**Typ:** derivation (1-loop QFT na curved background)

### N2: QCD trace anomaly (gluon condensate)

**Status:** **CLOSED** (2026-05-11) przez dedicated cycle
[[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/]] (STRUCTURAL_DERIVED,
24/24 sympy PASS, 6/6 P-requirements RESOLVED).

**Problem:** W QCD: `T^μ_μ_QCD ~ Λ_QCD⁴` (non-perturbative gluon
condensate from trace anomaly). To **DOES** sprzęga z Φ przez
`ρ_QCD ~ Λ_QCD⁴/c_0² ~ (217 MeV)⁴`.

**Status post-2026-05-11 — closure summary:**

- **β-function 1-loop QCD LOCK:** `β_QCD(g) = -(b_0/(16π²))·g³` z `b_0 = (11/3)N_c - (2/3)N_f`;
  N_c=3, N_f=6 daje **b_0 = 7** (asymptotic freedom, opposite to QED).
- **Trace anomaly explicit:** `T^μ_μ_QCD,1-loop = (β(g)/(2g))·G² = -(b_0 α_s)/(8π)·G²`.
- **Gluon condensate:** `⟨α_s G²/π⟩_0 ≈ 0.012 GeV⁴` (SVZ-1979 + lattice 2018+ external);
  ρ_QCD_vacuum equivalent ~ 2.8·10¹⁸ kg/m³.
- **Decomposition:** ρ_QCD(T) = ρ_QCD_vacuum (constant SVZ) + ρ_QCD_thermal(T)
  (transient z Δ(T) interaction measure z lattice EoS, peak Δ_max ≈ 4 near
  T_c = 156 ± 9 MeV crossover).
- **Friedmann modyfikacja:** ρ_QCD_anomaly/ε_total ≈ 26-35% w narrow T-window wokół T_c
  (FWHM ~50 MeV); t_H(T_c) ~10⁻⁵ s consistent z standard cosmology QCD epoch.
- **IR limit (T<<Λ_QCD):** Δ(T) → 0 strukturalnie; today's ρ_QCD ≡ ρ_QCD_vacuum.
- **Q2 F1 KONSTRUKTYWNIE VERIFIED dla QCD sektora:** matter-vacuum decoupling
  potwierdzone z dedicated derivation; T-Λ ratio empirical 1.020 ± 0.02 preserved.
- **BBN ⁴He Y_p:** TGP = ΛCDM = 0.247 ± 0.001 vs PDG 2024 obs 0.245 ± 0.003 — within **0.55σ**.
- **D/H:** TGP = ΛCDM = 2.5·10⁻⁵ vs PDG 2024 obs 2.527·10⁻⁵ — within **0.26σ**.
- **CMB Planck 2018:** ω_b=0.02237, ω_m=0.1430, Ω_Λ=0.6889 wszystkie preserved automatic.
- **PTA NANOGrav 15-yr:** crossover (NIE first-order) → no dominant QCD signal →
  consistent z SMBHB consensus origin. **R6 closed.**
- **R1-R7 risks:** 6 fully closed + R2 honestly documented (lattice external).

Patrz [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase_FINAL_close.md]].

**Closure:** dedicated cycle `op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11`
(STRUCTURAL_DERIVED).

**Open (deferred precision, NOT blockers):**
- Multi-loop β_QCD precision (b_1, b_2, b_3) + sub-percent lattice EoS — extension cycle
- N4 Higgs sector extension (separate cycle)
- N5 EW gauge anomaly extension (separate cycle)

**Typ:** derivation (non-perturbative QCD + cosmology integration)

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

**Status:** **CLOSED** (2026-05-11) przez dedicated compact cycle
[[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] (STRUCTURAL_DERIVED, 8/8 sympy
PASS, 6/6 P-requirements RESOLVED). **ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0²
verified do <10⁻⁶ precision** (target 1%, achieved 6 OOM below). R1 (double-counting)
closed strukturalnie; R2 honestly documented (galactic-disk regime scope).

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

**Status:** ✅ **CLOSED 2026-05-11** — STRUCTURAL_DERIVED via dedicated cycle.

**Closure source:** [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/Phase_FINAL_close.md]]
(24/24 sympy PASS, 6/6 P-requirements RESOLVED, 5/6 risks closed + 1 deferred R3).

**Resolution (post-N4 2026-05-11):**

1. **Classical SSB cancellation explicit:** bare T^μ_μ_vac_bare = -λv⁴ ≠ 0;
   post-renormalization T_vac_renorm = 0 (Q2 F1 substrate-decoupling enforced)
2. **1-loop trace anomaly form:** `T^μ_μ_quantum = (1+γ_m)m_H²h² + (β_λ/4)h⁴
   + curvature·h² + Riegert local`
3. **β_λ ≈ -0.033** (top Yukawa dominant negative); **γ_m ≈ -0.027 ÷ -0.035**
4. **Q2 F1 konstruktywnie verified dla Higgs sektora**: OOM gap 55.3 (bare
   ρ_Higgs_vac vs Λ_obs) substrate-decoupled per single-Φ + substrate-vacuum
   identification (sympy T5 Phase 2)
5. **R5 LOCK**: m_H=125.25 GeV > endpoint 80 GeV (KLRS 1996, DRR 2014) ⇒
   crossover (NOT first-order); LISA Ω_GW^EW = 0 strukturalnie (falsifiable post-2035)
6. **Hipoteza zweryfikowana**: Higgs vacuum TGP-decoupled (Q2 F1) jak
   hipotezowano; 1-loop daje m_H²·h² term ze γ_m anomalous dim; cosmological
   transient source ρ_Higgs_thermal(T_EW) ≈ 0.83% radiation already in g_*(T_EW)

**Empirical commitments (post-N4 falsifiable):**
- LHC m_H stability (Run 3, HL-LHC, FCC-ee)
- LISA Ω_GW^EW = 0 (R5 LOCK, post-2035, double-falsification jeśli detected)
- HL-LHC λ_HHH null test (Δλ/λ_SM ≈ 0, ±50% precision 2030+)
- FCC-ee λ_HHH null test (Δλ/λ_SM ≈ 0, ±5% precision 2045+)

**Deferred R3 (hierarchy problem):** Q2 F1 + S05 *strengthens consistency*
z m_H stability, ale full theoretical resolution byłaby revolutionary
breakthrough → outside cycle scope.

**Typ:** derivation (1-loop Higgs sektor + EW cosmology + phenomenology bounds)

### N5: Gauge field anomaly w SU(2) i U(1) (electroweak)

**Status:** ✅ **CLOSED 2026-05-11** — STRUCTURAL_DERIVED via dedicated cycle
(LAST L01 N-need closed; full SM matter+gauge coverage milestone).

**Closure source:** [[../op-L01-N5-EW-gauge-anomaly-2026-05-11/Phase_FINAL_close.md]]
(8/8 sympy PASS, 6/6 P-requirements RESOLVED, 6/6 risks closed konstruktywnie/inherited).

**Resolution (post-N5 2026-05-11):**

1. **β_SU(2) asymptotic freedom**: b₀=19/6 > 0; β(M_Z) ≈ -5.56·10⁻³
2. **β_U(1) Landau pole**: b₀=41/6 > 0; β(M_Z) ≈ +1.97·10⁻³; μ_LP ~10⁴² GeV >> M_Pl
3. **Trace anomaly form** explicit dla obu: T^μ_μ = (β/2g)·F² (Adler-Collins-Duncan 1977)
4. **Q2 F1 konstruktywnie verified dla gauge sektora**: OOM gap 54.9 (bare ρ_gauge_vac
   ≈ 1.94·10⁴⁴ eV⁴ vs Λ_obs); analog N4 Higgs pattern
5. **EW cosmology N4-inheritance**: R5 LOCK crossover; gauge bosons freeze-out
   T~3.2 GeV >> T_BBN; LISA Ω_GW^EW=0 inherited
6. **PDG 2024 EW precision** preserved exactly (M_W, M_Z, sin²θ_W tree+loop)
7. **Cross-cycle 5-fold SM coverage** post-N5: EM + QCD + SPARC + Higgs + EW gauge

**LAST L01 N-NEED CLOSED — full SM matter+gauge sektor coverage milestone 2026-05-11.**

**Typ:** derivation (compact single-session z silnej architecture inheritance N1+N2+N4)

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

  **Status post-2026-05-10:** **CLOSED** z native-first operator class analysis.
  - ψ.1.v3 canonical basis (Phase 7 Hilbert-series-style enumeration):
    `B_ψ.1.v3^dim-6 = {L₅'_a, L₅'_b}` — wszystkie operatory mają **≥1 ∂lnX leg**
    (cross-sector class)
  - Quantum trace anomaly EM `(α/(3π))·F²` jest **pure-photon** dim-4 operator
    (renormalizacja gauge coupling), **bez** ∂lnX leg
  - Operator classes są **disjoint** — pure-photon vs cross-sector (Phase 7 T7.1
    invariance filter explicite filtruje pure-photon sektor)
  - **Wniosek:** ψ.1.v3 dim-6 EFT NIE pokrywa quantum trace anomaly. **N1 dedicated
    cycle** [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] zrealizowany 2026-05-11
    (STRUCTURAL_DERIVED) — Theorem 2.1 (Disjointness) potwierdzona *konstruktywnie*
    z dedicated 1-loop QED derivation w obecności emergent-metric `g_eff[{Φ_i}]`.

  Patrz [[../op-psi1-substrate-light-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §3
  + [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase2_results.md]] §2.

- **Q2**: Jak `ρ_QCD = Λ_QCD⁴` interaguje z T-Λ closure (ρ_vac,TGP =
  M_Pl²·H₀²/12)? Czy to są te same gęstości energii vacuum, czy oddzielne?

  **Status post-2026-05-10:** **CLOSED** przez dedicated synthesis cycle
  [[../op-Q2-vacuum-budget-2026-05-10]] (STRUCTURAL DERIVED).
  - Native answer: **OSOBNE.** SM matter sector vacua (`ρ_QCD`, `ρ_Higgs_SSB`,
    `ρ_EW`) **NIE są additive** do `ρ_vac_TGP` — strukturalna konsekwencja
    single-Φ axiom + substrate-vacuum identification (Φ_eq=H₀, γ=M_Pl²·g̃)
  - Mechanism: matter sector vacua są *transient sources* during phase transitions
    (T~Λ_QCD, Λ_EW), NIE contribute do *today's* Λ
  - Empirical test: T-Λ 7/7 PASS + ratio 1.020 = direct evidence dla argumentu
  - Cosmological constant problem rozszerzony structural resolution z T-Λ
    (ŝ-quanta) na ALL SM matter sectors (Q2 closure)

  Patrz [[../op-Q2-vacuum-budget-2026-05-10/Phase_FINAL_close.md]].

- **Q3**: Dla magnetarów (B ~ 10¹¹ T), oczekiwany observable Φ-shift
  od `ρ_EM_quantum` jest jaki rząd? — to było pierwotnie framowane jako
  *concrete falsifiability test* dla magnetar polar shift TT10 z τ.3 cyclu.

  **Status post-2026-05-10:** **CLOSED** z native-first re-analysis.
  - Native estimate (L01 ADDENDUM §3.2): `ρ_EM_quantum/ρ_NS ~ 10⁻¹²` w typowym magnetar
  - Mechanism decoupling (τ.3 ADDENDUM §2): Phase3.TT10 testuje **L4 gradient-coupled
    mass mechanism**, NIE `ρ_EM_quantum`. Dwa mechanizmy decoupled przez **10 OOM**.
  - **Wniosek:** TT10 NIE jest falsifikatorem dla L01 N1; dedicated test dla
    `ρ_EM_quantum` wymaga lab Schwinger-class field w macroscopic volume (beyond 2030+).

  Patrz [[../op-tau3-substrate-clock-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §2.

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | N1 (EM trace anomaly) | ~~dedicated 1-loop QFT cycle~~ **CLOSED 2026-05-11** | [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] |
| B2 | N2 (QCD trace anomaly) | ~~dedicated QCD + cosmology cycle~~ **CLOSED 2026-05-11** | [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/]] (STRUCTURAL_DERIVED, 24/24 sympy PASS, 6/6 P-requirements) |
| B3 | N3 (SPARC consistency) | ~~krótki re-derivation w nbody/~~ **CLOSED 2026-05-11** | [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] (STRUCTURAL_DERIVED, 8/8 sympy PASS) |
| B4 | N4 (Higgs sector) | ~~dedicated 1-loop Higgs + EW cosmology cycle~~ **CLOSED 2026-05-11** | [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/]] (STRUCTURAL_DERIVED, 24/24 sympy PASS, 6/6 P-requirements, R5 LOCK + Q2 F1 Higgs konstruktywnie verified) |
| B5 | N5 (EW anomaly) | ~~część N2 cycle / dedicated EW gauge cycle~~ **CLOSED 2026-05-11** | [[../op-L01-N5-EW-gauge-anomaly-2026-05-11/]] (STRUCTURAL_DERIVED, 8/8 sympy PASS, 6/6 P-requirements; LAST L01 N-need closed; full SM matter+gauge coverage) |

## Closed needs (po L01 cyklu 2026-05-04)

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| Formal kowariantna definicja ρ | ρ ≡ -T^μ_μ/c_0² (formal_definition.md) | L01 cykl | 2026-05-04 |
| Mapping SM sektorów | 5 sektorów (Dirac, scalar, EM, YM, fluid) | SM_sector_mapping.md | 2026-05-04 |
| Photon treatment (ρ_EM=0) | analiza 3 implikacji + 1 open (quantum) | photon_treatment.md | 2026-05-04 |
| Spójność z ax:metric-coupling | Option-2 audytu *DERIVED* nie *postulowane* | formal_definition.md §5 | 2026-05-04 |
| S04 N1 (formal kowariantna derywacja L_mat) | wyprowadzenie z ax:metric-coupling + perturbation | formal_definition.md §4 | 2026-05-04 |
| Q3 (magnetar polar shift status) | native re-analysis: `ρ_EM_quantum/ρ_NS ~ 10⁻¹²`; TT10 decoupled od L01 N1 (testuje L4, nie trace anomaly) | ADDENDUM_2026-05-10 §3.2 + τ.3 ADDENDUM §2 | 2026-05-10 |
| Q1 (quantum trace anomaly EM coverage przez ψ.1.v2 dim-6 EFT) | operator class disjoint: ψ.1.v3 basis ma ≥1 ∂lnX leg, quantum trace anomaly to pure-photon dim-4; N1 dedicated cycle wymagany strukturalnie | ψ.1 ADDENDUM §3 (Phase 7 enumeration) | 2026-05-10 |
| Q1 (post-N1) **konstruktywna verification** Theorem 2.1 (Disjointness) | dedicated 1-loop QED derivation w g_eff[{Φ_i}] potwierdza disjointness explicit; trace anomaly TGP-reduced operator classes ⊥ ψ.1.v3 basis B={L₅'_a, L₅'_b} | [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase2_results.md]] §2 (sympy T4 LOCK) | 2026-05-11 |
| **N1 (Quantum trace anomaly EM dedicated cycle)** | STRUCTURAL_DERIVED: ρ_EM_quantum[{Φ_i}] explicit z prefactor α/(3π)≈7.74·10⁻⁴; Theorem 2.1 disjointness; GW170817 Δc/c~10⁻⁸⁰; MICROSCOPE η_TGP_EM_quantum=0 strukturalnie z S05 universal coupling; R5 documented (B≪B_QED); R6 (Asorey-2015 QEP) closed strukturalnie. **L01 ADDENDUM §3.2 typo correction:** `α²/(3π)≈7.7·10⁻⁷` → `α/(3π)≈7.74·10⁻⁴` (factor ~1000 OOM correction) | [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] (16/16 sympy PASS, 6/6 P-requirements, STRUCTURAL_DERIVED) | 2026-05-11 |
| Q2 (ρ_QCD vs T-Λ closure) | OSOBNE: SM matter sector vacua decoupled od bare Λ przez single-Φ axiom + substrate-vacuum identification; cosmological constant problem rozszerzony structural resolution | [[../op-Q2-vacuum-budget-2026-05-10]] cycle (STRUCTURAL DERIVED) | 2026-05-10 |

<a id="T"></a>

## §T — Three-layer specification per NEEDS (post-2026-05-10)

Aplikacja [[../../meta/PPN_AS_PROJECTION.md]] §3.1 mandatory three-layer presentation
do otwartych potrzeb L01. Pełne uzasadnienie i kontekst w
[[ADDENDUM_2026-05-10_native_observables_first.md]].

Każdy N_i jest *deferred ρ-correction*. Tabele poniżej specyfikują dla każdego N_i:
- **L1 (native predictions):** *jaką obserwowalną* korekcja ρ_X modyfikuje, liczoną
  bezpośrednio z `g_eff[Φ]` + Φ-EOM po podstawieniu poprawionego źródła
- **L2 (projection chart):** *na który* PPN/ppE/cosmology parametr to się projektuje
- **L3 (falsification map):** *który* observational bound limituje native ρ_X-coef,
  i w jakim window

### §T.1 — N1 (Quantum trace anomaly EM, 1-loop QED na curved background) — **CLOSED 2026-05-11**

**Native ρ_X correction (z dedicated cycle [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] Phase 1+2):**
```
ρ_EM_quantum[{Φ_i}](x) = -[
   (α/(3π)) · F²(x)
   + γ_1 · (∇²ψ) · F²
   + γ_2 · (∂_μ ∂_ν ψ) · F^{μρ} F^ν_ρ
   + γ_3 · σ_ab · F²
   + γ_4 · □F²
   + Riegert local with σ_eff = function(ψ)
] / c_0²
```

z **β(α)/(2α) = α/(3π) ≈ 7.74·10⁻⁴** (sympy LOCK; **NIE** 7.7·10⁻⁷ jak początkowo cytowane —
typo correction propagated 2026-05-11).

| Layer | Specification (post-N1 closure) |
|---|---|
| **L1 native** | (a) Lab regime (B=1 T): \|ρ_EM_quantum\| ~ 7·10⁻¹⁵ kg/m³ inside electromagnet; niewykrywalne acceleration. (b) Magnetar typical (B=10¹⁰ T): ρ_EM_quantum/ρ_NS ~ 1.7·10⁻¹². (c) Magnetar extreme (B=10¹¹ T, B/B_QED~23): ratio ~1.7·10⁻¹⁰ z corrected α/(3π) prefactor (NIE 10⁻¹²). (d) Solar system: deflekcja od ρ_EM_quantum suppressed by F²/(m_p²c²)·(α/(3π)) ~ 10⁻⁴⁰; insignificant. (e) Lab Schwinger Δclock (E~10¹⁵ V/m + B~10⁹ T macroscopic): future test, zasięg 2030+. |
| **L2 projection** | (a) **NIE wprowadza** nowych ppE/PPN free parameters; modyfikuje *precision* native predictions (deferred precision Wilson coefs γ_i, NIE new free coefs). (b) Theorem 2.1 (Disjointness, dedicated cycle Phase 2 §2): trace anomaly TGP-reduced operator classes ⊥ ψ.1.v3 basis B={L₅'_a, L₅'_b}. (c) Cosmology: w_eff(z) chart unchanged (radiation era T^μ_μ_radiation=0 automatic). |
| **L3 falsifikator** | (a) **MICROSCOPE η ≤ 1.1·10⁻¹⁵** (Pt vs Ti): η_TGP_EM_quantum = **0 strukturalnie** z S05 universal coupling; PASS automatic z B9 baseline 1.32·10⁻²⁶ (~11 OOM margin). MICROSCOPE-2 2027+ projection η ≤ 10⁻¹⁷: PASS automatic. (b) **GW170817 \|c_GW/c_EM-1\| ≤ 9·10⁻²²**: Δc/c ~ 10⁻⁸⁰ (~58 OOM margin). (c) **Eöt-Wash + LLR Nordtvedt**: PASS automatic z universal coupling. (d) **τ.3 TT10 magnetar X-ray timing** (XMM-Newton/Chandra): trace anomaly **NIE jest falsifierem** dla τ.3 (mechanism decoupled — testuje L4 gradient-coupled mass mechanism czysto; ratio 10⁻¹² typical, 10⁻¹⁰ extreme ≪ 1). (e) **R6 (Asorey-2015 QEP violations)**: closed strukturalnie z TGP universal coupling immunity. |

**Closure:** dedicated cycle `op-L01-N1-EM-trace-anomaly-TGP-2026-05-11`
([[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase_FINAL_close.md]]) — STRUCTURAL_DERIVED,
16/16 sympy PASS, 6/6 P-requirements RESOLVED.

**Open (deferred precision, NOT blockers):**
- Wilson coefs γ_1,γ_2,γ_3,γ_4 numerical pinning (multi-loop QED extension cycle)
- B ≳ B_QED magnetar regime non-perturbative analysis (Schwinger pair production effects)
- Schwinger-class lab predictions (2030+ experimental capability)

### §T.2 — N2 (QCD trace anomaly, gluon condensate) — **CLOSED 2026-05-11**

**Native ρ_X correction (z dedicated cycle [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/]] Phase 1+2):**
```
ρ_QCD(T)[{Φ_i}] = ρ_QCD_vacuum + ρ_QCD_thermal(T)
   ρ_QCD_vacuum  = -(b_0/8) · ⟨α_s G²/π⟩ / c_0²
                 ≈ 2.8·10¹⁸ kg/m³  (constant SVZ; substrate-decoupled per Q2 F1)
   ρ_QCD_thermal(T) = -Δ(T) · T⁴ / c_0²
                   Δ_max ≈ 4 near T_c = 156 ± 9 MeV (HotQCD lattice crossover)
                   Δ(T → 0) → 0 strukturalnie
                   Δ(T → ∞) → 0 (Stefan-Boltzmann conformal)
```

z **β_QCD(g) = -(b_0/(16π²))·g³**, b_0 = (11/3)N_c - (2/3)N_f = 7 (N_c=3, N_f=6 high-T).
**Λ_QCD = 217 ± 8 MeV** (PDG 2024 MS-bar N_f=5).

| Layer | Specification (post-N2 closure) |
|---|---|
| **L1 native** | (a) ρ_QCD_vacuum equivalent ~2.8·10¹⁸ kg/m³ (constant, substrate-decoupled). (b) Δ_max≈4 peak near T_c=156 MeV (lattice 2+1 flavor crossover). (c) ρ_QCD_anomaly/ε_total ≈ 26-35% transient at T_c (FWHM ~50 MeV). (d) Hubble t_H(T_c) ~ 10⁻⁵ s consistent z standard cosmology QCD epoch. (e) ρ_QCD_thermal(T<<Λ_QCD) → 0 strukturalnie (Q2 F1 verification). |
| **L2 projection** | (a) **NIE wprowadza** nowych ppE/PPN/cosmology free parametrów (analog do N1); modyfikuje *precision* w QCD epoce (deferred precision lattice inputs, NIE new free coefs). (b) w_eff(z) unchanged: matter-decoupling preserved. (c) Q2 F1 mechanism KONSTRUKTYWNIE verified dla QCD sektora (analog do N1 Theorem 2.1 dla EM). (d) Σ(z,k), μ(z,k) modified gravity charts unchanged. |
| **L3 falsifikator** | (a) **BBN ⁴He Y_p ≤ ±5σ** (PDG 2024 0.245±0.003): TGP=ΛCDM=0.247±0.001 within **0.55σ** — PASS. (b) **D/H** (PDG 2024 2.527·10⁻⁵): TGP=2.5·10⁻⁵ within **0.26σ** — PASS. (c) **Planck 2018 ω_b, ω_m, Ω_Λ**: preserved automatic — PASS exactly. (d) **T-Λ ratio 1.020 ± 0.02** (closure_2026-04-26): preserved — direct evidence dla Q2 F1. (e) **PTA NANOGrav 15-yr (2023)**: lattice 2+1 flavor crossover (NOT first-order) → no dominant QCD signal → consistent z SMBHB consensus origin — PASS. (f) **Hypothetical naive-additive scenario** daje ratio ~10⁷⁷ (catastrophe) → empirical 1.020 *direct evidence* dla Q2 F1 mechanism. |

**Closure:** dedicated cycle `op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11`
([[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase_FINAL_close.md]]) —
STRUCTURAL_DERIVED, **24/24 sympy PASS** (Phase 1+2+3), **6/6 P-requirements RESOLVED**.

**Open (deferred precision, NOT blockers):**
- Multi-loop β_QCD precision (b_1, b_2, b_3 coefs) + sub-percent lattice EoS
- Higher-order Wilson coefficients dla curvature × G² terms (analogous N1 γ_i)
- N4 Higgs sector + N5 EW gauge anomaly extensions (separate cycles)

**Deferred to:** ~~dedicated cycle `op-QCD-trace-anomaly-cosmology`~~ **CLOSED**
[[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/]]; precision extension cycles deferred (4-6 tyg., wymaga
lattice QCD inputs).

### §T.3 — N3 (SPARC fits consistency check) — **CLOSED 2026-05-11**

**Native ρ_X correction (z dedicated cycle [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] Phase 1):**
```
ρ_SPARC(r) ≡ ρ_baryon(r) = ρ_HI(r) + ρ_stars(r) + ρ_bulge(r)
            = -T^μ_μ_dust(r) / c_0²    [non-relativistic dust limit, verified]

Verified precision:
  Galactic stars (v ~ 200 km/s): ρ_TGP/ρ_rest deviation ~ 2.2·10⁻⁷ (= 2.2·10⁻⁵ %)
  HI gas thermal (v ~ 1 km/s):  deviation ~ 5.6·10⁻¹² (utterly negligible)
  Target precision: 1%
  Achieved precision: ~10⁻⁶ (6 OOM below target)
```

| Layer | Specification (post-N3 closure) |
|---|---|
| **L1 native** | (a) Rotation curve `v(r)` liczona bezpośrednio z `g_eff[{Φ_i}]` (multi-source emergent gravity) z source `ρ_baryon` (HI + stars + bulge). (b) **Dust-limit verification:** for galactic regime v²/c² ~ 10⁻⁷, ρ_TGP = ρ_rest exact do 10⁻⁶ precision (sympy Phase 1 T1-T4). (c) NO separate ρ_DM matter component (R1 closed). |
| **L2 projection** | (a) Galactic-scale modified gravity charts: effective surface density Σ_eff, extra acceleration at low-acceleration limit (MOND-like phenomenology). (b) **Modyfikacja v(r) flat curves z g_eff[Φ̄] mechanism (gravitational), NIE z ρ_DM (matter); S05 enforced.** (c) BTFR slope ≈ 4 + RAR Lelli+2017 match to ~15% (galaxy_scaling cycles gs10-gs61). |
| **L3 falsifikator** | (a) **SPARC residuals**: 175 galaxy rotation curve fits, chi²_red competitive z MOND simple (galaxy_scaling CLOSURE_2026-04-19). (b) **Solar System Cassini γ**: konsystencja `ρ_baryon` z jednostkową coupling Newton-limit on solar scale (preserved, separate cycle). (c) **Strong/weak lensing convergence**: cross-check że `ρ_SPARC ≡ ρ_baryon` tłumaczy soczewkowanie bez nadwyżkowych komponentów. (d) **Cluster-scale ~35% mass deficit** — separate issue (gs13-gs55, requires ~2 eV sterile ν), NIE part of N3 scope. |

**Closure:** dedicated compact cycle `op-L01-N3-SPARC-rho-consistency-2026-05-11`
([[../op-L01-N3-SPARC-rho-consistency-2026-05-11/Phase_FINAL_close.md]]) —
STRUCTURAL_DERIVED, **8/8 sympy PASS**, **6/6 P-requirements RESOLVED**.

**R1 (double-counting) closed strukturalnie:** TGP-emergent DM (z `g_eff[Φ̄]`
background) **NIE może** być dodawane jako oddzielny ρ_DM source — to byłoby
naruszenie S05 (single fundamental field). Galaxy_scaling cycles (gs10-gs61)
strukturalnie unikają tego — używają tylko ρ_baryon jako matter source.

**R2 honestly documented:** SPARC scope = galactic-disk regime (R~1-50 kpc, NR
valid); near-SMBH ISCO regime wymaga full GR / TGP-PPN treatment (poza zasięgiem
N3 cycle).

**Open (cosmetic, deferred):**
- galaxy_scaling/README + nbody/README small documentation note z explicit T^μ_μ → ρ mapping

### §T.4 — N4 (Higgs sector)

**Native ρ_X correction:**
```
ρ_Higgs_quantum ≈ ⟨(∂h)² - m_H²·h²⟩_1-loop / c_0²   [quantum fluctuations around vacuum]
ρ_Higgs|_classical_vacuum = 0  (SSB cancellation, F10)
```

| Layer | Specification |
|---|---|
| **L1 native** | (a) Vacuum energy contribution z 1-loop Higgs fluctuations: `δρ_vac_Higgs ~ m_H⁴/(64π²)` (bare estimate, requires renormalization). (b) Modyfikacja Φ-EOM w EW phase transition era (T~100 GeV). |
| **L2 projection** | (a) Contribution do T-Λ closure `ρ_vac_TGP`: czy renormalized `δρ_vac_Higgs` jest częścią observed Λ, czy structurally cut przez vacuum subtraction (`⟨T^μ_μ⟩_0 = 0`)? (b) Potencjalne `δw(z)` w late-time cosmology. |
| **L3 falsifikator** | (a) **SN+BAO Λ measurement** (1% precision na Ω_Λ): konstrainuje ρ_vac_TGP total. (b) **Planck 2018 ω_Λ**: ω_Λ = 0.6889 ± 0.0056. (c) **CMB-S4 future projections**: σ(w) ~ 10⁻³ — może detekować quantum Higgs contribution. |

**CLOSURE 2026-05-11:** **CLOSED konstruktywnie** w dedicated cycle
[[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/]] (STRUCTURAL_DERIVED, 24/24
sympy PASS, 6/6 P-requirements RESOLVED). Higgs sector trace anomaly +
EW phase transition cosmology + Q2 F1 konstruktywna verification + R5 LOCK
(crossover) + LISA Ω_GW^EW=0 falsifiable prediction.

### §T.5 — N5 (EW gauge anomaly, SU(2)×U(1))

**Native ρ_X correction:**
```
ρ_EW(T) ≈ Λ_EW⁴(T) / c_0²   [analogous do QCD, ale dla SU(2)×U(1)]
β-functions: SU(2) b₀ = -19/12, U(1) b₀ = +41/12
```

| Layer | Specification |
|---|---|
| **L1 native** | (a) Modyfikacja `H(z)` w EW epoce (T~100 GeV, z~10¹⁵). (b) Stochastic GW background z first-order EW phase transition (jeśli applicable w SM extension). |
| **L2 projection** | (a) `w_eff(z~10¹⁵)` cosmology chart. (b) Contribution do gravitational wave background spectrum. |
| **L3 falsifikator** | (a) **BBN ⁴He, D/H**: niewrażliwe bezpośrednio (T_EW >> T_BBN). (b) **LISA stochastic GW background**: detekcja first-order phase transition. (c) **PTA stochastic GW**: niewrażliwe na EW (frequency mismatch). |

**Deferred to:** część cyklu N2 (`op-QCD-trace-anomaly-cosmology` rozszerzony na EW
sector + cosmological phase transitions).

### §T.6 — Native parameter audit (post-N1-N5 closure)

Per `meta/PPN_AS_PROJECTION.md §3.3`:

```
Independent native source structures fixed by L01 (current closure):
  ρ ≡ -T^μ_μ/c_0²  (1 derived definition)
  5 SM sector mappings: Dirac, scalar, EM_classical, YM_classical, fluid
  q (coupling const) → G_0 mapping
  φ-prefactor in L_mat (derived from √(-g_eff))

Pending additions from N1-N5 closure (NO new free params, only renormalized refinements):
  N1 → ρ_EM_quantum coefficient (β(α)/(2α))   [renormalization-fixed, not free] ✓ CLOSED 2026-05-11
  N2 → ρ_QCD(T) profile (lattice input)        [lattice-fixed, not free]        ✓ CLOSED 2026-05-11
  N3 → ρ_SPARC ≡ ρ_baryon verification         [no new param, consistency check] ✓ CLOSED 2026-05-11
  N4 → ρ_Higgs_quantum vacuum subtraction      [renormalization scheme dependence] ✓ CLOSED 2026-05-11
  N5 → ρ_EW(T) profile (β-function input)      [perturbatively fixed, not free]  ✓ CLOSED 2026-05-11

Forced from substrate symmetry (unchanged by N1-N5):
  GW170817 c_GW = c_EM ≡ exact
  WEP universality (B9 verified)
  Radiation era non-source for Φ
  ζ_i ≡ 0 (energy/momentum conservation)
  α_i ≡ 0 (preferred frame, from Lorentz-invariance Φ-substrate)

Conclusion: L01 closure of N1-N5 increases PRECISION of native predictions
            in extreme regimes; does NOT add new free parameters.
```

## Cross-references

- [[README.md]] — cykl L01 indeks
- [[FINDINGS.md]] — eksportowalne wyniki
- [[ADDENDUM_2026-05-10_native_observables_first.md]] — native-first methodology overlay (binding 2026-05-10+)
- [[../../meta/PPN_AS_PROJECTION.md]] — parent methodology doc
- [[../../audyt/L01_rho_operational/]]
- [[../../audyt/S04_metric_coupling_axiom/]] (powiązany — N1 formal derivation)
- [[../op-psi1-substrate-light-acceleration]] (ψ.1.v2 dim-6 EFT)
- [[../op-tau3-substrate-clock-acceleration]] (τ.3 magnetar polar)
- [[../op-omega1-substrate-em-coupling]] (ω.1 EM coupling)
- [[../op-emergent-metric-from-interaction-2026-05-09]] (g_eff[Φ] funkcjonał, L01 jest L1 input)
- [[../closure_2026-04-26/Lambda_from_Phi0]] (T-Λ vacuum energy, Q2 territory)
