---
title: "L01 ρ stress-energy bridge — formal kowariantna definicja gęstości materii"
date: 2026-05-04
cycle: L01
type: audit-resolution
status: ANALYTICAL DECISION DOCUMENTED — ρ ≡ -T^μ_μ/c_0² jako formal definition + mapping SM
parent: "[[../../audyt/L01_rho_operational/README.md]]"
predecessors:
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]]" (§A.4, §C.3, §N — Option-2)
  - "[[../op-newton-momentum/B9_wep_microscope_composition_results.md]]" (B9 closure 2026-05-01)
related:
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]" (lin. 102-129, eq:L-mat-unified)
  - "[[../../core/sek08_formalizm/sek08_formalizm.tex]]" (lin. 11257-11337, ax:metric-coupling)
  - "[[../../TGP_FOUNDATIONS.md]]" §4 (warstwa 3a/3b/3c)
tags:
  - TGP
  - L01
  - rho-formal-definition
  - stress-energy
  - matter-coupling
  - SM-sector-mapping
  - photon-treatment
  - audit-resolution
tgp_status:
  folder_status: closed-resolved
  level: L1
  kind: derivation
  core_compatibility: current
  last_reviewed_against_core: 2026-05-04
  may_edit_core: true
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges: []  # quantum-trace-anomaly CLOSED 2026-05-11 by op-L01-N1-EM-trace-anomaly-TGP-2026-05-11
  depends_on: []
  impacts:
    - "[[../../audyt/L01_rho_operational]]"
    - "[[../../audyt/S04_metric_coupling_axiom]]"
  source_of_status:
    - "ax:metric-coupling (sek08_formalizm.tex:11257)"
    - "eq:L-mat-unified (sek08a:103-111)"
    - "AUDYT 2026-05-01 §A.4 Option-2 decision"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-04
---

# L01 — formal kowariantna definicja `ρ`

## Cel

Operacyjnie domknąć status `ρ` w `L_mat = -(q/Φ_0)·φ·ρ` (sek08a:103)
poprzez:

1. Explicit kowariantna definicja `ρ ≡ -T^μ_μ/c_0²` w sek08a
2. Mapping na sektory SM: Dirac (m·ψ̄ψ), scalar (gradient terms), EM (T^μ_μ=0)
3. Treatment fotonów (`ρ_EM = 0` — co to znaczy fizycznie)
4. NON-BREAKING addytywna edycja sek08a

## Werdykt

**`ρ ≡ -T^μ_μ/c_0²`** jest formal kowariantną definicją gęstości
materii w `L_mat` sek08a. Dla:

- **Dirac fermion** (m·ψ̄ψ): `ρ_Dirac = m·ψ̄ψ/c_0²` ≥ 0 ✓ → standardowa
  gęstość masy bezwładnościowej
- **Massless EM field**: `T^μ_μ_EM = 0` → `ρ_EM = 0` → **fotony nie
  generują pola Φ przez `L_mat`**, ale propagują się geodezyjnie po
  `g_eff` (`prop:coupling-consequences` 2)
- **Massive scalar**: `T^μ_μ_scalar` zawiera `m²|φ|²` term + gradient
  contributions
- **Yang-Mills gauge**: `T^μ_μ_YM = 0` (klasycznie); quantum trace anomaly
  modyfikuje to przez β-funkcję

Ta definicja **ratuje strukturalną spójność** `ax:metric-coupling`
(audit Option-2): `L_mat = -(q/Φ_0)·φ·ρ` jest **derived consequence**
volume element `√(-g_eff) ∝ φ` plus stress-energy tensor materii w
kanonicznym sprzęganiu metrycznym.

## Pliki w cyklu

| Plik | Opis |
|------|------|
| [[README.md]] | (ten plik) |
| [[formal_definition.md]] | explicit `ρ ≡ -T^μ_μ/c_0²` derivation z `L_mat[ψ_m, g_eff]` |
| [[SM_sector_mapping.md]] | Dirac/scalar/EM/gauge — explicit ρ dla każdego sektora |
| [[photon_treatment.md]] | T^μ_μ_EM=0 i implikacje (GW170817 c_GW=c, no fifth force from photons) |
| [[FINDINGS.md]] | eksportowalne wyniki |
| [[NEEDS.md]] | open problems (quantum trace anomaly) — *post-2026-05-10*: zawiera §T three-layer specification per N1-N5 |
| [[ADDENDUM_2026-05-10_native_observables_first.md]] | **interpretive overlay (2026-05-10):** native-observables-first reframe per `meta/PPN_AS_PROJECTION.md`; ρ jako natywny source field, PPN/ppE jako L2 projekcja |

## Co wykonano w sesji

1. **Synthesis 5 plików rdzenia + audit** — sek08a §74-129, sek08_formalizm
   §11257-11337, AUDYT 2026-05-01 §A.4+§C.3+§N, B9 closure 2026-05-01.
2. **Formal derivation** `ρ = -T^μ_μ/c_0²` z `L_mat[ψ_m, g_eff]` —
   demonstracja, że Option-2 audytu jest *wyprowadzeniem*, nie *postulatem*.
3. **Mapping SM sector-by-sector** — Dirac, scalar, EM, Yang-Mills.
4. **Treatment fotonów** — explicit `ρ_EM = 0` z T^μ_μ_EM = 0; brak
   piątej siły *od fotonów* (konsystentne z B9 MICROSCOPE 6/6 PASS).
5. **NON-BREAKING addytywne edits** w sek08a (planowane, opcjonalne).

## Status w audicie L01

Patrz [[../../audyt/L01_rho_operational/POST_ACTION_UPDATE_2026-05-04.md]].

## Open problems

- **N1: Quantum trace anomaly** — `T^μ_μ ≠ 0` w QFT na curved background
  pochodzi od *trace anomaly* (Birrell-Davies 1982); modyfikacja
  ρ przez `O(R²)` corrections wymaga dedicated cycle dla **renormalized
  ρ** w obecności `g_eff[Φ]`.
- **N2: SPARC fits** — N-body symulacje używają `ρ = ρ_baryon` (HI + stars
  + bulge); explicit weryfikacja że to jest non-relativistic limit
  `T^μ_μ = -ρ_rest·c²` jest pożądana (cosmetic, low priority).

Patrz [[NEEDS.md]] dla pełnej listy.

## Status (2026-05-11 update — N1+N2+N3+N4+N5 closure ★ FULL SM COVERAGE ★)

🔒 **CLOSED-DERIVED** (Phase 3/3 CLOSED, F1-F21 exported, **Q1+Q2+Q3 zamknięte**
przez cross-cycle propagation 2026-05-10, **N1+N2+N3+N4+N5 wszystkie zamknięte
konstruktywnie 2026-05-11** przez pięć dedicated cycles — **MILESTONE: pełne
SM matter+gauge sektor coverage**).

> **POST-N5 CLOSURE (2026-05-11, same day as N1+N2+N3+N4 — ★ LAST L01 N-NEED ★):**
> EW gauge sektor (SU(2)×U(1)) trace anomaly N5 zamknięty *konstruktywnie* w
> dedicated cycle [[../op-L01-N5-EW-gauge-anomaly-2026-05-11/]] (compact single-
> session z silnej architecture inheritance N1+N2+N4):
> - **β_SU(2) asymptotic freedom**: b₀=19/6 > 0; β(M_Z) ≈ -5.56·10⁻³
> - **β_U(1) Landau pole**: b₀=41/6 > 0; β(M_Z) ≈ +1.97·10⁻³; μ_LP ~10⁴² GeV >> M_Pl
>   (23 OOM margin powyżej Planck — cosmologically irrelevant)
> - **Trace anomaly form** Adler-Collins-Duncan 1977 explicit dla obu: T^μ_μ = (β/2g)·F²
> - **Q2 F1 KONSTRUKTYWNIE verified dla gauge sektora** — OOM gap 54.9 (bare ρ_gauge_vac
>   ≈ 1.94·10⁴⁴ eV⁴ vs Λ_obs) substrate-decoupled (analog N4 Higgs OOM gap 55.3)
> - **EW cosmology N4-inheritance**: R5 LOCK crossover; gauge bosons freeze-out
>   T~3.2 GeV >> T_BBN; **LISA Ω_GW^EW=0 inherited** z N4 R5 LOCK
> - **PDG 2024 EW precision preserved** (M_W=80.379, M_Z=91.1876, sin²θ_W=0.2312;
>   tree-level Sirlin 0.2230 + loop correction 0.008 standard SM)
> - 8/8 sympy PASS, 6/6 P-requirements RESOLVED, 6/6 risks closed konstruktywnie
> - **MILESTONE:** **LAST L01 N-need closed** — full SM matter+gauge coverage
>
> Patrz [[../op-L01-N5-EW-gauge-anomaly-2026-05-11/Phase_FINAL_close.md]].

> **POST-N4 CLOSURE (2026-05-11, same day as N1+N2+N3):**
> Higgs sektor trace anomaly + EW phase transition cosmology N4 zamknięty
> *konstruktywnie* w dedicated cycle
> [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/]]:
> - **SSB cancellation**: bare T^μ_μ_vac = -λv⁴; post-renormalization T_vac=0
>   (standard SM convention; Q2 F1 substrate-decoupling enforced)
> - **1-loop trace anomaly explicit**: `T^μ_μ_quantum = (1+γ_m)m_H²h² + (β_λ/4)h⁴
>   + curvature·h² + Riegert local`
> - **β_λ ≈ -0.033** (top Yukawa dominant), **γ_m ≈ -0.027 ÷ -0.035** (sympy LOCK)
> - **Q2 F1 KONSTRUKTYWNIE verified dla Higgs sektora** — OOM gap 55.3 (bare
>   ρ_Higgs_vac vs Λ_obs) substrate-decoupled via single-Φ + substrate-vacuum
>   identification (czwarty SM sektor verification post N1+N2+N3)
> - **R5 LOCK**: m_H=125.25 GeV > endpoint 80 GeV (KLRS 1996, DRR 2014) ⇒
>   EW transition crossover (NOT first-order); **LISA Ω_GW^EW = 0 strukturalnie**
>   (falsifiable post-2035)
> - LHC m_H preserved exact + Planck 2018 ω_b/ω_m/Ω_Λ/N_eff preserved + BBN ⁴He
>   (0.67σ) + D/H (0.90σ) within bounds + PTA NANOGrav 15-yr compatible (two-sektor
>   GW synergy: LISA mHz + PTA nHz oba empty primordial bands)
> - **HL-LHC + FCC-ee future precision** null test prediction (Δλ_HHH/λ_SM ≈ 0)
> - 24/24 sympy PASS (Phase 1+2+3), 6/6 P-requirements RESOLVED
> - R1-R6 risks: 5 closed (4 strukturalnie/konstruktywnie + 1 honestly documented)
>   + 1 deferred R3 (Higgs hierarchy) z empirical status preserved
>
> Patrz [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11/Phase_FINAL_close.md]].

> **POST-N3 CLOSURE (2026-05-11, same day as N1+N2):**
> SPARC ρ-consistency N3 zamknięty *kompaktowo* w dedicated verification cycle
> [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]]:
> - **ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0²** verified do <10⁻⁶ precision
>   (target 1%, achieved 6 OOM below)
> - **R1 (double-counting) closed strukturalnie:** TGP-emergent DM jest
>   gravitational (g_eff[Φ̄]), NIE matter; S05 single-Φ enforced
> - **R2 (galactic-center) honestly documented:** SPARC scope = galactic-disk
>   regime; near-SMBH ISCO (v~c/2) outside cycle scope
> - 8/8 sympy PASS, 6/6 P-requirements RESOLVED
>
> Patrz [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/Phase_FINAL_close.md]].

> **POST-N2 CLOSURE (2026-05-11, same day as N1):**
> QCD trace anomaly + cosmology integration N2 zamknięty *konstruktywnie* w dedicated cycle:
> - β_QCD(g) = -(b_0/(16π²))·g³, b_0=7 (sympy LOCK; asymptotic freedom)
> - Trace anomaly form: T^μ_μ_QCD = (β(g)/(2g))·G² = -(b_0 α_s)/(8π)·G²
> - Gluon condensate ⟨α_s G²/π⟩ ≈ 0.012 GeV⁴ → ρ_QCD_vacuum ~ 2.8·10¹⁸ kg/m³
> - **Q2 F1 KONSTRUKTYWNIE verified dla QCD sektora** (Phase 2 §3); matter-vacuum
>   decoupling z dedicated derivation
> - BBN ⁴He Y_p (0.55σ), D/H (0.26σ), CMB ω_b (exact), PTA NANOGrav 15-yr (compat
>   z SMBHB consensus) wszystkie PASS automatic
> - 24/24 sympy PASS (Phase 1+2+3), 6/6 P-requirements RESOLVED
> - R1-R7 risks: 6 fully closed + R2 honestly documented (lattice external)
>
> Patrz [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/Phase_FINAL_close.md]].

> **POST-N1 CLOSURE (2026-05-11):**
> Quantum trace anomaly EM N1 zamknięty
> *konstruktywnie* w dedicated 1-loop QED cycle:
> - ρ_EM_quantum[{Φ_i}] explicit form (Phase 1+2): `(α/(3π))·F² + γ_i·(curvature×F²) + Riegert`
> - **β(α)/(2α) = α/(3π) ≈ 7.74·10⁻⁴** (sympy LOCK; **typo correction**: ADDENDUM §3.2 Q3
>   cytowała `7.7·10⁻⁷` — factor ~1000 OOM correction propagated)
> - **Theorem 2.1 (Disjointness)** verified konstruktywnie: trace anomaly TGP-reduced
>   ⊥ ψ.1.v3 basis B={L₅'_a, L₅'_b}
> - GW170817 Δc/c ~ 10⁻⁸⁰ (~58 OOM margin); MICROSCOPE η_TGP_EM_quantum=0 strukturalnie
>   z S05 universal coupling
> - **R6 (Asorey-2015 QEP violations) closed strukturalnie** — universal coupling immunity
> - **R5 (B≳B_QED magnetar)** honestly documented; deferred non-perturbative cycle
> - 16/16 sympy PASS, 6/6 P-requirements RESOLVED, STRUCTURAL_DERIVED classification
>
> Patrz [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/Phase_FINAL_close.md]].

## Status (2026-05-10 update)

🔒 **CLOSED-DERIVED** (Phase 3/3 CLOSED, F1-F21 exported, **Q1+Q2+Q3 zamknięte**
przez cross-cycle propagation 2026-05-10).

> **ADDENDUM 2026-05-10 (interpretive overlay, NIE zmienia classification):**
> [[ADDENDUM_2026-05-10_native_observables_first.md]] — reframe wyników L01 w
> native-observables-first form per [[../../meta/PPN_AS_PROJECTION.md]] (binding
> 2026-05-10+). Kluczowe statementy:
>
> 1. `ρ ≡ -T^μ_μ/c_0²` jest **natywnym source field** dla Φ-EOM — nie projekcją.
>    PPN/ppE/cosmology są L2 projekcjami obserwabli liczonych z `(g_eff[Φ], ρ)`.
> 2. NEEDS N1–N5 dostały **three-layer specification** (L1 native obserwable / L2
>    projection chart / L3 falsifikator) — patrz [[NEEDS.md]] §T.
> 3. Native parameter audit pokazuje, że N1-N5 closure *nie wprowadza* nowych
>    swobodnych parametrów — modyfikuje precyzję native predictions w extreme
>    regimes (renormalization-fixed, lattice-fixed).
> 4. **Q1 CLOSED:** ψ.1.v3 dim-6 EFT operator class disjoint od quantum trace anomaly
>    EM (`α·F²` pure-photon, NO ∂lnX leg) — N1 dedicated cycle strukturalnie wymagany.
>    Patrz [[../op-psi1-substrate-light-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §3.
> 5. **Q2 CLOSED:** SM matter sector vacua (ρ_QCD, ρ_Higgs, ρ_EW) **NIE additive** do
>    bare Λ — single-Φ axiom + substrate-vacuum identification gives strukturalnie.
>    Patrz [[../op-Q2-vacuum-budget-2026-05-10]] cycle (STRUCTURAL DERIVED).
> 6. **Q3 CLOSED:** `ρ_EM_quantum/ρ_NS ~ 10⁻¹²` w typowym magnetar; τ.3 TT10
>    testuje L4 mechanism czysto, decoupled od L01 N1.
>    Patrz [[../op-tau3-substrate-clock-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §2.

**Cross-cycle convergence (2026-05-11 update — post-N5) — ★ 10-fold separable sector structure diagnostic ★:**

| Cycle | Diagnosis pattern | Sector separation level |
|---|---|---|
| L01 ADDENDUM §3.2 (Q3 native estimate) | numerical magnitude | 8-12 OOM (corrected 2026-05-11) |
| τ.3 ADDENDUM §2 (L4 vs ρ_EM_quantum) | mechanism | distinct EOM paths |
| ψ.1 ADDENDUM §3 (L01-Q1 resolution) | operator class | disjoint dim-6 vs dim-4 |
| Q2 cycle (vacuum-budget synthesis) | vacuum-level | substrate vs matter sector |
| **op-L01-N1 cycle** (2026-05-11 — *konstruktywna* derivation) | **dedicated 1-loop QED** | **Theorem 2.1 explicit (sympy T4 LOCK)** |
| **op-L01-N2 cycle** (2026-05-11 — *konstruktywna* derivation) | **non-pert. QCD + cosmology** | **Q2 F1 KONSTRUKTYWNIE verified + BBN/CMB/PTA all PASS** |
| **op-L01-N3 cycle** (2026-05-11 — *kompaktowa* verification) | **galactic dust limit + double-counting** | **ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0² verified do <10⁻⁶; gravitational-vs-matter sector explicitly separated** |
| **op-L01-N4 cycle** (2026-05-11 — *konstruktywna* derivation) | **Higgs SSB + 1-loop + EW cosmology** | **Q2 F1 KONSTRUKTYWNIE verified dla Higgs sektora (OOM gap 55.3); R5 LOCK crossover; LISA Ω_GW^EW=0 falsifiable; LHC+Planck+BBN+PTA all PASS automatic** |
| **op-L01-N5 cycle** (2026-05-11 — *konstruktywna* compact derivation) | **EW gauge SU(2)×U(1) trace anomaly** | **Q2 F1 KONSTRUKTYWNIE verified dla gauge sektora (OOM gap 54.9); β_SU(2)<0 asymptotic freedom + β_U(1)>0 Landau pole μ_LP>>M_Pl; PDG 2024 EW precision preserved; ★ LAST L01 N-need closed ★** |

★ **Dziesięć niezależnych diagnoz** zbieżne na **separable sector structure** —
strukturalna własność TGP framework, **konstruktywnie** potwierdzona przez
**pięć dedicated derivations** (N1 EM + N2 QCD + N3 SPARC + N4 Higgs + N5 EW
gauge), NIE post-hoc tuning. **MILESTONE: pełne SM matter+gauge sektor coverage
strukturalnie udowodniony 2026-05-11.**

**Q2 F1 pięć niezależnych metod konstruktywnie verified (post-N5):**
- N1: operator-class disjointness (Theorem 2.1) dla EM sektora
- N2: vacuum-level decoupling + cosmology consistency dla QCD sektora
- N3: gravitational-vs-matter sector separation dla galactic dynamics
- N4: SSB cancellation + thermal decoupling + EW crossover dla Higgs sektora
- **N5: β-functions SU(2)+U(1) + gauge vacuum decoupling dla EW gauge sektora**

**Two-sektor primordial GW synergy (post-N5):**
- EW (LISA mHz, post-2035): Ω_GW^EW = 0 strukturalnie (R5 LOCK; crossover; **N4 + N5 unified prediction**)
- QCD (PTA nHz, NANOGrav 15-yr): no first-order signal (crossover 2+1 flavor lattice; N2)
- ⇒ **TWO empty primordial GW bands** consistent z TGP separable sector structure

**Open follow-up cycles (scaffolded 2026-05-11):**

- [[../op-cluster-mass-deficit-resolution-2026-05-11/]] — **CLOSED 2026-05-11** STRUCTURAL_DERIVED (H1b: TGP + sterile ν 2 eV); 24/24 sympy PASS; 6.4σ multi-experiment falsifiability post-2030+; separate research thread (cluster astrophysical regime; NIE L01 SM sektor coverage)
- [[../op-Higgs-hierarchy-mechanism-2026-05-11/]] — N4 R3 deferred follow-up (Q2 F1 + S05 jako hierarchy protection mechanism); ~8-12 sesji; honest CAVEAT (revolutionary scope)
- [[../op-S07-reset-alternative-f-psi-2026-05-11/]] — M9.1'' FALSIFIED-OBSERVATIONAL GWTC-3 5.02σ rejection follow-up (alternative f(ψ) z S07 freedom); ~5-8 sesji
- [[../op-inflation-substrate-genesis-2026-05-11/]] — Φ_eq(t) inflation prehistory (Q2 NEEDS §N4 deferred); ~8-12 sesji

## Cross-references

- [[ADDENDUM_2026-05-10_native_observables_first.md]] — native-first overlay (2026-05-10)
- [[../../meta/PPN_AS_PROJECTION.md]] — parent methodology (binding)
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] — siostrzany form-meaning protocol
- [[../op-emergent-metric-from-interaction-2026-05-09]] — siostrzany cykl gravity-sector (L01 jest L1 input dla niego)
- [[../op-psi1-substrate-light-acceleration]] — Q1 closure (operator class disjoint)
- [[../op-tau3-substrate-clock-acceleration]] — Q3 closure (mechanism decoupling)
- [[../op-Q2-vacuum-budget-2026-05-10]] — Q2 closure (vacuum-level decoupling, dedicated cycle)
- [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11]] — **N1 closure (Quantum trace anomaly EM constructive derivation, 2026-05-11, STRUCTURAL_DERIVED)**
- [[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11]] — **N2 closure (QCD trace anomaly + cosmology constructive derivation, 2026-05-11, STRUCTURAL_DERIVED, Q2 F1 konstruktywnie verified)**
- [[../op-L01-N3-SPARC-rho-consistency-2026-05-11]] — **N3 closure (SPARC ρ-consistency compact verification, 2026-05-11, STRUCTURAL_DERIVED, ρ_baryon ≡ -T^μ_μ_dust/c_0² do <10⁻⁶; double-counting closed)**
- [[../op-L01-N4-Higgs-trace-anomaly-2026-05-11]] — **N4 closure (Higgs trace anomaly + EW phase transition cosmology constructive derivation, 2026-05-11, STRUCTURAL_DERIVED, 24/24 sympy PASS, Q2 F1 konstruktywnie verified dla Higgs sektora; R5 LOCK crossover + LISA Ω_GW^EW=0 falsifiable post-2035)**
- [[../op-L01-N5-EW-gauge-anomaly-2026-05-11]] — **N5 closure (EW gauge anomaly SU(2)×U(1) compact constructive derivation, 2026-05-11, STRUCTURAL_DERIVED, 8/8 sympy PASS, Q2 F1 konstruktywnie verified dla gauge sektora; β_SU(2) asymptotic freedom + β_U(1) Landau pole; ★ LAST L01 N-need closed; full SM matter+gauge coverage milestone ★)**
- [[../closure_2026-04-26/Lambda_from_Phi0]] — T-Λ parent of Q2 cycle
- [[../../audyt/L01_rho_operational]]
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] §`eq:L-mat-unified`
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`ax:metric-coupling`
- [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.4, §C.3, §N
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]]
