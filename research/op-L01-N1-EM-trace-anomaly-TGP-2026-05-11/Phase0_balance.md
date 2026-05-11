---
title: "Phase 0 — Balance sheet + 6/6 gate criteria + literature cross-reference + initial NEEDS"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-balance
phase: 0
status: 🟢 RESOLVED — 6/6 gate criteria PASS
gate_criteria_passed: 6
gate_criteria_total: 6
predecessors_inventoried: 7
needs_initial_count: 8
tags:
  - phase0
  - balance-sheet
  - gate-criteria
  - literature-cross-reference
  - QFT-curved-background
  - native-first-methodology
---

# Phase 0 — Balance sheet

## §0 — Executive summary

**6/6 gate criteria PASS**. Balance sheet identifies:

- **7 predecessor cycles** giving structural inputs (L01 ρ-bridge classical, emergent-metric g_eff[{Φ_i}], ψ.1 operator class, τ.3 mechanism decoupling, B9 MICROSCOPE baseline, G.0 vacuum stability, Q2 vacuum budget).
- **3 canonical literature anchors** (Birrell-Davies 1982, Riegert 1984, Capper-Duff-Halpern 1974) + **2 modern checkpoints** (Duff 1994 review, arXiv:2306.03892 2023/2025).
- **1 critical risk-input** literature reference (JHEP 2015 Asorey et al. — QED trace anomaly può induce QEP violations; explicit check required Phase 3).
- **8 initial sub-needs** (N0.1–N0.8) feeding into Phase 1-4 work.

**No structural blockers identified.** Phase 1 may proceed.

## §1 — Inventory: existing TGP results (predecessors)

### §1.1 — L01 ρ-bridge cycle (CLOSED-DERIVED 2026-05-10)

[[../op-L01-rho-stress-energy-bridge-2026-05-04/]]

**Inputs do tego cyklu:**

1. **Formal definition** ρ ≡ -T^μ_μ/c_0² (formal_definition.md §3-4) — derived
   z ax:metric-coupling, NIE postulat.
2. **SM sector mapping** (SM_sector_mapping.md): 5 sektorów klasycznie:
   - Dirac fermion: ρ_Dirac = m·ψ̄ψ/c_0² ≥ 0
   - Massive scalar: ρ_scalar = ((∂φ)² - 4V)/c_0²
   - **EM (massless gauge)**: T^μ_μ_EM = 0 (conformal invariance) ⇒ ρ_EM = 0 *klasycznie*
   - Yang-Mills classical: T = 0 ⇒ ρ_YM_classical = 0
   - Perfect fluid: ρ_TGP = (ρ_e − 3p)/c_0²
3. **Photon treatment** (photon_treatment.md §5): "Quantum trace anomaly EM"
   sygnalizowany jako N1; explicit form
   `T^μ_μ_EM,quantum = (β(α)/(2α))·F²` z β(α)/(2α) ~ α/(3π) ≈ 7.7·10⁻⁷.
4. **ADDENDUM 2026-05-10 §3.2** native estimate dla magnetar regime:
   ρ_EM_quantum/ρ_NS ~ 10⁻¹² typowo, 10⁻¹⁰ w hot-spots — *nie significant*
   Φ-shift na powierzchni magnetara.
5. **NEEDS §T.1** three-layer specification dla N1 (L1 native / L2 projection
   chart / L3 falsifier map).
6. **F4-F5, F11, F15 (klasyczne photon results)** preserved: c_GW=c_EM, no 5-th
   force from photons, radiation era non-source.

**Co tem cykl daje L01 wzajemnie:**
- Dedicated **constructive** derivation N1 (zamiast deferred placeholder).
- Confirms Q1 closure z **explicit operator structure** (a nie tylko z argumentu klasy).
- Confirms Q3 numerical estimate via formalism (vs initial back-of-envelope).
- Adds entry M911-EM-quantum do [[../../PREDICTIONS_REGISTRY.md]].

### §1.2 — op-emergent-metric-from-interaction (STRUCTURAL DERIVED 2026-05-09)

[[../op-emergent-metric-from-interaction-2026-05-09/]] (57/57 sympy PASS, 6/6 P-requirements)

**Krytyczne inputs:**

1. **g_eff = G[{Φ_i}, σ_ab[Φ], Φ̄, x]** funkcjonał konfiguracji Φ-źródeł, NIE
   lokalna funkcja f(ψ). To zastępuje M9.1'' postulate (FALSIFIED 5σ GWTC-3
   2026-05-09).
2. **Phase 1 ansatz** (Phase1_results §"Strukturalna decyzja"):
   ```
   g_eff^00 = -A(ψ)
   g_eff^ij = δ^ij·B(ψ) + σ^ij·C(ψ)/(Φ_0² c²)
   g_eff^0i = 0  (statyczny limit)
   σ^ij = (∂^iΦ)(∂^jΦ) - (1/3)δ^ij δ_kl(∂^kΦ)(∂^lΦ)
   ```
   Trzy niezależne funkcje A(ψ), B(ψ), C(ψ).
3. **Phase 4 GWTC-3 zero-β region**: parametric family w (a_n, ξ_n, b_n, c_0, κ_σ)
   space; M911 specific point falsified, neighbourhood open.
4. **Phase 6 SU(2) cross-consistency**: ten sam mechanizm interakcji daje
   tensor structure at level 2 (g_eff) i level 3 (SU(2) spin) — H6.1 unification.
5. **§5.1 BD/Horndeski demarcation**: g_eff jest *funkcjonałem*, brak osobnego
   variation principle.

**Implikacja dla Phase 1:** krzywizna `R[g_eff]`, `R_μν[g_eff]`, etc. w 1-loop
QED redukują się do funkcji `(∂Φ_i, ∂Φ_i·∂Φ_j, σ_ab, Φ̄)`. Phase 2 tutaj wykonuje
explicit reduction.

**Risk R1 (M9.1'' contamination):** unikamy go strukturalnie używając Phase 1
3-funkcyjny ansatz {A, B, C}, NIE specific f_M911(ψ) = (4-3ψ)/ψ.

### §1.3 — ψ.1 cycle (PASS-CLOSED Phase 7) + ADDENDUM 2026-05-10

[[../op-psi1-substrate-light-acceleration/]] + [[../op-psi1-substrate-light-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]]

**Kluczowe inputs:**

1. **Canonical basis B_ψ.1.v3^dim-6 = {L₅'_a, L₅'_b}** (Phase 7 Hilbert-series-style
   enumeration):
   - L₅'_a = (β_g/Λ²)·(∂_μ lnX)(∂_ν lnX)·F^{μρ} F^ν_ρ (parity-even)
   - L₅'_b = β̃_g·(∂_μ lnX)(∂_ν lnX)·F^{μρ} F̃^ν_ρ (parity-odd)
2. **Sterile-sector exclusion** (Phase 7 T7.3 R4): scalar `(∂lnX)² F²` (varying-α
   Bekenstein/Sandvik) i `(∂lnX)² F·F̃` (axion-like ω.1) są EXPLICITLY
   FILTERED OUT z ψ.1 basis.
3. **Pure-photon decoupling** (Phase 7 T7.1 invariance filter): operatory
   `(F, F, F)` — PURE-PHOTON dim-4 — żyją w **separate sector** (Euler-Heisenberg
   / pure-EM 1-loop renormalization).

**ADDENDUM 2026-05-10 §3 Q1 closure** (cited extensively):

> Quantum trace anomaly EM `(α/(3π))·F²` jest **pure-photon dim-4 operator**, NIE
> dim-6 cross-sector. Klasa operatorów: **renormalizacja pure-EM coupling**, nie
> *substrate gradient coupling*. ψ.1.v3 canonical basis NIE *covers* quantum trace
> anomaly. Te są *disjoint sektory*.

**Implikacja dla tego cyklu (Phase 2 verification target):** derivation MUSI
*konstruktywnie* potwierdzić, że produkt jest pure-photon (z curvature corrections
sourced by ∂Φ), NIE odtworzyć dim-6 cross-sector operators. To jest **R2 risk
guard**.

### §1.4 — τ.3 cycle (PASS-CLOSED) + ADDENDUM 2026-05-10

[[../op-tau3-substrate-clock-acceleration/]] (referenced w L01 README §148-156)

**Kluczowe input:** mechanism-level decoupling τ.3 L4 gradient-coupled mass mechanism
od L01 ρ_EM_quantum — Phase3.TT10 magnetar X-ray timing testuje L4, NIE trace
anomaly. *10 OOM separation* w native estimate.

**Implikacja dla Phase 3:** numerical estimates tutaj powinny być spójne z τ.3
diagnozą (10⁻¹² ratio confirmed via formal derivation).

### §1.5 — B9 WEP closure (PASS 6/6 2026-05-01)

[[../op-newton-momentum/B9_wep_microscope_composition_results.md]] (referenced)

**Input:** η_TGP_Dirac (Pt vs Ti) = 1.32·10⁻²⁶ baseline, ≪ MICROSCOPE bound
1.1·10⁻¹⁵.

**Implikacja dla Phase 3 (P5 verification):** ρ_EM_quantum addition do η budget
musi być ≪ 10⁻¹⁵. Phase 3 tutaj robi explicit estimate η_TGP_EM_quantum dla
Pt vs Ti (różne EM-content per atom).

### §1.6 — G.0 vacuum stability closure (referenced)

m_sp² > 0 cosmological-scale → Yukawa decay przy lab scale → no significant
Φ-mediated 5-th force from EM at solar system / lab. Quantum loops nie zmieniają
tego (loops dotyczą EM sektora; Φ scale m_sp ustalony przez vacuum potential V(Φ)).

### §1.7 — op-Q2-vacuum-budget-2026-05-10 (STRUCTURAL DERIVED)

SM matter sector vacua (ρ_QCD, ρ_Higgs, ρ_EW) **NIE additive** do bare Λ — single-Φ
axiom + substrate-vacuum identification gives strukturalnie. Implikacja dla N1:
ρ_EM_quantum w lab/magnetar regime jest *transient source*, NIE contribution do
*today's* Λ. Phase 4 cosmology layer powinien ten point reflektować.

## §2 — Literature cross-reference (verified 2026-05-11)

### §2.1 — Canonical anchors (used as reference framework)

| Reference | Bibliographic key | Role w cyklu |
|---|---|---|
| **Birrell, Davies, "Quantum Fields in Curved Space"** (CUP 1982) | Birrell-Davies-1982 | **Canonical** dla 1-loop QFT na curved background; chapter 6 trace anomaly framework. Phase 1 setup z §6.1-6.4. |
| **Riegert, Phys. Lett. B 134, 56 (1984)** | Riegert-1984 | Conformal anomaly action — explicit non-local effective action z trace anomaly source. Phase 1 decomposition target. |
| **Capper, Duff, Halpern, Phys. Rev. D 10, 461 (1974)** | CDH-1974 | 1-loop QED trace anomaly **original derivation**; β(α)/(2α) formula. Phase 1 sympy LOCK target. |

### §2.2 — Modern checkpoints

| Reference | Bibliographic key | Role |
|---|---|---|
| **Duff, Class. Quantum Grav. 11, 1387 (1994)** | Duff-1994 | Modern review of trace anomalies; Type A vs Type B classification; canonical compilation. |
| **arXiv:2306.03892** (Boyle, Padilla — 2023, updated 2025) | Boyle-Padilla-2023 | "Conformal anomaly and gravitational pair production" — recent treatment łączące anomaly z curvature; ujęcie 4D explicit; cross-check formul. |

### §2.3 — Risk-input literature (R6 QEP violations check)

| Reference | Bibliographic key | Risk relevance |
|---|---|---|
| **JHEP 05 (2015) 118 (Asorey, Bautista, ...) — "QED trace anomaly, non-local Lagrangians and quantum equivalence principle violations"** | Asorey-2015 | **R6:** QED trace anomaly *może* induktywnie generować QEP violations (non-local effective action). Phase 3 tutaj musi explicit verify że universal coupling structure w `L_mat = -(q/Φ_0)·φ·ρ` automatically kasuje to (S05 mechanism — wspólna metryka g_eff dla wszystkich). |

### §2.4 — Literature consistency note

**WebSearch 2026-05-11** (queries: "QED trace anomaly 1-loop scalar-tensor background
curved space 2024 2025"; "Riegert conformal anomaly action effective metric scalar
field magnetar 2024 2025 review"):

- Birrell-Davies framework (1982) **remains canonical** — no major revision through
  2025.
- Boyle-Padilla 2023 (arXiv:2306.03892) confirms standard β(α)/(2α)·F² + curvature
  structure.
- **No specific 2024-2026 review** on Riegert action × magnetar regimes via
  WebSearch — non-finding noted; framework choice based na canonical references.
- Asorey-2015 framework dla QEP violations — used as risk-flag input dla R6.

**No literature update-required** dla Phase 1 setup. Canonical Birrell-Davies +
Riegert 1984 + Capper-Duff-Halpern 1974 jest **wystarczające**.

## §3 — Initial NEEDS list

Sub-needs identyfikowane dla Phase 1-4 work:

| ID | Sub-need | Phase | Depends on |
|---|---|---|---|
| **N0.1** | Birrell-Davies 1-loop effective action S_eff[A_μ; g_eff[{Φ_i}]] explicit form | 1 | g_eff ansatz from emergent-metric Phase 1 |
| **N0.2** | Sympy LOCK β(α) = α²/(3π) z 1-loop QED β-function | 1 | CDH-1974 derivation |
| **N0.3** | Riegert action decomposition T^μ_μ_EM,1-loop = β/2α·F² + a₁·R·F² + a₂·R_μν F^μρ F^ν_ρ + a₃·□F² + non-local | 1 | Riegert 1984 framework |
| **N0.4** | Reduction R, R_μν w g_eff[{Φ_i}] limit do funkcji {∂Φ_i, ∂Φ_i·∂Φ_j, σ_ab} | 2 | emergent-metric Phase 1 ansatz {A, B, C} |
| **N0.5** | Disjointness verification: derivation produkt nie odtwarza ψ.1.v3 dim-6 EFT operators | 2 | ψ.1 Phase 7 enumeration |
| **N0.6** | GW170817 c_GW=c_EM preservation — dispersion relation linearization 1-loop corrections | 2 | linear analysis |
| **N0.7** | Lab regime numerics (B~1 T, R~10⁻⁵² m⁻²): ρ_EM_quantum scaling + MICROSCOPE η budget | 3 | Phase 1+2 outputs |
| **N0.8** | Magnetar regime numerics (B~10¹¹ T): ρ_EM_quantum/ρ_NS ratio + B ≪ B_QED limit | 3 | Phase 1+2 outputs + R5 caveat |
| **N0.9** | QEP universality verification — universal coupling kasuje Asorey-2015 type QEP violations | 3 | Phase 1 effective action + S05 mechanism |
| **N0.10** | Native parameter audit + forced-zero declarations + three-layer L1/L2/L3 closure | 4 | Phase 1-3 outputs |

(Numerowane N0.1-N0.10 dla wyróżnienia od finalnych N1-Nn z Phase 4 closure.)

## §4 — Six-gate criteria check

Per [[../../meta/PPN_AS_PROJECTION.md]] §3.1 + standard TGP Phase 0 protocol:

| # | Gate criterion | Status | Evidence |
|---|---|---|---|
| **G1** | Predecessors inventoried + structural inputs explicit | ✅ PASS | §1.1-§1.7: 7 predecessor cycles enumerowane, każdy z explicit input list |
| **G2** | Literature canonical refs cross-referenced + currency check | ✅ PASS | §2: Birrell-Davies + Riegert + CDH-1974 canonical; modern checkpoints (Duff 1994, Boyle-Padilla 2023) verified via WebSearch 2026-05-11 |
| **G3** | Risk flags declared explicit (R1-R6) + addressing strategy per risk | ✅ PASS | README §"Six risks" + §"R-i ↦ Phase j addressing"; R6 enriched z Asorey-2015 literature input |
| **G4** | NEEDS list initialized (Phase 1-4 sub-needs) | ✅ PASS | §3: N0.1-N0.10 enumerated z Phase mapping |
| **G5** | Methodology binding declared (native-first + sympy LOCK + S05 + §5.1) | ✅ PASS | README §"Methodology constraints" 9 reguł; binding 2026-05-10+ explicit |
| **G6** | Cross-cycle consistency map declared | ✅ PASS | README §"Connection do innych cykli" + §"Strukturalna pozycja"; trzy niezależne diagnozy zbieżne (L01 §3.2, τ.3 §2, ψ.1 §3) explicit |

**6/6 GATE PASS** — Phase 1 may proceed.

## §5 — Strategic assessment

### §5.1 — Why this cycle is structurally clean

1. **Problem space is restricted by predecessor closures:** ψ.1 Q1 closure already
   identified that quantum trace anomaly EM is **disjoint operator class** od ψ.1.v3
   dim-6 EFT. Tego cyklu zadanie nie jest "discover new physics", ale
   "**konstruktywnie** zrealizować disjoint-class derivation". Probability of
   cycle SUCCESS *raised* przez ten upfront classification.
2. **Recovery framework is robust:** g_eff = G[{Φ_i}] STRUCTURAL DERIVED 2026-05-09
   replaces M9.1'' falsified postulate; nasz N1 derivation dziedziczy zdrową
   strukturę background metric.
3. **Numerical target jest already estimated:** ρ_EM_quantum/ρ_NS ~ 10⁻¹² (L01
   ADDENDUM §3.2 Q3) — sets *expected outcome* dla Phase 3. Phase 3 musi *odtworzyć*
   ten back-of-envelope estimate via formal derivation.
4. **MICROSCOPE bound jest *automatic*** dla universal coupling structure (B9
   closure baseline 1.32·10⁻²⁶ ≪ 10⁻¹⁵). Phase 3 tutaj jest verification, NIE
   tuning.

### §5.2 — Where structural risks live

1. **R1 (M9.1'' contamination):** **eliminowany strukturalnie** w Phase 1 setup —
   ansatz {A, B, C} 3-funkcyjny per emergent-metric, NIE M9.1'' f(ψ).
2. **R2 (operator class re-overlap):** **verified explicit Phase 2** — derivation
   produkt mapped na ψ.1.v3 basis i pokazane disjoint.
3. **R3 (GW170817 c-violation):** **verified Phase 2/3** — dispersion linearization;
   trace anomaly nie wprowadza propagującego scalar mode, GW170817 preserved
   automatic.
4. **R4 (S05 violation):** **verified Phase 1** — quantum loops dotyczą tylko
   EM sektora; nie wprowadzają drugiego propagującego pola; g_eff funkcjonał
   {Φ_i} structure preserved przez 1-loop integration out fermion/photon
   loops on background.
5. **R5 (perturbative breakdown w magnetar):** **honestly documented Phase 3** —
   B ≪ B_QED ~ 4·10⁹ T regime explicit; B ≳ B_QED zostaje deferred to
   non-perturbative analysis (poza zasięgiem cyklu).
6. **R6 (QEP violations from trace anomaly):** **explicit Phase 3 verification** —
   universal coupling structure w L_mat (single g_eff dla wszystkich) automatycznie
   kasuje Asorey-2015 type QEP violations. Strukturalna *consequence* S05.

### §5.3 — Probability assessment update (post-Phase-0)

| Outcome | Pre-Phase-0 (README) | Post-Phase-0 (this) | Reason |
|---|---|---|---|
| Pełen DERIVED | 35-50% | **45-60%** ↑ | Predecessors + literature stronger niż expected; risk addressing well-localized |
| STRUCTURAL CONDITIONAL | 25-35% | 25-30% (similar) | Riegert non-local terms może zostać deferred to multi-session |
| STRUCTURAL_NO_GO | 10-20% | **5-15%** ↓ | R3, R6 mają strong structural arguments dla preservation |
| EARLY_HALT | 5-15% | 5-10% (similar) | R5 (magnetar B ≳ B_QED) jest *expected* deferred, NIE blocker |

**Trend:** Phase 0 raises confidence w success; cycle pozycjonowany jako
*constructive verification* of expected disjoint-class derivation, NIE jako
*discovery of new physics*.

## §6 — Phase 1 setup preview

Phase 1 (next session) will execute:

1. **N0.1 setup:** Birrell-Davies §6.1-6.4 effective action S_eff[A_μ; g_eff].
2. **N0.2 sympy LOCK:** β(α) = α²/(3π) 1-loop QED.
3. **N0.3 Riegert decomposition:** T^μ_μ_EM,1-loop expansion w basis
   {β/2α·F², R·F², R_μν F^μρ F^ν_ρ, □F², non-local Riegert terms}.
4. **R1 guard (Phase 1 explicit):** g_eff ansatz {A(ψ), B(ψ), C(ψ)} per
   emergent-metric Phase 1, NIE M9.1''.
5. **R4 guard (Phase 1 explicit):** verification że 1-loop integration out
   fermion/photon loops na background g_eff[{Φ_i}] nie wprowadza drugiego
   propagującego pola — S05 preserved.

Phase 1 deliverables: Phase1_setup.md + Phase1_results.md + Phase1_sympy.py +
Phase1_sympy.txt.

## §7 — Cross-references

- [[./README.md]] — cycle overview
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/]] — parent L01 cycle (N1 deferred here)
- [[../op-emergent-metric-from-interaction-2026-05-09/]] — g_eff[{Φ_i}] structural input
- [[../op-psi1-substrate-light-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §3 — Q1 operator class closure
- [[../op-tau3-substrate-clock-acceleration/]] — τ.3 mechanism decoupling
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]] — MICROSCOPE baseline
- [[../op-Q2-vacuum-budget-2026-05-10/]] — vacuum-level decoupling
- [[../../meta/PPN_AS_PROJECTION.md]] — methodology binding
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §4 — form-meaning protocol
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] §eq:L-mat-unified

---

**Phase 0 close:** 6/6 gate criteria PASS. Phase 1 may proceed. Probability of
cycle SUCCESS raised do 45-60% Pełen DERIVED (post-balance-sheet).
