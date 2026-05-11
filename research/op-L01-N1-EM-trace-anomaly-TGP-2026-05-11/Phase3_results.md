---
title: "Phase 3 results — Phenomenology numerics + MICROSCOPE/Eöt-Wash/LLR bound checks + R6 QEP closure"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-results
phase: 3
status: 🟢 RESOLVED — phenomenology consistent z bounds; R5 documented; R6 closed
sub_needs_resolved: [N0.7, N0.8, N0.9]
risks_addressed: [R5 documented, R6 closed]
predecessor: "[[./Phase2_results.md]] (8/8 PASS)"
tags:
  - phase3
  - phenomenology-results
  - lab-magnetar-numerics
  - microscope-bounded
  - QEP-universal-immunity
  - L01-typo-correction-propagated
---

# Phase 3 results

## §0 — Executive summary

**Phenomenology consistent z all observational bounds.** Phase 3 establishes:

1. **Lab regime (B=1 T):** ρ_EM_quantum_lab ~ 7·10⁻¹⁵ kg/m³ — local *inside*
   electromagnet; nieistotne dla gravitational dynamics solar system.
2. **Magnetar regime (B=10¹¹ T, extreme):** ρ_EM_quantum/ρ_NS ~ **1.7·10⁻¹⁰** z
   *corrected* α/(3π)≈7.74·10⁻⁴ prefactor (NIE 10⁻¹² jak L01 ADDENDUM cytuje
   z typo). 2 OOM korekta propagated to Phase 4 typo correction list.
3. **MICROSCOPE η bound:** η_TGP_EM_quantum (Pt vs Ti) = **0 strukturalnie** z
   universal coupling structure (NIE numerically suppressed; *exactly zero*).
4. **Eöt-Wash, LLR Nordtvedt bounds:** wszystkie *automatic* z universal coupling.
5. **R5 (perturbative breakdown):** documented — cycle valid B ≪ B_QED ≈ 4.4·10⁹ T;
   B ≳ B_QED extreme magnetar regions deferred do non-perturbative analysis.
6. **R6 (QEP violations from Asorey-2015 type non-local Lagrangians):**
   **structurally closed** — TGP universal coupling z S05 immune do takiego
   QEP violation mechanism.
7. **TT10 magnetar test status preserved:** τ.3 mechanism decoupling intact
   (8-10 OOM separation vs L01 ρ_EM_quantum).

## §1 — Lab regime numerics (B = 1 T)

### §1.1 — Energy density + ρ_EM_quantum

```
B = 1 T  →  u_mag = B²/(2μ_0) = 1/(2·4π·10⁻⁷) ≈ 4·10⁵ J/m³
F² = 2·(B² - E²/c²) ≈ 2 B² in SI · (μ_0/μ_0)  →  F²·F² density = 2·u_mag
T^μ_μ_EM,quantum_lab = (α/(3π)) · F² ≈ 7.74·10⁻⁴ · 2 · 4·10⁵
                      ≈ 6.2·10² J/m³
ρ_EM_quantum_lab = -T^μ_μ / c_0² ≈ -6.2·10² / (9·10¹⁶) ≈ -7·10⁻¹⁵ kg/m³
```

|ρ_EM_quantum_lab| ~ **7·10⁻¹⁵ kg/m³** (B=1 T inside magnet volume).

### §1.2 — Implikacje (lab)

- Local energy density inside electromagnet: ~7·10⁻¹⁵ kg/m³, ~6 OOM więcej
  niż interstellar ρ_baryon (10⁻²¹ kg/m³).
- **Φ-shift kontrybucja** w obrębie elektromagnetu: integrating across magnet
  volume V_magnet ~ (1 m)³ → m_eq ~ 7·10⁻¹⁵ kg.
- Φ-mediated 5-th force from this mass: G·m_eq/r² ~ 6.67·10⁻¹¹ · 7·10⁻¹⁵ / r²
  ≈ 5·10⁻²⁵ / r² m/s² (for r~m: 5·10⁻²⁵ m/s², ABSURDNIE małe).

⇒ **Niewykrywalne** z lab acceleration measurements (current sensitivity ~10⁻¹⁵
m/s² with atom interferometry, far below 10⁻²⁵).

## §2 — Magnetar regime numerics (B = 10⁹ to 10¹¹ T)

### §2.1 — Estimate dla B=10¹⁰ T (typical magnetar)

```
B = 10¹⁰ T  →  u_mag = 10²⁰ / 2.51·10⁻⁶ ≈ 4·10²⁵ J/m³
T^μ_μ_quantum = 7.74·10⁻⁴ · 2 · 4·10²⁵ ≈ 6·10²² J/m³
ρ_EM_quantum_magnetar ≈ -6·10²² / 9·10¹⁶ ≈ -7·10⁵ kg/m³
ρ_EM_quantum/ρ_NS ~ 7·10⁵ / 4·10¹⁷ ≈ 1.7·10⁻¹²
```

### §2.2 — Estimate dla B=10¹¹ T (extreme magnetar)

```
B = 10¹¹ T  →  u_mag ≈ 4·10²⁷ J/m³
T^μ_μ_quantum ≈ 7.74·10⁻⁴ · 8·10²⁷ ≈ 6·10²⁴ J/m³
ρ_EM_quantum_magnetar ≈ -7·10⁷ kg/m³
ρ_EM_quantum/ρ_NS ~ 1.7·10⁻¹⁰
```

**WAŻNE:** B/B_QED ≈ 23 dla B=10¹¹ T → **non-perturbative regime**, linear
extrapolation jest *upper bound estimate*. True non-perturbative result mógłby być
kilka rzędów wielkości mniejszy (przez Schwinger pair production back-reaction
i higher-order QED corrections).

### §2.3 — L01 ADDENDUM §3.2 Q3 typo correction

L01 ADDENDUM cytuje:
> ρ_EM_quantum/ρ_NS ~ 10⁻¹² w typowym magnetar (B=10¹¹ T)

z calculation:
- B²/(2μ_0) ~ 4·10²⁸ J/m³ (close to my 4·10²⁷, factor 10 due to definition)
- α²/(3π) ≈ 7.7·10⁻⁷  ← **typo**, should be α/(3π) ≈ 7.74·10⁻⁴

**Corrected ratio:** ~10⁻¹⁰ dla B=10¹¹ T (z proper α/(3π) prefactor).

**Action item Phase 4:** propagate correction do:
- L01 ADDENDUM §3.2 Q3
- L01 NEEDS.md §T.1 (gdzie typo jest cytowane)
- τ.3 ADDENDUM §2 (gdzie cross-cycle convergence diagnostic używa "10 OOM
  separation"; correction: **8 OOM separation** dla B=10¹¹ T extreme magnetar;
  10 OOM separation dla B=10¹⁰ T typical magnetar)

### §2.4 — Konsekwencja: TT10 status

τ.3 ADDENDUM §2 conclusion: TT10 magnetar X-ray timing test **insensitive** to
ρ_EM_quantum dla typical magnetar (10⁻¹² ratio) — testuje L4 gradient-coupled mass
mechanism czysto.

Po typo correction:
- typical magnetar B=10¹⁰ T: ratio ~1.7·10⁻¹², separation 12 OOM (similar)
- extreme magnetar B=10¹¹ T: ratio ~1.7·10⁻¹⁰, separation 10 OOM (still strong)

⇒ **TT10 mechanism decoupling preserved**. Trace anomaly contribution wciąż
≪ surface ρ_NS, NIE psuje L4 mechanism test.

## §3 — MICROSCOPE / Eöt-Wash / LLR bound checks

### §3.1 — η_TGP_EM_quantum (Pt vs Ti) — strukturalnie zero

**Universal coupling structure** w TGP:
```
L_mat = -(q/Φ_0) · φ · ρ
```

ρ ≡ -T^μ_μ/c_0² zbiera **wszystkie** contributions z atomic stress-energy:
1. rest mass elektronów + protonów + neutronów
2. nuclear binding energy
3. atomic EM binding energy
4. Lamb shifts + radiative corrections
5. **quantum trace anomaly contributions**

Wszystkie te contributions enter z **tym samym** prefactor `q/Φ_0`. Brak
sektor-specific coupling.

**Konsekwencja:** Pt atom z N_Pt elektronów ma swój ρ_Pt ≈ M_Pt/V_Pt z
*kompletnym* atomic mass-energy budget (włączając quantum corrections).
Ti atom analogicznie ma ρ_Ti = M_Ti/V_Ti.

Φ-mediated acceleration:
```
a_Pt = -(q/Φ_0) · ∇φ · (ρ_Pt/M_Pt) = -(q/Φ_0)·∇φ·(1/V_Pt)·...   for Pt-block
```

Dla *atom-level* WEP test (Pt atom vs Ti atom w gravitational field), free-fall
acceleration:
```
a_atom = -∇U_grav + (Φ-correction)·M_atom
```

z (Φ-correction) **uniwersalną** dla wszystkich atomów. Zatem:
```
a_Pt - a_Ti = 0 + 0 + 0 = 0    [structurally]
```

η_TGP_EM_quantum (Pt vs Ti) = **0 strukturalnie**.

### §3.2 — MICROSCOPE bound check

Bound: η_MICROSCOPE ≤ 1.1·10⁻¹⁵ (Pt vs Ti, 2017; tighter 2027+ projection).

TGP prediction: η_TGP_total = η_TGP_Dirac + η_TGP_EM_quantum + ... = 1.32·10⁻²⁶
+ 0 + ... ≈ **1.32·10⁻²⁶** (B9 baseline preserved).

**Margin:** ~11 OOM below MICROSCOPE bound. **PASS.**

Quantum trace anomaly **nie zmienia** B9 closure verdict.

### §3.3 — Eöt-Wash bound check

Eöt-Wash η ≤ 5·10⁻¹⁴ (terrestrial torsion balance, various test masses).

Same structural argument: universal coupling → η_TGP_EM_quantum = 0 →
**Eöt-Wash automatic PASS**.

### §3.4 — LLR Nordtvedt bound check

Lunar Laser Ranging η_N ≤ 4.4·10⁻⁴ (Earth vs Moon as different composition).

Same structural argument: universal coupling → η_N_TGP_EM_quantum = 0 →
**LLR automatic PASS**.

## §4 — R6 closure: QEP universality verification

### §4.1 — Asorey-2015 mechanism

JHEP 05 (2015) 118: "QED trace anomaly, non-local Lagrangians and quantum
equivalence principle violations".

Argument:
- 1-loop QED trace anomaly produces *non-local* effective action (Riegert action).
- W teoriach z **non-universal** coupling do gravity (np. Brans-Dicke z
  matter-dependent dilaton coupling), non-local Lagrangians can break QEP.
- Concrete falsifier: QEP-violating contributions can be bounded by MICROSCOPE.

### §4.2 — TGP immunity argument

**Strukturalna własność TGP** (S05 + ax:metric-coupling):

1. **S05** gwarantuje *single fundamental field substratu* (Φ).
2. **ax:metric-coupling** nakazuje, że materia sprzęga się **wyłącznie** przez
   `g_eff[Φ]` minimal coupling.
3. **L_mat = -(q/Φ_0)·φ·ρ** z `ρ ≡ -T^μ_μ/c_0²` — universal across all matter
   sectors.
4. **Riegert localization** (Phase 1 §2.2): σ_eff = funkcja Φ, NIE niezależny
   field.

⇒ Quantum trace anomaly w TGP NIE generuje non-universal coupling. Wszystkie
QED corrections (włączając non-local Riegert action) widzą **ten sam Φ-coupling**
przez `g_eff[Φ]`.

**Konsekwencja:** Asorey-2015 type QEP violations **nie istnieją** strukturalnie
w TGP. R6 risk closed.

### §4.3 — Konsekwencja dla MICROSCOPE constraint interpretation

W standard QFT-on-curved-space framework (Asorey 2015): MICROSCOPE η bound
constrains *non-universal coupling parameters*.

W TGP: MICROSCOPE constrains *η_TGP_total* z all sectors, ale TGP framework
*structurally* gives η = 1.32·10⁻²⁶ z B9 closure (mass-coupling structure),
unaffected by quantum corrections.

⇒ **MICROSCOPE PASS automatic** dla TGP, *independent* of post-2017 measurement
precision.

## §5 — R5 documentation: B ≪ B_QED limitation

### §5.1 — Schwinger critical field

```
B_QED = m_e² c² / (e ℏ) ≈ 4.41·10⁹ T
E_QED = m_e² c³ / (e ℏ) ≈ 1.32·10¹⁸ V/m
```

For B ≳ B_QED:
- Pair production (Schwinger effect) becomes appreciable
- Vacuum birefringence + photon splitting become observable
- Heisenberg-Euler effective Lagrangian (4-loop+ QED) needed
- Linear extrapolation z 1-loop fails

### §5.2 — Cycle scope

**Cycle valid:**
- All lab conditions (B ≪ 10⁹ T)
- Most magnetar interior + neutron-star surface (B ≲ B_QED locally)
- Cosmological + galactic (B ≪ B_QED everywhere)

**Cycle deferred (non-perturbative analysis):**
- Magnetar atmosphere/surface hot-spots (B ~ 10¹¹ T → B/B_QED ≈ 23)
- Hypothetical Schwinger-class lab fields (E ~ 10¹⁸ V/m, B ~ 10⁹ T macroscopic;
  beyond 2030+ technology)

### §5.3 — Honest caveat

**Phase 3 numerical estimates dla B=10¹¹ T są upper bound** — true
non-perturbative result mógłby być znacznie mniejszy (cancellation between
positive and negative contributions, vacuum back-reaction).

**Action item:** future cycle `op-EM-trace-anomaly-Schwinger-extension` (z
non-perturbative QED) byłby continuation tego cyklu dla extreme magnetar regimes.
**NIE blocker** dla obecnego N1 closure — R5 honestly documented.

## §6 — Findings (exportable)

| ID | Finding | Source |
|---|---|---|
| **F3.1** | Lab regime ρ_EM_quantum (B=1 T): ~7·10⁻¹⁵ kg/m³ inside magnet; niewykrywalne acceleration | §1.1 |
| **F3.2** | Magnetar typical (B=10¹⁰ T): ρ_EM_quantum/ρ_NS ~ 1.7·10⁻¹² | §2.1 |
| **F3.3** | Magnetar extreme (B=10¹¹ T, B/B_QED~23): ratio ~1.7·10⁻¹⁰ z corrected α/(3π) | §2.2 |
| **F3.4** | L01 ADDENDUM §3.2 Q3 typo: "α²/(3π)≈7.7·10⁻⁷" → corrected α/(3π)≈7.74·10⁻⁴; 2-3 OOM correction | §2.3, F1.3 |
| **F3.5** | η_TGP_EM_quantum (Pt vs Ti) = **0 strukturalnie** z universal coupling structure (S05) | §3.1 |
| **F3.6** | MICROSCOPE bound automatic PASS: η_TGP = 1.32·10⁻²⁶ (B9) ≪ 10⁻¹⁵, ~11 OOM margin | §3.2 |
| **F3.7** | Eöt-Wash + LLR Nordtvedt: automatic PASS z universal coupling | §3.3-§3.4 |
| **F3.8** | **R6 (QEP violations from Asorey-2015) closed strukturalnie** — TGP universal coupling immune | §4 |
| **F3.9** | R5 documented: cycle valid B ≪ B_QED; non-perturbative regime deferred | §5 |
| **F3.10** | TT10 magnetar test mechanism decoupling preserved: trace anomaly ≪ surface ρ_NS dla all magnetar regimes | §2.4 |

## §7 — Cross-cycle propagation list (Phase 4)

Phase 4 closure musi propagate następujące corrections / updates:

1. **L01 ADDENDUM §3.2 Q3:** typo correction `α²/(3π)≈7.7·10⁻⁷` → `α/(3π)≈7.74·10⁻⁴`;
   numerical estimates updated.
2. **L01 NEEDS.md §T.1 (N1):** numerical estimates updated (corrected prefactor).
3. **τ.3 ADDENDUM §2:** cross-cycle convergence diagnostic update (8-12 OOM
   separation depending on magnetar B; mechanism decoupling preserved).
4. **PREDICTIONS_REGISTRY.md:** new entry M911-EM-quantum z corrected estimates.
5. **L01 README.md:** Q1 closure update — *konstruktywnie* potwierdzony przez N1
   cycle dedicated derivation (Theorem 2.1).

## §8 — Phase 3 → Phase 4 handoff

### §8.1 — Co Phase 3 dało

1. **Lab + magnetar numerics** z corrected α/(3π) prefactor.
2. **MICROSCOPE/Eöt-Wash/LLR bounds:** automatic PASS strukturalnie.
3. **R5 documented** (B ≪ B_QED limitation).
4. **R6 closed strukturalnie** (TGP universal coupling immune do
   Asorey-2015 type QEP violations).
5. **L01 ADDENDUM typo identified** (factor ~1000) — propagation list.

### §8.2 — Co Phase 4 musi dostać

1. **Three-layer L1/L2/L3 closure** per
   [[../../meta/PPN_AS_PROJECTION.md]] §3.1 mandatory format.
2. **Native parameter audit** + forced-zero declarations.
3. **6/6 P-requirements verification** + final classification (STRUCTURAL_DERIVED
   target).
4. **Cross-cycle propagation** updates (§7).

## §9 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]] §3 NEEDS list
- [[./Phase1_results.md]] (β/(2α) prefactor LOCK)
- [[./Phase2_results.md]] (operator structure + GW170817)
- [[./Phase3_setup.md]]
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/ADDENDUM_2026-05-10_native_observables_first.md]] §3.2 Q3 (typo source — to be corrected)
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]] (η_TGP_Dirac baseline 1.32·10⁻²⁶)
- [[../op-tau3-substrate-clock-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §2
- JHEP 05 (2015) 118 (Asorey et al.) — Q6 risk-input

---

**Phase 3 close:** phenomenology consistent z bounds. R5 documented, R6 closed.
Phase 4 may proceed.
