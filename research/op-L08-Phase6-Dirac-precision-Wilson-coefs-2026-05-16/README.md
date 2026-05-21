---
title: "op-L08-Phase6-Dirac-precision-Wilson-coefs — Wilson coefs Φ-dependent corrections + a_e lab-scale estimate"
date: 2026-05-16
pre_registration_date: 2026-05-16
parent: "[[../op-L08-Phase6-Dirac-propagator-2026-05-16/Phase_FINAL_close.md]]"
cycle: L08 Phase 6 precision extension (Wilson coefs)
status: CLOSED-RESOLVED — A− single-session sprint 2026-05-16
folder_status: closed-resolved
claim_status: A-
closure_date: 2026-05-16
may_edit_core: false
audit_source: "[[../op-L08-Phase6-Dirac-propagator-2026-05-16/Phase_FINAL_close.md]] §10.1 (Phase 2 extension candidates)"
priority: P2 (extension precision)
related:
  - "[[../op-L08-Phase6-Dirac-propagator-2026-05-16/]] (predecessor cycle A−; S_F^TGP construction)"
  - "[[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] (Wilson coef framework reference; α/(3π) prefactor LIVE)"
  - "[[../op-newton-momentum/B9_wep_microscope_composition_results.md]] (B9 MICROSCOPE η_TGP_lab = 1.32·10⁻²⁶ baseline LIVE)"
  - "[[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] (m_obs LIVE)"
  - "[[../../audyt/L08_kink_fermion_closure/README.md]] (cluster context)"
  - "[[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] BINDING contract spec"
tgp_status:
  level: L1
  kind: derivation
  output_type: observable  # a_e prediction → PDG comparison
  core_compatibility: extension
  may_edit_core: false
  exports_findings: true
  depends_on:
    - "L08-Dirac S_F^TGP construction (A− 2026-05-16)"
    - "L01-N1 quantum trace anomaly prefactor α/(3π) (A− 2026-05-11)"
    - "B9 MICROSCOPE WEP universal coupling (8.3·10¹⁰× safe)"
    - "L05 m_obs LIVE"
    - "S05 single-Φ axiom"
  impacts:
    - "L08 audit warstwa 3c partial-(D) STRENGTHENED z lab-scale falsifiability evidence"
    - "TGP fenomenologia lab-scale (QED a_e match)"
    - "Future cycles: quark/neutrino/boson sector precision"
tags:
  - TGP
  - L08
  - Phase6
  - precision
  - Wilson-coefs
  - emergent-Dirac
  - anomalous-magnetic-moment
  - a_e
  - lab-scale
  - PDG-comparison
  - pre-registration-locked
---

# op-L08-Phase6-Dirac-precision-Wilson-coefs — BINDING contract

> **Cykl scaffold:** 2026-05-16 (post-L08-Dirac-propagator A− closure)
> **Status:** PARKING (Phase 0 scaffold ready)
> **claim_status:** TBD pending Phase 1 execution
> **Pre-registration date:** 2026-05-16 (immutable, BINDING)

## §0 — BINDING contract (pre-Phase-1 LOCKED)

### §0.1 — Cycle goal (binding)

**Wyprowadzić Wilson coef framework** dla S_F^TGP corrections w obecności background
Φ ≠ Φ_0 + **estymować leading-order a_e (electron anomalous magnetic moment) prediction**
w TGP-substrate.

**Inheritance starting point:** L08-Dirac-propagator A− closure provides:
```
S_F^TGP(p)|_{Φ_0} = i(γ·p + m_obs)/(p² − m_obs² + iε)
                 = S_F^Dirac canonical (Φ_0 limit verified Phase 1 T10)
```

**This cycle extends to Φ ≠ Φ_0:**
```
S_F^TGP(p; Φ) = S_F^TGP|_{Φ_0} + ζ_1 · (δΦ/Φ_0) · S_F^(1) + ζ_2 · (δΦ/Φ_0)² · S_F^(2) + O((δΦ)³)
```

z Wilson coefs ζ_1, ζ_2 constrained by:
- Universal coupling g_eff[Φ] (S04 + B9 WEP)
- L01 ρ ≡ -T^μ_μ/c_0² formal definition
- M9.1'' geometric structure (Lorentz signature preservation)

**Lab-scale prediction target:** a_e^TGP correction magnitude, comparison z PDG bound.

**3-warstwowa specyfikacja:**

- **L1 (native predictions):** Wilson coefs ζ_i derived z g_eff[Φ] structure; lab-scale
  (δΦ/Φ_0)_lab estimate z B9 MICROSCOPE inheritance; a_e^TGP correction
- **L2 (standard QED projection):** a_e^QED = α/(2π) + O(α²) + ... (LIT-anchored); TGP
  correction must be ≪ QED at lab scale
- **L3 (falsifikator):** PDG 2024 a_e = 0.00115965218073(28) z precision 0.28 ppb

### §0.2 — Pre-registered falsification rule (LOCKED, immutable)

> *Jeśli Phase 1 estymata a_e^TGP correction > 1 ppb (above experimental precision threshold),
> TGP emergent Dirac propagator z S05 universal coupling jest niezgodny z PDG a_e
> measurement → strukturalna konieczność reframing universal coupling lub Wilson coef
> structure. Komplementarnie: jeśli estymata < 0.001 ppb, TGP fenomenologicznie
> indistinguishable z QED at lab scale (success scenario: no detectable lab-scale
> deviation, structurally guaranteed S05 + B9 safety).*

**Recovery scope (anti-Lakatos LOCKED):**

```
allowed:
  - Wilson coef precision z multi-loop QED extension (post-Phase 1)
  - High-curvature regime (astrophysical) tests post-lab-scale verification
  - Cross-check z N1 EM trace anomaly α/(3π) prefactor LIVE
forbidden:
  - Post-hoc tuning ζ_1, ζ_2 dla PDG match
  - Universal coupling violation (S04 + B9 LOCKED preserve)
  - S05 single-Φ violation
if_recovery_exhausted:
  H1c: TGP Wilson coef framework requires non-universal coupling
       → extension narusza B9 MICROSCOPE; revisit cycle z explicit re-derivation
```

### §0.3 — Six P-requirements (declarative)

| # | P-requirement | Resolution mechanism |
|---|---|---|
| **P1** | Wilson coef expansion S_F^TGP[Φ] explicit | Phase 1 T1-T3: expansion ze background g_eff[Φ] |
| **P2** | Wilson coefs ζ_i constrained by g_eff structure | Phase 1 T4-T5: ζ_i derived (not free; analog L01-N1 prefactor) |
| **P3** | Lab-scale (δΦ/Φ_0)_lab estimate | Phase 1 T6: z B9 MICROSCOPE η_lab = 1.32·10⁻²⁶ inheritance |
| **P4** | a_e^TGP correction magnitude | Phase 1 T7-T8: combination ζ_i × (δΦ/Φ_0)_lab |
| **P5** | PDG comparison + falsification check | Phase 1 T9-T10: a_e^TGP_corr vs 0.28 ppb tolerance |
| **P6** | S05 + B9 + universal coupling preservation | Phase 1 T11 DEC: declarative inheritance check |

### §0.4 — Inheritance ledger

| Source | Element | Status |
|---|---|---|
| L08-Dirac-propagator | S_F^TGP|_{Φ_0} canonical form | ✅ LIVE A− |
| L08-Dirac-propagator | γ^μ Cl(1,3) operators | ✅ LIVE A− |
| L08-Dirac-propagator | m_obs pole z L05 inheritance | ✅ LIVE A− |
| L01-N1 EM trace anomaly | β(α)/(2α) = α/(3π) ≈ 7.74·10⁻⁴ prefactor structure | ✅ LIVE A− |
| B9 MICROSCOPE | η_TGP_lab = 1.32·10⁻²⁶ (8.3·10¹⁰× safe) | ✅ LIVE LOCKED |
| L05 | m_obs LIVE; k_obs(α=2,d=3)=3 | ✅ LIVE A− |
| Q2 F1 | Φ_eq(today) = H_0 anchor | ✅ LIVE LOCKED |
| L01 | ρ ≡ -T^μ_μ/c_0² | ✅ LOCKED |
| S05 single-Φ axiom | no new fundamental fields | ✅ LOCKED |
| PDG 2024 a_e | 0.00115965218073(28) | external LIT |

### §0.5 — Output type (binding)

`output_type: observable` — a_e correction magnitude IS native observable (compared
z PDG 2024 a_e measurement at 0.28 ppb precision).

**Initial target:** **A−** (STRUCTURAL_DERIVED_NATIVE_PARTIAL); upgrade do **A**
possible if Phase 2 includes multi-loop QED comparison.

### §0.6 — Phase plan

| Phase | Goal | Effort | Tests |
|---|---|---|---|
| **Phase 0** | Balance sheet 8/8 ☑ | ~15 min | N/A |
| **Phase 1** | Wilson coef expansion + a_e estimate | ~1-2 sesje | 11 tests target |
| **Phase FINAL** | Closure ceremony | ~15 min | N/A |

**Single-session ambition:** Phase 0 + Phase 1 + Phase FINAL → A−.

## §1 — Activation pre-conditions

1. ✅ Phase 0 balance sheet (per Phase0_balance.md)
2. ⏸ Authorization "Wilson coefs Phase 1 start"
3. ⏸ WIP slot

## Status

🟢 **CLOSED-RESOLVED** — single-session sprint 2026-05-16 (user "działaj dalej z cyklem").

**Phase 1 verdict:** **11/11 sympy PASS** — 8 FP (72.7%) + 2 LIT + 1 DEC; 0 hardcoded.
**6/6 P-requirements RESOLVED.**

**claim_status A−** (STRUCTURAL_DERIVED_NATIVE_PARTIAL).

**Kluczowy wynik fizyczny:**
- |δa_e^TGP| ~ 5.66·10⁻²⁹ (Wilson coef estimate)
- PDG 2024 a_e precision: 2.8·10⁻¹³ (0.28 ppb)
- Ratio TGP/PDG ~ 2.02·10⁻¹⁶ → **16 OOM safety margin**

**Strukturalna gwarancja PDG compliance przez S05 + B9 universal coupling.** TGP
fenomenologicznie indistinguishable z QED at lab scale — NIE przez fine-tuning,
ale przez fundamental inheritance.

**Full closure documentation:** [[./Phase_FINAL_close.md]].

## Cross-references

- [[../op-L08-Phase6-Dirac-propagator-2026-05-16/]] (predecessor A−)
- [[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] (Wilson framework reference)
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]] (B9 baseline)
- [[../op-L05-mass-exponent-k-alpha-d-2026-05-16/]] (m_obs LIVE)
- [[../../audyt/L08_kink_fermion_closure/README.md]] (cluster context)
- [[../../meta/CYCLE_LIFECYCLE.md]] (claim_status taxonomy)
- [[../../meta/CALIBRATION_PROTOCOL.md]] (Phase 6 ABSOLUTE BINDING)
