---
title: "Phase 0 balance sheet retrofit — op-sigma1-substrate-light-dispersion (σ.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: op-sigma1-substrate-light-dispersion
cycle_path: "[[../op-sigma1-substrate-light-dispersion/README.md]]"
auditor: Claudian
classification: STRUCTURAL
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - retrospective
  - sigma1
  - phase3-medium-risk
  - positive-example
  - self-correction
related:
  - "[[../op-sigma1-substrate-light-dispersion/Phase3_results.md]]"
  - "[[retrofit_op-omega1_2026-05-06.md]]"
---

# Phase 0 balance sheet retrofit — σ.1 (op-sigma1-substrate-light-dispersion)

## Metadata cyklu

- **Cykl:** [[../op-sigma1-substrate-light-dispersion/README.md]]
- **Data oryginalnego closure:** 2026-04-30 (Phase 3 PASS, "FULL CONVERGENCE 4/4")
- **Self-correction:** 2026-05-01 (sign flip + scope softening + LIVE PARTIAL downgrade)
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian (M03 Phase 3, medium-risk #10 — chain ω.1→σ.1)
- **Klasyfikacja końcowa:** **STRUCTURAL** ★ (analog ω.1 pattern)

## 1. Co cykl twierdzi że robi

Z [[../op-sigma1-substrate-light-dispersion/Phase3_results.md]] (18/18
PASS, "FULL CONVERGENCE 4/4", post-2026-05-01 patch):

> "σ.1 c-mechanism dispersion ω_±² = k² ± g·k·n_∥ generuje 4 falsifiable
> channels. 3 alt-dispersions cross-channel FALSIFIED (scalar c(X),
> tensor Bumblebee, Lorentz-violating). σ.1 closes c-mechanism question:
> prędkość światła zależy od substrate **polarization-dependent**, NIE
> skalarnie."

**Self-correction 2026-05-01:**
- Dispersion form rewritten dimensionally-explicit
- Group velocity sign **corrected to POSITIVE** v_g,± = 1 + (g·n_∥)²/(8k²)
- "Effective optical metric" softened to "helicity-dependent optical cone"
- W3.3 CMB downgraded "PARTIAL POST-CONFIRM" → "LIVE PARTIAL candidate"

Główne claims:
- **C1**: Dispersion ω_±² = k² ± g·k·n_∥ helicity-dependent (NOT scalar c(X))
- **C2**: Phase velocity v_φ,± = 1 ± g·n_∥/(2k) linear birefringence
- **C3**: Group velocity v_g,± = 1 + (g·n_∥)²/(8k²) POSITIVE quadratic (sign-corrected 2026-05-01)
- **C4**: 3 alt-dispersions FALSIFIED (scalar c(X), tensor Bumblebee, Lorentz-violating)
- **C5**: CMB anisotropic birefringence beyond ω.1 isotropic — LIVE PARTIAL
- **C6**: Atomic clock orthogonal cross-check (no scalar α_em variation)

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs

```
- Webb/Murphy QSO Δα/α NULL 1e-7                   [Webb 2003-17, Murphy 2022]
- Planck PR4 + ACT 2024 β = 0.34±0.09° (3.8σ)      [Eskilt 2024 — shared z ω.1]
- Planck consistent kinematic dipole               [tensor Bumblebee falsifier]
- Fermi LAT GRB 090510 NULL energy-dep TOA          [Lorentz-violating falsifier]
- LIGO O5/O6 + cold-atom 2030+++                    [lab Mach-Zehnder]
- FAST/SKA-2 PTA 2030+ pulsar polarized timing      [pulsar dispersion]
- SO 2027+ + LiteBIRD 2029+ CMB 5σ POST-CONFIRM     [future]
```

### 2.2 Structural axioms (TGP-internal LOCKED)

```
- ω.1 axion-like (ln X)F·F̃ coupling                 [op-omega1, sympy LOCK]
- ω.1 EOM □(ln X) = -(g/f_X²) E·B                    [op-omega1]
- φ.1 substrate-action L=½(∂ ln X)²                  [op-phi1, AXIOM]
- WKB special case (n_i ∥ k_i)                       [stated approximation]
- Polarization-dependent (NOT scalar) c-mechanism    [σ.1 derivation]
```

### 2.3 Derived outputs

```
- O1: ω_±² = k² ± g·k·n_∥ (helicity-dependent dispersion)   (sympy Phase 1)
- O2: v_φ,± = 1 ± g·n_∥/(2k) linear birefringence            (Phase 2)
- O3: v_g,± = 1 + (g·n_∥)²/(8k²) POSITIVE quadratic         (sign-corrected)
- O4: NO scalar c(X) at leading order                        (consistent z Webb/Murphy)
- O5: CMB anisotropic + isotropic birefringence              (LIVE PARTIAL ~3.8σ)
- O6: Atomic clock NO scalar α_em variation                  (orthogonal cross-check)
- O7: 3 alt-dispersions FALSIFIED via 3 niezależnych obs     (multi-domain)
```

### 2.4 Tautology test (CRITICAL)

**O1-O3 (dispersion form):** Direct from ω.1 axion-like coupling z WKB
approximation. sympy LOCK Phase 2 7/7 PASS. NIE tautology — substantive
modified Maxwell w gradient. Sign correction 2026-05-01 v_g POSITIVE
(error fixed). **PASS**.

**O7 (3 alt-dispersions FALSIFIED):**
- scalar c(X) (dilaton): Webb/Murphy QSO NULL 1e-7 → uniform α_em(z)
  drift ✗
- tensor c(X) (Bumblebee): Planck CMB consistent z kinematic dipole
  → direction-dependent residual ✗
- Lorentz-violating (E_QG): Fermi LAT GRB 090510 NULL → energy-dep TOA ✗

**3 niezależne empirical FALSIFICATIONS** w różnych domains
(cosmological, CMB, GRB). σ.1 axion-induced UNIQUE survivor. **PASS
strong**.

**Werdykt tautology test:** PASS — dispersion substantive z sympy
LOCK; alt-falsification multi-domain.

### 2.5 Falsifiability test (CRITICAL)

**Multi-channel concrete falsifiers:**
- Lab Mach-Zehnder + parallel E·B: LIGO O5/O6 + cold-atom 2030+++
- PTA pulsar polarized timing: SKA-2 2030+ distinct from Faraday `(1+z)^{-2}`
- CMB anisotropic β(θ,φ): SO 2027+ + LiteBIRD 2029+ targets <0.05° for 5σ
- Atomic clock: NO scalar drift (orthogonal — falsifier: detected scalar drift > 1e-22/yr 2035+)

**LIVE PARTIAL CMB:** Planck PR4 + ACT 2024 β=0.34±0.09° at 3.8σ —
**partial detection NIE confirmed**, awaits SO/LiteBIRD 2027+
corroboration. Honest "LIVE PARTIAL" labeling (analog ω.1, consistent
self-correction 2026-05-01).

**Werdykt falsifiability test:** PASS — multi-channel forward + 1 LIVE
PARTIAL.

### 2.6 Independent-path cross-validation

**O1-O3 (dispersion form):** 1 main TGP path z ω.1 EOM modified
Maxwell. Single derivation, sympy LOCKED.

**Cross-channel falsification map (W3.5):** **3 alt-dispersions
FALSIFIED z 3 niezależnych empirical domains** (Webb/Murphy QSO,
Planck CMB dipole, Fermi LAT GRB). To jest **substantive multi-domain
falsification**, NIE pojedynczy obs.

**4-channel "FULL CONVERGENCE":** 1) plane-wave dispersion
POST-DERIVED, 2) v_φ/v_g sympy LOCK POST-DERIVED, 3) optical metric
POST-DERIVED, 4) CMB E/B chirality LIVE PARTIAL. Channels 1+2+3 są
**internal mathematical checks** dla wybranej dispersion form;
Channel 4 = single LIVE PARTIAL observational.

**Identical pattern do ω.1, φ.1:** "FULL CONVERGENCE 4/4" framing
borderline promotional dla 3 internal + 1 LIVE PARTIAL.

**Werdykt independent-path:** STRUCTURAL — 1 main TGP path +
substantive cross-channel alt-falsification (3 niezależnych obs);
NOT multi-domain niezależna DERIVATION.

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists (this file)
☑ Tautology test PASS (sympy LOCK + sign correction explicit)
☑ Falsifiability test PASS (4-channel + 3 alt-dispersions FALSIFIED multi-domain)
☐ Independent-path cross-validation PARTIAL — 1 main path; cross-channel falsification multi-obs
☑ Alt-scan ≥4 candidates (3 alt-dispersions explicit + 1 σ.1 winner = 4 forms tested)
☑ NIE used post-hoc structural motivations (modified Maxwell from ω.1 axion-like)
☑ NIE circular anchor (uses ω.1 + φ.1 axioms, no self-reference)
☑ NIE inheriting drift > parent × 5× (sign correction 2026-05-01 explicit, drift 0)
```

**7/8 ☑ + 1 ☐** — independent-path PARTIAL → max status STRUCTURAL.

**Self-correction 2026-05-01 substantive:**
- Group velocity sign corrected POSITIVE (was negative)
- Optical metric scope softened ("helicity-dependent cone", NOT full pseudo-Riemannian)
- W3.3 CMB downgraded LIVE PARTIAL (consistent ω.1)
- WKB approximation explicit (full CFJ-class beyond σ.1)

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? |
|-------|----------|
| DERIVED FULL | NO — single TGP dispersion path; multi-channel testing not multi-derivation |
| DERIVED CONDITIONAL | NO — wymaga ≥2 niezależne paths convergent; 1 main path z 4 alt-falsified |
| **STRUCTURAL** | **YES** — dispersion form sympy-LOCKED z ω.1 EOM; 3 alt-dispersions multi-domain FALSIFIED; LIVE PARTIAL CMB explicit |
| ANSATZ | NO — concrete falsifiers + multi-channel test |
| NUMEROLOGICAL | NO — alt-dispersions FALSIFIED przez real obs constraints (NIE constructed criterion) |
| TAUTOLOGY | NO — substantive dispersion w gradient |

**Final verdict:** **STRUCTURAL** ★

**Strukturalne cechy positive example:**

1. **Self-correction 2026-05-01 substantive:** sign flip v_g + scope
   softening + LIVE PARTIAL downgrade — analog ω.1, ψ.1 self-correction
   discipline pattern.
2. **3 alt-dispersions FALSIFIED z 3 niezależnych empirical domains**
   (Webb/Murphy + Planck + Fermi LAT) — substantive multi-domain
   falsification.
3. **WKB approximation explicit** (full CFJ-class beyond σ.1) —
   honest scope limitation.
4. **Sign correction explicit** (group velocity POSITIVE) —
   mathematical correctness over rhetorical convenience.
5. **Atomic clock orthogonal cross-check** — TGP nie ingeruje w
   scalar α_em (consistent z τ.2 scale-protection theorem).

**Phase 6 gate compliance:**
1. ✓ Phase0_balance.md exists (this file)
2. ✓ "PARTIAL POST-CONFIRM" → "LIVE PARTIAL" downgrade discipline
3. ✓ Brak constructed criterion (alt-dispersions FALSIFIED przez real obs)
4. ⚠ "FULL CONVERGENCE 4/4" framing borderline (analog ω.1, φ.1)
5. ✓ Brak sympy-rationalization-as-DERIVED (sympy LOCK = direct EL of ω.1)

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status YAML | `status: PASS`, "FULL CONVERGENCE 4/4", 18/18 PASS | STRUCTURAL — dispersion sympy-LOCKED + multi-domain falsification + LIVE PARTIAL explicit |
| Counter | "Cumulative ledger: 697 + 18 = 715" | Stays w sensie sub-tests substantive |
| Sub-tests | 5+7+6 = 18/18 PASS | Substantive (Phase 1 dispersion + Phase 2 velocity sympy LOCK + Phase 3 falsification map) |
| Independence | "FULL CONVERGENCE 4/4" | Re-characterized: 3 internal mathematical checks (POST-DERIVED) + 1 LIVE PARTIAL observational + 3 alt-dispersions empirical FALSIFIED |

## 6. Recommended action

- [x] **NO-OP klasy** — STRUCTURAL z honest reporting + self-correction
- [x] **Phase 5 PREDICTIONS_REGISTRY annotation:** "FULL CONVERGENCE
      4/4" → "1 dispersion form sympy-LOCKED + 3 alt-dispersions
      multi-domain FALSIFIED + 1 LIVE PARTIAL CMB"
- [ ] CRITIQUE — nie wymaga (cycle ALREADY self-corrected 2026-05-01)
- [ ] CASCADE_AUDIT — none (downstream τ.2 + τ.3 already audited)
- [ ] CORE_IMPACT — none

## 7. Notes

**σ.1 jako member chain φ.1→ω.1→σ.1→τ.2→τ.3→ψ.1:**

σ.1 jest **dispersion-channel** member 6-cycle TGP-substrate light/clock
chain. ω.1 dostarcza axion-like coupling, σ.1 dostarcza dispersion
relations, τ.2 scale-protection theorem, τ.3 sub-leading mass
acceleration, ψ.1 photon kinetic.

**Cross-cycle consistency:**
- σ.1 NO scalar c(X) ↔ τ.2 NO scalar α_em variation (consistent)
- σ.1 helicity-dependent ↔ ω.1 axion-like F·F̃ (parity-odd shared)
- σ.1 CMB LIVE PARTIAL ↔ ω.1 same Planck β=0.34±0.09° hint (shared signature)

**Self-correction discipline pattern (2026-05-01 cluster):**

| Cykl | 2026-05-01 correction |
|------|------------------------|
| ω.1 | "POST-CONFIRM" → "LIVE PARTIAL" + "B²" → "F·F̃ ∝ E·B" |
| **σ.1** | **v_g sign POSITIVE + WKB scope + LIVE PARTIAL** |
| ψ.1 | v1 INVALIDATED + v2 corrected + v3 Hilbert series |
| τ.3 | A5 audit patch (1/m_e factor inheritance flagged) |

4 cykle z **same-day self-correction discipline** 2026-05-01 — wzór
**spontaneous adoption** CALIBRATION_PROTOCOL §4 (pre-binding
2026-05-04).

## 8. Cross-references

- [[../op-sigma1-substrate-light-dispersion/README.md]] — cykl audited
- [[../op-sigma1-substrate-light-dispersion/Phase3_results.md]] — main claims
- [[../op-omega1-substrate-em-coupling/]] — ω.1 axion-like (parent)
- [[retrofit_op-omega1_2026-05-06.md]] — ω.1 retrofit (consistent)
- [[../op-phi1-substrate-action-variational/]] — φ.1 axiom (input)
- [[../op-tau2-substrate-time-coupling/]] — τ.2 (downstream pair)
- [[../op-psi1-substrate-light-acceleration/]] — ψ.1 (sub-leading L₅ pair)
- [[README.md]] — M03 master plan
- [[audit_log.md]] — appended 2026-05-06 (Phase 3 #10)
- [[tracker.md]] — status updated to DONE_STRUCTURAL
- [[../../meta/CALIBRATION_PROTOCOL.md]] — protocol source
- [[../../meta/CALIBRATION_GATE_ENFORCEMENT.md]] — Phase 6 gate
