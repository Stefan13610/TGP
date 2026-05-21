---
title: "Phase 0 — Balance sheet RP² extension cyklu"
date: 2026-05-17
parent: "[[./README.md]]"
type: phase-zero-balance
phase: 0
status: 🟢 ACTIVE
sympy_substance_ratio: "6 FP / 1 LIT / 1 DEC = 75% FP"
hardcoded_T_pass: 0
---

# Phase 0 — Balance sheet RP² extension

## §1 — Inputs (LIVE inheritance)

### §1.1 — Z β-task (predecessor)

| Element | Source | Status |
|---|---|---|
| Source S = (2e/f_0)·(∂_μf_0)·A^μ | β-task T1 (2026-05-17) | LIVE LOCK |
| Linear-in-v scaling S ∝ v·B | β-task T3 | LIVE LOCK |
| Amplitude δθ_wake ~ e·B·v·L_kink | β-task T4 dimensional | LIVE |
| Gauge invariance U(1) | β-task T7 | LIVE LOCK |

### §1.2 — Z RP² geometry (PHASE3 why_n3)

| Element | Source | Status |
|---|---|---|
| Hedgehog orientation n(x)=x̂ | PHASE3 §2.1 (deg=1 na S²) | LIVE |
| RP² = S²/Z₂ antypodal | PHASE3 §2.2 (Q_eff=1/2) | LIVE |
| π₁(RP²) = Z₂ | PHASE3 §2.3 | LIVE |
| Berry phase γ=π pod 2π rotation | PHASE3 §2.4 | LIVE LOCK |
| Spin-1/2 emergence | PHASE3 §2.5 | LIVE |

### §1.3 — J_amp/J_phase split (lambda1)

| Element | Source | Status |
|---|---|---|
| Φ = |Φ|·exp(iθ) decomposition | lambda1 phase1L5 | LIVE |
| Independent amplitude/phase sectors | lambda1 | LIVE |

### §1.4 — Literature

| Element | Source | Use |
|---|---|---|
| Skyrme model hedgehog | Skyrme 1961 Proc Roy Soc A 260 | Analog reference (T6) |
| Berry phase formula | Berry 1984 Proc Roy Soc A 392 | γ_Berry computation (T6) |

## §2 — Outputs (intended deliverables)

**Phase 1 deliverables:**
- T1: f_0(r) magnitude profile spherical preservation pod RP² hedgehog
- T2: Source S structure identical z β-task (no modification w magnitude sector)
- T3: Linear-in-v preservation (β-task PASS robust)
- T4: n=0 winding consistency dla neutrino (∂θ_static=0 preserved)
- T5: Berry phase × motion → spinor-mediated coupling channel (NEW finding)
- T6: γ_Berry = π exact preservation
- T7: Gauge invariance pod A→A+∂λ
- T8: Structural equivalence theorem (β-task spherical ⇔ RP² hedgehog magnitude part)

**Phase FINAL deliverables:**
- R3 disposition (CLOSED structurally / refined / falsified)
- β-task verdict update (robust / refined / revised)
- Downstream impact: spinor channel candidacy dla μ_ν^TGP quantitative

## §3 — Risk register

### R1 — Hedgehog asymmetry (medium)

**Concern:** Hedgehog z orientation n(x)=x̂ może wprowadzać directional asymmetry
w f_0. W rzeczywistości jednak: spherical symmetry of hedgehog **preserves**
spherical f_0(r) bo orientation jest radial (każdy point on S² ma equivalent f_0).

**Mitigation:** T1 explicit sympy check że f_0 może być spherical w hedgehog
ansatz (decomposition Φ=f_0(r)·U(n(x))).

### R2 — Spinor-A coupling jest loop-level (medium)

**Concern:** Bezpośrednie coupling spinora do A_μ przebiega przez ψ̄γ^μA_μψ Dirac
term, NIE przez nasz scalar Lagrangian L = (∂|Φ|)² + |Φ|²(∂θ-eA)². Spinor jest
**emergent** z RP² topologii, więc jego coupling NIE jest fundamentalny.

**Mitigation:** T5 identifies channel heuristically (structural existence, NIE quantitative
loop computation). Quantitative deferred do dedicated cycle (post-W/Z sector).

### R3 — Skyrme-type non-minimal coupling (low)

**Concern:** Pełny TGP Lagrangian może zawierać (∇·n)², (n × ∂n)² terms
w analogii do Skyrme model. To może modyfikować source dla RP² geometry.

**Mitigation:** OUT OF SCOPE. Minimal U(1) coupling preserved per β-task framework;
Skyrme terms deferred to dedicated cycle.

## §4 — Sympy substance plan

| Test | Klasa | Substance | Hardcoded? |
|---|---|---|---|
| T1 | FP | Decomposition algebra symbolic | 0 |
| T2 | FP | Symbolic source comparison | 0 |
| T3 | FP | Sympy series + limit (analog β-task T3) | 0 |
| T4 | FP | Symbolic n=0 winding check | 0 |
| T5 | FP | Heuristic Berry-motion cross-term | 0 |
| T6 | LIT | γ_Berry exact computation | 0 |
| T7 | DEC | Gauge transformation algebra | 0 |
| T8 | FP | Structural theorem comparison | 0 |

**FP: 6/8 = 75% ✓.** Hardcoded T_pass=True: **0**.

## §5 — 8/8 gate

- [x] G1 — Contract WRITTEN BEFORE Phase 1 ✓
- [x] G2 — Pre-registered falsification (README §0.2) ✓
- [x] G3 — Native-first methodology (README §0.3) ✓
- [x] G4 — Pre-flight read confirmation (README §0.4) ✓
- [x] G5 — Sympy substance plan ≥75% FP, 0 hardcoded ✓
- [x] G6 — Risk register (3 risks) ✓
- [x] G7 — Decision tree explicit (robust/refined/revised) ✓
- [x] G8 — Downstream impact identified ✓

**Verdict:** 🟢 **8/8 PASS. Phase 1 sympy authorized.**

## §6 — Cross-references

- [[./README.md]] — scope + contract
- [[../op-neutrino-omega-motion-wake-2026-05-17/Phase_FINAL_close.md]] — β-task predecessor
- [[../why_n3/PHASE3_RP2_defect_quantization.md]] — RP² geometry source
- [[../../meta/CALIBRATION_PROTOCOL.md]] — 8/8 gate

---

**Sign-off:** Claudian @ 2026-05-17 (2nd cycle sesji). **Gate 8/8 PASS.**
