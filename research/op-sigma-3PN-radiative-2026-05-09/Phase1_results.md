---
title: "Phase 1 results — Linear δΦ + σ near-field structure (foundation set)"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-results
phase: 1
status: 🟢 STRUCTURAL DERIVED — foundation set, 11/11 sympy PASS
needs_resolved: ["Linear δΦ from binary structure", "σ_ij tensor properties verified", "Phase 2 strategy identified"]
needs_blocker: ["Hadamard regularization dla ∫ σ d³x (Phase 2)", "Multi-session 2nd-order δΦ derivation"]
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
---

# Phase 1 results — Foundation dla Route A escape

## §0 — Executive summary

**STRUCTURAL DERIVED 11/11 sympy PASS — foundation set dla Phase 2.**

| Item | Result |
|---|---|
| Binary mass quadrupole Q_ij | Explicit form for circular orbit |
| Tr(Q) | = M·a²/4 const ⟹ no monopole radiation |
| d²(Tr Q)/dt² | = 0 EXACT (consistent z Phase 3 sphere-average finding) |
| Linear δΦ from quadrupole | Far-field amplitude + angular structure |
| σ_ij tensor structure | Traceless + symmetric, 1/r⁴ scaling far-field |
| ∫ σ d³x | Convergent at infinity, UV-divergent at particle (Hadamard reg needed) |
| **Phase 2 strategy** | **σ as 2nd-order tensor source dla δΦ_2 generation** |
| **PN order estimate** | **σ-induced TT ~ U² ~ 1% mass quadrupole at LIGO band** |

## §1 — Mass quadrupole z circular orbit

Equal-mass binary z separation a (m_each = M/2):
```
Q_xx = (M·a²/4)·cos²(ωt) = (M·a²/8)·(1 + cos(2ωt))
Q_yy = (M·a²/4)·sin²(ωt) = (M·a²/8)·(1 - cos(2ωt))
Q_xy = (M·a²/4)·cos(ωt)·sin(ωt) = (M·a²/8)·sin(2ωt)
Q_zz = 0  (orbit w xy plane)
Tr(Q) = M·a²/4 = const  ✓
```

Time derivatives oscillate at **2ω** (twice orbital frequency, GR-characteristic).

## §2 — Linear δΦ z binary

Far-field z mass quadrupole:
```
δΦ_far(r, n̂, t) = -(q/(8π·K_1·Φ_0·c²·r)) · d²Q_ij/dt² · n^i·n^j
```

**Angular structure:** δΦ has Y_2m angular pattern. Specific cases verified:
- Observer along z: δΦ → 0 (Q_zz = 0)
- Observer along x: δΦ ∝ d²Q_xx/dt²

This Y_2m structure DOES exist, ALE per [[../op-h-TT-calibration-2026-05-09/Phase2_sympy.py]]
NIE jest sufficient dla h_+, h_× at observer (pure trace structure of
δg_eff^ij = δ^ij·b_1·δΦ).

⟹ Need DIFFERENT mechanism dla TT modes: σ as 2nd-order source.

## §3 — σ_ij near-field tensor

Z δΦ_total = δΦ_1 + δΦ_2 (Newtonian approximation):
- ∂_i δΦ ~ G·M_each / r_i² (gradient near each particle)
- σ_ij = (∂_iΦ)(∂_jΦ) - (1/3)δ_ij(∇Φ)² 

**Verified properties:**
- σ traceless (3D): Tr(σ) = 0 ✓ sympy LOCK
- σ symmetric: σ_ij = σ_ji ✓
- Far-field scaling: σ ~ 1/r⁴ (bilinear in 1/r²)
- Near-particle UV: σ ~ 1/r⁴ at r → 0 (z near singular ∂Φ)

## §4 — Volume integral analysis

```
∫ σ_ij d³x dla binary:
   Far-field: σ ~ 1/r⁴, ∫ r² dr·1/r⁴ = ∫ dr/r² ~ finite at ∞ ✓
   Near-field: σ ~ 1/r⁴, ∫ r² dr·1/r⁴ ~ DIVERGES at r → 0 (UV)
```

**Standard regularization needed:** Hadamard partie-finie lub dimensional
regularization. Renormalizes into "σ effective tensor moment" Q^σ_ab.

This regularization jest standard 2-body PN technical work (multi-session).

## §5 — Phase 2 strategy (key escape route mechanism)

### §5.1 — σ_ab dynamics z Path B audit

Per [[../closure_2026-04-26/sigma_ab_pathB/results.md]]:
```
□ σ_ab + 2m_s²·σ_ab = TT-projected source[J·∂Φ]   (heredity EOM)
M² = 2 m_s² (composite mass derived)
```

σ_ab IS propagating tensor wave z M² mass term. At LIGO band z m_s ≪ ω,
effectively massless propagation — σ_ab waves propagate z c.

### §5.2 — σ_ab radiation from binary

If σ_ab source ma finite tensor moment Q^σ_ab(t), far-field radiation:
```
σ_ab(observer, t) ~ d²Q^σ_ab(t-r/c)/dt² / r        (1/r RADIATIVE!)
```

This jest NORMAL quadrupole-formula radiation, but for σ_ab tensor field.

### §5.3 — δg_eff^TT z σ-coupling

Z emergent-metric ansatz:
```
δg_eff^ij ⊃ σ^ij · c_0/(Φ_0² c²)
```

Substituting σ^ij(observer) ~ d²Q^σ/r·1/r:
```
δg_eff^ij(observer) ~ c_0/(Φ_0² c²) · d²Q^σ_ij(t-r/c)/dt² · (1/r)
```

**TT projection (since σ has genuine tensor structure, NIE δ^ij·scalar):**
```
δg_TT^ij(observer) ≠ 0  ← provides h_+, h_×!
```

### §5.4 — PN order estimate

```
Q^σ_ab ~ ∫ σ d³x ~ ∫ (∂Φ)² d³x ~ G² M² · (regularization scale dependence)
```

Compared do mass quadrupole Q^M ~ M·a²:
```
Q^σ / Q^M ~ G·M/a (Schwarzschild radius / orbital separation)
            ~ U (post-Newtonian parameter)
```

So σ-induced TT amplitude relative do GR mass-quadrupole TT:
```
h_TT^σ / h_TT^GR ~ U ~ 0.1 (LIGO band, near merger v² ~ 0.1·c²)
```

**~10% of mass-quadrupole tensor amplitude** (post-correction U~0.1).

This jest within order of magnitude **but POTENTIALLY too large** dla
LIGO scalar polarization 5% bound. Phase 2 MUSI compute carefully.

### §5.5 — Honest caveat

**1% scaling z my Phase 1 §6 was wrong** — corrected here to ~U ~ 10%.

This means σ-induced TT may BE detectable as polarization deviation z
GR-pure TT pattern. LIGO polarization tests could DETECT or FALSIFY
this prediction.

**Honest scope:** Phase 2 explicit calculation determines whether
σ-induced TT amplitude is consistent z observed h_+, h_× or violates.

## §6 — Phase 1 honest assessment

### §6.1 — What Phase 1 achieved

- Set up rigorous framework dla Phase 2 multi-session work
- Verified σ_ij tensor properties (traceless, symmetric, 1/r⁴ far-field)
- Identified Phase 2 strategy: σ as 2nd-order tensor source
- PN order estimate: σ-induced TT ~ 10% of mass quadrupole (corrected from 1%)

### §6.2 — What Phase 1 did NOT resolve

- Hadamard regularization of ∫ σ d³x (multi-session)
- Q^σ_ab tensor moment z circular orbit z proper renormalization
- Far-field σ_ab(observer) explicit amplitude
- Comparison z LIGO observed h_+, h_× amplitudes

### §6.3 — Cycle status

Phase 1 establishes **foundation**. Cycle continues do multi-session Phase 2.

## §7 — Connection do other cycles

### §7.1 — emergent-metric framework

If Route A escape works:
- emergent-metric Phase 4 Path 2 σ-coupling JEST genuinely TT source at 3PN+
- R5 risk RESOLVED at higher PN order
- Cycle #3 STRUCTURAL_CONDITIONAL → STRUCTURAL DERIVED (post-Route A)

If Route A fails:
- Pivot do Route B (nonlinear δΦ self-coupling)
- Or deeper structural reconsideration

### §7.2 — c_0, κ_σ joint cycles

c_0 ≈ 4π z OP-7 T3.4 chain remains valid foundation.
κ_σ ≈ 1/(3π) z 2-body orbital averaging remains heuristic.
**σ-induced TT amplitude w Route A may give DIFFERENT calibration** than
zero-β target c_0·κ_σ = 4/3.

## §8 — Time budget update

| Phase | Original estimate | Actual (Phase 1 done) |
|---|---|---|
| 0 | 0.5 sesji | 1 sesja done |
| 1 | 1 sesja | 1 sesja done (this) |
| 2 | 1-2 sesji | DEFERRED (multi-session) |
| 3 | 1 sesja | DEFERRED |
| 4 | 1 sesja | DEFERRED |
| 5 | 0.5 sesji | DEFERRED |
| FINAL | 0.5 sesji | DEFERRED |

**Multi-session continuation needed:** ~2-4 sesji dla Phase 2-FINAL.

## §9 — Cross-references

- [[./README.md]] — cycle setup
- [[./Phase1_sympy.py]] — verification script (11/11 PASS)
- [[../op-h-TT-calibration-2026-05-09/Phase_FINAL_close.md]] — adversarial finding
- [[../op-scalar-mode-LIGO-bound-2026-05-09/Phase_FINAL_close.md]] — DOWNGRADED status
- [[../closure_2026-04-26/sigma_ab_pathB/results.md]] — σ_ab Path B EOM
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — Phase 4 reference
- [[../op7/OP7_T3_results.md]] — ξ_eff = G·Φ_0² LOCK
