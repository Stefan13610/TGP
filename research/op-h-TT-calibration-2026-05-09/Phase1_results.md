---
title: "Phase 1 results — Calibration cycle reveals subtle Phase 3 cycle #3 issue"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-results
phase: 1
status: 🟠 STRUCTURAL_CONDITIONAL_HALT — multi-session needed, also affects cycle #3 verdict
needs_resolved: ["4√π factor source identified", "polarization decomposition ambiguity flagged"]
needs_blocker: ["Multi-session proper-TT-projection analysis"]
sympy_script: "[[./Phase1_sympy.py]]"
sympy_output: "[[./Phase1_sympy.txt]]"
---

# Phase 1 results — Calibration cycle: SUBTLE ISSUE IDENTIFIED

## §0 — Executive summary

**STRUCTURAL_CONDITIONAL_HALT 6/6 sympy PASS.** Calibration single-session
goal NOT achieved. ALE cycle ujawniła **subtelny błąd w Phase 3 cycle #3
verdict**, który może invalidate R5 mitigation claim.

| Item | Original Phase 3 cycle #3 | Phase 1 cycle calibration (corrected) |
|---|---|---|
| h_+, h_× from δ^ij·b_1·δΦ | "non-zero z l=2 multipole" | **ZERO at fixed observer** (isotropic spatial) |
| h_S at observer | "0 EXACT (sphere-averaged)" | **NON-ZERO at given observer angle** |
| Polarization pattern | TT-like | **Pure trace at observer** |
| LIGO compatibility | R5 mitigated | **R5 RISK RESTORED** |

## §1 — Phase 3 cycle #3 error detail

Phase 3 cycle #3 §5 claimed:
> "h_+ ~ d²(Q_xx - Q_yy)/dt²"

ALE w TGP linearized:
```
δg_eff^xx = b_1·δΦ(observer)
δg_eff^yy = b_1·δΦ(observer)
δg_eff^zz = b_1·δΦ(observer)
```

ALL THREE EQUAL — ISOTROPIC SPATIAL. Therefore:
```
h_+ = δg^xx - δg^yy = 0 IDENTICALLY
h_× = 2·δg^xy = 0 IDENTICALLY
```

**At ANY GIVEN OBSERVER POSITION**, δg_eff^ij from δ^ij·b_1·δΦ has NO
TT structure. δΦ may depend on observer angle (different value for different θ),
but at given θ, all three diagonal components are EQUAL.

⟹ **TGP linearized z δ^ij·b_1·δΦ ansatz alone gives NO h_+, h_×.**

## §2 — Sphere-averaged argument re-examined

Phase 3 cycle #3 §3 claimed:
```
h_S ~ ⟨δΦ⟩_sphere = (1/3)·Tr(d²Q/dt²) = 0 dla circular orbit
```

**Re-examination:** sphere-average IS zero, ALE this doesn't make h_S = 0
at given observer.

**Standard polarization decomposition in GW context:**
- h_S (scalar polarization at observer) = (1/3)·tr(δg^ij at observer)
- For TGP: tr(δg^ij at observer) = 3·b_1·δΦ(observer)
- ⟹ h_S(observer) = b_1·δΦ(observer) ≠ 0 generically

What sphere-average tells us:
- ∫|h_S(θ,φ)|² dΩ over sphere = total scalar polarization POWER radiated
- For binary, this MIGHT be reduced by trace cancellation (Q_xx + Q_yy = const)
- ALE at INDIVIDUAL observer, h_S ≠ 0

So Phase 3 sphere-average argument was about TOTAL POWER, NIE about
scalar mode amplitude AT GIVEN DETECTOR.

LIGO measures at specific detector position. If h_S at LIGO position is
non-zero, scalar polarization will be detected.

## §3 — σ-coupling jako TT source — corrected analysis

σ_ij = (∂_iΦ)(∂_jΦ) - (1/3)δ_ij(∇Φ)² IS angular tensor (NOT isotropic).

At observer position, σ_ij carries genuine tensor structure. Different
diagonal components σ_xx ≠ σ_yy generally.

ALE: σ_ij ~ (∂Φ)² ~ 1/r² (NEAR-FIELD), NIE 1/r (RADIATION).

⟹ At LIGO detector (large r_dist), σ contribution suppressed by
r_orbit/r_dist ~ 10⁻²⁰. **NEGLIGIBLE.**

## §4 — Cycle #3 verdict status update

### §4.1 — Phase 2 verdict RESTORED

```
Phase 2 cycle #3 (STRUCTURAL_NO_GO at linearized): RESTORED CORRECT
Phase 3 cycle #3 (R5 mitigated via multipole): IDENTIFIED INCORRECT
```

Phase 1 calibration cycle UJAWNIŁA subtle error w Phase 3 cycle #3 §5.

### §4.2 — R5 risk RESTORED

| Aspect | Phase 3 claim | Phase 1 calibration (corrected) |
|---|---|---|
| h_S at observer | 0 | b_1·δΦ(obs) ≠ 0 |
| h_+, h_× | from l=2 multipole | 0 z δ^ij·b_1·h alone; σ at 1/r² → negligible |
| LIGO bound | trivially satisfied | **VIOLATED** by ~factor 70× |

### §4.3 — Implications dla framework

**emergent-metric Phase 4 Path 2 status:**
- Cycle #3 Phase 3 verdict was "R5 mitigated" → status validated
- Phase 1 calibration verdict: Phase 3 was incorrect → status threatened

**Honest assessment:** R5 risk RESTORED. emergent-metric framework
linearized predicts SCALAR-DOMINANT GW polarization, sprzeczne z
observed TT-dominant pattern.

## §5 — Possible resolutions (multi-session work)

### §5.1 — σ-coupling at higher PN orders

W 3PN+ orbit, σ may couple via retarded Green function differently,
possibly providing TT contribution at 1/r far-field.

**Effort:** 3-5 sesji. Most likely resolution if framework valid.

### §5.2 — Nonlinear δΦ self-coupling

V_M9.1''(ψ) = -γψ²(4-3ψ)²/12 has nonlinear terms. Nonlinear δΦ products
may generate proper TT radiation at higher PN.

**Effort:** 3-5 sesji.

### §5.3 — Reformulation of g_eff functional dependence

g_eff^ij = δ^ij·B(ψ) + σ^ij·C(ψ) ansatz might be incomplete. Could
include velocity-dependent coupling V_ij(ψ, ∂ψ/∂t) z genuine tensor
structure at LINEAR (NIE bilinear) order.

**Effort:** 5-10 sesji (framework revision).

### §5.4 — Honest acknowledgment: emergent-metric needs multi-Φ extension?

Single-Φ S05 framework w linearized linearity gives only scalar
polarization at observer. To match GR's TT modes, need separate
tensor degrees of freedom — possible only z multi-Φ extension OR
S05 axiom modification.

**Effort:** 5-10 sesji + structural reconsideration.

## §6 — Honest CALIBRATION_PROTOCOL check

| Pattern | Status |
|---|---|
| Multi-candidate fit | NOT used |
| Constructed criterion post-hoc | NOT used |
| Drift hardening | NOT used |
| Algebraic re-arrangement | NOT used |
| Definitional tautology | NOT used |
| Sympy "DERIVED" without first-principles | **CRITICAL CAVEAT**: previous Phase 3 cycle #3 verdict had subtle error; this cycle CORRECTS by identifying the error |

**Honest reporting:** Phase 1 calibration cycle did its JOB by identifying
ambiguity in Phase 3 cycle #3 verdict. This jest GOOD adversarial analysis.

## §7 — Update needed dla previous cycles

### §7.1 — Cycle #3 status correction

Cycle #3 verdict should be revised:
- Phase 2 STRUCTURAL_NO_GO: CORRECT verdict (at linearized)
- Phase 3 R5 mitigation: INCORRECT — δ^ij·b_1·δΦ doesn't give h_+, h_×
- Cycle #3 final classification should be STRUCTURAL_CONDITIONAL (depends
  on multi-session resolution per §5)

### §7.2 — TGP_FOUNDATIONS §3.6.10.4 should be amended

§3.6.10.4 currently states "h_S = 0 strukturalnie via multipole".
This is INCORRECT for given-observer h_S.

Corrected: sphere-AVERAGE h_S = 0, but at GIVEN observer h_S ≠ 0.
Needs honest update.

### §7.3 — PREDICTIONS_REGISTRY similarly needs correction

Banner statement "h_S = 0 EXACT" should be "sphere-averaged h_S = 0,
local h_S ≠ 0 at observer".

## §8 — Verdict

**Phase 1 calibration cycle: 6/6 sympy PASS, STRUCTURAL_CONDITIONAL_HALT.**

Cycle did NOT achieve original goal (resolve 4√π factor). ALE achieved
SOMETHING MORE IMPORTANT: identified subtle error in cycle #3 Phase 3
that affects N14 R5 risk status.

**Adversarial analysis worked.** Honest reporting REQUIRED to update
upstream cycle results.

## §9 — Cross-references

- [[./README.md]] — cycle setup
- [[./Phase1_sympy.py]] — verification + error identification
- [[../op-scalar-mode-LIGO-bound-2026-05-09/Phase3_results.md]] — Phase 3 (now identified as incorrect)
- [[../op-scalar-mode-LIGO-bound-2026-05-09/Phase2_results.md]] — Phase 2 (now restored correct)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — N14 deferral
- [[../../TGP_FOUNDATIONS.md]] §3.6.10 — needs amendment
- [[../../PREDICTIONS_REGISTRY.md]] — needs amendment
