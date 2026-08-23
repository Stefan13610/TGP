---
title: "Phase 2 results — op-CE-H-3D-native-interaction-2026-05-22 — Two-vortex V_int(L) interaction"
type: phase_results
status: COMPLETE_HONEST_PARTIAL
phase: 2
parent_cycle: op-CE-H-3D-native-interaction-2026-05-22
date_completed: 2026-05-23
sympy_total: "4/4 PASS execution (3 substantive PASS + 1 HONEST FAIL)"
substance_metrics: "3/4 FP substantive PASS (75%); 0 hardcoded T_pass=True (strict cycle 1/2/7 preserved); 0/1 DEC budget used (preserved); 1 honest FAIL (sign convention pre-registration error)"
verdict: PROCEED_TO_PHASE_3_WITH_HONEST_CAVEAT
honest_fail_pattern: "Same class as R1-1 (CE-H T_P3_2 + FFS-3 q=1 implicit). NEW META-IRONY: violation of §3.6 BINDING within first cycle post-§3.6 enforcement."
---

# Phase 2 results — Two-vortex V_int(L) interaction (HONEST PARTIAL)

**Status:** COMPLETE_HONEST_PARTIAL 2026-05-23
**Sympy:** [[./Phase2_sympy.py]] + [[./Phase2_sympy.txt]]

---

## §0 — Verdict + summary

```
████████████████████████████████████████████████████████████████████
█  op-CE-H-3D-native-interaction-2026-05-22 PHASE 2                █
█                                                                  █
█  PHASE 2 SYMPY: 4/4 PASS execution; 3/4 substantive PASS         █
█  SUBSTANCE METRIC: 3/4 FP substantive (75%); 0 hardcoded         █
█  STRICT CYCLE 1/2/7 PATTERN PRESERVED                            █
█                                                                  █
█  KEY RESULTS:                                                    █
█    Single-vortex |∇θ|² = n²/r_⊥² (sympy exact) ✓                █
█    2D Goldstone Green function G_2D = -log/(2π) ✓                █
█    Two-vortex log behavior CONFIRMED (99% magnitude match):     █
█      Numerical |slope| ≈ 6.22 vs analytical 2π ≈ 6.28           █
█    Dimensional consistency ✓                                     █
█                                                                  █
█  HONEST FAIL T_P2_4:                                             █
█    Pre-registered: slope = +2π (POSITIVE)                       █
█    Numerical: slope = -6.22 (NEGATIVE)                          █
█    Magnitude perfect match; SIGN convention pre-reg ERROR        █
█    Same class as T_P3_2 + FFS-3 + R1-1 pattern                  █
█                                                                  █
█  META-IRONY: Phase 0 §8 violated §3.6 BINDING (heuristic sign)  █
█    despite §3.6 being first-cycle-post-BINDING enforcement      █
█                                                                  █
█  F-γ-1 SUBSTANCE verdict (preliminary):                          █
█    Log behavior CONFIRMED, CLEARLY distinguishable from exp ✓    █
█    Phase 3 sign-agnostic differential test will finalize         █
████████████████████████████████████████████████████████████████████
```

---

## §1 — Pre-registered tests verdict table

| Test | Type | Pre-registered claim | Status | Substantive |
|------|------|---------------------|--------|-------------|
| T_P2_1 | LIT | (Phase 1 LIT covers cycle) | N/A | — |
| T_P2_2 | FP | Single-vortex \|∇θ\|² = n²/r_⊥² | PASS | ✓ |
| T_P2_3 | FP | G_2D = -log(r/r_0)/(2π) | PASS | ✓ |
| T_P2_4 | FP | Slope = +2π for V_int log scaling | **HONEST FAIL** | sign pre-reg error |
| T_P2_5 | FP | Dimensional consistency | PASS | ✓ |

**Substantive FP: 3/4 PASS (75%).**

---

## §2 — Per-test detailed findings

### §2.1 T_P2_2 [FP] Single-vortex |∇θ|²

**Computation:**
$$\theta(x,y) = n \cdot \arctan(y/x)$$
$$|\nabla\theta|^2 = \left(\frac{\partial\theta}{\partial x}\right)^2 + \left(\frac{\partial\theta}{\partial y}\right)^2 = \frac{n^2}{x^2+y^2} = \frac{n^2}{r_\perp^2}$$

**Sympy verification:** $|\nabla\theta|^2 - n^2/r_\perp^2 = 0$ exact.

**Verdict:** ✓ PASS

### §2.2 T_P2_3 [FP] 2D Goldstone Green function

**Computation:** G_2D(r_⊥) = -log(r_⊥/r_0)/(2π).

**Sympy verification:** $\nabla^2 G_{2D} = 0$ away from origin (Laplace eqn satisfied).

**Physical interpretation:** Mediates 2D Coulomb interaction (logarithmic, NIE exponential). Reflects 3D massless Goldstone in z-translation-invariant geometry.

**Verdict:** ✓ PASS

### §2.3 T_P2_4 [FP] Two-vortex V_int(L) — HONEST FAIL (sign convention)

**Pre-registered claim (Phase 0 §8.6):**
$$V_{int}(L)/L_z \approx +2\pi v^2 n_1 n_2 \log(L/r_0)$$

with **POSITIVE** slope coefficient $+2\pi n_1 n_2$ (for n_1=n_2=1: +2π ≈ +6.2832).

**Numerical computation (sympy Phase 2):**

Two vortices at $\vec{x}_{1,\perp} = (-L/2, 0)$ and $\vec{x}_{2,\perp} = (+L/2, 0)$, windings $n_1 = n_2 = 1$:

| L | V_int(L) numeric |
|---|------------------|
| 1.0 | 24.6676 |
| 2.0 | 20.4723 |
| 4.0 | 15.8797 |
| 8.0 | 11.5382 |
| 16.0 | 7.2423 |
| 32.0 | 3.1137 |

**Differential analysis** $\Delta V/\Delta \log(L)$:

| L_1 - L_2 | slope |
|-----------|-------|
| 1-2 | -6.0525 |
| 2-4 | -6.6258 |
| 4-8 | -6.2635 |
| 8-16 | -6.1977 |
| 16-32 | -5.9562 |

**Average slope:** -6.2191

**Magnitude analysis:**
- |numerical| / |analytical| = 6.2191 / 6.2832 = **0.9898 (99% match)**
- Magnitude pre-derivation correct to 1%

**Sign analysis:**
- Pre-registered: +2π
- Numerical: -2π
- **Sign opposite** of pre-registered

**Anti-Lakatos discipline LOCKED:**
- ✗ NIE modified threshold ex post (forbidden #1)
- ✗ NIE re-defined slope sign convention to fit pre-reg
- ✗ NIE hidden the FAIL
- ✓ Reported T_P2_4 as **HONEST FAIL** per literal threshold
- ✓ Documented sign convention error w pre-derivation

**Verdict:** ❌ **HONEST FAIL** per literal pre-registered threshold.

### §2.4 T_P2_4 honest declaration — root cause

**Source of pre-registration sign error (Phase 0 §8.6):**

Phase 0 §8.6 wrote:
> "For 2D Laplacian: $\nabla^2 \phi_i = 2\pi n_i \delta^2(\vec{r}_\perp - \vec{x}_{i,\perp})$, and 2D Green function $G_{2D}(\vec{r}_\perp) = -\frac{1}{2\pi}\log(|\vec{r}_\perp|)$.
>
> $\frac{E_{int}}{L_z}(L) \approx 2\pi v^2 n_1 n_2 \cdot \log(L/r_0)$"

**Error:** The coupling of vortex sources via 2D Coulomb Green function gives:

$$V_{int}(L)/L_z = -(q_1 \cdot q_2)/(2\pi) \cdot \log(L/r_0) \cdot \text{factor}$$

with effective "charges" $q_i \propto 2\pi v n_i$. The **negative sign** comes from standard convention dla scalar Goldstone field (analog 2D electrostatic potential of "charges").

**Physical check (correct sign):**

Same-sign vortices ($n_1 \cdot n_2 > 0$) should **REPEL** (like-charges in 2D Coulomb). This means:
- V_int(L) HIGH at small L (close vortices)
- V_int(L) LOW at large L (far vortices)
- V_int **DECREASES** with L → **NEGATIVE slope**

Numerical confirms: V_int(1) = 24.7 (close) > V_int(32) = 3.1 (far). REPULSION ✓.

**Conclusion:** Numerical result physically CORRECT (like-charges repel). Pre-registered sign was based on **heuristic intuition** about "interaction increasing with separation" — same class methodology error as R1-1.

### §2.5 T_P2_5 [FP] Dimensional consistency

**Dimensional analysis (natural units ħ=c=1):**

| Quantity | Dimension |
|----------|-----------|
| v (field VEV) | mass = 1/length |
| n (winding) | dimensionless |
| L, r_0 | length |
| log(L/r_0) | dimensionless |
| **v²·n²·log** | **mass²** |
| **E/L_z = energy/length** | **mass²** |

Match: ✓ true.

**Verdict:** ✓ PASS

---

## §3 — META-IRONY: §3.6 BINDING violation in first cycle post-BINDING

**Critical observation:**

CALIBRATION_PROTOCOL §3.6 BINDING (2026-05-22) was formalized SPECIFICALLY do prevent pre-registration analytical pre-derivation gaps (R1-1 pattern). This cycle (Poziom γ-1) is the FIRST cycle post-§3.6 BINDING.

Phase 0 §8 explicitly claimed §3.6 BINDING compliance via analytical pre-derivation section. **Yet sign convention error slipped through.**

**Why this happened:**
- §8.4 derived Green functions correctly (analytical)
- §8.5-§8.6 jumped to interaction formula z "≈" notation
- Sign convention was IMPLICIT (heuristic intuition "log grows positive")
- Anti-Lakatos check w §8.9 caught coupling form linearity + coupling constant C as heuristic, BUT MISSED sign convention as heuristic

**Methodology lesson dla future:**
- §3.6.2 forbids heuristic shortcuts dla DERIVED VALUES (m vs m·√2)
- **EXTENSION needed:** §3.6.2 also forbids heuristic shortcuts dla SIGN CONVENTIONS
- Phase 0 §8 analytical pre-derivation MUST include explicit sign verification (positive/negative pre-derivation derived from physics, NIE assumed)

**R1 flag created (NEW R1 dla future R2 audit):**

R1-2: **§3.6 BINDING extension required** dla sign conventions
- Pattern: T_P2_4 sign error (this cycle, 2026-05-23) + R1-1 m vs m·√2 (CE-H 2026-05-21) + FFS-3 q=1 implicit (FFS 2026-05-20)
- THREE instances of analytical pre-derivation gap (different aspects)
- Mitigation candidate: §3.6.2 explicit "sign conventions MUST be derived, NIE assumed" addition
- Authority: Future R2 audit cycle decides whether to extend §3.6

---

## §4 — F-γ-1 substantive status (preliminary)

**Question dla F-γ-1 CRUCIAL TEST:** Czy native 3D interaction tail jest power-law/log lub pure exponential?

**Phase 2 substantive evidence:**

1. **Magnitude:** Numerical |slope| matches analytical 2π at 99% accuracy.
2. **Scaling:** V_int(L) varies as **log(L)** (not exp(-mL)) — verified across factor 32 range.
3. **Physics:** Like-sign vortices REPEL via logarithmic interaction (correct cosmic string global-vortex behavior).
4. **Discrimination:** At L=32, V_int/V_int(L=2) = 3.1/20.5 = 0.152. For pure exponential exp(-m·30) z m=1: ratio = 10⁻¹³. **Clearly distinguishable.**

**Substantive F-γ-1 status: PRELIMINARILY PASS** (sign-agnostic test in Phase 3 will finalize).

**Honest caveat:** T_P2_4 sign FAIL is methodology issue, NIE substance issue. **Substance (log scaling) confirmed**.

---

## §5 — Anti-Lakatos discipline check

- ✅ T_P2_4 HONEST FAIL reported per literal threshold (NIE modified)
- ✅ Substance verification 99% accuracy noted, NIE adopted as PASS rescue
- ✅ Phase 0 §8 sign convention error documented explicit (NIE hidden)
- ✅ Pattern detection (3rd instance) honest declared
- ✅ Methodology lesson formalized (R1 flag dla R2 extension)
- ✅ NO Phase 0 ex post modifications
- ✅ Sympy results IMMUTABLE per audit trail

---

## §6 — Implications dla Phase 3

**Phase 3 sympy goal:** Differential test V_int(L) ~ log(L) vs exp(-mL).

**Phase 3 approach (sign-agnostic):**
- Fit |V_int(L)| vs log(L) and exp(-m_σ L)
- Compute R²_log vs R²_exp
- F-γ-1 PASS criterion: R²_log > 0.95 AND R²_log > R²_exp + 0.02 (per Phase 0 §8.10)

**Expected outcome:** R²_log ≈ 0.99+ (excellent log fit per Phase 2 data); R²_exp << R²_log (poor exponential fit because data spans 30+ in L with magnitude ratio ~8x — not exponential).

**Phase 3 will be sign-agnostic** to bypass T_P2_4 sign convention issue. F-γ-1 verdict will be based on FORM of decay (log vs exp), NIE sign.

---

**END OF Phase 2 results — PROCEED_TO_PHASE_3_WITH_HONEST_CAVEAT LOCKED 2026-05-23**

**Aggregate Phase 2:** 3/4 substantive PASS + 1 HONEST FAIL (sign convention pre-registration error). Substance F-γ-1 PRELIMINARILY PASS pending Phase 3 sign-agnostic differential test. Pattern detection: 3rd instance analytical pre-derivation gap → R1 flag dla §3.6 extension.
