---
title: "Phase 1 results — op-CE-H-3D-native-interaction-2026-05-22 — Vortex ansatz EL + mass spectrum"
type: phase_results
status: COMPLETE
phase: 1
parent_cycle: op-CE-H-3D-native-interaction-2026-05-22
date_completed: 2026-05-23
sympy_total: "5/5 PASS execution (4 FP + 1 LIT)"
substance_metrics: "4/4 FP substantive PASS (100%); 0 hardcoded T_pass=True (strict cycle 1/2/7 preserved); 0/1 DEC budget used (vortex default, NOT DEC switch)"
verdict: PROCEED_TO_PHASE_2
---

# Phase 1 results — Vortex ansatz EL + mass spectrum + RP²/Z₂ compatibility

**Status:** COMPLETE 2026-05-23 (Phase 1 closed PROCEED_TO_PHASE_2)
**Sympy:** [[./Phase1_sympy.py]] + [[./Phase1_sympy.txt]]

---

## §0 — Verdict + summary

```
████████████████████████████████████████████████████████████████████
█  op-CE-H-3D-native-interaction-2026-05-22 PHASE 1                █
█                                                                  █
█  PHASE 1 SYMPY: 5/5 PASS (4 FP + 1 LIT)                         █
█  SUBSTANCE METRIC: 4/4 FP (100%); 0 hardcoded; 0/1 DEC budget    █
█  STRICT CYCLE 1/2/7 PATTERN PRESERVED                            █
█                                                                  █
█  KEY RESULTS:                                                    █
█    Vortex EL equation derived from TGP Lagrangian                █
█    Mass spectrum analytically verified:                          █
█      m_σ = v·√(2λ) (Higgs massive)                              █
█      m_π = 0 (Goldstone massless, dla global U(1))              █
█    Far-field structure: power-law -n²/(2λv·r²) + exp(-m_σr)/√r  █
█    Z₂ × RP² compatibility preserved                              █
█                                                                  █
█  AGGREGATE VERDICT: PROCEED_TO_PHASE_2                           █
████████████████████████████████████████████████████████████████████
```

---

## §1 — Pre-registered tests verdict table

| Test | Type | Pre-registered claim | Status | Substantive |
|------|------|---------------------|--------|-------------|
| T_P1_1 | LIT | 5 literature anchors | LIT_PASS | Informational |
| T_P1_2 | FP | Vortex EL derived from TGP L | PASS | ✓ |
| T_P1_3 | FP | Far-field ρ→v z m_σ rate | PASS | ✓ |
| T_P1_4 | FP | Mass spectrum (Higgs massive, Goldstone massless) | PASS | ✓ |
| T_P1_5 | FP | Z₂ × RP² compatibility | PASS | ✓ |

**Substantive FP: 4/4 PASS (100%).**

---

## §2 — Per-test detailed findings

### §2.1 T_P1_1 [LIT] Literature anchors

Literature checkpoint: 5/5 anchors verified (Nielsen-Olesen, Vilenkin-Shellard, Manton-Sutcliffe, Goldstone theorem references, Peskin-Schroeder).

**Critical methodological note:** 4 of 5 anchors use Abelian Higgs framework (LOCAL gauge). TGP S05 in basic axiomatization has GLOBAL U(1) (gauge emergence post Foundations §3.4). **Different mass spectrum:** Goldstone remains massless (NOT eaten by gauge field).

This is the **key analytical reason** that 3D long-range structure exists in TGP (per Phase 0 §8 pre-derivation), unlike pure Abelian Higgs which has exponential screening.

### §2.2 T_P1_2 [FP] Vortex EL equation

**Lagrangian:** $\mathcal{L} = \frac{1}{2}|\nabla\Phi|^2 - \frac{\lambda}{4}(|\Phi|^2 - v^2)^2$

**Ansatz:** $\Phi_{vortex}(r, \phi) = \rho(r) \cdot e^{i n \phi}$ (z-translation invariant)

**Derived EL equation:**

$$\boxed{\rho''(r) + \frac{\rho'(r)}{r} - \frac{n^2 \rho(r)}{r^2} + \lambda \rho(r) (\rho(r)^2 - v^2) = 0}$$

**Boundary conditions:**
- $\rho(0) = 0$ (core regularity)
- $\rho(r \to \infty) = v$ (mexican hat VEV)

**Sympy verification:** $EL_{lhs} - EL_{expected} = 0$ exactly.

**Verdict:** ✓ PASS — well-posed BVP (boundary value problem).

### §2.3 T_P1_3 [FP] Far-field asymptotic

**Linearization** $\rho = v + \delta(r)$:

Vacuum residual at O(δ⁰): $-n^2 v/r^2$ (power-law source from angular winding)

First-order EL: $\delta'' + \delta'/r + (2\lambda v^2 - n^2/r^2)\delta - n^2 v/r^2 = 0$

**At large r (n²/r² → 0):**
- Homogeneous: $\delta'' + \delta'/r - m_\sigma^2 \delta = 0$ → $\delta_{hom}(r) \sim K_0(m_\sigma r) \sim e^{-m_\sigma r}/\sqrt{r}$
- Particular: $\delta_{part}(r) \approx -n^2/(2\lambda v r^2)$ (power-law)

**Full far-field structure:**

$$\rho(r) \approx v - \frac{n^2}{2\lambda v r^2} + A \cdot K_0(m_\sigma r) + O\left(\frac{1}{r^4}\right)$$

**Critical observation:** Far-field has BOTH power-law correction (from n² source) AND exponential correction (from massive Higgs mode). The presence of power-law structure in field profile **at single-vortex level** indicates 3D structure differs fundamentally from 1D Z2 (pure exp(-m·√2·L)).

**Verdict:** ✓ PASS — vacuum residual exactly $-n^2 v/r^2$; m_σ confirmed.

### §2.4 T_P1_4 [FP] Mass spectrum

**Decomposition:** $\Phi = (v + \sigma) \cdot e^{i\pi/v}$

Quadratic Lagrangian:

$$\mathcal{L}_{quad} = \frac{1}{2}(\partial\sigma)^2 - \lambda v^2 \sigma^2 + \frac{1}{2}(\partial\pi)^2$$

**Mass spectrum:**

| Mode | Mass² | Mass |
|------|-------|------|
| Higgs (σ) | $2\lambda v^2$ | $v\sqrt{2\lambda}$ |
| Goldstone (π) | 0 | **MASSLESS** |

**Goldstone theorem application:**
- Continuous global U(1) broken by ⟨Φ⟩ = v
- Expected Goldstone modes: 1
- Computed: 1 ✓

**Positivity:** m_σ² = 2λv² > 0 for λ, v > 0.

**Verdict:** ✓ PASS — mass spectrum exactly as analytically predicted (Phase 0 §8.3).

**CRITICAL implication dla Phase 2/3:** Goldstone π is genuinely massless w pełnym 3D TGP S05+Z₂. To **mediator** long-range interaction between defects → expected 1/L (point) lub log(L) (vortex) at large separation.

### §2.5 T_P1_5 [FP] Z₂ × RP² compatibility

**Z₂ check:** Under Φ → -Φ: vortex ρ·exp(inφ) → ρ·exp(inφ + iπ). Winding number n unchanged. **Z₂ compatible** ✓

**RP² check:** RP² topology relevant for σ_ab orientation field; pure Φ-vortex independent. **Default vortex compatible** ✓

**Note:** Hedgehog point defect (DEC alternative) WOULD use RP² explicitly. Vortex line jest simpler default.

**Verdict:** ✓ PASS — both symmetries preserved.

---

## §3 — §3.6 BINDING compliance check

**Per CALIBRATION §3.6:** Phase 0 included analytical pre-derivation; Phase 1 sympy VERIFIED predictions.

**Pre-registered predictions vs Phase 1 sympy results:**

| Prediction (Phase 0 §8) | Phase 1 sympy result | Match? |
|-------------------------|----------------------|--------|
| EL equation Nielsen-Olesen form | $\rho'' + \rho'/r - n^2\rho/r^2 + \lambda\rho(\rho^2-v^2)=0$ derived | ✓ exact |
| Mass spectrum $m_\sigma = v\sqrt{2\lambda}$ | $\sqrt{2\lambda} v$ derived | ✓ exact |
| Goldstone massless | $m_\pi = 0$ derived | ✓ exact |
| Far-field power-law correction $-n^2/(2\lambda v r^2)$ | Vacuum residual $-n^2 v/r^2$ → particular sol $-n^2/(2\lambda v r^2)$ | ✓ exact |
| Far-field exponential correction $e^{-m_\sigma r}/\sqrt{r}$ | $K_0(m_\sigma r)$ asymptotic $\sim e^{-m_\sigma r}/\sqrt{r}$ | ✓ exact |

**All analytical predictions verified.** §3.6 BINDING enforcement successful: NO discrepancy between heuristic pre-derivation and sympy execution.

---

## §4 — Anti-Lakatos discipline check

- ✅ Pre-registered analytical predictions LOCKED Phase 0 (§3.6 BINDING)
- ✅ Each test reported per pre-registered claim (NO ex post threshold modifications)
- ✅ Sympy results match analytical predictions EXACTLY (no rescue needed)
- ✅ 0 hardcoded T_pass=True (all 4 FP tests compute-then-compare via symbolic verification)
- ✅ DEC budget preserved 0/1 (vortex default, NOT DEC switch)
- ✅ Z₂ × RP² compatibility verified honestly (NIE force PASS by relaxing definition)

---

## §5 — Implications dla Phase 2

**Key Phase 1 outputs LOCKED:**
- Mass spectrum confirmed: $m_\sigma = v\sqrt{2\lambda}$ (massive), $m_\pi = 0$ (massless)
- Single-vortex profile structurally well-defined
- Goldstone π is genuine long-range mediator (3D propagator $G_\pi(r) = 1/(4\pi r)$)

**Phase 2 task:** Two-vortex interaction V_int(L) at large separation, with Goldstone mediating long-range component.

**Pre-registered expectation:** V_int(L)/L_z ∝ v²·n_1·n_2·log(L/r_0) (cosmic-string-style 2D log) **for parallel vortex lines** (z-translation invariant ansatz).

---

**END OF Phase 1 results — PROCEED_TO_PHASE_2 LOCKED 2026-05-23**
