---
title: "Phase 2 Results -- Parameter scan extended"
type: phase_results
status: CLOSED
phase: 2
parent_cycle: op-CE-H-two-particle-equilibrium-2026-05-21
date_completed: 2026-05-21
result: 5/5 PASS (after numerical seed hygiene fix), F-beta-3 + F-beta-4 fully verified
---

# Phase 2 Results -- Parameter scan extended

**Status:** CLOSED 2026-05-21
**Result:** 5/5 PASS (5/5 substantive FP)
**F-β-3, F-β-4 verification:** EXTENDED CONFIRMED across full (α, D, m) grid

---

## §1 — Test verdict table

| Test | Description | Result |
|------|-------------|--------|
| T_P2_1 | α=0.5: equilibrium exists, stable + unstable branches identified | PASS |
| T_P2_2 | α=2: equilibrium exists, stable + unstable branches identified | PASS |
| T_P2_3 | α=3: equilibrium exists, stable + unstable branches identified | PASS |
| T_P2_4 | L*(m) ∝ 1/m scaling (dimensional analysis) | PASS |
| T_P2_5 | F-β-4 robust across (α, D) 4×5=20 grid | PASS (20/20) |

**Substantive metrics:**
- 5/5 substantive FP PASS (100%)
- 0 hardcoded T_pass=True (strict cycle 1/2/7 preserved)
- 0/1 DEC budget used (cumulative Poziom β preserved unused)

---

## §2 — General analytical results (for arbitrary α)

### §2.1 Stationarity equation

For E_total(L) = 2·E_K - A·exp(-mL) + D/L^α:

$$\frac{dE}{dL} = A m e^{-mL} - \frac{\alpha D}{L^{\alpha+1}} = 0$$

Substituting u = m·L:

$$u^{\alpha+1} e^{-u} = \frac{\alpha D}{A m^{-\alpha}} \cdot \frac{1}{m} = \frac{\alpha D m^\alpha}{A}$$

### §2.2 Critical D_critical

Function g(u; α) = u^(α+1)·e^(-u) has interior maximum at u = α+1:

$$g_{max}(\alpha) = (\alpha+1)^{\alpha+1} \cdot e^{-(\alpha+1)}$$

Therefore:

$$D_{critical}(\alpha) = \frac{A \cdot (\alpha+1)^{\alpha+1} \cdot e^{-(\alpha+1)}}{\alpha \cdot m^\alpha}$$

### §2.3 Stability criterion

$$\frac{d^2 E}{dL^2}\bigg|_{L^*} = \frac{\alpha D}{L^{*(\alpha+2)}} \cdot \big[(\alpha+1) - m L^*\big]$$

**Stable iff m·L* < α+1**, i.e. u_stable < α+1.

### §2.4 Stability boundary verification (numerical)

| α | α+1 | u_stable@D=0.5·D_crit | u_unstable@D=0.5·D_crit | Stab indicator (stable) | Stab indicator (unstable) |
|---|-----|---|---|---|---|
| 0.5 | 1.5 | 0.4781 | 3.4367 | +1.0219 | -1.9367 |
| 1.0 | 2.0 | 0.6053 (Phase 1b) | 4.7079 (Phase 1b) | +1.395 | -2.708 |
| 2.0 | 3.0 | 1.3941 | 5.5254 | +1.6059 | -2.5254 |
| 3.0 | 4.0 | 2.0828 | 6.8379 | +1.9172 | -2.8379 |

All α values: u_stable < α+1 ✓ stable; u_unstable > α+1 ✓ unstable.

---

## §3 — Mass scaling (T_P2_4)

### §3.1 Dimensional analysis

Stationarity at fixed dimensionless ratio D·m^α/A:
- u_stable depends ONLY on this dimensionless combination
- L_stable = u_stable / m (natural length scale 1/m)

**Verified numerically:**

| m | u_stable (D·m/A = 0.3) | L* = u_stable / m |
|---|---|---|
| 0.5 | 0.8291 | 1.6581 |
| 1.0 | 0.8291 | 0.8291 |
| 2.0 | 0.8291 | 0.4145 |
| 4.0 | 0.8291 | 0.2073 |
| 8.0 | 0.8291 | 0.1036 |

**u_stable constant across m** (factor 16 range): ✓
**L* * m = u_stable** (perfect scaling): ✓

### §3.2 Physical interpretation

The natural length scale of equilibrium is 1/m (Compton-like scale of Phi fluctuations). For TGP at hadronic scale (m ~ Phi_0 ~ Λ_QCD ~ 200 MeV), L* ~ 1 fm (right order of magnitude for hadron size).

**Honest caveat:** this is **dimensional consistency**, not quantitative prediction. Mapping m → physical TGP scale requires Poziom γ.

---

## §4 — Full (α, D) grid robustness (T_P2_5)

### §4.1 Grid coverage

- α ∈ {0.5, 1.0, 2.0, 3.0} (factor 6 range)
- D/D_critical ∈ {0.1, 0.3, 0.5, 0.7, 0.9} (factor 9 range)
- Total: 20 (α, D) pairs

### §4.2 All 20 pass

After numerical seed hygiene (adaptive seed selection for nsolve):

**All 20 pairs:**
- Stable branch exists with u_s < α+1
- Unstable branch exists with u_u > α+1
- d²E/dL² sign correct on each branch

**F-β-4 ROBUST CONFIRMED:** equilibrium exists across full grid without fine-tuning.

### §4.3 Numerical hygiene note (HONEST METHODOLOGICAL DECLARATION)

**Initial run:** 18/20 PASS. Two failures at (α=0.5, D/D_crit ∈ {0.1, 0.3}).

**Root cause:** nsolve TypeError z fixed seed (α+1)/2 = 0.75 dla małego u_stable case. Function u^1.5 ma derivative-singular behavior near u=0 dla α=0.5, co utrudnia convergence Newton's method.

**Fix:** adaptive seed selection — try multiple seeds, pick first that converges to correct branch:
- Stable branch seeds: (α+1)·sqrt(D_fraction), (α+1)·D_fraction, (α+1)·0.5, (α+1)·0.1, 0.01
- Unstable branch seeds: (α+1)·1.5+1, (α+1)·2, (α+1)·3, (α+1)·(2-D_fraction)

**Anti-Lakatos check:**
- ✗ NOT modifying F-β-4 threshold (still 20/20 required)
- ✗ NOT modifying pre-registered prediction
- ✓ Fixing NUMERICAL METHOD (better seeds), NOT what we're testing
- ✓ Analytical structure unchanged: equilibrium PROVABLY exists for D < D_crit at any α > 0

**Honest declaration:** initial 18/20 failure was a **numerical artifact**, not structural. After standard numerical hygiene (multi-seed retry), 20/20 PASS.

---

## §5 — What Phase 2 adds beyond Phase 1a/1b

### §5.1 Robustness of CE-H mechanism

Phase 1b verified equilibrium for **one** (α=1, m=1, A_int=1) configuration. Phase 2 extends:
- Multiple α (Coulomb-like 1, sub-Coulomb 0.5, super-Coulomb 2-3)
- Dimensional scaling verified
- Grid coverage 20 (α, D) combinations all pass

**Implication:** CE-H mechanism is **not fragile** to specific functional form choice of D/L^α. The structural feature (bg-stabilized equilibrium) persists across broad family of long-range repulsion forms.

### §5.2 Critical confinement-deconfinement boundary

For each α, D_critical = A·(α+1)^(α+1)·e^(-(α+1)) / (α·m^α).

| α | D_crit/A·m^α | u at critical (α+1) |
|---|-----|-----|
| 0.5 | 0.41 | 1.5 |
| 1.0 | 0.54 | 2.0 |
| 2.0 | 1.34 | 3.0 |
| 3.0 | 4.69 | 4.0 |

**Larger α (steeper repulsion) → higher D_critical** (more bg coupling tolerated before deconfinement).

**Physical implication for Poziom γ:** in full TGP with α derived from substrate Lagrangian, D_crit determines confinement/deconfinement transition. This could connect to QCD phase diagram (T_c) at full cosmological scale.

**Not pursued in Poziom β** — noted for Poziom γ extension.

---

## §6 — F-β verifications summary (Phase 1a + 1b + 2)

| FP | Phase 1a | Phase 1b | Phase 2 |
|----|----------|----------|---------|
| F-β-1 NULL isolation | ✓ CONFIRMED | n/a | n/a |
| F-β-2 POSITIVE with bg | n/a | ✓ CONFIRMED (α=1) | ✓ EXTENDED to α∈{0.5,1,2,3} |
| F-β-3 monotonic L*(D) | n/a | ✓ CONFIRMED (factor 10) | ✓ EXTENDED (factor 9 × 4 α) |
| F-β-4 no fine-tuning | n/a | ✓ CONFIRMED (factor 10) | ✓ EXTENDED (20/20 grid) |
| F-β-5 self-consistency closure | deferred | deferred | deferred → Phase 3 |

---

## §7 — Honest caveats persistent

### §7.1 α is exogenous

- α value chosen by hand (Phase 1b: α=1; Phase 2: α ∈ {0.5, 1, 2, 3})
- Full TGP-native derivation of α from Phi-substrate Lagrangian = **Poziom γ scope**
- For Poziom β toy: any reasonable α gives equilibrium, so result is **robust** to α choice within tested range

### §7.2 D/L^α form remains constructed

- The repulsive bg contribution is **modeled**, not derived from substrate self-consistency
- Phase 3 will partially address via (EQ-1)↔(EQ-2) self-consistency check
- Full derivation remains Poziom γ scope

### §7.3 1D Z2 limitations persist

- Full TGP is 3D with U(1) + RP² topology
- Vortex/hedgehog dynamics in 3D may have additional structure not captured by 1D kinks
- Poziom β = structural proof-of-principle at simplest level

---

## §8 — Discipline status

### §8.1 Anti-Lakatos

- ✅ Initial 18/20 result reported honestly
- ✅ Investigation identified numerical (not structural) cause
- ✅ Fix is NUMERICAL METHOD (adaptive seeds), NOT threshold modification
- ✅ Post-fix 20/20 result transparent, fully documented

### §8.2 Strict cycle 1/2/7

- ✅ 0 hardcoded T_pass=True maintained
- ✅ All substantive FP use compute-then-compare
- ✅ DEC budget 0/1 still preserved

### §8.3 Native equations

- ✅ Built on Phase 1a/1b results
- ✅ No fitting to QCD/SM (despite tempting QCD confinement analogy noted for Poziom γ)

---

## §9 — Status końcowy Phase 2

- ✅ 5/5 tests PASS (5/5 substantive)
- ✅ F-β-3, F-β-4 extended across (α, D, m) grid
- ✅ Cumulative Poziom β: 14/14 substantive PASS (100%)
- ✅ Discipline LOCKED across Phase 1a + 1b + 2

**Phase 2 CLOSED 2026-05-21. Phase 3 (self-consistency closure F-β-5) authorized for execution.**
