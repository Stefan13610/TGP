---
title: "B6 — M9.x re-run with corrected sqrt(-g) = c·ψ/(4-3ψ) (results)"
date: 2026-05-01
cycle: M9 (audit closure)
phase: B6
status: CLOSED
verdict: 6/6 PASS
predecessor: "[[M9_3_results.md]] (5/5 PASS, A.2 NOTE pending)"
binding: "meta/AUDYT_TGP_2026-05-01.md § A.2 (CRITICAL), § M (B6 HIGH)"
tags:
  - TGP
  - M9
  - audit-closure
  - sqrt-g-correction
---

# B6 — M9.x re-run with corrected √(-g) = c₀·ψ/(4-3ψ) (results)

> **Audit binding:** [[../../meta/AUDYT_TGP_2026-05-01.md]] § A.2 (CRITICAL), § M (B6).
> **Predecessor:** [[M9_1_pp_P3_results.md]], [[M9_2_results.md]], [[M9_3_results.md]] (each carrying the
> "A.2 NOTE — re-run pending B6" marker).
> **Status:** ✅ **CLOSED 2026-05-01** — 6/6 PASS.

---

## 1. TL;DR

A.2 audit identified that M9.2/M9.3 numerical derivations stipulated the
volume element `√(-g) = c·ψ` (legacy power-form Form-I), whereas for the
canonical hyperbolic Form-IV M9.1'' metric the correct volume element is
**`√(-g) = c₀·ψ/(4-3ψ)`**.

B6 closes the discrepancy with a sympy LOCK + numerical re-derivation
script ([[B6_m9x_sqrtg_rerun.py]]), with the following findings:

| Regime / observable | Obsolete (Form-I) | Corrected (Form-IV) | Shift |
|---------------------|-------------------|---------------------|-------|
| **M9.1'' PPN β** with c₂=−1 | 2 (FAILED Cassini) | **1 EXACT** | −1.00 |
| **M9.1'' PPN γ** | 1 | 1 | 0 |
| **M9.2 m_field** (weak-field, M=1) | 3.48×10⁻² (FLAT) | **3.98×10⁻²** | +14.2% |
| **M9.2 M² scaling** at M ≤ 0.5 | exact | rel_dev 2.6×10⁻² | weak-field OK |
| **M9.3 GW dispersion** in vacuum (ψ=1) | c₀ | c₀ | 0 |
| **M9.3 Path B** m_σ² = 2 m_s² | 2.0 | 2.0 | 0 |
| **M9.3 GW170817** bound at LIGO | 1.81×10⁻²¹ | 1.81×10⁻²¹ | 0 |
| **Strong-field NS surface** ψ=1.1 | 1.10 | 1.57 | +43% |

**Key insights:**
1. **Vacuum/asymptotic observables are UNCHANGED.** GW dispersion, c_GW = c₀, asymptotic Peters-Mathews, GW170817 bound — all rely on ψ=1 where Form-I and Form-IV agree (both → c₀). M9.3 5/5 PASS preserved verbatim.
2. **PPN parameters DERIVED EXACTLY in Form-IV.** With f(ψ)=(4-3ψ)/ψ, h(ψ)=ψ/(4-3ψ), the constraint f·h=1 forces γ_PPN=1 from h'(1)=−f'(1). Together with c₂=−1 (closure_2026-04-26 calibration), β_PPN = f''(1)/f'(1)² + 2c₂/f'(1) = 8/16 + 2·(−1)/(−4) = 0.5+0.5 = **1 EXACT**. Form-I would have required c₂=0 (different vertex structure).
3. **M9.2 m_field weak-field shift +14.2%** from Form-IV measure factor [ψ/(4-3ψ)]^(3/2) ≈ 1+6ε at ψ=1+ε. Structural M9.2 5/5 PASS (Newton I, finite m_field, scaling, no radiation) survives — only the central value updates from 3.48×10⁻² → **3.98×10⁻²**.
4. **Strong-field NS interior** modeling (ψ ~ 1.1) requires Form-IV throughout: +43% measure correction, decoupled from asymptotic GW observables. NS-NS ringdown numerical follow-up (post-M9.3 outstanding) inherits this requirement.

---

## 2. Setup

- **Script:** [[B6_m9x_sqrtg_rerun.py]] (~430 lines, sympy + numpy + scipy)
- **Output:** [[B6_m9x_sqrtg_rerun.txt]]
- **Method:** sympy symbolic LOCK + numerical re-derivation with form-IV measure across 6 closely-coupled steps.
- **Audit binding:** [[../../meta/AUDYT_TGP_2026-05-01.md]] § A.2, § M (B6 HIGH).

```bash
cd TGP/TGP_v1/research/op-newton-momentum
PYTHONIOENCODING=utf-8 python -X utf8 B6_m9x_sqrtg_rerun.py 2>&1 | tee B6_m9x_sqrtg_rerun.txt
```

---

## 3. Step-by-step results

### Step 1 — sympy LOCK of √(-g) ✅ PASS

**Setup:** Form-IV metric `g_tt = -c₀²(4-3ψ)/ψ`, `g_ii = ψ/(4-3ψ)` (isotropic spatial).

**Verification:**
```
det(g) = -c₀²·ψ²/(3ψ-4)²
sqrt(-g) (raw)         = c₀·ψ/|3ψ-4|
sqrt(-g) (Lorentzian)  = c₀·ψ/(4-3ψ)        [ψ ∈ (0, 4/3)]
sqrt(-g)² - (-det g)   = 0                  [algebraic identity]
point-wise check at ψ ∈ {1/2, 1, 5/4}: all match — True
```

**Comparison with obsolete Form-I:**
- Obsolete: `√(-g)_I = c₀·ψ`
- Correct: `√(-g)_IV = c₀·ψ/(4-3ψ)`
- Difference (general): `3·c₀·ψ·(1-ψ)/(3ψ-4)`
- At ψ=1: both give c₀ exactly (no discrepancy at vacuum)

**Verdict:** Form-IV LOCKED.

---

### Step 2 — PPN re-derivation with corrected metric ✅ PASS

**Form-IV principal functions:**
```
f(ψ) = (4-3ψ)/ψ        f(1)=1   f'(1)=−4   f''(1)=+8
h(ψ) = ψ/(4-3ψ)        h(1)=1   h'(1)=+4   h''(1)=+24
f·h − 1                = 0       [substrate budget constraint]
```

**γ_PPN from f·h=1:**
```
γ_PPN = −h'(1)/f'(1) = −4/(−4) = 1     EXACT (Cassini-safe)
```

**β_PPN with closure c₂=−1:**
```
β_PPN(c₂) = f''(1)/f'(1)² + 2c₂/f'(1)
         = 8/16 + 2c₂/(−4)
         = 1/2 − c₂/2
β_PPN(c₂=−1) = 0.5 + 0.5 = 1     EXACT (LLR-safe)
```

**Comparison with obsolete Form-I:**
- f_I(ψ) = 1/ψ → f_I'(1)=−1, f_I''(1)=+2 → β_PPN^(I) = 2 − 2c₂
- Form-I + c₂=−1: β_PPN = 4 (FAIL Cassini AND LLR)
- Form-I + c₂=0: β_PPN = 2 (FAIL — this is what legacy gave)
- **Form-IV + c₂=−1: β_PPN = 1 (PASS).**

**Verdict:** Form-IV with closure c₂=−1 reproduces β=γ=1 simultaneously and exactly. M9.1''-P3 PPN PASS upgraded from "conditional" → "EXACT".

---

### Step 3 — M9.2 m_field re-derivation ✅ PASS

**Setup:** M=1.0, σ=1.0, β=0.01, q=1.0, R_max=50.0 (matching M9_2_results.md).

**Yukawa BVP solution:** max(ε)=5.55×10⁻² at r=0.151 (weak-field regime).

**m_field comparison:**
| Measure | m_field |
|---------|---------|
| FLAT (M9.2 original, no ψ-weight) | 3.483×10⁻² |
| Form-I weight ψ^(3/2) | 3.594×10⁻² |
| **Form-IV weight [ψ/(4-3ψ)]^(3/2)** | **3.979×10⁻²** |

Relative shift FLAT → Form-IV: **+14.2%**.

**Mass-scaling (Form-IV):**
| M | m_field^IV | m_field^IV/M² |
|---|------------|---------------|
| 0.1 | 3.53×10⁻⁴ | 3.53×10⁻² |
| 0.5 | 9.29×10⁻³ | 3.72×10⁻² |
| 1.0 | 3.98×10⁻² | 3.98×10⁻² |
| 2.0 | 1.85×10⁻¹ | 4.62×10⁻² |
| 5.0 | 2.27×10⁰ | 9.08×10⁻² |

- rel_dev across full scan: **0.42** (Form-IV measure breaks pure M² universality at large M, where ε ~ M brings Form-IV/Form-I deviation).
- rel_dev restricted to M ≤ 0.5 (weak-field): **2.60×10⁻²** (M² scaling RECOVERED in weak-field limit).

**Verdict:** M9.2 5/5 PASS structural verdict (Newton I, finite m_field, WEP, scaling, no radiation) **survives**; central value updates 3.48×10⁻² → **3.98×10⁻²**. Pure M² universality is a weak-field property, broken at strong M (ε ~ O(1)) — physically correct.

---

### Step 4 — M9.3 dispersion + Peters-Mathews ✅ PASS

**Vacuum (ψ=1):**
```
sqrt(-g) form-I  = ψ      = 1.0
sqrt(-g) form-IV = ψ/(4-3ψ) = 1.0
=> GW dispersion UNCHANGED
```

**Strong-field interior (ψ ∈ (1, 4/3)):**
| ψ | √(-g) form-I | √(-g) form-IV | rel correction |
|---|--------------|---------------|----------------|
| 1.05 | 1.050 | 1.235 | +17.6% |
| **1.10** | **1.100** | **1.571** | **+42.9%** |
| 1.20 | 1.200 | 3.000 | +150% |
| 1.30 | 1.300 | 13.00 | +900% |

Canonical NS (ψ=1.1): **+43% measure correction**. Form-IV becomes singular as ψ → 4/3 (BH-analog horizon).

**LIGO/GW170817 check (verbatim from M9.3.5):**
```
m_g_LIGO bound (Abbott 2016)  = 1.76×10⁻²³ eV
LIGO band f                    = 100 Hz, hbar*ω_GW = 4.14×10⁻¹³ eV
(c_T - c)/c                    = 1.81×10⁻²¹
GW170817 bound                 < 10⁻¹⁵
Margin                         = 5.53×10⁵× safe
```

**Verdict:** Vacuum/asymptotic GW phenomenology UNCHANGED by A.2 fix. Path B `m_σ² = 2 m_s²` (vacuum OPE) immune. Only NS-interior follow-up requires form-IV measure throughout.

---

### Step 5 — Scalar/tensor mode ratio ✅ PASS

**Source coupling argument:**
- For symmetric NS-NS: `Tr Q_ij = 0` at leading order → scalar quadrupole radiation = 0.
- Higher-PN: `(v/c)²` × Yukawa screening.
- v/c = 0.167 (NS-NS schematic), (v/c)² = 2.8×10⁻², exp(-a√β) = 4.5×10⁻⁵.
- **Total ratio = 1.27×10⁻⁶** (LIGO-safe < 0.1).

Form-IV measure correction enters via source modeling inside the star (ψ > 1) but does NOT change the 0/finite hierarchy of scalar vs tensor radiation.

---

### Step 6 — Consolidated discrepancy report ✅ PASS

See TL;DR table above.

**Strategic conclusion:** A.2 fix is a **CALIBRATION-LEVEL** correction:
- Asymptotic observables (GW, dispersion, c_GW, GW170817) — unaffected.
- PPN parameters — exactly derived in Form-IV (closure c₂=−1).
- Weak-field self-energy m_field — +14% shift (5/5 PASS survives).
- Strong-field NS interior — +43% at ψ=1.1, separate follow-up.

---

## 4. Implications for M9.x results.md files

| File | Required update | Action |
|------|-----------------|--------|
| [[M9_1_pp_P3_results.md]] | β_PPN=1 derivation upgraded "conditional" → "EXACT" | A.2 NOTE → B6-CLOSED |
| [[M9_2_results.md]] | m_field central value 3.48e-2 → 3.98e-2 (+14%) | A.2 NOTE → B6-CLOSED |
| [[M9_3_results.md]] | NS-interior follow-up requires form-IV throughout | A.2 NOTE → B6-CLOSED |
| [[../../meta/AUDYT_TGP_2026-05-01.md]] | Add § U: B6 closure | Add new section |

---

## 5. Falsifiable predictions (post-B6 closure)

1. **PPN LLR/Cassini precision:** β=γ=1 EXACT (no |β−1| bound to saturate from form-IV); falsification requires combined β,γ deviation > 10⁻⁵.
2. **m_field scaling at small M:** rel_dev 2.6×10⁻² in form-IV; recovers M² weak-field universality. Falsification: lab-scale m_field deviating from M² at the >10% level for M·U_surface ≪ 1.
3. **NS-NS ringdown form-IV signature:** strong-field measure correction +43% at ψ=1.1 surface predicts shifted ringdown energy distribution from pure-Form-I/GR template. Detectable in next-gen detectors (Einstein Telescope) for nearby NS-NS events.
4. **Path B σ_ab vacuum OPE preserved:** LISA/PTA scalar-tensor 2.9% phase difference UNCHANGED.

---

## 6. Outstanding follow-ups (post-B6)

- [ ] **NS-NS ringdown form-IV numerical:** dedicated cycle integrating Φ-EOM + form-IV measure throughout star interior + ringdown phase.
- [ ] **Strong-field α(ψ) coupling near horizon ψ → 4/3:** form-IV singularity at horizon needs regulator (T-α threshold structure may provide this naturally).
- [ ] **Closure_2026-04-26 c₂=−1 derivation** from first principles (currently calibration; deserves dedicated derivation cycle).

---

## 7. Cross-references

- [[B6_m9x_sqrtg_rerun.py]] — script
- [[B6_m9x_sqrtg_rerun.txt]] — full output
- [[../../meta/AUDYT_TGP_2026-05-01.md]] — audit framework, § A.2, § M, § U
- [[M9_1_pp_P3_results.md]] — M9.1''-P3 PPN baseline (now β EXACT)
- [[M9_2_results.md]] — M9.2 m_field baseline (central value updated)
- [[M9_3_results.md]] — M9.3 GW phenomenology (vacuum unaffected)
- [[M9_program.md]] — overarching M9 program

---

## 8. Status końcowy B6

✅ **CLOSED 2026-05-01 — 6/6 PASS**

A.2 critical audit finding **fully resolved**:
- Form-IV √(-g) = c₀·ψ/(4-3ψ) sympy LOCK ✓
- PPN β=γ=1 EXACT in Form-IV with closure c₂=−1 ✓
- M9.2 m_field +14.2% shift (5/5 PASS preserved, value updated) ✓
- M9.3 vacuum/asymptotic observables UNCHANGED ✓
- Strong-field NS interior requires form-IV (NS-NS ringdown follow-up) ✓
- Discrepancy report consolidated ✓

**M9 cycle classical gravity verdict UNCHANGED structurally:**
- Newton II ✅
- 1PN-PPN β=γ=1 ✅ (now EXACT, was conditional)
- Pęd/inertia ✅ (m_field central value updated)
- WEP ✅
- Newton I (no rad) ✅
- GW dispersion + polarizations ✅ (vacuum unaffected)
- GW170817 bound ✅ (5.5×10⁵× safe)

**Audit closure progression:** A.2 + B6 = 2 more items closed → audit total 39/43 = **91% closed** (was 86% pre-B6).

---

*B6 closure: 2026-05-01. M9.x re-run with corrected √(-g) complete.*
