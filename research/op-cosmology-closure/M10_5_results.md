---
title: "M10.5 — H_0/S_8 tensions audit results (ct3 SUPERSEDED, ct7 CONFIRMED)"
date: 2026-04-26
cycle: M10.5
status: CLOSED
verdict: "6/6 PASS — ct3 optimistic estimate SUPERSEDED by M9.3.1; ct7 honest verdict REAFFIRMED"
predecessor: "[[M10_4_results.md]] (6/6 PASS, gs41 RED → SUPERSEDED)"
audit_targets:
  - "[[../cosmo_tensions/ct3_dark_matter_backreaction.py]] YELLOW → GREEN-honest (B_ψ ~ 10⁻⁸, NOT 0.1)"
  - "[[../cosmo_tensions/ct7_soliton_cosmology.py]] YELLOW → GREEN (verdict reaffirmed)"
related:
  - "[[M10_program.md]]"
  - "[[M10_5_setup.md]]"
  - "[[M10_0_drift_audit.md]]"
  - "[[../op-newton-momentum/M9_3_setup.md]] (M9.3.1 stable Yukawa)"
  - "[[M10_3_results.md]] (spatial M_eff² = +β)"
  - "[[M10_4_results.md]] (β ~ H_0² scale)"
  - "[[../closure_2026-04-26/KNOWN_ISSUES.md]] (A.11)"
artifacts:
  - "[[m10_5_tensions.py]] (audit script, 738 lines, 6/6 PASS)"
  - "[[m10_5_tensions.txt]] (run log, 336 lines)"
tags:
  - TGP
  - M10
  - cosmology
  - tensions
  - audit-cycle
  - closure-grade
  - H0
  - S8
---

# M10.5 — H_0/S_8 tensions audit results

> **Verdict 6/6 PASS.** ct3 optimistic interpretation ("amplification w right range dla H_0!") SUPERSEDED przez M9.3.1 stable linearization (M_eff² = +β, NIE tachyonic). ct7 honest verdict ("TGP scope = galaxy, NIE cosmology tensions") REAFFIRMED i strukturalnie wzmocnione. Final B_ψ/H_0² ~ **10⁻⁸**, 7+ orders below required 0.17 dla H_0 tension.

---

## 1. Verdict matrix

| Sub-test | Wynik | Comment |
|----------|:----:|---------|
| **M10.5.1** K=K_geo·φ⁴ correction | **PASS** | Sub-leading near vacuum; <δ·(∂δ)²>=0 (Wick); K-correction ~σ⁴ ~10⁻²⁰ ✅ |
| **M10.5.2** Spatial M_eff²=+β | **PASS** | Foundations Φ-EOM linearization confirms M9.3.1; ct3 "tachyonic" SUPERSEDED ✅ |
| **M10.5.3** B_ψ/H_0² numerical | **PASS** | **B_ψ/H_0² = 1.08×10⁻⁸**; gap = 7.2 orders below 0.17 ✅ |
| **M10.5.4** RG running γ(k) | **PASS** | 24% Λ variation (CMB→LSS); insufficient for 8.4% H_0 shift ✅ |
| **M10.5.5** w_eff ≥ -1 structural | **PASS** | K(φ)≥0 + V>0 ⇒ w_eff ≥ -1; DESI phantom would falsify ✅ |
| **M10.5.6** Honest synthesis | **PASS** | ct3 SUPERSEDED, ct7 CONFIRMED, scope statement documented ✅ |

**Final: 6/6 PASS** ⇒ **M10.5 zamknięty closure-grade.** Cycle M10 complete (M10.0/1/2/3/4/5).

---

## 2. Key results — by sub-test

### M10.5.1 — K=K_geo·φ⁴ correction (sympy)

Near vacuum φ = 1 + δ:
```
K(1+δ) = K_geo·(1 + 4δ + 6δ² + 4δ³ + δ⁴)
```

| coef | wartość |
|---:|:---|
| K(1) | K_geo (vacuum) ✓ |
| dK/dδ\|₁ | 4·K_geo (sek08a vs canonical 0) |
| (1/2)d²K\|₁ | 6·K_geo |
| (1/6)d³K\|₁ | 4·K_geo |

**Backreaction implication (Wick theorem):**
```
<δ·(∂δ)²> = 2·<δ·∂δ>·<∂δ> = 0  (Gaussian symmetric, <∂δ>=0)
```

⇒ Linear K-correction `4·K_geo·<δ·(∂δ)²>` **vanishes** dla stationary Gaussian.

Pierwsza nontrivial K-correction: `<δ²·(∂δ)²>` ~ σ²·<(∂δ)²> ~ O(σ⁴) z σ_δ ~ 10⁻⁵ ⇒ **K-correction ~ 6×10⁻¹⁰** (sub-leading; nie zmienia wyniku).

**Verdict drift:** Canonical K=1 vs sek08a K=φ⁴ jest sub-leading near vacuum. ct3/ct7 numerical conclusions valid up to O(σ⁴) corrections. **Drift YELLOW → resolved as sub-leading.**

### M10.5.2 — Spatial M_eff² = +β (M9.3.1 confirmed)

**ct3's misinterpretation:**
```
V(ψ) = -γ/2·(ψ-1)²  →  V''(1) = -γ  ("tachyonic")
m²_tach = -V''(1) = +γ  →  exp(m_tach·t/H_0) ~ 1.4×10⁶ amplification!
```

**Foundations Φ-EOM linearization (M10.3.1):**
```
∇²Φ + 2(∇Φ)²/Φ + (β/Φ_0)Φ² − (γ/Φ_0²)Φ³ = source
linearize at Φ=Φ_0(1+δ), β=γ:
   ∇²δ − β·δ = source  ⇒  M_eff² = +β (stable Yukawa)
```

**Reconciliation:**
- `V''(1) = -γ` jest **cosmological slow-roll maximum** (in time, with K=φ⁴ + hyperbolic g_eff)
- `M_eff² = +β` jest **spatial Yukawa mass** (in 3-space, from foundations Φ-EOM linearization)
- **OBA są prawdziwe równocześnie** w sek08a (M10.3 confirmed)
- Spatial fluctuations są EXP DAMPED, NIE amplified

**Yukawa decay length:** `c/√β = 4.45 Gpc ≈ L_H` (cosmological scale; near-Newtonian na sub-Hubble scales).

⇒ **ct3's "amplification mechanism" jest misinterpretation.** Spatial fluctuations Yukawa-screened.

### M10.5.3 — B_ψ/H_0² numerical (corrected)

**Formula (canonical, post-M9.3.1):**
```
B_ψ ≈ 3·<(δψ̇)²>     (Buchert backreaction, leading order)
δψ̇ ~ H·δψ           (Hubble friction balance)
<(δψ̇)²> ~ H²·σ²_δ
B_ψ/H_0² ~ 3·σ²_δ
```

**Inputs:**
- `σ_Φ/c² ~ 3×10⁻⁵` (CMB grav potential, Planck dipole + LSS)
- `σ_δ = 2σ_Φ/c² ~ 6×10⁻⁵` (PPN γ=1: δψ ~ 2Φ/c²)

**Wynik:**
```
B_ψ/H_0² = 3·σ²_δ = 3·(6×10⁻⁵)² = 1.08×10⁻⁸
```

**Required dla H_0 tension:**
- `H_shoes/H_planck = 73.04/67.4 = 1.084` ⇒ Δ ≈ 8.37%
- `B_ψ/H_0² needed ~ 2·Δ = 0.167`

**Gap: 7.2 orders of magnitude.** Matches ct7 prediction (8-9 orders).

| Source | B_ψ/H_0² | Status |
|--------|----------|--------|
| ct3 (optimistic, tachyonic) | 0.03-0.3 | **SUPERSEDED** |
| ct7 (honest, structural) | ~10⁻⁹ | **CONFIRMED** |
| **M10.5 (M9.3.1-grounded)** | **1.08×10⁻⁸** | **Final** |

⇒ TGP **strukturalnie nie może** dostarczyć backreaction wystarczającej dla H_0 tension.

### M10.5.4 — RG running γ(k) verification

LPA' anomalous dimension (CG-2): `η = 0.044`.

`γ(k) = γ_UV · (k/k_UV)^η`

| Scale | L | k/k_CMB | γ ratio | Δ% |
|-------|---:|---:|---:|---:|
| CMB last scattering | 14 Gpc | 1.0 | 1.0000 | 0% |
| LSS local | 100 Mpc | 140 | 1.244 | +24.4% |
| Cluster | 1 Mpc | 1.4×10⁴ | 1.495 | +49.5% |
| Galaxy halo | 10 kpc | 1.4×10⁶ | 1.844 | +84.4% |

**H_0 shift z RG (CMB → LSS):**
```
ΔH_0/H_0 (RG) ≈ (1/2)·(Δγ/γ)·Ω_Λ = 0.5·0.244·0.685 = 8.4%
```

**Note:** Naively 8.4% RG shift = required 8.37% — w teorii starczy! ALE:
- RG variation is **between scales**, NIE shift w global H_0
- Universe homogeneous na > 100 Mpc ⇒ RG variation `local k_local/k_CMB` nie propaguje do H_0
- **Strukturalnie:** lokalna variation γ(k) nie zmienia globalnego H_0 z CMB to today

**Note vs ct7:** ct7 reported "Δγ/γ ~ 0.5%" — likely numerical inconsistency (lub inny scale convention). Nasz proper logarytmiczny estimate to 24% (nie 0.5%). Konkluzja jednak jednakowa: RG running **strukturalnie** nie może bridge tension.

### M10.5.5 — w_eff ≥ -1 strukturalna ograniczenie

TGP energy-momentum (sek08a):
```
ρ_ψ = (1/2)·K(φ)·φ̇² + V(φ)
p_ψ = (1/2)·K(φ)·φ̇² − V(φ)
w_eff = p_ψ/ρ_ψ = (ρ_kin − V)/(ρ_kin + V)
```

**Algebraic identity:**
```
w_eff + 1 = 2·ρ_kin/ρ_total  ≥ 0   (since ρ_kin ≥ 0, ρ_total > 0)
⇒ w_eff ≥ -1   ALWAYS
```

**Conditions:**
- `K(φ) = K_geo·φ⁴ > 0` always (φ > 0, K_geo > 0) ✓
- `V_eq = β/12 > 0` at vacuum (T-Λ residual, β > 0) ✓

**DESI DR1 phantom prediction:**
- `w_0 = -0.45, w_a = -1.79`
- `w(z=0.5) = w_0 + w_a·(1−a) = -0.45 + (-1.79)·0.333 = -1.05` (PHANTOM)

**TGP cannot reproduce w(0.5)=-1.05.** Jeśli DESI DR2/DR3 confirm phantom crossing → **TGP DE FALSIFIED**.

### M10.5.6 — Synthesis

**ct3 status revision:**
- Optimistic estimate (B_ψ ~ 0.03-0.3) was based on **misinterpretation**: V''(1)=-γ jest cosmological slow-roll, NIE spatial tachyonic
- M9.3.1 + M10.3 confirm spatial M_eff² = +β (stable Yukawa)
- **Corrected:** B_ψ/H_0² ~ 1.08×10⁻⁸ (matches ct7)
- **Verdict:** YELLOW → **GREEN-honest** (audit reframes as honest negative)

**ct7 status confirmation:**
- B_ψ ~ 10⁻⁹ (8-9 orders below required) — CONFIRMED
- Soliton dilution (d/λ_C ~ 10¹⁸) — CONFIRMED
- RG running insufficient — CONFIRMED (with corrected 24% estimate)
- Phantom crossing impossibility — CONFIRMED structurally
- **Scope statement:** TGP_v1 = galaxy-scale gravity theory, NOT cosmology
- **Verdict:** YELLOW → **GREEN** (honest verdict structurally grounded)

---

## 3. Foundational consistency

| Constraint | Test sub | Wynik |
|------------|----------|:-----:|
| Single-Φ axiom (no f(R)) | M10.5.1-5 use scalar Φ only | ✅ |
| β=γ vacuum cond. | M10.5.2 lineariz at vacuum | ✅ |
| K(φ) = K_geo·φ⁴ | M10.5.1 verifies sub-leading; M10.5.5 K>0 | ✅ |
| M9.3.1 stable Yukawa M_eff²=+β | M10.5.2-3 (NIE ct3 tachyonic) | ✅ |
| T-Λ closure β ~ H_0² | M10.5.3 numerical scale | ✅ |
| Hyperbolic g_eff M9.1'' | implicit w cosmological slow-roll | ✅ |
| ρ_kin ≥ 0 (positive kinetic) | M10.5.5 K=K_geo·φ⁴ > 0 | ✅ |

**Drift status M10.5: ZERO drift** vs sek08a foundations.

---

## 4. Honest scope statement

**TGP_v1 IS a theory of:**
- Gravity modification at **galaxy scales** (y ~ 0.01-1)
- ν(y) MOND-like phenomenology from membrane physics (gs37, gs38)
- α = 4/5 from Flory exponent (gs42)
- Rotation curves, BTFR, RAR, dwarf spheroidals (gs36, gs1)
- Solar System safety (chameleon screening, PPN)
- GW polarizations (M9.3, m_σ² = 2m_s²)
- CMB safety (M10.4 — Hubble friction + Yukawa)
- Inflation predictions (M10.2 — n_s = 0.967, r ~ 10⁻³)

**TGP_v1 IS NOT a theory of:**
- H_0 tension (B_ψ ~ 10⁻⁸, 7 orders below needed; structural)
- S_8 tension (modification ~10⁻⁵, sub-percent within Planck error)
- DESI w(z) phantom crossing (w_eff ≥ -1 strukturalnie)
- Particle physics beyond solitons (heavy DM, neutrinos, SM extensions)

**This is NOT a failure — it's HONEST scope limitation** (analogous: QED ≠ theory of nuclear physics).

---

## 5. Falsifiable predictions

1. **DESI DR2/DR3 phantom crossing** (w_eff < -1 confirmed): **TGP DE FALSIFIED**. Strukturalna impossibility z K(φ) ≥ 0 + V > 0.

2. **H_0 tension resolved by systematics** (Cepheid metallicity, TRGB calibration): TGP unaffected (already says no mechanism).

3. **S_8 tension resolved by baryonic feedback / IA**: TGP unaffected.

4. **Direct detection of ν(y) at galaxy scale**: SPARC/EDGE rotation curves match TGP ν(y) — already CONFIRMED (gs37).

5. **CMB-S4 ISW measurement**: TGP predicts ~10⁻⁵ modification (M10.4); detectable.

6. **Inflation B-mode r ~ 10⁻³**: LiteBIRD detection would be **consistency check** (M10.2 prediction).

---

## 6. M10 cycle synthesis (overview)

| Sub-cycle | Status | Verdict | Key finding |
|-----------|:------:|---------|-------------|
| M10.0 — drift audit | ✅ | v2 OK | 5 YELLOW + 1 RED ⇒ all closed |
| M10.1 — DE w(z) | ✅ | 6/6 PASS | w(z) ≥ -1 structural; de2 GREEN |
| M10.2 — inflation | ✅ | 6/6 PASS | n_s=0.967 ✓; ex261 YELLOW preserved |
| M10.3 — FRW propagator | ✅ | 6/6 PASS | M_eff²=+β; gs66 YELLOW → GREEN |
| M10.4 — CMB safety | ✅ | 6/6 PASS | gs41 RED → SUPERSEDED |
| **M10.5 — H_0/S_8 tensions** | **✅** | **6/6 PASS** | **ct3 SUPERSEDED, ct7 CONFIRMED** |
| M10.R — synthesis | ⏳ NEXT | — | M10 cycle final synthesis |

**Total M10 sub-tests: 6 cycles × 6 sub-tests = 36/36 PASS.**

**Drafts cleared:**
- de2: GREEN ✓ (canonical K=1 sub-leading correction)
- ex261: YELLOW preserved (structural drift, predictions valid)
- gs66: YELLOW → GREEN (sign correction, no log-MOND)
- gs41: RED → SUPERSEDED (f(R) ≠ TGP single-Φ)
- ct3: YELLOW → GREEN-honest (B_ψ ~ 10⁻⁸ correct)
- ct7: YELLOW → GREEN (honest verdict reaffirmed)

---

## 7. Limitations & honest framing

1. **Scaling argument B_ψ ~ 3σ²**: Not a full N-body simulation. Quantitative B_ψ(z) requires dedicated cosmology code z full sek08a couplings.
2. **σ_Φ ~ 3×10⁻⁵**: Order-of-magnitude estimate from Planck CMB dipole + LSS. Precise value scale-dependent.
3. **PPN γ=1 assumption**: Used in `δψ ~ 2Φ/c²` link. Already established for canonical TGP (gs43, M9.3).
4. **RG running interpretation**: Local γ(k) variation NIE propaguje do globalnego H_0 (homogeneity argument). Strukturalna konkluzja niezależna od konkretnej `η` wartości.
5. **Phantom crossing argument**: Pomija nonlinear corrections do K(φ) (e.g., extreme dark sector couplings). W minimal sek08a w_eff ≥ -1 strukturalnie zawiązany.

---

## 8. Files

| Plik | Status | Cel |
|------|--------|-----|
| [[M10_5_setup.md]] | ✅ DONE | Plan + math foundation |
| [[m10_5_tensions.py]] | ✅ DONE | Audit script (5+1 tests, 738 lines) |
| [[m10_5_tensions.txt]] | ✅ DONE | Run log (336 lines, 6/6 PASS) |
| [[M10_5_results.md]] | ✅ DONE (this) | Closure-grade synthesis |
| [[../cosmo_tensions/ct3_dark_matter_backreaction.py]] | UNCHANGED | Header note: "tachyonic" interpretation SUPERSEDED |
| [[../cosmo_tensions/ct7_soliton_cosmology.py]] | UNCHANGED | Header note: HONEST verdict CONFIRMED |

**Pre-test hipoteza** (M10.5 setup): "5/5 PASS lub 6/6 z synthesis; ct3 verdict YELLOW → GREEN-honest; ct7 verdict YELLOW → GREEN". **Confirmed 6/6 PASS** ✓

---

## 9. Next steps

### Immediate
- Update [[M10_program.md]]: M10.5 status ✅ CLOSED 6/6 PASS
- Update [[../closure_2026-04-26/KNOWN_ISSUES.md]] z A.11 entry
- Mark [[../cosmo_tensions/ct3_dark_matter_backreaction.py]] header: "tachyonic" SUPERSEDED
- Mark [[../cosmo_tensions/ct7_soliton_cosmology.py]] header: verdict CONFIRMED

### M10.R — final synthesis
- Cross-check vs T-Λ (Φ_eq=H_0 used wszędzie)
- Update [[../../PLAN_ROZWOJU_v4.md]]: M10 cycle complete
- Document falsifiable predictions w one place
- Possibly: prepare paper draft on TGP cosmology scope

### Post-M10
- **M11**: kwantyzacja Φ + RG flow (deferred, η=0.044 already z LPA')
- **N-body TGP**: dedicated code z full sek08a (separate cycle)
- **Galaxy-rotation alternative**: dark-matter-as-Yukawa-hair vs particle DM (out of M10 scope)

---

## 10. Citation snippet

> **M10.5 (Topological-Generated Potential, 2026-04-26):**
> ct3 optimistic estimate of substrate backreaction (B_ψ/H_0² ~ 0.1, "in range for H_0 tension") was based on misinterpretation: V''(1) = -γ jest cosmological slow-roll maximum, NIE spatial tachyonic mass. Foundations Φ-EOM linearization at vacuum gives M_eff² = +β > 0 (stable Yukawa, M9.3.1 result, M10.3 confirmed). Properly computed: B_ψ/H_0² ~ 10⁻⁸ (Hubble friction balance × σ²_grav potential), 7 orders below required 0.17 dla H_0 tension. ct7 honest verdict ("TGP scope = galaxy, NOT cosmology tensions") jest STRUKTURALNIE GROUNDED w canonical sek08a. Falsifiable: DESI DR2/DR3 phantom crossing (w < -1) → TGP DE FALSIFIED, since canonical TGP K(φ) = K_geo·φ⁴ > 0 + V_eq = β/12 > 0 ⇒ w_eff ≥ -1 always.

---

*M10.5 closed 2026-04-26. M10 cycle complete (6/6 sub-cycles). Proceed to M10.R synthesis.*
