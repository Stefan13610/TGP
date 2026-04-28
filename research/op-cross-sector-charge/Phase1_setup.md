---
title: "XS.1.Phase1 setup — dimensional + structural audit of √α₀ = κ_TGP"
date: 2026-04-28
cycle: XS.1.Phase1
status: PRE-EXECUTION
predecessor: "[[../op-bh-alpha-threshold/Phase3_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - cross-sector
  - alpha-0
  - kappa-TGP
  - dimensional-audit
  - falsification
---

# XS.1.Phase1 — Dimensional + structural audit of √α₀ = κ_TGP

> **Cel:** Sprawdzić, czy hipoteza tożsamości √α₀ = κ_TGP jest
> wewnętrznie spójna na poziomie wymiarów, danych źródłowych i
> statystycznego priora. Phase 1 nie udowadnia tożsamości — sprawdza
> czy ma sens **w ogóle ją badać** w Phase 2.

---

## Hipotezy

**H₀ (null, "coincidence"):** α₀ i κ_TGP² to dwie niezależne stałe; ich
zgodność na 0.75% to przypadkowy fit dwóch O(1) liczb.

**H₁ (TGP charge identity):** √α₀ = κ_TGP jest fundamentalną tożsamością
TGP — α₀ (BH photon-ring) i κ_TGP (SC pair-breaking) wynikają ze
wspólnego substrate-field generator pod TGP single-Φ axiom.

Phase 1 testuje **konieczne warunki** dla H₁ (dimension, data
independence, multi-anchor robustness, prior odds). Wszystkie 5
warunków musi PASS, by Phase 2 (substrate-action derivation) miała
strukturalne podstawy.

---

## Pre-Phase 1 stan

| Pole | BH side | SC side |
|---|---|---|
| Stała | **α₀ ≈ 4.018** (Phase 2 strict, ξ=1) | **κ_TGP = 2.012** (TGP-SC v2) |
| Pochodzenie | photon-ring shift target_shift = (1/2)(1 − 3/3.88) | T_c calibration V/Nb/Ta/Mo/Pd ensemble |
| Dataset | M9.2-D + Phase 1.B.2 + closure_2026-04-26/alpha_psi_threshold | TGP-SC v2 paper (5 BCS metals) |
| Precision | ~1% (ξ unresolved + 0.114 round-off) | ~0.5% (V/Nb/Ta/Mo/Pd RMS) |
| Status | PARTIALLY DERIVED (Phase 2.T2.4) | LOCKED (TGP-SC v2 calibrated) |

**Numerical match (BH.1.Phase2.T2.5):**
```
α₀ = 4.0179
κ_TGP² = 4.0481
|Δ|/κ_TGP² = 0.7464%
```

---

## 5 sub-testów Phase 1

### T1.1 — Dimensional audit α₀ vs κ_TGP² w TGP natural units

**Cel:** Sprawdzić, czy α₀ i κ_TGP² są w **tych samych** jednostkach
(obie dimensionless, lub obie z konkretnym wymiarem).

**Test:**
- α₀ enters α(ψ) = α₀(ψ−1)²; ψ dimensionless [GM/c²r ratio];
  α(ψ) = photon-ring fractional shift dimensionless → **α₀ DIMENSIONLESS**
- κ_TGP enters T_c formula (TGP-SC v2): T_c = T_c^base · F(λ_sf, κ_TGP, β);
  κ_TGP = pair-breaking coupling normalization → **κ_TGP DIMENSIONLESS**
- Therefore both α₀ and κ_TGP² are dimensionless O(1) constants

**Falsification:** if dimensions disagree (one carries [length]ⁿ or
[mass]ⁿ where the other does not) → **identity impossible**.

### T1.2 — Numerical match z uncertainty propagation

**Cel:** Sprawdzić, czy α₀ i κ_TGP² są zgodne **w obrębie złożonej
niepewności**.

**Test:**
- α₀ = 4.0179 ± 0.04 (1% precision: ξ uncertainty + Phase 2 round-off)
- κ_TGP² = 4.0481 ± 0.04 (0.5% on κ_TGP propagated)
- |Δ| = 0.0302
- Combined 1σ ≈ √(0.04² + 0.04²) ≈ 0.057
- Distance: |Δ|/σ_combined ≈ 0.53σ → **within 1σ agreement**

**Falsification:** if |Δ|/σ > 3 → identity rejected statistically.

### T1.3 — Multi-anchor consistency κ_TGP

**Cel:** Sprawdzić, czy κ_TGP jest **stabilny** across różnych anchor
materials (V/Nb/Ta/Mo/Pd) — niestabilność dyskwalifikuje identity test.

**Test:**
- TGP-SC v2 calibration: κ_TGP from each of V, Nb, Ta, Mo, Pd individually
- Compute mean ± std; require std/mean < 1% (sub-percent stability)
- For each anchor, also compute √α₀ implied (= κ_TGP from that anchor)
  and check consistency with α₀ ≈ 4.018

**Falsification:** if RMS spread of κ_TGP across anchors > 5% → κ_TGP
is anchor-dependent fit, not stable physical constant; identity
hypothesis meaningless.

### T1.4 — Data independence (BH vs SC datasets disjoint)

**Cel:** Sprawdzić, że α₀ (BH) i κ_TGP (SC) zostały skalibrowane na
**rozłącznych** datasets — brak shared fit-parameters wykluczający
"match by construction".

**Test:**
- α₀: derived from M9.2-D (photon-ring rate) + closure_2026-04-26/T-α
  (geometric scenario (e), shadow shift target 0.114)
- κ_TGP: calibrated from V/Nb/Ta/Mo/Pd T_c at standard pressure
- Datasets: SMBH photon-ring imaging (EHT) vs ambient-pressure
  superconductor T_c (NIST/PDG) — **completely disjoint physical sectors**
- No shared free parameter between BH and SC fits

**Falsification:** if any common fit parameter was used in both
calibrations → match is artifact, not structural prediction.

### T1.5 — Coincidence null hypothesis (Bayesian prior)

**Cel:** Estymować prior odds H₁ vs H₀ z czystej statystyki: jakie jest
prawdopodobieństwo że dwie niezależne O(1) stałe TGP zgadzają się w
0.75% bez strukturalnego powodu?

**Test:**
- Define "agreement at level p": |α₀ − κ_TGP²|/κ_TGP² < p
- Reference range for both constants: O(1)–O(10) (typical TGP
  structural constants like β=2.527, γ=π²/8, κ²=4.05, K_geo=1)
- Uniform prior over [1, 10] for each: P(agreement at 0.75%) for two
  independent draws ≈ 2 × 0.0075 / (10−1) ≈ 0.17%
- Better baseline: take 5 known TGP O(1) constants from registry
  (β, κ², γ, K_geo, m_σ²/m_s²=2) → 5C2 = 10 pairwise comparisons;
  count how many agree at 0.75% level
- Expected: < 1 pair (since 10 × 0.0017 ≈ 0.02 expected by chance)
- Observed: 1 pair (α₀ vs κ_TGP²) → **prior odds 50:1 in favor of
  structural origin**

**Falsification:** if Bayesian prior odds favor coincidence (P > 30%) →
match is mundane and not worth Phase 2 derivation.

---

## Verdict gate

5/5 PASS → identity hypothesis viable, **proceed to Phase 2**
(substrate-action derivation).

≤4/5 → revise hypothesis or close XS.1 as coincidence with documented
audit.

---

## Materiał wykonawczy

- **Skrypt:** [`phase1_cross_sector_audit.py`](phase1_cross_sector_audit.py) (5 sub-tests)
- **Output:** `phase1_cross_sector_audit.txt`
- **Memo:** `Phase1_results.md` (closure + Phase 2 decision)

## Środowisko

```bash
cd TGP/TGP_v1/research/op-cross-sector-charge
PYTHONIOENCODING=utf-8 python -X utf8 phase1_cross_sector_audit.py 2>&1 | tee phase1_cross_sector_audit.txt
```

---

## Constants used

```
alpha_0    = 4.0179      (BH.1.Phase2 strict, xi=1)
alpha_0_err = 0.04       (~1% precision)

kappa_TGP   = 2.0120     (TGP-SC v2, V/Nb/Ta/Mo/Pd RMS calibration)
kappa_TGP_err = 0.010    (~0.5% precision)
kappa_TGP_sq = 4.0481

beta = 2.527             (TGP-SC v2)
gamma = pi^2/8 = 1.2337  (TGP core)
K_geo = 1                (Phase 1.A K(phi) = K_geo * phi^4)
m_sigma_sq_over_m_s_sq = 2  (sympy-exact, GW4)

|Delta|/kappa_TGP^2 = 0.7464%
```

---

## Cross-references

- [`program.md`](program.md) — overall 3-phase XS.1 plan
- [`../op-bh-alpha-threshold/Phase2_results.md`](../op-bh-alpha-threshold/Phase2_results.md) — α₀ origin
- [`../op-sc-alpha-origin/Phase3_results.md`](../op-sc-alpha-origin/Phase3_results.md) — κ_TGP origin
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — BH8 entry (the same identity registered as falsification target)
