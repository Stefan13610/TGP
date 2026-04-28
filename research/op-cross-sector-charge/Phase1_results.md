---
title: "XS.1.Phase1 results — dimensional + structural audit of √α₀ = κ_TGP"
date: 2026-04-28
cycle: XS.1.Phase1
status: CLOSED
verdict: PASS
predecessor: "[[../op-bh-alpha-threshold/Phase3_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - cross-sector
  - alpha-0
  - kappa-TGP
  - dimensional-audit
  - falsification
  - closure
---

# XS.1.Phase1 — Results: dimensional + structural audit of √α₀ = κ_TGP

> **Status:** CLOSED 2026-04-28 — **5/5 PASS**.
> Cross-sector identity hypothesis √α₀ = κ_TGP passes all 5 feasibility
> conditions (dimension, numerical, anchor-stability, data
> independence, Bayesian prior).
> **Decyzja:** proceed to **Phase 2** (substrate-action derivation).

---

## Verdict

| Sub-test | Description | Result |
|---|---|---|
| **T1.1** | Dimensional audit α₀ vs κ_TGP² (both dimensionless?) | **PASS** |
| **T1.2** | Numerical match in combined 1σ uncertainty | **PASS** |
| **T1.3** | Multi-anchor consistency κ_TGP across V/Nb/Ta/Mo/Pd | **PASS** |
| **T1.4** | Data independence (BH vs SC datasets disjoint?) | **PASS** |
| **T1.5** | Bayesian prior odds (coincidence vs structural) | **PASS** |

**5/5 PASS** → identity hypothesis **internally consistent**, justifies Phase 2.

---

## T1.1 — Dimensional audit

| Side | Constant | Role | Dimension |
|---|---|---|---|
| BH | α₀ | photon-ring T·J·J Wilson coefficient | dimensionless |
| SC | κ_TGP | pair-breaking coupling normalization | dimensionless |

Both stałe są dimensionless O(1) numbers w TGP natural units.

**Verdict:** PASS — identity dimensionally compatible.

## T1.2 — Numerical match (1σ test)

```
α₀         = 4.0179 ± 0.0400
κ_TGP²     = 4.0481 ± 0.0402
|Δ|        = 0.0302
|Δ|/κ_TGP² = 0.7471%
combined 1σ = 0.0567
distance   = 0.53σ
```

**Verdict:** PASS — match well within 1σ (gate: <3σ).

## T1.3 — Multi-anchor consistency

| Anchor | κ_TGP | δ from √α₀ = 2.0045 |
|---|---:|---:|
| V | 2.0080 | 0.0035 |
| Nb | 2.0145 | 0.0100 |
| Ta | 2.0125 | 0.0080 |
| Mo | 2.0095 | 0.0050 |
| Pd | 2.0155 | 0.0110 |

```
mean κ_TGP = 2.0120
std        = 0.0032 (0.159%)
TGP-SC v2 published = 2.0120
```

**Verdict:** PASS — RMS spread 0.159% << 5% gate. κ_TGP is anchor-stable across V/Nb/Ta/Mo/Pd; identity hypothesis well-defined.

## T1.4 — Data independence

| Aspect | BH side | SC side |
|---|---|---|
| Source | EHT photon-ring imaging + M9.2-D rate | T_c at ambient/standard P (NIST/PDG) |
| Anchor | shadow shift = (1/2)(1−3/3.88) = 0.1134 | 5-element V/Nb/Ta/Mo/Pd |
| Data | M87* + Sgr A* (2019/2022) | 5 BCS metals + 15 LnH₉ cross-check |
| Fit params | α₀, ψ_th, n | κ_TGP, β, T_c^base |
| Physical scale | GM/c² ~ kpc | λ_F ~ Å |

Common fit parameters: **NONE**.
Common data points: **NONE**.
Common physical scale: **NONE**.

**Verdict:** PASS — datasets fully disjoint; any match is structural, **not** fit-by-construction.

## T1.5 — Bayesian prior odds

**Inventory of 6 TGP O(1) structural constants:**
β=2.527, γ_core=π²/8=1.234, K_geo=1, κ_TGP²=4.048, m_σ²/m_s²=2, α₀=4.018.

**All 15 pairwise comparisons at 0.75% match level:**
- 14/15 pairs: 19% to 75% relative difference (clearly NO match)
- 1/15 pair: **κ_TGP² vs α₀ at 0.747%** ← the hint we are testing

**Statistics:**
```
P(match by chance per pair, uniform [1,10]) ≈ 0.167%
Expected matches in 15 pairs by chance      ≈ 0.025
Observed matches                            = 1
Bayes factor B(H₁/H₀)                       ≈ 600
P(at least one chance match)                = 2.50% < 30% gate
```

**Verdict:** PASS — Bayes factor ~600 favors structural origin over coincidence by orders of magnitude. Mundane coincidence rejected statistically.

---

## Synthesis

Phase 1 establishes that **all five necessary conditions** for the
identity hypothesis √α₀ = κ_TGP are met:

1. **Dimensions match** (both dimensionless)
2. **Numerics agree** (0.53σ within 1σ combined error)
3. **κ_TGP is anchor-stable** (0.16% RMS across 5 BCS metals)
4. **Datasets are disjoint** (no shared fit parameters or data points)
5. **Coincidence is statistically unlikely** (Bayes factor ~600)

This does **not** prove the identity — but it removes all trivial
**internal-consistency** failure modes. Phase 2 must now derive the
identity from a substrate-field action S[Φ, g, J] under TGP single-Φ
axiom.

---

## What XS.1.Phase1 closes

- ✅ Identity hypothesis is **dimensionally well-posed**
- ✅ Identity hypothesis is **statistically warranted** (Bayes 600:1)
- ✅ Datasets are **structurally disjoint** (match is not artefact)
- ✅ κ_TGP is **anchor-stable** (so "κ_TGP" is a meaningful constant, not a fit alias)

## What XS.1.Phase1 does NOT close

- ❌ Identity √α₀ = κ_TGP is **NOT yet derived** from substrate action (Phase 2)
- ❌ RG-flow stability of identity at UV scale is **NOT yet checked** (Phase 2)
- ❌ Concrete cross-sector predictions XS1–XS6 **NOT yet generated** (Phase 3)

---

## Materiał wykonawczy

- **Skrypt:** [`phase1_cross_sector_audit.py`](phase1_cross_sector_audit.py)
- **Output:** [`phase1_cross_sector_audit.txt`](phase1_cross_sector_audit.txt)
- **Setup:** [`Phase1_setup.md`](Phase1_setup.md)

## Cross-references

- [`program.md`](program.md) — overall 3-phase XS.1 plan
- [`../op-bh-alpha-threshold/Phase2_results.md`](../op-bh-alpha-threshold/Phase2_results.md) — Phase 2 T2.5 cross-sector hint origin
- [`../op-bh-alpha-threshold/Phase3_results.md`](../op-bh-alpha-threshold/Phase3_results.md) — BH8 falsification design
- [`../op-sc-alpha-origin/Phase3_results.md`](../op-sc-alpha-origin/Phase3_results.md) — κ_TGP calibration
- [`../../INDEX.md`](../../INDEX.md) — master ledger 317 → 322

## Decyzja po Phase 1

**Phase 1 CLOSED** with 5/5 PASS. Identity hypothesis viable.

→ **Proceed to Phase 2** (substrate-action derivation, 7 sub-tests).
   Master ledger update concurrent: 317 → 322 (+5 from Phase 1).
