---
title: "UV.1.Phase1 setup — AS NGFP foundational audit"
date: 2026-04-29
cycle: UV.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
tags:
  - TGP
  - uv-completion
  - asymptotic-safety
  - NGFP
  - audit
  - falsification
---

# UV.1.Phase1 — AS NGFP foundational audit

> **Cel:** Sprawdzić, czy 5 fundamental AS NGFP constants są LOCKED w
> istniejących Phase 3.A / Phase 3.E closurach przed Phase 2 N_A
> first-principles derivation.

---

## Hipoteza Phase 1 (H_UV1)

Pod TGP Phase 3.A KEYSTONE (4-of-4 UV completion structural compatibility)
i Phase 3.E (B.4/B.6/Δ_target) wszystkie 5 fundamental AS NGFP inputs do
N_A derivation są **already LOCKED** z **zero free parameters**:

```
g*·λ* (Litim invariant) = 0.1349       (Reuter 1998 + Phase 3.A.1)
η_N* (anomalous dim)    = −2           (Reuter NGFP fixed point)
m_Φ/Λ_EFT scale sep      > 50 dex      (Phase 3.A.6: 60.93 dex)
T-FP IR consistency     = 12/12 POS    (Phase 3.A.5)
heat-kernel a₂ → α₀     = 0.0% drift   (Phase 3.E.4)
```

---

## 5 sub-testów Phase 1

### UV1.1 — Litim invariant g*·λ* = 0.1349 LOCKED

**Cel:** Verify Litim invariant z Phase 3.A.1 jest stable Reuter
reference (within 5%).

**Test:**
- Reuter 1998 reference: g*·λ* = 0.135 (Type IIa cutoff)
- Phase 3.A.1 closure: g*·λ* = 0.1349 (drift 0.07%)
- Independence: drift < 5% gate (Phase 3 R-final)

**Falsification:** drift > 5% z Reuter reference → AS NGFP foundation
broken; UV.1 cannot proceed.

### UV1.2 — η_N* anomalous dimension = −2 LOCKED

**Cel:** Verify NGFP-fixed-point η_N* = −2 (Reuter NGFP characterization).

**Test:**
- Reuter 1998: η_N* = −2 (FRG fixed point in Einstein-Hilbert truncation)
- Heat-kernel correction factor under NGFP: (1 + η_N*/2) = 0
- This means a₂ heat-kernel coefficients have peculiar UV scaling: nominally
  marginal but with η-induced log running

**Falsification:** η_N* drift > 10% z Reuter → NGFP fixed-point ambiguity.

### UV1.3 — Scale separation m_Φ/Λ_EFT = 60.93 dex > 50 dex

**Cel:** Confirm EFT-NGFP bridge is wide enough to support EFT closure.

**Test:**
- M9.1″: m_Φ ~ H₀ ~ 10⁻⁶¹ Pl
- Λ_EFT ~ 10⁻¹·⁵ Pl (UV cutoff for substrate EFT)
- Scale separation: log₁₀(Λ_EFT/m_Φ) = 60.93 dex
- 50 dex gate: Phase 3 EFT compatibility

**Falsification:** dex < 50 → EFT validity bridge collapsed; NGFP UV
prescription cannot reach m_Φ scale.

### UV1.4 — T-FP IR consistency 12/12 POSITIVE

**Cel:** Verify temperature-FP IR fixed-point analysis is fully consistent
across 12 thermodynamic checks (Phase 3.A.5).

**Test:**
- Phase 3.A.5: T-FP analysis 12/12 POSITIVE (Phase 3 R-final 60/60)
- IR fixed point sub-extensive: confirmed entropy area-law gravity behavior
- Cosmological IR flows compatible z NGFP UV: 12/12 sub-tests

**Falsification:** < 12/12 POSITIVE → IR-UV bridge inconsistent; NGFP
extrapolation z UV do IR broken.

### UV1.5 — Heat-kernel a₂ → α₀ reproducibility 0.0% drift

**Cel:** Verify a₂ heat-kernel framework correctly reproduces F4 rational
α₀ chain (Phase 3.E.4).

**Test:**
- α₀_F4 = 1069833/264500 = 4.04472 (sympy exact)
- α₀_repro from a₂ frame: 4 + correction_factor·Δ_target
- Phase 3.E.4: drift 0.0000% < 5% gate
- This confirms a₂ heat-kernel frame is consistent w F4 chain

**Falsification:** drift > 0.1% → a₂ frame broken; UV.1 Phase 2 derivation
cannot use heat-kernel approach.

---

## Verdict gate

**5/5 PASS** → all 5 fundamental NGFP inputs LOCKED, proceed Phase 2
N_A derivation z zero free parameters w premise.

**≤ 4/5 PASS** → audit gap, Phase 2 deferred until gap zaadresowany.

---

## Materiał wykonawczy

- **Skrypt:** [`phase1_ngfp_audit.py`](phase1_ngfp_audit.py) (5 sub-tests, sympy + numeric)
- **Output:** `phase1_ngfp_audit.txt`
- **Memo:** `Phase1_results.md`

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-uv-as-ngfp/phase1_ngfp_audit.py 2>&1 | tee research/op-uv-as-ngfp/phase1_ngfp_audit.txt
```

---

## Constants used

```
g* (NGFP Newton)        = 0.71              (Reuter 1998)
λ* (NGFP cosmological)  = 0.19              (Reuter 1998)
g*·λ* (Litim invariant) = 0.1349            (Reuter ref 0.135)
η_N* (anomalous dim)    = -2                (Reuter NGFP)
m_Φ/Λ_EFT               = 60.93 dex         (Phase 3.A.6)
T-FP                    = 12/12 POSITIVE    (Phase 3.A.5)
α₀_F4                    = 1069833/264500    (Phase 2.B.3)
N_A target              = 500/57 = 8.7719   (ξ.1.Phase2)
```

---

## Cross-references

- [`program.md`](program.md) — overall UV.1 plan
- [`../op-phase3-uv-completion/Phase3_A_results.md`](../op-phase3-uv-completion/Phase3_A_results.md) — AS KEYSTONE
- [`../op-phase3-uv-completion/Phase3_E_results.md`](../op-phase3-uv-completion/Phase3_E_results.md) — UV7 STRUCTURAL-DERIVED
- [`../op-xi-photon-ring/Phase3_results.md`](../op-xi-photon-ring/Phase3_results.md) — ξ.1 program END (N_A target)
- [`../../INDEX.md`](../../INDEX.md) — master ledger 355
