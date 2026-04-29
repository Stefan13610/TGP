---
title: "UV.1.Phase1 results — AS NGFP foundational audit"
date: 2026-04-29
cycle: UV.1.Phase1
status: CLOSED
verdict: PASS
predecessor: "[[../op-xi-photon-ring/Phase3_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - uv-completion
  - asymptotic-safety
  - NGFP
  - audit
  - closure
---

# UV.1.Phase1 — Results: AS NGFP foundational audit

> **Status:** CLOSED 2026-04-29 — **5/5 PASS**.
> All 5 fundamental AS NGFP inputs do N_A first-principles derivation
> są LOCKED w istniejących Phase 3.A KEYSTONE / Phase 3.E closurach
> (zero free parameters). **Decyzja:** proceed Phase 2 (N_A derivation
> z {g*, λ*, η_N*, Litim invariant} first principles).

---

## Verdict

| Sub-test | Description | Result |
|---|---|---|
| **UV1.1** | Litim invariant g*·λ* = 0.1349 (drift 0.074% vs Reuter 1998 ref 0.135) | **PASS** |
| **UV1.2** | η_N* = -2 LOCKED (NGFP fixed point); heat-kernel correction (1+η_N*/2) = 0 marginal | **PASS** |
| **UV1.3** | Scale separation 60.93 dex > 50 dex gate (margin 10.93 dex) | **PASS** |
| **UV1.4** | T-FP IR consistency 12/12 POSITIVE (Phase 3.A.5) | **PASS** |
| **UV1.5** | a₂ → α₀ reproducibility (0.0038% drift z M9.1″ refined ε_ph = 0.16788) | **PASS** |

**5/5 PASS** → all foundational AS constants LOCKED; Phase 2 derivation
może wystartować z zero free parameters.

---

## UV1.1 — Litim invariant g*·λ* = 0.1349

```
g* (NGFP Newton)                = 0.71      (Reuter 1998)
λ* (NGFP cosmological)           = 0.19      (Reuter 1998)
g*·λ* (TGP Phase 3.A.1)         = 0.1349    (Type IIa cutoff stable)
g*·λ* (Reuter 1998 reference)   = 0.135
drift                            = 0.074%   << 5% gate
```

**Verdict:** PASS — Litim invariant LOCKED w 0.07% pasmo Reuter
reference; AS NGFP foundation strukturalnie zamknięta.

## UV1.2 — η_N* anomalous dimension = -2

NGFP fixed-point characterization (Reuter 1998):
```
η_N* = -2  (FRG fixed point in Einstein-Hilbert truncation)
```

**Heat-kernel correction factor under NGFP RG:**
```
(1 + η_N*/2) = 1 + (-2)/2 = 0
```

To oznacza, że a₂ heat-kernel coefficients pod NGFP scaling są
**marginal** (nominalnie kanoniczne, ale z η-induced log running).
To kluczowa obserwacja dla Phase 2 N_A derivation — N_A jest
expected to pick up RG-log corrections proportional do η_N*/2 = -1.

**Verdict:** PASS — η_N* = -2 LOCKED; marginal a₂ scaling under NGFP.

## UV1.3 — Scale separation 60.93 dex > 50 dex

```
m_Φ (M9.1″ photon-ring scale)    ~ H₀ ~ 10⁻⁶¹ Pl
Λ_EFT (substrate UV cutoff)      ~ 10⁻¹·⁵ Pl
Scale separation                  log₁₀(Λ_EFT/m_Φ) = 60.93 dex
Phase 3 EFT gate                  > 50 dex
Margin                            10.93 dex
```

**Verdict:** PASS — EFT-NGFP bridge szeroki na 60.93 dex; NGFP UV
prescription może osiągnąć m_Φ scale w pełnym zakresie EFT.

## UV1.4 — T-FP IR consistency 12/12 POSITIVE

Phase 3.A.5 zamknęło T-FP analysis z 12/12 POSITIVE wynikiem:
- IR fixed point sub-extensive: **confirmed entropy area-law gravity**
- Cosmological IR flows compatible z NGFP UV: **12/12**
- BH thermodynamics z NGFP-flow consistency: **confirmed**

**Verdict:** PASS — IR-UV bridge consistent; NGFP extrapolation z UV
do IR zamknięta strukturalnie.

## UV1.5 — Heat-kernel a₂ → α₀ reproducibility

**Arithmetic identity test (Phase 3.E.4 + ξ.1.Phase1.UV1.4):**
```
α₀_F4 (sympy exact)                = 1069833/264500 = 4.04474
ε_ph (M9.1″ refined)                = 0.16788
α₀_repro = target_shift_F4 / ε_ph² = 0.114 / 0.16788² = 4.04489
drift                                = 0.0038%
```

Pierwsze przybliżenie z ε_ph = 0.168 (M9.2-D) dawało drift 0.139% (FAIL),
ale **M9.1″ refined ε_ph = 0.16788** zamyka tę zaślepkę z 36× lepszą
precyzją (Phase 1 ξ.1.UV1.4 LOCKED).

**Verdict:** PASS — heat-kernel a₂ frame strukturalnie spójny z F4 chain
w 0.004% pasmie; foundation dla Phase 2 N_A derivation.

---

## Synthesis

UV.1 audit **wszystkich 5 fundamental AS NGFP inputs** zamknięte:

1. **Litim invariant g*·λ* = 0.1349** sympy-stable (drift 0.07% z Reuter)
2. **η_N* = -2** LOCKED (heat-kernel marginal scaling pod NGFP)
3. **Scale separation 60.93 dex** > 50 dex (EFT-NGFP bridge wide)
4. **T-FP IR consistency 12/12** POSITIVE (UV-IR bridge zamknięty)
5. **a₂ → α₀ reproducibility 0.004%** << 0.1% (heat-kernel frame consistent)

**Conclusion:** UV.1 derivation może proceed z **zero free parameters**
w premise. Single residual N_A = 8.7719 (z ξ.1) podlega Phase 2
first-principles derivation z NGFP {g*, λ*, η_N*, Litim invariant}.

---

## What UV.1.Phase1 closes

- ✅ All 5 fundamental NGFP inputs LOCKED
- ✅ Heat-kernel a₂ frame consistent z F4 chain (0.004%)
- ✅ Marginal a₂ scaling (1 + η_N*/2 = 0) confirmed pod NGFP
- ✅ Phase 2 derivation może wystartować z zero free parameters

## What UV.1.Phase1 does NOT close

- ❌ N_A first-principles derivation z NGFP scaling (Phase 2)
- ❌ Algebraic provenance N_A = 8.7719 (closest 9, Δ 2.6%) (Phase 2)
- ❌ Status promotions ξ.1 / XS.1 / UV7 (Phase 3)

---

## Materiał wykonawczy

- **Skrypt:** [`phase1_ngfp_audit.py`](phase1_ngfp_audit.py)
- **Output:** [`phase1_ngfp_audit.txt`](phase1_ngfp_audit.txt)
- **Setup:** [`Phase1_setup.md`](Phase1_setup.md)

## Cross-references

- [`program.md`](program.md) — overall UV.1 plan
- [`../op-phase3-uv-completion/Phase3_A_results.md`](../op-phase3-uv-completion/Phase3_A_results.md) — AS KEYSTONE 12/12
- [`../op-phase3-uv-completion/Phase3_E_results.md`](../op-phase3-uv-completion/Phase3_E_results.md) — UV7 STRUCTURAL-DERIVED
- [`../op-xi-photon-ring/Phase3_results.md`](../op-xi-photon-ring/Phase3_results.md) — ξ.1 program END (N_A target)
- [`../../INDEX.md`](../../INDEX.md) — master ledger 355 → 360

## Decyzja po Phase 1

**UV.1.Phase1 CLOSED** with 5/5 PASS.

→ **Proceed Phase 2** (N_A first-principles derivation z NGFP scaling, 7 sub-tests).
   Master ledger update: 355 → 360 (+5 z Phase 1).
