---
title: "ξ.1.Phase1 results — EFT photon-ring frame audit"
date: 2026-04-29
cycle: ξ.1.Phase1
status: CLOSED
verdict: PASS
predecessor: "[[../op-cross-sector-charge/Phase3_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - xi-factor
  - photon-ring
  - audit
  - closure
---

# ξ.1.Phase1 — Results: EFT photon-ring frame audit

> **Status:** CLOSED 2026-04-29 — **5/5 PASS**.
> All 5 fundamental inputs do a₂ heat-kernel derivation są LOCKED w
> istniejących closurach (zero free parameters). **Decyzja:** proceed
> Phase 2 (a₂ derivation).

---

## Verdict

| Sub-test | Description | Result |
|---|---|---|
| **ξ1.1** | ξ_geom = 1 LOCKED (M9.1″ vacuum / Phase 2.B.2) | **PASS** |
| **ξ1.2** | α(α−1) = 2 LOCKED (Phase 1.A.1 Theorem alpha2; α=2 unique sympy solution) | **PASS** |
| **ξ1.3** | ψ_ph − 1 = 0.168 LOCKED (M9.2-D vs M9.1″ drift 0.010%) | **PASS** |
| **ξ1.4** | F4 rational 1069833/264500 reconstructed from arithmetic identity (within 0.004%) | **PASS** |
| **ξ1.5** | Phase 2 strict (1/2)(1 − 3/3.88) sympy-exact = 11/97 = 0.11340206 | **PASS** |

**5/5 PASS** → premise inputs strukturalnie zamknięte; Phase 2
derivation może wystartować z zero free parameters.

---

## ξ1.1 — ξ_geom = 1 derivation

```
ξ_geom = α(α-1)/2 pod F2 (α=2) + F3 (φ_eq=1) + vacuum
       = 2·1/2
       = 1   (sympy exact integer)
```

**Verdict:** PASS — ξ_geom LOCKED na 1 sympy-exact, derived z F2 + F3.

## ξ1.2 — α(α−1) = 2 derivation

Sympy solve α(α−1) − 2 = 0 yields **unique positive integer solution α = 2**.

| α | α(α−1) | Constraint violated |
|---|---:|---|
| 1 | 0 | C1 stability (V'' degenerate) |
| 2 | 2 | none ✓ (Theorem alpha2) |
| 3 | 6 | C2 positivity (V'' overconfining) |
| 4+ | ≥ 12 | C1 + Phase 1 root not exists |

**Verdict:** PASS — α=2 unique sympy positive-integer solution; F2 LOCKED.

## ξ1.3 — ψ_ph − 1 = 0.168

```
ψ_ph (M9.2-D)        = 1.168
ψ_ph (M9.1″ refined) = 1.16788
|drift|              = 0.0103%
drift gate           = 0.5%
ε_ph² (sympy)        = 168²/1000² = 441/15625 = 0.028224
```

**Verdict:** PASS — M9.2-D drift 0.0103% << 0.5% gate; ε_ph² sympy-exact.

## ξ1.4 — F4 rational 1069833/264500 provenance

**Arithmetic identity (Phase 2.B.3):**
```
α₀_F4 = target_shift_F4 / [(ψ_ph_M9_1″ − 1)² · ξ_geom]
      = 0.114 / (0.16788² · 1.0)
      = 4.04489
```

Compare with F4 sympy LOCKED rational:
```
F4 sympy = 1069833/264500 = 4.04474
|reconstructed − sympy| / sympy = 0.0038%
```

Sub-percent agreement. F4 rational LOCKED w `closure_2026-04-26`.

**Verdict:** PASS — F4 rational z arithmetic identity (Phase 2.B.3 chain).

## ξ1.5 — Phase 2 strict provenance

```
target_shift_strict = (1/2)(1 − r_ph^GR/r_ph^TGP)
                    = (1/2)(1 − 3/3.88)
                    = 11/97 (sympy exact)
                    = 0.11340206
α₀_strict = target_shift_strict / ε_ph²
          = 0.11340206 / 0.028224
          = 4.01793
```

**Split z F4 frame:**
```
0.114 (F4 literal) − 0.1134 (strict) = 0.000598
relative split = 0.527%
```

Strict form jest **direct geometric photon-ring shift** (bare form, no a₂
correction). F4 form **proposed** jako post-a₂-corrected — Phase 2 do test.

**Verdict:** PASS — strict form sympy-exact 11/97 (BH.1.Phase2 LOCKED).

---

## Synthesis

Phase 1 audit **wszystkich 5 fundamental inputs** zamknięte:

1. **ξ_geom = 1** sympy-exact, z F2 + F3 + vacuum
2. **α(α−1) = 2** sympy-exact, z F2 (α=2 unique pod (C1)–(C3))
3. **ψ_ph − 1 = 0.168** stable across M9.2-D + M9.1″ (drift 0.010%)
4. **F4 rational 1069833/264500** reconstructed z arithmetic identity (0.004% drift)
5. **Phase 2 strict 11/97 = 0.1134** sympy-exact z direct geometric formula

**Conclusion:** ξ.1 derivation może proceed z **zero free parameters**
w premise. The 0.527% split between F4 (0.114) i strict (0.1134) form
jest **the** ξ-factor podlegające Phase 2 derivation.

---

## What ξ.1.Phase1 closes

- ✅ Premise audit: all 5 inputs LOCKED w istniejących closurach
- ✅ Frame distinction explicit: F4 (post-a₂-corrected proposed) vs strict (bare-form)
- ✅ Split quantified: 0.527% (= F4 literal 0.114 − strict 11/97)

## What ξ.1.Phase1 does NOT close

- ❌ a₂ heat-kernel derivation (Phase 2)
- ❌ Frame A vs Frame B selection (Phase 2)
- ❌ ξ.1 classification (Phase 2)

---

## Materiał wykonawczy

- **Skrypt:** [`phase1_frame_audit.py`](phase1_frame_audit.py)
- **Output:** [`phase1_frame_audit.txt`](phase1_frame_audit.txt)
- **Setup:** [`Phase1_setup.md`](Phase1_setup.md)

## Cross-references

- [`program.md`](program.md) — overall ξ.1 plan
- [`../op-cross-sector-charge/Phase2_results.md`](../op-cross-sector-charge/Phase2_results.md) — XS.1 ξ unresolution
- [`../op-bh-alpha-threshold/Phase2_results.md`](../op-bh-alpha-threshold/Phase2_results.md) — Phase 2 strict provenance
- [`../op-phase2-quantum-gravity/Phase2_B_results.md`](../op-phase2-quantum-gravity/Phase2_B_results.md) — F4 rational chain
- [`../../INDEX.md`](../../INDEX.md) — master ledger 336 → 341

## Decyzja po Phase 1

**ξ.1.Phase1 CLOSED** with 5/5 PASS.

→ **Proceed Phase 2** (a₂ heat-kernel first-principles derivation, 7 sub-tests).
   Master ledger update: 336 → 341 (+5 z Phase 1).
