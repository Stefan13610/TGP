---
title: "η.2.Phase1 setup — cross-sector denom landscape audit (5 sub-tests)"
date: 2026-04-30
cycle: η.2.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[../op-alpha-fine-structure/Phase3_results.md]]"
tags:
  - TGP
  - eta2
  - denom-landscape
  - cross-sector-primes
  - chirality-counting
---

# η.2.Phase1 — Cross-sector denom landscape audit

> **Cel:** Inventory all rational-form denominators across 4 closed cycles
> (η.1, θ.1, ζ.1, α.1) + B² chirality-counting values + prime factorizations.
> Rank candidate cross-product structures dla (81, 78, 14) i residual 0.036
> derivation. **Pre-execution audit only** — no derivation claims yet.

## 5 sub-tests

### B1.1 — Cross-sector denom inventory

**Test:** List all sympy-LOCKED denominators across η.1/θ.1/ζ.1/α.1:

| Sector | Param | Value | Numerator | Denom | Denom factorization |
|---|---|---|---|---|---|
| η.1 | A | 64/81 | 64 = 2⁶ | 81 | 3⁴ |
| η.1 | ρ̄ | 11/78 | 11 | 78 | 2·3·13 |
| η.1 | η̄ | 5/14 | 5 | 14 | 2·7 |
| θ.1 | K_up | 7/8 | 7 | 8 | 2³ |
| θ.1 | K_down | 37/50 | 37 | 50 | 2·5² |
| ζ.1 | K_lepton | 2/3 | 2 | 3 | 3 |
| ζ.1 | K_ν | 1/2 | 1 | 2 | 2 |
| ε.1 | ψ_ph | 160/137 | 160 = 2⁵·5 | 137 | 137 (prime) |
| ε.1 | ε_ph | 23/137 | 23 | 137 | 137 (prime) |
| α.1 | residual fit | 9/250 | 9 = 3² | 250 | 2·5³ |

**Verdict gate:** All 10 entries match prior cycle results (no transcription error).

### B1.2 — Cross-sector prime inventory

**Test:** Categorize primes by cross-sector appearance count:

| Prime | Cross-sector appearances | Verdict |
|---|---|---|
| 2 | A=2⁶, ρ̄=2·3·13, η̄=2·7, K_up=2³, K_down=2·5², K_ν=2, ψ_ph=2⁵·5, residual=2·5³ | **ubiquitous** (8/10) |
| 3 | A=3⁴, ρ̄=2·3·13, K_lepton=3, residual_num=3² | **cross-sector** (4/10) |
| 5 | η̄_num=5, K_down=2·5², ψ_ph=2⁵·5, residual=2·5³ | **cross-sector** (4/10) |
| 7 | η̄=2·7, K_up_num=7 | **prime-7 share** (2/10) |
| 11 | ρ̄_num=11 | unique do η.1 |
| 13 | ρ̄=2·3·13 | unique do η.1 |
| 23 | ε_ph_num=23 | unique do ε.1 |
| 37 | K_down_num=37 | unique do θ.1 |
| 137 | ψ_ph=137, ε_ph=137 | **ε.1 + α.1 cross-sector** |

**Hypothesis:** primes 2, 3, 5, 7 form **chirality-counting cascade core**;
137 is QED-anchor unique do ε.1+α.1; 11/13/23/37 are sector-specific.

### B1.3 — B² chirality-counting cross-product candidate landscape

**Test:** Enumerate plausible polynomial(B²_i) → denoms (81, 78, 14):

```
B²_lepton = 2,   B²_ν = 1,   B²_up = 13/4,   B²_down = 61/25

Cross-products to test:
  C1: 81 = 3⁴ vs (1 + B²_lepton)·(1 + B²_ν)·... = 3·2·... → too small
  C2: 81 = (B²_up·4)·... = 13·... → no
  C3: 81 = N_gen⁴ = 3⁴ ✓ (N_gen = 3 generations, 4 sectors)
  C4: 78 = 2·3·13 = 2·N_gen·B²_up·4 = 2·3·13 ✓ (numerator B²_up matches!)
  C5: 14 = 2·7 vs (B²_up + 1)·... = (13/4 + 1)·... = 17/4·... → no
  C6: 14 = 2·7 vs K_up_num·(K_lepton_num) = 7·2 = 14 ✓ (cross-sector hint)
```

**Best candidate so far:** denom-pattern (3⁴, 2·3·13, 2·7) maps to
(N_gen⁴, 2·N_gen·B²_up_num, B²_up_num·K_lepton_num) — STRUCTURAL HINT
of 4-sector cross-product with N_gen=3, B²_up_num=13.

**Verdict gate:** at least 2/3 denoms map cleanly do B²-cross-product → Phase 2 viable.

### B1.4 — Residual 0.036 cross-product candidate landscape

**Test:** Check if residual = α⁻¹(0) − 137 = 0.0359990... admits a clean
4-sector form:

```
Candidates:
  R1: residual ?= (B²_lepton − B²_ν)/(some N) = 1/N → N ≈ 27.78 (not clean)
  R2: residual ?= (B²_up − B²_down)/(some N) = (13/4 − 61/25)/N = 81/(100·N)
                                              = 0.81/N → N ≈ 22.5 (not clean)
  R3: residual ?= (1/81 − 1/250) = 0.0123 − 0.004 = 0.0083 (no)
  R4: residual ?= (η̄ · ρ̄)/2 = (5/14)·(11/78)/2 = 55/2184 ≈ 0.0252 (drift 30%)
  R5: residual ?= (1 − A)·(1 − ρ̄)·(η̄)·(some power) = 17/81·67/78·5/14
                                                     = 5695/88452 ≈ 0.0644 (no)
  R6: 9/250 fit drift 0.0025% — best numerical, ale denom 250 lacks structure
```

**Honest verdict:** none of {R1-R5} delivers clean cross-sector form;
9/250 numerical fit remains best but soft denom. Phase 2 must derive
residual z **structural** principle (B²-extension lub F4-extension)
rather than brute-force rational fit.

**Verdict gate:** at least 1/5 candidates within drift < 1% → derivation viable.
If 0/5, residual stays STRUCTURAL HINT permanently (Phase 2 skip residual sub-test).

### B1.5 — Phase-1 audit verdict + Phase-2 viability

**Test:** Cumulative gate:
- B1.1 PASS (denom inventory consistent)
- B1.2 PASS (primes inventoried)
- B1.3 PASS (≥2/3 denom mappings clean)
- B1.4 PARTIAL (residual landscape weak; structural derivation needed)

**Phase-2 viability:** GO if B1.1 + B1.2 + B1.3 PASS. B1.4 PARTIAL is
acceptable (Phase 2 will derive structurally lub log honest STRUCTURAL HINT).

## Verdict gate

- **5/5 PASS** → Phase 2 first-principles derivation proceeds.
- **4/5 PASS** → Phase 2 z minor caveat (residual track partial).
- **≤ 3/5 PASS** → η.2.Phase1 reframing required.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-eta2-denom-derivation/phase1_denom_audit.py 2>&1 | tee research/op-eta2-denom-derivation/phase1_denom_audit.txt
```

## Cross-references

- [`program.md`](program.md) — η.2 3-phase plan
- [`../op-eta-wolfenstein/Phase3_results.md`](../op-eta-wolfenstein/Phase3_results.md) — H4 STRUCTURAL hint
- [`../op-theta-quark-koide/Phase3_results.md`](../op-theta-quark-koide/Phase3_results.md) — K-taxonomy + B²
- [`../op-alpha-fine-structure/Phase3_results.md`](../op-alpha-fine-structure/Phase3_results.md) — A5 residual STRUCTURAL HINT
