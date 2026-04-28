---
title: "XS.1.Phase3 results — multi-sector falsification map of √α₀ = κ_TGP"
date: 2026-04-28
cycle: XS.1.Phase3
status: CLOSED
verdict: PASS
predecessor: "[[Phase2_results.md]]"
parent: "[[program.md]]"
program_status: "XS.1 program END"
tags:
  - TGP
  - cross-sector
  - alpha-0
  - kappa-TGP
  - falsification
  - predictions
  - registry
  - closure
  - program-end
---

# XS.1.Phase3 — Results: multi-sector falsification map of √α₀ = κ_TGP

> **Status:** CLOSED 2026-04-28 — **7/7 PASS**.
> **XS.1 program END.** Identity √α₀ = κ_TGP **PARTIALLY DERIVED**
> (Phase 2) i teraz **falsyfikowalna w 6 niezależnych kanałach** w
> okresie 2027–2035 z 6 nowymi pre-registered predictions XS1–XS6
> w PREDICTIONS_REGISTRY.

---

## Verdict

| Sub-test | Description | Result |
|---|---|---|
| **T3.1** | XS1: ngEHT-α₀ × SC-κ_TGP² combined precision (2030+) | **PASS** |
| **T3.2** | XS2: cross-channel orthogonality vs g̃ ≈ 0.9803 (F5) | **PASS** |
| **T3.3** | XS3: lepton sector orthogonality (Koide / r_21 / r_31) | **PASS** |
| **T3.4** | XS4: QM sector orthogonality (Born n=2, CHSH 2√2) | **PASS** |
| **T3.5** | XS5: F4 rational anchor sub-percent consistency (0.084%) | **PASS** |
| **T3.6** | XS6: combined falsification roadmap (6 channels, 2027–2035) | **PASS** |
| **T3.7** | PREDICTIONS_REGISTRY entries XS1–XS6 (6 new) | **PASS** |

**Cumulative XS.1:** 5 (Phase 1) + 7 (Phase 2) + 7 (Phase 3) = **19 sub-tests**.

---

## T3.1 — XS1: ngEHT-α₀ × SC-κ_TGP² combined precision

**Setup:** Identity √α₀ = κ_TGP testowane przez **dwa niezależne
eksperymenty**:

| Side | Experiment | Anchor | Precision |
|---|---|---|---|
| BH (α₀) | ngEHT photon-ring imaging | Δb_crit shadow shift | ~5% on α₀ (2030+) |
| SC (κ_TGP²) | TGP-SC v2 + LnH₉ ensemble | Multi-Ln calibration | ~0.5% (current) → ~0.3% (post-LnH₉ 2030) |

**Combined precision (quadrature):**
```
σ_α₀(frac)             = 0.05 (5%)
σ_κ_TGP²(frac, current) = 2 × 0.005 = 0.010 (1%)
σ_κ_TGP²(frac, post)   = 2 × 0.003 = 0.006 (0.6%)
σ_combined(current)    = sqrt(0.05² + 0.010²)  ≈ 5.10%
σ_combined(post-LnH₉)  = sqrt(0.05² + 0.006²) ≈ 5.04%
```

**Current match:** |α₀ − κ_TGP²| / κ_TGP² = **0.747%** (Phase 2 strict);
**6.7× tighter** niż 5% falsification trigger.

**XS1 prediction:** by 2030+, identity testable at ~5% combined; identity
holds (current match 0.75% well below trigger).

**Verdict:** PASS — current match << trigger; ngEHT 2030+ resolves identity.

## T3.2 — XS2: cross-channel orthogonality vs g̃

**Setup:** F5 (g̃ ≈ 0.9803, drift 0.0306%) jest niezależnym O(1)
substrate identity z Phase 2 quantum-gravity entropy. Czy g̃ wpływa na
identity √α₀ = κ_TGP?

**Analysis:**
- κ_TGP² (bare) = 4.0481
- κ_TGP² · g̃² = 3.8902 (substrate-suppressed effective, Δ = 3.9%)
- Identity √α₀ = κ_TGP jest dla **bare** Wilson coefficients
- g̃ wchodzi do T-Λ entropy (vacuum prefactor) — niezależny sektor

**XS2 prediction:** identity √α₀ = κ_TGP i g̃ ≈ 0.9803 są
strukturalnie ortogonalne; deviacja g̃ poza [0.95, 1.05] nie falsyfikuje
XS1 (są niezależnymi O(1) substrate identities).

**Verdict:** PASS — g̃ within EFT window [0.95, 1.05]; identity is
bare-coefficient (independent of g̃ correction).

## T3.3 — XS3: lepton sector orthogonality

**Setup:** Czy lepton stałe (r_21=206.77, r_31=3477, K_koide=2/3)
zawierają κ_TGP/√α₀ jako factor?

**Tests (κ_TGP-factorization search):**

| Constant | / κ_TGP² | / κ_TGP³ | / κ_TGP⁴ | Integer match? |
|---|---:|---:|---:|---|
| r_21 = 206.77 | 51.08 | 25.39 | 12.62 | NO |
| r_31 = 3477 | 858.91 | 426.92 | 212.18 | NO |

Search range: N · κ_TGP^k for N ∈ {1,..,4}, k ∈ {1,..,4}; integer
match within 1%. **Result: NONE.**

**XS3 prediction:** lepton sector (r_21, r_31, K_koide) jest
**strukturalnie ortogonalny** do κ_TGP/α₀ identity. Żadna lepton stała
nie zawiera √α₀ ani κ_TGP jako factor.

**Verdict:** PASS — lepton-cross-sector orthogonality confirmed.

## T3.4 — XS4: QM sector orthogonality

**Setup:** QM Born n=2 (QM1) i α(ψ) n=2 (Phase 2 strict) — both n=2;
czy to cross-sector identity?

**Analysis:**
- Born n=2: integer (substrate-derived emergent 2)
- α(ψ) n=2: integer (Phase 2 strict, smoothness exponent)
- CHSH bound 2√2 ≈ 2.828; CHSH²/α₀ ≈ 1.99 (not κ_TGP-related)
- Both n=2 są **independent integer selections**, nie cross-sector identity

**XS4 prediction:** Born n=2 i α(ψ) n=2 są niezależnymi structural
ekspozycjami integer n=2 (substrate selecting integer 2 w różnych sektorach);
nie cross-sector √α₀ = κ_TGP identity.

**Verdict:** PASS — QM sector orthogonal; identity exclusive to BH/SC pair.

## T3.5 — XS5: F4 rational anchor sub-percent consistency

**Setup:** F4 LOCKED w registry: α₀ = 1069833/264500 ≈ 4.04472.
Compare with κ_TGP² = 4.0481.

```
α₀_F4         = 4.044737 (sympy 1069833/264500)
κ_TGP²        = 4.048144 (V/Nb/Ta/Mo/Pd RMS)
|Δ| / κ_TGP² = 0.0842%
              = 0.000842
```

**Phase 2 strict precision:** 0.747%
**F4 form precision:** 0.084%
**F4 vs strict precision factor:** 8.86× tighter w F4 frame.

F4 form jest **sub-percent (within 0.1%)** match — strukturalnie
najtighten potwierdza identity.

**XS5 prediction:** F4 sympy rational dla α₀ jest konsystentne
z κ_TGP² do 0.084%; precision potwierdza identity w F4 frame
(LOCKED-derivative status).

**Verdict:** PASS — sub-percent match w F4 frame.

## T3.6 — XS6: combined cross-sector falsification roadmap

**Decision matrix (6 independent channels):**

| Channel | Horizon | Falsification trigger |
|---|---|---|
| ngEHT photon-ring α₀ | 2030+ | \|α₀ − κ_TGP²\|/κ_TGP² > 5% |
| LnH₉ DAC κ_TGP | 2027–2030 | TGP-SC v2 multi-Ln drift > 1% |
| MICROSCOPE-2 η_TGP | 2027–2028 | η > 1·10⁻¹⁷ → α₀ ≠ 4.02 → identity broken |
| LIGO O5 QNM δf/f | 2027+ | δf null at 0.5% → α(ψ) ≠ 0.16 → α₀ shift |
| LISA SMBH ringdown | ~2035 | combined w/ ngEHT for tighter α₀ |
| Lepton g₀^τ precision | continuous | drift > 0.1% in r_21/r_31 (orthogonality check) |

**Falsification rule:** identity √α₀ = κ_TGP rejected if **≥ 2 niezależne
kanały** detect deviation > 5% (combined Bayes update).

**XS6 prediction:** identity będzie testowane przez **6 niezależnych
channels** w okresie 2027–2035; combined roadmap robust against
single-channel false positives.

**Verdict:** PASS — 6-channel falsification roadmap established.

## T3.7 — PREDICTIONS_REGISTRY entries XS1–XS6

**6 new pre-registered predictions:**

| Code | Sektor | Status | Anchor | Horizon |
|---|---|---|---|---|
| **XS1** | 8 (cross-link 3) | LIVE | \|α₀−κ_TGP²\|/κ_TGP² ≤ 3-5% by 2030+ | 2030+ ngEHT × SC v2 |
| **XS2** | 8 | STRUCTURAL | g̃ ≈ 0.9803 i √α₀=κ_TGP są niezależne | continuous |
| **XS3** | 8 (cross-link 6) | STRUCTURAL | lepton orthogonal do κ_TGP/α₀ | continuous |
| **XS4** | 8 (cross-link 9) | STRUCTURAL | Born n=2 niezależne od α(ψ) n=2 | continuous |
| **XS5** | 8 | LOCKED-derivative | F4 rational ≡ κ_TGP² do 0.084% | sympy-exact |
| **XS6** | 8 (combined) | LIVE | 6-channel falsification roadmap | 2027–2035 |

**Verdict:** PASS — 6/6 predictions ready for registry.

---

## Synthesis

**XS.1 program END** w 19 sub-testach (5+7+7):

1. **Phase 1 audit** (5/5): identity dimensionalna, statystycznie
   warranted (Bayes 600), datasets disjoint, anchor-stable.
2. **Phase 2 derivation** (7/7): substrate-action S[Φ, g, J, ψ_e]
   wymusza single g_TGP (F1 + F3); M_BH = M_SC = M_universal pod
   canonical normalization (K_geo=1, φ_eq=1); RG-invariant; **PARTIALLY
   DERIVED** (ξ unresolved between F4 rational i strict form).
3. **Phase 3 falsification** (7/7): identity testowalna w **6
   niezależnych channels** 2027–2035; XS1–XS6 registered.

**Cross-sector consequences:**
- BH photon-ring T·J·J coupling (α₀) i SC pair-breaking (κ_TGP²) są
  **same Wilson coefficient** of substrate-action under canonical
  normalization.
- Identity provides **falsification leverage**: jakakolwiek deviacja
  na ngEHT od α₀ = 4.02 ± 5% albo na LnH₉ od κ_TGP = 2.012 ± 1%
  zagraża identity AND obu Phase 2/3 cycles (BH.1 + SC.1 współspoinają
  się w XS.1).

---

## What XS.1 closes

- ✅ Identity √α₀ = κ_TGP **PARTIALLY DERIVED** (Phase 2)
- ✅ **F4 rational form** anchors identity at 0.084% precision (LOCKED)
- ✅ **6 niezależnych falsification channels** registered (XS1–XS6)
- ✅ **Cross-sector orthogonality maps:** lepton (XS3), QM (XS4), g̃ (XS2)
  are **separate** structural identities (not coincidence accumulation)
- ✅ Identity provides **combined falsification leverage** across BH+SC sectors

## What XS.1 does NOT close (long-term track)

- ❌ Full Phase 1 PLAN action-principles derivation of ξ-factor (target_shift = 0.114 vs 0.1134)
- ❌ UV completion route selection (which UV preserves identity)
- ❌ Lepton-sector deeper embedding (does √α₀ appear in g₀^τ first-principles?)
- ❌ Phase 0-grade derivation of identity z minimal axiomatic kernel

---

## Materiał wykonawczy

- **Skrypt:** [`phase3_falsification_map.py`](phase3_falsification_map.py)
- **Output:** [`phase3_falsification_map.txt`](phase3_falsification_map.txt)
- **Setup:** [`Phase3_setup.md`](Phase3_setup.md)

## Cross-references

- [`program.md`](program.md) — overall 3-phase XS.1 plan
- [`Phase1_results.md`](Phase1_results.md) — feasibility 5/5 PASS
- [`Phase2_results.md`](Phase2_results.md) — substrate-action 7/7 PASS, PARTIALLY DERIVED
- [`../op-bh-alpha-threshold/Phase3_results.md`](../op-bh-alpha-threshold/Phase3_results.md) — BH4–BH9 baseline
- [`../op-sc-alpha-origin/Phase3_results.md`](../op-sc-alpha-origin/Phase3_results.md) — SC4–SC7 baseline
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — XS1–XS6 entries (post-update)
- [`../../INDEX.md`](../../INDEX.md) — master ledger 329 → 336

## Decyzja po Phase 3

**XS.1.Phase3 CLOSED** with 7/7 PASS. **XS.1 program END.**

→ 6 new XS1–XS6 predictions LIVE/STRUCTURAL/LOCKED w registry.
→ Identity √α₀ = κ_TGP **PARTIALLY DERIVED** with **6-channel falsification** roadmap 2027–2035.
→ Master ledger update: 329 → 336 (+7 z Phase 3).
→ Cumulative XS.1: **19 sub-tests** (5+7+7).
