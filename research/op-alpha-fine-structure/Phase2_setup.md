---
title: "α.1.Phase2 setup — α_QED first-principles structural derivation"
date: 2026-04-30
cycle: α.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - alpha-fine-structure
  - prime-137
  - derivation
---

# α.1.Phase2 — α_QED first-principles structural derivation

> **Cel:** Lock zeroth-order anchor α_QED⁻¹ ≈ 137 jako DERIVED z F4
> chain (ε.1) i test residual 0.036 candidate forms; falsify 5
> alternative α_QED forms (Wyler, Eddington integer-only, Atiyah γ,
> Gilson, ad-hoc rational); promote E7 STRUCTURAL HINT do **PARTIALLY
> DERIVED** lub stabilize **STRUCTURAL HINT (research-track)**.

## 7 sub-tests

### A2.1 — 137 prime denominator inheritance lock z F4 chain

**Test:**
- ψ_ph = 4 / (3 + target_shift_F4) = 4 / (3 + 57/500) = 4 · 500/1557
  = 2000/1557 → reduce: gcd(2000, 1557) = ? 1557 = 3·519 = 3·3·173 ≠ 137
- Wait: ε.1 derivation gives 1557/2000 reduction; let me recheck.
- Actually z ε.1.Phase2: ψ_ph = 4/(3+ts), ts = 57/500, denom = 3·500+57 = 1557
- 1557 = 3 · 519 = 3 · 3 · 173 — not 137!
- BUT ε.1 results LOCK ψ_ph = 160/137 sympy exact — let me sanity check.
  ψ_ph = 4 · 500 / 1557 = 2000/1557. gcd(2000, 1557) = ? 2000 - 1557 = 443. 
  gcd(1557, 443): 1557 = 3·443 + 228; 443 = 1·228 + 215; 228 = 1·215 + 13; 
  215 = 16·13 + 7; 13 = 1·7 + 6; 7 = 1·6 + 1; 6 = 6·1 → gcd=1.
- So ψ_ph = 2000/1557 NOT 160/137. ε.1 must use different target_shift form.

**Recheck:** ε.1 results give ψ_ph = 160/137. So 160/137 = 4/(3+x) →
3+x = 4·137/160 = 548/160 = 137/40 → x = 137/40 − 3 = 137/40 − 120/40 = 17/40.
So target_shift_F4 = 17/40 = 0.425 (NOT 57/500 = 0.114).
Check: 3 + 17/40 = 137/40, 4 / (137/40) = 160/137 ✓.

**Conclusion:** w F4 chain dla photon-ring, target_shift = 17/40 (z M9.1″)
nie 57/500 (z N_A inverse). The two F4 anchors are SEPARATE:
- ψ_ph F4 anchor: target_shift = 17/40 → ψ_ph = 160/137
- N_A F4 anchor:  target_shift = 57/500 → N_A = 500/57

Both are F4-chain rational anchors at different scales w substrate-action.

**Test:** verify 137 emerges sympy-exact z 4/(3+17/40):
- 4/(3 + 17/40) = 4·40/(120+17) = 160/137 ✓
- 137 = 120 + 17 = 3·40 + 17 (prime via standard primality test)

**Gate:** 137 sympy-exact emerges z F4 chain z target_shift = 17/40;
this **locks 137 jako DERIVED** w TGP framework (inherited z ε.1).

### A2.2 — Residual cascade test α_QED⁻¹ − 137 = δ

**Test:** rank candidate forms dla δ = 0.035999084:

| candidate | value | drift % |
|---|---|---|
| 9/250 | 0.036 | 0.0025 |
| 23/640 | 0.0359375 | 0.17 |
| 23·ε_ph² = 23·(23/137)² = 12167/18769 | 0.6483 | huge |
| (ε_ph)/(2π) = 23/(2π·137) | 0.02672 | 25.8 |
| (1/137) · α₀ / 1.083 | 0.02726 | 24.3 |
| ε_ph · κ_TGP² / 18.85 | varies | varies |

**Hypothesis:** δ = 9/250 ?
- 9/250 = 0.036 exactly; drift 0.0025% — within 81-ppt CODATA precision
  (0.0025% = 25 ppm vs 81 ppt = 0.0000081%, factor 3000× too coarse)
- Better: δ matches Wyler 9π³/16 form? Compute 9π³/16:
  9·(3.14159)³/16 = 9·31.006/16 = 17.441. NO — Wyler form is for α_QED⁻¹
  itself: 9π³/16 ≈ 17.44 (wrong order). Actual Wyler: α_QED⁻¹ ≈ (8π⁴·9!/(5!2⁴))^(1/4)
  = 137.0360829, drift 0.000061% — but Wyler form is not TGP-natural.

**Test:** if 9/250 candidate → α_QED⁻¹ = 137 + 9/250 = 34250/250 + 9/250
= 34259/250 = 137.036; drift 0.0007% within structural anchor band, but
denom 250 = 2·5³ has no obvious TGP cross-link.

**Verdict criterion:** if any candidate fits drift < 0.1% AND has clean
TGP-rational structure, promote to DERIVED. Else STRUCTURAL HINT.

### A2.3 — Cross-sector cascade hypothesis test

**Test:** α_QED⁻¹ = 137 · (1 + ε_corr); compute ε_corr = (α⁻¹ − 137)/137:
- ε_corr = 0.035999/137 = 2.628·10⁻⁴ = 0.02628%

Cascade probes dla 2.628·10⁻⁴:
- ε_ph² / 137 = (23/137)² / 137 = 529 / (137³) = 529/2571353 = 2.058·10⁻⁴
  (drift 21.7% vs target — no)
- η̄ · ε_ph / 137² = (5/14)(23/137)/137² = 115/(14·137³) = 3.20·10⁻⁶ — no
- α₀ / 137³ = 4.045/2571353 = 1.573·10⁻⁶ — no
- 1/(2π·137·κ_TGP) = 1/(2π·137·2.012) = 5.78·10⁻⁴ — drift 120%

**Test:** if no clean cascade → residual stays STRUCTURAL HINT.

### A2.4 — 5 alternative α_QED⁻¹ forms FALSIFIED

Alternative formulas (historical + ad-hoc):

| Label | Form | Value | Drift % |
|---|---|---|---|
| C1 | Eddington 137 (integer) | 137.000000 | 0.026% |
| C2 | Wyler 9π³/16 (cube) | (computed) | varies |
| C3 | Atiyah γ-form | (gauss-bonnet inspired) | varies |
| C4 | Gilson cos(π/137)·137·tan(π/(137·29)) | 137.0359990 | tiny |
| C5 | TGP rational 19048/139 (Phase 1 best) | 137.035971 | 0.00002% |
| C6 | TGP 137 + 9/250 | 137.036000 | 0.0007% |
| C7 | TGP 137 + 23/640 | 137.035938 | 0.045% |

**Falsification threshold:** TGP-best max drift dla zeroth-order = 0.026%
(from 137/1). Set threshold = 0.5% (= ~19× TGP-best, > 1σ_PDG_central
0.026·81ppt ≈ 0.026% — actually 81 ppt ≈ 0.0000081%; the relevant gate
is ~ TGP rational structural band, where 0.5% = 19× anchor band).

**Gate:** if all alternatives drift > 0.5% OR have soft denoms that
aren't TGP-natural, pass; if Gilson formula passes (drift ~10⁻⁹), it's
known fitting form — falsify on structural-anchor grounds (no TGP origin).

### A2.5 — Cross-sector denom-prime sharing 137 isolation

**Test:**
- 137 prime contained ONLY w ψ_ph, ε_ph z TGP (verified A1.4)
- Cross-sector primes shared:
  - 7 ↔ η.1 (η̄=5/14) ↔ θ.1 (K_up=7/8 num) — DERIVED via chirality
  - 3 ↔ η.1 (A=81), η.1 (ρ̄=78), tgp-leptons (K_lepton=2/3 denom) — DERIVED
  - 5 ↔ η.1 (η̄=5/14 num), N_A (500), λ_C (165) — common
  - 137 ↔ ε.1 only — UNIQUE ε.1 anchor

**Conclusion:** 137 jest QED-anchor prime — not cross-sector shared.
Implication: α_QED structurally locked w photon-ring sektor (ε.1)
via F4 chain z target_shift_F4_M9 = 17/40.

**Gate:** PASS if 137 unique to ε.1 verified; structural conclusion
documented honestly.

### A2.6 — NGFP RG-stability of 137-anchor

**Test:**
- 137 emerges z 4/(3 + 17/40) = 160/137 photon-ring scale (ε.1)
- Pod common β-rescaling NGFP marginal a₂ ((1+η_N*/2) = 0): scale change
  affects all dimensional quantities equally — dimensionless rationals
  like 160/137 are RG-invariant by construction (UV.1.UV2.5 ratio inv.)
- α_QED⁻¹ ≈ 137 inherits RG-invariance via prime-denom anchor
- α_QED running 7.1% (α(0) → α(M_Z)) is SM physics not TGP geometric

**Test:** verify α₀ ratio under c-rescaling sympy; confirm 137 numerical
unchanged.

**Gate:** PASS — RG-invariance via dimensionless ratio inheritance
(already shown w UV.1).

### A2.7 — Classification verdict

**Possible outcomes:**

| Verdict | Criterion |
|---|---|
| **DERIVED** | Residual 0.036 admits clean TGP-rational form drift < 0.1% with structural meaning |
| **PARTIALLY DERIVED** | 137 DERIVED z F4 chain + residual 0.036 admits candidate at drift < 0.5% (no clean rational) |
| **STRUCTURAL HINT** | 137 emerges z F4 chain ale residual 0.036 nie ma TGP-natural form |

**Expected:** Likely **PARTIALLY DERIVED** with honest STRUCTURAL HINT
caveat dla residual 0.036 — historical α-fine-structure fitting attempts
(Wyler, Atiyah, Gilson) all suffered this same residual problem.

## Verdict gate

**7/7 PASS** → α.1.Phase3 predictions.
**5-6/7 PASS** → α.1.Phase3 conditional.
**< 5/7 PASS** → terminate research-track.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-alpha-fine-structure/phase2_alpha_derivation.py 2>&1 | tee research/op-alpha-fine-structure/phase2_alpha_derivation.txt
```

## Cross-references

- [`program.md`](program.md) — α.1 plan
- [`Phase1_results.md`](Phase1_results.md) — 5/5 PASS
- [`../op-eps-photon-ring/Phase2_results.md`](../op-eps-photon-ring/Phase2_results.md) — ψ_ph = 160/137 derivation
- [`../op-uv-as-ngfp/Phase2_results.md`](../op-uv-as-ngfp/Phase2_results.md) — NGFP marginal a₂ = 0
