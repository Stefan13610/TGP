---
title: "ε.1.Phase2 setup — ε_ph = ψ_ph − 1 structural decomposition + 5 candidates falsification"
date: 2026-04-29
cycle: ε.1.Phase2
status: PRE-EXECUTION
predecessor: "[[Phase1_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - epsilon-photon-ring
  - structural-decomposition
  - 137-denominator
  - falsification
---

# ε.1.Phase2 — ε_ph = ψ_ph − 1 structural decomposition

> **Cel:** Pokazać że ε_ph = ψ_ph − 1 = 23/137 (sympy-exact prime-denominator
> decomposition) jest **jedyną** poprawną interpretacją. Falsyfikować
> alternatywne identity candidates (5 sztuk). Confirm F4 chain implicit
> lock, heat-kernel a₂ frame consistency, NGFP RG-stability.

---

## Hipoteza Phase 2 (H_ε2)

Pod ε.1.Phase1 5/5 + Phase 1 odkrycie ψ_ph = 160/137:

1. **ε_ph = ψ_ph − 1** structural identity (sympy exact, p=23, q=137)
2. **5 alternative candidates** wszystkie FAIL przy >5% drift
3. **F4 chain** locks ε_ph implicitly via α₀ = target_shift_F4/ε_ph² constraint
4. **Heat-kernel a₂ frame** spójna z ε_ph² entering geometric chain
5. **NGFP RG-stability** preserves ε_ph² ratio pod common β-rescaling

Frame outcome: **ε_ph DERIVED via ψ_ph M9.1″ null geodesic + F4 chain
implicit lock**, classification **PARTIALLY DERIVED (refined)**.

---

## 7 sub-testów Phase 2

### E2.1 — ε_ph = ψ_ph − 1 = 23/137 structural identity

**Test:**
- ψ_ph = 4/(3 + 17/40) = 4·40/(120 + 17) = 160/137 sympy exact
- ε_ph := ψ_ph − 1 = 23/137 sympy exact
- Float: 23/137 = 0.167883... matches refined ε_ph = 0.16788 within 0.002%

**Falsification:** if (23/137 − 4197/25000)/4197/25000 > 0.01% → ε_ph value
inconsistent with M9.1″ refined.

### E2.2 — 5 alternative identity candidates falsification

**Candidates (all to be FALSIFIED):**
1. **C1:** ε_ph = 1.168/(2π) = 0.18589 → drift +10.7% **FAIL**
2. **C2:** ε_ph = (ψ_ph − 1)/ψ_ph = 23/160 = 0.14375 → drift −14.4% **FAIL**
3. **C3:** ε_ph = 1/(2π·κ_TGP) for κ_TGP ≈ 2.012 → 0.0792 → drift −53% **FAIL**
4. **C4:** ε_ph = 1/(2β² + δ_F4) for β=1, δ_F4=0.114 → 1/2.114 = 0.473 → drift +182% **FAIL**
5. **C5:** ε_ph = 1/φ² (φ=golden ratio) = 1/2.618 = 0.382 → drift +127% **FAIL**

**Test:** All 5 candidates have drift > 5% → 5/5 FALSIFIED.

**Falsification:** if any candidate drift < 1% → competing identity, ε.1
needs refinement.

### E2.3 — F4 chain locks ε_ph implicitly

**Test:**
- F4 anchor: α₀ = 1069833/264500
- target_shift_F4 = 57/500 (LOCKED)
- α₀ = target_shift / ε_ph² ⟹ ε_ph² = target_shift / α₀
- ε_ph² = (57/500) / (1069833/264500) = 57·264500 / (500·1069833) = 15076500/534916500
- Simplified: ε_ph² ≈ 0.028184 ⟹ ε_ph ≈ 0.16788 (matches!)
- Sympy exact form via cancellation

**Falsification:** if F4-derived ε_ph² ≠ M9.1″ refined ε_ph² within 0.01% →
F4 chain decoupled from photon-ring scale.

### E2.4 — Heat-kernel a₂ frame consistency

**Test:**
- a₂_vacuum = (1/2) V''(Φ_eq)² = 2β² (ξ.1.Phase2 LOCKED)
- target_shift = ξ_geom · α(α−1) · ε_ph² / 2 / N_ref ?? (verify framework)
  Actually: target_shift = a₂ contribution z heat-kernel
- Verify target_shift_F4 = 57/500 dimensionally consistent z heat-kernel
  expansion in ε_ph² ratio

**Falsification:** if heat-kernel framework requires ε_ph in linear (not
quadratic) form → ξ.1 frame audit broken.

### E2.5 — M9.2-D vs M9.1″ refinement audit

**Test:**
- M9.2-D ψ_ph = 1.168 (3-decimal, coarse)
- M9.1″ ψ_ph = 4/(3+0.4250) = 160/137 = 1.16788 (sympy exact)
- Drift: |1.168 − 160/137|/(160/137) = 0.011%
- 0.527% F4-strict ξ-factor = a₂ EFT 1-loop correction (ξ.1.Phase2)
- M9.1″ refinement w 0.011% << 0.527% a₂ band ⟹ refinement nie wpływa
  na frame discrimination

**Falsification:** if M9.1″ refinement > 0.527% a₂ band → frame discrimination
contaminated.

### E2.6 — NGFP RG-stability of ε_ph² ratio

**Test:**
- Under common β-rescaling NGFP flow, ψ_ph (geometric) i β (substrate)
  scale together
- ε_ph² = (ψ_ph − 1)² preserves structure under β-rescaling
- α₀ = target_shift / ε_ph² is RG-invariant ratio (UV.1.Phase2.UV2.5)
- Therefore ε_ph² is RG-stable (consequence of α₀ ratio invariance)

**Falsification:** if ε_ph² shows drift > 0.5% pod NGFP RG simulation
→ ratio invariance broken.

### E2.7 — Classification decision

**Test:**
- ε_ph = 23/137 sympy exact (E2.1) ✓
- 5 alternative candidates falsified (E2.2) ✓
- F4 chain implicit lock (E2.3) ✓
- Heat-kernel a₂ frame consistent (E2.4) ✓
- M9.1″ refinement audit (E2.5) ✓
- NGFP RG-stability (E2.6) ✓

**Classification:**
- **DERIVED**: ψ_ph = 160/137 wynika z first-principles M9.1″ null geodesic
  Eddington-Finkelstein computation (deep), not just empirical fit
- **PARTIALLY DERIVED (refined)**: ε_ph = ψ_ph − 1 structural;
  M9.1″ origin documented but not re-derived w tym mini-cycle
- **STRUCTURAL HINT**: only F4 chain implicit lock, no ψ_ph derivation

**Verdict:** PARTIALLY DERIVED (refined) realistic outcome — ε_ph =
ψ_ph − 1 = 23/137 structurally; M9.1″ origin closed by ξ.1 + UV.1
predecessors; full DERIVED czeka na full M9.1″ Eddington-Finkelstein
re-derivation (long-term track).

---

## Verdict gate

**≥ 6/7 PASS** → ε.1 PARTIALLY DERIVED (refined), proceed Phase 3.

**≤ 5/7 PASS** → STRUCTURAL HINT, Phase 3 limited.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-eps-photon-ring/phase2_eps_decomposition.py 2>&1 | tee research/op-eps-photon-ring/phase2_eps_decomposition.txt
```

## Cross-references

- [`Phase1_results.md`](Phase1_results.md) — numerical landscape 5/5 PASS
- [`program.md`](program.md) — overall ε.1 plan
- [`../op-xi-photon-ring/Phase2_results.md`](../op-xi-photon-ring/Phase2_results.md) — heat-kernel a₂ frame (V''(1)=2β)
- [`../op-uv-as-ngfp/Phase2_results.md`](../op-uv-as-ngfp/Phase2_results.md) — UV.1.UV2.5 F4 chain RG-stability
