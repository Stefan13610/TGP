---
title: "ξ.1.Phase1 setup — EFT photon-ring frame audit"
date: 2026-04-29
cycle: ξ.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
tags:
  - TGP
  - xi-factor
  - photon-ring
  - audit
  - falsification
---

# ξ.1.Phase1 — EFT photon-ring frame audit

> **Cel:** Sprawdzić, czy 5 fundamental inputs do a₂ derivation
> (ξ_geom=1, α(α−1)=2, ψ_ph=1.168, F4 rational provenance, Phase 2 strict
> provenance) są wszystkie LOCKED w istniejących closurach. PASS = ξ.1
> derivation może wystartować z zero free parameters.

---

## Hipoteza Phase 1 (H_ξ1)

Wszystkie 5 inputów do photon-ring a₂ derivation są **already LOCKED**
w istniejących closurach (Phase 1.A, Phase 2.B, M9.1″, M9.2-D,
closure_2026-04-26). Brak free parametera w premise → ξ.1 Phase 2 może
proceed.

---

## 5 sub-testów Phase 1

### ξ1.1 — ξ_geom = 1 derivation z M9.1″ vacuum

**Cel:** Potwierdzić, że ξ_geom (vacuum geometric prefactor w
heat-kernel a₂ frame) jest LOCKED na 1.0 w M9.1″ FRW background.

**Test:**
- Phase 2.B.2 derivation: ξ_geom = 1 dla single-Φ + Z₂ vacuum
- Phase 1.F.3: V'(Φ_eq=1) = 0 exact (β=γ vacuum)
- Pod Φ_eq = 1 i K_geo = 1: ξ_geom = α(α−1)/2 = 2/2 = 1
- Falsification: jeżeli ξ_geom ma free parameter > 1% wariancji → frame nie zamknięty

**Anchor:** ξ_geom = 1.0 (Phase 2.B.2, vacuum-exact)

### ξ1.2 — α(α−1) = 2 derivation z Phase 1.A.1

**Cel:** Potwierdzić, że α=2 unique selection w K(φ) = K_geo·φ^α pod
(C1)–(C3) constraints.

**Test:**
- Phase 1.A.1 (Theorem alpha2): α=2 unique under (C1: stability) ∧ (C2: positivity) ∧ (C3: dimensional)
- α(α−1) = 2(1) = 2 exact integer
- Sympy: α=2 jest jedyne dopuszczalne rational w {1, 2, 3, 4, 5}
- Falsification: jeżeli inny α dopuszczony pod jakimkolwiek constraint relaxation → a₂ structure ambiguous

**Anchor:** α(α−1) = 2 (Theorem alpha2 LOCKED, F2 LOCKED)

### ξ1.3 — ψ_ph − 1 = 0.168 z M9.2-D photon-ring

**Cel:** Potwierdzić photon-ring deviation jest stable na 0.168
across all M9.2-D sub-cycles + M9.1″ derivation.

**Test:**
- M9.2-D: ψ_ph = 1.168 (universal photon-ring TGP)
- M9.1″: ψ_ph = 1.16788 (Phase 2.B.1, refined to 5 digits)
- ε_ph = ψ_ph − 1 = 0.168 (or 0.16788) — drift 0.07%
- ε_ph² = 0.028224
- Falsification: jeżeli ψ_ph drift > 0.5% w jakimkolwiek scenariuszu → photon-ring frame nie zamknięty

**Anchor:** ε_ph = 0.168 (M9.2-D LOCKED, drift < 0.07%)

### ξ1.4 — F4 rational provenance audit

**Cel:** Potwierdzić, że F4 sympy rational 1069833/264500 wynika z
arithmetic identity α₀ = Δ_target / [(ψ_ph−1)² · ξ_geom] z
Δ_target = 0.114, ψ_ph = 1.16788, ξ_geom = 1.0.

**Test:**
- Phase 2.B.3 chain: α₀ = 0.114 / (0.16788² · 1.0) = 0.114 / 0.028184 ≈ 4.0447
- Sympy rational: 1069833/264500 = 4.044737… (exact)
- 0.114 jest input (target_shift normalization w a₂ frame)
- 0.16788 jest M9.1″ photon-ring deviation
- 1.0 jest ξ_geom z vacuum
- Falsification: jeżeli F4 rational nie z tej arithmetic identity → trace inconsistency

**Anchor:** F4 = 1069833/264500 = α₀_F4 ≈ 4.04474 (closure_2026-04-26 LOCKED)

### ξ1.5 — Phase 2 strict (1/2)(1 − 3/3.88) provenance

**Cel:** Potwierdzić, że Phase 2 strict target_shift = 0.1134 wynika z
direct geometric photon-ring shift formula (1/2)(1 − r_ph^GR/r_ph^TGP).

**Test:**
- BH.1.Phase2.T2.4: target_shift = (1/2)(1 − 3/3.88) = 0.113402
- 3.88 = (1 + ψ_ph)·something — actually r_ph^TGP/M = 3.88, r_ph^GR/M = 3
- ε_ph² = 0.168² = 0.028224
- α₀_strict = 0.113402 / 0.028224 = 4.0179 (Phase 2 strict)
- 0.114 (F4 frame) − 0.1134 (strict frame) = 0.0006 (split 0.5%)
- Strict form NOT same jako F4 — różny calibration w prefactorze
- Falsification: jeżeli strict nie z direct geometric formula → strict not bare-form

**Anchor:** Phase 2 strict form (1/2)(1 − 3/3.88) (BH.1.Phase2 LOCKED)

---

## Verdict gate

5/5 PASS → all 5 inputs LOCKED → Phase 2 a₂ derivation może proceed
z zero free parameters w premise.

≤ 4/5 PASS → audit gap identifier; Phase 2 deferred do gap-fill.

---

## Materiał wykonawczy

- **Skrypt:** [`phase1_frame_audit.py`](phase1_frame_audit.py) (5 sub-tests, sympy rational + numeric)
- **Output:** `phase1_frame_audit.txt`
- **Memo:** `Phase1_results.md`

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 TGP/TGP_v1/research/op-xi-photon-ring/phase1_frame_audit.py 2>&1 | tee TGP/TGP_v1/research/op-xi-photon-ring/phase1_frame_audit.txt
```

---

## Constants used

```
xi_geom_locked    = 1.0           (Phase 2.B.2 vacuum LOCKED)
alpha_K           = 2             (F2 LOCKED, Theorem alpha2)
alpha_K_x_alpha_K_m_1 = 2         (= α(α−1) = 2·1 = 2)
psi_ph            = 1.168         (M9.2-D)
psi_ph_M9_1       = 1.16788       (M9.1″ refined)
eps_ph            = 0.168
eps_ph_sq         = 0.028224
target_shift_F4    = 0.114        (F4 rational input)
target_shift_strict = (1/2)(1 - 3/3.88) = 0.113402
F4_rational       = 1069833/264500 ≈ 4.04472
alpha_0_strict    = 4.0179
split_ratio       = (0.114 - 0.1134) / 0.1134 = 0.5%
```

---

## Cross-references

- [`program.md`](program.md) — overall ξ.1 plan
- [`../op-cross-sector-charge/Phase2_results.md`](../op-cross-sector-charge/Phase2_results.md) — XS.1 ξ unresolution
- [`../op-bh-alpha-threshold/Phase2_results.md`](../op-bh-alpha-threshold/Phase2_results.md) — Phase 2 strict provenance
- [`../op-phase2-quantum-gravity/Phase2_B_results.md`](../op-phase2-quantum-gravity/Phase2_B_results.md) — F4 rational chain
- [`../op-phase3-uv-completion/Phase3_E_results.md`](../op-phase3-uv-completion/Phase3_E_results.md) — UV7 a₂ frame
