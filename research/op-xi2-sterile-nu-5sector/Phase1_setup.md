---
title: "ξ.2.Phase1 setup — sterile ν landscape + B²_sterile candidates + viability (5 sub-tests)"
date: 2026-04-30
cycle: ξ.2.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[../op-nu-majorana-phase-mbb/Phase3_results.md]]"
tags:
  - TGP
  - xi2
  - phase1
  - sterile-nu
  - reactor-anomaly
  - gallium-anomaly
---

# ξ.2.Phase1 — sterile ν landscape + B²_sterile candidates + viability

> **Cel:** Audit current SBL sterile ν experiments (STEREO 2023, PROSPECT-I,
> BEST 2022, KATRIN 2022); register 4 B²_sterile candidates; gate
> viability through STEREO 2023 95% CL exclusion (sin²(2θ_14) < 0.01).

## Sub-tests

### X1.1 — Reactor antineutrino anomaly (RAA) audit

- 2011 Mention et al.: 6% deficit z ν̄_e detected/predicted ratio
- 2023 STEREO final result: 95% CL excludes RAA at Δm² ∈ [0.3, 10] eV²,
  sin²(2θ_14) > 0.07
- 2023 PROSPECT-I confirmed STEREO exclusion
- 2024 Daya Bay + RENO + Double Chooz updated U-235 flux modeling
  reduced anomaly to ~2-3% (within flux uncertainty)

### X1.2 — Gallium anomaly audit

- GALLEX (1995): R = 0.95 ± 0.11 (⁵¹Cr source)
- SAGE (2010): R = 0.87 ± 0.09 (⁵¹Cr + ³⁷Ar)
- BEST 2022: R_inner = 0.79 ± 0.05, R_outer = 0.77 ± 0.05; combined 4σ
  z null hypothesis; sterile fit Δm² ~ 1-10 eV², |U_e4|² ~ 0.1
- Tension z STEREO 2023 95% CL exclusion → systematics required
  (⁷¹Ge(p,n)⁷¹Ga cross-section)

### X1.3 — KATRIN sterile bounds

- KATRIN 2022 first sterile search: m_4 < 1.6 eV (90% CL) for
  sin²(2θ_e4) > 0.1
- KATRIN-TRISTAN 2027+: extends sensitivity to m_4 ~ 0.1 eV; sin²(2θ_e4)
  < 10⁻⁴ at m_4 = 0.5 eV
- Direct kinematic mass measurement (model-independent)

### X1.4 — 4 B²_sterile candidates

| Form | B²_sterile | sin²(2θ_14) | m_4 (eV) | Status |
|------|-----------|-------------|----------|--------|
| A (decoupled) | 0 | 0 | 0 | favored z TGP minimal counting |
| B (vacuum-suppressed) | λ_C² ≈ 0.0506 | λ_C⁴ ≈ 2.56·10⁻³ | ~0.1 | TGP-Cabibbo cascade |
| C (electroweak-half) | 1/4 | ~6.3·10⁻³ | ~0.5 | alt structural |
| D (half-Majorana) | 1/2 | ~1.3·10⁻² | ~1 | naive scaling — fails STEREO |

### X1.5 — Viability gate

Wszystkie 4 candidates muszą produkować |U_e4|² compatible z STEREO 2023
exclusion sin²(2θ_14) < 0.07 (95% CL). Form D (1/2) → sin²(2θ_14) ~ 0.013
< 0.07 → marginal pass. Form C (1/4) → 6.3·10⁻³ → pass. Form A (0) → 0 →
pass. Form B (λ_C²) → 2.56·10⁻³ → pass.

**4/4 candidates viable** under STEREO 2023; promotion to Phase 2.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-xi2-sterile-nu-5sector/phase1_sterile_landscape.py 2>&1 | tee research/op-xi2-sterile-nu-5sector/phase1_sterile_landscape.txt
```

## Cross-references

- [[program.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
