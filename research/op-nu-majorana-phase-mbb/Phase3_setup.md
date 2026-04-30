---
title: "ν.1.Phase3 setup — predictions + falsification convergence (6 sub-tests)"
date: 2026-04-30
cycle: ν.1.Phase3
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
tags:
  - TGP
  - nu1
  - phase3
  - predictions
  - falsification
  - mbb
---

# ν.1.Phase3 — predictions + falsification convergence (6 sub-tests)

> **Cel:** Register concrete 0νββ falsification predictions across 4
> next-/next-next-gen experiments (KamLAND-Zen 2027+, LEGEND-1000 2030+,
> nEXO 2030+, NEXT-HD 2030+); register ★ headline PMNS 8 free → 0 free
> + 2 Majorana phases DERIVED dual; outline ξ.2/ο.1 future research-track;
> close ν.1 z 7-channel falsification convergence.

## Sub-tests

### N3.1 — KamLAND-Zen 2027+ ~5 meV both forms < limit

m_ββ_TGP_A = 1.584 meV < 5 meV ✓ no falsification
m_ββ_TGP_B = 3.249 meV < 5 meV ✓ no falsification

Status: **2027+ generation cannot falsify** any TGP form (both well under
5 meV); but if KamLAND-Zen 2027+ excludes m_ββ < 5 meV (positive signal),
it would constrain Form B closer to the ceiling.

### N3.2 — LEGEND-1000 2030+ ~3 meV Form A evades, Form B at edge

m_ββ_TGP_A = 1.584 meV < 3 meV ✓ Form A evades
m_ββ_TGP_B = 3.249 meV ≈ 3 meV (slightly above) → Form B at edge,
**potential 1σ tension** if LEGEND-1000 reaches 2 meV sensitivity (mild
falsification trigger).

### N3.3 — nEXO + NEXT-HD 2030+ ~0.5 meV discriminates 3σ

Δm_ββ = m_ββ_B − m_ββ_A = 3.249 − 1.584 = **1.665 meV** = 3.33×
sensitivity → **3σ separation** between Form A vs Form B (gap factor 2.05).
nEXO + NEXT-HD 2030+ are the **decisive discriminator** for ν.1
hypothesis (B).

### N3.4 — ★ Headline: combined PMNS 6 free → 0 free + 2 Majorana DERIVED dual

Pre-program (post-ι.1 2026-04 CKM closure):
- 8 fundamental PMNS-related parameters: 3 angles + 1 δ_CP + 4 (m_lightest +
  Δm²₂₁ + |Δm²₃₁| + ordering)

Post-μ.1 (2026-04):
- 3 angles refined² (s²₁₂ = 5149/16800, s²₂₃ = 4/7, s²₁₃ = 13627867/624000000)
- δ_CP dual (Form A 205.36° + Form B 259.82°) — DERIVED
- m_lightest + masses fixed z ζ.1 NO ordering

Post-ν.1 (2026-04):
- α₂₁ dual: π/2 (Form A) ‖ 11π/13 (Form B) — DERIVED
- α₃₁ dual: 9π/26 (Form A) ‖ 12π/7 (Form B) — DERIVED

**Net: 8 fundamental → 0 free + 2 Majorana phases DERIVED dual structural
form.** Combined PMNS sector closure with structural dual ambiguity
discriminable z 0νββ 2030+ experiments.

### N3.5 — ξ.2/ο.1 future research-track outline

Future mini-cycles z ν.1 implications:
1. **ξ.2 — sterile ν 5-sector extension B²_sterile**: extend 4-sector
   chirality framework do 5-sector z B²_sterile (testable via short-baseline
   reactor anomaly + STEREO/PROSPECT data).
2. **ο.1 — cosmological Σm_ν tension hardening**: ν.1 ζ.1 m₁ = 0.76 meV
   gives Σm_ν = 59.01 meV ≈ DESI DR2 limit Σ < 72 meV (margin +22%);
   DESI DR3 2027+ refined limit will harden this constraint.
3. **π.1 — neutrinoless double beta decay TGP cross-checks** (Te-130 vs
   Xe-136 vs Ge-76 isotope-dependent NME corrections in TGP).

### N3.6 — 7-channel ν.1 falsification convergence

| # | Channel | Observable | Sensitivity | Date |
|---|---------|------------|-------------|------|
| 1 | KamLAND-Zen 2027+ | m_ββ < 5 meV | Form A/B both evade | 2027+ |
| 2 | LEGEND-1000 2030+ | m_ββ < 3 meV | Form A evades, Form B at edge | 2030+ |
| 3 | nEXO 2030+ | m_ββ < 0.5 meV | 3σ Form A vs B discrimination | 2030+ |
| 4 | NEXT-HD 2030+ | m_ββ < 0.5 meV | 3σ Form A vs B (Xe-136 HP) | 2030+ |
| 5 | DESI DR3 2027+ | Σm_ν = 59 meV | confirms ζ.1 NO ordering | 2027+ |
| 6 | DUNE δ_CP 2030+ | δ_CP Form A vs B | refines 205° vs 260° | 2030+ |
| 7 | T2K-II + HK 2027+ | δ_CP × m_ββ cross-check | combined Form A/B | 2027+ |

**Convergence:** ≥6/7 channels → ν.1 PASS; 7/7 → FULL CONVERGENCE.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-nu-majorana-phase-mbb/phase3_predictions.py 2>&1 | tee research/op-nu-majorana-phase-mbb/phase3_predictions.txt
```

## Cross-references

- [[program.md]]
- [[Phase2_results.md]]
- [[../op-mu-pmns-phase-hardening/Phase3_results.md]]
- [[../op-zeta-mass-spectrum/Phase3_results.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
