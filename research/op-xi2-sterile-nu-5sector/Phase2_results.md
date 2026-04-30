---
title: "ξ.2.Phase2 results — first-principles B²_sterile + |U_e4|² + Δm²_{41} (7/7 PASS)"
date: 2026-04-30
cycle: ξ.2.Phase2
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - xi2
  - phase2
  - sterile-nu
  - PASS
---

# ξ.2.Phase2 results — first-principles B²_sterile + |U_e4|² + Δm²_{41}

**Score: 7/7 PASS → Phase 3 viable; Form A LOCKED primary.**

> **Headline:** B²_sterile = 0 LOCKED (TGP 4-sector minimal counting);
> |U_e4|² = 0 LOCKED; m_4 = 0 LOCKED; RAA + Gallium 4σ anomalies → NOT
> sterile (flux + ⁷¹Ge systematics); 5 alternative SBL sterile fits
> FALSIFIED; thermalized sterile excluded 6.2σ z Planck N_eff;
> Form A primary z Σm_ν = 59.01 meV unchanged.

## Sub-test results

### X2.1 — Form A (B²_sterile=0) null structural ✓ PASS

- B²_sterile = 0 (sympy-exact)
- |U_e4|² = 0
- sin²(2θ_14) = 0
- m_4 = 0

### X2.2 — Form B (B²_sterile=λ_C²) sin²(2θ)=λ_C⁴ ✓ PASS

- B²_sterile = λ_C² = 81/1600 ≈ 0.05063
- sin²(2θ_14) = λ_C⁴ = 6561/2560000 ≈ 0.002563
- < STEREO 2023 95% CL = 0.07 ✓

### X2.3 — Δm²_{41} Form B « SBL scale ✓ PASS

- Form B: m_4 ~ λ_C · m_3 = 0.225 · 49.53 meV ≈ 0.0111 eV
- Δm²_{41} ~ 1.24·10⁻⁴ eV² « SBL ~1 eV² (factor 10⁵ short)
- → Form B CANNOT natively explain SBL anomalies

### X2.4 — Thermalized sterile >5σ N_eff excluded ✓ PASS

- Planck 2018: N_eff = 2.99 ± 0.17
- 1 thermalized sterile: N_eff = 4.05 → 6.2σ excluded
- Form A (m_4=0): N_eff = 3.046 SM consistent
- Form B at sub-eV: thermalized variant cosmologically excluded

### X2.5 — 5 alt SBL sterile fits FALSIFIED ✓ PASS

| Form | Δm² (eV²) | |U_e4|² | Conflict |
|------|-----------|--------|----------|
| (i) sterile-3+1 RAA | 1.7 | 0.02 | STEREO 2023 95% CL |
| (ii) gallium 3+1 | 4.0 | 0.10 | STEREO + KATRIN 2022 |
| (iii) eV-Dirac ν_R | 1.0 | 0.05 | KATRIN 2022 |
| (iv) flav-democratic light | 0.5 | 0.05 | Daya Bay/RENO flux 2024 |
| (v) heavy ν_R seesaw | 10¹⁵ | 0 | irrelevant SBL/KATRIN |

### X2.6 — Σm_ν Form A primary, Form B falsifiable ✓ PASS

- Form A: Σm_ν = 59.01 meV (z ζ.1, unchanged)
- Form B (m_4 ~0.1 eV thermalized): Σm_ν ~159 meV → DESI DR2 < 72 meV FAILS
- → Form A only TGP-consistent option

### X2.7 — Classification cascade 5 promotions Form A ✓ PASS

1. **B²_sterile = 0** LOCKED (TGP 4-sector minimal counting)
2. **|U_e4|² = 0** LOCKED (no sterile mixing)
3. **m_4 = 0** LOCKED (no sterile mass)
4. **RAA → flux systematics** LOCKED (NOT sterile)
5. **Gallium 4σ → ⁷¹Ge(p,n) cross-section** LOCKED (NOT sterile)

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../../INDEX.md]]
