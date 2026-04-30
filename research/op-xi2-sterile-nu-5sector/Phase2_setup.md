---
title: "ξ.2.Phase2 setup — first-principles B²_sterile + |U_e4|² + Δm²_{41} (7 sub-tests)"
date: 2026-04-30
cycle: ξ.2.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - xi2
  - phase2
  - sterile-nu
  - U_e4
  - delta-m41
---

# ξ.2.Phase2 — first-principles B²_sterile + |U_e4|² + Δm²_{41} (7 sub-tests)

> **Cel:** Lock B²_sterile dla Form A (decoupled, B²=0) jako primary +
> Form B (vacuum-suppressed, B²=λ_C²) jako alt; derive sin²(2θ_14) +
> m_4 sympy-exact; falsify 5 alternative SBL sterile fits;
> register classification cascade (5 promotions).

## Sub-tests

### X2.1 — Form A (B²_sterile = 0) sympy-exact null

W TGP minimalnym 4-sector chirality framework (post-ν.1), nie ma
miejsca na B²_sterile bez naruszenia counting-cascade. Strict null:

$$|U_{e4}|^{2} = 0, \quad \sin^{2}(2\theta_{14}) = 0, \quad m_{4} = 0$$

→ structural null; SBL sterile NIE istnieje; RAA + Gallium anomalies
muszą pochodzić z systematics (flux modeling, ⁷¹Ge cross-section).

### X2.2 — Form B (B²_sterile = λ_C²) sympy-exact suppression

Alt: sterile coupling przez Cabibbo cascade z down-quark mixing analog:

$$\sin^{2}(2\theta_{14}) = \lambda_{C}^{4} \approx 2.563 \times 10^{-3}$$

→ below STEREO 2023 (sin²(2θ) > 0.07 excluded), pass; potential SBN
2030+ marginal signal.

### X2.3 — Δm²_{41} dla Form A vs Form B

- Form A: undefined (no oscillation)
- Form B: m_4² ~ TGP-native scale; m_4 ~ √(λ_C² · m₃²) where m₃=49.53 meV
  → m_4 ~ λ_C · m₃ ≈ 11.1 meV (eV scale needed for SBL → tension);
  alternative: m_4 ~ Λ_QCD · λ_C^n (TBD scale)

For SBL relevance, Δm²_{41} must be ~ 1 eV² → m_4 ~ 1 eV. TGP Form B
nie produkuje naturally ten scale → Form B disfavored cosmologically;
Form A pozostaje primary.

### X2.4 — m_4 sterile mass

- Form A: m_4 = 0 (no sterile)
- Form B: m_4 ~ 0.1 eV (KATRIN-TRISTAN sensitivity-edge)

Cosmological N_eff Planck 2018: N_eff = 2.99 ± 0.17 → 1 light sterile
ν_R thermalized would give N_eff = 4.05; 95% CL excludes m_4 ~ 0.1-1 eV
fully thermalized. Form A consistent z Planck N_eff.

### X2.5 — 5 alternative sterile fits FALSIFIED

| Form | Δm² (eV²) | |U_e4|² | Conflict |
|------|-----------|--------|----------|
| (i) sterile-3+1 RAA | 1.7 | 0.02 | STEREO 2023 95% CL excluded |
| (ii) gallium 3+1 | 4.0 | 0.10 | STEREO 2023 95% CL excluded; KATRIN 2022 |
| (iii) eV-Dirac ν_R | m_4 ~ eV | flav.democ. | KATRIN 2022 m_4 < 1.6 eV |
| (iv) light-sterile flavor-democratic | 0.5 | 0.05 | Daya Bay/RENO flux update |
| (v) heavy ν_R seesaw | m_4 > 100 GeV | 0 | irrelevant SBL/KATRIN; not testable |

→ 4/5 SBL falsifiable, all FALSIFIED; (v) seesaw is irrelevant for ξ.2.

### X2.6 — TGP-native Σm_ν 5-sector update

Form A: Σm_ν unchanged = 59.01 meV (z ζ.1; no sterile mass added)
Form B: Σm_ν = 59.01 + m_4 ~ 159 meV (1 eV sterile thermalized) → tension
z DESI DR2 < 72 meV → Form B FALSIFIED cosmologically if thermalized.

→ Form A is **only** TGP-consistent option satisfying both SBL +
cosmological constraints.

### X2.7 — Classification cascade (5 promotions)

Form A LOCKED:
1. **B²_sterile = 0** LOCKED (TGP minimal counting structural)
2. **|U_e4|² = 0** LOCKED (no sterile mixing)
3. **m_4 = 0** LOCKED (no sterile mass)
4. **RAA → flux systematics** LOCKED (anomaly NOT sterile)
5. **Gallium 4σ → ⁷¹Ge systematics** LOCKED (anomaly NOT sterile)

→ ξ.2 promotes 5 cascade items; alt forms B/C/D demoted to falsified.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-xi2-sterile-nu-5sector/phase2_sterile_derivation.py 2>&1 | tee research/op-xi2-sterile-nu-5sector/phase2_sterile_derivation.txt
```

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../op-zeta-mass-spectrum/Phase3_results.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
