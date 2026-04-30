---
title: "ο.1 program — cosmological Σm_ν tension hardening (DESI DR2/DR3 vs TGP NH)"
date: 2026-04-30
cycle: ο.1
status: ACTIVE
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-xi2-sterile-nu-5sector/Phase3_results.md]]"
tags:
  - TGP
  - omicron1
  - program
  - cosmology
  - sigma-mnu
  - DESI
  - CMB-S4
---

# ο.1 program — cosmological Σm_ν tension hardening

> **Cel:** Harden ν.1-derived Σm_ν_TGP_A = 59.01 meV (NH, Form A,
> m_1 = 0) against DESI DR2 (2024–2025) cosmological bound
> Σm_ν < 72 meV (95% CL) + DESI DR3 2027+ tightening + CMB-S4 2030+
> ~20 meV ultimate sensitivity. Falsify 6 alt mass-spectrum fits
> (IH, degenerate, quasi-degenerate, sub-NH-floor, sterile-augmented,
> early-DE-modified). Close ο.1 z 7-channel cosmological falsification
> convergence including DESI DR3 / Simons / CMB-S4 / KATRIN / Euclid /
> LiteBIRD / Roman.

## Predecessors

- **ν.1** Phase 1–3 → Σm_ν_A = 59.01 meV (NH, m_1 = 0)
  derived z TGP NO mass spectrum + Δm²_21, Δm²_31 PDG anchors.
- **ξ.2** Phase 1–3 → no thermalized sterile (N_eff < 3.05);
  sterile sector decoupled → Σm_ν_cosmo = Σm_ν_active strict.
- **μ.1** Phase 1–3 → δ_CP and PMNS angles 8 free → 0 free; mass
  hierarchy ordering NH locked by chirality assignment.
- **κ.1** Phase 1–3 → λ_C = 0.225 numerical anchor for Form B fallback.

## Hipoteza centralna

**TGP 4-sector chirality framework + ν.1 Majorana NH structure forces
Σm_ν_A = 59.01 meV w wąskim oknie [55, 65] meV** ze:
- **m_1 = 0** (Form A: lightest neutrino massless from TGP chirality lock)
- **m_2 = √Δm²_21 ≈ 8.61 meV** (solar splitting)
- **m_3 = √Δm²_31 ≈ 50.50 meV** (atmospheric splitting)
- **Σm_ν_A = 0 + 8.61 + 50.50 ≈ 59.11 meV** (PDG anchors)
- **TGP-corrected: 59.01 meV** (z ν.1 substrate-action mass refinement)

**Form B alternative** (m_1 = ζ_TGP ≠ 0): Σm_ν_B ≈ 70–80 meV w pobliżu
DESI DR2 95% CL → tension trigger.

## DESI DR2 2024–2025 status

Z DESI DR2 + Planck CMB:

| Bound | Σm_ν (95% CL) | Status |
|-------|---------------|--------|
| Planck 2018 alone | < 0.12 eV | trivial PASS |
| Planck + BAO BOSS | < 0.15 eV | trivial PASS |
| Planck + DESI DR1 | < 0.082 eV | PASS w 23 meV margin |
| **Planck + DESI DR2** | **< 0.072 eV** | **PASS w 13 meV margin** |
| DESI DR2 + CMB + SN | < 0.064 eV | borderline (5 meV margin) |

→ TGP Form A 59.01 meV pass DESI DR2; Form B ~70+ meV at edge.

## Plan 3 fazowy (5+7+6 = 18 sub-tests)

### Phase 1 — landscape (5 sub-tests, **PASS bramka ≥4/5**)

- **O1.1** Σm_ν cosmological bounds inventory (Planck, BAO, DESI, eBOSS, SN)
- **O1.2** TGP Σm_ν_A derivation review (m_1=0, m_2, m_3 z PDG)
- **O1.3** NH vs IH likelihood ratio (DESI DR2 disfavors IH ~3σ)
- **O1.4** m_lightest constraint (Form A m_1=0 vs Form B m_1=ζ_TGP)
- **O1.5** viability gate Σm_ν_TGP_A < 72 meV (DESI DR2)

### Phase 2 — derivation hardening (7 sub-tests, **PASS bramka ≥6/7**)

- **O2.1** Σm_ν_A = m_1 + m_2 + m_3 closed form z TGP NO + PDG Δm²
- **O2.2** Form A m_1=0 vs Form B m_1=ζ_TGP comparison
- **O2.3** Cosmological active-only mass (post-ξ.2 sterile decouple)
- **O2.4** Σm_ν_TGP vs DESI DR2 < 72 meV → 13 meV margin
- **O2.5** Σm_ν_TGP vs Planck+BAO < 0.15 eV → trivial PASS
- **O2.6** Tension analysis: Form A vs DESI central ~ 1.4σ
- **O2.7** 6 alt mass-spectrum fits FALSIFIED + 4 promotions cascade
  (Σm_ν_A=59.01 meV, m_1=0, NH-LOCK, m_lightest=0 LOCKED)

### Phase 3 — predictions + falsification (6 sub-tests, **PASS bramka ≥5/6**)

- **O3.1** DESI DR3 2027+ ~50 meV projection: Form A 9 meV tension (~2σ)
- **O3.2** Simons Observatory 2025+ projection
- **O3.3** CMB-S4 2030+ ~20 meV ultimate: Form A decisive 5σ falsification
- **O3.4** KATRIN 2030+ direct mass m_β ~ 0.2 eV bound (orthogonal)
- **O3.5** Euclid + Roman 2027–2030+ weak lensing + galaxy clustering
- **O3.6** 7-channel ο.1 cosmological falsification convergence
  (DESI DR2 confirmed + DR3 + Simons + CMB-S4 + KATRIN + Euclid + LiteBIRD)

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-omicron1-sigmamnu-cosmo/phase1_sigmamnu_landscape.py 2>&1 | tee research/op-omicron1-sigmamnu-cosmo/phase1_sigmamnu_landscape.txt
PYTHONIOENCODING=utf-8 python -X utf8 research/op-omicron1-sigmamnu-cosmo/phase2_sigmamnu_derivation.py 2>&1 | tee research/op-omicron1-sigmamnu-cosmo/phase2_sigmamnu_derivation.txt
PYTHONIOENCODING=utf-8 python -X utf8 research/op-omicron1-sigmamnu-cosmo/phase3_sigmamnu_predictions.py 2>&1 | tee research/op-omicron1-sigmamnu-cosmo/phase3_sigmamnu_predictions.txt
```

## Cross-references

- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]]
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
