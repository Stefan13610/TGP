---
title: "π.1 program — 0νββ NME isotope cross-checks (Te-130 / Xe-136 / Ge-76)"
date: 2026-04-30
cycle: π.1
status: ACTIVE
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-omicron1-sigmamnu-cosmo/Phase3_results.md]]"
tags:
  - TGP
  - pi1
  - program
  - 0nubb
  - mbb
  - NME
  - cross-isotope
---

# π.1 program — 0νββ NME isotope cross-checks

> **Cel:** Harden ν.1-derived m_ββ_TGP_A = 1.584 meV ‖ m_ββ_TGP_B = 3.249 meV
> against full nuclear-matrix-element (NME) systematics across the three
> dominant 0νββ isotopes — Ge-76 (LEGEND), Te-130 (CUORE/CUPID),
> Xe-136 (KamLAND-Zen / nEXO / NEXT-HD). Cross-isotope T_{1/2} ratios
> cancel m_ββ and isolate the NME × phase-space-factor (PSF) systematic;
> require TGP cross-isotope consistency. Close π.1 z 7-channel
> NME-isotope falsification convergence (KZ-900 / LEGEND-1000 / nEXO /
> NEXT-HD / CUPID / SNO+ / theory NME-cross).

## Predecessors

- **ν.1** Phase 1–3 → m_ββ_TGP_A = 1.584 meV ‖ m_ββ_TGP_B = 3.249 meV
  derived z PMNS angles + Majorana phases dual structural form.
- **ο.1** Phase 1–3 → Σm_ν_A = 59.01 meV LOCKED w 18% DR2 margin;
  CMB-S4 5σ DECISIVE 2030+; constrains m_lightest = 0 (Form A NH).
- **ξ.2** Phase 1–3 → no thermalized sterile; sterile-channel decay
  to 0νββ rate negligible (|U_e4|² = 0).
- **κ.1** Phase 1–3 → λ_C = 0.225 anchor (Form B fallback m_lightest).

## Hipoteza centralna

**TGP m_ββ Form A (1.584 meV) and Form B (3.249 meV) generate a fully
specified isotope half-life landscape through the standard mass-mechanism**:

$$
\left[T_{1/2}^{0\nu}\right]^{-1}
\;=\; G_{0\nu}^{(A)} \cdot \big|M_{0\nu}^{(A)}\big|^2 \cdot
      \left(\frac{m_{\beta\beta}}{m_e}\right)^2
$$

Cross-isotope ratios cancel m_ββ entirely:

$$
\frac{T_{1/2}^{(A)}}{T_{1/2}^{(B)}}
\;=\;
\frac{G_{0\nu}^{(B)} \cdot |M_{0\nu}^{(B)}|^2}
     {G_{0\nu}^{(A)} \cdot |M_{0\nu}^{(A)}|^2}
$$

→ **NME-systematic-only test** independent of TGP m_ββ value.
The TGP Form A vs Form B 2.05× factor in m_ββ drops to a uniform 4.21×
factor in T_{1/2} **across all isotopes simultaneously** — falsifiable
via simultaneous detection across two isotopes.

## NME literature anchors (g_A = 1.27 unquenched)

| Isotope | QRPA | IBM-2 | NSM | EDF |
|---------|------|-------|-----|-----|
| Ge-76 (Z=32) | 5.0 | 4.6 | 3.0 | 4.6 |
| Te-130 (Z=52) | 4.2 | 4.1 | 1.9 | 5.1 |
| Xe-136 (Z=54) | 3.4 | 3.3 | 1.8 | 4.2 |

Phase-space factors G_{0ν} (10⁻¹⁵ yr⁻¹, Kotila & Iachello 2012):

| Isotope | G_{0ν} |
|---------|--------|
| Ge-76 | 2.36 |
| Te-130 | 14.22 |
| Xe-136 | 14.58 |

## Plan 3 fazowy (5+7+6 = 18 sub-tests)

### Phase 1 — landscape (5 sub-tests, **PASS bramka ≥4/5**)

- **P1.1** NME inventory across QRPA / IBM-2 / NSM / EDF for 3 isotopes
- **P1.2** Phase-space factor G_{0ν} inventory + g_A quenching status
- **P1.3** m_ββ Form A=1.584 / Form B=3.249 meV anchored z ν.1
- **P1.4** Method-spread (max/min) NME systematic σ inventory
- **P1.5** Viability gate — both forms above current 90% CL bounds for all
  isotopes × all NME methods (no rate-only falsification)

### Phase 2 — derivation hardening (7 sub-tests, **PASS bramka ≥6/7**)

- **P2.1** T_{1/2}(iso, Form, NME) closed form + 18-cell matrix
- **P2.2** Cross-isotope T_{1/2} ratios (cancel m_ββ → pure NME×PSF)
- **P2.3** Universal Form A/B T_{1/2} factor 4.21× (= 2.05²)
- **P2.4** TGP-native NME estimator via B²·Z·A^{1/3} cascade scaling
- **P2.5** NME method best-match vs TGP (which method survives)
- **P2.6** 6 alt fits FALSIFIED (g_A quench=0.7, single-method, IH m_lightest, etc.)
- **P2.7** 4-way LOCK promotions (m_ββ_A=1.584 LOCKED, m_ββ_B=3.249 LOCKED,
  cross-isotope universality 4.21× LOCKED, Form A/B-isotope discrimination)

### Phase 3 — predictions + falsification (6 sub-tests, **PASS bramka ≥5/6**)

- **P3.1** KamLAND-Zen 800 + KZ-900 2027+ Xe-136 sensitivity vs both forms
- **P3.2** LEGEND-1000 2030+ Ge-76 sensitivity vs both forms
- **P3.3** nEXO 2030+ Xe-136 + NEXT-HD 2030+ Xe-136 (HP gas)
- **P3.4** CUPID 2030+ Te-130 + SNO+ Te-130 cross-isotope
- **P3.5** Cross-isotope T_{1/2} ratio falsification (3-way Ge/Te/Xe consistency)
- **P3.6** 7-channel π.1 falsification convergence
  (KZ-900 + LEGEND-1000 + nEXO + NEXT-HD + CUPID + SNO+ + theory NME-cross)

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-pi1-bb0nu-nme-isotope/phase1_nme_landscape.py 2>&1 | tee research/op-pi1-bb0nu-nme-isotope/phase1_nme_landscape.txt
PYTHONIOENCODING=utf-8 python -X utf8 research/op-pi1-bb0nu-nme-isotope/phase2_nme_derivation.py 2>&1 | tee research/op-pi1-bb0nu-nme-isotope/phase2_nme_derivation.txt
PYTHONIOENCODING=utf-8 python -X utf8 research/op-pi1-bb0nu-nme-isotope/phase3_nme_predictions.py 2>&1 | tee research/op-pi1-bb0nu-nme-isotope/phase3_nme_predictions.txt
```

## Cross-references

- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
- [[../op-omicron1-sigmamnu-cosmo/Phase3_results.md]]
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]]
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
