---
title: "π.1.Phase2 setup — derivation hardening (7 sub-tests)"
date: 2026-04-30
cycle: π.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - pi1
  - phase2
  - NME
  - derivation
---

# π.1.Phase2 — derivation hardening (7 sub-tests)

> **Cel:** Build T_{1/2}(iso, Form, NME) closed-form 24-cell matrix;
> derive cross-isotope T_{1/2} ratios that cancel m_ββ; verify universal
> 4.21× Form A/B half-life factor; estimate TGP-native NME via
> B²·Z·A^{1/3} cascade scaling; FALSIFY 6 alt fits; LOCK 4 promotions.

## Sub-tests

### P2.1 — T_{1/2}(iso, Form, NME) closed form

24-cell matrix; verify Form B/Form A T_{1/2} ratio = 4.21 universally.

### P2.2 — Cross-isotope T_{1/2} ratios cancel m_ββ

R(iso1, iso2) = [G(iso2)·M²(iso2)] / [G(iso1)·M²(iso1)]

Form-independent; pure NME×PSF systematics test.

### P2.3 — Universal Form A/B 4.21× factor

Single number across all (iso, NME) cells = (m_ββ_B / m_ββ_A)² = 2.051².

### P2.4 — TGP-native NME estimator via closure 1/A^{1/3}

Closure approximation (NME ∝ <r²>^{−1/2}, r ∝ A^{1/3}):

M_TGP(iso) = M_ref · (A_ref / A)^{1/3}

with A_ref=76, M_ref=4.30 (Ge-76 lit mean). B²-cascade enters m_ββ
Form A/B (already established in ν.1). Pass criterion: TGP-native NME
within factor 2 of QRPA/IBM-2/NSM/EDF mean per isotope.

### P2.5 — NME method best-match vs TGP

Find which method (QRPA/IBM-2/NSM/EDF) minimizes χ² vs TGP-native B²·Z·A^{1/3}.

### P2.6 — 6 alt fits FALSIFIED

| # | Alt fit | Falsification |
|---|---------|---------------|
| 1 | g_A quenched 0.7 | T_{1/2} drops 4× → still above bound but breaks ratio universality |
| 2 | Single-method ALL-NSM | predicts Ge/Xe ratio 2.22 vs QRPA 2.86 |
| 3 | Single-method ALL-IBM-2 | predicts Te/Xe 0.66 vs NSM 0.92 |
| 4 | IH m_lightest > 50 meV | m_ββ ≥ 18 meV → T_{1/2} 100× lower |
| 5 | m_ββ = 5 meV (LEGEND target) | 4× higher than Form B |
| 6 | TGP m_ββ_A boost ζ_TGP | m_ββ_A → 8.5 meV (5.4×) — ν.1 falsified retroactively |

### P2.7 — 4-way LOCK promotions

| # | Item | Status |
|---|------|--------|
| 1 | m_ββ_A = 1.584 meV LOCKED | ν.1 + π.1 cross-check |
| 2 | m_ββ_B = 3.249 meV LOCKED | ν.1 + π.1 cross-check |
| 3 | T_{1/2,B}/T_{1/2,A} = 4.21 universal | LOCKED across iso/NME |
| 4 | Ge-76 NME-cleanest isotope | σ_NME=0.465 vs Te 0.837 |

## PASS bramka

- ≥6/7 PASS → π.1.Phase3 OK
- 7/7 PASS → FULL CASCADE

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-pi1-bb0nu-nme-isotope/phase2_nme_derivation.py 2>&1 | tee research/op-pi1-bb0nu-nme-isotope/phase2_nme_derivation.txt
```

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
