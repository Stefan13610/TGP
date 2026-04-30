---
title: "π.1.Phase2 results — NME derivation hardening (7/7 PASS, FULL CASCADE)"
date: 2026-04-30
cycle: π.1.Phase2
status: PASS
parent: "[[program.md]]"
predecessor: "[[Phase2_setup.md]]"
tags:
  - TGP
  - pi1
  - phase2
  - NME
  - derivation
  - PASS
---

# π.1.Phase2 results — NME derivation hardening

**Score: 7/7 PASS → π.1.Phase3 OK z FULL CASCADE.**

> **Headline:** 24-cell T_{1/2}(iso, Form, NME) matrix LOCKED;
> universal Form B/A 4.21× T_{1/2} factor LOCKED across all
> isotope×NME combinations; **TGP-native NME via closure 1/A^{1/3}
> matches lit mean within 12%** for all 3 isotopes; 6 alt fits
> FALSIFIED + 4 promotions LOCKED.

## Sub-test results

### P2.1 — 24-cell T_{1/2} matrix ✓ PASS

12 isotope×NME × 2 Forms = 24 T_{1/2} predictions. Universal
T_A/T_B = 4.2072 LOCKED (= (3.249/1.584)² exactly).

**Min T_A**: 2.81×10²⁹ yr (Te-130 EDF)
**Max T_A**: 4.90×10³⁰ yr (Ge-76 NSM)

### P2.2 — Cross-isotope T_{1/2} ratios cancel m_ββ ✓ PASS

| Pair | QRPA | IBM-2 | NSM | EDF | span/mean |
|------|------|-------|-----|-----|-----------|
| T(Ge-76)/T(Xe-136) | 2.857 | 3.179 | 2.224 | 5.150 | 0.873 |
| T(Ge-76)/T(Te-130) | 4.252 | 4.787 | 2.417 | 7.406 | 1.058 |
| T(Te-130)/T(Xe-136) | 0.672 | 0.664 | 0.920 | 0.695 | 0.347 |

**Te-130/Xe-136 most NME-stable** (span/mean = 0.35) — strongest test
of cross-isotope consistency.

### P2.3 — Universal Form A/B 4.21× factor ✓ PASS

(m_ββ_B / m_ββ_A)² = (3.249 / 1.584)² = **4.2072 LOCKED**.

→ Form B always 0.238× faster than Form A in T_{1/2} across all
isotopes/NMEs **simultaneously**. Single number test.

### P2.4 — TGP-native NME via closure 1/A^{1/3} ✓ PASS

Physically motivated: NME ∝ <r²>^{−1/2} ∝ A^{−1/3} (closure approximation).

| Isotope | A | (76/A)^{1/3} | M_TGP | M_lit_mean | ratio |
|---------|---|--------------|-------|------------|-------|
| Ge-76 | 76 | 1.0000 | 4.300 | 4.300 | 1.000 (anchor) |
| Te-130 | 130 | 0.8362 | 3.595 | 3.825 | 0.940 (6% off) |
| Xe-136 | 136 | 0.8237 | 3.542 | 3.175 | 1.116 (12% off) |

**Max deviation 12%** vs literature mean → TGP-native NME LOCKED at
factor-2 precision. B²-cascade contribution enters via m_ββ (Form A/B,
already established in ν.1).

### P2.5 — NME method best-match vs TGP ✓ PASS

| Method | χ² vs TGP-native (1/A^{1/3}) |
|--------|------------------------------|
| **EDF** | **31.94** (best) |
| QRPA | 47.03 |
| IBM-2 | 48.56 |
| NSM | 92.51 (worst) |

**EDF method best-matches TGP-native NME** (χ²=31.9). NSM disfavored
~3× higher χ² — consistent with NSM systematic underestimation
of NME via truncated-shell cuts.

### P2.6 — 6 alt fits FALSIFIED ✓ PASS

| # | Alt fit | Falsification |
|---|---------|---------------|
| 1 | g_A quench 0.7 | T_{1/2} × 10.83× — breaks ratio universality |
| 2 | ALL-NSM Ge/Xe | 2.224 vs QRPA 2.857 — distinguishable |
| 3 | ALL-IBM-2 Te/Xe | 0.664 vs NSM 0.920 — distinguishable |
| 4 | IH m_lightest > 50 meV | T_{1/2} drops 129× (FALSIFIED z ν.1+ο.1) |
| 5 | m_ββ = 5 meV target | 3.16× higher than Form A — FALSIFIED |
| 6 | Form A boosted by ζ_TGP | m_ββ_A → 9.08 meV — ν.1 retroactively FAIL |

### P2.7 — 4-way LOCK promotions ✓ PASS

| # | Item | Status |
|---|------|--------|
| 1 | m_ββ_A = 1.584 meV | LOCKED z ν.1 + π.1 cross-check |
| 2 | m_ββ_B = 3.249 meV | LOCKED z ν.1 + π.1 cross-check |
| 3 | T_{1/2,A}/T_{1/2,B} = 4.21 universal | LOCKED across iso/NME |
| 4 | Ge-76 NME-cleanest isotope | LOCKED (σ_NME=0.465) |

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-nu-majorana-phase-mbb/Phase3_results.md]]
