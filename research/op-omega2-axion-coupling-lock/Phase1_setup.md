---
title: "ω.2.Phase1 setup — anomaly-based g_bare derivation + alt falsification"
date: 2026-05-01
cycle: ω.2.Phase1
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - omega2
  - phase1
  - anomaly
  - chirality
  - setup
---

# ω.2.Phase1 setup

**Score gate:** ≥4/5 PASS = Phase 1 forward.

## Sub-tests (5)

- **W2.1.1** Triangle anomaly E_TGP from chirality B² (sympy-exact)
- **W2.1.2** Bare-coupling structural form g_bare = η_chir = 1−C6 = 19/24
- **W2.1.3** IR effective g_eff = α_em·E_TGP/(2π) — drift vs LOCK candidates
- **W2.1.4** Dimensional consistency + scale-symmetry preservation
- **W2.1.5** Alt-form falsification (5 alts: α_em / 1/(2π) / κ_TGP / C6 / 19/24)

## Inputs from prior cycles (LOCKED)

```
B²_up   = 13/4   (θ.1 — Dirac quarks + QCD)
B²_down = 61/25  (θ.1 derived — B²_up − 81/100 sympy)
B²_lep  = 2      (Dirac leptons)
B²_ν    = 1      (Majorana, halved — ν.1)
K_up    = 7/8    (K-taxonomy (2+B²_up)/(2·N_gen))
K_lep   = 2/3    (K-taxonomy)
K_ν     = 1/2    (K-taxonomy)
C6      = K_up − K_lep = 5/24  (ρ.1 LOCKED)
N_gen   = 3      (φ.1 triple-locked)
N_c     = 3      (color factor)
Q_u     = +2/3   (up-type EM charge)
Q_d     = −1/3   (down-type)
Q_l     = −1     (charged lepton)
α_em    = 1/137.036  (CODATA, anchored)
```

## Cross-references

- [[program.md]]
- [[../op-omega1-substrate-em-coupling/Phase3_results.md]]
- [[../op-rho1-71Ge-cross-section/Phase3_results.md]] — C6 = 5/24 LOCKED
- [[../op-theta-quark-koide/Phase3_results.md]] — B²_up, B²_down
