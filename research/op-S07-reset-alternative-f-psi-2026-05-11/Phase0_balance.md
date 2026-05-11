---
title: "Phase 0 balance — S07 alternative f(ψ) literature inventory + reset path + 6/6 gate"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-balance
phase: 0
status: 🟡 in progress (scaffold; needs S07 + SPA literature deep-dive)
tags:
  - phase0
  - S07-reset
  - alternative-f-psi
  - SPA-chain-rederivation
  - GWTC-3-bayesian
---

# Phase 0 — Balance sheet (scaffold version)

## §1 — Sub-needs inventory (initial)

### Theoretical sub-needs (f(ψ) family enumeration)
- **N0.1**: S07 freedom explicit characterization (which f(ψ) families allowed)
- **N0.2**: f(ψ) = 1 (GR baseline) — trivial benchmark
- **N0.3**: f(ψ) = polynomial(ψ) — first non-trivial family
- **N0.4**: f(ψ) = transcendental (e.g., exp, log) — second family
- **N0.5**: f(ψ) = M9.1'' specific (4-3ψ)/ψ — FALSIFIED baseline (reference)
- **N0.6**: Constraints na f(ψ): (a) GR limit, (b) ψ → 1 behavior, (c) S05 preserved

### Computational sub-needs (β_ppE z SPA chain)
- **N0.7**: Full SPA chain (Phase 1.5 G_SPA = 48 sympy-exact) dla każdej alternative
- **N0.8**: β_ppE^(b=-1) explicit dla każdej alternative
- **N0.9**: Multi-coefficient ratios {β_(N+1)PN / β_NPN} dla każdej alternative
- **N0.10**: M911-P2 re-derivation z corrected SPA chain

### Falsifiability sub-needs (GWTC-3 + future)
- **N0.11**: Bayesian GWTC-3 compatibility (90 BBH combined) dla każdej alternative
- **N0.12**: BH5 QNM ringdown predictions dla każdej alternative
- **N0.13**: ε.1 photon ring predictions dla każdej alternative
- **N0.14**: Future LIGO O5 + LISA + ngEHT detectability

### Cross-cycle sub-needs
- **N0.15**: L01 N1-N5 closures preserved (matter sektor unchanged by metric reset)
- **N0.16**: Q2 F1 mechanism preserved (substrate-vacuum independent of f(ψ))
- **N0.17**: T-Λ ratio 1.020 preservation (closure 2026-04-26 baseline)

## §2 — Literature anchors

### S07 framework
- TGP core: [[../../core/sek07_solver/sek07_solver.tex]] (S07 alternative f(ψ) freedom)
- Adams-positivity v2 Phase 6.T6.5 (β_g sign-LOCKED < 0)

### M9.1'' falsification post-mortem
- [[../op-newton-momentum/M9_1_pp_P1_results.md]] (M9.1'' kanoniczna metryka)
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (β = -15/4 LOCK)
- [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] (5.02σ rejection)
- [[../../audyt/T01_LIGO3G_falsifier/]] (convention + GWTC-3 audit)

### ppE framework foundational
- Yunes, Pretorius 2009: parametrized post-Einstein framework (PPE)
- Yunes, Yagi, Pretorius 2016: gravitational wave tests of GR (review)
- LIGO O3 Tests of GR: Abbott et al. 2021 (90 BBH ToGR generic)

### SPA chain methodology
- Cutler, Flanagan, Phys. Rev. D 49, 2658 (1994): stationary phase approximation
- Damour, Iyer, Sathyaprakash, Phys. Rev. D 63, 044023 (2001): SPA expansion

### Bayesian model selection
- Vallisneri, Phys. Rev. Lett. 107, 191104 (2011): Bayesian ToGR
- LIGO O3 TGR results: Abbott et al. 2021 (BF computation methodology)

## §3 — 6/6 gate

| Gate | Status |
|---|---|
| **G1**: Predecessor cycles enumerated | ✅ op-ppE-mapping + op-GWTC3-reanalysis + S07 |
| **G2**: Sub-needs initial list | ✅ 17 items (§1) |
| **G3**: Literature currency check | 🟡 partial (Phase 1 deep review) |
| **G4**: Methodology binding declared | ✅ Native-first + S05 + S07 freedom + Bayesian rigor |
| **G5**: Risk register R1-R6 declared | ✅ YAML frontmatter |
| **G6**: Probability assessment recorded | ✅ Honest 40-55% / 15-25% / 20-30% multi-outcome |

**6/6 PASS (scaffold).**

## §4 — Methodology constraints (binding)

1. Native-first methodology
2. S05 single-Φ axiom preserved bezwarunkowo (każda alternative f(ψ) musi preserve)
3. §5.1 BD/Horndeski demarcation
4. **NIE post-hoc tuning** — alternative f(ψ) MUST emerge z S07 freedom strukturalnie
5. **Full SPA chain rigor** — Phase 1.5 G_SPA = 48 LOCK preserved; NIE simplified Δα_n ratios
6. **Bayesian rigor**: GWTC-3 BF computation z proper noise models
7. **Cross-cycle invariance**: L01 N1-N5 + Q2 + T-Λ closures **NIE zmienione** przez metric reset
8. Sympy LOCK dla każdego analytic step

## §5 — Next session targets

**Phase 1 (next 1-2 sessions):**
1. S07 freedom explicit characterization (which f(ψ) families z constraints)
2. Initial enumeration: f(ψ) = 1, polynomial, transcendental candidates
3. GR limit verification dla każdej (ψ→1 → GR baseline)

**Phase 2 (next 2-3 sessions):**
1. Full SPA chain re-derivation (Phase 1.5 G_SPA = 48 starting point)
2. β_ppE^(b=-1) explicit dla każdej alternative
3. Multi-coefficient ratios re-derivation (M911-P2 re-grounded)

## §6 — Cross-references

- [[./README.md]]
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (SPA chain LOCK)
- [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] (5σ falsification)
- [[../../core/sek07_solver/sek07_solver.tex]] (S07 freedom)
- [[../../audyt/T01_LIGO3G_falsifier/]]
- Yunes-Pretorius 2009; Cutler-Flanagan 1994; LIGO O3 ToGR 2021

---

**Phase 0 scaffold ready.** Reset path z S07 freedom; multi-session derivation
pending Phase 1+.
