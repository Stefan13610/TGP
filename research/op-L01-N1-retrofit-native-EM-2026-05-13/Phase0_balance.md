---
title: "Phase 0 — Balance sheet (N1-EM retrofit)"
date: 2026-05-13
parent: "[[./README.md]]"
phase: 0
---

# Phase 0 — Balance sheet

## §1 — Predecessor disposition

| Cycle | Status | Rola |
|---|---|---|
| `op-L01-N1-EM-trace-anomaly-TGP-2026-05-11` | CLOSED-DOWNGRADED **D (ALGEBRAIC_MIMICRY)** | Reference scope; algebraic facts preserved; 11/16 TAUTOLOGY+HARDCODED testów ZASTĄPIONE first-principles |
| `op-L01-rho-stress-energy-bridge-2026-05-04` | CLOSED-DERIVED | ax:metric-coupling foundation |
| `op-emergent-metric-from-interaction-2026-05-09` | 57/57 PASS | g_eff[{Φ_i}] background |
| `op-L01-N3-retrofit-native-SPARC-2026-05-13` | CLOSED-RESOLVED A− | Sister retrofit; template pattern z 11/11 PASS i 9 FP |

## §2 — 6/6 gate

| # | Requirement | Sympy test |
|---|---|---|
| P1 | 1-loop QED β-function w g_eff curved background | T2 (FP) |
| P2 | T^μ_μ_EM = (β/(2g)) F² + Riegert | T1 + T3 (FP) |
| P3 | Disjointness dim-4 vs dim-6 (Theorem 2.1) | T4 (FP) |
| P4 | GW170817 c_GW = c_EM | T6 (FP) + T7 (LIT) |
| P5 | MICROSCOPE η universality | T9 (FP) + T8 (LIT) |
| P6 | S05 preserved | T5 (FP) + T10 (DEC) |

## §3 — Risk register

| Risk | Mitigation |
|---|---|
| R1 (β as literature) | T2 explicit symbolic derivation via dimensional reg + Pattern 2.1 linearization |
| R2 (M9.1'' contamination) | Use generic {A,B,C} ansatz z emergent-metric Phase 1 |
| R3 (declarative Theorem 2.1) | T4 explicit operator dimension counting + structural disjointness sympy |
| R4 (non-pert magnetar OUT) | T11 scope declaration |

## §4 — Scope

**IN:** 1-loop QED on g_eff[{Φ_i}] curved background; GW170817 dispersion bound; MICROSCOPE η; WEP.

**OUT:** Non-perturbative magnetar B >> B_QED; ψ.1.v3 dim-6 EFT (separate cycle); CMB photon
trace anomaly (separate `op-L01-N2-retrofit` for QCD trace which couples z cosmology).

## §5 — Estymata

- Phase 1: 1 sesja (sympy first-principles + symbolic Riegert)
- Phase 2: 0.5 sesji (observational bound check)
- Phase FINAL: 0.5 sesji

**Total:** ~2 sesji (skrócony estymata z RESEARCH_RESTART §1.3 dzięki SPARC pattern replication).

## §6 — Phase 0 sign-off

Gate OPEN. Next: `Phase1_sympy.py`.
