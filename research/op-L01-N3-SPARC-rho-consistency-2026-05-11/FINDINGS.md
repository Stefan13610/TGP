---
title: "FINDINGS — op-L01-N3-SPARC-rho-consistency-2026-05-11 (eksportowalne wyniki)"
date: 2026-05-11
parent: "[[./README.md]]"
type: findings
total_findings: 10
total_sympy_pass: "8/8 (Phase 1)"
six_requirements: "6/6 RESOLVED"
classification: STRUCTURAL_DERIVED
tags:
  - findings
  - eksportowalne
  - L01-N3-closure
  - SPARC-consistency-verified
  - dust-limit-LOCK
---

# FINDINGS — eksportowalne wyniki cyklu

## §1 — Phase 1: SPARC ρ-consistency verification

| ID | Finding | Source |
|---|---|---|
| **F1.1** | Dust limit (p=0): T^μ_μ_dust = -ρ_rest·c² (analytic LOCK; standard GR Wald 1984, MTW 1973) | sympy T1 |
| **F1.2** | ρ_TGP ≡ -T^μ_μ_dust/c_0² = ρ_rest exact w c_0=c convention | sympy T2 |
| **F1.3** | Galactic stars v ~ 200 km/s: ρ_TGP/ρ_rest deviation **2.2·10⁻⁷** (= 2.2·10⁻⁵ %) — 6 OOM below 1% target | sympy T3 |
| **F1.4** | HI gas thermal v ~ 1 km/s (T~100 K): deviation **5.6·10⁻¹²** — utterly negligible | sympy T4 |
| **F1.5** | **R1 closed strukturalnie:** TGP-emergent DM jest *gravitational* (g_eff[Φ̄] modification z multi-source `{Φ_i}` interaction), NIE *matter* sektor; NIE additive z ρ_baryon (S05 single-Φ enforced) | sympy T5+T7 |
| **F1.6** | SPARC framework (galaxy_scaling cycles gs10-gs61) używa ρ_baryon only — strukturalnie consistent z TGP framework; chi²_red competitive z MOND simple | sympy T6 |
| **F1.7** | S05 single-Φ preservation: matter source ρ_baryon → ρ_TGP via dust limit; gravitational z g_eff[{Φ_i}]; brak second fundamental DM field | sympy T7 |
| **F1.8** | **R2 honestly documented:** SPARC scope = galactic-disk regime (R~1-50 kpc, NR valid); near-SMBH ISCO regime (v~c/2) wymaga full GR/TGP-PPN (outside L01-N3 cycle scope) | sympy T8 |
| **F1.9** | **ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0²** verified do **<10⁻⁶ precision** (target 1%, achieved 6 OOM below) | §1 + sympy T3 |
| **F1.10** | Cluster-scale ~35% mass deficit (separate issue, requires ~2 eV sterile ν per gs13-gs55) — outside N3 scope | §3.2 |

## §2 — Risk register (final status)

| Risk | Status | Source |
|---|---|---|
| **R1** (double-counting vs separate ρ_DM) | **closed strukturalnie** | TGP-emergent DM gravitational (NOT matter); S05 enforced |
| **R2** (galactic-center relativistic limit) | **honestly documented** | SPARC scope clarified (galactic-disk regime; near-SMBH outside scope) |

## §3 — Sympy LOCK summary

| Phase | Tests | PASS |
|---|---|---|
| Phase 1 | 8 | 8/8 |
| **Total** | **8** | **8/8** |

## §4 — Predictions for PREDICTIONS_REGISTRY (1 entry)

| Entry ID | Prediction | Test | Status |
|---|---|---|---|
| **M911-SPARC-rho-consistency** | ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0² (HI + stars + bulge); galactic dust limit (v²/c² ~ 10⁻⁷) gives correction ≪ 10⁻⁶; brak separate ρ_DM matter component (TGP-emergent DM = gravitational mechanism z g_eff[Φ̄]) | SPARC database 175 galaxies (Lelli+2016); galaxy_scaling cycles gs10-gs61 chi²_red competitive z MOND simple | TESTED-PASS (verification cycle 2026-05-11; 8/8 sympy; deviation < 10⁻⁶) |

## §5 — Cross-cycle synthesis findings

| ID | Finding | Source |
|---|---|---|
| **F5.1** | **L01 N3 closure source** — cycle daje compact verification N3 z 8/8 sympy PASS; closes residual cosmetic L01 cleanup item | Phase_FINAL §11 |
| **F5.2** | **Cross-cycle convergence updated diagnostic** — siedem niezależnych diagnoz (L01 §3.2, τ.3 §2, ψ.1 §3, Q2, op-L01-N1, op-L01-N2, **op-L01-N3 this cycle**) zbieżnych na separable sector structure jako *strukturalna własność TGP*, NIE post-hoc tuning | Phase_FINAL §5 |
| **F5.3** | **Q2 F1 trzy niezależne konstruktywne weryfikacje** (post-2026-05-11): N1 operator-class disjointness + N2 vacuum-level + N3 gravitational-vs-matter sector separation | F5.2 |

---

**Findings exported.** 10 distinct findings + 2 risks + 8/8 sympy + 6/6 P-requirements +
1 PREDICTIONS_REGISTRY entry.
