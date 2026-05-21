---
title: "Phase 0 — Balance sheet + 6/6 gate (retrofit native SPARC)"
date: 2026-05-13
parent: "[[./README.md]]"
type: phase-balance
phase: 0
---

# Phase 0 — Balance sheet

## §1 — Predecessor disposition

| Cycle | Status | Rola w tym retrofit |
|---|---|---|
| `op-L01-N3-SPARC-rho-consistency-2026-05-11` | CLOSED-DOWNGRADED **D (ALGEBRAIC_MIMICRY)** | Reference dla scope; algebraic facts (dimensional consistency, v²/c² OOM) preserved; sympy 5/8 hardcoded True ZASTĄPIONE first-principles symbolic |
| `op-L01-rho-stress-energy-bridge-2026-05-04/formal_definition.md` | CLOSED-DERIVED 2026-05-10 | Foundation: ax:metric-coupling derivation `ρ ≡ -T^μ_μ/c_0²` |
| `op-emergent-metric-from-interaction-2026-05-09` | CLOSED-RESOLVED 57/57 PASS | g_eff[{Φ_i}] background dla galactic rotation |
| `op-cluster-mass-deficit-resolution-2026-05-11` | CLOSED-NULL EARLY_HALT_HONEST | **Scope clarification:** cluster ~Mpc scale OUTSIDE niniejszego galactic-disk regime |

## §2 — 6/6 gate (six P-requirements)

| # | Requirement | Phase 1 verification path |
|---|---|---|
| P1 | Perfect fluid T_μν decomposition z 4-velocity u^μ symbolic | sympy T1 |
| P2 | ρ ≡ -T^μ_μ/c_0² consistent dust limit (p=0) → ρ_TGP = ρ_rest EXACT | sympy T2+T3 |
| P3 | Non-relativistic correction (1 - v²/2c²) explicit Lorentz boost | sympy T4+T5 |
| P4 | SPARC ρ_baryon mapping consistent z ρ_TGP | sympy T6+T7 (structural + LIT-anchored) |
| P5 | No double-counting vs TGP-emergent DM (S05 single-Φ via g_eff[Φ̄] gravitational, ρ_baryon matter — separable sectors) | structural argument (T9 declarative) |
| P6 | S05 single-Φ axiom preserved bezwarunkowo | sympy T8 (ax:metric-coupling consistency) + T9 declarative |

**Gate result:** 6/6 P-requirements MAPPABLE do Phase 1 sympy plan (NIE pre-declared
PASS — verification w Phase 1).

## §3 — Risk register

| Risk | Mitigation |
|---|---|
| R1 (Lorentz boost expansion higher-order) | Use sympy.series() z explicit order=O(v⁴); verify next-to-leading is positive O(v⁴/c⁴) |
| R2 (Bullet Cluster scope) | Explicit scope: galactic-disk regime SPARC; cluster-scale POZA |
| R3 (ρ_baryon column independence) | Structural argument w T9 declarative + verify w SPARC database |
| R4 (MOND chi²_red benchmark) | Use Lelli+2017 documented chi²_red ≈ 2.0 jako external; TGP must reach competitive |

## §4 — Scope (in / out)

**IN scope (this cycle):**
- First-principles ρ ≡ -T^μ_μ/c_0² derivation symbolic z ax:metric-coupling
- Dust limit (p=0) → ρ_TGP = ρ_rest EXACT
- Non-relativistic Taylor expansion w v²/c²
- SPARC galactic-disk regime (v ~ 200 km/s stars; v ~ 1 km/s HI)
- L1 native v_rot(R) chi²_red specification dla SPARC

**OUT scope (deferred):**
- Cluster-scale (Mpc) — separate cycle per cluster cycle EARLY_HALT_HONEST
- Bullet Cluster lensing-vs-X-ray offset — separate cycle
- Full numerical SPARC fitting (175 galaxies) — separate cycle (galaxy_scaling)
- MOND comparison — L2 optional last stage
- Near-SMBH relativistic (Sgr A*) — outside SPARC regime

## §5 — Estymata Phase 1

- Phase 1 sympy + results: 1 sesja (this session)
- Phase 2 SPARC mapping + structural verification: 1 sesja optional
- Phase FINAL closure + adversarial audit: 1 sesja

**Total: 2-3 sesji per RESEARCH_RESTART §1.3 estymata.**

## §6 — Phase 0 sign-off

Phase 0 balance complete. Phase 1 sympy gate OPEN.

**Next:** uruchom `Phase1_sympy.py`.
