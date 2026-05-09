---
title: "Phase 6 setup — SU(2) cross-consistency (N11), CRITICAL"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-setup
phase: 6
status: 🟢 OPEN
needs: [N11]
predecessor: "[[./Phase5_results.md]] (10/10 PASS, m_inertial=m_grav AUTOMATIC)"
related_cycle: "[[../op-SPIN-SU2-substrate-derivation-2026-05-08/]] (CLOSED, 47/47 PASS)"
---

# Phase 6 setup — SU(2) cross-consistency (N11)

## §0 — Goal

Resolve N11 from [[./NEEDS.md]] (CRITICAL programmatic):
> "Mechanizm interakcji generujący g_eff (gradient cross-terms) musi być
> zgodny z mechanizmem generującym SU(2) z SPIN cycle (3 paths)."

**Strategic context:** SPIN-SU2 cycle ([[../op-SPIN-SU2-substrate-derivation-2026-05-08/]])
closed STRUCTURAL DERIVED w 2026-05-09 z 47/47 sympy PASS. SU(2) emergence
via THREE INDEPENDENT PATHS:

- **Path A (N17, N18)**: Dynamic-equilibrium 2-state bifurcation (zanik/ekspansja
  z V(φ) saddle → Pauli matrix algebra)
- **Path B (N21)**: M9.1'' horizon multipole l=1 dipole → SO(3) → SU(2) double cover
- **Path C (N19)**: External embedding (lean direction θ,φ → induced SO(3) action)

## §1 — Critical question

**Czy emergent-metric framework (Phase 1-5) jest STRUKTURALNIE ZGODNY z SPIN
cycle SU(2) derivation?** Konkretnie:

1. Czy σ-coupling C(ψ) (Phase 3 lock) **psuje którąkolwiek z 3 SU(2) paths**?
2. Czy 3PN parameter changes (Path 1 Phase 4: ξ_3 = (32-a_3)/32) preservują
   V(φ) saddle structure (Path A) i M9.1'' horizon (Path B)?
3. Czy mechanizm "dynamic equilibrium" w SPIN cycle (interakcja z Φ̄)
   **JEST tym samym** mechanizmem co σ-cross-terms w emergent-metric?
4. Czy c_0 (σ-coupling leading coefficient) jest **derivable** z SPIN
   cycle dynamic-equilibrium balance (analog SPIN N16 quintic equation)?

## §2 — Hipoteza H6.1: structural unification

**H6.1:** TGP ma JEDNĄ ZASADĘ generowania tensor structure z interactions:
- **Level 2 (metryka)**: g_eff^μν = G[{Φ_i}, σ_ab, Φ̄] z gradient cross-terms
- **Level 3 (spin)**: SU(2) z dynamic-equilibrium soliton-Φ̄ interaction

Both share:
- Φ̄ background (vacuum reference field)
- Localized δΦ source perturbations
- Self-energy + interaction balance
- Tensor structure emergent z interakcji

**Falsifier:** jeśli emergent-metric C(ψ) jest niespójne z SPIN cycle SU(2)
derivation, H6.1 fails → STRUCTURAL_NO_GO honest.

**Confirmation:** jeśli c_0 jest **derivable** z dynamic-equilibrium constraint
(analog N16 balance), albo if 3 SU(2) paths są INDEPENDENT od C(ψ),
H6.1 confirmed → cycle CLOSES STRUCTURAL DERIVED FULL.

## §3 — Path-by-path consistency analysis

### §3.1 — Path A (dynamic-equilibrium bifurcation): c_0-independent?

SPIN N17 mechanism: V(φ) = γ[φ³/3 - φ⁴/4] saddle structure → 2-state
bifurcation. **V(φ) is matter-sector potential (V_orig in dual-V framework).**

Emergent-metric C(ψ) modifies g_eff^ij (gravity sector). In dual-V framework
(per [[../op-dual-V-structure-clarification-2026-05-09/]]):
- V_grav (gravity sector) ≠ V_matter (matter sector)
- C(ψ) is gravity-sector coupling
- N17 bifurcation uses V_matter ⟹ **independent of C(ψ)**

**Hypothesis: Path A is c_0-independent (fully).**

Sympy verification: derive V_matter saddle structure given C(ψ) modification of
g_eff. If V_matter unchanged (per dual-V framework), bifurcation preserved.

### §3.2 — Path B (M9.1'' horizon multipole): SENSITIVE to ansatz change

SPIN N21 mechanism: M9.1'' canonical metric f(ψ) = (4-3ψ)/ψ has horizon at
ψ=4/3. l=1 dipole multipole structure of horizon → SO(3) → SU(2).

**PROBLEM**: M9.1'' specific FALSIFIED 5σ. Phase 3-4 emergent-metric uses
parametric family, NOT M9.1'' canonical.

Phase 4 Path 2 (keep M9.1'' parameters, add c_0): horizon structure preserved.
Phase 4 Path 1 (change ξ_3): horizon may shift or disappear.

**Hypothesis: Path B compatible only with Phase 4 Path 2 (keep M9.1'' params).**

This GIVES PREFERENCE for Path 2 over Path 1 — strukturalny argument za
σ-coupling jako mechanizm post-falsyfikacji recovery, NIE 3PN parameter change.

### §3.3 — Path C (external embedding): c_0-independent?

SPIN N19 mechanism: lean direction (θ, φ) → SO(3) external action → induced SU(2).

External embedding is geometric structure on R³ — NIEZALEŻNE od g_eff details
at vacuum. **Hypothesis: Path C is c_0-independent.**

## §4 — c_0 derivation from dynamic-equilibrium balance?

### §4.1 — SPIN N16 balance (reference)

```
E(λ) = T/λ + U/λ³ + S·λ²    (kinetic + potential + skin)
dE/dλ = 0 ⟹ 2S·λ⁵ - T·λ² - 3U = 0    (quintic)
λ_eq* ≈ 1.169 (sympy nroots dla T=U=S=1)
```

S = "skin coupling" parametr (interaction with Φ̄ background).

### §4.2 — Analog dla emergent-metric c_0?

Hypothesis: **c_0 emerges from analogous balance equation in 2-source binary.**

Conjecture: σ-cross-coupling at binary equilibrium scales as:
```
c_0 · κ_σ ≈ specific function of (E_self, E_inter, λ_eq)
```

If yes: c_0·κ_σ = 4/3 (Phase 4 Path 2 GR-match) should be derivable
z analog quintic equation w 2-source case.

**Realization:** explicit derivation requires:
1. Setup 2-source soliton energy budget (E_self_1 + E_self_2 + E_inter_12)
2. Variational extremum condition
3. Identification σ-coupling z structural parameter

This is **multi-session sympy work** — Phase 6 will scope it.

## §5 — Phase 6 strategy

### §5.1 — Tractable subset

Phase 6 will **NOT** attempt full c_0 derivation in single session. Scope:

1. **Path A consistency** (c_0-independent claim): sympy verification że dual-V
   framework preserves V(φ) bifurcation under c_0 modification.

2. **Path B sensitivity**: structural analysis identifying which Phase 4 family
   region preserves horizon structure.

3. **Path C consistency** (c_0-independent claim): sympy verification że external
   embedding is independent of c_0.

4. **c_0 derivability roadmap**: identify the explicit calculation needed dla
   c_0 first-principles in subsequent sessions.

### §5.2 — Decision tree

| Phase 6 outcome | Cycle classification |
|---|---|
| All 3 paths c_0-independent (Path A, C) + Path B preserved | STRUCTURAL DERIVED FULL |
| Path A + C c_0-independent, Path B sensitive | STRUCTURAL DERIVED (3 paths still work, 2 robust to c_0) |
| Some path BREAKS under c_0 modification | STRUCTURAL_CONDITIONAL — re-examine emergent-metric ansatz |
| All paths break | STRUCTURAL_NO_GO — radical rethink needed |

**Default expectation:** Path A and C robust; Path B sensitive but Phase 4 Path 2
preserves it. Cycle should classify STRUCTURAL DERIVED.

### §5.3 — c_0 derivation deferral (honest)

Even with full Phase 6 PASS, c_0 numerical pinning likely deferred to:
- Dedicated cycle: op-c0-derivation-from-substrate
- Or: H_Γ coarse-graining cycle
- Estimated: 5-10 sesji multi-session

This is acknowledged in handoff. Phase 6 SUCCESS = structural compatibility
proven, NOT canonical c_0 number.

## §6 — Phase 6 deliverables

- Phase6_sympy.py
- Phase6_sympy.txt
- Phase6_results.md

## §7 — Gate criteria

| # | Criterion |
|---|---|
| G1 | Path A V(φ) bifurcation preserved under c_0 modification (dual-V) |
| G2 | Path B sensitivity analysis (which Phase 4 family region preserves horizon) |
| G3 | Path C external embedding c_0-independent |
| G4 | At least 2 of 3 paths preserved → SU(2) emergence robust |
| G5 | c_0 derivation roadmap documented (deferred path forward) |
| G6 | Cross-consistency H6.1 verdict (CONFIRMED / SENSITIVE / NO_GO) |

**Pass:** ≥4/6.

## §8 — Anti-pattern compliance

- NIE post-hoc favored Phase 4 Path 2 (Path B preserved) — must justify ex ante
- NIE c_0 number fitting (acknowledge multi-session derivation needed)
- Honest reporting: if Path B is uniquely sensitive, document explicitly

## §9 — Cross-references

- [[./README.md]] — cycle overview
- [[./Phase5_results.md]] — predecessor (10/10 PASS, EP automatic)
- [[../op-SPIN-SU2-substrate-derivation-2026-05-08/Phase6_absolute_binding.md]] — SPIN closure
- [[../op-SPIN-SU2-substrate-derivation-2026-05-08/Phase1_N16_results.md]] — dynamic-equilibrium
- [[../op-SPIN-SU2-substrate-derivation-2026-05-08/Phase1_N17_results.md]] — bifurcation
- [[../op-dual-V-structure-clarification-2026-05-09/]] — dual-V framework
