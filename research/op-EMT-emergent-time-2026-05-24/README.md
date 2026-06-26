---
title: "op-EMT-emergent-time — TGP emergent-time formalism + (z, d_L) re-derivation"
type: research_program
status: DEFERRED
phase: pre-Phase 0 (skeleton)
folder_status: deferred
created_date: 2026-05-24
parent_motivation: "γ-5 §3.5 'time emerges from reconfiguration dynamics, rate = c' + sek08c g_eff_00 metric form — formal τ(N) mapping never derived in TGP"
authorization: "User 2026-05-24: 'trzeba stworzyć wszystkie 4, ustawić im odpowiedni status i kolejność'"
deferred_reason: "Scope is research PROGRAM (4-6+ sesji multi-cycle), not single cycle. Requires formalism development before any falsifier test possible."
independent_of: "γ-3, γ-3', γ-5, γ-7 F8 cycles (different mechanism category: temporal frame analysis, NOT classical kinematic)"
expected_duration: "Multi-cycle research program (4-6+ sesji minimum)"
---

# op-EMT-emergent-time — Emergent time formalism

## Status: DEFERRED

This is a **research program**, not a single cycle. Activation requires:
- (a) Completion of B (PSR), and consideration of A (LAM) and D (G-substrate) results first
- (b) Explicit decision by user to commit to multi-cycle program
- (c) Phase 0 design with sub-cycles (formalism / observable derivation / data comparison)

## Scope (when/if activated)

**Primary objective:** Develop explicit TGP-native emergent-time formalism, then test whether cosmological observables (SN Ia d_L(z), CMB acoustic peak position, BAO ruler) require corrections from TGP-frame vs classical-t frame.

**Hypothesis** (user puzzle I, recorded verbatim in [[../../meta/F8_FORENSIC_2026-05-24.md]] §7 Puzzle I):
- γ-5 §3.5 quotes "time emerges from reconfiguration dynamics, rate = c"
- sek08c g_eff_00 = -c₀²/ψ → proper time observer rate depends on local ψ
- ψ_eq evolves through cosmological epochs (PR-011 chain inflation → ... → today)
- If τ(N) mapping varies with cosmic epoch, observed d_L(z) inferred in fixed-clock frame may differ from TGP-native d_L(z)
- Part of "observed ä > 0" in SN Ia may be clock-drift artifact

## Why DEFERRED

This is **not** a single cycle test for two reasons:

1. **Formalism gap.** TGP currently has:
   - g_eff_00 = -c₀²/ψ (sek08c) — gives proper time locally
   - ψ_eq(t) evolution sketch (PR-011 chain inflation)
   - γ-5 §3.5 quote on emergence rate
   - **BUT** explicit τ(N) global cosmological mapping is NOT in concept paper
   
   Building this formalism is itself a multi-sesja research item.

2. **Observable derivation.** Even after formalism is built, deriving TGP-native d_L(z) requires:
   - Integrating ψ_eq(t) through 6 cosmological epochs (PR-011 chain)
   - Re-mapping light propagation in evolving g_eff
   - Computing observables in TGP frame
   - Comparing to SN Ia data
   
   Each of these is significant work.

**Estimated total:** 4-6+ sesji minimum, possibly split into 2-3 sub-cycles.

## Anti-Lakatos status

**This document is research-program declaration, NOT a pre-registered falsifier.**

If/when activated:
- Will require fresh Phase 0 with sub-cycle decomposition
- Cannot inherit F8 four-cycle verdict structure
- Must NOT be framed as F8 rescue
- Standalone failure mode: if formalism fails to produce coherent observable corrections, document it honestly

**Independence:** This cycle tests temporal-frame effects on cosmological observables. F8 (dark energy acceleration) is one possible domain of impact, but **not the primary motivation**. Primary motivation is **gap in TGP formalism**: emergent time framework is incomplete in concept paper.

## Inheritance (LEGITIMATE when activated)

- γ-5 §3.5 emergence dynamics quote (concept paper)
- sek08c g_eff_00 metric (concept paper)
- PR-011 chain inflation epoch sequence (LOCKED)
- sek08a thm:einstein-emergence (FRW dynamics base)
- γ-3'/γ-5 LOCKED ℓ_P, E_P, c calibration

## Provisional sub-cycles (if program activated)

Pre-registration would split into independent cycles:

1. **EMT-1 (formalism):** Define τ(N) mapping rigorously from sek08c + γ-5 §3.5
2. **EMT-2 (observable derivation):** Compute d_L(z), CMB peak, BAO in TGP frame
3. **EMT-3 (data comparison):** Compare to observations, test for clock-drift signatures

Each sub-cycle independent pre-registration. Cumulative outcome: whether emergent-time matters cosmologically.

## Files (when activated)

(NONE yet — deferred status)

Future:
- `Phase0_balance.md` — program scope declaration
- Sub-cycle subfolders: `EMT-1/`, `EMT-2/`, `EMT-3/`

## Status summary

| Field | Value |
|-------|-------|
| Status | DEFERRED |
| Reason | Multi-cycle research program; single Phase 0 insufficient |
| Activation requires | User explicit re-prioritization after B/A/D completion |
| F8 connection | Indirect (one possible application domain) |
| Anti-Lakatos | COMPLIANT (no inheritance from F8 FAILs) |
