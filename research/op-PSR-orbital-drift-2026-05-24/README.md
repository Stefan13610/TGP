---
title: "op-PSR-orbital-drift — TGP native O(U³) orbital-drift falsifier (NS surface)"
type: research_cycle
status: CLOSED-RESOLVED
phase: FINAL (4-phase single-session sprint complete)
folder_status: closed-resolved
created_date: 2026-05-24
closed_date: 2026-05-24
claim_status: B+ (pre-observational consistency + future-test target)
parent_motivation: "core/sek08a_akcja_zunifikowana §3838-3840 (Schwarzschild O(U³) deviation, native TGP prediction)"
authorization: "User 2026-05-24: 'Możemy zacząć od B' → 'działaj Phase 2' → 'ok, faza Final, ale zaznacz gdzieś problem z sekcji 8'"
independent_of: "γ-3, γ-3', γ-5, γ-7 F8 cycles (different physical observable: NS surface, NOT cosmological acceleration)"
related_PR: "PR-017 LOCKED-PENDING-FUTURE-PRECISION 2026-05-24"
R1_flags: "R1 #18 NEW (sek08a §3840 gauge ambiguity — substantive finding, future sek08c v3.0 scope)"
duration_actual: "1 sesja (4 phases batched)"
---

# op-PSR-orbital-drift — TGP O(U³) orbital-drift native falsifier

## Cycle scope (declared in Phase 0)

**Primary observable:** binary pulsar orbital dynamics (Hulse-Taylor B1913+16 + J0737-3039 double pulsar).

**Mechanism tested:** native TGP prediction of O(U³) deviation from Schwarzschild metric, per sek08a §3838-3840:

$$\Delta g_{00} = -\tfrac{1}{6}U^3 + \tfrac{1}{3}U^4 + \ldots$$
$$\Delta g_{rr} = \tfrac{1}{2}U^2 + \tfrac{5}{6}U^3 + \ldots$$

with U = GM/(c²r) Newtonian potential.

**Hypothesis:** TGP O(U³) deviation produces measurable residual in binary pulsar timing observables (periastron advance, orbital decay, Shapiro delay) beyond standard GR PPN to 2PN.

**OUT OF SCOPE** (anti-Lakatos discipline):
- F8 cosmological acceleration (separate cycle category — see [[../../meta/F8_FORENSIC_2026-05-24.md]])
- BH shadow predictions (already PR-012 / op-bh-alpha-threshold)
- LIGO/CE inspiral phase residuals (already PR-002)
- Galaxy rotation curves (already PR-004 SPARC)

## Status

**CLOSED-RESOLVED** (2026-05-24) — claim_status **B+** — PR-017 LOCKED-PENDING-FUTURE-PRECISION.

## Cycle outcome summary

| Falsifier | Verdict |
|-----------|---------|
| F-PSR-A | PASS (magnitude derivation procedure; Δz_poly(α, U) parametric) |
| F-PSR-B | PASS_CONSISTENT (NICER allows |α| ≤ 7.6 ≫ S07 0.832) |
| F-PSR-C | PASS (cross-system independence) |

**11 substantive items**: 9 PASS + 2 GAUGE_FINDINGS (R1 #18 sek08a §3840 gauge ambiguity).
**0 hardcoded T_pass=True ✓** | **Anti-Lakatos COMPLIANT ✓**

## Files

- [[Phase0_balance.md]] — pre-registration LOCKED
- [[Phase1_sympy.py]] + [[Phase1_derivation.md]] — F-PSR-A magnitude derivation
- [[Phase2_sympy.py]] + [[Phase2_results.md]] — F-PSR-B/C observational comparison
- [[Phase_FINAL_close.md]] — closure ceremony + R1 #18 registration

## Inheritance (LEGITIMATE)

- sek08a thm:einstein-emergence + thm:general-emergence (concept paper)
- sek08c metric form (V_M911, PPN γ=β=1)
- γ-5 Phase 3 G_eff = c³·ℓ_P²/ℏ (LOCKED)
- γ-7 Phase 1 q = √(4πG)·m Newton matching (LOCKED, PASS verdict)
- Appendix E eq. 353 m_sp scale (concept paper)

## Forbidden (anti-Lakatos)

- ❌ Cite F8 four-cycle FAILs as motivation
- ❌ Cite F8_FORENSIC_2026-05-24.md or E1/E2 explorations as positive evidence
- ❌ Inherit factor-10 threshold from γ-7 (use observational precision instead)
- ❌ Frame as "γ-8" or continuation of cosmology cycles
