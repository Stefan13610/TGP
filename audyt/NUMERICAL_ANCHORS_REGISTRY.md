---
title: "NUMERICAL ANCHORS REGISTRY — centralized tracking of TGP numerical coincidences"
date: 2026-05-16
type: registry
status: ACTIVE
purpose: "Track all numerical-anchor-class observations (NOT structural derivations) z honest classification per PHASE6 §11 pattern (CLOSED-NEGATIVE 2026-05-01)."
last_update: 2026-05-16
tags:
  - registry
  - numerical-anchor
  - structural-vs-numerical
  - PHASE6
  - honest-classification
parent: "[[./README.md]]"
related:
  - "[[./AUDIT_REPORT_2026-05-16_8-cycle_integration.md]] §10"
  - "[[../research/why_n3/PHASE6_alpha_em_connection.md]]"
---

# NUMERICAL ANCHORS REGISTRY

> **Purpose:** Centralized tracking of TGP numerical-anchor observations — OOM-close
> dimensional coincidences that **DO NOT** rise to structural derivation status
> (typically 10% precision threshold). Honest classification per PHASE6 §11 pattern
> (CLOSED-NEGATIVE 2026-05-01).
>
> **Distinction:**
> - **STRUCTURAL DERIVATION:** explicit mechanism + precision within ~10%
> - **NUMERICAL ANCHOR:** OOM-close (~factor 3) coincidence, no known mechanism
> - **STRUCTURAL HINT:** intermediate (~5-15% off, mechanism candidate present)
>
> This registry tracks **NUMERICAL ANCHORS only** — honest scientific reporting that
> dimensional coincidences exist but **do not constitute derivations**.

## §0 — Classification framework

### §0.1 — Anchor vs derivation criteria

| Criterion | DERIVATION | STRUCTURAL HINT | NUMERICAL ANCHOR |
|---|---|---|---|
| Precision | within ≤10% | 5-15% off | within ~OOM (factor 3) |
| Mechanism | explicit derivation | mechanism candidate | NO known mechanism |
| TGP-natural | yes (single-Φ + Z₂) | partial | dimensional only |
| Status | A− eligible | B+ partial closure | B+ partial / extension cycle |
| Action | promote to LOCKED | investigate further | document; defer to extension |

### §0.2 — Statistical caveat

For typical enumeration of K dimensional combinations across N OOM-wide span:
- Expected anchor hits w ±0.5 OOM window: `K × (1/N)` (roughly)
- Example: 12 combinations × 80-OOM span → ~0.15 hits expected by chance
- **Single anchor hit** in such enumeration: weak evidence, easily attributable to chance
- **Multiple aligned anchors** in cross-sector pattern: stronger evidence (still NOT derivation)

## §1 — Active numerical anchors (2026-05-16 baseline)

### §1.1 — Anchor #1: L08 e_Euler² in mass formula

**Source cycle:** [[../research/op-L08-Phase6-e2-derivation-2026-05-16/]]

**Anchor:**
```
m_obs = c · A_tail^p · g₀^(e²(1−α/4))
where e_Euler² ≈ 7.389056...
Best 0.02% fit dla lepton mass ratios (m_μ/m_e, m_τ/m_e)
```

**Status:** **NUMERICAL ANCHOR** (NOT structural derivation)

**History:**
- 2026-05-01: PHASE6_alpha_em_connection.md CLOSED-NEGATIVE — "X = e²/4 jest EMPIRICAL FIT
  w R3 amplitude sector z e_Euler statystycznym anchor"
- 2026-05-16 L08-e²-derivation cycle (B+ partial): algebraic reconciliation F1↔F2 derived;
  e_Euler² structural origin OPEN
- 2026-05-16 L08-RG-flow cycle (HALT-B): RG flow path 1+3 OBSTRUCTED; TGP α=2 = free
  massive field in canonical ψ=φ²; literature O(0.01-0.1) ≠ e²/2 ≈ 3.69 (factor 50-100 mismatch)
- **Status REINFORCED:** numerical anchor classification preserved post-2026-05-16

**Candidate structural mechanisms (all OBSTRUCTED or open):**
- (a) Yukawa tail integration ∫ exp(-2mr) — e appears, coefficient not natural
- (b) RG flow Z_φ(μ) at AS NGFP — **OBSTRUCTED** (L08-RG-flow HALT-B)
- (c) Partition function evaluation at S=−2 — arbitrary anchor
- (d) Topological winding × Berry phase — gives π, NIE e_Euler
- (e) Numerical coincidence — **currently most defensible**

**Future paths:**
- Hobart-Derrick α=4 generalization (multi-session)
- Statistical interpretation X = 1.847 ± δ (alternative anchor)
- Lattice substrate calculation (multi-session)

### §1.2 — Anchor #2: L06 (M_Pl²·H_0)^(1/3) in m_X

**Source cycle:** [[../research/op-L06-axion-mass-derivation-2026-05-16/]]

**Anchor:**
```
(M_Pl² · H_0)^(1/3) = ((1.22·10²⁸)² × 1.5·10⁻³³)^(1/3) ≈ 6.07·10⁷ eV ≈ 60 MeV
vs target m_X = 100 MeV (phenomenological from ψ.1, ω.2)
Distance: factor 1.7 (40% off); within ±0.5 OOM anchor window
```

**Status:** **NUMERICAL ANCHOR** (NOT structural derivation)

**History:**
- 2026-05-16 L06-axion-mass-derivation cycle (B+ partial): 4 derivation paths (A-D)
  OBSTRUCTED z explicit proofs; Path E (FREE PARAMETER) confirmed strukturalnie
- Dimensional enumeration: 12 combinations tested w 80-OOM span;
  `(M_Pl²·H_0)^(1/3)` = ONLY hit within ±0.5 OOM tolerance (1 of 12 — statistically possible by chance)

**Candidate structural mechanisms (none known):**
- (a) Gravity-mediated SUSY-breaking analog: m_soft³ ~ M_Pl²·H_SUSY (requires SUSY, not in TGP)
- (b) Quintessence tracker scale (Hill-Steinhardt 1988): tracking H during slow-roll
- (c) See-saw analog z some hidden mass scale: requires additional assumption
- (d) Numerical coincidence — **most likely** given:
  - 40% off pre-registered 10% precision threshold
  - NO TGP-native mechanism connecting m_X to (M_Pl²·H_0)^(1/3)
  - Dimensional power 1/3 NIE natural in standard QFT (usually 1/2 from kinetic terms)

**Cross-cycle context:**
- ψ.1 uses m_X = 100 MeV (Yukawa range 1.97 μm SNR optimization)
- τ.3 uses m_X = g·f_X = 0.83 MeV (anomaly coupling product, f_X phenomenological)
- L06 confirms m_X = FREE PARAMETER strukturalnie per audit § A.7 option 2 + ω.3 endorsement
- (M_Pl²·H_0)^(1/3) ≈ 60 MeV does NOT match either phenomenological value

**Future paths:**
- Gravity-mediated mass mechanism investigation (extension cycle if motivated)
- Possible cross-pattern: ratio (M_Pl²·H_0)^(1/3) / Λ_QCD = 60/217 ≈ 0.28 — random?

## §2 — Pattern observations

### §2.1 — Cross-anchor consistency check

| Property | Anchor #1 (e_Euler²) | Anchor #2 ((M_Pl²·H_0)^(1/3)) |
|---|---|---|
| Cycle | L08-e² + L08-RG | L06 |
| Discovered | 2026-05-01 / reinforced 2026-05-16 | 2026-05-16 |
| Distance from target | 0.02% (best fit) | 40% (factor 1.7) |
| Mathematical structure | transcendental (e) | algebraic (cube root of product) |
| Sector | mass exponent (leptons) | axion mass (ALP) |
| Cross-sector pattern? | NO | NO |
| Status | NUMERICAL ANCHOR | NUMERICAL ANCHOR |

**Observation:** Two anchors with very different mathematical structure and very different
distance-from-target. No cross-sector pattern. Both honestly classified as numerical
coincidences.

### §2.2 — Disposition policy

**For each numerical anchor identified:**
1. Document in this registry z full context
2. Test against pre-registered structural-derivation threshold (typically 10%)
3. If fails: classify as NUMERICAL ANCHOR; document candidate mechanisms
4. NO post-hoc fitting allowed (forbidden direction per cycle convention)
5. Extension cycle ONLY if mechanism candidate becomes available (research-track open)

**Registry updates:**
- Add new entries when cycles discover OOM-close coincidences
- Promote to "STRUCTURAL HINT" or "DERIVED" if mechanism subsequently found
- Demote / archive entries if mechanism FALSIFIED (analog L08-RG-flow HALT-B for e_Euler²)

## §3 — Recommended future actions

### §3.1 — Pursue / monitor (low priority)

- **e_Euler²:** Hobart-Derrick α=4 generalization (multi-session); statistical
  interpretation X = 1.847 ± δ pursuit
- **(M_Pl²·H_0)^(1/3):** gravity-mediated mass mechanism investigation if motivated;
  ratio analysis vs other TGP scales

### §3.2 — Do NOT pursue

- Post-hoc reverse-engineering mechanism to fit anchor (analog L08 audit §11 caveat)
- Free-parameter introduction to bridge anchor↔derivation gap (S05 violation)
- Promotion to DERIVED without explicit mechanism + 10% precision

### §3.3 — Honest reporting requirements

Every cycle finding numerical anchor must:
- Pre-register precision threshold (typically 10% dla derivation; ±0.5 OOM dla anchor)
- Document explicit obstruction proofs dla candidate mechanisms tested
- Add entry do tej registry z full context
- Cross-link bidirectionally (cycle Phase_FINAL_close ↔ registry entry)

## §4 — Statistical baseline

**Empirical pattern (sesja 2026-05-16):**
- 8 derivation cycles + 1 housekeeping
- 2 numerical anchors documented (25% rate per cycle for anchor-finding)
- 9 explicit obstruction proofs documented (1.1 per cycle average)
- Pattern suggests TGP framework has **clear structural-vs-numerical distinction**
  that should be tracked systematically

**Future projection:**
- Expected ~2-3 numerical anchors per 10 derivation cycles (rough estimate)
- Registry should be reviewed quarterly; archive obsolete entries; promote confirmed mechanisms

## §5 — Cross-references

- [[./README.md]] — audyt index
- [[./AUDIT_REPORT_2026-05-16_8-cycle_integration.md]] §10 — pattern recognition that
  motivated creation of this registry
- [[../research/why_n3/PHASE6_alpha_em_connection.md]] — CLOSED-NEGATIVE 2026-05-01;
  original numerical anchor classification pattern
- [[../research/op-L08-Phase6-e2-derivation-2026-05-16/Phase_FINAL_close.md]] — Anchor #1 source
- [[../research/op-L08-Phase6-RG-flow-Z-phi-asymptotic-2026-05-16/Phase_FINAL_close.md]] —
  Anchor #1 mechanism obstruction proof
- [[../research/op-L06-axion-mass-derivation-2026-05-16/Phase_FINAL_close.md]] — Anchor #2 source
- [[../PREDICTIONS_REGISTRY.md]] §"Updated 2026-05-16" — sesja-level cross-reference

---

**Maintainer:** Add new anchors as cycles discover them; review quarterly; archive obsolete.
