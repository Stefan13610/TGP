# OP-M92 — Decision tree synthesis (T-M92.6)

**Data:** 2026-04-25
**Status:** OPEN — readiness package for ngEHT 2030+ verdict
**Parent:** [[research/op-m92/OP_M92_setup.md]]
**Cross-refs:** [[research/op-eht-A/OP_EHT_A_final_verdict.md]]

---

## TL;DR

> Po OP-EHT-A NEGATIVE, M9.1'' standalone falsifiable w 4-6 lat (ngEHT 2030+).
> OP-M92 to **deferred conditional research program** z pre-analyzed
> kandydatami (T-M92.1..T-M92.5 done jako readiness package).
>
> **Główny wniosek strategiczny:** wybór M9.2 jest data-driven po ngEHT
> verdict, ale ranking kandydatów już zrobiony:
>
> 1. **Candidate D (momentum back-reaction):** PROMISING — natural strong-field
>    activation, auto-PASS weak-field PPN, physically motivated (Lenz law).
> 2. **Candidate A (dual-field):** VIABLE z screening (Vainshtein/chameleon).
> 3. **Candidate B (conformal frame):** VIABLE ale tightly constrained.
> 4. **Candidate C (q-flow):** NOT VIABLE bez ψ-threshold mechanism.
>
> **Estimated response time post-verdict:** 2-4 weeks dla candidate selection
> i 6-12 months dla full M9.2 paper, vs naive 2-3 lat from-scratch.

---

## 1. Decision tree (formal)

```
                           ngEHT 2030-2032 verdict
                                     |
                  ┌──────────────────┼──────────────────┐
                  |                  |                  |
        b/b_GR ∈ [0.99, 1.01]   ∈ [0.90, 0.99]    ∈ [1.10, 1.20]
            "GR shadow"        "Intermediate"    "TGP-like"
                  |                  |                  |
            M9.2 PIVOT        M9.1''' refinement   M9.1'' VINDICATED
            MANDATORY         (intermediate pivot)
                  |                  |                  |
        ┌─────────┴───┐               |                  |
        |             |              |                  |
    Candidate D   Candidate A/B  Modify ψ_ph        OP-EHT upgrades
    (priority)    (fallback)     structure          CONDITIONAL → POSITIVE
        |             |              |                  |
    Phase 1:      Phase 2:       6-12 month         M9.1'' becomes
    structural    if D fails     research            standard theory
    derivation    Vainshtein     program             with EHT validation
    + 1PN check   screening etc.
        |             |
    OP-M92        OP-M92
    closure       closure
    POSITIVE      POSITIVE
```

**Decision thresholds** (numerical):
- "GR shadow": b_crit_observed/b_crit_GR - 1 ∈ [-0.01, +0.01]
- "Intermediate": ∈ [-0.10, -0.01] ∪ [+0.01, +0.10]
- "TGP-like": ∈ [+0.10, +0.20] (M9.1'' standalone predicts +0.1456)

## 2. Probability assessment (Bayesian, post-OP-EHT-A)

| Outcome | Prior | Post-OP-EHT-A | Reasoning |
|---------|-------|---------------|-----------|
| GR shadow at 1% | 50% | 50% | OP-EHT-A doesn't shift this (it's about TGP not GR) |
| Intermediate band | 30% | 30% | Allows partial M9.1''' refinement |
| TGP-like deviation | 20% | 20% | M9.1'' standalone success path |

**Most likely outcome: GR shadow** → M9.2 pivot mandatory in 2030-2032.

## 3. Candidate D (momentum back-reaction) — preferred path

**Why preferred:**

1. **Strong-field activation natural**: T^μν ~ U² in weak-field, ~1 in strong-field
   → α × T·J·J term grows non-linearly z curvature (no fine-tuning).

2. **Weak-field auto-pass**: O(U^4) suppression w PN expansion → γ_PPN, β_PPN
   = 1 at 1PN bez constraints na α coupling.

3. **Physically motivated**: substrate physics analog Lenz law w EM —
   moving charge induces back-reactive field opposing motion. W substrate
   physics, "moving" matter induces back-reactive Φ-flow opposing curvature growth.

4. **First-principles feasible**: action S = ∫ T^μν J_μ J_ν √-g d⁴x ma
   natural variational derivation z higher-curvature gravity literature.

**Open challenges (T-M92.5 scope):**

- **Ostrogradsky instability**: higher-derivative terms after IBP. Mitigation:
  use degenerate Lagrangian construction (à la Horndeski/DHOST).
- **Ghost modes**: requires careful sign analysis of α coupling.
- **Quantization**: classical-only theory until renormalizability checked.

**Phase 1 plan (post-verdict):**

1. Variational derivation T^μν J_μ J_ν → modified Φ-EOM
2. Photon ring analysis: does α tune to scenario (e) at ψ_ph?
3. PN expansion to O(U^4): verify weak-field auto-pass
4. Stability: linearize around Schwarzschild, check dispersion relation
5. OP-7 c_GW = c_0 cross-check (likely auto-pass, σ_ab kinetic invariant)

**Phase 2 plan (if Phase 1 PASS):**

1. Full M9.2 axiom paper draft
2. Cosmological consequences (does back-reaction affect H_0?)
3. Falsifiable predictions distinct from M9.1'' (e.g., black hole secondary
   photon ring, gravitational lensing tail)
4. OP-9 closure under M9.2

## 4. Fallback path (if Candidate D fails Phase 1)

**Candidate A (dual-field) Phase 1:**

1. Introduce auxiliary scalar σ(x) z action α(∇σ)² + V(Φ, σ)
2. Vainshtein screening derivation (à la massive gravity/galileon)
3. PPN audit: γ_PPN = 1 + O(1/ω_BD) — check Cassini bound
4. Photon ring analysis with σ-dressed metric

**Candidate B (conformal frame) Phase 1:**

1. Determine Ω(ψ) profile s.t. Ω²(1)=1, Ω²(1.168)=0.785
2. Equivalence principle: composition-dependent acceleration test
3. PN expansion: γ_PPN, β_PPN bounds
4. Cosmological viability (does Ω-flow break ΛCDM at late times?)

## 5. Timeline (formal)

| Period | Phase | Deliverable |
|--------|-------|-------------|
| 2026-04-25 | OP-M92 OPEN | Setup + diagnostic T1 + candidate map (DONE) |
| 2026 Q3-Q4 | Phase 0 | Candidate D variational derivation sketch |
| 2027 H1 | Phase 0 | Candidate D PN expansion + stability |
| 2027 H2 | Phase 0 | Candidate A/B fallback sketches |
| 2028-2029 | Phase 0 | OP-M92 paper draft (deferred conditional) |
| 2030 | TRIGGER 2 EVAL | ngEHT first Sgr A* results — verdict pending |
| 2031-2032 | RESPONSE | Decision tree execution; selected M9.2 in 2-4 weeks |
| 2032+ | Phase 1-2 | Full M9.2 paper / M9.1'' vindication paper |

## 6. Why this is "deferred conditional research"

OP-M92 jest **deferred** — żadna critical decision pre-2030. Powód:
- Bez ngEHT data, nie wiadomo czy M9.2 jest needed (50% probability of M9.1''' or vindication paths)
- M9.1'' standalone pozostaje active baseline z weak-field PASS
- 4 lat to luxury w teorii fundamental physics — pozwala na deeper structural analysis
- Pre-analysis kandydatów ELIMINUJE response delay 6-12 months → 2-4 weeks

OP-M92 jest **conditional** — outcome zależy od external observation (ngEHT):
- Trigger 1 satisfied (Track A failed) ✓
- Trigger 2 conditional na ngEHT verdict
- Multiple branches: pivot / refinement / vindication

OP-M92 jest **research program**, nie single closure:
- T-M92.1 done (diagnostic)
- T-M92.2..T-M92.5 done (candidate map)
- T-M92.6 done (this decision tree)
- Full closure deferred do 2030+

## 7. Cross-program implications

### F4 falsifiability hardening (current state, post-OP-M92 OPEN)

> **TGP M9.1'' standalone predykuje +14.56% strong-field shadow deviation.
> ngEHT 2030+ at 1% Sgr A* precision daje >10σ verdict. Jeśli ngEHT confirms
> GR → M9.2 pivot z pre-analyzed Candidate D (momentum back-reaction)
> as priority, response time 2-4 weeks. Jeśli ngEHT confirms TGP-like →
> M9.1'' VINDICATED, OP-EHT upgrades CONDITIONAL → POSITIVE.**

### Resource allocation (post-OP-M92)

- **75% effort:** continue OP-DESI, OP-Hubble (cosmology — independent)
- **20% effort:** OP-M92 Phase 0 (candidate readiness for 2030+)
- **5% effort:** OP-EHT continuous monitoring (EHT 2026-2030 incremental data)

## 8. Bottom line

OP-M92 nie jest crisis response. Jest pre-emptive readiness program
z explicit timeline 4-6 years. Po OP-EHT-A NEGATIVE strukturalnie:

- M9.1'' alive standalone z weak-field PASS
- M9.2 pivot mandate **conditional** na ngEHT 2030+ verdict
- Pre-analyzed candidates → response time 2-4 weeks (vs 6-12 months)
- Candidate D (momentum back-reaction) leads as natural strong-field
  rescue z auto-preserved weak-field PPN

**Falsifiability F4 staje się jeszcze stronger:** TGP M9.1'' już *dziś*
ma explicit decision tree z named candidates dla ngEHT response, co
demonstrates że teoria jest nie tylko falsifiable ale **operationally
prepared** dla obu outcomes.
