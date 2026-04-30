---
title: "ω.1 program — substrate ↔ EM back-reaction coupling (axion-like topological)"
date: 2026-04-30
cycle: ω.1
status: ACTIVE
parent: "[[../op-phi1-substrate-action-variational/Phase3_results.md]]"
tags:
  - TGP
  - omega1
  - substrate-em-coupling
  - axion-like
  - topological
  - back-reaction
  - structural
---

# ω.1 — substrate ↔ EM back-reaction coupling

**Goal:** Determine whether U(1) electromagnetism can structurally back-react on
the substrate field X (closing the G3 gap from `em_from_substrate/PLAN.md`),
identify the **unique non-trivial coupling form** compatible z φ.1 substrate-action
+ X→λX scale invariance, and propose falsifiable experimental signatures.

This addresses the user-flagged research question: *"czy istnieje strukturalne
połączenie między elektromagnetyzmem a polem generowanej przestrzeni; czy
istnieje szansa na oddziaływanie na samą przestrzeń za pomocą elektromagnetyzmu."*

## Hypothesis

The minimal extension of φ.1 substrate-action that admits non-trivial
EM ↔ substrate coupling z preserwacja gauge invariance + Noether scale-symmetry
ma postać **axion-like topological**:

$$\mathcal{L}_{\omega.1} = -\tfrac{1}{4}F_{\mu\nu}F^{\mu\nu}
   + \tfrac{1}{2} f_X^2 (\partial_\mu \ln X)(\partial^\mu \ln X)
   + \tfrac{g}{4} (\ln X) \, F_{\mu\nu}\tilde{F}^{\mu\nu}$$

gdzie $\tilde{F}^{\mu\nu} = \tfrac{1}{2}\epsilon^{\mu\nu\rho\sigma} F_{\rho\sigma}$.

**Critical math correction (vs naive minimal coupling):**

Naive minimal coupling $\mathcal{L}_{\text{minimal}} = e A^\mu J_\mu$
z $J_\mu = \partial_\mu \ln X$ jest **gauge-trivial**:

$$\int A^\mu \partial_\mu (\ln X) \, d^4x
 = -\int (\ln X) \partial_\mu A^\mu \, d^4x + \text{boundary}
 \xrightarrow{\text{Lorenz gauge}} 0$$

Stąd jedyny **strukturalnie nietrywialny** coupling kompatybilny z:
1. φ.1 scale-symmetry (X → λX → ln X → ln X + ln λ, F·F̃ niezmiennicze)
2. U(1) gauge invariance (F·F̃ jest total divergence = topological density)
3. Lorentz invariance + parity-odd channel (F·F̃ jest pseudoscalar)

to **dokładnie axion-like coupling** $\propto (\ln X) F\tilde{F}$.

## Reference frame

| Source | Status quo | ω.1 lift |
|---|---|---|
| em_from_substrate U(1) gauge | EM emerges FROM substrate (one-way) | ω.1 closes back-reaction |
| em_from_substrate G3 gap | "J_phase ↔ J_amp" unanswered | ω.1 explicit Lagrangian |
| φ.1 substrate-action | scalar EL `□(ln X) = 0` | ω.1: source term `-(g/4f_X²) F·F̃` |
| XS.1 κ_TGP ≈ 2.012 | cross-sector SC charge | g LOCK candidate |
| α_em ≈ 1/137.036 | EM coupling anchor | alternate g LOCK candidate |

## Phase plan (5 + 7 + 6 = 18 sub-tests)

### Phase 1 — coupling-form scan + structural invariance (5 sub-tests)

- **W1.1** — Coupling-form scan (4 candidates: minimal A·J trivial, dilaton (ln X)F², axion (ln X)F·F̃, gradient (∂u)²F²)
- **W1.2** — Gauge invariance check (axion form: ∫F·F̃ = total divergence)
- **W1.3** — Scale invariance X → λX (axion form: invariant mod boundary)
- **W1.4** — Coupling constant g LOCK candidates (κ_TGP, α_em, α₀, η_chir)
- **W1.5** — Phase 1 gate (≥4/5)

### Phase 2 — sympy LOCK + modified EOMs (7 sub-tests)

- **W2.1** — sympy LOCK modified Maxwell: `∂_ν F^νμ = -g F̃^μν ∂_ν(ln X)`
- **W2.2** — sympy LOCK modified substrate EOM: `f_X² □(ln X) = -(g/4) F·F̃`
- **W2.3** — Bianchi identity preservation: `∂_μ F̃^μν = 0` unchanged
- **W2.4** — Alt-coupling falsification (3 alts: dilaton, gradient, minimal)
- **W2.5** — Lorenz-gauge consistency check
- **W2.6** — Stress-energy T^μν derivation + scale-current modification
- **W2.7** — Phase 2 gate (≥6/7)

### Phase 3 — predictions + 4-channel convergence (6 sub-tests)

- **W3.1** — Vacuum birefringence prediction (PVLAS/OSQAR lab limits)
- **W3.2** — Magnetar B² substrate sourcing (∼10¹⁵ G fields, B·E channel)
- **W3.3** — CMB cosmological birefringence Δχ (Planck/SO/LB constraints)
- **W3.4** — Quasar polarization rotation vs redshift z (radio + optical)
- **W3.5** — Alt-coupling cross-channel falsification (dilaton/minimal predictions FAIL)
- **W3.6** — 4-channel ω.1 convergence (gauge + scale + EOM + experimental)

## Falsification path

ω.1 jest falsyfikowane przez:
1. Pomiar g = 0 z dokładnością wykluczającą topological coupling (PVLAS-like ground)
2. CMB/quasar birefringence Δχ inconsistent z g·(ln X) sourcing
3. Alt-coupling form (dilaton or minimal) better-fit niż axion

## Cross-references

- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — φ.1 substrate-action AXIOM
- [[../op-upsilon1-closure-cross-family/Phase3_results.md]] — υ.1 closure law DERIVED
- [[../../em_from_substrate/PLAN.md]] — G3 gap (J_phase ↔ J_amp)
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
