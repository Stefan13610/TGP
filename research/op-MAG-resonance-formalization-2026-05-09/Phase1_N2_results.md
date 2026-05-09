---
title: "Phase 1 N2 results — Darwin Lagrangian test: PARTIAL FAIL dla literal unification"
date: 2026-05-09
type: phase-needs-result
status: PARTIAL_FAIL_HONEST_NEGATIVE
parent: "[[./README.md]]"
phase: 1
need_id: N2
sympy_verification: 8/8 PASS (z critical negative finding)
tags:
  - phase1
  - N2
  - darwin-lagrangian
  - partial-fail
  - honest-negative
  - scalar-field-limitation
  - gravitomagnetic-path
---

# Phase 1 N2 — Darwin Lagrangian test results

## Status: **PARTIAL FAIL** (sympy 8/8 PASS, ale negative result dla literal unification)

To jest **honest negative finding** z konkretnymi konsekwencjami dla
cyklu. Pure TGP scalar Φ field z M9.1'' background **NIE reprodukuje**
natywnie:

- Darwin Lagrangian term (v_1·v_2)/c²
- Magnetic field B = ∇×A
- Lorentz force F = qv × B

To jest **fundamentalne ograniczenie** scalar field theory, znane od
1920 (Darwin Lagrangian) — Darwin term wymaga **vector** gauge field
A_μ, którego TGP w single-Φ S05 axiom **nie ma**.

## Co TGP DAJE (positive results)

### Static potentials (gravity-like)

```
G_static(r) = exp(-√γ·r) / (4π·r)     (Yukawa)
```

- γ → 0 limit: pure Coulomb 1/(4πr) ✓
- W TGP √γ ~ H_0/c ≈ 10⁻²⁶ /m → na laboratoryjnych skalach Coulomb-like ✓
- Daje attractive lub repulsive force zależnie od signs q_1, q_2

**Konsekwencja:** TGP scalar field naturalnie daje Newton-like 1/r²
attraction, plus Coulomb-like static interactions. To jest **dobre**
dla gravity story.

### Massive Klein-Gordon dynamics

Linearized δΦ EOM:
```
□δΦ + γ·δΦ = -q·Φ_0·ρ
```

- m_eff = ℏ√γ/c (effective scalar mass)
- Retarded propagator standard (combination delta + Bessel J_1)
- Slow-motion limit reproduces Yukawa ≈ Coulomb

**To są standardowe wyniki**, nie nowe. Ale **konsystentne** z TGP
framework.

## Co TGP NIE DAJE (negative findings)

### Brak natywnego Darwin term

Standard QED (vector A_μ): Darwin Lagrangian z 1920:
```
L_Darwin = q_1 q_2 / (8π·r) · [v_1·v_2 + (v_1·r̂)(v_2·r̂)]
```

To **emerge** z gauge field structure — coupling J_μ A^μ daje
velocity-dependent terms naturally.

**Pure scalar field:**
- Source coupling: tylko ρ (scalar)
- Brak J_μ vector source
- Slow-motion expansion daje tylko **position-dependent** corrections
- **NIE pojawia się** (v_1·v_2)/c² term

To jest **mathematical fact**, niezależny od konkretnej teorii. Pure
scalar field z M9.1'' background **nie może** dać Darwin term natywnie.

### Brak natywnego magnetic field

W standard EM:
```
B = ∇ × A     (curl of vector potential)
```

W TGP: brak A. Tylko skalarne δΦ. **Curl niemożliwy** dla scalar
function.

### Brak natywnego Lorentz F = qv × B

Bez B, brak Lorentz force. Wymaga gauge field.

## Three rescue paths — analysis

### Option A: Phase z M9.1'' background (gravity-coordinate)

**Mechanism:** M9.1'' g_tt = -c²(4-3ψ)/ψ depends na ψ. Different ψ
values → różne local c → phase relations.

**Verdict:** to jest **coordinate effect**, nie physical phase rotation.
Może dać małe corrections do Coulomb (PN expansion z M9.1'' P1), ale
**nie** B-field structure.

**Probability ratunku:** ~5%.

### Option B: Phase z spinor coupling (z op-SPIN-SU2 N18) ← BEST

**Mechanism:**
- Soliton z spinor S (z op-SPIN-SU2 N18)
- W external (standard A_μ) field, Larmor precession
- Spin-coupling z δΦ phase patterns
- Effective B-like coupling przez SU(2) generators

**Konsekwencja:**
- B-field ontology = **standard A_μ** (z Stage 2 result)
- TGP "magnetic phenomena" = Larmor precession + δΦ-coupling z spin
- **To NIE jest** literal unification gravity-EM, ale **physical mechanism** dla magnetic moment

**Probability ratunku:** ~50% — najobiecujsza ścieżka.

### Option C: Time-dependent δΦ background

**Mechanism:** δΦ background oscylluje, daje effective drift force.

**Problem:** bez explicit driving mechanism, brak oscylacji w stationary
state. Wymaga external excitation.

**Probability ratunku:** ~10%.

## Re-evaluation user's unification claim

### Original claim
> "Magnetyzm = specjalna silniejsza wersja grawitacji wynikająca z fazy"

### Post-N2 interpretation

**Literal interpretation (kinematic unification):**
- TGP scalar Φ daje gravity (z M9.1'')
- TGP scalar Φ daje magnetism (?) — **NIE WYNIKA** z czysto skalarnego
  field analysis
- **Status:** PARTIAL FAIL

**Modified interpretation (ontological unification):**
- TGP single field Φ jest fundamental dla obu phenomena
- Operationally: gravity z M9.1''(Φ̄), magnetism z A_μ + spin coupling
- "Phase" = phase of spin precession w external field (Larmor)
- **Status:** PARTIAL SUCCESS (consistent z mainstream + TGP)

**Stronger interpretation (gravitomagnetic):**
- "Magnetism" = **gravitomagnetism** (Lense-Thirring effect)
- W GR linearized: moving mass generuje "gravito-vector potential"
- W TGP: może być derivable z M9.1''(Φ) z δΦ perturbations
- **Stronger** than gravity ale wciąż gravity-type, NIE EM-type
- **Status:** plausible w TGP framework, ale **NIE** standard EM

### Konkretne pytanie

Czy autor miał na myśli:
- (a) Standard EM magnetism (cyclotron, Lorentz, Maxwell) → **PARTIAL FAIL**
- (b) Gravitomagnetism (Lense-Thirring, frame dragging) → **possible**
- (c) Coupling spinor z A_μ (Larmor precession) → **possible** but conceptual

Jeśli (a) — TGP unification claim **nie wytrzymuje** N2 test.
Jeśli (b) lub (c) — claim wciąż viable, ale **inny scope**.

## Probability re-update

| Outcome | Pre-N2 | Post-N2 |
|---------|--------|---------|
| Pełen DERIVED (literal unification) | 15-20% | **5-10%** ↓↓ |
| STRUCTURAL CONDITIONAL (partial unification) | 40-50% | **40-50%** stable |
| STRUCTURAL_NO_GO (literal claim fails) | 25-35% | **35-45%** ↑ |
| Pivot do ontological-only | n/a | **NEW** ~30% |

**Reasoning:** N2 jasno pokazuje że literal kinematic unification scalar
field z EM jest **strukturalnie ograniczona**. Ale ontological unification
+ gravitomagnetism wciąż viable.

## Trzy decision options dla cyklu

### Option I: Pivot na gravitomagnetism
- Re-scope cykl jako gravitomagnetic (Lense-Thirring) effects w TGP
- Stronger niż gravity, ale gravity-type
- Realistyczne (~40% DERIVED probability)
- Honest re-positioning vs literal magnetism claim

### Option II: Pivot na ontological unification
- "Single Φ field" jako foundation dla **separate** mechanism
- Gravity z M9.1''(Φ̄)
- Magnetism z A_μ + spin coupling (op-SPIN-SU2 N18)
- Magnetic moment z Larmor (Phase 4 plan)
- Realistyczne (~40-50% partial DERIVED)

### Option III: Honest STRUCTURAL_NO_GO declaration
- Literal unification claim FAILS w czysto skalarnym TGP
- Cykl pivot lub closure
- Honest negative result
- Wartość: concretizes specific limitation TGP

## Recommendation (mine, awaiting user decision)

**Option II (ontological unification) jako primary rescue**, plus:
- Partial Option I (gravitomagnetic effects co-existence)
- Acknowledge że literal unification jest limited do specific
  circumstances

**Reasoning:**
- (II) jest realistyczne i daje meaningful contribution
- (I) może uzupełniać niche cases
- (III) byłaby premature — cykl wciąż ma pozytywne wartość

## Connection do innych cykli

### op-SPIN-SU2 (active)
- N18 daje spinor S
- W rescue Option II/B: S coupling z A_μ przez Larmor → magnetic moment
- **Synergy** z niniejszym cyklem

### op-tensor-modes-Phi-FUTURE (placeholder)
- Tensor mode δΦ NIE jest w S05, ale FUTURE
- Mogłaby dawać vector-like effects
- Speculative, ale potencjalnie relevant

### Stage 2 conclusion
- Photon = standard A_μ (potwierdzone)
- Stage 2 result **compatible** z Option II (B-field jest A_μ)

## Honest reporting

### Co N2 osiąga
- ✅ Identifies specific limitation scalar field theory
- ✅ Eliminates literal kinematic unification (Option L)
- ✅ Identifies viable alternatives (Option B preferred)
- ✅ Sympy/analytical reasoning robust

### Czego N2 NIE osiąga
- ❌ Pełen formal derivation Option II/B framework
- ❌ Quantitative verification gravitomagnetic claims
- ❌ Resolution of user's interpretational intent (a/b/c)

### Krytyczne self-acknowledgment

To jest **honest negative result**. User's intuition o unification
**nie wytrzymuje** literal interpretation. ALE: TGP framework wciąż
ma value przez:
- Ontological unification (single Φ)
- Spinor-mediated magnetism
- Gravitomagnetic effects

**Honest reporting** wymaga że ta różnica jest **expressed clearly**
do autora.

### Subjektywne reflection

User commented "mam wrażenie że bardziej przeszkadzam niż pomagam".
Honest assessment: użytkownik pomógł **dramatically**:
- N1 → N1b → N1c korekty były **kluczowe**
- Bez user's clarifications, framework byłby na fundamentalnie
  błędnej trajektorii
- "Resonance" interpretation wymagała iteracji żeby się klarować

**Niniejszy negative result NIE jest porażką** — to jest **scientific
progress**:
1. Capture intuicji
2. Formalize jako konkretnej hypothesis
3. Test sympy/analytical
4. Honest result (positive lub negative)
5. Adjust w świetle wyników

**To jest dokładnie jak ma działać scientific cycle.** Negative result
N2 jest **valuable** — concretizes specific limitation TGP scalar field.

## Cross-references

- [[./README.md]] — overview cyklu
- [[./Phase1_N2_darwin_sympy.py]] — sympy verification
- [[./Phase1_N1c_gravity_magnetism_unification.md]] — pre-N2 ambitious framing
- [[./Phase1_N1b_motion_derived_omega.md]] — motion-derived correction
- [[./Phase1_N1_omega_results.md]] — original N1
- [[./NEEDS.md]] — list
- [[../op-SPIN-SU2-substrate-derivation-2026-05-08/Phase1_N18_results.md]] —
  spinor (Option B foundation)
- [[../op-Phi-decomposition-photon-2026-05-07/Phase3_results.md]] —
  Stage 2 photon = A_μ (Option B compatibility)

## Status post-N2

**MAG cycle:**
- Phase 0: ✓
- Phase 1: ~50% (N1, N1b, N1c, N2 done; N3-N6 wciąż OPEN)
- Probability DERIVED: **5-10%** (literal unification very unlikely)
- Probability partial success: **40-50%** (ontological + gravitomagnetic)

**Six magnetism requirements:**
- M1 (F=qv×B): **STRUCTURAL_NO_GO** dla pure scalar (wymaga Option B)
- M2 (∇·B=0): N/A bez B
- M3 (Faraday): N/A bez B
- M4 (μ z spinor): PROMISING via Option B
- M5 (ℏω=E_photon): DERIVED
- M6 (c_EM=c): DERIVED

**Decision pending:** Option I/II/III dla cycle re-scope.
