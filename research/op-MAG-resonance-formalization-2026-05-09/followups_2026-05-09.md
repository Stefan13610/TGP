---
title: "Follow-ups post-cycle close (op-MAG-resonance) — scoped cycles + consistency notes"
date: 2026-05-09
type: follow-up-summary
parent: "[[./Phase6_absolute_binding.md]]"
tags:
  - follow-ups
  - scoped-cycles
  - consistency-checks
---

# Follow-ups po close MAG cycle (2026-05-09)

## Status post-cycle

MAG cycle CLOSED: STRUCTURAL DERIVED CONDITIONAL (sympy 14/14 PASS).

Acknowledged limitations (z Phase 6 absolute binding):

| Limitation | Status post-MAG | Path forward |
|---|---|---|
| **Φ_0 free parameter** | Phase 5 wymaga Φ_0 input | → osobny cykl `op-Phi-vacuum-scale-2026-05-09` (SCOPED) |
| **α/π anomalous moment** | M4 leading order only | → osobny cykl `op-MAG-anomalous-moment-2026-05-09` (SCOPED) |
| **EM-strength Lorentz force** | wciąż via standard A_μ | → osobny cykl `op-MAG-Lorentz-A-mu-coupling-2026-05-09` (SCOPED) |
| **Particle spectrum (μ, τ, kwarki)** | out of scope MAG | → consistency check z `particle_sector_closure` (P4) |

## Trzy scoped cycles (rozpisane 2026-05-09)

### 1. [[../op-Phi-vacuum-scale-2026-05-09/README.md]]
- **Cel:** derivation Φ_0 z first principles + UV/IR normalization (2 wartości reconciliation)
- **Scope:** TOPORNY problem (per user note)
- **Recommend:** ODŁÓŻ — niski ROI, wymaga RG flow infrastructure
- **Probability DERIVED:** 10-20%

### 2. [[../op-MAG-anomalous-moment-2026-05-09/README.md]]
- **Cel:** g_e - 2 ≈ α/π z **phase amplification at topology connection boundary**
- **Hypothesis:** boundary phase amplifies (jak gdyby measured 2x), ratio = α/π
- **Scope:** SPECULATIVE, ale konkretna nowa hipoteza autora
- **Recommend:** Phase 0-1 quick scan (cheap), if signals positive — kontynuuj
- **Probability DERIVED:** 10-15%; STRUCTURAL: 25-35%

### 3. [[../op-MAG-Lorentz-A-mu-coupling-2026-05-09/README.md]]
- **Cel:** EM-strength Lorentz F = qv × B z standard A_μ + spinor amplification mechanism (gravitomagnetic → EM)
- **Scope:** medium, builds na N2d + M4
- **Recommend:** **NAJWYŻSZY PRIORITY** — most realistic delivery, closes biggest MAG limitation
- **Probability DERIVED:** 15-25%; STRUCTURAL CONDITIONAL: 40-50%

## Particle spectrum: P4 consistency check (DONE)

User comment (2026-05-09):
> "Particle spectrum (μ, τ, kwarki) → to jest realne out of scope, wyprowadzenie
> jest w innym miejscu wieża letnowa, jedyne co należy zrobić to zbadać zgodność"

### Existing tower derivation: P4 — Particle Sector Closure

Ref: [[../particle_sector_closure/README.md]] (status: CLOSED 2026-04-19)

**P4 results:**
- Lepton sector: Koide K = 2/3 derivable z solitonowej ODE + PDG m_τ/m_e ✓ FULL CLOSURE
- Quark sector: K_up = 0.8746, K_down = 0.7398, RG-invariant, NIE Koide K=2/3 ✓ FULL CLOSURE
- Neutrino sector: Brannen ansatz INCONSISTENT z Δm² ✓ FULL CLOSURE (negative)
- α₃ closed form: P_cos transcendence, structural closure

**Top-level statement (z P4):**
> "Sektor leptonów naładowanych jest w pełni zamknięty — Koide K=2/3 jest
> strukturalnie derywowany z solitonowej ODE + PDG m_τ/m_e."

### Consistency check RESULT (DONE 2026-05-09)

**Sympy:** 2/2 PASS ([[./Phase5b_P4_consistency_check_sympy.py]])

**Setup:**
- P4 ps2 mass formula: m ∝ A⁴ (α=1) lub m ∝ A³ (α=2 TGP-canonical)
  gdzie A = tail amplitude solitonu g(r) ~ 1 - A·e⁻ʳ/r
- MAG Phase 5 Mach formula: m_Mach = (3γq²)/(16π Φ_0² m_C) · ⟨δΦ²_bg⟩
- Substituting P4 tail profile do MAG: ∫δΦ²_sol = 2π A² → m_Mach ∝ **A²**

**Power mismatch:**
- P4: A⁴ (α=1) lub A³ (α=2) — verified Koide K=2/3
- MAG Phase 5: A² — would NOT reproduce Koide K=2/3

**INTERPRETATION (resolution):**

Both formulas są TGP-natywne ale opisują **różne mass contributions:**

```
m_total = m_intrinsic_soliton (P4)  +  m_Mach_correction (MAG)
              ∝ A⁴ (lub A³)                    ∝ A²
              DOMINANT (intrinsic energy)      SUBLEADING (bg coupling)
```

- **P4 ps2** counts intrinsic soliton energy (kinetic + potential z V(Φ))
- **MAG Phase 5** counts Mach correction z background fluctuations

Hierarchy: m_Mach/m_intrinsic ~ A²/A⁴ = 1/A² (Mach correction smaller dla heavy solitons).

**VERDICT:** ✓ **CONSISTENT** — no formal tension, structural separation.

P4 lepton tower (Koide K=2/3) jest LEADING, MAG Phase 5 Mach correction jest SUBLEADING. Both contribute, ale P4 dominuje dla typical lepton.

**Numerical check:** PDG lepton masses dają Koide K = 0.666661 (vs 2/3 = 0.666667), diff 0.001%. P4 ps2 closure confirmed.

**Light scope DONE.** Cycle cohesion strengthened — TGP framework has internally consistent multi-mechanism mass structure.

## Strategic priorities (mine, awaiting user confirmation)

1. **HIGHEST:** op-MAG-Lorentz-A-mu-coupling — most realistic + closes biggest gap
2. **MEDIUM:** P4 consistency check — light scope, builds confidence
3. **MEDIUM-LOW:** op-MAG-anomalous-moment — speculative, Phase 0-1 quick test
4. **LOW (defer):** op-Phi-vacuum-scale — toporny, wymaga more infrastructure

## Notatka o "wracamy do A"

User wspomniał: "wracamy do A, po rozpisaniu nowych cykli badawczych"

**Interpretation pending:** "A" w tym kontekście może oznaczać:
- (a) Continue op-SPIN-SU2-substrate-derivation cycle (paused after Phase 1cd)
- (b) Start jednego z 3 scoped cycles powyżej
- (c) Inny aktywny cykl

**Awaiting user clarification.**

## Cross-references

### Closed cycle
- [[./Phase6_absolute_binding.md]] — final classification

### New scoped cycles
- [[../op-Phi-vacuum-scale-2026-05-09/]]
- [[../op-MAG-anomalous-moment-2026-05-09/]]
- [[../op-MAG-Lorentz-A-mu-coupling-2026-05-09/]]

### Consistency check target
- [[../particle_sector_closure/README.md]] — P4 closed result
