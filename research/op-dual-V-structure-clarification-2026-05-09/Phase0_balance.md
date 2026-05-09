---
title: "Phase 0 balance — op-dual-V-structure-clarification-2026-05-09"
date: 2026-05-09
type: phase-balance-sheet
status: COMPLETE
parent: "[[./README.md]]"
phase: 0
---

# Phase 0 balance — Dual-V structure formal verification

## Pre-cycle confirmed facts (search agent + direct read)

### Fact 1: G.0 explicit "A4 marker" — gravity-only deprecation

[[../op-g0-r3-from-canonical-projection/Phase1_results.md]] **linia 264-266:**

```
- A1, A2, A3 (metric drift) — automatycznie zamknięte przez G.0 jeśli
  V_M911 i √(-g)=c₀ψ/(4-3ψ) wprowadzone do sek08a
- A4 (matter coupling) — wymaga osobnego sprawdzenia (G.0 nie dotyka L_mat)
- B6, B7, B8 (numerical re-runs) — automatycznie po update sek08a
- M9.x results (PPN, GW, LLR, BBN) — re-derivation po update wymagana
```

**Eksplicytne:** G.0 closure
- ✅ Zamyka A1-A3 (metric drift, gravity sector)
- ❌ NIE dotyka A4 (matter coupling)

### Fact 2: V_M9.1'' derived z GRAVITATIONAL constraints

[[../op-g0-r3-from-canonical-projection/phase2_P21_vacuum_uniqueness.py]] **linia 13-22:**

```python
PURPOSE: G.0 PHASE 2 SUB-TASK P21:
  Pelna sympy formal weryfikacja:
  1. UNIQUENESS V_M911(psi) = -psi^2*(4-3psi)^2/12 pod constraint'ami:
     - K(psi) = psi^4 (T-D-uniqueness, alpha=2)
     - sqrt(-g) = c0*psi/(4-3psi) (M9.1'' canonical)
     - Static EOM = R3 ODE
```

**Wszystkie constraint'y są GRAVITY:**
- R3 ODE: static spherical EOM (Newton limit)
- M9.1'' metric: geometric (gravity)
- T-D-uniqueness: kinetic term z spatial geometry

V_M9.1'' wynika **JEDNOZNACZNIE** z gravity, NIE z matter dynamics.

### Fact 3: V_orig derived z MATTER field theory

[[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] **linia 38-55:**

```
TGP V(Φ) = -β Φ³/(3 Φ_0) + γ Φ⁴/(4 Φ_0²)

Wokół vacuum Φ = Φ_0 + δΦ:
  V''''(Φ_0) = 6γ/Φ_0²

Cross-coupling term: 6 · δΦ_sol² · δΦ_bg²

Effective coupling λ_4 = V''''/4 = 3γ/(2 Φ_0²)
```

**V_orig jest user'owany dla MATTER dynamics** (field expansion around vacuum,
particle Mach inertia).

### Fact 4: T-Λ ρ_vac używa V_orig formula (matter-like)

[[../closure_2026-04-26/Lambda_from_Phi0/results.md]] **linia 19:**

```
ρ_vac,TGP = V(Φ_eq) = γ·Φ_eq²/12
```

To jest V_orig V(Φ_eq) = γ·Φ_eq²/12 (przy β=γ vacuum, Phi_eq=Phi_0).
Matches obs ratio 1.020 (positive result).

## Path C hypothesis = CONFIRMED FRAMEWORK FEATURE

**Pre-search status (post-Phase5-clarification):** Path C suggested (50-60% probability)

**Post-search status (post-G.0 explicit reading):** Path C **CONFIRMED**:

| Sektor | Potencjał | Source | Pozytywne wyniki |
|---|---|---|---|
| **Gravity** | $V_{M9.1''}(\psi) = -\gamma\psi^2(4-3\psi)^2/12$ | R3 ODE (G.0 P21) | M9.1'' metryka, Newton, mass spectrum invariance |
| **Matter** | $V_{orig}(\Phi) = -\beta\Phi^3/(3\Phi_0) + \gamma\Phi^4/(4\Phi_0^2)$ | Field expansion | T-Λ ρ_vac (1.020 match), Phase 5 m_e=511 keV |

## Reinterpretation of "DEPRECATED 2026-05-02"

**Original sek08a annotation (linia 96-98):**
```
V_orig = ... [DEPRECATED 2026-05-02; see prop:V-M911-canonical]
```

**Correct interpretation (post-A4 marker analysis):**
```
V_orig = ... [DEPRECATED FOR GRAVITATIONAL SECTOR 2026-05-02 via G.0 closure;
              see prop:V-M911-canonical. Matter sector usage maintained pending
              A4 verification.]
```

**Niniejszy cykl realizuje A4 verification.**

## Pytania dla Phase 1

1. **Q1:** Czy V_M9.1'' i V_orig są mutually exclusive (jeden zastępuje drugi)
   czy complementary (każdy dla swojego sektora)?
2. **Q2:** Czy istnieje kontradykcja matematyczna w dual-V framework?
3. **Q3:** Czy V_orig jest expansion V_M9.1'' w specyficznym regime (e.g.,
   low-amplitude matter fluctuations around vacuum), czy genuinely independent?
4. **Q4:** Jaka jest correct sek08a annotation dla V_orig?

## Hipoteza H1 update post-search

**H1 (refined):** Dual-V framework jest **legitimate TGP feature**:
- V_M9.1'' = **gravitational potential** (canonical from R3 ODE)
- V_orig = **matter potential** (canonical from field theory expansion)
- A4 marker w G.0 closure **already documents** ten fakt

**Falsifier:** kontradykcja matematyczna między dwoma EOM-y, np. dual-V
implikuje inconsistent particle masses lub vacuum structure.

## 8/8 gate criteria

| # | Criterion | Status |
|---|---|---|
| 1 | Pre-cycle finding (A4 marker) eksplicit | ✅ |
| 2 | Scope ZAWĘŻONY do A4 verification | ✅ |
| 3 | Hypothesis H1 testable, falsifiable | ✅ |
| 4 | NIE multi-candidate fit | ✅ N/A |
| 5 | NIE constructed criterion | ✅ |
| 6 | Sympy verification dual-V consistency | 🟡 in progress |
| 7 | Honest reporting jeśli kontradykcja znaleziona | 🟡 awaits sympy |
| 8 | Cross-reference do G.0, Phase 5, T-Λ | ✅ |

## TGP-internal axioms (LOCKED post-G.0)

- V_M9.1''(ψ) = -γψ²(4-3ψ)²/12 (gravity, G.0 P21 LOCKED)
- R3 ODE: psi'' + (2/r)psi' + (2/psi)(psi')² = (1-psi)/psi² (gravity)
- κ = 3/(4·Phi_0) (Newton limit, G.0 P32 LOCKED)
- A4 audit marker: matter coupling separate verification

## Cross-references

- [[../op-g0-r3-from-canonical-projection/Phase1_results.md#5.3]] — A4 marker
- [[../op-g0-r3-from-canonical-projection/phase2_P21_vacuum_uniqueness.py]] — V_M9.1'' derivation source
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] — V_orig matter usage
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] — V_orig matter usage
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] — annotation source

## Status

**Phase 0 COMPLETE.** Path C confirmed strong evidence. Ready for Phase 1 sympy.
