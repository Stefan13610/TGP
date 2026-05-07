---
title: "Phase 1 setup — Φ-EOM w FRW background z varying c, ℏ, G"
date: 2026-05-06
parent: "[[README.md]]"
status: SETUP — gotowy do execution
phase: 1
verdict_gate: "≥4/5 PASS dla Phase 2 enable"
predecessor: "[[Phase0_balance.md]] (8/8 ☑ PASS)"
tags:
  - TGP
  - EXT-1
  - Phase1
  - FRW
  - Phi-EOM
  - varying-constants
---

# Phase 1 setup — Φ-EOM w FRW z varying c, ℏ, G

## Cel

Wyprowadzić **H_TGP(t, ψ(t))** z Friedmann equation w TGP z explicit
ax:c-ax:G integration. Sprawdzić czy Φ-EOM w FRW jest **well-defined**
w erze radiacyjnej (z > z_eq).

## Pre-requisites

- ✅ [[Phase0_balance.md]] 8/8 ☑ PASS (mandatory per Phase 6 gate)
- ✅ L01 formal definition EXECUTED 2026-05-04
- ✅ ax:c-ax:G axiomy w sek04_stale.tex
- ✅ M9.1'' canonical form (z S07 caveat — trzecia iteracja)
- ✅ closure_2026-04-26 T-Λ closure

## Plan sub-tests (5 tests, gate ≥4/5)

### F1.1 — FRW background setup z ψ(t)

**Cel:** Setup metryki FRW w TGP z explicit ψ(t) field evolution.

**Setup:**
```
Metric: ds² = -c(Φ)² dt² + a(t)²[dr² + r²(dθ² + sin²θ dφ²)]
       = -c_0²·(Φ_0/Φ)·dt² + a(t)²·[...]   (z ax:c)

ψ(t) field evolution: cosmological scale background
Substrate field: Φ(t) = Φ_0·ψ(t)        (definicja w TGP)
```

**Test:**
- Sympy LOCK: czy FRW metryka z varying c(Φ) jest dobrze zdefiniowana
  (signature -+++, smooth)?
- Czy proper time τ_proper = ∫dt·c(Φ(t))/c_0 jest spójne?

**Expected:** PASS (FRW jest standard background, varying c jest
dimensional rescaling).

### F1.2 — Friedmann equation w TGP

**Cel:** Wyprowadzić H_TGP² = f(ρ, ψ, dψ/dt, V(ψ)) z action
TGP wariacja po a(t).

**Setup:**
```
S_TGP = ∫ √(-g) [R/(16πG(Φ)) + ½K(ψ)·g^μν ∂_μψ ∂_νψ - V(ψ) - L_mat] d⁴x

Variation δS/δg_tt → Friedmann eq w TGP:
H_TGP² = (8πG(Φ)/3) [ρ_matter + ½K(ψ)·(dψ/dt)²·c(Φ)⁻² + V(ψ)/c(Φ)²]
       + (curvature terms if k ≠ 0; we assume k = 0 flat)

Substytuja ax:c-ax:G:
G(Φ) = G_0·Φ_0/Φ
c(Φ) = c_0·√(Φ_0/Φ)

H_TGP² = (8πG_0/3)·(Φ_0/Φ) · [ρ_matter + ½K(ψ)·(dψ/dt)²·(Φ/Φ_0)/c_0² + V(ψ)·(Φ/Φ_0)/c_0²]
```

**Test:**
- Sympy LOCK Friedmann eq w TGP (sympy sympy.diff variational)
- Recovery w limicie ψ → 1: H_TGP² → H_GR² standard?
- Cross-check vs M10.1 (FRW DE w(z)) sub-cycle z M10 cosmology

**Expected:** PASS (Φ-EOM w FRW jest standard derivation; recovery
w ψ=1 limit jest necessary check).

### F1.3 — ax:c-ax:G integration

**Cel:** Substytuować ax:c, ax:ℏ, ax:G w wszystkich pojawiających się
wyrażeniach w Phase 1 i sprawdzić consistency.

**Setup:**
```
- ρ_matter = ρ_rest_mass·c(Φ)²/c_0² (energy density z mass)
- Energy momentum T^μν komponentu materii
- Quantum nuclear reactions (Saha, Coulomb barrier) — varying ℏ effects
```

**Test:**
- Sympy substitution wszędzie gdzie c, ℏ, G się pojawia
- Sprawdzenie czy redefiniowane "stałe" są konsistentne dimensionally
  (units check)
- Recovery w obecnej epoce: c(Φ_0) = c_0, ℏ(Φ_0) = ℏ_0, G(Φ_0) = G_0

**Expected:** PASS (axiomy są pre-defined; substitution jest mechaniczna).

### F1.4 — Limity asymptotyczne (radiation-dominated era)

**Cel:** W limicie z >> z_eq (radiation-dominated GR era), znaleźć ψ(z)
asymptotic behavior.

**Setup:**
```
GR: H_GR(z) ∝ a(z)⁻² ∝ (1+z)² (radiation-dominated)

TGP: H_TGP(z) = ? z Φ-EOM w FRW

Question: czy istnieje sensowne asymptotic ψ(z) takie że H_TGP(z)
jest finite + smooth dla z → ∞?
```

**Test:**
- Asymptotic analysis Φ-EOM dla z >> z_eq
- Sprawdzenie singularności (Big Bang singularity w TGP?)
- Sprawdzenie czy ψ(z=10¹⁰) jest realistic (np. ψ ∈ [0.1, 10] albo
  ψ → 0 lub ∞)

**KRYTYCZNY:** jeśli ψ(z) → 0 lub ∞ w erze radiacyjnej, M9.1''
założenia perturbacyjne wokół ψ=1 łamią się. To jest **R1** risk
z Phase0_balance §7.

**Expected:** PASS (TGP ma sensowne asymptotic) lub FAIL (singularność
→ Phase 1 FAIL → cykl ABANDONED → ścieżka E).

### F1.5 — Phase 1 GATE ≥4/5 PASS

**Cel:** Aggregated gate criterion.

**Test:**
- F1.1+F1.2+F1.3+F1.4 sum
- Minimum 4/5 PASS dla Phase 2 enable

**Decision:**
- 5/5 PASS → Phase 2 ENABLED, optimistic outlook
- 4/5 PASS → Phase 2 ENABLED with caveat
- 3/5 PASS → Phase 1 FAILED, considera ścieżki D/E
- ≤2/5 PASS → Cykl ABANDONED, raport do user-a

## Materiał do execution

**Skrypty proponowane (do napisania):**

```
phase1_FRW_setup.py            — F1.1 FRW metric setup
phase1_friedmann_TGP.py         — F1.2 Friedmann eq derivation
phase1_axcGh_substitution.py    — F1.3 ax:c-ax:G integration
phase1_asymptotic_radiation.py  — F1.4 asymptotic z >> z_eq
phase1_gate_summary.py          — F1.5 PASS/FAIL aggregator
```

**Tools:**
- sympy (variational, asymptotic analysis)
- scipy.integrate (ODE solver dla Φ-EOM jeśli numerical)
- numpy (asymptotic checks)

## Cross-references

- [[Phase0_balance.md]] — pre-Phase-1 mandatory balance sheet
- [[README.md]] — overall program
- [[../../core/sek04_stale/sek04_stale.tex]] — ax:c-ax:G axiomy
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
  — TGP action
- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]] — sek05 cosmology
- [[../op-cosmology-closure/M10_1_results.md]] — M10.1 FRW DE w(z)
- [[../op-cosmology-closure/M10_3_results.md]] — M10.3 FRW propagator
- [[../closure_2026-04-26/]] — T-Λ closure source
