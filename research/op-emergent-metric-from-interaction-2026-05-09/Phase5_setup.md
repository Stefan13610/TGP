---
title: "Phase 5 setup — Lenz back-reaction (m_inertial), N9 + N10"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-setup
phase: 5
status: 🟢 OPEN
needs: [N9, N10]
predecessor: "[[./Phase4_results.md]] (8/8 PASS)"
---

# Phase 5 setup — Lenz back-reaction + equivalence principle

## §0 — Goal

Resolve N9 + N10 from [[./NEEDS.md]]:
- **N9**: Lenz back-reaction formalization → m_inertial obliczalny sympy
- **N10**: Zasada równoważności m_inertial = m_grav AUTOMATYCZNIE z S05

**Reference:** [[../../TGP_FOUNDATIONS.md]] §6 — Pęd i bezwładność (Lenz-podobny obraz):

> "Spoczynek = lokalna równowaga pola Φ wokół źródła ρ. Przyspieszenie
> źródła łamie tę równowagę → pole reaguje zwrotnie → siła back-reakcji
> opiera się zmianie konfiguracji → manifestuje się jako bezwładność."

## §1 — Strategy

### §1.1 — Linearized Φ-EOM with point source

Z S_TGP (FOUNDATIONS §3 canonical form):
```
S_TGP = ∫ d⁴x [-(K/2)(∂Φ)² - V(Φ) - q ρ Φ / Phi_0]
```

Linearize wokół vacuum Φ → Phi_0 + δΦ:
```
(∂_t² - c²∇² + m_sp²) δΦ = -q ρ / (Phi_0 K(1))
```

m_sp² > 0 (G.0 P21 LOCK), efektywnie cosmological scale → effectively massless
at solar-system distances.

### §1.2 — Static source solution (massless limit)

Dla ρ_static(x) = M δ³(x) i m_sp → 0:
```
Φ_eq(r) = Phi_0 - q M / (4π Phi_0 K(1) r)   = Phi_0 + δΦ_eq(r)
```

Newton-like 1/r tail. The amplitude q · M / (4π Phi_0 K(1)) defines
the "gravitational charge" of the source.

### §1.3 — Energy stored in static configuration

```
E_static = (1/2) ∫ d³x [K(1)(∇δΦ_eq)² + m_sp² δΦ_eq²]
```

For massless limit and point source:
```
E_static = (1/2) K(1) ∫ d³x (∇δΦ_eq)²
         = (q²M²) / (8π Phi_0² K(1)) · (1/r_min)   [diverges at point]
```

Renormalize: E_static = m_field · c² (mass-energy of field configuration).
This is the SELF-ENERGY of the static profile — it is the MASS associated
with the source field configuration.

### §1.4 — Boosted source (Galilean Lenz I — no back-reaction)

Source moves with constant velocity v: ρ(x - vt). Field follows Galilean
transformation:
```
Φ(x, t) = Phi_0 + δΦ_eq(x - vt)
```

This is a SOLUTION of the EOM (Galilean invariance for non-relativistic limit).
**No back-reaction** — uniform motion is preserved (Newton I from translational
invariance + minimum equilibrium configuration).

### §1.5 — Accelerating source (Lenz II — back-reaction = inertia)

Source accelerates with a: x_source(t) = X_0 + (1/2) a t². Now Φ ≠ Φ_eq translated;
field cannot follow infinitely fast → back-reaction.

Linear order in a: δΦ = δΦ_eq(x - x_source(t)) + δΦ_BR
where δΦ_BR satisfies (∂_t² - c²∇² + m_sp²) δΦ_BR = -ä · ∂_x δΦ_eq + O(a²)

Back-reaction force on source:
```
F_BR^i = -∫ d³x ρ · ∂^i Φ = -q M ∂^i Φ |_at source
       = -m_inertial · a^i
```

After regularization (analog of EM self-energy), m_inertial emerges
from the field configuration:
```
m_inertial = (1/c²) · E_static = (1/c²) · (q²M²)/(8π Phi_0² K(1)) · (1/r_reg)
```

### §1.6 — Equivalence principle

m_grav (1PN coupling, Phase 2 N4): coupling of source M to Φ-field
through L_mat ∝ q M Φ / Phi_0. Effective Newton-like constant:
```
G_eff = q² / (4π Phi_0² K(1))
m_grav = M  (source mass enters geodesic motion through M·δ³(x) coupling)
```

**Key claim:** m_inertial / m_grav = 1 if both derive from same q, M, K(1)
in the TGP single-field framework (S05 lock).

**Sympy verification:** show this ratio = 1 algebraically, without numerical
specifications.

## §2 — Phase 5 sympy strategy

```
1. Setup linearized Φ-EOM
2. Solve static point source → Φ_eq(r)
3. Compute m_grav from coupling structure (1PN identification)
4. Compute m_inertial from back-reaction integral (regularized)
5. Show ratio m_inertial / m_grav = 1 (sympy LOCK)
6. Verify Galilean invariance for constant velocity (Newton I)
7. Verify F_BR = -m·a structural form (Newton II)
```

## §3 — Anti-pattern compliance

- NIE post-hoc fitting m_inertial / m_grav = 1 (must emerge structurally)
- NIE drift hardening (do NOT introduce free parameters to make ratio = 1)
- Honest reporting: if regularization scheme matters, document explicitly

## §4 — Phase 5 deliverables

- Phase5_sympy.py
- Phase5_sympy.txt
- Phase5_results.md

## §5 — Gate criteria

| # | Criterion |
|---|---|
| G1 | Linearized Φ-EOM derived (sympy) |
| G2 | Φ_eq(r) static solution computed |
| G3 | m_grav structural form identified |
| G4 | m_inertial structural form identified (via back-reaction) |
| G5 | m_inertial / m_grav = 1 sympy LOCK |
| G6 | Newton I + II structural verification |

**Pass:** ≥5/6 with sympy verification.

## §6 — Cross-references

- [[../../TGP_FOUNDATIONS.md]] §6 — Lenz-podobny obraz (canonical TGP)
- [[./README.md]] — cycle overview
- [[./Phase4_results.md]] — predecessor (8/8 PASS)
- [[../op-g0-r3-from-canonical-projection/]] — V_M911 + m_sp² LOCK
