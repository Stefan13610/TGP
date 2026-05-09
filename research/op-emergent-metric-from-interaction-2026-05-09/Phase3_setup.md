---
title: "Phase 3 setup — 2.5PN binary inspiral β_ppE^new via σ-cross-terms"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-setup
phase: 3
status: 🟢 OPEN
needs: [N6, N7, N8]
predecessor: "[[./Phase2_results.md]] (7/7 sympy PASS, σ-coupling free at 2PN)"
---

# Phase 3 setup — 2.5PN binary inspiral, β_ppE^new derivation

## §0 — Goal

Resolve N6, N7, N8 from [[./NEEDS.md]]:
- **N6**: 2-source case formalization (gradient cross-terms σ_cross)
- **N7**: Effective phase modification δφ(f) via SPA chain
- **N8**: β_ppE^new vs M9.1'' single-source -15/4

**Strategic context:** Phase 2 established γ=β=1 EXACT (1PN/2PN) with
**σ-coupling C(ψ) FREE in solar-system regime.** The leading coefficient
c_0 = C(1) enters at 2.5PN binary inspiral via gradient cross-terms
∂_μΦ_1·∂_νΦ_2 — *structurally absent* in single-source M9.1''.

This Phase 3 derives β_ppE^new(c_0): a parametric formula for the ppE
coefficient as function of c_0 (and other free coupling parameters).

## §1 — Phase 3 strategy

### §1.1 — Decomposition of contributions to β_ppE^new

Using SPA chain (G_SPA = 48 sympy-exact LOCK from [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]]):

```
β_ppE^(b=-1) = -(3/(128η)) · Δα_4

Δα_4 = (Δα_4)_diag + (Δα_4)_σ-cross
     [single-source]   [2-source, NEW]
```

**M9.1'' (single-source, c_0 = 0):**
- (Δα_4)_diag = -40 (from b_n single-source coefficients) — **FALSIFIED 5.02σ**
- (Δα_4)_σ-cross = 0 (c_0 = 0 turns off σ contribution)
- β_ppE^M911 = -(3/32)·(-40) = 15/4 = 3.75

**Phase 3 derivation:**
- (Δα_4)_diag depends on (a_n, b_n) Taylor coefs of A, B
- (Δα_4)_σ-cross depends on c_0 + 2-body geometry
- β_ppE^new = β_ppE_diag + κ_σ · c_0 (linear in c_0 leading order)

### §1.2 — Required derivation steps (per N6/N7/N8)

| Step | Sub-task | Sympy effort |
|---|---|---|
| §3.1 | σ_ij^cross in 2-body geometry | Medium |
| §3.2 | g_eff^ij correction from σ-coupling | Light |
| §3.3 | Binding energy E_orb modification at 2PN | Heavy |
| §3.4 | Radiation flux F(v) check (still GR-form?) | Light |
| §3.5 | SPA chain → α_4 with σ contribution | Medium |
| §3.6 | β_ppE^new(c_0) parametric formula | Light |

## §2 — N6 setup: 2-source binary geometry

### §2.1 — Configuration

Binary in COM frame with separation r_12, mass ratio η = m_1·m_2/(m_1+m_2)²:
- Particle 1 at x_1 = -m_2 r_12 / M_total, x_2 = +m_1 r_12 / M_total
- For equal mass η = 1/4: x_1 = -r_12/2, x_2 = +r_12/2

Newtonian-potential field:
```
δΦ_i(x) = -G M_i / |x - x_i|
∂_j δΦ_i = G M_i (x - x_i)_j / |x - x_i|³
```

### §2.2 — σ_ij = (∂_iΦ)(∂_jΦ) - (1/3)δ_ij(∇Φ)²

Cross-terms (the new contribution):
```
σ_ij^cross = (∂_iΦ_1)(∂_jΦ_2) + (∂_iΦ_2)(∂_jΦ_1) - (2/3)δ_ij(∇Φ_1·∇Φ_2)
```

This is **traceless 3D tensor** with anisotropy along the separation
axis. Magnitude order: ∂_jΦ_i ~ GM_i/r_i², so σ_ij^cross ~ G²M_1M_2/(r_1²r_2²).

### §2.3 — At particle position (test particle limit)

Evaluate σ_ij^cross at x = x_1 (position of particle 1):
- r_1 → 0 (singular self-term — handled via regularization)
- r_2 = r_12 (separation)
- Contribution from ∂_iΦ_2 evaluated at x_1: ~ GM_2 / r_12²

Cross product (∂_iΦ_1·∂_jΦ_2)|_{x=x_1} requires careful regularization
(Hadamard partie-finie or dimensional reg).

## §3 — N7 setup: Phase modification via SPA

### §3.1 — Modified circular orbit binding energy

E_orb = -(η v²/2) · [1 + e_1 v² + e_2 v⁴ + ...]

In modified g_eff with σ-coupling, e_2 receives correction:
```
e_2^new = e_2^GR + Δe_2^diag + Δe_2^σ(c_0)
       = -27/8 + Δe_2^diag + (function of c_0, η)
```

### §3.2 — Radiation flux F(v) — GR-form preserved

Per Phase 2 N4c: σ-coupling at O(σ²) only at 3PN+ in flux. So
F(v) = F_GR(v) at leading 2PN-orbital. (GW1 + GW2 cross-channel
consistency from Phase 1.5 still holds.)

### §3.3 — SPA chain → α_4

α_4 = 30·e_2 - 20·e_1·p_1 + 10·p_1² - 10·p_2

With e_2 modified by σ-coupling: Δα_4 has contribution 30·Δe_2^σ(c_0).

```
β_ppE^new = -(3/(128η)) · [Δα_4_diag + 30·Δe_2^σ(c_0)]

At η=1/4:
β_ppE^new = -(3/32) · [Δα_4_diag + 30·Δe_2^σ(c_0)]
```

## §4 — N8 setup: comparison with M9.1''

### §4.1 — Recovery of single-source limit

In single-source limit (one mass, no binary): σ_ij^cross → 0 trivially
(no second source). Then:
- β_ppE^new → β_ppE_diag (no σ contribution)
- If A=ψ/(4-3ψ), B=(4-3ψ)/ψ, c_0 free → recovers M9.1'' single-source -15/4

This is the BACKWARDS COMPATIBILITY: single-source M9.1'' is special case.

### §4.2 — Binary deviation

For binary η=1/4, σ-cross-term contribution is structurally NEW:
- M9.1'' didn't have σ-coupling (c_0 = 0 implicitly)
- New ansatz allows c_0 ≠ 0 → genuine 2-body modification

### §4.3 — Structural relation

```
β_ppE^new(c_0) = β_ppE_diag + κ_σ(η) · c_0

where κ_σ(η) is the "structural sensitivity" of β to σ-coupling.
```

## §5 — Phase 3 deliverables

- Phase3_sympy.py — full derivation
- Phase3_sympy.txt — output with X/Y PASS
- Phase3_results.md — synthesis + β_ppE^new(c_0) formula

## §6 — Phase 3 gate criteria

| # | Criterion |
|---|---|
| G1 | σ_ij^cross derived in 2-body sympy |
| G2 | (Δα_4)_σ-cross derived via SPA chain |
| G3 | β_ppE^new(c_0) parametric formula |
| G4 | Single-source recovery: c_0=0 → β = β_diag |
| G5 | Verification: M9.1'' canonical recovers (a, b) → β = -15/4 |
| G6 | c_0 status documented (derivable / free / framework-fixed) |

**Pass:** ≥5/6 with sympy verification.

## §7 — Critical question: c_0 status

**THE central question of Phase 3:**

Is c_0 (leading σ-coupling coefficient) **derivable from TGP framework**,
or is it a **free parameter** to be fitted?

| Option | Implication |
|---|---|
| (A) Derivable from σ_ab (OP-7 T2) coupling structure | Best — first-principles derivation |
| (B) Constrained by SU(2) cross-consistency (N11) | Good — programmatic unification |
| (C) Free parameter | STRUCTURAL_CONDITIONAL — reduces TGP predictivity |
| (D) Fixed by Φ_0 EFT scale-dependence | Connects to op-Phi-vacuum-scale findings |

This question must be addressed in Phase 3 §6 (G6 criterion).

## §8 — Anti-pattern compliance

- NIE multi-c_0 fit z minimum drift (CALIBRATION_PROTOCOL anti-pattern 1)
- NIE drift hardening (anti-pattern 3): if c_0 must be specific value to
  pass GWTC-3, this MUST be documented as TENTATIVE (pending §7 resolution)
- Honest reporting: structural form first, fitting only if §7 resolves to (C)
