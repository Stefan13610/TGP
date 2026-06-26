---
title: "Phase 1 plan — Cosmological ansatz + (EQ-1)-(EQ-6) self-consistency setup"
type: phase_plan
status: PRE_REGISTERED_LOCKED
phase: 1
parent_cycle: op-CE-H-gamma-3-cosmological-2026-05-23
pre_registration_date: 2026-05-23
authorization: "User explicit: 'Authorize γ-3 Phase 1+ multi-session (full commitment)'"
---

# Phase 1 plan — Cosmological ansatz + (EQ-1)-(EQ-6) self-consistency setup

**Status:** PRE_REGISTERED_LOCKED 2026-05-23.
**Methodology:** Strict cycle 1/2/7 + §3.6.6-3.6.10 BINDING + native equations FIRST (NIE Friedmann a priori).
**Scope:** Apply spatially homogeneous isotropic ansatz do (EQ-1)-(EQ-6); derive emergent equations form; test (EQ-3) self-consistency at cosmological scale; identify mapping do Friedmann post-derivation (NIE input).

---

## §1 — DEC declarations (pre-registered)

### DEC 1 (BUDGETED, per README §5)

**Decision:** Cosmological symmetry ansatz = **spatially homogeneous + isotropic** (Friedmann-Robertson-Walker-like emergent).

**Justification:**
- Concept paper §3.4: "NIE zakłada Robertson-Walker symmetry a priori (ma wynikać z self-consistency, nie być przyjęta)"
- W Phase 1: ansatz jest **test hypothesis**, NIE assumption: czy spatially homogeneous + isotropic Phi-substrate JEST self-consistent solution (EQ-3)?
- Alternative anisotropic ansatz (Bianchi-type) reserved dla DEC 3 jeśli homogeneous fails self-consistency

**DEC 1 expenditure:** 1 of 3 pre-registered DEC budget. Remaining: 2.

---

## §2 — Pre-registered FP (Fixed-point predictions) — strict cycle 1/2/7

**0 hardcoded T_pass=True.** Każdy FP = compute-then-compare against pre-registered threshold.

### T_P1_1 — Spatial homogeneity ansatz applied do (EQ-2)

**Hypothesis:** Pod ansatzem ⟨Φ⟩(x,t) = ⟨Φ⟩(t) (spatial homogeneity), (EQ-2) reduces do integral form:

$$\langle\Phi\rangle(t) = \langle\Phi\rangle_{frontier}(t) + \int d^3 x' \, G(\vec{x} - \vec{x}')\, \rho_{avg}(t)$$

gdzie ρ_avg(t) = uniform mean density.

**Pre-registered (strict literal):**
- Spatial derivative ∂⟨Φ⟩/∂x = 0 musi być **konsystentny** z (EQ-2) right-hand side spatial dependence
- Po homogeneity, równanie redukuje się do **algebraicznego** związku ⟨Φ⟩(t) jako function ρ_avg(t) + boundary
- PASS jeśli: redukcja explicit + RHS spatial dependence eliminuje się pod uniform ρ assumption
- FAIL jeśli: residual spatial term NIE eliminuje się → spatially homogeneous solution NIE istnieje

**Anti-Lakatos:** PASS threshold = symbolic identity check (LHS ≡ RHS spatially uniform).

### T_P1_2 — (EQ-5) cosmological ODE form derivation

**Hypothesis:** Pod ansatzem, (EQ-5) staje się **pierwszego rzędu ODE** dla ⟨Φ⟩(t):

$$\frac{d\langle\Phi\rangle}{dt} + 3H(t)\langle\Phi\rangle = S_{creation}(t)$$

gdzie H(t) jest external function (do be derived w (EQ-6)).

**Pre-registered (strict literal):**
- ODE structural form: linear pierwszego rzędu z **damping term** 3H⟨Φ⟩ + source S_creation
- Stationary fixed point: d⟨Φ⟩/dt = 0 → ⟨Φ⟩_∞ = S_creation / (3H)
- PASS jeśli: ODE symbolic derivation + stationary FP solution explicit
- FAIL jeśli: form nie zgadza się z concept paper §4.2 (EQ-5)

**Anti-Lakatos:** PASS = symbolic verification stationary FP równanie.

### T_P1_3 — (EQ-3) self-consistency at cosmological scale

**Hypothesis:** Self-consistency (EQ-1) ∘ (EQ-2) z ansatzem homogeneous → fixed-point ⟨Φ⟩*(t) such that:

$$\rho_i^* = \rho_i^*[\langle\Phi\rangle^*(t)] \quad \text{AND} \quad \langle\Phi\rangle^*(t) = \langle\Phi\rangle_{frontier}(t) + G^{(0)}\int d^3x' \rho^*(t)$$

gdzie G^(0) = ∫G(x-x')d³x' jest spatial integral Green function.

**Pre-registered (strict literal):**
- Fixed point musi **istnieć** dla rozumnych physical parameter values (NIE diverging)
- W limicie ρ → ρ_avg uniform, G^(0) musi być finite (NIE infinity per Green function regularization)
- PASS jeśli: G^(0) symbolic finite OR explicit IR regularization scheme
- FAIL jeśli: G^(0) diverging → spatially homogeneous Phi-substrate nie ma fixed point → DEC 3 anisotropic alternative wymagana

**Anti-Lakatos:** PASS = symbolic finiteness check G^(0).

### T_P1_4 — (EQ-6) Hubble equation functional form

**Hypothesis:** (EQ-6) H² = H[ρ_i, ⟨Φ⟩, S_creation] pod ansatzem reduces do:

$$H^2 = f(\rho_{avg}, \langle\Phi\rangle, S_{creation}, \lambda, v)$$

gdzie f jest **symbolic function** TGP parameters.

**Pre-registered (strict literal):**
- Functional form **derived** (NIE postulated)
- Late-time limit (z << 1, per Phase 0 §4.3 limit choice) gives algebraic relation H² ~ kombinacja parameter
- Mapping do Friedmann form H² = (8πG/3)ρ_total: **post-derivation comparison** (NIE input)
- PASS jeśli: f symbolic derivation explicit
- PARTIAL jeśli: derivation requires DEC 2 (numerical method) — defer to Phase 2
- FAIL jeśli: NIE derive z (EQ-1)-(EQ-5) under ansatz → fundamental gap

**Anti-Lakatos:** PASS = symbolic functional form, NIE Friedmann substitution. Mapping post-derivation comparison only.

---

## §3 — Cycle 1/2/7 budget

- 4 substantive FP pre-registered
- 0 hardcoded T_pass=True
- 1 PARTIAL allowed (T_P1_4 may be PARTIAL z odsyłaniem do Phase 2 numerical)
- 0-1 HONEST FAIL acceptable

---

## §4 — §3.6 extension compliance

### §3.6.6 Sign convention
- H > 0 expansion (Phase 0 §4.1)
- ρ > 0 matter (Phase 0 §4.1)
- d⟨Φ⟩/dt sign dependent on whether bulk underdense lub saturated (TBD Phase 2)

### §3.6.7 Fit DoF equalization
- Phase 1 jest **derivation only**, NIE fitting
- Mapping comparisons Phase 3+ z parameter count

### §3.6.8 Implicit assumptions
- Cosmological principle = TESTED w T_P1_3, NIE assumed
- Linearization around E2 equilibrium (Phase 0 §4.3)
- Late-time z << 1 limit
- Mean-field Hartree-Fock-like

### §3.6.9 Numerical precision
- Phase 1 = symbolic; Phase 3+ numerical
- Precision validation deferred do Phase 3

### §3.6.10 Methodology evolution
- Niniejszy cykl = third practical application §3.6 extension (post γ-1 retry + Phase 0)
- Pattern detection ongoing

---

## §5 — Forbidden moves (inherited from README §2)

1-10 standard anti-Lakatos
11. NIE ΛCDM fitting (native-equations-first)
12. NIE Hoyle-Bondi steady-state reuse
13. NIE ad-hoc Λ cosmological constant

---

## §6 — Computational plan (Phase1_sympy.py)

**Tool:** sympy symbolic.

**Sections:**
1. Symbol declarations (t, x, ⟨Φ⟩(t), ρ_avg(t), H(t), S_creation(t), λ, v, G_0)
2. T_P1_1 spatial homogeneity reduction: apply ⟨Φ⟩(x,t) → ⟨Φ⟩(t) do (EQ-2)
3. T_P1_2 (EQ-5) ODE form: stationary FP derivation
4. T_P1_3 (EQ-3) self-consistency: G_0 finiteness check + cosmological FP existence
5. T_P1_4 (EQ-6) functional form derivation z (EQ-1)+(EQ-2)+(EQ-5)

**Output:** Pass/Fail per FP per literal threshold.

**Result file:** Phase1_results.md.

---

## §7 — Cross-validation (Path A-D per Phase 0 §8)

- **Path A analytical:** Phase 1 symbolic derivation
- **Path B numerical:** Deferred Phase 2-3
- **Path C limit checks:** Late-time z << 1 within Phase 1; weak-coupling limit Phase 2
- **Path D observational anchors:** Deferred Phase 3+

---

## §8 — Risk register Phase 1 specific

| ID | Risk | Mitigation | Severity |
|----|------|-----------|----------|
| RP1-1 | G_0 spatial integral diverges (R3 IR divergence) | T_P1_3 explicit check; honest FAIL acceptable | MEDIUM |
| RP1-2 | Anisotropic Bianchi ansatz wymagany (homogeneous fails) | DEC 3 reserve; honest pivot acceptable | MEDIUM |
| RP1-3 | (EQ-6) NIE derives symbolically (calculational hell) | PARTIAL acceptable z odsyłaniem do Phase 2 | HIGH |
| RP1-4 | ΛCDM accidental ansatz drift | Anti-Lakatos forbidden #11 LOCK; native equations FIRST | LOW |

---

**Phase 1 plan LOCKED 2026-05-23. Ready for Phase1_sympy.py execution.**
