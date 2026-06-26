---
title: "Phase 2 plan — Frontier creation rate S_creation(t) derivation"
type: phase_plan
status: PRE_REGISTERED_LOCKED
phase: 2
parent_cycle: op-CE-H-gamma-3-cosmological-2026-05-23
pre_registration_date: 2026-05-23
---

# Phase 2 plan — Frontier creation rate S_creation(t) derivation

**Status:** PRE_REGISTERED_LOCKED 2026-05-23.
**Methodology:** Strict cycle 1/2/7 + §3.6.6-3.6.10 BINDING + native equations FIRST.
**Scope:** Derive functional form S_creation z (EQ-5) + frontier creation mechanism (concept paper §2.2-2.4); identify Hubble equation closure structure.

---

## §1 — Pre-registered FP (Fixed-point predictions)

**0 hardcoded T_pass=True.** Każdy FP = compute-then-compare.

### T_P2_1 — Frontier surface area per unit volume

**Hypothesis:** Dla spherically expanding cosmological region radius R(t):

$$\frac{A_{frontier}}{V_{universe}} = \frac{4\pi R^2}{(4/3)\pi R^3} = \frac{3}{R}$$

If H = c/R (geometric Hubble definition):

$$\frac{A}{V} = 3H/c = 3H \text{ (natural units c=1)}$$

**Pre-registered (literal):**
- A/V symbolic derivation explicit
- PASS: A/V = 3H (or 3H/c if c kept explicit)
- FAIL: NIE derive z geometry

### T_P2_2 — Per unit volume rate of new "saturated" volume

**Hypothesis:** Frontier moves at c (maximum signal speed in TGP, per Lorentz invariance of Phi-substrate); per unit cosmological volume rate of new "saturated" volume:

$$\dot{V}_{new}/V = (A/V) \cdot v_{frontier} = 3H \cdot c = 3H \text{ (natural)}$$

**Pre-registered (literal):**
- PASS: Symbolic derivation V̇_new/V = 3H
- FAIL: NIE derive

### T_P2_3 — S_creation = 3Hv (saturation efficiency η_sat = 1)

**Hypothesis:** Newly created volume saturates to ⟨Φ⟩ = v via σ-mode relaxation (timescale τ_σ = 1/m_σ). Pod assumption m_σ >> H (fast relaxation, late-time limit):

$$S_{creation} = (\dot{V}_{new}/V) \cdot v \cdot \eta_{sat} = 3H \cdot v$$

(η_sat = 1 exact saturation).

**Pre-registered (literal):**
- m_σ >> H limit verification (m_σ ~ 200 MeV >> H_0 ~ 10⁻³³ eV: factor ~ 10⁴¹) ✓
- Functional form S_creation = 3Hv derived
- PASS: S_creation symbolic explicit + η_sat = 1 justified
- FAIL: NIE derive OR factor of magnitude wrong

### T_P2_4 — Closure structure (EQ-5)+(EQ-6) and H determinacy

**Hypothesis:** Substituting S_creation = 3Hv into (EQ-5) stationary at ⟨Φ⟩ = v:

$$3H \cdot v = 3H \cdot v$$

**This is tautological.** H is NIE determined by (EQ-5)+(EQ-6) alone.

**Structural implication:** TGP-native cosmology has 1-parameter family parametrized by t_universe (initial condition / boundary condition); H_0 = 1/t_universe (geometric).

**Pre-registered (literal):**
- Tautology of (EQ-5)+(EQ-6) under stationary v explicit
- 1-parameter family identification
- PASS: structural insight documented + geometric H formula explicit
- PARTIAL: tautology confirmed but functional structure unclear
- FAIL: NIE tautological (i.e., system DOES determine H from v, λ alone)

**Anti-Lakatos honesty:** This finding pre-registered as one of 4 possible resolutions w Phase 0 §4.4 ("Frontier creation rate S_creation NIE jest direct v² scale (geometric factors)"). NIE post-hoc move.

---

## §2 — Cycle 1/2/7 budget

- 4 substantive FP pre-registered
- 0 hardcoded T_pass=True
- 1 PARTIAL allowed (T_P2_4 explicitly may be PARTIAL — pre-registered)

---

## §3 — DEC tracking

- **DEC 2 NIE expended w Phase 2** (numerical method deferred to Phase 3)
- Remaining: 2 of 3

---

## §4 — Honest disposition (pre-registered)

**Key structural finding expected:** TGP-native cosmology gives **geometric H_0 = 1/t_universe** in stationary E2 equilibrium.

**Consequences dla F-γ-3:**
- H_0 NIE pinned by (v, λ) alone — t_universe additional input
- Observed t_universe ≈ 14 Gyr → predicted H_0 ≈ 70 km/s/Mpc ✓ (in observed range [67, 73])
- But this is **PARTIAL DERIVATION** — t_universe is initial condition, not derived

**Pre-registered scenarios dla F-γ-3 verdict (Phase 3+):**
- **PASS-A**: Geometric H_0 formula derived; t_universe ≈ 14 Gyr independently observed → factor 2 PASS
- **PARTIAL**: Formula correct but t_universe not derived from TGP parameters → conditional PASS
- **FAIL**: Numerical mismatch z observed by factor > 2 → HALT-B

---

## §5 — §3.6 extension compliance

- §3.6.6 Sign: H > 0, ρ > 0, S > 0 (per Phase 0 §4.1)
- §3.6.7 DoF: Phase 2 derivation only (no fits)
- §3.6.8 Implicit: m_σ >> H late-time limit explicit
- §3.6.9 Precision: symbolic Phase 2; numerical Phase 3
- §3.6.10 Pattern monitoring: ongoing

---

## §6 — Pattern detection

Pre-registered: tautology finding might trigger R1 flag jeśli new methodology pattern emerges (e.g., "geometric vs functional H" issue NOT pre-anticipated in §3.6 extension).

If R1 flag → defer R2 audit cycle until Phase FINAL OR mid-cycle pause for user authorization.

---

## §7 — Computational plan (Phase2_sympy.py)

**Tool:** sympy symbolic.

**Sections:**
1. T_P2_1: A/V symbolic derivation (geometric)
2. T_P2_2: V̇/V derivation
3. T_P2_3: S_creation = 3Hv with η_sat = 1 justification
4. T_P2_4: (EQ-5)+(EQ-6) tautology check + geometric H formula

**Result file:** Phase2_results.md.

---

**Phase 2 plan LOCKED 2026-05-23. Ready dla Phase2_sympy.py execution.**
