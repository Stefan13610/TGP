---
title: "ξ.1.Phase2 setup — heat-kernel a₂ first-principles derivation"
date: 2026-04-29
cycle: ξ.1.Phase2
status: PRE-EXECUTION
predecessor: "[[Phase1_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - xi-factor
  - heat-kernel
  - a2-coefficient
  - photon-ring
  - falsification
---

# ξ.1.Phase2 — Heat-kernel a₂ first-principles derivation

> **Cel:** Wyprowadzić target_shift z heat-kernel a₂ coefficient pod
> TGP substrate inputs (F1 + F2 + F3 + F4). Decide czy F4 rational
> (Frame A), Phase 2 strict (Frame B), albo trzecia forma jest true
> a₂-derived value. Na PASS: ξ.1 promotion z STRUCTURAL POSTULATE → DERIVED
> lub PARTIALLY DERIVED.

---

## Hipoteza Phase 2 (H_ξ2)

Pod TGP substrate (F1 single-Φ ∧ F2 K_geo·φ⁴ ∧ F3 φ_eq=1) na M9.1″ FRW
background, heat-kernel a₂ Avramidi/Birrell-Davies expansion **exclusively
selects** jeden z trzech frames:

```
Frame A:   target_shift = 0.114      (F4 sympy rational)
Frame B:   target_shift = 11/97      (Phase 2 strict)
Frame C:   target_shift = ξ-degenerate w O(1%) precision
```

a₂ struktura na vacuum (M9.1″ Ricci R suppressed 10⁻¹²² dex):
```
a₂ ⊃ (1/2) V''(Φ_eq)² + (1/180) R² − (1/180) R_μν R^μν + …
   = (1/2) V''(1)²  (Ricci negligible na FRW)
```

V''(Φ_eq) z F2 lock K(φ) = K_geo·φ⁴, F3 lock φ_eq = 1, normalization F4
α₀ anchor = compute V''(1) → a₂(1) → target_shift via heat-kernel
expansion z photon-ring background.

---

## 7 sub-testów Phase 2

### ξ2.1 — Heat-kernel a₂ structure z Birrell-Davies / Avramidi

**Cel:** Skonstruować symbolicznie a₂ coefficient na FRW background w
single-scalar-field theory.

**Test:**
- Standard form (Birrell-Davies 1982, Eq. 6.45):
  ```
  a₂(x) = (1/2) m²(x)² + (1/6) R · m²(x) + (1/180)(R² − R_μν R^μν)
  ```
- Single-Φ + K(φ)·(∂Φ)²: m²(Φ) → V''(Φ)
- Pod F2: V''(Φ) jest wynikiem F2 + F3 substrate stability conditions

**Falsification:** a₂ structure deviates z standard convention (e.g.,
DeWitt-Schwinger sign ambiguities) → frame nie zamknięty.

### ξ2.2 — TGP V''(Φ) evaluation pod F2 + F3

**Cel:** Wyliczyć V''(Φ_eq=1) symbolicznie z TGP substrate Lagrangian.

**Test:**
- Phase 1.A: K(φ) = K_geo · φ⁴ (Theorem alpha2)
- Phase 1.F.3: Φ_eq = 1, V'(Φ_eq) = 0 exact
- Effective potential z Z₂: V(Φ) = (β/4) Φ⁴ − (γ/2) Φ² + const
- Pod β = γ vacuum: V''(Φ_eq) = (3β·Φ²−γ)|_{Φ=1} = 3β − γ = 2β = 2γ
- Pod F4 anchor: β = γ = wartość z α₀ chain

**Falsification:** V''(Φ_eq) ma free parameter not LOCKED.

### ξ2.3 — M9.1″ FRW Ricci suppression

**Cel:** Confirm R, R_μν R^μν contributions do a₂ są subdominant na
vacuum (M9.1″ FRW background).

**Test:**
- M9.1″ FRW: R = 12·H₀² ≈ 12 · (10⁻⁶¹ Pl)² ~ 10⁻¹²² Pl⁻²
- R²/m_Pl⁴ ~ 10⁻²⁴⁴ Pl⁰ — strongly negligible
- (1/2) V''(1)² dominuje w a₂(1) na vacuum

**Falsification:** Ricci R term NOT subdominant w a₂ → background-dependent ξ.

### ξ2.4 — a₂ → target_shift conversion via photon-ring back-reaction

**Cel:** Wykonać conversion z a₂(1) na target_shift dla photon-ring.

**Test:**
- Heat-kernel formula: target_shift = ξ_geom · α(α−1) · a₂_normalized_ratio
- a₂_normalized_ratio = (V''(1)²/2) / [reference scale]
- Reference scale: M9.1″ photon-ring perturbation amplitude
- ξ_geom = 1 (ξ1.1 LOCKED)
- α(α−1) = 2 (ξ1.2 LOCKED)
- target_shift = 2 · a₂_ratio

**Falsification:** conversion has free O(1) parameter not pinnable
z LOCKED inputs.

### ξ2.5 — Frame A test: a₂ yields F4 sympy rational

**Cel:** Sprawdzić, czy a₂ derivation z TGP inputs daje target_shift
zgodne z F4 sympy rational 1069833/264500 / [(0.16788)² · 1.0] ≈ 4.04474.

**Test:**
- Pod β = γ = α₀-derived vacuum choice
- V''(1) = 2β such that V''(1)²/2 = 2β² yields:
  - target_shift_A = ξ_geom · α(α−1) · (2β² / reference)
- Frame A holds if target_shift_A = 0.114 (within 0.1% precision)

**Falsification:** a₂ NOT 4.0447 sympy → Frame A FALSIFIED.

### ξ2.6 — Frame B test: a₂ yields Phase 2 strict

**Cel:** Sprawdzić, czy alternative frame yields target_shift = 11/97
≈ 0.11340.

**Test:**
- Pod alternate normalization (e.g., bare-form direct geometric)
- target_shift_B = (1/2)(1 − r_ph^GR/r_ph^TGP) = 11/97
- Frame B holds if a₂ derivation matches strict bare form (within 0.1%)

**Falsification:** a₂ NOT 4.0179 → Frame B FALSIFIED.

### ξ2.7 — Classification decision

**Cel:** Skompletować algebraic chain i przyporządkować klasyfikację.

**Test:**
- Pod consistency wszystkich 6 LOCKED inputs (F1, F2, F3, F4, ξ_geom,
  α(α−1)) + a₂ derivation pod canonical normalization
- Status:
  - **DERIVED**: jeden frame uniquely selected, sub-percent (<0.1%) match
  - **PARTIALLY DERIVED**: oba frames consistent within EFT 1% precision (frame degeneracy resolved by UV completion)
  - **STRUCTURAL HINT**: oba frames possible, no a₂-unique selection

**Falsification:** chain has free parameter NOT pinnable z TGP inputs.

---

## Verdict gate

**≥ 6/7 PASS** → ξ.1 DERIVED or PARTIALLY DERIVED, proceed Phase 3.

**≤ 5/7 PASS** → ξ.1 STRUCTURAL HINT, Phase 3 still proceeds (limited
registry update).

---

## Materiał wykonawczy

- **Skrypt:** [`phase2_a2_derivation.py`](phase2_a2_derivation.py) (7 sub-tests, sympy-symbolic + numeric)
- **Output:** `phase2_a2_derivation.txt`
- **Memo:** `Phase2_results.md`

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 TGP/TGP_v1/research/op-xi-photon-ring/phase2_a2_derivation.py 2>&1 | tee TGP/TGP_v1/research/op-xi-photon-ring/phase2_a2_derivation.txt
```

---

## Constants used

```
xi_geom              = 1                        (ξ1.1)
alpha_K              = 2                        (ξ1.2)
alpha_K_x_alpha_K_m_1 = 2                        (= 2·1)
psi_ph               = 1.16788                  (ξ1.3, M9.1″ refined)
eps_ph               = 0.16788                  (= ψ_ph - 1)
eps_ph_sq            = 0.0281836944             (M9.1″ refined)
target_shift_F4      = 0.114                    (F4 anchor)
target_shift_strict  = 11/97 ≈ 0.11340206       (sympy exact)
F4_rational          = 1069833/264500 ≈ 4.04474
alpha_0_strict       = 4.01793
phi_eq               = 1                        (F3)
K_geo                = 1                        (F2)
beta = gamma         = (α₀-derived vacuum)      (F4 chain)

a2_BD_pure_vacuum    = (1/2) V''(Φ_eq)²
Ricci_suppression    = 10⁻²⁴⁴ Pl⁰ (negligible)

split_F4_strict      = 0.527%                   (Phase 1 quantified)
EFT_precision_band   = 1%                       (a₂ EFT precision)
```

---

## Cross-references

- [`Phase1_results.md`](Phase1_results.md) — premise audit 5/5 PASS
- [`program.md`](program.md) — overall ξ.1 plan
- [`../op-cross-sector-charge/Phase2_results.md`](../op-cross-sector-charge/Phase2_results.md) — XS.1 PARTIALLY DERIVED
- [`../op-phase3-uv-completion/Phase3_E_results.md`](../op-phase3-uv-completion/Phase3_E_results.md) — UV7 a₂ frame
- [`../op-phase2-quantum-gravity/Phase2_B_results.md`](../op-phase2-quantum-gravity/Phase2_B_results.md) — F4 chain
