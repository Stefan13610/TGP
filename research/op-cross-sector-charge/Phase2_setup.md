---
title: "XS.1.Phase2 setup — substrate-action derivation of √α₀ = κ_TGP"
date: 2026-04-28
cycle: XS.1.Phase2
status: PRE-EXECUTION
predecessor: "[[Phase1_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - cross-sector
  - alpha-0
  - kappa-TGP
  - substrate-action
  - falsification
---

# XS.1.Phase2 — Substrate-action derivation of √α₀ = κ_TGP

> **Cel:** Wyprowadzić identyczność √α₀ = κ_TGP z **wspólnego
> substrate-field action** S[Φ, g_μν, J^μ, ψ_e] pod TGP single-Φ axiom
> (F1) + Z₂ vacuum (F3) + K(φ) = K_geo·φ⁴ (F2).

---

## Hipoteza Phase 2 (H_xs)

Pod TGP single-Φ axiom, oba sektory (BH photon-ring T·J·J, SC pair-
breaking λ_sf) mają **bilinearną strukturę** względem substrate Φ. Pod
canonical normalization (g_TGP = uniwersalna "TGP charge"), Wilson
coefficient każdego operatora to:

```
W_BH (photon-ring T·J·J) = g_TGP² · M_BH(K(φ_eq), φ_eq, geometry)
W_SC (pair-breaking)     = g_TGP² · M_SC(K(φ_eq), φ_eq, geometry)
```

Pod TGP single-Φ + Z₂ + β=γ vacuum (φ_eq = 1):
- **K(φ_eq) = K_geo = 1** (F2 lock)
- **φ_eq = 1** (F3 lock)
- **g_TGP = uniwersalne** (single-Φ — jedna stała sprzężenia substrate)

Pod canonical normalization geometric matrix elements **konwergują**:
M_BH = M_SC = M (universal substrate matrix element).

Stąd:
```
α₀ ≡ W_BH = g_TGP² · M
κ_TGP² ≡ W_SC = g_TGP² · M
→ α₀ = κ_TGP²   ↔   √α₀ = κ_TGP
```

---

## 7 sub-testów Phase 2

### T2.1 — TGP action operator content

**Cel:** Zidentyfikować minimalny zestaw operatorów Φ-bilinear w TGP
Lagrangian, który generuje **oba** sektory (BH T·J·J, SC pair-breaking)
z **jednego** sprzężenia g_TGP.

**Test:**
- Skonstruować symbolicznie:
  ```
  L_TGP = (1/2) K(φ) ∂_μΦ ∂^μΦ - V(Φ)
        + g_TGP · Φ · O_BH (T_μν J^μ J^ν Wilson operator)
        + g_TGP · Φ · O_SC (Cooper-pair bilinear)
  ```
- Pod single-Φ axiom: tylko jedno Φ → tylko jedna stała g_TGP
- Pod Z₂: tylko parzyste potęgi Φ w L → bilinear (Φ²) jest najniższy nietrywialny
- Brak alternatywnych single-Φ Lagrangianów konsystentnych z Z₂ w bilinear

**Falsification:** istnieje single-Φ Lagrangian z **dwoma** różnymi
g_TGP coupling do BH i SC sektora bez łamania Z₂.

### T2.2 — BH-channel Wilson coefficient α₀

**Cel:** Wyprowadzić α₀ z T·J·J Wilson coefficient na photon-ring scale.

**Test:**
- Phase 1.B.3 / closure_2026-04-26 / F4: α₀ = 0.114/0.168² ≈ 4.045 (sympy
  rational 1069833/264500); strict Phase 2 BH.1 form: α₀ = (1/2)(1−3/3.88)/0.168² = 4.018
- Wilson coefficient: α₀ = g_TGP² · M_BH
- M_BH = (1/2)(1−3/3.88)/0.168² · (1/g_TGP²) — universal photon-ring matrix element
- Pod canonical normalization g_TGP = 1: α₀ = M_BH ≈ 4.018

**Falsification:** Wilson coefficient α₀ nie wyraża się jako
g_TGP² · M_BH z M_BH dependent only na (K(φ_eq), φ_eq, geometric).

### T2.3 — SC-channel Wilson coefficient κ_TGP

**Cel:** Wyprowadzić κ_TGP² z λ_sf normalization na Fermi surface.

**Test:**
- TGP-SC v2: λ_sf = κ_TGP² · g_J · μ_eff² (de Gennes scaling)
- κ_TGP² = pair-breaking Wilson coefficient: κ_TGP² = g_TGP² · M_SC
- M_SC = κ_TGP² · (1/g_TGP²) — universal substrate matrix element
- Pod canonical normalization g_TGP = 1: κ_TGP² = M_SC ≈ 4.048

**Falsification:** κ_TGP² nie wyraża się jako g_TGP² · M_SC z M_SC
dependent only na (K(φ_eq), φ_eq, geometric).

### T2.4 — Common-generator test: M_BH = M_SC

**Cel:** Sprawdzić, czy substrate matrix elements **konwergują** pod TGP
single-Φ + Z₂ + β=γ vacuum + Phase 1 covariant K(φ) = K_geo·φ⁴.

**Test:**
- Phase 1.A: K(φ) = K_geo · φ⁴ (Theorem alpha2)
- Phase 1.F.3: φ_eq = 1 (β=γ vacuum, V'(Φ_eq=1)=0 exact)
- F4 / closure_2026-04-26: K_geo = 1 (substrate-anchor calibration)
- Pod (K_geo=1, φ_eq=1): K(φ_eq) = 1
- M_BH and M_SC oba zależą od **tylko** K(φ_eq) i φ_eq (no other free parameters under TGP single-Φ)
- → M_BH = M_SC = M_universal ≈ 4.03 ± 0.04 (combined precision)

**Falsification:** if M_BH ≠ M_SC pod canonical normalization → identity
nie strukturalna (lokalne korekty O(1) sektor-zależne).

### T2.5 — RG flow audit: ratio α₀/κ_TGP² invariant pod RG

**Cel:** Sprawdzić, czy identity √α₀ = κ_TGP zachowuje się pod RG flow z
IR (skala observation) do UV (Planck scale).

**Test:**
- Pod TGP single-Φ axiom, oba operators (T·J·J, Cooper-pair) mają
  **wspólną anomalous dimension γ_an** (because they share same Φ-coupling)
- Ratio α₀/κ_TGP² jest RG-invariant: d(α₀/κ_TGP²)/dlnμ = 0 pod common γ_an
- Numerical check: pod one-loop: g_TGP(μ)² evolves uniformly w obu kanałach
- Identity preserved at all scales between observation and Planck

**Falsification:** if RG flow generates differential running α₀ ≠ κ_TGP²
at some scale → identity is IR-only coincidence, broken at UV.

### T2.6 — Phase 1 covariant 4D embedding

**Cel:** Sprawdzić, czy XS.1 identity konsystentnie embeduje się w
Phase 1 covariant 4D Lagrangian (existing 50-test closure).

**Test:**
- Phase 1.A.1 confirmed K(φ) = K_geo·φ⁴ pod (C1)–(C3) constraints
- Phase 1.B.2 anchored α₀ ≈ 4.04 pod target_shift = 0.114 (closure_2026-04-26)
- Phase 1.F.3 confirmed φ_eq = 1 (β=γ vacuum)
- F2/F3/F4 są LOCKED w PREDICTIONS_REGISTRY → wszystkie inputy do XS.1
  identity są strukturalnie pinned
- Identity √α₀ = κ_TGP follows from F2 ∧ F3 ∧ F4 + Phase 1.A.1 +
  Phase 1.B.2 + canonical normalization (no new postulate)

**Falsification:** if XS.1 identity introduces a constraint that
contradicts Phase 1 covariant 50-test closure → identity wymaga
modyfikacji axiomatic framework.

### T2.7 — Sympy-exact closure: status √α₀ = κ_TGP

**Cel:** Skompletować algebraic chain i przyporządkować klasyfikację:
DERIVED, PARTIALLY DERIVED, czy STRUCTURAL HINT.

**Test:**
- Algebraic chain: F1 (single-Φ) ∧ F2 (K_geo) ∧ F3 (β=γ) ∧ F4 (α₀ rational)
  + Phase 1.B.2 (target_shift) + canonical normalization
  → α₀ = M_universal = κ_TGP²
- All inputs are LOCKED w registry; chain has no free parameters
- Numerical match: 0.75% — within combined uncertainty (~1.4% combined ξ)
- Status:
  - **DERIVED**: if algebraic chain has zero ξ uncertainty
  - **PARTIALLY DERIVED**: if ξ unresolved factor remains (current Phase 2.BH.1.T2.4 status)
  - **STRUCTURAL HINT**: if only numerical match, no algebraic chain

**Falsification:** if algebraic chain has free parameter that is
sector-dependent (not universal) → identity remains STRUCTURAL HINT.

---

## Verdict gate

≥ 6/7 PASS → identity DERIVED or PARTIALLY DERIVED → **proceed to
Phase 3** (multi-sector falsification map).

≤ 5/7 PASS → identity remains STRUCTURAL HINT → Phase 3 still proceeds
but with weaker classification.

---

## Materiał wykonawczy

- **Skrypt:** [`phase2_substrate_action.py`](phase2_substrate_action.py) (7 sub-tests, sympy-symbolic + numerical)
- **Output:** `phase2_substrate_action.txt`
- **Memo:** `Phase2_results.md` (closure + classification + Phase 3 decision)

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 TGP/TGP_v1/research/op-cross-sector-charge/phase2_substrate_action.py 2>&1 | tee TGP/TGP_v1/research/op-cross-sector-charge/phase2_substrate_action.txt
```

---

## Constants used

```
alpha_0_rational = 0.114 / 0.168^2 = 4.04472 (F4 sympy rational 1069833/264500)
alpha_0_strict   = (1/2)(1 - 3/3.88) / 0.168^2 = 4.0179 (BH.1.Phase2)
kappa_TGP        = 2.0120 (TGP-SC v2; V/Nb/Ta/Mo/Pd RMS)
kappa_TGP_sq     = 4.0481

K_geo            = 1   (F2 LOCKED)
phi_eq           = 1   (F3 LOCKED, beta=gamma vacuum)
target_shift_BH  = (1/2)(1 - 3/3.88) = 0.1134 (Phase 1.B.3 + closure 2026-04-26)
target_shift_F4  = 0.114 (rational anchor sympy 1069833/264500)
psi_ph - 1       = 0.168 (M9.2-D universal photon ring)
g_TGP            = 1   (canonical normalization)
M_universal      = α₀ = κ_TGP² ≈ 4.03 (target match)

|alpha_0_rational - kappa_TGP^2| / kappa_TGP^2 = 0.082%  (F4 form)
|alpha_0_strict   - kappa_TGP^2| / kappa_TGP^2 = 0.747%  (Phase 2 strict)
```

---

## Cross-references

- [`Phase1_results.md`](Phase1_results.md) — feasibility audit 5/5 PASS
- [`program.md`](program.md) — XS.1 3-phase plan
- [`../op-bh-alpha-threshold/Phase2_results.md`](../op-bh-alpha-threshold/Phase2_results.md) — α₀ provenance + cross-sector hint origin
- [`../op-sc-alpha-origin/`](../op-sc-alpha-origin/) — κ_TGP provenance
- [`../op-phase1-covariant/`](../op-phase1-covariant/) — F2/F3/F4 LOCKED inputs
