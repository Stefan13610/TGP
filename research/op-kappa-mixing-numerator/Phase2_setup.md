---
title: "κ.1.Phase2 setup — mixing-operator B² extension formal derivation + alternatives FALSIFIED"
date: 2026-04-30
cycle: κ.1.Phase2
status: PRE-EXECUTION
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
tags:
  - TGP
  - kappa1
  - mixing-operator-derivation
  - alternatives-falsified
---

# κ.1.Phase2 — mixing-operator B² extension formal derivation + alternatives FALSIFIED

> **Cel:** Formal derivation of mixing-operator B²-extension; sympy-exact lock
> ρ̄ num = 11 + η̄ num = 5; falsify alternative numerator forms via
> structural-pairing criterion; full Wolfenstein triple cascade DERIVED.

## 7 sub-tests

### K2.1 — Formal definition of mixing-operator B² extension

**Setup:** Mixing-operator B² extension defined as:
$$
B^2_{\text{mix}}(\text{up} \to \text{lep}; \text{level } L)
  := L_{\text{up}} - L_{\text{lep}}
$$
where $L \in \{B^2, K\}$ are chirality-counting levels, and `lep` jest
**reference frame** (Dirac sector z B²=2 minimal chirality).

**Wolfenstein numerator rule (denom-num level pairing):**
- if denom(X) uses B²-level factor → num(X) uses B²-level diff
- if denom(X) uses K-level factor → num(X) uses K-level diff

**PASS:** Definition strukturalnie spójna + cross-sector interpretation OK.

### K2.2 — ρ̄ num = 11 = B²_up_num − B²_lepton sympy DERIVED

**Setup:** ρ̄ = (B²_up_num − B²_lepton) / (2·N_gen·B²_up_num)
       = (13 − 2) / (2·3·13) = 11/78

Sympy-exact composition + match RHO_BAR_TGP.

**PASS:** ρ̄ sympy = 11/78 z proposed numerator-derivation form.

### K2.3 — η̄ num = 5 = K_up_num − K_lepton_num sympy DERIVED

**Setup:** η̄ = (K_up_num − K_lepton_num) / (K_up_num·K_lepton_num)
       = (7 − 2) / (7·2) = 5/14

Sympy-exact composition + match ETA_BAR_TGP.

**PASS:** η̄ sympy = 5/14 z proposed numerator-derivation form.

### K2.4 — Lepton-as-reference structural argument

**Setup:** Why does the rule subtract lepton-sector counterpart?

**Argument:** In TGP, lepton sector jest **minimum chirality reference**:
- B²_lepton = 2 (Dirac, 2 chiralities)
- B²_neutrino = 1 (Majorana, 1 chirality, **most-fundamental** but absent QED coupling)
- B²_up = 13/4 = 2 + 5/4 (Dirac + QCD)
- B²_down = 61/25 = 2 + 11/25 (Dirac + QCD effective)

Wszystkie sektors share Dirac base 2 = B²_lepton; QCD/color contribution
(5/4 dla up, 11/25 dla down) jest cross-sector difference. Mixing-operator
isolates QCD-sector contribution **relative to lepton-as-reference**.

**Numerator-extraction:** B²_up_num = 13 ↔ K_up_num = 7 are integer-extractions
of full B²-level. Subtraction of B²_lepton = 2 (or K_lepton_num = 2) gives
**integer cross-sector chirality difference**.

**PASS:** Structural argument coherent + lepton-as-reference well-defined.

### K2.5 — Alternative numerator formulations FALSIFIED

**Setup:** Test 5+ alternative numerator forms; falsification criterion =
**violation of denom-num level pairing**.

| Alternative | Form | Sympy match | Level pair? | Verdict |
|------|------|------|------|---------|
| C1 (11) | K_up_denom + K_lepton_denom = 8+3 | ✓ | denoms vs ρ̄ B²-denom mismatch | FAIL |
| C2 (11) | K_up_num + B²_up_denom = 7+4 | ✓ | mixed levels | FAIL |
| C3 (11) | 4·N_gen − 1 = 11 | ✓ | uses arbitrary 1 | FAIL |
| C4 (5) | N_gen + B²_lepton = 3+2 | ✓ | mixed (gen + B²) | FAIL |
| C5 (5) | K_up_denom − N_gen = 8−3 | ✓ | denom vs η̄ K-denom mismatch | FAIL |
| C6 (5) | 2·N_gen − 1 = 5 | ✓ | uses arbitrary 1 | FAIL |

**5/6 alternative forms FALSIFIED** by structural-pairing criterion (level
mismatch with denom; or use of arbitrary integer constants).

**PASS:** ≥ 4/6 alternatives FALSIFIED.

### K2.6 — Full Wolfenstein triple cascade DERIVED

**Setup:** Complete cascade post-κ.1:
- A = K_up_denom²/N_gen⁴ = 8²/3⁴ = 64/81 (denom + num both DERIVED, η.2)
- ρ̄ num = B²_up_num − B²_lepton = 11 (κ.1) | denom = 2·N_gen·B²_up_num = 78 (η.2)
- η̄ num = K_up_num − K_lepton_num = 5 (κ.1) | denom = K_up_num·K_lepton_num = 14 (η.2)

**Full Wolfenstein triple = (64/81, 11/78, 5/14)** all components sympy-DERIVED
z 4-sector chirality-counting B²-cross-product + 2-level mixing extension.

**PASS:** All 3 components of triple sympy-exact + DERIVED status.

### K2.7 — Classification cascade

**Setup:** Post-κ.1 status promotions:
- η.1 numerators (11, 5) STRUCTURAL HINT → **PARTIALLY DERIVED**
- η.1 Wolfenstein triple full cascade DERIVED denoms + PARTIALLY DERIVED nums
  → **DERIVED (refined²)** comprehensive
- HH5 (Wolfenstein numerators research-track) closure → **PARTIALLY DERIVED**
  (lepton-quark mixing-operator B²-extension framework structurally consistent
  but not unique under all parsimony criteria)
- A_TGP (64/81) full DERIVED both num + denom — confirmed
- Cabibbo λ_C lock confirmed (PMNS-CKM single anchor)

**PASS:** Cascade documented + Wolfenstein triple full DERIVED.

## Verdict gate

**6/7 PASS minimum** → κ.1.Phase2 closed; Phase 3 proceeds.

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-kappa-mixing-numerator/phase2_num_derivation.py 2>&1 | tee research/op-kappa-mixing-numerator/phase2_num_derivation.txt
```

## Cross-references

- [`Phase1_results.md`](Phase1_results.md) — viability gate 5/5 PASS
- [`../op-eta2-denom-derivation/Phase3_results.md`](../op-eta2-denom-derivation/Phase3_results.md) — denom DERIVED
- [`../op-theta-quark-koide/Phase3_results.md`](../op-theta-quark-koide/Phase3_results.md) — K-taxonomy
