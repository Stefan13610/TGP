---
title: "κ.1 program — Wolfenstein numerator derivation via lepton-quark mixing-operator B²-extension"
date: 2026-04-30
cycle: κ.1
status: PRE-EXECUTION
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-eta2-denom-derivation/Phase3_results.md]]"
tags:
  - TGP
  - kappa1
  - wolfenstein
  - numerator-derivation
  - mixing-operator
  - chirality-counting
---

# κ.1 program — Wolfenstein numerator derivation via lepton-quark mixing-operator B²-extension

> **Cel:** Lift the (11, 5) STRUCTURAL HINT for ρ̄/η̄ numerators
> (post-η.2 denom-DERIVED, num-STRUCTURAL HINT) do PARTIALLY DERIVED via
> rigorous derivation z lepton-quark mixing-operator B²-extension.

## Hypothesis κ.1

**Unified rule:** For each Wolfenstein cross-sector parameter X = a/b, where
the denominator b uses an up-quark sector factor F_up, the numerator a is the
**difference** between F_up and its lepton-sector counterpart F_lep.

**ρ̄ numerator (B²-level difference):**
- denom(ρ̄) = 78 = 2·N_gen·**B²_up_num** = 2·3·13
- num(ρ̄) = 11 = **B²_up_num − B²_lepton** = 13 − 2

**η̄ numerator (K-level difference):**
- denom(η̄) = 14 = **K_up_num**·K_lepton_num = 7·2
- num(η̄) = 5 = **K_up_num − K_lepton_num** = 7 − 2

**Structural interpretation:** the apex (ρ̄, η̄) of the unitarity triangle
encodes a 2-level chirality-counting cross-sector difference between up-quark
and lepton sectors:
- Real axis (ρ̄) → B²-level mixing (B²_up vs B²_lepton)
- Imaginary axis (η̄) → K-level mixing (K_up vs K_lepton)

**Identity check** (since K_lepton_num = 2 = B²_lepton):
- ρ̄ num: 11 = B²_up_num − B²_lepton (sympy exact)
- η̄ num: 5 = K_up_num − B²_lepton = K_up_num − K_lepton_num (sympy exact, dual form)

## 3-phase plan

**Phase 1 (5 sub-tests):** numerator landscape + alternatives ranking
- K1.1: Inventory of all sympy-LOCKED numerators across closed cycles
- K1.2: Sector-unique primes inventory (11 ρ̄, 5 η̄, 23 ε_ph, 37 K_down, 13 B²_up)
- K1.3: Mixing-operator hypothesis sympy test (11 = B²_up_num − B²_lepton; 5 = K_up_num − K_lepton_num)
- K1.4: Alternative numerator forms ranked z TGP-natural primitives
- K1.5: Viability gate dla mixing-operator B²-extension framework

**Phase 2 (7 sub-tests):** rigorous derivation + alternatives FALSIFIED
- K2.1: Formal definition of mixing-operator B² extension
- K2.2: ρ̄ num = 11 = B²_up_num − B²_lepton sympy DERIVED z full Wolfenstein composition
- K2.3: η̄ num = 5 = K_up_num − K_lepton_num sympy DERIVED z full Wolfenstein composition
- K2.4: Why subtraction of lepton sector? Lepton-as-reference structural argument
- K2.5: 5+ alternative numerator formulations FALSIFIED via parsimony criterion
- K2.6: Cross-check: full Wolfenstein triple (A, ρ̄, η̄) = (64/81, 11/78, 5/14) all components DERIVED + FULLY-DERIVED status
- K2.7: Classification cascade — η.1 numerators STRUCTURAL HINT → PARTIALLY DERIVED; Wolfenstein full cascade DERIVED

**Phase 3 (6 sub-tests):** predictions KK1-KK6
- K3.1 (KK1): Belle II 2027+ |V_ub|_κ.1 sharper window post-numerator DERIVATION
- K3.2 (KK2): Unitarity triangle β post-κ.1 hardened sin(2β) = 0.7090
- K3.3 (KK3): Cross-sector mixing-operator B² extension uniqueness test
- K3.4 (KK4): Wolfenstein full closure — NO new free parameters in CKM matrix
- K3.5 (KK5): Future ι.1 cycle hint — extend mixing-operator do (ν, up), (lep, down), (lep, ν) sektor pairs
- K3.6 (KK6): 5-channel κ.1 falsification convergence

## Verdict gate

**Phase 1:** ≥ 4/5 PASS → Phase 2 viable
**Phase 2:** ≥ 6/7 PASS → Phase 3 proceeds; classification PARTIALLY DERIVED
**Phase 3:** ≥ 5/6 PASS → κ.1 program END

**Total target:** 18/18 sub-tests, ledger 481 → 499.

## Falsification criteria

- If alternative numerator forms (e.g. 11 = K_up_denom + K_lepton_denom = 8+3,
  5 = N_gen + B²_lepton = 3+2) admit **fewer or equal** TGP-natural framework
  constants in their construction, classification stays STRUCTURAL HINT.
- If proposed forms break under cross-sector parsimony test, demote → STRUCTURAL HINT.

## Cross-references

- [[../op-eta-wolfenstein/Phase3_results.md]] — η.1 Wolfenstein triple LOCKED z (64/81, 11/78, 5/14)
- [[../op-eta2-denom-derivation/Phase3_results.md]] — η.2 denoms DERIVED + ρ̄/η̄ num STRUCTURAL HINT
- [[../op-theta-quark-koide/Phase3_results.md]] — K-taxonomy 4-sector: K_up=7/8, K_lepton=2/3, K_ν=1/2, K_down=37/50
- [[../../INDEX.md]] — master ledger
- [[../../PREDICTIONS_REGISTRY.md]] — HH5 entry promotion target

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-kappa-mixing-numerator/phase1_num_landscape.py 2>&1 | tee research/op-kappa-mixing-numerator/phase1_num_landscape.txt
PYTHONIOENCODING=utf-8 python -X utf8 research/op-kappa-mixing-numerator/phase2_num_derivation.py 2>&1 | tee research/op-kappa-mixing-numerator/phase2_num_derivation.txt
PYTHONIOENCODING=utf-8 python -X utf8 research/op-kappa-mixing-numerator/phase3_kappa_predictions.py 2>&1 | tee research/op-kappa-mixing-numerator/phase3_kappa_predictions.txt
```
