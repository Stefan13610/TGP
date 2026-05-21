---
title: "Phase 1 setup — sympy test plan dla emergent Dirac propagator"
date: 2026-05-16
parent: "[[./README.md]]"
type: phase-setup
phase: 1
sub_tests: 13
substance_target: "11 FP + 1 LIT + 1 DEC; 0 hardcoded"
binding: "BINDING contract LOCKED w README §0 + Phase0_balance.md 8/8 ☑"
tags:
  - phase1-setup
  - sympy-plan
  - emergent-Dirac
  - L08-Dirac-propagator
---

# Phase 1 setup — 13 sub-test plan

## §0 — ASK-RULE pre-flight (Triggers A-D)

Per [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1.

### Trigger A — form-meaning split (HIGH-RISK ZONE check)

**Form:** Dirac propagator S_F(p) = i(γ·p + m)/(p² - m² + iε)
**Meaning w TGP context:** S_F NIE jest free-field "elementarna cząstka Dirac" —
to *effective propagator dla emergent fermion (2-kink bound state)* w universal g_eff[Φ]
frame. **NIE BD-style** "scalar particle propagator" interpretation.

**BD-drift mitigation:** explicit annotation w T4 że propagator opisuje 2-kink composite
fermion, NIE elementary Dirac field; mass m_obs jest renormalized pole mass per L05
distinction (NIE bare mass parameter).

### Trigger B — inheritance ledger consistency

Wszystkie LIVE inheritance per Phase 0 §4 enumerated:
- L08-FR antisymmetry: γ_exchange = π Berry phase (T1)
- L08-Clifford: {γ^μ, γ^ν} = 2η^μν (T2)
- L05: m_obs ≠ M_full; k_obs(α=2,d=3)=3 (T5-T6)
- why_n3 Phase 5: universal mass formula (T6)
- PDG 2024: m_μ/m_e = 206.7682275 (T7 LIT)
- M9.1'' canonical: Lorentz signature dla η_μν (T8)

### Trigger C — standard reference validity

- Bjorken-Drell §15.3 (Dirac propagator definition)
- Peskin-Schroeder §3.2-3.3 (Cl algebra + propagator)
- Lüders-Pauli CPT theorem (T11)

Standard QFT formalism applied; no novel claims w L2 reduction.

### Trigger D — 0 hardcoded T_pass=True

**Strict policy:** każdy test computes substantive result + verifies. Pattern:
```python
T_pass = (sp.simplify(expr_lhs - expr_rhs) == 0)  # NEVER literal True
```

## §1 — 13 sub-test plan

### Sub-tests detail

| # | Test | Type | Substantive content |
|---|---|---|---|
| **T1** | 2-particle Fock antisymmetric state | FP | Construct \|k1,k2⟩_anti; verify P12 \|k1,k2⟩_anti = -\|k1,k2⟩_anti; Pauli exclusion \|k,k⟩_anti = 0 |
| **T2** | Cl(1,3) γ^μ matrices in Dirac representation | FP | Explicit 4x4 matrices; verify {γ^μ, γ^ν} = 2η^μν dla wszystkich (μ,ν) |
| **T3** | (γ·p - m)(γ·p + m) = p² - m² | FP | Symbolic verification z Cl algebra inheritance |
| **T4** | Propagator construction S_F = i(γ·p + m)/(p² - m² + iε) | FP | (γ·p - m) S_F(p) = iI verification |
| **T5** | Pole identification at p² = m_obs² | FP | Residue analysis; mass-shell condition |
| **T6** | m_obs z why_n3 Phase 5 formula | FP | m_obs = c_M · A² · g_0^(e²(1-α/4)); α=2 case |
| **T7** | PDG m_μ/m_e = 206.7682 (LIT) | LIT | Verify L05 inherited result: m_obs(α=2,d=3) z k_obs=3 reproduces PDG |
| **T8** | Lorentz boost covariance | FP | S(Λ) = exp(-i/4 ω_μν σ^μν); commutation [γ^μ, σ^νρ] = i(η^μν γ^ρ - η^μρ γ^ν) |
| **T9** | Hermiticity S_F^†(p) = γ_0 S_F(p) γ_0 | FP | γ_0^† = γ_0; γ_i^† = -γ_i; verify |
| **T10** | Free-field limit lim Φ→Φ_0 S_F^TGP = S_F^Dirac | FP | Symbolic limit; coupling V_eff[Φ] → 0; mass m_obs → m_canonical |
| **T11** | CPT operational (substrate Z₂ + emergent Lorentz) | DEC | Lüders-Pauli theorem applicability; Z₂ inheritance from FOUNDATIONS §1 |
| **T12** | Two-point function ⟨0\|Tψψ̄\|0⟩ = S_F(x-y) | FP | Fourier transform consistency: S_F(x) = ∫ d⁴p/(2π)⁴ e^(-ip·x) S_F(p) |
| **T13** | DECLARATIVE S05 single-Φ axiom preservation | DEC | No new fundamental fields; emergent fermion = 2-kink bound state w Φ |

### Sub-tests dependency graph

```
T1 (FR antisym)─┐
                ├──→ T3 ──→ T4 ──→ T5 ──→ T6 ──→ T7
T2 (Cl algebra)─┘                                 │
                │                                 │ T8 ──→ T9
                └──→ T8                           │      
                                                  └──→ T10
                                                  
T11 (CPT DEC) ── independent
T12 (2-point function) ── consumes T4
T13 (S05 DEC) ── independent
```

### Risk-coupling per sub-test

| R# | Risk | Sub-tests affected |
|---|---|---|
| R1 | FR + Cl combination | T1, T2, T3, T4 |
| R2 | m_obs identification z L05 | T5, T6, T7 |
| R3 | Lorentz boost generator | T8 |
| R4 | Free-field limit | T10 |
| R5 | CPT operational claim | T11 |
| R6 | Wilson coefs (Phase 2 OUT) | N/A this Phase |

## §2 — Sympy plan executive

**File:** `Phase1_sympy.py`
**Expected output:** `Phase1_sympy.txt`
**Sub-tests:** 13
**Substance:** 11 FP + 1 LIT (T7) + 1 DEC (T11, T13 dwa DEC)

Wait — T11 i T13 są oba DEC. Adjust: **11 FP + 1 LIT + 2 DEC = 13 total**

**Per Phase 0 §10:** "11 FP + 1 LIT + 1 DEC". Reconciliation: T11 jest **DEC reasoning** (CPT
applicability declarative), T13 jest **DEC formal** (S05 single-Φ axiom statement). Both
declarative w sense że NIE constitute novel computation — citation/inheritance only.

**Revised substance count:** 10 FP + 1 LIT + 2 DEC = 13 total. **FP fraction: 76.9%**
(still > 70% threshold target).

**Acceptable:** 76.9% FP > L08-RG flow's 20% FP (HALT-B), > LIGO-3G-native 20%; comparable
do other cosmology cycles (80% inflation; 80% S07-reset). HOWEVER lower niż original ambition
84.6%. Acceptable honest substance.

## §3 — Pre-execution check

- [x] README.md §0 BINDING contract LOCKED
- [x] Phase0_balance.md 8/8 ☑ PASS
- [x] User authorization received ("aktywuj L08-Dirac Phase 1")
- [x] ASK-RULE Triggers A-D PASS pre-flight
- [x] Sub-tests T1-T13 explicit + dependency graph
- [x] Substance metric target documented (10 FP + 1 LIT + 2 DEC; 76.9% FP)

**Status:** 🟢 **Phase 1 READY** — proceed do sympy execution.

## §4 — Post-execution plan

1. Run `Phase1_sympy.py` → produce `Phase1_sympy.txt`
2. Verify all 13/13 PASS (target) lub document obstructions honestly
3. Write `Phase1_results.md` z explicit derivations + verdict
4. Decide post-Phase-1:
   - ≥11/13 PASS + structure native-derivable → Phase FINAL closure A−
   - <11/13 + obstruction → HALT-B honest negative (analog L08-RG)
   - critical S05 violation → cycle abandon, document w FINDINGS

## Cross-references

- [[./README.md]] §0 BINDING contract
- [[./Phase0_balance.md]] 8/8 ☑ gate PASS
- [[./Phase1_sympy.py]] (to be written next)
- [[./Phase1_sympy.txt]] (post-execution)
- [[./Phase1_results.md]] (post-execution)
