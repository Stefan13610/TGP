---
title: "Phase 0 balance — op-MAG-Phase5-V-reference-clarification-2026-05-09"
date: 2026-05-09
type: phase-balance-sheet
status: COMPLETE
parent: "[[./README.md]]"
phase: 0
---

# Phase 0 balance — Phase 5 V_orig usage clarification

## Confirmed facts

### Fact 1: Phase 5 explicit cytuje V_orig (DEPRECATED)

**Phase 5 results (markdown), line 38:**
```
V(Φ) = -β Φ³/(3 Φ_0) + γ Φ⁴/(4 Φ_0²)
```

**Phase 5 sympy, line 63:**
```python
V_Phi = -beta_p * Phi**3 / (3 * Phi_0) + gamma_p * Phi**4 / (4 * Phi_0**2)
```

To **literalna kopia** V_orig (sek08a linie 95-110), DEPRECATED 2026-05-02.

### Fact 2: Phase 5 closed AFTER V_M9.1'' canonical lock

| Cykl | Closure date |
|---|---|
| G.0 closure (V_M9.1'' canonical lock) | **2026-05-02** |
| Phase 5 MAG (V_orig usage) | **2026-05-09** |

Phase 5 closed **7 dni PO** V_M9.1'' canonical lock, ale używa deprecated formuły.
To **NIE** "pre-G.0 historical artifact" (jak T-Λ closure 2026-04-26).

### Fact 3: Phase 5 derivation jest STRUCTURE-DEPENDENT na V_orig

Phase 5 derivation steps:
1. Expansion around vacuum Φ = Φ_0 (V_orig vacuum z β=γ)
2. V''''(Φ_0) = 6γ/Φ_0² (V_orig 4th derivative)
3. λ_4 = V''''/4 = 3γ/(2 Φ_0²) (positive)
4. m_Mach = λ_4 · ⟨δΦ²_bg⟩ · ∫δΦ²_sol

**Każdy krok zależy od V_orig form**. Re-derivation z V_M9.1'' wymaga
recompute wszystkich derivatives przy V_M9.1'' minimum (ψ=2/3 lub ψ=1).

## Key question dla Phase 1

**Czy V''''/4 (λ_4) zachowuje znak i magnitude pod V_orig→V_M9.1'' canonical?**

Pre-cycle expectation:
- V_orig λ_4 = +3γ/(2 Phi_0²) (positive, dimensional)
- V_M9.1'' λ_4 = ? (sign + magnitude TBD by sympy)

Jeśli λ_4 zmienia znak (positive→negative), m_Mach zmienia znak — **physical
implication poważna**.

## Three resolution paths (z README)

| Path | Description | Workload |
|---|---|---|
| A (lightweight) | Re-interpretation Phi_0_Phase5 = (2/3)·Phi_0_V_M911 | 1 cykl, sympy verify |
| B (medium) | Re-derivation Phase 5 around V_M9.1'' true minimum | 1-2 cykle, full re-do |
| C (heavy) | Multi-vacuum identification (Phase 5 around different ψ) | osobny cykl |

## 8/8 gate criteria

| # | Criterion | Status |
|---|---|---|
| 1 | Audit cycle, lightweight scope | ✅ |
| 2 | Phase 5 V_orig usage explicit verified (Fact 1) | ✅ |
| 3 | Pre-cycle hipoteza H1 (λ_4 sign change) | ✅ |
| 4 | NIE multi-candidate fit / drift selection | ✅ N/A |
| 5 | NIE constructed criterion | ✅ |
| 6 | Sympy verification dla Phase 1 sign comparison | 🟡 in progress |
| 7 | Honest reporting jeśli Path B/C wymagane | 🟡 awaits sympy |
| 8 | Cross-reference do parent audit + Phase 5 | ✅ |

## External anchors

- m_e = 0.511 MeV (PDG)
- v_EW = 246 GeV (Higgs VEV scale, Phase 5 best scenario)
- m_C ~ 1.5×10⁻³³ eV (Phase 5 Compton/cosmological scale)

## TGP-internal axioms

- V_M9.1''(ψ) = -γψ²(4-3ψ)²/12 (canonical)
- V_M9.1'' minimum: ψ=2/3, V=-4γ/27
- V_orig (DEPRECATED): -βΦ³/(3Φ_0) + γΦ⁴/(4Φ_0²)

## Cross-references

- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] linia 38
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]] linia 63
- [[../op-V-canonical-consistency-audit-2026-05-09/Phase1_audit_results.md]] §1.2

## Status

**Phase 0 COMPLETE.** Ready for Phase 1 sympy.
