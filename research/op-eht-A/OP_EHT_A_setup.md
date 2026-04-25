# OP-EHT-A — Track A: derive scenario (e) z M9.1'' first principles

**Data otwarcia:** 2026-04-25
**Status:** **OPEN** — Track A z OP-EHT closure (CONDITIONAL POSITIVE 13/18)
**Cel:** udowodnić że q-renormalization scenario (e) coupling
f(ψ) = √(g_tt^GR/g_tt^TGP) emerge naturalnie z M9.1'' Lagrangian,
rozwiązując Sgr A* tension (4σ → <1σ).

**Następnik:** OP-EHT (CLOSED 2026-04-25, 13/18 = 72%).
**Cross-references:** [[research/op-eht/OP_EHT_final_verdict.md]],
[[research/op-eht/OP_EHT_T3_results.md]] (scenario (e) ansatz),
[[paper/tgp_core.tex]] §matter coupling, F2 falsifiability.

---

## 1. Cel OP-EHT-A

OP-EHT closure pokazał że q-renormalization scenario (e) z fitted ansatz
f(ψ) = √(g_tt^GR/g_tt^TGP) **mathematically reduces** photon ring
deviation z +14.56% do +1.46% (within 5%, well within EHT envelope).
Track A pyta:

> **Czy taki coupling f(ψ) wynika naturalnie z TGP first principles
> (proper-time invariance, Z_2 symmetry, variational action)?**

Jeśli TAK → M9.1'' jest **fully reconciled** ze strong-field EHT, OP-EHT
closure upgrades z CONDITIONAL POSITIVE do **POSITIVE**, F4 falsification
window 2030+ zachowuje się jako wzmacnienie a nie test krytyczny.

Jeśli NIE → M9.2 conditional pivot wraca jako Track B (deferred).

## 2. Hipoteza centralna — proper-time matter coupling

W M9.1'' v2 axiom paper, matter coupling jest:
```
L_mat = -q × φ × ρ / Φ_0
```
gdzie φ = Φ - Φ_0 jest perturbation, q coupling constant, ρ rest density.

**Hipoteza Track A:** Naturalny coupling powinien być sourced w proper-time
frame matter, nie coordinate frame. W strong-field gdzie g_tt znacznie
różni się od -c_0², to daje dodatkowy faktor √(|g_tt|/c_0²):

```
L_mat^A = -q × φ × ρ × √(|g_tt(ψ)|/c_0²) / Φ_0
       = -q × φ × ρ × √((4-3ψ)/ψ) / Φ_0
```

Przy 1PN (ψ → 1), √((4-3ψ)/ψ) → 1 → standard matching γ_PPN = β_PPN = 1
zachowane.

W strong-field, effective source dla Φ-EOM jest renormalizowane przez
√(|g_tt|/c_0²), co modyfikuje effective A za który propaguje fotonowy
ring matching. Łatwy rachunek pokazuje że w photon ring (ψ=1.168):

```
A_eff = A × √(|g_tt^TGP|/|g_tt^GR|) ?
```

co odwraca scenario (e) — wymaga careful derivation.

## 3. Dlaczego proper-time coupling jest naturalny

1. **Lorentz invariance**: rest energy of point particle is m c² = m × √(-g_tt)
   in coord frame; identifies natural √|g_tt| weight in matter action.
2. **Z_2 symmetry**: ψ → ψ jest invariant; √|g_tt(ψ)| jest funkcją tylko ψ,
   więc Z_2 preserved.
3. **Variational consistency**: action principle naturally introduces
   √-g_full = √|g_tt| × √g_spatial w covariant volume measure.
4. **Compatibility with f·h = 1 substrate budget**: spatial g_ii follows
   from f_h(ψ) = ψ^something; product fh = 1 lock ⇒ √|g_tt| × g_spatial
   gives uniform substrate measure.

## 4. Plan testów OP-EHT-A

| # | Test | Cel | Metoda | Status |
|---|------|-----|--------|--------|
| **T-A1** | Derive proper-time coupling | Pokazać że L_mat^A wynika z natural action principle | sympy: variational derivation | **OPEN** |
| **T-A2** | 1PN matching | Zachowanie γ_PPN = β_PPN = 1 przy ψ → 1 | sympy: expand to 1PN, check coefficients | **OPEN** |
| **T-A3** | Strong-field photon ring | Numerical b_crit z renormalized coupling, sprawdzić deviation ≤ 5% | numpy: re-run photon ring | **OPEN** |
| **T-A4** | M9.1'' P3 audit re-validation | 3 PASS weak-field zachowane | re-run [[research/op-newton-momentum/M9_1_pp_P3_results.md]] | **OPEN** |
| **T-A5** | OP-7 c_GW = c_0 cross-check | GW propagation pod nowym coupling | analytic linearization → check σ_ab EOM | **OPEN** |
| **T-A6** | Synthesis + verdict | OP-EHT-A POSITIVE / NEGATIVE / INCONCLUSIVE | meta | **OPEN** |

## 5. Kryteria PASS

**OP-EHT-A POSITIVE** (M9.1'' rescue success):
- T-A1 PASS: derivation jest natural i unique (z principle)
- T-A2 PASS: 1PN γ_PPN = β_PPN = 1 zachowane
- T-A3 PASS: strong-field deviation ≤ 5% (within EHT envelope)
- T-A4 PASS: M9.1'' P3 audit (Mercury, Cassini, LLR) zachowany
- T-A5 PASS: c_GW = c_0 zachowane (OP-7 closure stays valid)
- ⇒ M9.1'' fully reconciled; OP-EHT upgrades CONDITIONAL → POSITIVE

**OP-EHT-A NEGATIVE** (M9.2 pivot mandatory):
- Jeśli któryś z T-A1..T-A5 FAIL → coupling nie jest naturalny w M9.1''
- ⇒ M9.2 conditional opening becomes mandatory, deferred do ngEHT 2030+

**OP-EHT-A INCONCLUSIVE**:
- T-A1..T-A3 PASS but T-A4 lub T-A5 FAIL → strong-field rescue ale weak-field
  cost; może wymagać partial M9.1''→M9.1''' refinement (intermediate pivot)

## 6. Pliki OP-EHT-A (planowane)

- `OP_EHT_A_setup.md` — ten plik
- `op_eht_a_T1_proper_time_derivation.py` — T-A1 (sympy: action variation)
- `op_eht_a_T1_proper_time_derivation.txt` — raw output T-A1
- `OP_EHT_A_T1_results.md` — werdykt T-A1
- `op_eht_a_T2_1pn_matching.py` — T-A2 (sympy: 1PN expansion)
- `op_eht_a_T2_1pn_matching.txt` — raw output T-A2
- `OP_EHT_A_T2_results.md` — werdykt T-A2
- `op_eht_a_T3_photon_ring.py` — T-A3 (numpy: re-run photon ring)
- `op_eht_a_T3_photon_ring.txt` — raw output T-A3
- `OP_EHT_A_T3_results.md` — werdykt T-A3
- `op_eht_a_T4_p3_revalidation.py` — T-A4 (Mercury/Cassini/LLR)
- `op_eht_a_T5_gw_propagation.py` — T-A5 (c_GW check)
- `OP_EHT_A_final_verdict.md` — synthesis

## 7. Bottom line — OP-EHT-A priority

OP-EHT-A jest **highest-priority follow-up** po OP-EHT closure:
- Uratowuje M9.1'' jeśli sukces
- Definiuje M9.2 pivot scope jeśli porażka
- Zamykany analitycznie + numerically w 1-2 dni
- **Single binary verdict**: M9.1'' alive lub M9.2 mandatory
