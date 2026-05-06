---
title: "Phase 0 balance sheet retrofit — op-eps-photon-ring (ε.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: ε.1 (op-eps-photon-ring)
cycle_path: "[[../op-eps-photon-ring/README.md]]"
auditor: Claudian
classification: STRUCTURAL
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - eps1
  - photon-ring
  - 137-anchor
  - high-risk-resolved
---

# Phase 0 balance sheet retrofit — ε.1 (op-eps-photon-ring)

## Metadata cyklu

- **Cykl:** [[../op-eps-photon-ring/README.md]]
- **Data oryginalnego closure:** 2026-04-29 (ε.1.Phase3 PASS 6/6)
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian
- **Klasyfikacja końcowa:** **STRUCTURAL** (downgrade z claim "PARTIALLY DERIVED")

## 1. Co cykl twierdzi że robi

Z [[../op-eps-photon-ring/Phase2_results.md]] werdykt:

> ε_ph = ψ_ph − 1 = 23/137 sympy-exact structural decomposition LOCKED;
> 5 alternative identity candidates falsified (drifts 10.7%, 14.4%, 52.9%,
> 181.8%, 127.5%); F4 chain implicit lock w 0.0019%.
> Classification: PARTIALLY DERIVED (refined).

Główne claims:

- **C1:** `ψ_ph = 160/137` photon-ring scale w M9.1'' Lorentzian metric
- **C2:** `ε_ph = ψ_ph - 1 = 23/137` (algebraic identity)
- **C3:** Prime-137 denominator "matches QED α_fine signature"
- **C4:** F4 chain implicit lock `ε_ph² = 529/18769 = (23/137)²`
- **C5:** EHT 2030+ (ngEHT) prediction `r_ph = (160/137)·r_g` z 0.1% precision

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs (PDG, CODATA, observational)

```
- CODATA α⁻¹(0) = 137.035999084(21)        [9 sig figs]
- M87* shadow size (EHT 2019)              [observational anchor for r_ph]
- ngEHT 2030+ projected precision 0.1%     [future falsification target]
```

### 2.2 Structural axioms (TGP-internal LOCKED)

```
- M9.1'' Lorentzian metric g_tt = -c²(4-3ψ)/ψ   [G.0 closure 2026-05-02]
- 4 = N_gen + 1                                  [from R3 ODE topology]
- 3 = N_gen                                      [why_n3 + G.0]
- ψ_horizon = 4/3                                [BH horizon w M9.1'']
- α₀ = 1069833/264500 ≈ 4.04474                  [F4 chain, UV.1 source]
```

**Critical anchor:**
- **137 prime** — NOT derived w ε.1; "matches α_fine signature" (claim z
  Phase1_results.md E1.2). Anchor borrowed from CODATA α_fine = 1/137.036.

### 2.3 Derived outputs (the cycle claims)

```
- Output 1: ψ_ph = 160/137                  (Phase1_results E1.2)
- Output 2: ε_ph = ψ_ph - 1 = 23/137         (algebraic, Phase1 E1.2)
- Output 3: r_ph = (160/137)·r_g             (ngEHT 2030+, Phase3 E3.1)
- Output 4: ε_ph² = 529/18769                (Phase2 E2.3)
```

### 2.4 Tautology test (CRITICAL)

**Sympy substitution dla ψ_ph = 160/137:**

```
ψ_ph = 4 / (3 + δ_ph)
     gdzie δ_ph = 17/40 = 0.4250  (claim z Phase1 E1.2)
     ψ_ph = 4 / (3 + 17/40) = 4 / (137/40) = 160/137
```

**Krytyczne pytanie:** skąd `δ_ph = 17/40 = 0.425`?

Zostało **odwrócone**: jeśli `ψ_ph = 160/137` jest "input" z M9.1'' photon-ring
calculation, to `δ_ph = 4/ψ_ph - 3 = 4·137/160 - 3 = 0.425` jest *re-expression*,
nie *derivation*.

**Realna kolejność:**
1. M9.1'' geodesic equation daje numerical photon-ring radius r_ph
2. r_ph normalized do r_g daje ψ_ph numerical (0.16788...)
3. **Sympy rationalization** dopasowuje `ψ_ph = 160/137` (sympy rational w
   simple form, denominator < 10⁵)
4. Prime decomposition: 160 = 2⁵·5, 137 = prime
5. **Observation:** 137 = α_fine⁻¹ denominator → "deep structural connection"

**Czy outputs kasują się tożsamościowo?**
- Częściowo: `ε_ph = ψ_ph - 1` to **trywialna algebra** (subtraction by 1)
- Output `23/137` ma niezależną informację od `137` jeśli 23 nie jest fitted
- Ale: 137 jest pre-fixed jako anchor, więc 23 = 160 - 137 = 23 jest **forced**
  by `ψ_ph = 160/137`. Trywialna konsekwencja.

**Real informational content:** ψ_ph = 160/137 (gdzie obie liczby pochodzą z
sympy rationalization 0.16788... numerical wartości z M9.1'' geometrii).

**Werdykt tautology test:** ⚠ **PASS QUALIFIED** — output `23/137` jest
trywialną konsekwencją ψ_ph; cała informacja w `160/137` rationalization.
Nie tautologia w sensie χ.1 (no canceling), ale "weak derivation" — sympy
rationalization to nie first-principles.

### 2.5 Falsifiability test (CRITICAL)

**Konkretny falsifier (Phase3 E3.1):**

```
Jeśli ngEHT 2030+ wykryje r_ph poza [0.99, 1.01]·(160/137)·r_g (0.1% band),
ε.1 FAILS.
```

**Band check:**
- TGP claim: ψ_ph = 160/137 = 1.16788
- ngEHT 2030+ projected band: 0.1%
- 5 alternative candidates falsified at drift > 5% (Phase2 E2.2)

**Critical observation:** Falsifier istnieje **operacyjnie** (ngEHT 2030+
przyniesie precision 0.1%, czyli 50× tighter niż obecne EHT). To dobry test.

**Ale:** 5 alternative candidates falsified są wszystkie **algebraic identities**
testowane ad-hoc (1.168/(2π), 23/160, 1/(2π·κ), etc.) — to jest *post-hoc
elimination* (autor próbował alternatives po fact). To **nie jest** alt-scan
≥4 first-principles candidates z falsifiability.

**Werdykt falsifiability test:** ⚠ **PASS QUALIFIED** — operacyjny falsifier
istnieje (ngEHT 2030+ 0.1%); ale alt-scan jest *algebraic* nie *physical*.

### 2.6 Independent-path cross-validation (CRITICAL for DERIVED)

**Path 1:** ψ_ph = 160/137 z M9.1'' geodesic equation (numerical → sympy rational)  
**Path 2:** F4 chain implicit lock `ε_ph² = 529/18769 = (23/137)²` (Phase2 E2.3)

**Convergence:** sympy-exact (529 = 23², 18769 = 137²)

**Krytyczne:** czy Path 2 jest *naprawdę niezależny*?

F4 chain mówi: `α₀ = target_shift / ε_ph²`, gdzie:
- `target_shift = 57/500 = 0.114` (heat-kernel a₂)
- `α₀ = 1069833/264500 ≈ 4.04474` (UV.1 source)

Substytucja:
```
ε_ph² = target_shift / α₀ = (57/500) / (1069833/264500)
      = 57·264500 / (500·1069833)
      = 15076500 / 534916500
      = 529/18769
      = (23/137)²    ✓ sympy-exact
```

**Czy F4 chain jest niezależny od M9.1''?**

UV.1 (cykl predecessor) używa `α₀ = 1069833/264500` jako anchor. Skąd to?
Sprawdzenie wymagałoby audit UV.1 (op-uv-as-ngfp). Jeśli α₀ jest **derived
z M9.1''**, to Path 2 nie jest niezależny od Path 1.

**Tymczasowo:** treat as 1 path (M9.1'' geometry); F4 chain consistency to
*algebraic re-arrangement*, nie nowa fizyczna ścieżka.

**Werdykt independent-path:** ⚠ **PARTIAL** — 1 niezależna ścieżka (M9.1''),
F4 chain to algebraic consistency check.

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists (this file)
⚠ Tautology test PASS (qualified — sympy rationalization to nie first-principles)
⚠ Falsifiability test PASS (qualified — alt-scan jest algebraic nie physical)
⚠ Independent-path cross-validation PARTIAL (1 path, F4 to algebraic check)
☐ Alt-scan ≥4 candidates with ≥3σ discrimination
   (5 alt-candidates są algebraic identities; nie ≥4 first-principles physical models)
☑ NIE used post-hoc structural motivations (sympy rationalization legit)
☐ NIE circular anchor — 137 anchor borrowed z α_fine, nie self-derived
☑ NIE inheriting drift > parent cycle drift × 5×
```

**3 ⚠ + 2 ☐ z 8 ⇒ status max STRUCTURAL.**

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? | Uzasadnienie |
|-------|----------|--------------|
| DERIVED FULL | NO | only 1 niezależna ścieżka, sympy rationalization, 137 anchor borrowed |
| DERIVED CONDITIONAL | NO | conditional na M9.1'' geometry only — nie 2 paths |
| **STRUCTURAL** | **YES** | algebraic structural constraint (ψ_ph = 160/137 z M9.1''); falsifier istnieje |
| ANSATZ | NO | jest stronger niż ansatz (concrete ngEHT prediction) |
| NUMEROLOGICAL | NO | nie pure coincidence — M9.1'' geometric derivation |
| TAUTOLOGY | NO | ε_ph = ψ_ph - 1 to weak derivation, ale ψ_ph ma niezależne info |

**Final verdict:** **STRUCTURAL** (downgrade z original "PARTIALLY DERIVED").

**Promotion path do DERIVED:** wymaga independent path 2 (np. ε_ph z innej
geometrii, NIE F4 chain re-arrangement). Lub: explicit derivation 137 anchor
z TGP topology (nie borrowed z α_fine).

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status | "PARTIALLY DERIVED (refined)" Phase2 | **STRUCTURAL** (downgrade) |
| Counter PREDICTIONS_REGISTRY | +18 | sub-tests mechanically OK; status downgrade w Phase 5 refactor |
| Sub-tests | 18/18 PASS | mechanically PASS, ale "PARTIALLY DERIVED" claim wymaga 2 paths |
| Independence | claim ψ_ph = 160/137 unique | weryfikacja: sympy rationalization, NIE first-principles unique |

## 6. Recommended action

- ☐ NO-OP — wymaga downgrade
- ☑ **DOWNGRADE w PREDICTIONS_REGISTRY**: status z `PARTIALLY DERIVED` na
  `STRUCTURAL` w Phase 5 registry refactor
- ☐ CRITIQUE — niepotrzebne (nie tautology / circular)
- ☑ **CASCADE_AUDIT**: η.2 zależy od ε.1 (137 anchor) — wymaga reconfirm
  η.2 conditional na ε.1 STRUCTURAL status
- ☐ CORE_IMPACT — brak

## 7. Notes

**Subtle issue:** ε.1 jest *algebraicznie* clean (sympy-exact 23/137,
529/18769) i ma operacyjny falsifier (ngEHT 2030+ 0.1%). Ale:
- "PARTIALLY DERIVED" claim wymaga **first-principles derivation** number
  137 z TGP topology — to jest **OPEN**
- 5 alt-candidates falsified są algebraic identities, nie physical models —
  to *nie jest* full alt-scan
- F4 chain "second path" jest re-arrangement, nie niezależna ścieżka

**Pozytywne:**
- M9.1'' geometric derivation jest legit (post-G.0 closure)
- 36.3× refinement gain coarse → refined wskazuje structural improvement
- ngEHT 2030+ falsifier jest concrete

**Conclusion:** ε.1 jest **STRUCTURAL** (algebraic constraint OK, ale claim
"DERIVED" przedwczesny bez first-principles 137 origin).

**Konsekwencja dla η.2:** η.2 prerequisite na 137 jest spełniony tylko
jeśli ε.1 jest accepted jako STRUCTURAL anchor. To znaczy η.2 sam jest
**STRUCTURAL+** (jego derivation cleaner niż ε.1, ale prerequisite chain).

## 8. Cross-references

- [[../op-eps-photon-ring/README.md]] — cykl audited
- [[../op-eps-photon-ring/Phase1_results.md]] — sympy rationalization
- [[../op-eps-photon-ring/Phase2_results.md]] — structural decomposition
- [[../op-eps-photon-ring/Phase3_results.md]] — ngEHT predictions
- [[../op-uv-as-ngfp/Phase3_results.md]] — predecessor (UV.1)
- [[retrofit_op-eta2_2026-05-06.md]] — η.2 retrofit (ε.1 prerequisite)
- [[README.md]] — M03 master plan
- [[audit_log.md]] — running log (added 2026-05-06)
- [[tracker.md]] — status updated do DONE_STRUCTURAL
- [[../../meta/CALIBRATION_PROTOCOL.md]] — protocol source
