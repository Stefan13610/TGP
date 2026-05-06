---
title: "Phase 0 balance sheet retrofit — op-iota-charge-pmns-unification (ι.1)"
date: 2026-05-06
parent: "[[README.md]]"
type: balance-sheet-retrofit
cycle_audited: ι.1 (op-iota-charge-pmns-unification)
cycle_path: "[[../op-iota-charge-pmns-unification/README.md]]"
auditor: Claudian
classification: ANSATZ
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - phase0
  - balance-sheet-retrofit
  - iota1
  - PMNS
  - mixing-operator
  - kappa-contagion
  - critical-downgrade
---

# Phase 0 balance sheet retrofit — ι.1 (op-iota-charge-pmns-unification)

## Metadata cyklu

- **Cykl:** [[../op-iota-charge-pmns-unification/README.md]]
- **Data oryginalnego closure:** 2026-04-30 (ι.1.Phase2 PASS 7/7, FULL CASCADE)
- **Data retrofit:** 2026-05-06
- **Auditor:** Claudian
- **Klasyfikacja końcowa:** **ANSATZ** (severe downgrade z claim "DERIVED FULL CASCADE")

## 1. Co cykl twierdzi że robi

Z [[../op-iota-charge-pmns-unification/Phase2_results.md]] werdykt:

> Verdict: 7/7 PASS — FULL CASCADE. PMNS angles DERIVED via mixing-operator
> framework (analog do κ.1 CKM closure); charge-sector unification PARTIALLY
> DERIVED via 3-way cascade lock; 5+ alternative PMNS forms FALSIFIED.

Główne claims:

- **C1:** sin²θ₁₂ = 1/N_gen = 1/3 (trimaximal)
- **C2:** sin²θ₁₃ = K_ν·λ_C² = (1/2)·λ_C² ≈ 0.0254
- **C3:** sin²θ₂₃ = K_ν = 1/2 (maximal)
- **C4:** 3-way charge cascade lock (|Q_u|²-|Q_d|² = 1/N_gen ↔ B²-diff ↔ α-residual)
- **C5:** Mixing-operator framework analog do κ.1
- **C6:** "ζ.1 PMNS angles PARTIALLY DERIVED → DERIVED (mixing-operator)"

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs (NuFit 5.3, PDG)

```
- NuFit 5.3 sin²θ₁₂ = 0.307 ± 0.013   [3 sig figs, 4.2% band]
- NuFit 5.3 sin²θ₁₃ = 0.022 ± 0.0007  [2 sig figs, 3.2% band]
- NuFit 5.3 sin²θ₂₃ = 0.572 ± 0.024   [3 sig figs, 4.2% band]
- PDG δ_CP = -1.97 ± 0.6 rad           [pending, JUNO/DUNE 2030+]
```

### 2.2 Structural axioms (TGP-internal)

```
- N_gen = 3                                    [why_n3 + R3 ODE]
- K_ν = 1/2 (Majorana, 1 chirality)            [θ.1 K-taxonomy + ζ.1]
- B²_up = 13/4, B²_down = 61/25                [θ.1 — NUMEROLOGICAL z M03 retrofit!]
- λ_C = 451/2000 = 0.2255                      [ζ.1; STATUS: NIE-AUDYTOWANE]
- Mixing-operator framework                    [κ.1 — NUMEROLOGICAL z M03 retrofit!]
```

**Critical:** ι.1 zależy od **κ.1 (NUMEROLOGICAL)** + **θ.1 K_down (NUMEROLOGICAL)**
+ ζ.1 (nie-audytowane). Cascade contamination heavy.

### 2.3 Derived outputs

```
- sin²θ₁₂ = 1/3 = 0.333         (NuFit 0.307, drift 8.58%)
- sin²θ₁₃ = 0.0254              (NuFit 0.022, drift 15.57%)
- sin²θ₂₃ = 1/2 = 0.500          (NuFit 0.572, drift 12.59%)
- 3-way cascade identities      (algebraic, all from B²_up - B²_down = 81/100)
```

### 2.4 Tautology test (CRITICAL)

#### Output 1: sin²θ₁₂ = 1/3

```
sin²θ₁₂ = 1/N_gen = 1/3 (trimaximal pattern, S₃ democratic permutation)
```

Z TGP axioms (N_gen=3) — algebraicznie OK. Pattern strukturalny.

**Werdykt tautology K_θ₁₂:** ✅ **PASS** (algebraic identity z N_gen).

#### Output 2: sin²θ₁₃ = (1/2)·λ_C²

```
sin²θ₁₃ = K_ν · λ_C² = (1/2)·(0.2255)² = 0.02543
NuFit measurement: 0.022 ± 0.0007
Drift: 15.57% — wykracza znacznie poza 3.2% PDG band
```

**Krytyczne:**
- Output zależy od λ_C (ζ.1, NIE-audytowane)
- TGP claim 0.0254 jest **outside** NuFit 1σ band (0.0213-0.0227)
- Drift 15.57% / band 3.2% = ~5σ tension
- Author używa "zeroth-order gate <25%" — to **non-standard accommodating gate**

**Werdykt tautology K_θ₁₃:** ❌ **FAIL** — claim "DERIVED" przy 5σ tension
vs eksperyment to **NIE jest derivation**.

#### Output 3: sin²θ₂₃ = 1/2

```
sin²θ₂₃ = K_ν = 1/2 (maximal mixing)
NuFit measurement: 0.572 ± 0.024
Drift: 12.59% — wykracza znacznie poza 4.2% band
```

**Krytyczne:**
- TGP claim 0.5 jest na **lower edge of NuFit 95% CL** (typically [0.41, 0.62])
- Drift 12.59% / 1σ band 4.2% = ~3σ tension
- Maximal mixing była *historic prediction* przed precision measurements;
  NuFit 5.3 + T2K + NOvA disfavor it

**Werdykt tautology K_θ₂₃:** ❌ **FAIL** — claim "DERIVED" przy 3σ tension.

#### Output 4: 3-way charge cascade

```
1. |Q_u|² - |Q_d|² = (2/3)² - (1/3)² = 4/9 - 1/9 = 3/9 = 1/3 = 1/N_gen ✓
2. B²_up - B²_down = 13/4 - 61/25 = 81/100 ✓
3. 2·(81/100)/(9·5) = 162/4500 = 9/250 ✓
```

Identity 1 jest **trywialna algebra** dla SM charges + N_gen=3.
Identity 2 zależy od B²_up=13/4 (forced by θ.1 K_up STRUCTURAL) +
**B²_down = 61/25 (forced by θ.1 K_down NUMEROLOGICAL)**.
Identity 3 jest η.2 derivation Form B (Form A independent).

**Cascade contamination:** identity 2 inherits θ.1 K_down NUMEROLOGICAL status.

**Werdykt tautology cascade:** ⚠ **PARTIAL** — identities sympy-exact, ale
inputs contaminated.

### 2.5 Falsifiability test (CRITICAL)

**Wszystkie 3 PMNS angles już *outside* NuFit 1σ band (drifts 8-15%):**

| Angle | TGP claim | NuFit 5.3 (1σ) | Drift / σ |
|-------|-----------|-----------------|-----------|
| sin²θ₁₂ | 0.333 | [0.294, 0.320] | **+2.0σ** outside upper edge |
| sin²θ₁₃ | 0.0254 | [0.0213, 0.0227] | **+5.4σ** outside upper edge |
| sin²θ₂₃ | 0.500 | [0.548, 0.596] | **-3.0σ** outside lower edge |

**Author's "zeroth-order gate <25%"** mówi: "drift <25% to OK". To jest:
- **Non-standard physics gate** (typowo 1σ-3σ to threshold)
- **Accommodating** — pozwala on TGP claim "DERIVED" mimo 3-5σ tensions
- **Non-falsifiable** — żaden NuFit precision improvement nie wykluczy
  TGP claim jeśli drift wciąż <25%

**Realizm:** PMNS angles są **już falsified** by NuFit 5.3 z 1σ standard.
Future JUNO/DUNE 2030+ z σ improvement do 1-2% **definitivně wykluczy**
TGP claims.

**Werdykt falsifiability test:** ❌ **FAIL** — wszystkie 3 angles już outside
1σ band; "zeroth-order gate" jest accommodating, nie discriminating.

### 2.6 Independent-path cross-validation

**Path 1:** Mixing-operator framework (analog do **κ.1 NUMEROLOGICAL!**)
**Path 2:** ζ.1 PMNS angles "PARTIALLY DERIVED" — ale ζ.1 nie audited
**Path 3:** brak

**Convergence:** Path 1 i Path 2 oba używają tych samych PMNS angle wartości
(1/3, 1/2, K_ν·λ_C²) — *circular* (ζ.1 source dla wartości; ι.1 "explains" via mixing-operator).

**Werdykt independent-path:** ❌ **FAIL** — mixing-operator framework jest
**bezpośrednio dziedziczone z κ.1 NUMEROLOGICAL**; brak independent path.

## 3. Audit gate checklist

```
☑ Phase 0 balance sheet exists
❌ Tautology test FAIL (PMNS outputs nie-derived; charge cascade z θ.1 K_down NUMEROLOGICAL)
❌ Falsifiability test FAIL (3-5σ tensions vs NuFit 5.3; "zeroth-order gate" accommodating)
❌ Independent-path FAIL (mixing-operator framework z κ.1 NUMEROLOGICAL)
❌ Alt-scan 5 candidates falsified — ale alternatives wybrane dla illustrative falsification, nie pełen first-principles physical scan
❌ POST-HOC structural motivation (mixing-operator framework constructed w κ.1, dziedziczone w ι.1)
☑ NIE circular anchor
⚠ Inheriting drift > parent × 5×: drifts 8-15% vs PDG bands 3-4% (heavily inflated)
```

**5 ❌ + 1 ⚠ z 8 ⇒ status max ANSATZ.**

## 4. Klasyfikacja końcowa

| Klasa | Spełnia? | Uzasadnienie |
|-------|----------|--------------|
| DERIVED FULL | NO | 3-5σ tension vs experiment, mixing-operator from κ.1 NUMEROLOGICAL |
| DERIVED CONDITIONAL | NO | conditional na κ.1+θ.1 oba NUMEROLOGICAL |
| STRUCTURAL | NO | claim "DERIVED" untenable z 3-5σ tensions |
| **ANSATZ** | **YES** | structural pattern (mixing-operator), niezweryfikowany field-theoretycznie i już falsified by experiment |
| NUMEROLOGICAL | partial | charge cascade ↔ B²-diff jest algebraicznie OK z forced inputs; ale claims o PMNS są ANSATZ-level |
| TAUTOLOGY | NO | nie definicyjna kasacja |

**Final verdict:** **ANSATZ** (severe downgrade z claim "DERIVED FULL CASCADE 7/7").

**Severity:** **CRITICAL** — ι.1 ma:
- **3-5σ tension w 3 PMNS angles** vs NuFit 5.3 (już falsified by 1σ standard)
- **"Zeroth-order gate <25%"** non-standard accommodating
- **κ.1 NUMEROLOGICAL contagion** (mixing-operator framework dziedziczone)
- **θ.1 K_down NUMEROLOGICAL contagion** (B²_down=61/25 w 3-way cascade)

Reklasyfikacja: ANSATZ (research-track only, NOT registry-DERIVED).

## 5. Comparison ze status oryginalnym

| Element | Original claim | Retrofit verdict |
|---------|----------------|------------------|
| Status | "DERIVED FULL CASCADE 7/7 PASS" | **ANSATZ** (severe downgrade) |
| Counter PREDICTIONS_REGISTRY | +18 (FULL CASCADE) | sub-tests mechanically PASS via accommodating gate; claims "DERIVED" untenable |
| Sub-tests | 7/7 PASS | mechanically PASS only via "zeroth-order gate <25%"; standard 1σ test → 3 of 3 FAIL |
| ζ.1 cascade promotion | "PARTIALLY DERIVED → DERIVED (mixing-operator)" | **REVERSE**: ζ.1 PMNS angles status questionable (zeroth-order accuracy only) |

## 6. Recommended action

- ☐ NO-OP — wymaga severe downgrade
- ☑ **DOWNGRADE w PREDICTIONS_REGISTRY**:
  - ι.1: "DERIVED FULL CASCADE" → **ANSATZ**
  - ζ.1 PMNS: "PARTIALLY DERIVED → DERIVED" rollback do **ANSATZ-level** (zeroth-order)
- ☑ **CRITIQUE w ι.1 folderze**: utworzyć `CRITIQUE_zeroth_order_gate_2026-05-06.md`:
  - Explicit comparison TGP claim vs NuFit 5.3 1σ bands (3-5σ tensions)
  - "Zeroth-order gate <25%" jest accommodating, nie standard physics test
  - Mixing-operator framework dziedziczenie z κ.1 NUMEROLOGICAL
- ☑ **CASCADE_AUDIT**:
  - **ζ.1** — PMNS source i λ_C source — wymaga audit (high priority)
  - **μ.1** PMNS-phase — wymaga audit (analogiczny pattern, juz w high-risk queue)
  - **ν.1** Majorana phase — może mieć ten sam wzorzec
- ☐ CORE_IMPACT — brak (ι.1 wewnętrzny w research/)

## 7. Notes

**Krytyczne odkrycie:**

ι.1 reprezentuje **gorszy** wariant κ.1 wzorca:
- κ.1 produkuje sympy-exact match dla η.1 numerators (multi-candidate selection)
- ι.1 produkuje **values outside experimental band** ale claim "DERIVED" via accommodating gate

To jest **classic ansatz pattern**: structural framework wprowadzona dla
elegant pattern, ale **nie matches precise experimental measurements**.

**Pattern observed (NEW):** "Zeroth-order accuracy + accommodating gate"

| Cykl | Drift vs PDG/NuFit | Author gate | Standard test | Status |
|------|---------------------|-------------|---------------|--------|
| η.2 (α-residual) | 0.000001% | n/a | <0.001% | DERIVED ✓ |
| η.1 (Wolfenstein) | 0.04% | <0.5% | <1σ (=14% for ρ̄) | STRUCTURAL |
| **ι.1 (PMNS)** | **8-15%** | **<25%** | **<1σ (=3-4%)** | **ANSATZ ❌** |

ι.1 używa **najbardziej accommodating gate** — explicit non-standard "zeroth-order gate"
który allows 25% drift. To jest **non-falsifiability w action**.

**Cascade reverse impact:**

ι.1 promotion ζ.1 → "DERIVED (mixing-operator)" jest untenable po retrofit:
- ζ.1 PMNS status questionable
- λ_C anchor (z ζ.1) używane w 3+ cyklach (η.1, ι.1, η.2)
- ζ.1 musi być audited next priority

**Author's caveat acknowledged:**

W oryginalnym Phase2_results.md author **sam pisze**:
> "Caveat preserved: Status PMNS DERIVED via mixing-operator framework
> structurally, but drifts 8.58–15.57% vs NuFit 5.3 inherited z ζ.1
> (zeroth-order). δ_CP phase otwarte JUNO/DUNE 2030+. **Full DERIVED w
> sensie CKM (z numerical lock < 0.05%) czeka na future μ.1/ν.1 cycle**
> z higher-order corrections."

To jest **explicit acknowledgment** że current status NIE jest DERIVED w
standardowym sensie. M03 retrofit formalizes to: **ANSATZ** w taxonomii
CALIBRATION_PROTOCOL.

## 8. Cross-references

- [[../op-iota-charge-pmns-unification/README.md]] — cykl audited
- [[../op-iota-charge-pmns-unification/Phase2_results.md]] — main derivation
- [[retrofit_op-kappa_2026-05-06.md]] — κ.1 retrofit (mixing-operator source)
- [[retrofit_op-theta_2026-05-06.md]] — θ.1 retrofit (B²_down NUMEROLOGICAL)
- [[retrofit_op-eta-wolfenstein_2026-05-06.md]] — η.1 retrofit
- [[../op-zeta-mass-spectrum]] — ζ.1 source dla PMNS i λ_C (NIE-audytowane, high priority)
- [[README.md]] — M03 master plan
- [[audit_log.md]] — running log (added 2026-05-06)
- [[tracker.md]] — status updated do DONE_ANSATZ
- [[../../meta/CALIBRATION_PROTOCOL.md]]
