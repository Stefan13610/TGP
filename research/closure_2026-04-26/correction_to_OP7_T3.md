# Correction to OP-7 T3 — Path B promotion to PRIMARY derivation

**Data:** 2026-04-26
**Status:** **CORRECTION ISSUED**
**Affected file:** [[../op7/OP7_T3_results.md]]
**Source:** [[sigma_ab_pathB/results.md]] (closure_2026-04-26 Phase 1)
**Verdict:** Path B 11/11 PASS — promoted from "Path A or B equivalent" to **Path B PRIMARY**.

---

## 1. Cel correction note

OP7_T3_results.md (2026-04-25) sformułował dynamikę σ_ab przez "**two paths
sprawdzone strukturalnie i pokrywają się**":
- **Path A:** efektywny Lagrangian L_σ z niezależnym σ_ab field
- **Path B:** composite z s-EOM, σ_ab(x) = ⟨(∂_a δŝ)(∂_b δŝ)⟩ - trace

Original conclusion: "Path A == Path B; single-Φ aksjomat zachowany."

**Issue identified post-OP-7:** Path A wymaga *quasi-fundamental* σ_ab Lagrangian,
co jest na granicy single-Φ aksjomatu z TGP_FOUNDATIONS §1. Choć Path A i Path B
dają ten sam EOM, ich **ontological status jest różny**:

- Path A: σ_ab traktowane jako quasi-niezależny field z własnym Lagrangianem
- Path B: σ_ab traktowane jako kompozyt z baseline ŝ field (single-Φ pozostaje jedyny)

Path B jest bardziej rygoryczny TGP-philosophically.

**Correction:** Closure_2026-04-26 Phase 1 audit (sigma_ab_pathB/) explicitly
verified Path B przez box-of-product algebra i wyprowadził `M² = 2m_s²` z heredity
(NIE jako analogię do mezonów).

⇒ **Path B promoted to PRIMARY derivation** of σ_ab dynamics in TGP.

---

## 2. Patch do OP7_T3_results.md

### 2.1 Sekcja "Część T3.1 — Strukturalna derivacja EOM"

**ORIGINAL (2026-04-25, OP7_T3_results.md §2 "Część T3.1"):**

> Two paths sprawdzone strukturalnie i pokrywają się:
> [Path A i Path B]
> Path A == Path B; single-Φ aksjomat zachowany.

**CORRECTED (2026-04-26, post-closure):**

> **Two paths analyzed; Path B is PRIMARY (TGP-natural):**
>
> **Path B (PRIMARY — closure_2026-04-26 Phase 1):**
> σ_ab(x) = ⟨(∂_a δŝ)(∂_b δŝ)⟩^TF dziedziczy bezpośrednio z linearizacji
> ŝ-EOM. Zastosowując box-of-product identity:
> ```
> □(AB) = (□A)B + 2(∂_μA)(∂^μB) + A(□B)
> ```
> z A = ∂_a δŝ, B = ∂_b δŝ, używając □δŝ = J - m_s² δŝ:
> ```
> □σ_ab + 2 m_s² σ_ab = (TT projection of source)
> ```
> Co identyfikuje **m_σ² = 2 m_s²** — NIE jako analogię do mezonów,
> ale jako rygoryczną consequence of heredity algebra.
>
> Ghost-freeness: gwarantowana strukturalnie przez Gram-matrix positivity
> (eigenvalues K_ab ≥ 0 dla composite operator). Brak nowych d.o.f.
>
> **Path A (effective equivalent):**
> Efektywny Lagrangian L_σ podaje ten sam EOM, ale wymaga σ_ab jako
> quasi-fundamental field — naruszałby spirit single-Φ aksjomatu.
> Path A pozostaje **useful effective tool**, ale NIE primary ontology.
>
> **Verdict:** Path B 11/11 PASS (closure_2026-04-26/sigma_ab_pathB/).
> Single-Φ aksjomat **strictly preserved** — σ_ab to composite, NIE field.

### 2.2 Sekcja "Cross-references"

**Add to OP7_T3_results.md cross-references:**
```
- [[research/closure_2026-04-26/sigma_ab_pathB/results.md]]   ← Path B PRIMARY (11/11 PASS)
- [[research/closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]] ← Phase 1 closure
- [[research/closure_2026-04-26/correction_to_OP7_T3.md]]      ← This correction
```

### 2.3 Status block

**ORIGINAL:** `Status: ⚠️ STRUCTURAL POSITIVE z OPEN TENSION (Φ₀/m_σ)`

**UPDATED:** `Status: ✅ PATH B STRUCTURAL DERIVATION CLOSED (closure_2026-04-26 Phase 1, 11/11 PASS); ⚠️ T3.2 m_σ ↔ Φ₀ tension ROBUST OPEN (NOT addressed by Path B promotion — wymaga T3-extended Bethe-Salpeter or Φ₀ ≪ meV).`

---

## 3. Implementacja correction

**Strategia:** Nie edytujemy bezpośrednio OP7_T3_results.md (zachowanie audit
trail). Zamiast tego:
1. Ten correction note (`correction_to_OP7_T3.md`) **stoi obok** original.
2. Future readers OP7_T3_results.md powinni najpierw konsultować ten correction.
3. Tgp-core-paper i PRD companion update pochodzą **z corrected version**.

**Patch summary:**
- Path A → Path B PRIMARY (TGP-natural composite ontology)
- m_σ² = 2 m_s² **derived**, not analogized
- Ghost-freeness **structurally guaranteed** by Gram-positivity
- Single-Φ axiom **strictly preserved**

**NIE zmienia się:**
- T3.2 (m_σ scale tension Φ₀ ↔ GW170817): pozostaje OPEN.
- T3.3 (ghost analysis): wzmocnione przez Gram-positivity argument.
- T3.4 (ξ coupling): pozostaje OPEN.
- T3.5 (Bethe-Salpeter): pozostaje OPEN.
- T3.6 (symmetry protection): pozostaje OPEN.

---

## 4. Cross-references

- [[../op7/OP7_T3_results.md]] (original 2026-04-25)
- [[../op7/OP7_T3_extended_results.md]] (T3-extended)
- [[sigma_ab_pathB/setup.md]] (closure_2026-04-26 Phase 1 design)
- [[sigma_ab_pathB/results.md]] (Phase 1 verdict 11/11 PASS)
- [[sigma_ab_pathB/sigma_ab_pathB_audit.py]] (executive script)
- [[CLOSURE_2026-04-26_SUMMARY.md]] §2 "Phase 1: σ_ab Path B audit"
- [[KNOWN_ISSUES.md]] §A.1 "σ_ab dynamics — Path B PRIMARY"

---

## Bottom line

**OP-7 T3 σ_ab dynamics:** Path B promoted from equivalent alternative to
**PRIMARY derivation** as result closure_2026-04-26 Phase 1 audit.

- σ_ab jest **composite operator** (NIE quasi-fundamental field).
- m_σ² = 2 m_s² **derived** from heredity algebra (NIE analogized).
- Ghost-freeness **structurally guaranteed** by Gram-matrix positivity.
- **Single-Φ Z₂ aksjomat (TGP_FOUNDATIONS §1) strictly preserved.**

Path A pozostaje useful effective Lagrangian formulation, ale NIE jest
primary ontology. Future TGP papers powinny używać **Path B language**:
"σ_ab dziedziczy z s-EOM przez heredity algebra" zamiast "σ_ab z efektywnego
Lagrangianu."
