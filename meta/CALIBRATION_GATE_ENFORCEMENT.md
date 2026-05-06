---
title: "CALIBRATION_GATE_ENFORCEMENT — operational gate dla future cykli (Phase 6)"
date: 2026-05-06
type: operational-protocol
status: ABSOLUTE BINDING for ALL new cycles 2026-05-06+
parent: "[[CALIBRATION_PROTOCOL.md]]"
related:
  - "[[research/AGENT_PROTOCOL.md]]"
  - "[[../research/op-M03-balance-sheet-retrofit-2026-05-06/template_Phase0_balance.md]]"
  - "[[../research/op-M03-balance-sheet-retrofit-2026-05-06/Phase5_registry_refactor_draft.md]]"
tags:
  - meta
  - calibration
  - gate
  - enforcement
  - phase6
  - operational
---

# CALIBRATION_GATE_ENFORCEMENT — operational guide

> **Status:** ABSOLUTE BINDING for ALL new cycles 2026-05-06+
> (Phase 6 of M03 retrofit framework).
>
> **Purpose:** Operational guide dla future cycles + future agents. After
> M03 retrofit revealed systemic over-claiming pattern w 9 cyklach
> (4 z 74394a8 + 5 pre-74394a8), Phase 6 establishes **mandatory pre-commit gate**
> that future cycles MUST pass.

## Executive summary

**Trzy kluczowe zasady (post-2026-05-06):**

1. **Phase0_balance.md PRZED commit** — żaden new cycle claim status bez
   completed Phase 0 balance sheet w folderze cyklu
2. **Honest classification > Optimistic claim** — STRUCTURAL z disclosure
   limitations > "DERIVED FULL" z cascade dependencies
3. **Independent path test PHYSICAL, NIE algebraic** — algebraic
   re-arrangement same equations ≠ independent derivation

## Mandatory workflow dla future cycles

### Krok 1: Folder creation

Każdy nowy cykl `research/op-X-Y-Z/` MUST create:

```
research/op-X-Y-Z/
├── README.md                # cycle description (existing convention)
├── Phase0_balance.md        # MANDATORY (Phase 6 enforcement)
├── Phase1_setup.md          # methodology
├── Phase1_results.md        # results
├── (Phase 2, 3 if applicable)
├── FINDINGS.md              # exportable
├── NEEDS.md                 # remaining open items
└── CRITIQUE_*.md            # if self-correction needed
```

**`Phase0_balance.md` jest MANDATORY**, nie optional.

### Krok 2: Phase0_balance.md template

Skopiować z [[../research/op-M03-balance-sheet-retrofit-2026-05-06/template_Phase0_balance.md]]:

Zawartość MANDATORY 6 sections:

1. **External inputs** — PDG/CODATA/observational (z sig figs + bands)
2. **Structural axioms** — TGP-internal LOCKED inputs (z source cycle)
3. **Derived outputs** — co cykl twierdzi że wyprowadza
4. **Tautology test** — sympy substitution → outputs nie kasują się
5. **Falsifiability test** — concrete falsifier + band check
6. **Independent-path cross-validation** — ≥2 PHYSICAL paths convergent

### Krok 3: Pre-commit gate checklist

Przed jakikolwiek registry write (`PREDICTIONS_REGISTRY.md` update):

```
☐ Phase0_balance.md exists (in cycle folder)
☐ Tautology test PASS (sympy substitution → outputs nie kasują się definicyjnie)
☐ Falsifiability test PASS (existing experimental band > 5× drift claim)
☐ Independent-path cross-validation PASS (≥2 physical paths convergent)
☐ Alt-scan ≥4 first-principles candidates with ≥3σ discrimination
   (NIE algebraic identities z TGP axioms — wymaga physical models)
☐ NIE post-hoc structural motivations
☐ NIE constructed criterion by select winner
☐ NIE circular anchor (output jako funkcja samego siebie)
☐ NIE inheriting drift > parent cycle drift × 5×
☐ Custom falsifiability gate (>1σ) explicit JUSTIFIED
```

**Brak choćby jednego ☐ → max status STRUCTURAL.**

**Specific FAIL conditions → automatic max status:**

| Failed criterion | Max possible status |
|------------------|---------------------|
| Tautology test FAIL | ANSATZ |
| Falsifiability test FAIL | NUMEROLOGICAL OBSERVATION |
| Independent-path FAIL | STRUCTURAL |
| Alt-scan inadequate (algebraic only) | STRUCTURAL |
| Constructed criterion | NUMEROLOGICAL |
| Multi-candidate selection | NUMEROLOGICAL (UV.2 wzorzec) |
| 3+ σ tension vs experiment | ANSATZ |
| Cascade z NUMEROLOGICAL source | inherits NUMEROLOGICAL |

### Krok 4: PREDICTIONS_REGISTRY entry

Per CALIBRATION_PROTOCOL §3, registry entry musi include:

```markdown
**Updated <date> (cycle X program END):** <counter> cumulative
(+ X.Phase1 N + X.Phase2 N + X.Phase3 N = total).
**Phase0 balance sheet:** [link to Phase0_balance.md]
**Classification:** <DERIVED FULL | DERIVED CONDITIONAL | STRUCTURAL | ANSATZ | NUMEROLOGICAL>
**Justification:** <1-2 sentences why this classification>
**Cascade dependencies:** <list of prerequisite cycles>
**Falsifier:** <concrete experimental test + horizon>
```

## Promotion path (post-Phase 6)

Status promotions wymagają **explicit gate re-pass**:

```
research-track ANSATZ
    ↓ [requires: independent structural motivation]
NUMEROLOGICAL OBSERVATION
    ↓ [requires: field-theoretic test]
STRUCTURAL
    ↓ [requires: independent-path cross-validation (2+ physical paths)]
DERIVED CONDITIONAL
    ↓ [requires: prerequisite audits OK + sympy LOCK]
DERIVED FULL
```

**KAŻDA promocja wymaga osobnego balance sheet review.**

**Auto-rejection scenarios:**

- Status promotion via "cascade" (e.g., "post-X.1 → DERIVED refined²") bez
  audit X.1 prerequisite — REJECTED
- Status promotion na podstawie sub-tests count (e.g., "7/7 PASS") bez
  Phase 0 balance sheet — REJECTED
- "FULL CASCADE" claim z multi-candidate selection — REJECTED

## Pattern recognition reference

### Negative examples (FAIL gate, status downgrade)

**κ.1-style "constructed criterion":**
```
Detect: cycle introduces new "rule" by select unique value z multi-candidate
        sympy-exact set
Action: AUTOMATIC FAIL alt-scan criterion → max NUMEROLOGICAL
```

**ι.1-style "accommodating gate":**
```
Detect: cycle uses "zeroth-order gate >5%" or similar non-standard band
        AND values outside experimental 1σ
Action: AUTOMATIC FAIL falsifiability → max ANSATZ
```

**μ.1-style "cascade fitting":**
```
Detect: cycle introduces correction (1-X) lub similar that exactly compensates
        upstream zeroth-order drift
Action: AUTOMATIC FAIL independent-path → cascade contagion
```

**θ.1 K_down-style "multi-candidate within band":**
```
Detect: ≥3 sympy candidates within experimental band, author selects
        minimum-drift winner without structural argument for selection
Action: AUTOMATIC FAIL falsifiability → max NUMEROLOGICAL
```

**ε.1-style "borrowed anchor":**
```
Detect: anchor borrowed z external constant (e.g., 137 z α_fine PDG) bez
        first-principles TGP derivation
Action: status max STRUCTURAL (NIE DERIVED — wymaga independent derivation
        anchor)
```

**χ.1-style "definitional tautology":**
```
Detect: sympy substitution of axioms cancels output do trivial identity
        (e.g., G_N = 1/M_Pl² po substitucji M_TGP definitions)
Action: AUTOMATIC FAIL tautology → max ANSATZ
```

### Positive examples (PASS gate)

**δ.x-style "honest PARTIAL POSITIVE":**
- Author explicit acknowledges incomplete derivation
- Multi-candidate scenario disclosed
- Status reported as "PARTIAL POSITIVE" or "Level B"
- → PASS as STRUCTURAL or STRUCTURAL_CONDITIONAL

**γ.1-style "multi-anchor reality":**
- Author identifies multiple structural alternatives
- Trade-off acknowledged (e.g., Ω_Λ↔α_s)
- Pure structural form distinguished z phenomenological corrections
- → PASS as STRUCTURAL

**η.2-style "two physical independent paths":**
- 2 niezależne forms (Form A, Form B) sympy-exact equivalent
- Each form uses different axioms
- Concrete falsifier + future experiment
- → PASS as DERIVED_CONDITIONAL

## Self-correction protocol (post-promotion FAIL)

Cykl który po-promocyjnie ujawnia tautologię/circular anchor/post-hoc
fitting wykonuje **mark-as-unproven** w 3 krokach:

1. **CRITIQUE_<issue>_<date>.md** w folderze cyklu z explicit algebraic
   investigation
2. **Phase3_results.md** — verdict downgrade w YAML + opening blockquote
   z linkiem do CRITIQUE
3. **PREDICTIONS_REGISTRY.md / INDEX.md** — REVISION block + per-row
   epistemic status table; counter mark `effective uncontested`

**Sub-tests PASS NIE są usuwane** — one są mechanicznie poprawne jako
algebraic identities. Tylko **interpretacja statusu** downgraded.

→ Reference: `research/op-chi1-newton-constant-derivation/CRITIQUE_circular_anchor_2026-05-02.md`,
`research/op-uv2-mtgp-absolute-scale/CRITIQUE_repackaged_circularity_2026-05-02.md`

## Cascade-aware classification

**Cascade rule:** Cykl X używający anchor/output Y ma **maksimum status**
zależny od status Y:

```
Y = TAUTOLOGY     → X max ANSATZ
Y = NUMEROLOGICAL → X max ANSATZ
Y = ANSATZ        → X max STRUCTURAL
Y = STRUCTURAL    → X max DERIVED_CONDITIONAL
Y = DERIVED_*     → X may achieve DERIVED_FULL
```

**Implementation:** Phase0_balance.md §"2.2 Structural axioms" MUST list
each input z status of source cycle. Audit gate cross-references status.

## Verification checklist dla future M03 audits (Phase 3-4)

Future agents continuing M03 retrofit (Phase 3 medium-risk + Phase 4
low-risk) MUST follow:

### Per-cycle audit checklist

```
1. ☐ Read [[research/op-M03-balance-sheet-retrofit-2026-05-06/resume_protocol.md]]
2. ☐ Check tracker.md status of cycle (PENDING / IN_PROGRESS / DONE)
3. ☐ Read cycle Phase1-3 results
4. ☐ Apply 9-pattern recognition (z patterns observed sekcja audit_log)
5. ☐ Apply 8-criteria audit gate checklist
6. ☐ Classify: DERIVED FULL / CONDITIONAL / STRUCTURAL / ANSATZ / NUMEROLOGICAL / TAUTOLOGY
7. ☐ Check cascade dependencies (prerequisite cycles)
8. ☐ Document w retrofit_<cycle>_<date>.md (use template_Phase0_balance.md)
9. ☐ Update tracker.md status
10. ☐ Update audit_log.md z entry + new patterns observed
```

### Per-session checklist

```
1. ☐ Audit log SESSION_START entry
2. ☐ Audit 1-3 cycles per session (avoid context overload)
3. ☐ Update tracker incrementally
4. ☐ Audit log SESSION_END entry z summary + next priorities
5. ☐ Don't proceed to Phase 5 (registry refactor) until Phase 3+4 complete
```

## Hooks integration (optional, future)

For automated enforcement, settings.json or git pre-commit hooks could
verify:

```bash
# Pseudo-code dla pre-commit hook
for new_cycle in research/op-*; do
  if [[ ! -f "$new_cycle/Phase0_balance.md" ]]; then
    echo "ERROR: $new_cycle missing Phase0_balance.md (CALIBRATION_PROTOCOL §6)"
    exit 1
  fi
done
```

→ Implementation deferred (manual enforcement working post-2026-05-06).

## Reference: M03 retrofit lessons learned (12 cycles, 5 honest + 7 over-claiming)

| Cykl | Original | M03 verdict | Pattern |
|------|----------|-------------|---------|
| η.2 | "FULL CASCADE 18/18" | DERIVED_CONDITIONAL | Honest cascade |
| ε.1 | "PARTIALLY DERIVED" | STRUCTURAL | Borrowed 137 anchor |
| θ.1 K_up | "PARTIALLY DERIVED" | STRUCTURAL | Universal pattern OK |
| θ.1 K_down | "STRUCTURAL refined" | NUMEROLOGICAL | Multi-candidate (UV.2 wzorzec) |
| η.1 | "DERIVED (refined²)" | STRUCTURAL | Numerator κ.1 contagion |
| **κ.1** | "DERIVED FULL CASCADE 7/7" | **NUMEROLOGICAL** | Constructed criterion |
| **ι.1** | "DERIVED FULL CASCADE 7/7" | **ANSATZ** | 3-5σ + accommodating gate |
| δ.1 ★ | "PARTIAL POSITIVE" | STRUCTURAL | Honest reporting |
| δ.2 ★ | "Level B PARTIAL" | STRUCTURAL_COND | Honest reporting |
| γ.1 ★ | "POSITIVE z H5" | STRUCTURAL | Multi-anchor honest |
| XS.1 ★ | "PARTIALLY DERIVED" | STRUCTURAL | Honest cross-sector |
| **μ.1** | "DERIVED FULL, 8 free→0 free" | **NUMEROLOGICAL** | Drift hardening fitting |
| ζ.1 ★ | "PARTIALLY DERIVED" | STRUCTURAL | Honest zeroth-order |

**Pattern:** 5 z 12 (42%) honest = pass gate. 3 z 12 (25%) "DERIVED FULL CASCADE" claims → NUMEROLOGICAL/ANSATZ severe downgrade.

**Lesson:** Honest reporting prevents over-claiming. CALIBRATION_PROTOCOL §2.4-2.6
tests are anti-overclaim mechanism. Phase 6 enforcement makes this MANDATORY.

## Cross-references

- [[CALIBRATION_PROTOCOL.md]] — original protocol (now ABSOLUTE BINDING post-Phase 6)
- [[../research/op-M03-balance-sheet-retrofit-2026-05-06/]] — retrofit framework + 12 examples
- [[../research/op-M03-balance-sheet-retrofit-2026-05-06/template_Phase0_balance.md]] — mandatory template
- [[../research/op-M03-balance-sheet-retrofit-2026-05-06/resume_protocol.md]] — anti-duplication guide
- [[SUBAGENT_AUDIT_74394a8_2026-05-02.md]] — root-cause exemplar (chi.1, UV.2)
- [[research/AGENT_PROTOCOL.md]] — research workflow (Phase 6 update)
- [[../audyt/M03_balance_sheet_missing/POST_ACTION_UPDATE_2026-05-06.md]]
