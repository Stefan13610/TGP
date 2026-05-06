---
title: "CALIBRATION_PROTOCOL — anti-overclaim discipline for new cycles"
date: 2026-05-04
last_updated: 2026-05-06
type: protocol
status: ABSOLUTE BINDING for ALL new cycles claiming any epistemic class (post-2026-05-06)
parent: "[[AUDYT_TGP_2026-05-01.md]]"
related:
  - "[[SUBAGENT_AUDIT_74394a8_2026-05-02.md]]"
  - "[[research/AGENT_PROTOCOL.md]]"
  - "[[PLAN_RESEARCH_WORKFLOW_v1.md]]"
  - "[[CALIBRATION_GATE_ENFORCEMENT.md]]"
  - "[[../research/op-M03-balance-sheet-retrofit-2026-05-06/]]"
tags:
  - meta
  - calibration
  - anti-overclaim
  - protocol
  - audit-gate
  - phase6-enforced
---

# CALIBRATION_PROTOCOL — anti-overclaim discipline

> **Cel:** zapobiec systemic over-claiming wzorca λ.1 / χ.1 / UV.2 (per
> [[SUBAGENT_AUDIT_74394a8_2026-05-02.md]] §3) **i wzorca mixing-operator
> family κ.1 / ι.1 / μ.1** (per M03 retrofit 2026-05-06).
>
> **Trigger:** każdy nowy cykl claiming jakąkolwiek klasę epistemiczną
> (`DERIVED FULL/CONDITIONAL`, `STRUCTURAL`, `ANSATZ`, `NUMEROLOGICAL`,
> `FULL CONVERGENCE`, `KEYSTONE`) musi przejść 2-page balance sheet
> review **przed** committed do PREDICTIONS_REGISTRY/INDEX master ledger.
>
> **Status:** **ABSOLUTE BINDING for ALL new cycles 2026-05-06+** (Phase 6
> enforcement). Previous cycles (greckie litery ε…ψ, M9, M10, M11,
> closure_2026-04-26) były audited retroactively in
> [[AUDYT_TGP_2026-05-01.md]] + [[../research/op-M03-balance-sheet-retrofit-2026-05-06/]]
> M03 retrofit framework.

## ⚠ Phase 6 enforcement update 2026-05-06

Po Phase 1-5 M03 retrofit framework (12 cykli audited, 9 systemic
over-claiming instances confirmed, 5 honest reporting positive examples),
CALIBRATION_PROTOCOL upgrade do **ABSOLUTE BINDING**:

### Rules (post-Phase 6)

1. **NIE może być nowy DERIVED claim bez `Phase0_balance.md` w folderze cyklu.**
   - File MUST exist PRZED commit do master ledger
   - Template: [[../research/op-M03-balance-sheet-retrofit-2026-05-06/template_Phase0_balance.md]]
2. **NIE może być promocja statusu (e.g., PARTIALLY DERIVED → DERIVED) bez explicit cascade audit.**
   - Promotion claims wymagają audit każdego prerequisite cyklu
   - Cascade-aware classification mandatory
3. **NIE może być "constructed criterion" by select winner z multi-candidate set.**
   - Criterion MUST be pre-derived w osobnym axiomatic cycle
   - Mixing-operator post-hoc construction (κ.1 wzorzec) **AUTOMATIC FAIL**
4. **NIE może być "accommodating gate" niestandardowy.**
   - Default falsifiability test: 1σ vs experimental band
   - Custom gate (e.g., "20%" or "25%") wymaga **explicit justification**
   - Brak justification = AUTOMATIC NUMEROLOGICAL classification
5. **NIE może być sympy-rationalization claim "DERIVED" without first-principles.**
   - "X = p/q sympy-exact" wymaga sprawdzenia: forced by axioms vs sympy fitted?
   - Sympy fitting (post-numerical) ≠ first-principles derivation

→ Patrz [[CALIBRATION_GATE_ENFORCEMENT.md]] dla pełnej operational guide.

### Pattern recognition z M03 (negative examples — automatic FAIL)

| Pattern | Examples | M03 verdict |
|---------|----------|-------------|
| Multi-candidate fit z minimum drift selection | UV.2 K_struct, θ.1 K_down | NUMEROLOGICAL |
| Constructed criterion by select C0 | κ.1 mixing-operator denom-num pairing | NUMEROLOGICAL |
| 3-5σ tensions + accommodating "zeroth-order gate" | ι.1 PMNS angles | ANSATZ |
| Drift hardening via fitted corrections z NUMEROLOGICAL source | μ.1 lift factors | NUMEROLOGICAL |
| Anchor borrowed from external (NIE first-principles) | ε.1 137 anchor | STRUCTURAL (limit) |
| Algebraic re-arrangement masquerading as second path | ε.1 F4 chain | NOT independent path |
| Definitional tautology (output kasuje się definicyjnie) | χ.1 G_N | TAUTOLOGY |

### Pattern recognition z M03 (positive examples — pass gate)

| Pattern | Examples | M03 verdict |
|---------|----------|-------------|
| Honest "PARTIAL POSITIVE" + acknowledged limitations | δ.1, δ.2 | STRUCTURAL ✓ |
| Multi-anchor reality acknowledgment | γ.1 (Ω_Λ↔α_s trade-off) | STRUCTURAL ✓ |
| Honest cascade conditionality | η.2 (conditional na ε.1+θ.1) | DERIVED_CONDITIONAL ✓ |
| Honest "PARTIALLY DERIVED" + explicit zeroth-order gate | ζ.1, XS.1 | STRUCTURAL ✓ |
| 2 independent paths z sympy-exact equivalence | η.2 Form A ≡ Form B | DERIVED_CONDITIONAL ✓ |

---

## 1. Trzy klasy epistemiczne

Każdy claim TGP wpisuje się w jedną z czterech klas:

| Klasa | Co znaczy | Wymagana evidence | Promotion gate |
|---|---|---|---|
| **DERIVED** | output wynika **first-principles** z {axioms, prior-LOCKED cycles}, bez post-hoc fittingu i bez circular anchors | Phase 0 balance sheet + sympy LOCK + alt-scan ≥4 falsified ≥3σ + falsifier identified | balance sheet + 5/5 + 7/7 + 6/6 + cross-validate from independent path |
| **STRUCTURAL** | output spełnia algebraic/structural constraint, ale wymaga external anchor lub jest jeden z multi-candidate winnerów | sympy LOCK + alt-scan ≥3 falsified | 5/5 + 6/7 + 5/6 |
| **ANSATZ** | output to wzorzec strukturalny (np. log-conformal mode, threshold form), niezweryfikowany field-theoretycznie | minimal — phase 1 setup OK | (research-track only, NOT registry-DERIVED) |
| **NUMEROLOGICAL OBSERVATION** | output to numeryczna koincydencja w teoretycznym paśmie ≥10× większym niż drift, bez first-principles motywacji | rzetelne reportowanie (drift, band, alt-scan) | (research-track only, NOT registry-DERIVED) |

**Zasady promocji:**
- **NUMEROLOGICAL → ANSATZ:** wymagane independent structural motivation
  (e.g. `2π² = vol(S³)` musi być zaderywowane z TGP ontology, NIE post-hoc).
- **ANSATZ → STRUCTURAL:** wymagane field-theoretic test (sympy LOCK
  z first-principles inputs, NIE tylko algebraic consistency).
- **STRUCTURAL → DERIVED:** wymagane independent-path cross-validation
  (e.g. UV.x M_TGP NIE z M_Pl_PDG anchor, ale z FRG flow + AS NGFP fixed-point).

---

## 2. Phase 0 balance sheet (binding template)

Każdy cykl claiming `DERIVED` musi posiadać `Phase0_balance.md` z polami:

### 2.1 External inputs (PDG, CODATA, observational)

Lista per-cycle z ilością cyfr znaczących i band:
```
- PDG α_em^-1 = 137.035999084(21)         [9 sig figs, 0.15 ppb band]
- PDG M_Z = 91.1876(21) GeV                [4 sig figs, 23 ppm band]
- DESI DR2 w_0 = -0.75 ± 0.10              [DR2 2025-03 arXiv:2503.14738]
- ...
```

### 2.2 Structural axioms (TGP-internal LOCKED)

Lista anchorów które mają **independent LOCK** (sympy diff=0 z innego
cyklu, nie self-reference):
```
- g* = ?? (UV.1 NGFP fixed-point)
- N_A = 500/57 (ξ.1 photon-ring count)
- B²_up = 13/4, B²_down = 61/25 (θ.1 quark Koide)
- E_TGP = 536/75 (ω.2, mechanically from θ.1)
- ...
```

### 2.3 Derived outputs (the cycle claims)

Lista co cykl twierdzi że wyprowadza:
```
- Output 1: G_N(SI)
- Output 2: M_Pl
- Output 3: M_TGP
- ...
```

### 2.4 Tautology test (CRITICAL)

Dla każdego output:
- **Pytanie:** czy output jest wyrażalny jako funkcja wyłącznie
  external inputs i axiomów, **bez** redukcji do tożsamości jednostkowej?
- **Sympy substitution:** podstawić wszystkie axiom relations w wzór output.
  Czy outputs **kasują się tożsamościowo**? Jeśli tak → **TAUTOLOGY**.
- **Przykład χ.1 (FAILED test):** `G_N = g*/(M_TGP²·ξ_grav)` z
  `M_TGP = M_Pl·√(g*/ξ_grav)` → `G_N = g*/(M_Pl²·g*) = 1/M_Pl²`.
  g* i ξ_grav się kasują. **TAUTOLOGY** → status max **ANSATZ**, NIE DERIVED.

### 2.5 Falsifiability test (CRITICAL)

Dla każdego output:
- **Pytanie:** czy istnieje wartość axiomu lub external input która
  **wykluczyłaby** match? Jeśli każdy axiom redukuje się do "fitting noise"
  w experimental band, output jest **non-falsifiable** → status max
  **NUMEROLOGICAL OBSERVATION**.
- **Przykład UV.2 (FAILED test):** drift `0.29%` < theoretical band `10–30%`
  M_GUT 2-loop. K_struct mógłby być 70-130% N_A·2π² i wciąż drift OK.
  **NON-FALSIFIABLE** → status max **NUMEROLOGICAL**, NIE DERIVED.
- **Przykład η.2 (PASSED test):** α-residual = 9/250 sympy-exact diff=0
  vs CODATA `α^-1 = 137.036` 9 sig figs. Falsifier: jeśli α^-1_CODATA
  ≠ 137 + 9/250, η.2 FAILS. **FALSIFIABLE** → DERIVED OK.

### 2.6 Independent-path cross-validation (CRITICAL for DERIVED)

- **Pytanie:** czy istnieje **niezależna ścieżka** od axiomów do output
  która daje ten sam result?
- **Przykład M9.1''**: β_PPN = 1 z (a) "master formula" + c₂=-1
  + (b) closure_2026-04-26 calibration via T-FP. **Two paths OK** → DERIVED.
- **Przykład UV.2** (FAILED): jedyna ścieżka to `K_struct = (M_Pl/M_GUT)·√(g*/N_A)`
  empiryczny. **NIE niezależna ścieżka** → max NUMEROLOGICAL.

---

## 3. Audit gate — checklist dla agenta przed registry write

Przed dopisaniem `DERIVED FULL` / `FULL CONVERGENCE` / `KEYSTONE` w
PREDICTIONS_REGISTRY:

```
☐ Phase 0 balance sheet exists (Phase0_balance.md w folderze cyklu)
☐ Tautology test PASS (sympy substitution → outputs nie kasują się)
☐ Falsifiability test PASS (existing experimental band < 5× drift claim)
☐ Independent-path cross-validation PASS (≥2 paths convergent)
☐ Alt-scan ≥4 candidates with ≥3σ discrimination (NIE drift-min aesthetic)
☐ NIE używane post-hoc structural motivations (np. "2π² = vol(S³)" musi
   być pre-derived)
☐ NIE circular anchor (output jako funkcja samego siebie po substitution)
☐ NIE inheriting drift > parent cycle drift × 5× (cascade discipline)
```

**Brak choćby jednego ☐ → max status STRUCTURAL.** **Tautology lub
falsifiability FAIL → max status ANSATZ lub NUMEROLOGICAL OBSERVATION.**

---

## 4. Self-correction discipline

Cykl który po-promocyjnie ujawnia tautologię/circular anchor/post-hoc
fitting wykonuje **mark-as-unproven** (NIE rollback, NIE delete) w
3 krokach:

1. **CRITIQUE_<issue>_<date>.md** w folderze cyklu z explicit algebraic
   investigation (sympy substitution showing tautology, lub fitting
   evidence) — wzór: [[SUBAGENT_AUDIT_74394a8_2026-05-02.md]] sec 2.1, 2.2.
2. **Phase3_results.md** — verdict downgrade w YAML + opening blockquote
   z linkiem do CRITIQUE; sub-tests PASS preserved.
3. **PREDICTIONS_REGISTRY.md / INDEX.md** — REVISION block + per-row
   epistemic status table; counter mark `effective uncontested`.

**Sub-tests PASS NIE są usuwane** — one są mechanicznie poprawne jako
algebraic identities. Tylko **interpretacja statusu** downgraded.

---

## 5. Cross-references

- [[SUBAGENT_AUDIT_74394a8_2026-05-02.md]] — root-cause exemplar
  (chi.1 G_N tautology, UV.2 K_struct numerologia)
- [[research/op-chi1-newton-constant-derivation/CRITIQUE_circular_anchor_2026-05-02.md]]
- [[research/op-uv2-mtgp-absolute-scale/CRITIQUE_repackaged_circularity_2026-05-02.md]]
- [[research/op-omega2-axion-coupling-lock/AUDIT_omega2_2026-05-04.md]]
  (positive example — cascade z θ.1 OK, no critique)
- [[research/op-omega3-axion-decay-constant/AUDIT_omega3_2026-05-04.md]]
  (cascade-conditional example — algebra OK, magnitude conditional)

---

**Status:** BINDING 2026-05-04+. **Apply przed each new cycle's promotion claim.**
