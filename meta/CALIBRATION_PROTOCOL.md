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

## 4.4 BD-drift audit (added 2026-05-10, BINDING)

### 4.4.1 — Trigger i scope

**Każdy cycle Phase FINAL** (post-2026-05-10 cycles, plus retroactive audit cycles dla
pre-2026-05-10 cykli z BD-form formulami) MUSI spawn `bd-drift-audit` subagent przed wystawieniem
verdict `STRUCTURAL DERIVED`.

**Cel:** wykryć systemic BD-drift (Brans-Dicke / scalar-tensor / Horndeski translation) w
cycle outputs zanim verdict propaguje cascadowo do downstream cykli.

**Scope:**
- Cykle dotyczące gravity / inertia / momentum / mass / GW sektora
- Cykle inheriting LOCKs z pre-2026-05-10 cykli (heavily suspect)
- Cykle używające `m_Φ`, `g_eff`, `Φ-EOM`, `T^μν`, `δΦ propagator`, `Yukawa`, `BD ω`,
  `scalar-tensor`, `Horndeski` (search literally tych terms)

**Wyłączenia:**
- Pure mathematical sympy proofs (algebraic identities, symbol manipulation)
- Documentation-only updates (no new derivations)
- Pure numerical simulations z jasno zdefiniowaną TGP-native equation

### 4.4.2 — Subagent prompt template

```
Subagent task: BD-drift audit dla cycle <CYCLE_NAME>

Context: Cycle <CYCLE_NAME> Phase <X> outputs (read these files):
- <CYCLE_PATH>/Phase<X>_results.md
- <CYCLE_PATH>/Phase<X>_sympy.py
- <CYCLE_PATH>/Phase<X>_sympy.txt
- <CYCLE_PATH>/README.md

Your task:
1. Read TGP_NATIVE_COMPUTATIONAL_PATTERNS.md §1 (ASK-RULE), §2 (patterns), §3 (red flags),
   §4 (form-meaning mapping)
2. Audit cycle outputs dla:
   (a) §3 red flags — list every detected red flag z line numbers
   (b) §4 form-meaning mismatch — every BD-form formula bez explicit annotation
   (c) §1 ASK-RULE triggers — moments gdy agent powinien był pytać user'a but didn't
   (d) §2 missing patterns — gdzie cycle używa std-physics derivation zamiast TGP-native pattern
3. Generate report `BD_DRIFT_AUDIT_<DATE>.md` z:
   - Section 1: detected drifts (list)
   - Section 2: severity classification (LOW/MEDIUM/HIGH/CRITICAL per drift)
   - Section 3: recommendation (PASS / CONDITIONAL / RECOMMEND_AMENDMENT / RECOMMEND_HALT)
   - Section 4: specific fixes per drift

Constraints:
- Read-only, NIE modify cycle files
- Max 3-4 file reads
- Concise report (under 800 words)
- If no drifts found, explicit state "NO BD-DRIFT DETECTED" w Section 1

Return: report file path + 2-sentence summary.
```

### 4.4.3 — Verdict consequences

| Drift severity | Subagent recommendation | Cycle verdict consequence |
|---|---|---|
| NONE | PASS | Cycle proceeds z planned verdict |
| LOW (1-2 minor BD-form formulas, easily annotated) | CONDITIONAL | Cycle adds explicit annotations + proceeds |
| MEDIUM (multiple §3 red flags, missing §2 patterns) | RECOMMEND_AMENDMENT | Cycle Phase X reopened dla TGP-native re-derivation; OR explicit downgrade do `STRUCTURAL_CONDITIONAL z BD-drift flag` |
| HIGH (systemic Φ-quantum carrier framing, fixed m_Φ, postulated mechanisms) | RECOMMEND_HALT | Cycle classification → `STRUCTURAL_CONDITIONAL_HALT z BD-drift`; spawn dedicated TGP-native re-derivation cycle |
| CRITICAL (verdict fundamentally relies on BD-mode physics) | RECOMMEND_HALT + cascade audit | Cycle DOWNGRADED; spawn audit cykli wszystkich downstream cykli inheriting LOCKs |

### 4.4.4 — Adversarial pattern continuation

BD-drift audit jest **rozszerzeniem** existing adversarial verification protocol (§4.3).
Pattern adversarial = każda sesja ma at least 1 audit subagent działający niezalezne od
primary cycle work.

**Demonstrated value (do 2026-05-10):** 5× pre-2026-05-10 catch (sphere-avg error, factor-4
ξ_eff gap, Channel B Yukawa flag, Yukawa audit verdict, m_Φ verification ruling out mech iii)
+ 1× meta-layer catch (BD-drift discovered post-recovery-V Phase 1 → triggered TGP_NATIVE_COMPUTATIONAL_PATTERNS).

### 4.4.5 — Fallback dla absence of subagent capability

Jeśli session bez ability spawn subagent (np. resource constraints, user-request limits):
agent MUSI **manually** wykonać §4.4.2 audit checklist on own outputs przed Phase FINAL,
i explicit document audit results w `Phase<FINAL>_results.md` §X.X "Self-audit BD-drift".

**Honest fallback:** "Self-audit jest weaker niż independent subagent audit." Future session
SHOULD re-run independent audit gdy capability available.

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
- [[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] — anti-BD-drift binding protocol (§1 ASK-RULE,
  §2 patterns, §3 red flags, §4 form-meaning mapping, §5 pre-flight checklist) referenced
  by §4.4 BD-drift audit

---

## 6. R1+R2+R3 two-tier discipline addendum (LOCKED 2026-05-22, post R2_PASS)

**Authorization:** Methodology innovation R1+R2+R3 introduced FFS pre-screening 2026-05-19 §6;
operationally validated dwukrotnie:
- **First operational test (FFS Phase 4 2026-05-20):** R3 trigger active z 1/3 evidence lines
  → axiom NOT accepted (rejection working correctly)
- **Second operational test (CE-H Poziom β Phase FINAL 2026-05-21):** R3 trigger 3/3 evidence
  lines confirmed → CE-H accepted as structural feature (acceptance with caveats)
- **R2 audit cycle (2026-05-22):** R2_PASS (8/9 CLOSED + 1/9 DEFERRED + 0 ESCALATED)
  → propagation z R1+R2+R3 §6 addendum candidate confirmed for CALIBRATION_PROTOCOL §3
  inclusion

### 6.1 — Rule R1 (research-tier permissive flagging)

W trakcie cyklu badawczego nowe strukturalne elementy NIE są zakazane z góry. **Każda nowa
struktura jest *flagged*** w inventory test z dokładnym opisem:
- **Derived:** logicznie z istniejących foundations bez nowych postulatów
- **Reinterpreted:** istniejąca struktura w nowej roli (ta sama mathematical object)
- **Flagged-new:** nowa struktura wymagająca uzasadnienia w integration audit

**Uzasadnienie:** Top-down exploration czasami niezbędna — bottom-up derivation z minimal
axioms może być impossible nawet jeśli mechanizm fizycznie real. Rygor R1 zapobiega
self-strangulation.

### 6.2 — Rule R2 (integration audit gate)

Przed włączeniem jakiegokolwiek flagged-new wyniku do rdzenia TGP (TGP_FOUNDATIONS,
STATE.md as binding), **każda flagged-new structure** MUSI przejść osobny *integration audit*:

- **Necessity test:** czy struktura jest *absolutnie niezbędna* dla obserwowalnej prediction?
- **Derivability check:** czy w międzyczasie pojawił się sposób derive jej z foundations?
- **Alternative path search:** czy istnieje równoważna formulation NIE wymagająca tej struktury?

**Audit cycle pattern:** osobny doc post-cycle `op-R2-integration-audit-XX/`; NIE część
cyklu badawczego — osobny etap.

### 6.3 — Rule R3 (multi-line convergence ≥3 threshold)

Jeśli flagged-new structure NIE może być derived ani eliminated po R2 audit, i hipoteza
włączenia do rdzenia jako *new axiom* jest rozważana, MUSI spełniać:

**R3 threshold:** **≥3 niezależne lines of evidence** (każda **pre-zarejestrowana** *przed*
observation) konwergujące na tej samej strukturze.

**Operacyjna walidacja R3 (dwukrotnie):**
- **FFS Phase 4 (R3 rejection):** 1/3 lines → axiom NOT accepted ✓
- **CE-H Poziom β Phase FINAL (R3 acceptance):** 3/3 lines → CE-H accepted as structural
  feature (NIE nowy axiom — konsekwencja S05+Z₂+U(1)+RP² ontologii) ✓

### 6.4 — Anti-Lakatos preservation pod R1+R2+R3

Trzy bezpieczniki:
1. R1 wymaga flagging każdej nowej struktury → no hidden additions
2. R2 integration audit jest gate przed core inclusion → każda flagged structure
   przechodzi necessity check
3. R3 multi-line convergence ≥3 niezależne pre-registered evidence lines konwergujące →
   no post-hoc convergence

**Krytyczne:** R3 acceptance derivation chain MUST be explicit verified during R2 audit
(NIE implicit axiom drift). R2 audit cycle 2026-05-22 §8 Phase 0 pre-registered + Phase 2
verified CE-H R3 acceptance chain: S05+Z₂ ontology → kink stability → asymmetric configuration
requires bulk Phi ≠ 0 → ⟨Phi⟩_bg structural feature.

---

## 3.6 — Analytical pre-derivation step (LOCKED 2026-05-22, BINDING per R1-1 audit verdict)

**Source:** R1-1 audit verdict (op-R2-integration-audit-CE-H-FFS-2026-05-22 Phase 3) —
pattern observed twice (CE-H Phase 3 T_P3_2 + FFS pre-screening T7 q=1 implicit) where
pre-registration used **heuristic intuition** instead of **analytical derivation**, leading
to formally failed pre-registered thresholds despite substantive 1% accuracy of actual results.

### 3.6.1 — Requirement (BINDING for falsifiers z numerical thresholds)

Dla każdego pre-registered falsifier specifying numerical threshold for derived quantity
(e.g., "fitted decay rate match m within 10%", "ratio within factor 2"), pre-registration
MUST include:

(a) **Analytical derivation** of expected value from underlying ansatz/field theory
(b) **Symbolic computation** demonstrating expected numerical value
(c) **Documentation w Phase 0** alongside numerical threshold

### 3.6.2 — Forbidden shortcut

❌ Heuristic intuition ("interaction should decay z mass scale m") without explicit analytical
derivation.

✅ Analytical derivation z asymptotic expansion ("kink tail tanh(mx/√2) decays as
e^{-m√2·x}; expected decay rate m·√2 ≈ 1.4142").

### 3.6.3 — Anti-Lakatos preservation

Analytical pre-derivation at pre-registration time = honest expectation; failure → real
structural finding. **PRE-EMPTS** post-hoc threshold modification temptation.

**Permitted relaxation:** Heuristic threshold acceptable IF marked explicitly as
"informational only" (LIT/INVENTORY class), NIE FP class.

### 3.6.4 — Retrospective application policy

NO retroactive re-evaluation of closed cycles. Closed cycles preserve pre-registered LOCKs.
Going forward, new cycles use the addendum.

**Cross-cycle annotation:** existing closed cycles MAY annotate pre-registration patterns
as "methodology lesson" (e.g., FFS-3 q=1 implicit, CE-H T_P3_2 m vs m·√2), NIE re-evaluate
verdicts.

### 3.6.5 — Phase 0 template requirement

Phase 0 balance sheet MUST include analytical pre-derivation slot (explicit subsection)
for each FP-class falsifier. Sympy symbolic computation lub LaTeX derivation acceptable.

---

## 3.6 EXTENSION (BINDING 2026-05-23) — R2 audit op-R2-audit-3-6-extension-2026-05-23

**Source:** Second R2 audit cycle (2026-05-23), triggered by R1-2 flag from γ-1 closure
(op-CE-H-3D-native-interaction-2026-05-22). 4 pattern instances revealed §3.6.1-§3.6.5
coverage insufficient dla: sign conventions, fit DoF asymmetry, implicit assumptions,
numerical precision validation.

**Pattern instances mapped:**
- Instance 1 (CE-H T_P3_2): → §3.6.9 (numerical precision validation)
- Instance 2 (FFS T_P4_3 q=1 implicit): → §3.6.8 (implicit assumption enumeration)
- Instance 3 (γ-1 T_P2_4 sign): → §3.6.6 (sign convention derivation)
- Instance 4 (γ-1 T_P3_3 DoF): → §3.6.7 (fit DoF equalization)

### 3.6.6 — Sign convention derivation (BINDING 2026-05-23)

For each pre-registered formula involving signs (slope coefficients, attraction/repulsion,
growth/decay rates), Phase 0 analytical pre-derivation MUST include explicit sign verification:

(a) **Physical principle:** State which physical principle determines sign (e.g., "like-charges
   repel in 2D Coulomb", "kink-antikink attract via Casimir-like long-range force")

(b) **Limiting case verification:** Verify sign in known limit (e.g., "at L→r_0, V_int should
   be high for like-sign because charges close" → positive in standard convention)

(c) **Convention statement:** Explicit convention chosen (e.g., "V_int = energy increase when
   bringing charges together; positive convention dla repulsion")

**Forbidden shortcut:** Assuming sign from form mathematical analogy without physical
principle derivation.

**Example (γ-1 Phase 2 lesson):**
- ❌ "V_int ~ 2π log(L/r_0) → assume positive coefficient"
- ✅ "Same-sign 2D Coulomb charges REPEL → V_int(L) high at small L, low at large L
   → coefficient -2π for log(L/r_0)"

### 3.6.7 — Fit parameter DoF equalization (BINDING 2026-05-23)

For each pre-registered FP-class falsifier involving comparison of fit forms (e.g., R²_log vs
R²_exp), comparison MUST use:

(a) **Equal parameter count:** Both fit forms with same number of free parameters
   (e.g., 2-param log vs 2-param exp, NIE 2-param log vs 3-param exp+offset)

(b) **OR adjusted R² (Akaike/BIC):** Information criteria accounting for parameter count if
   fits MUST have different DoF dla theoretical reasons

(c) **Phase 0 explicit declaration:** Number of parameters per fit form documented; jeśli
   unequal, justification + AIC/BIC adjustment specified

**Forbidden shortcut:** Comparing R² of fits with asymmetric parameter counts without
adjustment.

**Example (γ-1 Phase 3 lesson):**
- ❌ "R²_log (2-param) vs R²_exp+offset (3-param), threshold R²_log > R²_exp + 0.02"
- ✅ "R²_log (2-param: A, B) vs R²_exp (2-param: C, m), threshold R²_log > R²_exp + 0.02"
- ✅ "ΔAIC_log < ΔAIC_exp by margin 4 (significantly preferring log)"

### 3.6.8 — Implicit assumption enumeration (BINDING 2026-05-23)

For each pre-registered formula z derived value, Phase 0 analytical pre-derivation MUST
enumerate explicit:

(a) **Background assumptions:** What is assumed about field configuration vs general case
   (e.g., "spherically symmetric ansatz", "linearized about VEV")

(b) **Normalization conventions:** Field units, dimensional choices (e.g., "Phi has dimensions
   of mass", "ρ scaled by v")

(c) **Limit choices:** Which limits taken (e.g., "small-coupling λ << 1", "non-relativistic",
   "static")

(d) **Effective parameter substitutions:** Any parameter set to convenient value
   (e.g., "q=1 effective in pre-screening", "v=1 numerical normalization")

(e) **Implicit symmetries:** Symmetries imposed for tractability (e.g., "z-translation invariant")

**Phase 0 BALANCE SHEET MUST include explicit "assumptions enumeration" subsection per
FP-class falsifier.**

**Forbidden shortcut:** Using convenient pre-screening or simplified formulas without explicit
assumption documentation.

**Example (FFS T7 q=1 implicit lesson):**
- ❌ "σ = π·v² (string tension formula)" — implicit q=1
- ✅ "σ = π·q²·v² (Nielsen-Olesen z winding q); pre-screening uses q=1 effective for
   order-of-magnitude estimate; strict q=1/3 evaluated separately in full cycle"

### 3.6.9 — Numerical precision validation (BINDING 2026-05-23)

Pre-registered numerical thresholds dla FP-class falsifiers MUST be validated against
analytical computation z **explicit precision standard**:

(a) **5% accuracy standard:** Analytical pre-derivation must produce expected value w
   **±5% accuracy** of subsequent sympy execution (unless explicit looser tolerance
   justified physically)

(b) **Precision documentation:** Phase 0 documents:
   - Analytical expected value
   - Analytical derivation steps (verifying each)
   - Expected precision basis (5% default; relaxation justified)
   - Sympy verification plan

(c) **Mismatch handling:** Jeśli analytical pre-derivation gives value differing from sympy
   by > 5% (beyond expected numerical noise), Phase 0 MUST be RE-AUDITED before LOCK
   (anti-Lakatos discipline still preserved: NIE modification post-LOCK)

**Forbidden shortcut:** Pre-registering numerical thresholds based on heuristic intuition
without analytical precision validation.

**Example (CE-H T_P3_2 lesson):**
- ❌ "Expected m = 1.0 (heuristic intuition)" → fitted 1.40, 40% off
- ✅ "Expected m·√2 ≈ 1.4142 (analytical derivation z tanh tail expansion); threshold
   tolerance 10% → fitted 1.40 within 1% of analytical → PASS"

### 3.6.10 — Methodology evolution acknowledgment

§3.6.1-§3.6.5 (BINDING 2026-05-22) + §3.6.6-§3.6.9 (BINDING 2026-05-23) covered 4 pattern
instances observed across 4 cycles. Methodology evolution is **expected**: future cycles may
reveal additional pattern aspects not yet covered. Such discoveries → R1 flag → R2 audit
cycle → further §3.6 extension. **NIE infinite regress concern**: anti-Lakatos discipline
preserved (no retroactive cycle modification; closed cycles preserve LOCKs).

---

## 3.6 EXTENSION 2 (BINDING 2026-05-24) — R2 audit op-R2-audit-3-6-extension-2-2026-05-24

**Origin:** γ-3 cosmological cycle (op-CE-H-gamma-3-cosmological-2026-05-23, B+ verdict)
revealed 3 R1 flag CANDIDATES not covered by §3.6.1-§3.6.10:

- Instance 5 (γ-3 Phase 4 PARTIALs over budget): → §3.6.11 (PARTIAL taxonomy)
- Instance 6 (γ-3 Phase 5 §5 conflation): → §3.6.12 (Concept paper rigor classification)
- Instance 7 (γ-3 Phase 3 c=const audit gap, user identified): → §3.6.13 (Constants identification, HIGH priority)

### 3.6.11 — PARTIAL taxonomy + budget refinement (BINDING 2026-05-24)

**Two distinct PARTIAL categories observed across cycles must be classified separately:**

**(a) `PARTIAL_compute`** — Computation executed; result diverges od threshold by quantifiable
margin. Examples: sign convention mismatch z 99% magnitude match (γ-1 T_P2_4); R² difference
0.013 vs threshold 0.02 (γ-1 T_P3_3).

**(b) `PARTIAL_concept_mismatch`** — TGP-native framework structurally lacks direct counterpart
for observational test. Examples: Ω_m needs GR-equivalent ρ_critical (γ-3 T_P4_1); T_CMB = 2.725 K
is observational input nie TGP-derived (γ-3 T_P4_2).

**Budget:**
- `PARTIAL_compute`: max 1 per cycle (preserves cycle 1/2/7 numerical discipline)
- `PARTIAL_concept_mismatch`: NO hard limit; each must explicitly document why TGP-native
  cannot test (concept paper § / declared limit / scope statement)

**Pre-registration requirement:** Each FP musi specify outcome set:
PASS / PARTIAL_compute / PARTIAL_concept_mismatch / DEFERRED / FAIL z explicit thresholds.

### 3.6.12 — Concept paper claim rigor classification (BINDING 2026-05-24)

**Każdy concept paper claim referenced jako pre-registration foundation MUST be classified
pre-LOCK as jeden z trzech:**

**(I) `DERIVED`** — Rigorously derived from minimal axioms; reproducible w sympy/symbolic
computation. Concept paper sekcja musi cite derivation OR link to research cycle z derivation.

**(II) `STRUCTURAL_PLAUSIBLE`** — Plausible structural argument; not yet rigorously derived but
reasonable inference from axioms. Concept paper sekcja musi explicitly flag jako "do verification
w future cycle".

**(III) `QUALITATIVE`** — Intuitive/qualitative claim; NIE technical derivation. Concept paper
sekcja musi explicitly flag jako QUALITATIVE z note: "downstream cycles SHOULD NOT rely on this
claim as pre-registration without independent rigor."

**Pre-LOCK audit:** PRZED LOCKing concept paper (or any meta document) z pre-registered claims,
all referenced claims musi mieć (I/II/III) classification explicit.

**Existing LOCKED documents:** Retroactive R2 audit for classifications; identified flaws (e.g.,
TGP_GENERATED_SPACE_COSMOLOGY §5 acceleration conflation) flagged for future revision PRZED new
cycle pre-registration depends on them. **NIE retroactive verdict modification** for already-closed cycles.

### 3.6.13 — Fundamental constants identification BINDING (HIGH priority, BINDING 2026-05-24)

**Każdy cycle pre-registration MUST explicitly enumerate fundamental constants i numerical scales
used w derivations, z classification:**

**(α) `TGP_FUNDAMENTAL`** — Constant appears w minimal TGP Lagrangian/axioms as parameter (e.g.,
m_σ derived z λ, v). Justified by axiom set.

**(β) `EMERGENT_FROM_PHI`** — Constant emerges z Phi-substrate dynamics; may vary z background
Phi configuration. Examples: signal speed c (per TGP_GENERATED_SPACE_COSMOLOGY §1.1 ontology
"przestrzeń emergent z Phi"), effective metric g_μν.

**(γ) `OBSERVATIONAL_ANCHOR`** — External observational value used as input (e.g., t_universe,
T_CMB). Must be honestly declared as external.

**(δ) `APPROXIMATION_LIMIT`** — Constant approximation valid only w specific regime. Must declare
regime explicitly + flag for full derivation w future cycles.

**Pre-registration requirement:**
- Each derivation MUST list every constant used + classification
- If any constant classified (β) (emergent), cycle MUST verify whether constant actually behaves
  as constant w investigated regime, OR derive functional form explicitly
- If classified (δ) (approximation), regime of validity MUST be stated PRZED derivation

**Specific TGP-native constants requiring classification:**
- **c (signal speed):** (β) EMERGENT_FROM_PHI per concept paper §1.1 + §3.4 ontology. Cycles
  must justify c=c_0 explicitly OR derive c(Φ) functional form.
- **m_σ (sigma mode mass):** (α) TGP_FUNDAMENTAL via m_σ² = 2λv².
- **G_Newton (gravitational coupling):** NOT in current TGP Lagrangian → undefined w current
  TGP scope. Per AX-DECL declared limit.
- **t_universe (cosmological age):** (γ) OBSERVATIONAL_ANCHOR (stellar ages).
- **T_CMB (2.725 K):** (γ) OBSERVATIONAL_ANCHOR per γ-3 Phase 4 honest declaration.

**Critical:** γ-3 cycle (B+ LOCKED) used c=c_0 implicit without (β) classification. R2 audit
ujawniło ten gap. γ-3' cycle (revisit, post 2026-05-24) MUST apply §3.6.13 BINDING with c
classified (β) z explicit derivation c(Φ_frontier) PRZED any R(t) computation.

### 3.6.14 — Methodology evolution acknowledgment (extended)

§3.6.1-§3.6.5 (BINDING 2026-05-22) + §3.6.6-§3.6.10 (BINDING 2026-05-23) + §3.6.11-§3.6.13
(BINDING 2026-05-24) covered 7 pattern instances observed across 5 cycles. Methodology evolution
is **expected** + LEGITIMATE. Anti-Lakatos discipline preserved:
- No retroactive cycle modification (closed cycles preserve LOCKs)
- New sub-rules apply going FORWARD only
- R1 flag → R2 audit cycle pattern documented and proven (twice: 2026-05-23, 2026-05-24)

---

**Status:** BINDING 2026-05-04+ (§1-§5); §6 R1+R2+R3 BINDING 2026-05-22+; §3.6 BINDING
2026-05-22+ (§3.6.1-§3.6.5); §3.6 EXTENSION BINDING 2026-05-23+ (§3.6.6-§3.6.10); §3.6
EXTENSION 2 BINDING 2026-05-24+ (§3.6.11-§3.6.13). **Apply przed each new cycle's promotion claim.**
