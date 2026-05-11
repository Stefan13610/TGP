---
title: "PRE_REGISTERED_FALSIFIERS — append-only registry of pre-registered decision rules"
date: 2026-05-10
type: meta-registry
status: 🟢 ACTIVE — append-only; immutable timestamps
binding_scope: "Każdy cykl z falsifiable claim; mandatory before Phase 1 sympy"
related:
  - "[[CYCLE_KICKOFF_TEMPLATE.md]] §1, §2.3"
  - "[[VALIDATION_TRANSFERS.md]]"
  - "[[../PREDICTIONS_REGISTRY.md]]"
parent: "[[README.md]]"
tags:
  - meta
  - registry
  - pre-registration
  - falsification
  - immutable-timestamp
  - anti-moving-goalposts
---

# PRE_REGISTERED_FALSIFIERS — append-only registry

## §0 — Po co ten plik

### §0.1 — Diagnoza ryzyka (2026-05-10)

Critique od Claudian (post-mPhi-verification cascade analysis):

> "Recovery V parametric family OPEN" + "framework extension multi-session" + "specific point
> falsified, neighbourhood otwarte" — to jest classical degenerative research programme
> pattern (Lakatos): każda falsyfikacja → otwarcie nowej recovery space.

**Anti-pattern:** falsification observation → "ten konkretny point excluded, neighbourhood
otwarte" → recovery cycle → następna falsification → "ten conkretny shifted point excluded"
→ nieskończona regresja recovery spaces.

**Remediation:** **pre-registration** decision rule **PRZED** observation. Po observation
można tylko apply rule, nie redefine rule.

### §0.2 — Format pre-registration

Każdy falsifiable cycle MUSI mieć w opening commit:

```yaml
contract:
  L1_native:
    falsification_rule: "<exact decision rule>"
    pre_registration_date: <YYYY-MM-DD>      # IMMUTABLE
    pre_registration_hash: <git-SHA-of-this-commit>  # cryptographic seal
```

### §0.3 — Append-only invariant

Ten plik jest **append-only**:

- Wpisy NIGDY nie są removed lub modified
- Updates pojawiają się jako nowe wpisy `## §N+1 — Update of PR-### YYYY-MM-DD`
- Każdy update wymaga explicit reason + adversarial check
- Delete operations są forbidden (git history zachowuje audit trail)

**Hard rule:** post-observation revision rule (np. "po widzeniu 5σ falsification, redefinujemy
acceptance window") jest *forbidden* — każda revision musi być pre-registered z nową
timestampą i osobnym entry.

---

## §1 — Format entries

```markdown
### PR-<NUM>: <short cycle title>

- **Cycle:** [[../research/op-NAME/]]
- **Pre-registration date:** YYYY-MM-DD HH:MM (UTC if known)
- **Pre-registration commit:** <git SHA>
- **Native observable:** <observable in physical units>
- **Decision rule (immutable):**
  > <exact text decision rule; verbatim from kickoff commit>
- **Falsification target:** <which native coefs / framework aspect would be ruled out>
- **Confidence threshold:** <e.g., 5σ, 95% CL, ...>
- **Recovery scope (if any):** <pre-declared recovery directions, NIE post-hoc shifted points>
- **Status:** PENDING | TRIGGERED-FALSIFIED | TRIGGERED-CONFIRMED | EXPIRED
- **Result entry (post-trigger):**
  - Date observed: YYYY-MM-DD
  - Source: <observation source>
  - Outcome: <falsified | confirmed | inconclusive>
  - Reference: [[../research/op-NAME/Phase_X_results.md]]
- **Notes:** <optional>
```

---

## §2 — Initial entries (post-2026-05-10 cycles only)

> ⚠️ **Note:** Pre-2026-05-10 cycles **nie mają** pre-registration timestamps — moving-goalposts
> ryzyko unaddressed dla starszych cycles. Audit retrospective: każdy claimed-falsifiable
> result z pre-2026-05-10 cycles wymaga explicit annotation "no pre-registration; classical
> mode" w PREDICTIONS_REGISTRY.

### PR-001 (RETROACTIVE LOG): GWTC-3 RE-RUN M9.1'' falsification

- **Cycle:** [[../research/op-GWTC3-reanalysis/]]
- **Pre-registration date:** ⚠️ **NOT PRE-REGISTERED** — retrospective log only
- **Pre-registration commit:** N/A
- **Native observable:** β_ppE^TGP_(b=-1) projection na ppE chart (NOTE: not native; this is
  itself a projection-cycle, see CYCLE_KICKOFF §4 intentional projection)
- **Decision rule (retrospective):**
  > "If β_ppE^TGP prior z M9.1'' anchor falls outside GWTC-3 5σ window, M9.1'' specific
  > Taylor expansion form is excluded."
- **Falsification target:** M9.1'' specific f(ψ) = (4-3ψ)/ψ form
- **Confidence threshold:** 5σ
- **Recovery scope:** EX POST FACTO declared via emergent-metric framework (Path 1 c_0=0,
  Path 2 c_0·κ_σ=4/3) — **NOT pre-registered, declared after falsification**
- **Status:** TRIGGERED-FALSIFIED
- **Result entry:**
  - Date observed: 2026-05-09 (Phase 2 RE-RUN)
  - Source: GWTC-3 combined ~90 BBH posterior (LIGO/Virgo/KAGRA Collaboration)
  - Outcome: BF_TGP/GR = 3.5·10⁻⁶, log10 BF = -5.45, σ-level = 5.02σ FALSIFIED
  - Reference: [[../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]
- **Notes:**
  - Recovery declared *after* falsification observed — anti-pattern flag.
  - Future ET-D / CE test of recovered point should be pre-registered NOW (before
    recovery cycle proceeds).
  - Retrofit native-first version: `op-LIGO-3G-deviation` cycle in observable form
    (φ(f) function, not β_ppE parameter) — pre-registration NEXT.

### PR-002 (PROPOSED, TO-BE-FILLED): ET-D / CE recovered-point falsification

- **Cycle:** [[../research/op-LIGO-3G-deviation/]] (post-retrofit)
- **Pre-registration date:** PENDING — **MUST FILL BEFORE retrofit cycle Phase 1**
- **Pre-registration commit:** PENDING
- **Native observable:** GW inspiral phase φ(f) deviation residual przeciwko GR template,
  in radians as function of f (Hz), for BBH events with M_total > Z M_⊙ at SNR > W
- **Decision rule (PROPOSED — needs author lock):**
  > "If ET-D + CE 2027-2035 stack residual on φ(f) - φ(f)_GR < δ_threshold (TBD by Fisher
  > matrix calculation) at 5σ for any sub-window of frequency band, recovered-point
  > emergent-metric framework is structurally excluded — NOT just specific (a_n, ξ_n) point
  > but full Path 2 (σ-coupling) family."
- **Falsification target:** Emergent-metric Phase 4 Path 2 (canonical σ-coupling recovery)
- **Confidence threshold:** 5σ stack residual
- **Recovery scope:** **NONE** — pre-declared: if ET-D + CE 5σ excludes Path 2, framework
  needs deeper amendment (mechanism (v): framework extension), NOT shifted (a_n, ξ_n) point
  within same family
- **Status:** PENDING (cycle not yet retrofitted; pre-registration not yet committed)
- **Notes:** **HARD RULE for this cycle:** No recovery cycle on Path 2 falsification without
  separate explicit author authorization + new PR-### entry. Anti-Lakatos clause.

### PR-003 (PROPOSED, RECOMMENDED): TGP-native predictions time capsule

- **Cycle:** Cross-cycle (foundational meta)
- **Pre-registration date:** PENDING — **HIGH PRIORITY**
- **Pre-registration commit:** PENDING
- **Native observable:** Top-N TGP-native observable predictions in observational language
  (arcsec, ms, Hz, strain), even if exact numerical values not yet computed
- **Decision rule (PROPOSED):**
  > "Time capsule predictions sealed with cryptographic timestamp 2026-05-10. Future data
  > releases (CMB-S4, ET-D, CE, JWST cosmology, BBN refinements) compared against capsule
  > predictions. Silent revision between capsule and submission FORBIDDEN — any update
  > requires new PR-### entry with explicit reason + adversarial review."
- **Falsification target:** Anti-"we always said X" retrofit
- **Confidence threshold:** N/A — meta-protocol, not single test
- **Recovery scope:** N/A
- **Status:** PROPOSED — author authorization pending
- **Notes:** Capsule format: each prediction lists (a) observable + units, (b) TGP value
  range, (c) cycle reference, (d) measurement instrument. Sealed git tag.

---

## §3 — Anti-patterns

### §3.1 — Post-hoc rule revision

**Anti-pattern:**

```
T0: pre-registered "if β > 0.1, falsified"
T1: observation: β = 0.15
T2: revise rule: "if β > 0.2 in this specific BBH mass window, falsified"
T3: claim: "rule passed"
```

**Why bad:** Rule wasn't pre-registered with mass window restriction; restriction added after
seeing data.

**Remediation:** Hard rule: post-observation revision FORBIDDEN. Any revision = new PR-###
entry with new pre-registration timestamp + explanation why original rule was inadequate.
Original rule + result remain in registry (append-only).

### §3.2 — Underspecified decision rule

**Anti-pattern:**

```
falsification_rule: "if observation disagrees with TGP, framework is wrong"
```

**Why bad:** No specific threshold, observable, instrument, or window. Cannot trigger or fail
deterministically.

**Remediation:** Decision rule must be operationally testable: specific instrument, specific
observable, specific threshold, specific confidence level. Format: "if <instrument> measures
<observable> outside <window> at <CL>, <specific framework aspect> excluded."

### §3.3 — Unbounded recovery space

**Anti-pattern:**

```
falsification_rule: "if X exceeded, M9.1'' specific point excluded but recovery space open"
```

**Why bad:** "Recovery space open" without pre-declared bounds = degenerative research
programme. Each falsification just opens new recovery space → infinite regress.

**Remediation:** Pre-declare recovery scope in entry. Format:

```
recovery_scope:
  allowed_directions: ["σ-coupling addition with c_0·κ_σ in [3/2, 5/4]", "shift a_3 in [-1, 1]"]
  forbidden_directions: ["new free Taylor coefs beyond a_5/ξ_5", "modification of S05 axiom"]
  if_recovery_exhausted: "framework needs structural amendment (mechanism v); NOT continued
                          recovery cycles"
```

If observation falsifies and recovery_scope is exhausted, framework must enter
"structural amendment" mode (deeper change) or be acknowledged as failed.

### §3.4 — Cycle without pre-registration claiming falsifiable result

**Anti-pattern:** Cycle published as "STRUCTURAL_DERIVED falsifiable prediction" but no
PR-### entry exists.

**Why bad:** Without immutable timestamp, cycle effectively could revise rule post-observation.

**Remediation:** Hard rule: claim status `STRUCTURAL_DERIVED_NATIVE` (A-/A/A+) requires
linked PR-### entry. Without entry: max status `STRUCTURAL_VERIFIED` (C, internal consistency
only).

---

## §4 — Update protocol

When pre-registered rule needs revision (legitimate cases only):

1. **Open new PR-<NUM+1> entry** linking to original PR-<NUM>
2. **State explicit reason** (acceptable: "Phase 0 scope refinement before any data observed";
   unacceptable: "data didn't fit original rule")
3. **Adversarial review** of revision (separate agent checks revision is genuine scope change,
   not goal-post movement)
4. **Original entry preserved** — registry is append-only
5. **PREDICTIONS_REGISTRY entry updated** with reference do BOTH original i revised PR-###

**Audit trail invariant:** any future reader can reconstruct: "what was the rule at time T?"
by reading registry up to date T.

---

## §5 — Sign-off

**Doc authored:** 2026-05-10 (post-conversation autor + Claudian o pre-registration jako
anti-Lakatos clause).

**Status:** ACTIVE registry. Bootstrap entries §2 PR-001 (retroactive log) + PR-002 / PR-003
(proposed, pending author authorization).

**Insight credit:** Claudian (Lakatos diagnosis); autor (acceptance kalibracji "analytical
reduction OK, recovery without bound NOT OK").

**Mandatory next steps:**

1. Author lock decision rule for PR-002 (ET-D / CE retrofit cycle)
2. Author authorization for PR-003 (time capsule format)
3. Every new cycle post-2026-05-10 with `falsification_rule` MUST submit PR-### entry przed
   Phase 1 sympy commit
