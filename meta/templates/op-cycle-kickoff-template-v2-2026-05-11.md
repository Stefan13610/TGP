---
# ============================================================================
# MINIMAL VIABLE BOILERPLATE — copy this to research/op-<NAME>-<YYYY-MM-DD>/README.md
# ============================================================================
#
# Wszystkie pola oznaczone <<FILL>> MUSZĄ być wypełnione przed pierwszym commit.
# Pola oznaczone <<OPTIONAL>> mogą być empty list/string ale field MUSI istnieć.
#
# Po wypełnieniu, sprawdź lokalnie:
#   python tooling/validate_kickoff.py research/op-<NAME>-<YYYY-MM-DD>/README.md
#
# Validator MUSI zwrócić PASS przed Phase 0 commit.
#
# Reference: meta/CYCLE_KICKOFF_TEMPLATE.md (BINDING dla cykli post-2026-05-10)
# Validator: tooling/validate_kickoff.py
# Restart context: meta/RESEARCH_RESTART_2026-05-11.md
# ============================================================================

title: "<<FILL: krótki cycle title, ~10 słów>>"
date: <<FILL: YYYY-MM-DD>>
type: research-cycle
folder_status: parking   # default; → active po explicit user decision + WIP slot wolny
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (BINDING post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: "<<FILL: jednostki fizyczne (arcsec, Hz, ms, dimensionless strain ratio, etc.) — NIE β_ppE/γ_PPN/ξ_2>>"
    measurement_instrument: "<<FILL: Cassini | LLR | LIGO | BBN | JWST | etc.>>"
    native_coefs_constrained:
      - "<<OPTIONAL: lista native Taylor coefs z g_eff[Φ] / Φ-EOM constrained tym cyklem>>"
    falsification_rule: "<<FILL: 'jeśli <instrument> measures <observable> outside <window> at <CL>, <framework aspect> excluded'>>"
    pre_registration_date: "<<FILL: YYYY-MM-DD timestamp PRZED jakąkolwiek calculation; immutable po commit>>"

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks: []     # OPTIONAL: [Newton-limit, GR-weak-field, PPN-1PN, ppE-2PN-phase, BBN-Friedmann, ...]
    reduction_type: "not-attempted"  # MUST: analytical-exact | analytical-approximate | numerical-agreement | mapping-failed | not-attempted
    validation_transfer: ""   # OPTIONAL: jeśli analytical-exact, jakie bounds dziedziczą
    failure_disposition: "L1-stands"  # DEFAULT: failure of L2 ≠ failure of cycle

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map: []    # OPTIONAL: lista observational bounds + native coefs constrained
                              # Format: { bound: "<source>", constrains: "<coef>", window: "<value> at <CL>", status: "<pass/fail/pending>" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: T0                     # <<FILL: T0|T1|T2|T3|T4 per M9_RESTRUCTURE_NOTE §2>>
  kind: derivation              # <<FILL: derivation|audit|consistency-check|recovery|...>>
  output_type: structural       # MANDATORY post-2026-05-10: observable | projection | structural
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status: []

predecessors: []   # OPTIONAL: lista wikilinków do parent cykli
classification: SCAFFOLD  # SCAFFOLD | DERIVATION | AUDIT | RECOVERY | ...
goal: "<<FILL: 1-2 zdania głównego celu cyklu>>"
target_window: "<<FILL: zakres obserwabli + horizon czasowy>>"

six_requirements_target:
  - "P1: <<FILL>>"
  - "P2: <<FILL>>"
  - "P3: <<FILL>>"
  - "P4: <<FILL>>"
  - "P5: <<FILL>>"
  - "P6: <<FILL: S05 single-Φ axiom preserved bezwarunkowo>>"

risk_flags:
  - "R1: <<FILL>>"
  - "R2: <<FILL>>"

phase_plan:
  Phase_0: "Balance sheet + pre-flight methodology read confirmation"
  Phase_1: "<<FILL: native derivation core>>"
  Phase_2: "<<FILL>>"
  Phase_3: "<<FILL>>"
  Phase_FINAL: "Closure + optional L2 framework reduction + L3 falsification map"

tags:
  - <<FILL: relevant-tag-1>>
  - <<FILL: relevant-tag-2>>
  - cycle-scaffold-<<FILL: YYYY-MM-DD>>
---

# <<FILL: cycle title>>

> **Cel:** <<FILL: 1-2 zdania>>

## §0 — Cel + native-first contract

[CITE: `meta/CYCLE_KICKOFF_TEMPLATE.md` §1; `meta/PPN_AS_PROJECTION.md` §3.1; `meta/M9_RESTRUCTURE_NOTE.md` §2]

### §0.1 — Native observable target

<<FILL: co fizycznie liczymy, w jakich jednostkach, jaki instrument detekuje>>

### §0.2 — Pre-registered falsification rule

<<FILL: decision rule WRITTEN BEFORE any calculation; immutable timestamp>>

```
pre_registration_date: <<FILL: YYYY-MM-DD>>
pre_registration_hash: <<auto-set by git commit SHA>>
recovery_scope:
  allowed_directions: []     # explicit bounds (per PRE_REGISTERED_FALSIFIERS §3.3)
  forbidden_directions: []
  if_recovery_exhausted: "framework needs structural amendment, NOT continued recovery"
```

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

Q1-Q8 checklist per [[../../meta/CYCLE_LIFECYCLE.md]] §Phase 0 README template.

- [ ] **Q1 (Pattern coverage):** Reviewed `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §2?
- [ ] **Q2 (Red flags):** Zidentyfikowane §3 red flags?
- [ ] **Q3 (Inherited LOCKs):** Mapping w `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §4?
- [ ] **Q4 (Standard-physics tools):** Jeśli używam std tools — explicit justify czemu NOT BD-translation?
- [ ] **Q5 (m_Φ usage):** Universal vs environment-dependent (Pattern 2.5)?
- [ ] **Q6 (GR limit framing):** Distinguishes "TGP gives same number" vs "TGP IS GR"?
- [ ] **Q7 (ASK-RULE self-check):** Niewyjaśnione gaps gdzie sięgam po std physics?
- [ ] **Q8 (BD-drift audit plan):** Phase FINAL spawn `bd-drift-audit` subagent?

### §0.4 — Pre-flight methodology read confirmation

**BINDING per `meta/CYCLE_KICKOFF_TEMPLATE.md` §2.6:**

- [ ] Przeczytano [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3
- [ ] Przeczytano [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift
- [ ] Przeczytano [[../../meta/M9_RESTRUCTURE_NOTE.md]] §1.4 + §3 — M9.1'' jako anchor, NIE framework
- [ ] Przeczytano [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 — kickoff contract

**Sign-off:** <<FILL: agent name | author>> @ <<FILL: timestamp>>

### §0.5 — Sympy substance plan (NEW: post-2026-05-11 audit lesson)

Per [[../../meta/AUDIT_2026-05-11_sympy_substance.md]] §4: 24/104 testów w 2026-05-11 cohort
to literal `T_pass = True`. To jest anti-pattern. Plan sympy dla tego cyklu:

- [ ] Każdy test sympy ma **explicit pytanie fizyczne** które weryfikuje
- [ ] **Co najmniej 50% testów** to non-trivial symbolic manipulation (NIE `T_pass = True`)
- [ ] **Co najmniej 1 test** wykonuje first-principles derivation z TGP axioms (S05 / single-Φ / Φ-EOM / substrate-vacuum identification)
- [ ] Structural declarations (S05 preservation, scope documentation) są raportowane **osobno** od sympy tests, NIE counted jako "8/8 sympy PASS"

Plan:
- Phase 1 sympy first-principles derivations: <<FILL: lista konkretnych derivations>>
- Phase 1 literature-anchored consistency checks: <<FILL: lista>>
- Structural declarations: <<FILL: lista>> (NIE counted w sympy PASS total)

## §1 — Phase 0: scope mapping + balance sheet

[<<FILL: balance sheet content>>]

## §2 — Phase 1: native derivation

[<<FILL: po user authorization "active">>]

## §FINAL — Optional L2 framework reduction

[OPTIONAL — last stage; failure here does NOT invalidate Phase 1-N native results]

## §FINAL+1 — L3 falsification map check

[<<FILL>>]

---

## Status

🟡 **PARKING — scaffold opened <<FILL: YYYY-MM-DD>>**. Pre-flight methodology read
confirmation: <<FILL: complete | pending>>. Validator status: <<FILL: PASS | FAIL>>.

Phase 0 commit gate:
1. `python tooling/validate_kickoff.py research/<this-dir>/README.md` → must PASS
2. PR-### entry w `meta/PRE_REGISTERED_FALSIFIERS.md` z immutable timestamp (jeśli falsifiable claim)
3. User authorization "active" + WIP slot wolny

Aż wszystkie 3 gate'y PASS — cycle pozostaje w `parking`. **Bez wypełnienia BINDING contract::,
cycle NIE może aspirować do statusu wyższego niż `PROJECTION_VERIFIED`** per CYCLE_KICKOFF_TEMPLATE §0.2.

---

**Cycle scaffolded:** <<FILL: YYYY-MM-DD>> (<<FILL: author/agent>>, restart per
`meta/RESEARCH_RESTART_2026-05-11.md` clean kickoff schema).

**Cross-references:**
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] §1-§2 (BINDING contract)
- [[../../meta/CYCLE_LIFECYCLE.md]] §Claim status taxonomy
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] §3.3 (recovery_scope), §3.4 (PR-### entry)
- [[../../meta/AUDIT_2026-05-11_sympy_substance.md]] §4 (sympy substance lessons)
- [[../../tooling/validate_kickoff.py]] (technical enforcement gate)
- [[../../meta/RESEARCH_RESTART_2026-05-11.md]] (operational restart context)
