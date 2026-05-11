---
title: "CYCLE_KICKOFF_TEMPLATE — mandatory contract dla nowych cykli (post-2026-05-10)"
date: 2026-05-10
type: meta-policy
status: 🟢 ACTIVE — BINDING dla wszystkich cykli otwieranych post-2026-05-10
binding_scope: "Każdy nowy `research/op-*` cykl, szczególnie gravity/cosmology sectors"
related:
  - "[[CYCLE_LIFECYCLE.md]]"
  - "[[PPN_AS_PROJECTION.md]]"
  - "[[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]]"
  - "[[M9_RESTRUCTURE_NOTE.md]]"
  - "[[CALIBRATION_PROTOCOL.md]]"
  - "[[VALIDATION_TRANSFERS.md]]"
  - "[[PRE_REGISTERED_FALSIFIERS.md]]"
parent: "[[README.md]]"
tags:
  - meta
  - kickoff
  - native-first
  - anti-drift
  - mandatory-contract
---

# CYCLE_KICKOFF_TEMPLATE — mandatory contract dla nowych cykli

## §0 — Po co ten dokument

### §0.1 — Diagnoza weekendowa autora (2026-05-10)

> "Agenci pracowali autonomicznie w PPN/ppE-projection mode bo (a) literatura beyond-GR jest
> przytłaczająco w tym języku, (b) brak explicit kontraktu startowego cyklu o native-first
> → agenci defaultowo szukają 'compatibility layer' z istniejącymi narzędziami. Naiwne
> założenie 'agenci nie dryfują' zawiodło."

**Przyczyna driftu:** Brak kontraktu w opening commit cyklu pozwalał agentom inferować scope
z literatury (która jest w PPN/ppE basis). Methodology meta-docs istniały, ale nie były
*egzekwowane* na poziomie kickoff.

### §0.2 — Co ten template robi

Wymusza explicit contract w opening commit każdego nowego cyklu. Bez wypełnienia wszystkich
sekcji L1/L2/L3 kontraktu, cykl **nie może aspirować** do statusu wyższego niż
`PROJECTION_VERIFIED` (najsłabszy claim level).

### §0.3 — Filozofia: native-first, mapping last (opcjonalny)

**Kalibracja autora 2026-05-10 (po dyskusji z Claudian):**

> "Czasem najlepszym dowodem jest jednak analityka, bo matematycznie nie da się pokryć
> wszystkiego. Analityczne wyprowadzenie rdzeniowych wzorów TGP czasem powinno dać się
> zmapować na klasyczne frameworki, tylko to nie powinien być główny cel. Raczej ostatni
> opcjonalny etap cyklu którego fail nie determinuje błędu."

**Trzy tryby pracy** (rozróżnienie krytyczne):

| Tryb | Status | Action |
|---|---|---|
| **Mimicry mode** (drift) | Output = parameter w obcym frameworku jako *primary* (np. "β_ppE = -15/4 jest predykcją TGP") | ❌ REJECTED — nie pass kickoff |
| **Native-with-mapping** | Output = native equations + observables; ostatnia faza: analytical reduction do GR/Newton/PPN/ppE potwierdzona | ✅ NAJSILNIEJSZY — validation transfer aktywny |
| **Native-only** | Output = native equations + observables; mapping nie próbowany lub failed | ✅ VALID — brak transferu, ale brak błędu |

**Hard rule:** L1 musi być wypełnione i kompletne **zanim** L2 jest próbowane. Nie wolno
zacząć od L2.

---

## §1 — Mandatory contract template

**Skopiuj poniższe do `research/op-<NAME>-<YYYY-MM-DD>/README.md` przed Phase 0:**

```yaml
---
title: "<short cycle title>"
date: <YYYY-MM-DD>
type: research-cycle
folder_status: parking   # default; → active po explicit decyzji + WIP slot
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (mandatory post-2026-05-10) ==============
contract:
  # --- L1: Native (MANDATORY) ---
  L1_native:
    output_observable: ""           # MUST: jednostki fizyczne (arcsec, Hz, ms, dimensionless strain ratio, etc.)
                                    # NIE: bezwymiarowe parametry phenomenological frameworków (β_ppE, γ_PPN, ξ_2)
    measurement_instrument: ""      # Jaki realny pomiar dotyka tej obserwabli? (Cassini, LLR, LIGO, BBN, JWST, ...)
    native_coefs_constrained: []    # Lista native Taylor coefs z g_eff[Φ] / Φ-EOM constrained tym cyklem
    falsification_rule: ""          # Pre-registered decision rule: "jeśli observed X > Y at confidence Z, native coefs constrained do W"
    pre_registration_date: ""       # YYYY-MM-DD timestamp PRZED jakąkolwiek calculation; immutable

  # --- L2: Cross-framework reduction (OPTIONAL — last stage) ---
  L2_framework_reduction:
    target_frameworks: []           # [Newton-limit, GR-weak-field-Schwarzschild, PPN-1PN, ppE-2PN-phase, BBN-Friedmann, ...]
    reduction_type: ""              # analytical-exact | analytical-approximate | numerical-agreement | mapping-failed | not-attempted
    validation_transfer: ""         # Jeśli analytical-exact: jakie observational bounds dziedziczą do TGP w tym regime
    failure_disposition: "L1-stands"  # DEFAULT: failure of L2 mapping ≠ failure of cycle

  # --- L3: Falsification map (consistency) ---
  L3_falsification_map: []          # Lista observational bounds (Cassini, LLR, GWTC, ...) i które native coefs constrain
                                    # Format: { bound: "<source>", constrains: "<coef>", window: "<value> at <CL>", status: "<pass/fail/pending>" }

# ============== END KICKOFF CONTRACT ==============

tgp_status:
  level: <T0|T1|T2|T3|T4>           # per M9_RESTRUCTURE_NOTE §2 tier hierarchy
  kind: <derivation|audit|consistency-check|recovery|...>
  output_type: <observable|projection|structural>   # MANDATORY post-2026-05-10
  core_compatibility: <current|review-only>
  may_edit_core: <true|false>
  has_needs_file: <true|false>
  has_findings_file: <true|false>
  exports_findings: <true|false>
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status: []
---

# <cycle title>

## §0 — Cel + native-first contract

[CITE: meta/CYCLE_KICKOFF_TEMPLATE.md §1; meta/PPN_AS_PROJECTION.md §3.1; meta/M9_RESTRUCTURE_NOTE.md §2]

### §0.1 — Native observable target

<co fizycznie liczymy, w jakich jednostkach, jaki instrument detekuje>

### §0.2 — Pre-registered falsification rule

<decision rule WRITTEN BEFORE any calculation; immutable timestamp>

### §0.3 — TGP-native check (mandatory, pre-Phase-1)

[wymagane checklist Q1-Q8 per CYCLE_LIFECYCLE.md §"Phase 0 README template"]

## §1 — Phase 0: scope mapping
...

## §2 — Phase 1: native derivation
...

## §3 — Phase 2-N: native obserwable computation
...

## §FINAL — Optional L2 framework reduction

[OPTIONAL — last stage; failure here does NOT invalidate Phase 1-N native results]

## §FINAL+1 — L3 falsification map check
...
```

---

## §2 — Reguły egzekwujące

### §2.1 — Hard rule: ordering enforcement

**L1 musi być wypełnione **kompletnie** zanim L2 jest próbowane.** Konkretnie:

1. **Phase 1-(N-1)** to native derivation z observable output. Wszystkie sympy verification
   na poziomie L1.
2. **Phase N (penultimate)** może próbować L2 framework reduction — ale **tylko jeśli L1
   results są stabilne i sympy-PASS**.
3. **Phase FINAL** = closure z L3 falsification map check + optional L2 mapping summary.

**Anti-pattern:** "zacząć od policzenia β_ppE i wstecznie zrekonstruować native coefs". To
jest mimicry mode (drift). Auto-reject.

### §2.2 — Hard rule: output_type: observable for highest claim status

Status mapping (extension CYCLE_LIFECYCLE):

| `output_type` | Max claim status |
|---|---|
| `observable` (jednostki fizyczne) | **A+** (jeśli L2 transfer aktywny), **A** (L2 attempted-failed), **A-** (L2 not-attempted) |
| `structural` (algebra konsystentna, brak observable) | **C** |
| `projection` (parameter w bazie phenomenological framework jako primary) | **B** (drift-deprecated) |

**Hard rule:** Brak `output_observable` → `output_type: structural` lub `projection`. Cykl
NIE może aspirować do A+/A/A-.

### §2.3 — Soft rule: pre-registration timestamp

`pre_registration_date` musi być w opening commit cyklu. Po Phase 1+ rozpoczęciu, zmiana
`falsification_rule` wymaga explicit annotation `## §X — Falsification rule revision`
z uzasadnieniem (akceptowalne uzasadnienia: scope refinement post-Phase-0; nieakceptowalne:
"wynik nie pasował do oryginalnej reguły, zmieniam regułę").

### §2.4 — Soft rule: failure disposition declaration

`L2_framework_reduction.failure_disposition` musi być jednym z:

- `L1-stands` (DEFAULT) — failure of L2 mapping doesn't invalidate native L1 result
- `L1-and-L2-required` (RZADKO) — gdy cykl explicit potrzebuje L2 success dla claim
  (np. "cykl jest ABOUT pokazaniem GR limit"); musi być uzasadnione w §0.1
- `not-applicable` — cykl nie próbuje L2 mapping (cosmology pure native)

**Hard rule:** `L1-and-L2-required` requires explicit author authorization w §0.1.

### §2.5 — Hard rule: forced-zero parameters: declare, don't sympy-verify

Per [[PPN_AS_PROJECTION.md §3.2]]:

- α_i ≡ 0 (preferred-frame breaking) — strukturalnie z Lorentz-invariance Φ-substratu
- ζ_i ≡ 0 (energy/momentum non-conservation) — strukturalnie z covariant Φ-EOM

Cykl który "weryfikuje α_i = 0 sympy 5/5 PASS" robi nic merytorycznego — to tożsamość. Format:
jedna linijka deklaracji + cite S05 + §5.1 FOUNDATIONS.

### §2.6 — Cross-cycle pre-flight check

Przed Phase 1 sympy, agent MUST przeczytać:

1. [[meta/PPN_AS_PROJECTION.md]] §3.1 (three-layer L1/L2/L3)
2. [[meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 (anti-BD-drift)
3. [[meta/M9_RESTRUCTURE_NOTE.md]] §1.4 + §3 (M9.1'' inheritance reading rules)
4. [[meta/CYCLE_KICKOFF_TEMPLATE.md]] (this doc) §1-§2

**Hard rule:** kickoff commit **MUST** zawierać confirmation block:

```markdown
## §0.4 — Pre-flight methodology read confirmation

- [ ] Przeczytano PPN_AS_PROJECTION §3.1 — three-layer L1/L2/L3
- [ ] Przeczytano TGP_NATIVE_COMPUTATIONAL_PATTERNS §1-§4 — anti-BD-drift
- [ ] Przeczytano M9_RESTRUCTURE_NOTE §1.4 + §3 — M9.1'' jako anchor, nie framework
- [ ] Przeczytano CYCLE_KICKOFF_TEMPLATE §1-§2 (this doc) — kickoff contract

**Sign-off:** <agent name> @ <timestamp>
```

---

## §3 — Adversarial pre-flight check

Przed otwarciem nowego cyklu (przed pierwszym Phase 0 commit), opcjonalnie:

```bash
# Adversarial agent reviews kickoff for drift markers:
agent dispatch --role adversarial \
  --input research/op-NAME/README.md \
  --check-for "drift markers per CYCLE_KICKOFF §3.1"
```

**Drift markers (auto-reject z konkretnym powodem):**

| Marker | Powód | Remediation |
|---|---|---|
| `output_observable: β_ppE` lub podobne | Bezwymiarowy projection param jako primary | Zamień na native observable (φ(f), strain, arcsec, etc.) |
| `output_observable: ""` (puste) | Brak observable target | Wypełnij L1 lub mark `output_type: structural` |
| Phase 1 plan zaczyna od "policzymy β_X" | Mimicry ordering | Reorganizuj: Phase 1 = native; Phase FINAL = optional L2 |
| Brak `pre_registration_date` | Bez immutable timestamp moving-goalposts ryzyko | Add timestamp przed Phase 1 rozpoczęciem |
| `failure_disposition: L1-and-L2-required` bez §0.1 justification | Niejawne wymuszanie mapping success | Justify lub zmień na L1-stands |
| Brak `## §0.4 Pre-flight methodology read confirmation` | Bypass methodology read | Wymagana sekcja przed Phase 1 |

**Adversarial pre-flight result:**

- **PASS** → cykl może otrzymać `folder_status: active` przy WIP slot wolnym
- **PASS-WITH-FLAGS** → fix flags first, potem `active`
- **REJECT** → reorganizuj kickoff lub świadoma decyzja (akceptowalne też: "mapping cycle,
  nie native cycle" — wtedy explicit `output_type: projection`, max status B)

---

## §4 — Wzorzec dla cykli "intentional-projection"

Niektóre cykle są **intencjonalnie** projection-only — np. mapping TGP na ppE basis dla
external compatibility. To nie jest drift — to jest świadomy translation cycle. Format:

```yaml
contract:
  L1_native:
    output_observable: "N/A — intentional projection cycle"
    measurement_instrument: "N/A"
    native_coefs_constrained: []
    falsification_rule: "N/A — translation-only, no falsifiable claim"
    pre_registration_date: "<YYYY-MM-DD>"

  L2_framework_reduction:
    target_frameworks: ["ppE", "PPN", "BD"]
    reduction_type: "intentional-projection-this-IS-the-output"
    validation_transfer: "N/A — output is the projection, not native physics"
    failure_disposition: "L1-not-applicable"

  L3_falsification_map: []  # NIE robi falsification claims; tylko translation table

tgp_status:
  output_type: projection      # explicit, intentional
  kind: framework-translation  # explicit kind
```

**Hard rule:** intentional projection cycle **NIE może** aspirować do statusu wyższego niż
`PROJECTION_VERIFIED`. Może być cytowane jako *consistency check* dla native cycles, ALE
NIE może być cytowane jako *falsifiable prediction*.

**Przykład legitymny:** `op-GWTC3-reanalysis` — musiał użyć β-prior bo LIGO data są tak
filtrowane. To jest necessary translation, nie drift.

**Przykład illegitymny (drift):** `op-ppE-mapping` Phase 1 — produkował β_ppE jako "TGP
prediction", brak L1 native. Powinien być `output_type: projection`, max status B.

---

## §5 — Migration policy dla starszych cykli

**Pre-2026-05-10 cykle:** ten template **nie jest** wymagany retroaktywnie. ALE:

- Inheriting LOCKs z starszych cykli wymaga §0.4 confirmation per [[CYCLE_LIFECYCLE.md]] §
  TGP-native check (post-2026-05-10 binding)
- Cykle wybrane do retrofit (per Phase 5 retrofit case study) są przepisywane w nowym formacie
- Audit cycle `op-projection-cycles-triage-2026-05-10` (post-Phase-0 plan) klasyfikuje stare
  cykle według §2.2 status mapping

---

## §6 — Sign-off

**Doc authored:** 2026-05-10 (post-conversation autor + Claudian o "ok działaj" command).

**Status:** ACTIVE binding dla cykli post-2026-05-10. Pre-2026-05-10 cykle handled przez
migration policy §5.

**Insight credit:** autor TGP weekend audit 2026-05-10 (drift diagnosis); Claudian (template
formalization).

**Cross-references:**

- [[CYCLE_LIFECYCLE.md]] — claim status taxonomy (extended this doc)
- [[PPN_AS_PROJECTION.md]] — three-layer L1/L2/L3 methodology source
- [[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] — anti-BD-drift methodology source
- [[M9_RESTRUCTURE_NOTE.md]] — M9.1'' inheritance reading rules
- [[CALIBRATION_PROTOCOL.md]] — adversarial verification protocol
- [[VALIDATION_TRANSFERS.md]] — registry where L2-success transfers happen
- [[PRE_REGISTERED_FALSIFIERS.md]] — append-only registry of pre-registered decision rules
