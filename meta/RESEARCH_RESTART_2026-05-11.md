---
title: "RESEARCH_RESTART_2026-05-11 — operational restart na nowym czystym schemacie (post-external-review)"
date: 2026-05-11
type: meta-operational
status: 🔴 BINDING — wszystkie nowe cykle post-2026-05-11 MUSZĄ stosować się do tego protokołu
parent: "[[README.md]]"
authorization: "autor projektu, conversation 2026-05-11, options A+B+F+K+L (Rec 1+3+F+2+4)"
related:
  - "[[CYCLE_KICKOFF_TEMPLATE.md]] (BINDING contract spec)"
  - "[[CYCLE_LIFECYCLE.md]] §Claim status taxonomy"
  - "[[PRE_REGISTERED_FALSIFIERS.md]] §3.3 (recovery_scope), §3.4 (PR-### entry)"
  - "[[AUDIT_2026-05-11_sympy_substance.md]] (empirical baseline z 112 testów)"
  - "[[../tooling/validate_kickoff.py]] (technical enforcement gate)"
  - "[[templates/op-cycle-kickoff-template-v2-2026-05-11.md]] (minimal viable boilerplate)"
tags:
  - meta
  - operational
  - restart
  - post-audit
  - halt-mechanism
  - enforcement
  - clean-schema
---

# RESEARCH_RESTART_2026-05-11 — operational restart na nowym czystym schemacie

## §0 — Co to za plik

### §0.1 — Kontekst (skondensowany)

Zewnętrzna recenzja autora projektu 2026-05-11 zidentyfikowała procedural + substantive
drift w cohort 2026-05-11 (N1+N2+N3+N4+N5+cluster+hierarchy):

- **0/7 cykli** miało BINDING `contract::` blok per `meta/CYCLE_KICKOFF_TEMPLATE.md`
- **0/112 testów sympy** wykonywało first-principles derivation z TGP axioms
- **24/104 testów** to literal `T_pass = True` (algebraic mimicry pattern)
- **Cluster cycle** miał Lakatos OR-clause verdict-logic (H1a → H1b bez pre-bounded recovery_scope)

Wykonane (Rec 1-3, options A+B+F+K):
- ✅ **Rec 1** — reklasyfikacja 6 cykli z STRUCTURAL_DERIVED → STRUCTURAL_VERIFIED (C)
- ✅ **Rec 3** — adversarial audit 112 testów (decydowalne dane, append-only audit record)
- ✅ **Rec 3+F** — differential downgrade: N1, N3 → D (ALGEBRAIC_MIMICRY); N2/N4/N5/Cluster preserve C
- ✅ **Rec 2** — cluster cycle → EARLY_HALT_HONEST (`closed-NULL`)
- 🔄 **Rec 4** (TEN PLIK) — halt + enforcement mechanism + clean restart schema

### §0.2 — Cel tego pliku

Operational instructions dla agentów i autora **jak zacząć nowy cykl po restart**:
1. Halt mechanism (jakie cykle są obecnie zamrożone)
2. Technical enforcement gate (validator + jak go używać)
3. Clean kickoff workflow (krok-po-kroku od zera do active cycle)
4. Anti-drift checklist (co NIE robić)

### §0.3 — Co restart oznacza, a co NIE oznacza

**RESTART oznacza:**
- Nowe cykle MUSZĄ stosować nowy BINDING workflow
- Open scaffolds bez BINDING contract:: są zamrożone do parking-pending-new-kickoff
- Każdy nowy cycle requires validator PASS przed Phase 0 commit

**RESTART NIE oznacza:**
- Usuwanie istniejących cykli (audit trail invariant — append-only)
- Inwalidacja realnych findings z poprzednich cykli (np. L01 ADDENDUM α/(3π) typo correction; ROFM tools; emergent-metric Phase 1 foundation)
- Restart całej TGP framework (gravity sektor Tier 1 {A,B,C} foundation per M9_RESTRUCTURE §2 confirmed clean L1-native per #6 audit row)

---

## §1 — Halt mechanism (aktualnie zamrożone cykle)

### §1.1 — Scaffoldy 2026-05-11 status update

| Cycle | Pre-restart status | **Post-restart status** | Reason |
|---|---|---|---|
| `op-S07-reset-alternative-f-psi-2026-05-11` | OPEN scaffold (parking) | **parking-pending-new-kickoff** | Brak `contract::` blok; halt do new kickoff z BINDING template |
| `op-inflation-substrate-genesis-2026-05-11` | OPEN scaffold (parking) | **parking-pending-new-kickoff** | Brak `contract::` blok; halt do new kickoff z BINDING template |
| `op-LIGO-3G-native-phase-residual-2026-05-11` | OPEN parking (z `contract::`) | **parking-validator-PASS** | Jedyny cycle 2026-05-11 z BINDING template; validator PASS; ready dla user activation pending WIP slot |

### §1.2 — Re-activation procedure dla halted scaffolds

Aby przywrócić halted scaffold do active:

1. **Rewrite README.md** używając `meta/templates/op-cycle-kickoff-template-v2-2026-05-11.md` boilerplate
2. **Wypełnij wszystkie `<<FILL>>` placeholders** (output_observable, falsification_rule, pre_registration_date, etc.)
3. **Run validator:** `python tooling/validate_kickoff.py research/<scaffold-dir>/README.md`
4. **Validator MUSI zwrócić PASS** przed Phase 0 commit
5. **Submit PR-### entry** w `meta/PRE_REGISTERED_FALSIFIERS.md` jeśli falsifiable claim
6. **User authorization** "active" + WIP slot wolny per `meta/CYCLE_LIFECYCLE.md` WIP-limit

### §1.3 — Cykle z 2026-05-11 cohort (post-Rec-1+2+3+F status)

NIE są re-activated bez explicit user decision:

| Cycle | Final claim_status | Retrofit candidate |
|---|---|---|
| N1 EM | D (ALGEBRAIC_MIMICRY) | `op-L01-N1-retrofit-native` ~3-5 sesji |
| N2 QCD | C (LITERATURE_ANCHORED) | `op-L01-N2-retrofit-native` ~4-6 sesji |
| N3 SPARC | D (ALGEBRAIC_MIMICRY, 5/8 hardcoded) | `op-L01-N3-retrofit-native-SPARC` ~2-3 sesji |
| N4 Higgs | C (MIXED, Phase 1 substantive) | `op-L01-N4-retrofit-native-Higgs` ~5-8 sesji |
| N5 EW gauge | C (LITERATURE_ANCHORED) | `op-L01-N5-retrofit-native-EW` ~4-6 sesji |
| Cluster | EARLY_HALT_HONEST (`closed-NULL`) | `op-cluster-sterile-nu-prediction-2026-XX` (separate cycle z pre-bounded recovery_scope) |
| Hierarchy | STRUCTURAL_NO_GO (honest, preserved) | `op-composite-Higgs-substrate-TGP` (deferred sub-case) |

**Retrofit cycle path:** każdy retrofit jest **nowym cyklem** z BINDING template, NIE
in-place modification poprzedniego cyklu (audit trail invariant). Poprzedni cycle
pozostaje `closed-resolved` lub `closed-NULL` jako historical record.

---

## §2 — Technical enforcement gate

### §2.1 — validate_kickoff.py — opis

`tooling/validate_kickoff.py` jest **pure-stdlib** Python script (no PyYAML dependency)
który sprawdza README.md cyklu przeciwko BINDING `meta/CYCLE_KICKOFF_TEMPLATE.md` §1+§2:

**Sprawdzane pola:**
- `contract::` blok present (top-level YAML)
- `contract.L1_native.output_observable` non-empty
- `contract.L1_native.measurement_instrument` non-empty
- `contract.L1_native.falsification_rule` non-empty
- `contract.L1_native.pre_registration_date` non-empty (immutable timestamp)
- `contract.L2_framework_reduction.reduction_type` non-empty
- `output_type` field present (anywhere — top-level lub tgp_status nested)
- BODY: `## §0.4 — Pre-flight methodology read confirmation` section present
- BODY: §0.4 sekcja zawiera references do wszystkich 4 methodology files
  (PPN_AS_PROJECTION, TGP_NATIVE_COMPUTATIONAL_PATTERNS, M9_RESTRUCTURE_NOTE,
  CYCLE_KICKOFF_TEMPLATE)

**Cutoff:** validator SKIPuje cykle z `date:` < 2026-05-10 (pre-BINDING).

### §2.2 — Usage modes

```bash
# Single file (use case: dev'ing new cycle)
python tooling/validate_kickoff.py research/op-NAME-YYYY-MM-DD/README.md

# All post-cutoff cycles (use case: full audit; restart sanity check)
python tooling/validate_kickoff.py --all-new

# Pre-commit hook (use case: automated gate)
python tooling/validate_kickoff.py --pre-commit  # reads paths from stdin

# Quiet mode (only FAIL details + summary)
python tooling/validate_kickoff.py --all-new --quiet
```

**Exit codes:** 0 = PASS, 1 = FAIL, 2 = FILE_ERROR.

### §2.3 — Baseline test results (2026-05-11)

Wykonano `python tooling/validate_kickoff.py --all-new` po restart implementation:

- **Total cycles post-cutoff:** 18
- **PASS:** 1 — `op-LIGO-3G-native-phase-residual-2026-05-11` (jedyny ze BINDING template)
- **FAIL:** 17 — wszystkie pozostałe cykle 2026-05-11 cohort + scaffolds

To **empirycznie potwierdza recenzję zewnętrzną** (16/17 missing 3 fields:
`contract::`, `output_type`, `§0.4`).

### §2.4 — Git pre-commit hook (optional, recommended)

**Instalacja pre-commit hook (POSIX):**

```bash
cat > .git/hooks/pre-commit << 'EOF'
#!/bin/sh
# Pre-commit hook: validate kickoff contract for new/modified op-* README.md files
git diff --cached --name-only --diff-filter=ACMR | \
    grep -E 'research/op-[^/]+-[0-9]{4}-[0-9]{2}-[0-9]{2}/README\.md$' | \
    python TGP/TGP_v1/tooling/validate_kickoff.py --pre-commit
EOF
chmod +x .git/hooks/pre-commit
```

**Windows (PowerShell):** użyj `git config core.hooksPath` z osobnym shell script (cygwin/git-bash).

Hook blokuje commit jeśli jakikolwiek staged `research/op-*/README.md` nie passuje
validatora. Lokalne `--no-verify` bypass disabled by policy — patrz CLAUDE.md path rules.

---

## §3 — Clean kickoff workflow (krok-po-kroku)

### §3.1 — Step 1: scaffold directory

```bash
mkdir research/op-<DESCRIPTIVE-NAME>-$(date +%Y-%m-%d)
cd research/op-<DESCRIPTIVE-NAME>-$(date +%Y-%m-%d)
```

Naming conventions:
- `<DESCRIPTIVE-NAME>` = hyphenated short title (e.g., `cluster-sterile-nu-prediction`)
- `YYYY-MM-DD` = today (used by validator for cutoff check)

### §3.2 — Step 2: copy template

```bash
cp meta/templates/op-cycle-kickoff-template-v2-2026-05-11.md README.md
```

### §3.3 — Step 3: fill placeholders

Open `README.md`, replace wszystkie `<<FILL: ...>>` z konkretnymi wartościami:

**Mandatory (BINDING):**
- `title`, `date`
- `contract.L1_native.output_observable` (fizyczne jednostki!)
- `contract.L1_native.measurement_instrument`
- `contract.L1_native.falsification_rule`
- `contract.L1_native.pre_registration_date`
- `contract.L2_framework_reduction.reduction_type`
- `tgp_status.output_type` (observable | projection | structural)
- §0.4 Pre-flight methodology read confirmation checklist
- §0.5 Sympy substance plan

**Optional (recommended):**
- `predecessors`, `goal`, `target_window`, `six_requirements_target`, `risk_flags`
- §0.1, §0.2, §0.3 (Q1-Q8 checklist)
- §0.5 sympy substance plan (first-principles vs literature-anchored split)

### §3.4 — Step 4: run validator

```bash
python tooling/validate_kickoff.py research/op-<NAME>-<DATE>/README.md
```

**MUST return PASS** przed Phase 0 commit.

Typical FAIL output (do diagnozy):
```
FAIL — 3 missing field(s):
    - YAML.contract: (BINDING block per CYCLE_KICKOFF_TEMPLATE §1)
    - YAML.contract.L1_native.pre_registration_date
    - BODY: ## §0.4 — Pre-flight methodology read confirmation (BINDING §2.6)
```

Działaj na każdy FAIL until PASS.

### §3.5 — Step 5: PR-### entry (jeśli falsifiable claim)

Edit `meta/PRE_REGISTERED_FALSIFIERS.md`, dodaj entry:

```markdown
### PR-<NEXT-NUM>: <cycle title>

- **Cycle:** [[../research/op-<NAME>-<DATE>/]]
- **Pre-registration date:** <YYYY-MM-DD HH:MM>
- **Pre-registration commit:** <git SHA — filled after commit>
- **Native observable:** <z contract.L1_native.output_observable>
- **Decision rule (immutable):**
  > <verbatim z contract.L1_native.falsification_rule>
- **Falsification target:** <which native coefs / framework aspect>
- **Confidence threshold:** <5σ, 95% CL, ...>
- **Recovery scope (if any):**
    allowed_directions: [...]
    forbidden_directions: [...]
- **Status:** PENDING
```

### §3.6 — Step 6: user authorization "active"

Cycle pozostaje `folder_status: parking` aż:

1. Validator PASS
2. PR-### entry committed
3. User authorization explicit "active" + WIP slot wolny per CYCLE_LIFECYCLE WIP-limit (max 5)

Zmień YAML `folder_status: parking` → `folder_status: active` w osobnym commit z message:
```
op-<NAME>: parking → active (user-authorized, WIP slot X of 5)
```

### §3.7 — Step 7: Phase 0 balance sheet

Standard Phase 0 plan per CYCLE_LIFECYCLE §Phase 0 README template (Q1-Q8 checklist).

### §3.8 — Step 8: Phase 1 sympy z first-principles requirement

**NEW post-2026-05-11 audit lesson:** §0.5 sympy substance plan jest BINDING.

Minimum requirements dla Phase 1 sympy:
- Każdy test ma explicit pytanie fizyczne
- **≥50% testów = non-trivial symbolic manipulation** (NIE `T_pass = True`)
- **≥1 test = first-principles derivation** z TGP axioms (S05 / single-Φ / Φ-EOM /
  substrate-vacuum identification) — NIE substytucja literatury

Structural declarations (S05 preservation, scope documentation) raportowane **osobno**
od sympy tests w Phase results, NIE counted jako "8/8 sympy PASS".

---

## §4 — Anti-drift checklist (czego NIE robić)

Po wykonaniu Rec 1-4, lista anti-patterns które są **explicitly forbidden** w nowym
restart workflow:

### §4.1 — Anti-pattern: skip kickoff template

**Forbidden:**
- Otwieranie nowego cyklu bez `contract::` blok
- Otwieranie cyklu z `pre_registration_date: <pusty string>` lub `: TBD`
- Phase 1 sympy commit przed validator PASS

**Why bad:** powtórzenie 2026-05-11 cohort drift pattern (0/7 BINDING compliance).

### §4.2 — Anti-pattern: hardcoded `T_pass = True` jako sympy verification

**Forbidden:**
```python
T5_pass = True  # dimensional analysis is structural, no symbolic test
T5_pass = True  # honest documentation is structural
T_disjoint = True  # operator class enumeration is by construction
```

**Why bad:** to NIE jest sympy verification. To jest structural declaration.

**Remediation:**
- Either wykonaj rzeczywistą sympy weryfikację (np. sprawdź że symbolic difference == 0
  po explicit derivation)
- Or raportuj jako "structural declaration" osobno od sympy PASS count
- Or oznacz test jako `T5_status = "DECLARATIVE"` (NIE `T_pass = True`) i licz osobno

Per `meta/AUDIT_2026-05-11_sympy_substance.md` §4.3: 24/104 testów cohort = literal True.
Limit dla nowych cykli: **<25% testów może być `_pass = True`** (recommend: <10%).

### §4.3 — Anti-pattern: literature value cited as "TGP derivation"

**Forbidden:**
```python
# β-function from CDH 1974
beta_one_loop = sp.Rational(2, 3) * N_f * alpha**2 / pi
# ... potem ...
T1_pass = (beta - alpha**2 * 2/(3*pi) == 0)  # SAME expression substituted from itself
```

**Why bad:** weryfikujesz że wpisałeś poprawnie literacką stałą. To jest LITERATURE_ANCHORED,
NIE FIRST_PRINCIPLES.

**Remediation:**
- Either honestly classify cycle jako "literature-anchored verification" (output_type:
  structural; max claim status C)
- Or **wyprowadź** β-function z TGP axioms (Φ-EOM kowariantnej, 1-loop integral w
  `g_eff[{Φ_i}]` background, S05 single-Φ source) — to jest FIRST_PRINCIPLES

### §4.4 — Anti-pattern: Lakatos OR-clause werdykty

**Forbidden:**
```yaml
hypotheses:
  H1a: TGP-pure works
  H1b: TGP + sterile ν works (if H1a fails)  # ← Lakatos backstop
  H1c: STRUCTURAL_NO_GO
```

**Why bad:** decision tree z gwarantowanym ALL-PASS; H1b adopts post-hoc parameters
bez pre-bounded recovery_scope per PRE_REGISTERED_FALSIFIERS §3.3.

**Remediation:**
- H1a sole hypothesis; if H1a fails → EARLY_HALT_HONEST verdict (analog do
  op-MAG-anomalous-moment + cluster cycle Rec 2 precedent)
- H1b → **separate dedicated cycle** z explicit pre-bounded `recovery_scope.allowed_directions`
- PR-### entry dla H1b cycle przed Phase 1 sympy fit, NIE post-hoc

### §4.5 — Anti-pattern: "6/6 P-requirements RESOLVED" jako automatic verdict

**Forbidden:** declarative "6/6 RESOLVED" bez explicit per-P evidence i bez external
adversarial check.

**Why bad:** cohort 2026-05-11 wszystkie cykle uzyskały "6/6 RESOLVED" tego samego dnia —
to jest statistical anomaly, NIE positive result.

**Remediation:**
- Per-P evidence: każdy P-requirement linked do konkretnego sympy test lub Phase results section
- Honest verdict spectrum: 6/6, 5/6 + 1 documented, 4/6 + 2 deferred, 0/6 EARLY_HALT_HONEST
- Adversarial check przed Phase FINAL closure (spawn subagent z decydowalnym pytaniem)

### §4.6 — Anti-pattern: STRUCTURAL_DERIVED dla cycle bez observable target

**Forbidden:** claim `STRUCTURAL_DERIVED` (A−/A/A+) dla cycle z `output_type: structural`
(brak observable).

**Why bad:** per CYCLE_LIFECYCLE §taxonomy, `output_type: structural` → max claim_status
**C (STRUCTURAL_VERIFIED)**. `STRUCTURAL_DERIVED` impliuje falsifiable native observable.

**Remediation:** explicit observable target w `contract.L1_native.output_observable` lub
honest classification jako C.

---

## §5 — Recommended cycle order po restart

Per backlog z post-Rec-1+2+3+F outcome:

### §5.1 — Wysokie priorytety (autor decyzja, NIE auto-spawned)

1. **`op-LIGO-3G-native-phase-residual-2026-05-11`** — already validator PASS, ready dla
   user activation pending WIP slot
2. **Retrofit cycles N1-N5** — odpowiednio do retrofit candidates §1.3 tabeli; każdy
   jest nowym cyklem z BINDING template
3. **`op-cluster-sterile-nu-prediction-2026-XX`** — separate cycle z pre-bounded
   recovery_scope (m_ν, sin²2θ, ΔN_eff explicit bounds)

### §5.2 — Niskie priorytety (foundational, deferred)

4. **S07 reset cycle** (z explicit kickoff zgodny z new template) — post M9.1'' GWTC-3
   falsification, alternative f(ψ) families enumeration
5. **Inflation substrate genesis** — Φ_eq(t) prehistory; bardzo long-term (8-12 sesji)
6. **Composite Higgs sub-case** — hierarchy NO_GO follow-up

### §5.3 — Recommended workflow rate

Per CYCLE_LIFECYCLE WIP-limit (max 5 active): jeden cycle Phase 0 → Phase 1 → Phase
FINAL closure w 1-2 sesjach (compact) lub 5-8 sesjach (multi-phase), NIE 6 cykli
zamykanych w tym samym dniu (powtórzenie 2026-05-11 cohort anomaly).

---

## §6 — Sign-off

**Restart authorization:** autor projektu, conversation 2026-05-11, option (L)
"autoryzuję L, chciałbym zakończyć audyty i zacząć działać od nowa z tymi projektami
badawczymi, według nowego czystego schematu".

**Status:** BINDING dla wszystkich cykli post-2026-05-11.

**Wykonane deliverables:**

1. `tooling/validate_kickoff.py` — pure-stdlib Python validator (technical enforcement)
2. `meta/templates/op-cycle-kickoff-template-v2-2026-05-11.md` — minimal viable boilerplate
3. `meta/RESEARCH_RESTART_2026-05-11.md` (ten plik) — operational guidance
4. Scaffold #4 + #5 → `parking-pending-new-kickoff`
5. STATE.md update z halt notice + enforcement reference

**Pełna historia audytu z którego wynika restart:**

- [[AUDIT_2026-05-11_sympy_substance.md]] — niezależny adversarial audit 112 testów
- [[../PREDICTIONS_REGISTRY.md]] — STATUS DOWNGRADE PREAMBLE + Rec 2 + Rec 3+F updates
- Per-cycle §RETROACTIVE sections w 7 closure files (N1-N5 + cluster + hierarchy)

**Update protocol:** ten plik jest BINDING reference; updates wymagają explicit author
authorization + change log w §7 (nie yet created — to-be-added przy pierwszej revision).

**Next action expected:** user authorization dla pierwszego nowego cyklu na restart
schemacie. Recommended pierwszy candidate: `op-LIGO-3G-native-phase-residual-2026-05-11`
activation (already validator PASS).
