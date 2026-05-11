---
title: "CYCLE_LIFECYCLE — reguły statusu cykli badawczych research/op-*"
date: 2026-05-09
last_updated: 2026-05-10
type: meta-policy
status: ACTIVE (binding od 2026-05-09; extended 2026-05-10 z claim status taxonomy)
parent: "[[README.md]]"
related:
  - "[[../STATE.md]]"
  - "[[PLAN_RESEARCH_WORKFLOW_v1.md]]"
  - "[[CALIBRATION_PROTOCOL.md]]"
  - "[[CYCLE_KICKOFF_TEMPLATE.md]] (BINDING dla cykli post-2026-05-10)"
  - "[[PPN_AS_PROJECTION.md]]"
  - "[[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]]"
  - "[[M9_RESTRUCTURE_NOTE.md]]"
  - "[[VALIDATION_TRANSFERS.md]]"
  - "[[PRE_REGISTERED_FALSIFIERS.md]]"
supersedes: null
revisions:
  - 2026-05-09: Initial WIP-limit + folder_status słownik + Phase 0 README template
  - 2026-05-10: Extended z (a) claim status taxonomy A+/A/A-/B/C/D, (b) output_type
    field requirement, (c) observable-check anti-pattern #7, (d) kickoff template binding
---

# CYCLE_LIFECYCLE — reguły statusu cykli `research/op-*`

## Po co ten dokument

Diagnoza chaosu 2026-05-09 zidentyfikowała:

> **80 cykli z `folder_status: active`, realnie WIP ~5–10.**

`folder_status: active` straciło znaczenie operacyjne — oznaczało jedynie
"otwarty, nieformalnie nigdy nie zamknięty". Stąd reguły poniżej.

Ten dokument **nie zastępuje** `PLAN_RESEARCH_WORKFLOW_v1.md`
(który definiuje **fazy wewnątrz** cyklu: Phase 0..6+) — definiuje
**status zewnętrzny cyklu** (ile cykli równolegle żyje i w jakiej formie).

## Dwa poziomy statusu

| Poziom | Co opisuje | Gdzie | Kto aktualizuje |
|---|---|---|---|
| **Cycle status** (ten dokument) | Czy cykl żyje, śpi, czeka, jest zamknięty | `folder_status` w README.md cyklu + `STATE.md` | Po każdej sesji która dotyka cyklu |
| **Phase status** (`PLAN_RESEARCH_WORKFLOW_v1`) | Postęp wewnątrz cyklu (Phase 0..6) | `phase` w README.md cyklu + Phase\*_results.md | Po każdej zakończonej fazie |

## Słownik statusów (`folder_status` values)

| Status | Znaczenie | Warunek wejścia | Warunek wyjścia | Liczy się do WIP-limit? |
|---|---|---|---|---|
| `active` | Realnie pracuję nad nim w tej / następnej sesji | Critical path lub WIP slot wolny + decyzja user'a | Phase FINAL closed → `closed-*`, lub świadomy pivot → `paused` | **TAK** |
| `paused` | Świadomie zamrożony; blocker udokumentowany w README §Status | Decyzja "wracam później", blocker spisany | Blocker rozwiązany + WIP slot wolny → `active` | NIE |
| `needs-bridge` | Czeka na zewnętrzny precedens (inny cykl, decyzja, dane) | Zależność od `op-X CLOSED` lub external data | Bridge dostarczony → `active` (wymaga WIP slot) | NIE |
| `parking` | Pomysł / propozycja zarejestrowana, jeszcze nie startujemy | Otwarte przez user'a / zewnętrzną recenzję, brak Phase 0 | User decyzja "start" → `active` (wymaga WIP slot) | NIE |
| `closed-resolved` | Phase FINAL z verdict DERIVED / STRUCTURAL_DERIVED / STRUCTURAL_CONDITIONAL | Phase 6 ABSOLUTE BINDING gate PASS | Terminal status — nie zmienia się | NIE |
| `closed-NULL` | Phase FINAL z verdict EARLY_HALT (brak claimu, honest acknowledgment) | EARLY_HALT po Phase 0-2 z dobrym scope mapping | Terminal status | NIE |
| `closed-superseded` | Inny cykl objął zakres / zastąpił to ujęcie | Rozproszony zakres, nowy cykl owns problem | Terminal status; link do następcy MANDATORY w README §Status | NIE |
| `closed-FALSIFIED` | Hipoteza obalona empirycznie / strukturalnie z honest reporting | Phase 4-5 z observational ruling out | Terminal status | NIE |
| `archive` | Eksperymentalna struktura, system folder | Specjalne (np. `_archive`, `_sandbox`) | — | NIE |

## WIP-limit

> **Maksymalnie 5 cykli z `folder_status: active` w danym czasie.**

Wyjątki:

- **Critical path cycle** (oznaczony w STATE.md sekcja "Critical path") nie liczy
  się do WIP-limit, ale automatycznie blokuje 1 slot. Zwykle 0–1 critical path.
- **Same-day spawn z closure**: jeśli zamknięty cykl X spawnuje 1–4 sub-cykle Y_i
  jako natural continuation, mogą wszystkie być `active` przez 1–2 sesje (jak
  dual-V chain 2026-05-09: Phi-vacuum-scale spawned V-canonical + MAG-Phase5 +
  dual-V-clarification + Phase5-MAG-erratum w jednej fali). Po fali — zwykła reguła.

**Procedura zwiększenia WIP poza 5:**

1. Spisz w STATE.md uzasadnienie (np. cascade closure, external deadline).
2. Określ datę powrotu do limit-5.
3. Po tej dacie — przesuń najmniej krytyczne `active` do `paused`.

## Stale-detection (auto-pause candidate)

Jeśli cykl ma `folder_status: active` ale **brak commita w jego folderze >30 dni**:

- Domyślnie kandydat do `paused` (wymaga decyzji user'a, nie automatyczne).
- W STATE.md outstanding-debt sekcji odnotować jako "X stale-active cycles".
- Skrypt `tooling/check_stale_cycles.py` (do napisania) może co tydzień wypisywać
  listę kandydatów (read-only, nie modyfikuje YAML).

## Przykładowe transitions (z dual-V chain 2026-05-09)

| Cykl | Przejście | Powód |
|---|---|---|
| `op-Phi-vacuum-scale-2026-05-09` | (otwarcie) `parking` → `active` | User user iteration "kolektywny Schwarzschild" |
| `op-Phi-vacuum-scale-2026-05-09` | `active` → `closed-resolved` | Phase FINAL: STRUCTURAL_DERIVED_CONDITIONAL_HALT, 84/88 PASS |
| `op-Phase5-MAG-erratum-2026-05-09` | (spawn) `parking` → `active` (cascade exception) | Spawned z dual-V chain |
| `op-Phase5-MAG-erratum-2026-05-09` | `active` → `closed-resolved` | Phase 5 erratum applied (γ = m_C²), 5/5 PASS |
| `op-S07-alternative-f-psi-derivation-2026-05-09` | (otwarcie) — | Critical path; nie liczy się do WIP, ale blokuje 1 slot |

## Co aktualizować przy zmianie statusu

Każda zmiana `folder_status` cyklu wymaga:

1. **README.md cyklu** — pole `folder_status:` w YAML + sekcja §Status z datą i powodem.
2. **STATE.md** (root) — sekcja "Active WIP" lub "Recent closures" lub
   "Outstanding meta-debt" jeśli stale.
3. **PRIORITY_MATRIX.md** (audyt/) — jeśli cykl ma związek z S/L/D/M/T/EXT issue.
4. (opcjonalne) **PREDICTIONS_REGISTRY.md** — jeśli closure dotyka predykcji
   (FALSIFIED/PASS/PENDING/WITHDRAWN).

## Claim status taxonomy (post-2026-05-10)

**Cycle status (`folder_status`)** mówi *czy cykl żyje*. **Claim status (`claim_status`)**
mówi *jak silne jest twierdzenie cyklu*. Te są niezależne wymiary.

### Słownik `claim_status` (binding post-2026-05-10)

| Level | Tag | Definition | Wymagania |
|---|---|---|---|
| **A+** | `STRUCTURAL_DERIVED_NATIVE_WITH_TRANSFER` | Native L1 derivation + L2 framework reduction succeeded analitycznie → validation transfer aktywny | (1) `output_type: observable`, (2) full L1 native derivation, (3) L2 reduction analytical-exact lub analytical-approximate, (4) entry w [[VALIDATION_TRANSFERS.md]], (5) [[PRE_REGISTERED_FALSIFIERS.md]] entry |
| **A** | `STRUCTURAL_DERIVED_NATIVE` | Native L1 derivation + L2 mapping attempted i failed analytically, ale observables agree numerycznie (functional equivalence) | (1) `output_type: observable`, (2) full L1 native derivation, (3) L2 attempted z explicit failure note, (4) numerical agreement documented, (5) PR-### entry |
| **A−** | `STRUCTURAL_DERIVED_NATIVE_PARTIAL` | Native L1 derivation + L2 mapping not yet attempted (legitimate retrofit candidate) | (1) `output_type: observable`, (2) full L1 native derivation, (3) PR-### entry |
| **B** | `PROJECTION_VERIFIED` | L2-only output (parameter w obcym frameworku jako primary) — drift mode lub intentional translation | (1) `output_type: projection`, (2) explicit declaration: drift-deprecated lub intentional-translation; cykl NIE może być cytowany jako falsifiable native prediction |
| **C** | `STRUCTURAL_VERIFIED` | Internal algebraic consistency verified, brak observable target | `output_type: structural`; sympy LOCK ale brak L1 obserwabli |
| **D** | `SPECULATIVE_PARTIAL` | Phase 0-2 work-in-progress, nie aspirujący do closure claim | n/a — nie closing status |

**Hard rule:** `output_type: observable` jest wymagane dla A+/A/A-. Brak observable target lub
`output_type: projection` ogranicza max claim status do B.

**Hard rule:** PR-### entry w [[PRE_REGISTERED_FALSIFIERS.md]] jest wymagane dla A+/A/A-.
Bez immutable pre-registration timestamp → max C.

### Mapping claim_status ↔ verdict (legacy → taxonomy)

| Legacy verdict | Mapping na claim_status |
|---|---|
| `STRUCTURAL_DERIVED` (z observable + L2 transfer) | A+ |
| `STRUCTURAL_DERIVED` (z observable, no L2 lub L2 fail) | A or A− |
| `STRUCTURAL_DERIVED` (output = parameter, no observable) | **DOWNGRADE → B** (drift) lub **B** (intentional translation) |
| `STRUCTURAL_CONDITIONAL` (z observable, conditional path) | A− pending; or B if condition is L2-mapping match |
| `STRUCTURAL_CONDITIONAL_HALT` (conditional, blocked) | claim status uncommitted; closure as `closed-resolved` z `claim_status: pending-bridge` |
| `EARLY_HALT_HONEST` | C (no observable claim) lub `closed-NULL` |
| `STRUCTURAL_VERIFIED` (algebra) | C |

**Migration:** pre-2026-05-10 cykle zachowują legacy verdict tag w README; **dodanie**
`claim_status` field per audit (Phase 0-1 retrofit plan). Legacy verdict NIE jest usuwane —
przybywa nowe pole.

### Output type field (mandatory post-2026-05-10)

Każdy nowy cykl musi mieć w YAML:

```yaml
tgp_status:
  output_type: <observable|projection|structural>
```

- `observable` — output ma jednostki fizyczne (arcsec, Hz, ms, dimensionless strain ratio)
- `projection` — output to parameter w obcym frameworku (intentional translation lub drift)
- `structural` — output to algebraic consistency (np. axiom verification, sympy LOCK na
  identity)

**Pre-2026-05-10 cykle:** field dodawany retroactively per audit; default `pending-audit`
do triage decyzji.

---

## Anti-patterns (do unikania)

| # | Anti-pattern | Co zamiast |
|---|---|---|
| 1 | "Otwieram cykl ale w sumie nie wiem czy będę robić" | Otwórz w `parking`. Nie marnuj WIP slot. |
| 2 | Cykl `active` od 60+ dni bez commita | `paused` z spisanym blocker'em. |
| 3 | 5+ cykli `active` od jednego user'a w jednej sesji | Skup się na critical path; resztę → `parking`. |
| 4 | Zamknięty cykl ale `folder_status` dalej `active` | Zaktualizuj na `closed-*` przy commicie zamknięcia. Closure marker + status muszą się zgadzać. |
| 5 | Handoff w roocie zamiast cyklu | Każdy handoff >1 sesja = cykl `op-*`. (Patrz STATE.md migration log 2026-05-09.) |
| 6 | Spawn 5+ sub-cykli równocześnie bez cascade-exception entry w STATE.md | Spisz cascade exception przed spawningiem. |
| **7** | **L2-only output without L1 native derivation (drift mode)** — cykl produkuje *parameter w obcym frameworku* (β_ppE, β_PPN, γ_PPN, ξ_2, ...) jako primary, brak observable target z fizycznymi jednostkami | **Auto-downgrade do `claim_status: B (PROJECTION_VERIFIED)`**. Per [[CYCLE_KICKOFF_TEMPLATE.md]] §2.2. Cykl nie może aspirować do A+/A/A-. Jeśli intencjonalna projekcja: explicit `output_type: projection` + `kind: framework-translation`. Jeśli drift: retrofit do native form lub archive jako MIMICRY_DEPRECATED. |
| **8** | **Brak `pre_registration_date` dla falsifiable claim** | Per [[PRE_REGISTERED_FALSIFIERS.md]] §3.4: bez immutable timestamp moving-goalposts ryzyko aktywne. Cykl claim status max C (internal consistency). Z entry: A+/A/A-. |
| **9** | **L2-projection presented jako "TGP prediction"** zamiast jako consistency check | Reframe: native L1 jest predykcją; L2 to chart compatibility. Per [[PPN_AS_PROJECTION.md]] §3.1 three-layer mandatory. |

## Phase 0 README template (BINDING post-2026-05-10)

**Wszystkie cykle z 2026-05-10+ otwierane MUSZĄ mieć w README.md sekcję `## §X — TGP-native
check (mandatory, pre-Phase-1)`** zawierającą poniższy checklist (per
[[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §5):

```markdown
## §X — TGP-native check (mandatory, pre-Phase-1)

[ ] **Q1 (Pattern coverage):** Czy reviewuję §2 patterns z TGP_NATIVE_COMPUTATIONAL_PATTERNS.md
        relevant dla mojego problemu?
        Lista relevantnych patterns: __________________

[ ] **Q2 (Red flags):** Czy zidentyfikowałem §3 red flags w moim Phase 1 plan?
        Lista znalezionych red flags: __________________
        Jeśli brak: explicit explanation czemu plan jest red-flag-free: __________________

[ ] **Q3 (Inherited LOCKs §4 mapping):** Lista inherited LOCKs i ich §4 entry status:
        - LOCK 1: ____________  →  §4 entry F#: ____  status: ✅/⚠️/❌
        - LOCK 2: ____________  →  §4 entry F#: ____  status: ✅/⚠️/❌
        Jeśli ANY ❌ lub brak entry: WYKONAJ ASK-RULE Trigger B.

[ ] **Q4 (Standard-physics tools):** Czy mój Phase 1 plan używa standardowych tools
        (Yukawa, BD ω, GR derivation, propagator, exchange)?
        Jeśli TAK: explicit justify czemu NOT BD-translation
        (link do §2 pattern lub §4 mapping): __________________

[ ] **Q5 (m_Φ usage):** Czy m_Φ używam jako universal stałą czy environment-dependent
        (Pattern 2.5)?
        Jeśli universal: justify or use Pattern 2.5 environment-dependent: __________________

[ ] **Q6 (GR limit framing):** Czy plan distinguishes "TGP gives same number as GR
        (via TGP mechanism)" vs "TGP IS GR by translation"?
        Mechanism cite: __________________

[ ] **Q7 (ASK-RULE self-check):** Czy są niewyjaśnione gaps gdzie sięgam po std physics?
        Jeśli TAK: WYKONAJ ASK-RULE before Phase 1 sympy.
        Triggers fired: __________________
        Asked user: yes / no  (jeśli no: explain czemu wyjątki §1.4 applies)

[ ] **Q8 (BD-drift audit plan):** Phase FINAL będzie spawn `bd-drift-audit` subagent
        per CALIBRATION_PROTOCOL §4.4? (yes / no — jeśli no, justify per §4.4.5 fallback)
```

**Cykle bez §X TGP-native check:** automatically `STRUCTURAL_CONDITIONAL z BD-drift-audit-pending`
(per CALIBRATION_PROTOCOL §4.4.3 verdict consequences).

**Pre-2026-05-10 cykle:** zachowane bez §X (nie retroaktywnie wymuszane), ALE inheriting
LOCKs z tych cykli WYMAGA ASK-RULE Trigger B sprawdzenia (per
[[TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1.1 Trigger B definicja).

---

## Mapowanie z legacy `folder_status` values (sprzed 2026-05-09)

Dla cykli istniejących przed wprowadzeniem tej polityki:

| Legacy value | Nowy value |
|---|---|
| `active` (z commitem <14 dni) | sprawdź czy w WIP-5 → `active`; reszta → `paused` |
| `active` (>30 dni bez commita, brak Phase FINAL) | `paused` |
| `active` (z Phase_FINAL/Phase6_absolute_binding marker) | `closed-resolved` lub `closed-NULL` (zależne od verdict) |
| `needs-bridge` | bez zmian, ale verify że bridge dependency dalej istnieje |
| `research`, `research-active` | `active` (do WIP-5) lub `paused` |
| `audit` (wewnątrz research/) | rozważyć przeniesienie do `audyt/` |
| `archive`, `sandbox`, `review` | bez zmian (system) |
| null (brak statusu) | `parking` lub `paused` zależnie od kontekstu; flag dla manual review |

## Migracja 2026-05-09 (one-time)

Inwentaryzacja 2026-05-09 (z `STATE.md` outstanding-debt #4): 80 cykli oznaczonych
`active`, realnie ~5–10. Pełna reklasyfikacja: osobna sesja. Plan:

- Bucket A (19 cykli z commitem ostatnie 14 dni): triage do WIP-5 + reszta `paused`
- Bucket B (3 cykle z Phase6 marker, dalej oznaczone `active`): `closed-resolved`
- Bucket C (91 stale-active): zbiorczo → `paused`
- Bucket D (6 needs-bridge): verify dependencies
- Bucket E (10 unknown): manual triage

**Reguła:** mass-rewrite YAML 91 plików wymaga explicit user authorization
(blast radius — można uszkodzić strukturę cyklu). Preferowane: skrypt
`tooling/reclassify_cycles_2026-05-09.py` z dry-run + diff preview.
