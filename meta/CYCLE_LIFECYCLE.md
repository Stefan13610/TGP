---
title: "CYCLE_LIFECYCLE — reguły statusu cykli badawczych research/op-*"
date: 2026-05-09
type: meta-policy
status: ACTIVE (binding od 2026-05-09)
parent: "[[README.md]]"
related:
  - "[[../STATE.md]]"
  - "[[PLAN_RESEARCH_WORKFLOW_v1.md]]"
  - "[[CALIBRATION_PROTOCOL.md]]"
supersedes: null
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

## Anti-patterns (do unikania)

| # | Anti-pattern | Co zamiast |
|---|---|---|
| 1 | "Otwieram cykl ale w sumie nie wiem czy będę robić" | Otwórz w `parking`. Nie marnuj WIP slot. |
| 2 | Cykl `active` od 60+ dni bez commita | `paused` z spisanym blocker'em. |
| 3 | 5+ cykli `active` od jednego user'a w jednej sesji | Skup się na critical path; resztę → `parking`. |
| 4 | Zamknięty cykl ale `folder_status` dalej `active` | Zaktualizuj na `closed-*` przy commicie zamknięcia. Closure marker + status muszą się zgadzać. |
| 5 | Handoff w roocie zamiast cyklu | Każdy handoff >1 sesja = cykl `op-*`. (Patrz STATE.md migration log 2026-05-09.) |
| 6 | Spawn 5+ sub-cykli równocześnie bez cascade-exception entry w STATE.md | Spisz cascade exception przed spawningiem. |

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
