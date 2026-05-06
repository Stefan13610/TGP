---
title: "RESUME PROTOCOL — instrukcja dla agenta w przyszłej sesji M03"
date: 2026-05-06
parent: "[[README.md]]"
type: protocol
tgp_owner: research/op-M03-balance-sheet-retrofit-2026-05-06
tags:
  - protocol
  - resume
  - future-agent
  - anti-duplication
---

# RESUME PROTOCOL — M03 multi-session continuation

## Cel tego dokumentu

Instrukcja dla **przyszłego agenta** kontynuującego pracę M03 w kolejnej
sesji. Cel: **brak dublowania pracy** + spójna kontynuacja.

## Krok 1: Przeczytaj te pliki w tej kolejności

```
1. research/op-M03-balance-sheet-retrofit-2026-05-06/README.md       — master plan
2. research/op-M03-balance-sheet-retrofit-2026-05-06/tracker.md      — status każdego cyklu
3. research/op-M03-balance-sheet-retrofit-2026-05-06/audit_log.md    — co już zrobione
4. research/op-M03-balance-sheet-retrofit-2026-05-06/high_risk_queue.md — następna kolejka
```

## Krok 2: Sprawdź tracker

W [[tracker.md]] każdy cykl ma status:

- `PENDING` — nie tknięty, dostępny do retrofit
- `IN_PROGRESS` — sesja zaczęła ale nie skończyła (sprawdź audit_log dla powodu)
- `DONE_<class>` — zaaudytowany, sklasyfikowany (DERIVED_FULL / CONDITIONAL / STRUCTURAL / ANSATZ / NUMEROLOGICAL / TAUTOLOGY)
- `BLOCKED_<reason>` — np. czeka na decyzję autora lub inny cykl

## Krok 3: Wybierz następny cykl

**Kolejność priorytetów:**

1. **High-risk pending** ([[high_risk_queue.md]]) — pattern mixing-operator + cascade
2. **IN_PROGRESS resume** — dokończ sesję poprzedniego agenta (sprawdź audit_log
   ostatni wpis)
3. **Medium-risk** — claim DERIVED z niezależną fizyką
4. **Low-risk** — już post-audit downgraded

## Krok 4: Wykonaj retrofit per cykl

Per cykl `op-X`:

1. **Skopiuj template:** `template_Phase0_balance.md` →
   `retrofit_op-X_<dzisiaj>.md`
2. **Wypełnij sekcje 2.1-2.6** zgodnie z [[../../meta/CALIBRATION_PROTOCOL.md]] §2
3. **Wykonaj 3 testy:**
   - Tautology test (sympy substitution → outputs się kasują definicyjnie?)
   - Falsifiability test (axiom value które wykluczyłoby match?)
   - Independent-path cross-validation (≥2 paths convergent?)
4. **Klasyfikuj** zgodnie z tabelą w [[README.md]] §"Klasyfikacja klas epistemicznych"
5. **Update tracker.md** — zmień status cyklu na `DONE_<class>`
6. **Update audit_log.md** — dodaj wpis z datą, klasyfikacją, krótkim rationale

## Krok 5: Co ROBIĆ a czego NIE

### ROBIĆ ✅

- Tworzyć nowe pliki `retrofit_<cykl>_<data>.md` w
  `research/op-M03-balance-sheet-retrofit-2026-05-06/`
- Update `tracker.md` (modyfikuj status, datę, klasyfikację)
- Update `audit_log.md` (dodawaj nowe wpisy chronologicznie)
- Tworzyć `CRITIQUE_<issue>_<date>.md` w cyklu jeśli wykryjesz tautologię
  (zgodnie z [[../../meta/CALIBRATION_PROTOCOL.md]] §4)

### NIE ROBIĆ ❌

- **NIE** modyfikuj plików oryginalnego cyklu (`research/op-X/README.md`,
  `phase*_results.md`, etc.) — pozostają immutable jako historic record
- **NIE** zmieniaj statusu cyklu w jego YAML frontmatter (to robi Phase 5
  registry refactor po pełnym retrofit)
- **NIE** usuwaj entries z `PREDICTIONS_REGISTRY.md` (nawet jeśli okażą się
  TAUTOLOGY — zachowane jako historic + downgraded w post-Phase 5 refactor)
- **NIE** dubluj retrofit dla cyklu z `tracker.md` status `DONE_*`
- **NIE** rozpoczynaj nowego cyklu jeśli jest IN_PROGRESS od poprzedniej
  sesji (najpierw dokończ to)

## Krok 6: Kiedy zakończyć sesję

Sesja M03 może być długa. Sygnały by zakończyć:

- Kontekst > 70% pełny → save state, append do audit_log "session pause"
- Audyt 1 cyklu zakończony → naturalny moment by przerwać
- Wykryta poważna tautologia → zatrzymaj się, dokumentuj w `CRITIQUE_*`

**Always:** przed końcem sesji update `tracker.md` i `audit_log.md`.

## Krok 7: Sygnał ukończenia M03

M03 jest UKOŃCZONY gdy:

1. Wszystkie cykle z `tracker.md` mają status `DONE_*` (nie ma PENDING/IN_PROGRESS)
2. Phase 5 registry refactor zakończony
3. CALIBRATION_PROTOCOL absolute-binding gate enforced

Po tym: `audyt/M03_balance_sheet_missing/POST_ACTION_FINAL_<data>.md` +
`audyt/PRIORITY_MATRIX.md` update do `EXECUTED ALL PHASES`.

## Anti-pattern: czego unikać

### Anti-pattern 1: "I'll re-audit cycle X"

Jeśli `tracker.md` mówi `DONE_DERIVED_FULL`, **NIE rób re-audit** chyba że
masz konkretny dowód błędu w retrofit_<X>_<data>.md. Zamiast tego:
- Dodaj `CRITIQUE_<issue>_<data>.md` w folderze cyklu z explicit problem
- Update `tracker.md` z `DONE_X → DONE_X_DISPUTED_<data>`

### Anti-pattern 2: "Wszystkie cykle są podobne"

NIE używaj template'u na "Easy mode" bez explicit verification. Każdy cykl
wymaga **własnego** tautology test, falsifiability test, independent paths.
Mechaniczne wypełnianie template = spaprować retrofit.

### Anti-pattern 3: "Skipuję low-risk"

NIE skipuj low-risk cykli. Jeśli były zaklasyfikowane jako low-risk w
[[high_risk_queue.md]] to dlatego że dotychczas wyglądały OK — ale formal
balance sheet może ujawnić niespodzianki. Każdy cykl claiming LOCKED/DERIVED
**musi** mieć retrofit.

### Anti-pattern 4: Modyfikacja core LaTeX

M03 jest **meta-audit**, NIE core edit. Jeśli retrofit ujawnia że cykl X
miał błąd który prevoz wpłynął na core LaTeX (rzadkie ale możliwe):
- Dokumentuj w `retrofit_<X>_<data>.md` §"Core impact"
- NIE edytuj core (poza dodaniem comment-block)
- Eskaluj do osobnego cyklu naprawczego

## Konkretny workflow dla nowej sesji

```
1. cat research/op-M03-balance-sheet-retrofit-2026-05-06/audit_log.md | tail -20
   → ostatnie 5-10 wpisów, zrozum gdzie poprzednia sesja skończyła
2. Read tracker.md → identyfikuj IN_PROGRESS lub PENDING
3. Wybierz cykl (priorytet z high_risk_queue.md)
4. Read research/op-X/README.md + Phase*_results.md (~5-10 min)
5. Wypełnij retrofit_op-X_<dzisiaj>.md (~30-60 min)
6. Update tracker.md (1 min)
7. Update audit_log.md (1 min)
8. Repeat 3-7 dla kolejnego cyklu (jeśli czas)
9. Save state → naturalne zakończenie sesji
```

## Cross-references

- [[README.md]] — master plan
- [[tracker.md]] — status tracking
- [[template_Phase0_balance.md]] — template
- [[audit_log.md]] — log
- [[high_risk_queue.md]] — priority queue
- [[../../meta/CALIBRATION_PROTOCOL.md]]
