---
title: "Patch log — CORE-CLEANUP-B (2026-05-04)"
date: 2026-05-04
parent: "[[README.md]]"
type: patch-log
tgp_owner: research/op-CORE-CLEANUP-B-2026-05-04
tags:
  - patch-log
  - core-cleanup
  - non-breaking
---

# Patch log — CORE-CLEANUP-B 2026-05-04

## Sesja 1 — 2026-05-04 (Phase 0 + Phase 1 + Phase 2 partial)

### Phase 0: Setup cyklu

| Akcja | Plik | Linie zmienione |
|-------|------|------------------|
| Utworzenie folder cyklu | `research/op-CORE-CLEANUP-B-2026-05-04/` | NEW |
| README.md | `research/op-CORE-CLEANUP-B-2026-05-04/README.md` | NEW |

### Phase 1: Skeleton dodatek_pivot_history.tex

| Akcja | Plik | Linie zmienione |
|-------|------|------------------|
| Utworzenie skeleton | `core/formalizm/dodatek_pivot_history.tex` | NEW (~200 lin) |
| Include do main.tex | `main.tex` | +5 lin (po dodatekZ) |

**Compile risk:** brak — nowy plik, samodzielnie compile clean.

### Phase 2 step 1: sek08c forma (III) antypodalna → appendix

| Akcja | Plik | Linie |
|-------|------|--------|
| Wyciąć thm:antipodal-uniqueness | `core/sek08c_metryka_z_substratu.tex` | usunięto lin. 491-586 (oryginalne) |
| Zastąp short remark | `core/sek08c_metryka_z_substratu.tex` | dodano lin. 491-525 (rem:form-III-relocated) |
| Skopiować pełną treść do appendix | `core/formalizm/dodatek_pivot_history.tex` | rozszerzono ~150 lin |

**Effect:**

- `sek08c_metryka_z_substratu.tex`: 592 → 544 linii (−48 lin)
- `dodatek_pivot_history.tex`: 200 → 356 linii (+150 lin pełnej treści twierdzenia)

**Cross-references zachowane** (wszystkie labels w appendix):

- `thm:antipodal-uniqueness` ✓
- `rem:antipodal-implications` ✓
- `eq:antipodal-metric-ansatz` ✓
- `eq:antipodal-unique-solution` ✓
- `eq:anti-c1`, `eq:anti-c2`, `eq:anti-general-solution`, `eq:anti-gamma-of-s` ✓

**Compile risk:** zostaje pdflatex check w Phase 5.

## Status post-Sesja 1

| Phase | Plan | Status |
|-------|------|--------|
| 0 — Setup | utworzenie cyklu | **DONE** ✓ |
| 1 — Skeleton dodatek | skeleton + main.tex include | **DONE** ✓ |
| 2 — sek08c cleanup | formy I, II, III → appendix | **PARTIAL** (III done; I+II pending) |
| 3 — sek08a cleanup | v1.x → appendix; body → v2.0 | NEXT SESSION |
| 4 — Cross-ref audit | ~30 P33 files | NEXT SESSION |
| 5 — Compile check | pdflatex full | NEXT SESSION |

## Pozostałe deprecated content do refactoru (NEXT SESSIONS)

### sek08c forma (II) eksponencjalna

**Lokalizacja:** lin. 326-355 (deklaracja `f=e^{-2U}, h=e^{+2U}`),
lin. 357-396 (Krok 3 derivation z formy eksponencjalnej).

**Labels:** brak (forma II nie ma własnego label, jest in-line w narracji).

**⚠ KOMPLIKACJA wykryta 2026-05-04:** narracja formy (II) jest **wpleciona
w centralny proof** `thm:metric-from-substrate-full` (lin. 305-397). Cały
proof używa łańcucha: forma (I) → β_PPN=2 → forma (II) eksponencjalna →
β_PPN=1. Wycięcie formy (II) wymaga **przepisania całego proofu** z
M9.1'' canonical jako jedyną formą.

**Strategia next session:** dwuetapowa:

1. **Etap A**: Skopiować pełny `thm:metric-from-substrate-full` (lin. 262-397)
   do `dodatek_pivot_history.tex` Section 2 jako `thm:metric-from-substrate-full-HISTORICAL`
   z explicit STATUS marker + zachowane wszystkie labels
   (`thm:metric-from-substrate-full`, `eq:metric-full-derived`, etc.)
2. **Etap B**: Zastąpić body sek08c (lin. 262-397) NEW theorem
   `thm:metric-from-budget-M911-canonical` z M9.1'' canonical derivation:
   - Krok 1: bilans budżetowy `fh=1` (zachowane bez zmian, ~poprawne dla M9.1'')
   - Krok 2: M9.1'' jako kanoniczna forma `g_tt = -c²(4-3ψ)/ψ`
   - Krok 3: PPN matching z master formula (β=γ=1 EXACT)
   - Krok 4: `√(-g) = c·ψ/(4-3ψ)` (M9.1'' canonical)

**Estymata:** ~2-3h refactor + compile check (większa niż wcześniej zakładana).

### sek08c forma (I) potęgowa

**Lokalizacja:** lin. 131-200 (preliminary proposition `prop:antipodal-from-budget`,
**zachowane** — jest podstawą fizyczną), lin. 262-303 (`thm:metric-from-substrate-full`
z `eq:metric-full-derived`), lin. 432-488 (volume element `eq:vol-element-derived`).

**Labels:** `eq:metric-full-derived`, `prop:antipodal-from-budget`,
`thm:metric-from-substrate-full`, `eq:vol-element-derived`,
`rem:metric-action-consistency`.

**Strategia next session:** trudniejsze — `thm:metric-from-substrate-full`
jest podstawowym mostem substrat→metryka. Można:

- (a) **Refactor** twierdzenia z M9.1'' (zachować label, zmienić treść)
- (b) **Relocate** do appendix + nowy theorem dla M9.1'' canonical

**Strategia preferowana:** (b) — relocate, zachować label, dodać nowe
canonical theorem w body.

**Estymata:** ~1h refactor + compile check.

### sek08a deprecated v1.x propositions

**Lokalizacja:** lin. 1-805 (cała sekcja v1.x).

**Labels:** `eq:V-selfinterference`, `eq:sqrt-g-eff`, `prop:vacuum-condition`,
`prop:kappa-corrected`, `eq:full-field-alt`, prawdopodobnie kilka więcej.

**Strategia next session:** najtrudniejsze — 805 linii z dużą ilością
propositions. Większość ma kanoniczne replacements w v2.0 ADDENDUM.
Strategia:

1. Inwentaryzować pełną listę v1.x labels (grep)
2. Skopiować wszystkie do `dodatek_pivot_history.tex` Section 1
3. Zastąpić v1.x body short remark referencing appendix + v2.0
4. Compile check

**Estymata:** ~3-4h refactor + extensive compile check.

## Honest scope assessment

**Co wykonano w sesji 2026-05-04:**

- ✅ Phase 0: setup cyklu + folder
- ✅ Phase 1: skeleton `dodatek_pivot_history.tex` + main.tex include
- ✅ Phase 2 step 1: forma (III) antypodalna → appendix (z labels)
- ✅ Patch log + scope assessment

**Liczbowo (po sesji 1, faza 2.1):**

- `sek08c.tex`: 592 → 544 linii (−48 linii, forma III usunięta z body)
- `dodatek_pivot_history.tex`: 0 → 356 linii (utworzony z 3 sekcjami,
  forma III pełna treść skopiowana)
- `main.tex`: +5 linii (include nowego dodatku)

**Liczbowo (po sesji 1, faza 3.1 — Phase 3 step 1+2 wykonane):**

- `sek08a.tex`: 1021 → **858 linii (−163 lin)**, proof prop:kappa-corrected
  (~210 linii deprecated proof) relocated do appendix; zachowane:
  declaration prop:kappa-corrected + BOXED `eq:kappa-corrected` value
- `dodatek_pivot_history.tex`: 356 → **485 linii (+129 lin)** — pełen proof
  prop:kappa-corrected z 12 labels skopiowany
- `main.tex`: bez zmian (już zawiera include z fazy 1)

**Łączne zyski cleanup body:**

- sek08a + sek08c: **−211 linii deprecated content z body** (prop:kappa-corrected proof + forma III antypodalna)
- dodatek_pivot_history: +540 linii preserved historical record
- **Wszystkie cross-references zachowane** (12+ labels w appendix bez duplikatów; explicit avoidance multiply-defined errors)

**Dodatkowe edycje sek08a (sesja 1, faza 3.2):**

- prop:field-eq-from-action proof (lin. 261-368, ~108 linii): dodane explicit
  STATUS marker post-proof + canonical replacement reference (`prop:psi-EOM-R3`,
  `eq:R3-ODE`). Body proof zachowany dla cross-ref bezpieczeństwa; appendix
  ma structural summary bez duplicate labels.

**Liczbowy stan finalny po sesji 1+2:**

| Plik | Pre-sesja | Post-sesja 1+2 | Δ |
|------|-----------|------------|---|
| `sek08a.tex` | 1021 | **876** | **−145 lin** |
| `sek08c.tex` | 592 | **544** | **−48 lin** |
| `dodatek_pivot_history.tex` | 0 | **540** | **+540 lin** (NEW) |
| `main.tex` | (orig) | +5 | +5 lin (include) |

**Liczbowy stan finalny po sesji 3:**

| Plik | Post-sesja 1+2 | Post-sesja 3 | Δ sesja 3 |
|------|------------|----------|---|
| `sek08a.tex` | 876 | **876** | bez zmian |
| `sek08c.tex` | 544 | **682** | +138 (NEW canonical M9.1'' 76 lin + stary block w `\iffalse` 150 lin + separator 8 lin) |
| `dodatek_pivot_history.tex` | 540 | **641** | +101 (formy I+II expanded z pełną treścią) |

**Liczbowy stan finalny po sesji 4:**

| Plik | Post-sesja 3 | Post-sesja 4 | Δ sesja 4 |
|------|------------|----------|---|
| `sek08a.tex` | 876 | **981** | +105 (NEW `rem:hyp-unified-action-M911-canonical` z `eq:S-TGP-unified-M911-canonical` ~75 lin + DEPRECATED markers przy 3 equations ~30 lin) |
| `sek08c.tex` | 682 | 682 | bez zmian |
| `dodatek_pivot_history.tex` | 641 | 641 | bez zmian |

**Strategia sesji 4:** annotation-level cleanup z explicit `[DEPRECATED]`
markers w-place przy 3 equations + NEW canonical M9.1''-spojny hipotezę
jako `rem:hyp-unified-action-M911-canonical` (z BOXED
`eq:S-TGP-unified-M911-canonical`). Body cross-refs zachowane (5 \eqref
użyć w sek08a body + 1 w op-phase1-covariant); literal removal odroczone
do Phase 4 cross-ref audit (jeśli pewne że unique cross-refs).

## Sesja 5 — Cross-ref audit (2026-05-04, kontynuacja)

### Globalny grep: 30 plików zawierających deprecated labels

Inwentaryzacja `eq:V-selfinterference`, `eq:g-eff-unified`,
`eq:sqrt-g-eff`, `prop:field-eq-from-action`, `eq:field-eq-reproduced`,
`prop:kappa-corrected`, `eq:kappa-corrected`,
`thm:metric-from-substrate-full`, `eq:metric-full-derived`,
`thm:antipodal-uniqueness`, `rem:antipodal-implications`,
`eq:vol-element-derived`, `rem:metric-action-consistency`,
`rem:step3-elimination` w całym TGP_v1.

**Wyniki:**

| Plik | Typ | Kategoria | Status |
|------|-----|-----------|--------|
| `sek08a.tex` | core | host of `eq:sqrt-g-eff`, `prop:kappa-corrected`, etc. | OK (active labels) |
| `sek08c.tex` | core | host of `eq:metric-full-derived`, `eq:vol-element-derived` w `\iffalse` | OK (compile-skipped) |
| `dodatek_pivot_history.tex` | core | host of all relocated labels | OK (active labels) |
| `sek08_formalizm.tex` | core | 2× cross-ref do `prop:kappa-corrected` | OK (label aktywny w sek08a) |
| `status_map.tex` | core | 1× cross-ref do `prop:kappa-corrected` + `prop:kappa-corrected-G0` | OK |
| `sek00_summary.tex` | core | (informacyjne wzmianki) | OK |
| `sek09_cechowanie.tex` | core | (cross-refs do canonical sek08a) | OK |
| `dodatekH_lancuch_wyprowadzen.tex` | core | (cross-refs) | OK |
| `dodatekO_u1_formalizacja.tex` | core | (cross-refs) | OK |
| `op-phase1-covariant/Phase1_program.md` | research/.md | 1× `\eqref{eq:V-selfinterference}` + `\eqref{eq:g-eff-unified}` | OK info-level (md, nie LaTeX cross-ref) |
| `op-g0-r3-from-canonical-projection/*.md` | research/.md | meta-references | OK info-level |
| Audit files (`audyt/*.md`) | meta | meta-references do labels | OK info-level |

### Active LaTeX cross-refs do deprecated labels w body

**Compile-binding (`\eqref{}` lub `\ref{}` w `.tex` poza `\iffalse`):**

1. **`prop:kappa-corrected`** — 4+ active cross-refs:
   - sek08a body: lin. 319 (`stw.~\ref{prop:kappa-corrected}`),
     567 (short proof remark), 569, 592, 629
   - `sek08_formalizm.tex`: 2× (lin. 1810, 4132)
   - `status_map.tex`: 1× (lin. 254)
   → **Label MUSI być zachowany w body sek08a** (lin. 507).
   ✅ Status: zachowany.

2. **`eq:kappa-corrected`** — używane w sek08a body lin. 569 short proof
   remark + sek08c lin. 451 i 482 wewnątrz `\iffalse` (compile-skipped).
   → **Label aktywny w body sek08a** (BOXED VALUE INVARIANT, lin. 519).
   ✅ Status: zachowany.

3. **`eq:sqrt-g-eff`** — 5 active cross-refs w sek08a body:
   - lin. 294 (`rem:volume-element` opening)
   - lin. 313 (rem:volume-element comparison)
   - lin. 379 (`prop:field-eq-from-action` body)
   - lin. 517 (prop:kappa-corrected new short remark)
   - lin. 653 (v2.0 ADDENDUM porównanie z eksponencjalną)
   → **Label MUSI być zachowany** w body sek08a (lin. 179, z explicit DEPRECATED marker).
   ✅ Status: zachowany.

4. **`eq:V-selfinterference`** — 0 active LaTeX cross-refs w `.tex`;
   1× w `op-phase1-covariant/Phase1_program.md` (info-level, md).
   → **Bezpieczne kandydat dla literal removal w przyszłości**.
   Status: zachowane z explicit DEPRECATED marker (sesja 4).

5. **`eq:g-eff-unified`** — 0 active LaTeX cross-refs w `.tex`;
   1× w `op-phase1-covariant/Phase1_program.md` (info-level, md).
   → **Bezpieczne kandydat dla literal removal w przyszłości**.
   Status: zachowane z explicit DEPRECATED marker (sesja 4).

6. **`prop:field-eq-from-action`, `eq:field-eq-reproduced`,
   `eq:S-unified-psi`, `eq:EL-unified-raw`, `eq:EL-unified-divided`,
   `eq:S-static-correct`** — wszystkie aktywne w body sek08a (lin. 247,
   253, 264, 314, 326, 355). Cross-refs głównie wewnętrzne sek08a body.
   ✅ Status: zachowane (proof body z STATUS marker, sesja 2).

7. **`thm:metric-from-substrate-full`, `eq:metric-full-derived`,
   `rem:step3-elimination`, `eq:vol-element-derived`,
   `rem:metric-action-consistency`** — wszystkie cross-refs w sek08c
   body są wewnątrz `\iffalse` (lin. 389-624, compile-skipped).
   Active labels są tylko w appendix `dodatek_pivot_history.tex`.
   ✅ Status: compile-safe (no duplicates).

8. **`thm:antipodal-uniqueness`, `rem:antipodal-implications`**,
   plus 6 equations (`eq:antipodal-metric-ansatz`,
   `eq:antipodal-unique-solution`, `eq:anti-c1`, `eq:anti-c2`,
   `eq:anti-general-solution`, `eq:anti-gamma-of-s`) — wszystkie tylko
   w appendix (sesja 1 relocated).
   ✅ Status: compile-safe.

9. **12 equations z proof prop:kappa-corrected** (`eq:volume-FRW`,
   `eq:Lkin-FRW`, `eq:S-FRW-reduced`, `eq:ddt-pdotpsi`,
   `eq:dpsi-direct`, `eq:EL-FRW-full`, `eq:EL-FRW-simplified`,
   `eq:EL-FRW-divided`, `eq:psi-eq-unified`,
   `eq:cosmo-linearized-unified`, `eq:poisson-frw`,
   `eq:kappa-def-operational`) — wszystkie tylko w appendix (sesja 2 relocated).
   ✅ Status: compile-safe.

### Wnioski sesji 5

**Compile-safety verified:** brak duplicate labels, cross-refs prawidłowo
działają.

**Kandydaci dla literal removal w przyszłości** (po dodatkowych audit):
- `eq:V-selfinterference` (0 active LaTeX cross-refs)
- `eq:g-eff-unified` (0 active LaTeX cross-refs)

**Labels MUSZĄ być zachowane** (active cross-refs):
- `prop:kappa-corrected`, `eq:kappa-corrected`, `eq:sqrt-g-eff`,
  `prop:field-eq-from-action`, `eq:field-eq-reproduced`

**Pozostałe deprecated labels** są aktywne tylko w appendix
(compile-safe, no duplicates).

### Zmiany w sesji 5

**Brak edycji rdzenia** — sesja 5 to comprehensive cross-ref audit z
verification compile-safety. Wszystkie cross-refs działają prawidłowo
przez kombinację:
- Aktywnych labels w body sek08a (5 deprecated labels zachowane)
- `\iffalse` block w body sek08c (5 deprecated labels compile-skipped)
- Aktywnych labels w appendix (25 relocated labels)

### Status post-sesja 5

| Klastr | Status |
|--------|--------|
| Body sek08a active labels | ✅ 5 deprecated zachowane (active cross-refs) |
| Body sek08c `\iffalse` block | ✅ 5 deprecated compile-skipped |
| Appendix labels | ✅ 25 deprecated active w `dodatek_pivot_history.tex` |
| Cross-refs compile-safety | ✅ Verified (no duplicates) |
| Phase 6: pdflatex compile clean check | ✅ **PASS sesja 6** (547 stron, 0 nowych errors) |

## Sesja 6 — pdflatex compile clean check (2026-05-04, finalne)

### Wynik

**Status:** ✅ **PASS** (patrz [[POST_PATCH_VERIFICATION.md]] dla pełnych szczegółów)

- **PDF wygenerowany:** 547 stron (G.0 baseline 537 + 10 stron z cleanup)
- **0 nowych compile errors** wprowadzonych przez cleanup
- **3 multiply-defined labels** — wszystkie pre-existing identyczne z G.0 closure
  (`eq:MS-TGP`, `rem:formulation-dictionary`, `ssec:epistemic-table`)
- **8 undefined references** — wszystkie pre-existing identyczne z G.0 closure
  (`ax:substrat`, `app:A-aksjomaty`, etc.)

### Wykryta i naprawiona literówka w sesji 6

W `dodatek_pivot_history.tex` lin. 288: `\end{enumerate>` → `\end{enumerate}`
(błędny char `>` zamiast `}`). Naprawione.

### `\iffalse` block w sek08c działa poprawnie

Verification: 4 duplicate labels w body sek08c (`eq:metric-full-derived`,
`rem:step3-elimination`, `eq:vol-element-derived`,
`rem:metric-action-consistency`) wewnątrz `\iffalse...\fi` (lin. 389-622)
**są skipped przez TeX engine** — pdflatex nie zgłasza ich jako duplicates.
Mojej cleanup-introduced duplicates: **ZERO**.

## Cykl CORE-CLEANUP-B — ZAMKNIĘTY

**Wszystkie 6 faz wykonane:**

1. ✅ Phase 0: Setup cyklu (sesja 1)
2. ✅ Phase 1: Skeleton dodatek_pivot_history.tex + main.tex include (sesja 1)
3. ✅ Phase 2: sek08c forma III antypodalna → appendix (sesja 1) + formy I+II → appendix przez `\iffalse` (sesja 3)
4. ✅ Phase 3: sek08a v1.x cleanup — proof prop:kappa-corrected → appendix (sesja 2) + DEPRECATED markers (sesja 4)
5. ✅ Phase 4: Cross-ref audit ~30 plików P33 (sesja 5)
6. ✅ Phase 5/6: pdflatex compile clean check + POST_PATCH_VERIFICATION (sesja 6)

**Cleanup operations summary:**

- 6 sesji × ~2-3h każda = ~12-15h pracy (zgodne z pierwotną estymatą)
- 4 pliki rdzenia zmodyfikowane: `sek08a.tex`, `sek08c.tex`, `dodatek_pivot_history.tex` (NEW), `main.tex`
- 30+ deprecated labels relocated do appendix lub `\iffalse` block
- 6 NEW canonical labels dodane (M9.1''-spojne)
- pdflatex compile clean (547 stron, 0 nowych errors)

**Strategia sesji 3:** użycie `\iffalse...\fi` block do compile-bezpiecznego
wyłączenia starych labels w body sek08c (formy I+II w bloku 150 linii).
Labels (`thm:metric-from-substrate-full`, `eq:metric-full-derived`,
`rem:step3-elimination`, `eq:vol-element-derived`,
`rem:metric-action-consistency`) **są aktywne tylko w appendix** — body
zawiera kopie wewnątrz `\iffalse` (compile-skipped, brak duplicate-defined).

**Body cleanup post-sesja 3 (compile-effective):**

- sek08a body active: ~876 linii (155 deprecated content w komentarzach
  + canonical v2.0 ADDENDUM od lin. 677)
- sek08c body active: ~530 linii (compile-active, bez `\iffalse` block;
  76 lin NEW canonical M9.1'' theorem + 460 lin oryginalna struktura
  budżet+volume+forma III remarks)
- Appendix retencja: **641 linii historical record** z labels

**Body sek08c po `\iffalse` cleanup**:

- ~530 linii compile-active (NEW canonical M9.1'' + zachowane budżet substratowy
  prop:antipodal-from-budget, sssec:volume-element-check zastąpione
  rem:vol-element-M911, sssec:antipodal-uniqueness short remark formy III)
- 150 linii w `\iffalse...\fi` block (deprecated formy I+II, compile-skipped)
- Łącznie: 682 linii fizycznie w pliku, ~530 active, ~150 wycięte

**Body sek08a aktywne (post-sesja 1+2 cleanup):**

- v1.x section: zachowuje hyp:unified-action + ax:metric-coupling structure
  (zachowane, fundamental); proof prop:kappa-corrected i prop:field-eq-from-action
  zachowane jako short remark + canonical reference (V2.0 ADDENDUM canonical)
- v2.0 ADDENDUM: lin. 677-876 (kanoniczne G.0-closed, **CANONICAL**)

**Cross-references status:**

- ✅ 12+ labels w appendix: `thm:antipodal-uniqueness`,
  `rem:antipodal-implications`, `eq:antipodal-metric-ansatz`,
  `eq:antipodal-unique-solution`, `eq:anti-c1`, `eq:anti-c2`,
  `eq:anti-general-solution`, `eq:anti-gamma-of-s`, `eq:volume-FRW`,
  `eq:Lkin-FRW`, `eq:S-FRW-reduced`, `eq:ddt-pdotpsi`, `eq:dpsi-direct`,
  `eq:EL-FRW-full`, `eq:EL-FRW-simplified`, `eq:EL-FRW-divided`,
  `eq:psi-eq-unified`, `eq:cosmo-linearized-unified`, `eq:poisson-frw`,
  `eq:kappa-def-operational`
- ✅ Body sek08a zachowuje: `prop:field-eq-from-action`,
  `eq:field-eq-reproduced`, `eq:S-unified-psi`, `eq:EL-unified-raw`,
  `eq:EL-unified-divided`, `eq:S-static-correct`,
  `prop:kappa-corrected`, `eq:kappa-corrected` (BOXED VALUE INVARIANT),
  `prop:N07-resolved`, `eq:Gdot-corrected` (κ value invariant)
- ✅ NO duplicate labels (body vs appendix rozłączne sets)

**Status post-Sesja 1:**

| Klastr cleanup | Status |
|----------------|--------|
| `dodatek_pivot_history.tex` skeleton | ✅ DONE |
| sek08c forma (III) antypodalna | ✅ RELOCATED to appendix |
| sek08a `prop:kappa-corrected` proof (12 labels) | ✅ **RELOCATED to appendix** (sesja 1, faza 3.1) |
| sek08c forma (II) eksponencjalna | ✅ **RELOCATED to appendix** (sesja 3) |
| sek08c forma (I) potęgowa | ✅ **RELOCATED to appendix** (sesja 3, łącznie z formą II w `\iffalse` block) |
| sek08c NEW canonical M9.1'' theorem | ✅ **DODANE** (sesja 3): `thm:metric-from-budget-M911-canonical` + `eq:metric-M911-canonical` + `eq:vol-element-M911` + `rem:vol-element-M911` |
| sek08a `prop:field-eq-from-action` proof (eq:S-unified-psi, eq:EL-unified-*, eq:S-static-correct) | 🟡 PARTIAL — proof body zachowany z explicit STATUS marker; canonical replacement reference dodane; appendix summary stworzony. Phase 4 cross-ref audit zdecyduje czy literal removal body labels możliwe |
| sek08a `eq:V-selfinterference`, `eq:sqrt-g-eff`, `eq:g-eff-unified` w hyp:unified-action | ✅ **DEPRECATED MARKERS + canonical ref** (sesja 4): annotation-level cleanup. Body labels zachowane (5 \eqref użyć w body sek08a; 1 w op-phase1-covariant); explicit `[DEPRECATED]` flags inline; NEW `rem:hyp-unified-action-M911-canonical` z `eq:S-TGP-unified-M911-canonical` dodany |
| sek08a `prop:N07-resolved` (lin. 464+) — używa κ value | ✓ KEEP IN BODY (κ value invariant) |
| Cross-ref audit ~30 plików | ✅ **DONE sesja 5**: comprehensive globalny grep, klasyfikacja, compile-safety verification |
| pdflatex compile check | ⚠ NEXT SESSION (po pełnym cleanup) |

**Realistyczny czas pełnego Wariantu B:** 4-6 sesji × 2-3h każda
(podstawowa estymata 1-2 tygodni potwierdzona).

**Decyzja autora wymagana dla kolejnej sesji:** czy kontynuować
Phase 2 step 2 (sek08c forma II + I = przepisanie centralnego twierdzenia)
czy przeskoczyć do Phase 3 (sek08a v1.x deprecated, większa skala)?

## Cross-references

- [[README.md]] — plan cyklu
- [[../../audyt/CLEANUP_INVENTORY_2026-05-04.md]] — inwentaryzacja
- [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]
- [[../../core/formalizm/dodatek_pivot_history.tex]]
- [[../../main.tex]]
