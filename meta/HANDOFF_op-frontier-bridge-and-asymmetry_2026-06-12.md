---
title: "HANDOFF — prompt dla agenta wykonującego cykl op-frontier-bridge-and-asymmetry"
date: 2026-06-12
type: meta-handoff
status: READY (cykl REGISTERED-QUEUED; aktywacja = user 'działaj')
cycle: "[[../research/op-frontier-bridge-and-asymmetry-2026-06-12/]]"
authored_by: "Claudian, sesja #16 FINAL post-QA (user: 'zarejestruj przyszły cykl i rozpisz prompt dla nowego agenta')"
---

# HANDOFF PROMPT — skopiuj poniższy blok jako zadanie nowego agenta

---

Jesteś ekspertem fizyki teoretycznej pracującym w ramach TGP (Teoria Generowanej Przestrzeni,
projekt `TGP/TGP_v1/`). Twoim zadaniem jest wykonanie cyklu badawczego
**`op-frontier-bridge-and-asymmetry-2026-06-12`** (REGISTERED-QUEUED → aktywacja po Phase 0
+ autoryzacji użytkownika).

## 0. Tożsamość metodologiczna (nienegocjowalne)

Pracujesz w reżimie anty-Lakatosowskim: pre-rejestracja PRZED rachunkiem; falsyfikatory
z progami liczbowymi LOCKED w Phase 0 i nigdy nie korygowane ex post; HONEST_NEGATIVE jest
pełnoprawnym wynikiem; każda faza wymaga osobnej autoryzacji użytkownika („działaj");
0 hardcoded T_pass w sympy; 0 nowych stałych fundamentalnych bez deklaracji; werdykty
poprzedników LOCKED — nigdy nie modyfikowane wstecznie.

## 1. Obowiązkowe lektury (W TEJ KOLEJNOŚCI, przed Phase 0)

1. `TGP/TGP_v1/STATE.md` — wpisy sesji #15 (FCR) i #16 (FM) — stan łańcucha
2. `research/op-frontier-microphysics-2026-06-11/` — CAŁY cykl FM: Phase0_balance (wzorzec
   pre-rejestracji), Phase1-3_derivation + sympy, Phase_FINAL_close (**§4 GAP REGISTER —
   to jest Twoja lista zadań**)
3. `research/op-frontier-creation-rate-derivation-2026-06-11/` — FCR (marginalność, derived
   flow, C-DERIVED form)
4. `meta/SCOPING_op-frontier-bridge-and-asymmetry_2026-06-12.md` — analiza Q1-Q3 + hipotezy
   robocze modułu B (H-SORT/H-CP) i ich sygnatury falsyfikowalne
5. `meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md` — concept paper (ontologia §1-§2,
   EQ-1..6 §4, falsyfikatory F4-F9 §7, ryzyka §10)
6. `meta/CALIBRATION_PROTOCOL.md` **§3.6 + rozszerzenia §3.6.6-§3.6.13** — WIĄŻĄCE
   (analytical pre-derivation, sign conventions, assumptions enumeration, precision ±5%,
   klasyfikacja stałych)
7. `meta/CYCLE_KICKOFF_TEMPLATE.md` §1-§3 — kontrakt kickoff
8. Dla modułu B: `research/op-CE-H-two-particle-equilibrium-2026-05-21/` (kink–antykink:
   anihilacja w izolacji LOCKED; V_int ∝ exp(−m√2·L)) + FFS pre-screening (antymateria =
   operacja C: Φ → Φ*)

## 2. Cel cyklu

**Moduł A (bridge completion):** domknąć GAP-1..5 z FM Phase_FINAL §4:
- GAP-1: most polowy EQ-5 (S_Φ ↔ S_matter w jednostkach pola; w koncepcie schematyczny)
- GAP-2: bottom-up rate nukleacji J_source z mikrofizyki ściany (cross-check top-down
  Ṁ = 2c³/9G z marginalności — top-down NIE wolno modyfikować)
- GAP-3: A-iv (monochromatyczność wejścia) z dynamiki ściany
- GAP-4: uniqueness samouzgodnionego rozwiązania jednorodnego (existence+exactness już są)
- GAP-5: sferyczność locusu frontu
**Sukces modułu A ⇒ PR-022 append:** predykcja jednopunktowa bezparametrowa log₁₀G = 2.025
(p = 2/3 EdS EXACT; v_c = 2c/3) vs observed 3.0 — z OBOWIĄZKOWĄ honest-physics note (0.97 dex
poniżej; PASS_BAND krawędziowo 0.025 dex; dyskusja rozbieżności otwarta, nie rescue).

**Moduł B (creation selectivity, GAP-6):** czy ściana rozróżnia Φ od Φ*? Zbiór hipotez
CLOSED do pre-rejestracji w Phase 0:
- **H-SORT:** amplituda C-symetryczna (pary) + sortowanie orientacją ściany (soliton zgodny
  z gradientem → bulk; C-partner → sektor frontowy, ukryty za marginalnie nieosiągalnym
  frontem). Test mikrofizyczny: energia osadzenia kinku vs antykinku w tle gradientowym
  profilu ściany (machinery 1D z op-CE-H — reuse).
- **H-CP:** asymetria w samej amplitudzie kreacji (sprzężenie fazy U(1)/RP² z orientacją
  ściany; Sakharov-analog: out-of-eq ✓ ściana, B-violation ✓ konwersja, C/CP = do wyprowadzenia).
- **HONEST_NEGATIVE:** ściana nie rozróżnia ⇒ GAP-6 eskaluje do falsyfikatora ramy CE-H
  (kreacja frontowa netto-materii bez mechanizmu) — zapisz bez łagodzenia.

Pre-rejestrowalne sygnatury ilościowe H-SORT (Phase 0 ma je przekuć w falsyfikatory z progami):
(1) budżet popytu ×2 ⇒ t_* → √2·t_*; (2) frakcja leakage antysolitonów do bulku ⇒ przewidywane
rozproszone tło γ (porównanie z obserwacjami = comparison-only); (3) wyścig sortowanie vs
anihilacja: siła sortująca musi rozdzielić parę szybciej niż anihilacja (szybkość z LOCKED
V_int ∝ exp(−m√2·L)) — warunek ilościowy istnienia mechanizmu.

## 3. Phase 0 — wymagania twarde (przed jakimkolwiek rachunkiem)

1. Falsyfikatory F-BA-* per moduł, z klasami werdyktów (DERIVED/PARTIAL/GAP; dla B:
   H-SORT_DERIVED / H-CP_DERIVED / TWO_HYPOTHESES_PERSIST / HONEST_NEGATIVE) — CLOSED sets
2. Analytical pre-derivation (§3.6.1/.6/.8/.9): wartości oczekiwane policzalne symbolicznie
   PRZED sympy; konwencje znaków; enumeracja założeń per falsyfikator
3. Bands: log₁₀G_obs = 3.0 comparison-only (inherited LOCKED); dla modułu B zadeklaruj
   targety porównawcze (np. ograniczenia na frakcję antymaterii / tło γ) jako
   OBSERVATIONAL_ANCHOR comparison-only — **η_B ani żadna obserwowana asymetria NIGDY nie
   jest inputem wyprowadzenia** (to analog zakazu G_obs — circularity guard FP w każdym skrypcie)
4. Forbidden moves (min. 10; odziedzicz z FM Phase 0 §3 + dodaj: zakaz modyfikacji top-down
   rate; zakaz selekcji hipotezy B wynikiem; zakaz „ukrywania" antymaterii bez ilościowego
   warunku 3; zakaz fizyki ściany spoza TGP Lagrangianu — luka = GAP, nie improwizacja)
5. Stałe: budżet nowych fundamentalnych = 0; λ, Φ₀ symboliczne; każdy anchor zadeklarowany
6. Risk register (min.: §10.1 calculational hell; wariant „hemisfery przestrzenne" ODRZUCONY
   na starcie — Zel'dovich + izotropia CMB, odnotuj w §forbidden/risk; ryzyko
   niefalsyfikowalności H-SORT bez sygnatur 1-3 — dlatego są obowiązkowe)
7. Anticipated outcomes (INFORMATIONAL) — zapisz kierunek ciągnięcia kryteriów PRZED rachunkiem

## 4. Protokół wykonania

- Czekaj na user „działaj" przed Phase 0 LOCK i przed każdą kolejną fazą; po każdej fazie
  raportuj z decision menu
- Każda faza: dokument PhaseN_derivation.md + PhaseN_sympy.py/.txt (PASS/FAIL per FP,
  0 hardcoded, circularity-guard FP obowiązkowy)
- Kolejność rekomendowana: Phase 0 → moduł B test mikrofizyczny 1D (kink w tle gradientowym —
  najtańszy decydujący rachunek, reuse op-CE-H) → moduł A GAP-2/GAP-3 (mikrofizyka ściany,
  wspólna machinery) → GAP-1/GAP-4/GAP-5 → Phase FINAL (agregat; decyzje PR — wyłącznie user)
- STATE.md: wpis po każdej fazie (wzorzec: sesje #15-#16); README cyklu: flip parking → active
  przy aktywacji
- Multi-sesja dozwolona; jeśli rachunek nieprzejezdny (§10.1) — HONEST_NEGATIVE/INDETERMINATE
  zamiast przybliżeń bez deklaracji; każde przybliżenie (thin-wall, mean-field, 1D-proxy)
  jawnie zadeklarowane per-use z dyspozycją „1D-proxy ≠ 3D claim"

## 5. Czego absolutnie nie wolno

- Modyfikować LOCKED: FM/FCR/R17/P25/γ-*/F8/PR-017..020/CE-H two-particle
- Używać G_obs, η_B_obs, jakiejkolwiek obserwowanej asymetrii jako inputu wyprowadzeń
- Appendować PR-022 bez domknięcia GAP-1..5 ORAZ bez decyzji użytkownika
- Cytować archiwalnych ex263/ex279 (leptogeneza pre-restart) jako LIVE źródła — wolno
  wyłącznie jako kontekst historyczny z etykietą NIE-LIVE
- Zamykać modułu B „na miękko": werdykt musi być jednym z CLOSED set

---

**Koniec promptu.** Stan rejestracji: cykl REGISTERED-QUEUED (folder + README istnieją);
aktywacja wymaga user „działaj" → Phase 0 LOCK.
