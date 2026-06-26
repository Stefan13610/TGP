---
title: "Phase 2 — moduł A: F-BA-2 (GAP-2) DERIVED — bottom-up J_source z regulatorem marginalnościowym (Ṁ_bu = 2c³/9G EXACT; η = (t_*/t)² EXACT) · F-BA-3 (GAP-3) SUPPORTED_PARTIAL — dyspersja ścienna δ/(ct) EXACT (wykładnik 1), kanał odrzutu deklarowany"
type: phase_result
status: PHASE2_COMPLETE
phase: 2
cycle: op-frontier-bridge-and-asymmetry-2026-06-12
created_date: 2026-06-12
authorization: "User 2026-06-12: 'Phase 2'"
sympy_script: "[[./Phase2_sympy.py]]"
sympy_output: "[[./Phase2_sympy.txt]]"
sympy: "11/11 PASS; 0 hardcoded T_pass; circularity guard FP11"
falsifier_resolved: "F-BA-2 = DERIVED (regulator zidentyfikowany: marginalność/trychotomia na brzegu wejścia; cross-check EXACT; poziom ledger/hydrodynamiczny — kwalifikacja jawna) · F-BA-3 = SUPPORTED_PARTIAL (kanały ścienne DERIVED z wykładnikiem 1; kanał odrzutu kreacyjnego niewyprowadzony — deklarowany; atraktor pyłowy chroni machinery wzrostu)"
anti_lakatos_lock: PRESERVED
---

# Phase 2 — moduł A: bottom-up rate (GAP-2) + monochromatyczność A-iv (GAP-3)

## §0 — Verdict at a glance

| Falsyfikator | Werdykt | Esencja |
|---|---|---|
| **F-BA-2 (GAP-2)** | **DERIVED** (kwalifikacja §1.4) | Bottom-up funkcjonał J = ρ_e·(Ṙ − v_c) per area; **regulator = trychotomia marginalności zastosowana jako warunek brzegowy wejścia** (jedyne ρ_e = ρ̄; nad-depozycja → gałąź runaway, niedo-depozycja → gałąź blocked — obie wykluczone LOCKED); **Ṁ_bu = 2c³/9G EXACT == top-down**; wykładnik t⁰ vs capacity t²; η(t) = (t_*/t)² EXACT |
| **F-BA-3 (GAP-3)** | **SUPPORTED_PARTIAL** | Kanały kontrolowane ścianą **DERIVED**: dyspersja geometryczna i czasowa Δv/v_c = δ/(ct) EXACT (wykładnik 1 — warunek pre-rejestrowany spełniony dla tych kanałów; → 0 thin-wall); **kanał odrzutu kreacyjnego NIE wyprowadzony** (kinematyka kreacji nieznana — deklarowane); mitygacje DERIVED: v_pec ∝ 1/a_m, w_eff ∝ t^(−4/3) → 0 (atraktor pyłowy) |

Top-down Ṁ, v_c, marginalność, ρ̄ — użyte READ-ONLY (forbidden move #4 respektowany: zero modyfikacji).

## §1 — F-BA-2: bottom-up J_source z regulatorem (GAP-2)

### §1.1 Funkcjonał bottom-up (kinematyka kreacji brzegowej)

Kreacja jest boundary-localized (A-i LOCKED); materia zdeponowana na wewnętrznej krawędzi
porusza się z v_c = 2c/3, a krawędź ucieka z Ṙ = c ⇒ strumień masy pozostawianej w bulku:

```
J_source = ρ_e·(Ṙ − v_c)   [per area]   ⇒   Ṁ_bu = 4πR²·ρ_e·(c − 2c/3)
```

To jest brakująca w FM forma bottom-up: rate jako iloczyn gęstości depozycji × względnej
prędkości recesji krawędzi (FP1).

### §1.2 Regulator: trychotomia marginalności jako warunek brzegowy (FP2-FP3)

Pytanie GAP-2 brzmiało: czemu rate ∝ t⁰, skoro capacity ściany ∝ t²? Odpowiedź wyprowadzona:
**ρ_e nie jest swobodne** — księgowość mechaniczna na powierzchni wejścia musi być marginalna
(LOCKED zasada FCR P3, wyprowadzona tam trychotomią z M ∝ t i R = ct):

```
koszt(ρ_e) = ½v_c² − GM(ρ_e)/(ct),   M(ρ_e) = (4π/3)ρ_e(ct)³
koszt(ρ̄) = 0 EXACT;  ∂koszt/∂ρ_e < 0
ρ_e > ρ̄ ⇒ koszt < 0 (kreacja „z górki" → gałąź runaway — wykluczona przez R = ct)
ρ_e < ρ̄ ⇒ koszt > 0 (kreacja blokowana → gałąź blocked — wykluczona przez M ∝ t)
⇒ ρ_e = 1/(6πGt²) JEDYNE (FP2: jedyny pierwiastek)
```

**Regulator = marginalność czytana lokalnie na brzegu:** ściana konwertuje dokładnie tyle
offsetu ΔV, ile marginalne wchłanianie dopuszcza; rate jest demand-limited, nie
supply-limited. Wstawienie ρ_e do §1.1: **Ṁ_bu = 2c³/9G EXACT** (FP1) — cross-check
bottom-up↔top-down domknięty; wykładniki: rate t⁰, capacity t² (FP5 audit).

### §1.3 Efektywność konwersji i dyspozycja nadwyżki (FP4)

η(t) = popyt/podaż = (Ṁc²)/(ΔV·4πR²c) = **(t_*/t)² EXACT**, z t_* dokładnie w formie FM P3
LOCKED (tożsamość symboliczna — cross-check niezależny). η(t_*) = 1; nadwyżka 1−η → 1 idzie
w kinetykę ściany (γv rośnie, v → c — spójność z FM P1 FP5/R-FM-3, read-only). Epoka t < t_*:
reżim supply-limited — deficyt wczesnoepokowy odziedziczony bez zmian (INFORMATIONAL FM;
moduł B SIG-1 przesuwa próg do √2·t_* przy kreacji parowej).

### §1.4 Kwalifikacja werdyktu (uczciwa, jawna)

- **Poziom wyprowadzenia: ledger/hydrodynamiczny (consistency-closure).** Łańcuch
  Ṁ→ρ̄→Ṁ przechodzi przez konfigurację LOCKED (mapa powłok FM COR-1) — charakter punktu
  stałego nazwany wprost (precedens: FM COR-1 DERIVED_SELF_CONSISTENT). Niezależna treść
  bottom-up: (i) identyfikacja regulatora (trychotomia jako warunek brzegowy — FP3 z jawnymi
  kierunkami wykluczeń), (ii) forma funkcjonału J (FP1), (iii) η = (t_*/t)² łącząca ledger
  ściany z popytem (FP4), (iv) audit wykładników (FP5).
- **Czego NIE wyprowadzono (poza zakres pre-rejestracji, flagowane):** statystyczny rate
  nukleacji z fluktuacji wokół profilu ściany (prefaktor klasy Langera/Colemana — wymaga
  machinery fluktuacyjnej, Appendix E territory; teren §10.1). Pre-rejestrowana klasa DERIVED
  żądała: funkcjonał + regulator + zgodność wykładnika — wszystkie trzy spełnione; brak
  prefaktora statystycznego odnotowany jako kandydat przyszłej rafinacji, NIE ukryty.

**F-BA-2 = DERIVED.**

## §2 — F-BA-3: A-iv z dynamiki ściany (GAP-3)

### §2.1 Kanały kontrolowane ścianą — DERIVED (FP6-FP8)

- **Powierzchnia wejścia jest JEDYNA i ostra (FP6):** profil monotoniczny ⇒ pojedyncze
  przecięcie poziomicy |Φ| = Φ₀/√3; głębokość x* = δ·atanh(1/√3), m_Φx* = ln(2+√3) EXACT
  (zgodność z LOCKED FM P1 FP6) — zdarzenie wejścia jest dobrze określoną powierzchnią,
  nie pasmem.
- **Dyspersja geometryczna (FP7):** depozycja w podwarstwie grubości ≤ δ; flow-matching
  u = (2/3)x/t ⇒ **Δv/v_c = δ/(ct) EXACT — wykładnik 1 w (δ/ct)** (warunek pre-rejestrowany
  „≥ 1" spełniony); → 0 w limicie thin-wall.
- **Dyspersja czasowa (FP8):** czas formacji δ/c × |∂u/∂t| ⇒ ten sam rząd δ/(ct), wykładnik 1.

### §2.2 Kanał odrzutu kreacyjnego — NIE wyprowadzony (deklarowane) + mitygacje DERIVED

Rozkład pędów pary przy kreacji wymaga kinematyki zdarzenia kreacyjnego — niedostępnej
w statycznej machinery (ta sama deklaracja co Phase 1 §5.3). Bez tego **pełne** „zerowa
dyspersja wejścia" nie jest twierdzeniem wyprowadzonym. Mitygacje wyprowadzone:

- **Zanik dyspersji (FP9):** dv_pec/dt = −(∂u/∂x)v_pec ⇒ v_pec ∝ t^(−2/3) = 1/a_m EXACT
  (zgodne z LOCKED wykładnikami trajektorii FM P2 FP7 {2/3, 1/3}) — każdy odrzut początkowy
  jest transientem.
- **Atraktor pyłowy (FP10):** w_eff = ⟨v_pec²⟩/c² ∝ t^(−4/3) → 0 — człon ciśnieniowy
  w machinery wzrostu zanika szybciej niż człony pyłowe ⇒ **C-DERIVED form (LOCKED) jest
  chroniona asymptotycznie nawet przy niezerowej dyspersji wejścia**. Rola A-iv w K2/FM jest
  więc strukturalnie zabezpieczona, choć samo A-iv przy wejściu pozostaje częściowo otwarte.
- INFORMATIONAL (nie użyte do oceny): blisko-krytyczność wyścigu z Phase 1 sugeruje, że
  populacja posortowana ucieka z energią resztkową ≈ 0 (escape na progu) — kierunek
  „mały odrzut naturalny"; bez statusu dowodowego.

**F-BA-3 = SUPPORTED_PARTIAL** (kanały ścienne DERIVED; kanał odrzutu = luka deklarowana;
machinery wzrostu chroniona atraktorem). Ocena strict — bez łagodzenia: warunek
pre-rejestrowany dotyczył pełnej dyspersji wejścia, a kanał odrzutu nie jest kontrolowany
parametrem δ/ct.

## §3 — Status warunków PR-022 / F-BA-D po Phase 2

| GAP | Status |
|---|---|
| GAP-1 (EQ-5 field-level) | OPEN → Phase 3 |
| **GAP-2 (bottom-up J_source)** | **DERIVED** (§1) |
| **GAP-3 (A-iv)** | **SUPPORTED_PARTIAL** (§2) |
| GAP-4 (uniqueness) | OPEN → Phase 3 |
| GAP-5 (sphericity) | OPEN → Phase 3 |

**Uczciwa flaga progu (już teraz, nie w FINAL):** F-BA-D BRIDGE_COMPLETE wymaga F-BA-1..5
wszystkich DERIVED(-IN-CLASS). SUPPORTED_PARTIAL na GAP-3 w ścisłym odczytaniu **blokuje
PR-022 append** — chyba że użytkownik w Phase FINAL uzna próg za spełniony przy atraktorze
pyłowym (analogia: decyzja progowa FM (iv) PARTIAL — wtedy odmówiono obniżenia poprzeczki).
Decyzja należy wyłącznie do użytkownika; tu tylko flagujemy bez rekomendacji łagodzącej.
**NO PR-022 w tej fazie** (forbidden move #8).

## §4 — Anti-Lakatos (Phase 2): COMPLIANT ✓

0/11 hardcoded; top-down Ṁ/v_c/marginalność/ρ̄ read-only (forbidden move #4 ✓ — bottom-up
wyłącznie cross-check); F-BA-3 oceniony SUPPORTED_PARTIAL wbrew pokusie DERIVED (kanał
odrzutu nie został zamieciony pod atraktor); kwalifikacja poziomu F-BA-2 jawna (consistency-
closure; prefaktor statystyczny flagowany); naprawa plumbingu FP6 (exp-rewrite tożsamości
atanh/log — klasa identyczna jak Phase 1 §5.6) bez zmiany żadnego progu; circularity guard
FP11 czysty; λ, Φ₀ symboliczne; 0 nowych stałych; 0 predecessor verdicts modified.
