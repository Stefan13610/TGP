---
title: "Phase 3 — target #1: frontier marginality condition — PRINCIPLE_DERIVED + COEFFICIENT_TWO_POINT"
type: phase_result
status: PHASE3_COMPLETE
phase: 3
cycle: op-frontier-creation-rate-derivation-2026-06-11
created_date: 2026-06-11
authorization: "User 2026-06-11: 'ok wszystko jasne działaj' (po klaryfikacji semantyki 'zero-energy' → marginalność TGP-natywna)"
sympy_script: "[[./Phase3_sympy.py]]"
sympy_output: "[[./Phase3_sympy.txt]]"
sympy: "7/7 PASS; 0 hardcoded T_pass=True"
falsifier_resolved: "F-ZE = PRINCIPLE_DERIVED + COEFFICIENT_TWO_POINT {ε ∈ {2/3, 3/2}} + TIEBREAKER_OPEN; B1 upgrade STRUCTURAL_POSTULATE → MARGINALITY-DERIVED (two-point); F-FCR-D STRUCTURAL_CONDITIONAL (SHARPENED); NO PR-022"
anti_lakatos_lock: PRESERVED
---

# Phase 3 — frontier marginality condition (dawne „zero-energy" B1)

## §0 — Verdict at a glance

**F-ZE = PRINCIPLE_DERIVED + COEFFICIENT_TWO_POINT + TIEBREAKER_OPEN.**

Kluczowa zmiana statusu: postulat B1 („zero-energy") przestaje być arbitralny — jest
**konsekwencją stabilności stacjonarnej kreacji frontowej**, z dokładnym dwupunktowym zbiorem
współczynnika. Semantyka TGP-natywna (klaryfikacja użytkownika): „zero" = **zerowy koszt NETTO
mechaniczny** kreacji względem nasyconej próżni E2 (substrat = punkt odniesienia, nie nieobecność).

## §1 — Zasada (FP1): trychotomia stabilności wymusza marginalność

Stacjonarny front (γ-3: R = ct, LOCKED) + stacjonarna kreacja (M ∝ t, Phase 1 skeleton):
- koszt netto > 0 → kreacja zablokowana (jak w bulku, E2 saturation) → M ∝ t sprzeczne ✗
- koszt netto < 0 → kreacja lawinowa → R = ct steady sprzeczne ✗
- **koszt netto = 0 → marginalność → jedyna spójna opcja ✓**

Warunek marginalnego związania elementu kreowanego na froncie: (1/2)v_c² = GM/(ct) —
**dokładnie ta sama semantyka, którą standardowa kosmologia definiuje gęstość krytyczną**.
Unikalne rozwiązanie M(t) ∝ t (odtwarza szkielet Phase 1). 

## §2 — Współczynnik (FP2-FP3): ξ = v_c²/(2c²) EXACT; filtry zbioru bookkeepingów

| Bookkeeping | ξ | Status |
|---|---|---|
| B-k1 rest-energy (c²dM) | 1 | **EXCLUDED** — koszt spoczynkowy to konwersja wewnątrz-substratowa (sama kreacja), nie wiązanie mechaniczne; poza zakresem zasady |
| B-k2 uniform-sphere total (3/5·GM²/R) | 5/3 | **EXCLUDED** — warunek globalny, nie marginalny-na-froncie |
| **B-k3 v_c = c** (prędkość frontu, γ-3) | **1/2** | SURVIVES — **identyfikacja: to jest dokładnie warunek Schwarzschilda R = 2GM/c²** (stąd nazwa „zero-energy" w Phase 1) |
| **B-k4 v_c = 2c/3** (wyprowadzony przepływ Phase 2 przy x = ct) | **2/9** | SURVIVES — spójność kinematyczna z derived flow |

**ε_G = (3/2)(v_c/c)² EXACT ⇒ ε ∈ {3/2, 2/3}** — dwupunktowy zbiór dokładny, zero strojenia.

## §3 — Predykcja dwupunktowa (FP4): OBA punkty w paśmie PASS

| Kandydat | p (EXACT) | log₁₀G | Verdict |
|---|---|---|---|
| B-k3 (ε = 3/2) | (√55−1)/6 = 1.0694 | **3.249** | PASS_BAND (0.25 dex od obs.) |
| B-k4 (ε = 2/3) | **2/3 (= EdS dokładnie!)** | **2.025** | PASS_BAND (na krawędzi: 0.025 dex — uczciwa nota) |

Obserwowane log₁₀G = 3.0 leży **między** punktami, bliżej B-k3. Ciekawostka strukturalna: B-k4
odtwarza wykładnik wzrostu Einsteina–de Sittera z zupełnie innej ontologii (kreacja frontowa
zamiast materii dominującej w FRW).

## §4 — Tiebreaker (FP6): OPEN — i to jest uczciwa granica cyklu

Rozstrzygnięcie B-k3 vs B-k4 = pytanie „z jaką prędkością wchodzi materia kreowana na froncie?"
= mikrofizyka frontu = **concept paper §10.6 Q4 (»czym dokładnie jest frontier?«) — OPEN**.
Żaden element LOCKED nie rozstrzyga; cherry-pick zakazany. TIEBREAKER_OPEN zadeklarowany.

## §5 — Aggregate update (mechanical): STRUCTURAL_CONDITIONAL (SHARPENED)

- B1: STRUCTURAL_POSTULATE → **MARGINALITY-DERIVED (two-point)** ⭐
- PR-022: **WITHHELD** (tiebreaker + A-ii/C-2/A2 caveats) — ale kandydat-PR jest już formułowalny
  jako **dwupunktowa predykcja bezparametrowa**: log₁₀G_TGP ∈ {2.025, 3.249} vs obs. 3.0
- Łańcuch braków do PR-022 (zaktualizowany): (1) ~~zero-energy~~ → **tiebreaker mikrofizyki frontu
  (§10.6 Q4)** — JEDYNY brak główny; (2) ~~bulk transport~~ DONE; korolaria: A-ii, C-2, A2

## §6 — Anti-Lakatos (Phase 3): COMPLIANT ✓

0/7 hardcoded; zbiór bookkeepingów CLOSED pre-declared z filtrami zasadniczymi (nie wynikowymi —
B-k1/B-k2 wykluczone semantyką zasady, nie wartością wzrostu; FP6 guard: G_obs nieobecne);
oba surviving punkty raportowane (bez selekcji); krawędziowość B-k4 (0.025 dex) ujawniona;
0 predecessor verdicts modified; doc-cleanup note: rename „zero-energy" → „frontier marginality
condition" w przyszłych dokumentach.
