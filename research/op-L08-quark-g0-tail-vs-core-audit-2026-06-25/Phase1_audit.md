---
title: "Phase 1 FAST-AUDIT — werdykt NORM-OVERLOAD (9/9 FP PASS, wyliczony): zakres sek08b:529 [0,817;0,891] to PASMO BAZOWE g₀^e (≈ g₀ elektronu / φ-FP anchor; pokrywa się z pasmem g_min ogona), NIE pełna domena rdzeniowego g₀; sufit HALT-B T11=2,68× VOIDED przez wewnętrzną sprzeczność z hierarchią leptonową"
type: phase1_results
status: COMPLETE
date: 2026-06-25
cycle: op-L08-quark-g0-tail-vs-core-audit-2026-06-25
parent: "[[./README.md]]"
sympy_pass: "9/9 FP"
verdict: NORM-OVERLOAD
hardcoded: 0
---

# Phase 1 FAST-AUDIT — wyniki

> **Werdykt (WYLICZONY z flag T1–T8, NIE wybrany): NORM-OVERLOAD.**
> Zakres `g₀ ∈ [0,817; 0,891]` (sek08b:529) **NIE jest** pełną domeną rdzeniowego
> sprzężenia g₀ kwarków, jak założył sufit strukturalny T11=2,68× cyklu HALT-B.
> To **wąskie pasmo wokół wartości bazowej g₀^e ≈ 0,8694** (elektron / kotwica φ-FP),
> które dodatkowo **pokrywa się z pasmem minimów ogona** g_min∈[0,742; 0,898]. Sufit
> T11 zostaje **VOIDED** (nie sfalsyfikowany — *nieważny dla rdzeniowego g₀*). HALT-B
> → **reopened**; licencja na osobny cykl rescue-test.

## §1 — Co rozstrzygnęła ekstrakcja (dodatekJ_ogon_masy.tex)

### §1.1 — Trzy rozłączne wielkości „g" (T1 PASS)

`dodatekJ` definiuje jednoznacznie (eq:J-ode, eq:J-tail):

| Symbol | Definicja | Rola |
|---|---|---|
| **g₀** | warunek brzegowy `g(0)=g₀` (parametr „strzału") | **WEJŚCIE ODE** — sprzężenie rdzeniowe definiujące soliton |
| **g(r)** | rozwiązanie profilu | wynik |
| **g_min** | ekstremum oscylacyjnego ogona (eq:J-tail) | **WARTOŚĆ PROFILU** mierzona w polu dalekim (zewnętrzna) |

To są **różne kategorie** — wprost rozróżnienie ogon/rdzeń z uwagi użytkownika
(g_min = wartość mierzona z zewnątrz, ≠ g₀ rdzeniowe „w układzie własnym").

### §1.2 — Rdzeniowe g₀ leptonów (ODE substratowe, kanoniczne; tab. l.435–449, ex157)

```
e:   g₀ = 0,869470   (baza)
μ:   g₀ = 1,406833   (= φ·g₀^e)
τ:   g₀ = 1,729615
```
⇒ leptony rozpinają **rdzeniowe g₀ ∈ [0,869; 1,730]** (rozstęp ε_core ≈ 0,66; ratio ≈ 1,99×).

### §1.3 — Lokalizacja przedziału [0,817; 0,891] (T2, T3 PASS)

- **T2:** [0,817; 0,891] zawiera **tylko** g₀ elektronu (0,8694); **wyklucza** mion (1,407)
  i tauon (1,729). Czyli NIE jest to domena, na której żyje hierarchia — to pasmo **bazy**.
- **T3:** [0,817; 0,891] **leży wewnątrz** pasma minimów ogona g_min∈[0,742; 0,898]
  (μ:0,898, τ:0,742). Numerycznie pasmo bazowe g₀^e i pasmo g_min **zbiegają się** w tym
  samym oknie <1 (bo g₀^e≈0,87 i g_min leptonów oscylują wokół ~0,8) — stąd przeciążenie.

## §2 — Dlaczego sufit T11 jest VOIDED (T5, T6 — rdzeń werdyktu)

### §2.1 — Rekonstrukcja sufitu (T5 PASS)

Argument HALT-B: `max(m)/min(m) ≤ (A_max/A_min)²·(1+ε)^(e²/2)`, z ε wziętym z [0,817;0,891]:
- ε_I = 0,0867; `(1+ε_I)^(e²/2) = 1,359`; `(A_max/A_min)² ≈ 1,93` ⇒ **ceiling ≈ 2,62** (≈ raportowane 2,68). Rekonstrukcja potwierdzona.

### §2.2 — Test wewnętrznej spójności = nóż (T6 PASS, CENTRALNY)

**Ten sam mechanizm `m ∝ A_tail⁴` zastosowany do LEPTONÓW:**
- leptony używają rdzeniowego g₀ ∈ [0,869; 1,730] (NIE mieszczą się w [0,817; 0,891]);
- dają **m_τ/m_e = 3477** (odtworzone mechanizmem do 0,006%, dodatekJ ex157).

`3477 ≫ 2,62` (o czynnik ~1327×). Zatem **gdyby sufit 2,68× był ważnym ograniczeniem
strukturalnym mechanizmu, zakazywałby też ustalonej hierarchii leptonowej** — co jest
absurdem (leptony są fundamentem kalibracyjnym why_n3, LIVE). Wniosek wyliczony:

> Sufit 2,68× nie ogranicza mechanizmu — ogranicza jedynie **sztucznie zawężoną domenę
> [0,817; 0,891]**, która nie zawiera nawet mionowego/tauonowego rdzeniowego g₀.
> Cykl HALT-B przetestował **złą zmienną** (pasmo bazowe/ogona zamiast dynamicznej domeny
> rdzeniowego g₀). ⇒ **NORM-OVERLOAD**.

### §2.3 — Skala błędu kategorii (T7 PASS; INFORMATIONAL, NIE rescue)

Pod właściwą domeną rdzeniową leptonów ε_core ≈ 0,66 (>5× ε_I); `(1+ε_core)^(e²/2) ≈ 6,5`.
**To NIE jest dowód reprodukowalności kwarków** (R4/forbidden #10) — pokazuje wyłącznie, że
domena użyta w T11 była o rząd wielkości za wąska. Rekalkulacja stosunków = osobny rescue-test.

## §3 — Circularity guard (T8 PASS)

Werdykt kategoryzacji użył wyłącznie: rdzeniowych g₀ leptonów (ODE), g_min ogona, przedziału
[0,817;0,891] oraz **leptonowego** anchora m_τ/m_e. **Żadnej masy KWARKA PDG** — kategoria g
NIE została wybrana przez to, „co ratuje" 80000×.

## §4 — Werdykt + flagi

| Test | Pytanie | Wynik |
|---|---|---|
| T1 | g₀ (wejście) ≠ g_min (profil)? | PASS |
| T2 | I wyklucza μ/τ rdzeniowe, zawiera tylko e? | PASS |
| T3 | I ⊂ pasmo g_min ogona [0,742;0,898]? | PASS |
| T5 | rekonstrukcja sufitu T11 ≈ 2,68×? | PASS |
| T6 | **werdykt: leptony łamią sufit ⇒ NORM-OVERLOAD** | PASS (wyliczony) |
| T7 | ε_core ≫ ε_I (skala błędu)? | PASS |
| T8 | brak masy kwarka PDG w kategoryzacji? | PASS |

**9/9 FP PASS (T1,T2×3,T3,T5,T6,T7,T8); T9-DEC.** 0 hardcoded. Werdykt **NORM-OVERLOAD**.

## §5 — DOUBTS register

- **D1 (drugorzędna niespójność wzoru — POZA zakresem, do flagowania):** cykl HALT-B użył
  `m = c_M·A_tail²·g₀^(e²/2)`, podczas gdy `dodatekJ` (eq:J-mass-Atail4-analytic) ma
  `M = c_M·A_tail⁴`. To dwa różne wzory (A_tail² vs A_tail⁴ + dodatkowy g₀^3,69). Należy
  rozstrzygnąć w cyklu rescue-test, który wzór jest kanoniczny — wpływa na rachunek stosunków.
- **D2:** g₀^e w trzech wariantach ODE: 0,8339 / 0,8694 / 0,8993 — [0,817;0,891] to mniej
  więcej to rozmycie wariantu obliczeniowego wokół bazy, NIE okno fizyczne kwarków.
- **D3:** NORM-OVERLOAD zdejmuje błędne no-go; rescue-test może dać HALT-B z INNEGO powodu
  (np. bariera ducha w Formulacji B obcina spektrum przy g₀≳1,63 — l.398–410; w Formulacji A
  substratowej brak bariery, ale spektrum kwarkowe wymaga g₀ rozpiętego znacznie szerzej niż
  leptonowe 1,73 → realne pytanie rescue-testu).

## §6 — Dyspozycja (do Phase FINAL)

NORM-OVERLOAD ⇒ (per reguła §0.2 LOCKED):
1. **Sufit T11 = 2,68× VOIDED** (nieważny dla rdzeniowego g₀).
2. **HALT-B reopened** — status: `STRUCTURAL_INSUFFICIENCY` → `INDETERMINATE-PENDING-RESCUE`
   (werdykt poprzednika NIETYKALNY jako test hipotezy „I = domena core g₀", która okazała się
   błędną kategorią).
3. **Licencja** na osobny cykl `op-quark-mass-core-g0-rescue-test` (nowy PR, własny Phase 0).
4. **Korekta sek08b:529** (housekeeping, user-gated): „g₀ ∈ [0,817;0,891]" myląco sugeruje
   uniwersalną domenę; faktycznie to pasmo bazowe g₀^e — przeformułować na „baza g₀^e ≈ 0,869
   (±wariant ODE); hierarchia generowana przez rdzeniowe g₀ ∈ [g₀^e, g₀^τ]".
5. **PR-014 NIE formalizowany** (był warunkowy na NORM-COHERENT).
