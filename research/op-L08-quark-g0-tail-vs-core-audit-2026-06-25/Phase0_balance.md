---
title: "Phase 0 — pre-registration LOCKED: op-L08-quark-g0-tail-vs-core-audit — kategoria wielkości g∈[0,817;0,891] (rdzeń g₀ vs ogon g_min<1) determinująca ważność sufitu HALT-B 2,68×; falsyfikatory NORM-OVERLOAD/COHERENT/INDETERMINATE LOCKED"
type: phase0_balance
status: LOCKED
locked_date: 2026-06-25
cycle: op-L08-quark-g0-tail-vs-core-audit-2026-06-25
authorization: "User 2026-06-25: 'ok działaj z tym, rozpisz nowy cykl' + uwaga merytoryczna: 'w ramach TGP solitony (więc możliwe że kwarki też) mają strukturę wewnętrzną i ogon; ogon to wartość mierzalna z zewnątrz, różna od wartości mierzonych we własnym układzie odniesienia kwarka'. Fraza aktywacyjna pokrywa Phase 0 LOCK (precedens domu #15/#16/#17/#19); Phase 1+ = osobne 'działaj'."
methodology_binding: "CALIBRATION_PROTOCOL §3.6 (analytical pre-derivation + 8/8 gate) BINDING; CYCLE_KICKOFF_TEMPLATE §1-§3 (L1 observable-first); werdykt HALT-B op-L08-Phase6-quark-sector-mass-formula IMMUTABLE/NIETYKALNY (ten cykl bada ZAKRES WAŻNOŚCI argumentu T11, NIE jego poprawność dla testowanej hipotezy); PRE_REGISTERED_FALSIFIERS §0.3 append-only"
anti_lakatos_lock: PRESERVED
PR_candidate: "PR-014 (formalizacja TYLKO przy NORM-COHERENT; LOCK rejestru = user w Phase FINAL). Rescue-test = osobny przyszły PR (NIE rezerwowany tu)."
queue_next_after_cycle: "Zależne od werdyktu: NORM-OVERLOAD ⇒ user-decyzja o op-quark-mass-core-g0-rescue-test; NORM-COHERENT ⇒ housekeeping demotacji; NORM-INDETERMINATE ⇒ housekeeping definicji app:ogon-masy. Brak scope-creep w tej sesji."
---

# Phase 0 — pre-registration LOCKED (przed jakąkolwiek ekstrakcją/rachunkiem)

## §0 — Tożsamość cyklu (nienegocjowalne)

1. **AUDYT EPISTEMICZNY, NIE RATUNEK.** Werdykt HALT-B cyklu
   `op-L08-Phase6-quark-sector-mass-formula-2026-05-16` (0/5 stosunków, sufit T11=2,68×
   vs wymagane 80000×) jest **IMMUTABLE**. Pozostaje na zawsze poprawny **dla hipotezy,
   którą testował**: „kwarki mają rdzeniowe g₀ ∈ [0,817; 0,891]". Ten cykl pyta wyłącznie:
   **czy ta hipoteza dotyczyła właściwej kategorii fizycznej** — czy też testowała
   *wartość ogonową g_min* opisaną symbolem „g₀" przez przeciążenie notacji w sek08b:529.
2. **Wszystkie trzy werdykty = sukces.** NORM-COHERENT (sufit ważny ⇒ HALT-B potwierdzony,
   uczciwa demotacja predykcji kwarkowej) NIE jest „gorszy" od NORM-OVERLOAD. Pokusa
   ratunku (uwaga usera o ogonie czyni OVERLOAD atrakcyjnym) jest **najsilniejszym
   ryzykiem META** i dlatego reguła werdyktu jest LOCKED PRZED rachunkiem.
3. **Rescue-test poza zakresem.** Faktyczna rekalkulacja stosunków mas z rdzeniowym g₀ =
   OSOBNY cykl + osobny PR (forbidden #1). Nawet przy OVERLOAD ten audyt zatrzymuje się na
   „sufit VOIDED ⇒ HALT-B reopened + licencja", NIE „kwarki działają".
4. Każda faza = osobne „działaj"; raport + decision menu + wpis STATE.md po każdej.

### §0.4 — Pre-flight methodology read confirmation (per KICKOFF §2.6)

- [x] STATE.md #36-#39 (α=2/c₀ aksjomatyczne; #39 WIP: quark-mass HALT-B = następna rzecz)
- [x] CAŁY README + Phase_FINAL_close cyklu HALT-B (T11 ceiling logic; R3 normalization flag; §9.1 Path α)
- [x] op-FFS-quark-object Phase_FINAL (kwark = obiekt FFS z ogonem; Φ_0_local relacyjne)
- [x] sek08b_ghost_resolution.tex l.34-62 + l.326-374 + l.526-531 (g₀ rdzeniowe / ogon g_min / claim)
- [x] PRE_REGISTERED_FALSIFIERS §0-§2 (append-only; PR-014 candidate)
- [x] CALIBRATION_PROTOCOL §3.6 (BINDING) + CYCLE_KICKOFF_TEMPLATE §1-§3
- [ ] **Phase 1 (merytoryczne, NIE pre-flight):** dodatek app:ogon-masy (eq:J-ode, eq:J-tail, hip:J-mass-Atail4, A_tail(g₀)) + why_n3 PHASE2 (uniwersalny wzór + kalibracja)

**Sign-off:** Claudian @ 2026-06-25.

## §1 — Kontrakt kickoff L1 (observable-first)

```yaml
L1_native:
  output_observable: "Kategoria wielkości g w sek08b:529 [0,817;0,891] ∈ {A=rdzeniowe g₀, B=ogonowe g_min, C=inne}; wielkość pochodna: ważność sufitu T11=2,68×"
  measurement_instrument: "LIVE tekst rdzenia + dodatku ogon-masy + why_n3 + cykl-poprzednik (read-only); ZERO danych PDG kwarków w kategoryzacji"
  native_coefs_constrained: ["definicja g₀ (wejście eq:J-ode)", "g_min (ekstremum eq:J-tail)", "domena ε sufitu T11"]
  falsification_rule: "§0.2 README — NORM-OVERLOAD / NORM-COHERENT / NORM-INDETERMINATE, LOCKED 2026-06-25"
  pre_registration_date: "2026-06-25"
L2_framework_reduction:
  target_frameworks: ["spójność wewnętrzna sek08b ↔ app:ogon-masy ↔ why_n3 ↔ cykl HALT-B"]
  reduction_type: "internal-consistency-audit"
  failure_disposition: "L1-stands"
L3_falsification_map:
  - { bound: "g₀ rdzeniowe leptonów {1,0;2,0;2,34} (sek08b:346)", constrains: "[0,817;0,891] ⊂ {1..2,34}? — prima facie NIE", status: "pending Phase 1" }
  - { bound: "ogon g_min<1 (sek08b:60,328-329)", constrains: "[0,817;0,891] ⊂ region g<1? — prima facie TAK", status: "pending Phase 1" }
  - { bound: "T11 ceiling 2,68× (ε≈0,091 z [0,817;0,891])", constrains: "ważność sufitu zależna od kategorii", status: "pending Phase 1" }
  - { bound: "HALT-B verdict IMMUTABLE", constrains: "anti-Lakatos invariant — bez modyfikacji retroaktywnej", status: "LOCKED" }
```

## §2 — Pytania (zbiory CLOSED; kolejność WYMUSZONA)

**Q1 (kategoria — fast-audit, Phase 1):** czym JEST liczba g∈[0,817;0,891] w sek08b:529?
Zbiór kategorii **CLOSED**:

- **(A) RDZENIOWE g₀** — argument wejściowy eq:J-ode; dla leptonów {e≈1; μ≈2,0; τ≈2,34}.
- **(B) OGONOWE g_min** — ekstremum oscylacyjnego ogona g(r) (eq:J-tail), g_min<1, reżim
  przypróżniowy; wielkość **mierzona z zewnątrz** (≠ rest-frame solitonu) — wprost uwaga usera.
- **(C) INNA** — g* (próg duchowy), parametr podgrupy chiralnej, lub inna normalizacja;
  wymaga jawnego cytatu z app:ogon-masy. Brak rozstrzygnięcia ⇒ NORM-INDETERMINATE.

**Q2 (ważność sufitu — TYLKO po Q1):** rekonstrukcja, którą wielkość argument T11 wstawił
do `max(m)/min(m) ≤ (A_max/A_min)²·(1+ε)^(e²/2)`. Jeśli Q1=B, to ε=0,091 wzięte z `I`
jest **rozstępem wartości ogonowych**, a NIE rdzeniowych ⇒ sufit liczony na błędnej domenie.
Liczba kontrolna (strukturalna, BEZ fitu): rdzeniowy rozstęp leptonów {1;2,0;2,34} daje
ε_core ≈ 1,34 (vs 0,091) ⇒ `(1+ε_core)^(e²/2)` o rzędy większe — to pokazuje SKALĘ błędu
kategorii, nie ratuje kwarków (R4).

## §3 — Falsyfikatory (zbiory CLOSED; reguły LOCKED 2026-06-25)

### F-NORM-A — kategoria (Q1; werdykt Phase 1)

Klasy (CLOSED): **OVERLOAD / COHERENT / INDETERMINATE** (kryteria verbatim §0.2 README).

**NORM-OVERLOAD wymaga WSZYSTKICH (mechanicznie):**
1. eq:J-ode/eq:J-tail wykazują, że `I` indeksuje wartość profilu/ogona (kat. B), NIE wejście ODE;
2. sek08b:529 wiąże symbol „$g_0$" z tą wartością ogonową — jawny cytat (przeciążenie udokumentowane);
3. warunek wspierający: `I ⊄ conv(D_core)` (cały <1) — NIE wystarczający sam (R2).

**NORM-COHERENT wymaga:** eq:J-ode wykazuje, że `I` to faktyczny argument wejściowy ODE
dla kwarków (kat. A), różny od leptonowego {1..2,34} — wtedy sufit T11 ważny.

**NORM-INDETERMINATE:** dokumenty nie wystarczają / sprzeczne ⇒ GAP.

### F-NORM-B — ważność sufitu (Q2; werdykt Phase 1, sprzężony z A)

| F-NORM-A | Sufit T11 | Dyspozycja HALT-B |
|---|---|---|
| OVERLOAD | **VOIDED** (błędna kategoria domeny) | **reopened** → INDETERMINATE-PENDING-RESCUE + licencja osobnego cyklu |
| COHERENT | **VALID** | **CONFIRMED-STRENGTHENED** + PR-014 formalizacja + demotacja predykcji |
| INDETERMINATE | nierozstrzygnięty | **status quo** + housekeeping definicji |

**Protokół value-blind:** werdykt F-NORM-A jest WYLICZONY z T1-T6 (definicje tekstu),
LOCKOWANY PRZED jakimkolwiek obliczeniem dotyczącym stosunków mas kwarków. Żaden wybór
kategorii nie może być motywowany tym, która „ratuje" 80000× (forbidden #3).

## §4 — Analytical pre-derivation (§3.6.1-§3.6.9; PRZED ekstrakcją)

### §4.1 — Co tekst MÓWI prima facie (przed pełną ekstrakcją app:ogon-masy)

Z sek08b (lektura §0.4):
- **Rdzeniowe g₀** (rem:heavy-solitons l.327, l.346): „Profil solitonu z dużym $g_0 \gg 1$
  (ciężkie cząstki: mion, tauon)"; „Mion ($g_0 \approx 2{,}0$) i tauon ($g_0 \approx 2{,}34$)".
  ⇒ rdzeniowe g₀ leptonów: e≈1, μ≈2,0, τ≈2,34 (rozstęp ≈ 2,34×).
- **Ogon** (l.60, l.328-329, l.347-349): „oscylacyjny ogon, g_min < 1"; „Ich ogon
  oscylacyjny ($g \approx 1$ ... eq:J-tail) pozostaje w pełni w reżimie próżniowym".
  ⇒ region ogonowy: g<1, blisko g*.
- **Claim universalności** (l.528-529): „ten sam ODE działa na leptony i kwarki
  ($g_0 \in [0{,}817;\, 0{,}891]$)".

**Napięcie prima facie:** [0,817;0,891] (cały <1) NIE pokrywa się z rdzeniowym zbiorem
leptonów {1;2,0;2,34} (cały ≥1), ALE leży dokładnie w regionie ogonowym g_min<1. To
**wspiera hipotezę B (OVERLOAD)**, lecz NIE rozstrzyga — rozstrzygnięcie wymaga eq:J-ode
(co jest argumentem) + jawnej semantyki „$g_0$" w l.529 (Phase 1). Możliwe też, że l.529
opisuje **odrębne** rdzeniowe g₀ kwarków akurat w oknie <1 (kat. A, COHERENT) — wtedy
kwarki byłyby „lekkimi" solitonami (g₀<1) o słabszym sprzężeniu niż elektron, co jest
fizycznie nieoczywiste i samo w sobie wymaga audytu.

### §4.2 — Dlaczego to ma znaczenie dla sufitu (logika T11)

Sufit poprzednika: `max(m)/min(m) ≤ (A_max/A_min)²·(1+ε)^(e²/2)`, ε≈0,091 z `I`.
- Jeśli `I` = rdzeniowe g₀ (A): ε≈0,091 poprawne ⇒ sufit 2,68× ważny ⇒ COHERENT.
- Jeśli `I` = ogonowe g_min (B): rdzeniowe g₀ kwarków jest WIELKOŚCIĄ NIEOGRANICZONĄ tym
  przedziałem; analog leptonowy (rdzeń 1→2,34 przy ogonie ~1) sugeruje, że rdzeniowy
  rozstęp może być znacznie szerszy ⇒ ε_core ≫ 0,091 ⇒ sufit dużo wyższy ⇒ argument 2,68×
  **nieważny dla rdzeniowego g₀** (VOIDED). UWAGA (R4): to NIE dowodzi reprodukowalności —
  dowodzi tylko, że no-go policzono na złej zmiennej.

### §4.3 — Enumeracja założeń (§3.6.8)

(a) tekst rdzenia LIVE jest źródłem prawdy o definicjach (read-only; zero edycji w tym cyklu);
(b) eq:J-ode definiuje g₀ jako parametr wejściowy/warunek brzegowy; eq:J-tail definiuje
asymptotykę profilu g(r) — rozróżnienie do potwierdzenia Phase 1; (c) kalibracja leptonowa
why_n3 jest LIVE i niezmieniana; (d) żadne podstawienie wartości PDG kwarków w fazie
kategoryzacji (circularity guard FP T8); (e) e²/2 ≈ 3,69 (L08-e² LIVE) — niezmieniane.

### §4.4 — Klasyfikacja stałych (§3.6.13)

| Wielkość | Klasa | Nota |
|---|---|---|
| g₀ rdzeniowe | (α) TGP_FUNDAMENTAL (wejście ODE) | przedmiot rozróżnienia |
| g(r), g_min | (α) pochodna profilu | ekstremum ogona, eq:J-tail |
| A_tail(g₀) | (α) pochodna ODE | hip:J-mass-Atail4 |
| e²/2 ≈ 3,69 | (α) LOCKED (L08-e²) | wykładnik, niezmieniany |
| {1,0; 2,0; 2,34} | LITERATURE (sek08b:346) | rdzeniowe g₀ leptonów |
| [0,817; 0,891] | przedmiot audytu | sek08b:529 |
| PDG quark masses | comparison-only | ZAKAZANE w kategoryzacji (guard FP) |
| Nowe stałe | budżet: 0 | luka = GAP |

## §5 — Plan faz (fast-audit first)

| Faza | Zakres | Werdykt | Gate |
|---|---|---|---|
| **Phase 1 — FAST-AUDIT (Q1+Q2)** | ekstrakcja eq:J-ode/eq:J-tail/A_tail (T1) → przynależność `I` do D_core (T2) i R_tail (T3) → cytat semantyki „$g_0$" l.529 (T4) → rekonstrukcja domeny T11 (T5) → werdykt kategorii (T6) → sufit pod A vs B (T7) → circularity guard (T8) → S05/HALT-B IMMUTABLE (T9) | F-NORM-A + F-NORM-B | „działaj" |
| Phase FINAL | agregat; dyspozycja HALT-B; korekta sek08b:529 (jeśli OVERLOAD); PR-014 (jeśli COHERENT) lub licencja rescue-test; DOUBTS register; propagacja | — | „działaj" |

Każda faza: PhaseN_audit.md + PhaseN_sympy.py/.txt (PASS/FAIL per FP; 0 hardcoded;
circularity-guard FP: free-symbols audit na PDG quark masses w każdym skrypcie).

## §6 — Forbidden moves (LOCKED; 12)

1. Rekalkulacja 5 stosunków mas kwarków z rdzeniowym g₀ (= rescue-test; OSOBNY cykl/PR).
2. Modyfikacja/retraktacja werdyktu HALT-B poprzednika (IMMUTABLE; oraz PR-004/PR-010 read-only).
3. Wybór kategorii g motywowany tym, która „ratuje" 80000× (value-blind violation).
4. Użycie JAKIEJKOLWIEK danej PDG kwarka w fazie kategoryzacji (guard FP w każdym sympy).
5. Post-hoc redefinicja przedziału [0,817;0,891] lub progów werdyktu.
6. Edycja tekstu rdzenia w tym cyklu (read-only; korekty = osobny housekeeping po werdykcie).
7. Nowe stałe/pola (S05; budżet 0; luka = GAP).
8. Hardcoded T_pass.
9. Miękkie domknięcie / przeciąganie po jednoznacznym werdykcie (każdy z trzech = zamknięcie).
10. Nadinterpretacja „sufit VOIDED ⇒ kwarki reprodukowalne" (R4; OVERLOAD zdejmuje no-go, nie dowodzi pozytywu).
11. Scope-creep do rescue-test lub innego cyklu w tej sesji.
12. Append/LOCK PR-014 w rejestrze przed Phase FINAL i bez decyzji user.

## §7 — Risk register

| ID | Ryzyko | Severity | Mitygacja |
|---|---|---|---|
| R-1 | Pokusa ratunku (uwaga usera o ogonie czyni OVERLOAD atrakcyjnym) | **META/HIGH** | reguła §0.2 LOCKED przed rachunkiem; werdykt WYLICZONY z definicji; rescue = osobny PR |
| R-2 | „$g_0$" w rdzeniu faktycznie przeciążone | MED | wtedy OVERLOAD UCZCIWY; wymaga jawnego cytatu eq:J-ode (T4) — nie domysłu |
| R-3 | app:ogon-masy definiuje `I` jako jeszcze inną wielkość (g*, podgrupa chiralna) | MED | kategoria C + NORM-INDETERMINATE; zbiór CLOSED pokrywa |
| R-4 | Nadinterpretacja OVERLOAD jako „kwarki działają" | HIGH | forbidden #10; werdykt zatrzymuje się na VOIDED+licencja |
| R-5 | BD-drift | LOW | brak propagatora; profil klasyczny; self-audit Phase FINAL |
| R-6 | COHERENT odebrane jako „porażka" zamiast wynik | MED | §0.0 pkt 2: demotacja predykcji = pełnoprawny uczciwy wynik |

## §8 — Anticipated outcomes (INFORMATIONAL; kierunki pokusy zapisane PRZED rachunkiem)

1. **Prima facie sygnał wskazuje NORM-OVERLOAD** ([0,817;0,891] cały <1, w regionie
   ogona, poza rdzeniowym {1..2,34}). ALE: zapisane jest, że to NIE wystarcza — wymaga
   eq:J-ode pokazującego, że `I` to wartość profilu, ORAZ cytatu semantyki „$g_0$" l.529.
   Możliwy zwrot: l.529 może opisywać odrębne rdzeniowe g₀<1 kwarków (COHERENT) — wtedy
   sufit stoi i HALT-B potwierdzony.
2. **Pokusa #1 (META):** „uwaga usera potwierdza OVERLOAD" — antidotum: user dał HIPOTEZĘ
   fizyczną (ogon ≠ rest-frame), nie werdykt; rozstrzyga tekst eq:J-ode, nie autorytet.
3. **Pokusa #2:** przy OVERLOAD — policzenie sufitu rdzeniowego i ogłoszenie „kwarki
   reprodukowalne". Antidotum: forbidden #10; to osobny rescue-test z własnym falsyfikatorem.
4. **NORM-COHERENT jest realnym i wartościowym wyjściem:** czyściej domyka headline
   (lepton derived / quark (H)) i formalizuje PR-014.

## §9 — Anti-Lakatos compliance (Phase 0)

Zbiór kategorii CLOSED przed rachunkiem ✓ · reguła werdyktu trójdzielna LOCKED z wyliczalnymi
kryteriami ✓ · werdykt HALT-B poprzednika NIETYKALNY (audyt zakresu ważności, nie poprawności) ✓ ·
PDG quark masses comparison-only z guardem FP ✓ · rescue-test wykluczony z sesji (osobny PR) ✓ ·
nadinterpretacja OVERLOAD zakazana (R4/#10) ✓ · 0 nowych stałych ✓ · 0 edycji rdzenia ✓ ·
anticipated outcomes z kierunkami pokusy (w tym META: autorytet usera ≠ dowód) ✓ ·
COHERENT zadeklarowany jako pełnoprawny sukces (anty-bias ku ratunkowi) ✓.

**Phase 0 LOCKED 2026-06-25. Następny krok: Phase 1 FAST-AUDIT — wymaga user „działaj".**
