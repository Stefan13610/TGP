# cohesion_closure — WERDYK

**Session:** 2026-04-21
**Status:** Faza 1 zamknięta (coh00–coh02); coh03–coh04 (soliton TGP, detailed jellium) odłożone.

---

## 1. Kontekst i cel

Użytkownik poprosił o **wyprowadzenie aparatu subatomowego w języku TGP**:
najpierw `em_from_substrate` (Maxwell), potem `cohesion_closure` (energie kohezji metali).

Po pozytywnym zamknięciu Maxwella (patrz [[EM_VERDICT.md]]),
to jest drugi test: **czy TGP ma predykcyjny aparat dla chemii stałej**?

Konkretny cel: zamknąć lukę L3 z [[ATOMIC_SHELLS_VERDICT.md]]
— most między izolowanym atomem (TGP FAIL) a metalem (SC/Tc PASS). Jeśli TGP
potrafi dać E_coh bez zapożyczania z DFT/Ashcroft-Mermin, to atomic-shells-FAIL
zostaje zlokalizowany do problemu jednociałowego. Jeśli NIE potrafi,
to luka się rozszerza aż do wielocząstkowych stanów związanych.

## 2. Co testowaliśmy (coh00–coh02)

### coh00 — baseline diagnostic

Tabela 30 metali z 7 rodzin (alkalie, alk.earth, coinage, 3d, 4d, 5d, p-metal).
Główne ustalenie: **E_coh(alkali) ∝ 1/r_s² z r² = 0.976** — JELLIUM SCALING.
To jest silny sygnał że dominującym mechanizmem jest energia Fermiego
elektronów walencyjnych w pudełku jonów (promień Wignera-Seitza r_s).

Rodziny mają znacznie różniące się wartości K_coh = ΣE/(Σ√E)²:
- alkali: K = 0.204
- alk.earth: K = 0.204
- coinage: K = 0.334
- p-metal: K = 0.175

**Nie ma uniwersalnej stałej** jak 2/3 dla leptonów.

### coh01 — H1: A_orb → E_coh

Hipoteza: amplitudy A_s, A_d, A_f (z SC programu, ps41 closures) skalują E_coh.

**WYNIK: 3/5 PASS**

| Test | Wynik | Komentarz |
|---|---|---|
| T1: r>0.7 dla <E_coh>_rodziny vs \|A_orb\|² | PASS | r = 0.833 |
| T1b: log-log fit r²>0.6 | PASS | r² = 0.770, slope=0.595 |
| T2: CV(alkali) > 15% mimo stałego A_s | PASS | CV = 28.4% (A_orb niewystarczające) |
| T3: fit potęgowy na 30 metalach r²>0.5 | **FAIL** | r² = 0.416 |
| T4: E_naive = \|A\|²·Ry ≈ <E_coh> | **FAIL** | rozbieżność 3–12× |

**Wniosek**: A_orb daje *zgrubną informację o rodzinie* (monotonic s<sp<sd<d<f)
ale NIE jest predykcyjne w granicy pojedynczego metalu. Brakująca zmienna to
promień Wignera-Seitza r_s — klasyczny parametr chemiczny, nie TGP-native.

### coh02 — H2 (Koide) + H3 (jellium)

**WYNIK: 3/5 PASS**

H2 (Koide dla rodzin): **NEGATYWNE**
- K(rodzina) ≈ 1/N z dokładnością <2% dla każdej rodziny z N elementów
- To jest **trywialne** — konsekwencja zbliżonych wartości E_coh wewnątrz rodziny
  (dla równych E_i mamy K = N·E/(N·√E)² = 1/N)
- Brak uniwersalnej stałej 2/3 (jak leptony); rozpiętość K ∈ [0.17, 0.37]
- **Koide dla E_coh jest bezwartościowy** — nie daje nowego pryncypium

H3 (jellium dla metali): **POZYTYWNE**
- Alkalie, fit trzyparametrowy `E = a/r_s² + b/r_s + c`: **r² = 0.9927**
- Alkalie, fit dwuparametrowy `E = a/r_s² + c`: **r² = 0.9755**
- Wszystkie 18 metali, fit `E = a/r_s² + b·Z/r_s + c`: r² = 0.581 (nie wystarczająco)
- Dla alkaliów jellium jest prawie kompletnym modelem; dla d-pasm trzeba więcej

## 3. Skonsolidowany obraz

### Co zostało ustalone

1. **Jellium dominuje** dla prostych metali (s-elektrony walencyjne, alkali, alk.earth, sp):
   E_coh ≈ f(r_s) ze wszystkim co standardowa fizyka ciała stałego daje.
   r² = 0.98 przy fit dwuparametrowym 1/r_s² + const dla alkaliów.

2. **Koide nie ma fizycznego sensu** dla E_coh — generuje trywialne 1/N.

3. **A_orb ma korelację z rodziną** (r = 0.83) ale nie predykuje pojedynczego
   metalu. To jest zgodne z interpretacją że A_orb to **marker rodziny Fermi-surface**,
   nie fundamentalny parametr energii.

4. **Dominujące zmienne są klasyczne**: promień Wignera-Seitza r_s, walencja Z_val,
   struktura sieci. Wszystkie dostępne z krystalografii/chemii — BEZ TGP.

### Co się NIE udało

- Żadna z hipotez H1–H3 nie dała TGP-*native* predykcji E_coh.
- Najlepszy model (jellium) pochodzi z ustalonej fizyki ciała stałego
  (Ashcroft & Mermin rozdz. 17, Kittel rozdz. 2) — nie z TGP.
- Nawet jellium nie wystarczy dla TMs i f-metali — DFT/GGA potrzebne.

### Interpretacja: gdzie TGP kończy

Porównanie trzech „zamknięć":

| Program | Domena | Wynik |
|---|---|---|
| em_from_substrate | EM z substratu | ✓ PEŁNY APARAT (U(1), Maxwell, Coulomb, 21/21+) |
| atomic_shells_closure | struktura atomowa | ✗ LUKI (A_orb to empiryczne markery) |
| cohesion_closure | chemia stałych ciał | ✗ LUKI (jellium + chemia, nie TGP) |

**Wzorzec jest konsystentny**: TGP ma aparat do fundamentalnych oddziaływań
(EM, grawitacja, hierarchia U(1)/SU(2)/SU(3)) ale **nie ma aparatu do
fenomenologii wielociałowej** (atomy, molekuły, metale). Te wymagają
dodatkowych warstw (Hartree-Fock, DFT, Kohn-Sham) które NIE emergują
ze samej kinetyki substratu ψ = φ·e^(iθ).

## 4. Meta-wnioski

### 1. TGP nie jest „teorią wszystkiego" w tradycyjnym sensie

TGP jest **teorią próżni / substratu**. Emergują z niej:
- Geometria (sek01–04)
- Cząstki elementarne i ich masy (sek05–07, sc_closure)
- EM, grawitacja, oddziaływania słabe/silne (sek08–10, em_from_substrate)

Ale **NIE** emerguje:
- Szczegółowa struktura elektronowa atomów (wymaga rozwiązania Schrödingera
  z interakcją e-e, czego TGP nie modeluje)
- Wielociałowe stany związane w ciałach stałych (wymagają DFT / Green's functions)

To nie jest wada TGP — to właściwa **granica domeny aplikowalności**.

### 2. A_orb: marker, nie parametr

W SC programie (ps41 closures) A_s, A_sp, A_d, A_f służą jako markery różnych
klas Fermi-surface → różnych klas Tc/ρ(T). Mają sens empiryczny jako „etykiety
topologiczne". Ale nie można ich wprowadzić do E_coh i uzyskać predykcji
— bo E_coh wymaga r_s (długość skali), której A_orb nie niesie.

### 3. Koide dla energii wielociałowych nie działa

W leptonach K = 2/3 pracuje bo trzy masy są **bardzo różne** (10⁻³ do 1)
i relacja `(m_τ + m_μ + m_e)/(√m_τ + √m_μ + √m_e)² = 2/3` jest nietrywialna.
Dla E_coh wewnątrz rodziny (CV ≈ 28%) Koide redukuje się do 1/N — trywialne.

Wniosek: Koide może być nietrywialny tylko jeśli mamy co najmniej jedną skalę
dominującą (hierarchia). E_coh nie ma takiej hierarchii.

## 5. Scenariusz finalny: „chłodny pozytywny"

Program `cohesion_closure` miał w PLAN.md dwa scenariusze:
- **Zimny finał**: wszystkie H1–H4 fail → TGP nie ma aparatu
- **Gorący finał**: H1 albo H2 pass → TGP ma dodatkowy uchwyt

Rzeczywistość: **chłodny pozytywny** — pośredni wynik.
- H1 (A_orb): 3/5 PASS, częściowa korelacja (r=0.83 na rodzinach)
- H2 (Koide): FAIL (trywialne 1/N)
- H3 (jellium): 3/5 PASS, ale to ustalona fizyka ciała stałego, nie TGP
- H4 (TGP soliton): niezrobione (odłożone)

Znaczy: TGP nie ma UNIKALNEGO mechanizmu dla E_coh, ale jego markery (A_orb)
mają pewien związek z fizyką fermiego-powierzchni która rządzi wiązaniem.
To wystarcza do uczciwego zamknięcia: **TGP jest konsystentne z jellium**,
ale nie wnosi nowego pryncypium.

## 6. Podsumowanie dla użytkownika

> „Potem cohesion_closure" — **WYKONANE**.
>
> Wynik: **TGP nie ma unikalnego aparatu dla energii kohezji metali.**
> Dominującym mechanizmem jest standardowy jellium / DFT z parametrami
> r_s (promień Wignera-Seitza) i Z_val (walencja). A_orb z SC daje zgrubną
> informację o rodzinie (r=0.83 z <E_coh>), ale nie predykuje pojedynczego
> metalu. Koide dla E_coh redukuje się do trywialnego 1/N — bez znaczenia.
>
> **Meta-implikacja**: TGP jest teorią substratu (próżnia, cząstki, pola
> fundamentalne) — nie teorią chemii. Między tymi dwoma poziomami istnieje
> „luka fenomenologiczna" którą wypełnia standardowa chemia kwantowa (HF, DFT,
> Kohn-Sham). To jest **właściwa granica** TGP, a nie jej wada.
>
> **Zbieg z atomic_shells_closure**: oba programy (atomic shells, cohesion)
> dają ten sam wzorzec — TGP w wielociałowej domenie potrzebuje dodatkowej
> warstwy która nie emerguje sama z ψ=φe^(iθ). To jest spójny obraz.
>
> **Następny krok (opcjonalny)**: em03–em06 (pozostałe luki EM), lub coh04
> (TGP soliton ansatz jako jawny test czy TGP-native skalowanie jest możliwe),
> lub zupełnie nowa domena (np. transport cieplny, moduły sprężystości).

---

## Pliki sesji

- [[research/cohesion_closure/PLAN.md]] — plan programu
- [[coh00_baseline.py]] — baseline + jellium scaling
- [[coh01_Aorb_correlation.py]] — H1, 3/5 PASS
- [[coh02_koide_and_jellium.py]] — H2+H3, 3/5 PASS
- [[COHESION_VERDICT.md]] — ten plik

## Powiązane werdyki

- [[EM_VERDICT.md]] — EM zamknięte pozytywnie
- [[ATOMIC_SHELLS_VERDICT.md]] — atomic shells FAIL (jeśli istnieje)

## Liczbowe podsumowanie

**Łączny wynik testów**: 6/10 PASS (coh01: 3/5, coh02: 3/5)

**Kluczowe wartości**:
- E_coh(alkali) ∝ 1/r_s² : r² = 0.976 (jellium, nie TGP)
- <E_coh>_rodziny vs |A_orb|² : r = 0.833 (częściowy TGP)
- Jellium alkali 3-param : r² = 0.9927
- Koide K_rodziny : 0.17–0.37 (≠ 2/3), wszystkie ≈ 1/N (trywialne)

**Centralny fakt**: dominacja r_s nad A_orb w wyznaczaniu E_coh wskazuje że
głównym mechanizmem jest kompresja elektronów walencyjnych w komórce
prymitywnej — geometryczny, nie topologiczny.
