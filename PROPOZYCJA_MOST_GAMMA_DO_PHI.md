# Propozycja badawcza: most `Γ -> Φ` (substrat -> pole ciągłe)

## Cel

Celem tego dokumentu jest sformułowanie **minimalnego, ale publikowalnego programu domknięcia** przejścia

`substrat Γ -> pole blokowe Φ_B -> pole ciągłe Φ`.

Nie chodzi tu o pełny dowód wszystkich elementów ERG i operatora Banacha, lecz o takie zamknięcie, które:

- daje ścisły sens przejściu od mikrodynamiki do pola ciągłego,
- uzasadnia lokalny funkcjonał efektywny dla `Φ`,
- pokazuje, skąd bierze się struktura równania pola TGP,
- oddziela to, co już jest policzone numerycznie, od tego, co nadal pozostaje hipotezą.

## Diagnoza stanu obecnego

Z materiału w repo wynika spójny obraz:

- `Φ` jest już naturalnie zdefiniowane jako obserwabla blokowa substratu:
  `Φ_B(x) = N_B^{-1} Σ_{i in B(x)} φ_i^2`.
- warunek skali dla sensownej granicy continuum jest już zapisany:
  `a_sub << L_B << ξ`, gdzie `L_B = b a_sub`, a `ξ = 1/sqrt(γ)`,
- istnieją mocne argumenty wspierające zbieżność:
  ERG/LPA', de Giorgi, homogenizacja,
- ale nadal nie ma pełnego mostu:
  `Γ -> Φ` w sensie twierdzenia o zbieżności i identyfikacji operatora efektywnego.

Najważniejsza luka nie dotyczy więc pojedynczego współczynnika, lecz braku formalnego twierdzenia typu:

> z blokowania substratu istnieje granica ciągła i jej dynamika ma postać lokalnego równania TGP.

## Proponowany cel minimum

Najrozsądniejszym celem nie jest na dziś pełne zamknięcie `CG-1`, `CG-3`, `CG-4` w najmocniejszej wersji, lecz udowodnienie słabszego twierdzenia:

## Teza robocza

Niech

`Φ_B(x) = N_B^{-1} Σ_{i in B(x)} φ_i^2`

będzie polem blokowym na substracie `Γ` z lokalnym hamiltonianem i długością korelacji `ξ_corr`.

Jeżeli:

1. `sup_B ||Φ_B||_{L^2(K)} < ∞` dla każdego zwartego `K`,
2. energia swobodna po blokowaniu spełnia dolny szacunek gradientowy
   `F_B[Φ_B] >= c_* ||∇Φ_B||^2_{L^2(K)} - C_K`, gdzie `c_* > 0`,
3. korelacje dwupunktowe maleją szybciej niż `exp(-r/ξ_corr)`,
4. zachodzi separacja skal `a_sub << L_B << ξ_corr`,

to rodzina `{Φ_B}` jest prezwarta słabo w `H^1_loc` i mocno w `L^2_loc`, a punkty skupienia spełniają równanie Eulera-Lagrange'a pewnego **lokalnego funkcjonału efektywnego**.

To byłoby już realne zamknięcie mostu w wersji:

`weak continuum theorem`.

## Dlaczego ten cel jest właściwy

Taki rezultat jest wystarczająco mocny, bo:

- przenosi TGP z poziomu intuicji coarse-grainingu na poziom analizy funkcjonalnej,
- nie wymaga od razu dowodu, że operator blokowania jest kontrakcją Banacha,
- wystarcza do legitymizacji dalszych przejść `Φ -> g_{μν}`,
- pozwala traktować ciągłe równanie pola jako wynik granicy efektywnej, a nie tylko motywowany ansatz.

Innymi słowy: nie trzeba najpierw udowodnić wszystkiego. Trzeba najpierw udowodnić **najkrótszy most, który naprawdę niesie fizykę**.

## Proponowany łańcuch wyprowadzenia

### Krok 1. Ustalić ontologię zmiennej blokowej

Za punkt wyjścia przyjąć:

- zmienne mikroskopowe `φ_i`,
- hamiltonian lokalny na `Γ`,
- pole blokowe `Φ_B = <φ^2>_B`,
- próżnię jako stan termalny substratu w fazie `T < T_c`.

To jest ważne, ponieważ wtedy `Φ` nie jest nowym obiektem dodanym z zewnątrz, tylko **obserwablą substratu**.

### Krok 2. Ustalić domenę stosowalności granicy continuum

Formalnie używać tylko reżimu:

`a_sub << L_B << ξ`.

Znaczenie:

- `L_B >> a_sub` daje prawo dużych liczb dla średniej blokowej,
- `L_B << ξ` uzasadnia rozwinięcie gradientowe,
- continuum ma sens tylko blisko przejścia fazowego, gdzie `ξ/a_sub` jest duże.

To powinno być zapisane nie jako uwaga heurystyczna, ale jako jawny warunek twierdzenia.

### Krok 3. Wyprowadzić lokalny funkcjonał efektywny po blokowaniu

Najbardziej naturalny kandydat ma postać:

```text
F_eff[Φ] = ∫ d^3x [ K1(Φ) (∇Φ)^2 + U(Φ) - q ρ Φ ].
```

W wersji minimalnej nie trzeba od razu znać pełnych zamkniętych wzorów na `K1` i `U`.
Wystarczy wykazać:

- lokalność,
- dodatniość sektora kinetycznego dla `Φ > 0`,
- istnienie stabilnej próżni `Φ = Φ_0`,
- poprawny limit liniowy wokół `Φ_0`.

To jest dokładnie moment, w którym mikrodynamika zaczyna przechodzić w pole efektywne.

### Krok 4. Użyć kompaktowości zamiast pełnego punktu stałego

Najbardziej realistyczna ścieżka analityczna:

- ograniczenie `L^2` dla `Φ_B`,
- dolny szacunek energii dający kontrolę `∇Φ_B` w `L^2`,
- Rellich-Kondrachov,
- zbieżność podsekwencji w `L^2_loc`,
- słabe przejście do granicy w równaniu Eulera-Lagrange'a.

To pozwala zamienić problem

`czy istnieje pełna kontrakcja operatora blokowania?`

na słabszy i bardziej dostępny problem

`czy istnieje granica ciągła podsekwencji z dobrze określonym funkcjonałem efektywnym?`

### Krok 5. Zidentyfikować operator efektywny z równaniem TGP

Po wariacji funkcjonału dostajemy równanie typu:

```text
2 K1(Φ) ∇²Φ + K1'(Φ) (∇Φ)² - U'(Φ) + qρ = 0.
```

Po przeskalowaniu celem jest otrzymanie struktury:

```text
∇²Φ + α_eff (∇Φ)² / Φ + β_eff Φ²/Φ_0 - γ_eff Φ³/Φ_0² = - q_eff Φ_0 ρ.
```

To jest punkt, w którym trzeba pokazać nie tylko istnienie PDE, ale zgodność jego struktury z TGP.

### Krok 6. Uzasadnić `α_eff ~ 2` z przejścia `φ -> Φ = φ²`

To jest bardzo ważne, bo właśnie tu teoria może pokazać, że współczynnik nieliniowy nie jest dobierany ręcznie.

Jeżeli `Φ = φ²`, to:

```text
∇Φ = 2φ ∇φ,
(∇Φ)² = 4Φ (∇φ)²,
(∇φ)² = (∇Φ)² / (4Φ).
```

Każdy regularny człon gradientowy dla amplitudy `φ` po przejściu do zmiennej `Φ` generuje więc naturalnie składnik

`(∇Φ)² / Φ`.

To jest najlepszy kandydat na wyprowadzenie `α = 2` jako konsekwencji geometrii zmiennej, a nie osobnego postulatu.

## Co dokładnie trzeba zbadać

### 1. Dodatniość efektywnej sztywności

To jest najważniejszy punkt techniczny.
Bez `K1(Φ) > 0` dla `Φ > 0` nie ma ani eliptyczności, ani kompaktowości, ani kontroli nad limitem continuum.

Minimalny cel:

- pokazać analitycznie albo numerycznie, że po blokowaniu sektor gradientowy pozostaje dodatni,
- oszacować dolne ograniczenie `c_*`.

### 2. Lokalność po blokowaniu

Trzeba sprawdzić, czy nielokalne ogony po coarse-grainingu są tłumione jako

`O(L_B/ξ)` albo szybciej.

Jeżeli tak, można uczciwie przejść do lokalnego funkcjonału efektywnego.

### 3. Identyfikacja potencjału efektywnego

Należy wyznaczyć, czy histogramowe lub RG-owe `V_eff(Φ)` rzeczywiście daje rozwinięcie zgodne z:

- próżnią przy `Φ = Φ_0`,
- dodatnią masą `m_sp² > 0`,
- strukturą `β_eff`, `γ_eff` zgodną z istniejącym równaniem pola.

### 4. Kontrola błędu continuum

Most `Γ -> Φ` nie powinien kończyć się na samym istnieniu granicy.
Potrzebny jest jeszcze błąd:

`||Φ_B - Φ|| <= ε(L_B/ξ, a_sub/L_B)`.

Nawet oszacowanie jakościowe byłoby bardzo cenne, bo wyznacza domenę stosowalności formalizmu ciągłego.

## Program badawczy: wersja analityczna

### Etap A1. Lemat kompaktowości

Sformalizować słabe twierdzenie continuum:

- ograniczenie `Φ_B` w `L^2_loc`,
- ograniczenie `∇Φ_B` przez energię,
- prezwartość w `H^1_loc`/`L^2_loc`.

To jest najważniejszy pierwszy rezultat publikowalny.

### Etap A2. Lemat lokalnego funkcjonału

Pokazać, że energia swobodna po blokowaniu ma w granicy postać lokalną drugiego rzędu:

```text
F_eff[Φ] = ∫ [K1(Φ)(∇Φ)² + U(Φ)] d^3x + reszta,
```

gdzie `reszta -> 0` dla `L_B/ξ -> 0`.

### Etap A3. Lemat `α_eff`

Udowodnić, że przejście z amplitudy do jej kwadratu wymusza pojawienie się członu

`(∇Φ)² / Φ`,

a jego kanoniczna normalizacja daje `α_eff = 2 + O(η)`.

### Etap A4. Identyfikacja współczynników

Połączyć:

- parametry z coarse-grainingu,
- wynik LPA'/ERG,
- rozwinięcie wokół `Φ_0`,

aby zidentyfikować `β_eff`, `γ_eff`, `m_sp`.

### Etap A5. Twierdzenie końcowe

Złożyć wynik w formie:

> Dla substratu `Γ` spełniającego warunki regularności i separacji skal, pole blokowe `Φ_B` posiada granicę ciągłą `Φ` w sensie słabym, a `Φ` spełnia efektywne równanie TGP z błędem kontrolowanym przez `L_B/ξ` oraz `a_sub/L_B`.

## Program badawczy: wersja numeryczna

Repo już ma dobry punkt startowy do programu wspierającego dowód.

### Etap N1. Monte Carlo + block averaging

W praktyce mierzyć:

- `Φ_B` dla wielu rozmiarów bloków `b`,
- histogramy `Φ_B`,
- korelatory gradientowe,
- długość korelacji `ξ`.

### Etap N2. Ekstrakcja `V_eff(Φ)`

Z rozkładu histogramowego wyznaczyć:

- minimum `Φ_min`,
- krzywiznę wokół minimum,
- współczynniki sześcienne i czwartorzędowe,
- zgodność z kształtem potencjału TGP.

### Etap N3. Ekstrakcja sektora kinetycznego

Z korelatorów różnic blokowych wyznaczyć efektywny `C_kin` i sprawdzić:

- dodatniość,
- stabilność przy zmianie `L`,
- zgodność z przewidywaniem continuum.

### Etap N4. Finite-size scaling

Nie wystarczy pojedyncza symulacja.
Potrzebna jest ekstrapolacja:

- `L -> ∞`,
- kilka bloków `b`,
- reżim blisko `T_c`,
- porównanie z `ξ`.

### Etap N5. Test zgodności operatora

Najmocniejszy test numeryczny:

- podstawiać zmierzone `Φ_B` do efektywnego PDE,
- sprawdzać rezydua,
- badać, czy rezydua maleją wraz z `a_sub/L_B` i `L_B/ξ`.

To byłby numeryczny odpowiednik identyfikacji `K_hom = K_TGP`.

## Kryterium uznania mostu za zamknięty

Most `Γ -> Φ` należałoby uznać za praktycznie domknięty, jeśli równocześnie będą spełnione cztery warunki:

1. istnieje formalne słabe twierdzenie continuum dla `Φ_B`,
2. lokalny funkcjonał efektywny jest wyprowadzony lub bardzo silnie uzasadniony,
3. struktura operatora TGP jest odzyskana z coarse-grainingu,
4. istnieje niezależne wsparcie numeryczne dla znaków i skali współczynników.

Bez punktu 3 teoria nadal ma tylko przejście `substrat -> jakieś pole`.
Dopiero punkt 3 daje przejście `substrat -> właśnie pole TGP`.

## Najważniejsze ryzyko

Największe ryzyko polega na tym, że pełna identyfikacja współczynników może okazać się zbyt mocna na pierwszy etap.

Dlatego rekomendowana strategia to:

- najpierw domknąć **istnienie granicy i lokalnego PDE**,
- potem domknąć **dokładną identyfikację współczynników**,
- dopiero na końcu próbować najmocniejszej wersji z pełnym punktem stałym operatora blokowania.

To minimalizuje ryzyko utknięcia na najbardziej ambitnym, ale niekoniecznym technicznie kroku.

## Rekomendacja końcowa

Najbardziej sensowna wersja robocza mostu `Γ -> Φ` brzmi:

> Pole `Φ` należy traktować jako granicę blokowej obserwabli `Φ_B = <φ²>_B` w reżimie `a_sub << L_B << ξ`, a celem formalnym powinno być słabe twierdzenie continuum prowadzące do lokalnego funkcjonału efektywnego i równania Eulera-Lagrange'a o strukturze zgodnej z TGP.

To jest dziś najkrótsza droga z:

`heurystyka coarse-grainingu`

do:

`publikowalnego wyprowadzenia substrat -> pole ciągłe`.

## Proponowane następne artefakty

Jeżeli ten kierunek ma być dalej rozwijany, sensowne byłyby kolejne dwa dokumenty:

1. `MOST_GAMMA_DO_PHI_LEMATY.tex`
   Zbiór lematów A1-A5 w stylu manuskryptu.

2. `PLAN_NUMERYCZNY_CG3_CG4.md`
   Operacyjny plan domknięcia części MC/RG dla `CG-3` i `CG-4`.

## Status realizacji (2026-04-12)

Oba artefakty zostały utworzone:

1. ✅ `dodatekQ2_most_gamma_phi_lematy.tex` — Lematy A1–A5 w pełnej formie LaTeX
   - A1 (kompaktowość): Szkic z dowodem via Rellich-Kondrachov
   - A2 (lokalność): Szkic z argumentem tłumienia ogonów
   - A3 (α=2): Propozycja z pełnym dowodem algebraicznym
   - A4 (współczynniki): Program z identyfikacją β=γ
   - A5 (twierdzenie): Szkic kompozytowy z kontrolą błędu
   - Włączony do main.tex jako Dodatek Q2

2. ✅ `PLAN_NUMERYCZNY_CG3_CG4.md` — Etapy N1–N5 z kryteriami PASS/FAIL
   - N1: MC + block averaging (bazowy, 2-3 dni)
   - N2: Ekstrakcja V_eff(Φ) (identyfikacja β, γ)
   - N3: Sektor kinetyczny (c* > 0, α = 2)
   - N4: Finite-size scaling (ekstrapolacja L → ∞)
   - N5: Test zgodności operatora (residua PDE)
   - Szacowany czas: 8–12 dni roboczych

### Następne priorytety

1. **Zamknąć A1**: uzasadnić c* > 0 (dolny szacunek gradientowy)
   - Droga: ERG (CG-2 daje K_IR/K_UV = 1.000) → przenieść na szacunek blokowy
   - Alternatywa: N3 (numeryczny pomiar c*)
2. **Napisać skrypt N1**: substrate_mc_cg3.py (bazowy MC)
3. **Napisać skrypt N3**: substrate_kinetic_sector.py (test c* > 0)
