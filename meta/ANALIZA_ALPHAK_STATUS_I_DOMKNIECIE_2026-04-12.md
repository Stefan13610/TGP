# Analiza luki `alpha_K ~= 8.56` i plan domkniecia

Data: 2026-04-12

## Cel

Ten dokument porzadkuje status parametru `alpha_K ~= 8.56` po sesjach `v45-v47b` i rozdziela dwa logicznie rozne problemy:

1. bezposrednie wyprowadzenie `alpha_K` z Warstwy I-II,
2. formalna supersesja starego `alpha_K` przez kanoniczny jezyk `g0* / phi-FP`.

Najwazniejszy wniosek jest prosty:

`alpha_K` nie jest juz pojedyncza "brakujaca liczba do policzenia", tylko miejscem styku dwoch formulacji sektora leptonowego.

## 1. Status kanoniczny

Obecnie w repo wspolistnieja dwa jezyki opisu.

### Jezyk legacy

- `alpha_K_old ~= 8.56`
- relacje typu `r21 <-> alpha_K <-> a_Gamma`
- stare skany bifurkacyjne i energetyczne
- interpretacja: parametr efektywny starej formulacji ODE

### Jezyk kanoniczny

- `g0^e`, `g0^mu`, `g0^tau`
- punkt staly `g0*`
- regula `phi-FP`
- relacje `A_tail`, `Q_K = 3/2`, `r21`, `r31`

Interpretacja:

- jezyk kanoniczny jest obecnie silniejszy formalnie,
- jezyk legacy nadal jest przydatny jako mapa translacyjna,
- ale sam `alpha_K_old` nie ma juz tak mocnego statusu jak `g0*` i `phi-FP`.

## 2. Audyt starych sciezek OP-3

Ponizej stan najwazniejszych historycznych sciezek selekcji `alpha_K`.

| Sciezka | Idea | Status | Wniosek |
|---|---|---|---|
| 1 | `argmin E*(alpha)` | zamknieta negatywnie | minimum bylo artefaktem `r_max`; brak fizycznego minimum |
| 2 | `a_c(alpha_K) = a_Gamma` | zamknieta negatywnie | numerycznie `a_c != a_Gamma` w punkcie fizycznym |
| 3 | `alpha_K/alpha_SN ~ n_s` | zamknieta negatywnie | koincydencja okazala sie artefaktem |
| 4 | `M_conf^B / M_conf^A = 206.77` | zamknieta negatywnie | brak naturalnego `p` dajacego fizyczny wynik |
| 5 | `argmax K2/K1` | zamknieta negatywnie | maksimum nie wypada przy `alpha_K ~= 8.56` |
| 6/7 | `E_kin^B / E_tot^B = 206.77` | zamknieta negatywnie | crossing istnieje, ale przy zlym `alpha`; nie odtwarza mas |
| 8+ | modyfikacje `V_mod` starego typu | brak domkniecia | nie daly stabilnego mechanizmu |

Wniosek audytowy:

- stare sciezki nie zamknely OP-3,
- problem nie wyglada juz jak "trzeba mocniej przeskanowac",
- potrzebny jest inny typ mostu niz czysty skan parametru.

## 3. Co pozostaje sensowne po audycie

Po odrzuceniu starych sciezek zostaja tylko dwa sensowne kierunki.

### Kierunek A: bezposredni

Nie szukac `alpha_K` jako "ladnej liczby" z krajobrazu ODE, tylko jako parametru emergentnego z relacji implicitnej:

`g0(K) + alpha * g1(K) = 0`

To jest jedyny kandydat, ktory:

- zachowuje lacznosc z analityka Koide,
- oddziela czesc geometryczna od czesci sprzezeniowej,
- daje szanse zapisac `alpha` jako funkcje obiektow nizzszej warstwy.

### Kierunek B: reframing

Uznac, ze sektor leptonowy jest juz kanonicznie kodowany przez:

`substrat -> g0^e -> phi-FP -> (r21, Q_K, m_tau)`

a `alpha_K_old` ma status:

`legacy recovered parameter`

czyli parametru odtwarzanego z nowego jezyka, ale niekoniecznie fundamentalnego.

## 4. Ocena bezposredniej sciezki `g0 + alpha g1 = 0`

To jest jedyna bezposrednia sciezka, ktora warto dalej rozwijac.

### Dlaczego akurat ta

W przeciwienstwie do dawnych prob:

- nie opiera sie na przypadkowym ekstremum,
- nie wymaga identyfikacji "fizycznego" crossing point po fakcie,
- nadaje sie do rozdzielenia na czesc substratowa i efektywna.

### Czego jeszcze brakuje

Zeby ta sciezka zamknela luke, trzeba jeszcze pokazac:

1. jaka jest kanoniczna postac `g0(K)` i `g1(K)`,
2. w jakiej formulacji rownanie jest wazne,
3. jak `alpha` zalezy od obiektow Warstwy I-II, a nie od `r21`,
4. jaki warunek wybiera dokladnie sektor leptonowy.

### Obecny werdykt

Ta sciezka jest:

- sensowna,
- technicznie obiecujaca,
- ale jeszcze nie domknieta.

Jej status na dzis:

`jedyna zywa sciezka bezposrednia`.

## 5. Ocena sciezki `g0* / phi-FP`

Ta sciezka jest obecnie najsilniejsza epistemicznie.

### Co juz daje

- bardzo mocne odtworzenie `r21`,
- formalne domkniecie mechanizmu `Q_K = 3/2`,
- bardzo mocne wyniki dla `m_tau`,
- sensowny kandydat na kanoniczny jezyk sektora leptonowego.

### Co jeszcze brakuje

Nie jest jeszcze zamkniety most:

`substrat -> g0^e`

oraz nie jest jeszcze spisana pelna mapa:

`alpha_K_old <-> (g0*, g0^e, phi-FP)`

### Obecny werdykt

Ta sciezka jest:

- mocniejsza od starego jezyka `alpha_K`,
- bardziej zgodna z najnowszym stanem repo,
- najlepszym kandydatem na kanoniczna formulacje sektora leptonowego.

## 6. Rozbicie OP-3 na dwa podproblemy

Najuczciwsze uporzadkowanie luki to:

### OP-3A

`Czy alpha_K_old wynika z Warstwy I-II bez uzycia r21?`

To jest problem bezposredniego wyprowadzenia.

### OP-3B

`Czy istnieje jawna i jednoznaczna mapa alpha_K_old <-> g0* / phi-FP?`

To jest problem formalnej supersesji starego jezyka przez nowy.

W praktyce:

- OP-3A jest nadal otwarty,
- OP-3B jest czesciowo zamkniety,
- i to OP-3B jest dzis bardziej realistyczna droga domkniecia luki.

## 7. Kryterium sukcesu

Luka `alpha_K` moze zostac uznana za zamknieta na dwa rozne sposoby.

### Wariant A: domkniecie silne

Mozna napisac:

`alpha_K ~= 8.56` wynika z Warstwy I-II bez uzycia `r21`.

Minimalne warunki:

1. istnieje jawny lancuch
   `Warstwa I-II -> efektywne rownanie solitonu -> alpha_K -> r21`
2. `r21` nie sluzy jako input kalibracyjny,
3. liczba `8.56` wynika z mechanizmu, nie z dopasowania po fakcie.

### Wariant B: domkniecie przez supersesje

Mozna napisac:

`alpha_K_old` nie jest juz parametrem fundamentalnym; staje sie parametrem recovered/legacy, a kanoniczny most biegnie przez `g0* / phi-FP`.

Minimalne warunki:

1. istnieje jawna mapa `alpha_K_old <-> g0* / phi-FP`,
2. nowy jezyk reprodukuje wszystkie kluczowe wyniki starego,
3. nowy jezyk ma mniejsza liczbe wejsc i czystszy podzial `Input / Derived`.

## 8. Rekomendacja robocza

Najbardziej realistyczna strategia jest dwuetapowa:

1. formalnie zamknac OP-3B,
2. utrzymac OP-3A jako program badawczy wysokiego ryzyka.

Powod:

- OP-3B porzadkuje dokumentacje i status teorii juz teraz,
- OP-3A pozostawia otwarta mozliwosc prawdziwego wyprowadzenia `alpha_K`,
- ale nie blokuje calego sektora leptonowego na parametrze, ktorego wszystkie dawne sciezki selekcji byly negatywne.

## 9. Minimalne nastepne kroki

Jesli celem jest dalsze domykanie luki, najkrotszy sensowny program to:

1. spisac jawny slownik przejsc `alpha_K_old <-> g0* / g0^e / phi-FP`,
2. wydzielic jedna kanoniczna postac rownania `g0 + alpha g1 = 0`,
3. sprawdzic, czy `alpha = F(a_Gamma, Phi_0, xi, g0*)` daje sie zapisac bez `r21`,
4. na koncu zdecydowac, czy `alpha_K` zostaje jako parametr fundamentalny, czy jako parametr legacy.

## Konkluzja

Na dzis najuczciwszy opis brzmi:

`alpha_K ~= 8.56` nie jest jeszcze wyprowadzone z pierwszych zasad.

Ale jednoczesnie nie jest juz jedynym nosnikiem sektora leptonowego.
Najsilniejszy obecnie opis biegnie przez `g0* / phi-FP`, a prawdziwa luka dotyczy nie tyle samej liczby `8.56`, ile relacji miedzy starym i nowym jezykiem teorii.
