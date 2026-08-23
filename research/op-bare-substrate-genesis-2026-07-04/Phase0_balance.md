---
title: "Phase0_balance — LOCK: geneza samopodtrzymujacych solitonow z golego substratu"
date: 2026-07-04
type: phase0-lock
tgp_owner: research/op-bare-substrate-genesis-2026-07-04
status: LOCKED
anti_lakatos_lock: PRESERVED
tags: [substrate-genesis, soliton, lock, self-sustaining, nucleation, gravity-bridge]
related:
  - "[[../op-blocked-soliton-bang-2026-07-04/README.md]]"
  - "[[../op-blocked-soliton-bang-2026-07-04/POST_CLOSE_BRAINSTORM_UPDATE.md]]"
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
  - "[[../../core/sek01_ontologia/sek01_ontologia.tex]]"
  - "[[../../axioms/roznica_N0/dodatek0_aksjomatyka_roznicy.tex]]"
---

# Phase 0 — LOCK: Geneza samopodtrzymujacych solitonow z golego substratu

## 0. Hipoteza robocza

Ten cykl testuje poprawiona hipoteze autora (2026-07-04, po zamknieciu
`op-blocked-soliton-bang`):

> Gola faza substratu (Phi = 0) jest niemetryczna i strukturalnie
> nietrwala. Male zaburzenie s_i != 0 wlacza lokalny zalazek metrycznosci
> (Phi_i = s_i^2 > 0). Pojedyncze zaburzenie samo sie relaksuje (rozlewa /
> wraca do Phi = 0 — to jego wlasna "kreacja przestrzeni" / entropia).
> Ale odpowiednio gesta konfiguracja nakladajacych sie zaburzen podtrzymuje
> sie NAWZAJEM: wchodzi w lock, z ktorego nie da sie zrelaksowac. Jesli tak,
> ten samopodtrzymujacy lock jest kandydatem na fundament mostu do
> grawitacji. Jesli nie — trzeba szukac dalej.

Kluczowa poprawka wzgledem poprzedniego cyklu: **przezywalnosc struktury
musi wynikac z jednej dynamiki pola substratu, a NIE z osobnej reguly
amplitud/pozycji odsprzezonej od pola.** W poprzednim opie `overlap_i`
liczyl tylko geometrie sasiadow i pole `phi` nie wchodzilo do losu
solitonu — to dekretowanie wyniku, nie jego wyprowadzenie. Tego bledu ten
cykl nie powtarza.

Wersja minimalna NIE twierdzi, ze policzono grawitacje. Testuje tylko, czy
substrat potrafi SAM (bez wrzucania gotowych solitonow) wygenerowac
lokalna, samopodtrzymujaca sie strukture Phi > 0 dzieki kolektywnemu
lockowi.

## 1. Model zalockowany przed kodem

### 1.1 Substrat

Substrat `Gamma` reprezentujemy jako siatke 2D (graf regularny `V,E`).
Na kazdym wezle jedna zmienna substratowa `s(x)` w R. Definiujemy:

```text
Phi(x) = s(x)^2 >= 0
```

Interpretacja zgodna z sek01_ontologia:

- `Phi = 0` = faza niemetryczna substratu (nie "pusta przestrzen").
- `Phi > 0` = faza metryczna (propagacja / geometria maja sens).

Symetria `s -> -s` (Z2) jest zachowana; obserwablem fizycznym jest `Phi`,
nie znak `s`. Zgodne z lancuchem `N0 -> Z2 -> s_i -> Gamma -> Phi -> g`.

### 1.2 Dynamika (jedna, spojna, bez osobnej reguly przezycia)

Uzywamy przeplywu gradientowego (relaksacji) w parametrze selekcji `tau`:

```text
ds/dtau = - dH_Gamma / ds
        = kappa * Lap(s) - V'(s)
```

gdzie `H_Gamma[s] = sum_x [ 0.5 * kappa * |grad s|^2 + V(s) ]`.

`tau` jest parametrem relaksacji/selekcji, NIE czasem fizycznym: przed
faza metryczna czas fizyczny nie jest zdefiniowany (zgodne z
dodatek0_aksjomatyka_roznicy: niestabilnosc `N0/Phi=0` jest
strukturalna/logiczna/miarowa, nie procesem w czasie).

### 1.3 Potencjal: Phi = 0 metastabilne, Phi > 0 glebsze, z bariera

Zeby "single decays / many lock" moglo sie WYLONIC, a nie byc wpisane
regula, potencjal ma miec **metastabilna falszywa prozne w s = 0** i
glebsza prawdziwa prozne w `s = +-s*`, oddzielone bariera (przejscie
pierwszego rodzaju / teoria nukleacji):

```text
V(s) = 0.5 * a * s^2  -  (b/3) * |s|^3  +  0.25 * c * s^4
V'(s) = a*s - b*|s|*s + c*s^3 = s * (a - b*|s| + c*s^2)
```

Zalockowane wspolczynniki (Z2-symetryczne przez |s|):

```text
a = 0.50   (dodatnia masa w s=0  -> s=0 to LOKALNE minimum, metastabilne)
b = 1.60
c = 1.00
```

Konsekwencje (policzone analitycznie przed kodem):

- `s = 0` jest lokalnym minimum (V''(0) = a > 0) -> gola faza `Phi = 0`
  jest metastabilna, nie rozpada sie samorzutnie bez seedu.
- ekstrema niezerowe: `|s| = (b +- sqrt(b^2 - 4ac)) / (2c)`,
  `b^2 - 4ac = 0.56 > 0` -> istnieja dwa dodatnie pierwiastki:
  bariera przy `|s| ~ 0.426` (Phi_bar ~ 0.181),
  prawdziwa proznia przy `s* ~ 1.174` (Phi* ~ 1.378).
- `V(s*) ~ -0.044 < 0 = V(0)` -> faza metryczna jest energetycznie glebsza.

Fizyka lock: podkrytyczna kropla `Phi > 0` kurczy sie (relaksuje do
`Phi = 0`); nadkrytyczna nukleuje i utrzymuje sie. Klaster podkrytycznych
seedow moze polaczyc sie w nadkrytyczny obszar -> **kolektywny lock**.
Ten mechanizm wychodzi z samego `H_Gamma`, bez mobility-clampu i bez
reguly amplitud.

### 1.4 Brak ukrytych podporek

- BRAK osobnej ewolucji amplitud/pozycji. Los kazdej struktury czytamy
  wylacznie z pola `s`.
- BRAK mobility zerujacej relaksacje (to byl clamp-ryzykowny chwyt).
  Lock ma pochodzic z prawdziwej prozni `V`, nie z zamrozenia solvera.
- BRAK reguly nukleacji/spawn. Nowe struktury moga powstawac tylko przez
  dynamike pola.

## 2. Parametry Phase 1 (LOCKED)

Jednostki bezwymiarowe:

```text
grid:        N = 128,  L = 64.0   (dx = 0.5)
kappa = 0.50
a = 0.50, b = 1.60, c = 1.00
dt = 0.02     (stabilnosc jawnego Eulera: dt < dx^2/(4*kappa) = 0.125)
steps_decay = 6000     (dlugi run testu zaniku/locku)
steps_probe = 3000     (krotszy run pomocniczy)
seed profil: s(x) = A0 * exp(-r^2 / (2*w^2)),  w = 1.50, dodawany do s=0
seed amplitude scan: A0 in {0.4, 0.6, 0.8, 1.0, 1.2, 1.4}
cluster: K = 7 seedow (heks pierscien + centrum), separacje D in {2, 3, 4}
seed smoke (dla stochastycznych wariantow): 20260704
Phi_metric_threshold: eps = 0.30   (wyraznie powyzej bariery Phi_bar~0.18)
metric_area = frakcja wezlow z Phi > eps
A_min = minimalna metric_area uznawana za utrzymana strukture = 4 wezly / N^2

--- doprecyzowania LOCK (poprawka pre-code, patrz sekcja 11) ---
brzegi:      periodyczne (torus), Laplasjan 5-punktowy
bare noise:  s0(x) ~ U(-0.05, +0.05) iid na kazdym wezle
             (amplituda 0.05 << bariera |s_bar| ~ 0.426), seed 20260704
tau_mid:     krok 3000 = steps_decay/2 (mianownik persistence)
persistence_tail = metric_area(krok 6000) / metric_area(krok 5400)
             (ostatnie 10% runu; odroznia lock od powolnego zaniku)
delta_bump:  gaussowski, A_bump = 0.10, w_bump = 1.0, dodawany do
             biezacego pola s (nie do Phi)
localized:   zaden wezel w odleglosci Chebysheva > 0.4*L od centrum
             siatki nie ma Phi > eps; jesli front dotknie tej strefy,
             run dostaje flage boundary_contact i dalsze kroki nie
             wchodza do ocen G1-G6
```

Uwaga metodologiczna: jesli jawny Euler przy tych `dt/dx` okaze sie
numerycznie niestabilny, dozwolony jest JEDEN techniczny smoke-run
kalibracyjny obnizajacy `dt` (stabilnosc numeryczna != strojenie pod wynik).
Kryteria G1-G6 i procedura pozostaja niezmienione.

## 3. Scenariusze

1. `bare`: brak seedu (albo szum o amplitudzie << bariera). Test, czy gola
   faza pozostaje niemetryczna.
2. `single(A0)`: jeden seed w centrum, skan `A0`. Wyznacza operacyjnie
   prog krytyczny `A_c` (najwiekszy `A0`, przy ktorym pojedynczy seed
   zanika).
3. `cluster(A0, D)`: K = 7 seedow o TAKIM SAMYM `A0` jak podkrytyczny
   single, separacja `D`. Test kolektywnego locku.
4. `create_cost`: pomiar lokalnego kosztu utworzenia dodatkowego malego
   zaburzenia `DeltaE_create(x) = H[s + delta_bump@x] - H[s]` w regionie
   `Phi = 0` oraz w regionie zablokowanym `Phi > 0`.

## 4. Obserwable

```text
- Phi_max(tau), Phi_mean(tau)
- metric_area(tau) = frakcja wezlow z Phi > eps
- localized: czy region metryczny ma skonczony nosnik (Phi -> 0 przy brzegu)
- persistence = metric_area(tau_final) / metric_area(tau_mid)
- A_c = prog krytyczny amplitudy dla pojedynczego seedu
- DeltaE_create w Phi=0  vs  DeltaE_create w Phi>0
- H_Gamma(tau) (monotoniczny spadek = poprawnosc przeplywu gradientowego)
- finite / runaway flags
```

## 5. Kryteria sukcesu i falsyfikacji (LOCKED przed kodem)

### G1 — deterministycznosc i skonczonosc (techniczne)
Ten sam seed i parametry daja identyczny wynik. Brak NaN/Inf, brak runaway,
`H_Gamma(tau)` nierosnace (z dokladnoscia numeryczna). = PASS techniczny.

### G2 — gola faza pozostaje niemetryczna
Scenariusz `bare` po `steps_decay`: `metric_area < A_min`. Potwierdza, ze
`Phi = 0` jest metastabilne (nie wypelnia sie samorzutnie).
FALSYFIKATOR: gola faza sama nukleuje bez seedu -> `Phi = 0` nie jest
metastabilne, model zly.

### G3 — pojedynczy podkrytyczny seed zanika
Istnieje `A0` (podkrytyczny), przy ktorym `single(A0)` po `steps_decay`
ma `metric_area < A_min` (Phi_max < eps). "single decays".
FALSYFIKATOR: kazdy `A0` w skanie albo od razu lockuje, albo eksploduje —
brak reżimu, w ktorym pojedynczy seed realnie zanika.

### G4 — kolektywny lock (RDZEN HIPOTEZY)
Dla NAJWIEKSZEGO podkrytycznego `A0` z G3, scenariusz `cluster(A0, D)`
przy co najmniej jednej separacji `D` po `steps_decay` ma:

```text
metric_area >= A_min
ORAZ persistence >= 0.9          (utrzymuje sie lub rosnie)
ORAZ persistence_tail >= 0.99    (ogon nie maleje: lock, nie powolny zanik)
```

czyli klaster tych samych seedow, ktore pojedynczo zanikaja, podtrzymuje
sie nawzajem. "many lock".

Obowiazkowa flaga interpretacyjna (superpozycja vs oddzialywanie):
dla kazdego `cluster(A0, D)` raportowac `Phi_max(tau=0)`. Jesli poczatkowa
amplituda centralna klastra `s_max(0)` przekracza najmniejsza NADkrytyczna
amplitude pojedynczego seedu z G3, klaster oznaczyc `superposition-trivial`
(lock wynika z liniowego zsumowania seedow w jeden nadkrytyczny blob od
tau=0, np. przy D=2 i w=1.5 centrum dostaje ~3.5*A0). Lock uznany za
KOLEKTYWNY w mocnym sensie tylko, gdy zachodzi rowniez dla przypadku bez
tej flagi. Oba warianty raportowac wprost — flaga nie zmienia PASS/FAIL
G4, zmienia jego interpretacje w raporcie koncowym.
FALSYFIKATOR: klaster zanika tak samo jak single -> **substrat NIE
generuje samopodtrzymujacych solitonow przez kolektywny lock** -> hipoteza
falsyfikowana, trzeba szukac dalej.

### G5 — lock jest prawdziwa bifurkacja, nie artefaktem
Skan pokazuje OSTRY prog (po `A0` lub `D`) rozdzielajacy zanik od locku,
i lock utrzymuje sie przy PELNEJ relaksacji (brak clampu, `H_Gamma` nadal
spada). = istnieje krytyczna gestosc/rozmiar.
FALSYFIKATOR: lock pojawia sie tylko po wprowadzeniu jawnego clampu /
zamrozenia -> artefakt algorytmiczny, nie fizyka.

Kontrola pinningu siatki (zalockowana pre-code): szerokosc scianki
domenowej to ~sqrt(2*kappa/V''(s*)) ~ 1.07 ~ 2*dx, wiec front moze sie
"przypiac" do dyskretnej siatki i udawac lock. Dla przypadkow granicznych
G4/G5 (najmniejsze A0/D dajace lock oraz najwieksze dajace zanik)
dozwolony jest JEDEN kontrolny re-run na drobniejszej siatce:
`N = 256, L = 64 (dx = 0.25), dt = 0.005`, te same parametry fizyczne
i ten sam tau_final. Wymog: klasyfikacja lock/zanik musi sie jakosciowo
zgadzac z siatka bazowa. FALSYFIKATOR pinningu: lock znika na drobniejszej
siatce -> artefakt dyskretyzacji, G4/G5 dla tego punktu = FAIL.

### G6 — koszt kreacji spada w fazie metrycznej (kluczowe wg autora)

Pomiar `DeltaE_create(x) = H[s + delta_bump@x] - H[s]` (bump z sekcji 2)
wykonywany w TRZECH zalockowanych lokalizacjach na skonfigurowanym stanie
z lockiem:

```text
(i)   FAR:   gleboko w Phi = 0, odleglosc > 10*w od struktury
(ii)  CORE:  wnetrze regionu zablokowanego (s ~ s*), maksimum Phi
(iii) FRONT: interfejs regionu (wezel o Phi najblizszym eps na
             zewnetrznym zboczu struktury)
```

Predykcja analityczna zapisana PRZED kodem: dla malego bumpa w CORE
`DeltaE ~ 0.5*V''(s*)*delta^2 > 0` z `V''(s*) ~ 0.88 > V''(0) = a = 0.5`,
wiec CORE moze byc DROZSZY niz FAR — to NIE falsyfikuje G6. Fizycznym
nosnikiem tezy "kreacja na granicy metryki" jest FRONT, gdzie poszerzenie
domeny zamienia falszywa proznie (V=0) na prawdziwa (V(s*)<0).

Kryterium PASS:

```text
DeltaE_create(FRONT)  <  DeltaE_create(FAR)
```

i najlepiej `DeltaE_create(FRONT) <= 0` (kreacja kolejnego zaburzenia na
granicy regionu metrycznego jest bezkosztowa/samoczynna).
Wartosc CORE raportowac wprost obok FRONT i FAR (bez wplywu na PASS/FAIL).
FALSYFIKATOR: koszt kreacji na FRONT jest taki sam lub wyzszy niz w FAR ->
brak "bezkosztowej kreacji", most do "wielki wybuch trwa na granicy
metryki" nie ma tu podstawy.

## 6. Regula decyzyjna calego cyklu

```text
G2 + G3 + G4 + G5 + G6 PASS
  -> substrat POTRAFI wygenerowac samopodtrzymujace, lockowane struktury
     Phi>0 -> fundament mostu do grawitacji STOI, przejdz do nastepnego
     cyklu (uwiezienie powlokowe / dynamika scian / most do g_mu_nu).

G4 FAIL (klaster zanika jak single)
  -> substrat w tym minimalnym modelu NIE podtrzymuje sie kolektywnie
     -> hipoteza falsyfikowana w tej klasie modeli, szukaj dalej
     (inny H_Gamma, inne sprzezenie, wiecej pol).

Wynik mieszany raportowac wprost, bez strojenia progow.
```

## 7. Forbidden moves

1. Nie startowac od gotowych solitonow. Wejscie to `Gamma + s(x)`, koniec.
2. Nie dodawac osobnej reguly amplitud/pozycji ani spawn-rule. Los
   struktur czytany wylacznie z pola `s`.
3. Nie uzywac mobility -> 0 ani clampu jako mechanizmu locku. Clamp tylko
   jako guard przeciw overflow; dobicie do guardu = run oznaczony runaway.
4. Nie zmieniac progow G1-G6 ani wspolczynnikow `a,b,c` po pierwszym runie
   (poza jednym technicznym obnizeniem `dt` dla stabilnosci oraz jednym
   kontrolnym re-runem na drobniejszej siatce z sekcji G5 — oba to guardy
   numeryczne, nie strojenie pod wynik).
5. Nie interpretowac wyniku jako dowodu GR, MOND ani kosmologii.
6. Nie promowac wyniku do `core/` ani rejestrow bez osobnej autoryzacji.
7. Raportowac wyniki negatywne wprost. G4 FAIL jest pelnoprawnym wynikiem
   cyklu ("szukaj dalej"), nie porazka do ukrycia.
8. Nie mylic `tau` z czasem fizycznym w interpretacji.

## 8. Czego ten cykl NIE robi

- Nie stabilizuje solitonu do skonczonego rozmiaru: nadkrytyczny obszar
  moze rosnac (front domenowy). Rozroznienie "lokalny soliton vs rosnaca
  domena/sciana" jest osobnym pytaniem (por. [[../op-wall-dynamics-2026-07-03/README.md]]).
  Tu mierzymy tylko: czy struktura sie PODTRZYMUJE i czy pochodzi z locku.
- Nie liczy metryki `g_mu_nu` ani geodezyjnych.
- Nie liczy uwiezienia powlokowego (sila przywracajaca na srodku) — to
  nastepny cykl, jesli G4 PASS.
- Nie dotyka `core/` ani rejestrow predykcji.

## 9. Deliverables

- `Phase0_balance.md` (ten plik)
- `Phase0_output.txt` (numeryczna weryfikacja rachunkow analitycznych z sekcji 1.3, pre-code)
- `Phase1_substrate_genesis.py`
- `Phase1_output.txt`
- `Phase2_threshold_and_create_cost.py`  (skan progu + DeltaE_create)
- `Phase2_output.txt`
- `Phase_FINAL_close.md`
- `README.md`

## 10. Deklaracja LOCK

Model (`Gamma`, `s`, `Phi = s^2`, przeplyw gradientowy w `tau`, potencjal z
metastabilnym `Phi = 0`), wspolczynniki `a,b,c`, scenariusze i kryteria
G1-G6 oraz regula decyzyjna zapisano PRZED napisaniem kodu Phase 1. Wynik
dodatni, mieszany albo negatywny bedzie raportowany bez zmiany progow.
Rdzen: przezywalnosc czytana z jednej dynamiki pola, lock z prawdziwej
prozni potencjalu, zero osobnej reguly przezycia.

## 11. Log poprawek LOCK (wszystkie PRE-CODE, 2026-07-04)

Ponizsze poprawki naniesiono po review planu, ale PRZED napisaniem
pierwszej linijki kodu Phase 1 — anti-Lakatos zachowany (doprecyzowanie
locka przed eksperymentem, nie strojenie po wyniku):

1. **Sekcja 2**: zalockowano brakujace elementy: warunki brzegowe
   (periodyczne, torus), stencil Laplasjanu (5-punktowy), amplituda szumu
   `bare` (U(-0.05, +0.05)), definicja `tau_mid` (krok 3000), parametry
   `delta_bump` (A_bump = 0.10, w_bump = 1.0), operacyjna definicja
   `localized` + flaga `boundary_contact`.
2. **G4**: dodano `persistence_tail >= 0.99` (ogon ostatnich 10% runu nie
   maleje — odroznia lock od powolnego zaniku) oraz obowiazkowa flage
   `superposition-trivial` (lock przez liniowe zsumowanie seedow w jeden
   nadkrytyczny blob od tau=0 vs lock przez oddzialywanie frontow).
   Flaga nie zmienia PASS/FAIL, zmienia interpretacje.
3. **G5**: dodano kontrole pinningu siatki (scianka ~2*dx moze sie
   przypiac do dyskretnej siatki i udawac lock): jeden kontrolny re-run
   przypadkow granicznych na `N=256, dx=0.25, dt=0.005`.
4. **G6**: doprecyzowano lokalizacje pomiaru (FAR / CORE / FRONT) i
   przeniesiono kryterium PASS na FRONT, z zapisana przed kodem predykcja
   analityczna, ze CORE moze byc drozszy od FAR (`V''(s*) ~ 0.88 >
   V''(0) = 0.5`) — pierwotna definicja ("w regionie Phi>0", bez
   lokalizacji) FAILowalaby z powodow trywialnych, nie fizycznych.
5. **Forbidden move 4**: rozszerzono wyjatek techniczny o kontrolny
   re-run z G5.
6. **Faza 0**: rachunki analityczne z sekcji 1.3 zweryfikowano
   numerycznie przed kodem (pierwiastki, V(s*), krzywizny, szerokosc
   scianki, ograniczenia stabilnosci) — wynik w `Phase0_output.txt`.
