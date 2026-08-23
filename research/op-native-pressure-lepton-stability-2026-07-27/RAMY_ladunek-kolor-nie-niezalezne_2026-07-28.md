# Ładunek i kolor NIE są niezależne — ramy, warunki i stan

**Data:** 2026-07-28
**Typ:** dokument-program (RAMY / CONSTRAINTS), nie dokument-wyniki
**Status uczciwy:** *nie mamy jeszcze nic wyprowadzonego.* To jest **scaffold** —
zbiór warunków, które będą musiały spełnić obiekty, które spróbujemy wyprowadzić.
**Autor sesji:** dialog + audyt (Claudian, tryb Opus)

---

## 0. Teza centralna (oś całego programu)

> **Ładunek i kolor nie są od siebie niezależne. Muszą wypływać z JEDNEJ
> struktury — topologii/kształtu własnej przestrzeni obiektu.**

To jest to, czego najbardziej brakowało. W Modelu Standardowym ta zależność
istnieje (lock color↔charge siedzi w ilorazie cechowania ℤ₃), ale jest
**postulowana** przez grupę cechowania. Celem TGP-natywnym jest **wyprowadzić
tę zależność z kształtu**, nie wsadzić jej ręcznie.

Wszystko poniżej to warunki brzegowe dla tego celu.

---

## 1. Czego się pozbyliśmy w tej sesji (ślepe zaułki — nie powtarzać)

Uczciwy rejestr obalonych/fittowanych rzeczy, żeby nie wracać:

| Co | Werdykt | Dowód |
|---|---|---|
| **Bounce-hierarchy** (Tier-2 „native mechanism") | **artefakt** | N_neg skaluje z pudłem R (12/19/25 dla R=40/60/80), próżnia F-S ma tyle samo |
| Solitony crown w kanonicznej F-A | **nie istnieją** | ODE ucieka do ∞ dla wszystkich g₀ (runaway); F-A ma stabilną próżnię W''(1)=+1 |
| **Bariera = horyzont metryki** (why_n3/M9.1'') | **cyrkularne** | most ψ=0.3814g+0.6186 to fit 2-punktowy, gdzie „coincidence" jest jednym z warunków; dokument sam pisze „(definitionally)" |
| Multiplikatywności z twistów (Tier-A) | **nie pasują** do SM | Q=±⅓ daje 6, nie 3 |
| Parametr-free „3" stabilnych warkoczy | **realne, ale zła struktura** | 3 klasy = 1 order-2 + 2 order-3 (różne QN), nie 3 identyczne generacje |
| Lock color↔charge z gołych warkoczy | **no-go** | abelianizacja B₃ = ℤ (writhe); nieabelowa triality nie da się śledzić abelowym ładunkiem |

**Wzorzec (zgodny ze znaną fizyką):** struktura SM (locki, multiplikatywności,
generacje) **nie wypada z kinematyki za darmo** — wymaga warunków spójności.
To jest treść twierdzeń no-go dla preonów ('t Hooft anomaly matching).
Te więzy **SĄ** fizyką.

---

## 2. Jedyny twardy filar (wzorzec do naśladowania)

**Spin-½ wyprowadzony z topologii kształtu — bez zewnętrznego pola.** To jest
jedyna rzecz w tym obszarze, która jest realnie DERIVED (i to podwójnie):

- **Ścieżka A (Q5):** vielbein hedgehog → mapa S³→SU(2)~S³, π₃(S³)=ℤ,
  nawinięcie **B=1** (topologiczne, niezależne od profilu). Kwantyzacja
  Finkelsteina-Rubinsteina: obrót 2π → (−1)^B = −1 ⟹ **fermion, spin-½**.
  Bonus: **g=2** z „prąd EM sprzęga się z prądem topologicznym".
- **Ścieżka B (why_n3 Faza 3):** RP²=S²/ℤ₂, faza Berry'ego γ=π pod 2π ⟹
  Ψ(2π)=−Ψ ⟹ **spin-½**.

Dwie różne topologie → ten sam wynik = sygnatura realnej struktury.

**To jest wzorzec:** sztywna liczba kwantowa z topologii własnej przestrzeni
obiektu, zero zewnętrznego pola, zero dostrajania. **Ładunek musi wyjść
analogicznie.**

---

## 3. Wyłaniające się ramy (obraz jakościowy — warunki, nie wyniki)

Ontologia (potwierdzona w dialogu):

- **Obiekt = kształt** (soliton topologiczny) generujący własną przestrzeń
  („balon z ziarnem piasku": ziarno = soliton, balon = jego przestrzeń).
- **Tło = zbiorowa przestrzeń innych obiektów** (NIE osobne pole, NIE substrat
  jako pole). Każdy obiekt porusza się przez przestrzeń zrobioną przez resztę.
- **Brak zewnętrznego pola EM na wejściu.** EM jest *obserwacją*, nie tłem.

Warstwy (co czym jest):

| Wielkość | Czym jest w ramach | Status |
|---|---|---|
| **kolor** | struktura/kształt (permutacja/nawinięcie); triality = ℤ₃ = A₃ na 3 pasmach | koncept |
| **spin** | nawinięcie B | ✅ **wyprowadzone** |
| **ładunek** | **selektywny przekaz pędu kształtu**, ujawniany tylko przez pomiar zmieniający pęd | operacyjnie zdefiniowane |
| **EM** | emergentne z rezonansu kształtów (analog skyrmion emergent-EM: ruch tekstury → emergentne E, nawinięcie → emergentne B, siła Lorentza) | koncept, ugruntowany w literaturze |
| **{+, −, 0}** | znak nawinięcia / superpozycja (0 ≈ nałożenie +/−) | hipoteza |

Kluczowa rzecz o naturze ładunku (doprecyzowanie z dialogu):
- Ładunek jest **intrinsyczny dla kształtu** — istnieje też w ruchu jednostajnym.
- Ruch jednostajny daje **pole** (bezstratne, odwracalne).
- Ale **pomiar wymaga zmiany pędu** (interferencja). Bez zmiany pędu nie ma
  odczytu. To jest **operacyjna definicja ładunku = odpowiedź kształtu na
  przekaz pędu** (form factor / amplituda rozpraszania).

---

## 4. Falsyfikator: eksperyment e–p (operacyjna definicja ładunku)

```
Stan:   e(−) z Pe,  p(+) z Pp,  reszta = 0,  wszystkie mają jakiś pęd
Akcja:  dodaj Pd do e, przepchnij e przez cały basen
Predykcja:
   e_końcowe = (Pe + Pd) ∓ δ
   p_końcowe = Pp ± δ
   0-balony: pęd NETTO niezmieniony (ruszają się przelotnie, ale oddają — brak rezonansu)
   zachowanie pędu: δ przekazane WYŁĄCZNIE e ↔ p
```

To jest rozpraszanie Coulomba jako operacyjna definicja ładunku, **bez
postulowania EM**. Znak i wielkość δ = ładunek. Neutralne = nie rozpraszają.

**Warunek nie-cyrkularności:** selektywność (tylko e↔p) i znak δ muszą
**WYPAŚĆ z dynamiki kształtów**, nie być wstawione ręcznie. Inaczej — piąty fitting.

---

## 5. WARUNKI, które musi spełnić każdy wyprowadzony obiekt (checklista)

To jest sedno tego dokumentu — ramy dla przyszłej derywacji:

1. **Metastabilny** — z topologii/reguł, NIE z numerycznego hacka (jak „bounces").
   Metastabilność = bariera kombinatoryczna (niezmienniczość pod ruchami
   substratu), **nie** basen atraktora (bo singletony to saddle, nie minima).
2. **Symetryczny/składalny** — grupuje się w większe struktury (reguły dopasowania).
3. **Minimalny** — najmniej reguł/generatorów.
4. **Samoreferencyjny / background-free** — obiekt współdefiniuje swoją przestrzeń;
   żadnego zewnętrznego pola jako wejścia. Cząstka = punkt stały: struktura
   zgodna z geometrią, którą sama generuje.
5. **⟵ CENTRALNE: kolor i ładunek współ-wyprowadzone z JEDNEGO kształtu.**
   Nie mogą być niezależnymi stopniami swobody. Ich zależność (lock) ma
   **wypaść**, nie być postulowana.
6. **Spin z nawinięcia** — już spełnione przez wzorzec (§2); nowy obiekt musi
   być z nim zgodny.
7. **Ładunek = selektywny przekaz pędu** — sygnatura eksperymentu e–p (§4);
   skwantowany {+,−,0} ma pochodzić z topologii (jak B przy spinie), nie z
   ciągłego widma.
8. **EM emergentne, nie wejściowe** — ma wyjść jako rezonans/sprzężenie ruchu
   kształtów, z Coulombowską strukturą (znak: + i − się przyciągają, ++ odpycha).
9. **Sektor:** naturalny dom to **kwarki/kolor** (ładunek w trzecich + ℤ₃), NIE
   leptony (całkowity ładunek, brak koloru). Leptony były dopasowaniem
   odziedziczonym z fittowanej pracy why_n3 (bo pasują do Koide).
10. **„3":** pojawia się jako **kolor** (3 pasma = 3 kolory, sztywne z PSL(2,ℤ)),
    ale **liczba generacji pozostaje OTWARTA** (nie mylić tych dwóch trójek).

---

## 6. Co jest naprawdę otwarte (uczciwie: prawie wszystko)

- Jakie **konkretnie kształty** = +, −, 0? (kandydat: przeciwne nawinięcie / superpozycja)
- Czy lock **kolor↔ładunek WYPADA** z kształtu (§0)? — **nierozstrzygnięte, to główny cel**
- Czy EM emerguje z **Coulombowską** strukturą (nie tylko „jakieś sprzężenie")?
- Skąd **3 generacje** (osobne od 3 kolorów)?
- Minimalny zestaw **ruchów substratu Φ** (odkładany — to decyzja teoretyczna użytkownika)
- Kryterium **rezonansu** przygwożdżone do konkretnego warunku (dyspersja/faza)

---

## 7. Następny projekt (zscope'owany)

**Wierna symulacja rozpraszania e–p** (nie szybki toy — bo toy = szósty fit/artefakt):
1. Wybrać teorię pola z solitonami przeciwnego nawinięcia (+/−) i neutralnymi (0).
2. Bezstratna (hamiltonowska) dynamika nakładających się kształtów.
3. Puścić e z Pd przez „basen", zmierzyć pędy końcowe wszystkich.
4. Test: czy δ przekazuje się **selektywnie** e↔p, ze znakiem Coulomba, z
   neutralnymi transparentnymi — i czy to **wypada**, a nie jest wstawione.

Arena: **solitony topologiczne (skyrmiony)** — czysta warstwa, gdzie spin/B=1/g=2
już działają. NIE profile F-S (obalone).

---

## 8. Jednozdaniowe podsumowanie stanu

> Nie mamy jeszcze żadnego wyprowadzonego obiektu. Mamy: jeden twardy filar
> (spin z nawinięcia), rejestr obalonych ślepych zaułków, i spójny zestaw
> warunków — z których centralny brzmi: **ładunek i kolor muszą wypłynąć
> z jednej topologii kształtu, nie być niezależne.** To jest scaffold pod
> derywację, nie derywacja.

---

**Powiązane pliki (ta sesja):**
- `AUDYT_KRYTYCZNY_2026-07-28.md` (obalenie bounce-hierarchy + test F-A)
- `TIERA_braid_bruteforce.py`, `TIERB_stable_braids.py`, `TIERB_color_charge_lock.py`
- `AUDIT_diagnostic.py`, `AUDIT_FA_dual_test.py`

**Powiązane (rdzeń/istniejące):**
- `research/why_n3/PHASE3_RP2_defect_quantization.md` (spin-½, wzorzec)
- `research/qm_spin/README.md` (spin z π₃(S³), g=2)
- `core/sek08c_metryka_z_substratu.tex` (metryka z budżetu — z zastrzeżeniami)
