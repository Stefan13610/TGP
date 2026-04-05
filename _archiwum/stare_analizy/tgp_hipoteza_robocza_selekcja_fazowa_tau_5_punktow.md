# Hipoteza robocza: selekcja fazowa leptonów — 5 punktów do dalszej analizy

**Status:** robocze, do sprawdzenia  
**Cel:** uchwycić, gdzie w obecnym formalizmie może siedzieć brakująca informacja dla `tau`, oraz wskazać miejsca, w których analiza mogła wejść w ślepy zaułek.

---

## 1. Być może za wcześnie zredukowałem rodzinę rozwiązań do jednego parametru `g_0`

Obecna analiza w praktyce traktuje stan leptonowy jako punkt na osi `g_0`, a następnie czyta z niego amplitudę ogona `A_tail`. To mogło wystarczyć dla elektronu i mionu, ale być zbyt ubogie dla trzeciej generacji.

**Hipoteza robocza:** selekcja generacji nie zachodzi w przestrzeni jednowymiarowej `A_tail(g_0)`, tylko w przestrzeni co najmniej dwuwymiarowej, np. `(A_tail, faza)` albo `(B_tail, C_tail)`.

**Co mogłem przeoczyć:**
- dodatkowy stopień swobody ukryty w asymptotyce,
- rozdzielenie stanów przez fazę, a nie samą amplitudę,
- fakt, że `tau` nie jest „trzecim takim samym punktem”, tylko stanem innego typu.

**Gdzie mogłem wejść w kozi róg:**
- próbując na siłę wydobyć pełną triadę `e, μ, τ` z jednej funkcji skalarnej,
- zakładając, że jedna oś parametryzacji musi wystarczyć dla wszystkich generacji.

---

## 2. Faza ogona wygląda na drugi inwariant i to właśnie tam może siedzieć selekcja

Dla asymptotyki

```math
(g(r)-1)\,r \approx B_{\rm tail}\cos r + C_{\rm tail}\sin r
```

naturalnym rozszerzeniem amplitudy

```math
A_{\rm tail}=\sqrt{B_{\rm tail}^2 + C_{\rm tail}^2}
```

jest faza zakodowana przez relację między `B_tail` i `C_tail`.

**Hipoteza robocza:** generacje są wybierane przez wyróżnione warunki fazowe, np. punkty typu `B_tail = C_tail` albo `C_tail = 0`, a nie tylko przez samą wartość `A_tail`.

**Co mogłem przeoczyć:**
- że `μ` i ciężki stan `tau`-podobny są bardziej naturalnie opisane przez klasę fazową niż przez samą masę,
- że Koide może być cieniem geometrii fazowej rodziny rozwiązań.

**Gdzie mogłem wejść w kozi róg:**
- traktując `B_tail` i `C_tail` wyłącznie jako pomocniczy fit techniczny,
- koncentrując się tylko na `A_tail^4`, mimo że faza może wybierać klasę stanu.

---

## 3. Problem z `tau` może nie dotyczyć selekcji stanu, tylko mapy `stan -> masa`

Wyniki sugerują, że ciężki stan jest widoczny strukturalnie, ale jego masa nadal nie doskakuje do fizycznego `m_τ` w gołym modelu `m \propto A_tail^4`.

**Hipoteza robocza:** obecne ODE może poprawnie wskazywać typ stanu `tau`, ale masa wymaga jeszcze dodatkowego czynnika dynamicznego, np.

```math
m \propto A_{\rm tail}^4 \cdot F(B_{\rm tail}, C_{\rm tail})
```

lub innej poprawki zależnej od fazy, rdzenia albo stabilności.

**Co mogłem przeoczyć:**
- że reguła `A_tail^4` jest tylko pierwszym przybliżeniem,
- że ciężki stan wymaga korekty fazowej, topologicznej albo energetycznej,
- że brakujące ~kilka procent / rząd skali nie musi obalać selekcji, tylko wskazywać brakujący człon w mapie masowej.

**Gdzie mogłem wejść w kozi róg:**
- utożsamiając niepowodzenie prostej mapy masowej z nieistnieniem samego stanu `tau`,
- próbując ratować `tau` przez relacje typu `g_0^τ \approx 2 g_0^e`, zamiast sprawdzić, czy błąd siedzi w formule na masę.

---

## 4. Koide może być warunkiem globalnym dla trypletu, a nie lokalną własnością pojedynczego rozwiązania

Najbardziej uporczywy problem polega na tym, że obecne jednosolitonowe ODE nie generuje jeszcze czysto fizycznego Koidego dla `τ`, mimo że sama algebra Koidego bardzo dobrze domyka tryplet.

**Hipoteza robocza:** warunek `Q = 3/2` nie wynika z pojedynczego profilu `g(r)`, lecz z relacji między trzema dozwolonymi stanami w rodzinie rozwiązań.

**Co mogłem przeoczyć:**
- że Koide jest własnością geometrii przestrzeni stanów, a nie lokalnego ODE,
- że powinienem minimalizować funkcjonał dla całej triady, a nie dla jednego leptona naraz,
- że potrzebny jest warunek domknięcia na zbiorze trzech rozwiązań.

**Gdzie mogłem wejść w kozi róg:**
- próbując wyciągnąć Koidego bezpośrednio z mapy `g_0 -> A_tail`,
- interpretując brak lokalnego crossingu jako dowód, że „tam nic nie ma”, zamiast sprawdzić strukturę globalną.

---

## 5. `Tau` może być stanem z innej gałęzi lub innego typu niż `e` i `μ`

To, że elektron i mion wychodzą relatywnie naturalnie, a `tau` stawia wyraźny opór, może znaczyć nie tyle „brak trzeciej generacji”, ile raczej „zła klasa rozwiązań została użyta do jej szukania”.

**Hipoteza robocza:** `tau` może odpowiadać stanowi:
- z inną liczbą węzłów,
- z inną gałęzią stabilności,
- z inną topologią rdzeń–ogon,
- albo stanowi hybrydowemu, którego nie da się opisać tym samym ansatzem co `e` i `μ`.

**Co mogłem przeoczyć:**
- że `tau` nie powinno być szukane jako „kolejny punkt” na tej samej gałęzi,
- że obecny ansatz jednosolitonowy może być zbyt sztywny dla trzeciej generacji,
- że trzecia generacja może wymagać dodatkowego indeksu dyskretnego.

**Gdzie mogłem wejść w kozi róg:**
- zakładając z góry, że wszystkie trzy leptony muszą być opisane dokładnie tym samym typem rozwiązania,
- ignorując możliwość, że sukces dla `e` i `μ` nie gwarantuje poprawnej klasy dla `τ`.

---

## Wersja syntetyczna

Najkrótsza robocza intuicja jest taka:

> To, co obecnie wychodzi jako elektron i mion, może być prawdziwym cieniem głębszej struktury.  
> Problem `tau` najpewniej nie polega na tym, że teoria nie ma trzeciej generacji, tylko na tym, że aktualny opis spłaszcza rodzinę rozwiązań zbyt mocno albo używa zbyt ubogiej mapy `stan -> masa`.

---

## Co dalej sprawdzać

Najbardziej sensowne następne testy:

1. policzyć stabilnie mapę `(B_tail, C_tail, A_tail)` dla szerokiej siatki rozwiązań,  
2. sprawdzić, czy warunki fazowe są stabilne numerycznie i niezależne od sposobu fitu,  
3. przetestować poprawkę do prawa masowego zależną od fazy,  
4. spróbować sformułować Koide jako warunek globalny na triadzie stanów,  
5. sprawdzić, czy `tau` nie siedzi na innej gałęzi / w innym ansatzu niż `e` i `μ`.
