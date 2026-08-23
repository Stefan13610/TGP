# Ciśnienie sąsiadów a lokalizacja — test tezy „pojedynczy soliton nie może istnieć"

**Data:** 2026-07-31
**Typ:** WYNIKI (warstwa dynamiczna/lokalizacja) — **wynik NEGATYWNY, mocny**
**Weryfikacja:** `PRESSURE_wigner_seitz.py` (uruchomione; analityka + numeryka; oba sektory)
**Falsyfikator zadeklarowany PRZED liczeniem:** brak minimum E(L) ⟹ teza pada w danym sektorze.

---

## 0. Teza i jej realizacja

Teza (użytkownik, zgodna z RAMY §3 „tło = zbiorowa przestrzeń innych obiektów"):
**pojedynczy soliton nie może istnieć; istnieje tylko dzięki ciśnieniu innych.**
Ściana = dowód, że konfiguracja wielu obiektów jest korzystniejsza energetycznie niż rozpad.

Realizacja: **komórka Wignera-Seitza** — obiekt w kuli promienia L, na brzegu **u'(L)=0**
(symetria z sąsiadami) zamiast **u(∞)=1** (próżnia). To celowo usuwa premisę testu izolowanego:
runaway z AUDYTU (DODATEK A) liczony był właśnie z u(∞)→1, więc dotyczył **tylko obiektu w próżni**.

## 1. Wynik decydujący (analityczny, nie numeryczny)

Energia statyczna w sektorze kanonicznym:
```
E[u] = 4π ∫ [ ½u⁴u'² + (W(u) − W(1)) ] r² dr
```
- człon kinetyczny **½u⁴u'² ≥ 0** (dodatnio określony),
- **W(u) − W(1) ≥ 0**, bo W(u)=u⁸/8−u⁷/7 ma **JEDNO globalne minimum** przy u=1
  (W'(u)=u⁶(u−1) → jedyne zero na u>0; zweryfikowane numerycznie).

⟹ **E[u] ≥ 0, z równością TYLKO dla u ≡ 1.**

> **Żadna niejednorodna konfiguracja nie może być korzystniejsza energetycznie od jednorodnej —
> przy ŻADNEJ gęstości.** Ciśnienie sąsiadów nie ma czego stabilizować, bo nie ma nic, co
> powstrzymywałoby relaksację do próżni jednorodnej.

## 2. Potwierdzenie numeryczne (komórka W-S)

Strzelanie u(0)=g₀, u'(0)=0 → u'(L) dla L∈{1,2,5,10}, g₀∈{0.5…2.0}:
**żadne (g₀,L) nie spełnia warunku Neumanna u'(L)=0** (poza trywialnym u≡1); dla większych L
wszystko ucieka (runaway) do 0 lub ∞.

Powód mechaniczny: równanie = cząstka w potencjale **−W** z tarciem 2u'/r. Skoro u=1 jest
**minimum W**, to jest **maksimum −W** ⟹ równowaga niestabilna: cząstka puszczona w spoczynku
przy u≠1 stacza się **monotonicznie** i nigdy nie wraca do u'=0.

**E(L) = 0 identycznie dla każdego L** (jedyne rozwiązanie to u≡1) ⟹ **brak minimum przy
skończonym L** ⟹ **zadeklarowany falsyfikator zadziałał: teza PADA w tym sektorze.**

## 3. Kontrola drugiego sektora (żeby nie orzekać z jednego)

Sektor grawitacyjny M9.1'' z elementem objętościowym √−g = c₀ψ/(4−3ψ):
V_eff = √−g · V_M911 = **−c₀γψ³(4−3ψ)/12**, V_eff′ = −c₀γψ²(1−ψ) → krytyczne ψ=1;
**również JEDNO globalne minimum przy ψ=1** (zweryfikowane numerycznie).
⟹ twierdzenie z §1 stosuje się do **OBU** kanonicznych sektorów.

## 4. Retroaktywna diagnoza (dlaczego stare „pressure" MUSIAŁO być fitem)

AUDYT §4.3 zaklasyfikował `pressure+loops "111%"` jako **overfitting** (jawny skan parametru
do celu). Teraz wiemy **strukturalnie DLACZEGO to nie mogło zadziałać**: w kanonicznym sektorze
skalarnym jednorodna próżnia jest **globalnym** minimum energii, więc żadne „ciśnienie" nie
ustabilizuje grudki. Skan parametru zaklejał **niemożliwość strukturalną**. To zamyka tamten
wątek definitywnie — nie jako „źle policzone", lecz „niewykonalne w tym sektorze".

## 5. Czego to NIE obala (ważne)

**NIE obala ontologii użytkownika.** Obala jej realizację **w sektorze skalarnym**. Wynik wprost
wskazuje, czego potrzeba: **więzu zabraniającego u≡1**. Kandydaci:

1. **Sektor topologiczny (wstążki 2T)** — i tu jest zbieżność: udowodniliśmy już
   (`WYNIKI_domykanie-solitony`), że **pojedyncza nietrywialna wstążka NIE może się domknąć**
   (wszystkie 6 klas: `False`) i że domknięcie wymaga **partnerów** o sumie trialności 0.
   Topologia **zabrania** relaksacji do próżni — dokładnie brakujący więz. I niezależnie
   odtwarza tezę „pojedynczy obiekt nie istnieje sam".
2. **Źródło materii ρ** (L_mat = −(q/Φ₀)ψρ) — ale wtedy cząstka jest **wstawiona ręcznie**,
   więc nie nadaje się do wyprowadzania cząstek (ryzyko cyrkularności).

## 6. Klasyfikacja audytowa

- **DERIVED (mocne, analityczne):** E[u]≥0 z równością iff u≡1 w obu kanonicznych sektorach;
  brak niejednorodnego rozwiązania w komórce W-S; E(L)=const ⟹ brak minimum.
- **DERIVED (retroaktywnie):** strukturalny powód, dla którego stare „pressure+loops" musiało być fitem.
- **OTWARTE:** energetyka w sektorze topologicznym (brak natywnego funkcjonału — Skyrme
  zdyskwalifikowany); czy więz topologiczny + ciśnienie dają minimum przy skończonym L.
- **NIE twierdzę:** że teza jest obalona globalnie; że sektor topologiczny ją uratuje (nie policzone).

## 7. Jednozdaniowo

> W obu kanonicznych sektorach skalarnych jednorodna próżnia jest **globalnym** minimum energii
> (E≥0, równość iff u≡1), więc ciśnienie sąsiadów nie ratuje lokalizacji i zadeklarowany
> falsyfikator zadziałał — ale to wskazuje precyzyjnie brakujący element: **więz topologiczny
> zabraniający próżni**, który w sektorze wstążek 2T **już mamy udowodniony** (pojedyncza wstążka
> nie domyka się, wymaga partnerów).

---

**Plik:** `PRESSURE_wigner_seitz.py`
**Powiązane:** `WYNIKI_domykanie-solitony_2026-07-31.md` (wstążka nie domyka się sama),
`AUDYT_KRYTYCZNY_2026-07-28.md` (DODATEK A: runaway z u(∞)→1; §4.3 pressure=overfitting),
`STATE_lancuch-substrat-kwarki_2026-07-31.md`
