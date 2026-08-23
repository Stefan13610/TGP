# Czy pierwotna intuicja gradientowa przetrwała? TAK — §3 rdzenia. Plus audyt: co w niej stoi, a co nie

**Data:** 2026-08-10
**Typ:** ARCHEOLOGIA RDZENIA + AUDYT (weryfikacja algebry i zakresu ważności)
**Pytanie użytkownika:** czy wersja „przestrzeń = pole generowane przez zaburzenie z rozkładem
gradientowym; przyciąganie / odpychanie na bliskich / lock" przetrwała rozwój teorii?
**Skrypt:** `REZIMY_audit_trzy-rezimy.py`

---

## 0. Odpowiedź krótka

> **Przetrwała — i to jako jedna z głównych sekcji rdzenia (`sek03_rezimy`), nazwana wprost
> „centralną predykcją TGP".** Wszystkie trzy Twoje podpisy są tam, z tymi samymi nazwami:
> przyciąganie (grawitacja) / odpychanie / studnia (lock). Algebra warunku istnienia
> trzech reżimów jest **poprawna** (sprawdziłem).
> **ALE:** „studnia" — czyli właśnie ten lock — jest w rdzeniu **bezdenna** (kolaps, nie confinement),
> **leży poza zakresem ważności wzoru**, którym jest wyprowadzana, a **etykiety d_rep/d_well
> są zamienione** względem tabeli reżimów.

## 1. Co dokładnie przetrwało — cytaty

`sek03_rezimy` lin. 6–10:
> „**Centralną predykcją TGP** jest istnienie **trzech odrębnych reżimów** oddziaływania,
> wynikających z **interferencji pól Φ generowanych przez różne źródła**. Nie są to trzy osobne
> siły — jest to **jedna dynamika** pola Φ, która na różnych skalach daje różny znak siły."

**Mechanizm = dokładnie Twój** (lin. 15–24): **F ∝ −∇V_eff(Φ)**, i „V_eff **nie jest monotoniczną**
funkcją odległości od źródła".

**Grawitacja** (lin. 64–72):
> „każdy obiekt generuje »wygórowanie« Φ wokół siebie. Dwa odległe obiekty mają między sobą
> obszar, gdzie ich pola nakładają się **konstruktywnie**… To jest grawitacja w TGP:
> **nie krzywizna gotowej geometrii, lecz siła wynikająca z konstruktywnej interferencji
> wygenerowanej przestrzeni.**"

**Odpychanie** (lin. 94–113) — z „stromego wzrostu Φ" blisko źródła; przypisane mu: stabilność
materii, **zasada nieoznaczoności** (prop:uncertainty-tgp), ciśnienie degeneracji.

**Lock** (lin. 259–281, `prop:well`) — „dwa źródła są »sklejone« we wspólnym dołku pola Φ";
przypisane: **confinement kwarków**, stany hadronowe, asymptotyczna swoboda.

⟹ **Twoja pamięć jest trafna. Nic z tej intuicji nie zostało po drodze wyrzucone.**
Co więcej — to §3 jest źródłem „reżimu III / studni przestrzennej", na którą powołuje się
`sek09` przy konfinowaniu gluonów.

## 2. Co jest poprawne (zweryfikowane)

`prop:trzy-rezimy-beta-gamma`, przy β=γ:
```
V_eff(d) = −4πC²/d + 8πβC²/d² − 24πβC³/d³
F = 0  ⟺  d² − 4βd + 18βC = 0  ⟺  β > 9C/2
```
**Algebra przeliczona ręcznie i numerycznie — zgadza się.** Pierwiastki zerują F do ~1e−16.
Struktura znaków potwierdzona: **przyciąganie → odpychanie → przyciąganie** (3 reżimy).

Ładna i nietrywialna konsekwencja (`rem:trzy-rezimy-fizyczne`): warunek β > 9C/2 jest
spełniony dla **cząstek** (C≪1), a **złamany dla obiektów makroskopowych** (C≫1) — wtedy
zostaje **sam reżim I (grawitacja)**. To jest zgodne z obserwacją i **nie jest dopasowane**.

## 3. Trzy problemy, których rdzeń nie odnotowuje

### (a) ⛔ „Studnia" jest BEZDENNA — to kolaps, nie confinement

| d | V(d) | F(d) |
|---|---|---|
| 0.200 | −3.77e+00 | −8.17e+01 |
| 0.100 | −5.15e+01 | −1.77e+03 |
| 0.010 | −7.29e+04 | −2.21e+07 |
| 0.001 | −7.51e+07 | −2.26e+11 |

V ~ −24πβC³/d³ → **−∞**, F → **−∞**. Analiza ekstremów:

| pierwiastek | d | V'' | typ |
|---|---|---|---|
| d₁ (mniejszy) | 0.5168 | **−5.23** | **MAKSIMUM (niestabilne)** |
| d₂ (większy) | 3.4832 | **+2.52e−03** | **MINIMUM (stabilne)** |

**Jedyne stabilne minimum leży przy WIĘKSZYM pierwiastku** — na granicy odpychanie/grawitacja.
„Studnia" nie jest dołkiem, tylko **zboczem do −∞**.

Rdzeń (`ssec:studnia`) twierdzi: „asymptotyczna swoboda: wewnątrz studni gradienty są małe
(**płaskie dno**)". **W eq:Veff-beta-eq-gamma płaskiego dna NIE MA.** To twierdzenie
nie ma pokrycia we wzorze, do którego się odwołuje.

### (b) ⛔ Reżim III leży POZA zakresem ważności swojego wzoru

Szereg −C²/d + βC²/d² − βC³/d³ to rozwinięcie w potęgach **1/d** (duże d). Stosunki członów:

| d | E_β/E_lin | E_γ/E_β |
|---|---|---|
| 10.0 | 0.20 | 0.03 |
| 2.0 | **1.00** | 0.15 |
| 1.0 | **2.00** | 0.30 |
| 0.2 | **10.00** | **1.50** |

Przy **d < 2β kolejne człony przewyższają poprzednie** — szereg przestaje być asymptotyczny.
A d₁ = 0.517 (granica studni) leży **głęboko w tym obszarze** (E_β/E_lin ≈ 3.9).

⟹ **Studnia jest wyprowadzona przez ekstrapolację rozbieżnego szeregu poza jego zakres.**
To nie znaczy, że studni nie ma — znaczy, że **ten rachunek jej nie pokazuje**.

### (c) ⚠ Etykiety d_rep / d_well są zamienione

`eq:d-rep` = 2β − √(...) = **mniejszy** pierwiastek; `eq:d-well` = 2β + √(...) = **większy**.
Ale tabela w `ssec:diagram` mówi **r_well < r_rep**. **Sprzeczność nazw** — kosmetyczna,
ale każdy rachunek numeryczny na tych wzorach wyjdzie odwrotnie.

### (d) ⚠ §3 jest napisana w Formulacji B, a kanoniczna jest A

Przypis w lin. 97–101 mówi wprost:
> „W sformułowaniu kanonicznym (**Formulacja A**): g = Φ/Φ₀, **g < 1 blisko materii
> (pole MALEJE, nie rośnie)**… Opis w konwencji Φ **rosnącego** blisko materii odpowiada
> **Formulacji B**."

Twoja pierwotna intuicja („im bliżej źródła, tym większe natężenie pola") to **Formulacja B**.
Cała narracja §3 jest w B. **Kanoniczna jest A.** To nie unieważnia struktury trzech reżimów
(znaki siły są konwencjo-niezależne), ale znaczy, że **sekcja niosąca Twoją intuicję jest
napisana w nie-kanonicznej konwencji** — i przy przepisywaniu na A trzeba to zrobić uważnie.

## 4. Co z tego jest DOBRĄ wiadomością — i to konkretną

**Zewnętrzna struktura potencjału jest wiarygodna, wewnętrzna nie.** Sprawdziłem, gdzie
szereg jeszcze działa:

- przy **d₂ = 3.483**: E_β/E_lin = 0.57, E_γ/E_β = 0.086 ⟹ **szereg zbieżny, d₂ jest w zakresie**
- przy **d₁ = 0.517**: E_β/E_lin = 3.87 ⟹ **poza zakresem**

⟹ **Stabilne minimum d₂ jest wynikiem, któremu można ufać.**

I to jest bezpośrednio to, czego szuka obraz balonów: **d₂ to naturalna odległość równowagowa
między dwoma obiektami** — czyli **rozmiar balonu / stała upakowania**, wyprowadzona z (β, C),
a nie dopasowana. To jest szczebel wieży, który **stoi**.

## 5. Odpowiedź na „czy to najcięższe obliczeniowo"

**Nie — okazuje się najlżejsze z dostępnych dróg**, bo rdzeń ma już zamknięty rachunek
dwóch źródeł (`ssec:dwa-zrodla`, eq:Eint-decomp–eq:Egamma). Ciężkie jest dopiero to,
co **poza** nim: konfiguracja **wielu** źródeł z nakładaniem (Twoje przenikające się balony),
gdzie superpozycja przestaje działać, bo dynamika jest nieliniowa.

## 6. Czego NIE twierdzę

- Nie twierdzę, że confinement w TGP nie działa — twierdzę, że **`eq:Veff-beta-eq-gamma`
  go nie pokazuje**, bo rozbiega się do −∞ i leży poza zakresem ważności.
- Nie sprawdzałem `ssec:dwa-zrodla` (pełna dekompozycja E_int) — możliwe, że tam
  regularyzacja małych d jest zrobiona. **To jest następny check, nie wniosek.**
- Nie twierdzę, że d₂ to „rozmiar cząstki" — to **odległość równowagowa pary**. Przejście
  od pary do upakowania wielu obiektów jest **niezrobione**.

---

**Źródła (u źródła):** `core/sek03_rezimy/sek03_rezimy.tex` lin. 6–10, 15–51, 56–82 (reżim I),
87–123 (reżim II + nieoznaczoność), 259–281 (reżim III), 289–304 (tabela), 415–489
(prop:trzy-rezimy-beta-gamma + dowód + interpretacja), przypis lin. 97–101 (Formulacja A vs B);
`core/sek09_cechowanie` lin. 591 (konfinowanie z reżimu III)
**Rachunek:** `REZIMY_audit_trzy-rezimy.py`
