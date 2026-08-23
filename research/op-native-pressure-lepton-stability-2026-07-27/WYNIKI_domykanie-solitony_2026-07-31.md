# Domykanie wstążek w solitony — co wpada automatycznie, a co nie

**Data:** 2026-07-31
**Typ:** WYNIKI (strzałka kwarki→solitony, warstwa kombinatoryczna)
**Weryfikacja:** `RIBBONS_closure_solitons.py` (uruchomione, czysta teoria grup)
**Reframe (decyzja użytkownika):** ŁADUNEK zdjęty z tej warstwy — jest **relacyjny**, pojawia się
dopiero w relacji wielu obiektów (zgodne z RAMY §3: „ładunek = selektywny przekaz pędu,
ujawniany przez pomiar zmieniający pęd" i §4: e–p = wiele ciał). To usuwa błąd kategorii
z poprzedniego kroku (audyt: lock via holonomia jednego obiektu = tautologia).

---

## 0. Postawienie pytania (bez ładunku)

Kolekcja wstążek o klasach C₁…Cₙ **domyka się** (anihiluje do próżni) ⟺ istnieją reprezentanci
gᵢ∈Cᵢ z **g₁g₂…gₙ = e** w π₁=2T. (Dla defektów nieabelowych klasa jest określona do sprzężenia,
więc pytamy o ISTNIENIE reprezentantów.) Jeśli iloczyn = **−1**, konfiguracja domyka kolor,
ale zostaje centralny element spinowy — to nie jest próżniowe domknięcie.

## 1. Co WPADA AUTOMATYCZNIE (DERIVED, zero fitu, zero parametrów)

**(a) Konfinement topologiczny pojedynczej wstążki.**
Żadna nietrywialna wstążka nie domyka się sama — dla WSZYSTKICH 6 nietrywialnych klas
(C1…C6) wynik `domyka się? False`. **Każda nietrywialna wstążka jest topologicznie związana:**
nie może zniknąć w izolacji, musi kończyć się na innych. To automatyczna, kombinatoryczna
reguła — nie hack, nie energetyka.

**(b) Trialność-0 jest KONIECZNA.**
Wszystkie domykające się pary mają sumę trialności 0 (5/5). Domknięcie ⟹ singlet kolorowy.
Kandydat „mezon" = **C4+C5** (kolor+antykolor, rząd 3) domyka się do e.

**(c) Struktura par i trójek.**
- Pary domykające: C0+C0, C1+C1, C2+C6, C3+C3, **C4+C5** — wszystkie trialność-0.
- Trójki trialności (1,1,1) (kandydat „barion"): wszystkie 4 kombinacje C2/C4 **mogą** się domknąć.

**(d) Reszta bezbarwna = sektor Q₈.**
Suma trialności 0 ⟺ iloczyn całkowity leży w Q₈ (rdzeń bezbarwny), który **zawiera −1** (spin).
Domknięcie koloru NIE usuwa etykiety spinowej.

## 2. Co NIE WPADA automatycznie (uczciwie — i to jest ważniejsze)

**(a) Trialność-0 NIE WYSTARCZA.** Z 30 trójek o sumie trialności 0 domyka się do e tylko **21**.
Istnieje **dodatkowa przeszkoda w Q₈** poza kolorem. To jest **MOCNIEJSZE niż reguła singletu
kolorowego SM** — i nie wiadomo jeszcze, czy to zaleta (dodatkowa predykcja), czy problem
(nadmiar restrykcji). **Charakteryzacja tej przeszkody = OTWARTE.**

**(b) „Barion = fermion" NIE jest wymuszone.** Dla wszystkich trójek (1,1,1) osiągalne jest
**ZARÓWNO e, JAK I −1**. Czyli statystyka złożenia **nie wypada** z samej kombinatoryki —
wymaga dodatkowej informacji (framing/ścieżki bazowe). **Nie ogłaszam „barion=fermion".**

**(c) SOLITON — czyli to, o co pytano — NIE wpada automatycznie.**
Rachunek daje **regułę selekcji** (które konfiguracje MOGĄ się domknąć), **nie dynamikę**
(że MUSZĄ się związać, że powstaje zlokalizowany obiekt o skończonym rozmiarze). Konfinement
dynamiczny, energia wiązania i rozmiar to **ta sama otwarta warstwa lokalizacji**, która
zablokowała krok 2 (kanoniczny F-A jest runaway; płaski Skyrme zdyskwalifikowany jako nie-natywny).

> **Sedno:** „domknięcie" mamy jako **topologiczną regułę selekcji**; „soliton" wymaga
> **dynamiki**, której nie mamy. Kombinatoryka mówi CO może się domknąć — nie mówi, że powstaje
> cząstka.

## 3. Klasyfikacja audytowa

- **DERIVED:** konfinement topologiczny pojedynczej wstążki; trialność-0 konieczna; tabela par
  i trójek; reszta w Q₈; trialność-0 niewystarczająca (21/30).
- **OTWARTE:** natura dodatkowej przeszkody Q₈; statystyka złożenia (e vs −1) — niewymuszona;
  **dynamika/lokalizacja (soliton jako obiekt)**; ładunek (świadomie odłożony do warstwy relacyjnej).
- **NIE twierdzę:** że barion jest fermionem; że to jest konfinement QCD; że powstaje soliton.

## 4. Jednozdaniowo

> Domykanie wpada automatycznie **jako reguła kombinatoryczna** (żadna nietrywialna wstążka nie
> zniknie sama; domknięcie wymaga trialności 0 — a nawet więcej, bo trialność-0 nie wystarcza),
> ale **soliton nie wpada** — statystyka złożenia nie jest wymuszona, a zlokalizowany obiekt
> wymaga dynamiki, której wciąż nie mamy.

---

**Plik:** `RIBBONS_closure_solitons.py`
**Powiązane:** `STATE_lancuch-substrat-kwarki_2026-07-31.md`,
`WYNIKI_s54_domkniecie_frame-native_2026-07-30.md` (2T=Q₈⋊ℤ₃),
`WYNIKI_charge-lock_Z3_2026-07-30.md` (korekta po audycie — ładunek zdjęty z tej warstwy)
