# Substrat → wstążki: target Φ-2 wymuszony na 2T (binarna tetraedralna)

**Data:** 2026-07-30
**Typ:** WYNIKI (pierwsza strzałka łańcucha) — poziom topologiczno-grupowy
**Weryfikacja:** dwa niezależne skrypty, uruchomione u źródła (nie z narracji):
`RIBBONS_group_compare.py`, `RIBBONS_uniqueness_2T.py`
**Status uczciwy:** strzałka substrat→wstążki jest wyprowadzona **tylko na poziomie
klasyfikacji topologicznej**. Warstwa dynamiczna (realizacja jako metastabilne obiekty,
background-free) pozostaje OTWARTA. To NIE zamyka strzałki — ustala jej szkielet i target.

---

## 0. Kontekst i decyzje

- Arena substratu: **Φ-2 (nieabelowa)** — decyzja użytkownika. Motyw: uciec od
  obalonego no-go (abelianizacja B₃=ℤ / writhe nie śledzi trialności).
- Ontologia: substrat = pole parametru porządku Φ: ℝ³ → M = SO(3)/H (nematyk-podobny).
  **Wstążka = linia dysklinacji** (framed line): rdzeń = defekt liniowy, framing = obrót 2π.
- Klasyfikacja: linie defektów ⟷ **π₁(M) = H̃** = preobraz H w SU(2) (grupa binarna).
  To standardowy wynik teorii defektów (Volovik–Mineev, Mermin).
- Dopuszczalne targety = **wszystkie skończone podgrupy SU(2)** (klasyfikacja ADE):
  cykliczne ℤ_n, binarne dwuścienne Dic_n, oraz 2T, 2O, 2I.

## 1. Warunek konieczny na kolor (nie postulat — wymóg §0)

Żeby kolor (trialność = ℤ₃ = A₃) mógł WYPAŚĆ ze wstążki, a nie być wsadzony:
1. π₁(M) **nieabelowe** (inaczej = obalony no-go writhe→ℤ),
2. **ℤ₃ przeżywa abelianizację** (3 | |G/G'|) — czyli trialność jest widoczna na
   poziomie obserwowalnym (abelowym), a nie ginie jak writhe.

## 2. Wynik rachunku (zweryfikowany kodem)

Przegląd WSZYSTKICH skończonych podgrup SU(2), abelianizacja |G/G'|:

| Grupa | |G| | nieabelowa? | \|G/G'\| | 3 \| ab? | dopuszczalna? |
|---|---|---|---|---|---|
| ℤ_n (bin. cykliczne) | n | NIE | ℤ_n | — | nie (abelowa) |
| Dic_1..6 (bin. dwuścienne) | 4n | tak (n≥2) | rzędu 4 | **nie** | **nie** |
| **2T (bin. tetraedralna)** | **24** | **tak** | **ℤ₃** | **TAK** | **✅ JEDYNA** |
| 2O (bin. oktaedralna) | 48 | tak | ℤ₂ | nie | nie |
| 2I (bin. ikozaedralna) | 120 | tak (doskonała) | trywialna | nie | nie |

> **2T jest JEDYNYM targetem Φ-2, w którym kolor (ℤ₃) przeżywa abelianizację.**
> Kolor nie jest wybrany — jest **wymuszony** przez warunek konieczny. To jest duch §0
> (lock/struktura ma wypaść, nie być postulowana), zrealizowany na poziomie substratu.

Uwaga uczciwości: „rzędu 4" dla Dic_n to ℤ₄ (n nieparz.) lub ℤ₂×ℤ₂ (n parz.) — nieistotne,
bo 3 ∤ 4 w obu przypadkach. Cała rodzina Dic_n (w tym Q₁₂ i biaxialny Q₈) odpada tym samym
mechanizmem, który obalił „lock z warkoczy": abelianizacja gubi ℤ₃.

## 3. Struktura wstążek dla 2T (= SL(2,3))

Z `RIBBONS_group_compare.py` (zweryfikowane):

- **7 typów wstążek** = 7 klas sprzężoności. Rozmiary/rzędy:
  (1,ord1)=próżnia, (1,ord2)=centralny −1, (6,ord4), dwie (4,ord3), dwie (4,ord6).
- **Rozkład z JEDNEJ grupy** (to jest sedno):
  - **centrum ℤ₂ = {1, −1}** → element −1 (rząd 2) = generator obrotu 2π =
    **spin-½** (Finkelstein–Rubinstein / Berry π). Zgodne z filarem §2, spełnia §5.6.
  - **abelianizacja ℤ₃** (trzy 1-wym. reprezentacje) = **kandydat na 3 kolory** (trialność).
  - **komutator G' = Q₈** (rząd 8) = nieabelowy rdzeń (dawny biaxialny Q₈ — tu jako
    podstruktura, nie jako target).
  - **dwie klasy rzędu 3** = naturalny kandydat na **kolor / antykolor**.

## 4. Bilans wg reguł RAMY §5 (uczciwie)

| Reguła | Status na tym etapie |
|---|---|
| §5.1 metastabilny z topologii (nie hack) | ✅ klasa nietryw. π₁ = bariera topologiczna (textbook), NIE basen/hack |
| §5.2 składalny | ✅ iloczyn grupowy; nieabelowość → braiding fizyczny |
| §5.3 minimalny | ✅ i mocniej: **uniqueness** wśród π₁(SO(3)/H) |
| §5.4 samoreferencyjny / background-free | ⚠ OTWARTE — OP-space nie zakłada pola EM, ale formalizacja background-free niezrobiona |
| §5.5 kolor+ładunek z JEDNEJ topologii | ⚠ CZĘŚCIOWO — kolor-kandydat i spin z jednej grupy ✅; **ładunek jeszcze niezdefiniowany** na wstążce; lock NIEROZSTRZYGNIĘTY |
| §5.6 spin z nawinięcia | ✅ centralny −1 = generator FR/Berry |

## 5. Co jest DERIVED / co OTWARTE (klasyfikacja audytowa)

**DERIVED (teoria grup, zweryfikowane kodem):**
- Wstążki (linie dysklinacji) w Φ-2 klasyfikowane przez π₁ = grupa binarna.
- Target **wymuszony na 2T** przez warunek konieczny na trialność (uniqueness).
- 2T niesie spin (centrum ℤ₂), kolor-kandydat (abelianizacja ℤ₃), rdzeń Q₈ — z jednej grupy.
- Metastabilność i składalność wstążek z topologii (nie hack).

**OTWARTE (nie rości sobie — następne warstwy):**
- **Ładunek** na wstążce niezdefiniowany; lock kolor↔ładunek (§0) NIEROZSTRZYGNIĘTY. To jest
  właściwa treść NASTĘPNYCH strzałek, nie tej.
- **Dynamiczna realizacja:** czy istnieją energetycznie metastabilne obiekty (nie tylko klasy
  topologiczne)? Czy substrat Φ→SO(3)/T jest realizowalny background-free (§5.4)?
- **Mapa 3 kolory / 3 antykolory:** dwie klasy rzędu 3 + ℤ₃-abelianizacja trzeba dokładnie
  zmapować na strukturę koloru (to już zahacza o strzałkę wstążki→kwarki — NIE ruszać teraz).
- **NIE mylić** trialności (3 kolory, tu wymuszone) z liczbą generacji (pozostaje OTWARTA, §5.10).

## 6. Jednozdaniowo

> W arenie Φ-2 wstążki = linie dysklinacji klasyfikowane przez π₁(SO(3)/H); warunek, by kolor
> (ℤ₃) przeżył abelianizację, **wymusza jednoznacznie target 2T**, który z jednej grupy niesie
> spin (centrum ℤ₂) i kandydat-kolor (abelianizacja ℤ₃). To wyprowadza **szkielet topologiczny**
> strzałki substrat→wstążki; ładunek, lock i dynamika pozostają OTWARTE.

---

**Pliki tej sesji:**
- `RIBBONS_group_compare.py` — struktura Q₁₂ vs 2T (klasy, centrum, abelianizacja)
- `RIBBONS_uniqueness_2T.py` — uniqueness 2T wśród skończonych podgrup SU(2)

**Powiązane:**
- `RAMY_ladunek-kolor-nie-niezalezne_2026-07-28.md` (reguły §5, oś §0)
- `AUDYT_KRYTYCZNY_2026-07-28.md` (no-go writhe→ℤ, wzorzec obalania)
- `research/qm_spin/README.md` (filar spin-½: centralny −1 / FR)
