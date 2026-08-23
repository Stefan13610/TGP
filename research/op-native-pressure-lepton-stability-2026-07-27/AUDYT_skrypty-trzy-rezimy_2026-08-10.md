# Audyt trzech skryptów rdzenia liczących trzy reżimy — skrypty przeczą własnym podsumowaniom

**Data:** 2026-08-10
**Typ:** AUDYT WYKONAWCZY (uruchomienie skryptów rdzenia + weryfikacja ich twierdzeń)
**Skrypty rdzenia:** `tooling/scripts/gravity/two_body_Veff.py`,
`tooling/scripts/stability/three_regimes_quantitative.py`
**Weryfikacja:** `REZIMY_audit_skrypty_rdzenia.py`

---

## 0. Werdykt

> Skrypty **nie ratują** studni — **pogłębiają problem**. Znalazłem:
> **(1)** skrypt raportuje **niestabilny** punkt równowagi jako `(stable)`;
> **(2)** skrypt drukuje **triumfalne podsumowanie sprzeczne z własnym outputem**;
> **(3)** skale fizyczne są **katastrofalne** — „studnia" protonu leży **17 rzędów poniżej
> długości Plancka**, a „grawitacja" zaczyna się przy **0,80 promienia Hubble'a**;
> **(4)** twierdzenie wyjaśniające „dlaczego makroskopowo widzimy tylko grawitację"
> **łamie się o 20+ rzędów wielkości**.
> **Za to** stabilne minimum d\* = 4β potwierdza się **dokładnie** — także w skrypcie rdzenia.

## 1. ⛔ Niestabilny punkt równowagi zaraportowany jako stabilny

`three_regimes_quantitative.py`, output dla protonu:
```
EQUILIBRIUM POINTS (proton):
  d = 1.044502e-77 r0 = 2.772411e-52 m  (stable)
  d = 4.000000e+00 r0 = 1.061716e+26 m  (stable)
```

Weryfikacja (β=1, C=0.1 — przypadek liczbowo czysty):

| pierwiastek | V'' | F tuż poniżej | F tuż powyżej | werdykt |
|---|---|---|---|---|
| **d₁ = 0.5168** | **−5.23** | **−2.71e−03** (przyciąganie) | **+2.69e−03** (odpychanie) | **NIESTABILNY** |
| d₂ = 3.4832 | +2.53e−03 | +8.85e−06 (odpychanie) | −8.80e−06 (przyciąganie) | STABILNY ✅ |

**Poniżej d₁ siła jest przyciągająca** ⟹ układ spada do d=0 (kolaps).
**Powyżej d₁ odpychająca** ⟹ ucieka do d₂. To jest **maksimum potencjału — BARIERA, nie studnia.**

⟹ **Punkt, który rdzeń nazywa „studnią" i przypisuje mu confinement kwarków,
jest niestabilnym maksimum.** Skrypt oznacza go `(stable)`. To błąd w skrypcie,
a nie tylko w narracji.

## 2. ⛔ Podsumowanie sprzeczne z własnym outputem

`two_body_Veff.py`, sekcja [2] — **rzeczywisty output**:
```
[2] Numerical V_eff from single-source profile superposition
  Effective C = 0.00538
  Force zero crossings at d = ['14.92']        <-- JEDNO przejście przez zero
```
Ten sam skrypt, `SUMMARY`:
```
2. Numerical superposition of single-source profiles reproduces
   the three-regime structure in E_int(d) decomposition.
```

**Jedno przejście przez zero = DWA reżimy, nie trzy.** Numeryczna część skryptu
**nie odtwarza** struktury trzech reżimów, a podsumowanie twierdzi, że odtwarza.

> To jest **dokładnie ten sam wzorzec**, który złapałem u siebie przy Skyrme v1:
> zakodowany na sztywno triumfalny print, niezależny od tego, co faktycznie wyszło.
> Tam był mój; tu jest w rdzeniu.

## 3. ⛔ Skale fizyczne — wynik katastrofalny

Z `three_regimes_quantitative.py`: **r₀ = 2,654e+25 m** (dla porównania R_Hubble = 1,322e+26 m,
czyli **r₀ ≈ 0,2 R_Hubble**).

| obiekt | d_well [m] | d_grav [m] | d_well / l_Planck | d_grav / R_Hubble |
|---|---|---|---|---|
| **proton** | **2,77e−52** | **1,06e+26** | **1,7e−17** | **0,803** |
| elektron | 1,51e−55 | 1,06e+26 | 9,3e−21 | 0,803 |
| Ziemia | 9,90e−01 | 1,06e+26 | 6,1e+34 | 0,803 |

**Konsekwencja dla dwóch protonów w odległości 1 fm:** 2,77e−52 ≪ 1e−15 ≪ 1,06e+26,
czyli **reżim II — ODPYCHANIE**. Nie grawitacja, nie studnia.

- **„Studnia" protonu leży 17 rzędów wielkości PONIŻEJ długości Plancka** — poza jakąkolwiek
  dziedziną stosowalności teorii kontinualnej (`rem` w sek08 ogranicza opis do ψ ∈ (ε_sub/Φ₀, 4/3)).
- **Reżim I (grawitacja) zaczyna się dopiero powyżej ~10²⁶ m.**

Rdzeń (`rem:trzy-rezimy-fizyczne`) twierdzi:
> „Trzy reżimy (grawitacja, odpychanie, uwięzienie) występują na skalach **subatomowych**."

**To jest ilościowo fałszywe według własnego skryptu rdzenia** — o ~37 rzędów w dół
i ~41 rzędów w górę.

## 4. ⛔ Wyjaśnienie „dlaczego makro widzimy tylko grawitację" nie działa

Skrypt: **M_crit = 1,601e+50 kg = 8,05e+19 M_☉**.

| obiekt | M [kg] | M < M_crit? | wg rdzenia |
|---|---|---|---|
| Ziemia | 5,97e+24 | ✅ | **trzy reżimy** |
| Słońce | 1,99e+30 | ✅ | **trzy reżimy** |
| Droga Mleczna | 2,98e+42 | ✅ | **trzy reżimy** |
| gromada galaktyk | 1,99e+45 | ✅ | **trzy reżimy** |
| obserwowalny Wszechświat | 1,50e+53 | ❌ | tylko grawitacja |

Rdzeń twierdzi, że warunek jest „**złamany dla obiektów makroskopowych** (C≫1), co tłumaczy,
dlaczego na skalach makroskopowych obserwujemy **wyłącznie** grawitację".

**Wg skryptu łamie się dopiero powyżej 8×10¹⁹ mas Słońca** — czyli **dla niczego
we Wszechświecie poza całym obserwowalnym Wszechświatem**. Ziemia, Słońce, galaktyki
i gromady **wszystkie** mają „trzy reżimy" i powinny się odpychać aż do ~10²⁶ m.

**Wyjaśnienie zawodzi o ponad 20 rzędów wielkości i jest sprzeczne z obserwacją.**

## 5. ⚠ Rozbieżność skrypt ↔ tekst rdzenia (na korzyść skryptu)

`two_body_Veff.py` lin. 159, 175–176, 186 liczy całki na **δΦ = Φ − Φ₀** (zaburzeniu):
```python
f1 = f_single(r) - 1.0   # δf₁ = f - 1 (perturbation from vacuum)
I2 = ...  # ∫ δf₁·δf₂ d³x  (for E_β)
```
ale `eq:Ebeta`/`eq:Egamma` w rdzeniu są zapisane na **pełnym Φ₁Φ₂**.

**Skrypt ma rację, tekst nie.** Energia oddziaływania musi być zbudowana z zaburzeń —
tło Φ₀ kasuje się. Co więcej, `eq:Eint` = E[Φ₁+Φ₂] − E[Φ₁] − E[Φ₂] jest **dosłownie
niespójne**: jeśli Φ₁,Φ₂ → Φ₀ w nieskończoności, to Φ₁+Φ₂ → **2Φ₀**, czyli zła próżnia.
Poprawna postać to Φ = Φ₀ + δΦ₁ + δΦ₂. **Skrypt cicho poprawia błąd w równaniach rdzenia.**

## 6. ✅ Co się POTWIERDZA — i to dokładnie

Mój wynik z poprzedniego audytu: stabilne minimum leży przy **d\* = 4β**, niezależnie
od wątpliwego członu kwartycznego. Skrypt rdzenia daje dla protonu:

```
d = 4.000000e+00 r0     <-- dokładnie 4β
```
oraz d₁ = 1,0445e−77 r₀ = **4,5·C** (C_proton = 2,321e−78) — obie asymptotyki potwierdzone
co do cyfry. **Odległość równowagowa pary jest wynikiem odpornym.**

**Ale ma cenę, której nie zamiatam:** przy tej kalibracji d\* = 4r₀ = **1,06e+26 m**.
Czyli „rozmiar balonu" wychodzi **kosmologiczny**, nie cząstkowy — bo r₀ i β są
kalibrowane **ciemną energią** (Φ₀ ≈ 25 z Λ_obs, γ = Φ₀H₀²/c₀²).

To jest **ten sam problem**, który znalazłem w `ANALIZA_miara-ze-skali_wieza-v2`:
> „przy kosmologicznej kalibracji defekty Φ są obiektami GALAKTYCZNYMI, nie cząstkami"

**Struktura wyniku jest poprawna; skala jest odziedziczona po kalibracji kosmologicznej.**
Dopóki hierarchia Φ₀ między domenami jest wolnym parametrem, d\* nie jest predykcją
rozmiaru cząstki.

## 7. Bilans dla trzech reżimów po tym audycie

| | Status |
|---|---|
| **Reżim I — grawitacja** | ✅ **stoi** (E_lin, Yukawa pomijalna do 50 rzędów) |
| **Reżim II — odpychanie** | ⚠ **znak pewny**, skalowanie (1/d)ln(d/r₀) bez pokrycia, **zasięg absurdalny** (do 10²⁶ m dla protonu) |
| **Reżim III — studnia/confinement** | ⛔ **UPADA**: poza zakresem ważności funkcjonału, dwa sprzeczne wzory na skalę, punkt **niestabilny** (bariera, nie studnia), położenie **17 rzędów pod Planckiem** |
| **d\* = 4β (odległość równowagowa)** | ✅ **odporne**, potwierdzone skryptem rdzenia — ale skala kosmologiczna |
| **„makro ⟹ tylko grawitacja"** | ⛔ **fałszywe** o 20+ rzędów (M_crit = 8e19 M_☉) |

## 8. Czego NIE twierdzę

- Nie twierdzę, że TGP nie ma confinementu — twierdzę, że **cała ścieżka §3 → ssec:dwa-zrodla
  → te trzy skrypty go nie ustanawia**, a w kilku miejscach twierdzi coś przeciwnego do
  własnych liczb.
- Nie sprawdziłem `two_source_potential.py` ani `effective_potential.py` — możliwe,
  że któryś używa innej regularyzacji. To zostaje otwarte.
- Nie twierdzę, że kalibracja Φ₀≈25 jest jedyna — rdzeń **sam** deklaruje Φ₀ jako
  „EFT scale-dependent free parameter". Przy innym Φ₀ skale się przesuną. **Ale wtedy
  d\* nie jest predykcją.**

---

**Uruchomione:** `tooling/scripts/gravity/two_body_Veff.py`,
`tooling/scripts/stability/three_regimes_quantitative.py`
**Weryfikacja:** `REZIMY_audit_skrypty_rdzenia.py`
**Powiązane:** `AUDYT_ssec-dwa-zrodla_2026-08-10.md`, `ANALIZA_trzy-rezimy_intuicja-gradientowa_2026-08-10.md`
