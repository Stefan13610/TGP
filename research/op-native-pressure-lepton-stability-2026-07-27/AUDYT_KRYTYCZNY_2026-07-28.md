# Audyt Krytyczny — Tier 2 / Bounce-Hierarchy

**Data:** 2026-07-28
**Typ:** Niezależny audyt techniczny (weryfikacja u źródła, nie z narracji)
**Metoda:** Ponowne przeczytanie kodu solvera + własny test diagnostyczny (`AUDIT_diagnostic.py`)
**Wniosek jednozdaniowy:** Narracja "TIER 2 SOLVED / 95% publication-ready" **nie jest poparta** tym, co faktycznie policzono. Istnieje jednak realny, wąski sygnał fizyczny, który został przysłonięty przez pomiar niewłaściwej wielkości.

---

## 0. Jak czytać ten dokument

To nie jest dokument "obalający TGP". To audyt **jednego wątku** (Tier 2 — hierarchia generacji leptonów przez bounce-hierarchy) prowadzonego w sesjach #64–#65. Oddzielam trzy rzeczy, które wcześniejsze dokumenty mieszały:

1. **Narzędzie** (solver CP-7) — solidne, przetestowane, zbieżne. ✅
2. **Realny sygnał** (mody poniżej krawędzi kontinuum) — istnieje i jest odtwarzalny. ✅ (ale ograniczony)
3. **Interpretacja i wnioski** (bounce → N_neg, "native mechanism", gotowość do publikacji) — **niepoparte / błędne**. ❌

---

## 1. USTALENIE DECYDUJĄCE: „N_neg" mierzyło artefakt pudła

### 1.1 Co pokazał test box-size

Uruchomiłem ten sam solver dla R = 40, 60, 80 (proporcjonalnie N). Wynik:

```
g0        R     N(<0)   N(<-1 edge)   lam_min
e   1.249  40     12        0         -0.995
e   1.249  60     19        0         -0.998
e   1.249  80     25        0         -0.999
mu  2.021  40     13        2         -1.282
mu  2.021  60     19        2         -1.282
mu  2.021  80     25        3         -1.282
tau 3.189  40     12        3         -4.216
tau 3.189  60     19        3         -4.216
tau 3.189  80     25        3         -4.216
```

**Vacuum F-S (puste "pole", ZERO solitonu):**
```
VACUUM R=40: N(<0)=12   N(<-1 edge)=0
VACUUM R=60: N(<0)=19   N(<-1 edge)=0
VACUUM R=80: N(<0)=25   N(<-1 edge)=0
```

### 1.2 Interpretacja (twarda)

- **N(<0) = 12 / 19 / 25 skaluje się liniowo z rozmiarem pudła R** — identycznie dla e, μ, τ **oraz dla próżni**. To jest podręcznikowa sygnatura **dyskretyzowanego kontinuum** (gęstość stanów poniżej progu rośnie z objętością pudła). To **nie są mody fizyczne**.
- Formuła F-S ma **krawędź kontinuum przy λ ≈ −1**, nie przy 0 (kod to jawnie mówi: `[C5] B4 vacuum F-S ... expect edge = -1`, werdykt `CONFIRMED-NEGATIVE (as predicted analytically)`). Wszystkie ~19 "modów ujemnych" to po prostu stany kontinuum leżące między −1 a 0.
- Zatem **liczba "N_neg" ≈ 19 nie ma nic wspólnego z solitonem ani z generacją.** Próżnia ma tyle samo.

### 1.3 Konsekwencja dla tej sesji

Skrypty uruchomione dzisiaj (`BOUNCE_NNEG_INTERPOLATION.py`, `BOUNCE_DETAILED_INTERP_QUICK.py`) liczyły **N_neg = mody poniżej 0** — czyli dokładnie ten artefakt. Stąd "determinizm":

```
bounces=0 → N_neg=19
bounces=1 → N_neg=19
bounces=2 → N_neg=19
bounces=3 → N_neg=19
```

To nie jest "N_neg jest deterministyczną funkcją bounces". To jest "N_neg ≈ const ≈ 19, bo to rozmiar pudła". **Determinizm został ogłoszony na stałej liczbie będącej artefaktem.** Wniosek "✅ DETERMINISM CONFIRMED" jest nieuprawniony.

---

## 2. CO JEST REALNE: sygnał 0/2/3 poniżej krawędzi

Kiedy policzyłem **właściwą** wielkość — mody poniżej krawędzi kontinuum F-S (−1), czyli prawdziwe stany związane / niestabilne względem kontinuum:

```
e   (g0=1.249, bounces=0): 0 modów poniżej krawędzi -> []
mu  (g0=2.021, bounces=1): 2 mody poniżej krawędzi -> [-1.282, -1.057]
tau (g0=3.189, bounces=3): 3 mody poniżej krawędzi -> [-4.216, -1.114, -1.010]
```

**To jest odtwarzalne i (w większości) niezależne od R** (test 1: kolumna `N(<-1 edge)` = 0/0/0, 2/2/3, 3/3/3 — jedno drgnięcie μ przy R=80 to mod tuż przy krawędzi).

### Ważna uczciwość:
- Wzorzec **0/2/3, który podałeś na samym początku** ("e: 0, μ: 2, τ: 3 lokalne mody ujemne") — **jest poprawny**, ale jako **licznik modów poniżej krawędzi −1**, a nie jako "N_neg".
- **Twoja intuicja fizyczna była lepsza niż skrypty, które potem napisałem.** Skrypty mierzyły złą wielkość (N_neg<0 = artefakt) i na niej ogłosiły sukces. To jest błąd metodologiczny po mojej stronie z wcześniejszych sesji, który dziś wychodzi w audycie.

---

## 3. CO POZOSTAJE OTWARTE (a było ogłoszone jako zamknięte)

Realny sygnał 0/2/3 istnieje. Ale **kluczowe pytanie fizyczne jest NIEROZSTRZYGNIĘTE**:

> Czy mody poniżej krawędzi F-S to **fizyczne niestabilności**, czy **artefakt formulacji F-S**?

Dowody, że to może być **artefakt formulacji**:
- Formulacja **F-A (dual grawitacyjny)** daje **czyste widmo** (wszystkie λ ≥ 0, brak modów poniżej krawędzi). Kod: `[C2] ... C2 verdict: PASS (N_neg=0)`.
- Jeśli F-A i F-S opisują **tę samą fizykę** (teza L04 duality), to niestabilność widoczna tylko w F-S jest **zależna od formulacji** — czyli podejrzana jako artefakt (analogia: mody-duchy w złym wyborze zmiennych/cechowania).
- Formulacja F-S ma odwrócony znak w strukturze kinetycznej (stąd krawędź przy −1, stąd nazwa w kodzie "dualism control"). To środowisko sprzyjające artefaktom.

Dopóki to nie jest rozstrzygnięte, **nie można twierdzić**, że "saddle points są natywnymi cechami kodującymi generacje". To jest **hipoteza**, nie wynik.

---

## 4. ŁAŃCUCH POCHODNY — problemy odziedziczone

### 4.1 Ekstrakcja ładunków (`Phase2_charge_extraction_v2.py`) — SŁABA METODA

- Ładunek liczony przez dopasowanie ogona `u(r)=(g−1)·r` do `B·cos(r)+C·sin(r)`, czyli **założona oscylacja o liczbie falowej dokładnie k=1** — **bez uzasadnienia**. Statyczny soliton nie ma powodu mieć ogona ~sin(r); pole masywne → ogon eksponencjalny, Goldstone bezmasowy → ogon potęgowy 1/r. Oscylacja k=1 to prawdopodobnie **ripple numeryczny / odbicie od brzegu**.
- Okno dopasowania **r ∈ [45, 58]**, tuż przy ścianie pudła **R=60** (Dirichlet + przycięcie). Dopasowywanie oscylacji przy sztucznym brzegu = dopasowywanie artefaktu brzegowego.
- Wyniki q_e=1.200, q_μ=1.107, q_τ=1.049 to liczby rzędu 1 zależne od okna — **nie ma dowodu, że to ładunki topologiczne** (te byłyby kwantowane / miałyby jasny wzór).

### 4.2 Self-consistency / pressure (`Phase3...`, `Phase4...`) — NA PIASKU

- Model ciśnienia `E = ½ Σ qᵢqⱼ G(rᵢⱼ)` z `G = −ln(r)/2π` (propagator **2D** Goldstone) **narzucony**, nie wyprowadzony, dla problemu radialnego 3D.
- Wejściem są ładunki z 4.1 (wątpliwe). Odległości równowagi ~10⁷ (po przeskalowaniu 1e7) — liczby wyglądające arbitralnie.

### 4.3 „Pressure + loops = 111%" — OVERFITTING

- `Phase4_extended_v3_scaled.py` **jawnie skanuje** scale ∈ {1, 1.5, 2, 3, 5, 10} szukając wartości trafiającej w cel.
- `N4d_N4c_HYBRID` skanuje λ_loop aż do trafienia.
- To jest **dopasowywanie parametru do z góry znanego wyniku** — nie predykcja. A do tego **naprawiało mody, które wg §1 są artefaktem pudła** (nie wymagają naprawy).

### 4.4 Błąd numeryczny do odnotowania
- F-Sp (substrate α=1) dla τ: **grid pts = 199** (z 4000) — całkowanie profilu się załamało (g spadło do 0.16, pudło przycięte). Wynik `N_neg=1, lam_min=−2.125` dla `sub_tau` jest **numerycznie nieważny** i nie powinien być cytowany.

---

## 5. CO JEST SOLIDNE (bilans nie jest zerowy)

- **Solver CP-7 (`Phase2_bvp_spectrum.py`)** — dobrze zbudowany: self-adjoint flux form, staggered grid, symetryzacja B^(−1/2)MB^(−1/2), harness testowy z werdyktami PASS/FAIL, zbieżność potwierdzona (N=2000/4000/8000 dają stabilne lam_min). To jest **wartościowe, wielokrotnego użytku narzędzie**.
- **Ściana-duch przy G_GHOST = e^(−1/4)** — realna cecha F-S (F_S(g)=1+4ln(g)=0). Bounces to realne zjawisko numeryczne.
- **Sygnał 0/2/3 poniżej krawędzi** — realny, odtwarzalny, w większości odporny na R. Wart dalszego badania.
- **Krawędź −1 jest poprawnie udokumentowana** w kodzie jako "dualism control" — autor solvera (Ty/wcześniejsza sesja) **wiedział**, że F-S siedzi na −1. Ta wiedza po prostu nie przeszła do skryptów weryfikacyjnych tej sesji.

---

## 6. Klasyfikacja wg wagi

| # | Ustalenie | Waga | Status |
|---|-----------|------|--------|
| 1 | "N_neg = f(bounces)" mierzy artefakt pudła (N_neg skaluje z R, próżnia = 19) | **KRYTYCZNE** | Potwierdzone testem |
| 2 | "95% publication-ready / TIER 2 SOLVED / paradigm shift" — nieuprawnione | **KRYTYCZNE** | Narracja > dowody |
| 3 | Korelacja r=−0.9065 to bounces↔lam_min, obie napędzane przez g₀ (wspólna przyczyna) | **WYSOKIE** | Nie dowodzi przyczynowości |
| 4 | Pressure+loops "111%" = skan parametru do celu (overfitting) | **WYSOKIE** | Z kodu |
| 5 | Ekstrakcja ładunków: ansatz sin(r) k=1 przy brzegu pudła | **WYSOKIE** | Metodologia |
| 6 | Status fizyczny modów 0/2/3 (fizyka vs artefakt F-S) — nierozstrzygnięty | **ŚREDNIE** | Otwarte (F-A dual czyste) |
| 7 | sub_tau F-Sp: całkowanie załamane (199 pkt) | **ŚREDNIE** | Nie cytować |
| — | Solver CP-7, ściana-duch, sygnał 0/2/3, zbieżność | **POZYTYW** | Solidne |

---

## 7. Co realnie zostało udowodnione — jedno zdanie

> W formulacji F-S profile μ i τ rozwijają odpowiednio 2 i 3 mody fluktuacji poniżej krawędzi kontinuum (−1), a e żadnego; głębokość rośnie z g₀. Czy to fizyka, czy artefakt formulacji F-S — **nie wiadomo** (dual F-A jest czysty). Wszystko powyżej tego (determinizm N_neg, mechanizm natywny, ciśnienie/pętle stabilizujące, gotowość do publikacji) **nie jest ustalone**.

---

## 8. Skorygowana ścieżka (jeśli chcesz to ratować — a jest co)

Kolejność wg wartości poznawczej:

1. **Rozstrzygnij artefakt vs fizyka (NAJWAŻNIEJSZE).**
   Weź mod μ (λ=−1.282) i τ (najgłębszy λ=−4.216) i policz **to samo w formulacji F-A** dla dualnie odpowiadających tła. Jeśli w F-A te mody znikają → **to artefakt F-S**, cała bounce-hierarchia upada jako "mechanizm fizyczny". Jeśli przeżywają → masz realny wynik. **To jedno pytanie decyduje o wszystkim.**

2. **Popraw licznik.** Wszędzie licz mody **poniżej krawędzi kontinuum** (−1 dla F-S), nie poniżej 0. Przelicz interpolację e→τ tą wielkością — dopiero wtedy "0/2/3 vs g₀" ma sens.

3. **Sprawdź niezależność od g₀ vs bounces.** Zrób punkty gdzie bounces rośnie, ale g₀ prawie stałe (lub odwrotnie), żeby rozdzielić, co naprawdę steruje liczbą modów. Obecna "korelacja" tego nie rozdziela.

4. **Wyrzuć lub przeprojektuj ekstrakcję ładunków.** Ogon fituj z dala od brzegu (r < R/2), z uzasadnioną asymptotyką (najpierw sprawdź, czy ogon jest eksponencjalny czy potęgowy — plot log|g−1| vs r).

5. **Nie wracaj do pressure+loops** dopóki punkt 1 nie rozstrzygnie, że jest CO stabilizować.

---

## 9. Rekomendacja końcowa

**Nie publikować niczego z Tier 2 w obecnej formie.** Twój własny instynkt z sesji #65 ("spokojnie z tą publikacją xD", "czy to natywne efekty czy dopasowania znanych teorii") był trafny — trafniejszy niż dokumenty, które potem wygenerowałem.

Wartość tej pracy jest **realna, ale wąska**: dobrze zbudowany solver + jedno konkretne, odtwarzalne pytanie (mody 0/2/3 poniżej krawędzi F-S), którego status fizyczny wymaga **jednego rozstrzygającego testu** (punkt 8.1). To jest dobra pozycja wyjściowa do rzetelnej pracy — po prostu nie jest to "rozwiązany Tier 2".

Dokumenty do oznaczenia jako **SUPERSEDED / narracyjnie zawyżone**:
- `TIER2_NATIVE_MECHANISM_SUMMARY.md` (twierdzi 95% publication-ready)
- `BOUNCE_HIERARCHY_COMPLETE_EXPOSITION.md` (twierdzi determinizm N_neg)
- `TIER2_COMPLETE_MECHANISM_EXPOSITION.md` (pressure+loops jako mechanizm)

Zachować jako narzędzia/dane:
- `Phase2_bvp_spectrum.py` (solver) — solidny
- `AUDIT_diagnostic.py` (ten test) — punkt odniesienia
- sygnał 0/2/3 poniżej krawędzi — jedyny twardy wynik

---

**Przeprowadził:** Claudian (tryb audytu, Opus)
**Weryfikacja:** własny `AUDIT_diagnostic.py`, uruchomiony na tym samym solverze; liczby w §1–§2 pochodzą z faktycznego przebiegu, nie z wcześniejszej narracji.

---

# DODATEK A — TEST ROZSTRZYGAJĄCY (F-A kanoniczna) — 2026-07-28

**Pytanie z §8.1:** czy mody 0/2/3 to fizyka, czy artefakt formulacji F-S?
**Odpowiedź: ARTEFAKT reprezentacji niekanonicznej. Jednoznacznie.**

## A.1 Rama teoretyczna (z L04)

Cykl `op-L04-ODE-canonicalization-2026-05-04` rozstrzyga autorytatywnie:
> „TGP-canonical α=2 (K = K_geo·φ⁴) JEST jedyną poprawną kanoniczną formulacją.
>  Dualizm jest **pozorny**."

K = φ⁴ ⟹ `F(u) = u⁴` = **F-A**. Formulacja log `F-S` (1+4ln g) jest **niekanoniczna**, a rdzeń ma osobną sekcję `sek08b_ghost_resolution.tex` — czyli ściana-duch to znany problem reprezentacji „do rozwiązania", nie fizyka.

## A.2 Test numeryczny (`AUDIT_FA_dual_test.py`)

Wyprowadziłem ODE solitonu F-A z **tej samej akcji** co F-S
(`F(u)(u''+2u'/r) + ½F'(u)u'² = W'(u)`, F=u⁴, W'=u⁷−u⁶):

```
u'' + 2u'/r = u³ − u² − 2u'²/u,   u(0)=g₀, u'(0)=0, u(∞)→1
```

**Wynik — F-A NIE MA solitonów typu crown w ogóle:**

```
g₀       runaway?   u_end
1.2491 (e)   TRUE    11.56   ← ucieka do nieskończoności
2.0212 (μ)   TRUE    11.49   ← ucieka
3.1891 (τ)   TRUE    11.67   ← ucieka
skan 1.01…3.0: WSZYSTKIE uciekają (runaway=True)
```

Żadne g₀ nie daje związanego solitonu w F-A. Widmo nie istnieje, bo obiekt nie istnieje.

## A.3 Dlaczego (analitycznie, nie tylko numerycznie)

Potencjał F-A: W(u)=u⁸/8 − u⁷/7, W'(u)=u⁶(u−1), **W''(1)=+1 > 0** → próżnia u=1 to **stabilne minimum**. W obrazie mechanicznym soliton rolluje po −W; dla u>1 mamy −W → −∞, więc cząstka startująca w spoczynku przy g₀>1 stacza się do +∞. **Brak rozwiązania powracającego = brak solitonu crown.** To jest strukturalne, nie numeryczne.

Kontrast: F-S ma **W''(1) = −1 < 0** → próżnia **tachionowa**. Stąd krawędź −1, stąd „mody ujemne", stąd w ogóle istnienie crown-lump — ale tylko dzięki **ręcznemu hackowi odbicia** przy ścianie-duchu (kod `soliton_profile` odwraca znak g' gdy g osiąga G_GHOST+0.005; to są właśnie „bounces"). Genuinny gładki soliton nie potrzebuje ręcznych odbić.

## A.4 Werdykt

| | F-A (kanoniczna, α=2 φ⁴) | F-S (log, niekanoniczna) |
|---|---|---|
| Próżnia | stabilna (W''=+1) | tachionowa (W''=−1) |
| Ściana-duch | brak | tak (G_GHOST=e^−¼) |
| Soliton crown g₀>1 | **nie istnieje (runaway)** | istnieje tylko przez hack odbicia |
| Mody 0/2/3 pod krawędzią | N/A (brak obiektu) | tak — ale to widmo hacka |

**Cała bounce-hierarchia (bounces, mody 0/2/3, „hierarchia generacji z odbić") jest własnością obiektu, który istnieje wyłącznie w niekanonicznej reprezentacji F-S i tylko dzięki ręcznej procedurze odbicia od osobliwości log-formy. W kanonicznej formulacji TGP obiekt ten nie istnieje. Mechanizm NIE jest fizyczny.**

## A.5 Zakres (czego to NIE przekreśla)

- **Nie przekreśla** pracy nad **masami leptonów** (`m_obs = c·A²·g₀^[e²(1−α/4)]`, <0.1% PDG). Ta formuła używa g₀ i A_tail (projekcja asymptotyczna) i **nie zależy** od struktury saddle/spektralnej. Może być poprawna niezależnie od tego wyniku.
- **Nie przekreśla** solvera CP-7 ani całego TGP.
- **Przekreśla** konkretnie: Tier-2 „native bounce-hierarchy mechanism" jako twierdzenie fizyczne.

**Wartość wyniku:** negatywny, ale mocny — zamyka fałszywy trop definitywnie, zamiast zostawiać go jako „95% publication-ready". To jest dobry wynik audytu.
