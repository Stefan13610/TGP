---
title: "Phase_method_decisions — decyzje implementacyjne ZAMROŻONE przed startem obliczeń (forma EL 1D + rozjazd z #63, relaksacja Newton na pół-komórce, dyskretyzacja Blocha, solver, metryki błędu, reguła modu translacyjnego)"
date: 2026-08-31
type: method-decisions
tgp_owner: research/op-bloch-chain-stability-2026-08-31
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase3_two_sectors.py]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem)

Wszystkie poniższe decyzje zapisane przed uruchomieniem pierwszego skryptu.
Kryteria/progi/punkty LOCKa (Phase0_balance.md) NIETKNIĘTE — poniżej wyłącznie
doprecyzowania implementacyjne, które LOCK jawnie delegował, oraz definicje
metryk tam, gdzie LOCK ich nie domknął (każda z osiągalnym FAIL).

## 1. Forma EL 1D (M-C) i ROZJAZD z #63 (`Phase3_nonlinear_evolution.py`)

**M-C (LOCK, dosłownie):** `g″ + (2/g)g′² = g²(1−g)`, β=γ=1.

**Weryfikacja wariacyjna (algebra ręczna, jeden krok):** M-C jest DOKŁADNIE
równaniem EL akcji kanonicznej Noty kanonicznej 2026-04-07 (ANALIZA_N2),
zredukowanej do 1D:

```
S[g] = ∫ dt dx [ ½K(g)(g_t² − g_x²) − U(g) ],  K = g⁴,  U = βg⁷/7 − γg⁸/8
EL:   K g_tt + ½K′g_t² = K g_xx + ½K′g_x² − U′
statyka: K g″ + ½K′g′² = U′  ⟹  g″ + (2/g)g′² = U′/K = βg² − γg³ = g²(1−g) ✓
```

**ROZJAZD z #63 (udokumentowany, wymagany przez LOCK M-C):**
`Phase3_nonlinear_evolution.py` (#63, M0-f_ε) używa INNEGO reprezentanta
sektora: kinetyka F = f_ε(g) z f = 1+2α·ln g (α=2), potencjał W = g³/3−g⁴/4
z wagą 1 (W′ = g²(1−g)). Jego statyka to `f g″ + ½f′g′² = g²(1−g)`, co NIE
jest formą M-C (kinetyka f_log vs K=g⁴; W waga-1 vs U′ = K·g²(1−g)).
Obie formy mają IDENTYCZNĄ linearyzację wokół próżni g=1 (m² = −1,
ω²(k) = k²−1; sprawdzone: W″(1) = U″(1)/K(1) = −1) — rozjazd dotyczy
wyłącznie skończonych amplitud.

**DECYZJA (frozen):** cały cykl pracuje w JEDNYM spójnym modelu wariacyjnym —
akcji kanonicznej powyżej (źródle M-C). Powód: tło z Phase 2 jest wtedy
DOKŁADNIE stacjonarnym punktem dynamiki Phase 4 i linearyzacji Phase 3
(centralna lekcja Q1 poprzednika: tło niestacjonarne w danej dynamice
„pada samo z siebie" i uniemożliwia separację efektu). Mieszanie tła M-C
z dynamiką dosłownego M0-f_ε (#63) byłoby niespójne wariacyjnie — to jest
ten rozjazd i tak go rozstrzygamy.

**Rola f_ε (LOCK M-D: „ta sama regularyzacja f_ε, ε=0.2; kontrola 0.1"):**
przenosimy z #63 (a) strukturę numeryczną ewolucji (semi-dyskretny układ
hamiltonowski, RK4, gate energii, detektor breakdown g≤0), (b) samą MAPĘ
regularyzującą x ↦ ½(x + √(x²+ε²)), zastosowaną do funkcji kinetycznej TEGO
sektora: `K_ε(g) = ½(K + √(K²+ε²))`, K = g⁴; ε = 0.2 (kontrola 0.1) [INPUT
z #63]. K_ε użyte WYŁĄCZNIE w ewolucji nieliniowej (P1b, Phase 4) — pełni
tę samą rolę co f_ε w #63 (skończona kinetyka przy g→0; dla g≈1 odchyłka
K_ε−K ≤ ε²/8 ≈ 5e−3 bez wpływu na detekcję ucieczki). W operatorach
SPEKTRALNYCH (P1a, Phase 3) K bez regularyzacji: K = g⁴ > 0 na tłach
(regularyzacja zbędna), a K_ε(1) = 1.00995 przesuwałoby ω² próżni o ~1%
i sztucznie łamało bramkę P1a 1e−3 — spektra liczone dla ścisłej akcji.

## 2. Operator fluktuacji (druga wariacja; do Phase 1/3)

Linearyzacja g = g₀(x) + φ(x)e^{−iωt} wokół statycznego g₀ (rachunek ręczny,
niżej zweryfikowany bramką P1a i kontrolami P3b/P3c):

```
−(K φ′)′ + Q φ = ω² K φ,      K = g₀⁴  (waga; problem uogólniony)
Q = U″(g₀) − K′(g₀)g₀″ − ½K″(g₀)g₀′²
  = K·(2g₀ − 3g₀² + 2g₀′²/g₀²)   [po podstawieniu EL za g₀″; forma używana]
```

Próżnia: Q/K = −1 ⟹ ω² = k²−1 (tachion) ✓. Struktura identyczna ze
zwalidowanym operatorem Q2 poprzednika (waga w=ψ⁴, tam sektor stabilny).

**Cross-check kanoniczny (dodatkowa kontrola wewnętrzna, nieusuwana):**
podstawienie u = g³/3 kanonizuje kinetykę (½K g′² = ½u′², ½K g_t² = ½u_t²)
i daje równoważny problem wagi 1: `−χ″ + (4g₀ − 5g₀²)χ = ω²χ`, χ = g₀²φ
(próżnia: 4−5 = −1 ✓). W Phase 3 dla każdego d porównujemy λ_min obu form
w k∈{0, π/d} (raport zgodności; rozjazd rzędu >O(h²) = czerwona flaga
implementacji).

**Mod translacyjny (zapisane PRZED obliczeniami):** tło M-C nie ma
zewnętrznego źródła, więc translacja jest symetrią spontanicznie złamaną —
φ = g₀′ jest DOKŁADNYM modem zerowym w k=0 (kontinuum). Konsekwencja
strukturalna: min_k λ_min ≤ 0 analitycznie dla każdego tła niestałego.
Reguła pre-rejestrowana: w k=0 identyfikujemy mod translacyjny przez overlap
B-ważony z g₀′ (|⟨φ, g₀′⟩_B|/(‖φ‖_B‖g₀′‖_B) ≥ 0.9); jego λ raportujemy
osobno (oczekiwane |λ_trans| = O(h²)). Kryterium LOCKa ω²_min(d) > 0
stosujemy LITERALNIE do pełnego minimum (z modem translacyjnym włącznie);
dodatkowo raportujemy deskryptywnie (nie-kryterialnie) minimum po widmie
z wyłączonym modem translacyjnym. Żadnego progu nie zmieniamy.

## 3. Schemat relaksacji (Phase 2): NEWTON tłumiony na pół-komórce

- **Wybór: Newton** z line-search (wzorzec `Phase3_two_sectors.py`
  poprzednika), nie gradient flow — szybsza i ostrzejsza zbieżność do
  ‖R‖∞ ≤ 1e−10.
- **Pół-komórka + parzystość:** brak źródła ⟹ Jacobian pełnej komórki
  periodycznej ma dokładne zero translacyjne (Newton osobliwy). Relaksujemy
  na [0, d/2] w klasie parzystej (węzły x_j = j·h, j = 0..N/2, h = d/N;
  odbicie lustrzane na obu końcach: g₋₁ = g₁, g_{N/2+1} = g_{N/2−1}),
  bump wycentrowany w d/2. Po zbiegnięciu lustrzane rozwinięcie do pełnej
  komórki (N węzłów) i weryfikacja kryterium LOCKa ‖residuum EL‖∞ ≤ 1e−10
  NA PEŁNEJ KOMÓRCE periodycznej — to jest liczba raportowana.
- **Residuum:** R = D2(g) + 2·D1(g)²/g − g²(1−g), różnice centralne rzędu 2.
  Jacobian trójdiagonalny: diag = −2/h² − 2D1²/g² − (2g−3g²),
  up/lo = 1/h² ± 2D1/(g·h); na końcach lustrzanych podwojone sprzężenie
  2/h². Solver: `scipy.linalg.solve_banded`.
- **Podłoga arytmetyki** (por. korekta 2 poprzednika): szacunek szumu
  zaokrągleń członu D2 to ~4·eps·max(g)/h²; dla najgęstszej siatki
  nietrywialnej (d=3π, N=800) ≈ 2e−11 < 1e−10 — próg LOCKa osiągalny.
  Raportujemy osiągnięte residuum; stagnację powyżej 1e−10 raportujemy
  jako FAIL punktu (bez zmiany progu).
- **Starty (zalockowane 2):** g₀(x) = 1 − A·exp(−(x−d/2)²/(2σ²)),
  A ∈ {0.7 „deep", 0.3 „shallow"}; **σ = d/8** (decyzja frozen: szerokość
  skalowana z komórką). Oba starty raportowane zawsze.
- **Start pomocniczy S3 (dozwolony dodatek, uruchamiany TYLKO gdy oba
  zalockowane starty nie dadzą rozwiązania niestałego dla danego d):**
  profil z kwadratury pierwszej całki E_mech = ½g⁴g′² − U(g) = const
  (bisekcja E do okresu d). Starty zalockowane NIE są usuwane z raportu.
- **Siatki:** N ∈ {400, 800} (pełna komórka). Transfer: rozwiązanie N=400
  interpolowane (CubicSpline) na siatkę N=800 jako start Newtona.
  Zbieżność siatkowa ‖g_N − g_2N‖∞ liczona NA WĘZŁACH WSPÓLNYCH
  (węzły N=400 ⊂ N=800; siatka węzłowa wybrana m.in. po to).
- **Klasyfikacja punktu:** NIESTAŁE-ZBIEŻNE (‖g−1‖∞ ≥ 0.05,
  ‖g_N−g_2N‖∞ ≤ 1e−4, ‖R‖∞ ≤ 1e−10) / KOLAPS-DO-PRÓŻNI (‖g−1‖∞ < 0.05) /
  UCIECZKA-g→0 / NIEZBIEŻNE. Diagnostyka raportowana: g_min, g_max,
  liczba minimów lokalnych, rozrzut E_mech wzdłuż profilu.
- Uwaga strukturalna (zapisana przed startem, do uczciwej interpretacji):
  pierwsza całka implikuje, że okres małych oscylacji wokół g=1 wynosi
  dokładnie 2π i rośnie z amplitudą (dywergencja przy separatrysie g→0);
  możliwy więc brak rozwiązań niestałych dla d ≤ 2π — będzie to raportowane
  jako wynik, nie artefakt.

## 4. Dyskretyzacja Blocha (Phase 1/3) i solver własny

- Pełna komórka, węzły x_j = j·h, j = 0..N−1, h = d/N. Warunek Blocha
  φ(x+d) = e^{ikd}φ(x) ⟹ macierz Hermitowska z fazą na sprzęgle brzegowym:

```
A(k)[j,j]     = (K_{j+½} + K_{j−½})/h² + Q_j
A(k)[j,j+1]   = −K_{j+½}/h²                    (j = 0..N−2)
A(k)[N−1,0]   = −K_{N−½}/h² · e^{+ikd}         (narożnik; h.c. w [0,N−1])
B             = diag(K_j)                       (waga)
K_{j+½} = ½(K_j + K_{j+1 mod N})  (średnia węzłowa, jak w Q2 poprzednika)
Q_j = K_j·(2g_j − 3g_j² + 2·D1_j²/g_j²),  D1 = różnica centralna periodyczna
```

- **Problem uogólniony** A(k)φ = ω²Bφ symetryzowany: H(k) = B^{−½}A(k)B^{−½}
  (Hermitowska, zespolona). **Solver własny: `scipy.linalg.eigh`**
  (dense, LAPACK zheevr; N ≤ 800 — koszt pomijalny; scipy dostępny:
  sprawdzono importem, 1.17.1).
- Siatka k: 16 punktów równomiernie z KRAŃCAMI w [0, π/d] (kontrola 32) —
  jak zalockowano.
- Sektor stabilny (P1a-kontrola, P3c): operator poprzednika verbatim
  (waga w=ψ⁴, potencjał −w·(2ψ−3ψ²−2ψ′²/ψ²); próżnia: ω² = k²+1) —
  z tą samą fazą Blocha.

## 5. Metryki błędu i bramki (doprecyzowania frozen; progi LOCKa bez zmian)

- **P1a:** komórka d = 2π, siatki N ∈ {400, 800} (h i h/2); 16 punktów k
  w [0, π/d]; porównanie 3 najniższych gałęzi z dokładnymi zwiniętymi
  parabolami ω²_ex = (k+2πn/d)² + m², m² = −1 (tachion) i +1 (kontrola
  stabilna). **Błąd względny := |ω²_num − ω²_ex| / max(|ω²_ex|, 1)**
  (podłoga 1 = |m²|, naturalna skala; konieczna, bo gałąź tachionowa
  przechodzi przez ω² = 0 — bez podłogi błąd względny w zerze jest
  niezdefiniowany). Gate: max błąd ≤ 1e−3 na siatce bazowej (N=400)
  ORAZ stosunek maxbłędu(N=400)/maxbłędu(N=800) ∈ [3, 5] („~×4", rząd 2).
  FAIL osiągalny: zły znak/operator/faza ⟹ błąd O(1).
- **P3b:** ta sama metryka błędu i próg (≤1e−3) dla 8 najniższych gałęzi
  próżni w superkomórce 4d, siatka k: 8 punktów w [0, π/(4d)], dla każdego
  zalockowanego d.
- **P1b:** komórka d=2π, N=400, K_ε (ε=0.2 oraz kontrola 0.1);
  start g = 1 + a(0.5 + cos x), a = 1e−3; t_end = 4.0; RK4,
  dt ∈ {0.004, 0.002}. Gate energii |ΔE|/E ≤ 1e−6 przez cały bieg
  (H semi-dyskretne zachowane dokładnie przez ODE ⟹ gate mierzy czysty
  błąd integratora — jak #63). Zbieżność dt (×2): ‖g_dt(t_end) −
  g_dt/2(t_end)‖∞ ≤ 1e−6 (pre-rejestrowane; FAIL osiągalny).
- **Phase 3 zbieżność (LOCK):** wartość główna ω²_min := (N=800, k=32);
  Δ_siatka = |ω²_min(N400,k32) − ω²_min(N800,k32)|,
  Δ_k = |ω²_min(N800,k16) − ω²_min(N800,k32)|;
  ZBIEŻNE ⟺ obie ≤ 0.01·max(|ω²_min|, 0.1). Jeśli dwa zalockowane starty
  Phase 2 dają RÓŻNE tła dla tego samego d — spektrum liczone i raportowane
  dla obu (zakaz wyboru post-hoc).
- **P3c:** konstrukcja Q2 poprzednika verbatim: sektor stabilny, komórka
  d_s = 4.0, σ_s = 0.5, q = 1.0, Newton periodyczny (źródło łamie
  translację — brak modu zerowego), N ∈ {400, 800}, Bloch 16 k. Gate:
  ω²_min > 0. Kotwica krzyżowa (raport): λ_min(k=0) vs +1.88310
  (Phase_FINAL_close poprzednika).

## 6. Phase 4 (warunkowa) — doprecyzowania frozen

- Superkomórka 4d, N_sc = 4·400 = 1600 (h = d/400); tło = tło Phase 2
  (N=400) powielone ×4.
- **Mod perturbacji:** mod minimalny z Phase 3 w najbliższym argmin
  punkcie k WSPÓŁMIERNYM z superkomórką (k ∈ {0, π/(2d), π/d});
  re-diagonalizacja w tym k; perturbacja δg = ±a·Re(φ_k rozszerzone
  fazowo na 4 komórki), znormalizowane ‖Re φ‖∞ = 1. Jeżeli argmin
  jest modem translacyjnym (reguła §2) — używamy modu najniższej
  wartości NIE-translacyjnej (perturbacja wzdłuż czystej translacji
  jest trywialna: przesuwa tło).
- **Amplituda:** a = 0.01·‖tło‖, ‖tło‖ := max_x|g₀(x) − 1| (jak bg_amp
  w #63).
- **Ewolucja:** RK4 na semi-dyskretnym H z K_ε; ε = 0.2 oba znaki ±;
  kontrola ε = 0.1 przy znaku +; kontrola dt/2 przy (+, ε=0.2);
  dt = 0.004; t_end = 3·t*_ref = 10.86 (t*_ref = 3.62 [INPUT #63]).
- **Ucieczka :=** (g ≤ 0 lub niefinityczne) LUB max_x|g(x,t) − g₀(x)| >
  1.0·‖tło‖ (100% tła). Raport t_escape vs 2·t*_ref = 7.24
  i 3·t*_ref = 10.86. Gate energii raportowany (overflow po breakdown
  raportowany wprost, jak u poprzednika).
- Phase 4 uruchamiana dla KAŻDEGO d ze zbieżnym Phase 3 (warunek LOCKa
  to zbieżność, nie znak).

## 7. Higiena wykonania

- Python: `python` (CPython 3.14.2; numpy 2.4.3, scipy 1.17.1 —
  sprawdzone importem PRZED startem).
- Zero zmiany cwd; wszystkie ścieżki pełne; po każdym zapisie weryfikacja
  `ls`. Outputy: `Phase*_output.txt` (stdout skryptu przekierowany).
- Tła Phase 2 zapisywane do `Phase2_backgrounds.npz` (artefakt
  obliczeniowy dla Phase 3/4; wymienione w close-nocie).
- Rdzeń `.tex`, STATE.md, git — NIETYKANE. Katalogi innych cykli —
  NIETYKANE.

**FROZEN 2026-08-31, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
