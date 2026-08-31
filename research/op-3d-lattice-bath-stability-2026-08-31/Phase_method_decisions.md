---
title: "Phase_method_decisions — decyzje implementacyjne ZAMROŻONE przed startem obliczeń (rozdział modeli P1 vs P2–4, operator 3D kartezjański, Newton–Krylov, Bloch 3D, eigsh 'SA', interpolacja radialna→3D, metryki bramek)"
date: 2026-08-31
type: method-decisions
tgp_owner: research/op-3d-lattice-bath-stability-2026-08-31
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_method_decisions.md]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase2_bath_runaway.py]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem cyklu)

Wszystkie decyzje poniżej zapisane przed uruchomieniem pierwszego skryptu.
Kryteria/progi/punkty LOCKa (Phase0_balance.md) NIETKNIĘTE — poniżej wyłącznie
doprecyzowania, które LOCK jawnie delegował („wybierz i zamroź"), oraz
definicje metryk tam, gdzie LOCK ich nie domknął (każda z osiągalnym FAIL).

**Rejestr WEJŚĆ (flagowany w każdym zależnym wyniku):** g₀=2.02117
(kalibracja μ #63), ε=0.2 (kontrola 0.1) [INPUT #63], d*₁(μμ, M-P)=3.0790
[INPUT bath-two-sectors Phase 1], β=γ=1.

## 0. Rozdział modeli: bramka P1 (#63 M0-f_ε) vs rachunek P2–P4 (akcja kanoniczna)

Kotwice LOCKa P1b/P1c (λ_min(w1)=−1.3896, t*=3.62) istnieją WYŁĄCZNIE
w modelu #63 M0-f_ε (kinetyka F=f_ε∘f_log, f=1+2α·ln g, α=2; potencjał
W=g³/3−g⁴/4 z wagą 1). Dlatego:

- **P1b (kotwica radialna) i P1c (most radialny→kartezjański) liczone
  w DOSŁOWNYM modelu #63 M0-f_ε** — bramka sprawdza adekwatność
  dyskretyzacji 3D przez reprodukcję kotwic w ich własnym modelu
  (profil, operator wagi-1, dynamika hamiltonowska — konwencje verbatim
  z `Phase3_nonlinear_evolution.py`/`Phase2_bath_runaway.py`, w tym
  korekta 1 poprzednika: detektor breakdown g ≤ 0.01 co krok RK4).
- **P1a, Phase 2, Phase 3, Phase 4 — w akcji kanonicznej** (LOCK §2 M-3D):
  S = ∫dt d³x [½K(g)(g_t²−|∇g|²) − U(g)], K=g⁴, U=βg⁷/7−γg⁸/8, β=γ=1.
  Statyka EL: ∇²g + (2/g)|∇g|² = g²(1−g). Uzasadnienie spójności
  tło↔linearyzacja↔dynamika: bloch-chain method_decisions §1 —
  **przyjęte przez LOCK dosłownie** („uzasadnienie jak w bloch-chain
  method_decisions — przyjęte"). Regularyzacja: mapa x↦½(x+√(x²+ε²))
  zastosowana do K (K_ε), WYŁĄCZNIE w ewolucji nieliniowej (P1d-wariant
  kanoniczny, Phase 4); spektra bez regularyzacji (LOCK §2).
- Obie rodziny mają IDENTYCZNĄ linearyzację próżni g=1 w granicy ε→0
  (ω²=|k|²−γ); w modelu #63 przy ε=0.2 waga próżni F_ε(1)=1.00990
  (przesunięcie ~1% dyspersji) — dlatego bramka P1a (exact) liczona jest
  dla operatora kanonicznego (K(1)=1 ściśle), a most P1c dla operatora #63
  (kotwica −1.3896 zawiera F_ε). Zapisane tu PRZED obliczeniem.

## 1. Operator fluktuacji 3D (forma FROZEN)

**Sektor tachionowy (kanoniczny; Phase 1a/3):** druga wariacja wokół
statycznego g₀ (algebra jak bloch-chain §2, człony 1D→3D: g′²→|∇g₀|²,
g″→∇²g₀; po podstawieniu EL):

```
−∇·(K∇φ) + Qφ = ω² K φ,   K = g₀⁴  (problem uogólniony, waga B=K)
Q = K·(2βg₀ − 3γg₀² + 2|∇g₀|²/g₀²)
```

Próżnia: Q/K = −1 ⟹ ω² = |k|²−γ ✓ (bramka P1a).

**Cross-check kanoniczny (kontrola wewnętrzna, nieusuwana):** u=g³/3,
χ=g₀²φ ⟹ problem wagi 1: `−∇²χ + (4βg₀−5γg₀²)χ = ω²χ` (próżnia 4−5=−1 ✓).
Porównanie λ_min obu form w k∈{Γ, R} dla każdego istniejącego tła;
rozjazd ≫O(h²) = czerwona flaga implementacji.

**Sektor stabilny (P3c; konstrukcja Q2 poprzednika verbatim → 3D):**
w=ψ⁴, pot = −w·(2βψ − 3γψ² − 2|∇ψ|²/ψ²); próżnia: ω² = |k|²+γ ✓.

**Operator #63-model 3D (P1c):** waga 1, `−∇·(F_ε∇v) + Q₆₃ v = λ v`,
Q₆₃ = W″(g) − ½F_ε″(g)|∇g|² − F_ε′(g)·∇²g (dokładnie #63 build_tridiag_gen
w 3D; ε=0.2 jak kotwica). Dla profilu radialnego |∇g|²=g′(r)²,
∇²g = g″+2g′/r liczone ANALITYCZNIE ze splajnów radialnych (nie z różnic
na siatce 3D) — profil jest dokładnie radialny; w r=0: ∇²g = 3g″(0).

**Dyskretyzacja (wszystkie operatory):** siatka węzłowa x_j = j·h,
j=0..N−1 na wymiar, h=d/N (soliton/źródło w centrum d/2 — węzeł dokładny
dla N parzystych); stencil 7-punktowy; współczynnik kinetyczny na linkach:
K_{j+½} = K(½(g_j+g_{j+1}))(średnia węzłowa pola, jak F_mid #63 i jak Q2);
diag = Σ_kier (K_{+½}+K_{−½})/h² + Q; off = −K_{±½}/h². |∇g₀|² i ∇²g₀
w Q z centralnych różnic periodycznych (Phase 3; tła dyskretne) albo
z radialnych splajnów (P1c; profil analityczny).

## 2. Bloch 3D i solver spektralny (FROZEN)

- Warunek Blocha φ(x+d·e_α) = e^{ik_α d}φ(x): faza na sprzęgle owijającym
  w kierunku α. Zalockowany zbiór k: Γ=(0,0,0), X=(π/d,0,0), M=(π/d,π/d,0),
  R=(π/d,π/d,π/d) — we WSZYSTKICH tych punktach e^{ik_α d} ∈ {+1,−1}
  ⟹ **macierze RZECZYWISTE symetryczne** (float64; implementacja fazy
  jako ±1 na wrapie). Problem uogólniony symetryzowany:
  H(k) = B^{−½}A(k)B^{−½} (skalowanie wierszy/kolumn macierzy rzadkiej csr).
- **Solver: `scipy.sparse.linalg.eigsh(H, k=10, which='SA')`** (FROZEN;
  bez shift-invert — faktoryzacja LU 3D przy N=48..100 na wymiar
  przekracza pamięć). Parametry FROZEN: tol=1e−6, ncv=80,
  maxiter=200000, v0 = deterministyczny wektor pseudolosowy
  (np.random.default_rng(20260831), rozkład normalny, znormalizowany) —
  celowo NIE symetryczny: start symetryczny leżałby w podprzestrzeni
  niezmienniczej operatora i mógłby ominąć mody antysymetryczne
  (w tym translacyjne, które Phase 3 MUSI zidentyfikować).
  Brak zbieżności ARPACK → punkt raportowany INCOMPLETE
  (bez zmiany siatek).
- ≥10 najniższych wartości własnych wszędzie (LOCK).

## 3. Interpolacja radialna → 3D (FROZEN)

- Profil radialny #63: `profile_soft(g₀=2.02117, R=60, N=4000)` verbatim
  (siatka przesunięta r_j=(j+½)h, h=0.015); tablice g, g′, g″ z DOP853.
- CubicSpline osobno dla (g−1), g′, g″ na siatce radialnej; ewaluacja
  w r(x) = odległość minimalnego obrazu od centrum komórki/pudła;
  dla r > r_end: g:=1, g′:=g″:=0 (zero-padding do próżni — LOCK P1c);
  dla r < r₀=0.0075: ekstrapolacja splajnu (profil gładki, g′(0)=0),
  w węźle centralnym r=0: ∇²g := 3·g″(0⁺).
- **Start relaksacji Phase 2 (PRIMARY, zalockowany):** superpozycja
  periodyzowana g_start(x) = 1 + Σ_{m∈{−1,0,1}³} [g_rad(|x−c−m·d|) − 1]
  (27 obrazów; ogony 1/r → suma obrazów 3×3×3 wystarczająca — wskazówka
  HANDOFF; max odległość ≈ 2.6·d ≤ 33 < 60 = zasięg profilu ✓).
  **Start 2 (zalockowany):** 1 + 0.7·(ta sama suma) (0.7× amplitudy).
  Oba starty dokładnie symetryczne względem odbić x_α → d−x_α
  (węzeł j → (N−j) mod N).

## 4. Schemat relaksacji Phase 2 (FROZEN): scipy.optimize.newton_krylov

- **Wybór: `scipy.optimize.newton_krylov`** (nie gradient flow — tło jest
  punktem SIODŁOWYM energii w sektorze tachionowym: U″(1)=−γ<0; flow
  gradientowy z definicji ucieka z siodła; zapisane przed startem).
- Residuum (forma FROZEN, jak bloch-chain 1D→3D):
  `R(g) = Σ_α D2_α g + (2/g)·Σ_α (D1_α g)² − βg²(1−g)·(…)` dokładnie:
  R = ∇²_h g + (2/g)|∇_h g|² − g²(β−γg) z β=γ=1, różnice centralne
  rzędu 2, periodyczne (np.roll).
- Parametry FROZEN: f_tol=1e−9 (max-norm wewnątrz scipy), method='lgmres',
  inner_maxiter=200, maxiter=300; preconditioner inner_M = spilu
  macierzy (∇²_h − I) (csc; drop_tol=1e−5, fill_factor=15) — SPD po
  negacji, budowany raz per (d,N).
- Po zbiegnięciu: **weryfikacja kryterium LOCKa ‖R‖∞ ≤ 1e−8 na pełnej
  komórce** (liczba raportowana) + raport asymetrii
  max|g − g_odbite| (diagnostyka; translacyjne zero łamane siatkowo
  do O(h²), start symetryczny — dryf raportowany, nie korygowany).
- Podłoga arytmetyki: szum D2 ≈ 4·eps·max(g)/h²; najgęstsza siatka
  d=π, N=48 (h=0.0654): ≈ 4e−13 ≪ 1e−8 — próg LOCKa osiągalny.
- **Transfer siatek:** N=32 relaksowane ze startu; N=48 startowane
  z rozwiązania N=32 interpolowanego trygonometrycznie (zero-padding FFT;
  pole periodyczne gładkie). Gdy N=32 dał KOLAPS — N=48 startowany
  z surowego startu na własnej siatce (punkt raportowany niezależnie).
- **Zbieżność siatkowa (LOCK ≤5e−3):** ‖g₃₂−g₄₈‖∞ na WSPÓLNEJ podsiatce
  węzłowej 16³ (węzły x = m·d/16 należą dokładnie do obu siatek —
  bez interpolacji; to jest „wspólna siatka" LOCKa).
- **Klasyfikacja punktu (d, start):** NIESTAŁE-ZBIEŻNE (‖g−1‖∞ ≥ 0.05
  ∧ ‖R‖∞ ≤ 1e−8 ∧ ‖g₃₂−g₄₈‖∞ ≤ 5e−3) / KOLAPS-DO-PRÓŻNI (‖g−1‖∞ < 0.05) /
  UCIECZKA-g→0 (min g < 0.01) / NIEZBIEŻNE (reszta). Dedup teł:
  ‖Δ‖∞ < 1e−8. Diagnostyka: g_min, g_max, ‖g−1‖∞, asymetria.

## 5. Metryki bramek Phase 1 (doprecyzowania FROZEN; progi LOCKa bez zmian)

- **P1a (dyspersja próżni 3D exact):** pudło L=2π (dozwolone przez LOCK),
  N ∈ {32, 64}; operator kanoniczny (i kontrola stabilna) na g≡1;
  zalockowany zbiór k (d=L=2π ⟹ X=(½,0,0) itd.). **Bramka liczona na
  NAJNIŻSZEJ gałęzi w każdym punkcie k** (wielkość zasilająca Phase 3 to
  ω²_min; LOCK nie domyka liczby gałęzi — definicja FROZEN TERAZ):
  błąd := |ω²_num − ω²_ex|/max(|ω²_ex|,1), ω²_ex = min_G |k+G|² ∓ γ
  (podłoga 1 = |m²| jak u poprzednika — gałąź przechodzi przez 0).
  Gate: max po 4 punktach k i obu sektorach ≤ 1e−3 przy N=32 ORAZ
  stosunek maxerr(N=32)/maxerr(N=64) ∈ [3,5] na sektor (rząd 2; w Γ
  najniższa gałąź jest reprezentowana dokładnie — stąd stosunek liczony
  na maksimum po zbiorze k, nie per-punkt). DODATKOWO raport deskryptywny:
  10 najniższych gałęzi w każdym k z tą samą metryką (informacyjnie).
  FAIL osiągalny: zły znak/stencil/faza ⟹ błąd O(1).
- **P1b (kotwica radialna, model #63 verbatim):** R=60, N=4000;
  λ_min(w1) gate −1.3896 ± 1e−3; ewolucja radialna a=±0.01,
  dt ∈ {0.004, 0.002}, gate energii ≤1e−6, detektor g≤0.01 co krok;
  KAŻDY bieg musi dać t* ∈ 3.62 ± 0.05; t*_radial := min po biegach.
- **P1c (most, model #63):** pudło L=30 periodyczne (LOCK L≥30; przekątna
  półpudła 26 < 60 ⟹ zero-padding aktywny tylko formalnie), soliton
  w centrum. Siatki: N=76 (h=0.3947 ≈ 0.4 — „N odpowiadający h≈0.4")
  oraz zagęszczenie N=100 (h=0.30). Gate: |λ_min(3D, N=76) − λ_rad|
  ≤ 0.05·|λ_rad| (λ_rad z P1b) ORAZ |λ_min(3D, N=100) − λ_rad| <
  |λ_min(3D, N=76) − λ_rad| (poprawa z zagęszczeniem — literalnie).
  t*_izo(3D): ewolucja #63-3D na siatce N=76, start g₃D + a·v₃D
  (v₃D = mod najniższy z eigsh, znormalizowany ‖v‖_{L²(h³)}=1, znak:
  max|v| dodatni — konwencja #63), a=±0.01, dt ∈ {0.02, 0.01}
  (CFL: dt ≪ h=0.395), detektor g≤0.01 co krok; KAŻDY bieg
  t* ∈ 3.62 ± 15%; **t*_izo(3D) := min po biegach** (konwencja P2a) —
  wielkość zasilająca Phase 4.
- **P1d (gate energii 3D):** pudło L=2π, N=32, perturbacja
  g = 1 + a(0.5 + cos x + cos y + cos z), a=1e−3, t_end=4.0, RK4,
  dt ∈ {0.004, 0.002}; OBIE dynamiki (kanoniczna K_ε ε∈{0.2,0.1}
  i #63-owa F_ε ε=0.2): gate |ΔE|/E ≤ 1e−6 przez cały bieg;
  zbieżność dt: ‖g_dt(t_end) − g_{dt/2}(t_end)‖∞ ≤ 1e−6 (FAIL osiągalny).
- Dynamika 3D (obie rodziny): semi-dyskretny hamiltonian
  H = Σ_j h³[π_j²/(2Φ(g_j)) + P(g_j)] + ½Σ_links h³ Φ(g_mid)(Δg/h)²,
  (Φ,P) = (F_ε, W) dla #63 lub (K_ε, U) dla kanonicznej; g_mid = średnia
  węzłowa na linku; periodycznie; RK4; dH/dt=0 dokładnie dla ODE ⟹ gate
  mierzy czysty błąd integratora (wzorzec #63).

## 6. Phase 3 — doprecyzowania FROZEN

- Dla każdego istniejącego tła (po dedup): eigsh w 4 zalockowanych
  punktach k, N ∈ {32,48}, 10 wartości + wektory.
- **Identyfikacja modów translacyjnych (kontrola nieusuwalna, PRZED
  interpretacją):** T = span{D1_x g₀, D1_y g₀, D1_z g₀} ortonormalizowany
  w iloczynie B=K (Gram); mod φ uznany za translacyjny ⟺
  ‖P_T φ‖_B/‖φ‖_B ≥ 0.9 (reguła FROZEN). Raport: λ modów translacyjnych
  (oczekiwane O(h²) w Γ), overlapy wszystkich 10 modów.
- **ω²_min(d) := najniższa wartość własna po 4 punktach k PO ODJĘCIU
  modów zidentyfikowanych jako translacyjne** (litera LOCKa §3 Phase 3);
  pełne minimum (z translacjami) raportowane obok deskryptywnie.
- Zbieżność (LOCK): |ω²_min(N32) − ω²_min(N48)| ≤ 0.05·max(|ω²_min|,0.1),
  wartość główna := N=48.
- Cross-check u-formy w {Γ, R} dla każdego tła (N=32 — koszt), raport |Δ|.
- **P3b:** próżnia g≡1 w TEJ SAMEJ komórce [0,d)³, N ∈ {32,48}, 4 punkty k,
  10 najniższych gałęzi vs posortowane |k+G|²−γ (G = (2π/d)·Z³, zakres
  |G_α/(2π/d)| ≤ 6 wystarczający dla 10 gałęzi); metryka
  |Δ|/max(|ex|,1) ≤ 1e−2 (próg LOCKa; podłoga jak w P1a).
- **P3c:** sektor stabilny 3D ze źródłem gaussowskim: komórka d_s=4.0,
  σ_s=0.5, q=1.0 (konstrukcja Q2/bloch-chain → 3D), źródło periodyzowane
  ρ = Σ_{m∈{−2..2}³} exp(−|x−c−m·d_s|²/(2σ_s²)) / ((2π)^{3/2}σ_s³);
  EL: ∇²ψ + (2/ψ)|∇ψ|² + βψ² − γψ³ + qρ = 0 (1D-wzorzec newton_stab
  verbatim → 3D), relaksacja newton_krylov (źródło łamie translację),
  kontynuacja q we frakcjach {0.25,0.5,0.75,1.0}; N ∈ {32,48}; Bloch
  4 punkty k, operator sektora stabilnego (§1). Gate: ω²_min > 0 na obu
  siatkach (osiągalny FAIL maszynerii).

## 7. Phase 4 (warunkowa) — doprecyzowania FROZEN

- Superkomórka 2×2×2 (LOCK), siatka N_sc = 2·48 = 96 na wymiar
  (tło N=48 powielone ×2 w każdym kierunku).
- **Mod perturbacji:** minimalny NIE-translacyjny mod z Phase 3 (N=48)
  w argmin k po 4 zalockowanych punktach (wszystkie są współmierne
  z superkomórką 2×2×2: e^{ik_α d}=±1); rozszerzenie na superkomórkę
  φ_sc(x+m·d) = φ(x)·Π_α(±1)^{m_α} (rzeczywiste); normalizacja
  ‖φ_sc‖∞ = 1. Perturbacja wzdłuż czystej translacji jest trywialna
  (przesuwa tło) — wykluczenie translacji jak w §6.
- Amplituda: a = 0.01·‖tło‖, ‖tło‖ := max|g₀−1| (LOCK; konwencja #63).
- Ewolucja: kanoniczna K_ε, RK4, dt=0.01 (CFL: najgęstsza siatka
  h=π/48=0.0654 ⟹ limit RK4 ≈ 0.075), kontrole: dt/2 przy (+, ε=0.2)
  oraz ε=0.1 przy (+); oba znaki ± przy ε=0.2; t_end = 3·t*_izo(3D z P1c).
- **Ucieczka :=** (g ≤ 0.01 lub niefinityczne — co krok) LUB
  max|g−g₀| > 1.0·‖tło‖ (100% tła — co próbkę 0.02). Raport t_esc vs
  2·t*_izo(3D) i 3·t*_izo(3D); gate energii ≤ 1e−6 do momentu ucieczki.
- Phase 4 dla KAŻDEGO tła ze zbieżnym Phase 3 (warunek = zbieżność,
  nie znak). Ruling kwantyfikatora d przy istnieniu częściowym: zapisany
  w `Phase3_verdict_ruling.md` PO Phase 3, PRZED Phase 4 (wzorzec
  bloch-chain: PRIMARY po tłach istniejących, strict obok).

## 8. Higiena wykonania

- Python: `python` (CPython 3.14.2; numpy 2.4.3, scipy 1.17.1 —
  sprawdzone importem PRZED startem; sympy niepotrzebny).
- Zero zmiany cwd; ścieżki pełne; po każdym zapisie weryfikacja `ls`
  (znany artefakt zagnieżdżonych katalogów). Outputy: `Phase*_output.txt`.
- Tła Phase 2 → `Phase2_backgrounds3d.npz`; wyniki Phase 3 →
  `Phase3_results3d.json`.
- Runy >10 min uruchamiane w tle; zalockowane siatki NIE zmniejszane;
  pojedynczy punkt wolno raportować INCOMPLETE z przyczyną.
- Rdzeń `.tex`, STATE.md, git — NIETYKANE. Katalogi innych cykli —
  NIETYKANE (odczyt wzorców dozwolony).

**FROZEN 2026-08-31, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
