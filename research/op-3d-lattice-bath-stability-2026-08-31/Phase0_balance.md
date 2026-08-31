---
title: "Phase0_balance — LOCK: stabilność w kąpieli w 3D — sieć periodyczna solitonów, samouzgodnione tło + Bloch (ostatnia droga rachunkowa hipotezy stabilizacji gęstością)"
date: 2026-08-31
type: phase0-lock
tgp_owner: research/op-3d-lattice-bath-stability-2026-08-31
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-08-31: 'ok 1 Cykl 3D dla stabilności w kąpieli' (po Q-FAIL 1D w trzech klasach dynamiki: bloch-chain II rzędu, gradient flow, symplectic-Jspectrum)"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-symplectic-Jspectrum-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-nonlinear-charge-constraint-2026-07-03/Phase3_nonlinear_evolution.py]]"
---

# Phase 0 — LOCK cyklu `op-3d-lattice-bath-stability`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu. Kryteria,
punkty skanu, kontrole, progi i forbidden moves zamknięte PRZED kodem.**

---

## 0. Pytanie binarne (Q) i dlaczego 3D jest inne niż 1D

**Czy istnieje separacja d sieci prostej kubicznej solitonów w sektorze
tachionowym, przy której najniższa gałąź pasmowa fluktuacji wokół
SAMOUZGODNIONEGO tła 3D jest dodatnia: ω²_min(d) = min_k λ_min(L̂_d(k)) > 0
(po odjęciu modów translacyjnych/Goldstone)?**

Dlaczego negatyw 1D się nie przenosi (pre-rejestrowane, żeby werdykt nie
był interpretowany poza zakresem):
1. Argument węzłowy Hilla (1D: stan podstawowy pod węzłowym modem
   translacyjnym) NIE obowiązuje w 3D (mody translacyjne mają strukturę
   wektorową; twierdzenie węzłowe nie porządkuje widma tak samo).
2. Ogony 3D zanikają ~cos(r−δ)/r (izolowany soliton istnieje; #63),
   w 1D nie zanikały wcale — tło 3D jest bliżej próżni między węzłami,
   ale sprzężenie sieci jest słabsze; kierunek efektu netto = właśnie
   pytanie Q.
3. Wynik 1D (kąpiel SPOWALNIA wzrost: bloch −1.22 vs próżnia −1;
   symplectic 0.14 vs 0.25) dopuszcza monotoniczny trend, który w 3D
   może przekroczyć zero — albo nie. Oba wyniki pełnoprawne.

To jest RACHUNEK CENTRALNY programu (ω²(n) z retrospektywy 2026-08-16),
w ostatniej niewykluczonych geometrii.

## 1. Kontekst dziedziczony (NIE powtarzać)

- **Baseline izolowany #63 V3** (reprodukowany w bath-two-sectors P2a):
  soliton μ, g₀=2.02117, R=60, model M0-f_ε; λ_min(w1)=−1.3896, t*=3.62.
  Dokładna forma modelu: `../op-nonlinear-charge-constraint-2026-07-03/
  Phase3_nonlinear_evolution.py` + `../op-bath-two-sectors-2026-08-23/
  Phase2_bath_runaway.py` — implementator MUSI odtworzyć i zamrozić
  w `Phase_method_decisions.md` PRZED startem.
- **Metoda zwalidowana:** samouzgodnione tło + periodyczne BC + Bloch
  (bloch-chain: pełna zbieżność w 1D; przenieść na 3D kartezjańskie).
- **Zmierzona drabina:** d*₁(μμ, M-P) = 3.0790 [INPUT z bath-two-sectors
  Phase 1]; długość fali ogona 2π (ω_tail=1).

## 2. Model ZAMKNIĘTY

- **M-3D:** sektor tachionowy w akcji kanonicznej (K=g⁴, U=g⁷/7−g⁸/8;
  statyka radialna = maszyneria 2), siatka kartezjańska 3D, komórka
  [0,d)³ z periodycznymi BC, jeden soliton na komórkę (sieć prosta
  kubiczna). Regularyzacja f_ε (ε=0.2, kontrola 0.1 [INPUT #63])
  wyłącznie w ewolucji nieliniowej; spektra bez regularyzacji
  (uzasadnienie jak w bloch-chain method_decisions — przyjęte).
- **Skan d (zalockowany, 5 punktów):** d ∈ {π, 2π, 3.0790 [INPUT d*₁μμ],
  3π, 4π}; wolno dodać, NIE usuwać.
- **Starty relaksacji (zalockowane 2):** superpozycja profili izolowanych
  (PRIMARY; profil radialny z Phase 1) + wariant przeskalowany 0.7×
  amplitudy; wolno dodać, NIE usuwać.
- **Siatki (zalockowane):** N ∈ {32, 48} punktów na wymiar komórki
  (zbieżność = porównanie N; wolno dodać N=64, NIE zastępować pary).
- **Zbiór k (zalockowany):** punkty wysokiej symetrii strefy Brillouina
  sc: Γ=(0,0,0), X=(π/d,0,0), M=(π/d,π/d,0), R=(π/d,π/d,π/d);
  wolno dodać punkty, NIE usuwać.
- **Solver spektralny:** operator rzadki, scipy.sparse.linalg.eigsh
  (najniższe ≥10 wartości własnych, tryb 'SA' lub shift-invert —
  wybór FROZEN w method_decisions).
- **Rejestr WEJŚĆ:** g₀=2.02117 (kalibracja μ #63), ε, d*₁=3.0790,
  β=γ=1 — flagować INPUT w każdym zależnym wyniku.

## 3. Fazy i kryteria (zalockowane)

### Phase 1 — bramka maszynerii 3D (osiągalne FAIL-e; TANIA — przed ciężkim rachunkiem)
- **P1a (dyspersja próżni 3D, exact):** Bloch na g≡1: sektor tachionowy
  MUSI dać ω²(k)=|k|²−γ, kontrola stabilna |k|²+γ; gate: błąd wzgl.
  ≤ 1e−3 na zalockowanym zbiorze k przy N=32 ORAZ redukcja ~×4 przy
  N=64 (rząd 2; dla bramki dopuszczone małe pudło L=2π zamiast d).
- **P1b (kotwica radialna):** kod radialny 1D (tani): reprodukcja
  λ_min(w1) = −1.3896 ± 1e−3 (kotwica #63/P2a) oraz t* = 3.62 ± 0.05
  (dt ×2). FAIL ⟹ STOP.
- **P1c (most radialny→kartezjański 3D, KLUCZOWA bramka):** izolowany
  soliton μ interpolowany na siatkę 3D w dużym pudle (L ≥ 30, zero-padding
  do próżni): λ_min(3D) MUSI zgodzić się z λ_min radialnym w granicach
  ±5% przy N odpowiadającym h≈0.4 i poprawiać się z zagęszczeniem;
  t*_izo(3D) zmierzone (do kryteriów Phase 4) i zgodne z 3.62 ± 15%.
  FAIL ⟹ STOP (dyskretyzacja 3D nieadekwatna — raport bez Phase 2–4).
- **P1d (gate energii):** ewolucja 3D na próżni+perturbacja:
  |ΔE|/E ≤ 1e−6, zbieżność dt ×2.

### Phase 2 — istnienie tła sieci 3D
- Relaksacja (Newton-Krylov lub tłumiony gradient flow — FROZEN
  w method_decisions) do ‖residuum EL‖∞ ≤ 1e−8 (3D, luźniejszy niż 1D
  z powodu kosztu — zalockowane TERAZ), dla 5 d × 2 starty × N ∈ {32,48}.
- Kryteria: niestałość ‖g−1‖∞ ≥ 0.05; zbieżność siatkowa
  ‖g₃₂−g₄₈‖∞ ≤ 5e−3 (interpolacja na wspólną siatkę).
- Brak istnienia dla WSZYSTKICH d ⟹ **CLOSED-GATE-STOP** (pełnoprawny
  negatyw: sieć 3D nie istnieje w klasie zbadanej ⟹ hipoteza drabiny
  wymaga innej geometrii — raport wprost).

### Phase 3 — RACHUNEK CENTRALNY
- Dla każdego istniejącego tła: λ_min(L̂_d(k)) na zalockowanym zbiorze k,
  N ∈ {32,48}; **identyfikacja modów translacyjnych (3 sztuki,
  λ ~ O(h²)) i ewentualnych modów cechowania PRZED interpretacją —
  kontrola nieusuwalna** (lekcja Goldstone 1D); ω²_min = najniższa
  wartość własna PO odjęciu zidentyfikowanych modów zerowych.
- Zbieżność: |Δω²_min(N₃₂→N₄₈)| ≤ 0.05·max(|ω²_min|, 0.1)
  (próg luźniejszy niż 1D — koszt 3D; zalockowany TERAZ).
- **P3b (kontrola artefaktu, nieusuwalna):** spektrum próżni g≡1 w tej
  samej komórce/superkomórce: gałęzie |k+G|²−γ odtworzone ≤ 1e−2.
- **P3c (kontrola sektora stabilnego, nieusuwalna):** ta sama procedura
  (relaksacja + Bloch) w sektorze m²=+γ ze źródłem gaussowskim 3D:
  ω²_min > 0 wymagane (osiągalny FAIL maszynerii).

### Phase 4 (warunkowa — tylko przy zbieżnym Phase 3) — test nieliniowy
- Superkomórka 2×2×2 (8 komórek; koszt!), perturbacja ± wzdłuż modu
  minimalnego (0.01·‖tło‖), ewolucja M0-f_ε do t = 3·t*_izo(3D z P1c);
  gate energii ≤ 1e−6.
- **Q-PASS:** istnieje d z ω²_min(d) > 0 (zbieżnie) ORAZ brak ucieczki
  (g→0 lub ‖δg‖>100% tła) do 3·t*_izo przy obu znakach perturbacji.
- **Q-FAIL:** ω²_min(d) < 0 dla wszystkich ISTNIEJĄCYCH teł (zbieżnie)
  ORAZ ucieczka ≤ 2·t*_izo. (Ruling kwantyfikatora przy istnieniu
  częściowym: jak w bloch-chain — po tłach istniejących, odczyt strict
  raportowany obok; ruling do zapisania PRZED Phase 4.)
- **Q-INCONCLUSIVE:** pozostałe (w tym niezbieżność) — wprost.
- Deskryptywnie (bez progu, do close): trend ω²_min(d) vs 1D
  (czy 3D podnosi względem −1.22) — nawet przy Q-FAIL wartościowa
  charakterystyka kierunku.

## 4. Forbidden moves
1. Zero zmian kryteriów/progów/punktów d/k/N/startów po starcie obliczeń;
   korekty tylko dla udokumentowanego błędu implementacji
   (correction_note PRZED użyciem wyniku; pierwotne outputy zachowane).
2. Kontrole P1a/P1b/P1c/P3b/P3c nieusuwalne; każdy test z osiągalnym FAIL.
3. Zakaz wyboru „lepszej" siatki/startu/punktu k post-hoc; wszystko
   raportowane.
4. Wyniki negatywne wprost; niezbieżność jako niezbieżność; zakaz
   skanowania do celu.
5. Rdzeń `.tex` NIETKNIĘTY; STATE.md nieedytowany; git nieużywany
   (commit robi sesja główna); katalogi innych cykli nietknięte.
6. Decyzje implementacyjne FROZEN w `Phase_method_decisions.md` PRZED
   startem obliczeń danej fazy; budżet czasowy: pojedynczy run > 10 min
   → uruchamiać w tle i raportować częściowo, NIE zmniejszać zalockowanych
   siatek (wolno raportować INCOMPLETE punktu).
7. Higiena ścieżek: bez `cd`; weryfikacja materializacji po każdym zapisie.

## 5. Deliverables
`Phase_method_decisions.md`, `Phase1_gate3d.py`, `Phase2_relax_lattice3d.py`
(+ `Phase2_backgrounds3d.npz`), `Phase3_bloch3d.py` (+ json),
warunkowo `Phase4_nonlinear3d.py`, outputy `Phase*_output.txt`,
`Phase_FINAL_close.md`, `NEEDS.md` (user-gated), `README.md` (log faz).

## 6. Drzewo decyzyjne
```text
P1b/P1c FAIL → STOP (maszyneria; raport, bez Phase 2–4)
P2 brak istnienia wszędzie → CLOSED-GATE-STOP (negatyw istnienia sieci sc;
   NEEDS: inne sieci (bcc/fcc) / inna geometria — user-gate)
Q-PASS → NEEDS: „stabilność = własność konfiguracji 3D o skończonej
   gęstości" (user-gate: dopisek retrospektywa + sek08b) + charakterystyka
   ω²(d) + PILNE: konsekwencje dla drabiny mas
Q-FAIL → NEEDS: Limitations — hipoteza stabilizacji gęstością zamknięta
   także w 3D sc (ostatnia geometria klasy); decyzja aksjomatyczna
   o znaku W zostaje bez podpory rachunkowej
Q-INCONCLUSIVE → NEEDS: wniosek metodologiczny (co ograniczyło zbieżność
   3D), bez ontologii
```

---

**LOCK ZAMKNIĘTY 2026-08-31. Zmiany poniżej tej linii po starcie obliczeń
= forbidden move (poza datami realizacji faz).**
