---
title: "Phase_FINAL_close — zamknięcie: Q-FAIL (odczyt PRIMARY wg rulingu; strict-literal: Q-INCONCLUSIVE — oba raportowane): samouzgodniony łańcuch periodyczny w sektorze tachionowym NIE stabilizuje — ω²_min(d) = −1.222…−1.230 < 0 zbieżnie dla wszystkich istniejących teł (d∈{3π,4π,6π}; dla d∈{π,2π} tło niestałe NIE istnieje), ucieczka nieliniowa t_esc=4.58–4.92 ≤ 2t*_ref przy obu znakach; kontrole P1a/P3b/P3c czyste"
date: 2026-08-31
type: phase-final-close
tgp_owner: research/op-bloch-chain-stability-2026-08-31
status: CLOSED
verdict: "GATE Phase 1: PASS (P1a tach/stab maxerr 8.3e−5 / 3.2e−5, ratio 4.000/4.000; P1b |ΔE|/E ≤ 6.9e−15, dt-konsystencja ≤ 1.3e−13). Phase 2 (istnienie): d=π i d=2π — KOLAPS-DO-PRÓŻNI z obu zalockowanych startów (S3-strzelanie potwierdza brak orbity o okresie ≤2π; pierwsza całka: T ≥ 2π z równością tylko przy zerowej amplitudzie); d=3π istnieje (g∈[0.5098,1.1406], oba starty → identyczne tło, ‖Δ‖=3.3e−16); d=4π istnieje (A=0.7; g∈[0.3312,1.1427]; A=0.3 kolaps); d=6π: start zalockowany → tło 2-garbne (powielona orbita 3π), dodany start S3single → tło jednogarbne (g∈[0.1675,1.1429]); wszystkie istniejące: ‖R‖∞ ≤ 5.7e−12 (próg 1e−10), ‖g_N−g_2N‖∞ ≤ 4.4e−5 (próg 1e−4). Phase 3 (centralna): ω²_min(d) ZBIEŻNE wszędzie (Δ_siatka ≤ 7.4e−5, Δ_k = 0 przy progu ~1.2e−2): 3π: −1.222191; 4π: −1.229340; 6π(2-garb): −1.222209; 6π(1-garb): −1.229829 — WSZYSTKIE UJEMNE, argmin k=0, mod minimalny NIE-translacyjny (Goldstone zidentyfikowany osobno: λ_trans = −2.4e−6…−1.3e−5 ~ O(h²)); cross-check u-formy zgodny do 8.8e−6; P3b (próżnia w superkomórce 4d, 5/5 d): PASS maxerr ≤ 2.7e−5; P3c (sektor stabilny, konstrukcja Q2): PASS ω²_min=+1.88310/+1.88311 (kotwica poprzednika +1.88310 odtworzona). Phase 4: ucieczka ‖δg‖>100% tła w t_esc = 4.58–4.92 ≤ 2t*_ref=7.24 we WSZYSTKICH 16 biegach (oba znaki, ε=0.2/0.1, dt/2 — t_esc stabilne do 0.02–0.14), σ_fit = 1.077–1.113 vs σ_lin=√|ω²_min|=1.105–1.109 (odchył ≤2.6%), gate energii ≤ 1.6e−7 PASS wszędzie. WERDYKT: Q-FAIL (PRIMARY, kwantyfikator po tłach istniejących wg Phase3_verdict_ruling.md); strict-literal (kwantyfikator po 5 zalockowanych d, z których 2 nie mają tła): Q-INCONCLUSIVE — oba odczyty raportowane, progi nietknięte."
anti_lakatos_lock: PRESERVED
tags: [bloch-chain, tachyonic-sector, band-structure, self-consistent-background, negative-result, goldstone-mode, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase3_verdict_ruling.md]]"
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# Phase FINAL — zamknięcie cyklu op-bloch-chain-stability

Obliczenia i zamknięcie 2026-08-31 (jedna sesja implementatora). Kryteria
LOCKa (Phase0_balance.md §0, §3, §6) stosowane DOSŁOWNIE; zero zmian
progów/punktów/siatek po starcie. Jedyna niedomkniętość literału LOCKa
(kwantyfikator „dla WSZYSTKICH zalockowanych d" przy istnieniu częściowym)
rozstrzygnięta rulingiem zapisanym PRZED Phase 4
([[Phase3_verdict_ruling.md]]) — oba odczyty raportowane.

---

## 0. Werdykt

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| Gate Phase 1 | **PASS** | dyspersja próżni exact w obu sektorach (rząd 2, ratio 4.000), gate energii 6.9e−15 |
| Istnienie (Phase 2) | **CZĘŚCIOWE** | łańcuch istnieje dla d∈{3π,4π,6π}; dla d∈{π,2π} NIE istnieje (kolaps do próżni; T_orbit ≥ 2π) |
| **Q** (ω²_min(d) > 0 dla pewnego d?) | **Q-FAIL** (PRIMARY; strict: Q-INCONCLUSIVE) | ω²_min ∈ [−1.230, −1.222] < 0 zbieżnie dla WSZYSTKICH istniejących teł + ucieczka nieliniowa ≤ 2t*_ref przy obu znakach |

**Sens fizyczny (deskryptywnie):** skończona gęstość struktur w 1D NIE
odwraca znaku najniższej gałęzi — najniższy mod łańcucha (amplitudowy,
k=0, nie-translacyjny) jest nawet NIECO głębiej ujemny niż tachion próżni
(−1.22…−1.23 vs −1.0) i pogłębia się z amplitudą tła (3π: −1.2222 →
6π jednogarbne: −1.2298). Periodyczna kąpiel samouzgodniona w tej klasie
(1D, sektor EL Noty kanonicznej) nie stabilizuje; ostatnia 1D-wersja
hipotezy stabilizacji gęstością zamknięta negatywnie — pozostaje 3D
i warianty poza-klasowe (drzewo LOCKa §6, gałąź Q-FAIL).

## 1. Phase 1 — bramka maszynerii (Phase1_output.txt): PASS

- **P1a (dyspersja próżni, exact):** komórka d=2π, metryka
  |Δω²|/max(|ω²_ex|,1) (podłoga zdefiniowana PRZED startem — gałąź
  przechodzi przez 0): sektor tachionowy maxerr = 8.327e−5 (N=400),
  ratio(N400/N800) = 4.000; sektor stabilny 3.203e−5, ratio 4.000 —
  gate ≤1e−3 + rząd 2: **PASS/PASS**. Zwinięcie gałęzi poprawne
  (k=0: −1, 0, 0; k=1/2: −3/4 podwójne).
- **P1b (gate energii):** RK4 na dokładnym semi-dyskretnym hamiltonianie
  akcji kanonicznej z K_ε; |ΔE|/E ≤ 6.93e−15 (ε=0.2 i 0.1; dt=0.004
  i 0.002), konsystencja dt ‖Δg‖∞ ≤ 1.3e−13 — **PASS**. Wzrost tachionowy
  perturbacji ~e^t obecny (kontrola sensowności).

## 2. Phase 2 — istnienie tła (Phase2_output.txt + addendum)

| d | start A=0.7 | start A=0.3 | tło do Phase 3 |
|---|---|---|---|
| π | KOLAPS-DO-PRÓŻNI | KOLAPS-DO-PRÓŻNI | — (S3: brak orbity) |
| 2π | KOLAPS (amp. resztkowa 8e−4→0 z N) | KOLAPS (j.w.) | — (S3: brak orbity) |
| 3π | NIESTAŁE-ZBIEŻNE | NIESTAŁE-ZBIEŻNE (identyczne, ‖Δ‖=3.3e−16) | g∈[0.509826, 1.140602] |
| 4π | NIESTAŁE-ZBIEŻNE | KOLAPS-DO-PRÓŻNI | g∈[0.331188, 1.142721] |
| 6π | NIESTAŁE-ZBIEŻNE, ale 2-GARBNE (powielona orbita 3π) | KOLAPS | 2-garbne + DODANY start S3single: 1-garbne g∈[0.167506, 1.142862] |

- Wszystkie istniejące tła: ‖residuum EL‖∞ ≤ 5.7e−12 na pełnej komórce
  (próg 1e−10), ‖g_400−g_800‖∞ ≤ 4.4e−5 na węzłach wspólnych (próg 1e−4),
  ‖g−1‖∞ = 0.49–0.83 (próg 0.05).
- **Nieistnienie dla d ≤ 2π jest strukturalne**, nie numeryczne: pierwsza
  całka E = ½g⁴g′² − U(g) daje okres orbit T(amplituda) ≥ 2π (równość
  w granicy zerowej amplitudy; T→∞ przy separatrysie g→0) — potwierdzone
  strzelaniem S3 (brak orbity o okresie ≤ 2π). Amplituda tła rośnie z d
  (raportowana, nie strojona — LOCK §2 Rejestr WEJŚĆ).
- **Incydent 1 (rozwiązany dodatkiem, nie korektą):** dla d=6π zalockowany
  start głęboki zbiegł do tła 2-garbnego (2 obiekty na komórkę wbrew
  specyfikacji M-C „jeden obiekt"). Zgodnie z LOCK §2 („wolno dodać starty,
  NIE usuwać") dodano start S3single (strzelanie na półokres 3π):
  tło jednogarbne, zbieżne. Wynik startu zalockowanego zachowany
  i policzony w Phase 3/4 równolegle (zakaz wyboru post-hoc).

## 3. Phase 3 — spektrum pasmowe (Phase3_output.txt): RACHUNEK CENTRALNY

**Tabela ω²_min(d) (wartość główna N=800, k=32; wszystkie ZBIEŻNE:**
Δ_siatka ≤ 7.4e−5, Δ_k = 0.0 przy progu 0.01·max(|ω²_min|,0.1) ≈ 1.2e−2):

| tło | d | ω²_min | argmin k | status |
|---|---|---|---|---|
| 3π (1-garbne) | 9.4248 | **−1.222191** | 0 | ZBIEŻNE |
| 4π (1-garbne) | 12.5664 | **−1.229340** | 0 | ZBIEŻNE |
| 6π (2-garbne, start zalockowany) | 18.8496 | **−1.222209** | 0 | ZBIEŻNE |
| 6π (1-garbne, S3single) | 18.8496 | **−1.229829** | 0 | ZBIEŻNE |

- **Mod minimalny NIE jest modem translacyjnym:** overlap B-ważony
  z g₀′ = 0.000; Goldstone zidentyfikowany osobno (reguła FROZEN
  w method_decisions §2): λ_trans = −2.4e−6…−1.3e−5 ~ O(h²) ≈ 0 ✓.
  Struktura k=0 (np. 3π): {−1.2222 (amplitudowy), −0.0000 (translacja),
  +0.632, +1.615}.
- **Cross-check kanoniczny (u=g³/3, waga 1, potencjał 4g−5g²):** zgodność
  z formą ważoną K=g⁴ do |Δ| ≤ 8.8e−6 przy k∈{0, π/d} — operator
  poprawny niezależnie od reprezentacji.
- **P3b (nieusuwalna):** próżnia w superkomórce 4d dla WSZYSTKICH 5
  zalockowanych d: 8 gałęzi odtwarza k²−1 zwinięte, maxerr ≤ 2.742e−5
  (gate 1e−3) — **PASS 5/5**.
- **P3c (nieusuwalna):** pełna procedura w sektorze stabilnym
  (konstrukcja Q2 poprzednika, d_s=4, σ_s=0.5, q=1): ω²_min = +1.88310
  (N=400) / +1.88311 (N=800) > 0, argmin k=0; kotwica poprzednika
  +1.88310 odtworzona do 1e−5 — **PASS** (osiągalny FAIL maszynerii
  nie zaszedł).
- Deskryptywnie: ω²_min prawie nie zależy od d (rozstęp 0.6%), jest
  BLISKIE tachionowi próżni (−1) i lekko GŁĘBSZE; pogłębia się
  z amplitudą tła (kolejność: 3π −1.2222 < 4π −1.2293 < 6π-1garb
  −1.2298). Kierunek PRZECIWNY do hipotezy stabilizacji gęstością.

## 4. Phase 4 — test nieliniowy (Phase4_output.txt)

Superkomórka 4d, perturbacja ±0.01·‖tło‖ wzdłuż modu minimalnego (k=0,
nie-translacyjny), K_ε (ε=0.2; kontrole ε=0.1 i dt/2), t*_ref=3.62:

| tło | t_esc (znak +) | t_esc (znak −) | σ_fit vs σ_lin |
|---|---|---|---|
| 3π | 4.64 | 4.92 | 1.1056 vs 1.1055 (+0.0%) / −2.6% (znak −) |
| 4π | 4.58 | 4.84 | 1.1126 vs 1.1088 (+0.3%) / −2.5% |
| 6π 2-garb | 4.64 | 4.90 | 1.1053 vs 1.1055 (−0.0%) / −2.5% |
| 6π 1-garb | 4.60 | 4.78 | 1.1108 vs 1.1090 (+0.2%) / −2.1% |

- **UCIECZKA (‖δg‖∞ > 100% tła) we wszystkich 16 biegach w t = 4.58–4.92
  ≤ 2·t*_ref = 7.24**; kontrole: ε=0.1 przesuwa t_esc o ≤0.14; dt/2 nie
  zmienia t_esc (≤0.01); gate energii ≤ 1.6e−7 **PASS wszędzie** (do
  momentu ucieczki; bez overflow — detekcja przed breakdownem g≤0).
- Wzrost wykładniczy zgodny z przewidywaniem liniowym z Phase 3
  (σ_fit/σ_lin ∈ [0.974, 1.003]) — spektrum Blocha i dynamika nieliniowa
  wzajemnie się potwierdzają.

## 5. Złożenie werdyktu (litera LOCKa §3 Phase 4 + ruling)

- **Q-PASS:** wymaga d z ω²_min(d) > 0 (zbieżnie) — NIE ISTNIEJE
  (kwantyfikator egzystencjalny; rozstrzygalne identycznie w obu
  odczytach rulingu). Nie zaszło.
- **Q-FAIL:** „ω²_min(d) < 0 dla WSZYSTKICH zalockowanych d (zbieżnie)
  ORAZ ucieczka w t ≤ 2·t*_ref":
  - odczyt PRIMARY (po tłach istniejących, [[Phase3_verdict_ruling.md]]):
    4/4 tła ujemne zbieżnie + ucieczka 16/16 ≤ 2t* ⟹ **Q-FAIL**.
  - odczyt strict-literal (po 5 zalockowanych d): dla {π, 2π} ω²_min
    niezdefiniowane (tło nie istnieje) ⟹ warunek nieewaluowalny ⟹
    formalnie **Q-INCONCLUSIVE** („pozostałe przypadki").
- **Werdykt nagłówkowy: Q-FAIL** (PRIMARY), z jawną flagą odczytu strict.
  Wynik negatywny pełnoprawny (LOCK §0): „kąpiel periodyczna nie
  stabilizuje w klasie zbadanej — ostatnia wersja hipotezy stabilizacji
  gęstością w tym sektorze zamknięta negatywnie".

**Obserwacja analityczna post-hoc (nie kryterium; korroboruje wynik):**
w formie kanonicznej problem fluktuacji to równanie Hilla
−χ″ + Ū″(u₀(x))χ = ω²χ. Mod translacyjny χ = u₀′ ma węzły (2 na komórkę),
a stan podstawowy Hilla jest bezwęzłowy — stąd λ₀ < λ(u₀′) = 0 ŚCIŚLE dla
KAŻDEGO niestałego tła periodycznego. Ujemność ω²_min w 1D skalarnej
klasie M-C jest więc strukturalna (twierdzenie oscylacyjne), a numeryka
kwantyfikuje głębokość (≈ −1.22, mod amplitudowy). To wzmacnia gałąź
NEEDS „pozostaje 3D/poza-klasowe": w 3D argument węzłowy w tej formie
nie obowiązuje.

## 6. Korekty, dodatki, incydenty (pełna lista)

1. **Rozjazd form EL z #63** udokumentowany PRZED startem
   (Phase_method_decisions.md §1): #63 M0-f_ε (F=f_ε, W waga-1) vs akcja
   kanoniczna (K=g⁴, U′=K·g²(1−g)); identyczna linearyzacja próżni;
   cykl w całości w akcji kanonicznej (spójność tło↔dynamika);
   regularyzacja f_ε (mapa ½(x+√(x²+ε²)), ε=0.2/0.1) zastosowana do K
   w ewolucji nieliniowej, spektra bez regularyzacji (uzasadnienie tamże).
2. **Dodany start S3single dla d=6π** (LOCK §2 dozwala dodawać starty):
   `Phase2_s3_6pi_addendum.py` + `Phase2_output_addendum.txt`; wyniki
   startów zalockowanych zachowane w `Phase2_output.txt` i w npz;
   tło 2-garbne policzone równolegle w Phase 3/4.
3. **Ruling kwantyfikatora d** ([[Phase3_verdict_ruling.md]]) zapisany
   po Phase 3, PRZED Phase 4; oba odczyty raportowane (wzorzec:
   Phase1_gate_ruling poprzednika).
4. Zero plików `Phase*_correction_note.md` — żadna kontrola nie wykryła
   błędu implementacji wymagającego korekty; wszystkie bramki przeszły
   w pierwszym biegu (P1a/P1b/P3b/P3c), pierwotne outputy nienaruszone.
5. Incydent środowiskowy bez wpływu na wyniki: sandbox blokuje cwd poza
   vaultem — wszystkie wywołania z pełnymi ścieżkami (zgodnie z higieną
   ścieżek LOCKa §4 p.8 i tak obowiązującą).

## 7. Pliki cyklu

Obliczeniowe: `Phase1_machinery_gate.py` + `Phase1_output.txt`;
`Phase2_relax_chain.py` + `Phase2_output.txt`; `Phase2_s3_6pi_addendum.py`
+ `Phase2_output_addendum.txt`; `Phase3_bloch_spectrum.py`
+ `Phase3_output.txt` + `Phase3_results.json`; `Phase4_nonlinear.py`
+ `Phase4_output.txt`; artefakt danych: `Phase2_backgrounds.npz`.
Metodologiczne: `Phase_method_decisions.md`, `Phase3_verdict_ruling.md`.
Zamykające: ten plik, `NEEDS.md` (user-gated), `README.md`
(zaktualizowany). Rdzeń `.tex`, STATE.md, git — nietknięte.

## 8. Mapowanie na drzewo decyzyjne LOCKa §6

Q-FAIL → NEEDS: „Limitations — periodyczna kąpiel samouzgodniona nie
stabilizuje w 1D (klasa zbadana); pozostaje 3D/poza-klasowe" — dosłownie
wg drzewa; plus (z Phase 2) negatyw istnienia dla d ≤ 2π i (deskryptywnie)
argument węzłowy Hilla wskazujący, że 1D-klasa była strukturalnie
przesądzona — co podnosi priorytet wariantu 3D jako JEDYNEGO otwartego
w tej linii.
