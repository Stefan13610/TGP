---
title: "Phase_FINAL_close — zamknięcie: Q1-INCONCLUSIVE (niezbieżność siatka×komórka we wszystkich 15 punktach; komparator pokazuje breakdown także BEZ kąpieli), Q2-FAIL (gęstość NIE zmienia znaku w akcji stabilnej — ω²_min>0 i ROŚNIE z n; kontrola d=∞ czysta)"
date: 2026-08-29
type: phase-final-close
tgp_owner: research/op-bath-two-sectors-2026-08-23
status: CLOSED
verdict: "GATE Phase 1: PASS wg rulingu zakresu (Phase1_gate_ruling.md; strict-reading z opcjonalnym τ: FAIL — oba odczyty raportowane). Zmierzone (A,δ) PRIMARY [120,260]: M-P e (0.093576, −75.34°), M-P μ (0.615899, +97.20°), M-L e (0.095922, −81.43°), M-L μ (0.363732, +38.58°) — Δ_ML(e→μ)=120.01° reprodukuje N1 co do setnej; drabiny 2π±5% PASS 6/6 (odchył 0.96–3.34%); P1c czysta 6/6; κ_eff≤5.4e−6→κ=0. Q1: INCONCLUSIVE — P2a baseline PASS (λ_min(w1)=−1.3896 vs #63 −1.389; t*=3.62 stabilne w dt; gate ≤1.64e−8), ale ω²_min niezbieżne we WSZYSTKICH 15 punktach (rozrzut 0.28–274.9 przy progu ≤0.01·max(|ω²|,0.1)); znak ω² podąża za rozmiarem komórki (R<π: +, R>π: −), a komparator izolacja-w-komórce (bez kąpieli) daje niemal identyczne ω² i t* — konstrukcja komórkowa nie separuje efektu kąpieli. Q2: FAIL — kontrola d=∞ PASS (m²=0.99997/0.99999, gate γ±5%; ω²_min=+1.00687≈+γ), ω²_min(d)= +1.34464 (d=8), +1.56855 (d=6), +1.88310 (d=4), +2.46974 (d=2) — wszystkie ZBIEŻNE (|dv|≤5.7e−6), wszystkie DODATNIE i ROSNĄCE z gęstością; odpowiedź statyczna monotoniczna/Yukawa (0 zmian znaku) dla wszystkich d; wrażliwość q=0.3 (d=4): +1.28824. Hipoteza dwóch sektorów OBALONA w klasie zbadanej: rozdwojenie znaku W pozostaje problemem aksjomatycznym."
anti_lakatos_lock: PRESERVED
tags: [bath-stability, two-sectors-hypothesis, measured-phases, cell-spectrum, negative-result, inconclusive-q1, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[README.md]]"
  - "[[Phase1_gate_ruling.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase2_correction_note.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# Phase FINAL — zamknięcie cyklu op-bath-two-sectors

Obliczenia wykonane 2026-08-23 (sesja implementatora, przerwana przed
dokumentami zamykającymi); zamknięcie dokumentacyjne 2026-08-29.
Kryteria LOCKa (Phase0_balance.md §3, §6) stosowane DOSŁOWNIE; zero zmian
progów/punktów/okien po starcie (korekty implementacyjne: §4 poniżej).

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| Gate Phase 1 | **PASS** (ruling zakresu; strict-reading: FAIL — oba odczyty) | wszystkie wielkości zasilające Phase 2 przeszły; jedyny FAIL to opcjonalny τ w oknie kontrolnym [50,150] |
| **Q1** (runaway w kąpieli) | **Q1-INCONCLUSIVE** | ω²_min niezbieżne siatka×komórka we wszystkich 15 punktach; LOCK §3: „niezbieżność raportować JAKO niezbieżność" |
| **Q2** (dwa sektory jednej akcji) | **Q2-FAIL** | gęstość źródeł w akcji STABILNEJ nie generuje znaku tachionowego — ω²_min(d)>0 i ROŚNIE z gęstością; hipoteza autora obalona w klasie zbadanej |

Oba wyniki są pełnoprawnymi rozstrzygnięciami (LOCK §0). Q2-FAIL jest
czystszy poznawczo niż Q1-INCONCLUSIVE: procedura Q2 miała czystą kontrolę
d=∞ i pełną zbieżność — wynik negatywny jest wiarygodny.

## 1. Phase 1 — zmierzone (A, δ), drabina, gate (Phase1_output.txt)

Okno PRIMARY [120,260], R=450; R-kontrola (300 vs 450) idealna
(|Δδ|=0.0000° przy gate ≤3.6°); κ_eff ≤ 5.4e−6 → κ:=0 wszędzie.

| Model | Gatunek | A | δ | R² [120,260] | R² [50,150] |
|---|---|---|---|---|---|
| M-P | e | 0.093576 | −75.34° | 0.9999978 | 0.9999906 |
| M-P | μ | 0.615899 | +97.20° | 0.9998993 | 0.9995755 |
| M-P | τ (INPUT Q_K=3/2) | 1.006011 | −171.94° | 0.9997182 | **0.9987829 FAIL** |
| M-L | e | 0.095922 | −81.43° | 0.9999994 | 0.9999974 |
| M-L | μ | 0.363732 | +38.58° | 0.9999908 | 0.9999615 |

- **Δ_ML(e→μ) = 120.01°** — dokumentacyjna różnica faz dodatku H
  reprodukowana co do setnej w układzie logarytmicznym (potwierdzenie N1
  na niezależnie zmierzonych fazach). W M-P: Δ=172.54° (inny układ, inne
  fazy — spójne z werdyktem A3 poprzedniego cyklu).
- **Drabiny d\*** (odstęp 2π±5%, predykcja zalockowana): MP ee 3.3644
  (odchył 2.99%), MP eμ 6.5124 (1.19%), MP μμ 3.0790 (3.34%), ML ee 3.1315
  (3.27%), ML eμ 5.3504 (1.59%), ML μμ 7.4972 (0.96%) — **PASS 6/6**,
  korekta malejąca z d zgodnie z predykcją.
- **P1c (Yukawa bez cos): 0 minimów w 6/6 par×modeli — CZYSTA.**
- **τ:** FAIL deskryptywny w [50,150] (flaga TAU-NEAR-THRESHOLD; 1.9% pod
  progiem 8/5); wrażliwość +2% → **KOLAPS** (r_end=5.9) — zgodnie z progiem
  H7; −2% → PASS. Żadna liczba τ nie weszła do Phase 2.
- **Gate:** ruling zakresu zapisany PRZED Phase 2 ([[Phase1_gate_ruling.md]]):
  gate obejmuje wielkości zasilające Phase 2 → **PASS**; strict-reading
  (z opcjonalnym τ, oknem kontrolnym) → FAIL. Oba odczyty raportowane;
  progi nietknięte.

## 2. Q1 — runaway w kąpieli (Phase2_output.txt + komparator)

### P2a baseline (gate kodu)

Reprodukcja #63 V3 (soliton μ, g₀=2.02117, R=60, N=4000):
λ_min(w1)=−1.3896 (#63: −1.389 ✓), λ_min(F)=−7.8579; BREAKDOWN t\*=3.62
identyczny dla dt=0.004 i 0.002 oraz a=±0.01; gate |ΔE|/E ≤ 1.64e−8
(próg 1e−6) — **P2a PASS**. Komparator e (dodatek): t\*_izo(e)=7.08;
gate e FAIL przy a=+0.01 (artefakt overflow po przejściu g→0, korekta 1
w [[Phase2_correction_note.md]] — poza zakresem gate'u P2a).

### P2b/P2c skan (3 konfiguracje × 5 zalockowanych d)

Wynik strukturalny, identyczny we wszystkich 15 punktach:

1. **Niezbieżność siatka×komórka:** rozrzut ω²_clean 0.28–274.9 przy
   progu ≤0.01·max(|ω²|,0.1) — **0/15 punktów zbieżnych**. Znak ω²_min
   podąża za ROZMIAREM KOMÓRKI, nie za gęstością: komórki R<π dają ω²>0
   (brak modów komórkowych; np. ee d=1.682, R=0.84: +23.9), komórki R>π
   dają ω²<0 (np. μμ d=4.618, R=4.62: −272.7). P2c identyfikuje mody
   próżni (N_neg=floor(R/π) potwierdzone we wszystkich komórkach), ale po
   odjęciu modów próżni (w2_clean) rozrzut pozostaje.
2. **Mod runaway przeżywa w dużych komórkach:** dla centralnego μ
   ω²_min ≈ −7.8…−9.5 ≈ λ_min(F) izolacji (−7.8579) — kąpiel przy
   zalockowanych d przesuwa tło o c_bath ∈ [−0.03, +0.45] (ladder points),
   za mało by odwrócić znak. Duży |c_bath| tylko w punktach 0.5·d\*₁
   (μμ: +4.74, ee: −0.66) — tam spektra najbardziej rozjechane.
3. **Ewolucje:** BREAKDOWN we wszystkich 90 przebiegach
   (t\*/t\*_izo = 0.03–1.81), w tym przy amp=0 (sloshing niestacjonarnego
   tła — kontrola dodana w method_decisions p.8). Warunek Q1-PASS
   („brak runaway do 3·t\*_izo") nie zaszedł nigdzie; warunek Q1-FAIL
   („ω²<0 wszędzie ORAZ ucieczka ≤2·t\*_izo") też nie (punkty z ω²>0
   w małych komórkach; t\* rozrzucone po obu stronach 2·t\*_izo... nie —
   wszystkie t\* ≤ 1.81·t\*_izo, ale ω² NIE jest <0 wszędzie).
4. **Komparator izolacja-w-komórce (bez kąpieli, amp=0):** ω² i t\*
   niemal identyczne jak z kąpielą (np. eμ d=12.87 R=12.87: izol.
   ω²=−7.858 vs kąpiel −7.79/−8.07; breakdowny w tych samych oknach).
   **Diagnoza: obcięty profil + zero-flux BC w komórce jest
   niestacjonarny i pada sam z siebie — konstrukcja nie separuje efektu
   kąpieli od artefaktu komórki.**

**Werdykt: Q1-INCONCLUSIVE** (LOCK §3, litera: niezbieżność siatek/komórek
raportowana jako niezbieżność; zakaz wyboru „lepszej" siatki post-hoc).
Pytanie N3 (czy kąpiel stabilizuje runaway) pozostaje OTWARTE — ale cykl
ustala deskryptywnie, że komórkowy wariant M-B z obciętym tłem jest
NIEZDOLNY rozstrzygnąć Q1 (wniosek metodologiczny → NEEDS).

## 3. Q2 — dwa sektory jednej akcji (Phase3_output.txt)

Wariant wymiarowy FROZEN przed startem: komórka 1D periodyczna, σ_s=0.5,
q=1.0 [INPUT] (Phase_method_decisions.md Phase 3).

- **P3a kontrola d=∞ (osiągalny FAIL procedury): PASS.**
  N=4000/8000: fit m=0.99998/1.00000 → m²=0.99997/0.99999 (gate γ±5%,
  odchył 0.00%); zanik monotoniczny; ω²_min=+1.00687 ≈ +γ ✓.
- **ω²_min(d), oba N zbieżne (|dv| ≤ 5.7e−6):**

  | d | ω²_min | odpowiedź statyczna |
  |---|---|---|
  | ∞ | +1.00687 | monotoniczna (0 zmian znaku) |
  | 8 | +1.34464 | monotoniczna |
  | 6 | +1.56855 | monotoniczna |
  | 4 | +1.88310 | monotoniczna |
  | 2 | +2.46974 | monotoniczna |

- Wrażliwość q=0.3 (d=4): ω²_min=+1.28824 — również dodatnie.
- Klasyfikacja odpowiedzi: superkomórka 10 okresów, L̂χ=−δ_h:
  0 zmian znaku poza 2σ_s dla WSZYSTKICH d; |χ| zanika.

**Werdykt: Q2-FAIL.** Gęstość źródeł w akcji stabilnej nie tylko NIE
generuje znaku tachionowego — **usztywnia** potencjał fluktuacji
(ω²_min rośnie monotonicznie z n; mechanizm: tło ψ_n rośnie z gęstością,
a wokół ψ>1 krzywizna 3γψ²−2βψ rośnie). Hipoteza „pojedynczy obiekt vs
obiekt w kąpieli = dwa sektory jednej akcji" jest OBALONA w klasie
zbadanej (1D periodyczne, źródła gaussowskie, q∈{0.3,1.0}). Zgodnie
z drzewem LOCKa §6: **rozdwojenie znaku W (N2) pozostaje otwartym
problemem AKSJOMATYCZNYM — decyzja ontologiczna autora, nie numeryka.**

## 4. Korekty implementacyjne i odstępstwa (pełna lista)

1. **Ruling zakresu gate'u Phase 1** ([[Phase1_gate_ruling.md]], zapisany
   przed Phase 2): LOCK nie domknął, czy gate agreguje opcjonalny τ i okno
   kontrolne. Oba odczyty raportowane (PASS pipeline / FAIL strict).
2. **Korekta 1** (detektor breakdown co krok RK4 zamiast co próbkę;
   [[Phase2_correction_note.md]]) — udokumentowana przed użyciem wyników;
   kryteria nietknięte.
3. **Korekta 2** (próg Newtona skalowany podłogą arytmetyki double;
   tamże) — j.w.
4. **Dodatki dozwolone:** komparator izolacja-w-komórce
   (`Phase2_isolation_cell_comparator.py`), przebiegi amp=0, komparator e
   w P2a. Punkty zalockowane nieusunięte.
5. **Gate energii w ewolucjach komórkowych:** masowo FAIL (overflow po
   breakdown) — raportowane wprost; gate zalockowany dotyczył P2a
   (reprodukcji #63) i tam PASS.
6. **Wzorzec frontmattera:** HANDOFF wskazywał
   `op-wall-dynamics-2026-07-03/Phase_FINAL_close.md`, który NIE ISTNIEJE
   (tamten cykl nie ma close-nota); użyto wzorca
   `op-lattice-bath-runaway-2026-08-23/Phase_FINAL_close.md`.
7. Sesja implementatora 2026-08-23 urwała się po Phase 3, przed
   dokumentami zamykającymi; niniejsze zamknięcie 2026-08-29 wyłącznie
   dokumentuje istniejące outputy (zero nowych obliczeń).

## 5. Pliki cyklu

Obliczeniowe: `Phase1_tail_measured.py` + `Phase1_output.txt`
+ `Phase1_results.json`; `Phase2_bath_runaway.py` + `Phase2_output.txt`;
`Phase2_isolation_cell_comparator.py` + `Phase2_output_comparator.txt`;
`Phase3_two_sectors.py` + `Phase3_output.txt`.
Metodologiczne: `Phase_method_decisions.md`, `Phase1_gate_ruling.md`,
`Phase2_correction_note.md`. Zamykające: ten plik, `NEEDS.md` (user-gated),
`README.md` (zaktualizowany).

## 6. Mapowanie na drzewo decyzyjne LOCKa §6

- Q1-INCONCLUSIVE: drzewo nie przewiduje gałęzi (tylko PASS/FAIL) —
  do NEEDS trafia wniosek METODOLOGICZNY (niezdolność wariantu
  komórkowego), nie ontologiczny; N3 pozostaje OTWARTE.
- Q2-FAIL → NEEDS: „znak W = otwarty problem aksjomatyczny (decyzja
  ontologiczna autora, nie numeryka)" — dosłownie wg drzewa.
