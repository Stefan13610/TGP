---
title: "NEEDS — user-gated po Q-I1-INCONCLUSIVE + Q-I2-FAIL (drzewo LOCK §5, gałąź „INCONCLUSIVE → NEEDS metodologiczny"): stan podprogowy indukowany materią ISTNIEJE na siatce produkcyjnej (ψ̄(0)=0.773/0.808 < 5/6), ale jest ZALEŻNY OD SIATKI tuż pod λ̃_crit=0.1077 i — co rozstrzygające — NIE PRZEŻYWA wygaszenia źródła (4/4 RETURN-TO-VACUUM): materia statyczna deformuje albo niszczy, nie kreuje"
date: 2026-09-15
type: needs
tgp_owner: research/op-matter-induced-creation-2026-09-15
status: USER-GATED
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/NEEDS.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# NEEDS (wszystko user-gated; ZERO samowolnych dopisków do core)

Werdykty cyklu: **Q-I1-INCONCLUSIVE** + **Q-I2-FAIL** (litera LOCK §4).
Drzewo LOCK §5 aktywuje gałąź ostatnią: „INCONCLUSIVE → NEEDS
metodologiczny (t_on, Δ rampy, siatka λ̃)". Gałęzie „czyste"
(Q-I1-PASS ∧ Q-I2-FAIL oraz Q-I1-FAIL ∧ Q-I2-FAIL) nie zostały
aktywowane, bo litera Q-I1 nie jest rozstrzygnięta — ich skutki
(domknięcie gałęzi „kreacja z materii statycznej") wymagają decyzji
autora, patrz N2.

## N1 (z drzewa, METODOLOGICZNY — główny) — zbieżność siatkowa przy λ̃_crit i reguła potwierdzenia

Fakt: na siatce produkcyjnej h=0.05 okno podprogowe jest niepuste
(λ̃=0.08: ψ̄(0)=0.808031; λ̃=0.10: ψ̄(0)=0.772903; próg M911 5/6),
a λ̃_crit=0.107656±0.000156 (bisekcja 6 kroków, FROZEN). LOCK §3
przewidział JEDNO potwierdzenie h=0.025 — dla **najgłębszego**
SETTLED-SUB, czyli dla biegu leżącego 7% pod progiem (λ̃=0.10);
ten bieg na h=0.025 **kolabuje** (t_end=0.95, pas górny), co wg
litery LOCK §4 blokuje zarówno PASS, jak i FAIL.

Probe deskryptywny (`Phase3_desc_gridconv_output.txt`, uruchomiony PO
zapisaniu werdyktu, bez wpływu na werdykt) pokazuje, że rozjazd jest
**lokalny wokół progu**, nie systemowy: λ̃=0.06 i λ̃=0.08 są zbieżne
siatkowo (Δψ̄(0) ≤ 2.7e−5, ta sama kategoria), a λ̃=0.10 przy dt/2
odtwarza SETTLED-SUB z ψ̄(0) zgodnym do 1.6e−10 (więc to efekt
SIATKI, nie kroku czasowego). Innymi słowy λ̃_crit(h) maleje przy
zagęszczaniu siatki (h=0.025 ⟹ λ̃_crit<0.10) — okno podprogowe się
zawęża, ale nie znika.

Kandydat re-locku (wyłącznie decyzją autora — nie wolno tego
„docenić" post hoc w obecnym cyklu):
1. reguła potwierdzenia siatkowego: **najgłębszy SETTLED-SUB ORAZ
   co najmniej jeden SETTLED-SUB odsunięty od progu** (np. λ̃ ≤
   0.8·λ̃_crit), z jawnym kryterium „która kategoria decyduje przy
   rozjeździe";
2. ekstrapolacja λ̃_crit(h) (h ∈ {0.05, 0.025, 0.0125}) — czy próg
   ma granicę h→0, czy zbiega do λ̃_turn 0D = 1/12 = 0.08333;
3. detektor progu odporny na pojedynczy overshoot (zdarzenie brzegowe
   zachodzi w jednym kroku — dziedziczona obserwacja poprzednika).
[USER-GATE, re-lock metodologiczny]

## N2 (interpretacyjny — domknięcie gałęzi M911-N2) — „kreacja z materii statycznej"

Q-I2-FAIL jest wynikiem mocnym i **niezależnym od losu Q-I1**: KAŻDY
stan osiadły, jaki statyczne źródło korpusowe zdołało wytworzyć —
łącznie z obydwoma stanami **podprogowymi** (ψ̄(0)=0.773 i 0.808,
detektor M911 przekroczony) — po adiabatycznym wygaszeniu źródła
wraca do próżni (4/4 RETURN-TO-VACUUM; kontrola negatywu h=0.025
zgodna; kontrola czystości maszynerii λ̃=0.01 PASS). Zero histerezy:
już przy t=700 max_{r≤40}|ψ−1| spada z 0.227 do 9.8e−6, a praca
wygaszania oddaje polu dokładnie energię studni (ΔE = +9.1488 przy
λ̃=0.10).

Pytanie do autora: czy to wystarcza, by **formalnie domknąć
negatywnie** gałąź „sprzężenie z materią" drzewa hipotezy kreacji
(M911-NEEDS N2), mimo że litera Q-I1 została INCONCLUSIVE? Warstwa
faktu domyka ją w całości (obiekt indukowany jest cieniem źródła);
warstwa litery zostawia otwarte tylko pytanie o zbieżność KATEGORII
stanu ze źródłem WŁĄCZONYM (N1), które dla kreacji jest wtórne.
[USER-GATE — domknięcie gałęzi drzewa nie jest decyzją
implementatora]

## N3 (metodologiczny, z wykonania) — kryterium RETURN-TO-VACUUM jest puste przy rampie adiabatycznej

Zamrożona litera RETURN-TO-VACUUM miała dwie klauzule:
E_core(1700) < 0.05·E_core(700) **lub** max_{r≤40}|ψ−1| < 1e−3
w oknie końcowym. Zadziałała wyłącznie klauzula amplitudowa, bo
przy Δ=100 rampa jest dla tego układu tak adiabatyczna, że
**E_ref^off = E_core(700) jest już wielkością szumową** (≈4e−4, pięć
rzędów poniżej |E_core(600)| ≈ 9.15) — iloraz E(1700)/E(700) ≈ 0.75
nie niesie informacji o „powrocie", tylko o powolnym rozpływie
resztki masywnego pola. Kandydat do wspólnej metodyki: w LOCKach
z rampą definiować progi energetyczne względem **stanu ze źródłem
włączonym** (E_core(t_off)) albo dobierać Δ tak, by E_ref^off było
mierzalne; oraz raportować parę (amplituda, energia) zawsze łącznie.
[USER-GATE przy kolejnym LOCKu]

## N4 (deskryptywny, kandydat nowego LOCKa) — próg jest DYNAMICZNY, nie statyczny

Pre-rejestrowana predykcja P1-I2 wskazywała fold 0D
(λ̃_fold = 0.285769734239) jako mechanizm λ̃_crit. Zmierzone
λ̃_crit = 0.1077 leży **2.65× niżej**; kolapsy zachodzą w PIERWSZYM
overshoocie nagłego załączenia źródła (t_end ∈ [0.60, 0.83]), a
znacznie lepszym deskryptorem jest 0D punkt zwrotu ze startu ψ≡1,
ψ̇=0: λ̃_turn = 1/12 = 0.0833 (zaniża próg o 23%, bo pomija ucieczkę
energii do fal). Predykcja raportowana bez reinterpretacji (LOCK §6)
— jej niepotwierdzenie jest wynikiem.

Konsekwencja, do decyzji autora: cała mapa λ̃ zbadana w tym i w
poprzednim cyklu dotyczy wyłącznie **rodziny startów „źródło
włączone skokowo w t=0"** (LOCK §3: „start zawsze: ψ≡1, π₀=0,
źródło od t=0"). Reżim λ̃ ∈ (0.1077, 0.2858) — gdzie 0D minimum
statyczne NADAL ISTNIEJE, a nasz protokół daje wyłącznie COLLAPSE —
jest nietknięty. Kandydat LOCKa: **adiabatyczne ZAŁĄCZANIE** źródła
(rampa w górę, ta sama rodzina S, Δ), pytanie binarne: czy stany
osiadłe istnieją aż do λ̃_fold i czy któryś z nich przeżywa
wygaszenie (powtórka Q-I2 na głębszych stanach). UWAGA: to pytanie
o INDUKCJĘ, nie o kreację — przy obecnym wyniku Q-I2 (zero
histerezy przy wygaszeniu, energia oddawana co do joule'a)
oczekiwanie pozostaje negatywne. [USER-GATE, kandydat nowego LOCKa]

## N5 (dotyczy core — podtrzymanie N4 poprzednika) — zakres ważności 𝒰_mat przy λ̃ ≳ 0.1

Cykl potwierdza obserwację poprzednika w ostrzejszej formie: przy
λ̃ ≳ 0.108 (ρ₀ rzędu Φ₀/qc₀ ·0.1) forma liniowa w ρ prowadzi do
ucieczki z dziedziny w czasie t < 1 j.cz., a kierunek ucieczki
(sufit vs podłoga) jest kierunkiem overshootu, nie własnością
fizyczną (λ̃=0.11 → podłoga, λ̃=0.12–0.18 → sufit, λ̃=0.20 →
BREAKDOWN). Pytanie do autora pozostaje rdzeniowe (sek08a):
czy program traktuje ten reżim jako fizyczny zakres L_mat=−(q/Φ₀)ψρ,
czy wymaga członów wyższego rzędu / samouzgodnienia ρ(ψ) (S05/L01 —
decyzja rdzeniowa, nie implementatorska). [USER-GATE]

## N6 (higiena, drobne) — co zadziałało bez zarzutu

Do powtórzenia w kolejnych LOCKach: (a) przejęcie estymatora dryfu
sekularnego poprzednika i **zamrożenie go w MD PRZED pierwszym
biegiem** — P2c PASS za pierwszym razem, zero correction notes;
(b) kotwice regresyjne z poprzednika odtworzone co do cyfry
(ψ̄(0)=0.865982 @λ̃=0.05, t_end=0.3750 @λ̃=0.5, t_end=0.600 @λ̃=0.20)
— tani i mocny dowód, że kopia silnika jest wierna; (c) bisekcja
w protokole IDENTYCZNYM z biegami głównymi (bez „krótkiego" wariantu)
— zero nowych kryteriów przy okazji. [bez user-gate; notatka
metodyczna]
