---
title: "NEEDS — op-oscillon-small-amplitude: Q-G-INCONCLUSIVE ⟹ NEEDS metodologiczny (drzewo LOCK §6) + user-gate statusu łańcucha leptonowego (negatyw merytorycznie mocny w OBU klasach amplitud)"
date: 2026-09-15
type: needs
tgp_owner: research/op-oscillon-small-amplitude-2026-09-14
status: USER-GATED
related:
  - "[[Phase_FINAL_close.md]]"
  - "[[Phase0_balance.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/NEEDS.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/README.md]]"
---

# NEEDS (wszystkie pozycje user-gated; zero samowolnych działań)

Werdykt cyklu: **Q-G-INCONCLUSIVE** (0×OSCILLON, 0×OSCILLON-WEAK,
9×RADIATED zbieżnie w zakresie potwierdzeń, 2×COLLAPSE zbieżnie dt/2,
1×INCONCLUSIVE-RUN). Wg drzewa LOCKa §6: „NEEDS metodologiczny
(t_max/staging/detektor)". Maszyneria NIE jest źródłem INCONCLUSIVE
(P2 PASS 6/6, w tym regresja τ=206.9 z odchyłką 0.000; kontrola
negatywu h/2: τ=418.2 identyczne, E_end/E_ref zgodne do 6 cyfr;
COLLAPSE zbieżne dt/2 do ≤0.06 j.cz.). Warunek eskalacji do
cyklu-bliźniaka (dominacja COLLAPSE) NIE zachodzi (2/12).

## N1 (metodologiczny — litera drzewa): szczelina kategorii przy (a=0.05, σ=10)

Jedyny bieg INCONCLUSIVE-RUN spełnia warunki 1–2 detektora
(podtrzymanie 880 > 628 j.cz., 281 przejść), ale: E_core(2000) =
7.8%·E_ref leży w szczelinie między progiem RADIATED (5%) a progiem
plateau (50%); brak quasi-stacjonarności (spadek ~4.7× na ostatnich
1000 j.cz. — to gasnący ogon, nie plateau); ω_peak niemierzalne wg
zamrożonej reguły (segment stabilny ≥1000 j.cz. nie istnieje przy
podtrzymaniu 880). Propozycje do EWENTUALNEGO nowego LOCKa (zmiana
definicji możliwa tylko tam): (a) klasyfikacja biegów z podtrzymaniem
>100 T₀ dopiero z pełnego t_max (bez uśmiercania w triage t=2000 —
tu triage 0.2·E_ref i tak by go uśmiercił: 0.078<0.2, więc zmiana
dotyczyłaby progu triage), (b) domknięcie szczeliny kategorii
(RADIATED do 0.5·E_ref przy monotonicznym spadku), (c) segment FFT
dopasowany do okna podtrzymania. **Decyzja: user.**

## N2 (statusowy, najważniejszy merytorycznie): łańcuch leptonowy / hipoteza ratunkowa po DWÓCH cyklach

Predykcja P1b (ω₂=−139/24<0) została skonfrontowana w swojej domenie
naturalnej (a≤0.10, t≤10⁴, 10×dłużej niż poprzednik): **zero
oscylonów, zero stanów quasi-stacjonarnych**; wszystkie mierzalne
częstości rdzenia = 1.0013–1.0016 (próg kontinuum, płasko w a) —
deskryptywnie zero śladu zmiękczenia LP w pełnej dynamice, mimo że
przewidywane czasy formowania (~17–430 j.cz.) mieszczą się głęboko
w oknie biegów. Formalnie NIE jest to falsyfikacja P1b (litera FAIL
niespełniona — patrz werdykt), ale łącznie z Q-E-INCONCLUSIVE
poprzednika: **gałąź zdrowa nie wykazała nośnika oscylonowego
w ŻADNEJ z dwóch zbadanych klas amplitud** (|a|∈[0.15,0.30]
i a∈[0.02,0.10]; kształty gauss σ∈[3,10] + quasi-R3). Pytanie
user-gated (drzewo §6 gałąź FAIL, tu adekwatne co do treści przy
literze INCONCLUSIVE): status hipotezy ratunkowej „profile R3 jako
stany gałęzi zdrowej" i łańcucha leptonowego; ewentualny dopisek core
`rem:psi-EOM-R3-branch-status` — **wyłącznie decyzją usera** (rdzeń
.tex nietknięty przez ten cykl).

## N3 (deskryptywny, transfer do bliźniaka): kolaps szerokich startów małej amplitudy

Nowa obserwacja: nawet przy amplitudzie a=0.08–0.10 start dostatecznie
SZEROKI (σ=10; E_ref≈19.2 i 33) kolabuje do GÓRNEJ granicy dziedziny
(BOUNDARY-UPPER, t≈40–55, zbieżnie dt/2; przy dt/2 raz NONFINITE
w oknie 0.005 j.cz. — nadkategoria COLLAPSE z N4 poprzednika działa
zgodnie z projektem). Sferyczna implozja wzmacnia ψ w rdzeniu ponad
4/3 niezależnie od małości a, gdy całkowita energia jest duża.
Materiał wejściowy dla cyklu-bliźniaka
`op-collapse-matter-source-2026-09-14` (Q-H2) — bez eskalacji
(COLLAPSE nie dominują: 2/12). **Decyzja: user.**

## N4 (techniczny, drobny): budżet detektora ω dla krótkich podtrzymań

Zamrożony wymóg segmentu stabilnego ≥1000 j.cz. czyni warunek 3
detektora niespełnialnym dla obiektów żyjących 628–1050 j.cz. —
w tym cyklu bez konsekwencji (żaden bieg nie miał plateau), ale
w przyszłym LOCKu warto związać długość segmentu z oknem podtrzymania
(np. min(1000, 0.8·τ) z odpowiednio raportowanym Δω). **Decyzja:
user (tylko nowy LOCK).**

## Dotrzymane zakazy

ZAKAZ claimów o masach leptonów i o dyskretności rodzin — dotrzymany;
INCONCLUSIVE ≠ pozytyw — dotrzymane (werdykt wg litery, bez
reinterpretacji); mapa ω(a) pozostała miękka (konfrontacja wyłącznie
deskryptywna); starty/detektor/progi/triage/sponge niezmienione po
pierwszym biegu produkcyjnym; korekta 1 wyłącznie warstwy raportu
(nota PRZED użyciem, pierwotny output pośredni zachowany).
