---
title: "Phase_FINAL_close — zamknięcie: Q-G-INCONCLUSIVE (zero oscylonów małej amplitudy w domenie naturalnej predykcji P1b: 9×RADIATED [τ=174–699, w tym kontrola zbieżności negatywu a=0.05 σ=6: τ=418.2 IDENTYCZNE na h i h/2], 2×COLLAPSE zbieżnie dt/2 [σ=10, a≥0.08: BOUNDARY-UPPER t≈40–55, Δt zdarzeń ≤0.06], 1×INCONCLUSIVE-RUN [a=0.05 σ=10: podtrzymanie 880>628 j.cz. i 281 przejść, ale E_end=7.8%·E_ref w szczelinie kategorii 5–50% bez plateau i bez mierzalnego ω wg zamrożonej reguły segmentu ≥1000]); deskryptywnie ω_desc≈1.0013–1.0016±0.007–0.016 na WSZYSTKICH mierzalnych biegach — rdzeń dzwoni na progu kontinuum ω≈1, ZERO śladu zmiękczenia LP (mapa przewidywała 0.9977→0.9421) i zero monotonii w a; P1a′ PASS 24/24, P2 PASS 6/6 (regresja τ=206.9 dokładnie), triage: 12/12 startów martwych w Etapie A"
date: 2026-09-15
type: phase-final-close
tgp_owner: research/op-oscillon-small-amplitude-2026-09-14
status: CLOSED
verdict: "Q-G: INCONCLUSIVE wg litery (LOCK §5: PASS wymaga ≥1 OSCILLON potwierdzonego h/2 i dt/2 — jest 0; FAIL wymaga WSZYSTKIE 12 RADIATED zbieżnie ORAZ zero OSCILLON-WEAK — jest 9 RADIATED [tam gdzie kontrole wg §3: a=0.05 σ=6 zbieżnie τ=418.2/418.2, E_end/E_ref zgodne do 6 cyfr] + 2 COLLAPSE zbieżne dt/2 [a=0.08 σ=10: BOUNDARY-UPPER 55.240 vs NONFINITE 55.245; a=0.10 σ=10: BOUNDARY-UPPER 40.010/40.065 — nadkategoria działa zgodnie z projektem N4 poprzednika] + 1 INCONCLUSIVE-RUN [a=0.05 σ=10]). Deskryptywnie (obowiązkowe, miękkie): żaden start nie utrzymał E_core≥0.5·E_ref przez 100 T₀ z mierzalnym ω_peak≤0.99; wszystkie mierzalne ω_desc = 1.0013–1.0016 (± 0.007–0.016) — dokładnie próg kontinuum m=1, BEZ obniżenia częstości przewidywanego mapą LP (0.9977/0.9855/0.9629/0.9421 dla a=0.02/0.05/0.08/0.10) i BEZ monotonii w a ⟹ w klasie zbadanej (gauss a≤0.10, σ≤10, t≤10⁴) pole NIE samopułapkuje się: rdzeń wypromieniowuje na częstości liniowej. INCONCLUSIVE ≠ pozytyw; falsyfikacja P1b NIE jest orzeczona (litera FAIL niespełniona przez 2 kolapsy szerokich startów i 1 bieg w szczelinie kategorii); ZAKAZ claimów o masach leptonów i dyskretności rodzin dotrzymany."
anti_lakatos_lock: PRESERVED
tags: [oscillon-small-amplitude, lindstedt-poincare, soft-nonlinearity, second-order-dynamics, healthy-branch, radiated, collapse, inconclusive-qg, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase_correction_note_1_fft_report.md]]"
  - "[[NEEDS.md]]"
  - "[[README.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_FINAL_close.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/NEEDS.md]]"
---

# Phase FINAL — zamknięcie cyklu op-oscillon-small-amplitude

**Status: CLOSED-EXECUTED (2026-09-15, jedna sesja: MD FROZEN →
Phase 1 → Phase 2 PASS 6/6 → Phase 3 Etap A (13 biegów) → triage
FROZEN → Etap B (vac) + potwierdzenia (h/2 obowiązkowe dla a=0.05 σ=6;
dt/2 dla obu COLLAPSE) → zamknięcie).** Kryteria LOCKa stosowane
DOSŁOWNIE; zero zmian kryteriów/progów/detektora/rodzin startów/
triage/sponge po pierwszym biegu produkcyjnym; jedna korekta —
wyłącznie warstwy raportującej (nota 1).

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| **Q-G** (oscylony małej amplitudy — RACHUNEK CENTRALNY) | **Q-G-INCONCLUSIVE** (litera §5) | 0×OSCILLON, 0×OSCILLON-WEAK, 9×RADIATED (zbieżnie tam, gdzie kontrole wg §3), 2×COLLAPSE (zbieżnie dt/2), 1×INCONCLUSIVE-RUN — litera FAIL wymaga „WSZYSTKIE 12 RADIATED", więc INCONCLUSIVE |
| P1a′ (gate cytatów form) | **PASS 24/24** | sympy(bin double) vs float ≤1.3e−15 (próg 1e−12), ψ∈{0.9,1,1.1}, 8 form + 4 tożsamości simplify=0 (w tym dU bez kancelacji — korekta 2 poprzednika dziedziczona) |
| P1b′ (cytat + mapa) | **wpisane PRZED Phase 3** | ω₂=−139/24 (CYTAT z poprzednika, bez ponownego wyprowadzania); tabela miękka ω(a): 0.9977/0.9855/0.9629/0.9421 |
| P2 (bramka maszynerii) | **PASS 6/6** | próżnia 0.0 dokładnie (obie siatki); regresja τ=206.9 (odchyłka 0.000 od baseline poprzednika mimo R=200→400 i sponge [160,200]→[320,400]); dryf energii 9.65e−8 ≤1e−6 (obie siatki, wartość zbieżna między siatkami); odbicie sponge 1.13e−7 ≤1e−3 |

**Uwagi zalockowane, stosowane dosłownie:** INCONCLUSIVE ≠ pozytyw;
mapa ω(a) pozostaje miękka (bez rangi werdyktu); ZAKAZ claimów
o masach leptonów i o dyskretności rodzin — dotrzymane.

## 1. Wejścia i maszyneria (rejestr; MD §1–§4, §9)

- Formy CYTAT (bez modyfikacji): M=ψ⁶/(4−3ψ)², M′=12ψ⁵(2−ψ)/(4−3ψ)³,
  𝒦=ψ⁴, 𝒰=ψ⁴/4−ψ³/3; K_geo=γ=c₀=1; dziedzina (0,4/3); tożsamość
  energii 𝒰(ψ)−𝒰(1)=(ψ−1)²(3ψ²+2ψ+1)/12 (korekta 2 poprzednika —
  dziedziczona).
- Silnik: KOPIA engine_core.py poprzednika (uogólniony Störmer–Verlet
  na (ψ,π), punkt stały do stagnacji maszynowej, siatka przesunięta,
  zero podłóg/barier); zmiana WYŁĄCZNIE parametryczna: sponge
  smootherstep γ₀=1.0 na [320,400]; R=400; t_max=10000; dt=0.005.
- Starty (FROZEN, π₀=0, brak seeda): 12 gaussów a∈{0.02,0.05,0.08,
  0.10}×σ∈{3,6,10} + próżnia; staging: Etap A t=2000 → triage
  E_core(2000)≥0.2·E_ref i zero zdarzeń → Etap B t=10⁴.
- Detektor (FROZEN): E_core(r≤80)≥0.5·E_ref przez ≥100 T₀ (okno od
  ≥t=50) + ≥50 przejść + ω_peak≤0.99 (FFT Hanna, segment stabilny
  ≥1000 j.cz., interpolacja paraboliczna, Δω=2π/L) + potwierdzenia
  h/2 i dt/2 (τ ±10%, ω ±2%).
- Monitoring CFL (deskryptywny): max lokalnego c=(4−3ψ)/ψ w biegach
  bez kolapsu ≤1.015 ⟹ dt·c/h ≈ 0.10 — margines bezpieczny.

## 2. Phase 3 — tabela biegów (pełna: Phase3_output.txt, Phase3_results/)

| Start | Siatka | Klasa/status | τ / t_zdarzenia | ω (deskryptywnie*) | E_end/E_ref |
|---|---|---|---|---|---|
| a=0.02 σ=3 | h05 | RADIATED | τ=174.3 (55 przejść) | – (okno 174<200) | 1.10e−3 |
| a=0.02 σ=6 | h05 | RADIATED | τ=393.9 (126) | 1.0013±0.0159* | 9.84e−3 |
| a=0.02 σ=10 | h05 | RADIATED | τ=699.2 (223) | 1.0016±0.0090* | 4.61e−2 |
| a=0.05 σ=3 | h05 | RADIATED | τ=174.9 (55) | – | 1.12e−3 |
| a=0.05 σ=6 | h05 | RADIATED | τ=418.2 (134) | 1.0013±0.0150* | 1.130e−2 |
| a=0.05 σ=6 | **h025 (kontrola obowiązkowa)** | RADIATED | **τ=418.2 (identyczne)** | 1.0013±0.0150* | 1.130e−2 (zgodne do 6 cyfr) |
| a=0.05 σ=10 | h05 | **INCONCLUSIVE-RUN** | podtrzymanie [50,930.3] = 880>628; 281 przejść | 1.0015±0.0071* | 7.78e−2 (szczelina 5–50%, bez plateau: qs=3.7) |
| a=0.08 σ=3 | h05 | RADIATED | τ=178.6 (56) | – | 1.18e−3 |
| a=0.08 σ=6 | h05 | RADIATED | τ=487.1 (156) | 1.0013±0.0129* | 1.62e−2 |
| a=0.08 σ=10 | h05 | **COLLAPSE** | t=55.240 BOUNDARY-UPPER | – | – |
| a=0.08 σ=10 | h05_dt2 | COLLAPSE (kontrola) | t=55.245 NONFINITE (Δt=0.005≤1) | – | – |
| a=0.10 σ=3 | h05 | RADIATED | τ=183.5 (58) | – | 1.25e−3 |
| a=0.10 σ=6 | h05 | RADIATED | τ=637.9 (203) | 1.0016±0.0098* | 3.07e−2 |
| a=0.10 σ=10 | h05 | **COLLAPSE** | t=40.010 BOUNDARY-UPPER | – | – |
| a=0.10 σ=10 | h05_dt2 | COLLAPSE (kontrola) | t=40.065 BOUNDARY-UPPER (Δt=0.055≤1) | – | – |
| vac (kontrola) | h05 | OK | τ cenzurowane ≥9950; 0 przejść; kandydat=False | – (sygnał zerowy) | – |

\* ω_desc: pik FFT na całym oknie podtrzymania (krótszym niż wymóg
detektora 1000 j.cz.) — deskryptywne, NIE wchodzi do detektora
(correction note 1b); ω_peak detektorowe niemierzalne dla wszystkich
biegów (żaden segment stabilny ≥1000 j.cz. nie istnieje).

**Triage (FROZEN):** WSZYSTKIE 12 startów martwe w Etapie A
(E_core(2000)/E_ref = 0.0011–0.078 < 0.2 albo zdarzenie kolapsowe);
do Etapu B przeszła wyłącznie próżnia kontrolna (litera: 0≥0.2·0,
zero zdarzeń) — dobiegła do t=10⁴ z ψ≡1 dokładnie.

**WERDYKT Q-G: INCONCLUSIVE** — litera §5: PASS wymaga ≥1 OSCILLON
(jest 0); FAIL wymaga „WSZYSTKIE 12 startów RADIATED zbieżnie ORAZ
zero OSCILLON-WEAK" (jest 9 RADIATED + 2 COLLAPSE + 1 INCONCLUSIVE-RUN).

## 3. Konfrontacja z mapą LP (deskryptywna, miękka — bez progu)

| a | ω_LP = 1−(139/24)a² | ω_desc zmierzone (σ=6 / σ=10) |
|---|---|---|
| 0.02 | 0.9977 | 1.0013±0.0159 / 1.0016±0.0090 |
| 0.05 | 0.9855 | 1.0013±0.0150 / 1.0015±0.0071 |
| 0.08 | 0.9629 | 1.0013±0.0129 / – (kolaps) |
| 0.10 | 0.9421 | 1.0016±0.0098 / – (kolaps) |

Odczyt deskryptywny: wszystkie mierzalne częstości siedzą NA progu
kontinuum ω≈1.001 (w granicach Δω zgodne z m=1), bez jakiejkolwiek
monotonii w a — to sygnatura liniowego dzwonienia dyspersyjnego
gasnącego rdzenia, nie samopułapkowania. Dla a=0.08–0.10 (σ=6)
przewidywane zmiękczenie (−3.7…−5.8%) leży daleko poza słupkiem Δω
(±1.0–1.3%) — zmierzona częstość NIE podąża za mapą LP. Zastrzeżenie
zalockowane: pomiar wykonany na gasnącym (RADIATED) rdzeniu, nie na
plateau oscylonowym — mapa LP dotyczy oscylonu, który w tej klasie
NIE powstał; konfrontacja pozostaje miękka i bez rangi werdyktu.

## 4. Mapowanie na drzewo decyzyjne (LOCK §6)

**Q-G-INCONCLUSIVE → „NEEDS metodologiczny (t_max/staging/detektor)"**
— konkretyzacja w [[NEEDS.md]]. COLLAPSE NIE dominują (2/12, wyłącznie
najszersze/najcięższe starty σ=10, a≥0.08) ⟹ warunek eskalacji do
cyklu-bliźniaka `op-collapse-matter-source` (dominacja) NIE zachodzi;
wyniki kolapsowe przekazane tam deskryptywnie (NEEDS N3). Źródłem
INCONCLUSIVE nie jest maszyneria (P2 PASS 6/6; kontrola negatywu
h/2 idealnie zbieżna; oba COLLAPSE zbieżne dt/2 ≤0.06 j.cz.), lecz
litera FAIL wymagająca 12/12 RADIATED przy: (i) 2 kolapsach szerokich
startów, (ii) 1 biegu (a=0.05 σ=10) w szczelinie kategorii (E_end
7.8% ∈ (5%,50%), podtrzymanie 880>628 j.cz., ale bez plateau —
energia w ostatnich 1000 j.cz. spada ~4.7× — i bez mierzalnego
ω_peak wg zamrożonej reguły segmentu ≥1000 j.cz.).

## 5. Korekty / incydenty / higiena (anti-Lakatos)

- ✓ LOCK przeczytany w całości przed wszystkim; MD FROZEN przed
  jakimkolwiek kodem (computations_performed: ZERO); tabela ω(a)
  zapisana w MD §6 i Phase1_output PRZED Phase 3; Phase 3 wykonana
  W CAŁOŚCI po PASS Phase 2.
- ✓ **Korekta 1** ([[Phase_correction_note_1_fft_report.md]], PRZED
  użyciem wyników; output pośredni zachowany:
  `Phase3_output_interim_stageA.txt` + `verdict_interim_stageA.json`):
  (a) fft_peak zwracał pozorny pik (moc=0) dla tożsamościowo zerowego
  sygnału próżni — po korekcie „brak piku"; ZERO wpływu na
  trajektorie/detektor/kategorie/triage; (b) addendum deskryptywne
  ω_desc (kolumna obowiązkowej tabeli konfrontacyjnej) — nie wchodzi
  do żadnego warunku detektora.
- ✓ Incydent proceduralny (udokumentowany na żądanie koordynatora):
  krótki benchmark wydajności kroku (2000 kroków, 3 s, start a=0.05
  σ=6) wykonany PO zamrożeniu MD, w trakcie biegu Phase 2, PRZED
  Phase 3; wynik użyty wyłącznie do planowania batchy (1.43 ms/krok)
  — zero wpływu na kryteria i wyniki; żaden bieg produkcyjny Phase 3
  nie wystartował przed PASS Phase 2.
- ✓ Zakaz podłóg/barier dotrzymany (granice tylko klasyfikowane);
  starty/detektor/progi/triage/sponge niezmienione po pierwszym biegu;
  INCONCLUSIVE nie reinterpretowane; mapa ω(a) miękka; zakaz claimów
  o masach i dyskretności dotrzymany.
- ✓ Rdzeń `.tex`/STATE.md/git NIETKNIĘTE; katalogi innych cykli tylko
  odczyt; pełne ścieżki bez `cd`; weryfikacja `ls` po zapisach (zero
  artefaktów zagnieżdżonych); bez /dev/null i heredoc.
- ✓ **Integralność:** SHA256 (integrity_snapshot.txt) — Phase0_balance
  .md, Phase_method_decisions.md, engine_core.py UNCHANGED przy
  zamknięciu (weryfikacja 2026-09-15, dopisana do snapshotu).
- Środowisko: CPython 3.14.2, numpy 2.4.3, scipy 1.17.1, sympy 1.14.0
  (identyczne z poprzednikiem).

## 6. Odczyt (deskryptywnie, bez claimów poza klasą zbadaną)

1. **Domena naturalna predykcji P1b została wreszcie zsondowana**
   (a∈[0.02,0.10], t do 10⁴ = 10× dłużej niż poprzednik, R=400)
   i w całej siatce 12 startów gaussowskich **nie powstał ani jeden
   oscylon ani stan quasi-stacjonarny**: 9 startów wypromieniowuje
   rdzeń w τ=174–699 j.cz. (wszystkie < progu 100 T₀ liczonego jako
   podtrzymanie od t=50; jedyny wyjątek 880 j.cz. przy a=0.05 σ=10
   również gaśnie bez plateau), a 2 najszersze kolabują do górnej
   granicy dziedziny w t≈40–55.
2. **Miękka nieliniowość LP nie zostawia śladu w dynamice pełnej:**
   wszystkie mierzalne częstości rdzenia = 1.0013–1.0016 (próg
   kontinuum), płasko w a — pole małej amplitudy zachowuje się
   liniowo-dyspersyjnie i odpływa, zanim nieliniowość zdąży związać
   rdzeń (przewidywane czasy formowania ~1/(|ω₂|a²) ≈ 17–430 j.cz.
   mieszczą się w oknie biegu, więc „za krótki bieg" nie tłumaczy
   negatywu dla a≥0.05 — deskryptywnie).
3. **Werdykt formalny to INCONCLUSIVE, nie FAIL** — o literę: 2
   kolapsy (σ=10) i 1 bieg w szczelinie kategorii blokują „wszystkie
   12 RADIATED". Merytoryczna zawartość negatywu jest jednak mocna
   i zbieżna (h/2 idealnie; dt/2 ≤0.06 j.cz.), a łącznie
   z Q-E-INCONCLUSIVE poprzednika: **w OBU klasach amplitud
   (|a|∈[0.15,0.3] i a∈[0.02,0.10]) gałąź zdrowa nie wytworzyła
   żadnego nośnika oscylonowego** — status hipotezy ratunkowej
   i łańcucha leptonowego: user-gate (NEEDS N2).

## 7. Pliki cyklu

`Phase0_balance.md` (LOCK) · `HANDOFF_PROMPT.md` ·
`Phase_method_decisions.md` (FROZEN) · `engine_core.py` (kopia;
zmiana parametryczna sponge [320,400]) · `Phase1_analytic.py` →
`Phase1_output.txt` · `Phase2_gate.py` → `Phase2_output.txt`
(+ `Phase2_progress.log`) · `Phase3_evolve.py` + `Phase3_batchA.sh`
→ `Phase3_output.txt` (+ `Phase3_output_interim_stageA.txt`,
`Phase3_batchA.log`, `Phase3_batchB.log`, `Phase3_progress.log`)
+ `Phase3_results/` (16 katalogów biegów: state.npz [ψ,π,serie
ψ(0,t),E_core(t) co 0.1, profile co 100 j.cz.] + meta.json +
analysis.json; `verdict.json` + `verdict_interim_stageA.json`) ·
`Phase_correction_note_1_fft_report.md` · `integrity_snapshot.txt` ·
`NEEDS.md` · `README.md` (log).
