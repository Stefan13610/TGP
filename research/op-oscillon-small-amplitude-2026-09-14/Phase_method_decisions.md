---
title: "Phase_method_decisions — decyzje metodyczne FROZEN cyklu op-oscillon-small-amplitude (silnik dziedziczony 1:1 z op-r3-stationary-states; zmiany TYLKO parametryczne R=400/t_max=10⁴/sponge [320,400] + staging; detektor OSCILLON z kryterium ω_peak≤0.99; tabela przewidywanych ω(a) ZAPISANA PRZED Phase 3)"
date: 2026-09-15
type: phase-method-decisions
tgp_owner: research/op-oscillon-small-amplitude-2026-09-14
status: FROZEN
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_method_decisions.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase1_output.txt]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_correction_note_2_energy_eval.md]]"
---

# Phase_method_decisions (FROZEN przed jakimkolwiek kodem cyklu)

**Status: FROZEN 2026-09-15, ZERO obliczeń cyklu przed zapisem tego
dokumentu.** Zmiana po pierwszym biegu produkcyjnym = forbidden move
(LOCK §7). Pozycje wykraczające poza literę LOCKa oznaczone
**[INPUT-MD]**. Wszystko, co nie jest niżej jawnie zmienione
parametrycznie, jest DZIEDZICZONE 1:1 z
`../op-r3-stationary-states-2026-09-14/Phase_method_decisions.md`
(dalej: MD-poprzednika).

---

## 1. Formy modelu (CYTATY — dziedziczone bez modyfikacji; LOCK §1)

Źródło pierwotne: `../op-action-audit-spectrum-insert-2026-09-13/
Phase1_output.txt` (Q-D1-PASS, konwencja |g^tt|); cytat operacyjny za
`../op-r3-stationary-states-2026-09-14/Phase1_output.txt` (linie
„Formy CYTAT (MD sec.1)…"):

- **M(ψ) = ψ⁶/(4−3ψ)²**; **M′(ψ) = 12ψ⁵(2−ψ)/(4−3ψ)³**
- **𝒦(ψ) = ψ⁴**, 𝒦′ = 4ψ³
- **𝒰(ψ) = ψ⁴/4 − ψ³/3**; 𝒰′ = ψ²(ψ−1); 𝒰″ = 3ψ²−2ψ, 𝒰″(1)=1=m²
- π = M(ψ)ψ̇; EOM radialnie 3D: Mψ̈ + ½M′ψ̇² = (1/r²)(r²𝒦ψ′)′
  − ½𝒦′ψ′² − 𝒰′
- **K_geo = γ = c₀ = 1 [LOCK]**; dziedzina ψ∈(0,4/3); próżnia ψ*=1;
  T₀=2π; 100 T₀ = 628.319; ZERO podłóg/barier.
- **Gęstość energii potencjalnej (dziedziczona korekta 2 poprzednika,
  `Phase_correction_note_2_energy_eval.md`):** tożsamość DOKŁADNA bez
  kancelacji **𝒰(ψ)−𝒰(1) = (ψ−1)²(3ψ²+2ψ+1)/12** — używana we
  WSZYSTKICH ewaluacjach energii (E, E_core) tego cyklu.

## 2. Silnik (kopia; zmiany TYLKO parametryczne — LOCK §1)

`engine_core.py` tego cyklu = **kopia**
`../op-r3-stationary-states-2026-09-14/engine_core.py` (uogólniony
Störmer–Verlet na (ψ,π), iteracje punktu stałego do stagnacji
maszynowej z awaryjnym progiem 1e−12 po FP_MAXIT=200; siatka
przesunięta r_i=(i+½)h; zerowy strumień w r=0 [ψ′(0)=0, l'Hôpital]
i r=R; siła dyskretnie wariacyjna F=−(1/(h r²))∂E_sp/∂ψ; energia
próżniowo odjęta z czynnikiem 4π i tożsamością dUfun j.w.).
Zmiany względem oryginału (wyłącznie parametryczne, LOCK §1):

- **okno sponge: R_SP0=320, R_SP1=400** (oryginał: 160/200);
  smootherstep S(x)=6x⁵−15x⁴+10x³, **γ₀=1.0** (dziedziczone), sponge
  działa wyłącznie na π (człon −γ_sp(r)π).
- Nic więcej: formuły kroku, siły, energii, pas klasyfikacyjny
  (ψ≥4/3−1e−6, ψ≤1e−6), start_gauss, start_vacuum — BEZ ZMIAN.

**Parametry produkcyjne [LOCK §3]:** h=0.05 (główna), h=0.025
(potwierdzenia); dt=0.005 (dt/2=0.0025 w potwierdzeniach); **R=400**;
**t_max=10000**; E_core: r≤80 (człony gradientowe r_{i+½}≤80);
E_ref ≔ E_core(t=50); zapis ψ(0,t)≔ψ(r₀=h/2,t) i E_core(t) co
dt_out=0.1; pełny profil co 100 j.cz.; checkpointy npz co 500 j.cz.
Przy każdym checkpoincie log: max|ψ−1| i max lokalnego c=(4−3ψ)/ψ
(monitoring CFL, deskryptywny) **[INPUT-MD]**.

## 3. Rodzina startów i staging (FROZEN [LOCK §3])

Starty deterministyczne, π₀≡0, brak seeda:
ψ(r,0)=1+a·exp(−r²/(2σ²)), **a∈{0.02,0.05,0.08,0.10} ×
σ∈{3,6,10}** (12 startów) + próżnia kontrolna ψ≡1 (13 biegów Etapu A).

- **Etap A:** wszystkie 13 biegów, h=0.05, dt=0.005, do t=2000.
- **Triage (FROZEN):** do Etapu B przechodzi bieg ŻYWY w t=2000:
  **E_core(2000) ≥ 0.2·E_ref ORAZ zero zdarzeń brzegowych** (zdarzenie
  brzegowe = wejście w pas klasyfikacyjny lub niefinityczność).
  Biegi martwe klasyfikowane z Etapu A (t_end=2000).
- **Etap B:** kontynuacja żywych do t_max=10000 z checkpointów npz.

## 4. Detektor OSCILLON i kategorie (FROZEN [LOCK §3])

Definicje pomocnicze **[INPUT-MD — operacjonalizacje]**:

- **Okna podtrzymania:** maksymalne przedziały kolejnych próbek
  t∈[50, t_end] z E_core(t) ≥ 0.5·E_ref. Okno detektora = najdłuższe
  z nich (przy remisie: najwcześniejsze). LOCKowe „okno zaczynające
  się ≥ t=50" realizowane przez skan WSZYSTKICH przedziałów w
  [50, t_end], nie tylko startującego w t=50.
- **τ (czas zejścia pod 0.5·E_ref):** τ ≔ t_hold − 50, gdzie t_hold =
  koniec pierwszego okna podtrzymania zaczynającego się w t=50
  (identycznie jak u poprzednika — wymusza to bramka regresyjna P2b,
  która żąda zgodności liczbowej z baseline τ=206.9 poprzednika
  liczonym właśnie jako t_hold−50). Jeżeli t_hold = t_end biegu
  pełnego (t_max=10000) ⟹ τ cenzurowane („≥9950").
- **Przejścia:** zmiany znaku ψ(0,t)−1 po próbkach dt_out w oknie
  detektora.
- **FFT (pomiar ω):** ψ(0,t)−⟨ψ(0,t)⟩ na segmencie stabilnym
  [max(50, T−2000), T], gdzie T = t_hold (gdy okno podtrzymania
  istnieje) albo t_end biegu (pomiar deskryptywny dla żywych);
  wymagana długość segmentu ≥1000 j.cz. (inaczej ω niemierzone ⟹
  warunek 3 niespełnialny); okno Hanna; pik dominujący z interpolacją
  paraboliczną (log-moc, 3 biny); **Δω ≔ 2π/L_seg** (szerokość binu)
  raportowana obok ω_peak; harmoniki (piki o mocy ≥1% dominującego)
  deskryptywnie.
- **Quasi-stacjonarność (dla OSCILLON-WEAK):** |ΔE_core|/E_core ≤
  1e−3 na ostatnich 1000 j.cz. operacjonalizowane jako
  |⟨E_core⟩_[t_end−100,t_end] − ⟨E_core⟩_[t_end−1000,t_end−900]| /
  ⟨E_core⟩_[t_end−100,t_end] ≤ 1e−3 (średnie okienkowe eliminują
  oscylację przepływu energii przez brzeg rdzenia).

**Detektor OSCILLON (litera LOCKa §3):**
1. E_core(t) ≥ 0.5·E_ref nieprzerwanie przez ≥100 T₀ = 628.319 j.cz.
   w oknie zaczynającym się ≥ t=50;
2. ≥50 przejść ψ(0,t) przez 1 w tym oknie;
3. **ω_peak ≤ 0.99**;
4. potwierdzenie: h=0.025 ORAZ dt/2 — czas życia ±10% (reguła
   cenzurowania: jedna wartość cenzurowana ⟹ zgodność ⟺ druga
   ≥ 0.9·9950 = 8955; obie cenzurowane ⟹ zgodne **[INPUT-MD]**),
   ω_peak ±2%.

**Kategorie (litera LOCKa §3):**
- **OSCILLON** — detektor w całości (1–4).
- **OSCILLON-WEAK** (deskryptywna, NIE pozytyw Q-G): plateau E_core
  poniżej 0.5·E_ref, ale quasi-stacjonarne (def. wyżej) ORAZ
  ω_peak≤0.99; te same potwierdzenia co OSCILLON; nie zmienia Q-G.
- **RADIATED** — E_core(t_end) < 0.05·E_ref, spadek bez plateau;
  τ j.w.
- **COLLAPSE** (nadkategoria): wejście w pas graniczny (ψ>4/3−1e−6
  lub ψ<1e−6) LUB niefinityczność, jeżeli oba zdarzenia typu
  kolapsowego mieszczą się w oknie ≤1 j.cz. między siatkami
  (operacjonalizacja **[INPUT-MD]**: bieg bazowy h=0.05 i KAŻDY
  wykonany bieg kontrolny [dt/2 obowiązkowy; h/2 jeśli wykonany]
  wykazują zdarzenie kolapsowe, a czasy zdarzeń różnią się parami
  ≤1 j.cz.); podtyp BOUNDARY-UPPER / BOUNDARY-LOWER / NONFINITE
  raportowany deskryptywnie.
- **INCONCLUSIVE-RUN** — kategoria niezbieżna między siatkami (poza
  regułą nadkategorii COLLAPSE) lub bieg niespełniający żadnej
  z powyższych liter.

**Potwierdzenia (reguła FROZEN [LOCK §3]):** (i) każdy kandydat
OSCILLON/OSCILLON-WEAK — pełne potwierdzenie h=0.025 ORAZ dt/2;
(ii) obowiązkowo start (a=0.05, σ=6) — h=0.025 niezależnie od klasy
(kontrola zbieżności negatywu); (iii) dt/2 dla kandydatów i przy
zdarzeniach COLLAPSE. Horyzont biegu kontrolnego **[INPUT-MD]**:
lustrzany do bazowego (ten sam t_end wg stagingu zastosowanego do
samego biegu kontrolnego; dla kontroli COLLAPSE wystarcza
t ≤ min(t_end bazowego, t_zdarzenia+100)). dt/2 wykonywane na h=0.05
(jak u poprzednika **[INPUT-MD]**).

## 5. Werdykt Q-G (litera LOCKa §5 — kopia bez zmian)

- **Q-G-PASS:** ≥1 start z klasą OSCILLON potwierdzoną (h/2 i dt/2).
- **Q-G-FAIL:** WSZYSTKIE 12 startów RADIATED zbieżnie (kategoria
  i τ ±10% na h i h/2 w zakresie potwierdzeń §3) ORAZ zero
  OSCILLON-WEAK. = falsyfikacja P1b w klasie zbadanej (a≤0.10, t≤10⁴).
- **Q-G-INCONCLUSIVE:** każdy inny rozkład kategorii.

## 6. Tabela przewidywanych ω(a) (MIĘKKA — zapisana PRZED Phase 3; LOCK §2)

Cytat P1b (poprzednik, `Phase1_output.txt`): **ω₂ = −139/24 ≈
−5.791667 < 0** (bez ponownego wyprowadzania). Mapa dla modu
jednorodnego ω(a) ≈ 1 + ω₂a² (deskryptywna, BEZ progu PASS/FAIL;
pre-rejestrowany jest ZNAK ω_peak<1 i monotonia w a):

| a | ω(a) = 1 − (139/24)a² |
|---|---|
| 0.02 | 0.997683 |
| 0.05 | 0.985521 |
| 0.08 | 0.962933 |
| 0.10 | 0.942083 |

Uwaga deskryptywna (nie zmienia progów): dla a=0.02 mapa LP daje
ω=0.9977 > 0.99, tj. powyżej zamrożonego progu detektora (warunek 3);
próg pozostaje literą LOCKa.

## 7. Phase 1 (analityczna, przed numeryką; LOCK §2)

- **P1a′:** gate cytatów form — sympy (na dokładnej binarnej
  reprezentacji double wejścia, lekcja poprzednika) vs float dla
  M, M′, 𝒦, 𝒦′, 𝒰, 𝒰′, 𝒰″ oraz tożsamości energetycznej
  dU=(ψ−1)²(3ψ²+2ψ+1)/12 w ψ∈{0.9,1,1.1}; próg |Δ| ≤
  1e−12·max(1,|wartość|); dodatkowo sympy simplify=0 dla
  dU−(𝒰(ψ)−𝒰(1)) i dla cytatu M′ **[INPUT-MD — zestaw form]**.
- **P1b′:** wpis ω₂=−139/24 jako CYTAT (źródło: poprzednik
  Phase1_output.txt) + tabela §6 do outputu.

## 8. Phase 2 — bramka maszynerii (FROZEN; FAIL ⟹ STOP; LOCK §4)

- **P2a (próżnia):** konfiguracja produkcyjna (R=400, sponge ON),
  ψ≡1, π≡0, 100 T₀=628.319, obie siatki (h=0.05, h=0.025); gate:
  ‖ψ−1‖∞ ≤ 1e−10 przez cały bieg; detektor zero alarmów.
- **P2b (regresja):** start (a=+0.15, σ=3), R=400, h=0.05, sponge ON,
  do t=300; gate: **τ (= t_hold−50) zgodne z baseline poprzednika
  206.9 do ±5%** (dopuszczalna różnica od R=200→400 i innego okna
  sponge).
- **P2c-energia:** start (a=0.05, σ=6) — start produkcyjny — pudło
  zamknięte (**sponge OFF [INPUT-MD]**, R=400), do t=700, obie
  siatki; dryf ≔ |⟨E⟩_{t∈[90T₀,100T₀]} − ⟨E⟩_{t∈[0,10T₀]}| /
  ⟨E⟩_{[0,10T₀]} ≤ **1e−6** (operacjonalizacja średnimi okienkowymi
  dziedziczona z MD-poprzednika §7); deskryptywnie max|E−E(0)|/E(0).
- **P2c-sponge (odbicie):** bieg A: R=400 sponge ON; bieg B
  (referencja): R=800 sponge OFF; h=0.05, t=350; start wspólny
  **[INPUT-MD — skalowanie geometrii poprzednika ×2]:**
  ψ = 1+1e−3·exp(−(r−200)²/(2·5²)); w zmiennej u=r(ψ−1):
  gate: max_{r≤240,t≤350}|u_A−u_B| / max_{r∈[240,320],t≤350}|u_B|
  ≤ **1e−3** (mianownik = amplituda padająca na wejściu warstwy
  sponge [320,400]; ściana R=800 nie zawraca sygnału do r≤240 przed
  t=350).
- FAIL któregokolwiek ⟹ STOP (LOCK §4).

## 9. Rejestr WEJŚĆ (flagowane)

[LOCK]: K_geo=γ=c₀=1; formy §1 (cytaty); rodzina startów §3;
h∈{0.05,0.025}; R=400; sponge [320,400] γ₀ smootherstep; dt=0.005
(dt/2=0.0025); t_max=10000; staging Etap A t=2000, triage
0.2·E_ref + zero zdarzeń brzegowych; E_core r≤80; E_ref=E_core(50);
detektor: 0.5·E_ref, 100 T₀, ≥50 przejść, ω_peak≤0.99, potwierdzenia
h/2+dt/2 (τ ±10%, ω ±2%); kategorie OSCILLON/OSCILLON-WEAK
(1e−3/1000 j.cz.)/RADIATED (0.05·E_ref)/COLLAPSE (pas 4/3−1e−6,
1e−6, okno ≤1 j.cz.)/INCONCLUSIVE-RUN; P2a 1e−10; P2b 206.9±5%;
P2c 1e−6 i 1e−3 (vs R=800, t=350); dt_out=0.1; profil co 100 j.cz.;
checkpointy co 500 j.cz.; brak seeda; werdykty §5.
[INPUT-MD]: γ₀=1.0 (dziedziczone); sponge OFF w P2c-energia i biegu
referencyjnym P2c-sponge; puls odbiciowy a=1e−3 σ=5 @ r=200, strefy
r≤240 / [240,320]; τ ≔ t_hold−50 (konwencja poprzednika, wymuszona
przez P2b); skan wszystkich okien podtrzymania w [50,t_end];
reguła cenzurowania (≥9950, próg 8955); segment FFT
[max(50,T−2000),T], długość ≥1000, interpolacja paraboliczna,
Δω=2π/L_seg; quasi-stacjonarność średnimi okienkowymi 100 j.cz.;
operacjonalizacja okna ≤1 j.cz. COLLAPSE (pary bazowy–kontrolne);
horyzont biegów kontrolnych lustrzany (COLLAPSE: t_zdarzenia+100);
dt/2 na h=0.05; ψ(0,t)≔ψ(h/2,t); log max|ψ−1| i max c przy
checkpointach; tolerancja punktu stałego = stagnacja maszynowa
(FP_MAXIT=200, próg awaryjny 1e−12 — dziedziczone z korekty 1(b)
poprzednika); P1a′ zestaw form z tożsamością dU.

**FROZEN. Zmiany poniżej po starcie obliczeń = forbidden move.**
