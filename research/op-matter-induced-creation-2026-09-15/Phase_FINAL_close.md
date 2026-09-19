---
title: "Phase_FINAL_close — zamknięcie: Q-I1-INCONCLUSIVE (stan PODPROGOWY osiadły ISTNIEJE na siatce produkcyjnej — λ̃=0.08: ψ̄(0)=0.80803, λ̃=0.10: ψ̄(0)=0.77290, oba <5/6 — ale LOCKowe potwierdzenie h=0.025 dla NAJGŁĘBSZEGO przypadku (λ̃=0.10) ROZJECHAŁO SIĘ [COLLAPSE t=0.95]; λ̃_crit(h=0.05)=0.107656±0.000156) + Q-I2-FAIL (KREACJI NIE MA: wszystkie 4 stany osiadłe po adiabatycznym wygaszeniu źródła dają RETURN-TO-VACUUM, kontrola negatywu h=0.025 zgodna, kontrola czystości λ̃=0.01 PASS; predykcja pre-rejestrowana TRAFIONA — obiekt indukowany jest CIENIEM źródła, zero histerezy)"
date: 2026-09-15
type: phase-final-close
tgp_owner: research/op-matter-induced-creation-2026-09-15
status: CLOSED
verdict: "Q-I1: INCONCLUSIVE wg litery (LOCK §4: PASS wymaga ≥1 λ̃ z listy z SETTLED-SUB ZBIEŻNIE — „kategoria na h i h=0.025 dla najgłębszego przypadku"; najgłębszy SETTLED-SUB to λ̃=0.10 [ψ̄(0)=0.772903 na h=0.05], a jego potwierdzenie h=0.025 dało COLLAPSE t_end=0.95 ⟹ brak zbieżnego SETTLED-SUB; FAIL też nie zachodzi, bo lista NIE jest w całości SETTLED-DEF/COLLAPSE: 0.08 i 0.10 to SETTLED-SUB na siatce produkcyjnej. Obraz fazy włączonej: λ̃∈{0.01,0.05,0.06} SETTLED-DEF, {0.08,0.10} SETTLED-SUB, {0.12,0.14,0.16,0.18,0.20} COLLAPSE w pierwszym overshoocie [t_end 0.83/0.75/0.695/0.65/0.60 — kotwica 0.20 odtwarza poprzednika co do cyfry]; bisekcja 6 kroków: λ̃_crit=0.107656±0.000156 na h=0.05, przy h=0.025 λ̃_crit<0.10 [rozjazd = zawężenie okna, nie jego zniknięcie]). Q-I2: FAIL wg litery (LOCK §4: wszystkie biegi SETTLED-* po wygaszeniu dają RETURN-TO-VACUUM lub COLLAPSE zbieżnie — 4/4 RETURN-TO-VACUUM [λ̃=0.05,0.06,0.08,0.10], obowiązkowa kontrola negatywu λ̃=0.05 na h=0.025 ZGODNA, kontrola czystości wygaszania λ̃=0.01 PASS [max|ψ−1|=1.38e−4 < 1e−3]; PERSISTENT-OBJECT = 0). Mechanizm deskryptywny: rampa smootherstep Δ=100 jest dla tego układu adiabatyczna — już w t=700 max_{r≤40}|ψ−1| spada z 0.227 do 9.8e−6, a E_core z −9.15 do +4.1e−4; po wygaszeniu nie zostaje NIC zlokalizowanego (max|ψ−1|≤2.5e−4 do t=1700). Predykcja pre-rejestrowana Q-I2 („RETURN-TO-VACUUM lub COLLAPSE") TRAFIONA. Predykcja mechanizmu P1-I2 (λ̃_fold jako mechanizm λ̃_crit) NIE potwierdzona deskryptywnie: λ̃_fold(0D)=0.285770 ≫ λ̃_crit(dyn)=0.1077 — próg jest dynamiczny (overshoot nagłego załączenia), nie statyczny; ψ̄(0) zmierzone leży systematycznie POWYŻEJ ψ_min(0D) o +0.0075…+0.0231 (sztywność gradientowa). P2 PASS 3/3 bez korekt (ψ̄(0)@0.05=0.865982 i t_end@0.5=0.3750 — kotwice poprzednika co do cyfry). INCONCLUSIVE ≠ pozytyw; zakazy claimów (masy leptonów, oscylony, ρ(ψ)) dotrzymane."
anti_lakatos_lock: PRESERVED
claim_status: "B"
tags: [matter-induced-creation, matter-coupling, L-mat-unified, subthreshold-state, settled-sub, ramp-off, return-to-vacuum, no-creation, lambda-crit, grid-divergence, second-order-dynamics, healthy-branch, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[NEEDS.md]]"
  - "[[README.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# Phase FINAL — zamknięcie cyklu op-matter-induced-creation

**Status: CLOSED-EXECUTED (2026-09-15, jedna sesja: MD FROZEN →
integrity_snapshot → Phase 1 (sympy/0D) → Phase 2 (bramka, PASS 3/3,
ZERO korekt) → Phase 3 Q-I1 (10 biegów + 6 kroków bisekcji +
1 potwierdzenie h/2) → Phase 3 Q-I2 (5 biegów wygaszania +
1 potwierdzenie h/2) → probe deskryptywny zbieżności → zamknięcie).**
Kryteria LOCKa stosowane DOSŁOWNIE; **zero zmian kryteriów / progów /
detektorów / listy λ̃ / rampy / okien / sponge po pierwszym biegu
produkcyjnym; zero correction notes** (żaden gate nie wymagał
diagnostyki).

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| **Q-I1** (indukowany obiekt: czy istnieje trwały stan podprogowy ze źródłem WŁĄCZONYM) | **Q-I1-INCONCLUSIVE** (litera LOCK §4) | Stan podprogowy osiadły ISTNIEJE na siatce produkcyjnej (λ̃=0.08: ψ̄(0)=0.80803; λ̃=0.10: ψ̄(0)=0.77290; próg 5/6=0.83333), ale LOCKowe potwierdzenie h=0.025 dotyczy **najgłębszego** przypadku (λ̃=0.10) i **rozjechało się** (COLLAPSE t=0.95) ⟹ brak zbieżnego SETTLED-SUB; FAIL też nie zachodzi |
| **Q-I2** (KREACJA — centralne: czy obiekt przeżywa wygaszenie źródła) | **Q-I2-FAIL** (litera LOCK §4) | 4/4 stany SETTLED-* → **RETURN-TO-VACUUM**; PERSISTENT-OBJECT = 0; kontrola negatywu (λ̃=0.05, h=0.025) zgodna; **kreacji nie ma — obiekt jest cieniem źródła** |
| P1-I1 (gate cytatów form + rampy) | **PASS** | 6/6 tożsamości simplify=0; 27/27 kontroli numerycznych ≤1.4e−15 (próg 1e−12); λ̃(t)∈{λ̃,λ̃,λ̃/2,0,0} w t∈{0,600,650,700,1700} do 1e−15; λ̃=0 **dokładnie** od t=700 |
| P1-I2 (PRE-REJESTRACJA ilościowa 0D) | **ZAPISANE PRZED Phase 3** | ψ_min(λ̃) dla całej listy + λ̃_fold=**0.285769734239** (ψ_fold=0.327917), sympy=numeryka do 5.6e−17 |
| P1-I3 (kontekst predykcji Q-I2) | **PASS** (cytat) | bez źródła jedynym znanym stanem trwałym jest ψ≡1 (Q-A-PASS + Q-B-FAIL M911) |
| P2 (bramka) | **PASS 3/3**, bez korekt | próżnia 0.000e+00 (≤1e−10, obie siatki); regresja ψ̄(0)=0.865982 (kotwica 0.865982, odchył 0.000%); regresja COLLAPSE t_end=0.3750 (kotwica 0.375); dryf sekularny 8.80e−09 / 8.99e−09 (≤1e−6/100T₀) |
| Kontrola czystości wygaszania (λ̃=0.01) | **PASS** | RETURN-TO-VACUUM; max_{r≤40}|ψ−1|(1700)=1.38e−4 < 1e−3 |
| Predykcja pre-rejestrowana Q-I2 | **TRAFIONA** | „RETURN-TO-VACUUM lub COLLAPSE" — zrealizowane 4/4 jako RETURN-TO-VACUUM |
| Predykcja mechanizmu P1-I2 (λ̃_fold ⟶ λ̃_crit) | **NIE POTWIERDZONA** (deskryptywnie) | λ̃_fold=0.2858 ≫ λ̃_crit=0.1077; próg jest **dynamiczny**, nie statyczny (raport bez reinterpretacji — LOCK §6) |

**Uwagi zalockowane, stosowane dosłownie:** INCONCLUSIVE ≠ pozytyw
(Q-I1 **nie jest** „kreacją indukowaną" — jest niezbieżnym kandydatem);
ZAKAZ claimów o masach leptonów i oscylonach dotrzymany; S05
dotrzymane (ρ̂ statyczne, λ̃(t) wyłącznie wg zamrożonej rampy, zero
pól dynamicznych, zero ρ(ψ)); rdzeń `.tex`, STATE.md i git nietknięte.

## 1. Wejścia i formy (rejestr [INPUT]; MD §1–§3, §9)

- **Sektor pola [CYTAT]** (za `op-collapse-matter-source` ⟵
  `op-r3-stationary-states`, Q-D1-PASS, konwencja |g^tt|):
  M=ψ⁶/(4−3ψ)², M′=12ψ⁵(2−ψ)/(4−3ψ)³, 𝒦=ψ⁴, 𝒦′=4ψ³,
  𝒰=ψ⁴/4−ψ³/3, 𝒰′=ψ²(ψ−1); K_geo=γ=c₀=1; dziedzina (0,4/3);
  pas 4/3−1e−6 / 1e−6; ZERO podłóg i barier;
  𝒰−𝒰(1)=(ψ−1)²(3ψ²+2ψ+1)/12.
- **Człon materii [CYTAT wyniku P1-H1 poprzednika]** (nie wyprowadzany
  ponownie): 𝒰_mat=λ̃(t)ρ̂(r)ψ²/(4−3ψ);
  ∂𝒰_mat/∂ψ=λ̃ρ̂ψ(8−3ψ)/(4−3ψ)²; 𝒰_mat−𝒰_mat(1)=λ̃ρ̂(ψ−1)(ψ+4)/(4−3ψ);
  ρ̂=exp(−r²/18). Geneza korpusowa: `eq:L-mat-unified`
  L_mat=−(q/Φ₀)ψρ × √−g=c₀ψ/(4−3ψ) (rdzeń, TYLKO ODCZYT).
- **Rampa [LOCK, JEDYNA zmiana merytoryczna silnika]:**
  λ̃(t)=λ̃·S((t_off+Δ−t)/Δ), S=smootherstep, Δ=100, t_off=600 ⟹
  λ̃(t)=λ̃ dla t≤600, λ̃=0 **dokładnie** dla t≥700. Implementacja:
  F₁ liczone z λ̃(t), F₂ z λ̃(t+dt); dla λ̃=const krok jest bitowo
  identyczny z silnikiem poprzednika (dowód empiryczny: P2b odtwarza
  obie kotwice **co do ostatniej podanej cyfry**).
- **Protokół [LOCK]:** h=0.05 (potwierdzenia 0.025), R=200, dt=0.005
  (dt/2 w potwierdzeniach), sponge smootherstep γ₀=1 na [160,200],
  E_core r≤80, dt_out=0.1, profil co 5 j.cz., checkpoint npz co 100
  j.cz.; start ZAWSZE ψ≡1, π₀=0, źródło od t=0; t_on=600, t_max=1700.
- **Listy [LOCK]:** Q-I1 λ̃∈{0.06,0.08,0.10,0.12,0.14,0.16,0.18}
  + kotwice {0.05,0.20}; kontrola wygaszania λ̃=0.01 (poza werdyktami).

## 2. Phase 1 — analityka pre-rejestrowana (`Phase1_analytic.py` → `Phase1_output.txt`)

**P1-I1 PASS.** Tożsamości symboliczne (6/6, simplify=0) i kontrola
numeryczna sympy↔silnik dla ψ∈{0.9,1,1.1}: maks. odchył 1.33e−15
przy progu 1e−12. Gate rampy: λ̃(0)=λ̃(600)=λ̃, λ̃(650)=λ̃/2 (dokładnie
0.5 — S(0.5)=1/2), λ̃(700)=λ̃(1700)=0 **dokładnie**; S′(0)=S′(1)=0.

**P1-I2 (PRE-REJESTRACJA, zapisana PRZED Phase 3).** Minimum
algebraiczne f(ψ)=𝒰(ψ)+λ̃ψ²/(4−3ψ) przy ρ̂=1; gałąź f′=0 daje
**λ̃ = LAM(ψ) = ψ(1−ψ)(4−3ψ)²/(8−3ψ)** (gate simplify=0), więc fold
(f′=f″=0) = maksimum LAM na (0,1):

> **ψ_fold = 0.327916542814, λ̃_fold = 0.285769734239** (sympy vs
> numeryka: |Δ| = 5.6e−17). Dla λ̃ > λ̃_fold lokalne minimum 0D nie
> istnieje.

| λ̃ | 0.06 | 0.08 | 0.10 | 0.12 | 0.14 | 0.16 | 0.18 | *(0.01)* | *(0.05)* | *(0.20)* |
|---|---|---|---|---|---|---|---|---|---|---|
| ψ_min(0D) | 0.826612 | 0.786726 | 0.749756 | 0.714585 | 0.680435 | 0.646671 | 0.612685 | 0.957825 | 0.848184 | 0.577787 |

(kursywą kotwice — deskryptywne rozszerzenie tabeli [INPUT-MD]).
Cała lista LOCKa ma ψ_min < 5/6, tj. 0D **przewidywało** stan
podprogowy w całym oknie.

Dodatek deskryptywny [INPUT-MD, dopisany PRZED Phase 3]: 0D
zachowawczy punkt zwrotu ze startu ψ≡1, ψ̇=0 (ΔV=0) —
λ̃_turn(ψ→0⁺)=**1/12=0.083333**, tj. w czystym 0D każda λ̃>0.0833
sięga podłogi w pierwszym overshoocie.

**P1-I3 PASS** (cytat): bez źródła jedynym znanym stanem trwałym
gałęzi zdrowej jest ψ≡1 (M911: Q-A-PASS + Q-B-FAIL; brak oscylonów
Q-E/Q-G) — kontekst predykcji Q-I2.

## 3. Phase 2 — bramka (`Phase2_gate.py` → `Phase2_output.txt`) — PASS 3/3

| gate | wynik | próg | status |
|---|---|---|---|
| P2a próżnia h=0.05 / h=0.025 | 0.000e+00 / 0.000e+00 | ≤1e−10 | PASS |
| P2b-1 regresja λ̃=0.05: ψ̄(0) | **0.865982** (δψ(0)=−0.134018) | 0.865982 ±1% | PASS (odchył 0.000%) |
| P2b-2 regresja λ̃=0.5: t_end | **0.3750** (BREAKDOWN-BOUNDARY) | 0.375 ±5% | PASS |
| P2c dryf sekularny h=0.05 / h=0.025 | 8.801e−09 / 8.993e−09 | ≤1e−6/100T₀ | PASS |

Deskryptywnie (P2c): stary estymator offsetu ΔC₂ = 2.282e−06/2.298e−06
(dokładnie odtwarza wartości poprzednika — potwierdzenie, że korekta 1
poprzednika była poprawną diagnozą, a nie dopasowaniem); fit LSQ
[10,100]T₀ = 1.61e−07; max|E−E(0)| = 3.38e−04. **Żadna korekta nie
była potrzebna w tym cyklu** (estymator sekularny przejęty i zamrożony
w MD §5 PRZED pierwszym biegiem).

## 4. Phase 3 / Q-I1 — faza włączona (`Phase3_qi1_settle.py` → `Phase3_qi1_output.txt`)

### 4.1 Tabela główna (h=0.05, dt=0.005, okno [500,600])

| λ̃ | rola | klasa | ψ̄(0) zmierz. | ψ_min(0D) P1-I2 | Δ(meas−0D) | min ψ w transjencie | t_end |
|---|---|---|---|---|---|---|---|
| 0.01 | kontrola | SETTLED-DEF | 0.965353 | 0.957825 | +0.007528 | 0.93135 | — |
| 0.05 | kotwica | SETTLED-DEF | 0.865982 | 0.848184 | +0.017798 | 0.70638 | — |
| 0.06 | **lista** | SETTLED-DEF | 0.845733 | 0.826612 | +0.019121 | 0.65244 | — |
| 0.08 | **lista** | **SETTLED-SUB** | **0.808031** | 0.786726 | +0.021305 | 0.52822 | — |
| 0.10 | **lista** | **SETTLED-SUB** | **0.772903** | 0.749756 | +0.023148 | 0.42689 | — |
| 0.12 | **lista** | COLLAPSE (BB-górny) | n/a | 0.714585 | — | — | 0.8300 |
| 0.14 | **lista** | COLLAPSE (BB-górny) | n/a | 0.680435 | — | — | 0.7500 |
| 0.16 | **lista** | COLLAPSE (BB-górny) | n/a | 0.646671 | — | — | 0.6950 |
| 0.18 | **lista** | COLLAPSE (BB-górny) | n/a | 0.612685 | — | — | 0.6500 |
| 0.20 | kotwica | COLLAPSE (BREAKDOWN) | n/a | 0.577787 | — | — | **0.6000** |

Osiadłość wszystkich SETTLED-* (biegi główne): V ∈ [2.5e−5, 3.6e−5] przy progu
0.01·D (D ∈ [0.035, 0.227]) — margines 2–3 rzędy wielkości.
Kotwica λ̃=0.20 odtwarza poprzednika (Q-H1: „lam=0.2 … COLLAPSE
t_end=0.600") **co do cyfry**.

### 4.2 Bisekcja λ̃_crit (6 kroków, protokół identyczny, FROZEN)

Przedział startowy [0.100, 0.120] (najwyższe bez COLLAPSE / najniższe
z COLLAPSE):

| krok | λ̃ | wynik | przedział po kroku |
|---|---|---|---|
| 1 | 0.1100000 | COLLAPSE (BB-dolny) | [0.100000, 0.110000] |
| 2 | 0.1050000 | SETTLED-SUB (ψ̄(0)=0.764405) | [0.105000, 0.110000] |
| 3 | 0.1075000 | **UNSETTLED** (min ψ=0.2747) | [0.107500, 0.110000] |
| 4 | 0.1087500 | COLLAPSE (BREAKDOWN) | [0.107500, 0.108750] |
| 5 | 0.1081250 | COLLAPSE (BB-górny) | [0.107500, 0.108125] |
| 6 | 0.1078125 | COLLAPSE (BB-górny) | [0.107500, 0.107813] |

> **λ̃_crit (h=0.05) = 0.1076563 ± 0.0001562** (deskryptywnie).

Zgodne z dziedziczonym oknem λ̃_crit∈(0.05, 0.2] i **zawężające je
~64×**. Uwaga deskryptywna: krok 3 (λ̃=0.1075) dał UNSETTLED — tuż
pod progiem transjent nie osiada do t=600 (kategoria niezerowa
w regule bisekcji „COLLAPSE vs nie-COLLAPSE", bez wpływu na regułę).

### 4.3 Potwierdzenie siatkowe (LOCK §3) i litera werdyktu

LOCK §3 nakazuje potwierdzenie h=0.025 dla **najgłębszego
SETTLED-SUB** ⟹ λ̃=0.10. Wynik: **COLLAPSE, t_end=0.95
(BREAKDOWN-BOUNDARY górny)** — kategoria **ROZJECHANA** względem
h=0.05 (SETTLED-SUB). Zatem:

- **Q-I1-PASS nie zachodzi** (LOCK §4 wymaga SETTLED-SUB „zbieżnie —
  kategoria na h i h=0.025 dla najgłębszego przypadku");
- **Q-I1-FAIL nie zachodzi** (lista nie jest w całości
  SETTLED-DEF/COLLAPSE: 0.08 i 0.10 to SETTLED-SUB na siatce
  produkcyjnej);
- ⟹ **Q-I1-INCONCLUSIVE**.

### 4.4 Konfrontacja deskryptywna z P1-I2 (bez reinterpretacji)

- ψ̄(0) zmierzone leży **systematycznie POWYŻEJ** ψ_min(0D), z luką
  rosnącą monotonicznie z λ̃: +0.0075 (0.01) → +0.0231 (0.10).
  Odczyt: skończona szerokość źródła (σ_ρ=3) i sztywność gradientowa
  (𝒦ψ′²) spłycają studnię względem granicy 0D.
- **Mechanizm przewidziany przez P1-I2 (fold) NIE jest mechanizmem
  λ̃_crit:** λ̃_fold=0.285770 vs λ̃_crit(dyn)=0.107656 — próg leży
  2.65× niżej. Kolapsy zachodzą w PIERWSZYM overshoocie nagłego
  załączenia (t_end ∈ [0.60, 0.83]), zanim jakiekolwiek minimum
  statyczne mogłoby być testowane. Deskryptywnie znacznie bliższa
  jest granica **0D zachowawcza** λ̃_turn=1/12=0.0833 (zaniża próg
  o 23%, bo pomija ucieczkę energii do fal — promieniowanie ratuje
  biegi 0.08–0.107 przed podłogą). **Predykcja P1-I2 raportowana bez
  przeformułowania (LOCK §6); jej niepotwierdzenie jest wynikiem.**

## 5. Phase 3 / Q-I2 — wygaszanie źródła (`Phase3_qi2_rampoff.py` → `Phase3_qi2_output.txt`)

Do wygaszenia weszły **wszystkie** biegi SETTLED-* (lista + kotwice):
λ̃ ∈ {0.05, 0.06, 0.08, 0.10} + kontrola λ̃=0.01. Kontynuacja ze stanu
t=600, rampa 600→700, ewolucja swobodna do t=1700.

| λ̃ | rola | klasa po rampie | E_core(700) | E_core(1700) | E(1700)/E(700) | max_{r≤40}|ψ−1|(1700) | τ_obj |
|---|---|---|---|---|---|---|---|
| 0.05 | kotwica | **RETURN-TO-VACUUM** | +3.7307e−04 | +2.7856e−04 | 0.7467 | 2.02e−04 | 0.0 |
| 0.06 | lista | **RETURN-TO-VACUUM** | +3.8994e−04 | +2.9378e−04 | 0.7534 | 2.06e−04 | 0.0 |
| 0.08 | lista | **RETURN-TO-VACUUM** | +4.0737e−04 | +3.1107e−04 | 0.7636 | 2.10e−04 | 0.0 |
| 0.10 | lista | **RETURN-TO-VACUUM** | +4.1061e−04 | +3.1698e−04 | 0.7720 | 2.10e−04 | 0.0 |
| 0.01 | kontrola | **RETURN-TO-VACUUM** | +2.0759e−04 | +1.4174e−04 | 0.6828 | 1.38e−04 | 0.0 |

**Przebieg wygaszania (deskryptywnie; max_{r≤40}|ψ−1| i E_core):**

| λ̃ | t=600 | t=700 | t=1000 | t=1328.3 (koniec W_P) | t=1700 |
|---|---|---|---|---|---|
| 0.10 | 0.227 / −9.148 | 9.8e−6 / +4.1e−4 | 7.1e−5 / +2.2e−4 | 1.8e−4 / +3.2e−4 | 2.1e−4 / +3.2e−4 |
| 0.05 | 0.134 / −2.836 | 9.6e−6 / +3.7e−4 | 2.8e−5 / +1.9e−4 | 1.4e−4 / +2.8e−4 | 2.0e−4 / +2.8e−4 |

**Klasyfikacja zadziałała klauzulą amplitudową**
(max_{r≤40}|ψ−1| < 1e−3 w oknie końcowym W_F=[1600,1700]), nie
energetyczną: rampa jest dla tego układu tak adiabatyczna, że
**E_ref^off = E_core(700) jest już wielkością szumową** (≈4e−4, pięć
rzędów poniżej |E_core(600)|), więc iloraz E(1700)/E(700)≈0.75 nie
niesie informacji. Odnotowane jako fakt metodologiczny (NEEDS N3).

**Bilans energii (deskryptywnie; podczas rampy energia NIE jest
zachowana — praca źródła, LOCK: nie bramkować):**

| λ̃ | E_core(600) | E_core(700) | ΔE rampy (praca źródła) |
|---|---|---|---|
| 0.01 | −1.47591e−01 | +2.07592e−04 | +1.47798e−01 |
| 0.05 | −2.83625e+00 | +3.73070e−04 | +2.83662e+00 |
| 0.06 | −3.88556e+00 | +3.89936e−04 | +3.88595e+00 |
| 0.08 | −6.32327e+00 | +4.07374e−04 | +6.32368e+00 |
| 0.10 | −9.14841e+00 | +4.10611e−04 | +9.14882e+00 |

Czyli: wygaszanie źródła **oddaje polu dokładnie tyle energii, ile
studnia materii trzymała** — pole wraca do próżni bez pozostałości,
zero histerezy. (Deskryptywnie przy λ̃=0.10: E_core=−9.148 to suma
członu materii i **dodatniego** kosztu deformacji pola
E_core^field=+5.842 — stan osiadły jest wyłącznie równowagą
wymuszoną.)

**Potwierdzenia (LOCK §3):**
- PERSISTENT-OBJECT: **brak** ⟹ potwierdzenia h/2+dt/2
  bezprzedmiotowe.
- Obowiązkowa kontrola negatywu (najniższe λ̃ wśród SETTLED = 0.05),
  bieg pełny 0→1700 na h=0.025: klasa fazy włączonej SETTLED-DEF
  (zgodna), klasa po wygaszeniu **RETURN-TO-VACUUM — ZGODNE**.
- Kontrola czystości maszynerii wygaszania (λ̃=0.01):
  RETURN-TO-VACUUM, max_{r≤40}|ψ−1|(1700)=**1.382e−04 < 1e−3** —
  **PASS**.

⟹ **Q-I2-FAIL wg litery** (wszystkie biegi SETTLED-* dają
RETURN-TO-VACUUM zbieżnie).

**Konfrontacja z predykcją pre-rejestrowaną (bez reinterpretacji):**
LOCK §0 przewidywał „RETURN-TO-VACUUM lub COLLAPSE" — **predykcja
TRAFIONA**, 4/4 jako RETURN-TO-VACUUM. Kreacji w konwencji
kanonicznej NIE ma; stan podprogowy indukowany przez materię jest
**cieniem źródła**, nie obiektem.

## 6. Probe deskryptywny (poza protokołem werdyktowym)

`Phase3_desc_gridconv.py` → `Phase3_desc_gridconv_output.txt`,
uruchomiony **po** zapisaniu werdyktu Q-I1, wyłącznie jako materiał
do NEEDS (werdykt NIEZMIENIONY):

| bieg | klasa | ψ̄(0) | odniesienie h=0.05, dt=0.005 |
|---|---|---|---|
| λ̃=0.06, h=0.025 | SETTLED-DEF | 0.845711 | SETTLED-DEF 0.845733 (Δ=2.2e−5) |
| λ̃=0.08, h=0.025 | **SETTLED-SUB** | 0.808004 | SETTLED-SUB 0.808031 (Δ=2.7e−5) |
| λ̃=0.10, h=0.05, dt/2 | SETTLED-SUB | 0.772903 | SETTLED-SUB 0.772903 (Δ=1.6e−10) |

Odczyt deskryptywny: rozjazd **nie jest** efektem kroku czasowego
(dt/2 przy λ̃=0.10 daje tę samą kategorię i tę samą ψ̄(0)) i **nie
obejmuje** wnętrza okna (λ̃=0.06 i 0.08 są zbieżne siatkowo do 3e−5).
Rozjazd dotyczy **wyłącznie** biegu najbliższego progowi
(λ̃=0.10 vs λ̃_crit(h=0.05)=0.1077, tj. 7% poniżej progu): siatka
h=0.025 rozdziela ostrzejszy pik centralny w pierwszym overshoocie
(min ψ = 0.427 przy h=0.05) i tam ucieka z dziedziny. Innymi słowy:
**λ̃_crit jest zależne od siatki** — h=0.025 przesuwa próg poniżej
0.10, okno podprogowe się **zawęża, ale nie znika** (λ̃=0.08
pozostaje SETTLED-SUB na obu siatkach).

**To NIE zmienia werdyktu Q-I1.** LOCK §3 wyznaczył potwierdzenie
h=0.025 dla najgłębszego SETTLED-SUB; wskazanie po fakcie innego,
wygodniejszego przypadku jako „właściwego potwierdzenia" byłoby
forbidden move (anti-Lakatos). Wniosek idzie do NEEDS jako re-lock
metodologiczny (N1).

## 7. Higiena, incydenty, artefakty

- **Zero correction notes** — żaden gate nie wymagał diagnostyki;
  estymator dryfu sekularnego (correction note 1 poprzednika) został
  **przejęty i zamrożony w MD §5 PRZED pierwszym biegiem** tego cyklu
  (nie jest korektą tego cyklu; flagowany [INPUT-MD]).
- **Zero zmian** kryteriów / progów / detektorów / listy λ̃ / rampy /
  okien / sponge po pierwszym biegu produkcyjnym. Jedyny dopisek
  deskryptywny (0D punkt zwrotu, kotwice tabeli P1-I2) został
  wykonany **przed** Phase 3 i wyraźnie oflagowany w outputach.
- **Incydenty: brak.** Wszystkie biegi zakończone kodem 0; żadnego
  NonConvergence poza kategorią COLLAPSE (tam: 2 biegi BREAKDOWN,
  reszta pas graniczny).
- Integralność: `integrity_snapshot.txt` (SHA256 LOCKa, MD, silnika
  po FROZEN) zweryfikowany przy zamknięciu — **hashe niezmienione**.
- Artefakty: `Phase3_results/` — 20 biegów fazy włączonej
  (10 głównych + 6 bisekcji + 1 potwierdzenie h/2 + 3 probe) +
  5 biegów wygaszania + 1 potwierdzenie pełne 0→1700 (json + npz serii +
  checkpointy co 100 j.cz. + stany t=600), `qi1_summary.json`,
  `qi2_summary.json`, `verdict.json`.
- Budżet faktyczny: Phase 2 356 s; Q-I1 478 s; Q-I2 1008 s; probe
  326 s (zrównoleglenie po biegach; etapy ≪50 min).

## 8. Co z tego wynika (drzewo LOCK §5 — litera)

Zestaw werdyktów **Q-I1-INCONCLUSIVE ∧ Q-I2-FAIL** nie trafia
w żadną z trzech gałęzi „czystych" drzewa LOCK §5 (te wymagają
Q-I1-PASS albo Q-I1-FAIL); aktywowana jest gałąź ostatnia:
**„INCONCLUSIVE → NEEDS metodologiczny (t_on, Δ rampy, siatka λ̃)"**
— rozwinięta w [[NEEDS.md]] (N1).

Odczyt deskryptywny, **do decyzji autora** (nie werdykt):

1. **Q-I2-FAIL jest wynikiem mocnym i niezależnym od losu Q-I1.**
   Każdy stan osiadły, jaki udało się wytworzyć statyczną materią
   korpusową — również oba stany **podprogowe** (ψ̄(0)=0.773 i 0.808,
   detektor M911 przekroczony) — znika bez śladu po adiabatycznym
   wygaszeniu źródła. Hipoteza „materia statyczna kreuje obiekt"
   nie ma nośnika numerycznego: **materia deformuje albo niszczy,
   nie kreuje**.
2. Gałąź „kreacja z materii statycznej" jest tym samym **domknięta
   negatywnie w warstwie faktu** (predykcja pre-rejestrowana
   trafiona); formalne domknięcie gałęzi drzewa M911-N2 wymaga
   jednak user-gate, bo litera Q-I1 nie jest rozstrzygnięta.
3. Otwarte pozostaje wyłącznie to, co LOCK §5 kierował do osobnego
   LOCKa: **samouzgodnione ρ(ψ)** („lepton = stan związany ze
   źródłem") — tu nietykalne (S05, zakaz cyklu).

**Zakazy dotrzymane:** zero claimów o masach leptonów, zero claimów
o oscylonach, zero samouzgodnienia ρ(ψ), zero dopisków do `core/`,
`STATE.md` i gita.
