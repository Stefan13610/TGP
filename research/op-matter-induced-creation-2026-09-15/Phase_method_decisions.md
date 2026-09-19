---
title: "Phase_method_decisions — decyzje metodyczne FROZEN cyklu op-matter-induced-creation (formy CYTAT, 𝒰_mat CYTAT z P1-H1 poprzednika, JAWNA rampa λ̃(t), klasyfikatory fazy włączonej SETTLED-SUB/DEF, klasyfikatory po wygaszeniu PERSISTENT-OBJECT/RETURN-TO-VACUUM, bisekcja λ̃_crit, okna pomiarowe, gate'y P2)"
date: 2026-09-15
type: phase-method-decisions
tgp_owner: research/op-matter-induced-creation-2026-09-15
status: FROZEN
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase_method_decisions.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase1_output.txt]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase3_qh1_output.txt]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase_correction_note_1_p2c_estimator.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# Phase_method_decisions (FROZEN przed jakimkolwiek kodem i biegiem cyklu)

**Status: FROZEN 2026-09-15, ZERO obliczeń cyklu przed zapisem tego
dokumentu** (również zero benchmarków silnika). Wszystkie decyzje
poniżej są zamrożone; zmiana po pierwszym biegu produkcyjnym =
forbidden move (LOCK §6). Pozycje wykraczające poza literę LOCKa
oznaczone **[INPUT-MD]**; wejścia zamrożone LOCKiem — **[LOCK]**;
wielkości przejęte z zamkniętych cykli — **[CYTAT]**.

---

## 1. Formy sektora pola (CYTATY — bez modyfikacji; LOCK §1)

Źródło: `../op-collapse-matter-source-2026-09-14/Phase_method_decisions.md`
§1 (za `../op-r3-stationary-states-2026-09-14/Phase1_output.txt`,
Q-D1-PASS, konwencja |g^tt|):

> „M(psi) = psi^6/(4-3psi)^2 … M' = 12 psi^5 (2-psi)/(4-3psi)^3
>  K(psi) = psi^4 … K' = 4 psi^3
>  U(psi) = psi^4/4 - psi^3/3 … U' = psi^2(psi-1)
>  K_geo = gamma = c0 = 1 [LOCK]. ZERO podlog/barier."

⟹ **M(ψ)=ψ⁶/(4−3ψ)²**, **M′=12ψ⁵(2−ψ)/(4−3ψ)³**, **𝒦(ψ)=ψ⁴**,
**𝒦′=4ψ³**, **𝒰(ψ)=ψ⁴/4−ψ³/3**, **𝒰′=ψ²(ψ−1)**; K_geo=γ=c₀=1;
dziedzina ψ∈(0,4/3); pas klasyfikacyjny ψ≥4/3−1e−6 (górny) /
ψ≤1e−6 (dolny); **ZERO podłóg i barier poza pasem**.

Tożsamość bez kancelacji do ewaluatora energii (correction note 2
cyklu `op-r3-stationary-states`, dziedziczona):
**𝒰(ψ)−𝒰(1) = (ψ−1)²(3ψ²+2ψ+1)/12**.

EOM (LOCK §1): Mψ̈ + ½M′ψ̇² = (1/r²)(r²𝒦ψ′)′ − ½𝒦′ψ′² − 𝒰′
− ∂𝒰_mat/∂ψ; π ≔ Mψ̇.

## 2. Człon materii — CYTAT wyniku P1-H1 poprzednika (bez ponownego wyprowadzania)

Źródło: `../op-collapse-matter-source-2026-09-14/Phase1_output.txt`
(wszystkie gate'y PASS: simplify=0, kontrola numeryczna ≤3.55e−14):

> „U_mat(psi,r) = lam * rho_hat(r) * psi^2/(4-3psi)   [WYPROWADZONE]
>  dU_mat/dpsi  = lam * rho_hat(r) * psi(8-3psi)/(4-3psi)^2
>  U_mat - U_mat(1) = lam * rho_hat * (psi-1)(psi+4)/(4-3psi)
>  dpsi_lin(r) = -(5 lam/r) int rp rho_hat [e^{-|r-rp|}-e^{-(r+rp)}]/2 drp ;
>  dpsi_lin(0) = -5 lam I0, I0=0.776061935"

⟹ **𝒰_mat(ψ,r) = λ̃(t)·ρ̂(r)·ψ²/(4−3ψ)**,
**∂𝒰_mat/∂ψ = λ̃(t)ρ̂(r)·ψ(8−3ψ)/(4−3ψ)²**,
**𝒰_mat−𝒰_mat(1) = λ̃(t)ρ̂(r)·(ψ−1)(ψ+4)/(4−3ψ)**;
**ρ̂(r)=exp(−r²/18)** [LOCK, FROZEN], statyczne (S05: ρ = źródło
zewnętrzne, ZERO pól dynamicznych, ZAKAZ ρ(ψ)).
λ̃ ≡ (q c₀/Φ₀)ρ₀ ≥ 0. Formy FROZEN — ponowne wyprowadzanie NIE jest
wykonywane; P1-I1 sprawdza wyłącznie **gate cytatów** (§4).

## 3. Rampa λ̃(t) — JEDYNA zmiana merytoryczna silnika (LOCK §1, §3)

Zapis JAWNY (LOCK §3, dosłownie):

> **λ̃(t) = λ̃ · S( (t_off + Δ − t)/Δ )**, S = smootherstep
> S(x) = 0 dla x≤0, S(x) = 6x⁵ − 15x⁴ + 10x³ dla 0<x<1, S(x) = 1
> dla x≥1; **Δ = 100**, **t_off = 600**.

Konsekwencje (arytmetyka zamrożona, nie decyzja):
- t ≤ 600 ⟹ (700−t)/100 ≥ 1 ⟹ **λ̃(t) = λ̃** (faza włączona, stała);
- 600 < t < 700 ⟹ zjazd smootherstep (S′=S″=0 na obu końcach);
- **t ≥ 700 ⟹ λ̃(t) = 0 DOKŁADNIE** (ewolucja swobodna).

Implementacja [INPUT-MD — operacjonalizacja, bez wpływu na definicję]:
w silniku `self.lrho = λ̃(t)·ρ̂(r_i)`; w kroku Störmera–Verleta siła
F₁ liczona z λ̃(t), siła F₂ z λ̃(t+dt) (naturalne uogólnienie schematu
symetrycznego na siłę zależną od czasu; dla λ̃=const krok jest
**bitowo identyczny** z silnikiem poprzednika — podstawiana jest ta
sama tablica). Ewaluator energii używa λ̃(t) bieżącego.
Podczas rampy energia NIE jest zachowana (praca źródła) — raport
deskryptywny, **bez bramkowania** (LOCK §3 / handoff).

## 4. Phase 1 — analityka pre-rejestrowana (sympy/numeryka; PRZED Phase 3)

- **P1-I1 (gate cytatów, próg 1e−12):** dla ψ∈{0.9, 1, 1.1} [LOCK §2]
  i λ̃ρ̂=1: |wartość sympy(binarna) − wartość float silnika| ≤
  1e−12·max(1,|v|) dla M, M′, 𝒦, 𝒦′, 𝒰, 𝒰′, 𝒰−𝒰(1), 𝒰_mat,
  ∂𝒰_mat/∂ψ, 𝒰_mat−𝒰_mat(1); dodatkowo tożsamości symboliczne
  simplify(·)=0: ∂𝒰_mat/∂ψ − λ̃ρ̂ψ(8−3ψ)/(4−3ψ)²,
  𝒰_mat−𝒰_mat(1) − λ̃ρ̂(ψ−1)(ψ+4)/(4−3ψ), 𝒰−𝒰(1) −
  (ψ−1)²(3ψ²+2ψ+1)/12, M′−dM/dψ, 𝒦′−d𝒦/dψ, 𝒰′−d𝒰/dψ.
  Dodatkowo (LOCK §2) gate rampy: λ̃(t) w t∈{0, 600, 650, 700, 1700}
  = {λ̃, λ̃, λ̃/2, 0, 0} do 1e−15 oraz S′(0)=S′(1)=0.
- **P1-I2 (PRE-REJESTRACJA ILOŚCIOWA; LOCK §2):** statyczna odpowiedź
  nieliniowa 0D: f(ψ) ≔ 𝒰(ψ) + λ̃·ψ²/(4−3ψ) przy **ρ̂=1** (rdzeń
  źródła); **ψ_min(λ̃)** = położenie lokalnego minimum f w (0,4/3)
  spełniające f′=0, f″>0, ψ_min<1 (gałąź ciągła od ψ_min(0)=1);
  **λ̃_fold** = wartość λ̃, przy której lokalne minimum ZNIKA
  (saddle-node: f′=f″=0) — wyznaczana symbolicznie (rugownik/solve
  sympy) i numerycznie (brentq/maksimum gałęzi λ̃(ψ) z f′=0).
  **TABELA ψ_min(λ̃) dla λ̃∈{0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18}
  (lista LOCKa) — zapisana w `Phase1_output.txt` PRZED Phase 3**;
  deskryptywnie dopisane kotwice λ̃∈{0.01, 0.05, 0.20} [INPUT-MD —
  wyłącznie rozszerzenie deskryptywne tabeli, bez zmiany listy].
  Predykcja pre-rejestrowana (LOCK §2): λ̃_fold jest przewidywanym
  MECHANIZMEM λ̃_crit; konfrontacja z ψ̄(0) osiadłym = **deskryptywna**
  (nie bramkuje żadnego werdyktu; LOCK §6: predykcja P1-I2
  nienaruszalna — raportowana bez reinterpretacji).
- **P1-I3 (fakt, cytat kontekstu):** bez źródła (λ̃=0) jedynym znanym
  stanem trwałym gałęzi zdrowej jest ψ≡1 — cytat z
  `../op-metric-pair-M911-2026-09-02/NEEDS.md` (Q-A-PASS + Q-B-FAIL)
  oraz z LOCK §0 (brak oscylonów: Q-E/Q-G). Kontekst predykcji Q-I2.

## 5. Phase 2 — bramka (FROZEN; FAIL ⟹ STOP; LOCK §3)

Konfiguracja produkcyjna: h=0.05, R=200, dt=0.005, sponge ON
(smootherstep γ₀=1.0 na [160,200]), start ψ≡1, π₀=0, źródło od t=0.

- **P2a (próżnia):** λ̃=0, ψ≡1, π≡0, 100·T₀ (T₀=2π), **obie siatki**
  h∈{0.05, 0.025} [INPUT-MD — obie siatki, dziedziczone];
  gate **‖ψ−1‖∞ ≤ 1e−10** w całym biegu.
- **P2b (regresje, kotwice z `Phase3_qh1_output.txt` poprzednika):**
  1. λ̃=0.05 stały, h=0.05, t_max=300, okno [250,300]: gate
     **ψ̄(0) = 1 − 0.134018 = 0.865982 ± 1%** [LOCK];
  2. λ̃=0.5 stały, h=0.05: gate **COLLAPSE z t_end = 0.375 ± 5%**
     [LOCK]; bieg prowadzony do t=30 lub do zdarzenia.
  Okno/próbkowanie P2b-1 identyczne z §6 (profil co 5 j.cz. ⟹ 11
  próbek w [250,300]) — operacjonalizacja dziedziczona [INPUT-MD].
- **P2c (dryf energii przy λ̃ STAŁYM):** λ̃ = const = 0.05, start
  gauss a=+0.05 σ=3 (konfiguracja P2c poprzednika, [INPUT-MD
  dziedziczone]), **sponge OFF** (układ zamknięty) [INPUT-MD
  dziedziczone], t=700, **obie siatki**; estymator dryfu SEKULARNEGO
  przejęty z `Phase_correction_note_1_p2c_estimator.md` poprzednika
  (**przejęcie ZAMROŻONE TU, przed pierwszym biegiem** — nie jest
  korektą tego cyklu) [INPUT-MD]:
  **dryf ≔ (100/80)·|⟨E⟩_[90T₀,100T₀] − ⟨E⟩_[10T₀,20T₀]| /
  |⟨E⟩_[0,10T₀]| ≤ 1e−6 / 100T₀**;
  deskryptywnie: stary estymator offsetu ΔC2, fit LSQ na [10T₀,100T₀],
  max|E−E(0)|.

**FAIL któregokolwiek z P2a/P2b/P2c ⟹ STOP cyklu** (LOCK §3).

## 6. Phase 3 / Q-I1 — faza włączona: okna, detektory, klasyfikatory (FROZEN)

Biegi (LOCK §3): start ψ≡1, π₀=0, źródło od t=0, λ̃ STAŁE do t=600;
**lista λ̃ = {0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18}** (7, FROZEN)
**+ kotwice {0.05, 0.20}**; h=0.05, dt=0.005, t_on=600 (bieg fazy
włączonej prowadzony do t=600 = koniec plateau λ̃(t)).
Dodatkowo (LOCK §3, kontrola czystości wygaszania Q-I2): bieg
**λ̃=0.01** w identycznym protokole — NIE wchodzi do listy Q-I1
i NIE uczestniczy w werdykcie Q-I1 (wyłącznie kontrola maszynerii).

Definicje pomiarowe [INPUT-MD — operacjonalizacje, dziedziczone z MD
§6 poprzednika]:
- **ψ(0,t) ≔ ψ(r=h/2, t)**; próbkowanie skalarów co **dt_out=0.1**;
  pełny profil co **5 j.cz.**; checkpoint `.npz` co 100 j.cz.
- **ψ̄(r) ≔ średnia próbek profilowych z okna t∈[500,600]** (21
  próbek); **ψ̄(0) ≔ ψ̄(h/2)**; δψ(0) ≔ ψ̄(0)−1.
- **Osiadłość:** V ≔ max_{t∈[500,600]} max_{r≤80} |ψ̇(r,t)|
  (ψ̇=π/M, próbki co dt_out); D ≔ max_{r≤80}|ψ̄(r)−1|;
  **osiadły ⟺ V ≤ 0.01·max(D, 1e−12)** („‖ψ̇‖∞ → poziom szumu").
- **COLLAPSE** (nadkategoria, priorytet 1): ψ ≥ 4/3−1e−6 lub
  ψ ≤ 1e−6 w dowolnym węźle, albo niefinityczność/NonConvergence,
  w dowolnym t; zdarzenie rejestrowane w kroku wystąpienia, t_end
  z rozdzielczością dt=0.005 ≪ 1 j.cz. (spełnia „okno ≤1 j.cz."
  LOCKa); podtyp (sufit/podłoga/BREAKDOWN) raportowany deskryptywnie.
- **SETTLED-SUB** (priorytet 2): osiadły ORAZ **ψ̄(0) < 5/6**.
- **SETTLED-DEF** (priorytet 3): osiadły ORAZ **ψ̄(0) ≥ 5/6**.
- **UNSETTLED** (priorytet 4): nie osiadł do t=600 (V > 0.01·max(D,1e−12))
  bez zdarzenia brzegowego.
- **Potwierdzenie siatkowe (LOCK §3):** h=0.025 dla **najgłębszego
  SETTLED-SUB** (najmniejsze ψ̄(0)); jeśli SETTLED-SUB nie wystąpi —
  dla **najgłębszego SETTLED-DEF**. **„Zbieżnie" [INPUT-MD,
  dziedziczone z MD §6 poprzednika]** ⟺ kategoria na h=0.025 równa
  kategorii na h=0.05 tam, gdzie LOCK przewiduje drugą siatkę;
  dla pozostałych biegów kategoria = kategoria z h=0.05 (LOCK nie
  przewiduje tam drugiej siatki).
- **Bisekcja λ̃_crit (LOCK §3, FROZEN):** przedział startowy
  [λ̃_lo, λ̃_hi] = [najwyższe λ̃ z listy+kotwic BEZ COLLAPSE,
  najniższe λ̃ z listy+kotwic Z COLLAPSE]; **dokładnie 6 kroków**
  bisekcji; w każdym kroku bieg w **identycznym protokole**
  (h=0.05, dt=0.005, t_max=600, ten sam detektor COLLAPSE) w punkcie
  środkowym; COLLAPSE ⟹ λ̃_hi ← środek, brak COLLAPSE ⟹ λ̃_lo ←
  środek. Wynik **deskryptywny**: λ̃_crit = (λ̃_lo+λ̃_hi)/2 ±
  (λ̃_hi−λ̃_lo)/2. Protokół biegów bisekcji (t_max=600, ta sama
  siatka i detektor) [INPUT-MD — LOCK mówi „6 kroków", nie podaje
  t_max; wybrano protokół IDENTYCZNY z biegami głównymi, aby nie
  wprowadzać nowego kryterium].
- Deskryptywnie obowiązkowo: ψ̄(0) zmierzone vs tabela P1-I2
  (ψ_min(λ̃), 0D) — **bez bramkowania**.
- **E_core** ≔ energia próżniowo odjęta (czynnik 4π) z **gęstością
  PEŁNĄ** (z członem 𝒰_mat, próżniowo odjętym) dla r_i ≤ 80
  (człony gradientowe: r_{i+½} ≤ 80); deskryptywnie E_core^field
  (bez członu materii) [INPUT-MD, dziedziczone].

**Werdykty Q-I1 (litera LOCK §4, przywołanie):** Q-I1-PASS ⟺ ≥1 λ̃
**z listy** daje SETTLED-SUB zbieżnie; Q-I1-FAIL ⟺ wszystkie λ̃
z listy dają SETTLED-DEF lub COLLAPSE zbieżnie; inaczej
Q-I1-INCONCLUSIVE. (Kotwice 0.05/0.20 i bieg 0.01 są kontrolami,
nie wchodzą do kwantyfikatora „wszystkie λ̃ z listy".)

## 7. Phase 3 / Q-I2 — faza wygaszania: okna i klasyfikatory (FROZEN)

Biegi: dla **KAŻDEGO** biegu SETTLED-* z §6 (lista + kotwice)
kontynuacja z checkpointu t=600 (stan ψ,π zapisany na dysk),
rampa §3 (600→700), dalej ewolucja swobodna do **t_max=1700**;
h=0.05, dt=0.005, sponge ON. Plus kontrola czystości λ̃=0.01.

- **E_ref^off ≔ E_core(700)** (λ̃(700)=0 dokładnie ⟹ gęstość pełna =
  gęstość pola). **Okno trwałości W_P ≔ [700, 700+100·T₀]**,
  T₀=2π ⟹ W_P = [700, 1328.3185…]; **okno końcowe W_F ≔ [1600,1700]**
  [INPUT-MD — „okno końcowe" nie jest liczbowo w LOCKu; wybrane
  100 j.cz. przed t_max, symetrycznie do pozostałych okien cyklu].
- **PERSISTENT-OBJECT** (priorytet 2, po COLLAPSE): E_core(t) ≥
  0.5·E_ref^off dla WSZYSTKICH próbek t∈W_P **ORAZ**
  max_{r≤40}|ψ(r,t)−1| ≥ 0.02 dla WSZYSTKICH próbek t∈W_P.
  (Gdy E_ref^off ≤ 0 — kategoria niezbieżna ⟹ INCONCLUSIVE-RUN
  [INPUT-MD — zabezpieczenie dziedziczone; przy λ̃(700)=0 gęstość
  pola jest dodatnio określona, więc przypadek nie jest oczekiwany.)
- **RETURN-TO-VACUUM** (priorytet 3): E_core(1700) < 0.05·E_ref^off
  **lub** max_{r≤40}|ψ−1| < 1e−3 dla wszystkich próbek t∈W_F.
- **COLLAPSE** (nadkategoria, priorytet 1): jak §6, w dowolnym
  t ∈ (600, 1700].
- **INCONCLUSIVE-RUN** (priorytet 4): każdy inny przypadek.
- **τ_obj (czas życia do potwierdzeń ±10%)** [INPUT-MD,
  dziedziczone]: τ_obj ≔ największe T takie, że OBA warunki
  PERSISTENT-OBJECT zachodzą nieprzerwanie na [700, 700+T];
  τ_obj = 1000 ⟹ cenzurowane („≥1000"); zgodność potwierdzenia ⟺
  ta sama kategoria ORAZ (dla PERSISTENT-OBJECT) τ_obj ±10%
  (cenzurowane: zgodność ⟺ drugie ≥900; obie cenzurowane ⟹ zgodne);
  dla COLLAPSE t_end ±10%.
- **Potwierdzenia (LOCK §3, FROZEN):** każdy PERSISTENT-OBJECT →
  bieg h=0.025 (dt=0.005) ORAZ bieg dt=0.0025 (h=0.05); obowiązkowo
  jeden RETURN-TO-VACUUM — **bieg o najniższym λ̃ wśród SETTLED-***
  → h=0.025. Biegi potwierdzające prowadzone **od t=0** w tym samym
  protokole (checkpoint z h=0.05 nie przenosi się na inną siatkę)
  [INPUT-MD — operacjonalizacja].
- **Kontrola czystości wygaszania (LOCK §3):** λ̃=0.01, identyczny
  profil rampy; **oczekiwane RETURN-TO-VACUUM** oraz
  **max_{r≤40}|ψ−1| < 1e−3 na końcu** (próbka t=1700). Gate
  maszynerii wygaszania: FAIL ⟹ raportowany jako incydent
  metodologiczny i unieważnia interpretację pozytywu Q-I2
  (nie zmienia litery werdyktów pozostałych biegów).
- Bilans energii deskryptywnie: E_core(600), E_core(700),
  E_core(1700), E_core(1700)/E_core(700), praca źródła podczas rampy
  (ΔE w [600,700]) — **bez bramkowania** (LOCK §3).

**Werdykty Q-I2 (litera LOCK §4, przywołanie):** Q-I2-PASS ⟺ ≥1
PERSISTENT-OBJECT potwierdzony (h/2 i dt/2); Q-I2-FAIL ⟺ wszystkie
biegi SETTLED-* po wygaszeniu dają RETURN-TO-VACUUM lub COLLAPSE
zbieżnie; inaczej Q-I2-INCONCLUSIVE. **Predykcja pre-rejestrowana
(LOCK §0): RETURN-TO-VACUUM lub COLLAPSE** — konfrontacja bez
reinterpretacji; INCONCLUSIVE ≠ pozytyw.

## 8. Higiena, artefakty, budżet

- Silnik: `engine_core.py` = **kopia** pliku poprzednika z cytatem
  w nagłówku; jedyna zmiana merytoryczna = λ̃(t) (§3). Katalogi
  innych cykli **tylko do odczytu**.
- Artefakty biegów: `Phase3_results/*.json` (skalary, klasyfikacja),
  `Phase3_results/*.npz` (checkpointy stanu co 100 j.cz. + stan
  t=600 + serie czasowe), `Phase3_results/verdict.json`.
- `integrity_snapshot.txt`: SHA256 LOCKa, MD i silnika po FROZEN;
  weryfikacja przy zamknięciu.
- Batch w tle + AKTYWNE czekanie (pętla sprawdzająca marker w pliku),
  etapy ≤50 min; zrównoleglenie po niezależnych biegach [INPUT-MD].
- Outputy tekstowe: `Phase1_output.txt`, `Phase2_output.txt`,
  `Phase3_qi1_output.txt`, `Phase3_qi2_output.txt`.

## 9. Rejestr WEJŚĆ (flagowany)

**[LOCK]:** formy M,𝒦,𝒰 i 𝒰_mat (cytaty §1–§2); ρ̂=exp(−r²/18);
rampa λ̃(t)=λ̃·S((t_off+Δ−t)/Δ), Δ=100, t_off=600; lista λ̃
{0.06,0.08,0.10,0.12,0.14,0.16,0.18} + kotwice {0.05,0.20} + kontrola
0.01; h∈{0.05,0.025}; dt=0.005 (dt/2 w potwierdzeniach); R=200;
sponge smootherstep γ₀=1 na [160,200]; E_core r≤80; t_on=600;
t_max=1700; okno [500,600]; klasyfikacja po wygaszeniu od t=700;
progi 5/6, 1e−6/4/3−1e−6, 0.5·E_core(700), 100 T₀, 0.02, 0.05·E_core(700),
1e−3; bisekcja 6 kroków; P2a 1e−10; P2b ψ̄(0)=0.865982±1% i
t_end=0.375±5%; P2c ≤1e−6/100T₀; potwierdzenia h/2 (+dt/2 dla
PERSISTENT-OBJECT), ±10%; start ψ≡1, π₀=0, źródło od t=0.
**[CYTAT]:** δψ_lin(0)=−5λ̃·I₀, I₀=0.776061934827; kotwica
δψ(0)=−0.134018 @λ̃=0.05; kotwica COLLAPSE t=0.375 @λ̃=0.5;
λ̃_crit∈(0.05,0.2] (kontekst); ψ≡1 jedyny stan trwały bez źródła.
**[INPUT-MD]:** γ₀=1.0; dt_out=0.1; profil co 5 j.cz.; checkpoint co
100 j.cz.; ψ(0,t)≔ψ(h/2,t); ψ̄ = średnia profili okna; osiadłość
V≤0.01·max(D,1e−12); E_core gęstość pełna + E_core^field
deskryptywnie; priorytety klasyfikacji; „zbieżnie" = kategoria na
siatkach przewidzianych LOCKiem; t_max=600 w biegach bisekcji;
okno końcowe W_F=[1600,1700]; τ_obj i reguła cenzurowania;
biegi potwierdzające od t=0; estymator dryfu sekularnego P2c
(przejęty z correction note 1 poprzednika, zamrożony tutaj PRZED
pierwszym biegiem); konfiguracja P2c (gauss a=+0.05 σ=3, sponge OFF,
t=700, obie siatki); rozszerzenie tabeli P1-I2 o kotwice deskryptywne;
tolerancja punktu stałego = stagnacja maszynowa (dziedziczona);
zrównoleglenie biegów.

**FROZEN. Jakakolwiek zmiana kryteriów/progów/detektorów/okien/listy
λ̃/rampy/sponge po pierwszym biegu produkcyjnym = forbidden move
(LOCK §6). Korekta wyłącznie dla udokumentowanego błędu implementacji:
`Phase_correction_note_*.md` PRZED użyciem poprawionego wyniku,
pierwotne outputy zachowane.**
