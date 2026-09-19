---
title: "Phase0_balance (LOCK) — op-matter-induced-creation: czy materia INDUKUJE kreację obiektu (stan podprogowy) i czy obiekt PRZEŻYWA wygaszenie źródła?"
date: 2026-09-15
type: phase0-balance
tgp_owner: research/op-matter-induced-creation-2026-09-15
status: LOCKED
anti_lakatos_lock: ACTIVE
related:
  - "[[HANDOFF_PROMPT.md]]"
  - "[[README.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# Phase 0 — LOCK (zero obliczeń przed zapisem tego pliku)

**Autoryzacja:** user-gate 2026-09-15 — wybór „Wątek kreacji: materia-indukcja
(M911-N1)" po zamknięciu `op-collapse-matter-source` (Q-H1-PULL).
Realizuje NEEDS N1 tamtego cyklu + gałąź „sprzężenie z materią" drzewa
hipotezy kreacji z [[../op-metric-pair-M911-2026-09-02/NEEDS.md]] (N2/Q-B-FAIL:
czysty sektor nie kreuje — czy kreuje materia?).

## 0. Kontekst i pytania

Poprzednik ustalił: λ̃≤0.05 → DEFORMATION (min ψ̄=0.866>5/6);
λ̃≥0.2 → COLLAPSE indukowany; **λ̃_crit∈(0.05,0.2]**. Deskryptywnie
zaobserwowano też osiadły stan PODPROGOWY ψ(0)→0.61<5/6 (przetrwaniec
Q-H2 przy λ̃=0.2 na tle startu gaussowskiego). Pytanie o KREACJĘ ma
dwie części, obie rozstrzygalne w oknie krytycznym:

- **Q-I1 (indukowany obiekt):** czy istnieje λ̃ w oknie (0.05,0.2),
  przy którym próżnia + źródło osiada w TRWAŁYM stanie podprogowym
  (ψ̄(0)<5/6 — detektor progowy M911) bez ucieczki z dziedziny?
  (= kreacja indukowana wg detektorów, ze źródłem WŁĄCZONYM.)
- **Q-I2 (prawdziwa kreacja — centralne):** czy JAKIKOLWIEK stan
  osiadły (podprogowy lub deformacja) PRZEŻYWA adiabatyczne wygaszenie
  źródła jako zlokalizowany obiekt? RETURN-TO-VACUUM = deformacja
  (obiekt jest cieniem źródła); PERSISTENT-OBJECT = kreacja.

**Predykcja pre-rejestrowana (z dotychczasowej wiedzy programu):**
w gałęzi zdrowej BEZ źródła nie ma statycznych solitonów (Q-B-FAIL)
ani oscylonów (Q-E/Q-G) ⟹ **oczekiwany wynik Q-I2: RETURN-TO-VACUUM
lub COLLAPSE** (pozytyw byłby PIERWSZĄ kreacją w konwencji kanonicznej
— dlatego wymaga najostrzejszych potwierdzeń). Pre-rejestrowany
pozytyw NIE jest oczekiwany; wynik zgodny z predykcją też jest
pełnoprawny (domyka gałąź „kreacja z materii statycznej").

**Zakres:** indukcja i trwałość obiektu. ZAKAZ claimów o masach
leptonów, oscylonach i o samouzgodnionym ρ(ψ) (osobne przyszłe LOCKi).

## 1. Model (formy FROZEN — CYTAT z zamkniętych poprzedników)

Jak w `op-collapse-matter-source` (P1-H1 wykonane — CYTAT, bez
ponownego wyprowadzania): M=ψ⁶/(4−3ψ)², 𝒦=ψ⁴, 𝒰=ψ⁴/4−ψ³/3,
**𝒰_mat=λ̃(t)·ρ̂(r)·ψ²/(4−3ψ)**, ρ̂(r)=exp(−r²/18); K_geo=γ=c₀=1;
π=Mψ̇; dziedzina (0,4/3); tożsamość energii bez kancelacji dziedziczona.
Silnik: kopia `../op-collapse-matter-source-2026-09-14/engine_core.py`
(z członem materii) — jedyna zmiana merytoryczna: **λ̃(t) zależne od
czasu** (profil wygaszania §3); zapis JAWNY w method_decisions.

## 2. Phase 1 — analityka (sympy; PRZED numeryką)

- P1-I1: gate cytatów form (1e−12, ψ∈{0.9,1,1.1}) + tożsamość
  ∂𝒰_mat/∂ψ=λ̃ρ̂ψ(8−3ψ)/(4−3ψ)².
- P1-I2 (pre-rejestracja ilościowa): statyczna odpowiedź nieliniowa
  0D (min. algebraiczne 𝒰+𝒰_mat po ψ przy ρ̂=1): krzywa ψ_min(λ̃)
  i λ̃_fold, przy którym lokalne minimum ZNIKA (fold/saddle-node) —
  przewidywany mechanizm λ̃_crit; tabela ψ_min dla λ̃∈{0.06…0.18}
  zapisana PRZED Phase 3 (konfrontacja deskryptywna z ψ̄(0) osiadłym).
- P1-I3 (fakt): po wygaszeniu (λ̃=0) jedyne znane stany trwałe to
  ψ≡1 (Q-B-FAIL, Q-A-PASS) — cytat, kontekst predykcji Q-I2.

## 3. Protokół numeryczny (FROZEN)

Siatka r_i=(i+½)h, **h=0.05** (potwierdzenia h=0.025), R=200, sponge
smootherstep γ₀=1 na [160,200], dt=0.005 (dt/2 w potwierdzeniach),
E_core r≤80; zapis ψ(0,t), min/max ψ, E_core co dt_out=0.1; start
zawsze: ψ≡1, π₀=0, źródło od t=0.

**Faza włączona (Q-I1):** λ̃∈{0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18}
(7 wartości, FROZEN) + kontrole {0.05, 0.20} (kotwice zgodności
z poprzednikiem); t_on=600. Klasyfikacja stanu w oknie [500,600]
(średnie po oknie): **SETTLED-SUB** (‖ψ̇‖∞→szum ORAZ ψ̄(0)<5/6),
**SETTLED-DEF** (osiadły, ψ̄(0)≥5/6), **COLLAPSE** (nadkategoria:
pas 1e−6/4/3−1e−6 lub niefinityczność, okno ≤1 j.cz.), **UNSETTLED**
(nie osiadł do t=600).
**Bisekcja λ̃_crit (FROZEN):** 6 kroków bisekcji między najwyższym λ̃
bez COLLAPSE a najniższym z COLLAPSE (z listy+kotwic); wynik
deskryptywny λ̃_crit±okno.

**Faza wygaszania (Q-I2):** dla KAŻDEGO biegu SETTLED-* : kontynuacja
z λ̃(t) = λ̃·S((t_off+Δ−t)/Δ) (smootherstep, Δ=100, t_off=600 →
λ̃=0 od t=700), dalej ewolucja swobodna do t_max=1700 (≥100 T₀ po
wygaszeniu + margines). Klasyfikacja po wygaszeniu (od t=700):
- **PERSISTENT-OBJECT:** E_core(t)≥0.5·E_core(700) nieprzerwanie
  przez ≥100 T₀ ORAZ max|ψ−1|(r≤40)≥0.02 w całym oknie (obiekt
  zlokalizowany, nie szum); wymaga potwierdzenia h=0.025 ORAZ dt/2
  (czas życia ±10%).
- **RETURN-TO-VACUUM:** E_core(1700)<0.05·E_core(700) lub
  max|ψ−1|(r≤40)<1e−3 w oknie końcowym.
- **COLLAPSE** (nadkategoria j.w.) / **INCONCLUSIVE-RUN** (inaczej
  lub kategoria niezbieżna).
**Potwierdzenia (FROZEN):** każdy PERSISTENT-OBJECT → h=0.025+dt/2;
obowiązkowo jeden RETURN-TO-VACUUM (najniższe λ̃ SETTLED) → h=0.025
(kontrola negatywu). Kontrola czystości wygaszania: bieg λ̃=0.01
(głęboka deformacja liniowa) z identycznym profilem wygaszenia —
oczekiwany RETURN-TO-VACUUM (gate maszynerii wygaszania, próg
max|ψ−1|<1e−3 na końcu).

**Phase 2 — bramka (FROZEN; FAIL ⟹ STOP):** P2a próżnia λ̃=0
‖ψ−1‖∞≤1e−10 (100 T₀); P2b regresja: λ̃=0.05 stała → osiadła
DEFORMATION z ψ̄(0) zgodnym z poprzednikiem (1−0.134018) do ±1%;
λ̃=0.5 stała → COLLAPSE t=0.375±5%; P2c dryf energii przy λ̃ STAŁYM
0.05 ≤1e−6/100 T₀ (energia z 𝒰_mat zachowana gdy λ̃=const; podczas
rampy NIE jest — odnotować, nie bramkować).

## 4. Werdykty (litera; INCONCLUSIVE ≠ pozytyw)

- **Q-I1-PASS:** ≥1 λ̃ z listy daje SETTLED-SUB zbieżnie (kategoria
  na h i h=0.025 dla najgłębszego przypadku).
- **Q-I1-FAIL:** wszystkie λ̃ z listy dają SETTLED-DEF lub COLLAPSE
  zbieżnie (okno krytyczne puste: deformacja przechodzi w kolaps bez
  stanu podprogowego).
- **Q-I1-INCONCLUSIVE:** inaczej.
- **Q-I2-PASS (KREACJA):** ≥1 PERSISTENT-OBJECT potwierdzony
  (h/2 i dt/2). (Wbrew pre-rejestrowanej predykcji — dopuszczalne,
  litera rozstrzyga.)
- **Q-I2-FAIL (zgodnie z predykcją):** wszystkie biegi SETTLED-* po
  wygaszeniu dają RETURN-TO-VACUUM lub COLLAPSE zbieżnie.
- **Q-I2-INCONCLUSIVE:** inaczej.

## 5. Drzewo decyzyjne (pre-rejestrowane)

- **Q-I2-PASS** → pierwsza kreacja w konwencji kanonicznej: eskalacja
  user-gate CORE (konfrontacja z Q-B-FAIL M911; charakteryzacja
  obiektu = osobny LOCK; ZAKAZ claimów o leptonach do tego czasu).
- **Q-I1-PASS ∧ Q-I2-FAIL** → obiekt indukowany istnieje, ale jest
  cieniem źródła: hipoteza „lepton = stan związany ze źródłem"
  wymaga samouzgodnienia ρ(ψ) — kandydat następnego LOCKa (user).
  Gałąź „kreacja z materii statycznej" DOMKNIĘTA negatywnie.
- **Q-I1-FAIL ∧ Q-I2-FAIL** → okno krytyczne puste i nic nie
  przeżywa: materia statyczna wyłącznie deformuje albo niszczy —
  wątek kreacji wraca do genezy poziomu 0 (Γ+s_i, M911-N1) — user.
- INCONCLUSIVE → NEEDS metodologiczny (t_on, Δ rampy, siatka λ̃).

## 6. Forbidden moves (egzekwowane)

Rdzeń `.tex`/STATE.md/git NIETYKANE; katalogi innych cykli tylko
odczyt; formy M,𝒦,𝒰,𝒰_mat bez modyfikacji (cytaty); lista λ̃,
profil rampy (smootherstep, Δ=100), progi, detektory, okna, sponge
niezmienialne po pierwszym biegu produkcyjnym; ZAKAZ podłóg/barier
poza pasem klasyfikacyjnym; ZAKAZ pól dynamicznych i samouzgodnienia
ρ(ψ) (S05; ρ̂ statyczne, λ̃(t) tylko wg zamrożonej rampy); ZAKAZ
claimów o masach leptonów/oscylonach; predykcja P1-I2 i predykcja
Q-I2 nienaruszalne (raport bez reinterpretacji); INCONCLUSIVE ≠
pozytyw; correction note tylko dla błędu implementacji, PRZED użyciem
wyniku, pierwotne outputy zachowane.

## 7. Deliverables

`Phase_method_decisions.md` (FROZEN) · `engine_core.py` (kopia+cytat
+λ̃(t)) · `Phase1_analytic.py`+output · `Phase2_gate.py`+output ·
`Phase3_qi1_settle.py`+output · `Phase3_qi2_rampoff.py`+output +
`Phase3_results/` (json+npz, verdict.json) · `Phase_FINAL_close.md` ·
`NEEDS.md` · dopis logu `README.md` (folder_status/verdict).
Cykl bez FINAL+NEEDS+README NIE jest zakończony.
