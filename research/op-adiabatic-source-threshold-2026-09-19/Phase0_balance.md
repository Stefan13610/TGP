---
title: "Phase0_balance (LOCK) — op-adiabatic-source-threshold: czy λ̃_crit jest własnością MODELU, czy PROTOKOŁU nagłego załączenia?"
date: 2026-09-19
type: phase0-balance
tgp_owner: research/op-adiabatic-source-threshold-2026-09-19
status: LOCKED
anti_lakatos_lock: ACTIVE
related:
  - "[[HANDOFF_PROMPT.md]]"
  - "[[README.md]]"
  - "[[../op-matter-induced-creation-2026-09-15/NEEDS.md]]"
  - "[[../op-matter-induced-creation-2026-09-15/Phase_FINAL_close.md]]"
  - "[[../op-collapse-matter-source-2026-09-14/Phase_FINAL_close.md]]"
---

# Phase 0 — LOCK (zero obliczeń przed zapisem tego pliku)

**Autoryzacja:** user-gate 2026-09-19 — wybór N4 z
[[../op-matter-induced-creation-2026-09-15/NEEDS.md]] („próg jest
DYNAMICZNY, nie statyczny") po zamknięciu Q-I1-INCONCLUSIVE + Q-I2-FAIL.
Cykl realizuje też **N1 tamtego cyklu** (reguła potwierdzenia siatkowego)
— wbudowaną tutaj jako poprawiony protokół, nie jako osobne pytanie.

## 0. Kontekst i pytania

Dwa cykle zbudowały mapę λ̃ **wyłącznie dla rodziny startów „źródło
włączone skokowo w t=0"** (oba LOCK-i §3: „start zawsze: ψ≡1, π₀=0,
źródło od t=0"). Zmierzony próg to λ̃_crit = 0.107656 ± 0.000156
(bisekcja 6 kroków, h=0.05). Pre-rejestrowana predykcja mechanizmu
— fold 0D, λ̃_fold = 0.285769734239 — **NIE trafiła: rozjazd 2.65×**.
Kolapsy zachodzą w PIERWSZYM overshoocie (t_end ∈ [0.60, 0.83]),
a lepszym deskryptorem jest punkt zwrotu 0D ze startu ψ≡1, ψ̇=0:
λ̃_turn = 1/12 = 0.08333.

Wniosek deskryptywny poprzednika: **reżim λ̃ ∈ (0.1077, 0.2858)
— gdzie statyczne minimum 0D NADAL ISTNIEJE, a nasz protokół daje
wyłącznie COLLAPSE — jest nietknięty.** Ten cykl go otwiera.

- **Q-J1 (centralne):** czy przy ADIABATYCZNYM załączaniu źródła
  (rampa w górę, ta sama rodzina smootherstep) istnieją stany osiadłe
  przy λ̃ **powyżej** progu nagłego załączenia (λ̃ > 0.1078)?
  (= czy λ̃_crit jest własnością protokołu, a nie modelu.)
- **Q-J2:** czy λ̃_crit(Δ_on) SATURUJE przy wydłużaniu rampy, czy
  dryfuje bez granicy? (= czy istnieje dobrze określona granica
  adiabatyczna, do której można się odwoływać.)
- **Q-J3 (warunkowe na Q-J1-PASS):** czy NAJGŁĘBSZY stan osiągnięty
  adiabatycznie — głębszy niż cokolwiek dotąd — przeżywa wygaszenie
  źródła? (= powtórka Q-I2 na głębszych stanach.)

**Predykcje pre-rejestrowane (dwie, rozbieżne co do znaku):**

1. **Q-J1: PASS oczekiwany** (uporządkowanie λ̃_crit^adiab >
   λ̃_crit^nagłe = 0.1077). Uzasadnienie: kolaps nagły jest
   przestrzeleniem energetycznym, a nie utratą równowagi; rampa
   odbiera nadwyżkę kinetyczną. **To jest pre-rejestracja POZYTYWNA
   — pierwsza w tej serii cykli.** Q-J1-FAIL byłby mocnym
   zaskoczeniem: oznaczałby, że próg jest odporny na protokół,
   czyli JEST własnością modelu mimo rozjazdu z foldem.
2. **Q-J3: RETURN-TO-VACUUM / COLLAPSE oczekiwane** (zgodnie
   z Q-I2-FAIL: zero histerezy, energia oddawana co do jednostki
   przy wygaszaniu). Pozytyw byłby pierwszą kreacją w konwencji
   kanonicznej ⟹ eskalacja user-gate CORE.

**Odniesienie ilościowe (pre-rejestrowane, NIE bramkujące):**
λ̃_fold(0D, ρ̂=1) = 0.285769734239 — cytat z P1-I2 poprzednika.
Hipoteza robocza: λ̃_crit(Δ_on→∞) → λ̃_fold. Zastrzeżenie zapisane
PRZED obliczeniami: fold policzono w 0D przy ρ̂=1, a bieg 3D ma
człony gradientowe i ρ̂(r) < 1 poza centrum — **zgodność co do
liczby NIE jest oczekiwana, oczekiwane jest uporządkowanie**
λ̃_crit(Δ_on) rosnące i ograniczone z góry przez λ̃_fold.

**Zakres:** zależność progu od protokołu załączania + trwałość
stanów głębokich. ZAKAZ claimów o masach leptonów, oscylonach
i o samouzgodnionym ρ(ψ) (osobne LOCK-i).

## 1. Model (formy FROZEN — CYTAT z zamkniętych poprzedników)

Bez zmian wobec `op-matter-induced-creation` (P1-I1/P1-H1 wykonane
— CYTAT, bez ponownego wyprowadzania):
M = ψ⁶/(4−3ψ)², 𝒦 = ψ⁴, 𝒰 = ψ⁴/4 − ψ³/3,
**𝒰_mat = λ̃(t)·ρ̂(r)·ψ²/(4−3ψ)**, ρ̂(r) = exp(−r²/18);
K_geo = γ = c₀ = 1; π = Mψ̇; dziedzina (0, 4/3); tożsamość energii
bez kancelacji dziedziczona.

Silnik: kopia `../op-matter-induced-creation-2026-09-15/engine_core.py`
(ma już λ̃(t)) — **jedyna zmiana merytoryczna: profil λ̃(t) z rampą
W GÓRĘ** (§3). Zapis JAWNY w `Phase_method_decisions.md`.

## 2. Phase 1 — analityka (sympy; PRZED numeryką)

- **P1-J1:** gate cytatów form (1e−12, ψ ∈ {0.9, 1, 1.1}) —
  identyczny z P1-I1 poprzednika; dodatkowo tożsamość ciągłości
  λ̃(t) i jej pochodnej na obu końcach rampy w górę (S′ = S″ = 0).
- **P1-J2 (cytat + kontrola):** krzywa ψ_min(λ̃) i λ̃_fold z minimum
  algebraicznego 𝒰 + λ̃ψ²/(4−3ψ) przy ρ̂ = 1 — **przeliczyć
  niezależnie i porównać z wartością poprzednika 0.285769734239
  (zgodność do 1e−9 = gate)**; tabela ψ_min dla całej listy λ̃ z §3
  zapisana PRZED Phase 3.
- **P1-J3 (pre-rejestracja kryterium adiabatyczności):** najkrótsza
  skala czasowa liniowej odpowiedzi wokół ψ = 1 to T₀ = 2π
  (m² = c_s² = 1, cytat Q-D1). Zapisać PRZED Phase 3:
  Δ_on ∈ {25, 100, 400} = {3.98, 15.9, 63.7} × T₀. Rampa jest
  adiabatyczna wtedy, gdy Δ_on ≫ T₀ — wszystkie trzy spełniają,
  z rozpiętością 16×, co czyni Q-J2 dobrze postawionym.
- **P1-J4 (fakt, cytat):** λ̃_crit^nagłe = 0.107656 ± 0.000156;
  λ̃_turn(0D) = 1/12; Q-I2-FAIL (4/4 RETURN-TO-VACUUM) — kontekst
  predykcji Q-J1 i Q-J3.

## 3. Protokół numeryczny (FROZEN)

**Siatka (dziedziczona bez zmian):** r_i = (i+½)h, **h = 0.05**
(potwierdzenia h = 0.025), R = 200, sponge smootherstep γ₀ = 1
na [160, 200], dt = 0.005 (dt/2 w potwierdzeniach), E_core r ≤ 80,
zapis ψ(0,t), min/max ψ, E_core co dt_out = 0.1, checkpoint `.npz`
co 100 j.cz. Start zawsze: **ψ ≡ 1, π₀ = 0**.

**Profil załączania (JEDYNA zmiana merytoryczna, FROZEN):**

> **λ̃(t) = λ̃ · S(t/Δ_on)**, S = smootherstep
> S(x) = 0 dla x ≤ 0; S(x) = 6x⁵ − 15x⁴ + 10x³ dla 0 < x < 1;
> S(x) = 1 dla x ≥ 1.
> Plateau: λ̃(t) = λ̃ dla t ≥ Δ_on.
> **t_on ≔ Δ_on + 600** (600 j.cz. plateau po rampie — tyle samo,
> ile miał poprzednik). Okno klasyfikacji: **[t_on − 100, t_on]**.

**Lista λ̃ (FROZEN, 8 wartości):**
λ̃ ∈ {0.10, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.30}
+ **kotwica regresyjna {0.05}** (musi odtworzyć ψ̄(0) = 0.865982
poprzednika — rampa nie może zmienić stanu osiadłego przy λ̃ głęboko
podprogowym) + **kotwica falsyfikacyjna {0.30}** (leży POWYŻEJ
λ̃_fold = 0.2858; jeśli osiądzie, hipoteza „fold ogranicza z góry"
jest obalona — odnotować, nie ratować).

**Siatka Δ_on (FROZEN, 3 wartości):** {25, 100, 400}.
- **Δ_on = 100 = bieg PODSTAWOWY:** pełna lista 8 λ̃ + kotwica 0.05.
- **Δ_on = 25 i Δ_on = 400:** wyłącznie bisekcja λ̃_crit (bez pełnej
  listy) — oszczędność budżetu, wystarcza do Q-J2.

**Bisekcja λ̃_crit(Δ_on) (FROZEN, dla każdego z 3 Δ_on):** 6 kroków
między najwyższym λ̃ bez COLLAPSE a najniższym z COLLAPSE; dla
Δ_on ∈ {25, 400} przedział startowy wyznaczyć biegiem zgrubnym po
{0.10, 0.15, 0.21, 0.27} (4 biegi sondujące, FROZEN). Wynik
deskryptywny λ̃_crit ± okno.

**Klasyfikatory (dziedziczone 1:1 z MD §6 poprzednika, FROZEN):**
ψ(0,t) ≔ ψ(h/2, t); ψ̄(r) ≔ średnia próbek profilowych z okna
klasyfikacji; ψ̄(0) ≔ ψ̄(h/2); V ≔ max po oknie max_{r≤80}|ψ̇|;
D ≔ max_{r≤80}|ψ̄(r)−1|; **osiadły ⟺ V ≤ 0.01·max(D, 1e−12)**.
Priorytety: **COLLAPSE** (ψ ≥ 4/3−1e−6 lub ψ ≤ 1e−6 lub
niefinityczność, w dowolnym t) > **SETTLED-SUB** (osiadły ∧
ψ̄(0) < 5/6) > **SETTLED-DEF** (osiadły ∧ ψ̄(0) ≥ 5/6) >
**UNSETTLED**.

**Reguła potwierdzenia siatkowego (POPRAWIONA — realizuje N1
poprzednika; FROZEN):** potwierdzeniu h = 0.025 podlegają **DWA**
biegi z Δ_on = 100:
  (a) **najgłębszy SETTLED-\*** z listy, ORAZ
  (b) **jeden SETTLED-\* odsunięty od progu**: największe λ̃ z listy
      spełniające λ̃ ≤ 0.8·λ̃_crit(Δ_on=100).
**Reguła rozjazdu (zapisana PRZED obliczeniami):** przy niezgodności
kategorii h vs h/2 decyduje **h = 0.025**, a bieg otrzymuje etykietę
`GRID-DIVERGENT` i **nie może** podpierać werdyktu PASS. Werdykt
PASS wymaga zbieżności biegu (b). Jeśli (b) nie istnieje (brak
SETTLED-\* poniżej 0.8·λ̃_crit), werdykt Q-J1 = INCONCLUSIVE
niezależnie od (a).

**Faza wygaszania (Q-J3, warunkowa na Q-J1-PASS):** dla
**najgłębszego zbieżnego SETTLED-\*** (Δ_on = 100) kontynuacja
z checkpointu t = t_on: λ̃(t) = λ̃·S((t_off + Δ − t)/Δ),
**Δ = 100, t_off = t_on** (profil zjazdu identyczny z poprzednikiem),
dalej ewolucja swobodna do **t_max = t_on + 1100**. Klasyfikacja
od t_off + Δ, progi **identyczne z Q-I2 poprzednika**:
- **PERSISTENT-OBJECT:** E_core(t) ≥ 0.5·E_core(t_off+Δ)
  nieprzerwanie przez ≥ 100 T₀ ORAZ max|ψ−1|(r ≤ 40) ≥ 0.02
  w całym oknie; wymaga potwierdzenia h = 0.025 ORAZ dt/2.
- **RETURN-TO-VACUUM:** max|ψ−1|(r ≤ 40) < 1e−3 w oknie końcowym
  [t_max−100, t_max]. **Uwaga (realizuje N3 poprzednika): człon
  energetyczny E(t_max) < 0.05·E(t_off+Δ) zostaje USUNIĘTY
  z kryterium** — u poprzednika był martwy (stosunek dwóch resztek
  = 0.75 przy progu 0.05). Raportować go deskryptywnie.
- **COLLAPSE** (nadkategoria j.w.) / **INCONCLUSIVE-RUN** (inaczej
  lub kategoria niezbieżna).

**Phase 2 — bramka (FROZEN; FAIL ⟹ STOP):**
- **P2a:** próżnia λ̃ = 0 z pełnym profilem rampy: ‖ψ−1‖∞ ≤ 1e−10
  przez 100 T₀ (rampa mnożona przez zero nie może nic wzbudzić).
- **P2b (regresja rampy w dół — dziedziczona):** λ̃ = 0.05 STAŁE
  (bez rampy) ⟹ ψ̄(0) = 0.865982 ± 1%; λ̃ = 0.5 STAŁE ⟹ COLLAPSE
  t = 0.375 ± 5%. Kotwice poprzednika co do cyfry.
- **P2c (NOWY gate rampy w górę):** λ̃ = 0.05 z rampą Δ_on = 100
  musi dać **ten sam stan osiadły** co λ̃ = 0.05 bez rampy:
  |ψ̄(0)_rampa − ψ̄(0)_nagłe| ≤ 1e−3. *To jest właściwy test tezy
  cyklu na przypadku kontrolnym: głęboko pod progiem protokół
  załączania NIE MOŻE mieć znaczenia. FAIL ⟹ STOP (rampa zmienia
  stan końcowy tam, gdzie nie powinna ⟹ błąd implementacji).*
- **P2d:** dryf energii przy λ̃ STAŁYM 0.05 ≤ 1e−6/100 T₀
  (estymator sekularny poprzednika, N6). Podczas rampy energia
  NIE jest zachowana (praca źródła) — odnotować, nie bramkować.

## 4. Werdykty (litera; INCONCLUSIVE ≠ pozytyw)

- **Q-J1-PASS:** ≥ 1 λ̃ **> 0.1078** (ściśle powyżej progu nagłego,
  górny koniec okna bisekcji poprzednika) daje SETTLED-\* przy
  Δ_on = 100, przy spełnionej regule potwierdzenia siatkowego §3
  (w tym zbieżność biegu (b)).
- **Q-J1-FAIL:** wszystkie λ̃ > 0.1078 z listy dają COLLAPSE lub
  UNSETTLED przy Δ_on = 100, zbieżnie. *(= próg odporny na protokół;
  wbrew pre-rejestrowanej predykcji.)*
- **Q-J1-INCONCLUSIVE:** inaczej (w tym: brak biegu (b)).
- **Q-J2-SATURATES:** λ̃_crit rośnie monotonicznie z Δ_on ORAZ
  I₂ ≤ 0.25·I₁, gdzie I₁ = λ̃_crit(100) − λ̃_crit(25),
  I₂ = λ̃_crit(400) − λ̃_crit(100), oba > okna bisekcji.
- **Q-J2-DRIFTS:** monotonicznie rosnące, ale I₂ > 0.25·I₁.
- **Q-J2-INCONCLUSIVE:** niemonotoniczne poza oknami bisekcji.
- **Q-J3-PASS (KREACJA):** ≥ 1 PERSISTENT-OBJECT potwierdzony
  (h/2 i dt/2). *(Wbrew predykcji — dopuszczalne, litera rozstrzyga.)*
- **Q-J3-FAIL (zgodnie z predykcją):** najgłębszy stan po wygaszeniu
  daje RETURN-TO-VACUUM lub COLLAPSE zbieżnie.
- **Q-J3-INCONCLUSIVE:** inaczej. **Q-J3-NIEURUCHOMIONE:** gdy
  Q-J1 ≠ PASS (warunkowość — nie jest brakiem wyniku).

## 5. Drzewo decyzyjne (pre-rejestrowane)

- **Q-J1-PASS ∧ Q-J2-SATURATES** → λ̃_crit^adiab jest dobrze
  określoną wielkością modelu; **cała mapa λ̃ obu poprzedników
  zostaje przeklasyfikowana jako własność protokołu** (kandydat
  dopisku do FINAL-i tamtych cykli — user-gate, bez zmiany ich
  werdyktów). Reżim (0.1077, λ̃_crit^adiab) staje się nowym
  obszarem roboczym dla stanów związanych ze źródłem.
- **Q-J1-PASS ∧ Q-J2-DRIFTS** → progu adiabatycznego NIE MA
  w zbadanym zakresie; λ̃_crit jest funkcją protokołu bez granicy
  — wniosek metodologiczny mocny (każdy próg w tym programie
  wymaga odtąd deklaracji protokołu). Kandydat dopisku meta.
- **Q-J1-FAIL** → próg jest własnością MODELU mimo rozjazdu 2.65×
  z foldem; mechanizmem jest coś innego niż overshoot — user-gate,
  kandydat na cykl o naturze granicy dziedziny (podejrzenie:
  metryczny defocusing c = (4−3ψ)/ψ, cytat dopisku core 2026-09-15).
- **Q-J3-PASS** → pierwsza kreacja w konwencji kanonicznej:
  eskalacja user-gate CORE (konfrontacja z Q-B-FAIL i Q-I2-FAIL);
  charakteryzacja obiektu = osobny LOCK; ZAKAZ claimów o leptonach.
- **Q-J3-FAIL** → brak histerezy potwierdzony na stanach GŁĘBSZYCH
  niż dotąd: gałąź „stan związany ze statycznym źródłem" domknięta
  mocniej; wątek kreacji → samouzgodnione ρ(ψ) (S05, user-gate)
  albo geneza Γ+s_i (M911-N1).
- INCONCLUSIVE → NEEDS metodologiczny (Δ_on, siatka λ̃, okno).

## 6. Forbidden moves (egzekwowane)

Rdzeń `.tex` / STATE.md / git NIETYKANE; katalogi innych cykli tylko
odczyt; formy M, 𝒦, 𝒰, 𝒰_mat bez modyfikacji (cytaty); **lista λ̃,
siatka Δ_on, profil rampy (smootherstep), t_on = Δ_on+600, progi,
detektory, okna, sponge, reguła potwierdzenia siatkowego i reguła
rozjazdu — niezmienialne po pierwszym biegu produkcyjnym**; ZAKAZ
podłóg/barier poza pasem klasyfikacyjnym; ZAKAZ pól dynamicznych
i samouzgodnienia ρ(ψ) (S05; ρ̂ statyczne, λ̃(t) wyłącznie wg
zamrożonej rampy); ZAKAZ claimów o masach leptonów i oscylonach;
**predykcje P1-J2, Q-J1 i Q-J3 nienaruszalne** (raport bez
reinterpretacji — w szczególności Q-J1-FAIL NIE wolno przedstawić
jako „częściowego potwierdzenia hipotezy protokołu");
INCONCLUSIVE ≠ pozytyw; correction note tylko dla błędu
implementacji, PRZED użyciem wyniku, pierwotne outputy zachowane;
rejestr WEJŚĆ flagowany [INPUT].

## 7. Deliverables

`Phase_method_decisions.md` (FROZEN) · `engine_core.py` (kopia +
cytat + rampa w górę) · `Phase1_analytic.py` + output ·
`Phase2_gate.py` + output · `Phase3_qj1_settle.py` + output ·
`Phase3_qj2_bisect.py` + output · `Phase3_qj3_rampoff.py` + output
(warunkowo) + `Phase3_results/` (json + npz, verdict.json) ·
`integrity_snapshot.txt` · `Phase_FINAL_close.md` · `NEEDS.md` ·
dopis logu `README.md` (folder_status / verdict).
**Cykl bez FINAL + NEEDS + README NIE jest zakończony.**
