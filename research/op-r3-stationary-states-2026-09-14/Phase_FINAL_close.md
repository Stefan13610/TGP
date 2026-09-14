---
title: "Phase_FINAL_close — zamknięcie: Q-E-INCONCLUSIVE (zero oscylonów w zamrożonej klasie startów: 2×RADIATED zbieżnie [a=±0.15, σ=3: τ≈207–211 ≪ 628, rdzeń wypromieniowuje], 6×BREAKDOWN-BOUNDARY zbieżnie [t≈2.3–19.8: pole ucieka do granicy dziedziny — w tym OBA starty quasi-R3 sin(r)/r], 2×INCONCLUSIVE [a=+0.25: kategoria niezbieżna między siatkami]); Q-F NIEURUCHOMIONE (warunkowe na Q-E-PASS); P1a PASS (κ²=ω²−1, κ=1⟺ω²=2 zgodne z R3), P1b PREDYKCJA pre-rejestrowana ω₂=−139/24<0 (miękka; domena małych amplitud POZA zbadaną rodziną startów |a|≥0.15 — bez konfrontacji rozstrzygającej), P2 PASS 6/6 po 2 korektach implementacyjnych harnessu (kwadranty FFT; kancelacja w ewaluatorze energii)"
date: 2026-09-14
type: phase-final-close
tgp_owner: research/op-r3-stationary-states-2026-09-14
status: CLOSED
verdict: "Q-E: INCONCLUSIVE wg litery (LOCK §2 Phase 3: PASS wymaga ≥1 OSCILLON zbieżnego — jest 0; FAIL wymaga WSZYSTKIE RADIATED zbieżnie — są 2/10; reszta: 6 BREAKDOWN-BOUNDARY zbieżnych kategorii i czasów [h, h/2 i dt/2; czasy zdarzeń zgodne do ~1%: qR3−0.20: 4.74/4.74/4.74; a−0.30σ3: 6.64/6.635/6.635; a−0.15σ6: 19.83/19.82/19.82; qR3+0.20: 2.275/2.27/2.2675] + 2 niezbieżne kategorie przy a=+0.25 [h05: BREAKDOWN-BOUNDARY(-LOWER) vs h025: BREAKDOWN niefinityczny]). Deskryptywnie (obowiązkowe): OBA starty quasi-R3 (a=±0.2, sin(r)/r, σ_w=15) NIE przeżywają jako stany oscylacyjne — kolaps do granicy dziedziny w t≈2.3–4.7 zbieżnie; kształt R3 w dynamice 2. rzędu gałęzi zdrowej natychmiast opuszcza dziedzinę. RADIATED: a=±0.15, σ=3: E_core spada monotonicznie do 1.4–1.5% E_ref w t=1000, τ=207–211 ≪ próg 628; ψ(0,t) po t≈300 to ogon dyspersyjny ~1e−3. Q-F: nieuruchomione wg litery (tylko przy Q-E-PASS). P1b: ω₂ = M′(1)²/16 + M′(1)𝒰‴(1)/8 − M″(1)/8 − 5𝒰‴(1)²/48 + 𝒰⁗(1)/16 = −139/24 ≈ −5.79 < 0 (miękka nieliniowość; warunek istnienia oscylonu małej amplitudy SPEŁNIONY jako predykcja) — zamrożona rodzina startów ma |a|≥0.15 i nie sonduje reżimu małych amplitud, więc wynik Phase 3 NIE konfrontuje predykcji w jej naturalnej domenie (raport bez reinterpretacji, forbidden move h dotrzymany). ZAKAZ claimów o masach leptonów dotrzymany; INCONCLUSIVE ≠ pozytyw."
anti_lakatos_lock: PRESERVED
tags: [r3-stationary-states, oscillon-search, second-order-dynamics, healthy-branch, breakdown-boundary, radiated, inconclusive-qe, lindstedt-poincare, soft-nonlinearity, closed]
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase_correction_note_P2_gates.md]]"
  - "[[Phase_correction_note_2_energy_eval.md]]"
  - "[[NEEDS.md]]"
  - "[[README.md]]"
  - "[[../op-action-audit-spectrum-insert-2026-09-13/Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md]]"
---

# Phase FINAL — zamknięcie cyklu op-r3-stationary-states

**Status: CLOSED-EXECUTED (2026-09-14, jedna sesja: LOCK → MD FROZEN
→ Phase 1 → Phase 2 (2 korekty implementacyjne, progi nietknięte)
→ Phase 3 (22 biegi główne + 8 kontroli dt/2) → zamknięcie).**
Kryteria LOCKa stosowane DOSŁOWNIE; zero zmian
kryteriów/progów/detektora/rodzin startów/sponge po pierwszym biegu
produkcyjnym.

---

## 0. Werdykty

| Pytanie | Werdykt | Jedno zdanie |
|---|---|---|
| **Q-E** (istnienie oscylonów — RACHUNEK CENTRALNY) | **Q-E-INCONCLUSIVE** (litera §2 Phase 3) | Zero OSCILLON; 2/10 RADIATED zbieżnie; 6/10 BREAKDOWN-BOUNDARY zbieżnie (kategoria deskryptywna); 2/10 niezbieżne — litera FAIL wymaga „wszystkie RADIATED", więc INCONCLUSIVE |
| **Q-F** (struktura R3) | **NIEURUCHOMIONE** | Warunkowe na Q-E-PASS (LOCK §2 Phase 4) |
| P1a (linearyzacja stacjonarna) | **PASS** | ∇²f=−κ²f, κ²=(M(1)ω²−𝒰″(1))/𝒦(1)=ω²−1 wyprowadzone sympy; κ=1⟺ω²=2 — zgodność z linearyzacją R3 |
| P1b (Lindstedt–Poincaré, PRE-REJESTRACJA) | **ω₂=−139/24<0** | Miękka nieliniowość ⟹ oscylony małej amplitudy oczekiwane; domena predykcji (a→0) POZA zamrożoną rodziną startów (|a|≥0.15) — bez rozstrzygnięcia |
| P1c (gate sympy vs float) | **PASS 27/27** | ≤1.3e−15 (próg 1e−12·max(1,|v|)), ψ∈{0.9,1,1.1}, 9 form |
| P2 (bramka maszynerii) | **PASS 6/6** | próżnia 0.0; dyspersja ω²=k²+1 do 1.5e−4–1.9e−3 (≤1%); dryf energii 1.9e−9 (≤1e−6); odbicie sponge 1.9e−4 (≤1e−3) — po korektach 1–2 |

**Uwagi zalockowane, stosowane dosłownie:** INCONCLUSIVE ≠ pozytyw;
Q-E bez PASS ⟹ konsekwencje R3 pozostają CONDITIONAL-ON-BRANCH bez
wykazanego nośnika w klasie zbadanej; ZAKAZ claimów o masach leptonów.

## 1. Wejścia i formy (rejestr; MD §1, §10)

- Formy CYTAT z zamkniętego poprzednika (op-action-audit
  Phase1_output.txt; Q-D1-PASS; konwencja |g^tt| — user-gate N2):
  M=ψ⁶/(4−3ψ)², M′=12ψ⁵(2−ψ)/(4−3ψ)³, 𝒦=ψ⁴, 𝒰=γ(ψ⁴/4−ψ³/3);
  π=Mψ̇; EOM: Mψ̈+½M′ψ̇²=(1/r²)(r²𝒦ψ′)′−½𝒦′ψ′²−𝒰′; K_geo=γ=c₀=1.
- Integrator: uogólniony Störmer–Verlet (symetryczny, 2. rzędu) na
  (ψ,π); dt=0.005 (dt/2 przy zdarzeniach); iteracje punktu stałego do
  stagnacji maszynowej (korekta 1b). Odwracalność zmierzona: powrót
  po 2×350 j.cz. do 5.1e−12 (diag2 T2).
- Siatka: radialna przesunięta r_i=(i+½)h, h∈{0.05,0.025}, R=200;
  sponge smootherstep γ₀=1.0 na [160,200]; E_core r≤80 (z korektą 2:
  gęstość 𝒰(ψ)−𝒰(1)=(ψ−1)²(3ψ²+2ψ+1)/12 bez kancelacji).
- Starty (LOCK, deterministyczne, π₀=0): 8 gaussów (a∈{−0.3,−0.15,
  +0.15,+0.25}×σ∈{3,6}), 2 quasi-R3 (±0.2·sin(r)/r·e^{−r²/450}),
  próżnia kontrolna; detektor FROZEN (0.5·E_ref przez ≥100T₀ + ≥50
  przejść; potwierdzenie h,h/2,dt/2 ≤10%); t_max=1000; brak seeda.

## 2. Phase 1 — analityka (Phase1_output.txt)

- **P1a:** pełna linearyzacja EOM na ψ=1+εe^{−iωt}f(r) (sympy,
  simplify=0): −ω²f + f − f″ − 2f′/r = 0 ⟹ κ² = ω²−1. Klasy:
  ω<1 zlokalizowane e^{−|κ|r}/r; ω>1 kontinuum sin(κr)/r.
  **Mapowanie na R3: κ=1 ⟺ ω²=2 — formy zgodne.**
- **P1b (PREDYKCJA pre-rejestrowana PRZED Phase 3):** LP do O(a²)
  dla modu jednorodnego: ω₂ = M₁²/16 + M₁U₃/8 − M₂/8 − 5U₃²/48
  + U₄/16, z M′(1)=12, M″(1)=156, 𝒰‴(1)=4, 𝒰⁗(1)=6 ⟹
  **ω₂ = −139/24 ≈ −5.792 < 0** (kontrola poprawności: dla M≡1
  redukuje się do standardu 3β/8−5α²/12). Miękka nieliniowość ⟹
  ω(a)<m — warunek konieczny oscylonów małej amplitudy SPEŁNIONY.
- **P1c: PASS 27/27** (sympy na binarnym double vs float; max
  |Δ|=1.3e−15).
- Deskryptywnie: tożsamość π-formulacji (π̇=π²M′/2M²+RHS ⟺ EOM)
  sympy simplify=0.

## 3. Phase 2 — bramka maszynerii (Phase2_output.txt; PASS 6/6)

| Gate | Wynik | Wartość |
|---|---|---|
| P2a próżnia (sponge ON, obie siatki, 100T₀) | PASS | ‖ψ−1‖∞ = 0.0 dokładnie; zero alarmów |
| P2b dyspersja (3 biny gate h=0.05) | PASS | k=0.597: 1.48e−4; k=1.005: 1.25e−3; k=1.414: 1.94e−3 (biny 1.79/2.20 deskr.: 2.4/2.6e−3; h=0.025 zgodne) — niezależna walidacja Q-D1 (ω²=k²+1) w pełnej nieliniowej maszynerii |
| P2c dryf energii (obie siatki) | PASS | 1.876e−9 / 1.990e−9 (deskr. max\|E−E0\|/E0 = 7.3e−6 = offset hamiltonianu-cienia ∝dt², nie dryf) |
| P2c odbicie sponge (różnicowe vs R=400) | PASS | 1.931e−4 |

## 4. Phase 3 — Q-E (Phase3_output.txt, Phase3_results/)

**Tabela klasyfikacji per start (h=0.05 i h=0.025 + dt/2 przy
zdarzeniach; czasy zdarzeń zbieżne do ~1%):**

| Start | Klasa | Szczegóły |
|---|---|---|
| g a=−0.30 σ=3 | BREAKDOWN-BOUNDARY (zbieżne h,h/2,dt/2) | t≈6.64/6.635/6.635 |
| g a=−0.30 σ=6 | BREAKDOWN-BOUNDARY (zbieżne) | t≈8.07/7.93/8.10 |
| g a=−0.15 σ=3 | **RADIATED** (zbieżnie) | τ=211.0 na OBU siatkach; E_end/E_ref=1.49e−2; 68 przejść |
| g a=−0.15 σ=6 | BREAKDOWN-BOUNDARY (zbieżne) | t≈19.83/19.82/19.82 |
| g a=+0.15 σ=3 | **RADIATED** (zbieżnie) | τ=206.9 na OBU siatkach; E_end/E_ref=1.42e−2; 66 przejść |
| g a=+0.15 σ=6 | BREAKDOWN-BOUNDARY (zbieżne) | t≈18.19/10.34/18.20 (kategoria zbieżna; czas różny między siatkami — raport) |
| g a=+0.25 σ=3 | INCONCLUSIVE | h05: BREAKDOWN-BOUNDARY-LOWER t≈3.85 (+dt/2 3.85); h025: BREAKDOWN (niefinityczność) t≈3.76 — kategoria niezbieżna |
| g a=+0.25 σ=6 | INCONCLUSIVE | h05: BREAKDOWN-BOUNDARY t≈4.86; h025 i dt/2: BREAKDOWN t≈4.67/4.78 |
| qR3 a=−0.20 | BREAKDOWN-BOUNDARY (zbieżne) | t=4.74/4.74/4.74 |
| qR3 a=+0.20 | BREAKDOWN-BOUNDARY (zbieżne) | t≈2.275/2.27/2.2675 |
| vac (kontrola) | OK | zero alarmów detektora (obie siatki), ψ≡1 dokładnie przez t=1000 |

**WERDYKT Q-E: INCONCLUSIVE** — litera: PASS wymaga ≥1 OSCILLON
(jest 0); FAIL wymaga „wszystkie RADIATED zbieżnie" (są 2/10).

**Deskryptywnie (obowiązkowe, LOCK §2): los startów quasi-R3** —
kształt sin(r)/r (przestrzenny profil R3 n≥1) w dynamice 2. rzędu
gałęzi zdrowej **NIE przeżywa jako stan oscylujący**: obie
polaryzacje (a=±0.2) kolabują do granicy dziedziny w t≈2.3–4.7
(zbieżnie w h, h/2, dt/2), zanim wykonają choć jedną pełną oscylację
T₀.

**Deskryptywna charakterystyka klas:**
- RADIATED (jedyne biegi, które przeżyły t_max): wąskie starty
  |a|=0.15, σ=3 — E_core spada monotonicznie (bez plateau
  oscylonowego): 100% → 37% (t=300) → 10% (t=500) → 1.4–1.5%
  (t=1000); rdzeń po t≈300 to dyspersyjny ogon małej amplitudy
  (~1e−3), ψ(0,t) oscyluje wokół 1 z malejącą obwiednią. Zero śladu
  stabilizacji częstości poniżej progu kontinuum w t≤1000.
- BREAKDOWN-BOUNDARY: sferyczna implozja/kompresja wypycha pole poza
  pas dziedziny (górny 4/3−1e−6 lub dolny 1e−6) w czasie t≈2–20;
  kategoria deskryptywna „pole wybiera granicę" (NIE pozytyw, NIE
  falsyfikacja) — zgodna z BREAKDOWN-BOUNDARY-LOWER poprzednika
  (op-action-audit, pin A=1.30).
- Przy a=+0.25 głęboki kolaps na siatce drobniejszej kończy się
  niefinitycznością (przekroczenie CFL przy ψ→0, c=(4−3ψ)/ψ→∞) tuż
  przed/przy wejściu w pas — kategoria niezbieżna ⟹ INCONCLUSIVE
  wg zamrożonej reguły (deskryptywnie: obie siatki zgodne co do
  KOLAPSU i czasu ~3.8–4.9, różnią się etykietą końcową).

**P1b vs Phase 3 (bez reinterpretacji, forbidden move h):**
predykcja ω₂<0 dotyczy istnienia oscylonów MAŁEJ amplitudy (a→0);
zamrożona rodzina startów zaczyna się od |a|=0.15 (a starty σ=6 /
|a|≥0.25 kolabują). Wynik Phase 3 nie potwierdza i nie obala
predykcji w jej naturalnej domenie — konfrontacja wymaga osobnego
cyklu z małymi amplitudami i dłuższym t_max (czasy życia oscylonów
małych amplitud skalują się ~1/(|ω₂|a²) ≫ 1000 przy a≲0.05).

## 5. Mapowanie na drzewo decyzyjne (LOCK §5)

**Q-E-INCONCLUSIVE → „NEEDS metodologiczny (t_max, sponge, siatki)"**
— dosłownie; konkretyzacja w [[NEEDS.md]]: źródłem INCONCLUSIVE nie
jest maszyneria (P2 PASS 6/6, zdarzenia zbieżne do ~1%), lecz (i)
rodzina startów omijająca domenę predykcji P1b (małe amplitudy),
(ii) dominacja kolapsu do granicy dziedziny (6/10) — pytanie
o status pasa granicznego w dynamice 2. rzędu jest user-gated,
(iii) niezbieżność KATEGORII (nie zjawiska) przy a=+0.25.
Hipoteza ratunkowa (profile R3 jako stany stacjonarne gałęzi
zdrowej) pozostaje CONDITIONAL-ON-BRANCH **bez wykazanego nośnika
w klasie zbadanej**; decyzje o dalszych krokach — user (NEEDS).

## 6. Korekty / incydenty / higiena (anti-Lakatos)

- ✓ LOCK przeczytany w całości przed wszystkim; MD FROZEN przed
  jakimkolwiek kodem; kryteria/progi/detektor/starty/sponge
  niezmienione po pierwszym biegu produkcyjnym; Phase 3 wykonane
  W CAŁOŚCI po zamknięciu korekt Phase 2.
- ✓ **Korekta 1** ([[Phase_correction_note_P2_gates.md]], PRZED
  użyciem wyników; pierwotny output `Phase2_output_pre_correction1.txt`):
  (a) harness P2b czytał kwadrant (+ω,+k) FFT 2D (składowa
  przychodząca) zamiast sumy kwadrantów — dla k=0.6 fałszywy pik
  (odbicie nie wraca w oknie); po korekcie 1.48e−4. (b) kryterium
  stopu iteracji → stagnacja maszynowa.
- ✓ **Korekta 2** ([[Phase_correction_note_2_energy_eval.md]], PRZED
  użyciem wyników; output pośredni `Phase2_output_pre_correction2.txt`):
  mechanizm z noty 1(b) REFUTOWANY dalszą diagnostyką (trajektoria
  identyczna) — rzeczywista przyczyna „dryfu" 7e−6: katastrofalna
  kancelacja U(ψ)−U(1) w EWALUATORZE energii (bias bezwzględny
  ~−6e−10, niezależny od dt/h/amplitudy — dowód: T4d, C1/E₀ ∝ 1/a²;
  odwracalność trajektorii 5e−12). Korekta: tożsamość
  𝒰(ψ)−𝒰(1)=(ψ−1)²(3ψ²+2ψ+1)/12 (sympy simplify=0). Po korekcie
  dryf 1.9e−9. Diagnostyki zachowane: `Phase2_diag_output.txt`,
  `Phase2_diag2_output.txt`, `Phase2_diag3_output.txt`.
- ✓ Zakaz podłóg/barier dotrzymany (jedyna obsługa granic: pas
  klasyfikacyjny; BREAKDOWN* klasyfikowane, nie korygowane);
  INCONCLUSIVE nie reinterpretowane; predykcja P1b nienaruszona po
  Phase 3; zakaz claimów o masach dotrzymany.
- ✓ Rdzeń `.tex`/STATE.md/git NIETKNIĘTE; katalogi innych cykli
  tylko odczyt; pełne ścieżki bez `cd`; weryfikacja `ls` po zapisach
  (zero artefaktów zagnieżdżonych ścieżek); bez /dev/null i heredoc.
- ✓ **Integralność**: hashe SHA256 plików zamrożonych rejestrowane
  w `integrity_snapshot.txt` (migawka przed obliczeniami + stan po
  korektach); Phase0_balance.md, Phase_method_decisions.md, Phase1_*
  NIEZMIENIONE przez cały cykl (weryfikacja przy zamknięciu; w tle
  vaulta działał niezależny agent porządkowy — zero ingerencji
  w pliki cyklu).
- Środowisko: CPython 3.14.2, numpy 2.4.3, scipy 1.17.1,
  sympy 1.14.0 (identyczne z poprzednikami).

## 7. Odczyt (deskryptywnie, bez claimów poza klasą zbadaną)

1. **Maszyneria dynamiki 2. rzędu działa i jest zwalidowana** (pierwsza
   klasa dynamiczna poza gradient flow w programie — N3 poprzednika
   zrealizowane): dyspersja pełnej nieliniowej ewolucji odtwarza
   ω²=k²+1 do 1.5e−4, energia zachowana do 1.9e−9/100T₀,
   odwracalność 5e−12, sponge 1.9e−4.
2. **W zbadanej klasie startów gałąź zdrowa nie wytwarza oscylonów:**
   umiarkowane wąskie pulsy (|a|=0.15, σ=3) w całości wypromieniowują
   (τ≈210 ≪ 628), a wszystkie szersze/głębsze starty — w tym oba
   „w kształcie R3" — kolabują do granicy dziedziny w kilka–kilkanaście
   jednostek czasu. Dziedzina (0,4/3) okazuje się dynamicznie
   „dziurawa" dla dużych zaburzeń w 2. rzędzie — to nowa, twarda
   obserwacja programu (w gradient flow poprzedników kolaps był
   wyjątkiem, tu jest regułą dla 6/10 startów).
3. **Napięcie z predykcją P1b pozostaje otwarte:** miękka nieliniowość
   (ω₂=−139/24) wciąż dopuszcza oscylony małej amplitudy — ale ich
   sondowanie wymaga startów a≲0.1 i t_max ≫ 1000, poza zamrożoną
   rodziną. To najbliższy dobrze postawiony test hipotezy ratunkowej.

## 8. Pliki cyklu

`Phase0_balance.md` (LOCK) · `HANDOFF_PROMPT.md` ·
`Phase_method_decisions.md` (FROZEN) · `engine_core.py` (silnik;
korekty 1b, 2 udokumentowane) · `Phase1_stationary.py` →
`Phase1_output.txt` · `Phase2_gate_dynamics.py` → `Phase2_output.txt`
(+ `Phase2_output_pre_correction1.txt`, `Phase2_output_pre_correction2.txt`)
· `Phase_correction_note_P2_gates.md` ·
`Phase_correction_note_2_energy_eval.md` · `Phase2_diag_gates.py` /
`Phase2_diag2_energy.py` / `Phase2_diag3_energy.py` →
`Phase2_diag_output.txt` / `Phase2_diag2_output.txt` /
`Phase2_diag3_output.txt` · `Phase3_evolve.py` → `Phase3_output.txt`
+ `Phase3_results/` (json+npz per bieg, `verdict.json`) ·
`Phase4_families.py` (warunkowy — nieaktywowany) → `Phase4_output.txt`
· `integrity_snapshot.txt` · `NEEDS.md` · `README.md` (log).
