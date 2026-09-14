---
title: "Phase0_balance — LOCK: R3 jako stany stacjonarne gałęzi zdrowej — czy dynamika 2. rzędu (M, 𝒦, 𝒰) z |g^tt| ma długożyciowe zlokalizowane stany oscylacyjne (oscylony) i czy ich struktura odtwarza dyskretność R3 (rodziny węzłowe, bariera)?"
date: 2026-09-14
type: phase0-lock
tgp_owner: research/op-r3-stationary-states-2026-09-14
status: PHASE0-LOCKED
computations_performed: ZERO
authorization: "User 2026-09-14: wybór «Sekwencja minimalnego ryzyka» (dopisek core: konwencja |g^tt| dla dynamiki + reklasyfikacja prop:psi-EOM-R3 + ten LOCK) po analizie N2 cyklu op-action-audit-spectrum-insert"
anti_lakatos_lock: ACTIVE
related:
  - "[[README.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-action-audit-spectrum-insert-2026-09-13/Phase_FINAL_close.md]]"
  - "[[../op-action-audit-spectrum-insert-2026-09-13/NEEDS.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Phase 0 — LOCK cyklu `op-r3-stationary-states`

**ZERO obliczeń wykonanych przed zapisaniem tego dokumentu.**

---

## 0. Pytanie i kontekst

Decyzja user-gate N2 (2026-09-14, dopisek rem:W-sign-axiomatic(iv)
+ rem:psi-EOM-R3-branch-status): konwencją kanoniczną DYNAMIKI jest
odczyt |g^tt| (gałąź stabilna). Cena: R3 ODE (fundament N=3 i stosunków
mas leptonów) przestaje być statyką konwencji kanonicznej. Hipoteza
ratunkowa (zapisana w core jako CONDITIONAL-ON-BRANCH): profile
oscylacyjne R3 to przestrzenne profile stanów STACJONARNYCH gałęzi
zdrowej — ψ=1+e^{−iωt}f(r), ∇²f=−κ²f, κ²=ω²−1 (linearyzacja R3:
κ=1 ⟺ ω²=2); masy = częstości quasi-stacjonarnych wzbudzeń
(oscylony/breathery), nie energie profili statycznych.

Kontekst dojrzałości: operator kinetyczny ZALOCKOWANY poprzednikiem
(Q-D1-PASS): M(ψ)=K_geoψ⁶/(c₀(4−3ψ)²), 𝒦=K_geoψ⁴, 𝒰=γ(ψ⁴/4−ψ³/3);
π=Mψ̇; EOM: Mψ̈+½M′ψ̇² = ∇·(𝒦∇ψ)−½𝒦′|∇ψ|²−𝒰′. K_geo=γ=c₀=1 [LOCK].
Ten cykl realizuje też N3 poprzednika (dynamika 2. rzędu) — pierwsza
klasa dynamiczna poza gradient flow w programie.

**Pytania binarne:**
- **Q-E (istnienie — RACHUNEK CENTRALNY):** czy dynamika 2. rzędu
  gałęzi zdrowej ma DŁUGOŻYCIOWE zlokalizowane stany oscylacyjne
  (oscylony): energia zlokalizowana utrzymująca się ≥100 okresów
  podstawowych (T₀=2π/m=2π) po odpromieniowaniu przejściowym,
  zbieżnie w siatce?
- **Q-F (struktura R3, warunkowe przy Q-E-PASS):** czy zbiór stanów
  długożyciowych organizuje się DYSKRETNIE: (a) skończona liczba
  rodzin różniących się liczbą węzłów obwiedni f(r), (b) częstości
  ω poniżej progu kontinuum (ω<m=1 dla rdzenia związanego lub
  dyskretne ω>1 rezonansów — raportować, która klasa realnie
  występuje), (c) istnienie bariery (rodziny wyżej-węzłowe zanikają
  szybko — analog N=3)?

**Poza zakresem (zalockowane):** ilościowe stosunki mas (206.766 itd.)
— wymagają osobnego cyklu z kalibracją; ten cykl testuje ISTNIENIE
i DYSKRETNOŚĆ. Zakaz claimów o masach leptonów z tego cyklu.

Uwaga interpretacyjna zalockowana: Q-E-FAIL (wszystko promieniuje do
próżni) = hipoteza ratunkowa upada w klasie zbadanej — konsekwencje
R3 pozostają CONDITIONAL-ON-BRANCH bez nośnika w gałęzi zdrowej
(decyzja o dalszym statusie łańcucha leptonowego: user). NIE oznacza
to falsyfikacji liczb R3 jako takich.

## 1. Model ZAMKNIĘTY

- **Równanie (dziedziczone, cytat w MD):** Mψ̈+½M′ψ̇² =
  (1/r²)(r²𝒦ψ′)′ − ½𝒦′ψ′² − 𝒰′, radialnie 3D; formy M,𝒦,𝒰 z cyklu
  op-action-audit (Phase1_output.txt) i sek08a (cytaty); dziedzina
  ψ∈(0,4/3), pas ψ>4/3−1e−6 i ψ<1e−6 = BREAKDOWN-BOUNDARY
  (klasyfikacja); ZERO podłóg/barier.
- **Integrator:** 2. rząd w czasie (leapfrog/Störmer–Verlet na parze
  (ψ, π=Mψ̇) lub równoważny; schemat FROZEN w MD przed startem);
  dt=0.005, kontrola dt/2 przy zdarzeniach; **gate zachowania
  energii** (dryf względny ≤1e−6 na 100 T₀ na czystej próżni
  z małym pulsem).
- **Brzeg:** pudło radialne R=200 z warstwą absorpcyjną (sponge)
  r∈[160,200] (tłumienie γ_sp(r) gładkie, FROZEN w MD); obiekt żyje
  w r≲40 — separacja 4×. Diagnostyka energii liczona w r≤80 (E_core).
- **Starty (zalockowane, deterministyczne):** (i) rodzina gaussów
  ψ=1+a·exp(−r²/2σ²), a∈{−0.3,−0.15,+0.15,+0.25}, σ∈{3,6} (8 biegów);
  (ii) profil quasi-R3: ψ=1+a·sin(r)/r·exp(−r²/2σ_w²), a∈{±0.2},
  σ_w=15 (2 biegi — start „w kształcie" R3 n≥1); (iii) kontrola:
  czysta próżnia (zero alarmów detektora). Siatki h∈{0.05,0.025}.
- **Detektor oscylonu (FROZEN):** E_core(t) po t_transient=50:
  jeżeli E_core(t)≥0.5·E_core(t_transient) przez Δt≥100·T₀ ORAZ
  ψ(0,t) oscyluje (≥50 przejść przez 1) ⟹ kandydat; POTWIERDZENIE:
  zbieżność siatkowa (h, h/2: czas życia różni się ≤10%) i dt/2.
  Pomiar ω: FFT ψ(0,t) w oknie stabilnym (pik dominujący + harmoniki
  raportowane).
- **t_max=1000** (≈159 T₀; wystarcza na próg 100 T₀ + transient).
- **Rejestr WEJŚĆ:** rodziny a,σ; R, sponge, h, dt, progi detektora,
  t_transient, t_max; K_geo=γ=c₀=1; brak seeda.

## 2. Fazy i kryteria (zalockowane)

### Phase 1 — analityka stanów stacjonarnych (sympy + shooting liniowy, zero ewolucji)
- P1a: linearyzacja stacjonarna: ∇²f=−κ²f, κ²=(M(1)ω²−𝒰″(1))/𝒦(1)
  = ω²−1 (wyprowadzić, cytat form); klasy: ω<1 ⟹ κ²<0 (f~e^{−|κ|r}/r,
  zlokalizowane), ω>1 ⟹ oscylacyjne sin(κr)/r (kontinuum); mapowanie
  na linearyzację R3 (κ=1 ⟺ ω²=2) — raport zgodności form.
- P1b: znak nieliniowego przesunięcia częstości (Lindstedt–Poincaré
  do O(a²) dla modu jednorodnego rdzenia): warunek istnienia oscylonu
  małej amplitudy (ω(a)<m — miękka nieliniowość) — wyprowadzenie
  symboliczne z 𝒰‴(1), 𝒰⁗(1), M′(1); wynik raportowany jako
  PREDYKCJA przed Phase 3 (pre-rejestracja kierunku).
- P1c (gate): sympy vs float 1e−12 w {0.9, 1, 1.1}.

### Phase 2 — bramka maszynerii dynamicznej
- P2a: próżnia zostaje (‖ψ−1‖∞≤1e−10 przez 100 T₀; obie siatki).
- P2b: **test dyspersji** — mały puls (a=1e−3): zmierzone ω(k) z FFT
  2D (r,t) vs k²+1 (zgodność ≤1% na ≥3 modach) — niezależna walidacja
  Q-D1 w pełnej nieliniowej maszynerii.
- P2c: gate energii (dryf ≤1e−6/100 T₀) + gate sponge (odbicie
  ≤1e−3 amplitudy padającej).
- FAIL ⟹ STOP.

### Phase 3 — Q-E: RACHUNEK CENTRALNY (ewolucje)
- 10 startów × 2 siatki (+dt/2 przy zdarzeniach); klasyfikacja per
  bieg: OSCILLON (detektor, zbieżnie) / RADIATED (E_core→0 przed
  progiem) / BREAKDOWN-BOUNDARY / INCONCLUSIVE.
- **Q-E-PASS:** ≥1 start daje OSCILLON zbieżny. **Q-E-FAIL:**
  wszystkie RADIATED zbieżnie. **Q-E-INCONCLUSIVE:** reszta.
- Deskryptywnie obowiązkowo: los startów quasi-R3 (ii) — czy kształt
  sin(r)/r przeżywa jako stan oscylujący.

### Phase 4 — Q-F (tylko przy Q-E-PASS)
- Dla każdego OSCILLON: ω (FFT), liczba węzłów obwiedni (średnia
  |ψ−1| po oknie), czas życia; grupowanie w rodziny (n, ω).
- **Q-F-PASS:** ≥2 rozróżnialne rodziny węzłowe ORAZ monotonia
  ω(n) ORAZ brak rodzin n≥3 o czasie życia ≥ progu (bariera-analog).
  **Q-F-PARTIAL:** tylko n=0. **Q-F-FAIL:** brak struktury dyskretnej.

## 3. Forbidden moves
Dziedziczone (op-action-audit §3) + (a) formy M,𝒦,𝒰 wyłącznie
z zamkniętego poprzednika/sek08a z cytatem; (b) zakaz podłóg/barier
i zakaz modyfikacji sponge po pierwszym biegu; (c) detektor oscylonu
i progi niezmienialne po pierwszym biegu; (d) rdzeń .tex/STATE.md/git
NIETYKANE (dopiski core zrobiła sesja główna PRZED tym LOCKiem);
(e) katalogi innych cykli tylko odczyt; (f) INCONCLUSIVE ≠ pozytyw;
(g) ZAKAZ claimów o masach leptonów (zakres §0); (h) P1b jest
pre-rejestracją — zakaz reinterpretacji po Phase 3.

## 4. Deliverables
`Phase_method_decisions.md` (FROZEN), `Phase1_stationary.py`+output,
`Phase2_gate_dynamics.py`+output, `Phase3_evolve.py`+output+npz/json
per bieg+log, `Phase4_families.py` (warunkowy)+output,
`Phase_FINAL_close.md`, `NEEDS.md`, `README.md`.

## 5. Drzewo decyzyjne
```text
P2 FAIL → STOP (maszyneria; raport)
Q-E-PASS + Q-F-PASS → hipoteza ratunkowa POTWIERDZONA strukturalnie:
   sukcesy R3 mają kandydata nośnika w gałęzi zdrowej; NEEDS (user):
   cykl ilościowy widma ω vs stosunki mas + dopisek core (zmiana
   statusu CONDITIONAL-ON-BRANCH)
Q-E-PASS + Q-F-PARTIAL → oscylony są, struktury rodzin brak w klasie
   zbadanej; NEEDS: szersza rodzina startów / dłuższe t_max (user)
Q-E-PASS + Q-F-FAIL → nośnik istnieje, ale nie odtwarza dyskretności
   R3 — raport wprost; user decyduje o statusie łańcucha leptonowego
Q-E-FAIL → hipoteza ratunkowa upada w klasie zbadanej; konsekwencje
   R3 pozostają CONDITIONAL-ON-BRANCH bez nośnika; user-gate
   (alternatywy: sprzężenie z materią N2-M911, geneza poziomu 0)
Q-E-INCONCLUSIVE → NEEDS metodologiczny (t_max, sponge, siatki)
```

---

**LOCK ZAMKNIĘTY 2026-09-14. Zmiany poniżej tej linii po starcie
obliczeń = forbidden move.**
