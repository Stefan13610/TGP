---
title: "Phase_correction_note_P2_gates — korekta 1: (a) odczyt kwadrantu FFT 2D w P2b (harness), (b) kryterium stopu iteracji punktu stałego w integratorze (dyssypacja numeryczna złapana przez gate energii). Progi/definicje LOCKa NIETKNIĘTE."
date: 2026-09-14
type: phase-correction-note
tgp_owner: research/op-r3-stationary-states-2026-09-14
status: APPLIED-BEFORE-USE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
  - "[[Phase2_output_pre_correction1.txt]]"
  - "[[Phase2_diag_output.txt]]"
---

# Correction note 1 (zapisana PRZED użyciem poprawionych wyników)

**Kontekst:** pierwszy przebieg Phase 2 (output pierwotny ZACHOWANY:
`Phase2_output_pre_correction1.txt`) dał P2a PASS, P2c-sponge PASS
(1.93e−4), ale P2b FAIL (bin k=0.6: 7.36e−2) i P2c-energia FAIL
(7.0e−6 > 1e−6). Diagnostyka deskryptywna: `Phase2_diag_gates.py`
→ `Phase2_diag_output.txt`. Obie usterki to błędy IMPLEMENTACJI
(harness pomiarowy / kryterium stopu iteracji), nie własności
maszynerii fizycznej. **Progi, definicje gate'ów, detektor, sponge,
rodziny startów, dt, siatki — NIETKNIĘTE.** Żaden bieg Phase 3 nie
istniał przed tą korektą (zero reinterpretacji wyników).

## Defekt (a) — P2b: odczyt złego kwadrantu FFT 2D

- **Objaw:** bin k=0.597: ω_meas=1.2504 vs teoria 1.1646 (7.4e−2),
  identycznie na obu siatkach (artefakt deterministyczny analizy,
  nie dyskretyzacji); biny k≥1.0 zgodne ~2.5e−3.
- **Przyczyna:** dla rzeczywistego u(r,t) widmo fft2 jest hermitowskie;
  fala WYCHODZĄCA e^{i(kr−ωt)} niesie moc w kwadrancie (+ω, −k)
  (indeks przestrzenny N−n_k), a harness czytał A[+ω, +n_k] =
  składową PRZYCHODZĄCĄ. Puls testowy jest wychodzący; odbicie od
  ściany dla k=0.6 (v_g=0.51) nie wraca w oknie T=400 (powrót ~t≈585)
  ⟹ bin (+,+) zawiera wyłącznie przeciek widmowy ⟹ fałszywy pik.
  Biny k≥1.0 przechodziły, bo ich szybsze odbicia częściowo wracają.
- **Dowód (diag D1):** kwadrant (+ω,−k) dla k=0.597: pik ω=1.164427
  (|d|/ω=1.49e−4), amplituda 4.30e+1 vs 1.17e−2 w (+,+) (stosunek
  ~3.7e+3); po sumie kwadrantów: k=0.6→1.5e−4, k=1.0→1.25e−3,
  k=1.4→1.94e−3.
- **Korekta:** odczyt s(ω) = A[:, n_k] + A[:, (N−n_k) mod N] (suma
  mocy obu znaków przestrzennych przy ω>0; poprawna dla fal
  wychodzących, przychodzących i stojących). Puls, okno, biny
  gate'owe, próg 1% — bez zmian.

## Defekt (b) — dyssypacja od kryterium stopu iteracji punktu stałego

- **Objaw:** „dryf" 7.0e−6 > 1e−6, niezależny od h (7.000e−6 vs
  6.984e−6).
- **Diagnostyka (D2, okna 10 T₀, dt vs dt/2):** błąd energii =
  C1(t) + C2·dt², gdzie C2·dt² = −3.68e−6 (dt=0.005) jest STAŁY
  między oknami [0,10] i [90,100] (offset hamiltonianu-cienia metody
  symetrycznej 2. rzędu — zgodny ze skalowaniem ×1/4 przy dt/2:
  −0.92e−6; NIE dryfuje: ΔC2 ≈ 2e−8), natomiast C1 narasta od
  ~−5e−8 do −7.03e−6 (saturacja ~t≈300) i jest NIEZALEŻNY od dt i h
  — to nie jest błąd obcięcia integratora, lecz systematyczna
  DYSSYPACJA od zamrożonego w MD §3 kryterium stopu iteracji punktu
  stałego: tol absolutny 1e−14 z podłogą max(1,‖·‖∞) przy skali
  π ~ 1e−4 dopuszcza residuum względne ~1e−10/krok o systematycznym
  znaku (iteracja kontrakcyjna zbiega monotonicznie), wzmacniane
  wagą r² przy kompresji fali na ścianie (t≈190–300), potem
  saturujące po rozproszeniu amplitud. Dekompozycja pasuje do OBU
  dt z dokładnością ~1%: −7.03−3.70=−10.73 (obs. −10.73);
  −7.03−0.92=−7.95 (obs. −7.96).
- **Korekta:** kryterium stopu iteracji punktu stałego = STAGNACJA
  maszynowa (iteruj aż ‖Δ‖∞=0 lub ‖Δ‖∞ przestaje maleć; cap 200,
  brak spadku poniżej 1e−12·max(1,skala) po cap ⟹ NonConvergence
  ⟹ BREAKDOWN jak w MD §3). Formuły schematu (uogólniony
  Störmer–Verlet, trzy kroki) — BEZ ZMIAN. Gate dryfu, jego okna
  i próg 1e−6 [LOCK] — BEZ ZMIAN (po eliminacji C1 oczekiwane
  ~e−8 w pierwotnej definicji okien z MD §7).

## Zakres zmian plików

- `engine_core.py`: tylko pętle stopu iteracji w `Engine.step`
  (kryterium stagnacji zamiast progu absolutnego).
- `Phase2_gate_dynamics.py`: tylko odczyt widma w `dispersion_run`
  (suma kwadrantów).
- Phase 2 zostanie wykonany PONOWNIE w całości (świeży
  `Phase2_output.txt`); pierwotny output zachowany.
- Hashe plików po korekcie: dopisane do `integrity_snapshot.txt`.
