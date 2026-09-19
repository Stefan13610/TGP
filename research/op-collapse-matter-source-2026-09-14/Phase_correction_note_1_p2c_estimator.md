---
title: "Phase_correction_note_1 — korekta estymatora dryfu sekularnego w P2c (okno [0,10T₀] zawiera dt²-owy offset stanu początkowego ΔC2 — estymator mierzył ograniczony efekt hamiltonianu-cienia, nie zmianę sekularną). Progi/konfiguracja/definicje LOCKa NIETKNIĘTE."
date: 2026-09-15
type: phase-correction-note
tgp_owner: research/op-collapse-matter-source-2026-09-14
status: APPLIED-BEFORE-USE
related:
  - "[[Phase2_diag_output.txt]]"
  - "[[Phase2_output_pre_correction1.txt]]"
  - "[[../op-r3-stationary-states-2026-09-14/Phase_correction_note_2_energy_eval.md]]"
---

# Correction note 1 (zapisana PRZED użyciem poprawionego wyniku)

## Zdarzenie

P2c (pierwszy przebieg, `Phase2_output_pre_correction1.txt`): dryf
okienny 2.282e−6 (h=0.05) / 2.298e−6 (h=0.025) > próg 1e−6 ⟹ FAIL.
Przed ogłoszeniem STOP wykonano diagnostykę deskryptywną
(`Phase2_diag_energy.py` → `Phase2_diag_output.txt`; progi/okna LOCKa
nietknięte podczas diagnostyki).

## Diagnostyka (rozstrzygająca)

Konfiguracja P2c, h=0.05, dt∈{0.01, 0.005, 0.0025}:

- **T1 (skalowanie dt):** dryf okienny 9.131e−6 → 2.282e−6 →
  5.705e−7; ilorazy **4.00, 4.00** — czyste skalowanie dt².
  Amplituda oscylacji E wokół fitu: 1.212e−3 → 3.031e−4 → 7.578e−5
  (również 4.00, 4.00). Zero komponentu stałego w dt (w odróżnieniu
  od biasu C1 poprzednika, który był dt-niezależny i był błędem
  ewaluatora).
- **T3 (struktura czasowa):** średnie E po kolejnych oknach 10T₀
  (odchyłki od okna 1, ×1e−6, dt=0.005): +0.00, **+2.29, +2.30,
  +2.29, +2.31, +2.29, +2.26, +1.87, +2.27, +2.28** — SKOK między
  oknem 1 (zawierającym stan początkowy t=0, π≡0) a stałym plateau
  okien 2–10; znaki różnic kolejnych okien „++-+---++" (brak trendu).
  Zmiana sekularna (liniowa w t) dawałaby narastanie ~+0.25, +0.5,
  …, +2.3 — nie występuje.

**Wniosek:** zmierzony „dryf" = dt²-owy offset hamiltonianu-cienia
między oknem [0,10T₀] (stan specjalny t=0) a oknami późnymi — składnik
ΔC2 zidentyfikowany już w diagnostyce poprzednika (correction note 2,
T1: „offset wczesnego okna skaluje się czysto ×dt²; hamiltonian-cień,
zgodny z teorią metod symetrycznych"). Jest to efekt OGRANICZONY
(plateau), nie sekularny; przy amplitudzie startu FROZEN a=0.05
(50× większej niż puls 1e−3 poprzednika) skaluje się ~amplituda² i
przekracza próg 1e−6, choć trajektoria jest zachowawcza.

## Błąd implementacji (udokumentowany)

Zamrożona definicja (MD §5 za MD §7 poprzednika): „średnie po oknach
eliminują ograniczoną oscylację energii właściwą metodom symetrycznym;
**dryf = zmiana sekularna — to ona jest przedmiotem gate'u LOCKa**".
Mój estymator (okna [0,10T₀] vs [90T₀,100T₀], [INPUT-MD]) tej definicji
NIE implementuje: różnica okien jest zdominowana przez ograniczony
offset ΔC2 stanu początkowego (dowód: T1+T3 wyżej), a nie przez zmianę
sekularną. To błąd POMIARU (operacjonalizacji [INPUT-MD]), nie
trajektorii — analogicznie do correction note 2 poprzednika.

## Korekta (wyłącznie estymator; próg i konfiguracja BEZ ZMIAN)

Estymator dryfu sekularnego w P2c: **okna [10T₀,20T₀] vs
[90T₀,100T₀]** (oba w reżimie plateau, poza stanem specjalnym t=0);
przeliczenie na 100T₀ (litera LOCKa „≤1e−6/100T₀"): odstęp środków
okien = 80T₀ ⟹

  dryf ≔ (100/80)·|⟨E⟩_[90T₀,100T₀] − ⟨E⟩_[10T₀,20T₀]| / |⟨E⟩_[0,10T₀]| ≤ 1e−6.

Deskryptywnie raportowane też: stary estymator (offset ΔC2) oraz
nachylenie fitu LSQ na [10T₀,100T₀]. Próg 1e−6, konfiguracja
(λ̃=0.05, gauss a=+0.05 σ=3, sponge OFF, t=700, obie siatki),
wszystkie pozostałe gate'y i definicje — BEZ ZMIAN.

## Przewidywanie (falsyfikowalne, PRZED przebiegiem poprawionym)

Z plateau T3 (dt=0.005): |⟨E⟩_w10 − ⟨E⟩_w2| ~ 1e−8 względnie ⟹
po korekcie dryf ~ 1e−8–1e−7 ≤ 1e−6 ⟹ P2c PASS na obu siatkach.
Jeśli po korekcie dryf > 1e−6 — bramka FAIL i STOP (bez dalszych
korekt estymatora).

## Higiena

- Output pierwotny zachowany: `Phase2_output_pre_correction1.txt`.
- Diagnostyka: `Phase2_diag_energy.py` → `Phase2_diag_output.txt`
  (+ serie E(t) w `Phase3_results/diag_p2c_E_dt*.npy`).
- Zero biegów Phase 3 przed tą korektą (zero reinterpretacji).
- Zmiana kodu: wyłącznie sekcja P2c w `Phase2_gate.py` (okna
  estymatora + raport deskryptywny).
