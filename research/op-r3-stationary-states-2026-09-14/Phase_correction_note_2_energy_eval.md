---
title: "Phase_correction_note_2 — korekta 2: kancelacja zmiennoprzecinkowa w ewaluatorze energii (U(ψ)−U(1) → tożsamość (ψ−1)²(3ψ²+2ψ+1)/12); rewizja błędnej diagnozy z correction note 1(b). Progi/definicje LOCKa NIETKNIĘTE."
date: 2026-09-14
type: phase-correction-note
tgp_owner: research/op-r3-stationary-states-2026-09-14
status: APPLIED-BEFORE-USE
related:
  - "[[Phase_correction_note_P2_gates.md]]"
  - "[[Phase2_diag2_output.txt]]"
  - "[[Phase2_diag3_output.txt]]"
---

# Correction note 2 (zapisana PRZED użyciem poprawionych wyników)

## Rewizja noty 1(b) — uczciwość diagnostyczna

Mechanizm podany w correction note 1(b) (dyssypacja od kryterium
stopu iteracji punktu stałego) okazał się BŁĘDNY: po zmianie
kryterium na stagnację maszynową trajektoria i wynik gate'u energii
nie zmieniły się (7.000e−6 — `Phase2_output_pre_correction2.txt`),
czyli stare iteracje były już zbieżne. Zmiana kryterium stopu
(nota 1(b)) pozostaje w mocy jako nieszkodliwe wzmocnienie; korekta
(a) noty 1 (kwadranty FFT) była trafna — P2b przechodzi (1.5e−4 /
1.3e−3 / 1.9e−3 ≤ 1%).

## Diagnostyka właściwej przyczyny (deskryptywna; diag2 + diag3)

- **T2 (odwracalność):** forward 350 j.cz., π→−π, back 350: powrót
  do ψ₀ z dokł. 5.1e−12, energia odzyskana do 9.0e−11 ⟹ integrator
  dokładnie odwracalny, ZERO dyssypacji trajektorii.
- **T1 (skalowanie dt):** offset wczesnego okna skaluje się czysto
  ×dt² (−1.476e−5 / −3.73e−6 / −2.80e−7 dla dt=0.01/0.005/0.00125)
  — hamiltonian-cień, zgodny z teorią metod symetrycznych; ale
  komponent C1 (różnica okien [0,10]→[90,100]) = −7.0e−6 STAŁY
  w dt (0.3%) i h (0.2%).
- **T4d (rząd amplitudy — rozstrzygający):** puls a=1e−4 (energia
  100× mniejsza): C1/E₀ = −6.9e−4, tj. 100× WIĘKSZY względnie ⟹
  **C1 bezwzględnie stały ≈ −6.1e−10, niezależny od amplitudy** ⟹
  nie jest własnością dynamiki pulsu.

## Przyczyna (zidentyfikowana)

Ewaluator energii liczył gęstość potencjalną jako `Ufun(psi) − U_VAC`
= U(ψ)−U(1): **katastrofalna kancelacja** (|U(1)|=1/12≈0.0833,
różnica ~½(ψ−1)² ~ 1e−12 przy rozproszonym polu ψ−1~1e−6): szum
zaokrągleń ~eps·|U(1)|≈1.2e−17 na punkt, ważony 4πh r_i² (do ~1.4e4
przy r≈200), o systematycznym znaku (zaokrąglenia U(ψ) dla ψ≈1
skorelowane), po ~4000 punktach ⟹ bias do ~1e−9 bezwzględnie.
Bias aktywuje się, gdy pole ROZPŁYWA SIĘ na obszar dużych r (stąd
narastanie t≈60–300 i saturacja — zgodnie z dekompozycją przestrzenną
T3), jest niezależny od dt/h/amplitudy i odwracalny (bias POMIARU,
nie trajektorii). **Przepływ był cały czas zachowawczy: siła używa
U′=ψ²(ψ−1) w formie sfaktoryzowanej (bez kancelacji), a F=−∇E
dokładnie (weryfikacja algebraiczna).**

## Korekta

W `engine_core.energy` (i wszędzie, gdzie liczona jest gęstość
próżniowo odjęta): zamiast `Ufun(p) − U_VAC` używa się DOKŁADNEJ
tożsamości wielomianowej:

  U(ψ) − U(1) = (ψ−1)²(3ψ²+2ψ+1)/12

(weryfikacja: (ψ−1)²(3ψ²+2ψ+1) = 3ψ⁴−4ψ³+1; /12 = ψ⁴/4−ψ³/3+1/12
= U(ψ)+1/12 ✓; tożsamość potwierdzona sympy w nagłówku poprawionego
przebiegu). Forma sfaktoryzowana eliminuje kancelację (~12 rzędów
zysku precyzji względnej przy ψ≈1). **Formy modelu, siła, integrator,
progi, okna gate'ów — BEZ ZMIAN; to korekta wyłącznie ewaluatora
diagnostycznego E** (używanego też przez detektor E_core w Phase 3 —
korekta PRZED pierwszym biegiem produkcyjnym).

Przewidywanie (falsyfikowalne, przed przebiegiem): po korekcie dryf
(okna MD §7) ≈ |ΔC2| ~ 3e−8 ≤ 1e−6 ⟹ P2c-energia PASS na obu
siatkach.

## Higiena

- Outputy pierwotne zachowane: `Phase2_output_pre_correction1.txt`
  (pierwszy przebieg), `Phase2_output_pre_correction2.txt` (po
  korekcie 1).
- Diagnostyki: `Phase2_diag_gates.py` → `Phase2_diag_output.txt`;
  `Phase2_diag2_energy.py` → `Phase2_diag2_output.txt`;
  `Phase2_diag3_energy.py` → `Phase2_diag3_output.txt`.
- Zero biegów Phase 3 przed tą korektą (zero reinterpretacji).
