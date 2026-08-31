---
title: "Correction note 1 (tego cyklu) — P3b d=4π N=48: eigsh(k=10) nie może zwrócić 8-krotnie zdegenerowanego poziomu obejmującego pozycje 3–10 (zwraca 7 kopii + następny poziom → przesunięcie pozycji, metryka 0.248); korekta: eskalacja k=32, ncv=240 dla tego punktu kontroli (zapisana PRZED użyciem wyniku)"
date: 2026-08-31
type: correction-note
tgp_owner: research/op-3d-canonical-lattice-2026-08-31
status: RECORDED-BEFORE-USE
related:
  - "[[Phase_method_decisions.md]]"
  - "[[Phase3_output.txt]]"
  - "[[Phase3_addendum_p3b_diag.py]]"
  - "[[Phase3_output_addendum_p3b.txt]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/Phase_correction_note_eigsh.md]]"
---

# Korekta 1 tego cyklu (udokumentowany błąd implementacji; kryteria/progi/siatki NIETKNIĘTE)

**Objaw:** `Phase3_output.txt` (zachowany bez zmian): P3b d=4π N=48:
maxerr = 2.482e−1 > 1e−2 → FAIL (pozostałe 9 punktów P3b: PASS ≤9.5e−3).

**Dowód przyczyny** (`Phase3_addendum_p3b_diag.py` +
`Phase3_output_addendum_p3b.txt`; próżnia ⟹ pełne spektrum FD znane
analitycznie z symbolu):

- k=X: FD-exact ma poziom −0.6879 o krotności 8 na pozycjach 3–10;
  `eigsh(k=10, ncv=80, tol=0)` zwraca 7 kopii i poziom następny (−0.4393)
  — porównanie posortowane rozjeżdża się PRZEZ PRZESUNIĘCIE POZYCJI:
  |eigsh−FD-exact|_max = 2.486e−1, dokładnie obserwowany maxerr.
- Operator poprawny: po eskalacji k=32, ncv=240:
  |eigsh−FD-exact|_max = 1.299e−13; czysty błąd dyskretyzacji
  |FD-exact−continuum| = 3.791e−4 ≪ 1e−2.
- Pozostałe 3 punkty k (Γ/M/R): |eigsh−FD-exact| ≤ 1.3e−13 już przy k=10
  (krotności ≤ pozycji 10) — bez artefaktu.

To ta sama rodzina błędu co korekta 1 poprzednika (tol=1e−6, lockout
kopii); tol=0 usunął wariant tolerancyjny, ale przy krotności 8 kończącej
się DOKŁADNIE na pozycji k=10 ARPACK z k=10 nie może z konstrukcji
reprezentować pełnej przestrzeni własnej. d=4π ma najgęstsze degeneracje
(2π/d najmniejsze) — stąd jedyny dotknięty punkt.

**Zakres skutków:**
- P3a (rachunek centralny, tła d=2π): NIEDOTKNIĘTY — λ_min izolowane,
  krotności ≤3 (translacje) w oknie 10; pokrycie translacji
  kontrolowane niezależnie (coverage=1.0, Rayleigh zgodny).
- P3c: NIEDOTKNIĘTY (poziomy w oknie 10 bez krotności sięgających
  pozycji 10; wynik +1.065 na obu siatkach).
- P3b: dotknięty WYŁĄCZNIE punkt d=4π N=48 (bez tła w Phase 2 —
  kontrola maszynerii, nie wynik fizyczny).

**Korekta (parametr solvera, nie kryterium):** dla punktu kontroli P3b
d=4π N=48: `eigsh(k=32, ncv=240, tol=0)`, 10 najniższych gałęzi z 32.
Wynik skorygowany (policzony w addendum PRZED niniejszym użyciem,
metryka LOCKa bez zmian): maxerr(4 pkt k) = **1.849e−3 ≤ 1e−2 → PASS**
(per-k: Γ 7.135e−4; X 3.791e−4 po eskalacji; M 4.014e−4; R 1.849e−3).
**P3b po korekcie: PASS 10/10 punktów.** Pierwotny `Phase3_output.txt`
zachowany bez modyfikacji.

Zapisano PO diagnostyce, PRZED użyciem skorygowanego P3b w jakimkolwiek
dokumencie wynikowym i PRZED uruchomieniem Phase 4.
