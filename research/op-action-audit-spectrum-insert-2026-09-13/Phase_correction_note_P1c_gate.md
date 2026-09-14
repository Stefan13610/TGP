---
title: "Correction note 1 — P1c: harness gate'u porównywał wartości w dwóch różnych punktach wejściowych (dokładne 13/10 vs binarny double 1.3) + naiwna forma float M traciła 2.1e−12 przez zaokrąglenie 3ψ wzmocnione 1/(4−3ψ)²; korekta: identyczne wejście binarne + math.fma dla (4−3ψ); progi/punkty/formy NIETKNIĘTE"
date: 2026-09-13
type: correction-note
tgp_owner: research/op-action-audit-spectrum-insert-2026-09-13
status: RECORDED-BEFORE-USE
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_method_decisions.md]]"
---

# Correction note 1 — gate P1c (zapisana PRZED użyciem skorygowanego wyniku)

**Zdarzenie:** pierwszy bieg `Phase1_canonical.py` (output pierwotny
zachowany: `Phase1_output_pre_correction1.txt`) dał P1c FAIL w punkcie
ψ=1.3 dla M: |Δ| = 3.47e−12 > 1e−12 (pozostałe 11/12 porównań PASS,
w tym wszystkie 𝒦 i 𝒰; P1a/P1b bez zastrzeżeń).

**Diagnoza (wykonana przed decyzją, wynik diagnostyki):**
1. **Błąd harnessu gate'u (dwa różne wejścia):** sympy liczone
   w DOKŁADNYM ψ=13/10, implementacja float w double(1.3)
   = 13/10 + 4.44e−17. Przy M′(1.3) ≈ 3.12e4 sam wkład reprezentacji
   wejścia = |M(13/10) − M(double(1.3))| = 1.36e−12 > progu 1e−12 —
   porównanie w dwóch różnych punktach NIE mierzy błędu implementacji,
   wbrew celowi gate'u (LOCK §2: „gate implementacji"; analogicznie
   gate poprzednika porównywał implementację silnika).
2. **Uwarunkowanie naiwnej formy float:** `p**6/(4-3*p)**2`
   w ψ=1.3: zaokrąglenie `fl(3*1.3)` (błąd ~2.2e−16) przechodzi
   w całości do (4−3ψ)=0.1 (błąd względny ~2.2e−15), po podniesieniu
   do kwadratu i pomnożeniu przez |M|≈482.7 daje ~2.1e−12 — zmierzone:
   |M_sympy(double 1.3) − M_float(double 1.3)| = 2.10e−12.

**Korekta (czysto implementacyjna; progi 1e−12, punkty {0.5, 1, 7/6,
1.3}, formy M/𝒦/𝒰 — NIETKNIĘTE):**
1. Harness P1c ewaluuje sympy w IDENTYCZNYM wejściu binarnym
   (`sp.Rational(float(pt))` — dokładna konwersja double), izolując
   błąd implementacji zgodnie z literą „gate implementacji".
   Wartości sympy w dokładnych punktach wymiernych raportowane
   RÓWNOLEGLE (deskryptywnie).
2. Referencyjna implementacja float M(ψ) używa `math.fma(-3.0, p, 4.0)`
   dla czynnika (4−3ψ) (pojedyncze zaokrąglenie; błąd |Δ| spada do
   ~1e−13 w najgorszym punkcie). Ta implementacja jest FROZEN dla
   każdego dalszego użycia M float. Uwaga: silnik Phase 3 używa
   wyłącznie 𝒦 i 𝒰 (statyka) — M nie występuje w żadnym rachunku
   numerycznym cyklu; korekta dotyczy tylko referencji gate'u.

**Wpływ na wyniki:** zero — żaden wynik zależny od P1c nie został
użyty przed korektą; P1a/P1b (symboliczne) niezmienione.
