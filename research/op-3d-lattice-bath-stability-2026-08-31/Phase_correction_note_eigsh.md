---
title: "Correction note 1 — eigsh tol=1e−6 gubi krotności zdegenerowanych wartości własnych (ARPACK lockout); korekta: tol=0 + kontrola pokrycia podprzestrzeni translacyjnej (zapisana PRZED Phase 3)"
date: 2026-08-31
type: correction-note
tgp_owner: research/op-3d-lattice-bath-stability-2026-08-31
status: RECORDED-BEFORE-USE
related:
  - "[[Phase_method_decisions.md]]"
  - "[[Phase0_balance.md]]"
  - "[[Phase1_output.txt]]"
---

# Korekta 1 (udokumentowany błąd implementacji; kryteria/progi/siatki NIETKNIĘTE)

**Objaw (dowód):** w P1a (Phase1_output.txt, zachowany bez zmian) deskryptywna
tabela 10 gałęzi próżni pokazuje błędy O(1) mimo poprawnej najniższej gałęzi
(np. Γ: maxerr(10 gałęzi)=9.94e−1 przy err(najniższa)=6e−14). Diagnostyka
bezpośrednia (operator próżni Γ, N=32, L=2π; multiset FD-symbolu znany
analitycznie):

```
tol=1e-6, ncv=80 : [-1.0000, -0.0032 x3, +0.9936 x3, +1.9904 x3]  (BŁĄD)
tol=0,    ncv=80 : [-1.0000, -0.0032 x6, +0.9936 x3]              (POPRAWNIE)
FD-exact         : [-1.0000, -0.0032 x6, +0.9936 x5, ...]
```

Przy tol=1e−6 ARPACK zwraca 3 z 6 kopii zdegenerowanej wartości −0.0032
(lockout restartów przy luźnej tolerancji — pojedynczy wektor startowy daje
jeden wektor na przestrzeń własną; kopie odzyskiwane są z szumu zaokrągleń,
na co luźny tol nie zostawia iteracji). Posortowane porównanie z gałęziami
dokładnymi rozjeżdża się PRZEZ PRZESUNIĘCIE POZYCJI, nie przez błąd operatora.

**Zakres skutków:**
- Bramka P1a: NIEDOTKNIĘTA (metryka FROZEN = gałąź najniższa,
  niezdegenerowana w każdym zalockowanym k; PASS ważny).
- P1c: NIEDOTKNIĘTE (λ_min izolowane, ekstremalne — zbiega poprawnie).
- P3b (10 gałęzi próżni, degeneracje kubiczne) i P3a (3-krotna degeneracja
  modów translacyjnych w Γ): BYŁYBY skażone — stąd korekta PRZED
  jakimkolwiek użyciem wyników Phase 3.

**Korekta (parametr solvera, nie kryterium):**
1. `eigsh(..., tol=0)` (pełna precyzja maszynowa) we WSZYSTKICH wywołaniach
   Phase 3 (P3a/P3b/P3c); pozostałe parametry FROZEN bez zmian
   (which='SA', k≥10, ncv=80, maxiter=200000, v0 deterministyczny
   pseudolosowy).
2. Kontrola pokrycia (nieusuwalna, dodana): w P3a dla każdego tła i k=Γ
   raportowane pokrycie podprzestrzeni translacyjnej przez zwrócone mody:
   c_i = ‖P_span(mody) ∂_i g‖_B/‖∂_i g‖_B, i=x,y,z; jeżeli c_i < 0.99
   dla któregoś i ⟹ eskalacja k=16, ncv=160 (ta sama macierz; punkt
   raportowany z flagą); dodatkowo iloraz Rayleigha λ_R(∂_i g) liczony
   BEZPOŚREDNIO (niezależnie od ARPACK) jako kontrola λ_trans ~ O(h²).
3. Pierwotne outputy zachowane: `Phase1_output.txt` bez modyfikacji
   (deskryptywna tabela 10 gałęzi P1a nosi ślad artefaktu — wyjaśniona
   niniejszą notą, nie przeliczana).

Zapisano PRZED uruchomieniem Phase2/Phase3 (Phase 1 w trakcie biegu:
P1a zakończone, P1c w toku — żaden wynik Phase 3 nie istnieje).
