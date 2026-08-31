---
title: "Phase3_verdict_ruling — ruling odczytu kryterium Q-FAIL dla d bez istniejącego tła (zapisany PO Phase 3, PRZED Phase 4; progi nietknięte)"
date: 2026-08-31
type: gate-ruling
tgp_owner: research/op-bloch-chain-stability-2026-08-31
status: FROZEN-PRE-PHASE4
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase2_output.txt]]"
  - "[[Phase3_output.txt]]"
  - "[[../op-bath-two-sectors-2026-08-23/Phase1_gate_ruling.md]]"
---

# Ruling: odczyt „ω²_min(d) < 0 dla WSZYSTKICH zalockowanych d"

**Problem (nie domknięty w LOCKu):** kryterium Q-FAIL (LOCK §3 Phase 4)
wymaga „ω²_min(d) < 0 dla WSZYSTKICH zalockowanych d (zbieżnie)".
Phase 2 wykazała, że dla d ∈ {π, 2π} tło niestałe NIE ISTNIEJE (kolaps do
próżni; zgodne z pierwszą całką: okres orbit ≥ 2π z równością tylko
w granicy zerowej amplitudy). Dla d ∈ {π, 2π} wielkość ω²_min(d) jest więc
NIEZDEFINIOWANA (nie ma łańcucha, którego stabilność można by mierzyć) —
LOCK nie przewidział istnienia częściowego (przewidział tylko brak istnienia
dla WSZYSTKICH d ⟹ CLOSED-GATE-STOP, co nie zaszło).

**Ruling (zapisany PRZED uruchomieniem Phase 4; wzorzec: Phase1_gate_ruling
poprzednika — oba odczyty raportowane, progi nietknięte):**

- **Odczyt PRIMARY (domain-restricted):** pytanie Q dotyczy stabilności
  łańcucha TAM, GDZIE ŁAŃCUCH ISTNIEJE. Kwantyfikator „dla wszystkich
  zalockowanych d" przebiega po d z istniejącym (niestałym, zbieżnym) tłem:
  {3π, 4π, 6π}. Nieistnienie dla {π, 2π} raportowane osobno jako negatyw
  istnienia (zawęża zakres, nie zmienia znaku wyniku — próżnia, do której
  te punkty kolabują, jest tachionowa z definicji sektora).
- **Odczyt STRICT-LITERAL:** kwantyfikator po wszystkich 5 zalockowanych d;
  dla {π, 2π} warunek nieewaluowalny ⟹ Q-FAIL formalnie nieosiągalny
  ⟹ werdykt spada do Q-INCONCLUSIVE („pozostałe przypadki").

Oba odczyty zostaną podane w Phase_FINAL_close.md. Werdykt nagłówkowy wg
odczytu PRIMARY, z jawną flagą odczytu strict. Kryteria Q-PASS
i Q-INCONCLUSIVE oraz wszystkie progi liczbowe — bez zmian; Q-PASS jest
rozstrzygalny identycznie w obu odczytach („istnieje d…" — kwantyfikator
egzystencjalny nie cierpi na nieistnienie punktów).

**Stan wejściowy Phase 4 (fakty z Phase 3, bez interpretacji):**
ω²_min zbieżne i ujemne dla wszystkich czterech policzonych teł:
3π/0.7: −1.222191; 4π/0.7: −1.229340; 6π/0.7 (2-garbne): −1.222209;
6π/S3single: −1.229829. Argmin k=0 wszędzie (współmierne z superkomórką 4d);
mod minimalny nie jest modem translacyjnym (overlap 0.000; Goldstone
zidentyfikowany osobno, λ_trans = −2.4e−6…−1.3e−5 ~ O(h²)).
