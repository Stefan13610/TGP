# OP-EHT T1 — PN truncation robustness — results

**Data:** 2026-04-25
**Status:** **4/4 PASS** — deviation +14.56% jest **ROBUST GENUINE**
**Skrypt:** `op_eht_t1_pn_truncation.py`
**Raw output:** `op_eht_t1_pn_truncation.txt`

## Cel

Sprawdzić czy +14.6% deviation w b_crit (impact parameter photon ring,
TGP M9.1'' vs GR Schwarzschild) jest artefaktem PN truncation.

## Metoda

1. Rozszerzono PN tail recurrence c_n do n=15 (analitycznie z M9.1'' ODE).
2. Dla każdego n_max ∈ {5, 6, 7, ..., 15} obliczono:
   - r_ph (areal coordinate of photon orbit)
   - b_crit = r_ph / A_t(r_ph), A_t = (4-3ψ)/ψ
3. Sprawdzono konwergencję: |b(n+1)-b(n)| → 0?

## Wyniki kluczowe

| n_max | r_ph (r_g) | b_crit (r_g) | Δb_crit | dev vs GR |
|-------|-----------|-------------|---------|-----------|
| 5 | 3.890 | 5.972 | — | +14.94% |
| 7 | 3.880 | 5.952 | -0.020 | +14.55% |
| 11 | 3.880 | 5.953 | +0.001 | +14.56% |
| 15 | 3.880 | 5.952 | -0.000 | +14.56% |

**Convergence ratio:** |Δb(n+1)/Δb(n)| ~ **0.269** → geometric series convergent.
**Truncation error at n=15:** < 1×10⁻⁴ (4 orders of magnitude below 14.56%).

## Kryteria PASS

- **T1.1** [PASS]: c_n recurrence extended successfully to n=15.
- **T1.2** [PASS]: Δb_crit(n=15→11) < 0.1% (truncation negligible).
- **T1.3** [PASS]: convergence pattern ROBUST (geometric ratio 0.269 < 1).
- **T1.4** [PASS]: b_crit ≈ 5.952 stable across all n=7..15.

## Werdykt

**T1 closes POSITIVE 4/4** — deviation +14.56% jest **genuine prediction**
M9.1'' hyperbolic metric, NIE artefakt PN truncation. Wszelkie próby
dalszego rozszerzenia PN tail (n>15) nie wniosą zmiany >0.001% w b_crit.

## Implikacje

- Path 1 z EHT_quick §6 (truncation artefact) **WYKLUCZONY**.
- Strong-field deviation jest realna; należy pytanie do T2 (calibration MC),
  T3 (q-renorm), T4 (environment), T5 (ngEHT window).
- M9.1'' jako exact prediction: b_crit_TGP/b_crit_GR = 1.1456 (lock).

## Cross-references

- [[OP_EHT_setup.md]]
- [[research/op7/EHT_quick_results.md]] — punkt wyjścia +14.6%.
- [[research/op-newton-momentum/M9_1_pp_p1_higher_pn.py]] — base recurrence.
