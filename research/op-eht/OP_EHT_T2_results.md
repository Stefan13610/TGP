# OP-EHT T2 — Mass calibration Monte Carlo — results

**Data:** 2026-04-25
**Status:** **1/3 PASS** — Sgr A* DISCRIMINATOR
**Skrypt:** `op_eht_t2_mass_calibration_mc.py`
**Raw output:** `op_eht_t2_mass_calibration_mc.txt`

## Cel

Sprawdzić czy effective EHT envelope (M_BH, D, ring modeling) absorbuje
+14.56% TGP deviation gdy wszystkie systematyki są propagowane przez
Monte Carlo.

## Metoda

- N=10⁰⁰⁰⁰ samples per source.
- Priors:
  - **M87***: M = (6.5 ± 0.7)×10⁹ M☉, D = 16.8 ± 0.8 Mpc, θ_obs = 42.0 ± 3.0 µas (~10%)
  - **Sgr A***: M = (4.297 ± 0.013)×10⁶ M☉ (GRAVITY tight), D = 8.277 ± 0.033 kpc, θ_obs = 51.8 ± 2.3 µas (~4.4%)
- Dla każdego sample obliczono predicted shadow diameter (GR i TGP).
- Residual w jednostkach σ_total: (θ_pred - θ_obs)/σ_θ.

## Wyniki kluczowe

| Source | GR predicted | TGP predicted | EHT obs | TGP within 95% CL |
|--------|------------- |---------------|---------|--------------------|
| M87* | 39.7 ± 4.5 µas | 45.5 ± 5.1 µas | 42.0 ± 3.0 µas | 59.1% [PASS] |
| Sgr A* | 54.1 ± 0.3 µas | 62.0 ± 0.4 µas | 51.8 ± 2.3 µas | 2.1% [FAIL] |

- **Sgr A* mean TGP residual: +4.01σ** (vs GR mean residual: +0.64σ)
- **Joint χ² (2 dof): TGP = 22.68, GR = 5.32** (95% CL threshold = 5.991)

## Kryteria PASS

- **T2.1** [PASS]: M87* TGP within 95% CL of EHT obs (59.1% > 50%).
- **T2.2** [FAIL]: Sgr A* TGP within 95% CL of EHT obs (2.1% << 50%).
- **T2.3** [FAIL]: Combined joint χ² < 5.991 (TGP = 22.68 → 4.4σ tension).

## Werdykt

**T2 closes 1/3 PASS** — **Sgr A* jest DISCRIMINATOR**.
- M87* prior jest na tyle szeroki (~10% systematic, ~10% mass uncertainty),
  że TGP wpada w 95% CL z marginalnym margin.
- Sgr A* GRAVITY 2019 prior jest **20× tighter** w masie, oraz EHT 2022
  ring at 4.4% precision daje 4.0σ tension dla TGP.
- Effective envelope NIE absorbuje +14.56% deviation gdy włączymy Sgr A*.

## Implikacje

- Path 2 z EHT_quick §6 (calibration absorpcja) **CZĘŚCIOWY** —
  fits dla M87* ale NIE dla Sgr A*.
- M9.1'' deviation pozostaje genuine; rozwiązanie musi przyjść z T3/T4/T5.

## Cross-references

- [[OP_EHT_T1_results.md]]
- [[research/op7/EHT_quick_results.md]]
- [[OP_EHT_setup.md]]
