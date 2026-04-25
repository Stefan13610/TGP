# OP-EHT T4 — Non-vacuum environment correction — results

**Data:** 2026-04-25
**Status:** **1/3 PASS** — M87* env OK, Sgr A* env INSUFFICIENT
**Skrypt:** `op_eht_t4_environment.py`
**Raw output:** `op_eht_t4_environment.txt`

## Cel

Sprawdzić czy accretion disk + plasma scattering może shifftować apparent
shadow diameter na tyle aby +14.56% TGP deviation był absorbowany przez
astrophysical environment.

## Metoda

Literature-based estimate σ_env / θ:
- emission ring vs photon ring offset (geometric thickness)
- synchrotron self-absorption (SSA)
- plasma scattering
- variability + flaring

Sources: Psaltis 2020, Younsi 2023, EHT 2019/2022 papers IV/V, Bower 2014.

## Wyniki kluczowe

| Source | M_dot (M☉/yr) | σ_env / θ | 95%-CL reach | residual TGP |
|--------|---------------|-----------|--------------|--------------|
| M87* | ~10⁻⁴ (high) | ±5% | ±9.8% | +8.3% |
| Sgr A* | ~10⁻⁸ (low) | ±4% | ±7.8% | +19.7% |

- **M87***: residual +8.3% < 9.8% reach → environment **CAN absorb**.
- **Sgr A***: residual +19.7% >> 7.8% reach → environment **CANNOT absorb**.
  - Gap remaining: +11.9% (~3.0σ)

## Kryteria PASS

- **T4.1** [PASS]: M87* env can absorb residual (9.8% > 8.3%).
- **T4.2** [FAIL]: Sgr A* env CANNOT absorb deviation (7.8% < 19.7%).
- **T4.3** [FAIL]: BOTH sources reconciled by env (Sgr A* gap 3σ).

## Werdykt

**T4 closes 1/3 PASS** — environment uratowuje M87* ale **Sgr A* tension
persists at ~3σ** even after maximum environment correction.

## Implikacje

- **Sgr A* jest cleanest test photon ring TGP M9.1''** — niska accretion
  rate eliminuje environment confounders, pozostawiając structural deviation.
- Dlatego ngEHT (T5) musi celować na Sgr A* jako PRIMARY falsification target.
- Combined T2+T3+T4 wniosek: jeśli q-renorm scenario (e) z T3 NIE BĘDZIE
  derivable z M9.1'' first principles, M9.2 pivot becomes mandatory after
  ngEHT 2030+ confirms GR shadow at 1% precision.

## Cross-references

- [[OP_EHT_T2_results.md]]
- [[OP_EHT_T3_results.md]]
- [[OP_EHT_setup.md]]
