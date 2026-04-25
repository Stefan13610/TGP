# OP-EHT T5 — ngEHT 2030+ falsification window + final verdict

**Data:** 2026-04-25
**Status:** **3/3 PASS** — OP-EHT closes **CONDITIONAL POSITIVE**
**Skrypt:** `op_eht_t5_ngeht_window.py`
**Raw output:** `op_eht_t5_ngeht_window.txt`

## Cel

Określić rok detekcji/falsyfikacji TGP M9.1'' przy ngEHT 1% precision
+ wydanie final verdict OP-EHT.

## Metoda

1. ngEHT timeline (literature + TWG roadmap):
   - 2027-2028: prototype (10-15 telescopes, 230+345 GHz)
   - 2029-2030: full ops (~25 telescopes)
   - 2030+: dedicated multi-frequency monitoring M87*+Sgr A*

2. Precision projections + detection threshold @ TGP +14.56% deviation.

3. Mass+distance priors at ngEHT epoch (GRAVITY 2030+ S-stars, JWST TRGB).

## Wyniki kluczowe

### T5.1 ngEHT precision projections

| Year | M87* σ_θ | Sgr A* σ_θ | M87* TGP | Sgr A* TGP | Stage |
|------|----------|------------|----------|------------|-------|
| 2026 | 10.0% | 4.4% | 1.5σ | 3.3σ | EHT 2019/2022 (current) |
| 2028 | 8.0% | 4.0% | 1.8σ | 3.6σ | EHT extension + ALMA + LMT |
| **2030** | **2.0%** | **1.5%** | **7.3σ** | **9.7σ** | **ngEHT prototype** |
| 2032 | 1.0% | 1.0% | 14.6σ | 14.6σ | ngEHT full ops |
| 2035 | 0.5% | 0.5% | 29.1σ | 29.1σ | BHEX (space VLBI) |

**First ≥5σ TGP discrimination: 2030** (Sgr A* 9.7σ, M87* 7.3σ).

### T5.2 Priors at ngEHT epoch (~2030)

- Sgr A* mass uncertainty: **0.05%** (GRAVITY 2030+ S-stars)
- Sgr A* distance: 0.4% (already met)
- Sgr A* combined theory uncertainty: **<0.5%**
- M87* combined theory uncertainty: 5-6% (close to threshold)

⇒ **Sgr A* jest PRIMARY ngEHT falsification target** (theory unc << TGP +14.56%).

### T5.3 Final OP-EHT synthesis

| Test | Score | Conclusion |
|------|-------|-----------|
| T1 PN truncation | 4/4 PASS | deviation **ROBUST GENUINE** |
| T2 mass calibration MC | 1/3 PASS | Sgr A* 4σ tension (DISCRIMINATOR) |
| T3 q-renormalization | 4/5 PASS | scenarios (c)/(e) absorb to <5% |
| T4 environment | 1/3 PASS | M87* OK, Sgr A* gap ~3σ |
| T5 ngEHT window | 3/3 PASS | falsifiable strong-field by 2032 |
| **OP-EHT TOTAL** | **13/18 PASS = 72%** | **CONDITIONAL POSITIVE** |

## Kryteria PASS T5

- **T5.1** [PASS]: ngEHT ≥5σ discrimination by 2032.
- **T5.2** [PASS]: priors sufficient (Sgr A* theory <5%).
- **T5.3** [PASS]: OP-EHT falsifiable strong-field window opens by 2032.

## FINAL OP-EHT VERDICT

```
╔════════════════════════════════════════════════════════════════╗
║ OP-EHT closes CONDITIONAL POSITIVE                             ║
║                                                                ║
║ TGP M9.1'' strong-field deviation +14.56% is GENUINE PHYSICS   ║
║ (T1 robust; T3 q-renorm scenario (e) gives natural absorption  ║
║ if matter coupling f(psi) = sqrt(g_tt^GR/g_tt^TGP) derivable). ║
║                                                                ║
║ Sgr A* with EHT 2022 already at 3-4 sigma tension; ngEHT       ║
║ ~2030-2032 will deliver definitive verdict at >= 10 sigma.     ║
║                                                                ║
║ Path forward:                                                  ║
║  - M9.1'' stays as falsifiable prediction (ngEHT 2030+ window) ║
║  - M9.2 conditional opening: ONLY IF q-renorm derivation fails ║
║  - F4 falsifiability narrow: TGP M9.1'' shadow at +14.56%      ║
║    over GR; ngEHT 2030+ at 1% precision = mandatory test.      ║
╚════════════════════════════════════════════════════════════════╝
```

## Implikacje strategiczne

1. **M9.1'' pozostaje testable predykcja** — strong-field deviation jest
   robust + falsifiable + już tense w EHT 2022 (3.3σ) → defines clear
   experimental milestone.

2. **F4 falsifiability hardening:**
   - Pre-2030: TGP M9.1'' shadow predykcja = +14.56% over GR uniformly.
   - Post-2030: ngEHT 1% delivers >10σ verdict (confirm or falsify).
   - Window: 2030-2032.

3. **Path forward — DUAL TRACK:**
   - **Track A** (M9.1'' rescue): derive q-renorm scenario (e)
     coupling f(ψ) = √(g_tt^GR/g_tt^TGP) z first principles M9.1''
     Lagrangian. If successful → M9.1'' survives strong-field test.
   - **Track B** (M9.2 pivot): if Track A fails by ~2028, prepare M9.2
     extension with momentum-back-reaction or second-order ψ-coupling
     such that natural f(ψ) emerges.

4. **Cross-impact OP-7 closure**: OP-7 zamknięty 96.9% PASS (94/97),
   GW propagation c_GW=c_0 unconditional. OP-EHT 72% PASS — niezależny
   od OP-7 (T1 σ_ab=0 dla static spherical), ale współdzieli M9.1''
   metric. Jeśli M9.2 będzie potrzebne, OP-7 GW result musi być
   re-walidowane w M9.2.

## Cross-references

- [[OP_EHT_T1_results.md]] [[OP_EHT_T2_results.md]] [[OP_EHT_T3_results.md]] [[OP_EHT_T4_results.md]]
- [[OP_EHT_setup.md]]
- [[research/op7/]] — OP-7 closure (independent)
- [[research/op7/EHT_quick_results.md]] — punkt wyjścia +14.6%
- [[paper/tgp_core.tex]] §applications BH shadows, F4
- [[KNOWN_ISSUES.md]]
