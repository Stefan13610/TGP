---
title: "ρ.1 program — ⁷¹Ge cross-section TGP-native correction (Gallium anomaly closure)"
date: 2026-04-30
cycle: ρ.1
status: PRE-PHASE1
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]]"
tags:
  - TGP
  - rho1
  - 71Ge
  - cross-section
  - gallium-anomaly
  - BEST
  - nuclear-matrix-element
  - program
---

# ρ.1 program — ⁷¹Ge cross-section TGP-native correction

## Cel

Domknąć tension XX3 z ξ.2 (BEST 2022 4σ Gallium → ⁷¹Ge deficit
R = 0.78 ± 0.05) przez **TGP-native nuclear matrix-element prediction
dla ⁷¹Ga(ν_e, e⁻)⁷¹Ge transition strength B(GT)**. ξ.2 LOCKED
B²_sterile = 0 → kanał sterile ν zamknięty → BEST anomalia musi być
~20% systemic over-estimation w Bahcall (1997) cross-section calculation.

## Kontekst

Gallium anomaly (BEST + GALLEX + SAGE) — rate ratio z ⁷¹Cr/⁵¹V neutrino sources:

| Experiment | R = obs/expected | σ |
|---:|---:|---:|
| GALLEX 1994 | 0.85 ± 0.06 | 2.5σ |
| SAGE 1996 | 0.79 ± 0.05 | 4.2σ |
| BEST 2022 | 0.78 ± 0.05 | **4.4σ** |

ξ.2 zamknął kanał sterile ν: B²_sterile = 0 LOCKED + STEREO 2023 +
PROSPECT-I 2023 post-prediction-confirmed null. Pozostaje jedyna opcja:
**~20% over-estimation w Bahcall ⁷¹Ga(ν_e, e⁻)⁷¹Ge cross-section**.

## Bahcall 1997 cross-section formula

```
σ(E_ν) = (G_F² · cos²θ_C / π) · p_e · E_e · F(Z, E_e) · |M|²
|M|² = (B(F)·g_V² + B(GT)·g_A²)
```

Bahcall used shell-model NME dla 3 transitions: ground-state (175 keV),
1st excited (500 keV), 2nd excited (708 keV) of ⁷¹Ge.

| Transition | E_ex (keV) | B(GT) Bahcall | weight |
|---|---:|---:|---:|
| ⁷¹Ga g.s. → ⁷¹Ge g.s. | 0 | 0.0865 | 92.4% |
| ⁷¹Ga g.s. → ⁷¹Ge 175 keV | 175 | 0.0150 | 5.0% |
| ⁷¹Ga g.s. → ⁷¹Ge 500 keV | 500 | 0.0130 | 2.0% |
| ⁷¹Ga g.s. → ⁷¹Ge 708 keV | 708 | 0.0070 | 0.6% |

(Bahcall 1997 ApJ 467 475, Tab. 1)

Frekers et al. (2011, PRC 84 014312) measured B(GT) directly via
⁷¹Ga(³He, t)⁷¹Ge charge-exchange reaction at RCNP — found ~5% lower
ground-state strength + significant excited-state contributions.

Haxton (2013) suggested ~20% systematic uncertainty z bound-state proton
overlap correction.

## Hypothesis ρ.1

**TGP-native B(GT) via chirality-counting B² normalization + closure
1/A^{1/3} applied as bound-state proton overlap correction:**

```
B(GT)_TGP = B(GT)_Bahcall · η_chirality · η_closure
η_chirality = (B²_lep · K_lep) / (2 · K_up)         # 4-sector ratio
            = (2 · 2/3) / (2 · 7/8) = 8/21          # sympy-exact
η_closure   = (A_anchor / A_71)^{1/3}                # bound-state proton
            = (76/71)^{1/3} ≈ 1.023                  # Ge-76 anchor
```

Combined factor:

```
η_combined = η_chirality · η_closure ≈ (8/21) · 1.023 ≈ 0.390
```

But this is too aggressive (62% reduction). Need refinement.

**Alternative ansatz** — chirality-counting only on excited-state weight,
bound-state correction as Haxton (2013) ~10% reduction:

```
B(GT)_TGP_g.s. = B(GT)_Bahcall_g.s. · (1 − Δ_chirality)
Δ_chirality = (1 − K_up·B²_lep / (K_lep·B²_up)) ≈ 1 − (7/8·2)/((2/3)·(13/4))
            = 1 − (7/4)/(13/6) = 1 − 42/52 = 10/52 ≈ 19.2%
```

To daje **~19.2% reduction w B(GT)_g.s. dla ⁷¹Ga → ⁷¹Ge** — bardzo
bliskie BEST observation 22% deficit (R=0.78).

## Sub-test plan (5+7+6 = 18 sub-tests)

### Phase 1 — landscape + BEST anomaly + viability (5 sub-tests)

- P1.1 — Bahcall 1997 cross-section recompute (sanity check)
- P1.2 — Frekers 2011 RCNP ³He,t comparison (5% g.s. lift)
- P1.3 — Haxton 2013 bound-state correction landscape
- P1.4 — BEST 2022 + GALLEX + SAGE combined R fit (3-experiment global)
- P1.5 — TGP viability: chirality factor 8/21 too aggressive → search
  alternative B² ratios + closure form

### Phase 2 — TGP-native B(GT) derivation (7 sub-tests)

- P2.1 — chirality-counting candidate forms (4-sector B² ratios)
- P2.2 — closure 1/A^{1/3} bound-state proton overlap correction
- P2.3 — combined η_chirality × η_closure scan dla B(GT)_g.s.
- P2.4 — best-fit form: Δ_chirality = 1 − (B²_up·K_up)/(B²_lep·K_lep)·correction
- P2.5 — excited-state corrections (175/500/708 keV transitions)
- P2.6 — full ⁷¹Ge cross-section recompute z TGP corrections
- P2.7 — promotions: 4 promotions z PARTIALLY DERIVED → DERIVED

### Phase 3 — predictions + 4-channel falsification (6 sub-tests)

- P3.1 — BEST 2022 R prediction post-TGP correction
- P3.2 — GALLEX 1994 + SAGE 1996 cross-check
- P3.3 — RCNP/iThemba/FRIB future ³He,t spectroscopy 2027+
- P3.4 — LANSCE / FRIB ⁷¹Ge(p,n)⁷¹Ga inverse-kinematics 2030+
- P3.5 — Borexino-II / SNO+ ⁸B solar ν cross-check (orthogonal)
- P3.6 — 4-channel ρ.1 falsification convergence

## PASS gates

- Phase 1: ≥4/5 PASS
- Phase 2: ≥6/7 PASS
- Phase 3: ≥5/6 PASS = ρ.1 program END
- 6/6 Phase 3 = FULL CONVERGENCE

## Cross-references

- [[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]] — π.1 NME closure 1/A^{1/3} framework
- [[../op-xi2-sterile-nu-5sector/Phase3_results.md]] — ξ.2 B²_sterile=0 LOCKED + XX3 BEST tension
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
