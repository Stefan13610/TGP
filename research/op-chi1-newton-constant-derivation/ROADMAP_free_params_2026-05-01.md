---
title: "Free-parameters derivation roadmap (post-υ.1/φ.1)"
date: 2026-05-01
status: PROPOSED
type: meta-roadmap
tags:
  - TGP
  - free-params
  - roadmap
  - postulates-to-derived
  - mini-cycles-queue
parent: "[[../../PREDICTIONS_REGISTRY.md]]"
---

# Free-parameters derivation roadmap

**Cel:** zdjąć status "naiwnie postulowane" z stałych wciąż anchorowanych
zewnętrznie (G_N, c₀, ℏ, m_e, M_TGP, Λ, f_a, g_ω.1, β_g_ψ.1.v2). Każdy z
parametrów otrzymuje własny **strukturalny mini-cykl** (3 fazy, 18 sub-tests,
≥4-channel convergence gate).

## Priority queue (ordered by infrastructure readiness × structural impact)

| Order | Cycle | Param | Infrastructure | Strategy | Status |
|:-----:|:-----:|---|---|---|---|
| **1** | **χ.1** | G_N + M_Pl (joint M_TGP) | UV.1, ξ.1, Phase 2.A, sek08 | AS NGFP threshold-match: `G_N = g*/(M_TGP²·ξ_grav)` | **PROPOSED** ✓ ([[program.md]]) |
| **2** | **ω.2** | g coupling ω.1 unique selection | ω.1 Phase 3 (4 LOCK candidates) | Cross-channel B² + CMB β + magnetar E·B precision pick | queue |
| **3** | **ψ.2** | β_g_ψ.1.v2 absolute value | ψ.1.v2 (Adams forces NEG sign), AS NGFP | UV match: β_g = O((m/Λ)²Q²/48π²) z heavy-mode 1-loop | queue |
| **4** | **UV.2** | M_TGP standalone (post-χ.1) | χ.1 partial M_TGP lock, UV.1 N_A=500/57 | ngEHT 2030+ + GUT scale + Λ_QCD orthogonal anchor | queue (post-χ.1) |
| **5** | **ζ.2** | m_e absolute scale | ζ.1 lepton-mass Koide, η.2 Wolfenstein, F4 | EW symmetry breaking + Yukawa + Koide ratios | queue (research-track) |
| **6** | **Λ.1** | Λ EFT cutoffs unification | per-cycle Λ_ψ, Λ_τ, Λ_σ, Λ_ω | cross-cycle consistency → single Λ_TGP | queue (post-UV.2) |
| **7** | **σ.2** | c₀ absolute scale vacuum-substrate | σ.1 dispersion, ψ.1.v2 anisotropy | substrate-action `S[X]` dim-embedding c₀=1 → c₀ = f(M_TGP, ℏ) | queue (post-χ.1+UV.2) |
| **8** | **φ.2** | ℏ absolute scale (quantum substrate-action) | φ.1 classical AXIOM | quantum φ.1 path-integral; ℏ from canonical commutator | queue (research-track far) |
| **9** | **ω.3** | f_a axion decay constant | ω.1 (PVLAS bound), χ.1 M_TGP | f_a = M_TGP / g_ω.1; post-ω.2 + χ.1 | queue (post-χ.1+ω.2) |

## Logical dependency graph

```
                    φ.1 AXIOM (substrate-action)
                           │
          ┌────────────────┼────────────────┐
          ▼                ▼                ▼
    UV.1 NGFP         ξ.1 photon       υ.1 closure
   {g*, η_N*=−2}      ring ξ-factor    (X_ref/X_obs)^{1/3}
          │                │                │
          └──────┬─────────┴────────────────┘
                 ▼
            ┌─── χ.1 (G_N + M_Pl, joint M_TGP partial) ───┐
            │                                              │
            ▼                                              ▼
       UV.2 (M_TGP standalone)                    ω.2 (g_ω.1 unique select)
            │                                              │
   ┌────────┼────────┐                            ┌────────┼────────┐
   ▼        ▼        ▼                            ▼        ▼        ▼
 σ.2     Λ.1      ω.3                          ψ.2      φ.2     [...]
 (c₀)  (Λ unify) (f_a)                       (β_g)    (ℏ quantum)
```

**Reading:** χ.1 jest `entry-gate` — większość downstream cycles depends
na M_TGP partial lock z χ.1. ω.2 / ψ.2 są równoległe (orthogonal sektory:
EM-axion, ψ.1.v2 photon-kinetic).

## Joint vs sequential strategy

**Joint locks (recommended):**
- **χ.1**: G_N + M_TGP partial (orthogonal anchors w X2.3 sub-test)
- **UV.2**: M_TGP standalone (post-χ.1 fixes G_N relation)
- **σ.2 + φ.2**: c₀ + ℏ (joint quantum substrate-action — c·ℏ dimensional pair)

**Sequential locks:**
- **ω.2** → **ω.3** (g_ω.1 selected → f_a = M_TGP/g_ω.1)
- **ζ.2** standalone (Yukawa + EW orthogonal do gravity sektor)

## Falsification overlap

Każdy mini-cykl musi mieć ≥3 channel-distinct falsifiers, w tym ≥1
post-confirm-able (existing data) lub ≥1 forward-LIVE (2027–2035 horizon).
Cross-cycle falsification napoły: jeśli ω.2 wybiera g_ω.1 = κ_TGP, to
χ.1 X1.5 alt (iv) `G_N = κ_TGP²/M_TGP²` musi pozostać FALSIFIED.

## Closure target post-roadmap (2027–2030 horizon)

| Status | Pre-roadmap | Post-roadmap target |
|---|---|---|
| Naiwnie postulowane | G_N, c₀, ℏ, M_TGP, m_e, Λ, f_a, M_Pl | (none — wszystkie DERIVED) |
| Częściowo derywowane | g_ω.1 (4 candidates), β_g_ψ.1.v2 (sign only) | full DERIVED post-ω.2/ψ.2 |
| AXIOM | φ.1, F1, F2, F3, N_gen=3 | (unchanged — fundament) |
| DERIVED (sympy-exact) | α_em, CKM, PMNS, Σm_ν, m_ββ, R_TGP, K's, B²'s, A, ρ̄, η̄, … | + G_N, M_Pl, M_TGP, c₀, ℏ, m_e, Λ, f_a, g_ω.1, β_g_ψ.1.v2 |

**Bilans po pełnym roadmap:** programy 1–9 zamykają **9 wolnych parametrów
fizyki fundamentalnej** strukturalnie do φ.1 + UV.1 + ξ.1 + υ.1 axiomatic
core (5 axioms, 0 free parameters).

## Cross-references

- [[program.md]] — χ.1 program (PROPOSED)
- [[../../PREDICTIONS_REGISTRY.md]] — F6 (κ STRUCTURAL → DERIVED post-χ.1)
- [[../../INDEX.md]] — master cycles ledger
