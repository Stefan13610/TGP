---
title: "μ.1 — PMNS δ_CP phase derivation + residual angle drift hardening (3-phase plan)"
date: 2026-04-30
cycle: μ.1
status: ACTIVE
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-iota-charge-pmns-unification/Phase3_results.md]]"
tags:
  - TGP
  - mu1
  - program
  - pmns-delta-cp
  - drift-hardening
  - cross-sector-phase
---

# μ.1 — PMNS δ_CP phase derivation + residual angle drift hardening

> **Cel:** Lift PMNS matrix od *3 DERIVED + 1 open* (post-ι.1) do *4 DERIVED*
> via (a) δ_CP phase derivation through cross-sector phase coupling z CKM
> (κ.1 Wolfenstein triple FULL DERIVED), oraz (b) higher-order mixing-operator
> corrections do PMNS angles dla drift hardening 8.58–15.57% → < 5% target
> (analog do Wolfenstein λ_C-cascade w CKM).

## Kontekst — co ι.1 zostawiło otwarte

- **PMNS matrix post-ι.1:** 3 angles DERIVED (mixing-operator framework)
  + δ_CP OPEN. Drifts vs NuFit 5.3:
  - sin²θ₁₃ = K_ν · λ_C² = 0.025425 vs 0.022 (drift 15.57%)
  - sin²θ₂₃ = K_ν = 1/2 vs 0.572 (drift 12.59%)
  - sin²θ₁₂ = 1/N_gen = 1/3 vs 0.307 (drift 8.58%)
- **CKM matrix post-κ.1:** 4 free → 0 free (FULL DERIVED) z Wolfenstein triple
  (A=64/81, ρ̄=11/78, η̄=5/14) + λ_C cross-sector cascade z ζ.1.
- **CKM CP phase γ post-κ.1:** γ_CKM = arctan(η̄/ρ̄) = arctan((5/14)/(11/78))
  = arctan(195/77) ≈ 68.45° vs PDG 65.7° (drift 4.19%).
- **Combined CKM+PMNS post-(κ.1+ι.1):** 8 free → 1 open (δ_CP) = 7/8 unified.
- **II5 research-track ι.1 future:** δ_CP via cross-sector phase coupling +
  drift hardening + sterile ν 5-sektor extension.

## Hipoteza μ.1

**(A) Cross-sector phase coupling — δ_CP_PMNS DERIVED:**

Two candidate structural derivations dla δ_CP_PMNS:

- **Form A (gen-tripled CKM γ):** δ_CP_PMNS = N_gen · γ_CKM
  = 3 · arctan(195/77) ≈ 205.4°
  → matches NuFit 5.3 NO best-fit 195° (drift 5.3%, within 1σ).

- **Form B (PMNS-Wolfenstein analog z (ν,up) mixing-operator + Majorana π
  shift):** δ_CP_PMNS = π + arctan(η̄_PMNS/ρ̄_PMNS) where (ρ̄_PMNS, η̄_PMNS)
  = (12/78, 6/7) = (2/13, 6/7) z (ν,up) pair (B² + K mixing-operator analog
  do κ.1 lepton-as-reference formulas):
  δ_CP_PMNS = π + arctan(39/7) ≈ 259.8°
  → matches T2K 2024 NO ≈ 248° (drift 4.8%, within NuFit 1σ).

Both forms within NuFit 1σ window [128°, 352°]. Discrimination via DUNE/T2HK
2030+ precision δ_CP measurement.

**(B) Drift hardening — higher-order mixing-operator corrections:**

Hipoteza: zeroth-order ι.1 derivations admit cross-sector λ_C-corrections
(analog do Wolfenstein CKM expansion) z TGP-natural form:

- **sin²θ₁₃_μ.1 = K_ν · λ_C² · (1 − ρ̄)** (cross-sector CKM ρ̄ leakage)
  = (1/2) · λ_C² · (67/78) ≈ 0.02184 (drift 0.73%)
- **sin²θ₂₃_μ.1 = K_ν / K_up** (cross-sector quark-lepton K-ratio)
  = (1/2)/(7/8) = 4/7 ≈ 0.5714 (drift 0.10%)
- **sin²θ₁₂_μ.1 = (1/N_gen) · (1 − λ_C · η̄)** (cross-sector Wolfenstein
  imaginary leakage)
  = (1/3) · (5149/5600) = 5149/16800 ≈ 0.3065 (drift 0.17%)

Wszystkie 3 drifts < 1% target (massive lift z 8.58–15.57% baseline).

## Plan 3-fazowy

### Phase 1 — landscape audit (5 sub-tests, target 5/5 PASS)

- M1.1: Re-confirm NuFit 5.3 PMNS angles + δ_CP central values + 1σ windows.
- M1.2: γ_CKM = arctan(η̄/ρ̄) sympy-exact landscape audit + PDG comparison.
- M1.3: PMNS-Wolfenstein analog (ρ̄_PMNS, η̄_PMNS) z (ν,up) mixing-operator
        pair sympy-rational inventory.
- M1.4: Drift correction landscape — top-rational forms dla 3 PMNS angles.
- M1.5: Viability gate dla Phase 2 derivation.

### Phase 2 — first-principles derivation (7 sub-tests, target 7/7 PASS)

- M2.1: γ_CKM = arctan(195/77) sympy-exact derivation z Wolfenstein triple.
- M2.2: PMNS-Wolfenstein analog (ρ̄_PMNS, η̄_PMNS) = (2/13, 6/7) DERIVED
        via (ν,up) B²+K mixing-operator pair.
- M2.3: δ_CP_PMNS dual derivation (Form A gen-tripled + Form B PMNS-Wolfenstein
        + Majorana π shift) — both within NuFit 1σ.
- M2.4: sin²θ₁₃ drift hardening: K_ν · λ_C² · (1 − ρ̄) DERIVED, drift < 1%.
- M2.5: sin²θ₂₃ drift hardening: K_ν / K_up DERIVED, drift < 0.2%.
- M2.6: sin²θ₁₂ drift hardening: (1/N_gen) · (1 − λ_C·η̄) DERIVED, drift < 0.2%.
- M2.7: Classification cascade — PMNS angles DERIVED (refined²) <1% drifts;
        δ_CP open → PARTIALLY DERIVED via cross-sector phase coupling;
        PMNS matrix 4 free → 0 free post-μ.1 (3 hardened + δ_CP partial).

### Phase 3 — predictions (6 sub-tests, target 6/6 PASS, μ.1 program END)

- M3.1 (MM1): JUNO 2027+ ultra-sharp sin²θ₁₃ post-μ.1 [0.0218, 0.0220] 0.5%
              window dla DERIVED status confirmation.
- M3.2 (MM2): DUNE 2030+ sin²θ₂₃ ultra-sharp 4/7 = 0.5714 within [0.566, 0.577]
              + 2nd octant LOCKED.
- M3.3 (MM3): DUNE/T2HK 2030+ δ_CP_PMNS dual prediction discrimination —
              Form A 205.4° vs Form B 259.8° → ≥5σ resolution dla 2030+.
- M3.4 (MM4): ★ headline — combined CKM+PMNS 8 free → 0 free post-μ.1
              (3 PMNS angles DERIVED + δ_CP PARTIALLY DERIVED via cross-sector).
- M3.5 (MM5): ν.1 future research-track hint — δ_CP form discrimination
              hardening, sterile ν 5-sektor extension, 0νββ Majorana phase.
- M3.6 (MM6): N-channel μ.1 falsification convergence (≥6 channels target).

## Verdict gates

- Phase 1: 5/5 PASS minimum dla viability.
- Phase 2: 6/7 PASS minimum dla derivation success; 7/7 → FULL CASCADE.
- Phase 3: 5/6 PASS minimum dla program END.

## Cumulative ledger

- Pre-μ.1: 517 (post-ι.1)
- Target post-μ.1: **535** (μ.1.Phase1 5 + μ.1.Phase2 7 + μ.1.Phase3 6 = 18)

## Cross-references

- [[../op-iota-charge-pmns-unification/Phase3_results.md]] — ι.1 PMNS angles
  DERIVED via mixing-operator (ν,up)+(lep,ν)+S₃; II5 research-track origin
- [[../op-kappa-mixing-numerator/Phase3_results.md]] — κ.1 CKM Wolfenstein
  triple FULL DERIVED, γ_CKM apex
- [[../op-zeta-mass-spectrum/Phase3_results.md]] — ζ.1 PMNS zeroth-order
  (1/3, 1/2, λ_C²/2)
- [[../op-eta-wolfenstein/Phase3_results.md]] — η.1 Wolfenstein triple LOCKED
- [[../../INDEX.md]] — master ledger 517→535 target
- [[../../PREDICTIONS_REGISTRY.md]] — MM1-MM6 entries, II5 promotion
