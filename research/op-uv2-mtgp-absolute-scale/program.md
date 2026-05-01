---
title: "UV.2 program — M_TGP absolute scale derivation (M_GUT anchor + TGP dimensionless ratio)"
date: 2026-05-01
cycle: UV.2
status: PROPOSED
parent: "[[../op-chi1-newton-constant-derivation/Phase3_results.md]]"
predecessors:
  - "[[../op-uv-as-ngfp/Phase3_results.md]]"
  - "[[../op-xi-photon-ring/Phase3_results.md]]"
  - "[[../op-chi1-newton-constant-derivation/Phase3_results.md]]"
tags:
  - TGP
  - UV2
  - M_TGP
  - absolute-scale
  - GUT-anchor
  - substrate-scale
  - PROPOSED
---

# UV.2 — M_TGP absolute scale derivation

**Goal:** Złamać tautologiczną cyrkularność z χ.1 (M_TGP anchored z M_Pl PDG)
poprzez **independent M_TGP derivation** używając **M_GUT (gauge-unification scale)
jako single dimensionful anchor** + TGP dimensionless invariants {g*, N_A, c_χ, ...}.
Cykl promotuje M_TGP **DERIVED PARTIAL → DERIVED FULL** i zamyka G3-scale gap.

## Hypothesis

χ.1 ustaliło:
$$M_{TGP} = M_{Pl}\sqrt{g^*/N_A} \approx 0.2845\, M_{Pl} \approx 3.473 \cdot 10^{18}\, \text{GeV}$$

Joint-lock z M_Pl PDG → M_TGP DERIVED PARTIAL (anchored z M_Pl).
Numerycznie: **M_TGP / M_GUT ≈ 173.7** (M_GUT ~ 2·10¹⁶ GeV).

**UV.2 hypothesis:** istnieje TGP-strukturalny dim-less factor `K_struct`
(funkcja {g*, N_A, π, ...}) taki że:
$$\boxed{\;M_{TGP} = K_{\text{struct}} \cdot M_{GUT}\;}$$

z `K_struct` derywowalnym z TGP dim-less inwariantów.

**Kandydat winner (Phase 1 scan):** K_struct = **N_A · 2π²** ≈ 173.18
(geometric: 2π² = vol(S³) unit; N_A z ξ.1 photon-ring inheritance).

**Equivalent dual-anchor relation (consistency):**
$$\frac{M_{Pl}}{M_{GUT}} = \frac{K_{\text{struct}}}{\sqrt{g^*/N_A}} = \frac{2\pi^2 \cdot N_A^{3/2}}{\sqrt{g^*}}$$

→ M_Pl/M_GUT predicted z {g*, N_A} alone — falsyfikowalne przez gauge unification scale + PDG M_Pl.

## Reference frame

| Source | Status quo | UV.2 lift |
|---|---|---|
| χ.1 | M_TGP DERIVED PARTIAL z M_Pl PDG | M_TGP DERIVED FULL z M_GUT + dim-less |
| UV.1 AS NGFP | g* = 71/100 marginal IR | TGP dim-less anchor input |
| ξ.1 | N_A = 500/57 photon-ring | TGP dim-less anchor input |
| Gauge unification | M_GUT ≈ 2·10¹⁶ GeV (sin²θ_W RG) | single dimensionful M-anchor |
| F1 Single-Φ | substrate ↔ all sectors | gauge-grav unification structural |

## Phase plan (5 + 7 + 6 = 18 sub-tests)

### Phase 1 — Structural setup + alt-K_struct falsification (5)

- **U1.1** M_TGP–M_GUT separation hypothesis: substrate scale > gauge unification
  scale by structural factor; gauge-grav unification ~ M_TGP joint-lock z UV.1
  AS NGFP marginal IR.
- **U1.2** K_struct candidates (4 dim-less form scan):
  - (a) K_a = N_A · 2π² (S³-volume × photon-ring invariant)
  - (b) K_b = N_A² · √(2π) (square-N_A correction)
  - (c) K_c = (4π) · N_A · κ_TGP (XS1 cross-sector)
  - (d) K_d = α₀ · 4π² · √N_A (F4 algebraic anchor)
- **U1.3** M_Pl/M_GUT cross-check: prediction vs PDG M_Pl + GUT scale measurement.
- **U1.4** Joint-system consistency: M_TGP–M_Pl–M_GUT triple-lock pod χ.1+UV.2.
- **U1.5** Alt-anchor ansatz falsification (5 alts):
  - (i) M_TGP = M_Pl postulate (circular; χ.1 tautology)
  - (ii) M_TGP = M_GUT direct (no separation; M_TGP/M_GUT = 1 fails)
  - (iii) M_TGP = √(M_Pl·M_GUT) geometric mean (XS dim-mix; structurally orthogonal)
  - (iv) M_TGP = M_GUT/g* (g*-only correction; misses N_A photon-ring)
  - (v) **M_TGP = N_A·2π²·M_GUT** (TEST candidate UV.2 winner)

**Score gate:** ≥4/5 PASS → Phase 2 forward

### Phase 2 — Sympy LOCK + numerical match (7)

- **U2.1** TGP dim-less invariant ledger: enumeruj {g*, λ*, η_N*, N_A, c_χ, α₀, g̃, κ}
  i ich algebraiczne kombinacje na poziomie dim-less.
- **U2.2** K_struct sympy form: derywować K = N_A · 2π² jako combination of
  ξ.1 photon-ring (N_A) × geometric (2π²) — uniqueness w 4-cand scan.
- **U2.3** M_Pl/M_GUT prediction: sympy `M_Pl/M_GUT = 2π²·N_A^(3/2)/√g*`
  vs observed PDG M_Pl / M_GUT(2-loop SM RG). **Drift target:** < 1%.
- **U2.4** M_TGP numerical reproduction: M_TGP_UV.2 = K·M_GUT vs χ.1 M_TGP
  (joint-lock z M_Pl PDG). **Drift target:** < 0.5%.
- **U2.5** G_N → M_Pl chain post-UV.2: G_N predicted z M_GUT + dim-less only;
  drift vs CODATA. **Drift target:** < 0.1%.
- **U2.6** UV-IR cascade self-consistency: AS NGFP UV cutoff Λ_AS → M_TGP IR
  scale matching consistent z marginal η_N* = -2.
- **U2.7** F-cluster post-UV.2: F4/F5/F6 untouched (algebraic + EFT-scaling
  + κ-DERIVED) ↔ XS1 √α₀ = κ_TGP preserved.

**Score gate:** ≥6/7 PASS → Phase 3 forward

### Phase 3 — Predictions + 4-channel convergence (6)

- **U3.1** **M_GUT precision forecast**: gauge unification 2-loop SM ≈ 2·10¹⁶ GeV
  z theoretical uncertainty ~10-30% (threshold corrections, susy assumption).
  UV.2 prediction tightens M_GUT z TGP-side: M_GUT = M_TGP/(N_A·2π²) ≈ 2.005·10¹⁶ GeV.
  **Drift gate:** < 5% vs SM RG central value.
- **U3.2** **M_Pl reproduction post-UV.2**: M_Pl = M_GUT · 2π²·N_A^(3/2)/√g*
  z M_GUT independent. **Drift target:** < 1% vs PDG.
- **U3.3** **G_N(SI) reproduction post-UV.2** (no M_Pl input):
  G_N = g*/(M_TGP²·N_A) = g*/(K²·M_GUT²·N_A). **Drift gate:** < 5·10⁻⁵ CODATA.
- **U3.4** **f_a axion decay constant** post-UV.2 prediction:
  f_a ~ M_TGP / (something) — opens ω.3 mini-cycle z fixed M_TGP scale.
- **U3.5** **gauge unification structural lock**: SM RG running α_1 = α_2 = α_3
  at M_GUT (within threshold uncertainty); UV.2 promotes M_GUT od observational
  → STRUCTURAL anchor z TGP-side derivation.
- **U3.6** **4-channel UV.2 convergence**:

| # | Channel | Form | Method | Target |
|---|---|---|---|---|
| 1 | UV-anchor (g*) | g* = 71/100 (UV.1 NGFP) | UV.1 NGFP | structural |
| 2 | Photon-ring (N_A) | N_A = 500/57 (ξ.1) | ξ.1 photon-ring | structural |
| 3 | M_GUT independent | gauge unification 2-loop SM | sin²θ_W RG | < 5% |
| 4 | M_Pl reproduction | M_Pl PDG | dim-less ratio | < 1% |

**Score gate:** ≥5/6 PASS → UV.2 program END (FULL CONVERGENCE)

## Promotions post-UV.2 (jeśli FULL CONVERGENCE)

- **M_TGP DERIVED FULL**: K_struct·M_GUT z dim-less ratio z TGP invariants
- **M_Pl/M_GUT structural prediction**: 2π²·N_A^(3/2)/√g* dim-less LOCKED
- **M_GUT** od observational → STRUCTURAL anchor (TGP-grav cross-link)
- **G3-scale gap CLOSED**: M_TGP no longer dependent na M_Pl PDG; circularity broken
- **f_a band fixed** (post-UV.2 → ω.3 enabled)
- **All Planck-scale physics** TGP-grounded: M_Pl, M_TGP, G_N, κ wszystkie z {g*, N_A, M_GUT}

## Falsification path

UV.2 jest falsyfikowane przez:
1. K_struct N_A·2π² drift > 1% vs χ.1 M_TGP (no other K candidate w 4-scan ratuje)
2. M_Pl/M_GUT predicted ratio drift > 1% vs PDG·M_GUT (struct mismatch)
3. M_GUT 2-loop SM running drift > 30% (threshold corrections push outside)
4. CODATA G_N drift post-UV.2 > 5·10⁻⁵ (joint-system inconsistent)
5. ngEHT/LIGO-grav 2030+ M_TGP indirect bound drift > 1% (cross-channel falsification)

## Open frontiers post-UV.2 (jeśli zamknięty)

- **f_a** axion decay constant (ω.3 follow-up; M_TGP fixes f_a band)
- **c₀** vacuum-substrate light speed (σ.2 + ψ.1.v2 fusion; orthogonal scale)
- **ℏ** quantum substrate (φ.2; quantum substrate-action)
- **Λ_TGP** EFT cutoffs unification (Λ.1; single Λ across cycles)
- **m_e** electron mass (ζ.2; electron Yukawa from substrate)

## Cross-references

- [[../op-chi1-newton-constant-derivation/Phase3_results.md]] — χ.1 M_TGP DERIVED PARTIAL
- [[../op-uv-as-ngfp/Phase3_results.md]] — UV.1 g* = 0.71, η_N* = −2
- [[../op-xi-photon-ring/Phase3_results.md]] — ξ.1 N_A = 500/57
- [[../op-phase2-quantum-gravity/Phase2_A_results.md]] — F6 STRUCTURAL→DERIVED
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]] — XX2 entry (M_TGP DERIVED PARTIAL → FULL post-UV.2)
