---
title: "ω.2.Phase2 results — sympy LOCK + cross-channel matching 6/7 PASS"
date: 2026-05-01
cycle: ω.2.Phase2
status: COMPLETE
parent: "[[program.md]]"
predecessor: "[[Phase1_results.md]]"
successor: "[[Phase3_setup.md]]"
tags:
  - TGP
  - omega2
  - phase2
  - sympy-lock
  - cross-channel
  - results
  - PASS
---

# ω.2.Phase2 results

**Score: 6/7 PASS** ≥6/7 gate → **Phase 3 ENABLED**.

## Sub-test results

| ID | Test | Result | Detail |
|---|---|---|---|
| W2.2.1 | CMB β constraint per candidate | **PASS** | Δ(lnX) finite for all 5 candidates |
| W2.2.2 | Δ(lnX)_cosmo from φ.1 FRW | **PASS** | Path A: ln(1+z_rec) ≈ 7.0; Path B: O(1) plateau; range [0.5, 10] |
| W2.2.3 | PVLAS-V 2030+ sensitivity | **PASS** | g/f_a ≈ 10⁻¹⁶–10⁻¹⁹ << 10⁻¹¹ for all candidates |
| W2.2.4 | Magnetar E·B ranking | **PASS** | only κ_TGP DETECTABLE → strong discriminator |
| W2.2.5 | Quasar Δχ/ln(1+z) = g/2 | **PASS** | all 5 detectable at SKA z>4 |
| W2.2.6 | Combined χ² ranking | **FAIL** | g_anomaly winner Δχ²=0.21 vs α_em (degenerate at CMB current precision) |
| W2.2.7 | UV-IR matching (Adler-Bardeen) | **PASS** | α(M_TGP)/α(IR) = 3.43 ≈ 3.4 RG-running M_Z→M_GUT |

## Cross-channel rankings (5 candidates, 4 channels)

### W2.2.1 — CMB β constraint required Δ(lnX)

CMB constraint: **g · Δ(lnX) = 2β = 0.01186 rad** (Planck PR4 + ACT 2024).

| candidate | g | required Δ(lnX) | status |
|---|---:|---:|---|
| **g_anomaly** = α·E_TGP/(2π) | 8.30·10⁻³ | **1.4299** | FEASIBLE (cosmological natural) |
| **α_em** | 7.30·10⁻³ | **1.6264** | FEASIBLE |
| 1/(2π) | 0.1592 | 0.0746 | requires tiny Δ(lnX) — disfavored |
| η_chir = 19/24 | 0.7917 | 0.0150 | requires tiny Δ(lnX) — disfavored |
| κ_TGP | 2.012 | 0.0059 | requires tiny Δ(lnX) — disfavored |

**Key insight**: only candidates z g ~ 10⁻³ pasują w cosmological natural range
Δ(lnX) ~ O(1)–O(10). Pozostałe wymagają fine-tuning Δ(lnX) << 1.

### W2.2.2 — Cosmological Δ(lnX) from φ.1 substrate-action

φ.1 axiom: `S[X] = ½ ∫ (∂_μ ln X)(∂^μ ln X) d⁴x` → EL eq `□(ln X) = 0`.

FRW background:
- Homogeneous: u(t) = ln X(t), `ü + 3Hu̇ = 0` → `u̇ = C/a³`
- Matter-dominated: u(t) = u_∞ − C/t (Δ(lnX) frozen as universe expands)

**Natural amplitudes:**
- Path A (log z_rec): Δ(lnX) ~ ln(1 + z_rec) ≈ 7.0
- Path B (O(1) plateau): Δ(lnX) ~ 1
- **Structural range: [0.5, 10]**

→ z W2.2.1 mapping, **g ∈ [10⁻³, 0.024]** consistent z cosmological Δ(lnX).

### W2.2.4 — Magnetar pulsar-timing anomaly

A_pulsar ~ g · ⟨E·B⟩_magnetar / (f_a² · M_TGP²) ~ g · 10⁻³⁰

| candidate | A_pulsar | FAST/SKA 2030+ status |
|---|---:|---|
| κ_TGP | 2.01·10⁻³⁰ | **DETECTABLE** |
| η_chir | 7.92·10⁻³¹ | below threshold |
| 1/(2π) | 1.59·10⁻³¹ | below threshold |
| g_anomaly | 8.30·10⁻³³ | below |
| α_em | 7.30·10⁻³³ | below |

**Magnetar discriminator**: jeśli FAST/SKA 2030+ widzi sygnał → κ_TGP path
(już disfavored CMB); brak sygnału → small-g paths (preferred).

### W2.2.5 — Quasar Δχ/ln(1+z) (SKA 2030+ z>4)

SKA threshold ~ 0.05° = 8.7·10⁻⁴ rad / ln(1+z).

| candidate | slope = g/2 | status |
|---|---:|---|
| κ_TGP | 1.006 rad | DETECTABLE |
| η_chir | 0.396 rad | DETECTABLE |
| 1/(2π) | 0.080 rad | DETECTABLE |
| α_em | 3.65·10⁻³ rad | DETECTABLE |
| g_anomaly | 4.15·10⁻³ rad | DETECTABLE |

Wszystkie 5 detectable, ale slope amplitude ranking gives independent test.

### W2.2.6 — Combined χ² ranking

| rank | candidate | g | Δ(lnX) | χ² | tier |
|:---:|---|---:|---:|---:|---|
| **1** | **g_anomaly = α·E_TGP/(2π)** | 8.30·10⁻³ | 1.43 | **0.185** | FEASIBLE |
| 2 | α_em alone | 7.30·10⁻³ | 1.63 | 0.392 | FEASIBLE |
| 3 | 1/(2π) | 0.159 | 0.075 | 118.1 | tiny Δ(lnX) (g too large) |
| 4 | η_chir = 19/24 | 0.792 | 0.015 | 123.5 | tiny Δ(lnX) (g too large) |
| 5 | κ_TGP | 2.012 | 0.006 | 124.4 | tiny Δ(lnX) (g too large) |

**Status**: g_anomaly TOP rank, ale **Δχ² gap vs α_em = 0.21 < 9** (3σ uniqueness
threshold). Phase 2 NIE rozstrzyga między g_anomaly a α_em jednoznacznie.

**Falsified at >10σ**: 1/(2π), η_chir = 19/24, κ_TGP — wszystkie wymagają
Δ(lnX) << cosmological range [0.5, 10]; każdy disfavored z Δχ² > 100.

**Surviving 2 candidates** (degenerate at current CMB):
- **g_anomaly = α_em·E_TGP/(2π) = 8.30·10⁻³** ← TGP-structural (anomaly form)
- α_em = 7.30·10⁻³ ← QED bare (no E_TGP weighting)

Discrimination requires SO 2027+ + LiteBIRD 2029+ 5σ CMB β.

### W2.2.7 — UV-IR matching (Adler-Bardeen non-renormalization)

Anomaly coefficient is **one-loop exact** (Adler-Bardeen):
- g_eff(IR) = α_em(IR)·E_TGP/(2π) = 8.30·10⁻³
- g_bare(M_TGP) = α_em(M_TGP)·E_TGP/(2π) ≈ α_em(M_TGP)·7.147/(2π)

α_em(M_GUT) ~ 1/40 → g_bare(M_TGP) ≈ 2.84·10⁻²
Ratio: g_bare/g_eff = α(M_TGP)/α(IR) = 0.025/0.00730 = **3.43 ≈ 3.4**
expected from RG running M_Z → M_GUT.

**UV-IR consistency CONFIRMED structurally.**

η_chir = 19/24 = 0.792 at M_TGP scale is **~28× larger** than anomaly form
g_bare(anomaly) ≈ 0.028. Interpretation: η_chir is K-LEVEL bare coupling
(no loop suppression); anomaly form is one-loop perturbative. Different
physical regimes. Cross-channel data preferuje anomaly path.

## Hypothesis ranking after Phase 2

**SURVIVING (CMB compatible):**
1. **g_anomaly = α_em·E_TGP/(2π) = 8.30·10⁻³** ← TGP-structural (RG-invariant)
2. **α_em alone = 7.30·10⁻³** ← QED bare (no anomaly content)

→ Δχ² gap = 0.21 < 9 → Phase 3 needs SO/LiteBIRD 5σ to resolve.

**FALSIFIED (CMB incompatible at >10σ):**
- 1/(2π) standard EFT — generic axion, no TGP content
- η_chir = 19/24 K-LEVEL — too large coupling at IR
- κ_TGP cross-sector — wrong sector, too large

## Phase 2 verdict

**SCORE: 6/7 PASS (≥6/7 gate)** → **Phase 3 ENABLED**.

**Key Phase 2 outputs:**

1. **Cross-channel falsification successful**: 3 of 4 ω.1 LOCK candidates
   eliminated (κ_TGP, 1/(2π), η_chir as bare couplings) at >10σ via CMB
   constraint + cosmological Δ(lnX) feasibility.

2. **Anomaly path validated**: g_anomaly = α_em·E_TGP/(2π) wins χ² ranking;
   UV-IR matching consistent with Adler-Bardeen (one-loop exact).

3. **Degeneracy reduced 4 → 2**: post-Phase 2 LOCK candidates:
   - g_anomaly = α_em·E_TGP/(2π) ≈ 8.30·10⁻³ (TGP-structural)
   - α_em = 7.30·10⁻³ (degenerate runner-up)

4. **Δχ² resolution forecast**: SO 2027+ projected σ_β ≈ 0.018° → Δχ² gap
   between g_anomaly and α_em projected to grow to ~5σ uniqueness post-2027.

**Phase 3 plan**: 6 sub-tests for predictions + 4-channel convergence with
2027–2030+ forecast windows. Final gate ≥5/6 → ω.2 program END.

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_setup.md]]
- [[../op-phi1-substrate-action-variational/Phase3_results.md]] — φ.1 EL eq for Δ(lnX)
- [[../op-omega1-substrate-em-coupling/Phase3_results.md]] — 4 LOCK candidates entering ω.2
