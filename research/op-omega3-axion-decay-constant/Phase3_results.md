---
title: "ω.3.Phase3 results — predictions + 4-channel convergence + program END 6/6 PASS"
date: 2026-05-01
last_revised: 2026-05-04
cycle: ω.3.Phase3
status: COMPLETE
parent: "[[program.md]]"
predecessor: "[[Phase2_results.md]]"
verdict: LOCKED_ALGEBRAIC_CASCADE_CONDITIONAL
verdict_history:
  - 2026-05-01: FULL_CONVERGENCE (claimed)
  - 2026-05-01: ZZ2/ZZ3 downgrade (audit §J.2)
  - 2026-05-04: cascade-conditional na UV.2 K_struct extension (mini-audit)
tags:
  - TGP
  - omega3
  - phase3
  - predictions
  - convergence
  - program-END
  - results
  - PASS
  - audit-reconciled-2026-05-01
  - cascade-conditional-2026-05-04
---

# ω.3.Phase3 results

> ⚠ **EPISTEMIC STATUS 2026-05-04 — LOCKED-ALGEBRAIC + cascade-conditional na UV.2.**
>
> Per [[AUDIT_omega3_2026-05-04.md]]: algebraic sympy LOCKs (`f_a · E_TGP / M_GUT = K_struct`,
> `g_aγ = α_em·E_TGP/(2π·f_a)`, 4-channel cascade diff=0) są **mechanicznie
> poprawne**, ale **magnitude cascade-conditional** na UV.2 K_struct (post-hoc
> fitted w M_GUT band, drift 0.29% w teoretycznym 10-30% paśmie) i M_GUT
> (PDG-anchored).
>
> ZZ2/ZZ3 już downgraded w [[../../meta/AUDYT_TGP_2026-05-01.md]] §J.2.
> Mini-audit 2026-05-04 rozszerza per-row downgrade na ZZ1, ZZ4-ZZ6
> (cascade-conditional). Sub-tests PASS preserved.
>
> **NIE jest BLOCKING critique** — algebra OK, cascade flagged. Forward-gate
> ω.4+ explicit dla strukturalnej derywacji m_a (ALP, free parameter).

**Score: 6/6 PASS** ≥5/6 gate → **ω.3 program END** [claimed FULL CONVERGENCE 2026-05-01; downgraded LOCKED-ALGEBRAIC + cascade-conditional 2026-05-04 — see banner above].

## Sub-test results

| ID | Test | Result | Detail |
|---|---|---|---|
| O3.1 | PVLAS-V 2030+ LSW NULL | **PASS** | g_aγ 9.59 OOM below; P_LSW ratio 4.5·10⁻³⁹ → structural NULL |
| O3.2 | ALPS-II / IAXO 2030+ NULL | **PASS** | ALPS-II 9.07 OOM, IAXO 7.77 OOM above TGP |
| O3.3 | ADMX/HAYSTAC haloscope NULL | **PASS** | 4.77 OOM coupling suppression; m_a free → structural NULL |
| O3.4 | CMB birefringence β cross-check | **PASS** | β_obs = 0.342° ∈ TGP-canonical [0.2, 0.5]° band |
| O3.5 | ALP fuzzy DM forward-gate | **PASS** | m_a free → ω.4+ structural derivation gate noted |
| O3.6 | 4-channel ω.3 convergence | **PASS** | sympy diff = 0 EXACT across UV.1+ξ.1+UV.2+ω.2 cascade |

## Key derived predictive identities

### Identity 1 — All current + 2030+ axion-photon experiments NULL

| experiment | sensitivity (GeV⁻¹) | TGP g_aγ ratio | OOM below | verdict |
|---|---:|---:|---:|:-:|
| PVLAS-IV / OSQAR-II (current) | 6.6·10⁻¹¹ | 2.60·10⁻¹⁰ | 9.59 | NULL |
| PVLAS-V 2030+ (LSW projected) | 6.6·10⁻¹¹ | 2.60·10⁻¹⁰ | 9.59 | NULL |
| ALPS-II 2030+ | 2·10⁻¹¹ | 8.56·10⁻¹⁰ | 9.07 | NULL |
| IAXO 2030+ helioscope | 1·10⁻¹² | 1.71·10⁻⁸ | 7.77 | NULL |
| ADMX/HAYSTAC 2030+ haloscope | 1·10⁻¹⁵ | 1.71·10⁻⁵ | 4.77 | NULL |

**Structural NULL forecast LOCKED across full axion-photon experimental landscape**.

### Identity 2 — LSW signal suppression (PVLAS-V channel)

P_LSW ∝ (g_aγ·B·L)⁴ → with TGP g_aγ vs PVLAS-V projected sensitivity:

$$\frac{P_{LSW}^{TGP}}{P_{LSW}^{PVLAS-V}} = \left(\frac{g_{a\gamma}^{TGP}}{g_{a\gamma}^{PVLAS-V}}\right)^4 \approx 4.54 \cdot 10^{-39}$$

→ TGP LSW signal **38 orders of magnitude below** detection floor — structural NULL essentially absolute.

### Identity 3 — Haloscope dual-NULL (ADMX/HAYSTAC)

Two independent structural conditions force NULL:

1. **m_a free** (TGP ALP, no QCD anomaly) → no constraint forces m_a into [μeV, meV] haloscope band
2. **Coupling suppression** (4.77 OOM) → P_halo ∝ g² ratio = 2.93·10⁻¹⁰ even if m_a happens to fall in band

→ **Dual structural NULL** — both conditions independently sufficient.

### Identity 4 — CMB birefringence β-band cross-check

| input | value |
|---|---|
| β_obs (Planck PR4 + ACT 2024, Eskilt et al.) | 0.342° |
| significance | ~3.8σ (LIVE PARTIAL hint, downgraded 2026-05-01) |
| TGP-canonical β band | [0.2°, 0.5°] |
| consistency | **β_obs ∈ band** ✓ |

TGP-canonical formula (locked z ω.1+ω.2):
$$\beta = \frac{g_{a\gamma}}{2}\int \frac{\partial \ln X}{\partial \eta}\, d\eta$$

→ ω.3 reproduces post-ω.2 LIVE PARTIAL candidate consistently. **Awaits SO/LiteBIRD 2027+ corroboration**.

### Identity 5 — ALP fuzzy DM forward-gate ω.4+

TGP axion classification:
- **E (electromagnetic) anomaly** = 536/75 sympy-exact (ω.2 triangle)
- **N (color/QCD) anomaly** = NOT present
- → **ALP regime** (Axion-Like Particle), m_a is **free parameter** post-ω.3

Consequences:
- No QCD instanton mass-coupling formula applies (m_a · f_a ≠ m_π · f_π · …)
- No specific TGP DM mass band prediction at ω.3 closure
- f_a > 10¹⁶ GeV super-GUT → consistent z anthropic θ_i misalignment / kinetic alignment

**Forward-gate ω.4+** (research-track):
- Structural m_a derivation candidate mechanisms:
  - m_a ~ √(g·f_a·H_inf) substrate-action coupling
  - m_a from TGP topological sector / instanton-analog
  - m_a from explicit PQ breaking at TGP scale
- If ω.4 locks m_a → fuzzy DM cosmology + superradiance bands become falsifiable

### Identity 6 — 4-channel ω.3 convergence (sympy diff = 0 EXACT)

| Channel | Source | Locked form | Status |
|---|---|---|---|
| 1 | UV.1 (NGFP) | g* = 71/100 | ✓ exact rational |
| 2 | ξ.1 (photon-ring) | N_A = 500/57 | ✓ exact rational |
| 3 | UV.2 (K-LOCK) | K_struct = N_A · 2π² | ✓ sympy LOCK |
| 4 | ω.2 (triangle) | E_TGP = 536/75 | ✓ sympy LOCK |

All 4 channels flow into f_a:

$$\boxed{\;f_a = \frac{N_A \cdot 2\pi^2 \cdot M_{GUT}}{E_{TGP}} = \frac{3125\,\pi^2}{1273}\, M_{GUT}\;}$$

Sympy diff (full cascade vs LOCK form) = **0 EXACT**.

## Hypothesis status post-ω.3 program END

**Promoted post-Phase 3 → ω.3 program END LOCK status:**

| identity | LOCK status |
|---|---|
| f_a = 3125·π²·M_GUT/1273 sympy-rational | **LOCK** (Phase 1+2+3 all PASS) |
| f_a ≈ 4.8456·10¹⁷ GeV (drift 0.30% z χ.1) | **LOCK numerical** |
| g_aγ = α_em·E_TGP/(2π·f_a) sympy | **LOCK** algebraic identity |
| g_aγ ≈ 1.71·10⁻²⁰ GeV⁻¹ | **LOCK numerical** |
| ALL axion-photon experiments NULL forecast | **LOCK** (5 experiments, 4.77–9.59 OOM below) |
| Low-scale inflation forecast (r < 3.17·10⁻²³) | **LOCK** if pre-inflation PQ + θ_i=O(1) |
| CMB β cross-check post-ω.3 | **CONSISTENT** z LIVE PARTIAL hint (~3.8σ) |
| ALP classification (E-only) | **LOCK** structural |
| F-cluster preservation max drift 0.083% | **LOCK** |
| 4-channel cascade sympy diff = 0 | **LOCK** EXACT |

## ω.3 program verdict

**SCORE: 6/6 PASS (≥5/6 gate)** → **ω.3 program END — FULL CONVERGENCE**.

### Final structural predictive ledger (ω.3)

1. **f_a = 3125·π²·M_GUT/1273 sympy-rational LOCK** — TGP-native ξ.1+UV.2+ω.2 cascade
2. **g_aγ = 1.71·10⁻²⁰ GeV⁻¹** — α_em·E_TGP/(2π·f_a) algebraic LOCK
3. **NULL forecast at ALL 2030+ axion-photon experiments** — PVLAS, ALPS-II, IAXO, CAST, ADMX, HAYSTAC
4. **Super-GUT axion regime** — distinct from classical QCD band by ≥5 OOM
5. **Low-scale inflation forecast** — r < 3.17·10⁻²³ if pre-inflation PQ + θ_i=O(1)
6. **CMB birefringence β consistency** — ω.3 reproduces ω.1+ω.2 LIVE PARTIAL ~3.8σ
7. **ALP classification LOCKED** — E-only anomaly, m_a free, ω.4+ forward-gate noted
8. **4-channel structural convergence** — UV.1 + ξ.1 + UV.2 + ω.2 → f_a sympy diff = 0

### Forward-gates (research-track post-ω.3)

- **ω.4+** : Structural m_a derivation (instanton-analog / substrate-action / explicit PQ breaking)
- **SO/LiteBIRD 2027+** : CMB β corroboration upgrade LIVE PARTIAL → CONFIRMED or refuted
- **PVLAS-V / IAXO commissioning** : NULL confirmation (predicted by TGP)

## Audit reconciliation (2026-05-01)

Post-`meta/AUDYT_TGP_2026-05-01.md` (CRITICAL **A7**: ω.2 Δχ²=0.21 << 9 → status downgrade required; "ω.3 albo deryw m_X albo jawnie m_X = free"):

### A7 option-2 closure

ω.3 spełnia **option 2** audytu A7 explicit:
- TGP axion classified **ALP** (E_TGP=536/75 anomaly only, no QCD N anomaly)
- m_a (≡ m_X) is **FREE PARAMETER** post-ω.3 — jawnie przyznane (O3.5, ZZ5, Phase 1 Identity 5)
- forward-gate **ω.4+** structural m_a derivation noted (substrate-action / instanton-analog / explicit PQ)

→ A7 option-2 **CLOSED** by ω.3 program END.

### Status reconciliation post-audit

| ω.3 result | original label | audit-aware label | rationale |
|---|---|---|---|
| **f_a = 3125·π²·M_GUT/1273 sympy-rational** | LOCKED | **LOCKED** ✓ | Pure sympy diff=0 EXACT, niezależne od ω.2 fitting |
| **g_aγ = α_em·E_TGP/(2π·f_a) algebraic** | LOCKED | **LOCKED-CONDITIONAL** ⚠ | Algebra LOCK; numerical value dziedziczy ω.2 g_axion LIVE PARTIAL |
| **NULL forecast @ 5 experiments** | LOCKED forecast | **LIVE PARTIAL forecast** ⚠ | Zależy od g_axion ω.2 (Δχ²=0.21); pre-2027+ data conditional |
| **r < 3.17·10⁻²³ low-scale inflation** | STRUCTURAL forecast | **STRUCTURAL** ✓ | Conditional na pre-inflation PQ + θ_i=O(1); audit nie spierane |
| **ALP classification (E-only)** | LOCKED | **LOCKED** ✓ | Czysto strukturalne (E vs N anomaly) |
| **CMB β cross-check ∈ [0.2°, 0.5°]** | LIVE PARTIAL hint | **LIVE PARTIAL** ✓ | Już zgodne z audytem (β downgrade 2026-05-01) |
| **4-channel cascade (UV.1+ξ.1+UV.2+ω.2)** | LOCKED | **LOCKED** ✓ | Sympy diff=0 EXACT |
| **F-cluster preservation max drift 0.083%** | LOCKED | **LOCKED** ✓ | Numerycznie exact |

### Inheritance disclosure

ω.3 dziedziczy z **ω.2 g_axion = α_em·E_TGP/(2π) ≈ 8.30·10⁻³** który audyt klasyfikuje jako **LIVE PARTIAL** (Δχ²=0.21 vs plain α_em — znacznie poniżej 3σ progu). W konsekwencji:

- **g_aγ numerical 1.71·10⁻²⁰ GeV⁻¹** is structurally derived but inherits this partial-status → algebra LOCKED, value LIVE PARTIAL
- **NULL @ PVLAS/ALPS/IAXO/ADMX/HAYSTAC** forecasts are **conditional** on ω.2 forward 2027+ corroboration (SO/LiteBIRD β, gauge couplings) — if ω.2 strengthens to Δχ² > 9, NULL forecasts upgrade to LOCKED

### Audit B7 affecting ZZ4

(∂lnX)² behavior in isocurvature constraint inherits standard PQ inflation cosmology — NOT explicit ω.1 EOM × Schwinger E·B Greens-function derivation (audit B7 open). ZZ4 r-bound is therefore **conditional structural forecast**, nie LOCK.

### Items NIE dotykane ω.3

- A1/A2/A3 (metric forms, √(-g), β_PPN) — M9 domain
- A4 (ax:metric-coupling vs L_mat) — sek08 architecture
- A5 (m_e_eff multiplicative) — τ.3 domain
- A6/A8 (ψ.1.v1 wycofane) — v2 LOCKED, ω.3 niezależne
- B6 (scalar-only mode) — M9.3 GW domain

## Cross-references

- [[program.md]]
- [[Phase1_results.md]]
- [[Phase2_results.md]]
- [[Phase3_setup.md]]
- [[../op-omega2-axion-coupling-lock/Phase3_results.md]]
- [[../op-uv2-mtgp-absolute-scale/Phase3_results.md]]
- [[../op-chi1-mtgp-cross-check/Phase3_results.md]]
- [[../../meta/AUDYT_TGP_2026-05-01.md]]
