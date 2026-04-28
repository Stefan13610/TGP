---
title: "Phase 2 setup — substrate-physics upgrade of (α₀, n, ψ_th) from postulate to derivation"
date: 2026-04-28
cycle: BH.1.Phase2
status: PRE-EXECUTION
predecessor: "[[Phase1_results.md]]"
parent: "[[program.md]]"
related:
  - "[[../closure_2026-04-26/alpha_psi_threshold/results.md]]"
tags:
  - TGP
  - BH
  - alpha-psi
  - substrate-physics
  - EFT
  - symmetry
  - cross-sector
---

# Phase 2 — Setup: substrate-physics upgrade of (α₀, n, ψ_th)

> **Cel:** Wykorzystując już istniejący T-α audit (closure_2026-04-26,
> 5/5 PASS), **upgrade'ować** parametry `α(ψ) = α₀(ψ-1)²Θ(ψ-1)` z
> "calibrated/motivated" do "derived" tam gdzie to możliwe, plus dodać
> cross-sector consistency check `α₀ ?= κ_TGP²` z TGP-SC.

---

## Stan pre-Phase 2 (z T-α closure 2026-04-26)

T-α audit **już** wyprowadził baseline parametry:

| Parameter | Value | Status pre-Phase 2 | Source |
|-----------|-------|--------------------|--------|
| ψ_th | **1** | postulate (vacuum motivation V'(Φ_eq)=0) | T-α §3.1 |
| n | **2** | structural minimality (C¹ + WEP + non-overkill) | T-α §3.2 |
| α₀ | **≈ 4.04** | calibrated from scenario (e) photon-ring shift | T-α §3.3 |

T-α flagged jako "open dla Phase 1":
- ❌ ψ_th = 1 first-principles derivation
- ❌ n = 2 mathematical (axiomatic) blocking of alternatives
- ❌ α₀ = 4 derived from RG/geometry (not calibrated)
- ❌ Pełna covariant action S[Φ, g, T_μν, J_μ]

BH.1.Phase2 NIE replikuje pełnego Phase 1 PLAN (15-month covariant
derivation), ale **częściowo upgraduje** te parametry o:
- EFT/symmetry rigorous arguments
- Numerical empirical constraints (WEP-MICROSCOPE-2 projected)
- Cross-sector consistency tests (α₀ ?= κ_TGP² SC)

---

## Hipoteza H_Tα-derived

**H₀ (T-α-baseline):** α(ψ) = α₀(ψ-1)²Θ(ψ-1) z parametrami z T-α closure
pozostają w statusie "calibrated/motivated", bez upgrade'u.

**H₁ (T-α-upgraded):** Co najmniej JEDEN z parametrów (α₀, n, ψ_th)
da się **strukturalnie zaostrzyć** do statusu "derived" w Phase 2:
- ψ_th = 1 z **Z₂ symmetry argument** (wokół vacuum point) — nie tylko vacuum motivation
- n = 2 z **WEP-MICROSCOPE-2 lower bound + non-overkill upper bound** — więcej niż T-α minimality
- α₀ z **cross-sector consistency** α₀ = κ_TGP² (SC) — niezależna calibracja

---

## 7 sub-testów Phase 2

### T2.1 — EFT/symmetry classification of α(ψ)

**Cel:** Pokazać że (ψ - ψ_th)² jest **unikalny** functional form pod
założeniami: (a) Lorentz invariance + diff-invariance (lokalna symetria),
(b) Z₂ reflection symmetry around ψ_th (vacuum point), (c) C¹ smoothness,
(d) minimum non-trivial power.

**Test:**
- Sympy expansion of α(ψ) around ψ = ψ_th
- Z₂ reflection ψ → 2·ψ_th - ψ → α(ψ) → α(ψ): wymaga even powers (ψ-ψ_th)^{2k}
- Smoothness: lowest 2k that's C¹ → 2k = 2 → n = 2
- Higher 2k = 4, 6, ... są subdominant w Taylor expansion blisko progu

**PASS criterion:** (ψ-ψ_th)² unique minimal Z₂-invariant form.

### T2.2 — n ≥ 2 lower bound from WEP-MICROSCOPE-2

**Cel:** Pokazać że WEP-MICROSCOPE-2 projected sensitivity η < 10⁻¹⁷
**wymaga** n ≥ 2 ściśle (n = 1 nie jest dopuszczalne dla projected).

**Test:**
- ψ_Earth - 1 ≈ 7·10⁻¹⁰ (Earth surface gravitational potential)
- α(ψ_Earth)/α₀ = (ψ_Earth - 1)^n
- WEP-induced η ~ α(ψ_Earth)/α(ψ_ph) ~ (7e-10/0.168)^n ~ 4·10⁻⁹·^n
- Required: η_TGP < 10⁻¹⁷ → 4·10⁻⁹·^n < 10⁻¹⁷ → n > 1.93
- Strict: n ≥ 2 minimal integer.

**PASS criterion:** n = 1 violates projected WEP-MICROSCOPE-2; n = 2 passes.

### T2.3 — n = 2 unique under non-overkill upper bound

**Cel:** Zaostrzyć: n = 3 daje suppression 10⁻²⁸ (16 orders of magnitude
beyond projected MICROSCOPE-2), nie ma fizycznej motywacji empirycznej.
n = 2 jest **unique** pod combined empirical constraint.

**Test:**
- η_TGP(n=2) = 4·10⁻¹⁸ (just at MICROSCOPE-2 threshold, safe)
- η_TGP(n=3) = 1.6·10⁻²⁶ (zbędne overkill)
- Argument minimal-empirical: n = 2 minimalna wartość zgodna z empirical bound; brak fizycznej motywacji do n > 2.

**PASS criterion:** combined T2.2 + T2.3 → n = 2 unique under empirical bracket.

### T2.4 — α₀ calibration audit z explicit ξ-tracking

**Cel:** Zaudytować T-α calibration α₀ ≈ 4.04 z explicit traceability
geometric vs sketch factors.

**Test:**
- target_shift = (1/2)·(1 - 3/3.88) = 0.1134 (geom units convention)
- ψ_ph - 1 = 0.168
- α₀ = 0.1134 / 0.168² / ξ = 4.018 / ξ
- T-α used ξ = 1 (Phase 0+ sketch O(1) factor)
- Report: α₀ = 4.018 if ξ = 1; α₀ uncertainty driven by ξ uncertainty

**PASS criterion:** α₀ = 4.02 ± δξ traceable (transparent calibration).

### T2.5 — Cross-sector consistency: α₀ ?= κ_TGP²

**Cel:** Test hipotezy że α₀ z BH sektora (T-α photon-ring calibration)
jest **strukturalnie identyczne** z κ_TGP² gdzie κ_TGP = 2.012 to
spin-fluctuation coupling z TGP-SC sektora.

**Test:**
- κ_TGP = 2.012 (TGP-SC v2 paper, eq:lambda_sf, calibrated on V/Nb/Ta/Mo/Pd)
- κ_TGP² = 4.0481
- α₀ z T-α (T2.4) = 4.018
- Difference: |α₀ - κ_TGP²| / κ_TGP² = |4.018 - 4.048| / 4.048 = 0.74%
- Threshold: < 1% → strong structural hint
- Implications: jeśli match, sugeruje że SC κ_TGP i BH α₀ to dwie strony tej samej fundamental coupling — substrate-matter interaction strength sqrt(α₀) reprezentuje fundamental "TGP charge".

**PASS criterion:** dopasowanie α₀ do κ_TGP² w obrębie 1% → STRUCTURAL HINT for cross-sector unification (Phase 1 PLAN może upgrade'ować to do strict identity).

### T2.6 — Multi-source ψ_ph universality (extended)

**Cel:** Replikować T-α.2 z **rozszerzonym** zakresem M_BH (LISA band,
NS-NS post-merger, intermediate mass BHs) — verify α(ψ_ph) pozostaje
universal.

**Test:**
- Extended source list: Sgr A*, M87*, NGC 1277 (M~1.7e10 M_⊙ supermassive),
  intermediate-mass BHs (M~10³–10⁴ M_⊙), GW200129 (M~88 M_⊙ final),
  NS-NS post-merger remnant (M~2.7 M_⊙ near collapse), Sgr A* progenitor seed (M~10² M_⊙ z early-universe)
- For each source, compute ψ_ph (= 1.168 universal w geom units).
- α(ψ_ph)/α₀ = 0.0282 dla wszystkich (universality check).
- Span M_BH covered: 10²–10¹⁰ M_⊙ (8 orders).

**PASS criterion:** max deviation in α(ψ_ph)/α₀ across all sources = 0 (machine precision).

### T2.7 — Phase 1 PLAN compatibility note

**Cel:** Po-test classification: które parametry Phase 2 zostawia jako
"derived" vs "calibrated" do pełnego Phase 1 PLAN.

**Test:**
- Tabela post-Phase-2 status każdego parametru z update'em z T2.1–T2.6
- Lista open issues dla full covariant Phase 1 PLAN
- Check: czy Phase 2 introducing nowe wymagania paper-rigor

**PASS criterion:** transparent classification map dostarczona; brak hidden assumptions.

---

## Materiał wykonawczy

- **Skrypt:** `phase2_alpha_psi_derivation.py` (sympy + numerical 7 sub-tests)
- **Output:** `phase2_alpha_psi_derivation.txt`
- **Memo:** `Phase2_results.md` (closure)

## Środowisko

```bash
cd TGP/TGP_v1/research/op-bh-alpha-threshold
PYTHONIOENCODING=utf-8 python -X utf8 phase2_alpha_psi_derivation.py 2>&1 | tee phase2_alpha_psi_derivation.txt
```

---

## Constants used

```
psi_th       = 1            (T-α postulate, vacuum point)
n            = 2            (T-α minimality argument)
psi_ph       = 1.168        (M9.2-D universal photon ring)
psi_Earth    = 1 + 7e-10    (Earth surface gravitational potential)
target_shift = 0.1134       ((1/2)·(1 - 3/3.88), geom convention)
kappa_TGP    = 2.012        (TGP-SC v2 calibrated, V/Nb/Ta/Mo/Pd)
WEP_MICR2    = 1e-17        (MICROSCOPE-2 projected sensitivity)
```

---

## Cross-references

- [`program.md`](program.md) — overall 3-phase plan
- [`Phase1_results.md`](Phase1_results.md) — H₀ rejected, Path E unique
- [`../closure_2026-04-26/alpha_psi_threshold/results.md`](../closure_2026-04-26/alpha_psi_threshold/results.md) — T-α baseline (5/5 PASS)
- `../../papers_external/tgp_sc_paper/tgp_sc.tex` — κ_TGP = 2.012 source
- [`../op-sc-alpha-origin/Phase2_results.md`](../op-sc-alpha-origin/Phase2_results.md) — wzorzec mini-cyklu (analog SC.1.Phase2)

---

## Decyzja po Phase 2

- **Optymalna ścieżka:** wszystkie 7 sub-testów PASS, T2.5 daje structural hint (α₀ ≈ κ_TGP²), Phase 2 closes z **upgraded T-α status** (3-of-3 parametrów ma teraz first-principles argument or empirical bracket lock-in).
- **Jeśli T2.5 fail (α₀ ≠ κ_TGP² > 1%):** cross-sector unification falsified; α₀ pozostaje calibrated; Phase 1 PLAN remains required for derivation.
- **Jeśli T2.1 fail (Z₂ argument niepoprawny):** ψ_th uniqueness traci structural support; sektor wraca do T-α status.
