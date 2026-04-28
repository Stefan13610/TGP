---
title: "Phase 3 setup — multi-source falsification map of α(ψ) threshold function"
date: 2026-04-28
cycle: BH.1.Phase3
status: PRE-EXECUTION
predecessor: "[[Phase2_results.md]]"
parent: "[[program.md]]"
tags:
  - TGP
  - BH
  - alpha-psi
  - falsification
  - ngEHT
  - LIGO
  - LISA
  - MICROSCOPE
  - NICER
  - Cassini
---

# Phase 3 — Setup: multi-source falsification map of α(ψ)

> **Cel:** Zarejestrować pełną falsifiable predictions map dla α(ψ) =
> α₀(ψ-1)²Θ(ψ-1) z parametrami z Phase 2 (α₀≈4.02, n=2, ψ_th=1) z
> konkretnymi horyzontami eksperymentalnymi i discrimination factors.
> Generate 6 nowych entries BH4–BH9 do PREDICTIONS_REGISTRY.md (analog
> SC.1.Phase3 → SC4–SC7).

---

## Stan pre-Phase 3

Phase 1+2 dostarczyły:
- Path E α(ψ) jako jedyna minimal extension Candidate D (Phase 1)
- ψ_th = 1 DERIVED, n = 2 DERIVED (Phase 2)
- α₀ ≈ 4.02 PARTIALLY DERIVED + STRUCTURAL HINT α₀ = κ_TGP² (0.75% match)
- Multi-source ψ_ph universality across 10 orders of M_BH

Phase 3 nie tworzy nowej teorii — **konsoliduje** Phase 1+2 i Phase 0+
heurystykę w fałsyfikowalne pre-registered predictions z konkretnymi
eksperymentami.

---

## 7 sub-testów Phase 3

### T3.1 — ngEHT photon ring multi-source map

**Cel:** Zarejestrować pełną listę SMBH źródeł dla ngEHT 2030+ z
predykcją uniwersalnego +14.56% shadow diameter shift.

**Sources (z masami z literatury):**
- SgrA* (4.3·10⁶ M_⊙)
- M87* (6.5·10⁹ M_⊙)
- NGC 1277 SMBH (1.7·10¹⁰ M_⊙)
- Cen A / NGC 5128 (5.5·10⁷ M_⊙)
- NGC 4258 (4·10⁷ M_⊙)
- Sombrero Galaxy / M104 (1·10⁹ M_⊙)
- IC 1101 (4·10¹⁰ M_⊙)
- M84 (1.5·10⁹ M_⊙)
- M81 (7·10⁷ M_⊙)
- TON 618 (6.6·10¹⁰ M_⊙ — most massive known)

**Test predykcji:**
- Universal r_ph^TGP/M = 3.88 (vs GR 3.0)
- Shadow diameter shift +14.56%: D_TGP = 1.1456 × D_GR
- ngEHT angular diameter resolution ~ 5-10 μas
- For each source: predict shadow size from M, distance, plus +14.56% TGP signature

**Falsification:** dwóch lub więcej sources daje **różne** shifty (>5% spread) → α(ψ_ph) nie jest universal → α(ψ) wymaga Phase 1 PLAN paper-rigor revision.

### T3.2 — LIGO O5 / LISA ringdown frequency shift

**Cel:** Predykcja phase-shift w QNM ringdown dla GW mergerów, dla M_f
od stellar (10 M_⊙) do supermassive (10⁷ M_⊙).

**Test:**
- ψ_ringdown ~ 1.20 (estimate, post-merger near-horizon strong-field)
- α(ψ_ringdown) = α₀ × (0.20)² = 4.02 × 0.04 = 0.161
- Frequency shift: δf/f ~ α(ψ) × O(geometry) ~ 1-5%
- LIGO O5 sensitivity (high-SNR events): ~0.5% on f, ~1% on damping τ
- LISA SMBH mergers: ~0.1% precision

**Predicted observables:**
- GW150914-class (M_f = 65 M_⊙): δf ~ 1-3% from GR
- LISA SMBH (M_f = 10⁶ M_⊙): δf ~ same fractional (universal w geom)

**Falsification:** post-O5 high-SNR ringdowns z δf consistent z 0 (within 0.1%) → α(ψ) coupling falsified at ringdown, OR ψ_ringdown estimate wrong.

### T3.3 — NICER pulsar M-R relation

**Cel:** Predykcja modyfikacji NS M-R curve od α(ψ) coupling przy
surface ψ ~ 1.3-1.4.

**Test:**
- NS surface ψ z compactness 2GM/(c²R):
  - 1.4 M_⊙ NS, R ≈ 12 km: 2M/R ≈ 0.41 → ψ ≈ 1.31
  - 2.08 M_⊙ NS (J0740), R ≈ 13.7 km: 2M/R ≈ 0.44 → ψ ≈ 1.34
- α(ψ_NS) = α₀ × (ψ_NS - 1)² ~ 0.4 (significant)
- Effective TOV equations modified; M-R curve shifted by ~1-3%

**NICER measurements:**
- J0030+0451: M = 1.44 M_⊙, R = 13.0 km
- J0740+6620: M = 2.08 M_⊙, R = 13.7 km

**Falsification:** TGP M-R deviation z GR przewyższa NICER precision (~3%) → albo TGP α(ψ) coupling przy NS-scale falsified, albo EOS uncertainty maskuje sygnał.

### T3.4 — MICROSCOPE-2 WEP test

**Cel:** Re-rejestracja Phase 2 T2.2 prediction jako konkretne pre-registered
falsification dla MICROSCOPE-2 mission (~2030).

**Test:**
- η_TGP = α₀ × (ψ_Earth - 1)² = 4.02 × (7·10⁻¹⁰)² = **1.97·10⁻¹⁸**
- MICROSCOPE-2 sensitivity ~ 10⁻¹⁷
- Margin 5.1× below detection threshold

**Falsification:** MICROSCOPE-2 wykrywa η > 10⁻¹⁷ → TGP n=2 falsified
(η > predicted by 5×); n=3 wymagane; lub Path E α(ψ) functional form
falsified.

### T3.5 — Solar PPN Cassini-class precision

**Cel:** Zarejestrować TGP α(ψ) prediction dla Solar PPN γ parameter w
Cassini bound (γ-1 < 2.3·10⁻⁵) i future precision (~10⁻¹⁰).

**Test:**
- ψ_Sun surface = 1 + GM_⊙/(c²R_⊙) = 1 + 2.12·10⁻⁶
- α(ψ_Sun) = α₀ × (2.12·10⁻⁶)² = 4.02 × 4.49·10⁻¹² = **1.81·10⁻¹¹**
- Cassini bound γ-1 < 2.3·10⁻⁵ → TGP at 1.8·10⁻¹¹ → margin **10⁶**

**Future precision missions** (LATOR, BEACON, post-Cassini PPN):
- Targeted precision γ-1 < 10⁻¹⁰ → TGP signature still 10× below

**Falsification:** future Solar PPN test wykrywa γ-1 > 10⁻¹⁰ niezgodne z TGP α(ψ_Sun) → coupling form czy ψ_Sun calibration falsified.

### T3.6 — Cross-sector consistency falsification design

**Cel:** Wykorzystując Phase 2 T2.5 hint (α₀ = κ_TGP² match 0.75%),
zaprojektować eksperymentalny test cross-sector identity z konkretnymi
precyzjami.

**Test design:**
- ngEHT precision na photon ring shift: 1-5% (high-quality SMBH targets)
- → constrain α₀ to 1-5% precision
- TGP-SC v2 calibration of κ_TGP: ~0.5% from V/Nb/Ta/Mo/Pd ensemble
- → κ_TGP² known to ~1% precision
- Combined: |α₀ - κ_TGP²| measurable to ~3% precision

**Falsification:** if combined measurement gives |α₀ - κ_TGP²|/κ_TGP² > 5%
→ cross-sector unification hypothesis rejected; α₀ pozostaje z BH sektora
calibrated; κ_TGP pozostaje z SC sektora calibrated.

**Confirmation:** if combined measurement gives match < 1% → strong
support dla √α₀ = κ_TGP jako fundamental TGP charge.

### T3.7 — PREDICTIONS_REGISTRY entries BH4–BH9

**Cel:** Generate 6 new BHx predictions z Phase 3 sub-tests, gotowe do
inklucji w PREDICTIONS_REGISTRY.md.

| Entry | Sector | Source | Prediction | Horizon |
|-------|--------|--------|------------|---------|
| **BH4** | Photon ring | ngEHT multi-SMBH | Universal +14.56% shadow diameter, 10 sources, 10.1 orders M_BH | 2030+ |
| **BH5** | Ringdown | LIGO O5 / LISA | δf/f ~ 1-5% w QNM frequency post-merger | 2027+ (LIGO) / 2035+ (LISA) |
| **BH6** | NS M-R | NICER | M-R curve shift ~1-3% from GR | 2027+ (NICER, NICER+) |
| **BH7** | WEP | MICROSCOPE-2 | η_TGP = (2±1)·10⁻¹⁸ | 2030+ |
| **BH8** | Cross-sector | ngEHT + TGP-SC | α₀ = κ_TGP² (0.75% match, 5% test) | 2030+ |
| **BH9** | Solar PPN | LATOR/BEACON | γ-1 ~ 10⁻¹¹ (10× below 10⁻¹⁰ precision) | 2035+ |

**Status każdego:** PRE-REGISTERED (Phase 3, 2026-04-28), LIVE.

---

## Materiał wykonawczy

- **Skrypt:** `phase3_multi_source_falsification.py` (7 sub-tests)
- **Output:** `phase3_multi_source_falsification.txt`
- **Memo:** `Phase3_results.md` (closure + verdict)

## Środowisko

```bash
cd TGP/TGP_v1/research/op-bh-alpha-threshold
PYTHONIOENCODING=utf-8 python -X utf8 phase3_multi_source_falsification.py 2>&1 | tee phase3_multi_source_falsification.txt
```

---

## Constants used

```
alpha_0    = 4.02       (Phase 2 strict, xi=1)
n          = 2          (Phase 2 derived)
psi_th     = 1          (Phase 2 derived)
psi_ph     = 1.168      (M9.2-D universal photon ring)
psi_ringdown ~ 1.20     (estimate, post-merger near-horizon)
psi_NS_surf ~ 1.31      (1.4 M_sun NS canonical, 12 km radius)
psi_Earth  = 1 + 7e-10  (Earth surface gravitational potential)
psi_Sun    = 1 + 2.12e-6 (Sun surface gravitational potential)
kappa_TGP  = 2.012      (TGP-SC v2)

Experimental sensitivities:
  ngEHT 2030+:        5-10 μas (shadow angular diameter)
  LIGO O5 (2027+):    ~0.5% f, ~1% τ ringdown
  LISA (2035+):       ~0.1% f, ~0.5% τ
  NICER+ (2027+):     ~3% on M-R curve
  MICROSCOPE-2:       ~10^-17 η_WEP sensitivity
  Cassini bound:      γ-1 < 2.3·10⁻⁵
  Future PPN:         γ-1 < 10⁻¹⁰
```

---

## Cross-references

- [`program.md`](program.md) — overall 3-phase plan
- [`Phase1_results.md`](Phase1_results.md) — H₀ rejected
- [`Phase2_results.md`](Phase2_results.md) — T-α upgraded, cross-sector hint
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — to be updated z BH4–BH9
- [`../op-sc-alpha-origin/Phase3_results.md`](../op-sc-alpha-origin/Phase3_results.md) — wzorzec mini-cyklu

---

## Decyzja po Phase 3

- Wszystkie 7 sub-testów PASS → Phase 3 CLOSED, 6 nowych BHx predictions
  zarejestrowanych, **BH.1 program (3 phases) END**.
- BH.1 program zamknięty paralelnie do SC.1: 5+7+7 = 19 cumulative tests.
- Open issues dla full Phase 1 PLAN: 5 (z Phase 2 T2.7) + experimental
  wait na ngEHT/LIGO O5+/MICROSCOPE-2/NICER (2027–2035 horizon).
