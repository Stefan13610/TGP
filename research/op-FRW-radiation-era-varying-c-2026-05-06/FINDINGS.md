---
title: "FINDINGS — op-FRW-radiation-era-varying-c (post-Phase-2 final)"
date: 2026-05-07
parent: "[[README.md]]"
type: findings
status: PHASE_2_COMPLETE_FAIL — ścieżka A FAILS confirmed numerycznie, decyzja D/E/F pending
classification: STRUCTURAL_NO_GO (analog μ.1', ο.2)
tgp_owner: research/op-FRW-radiation-era-varying-c-2026-05-06
tags:
  - findings
  - EXT-1
  - phase2-complete
  - sciezka-A-FAILS
  - decision-pending
  - STRUCTURAL_NO_GO
---

# FINDINGS — post-Phase-2 final

## Status FINAL

**EXT-1 ŚCIEŻKA A — DEFINITIVELY FAILS 2026-05-07** post Phase 2
short-cycle (2/2 PASS).

**Klasyfikacja:** **STRUCTURAL_NO_GO** (analog μ.1' substrate-log, ο.2
Hubble tension z M03 retrofit — Phase 6 honest reporting baseline
preserved).

## Phase 1 + Phase 2 cumulative

### ✅ Phase 1 (4/5 PASS): Analytical scaffold

- F1.1 FRW background z varying c(Φ): PASS
- F1.2 Friedmann eq w TGP (sympy LOCK): PASS
- F1.3 ax:c-ax:G consistency: PASS
- F1.4 Limity asymptotyczne: PARTIAL — binary outcome ujawniony
- F1.5 Phase 1 GATE: 4/5 PASS

### ✅ Phase 2 short-cycle (2/2 PASS): Numerical confirmation

- F2.1 Numerical Φ-EOM analytical: PASS — scenariusz (a) ψ FROZEN confirmed
- F2.2 H_TGP/H_GR scaling: PASS — **0.184% drift** confirmed analitycznie

## Critical findings

### 1. ψ-EOM analytical solution

W TGP framework w erze radiacyjnej **H >> m_eff** dla wszystkich z > 0
→ Hubble friction DOMINATE → ψ frozen ≈ 1 (overdamped harmonic oscillator).

**Source term analysis** ujawnia non-stable equilibrium: V'(ψ) max
przy ψ=2/3 z V'_max = 4γ/27, ale w radiation era source/K_geo
~ 3·10²⁶·γ >> 4γ/27 → brak stable equilibrium → ψ → 0 lub ψ → ∞.

W obu przypadkach M9.1'' łamie założenia + numerical solver diverges.

### 2. H_TGP/H_GR scaling

W scenariuszu (a) — ψ frozen, varying constants efektywnie wyłączone:

```
H_TGP² / H_GR² = ρ_matter / ρ_rad = (z_eq + 1) / (1 + z) ≈ 3400/(1+z)
Dla z = 10⁹: H_TGP / H_GR ≈ 0.184%
```

### 3. BBN ⁴He impact

```
T_freeze_TGP ≈ 0.20 MeV (vs GR 0.7 MeV)
n/p_TGP ≈ 0.00155 (vs GR 0.157)
Y_p_TGP ≈ 0.31% vs PDG 24.5%
```

**Drift ~99% << 5% BBN gate. CATASTROPHIC FAIL.**

## Probability evolution (subiektywna ocena)

| Outcome | Pre-Phase-1 (EXT-1 v2) | Post-Phase-1 | **Post-Phase-2 (FINAL)** |
|---------|------------------------|--------------|---------------------|
| P(ścieżka A → DERIVED) | 35-45% | 5-10% | **<1%** |
| P(ścieżka A → STRUCTURAL CONDITIONAL) | 30-40% | 10-15% | **<2%** |
| P(FAIL → pivot D/E/F) | 25-35% | 75-85% | **>97%** |

## Decision matrix

### Ścieżka E (RECOMMENDED short-term) — scope acknowledgment

**Co:** Update `TGP_FOUNDATIONS.md` § scope:
> "TGP_v1 is consistent with GR for z < z_recombination (≈ 1100). For
> earlier epochs (BBN, inflation), TGP defers to standard cosmology
> until structural extension (radiation coupling) is developed."

**Pros:**
- Honest acknowledgment of scope limitation
- Trivial implementation (~1 paragraph w FOUNDATIONS)
- Preserves S04 closure (B9 MICROSCOPE 6/6 PASS)
- Preserves M9.1'' weak-field PPN, GW170817, dark energy
- Peer-review rzetelność wzmocniona

**Cons:**
- Drastycznie obniża rangę TGP (nie pełna kosmologia)
- TGP staje się "theory of late-time gravity"

**Estymata:** <1 miesiąc.

### Ścieżka D (long-term research-track) — L_mat extension dla pól cechowania

**Co:** Dodać do L_mat sprzężenie dilatonowe dla pól cechowania:
```
L_mat = -(q/Φ_0)·φ·ρ + L_rad(φ, F_μν)
```

**Pros:**
- Potencjalnie ratuje phenomenologię BBN/CMB
- TGP może zachować claim pełnej kosmologii

**Cons:**
- **NARUSZA ax:metric-coupling** (S04 ZAMKNIĘTY 2026-05-04) — wymaga
  RE-OPEN S04
- B9 MICROSCOPE musi być re-verified post-extension
- Wymaga cyklu `op-Lmat-extension-S04-reopen/` (Phase 1+2+3)
- **30-50% probability success**

**Estymata:** 6-12 miesięcy.

### Ścieżka F (speculative) — pre-BBN inflation novel physics

**Co:** Dodać inflacyjną fazę z ψ_init << 1 → relaxes do ψ=1 today.

**Pros:** Może rozwiązać EXT-1 + dodać inflation physics
**Cons:** 12+ miesięcy, 10-25% probability, M9.1'' łamie założenia w ψ<<1
**Estymata:** 12+ miesięcy.

## Recommendation

**Krótkoterminowo:** Ścieżka **E** (scope acknowledgment) — natychmiastowa
honest implementation.

**Długoterminowo (parallel research-track):** Ścieżka **D** jako separate
cycle.

**Ścieżka F:** zarezerwowana jako "if all else fails".

## Status

✅ **Phase 2 short-cycle COMPLETE 2026-05-07.**

**EXT-1 cycle status:** **STRUCTURAL_NO_GO** (analog μ.1' + ο.2 z M03;
Phase 6 honest reporting baseline preserved).

**Decyzja D/E/F pending — w gestii autora.**

## Honest reporting note

Cykl wykonał Phase 1 (analytical scaffold) + Phase 2 short-cycle
(numerical confirmation), zamknął się jako STRUCTURAL_NO_GO **w 1
sesji** (2026-05-06 setup + 2026-05-07 execution + decision report).

**Wzór godny EXT-1 v2 spirit:** explicit acknowledgment niepowodzenia
ścieżki A jest **honest scientific output**, NIE failure. Phase 6
ABSOLUTE BINDING gate enforced.

## Cross-references

- [[Phase2_results.md]] — pełne F2.1+F2.2 z analytical derivation
- [[Phase1_results.md]] — Phase 1 4/5 PASS
- [[Phase0_balance.md]] — pre-derivation balance sheet
- [[README.md]] — program plan
- [[NEEDS.md]] — N1-N9 (większość RESOLVED post-Phase-2; N9 decision pending)
- [[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2
- [[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]
- [[../../audyt/S04_metric_coupling_axiom/]] — preserve w E lub re-open w D
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]] — B9
- [[../../TGP_FOUNDATIONS.md]] — scope statement target dla ścieżki E
