---
title: "FINDINGS — op-FRW-radiation-era-varying-c (post-Phase-1)"
date: 2026-05-07
parent: "[[README.md]]"
type: findings
status: PHASE_1_COMPLETE — 4/5 PASS, critical binary outcome ujawniony
tgp_owner: research/op-FRW-radiation-era-varying-c-2026-05-06
tags:
  - findings
  - EXT-1
  - phase1-complete
  - critical-finding
---

# FINDINGS — post-Phase-1

## Status

**PHASE 1 EXECUTED 2026-05-07** (4/5 PASS). Phase 2 ENABLED z caveat —
F1.4 ujawnia **binary outcome bifurcation**.

## Phase 1 highlights

### ✅ F1.1 — FRW background z varying c(Φ)

```
ds² = -c_0² · (Φ_0/Φ) · dt² + a(t)²[dr² + r²dΩ²]
√(-g) = c(Φ) · a³ · r² · sin θ
```

Signature (-,+,+,+) zachowane. Smooth dla Φ > 0.

### ✅ F1.2 — Friedmann eq w TGP (sympy LOCK)

```
H_TGP² = (8πG_0/3) · (Φ_0/Φ) · [ρ_matter/c_0² + (1/2)K(Φ)(dΦ/dτ)² + V(Φ)]
```

Recovery w ψ=1 limit: H_TGP² → H_GR² (sympy diff = 0).

Cross-check vs M10.1 V_0 = β/12 ≈ Ω_DE0 = 0.685 ✓.

### ✅ F1.3 — ax:c-ax:G consistency

**Krytyczne invariants:**
- G/c² = G_0/c_0² (Schwarzschild scale **invariant** pod Φ)
- G/(ℏc) = G_0/(ℏ_0 c_0) (constant)
- ℓ_Pl = √(ℏG/c³) ∝ ψ^(1/4) (varies as Φ-quarter-power)

### ⚠ F1.4 — Limity asymptotyczne (BIFURCATION CRITICAL)

**m_eff = √(γ/K_geo) ~ H_0** (z T-Λ closure 2026-04-26).

**W erze radiacyjnej:** H(z) >> m_eff dla wszystkich z > 0 → ψ frozen ≈ 1.

**Scenariusz (a) — DOMINANT:** ψ frozen → varying constants efektywnie
wyłączone w erze radiacyjnej. ALE TGP NIE ma ρ_rad strukturalnie
(T^μ_μ_EM = 0). Wynik:

```
H_TGP(z=10⁹) / H_GR(z=10⁹) ≈ √(z_eq / (1+z)) = √(3400 / 10⁹) ≈ 0.0018 = 0.18%
```

**Drift ~99.8% << 5% gate** → BBN ⁴He **CATASTROPHIC FAIL**.

**Scenariusz (b)** — wymaga ψ_BBN ≈ 3·10⁻⁶ (5 orders of magnitude poniżej
today). Wymaga **inicjalnego warunku** ψ << 1 z mechanizmu **przed**
TGP_v1 (np. inflacyjna faza). M9.1'' łamie założenia perturbacyjne.

### ✅ F1.5 — Phase 1 GATE 4/5 PASS (z caveat)

Phase 2 ENABLED. **Krytyczna rekomendacja:** numerical Φ-EOM solver
w Phase 2 powinien szybko potwierdzić/odrzucić scenariusz (a).

## Critical finding (Phase 1 main outcome)

> **EXT-1 ścieżka A jest PRAWDOPODOBNIE NIE DO URATOWANIA w obecnym
> TGP_v1 framework.**
>
> Najprawdopodobniejszy outcome (P ≈ 75-85%): scenariusz (a) confirmed
> w Phase 2 numerical → H_TGP(z=10⁹) ≈ 0.18% H_GR → BBN MASSIVE FAIL
> → ścieżka A FAILS.
>
> Wymagana decyzja autora (NEEDS N9):
> - **Ścieżka D** (L_mat extension dla pól cechowania) — narusza S04
>   ZAMKNIĘTY 2026-05-04 przez B9 MICROSCOPE 6/6 PASS. Re-open
>   S04 audit konieczny.
> - **Ścieżka E** (przyznanie zakresu post-recombination) — TGP
>   jako "GR w limicie post-recombination", drastycznie obniża rangę.

## Subiektywna ocena post-Phase-1

| Outcome | Pre-Phase-1 (z EXT-1 v2) | Post-Phase-1 |
|---------|--------------------------|--------------|
| P(ścieżka A → DERIVED) | 35-45% | **5-10%** |
| P(ścieżka A → STRUCTURAL CONDITIONAL) | 30-40% | **10-15%** |
| P(ścieżka A → FAIL → pivot D/E) | 25-35% | **75-85%** |

**Phase 1 obniżyła subiektywną ocenę szans ścieżki A** z 65-85%
(P(DERIVED) + P(STRUCTURAL) = 65-85%) do **15-25%**.

## Recommended next steps

1. **Phase 2 short-cycle** (1-2 sub-tests):
   - F2.1 Numerical Φ-EOM solver dla ψ(z), z ∈ [10³, 10¹⁰]
   - F2.2 Confirmation scenariusza (a) — ψ frozen ≈ 1 → H_TGP ≈ 0.18% H_GR
   - Jeśli (a) confirmed: STOP, Phase 2 verdict "ścieżka A FAILS"

2. **Decision report do user-a:**
   - Ścieżka D vs E vs novel approach
   - Implications dla S04 (B9 MICROSCOPE preserved jeśli E; re-open jeśli D)
   - Implications dla TGP rangę (E = "GR z limit post-rec")

3. **Alternative path opening:**
   - Jeśli decyzja: ścieżka D → otwórz cykl `op-Lmat-extension-S04-reopen/`
   - Jeśli decyzja: ścieżka E → update TGP_FOUNDATIONS § scope statement

## Cross-references

- [[Phase1_results.md]] — pełne sub-test results
- [[Phase0_balance.md]] — pre-derivation balance sheet
- [[README.md]] — program plan
- [[NEEDS.md]] — N1-N9 open questions (wszystkie nadal OPEN)
- [[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2
- [[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]
- [[../../audyt/S04_metric_coupling_axiom/]] — S04 audit (preservation
  check w ścieżce E vs re-open w ścieżce D)
