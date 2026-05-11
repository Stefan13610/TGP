---
title: "op-L01-N3-SPARC-rho-consistency — verification że ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0² w non-relativistic galactic limit (<1% precision)"
date: 2026-05-11
type: research-cycle
status: 🟡 OPEN — Phase 0+1 in single session (low-priority cosmetic)
parent: "[[../op-L01-rho-stress-energy-bridge-2026-05-04/README.md]]"
predecessors:
  - "[[../op-L01-rho-stress-energy-bridge-2026-05-04/]] (CLOSED-DERIVED 2026-05-10, defers N3 here)"
  - "[[../op-L01-N1-EM-trace-anomaly-TGP-2026-05-11/]] (sister N1 closed)"
  - "[[../op-L01-N2-QCD-trace-anomaly-cosmology-2026-05-11/]] (sister N2 closed)"
  - "[[../op-emergent-metric-from-interaction-2026-05-09/]] (g_eff[{Φ_i}] formalism)"
related:
  - "[[../galaxy_scaling/]] (SPARC fits w rotation curves)"
  - "[[../nbody/]] (N-body symulacje TGP)"
classification: VERIFICATION_CYCLE_LOW_PRIORITY
priority: cosmetic / cleanup
goal: "Verify że ρ_SPARC = ρ_baryon (HI + stars + bulge) ≡ -T^μ_μ_dust/c_0² w non-relativistic galactic limit do <1% precision; sanity check double-counting vs TGP-emergent DM (z g_eff[Φ̄] background); update nbody/ documentation z explicit T^μ_μ → ρ mapping."
estimated_effort: "1 sesja (compact cycle)"
target_window: "L1 native: ρ_TGP = -T^μ_μ_fluid/c_0² for dust limit (p=0) gives ρ_TGP = ρ_rest exactly; for galactic stars v² ~ 4·10⁻⁷ c² gives correction < 10⁻⁶ ≪ 1%."
six_requirements_target:
  - "P1: Sympy LOCK ρ_TGP = ρ_rest in dust limit (p=0)"
  - "P2: Numerical bound: galactic v ~ 200 km/s → v²/c² ~ 4·10⁻⁷; correction (1 - v²/(2c²)) < 10⁻⁶"
  - "P3: HI gas thermal v ~ 1 km/s → utterly negligible correction"
  - "P4: No double-counting vs TGP-emergent DM (g_eff[Φ̄] background source)"
  - "P5: nbody/ + galaxy_scaling/ documentation update z explicit T^μ_μ → ρ mapping"
  - "P6: SPARC residuals chi^2 unchanged (verification, not new fit)"
risk_flags:
  - "R1: Double-counting risk — gdyby SPARC fits used both ρ_baryon AND ρ_DM_separate, byłaby naruszenie S05 (single fundamental field)"
  - "R2: Relativistic correction at galactic centers — near SMBH (Sgr A*) v can approach c, but this is outside SPARC scope (galactic dynamics far from BH)"
tags:
  - L01
  - L01-N3
  - SPARC-consistency
  - dust-limit-verification
  - non-relativistic-galactic
  - low-priority-cosmetic
  - cycle-open-2026-05-11
---

# op-L01-N3-SPARC-rho-consistency-2026-05-11

> **Cel:** zamknąć **N3** z [[../op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]]
> — verify że SPARC fits (galaxy rotation curves) używają konsystentnej ρ
> w sensie L01 framework: `ρ_SPARC ≡ ρ_baryon ≡ -T^μ_μ_dust/c_0²` non-relativistic
> galactic limit do <1% precision.

## Geneza

L01 NEEDS §N3 (low-priority cosmetic):

> N-body symulacje TGP (SPARC galaxy rotation, `nbody/`) używają `ρ = ρ_baryon`
> (HI + stars + bulge) bez explicit T^μ_μ derivation. Pytanie: czy `ρ_baryon`
> w SPARC is the same as `ρ = -T^μ_μ/c_0²` in non-relativistic limit?

**Pre-existing analysis (L01 NEEDS §T.3):**
```
T^00 = ρ_rest · c² + (1/2)·ρ_rest·v²  ≈  ρ_rest·c²   (v² ≪ c²)
T^ii = ρ_rest·v²  ≈  0
T^μ_μ = T^00 - T^ii  ≈  -ρ_rest·c²
ρ = -T^μ_μ/c² = ρ_rest    ✓
```

Klasycznie konsystentne. Tego cyklu zadanie: **explicit numerical verification**
z bound check + double-counting sanity check + nbody documentation note.

## Centralna hipoteza H1

**H1:** Dla galactic regime (stars + HI gas + bulge), ρ_SPARC = ρ_baryon
identyfikuje się z `ρ_TGP = -T^μ_μ_dust/c_0²` z dokładnością far better than 1%
(typowo < 10⁻⁶ relativistic correction).

## Six P-requirements

| # | Requirement | Notes |
|---|-----------|-------|
| **P1** | Sympy LOCK ρ_TGP = ρ_rest in dust limit (p=0) | analytic |
| **P2** | Galactic stars v ~ 200 km/s: correction < 10⁻⁶ | numerical |
| **P3** | HI gas v ~ 1 km/s: correction utterly negligible | numerical |
| **P4** | No double-counting vs TGP-emergent DM | structural argument |
| **P5** | nbody/+galaxy_scaling/ docs updated z mapping | cosmetic |
| **P6** | SPARC residuals preserved (verification, NOT new fit) | structural |

## Risks

R1 (double-counting), R2 (galactic-center relativistic limit).

## Pliki w cyklu (compact)

| Plik | Status | Opis |
|------|--------|------|
| [[README.md]] | ✅ ten plik | overview |
| [[Phase0_balance.md]] | 🟡 next | balance + 6/6 gate |
| [[Phase1_sympy.py]] / [[Phase1_sympy.txt]] | next | dust-limit + bounds |
| [[Phase1_results.md]] | next | verification + double-counting |
| [[Phase_FINAL_close.md]] | next | sign-off |
| [[FINDINGS.md]] | next | exportable |
| [[NEEDS.md]] | next | residual |

## Probability assessment

| Outcome | Prob | Rationale |
|---|---|---|
| Pełen DERIVED | 80-90% | Trivial dust-limit derivation; existing L01 analysis already 90% complete |
| STRUCTURAL CONDITIONAL | 5-15% | gdyby double-counting issue zostałaby found w nbody code |
| STRUCTURAL_NO_GO | <5% | very unlikely |

## Connection do innych cykli

- **L01 ρ-bridge** (CLOSED-DERIVED): zamyka N3; uses §T.3 pre-existing analysis
- **op-emergent-metric** (STRUCTURAL_DERIVED): TGP-emergent DM source z `g_eff[Φ̄]`
  background — separate od ρ_baryon source
- **galaxy_scaling cycles** (gs10-gs61): SPARC fits używają ρ_baryon, w niniejszej
  weryfikacji potwierdzane consistent z TGP framework

## Status

🟡 **OPEN — Phase 0+1 in single session** (low-priority cosmetic per L01 NEEDS).

---

**Cycle opened:** 2026-05-11 (Claudian, post-N1+N2 closure same day; quickest
remaining L01 closure).
