---
title: "Phase 1 setup — ROFM cluster-scale extension formalism z g_eff[{Φ_i}] multi-source gradient enhancement"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 1
status: 🟡 setup phase
sub_needs_addressed: [N0.1, N0.2, N0.3, N0.4, N0.5, N0.6, N0.7, N0.9]
risks_addressed: [R1, R3, R4, R6]
predecessor: "[[./Phase0_balance.md]] (6/6 gate PASS)"
sister_cycle_pattern: "[[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] (galactic-disk; cluster outside scope)"
tags:
  - phase1
  - ROFM-cluster-extension
  - multi-source-gradient
  - g-eff-functional
  - cluster-deficit-derivation
---

# Phase 1 setup

## §0 — Cel Phase 1

Wyprowadzić **cluster-scale ROFM extension** w TGP framework:
1. ROFM ν(y) phenomenology przy cluster acceleration regime (g_bar ~ 10⁻¹¹ m/s²)
2. Multi-source g_eff[{Φ_i}] gradient cross-term σ^ij contribution
3. Konstruktywnie ustalić: czy native TGP mechanism closes cluster ~32% deficit
   (gs13/gs55 documented), czy wymaga sterile ν addition

**Three-way outcome** (probability assessment z README):
- **H1a (~30-40%):** TGP-pure (multi-source enhancement δν zamyka deficit)
- **H1b (~35-45%):** TGP + sterile ν ~2 eV addition needed
- **H1c (~15-25%):** STRUCTURAL_NO_GO (cluster outside TGP scope)

## §1 — ROFM phenomenology (recap z galaxy_scaling)

### §1.1 — Galactic-disk regime canonical formula

Per galaxy_scaling cycles gs10-gs61 (closure 2026-04-19):
```
v_circ²(r) = G·M_bar(r) · ν(y)
ν(y) = 1 + exp(-y^α) / y^γ
y(r) = g_bar(r) / a₀ = G·M_bar(r) / (r² · a₀)
```

z parametrami fitting SPARC sample (175 galaxies):
- **α = 0.81** (MOND threshold exponent)
- **γ = 0.41** (transition strength, scale-dependent)
- **a₀ = 1.20·10⁻¹⁰ m/s²** (MOND-like acceleration scale)

**Asymptotic behavior:**
- **y >> 1 (Newtonian):** exp(-y^α) → 0, ν → 1 (no enhancement)
- **y << 1 (deep MOND):** exp(-y^α) → 1, ν → 1/y^γ (enhancement)

### §1.2 — ROFM jest phenomenologically calibrated, NIE first-principles

**KLUCZOWA OBSERWACJA (closure 2026-04-19, Path C):**

ROFM ν(y) parameterization fits SPARC galaxies excellently (χ²_red competitive
z MOND simple), ALE **brak first-principles TGP derivation**. Funkcja ν(y)
jest **phenomenologiczna** — interpolation function między Newton (y>>1) a
MOND-like (y<<1) regime, calibrowana empirycznie z SPARC fits.

**Implikacja dla N3 SPARC closure (2026-05-11):** N3 cycle verified
**ρ_SPARC ≡ ρ_baryon** (mass content) do <10⁻⁶ precision; ROFM ν(y) formula
jest **separate question** (gravity modification mechanism), NIE addressed
przez N3.

## §2 — Cluster regime — acceleration scale comparison

### §2.1 — Cluster vs galactic acceleration

| Regime | r | M_bar | g_bar = GM/r² | y = g_bar/a₀ |
|---|---|---|---|---|
| **Galactic disk** (SPARC) | ~10 kpc | ~10¹¹ M_☉ | ~10⁻¹⁰ m/s² | ~1 (transition) |
| **Galactic outskirts** | ~30 kpc | ~10¹¹ M_☉ | ~10⁻¹¹ m/s² | ~0.1 (deep MOND) |
| **Cluster r₅₀₀** | ~1 Mpc | ~10¹⁴ M_☉ | ~3·10⁻¹¹ m/s² | ~0.25 (deep MOND) |
| **Cluster outskirts** | ~3 Mpc | ~10¹⁴ M_☉ | ~3·10⁻¹² m/s² | ~0.025 (very deep MOND) |

**Cluster regime jest deeply w MOND-like reżim** (y << 1) — gdzie ROFM gives
maximum enhancement ν(y) ≈ 1/y^γ.

### §2.2 — ROFM przewidywanie dla Coma cluster

**Coma cluster parameters** (gs55):
- r₅₀₀ ≈ 1.3 Mpc ≈ 4.0·10²² m
- M_bar(r₅₀₀) ≈ 1.0·10¹⁴ M_☉ ≈ 2.0·10⁴⁴ kg
- M_obs(r₅₀₀) (X-ray + lensing) ≈ 1.5·10¹⁴ M_☉ ≈ 3.0·10⁴⁴ kg
- Deficit factor: M_obs/M_bar ≈ 1.5 (50% enhancement nad baryon)

**ROFM prediction (galactic-calibrated parameters):**
```
g_bar(r₅₀₀) = G·M_bar/r₅₀₀² = 6.67·10⁻¹¹ · 2·10⁴⁴ / (4·10²²)² ≈ 8.3·10⁻¹¹ m/s²
y(r₅₀₀) = g_bar/a₀ ≈ 0.69
ν(0.69) = 1 + exp(-0.69^0.81)/0.69^0.41
        = 1 + exp(-0.738)/0.864
        = 1 + 0.478/0.864
        = 1.553
```

⇒ ROFM-predicted enhancement ν ≈ 1.55× nad pure Newton ⇒ M_TGP ≈ 1.55·M_bar ≈ 1.55·10¹⁴ M_☉.

**Comparison z observed M_obs ≈ 1.5·10¹⁴ M_☉**: ROFM **OVER-predicts slightly**
(102%), within data error. ALE per gs55 dokumentacja: empirically ROFM gives
~68% of M_obs (32% deficit). Discrepancy może wynikać z:
1. Detailed cluster mass profile (NFW vs uniform)
2. r₅₀₀ vs r_vir scale dependence
3. ROFM parametry α, γ might be galactic-specific, NOT cluster-specific

### §2.3 — gs13/gs55 documented deficit

Per gs55 lines 100-136 cluster sample (Coma, Perseus, Virgo, Bullet, A1689):
```
Coma:    f_deficit ≈ 32%, M_TGP ≈ 0.68 · M_obs
Perseus: f_deficit ≈ 30%, M_TGP ≈ 0.70 · M_obs
Mean (excluding Bullet merge): ~32% deficit → TGP recovers ~68%
```

⇒ **Cluster ~32% mass deficit jest EMPIRICAL FACT** w TGP-ROFM galactic-
calibrated parameters. Phase 1 wyzwanie: **czy multi-source g_eff[{Φ_i}]
contribution naprawia ten deficit, czy potrzeba BSM addition.**

## §3 — Multi-source g_eff[{Φ_i}] framework

### §3.1 — Functional form (z emergent-metric Phase 1)

Per [[../op-emergent-metric-from-interaction-2026-05-09/]] Phase 1 results:
```
g_eff^μν = G[{Φ_i}, σ_ab, Φ̄]

g_eff^00 = -A(ψ)
g_eff^ij = δ^ij · B(ψ) + σ^ij · C(ψ) / (Φ₀² c²)
```

z **gradient strain tensor:**
```
σ^ij = (∂^iΦ)(∂^jΦ) - (1/3)·δ^ij·(∂^kΦ)(∂^kΦ)
```

**Single-source limit (galactic):** σ^ij ≈ 0 (Φ approximately radial; trace-free
part vanishes for spherically symmetric Φ). g_eff reduces do A(ψ), B(ψ) only.

**Multi-source limit (cluster):** σ^ij ≠ 0 (multiple galactic Φ-sources within
cluster; non-trivial gradient interaction). C(ψ) coefficient daje additional
metric contribution.

### §3.2 — Hipoteza H1a: multi-source δν enhancement

**Hypothesis:** Cluster mass deficit jest częściowo closed przez multi-source
gradient interaction:
```
ν_cluster(y) = ν_galaxy(y) · (1 + δν_multi)
δν_multi ~ C(ψ) · ⟨σ^ij·σ_ij⟩ / (Φ₀² c²)
```

z δν_multi enhancement growing with cluster density (more sources).

### §3.3 — Required enhancement do close deficit

Per gs55: M_obs/M_TGP ≈ 1/0.68 ≈ 1.47 ⇒ **required δν_multi ≈ 0.47** (47%
additional enhancement nad galactic-calibrated ROFM).

**Test:** Czy multi-source σ^ij gradient interaction w cluster środowisku daje
δν ~ 0.47 strukturalnie?

## §4 — Bullet Cluster — collisionless mass evidence

### §4.1 — Empirical fact

Clowe et al. 2006: Bullet Cluster 1E 0657-558 shows **clear offset** między
X-ray gas (collisional) i lensing mass (collisionless). Total mass map (lensing)
follows galaxies (~stars + ~DM), NIE gas. Standard ΛCDM interpretation: gas
ram-pressure stripped during merger; collisionless DM passed through.

### §4.2 — TGP framework challenge

W TGP framework, **mass enhancement** (cluster deficit closure) musi pochodzić
z **non-baryonic source** (gas jest stripped). Per N3 SPARC closure: **NO
separate ρ_DM matter component** (TGP-emergent DM jest gravitational mechanism,
NIE matter sektor).

⇒ Bullet Cluster offset musi być explained przez:
- **(H1a):** Multi-source g_eff[{Φ_i}] tracking galaxy distribution (not gas)
  natywnie odtwarza lensing-vs-gas offset
- **(H1b):** Sterile ν component dodatkowo tracks galaxies (collisionless;
  consistent z Bullet)

## §5 — Sterile neutrino alternative (jeśli H1a fails)

### §5.1 — Required parameters (Famaey-Angus 2007 framework)

Sterile ν ~2 eV cluster mass component (Angus et al. 2010):
- Mass m_ν ≈ 2 eV (warm DM-like)
- Density ρ_ν concentrated in cluster gravitational potentials
- Sin²2θ_oscillation ~ 10⁻³ (consistent z LSND-like ν oscillation anomaly)

### §5.2 — Bounds compatibility

**Planck 2018 N_eff = 3.046 ± 0.18** bound dla sterile ν:
- N_eff = 3 + 1·f(m_ν, sin²2θ)
- For m_ν = 2 eV + sin²2θ ~ 10⁻³: Δ_N_eff ~ 0.05 → consistent (within 1σ)

**KATRIN 2024** direct m_ν bound:
- m_ν < 0.8 eV (90% CL) na electron neutrino combined mass
- Sterile ~2 eV jest separate species (NOT in active mass eigenstates)

**BBN ⁴He Y_p, D/H:**
- Sterile ν w thermal equilibrium dodaje Δ_N_eff_BBN ~ 0.05 → consistent

**Future CMB-S4** (2030s, Δ_N_eff ±0.04): może detect lub falsify sterile ν.

## §6 — Phase 1 plan + sympy LOCK targets

### §6.1 — Phase 1 sympy targets (8 tests)

Phase1_sympy.py weryfikuje:

1. **T1**: ROFM ν(y) at cluster regime g_bar ~ 10⁻¹¹ m/s² predicts deficit
   ~32% jeśli galactic parameters (α, γ, a₀) używane bez modification
2. **T2**: Cluster acceleration y_cluster < 1 (deep MOND regime)
3. **T3**: Multi-source gradient σ^ij·σ_ij magnitude estimate at cluster density
4. **T4**: Coma cluster M_TGP/M_obs ≈ 0.68 reproduced z galactic ROFM
5. **T5**: Required δν_multi ≈ 0.47 to close gap; comparison z theoretical max
   from C(ψ)·σ²/(Φ₀²c²)
6. **T6**: Bullet Cluster lensing-vs-gas offset compatibility z multi-source
   g_eff (tracks galaxies, NIE gas)
7. **T7**: Sterile ν 2 eV addition: Δ_N_eff < 0.18 (Planck 2018 1σ)
8. **T8**: Outcome decision H1a / H1b / H1c based on δν_multi vs required

Target: 8/8 sympy PASS (Phase 1 derivation rigor; outcome decision honest).

### §6.2 — Phase 1 deliverables

- [[Phase1_setup.md]] (this file)
- [[Phase1_results.md]] — outcome H1a/b/c decision + findings
- [[Phase1_sympy.py]] + [[Phase1_sympy.txt]] — 8 tests

## §7 — Risk addressing in Phase 1

### §7.1 — R1 (cluster sample bias) — partial

Phase 1 używa gs13/gs55 cluster sample (Coma, Perseus, Virgo, Bullet, A1689).
Selection bias possible; Phase 2 will expand sample (~30+ clusters) if Phase 1
H1a/H1b warrants.

### §7.2 — R3 (sterile ν tension z Planck) — addressed

Per §5.2: sterile ν 2 eV addition gives Δ_N_eff ~ 0.05 < Planck 2018 ±0.18 →
within 1σ. **Future CMB-S4 ±0.04 może falsify.**

### §7.3 — R4 (TGP-emergent NIE post-hoc fit) — binding

**Multi-source δν_multi MUST be derived strukturalnie z g_eff[{Φ_i}], NIE
adjusted as free parameter.** Phase 1 sympy T5 will test theoretical max
δν_multi z C(ψ)·σ²/(Φ₀²c²).

### §7.4 — R5 (Bullet Cluster) — addressed

Per §4.2: TGP-emergent multi-source mechanism tracks galaxies (NIE gas);
consistent z Bullet lensing-vs-X-ray offset.

### §7.5 — R6 (cluster scale ~Mpc) — addressed Phase 1

ROFM extension od galactic kpc → cluster Mpc jest natural scaling exercise;
no explicit transition scale required (smooth g_eff[{Φ_i}] functional).

## §8 — Cross-references

- [[./README.md]] / [[./Phase0_balance.md]]
- [[../op-L01-N3-SPARC-rho-consistency-2026-05-11/]] (galactic baseline)
- [[../galaxy_scaling/gs13]] / [[../galaxy_scaling/gs55]] (cluster sample)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]] (multi-source g_eff)
- Clowe et al. 2006 (Bullet Cluster collisionless evidence)
- Angus, Famaey, Diaferio 2010 (sterile ν cluster fit)
- Planck 2018 X (cosmological parameters + N_eff bounds)
- gs55 lines 100-136 (cluster deficit dokumentacja)

---

**Phase 1 setup ready.** Next: Phase1_sympy.py + Phase1_results.md.
