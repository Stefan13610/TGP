---
title: "ADDENDUM 2026-05-10 — B9 WEP MICROSCOPE composition test, native-first reframe"
date: 2026-05-10
parent: "[[B9_wep_microscope_composition_results.md]]"
type: addendum
addendum_kind: methodological-reframing
modifies_B9_status: NO
modifies_eta_TGP: NO
modifies_sympy_LOCK: NO
status: 🟢 ACTIVE — interpretive overlay (lightweight three-layer wrapper)
related:
  - "[[../../meta/PPN_AS_PROJECTION.md]] (parent methodology, binding 2026-05-10+)"
  - "[[B9_wep_microscope_composition_results.md]] (target of three-layer layering)"
  - "[[../op-L01-rho-stress-energy-bridge-2026-05-04]] (B9 jest L1 input dla L01 N3 SPARC consistency)"
  - "[[../op-L01-rho-stress-energy-bridge-2026-05-04/ADDENDUM_2026-05-10_native_observables_first.md]] §2.2 F8 (B9 6/6 PASS jako L3 falsifier dla ρ_Dirac universality)"
  - "[[M9_1_pp_P1_results.md]] (M9.1'' P1 c_n derivation)"
  - "[[M9_2_results.md]] (M9.2 single-source m_field universality)"
tags:
  - addendum
  - methodological-reframing
  - native-observables-first
  - WEP-test
  - MICROSCOPE
  - composition-test
  - B9
  - lightweight-wrapper
---

# ADDENDUM 2026-05-10 — B9 WEP composition native-first wrapper

## §0 — Status uwagi

Ten addendum **NIE zmienia** B9 CLOSED 6/6 PASS verdict, **NIE zmienia**
sympy LOCK `m_field ~ (qM)²/(4π·σ)`, **NIE zmienia** η_TGP_lab values
(1.32×10⁻²⁶ without T-α, 6.49×10⁻⁴⁵ with T-α). Co addendum **dodaje**:

1. **Three-layer presentation** (L1 / L2 / L3) per `meta/PPN_AS_PROJECTION.md`
   §3.1 mandatory binding 2026-05-10+. B9 results są technically detailed ale
   *implicit* w 3-layer terms — addendum czyni je explicit.
2. **L1 link upstream** do L01 ρ-bridge cycle (B9 jako **L3 falsifier**
   dla L01 F8 native ρ_Dirac universality claim).
3. **Native parameter audit** + forced-zero declarations explicit.
4. **Lightweight wrapper** — najmniejszy scope wśród propagation cycles
   2026-05-10, bo B9 sam w sobie jest tightly closed i numerically clean.

## §1 — Three-layer specification dla B9

### §1.1 — L1 (native predictions)

**Native source structure (z L01 ρ-bridge cycle F8):**

```
ρ_Dirac = m·|ψ|²/c_0²    [native, z T^μ_μ_Dirac = m·ψ̄ψ on-shell]
```

`ρ_Dirac` jest *uniwersalna* dla wszystkich masywnych fermionów (e, μ, τ, kwarki) —
**universality of free fall** w klasycznej granicy jest natywną własnością
TGP framework przez `ax:metric-coupling`.

**Native sub-leading correction (B9 main result):**

W TGP fenomenologii, deviation od WEP nie jest *zerowa* — jest **suppressed
przez (qM)²/σ scaling** plus dodatkowe T-α suppression przy lab regimes:

```
m_field ~ (qM)²/(4π σ)                              [B9 sympy LOCK, weak-field]
m_field/m_g ~ q · U_surface                         [linear w q i σ⁻¹]
η_AB ~ |δq/q| + |δσ/σ|                              [LEADING composition effect]
η_TGP_lab(Pt vs Ti) = 1.32·10⁻²⁶                    [no T-α suppression]
η_TGP_lab(Pt vs Ti) = 6.49·10⁻⁴⁵                    [with T-α suppression]
```

To są **native predykcje** — wynikają strukturalnie z M9.2 Yukawa equation +
two-component composition analysis na BVP.

### §1.2 — L2 (projection charts)

| Native quantity | Phenomenological chart | Status |
|---|---|---|
| `ρ_Dirac` universality | PPN η chart (Eötvös ratio) | η_TGP_native = 0 (universal) + η_geom from σ-mismatch |
| `m_field/m_g ~ q·U_surface` | 5-th force chart (Yukawa-like, scale 1/√β) | suppressed by σ-coupling, not Yukawa-like |
| `η_AB ~ |δq/q| + |δσ/σ|` | composition-dependent η chart | tested by Pt vs Ti pair |
| MICROSCOPE Pt vs Ti ratio | lab projection (U_lab ~ 2·10⁻²⁶) | full spectroscopic + geometric mass calc |

**Multi-coefficient signature pattern (TGP-specific):**

B9 daje *trzy* independent composition-tests (geometric, coupling, inhomogeneous ρ).
Wszystkie trzy mają specific TGP signature (linear w δq/q, quadratic w geometric
σ-mismatch, structure-sensitive dla ρ-inhomogeneity). To jest **unique TGP
signature** distinguishing od generic 5-th force theories.

### §1.3 — L3 (falsification map)

| Bound | Native coef constrained | Window | Status |
|---|---|---|---|
| **MICROSCOPE 2017** η ≤ 1.1·10⁻¹⁵ | ρ_Dirac universality + η_TGP_lab | 8.3·10¹⁰× safe (no T-α), 1.7·10²⁹× safe (with T-α) | **PASS** ✓ |
| **MICROSCOPE-2 future** η ≤ 10⁻¹⁷ | Same coefs, tighter | 7.6·10⁸× safe (no T-α), 1.5·10²⁷× safe (with T-α) | **PASS** projected |
| Eöt-Wash 5-th force | δq atomic-content variation | broad (composition tests) | PASS automatic |
| LLR Nordtvedt η < 4.4·10⁻⁴ | Universal coupling (η_TGP=0 in classical) | broad | PASS automatic |

**Kluczowy native test:** B9 6/6 PASS *demonstruje empirycznie* że TGP `ρ_Dirac
= m·|ψ|²/c_0²` universality jest preserved przy current sensitivity (10¹⁰-10²⁹×
safety margin). To jest **direct L3 falsifier** dla L01 F8 native claim.

## §2 — Upstream/downstream alignment

### §2.1 — B9 jako L3 falsifier dla L01 F8

L01 ADDENDUM §2.2 F8 (Dirac fermion ρ-mapping) explicite cytuje B9 jako
**L3 falsifier dla native ρ_Dirac universality**:

> (a) Eöt-Wash, MICROSCOPE: konstrukcjonalnie testują `δρ_Dirac` między atomami.
> PASS (B9 6/6).

Ten addendum zwraca cross-link explicite z B9 strony — B9 results są
**operacyjną demonstracją** L01 F8 native claim w lab regime.

### §2.2 — B9 jako input dla potencjalnego op-SPARC-rho-consistency cycle (N3)

L01 NEEDS N3 (SPARC consistency) wymaga weryfikacji `ρ_SPARC ≡ ρ_baryon ≡
-T^μ_μ_dust/c_0²` dla galactic rotation curves. B9 result `ρ_Dirac = m·|ψ|²/c_0²`
jest **L1 native source** dla SPARC analiza (galactic stars + HI gas są
non-relativistic Dirac fermions).

Native form alignment: B9 verifies non-relativistic Dirac ρ w *lab* regime;
SPARC tested w *galactic* regime; **identical** ρ formula expected. To powinien
być explicit w przyszłym N3 cycle.

### §2.3 — Multi-coefficient signature jako TGP-distinguishable

B9 daje **3 independent composition tests** (geometric σ, coupling δq, inhomogeneous
ρ). Per `meta/PPN_AS_PROJECTION.md §3.1` falsification map, multi-test pattern
distinguishes TGP od:

- **Standard scalar 5-th force** (single Yukawa coupling) — B9 ma specific
  scale `m_field/m_g ~ q·U_surface`, NIE Yukawa-like exp(-r/λ)
- **MOND / TeVeS** (acceleration-scale modifications) — B9 jest universal
  w acceleration, lokalna w composition
- **Brans-Dicke dilaton** (Bekenstein/Sandvik varying-α) — B9 NIE wprowadza
  varying gauge coupling (per ψ.1 ADDENDUM §3 sterile sector exclusion)

B9 + ψ.1 + L01 razem dają **TGP-specific multi-test signature** — pattern
constraints są *zwieloraźnie konstrainowane*, nie pojedynczo.

## §3 — Native parameter audit (B9)

Per `meta/PPN_AS_PROJECTION.md §3.3`:

```
Independent native parameters constrained by B9 cycle:
  - q (TGP coupling per mass)        [= G_0 mapping z Newton limit; structural]
  - σ (object size/density profile)   [composition observable; not free param]
  - ρ_Dirac universality (m·|ψ|²/c_0²) [native L01 F8, B9 demonstrates empirically]

Free coefs (deferred):
  - T-α suppression factor exact value     → depends on local Φ regime
  - Non-leading composition corrections     → higher-order (negligible at MICROSCOPE)

Forced from substrate symmetry:
  - WEP universality classical (η_TGP_native = 0 dla universal ρ_Dirac)
  - All atoms couple universally via ρ = m·|ψ|²/c_0² (FORCED z L01 F8)
  - GW170817 c_GW = c_EM (FORCED z single g_eff, irrelevantne dla B9 ale framework consistency)
  - α_i ≡ 0 (preferred frame, FORCED z Lorentz-invariance; irrelevantne dla B9)

Native parameter count for B9: 0 free (all derived z L01 + M9.2 sympy LOCK).
                              B9 sondzi *empirically* native universality, 
                              z η_TGP_lab numerical extrapolation.
```

## §4 — Co addendum NIE zmienia

- **B9 CLOSED 6/6 PASS verdict:** unchanged.
- **Sympy LOCK m_field ~ (qM)²/(4π σ):** unchanged.
- **η_TGP_lab(Pt vs Ti) = 1.32·10⁻²⁶ (no T-α), 6.49·10⁻⁴⁵ (with T-α):** unchanged.
- **MICROSCOPE 2017 + MICROSCOPE-2 future margins:** unchanged.
- **6/6 PASS test sequence (sympy LOCK, geometric, coupling, inhomogeneous,
  realistic, T-α):** unchanged.
- **Predecessor M9_2_results.md single-source verification:** unchanged.

Addendum służy jako:
- *lightweight three-layer wrapper* dla B9 technical results
- *cross-link explicit* z L01 F8 (downstream consumer L01 → upstream demonstrator B9)
- *multi-test signature documentation* dla TGP-distinguishable composition tests
- *binding methodology application* dla future cycles citing B9

## §5 — Sign-off

**Addendum authored:** 2026-05-10 (Claudian, kontynuacja propagacji native-first
methodology, lightweight wrapper jako trzeci w trio T01 + Φ-vacuum + B9).

**Status:** ACTIVE interpretive overlay. B9 classification preserved (CLOSED
6/6 PASS).

**Cross-references:**
- `meta/PPN_AS_PROJECTION.md` — parent methodology (binding 2026-05-10+)
- `op-L01-rho-stress-energy-bridge-2026-05-04/ADDENDUM_2026-05-10` §2.2 F8 —
  L01 native claim demonstrated empirically by B9
- `B9_wep_microscope_composition_results.md` — target of three-layer layering
- `M9_1_pp_P1_results.md` — M9.1'' P1 c_n derivation (parent of T01 framework)
- `M9_2_results.md` — single-source m_field universality (B9 predecessor)
- `audyt/T01_LIGO3G_falsifier/ADDENDUM_2026-05-10` — sibling addendum (gravity
  cycle, GW regime; B9 covers WEP regime)
- `op-Phi-vacuum-scale-2026-05-09/ADDENDUM_2026-05-10` — sibling addendum
  (vacuum sector; B9 covers matter sector via L01 F8)

**Insight credit:** autor TGP (γ vs β natywność insight 2026-05-10) + 5-fold
cross-cycle convergence diagnostic now extended to **gravity-WEP-matter triad**
(T01 + Φ-vacuum + B9 propagation 2026-05-10).
