---
title: "M10.0 — Drift audit kosmologicznych drafts (TGP_v1)"
date: 2026-04-26
cycle: M10
phase: M10.0
status: COMPLETED
verdict: 0 GREEN / 5 YELLOW / 1 RED (corrected v2)
parent: "[[M10_program.md]]"
revisions:
  - "v1 (initial): 2 GREEN / 1 YELLOW / 3 RED — error: ct3/ct7 misclassified as RED"
  - "v2 (corrected): 0 GREEN / 5 YELLOW / 1 RED — sek08a verified V''(1)=-γ slow-roll max"
tags:
  - TGP
  - M10
  - drift-audit
---

# M10.0 — Drift audit kosmologicznych drafts

> **Cel:** weryfikacja czy istniejące drafts kosmologiczne TGP_v1 używają **aktualnych** równań/założeń (post-closure_2026-04-26 + M9 cycle), zanim audytujemy je do closure-grade.
>
> **Metoda:** dla każdego draftu wyciągamy faktyczne `V(Φ)`, `K(Φ)`, source coupling, background metric, parametry, closures referenced — i porównujemy z aktualnymi fundamentami.

---

## Aktualny stan TGP_v1 (reference baseline)

### Single-Φ scalar action (sek08a)

```
S_TGP = ∫ d⁴x √(-g_eff) [ (1/2) K(φ) g_eff^μν ∂_μφ ∂_νφ - V(φ) - (q/Φ_0) φ ρ ]
K(φ) = K_geo φ⁴
V(φ) = (β/3)φ³ - (γ/4)φ⁴
β = γ (vacuum condition, sek08a prop:vacuum-condition)
```

### Φ-EOM (z sek08a + sek08c)

```
∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ_0 - γΦ³/Φ_0² = -qΦ_0 ρ
```

### Hyperbolic metric (M9.1''-P3)

```
ds² = -c² (4-3ψ)/ψ dt² + ψ/(4-3ψ) δ_ij dx^i dx^j
β_PPN = γ_PPN = 1
domain: ψ ∈ (0, 4/3)
```

### Linearyzacja wokół próżni `ψ=1`

**Spatial (M9.3.1):**
```
EOM linearized: ∇²δ - β δ = source
M_eff² = +β > 0    (Yukawa STABLE w sensie spatial range 1/√β)
m_σ² / m_s² = 2    (Path B closure)
```

**Temporal/cosmological (sek08a + de2 V form):**
```
V(ψ) = (β/3)ψ³ - (γ/4)ψ⁴,  β=γ  (sek08a line 56)
V(1)   = β/12          (residual vacuum energy = Λ source)
V'(1)  = 0             (vacuum cond)
V''(1) = -β            (SLOW-ROLL MAXIMUM at ψ=1)
```

**Reconciliation (CRITICAL):** spatial Yukawa stability i temporal slow-roll **KOEGZYSTUJĄ** w sek08a:
- Spatial: `∇²δ - β δ = source` → Yukawa screening, range `1/√β`.
- Temporal FRW: `δ_ddot + 3Hδ_dot - β δ = 0` → slow-roll growing mode (Hubble-damped).
- Differences come from `□` w (-,+,+,+) signature + non-canonical `K(φ)=K_geo·φ⁴` + `√(-g_eff)=c₀ψ`.
- Cosmologically: ψ slow-rolls, providing DE behavior z `w(z) ≥ -1`.

**Errata v1 → v2:**
- v1 incorrectly flagged ct3/ct7 jako RED for "wrong sign V''=-γ tachyonic".
- v2: `V''(1) = -γ` jest **POPRAWNY** w sek08a/de2 (slow-roll max). ct3/ct7 są YELLOW (kanoniczny K=1 vs sek08a K=φ⁴), nie RED.

**Real drift:**
- Drafts używają canonical kinetic `K=1` (de2 line 28-30 explicit); sek08a ma `K(φ)=K_geo·φ⁴` non-canonical.
- Near vacuum ψ=1: K≈K_geo=const, sub-leading. Daleko od vacuum (inflacja, nieliniowy reżim) — może być istotne.

### Closure_2026-04-26

- **T-Λ:** Φ_eq = H_0, ρ_vac = M_Pl² H_0² / 12
- **T-α:** α(ψ) = α_0 (ψ-1)² Θ(ψ-1), α_0 ≈ 4
- **Path B:** σ_ab = ⟨(∂_a δŝ)(∂_b δŝ)⟩^TF, m_σ² = 2 m_s²

### M9.2

- m_field = ∫(∇ε_eq)² d³x (kanoniczny field momentum)

---

## Audit per draft

### Draft 1: [[../desi_dark_energy/de2_tgp_frw_evolution.py]]

**Date:** 2026-04-06

**V(Φ):** `U(ψ) = (β/3)ψ³ - (γ/4)ψ⁴`, normalizowane jako `V(ψ) = 4ψ³ - 3ψ⁴` (V(1)=1, V'(1)=0, V''(1)=-12 w normalizacji `V₀`)

**K(Φ):** Canonical kinetic only: `ρ_ψ = (1/2)ψ̇² + V₀·V(ψ)`. UWAGA: brak K_geo·φ⁴ non-canonical kinetic (uproszczenie liniowe).

**Source coupling:** Brak explicit `qΦρ` — używa Friedmann equation `H² = (8πG/3)(ρ_m + ρ_r + ρ_ψ)` jako dynamicznego sprzężenia.

**Background metric:** Flat FRW (Friedmann), brak hyperbolic metric M9.1''.

**Key params:** `Ω_m = 0.315`, `Ω_r = 9.1×10⁻⁵`, `Ω_DE ≈ 0.685`, `δ ∈ [10⁻⁴, 10⁻²]` (initial displacement od ψ=1), shooting V₀.

**Closures referenced:** Brak explicit (ani T-Λ ani Path B). Ale Φ_eq=H_0 jest implicite przez `Ω_Λ=0.685` matching.

**Predictions:**
- `w(z) ≥ -1` strukturalnie (proven analytically: `w + 1 = ψ̇²/ρ_ψ ≥ 0`)
- `w_0 ≈ -0.9 do -1.0`, `w_a ≈ 0` (near-ΛCDM)
- Falsifiable: DESI `w_0 ≠ -1 ± 0.21` at 3σ

**Drift status:** **GREEN** (zgodny z aktualną akcją sek08a)

**Drift notes:**
- ✅ `V(ψ)` zgodne z β=γ formą.
- ⚠️ Canonical kinetic only (no K_geo·φ⁴) — to OK dla scalar field FRW, ale TGP-strict wymagałoby K(φ)=φ⁴. **Decyzja audit M10.1:** test z OBYDWIEMA formami, sprawdzić czy `w(z)≥-1` zachowane.
- ⚠️ Brak hyperbolic metric — używa standard FRW. Dla cosmology-scale OK (low-curvature limit), ale warto verify.
- 🔧 **M10.1 audit plan:** dodaj K_geo·φ⁴ kinetic + hyperbolic-metric-corrected Friedmann (jeśli zmienia wynik), inaczej GREEN.

---

### Draft 2: [[../nbody/examples/ex261_inflation_tgp.py]]

**Date:** 2026-04-06

**V(Φ):** `P(g) = (β/7)g⁷ - (γ/8)g⁸` z `β=γ` (vacuum cond.); conformal frame: `V_eff(g) = g⁴/8 - g³/7`; hilltop: `V_infl ≈ V₀(1 - (g/μ)³)`

**K(Φ):** `K(g) = g⁴` (non-canonical, rescale do canonical w conformal frame `K_eff = 1` po `g → φ = g³/3`)

**Source coupling:** Brak (inflation epoch — pre-matter)

**Background metric:** FRW during inflation, conformal `g_μν = g²η_μν`

**Key params:**
- `g_0_e = 0.86941` (end-of-inflation field value)
- `Ω_Λ = 0.6847` (current epoch)
- `N = 3` (number of generations)
- `|GL(3,F₂)| = 168`
- `n_s = 1 - 2/N_e ≈ 0.967` (Planck 0.9649)
- `r ~ 0.001-0.01`
- `N_e ~ 60` (e-folds)
- `E_infl ~ M_Pl/168^(1/4) ~ 7×10¹⁷ GeV`

**Closures referenced:** Pre-T-Λ (sprzed closure_2026-04-26), używa GL(3,F₂) heuristic z N=3 jako geometric origin.

**Predictions:**
- `n_s = 1 - 2/N_e` (Starobinsky-like!)
- `r ≪ 0.036` (BICEP/Keck-safe)
- `dn_s/dlnk = -2/N_e² ≈ -5.6×10⁻⁴`
- `T_reh ~ 10¹¹ GeV` (sufficient for leptogenesis)
- No monopoles (`π_2` trivial)
- Possible cosmic strings z Z_2

**Drift status:** **GREEN*** (z gwiazdką — wymaga weryfikacji `g⁷/g⁸` vs `(β/3)φ³ - (γ/4)φ⁴`)

**Drift notes:**
- 🔍 **CRUCIAL:** Standardowa TGP akcja sek08a ma `V = (β/3)φ³ - (γ/4)φ⁴` (potęgi 3-4). ex261 używa `(β/7)g⁷ - (γ/8)g⁸` (potęgi 7-8). Te formy są związane przez `g = φ^(?)` substitution? Trzeba zweryfikować.
- 💡 Hipoteza: jeśli `g = φ^(1/p)` z `p` od kaskady GL(N, F_2), to potęgi w V są skalowane. Dla N=3: `2N+1=7`, `2(N+1)=8` — pasuje strukturalnie.
- 🔧 **M10.2 audit plan:** verify map `g → φ` (substitution + Jacobian) i sprawdzić czy `n_s, r` są niezmiennikami transformacji. Jeśli tak — GREEN. Jeśli różnica — re-derive w canonical sek08a formie.

---

### Draft 3: [[../galaxy_scaling/gs66_frw_propagator.py]]

**Date:** Not visible w header

**V(Φ):** `U'(Φ) = -γ` (linearized form, gradient `-γ` używane jako effective force)

**K(Φ):** Canonical kinetic (linearized regime)

**Source coupling:** Linear Φ couples to matter (used numerically, source schematic)

**Background metric:** Flat FRW; linearization wokół `Φ = Φ_0` (vacuum)

**Key params:**
- `γ > 0` (cosmologically natural)
- `H(a)` from Friedmann
- 3 scenarios: `L_nat = 3 kpc`, `L_H` (cosmological), `H_0` (minimum)

**Closures referenced:** Brak

**Predictions:** **CRITICAL NEGATIVE RESULT (honest):**
- Propagator `G(r) = exp(-μr) / (4πr)` z `μ = √(-γ + 3iH²)`
- 1/r at small r, exponentially screened at large r — **never log(r)**
- Fourier-power proof: `log(r)` requires `G̃(k) ~ 1/k³` jako `k→0`, ale TGP daje `G̃(k) → const lub 1/k²`
- **Bridge (a) "MOND z TGP" FALSIFIED**

**Drift status:** **YELLOW** (eq. forma sub-leading, wynik szczery ale derivation z linearized U' nie z full Φ-EOM)

**Drift notes:**
- ⚠️ `U'(Φ) = -γ` to TYLKO leading-order linearization wokół próżni, nie pełne `βΦ²/Φ_0 - γΦ³/Φ_0²`. Próżnia `ψ=1` daje `U'(1) = β - γ = 0` (vacuum cond), a sub-leading expansion ma `U''(1) δ = -β δ` (tachyonic w force sense → stable Yukawa M_eff²=+β w EOM, sprawdzone w M9.3.1).
- 💡 Wynik (no log-MOND) prawdopodobnie zachowany w pełnej formie, bo bazuje na Fourier-power argument (uniwersalny).
- 🔧 **M10.3 audit plan:** re-derive w canonical Φ-EOM linearized do M_eff²=+β; verify że Yukawa-only conclusion się utrzymuje.

---

### Draft 4: [[../galaxy_scaling/gs41_cmb_compatibility.py]]

**Date:** Not visible

**V(Φ):** `f(R) = R + R_0^γ R^(1-γ) exp(-(R/R_0)^α)` z `α=4/5`, `γ=0.4` (disk default)

**K(Φ):** Not applicable (f(R) theory, **NIE scalar field**)

**Source coupling:** Modified Friedmann z f(R) derivative terms

**Background metric:** FRW; Friedmann modified by f'(R) terms

**Key params:**
- `a_0 = 1.2×10⁻¹⁰ m/s²`
- `R_0 ~ 10⁻⁴⁷ m⁻²`
- `H_0 = 67.4 km/s/Mpc`
- exponential suppression `exp(-(R/R_0)^0.8)`

**Closures referenced:** Brak

**Predictions:**
- Scalaron mass → ∞ at all cosmological epochs (R ≫ R_0)
- CMB primary anisotropies unchanged
- ISW ≈ 0
- Growth rate ≈ GR
- σ_8 enhancement negligible
- BBN safe
- **TGP f(R) cosmologically invisible by design** (chameleon screening)

**Drift status:** **RED** (używa f(R) framework, NIE single-Φ canonical TGP)

**Drift notes:**
- ❌ **MAJOR DRIFT:** ten draft używa f(R) model, który jest **inną teorią** — nie reprezentuje canonical TGP (single scalar Φ z S_TGP).
- 💡 Konkluzja "CMB-safe" prawdopodobnie zachowana w canonical TGP (Yukawa screening exp(-r/r_Y) działa analogicznie do chameleon), ale derivation musi być z scalar Φ-EOM, nie z f(R).
- 🔧 **M10.4 audit plan:** rebuild from scratch w canonical Φ-EOM. Test ISW z δΦ at recombination, σ_8 growth z linear FRW perturbations canonical, BBN compatibility (T~MeV epoch — Φ powinno być w vacuum). Oczekuję CMB-safe konkluzji, ale z poprawnej derivation.

---

### Draft 5: [[../cosmo_tensions/ct3_dark_matter_backreaction.py]]

**Date:** Not visible

**V(Φ):** `V(ψ) = -γ/2 · (ψ-1)²` (**TACHYONIC** — przeciwny znak vs M9.3.1)

**K(Φ):** Canonical kinetic + geometric friction `3ψ̇²/ψ` (TGP-specific)

**Source coupling:** Schematic via gravitational potential coupling

**Background metric:** FRW + perturbations `δψ`

**Key params:**
- `γ_eff = 12 × 3H_0² Ω_Λ ~ (25 km/s/Mpc)²`
- `f_collapsed = 0.7`
- `v_rms ~ 400 km/s`
- tachyonic timescale `~ 4 Gyr`
- amplification `exp(m_tach/H_0) ~ 1.4×10⁶`

**Closures referenced:** Brak

**Predictions:** **NEGATIVE:**
- Gap between needed `B_ψ/H_0² ≈ 0.174` and achieved `~ 10⁻⁹`: **8-9 orders**
- Tachyonic amplification + geometric friction nie wystarczające dla 5σ H_0 shift
- Multiple mechanisms (kinetic K=ψ⁴, volume √ψ, tachyonic, substrate budget) — wszystkie failują

**Drift status:** **YELLOW** (correct V''=-γ; uses canonical K=1 vs sek08a K=φ⁴)

**Drift notes (v2 corrected):**
- ✅ `V''(1) = -γ` MATCHES sek08a/de2 slow-roll maximum picture. **NIE jest drift sign error** (v1 incorrectly flagged this RED).
- ⚠️ Uses canonical kinetic K=1 (sek08a has K=K_geo·φ⁴ non-canonical). Near vacuum sub-leading; może wpływać na backreaction integral.
- 💡 Geometric friction `3ψ̇²/ψ` term — to dodatkowy przyczynek od non-canonical kinetic structure (good!).
- 🔧 **M10.5 audit plan:** verify backreaction `B_ψ/H_0² ~ 10⁻⁹` z full sek08a kinetic K=φ⁴; conclusion "8-9 order gap" prawdopodobnie się utrzyma (slow-roll timescale 1/√β ~ 1/H_0 zbyt wolny dla 5σ H_0 shift).

---

### Draft 6: [[../cosmo_tensions/ct7_soliton_cosmology.py]]

**Date:** Not visible

**V(Φ):** `V(ψ) = -γ/2(ψ-1)²` (**TACHYONIC** — same issue jak ct3)

**K(Φ):** Canonical + `3ψ̇²/ψ` friction

**Source coupling:** Soliton population on cosmological background

**Background metric:** FRW

**Key params:**
- `γ_eff` jak w ct3
- soliton mass variation `δm/m ~ 20-45%`
- RG running γ(k): only ~0.5% change w Λ between CMB i local (need 8% for H_0 tension)
- soliton inter-spacing `d ~ 10³⁸ × Compton wavelength`

**Closures referenced:** Brak

**Predictions:** **COMPREHENSIVE FALSIFICATION (honest):**
- Soliton population effects negligible
- RG running γ daje only ~0.5% Λ change (need 8%)
- Substrate phase transition frozen since `z >> 10¹⁰`
- Two-scale architecture (particle vs cosmological) unresolved
- **Conclusion: TGP scope = galaxy dynamics (y~0.01-1), NIE cosmology (y>>1)**

**Drift status:** **YELLOW** (same as ct3: correct V''=-γ; canonical K=1 vs sek08a K=φ⁴)

**Drift notes (v2 corrected):**
- ✅ V''(1) = -γ MATCHES sek08a/de2 slow-roll max (NIE drift sign error).
- ⚠️ Canonical K=1 vs sek08a K=φ⁴.
- 💡 RG running γ(k) z LPA' (η=0.044) — analytical tool, może być zachowane.
- 💡 Konkluzja "TGP scope = galaxy-scale (y~0.01-1), NIE cosmology (y>>1)" jest STRUKTURALNA — wynika z chameleon-screening Yukawa, niezależne od sub-leading K corrections.
- 🔧 **M10.5 audit plan:** złącz z M10.5 (ct3 + ct7 razem). Verify RG conclusions z poprawną kinetic; główna konkluzja prawdopodobnie się utrzyma.

---

## Summary table (v2 corrected)

| Draft | V(Φ) | K(Φ) | Background | Closures | Drift | M10 plan |
|-------|------|------|------------|----------|-------|----------|
| de2 | (β/3)ψ³-(γ/4)ψ⁴ ✓ | canonical (K=1) | FRW | implicit Ω_Λ | **YELLOW** | M10.1 verify K=φ⁴ correction |
| ex261 | (β/7)g⁷-(γ/8)g⁸ | g⁴ ✓ | FRW conformal | GL(3,F₂) | **YELLOW** | M10.2 verify g↔φ map |
| gs66 | linearized U'=-γ | canonical (K=1) | FRW | none | **YELLOW** | M10.3 full canonical Φ-EOM |
| gs41 | **f(R) chameleon** | n/a (NOT scalar) | FRW + f'(R) | none | **RED** | M10.4 rebuild w scalar Φ |
| ct3 | -(γ/2)(ψ-1)² ✓ | canonical (K=1) | FRW + δψ | none | **YELLOW** | M10.5 verify K=φ⁴ correction |
| ct7 | same -(γ/2)(ψ-1)² ✓ | canonical (K=1) | FRW soliton | none | **YELLOW** | M10.5 verify K=φ⁴ correction |

**Total v2:** 0 GREEN, 5 YELLOW, 1 RED

**Główny pattern drift:**
- 5 z 6 drafts używa kanonicznego kinetic K=1 zamiast sek08a `K(φ)=K_geo·φ⁴`. To uproszczenie OK przy ψ≈1 (sub-leading), ale wymaga weryfikacji.
- 1 draft (gs41) używa innego frameworka (f(R) zamiast single-Φ scalar) — RED structurally.

---

## Decisions for M10 closure cycle (v2 corrected)

1. **M10.1 (de2 audit):** Verify K=φ⁴ non-canonical kinetic correction do `w(z)`; expect bound `w≥-1` robust przy ψ≈1; compare DESI DR1.
2. **M10.2 (ex261 audit):** Verify `g↔φ` map (substitution + Jacobian); konfirmacja że n_s=1-2/N_e jest invariant; geometric N=3 origin.
3. **M10.3 (gs66 rework):** Re-derive FRW propagator z full sek08a Φ-EOM (with K=φ⁴); confirm Yukawa-only conclusion (no log-MOND); upewnić się że Fourier-power argument trzyma.
4. **M10.4 (gs41 rebuild):** Drop f(R) framework. Build CMB safety w scalar Φ-EOM perturbations. Test ISW, σ_8, BBN canonically.
5. **M10.5 (ct3+ct7 audit):** Verify backreaction `B_ψ/H_0²` z K=φ⁴ correction (NIE re-fix sign V''!); honest conclusion że TGP scope = galaxy-scale.

---

## Cross-check vs M9 closures

| M9 closure | Used in M10 sub-cycle |
|------------|------------------------|
| M9.1'' hyperbolic metric | M10.4 (CMB), M10.3 (FRW propagator) |
| M9.2 m_field inertia | M10.5 (backreaction calculation) |
| M9.3.1 stable Yukawa M_eff²=+β | M10.3, M10.5 (correct sign V'') |
| Path B σ_ab m_σ²=2m_s² | M10.4 (tensor mode in CMB perturbations) |
| T-Λ Φ_eq=H_0 | All sub-cycles (boundary condition) |
| T-α α(ψ_NS)≈0.65 | M10.4 (compact-object lensing in CMB?) — defer |

---

## Outstanding questions (post drift audit)

1. **g ↔ φ substitution w ex261:** czy to legalna re-parametrization, czy implicit nowy field? Wymaga careful sympy check.
2. **K_geo·φ⁴ effect on w(z):** czy non-canonical kinetic zmienia bound `w≥-1`? Hipoteza: zachowany strukturalnie (kinetic energy density positive).
3. **Hyperbolic metric in cosmology:** czy `ds² = -c²(4-3ψ)/ψ dt² + ...` redukuje się do FRW dla `ψ` cosmologically slow-varying? Potential M9.1'' extension target.
4. **M9.3 σ_ab tensor in CMB:** czy m_σ² = 2 m_s² produkuje detectable tensor signal w CMB B-mode? Cross-link M9.3 → M10.4.

---

*Drift audit completed 2026-04-26. Ready to proceed to M10.1 setup.*
