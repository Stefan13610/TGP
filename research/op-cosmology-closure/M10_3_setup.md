---
title: "M10.3 — FRW propagator audit (gs66 → canonical Φ-EOM)"
date: 2026-04-26
cycle: M10.3
status: SETUP
predecessor: "[[M10_2_results.md]] (6/6 PASS, ex261 YELLOW preserved)"
audit_target: "[[../galaxy_scaling/gs66_frw_propagator.py]]"
related:
  - "[[M10_program.md]]"
  - "[[M10_0_drift_audit.md]]"
  - "[[../op-newton-momentum/M9_3_setup.md]] (M9.3.1 spatial Yukawa)"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] (Φ-EOM)"
tags:
  - TGP
  - M10
  - FRW
  - propagator
  - audit-cycle
  - galaxy-scaling
---

# M10.3 — FRW propagator audit (canonical sek08a)

> **Cel:** zaudytować [[../galaxy_scaling/gs66_frw_propagator.py]] przeciwko **canonical sek08a Φ-EOM** (z foundatons §3.5) i M9.3.1 (stabilny Yukawa M_eff²=+β). Sprawdzić czy konkluzja **"no log-MOND z linear TGP"** jest robust przy poprawnych znakach.

---

## 1. Audit target — gs66 status v2

Z [[M10_0_drift_audit.md]]:

| Aspekt | gs66 | sek08a / M9.3.1 |
|--------|------|-----------------|
| Linearized force | `U''(1) = −γ` (tachyonic) | `M_eff² = +β` (stable Yukawa) |
| Kinetic | canonical `K=1` | `K(φ) = K_geo·φ⁴` |
| Source | schematic | `q·Φ_0·ρ` (sek08a coupling) |
| Background | flat FRW | hyperbolic M9.1'' |
| Conclusion | "no MOND" via Fourier-power | (claimed robust) |

**Drift v2:** YELLOW — formalna usterka znaków/równań, ale Fourier-power proof (`log(r)` requires `G̃(k) ~ 1/k³` at `k→0`) jest **uniwersalny** i prawdopodobnie zachowuje się przy korekcjach.

**Cel M10.3:** **udowodnić** że konkluzja gs66 jest robust przez **rebuild z canonical EOM** i nie zależy od konkretnych znaków.

---

## 2. Strukturalny punkt wyjścia — Foundations Φ-EOM

Z [[../op-newton-momentum/M9_3_setup.md]] §3.1 (sek08a `eq:field-eq-reproduced`):

```
∇²Φ + 2(∇Φ)²/Φ + (β/Φ_0)Φ² − (γ/Φ_0²)Φ³ = −q·Φ_0·ρ
```

**Wakuum:** `Φ = Φ_0 ⇒ ψ = Φ/Φ_0 = 1`, `(β−γ)Φ_0 = 0` (β=γ vacuum cond).

### 2.1 Linearyzacja statyczna (spatial M9.3.1 result)

Pivot `Φ = Φ_0(1+δ)`, `|δ| ≪ 1`:

- **Kinetic gradient term** `2(∇Φ)²/Φ`: dla statycznego baseline `Φ_eq = Φ_0` mamy `∇Φ_0 = 0`, więc bilinear `2·(2∇Φ_0·∇δΦ)/Φ_0 − 2(∇Φ_0)²δΦ/Φ_0² = 0` → znika w linear order
- **Potential force** `(β/Φ_0)Φ² − (γ/Φ_0²)Φ³` linearized:
  ```
  (β/Φ_0)(Φ_0² + 2Φ_0·δΦ) − (γ/Φ_0²)(Φ_0³ + 3Φ_0²·δΦ) + O(δΦ²)
  = (β−γ)Φ_0 + (2β−3γ)·δΦ
  = 0 + (2β−3γ)·δΦ                                     [β=γ]
  = −β·δΦ
  ```

**Statyczny linearized EOM:**
```
∇²δΦ − β·δΦ = −q·Φ_0·δρ
```

✓ Yukawa z **stabilnym** `M_eff² = +β`, dokładnie jak w M9.3.1 (i sympy `m9_3_gw.py:64-89`).

### 2.2 FRW extension

Dla niejednorodnego pola w FRW background `a(t)`:
```
δΦ̈ + 3H·δΦ̇ − (1/a²)∇²δΦ + β·δΦ = q·Φ_0·δρ      [sign convention M9.3.1]
```

Fourier `δΦ ~ e^{−iωt + i k·x}`:
```
D(ω,k) = −ω² − 3iHω + k²/a² + β
```

**Quasi-static** (`ω → 0`, `ω ~ H` from source slow evolution):
```
D(k) ≈ k²/a² + β + 3iH²
G̃(k) = 1/(k²/a² + β + 3iH²)
```

⚠ **Sign correction vs gs66:** gs66 ma `−γ + 3iH²` (tachyoniczne), my mamy `+β + 3iH²` (stabilne). Wynik Re(μ) i Im(μ) inny, ale Fourier-power argument (k→0) zachowany.

### 2.3 Real-space Green's function

```
G(r) = exp(−μ_eff·r) / (4π r),    μ_eff² = β + 3iH²
```

W limit `β ≫ H²` (Yukawa-dominant):
- `Re(μ_eff) ≈ √β` (damping length `1/√β`)
- `Im(μ_eff) ≈ 3H²/(2√β)` (slow oscillation)

**Stable Yukawa, exponential screening.** Brak log(r) at any scale.

### 2.4 Non-canonical kinetic K=K_geo·ψ⁴

Pełne sek08a ma `K(ψ) = K_geo·ψ⁴`. Linearyzacja przy ψ=1:
- `K(1) = K_geo`, `K'(1) = 4K_geo`
- Gradient term: `2(∇Φ)²/Φ` w foundations EOM **już zawiera** efekt non-canonical kinetic (bo Foundations EOM jest variational z full action sek08a)
- Linearized EOM (statyczny baseline): same form, M_eff² = +β unchanged near vacuum

Sub-leading correction przy `δ ≠ 0` jest `O(δ²)` (jak w M10.1 dla DE).

---

## 3. Fourier-power argument — universal no-log-MOND

Klucz **niezależny od znaków**: log(r) wymaga `G̃(k) ~ 1/k³` przy k→0 (z odwrotnej Fouriera w 3D):

```
∫ d³k/(2π)³ · e^{ikr} / k³ ~ log(r)   (regularized)
```

Dla LINEAR theory z polynomial dispersion `D(k) = D_0 + D_2·k² + O(k⁴)`:
- `k → 0`: `G̃(k) → 1/D_0 = const` (not 1/k³)
- `k → ∞`: `G̃(k) → 1/(D_2·k²)` (Newton 1/k², not 1/k³)

**Niezależnie od znaku M_eff²:** jeśli `D_0 ≠ 0`, far-field jest exp-suppressed (Yukawa) lub oscillating (Helmholtz), nigdy log(r). Jeśli `D_0 = 0`, far-field jest 1/r (Newton).

**Theorem (M10.3):** Dla każdej polynomial dispersion `D(k)` z full sek08a **canonical** linearization wokół vacuum: `G̃(k) → const` przy `k→0`. Bridge (a) MOND **strukturalnie wykluczony** przez Fourier-power, **niezależnie od konkretnych znaków potencjału**.

---

## 4. Sub-tests (planowane)

| ID | Cel | Metoda | PASS criterion |
|----|-----|--------|----------------|
| **M10.3.1** | Linearyzacja statyczna foundations Φ-EOM (sympy) | sympy series expansion `(β/Φ_0)Φ²−(γ/Φ_0²)Φ³` przy ψ=1, β=γ | linear coef = `−β·δΦ`, M_eff² = +β |
| **M10.3.2** | FRW quasi-static propagator (sympy) | sympy `D(ω,k)` derivation, quasi-static limit | `D(k) = k²/a² + β + 3iH²`, sign POSITIVE |
| **M10.3.3** | Real-space Green's function G(r) (numerical) | numpy: G(r) = exp(−μr)/(4πr) na grid, dla 3 scenariuszy β | stable Yukawa, brak log(r) na żadnej skali |
| **M10.3.4** | Sign correction vs gs66 (sympy) | porównanie gs66 `−γ` vs M9.3.1 `+β`; pokazać że oba dają NO log-MOND | gs66's tachyonic NIE produkuje log(r) ALE jest formalnie błędny; canonical +β robust |
| **M10.3.5** | Fourier-power universality (sympy) | symbolic proof: dla każdej polynomial D(k), `G̃(k→0) = 1/D(0) ≠ 1/k³` | uniwersalny no-log-MOND theorem |
| **M10.3.6** | Honest synthesis verdict | dokumentacja: gs66 conclusion preserved, znaki naprawione, robust theorem | PASS jeśli (a)-(e) wszystkie OK |

---

## 5. Foundational constraints check

| Constraint | Test sub | Expected |
|------------|----------|:--------:|
| Single-Φ axiom | M10.3.1-2 use scalar Φ only | ✅ |
| β=γ vacuum cond. | Linearization gives `2β−3γ = −β` correctly | ✅ |
| K(φ) = K_geo·φ⁴ | M_eff² = +β unchanged near vacuum (M10.1 pattern) | ✅ |
| Hyperbolic metric M9.1'' | Static baseline `Φ_eq = Φ_0` w `g_eff` z M9.1''; modyfikacja `2(∇Φ)²/Φ` term | ⏳ (M10.3 sub-leading) |
| M9.3.1 stable Yukawa | M_eff² = +β reproduced sympy | ✅ |
| Fourier-power universality | Theorem niezależny od konkretnych potencjałów | ✅ |

---

## 6. Expected verdict

**Hipoteza pre-test:** M10.3 zwróci **6/6 PASS** z verdictem:
- **gs66 conclusion** ("no MOND from linear TGP") **PRESERVED**
- **gs66 sign error** ujawniony i naprawiony: `M_eff² = +β > 0` (stable Yukawa), nie `−γ < 0` (tachyonic)
- **gs66 verdict upgrade:** YELLOW → GREEN (sign error documented, Fourier-power theorem robust)
- **Falsifiable prediction preserved:** TGP nie produkuje log(r) far-field z LINEAR theory; jakikolwiek jednoznaczny log-MOND signature falsifikuje TGP single-Φ

---

## 7. Files

| Plik | Status | Cel |
|------|--------|-----|
| [[M10_3_setup.md]] | NEW (this) | Plan + math foundation |
| [[m10_3_propagator.py]] | TO CREATE | Audit script (6 tests, sympy + numpy) |
| [[m10_3_propagator.txt]] | TO CREATE | Run log |
| [[M10_3_results.md]] | TO CREATE | Closure-grade synthesis |
| [[../galaxy_scaling/gs66_frw_propagator.py]] | UNCHANGED | Original draft, will retain (verdict update only) |

---

## 8. Następne (M10.4)

- **M10.4:** REBUILD gs41 CMB safety w canonical scalar Φ (currently uses f(R) — RED)
- **M10.5:** ct3+ct7 H₀/S₈ tensions audit (verify K=φ⁴ correction)
- **M10.R:** synthesis cyklu M10

---

*M10.3 setup ready 2026-04-26. Proceed with script.*
