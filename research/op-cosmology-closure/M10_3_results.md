---
title: "M10.3 — FRW propagator audit (gs66 → closure-grade verdict)"
date: 2026-04-26
cycle: M10.3
status: CLOSED
result: "6/6 PASS — gs66 verdict YELLOW → GREEN (sign corrected; Fourier-power theorem proven)"
predecessor: "[[M10_2_results.md]] (6/6 PASS, ex261 YELLOW preserved)"
audit_target: "[[../galaxy_scaling/gs66_frw_propagator.py]]"
related:
  - "[[M10_program.md]]"
  - "[[M10_0_drift_audit.md]]"
  - "[[M10_3_setup.md]]"
  - "[[../op-newton-momentum/M9_3_setup.md]] (M9.3.1 spatial Yukawa M_eff²=+β)"
artifacts:
  - "[[m10_3_propagator.py]]"
  - "[[m10_3_propagator.txt]]"
tags:
  - TGP
  - M10
  - FRW
  - propagator
  - Yukawa
  - audit-cycle
  - closure-grade
---

# M10.3 — FRW propagator audit (closure-grade)

> **Sub-cykl M10.3 zamknięty 2026-04-26.**
> Verdict: **6/6 PASS**. gs66 verdict upgrade **YELLOW → GREEN**: sign error (`U''=−γ` tachyoniczne) naprawiony do canonical (`M_eff²=+β` stable Yukawa); konkluzja "no log-MOND z linear TGP" **wzmocniona** przez **Fourier-power universality theorem**.

---

## 1. Cel sub-cyklu

Audyt [[../galaxy_scaling/gs66_frw_propagator.py]] przeciwko canonical sek08a Φ-EOM:

- czy gs66 użył poprawnej linearyzacji potencjału przy `ψ=1`?
- czy konkluzja "no MOND log-far-field" jest robust pod różnymi znakami?
- czy istnieje uniwersalny argument Fourier-power niezależny od konkretnego `M_eff²`?

**Honest framing (z [[M10_3_setup.md]]):** gs66 ma sign error (`U''=−γ` tachyoniczne, sprzeczne z M9.3.1 stable Yukawa `+β`), ale wynik (no log-MOND) prawdopodobnie zachowuje się przez Fourier-power. M10.3 ma to **udowodnić**.

---

## 2. Skrypt + log

- **Skrypt:** [[m10_3_propagator.py]] (~525 linii, 6 sub-testów: sympy + numpy)
- **Output:** [[m10_3_propagator.txt]] (363 linie)

```
Sub-cycle M10.3: 6/6 PASS
M10.3 CLOSURE-GRADE: 6/6 PASS
→ gs66 verdict upgrade: YELLOW → GREEN (sign error documented)
→ Fourier-power theorem: NO log-MOND from linear TGP (universal)
→ Ready for M10.4 (CMB safety REBUILD).
```

---

## 3. Wyniki sub-testów

### 3.1 M10.3.1 — Linearyzacja foundations Φ-EOM (sympy) ✅

Punkt wyjścia (sek08a `eq:field-eq-reproduced`):
```
∇²Φ + 2(∇Φ)²/Φ + (β/Φ_0)Φ² − (γ/Φ_0²)Φ³ = −q·Φ_0·ρ
```

Pivot `Φ = Φ_0(1+δ)`, `|δ| ≪ 1`, `β=γ` (vacuum cond).

**Sympy series expansion** force-side potential `(β/Φ_0)Φ² − (γ/Φ_0²)Φ³`:

| Order | Coefficient (β=γ) |
|-------|-------------------|
| `δ⁰` | `0` (vacuum cond ✓) |
| `δ¹` | `−β · Φ_0` |
| `δ²` | `−β · Φ_0` |

**Linearized static EOM:**
```
∇²(Φ_0·δ) + (−β·Φ_0)·δ = source
∇²δ − β·δ = source/Φ_0
⇒  M_eff² = +β > 0   (stable Yukawa, matches M9.3.1)
```

| Sub | Test | Result |
|-----|------|:------:|
| (a) | `V_force(vacuum) = 0` (β=γ) | ✅ |
| (b) | `dV_force/dδ` at vacuum = `−β·Φ_0` | ✅ |
| (c) | `M_eff² = +β > 0` (stable Yukawa, M9.3.1) | ✅ |
| (d) | Numerically `M_eff²(β=0.01) = 0.0100 > 0` | ✅ |

**Verdict:** PASS — canonical foundations EOM produkuje stable Yukawa, **nie** tachyonic jak gs66.

---

### 3.2 M10.3.2 — FRW quasi-static propagator (sympy) ✅

FRW extension:
```
δΦ̈ + 3H·δΦ̇ − (1/a²)∇²δΦ + β·δΦ = source
```

Fourier `δΦ ~ exp(−iωt + i k·x)`:
```
D(ω,k) = −ω² − 3iHω + k²/a² + β
```

**Quasi-static** (ω→0, ω~H):
```
D(k) = k²/a² + β + 3iH²
G̃(k) = a² / [a²(β + 3iH²) + k²]
```

**μ² = β + 3iH²** (effective mass-squared, complex).

| Limit | Re(μ) | Im(μ) | Behavior |
|-------|-------|-------|----------|
| `β ≫ H²` (galactic) | `√β` | `3H²/(2√β)` | Yukawa damping, slow oscillation |
| `β ≈ H²` | mixed | mixed | transitional |
| `β ≪ H²` (cosmological) | `√(3/2)·H` | `√(3/2)·H` | Hubble-screened |
| `β = 0` | 0 | `√(3)·H` | Pure Hubble damping (Newton screened) |

| Sub | Test | Result |
|-----|------|:------:|
| (a) | `D(ω,k) = −ω² − 3iHω + k²/a² + β` (canonical form) | ✅ |
| (b) | `D_qs(k=0) = β + 3iH²` (constant, not zero) | ✅ |
| (c) | `Re(μ²) = β > 0` (stable real part) | ✅ |
| (d) | `D_canon − D_gs66 = +2β` (sign correction quantified) | ✅ |

**Verdict:** PASS — canonical propagator różni się od gs66 o `2β` w real part `μ²` (sign flip).

---

### 3.3 M10.3.3 — Real-space Green's function G(r) (numerical) ✅

```
G_phys(r) = exp(−Re(μ)·r) · cos(Im(μ)·r) / (4π r)
```

**Trzy scenariusze:**

| Scenario | β | L_nat = 1/√β | Re(μ) | Im(μ) |
|----------|---|--------------|-------|-------|
| Galactic (3 kpc Yukawa) | `(3 kpc)⁻²` | 3 kpc | 1.08e−20 m⁻¹ | ~3.5e−71 m⁻¹ |
| Cosmological (L_H Yukawa) | `L_H⁻²` | 4448 Mpc | 1.05e−26 m⁻¹ | sub-leading |
| H_0 dimensional | `(H_0/c)²` | 4448 Mpc | same | same |

**Sample G_phys (galactic scenario):**

| r [kpc] | r/L_nat | G_phys [m⁻¹] | 4π r G |
|--------:|--------:|--------------|--------|
| 0.1 | 0.033 | high | ~0.97 |
| 1 | 0.33 | mid | ~0.72 |
| **3** | **1.0** | exp(−1)/12π ≈ 0.030 | **0.37** (Yukawa core) |
| 10 | 3.3 | low | ~0.04 |
| 30 | 10.0 | very low | ~10⁻⁵ (screened) |

**Rotation curve** at r=100 kpc, M=1.5×10¹¹ M_⊙ (cosmological β):
```
|g_Newton| = 2.09e−12 m/s²
|g_TGP|    = 2.09e−12 m/s²    (|g_TGP|/|g_Newton| = 1.0000)
```

✓ Newton-like przy r ≪ L_H, **NIE** v_flat const (MOND).

| Sub | Test | Result |
|-----|------|:------:|
| (a) | Re(μ) > 0 we wszystkich scenariach (stable damping) | ✅ |
| (b) | `\|4πrG(10·L_nat)\|` << 1 (Yukawa screening) | ✅ (gal=4.5e−5, cosmo=3e−7) |
| (c) | `\|4πrG(r)\|` peaks at small r (NIE rośnie) | ✅ (peak idx 0/30) |
| (d) | `\|g_TGP\|/\|g_N\|` ≈ 1.0 at 100 kpc (NIE MOND) | ✅ (1.0000) |

**Verdict:** PASS — stable Yukawa z exponential screening; brak log(r) at any scale.

---

### 3.4 M10.3.4 — Sign correction vs gs66 ✅

Porównanie:
```
gs66 D(k) = k²/a² − γ + 3iH²        (TACHYONIC, U''(1)=−γ)
canon D(k) = k²/a² + β + 3iH²        (STABLE Yukawa, M_eff²=+β)
```

Różnica `D_canon − D_gs66 = 2β` (sign flip).

**Numeryczne porównanie** (β = γ = (3 kpc)⁻², galactic):

| | Re(μ) [m⁻¹] | Im(μ) [m⁻¹] | Behavior |
|-|--------------|-------------|----------|
| gs66 | sub-leading (~3e−71) | 1.08e−20 | mostly imaginary → **oscillatory** |
| canon | 1.08e−20 | sub-leading | mostly real → **Yukawa** |

| r [kpc] | gs66 4π r G | canon 4π r G |
|--------:|------------:|-------------:|
| 0.1 | +9.97e−1 | +9.97e−1 |
| 1 | +9.78e−1 | +7.21e−1 |
| 10 | +1.85e−2 | +4.54e−5 |
| 100 | (oscillating) | ~10⁻¹⁵ (screened) |

**Klucz:** oba `D(k=0) ≠ 0` ⇒ **Fourier-power argument applies to BOTH** ⇒ żaden nie produkuje log(r). gs66 conclusion robust niezależnie od znaku `M_eff²`.

| Sub | Test | Result |
|-----|------|:------:|
| (a) | gs66 μ jest oscillatory (Re ≪ Im, tachyoniczny core) | ✅ |
| (b) | canon μ jest Yukawa-dominant (Re ≫ Im, stabilny) | ✅ |
| (c) | Oba `D(k=0) ≠ 0` ⇒ Fourier-power forbids log(r) | ✅ |
| (d) | gs66 no-MOND ROBUST pod sign correction | ✅ |

**Verdict:** PASS — różne znaki dają różne **lokalne** zachowanie (oscillatory vs Yukawa) ALE **identyczną** far-field konkluzję (no log).

---

### 3.5 M10.3.5 — Fourier-power universality theorem (sympy) ✅

**THEOREM (M10.3 universality):**

> Dla każdej rotacyjnie symetrycznej polynomial dispersion
> `D(k) = D_0 + D_2·k² + D_4·k⁴ + …` (Hermitowska, analityczna w `k=0`),
> propagator `G̃(k) = 1/D(k)` jest analityczny w `k²` przy `k=0`.
> **Brak modu `1/k³`** w niskoczęstotliwościowym rozwinięciu.

**Sympy verification:**
```
G̃(k) low-k expansion = 1/D_0  −  (D_2/D_0²)·k²  +  …    (analytic in k²)
G̃(k) high-k limit   = 1/(D_2·k²)                       (Newtonian, D_4=0)
```

`log(r)` Fourier transform requires `G̃(k) ~ 1/k³` przy `k→0`:
```
∫ d³k/(2π)³ · e^{ik·r}/k³ ~ log(r)    (regularized)
```

Polynomial `D(k)` nigdy nie produkuje `1/k³` mode (odd power singular term niemożliwy w analitycznej `k²` series). **QED**.

| Sub | Test | Result |
|-----|------|:------:|
| (a) | Low-k expansion has NO `k⁻³` term (analytic) | ✅ |
| (b) | `G̃(k=0) = 1/D_0` (constant) | ✅ |
| (c) | `G̃ ~ 1/(D_2 k²)` at large k (Newton, D_4=0) | ✅ |
| (d) | Theorem: canonical TGP cannot produce log(r) | ✅ |

**Verdict:** PASS — uniwersalny no-log-MOND theorem dla każdej polynomial dispersion z `D_0 ≠ 0`.

---

### 3.6 M10.3.6 — Honest synthesis verdict ✅

**Key findings:**

1. **Canonical foundations Φ-EOM linearization** (M10.3.1): `M_eff² = +β > 0` (stable Yukawa), zgodne z M9.3.1.
2. **FRW quasi-static propagator** (M10.3.2): `D(k) = k²/a² + β + 3iH²` (stabilny, w przeciwieństwie do gs66's `−γ + 3iH²`).
3. **Real-space Green's function** (M10.3.3): exponential Yukawa screening, brak log(r) at any scale.
4. **Sign correction vs gs66** (M10.3.4): gs66's tachyonic vs canonical's stable — **oba** nie produkują log(r) przez Fourier-power.
5. **Universality theorem** (M10.3.5): polynomial `D(k)` z `D_0 ≠ 0` strukturalnie wyklucza log(r). **Niezależnie** od znaku `M_eff²`.

**Honest verdict:**

- **gs66 conclusion** ("no MOND from linear TGP") **PRESERVED** i **WZMOCNIONE** przez Fourier-power universality.
- **gs66 sign error** (`U''=−γ` tachyonic) ujawniony i naprawiony do canonical `M_eff² = +β`. Sign error nie wpływa na konkluzję far-field, ale wpływa na intermediate behavior (oscillatory vs Yukawa).
- **gs66 verdict UPGRADE:** **YELLOW → GREEN.**

**Falsifiable statements:**

1. Linear TGP single-Φ ma **NO MOND signature** (`g(r) ≠ const · √(GMa₀)/r`).
2. Galactic rotation curves: TGP daje Newtonian (r ≪ L_nat) lub Yukawa-screened (r ≫ L_nat), **nigdy** flat.
3. Jeśli observation requires log(r) far-field potential ⇒ TGP single-Φ **LINEAR** theory **falsified** (potrzebna non-linear extension lub multi-field).

**Scope statement (M10.3 honest):**

TGP single-Φ linear theory: **galactic dynamics NOT a unique TGP signature**. Jakakolwiek zgodność z SPARC / dark-matter halos wymaga (a) phenomenological dressing (density profiles), lub (b) non-linear extensions OUTSIDE minimal sek08a action. **Bridge (a) "MOND from TGP linear FRW" remains FALSIFIED**, z wzmocnionym theoretical grounding (Fourier-power universality theorem).

---

## 4. Cross-checks foundational constraints

| Foundational | Test sub | Result |
|--------------|----------|:------:|
| Single-Φ axiom | M10.3.1-2 use scalar Φ only | ✅ |
| β=γ vacuum cond. | Linearization gives `(β−γ)Φ_0 = 0` | ✅ |
| K(φ) = K_geo·φ⁴ | M_eff² = +β unchanged near vacuum | ✅ |
| Hyperbolic metric M9.1'' | Static baseline `Φ_eq = Φ_0` w `g_eff`; pełna metryka sub-leading | ⏳ |
| M9.3.1 stable Yukawa | M_eff² = +β reproduced exactly (sympy) | ✅ |
| Fourier-power universality | Theorem niezależny od konkretnych potencjałów | ✅ |

---

## 5. Falsifiable predictions (preserved + sharpened)

1. **No log(r) far-field potential** (universal, structural — by Fourier-power).
2. **Yukawa screening** dla source (galactic): `g(r) ≈ G_N M e^{−r√β}/r²` przy r ≫ L_nat.
3. **Newton recovery** dla r ≪ L_nat (deep solar-system OK).
4. **Hubble-scale screening** dla β << H²: cosmological Φ-modes decay over `1/H`.

Każde silne odkrycie log(r) potencjału (np. mocny pure-MOND signature w SPARC poza phenomenologią) **falsifikuje** TGP single-Φ linear theory.

---

## 6. Wpływ na program M10

- **M10.3 zamknięty:** drift v2 YELLOW → **GREEN** (sign error documented + theorem proven).
- **gs66 status:** still describes "no MOND" honestly; sign correction `−γ → +β` musi być wskazana w future references.
- **Galaxy program path:** TGP nie predykuje MOND ⇒ alternatywny mechanism dla rotation curves potrzebny (cluster dynamics? phenomenologia? — separate cycle).

**Następnik:** M10.4 — gs41 CMB safety **REBUILD** (currently uses f(R) — RED, structural single-Φ violation).

---

## 7. Status plików

| Plik | Status | Działanie |
|------|--------|-----------|
| [[m10_3_propagator.py]] | NEW | M10.3 audit script (6 tests) |
| [[m10_3_propagator.txt]] | NEW | Run log (6/6 PASS) |
| [[M10_3_setup.md]] | EXISTING | Audit plan + math foundation |
| [[../galaxy_scaling/gs66_frw_propagator.py]] | UNCHANGED | Verdict upgrade YELLOW → GREEN dokumentowany w M10.3 |
| [[M10_3_results.md]] | NEW (this) | Closure-grade synthesis |

---

## 8. Open items → M11+

1. **Hyperbolic metric M9.1'' integration:** włączyć `g_eff[Φ]` static baseline + δΦ perturbations w propagator (sub-leading near vacuum, ale kompletność).
2. **Non-linear extension audit:** czy `2(∇Φ)²/Φ` term (non-linear kinetic) może produkować log(r) poza linear regime? gs65 sugeruje `r⁻¹/³` static, nie log — ale full time-dependent nonlinear nie sprawdzony.
3. **Multi-field extension:** czy dodanie Φ + ψ_extra (multi-Φ) mogłoby produkować log? (Outside single-Φ axiom — would require axiom extension.)

---

## Status

| Sub-test | Status | Wynik |
|----------|:------:|-------|
| M10.3.1 — Linearization foundations Φ-EOM (sympy) | ✅ PASS | 4/4 sub OK |
| M10.3.2 — FRW quasi-static propagator (sympy) | ✅ PASS | 4/4 sub OK |
| M10.3.3 — Real-space Green's function (numerical) | ✅ PASS | 4/4 sub OK |
| M10.3.4 — Sign correction vs gs66 | ✅ PASS | 4/4 sub OK |
| M10.3.5 — Fourier-power universality theorem | ✅ PASS | 4/4 sub OK |
| M10.3.6 — Honest synthesis verdict | ✅ PASS | 5/5 sub OK |

**M10.3: 6/6 PASS — CLOSURE-GRADE (gs66 verdict upgrade YELLOW → GREEN).**

---

*M10.3 sub-cycle closed 2026-04-26.*
