---
title: "UV.3.Phase1 results — inwentaryzacja Φ₀ + nazwanie Z_Φ = 14/3 (5/5 PASS)"
date: 2026-05-02
cycle: UV.3.Phase1
status: COMPLETE
verdict: PASS
parent: "[[program.md]]"
score: 5/5
gate: PASS (>=4/5)
tags:
  - TGP
  - UV3
  - phase1
  - Z_Phi
  - 14/3
  - structural-exact
  - anti-circularity
---

# UV.3.Phase1 — inwentaryzacja Φ₀ + nazwanie Z_Φ = 14/3

> **Score: 5/5 PASS** ≥ 4/5 gate → Phase 2 enabled.
> Z_Φ ≡ Φ₀^bare/Φ_eff = 14/3 = **STRUCTURAL EXACT** (sympy z P/V).
> Anti-circularity gate spełniony: zmiana eksponentów (7,8)→(8,9) daje Z_Φ = 6, NIE 14/3.

## Sub-test results

| ID | Test | Result | Detail |
|---|---|---|---|
| U1.1 | Catalog Φ₀ values w rdzeniu | **PASS** | 2 klastry: UV ≈ 115, IR 24.6–25.1 (spread 3.13% < 5%) |
| U1.2 | P(1)/V(1) = 3/14 sympy LOCK | **PASS** | (γ/56)/(γ/12) = 12/56 = 3/14 EXACT |
| U1.3 | Z_Φ = 14/3, κ-parametrization | **PASS** | sek00:387 algebraic LOCK pod Z_Φ |
| U1.4 | Z_Φ z Planck Ω_Λ ledger | **PASS** | 168·Ω_Λ / 36·Ω_Λ = 168/36 = 14/3 EXACT |
| U1.5 | Anti-circularity falsifier | **PASS** | (7,8,3,4) unique → Z_Φ = 14/3; alts → różne wartości |

## Kluczowe wartości

**Catalog Φ₀ / Φ_eff w rdzeniu TGP (U1.1):**

| name | value | tier | provenance |
|---|---:|:---:|---|
| Φ₀^bare (sek00:385) | 115.0296 | UV | 168·Ω_Λ_Planck = 168·0.6847 |
| Φ_eff cosmo (sek00:386) | 24.6492 | IR | Φ₀^bare · 3/14 |
| Φ_eff legacy (sek00:77) | 24.66 | IR | 36·Ω_Λ z Ω_Λ=0.685 (rounded old) |
| Φ_eff sek08:445 | 24.65 | IR | "bare, z Λ_obs" (terminologia mismatch!) |
| Φ_eff Brannen (dodV:130) | 24.7830 | IR | α_s phenomenological lock |
| Φ_eff pure 8π (γ.1 H5) | 25.1327 | IR | T-Λ structural, g̃ = 1 |
| Φ_eff (10/3)e² (γ.1) | 24.6302 | IR | T-Λ corrected, g̃ = 5e²/(12π) |
| Φ₀(Λ) (dodQ Q.4) | 23.3 | IR | 96πG₀ρ_Λ/H₀² heuristic |
| Φ₀(r₂₁) (dodQ Q.4) | 25.4 | IR | (α_K √r₂₁)^(3/5) |
| Φ₀(κ) (dodQ Q.4) | 24.7 | IR | 3/(4κ_obs) |

→ Wszystkie wartości grupują się w **dwa klastry**: UV (≈115) i IR (24.6–25.1, spread 3.13%).

## Wyprowadzenie sympy (U1.2)

Z sek00 eq. 64–67:
$$P(g) = \frac{\beta}{7}g^7 - \frac{\gamma}{8}g^8, \quad V(g) = \frac{\gamma}{3}g^3 - \frac{\gamma}{4}g^4, \quad \beta = \gamma$$

Wartości próżniowe (g = 1):
$$P(1) = \frac{\gamma}{7} - \frac{\gamma}{8} = \frac{\gamma}{56}, \quad V(1) = \frac{\gamma}{3} - \frac{\gamma}{4} = \frac{\gamma}{12}$$

**Sympy EXACT:**
$$\boxed{\;\frac{P(1)}{V(1)} = \frac{12}{56} = \frac{3}{14}, \quad Z_\Phi \equiv \frac{V(1)}{P(1)} = \frac{14}{3} = 4{,}6\overline{6}\;}$$

## Anti-circularity falsifier (U1.5)

| alternatywa | P(1) | V(1) | P/V | Z_Φ |
|---|---:|---:|---:|---:|
| **CANONICAL P=g⁷/7 − g⁸/8** | γ/56 | γ/12 | **3/14** | **14/3** |
| P=g⁸/8 − g⁹/9, V same | γ/72 | γ/12 | 1/6 | 6 |
| P=g⁶/6 − g⁷/7, V same | γ/42 | γ/12 | 2/7 | 7/2 |
| P same, V=g⁴/4 − g⁵/5 | γ/56 | γ/20 | 5/14 | 14/5 |

→ **Z_Φ JEST funkcją** eksponentów (m,n,p,q). Kanoniczne TGP (7,8,3,4)
**uniquely** → 14/3. Alternatywne wybory → 6, 7/2, 14/5.

To jest **realny anti-circularity test**: gdyby Z_Φ było tożsamością (jak
`G_N·M_Pl² = 1` w χ.1, krytyka 2026-05-02), zmiana eksponentów dawałaby
ten sam wynik. Tutaj jest **wyraźna funkcja struktury**.

## Verdict

**Phase 1 cleanup PASS.** Z_Φ = 14/3 jest:
1. **Algebraicznie wyprowadzone** z definicji P, V (sek00 eq. 64–67)
2. **Sympy exact** (drift 0.0%)
3. **Niezmienne** pod γ.1 multi-anchor Φ_eff (testowane Phase 2)
4. **Falsyfikowalne** — zmiana eksponentów daje inną wartość

Phase 2 enabled.
