---
title: "M11.S PoC — feasibility check (single soliton Φ_sol)"
date: 2026-04-26
cycle: M11
sub-cycle: M11.0+ (PoC, pre-execution)
status: PoC PASS (4/4)
predecessor: "[[M11_branch_strategy.md]]"
related:
  - "[[m11_S_soliton_poc.py]] (script)"
  - "[[m11_S_soliton_poc.txt]] (output)"
  - "[[../op-newton-momentum/M9_3_results.md]] (M9.3.1 source)"
tags:
  - TGP
  - M11
  - M11.S
  - PoC
  - soliton
---

# M11.S PoC — feasibility summary

> **Cel:** zanim formalnie commitujemy do M11.S audit (Branch I, soliton-based quantization), sprawdzić numerycznie czy klasyczny `Φ_sol(r)` istnieje, czy linearizacja zgadza się z M9.3.1 Yukawa, i czy s-wave spectrum jest stabilny.
> **Wynik:** ✅ **4/4 PASS** — feasibility CONFIRMED, **READY for full M11.S audit**.

---

## Wyniki

| Test | Cel | Wynik |
|---|---|---|
| **PoC.1** | Linear Yukawa benchmark — verify `(∇²-μ²)δφ = -source/K` z M9.3.1 sign | ✅ PASS — Yukawa decay rate exactly `e^(-μΔr)` |
| **PoC.2** | Nonlinear ODE `Φ_sol(r)` z regularized source | ✅ PASS — Φ_sol regular at origin, asymptotic Φ_0 |
| **PoC.3** | Linear vs nonlinear comparison w tail (r > 4 λ_C) | ✅ PASS — 4.21% max diff |
| **PoC.4** | Linearization s-wave spectrum (eigenvalues ω²) | ✅ PASS — 0 negative modes, lowest ω² = 1.069 |

### Detale numeryczne

Units: β = K_geo = Φ_0 = q/Φ_0 = 1; M_source = 1.0; a_source = 0.15 (Gaussian regularization).
Yukawa range λ_C = √(K_geo/β) = 1.0.

| r [λ_C] | Linear δφ_Yuk | Nonlinear Φ_sol−1 | rel diff |
|---:|---:|---:|---:|
| 0.001 | (singular) | +2.762 | (core, dominated by source) |
| 1.0 | +0.0293 | +0.0270 | 7.7% |
| 2.0 | +5.4×10⁻³ | +5.1×10⁻³ | 4.7% |
| 5.0 | +1.07×10⁻⁴ | +1.02×10⁻⁴ | 4.2% |
| 10.0 | +3.6×10⁻⁷ | +3.6×10⁻⁷ | <1% |

### Linearization spectrum (s-wave, around Φ_sol)

```
Lowest 5 ω² values:  1.069, 1.273, 1.613, 2.087, 2.695
Asymptotic ω² → β + k² = 1 + k²        (free Yukawa continuum, k = radial momentum)
0 negative modes (no instabilities)
0 zero modes in s-wave (zero modes are L=1 translational, NOT in s-wave sector)
```

Mass gap ω² ≈ β confirms M9.3.1: stable Yukawa regime around solitonowy background.

---

## Wnioski + finding

### Główne (potwierdza M11.S feasibility)

1. **H1 ✅ existence:** Klasyczny non-trivial Φ_sol(r) wynik numerycznego shootingu istnieje i jest jednoznaczny (przy danym M_source).

2. **H2 ✅ asymptotic + regular:** Φ_sol(r→r_max) → Φ_0=1 dokładnie; Φ_sol(r→0) finite (no naked singularity), zgodne z M9.3.1.

3. **H3 ✅ s-wave stability:** Spectrum fluktuacji wokół Φ_sol ma wszystkie ω² > 0 w sektorze s-fali. Mass gap ω₀² = 1.07 ≈ β, jak przewiduje M9.3.1 stable Yukawa.

4. **Linearization works in tail:** Liniowe Yukawa coincydujhe z pełnym nonlinear Φ_sol w 4% w tail r > 4·λ_C — dobra weryfikacja consistency.

### Finding NOWY (do dalszego audytu w M11.S)

**Φ_sol(r→0) = 3.76 — penetracja poza domenę psi ∈ (0, 4/3) z M9.1''.**

W naszych dimensionless units (q·M = 1) klasyczne pole rośnie powyżej 4/3 blisko źródła, czyli wychodzi poza zakres gdzie M9.1'' hyperbolic metric ma valid Lorentzian signature (M9.3.1 mówi: "psi ∈ (0, 4/3) — sygnatura Lorentzowska zachowana").

**Implikacje dla full M11.S:**
- Realistyczne TGP sources musi mieć słabe coupling: `q·M / (K_geo · Φ_0³) ≲ 1` żeby Φ_sol pozostało w (0, 4/3).
- Alternatywnie: regularized source musi mieć szerokość a_source ≳ λ_C/2 — wówczas response jest "rozmiyty" i nie penetruje wysokich Φ.
- Może to być **kanaliacja konstrukcyjna**: TGP nie dopuszcza arbitralnie skoncentrowanych klasycznych źródeł — fizyczne sources są rozmiyte z natury (matter solitony mają intrinsic core scale).

→ M11.S full audit musi sparametryzować przejście (Φ_sol exits valid domain) jako warunek dla "physical" sources.

### Implikacje dla M11.I + M11.G

- **M11.I:** ansatz Φ_2sol(r) = Φ_sol(r-r_1) + Φ_sol(r-r_2) - Φ_0 jest ważny TYLKO jeśli sources są wystarczająco oddzielone i słabe że Φ wszędzie w (0, 4/3). Trzeba sparametryzować `r_12_min` poniżej którego ansatz fails.
- **M11.G:** decompozycja Φ = Φ_cl[{r_i}] + δΦ_rad zakłada że Φ_cl jest nigdzie nie singular — PoC potwierdza to dla single source, multi-source jest do sprawdzenia w M11.I.

---

## Decyzja (post-PoC)

**M11.S feasibility = CONFIRMED.** Można commitować do full M11.S audit.

Rekomendacja kolejności:
1. **Update M11_program.md** (sub-cykle 0/S/I/G/1/2/3/4/R)
2. **Launch M11.S full audit:**
   - Pełna 1-loop analiza (zeta-regularization or dim reg)
   - δM_phys finite verification (H4)
   - Pełny spectrum (l=0,1,2 modes; zero modes from translations)
   - Domain-of-validity sweep (jaki q·M/Φ_0 zachowuje Φ ∈ (0, 4/3))
3. **M11.1 (Branch II) parallel:** m2b 1-loop audit niezależny od M11.S

---

## Files manifest

```
M11_branch_strategy.md         — strategy (this PoC sits BETWEEN strategy and full audit)
m11_S_soliton_poc.py           — PoC script (≈410 lines)
m11_S_soliton_poc.txt          — execution output
M11_S_PoC_summary.md           — this file
```

---

*M11.S PoC closed 2026-04-26. Feasibility: CONFIRMED 4/4 PASS. Awaiting M11_program.md restructure + M11.S full audit launch approval.*
