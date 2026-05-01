---
title: "ψ.1.v2.Phase5 setup — eikonal + dispersion + corrected Sagnac (correction phase)"
date: 2026-05-01
cycle: ψ.1.v2.Phase5
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - psi1
  - phase5
  - correction
  - eikonal
  - sagnac
  - setup
---

# ψ.1.v2.Phase 5 setup — eikonal + dispersion + corrected Sagnac (5 sub-tests)

## Reason

Phase 4 strukturalnie ustaliło tensorowy operator L₅'_a oraz efektywną metrykę
optyczną $g^{\mu\nu}_{eff} = \eta^{\mu\nu} + (\xi/\Lambda^2) n^\mu n^\nu$ z
$\xi = +2|\beta_g| > 0$ (Adams positivity). Phase 5 propaguje tę strukturę
przez eikonal limit do realnych obserwabli — anizotropowy $c_{eff}(\theta)$,
re-derivacja Sagnac fazowy z **kierunkowym** chopperem, Yukawa-tensor source
w lab, oraz cross-falsyfikacja przez alt-tensor candidates.

Cel: dostać realistyczny SNR Sagnac (poprawiony, z anizotropią) — oczekiwana
redukcja względem v1 ze SNR ~3×10⁴ do bardziej skromnego (~10²-10³)
po uwzględnieniu, że tylko składowa $\vec{k} \parallel \vec{\nabla}\ln X$
modyfikuje $c$.

## Sub-tests

### T5.1: Eikonal dispersion relation k_μ k_ν g^{μν}_eff = 0 sympy LOCK

Eikonal limit Lagrangian'u L₅'_a:
- $A_\mu = a_\mu \exp(i\phi/\hbar)$, $\hbar \to 0$
- Principal symbol of EOM: kontrakcja $g^{\mu\nu}_{eff} k_\mu k_\nu$
- Verify: dla $g^{\mu\nu}_{eff} = \eta^{\mu\nu} + (\xi/\Lambda^2) n^\mu n^\nu$
  dispersion = $k^2 + (\xi/\Lambda^2)(n\cdot k)^2 = 0$ (sympy)

PASS: sympy LOCK na dispersion relation $k^2 + (\xi/\Lambda^2)(n\cdot k)^2 = 0$.

### T5.2: Anisotropic c_local(θ) sympy LOCK formula

Z dispersion $\omega^2 = |\vec{k}|^2 - (\xi/\Lambda^2)(\omega - \vec{k}\cdot\vec{n})^2$
(dla statycznego $n^\mu = (0, \vec{n})$, $\vec{n} = \vec\nabla \ln X$):

- Solve dla $\omega^2/|\vec{k}|^2 \equiv c_{eff}^2$ z dokładnością do
  pierwszego rzędu w $\xi |\vec{n}|^2/\Lambda^2 \ll 1$
- Verify: $c_{eff}^2(\theta) = 1 - (\xi/\Lambda^2) |\vec{n}|^2 \cos^2\theta$
  gdzie $\theta = \angle(\vec{k}, \vec{n})$

PASS: sympy LOCK na formuł $c_{eff}^2(\theta)$.

### T5.3: Re-compute Sagnac fazowy z directional chopper

v1 błąd: izotropowy $\Delta c/c$ → ~3×10⁴ SNR był overestimate. v2 poprawka:
**chopper differential** między dwoma orientacjami pętli interferometru:

- **Configuration A**: pętla w płaszczyźnie zawierającej $\vec\nabla\ln X$
  (pełny anizotropowy efekt, oba ramiona różnie zmodyfikowane)
- **Configuration B**: pętla prostopadła do $\vec\nabla\ln X$
  (BRAK efektu, oba ramiona equally unaffected)

Differential signal = $\Delta\phi_A - \Delta\phi_B$ (eliminuje common-mode
systematics). Compute order-of-magnitude SNR z realistic params:
- $|\vec\nabla\ln X|^2 \cdot L^2 \sim$ from φ.1 substrate gradient
- $|\beta_g| \sim O(0.1)$ z UV matching (Adams forces $\beta_g < 0$)
- $\Lambda \sim 100$ TeV (EFT cutoff)

PASS: sympy/numeric formula dla $\Delta\phi_{A-B}$ z explicit dependence
na orientation; SNR ~10²-10³ (realistic, **NIE** 10⁴).

### T5.4: Yukawa Greens dla tensor source

Static lab E∥B configuration: induced tensor source
$T^{\mu\nu}_{induced} \sim (\partial \ln X)^2 F^{\mu\rho}F^\nu_{\;\rho}$.
Solve linearized EOM dla efektywnego potential — Yukawa-like decay
z mass scale $\Lambda$:

$\Phi_{eff}(r) \sim \frac{|\beta_g|}{\Lambda^2} \cdot \frac{e^{-r\Lambda}}{r} \cdot$ (tensor projector)

PASS: sympy verify Yukawa form-factor + tensor projector consistent z
T4.3 effective metric.

### T5.5: 4 alt-tensor-L₅' falsification matrix

Cross-validate L₅'_a vs alt-tensor candidates pod kątem realnych obserwabli:

| Operator | c_eff anisotropic | Sagnac chopper | TOF directional | Cherenkov |
|----------|:---:|:---:|:---:|:---:|
| **L₅'_a** $(\partial\ln X)(\partial\ln X)F\cdot F$ | YES (cos²θ) | YES (A-B) | YES | safe |
| L₅'_b parity-odd $F\cdot\tilde F$ | birefringence (helicity) | NO (helicity) | YES | safe |
| L₅'_c $(\partial\partial \ln X)F\cdot F$ | YES (reduces to L₅'_a) | YES | YES | safe |
| L₅'_d scalar $\Box(\ln X) F^2$ | NO (scalar — v1 mistake) | NO | NO | n/a |

PASS: matrix shows L₅'_a uniquely produces tensor anisotropy w Sagnac chopper
(L₅'_b daje birefringence ale NIE Sagnac differential; L₅'_d zerowe — to było
fałszywe v1).

## Phase verdict logic

5/5 PASS → Phase 6 forward (corrected predictions TT19-TT23)
≤4/5 → re-evaluate eikonal/Sagnac derivation
