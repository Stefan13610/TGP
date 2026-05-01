---
title: "ψ.1.v2.Phase4 setup — tensor L₅' structural derivation + causality (correction phase)"
date: 2026-05-01
cycle: ψ.1.v2.Phase4
status: SETUP
parent: "[[program.md]]"
tags:
  - TGP
  - psi1
  - phase4
  - correction
  - tensor-operator
  - setup
---

# ψ.1.v2.Phase 4 setup — tensor L₅' structural derivation (5 sub-tests)

## Reason for correction phase

**ψ.1.v1 (Phase 1-3) negative structural result**: skalarny operator
$L_5 = -(1/4)(\beta_g/\Lambda^2)(\partial \ln X)^2 F^2$ jest **wave-function
renormalization** typu Bekenstein/Sandvik dilaton-photon coupling. Po
redefinicji $A'_\mu = Z^{1/2}A_\mu$ struktura stożków świetlnych pozostaje
**niezmieniona** (Z>0). Skalarny F² → "varying $\alpha_{em}$" (przez
$e_{eff} = e/\sqrt{Z}$), **NIE** "varying c". External agent's critique
(2026-05-01) trafny — Phase 1 T1.3 + Phase 2 T2.1 sympy LOCK weryfikował
algebrę Taylor expansion, **nie fizyczne znaczenie**. Sagnac SNR ~3×10⁴
**artefakt** błędnej Δc/c.

**Naprawa**: tensorowy operator
$$L_5' = \frac{\beta_g}{\Lambda^2}\,(\partial_\mu \ln X)(\partial_\nu \ln X)\,F^{\mu\rho}F^\nu_{\;\rho}$$
wyróżnia kierunek $n^\mu = \partial^\mu \ln X$ → efektywna metryka optyczna
$g^{\mu\nu}_{eff} = \eta^{\mu\nu} + \xi (n^\mu n^\nu)/\Lambda^2$ → realna
modyfikacja stożka świetlnego (anizotropowa, zależna od kąta między $\vec{k}$
a $\vec{\nabla}\ln X$).

## Sub-tests

### T4.1: Tensor candidate scan + φ.1 X→λX scale-invariance

Test 4 candidate tensorowych operatorów dim-6 EFT:
- **L₅'_a** = $(\partial_\mu \ln X)(\partial_\nu \ln X) F^{\mu\rho}F^\nu_{\;\rho}$ — symmetric tensor coupling
- **L₅'_b** = $(\partial_\mu \ln X)(\partial_\nu \ln X) F^{\mu\rho}\tilde F^\nu_{\;\rho}$ — parity-odd tensor (helicity)
- **L₅'_c** = $(\partial_\mu \partial_\nu \ln X) F^{\mu\rho}F^\nu_{\;\rho}$ — second-derivative
- **L₅'_d** = $(\partial^2 \ln X) F^{\mu\nu}F_{\mu\nu}$ — scalar with 2nd derivative

Filter: φ.1 X→λX invariance + dim-6 EFT + irreducibility + tensor structure.

PASS criterion: at least 1 candidate uniquely identified as canonical.

### T4.2: Formal proof — scalar Z(x)F² fails to modify light cones

Sympy/symbolic proof:
- Lagrangian $L = -(1/4)Z(x)F^2$
- After field redefinition $A'_\mu = Z^{1/2}A_\mu$ → standard kinetic + induced terms involving $\partial Z$
- Principal symbol of EOM: characteristics determined by $\eta^{\mu\nu}$, NOT modified by $Z(x) > 0$
- Conclusion: **scalar F² coupling renormalizes $\alpha_{em}$ (varying-α), NOT c**

PASS: structural proof that v1 mechanism gives null Δc, with sympy LOCK on the redefinition.

### T4.3: Effective optical metric derivation from L₅'_a

For Lagrangian $L = -(1/4)F^2 + (\beta_g/\Lambda^2) S_{\mu\nu}F^{\mu\rho}F^\nu_{\;\rho}$
with $S_{\mu\nu} = (\partial_\mu \ln X)(\partial_\nu \ln X)$:

Eikonal limit $A_\mu = a_\mu \exp(i\phi/\hbar)$, $\hbar \to 0$:
- Principal symbol $D(k)$ from kinetic operator
- Dispersion relation: $g^{\mu\nu}_{eff} k_\mu k_\nu = 0$
- Derive $g^{\mu\nu}_{eff} = \eta^{\mu\nu} + \xi (n^\mu n^\nu)/\Lambda^2$ where $\xi = O(\beta_g)$

PASS: sympy LOCK on $g^{\mu\nu}_{eff}$ formula.

### T4.4: Causality consistency + positivity bounds

For anisotropic dispersion $\omega^2 = k^2[1 + \xi(\hat{n}\cdot\hat{k})^2 (\partial \ln X)^2/\Lambda^2]$:

- **Subluminal**: requires $\xi < 0$ (so $c_{eff} \leq c_0$ in all directions) — otherwise vacuum Cherenkov by SM particles
- **Closed timelike curves**: check no CTC formation in lab E∥B configurations
- **Effective metric Lorentzian**: ensure signature preserved
- **Cross-sector consistency**: σ.1 helicity birefringence + τ.3 mass shift unaffected

PASS: positivity bound $\beta_g < 0$ derived; no CTC; Lorentzian preserved at $|\beta_g|(\partial\ln X)^2/\Lambda^2 \ll 1$.

### T4.5: UV matching β_g sign for tensor operator (NEW Wilson coef)

Tensor $L_5'$ has DIFFERENT Wilson coefficient than scalar L₅. Re-do UV matching:
- **Channel A**: AS NGFP — sign of tensor coefficient at fixed point (different operator class)
- **Channel B**: heavy-mode 1-loop with tensor structure — generic sign analysis
- **Channel C**: cosmological/BBN consistency for anisotropic c

PASS: sign of tensor β_g consistent with positivity bound from T4.4 (subluminal).

## Phase verdict logic

5/5 PASS → Phase 5 forward
≤4/5 → re-evaluate operator choice or additional candidates
