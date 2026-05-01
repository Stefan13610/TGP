---
title: "ψ.1.v2.Phase4 results — tensor L₅' structural derivation 5/5 PASS"
date: 2026-05-01
cycle: ψ.1.v2.Phase4
status: PASS
parent: "[[program.md]]"
tags:
  - TGP
  - psi1
  - phase4
  - correction
  - tensor-operator
  - results
---

# ψ.1.v2.Phase4 results — 5/5 FULL CASCADE (correction phase)

## Sub-test outcomes

| ID | Test | Result |
|----|------|--------|
| **T4.1** | Tensor candidate scan + φ.1 X→λX scale-invariance | ✅ PASS |
| **T4.2** | Formal proof: scalar Z(x)F² fails to modify light cones | ✅ PASS |
| **T4.3** | Effective optical metric derivation from L₅'_a | ✅ PASS |
| **T4.4** | Causality consistency + positivity bounds | ✅ PASS |
| **T4.5** | UV matching β_g sign for tensor operator (NEW Wilson coef) | ✅ PASS |

**Score: 5/5 → Phase 5 forward**

## Key results

### T4.1: L₅'_a CANONICAL uniquely identified

| Operator | φ.1 inv | Tensor | Parity-even | Irreducible | Status |
|----------|:---:|:---:|:---:|:---:|--------|
| **L₅'_a** $(\partial_\mu \ln X)(\partial_\nu \ln X) F^{\mu\rho}F^\nu_{\;\rho}$ | ✓ | ✓ | ✓ | ✓ | **CANONICAL** |
| L₅'_b $(\partial_\mu \ln X)(\partial_\nu \ln X) F^{\mu\rho}\tilde F^\nu_{\;\rho}$ | ✓ | ✓ | ✗ | ✓ | parity-odd helicity-discriminator |
| L₅'_c $(\partial_\mu \partial_\nu \ln X) F^{\mu\rho}F^\nu_{\;\rho}$ | ✓ | ✓ | ✓ | ✗ | reduces to L₅'_a via parts |
| L₅'_d $(\Box \ln X) F^2$ | ✓ | ✗ | ✓ | ✗ | **SCALAR — same v1 pathology**, no Δc |

L₅'_a unique among tensorowych dim-6 EFT, derivative-only φ.1 X→λX scale-invariant.

### T4.2: Formal proof — scalar Z(x)F² inert (v1 was wrong)

Lagrangian $L = -\tfrac{1}{4}Z(x)F^2$ z $Z > 0$:
- EOM: $\partial_\mu(Z F^{\mu\nu}) = Z\,\partial_\mu F^{\mu\nu} + (\partial_\mu Z)F^{\mu\nu} = 0$
- Eikonal leading order: $Z \cdot k^2 = 0 \Rightarrow k^2 = 0$ (since Z > 0)
- **Null cones determined by $\eta^{\mu\nu}$, INDEPENDENT of Z(x)**

sympy verification: null surface from $\eta^{\mu\nu}k_\mu k_\nu = 0$ matches null surface from $Z\cdot\eta^{\mu\nu}k_\mu k_\nu = 0$ EXACTLY → scalar Z(x) does **NOT** shift light cone.

**Interpretation:** scalar F² coupling = wave-function renormalization typu Bekenstein/Sandvik dilaton-photon coupling → zmienia efektywne $\alpha_{em}$ (przez $e_{eff} = e/\sqrt{Z}$) ale **NIE c**. To wyjaśnia dlaczego Webb/Murphy mierzą $\Delta\alpha/\alpha$ a nie $\Delta c/c$ w tych modelach.

**ψ.1.v1 negative result strukturalnie udokumentowany.**

### T4.3: Effective optical metric $g^{\mu\nu}_{eff}$

Lagrangian $L_5'_a$ z eikonal limit daje:
$$\boxed{g^{\mu\nu}_{eff} = \eta^{\mu\nu} + \frac{\xi}{\Lambda^2}\,n^\mu n^\nu, \qquad n^\mu = \partial^\mu \ln X}$$

gdzie $\xi$ jest funkcją bare coupling $\beta_g$ (variational principle daje $\xi = -2\beta_g$ przy konwencji $L_5'_a = +(\beta_g/\Lambda^2)\,S_{\mu\nu}F^{\mu\rho}F^\nu_{\;\rho}$).

Limit $n \to 0$: $g^{\mu\nu}_{eff} = \eta^{\mu\nu}$ (standard cone recovered) ✓.

### T4.4: Causality + positivity bound

Statyczny gradient substratu $n^\mu = (0, \vec{n})$:
$$g^{\mu\nu}_{eff} k_\mu k_\nu = \omega^2 - |\vec{k}|^2\big[1 - \xi |\vec{n}|^2 \cos^2\theta/\Lambda^2\big] = 0$$

gdzie $\theta$ = kąt między $\vec{k}$ a $\vec{\nabla}\ln X$:

$$\boxed{c_{eff}^2(\theta) = 1 - \xi |\vec{n}|^2 \cos^2\theta / \Lambda^2}$$

| θ | $c_{eff}^2$ | Interpretacja |
|---|-------------|---------------|
| 0 (parallel) | $1 - \xi n^2/\Lambda^2$ | maksymalna modyfikacja |
| π/2 (perpendicular) | 1 | brak efektu (transversely) |

**Subluminal everywhere** (c_eff ≤ c_0) wymaga $\xi \geq 0$. W konwencji $\xi = -2\beta_g$ → **β_g ≤ 0**.

**Vacuum Cherenkov constraint**: gdyby $c_{eff} > c_0$, SM particles relatywistyczne ($v_{SM} \approx c_0$) wypromieniowywałyby photons → niestabilność próżni. Wymaganie subluminalności **wymusza** sign β_g.

**No closed timelike curves**: dla statycznego $n^\mu$ przestrzennego, sygnatura efektywnej metryki (+,-,-,-) zachowana przy $|\xi|n^2/\Lambda^2 \ll 1$. CTC nie powstają w lab E∥B configurations.

### T4.5: UV matching β_g sign (3 channels, Channel C decisive)

| Channel | Argument | Wynik |
|---------|----------|-------|
| **A: AS NGFP** | UV fixed-point Wilson coef dla tensor klas operatorów | sign UNDETERMINED bez explicit AS+matter-tensor calculation |
| **B: heavy-mode 1-loop** | $\beta_g \sim -(Q_f^2 m_f^2)/(48\pi^2 \Lambda^2)$ × tensor-projection | sign suggestywnie negative (zgodne z Channel C) |
| **C: Adams et al positivity** [arXiv:hep-th/0602178] | causality + analyticity forward amplitude → forced subluminality | **DECISIVE: β_g < 0 strict** |

**Synthesis:** Adams-Arkani-Hamed-Dubovsky-Nicolis-Rattazzi positivity bound (UV-independent) **wymusza** $\beta_g < 0$. Channels A, B konsystentne z tym (nie wymuszają niezależnie).

**Sign convention v2:** $|\beta_g| > 0$ z explicit Lagrangian:
$$L_5'_a = -\frac{|\beta_g|}{\Lambda^2}\,(\partial_\mu \ln X)(\partial_\nu \ln X)\,F^{\mu\rho}F^\nu_{\;\rho}$$

→ $\xi = +2|\beta_g| > 0$, subluminal $c_{eff} \leq c_0$.

## Phase verdict

**ψ.1.v2.Phase 4 PASS (FULL CASCADE 5/5) → Phase 5 forward**

Strukturalne wyniki Phase 4:
- L₅'_a $(\partial_\mu \ln X)(\partial_\nu \ln X)F^{\mu\rho}F^\nu_{\;\rho}$ CANONICAL uniquely identified (tensor + parity-EVEN + irreducible + φ.1-inv)
- **Formal proof**: skalarny Z(x)F² nie modyfikuje stożka świetlnego (sympy LOCK on null cone equality) — ψ.1.v1 NEGATIVE result formally documented
- Effective optical metric $g^{\mu\nu}_{eff} = \eta^{\mu\nu} + (\xi/\Lambda^2)n^\mu n^\nu$ derived
- **Anizotropowy c_local(θ)**: parallel-to-gradient maximally slowed, perpendicular unchanged
- **Positivity bound** β_g < 0 forced by Adams et al causality (UV-independent) — sign FIXED
- Light **SPOWALNIA** w kierunku gradientu (Δc/c < 0 for k∥∇lnX), brak vacuum Cherenkov

**Cross-cycle implication:** ψ.1 jest teraz **anizotropowym** sub-leading kanałem σ.1 (skalarny scalar c shift przy ∂lnX≠0 BYŁ fałszywy; rzeczywisty efekt jest **kierunkowo-zależny**). σ.1 helicity-dependent dispersion preserved leading order; ψ.1 tensor channel dodaje **kierunkową anizotropię scalar c shift** w sub-leading O((∂lnX)²).
