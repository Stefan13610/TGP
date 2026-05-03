---
title: "FINDINGS — op-psi1-substrate-light-acceleration"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-psi1-substrate-light-acceleration
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-psi1-substrate-light-acceleration

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `ψ.1.Phase1` | `INVALIDATED` | — | — | — |
| `Phase2_results.md` | `ψ.1.Phase2` | `INVALIDATED` | — | — | — |
| `Phase3_results.md` | `ψ.1.Phase3` | `INVALIDATED` | — | — | — |
| `Phase4_results.md` | `ψ.1.v2.Phase4` | `PASS` | — | — | — |
| `Phase5_results.md` | `ψ.1.v2.Phase5` | `PASS` | — | — | — |
| `Phase6_results.md` | `ψ.1.v2.Phase6` | `PASS` | — | — | — |
| `Phase7_results.md` | `ψ.1.v3.Phase 7` | `PASS — CLOSED` | — | — | — |

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — Key results

> ## Key results
> 
> ### T1.1: L₅_a CANONICAL uniquely identified
> 
> | Operator | φ.1 inv | Scalar | Parity-even | Status |
> |----------|:---:|:---:|:---:|--------|
> | **L₅_a (∂lnX)²·F²** | ✓ | ✓ | ✓ | **CANONICAL** |
> | L₅_b (∂lnX)²·F·F̃ | ✓ | ✗ | ✗ | parity-odd, helicity-discriminator |
> | L₅_c (□lnX)·F² | ✓ | ✓ | ✓ | reducible to L₅_a via parts |
> | L₅_d ln(X)·F² | ✗ | ✓ | ✓ | DILATON, breaks φ.1 X→λX |
> 
> L₅_a = $-(1/4)(\beta_g/\Lambda^2)(\partial_\mu \ln X)(\partial^\mu \ln X)\,F_{\nu\rho}F^{\nu\rho}$ is unique scale-invariant scalar irreducible candidate.
> 
> ### T1.2: β_g sign 3-channel agreement (β_g >...(truncated)

### `Phase4_results.md` — Key results

> ## Key results
> 
> ### T4.1: L₅'_a CANONICAL uniquely identified
> 
> | Operator | φ.1 inv | Tensor | Parity-even | Irreducible | Status |
> |----------|:---:|:---:|:---:|:---:|--------|
> | **L₅'_a** $(\partial_\mu \ln X)(\partial_\nu \ln X) F^{\mu\rho}F^\nu_{\;\rho}$ | ✓ | ✓ | ✓ | ✓ | **CANONICAL** |
> | L₅'_b $(\partial_\mu \ln X)(\partial_\nu \ln X) F^{\mu\rho}\tilde F^\nu_{\;\rho}$ | ✓ | ✓ | ✗ | ✓ | parity-odd helicity-discriminator |
> | L₅'_c $(\partial_\mu \partial_\nu \ln X) F^{\mu\rho}F^\nu_{\;\rho}$ | ✓ | ✓ | ✓ | ✗ | reduces to L₅'_a via parts |
> | L₅'_d $(\Box \ln X) F^2$ | ✓ | ✗ | ✓ | ✗ | **SCALA...(truncated)

### `Phase5_results.md` — Key results

> ## Key results
> 
> ### T5.1: Eikonal dispersion relation LOCKED
> 
> Z $g^{\mu\nu}_{eff} = \eta^{\mu\nu} + (\xi/\Lambda^2) n^\mu n^\nu$ i konwencji
> sygnaturalnej (+,-,-,-):
> 
> $$\boxed{g^{\mu\nu}_{eff}\, k_\mu k_\nu = (\omega^2 - |\vec{k}|^2) + \frac{\xi}{\Lambda^2}(n\cdot k)^2 = 0}$$
> 
> gdzie $n\cdot k = n^\mu k_\mu = n^0\omega - \vec{n}\cdot\vec{k}$.
> 
> sympy diff vs. expected = 0 (po fix Minkowski signature contraction —
> n·k = n^μ k_μ z metryką, nie Euclidean dot).
> 
> ### T5.2: Anisotropic $c_{local}(\theta)$ LOCKED
> 
> Statyczny gradient $n^\mu = (0, \vec{n})$, $\vec{n} = \vec\nabla\ln X$, $|\vec{n}| = n_{m...(truncated)

### `Phase6_results.md` — Key results

> ## Key results
> 
> ### TT19: Sagnac chopper differential
> 
> $$\boxed{\Delta\phi_{A-B} = \frac{4\pi L_{arm}}{\lambda_\gamma}\cdot|\beta_g|\cdot\frac{|\vec\nabla\ln X|^2}{\Lambda^2}}$$
> 
> - Lab realistic ($L=1$ m, $\lambda=1$ μm, $|\beta_g|=0.1$, $\Lambda=100$ TeV, $|\vec\nabla\ln X|=10^{-3}$ m⁻¹): $\Delta\phi \sim 5\times 10^{-36}$ rad
> - SNR (1 month, shot-floor 6×10⁻¹³ rad/√Hz): ~8×10⁻²⁴ (**sub-detection by ~23 OOM**)
> - **Falsifier**: positive lab Sagnac chopper SNR > 10⁻³ (without external $\vec\nabla\ln X$ amplification) FALSIFIES ψ.1.v2
> 
> **v1 corrected**: porównanie z fałszywym v1 SNR ~3×10⁴ → ψ.1...(truncated)

### `Phase7_results.md` — TL;DR

> ## TL;DR
> 
> **Audit C8** demanded systematic Hilbert-series-style enumeration of
> dim-6 EFT operators for the ψ.1 photon-substrate sector — replacing the
> manual 4-operator scan from ψ.1.v2.Phase 4 T4.1 with a structured catalog
> that includes IBP/Bianchi/EOM reduction relations.
> 
> **Phase 7** delivers exactly that. Result:
> 
> $$
> \boxed{\;
> \mathcal{B}_{\psi.1\text{-v3}}^{\dim\text{-}6}
> \;=\;
> \{
> \,L_5'^{(a)}{=}(\partial_\mu \ln X)(\partial_\nu \ln X)F^{\mu\rho}F^\nu{}_\rho\,,
> \;L_5'^{(b)}{=}(\partial_\mu \ln X)(\partial_\nu \ln X)F^{\mu\rho}\widetilde{F}^\nu{}_\rho\,
> \}\;}
> $$
> 
> Two-element canonical bas...(truncated)

## Eksportowalne formuły (boxed)

**`Phase1_results.md`:**

$$
\frac{\Delta c}{c_0} \;=\; -\frac{\beta_g}{2\Lambda^2}(\partial \ln X)^2 \;<\; 0 \qquad (\text{w obszarze gradientu})
$$

**`Phase1_results.md`:**

$$
L_{em} + L_5 = -\frac{1}{4}\Big[1 + \beta_g\frac{(\partial \ln X)^2}{\Lambda^2}\Big]F^2
$$

**`Phase1_results.md`:**

$$
\Rightarrow\quad c_{local} = \frac{c_0}{\sqrt{1 + \varepsilon}}, \quad \varepsilon \equiv \beta_g(\partial\ln X)^2/\Lambda^2
$$

**`Phase1_results.md`:**

$$
\boxed{\frac{\Delta c}{c_0} \;=\; -\frac{\beta_g}{2\Lambda^2}(\partial \ln X)^2}
$$

**`Phase2_results.md`:**

$$
L_{em} + L_5 = -\frac{1}{4}\Big[1 + \beta_g\frac{(\partial\ln X)^2}{\Lambda^2}\Big] F^2
$$

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa