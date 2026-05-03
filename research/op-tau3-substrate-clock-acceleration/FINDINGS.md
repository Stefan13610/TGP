---
title: "FINDINGS — op-tau3-substrate-clock-acceleration"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-tau3-substrate-clock-acceleration
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-tau3-substrate-clock-acceleration

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `τ.3.Phase1` | `PASS` | — | — | — |
| `Phase2_results.md` | `τ.3.Phase2` | `PASS` | — | — | — |
| `Phase3_results.md` | `τ.3.Phase3` | `PASS-A5-PATCHED` | — | — | — |
| `Phase4_results.md` | `τ.3.Phase 4` | `PASS — CLOSED` | — | — | — |

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — Key results

> ## Key results
> 
> ### T1.1: L4 candidate forms
> Of 4 candidate gradient/F·F̃ couplings, **3 are scale-invariant**:
> - **L4_a = m_0 + α_g (∂_μ ln X)(∂^μ ln X) / Λ²** ← canonical (lowest-dim, derivative-only)
> - **L4_b = m_0 + β F·F̃ / Λ⁴** (direct F·F̃ coupling)
> - **L4_d = m_0 + η (E·B)² / Λ⁶** (dim-10 EFT, more suppressed)
> - L4_c = m_0 + γ ln(F·F̃) — REJECTED (diverges in vacuum F·F̃ → 0)
> 
> Sign of α_g is structurally free under φ.1 — both signs allowed.
> 
> ### T1.2: UV matching predicts α_g > 0
> Three independent UV channels:
> 1. **AS NGFP**: substrate kinetic-positive → Wilson coef of (∂ ln X)² ψ̄ψ ge...(truncated)

### `Phase4_results.md` — TL;DR

> ## TL;DR
> 
> **Audit B1** flagged the original Phase-1 sign claim α_g > 0 as resting on
> 3 weak channels (AS NGFP "attractive" unverified, heavy-mode magnitude ≠
> sign, BBN constraint not derivation). B1 explicitly required:
> 
> > "Explicit forward-amplitude calc; tylko Adams positivity v2 robust"
> 
> **τ.3.Phase 4** delivers exactly that:
> 
> $$
> \boxed{\;\frac{d^2 A(s, t{=}0)}{ds^2}\bigg|_{s=0}
>   \;=\; \frac{2\,\alpha_g}{\Lambda^4}
>   \;=\; \frac{4}{\pi}\!\int_{4m^2}^{\infty}\!ds'\,\frac{\sigma_{\text{tot}}(s')}{s'^3}
>   \;>\;0
>   \;\;\Longrightarrow\;\;
>   \alpha_g\;>\;0\;\text{strict, UV-independent}\;}
> $$
> 
> ...(truncated)

## Eksportowalne formuły (boxed)

**`B7v2_results.md`:**

$$
|\partial \ln X|(d) = \frac{J}{m_X} e^{-m_X d}
$$

**`B7v2_results.md`:**

$$
(\partial \ln X)^2(d) = \frac{J^2}{m_X^2} e^{-2 m_X d}
$$

**`B7v2_results.md`:**

$$
\int_V (\partial \ln X)^2 dV = A_{\text{edge}} \cdot \int_0^\infty \frac{J^2}{m_X^2} e^{-2 m_X d}\,dd = \frac{J^2 \cdot A_{\text{edge}}}{2\,m_X^3}
$$

**`B7v2_results.md`:**

$$
\langle (\partial \ln X)^2 \rangle_{\text{avg}} = \frac{1}{V_{\text{clock}}} \int_V (\partial \ln X)^2 dV = \frac{J^2 \cdot A_{\text{edge}}}{2\,m_X^3\,V_{\text{clock}}}
$$

**`B7v2_results.md`:**

$$
\frac{\langle (\partial \ln X)^2 \rangle_{\text{real}}}{(\partial \ln X)^2_{\text{naive}}} = \frac{A_{\text{edge}}}{2\,m_X\,V_{\text{clock}}}
$$

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa