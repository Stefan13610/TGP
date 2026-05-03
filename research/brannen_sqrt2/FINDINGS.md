---
title: "FINDINGS — R6: B=√2 analitycznie z ODE solitonowego"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/brannen_sqrt2
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — R6: B=√2 analitycznie z ODE solitonowego

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## TL;DR / Wynik / Verdict sections (cytaty)

### `README.md` — Wyniki α₅ (sesja 2026-04-17)

> ## Wyniki α₅ (sesja 2026-04-17)
> 
> ### **α₅ = 0.02751 (5 cyfr znaczących) — ZWERYFIKOWANE PERTURBACYJNIE**
> 
> Dwa niezależne podejścia dają zgodne wyniki:
> 
> **Podejście 1: Pomiar z pełnego ODE** (R_max do 3000, Richardson):
> - R_max=2000 (deg 3): α₅ = 0.027508
> - R_max=3000 (deg 3): α₅ = 0.027513
> - Richardson 1/R²: α₅ = 0.02752
> 
> **Podejście 2: Perturbacyjna formuła** (ODE f₂-f₅ do R_max=5000):
> ```
> α₅ = B...(truncated)

### `README.md` — Status checklist

> ## Status checklist
> 
> - [x] Ścieżka 4: F(φ) = A(φg₀)/A(g₀) — NIE stałe
> - [x] Łańcuch algebraiczny: Z₃ → 120° → K=2/3 → B=√2 — KOMPLETNY
> - [x] Best-fit tau: B = 1.4142 z |K-2/3| = 6×10⁻¹⁰
> - [x] η(δ) korekcja: c₁ = 0.72538, perturbacyjnie zrozumiana
> - [x] K(g₀^e) skan: K zależy od g₀^e i g₀^τ
> - [x] g₀^τ candidates: 2·g₀^e najlepszy (0.55%), ale nie dokładny
> - [x] r₂₁ nie jest uniwersalne — zależy od ...(truncated)

## Eksportowalne formuły (boxed)

**`README.md`:**

$$
s\,(f')^2 = \frac{d}{ds}[s\,f\,f'] + s\,f^2 + \tfrac{1}{2}\tfrac{d}{ds}[f^2]
$$

**`README.md`:**

$$
K_c^{(II)} = -A_1 - A_2 + A_3 - A_4
$$

**`README.md`:**

$$
\Phi_1(t) = \tfrac{1}{2}[\text{Ci}(3t) - \text{Ci}(t)] - \tfrac{\sin(t)+\sin(3t)}{4t}
$$

**`README.md`:**

$$
\Phi_2(t) = \tfrac{1}{2}[\text{Si}(3t) - \text{Si}(t)] - \tfrac{\cos(t)-\cos(3t)}{4t}
$$

**`README.md`:**

$$
\Phi_3(t) = \tfrac{3\sin(t) - \sin(3t)}{8t} + \tfrac{3}{8}[\text{Ci}(3t)-\text{Ci}(t)] + \tfrac{\cos(3t)-\cos(t)}{8t^2}
$$

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa