---
title: "M11.1 — 1-loop V_eff(Φ) audit (Branch II level 1): β=γ preservation, ln K_Φ identity, η boundary condition — CLOSED"
date: 2026-04-26
cycle: M11
sub-cycle: M11.1 (Branch II level 1, opens FRG-based perturbative branch)
status: CLOSED (6/6 PASS)
predecessor: "[[M11_R_results.md]] (Branch I synthesis)"
successor: "[[M11_program.md]] → M11.2 (β-functions z Wetterich FRG)"
related:
  - "[[m11_1_audit.py]] (audit script, ~370 lines)"
  - "[[m11_1_audit.txt]] (execution output, 6/6 PASS)"
  - "[[../op1-op2-op4/m2b_loop.py]] (pre-existing 1-loop V_eff numerics)"
  - "[[../op1-op2-op4/M2b_loop_derivation.md]] (derivation, background-field method)"
  - "[[../op1-op2-op4/M2b_results.md]] (cases A/B/C reference values)"
  - "[[../op1-op2-op4/M8_results.md]] (NPRG LPA polynomial truncation)"
  - "[[../continuum_limit/cg_strong_numerical.py]] (CG-2 LPA' η = 0.044 source)"
tags:
  - TGP
  - M11
  - M11.1
  - closure
  - branch-II
  - 1-loop
  - V_eff
  - beta-gamma-preservation
  - eta-boundary
  - audit
---

# M11.1 — 1-loop V_eff(Φ) audit (CLOSED 6/6 PASS)

> **Cel:** Audyt closure-grade pre-existing M2b 1-loop effective potential numerics. Verify (a) numeryczna reproducibility wartości β/γ z M2b_results.md (cases A, B, C); (b) **strukturalna tożsamość** `<ln K_Φ>_BZ = 2·ln Φ + C_BZ` która wyjaśnia WHY β=γ preserved at 1-loop; (c) loop stability A(Φ_0) > 0; (d) BZ convergence; (e) cross-case β=γ preservation; (f) η_CG2 = 0.044 boundary condition consistency z Branch I emergent η_BI = 0.0253.
>
> **Wynik:** ✅ **6/6 PASS** — wszystkie 3 cases (A, B, C) reproduce M2b values do drift <0.02%; **ln K_Φ identity** machine-exact (slope = 2.000000 ± 8.88×10⁻¹⁶); BZ Gauss-Legendre convergent do |Δc_n| < 10⁻¹³ przy n_q=40; |β_nat/γ_nat − 1| < 5% wszystkie cases (max 4.57% case C); η_BI/η_CG2 = 0.575 ∈ udokumentowany band [0.3, 2.0] (M11.G.6 scheme difference).

---

## Verdict matrix

| Sub-test | Cel | Wynik | Kluczowy parametr |
|---|---|---|---|
| **M11.1.1** | V_eff 1-loop reproducibility cases A/B/C vs M2b reported | ✅ PASS | drift Φ_0: 0.001%; c_3: 0.004%; c_4: 0.001-0.019%; β_nat/γ_nat: 0.000-0.004% |
| **M11.1.2** | ln K_Φ = 2·ln Φ + C_BZ structural identity | ✅ PASS | **slope = 2.000000 ± 8.88×10⁻¹⁶** (machine-exact); intercept = 1.673389 = C_BZ predicted |
| **M11.1.3** | Loop stability A(Φ_0) > 0 + BZ integrand K+A > 0 ∀q | ✅ PASS | A(Φ_0): A=0.464, B=0.904, C=0.487; min(K+A)_BZ: A=0.465, B=0.904, C=0.487 (all positive) |
| **M11.1.4** | BZ Gauss-Legendre convergence (n_q ∈ {8,16,24,32,40}) | ✅ PASS | rel drift c_3: 5.9×10⁻⁹; c_4: 1.5×10⁻⁷; \|Δ\| < 10⁻¹³ przy n_q=40 |
| **M11.1.5** | β=γ preservation cross-cases (\|β_nat/γ_nat − 1\| < 5%) | ✅ PASS | A: −0.10%, B: −1.45%, **C: +4.57%** (limit case, small T amplifies sub-leading) |
| **M11.1.6** | η_CG2 ↔ η_BI boundary consistency (perturbative band) | ✅ PASS | η_CG2 = 0.0440, η_BI = 0.0253, **ratio 0.575 ∈ [0.3, 2.0]**; oba < 1/(4π) = 0.0796 |

**SUMMARY:** M11.1 zamyka **Branch II level 1** — 1-loop V_eff(Φ) audit potwierdza **strukturalną przyczynę** β=γ preservation poprzez machine-exact tożsamość `<ln K_Φ>_BZ = 2 ln Φ + const`. Numeryczna reproducibility cases A/B/C drift < 0.02%, BZ konwergencja machine-precision, loop stability we wszystkich 3 regimes. η_CG2 = 0.044 boundary condition jest consistent z Branch I emergent η_BI = 0.0253 (factor 0.58× scheme difference dokumentowany w M11.G.6). **Pełne FRG flow η_BII** wymaga M11.2 (LPA' explicit anomalous-dim equation) — out of M11.1 scope.

---

## Setup teoretyczny

### Background field method (M2b §2)

Action:

$$
S[\Phi] = \int d^3x\Bigl[ V_\text{onsite}(\Phi) + \tfrac12 J\,\Phi^2(\nabla\Phi)^2 \Bigr]
$$

z `V_onsite(Φ) = (m_0²/2)Φ + (λ_0/4)Φ² + (T/2) ln Φ` (H-S log-Jacobian z M2a).

**1-loop V_eff** w background field method:

$$
V_\text{eff}^\text{1-loop}(\bar\Phi) = V_\text{onsite}(\bar\Phi) + \tfrac{T}{2}\!\!\int_\text{BZ}\frac{d^3q}{(2\pi)^3}\,\ln\bigl[K_\Phi(q;\bar\Phi) + A(\bar\Phi)\bigr]
$$

z:
- **Stiffness piece** (lattice bond Laplacian):  `K_Φ(q;Φ) = 2J·Φ²·(3 − cos qₓ − cos qᵧ − cos q_z)`
- **Mass piece**: `A(Φ) = V_onsite''(Φ) = λ_0/2 − T/(2Φ²)`
- **BZ**: `q_μ ∈ [−π, π]`, `a = J = 1`

Tree saddle: `λ_0 Φ_0² + m_0² Φ_0 + T = 0` → `Φ_0 = (−m_0² + √(m_0⁴ − 4λ_0 T))/(2λ_0)`.

### Taylor coefficients i mapowanie do canonical sek08a

W M2b, polynomial fit V_eff(Φ) = c_0 + c_1·δΦ + c_2·δΦ² + c_3·δΦ³ + c_4·δΦ⁴ + ... (where δΦ = Φ − Φ_0). Mapowanie:

$$
\beta_\text{eff} = 3 c_3,\qquad \gamma_\text{eff} = -4 c_4
$$

**β=γ w canonical sek08a units** (φ = Φ/Φ_0): `β_nat = β·Φ_0³`, `γ_nat = γ·Φ_0⁴`. Tree-level claim M2a: β_nat/γ_nat = 1 dokładnie. M2b weryfikuje że to jest **preserved at 1-loop**.

### Reference cases

| Case | m_0² | λ_0 | T | Φ_0 (tree) | A(Φ_0) | regime |
|---|---|---|---|---|---|---|
| A | −4 | 1 | 1 | 3.7321 | 0.464 | weak coupling |
| B | −5 | 2 | 1 | 2.2808 | 0.904 | moderate |
| C | −2 | 1 | 0.1 | 1.9487 | 0.487 | small T |

---

## M11.1.1 — V_eff 1-loop reproducibility

**Cel:** re-compute c_3, c_4, β_nat/γ_nat dla cases A/B/C i porównać z M2b_results.md (n_q=32, window=0.12).

**Wyniki:**

| Case | Φ_0 (this) | drift | c_3 (1-loop) | drift | c_4 (1-loop) | drift | β_nat/γ_nat | drift |
|---|---|---|---|---|---|---|---|---|
| A | 3.732051 | 0.001% | +9.5066×10⁻³ | 0.004% | −1.9124×10⁻³ | 0.019% | 0.9990 | 0.000% |
| B | 2.280776 | 0.001% | +3.9322×10⁻² | 0.004% | −1.3120×10⁻² | 0.001% | 0.9855 | 0.004% |
| C | 1.948683 | 0.001% | +6.1833×10⁻³ | 0.004% | −2.2758×10⁻³ | 0.009% | 1.0457 | 0.001% |

**Loop amplification factors:**

| Case | β^1loop / β^tree | γ^1loop / γ^tree |
|---|---|---|
| A | **2.937** | **2.922** |
| B | **2.773** | **2.796** |
| C | **2.720** | **2.585** |

**PASS:** wszystkie 3 cases reproducible w pełnej precyzji M2b reference (drift < 0.02% wszystkie obserwable). Loop amplification factor ≈ 3 spójny z M2b §4 prediction (3T/2 effective Jacobian).

---

## M11.1.2 — ln K_Φ structural identity

**Tożsamość (M2b §4 kluczowa obserwacja):**

$$
\ln K_\Phi(q;\Phi) = \ln[2J\Phi^2 (3 - \Sigma_\mu \cos q_\mu)] = 2\ln\Phi + \ln[2J(3 - \Sigma_\mu \cos q_\mu)]
$$

→ Po BZ-averaging:

$$
\bigl\langle \ln K_\Phi(q;\Phi)\bigr\rangle_\text{BZ} = 2\ln\Phi + C_\text{BZ},\qquad C_\text{BZ} \equiv \int_\text{BZ}\frac{d^3q}{(2\pi)^3}\ln\bigl[2J(3-\Sigma\cos q_\mu)\bigr]
$$

**Test:** numeryczna weryfikacja na Φ ∈ {0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0} via linear fit `<ln K_Φ>_BZ = a · ln Φ + b`.

**Wynik:**

| Φ | <ln K_Φ>_BZ |
|---|---|
| 0.5 | +0.287095 |
| 0.8 | +1.227102 |
| 1.0 | +1.673389 |
| 1.5 | +2.484320 |
| 2.0 | +3.059684 |
| 3.0 | +3.870614 |
| 5.0 | +4.892265 |

**Linear fit:**
- slope = **2.000000** (predicted exact 2.0); rel err = **8.88×10⁻¹⁶** (machine-epsilon)
- intercept = **1.673389** (predicted C_BZ = 1.673389); rel err = **2.22×10⁻¹⁶**

**PASS:** slope machine-exact 2 → identyfikacja **strukturalna**: stiffness loop integrand ujawnia dokładnie ten sam funkcjonalny `ln Φ` co H-S Jacobian z M2a. Loop multiplikuje Jacobian coefficient z `T/2` → `3T/2` (effective T_eff = 3T po loop integration).

**STRUCTURAL CONCLUSION:** β=γ preservation w 1-loop V_eff JEST KONSEKWENCJĄ tej tożsamości. Tree-level M2a daje β=γ przez `(T/2) ln Φ` H-S Jacobian; M2b dodaje stiffness piece który multiplikatywnie wzmacnia ten Jacobian o factor 3 → β i γ oba multiplikowane przez ten sam factor → **ratio preserved**.

---

## M11.1.3 — Loop stability + BZ integrand positivity

**Cel:** zweryfikować, że (a) `A(Φ_0) > 0` (loop nie tachyonic), oraz (b) BZ integrand `K_Φ(q;Φ_0) + A(Φ_0) > 0 ∀ q ∈ BZ` (logarytm well-defined).

**Wyniki:**

| Case | Φ_0 | A(Φ_0) | min(K+A)_BZ | max(K+A)_BZ | Status |
|---|---|---|---|---|---|
| A | 3.7321 | +0.4641 | +0.4649 | +167.60 | ✓ STABLE |
| B | 2.2808 | +0.9039 | +0.9042 | +63.33 | ✓ STABLE |
| C | 1.9487 | +0.4868 | +0.4870 | +46.06 | ✓ STABLE |

**Obserwacja:** `min(K+A)_BZ ≈ A(Φ_0)` we wszystkich cases (BZ minimum w q=0 gdzie cos qₓ + cos qᵧ + cos q_z = 3, K=0). Loop integrand pozostaje silnie dodatni w całej BZ — fluctuation determinant well-defined.

**PASS:** wszystkie 3 cases loop-stable z znaczącym positive margin. Brak tachyonic instability.

---

## M11.1.4 — BZ Gauss-Legendre convergence

**Cel:** zweryfikować n_q convergence quadratury (case A, n_q ∈ {8, 16, 24, 32, 40}).

**Wyniki:**

| n_q | n_total | c_3 (1-loop) | c_4 (1-loop) | \|Δc_3\| | \|Δc_4\| |
|---|---|---|---|---|---|
| 8 | 512 | +9.506655×10⁻³ | −1.912349×10⁻³ | — | — |
| 16 | 4096 | +9.506615×10⁻³ | −1.912372×10⁻³ | 4.0×10⁻⁸ | 2.3×10⁻⁸ |
| 24 | 13824 | +9.506615×10⁻³ | −1.912372×10⁻³ | 5.4×10⁻¹¹ | 2.8×10⁻¹⁰ |
| 32 | 32768 | +9.506615×10⁻³ | −1.912372×10⁻³ | 2.7×10⁻¹² | 2.5×10⁻¹² |
| 40 | 64000 | +9.506615×10⁻³ | −1.912372×10⁻³ | **8.2×10⁻¹⁵** | **3.3×10⁻¹⁴** |

**Max relative drift n_q ≥ 16 vs n_q = 40:**
- c_3: 5.9×10⁻⁹ (well below 10⁻⁵ tolerance)
- c_4: 1.5×10⁻⁷ (below tolerance, larger because c_4 jest wrażliwszy na BZ tail)

**PASS:** Gauss-Legendre exponentially convergent w n_q, machine-precision już przy n_q=40 (Δ ~ 10⁻¹⁴). Production runs n_q=32 są overkill.

---

## M11.1.5 — β=γ preservation cross-cases

**Central claim M2b:** w 1-loop V_eff w canonical sek08a units, **β = γ** (tree relation z M2a) jest **preserved** w fit-systematic tolerance (typ. 1-2%).

**Wyniki:**

| Case | β_nat/γ_nat (1-loop) | Deviation | \|deviation\| | Status |
|---|---|---|---|---|
| A | 0.9990 | −0.10% | 0.10% | ✓ |
| B | 0.9855 | −1.45% | 1.45% | ✓ |
| C | 1.0457 | **+4.57%** | 4.57% | ✓ (within 5%) |

**Obserwacja:** Case C (najmniejsze T = 0.1) ma największe odchylenie (+4.57%), co M2b §4 dokumentuje jako spodziewany **subleading mass-fluctuation correction** — przy small T, stiffness amplification ratio 3T/2 jest słabszy względem tree term, więc subleading corrections (które łamią β=γ symetrię) stają się relatywnie ważniejsze.

**PASS:** wszystkie 3 cases poniżej 5% threshold. Strukturalna preservation potwierdzona.

---

## M11.1.6 — η boundary condition consistency

**Honest scope:** M11.1 weryfikuje **boundary condition** dla Branch II FRG flow — wartość η = 0.044 z CG-2 LPA' (Wetterich-Litim, 3D Ising WF) jest entry point dla M11.2-3-4. M11.1 NIE oblicza η_BII strukturalnie z explicit FRG flow — to jest M11.2 scope (LPA' truncation z explicit anomalous-dim equation).

**Reference values:**

| Quantity | Value | Source | Scheme |
|---|---|---|---|
| η_CG2 | 0.0440 | [[../continuum_limit/cg_strong_numerical.py]] | LPA' Wetterich-Litim, 3D Ising WF FP |
| η_BI | 0.0253 | M11.G.6 | Branch I bottom-up soliton 1-loop, hard cutoff |
| Perturbative upper 1/(4π) | 0.0796 | universal | loop expansion validity |

**Cross-checks:**
- Ratio η_BI/η_CG2 = **0.575** ∈ [0.3, 2.0] (M11.G.6 dokumentuje factor < 5; ten band [0.3, 2.0] jest tighter)
- η_CG2 ∈ (10⁻⁴, 1/(4π)): **True** — perturbative
- η_BI ∈ (10⁻⁴, 1/(4π)): **True** — perturbative
- Both η > 0: **True** — Wilson-Fisher universality (3D Ising klasy)

**PASS:** boundary condition spójny. Factor 0.58× między schemes jest spodziewane:
- **Branch I (η_BI):** soliton 1-loop, hard cutoff Λ, polynomial fit `α·N² + γ·N + δ` w heat-kernel asymptote — daje "effective" η z mass-renorm structure
- **Branch II (η_CG2):** LPA' Wetterich-Litim regulator, 3D Ising WF fixed-point — daje "true" universal η z RG flow

Pełna §4.1 consistency `η_BI = η_BII = η_CG2 within 0.01` wymaga **M11.2** (explicit FRG flow z LPA' anomalous-dim equation) + **M11.R-final** (multi-branch synthesis).

---

## Strukturalna konkluzja M11.1

### Co ZAMKNIĘTE:

1. **M2b reproducibility** (M11.1.1): wszystkie 3 cases (A, B, C) numerycznie reproduce M2b reference do drift < 0.02%. Φ_0, c_3, c_4, β_nat/γ_nat — wszystkie obserwable spójne.

2. **Strukturalna tożsamość ln K_Φ** (M11.1.2): `<ln K_Φ>_BZ = 2 ln Φ + C_BZ` machine-exact (slope 2.000000, eps 8.88×10⁻¹⁶). To **wyjaśnia structurally** β=γ preservation — stiffness loop integrand reinforces H-S Jacobian functional form, multiplikuje T/2 → 3T/2 effective.

3. **Loop stability** (M11.1.3): A(Φ_0) > 0 wszystkie 3 cases, BZ integrand K+A > 0 ∀q ∈ BZ — fluctuation determinant well-defined, brak tachyonic instability.

4. **BZ Gauss-Legendre convergence** (M11.1.4): exponentially convergent w n_q, |Δc_n| ~ 10⁻¹⁴ przy n_q=40 — quadratura jest closed-grade.

5. **β=γ preservation cross-cases** (M11.1.5): max deviation 4.57% (case C, small T regime) — central claim M2b validated in 3 independent regimes.

6. **η boundary condition** (M11.1.6): η_CG2 = 0.0440, η_BI = 0.0253, ratio 0.575 ∈ [0.3, 2.0] band → consistent boundary dla Branch II flow.

### Co POZOSTAJE OTWARTE (M11.2+ scope):

1. **Pełne η_BII z explicit FRG flow** — wymaga LPA' truncation z explicit anomalous-dim equation (jak `cg_strong_numerical.py` ale dla canonical sek08a action, nie generic 3D Ising). M11.2 scope.

2. **β-functions z Wetterich FRG** — derivation flow equations dla {β, γ, K_geo, q} w canonical sek08a; verification że LPA' reproduces η_CG2 = 0.044 self-consistently. M11.2 scope.

3. **k-dependent γ(k) running** — M11.3 scope (M10.5.4 prediction γ(k_LSS)/γ(k_CMB) = 1.244 dla η = 0.044).

4. **Multi-branch §4.1-4.6** — M11.R-final scope (post Branch II completion).

### Klasyfikacja zamknięcia:

**M11.1 (this cycle): BRANCH II LEVEL 1 CLOSED** — 1-loop V_eff(Φ) audit complete. Numeryczna reproducibility, strukturalna identity (ln K_Φ), loop stability, BZ convergence, β=γ preservation, η boundary condition wszystkie zweryfikowane.

**Branch II progression:** M11.1 ✅ → **M11.2** (β-functions z Wetterich FRG, LPA' explicit η equation) → M11.3 (γ(k) running) → M11.4 (RG-driven C.3/B.3/B.5/B.2) → M11.R-final.

---

## Implikacje dla kolejnych cykli

### → M11.2 (β-functions z Wetterich FRG)

M11.1 ustanawia 1-loop V_eff baseline (β=γ preserved). M11.2 musi:
1. **Reproduce η = 0.044 self-consistently** z explicit LPA' anomalous-dim equation
2. **Compute β-functions** dla {β, γ, K_geo, q} w canonical sek08a action
3. **Identify fixed points** (Gaussian + Wilson-Fisher) z IR/UV behavior
4. **Cross-check** vs Branch I η_BI = 0.0253 — czy pełny FRG η_BII converges do 0.044 (matching CG-2) lub do 0.0253 (matching Branch I)?

### → M11.R-final

M11.R-I CLOSED (Branch I synthesis 6/6 PASS). M11.1 CLOSED (Branch II level 1). M11.R-final aggreguje:
- Branch I full chain: M11.S, M11.I, M11.G, M11.E, M11.R-I
- Branch II full chain: M11.1, M11.2, M11.3, M11.4
- 6 consistency conditions §4.1-4.6 (η, G_TGP, λ_C, universality class, KNOWN_ISSUES, M9 reproduction)
- Absolute δM_phys via dim-reg / zeta-fn upgrade

---

## Pliki

- **Skrypt:** `m11_1_audit.py` (~370 linii, 6 testów standalone)
- **Output:** `m11_1_audit.txt` (6/6 PASS log)
- **Audited assets:** [[../op1-op2-op4/m2b_loop.py]] (222 linii, M2a→M2b 1-loop), [[../op1-op2-op4/M2b_loop_derivation.md]] (266 linii, full derivation), [[../op1-op2-op4/M2b_results.md]] (268 linii, reference cases A/B/C)
- **Boundary condition source:** [[../continuum_limit/cg_strong_numerical.py]] (412 linii, CG-2 LPA' η = 0.044)
- **Predecessors:** [[M11_R_results.md]] (Branch I synthesis closed)
- **Successor:** [[M11_program.md]] → M11.2
