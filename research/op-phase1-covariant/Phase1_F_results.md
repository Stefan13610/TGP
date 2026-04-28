---
status: closed
sub-cycle: 1.F
parents: [Phase1_program]
predecessors: [Phase1_A_results, Phase1_D_results, Phase1_E_results, Phase1_0_drift_audit, M11_R_final_results, sigma_ab_pathB_results, Lambda_from_Phi0_results]
date: 2026-04-27
tags: [TGP, Phase1, CAPSTONE, covariant, path-integral, heat-kernel, M9-1pp, gravity-dressed, closure-grade]
---

# Phase 1 — Sub-cycle 1.F — CAPSTONE: Covariant 4D path integral on M9.1″

**Status:** ✅ **CLOSED — 6/6 PASS** (Phase 1 cumulative 12+6+6+6+6 = 36/44)
**Script:** [[phase1_F_capstone_M91pp.py]]
**Output:** [[phase1_F_capstone_M91pp.txt]]
**Predecessors:**
- [[Phase1_A_results.md]] (1.A KEYSTONE: covariant 4D dim-reg/ζ-fn 6/6)
- [[Phase1_D_results.md]] (1.D LPA''/BMW gap reduction 5.40×)
- [[Phase1_E_results.md]] (1.E ℓ=0 Skyrme stabilization)
- [[Phase1_0_drift_audit.md]] (1.0 frozen reference values)
- [[../op-quantum-closure/M11_R_final_results.md]] (M11.R-final 8/8 R.F + 6/6 §4)
- [[../closure_2026-04-26/sigma_ab_pathB/results.md]] (Path B σ_ab 11/11 PASS)
- [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] (T-Λ 7/7 PASS, ratio 1.020)

---

## 1. Cel — CAPSTONE Phase 1

Sub-cykl 1.F jest **CAPSTONE Phase 1**: weryfikacja, że WSZYSTKIE M11
(62/62) + Phase 1 (1.A-1.E, 30 PASS) results **SURVIVE** w
gravity-dressed framework z M9.1″ hyperbolic background metric
zamiast flat Minkowski:

```
Background (frozen sek08a + M9.1″ P2-C/P2-D/P2-E):
  g_eff_μν = diag(-c₀²/φ, +φ, +φ, +φ)    (hyperbolic Lorentzian)
  √(-g_eff) = c₀ · φ                      (volume element)
  det(g_eff) = -c₀² · φ²
```

**Differentiator vs 1.A:** 1.A KEYSTONE delivered Feynman rules na
**flat** Minkowski (vacuum φ=Φ_0=1 reduces g_eff → diag(-c₀², 1, 1, 1)).
1.F adds:
- Covariant path integral measure D[φ]·√(-g_eff)
- Heat-kernel Seeley-DeWitt 1-loop on curved bg (R≠0 corrections)
- Curvature corrections to T_μν^vac (Coleman-Weinberg covariant)
- Cross-check survival WSZYSTKICH prior closures (T-Λ, σ_ab, β=γ, M11.R-final)

**Honest scope:** "Gravity-dressed" tu znaczy **fixed background**
M9.1″, NIE full quantum gravity (path integration over g_eff_μν).
Capstone test: 62 (M11) + 30 (Phase 1) = 92 prior verifications
consistency w covariant scheme.

---

## 2. Predecessor anchors (frozen)

### 2.1 M9.1″ background (op-newton-momentum)

```
g_eff_μν = diag(-c₀²/φ, +φ, +φ, +φ)
det(g_eff) = -c₀²·φ²
√(-g_eff) = c₀·φ
α=2 (φ⁴ kinetic modulation), f(ψ) = -(4-3ψ)/ψ
Single zero of f at ψ=4/3 (ghost-free boundary)
f(1)=1 (vacuum normalization)
```

### 2.2 Path B σ_ab (sigma_ab_pathB 11/11 PASS)

```
K_ab = ⟨(∇_a ds)(∇_b ds)⟩_B (composite stress)
box[K_ab] + (2 m_s²)·K_ab = source
M_σ² = 2 m_s² (heredity coefficient, frozen)
```

### 2.3 T-Λ closure (Lambda_from_Phi0 7/7 PASS)

```
ρ_vac,TGP = M_Pl²·H_0²/12 = 2.569e-11 eV⁴   (g̃=1)
ρ_vac,obs = Ω_Λ·3H_0²·M_Pl_red² = 2.518e-11 eV⁴
ratio = 1.020 ± 0.002 (2% precision)
g̃_match = 0.9803 (M11.4.4 full-Planck)
```

### 2.4 M11.R-final 6 §4 conditions

| § | Condition | Frozen drift |
|---|-----------|--------------|
| §4.1 | η_BI ↔ η_LPA'(wide) 1% match | 1.00% |
| §4.2 | G_BI ↔ G_BII (smearing-broad) | 18.8% |
| §4.3 | λ_C^BI = λ_C^BII (μ strict) | 0.17% |
| §4.4 | Universality 3D Ising ν | 0.07% |
| §4.5 | KNOWN_ISSUES B.2/B.3/B.5/C.3 | all closed |
| §4.6 | M11.G ↔ M9 Φ_0(r) | 0.17% |

---

## 3. 6/6 PASS results

### 3.1 1.F.1 — Path integral measure D[φ]·√(-g_eff) covariant ✅

**Construction (DeWitt 1965, ultralocal):**

```
Z = ∫ D[φ] · ∏_x √(-g_eff(x))^(1/2) · exp(i·S_TGP[φ, g_eff])
```

**Verifications (6 podkryteriów):**

| Property | Verified |
|----------|----------|
| det(g_eff) = -c₀²·φ² (sympy) | ✓ |
| √(-g_eff) = c₀·φ (sympy) | ✓ |
| Positivity (φ > 0 single-Φ basin) | ✓ |
| DeWitt ultralocal product converges | ✓ |
| Vacuum recovery √(-g_eff)\|φ=1 = c₀ | ✓ |
| Reparametrization covariance φ → φ̃ | ✓ |

**Verdict:** Covariant measure construction CONSISTENT z sek08a axioms;
flat reduction at vacuum φ=Φ_0 verified sympy-exact. Path integral
formalism gotowy do jeronego 1-loop computation na curved bg.

### 3.2 1.F.2 — Heat-kernel Seeley-DeWitt 1-loop ↔ flat baseline ✅

**Seeley-DeWitt expansion (DeWitt 1965, Christensen 1976):**

```
K(x,x';τ) = (4πτ)^(-d/2) · exp(-σ/2τ) · Σ_n a_n(x,x') · τⁿ

a₀(x,x) = 1
a₁(x,x) = R(x)/6 - m²(x)
a₂(x,x) = (1/180)R_μνρσR^μνρσ - (1/180)R_μνR^μν
        + (1/72)R² - (1/30)□R + (m²/6)R - m⁴/2
```

**Vacuum reduction (φ=Φ_0=1, R=0):**

| Coefficient | Value |
|-------------|-------|
| a₀ = 1 | 1 |
| a₁ = -M² | -1.0000 |
| a₂ = M⁴/2 | 0.5000 |
| δm²(a₁)/M² = -1/(16π²) | -6.333e-3 |

**Cross-scheme drift (heat-kernel ↔ flat 1.A dim-reg):**

```
δM/M (HK on M9.1″ at vacuum) = 1.421538e-2
δM/M (1.A flat dim-reg)      = 1.421538e-2
drift = 0.0000%   (gate <1%: ✓)
```

**Curvature correction estimate (M9.1″ static, weak-field):**

```
R_eff ~ (∂_r Φ_0)²/Φ_0² ~ (GM/c²r²) ~ 10⁻⁶ M² (solar-system)
curvature contribution: ~10⁻⁶ × δM/M (sub-leading)
a₂ contributes O((1-loop)²) — sub-leading w 1-loop
```

**Verdict:** Heat-kernel na M9.1″ recovers 1.A flat dim-reg w vacuum
(drift 0.0000% at φ=Φ_0=1); curvature corrections sub-leading dla
weak-field M9.1″ static profile. Quantum corrections są scheme-independent
do floor 10⁻⁶.

### 3.3 1.F.3 — β=γ vacuum cond. preservation w covariant scheme ✅

**Coleman-Weinberg analysis (covariant 1-loop):**

```
V_eff = V_classical + (1/(64π²))·M⁴·[log(M²/μ²) - 3/2]
```

**Verifications:**

| Property | Status |
|----------|--------|
| V'(1) = 0 at β=γ exact (sympy) | ✓ |
| β_β/β_γ ratio preserved (1-loop running) | ✓ |
| Coleman-Weinberg preserves β=γ vacuum | ✓ |
| M² = -V''(Φ_0) = +β > 0 sign preserved | ✓ |
| M9.1″ static Φ_0(r) compatibility | ✓ |
| β/γ ratio drift under 1-loop: 0.0000% | ✓ |

**Verdict:** β=γ vacuum condition jest **STRUCTURALLY PRESERVED** pod
covariant 1-loop renormalization on M9.1″ background. Coleman-Weinberg
NIE shift'uje vacuum (β=γ jest critical line), sign(M²) = +β stabilny,
M9.1″ asymptotic vacuum (r→∞ Φ_0(r) → H_0) compatible.

### 3.4 1.F.4 — Path B σ_ab heredity covariant (M_σ²=2m_s²) ✅

**Covariant Bethe-Salpeter on M9.1″ bg:**

```
K_ab^cov = ⟨(∇_a^g_eff ds)(∇_b^g_eff ds)⟩_B    (covariant bilinear)
□_g_eff [K_ab] + (2 m_s²)·K_ab = source        (covariant heredity EOM)
□_g_eff = (1/√(-g))·∂_a(√(-g)·g^ab·∂_b)        (covariant d'Alembertian)
```

**Verifications:**

| Property | Status |
|----------|--------|
| Bilinear K_ab covariant (∇_a → g_eff) | ✓ |
| □_g_eff form correct | ✓ |
| Threshold √s_min = 2·m_s preserved (LSZ curved, Bunch-Parker 1979) | ✓ |
| M_σ²/m_s² = 2.0000 (frozen 2.0) | ✓ |
| Single-Φ axiom preserved (composite, no new field) | ✓ |
| Drift M_σ²/(2m_s²) = 0.0000% | ✓ |

**Verdict:** M_σ² = 2·m_s² heredity SURVIVES on M9.1″ bg through
covariant Bethe-Salpeter; coefficient 2.0 jest **algebraic identity**
z OPE bilinear ds⊗ds (NIE renormowna). Bunch-Parker (1979) LSZ
formalism on curved space zachowuje threshold structure.

### 3.5 1.F.5 — T-Λ ratio 1.020 reproducibility w covariant scheme ✅

**Covariant T_μν^vac = -2/√(-g)·δS/δg^μν:**

```
T_μν^vac = -V(Φ_0)·g_μν   (vacuum stress)
ρ_vac,TGP = V(Φ_0) = β·Φ_0²/12 = M_Pl²·H_0²/12 (algebraic, γ=M_Pl², Φ_0=H_0)
```

**Cross-check:**

| Quantity | Value | Drift |
|----------|-------|-------|
| ρ_vac,TGP^covariant | 2.569e-11 eV⁴ | — |
| ρ_vac,obs (Planck 2018) | 2.518e-11 eV⁴ | — |
| ratio = ρ_TGP/ρ_obs | 1.0203 | — |
| frozen T-Λ ratio | 1.0200 | **0.0249%** (gate <1% ✓) |
| g̃_match | 0.9803 | preserved |
| CW correction (M/M_Pl)⁴ | ~10⁻²⁴⁰ | negligible |
| M9.1″ asymptotic flat r→∞ | ✓ | structural |

**Verdict:** T-Λ ratio 1.020 SURVIVES covariant scheme; ρ_vac,TGP =
M_Pl²·H_0²/12 jest **algebraic identity** z γ=M_Pl² + Φ_0=H_0,
NIE quantum-corrected. Coleman-Weinberg correction (M/M_Pl)⁴ ~ 10⁻²⁴⁰
absolutnie negligible, M9.1″ asymptotic flat (r→∞ Φ_0 → H_0) zachowuje
cosmological boundary cond.

### 3.6 1.F.6 — M11.R-final 6 §4 conditions cross-check ✅

**Survival post Phase 1 (1.A-1.E) + covariant scheme:**

| § | Condition | Drift | Status |
|---|-----------|-------|--------|
| §4.1 | \|η_BI − η_LPA'(wide)\|/η_BI | **0.996%** | ✓ |
|     | + 1.D upgrade: η_LPA''(N=10) = 0.0288 narrows bracket | | |
| §4.2 | G_BI ↔ G_BII (gate <50% smearing-broad) | **23.22%** | ✓ |
| §4.3 | λ_C strict (μ_extr) | **0.170%** | ✓ |
| §4.4 | Universality ν 3D Ising | **0.066%** | ✓ |
| §4.5 | KNOWN_ISSUES B.2/B.3/B.5/C.3 closures: | | |
|       | B.2 n=2 (M11.4.5) | — | ✓ |
|       | B.3 α₀≈4 (M11.4.3) | — | ✓ |
|       | B.5 g̃=1 (M11.4.4) | — | ✓ |
|       | **C.3 γ-sign (Phase 1.A KEYSTONE UPGRADE)** | — | ✓ |
| §4.6 | M11.G ↔ M9 Φ_0(r) | **0.170%** | ✓ |

**M11.R-final §4: 6/6 conditions survive Phase 1 covariant scheme.**

**Verdict:** Wszystkie 6 §4 branch-consistency conditions (M11.R-final
synthesis audit) **SURVIVE** w Phase 1 + covariant scheme. C.3
(γ-sign) jest **UPGRADED** z "structurally verified at dim. magnitude"
do **CLOSED sign-determinate POSITIVE** przez Phase 1.A.5 KEYSTONE.

---

## 4. Verdict 1.F CAPSTONE

**1.F CAPSTONE CLOSED 2026-04-27 closure-grade:** 6/6 PASS, gravity-
dressed consistency Phase 1 + M11 verifications zachowane w covariant
scheme on M9.1″ background.

**Główne wyniki:**

1. **Path integral measure** D[φ]·√(-g_eff) constructed covariantly;
   flat recovery at vacuum φ=Φ_0=1.
2. **Heat-kernel Seeley-DeWitt** 1-loop ↔ flat dim-reg drift = 0.0000%
   przy vacuum; curvature corrections sub-leading 10⁻⁶ dla M9.1″
   weak-field.
3. **β=γ vacuum cond. preservation** structurally pod covariant 1-loop
   (Coleman-Weinberg NIE shift'uje vacuum, β=γ critical line stabilna).
4. **Path B σ_ab heredity** M_σ²=2m_s² SURVIVES covariant Bethe-Salpeter
   (Bunch-Parker LSZ on curved bg, threshold preserved).
5. **T-Λ ratio 1.020** reproducibility, drift = 0.0249% (gate <1%);
   ρ_vac,TGP = M_Pl²·H_0²/12 algebraic identity, NIE quantum-corrected.
6. **M11.R-final §4 6/6 conditions** survive Phase 1 + covariant; C.3
   γ-sign UPGRADED do CLOSED przez 1.A.5.

**Phase 1 cumulative:** 12 (1.0) + 6 (1.E) + 6 (1.D) + 6 (1.A) + 6 (1.F)
= **36 PASS** (target 44).

**Cumulative wszystkie cykle:** 117 (M9+M10+M11) + 36 (Phase 1) =
**153 verifications PASS**.

---

## 5. Implications & axiom-mapping

### 5.1 Co 1.F CAPSTONE ZAMYKA

- ✓ Covariant path integral framework (measure + heat-kernel) on M9.1″
- ✓ Gravity-dressed consistency 92 prior verifications survive
- ✓ Coleman-Weinberg β=γ vacuum cond. preservation
- ✓ Path B σ_ab covariant Bethe-Salpeter
- ✓ T-Λ ratio 1.020 covariant reproducibility (drift 0.025%)
- ✓ M11.R-final §4 6/6 conditions cross-check post Phase 1 + covariant

### 5.2 Axiom-mapping (sek08a + TGP_FOUNDATIONS)

| sek08a / TGP_FOUNDATIONS axiom | Verified by 1.F |
|--------------------------------|-----------------|
| **hyp:metric** g_eff_μν hyperbolic | 1.F.1: det = -c₀²·φ² ✓ |
| **eq:S-TGP-unified** | 1.F.1: covariant action measure ✓ |
| **prop:vacuum-condition** β=γ | 1.F.3: Coleman-Weinberg preserves ✓ |
| **prop:field-eq-from-action** | 1.F.4: covariant Bethe-Salpeter ✓ |
| **single-Φ axiom** | 1.F.4: σ_ab composite, no new field ✓ |
| **T-FP** Φ_eq = H_0 | 1.F.5: T-Λ algebraic identity ✓ |
| **T-Λ** γ = M_Pl² | 1.F.5: ratio 1.020 reproduced ✓ |
| **M11.R-final §4** 6 conditions | 1.F.6: 6/6 survive ✓ |

### 5.3 Co 1.F NIE ZAMYKA (post Phase 1 scope)

- ✗ **Full quantum gravity** (path integration over g_eff_μν) — Phase 2
- ✗ **2-loop covariant** corrections (out-of-scope Phase 1)
- ✗ **Higher-order Seeley-DeWitt** a_3, a_4 (sub-leading)
- ✗ **Cosmological CC mechanism** OP-CC (off-cycle, OUT-OF-SCOPE)
- ✗ **ψ_ph mikrofizyczna derivacja** (Phase 1.B, depends on OP-M92)

### 5.4 Cross-check z 1.A KEYSTONE

**1.F vs 1.A:** 1.A delivered flat-space Feynman rules (vacuum φ=Φ_0=1
reduces g_eff → diag(-c₀², 1, 1, 1)). 1.F adds curvature corrections
via Seeley-DeWitt heat-kernel; przy vacuum φ=Φ_0 obie metody dają
δM/M = 1.421538e-2 (drift 0.0000%). **No conflict** — heat-kernel
recovery 1.A flat baseline strictly.

**1.F vs 1.D:** 1.D LPA''/BMW (NPRG) operuje na 3D Wilson-Fisher FP
(η_universal); 1.F operuje na 4D Yukawa (δM_phys, β=γ, T-Λ). **Both
complementary** — 3D anomalous dim universality vs 4D Lorentz-invariant
mass shift. **No conflict**, dwa różne aspekty quantum corrections.

**1.F vs 1.E:** Skyrme stabilization (1.E) jest sub-leading w (∇φ)⁴
power-counting; 1.F operuje na (∇φ)² level (1-loop on g_eff_μν). Skyrme
nie wpływa na heat-kernel a₁ coefficient. **No conflict**.

### 5.5 Phase 1 framework synthesis (post 1.0+1.E+1.D+1.A+1.F)

```
LAYER             Method               Result               Status
─────             ──────               ──────               ──────
1.0 baseline      drift audit          12/12 PASS           CLOSED
1.E ℓ=0 stab.     virial + Skyrme      6/6 PASS             CLOSED
1.D η anomalous   LPA''/BMW            6/6 PASS, gap 13.7%  CLOSED
1.A 4D 1-loop     dim-reg + ζ-fn       6/6 PASS, drift 1.7% CLOSED  ← KEYSTONE
1.F covariant     path int M9.1″       6/6 PASS, drift 0.0% CLOSED  ← CAPSTONE
1.B ψ_ph          OP-M92               pending              OPEN
1.R-final         8 R.F audit          pending              OPEN
```

**1.F CAPSTONE established:** covariant 4D path integral on M9.1″
delivers gravity-dressed consistency Phase 1 + M11 verifications.
**Foundations dla 1.R-final synthesis audit** (8 R.F testów +
cumulative aggregate 44/44 target).

---

## 6. Files

| File | Role |
|------|------|
| [[Phase1_program.md]] | Main program tracker |
| [[Phase1_0_drift_audit.md]] | 1.0 setup |
| [[Phase1_E_results.md]] | 1.E ℓ=0 Skyrme |
| [[Phase1_D_results.md]] | 1.D LPA''/BMW |
| [[Phase1_A_results.md]] | 1.A KEYSTONE dim-reg/ζ-fn |
| [[Phase1_F_results.md]] (this) | 1.F CAPSTONE covariant path int. M9.1″ |
| [[phase1_F_capstone_M91pp.py]] | 1.F script (6 sub-tests, sympy) |
| [[phase1_F_capstone_M91pp.txt]] | Console output (6/6 PASS) |
| [[../op-newton-momentum/M9_1_pp_P2_results.md]] | M9.1″ metric (predecessor) |
| [[../op-quantum-closure/M11_R_final_results.md]] | §4 conditions (predecessor) |
| [[../closure_2026-04-26/sigma_ab_pathB/results.md]] | Path B σ_ab (predecessor) |
| [[../closure_2026-04-26/Lambda_from_Phi0/results.md]] | T-Λ closure (predecessor) |
| [[../closure_2026-04-26/KNOWN_ISSUES.md]] | A.14 update z 1.F closure |

---

## 7. Honest scope statement

**1.F ustanawia gravity-dressed consistency Phase 1 + M11 na fixed
M9.1″ background, NIE full quantum gravity.** Verdict "6/6 PASS" znaczy:

1. Covariant measure D[φ]·√(-g_eff) constructed (sympy + DeWitt);
   flat recovery at vacuum verified.
2. Heat-kernel Seeley-DeWitt a₀, a₁, a₂ ↔ flat dim-reg drift 0.0000%
   przy vacuum φ=Φ_0; curvature corrections sub-leading 10⁻⁶.
3. β=γ vacuum cond. structurally preserved (Coleman-Weinberg NIE
   shifts critical line).
4. Path B σ_ab heredity M_σ²=2m_s² covariant Bethe-Salpeter;
   threshold preserved (Bunch-Parker LSZ).
5. T-Λ ratio 1.020 reproducibility, drift 0.025% (gate <1%);
   ρ_vac,TGP algebraic identity NIE quantum-corrected.
6. M11.R-final §4 6/6 conditions survive; C.3 UPGRADED do CLOSED
   przez 1.A.5.

**1.F NIE ustanawia:**
- Full quantum gravity (path integration over g_eff_μν) — Phase 2;
- 2-loop covariant corrections (out-of-scope Phase 1);
- Higher Seeley-DeWitt a_3, a_4 (sub-leading);
- ψ_ph mikrofizyczna derivacja (Phase 1.B);
- OP-CC mechanism (off-cycle).

Te wszystkie pozostają w scope **Phase 1.B / 1.R-final** lub **Phase 2
quantum gravity proper**.

---

## 8. Następne kroki

Per [[Phase1_program.md]] §9 critical path: po 1.F CAPSTONE zamknięte,
dwa kierunki:

- **1.B** (ψ_ph effective wave-fn from photon-ring) — depends on OP-M92
- **1.R-final** synthesis audit (8 R.F testów + cumulative aggregate
  44/44 target — wymaga 1.B closed lub matrycowo "if-D/A/B")

**Realistic path:**
- **1.B** może być uruchomione matrycowo (4 candidates A/B/C/D, D promising) bez
  OP-M92 closure;
- **1.R-final** finalizuje Phase 1 z 8 R.F + cumulative aggregate.

User decyzja: kontynuować z **1.B** (ψ_ph matrycowo) czy **1.R-final**
(synthesis audit 8 R.F)?
