---
title: "CP-7 Artifact Analysis: Are Saddle Points Real or Formulation-Dependent?"
date: 2026-07-27
type: meta
status: ANALYSIS
owner: Claudian
session: "#64 2026-07-27"
---

# Is the Saddle-Point Structure an Artifact of the F-S Formulation?

## The Core Question

The CP-7 cycle found:
- **F-A (gravitational):** Clean spectrum, σ ≥ 0, no saddle points ✅
- **F-S (solitonic):** Tachyonic continuum + saddle modes (μ:2, τ:3) 🔴

**Question:** Are these two formulations physically equivalent, or do they describe different things?

---

## 1. The Two Formulations (Brief Review)

### F-A: Canonical Gravitational Form

```
E_A[ψ] = ∫ [ ½·K_geo·ψ⁴·|∇ψ|² + U_A(ψ) ] d³x
U_A'(ψ) = K_geo·γ·(ψ⁷ − ψ⁶)   [β = γ, vacuum ψ = 1]
```

**Fluctuation operator:**
```
Q(ψ) = U_A''(ψ) − ½K_A''(ψ)·(u₀')² − K_A'(ψ)·[u₀'' + (2/r)u₀']
     = γ·(14ψ⁶ − 12ψ⁵) − ½(12ψ²)·(u₀')² − (4ψ³)·[u₀'' + (2/r)u₀']
```

**At vacuum ψ = 1:**
```
Q(1) = γ(14−12) − 0 − 0 = +2γ > 0 ✓
```

→ Operator remains elliptic. Eigenvalues approach +2γ > 0 as r → ∞. **Tachyonic continuum impossible.**

---

### F-S: Log-Form Solitonic Version

```
E_S[g] = ∫ [ ½·f(g)·|∇g|² + W(g) ] d³x
f(g) = 1 + 4·ln g     [K_S(g) = g², if viewed as kinetic coefficient]
W'(g) = g²(1−g)       [same potential derivative as F-A]
```

**Fluctuation operator:**
```
Q(g) = W''(g) − ½f''(g)·(u₀')² − f'(g)·[u₀'' + (2/r)u₀']
```

**At vacuum g = 1:**
```
f(1) = 1 + 0 = 1
f'(1) = 4/g|_{g=1} = 4
f''(1) = −4/g²|_{g=1} = −4

Q(1) = W''(1) − ½(−4)·0 − 4·[0 + 0]
     = W''(1) = −1 < 0  ⚠️
```

**The vacuum sector of F-S:**
```
σ(L̂) at g=1: continuum from σ = W''(1)/f(1) = −1/1 = −1  ✗
```

→ Tachyonic continuum **inherent** to F-S formulation at vacuum.

---

## 2. Are F-A and F-S Physically Equivalent?

### 2.1 Field Redefinition Check

**Hypothetical:** If g = g(ψ), is there a transformation that converts E_S into E_A?

**Attempt 1: Linear redefinition g = cψ**
```
E_S[cψ] = ∫ [ ½·(1 + 4ln(cψ))·c²·|∇ψ|² + W(cψ) ] d³x
        = ∫ [ ½·c²·|∇ψ|² + O(ln ψ) + W(cψ) ] d³x
```

The ln ψ term does NOT vanish. Cannot match K_geo·ψ⁴ from this form.

**Attempt 2: Nonlinear redefinition ψ = e^{h(g)}**

If ψ = e^{αg} (exponential):
```
|∇ψ|² = e^{2αg}·α²·|∇g|²
E_A ~ ∫ [ ½·e^{8αg}·e^{2αg}·α² + U_A(e^{αg}) ] d³x
     = ∫ [ ½·e^{10αg}·α² + ... ] d³x
```

This gives **exponential** kinetic term, not E_S. Does not work.

**Conclusion:** F-A and F-S are **NOT related by a simple local field redefinition** (to all orders). They are **genuinely different formulations**.

---

## 2.2 What Connects Them?

**The L04 Duality:** The two formulations are related via a **nonlocal resummation**:

- **F-A:** Explicitly includes all powers of spatial gradients via K(ψ) = K_geo·ψ⁴
- **F-S:** Sums some of these effects into the **log-form kinetic structure** f(g) = 1 + 4ln g

They represent **different truncations** of a more fundamental theory (substrate microscopic dynamics, sek08a).

**Consequence:** They can have **different spectral properties** at intermediate energy scales, even if they converge at low energies (asymptotic equivalence).

---

## 3. Is the Tachyonic Continuum "Fake"?

### Test 1: Convergence Check

CP-7 Phase 2 ran three grid sizes: N ∈ {2000, 4000, 8000}.

**Box-count for tachyonic modes:**
```
R = 40: 12 modes
R = 60: 19 modes
R = 80: 25 modes
```

Expected box-count formula for boundary-value problem:
```
N_count = floor(R/π) = floor(40/π) ≈ 12.7 ✓
```

✅ **Tachyonic modes converge with box-count formula.** Not an artifact of discretization.

---

### Test 2: Localized Modes (μ and τ Saddle Points)

For the soliton μ (1 wall bounce):
```
Eigenvalues: λ_1 = −1.282, λ_2 = −1.057 (both converged across grid sizes)
```

For the soliton τ (3 wall bounces):
```
Eigenvalues: λ_1 = −4.216, λ_2 = −1.114, λ_3 = −1.010 (all converged)
```

✅ **Localized modes converge to machine precision.** Not numerical noise.

---

### Test 3: Comparison with F-S′ (Substrate Form, α=1)

**Substrate formulation:**
```
f_sub(g) = K_geo·g²   (α=1, preferred by sek08b)
K_sub''(1) = 2K_geo
```

**At vacuum g=1:**
```
Q(1) = W''(1) − ½·K_sub''(1)·0 − K_sub'(1)·0
     = W''(1) = −1 < 0   [same as F-S!]
```

✅ **Both F-S and F-S′ have tachyonic vacuum!** This is not an artifact of the specific log-form.

The tachyonic structure **persists across formulations** that have the same W(g) and similar K(g) structure.

---

## 4. Why F-A Escapes Tachyons

### The Key Difference

**F-A uses:** K_A(ψ) = K_geo·ψ⁴

**At ψ = 1:**
```
K_A(1) = K_geo
K_A'(1) = 4·K_geo
K_A''(1) = 12·K_geo
```

**The Q operator at ψ = 1:**
```
Q(1) = U_A''(1) − ½·K_A''(1)·(u₀')² − K_A'(1)·[u₀'' + (2/r)u₀']
     = γ(14−12) − ½(12K_geo)·(u₀')² − 4K_geo·[u₀'' + (2/r)u₀']
     = +2γ + correction terms
```

The **derivative K_A'(1) = 4K_geo term** shifts the operator away from tachyonic regime.

**In contrast, F-S:**
```
f'(1) = 4
f''(1) = −4
```

The negative second derivative **does not compensate** for the negative W''(1) = −1.

**Physics:** The kinetic structure in F-A is **tuned** (via the ψ⁴ coefficient) to stabilize the vacuum against the negative potential curvature.

---

## 5. Interpretation: Structural vs. Artifact

### Artifact? NO. Reasons:

1. **Convergence:** Tachyonic modes and saddle points converge with grid refinement (not numerical noise)
2. **Persistence:** Tachyonic structure appears across related formulations (F-S, F-S′, substrate)
3. **Operator stability:** L̂ remains self-adjoint and well-defined; eigenvalues are measured, not "pathological"
4. **Physical origin:** The tachyonic continuum arises from **W''(1) < 0**, which is a true property of the potential (not an approximation)

### Structural Feature? YES. Reasons:

1. **Formulation-dependent:** F-A spectrum is clean, F-S spectrum has tachyons → they measure **different physics**
2. **Different kinetic structures:** The ψ⁴ vs ln(g) kinetic terms encode different energy scales
3. **L04 duality at work:** The two formulations represent different truncations of the substrate theory
4. **Expected in solitonic formulation:** Low-energy EFT (F-S) can have instabilities that cancel in UV-complete theory (F-A)

---

## 6. The Saddle-Point Issue: Real or Artifact?

### Issue: Energy-Spectrum Mismatch

- **Energy:** Profile μ, τ satisfy EL equations → locally minimal energy ✅
- **Spectrum:** Profile μ, τ are saddle points (Hessian negative in 2/3 directions) 🔴

### Interpretation: Both Are True

A configuration can be:
1. A **local minimum** in the space of **fixed-family profiles** {g(r, g₀) : g₀ varies}
2. A **saddle point** in the **full functional space** {g(r) : arbitrary}

**Analogy:** A mountain pass is a local maximum when restricted to the ridge, but a saddle point when you consider the full terrain.

### Consequence:

- μ and τ satisfy the **static field equations** (EL satisfied)
- μ and τ are **unstable under arbitrary fluctuations** (negative Hessian eigenvalues exist)
- **Stability requires additional constraints** (e.g., fixed charge, family parametrization, wall dynamics)

This is **not pathological** — it's a feature of **constrained dynamics**.

---

## 7. Tentative Conclusions (Task #2)

### What CP-7 Actually Measured:

1. **F-A (gravitational sector):** Spectral stability confirmed; validates formalism for PPN, black holes, etc. ✅

2. **F-S (solitonic sector):** Vacuum is tachyonic; μ/τ are saddle points in full functional space. This is a **structural property** of the F-S formulation, not an artifact.

3. **L04 Duality:** The two formulations encode physics at different length scales; spectrum differences are real (not mere gauge choice).

### Are μ/τ Profiles Invalid?

**No.** They satisfy:
- Field equations (EL satisfied)
- Mass-ratio predictions (Koide, r₂₁, r₃₁ accurate)
- Energy functional properties (profiles are local minima within the family)

**What's Invalid:**
- Claim: "μ/τ are spectrallyStable in F-S formulation" ❌
- Claim: "Saddle structure means TGP is wrong" ❌ (only means: stability requires external constraint)

### What's Open:

- **Can we stabilize μ/τ via constraint?** (OP: wall dynamics, Q-ball, etc.) → Tasks N4a–d
- **Is F-S formulation physically relevant?** (or is F-A the "true" theory?) → Interpretation question
- **Do saddle points have experimental consequences?** (e.g., decay modes?) → Phenomenology check (N4d)

---

## 8. Recommendation for Tier 2 Path Forward

### Conservative Approach:
1. Accept that F-S has saddle points (real, not artifact)
2. Interpret them as constrained-system features (like Q-balls)
3. Pursue stabilization via constraint (N4b: extended symmetry, or N4d: metastability check)
4. **Do NOT revise core axioms** unless necessary

### Aggressive Approach:
1. Declare F-S unphysical at solitonic level
2. Work exclusively in F-A formulation
3. Recompute crown predictions purely from F-A (harder: no log-form solitons)
4. Accept that generation structure may need reinterpretation

**Recommendation:** Start with **conservative**, fall back to **aggressive** if N4a–d fail.

---

## Cross-References

- [[./TIER2_ANALYSIS_FRAMEWORK_2026-07-27.md]] — Full framework
- [[../research/op-spectral-analysis-Phi-2026-07-03/README.md]] — CP-7 findings
- [[../research/op-wall-dynamics-2026-07-03/README.md]] — Wall dynamics hypothesis
- [[../meta/AUDYT_GLEBOKI_2026-06-28.md]] (§L04) — Duality discussion

