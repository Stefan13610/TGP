---
title: "TIER 2 Analysis Framework — Spectral Saddle-Point Problem"
date: 2026-07-27
type: meta
status: ACTIVE
owner: Claudian
purpose: "Systematic organization of the CP-7 spectral instability (tachyonic continuum + saddle modes μ/τ) and evaluation of four forward paths (N4a–d)"
session: "#64 2026-07-27"
---

# TIER 2 Analysis Framework — TGP Spectral Instability

> **Problem statement:** CP-7 numerical diagonalization (2026-07-03) revealed that the F-S solitonic sector has a tachyonic continuum (σ ∈ [−1, 0)) and saddle-point modes for μ (l=0: 2 modes) and τ (l=0: 3 modes). This contradicts the spectral-synthesis theorem (L03, 2026-05-06) which claimed σ ⊂ [0,∞) for "all physical backgrounds". Three stabilization paths (#62, #63) failed. **Question: Is this a model artifact or genuine physics?**

---

## 1. THE CORE PROBLEM: Understanding the Saddle Points

### 1.1 What CP-7 Found

| Property | Result | Implication |
|----------|--------|-------------|
| **Gravitational sector (F-A)** | σ ≥ 0, clean | ✅ OK — stable |
| **Solitonic vacuum (F-S at g=1)** | Tachyonic continuum from σ = −1 | ⚠️ Instability |
| **Electron (gen-1, e)** | 0 localized modes | ✅ Stable |
| **Muon (gen-2, μ)** | 2 saddle modes: λ ∈ {−1.282, −1.057} | 🔴 Saddle point (index=2) |
| **Tau (gen-3, τ)** | 3 saddle modes: λ ∈ {−4.216, −1.114, −1.010} | 🔴 Saddle point (index=3) |

**Interpretation:** μ and τ are saddle points (local minima in 2/3 directions, but unstable in 2/3 other directions).

### 1.2 Energy vs. Spectrum Contradiction

**The paradox:**
- **Energy functional E_S:** μ and τ profiles have minimal energy **locally** (they satisfy Euler–Lagrange equations)
- **Spectral analysis:** The same profiles are saddle points (Hessian has negative eigenvalues)

This is not impossible in principle (saddle points can have local minima in restricted subspaces), but it's **unusual and requires interpretation**.

### 1.3 Why the Synthesis Theorem Failed

The spectral-synthesis theorem from 2026-05-06 claimed:
> "σ(L̂) ⊂ [0,∞) for all physical backgrounds" 

**Diagnosis:**
- This is **true in the F-A gravitational formulation** (proven numerically, 10 PPN parameters exact)
- This is **false in the F-S solitonic formulation**

**Root cause:** The synthesis derivation assumed Q(φ) → +γ asymptotically. This is a property of **F-A only**. In F-S, the boundary condition yields Q(φ) → −1, which generates the tachyonic continuum.

**This is the L04 duality measured at the source.**

---

## 2. WHAT IS NOT RULED OUT

### 2.1 ✅ Profile Properties Remain Valid

- Lepton mass ratios (r₂₁ = μ/e, r₃₁ = τ/e) derived from **profile shapes**, not spectral stability
- Koide formula K = 2/3 from **algebraic identity** (Brannen parametrization)
- Crown selection mechanism (why N=3 generations) still functions

### 2.2 ✅ Why Saddle Points Don't Invalidate Masses

Saddle points are still **solutions to the field equations**. They satisfy:
```
δE/δφ = 0  ⟺  EOM satisfied
```

The profile μ(r) is a valid classical configuration, regardless of spectral stability. The mass ratios come from the **shape** of the profile, not from stability.

**Analogy:** A particle sitting on a saddle point of a mountain has definite position; the fact that it's unstable doesn't change its height.

---

## 3. FAILED STABILIZATION PATHS (Sessions #62, #63)

### 3.1 Session #62: `op-wall-dynamics` — Linear Constraints (FAILED)

**Hypothesis:** If we impose a linear constraint on fluctuations (e.g., ∫δΦ d³x = 0), can we remove the saddle-point modes?

**Method:**
- K1–K4: Four types of constraints (simple integral, metric-weighted, kinetic-weighted, core-budget)
- K_i ∧ K_j: Pairs of constraints
- Criterion: Project operator L̂ onto null space of constraint

**Results (W1):**
- Single constraints (K1–K3) remove only edge modes (−1.0098 for μ), leave deep modes (−1.282, −4.216) untouched
- K4 (core budget) removes deep modes but leaves residual mode at λ ≈ 0
- Pairs K_i ∧ K_j: μ index → 1, τ index → 1, but residual mode is not orthogonal to ∂g/∂g₀ (which defines the family of profiles)

**Verdict:** Linear constraints alone cannot stabilize μ/τ within the family of profiles. A different approach needed.

### 3.2 Session #62: Soft Wall (f_ε) — Model Sensitivity (FAILED)

**Hypothesis:** The ghost wall (reflecting boundary at φ → 0) might be an artifact of the hard cutoff. Can a smooth approximation f_ε = ½[f + √(f² + ε²)] stabilize μ/τ?

**Method:**
- Test ε ∈ {0.2, 0.1, 0.05, 0.02}
- Check whether λ_min(ε → 0) converges
- Measure mass drift r₂₁(ε) and r₃₁(ε)

**Results (W2):**
- **τ collapses for every ε** (profile terminates at r ≈ 2.6, like α=1 substrate)
- τ exists **only with hard-wall reflection** (non-EL regularization)
- λ_min(ε → 0) **does not converge** for τ
- Mass drift r₂₁: +1.9% to +23% depending on ε (not <0.1% as required)

**Verdict:** The wall is not a soft-cutoff artifact. It's a **structural feature** of the model. Smooth regularizations destroy μ/τ (especially τ).

### 3.3 Session #63: Charge Extension (U(1)) — Nonlinear Test (FAILED)

**Hypothesis:** Extend Z₂ → U(1)×Z₂, formulate Q-ball dynamics with Vakhitov–Kolokolov criterion. Does fixed-charge quantization stabilize μ/τ?

**Method:**
- Formulate canonical Q-ball EOM: L_S = ½f ġ² − ½f|∇g|² − W(g)
- Quantize conserved charge Q = ω∫fφ² r² dr (Noether current)
- Q-ball family φ_ω(r): continuous from CP-7 profiles as ω → 0
- Vakhitov–Kolokolov criterion: dQ_sol/dω < 0 ⟹ stabilization

**Results (V1–V3):**
- **V1 NEGATIVE (charge inventory):** No conserved charge C1–C5 survives Noether test (9/9 candidates fail, exact in sympy)
- **V2 NEGATIVE (VK criterion):** dQ_sol/dω > 0 everywhere (slope positive; only e@ω=0.05 is slope-negative, outside hypothesis)
- **V3 NEGATIVE (nonlinear dynamics):** μ field runs away to g=0 in time t* ≈ 3.6 (nonlinearity does not stabilize)

**Verdict:** The U(1) extension does not provide a stabilizing charge. The three smooth stabilization paths are now **all exhausted**.

---

## 4. FOUR REMAINING PATHS (N4a–d)

Each path attacks the problem from a different angle:

### 4.1 **N4a: Discrete Substrate Hypothesis**

**Question:** Does shifting from continuous ŝ(r) to lattice ŝ_i break the tachyonic mode structure?

**Rationale:**
- CP-7 assumes continuum field theory (BVP on R)
- If the substrate is discrete (lattice), continuum modes might not exist
- Lattice dispersion ω(k) might cut off or modify the tachyonic band

**Plan:**
1. Reformulate substrate Hamiltonian on lattice (e.g., 1D ring with periodic BCs)
2. Diagonalize nearest-neighbor Heisenberg coupling
3. Compare spectrum: does tachyonic continuum disappear?
4. Can we recover μ/τ with index=0?

**Risk:** Entire TGP ontology built on continuous substrate. This is a deep axiom change.

**Criterion:** If discrete substrate removes tachyonic continuum AND preserves mass ratios (r₂₁, r₃₁ within 0.1%), then N4a is viable.

---

### 4.2 **N4b: Alternative Symmetry (Z₂ → larger group)**

**Question:** Is Z₂ the right symmetry? Could extending to U(1)×Z₂, SU(2), or discrete non-abelian group stabilize μ/τ?

**Rationale:**
- #63 tested naive U(1) charge conservation; it failed
- But SU(2) or other discrete symmetries might have different consequences
- Lepton-sector structure (e, μ, τ = 3 generations) suggests a 3-fold structure

**Plan:**
1. Survey discrete/continuous symmetry extensions that preserve Z₂ as subgroup
2. For each, derive Noether charges (multiple currents, not just one)
3. Test whether fixed multi-charge quantization (e.g., Q₁, Q₂ simultaneously) stabilizes
4. Check whether mass ratios persist under symmetry change

**Risk:** Introducing new symmetries changes the entire gauge sector. Could conflict with phenomenology (no new massless bosons observed).

**Criterion:** If we find a symmetry with N ≥ 2 independent conserved charges, both stably quantizable, and E_min at (Q₁, Q₂) with index=0, then N4b is viable.

---

### 4.3 **N4c: Disformal Sector Coupling (F-A Radiative Corrections)**

**Question:** Does including radiative corrections from the gravitational sector (F-A) modify the effective potential for the soliton and remove the saddle points?

**Rationale:**
- CP-7 diagonalized L̂ for the **solitonic sector in isolation** (F-S)
- But the full theory couples F-S to F-A (metric-dependent couplings)
- Radiative loops in the gravitational sector might renormalize V_eff(φ)

**Plan:**
1. Compute one-loop effective potential V_eff^(1)(φ) including F-A box diagrams
2. Solve for renormalized vacuum φ_0(1-loop)
3. Diagonalize L̂ with V_eff^(1) instead of bare V
4. Check: do saddle-point modes vanish?

**Risk:** One-loop computations are technically demanding and scheme-dependent. Result might depend on renormalization scale.

**Criterion:** If V_eff^(1)(φ) has a simple minimum (V″ > 0) and gives stable spectrum, then N4c is viable.

---

### 4.4 **N4d: Metastability Hypothesis (Tunneling)**

**Question:** If saddle points are real, are μ/τ metastable? Can they decay via WKB tunneling to other minima?

**Rationale:**
- Saddle points can be long-lived if tunneling rates are suppressed
- If tunneling lifetime >> age of universe, configuration is effectively stable
- This is physically acceptable (e.g., proton decay in GUT models)

**Plan:**
1. Map the energy landscape around μ and τ profiles
2. Identify potential minima (other than the current profile)
3. Calculate tunneling rate T via WKB: Γ ~ e^{−S_bounce}
4. Estimate lifetime τ_decay = 1/Γ
5. Compare to 13.8 Gyr (age of universe)

**Risk:** WKB tunneling calculations are sensitive to the exact form of V(φ). If tunneling is fast, μ/τ decay, and TGP is ruled out experimentally (would show τ → other particles).

**Criterion:** If tunneling lifetime > 10^{18} s, then metastability is consistent with observations.

---

## 5. EVALUATION MATRIX

| Path | Mechanism | Difficulty | Physics Cost | Prediction Impact | Timeline |
|------|-----------|------------|--------------|-------------------|----------|
| **N4a** | Discrete substrate | Very High | Deep axiom change | Large (recalculate all masses) | 2 weeks |
| **N4b** | Extended symmetry | High | New gauge symmetry | Medium (mass ratios preserved?) | 1.5 weeks |
| **N4c** | Radiative corrections | High | Perturbative; scheme-dependent | Small (mass ratios preserved) | 1 week |
| **N4d** | Metastability | Medium | None (saddle points OK) | None (phenomenology unchanged) | 3 days |

---

## 6. DECISION TREE

```
Are the saddle points real physics?
│
├─ YES (they persist in all tests)
│   ├─ Can we live with them? (Metastability test)
│   │   ├─ τ_decay >> age of universe → ACCEPT (N4d)
│   │   └─ τ_decay < age of universe → FAIL (TGP ruled out)
│   │
│   └─ Must we fix them? (Axiom change)
│       ├─ Try discrete substrate (N4a)
│       ├─ Try extended symmetry (N4b)
│       └─ Try radiative corrections (N4c)
│
└─ NO (they vanish under axiom/model change)
    └─ Adopt new axiom, recalculate, check mass ratios
```

---

## 7. CURRENT STATUS

### ✅ Closed (Sessions #62–#63)
- Linear constraints insufficient (W1)
- Soft-wall regularization fails (W2)
- U(1) charge quantization fails (V1–V3)

### ⏳ Pending (Session #64 onwards)
- N4a evaluation (discrete substrate)
- N4b evaluation (alternative symmetry)
- N4c evaluation (radiative corrections)
- N4d evaluation (tunneling lifetime)

### 🎯 Deliverables This Session
1. Task #2 — Interpret whether saddle points are F-S artifacts or structural
2. Tasks #3–6 — Preliminary assessment of N4a–d feasibility
3. Task #7 — Decision on which path to pursue (or pivot)

---

## 8. WHAT SUCCESS LOOKS LIKE

**Success = One of:**
1. Saddle modes vanish under axiom/model change (N4a/b/c) → core revision + mass ratios recomputed + mass ratios preserved
2. Metastability proven (N4d) → tunneling lifetime >> universe age → TGP phenomenologically consistent despite saddle structure
3. Interpretation clarified → document why saddle structure is acceptable → update Limitations in lepton paper

**Failure = None of the above:**
→ Saddle-point instability fundamentally incompatible with TGP → model needs deeper revision

---

## Cross-References

- [[../research/op-spectral-analysis-Phi-2026-07-03/README.md]] — CP-7 findings
- [[../research/op-wall-dynamics-2026-07-03/Phase0_balance.md]] — Wall dynamics hypothesis
- [[../research/op-nonlinear-charge-constraint-2026-07-03/README.md]] — Q-ball test
- [[../audyt/L03_K_phi_stability/README.md]] — Ontological context
- [[../meta/AUDYT_GLEBOKI_2026-06-28.md]] — Tier 2 original scope

