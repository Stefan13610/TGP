---
title: "N4 Paths Feasibility Assessment — Four Routes Forward from CP-7 Saddle Points"
date: 2026-07-27
type: meta
status: FEASIBILITY_ANALYSIS
owner: Claudian
session: "#64 2026-07-27"
---

# N4a–d: Preliminary Feasibility Assessment

Each path addresses the CP-7 saddle-point problem from a different angle. Below is a systematic evaluation.

---

## N4a: Discrete Substrate Hypothesis

**Question:** Does a lattice substrate (instead of continuous field) eliminate tachyonic modes?

### Mechanism
```
Continuous substrate (current):
  H_Γ = Σ_i J_ij · ŝ_i ŝ_j  (GL Ginzburg–Landau bond, v2 2026-04-24)
  → Coarse-grain via renormalization group (RG)
  → Continuum field ŝ(r) with kinetic term K(φ)·|∇φ|²

Discrete substrate (proposed):
  H_lattice = Σ_i J_ij · ŝ_i ŝ_j  (same bond, but keep lattice index i, j ∈ 1..N_sites)
  → Diagonalize on finite lattice (e.g., 1D ring, 2D square, 3D cubic)
  → Spectrum ω(k) with k = 2πn/N, cutoff at k_max = π/a (Brillouin zone)
  → Check: do tachyonic modes disappear at continuum limit N → ∞?
```

### Pros
- ✅ **Preserves Z₂ axiom** — substrate is still ŝ with Z₂ symmetry
- ✅ **Keeps single-Φ ontology** — no new fields added
- ✅ **Existing substrate code** in `zusatz008c_metryka_z_substratu` could be adapted
- ✅ **Testable in principle** — compute spectrum for increasing N, check convergence

### Cons
- ❌ **Deep axiom change** — TGP built on continuum field theory (PDE, variational calculus)
- ❌ **Mass predictions recalculate** — lepton/quark masses depend on continuum soliton profiles; lattice profiles will differ
- ❌ **Computational explosion** — 3D lattice diagonalization: 10^3–10^6 sites needed for realistic Brillouin zone
- ❌ **Lose continuous translation** — discrete lattice breaks continuous SO(3) → discrete octahedral O_h symmetry; affects PPN, GW predictions
- ❌ **Risk of catastrophic change** — if mass ratios (r₂₁, r₃₁, Koide) shift by >0.1%, the entire model fails

### Technical Requirements
1. Substrate Hamiltonian on lattice
2. Python: sparse matrix diagonalization (scipy.sparse.eigh on 10^5–10^6 matrix)
3. Scan N_sites, extract spectrum, measure tachyonic band vs continuum limit
4. Recalculate soliton profiles on lattice (replacing BVP with eigenvector extraction)
5. Recompute mass ratios via new profiles

### Timeline
- **Research phase:** 2 weeks (theory + code setup)
- **Computation:** 1 week (parameter sweeps, convergence tests)
- **Result interpretation:** 3 days
- **Total:** ~3 weeks (if no major failures)

### Decision Criterion
**VIABILITY = YES if:**
- Tachyonic continuum vanishes by N → ∞
- **AND** mass ratios r₂₁, r₃₁ remain within 0.1% of CP-7 baseline (206.77, 3476.9)
- **AND** core .tex modifications are additive only (no removal of continuum formalism)

**Else:** Path N4a fails → pursue N4b/c/d

### Risk Level: 🔴 VERY HIGH

---

## N4b: Alternative Symmetry (Z₂ → Extended Group)

**Question:** Does extending the symmetry group beyond Z₂ provide new conserved charges for stabilization?

### Mechanism
```
Current:  Z₂ ≡ φ → −φ   (discrete, abelian)
          Noether charge: Q_Z₂ = 0 (Z₂ is discrete; no continuous parameter ω)

Extended: Z₂ × U(1)   → (φ, e^{iθ}·φ)
          OR Z₂ × Z₃  → (φ, e^{2πi/3·k}·φ)  [for 3 generations]
          OR SU(2)    (non-abelian, 3 generators)
          
Noether charges:
  Q_U(1) = ω ∫ |∂_tφ|² d³x
  Q_a (SU(2)) = ∫ φ^a d³x  (for a=1,2,3)
  
Fixed-charge quantization: E(Q₁, Q₂) minimum at physical charges → stable?
```

### Pros
- ✅ **Single-Φ preserved** — no new fundamental fields
- ✅ **Lepton structure suggestive** — three generations → 3-fold symmetry hint
- ✅ **Q-ball formalism applies directly** — extend #63 U(1) test to multi-charge
- ✅ **Moderate axiom change** — symmetry enlargement within Z₂×? structure
- ✅ **Existing physics analogs** — QCD has SU(3), electroweak has SU(2)×U(1)

### Cons
- ❌ **Gauge sector complications** — new symmetry requires new gauge bosons (massless unless Higgsed)
- ❌ **Why Z₂ × U(1) specifically?** — choice is somewhat arbitrary; many alternatives exist
- ❌ **Multi-charge quantization is hard** — how do ω₁, ω₂ couple? Are they independent or correlated?
- ❌ **Phenomenology risk** — new massless/light bosons (e.g., extra photons) would be observed; tension with experiments
- ❌ **Mass predictions sensitive** — extending symmetry changes interaction vertices; r₂₁, r₃₁ could shift significantly

### Technical Requirements
1. Extend Lagrangian to include new gauge sector
2. Derive multi-current Noether charges (e.g., Q₁ from U(1), Q₂ from another conserved current)
3. Solve Q-ball equations with dual-charge constraint: E = E(Q₁, Q₂)
4. Test Vakhitov–Kolokolov criterion for both charges: dQ_i/dω_i < 0
5. Check spectral stability of fixed-(Q₁, Q₂) family

### Timeline
- **Literature survey + design:** 1 week
- **Code development:** 1.5 weeks
- **Numerical tests + multi-charge Q-ball:** 1 week
- **Phenomenology check (new bosons):** 3 days
- **Total:** ~3.5 weeks

### Decision Criterion
**VIABILITY = YES if:**
- Multi-charge quantization is mathematically well-defined
- **AND** at least one charge pair (ω₁, ω₂) satisfies VK criterion
- **AND** no new light bosons required (or they decouple)
- **AND** mass ratios preserved within 0.1%
- **AND** no experimental tension

**Else:** Path N4b fails → pursue N4c/d

### Risk Level: 🟡 HIGH

---

## N4c: Disformal Sector Coupling (F-A Radiative Corrections)

**Question:** Do radiative corrections from the gravitational sector stabilize the solitonic saddle points?

### Mechanism
```
Isolated solitonic sector (CP-7):
  E_S[g] = ∫ [ ½f(g)|∇g|² + W(g) ] d³x
  W(g) = (β/3)g³ − (γ/4)g⁴,  β = γ
  → Tachyonic continuum from W''(1) = −1 < 0

Full theory (F-A + F-S coupled):
  E_total = E_A[ψ] + E_S[g] + E_coupl[ψ, g, σ_ab]
  
One-loop effective potential:
  V_eff^(1)(g) = W(g) + ΔV_1loop(g)
  ΔV_1loop(g) = (contributions from graviton/scalar box diagrams)
  
If ΔV_1loop(g) > W''(1) at g=1, then:
  V_eff''(1) = W''(1) + ΔV_1loop''(1) → possibly > 0
  → Tachyonic continuum removed ✓
```

### Pros
- ✅ **No axiom change** — coupling to F-A was always intended
- ✅ **Scheme-controlled** — perturbation theory, systematic order-by-order
- ✅ **Physical motivation** — gravity affects matter, especially in high-curvature regimes
- ✅ **Relatively contained** — only affects V_eff, not kinetic structure

### Cons
- ❌ **Technically demanding** — one-loop field theory integrals with Φ-dependent metrics
- ❌ **Scheme and scale dependent** — renormalization scale choice affects answer
- ❌ **UV divergences likely** — Φ⁴ theory is super-renormalizable, but couplings might need regulator
- ❌ **May not be enough** — why should one-loop radiative corrections >> tree tachyonic gap?
- ❌ **Convergence uncertain** — is one-loop sufficient? Or do multi-loop terms matter?

### Technical Requirements
1. Compute graviton propagator in curved background (based on M9.1'' metric)
2. Evaluate box diagram: graviton/scalar exchange contributions to V_eff(g)
3. Regulate divergences (dimensional regularization, MS-bar scheme)
4. Extract finite part; express as ΔV_1loop(g)
5. Test whether V_eff''(1) > 0 after renormalization

### Timeline
- **Feynman diagram setup + dimensional regularization:** 1 week
- **Loop integral computation (symbolic + numerical):** 1 week
- **Renormalization + scheme choice:** 3 days
- **Check V_eff spectral stability:** 3 days
- **Total:** ~2.5 weeks

### Decision Criterion
**VIABILITY = YES if:**
- One-loop ΔV_1loop can be computed in closed form (or accurately numerically)
- **AND** V_eff''(1) > 0 after renormalization
- **AND** choice of renormalization scale is justified (not ad-hoc)
- **AND** higher-loop estimates show one-loop is dominant

**Else:** Path N4c fails → pursue N4d

### Risk Level: 🟡 HIGH

---

## N4d: Metastability Hypothesis (Tunneling Dynamics)

**Question:** If μ/τ saddle points are real, can they be long-lived via WKB tunneling?

### Mechanism
```
Accept the saddle-point structure as physical. Then:
  1. Compute potential energy landscape around μ, τ profiles
  2. Identify nearest alternative minima (e.g., lower-energy configurations)
  3. Estimate tunneling rate Γ = A · exp(−S_bounce)
     S_bounce = ∫ 2[E_barrier − E_saddle] from saddle to barrier
  4. Lifetime: τ = 1/Γ
  5. If τ >> 13.8 Gyr (universe age), metastability is acceptable
     → TGP phenomenology unchanged; saddle structure is "hidden" by long lifetime
```

### Pros
- ✅ **No axiom change** — accepts saddle structure as-is
- ✅ **Fast to evaluate** — WKB tunneling is standard QFT tool
- ✅ **No core recalculation** — mass ratios, predictions unchanged
- ✅ **Phenomenologically conservative** — if tunneling is slow, μ/τ are effectively stable
- ✅ **Resolves paradox** — explains how saddle points can coexist with long-lived masses

### Cons
- ❌ **Interpretation shifts** — from "TGP predicts stable μ/τ" to "TGP predicts metastable μ/τ"
- ❌ **WKB accuracy** — tunneling rates are exponentially sensitive to barrier height; small errors → large Γ errors
- ❌ **Multi-dimensional tunneling** — μ/τ have 2/3 negative directions; which decay channel dominates?
- ❌ **Non-perturbative effects** — bouncing solutions might not exist in full theory; WKB is semi-classical approximation
- ❌ **Experimental signatures** — if τ decays (even slowly), there might be faint signals in rare-decay searches

### Technical Requirements
1. Scan potential energy landscape E_S[g(r, g₀)] around μ, τ
2. Find critical points (saddle point itself, nearby minima)
3. Compute tunneling path (bounce solution in imaginary time)
4. Integrate WKB exponent: S_bounce
5. Estimate tunneling rate from Γ ~ ω · e^{−S_bounce}
6. Compare τ decay with 13.8 Gyr and precision constraints

### Timeline
- **Landscape mapping:** 2 days
- **Bounce solution (analytical approximation or numerical):** 3 days
- **WKB exponent computation:** 2 days
- **Phenomenology assessment:** 2 days
- **Total:** ~1 week

### Decision Criterion
**VIABILITY = YES if:**
- Tunneling lifetime τ > 10^{20} seconds (>> 13.8 Gyr)
- **OR** tunneling is kinematically forbidden (no lower-energy minimum exists within profile family)
- **AND** decay channel (if any) is sufficiently suppressed

**Else:** Path N4d fails → TGP is ruled out by metastability

### Risk Level: 🟢 MEDIUM

---

## Summary Table

| Path | Mechanism | Difficulty | Risk | Timeline | Recommendation |
|------|-----------|-----------|------|----------|-----------------|
| **N4a** | Discrete substrate | Very High | 🔴 Very High | 3 weeks | **Last resort** — deep axiom change |
| **N4b** | Extended symmetry | High | 🟡 High | 3.5 weeks | **Pursue after N4c/d** — if alternative symmetries exist |
| **N4c** | Radiative corrections | Medium-High | 🟡 High | 2.5 weeks | **Pursue second** — addresses F-A coupling |
| **N4d** | Metastability | Low | 🟢 Medium | 1 week | **START HERE** — fastest, most conservative |

---

## Recommended Execution Order

### Phase 1: Fast check (N4d — Metastability)
1. **Goal:** Is there a nearby energy minimum? If yes, estimate tunneling lifetime.
2. **Time:** ~1 week
3. **Decision:** If τ > 10^{20} s → Continue Tier 2 with metastable interpretation
           If τ < 10^{15} s → Proceed to N4b/c

### Phase 2: Medium-term (N4c — Radiative Corrections)
1. **Goal:** Does F-A coupling stabilize F-S vacuum?
2. **Time:** ~2.5 weeks
3. **Decision:** If V_eff''(1) > 0 → Mass ratio check → Core revision
           If V_eff''(1) < 0 → Proceed to N4b

### Phase 3: Long-term (N4b — Symmetry Extension or N4a — Discrete Substrate)
1. **Goal:** Explore structural modifications
2. **Time:** 3–3.5 weeks each
3. **Decision:** If successful → Deep revision of TGP axioms
           If all fail → Admit TGP needs fundamental rethinking

---

## Critical Milestones

| Milestone | Trigger | Action |
|-----------|---------|--------|
| **N4d complete** | Day 7 | User decision: accept metastability, or escalate to N4c? |
| **N4c complete** | Day 18 | If radiative corrections fix spectrum → proceed to mass ratio check |
| **N4b/a complete** | Day 40–50 | Final synthesis: which path is viable? |
| **Tier 2 decision** | Day 50 | Go/no-go on revised TGP axioms |

---

## What Success Looks Like

**Tier 2 SUCCESS = At least one of:**
1. ✅ Metastability proven (N4d): tunneling lifetime >> 13.8 Gyr → Keep TGP, add "metastable interpretation"
2. ✅ Radiative stabilization (N4c): V_eff has minimum → Recompute masses, check r₂₁ within 0.1%
3. ✅ Symmetry extension (N4b): Multi-charge Q-ball stabilizes → Extend Lagrangian, update core
4. ✅ Discrete substrate (N4a): Tachyonic band vanishes at continuum limit → Reinterpret axiom I

**Tier 2 FAILURE = All paths fail:**
→ Admit CP-7 saddle structure is fundamental incompatibility
→ TGP needs deep revision (new field? new metric? new dynamics?)
→ Enter "restructuring" phase (not covered in Tier 2)

---

## Cross-References

- [[./TIER2_ANALYSIS_FRAMEWORK_2026-07-27.md]] — Overall framework
- [[./CP7_ARTIFACT_ANALYSIS_2026-07-27.md]] — Artifact vs. structural analysis
- [[../research/op-spectral-analysis-Phi-2026-07-03/README.md]] — CP-7 findings
- [[../research/op-wall-dynamics-2026-07-03/README.md]] — Wall dynamics (N6)
- [[../research/op-nonlinear-charge-constraint-2026-07-03/README.md]] — Q-ball test (context for N4b)

