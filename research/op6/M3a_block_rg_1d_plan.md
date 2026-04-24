# OP-6 / M3-a — Block-RG bond-form test of H₁ in 1D

**Status:** working, 2026-04-24.
**Predecessors:** `M2b_results.md`, `M3c_scaling_dimensions.md`.
**Constraint:** preserve Z₂-odd microscopic field ŝ; H₁ as in M1 §1.
**Question:** under explicit numerical block-RG of H₁ in 1D, does the
effective kinetic stiffness K_B(⟨Φ_B⟩) flow towards **K(Φ) ∝ Φ**
(p_B → +1, TGP target / M1-B) or towards **K(Φ) = const**
(p_B → 0, M3-c analytical prediction)?

---

## 1. Logic of the test

M2-b measured K_eff(⟨Φ⟩) at the lattice scale (B = 1, single-site
Φ_i = ŝ_i²) and found p = −0.28 ± 0.43, rejecting the TGP target
p = +1 at 3σ. That is a one-scale measurement; it does not exclude
the possibility that **flow** under block-RG drives the exponent
towards +1 at larger scales.

M3-c's scaling-dimension argument at 3D Ising WF predicts the
opposite: **K^{(0)} dominates K^{(1)} more strongly in the IR**,
ratio `K^{(1)}/K^{(0)} ~ k^{Δ_ε} = k^{1.41}` → 0. The analogous
statement in 1D (no WF fixed point, just Gaussian flow) is softer
but qualitatively the same: the dimension-zero bond stiffness
dominates higher-Φ-order bond structures at coarse scales.

M3-a tests this directly by measuring K_B(⟨Φ_B⟩) at a **ladder of
block sizes** B ∈ {1, 2, 4, 8, 16} from a single long MC simulation,
and fitting p_B for each B. Three possible outcomes:

| Flow pattern         | Prediction of           | Verdict for M1-B            |
|----------------------|-------------------------|------------------------------|
| p_B → +1 as B grows  | M1-B (H₁ → H_GL in IR)  | M1-B **revived**             |
| p_B stays ≈ −0.28/0  | M3-c (K^{(0)} dominates)| M1-B **closed on 2 grounds** |
| p_B bifurcates/chaos | (unexpected)            | Investigate                  |

## 2. Method

### 2.1 Simulation

- Hamiltonian: H₁ = Σ_i [(m²/2)ŝ² + (λ/4)ŝ⁴] − J Σ_⟨ij⟩ ŝ_iŝ_j .
- Parameters: J = 1, λ = 1 (matching M2-b Part C). Scan
  (β, m²) grid same as M2-b: β ∈ {0.3, 0.5, 0.7, 1.0, 1.5, 2.0,
  3.0}, m² ∈ {−1, 0, 1, 2, 3}. That is 7 × 5 = 35 runs.
- Lattice N = 4096 (2× M2-b's N = 1024, so B = 16 blocks give
  N/B = 256 block sites — still plenty for FFT).
- Thermalisation / measurement: n_therm = 3000, n_measure = 8000,
  measure every 2 sweeps. Same as M2-b.

### 2.2 Block averaging

For block size B, define **block-averaged Φ field**:

    Φ_B,I  =  (1/B) · Σ_{i ∈ block I}  ŝ_i² ,   I = 0, …, N/B − 1 .

(Block averaging rather than decimation: averaging is the standard
block-RG of a scalar field and preserves the sum rule for Φ.)

### 2.3 Block-level OZ fit

Same OZ protocol as M2-b §1, applied to the block field:

    1 / S_B(k̂_B)  ≈  1/χ_{Φ_B}  +  K_B · k̂_B² ,
    k̂_B²  =  2(1 − cos(2π j / (N/B)))  ,   j = 0, …, N/(2B) .

Fit the lowest n_k_fit = 4 non-zero modes. Extract K_B. Error via
jackknife (n_jk = 8), identical to M2-b.

### 2.4 Power-law fit

For each B, fit across the (β, m²) scan:

    K_B(⟨Φ_B⟩)  =  C_B · ⟨Φ_B⟩^{p_B} .

Report (p_B, stderr(p_B), n_B). Plot p_B vs B.

## 3. Predictions

### 3.1 Quantitative prediction — M3-c

At 1D (Gaussian flow), the scaling-dimension ratios don't apply
directly (no anomalous dimensions). But the *qualitative* statement
— that the Φ-independent bond stiffness dominates the Φ-dependent
one — is still expected. So:

    p_B  ≈  −0.3 to +0.1  for all B ∈ {1, 2, 4, 8, 16}  (M3-c-consistent)

### 3.2 Quantitative prediction — M1-B

If block-RG in 1D does generate K(Φ) ∝ Φ, we should see:

    p_B → +1  as B grows  (M1-B revival)

A visible monotonic drift p_1 → +1 over the block ladder would be
positive evidence. A flat p_B ≈ −0.3 would be negative evidence.

## 4. Cost and schedule

- Estimated wall time: ≈ 2–5 min (reuses M2-b's H₁ MC; block
  averaging is a free numpy reshape; extra cost is 5 OZ fits per
  parameter point instead of 1).
- Exit criterion: clear result for p_B(B) profile; decision on
  whether to proceed to M3-b (3D cluster MC) or reopen the axiom
  conversation.

## 5. Caveats

- 1D has no Wilson–Fisher; this test probes Gaussian-flow block-RG,
  not the WF fixed point. A 1D `flat p_B` result does not *prove*
  p_B stays flat at 3D WF — but combined with M3-c's analytical
  argument, it tightens the noose.
- OZ fit breaks down when ξ_{Φ_B} ≲ 1 in block units (i.e. ξ_Φ ≲ B
  in lattice units). For small ξ_Φ (M2-b had ξ_Φ ≤ 3), large-B
  fits will have correspondingly fewer usable modes. Expect B = 16
  measurements to be the noisiest; B = 1, 2 most reliable.
- The block field Φ_B is non-Gaussian (sum of B lattice Φ values,
  Φ_i ≥ 0). This does not invalidate the OZ expansion (still a
  second-derivative of the effective action at k = 0), but the
  log-log fit may acquire B-dependent finite-range corrections.

## 6. Exit criteria

- [ ] p_B computed for B ∈ {1, 2, 4, 8, 16} across the full scan.
- [ ] Tabulated (B, p_B, stderr, n_valid, χ²_log-log) table.
- [ ] Verdict recorded: *flat p_B* (M3-c confirmed) vs *drifting
      p_B → +1* (M1-B revived) vs *unexpected*.
- [ ] Recommendation on M3-b (3D MC) vs axiom-reopen.
