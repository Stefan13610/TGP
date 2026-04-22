# DESI w(z) Phantom-Crossing Falsification — Program P2 (TGP Redirect)

> **STATUS 2026-04-19: P2 ACTIVE — FORMAL DERIVATION + FALSIFICATION ROADMAP COMPLETE.**
> This is program **P2** from the [[REDIRECT_PROGRAM_2026-04-19.md|post-closure redirect]].
> TGP makes a **structural, derived, falsifiable** prediction for dark energy: **w(z) ≥ −1 for all z, no phantom crossing.** DESI DR2 (2026 Q3) is the decisive test point.

---

## 1. The TGP-unique prediction (derived, not fit)

From the sek05 potential `U(ψ) = (β/3)ψ³ − (γ/4)ψ⁴` with β = γ:
- `U(1) = γ/12`, `U'(1) = 0`, `U''(1) = −γ < 0` (slow-roll maximum)
- Action has **canonical kinetic term** (verified gs08b — no ghost degrees of freedom)

**Theorem (de2):** For a canonical scalar with positive kinetic term in FRW,
```
w + 1 = ψ̇² / ρ_ψ  ≥  0     (always)
```
Phantom crossing `w < −1` would require wrong-sign kinetic (ghost instability), which is **structurally forbidden** in sek05.

**Numerical verification (de2):** global `w_min = −1.00000000` across all trajectories (δ ∈ {1e-4, 1e-3, 1e-2}).

**CPL projection of TGP trajectory:**

| Parameter | TGP (sek05 derived) | ΛCDM | DESI DR1 DES-SN5YR |
|----------|---------------------|------|---------------------|
| w₀ | **−0.99998** | −1.00 | −0.727 ± 0.067 |
| wₐ | **−0.00003** | 0.00 | −1.050 ± 0.310 |

TGP is **indistinguishable from ΛCDM** within ~10⁻⁴ on (w₀, wₐ). Both sit at **~4.1σ tension** with DESI DR1 combined (DES-SN5YR).

---

## 2. Why P2 is a clean falsifier

1. **TGP-unique**: no adjustable parameters; w(z) follows from sek05 ansatz.
2. **Orthogonal to galaxy-MOND**: does not rely on the ν(y) claim retracted in [[gs65_formal_derivation.py|gs65]]/[[gs66_frw_propagator.py|gs66]].
3. **Cannot be rescued**: bridge-(b)/(c) substrate modifications affect structure growth, not the sign of kinetic energy; the w ≥ −1 bound is a theorem.
4. **Falsifier timeline**: DR2 in ~6 months, DR3 in ~18 months.

---

## 3. Program files

| File | Role | Status |
|------|------|--------|
| `de1_tgp_w_prediction.py` | Unified tension scan (H₀/S₈/w(z)); first look at DESI | ✔ ran 2026-04-18 |
| `de2_tgp_frw_evolution.py` | **Full FRW integration from sek05**; derived w₀, wₐ, w_min | ✔ ran 2026-04-19 |
| `de3_desi_dr1_comparison.py` | Mahalanobis distance to DESI DR1 posterior; AIC/BIC; Bayes | ✔ ran 2026-04-19 |
| `de4_falsification_roadmap.py` | DR2/DR3 projection, falsification criteria, timeline | ✔ ran 2026-04-19 |

Outputs: `de{1,2,3,4}_results.txt`

---

## 4. Consolidated results

### de2 — FRW derivation
- Integrated coupled ODE: `ψ̈ + 3Hψ̇ + V′(ψ) = 0`, `H² = Ω_r a⁻⁴ + Ω_m a⁻³ + ρ_ψ/ρ_crit`.
- Normalized `V(ψ) = 4ψ³ − 3ψ⁴`, shooting parameter V₀ fixed by Ω_DE0 = 0.685.
- Initial conditions at z=99: `ψ = 1 − δ`, `u = 0` (Hubble-frozen).
- **Result:** w(z) stays within 10⁻⁴ of −1 for all realistic δ; TGP ≈ ΛCDM on cosmology.

### de3 — DR1 comparison (combined posterior `(w₀, wₐ) = (−0.727 ± 0.067, −1.05 ± 0.31)`, ρ = −0.85)

| Model | w₀ | wₐ | χ² | n-σ | ΔAIC | ΔBIC |
|---|---|---|---|---|---|---|
| ΛCDM | −1.000 | 0.000 | 16.62 | 4.08σ | +12.62 | +15.24 |
| **TGP (de2)** | −0.99998 | −0.00003 | 16.62 | **4.08σ** | +12.62 | +15.24 |
| CPL best-fit | −0.727 | −1.050 | 0.00 | 0.00 | 0.00 | 0.00 |

Bayes factor log B(CPL/TGP) = 2.46 → "positive" for CPL (Jeffreys scale), but short of decisive. **SN-sample spread** (Pantheon+ 2.2σ, Union3 3.1σ, DES-SN5YR 4.1σ) is the dominant systematic.

### de4 — Falsification roadmap

**DR2 (~2026 Q3, √3 error reduction):** 

| Scenario | Projected n-σ (DES-SN5YR-like) | Verdict |
|---|---|---|
| (i) central migrates to LCDM (w₀=−0.9, wₐ=−0.3) | 2.77σ | TGP in tension, not falsified |
| (ii) DR1 persists (w₀=−0.727, wₐ=−1.05) | **7.06σ** | **TGP FALSIFIED** |
| (iii) drifts further (w₀=−0.6, wₐ=−1.5) | 10.4σ | TGP catastrophically falsified |

**DR3 (~2027 Q3, √5 error reduction):**

| Scenario | Projected n-σ | Verdict |
|---|---|---|
| (i) | 3.57σ | TGP marginally falsified |
| (ii) | **9.12σ** | **TGP falsified** |
| (iii) | 13.4σ | TGP catastrophically falsified |

**Even scenario (i)** — DESI central values migrating toward ΛCDM — triggers 3σ falsification at DR3, unless migration continues all the way to exactly (w₀, wₐ) = (−1, 0).

**Survival probability estimate (given DR1 tension):**
- P(TGP survives DR2) ≈ 30–50%
- P(TGP survives DR3) ≈ 15–30%

---

## 5. Specific falsification criteria

- **(F1) DR2 DES-SN5YR gives > 3.5σ + Pantheon+/Union3 both > 2.5σ** → TGP falsified at cosmological scale.
- **(F2) DR3 gives > 3σ across all three SN samples** → TGP definitively falsified.
- **(F3) Direct binned w(z) measurement detects w(z) < −1 at > 3σ in any bin** → immediate falsification.
- **(S1) DR2/DR3 converges to (w₀, wₐ) ≈ (−1, 0) within errors** → TGP supported (equivalent to ΛCDM at observational level).

---

## 6. Timeline

| Quarter | Milestone |
|---|---|
| 2026 Q2 | Program P2 scripts complete (de1–de4). **← current** |
| 2026 Q3 | **DESI DR2 release → decisive test point (F1).** |
| 2026 Q4 | DR2 + CMB + full SN joint re-analysis. |
| 2027 Q2–Q4 | **DESI DR3 full survey release → final test (F2).** |
| 2028+ | Euclid + LSST + CMB-S4 independent cross-check. |

---

## 7. What P2 does NOT address (scope)

- **H₀ tension**: TGP structurally cannot resolve (gs20 + de1; ~8 orders of magnitude too weak). *Not in scope.*
- **S₈ tension**: TGP substrate effect ~0.001%, 8500× too weak ([[TGP/TGP_v1/research/s8_tension/]]). *Not in scope.*
- **Galaxy-MOND**: retracted as first-principles derivation ([[CLOSURE_2026-04-19.md]]). *Separate program.*
- **Clusters**: handled by P1 (`research/cluster_anomaly/` (planowane) + sterile ν hybrid).

P2 is narrowly targeted on the **phantom-crossing falsification** — the cleanest, fastest, most TGP-unique cosmological test available.

---

## 8. Key references

- DESI Collaboration DR1 (arXiv:2404.03002, 2024)
- DESI Collaboration DR2 (expected 2026 Q3)
- Chevallier & Polarski 2001, Linder 2003 — CPL parametrization
- `research/galaxy_scaling/gs08b_ghost_free_verification.py` (planowane) — kinetic sign verification in sek02+sek05
- [[REDIRECT_PROGRAM_2026-04-19.md]] — parent redirect program
- [[TGP_STATUS_2026-04-19.md]] — overall TGP status post-closure

---

**Bottom line:** If DESI DR2 confirms DR1 central values at √3-shrunken errors, TGP is falsified on dark energy at ~7σ within 6 months. If DR2 migrates toward ΛCDM, TGP is supported. **Either way, by 2027 TGP cosmology has a definitive answer.**
