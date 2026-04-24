# M6 — GL-bond operator in MK-RG (verdict)

**Date:** 2026-04-25.
**Test designation:** Test C (P3.2) of
`external_review_2026-04-25/review_response_plan.md`.
**Verdict:** **CLOSURE-IN-PRINCIPLE only.** A crossing exists at
J_GL ≈ +5.89, but it is **outside the perturbative regime** in
which Track A's first-order-in-J_GL truncation is justified.
Track A alone cannot decide whether the GL bond closes OP-2b at
the non-perturbative level.

Setup: `M6_glbond_derivation.md`.
Code:  `mk_rg_glbond.py`.
Raw:   `mk_rg_glbond_results.txt`.

---

## 1. Summary of result

Adding the GL bond `J_GL · Σ A_ij Φ_i Φ_j (Φ_j − Φ_i)²` as a fixed
parameter (Track A) and treating it to first order in J_GL in the
on-site V update of the M3 MK-RG step:

| Regime | J_GL | r* | u* | B*/Γ* | 1/v*² | diff |
|---|---|---|---|---|---|---|
| M3 baseline | 0 | -2.45694 | +3.01611 | -0.56869 | +1.22759 | -1.79628 |
| Perturbative max | +1.000 | -2.46040 | +3.07160 | -0.39124 | +1.24841 | -1.63965 |
| Mid-range | +3.000 | -2.47176 | +3.18271 | +0.09019 | +1.28763 | -1.19744 |
| Mid-range | +5.000 | -2.49323 | +3.29974 | +0.85629 | +1.32348 | -0.46719 |
| **Crossing** | **+5.893** | **-2.5079** | **-** | **≈ +1.339** | **≈ +1.339** | **≈ 0** |
| Past crossing | +6.000 | -2.50993 | +3.36612 | +1.40361 | +1.34112 | +0.06249 |
| Far past | +10.000 | -2.61354 | +3.76329 | +4.51873 | +1.43992 | +3.07881 |

Negative J_GL pushes diff *more negative* (away from closure):
J_GL = -1 gives diff = -1.924, J_GL = -2 gives diff = -2.029,
and convergence degrades for J_GL ≲ -2.

**Crossing:** linear interpolation on the fine scan
J_GL ∈ [4.5, 7.0] step 0.05 places the crossing at
**J_GL_* = +5.8933 ± 0.05** (bracket [5.85, 5.90]).

## 2. Validation

`J_GL = 0` reproduces the M3 8-operator FP to 5 decimals:

```
M3 expects:  r* = -2.45694, u* = +3.01611, B*/Γ* = -0.5687
M6 J=0:      r* = -2.45694, u* = +3.01611, B*/Γ* = -0.56869   ✓
```

Tiny-J_GL continuity is smooth:

```
J_GL = +1e-3:  B*/Γ* = -0.56853   (shift +1.6e-4 vs J=0)
J_GL = -1e-3:  B*/Γ* = -0.56885   (shift -1.6e-4 vs J=0)
```

## 3. Why the crossing is outside the perturbative regime

Track A truncates the J_GL expansion at first order:

```
F(s_1, s_3) ≈ F_M3(K_eff(s_1+s_3)) − J_GL ⟨bond_12 + bond_23⟩
            + O(J_GL²).
```

For this truncation to be controlled, the per-step correction must
be small relative to F_M3. At the M3 fixed point:

- Typical scale: ⟨bond_12⟩(K_eff·s) ≈ ⟨s_2² s_1^a⟩ * (combinatorial)
  ~ O(1) at the WF FP for s_1 ~ s* ~ 1.
- Per-step correction at J_GL ≈ 6: J_GL · ⟨bond⟩ ~ 6, comparable to
  the M3 polynomial sector itself (|c_4*| ≈ 3, |c_6*| ≈ 5).
- Hence at the crossing J_GL ≈ 5.89, **first-order J_GL truncation
  is grossly violated**; J_GL² and higher terms are O(1) corrections
  to the answer, not small.

The **perturbative regime** for Track A is roughly |J_GL| ≲ 1.
Inside this regime:

- diff(0)   = -1.796
- diff(0.5) = -1.722  (4 % closure)
- diff(1.0) = -1.640  (9 % closure)

i.e., the GL bond at first order closes the M3 gap at the **few-%
per unit J_GL** rate. Reaching closure requires |J_GL| ~ 6, which
is outside Track A's domain of validity.

## 4. Decision against M6 §6 criteria

| Criterion | Outcome |
|---|---|
| (1) Closure with |J_GL_*| ≲ 1 (perturbative) | ✗ NOT satisfied |
| (2) Closure-in-principle: crossing at unphysical |J_GL_*| | ✓ at J_GL_* ≈ +5.89 |
| (3) No crossing in convergent regime | ✗ (crossing exists) |

So the result is **(2) closure-in-principle**. The GL bond moves the
flow in the right direction (B*/Γ* climbs from −0.57 through 0 to
+4.5 across J_GL ∈ [0, 10], crossing the target +1/v*² ≈ +1.34 at
J_GL ≈ 5.89). But Track A's first-order truncation cannot be trusted
near the crossing, so we cannot conclude OP-2b is closed.

## 5. Comparison with M3, M4, M5

| Test | Channel | Result | Min |diff| | Closure J/μ/η |
|---|---|---|---|---|
| M3 | bilinear bond only | gap | 1.796 (J=0) | — |
| M4 (Test A) | + H-S Jacobian μ ln(s² + ε²) | NEGATIVE | ≈ 1.62 at μ ≈ 0.32 | none |
| M5 (Test B) | + Z_Φ via η-deformation | NEGATIVE | 1.694 at η ≈ -0.55 | none |
| **M6 (Test C)** | **+ GL bond, Track A** | **CLOSURE-IN-PRINCIPLE** | **0.004** | **J_GL ≈ +5.89** |

M6 is qualitatively different from M4/M5: those two scans bound
|diff| from below (the gap never closes for any value of the
deformation parameter). M6 *does* close the gap, but only at
parameter values where Track A's perturbative assumption fails.

This is genuinely informative: the GL bond is the **only**
single-channel candidate from the M3 §6 list that drives the flow
through the target. The other two (H-S Jacobian, Z_Φ) are simply
incapable of closing the gap regardless of how large the deformation
parameter is taken.

## 6. What this does (and does not) prove

**It does prove:**
- The GL bond at finite J_GL > 0 monotonically reduces the M3
  gap |B*/Γ* − 1/v*²|, eventually crossing zero.
- At least at the level of single-site-projected MK-RG with Track-A
  truncation, the GL bond is the *unique* single-channel mechanism
  among {H-S Jacobian, Z_Φ, GL bond} that can close OP-2b.
- The closure direction is **physically correct**: positive J_GL
  (gradient bond favouring smooth Φ-configurations) pushes B/Γ
  toward the v² = β/γ Lorentz-locking value, exactly as expected
  on physical grounds.

**It does not prove:**
- That closure persists when J_GL² and higher terms are included
  (Track A is incomplete near the crossing).
- That the closure value J_GL_* ≈ 5.89 is the same as the bare
  Hamiltonian J_GL of the v2 substrate (since Track A is in the
  on-site projection, not the full bond flow).
- That the closure persists when J_GL flows under MK along with
  (r, u, B, Γ, …) — Track B (J_GL flowing) was deferred from this
  test.

## 7. What's next

Three options remain:

1. **Track B (J_GL flowing).** Extend `mk_rg_glbond.py` so that
   J_GL is updated under MK alongside (r, u, B, Γ, …). Operator-mixing
   contributions from `bond_12 ↔ V` decimation enter at the level
   of J_GL² (which Track A drops). Closure under Track B would be
   strong evidence for genuine GL-bond closure of OP-2b.

2. **NPRG (Wetterich) with full Z_Φ + GL kinetic ansatz (P3.4).**
   The non-perturbative gold standard. Resolves the perturbative-
   regime ambiguity by working at all orders in the couplings.

3. **Accept the present "closure-in-principle" status.** Document
   that:
   - All three single-channel candidates have been tested.
   - Two (H-S, Z_Φ) cannot close OP-2b at any value of their
     deformation parameter.
   - One (GL bond) does close OP-2b at finite J_GL > 0, but
     outside the perturbative regime where the test is fully
     trustworthy.
   - The gap is therefore **most likely** the GL bond, but a
     definitive proof requires either Track B or NPRG.

The recommendation: report this result honestly in the v3 paper
patch (status of OP-2b moves from "open" to "closed-in-principle
via GL bond, pending Track B or NPRG confirmation"). Schedule
NPRG (P3.4) as the next item.

## 8. Files

- `mk_rg_glbond.py` — Track A implementation (J_GL fixed, on-site
  projection, first-order truncation).
- `mk_rg_glbond_results.txt` — full scan output:
    - validation J_GL = 0 reproduces M3 ✓
    - tiny-J_GL continuity ✓
    - coarse scan J_GL ∈ [-2, 10] step 0.25
    - fine scan |J_GL| ≤ 1 step 0.05 (perturbative)
    - fine scan J_GL ∈ [4.5, 7] step 0.05 (around crossing)
    - crossing at J_GL = +5.8933 (linear interp)
    - min |diff| = 0.0039 at J_GL = +5.90
- `M6_glbond_derivation.md` — analytical setup.
- `M6_results.md` — this verdict.

## 9. Bottom line

**Track A:** the GL bond does cross the closure target
(B*/Γ* = 1/v*²) at J_GL_* ≈ +5.89, monotonically through the
M3 baseline, but at a value where first-order J_GL truncation is
no longer reliable.

**Status of OP-2b after M3, M4, M5, M6:**
- *Single-site MK + bilinear bond only* (M3): gap −1.80.
- *+ H-S Jacobian* (M4): gap unchanged (no closure).
- *+ Z_Φ via η-deformation* (M5): gap unchanged (no closure).
- *+ GL bond, Track A* (M6): closure at J_GL ≈ 5.89, **outside
  the perturbative regime**.

The GL-bond candidate is the **strongest** of the three; the
remaining work is to confirm or deny it at the non-perturbative
level (Track B or NPRG).
