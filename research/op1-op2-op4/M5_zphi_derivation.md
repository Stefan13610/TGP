# M5 — Z_Φ / η-deformation of MK-RG (analytical setup)

**Date:** 2026-04-25.
**Test designation:** "Test B" (P3.1) of
`external_review_2026-04-25/review_response_plan.md`.
**Scope:** analytical setup supporting `mk_rg_zphi.py`.
Conclusion in `M5_results.md`.

---

## 1. Hypothesis under test

M3 (`mk_rg_bgamma.py`) found `B*/Γ* ≈ -0.57` at the 3D-Ising WF FP,
with target `1/v*² ≈ +1.23`. M4 (Test A, `mk_rg_phi.py`) ruled out
the Hubbard–Stratonovich Jacobian as the missing physics.

This note tests the **second candidate** from the M3 §6 list:
**wave-function renormalisation** of the composite field `Φ = ŝ²`.

Standard NPRG-LPA' tracks both `V_k(Φ)` (potential) and `Z_k`
(field-strength renormalisation), with anomalous dimension

```
η_Φ = -d/dt ln Z_k|_FP        (t = ln(1/k), k = RG scale)
```

For 3D Ising bootstrap: `η_φ ≈ 0.036`, hence `η_Φ = η_{φ²} = 2 η_φ ≈ 0.07`
(at one-loop the composite picks up ≈ twice the elementary η).
This is small in standard FRG, but its effect on B*/Γ* could still
be O(1) because B and Γ couple to `Φ³` and `Φ⁴` respectively, and
thus are sensitive to `Z_k^3` and `Z_k^4` factors.

**Question:** does adding an η_Φ deformation to M3's RG step shift
the polynomial FP enough to recover `B*/Γ* = 1/v*²`?

If yes → η_Φ is the missing physics; OP-2b closes.
If no → η_Φ alone insufficient; the GL bond (P3.2) is the
remaining single-channel candidate.

## 2. Where η enters the MK-RG step

In single-site MK-RG with bilinear bond `−K Σ s_i s_j`, the basic
RG step is:

```
1. Bond move:    K → K_eff = b^{d-1} K
2. Decimate:     V'(s) = V(s) − 2 [F(s) − F_0]
3. Bar-rescale:  cb_{2k} = c_{2k} / K_new   with K_new = F_2
```

The third step implicitly fixes the field renormalisation by setting
`K_new ≡ 1` each iteration. In standard MK with no anomalous dim,
this corresponds to `Z_φ = 1` per step.

To introduce η_φ ≠ 0, we generalise step 1:

```
K_eff(η_φ) = b^{d-1+η_φ} K
           = b^{d-1} K · b^{η_φ}
           = (M3 K_eff) × b^{η_φ}.
```

**Interpretation.** Under MK coarse-graining with anomalous field
dim `η_φ`, the field rescales as `s → b^{−(d−2+η_φ)/2} s`. The
bilinear bond `−K s_i s_j` then absorbs a factor `b^{d−2+η_φ}` from
the field rescaling (combined with `b^{d-1}` from bond move and the
canonical `b^{−(d−2)}` from decimation cancellation), netting

```
K_eff_full = K · b^{d-1+η_φ}.
```

The factor `b^{η_φ}` is the explicit η-deformation. For η_φ = 0 we
recover M3 exactly.

**Caveat — composite vs. elementary.** The relevant anomalous dim
for the cubic (B) and quartic (Γ) operators in `Φ = ŝ²` is the
**composite η_Φ**, not the elementary `η_φ`. The relation
η_Φ = 2 η_φ holds at one-loop in φ⁴ theory but is corrected at
higher orders. In this test we treat η_φ as the **single deformation
parameter** (since it enters MK-RG through the bilinear-bond
rescaling) and report what value, if any, gives `B*/Γ* = 1/v*²`. The
result can then be compared to `2 × 0.036 ≈ 0.07` (3D Ising
bootstrap).

## 3. Implementation: a one-line modification of M3

In `mk_rg_phi.py` we already have `MigdalKadanoffPhi` extending
`MigdalKadanoffRGN`. Here we extend it once more by overriding only
the `K_eff` factor in `rg_step`:

```python
class MigdalKadanoffZPhi(MigdalKadanoffRGN):
    def __init__(self, n_ops, eta=0.0, ...):
        super().__init__(...)
        self.eta = float(eta)

    def rg_step(self, couplings, K):
        K_eff = (BD1 * K) * (B_RESCALE ** self.eta)
        ...   # rest identical to M3
```

η = 0 must reproduce M3 to 5 decimals (sanity check).

## 4. Decision criteria for OP-2b (Test B)

Define `R(η) ≡ B*(η)/Γ*(η)` and `T(η) ≡ 1/v*²(η)`.

1. **Closure (positive):** `R(η_*) = T(η_*)` for some η_* with
   `|η_*| ≲ 0.1` (i.e., physically reasonable, comparable to the
   bootstrap value `2 η_φ ≈ 0.07`). η_Φ is then THE missing physics.

2. **Closure-in-principle (intermediate):** `R(η_*) = T(η_*)` for some
   η_* but only at unphysically large `|η_*| > 0.5`. η_Φ contributes
   in the right direction but is not the dominant mechanism — the
   GL bond is still required.

3. **No closure (negative):** `R(η)` does not cross `T(η)` for any
   η in the convergent regime, OR signs of R(η) and T(η) never
   align (analogous to M4's outcome). Then η_Φ alone cannot close
   OP-2b; the GL bond is the remaining candidate.

## 5. Convergence regime

For η > 0 (positive anomalous dim), `K_eff = 4 · 2^η > 4`. Larger
bond strength means narrower distribution per site, easier
convergence. Numerically safe.

For η < 0, `K_eff < 4`. As η → −∞, K_eff → 0, the bond decouples,
and the FP iteration loses the WF branch (degenerates to single-site
Gaussian). The convergent regime is bounded below by some
`η_min ≈ −1` (`K_eff ≈ 2`).

We scan `η ∈ [−1, 2]`, with finer resolution around the bootstrap
value `η ≈ 0.07` and around any sign-flip of `R(η) − T(η)`.

## 6. Expected scaling

A simple analytic estimate: at the FP, the cumulants `κ_{2k}` scale
roughly as powers of `K_eff` (from M3's update rule
`F_{2k} = K_eff^{2k} κ_{2k}` and the bar-norm condition `K_new = 1`,
giving `κ_2 = K_eff^{-2}`). For higher cumulants, the FP balances
`c*_{2k} = 2 F*_{2k}/(2k-1)!`, so `c*_{2k} ∝ K_eff^{2k}`.

Therefore: `c*_{2k} ∝ K_eff^{2k} = (4 b^η)^{2k}`. In bar variables
this rescales as `cb*_{2k} = c*_{2k}/K* = c*_{2k}` (since K_new=1
in bar mode). So the η-dependence of cb*_{2k} is approximately
`cb*_{2k}(η) ≈ cb*_{2k}(0) · b^{2k η}`.

**Predictions**:
- `B*(η)/Γ*(η)` ≈ `B*(0)/Γ*(0) · b^{−2η}` (since Γ has 2k=8 vs B has 2k=6).
- `v*²(η) = |r*|/u*` ≈ `v*²(0) · b^{−2η}` (r has 2k=2, u has 2k=4).
- `1/v*²(η) ≈ 1/v*²(0) · b^{2η}`.

**At b=2, η=0.07 (bootstrap)**: `2^{0.14} ≈ 1.10`. So the η_Φ-induced
shift is at the 10% level, while the M3 gap is ~300%. This is a
**rough analytical estimate that η_Φ at the bootstrap value is
two orders of magnitude too small to close OP-2b.** The numerical
test will quantify this precisely.

The interesting question is whether **larger η** (un-bootstrap-like
but consistent with some lattice MK-RG anomaly) can close it. For
B*/Γ* to flip sign requires either dramatic FP shift or genuinely
new physics from the bond.

## 7. Files

- `mk_rg_zphi.py` — implementation (extension of `mk_rg_bgamma.py`).
- `mk_rg_zphi_results.txt` — raw output.
- `M5_results.md` — verdict on η_Φ closing OP-2b.
