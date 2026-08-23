# CHECKPOINT: Phase 4 Extended Implementation Plan

**Status:** Ready for Phase 4 Extended (Full Numerical Test)  
**Date Created:** 2026-07-27  
**Last Updated:** 2026-07-27  
**Location:** `op-native-pressure-lepton-stability-2026-07-27/`

---

## 🎯 Current Position

### Completed Work (Phases 1-4 Qualitative)
✅ Phase 1: Toy model validates pressure mechanism  
✅ Phase 2: Charges extracted from CP-7 profiles  
✅ Phase 3: Three-body equilibrium found and verified  
✅ Phase 4: Qualitative sign analysis confirms repulsive pressure  

### Next Phase (Phase 4 Extended - Numerical)
⏳ Integrate pressure term into CP-7 BVP solver  
⏳ Compute L̂_eff = L̂_TGP + V_bg(r)  
⏳ Diagonalize and extract exact spectral shifts Δλ  
⏳ Compare: isolated solitons vs. with pressure  

**Estimated effort:** 1-2 weeks  
**Critical files:** CP-7 Phase2_bvp_spectrum.py (spectral solver)

---

## 📋 Current State Summary

### Constants & Parameters (Verified)

```python
# CP-7 Formulacja B (log-form)
PHI = (1 + np.sqrt(5)) / 2
G0_E = 1.24915      # Electron reference
G0_MU = 2.02117     # Muon (φ * G0_E)
G0_TAU = 3.18912    # Tau (Koide-determined)

# Extracted Charges (Phase 2)
Q_E = 1.20009652
Q_MU = 1.10685470
Q_TAU = 1.04924277

# Equilibrium Configuration (Phase 3)
R_SCALE = 1e7  # Physical Planck units
R_EQ_EMU = 32472644.977 / R_SCALE  ≈ 3.247
R_EQ_ETAU = 15403380.105 / R_SCALE  ≈ 1.540
R_EQ_MUTAU = 21302469.877 / R_SCALE  ≈ 2.130

# Pressure Energy
E_PRESS_EQ = -5.046  # Binding energy
```

### Baseline Spectral Results (CP-7, Isolated)

```
Soliton e (g₀=1.24915):
  l=0: λ_min = -0.998,   N_neg = 0    (STABLE)
  l=1: λ_min = -0.995,   N_neg = 0

Soliton μ (g₀=2.02117):
  l=0: λ_min = -1.282,   N_neg = 2    (SADDLE)
  l=1: λ_min = -1.031,   N_neg = 1

Soliton τ (g₀=3.18912):
  l=0: λ_min = -4.216,   N_neg = 3    (SADDLE)
  l=1: λ_min = -1.108,   N_neg = 1
```

**Source:** CP-7 Phase2_bvp_spectrum.py output (Lines: [C4/C6a] section)

---

## 🔧 Phase 4 Extended: Detailed Implementation Roadmap

### Step 1: Import & Understand CP-7 Solver

**File Location:** `../op-spectral-analysis-Phi-2026-07-03/Phase2_bvp_spectrum.py`

**Key Functions to Study:**
```python
def soliton_profile(formulation, g0, R=60.0, N=4000):
    """
    Generates soliton profile g(r), g'(r), g''(r)
    Input: formulation ('F-S' or 'F-A'), initial g₀, domain R, grid points N
    Output: r, g, gp, d2, bounces, g_min
    """

def spectrum_on_background(formulation, g0_bg, g0_test, l_angular, mu_list, N_pts):
    """
    Computes spectrum of test soliton on background field
    Input: formulation, background g₀, test g₀, angular momentum l, eigenvalue guesses, grid
    Output: eigenvalues, eigenvectors, residuals
    """

def phase_shift_l0(g0, R_phase=60.0, N_pts=4000):
    """
    Computes scattering phase shift (used to count negative modes)
    For l=0, counts box spectrum intersection with tachyonic continuum
    """
```

**What You Need to Extract:**
1. BVP setup for computing eigenvalues (look at how λ are computed)
2. Operator matrix construction (L̂ discretization on radial grid)
3. Boundary conditions at r=0 and r=R_max

### Step 2: Compute V_bg(r) — Pressure Potential

**Theory:**
```
Background field from three solitons:
  ψ_bg(r) = ψ_e(r - r_e) + ψ_μ(r - r_μ) + ψ_τ(r - r_τ)

where ψ_i are interpolated from Phase 4 profiles:
  r_e = 0
  r_μ ≈ 3.247 (scaled units)
  r_τ ≈ 1.540

Pressure potential (qualitative):
  V_bg(r) ∝ q_i Σ_{j≠i} q_j / (2π r_ij) · [some function of ψ_bg]
```

**Numerical Approach:**
1. Generate three full profiles: ψ_e(r), ψ_μ(r), ψ_τ(r) on 0 ≤ r ≤ 60
2. Superpose shifted: ψ_bg(r) = ψ_e(r) + ψ_μ(|r - 3.247|) + ψ_τ(|r - 1.540|)
3. Estimate V_bg from functional derivative or direct coupling formula
4. Interpolate V_bg(r) to match spectral grid

**Implementation Sketch:**
```python
def compute_pressure_potential(r, r_eq, charges, profiles_interp):
    """
    Compute V_bg(r) at soliton locations from three-soliton configuration
    
    Args:
        r: radial grid (array)
        r_eq: dict with r_eμ, r_eτ, r_μτ (equilibrium distances)
        charges: dict with q_e, q_μ, q_τ
        profiles_interp: dict with interp1d objects for ψ_e, ψ_μ, ψ_τ
    
    Returns:
        V_bg: array of shape (len(r),), pressure potential values
    """
    q_e, q_mu, q_tau = charges['e'], charges['μ'], charges['τ']
    r_eμ, r_eτ, r_μτ = r_eq['eμ'], r_eq['eτ'], r_eq['μτ']
    
    # Coupling strength estimates
    V_eμ = q_e * q_mu / (2 * np.pi * r_eμ)  # e-μ coupling
    V_eτ = q_e * q_tau / (2 * np.pi * r_eτ)  # e-τ coupling
    V_μτ = q_mu * q_tau / (2 * np.pi * r_μτ)  # μ-τ coupling
    
    # For each radial point, estimate local pressure contribution
    # This requires careful geometric calculation of distances
    V_bg = np.zeros_like(r)
    
    for i, r_val in enumerate(r):
        # Distance from this point to each soliton center
        dist_to_e = r_val
        dist_to_μ = np.abs(r_val - r_eμ)
        dist_to_τ = np.abs(r_val - r_eτ)
        
        # Pressure from Goldstone coupling
        # V ~ Σ q_i q_j / (2π r_ij)  [rough estimate]
        V_bg[i] = (V_eμ * np.exp(-dist_to_μ / 2) +     # Decay away from μ
                   V_eτ * np.exp(-dist_to_τ / 2))      # Decay away from τ
    
    return V_bg
```

### Step 3: Modify CP-7 Operator to Include Pressure

**Original Operator (from F-S formulation):**
```
L̂ ψ = -ψ'' - (2/r) ψ' + W''(g) · ψ
```

where W(g) is the potential energy density.

**Modified Operator with Pressure:**
```
L̂_eff ψ = L̂ ψ + V_bg(r) · ψ
        = -ψ'' - (2/r) ψ' + W''(g) · ψ + V_bg(r) · ψ
        = -ψ'' - (2/r) ψ' + [W''(g) + V_bg(r)] · ψ
```

**Implementation in Solver:**
1. Discretize L̂ as before (BVP on radial grid)
2. Add V_bg to diagonal (or off-diagonal if coupling present)
3. Create matrix form: `L_eff_matrix = L_matrix + diag(V_bg)`
4. Diagonalize and extract eigenvalues

### Step 4: Create Phase 4 Extended Script

**File:** `Phase4_extended_spectral_test.py` (new)

**Structure:**
```python
import sys
sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')
from Phase2_bvp_spectrum import soliton_profile, spectrum_on_background

# Section 1: Load baseline (CP-7 results)
# → Print CP-7 λ_min values for e, μ, τ

# Section 2: Generate background configuration
# → Regenerate ψ_e, ψ_μ, ψ_τ at equilibrium r_eq

# Section 3: Compute V_bg(r)
# → For each soliton, superpose pressure from other two

# Section 4: Modify operator and diagonalize
# → L_eff = L_TGP + V_bg for each generation
# → Extract new λ_min values

# Section 5: Compare and analyze
# → Δλ = λ_new - λ_old
# → Assess: stabilization sufficient?

# Section 6: Detailed mode analysis
# → Which modes stabilize? Which remain negative?
# → Track mode energies and spatial profiles

# Section 7: Report results
# → Table: before/after spectral shifts
# → Verdict: does pressure mechanism work?
```

---

## 📊 Success Criteria (for Phase 4 Extended)

### Metric 1: Spectral Shifts
For each soliton, measure:
```
Δλ_i = λ_i^(with_pressure) - λ_i^(isolated)

Desired: Δλ > 0 for saddle modes
Target:  Δλ ≥ |λ_min^(isolated)| for full stabilization
Acceptable: 0 < Δλ < |λ_min| (partial/metastable)
Failure: Δλ ≤ 0
```

### Metric 2: Mode Count
```
Before:  e: N_neg=0,  μ: N_neg=2,  τ: N_neg=3
After:   e: N_neg=?,  μ: N_neg=?,  τ: N_neg=?

Success if: Δ(N_neg) ≥ number_of_saddle_modes
```

### Metric 3: Stability Interpretation
- **Full Stability:** λ_min > 0 for all solitons
- **Metastability:** λ_min < 0 but |λ_min| < scale of potential (tunneling decay >> age of universe)
- **Failure:** No shift or wrong sign

---

## 🗂️ Files & Resources

### Source Code (Read Only)
```
../op-spectral-analysis-Phi-2026-07-03/
├── Phase2_bvp_spectrum.py          ← Import soliton_profile, spectrum_on_background
├── Phase1_discretization_test.py
└── ...
```

### Current Session Results
```
./
├── Phase1b_sympy.py                 (Toy model, completed)
├── Phase2_charge_extraction_v2.py   (Charges: q_e, q_μ, q_τ)
├── Phase3_self_consistency_solver.py (Equilibrium r_eq)
├── Phase4_spectral_stability_test.py (Qualitative, completed)
├── TIER2_PHASE_SUMMARY_2026-07-27.md
└── CHECKPOINT_PHASE4_EXTENDED_2026-07-27.md (This file)
```

### Meta Documentation
```
../meta/
├── TIER2_SESSION64_AUDIT_SUMMARY_2026-07-27.md
└── TIER2_SESSION65_PHASE4_RESULTS_2026-07-27.md
```

---

## ⚠️ Known Issues & Workarounds

### Issue 1: Equilibrium Distances Very Large
```
r_eμ ≈ 32.5 × 10^7 Planck units (scaled to 3.247 for computation)
```
**Workaround:** Use R_SCALE = 1e7 to bring values into numerically stable range

### Issue 2: Three-Soliton Background Overlap
At equilibrium, μ and τ are at different radii. Spatial superposition requires:
```
ψ_bg(r) = ψ_e(r) + ψ_μ(|r - r_μ|) + ψ_τ(|r - r_τ|)
```
**Risk:** If r_μ, r_τ not on same grid as computation points, need interpolation

**Workaround:** Use scipy.interpolate.interp1d with kind='cubic', fill_value=1.0 (vacuum)

### Issue 3: Spectral Solver Expects Single Soliton
CP-7's `spectrum_on_background()` may assume only one soliton + perturbation. For three-soliton system:
**Option A:** Modify to accept three-soliton background  
**Option B:** Apply pressure perturbatively (compute shifts using first-order theory)  
**Option C:** Use numerical eigenvalue solver directly on modified matrix

---

## 📝 Resumption Instructions (If You Stop & Resume)

**To restart Phase 4 Extended from this checkpoint:**

1. **Load this file** to recall current position & constants
2. **Verify Phases 1-3 output** (check TIER2_PHASE_SUMMARY_2026-07-27.md)
3. **Confirm constants** in this checkpoint section
4. **Follow roadmap** (Step 1-7 above)
5. **Use template script** from Step 4
6. **Reference success criteria** (Section "Success Criteria")

**Key files to have ready:**
- `Phase3_self_consistency_solver.py` (for r_eq, E_press values)
- `Phase2_charge_extraction_v2.py` (for q values)
- `../op-spectral-analysis-Phi-2026-07-03/Phase2_bvp_spectrum.py` (spectral solver)

**Timeline Estimate:**
- Step 1 (Understand CP-7): 1-2 hours (reading code)
- Step 2 (Compute V_bg): 3-4 hours (implementation + testing)
- Step 3 (Modify operator): 4-6 hours (BVP modification)
- Step 4 (Create script): 2-3 hours (integration)
- Step 5-7 (Testing & reporting): 2-3 hours (debugging, output)
- **Total:** ~15-20 hours over 2-3 working days

---

## 🎯 Expected Outcomes (After Phase 4 Extended)

### Outcome A: Full Stabilization ✅
```
μ: Δλ ≥ 1.28  (N_neg: 2 → 0)
τ: Δλ ≥ 4.22  (N_neg: 3 → 0)
Result: Native pressure alone is sufficient
→ Path N4d: SUCCESS
→ Tier 2: COMPLETE
→ Publication: "Self-Consistent Pressure Stabilizes Lepton Solitons"
```

### Outcome B: Partial Stabilization ⚠️
```
μ: 0 < Δλ < 1.28  (N_neg: 2 → 1)
τ: 0 < Δλ < 4.22  (N_neg: 3 → 1-2)
Result: Pressure helps but insufficient
→ Path: N4d + N4c hybrid (pressure + radiative corrections)
→ Tier 2: EXTENDED (2-3 more weeks)
→ Publication: "Pressure and Radiative Effects in Lepton Stabilization"
```

### Outcome C: No Stabilization ❌
```
Δλ ≤ 0 or minimal shift
Result: Pressure mechanism doesn't work
→ Path: Investigate N4c (radiative), N4b (symmetry), N4a (substrate)
→ Tier 2: RESTART (pick alternative path)
→ Publication: "Pressure mechanism limitations; exploring alternatives"
```

---

## 💾 Backup & Archiving

**Current session files are saved at:**
```
TGP/TGP_v1/research/op-native-pressure-lepton-stability-2026-07-27/
```

**Checkpoint marker:** 2026-07-27 end-of-session  
**Next session should:** Read this file first, then proceed with Step 1 of roadmap

---

## 📞 Quick Reference

**Constants (Copy-Paste Ready):**
```python
G0_E, G0_MU, G0_TAU = 1.24915, 2.02117, 3.18912
Q_E, Q_MU, Q_TAU = 1.20009652, 1.10685470, 1.04924277
R_EQ = {'eμ': 3.247264, 'eτ': 1.540338, 'μτ': 2.130247}
E_PRESS = -5.046
```

**Baseline Spectra (Copy-Paste Ready):**
```python
BASELINE = {
    'e': {'l0': {'λ': -0.998, 'N_neg': 0}},
    'μ': {'l0': {'λ': -1.282, 'N_neg': 2}},
    'τ': {'l0': {'λ': -4.216, 'N_neg': 3}},
}
```

---

**Created:** 2026-07-27  
**Status:** Ready for Phase 4 Extended  
**Author:** Claudian (AI Assistant)  
**Next Step:** Implement Phase 4 Extended following roadmap above
