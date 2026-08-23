# 🚀 Phase 4 Extended: Quick Start

**Status:** Ready to begin implementation  
**Timeline:** 1-2 weeks  
**Parallel:** User research on topology  
**Date:** 2026-07-27

---

## 📍 Where You Are

✅ Phases 1-4 (Qualitative) complete
✅ All constants/equilibrium known
✅ Roadmap documented in CHECKPOINT_PHASE4_EXTENDED_2026-07-27.md
⏳ Ready to build L̂_eff and diagonalize

---

## 🎯 What Phase 4 Extended Does (1-Liner)

**Take CP-7's spectral solver + add pressure term V_bg(r) + diagonalize + measure spectral shifts Δλ**

---

## 🔴 START HERE: Step-by-Step

### STEP 0: Setup (Now)

```bash
# Already have:
# - Phase2_charge_extraction_v2.py (charges: q_e, q_μ, q_τ)
# - Phase3_self_consistency_solver.py (equilibrium: r_eq, E_press)
# - CP-7 at: ../op-spectral-analysis-Phi-2026-07-03/Phase2_bvp_spectrum.py

# Verify CP-7 can import:
python3 -c "import sys; sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03'); from Phase2_bvp_spectrum import soliton_profile; print('✓ CP-7 imports OK')"
```

### STEP 1: Understand CP-7 Operator (2-3 hours)

**Read:** `../op-spectral-analysis-Phi-2026-07-03/Phase2_bvp_spectrum.py`

**Find and understand:**
1. Function `spectrum_on_background()` - how does it compute eigenvalues?
2. How is L̂ discretized? (BVP matrix form)
3. What are boundary conditions?
4. How are eigenvalues extracted? (eigsh? eigh?)

**Key lines to locate:**
```python
# Look for:
# - Operator matrix construction (L_matrix = ...)
# - Diagonalization call (eigh(...) or eigsh(...))
# - Eigenvalue extraction (eigenvals = ...)
# - Boundary condition setup
```

**Deliverable:** Understand how CP-7 constructs and diagonalizes L̂

### STEP 2: Write V_bg Function (3-4 hours)

**Template from:** `Phase4_extended_spectral_test_TEMPLATE.py` (Section 3)

**Goal:** Function that computes V_bg(r) on spectral grid

**Implementation:**
```python
def compute_pressure_potential(r_grid, r_eq, charges, profiles_interp):
    """
    Create V_bg(r) from three-soliton pressure configuration
    
    Inputs:
      r_grid: radial grid from CP-7 (array, r ∈ [0, R_max])
      r_eq: equilibrium distances {eμ, eτ, μτ}
      charges: {q_e, q_μ, q_τ}
      profiles_interp: interp1d objects for ψ_e, ψ_μ, ψ_τ
    
    Output:
      V_bg: array, pressure potential on r_grid
    
    Physics:
      V_bg ~ Σ q_i q_j exp(-|r - r_ij| / scale)
      or smooth falloff from soliton centers
    """
    
    # Start simple: exponential decay from each center
    # Refine if needed
    pass
```

**Test:** `V_bg_test = compute_pressure_potential(test_r, R_EQ, Q, interp_g)`  
Check: V_bg ~ 10^-4 magnitude (like Phase 4 estimate)

### STEP 3: Modify CP-7 to Add V_bg (4-6 hours)

**Goal:** Create modified operator L̂_eff = L̂_TGP + V_bg

**Approach 1 (Simplest):**
```python
# Get L from CP-7
L_matrix, r_grid = extract_cp7_operator(g0, formulation, R, N)

# Add pressure
V_bg = compute_pressure_potential(r_grid, ...)
L_eff = L_matrix + np.diag(V_bg)

# Diagonalize
eigenvals, eigenvecs = eigh(L_eff)
```

**Approach 2 (If CP-7 uses sparse matrix):**
```python
# If L is sparse, convert/manipulate carefully
L_eff = L_matrix + scipy.sparse.diags(V_bg)
eigenvals = scipy.sparse.linalg.eigsh(L_eff, k=20, which='SM')[0]
```

**Deliverable:** Function that creates L̂_eff and returns eigenvalues

### STEP 4: Write Main Phase 4 Extended Script (2-3 hours)

**Template:** Use `Phase4_extended_spectral_test_TEMPLATE.py`

**Fill in:**
- Section 4: Diagonalization loop (call your L_eff function)
- Section 5: Comparison table (λ_old vs λ_new, compute Δλ)
- Section 6: Mode analysis (which modes shift?)
- Section 7: Verdict (does it stabilize?)

**Minimal version:**
```python
#!/usr/bin/env python3
import numpy as np
import sys
sys.path.insert(0, '../op-spectral-analysis-Phi-2026-07-03')
from Phase2_bvp_spectrum import soliton_profile

# Constants
Q_E, Q_MU, Q_TAU = 1.20009652, 1.10685470, 1.04924277
R_EQ_EMU, R_EQ_ETAU, R_EQ_MUTAU = 3.247264, 1.540338, 2.130247
BASELINE = {'e': -0.998, 'μ': -1.282, 'τ': -4.216}

# For each soliton
for gen in ['e', 'μ', 'τ']:
    g0 = {'e': 1.24915, 'μ': 2.02117, 'τ': 3.18912}[gen]
    
    # Get profiles
    r, g, gp, d2, *_ = soliton_profile('F-S', g0, R=60.0, N=4000)
    
    # Compute V_bg (from all three solitons)
    V_bg = compute_pressure_potential(r, r_eq, charges, ...)
    
    # Get CP-7 operator
    L_matrix = extract_cp7_operator(gen, g0, r)
    
    # Add pressure
    L_eff = L_matrix + np.diag(V_bg)
    
    # Diagonalize
    eigenvals = np.linalg.eigvalsh(L_eff)
    
    # Compare
    lambda_new = eigenvals[0]
    lambda_old = BASELINE[gen]
    delta_lambda = lambda_new - lambda_old
    
    print(f"{gen}: λ_old={lambda_old:7.3f} → λ_new={lambda_new:7.3f}, Δλ={delta_lambda:7.3f}")
```

### STEP 5: Debug & Refine (2-3 hours)

**Common Issues:**
1. V_bg magnitude wrong → adjust coupling constant
2. Matrix dimensions mismatch → check grid sizes
3. Eigenvalue not real → L_eff not symmetric → check V_bg
4. No change in λ → V_bg too weak → increase magnitude

**Testing:**
- Run on single soliton first (e)
- Verify Δλ ≠ 0
- Then extend to all three

### STEP 6: Produce Results (1-2 hours)

**Output Table:**

```
PHASE 4 EXTENDED RESULTS
═══════════════════════════════════════════════════════════

Soliton e (g₀=1.24915):
  Isolated:        λ_min = -0.998   N_neg = 0
  With pressure:   λ_min = ?        N_neg = ?
  Δλ = ?           Status: ?

Soliton μ (g₀=2.02117):
  Isolated:        λ_min = -1.282   N_neg = 2
  With pressure:   λ_min = ?        N_neg = ?
  Δλ = ?           Status: ? SADDLE POINT SOLVED?

Soliton τ (g₀=3.18912):
  Isolated:        λ_min = -4.216   N_neg = 3
  With pressure:   λ_min = ?        N_neg = ?
  Δλ = ?           Status: ? SADDLE POINT SOLVED?

VERDICT:
  [ ] Full Stabilization (Δλ ≥ |λ_min| for μ, τ)
  [ ] Partial Stabilization (0 < Δλ < |λ_min|)
  [ ] No Stabilization (Δλ ≤ 0)
  [ ] Unexpected (something else)
```

### STEP 7: Analysis & Interpretation (1-2 hours)

**Interpret based on result:**

**If Full Stabilization:**
```
✓ Native pressure alone is sufficient
✓ N4d: SUCCESS
✓ Hierarchy hypothesis: Not needed
→ Report findings, prepare publication
```

**If Partial Stabilization:**
```
⚠ Pressure helps but insufficient
→ Combine with N4c (radiative corrections)
→ Check if user's topology research offers alternative
→ Δλ_total = Δλ_pressure + Δλ_loops?
```

**If No Stabilization:**
```
❌ Pressure doesn't help (surprising!)
→ But mechanism is sound (sign is right)
→ Problem: magnitude, not concept
→ Pivot to N4c or investigate topology deeper
```

---

## 📋 Checklist

- [ ] STEP 0: CP-7 imports successfully
- [ ] STEP 1: Understand CP-7 operator (document findings)
- [ ] STEP 2: V_bg function works, correct magnitude
- [ ] STEP 3: L_eff creation works, eigenvalues real
- [ ] STEP 4: Main script runs without errors
- [ ] STEP 5: Results table populated
- [ ] STEP 6: Interpretation written
- [ ] STEP 7: Comparison with user's topology research (when ready)

---

## 🎲 Expected Results (For Reference)

**Optimistic (60% chance):**
```
Δλ_μ ~ 1.0 to 1.3    (saddle mode shifts to positive)
Δλ_τ ~ 3.0 to 4.5    (saddle modes shift to positive)
→ Full stabilization
```

**Pessimistic (25% chance):**
```
Δλ_μ ~ 0.1 to 0.5    (shifts but insufficient)
Δλ_τ ~ 0.5 to 2.0    (helps but doesn't fully stabilize)
→ Partial; need N4c
```

**Surprise (15% chance):**
```
Δλ ≤ 0 or very small
→ Mechanism doesn't work numerically (but math is sound)
→ Or topology is more important than we thought (user's research)
```

---

## 🔗 Integration with User's Topology Research

**When you (user) find winding numbers:**
```
Compare:
  W_i (winding) vs Δλ_i (spectral shift)
  
If they correlate:
  → Topology IS the real mechanism
  → Phase 4 Δλ is consequence, not cause
  → Mechanism is hierarchical+topological, not pressure
  
If they don't:
  → Pressure is real mechanism
  → Topology is independent feature
  → Both true: different levels of description
```

---

## 📞 Questions While Working

**If stuck:**
- Check CHECKPOINT_PHASE4_EXTENDED_2026-07-27.md (detailed roadmap)
- Check Phase4_extended_spectral_test_TEMPLATE.py (pseudocode)
- Check TIER2_PHASE_SUMMARY_2026-07-27.md (theory background)

---

## 🚀 GO!

Start with STEP 0 (verify imports), then 1-2 hours on STEP 1 (understand CP-7).

Report back when you hit first issue or first success.

Good luck! 🎯
