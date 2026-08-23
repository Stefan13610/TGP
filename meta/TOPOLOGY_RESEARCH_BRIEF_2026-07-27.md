# 🧬 Research Brief: Topological Structure in CP-7 Profiles

**Purpose:** Investigate whether topological features (winding numbers, topological degree) explain saddle points better than energy/symmetry arguments

**Parallel to:** Phase 4 Extended (numerical spectral test)  
**Timeline:** Independent, self-paced  
**Status:** OPEN RESEARCH

---

## 🎯 Core Questions

1. **Do CP-7 profiles encode winding numbers?**
   - Each soliton g_i(r) can be analyzed for topological degree
   - Is there an integer-valued quantity that's conserved?
   - Can we extract W_i (winding number) from g_i just like we extracted charges q_i?

2. **Is hierarchy topological?**
   - Does e have W_e = 1?
   - Does μ require W_μ = something different in presence of e?
   - Does τ configuration depend on (ψ_e + ψ_μ)?

3. **Are saddle points signature of topological constraint?**
   - Saddle in μ: means μ "wants to twist" to different configuration
   - Saddle in τ: means τ "wants to change winding"
   - But in hierarchy, winding is fixed → saddle becomes irrelevant?

4. **Charge vs Winding:**
   - Is q_i really just ∂W/∂g (derivative of winding)?
   - Or are they independent?

---

## 📊 Analysis Steps (Concrete)

### Step 1: Extract Winding Numbers from Profiles

**From Phase 2/3 code you already have:**
```python
# In Phase2_charge_extraction_v2.py or Phase3_self_consistency_solver.py
# You have: r, g(r) for each soliton (e, μ, τ)

# NEW: Compute winding number
def compute_winding_number(r, g, g_prime):
    """
    For each soliton, compute topological degree:
    
    W = (1/π) ∫ d(arctan(g'/g)) dr  [if g is phase-like]
    
    Or more careful:
    W = (1/2π) ∮ dθ/dφ  [on circle at infinity]
    
    For radial profiles, this becomes:
    W = phase_shift_at_infinity / π
    """
    
    # Compute phase from profile
    # Use CP-7's phase_shift_l0() concept
    # Extract: integer W
    pass
```

**What to look for:**
- W_e ≈ 1? (reference)
- W_μ ≠ 1? (different configuration)
- W_τ ≠ W_μ? (third level)

### Step 2: Check If Charges Correlate with Winding

```python
# Compare from Phase 2:
# q_e = 1.20009652
# q_μ = 1.10685470
# q_τ = 1.04924277

# With new:
# W_e = ?
# W_μ = ?
# W_τ = ?

# Hypothesis:
# q_i = f(W_i, g₀_i)  [charge depends on winding]
# 
# If true: charge is NOT fundamental, winding is
```

### Step 3: Analyze Eigenvectors in CP-7

**From CP-7 Phase2_bvp_spectrum.py output:**
```python
# CP-7 computes eigenvalues λ_i for each soliton
# But also eigenvectors v_i(r)

# NEW: For saddle modes (λ < 0), check:
# - Do they preserve topological structure?
# - Do they flip between g_+ and g_- components?
# - Is there a protected topological mode?

def analyze_saddle_mode_topology(r, v_saddle, g_ref):
    """
    Saddle mode v(r) = perturbation to g(r)
    
    Question: Does v(r) try to change the winding?
    
    If yes: saddle exists because winding is protected
            (in isolation, free to unwind)
    If no: saddle has different origin
    
    Test: Compute winding of (g + ε·v)
          Does it stay integer? Or try to branch cut?
    """
    pass
```

### Step 4: Check Hierarchy Hypothesis

**Prediction:**
```
If hierarchy is topological:
  
  μ alone: W_μ = ? (saddle points, λ < 0)
  μ in e's field: W_μ = ? (different? protected? stable?)
  
  τ alone: W_τ = ? (many saddle points, λ < 0)
  τ in (e+μ) field: W_τ = ? (changes?)
```

**Test:**
- Regenerate CP-7 spectra but with "μ ON ψ_e background" (not isolated)
- Does W_μ stay same or change?
- Do saddle modes persist or vanish?

### Step 5: Look for Chirality Structure

**Spinor Hypothesis:**
```
If g(r) = (g_+(r), g_-(r)), then:

q_i = ∫ (g_+(r) - g_-(r)) dr  [chirality imbalance]

or

W_i = (winding of g_+) - (winding of g_-)  [chiral winding difference]
```

**How to test:**
- Can you decompose CP-7 profiles into ± components?
- Do they have different winding?
- Does the difference correlate with observed charge?

---

## 🔬 Concrete Deliverables

### If You Find: Winding Numbers in Profiles

**Output:**
```
W_e = [integer]
W_μ = [integer]
W_τ = [integer]

Correlation with q_i:
  q_e = f(W_e, g₀_e)  → reveals function f
  q_μ = f(W_μ, g₀_μ)
  q_τ = f(W_τ, g₀_τ)
```

### If You Find: Hierarchy Dependence

**Output:**
```
Spectrum(μ alone)           vs  Spectrum(μ in ψ_e)
λ_sad: -1.282, N_neg: 2        λ_sad: ?, N_neg: ?

Spectrum(τ alone)           vs  Spectrum(τ in ψ_e+ψ_μ)
λ_sad: -4.216, N_neg: 3        λ_sad: ?, N_neg: ?

Interpretation:
  If saddle vanishes in hierarchy → it's topological artifact
  If saddle persists → independent instability
```

### If You Find: Topological Protection

**Output:**
```
Eigenvector v_i of saddle mode:
  - Tries to unwind? (flips topological index)
  - Or constrained? (protected by topology)

Implication:
  Protected → winding is fundamental
  Not protected → saddle is "free to move"
```

---

## 🧠 Intellectual Roadmap

### Week 1: Exploration
- Extract winding from profiles (if possible)
- Check correlation with known charges
- Look at eigenvector structure

### Week 2: Testing Hypothesis
- Modify CP-7 to compute "μ on ψ_e background"
- Compare spectra: isolated vs hierarchical
- Check if saddle structure changes

### Week 3: Synthesis
- If winding explains saddle → topological interpretation
- If hierarchy explains saddle → energetic interpretation changes
- Compare with Phase 4 Extended results

---

## 📚 Reference Material (in Vault)

**CP-7 files:**
- `op-spectral-analysis-Phi-2026-07-03/Phase2_bvp_spectrum.py`
  → Look for phase_shift_l0() function (already computes phases!)
  → Can extract winding from this

**Your own:**
- `op-native-pressure-lepton-stability-2026-07-27/Phase2_charge_extraction_v2.py`
  → Already extracts A_tail (far-field amplitude)
  → Can add winding extraction next to it

- `op-native-pressure-lepton-stability-2026-07-27/Phase3_self_consistency_solver.py`
  → Has profiles loaded
  → Can instrument to extract topological properties

---

## 💻 Minimal Working Example

```python
import numpy as np
from scipy.interpolate import interp1d

# From your Phase 2/3 code:
r, g, gp, d2 = profile_data  # (r, g(r), g'(r), g''(r))

# Compute phase angle (if g encodes phase information)
# or try: θ(r) = arctan(some combination of g, g')

# Extract winding by phase integral:
dtheta = np.diff(np.arctan2(gp, g))  # or similar
W = np.sum(dtheta) / (2 * np.pi)  # Total winding

print(f"Winding number W ≈ {W:.2f}")

# Compare with charge:
A_tail = 1.2  # (from Phase 2)
q = A_tail / (g[0] - 1)
print(f"Charge q ≈ {q:.4f}")
print(f"Ratio q/W = {q/W:.4f}  [if related, ratio should be meaningful]")
```

---

## 🎯 Success Criteria

### For This Research:

**You've Found Something Real If:**
1. Winding numbers are extractable and non-trivial (not all 0 or 1)
2. Charges correlate with winding in interpretable way
3. Hierarchy affects winding in predictable pattern
4. Saddle eigenvectors show topological structure

**You've Found a Blind Spot If:**
1. Phase 4 Extended shows Δλ > 0 (pressure works)
2. BUT winding analysis shows it's actually topological protection
3. → Then pressure is red herring, mechanism is deeper

---

## ⚡ Key Insight to Keep

> "Ładunek to konfiguracja przestrzeni"

If you can show:
```
Charge q_i = [topological descriptor of g_i]
           = winding / chiral imbalance / something integer-like
```

Then Session #63 failed for reason: **looked for wrong conserved quantity**
- Searched: U(1) phase symmetry
- Should have: topological winding number

And Session #65 might succeed for **different reason** than we think:
- We thought: pressure stabilizes
- Might be: pressure enforces hierarchy which protects topology

---

## 📝 Report Back When You Find

- Winding numbers (or "why they don't exist")
- Hierarchy effects on spectrum
- Topological structure in eigenvectors
- Correlation with charges

**Your findings + Phase 4 Extended results = Complete Picture** 

---

**Research Type:** Exploratory  
**Urgency:** None (parallel to Phase 4)  
**Direction:** Topological aspects of TGP lepton sector  
**Status:** Ready to start

Go dig. 🔍
