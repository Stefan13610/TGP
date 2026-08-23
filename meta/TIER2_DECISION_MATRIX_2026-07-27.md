---
title: "TIER 2 Decision Matrix — Session #64 Synthesis"
date: 2026-07-27
type: decision_matrix
status: READY_FOR_USER_DECISION
owner: Claudian
session: "#64 2026-07-27"
---

# TIER 2 Decision Matrix: Where Do We Go From Here?

## Executive Summary

CP-7 found that μ and τ solitons are **saddle points** (index 2/3) in the F-S formulation, with tachyonic continuum starting from σ = −1. Three stabilization paths (linear constraints, soft wall, U(1) charge) have failed. **Four forward options remain, each with distinct tradeoffs.**

---

## Context Recap (60 seconds)

| Session | Cycle | Finding | Verdict |
|---------|-------|---------|---------|
| #60 | CP-7 | Solitonic spectrum: tachyons + saddle points | Unexpected, structural |
| #62 | op-wall-dynamics | Linear constraints cannot stabilize | Failed (W1–W3) |
| #63 | op-nonlinear-charge | U(1) charge quantization fails | Failed (V1–V3) |
| **#64** | **Tier 2 Analysis** | **Four forward paths identified** | **Ready to choose** |

---

## The Four Paths: Quick Comparison

### 🟢 **N4d: Metastability (FASTEST, LOWEST RISK)**
**Ask:** Are μ/τ long-lived via tunneling?
- **What changes:** Interpretation only (accept saddle structure as metastable)
- **Timeline:** 1 week
- **Risk:** Medium (WKB accuracy, barrier height uncertainty)
- **Upside:** No axiom change; no mass ratio recalculation; fast answer
- **Downside:** Shifts narrative from "stable" → "metastable"; requires WKB verification

**✅ Recommendation:** **START HERE**

---

### 🟡 **N4c: Radiative Corrections (MEDIUM RISK)**
**Ask:** Do F-A loop corrections stabilize F-S vacuum?
- **What changes:** Effective potential V_eff gets radiative corrections
- **Timeline:** 2.5 weeks
- **Risk:** High (scheme-dependent, technically demanding)
- **Upside:** Addresses F-A/F-S coupling; principled QFT approach
- **Downside:** Convergence uncertain; might require multi-loop; renormalization scale choice

**✅ Recommendation:** **Pursue if N4d says "not metastable"**

---

### 🟡 **N4b: Symmetry Extension (MEDIUM-HIGH RISK)**
**Ask:** Does Z₂ → larger symmetry provide stabilizing charges?
- **What changes:** Gauge sector expanded (new symmetry, potentially new bosons)
- **Timeline:** 3.5 weeks
- **Risk:** High (multi-charge Q-ball complexity; phenomenology risk from new bosons)
- **Upside:** 3-fold structure hints at SU(2) or Z₃; natural for generation problem
- **Downside:** Many choices; risk of new particles conflicting with experiments

**✅ Recommendation:** **Pursue if N4c fails, or in parallel with N4c**

---

### 🔴 **N4a: Discrete Substrate (SLOWEST, HIGHEST RISK)**
**Ask:** Does lattice discretization eliminate tachyonic modes?
- **What changes:** Fundamental axiom (continuum → discrete substrate)
- **Timeline:** 3 weeks
- **Risk:** Very High (axiom change, mass ratio recalculation, computational load)
- **Upside:** Addresses tachyons at source; might recover Brillouin zone structure
- **Downside:** Entire TGP built on continuum; PPN/GW predictions could change; mass ratios uncertain

**✅ Recommendation:** **Last resort** — only if N4b/c both fail

---

## Decision Tree

```
Are you willing to accept metastability?
│
├─ YES (tunneling lifetime >> 13.8 Gyr)
│   └─ ✅ CONTINUE TIER 2 with metastable interpretation
│       (Run N4d for WKB verification; if viable, update Limitations paper)
│       
├─ MAYBE (need to check whether tunneling is slow)
│   └─ ✅ RUN N4d FIRST (1 week), then decide based on τ_decay
│       If τ > 10^{20} s → metastability OK → continue
│       If τ < 10^{15} s → metastability fails → escalate
│       
└─ NO (must fix saddle structure)
    ├─ TRY N4c (radiative corrections)
    │   If V_eff''(1) > 0 → mass ratio check → proceed
    │   If V_eff''(1) < 0 → proceed to N4b
    │
    ├─ TRY N4b (symmetry extension)
    │   If multi-charge Q-ball stabilizes → extend core → proceed
    │   If fails → proceed to N4a (if desperate)
    │
    └─ TRY N4a (discrete substrate, last resort)
        If tachyonic band vanishes → recalculate masses → proceed
        If fails → **restructure TGP** (beyond Tier 2)
```

---

## What Each Path Means for the Lepton Paper

### If N4d (Metastability) Succeeds
**Limitations addition:**
> "The μ and τ solitons are saddle points in the static variational sense, but may be long-lived via tunneling suppression (Appendix: WKB lifetime estimate). If tunneling timescale >> 13.8 Gyr, these configurations are effectively stable for cosmological purposes."

**Mass ratios:** Unchanged
**Core revision:** Additive only (remark in sek08b on metastability)

---

### If N4c (Radiative Corrections) Succeeds
**Limitations addition:**
> "The tree-level solitonic potential has tachyonic continuum. Including one-loop F-A (gravitational sector) corrections stabilizes the vacuum. Renormalized effective potential has V_eff''(1) > 0, restoring spectral stability of μ/τ."

**Mass ratios:** Recalculate (might shift by ~0.1–1%)
**Core revision:** Update sek08c with V_eff structure; additive remarks in sek08b

---

### If N4b (Symmetry Extension) Succeeds
**Limitations removal:**
> "T-OP4 fully resolved. Extended Z₂ → Z₂×U(1)×... symmetry with multi-charge quantization stabilizes the solitonic sector. Fixed-charge formulation recovers spectral stability."

**Mass ratios:** Recalculate (depends on new charges)
**Core revision:** Major — Lagrangian expanded, gauge sector modified, core sek08b rewritten

---

### If N4a (Discrete Substrate) Succeeds
**Implications:**
> "TGP is fundamentally discrete at the Planck scale; continuum limit emerges. Tachyonic modes are continuum artifact. Lattice spectrum converges to stable band structure by N → ∞."

**Mass ratios:** Recalculate (lattice solitons different from continuum)
**Core revision:** Catastrophic — axioms §1 rewritten, substrate description §0 revised, all energy calculations changed

**Realistic assessment:** If N4a becomes necessary, paper goes back to draft stage (2–3 months revision).

---

## User Decision Point: Three Options

### Option A: Optimize for Speed 🚀
**Execute order:** N4d → N4c → N4b → (if needed: N4a)
- **Timeline to decision:** 4–5 weeks (each phase ~1 week)
- **Risk:** Might miss optimal path if N4d says "not metastable" but N4c is viable
- **Outcome:** Fastest path to working solution

### Option B: Parallel Exploration ⚡
**Execute in parallel:** N4d + N4c (overlapping teams)
- **Timeline to decision:** 3 weeks (both done by day 21)
- **Risk:** Requires splitting focus; computational resources
- **Outcome:** Faster decision on N4d vs N4c; can pivot faster

### Option C: Deep Dive into Metastability 🔍
**Execute only:** N4d (comprehensive study, multi-barrier analysis)
- **Timeline to decision:** 2 weeks (detailed WKB for all decay channels)
- **Risk:** If metastability fails, wasted time; must then escalate to N4b/c anyway
- **Outcome:** Authoritative answer on whether metastability is viable

---

## My Recommendation

**Start with Option B (Parallel) — Day 1 onwards:**

**Path 1 (Fast, this session):**
- Begin **N4d** (metastability assessment)
- Identify critical points in potential landscape
- Estimate WKB exponents for μ/τ decay paths

**Path 2 (Medium, weeks 2–3):**
- Prepare **N4c** (radiative corrections)
- Set up Feynman diagrams; plan loop integrals
- Parallel with N4d completion

**Milestone (Day 7–14):**
- N4d results: "metastability viable" or "fails"
- If viable → Continue Tier 2 with metastable interpretation ✅
- If fails → Escalate to N4c/b

**Rationale:**
- N4d is fastest, lowest-risk, quick answer
- Can start N4c prep while waiting for N4d results (no conflict)
- If N4d succeeds, time saved; if not, N4c ready to start immediately
- Avoids "wasted week" if N4d is decision-killing

---

## Realistic Outcomes: Three Scenarios

### 🟢 Scenario 1: Metastability Sufficient (Probability: ~40%)
**What happens:**
- N4d: τ_decay > 10^{20} s ✓
- Metastable interpretation accepted
- Limitations paper updated with one paragraph
- **Tier 2 COMPLETE** (day 7)
- Proceed to publications

**Impact:** TGP narrative shifts slightly, but no deep changes

---

### 🟡 Scenario 2: Radiative Corrections Save the Day (Probability: ~30%)
**What happens:**
- N4d: τ_decay too short ✗
- N4c: V_eff''(1) > 0 after loop corrections ✓
- Mass ratios recalculated, r₂₁/r₃₁ preserved within 0.1% ✓
- Core updated with V_eff structure
- **Tier 2 COMPLETE** (day 18)
- Lepton paper: new section on radiative stabilization

**Impact:** TGP grounded in more complete QFT; more rigorous foundation

---

### 🔴 Scenario 3: All Paths Fail (Probability: ~30%)
**What happens:**
- N4d: Metastability fails ✗
- N4c: Radiative corrections insufficient ✗
- N4b: Multi-charge Q-ball doesn't work ✗
- N4a: Discrete substrate loses mass ratios ✗
- **Tier 2 FAILURE** (day 40–50)

**Next steps:**
- Admit CP-7 saddle structure is fundamental incompatibility
- Enter "restructuring phase" (6–12 weeks):
  - Reconsider Z₂ axiom?
  - Add new degree of freedom (scalar, spinor, metric)?
  - Revisit metric ansatz (new f(ψ))?

**Impact:** TGP needs substantial revision; publications delayed; new research cycle

---

## Checkpoints & Go/No-Go Criteria

| Checkpoint | Day | Go/No-Go Criterion |
|-----------|-----|-------------------|
| **N4d preliminaries** | 3 | Potential landscape mapped; decay paths identified → **GO** |
| **N4d WKB done** | 7 | τ_decay estimated; clear answer → **GO or ESCALATE** |
| **N4c setup** | 10 | Feynman diagrams ready; loop integral strategy clear → **GO** |
| **N4c result** | 18 | V_eff''(1) sign determined; clear verdict → **GO or NEXT** |
| **N4b/a decision** | 21 | Based on N4c: pursue N4b or N4a? → **USER DECISION** |
| **Tier 2 complete** | 50 | At least one path succeeded and mass ratios preserved → **TIER 2 END** |

---

## Questions for You (User)

Before I recommend which path to pursue first, please clarify:

### Q1: Philosophical stance on metastability
> "Would TGP work acceptably if μ/τ are metastable (lifetimes >> 13.8 Gyr) rather than stable?"

- [ ] Yes, acceptable (pursue N4d first)
- [ ] Maybe, need to see numbers (pursue N4d first)
- [ ] No, must have true stability (skip N4d, go to N4c/b)

### Q2: Tolerance for core axiom changes
> "Are you open to revising the TGP axioms (e.g., discrete substrate, extended symmetry) if N4c fails?"

- [ ] Yes, explore fully (N4a/b viable)
- [ ] Reluctantly, if necessary (N4b before N4a)
- [ ] No, must preserve core axioms (only N4c/d)

### Q3: Timeline / Resource constraints
> "How much time/computational resource can we spend on Tier 2?"

- [ ] Full speed, whatever it takes (Option B: parallel)
- [ ] ~4–5 weeks max (Option A: sequential)
- [ ] ~1–2 weeks max (Option C: metastability deep dive)

### Q4: Publication timeline
> "When do you want the lepton paper finalized?"

- [ ] ASAP (days, depends on Tier 2 outcome)
- [ ] In 4–6 weeks (can wait for N4 path)
- [ ] In 2–3 months (can do major revisions if needed)

---

## What I'm Ready to Do

**Immediately (Session #64 ongoing):**
- ✅ Create detailed N4d WKB framework
- ✅ Write pseudocode for metastability search (identify critical points, decay paths)
- ✅ Estimate computational cost for parallel runs

**Week 1:**
- ✅ Run N4d (WKB tunneling analysis)
- ✅ Generate preliminary answer: "metastability viable?" Y/N

**Week 2–3:**
- ✅ Based on N4d result, launch N4c or N4b as appropriate
- ✅ Generate mass-ratio recalculations if new physics found

**Week 4–5:**
- ✅ Synthesize Tier 2 results; prepare final decision
- ✅ Update lepton paper with new findings (Limitations + methods)

---

## Cross-References

- [[./TIER2_ANALYSIS_FRAMEWORK_2026-07-27.md]] — Complete framework
- [[./CP7_ARTIFACT_ANALYSIS_2026-07-27.md]] — Artifact analysis
- [[./N4_PATHS_FEASIBILITY_ASSESSMENT_2026-07-27.md]] — Path details
- [[../research/op-spectral-analysis-Phi-2026-07-03/README.md]] — CP-7 findings
- [[../paper_lepton_masses/tgp_lepton_masses.tex]] — Target paper for updates

---

## Next Step: Your Input

**Please provide:**
1. Answers to Q1–Q4 above (or similar guidance)
2. Any additional constraints or preferences not listed
3. Go-ahead signal to start N4d (metastability assessment)

**Then I will:**
- Launch N4d pipeline immediately
- Report preliminary results (day 3–7)
- Iterate based on findings

