---
title: "Phase 3 plan — H_0 numerical estimate z derived equations (F-γ-3 PRIMARY KILLER test)"
type: phase_plan
status: PRE_REGISTERED_LOCKED
phase: 3
parent_cycle: op-CE-H-gamma-3-cosmological-2026-05-23
pre_registration_date: 2026-05-23
---

# Phase 3 plan — H_0 numerical estimate z derived equations (F-γ-3 PRIMARY KILLER test)

**Status:** PRE_REGISTERED_LOCKED 2026-05-23.
**Pre-registration:** **PRE-NUMERICAL** — wszystkie thresholds + interpretations LOCKED przed sympy/numerical execution.
**Anti-Lakatos LOCK CRITICAL:** Po Phase 2 result (geometric H_0 = 1/t_universe), prediction wygląda na FAVORABLE. NIE wolno modyfikować F-γ-3 thresholds — bez względu na sympatie.

---

## §1 — Pre-registered FP (Fixed-point predictions)

**0 hardcoded T_pass=True.** Każdy FP = compute-then-compare against PRE-REGISTERED threshold.

### T_P3_1 — R(t) = c·t consistency (frontier dynamics)

**Hypothesis:** Per concept paper §2.4 + §3.3, frontier moves at c (relativistic Phi-substrate wave equation z (D3); maximum signal speed). Therefore R(t) = c · t + R_0; w late-time limit R_0 << c·t → R ≈ c·t.

**Pre-registered (literal):**
- Frontier velocity = c (consistency check z wave equation)
- R(t) = c·t late-time limit (R_0 negligible)
- PASS: derivation consistent z relativistic dispersion
- FAIL: NIE consistent

### T_P3_2 — H(t) = 1/t derived z R(t) = c·t

**Hypothesis:** H(t) = (1/R)(dR/dt) = c/(c·t) = 1/t

**Pre-registered (literal):**
- Symbolic derivation H(t) = 1/t
- PASS: sympy verification
- FAIL: NIE derive

### T_P3_3 — F-γ-3 PRIMARY KILLER verdict (numerical)

**Pre-registered threshold (LOCKED 2026-05-21; activated 2026-05-23):**
- TARGET: H_0 ∈ [67, 73] km/s/Mpc (observed Planck + SH0ES)
- TOLERANCE: H_0 ∈ [33.5, 146] km/s/Mpc (factor 2 anti-cherry-pick)
- FAIL trigger: H_0 outside [33.5, 146] → HALT-B

**Independent t_universe anchor (anti-circularity):**
- Stellar ages (globular clusters): t_universe ∈ [12.5, 14.0] Gyr
- This is INDEPENDENT of H_0 measurement; NIE circular
- (CMB-derived t_universe = 13.8 Gyr uses ΛCDM → potential circularity for TGP test)

**Predicted H_0 range:**
- H_0(t = 12.5 Gyr) = 1/12.5 Gyr → convert to km/s/Mpc
- H_0(t = 14.0 Gyr) = 1/14.0 Gyr → convert to km/s/Mpc

**Pre-registered verdict logic:**
- If predicted range ⊂ [33.5, 146]: **PASS factor 2** (basic)
- If predicted range ∩ [67, 73] ≠ ∅: **PASS target** (excellent)
- If predicted range outside [33.5, 146]: **FAIL → HALT-B**

### T_P3_4 — Anti-cherry-pick sensitivity check

**Hypothesis:** Verify result NIE depends on extreme t_universe choices.

**Pre-registered:**
- Compute H_0(t_universe ∈ [12.0, 14.5] Gyr) range
- All values must fit within [33.5, 146] for honest PASS
- If only extreme t values fit → PARTIAL z honest declaration

---

## §2 — Anti-Lakatos LOCK declarations

**CRITICAL after Phase 2 favorable finding:**

1. **NIE modify F-γ-3 threshold ex post.** [33.5, 146] LOCKED Phase 0; pre-registered factor 2 anti-cherry-pick. STAYS.

2. **NIE redefine "TGP-native derivation" mid-cycle.** Geometric derivation (H = 1/t_universe via frontier R = c·t) was pre-registered Phase 0 §4.4 resolution #2. STAYS.

3. **NIE cherry-pick t_universe.** Use INDEPENDENT stellar ages anchor; not CMB-derived (circularity risk).

4. **Hubble tension correlation = OBSERVATION ONLY.** NIE claim TGP resolves tension; że wynik jest closer SH0ES vs Planck = observation; for future investigation, NIE concluded here.

5. **Honest disposition pre-declared:**
   - Predicted H_0 ~ 70 km/s/Mpc → factor 2 PASS expected
   - PASS verdict = **NATIVE GEOMETRIC DERIVATION + observational input t_universe**
   - PARTIAL caveat (Interpretation B): t_universe-dependent prediction

---

## §3 — DEC tracking

- DEC 2 reserved (numerical method) — for Phase 3 simple geometric calculation; **likely NIE expended** (analytical sufficient)
- DEC 3 reserved (anisotropic Bianchi fallback) — likely unused

---

## §4 — Computational plan (Phase3_sympy.py)

**Tool:** sympy symbolic + numerical evaluation.

**Sections:**
1. T_P3_1: R(t) = c·t consistency check
2. T_P3_2: H(t) = 1/t symbolic derivation
3. T_P3_3: F-γ-3 numerical verdict z stellar age anchor
4. T_P3_4: Sensitivity sweep + cherry-pick check
5. F-γ-3 verdict declaration (PASS / PARTIAL / FAIL)

**Result file:** Phase3_results.md.

---

## §5 — §3.6 extension compliance

- §3.6.6 Sign: H > 0 ✓
- §3.6.7 DoF: Phase 3 NIE fit (one-parameter geometric formula); no fit comparison needed
- §3.6.8 Implicit: late-time + matter-dominated late epoch + frontier-at-c
- §3.6.9 Precision: numerical computation z literal threshold
- §3.6.10 Pattern monitoring: ongoing

---

## §6 — Phase 3 risk register

| ID | Risk | Mitigation | Severity |
|----|------|-----------|----------|
| RP3-1 | Geometric prediction matches "too well" → suspect ΛCDM-fitting drift | Anti-Lakatos check #3 (use independent anchor) | LOW |
| RP3-2 | Hubble tension speculation overreach | §2 declaration #4: observation only, NIE claim | LOW |
| RP3-3 | t_universe ambiguity (12.5 vs 14.0 Gyr) → wide H_0 range | T_P3_4 sensitivity sweep within reasonable range | LOW |
| RP3-4 | Interpretation A vs B ambiguity → unclear verdict | §2 declaration #5 pre-resolves | LOW |

Phase 3 RISK PROFILE: LOW (Phase 2 already gave numerical result; Phase 3 = verdict formalization).

---

**Phase 3 plan LOCKED 2026-05-23. Anti-Lakatos LOCK reinforced. Ready dla Phase3_sympy.py.**
