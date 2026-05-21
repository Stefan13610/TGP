---
title: "Phase 1 results — first-principles symbolic derivation native ρ ≡ -T^μ_μ/c_0² z TGP axiom ax:metric-coupling"
date: 2026-05-13
parent: "[[./README.md]]"
type: phase-results
phase: 1
sympy_total: "11/11 PASS (100%)"
substance_metrics: "9 FIRST_PRINCIPLES (81.8%) / 2 LITERATURE_ANCHORED (18.2%) / 0 HARDCODED True; 100% non-trivial"
substantive_correction: "Factor-2 error fixed: predecessor (1 - v²/(2c²)) → first-principles γ⁻² (1 - v²/c²) at O(v²/c²)"
---

# Phase 1 results — first-principles symbolic derivation

## §0 — Headline

```
████████████████████████████████████████████████████████████████████
█                                                                  █
█  op-L01-N3-retrofit-native-SPARC-2026-05-13                      █
█                                                                  █
█  Phase 1 sympy: 11/11 PASS (100%)                                █
█  Substance: 9 FP / 2 LIT / 0 hardcoded True                      █
█  Non-trivial: 100%                                               █
█                                                                  █
█  vs predecessor (D-ALGEBRAIC_MIMICRY): +9 FP, -5 hardcoded       █
█                                                                  █
█  Substantive correction caught: factor-2 in non-rel expansion    █
█                                                                  █
████████████████████████████████████████████████████████████████████
```

## §1 — Test results per ID

| Test | Klasa | Status | Substance |
|---|---|---|---|
| T1a | FIRST_PRINCIPLES | PASS | 4-velocity normalization u^μ u_μ = -c² (mostly-plus signature) |
| T1b | FIRST_PRINCIPLES | PASS | Static fluid (v=0): T_00 = ρ·c² (energy density at rest) |
| T2 | FIRST_PRINCIPLES | PASS | Trace T^μ_μ = -ρ·c² + 3p (general perfect fluid, mostly-plus) |
| T3a | FIRST_PRINCIPLES | PASS | Dust limit: T^μ_μ_dust = -ρ·c² (v-independent EXACT) |
| T3b | FIRST_PRINCIPLES | PASS | Dust limit: ρ_TGP ≡ -T^μ_μ/c² = ρ_rest EXACT |
| T4a | FIRST_PRINCIPLES | PASS | Lorentz boost u^μ z velocity w=v → fluid rest frame u'^μ = (1,0,0,0) |
| T4b | FIRST_PRINCIPLES | PASS | T^μ_μ Lorentz invariant (boosted = rest frame value) |
| T5 | FIRST_PRINCIPLES | PASS | Non-rel: ρ_rest = ρ_local·(1 - v²/c²) z γ⁻² (correct relativistic relationship) |
| T6 | LITERATURE_ANCHORED | PASS | SPARC stars v=200 km/s: v²/c² ≈ 4.45·10⁻⁷ (<< 1%) |
| T7 | LITERATURE_ANCHORED | PASS | HI gas v_thermal=1 km/s: v²/c² ≈ 1.11·10⁻¹¹ (utterly negligible) |
| T8 | FIRST_PRINCIPLES | PASS | ax:metric-coupling consistency: T_μν_dust via δS_mat/δg^μν matches perfect fluid form EXACT (zero matrix difference) |

**Structural declarations (separate from PASS total):**

- T9: S05 single-Φ preservation — algebraic declaration
- T10: Scope clarification (galactic-disk; cluster + near-SMBH outside)

## §2 — Substance vs predecessor

| Metric | Predecessor N3 (2026-05-11) | This retrofit (2026-05-13) | Delta |
|---|---|---|---|
| Total sympy tests | 8 | 11 | +3 |
| PASS rate | 8/8 (100%) | 11/11 (100%) | preserved |
| FIRST_PRINCIPLES | 0 | 9 | **+9** |
| LITERATURE_ANCHORED | 3 | 2 | -1 |
| HARDCODED `T_pass = True` | 5 | **0** | **-5** |
| Non-trivial substance | ~25% (only 2 unit checks + 1 trivial subst) | **100%** | **+75pp** |

**Per audit_2026-05-11 §4.3 baseline:**
- Cohort 2026-05-11: 0/112 FIRST_PRINCIPLES testów
- This retrofit: 9/11 FIRST_PRINCIPLES testów (**first cycle z substantive FP majority dla L01 family**)

## §3 — Key substantive correction caught

**Predecessor cycle claim (from L01 NEEDS §T.3 i Phase_FINAL §4.1):**

```
ρ_SPARC = ρ_baryon ≈ ρ_rest · (1 - v²/(2c²))
                            ^^^^^^^^^^^^^^^^
                            factor 1/2 — WRONG
```

**First-principles symbolic derivation (this retrofit, T5):**

Starting z T^00 = γ²·ρ_rest·c² + p·O(v²/c²) (mostly-plus, perfect fluid):

```
ρ_local ≡ T^00 / c² = γ² · ρ_rest + ...
ρ_rest = ρ_local · γ⁻² = ρ_local · (1 - v²/c² + O(v⁴/c⁴))
                                    ^^^^^^^^^
                                    NO factor 1/2
```

**Factor 2 difference w correction term.** Predecessor "(1 - v²/(2c²))" wynika z conflation
γ⁻¹ vs γ⁻² (kinetic energy density `(1/2)ρv²` w mechanics vs energy-momentum tensor T_00).

**Numerical consequence:**
- Predecessor claim galactic correction: v²/(2c²) ~ 2.2·10⁻⁷
- First-principles correction: v²/c² ~ 4.45·10⁻⁷ (**factor 2 wyższe**)

Oba są << 1% target (cycle scope preserved), ALE first-principles wartość jest correct
relativistically. To jest **substantive finding z retrofit**: predecessor algebraic mimicry
faktycznie zamaskował factor-2 systematic error, NIE tylko procedural drift.

## §4 — Six P-requirements final status (Phase 1)

| # | Requirement | Resolution |
|---|---|---|
| **P1** | Perfect fluid T_μν decomposition z 4-velocity u^μ symbolic | ✅ T1a + T1b (FP) |
| **P2** | ρ ≡ -T^μ_μ/c_0² consistent dust limit → ρ_TGP = ρ_rest EXACT | ✅ T2 + T3a + T3b (FP) |
| **P3** | Non-relativistic correction explicit Lorentz boost | ✅ T4a + T4b + T5 (FP) — factor-2 corrected |
| **P4** | SPARC ρ_baryon mapping consistent | ✅ T6 + T7 (LIT, SPARC OOM check) |
| **P5** | No double-counting vs TGP-emergent DM (S05) | ✅ T9 (DECLARATIVE — structural separability g_eff[Φ̄] gravitational vs ρ_baryon matter) |
| **P6** | S05 single-Φ axiom preserved | ✅ T8 (FP — ax:metric-coupling consistency) + T9 (DECLARATIVE) |

**6/6 RESOLVED with substance (verifiable evidence per P).**

## §5 — Risks final status (Phase 1)

| Risk | Resolution |
|---|---|
| R1 (Lorentz boost expansion higher-order) | T5 used `sympy.series(..., 4).removeO()` — explicit O(v⁴) cutoff documented |
| R2 (Bullet Cluster scope) | T10 declarative scope clarification |
| R3 (ρ_baryon column independence) | T9 declarative S05 preservation |
| R4 (MOND chi²_red benchmark) | DEFERRED do Phase 2 / Phase FINAL — Phase 1 dust-limit derivation orthogonal to MOND comparison |

3/4 closed Phase 1; R4 deferred do następnej fazy (post-Phase-1 nie blocker dla Phase 1 closure).

## §6 — L1 native observable status

**Wypełnione w Phase 1:**
- Perfect fluid T_μν symbolic decomposition: ✅
- ρ ≡ -T^μ_μ/c_0² definition consistent z ax:metric-coupling: ✅
- Dust limit ρ_TGP = ρ_rest EXACT: ✅
- Non-relativistic correction explicit: ✅
- SPARC galactic regime OOM check: ✅

**Pozostałe (Phase 2 lub Phase FINAL):**
- Faktyczne 175-galaxy SPARC chi²_red TGP computation: **deferred** — wymaga galaxy_scaling
  cycles + N-body symulacja (out-of-scope dla niniejszego retrofit; predecessor cycle scope też
  nie zawierał faktycznego fittingu)
- L2 framework reduction (Newton-limit + MOND comparison): **deferred do Phase FINAL**
- L3 falsification map status update: **deferred do Phase FINAL**

## §7 — Audit invariant check

Per `meta/CALIBRATION_PROTOCOL.md` §4.4 — adversarial verification dla Phase 1.

**Self-audit checklist (light-touch dla Phase 1):**

| § | Question | Answer |
|---|---|---|
| 4.4.2(a) | §3 red flags w tym sympy? | None detected. Pure GR symbolic; explicit mostly-plus convention; no BD-form |
| 4.4.2(b) | §4 form-meaning mismatch? | None — perfect fluid T_μν z u^μ jest standardowy GR (per Wald §4.3), z g_eff[Φ] coupling przez ax:metric-coupling |
| 4.4.2(c) | ASK-RULE triggers? | None fired Phase 1 |
| 4.4.2(d) | Missing §2 patterns? | Pattern 2.1 (g_eff linearization) relevant ale NIE used Phase 1 (dust limit nie wymaga); Pattern 2.4 (collective gradient strain) relevant for Phase FINAL L2 reduction (deferred) |

**Self-audit verdict:** ✅ NO BD-DRIFT detected. Sympy operates per TGP-native S05 +
ax:metric-coupling principles.

**Recommendation:** Phase FINAL closure powinien spawn independent adversarial subagent dla
final audit per CALIBRATION_PROTOCOL §4.4 binding (independent re-classification 11 tests).

## §8 — Status & next steps

**Phase 1 GATE: ✅ OPEN — 11/11 PASS, 9 FP / 2 LIT / 0 hardcoded; 6/6 P-requirements RESOLVED z substance.**

**Next steps (NIE w niniejszej sesji — pending user authorization "active"):**

1. **Phase FINAL** closure z:
   - L2 framework reduction: Newton-limit (ρ_baryon → 4πG·ρ Poisson) + MOND simple comparison
   - L3 falsification map check: SPARC chi²_red benchmark (defer faktyczne fitting do
     galaxy_scaling cycle)
   - Adversarial audit subagent per CALIBRATION_PROTOCOL §4.4
   - Final claim_status determination: target **A−** (STRUCTURAL_DERIVED_NATIVE z L2
     deferred) lub **A** (jeśli L2 closed Phase FINAL)

2. **Cross-cycle propagation:**
   - L01 NEEDS §N3 status update: D-downgraded → **retrofit-A−** (post-Phase-FINAL)
   - galaxy_scaling cycles cross-cycle consistency check
   - nbody/ documentation update z explicit T^μ_μ → ρ mapping (poprawiony do γ⁻²)

3. **Replication template:** Sympy substance pattern z tego retrofit (≥75% FP, no hidden True,
   classification explicit) jako reference dla N1-EM, N2-QCD, N4-Higgs, N5-EW retrofits.

---

**Phase 1 sympy verification COMPLETE.** L1 native first-principles derivation `ρ ≡
-T^μ_μ/c_0²` z TGP axiom ax:metric-coupling **substantywnie verified** — replaces D-downgraded
predecessor algebraic mimicry. Substantive factor-2 correction caught dla non-relativistic
expansion.

**Cycle status:** Phase 1 closed; folder_status pozostaje `parking` aż user authorization
"active" + WIP slot. Cycle gotowy dla Phase FINAL spawn.
