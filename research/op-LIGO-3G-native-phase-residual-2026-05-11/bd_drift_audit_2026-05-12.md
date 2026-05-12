---
title: "bd-drift-audit — op-LIGO-3G-native-phase-residual Phase 1+2+3 (mid-cycle, 2026-05-12)"
date: 2026-05-12
type: adversarial-audit
status: IMMUTABLE — append-only adversarial verification record
audit_protocol: meta/CALIBRATION_PROTOCOL.md §4.4
binding_scope: "Phase 1+2+3 sympy substance audit; mid-cycle (przed Phase 4 spawn)"
audit_method: independent-read-of-sympy-code-only (NO results.md or README pre-read)
parent: "[[./README.md]]"
---

## §1 — Audit scope

**Read (allowed):**
- `Phase1_sympy.py` (748 lines, 13 tests)
- `Phase2_sympy.py` (1103 lines, 14 tests)
- `Phase3_sympy.py` (722 lines, 9 tests)
- `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §3 red flags + Pattern 2.2/2.5 sections (reference)
- File listing of cycle directory

**Intentionally NOT read (independence requirement):**
- Phase1_results.md / Phase2_results.md / Phase3_results.md (cycle self-verdicts)
- README.md (cycle claims, P-requirements progress)
- STATE.md (cycle status row)
- Phase1_sympy.txt / Phase2_sympy.txt / Phase3_sympy.txt (run output — relied on assertion-based sympy.py code analysis instead)
- AUDIT_2026-05-11_sympy_substance.md, CYCLE_KICKOFF_TEMPLATE.md (referenced indirectly only; minimal grep)

Audit examined every test individually for: (a) classification accuracy, (b) tautology risk, (c) BD-drift markers per Pattern 2.2/2.5 §3 red flags, (d) `T_pass = True` literal hardcoding, (e) declarative-vs-derivative honesty.

---

## §2 — Per-phase analysis

### §2.1 — Phase 1 (13 tests)

**T1 (FP1) — Newton emergence z covariant Phi-EOM linearization** (Lines 72-149)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: PARTIAL-DOWNGRADE → MIXED FP/LIT**
- Evidence: Lines 92-118 derive V'(Phi_0), vacuum condition β=γ, V''(Phi_0)=β — these ARE genuine sympy derivations (FP-quality). HOWEVER, line 138 `T1_newton_coeff = sp.solve(G_identification, q_sym)[0]` solves the IDENTIFICATION `q - 4*pi*G/c^2 = 0` for q, then verifies `q - 4*pi*G/c^2 == 0` → **this specific sub-test (T1_coupling_pass) is tautological** (lines 135-139: identification written, then verified). The vacuum/mass/linearization parts are legitimately FP.
- Drift marker: Type A.3 (literature-anchored identification verified as identity — partial; sub-step only, not full test).
- Overall: T1 still earns FP status because the *substantive* steps (vacuum, mass, linearization expansion) are real derivations. Sub-tautology is minor.

**T2 (FP2) — Newton force via momentum-flux integral** (Lines 152-239)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT** (substantive concern)
- Evidence: Lines 194-202: compute `grad_phi_2 at body 1` from `phi_2 = -G*m_2/(r_12 - xi)` via `diff`. Line 205-213: **the "momentum-flux integral" is NEVER actually computed**. Author writes (line 207): "Gauss-theorem reduction ... F^i_1 = -m_1 * (∂^i phi_2)|_{body 1}". This is the *test particle in external field* formula, written by hand. No ∮T^{0i}dS_i is computed; the result F = G*m_1*m_2/r_12^2 is recovered because the formula F = -m·grad·φ was *inserted* (line 213).
- Drift marker: **Type A.3 (literature-as-derivation) + Type B.5 (BD propagator avoided in form, but momentum-flux NOT actually demonstrated).** Per Pattern 2.2 §2.2.2 Step 4-5, the unique TGP-native claim requires computing ∮T^{ij}n_j dA, NOT the geodesic force formula. The author acknowledges this on line 207 ("Gauss-theorem reduction") but doesn't carry it out.
- Severity: MEDIUM — does not break the test (Newton 3rd law verification is symbolic), but mis-classifies test as FP when it is structurally LIT (Newton form recovered, momentum-flux mechanism claimed but not symbolically verified).

**T3 (FP3) — m_Phi_eff(r) environment-dependent** (Lines 242-311)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: CONFIRM** (genuine FP)
- Evidence: Lines 268-282 actually substitute Phi_local = Phi_0*(1+h_r) into V'' and verify expansion `1 + 4h + 3h^2` matches symbolic computation. Limit r→∞, ∂/∂r ≠ 0, ∂/∂m_1 ≠ 0 are all real symbolic operations. Pattern 2.5 §2.5 implementation faithful.

**T4 — sigma_cross_12 anisotropic uniaxial inheritance** (Lines 314-381)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM**
- Evidence: Lines 332-348 build sigma_cross matrix from explicit 2-body Newtonian potentials, then numerically verify uniaxial pattern at 3 sample points. Honest inheritance verification.

**T5 — sigma_cross_12 single-source limit** (Lines 384-411)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM**
- Evidence: Substitutes m_2 → 0 in pre-computed sigma_cross, checks all 9 entries vanish via simplify. Real sympy work.

**T6 — V_orig dual-V structure** (Lines 414-453)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: PARTIAL DOWNGRADE-TO-TAUTOLOGY**
- Evidence: Lines 432-443: `V_of_Phi = -(beta/(3*Phi_0))*Phi^3 + (gamma/(4*Phi_0^2))*Phi^4` is *defined* at line 92, then T6 extracts coefficients and verifies they equal `-beta/(3*Phi_0)` and `gamma/(4*Phi_0^2)`. **This is tautology: extracting back what was assigned.** Same pattern as N1 T1 precedent in AUDIT §4.3.
- Drift marker: Type A.2 (tautology). Severity: LOW (the test functions as a CHECKSUM that author didn't typo the V definition; harmless but not LIT-grade).
- Recommendation: Re-classify as TAUTOLOGY-CHECKSUM (or DECLARATIVE).

**T7 — Retarded Green's function NOT Yukawa** (Lines 456-507)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT** (good content, wrong classification)
- Evidence: Lines 481-491 verify Laplacian of 1/r vanishes for r>0, and Yukawa satisfies (∇²-m²)G_Y=0. These are *textbook* Green's function facts (Laplace and Helmholtz equations). The test does explicitly check the anti-Yukawa criterion (NO exp factor in TGP-native form per Pattern 2.2 §2.2.3 Warning B), which is good defensive sympy. But: this is verification of *known mathematical identities*, not first-principles TGP derivation.
- Severity: LOW. Reclassify as LIT.

**T8 — Dimensional analysis** (Lines 510-537)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM** (correct, basic LIT).

**T9 — PN ordering O(v²/c²)** (Lines 540-576)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM** (PN convention check, honest LIT).

**T10 — c_0 * κ_σ = 4/3 EXACT** (Lines 579-602)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: PARTIAL DOWNGRADE-TO-TAUTOLOGY**
- Evidence: Lines 57-58: `c0_inh = 4*pi`, `kappa_sigma_inh = 1/(3*pi)`. Line 593: `product_inh = c0_inh * kappa_sigma_inh` = `4*pi * 1/(3*pi)` = `4/3`. Test verifies `4/3 - 4/3 == 0`. **This is an arithmetic checksum**, not even a substantive identity. The HONESTY here is the README §7.1 caveat (lines 59-61) explicitly flagging c_0/kappa_sigma as "C-heuristic-level locks (NIE FIRST_PRINCIPLES)". So author IS honest about LIT classification, but the test itself is a numerical equality of constants.
- Severity: LOW. Acceptable as LIT given the comment honesty, but it is mechanically a constant-equality check.

**T11 — g_eff ansatz {A,B,C} inheritance** (Lines 605-651)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: PARTIAL DOWNGRADE-TO-TAUTOLOGY**
- Evidence: Line 628 defines `A_psi = 1 + a1*delta_psi + ...`. Line 642 verifies that the O(δ³) Taylor coefficient equals `a_3/6`. **This is verifying Taylor expansion derivative rule** (d³/dδ³ of (a_3·δ³/6) = a_3, divided by 3! = a_3/6). The b_1 = -a_1 constraint check on line 637 substitutes `b_1 → -a_1` then verifies `(-a_1 + a_1) == 0` — pure substitution tautology.
- Drift marker: Type A.2 (tautology after definition). Severity: LOW.

**T12 — m_Phi_eff vacuum limit** (Lines 654-683)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT**
- Evidence: Lines 671-672: `m_Phi_eff_sq_h = beta*(1 + 4*h + 3*h^2)` — note this expression is RE-DEFINED from h (not inherited from T3's actual `m_Phi_eff_sq_dimless`). Line 673 takes limit h→0 → β. Line 675 verifies β - β = 0. This is taking a limit of `beta*(1+0+0) = beta`. **Tautology after fresh re-definition.**
- Severity: LOW. The CONTENT is fine (it's a consistency check) but classification is too generous. Should be LIT (or even TAUTOLOGY).

**T13 — S05 single-Phi declarative** (Lines 686-708)
- Author: DECLARATIVE
- **Adversarial verdict: CONFIRM**
- Evidence: Line 702: `T13_pass = True` — but explicitly flagged as DECLARATIVE on line 696-697 and label `T13_status = "DECLARATIVE"`. This is HONEST declarative usage per §0.5b budget. Comment on line 703 explicitly notes "STRUCTURAL declaration (NOT sympy-derivable per se)".

### §2.2 — Phase 2 (14 tests)

**T1 (FP4) — sigma_cross_12 uniaxial DERIVATION** (Lines 110-211)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: CONFIRM with caveat**
- Evidence: Lines 132-163 build T^{ij}_cross from explicit gradient products, trace-remove. Numerical verification at 4 sample points (lines 176-186). The "derivation" of uniaxial pattern from trace-removal + y↔z symmetry IS genuine. Caveat: same numerical-verification approach as Phase 1 T4 (which was LIT). The DERIVATIVE step here (T_cross construction from L_TGP variation) is mostly hand-waved in comments (lines 145-153) — actual code is identical to Phase 1 T4's K_cross. Difference: Phase 2 T1 calls it "derivation z S05" while Phase 1 T4 called same code "inheritance verification". Same code, different label.
- Drift marker: Type A.4 (re-classification of identical computation). Severity: LOW-MEDIUM. The FP label is defensible if one reads the comment narrative; structurally the sympy is the same as T4.

**T2 (FP5) — Delta_phi(f) full chain** (Lines 214-377)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: CONFIRM**
- Evidence: Lines 232-362 chain: define Delta_e_2_native expression → diff E_b w.r.t. v → divide by P_GW → integrate → substitute v(f). Each step is genuine symbolic operation. The "Delta_e_2_native" expression IS hand-typed (line 245) from parent cycle structure, BUT subsequent operations (differentiation, integration, substitution) are real sympy. End-to-end chain symbolically verified. Genuine FP.
- Caveat: Step 4 expected forms (lines 344, 353) are hand-written then checked against sympy output — could be tautology IF the comment-author was reverse-engineering. However, the sequence dE/dv → dt/dv → integrand → integration IS a substantive chain.

**T3 (FP6) — sigma-coupling 2.5PN contribution** (Lines 380-468)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT**
- Evidence: The "TT-projection of grad cross-terms" claimed in test title is NOT computed. Lines 401-415: comments describe TT projector P^ij = δ^ij - n^i n^j, but the actual sympy code (lines 417, 427) just checks: (a) sigma_xx_axis is non-zero numerically, (b) `Delta_e_2_sigma_only = c_0_test * kappa_sigma_test` is defined (line 428) then multiplied by 45/16 (line 432), and `beta_ppE - 45/16*c_0*kappa_sigma` is verified zero — pure substitution tautology. (c) Linearity check via d²/dc_0² of a linear-in-c_0 expression. (d) "no exp(-m*r) screening" → `not sigma_form.has(sp.exp)` — but sigma_form is built from G*m/r terms which obviously have no exp; this check would pass for ANY non-Yukawa form.
- Drift marker: Type A.2 (tautology) + Type A.3 (claim-not-derivation). The 45/16·c_0·κ_σ "match" is between hand-written expression and same hand-written expression.
- Severity: MEDIUM. The test does NOT verify TT-projection mechanism; it confirms arithmetic.

**T4 — Geodesic preserves S05** (Lines 471-535)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT**
- Evidence: Lines 491-502 compute Christoffel Gamma^x_xx symbolically from B(h) Taylor. Line 506: `gamma_free = Gamma_xxx_formal.free_symbols` then checks `g_eff_independent not in gamma_free`. But `g_eff_independent` is a symbol DEFINED on line 507 that was never inserted into Gamma_xxx_formal. **Trivially true by construction.** Line 511-513: vacuum limit h=0, h'=0 → Gamma→0 (any expression vanishes when its argument is 0). Line 518: leading order check `-1/2 * b_1 * hp` matches its own Taylor derivative.
- Drift marker: Type A.2 (tautology). The S05 claim is real but the sympy is verifying mechanical Taylor properties.
- Severity: LOW-MEDIUM.

**T5 — PN expansion to 2.5PN** (Lines 538-586)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: PARTIAL DOWNGRADE-TO-TAUTOLOGY**
- Evidence: Lines 562-572: defines E_b_test with e_1_test, e_2_test, Taylor-expands, then verifies the coefficients of v² are -η/2 (= what was put in), v⁴ = -η·e_1/2 (= what was put in), v⁶ = -η·e_2/2 (= what was put in). Pure echo of input.
- Severity: LOW.

**T6 — SPA consistency** (Lines 589-637)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: PARTIAL DOWNGRADE-TO-TAUTOLOGY**
- Evidence: Line 611-613: checks `Delta_phi_eta_qtr - (-15/4*Delta_e_2_native/(M*v)) == 0`. But `Delta_phi_eta_qtr` was computed in T2 to equal exactly `-15/4*Delta_e_2_native/(M*v)` (line 360-361). Re-verifying same expression. Line 627: `T6_consistency = True` literal (line 627: `T6_consistency = True  # structural form linear in Delta_e_2 verified above`). 
- Drift marker: **Type A.1 (literal `True` hardcoded)** on line 627. NOT flagged as DECLARATIVE; flows into `T6_pass`. This contradicts author's "0 hardcoded T_pass = True" claim.
- Severity: MEDIUM. Hidden hardcoded True passing as LIT classification.

**T7 — Kepler r_12(f)** (Lines 640-684)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM** (Kepler's law substitutions, honest LIT).

**T8 — Delta_phi at f_ISCO** (Lines 687-722)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM** (numerical sanity at Schwarzschild ISCO).

**T9 — sigma_TT decomposition h_+/h_×** (Lines 725-777)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM** (real TT decomposition; honest LIT).

**T10 — Native coefs sensitivity Fisher prep** (Lines 780-843)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT**
- Evidence: Lines 800-803 take ∂/∂a_3, ∂/∂xi_3, ∂/∂c_0, ∂/∂kappa of `Delta_phi_can`. Verifies each is non-zero at test point. The comment on lines 825-826 acknowledges: **"all proportional via Delta_e_2 at b=-1 level"** — i.e., the Fisher matrix IS degenerate. Test only checks `T10_all_nonzero` (each derivative ≠ 0), which is trivially true for any non-constant function with respect to each of its parameters. Does NOT actually verify Fisher non-degeneracy; the conclusion is the OPPOSITE (degeneracy admitted).
- Drift marker: Type A.3 (test name claims more than test verifies). Severity: LOW-MEDIUM. Author is HONEST in comment that real Fisher non-degeneracy needs higher PN; should be classified LIT (preparation, not verification).

**T11 — M9.1'' Path 2 anchor substitution** (Lines 846-895)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM**
- Evidence: Substitutes a_3=36, xi_3=5/24, c_0*kappa=4/3 into Delta_e_2_native, verifies result = 0. Real arithmetic (the cancellation: -4·(5/24) + 4 - 36/8 + 4/3 = -5/6+4-9/2+4/3 = 0). Honest LIT.

**T12 — GR limit multi-parameter** (Lines 898-959)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: CONFIRM** (real sympy solve, free symbols check, multi-D surface).

**T13 — Cross-cycle consistency parent Phase 4** (Lines 962-1012)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: PARTIAL DOWNGRADE-TO-TAUTOLOGY**
- Evidence: Lines 988-993 define `beta_ppE_phase2 = (45/16)*Delta_e_2_native` and `beta_ppE_phase4 = (45/16)*Delta_e_2_diag + (45/16)*c_0*kappa`. Both forms factor 45/16 are HAND-WRITTEN; Delta_e_2_native ≡ Delta_e_2_diag + c_0*kappa by construction. Verification is then arithmetic. The "cross-cycle match" is a definitional restatement.
- Severity: LOW. Honest LIT category given inheritance nature.

**T14 — S05 declarative through Phase 2** (Lines 1015-1043)
- Author: DECLARATIVE
- **Adversarial verdict: CONFIRM**
- Evidence: Line 1035: `T14_pass = True` BUT explicitly flagged DECLARATIVE on line 1027. Honest budget usage.

### §2.3 — Phase 3 (9 tests)

**T1 (FP7) — Analytical-exact SPA projection** (Lines 100-186)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT**
- Evidence: Each "step" is substitution of pre-known constants:
  - Step 1 (line 129): `(3/(128*eta))*delta_alpha_4 - (3/(128*eta))*delta_alpha_4 == 0` (tautology after definition).
  - Step 2 (line 139-144): `delta_alpha_4 = 30*Delta_e_2 - 20*Delta_e_1*p_1`, then substitute Delta_e_1=0, get 30*Delta_e_2. The CF Eq. 3.18 form is hand-typed from literature; the substitution is arithmetic.
  - Step 3 (line 148-151): `(3/128)*30/eta = 45/(64*eta)` — arithmetic.
  - Step 4 (line 155-157): substitute eta=1/4 → 45/16. Arithmetic.
  - Step 5 (line 162-166): substitute Delta_e_2_native_canonical. Arithmetic substitution.
  - Step 6 (line 170-171): re-verify Step 5. Same equality.
- Drift marker: Type A.3 (literature CF Eq. 3.18 → 30*Delta_e_2 hand-inserted, then "verified" as identity).
- Severity: MEDIUM. Test is a series of arithmetic checks on hand-typed literature formulas. Honest LIT, not FP.

**T2 (FP8) — Cross-cycle consistency parent Phase 4 LOCK** (Lines 189-274)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-LIT**
- Evidence: Lines 207-208 define `beta_ppE_phase4_lock = (45/16)*Delta_e_2_diag + (45/16)*c_0*kappa`. Line 218: `diff_with_parent = simplify(beta_ppE_native - beta_ppE_phase4_lock)` where beta_ppE_native = (45/16)*Delta_e_2_native_canonical from T1, and Delta_e_2_native_canonical = Delta_e_2_diag + c_0*kappa_sigma_sym (line 84-85). Thus the "cross-cycle identity" reduces to `(45/16)*(A+B) - (45/16)*A - (45/16)*B == 0` (arithmetic distribution). Coefficient checks (lines 224-237) re-extract coefficients of input definitions — tautology.
- Drift marker: Type A.2 + Type A.3. Severity: MEDIUM. This is the SAME form-match pattern as Phase 2 T13, also DOWNGRADE-TO-TAUTOLOGY.

**T3 (FP9) — VT-002 AF1 closure** (Lines 277-352)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: DOWNGRADE-TO-DECLARATIVE/LIT**
- Evidence: Criterion (a) uses `T1_step6_analytical_exact` from T1 (which we already downgraded). Criterion (b): 4 distinct substitution points → 3 distinct beta_ppE values. Verifies linear function takes different values for different inputs — trivially true for non-constant function. Criterion (c): re-derives 45/16 = 3·30/128/(1/4) — arithmetic. Criterion (d): same arithmetic. Criterion (e) is `T3_criterion_e = True` literal (line 338).
- Drift marker: **Type A.1 (literal `T_pass = True`** on line 338: `T3_criterion_e = True  # symbolic structure trivially dimensionless`). Hidden hardcoded True flowing into FP-classified test pass.
- Severity: MEDIUM-HIGH. This is exactly the antipattern §4.2 from AUDIT precedent.

**T4 — ppE basis convention** (Lines 355-394)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: DOWNGRADE-TO-TAUTOLOGY/DECLARATIVE**
- Evidence: Line 376: `phase_term_b_m1 = u_v ** b_val` with b_val=-1, then verifies it equals `1/v_pn`. Tautology. **Line 381: `T4_PN_order_correspondence = True` literal**. **Line 385: `T4_dimensional = True` literal**. Both flow into `T4_pass`.
- Drift marker: **Type A.1 (literal `True` × 2 hardcoded)** on lines 381, 385. NOT flagged as DECLARATIVE.
- Severity: MEDIUM-HIGH. Hidden hardcoded True in LIT-classified test.

**T5 — M9.1'' beta_ppE = -15/4** (Lines 397-442)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM**
- Evidence: Real arithmetic substitution to compute -4/3 → (45/16)*(-4/3) = -15/4. Plus inequality check against GWTC-3 bound. Honest LIT.

**T6 — GWTC-3 1σ window** (Lines 445-491)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: CONFIRM** (real solve for y-window, numerical comparison).

**T7 — Phase 2 vs Phase 3 SPA convention** (Lines 494-553)
- Author: LITERATURE_ANCHORED
- **Adversarial verdict: PARTIAL DOWNGRADE-TO-TAUTOLOGY**
- Evidence: Line 520: ratio = -(15/4)/(45/16) = -4/3 (arithmetic). **Line 533: `T7_convention_explained = True` literal** flows into `T7_pass`.
- Drift marker: **Type A.1 (literal True hardcoded)** on line 533.
- Severity: MEDIUM. Hidden hardcoded True.

**T8 — GR limit multi-parameter** (Lines 556-620)
- Author: FIRST_PRINCIPLES
- **Adversarial verdict: CONFIRM** (real solve, free_symbols check).

**T9 — S05 declarative** (Lines 623-650)
- Author: DECLARATIVE
- **Adversarial verdict: CONFIRM** (explicitly flagged DECLARATIVE).

---

## §3 — Drift markers found (cumulative)

| # | Marker type | Phase | Test | Severity | Evidence | Recommended action |
|---|-------------|-------|------|----------|----------|--------------------|
| 1 | A.1 literal `T_pass = True` (hidden) | Phase 2 | T6 line 627 | MEDIUM | `T6_consistency = True` flows into T6_pass; classified LIT but hardcoded | Reclassify component as DECLARATIVE-substep; document |
| 2 | A.1 literal `T_pass = True` (hidden) | Phase 3 | T3 line 338 | MEDIUM | `T3_criterion_e = True` flows into FP7 pass | Reclassify FP9 → LIT or split criterion |
| 3 | A.1 literal `T_pass = True` (hidden) | Phase 3 | T4 lines 381, 385 | MEDIUM | Two `True` literals flow into T4 LIT pass | Replace with real check or reclassify |
| 4 | A.1 literal `T_pass = True` (hidden) | Phase 3 | T7 line 533 | MEDIUM | `T7_convention_explained = True` flows into T7 LIT pass | Same |
| 5 | A.3 momentum-flux mechanism claimed, formula inserted | Phase 1 | T2 (FP2) | MEDIUM | Test recovers Newton via `F = -m·grad·phi` test-particle formula (line 213), NOT via ∮T^{0i}dS_i computation despite test name. Pattern 2.2 §2.2.2 Step 4-5 requires actual stress-energy surface integral. | Downgrade FP2 → LIT; OR add explicit T^{ij} surface integral computation |
| 6 | A.2 tautology | Phase 1 | T6 | LOW | V coefficients extracted from definition match definition | Reclassify as TAUTOLOGY-CHECKSUM |
| 7 | A.2 tautology | Phase 1 | T10, T11, T12 | LOW each | Arithmetic of constants / Taylor coefficient extraction / fresh re-definition limit | Reclassify as TAUTOLOGY-CHECKSUM |
| 8 | A.2/A.3 same-form match | Phase 2 | T6, T13 | LOW each | Forms hand-written then equality verified | Acknowledge LIT-level checksum |
| 9 | A.3 claim>verification | Phase 2 | T3 (FP6) | MEDIUM | "TT-projection mechanism" claimed; sympy only verifies `45/16*c_0*kappa - 45/16*c_0*kappa == 0` and `not has(exp)` | Downgrade FP6 → LIT |
| 10 | A.2 tautology by construction | Phase 2 | T4 | LOW-MEDIUM | `g_eff_independent not in Gamma_xxx_formal.free_symbols` — symbol never inserted | Reclassify FP supporting → LIT |
| 11 | A.3 test claims more than verified | Phase 2 | T10 | LOW | Fisher non-degeneracy claimed; only `each_partial != 0` checked; degeneracy ADMITTED in comment | Reclassify FP supporting → LIT |
| 12 | A.3 substitution chain literature-anchored | Phase 3 | T1 (FP7), T2 (FP8) | MEDIUM | CF Eq. 3.18 inserted (30*Delta_e_2); subsequent steps arithmetic substitution; "analytical-exact" reduces to (45/16)·(A+B)=(45/16)·A+(45/16)·B distribution | Downgrade FP7, FP8 → LIT |
| 13 | B.5 BD-Yukawa avoided syntactically | Phase 1 T7, Phase 2 T3/T4 | LOW | `not expr.has(sp.exp)` is trivially true for any non-Yukawa form built from G·m/r terms | Acceptable defensive check; recognize no real anti-BD content |

**HIGH-severity flags:** 0 (no fundamental BD-drift; the framework is TGP-consistent in narrative but implementation has Tautology layering).

**MEDIUM-severity flags:** 7 (mostly tautology/hidden-True patterns).

**LOW-severity flags:** 6.

---

## §4 — Substance metrics independent count

**Author's claim** (from cycle README §0.5b / Phase summary lines):
- TOTAL: 36 tests
- FIRST_PRINCIPLES: 15 (41.7%)
- LITERATURE_ANCHORED: 18 (50%)
- DECLARATIVE: 3 (8.3%)
- TAUTOLOGY: 0 (claimed 0% literal hardcoded)

**Adversarial re-classification:**

Phase 1 (13 tests):
- FP confirmed: T1 (partial — accept), T3 — **2 FP**
- LIT confirmed/downgraded-to: T2 (downgrade from FP), T4, T5, T7 (downgrade from FP), T8, T9, T10, T12 (downgrade from FP) — **8 LIT**
- TAUTOLOGY-CHECKSUM: T6, T11 — **2 TAUT**
- DECLARATIVE: T13 — **1 DEC**

Phase 2 (14 tests):
- FP confirmed: T1 (FP4 — accept), T2 (FP5), T12 — **3 FP**
- LIT: T3 (downgrade FP6→LIT), T4 (downgrade FP supp→LIT), T7, T8, T9, T10 (downgrade FP supp→LIT), T11, T13 — **8 LIT**
- TAUTOLOGY: T5, T6 (hidden True) — **2 TAUT**
- DECLARATIVE: T14 — **1 DEC**

Phase 3 (9 tests):
- FP confirmed: T8 — **1 FP**
- LIT: T1 (downgrade FP7→LIT), T2 (downgrade FP8→LIT), T3 (downgrade FP9→LIT, hidden True), T4 (hidden True), T5, T6, T7 (hidden True) — **7 LIT**
- DECLARATIVE: T9 — **1 DEC**
- TAUTOLOGY: 0 (those with hidden True classified as LIT-with-flag)

**Adversarial totals (36 tests):**
- FIRST_PRINCIPLES: **6 / 36 (16.7%)** vs. claimed 15 (41.7%)
- LITERATURE_ANCHORED: **23 / 36 (63.9%)** vs. claimed 18 (50%)
- DECLARATIVE: **3 / 36 (8.3%)** matches author
- TAUTOLOGY (newly identified): **4 / 36 (11.1%)** vs. claimed 0

**Reclassification rate:** 9 substantive reclassifications (out of 36 = 25%). This exceeds the 5% PASS threshold per protocol.

**Hidden literal `T_pass = True` count:** 4 instances (Phase 2 T6, Phase 3 T3/T4/T7) — author's claim of "0% literal hardcoded" is FALSIFIED for hidden-component variables that flow into pass conditions. The visible top-level `T_n_pass = True` are absent (correct), but sub-flags `T_x_criterion = True` exist and propagate.

---

## §5 — Verdict

**AMENDMENT NEEDED**

Reasoning:
- Reclassification rate 25% (9/36) substantively exceeds 5% threshold for PASS.
- 4 instances of hidden hardcoded `_criterion = True` / `_consistency = True` / `_convention_explained = True` / `_dimensional = True` that flow into pass conditions — falsifies the "0 literal True" claim.
- Phase 1 T2 (FP2) substantive concern: the Pattern 2.2 momentum-flux MECHANISM (the *unique* TGP-native distinguisher from BD) is NOT computed; the test recovers Newton via the geodesic test-particle formula F = -m·∇φ, which is correct physics but NOT the claimed derivation. This is the single biggest Pattern-2.2-compliance concern.
- Phase 3 T1/T2 (FP7/FP8) are arithmetic substitution chains anchored in Cutler-Flanagan literature (Eq. 3.18 hand-inserted). The "analytical-exact" claim reduces to distributing (45/16) over a sum that was DEFINED as a sum. This is LIT-quality, not FP.
- No HIGH-severity BD-drift markers found. NO Yukawa propagator usage. NO m_Phi universal slip (Phase 1 T3 genuinely environment-dependent). S05 declared honestly via flagged DEC budget. Forced-zero parameters not relevant here.
- Substance IS real: Phase 2 T2 (FP5) full chain Phi-EOM → geodesic → dE/dt → phase IS a genuine multi-step symbolic derivation. Phase 1 T3 (Pattern 2.5 m_Phi_eff) IS first-principles. Phase 2 T12 / Phase 3 T8 (multi-D GR surface) are real solve operations.

Net assessment: the cycle has **legitimate FP substance** (~6 genuine FP tests = 16.7%) but **the headline "41.7% FP" claim overstates by 2.5×** via reclassification leniency and the "0 hardcoded True" claim is falsified by hidden subflag literals.

---

## §6 — Recommendation dla Phase 4 spawn decision

Per pre-locked decision rule (cycle README §7.6 #3 + Q8):

**Phase 4 spawn BLOCKED pending amendment.**

Specific amendments required before Phase 4 spawn:

1. **Mandatory:** Audit & remove/reclassify hidden literal `True` assignments to pass-flag components:
   - Phase 2 line 627 (`T6_consistency = True`)
   - Phase 3 line 338 (`T3_criterion_e = True`)
   - Phase 3 lines 381, 385 (`T4_PN_order_correspondence = True`, `T4_dimensional = True`)
   - Phase 3 line 533 (`T7_convention_explained = True`)
   Either implement real checks, mark those individual criteria as DECLARATIVE-subflag, or split test into FP-substantive + DEC-substep.

2. **Mandatory:** Reclassify the following from FP → LIT in cycle results documentation:
   - Phase 1 T2 (FP2) → LIT until ∮T^{0i}dS_i sphere integral is symbolically computed
   - Phase 2 T3 (FP6) → LIT (TT-projection claim not symbolically realized)
   - Phase 2 T4 (FP supp) → LIT (trivial-by-construction symbol-absence check)
   - Phase 2 T10 (FP supp) → LIT (Fisher non-degeneracy not verified; degeneracy admitted)
   - Phase 3 T1 (FP7), T2 (FP8), T3 (FP9) → LIT (CF literature anchor + arithmetic substitution; no genuine FP step)
   
   This drops adversarial FP count to 6 (16.7%) which is below cycle target but accurately reflects substance.

3. **Recommended:** For Phase 1 T2 (FP2), if author wishes to retain FP classification, ADD a symbolic computation of T^{ij}_{cross} = ∂^i·phi_2 ⊗ ∂^j·phi_1 integrated over a sphere surface (n^j area element) and demonstrate the result equals -m_1·∇phi_2|_{body 1} via Gauss theorem APPLIED IN SYMPY (not just claimed in comments). This would convert the claim into a Pattern 2.2 §2.2.2 verification.

4. **Recommended:** Phase 3 FP-claim depends critically on the inserted CF Eq. 3.18 (`delta_alpha_4 = 30*Delta_e_2 - 20*Delta_e_1*p_1`). If this is to be FP, the test should DERIVE the prefactor 30 from the SPA integration of `(dt/dv)·omega_GW(v)` with explicit binding-energy coefficient — not insert it as literature.

5. **Not required but advisable:** Add a TAUTOLOGY classification category to the test registry (separate from DECLARATIVE) so checksums vs. genuine declarations are visible. Currently the "0% TAUT" claim is achieved by reclassifying as LIT.

After amendments, re-run the audit; if reclassifications are accepted and hidden-True flags removed/legitimized, the cycle would be PASS-WITH-FLAGS (LOW severity) and Phase 4 spawn authorized.

---

## §7 — Sign-off

Audit conducted: 2026-05-12 (adversarial subagent, autonomous read of 3 sympy files; protocol per CALIBRATION_PROTOCOL §4.4).

Files read: Phase1_sympy.py (full, 748 lines), Phase2_sympy.py (full, 1103 lines in two chunks), Phase3_sympy.py (full, 722 lines), meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md (relevant grep §3 Pattern 2.2/2.5). Did NOT read any results.md, README.md cycle claims, or STATE.md to preserve independence.

Methodology: line-by-line classification check against §3 red flags + §4 anti-patterns + Pattern 2.2/2.5 mechanism criteria.

Audit invariant: dokument jest IMMUTABLE; append-only.

---

## §8 — bd-drift-audit ITERATION II — Amendment verification (2026-05-12)

Independent verification of post-amendment state (Scope A, mandatory items §6.1+§6.2 only).
Method: read Phase{1,2,3}_sympy.py + .txt post-amendment; NOT read Phase{n}_results.md,
README.md §7.5b, STATE.md, or any amendment_2026-05-12.md self-report. Compared directly
against §6 mandatory items.

### §A — Mandatory item 1 verification (4 hidden Trues)

- **Phase 2 T6 line ~627 (`T6_consistency`):** PASS — now real `sp.diff` linearity check. Code (lines 646-652) introduces `_Delta_e2_var`, computes `sp.diff(_phase_form, _Delta_e2_var)` and verifies it equals `-(15/4)/(M*v)` AND second derivative is zero. `T6_consistency_status = "DERIVED"`. Substantive sympy, not literal.
- **Phase 3 T3 crit_e line ~338 (now ~369):** PASS — `T3_criterion_e_status = "DECLARATIVE"` explicit flag (line 368), literal `True` retained on line 369 with explicit comment that this is DECLARATIVE-substep (units not tracked in sympy).
- **Phase 3 T4 PN_order line ~381 (now ~423):** PASS — `T4_PN_order_correspondence_status = "DECLARATIVE"` (line 422), `T4_PN_order_correspondence = True` annotated as DECLARATIVE-substep.
- **Phase 3 T4 dimensional line ~385 (now ~428):** PASS — `T4_dimensional_status = "DECLARATIVE"` (line 427), `T4_dimensional = True` annotated as DECLARATIVE-substep.
- **Phase 3 T7 convention line ~533 (now ~582):** PASS — `T7_convention_explained_status = "DECLARATIVE"` (line 581), annotated as DECLARATIVE-substep narrative.

### §B — Mandatory item 2 verification (7 reclassifications FP → LIT)

- **Phase 1 T2:** now `T2_classification = "LITERATURE_ANCHORED"` (line 179, with RECLASSIFICATION NOTE block lines 156-164).
- **Phase 2 T3:** now `T3_classification = "LITERATURE_ANCHORED"` (line 405).
- **Phase 2 T4:** now `T4_classification = "LITERATURE_ANCHORED"` (line 500).
- **Phase 2 T10:** now `T10_classification = "LITERATURE_ANCHORED"` (line 826).
- **Phase 3 T1 (FP7):** now `T1_classification = "LITERATURE_ANCHORED"` (line 125).
- **Phase 3 T2 (FP8):** now `T2_classification = "LITERATURE_ANCHORED"` (line 220).
- **Phase 3 T3 (FP9):** now `T3_classification = "LITERATURE_ANCHORED"` (line 318).

All 7 reclassifications confirmed; each carries an explicit RECLASSIFICATION NOTE referencing bd-drift-audit 2026-05-12 §6 item 2.

### §C — Sympy still passes (re-run)

- Phase 1: **13/13 PASS** (txt line 144). 4 FP / 8 LIT / 1 DEC.
- Phase 2: **14/14 PASS** (txt line 174). 3 FP / 10 LIT / 1 DEC.
- Phase 3: **9/9 PASS** (txt line 125). 1 FP / 7 LIT / 1 DEC.

### §D — Amended substance metrics (independent recount)

- TOTAL: 36 tests
- FIRST_PRINCIPLES: **8 / 36 (22.2%)** — Phase 1 (T1, T3, T7, T12) + Phase 2 (T1, T2, T12) + Phase 3 (T8). Matches Scope A expected post-amendment count.
- LITERATURE_ANCHORED: **25 / 36 (69.4%)**.
- DECLARATIVE: **3 / 36 (8.3%)** — T13 (P1), T14 (P2), T9 (P3).
- Hidden literal `True` at substantive (non-DECLARATIVE) level: **0** (only legitimate uses remain: 3 DECLARATIVE `T_pass = True` per budget; 4 DECLARATIVE-substep `True` explicitly flagged via `_status = "DECLARATIVE"`; 1 loop-sentinel `T5_all_zero = True` in Phase 1 line 407 falsified by inner loop). No new hidden Trues introduced.

### §E — Verdict (per pre-locked decision rule)

**PASS — amendments successfully addressed all mandatory items.**

- All 4 hidden Trues from §6 item 1 resolved: 1 replaced with substantive `sp.diff` check (Phase 2 T6), 4 explicitly annotated with `_status = "DECLARATIVE"` flags (Phase 3 T3/T4/T7).
- All 7 FP→LIT reclassifications from §6 item 2 applied with explicit RECLASSIFICATION NOTE comments.
- No new HIGH-severity issues introduced; no new hidden Trues; no upward reclassifications LIT→FP. Phase 1 T7+T12 correctly remain FP per Scope A (recommended-only, outside mandatory list).
- Independent recount reproduces 8 FP / 25 LIT / 3 DEC. Reclassification rate post-amendment: 0% beyond Scope A.

### §F — Phase 4 spawn recommendation

**AUTHORIZED** — Phase 4 spawn unblocked. Cycle now honestly documents 22.2% FP substance with all sub-flag literals either substantively verified or explicitly DECLARATIVE-flagged.

Sign-off: bd-drift-audit ITERATION II 2026-05-12, independent verification subagent. Files read: Phase{1,2,3}_sympy.py + .txt only; results.md / README §7.5b / STATE.md / amendment_2026-05-12.md NOT read (independence preserved).

---

## §9 — Iter III FINAL audit (2026-05-12, pre-closure)

Method: independent read of Phase4_sympy.py + Phase5_sympy.py (NEW since Iter II) +
.txt run outputs + spot-check grep on Phase1-3_sympy.py post-amendment state. NOT read:
Phase{4,5}_results.md, README.md §7.5e/§7.5f, STATE.md, Phase1-3_amendment_2026-05-12.md.
Per CALIBRATION_PROTOCOL §4.4 (Phase FINAL trigger, mandatory pre-closure).

### §9.A — Phase 4 audit (T1-T9)

- **T1 (FP10) native Fisher rank-1:** **CONFIRM** as genuine FP. Lines 154-220 build
  `alpha_vec = sp.Matrix([alpha_a3, alpha_xi3, alpha_y])` from explicit `sp.diff` of
  `Delta_e_2_native` w.r.t. each native param (alpha = (-1/8, -4, +1) emerges, NOT
  hand-coded). `Gamma = W * alpha_vec * alpha_vec.T` constructed. Substantive operations:
  `.rank()` (line 162 → 1), `.eigenvals()` (line 170 → `{1089*W/64: 1, 0: 2}`), explicit
  eigenvector check `Gamma * alpha_vec` (line 187), TWO independent null mode verifications
  (kernel + orthogonality + rank-2 null_matrix). 17 sympy-substantive sub-checks all
  contribute to T1_pass. Rank-1 emerges structurally from Phase 2 FP5 chain (Delta_phi
  linear in single combination) — NOT postulated.
- **T2 (FP11) beta_ppE↔native equivalence:** **CONFIRM** as genuine FP. Substantive ops:
  `sp.diff(beta_ppE_as_fn, DeltaE2_sym)` (line 264), `sp.solve` for a_3 in 3D constraint
  (line 309), hyperplane verification at 3 distinct (a_3, xi_3, y) points proving same
  Delta_e_2 (lines 296-304). The (16/45) factor is DERIVED via chain rule from (45/16) LOCK,
  not postulated. GWTC-3 anchor cross-check (4.808 sigma) matches both projection paths.
- **T3-T8 LIT classifications:** **CONFIRM overall**. T3-T7 are honest numerical anchor
  projections + dimensional / symmetry / determinant checks (`Gamma.det() == 0`, `Gamma == Gamma.T`).
  T8 inheritance attribution is real `sp.simplify` equality verification against Phase 2/3 forms.
  Two flagged literals `T5_psd_structure = True` (line 468, `T5_psd_status="STRUCTURAL"`) and
  `T6_reparam_structure = True` (line 530, `T6_reparam_status="DERIVED"`) — **neither is
  used in T5_pass / T6_pass** (verified: T5_pass uses only T5_symmetric+T5_det_zero, T6_pass
  uses only T6_reparam_forward+T6_reparam_inverse). Status flags are honest annotations, not
  hidden pass-flow. **No hidden Trues.**
- **T9 DECLARATIVE 3PN deferral (FP12):** **PROPER FLAG**. T9_classification="DECLARATIVE",
  T9_status="DECLARATIVE", explicit recovery scope to Phase 4b documented (lines 681-684).
  1 DEC out of 9 tests = 11.1%, within ≤15% budget.

### §9.B — Phase 5 audit (T1-T10)

- **T5 (FP13) network Fisher quadrature DERIVED:** **CONFIRM** as genuine FP. Lines 342-380
  build `Gamma_ETD = W_ETD * alpha_vec * alpha_vec.T` and `Gamma_CE = W_CE * alpha_vec * alpha_vec.T`,
  sum them, verify the sum equals `(W_ETD + W_CE) * alpha_vec * alpha_vec.T` element-by-element
  via `sp.simplify` (lines 354-357). `.rank()` checks on each per-detector matrix AND on
  network matrix (all = 1). Eigenvector preservation: `Gamma_net * alpha_vec == net_eig * alpha_vec`
  (lines 365-370). **Crucial step**: quadrature rule `sigma_net^-2 = sigma_ETD^-2 + sigma_CE^-2`
  derived symbolically from `(W_ETD + W_CE) * ||alpha||^2 = W_ETD*||alpha||^2 + W_CE*||alpha||^2`
  (line 379) — emerges from rank-1 additivity + Fisher linearity, NOT postulated from
  quadrature rule. Numerical cross-check: quadrature → 0.00409 matches direct ET+CE
  network Phase 2 calibration 0.00409.
- **T1-T4, T6-T9 LIT:** **CONFIRM overall**. T1 SNR scaling via `sp.diff(log_A_sq, ...)*Mc`
  recovers 5/3 / -2 / -7/3 exponents — real derivative ops, not echo. T2-T4, T6, T7
  numerical detector projections via Phase 3 LOCK; T8 N-event stacking + 1-yr stack
  yields 1.292e-05 sigma_DeltaE2 (~10^5 sigma vs M9.1''); T9 microradian thresholds
  ordering check. All honest LIT, no FP-claim inflation. No hidden Trues (only legitimate
  T10_pass=True with T10_status="DECLARATIVE" flag).
- **T10 DECLARATIVE degeneracy_factor=5:** **PROPER**. Yagi-Yunes 2016 inheritance
  explicitly documented with assumption boundary (line 644-649). 1 DEC of 10 = 10%, at
  budget edge but compliant per §0.5b ≤10% target.

### §9.C — Phase 1-3 regression spot-check

- **Phase 2 T6 sp.diff replacement still substantive:** **YES** — lines 646-652 retain
  `_Delta_e2_var = sp.Symbol(...)`, `_phase_form = -(15/4)*var/(M*v)`, real
  `sp.diff(_phase_form, _Delta_e2_var)` + second-derivative check + `T6_consistency_status = "DERIVED"`.
- **Phase 3 T3/T4/T7 DECLARATIVE flags still present:** **YES** — `T3_criterion_e_status="DECLARATIVE"`
  (line 368), `T4_PN_order_correspondence_status="DECLARATIVE"` (line 422), `T4_dimensional_status="DECLARATIVE"`
  (line 427), `T7_convention_explained_status="DECLARATIVE"` (line 581). All bd-drift-audit
  2026-05-12 comments preserved verbatim.
- **Phase 1 T7/T12 preserved as FP per Scope A:** **YES** — `T7_classification="FIRST_PRINCIPLES"`
  (line 480), `T12_classification="FIRST_PRINCIPLES"` (line 675). T2 reclassified to LIT
  (line 179, with downgrade note). Correctly matches Scope A (mandatory list).
- **New hidden Trues post-amendment?** **NO** — Phase 4: 4 `=True` literals, all either
  (a) loop sentinels falsified inside loops (`all_detector_consistent`), (b) flagged
  STRUCTURAL/DERIVED but NOT in T_pass expression, or (c) explicit T9_pass=True as
  DECLARATIVE. Phase 5: 1 `=True` (T10_pass DECLARATIVE). Audit clean.

### §9.D — Independent count Phase 4+5 (Phase 1-3 from Iter II)

- **Phase 4:** 9 tests | **2 FP** (T1=FP10, T2=FP11) / **6 LIT** (T3-T8) / **1 DEC** (T9).
  Matches self-claim.
- **Phase 5:** 10 tests | **1 FP** (T5=FP13) / **8 LIT** (T1-T4, T6-T9) / **1 DEC** (T10).
  Matches self-claim.
- **Cumulative Phase 1-5 (Iter II Phase 1-3 totals 8 FP / 25 LIT / 3 DEC + Phase 4+5 independent):**
  Phase 1: 4 FP / 8 LIT / 1 DEC; Phase 2: 3 FP / 10 LIT / 1 DEC; Phase 3: 1 FP / 7 LIT / 1 DEC;
  Phase 4: 2 FP / 6 LIT / 1 DEC; Phase 5: 1 FP / 8 LIT / 1 DEC.
  **Totals: 55 tests | 11 FP (20.0%) / 39 LIT (70.9%) / 5 DEC (9.1%) | 90.9% non-trivial.**
  Matches self-report exactly.
  Note: Iter II §D had 8 FP / 25 LIT / 3 DEC for Phase 1-3 = 36 tests. Phase 4+5 add
  3 FP / 14 LIT / 2 DEC = 19 tests. Total 11/39/5. Consistent.
- **Hidden True count:** **0** at substantive level. Only legitimate DECLARATIVE-flagged
  True literals (4 in Phase 3 amendment, 1 in Phase 4 T9, 1 in Phase 5 T10) + 2
  status-flagged STRUCTURAL/DERIVED literals in Phase 4 T5/T6 that are NOT in pass flow.
- **Adversarial vs self-claim delta:** **0.0 pp** for FP count (11 = 11), 0.0 pp for LIT
  (39 = 39), 0.0 pp for DEC (5 = 5). Adversarial recount exactly matches self-claim.

### §9.E — VERDICT

**PASS — closure AUTHORIZED**

Reasoning:
- 0 hidden Trues at substantive level (Iter I had 4; Iter II amendments resolved).
- 0 reclassifications needed (Iter I had 9/36 = 25%; Iter II had 0; Iter III also 0).
- Phase 4 FP10 (rank-1 native Fisher) is substantive: real `eigenvals()`, real `.rank()`,
  alpha emerges from `sp.diff` of `Delta_e_2_native`, 17 sympy sub-checks. NOT postulated.
- Phase 4 FP11 (1D collapse) is substantive: `sp.solve` for native param, hyperplane
  proof via 3 distinct points, chain-rule (16/45) derivation NOT inserted.
- Phase 5 FP13 (network quadrature) is substantive: matrix sum + factorization +
  rank preservation under additivity + eigenvalue scalar additivity → quadrature rule
  DERIVED from Phase 4 rank-1 outer-product structure. NOT postulated.
- No HIGH-severity flags. No new MEDIUM flags. Pre-locked decision rule (HIGH OR ≥3
  reclassifications → AMENDMENT) NOT triggered.

### §9.F — Cycle integrity assessment

- **Substance integrity:** HIGH quality vs cohort 2026-05-11. 20.0% FP cumulative
  is honest (post-amendment), with the 3 NEW FP claims (FP10/FP11/FP13) each carrying
  genuine matrix operations (eigenvals, rank, hyperplane solve, outer-product additivity).
  Phase 4b 3PN extension explicitly deferred (not silently omitted) — anti-drift
  discipline maintained.
- **Audit trail invariant:** PRESERVED. §1-§8 untouched (read-verified line counts match
  prior state). §9 appended only.
- **Anti-drift compliance:** COMPLY. Zero hidden Trues. All DECLARATIVE literals have
  matching `_status` flags. Substance budgets per phase: P4 11.1% DEC (≤15% ✓), P5 10%
  DEC (≤10% ✓), cumulative 9.1% DEC (≤15% ✓), 90.9% non-trivial (≥60% ✓).
- **Pre-registration LOCK (PR-002):** Phase 5 outputs σ_Δφ thresholds in microradians
  at f=100Hz for 4 detector classes (LIGO-O5 / ET-D / CE / ET+CE network) — falsifiable
  observable units consistent z PR-002 spec. M9.1'' Path 2 falsifiability: ET-D/CE/network
  decisive (>5σ single-event), LIGO-O5 borderline (~15σ but below "decisive ≥20σ"). Lock
  intact.

### §9.G — Closure recommendation

**Recommended claim_status: A−** (STRUCTURAL_DERIVED_NATIVE, L2 not-attempted).

Per CYCLE_LIFECYCLE §Claim status taxonomy:
- Substance qualifies as STRUCTURAL_DERIVED_NATIVE: native (a_3, xi_3, c_0·κ_σ) parameter
  space → beta_ppE = (45/16)·Δe_2 LOCK derived from chain rule + rank-1 Fisher; Phase 4b
  3PN deferral is honest scope (rank-1 → rank-2 disambiguation deferred).
- L2 (literature transfer to standard ppE constraint paper) appears NOT-attempted in this
  cycle (Phase 6 BINDING gate still pending, no separate L2 doc spawned). Cycle is
  internally self-contained.
- A+ would require L2 active transfer to a published constraint paper.
- A would require L2 attempted (e.g., draft + failed integration) — not evident.
- **A− is honest: structurally derived native result, audit clean, L2 deferred.**

If Phase 6 BINDING gate triggers L2 spawn post-closure, cycle could be upgraded A− → A
on subsequent re-evaluation. Do NOT inflate to A or A+ pre-L2.

### §9.H — Sign-off

Audit conducted: 2026-05-12 (Iter III final).
Files audited: Phase4_sympy.py (771 lines, full) + Phase5_sympy.py (743 lines, full) +
Phase4/5_sympy.txt run outputs; spot-check Phase{1,2,3}_sympy.py amendment lines via grep.
NOT read: Phase{4,5}_results.md, README.md §7.5e/§7.5f, STATE.md, Phase1-3_amendment_2026-05-12.md.
Audit invariant: §1-§9 IMMUTABLE append-only.

Iter III subagent: independent adversarial verification per CALIBRATION_PROTOCOL §4.4
Phase FINAL pre-closure trigger.
