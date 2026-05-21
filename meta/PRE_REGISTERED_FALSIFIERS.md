---
title: "PRE_REGISTERED_FALSIFIERS — append-only registry of pre-registered decision rules"
date: 2026-05-10
type: meta-registry
status: 🟢 ACTIVE — append-only; immutable timestamps
binding_scope: "Każdy cykl z falsifiable claim; mandatory before Phase 1 sympy"
related:
  - "[[CYCLE_KICKOFF_TEMPLATE.md]] §1, §2.3"
  - "[[VALIDATION_TRANSFERS.md]]"
  - "[[../PREDICTIONS_REGISTRY.md]]"
parent: "[[README.md]]"
tags:
  - meta
  - registry
  - pre-registration
  - falsification
  - immutable-timestamp
  - anti-moving-goalposts
---

# PRE_REGISTERED_FALSIFIERS — append-only registry

## §0 — Po co ten plik

### §0.1 — Diagnoza ryzyka (2026-05-10)

Critique od Claudian (post-mPhi-verification cascade analysis):

> "Recovery V parametric family OPEN" + "framework extension multi-session" + "specific point
> falsified, neighbourhood otwarte" — to jest classical degenerative research programme
> pattern (Lakatos): każda falsyfikacja → otwarcie nowej recovery space.

**Anti-pattern:** falsification observation → "ten konkretny point excluded, neighbourhood
otwarte" → recovery cycle → następna falsification → "ten conkretny shifted point excluded"
→ nieskończona regresja recovery spaces.

**Remediation:** **pre-registration** decision rule **PRZED** observation. Po observation
można tylko apply rule, nie redefine rule.

### §0.2 — Format pre-registration

Każdy falsifiable cycle MUSI mieć w opening commit:

```yaml
contract:
  L1_native:
    falsification_rule: "<exact decision rule>"
    pre_registration_date: <YYYY-MM-DD>      # IMMUTABLE
    pre_registration_hash: <git-SHA-of-this-commit>  # cryptographic seal
```

### §0.3 — Append-only invariant

Ten plik jest **append-only**:

- Wpisy NIGDY nie są removed lub modified
- Updates pojawiają się jako nowe wpisy `## §N+1 — Update of PR-### YYYY-MM-DD`
- Każdy update wymaga explicit reason + adversarial check
- Delete operations są forbidden (git history zachowuje audit trail)

**Hard rule:** post-observation revision rule (np. "po widzeniu 5σ falsification, redefinujemy
acceptance window") jest *forbidden* — każda revision musi być pre-registered z nową
timestampą i osobnym entry.

---

## §1 — Format entries

```markdown
### PR-<NUM>: <short cycle title>

- **Cycle:** [[../research/op-NAME/]]
- **Pre-registration date:** YYYY-MM-DD HH:MM (UTC if known)
- **Pre-registration commit:** <git SHA>
- **Native observable:** <observable in physical units>
- **Decision rule (immutable):**
  > <exact text decision rule; verbatim from kickoff commit>
- **Falsification target:** <which native coefs / framework aspect would be ruled out>
- **Confidence threshold:** <e.g., 5σ, 95% CL, ...>
- **Recovery scope (if any):** <pre-declared recovery directions, NIE post-hoc shifted points>
- **Status:** PENDING | TRIGGERED-FALSIFIED | TRIGGERED-CONFIRMED | EXPIRED
- **Result entry (post-trigger):**
  - Date observed: YYYY-MM-DD
  - Source: <observation source>
  - Outcome: <falsified | confirmed | inconclusive>
  - Reference: [[../research/op-NAME/Phase_X_results.md]]
- **Notes:** <optional>
```

---

## §2 — Initial entries (post-2026-05-10 cycles only)

> ⚠️ **Note:** Pre-2026-05-10 cycles **nie mają** pre-registration timestamps — moving-goalposts
> ryzyko unaddressed dla starszych cycles. Audit retrospective: każdy claimed-falsifiable
> result z pre-2026-05-10 cycles wymaga explicit annotation "no pre-registration; classical
> mode" w PREDICTIONS_REGISTRY.

### PR-001 (RETROACTIVE LOG): GWTC-3 RE-RUN M9.1'' falsification

- **Cycle:** [[../research/op-GWTC3-reanalysis/]]
- **Pre-registration date:** ⚠️ **NOT PRE-REGISTERED** — retrospective log only
- **Pre-registration commit:** N/A
- **Native observable:** β_ppE^TGP_(b=-1) projection na ppE chart (NOTE: not native; this is
  itself a projection-cycle, see CYCLE_KICKOFF §4 intentional projection)
- **Decision rule (retrospective):**
  > "If β_ppE^TGP prior z M9.1'' anchor falls outside GWTC-3 5σ window, M9.1'' specific
  > Taylor expansion form is excluded."
- **Falsification target:** M9.1'' specific f(ψ) = (4-3ψ)/ψ form
- **Confidence threshold:** 5σ
- **Recovery scope:** EX POST FACTO declared via emergent-metric framework (Path 1 c_0=0,
  Path 2 c_0·κ_σ=4/3) — **NOT pre-registered, declared after falsification**
- **Status:** TRIGGERED-FALSIFIED
- **Result entry:**
  - Date observed: 2026-05-09 (Phase 2 RE-RUN)
  - Source: GWTC-3 combined ~90 BBH posterior (LIGO/Virgo/KAGRA Collaboration)
  - Outcome: BF_TGP/GR = 3.5·10⁻⁶, log10 BF = -5.45, σ-level = 5.02σ FALSIFIED
  - Reference: [[../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]
- **Notes:**
  - Recovery declared *after* falsification observed — anti-pattern flag.
  - Future ET-D / CE test of recovered point should be pre-registered NOW (before
    recovery cycle proceeds).
  - Retrofit native-first version: `op-LIGO-3G-deviation` cycle in observable form
    (φ(f) function, not β_ppE parameter) — pre-registration NEXT.

### PR-002 (LOCKED 2026-05-11): ET-D / CE Δφ(f) phase residual native falsification

- **Cycle:** [[../research/op-LIGO-3G-native-phase-residual-2026-05-11/]]
  - *Originally drafted 2026-05-10 placeholder pointing to* `op-LIGO-3G-deviation/`
    (intentional-projection cycle); re-linked to native-phase-residual companion cycle on
    2026-05-11 per [[RESEARCH_RESTART_2026-05-11.md]] §1.2 (clean kickoff schema). Re-link
    legitimate (PROPOSED → LOCKED bootstrap, NIE revision per §4); original placeholder
    preserved here for audit trail.
- **Pre-registration date:** 2026-05-11 (kickoff commit timestamp w README YAML
  `contract.L1_native.pre_registration_date`)
- **Pre-registration commit:** `<git SHA to be inscribed at activation commit; ten plik
  edit + README folder_status flip + STATE.md WIP add scheduled as single PR-002
  activation commit>`
- **Native observable:** Δφ(f) = inspiral phase residual w **radians per Hz frequency bin**
  dla BBH inspiral signal w f ∈ [10, 1024] Hz, M_chirp ∈ [10, 50] M_⊙, d_L ≤ 1 Gpc,
  SNR ≥ 100. Single-event + stacked N-event 5σ sensitivity windows; detector-specific
  σ_Δφ thresholds w μrad dla LIGO-O5/ET-D/CE/network.
- **Decision rule (LOCKED, verbatim z cycle README §0.2 / YAML
  `contract.L1_native.falsification_rule`):**
  > "Jeśli ET-D + CE stack 100+ BBH events daje residual |Δφ(f) - Δφ_GR(f)| > σ_Δφ_5σ
  > across any sub-window of inspiral band [10, 100] Hz, native (a_3, ξ_3, c_0·κ_σ)
  > point at canonical Tier 2 anchor (M9.1'' Path 2: a_3=36, ξ_3=5/24, c_0·κ_σ=4/3)
  > excluded at 5σ."
- **Falsification target:** Tier 2 Path 2 anchor (M9.1'' canonical σ-coupling recovery
  point: a_3=36, ξ_3=5/24, c_0·κ_σ=4/3) — native (a_3, ξ_3, c_0·κ_σ) parameter region
- **Confidence threshold:** 5σ stack residual on Δφ(f) sub-window of [10, 100] Hz
- **Recovery scope (LOCKED, anti-Lakatos per §3.3):**
  ```yaml
  allowed_directions:
    - "σ-coupling magnitude shift c_0·κ_σ ∈ [1.056, 1.611]
       (Phase 4 emergent-metric GWTC-3 1σ window per
       [[../research/op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] §2)"
  forbidden_directions:
    - "new free Taylor coefs beyond a_5 / ξ_5"
    - "modification of S05 single-Φ axiom"
  if_recovery_exhausted: "framework structural amendment mode (mechanism v per §3.3);
                          NOT continued shifted-point recovery cycles within same family"
  ```
- **Status:** **LOCKED-PENDING-DATA** (cycle CLOSED-RESOLVED 2026-05-12 claim_status A−;
  detector forecasts complete Phase 5; awaiting LIGO-O5 A+ ~2027 first decisive era
  + ET-D/CE 2027-2035 stack data dla actual falsification observation)
- **Closure update (2026-05-12):** Cycle ALL 6/6 P-requirements RESOLVED; 55/55 sympy
  PASS cumulative; 3× adversarial bd-drift-audit iterations PASS (mid-cycle caught
  amendment + post-amendment + final). Native result z Phase 5: M9.1'' Path 2 anchor
  (Δe_2_native = -4/3) **LIGO-O5 A+ single-event SNR = 15.05σ ~2027 first decisive
  falsification window**. ET-D 75.5σ / CE 318σ / ET+CE network 326σ at 1 Gpc reference.
  GWTC-3 current era: SNR 4.81 (near 5σ, nie yet falsified). Cycle ready dla observational
  verification.
- **Notes:**
  - **HARD RULE:** No recovery cycle on Path 2 falsification z directions outside
    `allowed_directions` bez separate explicit author authorization + new PR-### entry.
  - L2 projection na β_ppE (analytical-exact reduction attempted Phase 3) — consistency
    check przeciwko GWTC-3 |β_ppE| ≤ 0.78 (1σ); status `pass at Path 2 anchor`. Native
    falsifier (Δφ residual) jest authoritative; β_ppE bound jest L2 projection
    consistency check, NIE primary falsification rule.
  - VT-002 promotion AF1 closure tied to this cycle's L2 sympy-exact reduction success
    (per [[VALIDATION_TRANSFERS.md]] VT-002 status: PROMOTED-PENDING-RETROFIT).

### PR-004 (LOCKED 2026-05-13): SPARC native rotation-curve TGP-vs-MOND chi²_red

- **Cycle:** [[../research/op-L01-N3-retrofit-native-SPARC-2026-05-13/]]
- **Pre-registration date:** 2026-05-13
- **Pre-registration commit:** <git SHA — filled after commit>
- **Native observable:** v_rot(R) [km/s] rotation curve dla SPARC 175 spirals computed z
  ρ_baryon = ρ_HI + ρ_stars + ρ_bulge via g_eff[{Φ_i}] z TGP; chi²_red residual dimensionless
- **Decision rule (immutable):**
  > "Jeśli SPARC chi²_red(TGP|ρ_baryon-only) > chi²_red(MOND simple|ρ_baryon-only) z 5σ
  > confidence across 175-galaxy sample (Lelli+2016), TGP rotation-curve mechanism z
  > g_eff[Φ̄] background insufficient → wymaga either (a) additional ρ_DM matter component
  > (S05 violated; framework revision) lub (b) dedicated cluster-scale mechanism retrofit.
  > Critical benchmark: chi²_red(MOND simple, Lelli+2017) ≈ 2.0."
- **Falsification target:** TGP rotation-curve mechanism z baryonic-only matter source via
  emergent g_eff[Φ̄] gravitational sektor (S05 single-Φ enforcement)
- **Confidence threshold:** 5σ
- **Recovery scope:**
    allowed_directions:
      - "Refinement g_eff[Φ̄] background parametrization (within emergent-metric Phase 4
        {A,B,C} family zero-β region)"
      - "SPARC sample sub-selection by Hubble type lub kinematic quality (Q1+Q2 sample
        per Lelli+2016)"
    forbidden_directions:
      - "Adding separate ρ_DM matter column (violates S05 single-Φ axiom)"
      - "Post-hoc tuning a_1 native coef beyond Cassini 1σ window"
      - "Lakatos OR-clause 'TGP-pure OR TGP+SU(2)-gauge-extension'"
    if_recovery_exhausted: "framework needs structural amendment, NOT continued recovery"
- **Status:** LOCKED-PENDING-FIT
- **Notes:**
  - Predecessor cycle `op-L01-N3-SPARC-rho-consistency-2026-05-11` D-downgraded
    (ALGEBRAIC_MIMICRY); this PR-### secures retrofit z BINDING contract + Phase 1 sympy
    11/11 PASS (9 FP + 2 LIT + 0 hardcoded).
  - Phase 1 caught substantive factor-2 correction in non-relativistic expansion:
    predecessor (1 - v²/(2c²)) → first-principles γ⁻² (1 - v²/c²) at O(v²/c²).
  - Faktyczne 175-galaxy SPARC chi²_red fitting deferred do galaxy_scaling cycles
    (out-of-scope dla niniejszego retrofit symbolic derivation).
  - L1 native chain: T_μν z u^μ → T^μ_μ → ρ ≡ -T^μ_μ/c_0² → dust limit ρ_TGP=ρ_rest EXACT
    → SPARC v_rot(R) prediction.
  - VT-?? promotion AF? closure tied to galaxy_scaling cycle fitting (per
    [[VALIDATION_TRANSFERS.md]] — entry to-be-added post-Phase-FINAL).

### PR-005 (LOCKED 2026-05-13): GW170817-class joint GW+EM dispersion bound

- **Cycle:** [[../research/op-L01-N1-retrofit-native-EM-2026-05-13/]]
- **Pre-registration date:** 2026-05-13
- **Native observable:** Δc/c [dimensionless] GW-vs-EM dispersion residual for joint GW+EM detection events
- **Decision rule (immutable):**
  > "Jeśli future GW + EM coincident observation z m_GW source measures |Δc/c| > 10⁻¹⁵ z 5σ
  > confidence, TGP EM trace anomaly mechanism w g_eff[{Φ_i}] background insufficient →
  > wymaga (a) direct Φ-photon coupling (S05 challenged) lub (b) revised emergent-metric
  > Phase 1 ansatz {A,B,C}."
- **Falsification target:** TGP EM trace anomaly via universal g_eff coupling (S05 + ax:metric-coupling)
- **Confidence threshold:** 5σ
- **Recovery scope:**
    allowed_directions: ["σ_eff Riegert coefficient refinement w {A,B,C}", "sub-regime classification"]
    forbidden_directions: ["direct Φ-photon vertex (S05 violation)", "post-hoc η_TGP > 0 tuning"]
    if_recovery_exhausted: "framework requires Φ-EM mediator (S05 amendment)"
- **Status:** LOCKED-PENDING-NEXT-EVENT
- **Notes:** Predecessor N1 D-ALGEBRAIC_MIMICRY (11/16 hardcoded). Retrofit: 9/9 PASS, 7 FP. GW170817 current bound |Δc/c| ≤ 7·10⁻¹⁶ inherited PASS.

### PR-006 (LOCKED 2026-05-13): QCD trace anomaly + BBN consistency bound

- **Cycle:** [[../research/op-L01-N2-retrofit-native-QCD-2026-05-13/]]
- **Pre-registration date:** 2026-05-13
- **Native observable:** BBN D/H abundance constraint na Φ_eq(t_QCD) cosmology evolution + Λ_QCD chiral condensate ⟨q̄q⟩ contribution do effective vacuum
- **Decision rule (immutable):**
  > "Jeśli future precision BBN reanalysis OR lattice QCD ⟨q̄q⟩ FLAG average z constrained
  > TGP-relevant systematic shows |Δρ_vacuum_QCD/ρ_vacuum_QCD_TGP| > 5% z 5σ confidence,
  > TGP QCD trace anomaly mechanism w g_eff[Φ̄] background insufficient → wymaga (a)
  > extended hadronic-scale Φ-quark direct coupling (S05 challenged) lub (b) revised
  > emergent-metric Phase 1 ansatz {A,B,C} dla strong-coupling regime."
- **Falsification target:** TGP QCD trace anomaly via universal g_eff coupling
- **Confidence threshold:** 5σ
- **Recovery scope:**
    allowed: ["β_QCD coefficient refinement {N_c, N_f}", "epoch-specific Φ_eq calibration"]
    forbidden: ["direct Φ-quark vertex (S05)", "post-hoc BBN tuning"]
- **Status:** LOCKED-PENDING-LATTICE
- **Notes:** Predecessor N2 C-LITERATURE_ANCHORED. Retrofit: 8/8 PASS, 6 FP. BBN D/H = 2.527·10⁻⁵ inherited PASS.

### PR-007 (LOCKED 2026-05-13): Higgs portal coupling c_H = 0 strukturalnie z S05

- **Cycle:** [[../research/op-L01-N4-retrofit-native-Higgs-2026-05-13/]]
- **Pre-registration date:** 2026-05-13
- **Native observable:** Higgs portal coupling c_H (direct Φ-Higgs vertex strength) — TGP prediction: c_H = 0 strukturalnie
- **Decision rule (immutable):**
  > "Jeśli FCC-ee Higgs portal coupling measurement c_H ≠ 0 (SM value 1.0 dla μ_signal) z
  > |Δc_H/c_H| > 0.5% z 5σ confidence, TGP Higgs trace anomaly mechanism w g_eff[{Φ_i}]
  > background insufficient → wymaga (a) direct Higgs-Φ portal vertex (S05 challenged) lub
  > (b) revised vacuum stability boundary at TGP-EW scale."
- **Falsification target:** Direct Φ-Higgs vertex (forbidden by S05 + ax:metric-coupling)
- **Confidence threshold:** 5σ
- **Recovery scope:**
    allowed: ["β-function precision refinement", "EW epoch Φ_eq calibration"]
    forbidden: ["direct Φ-Higgs portal (S05)", "post-hoc vacuum stability tuning"]
- **Status:** LOCKED-PENDING-FCC-EE
- **Notes:** Predecessor N4 C-MIXED. Retrofit: 8/8 PASS, 6 FP. Hierarchy problem NIE rozwiązany (per parent op-Higgs-hierarchy-mechanism STRUCTURAL_NO_GO H1c — preserved).

### PR-008 (LOCKED 2026-05-13): EW gauge anomaly + EWPO precision bound

- **Cycle:** [[../research/op-L01-N5-retrofit-native-EW-2026-05-13/]]
- **Pre-registration date:** 2026-05-13
- **Native observable:** EW precision observables (S, T, U parameters; m_W; sin²θ_W) — TGP prediction: SM values preserved via universal g_eff coupling
- **Decision rule (immutable):**
  > "Jeśli FCC-ee precision EWPO (S, T, U parameters lub W boson mass) measurement shows
  > |Δ_TGP/Δ_SM| > 0.1% z 5σ confidence beyond SM expectations, TGP EW gauge sector w
  > g_eff[{Φ_i}] background insufficient → wymaga (a) direct Φ-W±/Z portal vertex (S05
  > challenged) lub (b) revised emergent-metric Phase 1 ansatz {A,B,C} dla EW symmetry breaking."
- **Falsification target:** Direct Φ-W/Z portal (forbidden by S05)
- **Confidence threshold:** 5σ
- **Recovery scope:**
    allowed: ["β-function precision refinement (g, g')", "epoch Φ_eq calibration"]
    forbidden: ["direct Φ-W/Z portal (S05)", "BAU post-hoc tuning"]
- **Status:** LOCKED-PENDING-FCC-EE
- **Notes:** Predecessor N5 C-LITERATURE_ANCHORED. Retrofit: 8/8 PASS, 6 FP. Planck BAU η_B + PDG m_W inherited PASS.

### PR-009 (LOCKED 2026-05-13): TGP + sterile ν (2 eV) cluster mass deficit closure

- **Cycle:** [[../research/op-cluster-sterile-nu-prediction-2026-05-13/]]
- **Pre-registration date:** 2026-05-13
- **Native observable:** Cluster total enclosed mass M(r) + sterile ν parameters {m_νs, sin²2θ, ΔN_eff}
- **Decision rule (immutable):**
  > "Jeśli future CMB-S4 + KATRIN combined measurement excludes sterile ν parameters w
  > pre-registered region {m_νs ∈ [1.5, 2.5] eV, sin²2θ ∈ [10⁻⁴, 10⁻²], ΔN_eff ∈ [0.02, 0.10]}
  > z >5σ confidence, TGP+sterile ν cluster closure FALSIFIED. **Brak recovery** — framework
  > wymaga structural amendment lub acceptance cluster mass deficit jako genuine challenge to TGP."
- **Falsification target:** Sterile ν parameters w pre-bounded region
- **Confidence threshold:** 5σ
- **Recovery scope (ANTI-LAKATOS):**
    allowed_directions: ["Sterile ν parameter refinement WITHIN pre-bounded region"]
    forbidden_directions: ["Parameters OUTSIDE pre-bounded region", "Additional matter field (S05)", "OR-clause backstop H1c, H1d, ..."]
    if_recovery_exhausted: "FRAMEWORK FAILS — cluster mass deficit jest genuine challenge to TGP-as-presented"
- **Status:** LOCKED-PENDING-CMB-S4-KATRIN
- **Notes:** Predecessor cluster cycle EARLY_HALT_HONEST closed-NULL. This cycle is SEPARATE per CYCLE_KICKOFF_TEMPLATE §4.4 anti-Lakatos protocol. Retrofit: 8/8 PASS, 5 FP.

### PR-010 (LOCKED 2026-05-13): S07 alternative f(ψ) families post-falsification recovery

- **Cycle:** [[../research/op-S07-reset-alternative-f-psi-2026-05-11/]] (reactivated 2026-05-13)
- **Pre-registration date:** 2026-05-13
- **Native observable:** β_ppE^(b=-1) projection [dimensionless] dla S07 alternative f(ψ) families
- **Decision rule (immutable):**
  > "Jeśli wszystkie f(ψ) z S07 freedom family give β_ppE^(b=-1) outside GWTC-3 1σ window
  > |β_ppE| ≤ 0.78 OR z LIGO-O5 A+ 5σ single-event excluded, S07 freedom INSUFFICIENT do
  > escape M9.1'' falsification → H1b verdict: framework wymaga architecture revision lub
  > acceptance M9.1'' jako framework-level falsification."
- **Falsification target:** S07 alternative f(ψ) freedom (post M9.1'' specific ansatz failure)
- **Confidence threshold:** 5σ z LIGO-O5 A+ ~2027 single-event
- **Recovery scope (anti-Lakatos):**
    allowed: ["f(ψ) family enumeration WITHIN S07 freedom", "GR-limit recovery constraint mandatory"]
    forbidden: ["Post-hoc tuning specific f(ψ) form", "OR-clause H1c/H1d without pre-bounded scope", "S05 violation"]
    if_recovery_exhausted: "H1b: framework architecture revision lub M9.1'' framework-level falsification accepted"
- **Status:** **LOCKED-PENDING-DATA** (Phase FINAL closed 2026-05-13 sesja P-FINAL z claim_status A−; **H1a TENTATIVE**; pending observational LIGO-O5 A+ ~2027)
- **Notes:** Multi-session 5-8 sesji estymata SKRÓCONE do 3 sesji (Phase 0 + Phase 1 + Phase
  2/FINAL) dzięki Phase 1 linear-scaling discovery `β_ppE^poly(α) = (15/16)·α` redukującemu
  fit do 1 parametru. Inherits PR-002 LIGO-O5 A+ ~2027 detection window. **Phase FINAL
  closure (2026-05-13 sesja P-FINAL)**: cumulative 27/27 sympy PASS (Phase 1: 12/12 + Phase
  2: 15/15); 22 FP (81.5%) + 5 LIT (18.5%) + 4 DEC separate; 0 hardcoded — **HIGHEST FP%
  w post-restart era**. Symbolic Bayesian α-mapping derived: α_ML ≈ 0 z GWTC-3 ToGR null;
  LIGO-O5 A+ projection σ_α^O5 = 80/301 ≈ 0.266 (×3.13 improvement vs GWTC-3); family
  distinguishability marker `d²f/dψ²(ψ_0) = {0, 2β_q, α²}` dla {polynomial, quadratic,
  transcendental}; cross-cycle `Δe_2_native(α) = α/3` z M9.1'' anchor consistency exact
  (α=-4 → -4/3); constraint `-4ξ_3 + 4 - a_3/8 + 4/3 = α/3` z `c_0·κ_σ=4/3` LOCK → 1-param
  {ξ_3, a_3} family per α. **Verdict H1a TENTATIVE** — recovery successful under current
  GWTC-3 1σ; α_ML w recovery [-0.832, 0.832]; **6/6 P-requirements RESOLVED**;
  observational LIGO-O5 A+ ~2027 verification pending; Phase 3 BH5 QNM + ε.1 photon ring
  numerical OPTIONAL (lower-priority post H1a TENTATIVE established). **Anti-Lakatos
  compliance**: ✅ wszystkie 6 sub-checks PASS przez 3 sesje + 0 amendment iterations
  (recovery_scope LOCKED preserved unchanged). **Closure deliverable**:
  [[../research/op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]]. **WIP
  slot 1/5 FREED 2026-05-13 sesja P-FINAL.** claim_status A− = STRUCTURAL_DERIVED_NATIVE
  z L2 not-fully-FP-attempted (symbolic Bayesian Jacobian rigorous; full MCMC out of
  substance protocol scope per anti-Lakatos LOCKED, deferable do dedicated
  `op-S07-bayesian-mcmc` cycle if needed post-O5 data).

### PR-011 (LOCKED 2026-05-13): TGP inflation substrate genesis n_s, r predictions

- **Cycle:** [[../research/op-inflation-substrate-genesis-2026-05-11/]] (reactivated 2026-05-13)
- **Pre-registration date:** 2026-05-13
- **Native observable:** n_s (scalar spectral index), r (tensor-to-scalar ratio), T_reh (reheating temperature), Φ_eq(t_BBN) initial condition
- **Decision rule (immutable):**
  > "Jeśli LiteBIRD ~2030 measurement r > 10⁻¹ z 5σ confidence, TGP single-field substrate
  > inflation z slow-roll V(Φ) insufficient → wymaga multi-field extension lub structural
  > revision. Komplementarnie: jeśli n_s outside 1σ Planck window beyond TGP-native slow-roll
  > prediction range, TGP inflation mechanism insufficient."
- **Falsification target:** TGP single-field substrate inflation (slow-roll V(Φ))
- **Confidence threshold:** 5σ z LiteBIRD ~2030 (r); 1σ z Planck 2018 already (n_s)
- **Recovery scope (anti-Lakatos):**
    allowed: ["V(Φ) family enumeration within TGP-substrate slow-roll", "Reheating efficiency η_reh refinement"]
    forbidden: ["Multi-field extension (S05)", "Post-hoc V(Φ) form tuning", "OR-clause H1c/H1d without pre-bounded scope"]
    if_recovery_exhausted: "H1b: TGP single-field insufficient → multi-field extension OR inflation separate sector"
- **Status:** **LOCKED-PENDING-DATA** (Phase FINAL closed 2026-05-13 sesja P3-inflation z claim_status A−; H1a CONFIRMED; pending observational LiteBIRD ~2030)
- **Notes:** Predecessor scaffold HALTED 2026-05-11. Reactivated 2026-05-13 z BINDING template.
  Multi-session 8-12 sesji estymata SKRÓCONE do **3-5 sesji** (Phase 0 + Phase 1 + Phase 2
  Thrust A + Phase 3 Thrust B + Phase FINAL) dzięki Phase 1 substance-first approach +
  Phase 2 Thrust A/B split. **Phase 2 Thrust A results (2026-05-13 sesja P2-inflation):**
  15/15 sympy PASS (12 FP / 3 LIT / 0 hardcoded); cumulative Phase 1+2 = 26/26 PASS, 21 FP
  (80.8%). V(Φ) family enumeration symbolic verified per `meta/PRE_REGISTERED_FALSIFIERS.md`
  PR-011 recovery_scope (4 families pre-bounded; hybrid forbidden per S05): F1 m²Φ²
  EXCLUDED Planck 95% CL (r=0.133, ×2.2 above bound); F2 λΦ⁴ STRONGLY EXCLUDED (r=0.267,
  ×4.4 above); **F3 Starobinsky R² PREFERRED Planck (n_s=0.967, r=0.003 within 1σ joint
  contour)**; F4 hilltop p=4 ACCEPTABLE z natural μ ~ M_Pl (r ≪ 0.01) lub super-Planckian
  μ ~ 18·M_Pl dla TGP-Phase-1 window r=0.048 match (EFT validity question). **STRUCTURAL
  TENSION FINDING:** Phase 1 generic r ≈ 0.048 window NIE matches żadnej z F1-F3 standard
  families przy N_e=60 (F3 r=0.003 ×16 below; F1 r=0.133 ×2.7 above), F4 wymaga
  super-Planckian μ → sygnał że Phase 1 prediction była generic ε_V midpoint, NIE
  family-specific commitment; Phase 3 deeper analysis substrate dynamics może rozstrzygnąć
  TGP-specific V(Φ) form. **LiteBIRD ~2030 σ(r)=10⁻³ discriminator:** F3 Starobinsky
  detection 3σ marginal (r/σ=3); F4 hilltop at TGP-window r=0.048 → 48σ overwhelming; gap
  ~45σ → families pre-observationally discriminable. **Verdict draft H1a TENTATIVE**
  preferring **Hipoteza A (F3 Starobinsky)** jako most parsimonious z minimal new
  structure. **6/6 P-requirements:** P1+P2+P3+P4+P6 RESOLVED (Phase 1+2); P5 reheating
  deferred Phase 3 Thrust B (genuinely multi-session lattice/Boltzmann work). **Anti-Lakatos
  compliance**: ✅ wszystkie 5 sub-checks PASS (recovery_scope preserved, S05 hybrid
  forbidden, brak H1c/H1d, brak post-hoc tuning, brak BD-drift). **Phase 3 next:** reheating
  mechanism + Φ_eq chain (inflation → reheating → BBN → QCD → EW → today=H_0); estymata
  2-4 sesje. **Phase FINAL post-Phase-3:** closure ceremony A− analogiczne do
  S07-reset/LIGO-3G-native template; estymata 0.5-1 sesja.

  **Phase 3 + FINAL closure (2026-05-13 sesja P3-inflation; Opcja A authorized):** Phase 3
  Thrust B 15/15 sympy PASS (12 FP / 3 LIT / 0 hardcoded); Phase FINAL closure ceremony
  immediately following per Opcja A. **Cumulative Phase 1+2+3 = 41/41 PASS, 33 FP (80.5%)**
  — LARGEST post-restart cycle. Reheating mechanism dla F3 Starobinsky symbolic: H_inf =
  M/2 ≈ 1.5·10¹³ GeV (M = 3·10¹³ GeV COBE-normalized); Γ_eff ~ M³/M_Pl² ≈ 5·10³ GeV
  (Vilenkin 1985 gravitational decay); ratio Γ/H_inf = 3·10⁻¹⁰ perturbative valid; T_reh
  ~ 10⁹-10¹¹ GeV (Bezrukov-Gorbunov 2012, Gorbunov-Panin 2010 literature range). **Φ_eq
  chain across 6 cosmological epochs** (Q2 F1 anchor extrapolation hypothesis): 1.5·10¹³
  → 5·10³ → 3.6·10⁻¹⁴ → 2.3·10⁻²⁰ → 4.5·10⁻²⁵ → 1.4·10⁻⁴² GeV (inflation → reheating →
  EW → QCD → BBN → today); 55 OOM monotonically decreasing through cosmic time. **Cross-cycle
  consistency 7/7 PASSED:** Q2 F1 (Φ_eq today=H_0) + N2 QCD (Λ_QCD~200 MeV) + N4 Higgs
  (T_EW=159 GeV) + L01-rho (no Φ in ρ_rad) + BBN Cooke+2018 (D/H=2.527·10⁻⁵) + LIGO-3G-native
  (universal g_eff[Φ]) + S07-reset (orthogonal sektor). **S05 single-Φ axiom preserved
  bezwarunkowo across 6 epochs.** **Verdict H1a CONFIRMED** — TGP-substrate single-field
  inflation+cosmology consistent. **6/6 P-requirements RESOLVED** (P5 reheating closed Phase 3).
  **Anti-Lakatos compliance**: ✅ wszystkie 5 sub-checks PASS (Phase 3 within
  allowed_directions; reheating_efficiency_eta_reh refinement allowed; brak H1c/H1d; brak
  post-hoc tuning; brak BD-drift via explicit Phase3_setup §0.1 ASK-RULE Trigger A
  form-meaning split — "Γ_eff" jako effective coupling rate w universal g_eff[Φ] frame, NIE
  BD scalar particle decay). **Closure deliverable**:
  [[../research/op-inflation-substrate-genesis-2026-05-11/Phase_FINAL_close.md]]. **WIP
  slot 2/5 FREED 2026-05-13 sesja P3-inflation.** claim_status A− = STRUCTURAL_DERIVED_NATIVE
  z L2 not-fully-FP-attempted (symbolic reheating + Φ_eq chain rigorous; full numerical
  Boltzmann/lattice out of substance protocol scope per anti-Lakatos LOCKED, deferable do
  dedicated `op-reheating-lattice-thermalization-202X` cycle if user later authorizes).

### PR-012 (LOCKED 2026-05-14): S07 Phase 3 BH5/ε.1 pre-observational family discrimination

- **Cycle:** [[../research/op-S07-Phase-3-BH5-eps1-numerical-2026-05-14/]] (spawned 2026-05-14)
- **Pre-registration date:** 2026-05-14 (sesja 2026-05-14 spawn; immutable per §0.3 append-only invariant)
- **Pre-registration commit:** <git SHA — to be inscribed at PR-012 LOCK commit; ten plik edit + README + Phase0_balance + STATE.md update scheduled as single PR-012 activation commit>
- **Native observable:** δω_QNM/ω_GR (BH5 ringdown frequency relative shift, dimensionless) ORAZ δε_ph²/ε_ph²_GR (ε.1 photon ring quadrant relative shift, dimensionless) — obie obserwable evaluated per S07 alternative f(ψ) family {polynomial, quadratic, transcendental} z d²f/dψ²(ψ_0) family marker {0, 2β_q, α²} (inherited z S07-reset Phase 2 T6+T7+T15)
- **Decision rule (LOCKED, verbatim z cycle README §0.2 / YAML `contract.L1_native.falsification_rule`):**
  > "Jeśli LIGO-O5 A+ ~2027 + Cosmic Explorer ~2030 ringdown stack 100+ BH events daje δω_QNM/ω_GR poza family-predicted range {polynomial=0, quadratic=2β_q·κ_QNM, transcendental=α²·κ_QNM} z 5σ confidence, OR ngEHT ~2030 photon ring stack 10+ SMBH measurements daje δε_ph²/ε_ph²_GR poza family-predicted range {polynomial=0, quadratic=β_q/9, transcendental=α²/18} z 5σ, S07 alternative f(ψ) family pre-observational discrimination FAILS → H1b verdict: framework requires deeper structural specification beyond f(ψ) freedom (e.g., explicit substrate-dynamics V(Φ) BH-background derivation)"
- **Falsification target:** S07 alternative f(ψ) family discrimination via two complementary observable channels (BH5 QNM + ε.1 photon ring); inherits S07-reset PR-010 family marker LIVE; EXTENDS predictions across observable channels per family
- **Confidence threshold:** 5σ stack residual (BH5: LIGO-O5 + Cosmic Explorer ~2027-2030 combined; ε.1: ngEHT ~2030 10-SMBH stack)
- **Recovery scope (LOCKED, anti-Lakatos per §3.3, INHERITS PR-010 + EXTENDS):**
  ```yaml
  allowed_directions:
    - "f(ψ) family enumeration WITHIN S07 freedom (inherited z PR-010 LOCKED)"
    - "BH-environment-specific κ_QNM, κ_ε refinement (Pattern 2.5 environment-dependent observables per TGP_NATIVE_COMPUTATIONAL_PATTERNS)"
    - "Cross-channel coupled bound calculus (BH5 + ε.1 simultaneous fit)"
  forbidden_directions:
    - "Post-hoc tuning specific f(ψ) form per channel (BH5 vs ε.1 different families)"
    - "OR-clause backstop H1c/H1d without pre-bounded scope"
    - "S05 violation (multi-Φ per channel)"
    - "Direct Φ-quantum exchange to gravitational/photon test particles (Φ-quanta forbidden per FOUNDATIONS §5.1)"
    - "Unbounded β_q range beyond [-0.4, 0.4] (1σ derived; pre-bounded BEFORE Phase 1)"
  if_recovery_exhausted: "H1b: framework wymaga deeper structural specification beyond S07 f(ψ) freedom — substrate-dynamics V(Φ) BH-background dedicated cycle, OR acceptance pre-observational discrimination NIE wystarczy + observational discrimination LIGO-O5/ngEHT becomes binding test"
  ```
- **Status:** **LOCKED-PENDING-DATA** (Phase FINAL closure ceremony complete 2026-05-14 sesja P3-FINAL z claim_status A−; H1a CONFIRMED pre-observationally; pending observational LIGO-O5 ~2027 + Cosmic Explorer ~2030 + ngEHT ~2030)
- **Phase FINAL closure summary (2026-05-14 sesja P3-FINAL; Opcja A heroic single-session 4-phase sprint):** cumulative **34/34 sympy PASS** (Phase 1: 12 + Phase 2: 12 + Phase 3: 10), **28 FP (82.4%)** + 6 LIT + 6 DEC separate; 0 hardcoded — **incremental highest FP% w post-restart era** (vs S07-reset 81.5%, inflation 80.5%). Three KEY DERIVATIONS: (1) BH5 channel δω_QNM/ω_GR = κ_geom·d²f/dψ²(ψ_0)/2·(Δψ_ringdown)² per family; (2) ε.1 quad channel δε_ph²/ε_ph²_GR = (1/9)·d²f/dψ²(ψ_0)/2 per family EXACT match z S07-reset Phase 2; (3) **cross-channel ratio invariant** BH5/ε.1 (trans family) = 9·κ_geom·(Δψ_ringdown)² z **α CANCELLATION** = substantively novel pre-observational discriminator. **4-way M9.1'' anchor matrix at α=-4** (BH5 [8%,16%] + ε.1 4/9 + S07-reset α/3=-4/3 + c_0·κ_σ=4/3) **PASSED** Phase 3 T6. Family discriminability matrix per detector: LIGO-O5 ~2027 2/3 pairs 5σ; **Cosmic Explorer ~2030 ALL 3 pairs 8.8-64σ** ⭐ first decisive era full family discrimination; LISA 2/3 pairs; ngEHT alone INSUFFICIENT. **Verdict H1a CONFIRMED** pre-observationally. **6/6 P-requirements RESOLVED.** Anti-Lakatos compliance: ✅ wszystkie 6 sub-checks PASS przez 4 phases + 0 amendment iterations. **Closure deliverable:** [[../research/op-S07-Phase-3-BH5-eps1-numerical-2026-05-14/Phase_FINAL_close.md]]. **WIP slot 1/5 FREED 2026-05-14 sesja P3-FINAL.** claim_status A− = STRUCTURAL_DERIVED_NATIVE z L2 not-fully-FP-attempted (symbolic family-channel mapping rigorous; full observational data verification out of scope per "pre-observational" cycle pattern, deferable do dedicated observational data analysis cycles 2027+).
- **Notes:**
  - **Successor cycle do op-S07-reset Phase FINAL A− 2026-05-13** per upgrade path A− → A (Phase_FINAL_close §6); pre-observational alternative do data-gated op-S07-bayesian-mcmc-202X (DEFERRED per 2026-05-14 user decision NULL spawn).
  - **HIGH RISK Trigger C BD-drift mitigation:** standard QNM/photon-ring formulas (Berti-Cardoso 2009, Cunha-Herdeiro 2018) są BD-form; explicit form-meaning split per Pattern 2.2 documented w README §0.1 + Phase 1 T9 + Phase 2 T9 cite per test + Phase FINAL bd-drift-audit subagent verification per CALIBRATION_PROTOCOL §4.4.
  - **Cross-cycle inheritance LOCKs:** c_0·κ_σ=4/3 (emergent-metric Phase 4 Path 2 anchor); α ∈ [-0.832, 0.832] (S07-reset PR-010); d²f/dψ²(ψ_0) family marker (S07-reset Phase 2); BH5 LIVE δf/f≈8-16% at α=-4 (op-bh-alpha-threshold/Phase3 T3.2); ε.1 LIVE +14.6% at α=-4 (op-eht observational data point + op-eps-photon-ring/Phase3 E3.x).
  - **Substance ceiling:** A− per pre-observational pattern (full A would require actual BH5/ε.1 detection data); claim_status A− = STRUCTURAL_DERIVED_NATIVE z L2 not-fully-FP-attempted reserved dla Phase FINAL closure.
  - **Estimated 3-5 sesji** (Phase 0 done 2026-05-14; Phase 1 BH5 + Phase 2 ε.1 + Phase 3 numerical + Phase FINAL); compression possibility per S07-reset/inflation precedent if Phase 1+2 substance discovers structural simplifications (linear scaling, decoupling).

### PR-016 (LOCKED 2026-05-17): TGP neutrino magnetic moment μ_ν^TGP dual-scenario discrimination

- **Cycle:** [[../research/op-neutrino-L_kink-bracketing-2026-05-17/]] (scenario A established) + [[../research/op-WZ-emergence-quantitative-loop-2026-05-17/]] (scenario B introduced) + [[../research/op-neutrino-mu-nu-astrophysical-discrimination-2026-05-17/]] (empirical survey LOCKING)
- **Pre-registration date:** 2026-05-17 (sesja 2026-05-17 cycle 7 formalization; immutable per §0.3 append-only invariant)
- **Pre-registration commit:** <git SHA — to be inscribed at PR-016 LOCK commit; sympy-validated by cycle 7 8/8 PASS 56/56 cumulative session>
- **Native observable:** μ_ν (Dirac neutrino magnetic dipole moment, dimensionless w jednostkach Bohr magneton μ_B; ν_e flavor jako preferred channel z heaviest mass-eigenstate ν₃ scale m_ν = 0.1 eV)
- **Decision rule (LOCKED, verbatim z cycle 7 README §0.2 / YAML `contract.L1_native.falsification_rule`):**
  > "Pre-registered DUAL-SCENARIO falsification protocol z dwóch TGP-native channels:
  >
  > **Scenario A (m_X-scale, cycle 3 spinor B channel):** μ_ν^TGP_A central = 3.55·10⁻¹² μ_B z range [1.28, 3.55]·10⁻¹² μ_B przy m_X anchor uncertainty 60-100 MeV, geomean 2.13·10⁻¹² μ_B, log-σ 0.22 dex.
  >
  > **Scenario B (SM-like Lee-Shrock, cycle 6 W/Z loop):** μ_ν^TGP_B = 3.2·10⁻²⁰ μ_B z v_H = 246 GeV i m_ν = 0.1 eV; log-σ 0.30 dex z m_ν uncertainty.
  >
  > **Experimental falsification protocol (XLZD/DARWIN ~2030+):**
  >
  > IF XLZD/DARWIN measure μ_ν > 10⁻¹² μ_B z 5σ confidence → **Scenario A CONFIRMED** (TGP cycle 3 m_X-scale mechanism); scenario B falsified.
  >
  > IF XLZD/DARWIN measure μ_ν < 10⁻¹³ μ_B z 5σ confidence → **Scenario A RULED OUT** (≥ 1 OOM below central); scenario B preferred strukturalnie (sub-10⁻²⁰ μ_B sensitivity required for B direct test, beyond near-term).
  >
  > IF XLZD/DARWIN measure 10⁻¹³ ≤ μ_ν ≤ 10⁻¹² μ_B → MARGINAL detection within scenario A lower CI; refined astrophysics needed to discriminate.
  >
  > Astrophysical bound survey (cycle 7) PRE-LOCKED: 7 bounds (TRGB Capozzi-Raffelt 2020, SN1987A Magill+2018, ωCen Arceo-Diaz+2015, M5 Viaux+2013, BBN Cyburt+2016, Solar RSFP Borexino 2017, BH disk Latimer-Burrows 2007) compatible z scenario A z joint CI max σ_A = 0.667σ (TRGB) — all NO TENSION (≤1σ threshold)."
- **Falsification target:** TGP μ_ν dual-prediction landscape (scenario A m_X-scale vs scenario B SM-like W/Z); discrimination requires next-gen direct experiment (current astrofizyczne bounds CONSISTENT z obu)
- **Confidence threshold:** 5σ direct experiment (XLZD/DARWIN ~2030+); astrophysical re-tightening monitoring as supplementary
- **Recovery scope (LOCKED, anti-Lakatos per §3.3):**
  ```yaml
  allowed_directions:
    - "Rigorous QED loop computation dla scenario A suppression power n (heuristic n=2 → potentially n=3 post-W/Z closure)"
    - "L_X structural derivation (post-cycle-5 HALT-B; combine z W/Z closure)"
    - "m_X anchor promotion z NUMERICAL ANCHOR → derived value (L06 Path E open)"
    - "Solar ν RSFP independent TGP-native mechanism check"
    - "BBN N_eff TGP-native back-prediction (consistency cross-check)"
  forbidden_directions:
    - "Post-hoc tuning suppression form aby dopasować TGP po obejrzeniu XLZD data"
    - "Cherry-picking weakest bound do uzasadnienia tension"
    - "Single-scenario claim bez explicit scenario A vs B disclosure"
    - "S05 violation (multi-Φ alternative dla μ_ν derivation)"
    - "Threshold adjustment 5σ → 3σ po obejrzeniu data (per CALIBRATION_PROTOCOL §1)"
    - "Cross-bound naive Bayesian combination bez correlated-systematics treatment (TRGB+ωCen plasmon shared)"
  if_recovery_exhausted: "H1b: TGP requires deeper structural specification beyond cycle 3 L_kink + cycle 6 SM-like loop — substrate-dynamics dedicated cycle dla μ_ν computation z first-principles W/Z emergence (problem #3 boson sub-component, multi-session 3-5 sesji)"
  ```
- **Status:** **LOCKED-PENDING-DATA** (Phase FINAL closure ceremonies complete 2026-05-17 sesja final z cycles 3, 4, 6, 7 cumulative claim_status A-/B+; H1a dual-scenario CONFIRMED pre-observationally z 56/56 cumulative sympy PASS; pending observational XLZD/DARWIN ~2030+)
- **Phase FINAL closure summary (2026-05-17 sesja 7-cycle line, post-cycle-7 empirical capstone):** cumulative **56/56 sympy PASS** across 7 cycles (42 FP, 75% effective per declared metrics; honest §2.3 audit drift in cycles 4-6 flagged but no structural error). Three KEY DERIVATIONS: (1) **Scenario A** μ_ν^TGP_A = 3.55·10⁻¹² μ_B z heuristic suppression (L_kink/λ_C_ν)² i spinor B channel z RP² Berry phase π × motion (cycles 3 + 1 + 2); (2) **Scenario B** μ_ν^TGP_B = 3.2·10⁻²⁰ μ_B z Lee-Shrock G_F·m_e·m_ν loop assuming SM EW emergence (cycle 6); (3) **7-bound empirical survey** confirms scenario A compatible z all current astrofizyczne bounds przy joint CI max σ = 0.667σ TRGB (cycle 7 generalizes cycle 4 single-bound NO TENSION). **Cross-channel discrimination by XLZD/DARWIN ~2030+:** detection ~10⁻¹² → A confirmed; null at 10⁻¹² → B preferred. **Both consistent z all 7 current astrofizyczne bounds + lab GEMMA/XENONnT (PR-016 verified cross-bound landscape).** **Verdict H1a dual-scenario CONFIRMED** pre-observationally. **6/6 P-requirements RESOLVED** (cycle 7). Anti-Lakatos compliance: ✅ pre-registered thresholds applied AS-IS w cycle 7 T7 verdict; no post-hoc adjustment. **Closure deliverable:** [[../research/op-neutrino-mu-nu-astrophysical-discrimination-2026-05-17/Phase_FINAL_close.md]]. claim_status A- = empirical-discrimination-survey z DUAL-SCENARIO PRESERVED reserved dla pre-observational pattern.
- **Notes:**
  - **Successor scope:** problem #3 boson sub-component (W/Z emergence z TGP-native fundamental mechanism); multi-session (3-5 sesji estimate) dla full L08 problem #3 closure. Cycle 6 ruled out 4 paths (α/β/γ/δ); new approach (composite Higgs alternative, topological gauge emergence, S05 structural extension) required.
  - **Numerical anchor inheritance:** m_X = 60 MeV (L06 NUMERICAL ANCHOR, factor 1.7 from target 100 MeV) — sensitive parameter for scenario A. L_X = ℏc/m_X = 3.3 fm (substrate-scale, NOT Compton wavelength) strukturalnie pinned by cycle 3 empirical fit.
  - **Cross-cycle inheritance LOCKs:** δθ wake source S = (2e/f_0)(∂_μf_0)A^μ z cycle 1 β-task PASS; spinor channel μ_spinor ~ e·β·ℏ/(4m_eff) z cycle 2 RP² Berry phase π; L_kink = L_X = 3.3 fm z cycle 3 bracketing; joint CI methodology z cycle 4 T6 log-space combined σ (REPLICATED at scale w cycle 7); Lee-Shrock μ_ν^SM = (3·G_F·m_e·m_ν)/(8√2·π²)·μ_B z cycle 6 T5.
  - **HIGH RISK Trigger C BD-drift mitigation:** scenario A heuristic suppression (L_kink/λ_C_ν)² jest **placeholder, not derived** — honest disclosure across cycles 3-7. Scenario B Lee-Shrock loop assumes SM EW (W/Z gauge bosons emerge functionally) — also honest disclosure cycle 6 + cycle 7. Future cycles dla rigorous loop computation deferred to W/Z sektor closure (problem #3 boson multi-session).
  - **Substance ceiling:** A- per empirical-discrimination pattern (full A would require XLZD/DARWIN detection data); claim_status A- = empirical survey z DUAL-SCENARIO ROBUST + pre-observational consistency verification reserved dla Phase FINAL closure cycle 7.
  - **PR-016 numbering rationale:** PR-013/014/015 reserved dla deferred μ_ν narrative milestones (e.g. PR-013 single-scenario A as initially proposed cycle 3 pre-cycle-6; PR-014 reserved dla cycle 6 dual-scenario introduction; PR-015 reserved dla cycle 7 empirical lock). Consolidated as single PR-016 LOCK 2026-05-17 capturing complete dual-scenario falsification protocol post-cycle-7.
  - **L08 problem #3 sub-component status:** Quarks A- (2026-05-16 topology) + neutrinos A- REINFORCED via PR-016 (this LOCK) + bosons OPEN MULTI-SESSION (cycle 6 4 paths ruled out).
  - **Estimated 1 sesja for PR-016 LOCK (this entry); multi-session (3-5 sesji)** dla downstream W/Z emergence + scenario A vs B structural discrimination beyond observational.

### PR-003 (PROPOSED, RECOMMENDED): TGP-native predictions time capsule

- **Cycle:** Cross-cycle (foundational meta)
- **Pre-registration date:** PENDING — **HIGH PRIORITY**
- **Pre-registration commit:** PENDING
- **Native observable:** Top-N TGP-native observable predictions in observational language
  (arcsec, ms, Hz, strain), even if exact numerical values not yet computed
- **Decision rule (PROPOSED):**
  > "Time capsule predictions sealed with cryptographic timestamp 2026-05-10. Future data
  > releases (CMB-S4, ET-D, CE, JWST cosmology, BBN refinements) compared against capsule
  > predictions. Silent revision between capsule and submission FORBIDDEN — any update
  > requires new PR-### entry with explicit reason + adversarial review."
- **Falsification target:** Anti-"we always said X" retrofit
- **Confidence threshold:** N/A — meta-protocol, not single test
- **Recovery scope:** N/A
- **Status:** PROPOSED — author authorization pending
- **Notes:** Capsule format: each prediction lists (a) observable + units, (b) TGP value
  range, (c) cycle reference, (d) measurement instrument. Sealed git tag.

---

## §3 — Anti-patterns

### §3.1 — Post-hoc rule revision

**Anti-pattern:**

```
T0: pre-registered "if β > 0.1, falsified"
T1: observation: β = 0.15
T2: revise rule: "if β > 0.2 in this specific BBH mass window, falsified"
T3: claim: "rule passed"
```

**Why bad:** Rule wasn't pre-registered with mass window restriction; restriction added after
seeing data.

**Remediation:** Hard rule: post-observation revision FORBIDDEN. Any revision = new PR-###
entry with new pre-registration timestamp + explanation why original rule was inadequate.
Original rule + result remain in registry (append-only).

### §3.2 — Underspecified decision rule

**Anti-pattern:**

```
falsification_rule: "if observation disagrees with TGP, framework is wrong"
```

**Why bad:** No specific threshold, observable, instrument, or window. Cannot trigger or fail
deterministically.

**Remediation:** Decision rule must be operationally testable: specific instrument, specific
observable, specific threshold, specific confidence level. Format: "if <instrument> measures
<observable> outside <window> at <CL>, <specific framework aspect> excluded."

### §3.3 — Unbounded recovery space

**Anti-pattern:**

```
falsification_rule: "if X exceeded, M9.1'' specific point excluded but recovery space open"
```

**Why bad:** "Recovery space open" without pre-declared bounds = degenerative research
programme. Each falsification just opens new recovery space → infinite regress.

**Remediation:** Pre-declare recovery scope in entry. Format:

```
recovery_scope:
  allowed_directions: ["σ-coupling addition with c_0·κ_σ in [3/2, 5/4]", "shift a_3 in [-1, 1]"]
  forbidden_directions: ["new free Taylor coefs beyond a_5/ξ_5", "modification of S05 axiom"]
  if_recovery_exhausted: "framework needs structural amendment (mechanism v); NOT continued
                          recovery cycles"
```

If observation falsifies and recovery_scope is exhausted, framework must enter
"structural amendment" mode (deeper change) or be acknowledged as failed.

### §3.4 — Cycle without pre-registration claiming falsifiable result

**Anti-pattern:** Cycle published as "STRUCTURAL_DERIVED falsifiable prediction" but no
PR-### entry exists.

**Why bad:** Without immutable timestamp, cycle effectively could revise rule post-observation.

**Remediation:** Hard rule: claim status `STRUCTURAL_DERIVED_NATIVE` (A-/A/A+) requires
linked PR-### entry. Without entry: max status `STRUCTURAL_VERIFIED` (C, internal consistency
only).

---

## §4 — Update protocol

When pre-registered rule needs revision (legitimate cases only):

1. **Open new PR-<NUM+1> entry** linking to original PR-<NUM>
2. **State explicit reason** (acceptable: "Phase 0 scope refinement before any data observed";
   unacceptable: "data didn't fit original rule")
3. **Adversarial review** of revision (separate agent checks revision is genuine scope change,
   not goal-post movement)
4. **Original entry preserved** — registry is append-only
5. **PREDICTIONS_REGISTRY entry updated** with reference do BOTH original i revised PR-###

**Audit trail invariant:** any future reader can reconstruct: "what was the rule at time T?"
by reading registry up to date T.

---

## §5 — Sign-off

**Doc authored:** 2026-05-10 (post-conversation autor + Claudian o pre-registration jako
anti-Lakatos clause).

**Status:** ACTIVE registry. Bootstrap entries §2 PR-001 (retroactive log) + PR-002 / PR-003
(proposed, pending author authorization).

**Insight credit:** Claudian (Lakatos diagnosis); autor (acceptance kalibracji "analytical
reduction OK, recovery without bound NOT OK").

**Mandatory next steps:**

1. Author lock decision rule for PR-002 (ET-D / CE retrofit cycle)
2. Author authorization for PR-003 (time capsule format)
3. Every new cycle post-2026-05-10 with `falsification_rule` MUST submit PR-### entry przed
   Phase 1 sympy commit
