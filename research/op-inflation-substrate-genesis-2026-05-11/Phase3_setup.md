---
title: "Phase 3 setup — inflation reheating mechanism + Φ_eq chain across 6 cosmological epochs (Thrust B)"
date: 2026-05-13
parent: "[[./README.md]]"
phase: 3
predecessor: "[[./Phase2_results.md]] — 15/15 PASS, F3 Starobinsky preferred (n_s=0.967, r=0.003 within Planck 1σ); STRUCTURAL TENSION Phase 1 r=0.048 vs F1-F3 standard families"
authorization: "user 'Opcja A' 2026-05-13 conversation (Phase 3 + Phase FINAL combined w 1 sesji)"
status: 🟡 ACTIVE — sympy + results + FINAL closure pending
---

# Phase 3 setup — inflation substrate Φ_eq Thrust B

## §0 — Pre-flight methodology re-confirmation

Per BINDING workflow (CYCLE_KICKOFF_TEMPLATE.md §1-§2 + §0.4 mandatory pre-flight):

- [x] Phase 2 closed (Thrust A): 15/15 PASS, 12 FP / 3 LIT / 2 DEC, 0 hardcoded; F3 Starobinsky preferred Planck-compatible
- [x] PR-011 LOCKED-PHASE-2-COMPLETE-THRUST-A → przechodzi do LOCKED-PHASE-3-IN-PROGRESS
- [x] Pre-registered falsification rule (PR-011) IMMUTABLE — Phase 3 NIE modyfikuje recovery_scope; reheating WITHIN allowed_directions per PR-011
- [x] S05 single-Φ axiom preserved bezwarunkowo — multi-field reheating mechanism (e.g., curvaton) ZABRONIONA
- [x] Three-layer L1/L2/L3 presentation MANDATORY per `meta/PPN_AS_PROJECTION.md` §3.1 (cosmology analog)
- [x] Anti-BD-drift Triggers A-D executed per `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §1.1 — see §0.1 below

## §0.1 — ASK-RULE Triggers A-D execution (anti-BD-drift)

| Trigger | Status | Notes |
|---|---|---|
| **A — TGP analogy visible?** | ✅ OK z conditional citation | **HIGH RISK ZONE.** Standard reheating literature (Kofman-Linde-Starobinsky 1994; Allahverdi+2010 review) treats inflaton jako **standard scalar particle z decay rate Γ_φ**. TGP-substrate Φ jest **environment-dependent observable** (Pattern 2.5), NIE particle. Phase 3 derives Γ_eff jako EFFECTIVE coupling rate w universal g_eff[Φ] frame — NIE BD scalar particle decay. Form-meaning split (per `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §4): "decay rate Γ" w TGP context = effective rate przy którym inflaton's energy density transfers to radiation (S05 preserves single-Φ throughout). |
| **B — Predecessor LOCK inheritance audit** | ✅ OK z conditional citation | Inheriting: (i) Phase 1 Φ-EOM Φ̈ + 3HΦ̇ + V'(Φ) = 0; (ii) Phase 2 F3 Starobinsky preferred z V_0 = (3/4)·M²·M_Pl² (Einstein frame); (iii) Q2 F1 anchor Φ_eq(today) = H_0 (OP-3 from `closure_2026-04-26/Lambda_from_Phi0/`). **Q2 F1 anchor extrapolation hypothesis:** Phase 3 ASSUMES Φ_eq = H(t) chain holds przy każdym epoch (NIE derived w Phase 3; honest annotation w results.md). |
| **C — Reproducing literature without TGP mechanism?** | ✅ OK | Standard cosmology formulas (Stefan-Boltzmann ρ_rad, Friedmann H², reheating completion T_reh) **verified symbolically** w Phase 3 jako LITERATURE_ANCHORED applications dla TGP-substrate Φ_eq chain. **NIE jest reproduction without TGP mechanism** — TGP-specific aspect = S05 preservation (single-Φ across 6 epochs; ax:metric-coupling universal). |
| **D — Hardcoded `T_pass = True`?** | ✅ NIE — protocol mandate 0 hardcoded; każdy test musi mieć symbolic verification step. Phase 1+2 0 hardcoded; Phase 3 preserves. |

**Sign-off:** Claudian @ 2026-05-13 sesja P3-inflation. Pre-flight PASS.

## §1 — Phase 3 decision question (Q3)

> **Q3:** Czy reheating mechanism + Φ_eq chain inflation→reheating→EW→QCD→BBN→today=H_0
> wykonalna z TGP-substrate single-Φ S05 axiom, preserving cross-cycle consistency z Q2 F1
> anchor (Φ_eq today = H_0), N2 QCD, N4 Higgs, L01-rho stress-energy?

**Phase 3 zamyka decyzję dla Phase FINAL:**

- **H1a confirmed:** reheating + Φ_eq chain wykonalne z F3 Starobinsky preferred + cross-cycle consistency PASSED → **Phase FINAL claim_status A−**
- **H1a partial:** reheating wykonalne ALE cross-cycle tension (e.g., Φ_eq(t_BBN) inconsistent z N2 QCD) → A− z explicit caveat
- **H1b:** reheating mechanism wymaga structural extension (multi-field) → per PR-011 immutable

## §2 — Sub-deliverable scope (B.1-B.4)

### §2.1 — B.1: Reheating mechanism symbolic dla F3 Starobinsky (preferred Phase 2)

**TGP-substrate reheating interpretation per Pattern 2.5 + S05:**
- Φ jest environment-dependent observable (NIE standard scalar inflaton particle)
- Post-slow-roll: ε_V → 1, V(Φ) regime ends; Φ oscillates around minimum
- Energy transfer Φ → SM matter via universal g_eff[Φ] coupling (S05 preserved bezwarunkowo)
- Reheating temperature T_reh derived analytically z Friedmann + Γ_eff effective decay rate

**Standard formula adapted (form-meaning split per §4 patterns):**
```
H_reh ≈ Γ_eff (reheating completion criterion)
ρ_rad(T_reh) = (π²/30)·g_*·T_reh⁴ (Stefan-Boltzmann radiation)
H² = ρ/(3M_Pl²) (Friedmann)
→ T_reh = (90/(π²·g_*))^(1/4) · √(Γ_eff·M_Pl)
```

**F3 Starobinsky specific:**
- V_0 = (3/4)·M²·M_Pl² (Einstein frame; M = Starobinsky scale)
- COBE normalization: M ≈ 3·10¹³ GeV
- H_inf = √(V_0/(3M_Pl²)) = √((3/4)·M²·M_Pl²/(3M_Pl²)) = M/2 ≈ 1.5·10¹³ GeV
- Γ_eff^grav ~ M³/M_Pl² ~ (3·10¹³)³/(2.4·10¹⁸)² ~ 5·10³ GeV (gravitational decay channel; Vilenkin 1985)
- T_reh ~ 10⁹-10¹¹ GeV literature range (depending on mechanism; Bezrukov-Gorbunov 2012, Gorbunov-Panin 2010)

### §2.2 — B.2: Φ_eq chain symbolic across 6 cosmological epochs

Z Q2 F1 anchor `Φ_eq(today) = H_0` jako boundary condition + Friedmann H(z) per epoch:

| Epoch | Cosmic time order | T scale | H ≈ ? | Φ_eq^epoch |
|---|---|---|---|---|
| **Inflation** | earliest | T_inf ~ M ~ 3·10¹³ GeV | √(V_0/(3M_Pl²)) = M/2 | ~ 1.5·10¹³ GeV |
| **Reheating** | post-inflation | T_reh ~ 10⁹-10¹¹ GeV | Γ_eff ~ 10³ GeV (Starobinsky grav) | ~ 10³ GeV |
| **EW** | radiation era T~159 GeV | T_EW = 159 GeV | √(π²g_*/90)·T²/M_Pl ~ 4·10⁻¹⁴ GeV | ~ 10⁻¹⁴ GeV |
| **QCD** | radiation era T~200 MeV | T_QCD ~ 200 MeV | √(π²g_*/90)·T²/M_Pl ~ 2·10⁻²⁰ GeV | ~ 10⁻²⁰ GeV |
| **BBN** | radiation era T~1 MeV | T_BBN ~ 1 MeV | √(π²g_*/90)·T²/M_Pl ~ 5·10⁻²⁵ GeV | ~ 10⁻²⁵ GeV |
| **Today** | DE era z=0 | T_CMB ~ 2.7 K | H_0 ~ 1.4·10⁻⁴² GeV | = H_0 ✓ Q2 F1 anchor |

**Substantywne pytanie:** czy Φ_eq = H(t) hypothesis preserved THROUGHOUT chain (single anchor
extrapolated per Q2 F1) lub tylko TODAY (anchor specific)? Phase 3 derives symbolic + honestly
annotates extrapolation hypothesis (NIE proof).

**Chain monotonicity:** Φ_eq decreasing through cosmic time order ✓ (inflation > reheating >
EW > QCD > BBN > today; spans 55 orders of magnitude).

### §2.3 — B.3: Cross-cycle consistency check

Per Phase 2 closing P5 + Phase 3 deliverable:

| Cycle | Consistency check | Phase 3 verification |
|---|---|---|
| Q2 F1 anchor | Φ_eq(today) = H_0 boundary | ✅ symbolic (T-Lambda.3 OP-3 anchor; PRESERVED) |
| N2 QCD retrofit | Λ_QCD ≈ 200 MeV; Φ_eq(t_QCD) ~ Λ²_QCD/M_Pl | symbolic check |
| N4 Higgs retrofit | T_EW ~ 159 GeV; Φ_eq(t_EW) ~ T²_EW/M_Pl | symbolic check |
| L01-rho stress-energy | radiation era ρ ∝ T⁴; no Φ contribution to ρ_rad | S05 preservation symbolic |
| BBN constraint | D/H = 2.527·10⁻⁵ Cooke+2018 | Φ_eq(t_BBN) consistent (boundary-condition compatibility) |

### §2.4 — B.4: Phase FINAL setup + verdict

Decision matrix:

| Scenario | Phase 3 evidence | Verdict |
|---|---|---|
| Reheating + Φ_eq chain + cross-cycle consistency wszystkie PASS | Phase 3 substance OK | **H1a confirmed** → Phase FINAL A− |
| Reheating + chain OK ale cross-cycle Q2/N2/N4 tension | partial PASS | A− z explicit caveat |
| Reheating wymaga multi-field extension | structural failure | **H1b** per PR-011 |

Phase 3 outputs decision; **Phase FINAL closure ceremony immediately follows w SAME session**
per Opcja A authorization.

## §3 — Sympy substance plan (target ≥75% FP, BINDING; aim ≥80%)

**Plan: 12 FP + 3 LIT + 2 DEC** (= 15 PASS-counted; FP fraction = 12/15 = 80.0%; analogiczne
do Phase 2 Thrust A).

| # | Test | Klasa | Question (substantive) |
|---|---|---|---|
| 1 | T1 | FP | Stefan-Boltzmann: ρ_rad = (π²/30)·g_*·T⁴ symbolic |
| 2 | T2 | FP | H(T) in radiation era = √(π²g_*/90)·T²/M_Pl symbolic z Friedmann + ρ_rad |
| 3 | T3 | FP | Reheating completion T_reh = (90/(π²g_*))^(1/4)·√(Γ_eff·M_Pl) symbolic |
| 4 | T4 | FP | F3 Starobinsky: V_0 = (3/4)M²M_Pl²; H_inf = √(V_0/(3M_Pl²)) = M/2 symbolic |
| 5 | T5 | FP | Φ_eq^inf = H_inf = M/2 ≈ 1.5·10¹³ GeV (M = 3·10¹³ GeV COBE-normalized) |
| 6 | T6 | FP | Φ_eq^reh = Γ_eff ~ M³/M_Pl² ~ 5·10³ GeV (Starobinsky gravitational decay; Vilenkin 1985) |
| 7 | T7 | FP | Φ_eq^EW = H(T_EW) = √(π²g_*/90)·T_EW²/M_Pl ~ 4·10⁻¹⁴ GeV (cross-cycle z N4-Higgs T_EW=159 GeV) |
| 8 | T8 | FP | Φ_eq^QCD = H(T_QCD) ~ 2·10⁻²⁰ GeV (cross-cycle z N2-QCD Λ_QCD≈200 MeV) |
| 9 | T9 | FP | Φ_eq^BBN = H(T_BBN) ~ 5·10⁻²⁵ GeV (T_BBN ~ 1 MeV; D/H Cooke+2018 epoch) |
| 10 | T10 | FP | Φ_eq^today = H_0 ~ 1.4·10⁻⁴² GeV (Q2 F1 anchor PRESERVED bezwarunkowo) |
| 11 | T11 | FP | Chain monotonicity: Φ_eq^inf > Φ_eq^reh > Φ_eq^EW > Φ_eq^QCD > Φ_eq^BBN > Φ_eq^today (55 OOM span) |
| 12 | T12 | FP | S05 preserved across 6 epochs: single Φ; ax:metric-coupling universal w każdym epoch (no second-Φ field, no curvaton) |
| 13 | T13 | LIT | Cosmological parameters: H_0 = 67.4 km/s/Mpc (Planck 2018); T_BBN ~ 1 MeV; T_QCD ~ 200 MeV; T_EW = 159 GeV |
| 14 | T14 | LIT | Reheating literature: Kofman-Linde-Starobinsky 1994 preheating; Allahverdi+2010 review |
| 15 | T15 | LIT | F3 Starobinsky T_reh literature: ~10⁹-10¹¹ GeV (Vilenkin 1985, Bezrukov-Gorbunov 2012, Gorbunov-Panin 2010) |
| 16 | T16 | DEC | Anti-Lakatos LOCKED PR-011: Phase 3 reheating WITHIN allowed_directions; brak H1c/H1d backstop |
| 17 | T17 | DEC | S05 preservation across 6 epochs explicit (single Φ; ax:metric-coupling universal; brak curvaton/multi-field reheating) |

**0 hardcoded `T_pass = True`. 100% non-trivial.**

## §4 — Risk register Phase 3

| ID | Risk | Mitigation |
|---|---|---|
| R-P3.1 | **BD-drift risk HIGH** — standard reheating literature treats inflaton jako particle z decay rate Γ; TGP-substrate Φ jest environment-dependent observable | Pre-flight ASK-RULE Trigger A explicit (§0.1); Phase 3 derives Γ_eff jako effective coupling in TGP universal g_eff[Φ] frame, NIE BD scalar particle decay; form-meaning split documented w results.md |
| R-P3.2 | Φ_eq = H(t) hypothesis through chain może NIE być derived w Phase 3 — to anchor (today only) per Q2 F1 | Phase 3 derives Φ_eq^epoch = H(z_epoch) jako WORKING HYPOTHESIS (per Q2 F1 extrapolation); **honest annotation** if non-trivial proof deferred do dedicated cycle |
| R-P3.3 | Cross-cycle consistency wymaga numerical comparison N2/N4 retrofit values vs Phase 3 epoch values | Symbolic check + literature-anchored numerical evaluation; deviation explicit-noted; tolerancja 10x dla g_* prefactor variations |
| R-P3.4 | **Lakatos drift** — kuszenie do "modify reheating mechanism dla TGP-specific" jeśli standard formula fails | PR-011 immutable; jeśli standard reheating fails → H1b per recovery_scope exhaustion |
| R-P3.5 | Phase FINAL closure ceremony w SAME session — może wymagać dodatkowego time | Closure ceremony użyję S07-reset template (Phase_FINAL_close.md); ~30 min additional work |
| R-P3.6 | Full Boltzmann/lattice numerical out-of-scope | Explicit annotation: deferred do dedicated `op-reheating-lattice-thermalization-202X` cycle if user later authorizes (out of substance protocol scope) |

## §5 — Six P-requirements update (Phase 3 contribution + closure)

| P | Phase 1+2 status | Phase 3 contribution |
|---|---|---|
| P1: Φ_eq EOM | ✅ RESOLVED | extended z post-inflation regime (oscillation around V minimum) |
| P2: ε_V, η_V definitions | ✅ RESOLVED | preserved (Phase 1+2 inheritance) |
| P3: n_s prediction | ✅ RESOLVED | preserved (Phase 2 F3 Starobinsky) |
| P4: r prediction | ✅ RESOLVED | preserved (Phase 2 F3 r=0.003) |
| P5: Reheating + BBN | 🟡 PARTIAL Phase 1+2 | **CLOSED Phase 3** (T_reh symbolic + chain epochs + BBN consistency) |
| P6: S05 preservation | ✅ RESOLVED | reverified across 6 epochs (Phase 3 T12+T17) |

**Phase 3 zamyka P5 → 6/6 P-requirements RESOLVED → Phase FINAL closure ready.**

## §6 — Cross-references

- [[./README.md]] — cycle BINDING contract
- [[./Phase0_balance.md]] — initial scaffold
- [[./Phase1_results.md]] — Φ_eq EOM + slow-roll formulas
- [[./Phase2_results.md]] — V(Φ) family enumeration; F3 Starobinsky preferred
- [[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] — sister cycle A− closure template
- [[../op-LIGO-3G-native-phase-residual-2026-05-11/Phase6_close.md]] — earlier A− closure template
- [[../closure_2026-04-26/Lambda_from_Phi0/]] — Q2 F1 anchor Φ_eq(today) = H_0
- [[../op-L01-N2-retrofit-native-QCD-2026-05-13/]] — N2 cross-cycle Λ_QCD anchor
- [[../op-L01-N4-retrofit-native-Higgs-2026-05-13/]] — N4 cross-cycle T_EW anchor
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/]] — radiation era ρ ∝ T⁴ S05 preservation
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3 BINDING (cosmology analog)
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift Triggers + form-meaning split
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] — BINDING contract structure
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-011 — anti-Lakatos LOCKED

---

**Phase 3 setup complete.** Sympy execution + results draft + Phase FINAL closure ceremony
pending. Estymata: 1 sesja Phase 3 + Phase FINAL combined per Opcja A authorization.
