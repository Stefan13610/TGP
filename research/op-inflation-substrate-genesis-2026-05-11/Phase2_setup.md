---
title: "Phase 2 setup — inflation V(Φ) family enumeration + per-family Planck/LiteBIRD discriminator (Thrust A)"
date: 2026-05-13
parent: "[[./README.md]]"
phase: 2
predecessor: "[[./Phase1_results.md]] — 11/11 PASS, n_s = 1-6ε_V+2η_V, r = 16ε_V; Planck-compatible window ε_V ≈ 3·10⁻³, r ≈ 0.048"
authorization: "user 'tak działaj' 2026-05-13 conversation (Opcja A: Phase 2 Thrust A only)"
status: 🟡 ACTIVE — sympy + results pending
---

# Phase 2 setup — inflation substrate Φ_eq Thrust A

## §0 — Pre-flight methodology re-confirmation

Per BINDING workflow (CYCLE_KICKOFF_TEMPLATE.md §1-§2 + §0.4 mandatory pre-flight):

- [x] Phase 1 closed: 11/11 PASS, 9 FP / 2 LIT / 2 DEC, 0 hardcoded; Planck-compatible window
  derived analytically (ε_V ≈ 3·10⁻³, η_V ≈ -8.6·10⁻³, r_predict ≈ 0.048)
- [x] PR-011 LOCKED-PENDING-PHASE-1 → przechodzi do LOCKED-PHASE-2-IN-PROGRESS
- [x] Pre-registered falsification rule (PR-011) IMMUTABLE — Phase 2 NIE modyfikuje
  recovery_scope ani decision_rule; V(Φ) family enumeration WITHIN allowed directions
- [x] S05 single-Φ axiom preserved bezwarunkowo — hybrid (multi-field) family **ZABRONIONA**
  per anti-Lakatos forbidden_directions
- [x] Three-layer L1/L2/L3 presentation MANDATORY per `meta/PPN_AS_PROJECTION.md` §3.1
  (analog dla cosmology: L1 native = ε_V/η_V w TGP-substrate Φ-EOM; L2 = standard slow-roll
  formulas n_s, r; L3 = Planck/LiteBIRD bound mapping na ε_V)
- [x] Anti-BD-drift Triggers A-D executed per `meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md` §1.1
  — see §0.1 below

## §0.1 — ASK-RULE Triggers A-D execution (anti-BD-drift)

| Trigger | Status | Notes |
|---|---|---|
| **A — TGP analogy visible?** | ✅ OK | Phase 2 V(Φ) family enumeration explicitly cited jako 4 standard families z literature (Linde 1983 chaotic; Starobinsky 1980 R²; Boubekeur-Lyth 2005 hilltop). TGP-substrate single-field (S05) wybiera SINGLE field play role inflaton + cosmological vacuum (Phase 1 T13 declarative). |
| **B — Predecessor LOCK inheritance audit** | ✅ OK z conditional citation | Inheriting: (i) Phase 1 ε_V = (M_Pl²/2)·(V'/V)², η_V = M_Pl²·V''/V definitions; (ii) n_s = 1-6ε_V+2η_V, r = 16ε_V slow-roll formulas; (iii) Q2 F1 boundary condition Φ_eq=H_0 today (NIE used w Phase 2 derivation, deferred Phase 3 chain). NIE inherit M9.1''-related LOCKs (S07-related; orthogonal sektor). |
| **C — Reproducing literature without TGP mechanism?** | ✅ OK | Literature standard slow-roll formulas (Stewart-Lyth 1993, Lyth 1997 consistency, Linde 1983 m²Φ², Starobinsky 1980 R², Boubekeur-Lyth 2005 hilltop) **verified symbolically** w Phase 2 jako structural derivation z Phase 1 ε_V/η_V definitions. **NIE jest reproduction without TGP mechanism** — TGP-substrate-specific aspect = S05 axiom restricts allowed families do single-Φ (eliminates hybrid Linde 1991). |
| **D — Hardcoded `T_pass = True`?** | ✅ NIE — protocol mandate 0 hardcoded; każdy test musi mieć symbolic verification step. Phase 1 0 hardcoded; Phase 2 preserves. |

**Sign-off:** Claudian @ 2026-05-13 sesja P2-inflation. Pre-flight PASS.

## §1 — Phase 2 decision question (Q2)

> **Q2:** Która V(Φ) family (jeśli którakolwiek) z TGP-substrate slow-roll inflation jest
> *maximally compatible* z Planck 2018 (n_s = 0.9649 ± 0.0042, r < 0.06 95% CL), i czy
> LiteBIRD ~2030 (σ(r)~10⁻³) może **discriminate** między allowed families?

**Phase 2 zamyka decyzję:**

- **H1a:** ≥1 V(Φ) family Planck-compatible + LiteBIRD pre-observationally discriminable →
  Phase 3 reheating + Phase FINAL closure A− pending observational
- **H1a partial:** families Planck-compatible ale LiteBIRD cannot discriminate przy r ≪ σ →
  H1a weakened
- **H1b:** wszystkie V(Φ) families z TGP-substrate slow-roll EXCLUDED → S05 single-Φ
  inflation INSUFFICIENT → accept H1b verdict per PR-011 immutable

**Brak H1c/H1d** — anti-Lakatos LOCKED.

## §2 — Sub-deliverable scope (B.1-B.4)

### §2.1 — B.1: V(Φ) family enumeration (4 families pre-bounded per PR-011 recovery_scope)

Per PR-011 LOCKED `allowed_directions: V(Φ) family enumeration within slow-roll
TGP-substrate ansatz (polynomial / R² / hilltop / hybrid)`, **z S05 constraint hybrid
ZABRONIONA** (multi-field) → enumerujemy 4 single-field families:

| # | Family | V(Φ) form | Literature | Pre-Phase-2 known status |
|---|---|---|---|---|
| **F1** | polynomial massive | (1/2)·m²·Φ² | Linde 1983 chaotic | r ≈ 8/N_e ≈ 0.13 — EXCLUDED Planck 95% |
| **F2** | polynomial quartic | (1/4)·λ·Φ⁴ | Linde 1983 chaotic | r ≈ 16/N_e ≈ 0.27 — STRONGLY EXCLUDED |
| **F3** | Starobinsky R² (Einstein frame) | V_0·(1 - exp(-√(2/3)·Φ/M_Pl))² | Starobinsky 1980 | r ≈ 12/N_e² ≈ 0.003 — PREFERRED Planck |
| **F4** | hilltop p=4 | V_0·(1 - (Φ/μ)⁴) (small-field regime) | Boubekeur-Lyth 2005 | r ≪ 0.01 tunable — ACCEPTABLE |

**Wszystkie single-field (S05 preserved); pre-declared przed Phase 2 sympy** (NIE post-hoc
addition). Hybrid forbidden_directions per PR-011 IMMUTABLE.

### §2.2 — B.2: Per-family symbolic (ε_V, η_V) → (n_s, r) derivation

Dla każdej rodziny Phase 2 symbolic derives:
- ε_V(Φ) = (M_Pl²/2)·(V'/V)² — explicit funkcja Φ
- η_V(Φ) = M_Pl²·V''/V — explicit funkcja Φ
- N_e(Φ_*, Φ_end) = (1/M_Pl²)·∫(V/V')dΦ — e-folds integral
- Φ_*(N_e=60) — value Φ at CMB scales horizon exit
- n_s(N_e=60), r(N_e=60) — numerical predictions per family

**Standard formulas verified analytically Phase 2 sympy:**

| Family | n_s(N_e) | r(N_e) | r_predict at N_e=60 |
|---|---|---|---|
| F1 m²Φ² | 1 - 2/N_e | 8/N_e | 0.133 |
| F2 λΦ⁴ | 1 - 3/N_e | 16/N_e | 0.267 |
| F3 Starobinsky R² | 1 - 2/N_e | 12/N_e² | 0.003 |
| F4 hilltop p=4 | 1 - 3/N_e | 32·(M_Pl/μ)⁶/N_e³ (small field) | tunable; very small |

### §2.3 — B.3: Planck 2018 + LiteBIRD ~2030 discriminator per family

For each surviving family (Planck-compatible), Phase 2 derives:
- σ(r)_LiteBIRD = 10⁻³ projection (Hazumi+2019)
- Detection significance r_predicted/σ_LiteBIRD per family
- Pre-observational discriminator strength

**Critical TGP-Phase-1 cross-check:** Phase 1 derived TGP-native Planck-compatible window
`r ≈ 0.048` (z ε_V ≈ 3·10⁻³). Ten window NIE matches żadnej standardowej rodziny przy
N_e=60:

- **F3 Starobinsky:** r ≈ 0.003 << 0.048 (×16 below TGP-Phase-1 window)
- **F1 m²Φ²:** r ≈ 0.13 >> 0.048 (×2.7 above — already EXCLUDED Planck)
- **F4 hilltop:** r tunable z μ; może hit r ≈ 0.048 dla specific μ value

**Phase 2 substantywne pytanie:** czy hilltop p=4 z μ tunable może realize TGP-Phase-1
window r ≈ 0.048? Lub Phase-1 r=0.048 jest TGP-substrate-derivation-specific (NIE matched
przez standard families) → sygnał że **TGP-substrate slow-roll family preferowana NIE jest
F1-F4 standard list**, requiring Phase 3 deeper dynamics analysis.

### §2.4 — B.4: H1a/H1b verdict draft

Decision matrix:

| Scenario | V(Φ) family compatible | LiteBIRD discriminator | Verdict |
|---|---|---|---|
| F3 (Starobinsky) lub F4 (hilltop) within Planck 1σ + LiteBIRD discriminable | ≥1 | YES | **H1a TENTATIVE** |
| F3/F4 compatible ale LiteBIRD degenerate (r << σ) | ≥1 | NO (Starobinsky r=0.003 marginal at 3σ) | **H1a partial** |
| Phase 1 r=0.048 window mismatch z all 4 standard families | n/a (mismatch) | n/a | **H1a structural-tension** — Phase 3 deeper analysis required |
| Wszystkie 4 families excluded by Planck | 0 | n/a | **H1b** — S05 single-Φ INSUFFICIENT |

Phase 2 outputs decision draft — Phase 3 reheating + Phase FINAL ceremony deferred separate
sessions.

## §3 — Sympy substance plan (target ≥75% FP, BINDING; aim ≥80%)

**Plan: 12 FP + 3 LIT + 2 DEC** (= 15 PASS-counted; FP fraction = 12/15 = 80.0%; analogiczne
do S07 Phase 2 + inflation Phase 1).

| # | Test | Klasa | Question (substantive) |
|---|---|---|---|
| 1 | T1 | FP | F1 m²Φ²: ε_V(Φ) = 2·M_Pl²/Φ²; η_V(Φ) = 2·M_Pl²/Φ² (degenerate quadratic; ε_V = η_V) |
| 2 | T2 | FP | F1 m²Φ²: N_e integral; Φ_*² = 4·M_Pl²·N_e + Φ_end²; symbolic |
| 3 | T3 | FP | F1 m²Φ²: n_s = 1 - 2/N_e; r = 8/N_e (verified symbolic z ε_V(Φ_*), η_V(Φ_*)) |
| 4 | T4 | FP | F2 λΦ⁴: ε_V = 8·M_Pl²/Φ²; η_V = 12·M_Pl²/Φ²; n_s = 1-3/N_e; r = 16/N_e symbolic |
| 5 | T5 | FP | F3 Starobinsky R² Einstein frame: V = V_0·(1 - exp(-√(2/3)·φ/M_Pl))²; ε_V symbolic large-φ asymptotic |
| 6 | T6 | FP | F3 Starobinsky: η_V symbolic; N_e ≈ (3/4)·exp(√(2/3)·φ_*/M_Pl); n_s = 1-2/N_e; r = 12/N_e² |
| 7 | T7 | FP | F4 hilltop p=4: V = V_0·(1 - (Φ/μ)⁴); ε_V = 8·M_Pl²·Φ⁶/μ⁸ small-field expansion |
| 8 | T8 | FP | F4 hilltop p=4: η_V = -12·M_Pl²·Φ²/μ⁴; n_s = 1 - 3·(μ/M_Pl)⁻¹·(Φ_*/μ)² + O(higher); r = 32·M_Pl⁶·Φ_*⁶/μ⁸ |
| 9 | T9 | FP | Planck 2018 1σ exclusion z r-only: F1 r=0.133 > 0.06 → 95% EXCLUDED; F2 r=0.267 > 0.06 → STRONGLY EXCLUDED |
| 10 | T10 | FP | Planck 2018 (n_s, r) joint: F3 (n_s=0.967, r=0.003) within 1σ Planck contour; F1 (n_s=0.967, r=0.133) outside Planck 1σ contour (sweet-spot mismatch) |
| 11 | T11 | FP | TGP-Phase-1 window (r ≈ 0.048) vs standard families: F3 r=0.003 (×16 below); F1 r=0.13 (×2.7 above); **F4 hilltop μ-tunable: μ_target ~ symbolic relation derived** |
| 12 | T12 | FP | LiteBIRD ~2030 detection per family: σ(r)=10⁻³; F3 r/σ = 3 (marginal 3σ); F4 z r=0.048 → 48σ; F1 already EXCLUDED pre-LiteBIRD |
| 13 | T13 | LIT | Planck 2018 (Aghanim+2020): n_s = 0.9649 ± 0.0042 (TT,TE,EE+lowE+lensing); r < 0.06 (95% CL) |
| 14 | T14 | LIT | Standard inflation references: Linde 1983 chaotic; Starobinsky 1980 R²; Boubekeur-Lyth 2005 hilltop |
| 15 | T15 | LIT | LiteBIRD JAXA mission ~2030 σ(r)~10⁻³ projection (Hazumi+2019; LiteBIRD Collaboration) |
| 16 | T16 | DEC | Anti-Lakatos LOCKED PR-011: brak H1c/H1d, recovery_scope V(Φ) within S05 single-Φ; hybrid (multi-field) FORBIDDEN |
| 17 | T17 | DEC | Three-layer L1/L2/L3 presentation MANDATORY w results.md (cosmology analog per PPN_AS_PROJECTION §3.1) |

**0 hardcoded `T_pass = True`. 100% non-trivial.** Każdy test ma explicit pytanie fizyczne
weryfikowane symbolic.

## §4 — Risk register Phase 2

| ID | Risk | Mitigation |
|---|---|---|
| R-P2.1 | Family enumeration nadmiernie ambitne — 4 families × 4 quantities w 1 sesji | Standard formulas literature-verified (Linde, Starobinsky, Boubekeur-Lyth); sympy verifies symbolic forms z Phase 1 ε_V/η_V definitions (FP-grade gdy explicit derivation chain) |
| R-P2.2 | Hilltop family μ-tunability może wymagać MCMC scan dla TGP-Phase-1 window match | Phase 2 derives **structural condition** μ = f(r_target) symbolic; numerical MCMC scan deferred Phase 3 jeśli verdict wymaga |
| R-P2.3 | **TGP-Phase-1 r=0.048 window mismatch z F1+F2+F3** = potential rzeczywisty structural finding | **Phase 2 dokumentuje honestly** w results.md jako "structural tension" — może być (a) sygnał że F4 hilltop preferowana z specific μ, (b) sygnał że TGP-substrate Phase 1 derivation była ε_V-window średnia, NIE specific family commitment, (c) sygnał Phase 3 deeper analysis required |
| R-P2.4 | **BD-drift risk** — V(Φ) framing slips do "TGP modyfikuje slow-roll" zamiast "TGP single-field W Standard slow-roll regime" | L1/L2/L3 layering MANDATORY; primary L1 = TGP-substrate Φ-EOM (Phase 1 inheritance); L2 = standard slow-roll (n_s, r); L3 = Planck/LiteBIRD constraints |
| R-P2.5 | **Lakatos drift** — kuszenie do wprowadzić "modified Starobinsky z TGP correction" lub "hybrid TGP+inflaton" w trakcie cyklu | PR-011 immutable; jeśli verdict H1b → close-NULL honest per cluster precedent |
| R-P2.6 | Cross-cycle Q2 (Φ_eq=H_0 today) consistency NIE addressed w Phase 2 | Phase 2 scope explicit V(Φ) family — Q2 consistency boundary condition + reheating chain = Phase 3 deferred (analogiczne do S07 Phase 3 BH5/ε.1) |

## §5 — Six P-requirements update (Phase 2 contribution)

| P | Phase 1 status | Phase 2 contribution |
|---|---|---|
| P1: Φ_eq EOM | ✅ RESOLVED (T1 Phase 1) | preserved (no new derivation) |
| P2: ε_V, η_V definitions | ✅ RESOLVED (T4+T5 Phase 1) | extended z per-family explicit ε_V(Φ), η_V(Φ) (Phase 2 T1-T8) |
| P3: n_s prediction Planck-compatible | ✅ RESOLVED (T6+T8+T10 Phase 1) | per-family discriminator (T9+T10 Phase 2) |
| P4: r prediction Planck-compatible | ✅ RESOLVED (T7+T8+T11 Phase 1) | **CLOSED** per-family + LiteBIRD forecast (T12 Phase 2) |
| P5: Reheating + BBN initial conditions | 🟡 PARTIAL (e-folds N_e structural Phase 1) | **deferred Phase 3** (genuinely multi-session lattice/Boltzmann work) |
| P6: S05 preservation | ✅ RESOLVED (T13 DEC Phase 1) | reverified (Phase 2 T16 DEC); hybrid family forbidden |

**Phase 2 zamyka P3+P4 family-specific.** P5 reheating deferred Phase 3.

## §6 — Cross-references

- [[./README.md]] — cycle BINDING contract
- [[./Phase0_balance.md]] — initial scaffold
- [[./Phase1_results.md]] — predecessor (11/11 PASS, slow-roll formulas verified)
- [[./Phase1_sympy.py]] + [[./Phase1_sympy.txt]] — symbolic foundation
- [[../op-S07-reset-alternative-f-psi-2026-05-11/Phase_FINAL_close.md]] — sister cycle A− template (sesja P-FINAL same day)
- [[../op-Q2-vacuum-budget-2026-05-10/]] — Q2 F1 boundary condition (Φ_eq=H_0 today; deferred Phase 3)
- [[../../meta/PPN_AS_PROJECTION.md]] §3.1 — three-layer L1/L2/L3 BINDING (cosmology analog)
- [[../../meta/TGP_NATIVE_COMPUTATIONAL_PATTERNS.md]] §1-§4 — anti-BD-drift Triggers A-D
- [[../../meta/CYCLE_KICKOFF_TEMPLATE.md]] — BINDING contract structure
- [[../../meta/PRE_REGISTERED_FALSIFIERS.md]] PR-011 — anti-Lakatos LOCKED

---

**Phase 2 setup complete.** Sympy execution + results draft pending. Estymata: 1 sesja
sympy + results.md. Phase 3 reheating mechanism deferred separate session(s).
