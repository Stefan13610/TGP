---
title: "op-CE-H-two-particle-equilibrium — Toy test of CE-H structural claim (Poziom β)"
type: research_cycle
status: PRE_REGISTERED_LOCKED
folder_status: active
pre_registration_date: 2026-05-21
parent_concept_paper: meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md
parent_cycle: op-FFS-quark-object-2026-05-20 (A- conditional closure 2026-05-20)
hypothesis_codes:
  - CE-H (Cosmic Equilibrium Hypothesis)
  - RE-H (Relational Energy Hypothesis)
test_level: TOY_MODEL (simplest possible test)
test_scope: 2-particle equilibrium in TGP Phi-substrate
methodology: native-equations FIRST (NIE derivation z ΛCDM/QCD)
claim_status_target: A | A- (depending on Phase 1a/1b dichotomy outcome)
authorization_chain:
  - "2026-05-21: concept paper Poziom α LOCKED"
  - "2026-05-21: explicit 'działaj' authorization by user for Poziom β"
discipline:
  - anti-Lakatos LOCKED
  - strict cycle 1/2/7 (0 hardcoded FP T_pass=True)
  - max 1 DEC budget per cycle
  - R1+R2+R3 two-tier discipline
  - pre-rejestracja PRZED any sympy
falsifiers_pre_registered: F-β-1 do F-β-5 (LOCKED 2026-05-21, see §3)
---

# op-CE-H-two-particle-equilibrium — BINDING contract

**Pre-registration date:** 2026-05-21
**Status:** LOCKED — żadnych modyfikacji L1/L2/L3 ex post bez HALT-B.

---

## §0 — Background i origin

Niniejszy cykl jest **Poziom β** roadmapy zdefiniowanej w concept paper [[../../meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md]] (Poziom α, LOCKED 2026-05-21).

**Cel:** najprostszy techniczny test strukturalnej hipotezy CE-H: *czy stabilne equilibrium w TGP wymaga cosmic <Phi>_bg background?*

**Strategia testu (dichotomia):**

- **Phase 1a — N=2 isolation:** dwa solitony w pustej przestrzeni (Phi → 0 at infinity). **Pre-registered prediction:** brak stabilnego L*. *Pozytywny null result wzmacnia CE-H (stabilność wymaga bg).*
- **Phase 1b — N=2 + bg:** dwa solitony + non-zero <Phi>_bg representing rest of universe (Phi → Phi_bg at infinity). **Pre-registered prediction:** stabilne L* istnieje.

**Dichotomia kluczowa:** *both* muszą wyjść zgodnie z prediction by CE-H przetrwało Poziom β. Jeśli (1a) ma equilibrium → CE-H jest niepotrzebne (background nie jest required). Jeśli (1b) NIE ma equilibrium → CE-H niepoprawne (background nie wystarczy).

---

## §1 — Pre-registered falsifiers (LOCKED 2026-05-21)

### F-β-1 — Isolation null result (NULL CONSISTENCY)
**Pre-registracja:** Phase 1a (N=2 isolation) MUSI dać brak stabilnego L* (lub L* = 0 collapse / L* = ∞ no binding).
**Tolerancja:** zero. Jeśli stable L* istnieje w isolation, CE-H falsyfikowane (background nie był required).
**Severity:** STRUCTURAL — failure → CE-H wrong.

### F-β-2 — Background-stabilized equilibrium (POSITIVE PREDICTION)
**Pre-registracja:** Phase 1b (N=2 + <Phi>_bg) MUSI dać stable L* > 0 finite, z d²E/dL² > 0.
**Tolerancja:** L* musi być w fizycznie sensownym zakresie (1 < L*/(Phi_0)^(-1) < 10⁶, tj. od microscopic do mesoscopic w jednostkach Phi-substrate scale).
**Severity:** STRUCTURAL — failure → CE-H wrong.

### F-β-3 — Monotonic L*(<Phi>_bg) dependence (POSITIVE PREDICTION)
**Pre-registracja:** L* MUSI być monotonic function <Phi>_bg (gęstsze bg → krótsze L* lub luźniejsze bg → dłuższe L*, jeden trend, nie oscillation).
**Tolerancja:** sign konsystentny across parameter range.
**Severity:** SECONDARY — failure suggests structural inconsistency.

### F-β-4 — No fine-tuning required (POSITIVE METHODOLOGY)
**Pre-registracja:** stable L* musi istnieć dla **rozsądnego zakresu** parametrów (Phi_0, λ, <Phi>_bg), nie tylko dla wąskiego fine-tuned punktu.
**Tolerancja:** zakres ≥ factor 10 w przynajmniej jednym z parameters by NOT-fine-tuned.
**Severity:** STRUCTURAL — fine-tuned solution = potential Lakatos red flag.

### F-β-5 — Self-consistency closure (POSITIVE STRUCTURAL)
**Pre-registracja:** rozwiązanie (EQ-2) z dwóch source particles MUSI być konsystentne z (EQ-1) bound state energy (no infinite regress, fixed point exists).
**Tolerancja:** convergence numeryczne lub analityczne demonstration.
**Severity:** STRUCTURAL — failure → fundamental math problem z (EQ-1)-(EQ-6) systemem.

---

## §2 — 10 forbidden post-hoc moves (inherited z FFS cycle)

Jakikolwiek z poniższych w trakcie Poziom β = **automatyczny HALT-B**:

1. Modyfikacja F-β-1...F-β-5 tolerancji ex post.
2. Renaming falsified prediction by avoid penalty.
3. Dodawanie "additional fields" by retain failed prediction.
4. Cherry-picking parameter range gdzie prediction zgodna.
5. Re-defining "equilibrium" by include unstable points.
6. Re-defining "stable" by include marginal stability.
7. Switching ansatz mid-cycle by avoid failure.
8. Hardcoding FP T_pass=True (strict cycle 1/2/7 violation).
9. Using DEC budget powyżej 1 (max budget exceeded).
10. Introducing new axioms by rescue failed prediction (R3 threshold violation).

---

## §3 — L1 / L2 / L3 falsification map

### L1 (native, cycle-local):
- F-β-1: isolation null result
- F-β-2: bg-stabilized equilibrium
- F-β-3: monotonic dependence
- F-β-4: no fine-tuning
- F-β-5: self-consistency closure

### L2 (framework targets, post-cycle):
- TGP-native (EQ-1)-(EQ-6) self-consistency in toy regime
- Bridge do FFS quark mass spectrum (Phi_0_local jako relational quantity)

### L3 (cross-cycle propagation, conditional on success):
- R3 trigger linia 3 evidence ✓ (jeśli oba 1a/1b zgodne z prediction)
- FFS Phase 4 C6 PARTIAL → potential RESOLVED_STRUCTURALLY
- Concept paper Poziom α → upgraded od "structural conjecture" do "structural feature verified at toy level"
- Poziom γ (N-body + continuum cosmologiczny) authorization gate **conditional** na sukcesie Poziom β

---

## §4 — Phase plan (5 faz)

### Phase 0 — Balance sheet
**Scope:** external inputs, LOCKED structural axioms, derived outputs claim, tautology + falsifiability checks, anti-BD-drift check.
**Deliverable:** `Phase0_balance.md`
**Estimated:** 1 dzień (today).

### Phase 1a — N=2 isolation (null test)
**Scope:** Two solitons in pure Phi-substrate (Phi → 0 at infinity). Compute E(L), find stationary points, check stability.
**Deliverable:** `Phase1a_sympy.py` + `Phase1a_results.md`
**Pre-registered prediction:** brak stable L* (NULL CONSISTENCY for F-β-1).
**Estimated:** 1-2 dni.

### Phase 1b — N=2 + background (positive test)
**Scope:** Two solitons + non-zero <Phi>_bg (Phi → Phi_bg at infinity, Phi_bg = Phi_0 - δ dla małego δ > 0). Compute E(L), find L*, check stability.
**Deliverable:** `Phase1b_sympy.py` + `Phase1b_results.md`
**Pre-registered prediction:** stable L* > 0 finite (POSITIVE for F-β-2).
**Estimated:** 1-2 dni.

### Phase 2 — Parameter scan + monotonicity
**Scope:** L*(<Phi>_bg, Phi_0, λ) scan. Check F-β-3 (monotonicity) i F-β-4 (no fine-tuning).
**Deliverable:** `Phase2_sympy.py` + `Phase2_results.md`
**Estimated:** 1 dzień.

### Phase 3 — Self-consistency closure check
**Scope:** F-β-5 verification. Verify że (EQ-2) z 2 source particles + own contributions converges to consistent (EQ-1) solution.
**Deliverable:** `Phase3_sympy.py` + `Phase3_results.md`
**Estimated:** 1 dzień.

### Phase FINAL — Closure
**Scope:** Verdict per falsifier, claim_status assignment, R3 trigger status update, cross-cycle propagation map (DEFERRED actual updates until R2 audit).
**Deliverable:** `Phase_FINAL_close.md`
**Estimated:** 0.5 dnia.

**Total estimated:** 5-7 dni.

---

## §5 — Discipline declarations (binding)

### Strict cycle 1/2/7 pattern
- 0 hardcoded FP T_pass=True
- LIT/INVENTORY tests informational only (status_class: INFORMATIONAL)
- Substantive FP tests must compute from inputs, then compare to pre-registered threshold

### Max DEC budget = 1
- Across całego cyklu, max 1 "decision" test (deterministic choice between alternatives)
- Current usage: 0/1

### R1+R2+R3 discipline
- R1 (research-tier): all 5 phases permissive flagging
- R2 (integration audit): scope = items aggregated z R1, deferred do post-cycle audit
- R3 (multi-line convergence): existing 3 lines z concept paper (Phase 4 FFS + Archimedean + structural argument). Poziom β = **verification** linii 3, nie new evidence (jeśli passes, R3 confirmed)

### Anti-BD-drift LOCK
- Phase 0 anti-BD-drift check explicit
- NO fitting do QCD, ΛCDM, lub innych frameworks
- Native equations FIRST; mapping post-hoc bonus only

### Anti-Lakatos LOCK
- Pre-registration LOCKED before any sympy
- Each phase reports honestly vs pre-registration
- Any tolerance modification ex post = HALT-B

---

## §6 — Risk register

| ID | Risk | Mitigation | Severity |
|----|------|-----------|----------|
| R1 | Choice of soliton ansatz (kink/vortex/hedgehog) | Phase 0 explicit declaration; pick simplest = 1D kink in Z2 toy | MEDIUM |
| R2 | Kink-antikink anihilują w isolation | Pre-registered prediction NULL for isolation; this is FEATURE not bug | LOW (handled) |
| R3 | Sympy może nie convergeować dla L>>1 lub L<<1 | Numeric fallback; document limit ranges | MEDIUM |
| R4 | "Equilibrium" definition ambiguity | Phase 0 explicit: dE/dL=0 AND d²E/dL²>0 | LOW |
| R5 | Calculational complexity dla 2D vortex | Phase 0 default to 1D kink; 2D/3D deferred to Poziom γ | LOW |
| R6 | Background <Phi>_bg parametrization | Phase 0 explicit: Phi_bg = Phi_0 · (1 - δ), δ → 0 limit checked | MEDIUM |
| R7 | Numerical vs analytical sympy results | Prefer analytical; fallback to numerical with seed reproducibility | LOW |
| R8 | Cross-cycle bridge nadinterpretacja sukcesu Poziom β | Poziom β = TOY test only; full claim wymaga Poziom γ; explicit declaration | HIGH |
| R9 | Scope creep do cosmologicznych predictions | LOCKED scope: tylko 2-particle toy; cosmologia (H_0, etc.) DEFERRED do Poziom γ | HIGH |
| R10 | R3 trigger over-claim | Poziom β potwierdza linię 3, NIE genertuje nowych linii; explicit | MEDIUM |

---

## §7 — Authorization gate dla Phase 1a

Po Phase 0 (balance sheet) zakończonej, **WYMAGANA jest osobna autoryzacja user** dla Phase 1a (pierwsza sympy implementacja).

**Sekwencja autoryzacji:**
- ✅ 2026-05-21: Poziom β scaffold (READMEs + Phase 0) authorized by "działaj"
- ⏳ Phase 1a sympy → wymaga osobnego potwierdzenia
- ⏳ Phase 1b sympy → wymaga osobnego potwierdzenia
- ⏳ Phase 2 sympy → wymaga osobnego potwierdzenia (lub batch auth)
- ⏳ Phase 3 sympy → wymaga osobnego potwierdzenia (lub batch auth)
- ⏳ Phase FINAL → wymaga osobnego potwierdzenia

**Alternative:** user może udzielić batch authorization dla all phases at once, ale każda phase nadal raportuje vs pre-registered falsifiers.

---

**END OF README — Poziom β cycle BINDING contract LOCKED 2026-05-21**
