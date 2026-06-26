---
title: "Phase 2 — CE-H items audit (op-R2-integration-audit-CE-H-FFS-2026-05-22)"
type: phase_audit
status: COMPLETE
date_completed: 2026-05-22
phase: 2
parent_cycle: op-R2-integration-audit-CE-H-FFS-2026-05-22
items_audited: 4 (CE-H-1, CE-H-2, CE-H-3, CE-H-4)
---

# Phase 2 — CE-H items audit

**Status:** COMPLETE 2026-05-22.
**Purpose:** Per-item structural audit dla 4 CE-H-source items (R2 scope §1.2 README).
**Audit methodology:** Path A (structural assessment) + Path B (alternative path search) + Path C (cross-cycle consistency) + Path D (pre-registration check).

---

## §0 — Verdict summary

| Item | Verdict | Justification (short) |
|------|---------|------------------------|
| CE-H-1 | ✅ **CLOSED** | Toy-model limitation explicit (D/L^α exogenous w 1D Z2); Poziom γ-1 scope LOCKED (F-γ-1 native 3D interaction test) |
| CE-H-2 | ✅ **CLOSED** | 1D Z2 substrate gives EXPONENTIAL (no native α); 3D U(1) Coulomb α=1 expected from standard propagator analysis; Poziom γ scope LOCKED |
| CE-H-3 | ✅ **CLOSED** | Dimensional analysis self-consistent within toy framework: [m]=1/L, [A]=E, [D]=E·L^α |
| CE-H-4 | ✅ **CLOSED** | Structural feature documented (boundary D_critical(α)); F-γ-4 pre-registered (LOCKED PENDING_POZIOM_GAMMA); NO over-claim QCD analog |

**Aggregate Phase 2:** 4 CLOSED + 0 DEFERRED + 0 ESCALATED (4/4 items).

---

## §1 — CE-H-1: D/L^α exogenous nature

### §1.1 Status assessment (Path A)

**Source:** CE-H Poziom β Phase 3 honest finding 2026-05-21 — native 1D Z2 substrate gives EXPONENTIAL exp(-m·√2·L), NOT power-law D/L^α.

**Structural facts (from Phase 3 results):**
- Variational two-soliton ansatz: $\Phi_{AK,K}(x; L) = -v \tanh(m(x+L/2)/\sqrt{2}) + v \tanh(m(x-L/2)/\sqrt{2}) + v$
- Numerical integration shows V_int → 0 as L → ∞ exponentially
- Fit results: V_int ~ -10.04 · exp(-1.40·L), R² = 0.9999 (exponential perfect fit)
- Power-law fit: R² = 0.967 (significantly worse)
- Native 1D Z2 interaction IS exponential, NIE power-law

**Phase 1b/2 D/L^α model status:**
- D/L^α was used dla bg-source representing "rest of universe" effect
- Demonstrated MECHANISM (bg can stabilize), NIE derivation z 1D Z2 substrate
- Phase 1b/2 success showed dichotomy: isolation (no equilibrium) vs bg (equilibrium exists)
- Mechanism verified; specific form of bg = exogenous modeling tool

### §1.2 Derivability check (Path A)

**Question:** Can D/L^α power-law be derived z 1D Z2 substrate?

**Answer:** NO. Native 1D Z2 substrate yields exponential interactions universally (confirmed Phase 3 to 1% accuracy: fitted decay m·√2 vs analytical m·√2).

**Reason:** In 1D, kink tail decay is exponential exp(-m·√2·x) because m is the field mass scale (V''(v) = 2λv²·tanh²... etc.). Power-law would require massless or long-range mode, which doesn't exist w 1D Z2.

### §1.3 Alternative path search (Path B)

**Alt-1: Multi-field 1D toy z effective long-range mode**
- Add Goldstone boson (continuous symmetry broken)
- 1D Z2 has discrete Z₂ symmetry — no continuous symmetry breaking → no Goldstone
- → Alt-1 violates 1D Z2 toy scope; would change AX-Z2

**Alt-2: Higher-dimensional native long-range emerges from 1D Z2 via integrate-out**
- Could 1D Z2 effective theory in higher-D give power-law?
- Standard QFT: integrating out heavy modes gives effective Lagrangian z higher-derivative terms
- These give modified short-distance behavior, NIE long-range power-law
- → Alt-2 doesn't work

**Verdict on alternatives:** Native 1D Z2 substrate cannot produce power-law. D/L^α IS exogenous w 1D Z2 — this is **structural feature**, NIE pre-registration error.

### §1.4 Cross-cycle consistency (Path C)

**Path forward LOCKED:**
- Concept paper Poziom α §10.6 open question Q1: "Jak konkretnie wygląda G(x-x') Greenowska funkcja Phi-substrate?" — answered for 1D (exponential), pending 3D
- Phase FINAL §6.2 F-γ-1 LOCKED: "native 3D U(1) interaction tail MUST be power-law lub logarithmic OR Coulomb-like"
- This is the structural argument: 3D U(1)+RP² PROPER framework should give natural long-range, complementing 1D Z2 exponential

**Cross-cycle integrity:** Phase 3 honest declaration + F-γ-1 pre-registration form coherent path: 1D Z2 limitation acknowledged, 3D U(1) future test specified.

### §1.5 Pre-registration check (Path D)

**Pre-registration:**
- Phase 0 EXT-7 explicit: "1D Z2 toy is methodological simplification; pełna 3D U(1)+RP² odłożona do Poziom γ"
- Phase 0 §2 AX-U1, AX-RP2: "NIE używana w 1D toy (EXT-7), zachowana dla Poziom γ"
- Phase FINAL §3.2 honest declaration: D/L^α was modeling tool, native derivation = Poziom γ scope
- → Pre-registration intact; toy limitation declared from Phase 0 onward

### §1.6 Verdict CE-H-1: CLOSED

**Justification:** Phase 3 honestly revealed native 1D Z2 substrate gives EXPONENTIAL exp(-m·√2·L), NIE power-law D/L^α. Phase 1b/2 D/L^α was **modeling tool** demonstrating mechanism (bg can stabilize), NIE derivation. Toy-model limitation **explicit from Phase 0**. Path forward LOCKED via F-γ-1 (3D U(1) native long-range test, Poziom γ-1 scope).

**Decision matrix mapping:** "Toy-model limitation explicit + Poziom γ-1 scope locked" → **CLOSED** per §3 README.

**R2 audit verdict:** CE-H-1 CLOSED. D/L^α exogeneity status structurally legitimate as toy limitation; Poziom γ-1 = next research direction for native 3D derivation.

---

## §2 — CE-H-2: α derivation gap

### §2.1 Status assessment (Path A)

**Source:** CE-H Phase 2 used exogenous α ∈ {0.5, 1, 2, 3} for parameter scan.

**Structural facts:**
- Phase 2 demonstrated CE-H mechanism robust across α range (factor 6 range)
- 4 distinct repulsion forms tested
- Equilibrium L*(α, D) exists in 20/20 grid cases
- α itself was scanned exogenously

### §2.2 Derivability check (Path A)

**Question:** Czy α derivable z TGP foundations?

**1D Z2 toy:**
- Native interaction is EXPONENTIAL (Phase 3 result)
- No native α scale in 1D Z2 substrate (α is power-law exponent which doesn't exist natively)
- → In 1D Z2, α is fundamentally exogenous

**3D U(1) expectation:**
- Standard scalar field in 3D: massless free propagator G(r) ~ 1/r → V_int ~ 1/r (Coulomb-like)
- Massive scalar: G(r) ~ e^(-mr)/r (Yukawa) → V_int ~ e^(-mr)/r
- 3D U(1) vortex-vortex: literature shows logarithmic interaction E ~ ln(L) for 2D vortices, modifications for 3D
- → 3D analog suggests α=1 (Coulomb) for massless modes, or logarithmic for vortex-vortex in 2D

**Pattern:** α derivation is intrinsically tied to dimensionality + field content:
- 1D massive Z2 → exponential (no power-law)
- 3D U(1) (massless) → α=1 (Coulomb)
- 2D U(1) vortex → logarithmic
- Higher-D modifications possible

**Conclusion:** α is **derivable** in 3D U(1) framework (Coulomb α=1 from standard propagator); **NOT derivable** in 1D Z2 (no power-law structure). This is dimensional dependence, NIE methodology gap.

### §2.3 Alternative path search (Path B)

**Alt-1: Effective α from RG flow in 1D Z2**
- RG analysis in 1D Z2 with critical regime
- 1D Z2 has phase transition at finite temperature; renormalization could induce power-law correlations
- But we're at T=0 ground state; RG flow doesn't add long-range mode
- → Alt-1 not applicable to T=0 1D Z2 toy

**Alt-2: α emerges from Φ-substrate self-consistency at finite source density**
- Possible: high source density (many particles) could induce effective long-range correlations
- But this requires N-body analysis (N≥3), not 2-body Phase 1b/2 scope
- Could be relevant for Poziom γ-2 (self-consistency closure z native bg)
- → Alt-2 = Poziom γ-2 scope

**Verdict on alternatives:** α in 1D Z2 fundamentally exogenous. In 3D U(1), α=1 expected from Coulomb-like propagator. Multi-particle self-consistency analysis = Poziom γ-2 scope.

### §2.4 Cross-cycle consistency (Path C)

**Inheritance:**
- Concept paper §4.2 EQ-2: $\langle\Phi\rangle(\vec{x}, t) = \langle\Phi\rangle_{frontier} + \sum_j G(\vec{x} - \vec{x}_j) \rho_j$
- G is Greenowska function — its form determines α
- 1D Z2: G ~ exp(-mr); 3D U(1) massless: G ~ 1/r
- → Concept paper explicit consistent z this verdict

### §2.5 Pre-registration check (Path D)

**Pre-registration:**
- Phase 2 README pre-registered scan range α ∈ {0.5, 1, 2, 3} explicit
- Phase 0 EXT-7 explicit 1D Z2 limitation
- F-γ-1 LOCKED: "3D U(1) native long-range power-law/logarithmic/Coulomb"
- → Pre-registration intact; α derivation deferred to Poziom γ

### §2.6 Verdict CE-H-2: CLOSED

**Justification:** α derivation gap is **dimensional artifact** of 1D Z2 toy (no power-law structure native), NIE methodology failure. 3D U(1) framework expected to provide α=1 (Coulomb) z standard massless scalar propagator analysis. Multi-particle self-consistency = Poziom γ-2 scope. Honest declaration consistent z Phase 0 1D limitation.

**Decision matrix mapping:** "α structural origin documented (1D Z2 limitation + 3D U(1) expected)" → **CLOSED** per §3 README.

**R2 audit verdict:** CE-H-2 CLOSED. α exogeneity in 1D Z2 is fundamentally dimensional; 3D Coulomb α=1 expected from standard propagator; explicit Poziom γ-1/γ-2 scope.

---

## §3 — CE-H-3: Dimensional structure verification

### §3.1 Status assessment (Path A)

**Source:** Phase 2 dimensional analysis (m, D, L) 2026-05-21.

**Structural facts (from Phase 2 results):**

Total energy: $E_{total}(L) = 2 E_K - A \cdot e^{-mL} + D/L^{\alpha}$

Dimensional analysis required:
- [L] = length
- [m] = 1/length (inverse length, kink mass scale)
- [mL] = dimensionless ✓ (argument of exp)
- [E] = energy
- [E_K] = energy ✓
- [A · e^{-mL}] = [A] · [dimensionless] = [A] = energy → [A] = E
- [D/L^α] = [D] · [length^(-α)] = E → [D] = E · L^α

For specific α:
- α=1: [D] = E·L (analog momentum, or coupling × length)
- α=2: [D] = E·L² (analog moment of inertia)
- General α: [D] depends on α

**Within toy framework, dimensional structure self-consistent.**

### §3.2 Derivability check (Path A)

**Question:** Czy dimensional structure jest internally consistent z TGP-native units?

**Within toy:**
- v (VEV) has dimension of field; in natural units (ħ=c=1), [v] = mass = 1/length
- m² = 2λv² → m = √(2λ)·v; if λ dimensionless, [m] = [v] = 1/length ✓
- [E_K] = (energy density × length) per unit transverse area; in 1D: [E_K] = energy ✓
- For exponential interaction A·exp(-mL): [A] = energy ✓ (energy times dimensionless)
- For power-law D/L^α: [D] = energy × L^α

**Consistency check:**
- All terms in E_total have dimension energy ✓
- mL is dimensionless ✓
- A, D have appropriate dimensions w toy units ✓

### §3.3 Alternative path search (Path B)

**Alt: Could dimensional inconsistency be hiding via numerical normalization?**

- Phase 1b/2 used numerical m=v=1 normalization
- This sets natural unit scale; A, D become dimensionless numbers w these units
- BUT: dimensional structure remains intact in unnormalized form
- Sympy verification (Phase 2) explicit kept dimensions; numerical evaluation was z m=v=1

**Verdict on alternatives:** No dimensional inconsistency. Numerical normalization preserves dimensional structure.

### §3.4 Cross-cycle consistency (Path C)

**Pattern 2.5 consistency:**
- TGP-native scales: Φ_0_local (mass scale), λ (dimensionless coupling), m_kink = √(2λ)·v
- These all have proper dimensions
- E_K and A, D in toy follow consistent dimensional scheme

**No conflicts z other cycles** (Phase 2 dimensional analysis is internal).

### §3.5 Pre-registration check (Path D)

**Pre-registration:**
- Phase 2 sympy preserved dimensions in symbolic form
- Numerical evaluation z m=v=1 normalized choice
- Dimensional analysis IS part of Phase 2 documentation
- → No pre-registration gap

### §3.6 Verdict CE-H-3: CLOSED

**Justification:** Dimensional structure of Phase 1b/2 energy functional E_total(L) = 2E_K - A·exp(-mL) + D/L^α is **internally self-consistent** within toy framework: [m]=1/L, [A]=E, [D]=E·L^α. Numerical normalization (m=v=1) preserves dimensional structure. NO inconsistency wykryta.

**Decision matrix mapping:** "Dimensional analysis confirmed self-consistent" → **CLOSED** per §3 README.

**R2 audit verdict:** CE-H-3 CLOSED. Dimensional structure clean.

---

## §4 — CE-H-4: Confinement/deconfinement boundary structural feature

### §4.1 Status assessment (Path A)

**Source:** CE-H Phase 2 unexpected observation D_critical(α) 2026-05-21 (declared NOT pre-registered, R1 noteworthy).

**Structural facts (from Phase 2 §2.2):**

For E_total(L) = 2·E_K - A·exp(-mL) + D/L^α:

Critical D value:
$$D_{critical}(\alpha) = \frac{A \cdot (\alpha+1)^{\alpha+1} \cdot e^{-(\alpha+1)}}{\alpha \cdot m^\alpha}$$

- D < D_critical(α): stable bound state exists
- D > D_critical(α): no equilibrium (deconfined regime)
- Larger α → higher D_critical (more bg coupling tolerated)

This was noted Phase 2 as structural observation analog to QCD phase diagram (confinement/deconfinement transition at finite T).

### §4.2 Derivability check (Path A)

**Question:** Czy D_critical(α) jest structural feature TGP, czy 1D Z2 toy artifact?

**Derivation:**

D_critical comes from analytical solution of dE/dL = 0:
- Function g(u; α) = u^(α+1)·e^(-u) has interior maximum at u = α+1
- $g_{max}(\alpha) = (\alpha+1)^{\alpha+1} \cdot e^{-(\alpha+1)}$
- D_critical = A·g_max(α) / (α·m^α)

This is **mathematical consequence** of exponential attraction + power-law repulsion balance. As long as both modes coexist, similar boundary will exist (potentially different functional form in 3D).

**Status:** D_critical IS structural feature of mixed exponential+power-law potentials. NOT 1D-specific artifact (the math applies to any exponential vs power-law balance).

### §4.3 Alternative path search (Path B)

**Alt-1: D_critical exists only in 1D Z2 toy, NIE in 3D**

- 3D: native interaction is power-law (Coulomb 1/r) or logarithmic
- Need bg form for analog test
- If bg gives different functional form (e.g., also Coulomb), maybe no boundary
- → Alt-1: structural feature might NOT survive transition to 3D

**Counter-argument:**
- Generic feature: any "attractive at short + repulsive at long" potential has equilibrium for some parameter range
- Boundary will exist at parameters where minimum touches inflection (D_critical equivalent)
- → Boundary IS generic structural feature; specific functional form may change

**Verdict:** D_critical formula 1D-specific, but EXISTENCE of confinement/deconfinement-like boundary is structural feature of mixed potentials.

### §4.4 Cross-cycle consistency (Path C)

**F-γ-4 pre-registration:**
- Phase FINAL §6.2 LOCKED: "F-γ-4 — Confinement/deconfinement match observed: jeśli D_critical analog observed QCD T_c, ratio musi być w factor 10"
- Severity: SECONDARY (consistency check)
- Note: "speculative — może być niepotwierdzalna w Poziom γ scope"

**Status:** F-γ-4 explicitly pre-registered as **speculative test**. R2 audit confirms this is correct positioning — NOT over-claim QCD analog.

### §4.5 Pre-registration check (Path D)

**Pre-registration:**
- Phase FINAL §2.2 explicit declaration: "Nie pre-registered w Poziom β, więc nie claim tutaj. Noteworthy dla Poziom γ extension"
- F-γ-4 LOCKED with explicit "speculative" annotation
- → Pre-registration honest; structural feature declared, NIE adopted as Poziom β result

### §4.6 Verdict CE-H-4: CLOSED

**Justification:** D_critical(α) boundary IS structural feature of mixed exponential+power-law potential (math applies generically). Specific functional form is 1D-specific, but boundary existence generalizes. Phase FINAL §2.2 honest declaration "noteworthy, NOT claim" + F-γ-4 LOCKED pre-registration (speculative) = honest scientific reporting. **NO over-claim QCD analog**. Status: structural feature candidate, awaiting Poziom γ-3 cosmological test (F-γ-4).

**Decision matrix mapping:** "Structural status documented + Poziom γ-3 F-γ-4 scope locked" → **CLOSED** per §3 README.

**R2 audit verdict:** CE-H-4 CLOSED. Confinement/deconfinement boundary jest legitimate structural feature candidate; F-γ-4 pre-registered (LOCKED PENDING_POZIOM_GAMMA) preserves honest scope.

---

## §5 — Phase 2 aggregate

### §5.1 Verdict count

| Verdict | Count |
|---------|-------|
| CLOSED | 4 (CE-H-1, CE-H-2, CE-H-3, CE-H-4) |
| DEFERRED | 0 |
| ESCALATED | 0 |

### §5.2 Implications dla CE-H Poziom β claim_status

**Pre-R2 status:** CE-H Poziom β A− conditional (16/17 substantive PASS + Warstwa 1 honest fail + Warstwa 2 D/L^α exogenous).

**Post-Phase 2 R2 audit (CE-H items only):**
- All 4 items CLOSED z structural justification
- CE-H-1: D/L^α exogenous nature → declared toy limitation, Poziom γ-1 scope
- CE-H-2: α derivation gap → 1D Z2 fundamental limitation, 3D Coulomb expected
- CE-H-3: dimensional structure → self-consistent
- CE-H-4: confinement/deconfinement boundary → structural feature candidate, F-γ-4 LOCKED

**Joint impact:** All 4 CE-H R2 audit items closed strengthens CE-H Poziom β structural proof-of-principle. **However**: Warstwa 1 honest fail (T_P3_2) NIE addressed by R2 (it's R1-1 item, Phase 3 audit scope).

**CE-H claim_status proposed update:**
- claim_status A− conditional **PRESERVED** until R1-1 closed AND Poziom γ-1 success
- A− → A upgrade trajectory: requires R1-1 CLOSED (Phase 3 audit) + Poziom γ-1 success
- Current R2 verdict consolidates "structural proof-of-principle z honest caveats"

### §5.3 Cross-cycle propagation implications (deferred do Phase 4)

| Doc target | Impact from Phase 2 |
|------------|---------------------|
| op-CE-H-two-particle-equilibrium-2026-05-21/Phase_FINAL_close.md §5 | Add R2 audit cycle reference + per-CE-H-item verdicts |
| meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md §13 | Add Poziom β R2 audit closure annotation |
| meta/TGP_W_Z_THEORETICAL_LIMIT.md §6.5 | Path η cosmology toy extension entry |
| meta/PRE_REGISTERED_FALSIFIERS.md | F-γ-1..4 formal entries pending Poziom γ-1 |

### §5.4 R3 acceptance derivation chain (per Phase 0 §8)

**Pre-registered chain (Phase 0 §8):**
1. AX-S05 single Phi field with U(1) phase
2. AX-Z2 discrete Z₂ symmetry → topological solutions
3. Kink/vortex/hedgehog stability requires localized energy density profile
4. Localized energy density requires asymmetric Phi configuration vs surrounding bulk
5. Asymmetric configuration requires bulk Phi ≠ 0
6. Bulk Phi ≠ 0 = effectively ⟨Phi⟩_bg > 0
7. Two-particle stability in bulk Phi ≠ 0 requires balance (Phase 1b verified)
8. → CE-H structural feature follows z S05+Z₂ ontologii

**Phase 2 audit verification:**
- Steps 1-2: foundational, LOCKED
- Steps 3-4: kink stability standard result; localized energy density structural
- **Step 5-6 critical check:** czy "asymmetric configuration requires bulk Phi ≠ 0"?
  - In 1D Z2 isolation (Phi → 0 at infinity): kink-antikink ATTRACT and ANIHILATE (no equilibrium)
  - Phase 1a verified this: F-β-1 NULL CONFIRMED, dE/dL > 0 everywhere
  - W bulk Phi ≠ 0 (Phase 1b): equilibrium L* exists
  - → Step 5-6 STRUCTURALLY VERIFIED z Phase 1a/1b dichotomy
- Step 7-8: follows from 5-6 + Phase 1b structural result

**Verdict:** R3 acceptance derivation chain VALID — CE-H structural feature genuinely derives z S05+Z₂+U(1)+RP² minimal axioms, NIE separately postulated.

**Anti-Lakatos discipline check:** No new axiom required. Minimal axioms PRESERVED.

---

## §6 — Anti-Lakatos discipline check

- ✅ Each item verdict reported per decision matrix §3 LOCKED
- ✅ No threshold modifications ex post
- ✅ No new items added mid-cycle (forbidden #1)
- ✅ CE-H-1 D/L^α exogeneity DOES NOT undermine Phase 1b/2 — mechanism vs functional form distinct
- ✅ CE-H-4 structural feature CANDIDATE, NIE adopted as result (F-γ-4 LOCKED preserves speculation)
- ✅ R3 derivation chain explicit verified (anti-implicit-axiom-drift)
- ✅ Cross-cycle inheritance preserved

**Self-audit:** Phase 2 audit clean, anti-Lakatos preserved.

---

**END OF PHASE 2 — CE-H items audit LOCKED 2026-05-22**

**Aggregate Phase 2:** 4 CLOSED + 0 DEFERRED + 0 ESCALATED. Ready dla Phase 3 (R1 flag audit).
