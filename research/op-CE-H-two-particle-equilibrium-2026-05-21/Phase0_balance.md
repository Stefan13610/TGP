---
title: "Phase 0 — Balance sheet (op-CE-H-two-particle-equilibrium)"
type: phase_balance
status: LOCKED
pre_registration_date: 2026-05-21
phase: 0
parent_cycle: op-CE-H-two-particle-equilibrium-2026-05-21
parent_concept_paper: meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md
---

# Phase 0 — Balance sheet

**Status:** LOCKED 2026-05-21.
**Purpose:** Explicit accounting wszystkich inputów (external + structural axioms) i outputów (derived claims) PRZED jakimkolwiek sympy. Anti-Lakatos discipline.

---

## §1 — External inputs (co przyjmujemy z zewnątrz)

| ID | Input | Source | Status |
|----|-------|--------|--------|
| EXT-1 | TGP Phi-substrate Lagrangian L = (1/2)\|∂Phi\|² - (λ/4)(\|Phi\|²-Phi_0²)² | meta/TGP_GENERATED_SPACE_COSMOLOGY_2026-05-21.md §3.2; FFS Phase 1/3 | DERIVED w previous cycles |
| EXT-2 | Mexican hat SSB form V(Phi) = (λ/4)(\|Phi\|²-v²)² | Pattern 2.5 §3.5.6; FFS Phase 3 T_P4_1 | DERIVED w FFS |
| EXT-3 | Single kink solution Phi(x) = v·tanh(x·m/√2) gdzie m² = 2λv² (w 1D Z2 limit) | Standard QFT result; Coleman 1985 | LITERATURE |
| EXT-4 | Single kink energy E_kink = (2√2/3)·m·v² | Standard QFT result | LITERATURE |
| EXT-5 | Kink-kink interaction energy V_int(L) ~ -e^(-mL) at large L | Manton-Sutcliffe 2004 "Topological Solitons" Ch. 5 | LITERATURE |
| EXT-6 | Boundary condition: Phi → Phi_∞ at infinity (Phi_∞ = 0 dla Phase 1a, Phi_∞ = Phi_bg ≠ 0 dla Phase 1b) | TGP-native CE-H structural ansatz | THIS CYCLE |
| EXT-7 | 1D toy reduction (Z2 only, no U(1) for simplicity) | Methodological simplification | THIS CYCLE DECISION |
| EXT-8 | Sympy 1.12+ symbolic math + numerical fallback | Tool dependency | INFRASTRUCTURE |

**Comment on EXT-7:** wybór 1D Z2 toy jest **metodologiczna simplifikacja** by uzyskać najprostszy test. U(1) phase + RP² topology (full TGP) dadzą bogatszą strukturę solitonów (vortices), ale 1D kink wystarczy dla **strukturalnego** testu CE-H (czy bg potrzebny dla stabilności). Pełna 3D U(1)+RP² odłożona do Poziom γ.

**Comment on EXT-3..EXT-5:** literature inputs są dla **single kink** w izolacji. Two-kink system: będziemy **derive** native, NIE assume Manton-Sutcliffe form (anti-BD-drift).

---

## §2 — LOCKED structural axioms (co bierzemy z TGP foundations)

| ID | Axiom | Source | Role w cyklu |
|----|-------|--------|-------------|
| AX-S05 | Single scalar field Phi with U(1) phase | Foundations §1 | Lagrangian structure |
| AX-Z2 | Discrete Z2 symmetry Phi → -Phi | Foundations §2 | Kink topological charge |
| AX-U1 | U(1) phase symmetry | Foundations §3 | Pełna U(1) NIE używana w 1D toy (EXT-7), zachowana dla Poziom γ |
| AX-RP2 | RP² topology dla full TGP | Foundations §4 | NIE używana w 1D toy (EXT-7), zachowana dla Poziom γ |
| AX-CE | CE-H structural ansatz: stabilność wymaga <Phi>_bg | concept paper Poziom α | **TO JEST TO CO TESTUJEMY** |

**Critical observation on AX-CE:** to NIE jest nowy axiom. Jest **strukturalną konsekwencją** AX-S05 + AX-Z2 (kink-kink dynamics w pure mexican hat). Phase 0 deklaruje to explicit, Phase 1a/1b weryfikuje.

---

## §3 — Derived outputs (co claimuje wyprodukować)

| ID | Output | Phase | Pre-registered prediction |
|----|--------|-------|--------------------------|
| OUT-1 | E(L) function for two-kink system in isolation | Phase 1a | Monotonic in L (attractive lub repulsive), no L* > 0 stationary point z d²E/dL² > 0 |
| OUT-2 | E(L) function for two-kink + <Phi>_bg | Phase 1b | Has stationary point L* > 0 z d²E/dL² > 0 |
| OUT-3 | L*(<Phi>_bg, Phi_0, λ) scan | Phase 2 | Monotonic w <Phi>_bg; smooth across factor 10 parameter range |
| OUT-4 | Self-consistency convergence check | Phase 3 | (EQ-2) z 2 source particles converges to (EQ-1) solution within iteration budget |

---

## §4 — Tautology test (anti-circular reasoning)

**Question:** Czy zakładamy to, co chcemy udowodnić?

**Pre-registered answer:**

- **Phase 1a (isolation null):** NIE zakładamy że isolation nie ma equilibrium. Wyprowadzamy E(L) z Lagrangianu i sprawdzamy. **Result is genuine prediction**, nie tautologia.
- **Phase 1b (bg-stabilized):** NIE zakładamy że <Phi>_bg stabilizuje. Wyprowadzamy E(L) z Lagrangianu + boundary condition i sprawdzamy. **Result is genuine prediction**, nie tautologia.

**Anti-tautology check:** w Phase 1a ansatz pomija <Phi>_bg explicit. W Phase 1b ansatz includes <Phi>_bg jako parametr. Różnica ansatzów jest **strukturalna**, nie ad-hoc. Każdy ansatz daje compute-then-compare result, nie assume-then-verify.

---

## §5 — Falsifiability test

**Question:** Czy każdy claim jest falsyfikowalny?

| Output | Falsifiable? | How |
|--------|--------------|-----|
| OUT-1 (isolation no L*) | YES | Find stable L* w isolation → F-β-1 fail → CE-H wrong |
| OUT-2 (bg has L*) | YES | No stable L* w bg → F-β-2 fail → CE-H wrong |
| OUT-3 (monotonic scan) | YES | Non-monotonic L*(<Phi>_bg) → F-β-3 fail |
| OUT-4 (self-consistency) | YES | Divergent iteration → F-β-5 fail |

All 4 outputs są falsyfikowalne z explicit threshold per pre-registered F-β-1...F-β-5.

---

## §6 — Anti-BD-drift check (CRITICAL)

**Question:** Czy w jakimkolwiek miejscu fitujemy do SM/QCD/ΛCDM frameworks?

**Pre-registered answer:** NIE.

- Phase 1a/1b: TGP Phi-substrate Lagrangian only. NO comparison z QCD potential, NO comparison z Higgs vacuum.
- Phase 2: parameter scan w jednostkach TGP-native (Phi_0, λ, <Phi>_bg). NO comparison z observed cosmological parameters (H_0, Omega_m). Cosmologia DEFERRED do Poziom γ.
- Phase 3: self-consistency check w TGP-native equations. NO fitting do empirycznych obserwacji.
- Phase FINAL: claim_status based on F-β-1...F-β-5 verdict only. NO claim "TGP reproduces QCD/ΛCDM".

**Methodology:** Native equations FIRST. Mapping do other frameworks DEFERRED do post-Poziom-β R2 audit lub Poziom γ implementations.

**User explicit (2026-05-21):**
> "tworzymy natywne równania TGP, sprawdzamy czy wyniki są zgodne z pomiarami, a opcjonalnie robimy mapowanie"

---

## §7 — Independent-path cross-validation

**Pre-registered cross-checks:**

- **Path A — Analytical:** explicit derivation E(L) z Lagrangianu, find dE/dL=0 analytically (where possible).
- **Path B — Numerical:** sympy nsolve dla dE/dL=0; check d²E/dL² sign numerically.
- **Path C — Limit checks:** L → 0 limit (collapse), L → ∞ limit (free particles); check energy makes sense.
- **Path D — Background limit:** <Phi>_bg → 0 limit dla Phase 1b should recover Phase 1a result (smooth limit).

If Path A vs B konfliktują → flag i investigate (R1 research-tier permissive).
If Path C/D limits fail → structural problem, escalate to R2.

---

## §8 — Choice of ansatz declaration

**Decision (locked in Phase 0, before sympy):**

**1D Z2 kink ansatz:**

Single kink:
$$\Phi_{kink}(x; x_0) = v \cdot \tanh\left(\frac{m(x - x_0)}{\sqrt{2}}\right), \quad m^2 = 2\lambda v^2$$

Single anti-kink:
$$\Phi_{antikink}(x; x_0) = -v \cdot \tanh\left(\frac{m(x - x_0)}{\sqrt{2}}\right)$$

**Two-soliton ansatz (Phase 1a — isolation):**

Trial: superposition + correction. Standard Manton method:
$$\Phi_{2K}(x; L) = \Phi_{kink}(x; -L/2) + \Phi_{antikink}(x; +L/2) - v_{\text{correction}}$$

gdzie v_correction zapewnia boundary condition Phi → 0 at infinity. **NIE assume Manton form a priori** — wyprowadzimy w Phase 1a sympy.

**Two-soliton ansatz (Phase 1b — bg):**

Modified boundary: Phi → Phi_bg at infinity, gdzie Phi_bg = v · (1 - δ), δ > 0 small parameter.

$$\Phi_{2K,bg}(x; L) = \Phi_{kink}(x; -L/2) + \Phi_{antikink}(x; +L/2) + \Phi_{bg}$$

z Phi_bg jako external boundary value (representing "rest of universe" contribution).

**Why kink+antikink (not kink+kink):** kink+kink ma topological charge ±2, anihiluje only przez quantum tunneling (long timescale). Kink+antikink ma total charge 0, classical anihilation channel. **Both physically relevant**, kink+antikink jest **prostszy test** bo finite-energy configuration trywialnie istnieje.

**Honest caveat:** kink+antikink w isolation **standardly** anihilują (collapse to vacuum), no equilibrium. To jest **wynikłe** prediction Phase 1a, nie assumption. Sympy weryfikuje.

---

## §9 — Tests planned dla Phase 1a/1b/2/3 (counts only)

| Phase | T_xxx tests planned | DEC budget | LIT/INVENTORY | Substantive FP |
|-------|---------------------|-----------|---------------|----------------|
| 1a | 5 (analytical E(L), stationary points, stability, limits) | 0 | 1 (Manton form check) | 4 |
| 1b | 5 (analytical E(L,bg), stationary points, stability, limits) | 0 | 0 | 5 |
| 2 | 4 (parameter scans, monotonicity, fine-tuning check) | ≤1 | 0 | 4 |
| 3 | 3 (self-consistency convergence check) | 0 | 0 | 3 |
| Total | 17 | ≤1 | 1 | 16 |

**Substantive FP ratio target:** ≥90% pass dla CE-H structural verification.
**Hardcoded T_pass=True count target:** 0 (strict cycle 1/2/7).

---

## §10 — Literature checkpoint (informational, NOT validation)

**Status:** INFORMATIONAL. Literature jest **kontekst**, nie target.

**Anchors:**
- Coleman 1985 "Aspects of Symmetry" — single kink classical solution
- Manton-Sutcliffe 2004 "Topological Solitons" — kink-kink interactions
- Rajaraman 1982 "Solitons and Instantons" — Mexican hat models
- Drazin-Johnson 1989 "Solitons" — 1D soliton equilibria

**Key observation:** literature **predicts** that pure isolation kink-antikink **lacks** stable L* > 0 equilibrium (kink-antikink attract, anihilate at L=0). To **jest** dokładnie nasza Phase 1a pre-registered prediction.

**This is feature, not bug:** literature confirms isolation null result, ALE literature NIE zawiera CE-H bg-stabilization mechanism. Phase 1b jest **novel** test.

**Anti-BD-drift discipline:** używamy literature dla **single kink** properties (EXT-3..EXT-5). Two-kink + bg system **wyprowadzamy native**, NIE adopting literature ansatz.

---

## §11 — Open questions (do rozwiązania w Phase 1+)

1. Dokładna forma ansatz dla Phi_2K,bg(x; L, δ) — analytical lub variational?
2. Czy E(L; δ) ma analytical closed form, czy wymaga numerics?
3. L*(δ) scaling — power law? exponential? other?
4. Limit δ → 0 — smooth transition do Phase 1a result?
5. Czy result jest sensitive do exact form of "boundary" Phi_bg (constant vs gradient)?

Wszystkie pre-registered jako open; addressed w odpowiednich fazach.

---

## §12 — Status końcowy Phase 0

- ✅ External inputs (EXT-1..EXT-8) inventoried
- ✅ LOCKED structural axioms (AX-S05, AX-Z2, AX-U1, AX-RP2, AX-CE) declared
- ✅ Derived outputs (OUT-1..OUT-4) z pre-registered predictions
- ✅ Tautology test passed (no circular reasoning)
- ✅ Falsifiability test passed (all outputs falsifiable)
- ✅ Anti-BD-drift check passed (no fitting do other frameworks)
- ✅ Independent-path cross-validation declared (Paths A/B/C/D)
- ✅ Ansatz declaration locked (1D Z2 kink + antikink)
- ✅ Test counts pre-registered (17 total, ≤1 DEC)
- ✅ Literature checkpoint informational (NOT validation target)
- ✅ Open questions identified

**Phase 0 LOCKED 2026-05-21. Ready for Phase 1a authorization.**

---

**END OF PHASE 0 — Balance sheet LOCKED 2026-05-21**
