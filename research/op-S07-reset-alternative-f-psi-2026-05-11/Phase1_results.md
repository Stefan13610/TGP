---
title: "Phase 1 results — S07 alternative f(ψ) family enumeration + GR-limit + β_ppE predictions"
date: 2026-05-13
parent: "[[./README.md]]"
phase: 1
sympy_total: "12/12 PASS (100%)"
substance_metrics: "10 FP (83.3%) / 2 LIT (16.7%) / 2 DEC (separate); 100% non-trivial"
phase1_outcome: "S07 alternative f(ψ) families enumerated; GR-limit verified; β_ppE^poly(α) = (15/16)·α linear scaling derived; GWTC-3 compatible range α ∈ [-0.832, 0.832]; M9.1'' α=-4 excluded; recovery region NON-TRIVIAL z c_0·κ_σ=4/3 LOCK"
---

# Phase 1 results — S07 alternative f(ψ)

## §0 — Headline

```
████████████████████████████████████████████████████████████████████
█  op-S07-reset-alternative-f-psi-2026-05-11                       █
█  Phase 1: 12/12 PASS — 10 FP / 2 LIT / 0 hardcoded               █
█                                                                  █
█  KEY FINDING:                                                    █
█    β_ppE^poly(α) = (15/16)·α  [linear scaling]                   █
█    α_M911 = -4  →  β = -15/4 (FALSIFIED 5σ)                      █
█    GWTC-3 1σ:  α ∈ [-0.832, 0.832]                               █
█                                                                  █
█  ZERO-β REGION NON-TRIVIAL z c_0·κ_σ = 4/3 LOCK                  █
█                                                                  █
████████████████████████████████████████████████████████████████████
```

## §1 — Test results

| Test | Klasa | Status | Substance |
|---|---|---|---|
| T1 | FP | PASS | Baseline f=1 trivial GR; β_ppE=0 |
| T2 | FP | PASS | Polynomial f(ψ)=1+α(ψ-ψ_0); GR-limit built in |
| T3a | FP | PASS | M9.1'' Taylor at ψ_0=1: df/dψ=-4 (α_M911=-4) |
| T3b | FP | PASS | **β_ppE^poly(α) = (15/16)·α linear scaling** |
| T4 | FP | PASS | GWTC-3 compat range: α ∈ [-0.832, 0.832]; M9.1'' excluded |
| T5a | FP | PASS | Quadratic family GR-limit preserved |
| T5b | FP | PASS | β_ppE^quad at 2.5PN leading = polynomial (higher PN suppressed) |
| T6a | FP | PASS | Transcendental f=exp(α(ψ-ψ_0)); GR-limit exp(0)=1 |
| T6b | FP | PASS | Transcendental β_ppE leading = polynomial (Taylor coef α) |
| T7 | FP | PASS | Cross-cycle c_0·κ_σ=4/3 LOCK (emergent-metric Phase 4) |
| T8 | LIT | PASS | GWTC-3 |β_ppE| ≤ 0.78 + M9.1'' rejection 5.02σ |
| T9 | LIT | PASS | LIGO-O5 A+ ~2027 SNR=15.05σ threshold (PR-002) |

**Declarative (separate):**
- T10: Anti-Lakatos commitment
- T11: S05 preservation across f(ψ) alternatives

## §2 — Six P-requirements (Phase 1 contribution)

| P | Resolution z evidence |
|---|---|
| P1 | f(ψ) family enumeration: trivial / polynomial / quadratic / transcendental (T1+T2+T5a+T6a FP) |
| P2 | β_ppE^poly(α) = (15/16)·α explicit z Δe_2_native chain (T3b FP) |
| P3 | GR-limit recovery dla wszystkich 3 alternative families (T1+T2+T5a+T6a FP) |
| P4 | GWTC-3 compat range α ∈ [-0.832, 0.832] (T4 FP); Phase 2 Bayesian fit deferred |
| P5 | Cross-cycle c_0·κ_σ=4/3 LOCK preserved (T7 FP) |
| P6 | S05 preserved across wszystkie alternatives (T11 DEC) |

**P1-P3, P5, P6 Phase 1 RESOLVED. P4 PARTIAL (analytical range; Bayesian fit Phase 2).**

## §3 — Substantive findings

### §3.1 — Linear scaling β_ppE^poly(α) = (15/16)·α

**Najistotniejszy wynik Phase 1**: dla polynomial f(ψ) family, β_ppE^(b=-1) skaluje się
LINIOWO z leading-order Taylor coefficient α:

```
β_ppE^poly(α) = (15/16) · α
```

**Konsekwencje:**

- **α = 0** (trivial GR): β_ppE = 0 — passes GWTC-3 by construction (no signal)
- **α = -4** (M9.1''): β_ppE = -15/4 — FALSIFIED 5.02σ (Phase 2 RERUN 2026-05-09)
- **|α| ≤ 0.832** (GWTC-3 1σ): recovery region

### §3.2 — GWTC-3 compatibility range

Bound |β_ppE| ≤ 0.78 (Abbott+2021):

```
|(15/16)·α| ≤ 0.78
|α| ≤ 0.78 · (16/15) = 0.832
```

**Recovery region:** α ∈ [-0.832, 0.832] — non-empty, signifikantna szerokość rzędu jednostki.

### §3.3 — Cross-cycle consistency z emergent-metric Phase 4

Wynik Phase 4 emergent-metric: zero-β region w {A,B,C} family is non-trivial.
Kombinacja c_0·κ_σ = 4/3 EXACT (joint product LOCK from independent cycles
`op-c0-derivation-from-substrate` + `op-kappa-sigma-2body-PN`) provides σ-coupling
compensation term in Δe_2_native:

```
Δe_2_native = -4·ξ_3 + 4 - a_3/8 + c_0·κ_σ
            = -4·ξ_3 + 4 - a_3/8 + 4/3
```

Zero-β region w {A,B,C} family corresponds to combinations w {α, ξ_3, a_3} space gdzie ten
wyrażenie kasuje się do zera. **NIE wymaga α = 0 trivially**; może być finite α
z compensating ξ_3, a_3 contributions.

### §3.4 — Family universality

3 enumerated families (polynomial / quadratic / transcendental) wszystkie give **SAME**
β_ppE^(2.5PN leading order) = (15/16)·α gdzie α = df/dψ|_{ψ_0}. Different families differ
TYLKO at higher PN orders (3PN+).

**Implikacja:** at 2.5PN GWTC-3 precision, **tylko one parameter α matters** for recovery
test. Family-specific differences (β_quad, exp series) appear przy BH ringdown / photon ring
channels (BH5 / ε.1 — separate cycles).

## §4 — Anti-Lakatos commitment status

Per T10 DEC + PR-010 LOCKED:

- **Brak H1c/H1d backstop**
- Pre-bounded recovery_scope: α ∈ [-0.832, 0.832] z mandatory GR-limit f(ψ_0)=1
- Jeśli LIGO-O5 A+ ~2027 single-event excludes all α in this range → **H1b verdict**:
  framework architecture revision OR M9.1'' framework-level falsification accepted

## §5 — Status & next steps

**Phase 1 GATE: ✅ OPEN — 12/12 PASS z substantive linear-scaling derivation.**

**Phase 2 plan (next session):**
- Bayesian GWTC-3 fit dla pre-bounded α range — najlikely α posterior z 90 BBH combined
- Higher-PN corrections from quadratic + transcendental families (BH5 + ε.1 channel
  predictions)
- Phase FINAL closure z dedicated H1a/H1b verdict + falsification map check

**Estymata pozostałych sesji:** 2-4 sesje (z 5-8 oryginalnej multi-session estymata
**SKRÓCONE** dzięki linear-scaling discovery upraszczającemu fit).

---

**Phase 1 close.** Linear scaling β_ppE^poly(α) = (15/16)·α + recovery region α ∈ [-0.832,
0.832] derived first-principles. Phase 2 Bayesian fit gotowy do uruchomienia.
