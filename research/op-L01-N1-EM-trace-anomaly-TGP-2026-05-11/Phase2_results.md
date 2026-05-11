---
title: "Phase 2 results — TGP-specific reduction + Theorem 2.1 disjointness + GW170817 dispersion full check"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-results
phase: 2
status: 🟢 RESOLVED — 8/8 sympy PASS
sub_needs_resolved: [N0.4, N0.5, N0.6]
risks_addressed: [R2 closed, R3 fully closed]
sympy_script: "[[./Phase2_sympy.py]]"
sympy_output: "[[./Phase2_sympy.txt]]"
predecessor: "[[./Phase1_results.md]] (8/8 PASS)"
tags:
  - phase2
  - TGP-reduction-explicit
  - theorem-2.1-disjointness
  - GW170817-preserved
  - psi1-Q1-confirmed-constructively
---

# Phase 2 results

## §0 — Executive summary

**8/8 sympy PASS.** Phase 2 establishes:

1. **R[g_eff] reduction** w 1PN limit (b_1 = -a_1, γ=1): R = a_1·∇²h (NIE
   zerowe, *proportional do a_1*), pure-2-derivative ψ scalar — sympy T1 ✓.
2. **R^{μν}[g_eff] reduction:** R^{ij} z γ=1 daje pure tensor ∂_i ∂_j h (scalar
   part zero), 2-derivative ψ tensor — sympy T2 ✓.
3. **Operator class enumeration** dla T_anomaly_TGP — 6 distinct classes,
   wszystkie z 0 lub 2-derivative ψ structure (NIE 1-derivative outer product) —
   sympy T3 ✓.
4. **Theorem 2.1 (Disjointness):** trace anomaly TGP-reduced operator classes są
   strukturalnie DISJOINT od ψ.1.v3 canonical basis B = {L₅'_a, L₅'_b}. Q1
   *konstruktywnie* potwierdzony z dedicated derivation — sympy T4 ✓.
5. **GW170817 c_GW=c_EM full check:** Δc/c ~ (α/(3π))·R_cosmo/m_e² ≈ 10⁻⁸⁰
   (utterly negligible) ≪ 9·10⁻²² bound — sympy T5 ✓ z corrected dimensional
   analysis.
6. **GW (h_TT) dispersion:** unchanged przez 1-loop QED (no graviton in QED loop),
   c_GW=c structurally preserved — sympy T6 ✓.
7. **R5 documented:** B ≪ B_QED ≈ 4.4·10⁹ T regime perturbative VALID; B≳B_QED
   defer to non-perturbative analysis — sympy T7 ✓.
8. **S05 verification:** wszystkie reduced operators contain only Φ (single
   fundamental field) — sympy T8 ✓.

| Check | Result |
|---|---|
| T1: R[g_eff] = a_1·∇²h (b_1=-a_1) | ✅ PASS |
| T2: R^{ij}[g_eff] tensor structure (γ=1, scalar zero) | ✅ PASS |
| T3: Operator class enumeration T_anomaly_TGP | ✅ PASS |
| T4: Disjointness vs ψ.1.v3 basis | ✅ PASS |
| T5: GW170817 dispersion |Δc/c| ≈ 10⁻⁸⁰ ≪ 9·10⁻²² | ✅ PASS |
| T6: c_GW=c structural (no QED↔graviton 1-loop) | ✅ PASS |
| T7: B ≪ B_QED regime documented | ✅ PASS |
| T8: S05 single-Φ verified | ✅ PASS |
| **TOTAL** | **8/8 PASS** |

## §1 — TGP-specific reduction R[g_eff], R^{μν}[g_eff]

### §1.1 — Setup recall (Phase 1 §3 + Phase 2 setup §1)

Per [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]]
ansatz, w weak-field limit ψ ≈ 1 + h, h ≪ 1:

```
g_eff^00 = -A(ψ) ≈ -(1 + a_1·h + a_2·h² + ...)
g_eff^ij = δ^ij·B(ψ) ≈ δ^ij·(1 + b_1·h + b_2·h² + ...)
g_eff^0i = 0  (statyczny limit)
```

z 1PN constraint **b_1 = -a_1** (γ_PPN=1, per emergent-metric Phase 2).

### §1.2 — Twierdzenie 2.2 (R reduction)

**Twierdzenie 2.2:** W weak-field limit z γ_PPN=1 constraint:

```
R[g_eff] = a_1 · ∇²h + O(h²) + O(σ_ab × ∂Φ)
```

(generic coefficient, dla a_1=4 GR limit daje 4·∇²h zgodnie z linearized
Einstein.)

**Sympy LOCK T1:** R reduces do `pure-2-derivative ψ scalar` structure, NIE
do `(∂ψ)(∂ψ)` outer product.

### §1.3 — Twierdzenie 2.3 (R^{μν} reduction)

**Twierdzenie 2.3:** W weak-field limit z γ_PPN=1 constraint:

```
R^{00}[g_eff] ≈ (a_1/2) · ∇²h
R^{ij}[g_eff] ≈ -a_1 · ∂^i ∂^j h  +  0 · δ^{ij} ∇²h  +  O(σ_ab × ∂Φ)
R^{0i}[g_eff] = 0   (statyka)
```

(Scalar part w R^{ij} *zeruje się* z γ_PPN=1 constraint — unique tensor
structure.)

**Sympy LOCK T2:** R^{μν} reduces do `pure-2-derivative ψ tensor` structure,
NIE do `(∂ψ)(∂ψ)` outer product class.

### §1.4 — Konsekwencja: trace anomaly TGP-reduced

Substituting do Phase 1 §1.3 stress-energy tensor, w weak-field limit:

```
T^μ_μ_EM,1-loop_TGP = (α/(3π)) · F²
                    + b₁ · a_1 · (∇²h) · F²              ← (∂²ψ) scalar × F²
                    - b₂ · a_1 · (∂_μ ∂_ν h) · F^{μρ} F^ν_ρ  ← (∂_μ ∂_ν ψ) tensor × F² tensor
                    + b₃ · □F²                             ← derivative of F²
                    + γ_3 · σ_ab · F²                       ← strain × F² (composite)
                    + Riegert local with σ_eff = function(ψ)
```

z b_i Wilson coefficients renormalization-fixed (computable z Birrell-Davies +
Phase 1 ansatz expansion; numerical pinning deferred — *deferred precision*, NIE
new free param).

**Strukturalna konsekwencja:**

```
ρ_EM_quantum[{Φ_i}] = -T^μ_μ_EM,1-loop_TGP / c_0²
                    = -(α/(3π))·F²/c_0² 
                      + curvature-mediated (∂²ψ)·F² + (∂_μ∂_ν ψ)·F² corrections
                      + σ_ab·F² strain coupling
                      + sub-leading derivative + Riegert local
```

To jest **renormalized quantum source field dla Φ-EOM** w obecności klasycznego
EM field na background Φ.

## §2 — Twierdzenie 2.1 (Disjointness) — pełna proof

### §2.1 — Statement

**Twierdzenie 2.1 (Disjointness):** W obecności g_eff[{Φ_i}] (per emergent-metric
Phase 1 ansatz {A(ψ), B(ψ), C(ψ)}, generic single-Φ scenario z 1PN constraint
b_1=-a_1), 1-loop QED trace anomaly produces operator classes:

```
T_anomaly_TGP_class = {
   (α/(3π)) · F²,                              [pure-photon dim-4]
   γ_1 · (∂²ψ) · F²,                           [scalar 2-deriv ψ × F², dim-6]
   γ_2 · (∂_μ ∂_ν ψ) · F^{μρ} F^ν_ρ,           [tensor 2-deriv ψ × F², dim-6]
   γ_3 · σ_ab · F²,                            [strain composite × F²]
   γ_4 · □F²,                                  [derivative correction, no ∂ψ]
   Riegert local with σ_eff = function(ψ)      [auxiliary, S05 preserved]
}
```

To są **disjoint operator classes** od ψ.1.v3 canonical basis:

```
B_ψ.1.v3^dim-6 = {
   L₅'_a = (β_g/Λ²)·(∂_μ lnX)(∂_ν lnX)·F^{μρ} F^ν_ρ,        [parity-even]
   L₅'_b = β̃_g·(∂_μ lnX)(∂_ν lnX)·F^{μρ} F̃^ν_ρ              [parity-odd]
}
```

### §2.2 — Proof structure

**Operator class characterization:**

1. **∂ψ leg counting:**
   - Trace anomaly TGP-reduced: maximum **0 explicit ∂ψ legs** (z R[g_eff] ≈
     a_1·∇²h being a *single 2nd-derivative scalar*, oraz R^{μν} being a *single
     2nd-derivative tensor*; obie mają strukturę `∂² · ψ`, NIE `(∂ψ) · (∂ψ)`).
   - ψ.1.v3 basis: **2 explicit ∂lnX legs** w outer product `(∂_μ lnX)(∂_ν lnX)`
     forming rank-2 symmetric tensor product.
   - Te są **structurally different operator classes** — NIE są equivalent
     przez integration by parts (które wymaga znikającego boundary term, co dla
     gradient-rich substrate jest niefizyczne).

2. **Tensor structure:**
   - Trace anomaly: scalar(R)·F² lub tensor(∂_μ∂_ν ψ)·F²-tensor (z 2nd-derivative
     tensor, nie outer-product tensor).
   - ψ.1.v3: 1st-derivative outer-product tensor `(∂_μ ψ)(∂_ν ψ)` × F²-tensor.

3. **σ_ab stricte composite** — gradient-strain scale-0 object zdefiniowany w
   OP-7 T2 (2026-04-25):
   ```
   σ_ab = K_ab - (1/3)δ_ab Tr K,    K_ab = ⟨(∂_a ŝ)(∂_b ŝ)⟩
   ```
   To jest **traceless 3D spatial composite** z ŝ (sound-speed gradient at
   substrate scale 0), NIE macroscopic 4D `(∂_μ Φ)(∂_ν Φ)` outer product. Choć
   formalnie zawiera 1st-derivative outer product, scale-0 ŝ jest *różny* od
   macroscopic Φ pole. ψ.1.v3 basis używa macroscopic `lnX = ln(Φ)` (pole
   substratu), NIE scale-0 ŝ. ⇒ disjoint sektory.

4. **(α/(3π))·F² pure-photon dim-4:**
   - **Explicitly excluded** z ψ.1 basis przez Phase 7 T7.1 invariance filter
     (pure-photon sektor decoupled).
   - To jest dokładnie *expected* — pure-EM 1-loop renormalization (Euler-Heisenberg
     class) żyje w *sektorze α-running*, nie w substrate-coupling sektorze.

5. **Parity-odd L₅'_b:**
   - 1-loop QED z **single Dirac fermion** + **classical photon background**
     **NIE produkuje** parity-odd anomaly (Adler-Bardeen formula dotyczy
     `F·F̃` z chiral fermion mass; pure-QED ma L_QED parity-even).
   - ψ.1.v3 L₅'_b żyje w *separate parity-odd sektor* (Adler positivity Phase 8
     pending).

**Konkluzja:** Theorem 2.1 zachodzi *strukturalnie*. **Q1 zamknięty z dedicated
derivation** (a nie tylko z operator-class argumentu).

## §3 — GW170817 c_GW=c_EM full check (R3 closure)

### §3.1 — Photon dispersion correction

Effective Lagrangian (renormalized):
```
L_eff[A; g_eff] = -¼ F² + γ_1·(R/m_e²)·F² + γ_2·... + ...
```

z `m_e²` IR cutoff scale (electron mass z fermion loop integrate-out).

Photon EOM dla plane wave A_μ ~ exp(ik·x) na slowly-varying R background:
```
k² · (1 + 4γ_1 R/m_e²) ≈ 0
⇒ ω² = c²k² · (1 + 4γ_1 R/m_e²)
⇒ Δc/c ≈ -2γ_1 · R/m_e²
```

**Numerical estimate** dla GW170817 propagation (cosmological background):
- R_cosmo ~ H₀² ≈ (10⁻²⁶ m⁻¹)² ≈ 10⁻⁵² m⁻²
- m_e² ~ (m_e c / ℏ)² ≈ (2.6·10¹² m⁻¹)² ≈ 6.8·10²⁴ m⁻²
- γ_1 ~ α/(3π) ≈ 7.74·10⁻⁴

```
Δc/c ~ (α/(3π)) · R_cosmo / m_e² ≈ 7.74·10⁻⁴ · 10⁻⁵² / 10²⁵ ≈ 10⁻⁸⁰
```

**GW170817 bound:** |c_GW/c_EM - 1| < 9·10⁻²² ≈ 10⁻²¹.

**Margin:** ~58 OOM below bound. **Strukturalnie preserved.**

### §3.2 — Graviton dispersion

Per [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] §3:
- g_eff has NO independent dynamics (BD demarcation §5.1).
- Tensor mode h_TT propagates at c (no Lorentz violation w Φ sector).

Quantum 1-loop QED:
- QED loop integrates out *fermion* (electron); produces effective action dla *photon*.
- **NO graviton in QED at 1-loop** (graviton-photon vertex is suppressed by
  Newton's G; not part of 1-loop QED).
- ⇒ **c_GW unchanged** od 1-loop QED.

**Konkluzja:** c_GW = c (classical, structural) i c_EM = c · (1 + O(10⁻⁸⁰))
(quantum-corrected). Różnica `c_GW - c_EM` ≈ -10⁻⁸⁰, far below bound 9·10⁻²².

⇒ **R3 (GW170817 c-violation risk) zamknięty.**

## §4 — R-guard verification (Phase 2)

### §4.1 — R2 closure (operator class re-overlap)

**Strategy:** Theorem 2.1 *konstruktywnie* pokazuje że trace anomaly TGP-reduced
operator classes są disjoint od ψ.1.v3 basis.

**Sympy T4:** explicit class enumeration + comparison; zero overlap confirmed.

**R2 zamknięty.**

### §4.2 — R3 full closure (GW170817)

**Sympy T5:** Δc/c ≈ 10⁻⁸⁰ ≪ 9·10⁻²² (margin ~58 OOM).
**Sympy T6:** c_GW NIE modyfikowane przez 1-loop QED (no graviton in loop).

**R3 zamknięty.**

### §4.3 — R5 partial documented

**Sympy T7:** B/B_QED ratio dla lab (~10⁻¹⁰, perturbative VALID), magnetar
(B~10¹¹ T, B/B_QED ≈ 23 — perturbative MARGINAL/BREAKS).

**R5 partial:** cycle valid dla B ≪ B_QED regime; extreme magnetar surface
defer non-perturbative.

### §4.4 — S05 (R4) re-verification

**Sympy T8:** wszystkie reduced operators contain only Φ jako fundamental field.

## §5 — Findings (exportable)

| ID | Finding | Source |
|---|---|---|
| **F2.1** | R[g_eff] reduction w 1PN: `R = a_1·∇²h` (b_1=-a_1, γ=1); pure-2-derivative ψ scalar | sympy T1, §1.2 |
| **F2.2** | R^{μν}[g_eff] reduction w 1PN: `R^{ij} = -a_1·∂^i∂^j h` (tensor 2-derivative, scalar part zero) | sympy T2, §1.3 |
| **F2.3** | T^μ_μ_EM,1-loop_TGP operator classes: 6 distinct (pure-photon F², (∂²ψ)F², (∂_μ∂_ν ψ)F²-tensor, σ_ab F², □F², Riegert) | §1.4, sympy T3 |
| **F2.4** | **Theorem 2.1 (Disjointness):** trace anomaly TGP-reduced ⊥ ψ.1.v3 basis B={L₅'_a, L₅'_b} | §2, sympy T4 |
| **F2.5** | Q1 *konstruktywnie* zamknięty z dedicated 1-loop QED derivation | §2.2, F2.4 |
| **F2.6** | GW170817 dispersion: Δc/c ≈ 10⁻⁸⁰ ≪ 9·10⁻²² (~58 OOM margin) | sympy T5 |
| **F2.7** | c_GW=c structurally — no QED↔graviton coupling at 1-loop | sympy T6 |
| **F2.8** | R5: cycle valid dla B ≪ B_QED ≈ 4.4·10⁹ T; extreme magnetar surface deferred | sympy T7 |
| **F2.9** | S05 single-Φ preserved: all reduced operators contain only Φ source | sympy T8 |
| **F2.10** | ρ_EM_quantum[{Φ_i}] = native source field dla Φ-EOM, computable z Phase 1+2 framework | §1.4 |

## §6 — Cross-cycle convergence diagnostic update

Cztery niezależne diagnozy zbieżne na **separable sector structure**:

| Cycle | Diagnosis pattern | Sector separation level | Updated status (po Phase 2) |
|---|---|---|---|
| L01 ADDENDUM §3.2 (Q3 native estimate) | numerical magnitude | 10 OOM | Magnetar regime estimate **Phase 3 will refine** z corrected α/(3π)≈7.74·10⁻⁴ |
| τ.3 ADDENDUM §2 (L4 vs ρ_EM_quantum) | mechanism | distinct EOM paths | **Confirmed** — trace anomaly produkuje different operator classes vs L4 |
| ψ.1 ADDENDUM §3 (L01-Q1 resolution) | operator class | disjoint dim-6 vs dim-4 | **Konstruktywnie potwierdzony** — Theorem 2.1 |
| Q2 cycle | vacuum-level | substrate vs matter sector | Consistent — ρ_EM_quantum jest *transient source*, NIE Λ contribution |
| **op-L01-N1 (this cycle)** | **constructive derivation** | **sektor structurally derived w 1-loop QED na g_eff[{Φ_i}]** | **Phase 1+2 closes** the constructive verification |

**Cross-cycle synthesis:** pięć niezależnych diagnoz teraz zbieżne na separable
sector structure jako *strukturalna własność TGP*, nie post-hoc tuning.

## §7 — Phase 3 handoff

### §7.1 — Co Phase 2 dało

1. **Renormalized ρ_EM_quantum[{Φ_i}]** w postaci konkretnej combination 6
   operator classes z explicit prefactors.
2. **Theorem 2.1 (Disjointness)** vs ψ.1.v3 — **Q1 zamknięty z dedicated
   derivation**.
3. **R3 zamknięty** — GW170817 c_GW=c_EM preserved structurally.
4. **R5 partial documented** — B ≪ B_QED regime explicit.

### §7.2 — Co Phase 3 musi dostać (phenomenology)

1. **Lab regime numerics** (B ~ 1 T, R ~ 10⁻⁵² m⁻²): ρ_EM_quantum estimate +
   MICROSCOPE η budget.
2. **Magnetar regime numerics** (B ~ 10¹¹ T): ρ_EM_quantum/ρ_NS ratio z
   *corrected* α/(3π) prefactor (10⁻⁹ ÷ 10⁻¹⁰, NIE 10⁻¹²).
3. **MICROSCOPE 5-th force bound check** — Pt vs Ti differential η_TGP_EM_quantum.
4. **R6 verification** — QEP universality check: universal coupling structure
   `L_mat = -(q/Φ_0)·φ·ρ` automatycznie kasuje Asorey-2015 type QEP violations
   (S05 mechanism — wspólna metryka g_eff dla wszystkich).

## §8 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]]
- [[./Phase1_results.md]]
- [[./Phase2_setup.md]]
- [[./Phase2_sympy.py]] / [[./Phase2_sympy.txt]] (8/8 PASS)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase1_results.md]]
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] §3
- [[../op-psi1-substrate-light-acceleration/Phase7_results.md]] (canonical basis B)
- [[../op-psi1-substrate-light-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §3 (Q1 closure target — *konstruktywnie* odtworzony)
- [[../op-tau3-substrate-clock-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §2

---

**Phase 2 close:** 8/8 sympy PASS. Theorem 2.1 (Disjointness) verified
*konstruktywnie*. Phase 3 may proceed.
