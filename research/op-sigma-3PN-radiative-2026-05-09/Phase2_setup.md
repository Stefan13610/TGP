---
title: "Phase 2 setup — PN-counting lock + dimensional analysis + Path A direct strategy"
date: 2026-05-09
parent: "[[./README.md]]"
type: phase-setup
phase: 2
predecessor: "[[./Phase1_results.md]]"
classification: SETUP_PRE_COMPUTATION
critical_decisions:
  - "Path A direct (cleaner) vs Path B Hadamard PF (HANDOFF original plan)"
  - "PN-order resolution: U vs U² scaling (Phase 1 §5.5 ambiguity)"
  - "Dimensional bookkeeping: c_0/(Φ_0²·c²)·σ explicit units"
---

# Phase 2 setup — Pre-computation locks

## §0 — Strategic decision: Path A direct

### §0.1 — HANDOFF original plan vs alternative

**HANDOFF §3 plan (Path B near-field):**
1. Hadamard partie-finie ∫σ d³x dla binary (UV-divergent at particle positions)
2. Renormalize Q^σ_ab(t) tensor moment
3. Far-field σ_ab ~ d²Q^σ/dt² · 1/r
4. TT-projection at observer

**Concerns z analizy ekspert §C.1-C.2:**
- C.1: Hadamard PF jest scheme-dependent powyżej 2PN (Blanchet 2002+ używa dim-reg)
- C.2: Finite part po PF zależy od cut-off scheme; może wprowadzić unphysical parameter
- C.3: σ jako (∂Φ)² jest **near-field bilinear**, nie evident źródło 1/r radiation
- PN-order ambiguity (U vs U²) Phase 1 §5.5 wynika właśnie z tego mętnego sourcing

### §0.2 — Path A direct (RECOMMENDED PRIMARY)

**Z OP-7 T3.1 LOCK** (closure_2026-04-26 i op-7/OP7_T3_results.md):
```
L_σ = -(1/4)(∂_μ σ_ab)(∂^μ σ^ab) - (1/2)m_σ² σ_ab σ^ab - (ξ_eff/2)·σ_ab·T^{ab,TT}
```

EOM (variational principle):
```
□ σ_ab + m_σ² σ_ab = -ξ_eff · T_ab^TT
```

**Co to daje:**
1. σ_ab jest polem effective z kanoniczną kinematyką (ghost-free, OP-7 T3.3)
2. **Source jest T_ab^TT** — TT-projection of matter stress-energy tensor (NIE (∂Φ)² composite)
3. Standard retarded Green function + far-field expansion daje **bezpośrednio** σ_ab(observer)
4. Zero Hadamard ambiguity: T_ab^TT dla matter point particles jest dobrze zdefiniowane

**Equivalence z Path B (OP-7 T3.1 + closure 2026-04-26):**
- Path A == Path B **structurally** (T-PB closure 11/11 PASS)
- m_σ² = 2m_s² **derived** (nie postulowane)
- ξ_eff = G·Φ_0² (T3.4 LOCK)
- σ_ab w Path A jest effective opis — same physics as Path B operatorowy

**Wniosek:** Path A daje ten sam wynik fizyczny z **kompletnie różną** (i czystszą)
strukturą obliczeniową. Phase 2 używa Path A jako primary; Path B Hadamard
jako cross-check sub-cycle (jeśli czas).

## §1 — PN-counting lock (resolution U vs U² ambiguity)

### §1.1 — Krytyczna obserwacja

W Path A EOM:
```
□σ_ab = -ξ_eff·T_ab^TT  (massless limit, m_σ ≪ ω)
```

Source `T_ab^TT` to **TT-projection mass-current quadrupole** z source matter:
```
T_ij^TT(point particles) = TT[Σ m_a v_a^i v_a^j δ³(x - x_a)]
```

Dla circular binary equal-mass:
```
T_xx = m·v²·cos²(ωt/2)·δ³(x-x_1) + ...   (oscillates at 2ω)
```

**Time-derivative structure:** d²T_ij^TT/dt² oscillates at 2ω (twice orbital frequency),
**identycznie jak w GR mass quadrupole**.

### §1.2 — PN order σ-induced TT — NOT 3PN, NOT U²

**Phase 1 §5.5 i HANDOFF claimed σ-induced TT entry at 3PN+ z scaling U lub U².**

**Path A correction:** EOM source `-ξ_eff·T_ab^TT` jest **LINEAR** w T_ab. Więc σ_ab
amplitude jest **leading order** z dispersion analog standard quadrupole formula:

```
σ_ab^far(observer) = -(ξ_eff/(4πr·c⁴)) · d²Q^M_ij^TT(t-r/c)/dt²
```

To jest **LEADING** quadrupole order (analog 2.5PN dla standard GR), **NIE**
suppressed U² ani U.

**PN counting RIGOROUS:**

| Quantity | PN order | GR analog |
|---|---|---|
| Mass quadrupole Q^M_ij | 0PN (leading) | Same |
| d²Q^M/dt² (radiation amplitude) | 0PN (leading) | Same |
| T_ij^TT source | 0PN | Same |
| σ_ab from □σ = -ξ T^TT (Path A) | 0PN amplitude | Equivalent |
| h_TT^σ at observer | 0PN amplitude | Equivalent |

**Wniosek:** σ-induced TT amplitude jest **same PN order as GR mass quadrupole**.
Różnica between TGP h_TT^σ i GR h_TT^GR jest w **prefactorze**, nie w PN suppression.

### §1.3 — Why Phase 1 §5.5 was misled

Phase 1 §5.4-5.5 worked z formula σ_ij = (∂_iΦ)(∂_jΦ) - trace, computed Q^σ ~ ∫(∂Φ)².
To daje (∂Φ)² ~ U² scaling **dla NEAR-FIELD bilinear**, ale to jest **fundamentally
different quantity** od σ_ab jako propagating field z source T_ab^TT.

**Path A vs Path B distinction:**
- Path B (operatorowy): σ_ab = ⟨(∂_aŝ)(∂_bŝ)⟩ — composite operator
- Path A (effective): σ_ab is propagating effective field z own EOM

Equivalent **fizycznie** ale Path A jest cleaner dla radiation problem.

**Critical insight:** the U² counting from (∂Φ)² applies do NEAR-FIELD bilinear,
but at far-field σ_ab solving Path A EOM has **leading-order** amplitude.

### §1.4 — Sphere-averaged ≠ observer amplitude (anti-pattern check)

§5 anti-pattern §5.1 z HANDOFF: "NIE conflate sphere-averaged z observer-amplitude".

**Phase 2 will compute σ_ab^TT(observer) at SPECIFIC angular position n**, NIE
sphere-integrated. Will use explicit angular structure z d²Q^M_ij^TT(n)·n^i·n^j.

## §2 — Dimensional bookkeeping (resolution C.4 concern)

### §2.1 — Units of relevant quantities

Working w SI-compatible units (z explicit c, G, Φ_0):

| Quantity | Dimension | Value/Status |
|---|---|---|
| G | [L³/(M·T²)] | Newton constant |
| c | [L/T] | speed of light |
| Φ_0 | [v²] = [L²/T²] | substrate background scalar (from T3.4) |
| K_1 | dimensionless | linearized normalization |
| q | dimensionless or [Φ_0]^(-1) | matter coupling charge |
| m_σ² | [1/T²]·[1/L²]·c² = [1/L²] (mass²) | composite mass squared |
| m_s² | [1/L²] | substrate mass squared |
| ξ_eff = G·Φ_0² | [L³/(M·T²)]·[L⁴/T⁴] = [L⁷/(M·T⁶)] | dimensionful coupling |

**T_ab^TT dimension:** [energy density] = [M/(L·T²)]

### §2.2 — Path A EOM dimensional check

```
□σ_ab + m_σ²σ_ab = -ξ_eff · T_ab^TT
```

LHS dimension: [σ]/[L²] (z □ = ∂²/∂x²)
RHS dimension: [ξ_eff]·[T] = [L⁷/(M·T⁶)]·[M/(L·T²)] = [L⁶/T⁸]

Therefore: [σ] = [L⁴/T⁶]·[L²] = [L⁶/T⁶]?? Hmm that's strange.

Actually z OP-7 T3.4 chain: ξ_eff = G·Φ_0² has structure [L³/(M·T²)]·[L⁴/T⁴] for SI units...

Actually let me redo this w **geometric units c=G=1** (standard dla GR PN):

In c=G=1:
- [Φ_0] = dimensionless (or [length]^-2 depending on convention)
- [σ_ab] = dimensionless metric perturbation analog
- [T_ab] = [length]^-2
- [m_σ²] = [length]^-2
- [□] = [length]^-2

Path A EOM dimensionally consistent if [ξ_eff]·[T] = [σ]·[length]^-2.
[ξ_eff·T] = [σ]/[length]^2 ⟹ [ξ_eff] = [σ]·[length]^0/[T] = dimensionless × [length]² (since T~[length]^-2)
Actually if σ dimensionless, then ξ_eff has [length²]. With ξ_eff ~ G·Φ_0² and Φ_0 ~ 1 (dimensionless), ξ_eff ~ G ~ [length²] ✓ in c=G=1 (since G has Mass scale via energy = mass).

Hmm dimensional analysis in GR units is subtle. **Let me just track factors symbolically in sympy** and verify dimensional consistency post-hoc.

### §2.3 — Emergent-metric coupling

```
δg_eff^ij ⊃ (c_0/Φ_0²) · σ^ij  (in c=1 units; c² factor reabsorbed)
```

This is the coupling that determines what LIGO detector measures.
**c_0 dimensionless** (geometric constant, ≈ 4π z OP-7 T3.4 chain).

### §2.4 — Target dimensionless ratio

**Final observable:** h_TT^σ(observer) / h_TT^GR(observer)

Both have units of dimensionless metric perturbation. Their ratio is dimensionless.

**Hipoteza dla Phase 2 sympy:**
```
h_TT^σ / h_TT^GR = (c_0 · ξ_eff) / (8π · G · Φ_0²) · (numerical factor from kinematics)
                 = c_0 / (8π) · (ξ_eff / (G·Φ_0²)) · κ_kin
                 = c_0 / (8π) · 1.06 · κ_kin   [GW150914 calibration ξ_eff/G·Φ_0² ≈ 1.06]
```

Z c_0 ≈ 4π: ratio = (4π)/(8π) · 1.06 · κ_kin = 0.53 · κ_kin

**Critical question:** what is κ_kin (additional kinematic factor z TT-projection,
orbital averaging, or κ_σ)? Phase 2 sympy MUST compute explicitly.

If κ_kin ≈ κ_σ ≈ 1/(3π): ratio ≈ 0.53/3π ≈ 0.056 (~5.6%, borderline LIGO).
If κ_kin ≈ 1: ratio ≈ 0.53 (~53%, **massively violates** LIGO 5%).
If κ_kin ≈ 1/(2π): ratio ≈ 0.084 (~8.4%, violates).

**This is THE key Phase 2 calculation.** Absolute dimensional answer determines
whether framework PASSES or FAILS LIGO bound.

## §3 — Sympy structure plan

### §3.1 — Sub-section 1: Path A EOM solution

Massless limit (m_σ ≪ ω, justified per OP-7 T6 decoupling):
```
□σ_ab = -ξ_eff · T_ab^TT
```

Retarded Green function (massless wave equation):
```
σ_ab(x,t) = (ξ_eff/(4π)) · ∫ d³y T_ab^TT(y, t-|x-y|/c) / |x-y|
```

Far-field expansion:
```
σ_ab^far(r→∞,t) = (ξ_eff/(4π·r·c⁴)) · ∫ d³y T_ab^TT(y, t_ret) · (1 + O(c·t_ret/r))
                ~ (ξ_eff/(4π·r·c⁴)) · d²Q^M_ij^TT(t-r/c)/dt²
```

Where Q^M_ij = ∫ ρ x_i x_j d³x = mass quadrupole.

**Sympy task:**
- Verify Green function solution dimensional consistency
- Verify multipole expansion formula
- Compute d²Q^M_TT/dt² for circular binary explicit

### §3.2 — Sub-section 2: TT-projection at observer

Observer direction n̂. TT-projector P^ij_kl applied to σ_ab^far at observer.

Critical: σ_ab is already symmetric and traceless (verified Phase 1 §3).
TT-projection further projects orthogonal to direction n̂.

```
σ_ab^TT(observer) = Λ^ij_kl · σ_kl^far(observer)
                  = Λ^ij_kl · (ξ_eff/(4π·r·c⁴))·d²Q^M_kl^TT/dt²
```

**Sympy task:**
- Verify TT-projection of σ_ab is non-zero (since σ has tensor structure)
- Compute explicit form for binary z observer along z (face-on)
- Compute h_+ = (1/2)(σ_xx^TT - σ_yy^TT), h_× = σ_xy^TT

### §3.3 — Sub-section 3: Emergent-metric coupling

```
δg_eff_ij(observer) = (c_0/Φ_0²) · σ_ij(observer)   (c=1 units)
                    = (c_0/(Φ_0²·c²)) · σ_ij(observer)   (SI units)
```

**Sympy task:**
- Track all c, G, Φ_0 factors
- Substitute c_0 = 4π (geometric LOCK)
- Substitute ξ_eff = G·Φ_0² (T3.4 LOCK)
- Verify dimensional consistency end-to-end

### §3.4 — Sub-section 4: GR comparison + LIGO test

GR mass quadrupole formula:
```
h_ij^GR_TT(observer) = (2G/(c⁴·r)) · d²Q^M_TT/dt²
```

TGP σ-induced:
```
h_ij^σ_TT(observer) = (c_0/(Φ_0²·c²)) · σ_ij^TT(observer)
                    = (c_0/(Φ_0²·c²)) · (ξ_eff/(4π·r·c⁴)) · d²Q^M_TT/dt²
                    = (c_0·G/(4π·c⁶·r)) · d²Q^M_TT/dt²   [using ξ_eff=G·Φ_0²]
```

Wait — this gives different c-power dependence than GR. Let me re-check.

Actually: GR formula [M/(c⁴·r)]·d²Q/dt² ~ [M·L²/T²]/(L³/T²·L) = [M/L²]??
Hmm dimensions of GR formula actually:
- [G/(c⁴·r)] = [L³/(M·T²)]/[L⁴/T⁴]·[1/L] = [L³·T⁴/(M·T²·L⁴·L)] = [T²/(M·L²)]
- [d²Q/dt²] = [M·L²/T²]
- product = [T²/(M·L²)]·[M·L²/T²] = dimensionless ✓

For TGP σ amplitude (above): [c_0·G/(c⁶·r)]·[d²Q/dt²]
- [G/(c⁶·r)] = [L³/(M·T²)]/[L⁶/T⁶]·[1/L] = [L³·T⁶/(M·T²·L⁶·L)] = [T⁴/(M·L⁴)]
- [d²Q/dt²] = [M·L²/T²]
- product = [T⁴/(M·L⁴)]·[M·L²/T²] = [T²/L²] = NOT dimensionless ✗

Dimensional inconsistency! Something wrong w naive substitution. **This is exactly the
issue Phase 2 must resolve carefully.**

Possible resolutions:
1. ξ_eff = G·Φ_0² has hidden c factors in OP-7 T3.4 derivation
2. δg_eff coupling C(ψ)/(Φ_0²·c²) has different power of c
3. T_ab^TT in Path A EOM has different units

**Phase 2 sympy MUST resolve this** before claiming any numerical ratio.

### §3.5 — Sub-section 5: Sanity check c_0·κ_σ = 4/3 connection

Existing LOCK: c_0·κ_σ = 4/3 (joint cycles #1+#2, β_ppE = 0).

This is constraint at 2.5PN PHASE level. The question:
**Does h_TT^σ ratio computation use c_0 or c_0·κ_σ effective combination?**

If h_TT^σ amplitude factor is c_0 (= 4π): ratio ~ 0.5 (massive violation)
If h_TT^σ amplitude factor is c_0·κ_σ (= 4/3): ratio ~ 0.053 (borderline LIGO)
If h_TT^σ amplitude factor is c_0·κ_σ² or c_0²·κ_σ: different again

**Phase 2 sympy MUST derive which combination appears**, NIE assume.

## §4 — Anti-pattern compliance

### §4.1 — Pre-declared methodology

Phase 2 plan declared **before** computation:
- Path A EOM as primary tool (per §0.2)
- Dim-reg if Path B sub-cycle needed (per §C.2 concern)
- Specific numerical comparison: h_TT^σ / h_TT^GR ratio
- LIGO comparison: O3 polarization tests + 5% scalar bound
- Verdict: NO predefined "pass band", computed value compared to observation

### §4.2 — Verdict reformulation (per §C.6 concern)

Original HANDOFF §3.5 had pass/fail framing that was confirmation-bias-prone.
**Reformulated for Phase 2:**

| Computed h_TT^σ / h_TT^GR | Implication |
|---|---|
| < 1% | Within LIGO O3 noise, framework consistent (TT-dominant prediction) |
| 1-5% | Borderline; specific TGP prediction testable z O5+ data |
| 5-20% | Tension z LIGO O3 bounds; framework challenged |
| > 20% | Framework FALSIFIED at LIGO sensitivity level |

Verdict NIE jest "pass/fail at 5%", lecz **specific quantitative prediction**
compared do **specific quantitative observation**.

### §4.3 — Adversarial verification commitment

Per HANDOFF §5.4 + ekspert §C.7: Phase 2 result will undergo **independent
verification** (separate sympy script lub independent agent) **before**
propagating to amendment cascade.

This protocol JUST PRZYNIÓSŁ wartość w op-h-TT-calibration cycle (caught
Phase 3 cycle #3 error). Same protocol for Phase 2.

## §5 — Phase 2 verdict criteria

### §5.1 — STRUCTURAL DERIVED requirement

To classify Phase 2 jako STRUCTURAL DERIVED:
1. ✓ Path A EOM solution explicit (sympy)
2. ✓ Far-field σ_ab amplitude derived (no scheme ambiguity)
3. ✓ Dimensional consistency end-to-end verified
4. ✓ Numerical ratio h_TT^σ / h_TT^GR computed z explicit c_0, κ_σ, ξ_eff values
5. ✓ Cross-check z c_0·κ_σ = 4/3 constraint (consistent or contradictory)
6. ✓ LIGO comparison z proper observational bound

### §5.2 — STRUCTURAL_CONDITIONAL requirement

If 1-3 ✓ but 4-6 reveal dimensional inconsistency (per §3.4 concern):
- Phase 2 will be STRUCTURAL_CONDITIONAL
- Sub-cycle needed to resolve dimensional structure (1-2 sesji)
- Multi-session continuation

### §5.3 — STRUCTURAL_NO_GO requirement

If h_TT^σ / h_TT^GR > 20% z robust dimensional analysis:
- Framework violates LIGO bounds even at qualitative level
- Route A FAILS
- Pivot do Route B (nonlinear δΦ self-coupling)

## §6 — Time budget

**Phase 2 estimated z this setup: 1 sesja for explicit calculation.**

Sub-tasks:
- Phase2_sympy.py: ~30-50 sympy checks tracking dimensional + numerical
- Phase2_results.md: synthesis + verdict
- Adversarial sub-cycle (jeśli czas): independent verification

## §7 — Cross-references

- [[./Phase1_results.md]] — predecessor phase
- [[./HANDOFF_NEXT_SESSION.md]] — original strategy plan
- [[../op7/OP7_T3_results.md]] — Path A EOM LOCK + ξ_eff = G·Φ_0²
- [[../closure_2026-04-26/sigma_ab_pathB/results.md]] — Path B equivalence
- [[../op-c0-derivation-from-substrate-2026-05-09/Phase_FINAL_close.md]] — c_0 = 4π
- [[../op-kappa-sigma-2body-PN-2026-05-09/Phase_FINAL_close.md]] — κ_σ = 1/(3π)
- [[../op-emergent-metric-from-interaction-2026-05-09/Phase4_results.md]] — Phase 4 ansatz
- [[../op-h-TT-calibration-2026-05-09/Phase2_sympy.py]] — TT-projection identity LOCK

---

**End Phase 2 setup.** Path A direct strategy chosen as primary; PN-counting clarified
(σ-induced TT enters at LEADING quadrupole order, NOT 3PN); dimensional bookkeeping
plan declared. Ready for explicit sympy computation.
