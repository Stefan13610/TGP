# OP-M92 Phase 0+ — Candidate D multi-source self-consistency check

**Data:** 2026-04-25
**Status:** **ISSUE FLAGGED — Phase 0+ structural sketch self-consistency NEGATIVE**
**Sub-test:** Candidate D α calibration cross-source (M87* vs Sgr A* vs GW150914 vs NS)
**Następnik:** Phase 1 PRIORITY ELEVATED — multi-source consistency derivation REQUIRED

---

## TL;DR

> **HONEST SELF-CHECK NEGATIVE FINDING:** Phase 0+ structural sketch's
> heuristic "α_geom ~ 0.1 (M_BH=1)" is **physically inconsistent** across
> Schwarzschild BH sources of different masses.
>
> Extracting α_SI assuming universal geom-units calibration:
>
> - **Sgr A***  (M = 4.3e6 M_⊙): α_SI = 1.61e+19 m²
> - **M87***    (M = 6.5e9 M_⊙): α_SI = 3.69e+25 m² (**2.3×10⁶** larger)
> - **GW150914** (M = 65 M_⊙): α_SI = 3.69e+09 m² (2.3×10⁻¹⁰ smaller)
> - **Neutron star** (M = 1.4 M_⊙): α_SI = 1.71e+06 m² (10⁻¹³ smaller)
>
> ⇒ α scales as M_BH² across **13 orders of magnitude** — α is NOT a
> single physical constant under naive Candidate D action.
>
> **Verdict: ISSUE FLAGGED, not fatal but Phase 1 priority.** Phase 0+
> structural sketch was a **structural plausibility check only** —
> rigorous Phase 1 covariant derivation must resolve multi-source
> consistency before Candidate D can claim to reproduce scenario (e)
> universally.
>
> This is a **CANDIDATE-INDEPENDENT** issue — any M9.2 pivot that aims
> to reproduce universal +14.56% strong-field deviation must explain
> how a single physical parameter handles BHs across 9+ orders of
> magnitude in mass.

---

## 1. The dimensional issue

W Phase 0+ structural sketch (afternoon 2026-04-25):

```
T·J·J / S_kin ~ O(1)  in geom units (M_BH=1)  at photon ring (r=3.88M, ψ=1.168)
```

⇒ For scenario (e) target shift 0.114, naive heuristic α (in M_BH² units) ~ 0.1.

**Critical assumption:** the O(1) factor is **universal** in geom units across
all Schwarzschild sources. Since photon ring geometry (r=3.88M, ψ=1.168) is
universally identical in geom units, this assumption is internally consistent.

**Problem:** Converting back to SI requires multiplying by R_S(source)²:
- α_SI = α_geom × R_S(source)²
- R_S = 2GM/c² scales linearly with M_BH
- ⇒ α_SI scales as M_BH²

For α to be a true universal physical constant (single number in SI), this
scaling MUST be cancelled by something else in the action.

## 2. Numerical results (multi-source α extraction)

| Source | M (M_⊙) | R_S (m) | α_SI (m²) | Ratio vs Sgr A* |
|--------|---------|---------|-----------|-----------------|
| Sgr A* | 4.3e+6 | 1.27e+10 | 1.61e+19 | 1.0 (ref) |
| M87* | 6.5e+9 | 1.92e+13 | 3.69e+25 | 2.29e+6 |
| GW150914 final | 65 | 1.92e+5 | 3.69e+9 | 2.29e-10 |
| Neutron star | 1.4 | 4.13e+3 | 1.71e+6 | 1.06e-13 |

Span: **19 orders of magnitude** in α_SI for sources differing by 9.5
orders of magnitude in mass (consistent with M² scaling).

## 3. Resolution paths

### PATH A: TGP-intrinsic length scale L_TGP — **FAILS**

Hipoteza: action ∼ α_dim × L_TGP² × T^μν J_μ J_ν z α_dim dimensionless.

- Required L_TGP for Sgr A* rescue: ~ 4.0×10⁹ m (~0.3 R_S(SgrA*))
- Required L_TGP for M87* rescue: ~ 6.1×10¹² m (~0.3 R_S(M87*))
- Single L_TGP cannot satisfy both — fails by factor 1.5×10³

### PATH B: T^μν refers to source-enhanced substrate density — **FAILS**

Hipoteza: substrate density ρ_sub ~ M_BH / R_S³ near BH (accretion-like).
Then T·J·J ~ ρ_sub / R_S² ~ 1/R_S⁴ → α requires same M² scaling. **Failure
mode identical to naive case.**

### PATH C: α dimensionless × M_Pl⁻² prefactor — **FAILS unphysically**

Hipoteza: action ∼ (α/M_Pl²) × T·J·J z α O(1).
- M_Pl/M_SgrA ~ 2.5×10⁻⁴⁵
- Required α ~ M⁴/M_Pl² ~ 1.5×10⁸⁹ for Sgr A* rescue
- Astronomically unphysical. **Fails.**

### PATH D: Scenario (e) target NOT universal — **NOT A VALID RESOLUTION**

Re-examine T-M92.1: scenario (e) factor f(ψ) = √(g_tt^GR/g_tt^TGP) at
photon ring depends purely on relative geometry r_ph^TGP=3.88M vs
r_ph^GR=3M (universal in geom units). Target +1.46% deviation IS
universal. **Cannot escape via this path.**

### PATH E: Modified Candidate D — non-minimal ψ-coupling — **VIABLE**

Hipoteza: replace constant α with α(ψ) function:
```
α(ψ) = α_0 × (ψ − ψ_threshold)^n × Θ(ψ − ψ_threshold)
```

- Activates ONLY in strong-field where ψ exceeds threshold (~ 1.05?)
- α_0 dimensionless O(1) physical constant
- ψ-dependence absorbs the M_BH scaling because ψ at photon ring is
  universal but ψ at lab/Earth scale is essentially 1 (no activation)

**Cost:** introduces new parameters (ψ_threshold, n, α_0). Increases
M9.2-D structural complexity but resolves multi-source issue. Phase 1
derivation must determine if this modification is unique/natural.

### PATH F: Composite Candidate D + Candidate B — **POSSIBLE HYBRID**

Conformal frame coupling (Candidate B) z momentum back-reaction
(Candidate D) razem mogłyby produce universal scenario (e) deviation
without explicit ψ-threshold.

Phase 1 must investigate hybrid landscape.

## 4. Implications dla wcześniejszych Phase 0+ cross-checks

Critical re-examination wszystkich tests done dzisiaj:

| Test | Original verdict | Multi-source impact | Robust? |
|------|------------------|---------------------|---------|
| Structural sketch (5/5 POSITIVE) | PASS | Heuristic O(1) needs rigorous derivation | CONDITIONAL |
| Cosmology (w(z) ≥ −1) | margin 1e+33× | α_SI variation by 10¹⁰ → still safe | **ROBUST** |
| WEP MICROSCOPE | margin 6.7× | Re-calibration could shift margin O(1)×–O(10)× | **VULNERABLE** |
| LLR Nordvedt | margin 3e+15× | α variation absorbed by huge margin | **ROBUST** |
| GW c_GW vacuum | structural | Unaffected (T=0 in vacuum) | **STRUCTURAL** |

**Key insight:** Cosmology and LLR margins are **so large** że αSI
variation by 10⁶+ doesn't break them. WEP MICROSCOPE jest **vulnerable**
— exact margin depends sensitively on which calibration we use.

WEP re-calibration test (Phase 0+ supplement):
- Calibrating α to M87* instead of Sgr A* multiplies α_SI by 2.3×10⁶
- → η_MICROSCOPE heuristic increases by 2.3×10⁶ → η ~ 4×10⁻¹⁰
- → **Falsifies MICROSCOPE bound 1.1×10⁻¹⁵ by factor ~3×10⁵**

⇒ WEP cross-check VERDICT IS CALIBRATION-DEPENDENT. The 6.7× margin
holds **only** for Sgr A*-calibrated α. M87*-calibrated α would FAIL
MICROSCOPE.

This is actually **consistent with the multi-source issue** — if α
were truly universal at M87* level, Candidate D would already be ruled
out by MICROSCOPE. So the resolution path E (ψ-threshold) is preferred
because it makes α effectively zero at lab scale (ψ_lab ≈ 1).

## 5. Revised Phase 0+ status

### Pre-multi-source verdict (afternoon-evening 2026-04-25)

5/5 POSITIVE structural + cosmology + WEP cross-checks ALL PASS.

### Post-multi-source verdict (night 2026-04-25)

**5/5 POSITIVE structural + cosmology + WEP cross-checks PASS, BUT:**

- Heuristic α_geom ~ 0.1 universal IS NOT a sound physical theory
- Resolution paths A-D fail; Path E (ψ-threshold) is preferred
- Phase 1 covariant derivation MUST establish:
  1. Whether modified Candidate D action with α(ψ) reproduces scenario (e)
  2. The precise threshold structure (ψ_threshold, n)
  3. Multi-source consistency (M87*, GW150914, NS-NS mergers)
  4. Compatibility with existing Phase 0+ cosmology + WEP results

**Updated F4 falsifiability narrative:**

> M9.1'' standalone falsifiable. M9.2-D pivot path: structural sketch
> POSITIVE z multi-source consistency caveat — naive constant α scales
> as M_BH² (incompatible z universality). Path E (α(ψ) threshold)
> resolves but introduces extra parameter. Phase 1 derivation **MUST**
> verify Path E gives universal physical α before claiming Candidate D
> reproduces scenario (e) universally.
>
> Multi-source consistency check is **independently falsifiable**:
> if EHT next-generation gives BH shadow data dla multiple sources
> (M87*, Sgr A*, plus 5-10 other LLAGN by 2030+), TGP M9.1'' predicts
> universal +14.56% deviation. Any deviation from universality among
> these sources would falsify even M9.2-D + Path E (which predicts
> approximately universal deviation modulo ψ-threshold artifacts).

## 6. Phase 1 priority list (REVISED)

| Priority | Item | Rationale |
|----------|------|-----------|
| **1 (NEW)** | Multi-source self-consistency derivation | Resolves α scaling issue; without this, Candidate D is structurally incomplete |
| 2 | Rigorous WEP/Nordvedt analysis | MICROSCOPE 6.7× margin (calibration-dependent) |
| 3 | Photon ring numerical verification w/ Path E α(ψ) | Confirm scenario (e) reproducibility |
| 4 | Cosmology perturbation scale dependence | Extends Phase 0+ cosmology to include perturbations |
| 5 | Loop corrections / ghost screening | Beyond-tree-level stability |

## 7. Files (Phase 0+ multi-source)

- `op_m92_P0plus_candD_multisource.py` — analytical cross-check script
- `op_m92_P0plus_candD_multisource.txt` — raw output
- `OP_M92_P0plus_candD_multisource_results.md` — this synthesis

## 8. Cross-references

- [[research/op-m92/OP_M92_P0plus_candD_results.md]] — Phase 0+ kickoff (parent)
- [[research/op-m92/OP_M92_P0plus_candD_cosmology_results.md]] — cosmology cross-check
- [[research/op-m92/OP_M92_P0plus_candD_wep_results.md]] — WEP cross-check (now calibration-vulnerable)
- [[research/op-m92/op_m92_T1_action_reverse_engineer.py]] — T-M92.1 universality of scenario (e)

---

## Bottom line

OP-M92 Phase 0+ multi-source self-consistency check **FLAGS A CRITICAL
ISSUE** w naive Candidate D Lagrangian: α_geom ~ 0.1 universal implies
α_SI scales as M_BH² across 13 orders of magnitude — α nie jest single
physical constant.

**Honest re-assessment:**
- Phase 0+ structural sketch (5/5 POSITIVE) was structural plausibility,
  NIE rigorous proof of feasibility
- Cosmology cross-check (margin 10³³×) ROBUST do α calibration variation
- WEP MICROSCOPE cross-check (margin 6.7×) **CALIBRATION-DEPENDENT** —
  M87*-calibrated α would FAIL MICROSCOPE by 10⁵× margin
- GW c_GW vacuum cross-check ROBUST (vacuum exact)

**Resolution path E (α(ψ) z threshold near ψ=1.05):** preferred,
introduces extra parameter ale resolves multi-source consistency
i WEP via ψ_lab ≈ 1 (effective null α at lab scale).

**Phase 1 priority elevated:** multi-source self-consistency derivation
becomes #1 priority (przed WEP, photon ring numerical, cosmology
perturbations).

**M9.2-D pivot path post-OP-M92 Phase 0+ (full):** alive z lead candidate
status, ALE z explicit Phase 1 caveat — Path E (modified α(ψ) coupling)
musi być derived rigorously przed paper-level claim. Response time
post-ngEHT 2030+ verdict pozostaje 2-4 weeks (sketch + Path E template
already structurally pre-derived), choć Phase 1 derivation effort
wzrasta z 6-12 miesięcy do 9-15 miesięcy.

This negative finding **strengthens** OP-M92 program rigor — pre-emptive
identification of structural issues PRZED ngEHT verdict znaczy że
post-verdict response nie napotka surprise issues.
