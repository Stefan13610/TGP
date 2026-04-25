# OP-M92 Phase 1 — Rigorous Derivation Plan (Candidate D + Path E)

**Data utworzenia:** 2026-04-25 (post Phase 0+ FULL closure)
**Status:** PLAN ONLY — execution scheduled to **2026 Q3-Q4** OR **post-ngEHT 2030+ verdict**
**Estimated effort:** 9-15 miesięcy part-time analysis
**Trigger conditions:** ngEHT 2030+ confirms GR shadow at 1% Sgr A* (50% Bayesian probability)

---

## 1. Cel końcowy (Phase 1 deliverable)

**Paper-level rigorous M9.2-D + Path E framework** ready for:
- Full M9.2 axiom paper drafting (Phase 2)
- Submission to PRD / CQG within 2-4 weeks of ngEHT verdict
- Independent verification by competent reviewers

**Output structure:**
1. Rigorous covariant action `S_{M9.2-D}` with α(ψ) coupling
2. Modified Φ-EOM (Klein-Gordon-like z back-reaction term)
3. Photon ring numerical reproduction of scenario (e) target +1.46% deviation
4. Multi-source consistency proof (M87*, Sgr A*, GW150914 final, NS)
5. WEP/Nordvedt rigorous derivation with MICROSCOPE prediction
6. Cosmology perturbation theory (CMB Planck constraints + DESI w(z) shifts)
7. Stability analysis (no-ghost, beyond-tree-level renormalization)
8. 5th force null prediction at lab scales (ψ_lab → α effective ≈ 0)

---

## 2. Phase 0+ inputs (already shipped 2026-04-25)

**Sketch-level inputs ready for Phase 1 hardening:**
- Action ansatz: `S = S_{M9.1''} + ∫ α(ψ) T^μν J_μ J_ν √-g d⁴x`
- Path E threshold form: `α(ψ) = α_0 (ψ - ψ_th)^n Θ(ψ - ψ_th)` z `ψ_th ≈ 1.05`
- 5/5 POSITIVE structural sketch (cosmology, WEP, GW propagation, stability tree-level)
- Multi-source ISSUE FLAGGED (naive constant α fails by factor 2.3×10⁶)
- Resolution paths A-F analyzed (Path E preferred, Path F hybrid possible)

---

## 3. Brainstorm — kluczowe wymiary problemu

### 3.1 Działania (action) — co rygorystycznie wybiera α(ψ)?

**Pytania:**
- Czy α(ψ) form jest uniquely determined by jakimś principle (uniqueness)?
- Czy threshold ψ_th = 1.05 jest fixed by photon ring geometry, czy free parameter?
- Czy power n w (ψ-ψ_th)^n jest determined (np. n=2 dla canonical kinetic-like)?
- Alternative ansatze: `α(ψ) = α_0 ψ^k`, `α(ψ) = α_0 / (ψ_max - ψ)`, `α(R, ψ)`?
- Czy α może zależeć od Ricci scalar R lub Weyl tensor C?

**Approaches:**
1. **Top-down:** postulate α(ψ) form, verify all constraints
2. **Bottom-up:** derive α(ψ) z TGP substrate microphysics (axiomatic)
3. **EFT:** write most general scalar-tensor action respecting symmetries (Z_2, gauge),
   identify minimal α(ψ) that resolves multi-source
4. **Phenomenological:** parameterize family, scan {α_0, ψ_th, n}, find consistent region

**Recommendation:** start z (1) Top-down + (4) Phenomenological scan, fall back to
(3) EFT jeśli (1) underdetermined. (2) Bottom-up jest desirable ale most expensive.

### 3.2 Modified Φ-EOM (variational derivation)

**Pytania:**
- Jak wygląda δS/δΦ = 0 z α(ψ) coupling?
- Czy wynik jest second-order (Ostrogradsky-free) czy higher-order?
- Czy istnieje ekwiwalentne perfect-fluid representation (auxiliary field trick)?

**Tools needed:**
- Sympy variational derivation (już mamy template z OP-7 work)
- Cross-check w Mathematica / xAct (GR specialized)

**Risk:** high-order EOM → ghost modes → Phase 1 dead-end. Mitigation: choose
α(ψ) form that preserves second-order structure (Horndeski-like discipline).

### 3.3 Photon ring self-consistency (numerical)

**Pytania:**
- Czy modified Φ-EOM produces ψ_ph = 1.168, r_ph = 3.88M w vacuum exterior?
- Czy scenario (e) target +1.46% deviation jest reproducible?
- Czy multi-source check passes (M87*, Sgr A*, NS-NS final)?

**Tools needed:**
- Numerical Φ-EOM integration (scipy ODE)
- Photon geodesic in modified metric
- Newman-Janis or similar machinery for rotating BHs (later)

**Risk:** photon ring shifts away from scenario (e) target → must adjust α(ψ) form,
iterate. **Branch point:** if no α(ψ) reproduces scenario (e) for all M_BH,
Path E fails → must pivot to Candidate A or B (or Path F hybrid).

### 3.4 Multi-source consistency (the critical Phase 0+ flag)

**Pytania:**
- Czy α(ψ) z Path E daje consistent shadow shift across BHs?
- Czy ψ_th = 1.05 jest TOO HIGH (some BHs miss threshold) lub TOO LOW (all activate)?
- Czy lab scale (ψ_lab ≈ 1) faktycznie suppresses α effectively?

**Method:** numerical multi-source scan:
- Sgr A*, M87*, GW150914 final, several stellar-mass BHs from LIGO/Virgo
- Compute predicted shadow deviation for each
- Check universal +14.56% deviation (M9.1'' baseline) vs +1.46% modification (M9.2-D)

**Branch point:** if consistency requires fine-tuning of ψ_th per source → Path E fails.

### 3.5 WEP / 5th force rigorous

**Pytania:**
- Jaki jest rigorous Nordvedt-like η dla MICROSCOPE Ti-Pt z Path E?
- Czy ψ_lab scaling (≈ 1 + 2U_Earth ≈ 1 + 1.4×10⁻⁹) jest below ψ_th = 1.05?
  - Yes! → α effectively zero at lab → η ≈ 0 → WEP MICROSCOPE comfortably PASS
- Composition-dependent contributions od kinetic vs binding energies w internal T^μν?

**Method:**
- Linear perturbation theory in weak field z α(ψ) ≈ 0 (lab)
- Higher-order: small finite α from rare ψ_lab fluctuations
- MICROSCOPE numerical η prediction

**Expected result:** Path E TIGHTENS WEP margin dramatycznie (z 6.7× do astronomical) — to dobre.

### 3.6 Cosmology perturbations (DESI w(z) + CMB)

**Pytania:**
- Czy α(ψ) coupling wpływa na δΦ perturbations w FRW?
- DESI DR3 w(z) prediction: wciąż ≥ −1?
- CMB Planck constraints: czy α(ψ) coupling wprowadza CMB anomalies?
- Sound horizon / H_0 cross-impact?

**Method:**
- Linear perturbation theory in FRW (already partial closure Phase 0+)
- CAMB / CLASS modifications dla M9.2-D + Path E
- Run vs Planck 2018 + DESI DR2

**Risk:** Path E may activate at high-z (cosmological ψ varies), modifying CMB.
Mitigation: ψ_cosmological ≈ 1 ± 10⁻⁵ → still below ψ_th = 1.05 → safe.

### 3.7 Stability (beyond tree level)

**Pytania:**
- 1-loop counterterms: czy α(ψ) coupling generates non-renormalizable divergences?
- Ghost screening dla α > 0 marginal cases?
- Vainshtein-like mechanism activation w strong-field?

**Method:**
- Background field method
- Compare to Horndeski / DHOST literature on screening
- Effective theory cutoff Λ analysis

**Risk:** non-renormalizable → effective theory only, treebly OK. Standard practice.

---

## 4. Dependency graph

```
[1. Action form α(ψ)]
        ↓
[2. Sympy Φ-EOM derivation] ← cross-check Mathematica
        ↓
[3a. Photon ring numerical]    [3b. Linear pert. theory]
        ↓                              ↓
[4a. Multi-source consistency] [4b. WEP MICROSCOPE]
        ↓                              ↓
        +══════════════════════════════+
                       ↓
[5. Cosmology perturbations]
        ↓
[6. Stability beyond tree-level]
        ↓
[7. Paper drafting]
```

**Critical path:** Action form → Φ-EOM → photon ring numerical → multi-source.
Jeśli multi-source fails → loop back to action form.

---

## 5. Milestones (15-month schedule)

### Quarter 1 (Months 1-3): Action + EOM
- **M1.1** Action form α(ψ) determination (top-down + EFT cross-check)
- **M1.2** Sympy variational Φ-EOM derivation
- **M1.3** Cross-check w Mathematica / xAct
- **Deliverable:** rigorous modified Φ-EOM dla M9.2-D + Path E
- **Branch point 1:** if EOM higher-order / unstable → restart with new α(ψ)

### Quarter 2 (Months 4-6): Photon Ring + Multi-Source
- **M2.1** Numerical Φ-EOM integration in Schwarzschild background
- **M2.2** Photon ring location + shadow shift dla Sgr A*
- **M2.3** Multi-source scan: M87*, GW150914, NS-NS finals
- **Deliverable:** scenario (e) +1.46% reproducibility verified universally
- **Branch point 2:** if multi-source fails → consider Path F (D+B hybrid) or
  pivot to Candidate A

### Quarter 3 (Months 7-9): WEP + Cosmology
- **M3.1** WEP rigorous Nordvedt derivation z Path E α(ψ_lab)
- **M3.2** MICROSCOPE Ti-Pt prediction
- **M3.3** Linear cosmology perturbations w/ Path E
- **M3.4** CAMB / CLASS run for M9.2-D vs Planck 2018 + DESI DR2
- **Deliverable:** WEP and cosmology cross-checks at paper rigor
- **Branch point 3:** if MICROSCOPE marginal → adjust ψ_th, re-iterate

### Quarter 4 (Months 10-12): Stability + Refinement
- **M4.1** 1-loop renormalization analysis
- **M4.2** Ghost mode screening
- **M4.3** Effective theory cutoff Λ
- **M4.4** Re-iterate Q1-Q3 z final action structure
- **Deliverable:** complete M9.2-D + Path E theoretical framework
- **Branch point 4:** if non-renormalizable beyond control → restart α(ψ) form

### Quarter 5 (Months 13-15): Paper Drafting
- **M5.1** Paper outline
- **M5.2** First draft (axioms, derivations, predictions)
- **M5.3** Peer feedback (independent reviewers, internal collaborators)
- **M5.4** Final draft
- **Deliverable:** PRD/CQG-ready M9.2 axiom paper

---

## 6. Risk register

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Path E α(ψ) form non-unique | Medium | Medium | EFT classification, scan family |
| EOM higher-order → ghost | Low | High | Choose α(ψ) preserving 2nd-order |
| Multi-source fails after rigorous derivation | Medium | High | Path F hybrid; or pivot to Candidate A |
| MICROSCOPE 6.7× margin shrinks z calibration | Medium | Medium | ψ_lab ≈ 1 → α ≈ 0 saves margin |
| Cosmology perturbations break Planck | Low | High | ψ_cosmological ≈ 1 → α ≈ 0 |
| Non-renormalizable | High | Low | Standard EFT cutoff treatment |
| Phase 1 takes 18+ months | Medium | Low | Prioritize critical path; defer secondary |
| ngEHT 2030+ confirms TGP-like (M9.1'' VINDICATED) | 20% | N/A | Phase 1 work obsolete; pivot to OP-EHT-A++ |
| ngEHT 2030+ delayed to 2034+ | High | Low | Phase 1 absorbs delay (extra time helpful) |

---

## 7. Decision points (branch logic)

**After M1.3 (rigorous EOM):**
- ✅ Second-order EOM, action well-defined → proceed M2
- ❌ Higher-order or ghost → loop back to M1.1 with new α(ψ)
- ❌ Cannot find any α(ψ) that works → Phase 1 EXIT, declare Candidate D
  fundamentally inconsistent, pivot to Candidate A

**After M2.3 (multi-source check):**
- ✅ Scenario (e) reproducible universally → proceed M3
- ⚠️ Reproducible only for some sources → consider Path F hybrid
- ❌ Cannot reproduce → Phase 1 EXIT, pivot to Candidate A

**After M3.4 (cosmology):**
- ✅ Planck + DESI DR2 consistent → proceed M4
- ⚠️ Marginal tension → adjust ψ_th, re-iterate
- ❌ Strong tension → re-examine action structure

**After M4.4 (stability):**
- ✅ Renormalizable / EFT-safe → proceed M5 paper
- ⚠️ Marginal → declare effective theory only
- ❌ Catastrophic instability → restart Phase 1 z new candidate

---

## 8. Resource requirements

**Computational:**
- Symbolic algebra (Sympy + Mathematica/xAct license)
- Numerical (Numpy/Scipy + CAMB/CLASS access)
- Cluster / HPC for cosmology runs (~ 1000 CPU-hours over 15 months)

**Personal:**
- 15-25% effort allocation over 15 months (compatible z OP-DESI 75% during 2026-2030)
- Periodic peer review check-ins (every 3 months)
- 1-2 collaborators z scalar-tensor / modified gravity expertise (recommended)

**Schedule windows:**
- Primary: 2026 Q3-Q4 to 2027 Q4 (5 quarters)
- Secondary: post-ngEHT 2030+ if trigger 2 satisfied (compressed 6-month sprint)

---

## 9. Success criteria

Phase 1 jest zakończona SUKCESEM jeśli:

1. ✅ Rigorous covariant action `S_{M9.2-D + Path E}` shipped
2. ✅ Modified Φ-EOM derived & verified (sympy + Mathematica cross-check)
3. ✅ Photon ring scenario (e) +1.46% reproduced numerically dla Sgr A* AND M87*
4. ✅ Multi-source consistency proven dla 5+ sources spanning 9+ orders of M_BH
5. ✅ MICROSCOPE Ti-Pt prediction comfortably below 1.1×10⁻¹⁵ bound
6. ✅ Planck 2018 CMB consistent (within 2σ)
7. ✅ DESI DR3 prediction within published bounds
8. ✅ Beyond-tree-level stability OK (or clean EFT cutoff identified)
9. ✅ Paper draft ready for submission

Phase 1 jest **PARTIAL SUCCESS** jeśli 6/9 satisfied → publication possible
z explicit caveats.

Phase 1 jest **EXIT-FAILURE** jeśli (1)-(4) niedosięgnięte → Candidate D + Path E
fundamentally inconsistent, must pivot to Candidate A or B as M9.2 path.

---

## 10. Phase 2 (post-Phase 1) preview

Po Phase 1 success:
- **Phase 2.1**: M9.2 axiom paper (PRD/CQG, ~ 6 months drafting)
- **Phase 2.2**: ngEHT collaboration outreach (joint analysis if trigger 2 hits)
- **Phase 2.3**: Followup theoretical work (rotating BHs, GW ringdown, etc.)

---

## 11. Cross-references

- [[research/op-m92/OP_M92_setup.md]] — formal scope (parent)
- [[research/op-m92/OP_M92_decision_tree.md]] — T-M92.6 decision tree
- [[research/op-m92/OP_M92_readiness_summary.md]] — Phase 0 closure
- [[research/op-m92/OP_M92_P0plus_candD_results.md]] — Phase 0+ kickoff
- [[research/op-m92/OP_M92_P0plus_candD_cosmology_results.md]] — cosmology cross-check
- [[research/op-m92/OP_M92_P0plus_candD_wep_results.md]] — WEP cross-check
- [[research/op-m92/OP_M92_P0plus_candD_multisource_results.md]] — multi-source ISSUE
- [[paper/tgp_core.tex]] §applications BH shadows + §7 OP-9 row + F4

---

## Bottom line

Phase 1 jest **9-15 miesięcy rigorous derivation work** zaplanowane do **2026 Q3-Q4
start** (lub post-ngEHT 2030+ accelerated sprint). Plan above provides:

- **Clear deliverables** (9 success criteria)
- **Branch points** (4 critical decisions)
- **Risk register** (9 identified risks z mitigation)
- **Schedule** (5 quarters z monthly milestones)

Phase 0+ readiness package (shipped 2026-04-25) provides **structural scaffolding**:
- Action ansatz pre-derived (sketch)
- Path E α(ψ) preferred resolution identified
- 4 cross-checks at sketch level (cosmology, WEP, GW, multi-source)
- Decision tree pre-prepared dla ngEHT outcomes

⇒ Phase 1 starts NOT from scratch ale z **structural template ready** — response time
post-ngEHT verdict 2-4 weeks (paper drafting + final adjustments) zamiast 2-3 lat.

**Recommended start:** 2026 Q3 (lipiec 2026) — gives 4+ years runway przed ngEHT
2030+ verdict, compatible z OP-DESI/OP-Hubble 75% effort allocation.

**Plan freezing date:** 2026-04-25 (this document). Re-evaluation co kwartał;
substantial revision tylko po ngEHT 2030+ verdict known.
