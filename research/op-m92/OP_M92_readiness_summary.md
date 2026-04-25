# OP-M92 — Readiness summary (interim closure)

**Data:** 2026-04-25
**Status:** **OPEN — PHASE 0 COMPLETE + PHASE 0+ KICKOFF DONE**
**Sub-tests:** T-M92.1 + T-M92.2..T-M92.5 + T-M92.6 + Phase 0+ Candidate D sketch done
**Deferred closure:** ngEHT 2030+ verdict

**Phase 0+ update (2026-04-25):** Candidate D (momentum back-reaction)
structural sketch shipped (5/5 POSITIVE). Single coupling α ~ 0.1 (geom
units) simultaneously tunes scenario (e) target i preserves weak-field PPN
z 9e+12× safety margin. No-ghost stable, c_GW = c_0 vacuum exact.
See [[research/op-m92/OP_M92_P0plus_candD_results.md]].

**Phase 0+ cosmology cross-check (2026-04-25 popołudnie):** Candidate D
preserves OP-DESI w(z) ≥ −1. Ratio 8πG α ρ_m ~ 8e-34 today, ~5e-8 at BBN,
phantom transition pushed to z ≈ 1e+11 (pre-BBN, unobservable régime).
de2 theorem strukturalnie zachowane.
See [[research/op-m92/OP_M92_P0plus_candD_cosmology_results.md]].

**Phase 0+ WEP cross-check (2026-04-25 wieczór):** Candidate D modyfikuje
TYLKO sektor grawitacyjny (matter EOM niezmodyfikowane) → UFF preserved
strukturalnie. Residual Nordvedt-like η_MICROSCOPE ~ 1.6e-16 vs bound
1.1e-15 (margin 6.7× TIGHT); η_LLR ~ 4.7e-29 vs 1.4e-13 (margin 3e+15×).
MICROSCOPE = tightest constraint, wymaga rigorous Phase 1 derivation.
See [[research/op-m92/OP_M92_P0plus_candD_wep_results.md]].

---

## TL;DR

> OP-M92 jest **deferred conditional research program**. Phase 0 (readiness
> package) closed 2026-04-25 z 4 pre-analyzed kandydatami M9.2 i decision
> tree gotowym do execution post-ngEHT 2030+.
>
> **Phase 0 deliverable:** zamiast 2-3 lat from-scratch po verdict,
> M9.2 paper response time 2-4 weeks po ngEHT decisive measurement.
>
> **Recommendation ranking (pre-verdict):**
> 1. Candidate D (momentum back-reaction) — PROMISING
> 2. Candidate A (dual-field) — VIABLE z screening
> 3. Candidate B (conformal frame) — VIABLE constrained
> 4. Candidate C (q-flow) — NOT VIABLE bez ψ-threshold
>
> M9.1'' standalone status: **alive, falsifiable in 4-6 lat.**

---

## Files (Phase 0)

- `OP_M92_setup.md` — formal scope (trigger conditions, candidates, plan)
- `op_m92_T1_action_reverse_engineer.py` + `.txt` — T-M92.1 diagnostic
- `op_m92_T2_T5_candidate_analysis.py` + `.txt` — T-M92.2..T-M92.5
- `OP_M92_decision_tree.md` — T-M92.6 decision tree synthesis
- `OP_M92_readiness_summary.md` — this file (interim closure)

## Key findings (Phase 0)

### Finding 1: Scenario (e) coupling is doubly non-local

T-M92.1 confirmed że f(ψ) = √(g_tt^GR/g_tt^TGP) wymaga porównania:
- TGP photon ring (ψ_ph=1.168, r_ph^TGP=3.88M)
- GR photon ring (r_ph^GR=3M)

To jest **doubly non-local** — różne r values dla dwóch metryk —
fundamentalnie niemożliwe do reprodukcji w single-field M9.1''.

### Finding 2: Naive 1PN expansion of scenario (e) breaks γ_PPN

Linear-U coefficient w f_e(ψ=1+2U) = 1 + 3U + ... → γ_PPN shift O(1)
naruszający Cassini bound 2.3e-5. Scenario (e) jest fitted ansatz tylko
przy photon ring, nie globally consistent matter coupling.

⇒ Każdy M9.2 kandydat musi **strukturalnie** odróżnić strong-field od
weak-field (np. ψ-threshold, T^μν-trigger, conformal screening).

### Finding 3: Candidate D (momentum back-reaction) leads ranking

Power: T·J·J term naturally O(U^4) w weak-field → γ_PPN, β_PPN auto-PASS
bez fine-tuning. W strong-field T^μν dominuje → activation. Single candidate
that satisfies (a)+(b) bez explicit screening.

### Finding 4: q-flow candidate eliminated

Power-law q(ψ) ~ ψ^(-3/4) matches photon ring (0.4% gap) ale 1PN linear-U
coefficient -3/2 dramatycznie naruszający Mercury γ_PPN. Bez ψ-threshold
mechanism, nie ma viable q-flow path.

## Status checklist

- [x] OP-EHT-A NEGATIVE → trigger 1 satisfied
- [x] Scope document (OP_M92_setup.md)
- [x] Diagnostic T1 (action principle requirements)
- [x] Candidate analysis T2-T5 (4 candidates ranked)
- [x] Decision tree T6 (ngEHT response plan)
- [x] Readiness summary (this file)
- [x] **Phase 0+ Candidate D sketch (2026-04-25 afternoon)** — 5/5 POSITIVE
- [x] **Phase 0+ Candidate D cosmology cross-check (2026-04-25)** — POSITIVE (w(z) ≥ −1 preserved)
- [x] **Phase 0+ Candidate D WEP cross-check (2026-04-25)** — POSITIVE (UFF structural; MICROSCOPE margin 6.7× TIGHT)
- [ ] ngEHT 2030+ verdict (DEFERRED)
- [ ] Candidate D Phase 1 derivation (scheduled 2026 Q3-Q4)
- [ ] Full M9.2 axiom paper (DEFERRED to 2032+)

## Bottom line

OP-M92 Phase 0 **complete**. Pre-emptive readiness reduces post-verdict
response time z 6-12 months do 2-4 weeks. M9.1'' status post-OP-M92:

- **Theory health:** HIGH (active baseline, falsifiable, operationally prepared)
- **Strong-field**: deferred conditional (ngEHT 2030+)
- **Weak-field**: PASS (Mercury, Cassini, LLR)
- **GW propagation**: PASS (OP-7 c_GW=c_0)
- **Cosmology**: separate independent scope (OP-DESI, OP-Hubble)

**Next operational priority** post-OP-M92 Phase 0:
- Continue OP-DESI / OP-Hubble cosmology work (75% effort)
- Phase 0+ Candidate D derivation sketch (20% effort, 2026 Q3-Q4)
- EHT incremental data monitoring (5% effort, 2026-2030)

OP-M92 stays OPEN with no immediate work blocked on it. **Phase 0 ships**
za 1 dzień (this session); next milestone Phase 1 trigger w 2030+.
