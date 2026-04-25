# OP-EHT — Final closure verdict

**Data zamknięcia:** 2026-04-25
**Status:** **CLOSED — CONDITIONAL POSITIVE** (13/18 PASS = 72%)
**Sub-tests:** T1 + T2 + T3 + T4 + T5

---

## TL;DR

> TGP M9.1'' strong-field shadow deviation +14.56% over GR jest **GENUINE
> PHYSICS** (T1 robust to PN truncation n=15). Effective EHT envelope
> NIE absorbuje deviation w obecnym pomiarze Sgr A* (T2 4σ tension; T4
> environment gap ~3σ). Q-renormalization scenario (e) w T3 mathematically
> reduces deviation < 5% z fitted matter coupling f(ψ) = √(g_tt^GR/g_tt^TGP)
> — wymaga first-principles derivation z M9.1'' Lagrangian.
>
> ngEHT 2030+ at 1% precision delivers definitive ≥10σ verdict on Sgr A*
> shadow. M9.1'' pozostaje **falsifiable testable prediction** w window
> 2030-2032.

---

## Summary table

| Test | Cel | Wynik | Score |
|------|-----|-------|-------|
| **T1** | PN truncation robustness | deviation +14.56% genuine, NOT artefact (n=15 stable, ratio 0.269) | 4/4 ROBUST |
| **T2** | Mass calibration MC | M87* OK (59% in 95% CL), Sgr A* 4σ tension | 1/3 |
| **T3** | Q-renormalization audit | G/c² invariant; scenario (e) absorbs to +1.46%; first-principles derivation OPEN | 4/5 |
| **T4** | Non-vacuum environment | M87* env OK (5%); Sgr A* env insufficient (~3σ gap) | 1/3 |
| **T5** | ngEHT window + verdict | ≥5σ discrimination by 2030 (Sgr A* 9.7σ); falsifiable by 2032 | 3/3 |
| **TOTAL** | OP-EHT close | CONDITIONAL POSITIVE | **13/18 = 72%** |

---

## Three resolution paths from EHT_quick §6 — closure

| Path | Test | Outcome |
|------|------|---------|
| 1. Truncation artefact | T1 | **EXCLUDED** — robust to n=15 |
| 2. Calibration absorption | T2 + T4 | **PARTIAL** — M87* OK, Sgr A* fails |
| 3. Q-renorm shift | T3 | **MATHEMATICALLY POSSIBLE** — derivation OPEN |

---

## Strategic implications

### A. M9.1'' status
- **Strong-field deviation +14.56% is genuine prediction** (T1 lock).
- **Falsifiable in 2030-2032 window** (T5 ngEHT 1% → 10σ Sgr A*).
- **Currently 3.3σ tense** with EHT 2022 Sgr A* — within current systematics
  but discriminator emerging.

### B. Q-renormalization scenario (e) — Track A
The fitted ansatz f(ψ) = √(g_tt^GR/g_tt^TGP) absorbs deviation to +1.46%.
For M9.1'' to survive strong-field test, this coupling must emerge naturally
from the action. **Open question:** is this consistent with TGP variational
principle?

If YES → M9.1'' fully reconciled with EHT (3-4σ tension dissolves).
If NO → M9.2 pivot mandatory after ngEHT 2030+ confirms GR.

### C. M9.2 pivot conditional — Track B
M9.2 would need:
- Natural f(ψ) coupling between matter and ψ-field that recovers GR
  in strong-field limit while preserving weak-field PPN.
- Compatibility with OP-7 closure (c_GW = c_0 in scalar/tensor sectors).
- Compatibility with M9.1'' weak-field success (PP P3 audit 3 PASS).

### D. F4 falsifiability hardening
Pre-OP-EHT: F4 was "TGP shadow within EHT 2024+ envelope".
Post-OP-EHT: **F4 narrowed to specific +14.56% deviation, falsifiable
by ngEHT 2030-2032 at >10σ on Sgr A***.

### E. Cross-impact OP-7
OP-7 closure (96.9% PASS, GW c_GW=c_0 unconditional) jest INDEPENDENT
od OP-EHT (T1 σ_ab=0 for static spherical). However:
- OP-7 confirmed σ_ab dynamics cannot rescue static-spherical shadow.
- If M9.2 pivot needed, OP-7 GW result must be re-validated in M9.2.

---

## Path forward — operational

### Immediate (2026):
1. **Update KNOWN_ISSUES.md** — OP-EHT closure entry with full T1-T5 detail.
2. **Update paper §7 OP-EHT row** — CONDITIONAL POSITIVE status, F4
   sharpening.
3. **F4 falsifiability section update** — explicit ngEHT 2030+ window,
   Sgr A* primary target, +14.56% TGP shadow as smoking gun.

### Medium-term (2026-2028):
1. **Track A research:** derive q-renorm scenario (e) coupling
   f(ψ) = √(g_tt^GR/g_tt^TGP) from M9.1'' first principles.
2. **OP-EHT-followup operational program:** if Track A fails, prepare
   M9.2 pivot proposal.
3. **EHT 2024+ analysis:** as soon as new EHT Sgr A* data published,
   re-run T2 with updated priors.

### Long-term (2029-2032):
1. **ngEHT prototype data (2029-2030):** Sgr A* shadow at 1.5% precision
   → 9.7σ discrimination of TGP vs GR.
2. **ngEHT full ops (2030-2032):** Sgr A* + M87* at 1% → 14.6σ verdict.
3. **OP-EHT FINAL closure:**
   - If GR confirmed at 1%: M9.1'' falsified strong-field → M9.2 pivot mandatory.
   - If TGP +14.56% confirmed: M9.1'' decisive triumph.
   - If intermediate (5-10% deviation): partial confirmation, refine theory.

---

## Files

### Setup + scripts
- [[OP_EHT_setup.md]] — formal plan T1-T5
- [[op_eht_t1_pn_truncation.py]] / [[op_eht_t1_pn_truncation.txt]]
- [[op_eht_t2_mass_calibration_mc.py]] / [[op_eht_t2_mass_calibration_mc.txt]]
- [[op_eht_t3_q_renormalization.py]] / [[op_eht_t3_q_renormalization.txt]]
- [[op_eht_t4_environment.py]] / [[op_eht_t4_environment.txt]]
- [[op_eht_t5_ngeht_window.py]] / [[op_eht_t5_ngeht_window.txt]]

### Sub-test verdicts
- [[OP_EHT_T1_results.md]] — 4/4 PASS (deviation ROBUST)
- [[OP_EHT_T2_results.md]] — 1/3 PASS (Sgr A* DISCRIMINATOR)
- [[OP_EHT_T3_results.md]] — 4/5 PASS (q-renorm scenario (e))
- [[OP_EHT_T4_results.md]] — 1/3 PASS (M87* env OK)
- [[OP_EHT_T5_results.md]] — 3/3 PASS (ngEHT 2030+ window)

### Cross-references
- [[research/op7/EHT_quick_results.md]] — punkt wyjścia +14.6%
- [[research/op7/OP7_T1_results.md]] — σ_ab=0 dla static spherical
- [[research/op-newton-momentum/M9_1_pp_P3_results.md]] — M9.1'' P3 audit
- [[paper/tgp_core.tex]] §applications + F4 falsifiability
- [[KNOWN_ISSUES.md]] — OP-EHT closure entry

---

## Bottom line

**OP-EHT CLOSED 2026-04-25 — CONDITIONAL POSITIVE 13/18 PASS (72%).**

TGP M9.1'' nie jest sfalsyfikowane w strong-field, ale **stand-alone defense**
opiera się na otwartej derivation q-renorm scenariusza (e). Realne falsifiable
window 2030-2032 z ngEHT 1% precision na Sgr A* shadow daje >10σ verdict.

M9.1'' status: **alive, testable, falsifiable** — exactly co M9.1'' powinno
być w post-OP-7 era.

M9.2 pivot: **conditional, deferred** — zaadresowany tylko jeśli (a) Track A
q-renorm derivation fails, oraz (b) ngEHT 2030+ confirms GR shadow at 1%.
