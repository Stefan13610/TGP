---
title: "op-omicron2-phi-mean-shift-cosmo — TGP Phi_0(t) tracking — NULL Hubble tension"
date: 2026-05-03
parent: "[[../INDEX.md]]"
type: research_program
status: STAGE_1_NULL_CLOSED_2026-05-03
verdict: |
  STAGE 1 NULL: Stage 0 magnitude (107% Hubble tension match) was structurally wrong.
  Proper D_A integration gives 0.6% (5/5 ICs) — 14× too small. TGP NIE rozwiązuje
  Hubble tension via Phi_0(t) tracking. M10.5 verdict effectively REINFORCED.
key_finding: |
  Pełna verification: matter source term w Phi-EOM strukturalnie istnieje (sek01
  ontology confirmed), psi tracks rho_bar(t) (~20% shift cosmologicznie), ALE
  Lambda_recomb / rho_total_recomb ~ 1e-5 (negligible), więc EDE-like H_recomb boost
  nie triggerуje. Sek05/sek08a internally consistent (NIE sprzeczne); minimal TGP
  coupling jest cosmologicznie tension-insufficient.
trigger: "User intuition 2026-05-03: Phi_0 should track global matter distribution(t)"
tgp_status:
  folder_status: paused
  level: L1
  kind: cosmology_research
  core_compatibility: "POTENTIAL_BREAKING — caveats noted, M10.x results need re-evaluation"
  last_reviewed_against_core: "2026-05-03 (Stage 0 only)"
  may_edit_core: false
  exports_findings: true
  has_needs_file: false
  has_findings_file: false
  open_bridges:
    - "M10.1 w(z) audit needs revision (was -1, now -0.93)"
    - "M10.5 H_0 tension verdict needs revision (was NO solver, now POTENTIAL solver)"
    - "M10.R falsification matrix may need updating"
    - "sek01 ↔ sek08a Phi_0 reconciliation"
  depends_on:
    - "[[../closure_2026-04-26/Lambda_from_Phi0/]] (T-Lambda)"
    - "[[../op-uv3-phi0-renormalization/]] (Z_Phi = 14/3)"
    - "[[../op-gamma1-phi-eff-anchor-resolution/]] (Phi_eff = 8 pi g_tilde)"
    - "[[../op-delta1-g-tilde-derivation/]] (g_tilde formula)"
    - "[[../op-delta2-Nf-derivation/]] (N_f = 5)"
  impacts:
    - "M10.1, M10.4, M10.5, M10.R results need re-evaluation"
    - "sek05 sek06 sek07 may need post-cascade addenda"
    - "DESI compatibility analysis changes"
    - "PUBLICATION blocked until Stage 7 verification"
  promoted_to_core: null
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
tags:
  - TGP
  - omicron2
  - cosmology
  - Hubble-tension
  - Phi-cosmological-evolution
  - sek01-sek08a-reconciliation
  - PRELIMINARY
  - publication-blocked
---

# op-omicron2 — TGP cosmological Φ_0(t) tracking matter via source term

## Cel

Zbadać czy hipoteza wynikająca z sek01 ontologii ("Φ_0 generowane przez kosmologiczną
materię") może rozwiązać Hubble tension przez **cosmological evolution Φ_0(t)** —
mechanism który M10.5 pominął.

## Trigger

User session 2026-05-03 z 5 turami:
1. Initial drift audit — zauważył że M10.5 daje gap 8 dekad do Hubble tension
2. Quick check (Opcja A) — confirmed M10.5 ~10⁻⁸
3. Audit (13 candidates) — kandydat K = matter source term
4. K-test naive — surprise 6 rzędów większy niż M10.5
5. **Z-test precise — 107% Hubble tension coverage z zero free parameters**

## Stage 0 → Stage 1 trajectory

### Stage 0 (initial — INCORRECT magnitude)
```
s_natural = 2.073×10⁻²  (precise z UV.3+γ.1+δ.1+δ.2+K_geo)
|ΔH/H| (Stage 0 formula): 8.93%  ← used `0.5·|dL/L|·Ω_L` [WRONG]
Coverage: 107%  ← would have been hubble tension solver
```

### Stage 1 (verification — NULL result)
```
5 ICs tested z proper V_0 shooting (Ω_DE0(today) = 0.685 enforced)
Direction: ALL EDE-like ✓ (Λ_recomb > Λ_today)
Magnitude: ALL ~0.6% via proper D_A integration
Required: 8.37%
Coverage: 7.3% (NOT 107%)
```

**Stage 1 NULL VERDICT — TGP nie rozwiązuje Hubble tension.**

**Detail:** zob. [[results.md]] (updated post-Stage_1)

## Status — verification cycle CLOSED-NULL

```
✅ Stage 0:  Initial calc DONE 2026-05-03 popołudnie [magnitude WRONG]
✅ Stage 1:  Re-shoot V_0 + D_A integration DONE 2026-05-03 wieczór [NULL]
🔴 Stage 2:  DESI w(z) — ABANDONED (premise nullified)
🔴 Stage 3:  CMB + BBN — ABANDONED
🔴 Stage 4:  Independent replication — ABANDONED
🔴 Stage 5:  Cross-check core — ABANDONED
🔴 Stage 6:  Comprehensive comparison — ABANDONED
🔴 Stage 7:  Snapshot repo + Zenodo — ABANDONED (no claim to publish)
```

**Detail:** zob. [[ROADMAP.md]] (updated post-Stage_1)

## Alternative paths open (post-NULL)

Choć Hubble tension claim nie poszedł, Stage 1 ujawnił:

### Path A: M10.1 w(z) revision (priority HIGH, ~4-6h)
- Matter source ON daje w_today ≈ -0.93 (NIE -1.0 jak M10.1 claimed z de2)
- DESI DR2 -0.75 ± 0.10 — TGP -0.93 may be compatible (closer than LCDM!)
- Worth pursuing as separate cycle

### Path B: Non-minimal coupling speculation (priority LOW)
- Could ax:metric-coupling extension fix tension? Speculative, deferred.

### Path C: Document NULL closure (priority HIGH, ~1h)
- Update sek05/M10.5 with explicit "matter source confirmed; tension-insufficient"
- This README + results.md serve as record

## Pliki

### Source code (Python)
- [[Z_precise_derivation.py]] — main derivation, Stage 0 ✅
- [[Z_output.txt]] — full output Stage 0
- [[omicron2_phi_mean_shift_check.py]] — Opcja A (Buchert variance check)
- [[omicron2_K_test.py]] — Opcja A.K (naive coupling)
- [[quick_test.py]] — fast numerical scan
- [[output.txt]] — Opcja A output
- [[K_test_output.txt]] — K-test output (partial)
- [[quick_test_output.txt]] — quick scan output

### Documentation
- [[README.md]] — this file
- [[results.md]] — Stage 0 results document (PRELIMINARY)
- [[AUDIT_blind_spots.md]] — 13 candidates analysis
- [[ROADMAP.md]] — staged verification + publication plan

## Cross-references — closure cascade input

| Closure | Date | Provides |
|---|---|---|
| T-Λ ([[../closure_2026-04-26/Lambda_from_Phi0/]]) | 2026-04-26 | γ = g̃·M_Pl²·H_0² |
| UV.3 ([[../op-uv3-phi0-renormalization/]]) | 2026-05-02 | Z_Φ = 14/3, Φ_eff vs Φ_0^bare |
| γ.1 ([[../op-gamma1-phi-eff-anchor-resolution/]]) | 2026-05-02 | Φ_eff = 8π·g̃ |
| δ.1 ([[../op-delta1-g-tilde-derivation/]]) | 2026-05-02 | g̃ = N_f·e²/(12π) |
| δ.2 ([[../op-delta2-Nf-derivation/]]) | 2026-05-02 | N_f = 5 derivable |

Bez WSZYSTKICH 5 closures, Z-test result nie jest powtarzalny. Cascade musi być
in-place przed independent replication.

## Cross-references — affected results

Po Z-test, następujące M10 results wymagają re-evaluation:

| Plik | M10 verdict pre-Z | Status post-Z |
|---|---|---|
| [[../op-cosmology-closure/M10_1_results.md]] | w(z) ≈ -1 | 🔴 needs revision (w today = -0.93) |
| [[../op-cosmology-closure/M10_5_results.md]] | TGP nie solver H_0 | 🔴 needs revision (z source term, IS solver candidate) |
| [[../op-cosmology-closure/M10_R_results.md]] | Honest scope: NIE cosmology tensions | 🔴 needs revision |
| [[../desi_dark_energy/README.md]] | TGP indistinguishable from ΛCDM | 🔴 needs revision (TGP quintessence-like) |
| [[../audyt_cosmology_drift_2026-05-03/README.md]] | "M10.5 verdict reinforced" | 🔴 needs **re-statement** (M10.5 verdict actually reversed) |

**Wszystkie te updates są BLOCKED do Stage 4 (independent replication) PASS.**

## Decision points

### Stage 0 → Stage 1
✅ Stage 0 result interesting enough → continue
⬜ User decision pending: which stage to run NOW?

### Po wszystkich Stages 1-6 PASS
→ Stage 7: separate repo `TGP_hubble_tension_2026/` w `C:\Users\Mateusz\Documents\ObsydnianMain\TGP\`
→ Plus Zenodo deposit

### Jeśli któryś stage FAIL
→ Document failure mode, fix, retry
→ Result MAY become "partial Hubble tension coverage" instead of "full solver"
→ Even partial result is valuable (vs M10.5 zero coverage)

## Honest scope post Stage 1 NULL

🔴 **TGP NIE jest Hubble tension solver** — Stage 1 confirms (0.6% vs required 8.4%)
🟢 **User intuition** (sek01 ontology) **strukturalnie POTWIERDZONA** —
   matter source term istnieje, ψ tracks ρ̄(t), sek01-sek08a internally consistent
🟢 **M10.5 verdict efektywnie REINFORCED** — different mechanism (D_A vs Buchert)
   but same conclusion

🟡 **Stage 0 magnitude error** — used wrong formula for tension impact:
   - Stage 0: `dH/H = 0.5·|dL/L|·Ω_L` (LDE-style today's Λ shift)
   - Stage 1: `dH/H_inferred = (1/D_A_ratio - 1)` (EDE-style H_recomb boost)
   - 14× discrepancy between the two — Stage 1 is the right one

🟡 **What WAS revealed (real physics):**
- ψ_today ≈ 0.78 (not 1.0 as M10.1 claimed)
- w_DE_today ≈ -0.93 (not -1.0)
- Lambda_recomb / Lambda_today = 35,319 (huge ratio, but absolute Λ_recomb still
  small vs ρ_m_recomb so doesn't shift H)
- sek01 ontology vs sek08a formalism: NOT contradictory after careful analysis

✅ **NULL result publication-not-needed; documented as honest negative finding**

## Bottom line — post Stage 1

**Pre-cascade TGP (M10.5 era):** TGP NIE rozwiązuje Hubble tension. ✓ confirmed
**Post-cascade TGP (Stage 1 era):** TGP **wciąż** NIE rozwiązuje Hubble tension.

Stage 0 enthusiasm (107% match) okazał się **strukturalnie wrong formula**.
Niemniej cycle wykrył:
1. Matter source term efektywnie jest w sek08a (potwierdza sek01 ontology)
2. ψ rzeczywiście tracks ρ̄(t) cosmologically
3. M10.1 w(z) ≈ -1 może wymagać revision z proper source ON (separate work)
4. M10.5 conclusion (no Hubble tension solver) reinforced przez niezależny mechanism

**User's instinct ("coś tu nie dostrzegam") was right** — Stage 0 was misleading.
**Verification rigor pay-off** — bug caught at Stage 1, NOT after publication.

**Plan publikacyjny ABANDONED.** Path A (w(z) revision) opens as next.
