---
title: "ROADMAP — staged verification (NULL after Stage 1) + alternative paths"
date: 2026-05-03
parent: "[[README.md]]"
type: roadmap
status: stage_1_NULL_stages_2-7_ABANDONED_alternative_paths_open
target: "ABANDONED: Hubble tension solver path. ALTERNATIVE: M10.1 w(z) revision pursue"
---

# ROADMAP — Status post Stage 1 NULL

> **🔴 STAGE 1 NULL VERDICT 2026-05-03 wieczór.**
>
> Stage 0 result (107% Hubble tension match) okazał się strukturalnie błędny
> w magnitude. Stage 1 (proper D_A integration) daje ~0.6% — 14× za mało.
>
> **Hubble tension publication path: ABANDONED** (Stages 2-7 nieistotne).
>
> **Alternative paths open:** M10.1 w(z) revision (separate cycle), non-minimal
> coupling speculation (lower priority).

## Faza 0 — DONE 2026-05-03 popołudnie

✅ Quick check (Opcja A) — confirms M10.5 ~10⁻⁸ (Buchert variance kanał)
✅ Audit (13 candidates) — identifies kandydat K (matter source term)
✅ K-test naive (s=2.5e-3) — surprise: **2.0e-3 ΔH/H** (6 rzędów > M10.5)
✅ Z-test precise (s=2.07e-2 z UV.3+γ.1+δ.1+δ.2+K_geo) — **8.93% match**
   ⚠️ **MAGNITUDE WRONG** — used `dH/H = 0.5·|dL/L|·Ω_L` which is LDE-style today's
       Lambda shift, NOT EDE-style H_recomb boost (Hubble tension's actual mechanism)
✅ results.md saved [[results.md]] (Stage 0 PRELIMINARY, later REVERSED)
✅ ROADMAP.md saved (this doc)

## Faza 1 — DONE 2026-05-03 wieczór — NULL VERDICT

✅ Re-shoot V_0 succeeds for 5 ICs (Ω_DE0(today) = 0.685 ± 0.001 enforced)
✅ Direction CORRECT — wszystkie ICs dają EDE-like (Λ_recomb > Λ_today)
🔴 Magnitude: **dH/H ~ 0.6% via D_A integration** (vs Stage 0 claim 8.9%)
🔴 14× error in Stage 0 — formula was structurally wrong for tension impact
🔴 **All 5 ICs match within 0.1% — robust NULL result**
🔴 Best match: 7.3% of required 8.37%

**Verdict:** TGP NIE rozwiązuje Hubble tension przez Φ_0(t) tracking mechanism.

**Strukturalny powód:** matter source term istnieje (sek01 ontology confirmed),
ALE coupling jest "active" tylko w late universe (z<1) gdzie ρ_m ~ ρ_DE.
Przy z=1100, ρ_m_recomb >> Λ_recomb_TGP (4.4×10⁻⁵ stosunek), więc TGP
nie boostuje H(recomb), nie zmniejsza r_s, nie shiftuje H_0_inferred.

**M10.5 verdict efektywnie REINFORCED przez różne tools (D_A vs Buchert variance).**

## Fazy 2-7 — ABANDONED

```
🔴 Stage 2:  DESI w(z) — ABANDONED (premise nullified by Stage 1)
🔴 Stage 3:  CMB + BBN — ABANDONED
🔴 Stage 4:  Independent replication — ABANDONED
🔴 Stage 5:  Cross-check core — ABANDONED (no Hubble tension claim to check)
🔴 Stage 6:  Comprehensive comparison — ABANDONED
🔴 Stage 7:  Snapshot repo + Zenodo — ABANDONED
```

Hubble tension publication BLOCKED permanently. Plan Stage 7 osobnego repo
+ Zenodo deposit zostaje cancelled.

## ALTERNATIVE PATHS — open after Stage 1

Choć Hubble tension nie wyszedł, Stage 1 ujawnił **realną fizykę** którą warto
zbadać oddzielnie:

### Path A: M10.1 w(z) revision (priority HIGH, ~4-6h)

**Motivation:** Stage 1 ujawnił że w(z) prediction z proper matter source term
nie jest -1, tylko **w_today ≈ -0.93** (zależnie od V_0 i ICs).

DESI DR2 (2025-03) shows w₀ ≈ -0.75 ± 0.10 — TGP -0.93 jest **bliżej DR2 niż
LCDM**. Może być compatible z DESI DR2 evolving DE preference.

**Skrypt:** `pathA_w_z_revision.py`
- FRW solver z matter source ON (s_precise)
- Compute (w₀, w_a) CPL fit z full FRW history
- Compare z DESI DR2 best-fit
- Update M10.1 audit verdict

**Pass criteria:** TGP CPL prediction w 1-2σ od DESI DR2.

### Path B: Non-minimal ρ-Φ coupling speculation (priority LOW, ~separate cycle)

**Speculative:** Could TGP with **non-minimal** ρ-Φ coupling (beyond
ax:metric-coupling) solve Hubble tension?

E.g., adding term `~ ρ²/Φ_0³` w action. This would amplify source effect.

**Risk:** modyfikuje TGP_FOUNDATIONS §4 axiom. **NIE jest minimal TGP.**

**Status:** Park as future research, not immediate.

### Path C: Document and close (priority HIGH, ~1h)

**Cel:** Honest closure of omicron2 cycle as **NULL result**.
- Update sek05/M10.5 z notą "matter source confirmed exists but Hubble-insufficient"
- Mark omicron2 as CLOSED-NULL
- Lessons learned documented

**Recommendation:** Path C + Path A pursue.

## Faza 1 — Re-shoot V_0 (priority HIGH, ~2-3h)

**Cel:** Sprawdzić czy z source term ON, można shoot V_0 by dać Ω_DE0 = 0.685
**dziś** (po 26% Λ-decay), oraz czy ΔH/H pozostaje ~8% z properly shooted V_0.

**Skrypt:** `stage1_reshoot_V0.py`

```python
# Shooting algorithm:
# while not converged:
#     V_0 = V_0_guess
#     run FRW with s_precise + V_0
#     measure rho_DE(today) = (1/2)u² + V_0·V(ψ_today)
#     if rho_DE(today) > Omega_DE0: V_0_guess down
#     else: up
# Once converged, measure |ΔH/H|
```

**Pass criteria:**
- Ω_DE0(today) = 0.685 ± 0.001 osiągnięte przez shooting
- |ΔH/H| ∈ [0.05, 0.12] (50-150% Hubble tension; some flexibility allowed)
- ψ̇(today) > 0 (consistent slow-roll)
- ψ_today ∈ [0.5, 1.0]

**Fail action:** Zaadresować dlaczego shoot V_0 nie działa. Może wymagać większej
elastyczności (np. allow V_0 > 0.685, lub start from bigger ψ_init).

## Faza 2 — DESI w(z) compatibility (priority HIGH, ~3-4h)

**Cel:** Compute (w₀, w_a) CPL fit z full FRW history po Stage 1. Compare
z DESI DR2 actual `(w₀, w_a) = (-0.75 ± 0.10, -0.90 ± 0.35)` przy 3.1σ.

**Skrypt:** `stage2_DESI_compatibility.py`

```python
# Compute w_eff(z) = (KE - PE) / (KE + PE) along FRW trajectory
# Linear lstsq fit w(a) = w_0 + w_a (1-a) on z ∈ [0, 2]
# Mahalanobis distance to DESI DR2 posterior
# Bayes factor TGP vs LCDM vs CPL best-fit
```

**Pass criteria:**
- TGP CPL `(w₀, w_a)` w 2σ od DESI DR2 best-fit (włącznie z covariance)
- TGP fitness lepsza lub porównywalna z LCDM (DR2 favors evolving DE)
- Mahalanobis dist < 2.5 σ

**Acceptable range:**
- w₀ ∈ [-0.95, -0.65]
- w_a ∈ [-1.3, -0.5]

**Fail action:** TGP w(z) niezgodne z DESI → nie jest to tylko Hubble tension solver,
ale też DESI-compatible. Należy sprawdzić, gdzie tracker model się gubi.

## Faza 3 — CMB + BBN check (priority CRITICAL, ~4-6h)

**Cel:** Sprawdzić czy `Λ_eff(recomb) = 1.27× today's value` psuje:
- Planck CMB acoustic peaks (sound horizon, ISW)
- BBN He⁴ abundance

**Skrypty:**
- `stage3a_CMB_check.py` — sound horizon r_d, acoustic peak l_A, Λ-effect
- `stage3b_BBN_check.py` — G_eff drift przy z=10⁹

**Pass criteria CMB:**
- r_d (sound horizon at recomb) shift < 1% vs LCDM
- l_A (acoustic angular scale) shift < 1%
- ISW at z<2 shift < 5×10⁻³
- Planck χ² nie rosnie więcej niż 5 vs LCDM (full likelihood would be bigger task)

**Pass criteria BBN:**
- G_eff(z=10⁹) drift < 1% vs G_Newton
- He⁴ abundance shift w obrębie obs. error ±0.005

**Fail action:** Jeśli CMB ostro psuje, scenario nie zamyka się na pełnym dataset.
Może wymagać modyfikacji (np. matter source becomes active later, after BBN+recomb).

## Faza 4 — Independent replication (priority CRITICAL, ~1-2 dni)

**Cel:** Niezależna derywacja s_natural przez subagenta (clean session),
niezależna FRW integration. Cross-verify s = 2.07e-2.

**Akcja:**
1. Spawn subagent z prompt: "Derive cosmological matter source coupling s_natural
   in TGP from sek08a action + T-Λ closure + Z_Φ = 14/3 (UV.3 closure). Then
   integrate FRW Φ-EOM in matter era and compute |ΔH/H| between recombination
   and today."
2. Compare independent s, ΔH/H z naszym wynikiem.
3. Investigate any discrepancy.

**Pass criteria:**
- Independent s ∈ [1.5e-2, 2.7e-2] (within ~30% of our 2.07e-2)
- Independent ΔH/H ∈ [0.06, 0.12] (within ~30% of our 0.0893)
- Major disagreement → investigate algebra; either we missed factor, or replicator did

**Fail action:** Independent disagreement reveals systematic error in our derivation.
Critical to identify before publication.

## Faza 5 — Cross-check core (priority HIGH, ~1 dzień)

**Cel:** Sprawdzić czy post-cascade scenario (z source term ON) NIE rozwala
istniejących TGP predictions.

**Pliki do sprawdzenia:**
- M10.1 w(z) audit — będzie shifted (already noted in caveats)
- M10.2 inflation — N_e, n_s, r — czy ψ-dynamics during inflation tę samą?
- M10.3 propagator — spatial Yukawa M_eff² = +β — independent of cosmology
- M10.4 CMB safety — covered w Stage 3
- sek02 sektor pole — depends on Φ-EOM stability w obrębie ψ ~ 1
- sek06 czarne dziury — szukać impact: ψ_local przy black hole ≠ ψ_cosmological
- sek07 predictions — większość unaffected (mass spectrum, fundamental constants)
- 40 quantitative predictions w README.md — które shift?

**Skrypt:** `stage5_core_crosscheck.py` (audit-style, list affected vs unaffected)

**Pass criteria:** ≤3 predictions affected non-trivially. Affected ones:
- explicitly listed
- NEW values computed
- consistency z observations checked

**Fail action:** >3 predictions broken → ostrzeżenie: post-cascade scenario
ma znaczące koszty. Decision: czy Hubble tension solution warta tego.

## Faza 6 — Comprehensive comparison (priority MEDIUM, ~3-5 dni)

**Cel:** Pełne porównanie TGP post-cascade vs ΛCDM vs DESI quintom-B vs
CPL vs early-DE vs interacting-DE z full likelihood approach.

**Datasets:**
- Planck 2018 CMB + lensing
- DESI DR2 BAO (arXiv:2503.14738)
- Pantheon+ / DES-SN5YR / Union3 SN
- KiDS-1000 weak lensing
- ACT DR6 small-scale CMB

**Methodology:**
- CosmoMC / Cobaya likelihood pipeline
- TGP Boltzmann solver (modify CLASS or CAMB) — TASK
- MCMC sampling
- Model comparison: ΔAIC, ΔBIC, log-Bayes factor

**Pass criteria:**
- TGP χ² ≤ ΛCDM + 5 (acceptable model, not strongly disfavored)
- ΔAIC vs ΛCDM ≤ 4 (or favored if ΔAIC negative)
- log-Bayes factor TGP/ΛCDM > -1 (no strong evidence against)

**Fail action:** TGP statistically disfavored przez full datasets → not viable
solver, even if Z-test alone matched.

## Faza 7 — Snapshot repo + Zenodo (priority CONTINGENT na Stages 1-6 PASS)

**Cel:** Stworzenie osobnego repository jako snapshot stanu Hubble-tension
research, z pełną reproducibility, plus publikacja Zenodo.

### 7.1 Folder structure (osobne repo)

```
C:\Users\Mateusz\Documents\ObsydnianMain\TGP\TGP_hubble_tension_2026/
├── README.md                  Project overview, abstract, key result
├── LICENSE                    All rights reserved (or CC-BY 4.0)
├── CITATION.cff               Citation metadata
├── .zenodo.json               Zenodo deposit metadata
├── .gitignore
├── derivation/
│   ├── Z_precise_derivation.py     (this script, pinned)
│   ├── stage1_reshoot_V0.py
│   ├── stage2_DESI_compatibility.py
│   ├── stage3a_CMB_check.py
│   ├── stage3b_BBN_check.py
│   ├── stage5_core_crosscheck.py
│   └── outputs/                  (all .txt outputs)
├── papers/
│   ├── tgp_hubble_tension.tex    (PRD-style paper)
│   ├── tgp_hubble_tension.bib
│   └── figures/                  (matplotlib outputs)
├── reproducibility/
│   ├── requirements.lock
│   ├── REPRODUCE.md              Step-by-step instructions
│   └── docker/                   (optional containerization)
└── snapshots/
    ├── results_stage_0.md
    ├── results_stage_1.md
    ├── ...
    └── results_final.md
```

### 7.2 Git workflow

```bash
# Initialize separate repo
cd C:/Users/Mateusz/Documents/ObsydnianMain/TGP
mkdir TGP_hubble_tension_2026
cd TGP_hubble_tension_2026
git init
# Copy artifacts from TGP_v1/research/op-omicron2-phi-mean-shift-cosmo/
# Plus full Stage 1-6 outputs
# Plus paper draft
git add .
git commit -m "Initial snapshot: TGP Hubble tension solver post-cascade Stage 0-6"
git branch -M main
# Optionally: GitHub remote
gh repo create Stefan13610/TGP_hubble_tension_2026 --private
git push -u origin main
```

### 7.3 Zenodo deposit

```yaml
# .zenodo.json
{
  "title": "Topology-Generated Potential resolves Hubble tension via cosmological Φ_0 evolution",
  "description": "Structural derivation of TGP cosmological matter source coupling, giving |ΔH_0/H_0| = 8.93% from zero free parameters; covers Hubble tension at 107%.",
  "creators": [{"name": "Mateusz Serafin"}],
  "keywords": ["TGP", "Hubble tension", "dark energy", "cosmology", "scalar field", "quintessence"],
  "access_right": "open",
  "license": "CC-BY-4.0",
  "communities": [{"identifier": "TGP"}]
}
```

### 7.4 Paper draft outline (PRD-style, ~10 pages)

1. **Abstract** — TGP solves Hubble tension structurally with 0 free parameters
2. **Introduction** — Hubble tension status (DESI Y3, SH0ES vs Planck)
3. **TGP framework** — sek08a action, T-Λ closure, UV.3 + γ.1 + δ.1 + δ.2 cascade
4. **Φ-EOM matter source derivation** — Z_precise step-by-step
5. **FRW integration** — numerical results
6. **Cross-checks** — DESI w(z), CMB, BBN
7. **Comparison** — vs ΛCDM, quintom-B, early-DE
8. **Falsifiability** — testing predictions (DESI DR3, CMB-S4, Euclid)
9. **Conclusions**

### 7.5 Pass criteria publikacji

- All Stages 1-6 PASS
- Independent replication CONFIRMED
- DESI compatibility CONFIRMED
- CMB+BBN safety CONFIRMED
- Core not broken (≤3 affected predictions, all OK observationally)
- Comprehensive comparison: TGP not statistically disfavored
- Paper drafted, internal review done
- All scripts reproducibility-tested

**Only after ALL these → Zenodo deposit + arXiv submission.**

## Failure modes (each stage)

| Stage | Failure mode | Action |
|---|---|---|
| 1 | V_0 shoot impossible / oscillation | Different initial conditions, larger ψ_init range |
| 2 | DESI w(z) > 3σ off | Investigate phase mismatch, modify s coupling form |
| 3a | CMB peaks shift > 1% | TGP needs late-time only mechanism, not full history |
| 3b | BBN G_eff > 1% | Source term must be suppressed at z > 10⁹ |
| 4 | Independent disagreement | Algebra error somewhere, line-by-line check |
| 5 | >3 core predictions broken | Trade-off analysis: H_0 worth it? |
| 6 | TGP statistically disfavored | Full likelihood reveals systematics |
| 7 | Repo / Zenodo blocked | wait until 1-6 fully pass |

## Timeline estimate

| Stage | Time | Cumulative |
|---|---|---|
| 0 | DONE | 2026-05-03 |
| 1 | 2-3h | +2-3h |
| 2 | 3-4h | +5-7h |
| 3 | 4-6h | +9-13h |
| 4 | 1-2 dni | +1-3 dni |
| 5 | 1 dzień | +2-4 dni |
| 6 | 3-5 dni | +5-9 dni |
| 7 | 1-2 dni (po pass) | +6-11 dni |

**Realistic total: 1-2 tygodni intensywnej pracy do Zenodo deposit.**

## Decision points

### Po Stage 1
- ✅ V_0 shoot OK i ΔH/H ~ 8% → **continue**
- ❌ V_0 shoot fails / ΔH/H rośnie do 30% / spada do 1% → **debug stage 0 derivation**

### Po Stage 2
- ✅ DESI compatible → **continue**
- ❌ DESI breaks → **either modify TGP, or accept partial solution**

### Po Stage 3
- ✅ CMB+BBN safe → **continue**
- ❌ CMB+BBN broken → **late-time-only modifications needed (Stage 0 has hidden assumption)**

### Po Stage 4
- ✅ Independent confirms → **continue z confidence**
- ❌ Disagreement → **CRITICAL: identify error before any publication**

### Po Stage 5
- ✅ Core safe (≤3 affected) → **continue**
- ❌ Core broken (>3 affected) → **trade-off review, decide if H_0 solver worth structural cost**

### Po Stage 6
- ✅ Statistical fit OK → **publication green light**
- ❌ TGP disfavored → **paper as "interesting candidate, future work" not "solver"**

### Po Stage 7
- ✅ All checks → **arXiv + Zenodo**
- ❌ Anything fails → **back to relevant stage**

---

## Checklist (mark when DONE)

```
Stage 0:  ✅ DONE 2026-05-03 (preliminary)
Stage 1:  ⬜ pending — re-shoot V_0
Stage 2:  ⬜ pending — DESI w(z) compatibility
Stage 3:  ⬜ pending — CMB + BBN check
Stage 4:  ⬜ pending — independent replication
Stage 5:  ⬜ pending — cross-check core
Stage 6:  ⬜ pending — comprehensive comparison
Stage 7:  ⬜ pending — snapshot repo + Zenodo (CONTINGENT)
```

---

*ROADMAP przygotowany 2026-05-03 post Z-test breakthrough. Każdy stage musi
PASS przed następnym. Publikacja blocked do Stage 7.*
