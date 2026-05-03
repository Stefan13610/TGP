---
title: "ω.2 mini-audit — E_TGP cascade verification + g_anomaly status"
date: 2026-05-04
cycle: ω.2
type: post-hoc-mini-audit
status: CLOSED — NO BLOCKING CRITIQUE
parent: "[[Phase3_results.md]]"
sibling-audits:
  - "[[../op-chi1-newton-constant-derivation/CRITIQUE_circular_anchor_2026-05-02.md]]"
  - "[[../op-uv2-mtgp-absolute-scale/CRITIQUE_repackaged_circularity_2026-05-02.md]]"
  - "[[AUDIT_omega3_2026-05-04.md]]" (sibling cascade)
binding:
  - "[[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] §5.1.1"
  - "[[../../meta/AUDYT_TGP_2026-05-01.md]] §A.7 (already-LIVE-PARTIAL)"
tags:
  - TGP
  - omega2
  - mini-audit
  - cascade-verification
  - anti-overclaim
  - already-LIVE-PARTIAL
---

# ω.2 mini-audit — E_TGP cascade verification

**Verdict:** ω.2 entry **NIE wymaga BLOCKING critique** analogicznej do χ.1/UV.2.
E_TGP = 536/75 jest mechanicznie wyprowadzone z pre-existing θ.1 anchors
(B²_up = 13/4, B²_down = 61/25, B²_lep = 2) przez standardową triangle-anomaly
formułę QED. Status `LIVE PARTIAL` (Δχ² = 0.21 vs α_em alone) **był już
poprawnie ustanowiony** w [[../../meta/AUDYT_TGP_2026-05-01.md]] §A.7.
Niniejszy mini-audit potwierdza tę decyzję; nie wprowadza nowych downgrade'ów.

---

## 1. Pytania §5.1.1 z subagent_audit

### 1.1 Czy E_TGP = 536/75 jest mechanicznie z θ.1/ρ.1/η.2 czy fitted internally?

**Werdykt: MECHANICZNIE z θ.1.** [[Phase1_results.md]] lin. 40-49:

```
B²_up = 13/4    (z θ.1 quark Koide)
B²_down = 61/25 (z θ.1 derived: B²_up − 81/100)
B²_lep = 2      (Dirac leptons, standard)
N_c = 3, Q_u = 2/3, Q_d = -1/3, Q_l = -1   (SM gauge charges)

E_TGP = N_c·[Q_u²·B²_up + Q_d²·B²_down] + 1·Q_l²·B²_lep
      = 3·[4/9·13/4 + 1/9·61/25] + 1·1·2
      = 3·[13/9 + 61/225] + 2
      = 3·386/225 + 2
      = 386/75 + 150/75
      = 536/75 ≈ 7.1467
```

**To jest standardowy QED triangle anomaly z TGP B² weighting** —
wszystkie inputs są pre-existing structural (θ.1 + Dirac). **NIE
post-hoc fitted.** Sympy LOCK rzeczywiście jest algebraiczną tożsamością.

### 1.2 Czy g_anomaly = α_em·E_TGP/(2π) ≈ 8.30·10⁻³ ma derivation czy pasuje do CMB birefringence post-hoc?

**Werdykt: STANDARDOWA forma EM anomaly one-loop.** Phase 1 lin. 63:

```
g_eff(IR) = α_em · E_TGP / (2π) = (1/137.036)·(536/75)/(2π) ≈ 8.300·10⁻³
```

To jest **standardowa formuła axion-photon coupling** z one-loop EM anomaly:
`g_aγ = α_em · (E/N) / (2π·f_a)` w zwykłej notacji KSVZ/DFSZ. TGP nie
"fituje" — wstawia własną wartość E/N (czyli E_TGP, bo brak QCD anomaly N).

**Phase 2 lin. 32 jawnie raportuje:** "Combined χ² ranking FAIL — Δχ²=0.21
vs α_em (degenerate at CMB current precision)". To jest **honest
nie-rozstrzygający wynik**, nie over-claim. Phase 3 status PARTIAL DERIVED
+ forward-gate SO 2027+/LiteBIRD 2029+ jest już zgodny z tą diagnozą.

### 1.3 Czy 4 alternative cross-channel candidates (κ_TGP, 1/(2π), η_chir) faktycznie wykluczane > 10σ czy fitted?

**Werdykt: alt-scan strukturalnie poprawny, ale dyskryminacja słaba.**
Phase 1 lin. 67-72 (cosmological feasibility Δ(lnX)):

| Candidate | g | Δ(lnX)_required | feasibility |
|---|---|---|---|
| α_em alone | 7.30·10⁻³ | 1.63 | OK (closest) |
| α_em·E_TGP/(2π) | 8.30·10⁻³ | 1.43 | OK |
| η_chir = 19/24 | 0.7917 | 0.0149 | borderline (very small) |
| κ_TGP | 2.012 | 0.0058 | borderline |

3/4 falsified by extreme Δ(lnX) requirements; α_em vs g_anomaly remain
**degenerate**. **NIE 10σ exclusion** — tylko ranking-by-feasibility.
g_anomaly winner przez najmniejszy χ² + structural inheritance, ale
gap < 3σ.

---

## 2. Co ω.2 naprawdę locks vs co claim

| Claim w prose-bloku 2026-05-01 | Reality post-mini-audit |
|---|---|
| "ω.2 g_axion = α_em·E_TGP/(2π) STRUCTURAL WINNER" | **Bayesian-favored, NOT 3σ unique** (Δχ²=0.21) |
| "3/4 candidates FALSIFIED >10σ" | **3/4 falsified by feasibility, not >10σ** |
| "FULL CONVERGENCE 18/18" | **18/18 sub-tests PASS, ale Phase 2 W2.2.6 explicit FAIL ranking 3σ** |
| "g_anomaly LOCKED" | **PARTIAL DERIVED — forward gate SO 2027+ 5σ pending** |
| "ω.2 program END" | **PROGRAM ENDED, lock CONDITIONAL** |

### 2.1 Co jest faktycznie strukturalnie poprawne (zachować)

- E_TGP = 536/75 sympy LOCK z θ.1/Dirac inheritance (Phase 1 W2.1.1) ✓
- η_chir = 19/24 dual-derivation Form A/B (Phase 1 W2.1.2) ✓
- 4-alt cross-channel falsification matrix (Phase 1 W2.1.5) ✓
  (struktura tak, magnitude rankingu warunkowy)
- Phase 2 χ² ranking honest — explicit Δχ²=0.21 < 9 reported ✓
- Forward-gate SO 2027+/LiteBIRD 2029+ realistic ✓

### 2.2 Co jest over-claimed (downgrade)

- "STRUCTURAL WINNER" sugeruje >3σ — w rzeczywistości 1.6σ now
- "FULL CONVERGENCE" mylący gdy W2.2.6 explicit FAIL z gap=0.21
- "LOCKED" → **PARTIAL DERIVED** (Phase 3 sam to potwierdza)

---

## 3. Status post-mini-audit

| audit item | werdict | rationale |
|---|---|---|
| ω.2 strukturalne sympy LOCK (E_TGP=536/75) | **CONFIRMED** | mechanicznie z θ.1 |
| ω.2 g_anomaly status | **LIVE PARTIAL** (already in §J of audit) | Δχ²=0.21 < 9 |
| ω.2 cascade-conditional na θ.1 | **OK** | θ.1 itself stably LOCKED |
| ω.2 nowe BLOCKING critique | **BRAK** | nie analogiczne do χ.1/UV.2 |

**Decyzja:** rejestrowy status `LIVE PARTIAL` z PREDICTIONS_REGISTRY pozostaje;
Phase3_results.md verdict `PARTIAL DERIVED` pozostaje. Tylko prose-blok
"Updated 2026-05-01" w registry jest oznaczany jako "PENDING audit"
w `REVISION 2026-05-04` (już zrobione).

---

## 4. WW7-WW12 per-row epistemic status

| entry | original 2026-05-01 | post-mini-audit 2026-05-04 |
|---|---|---|
| WW7 (E_TGP = 536/75) | LOCKED structural | **LOCKED-ALGEBRAIC** (mechanically from θ.1) |
| WW8 (g_anomaly cosmological) | LOCKED-CONDITIONAL | **LIVE PARTIAL** (already §J) |
| WW9 (η_chir dual-derivation) | LOCKED structural | **LOCKED-ALGEBRAIC** (Phase 1 W2.1.2 dual) |
| WW10 (g_anomaly Bayesian winner) | LOCKED | **LIVE PARTIAL** (Δχ²=0.21 honest reported) |
| WW11 (CMB β = 0.342° band) | LIVE PARTIAL | **LIVE PARTIAL** (already §B11 caveat) |
| WW12 (forward-gate SO 2027+) | LOCKED forecast | **STRUCTURAL** forecast (5σ post-2027) |

---

## 5. Wnioski

ω.2 jest **honest**, mimo że ledger prose 2026-05-01 używał optymistycznego
języka ("FULL CONVERGENCE", "STRUCTURAL WINNER"). Audit AUDYT 2026-05-01 §A.7
już to dostrzegł i downgrade'ował status do LIVE PARTIAL. Niniejszy
mini-audit:

1. **Confirmuje** że E_TGP=536/75 jest mechaniczna, NIE fitted.
2. **Confirmuje** że Phase 2 W2.2.6 explicit FAIL Δχ²=0.21 jest reported transparently.
3. **NIE wprowadza** nowego BLOCKING critique.
4. **Synchronizuje** ledger prose (PREDICTIONS_REGISTRY §"REVISION 2026-05-04")
   z już-existing §A.7 downgrade.

**Kontrast z χ.1 i UV.2:** chi.1 ma algebraiczną tautologię
`G_N · M_Pl² = 1`; UV.2 ma post-hoc K_struct fitting w M_GUT band 10-30%.
ω.2 ma uczciwą cascade z θ.1 + standardową EM anomaly formułę. **NIE jest
podobna patologia.**

---

**Autor:** ω.2 mini-audit per [[../../meta/SUBAGENT_AUDIT_74394a8_2026-05-02.md]] §5.1.1.
**Data:** 2026-05-04.
**Status:** CLOSED — no new critique. Confirms LIVE PARTIAL status (already audited 2026-05-01).
**Output:** WW7-WW12 per-row epistemic status table → [[../../PREDICTIONS_REGISTRY.md]] §"REVISION 2026-05-04".
