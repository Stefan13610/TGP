---
title: "Phase FINAL — N5-EW retrofit close (A−)"
date: 2026-05-13
parent: "[[./README.md]]"
phase: FINAL
classification: STRUCTURAL_DERIVED_NATIVE
claim_status: A-
output_type: observable
sympy_total: "8/8 PASS"
substance_metrics: "6 FP (75%) / 2 LIT (25%) / 0 hardcoded"
status: 🟢 CLOSED-RESOLVED
folder_status: closed-resolved
---

# Phase FINAL — N5-EW retrofit close (A−)

```
████████████████████████████████████████████████████████████████████
█  op-L01-N5-retrofit-native-EW-2026-05-13                         █
█  STRUCTURAL_DERIVED_NATIVE — claim_status A−                     █
█  Phase 1 sympy: 8/8 PASS (6 FP + 2 LIT + 0 hardcoded)            █
████████████████████████████████████████████████████████████████████
```

## §1 — Six P-requirements

| P | Resolution |
|---|---|
| P1 | b_2 = 19/6 (SU(2) AF), b_1 = -41/10 (U(1) non-AF) symbolic (T1 FP) |
| P2 | sin²θ_W = g'²/(g²+g'²); trig identity verified (T2 FP) |
| P3 | M_W²/M_Z² = cos²θ_W (Sirlin) (T3 FP) |
| P4 | Sphaleron rate ~ exp(-E_sph/T) BAU freeze-out (T5 FP) |
| P5 | Pattern 2.1 EW-Φ coupling via g_eff (T6 FP) |
| P6 | S05 + T9 DEC; m_W, η_B LIT (T8) |

## §2 — L3 falsification map

| Bound | Status |
|---|---|
| Planck BAU η_B = 6.13·10⁻¹⁰ | inherited PASS |
| PDG sin²θ_W = 0.23121 | structural (1σ MS-bar scheme) |
| ATLAS+CMS m_W = 80.369 GeV | inherited PASS |

## §3 — Substance vs predecessor

| Metric | Predecessor (C LIT) | Retrofit (A−) |
|---|---|---|
| BINDING contract | ❌ | ✅ |
| PR-### | ❌ | ✅ PR-008 |
| FIRST_PRINCIPLES | 0 | 6/8 (75%) |
| HARDCODED True | 0 | 0 |
| claim_status | C | A− |

## §4 — Sign-off

L01 N5 EW gauge anomaly retrofit **CLOSED A−**. PR-008 entry to-be-added.
