---
title: "CORE_INTEGRATION_PLAN — M03 retrofit annotations w core/ (Opcja C+ extension)"
date: 2026-05-06
parent: "[[POST_ACTION_FINAL_2026-05-06.md]]"
type: integration-plan
status: AWAITING_USER_CONFIRMATION
tgp_owner: audyt/M03_balance_sheet_missing
tags:
  - integration-plan
  - core-edits
  - opcja-c-extension
  - awaiting-confirmation
related:
  - "[[POST_ACTION_FINAL_2026-05-06.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../../research/op-M03-balance-sheet-retrofit-2026-05-06/audit_log.md]]"
---

# CORE_INTEGRATION_PLAN — M03 retrofit core annotations

## Status

**AWAITING USER CONFIRMATION** przed faktycznymi edycjami core LaTeX.

> Per AGENT_PROTOCOL §0 + M03 framework: nie modyfikujemy core/ bez
> explicit OK. Niniejszy dokument identyfikuje proposed edits.

## 1. Scan core/ — current state

### 1.1 Cycle naming references w core/

Po `grep`:
- `core/sek00_summary/sek00_summary.tex` — JEDYNY plik core używający
  cycle naming z M03 retrofit (κ.1, ι.1, μ.1, ζ.1)
- `core/_meta_latex/status_map.tex` — używa abstract framework
  (Aksjomat/Twierdzenie/Propozycja/Hipoteza), NIE cycle naming
- Brak referencji do UV.2 K_struct, XS.1, ν.1, Majorana phase, F4
  chain w sek00 ani innych core files

**Wniosek:** Core LaTeX **głównie operuje na meta-level structural
axioms** (M9, sek08, sek09, dodatekU itp.). Cycle naming używane TYLKO w
master tabeli sek00:320-348 jako reference do research-level cykli.

### 1.2 Annotations już dodane (sesja A 2026-05-06, Opcja C)

✅ **Już w core/sek00_summary/sek00_summary.tex:**

| Linia | Cycle | Annotation |
|-------|-------|------------|
| 323 | ζ.1 PMNS | "M03 retrofit 2026-05-06: STRUCTURAL — ι.1+μ.1 promotions DERIVED *withdrawn*" |
| 326 | μ.1 δ_PMNS Form A | "M03 retrofit 2026-05-06: NUMEROLOGICAL --- multi-form ambiguity z Form B 260°, 55° gap; 'fitted corrections' via κ.1 cascade" |
| 331 | Σm_ν Z1 | "B4-locked Z1 anchor 2026-05-01, D01 2026-05-06" |
| 341 | m_H | "D01-locked 2026-05-06" |

**Sesja A już pokryła główne mixing-operator family + D01 anchor M03
findings w core.**

## 2. Gaps post-sesja-J (proposed Opcja C+ extension)

### 2.1 Minimal extension (recommended, ~3-5 line edits)

Z M03 100% findings (sesje B-J) potencjalne extensions:

**E1 (sek00:331 Σm_ν line) — DODAĆ ο.1 7-channel forecast:**
```diff
- ($\Sigma=59.01$, B4-locked Z1 anchor 2026-05-01, D01 2026-05-06; pre-lock 59.6); IO wykluczone
+ ($\Sigma=59.01$, B4-locked Z1 anchor 2026-05-01, D01 2026-05-06; pre-lock 59.6;
+ \textbf{$\omicron$.1 Form A LOCKED}: 7-channel falsification, CMB-S4 2030+ 5.78$\sigma$ DECISIVE; M03 retrofit 2026-05-06: STRUCTURAL); IO wykluczone
```

**E2 (sek00:323 PMNS line) — DODAĆ ν.1 Majorana annotation:**
```diff
- M03 retrofit 2026-05-06: STRUCTURAL — $\iota$.1+$\mu$.1 promotions DERIVED \emph{withdrawn})
+ M03 retrofit 2026-05-06: STRUCTURAL — $\iota$.1+$\mu$.1 promotions DERIVED \emph{withdrawn};
+ $\nu$.1 Majorana phases: Form A $\alpha_{21}{=}\pi/2$ STRUCTURAL clean (chirality halving) + Form B $\alpha_{21}{=}11\pi/13$ NUMEROLOGICAL ($\mu$.1 cascade))
```

**E3 (sek00:348 Master Podsumowanie v2 line) — DODAĆ M03 closure annotation:**
```diff
- Podsumowanie v2: 35 skrypt\'ow, 307/348 (88.2\%); 12 $\star\star\star$; 40 predykcji; 0/15 kill
+ Podsumowanie v2: 35 skrypt\'ow, 307/348 (88.2\%); 12 $\star\star\star$; 40 predykcji; 0/15 kill;
+ \textbf{M03 retrofit 2026-05-06}: 40 cykli audited (33 $\star$ honest = 82.5\%), counter 856 $\to$ $\sim$712 effective uncontested (estimate)
```

**E4 (status_map.tex append section) — DODAĆ M03 closure reference:**
Add at end (after line ~1690):
```latex
\subsection*{M03 balance-sheet retrofit closure (2026-05-06)}

Per \texttt{audyt/M03\_balance\_sheet\_missing/POST\_ACTION\_FINAL\_2026-05-06.md}:
\begin{itemize}
  \item Phase 1-4: 100\% COMPLETE (45/45 effective scope audited)
  \item Honest reporting baseline: 33/40 = 82.5\%
  \item Mixing-operator family ISOLATED: $\kappa$.1+$\iota$.1+$\mu$.1+$\nu$.1 (NUMEROLOGICAL/ANSATZ)
  \item Counter recalc: 856 $\to$ $\sim$712 effective uncontested (estimate; full Phase 5 implementation deferred)
  \item Phase 6 ABSOLUTE BINDING gate: \texttt{meta/CALIBRATION\_PROTOCOL.md} + \texttt{CALIBRATION\_GATE\_ENFORCEMENT.md}
\end{itemize}
```

**Total minimal extension:** 4 edits w 2 plikach (sek00 + status_map).

### 2.2 Larger extensions (deferred, NIE recommended w 1 sesji)

| # | Extension | Scope | Comment |
|---|-----------|-------|---------|
| L1 | sek09 PMNS sektor explicit M03 annotations | sek09:1428+ multiple lines | Wymaga careful integration z mass formulas |
| L2 | dodatekU (UV completion) DEPRECATES UV.2 K_struct | dodatekU 5+ lines | UV.2 explicit referenced w UV taxonomy; potrzebuje DEPRECATION block |
| L3 | dodatekQ (coarse graining) M03 retrofit references | dodatekQ Q.4 a_Γ identification | Już corrected przez UV.3 retrofit (Φ_eff=24.65), ale NIE explicit referenced |
| L4 | tgp_companion.tex M03 closure summary | companion full integration | ~10+ lines, structural rewrite section |
| L5 | sek08a vacuum stability annotations | sek08a multiple sections | Phase 1 covariant + Phase 2 EFT closure references |

**Recommendation:** L1-L5 są **deferred** — wymagałyby osobnej sesji
2-3h dla careful structural integration. Phase 5 full PREDICTIONS_REGISTRY
(per-row epistemic class tagging) byłby naturalniejszym następnym
krokiem niż core LaTeX deep integration.

## 3. Decision matrix

| Opcja | Scope | Time | Risk |
|-------|-------|------|------|
| **Opcja A (minimal):** E1+E2+E3+E4 | 4 edits, 2 pliki | ~30 min | Low (annotations only, NIE structural) |
| **Opcja B (medium):** + L1 (sek09 PMNS) | +5 edits, sek09 | ~1.5h | Medium (PMNS sektor delicate) |
| **Opcja C (large):** + L1-L5 | +20 edits, 5 plików | ~3-4h | High (structural rewrite, careful review needed) |
| **Opcja D (skip):** No core edits | 0 edits | 0 | None — Phase 5 PREDICTIONS_REGISTRY priority instead |

## 4. Recommended action

**Opcja A (minimal)** — 4 line edits w sek00 + status_map dla M03 closure
post-sesja-J reflection. Niska risk, wysoka clarity.

Po Opcja A → user może ocenić czy continue z Opcja B/C lub przejść do
Phase 5 full implementation.

## 5. Awaiting user confirmation

Zatwierdzenie potrzebne przed faktycznymi edycjami core/. Format
zgody:
- "OK Opcja A" — execute E1+E2+E3+E4 (4 edits)
- "OK Opcja B" — execute Opcja A + L1 (sek09 PMNS)
- "OK Opcja C" — execute wszystko
- "Skip core, focus Phase 5" — Opcja D, fokus na PREDICTIONS_REGISTRY refactor
- Ad-hoc selection (wybierz E1/E2/E3/E4 osobno)

## 6. Cross-references

- [[POST_ACTION_FINAL_2026-05-06.md]] — M03 100% closure synthesis
- [[../../research/op-M03-balance-sheet-retrofit-2026-05-06/audit_log.md]] — 10-session log
- [[../../research/op-M03-balance-sheet-retrofit-2026-05-06/Phase5_registry_refactor_draft.md]] — Phase 5 alternative path
- [[../../core/sek00_summary/sek00_summary.tex]] — primary annotation target (lines 320-348)
- [[../../core/_meta_latex/status_map.tex]] — secondary annotation target (line ~1690)
- [[../../meta/CALIBRATION_PROTOCOL.md]] — ABSOLUTE BINDING reference
