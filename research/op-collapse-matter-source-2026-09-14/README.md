---
title: "README — op-collapse-matter-source: status kolapsu do granicy dziedziny — czy korpusowe sprzężenie z materią (eq:L-mat-unified) deformuje próżnię Yukawa-podobnie i stabilizuje kolaps?"
date: 2026-09-14
type: cycle-readme
tgp_owner: research/op-collapse-matter-source-2026-09-14
folder_status: closed
status: CLOSED
verdict: "Q-H1-PULL (źródło przy λ̃≥0.2 indukuje ucieczkę z dziedziny; λ̃_crit∈(0.05,0.2]; przy λ̃≤0.05 DEFORMATION Yukawa-podobna spłycona nieliniowo — gate liniowy 10.7%>5%) + Q-H2-INCONCLUSIVE (STABILIZED=0; COLLAPSE=8/12; RADIATED=1 potwierdzone h/2+dt/2 — kolaps baseline'u qR3−0.20 UCHYLONY przy λ̃=0.05, pole osiada w zdeformowanej próżni; INCONCLUSIVE-RUN=3, w tym 2 przetrwańcy t=1000 z E_ref≤0 i 1 rozjazd siatek)"
claim_status: "B"
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[NEEDS.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/NEEDS.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# op-collapse-matter-source (2026-09-14)

Realizacja NEEDS **N2** poprzednika (`op-r3-stationary-states`: kolaps
do granicy dziedziny = reguła 6/10 w dynamice 2. rzędu) ORAZ kandydata
LOCKa **N2-M911** (sprzężenie z materią, osobny LOCK). Materia wchodzi
JEDYNĄ formą z rdzenia: `eq:L-mat-unified` L_mat=−(q/Φ₀)ψρ przez
√−g=c₀ψ/(4−3ψ) ⟹ 𝒰_mat=λ̃ρ̂(r)ψ²/(4−3ψ) (wyprowadzane w P1-H1, nie
zakładane; S05 nienaruszone — ρ to statyczne źródło zewnętrzne).

- **Q-H1 (binarne M911-N2):** źródło na próżni — DEFORMATION
  (Yukawa-podobnie; gate liniowy δψ=−5λ̃·G_Yuk∗ρ̂ przy λ̃=0.01) czy
  THRESHOLD-PULL/COLLAPSE (indukcja przejść przez progi 5/6, 7/6)?
- **Q-H2 (centralne):** czy źródło STABILIZUJE 4 starty-reprezentanty
  kolapsu poprzednika (qR3 ±0.20, gauss −0.30σ3, gauss +0.15σ6)
  przy λ̃∈{0.05,0.2,0.5}?

Pre-rejestrowane (P1-H2/H3): znak deformacji **δψ<0**; asymetria
brzegowa (𝒰_mat→+∞ przy ψ→4/3, →0 przy ψ→0 ⟹ materia odpycha od
sufitu, słabo chroni podłogę). Zakres: odpowiedź na źródło +
stabilizacja; ZAKAZ claimów o masach/oscylonach. Cykl-bliźniak:
[[../op-oscillon-small-amplitude-2026-09-14/README.md]].

Status: **CLOSED** — werdykty wg litery LOCK §4:
[[Phase_FINAL_close.md]]; otwarte pytania: [[NEEDS.md]] (N1–N5).

## Log
- 2026-09-14 — LOCK zapisany (sesja główna; autoryzacja: wybór usera
  „N2: status kolapsu / sprzężenie z materią" po zamknięciu poprzednika
  Q-E-INCONCLUSIVE). Obliczeń zero.
- 2026-09-15 — realizacja pełna (agent-implementator): MD FROZEN →
  Phase 1 sympy (P1-H1/H2/H3 PASS; 𝒰_mat=λ̃ρ̂ψ²/(4−3ψ) WYPROWADZONE,
  δψ_lin=−5λ̃(G_Yuk∗ρ̂), I₀=0.776061935) → Phase 2 PASS 3/3 po
  korekcie 1 (estymator dryfu sekularnego; pierwotny FAIL 2.28e−6 =
  dt²-owy offset ΔC2, po korekcie 8.8e−9; correction note PRZED
  użyciem, pierwotny output zachowany) → Phase 3: Q-H1 6 biegów
  (DEFORMATION przy λ̃≤0.05 z saturacją nieliniową 0.893/0.691
  liniowej; COLLAPSE przy λ̃≥0.2, zbieżnie przy 0.5) i Q-H2 12+9
  biegów (8 COLLAPSE, 1 RADIATED potwierdzone [kolaps uchylony,
  qR3−0.20@0.05], 3 INCONCLUSIVE-RUN) → **Q-H1-PULL +
  Q-H2-INCONCLUSIVE**. Incydenty harnessu (2, bez wpływu na wyniki)
  i probe'y deskryptywne w Phase_FINAL_close §7 i §4–§5. Zamknięcie:
  FINAL + NEEDS + README; integrity_snapshot zweryfikowany.
