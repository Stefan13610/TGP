---
title: "Phase 1 — werdykt RESCUE-CONFIRMED (8/8 FP, wyliczony): niezależny solver odtworzył dodatekX (m_b 0,59%, m_t 0,77% pole); HALT-B testował strawmana (3/3 składników brak). Sektor kwarkowy = częściowo predyktywny (r_31 z r_21+𝒜), NIE structural insufficiency."
type: phase1_results
status: COMPLETE
date: 2026-06-25
cycle: op-quark-mass-core-g0-rescue-test-2026-06-25
parent: "[[./README.md]]"
sympy_pass: "8/8 FP"
verdict: RESCUE-CONFIRMED
hardcoded: 0
---

# Phase 1 — wyniki

> **Werdykt (WYLICZONY z flag T1–T5 per reguła §0.2 LOCKED): RESCUE-CONFIRMED.**
> Niezależny solver samozgodnego domknięcia `{m_0=𝒜·m_3/m_1 ; Q_K(m_i+m_0)=2/3}`
> **odtworzył** predykcję dodatekX: m_b 0,59%, m_t 0,77% (pole). Formuła HALT-B
> nie zawiera **żadnego** z 3 składników maszynerii TGP ⇒ **strawman**. Sektor
> kwarkowy jest **częściowo predyktywny** (r_31 predykowane z r_21 + uniwersalne
> 𝒜=a_Γ/φ=C_F²α_s²), a **nie** „structural insufficiency 0/5".

## §1 — Niezależna re-weryfikacja (T2, CENTRALNY)

Solver (sympy `nsolve`, zero hardcoded) rozwiązał dla każdego sektora samozgodnie
`m_0=𝒜·m_3/m_1` ∧ `Q_K(m_1+m_0,m_2+m_0,m_3+m_0)=2/3`, 𝒜=a_Γ/φ=0,02472:

| Sektor | wejście (m_1,m_2) | m_0^pred | **m_3^pred** | PDG | δ |
|---|---|---|---|---|---|
| down | (4,67; 93,4) | 22,3 MeV | **m_b = 4205 MeV** | 4180 | **0,59%** |
| up | (2,16; 1270) | 1962 MeV | **m_t = 171 435 MeV** | 172760 (pole) | **0,77%** |
| | | | | 162500 (MS-bar) | 5,50% |

Zgodne z dodatekX (m_b 0,24%, m_t 2,5%) w granicach metody. **m_0 odtworzone** (22,3≈21,9;
1962≈1981). **Predykcja m_t bliższa schematowi pole (0,77%) niż MS-bar (5,50%).**

## §2 — Pozostałe testy

- **T1 (PASS):** φ-FP anchors {down 0,817; lepton 0,870; up 0,891} = sek08b:529 [0,817;0,891]
  → **dosłownie potwierdza #40 NORM-OVERLOAD** (przedział to kotwice sektorów, nie domena).
- **T3 (PASS):** 𝒜=a_Γ/φ vs 𝒜_emp = 0,33% (≤2%); most 𝒜=C_F²α_s² → α_s=√𝒜/C_F=**0,11792**
  vs PDG 0,1179±0,0009 (**0,03σ**). [Warunkowe na K_geo·m_sp²=π·Φ₀², CG nie domknięte — R5.]
- **T4 (PASS):** formuła HALT-B `m=c_M·A_tail²·g₀^(e²/2)` (wszystkie g₀∈[0,817;0,891]) zawiera
  **0/3** składników dodatekX: brak {addytywne m_0, φ-FP per-sektor, shifted-Koide} ⇒ **strawman**.
- **T5 (PASS):** uczciwy licznik — per sektor **2 wejścia** (m_1, m_2) → m_3 **predykcja**;
  6 mas kwarków = **4 inputy + 2 predykcje** (m_b, m_t) + shared {a_Γ, φ z leptonów}. m_3 NIE
  użyty jako input do własnej predykcji. **r_31 = genuine predykcja, NIE „zero parametrów", NIE fit.**

## §3 — Honest caveats (anti-Lakatos; forbidden #5/#10/R4/R5)

1. **Werdykt CONFIRMED zgodnie z LOCKED regułą** (δ_t≤5% w ≥1 schemacie): m_t trafia w **pole**
   (0,77%), a w **MS-bar** wypada 5,50% (tuż za progiem). **Zgodność m_t jest scheme-zależna** —
   raportuję oba ZANIM ocena (forbidden #5 dotrzymane). Reguła nie została przesunięta ex post.
2. **To NIE „cały SM z 3 inputów"** (forbidden #10). To: **r_31 predykowane z r_21 + uniwersalne 𝒜**.
   Sektor kwarkowy: 4 inputy → 2 predykcje. Słabsze niż leptony (3 masy z 1 inputu via φ+Koide),
   bo r_21 kwarkowe nie jest predykowane (g₀^(1) strojone per sektor).
3. **𝒜=C_F²α_s² to most warunkowy** (CG-1/CG-3 nie domknięte) — nie pełna derywacja. Liczbowo
   uderzający (α_s 0,03σ), ale status [AN+NUM warunkowy].
4. **D1 rozstrzygnięty częściowo:** kanoniczną maszynerią masową kwarków jest **M∝A_tail⁴ +
   addytywne m_0** (dodatekJ/dodatekX), a NIE `A_tail²·g₀^(e²/2)` (formuła HALT-B). Formuła HALT-B
   była błędną reprezentacją — to część diagnozy strawmana.

## §4 — Tabela flag

| Test | Pytanie | Wynik |
|---|---|---|
| T1 | anchors=[0,817;0,891]? (#40) | PASS |
| T2 | m_b≤5% ∧ m_t≤5% (schemat)? | PASS (0,59% / 0,77% pole) |
| T3 | 𝒜 univ ≤2% ∧ α_s ≤1σ? | PASS (0,33% / 0,03σ) |
| T4 | strawman ≥2/3 brak? | PASS (3/3) |
| T5 | genuine predykcja, nie fit? | PASS |
| T6 | werdykt wyliczony | RESCUE-CONFIRMED |

**8/8 FP PASS. 0 hardcoded. Werdykt RESCUE-CONFIRMED (z caveatami §3).**

## §5 — Dyspozycja (Phase FINAL)
1. **HALT-B re-scoped:** testował strawmana (0/3 składników + błędna domena g₀ #40). Werdykt
   IMMUTABLE, ale jego **zakres**: falsyfikuje misformułowanie, NIE sektor kwarkowy TGP.
2. **STATE WIP korekta:** „quark-mass HALT-B (sufit 2,68× vs 80000×)" jako żywy status =
   **fałszywie pesymistyczny**; zastąpić „sektor kwarkowy: r_31 predykowane (m_b 0,6%, m_t 0,8% pole)
   przez φ-FP+m_0=a_Γ/φ+shifted-Koide; HALT-B = strawman".
3. **R12 (sek07):** zawężony — m_b/m_t recovery DZIAŁA; pozostaje otwarte pełne [AN] dla 𝒜 (CG).
4. **PR-025 candidate** (user-gated): predykcja m_b/m_t + α_s z 𝒜.
5. **NIE nadinterpretować** (#10): częściowa predykcja, scheme-zależny m_t, most 𝒜 warunkowy.
