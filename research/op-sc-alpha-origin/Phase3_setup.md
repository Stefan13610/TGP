---
title: "Phase 3 setup — multi-LnH₉ validation of α_PB scaling"
date: 2026-04-28
cycle: SC.1.Phase3
status: PRE-EXECUTION
predecessor: "[[Phase2_results.md]] (H_AG_PARTIAL; SmH₉/YbH₉ falsification targets identified)"
flask: "tgp-sc v2 (DOI 10.5281/zenodo.19670557)"
related:
  - "[[program.md]]"
  - "[[../../INDEX.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
tags:
  - TGP
  - SC
  - alpha-PB
  - multi-LnH9
  - falsification
  - SmH9
  - YbH9
  - PmH9
  - EuH9
  - TbH9
---

# Phase 3 — Setup: multi-LnH₉ walidacja TGP μ_eff² scaling

> **Cel:** Rozszerzyć Phase 2 dyskryminator dG-vs-μ_eff² na **całą rodzinę
> lantanowców** Ln³⁺ (Ln = La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu).
> Dla każdego: (a) przewidzieć T_c pod ~150 GPa według TGP μ_eff² scaling,
> (b) przewidzieć T_c według A-G + de Gennes, (c) zarejestrować predykcję
> w PREDICTIONS_REGISTRY z eksplicytnym horizontem i discriminator factor.

---

## Co Phase 3 ma zrobić (i czego NIE robi)

### Zakres Phase 3 (offline-doable):

1. **Pełna tabela 15 lantanowców** z g_J, J, μ_eff², dG (Hund GS).
2. **Predictions table:** dla każdego LnH₉ (15 mat.), TGP T_c i A-G T_c.
3. **Falsification map:** który Ln daje największy dyskryminator?
4. **Literatura update:** szukanie dotychczas opublikowanych pomiarów
   (SmH₉, EuH₉, TbH₉, etc.) — jeśli istnieją, użyte natychmiast jako test.
5. **PREDICTIONS_REGISTRY.md update**: rejestracja SC4 (SmH₉), SC5 (YbH₉)
   plus inne high-discrimination cele.
6. **RMS_log comparison:** TGP vs A-G+dG na 5 znanych materiałach (z SC v2 paper).

### NIE w Phase 3 (poza zakresem):

- DFT calculation N(0) i J_sf w high-P LnH₉ — wymaga DFT setup, queued.
- Eksperyment SmH₉/YbH₉ — wymaga grupy DAC (Eremets, Hemley, Prakapenka).
- Cyclic recheck T_c^base = 143 K dla każdego Ln — założenie 1st-order.

---

## Hipoteza testowa Phase 3

| Hipoteza | Twierdzenie | Falsyfikacja |
|---------|-------------|--------------|
| **H_TGP** (μ_eff² universal) | TGP fit α_PB = 0.2887 μ_B⁻² działa na wszystkie LnH₉ z RMS_log ≤ 0.5 | jakiekolwiek mat. z |T_c^pred / T_c^obs - 1| > factor 5 i kierunek pasujący do dG-scaling |
| **H_AG** (de Gennes universal) | A-G + dG fit działa na wszystkie LnH₉ | TGP scaling lepszy w RMS_log dla SmH₉ + YbH₉ jednocześnie |

**Critical discrimination:** TGP i A-G + de Gennes **muszą** dawać przeciwne predykcje co najmniej w **dwóch** materiałach. To dwukrotny test eliminujący degenerację 2-punktowego fitu Pr+Nd.

---

## Pełna tabela lantanowców Ln³⁺ (Hund GS)

| Z  | Ln | 4f^n | term | L | S | J | g_J | μ_eff² (μ_B²) | dG = (g_J-1)² J(J+1) |
|----|------|------|---------|---|-----|-----|------|-------|-------|
| 57 | La³⁺ | 4f⁰  | ¹S₀     | 0 | 0   | 0   | —    | 0     | 0     |
| 58 | Ce³⁺ | 4f¹  | ²F₅/₂   | 3 | 1/2 | 5/2 | 6/7  | 6.43  | 0.18  |
| 59 | Pr³⁺ | 4f²  | ³H₄     | 5 | 1   | 4   | 4/5  | 12.80 | 0.80  |
| 60 | Nd³⁺ | 4f³  | ⁴I₉/₂   | 6 | 3/2 | 9/2 | 8/11 | 13.09 | 1.84  |
| 61 | Pm³⁺ | 4f⁴  | ⁵I₄     | 6 | 2   | 4   | 3/5  | 14.40 | 3.20  |
| 62 | Sm³⁺ | 4f⁵  | ⁶H₅/₂   | 5 | 5/2 | 5/2 | 2/7  | 0.71  | 4.46  |
| 63 | Eu³⁺ | 4f⁶  | ⁷F₀     | 3 | 3   | 0   | —    | 0     | 0     |
| 64 | Gd³⁺ | 4f⁷  | ⁸S₇/₂   | 0 | 7/2 | 7/2 | 2    | 63.0  | 15.75 |
| 65 | Tb³⁺ | 4f⁸  | ⁷F₆     | 3 | 3   | 6   | 3/2  | 94.5  | 10.50 |
| 66 | Dy³⁺ | 4f⁹  | ⁶H₁₅/₂  | 5 | 5/2 | 15/2| 4/3  | 113.3 | 7.08  |
| 67 | Ho³⁺ | 4f¹⁰ | ⁵I₈     | 6 | 2   | 8   | 5/4  | 112.5 | 4.50  |
| 68 | Er³⁺ | 4f¹¹ | ⁴I₁₅/₂  | 6 | 3/2 | 15/2| 6/5  | 91.8  | 2.55  |
| 69 | Tm³⁺ | 4f¹² | ³H₆     | 5 | 1   | 6   | 7/6  | 57.2  | 1.17  |
| 70 | Yb³⁺ | 4f¹³ | ²F₇/₂   | 3 | 1/2 | 7/2 | 8/7  | 20.57 | 0.32  |
| 71 | Lu³⁺ | 4f¹⁴ | ¹S₀     | 0 | 0   | 0   | —    | 0     | 0     |

**Obserwacja:** μ_eff² i dG **dramatycznie się różnią** w środku rodziny:
- **Sm³⁺**: μ_eff² małe (0.71), dG bardzo duże (4.46) — ratio 0.16
- **Gd³⁺**: μ_eff² ogromne (63.0), dG ogromne (15.75) — ratio 4.0
- **Yb³⁺**: μ_eff² duże (20.6), dG małe (0.32) — ratio 64.0

Sm³⁺ i Yb³⁺ to **najczystsze** discriminatory (~2 rzędy wielkości w jednym kierunku).

**Eu³⁺ i Lu³⁺** mają J=0 → nie mają suppressowania pair-breaking — TGP-i-A-G zgadzają się: T_c ~ T_c^base = 143 K.

---

## Procedura testów (T3.1 – T3.7)

### T3.1 — Pełna Hund tabela 15 lantanowców

Sympy + analytical g_J, μ_eff², dG dla każdego Ln³⁺. Cross-check vs Jensen & Mackintosh.

**PASS** jeśli każdy Ln³⁺ matches do <0.05.

### T3.2 — TGP T_c predictions dla 15 LnH₉

Z Phase 2 fit slope c_TGP = 0.2631 μ_B⁻² (Pr+Nd 2-point):
```
T_c^TGP(Ln) = T_c^base · exp(-c_TGP · μ_eff²(Ln))
```

PASS:  Tabela 15 wartości generowana correctly.

### T3.3 — A-G + de Gennes T_c predictions dla 15 LnH₉

Z Phase 2 fit slope c_AG = 3.04 (Pr+Nd):
```
T_c^AG(Ln) = T_c^base · exp(-c_AG · dG(Ln))
```

PASS:  Tabela 15 wartości.

### T3.4 — Discrimination factor map

Dla każdego Ln³⁺ obliczyć:
```
discrim(Ln) = log10(T_c^TGP / T_c^AG)
```

Identyfikować top-5 discriminating Ln (największe |discrim|).

**PASS** jeśli top-5 zawiera SmH₉ i YbH₉ + co najmniej jeden inny clean test
(prawdopodobnie GdH₉ lub TbH₉).

### T3.5 — RMS_log comparison na 5 znanych materiałach

Z TGP SC v2 paper (eq:BPB context) używanych był 5 LnH₉ fit:
- LaH₁₀ (T_c ≈ 250 K) — μ_eff = 0
- CeH₉ (T_c ≈ 100 K) — outlier flagged
- PrH₉ (T_c ≈ 5 K)
- NdH₉ (T_c ≈ 4.5 K)
- ThH₉/UH₉/YbH₉ (jeszcze nie pomierzone w v2 — drop)

Compute:
```
RMS_log(TGP) = sqrt(sum_i [ln(T_c^TGP(i) / T_c^obs(i))]²) / sqrt(N)
RMS_log(AG)  = sqrt(sum_i [ln(T_c^AG(i)  / T_c^obs(i))]²) / sqrt(N)
```

W TGP SC v2 paper RMS_log = 0.316 (5 materiałów). Sprawdzić czy A-G + dG bije
to dla tych samych 4 materiałów (LaH₁₀ wykluczone, μ_eff=0 trywialne).

**PASS** jeśli oba scalingi mają RMS_log < 1.0; H_TGP supported jeśli RMS_log(TGP) < RMS_log(AG)+0.1.

### T3.6 — Literature search: czy SmH₉/YbH₉/inne LnH₉ pod high-P są opublikowane?

Web search trzy queries:
1. "SmH9 superconductivity high pressure"
2. "YbH9 superconductivity high pressure"
3. "lanthanide hydride T_c review 2024 2025"

Cele:
- Jeśli SmH₉ pomierzone: czy T_c bliskie 100 K (TGP) czy 0 K (A-G)?
- Jeśli YbH₉ pomierzone: czy T_c bliskie 0 K (TGP) czy 50 K (A-G)?
- Czy są inne LnH₉ (Tb, Dy, Ho, Er, Tm) pomierzone?

**Status:** "PASS" w sensie wykonano search. Wynik (positive/negative/null)
loguje się; nie blokuje closure.

### T3.7 — PREDICTIONS_REGISTRY.md update

Dodać entries:
- **SC4: SmH₉ T_c at ~150 GPa**
  - TGP: T_c ≈ 100 K (μ_eff² = 0.81 Hund) lub 96 K (μ_eff² = 1.5 Van Vleck)
  - A-G: T_c ≈ 0 K
  - Discriminator: factor 10⁵
  - Horizon: experiment 2027–2029 (Eremets/Hemley/Prakapenka groups)
- **SC5: YbH₉ T_c at ~150 GPa**
  - TGP: T_c ≈ 0.65 K
  - A-G: T_c ≈ 54 K
  - Discriminator: factor 80 (opposite direction)
  - Horizon: experiment 2027–2029
- **SC6 (optional): GdH₉ T_c at ~150 GPa**
  - TGP: extremely small (μ_eff² = 63)
  - A-G: extremely small (dG = 15.75)
  - Discriminator: small (both predict T_c → 0)
  - Horizon: not high-priority

**PASS** jeśli registry updated z 2 (Sm, Yb) lub 3 (Sm, Yb, Gd) entries.

---

## Materiał wykonawczy

- **Skrypt:** `phase3_multi_LnH9_validation.py` (sympy + numerical)
- **Output:** `phase3_multi_LnH9_validation.txt`
- **Wynik:** `Phase3_results.md`

## Środowisko

Python 3.x + sympy + numpy. Plus WebSearch (literature scan).

```bash
cd TGP/TGP_v1/research/op-sc-alpha-origin
python -X utf8 phase3_multi_LnH9_validation.py 2>&1 | tee phase3_multi_LnH9_validation.txt
```

---

## Możliwe wyniki Phase 3

1. **TGP wins clean (RMS_log(TGP) < RMS_log(AG))**: μ_eff² scaling preferowany
   na obecnych danych. Phase 3 closure CONFIRMED, czeka na SmH₉.

2. **A-G + dG wins clean**: TGP SC v2 prediction TENTATIVELY UNDERMINED.
   Wymaga update flask paper z addendum (Phase 2 + 3 corrections).

3. **Tied (both RMS_log similar)**: 2-point fit symmetric — dyskryminator
   tylko z **nowych** pomiarów (SmH₉/YbH₉). Phase 3 closure CONFIRMED.

W wszystkich 3 przypadkach: PREDICTIONS_REGISTRY zaktualizowane,
2-3 nowe formal predictions zarejestrowane z timestamps.

---

## Cross-references

- `program.md` — overall 3-phase plan
- `Phase1_results.md`, `Phase2_results.md` — predecessor cycles
- `phase2_alpha_PB_AG_derivation.py` — Phase 2 numerical machinery
- `../../PREDICTIONS_REGISTRY.md` — output target dla T3.7
- TGP SC v2 paper, Sec.5 (eq:BPB), figure 6 (5-LnH₉ fit) — RMS_log baseline
- Jensen & Mackintosh, "Rare Earth Magnetism" (1991), table 1.2 — Hund GS values
