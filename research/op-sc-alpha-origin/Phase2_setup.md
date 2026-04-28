---
title: "Phase 2 setup — α_PB Abrikosov–Gorkov first-principles derivation"
date: 2026-04-28
cycle: SC.1.Phase2
status: PRE-EXECUTION
predecessor: "[[Phase1_results.md]] (4/4 PASS for H1 — α_PB i α_0 strukturalnie różne)"
tags:
  - TGP
  - SC
  - alpha-PB
  - abrikosov-gorkov
  - de-gennes
  - lanthanide-pair-breaking
---

# Phase 2 — Setup: α_PB first-principles derivation z teorii Abrikosov–Gorkov

> **Cel:** Wyprowadzić α_PB ≈ 0.2887 μ_B⁻² z teorii Abrikosov–Gorkov (A-G) +
> standardowych parametrów ciała stałego (N(0), J_sf, g_J) **bez fitowania**
> na T_c PrH₉/NdH₉. Sprawdzić czy α_PB jest **uniwersalną stałą strukturalną
> TGP** czy też **fitowanym artefaktem** specyficznym dla pary PrH₉+NdH₉.

---

## Kontekst (Phase 1 → Phase 2)

Phase 1 zamknęła pytanie unit-bridge negatywnie:
- α_PB ma wymiar **[μ_B⁻²]** = [T²/J²] w SI; α_0 jest **bezwymiarowe**.
- Brak struktury TGP-core (κ_TGP, β, α_em) z naturalnym wymiarem μ_B⁻².
- T1.3 odkryło tylko **numeryczną zbieżność**: α_PB · ⟨μ_eff²⟩_(PrH9,NdH9) ≈ 3.74 ≈ α_0 = 4.04 (ratio 0.926, off ~7.5%).

Phase 2 odpowiada na pytanie: **czy α_PB jest derywowalne z core fizyki SC**
(A-G + de Gennes + literaturowe N(0)) **niezależnie** od fitu T_c?

---

## Hipoteza testowa

| Hipoteza | Twierdzenie | Falsyfikacja |
|---------|-------------|--------------|
| **H_AG** (A-G derivable) | α_PB = (π/4·T_c^base) · N(0) · J_sf² · K (gdzie K to scaling-factor de Gennes lub μ_eff²) z odchyłką < 30% od 0.2887 μ_B⁻² | jeśli α_PB^pred odbiega o > 30% lub scaling-factor jest niespójny dla całej rodziny LnH₉ |
| **H_fit** (independent fit) | α_PB jest **niezależnym fitowanym** parametrem; A-G + de Gennes nie redukuje SC sektor TGP do core | derywacja ZGADZA się — wtedy α_PB awansuje do "derived constant" |

**Threshold akceptacji H_AG:** odchyłka |α_PB^pred / α_PB^obs - 1| < 30%
ORAZ ten sam scaling-factor (de Gennes vs μ_eff²) wewnętrznie spójny dla
przynajmniej 3 z 5 lantanowców rodziny LnH₉.

---

## Ramka teoretyczna

### Standardowa formuła Abrikosov–Gorkov (1960)

Magnetyczne domieszki / lokalne momenty paramagnetyczne łamią pary Coopera
przez spin-flip scattering. Klasyczny wynik A-G:

```
ln(T_c^(0)/T_c) = ψ(1/2 + ρ) - ψ(1/2)
ρ = ℏ / (2π τ_sf k_B T_c) = Γ_sf / (2π k_B T_c)
Γ_sf = (π/ℏ) N(0) J_sf² · S(S+1)              (S(S+1) = "magnetic strength")
```

gdzie:
- T_c^(0) = baseline T_c bez domieszek (≡ T_c^base = 143 K w TGP SC v2)
- N(0) = density of states na poziomie Fermiego (per spin per atom)
- J_sf = exchange coupling f-electron ↔ conduction band (eV)
- ψ = digamma function
- Spin-flip strength **S(S+1)** dla pojedynczego momentu spin S

### Scaling factor — de Gennes vs μ_eff²

W rodzinie metali ziem rzadkich (RE) pair-breaking idzie z **de Gennes factor**:

```
de Gennes factor:  dG = (g_J - 1)² · J(J+1)
```

NIE z `μ_eff² = g_J² · J(J+1) μ_B²`. Różnica jest kluczowa:

| Ln³⁺ | g_J  | J   | dG = (g_J-1)² J(J+1) | μ_eff² (μ_B²) | dG/μ_eff² |
|------|------|-----|----------------------|---------------|-----------|
| La³⁺ | —    | 0   | 0.00                 | 0.00          | —         |
| Ce³⁺ | 6/7  | 5/2 | 0.18                 | 6.43          | 0.028     |
| Pr³⁺ | 4/5  | 4   | 0.80                 | 12.82         | 0.062     |
| Nd³⁺ | 8/11 | 9/2 | 1.84                 | 13.10         | 0.140     |
| Sm³⁺ | 2/7  | 5/2 | 4.46                 | 0.81          | 5.51      |
| Yb³⁺ | 8/7  | 7/2 | 2.57                 | 20.56         | 0.125     |

**Punkt krytyczny:**
- W A-G (z lokalnym RE momentem), pair-breaking ∝ **de Gennes** dG, nie μ_eff².
- W TGP SC v2 (eq:BPB), pair-breaking ∝ **μ_eff²** (przez czynnik exp(-α_PB μ_eff²)).
- Te dwie skale **są jakościowo różne** dla Sm i Yb (różnica nawet ×100).

To daje Phase 2 ostry test:
- **H_AG** poprawne → przewidywanie T_c z dG **różne** od TGP-derived μ_eff² scaling.
- Najwiekszy rozdźwięk: **SmH₉**.

### Identyfikacja TGP ↔ A-G

TGP SC v2 (eq:BPB):  B_PB = exp(-α_PB μ_eff²)
A-G weak-coupling limit (ρ ≪ 1):  T_c/T_c^(0) ≈ 1 - π² ρ /4   ⟹   ln(T_c^(0)/T_c) ≈ Γ_sf/(2k_B T_c)
W mocnym limicie (ρ ≫ 1, czyli T_c ≪ T_c^(0)):  ln(T_c^(0)/T_c) ≈ ln(2γ ρ) - constant

Dla PrH₉/NdH₉: T_c ≈ 5K, T_c^base = 143K → ratio = 0.035 → **strong limit**.

```
ln(143/5) = 3.35    ↔ TGP: α_PB · 12.82 = 3.70  (Pr)
ln(143/4.5) = 3.46  ↔ TGP: α_PB · 13.10 = 3.78  (Nd)
```

**Bardzo blisko**, ale TGP scaling = μ_eff², A-G scaling = dG (de Gennes).
W bezwymiarowej formie:

```
α_PB μ_eff²       (TGP)       ≈ ln(T_c^base / T_c)
π Γ_sf / (4 T_c^base) · K (A-G)  ≈ ln(T_c^base / T_c)
```

gdzie K to scaling-factor (1, dG, lub μ_eff²).

---

## Procedura testów (T2.1 – T2.6)

### T2.1 — A-G symbolic derivation w słabym i silnym limicie (sympy)

Wyprowadzić formuły digamma w obu limitach:
```
weak (ρ ≪ 1):     T_c/T_c^(0) = 1 - (π/4)·ρ + O(ρ²)
strong (ρ ≫ 1):   T_c/T_c^(0) = (π/(2γ))·exp(-π·ρ/2)   (γ = Euler 1.781)
```
Cel: ekstrakować bezwymiarową kombinację `π·Γ_sf/(2·T_c^base)` jako "α_PB · μ_eff²"-equivalent.

**PASS** jeśli sympy zgadza się z literaturą (Maki-de Gennes).

### T2.2 — de Gennes factor analytical dla La/Ce/Pr/Nd/Sm/Yb

Z tabeli Hund GS:
- Ce³⁺ ²F_5/2:  J=5/2, S=1/2, L=3, g_J=6/7, dG=0.18, μ_eff²=6.43
- Pr³⁺ ³H_4:   J=4,   S=1,   L=5, g_J=4/5, dG=0.80, μ_eff²=12.82
- Nd³⁺ ⁴I_9/2: J=9/2, S=3/2, L=6, g_J=8/11, dG=1.84, μ_eff²=13.10
- Sm³⁺ ⁶H_5/2: J=5/2, S=5/2, L=5, g_J=2/7, dG=4.46, μ_eff²=0.81
- Yb³⁺ ²F_7/2: J=7/2, S=1/2, L=3, g_J=8/7, dG=2.57, μ_eff²=20.56

**PASS** jeśli wartości zgadzają się z Ashcroft & Mermin / Jensen & Mackintosh
do 3 cyfr znaczących.

### T2.3 — μ_eff² analytical (Hund) — confirm sources w SC v2

μ_eff² = g_J² · J(J+1) μ_B² (Hund GS, free-ion approximation).

Sprawdzić czy dane z TGP SC v2 paper (eq:BPB context) zgadzają się
z Hund GS dla każdego LnH₉ w bazie 4 wartości literaturowych.

**PASS** jeśli max odchyłka < 5% (kryształowe pole CEF wymusza odchylenia
~1-3% w realistycznych materiałach).

### T2.4 — Scaling-factor ambiguity: dG vs μ_eff²

**Najważniejszy test**: czy fit α_PB · μ_eff² na PrH₉+NdH₉ daje dobre T_c
DLA SmH₉ i YbH₉?

```
Z fitu PrH₉+NdH₉:    α_PB ≈ 0.2887 μ_B⁻²
SmH₉ predict (μ_eff² scale):  ln(T_c^base/T_c) = 0.2887 · 0.81 = 0.234
                              → T_c = 143/exp(0.234) = 113 K   ★ TGP prediction
SmH₉ predict (dG scale):      ln(T_c^base/T_c) = (0.2887/0.80)·dG_Pr · 4.46/0.80
                              = c · 4.46 = 0.36·4.46 = 1.61
                              → T_c = 143/exp(1.61) = 28.7 K
```

WAIT — to jest tyko TGP scaling. Należy oddzielnie policzyć **A-G prediction
z dG** dla SmH₉ używając *tej samej* przeliczeniowej (Pr→Sm):

```
A-G + de Gennes:  Γ_sf ∝ dG
PrH₉ ln(143/5) = 3.35 = c_AG · dG_Pr = c_AG · 0.80   →  c_AG = 4.19
SmH₉:  ln(143/T_c) = 4.19 · 4.46 = 18.7   →   T_c = 143/exp(18.7) ≈ 9·10⁻⁷ K
```

**Critical falsification target:**
- TGP (μ_eff² scaling): SmH₉ T_c ≈ 100 K
- A-G + de Gennes:      SmH₉ T_c ≈ 0 K (catastrophically suppressed)

Eksperyment SmH₉ pod 150 GPa rozstrzyga jednoznacznie.

### T2.5 — Numerical Phase-2 derivation z literaturowych N(0) i J_sf

Inputs (literature):
- N(0) dla LaH₁₀ baseline ≈ 1.5 states/eV/spin/atom (DFT, e.g. Liu et al. 2017)
- J_sf 4f-conduction hybridization w lantanidach pod high-P ≈ 100-200 meV
  (Anderson model estimates, Nakamura et al. 2018)
- T_c^base = 143 K = 0.0123 eV (TGP SC v2)

Formula (A-G, μ_eff² scaling — TGP-style):
```
α_PB^pred = (π/4) · (1/T_c^base) · N(0) · J_sf² · (1/μ_B²)
         = (π/4) · (1/0.0123) · 1.5 · (0.15)² · (1/μ_B²)  [in eV units]
         = 0.785 · 81.3 · 1.5 · 0.0225 · μ_B⁻²
         = 2.15 μ_B⁻²
```

vs observed 0.2887 μ_B⁻². Ratio ≈ 7.5.

**PASS** (H_AG supported) jeśli ratio w [0.7, 1.3]; FAIL jeśli > 2× off.
Aktualne oszacowanie sugeruje **FAIL z A-G + literaturowych J_sf** —
α_PB byłaby raczej "small fraction" pełnego A-G coupling, co znaczy że
albo J_sf jest **istotnie mniejsze** w LnH₉ pod 150 GPa, albo TGP μ_eff²
scaling **nie jest** prosty A-G limit.

### T2.6 — Internal consistency: czy fit-α_PB jest universal w 5 LnH₉?

Z TGP SC v2 paper, fit RMS_log = 0.316 dla 5 LnH₉. To nie jest mały błąd
(~factor 2 w T_c). Sprawdzić czy odchyłki są **systematic** (np. korelują z dG)
czy **random** (true noise w fit).

**PASS** (TGP universal) jeśli no systematic correlation z dG.
**FAIL** (TGP fit-artifact) jeśli residua skorelowane z dG → A-G + de Gennes
byłby lepszy fit.

---

## Materiał wykonawczy

- **Skrypt:** `phase2_alpha_PB_AG_derivation.py` (sympy A-G + numerical)
- **Output:** `phase2_alpha_PB_AG_derivation.txt`
- **Wynik:** `Phase2_results.md`

## Środowisko

Python 3.x + sympy + numpy. Brak literaturowych danych pobieranych online.

```bash
cd TGP/TGP_v1/research/op-sc-alpha-origin
PYTHONIOENCODING=utf-8 python -X utf8 phase2_alpha_PB_AG_derivation.py 2>&1 | tee phase2_alpha_PB_AG_derivation.txt
```

---

## Możliwe wyniki Phase 2

1. **H_AG strongly supported** (ratio 0.8-1.2): α_PB awansuje do "derived
   constant" — TGP SC sektor redukuje się do core + literature inputs.
   Phase 3 staje się walidacją (multi-LnH₉) zamiast kompletnie nowego fitu.

2. **H_AG weakly supported / scaling-factor ambiguous** (ratio 1.5-3, no
   strong dG correlation): α_PB jest A-G-like ale wymaga TGP-specific
   modyfikacji (np. ψ-mediated coupling, nie standard A-G). Phase 3 jako
   krytyczny test SmH₉.

3. **H_AG rejected** (ratio > 5 lub residua silnie skorelowane z dG):
   α_PB pozostaje niezależnym fitowanym parametrem TGP SC sektora.
   Phase 3 staje się falsyfikacyjnym testem TGP μ_eff² scaling vs A-G.

W każdym przypadku **SmH₉ T_c experiment** rozstrzyga: TGP ≈ 100 K vs A-G ≈ 0 K.

---

## Cross-references

- `program.md` — overall 3-phase plan
- `Phase1_setup.md`, `Phase1_results.md` — unit-bridge audit (negative)
- `../../../tgp-sc-paper/paper/tgp_sc.tex` (lines 425-466) — α_PB definition
- A-G original: Abrikosov & Gorkov (1960), Sov. Phys. JETP 12, 1243
- de Gennes scaling: Maple (1976), Appl. Phys. 9, 179
- LnH₉ DFT: Liu, Naumov, Hoffmann et al. (2017), PNAS 114, 6990
