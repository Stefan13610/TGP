---
title: "op-sc-alpha-origin — Structural origin of α_PB"
date: 2026-04-28
cycle: SC.1
status: ACTIVE
predecessor: "[[../closure_2026-04-26/alpha_psi_threshold/results.md]] (T-α: α_0 ≈ 4.04 dimensionless)"
flask: "tgp-sc v2 (DOI 10.5281/zenodo.19670557)"
related:
  - "[[../../INDEX.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
tags:
  - TGP
  - SC
  - alpha-PB
  - structural-origin
  - lanthanide
---

# op-sc-alpha-origin — Structural origin of α_PB

> **Cel:** Sprawdzić, czy stała Pauli–Bohr α_PB ≈ 0.2887 μ_B⁻² (TGP SC v2, eq:BPB) jest **niezależnym fitowanym parametrem** czy też **derywowalna** z bardziej fundamentalnych wielkości — albo z TGP core (κ_TGP, β, α_0 z T-α), albo ze standardowego Abrikosov–Gorkov (Γ_sf, N(0), g_J).

---

## Tło

W TGP SC v2 (paper `eq:BPB`):

```
B_PB(μ_eff) = exp(-α_PB μ_eff²)
α_PB = 0.2887 μ_B⁻²    (fitted on PrH₉ + NdH₉)
```

α_PB jest jedną z **trzech "uniwersalnych stałych strukturalnych"** w T_c-formule (obok κ_TGP ≈ 2.012 i β ≈ 2.527), ale w przeciwieństwie do κ_TGP/β jej pochodzenie z TGP core nie jest spelled out — paper podaje tylko **2-punktowy fit** (PrH₉ T_c = 5 K, NdH₉ T_c = 4.5 K, bazowe T_c^base = 143 K).

W tym samym czasie z closure_2026-04-26 mamy **drugą α** w teorii: α_0 ≈ 4.04 (T-α threshold), bezwymiarową, opisującą substrate–matter coupling przy ψ = 1. Pytanie strategiczne: czy jest jakaś relacja między α_PB a α_0?

---

## 3-fazowy plan

### Phase 1 — Unit-bridge audit (TYLKO sprawdzenie wymiarowe)

**Cel:** Czy α_PB i α_0 to *ta sama* fizyczna α w różnych jednostkach (np. naturalne ↔ μ_B⁻²)?

**Hipoteza nullowa H₀:** TAK — istnieje czynnik konwersji jednostek `f` taki że `α_PB = f · α_0`.

**Hipoteza alternatywna H₁:** NIE — α_PB i α_0 to **strukturalnie różne** stałe (różne jednostki, różna fizyka, różne calibration).

**Testy:**
- T1.1 — analiza wymiarowa (sympy `dimsys_SI`)
- T1.2 — porównanie układów odniesienia (geom-units vs SI)
- T1.3 — sprawdzenie czy T-α residue ma reprezentację z μ_B⁻² jednostkach
- T1.4 — eksplicytne sprawdzenie bilansu: czy α_PB · ⟨μ_eff⟩² ↔ α_0 · ⟨ψ - 1⟩² dla jakiejś naturalnej kalibracji

**Status:** zamknij negatywny wynik w 1 cyklu jeśli H₀ odpada wymiarowo (oczekiwane). Phase 1 closure note → stwierdzenie że α_PB **nie jest** unit-cousin α_0 i wymaga osobnego strukturalnego pochodzenia.

### Phase 2 — Abrikosov–Gorkov first-principles derivation (jeśli Phase 1 = negatywne)

**Cel:** Wyprowadzić α_PB ze standardowej teorii Abrikosov–Gorkov + TGP core constants, **bez fitu na T_c**.

**Strategia:**
- Standard A-G: `ln(T_c^(0)/T_c) = ψ(1/2 + Γ_sf/(2π T_c)) - ψ(1/2)`
- W limicie Γ_sf >> T_c: T_c ≈ T_c^(0) · exp(-π Γ_sf / (4 T_c^(0)))
- Γ_sf ∝ μ_eff² · N(0) · τ_sf (standardowa Born approximation w spin-flip channel)
- Identyfikacja: α_PB · μ_eff² = π Γ_sf / (4 T_c^base) → α_PB = π · N(0) · τ_sf / (4 T_c^base) · const

**Predykcja:** α_PB powinno **wyjść z**:
- T_c^base = 143 K (już fixed by κ_TGP, β z core)
- N(0) — gęstość stanów na poziomie Fermiego dla rodziny LnH₉ (literaturowe dane DFT)
- τ_sf — czas spin-flip scattering (literaturowe NMR/μSR)
- g_J — Hund Landé factor (analytical)

**Falsyfikacja:** jeśli `α_PB^pred` z above formula odbiega o > 30% od `α_PB^fit = 0.2887`, to TGP **nie redukuje** SC do core constants → α_PB pozostaje niezależnym fitem.

### Phase 3 — Multi-material validation (jeśli Phase 2 daje predykcję)

**Cel:** Sprawdzić Phase-2 pochodzenie na **dodatkowych** lantanowych nadprzewodnikach poza PrH₉/NdH₉:
- LaH₁₀ (μ_eff = 0, expected T_c ≈ 200 K, baseline)
- CeH₉ (μ_eff = 2.54, mixed valence — outlier)
- SmH₉ (μ_eff = 1.55, named prediction T_c ≈ 100 K)
- YbH₉ (μ_eff = 0, expected T_c, η-limited)

**Sukces:** RMS_log Phase-2-derived ≤ RMS_log original fit (= 0.316 z 5 LnH₉).

**Falsyfikacja:** systematic outlier ≥ 50% w którymkolwiek niewykluczonym materiale.

---

## Co robimy najpierw

**Phase 1 — TERAZ.** Mała, deterministyczna analiza wymiarowa. Zamknięcie 1-cyklowe, niezależnie od wyniku.

Pozostałe fazy queued — wymagają literatury (DFT N(0) dla LnH₉) lub user-go-ahead.

---

## Cross-references

- `paper/tgp_sc.tex` — eq:BPB, def:BPB, fit description (lines 425-466)
- `closure_2026-04-26/alpha_psi_threshold/results.md` — α_0 ≈ 4.04 origin
- `INDEX.md` — flask deposits map
- `PREDICTIONS_REGISTRY.md` — SC sector predictions LuH10, Hg1245/SrTiO3, FeSe/BaTiO3, YbH9, SmH9

---

## Decyzja po Phase 1

- Jeśli Phase 1 = H₀ rejected (oczekiwane): pisz negatywny audyt, zamknij cykl SC.1.Ph1, propose Phase 2 do user-a.
- Jeśli Phase 1 = H₀ supported (mało prawdopodobne): rozszerz audyt o explicit derivation α_PB = f · α_0, raport do user-a.
