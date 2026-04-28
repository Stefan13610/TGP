---
title: "Phase 1 results — α_PB ↔ α_0 unit-bridge audit (negative result)"
date: 2026-04-28
cycle: SC.1.Phase1
status: CLOSED
verdict: 4/4 PASS for H1 (structurally distinct)
predecessor: "[[Phase1_setup.md]]"
related:
  - "[[program.md]]"
  - "[[phase1_unit_bridge_audit.py]]"
  - "[[phase1_unit_bridge_audit.txt]]"
  - "[[../closure_2026-04-26/alpha_psi_threshold/results.md]]"
tags:
  - TGP
  - SC
  - alpha-PB
  - dimensional-analysis
  - negative-result
  - closure
---

# Phase 1 results — α_PB ↔ α_0 unit-bridge audit

> **Status:** ✅ **CLOSED 2026-04-28** — **4/4 PASS for H1** (structurally distinct).
> α_PB ≈ 0.2887 μ_B⁻² (TGP SC v2 lanthanide pair-breaking) i α_0 ≈ 4.04
> (T-α substrate-matter coupling) **nie są** unit-cousins — to **strukturalnie różne** stałe wymagające osobnych derywacji.

---

## TL;DR

**Hipoteza H₀** (unit cousin): α_PB = f · α_0 dla pewnego pure-unit czynnika f → **REJECTED**.

**Hipoteza H₁** (structurally distinct): α_PB i α_0 to różne fizyczne stałe (różne sektory: SC vs gravity; różne wymiary: μ_B⁻² vs dimensionless; różne calibration sources: PrH₉/NdH₉ vs photon-ring) → **SUPPORTED 4/4**.

**Nieoczekiwana obserwacja (Phase 2 hook):** `α_PB · ⟨μ_eff²⟩_LnH₉ ≈ 3.74` zaskakująco blisko `α_0 ≈ 4.04` (off ~7.5%). To może być:
- (a) **numerical coincidence** — bo ⟨μ_eff²⟩ jest średnią z **tych samych** PrH₉/NdH₉ z których fitowano α_PB, więc nie jest niezależnym estymatorem.
- (b) **strukturalna zbieżność** — jeśli dla **innych** lantanowców (LaH₁₀, CeH₉, SmH₉, YbH₉) ten sam iloczyn `α_PB · ⟨μ_eff²⟩` mieści się ~5% wokół α_0, to jest realna relacja TGP-core ↔ SC.

Phase 2 (queued) zajmie się rozstrzygnięciem (a) vs (b).

---

## Wyniki test-by-test

### T1.1 — Pełna analiza wymiarowa (sympy `dimsys_SI`) ✅ H₀ REJECTED

```
dim(mu_B)     = {length: 2, current: 1}        [SI: A·m²]
dim(alpha_PB) = {length: -4, current: -2}      [SI: A⁻²·m⁻⁴]
dim(alpha_0)  = {}                              [dimensionless]
dim(alpha_PB) == dim(alpha_0)?  False
```

**Wniosek:** α_PB ma niezerowy wymiar SI (`A⁻²·m⁻⁴`); α_0 jest bezwymiarowe. **Żadna pure-unit konwersja nie może je zidentyfikować** — to fundamentalna nierównoważność wymiarowa.

---

### T1.2 — Geometric-unit equivalence check ✅ H₀ REJECTED

Test: gdyby α_PB i α_0 były unit cousins, w SI powinno być:

```
alpha_PB^pred = alpha_0 / mu_B² = 4.04 / (9.27e-24 J/T)² ≈ 4.70e+46 J⁻²T²
alpha_PB^observed_SI                                      = 3.36e+45 J⁻²T²
ratio (observed/pred)                                     = 0.072  (off factor 14)
```

Albo równoważnie, w jednostkach μ_B⁻² (= setting μ_B = 1):

```
alpha_0 (in mu_B⁻²)  = 4.04          [if alpha_0 -> alpha_0 · 1/mu_B²]
alpha_PB             = 0.2887        [observed]
ratio                = 0.0715        (off factor 14)
```

**Wniosek:** nawet po wymuszeniu identyfikacji unit-system (μ_B = 1), wartości różnią się o **czynnik 14**. Brak naturalnej kalibracji TGP daje α_PB z α_0.

---

### T1.3 — μ_eff calibration coincidence ⚠️ NUMERICAL_COINCIDENCE_NOT_STRUCTURAL

Średnia ⟨μ_eff²⟩ z PrH₉ + NdH₉:

```
mu_eff(PrH₉) = 3.58 μ_B,  μ_eff² = 12.82
mu_eff(NdH₉) = 3.62 μ_B,  μ_eff² = 13.10
avg(μ_eff²) = 12.96
```

Iloczyn:

```
alpha_PB · avg(mu_eff²) = 0.2887 · 12.96 = 3.74        (dimensionless)
alpha_0                  = 4.04                          (dimensionless)
ratio                    = 3.74 / 4.04 = 0.926           (off ~7.5%)
```

**Interpretacja:** liczba `3.74` jest blisko α_0 ≈ 4.04, ale `<μ_eff²>` jest średnią z **dokładnie tych** dwóch lantanowców, na których fitowano α_PB. To nie jest niezależny estymator. Konstrukcyjnie iloczyn `α_PB · μ_eff²` zostaje **fitowane** wraz z α_PB.

**Test który mógłby tę zbieżność strukturalnie potwierdzić** (Phase 3): czy ten sam iloczyn `α_PB · μ_eff²` dla **niezależnych** materiałów z podobnym zachowaniem A-G (z literaturowych Γ_sf, niezależnie od TGP-fit) zachowuje się jako stała ≈ α_0?

**Wniosek:** to **numerical coincidence** w obecnym audicie — strukturalnie nieuzasadniona.

---

### T1.4 — Structural identity criterion ✅ H₁ SUPPORTED

Próg akceptacji H₀ (structural identity):
- (a) closed-form `α_PB = f(TGP core constants only)`
- (b) deviation < 0.1% od observed value

W obecnej literaturze TGP **brak** closed-form relacji `α_PB = f(α_0, κ_TGP, β)`. α_PB jest **fitowany** (2-point, PrH₉ + NdH₉, z T_c^base = 143 K). Jego pochodzenie z core constants pozostaje **otwartym problemem** (Phase 2).

**Wniosek:** kryterium structural identity nie jest spełnione. **H₁ supported** — α_PB i α_0 są strukturalnie różne.

---

## Implikacje strategiczne

### Co Phase 1 USTALA (negatywny wynik konstruktywny)

1. **α_PB nie jest unit-aliasem α_0.** Każda przyszła praca która próbuje "zredukować" α_PB do α_0 przez konwersję jednostek jest skazana na porażkę — to jest niemożliwe wymiarowo.

2. **TGP ma minimum dwie niezależne α-stałe:** α_0 (gravitational substrate-matter coupling) z core, α_PB (lanthanide pair-breaking) z SC sector. Każda wymaga osobnej derywacji albo osobnego fitu.

3. **Audit-trail value:** ten audyt zamyka "low-hanging fruit" (czy oba αs to ten sam parametr) negatywnie, oszczędzając przyszłemu auditorowi cykli przemyśleń.

### Co Phase 1 NIE rozstrzyga (queued for Phase 2/3)

1. **Czy α_PB można wyprowadzić z TGP core + standard A-G** (κ_TGP, β, N(0), τ_sf, g_J)? Phase 2.
2. **Czy zbieżność z T1.3 (`α_PB · ⟨μ_eff²⟩ ≈ α_0`) jest strukturalna czy przypadkowa?** Phase 3 (multi-material extension).

### Phase 2 plan (queued)

W limicie Γ_sf >> T_c, Abrikosov-Gorkov daje:
```
T_c ≈ T_c^(0) · exp(-π Γ_sf / (4 T_c^(0)))
Γ_sf ∝ μ_eff² · N(0) · τ_sf       [Born approx, spin-flip channel]
```

Identyfikacja z TGP `B_PB = exp(-α_PB · μ_eff²)`:
```
α_PB = (π / 4 T_c^base) · N(0) · τ_sf · g²(z)
```

gdzie:
- `T_c^base = 143 K` z TGP core (już pinned by κ_TGP, β),
- `N(0)` z DFT (literatura, materiał-specyficzne),
- `τ_sf` z μSR/NMR (literatura),
- `g(z)` Lande g-factor (analytical Hund GS).

Predykcja: `α_PB^pred(LnH₉ family) = (π/4·143K) · ⟨N(0)·τ_sf⟩_LnH₉ · g²` — to jest fully closed-form **bez fitowania na T_c**.

**Falsyfikacja Phase 2:** jeśli `α_PB^pred / α_PB^observed = 1 ± 0.30`, structural derivation supported. Jeśli odbiega > 30%, α_PB pozostaje niezależnym fitem w SC sector.

---

## Cross-references

- `program.md` — full 3-phase plan
- `Phase1_setup.md` — pre-execution test plan
- `phase1_unit_bridge_audit.py` — sympy + numerical script
- `phase1_unit_bridge_audit.txt` — execution log (4/4 PASS)
- `../../INDEX.md` — flask deposits map
- `../closure_2026-04-26/alpha_psi_threshold/results.md` — α_0 origin
- `../../../tgp-sc-paper/paper/tgp_sc.tex` (lines 425-466) — α_PB definition

---

## Drift-check

| Kryterium | Wynik |
|-----------|-------|
| Single-Φ axiom zachowane | ✓ (T1.1 wymiarowo niezależna od fizycznej zawartości) |
| Brak nowych input parameters | ✓ (audyt tylko porównuje istniejące stałe) |
| TGP core constants nie modyfikowane | ✓ (κ_TGP, β, α_0 niezmienione) |
| SC v2 flask nie modyfikowany | ✓ (α_PB = 0.2887 μ_B⁻² pozostaje fitted; flask deposit at DOI 10.5281/zenodo.19670557 unchanged) |
| Audit completes without falsifying anything | ✓ (jedynie classifies α_PB ↔ α_0 as not-unit-cousins) |

**Brak drift-u od fundamentów TGP_v1.**

---

## Decyzja po Phase 1

**Phase 1 CLOSED 2026-04-28.** 4/4 PASS, hipoteza H₁ (structurally distinct) supported.

Następne kroki:
- **Phase 2 (queued):** spróbować closed-form derivation α_PB z A-G + TGP core. Wymaga literaturowych N(0)·τ_sf dla LnH₉ family (DFT + μSR/NMR). Czas: ~1-2 cykle.
- **Phase 3 (queued):** test Phase-2 derivation na poszerzonym 5-7-materiałowym sample (LaH₁₀, CeH₉, PrH₉, NdH₉, SmH₉, YbH₉). Wymaga eksperymentalnych T_c. Czas: ~1 cykl.

Phase 1 wniosek dodany do INDEX.md / PREDICTIONS_REGISTRY.md w kolejnym update zachowuje audit-trail.

---

*Phase 1 closure: 2026-04-28. Cykl op-sc-alpha-origin Phase 1 zamknięty z verdict 4/4 PASS dla H₁.*
