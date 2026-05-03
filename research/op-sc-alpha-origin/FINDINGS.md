---
title: "FINDINGS — op-sc-alpha-origin"
date: 2026-05-03
parent: "[[README.md]]"
type: findings
tgp_owner: research/op-sc-alpha-origin
source_session: S5 (auto-extraction; cite-only per AGENT_PROTOCOL §3)
tags:
  - findings
---

# FINDINGS — op-sc-alpha-origin

> **Sesja 5 auto-generated** (2026-05-03). Ekstrakcja cite-only z istniejących
> plików folderu. **Żadne treści nie są wymyślone** — każdy item ma `source:`
> cytujący plik. Manual review w Sesji 6 (RESEARCH_BUS broadcasts).

## Phase results — frontmatter verdicts

| Plik | Cycle | Status | Verdict | Score | Program |
|------|-------|--------|---------|-------|---------|
| `Phase1_results.md` | `SC.1.Phase1` | `CLOSED` | `4/4 PASS for H1 (structurally distinct)` | — | — |
| `Phase2_results.md` | `SC.1.Phase2` | `CLOSED` | `H_AG_PARTIAL — SmH₉/YbH₉ experimental falsification target identified` | — | — |
| `Phase3_results.md` | `SC.1.Phase3` | `CLOSED` | `7/7 PASS — multi-LnH₉ falsification matrix REGISTERED; SC.1 program END` | — | — |

## Numerical PASS counts (cited from .txt outputs)

- `phase1_unit_bridge_audit.txt` — 4/4 PASS — context: …OT_STRUCTURAL   T1.4  structural identity criterion  : H1_SUPPORTED    Verdict: 4/4 PASS for H1 (structurally distinct).            alpha_PB and alpha_0 are NO

## TL;DR / Wynik / Verdict sections (cytaty)

### `Phase1_results.md` — TL;DR

> ## TL;DR
> 
> **Hipoteza H₀** (unit cousin): α_PB = f · α_0 dla pewnego pure-unit czynnika f → **REJECTED**.
> 
> **Hipoteza H₁** (structurally distinct): α_PB i α_0 to różne fizyczne stałe (różne sektory: SC vs gravity; różne wymiary: μ_B⁻² vs dimensionless; różne calibration sources: PrH₉/NdH₉ vs photon-ring) → **SUPPORTED 4/4**.
> 
> **Nieoczekiwana obserwacja (Phase 2 hook):** `α_PB · ⟨μ_eff²⟩_LnH₉ ≈ 3.74` zaskakująco blisko `α_0 ≈ 4.04` (off ~7.5%). To może być:
> - (a) **numerical coincidence** — bo ⟨μ_eff²⟩ jest średnią z **tych samych** PrH₉/NdH₉ z których fitowano α_PB, więc nie jest niezależnym es...(truncated)

### `Phase1_results.md` — Wyniki test-by-test

> ## Wyniki test-by-test
> 
> ### T1.1 — Pełna analiza wymiarowa (sympy `dimsys_SI`) ✅ H₀ REJECTED
> 
> ```
> dim(mu_B)     = {length: 2, current: 1}        [SI: A·m²]
> dim(alpha_PB) = {length: -4, current: -2}      [SI: A⁻²·m⁻⁴]
> dim(alpha_0)  = {}                              [dimensionless]
> dim(alpha_PB) == dim(alpha_0)?  False
> ```
> 
> **Wniosek:** α_PB ma niezerowy wymiar SI (`A⁻²·m⁻⁴`); α_0 jest bezwymiarowe. **Żadna pure-unit konwersja nie może je zidentyfikować** — to fundamentalna nierównoważność wymiarowa.
> 
> ---
> 
> ### T1.2 — Geometric-unit equivalence check ✅ H₀ REJECTED
> 
> Test: gdyby α_PB i α_0 były uni...(truncated)

### `Phase2_results.md` — Wynik 6/6 sub-testów PASS

> ## Wynik 6/6 sub-testów PASS
> 
> | Sub-test | Wynik | Wniosek |
> |---------|-------|---------|
> | **T2.1** A-G symbolic derivation (sympy digamma, weak/strong limit) | **PASS** | A-G inversion stabilna; ρ_Pr = 4.0 (deep strong limit), Γ_sf(Pr) ≈ 10.84 meV |
> | **T2.2** de Gennes factor analytical (Hund GS, La/Ce/Pr/Nd/Sm/Yb) | **PASS** | Wszystkie 6 ionów zgadza się z Jensen & Mackintosh do <0.05 |
> | **T2.3** μ_eff² cross-check (Hund vs experimental) | **PASS** | Tylko Sm³⁺ enchanced przez Van Vleck/J-mixing (0.71 → 1.5), reszta clean |
> | **T2.4** Scaling-factor ambiguity (dG vs μ_eff²) | **PASS** |...(truncated)

### `Phase2_results.md` — Status α_PB w hierarchii TGP-SC po Phase 2

> ## Status α_PB w hierarchii TGP-SC po Phase 2
> 
> ```
>        BEFORE Phase 2          AFTER Phase 2
>        --------------          --------------
> α_PB:  fitted parameter   →    physically motivated A-G-style
>                                + experimental falsification path
>        (no derivation)         (Sm,Yb experiments resolve)
> ```
> 
> **α_PB nie jest jeszcze "derived"** w sensie wyprowadzenia z core
> TGP constants (κ_TGP, β, α_em). Pozostaje **fit z 2 punktów PrH₉/NdH₉**.
> Ale jego **fizyczna interpretacja** jako A-G pair-breaking constant jest teraz
> w pełni zrozumiała (T2.1, T2.5: "in the ballpark"...(truncated)

### `Phase3_results.md` — Wynik 7/7 sub-testów PASS

> ## Wynik 7/7 sub-testów PASS
> 
> | Sub-test | Wynik | Wniosek |
> |---------|-------|---------|
> | **T3.1** Hund GS table 15 Ln³⁺ | **PASS** | Wszystkie 15 lantanowców z analytical g_J, μ_eff², dG, zgodność z Jensen & Mackintosh |
> | **T3.2** TGP T_c predictions | **PASS** | c_TGP = 0.2631 μ_B⁻², tabela 15 wartości |
> | **T3.3** A-G + de Gennes T_c predictions | **PASS** | c_AG = 3.036, tabela 15 wartości |
> | **T3.4** Discrimination factor map | **PASS** | Top-3: Gd³⁺ (10¹³·⁶, both → 0), Er³⁺/Ho³⁺ (~10⁷, A-G > TGP), **Sm³⁺ (10⁵·⁸, TGP > A-G)** |
> | **T3.5** RMS_log on known data | **TGP_BETTER** | TGP ...(truncated)

### `Phase3_results.md` — Status α_PB w hierarchii TGP po Phase 3

> ## Status α_PB w hierarchii TGP po Phase 3
> 
> ```
> After Phase 1:  α_PB nie jest unit-cousin α_0  (H₀ rejected)
> After Phase 2:  α_PB JEST A-G-like, ale a priori J_sf ratio 2.59  (H_AG_PARTIAL)
> After Phase 3:  TGP μ_eff² scaling preferred over A-G+dG na istn. danych (RMS 0.42 vs 1.5);
>                 cleanest tests = SmH₉ (TGP wins) + YbH₉ (A-G wins) — pre-registered.
> ```
> 
> **SC.1 program (3 fazy) END.**
> - Phase 1: 4/4 PASS (unit-bridge audit, H₁ supported)
> - Phase 2: 6/6 PASS (A-G derivation audit, H_AG_PARTIAL)
> - Phase 3: 7/7 PASS (multi-LnH₉ falsification map, registered)
> - **Total: 17 sub-tests...(truncated)

---

## Cross-references

- [[README.md]] — opis folderu + YAML status
- [[NEEDS.md]] — otwarte luki tego folderu
- [[meta/research/RESEARCH_BUS.md]] — broadcast tych findings
- [[meta/research/FOLDER_STATUS_INDEX.md]] — globalna mapa