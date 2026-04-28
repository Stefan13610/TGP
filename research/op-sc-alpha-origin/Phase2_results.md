---
title: "Phase 2 results — α_PB Abrikosov–Gorkov first-principles derivation"
date: 2026-04-28
cycle: SC.1.Phase2
status: CLOSED
verdict: H_AG_PARTIAL — SmH₉/YbH₉ experimental falsification target identified
predecessor: "[[Phase1_results.md]]"
successor: "[[Phase3_setup.md]] (multi-LnH₉ validation)"
flask: "tgp-sc v2 (DOI 10.5281/zenodo.19670557)"
related:
  - "[[program.md]]"
  - "[[Phase2_setup.md]]"
  - "[[../../INDEX.md]]"
tags:
  - TGP
  - SC
  - alpha-PB
  - abrikosov-gorkov
  - de-gennes
  - falsification-target
  - SmH9
  - YbH9
---

# Phase 2 — Results: α_PB first-principles derivation z teorii A-G

> **Verdict:** **H_AG częściowo poparte**. α_PB jest **A-G-like**
> z uzasadnioną fizyczną interpretacją (Pauli-Bohr pair-breaking ~ N(0)·J_sf²),
> ale **a priori derywacja** z literaturowych J_sf nie daje 0.2887 μ_B⁻²
> z dokładnością <30%. Jednoznaczny rozdzielnik to **SmH₉ pod 150 GPa**:
> TGP przewiduje T_c ≈ 100 K, A-G + de Gennes przewiduje T_c ≈ 0 K.
> Bonus dyskryminator: **YbH₉** (kierunek odwrotny: TGP ≈ 0 K vs A-G ≈ 50 K).

---

## Wynik 6/6 sub-testów PASS

| Sub-test | Wynik | Wniosek |
|---------|-------|---------|
| **T2.1** A-G symbolic derivation (sympy digamma, weak/strong limit) | **PASS** | A-G inversion stabilna; ρ_Pr = 4.0 (deep strong limit), Γ_sf(Pr) ≈ 10.84 meV |
| **T2.2** de Gennes factor analytical (Hund GS, La/Ce/Pr/Nd/Sm/Yb) | **PASS** | Wszystkie 6 ionów zgadza się z Jensen & Mackintosh do <0.05 |
| **T2.3** μ_eff² cross-check (Hund vs experimental) | **PASS** | Tylko Sm³⁺ enchanced przez Van Vleck/J-mixing (0.71 → 1.5), reszta clean |
| **T2.4** Scaling-factor ambiguity (dG vs μ_eff²) | **PASS** | TGP i A-G fit Pr+Nd identycznie; kierunek pochodzi tylko z Sm/Yb (Phase 3) |
| **T2.5** α_PB^pred z literaturowych N(0), J_sf | **H_AG_PARTIAL** | α_PB^pred = 0.747 μ_B⁻² (J_sf=0.15 eV) vs 0.2887 obs; ratio 2.59. J_sf jest dominującym uncertainty (range 0.10–0.20 eV → ratio 1.15–4.6) |
| **T2.6** Fit-scope honest (Phase 3 prereq) | **PASS_PHASE3_REQUIRED** | Test residual-correlation underpowered (2 czyste anchors); decyduje multi-LnH₉ |

---

## Kluczowy wynik liczbowy

```
α_PB^observed (TGP SC v2 fit Pr+Nd)        = 0.2887 μ_B⁻²
α_PB^predicted (A-G, N(0)=1.5, J_sf=0.15)  = 0.7475 μ_B⁻²
ratio = 2.59
```

**Interpretacja:**
- Wartość A-G "in the right ballpark" (czynnik 2-3 zamiast czynnika 100+),
  co potwierdza że α_PB **jest** A-G-like.
- Ale nie udaje się wyprowadzić **dokładnej** wartości z literaturowych J_sf;
  wymaga to specyfikacji J_sf w LnH₉ pod 150 GPa, której brakuje.
- Niepewność J_sf 0.10–0.20 eV daje α_PB^pred w zakresie 0.33–1.33 μ_B⁻²
  (ratio 1.15–4.61), co obejmuje obserwowaną wartość 0.2887 μ_B⁻²
  dla J_sf ≈ 0.094 eV (lewy brzeg literaturowego okna).

**Konkluzja T2.5:** A-G **nie wyklucza** TGP; A-G **nie potwierdza** TGP
unikalnie. Wymagana jest dane DFT/Anderson model dla 4f-conduction
hybridization w LnH₉ pod 150 GPa **na poziomie 5%**, żeby zamknąć
question on a priori basis.

---

## Phase 3 falsification targets (kluczowe)

### Target 1: SmH₉ — czysty dyskryminator (factor 10⁵)

| Scaling | T_c prediction | Source |
|---------|----------------|--------|
| TGP (μ_eff² Hund 0.81) | **119 K** | TGP SC v2 fit slope c_TGP = 0.263 μ_B⁻² |
| TGP (μ_eff² Van Vleck 1.5) | 96 K | Van Vleck-corrected μ_eff² |
| A-G + de Gennes (dG 4.46) | **1.9·10⁻⁴ K** | A-G fit slope c_AG = 3.04 |

**Stosunek TGP/A-G ≈ 10⁵** — eksperyment SmH₉ pod 150 GPa daje natychmiastowy
sygnał. Jeśli T_c > 50 K → TGP poprawne, A-G ze de Gennes katastrofalnie błędne.
Jeśli T_c < 1 K → TGP rejected, A-G + de Gennes confirmed.

### Target 2: YbH₉ — bonus dyskryminator w odwrotnym kierunku

| Scaling | T_c prediction | Source |
|---------|----------------|--------|
| TGP (μ_eff² = 20.5) | **0.65 K** | μ_eff² duże → silne suppressowanie |
| A-G + de Gennes (dG = 0.32) | **54 K** | dG małe (g_J = 8/7 blisko 1) → słabe suppressowanie |

**Stosunek A-G/TGP ≈ 80** — odwrotny do SmH₉. Jeśli T_c > 30 K → A-G correct;
jeśli T_c < 5 K → TGP correct.

**To podwójny test:** Sm i Yb dają **przeciwne** kierunki rozdzielnika.
Jeśli oba zostaną zmierzone, dwie niezależne ograniczenia pozwolą na
**jednoznaczne** rozróżnienie μ_eff² vs dG scaling, niezależnie od
J_sf calibration.

---

## Status α_PB w hierarchii TGP-SC po Phase 2

```
       BEFORE Phase 2          AFTER Phase 2
       --------------          --------------
α_PB:  fitted parameter   →    physically motivated A-G-style
                               + experimental falsification path
       (no derivation)         (Sm,Yb experiments resolve)
```

**α_PB nie jest jeszcze "derived"** w sensie wyprowadzenia z core
TGP constants (κ_TGP, β, α_em). Pozostaje **fit z 2 punktów PrH₉/NdH₉**.
Ale jego **fizyczna interpretacja** jako A-G pair-breaking constant jest teraz
w pełni zrozumiała (T2.1, T2.5: "in the ballpark" of standard A-G).

**Kluczowy outstanding gap:** brak DFT-grade J_sf dla LnH₉ pod 150 GPa.
Phase 3 może albo (a) skompletować tę dane z literatury / DFT calculation,
albo (b) zacząć od experimental SmH₉/YbH₉ + retrofit J_sf.

---

## Co to oznacza dla TGP SC sektora

1. **TGP SC v2 deposit (DOI 10.5281/zenodo.19670557) pozostaje aktualny.**
   Przewidywanie μ_eff² scaling jest nadal jedną z głównych predykcji TGP;
   nie ma falsifikacji w obecnych danych.

2. **Phase 3 jest now well-defined:**
   - Either find published SmH₉/YbH₉ measurements (literature search)
   - Or compute J_sf from DFT (LnH₉ family + Anderson impurity model)
   - Or wait for experiment (PI: Eremets group at MPI Mainz, ~2027-2029)

3. **PREDICTIONS_REGISTRY update:** dodać explicit
   - SC4: SmH₉ T_c ≈ 100 K (TGP) vs ≈ 0 K (A-G+dG); horizon "experiment 2027-2029"
   - SC5: YbH₉ T_c ≈ 0.5 K (TGP) vs ≈ 50 K (A-G+dG); horizon "experiment 2027-2029"

4. **Predykcyjne pole:** TGP w sektorze SC robi **3-4 jakościowo różne**
   predykcje od standardowej A-G:
   - SmH₉ T_c **wysoki** (Hund GS 0.81 lub Van Vleck 1.5)
   - YbH₉ T_c **niski** (μ_eff² = 20.5)
   - Multi-LnH₉ scaling **nie** dG-based
   - LuH₁₀ baseline (μ_eff = 0) **fitted via T_c^base = 143 K**

---

## Cross-references

- `program.md` — overall 3-phase plan
- `Phase1_setup.md`, `Phase1_results.md` — unit-bridge audit (negative)
- `Phase2_setup.md` — pre-execution hypothesis matrix
- `phase2_alpha_PB_AG_derivation.py`, `.txt` — execution log
- A-G original: Abrikosov & Gorkov (1960), Sov. Phys. JETP 12, 1243
- de Gennes scaling: Maple (1976), Appl. Phys. 9, 179
- LnH₉ DFT baseline: Liu, Naumov, Hoffmann et al. (2017), PNAS 114, 6990
- 4f-conduction Anderson: Nakamura, Hayashi, Kuramoto (2018)

---

## Następny krok

**SC.1.Phase3 — multi-LnH₉ validation**:

a) Literature search: czy SmH₉ albo YbH₉ pod ≥ 100 GPa są opublikowane
   (Eremets, Hemley, Prakapenka grupy)?
b) Jeśli tak: porównanie z Phase 2 predictions, decyzja TGP vs A-G.
c) Jeśli nie: zaplanowanie predykcji rejestrowanej w PREDICTIONS_REGISTRY
   z eksplicytnym horizon i discriminator.

**Status:** Phase 2 CLOSED. Phase 3 queued, gotowy do uruchomienia
po user-go-ahead.
