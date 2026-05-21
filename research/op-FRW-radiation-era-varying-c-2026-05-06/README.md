---
title: "op-FRW-radiation-era-varying-c — H_TGP(z) z varying c, ℏ, G dla BBN/CMB"
date: 2026-05-06
parent: "[[../INDEX.md]]"
status: CLOSED-SUPERSEDED  # post-2026-05-16 closure ceremony per Phase_FINAL_close.md
closure_date: 2026-05-16
folder_status: closed-superseded
claim_status: STRUCTURAL_NO_GO  # Path A scope; REDIRECTED-RESOLUTION via Path F successor cycle
successor: "[[../op-inflation-substrate-genesis-2026-05-11/Phase_FINAL_close.md]]"  # Path F A− (2026-05-13)
cycle: EXT-1 (extension L01)
audit_source: "[[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2"
priority: P1 OTWARTE RYZYKO → CLOSED-SUPERSEDED
related:
  - "[[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]"
  - "[[../op-L01-rho-stress-energy-bridge-2026-05-04/]]"
  - "[[../../core/sek04_stale/sek04_stale.tex]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]]"
flask: "TBD (Zenodo deposit po Phase 3)"
tags:
  - TGP
  - EXT-1
  - L01-extension
  - FRW
  - radiation-era
  - varying-constants
  - BBN
  - CMB
  - cosmology
  - falsification-test
tgp_status:
  folder_status: active
  level: L1
  kind: derivation
  core_compatibility: extension
  last_reviewed_against_core: 2026-05-06
  may_edit_core: false
  exports_findings: true
  has_needs_file: true
  has_findings_file: true
  open_bridges:
    - "ax:c-ax:G integration"
    - "L01 NEEDS N1 quantum EM trace"
    - "L01 NEEDS N2 QCD vacuum"
    - "M9.1'' validity w reżimie ψ << 1"
  depends_on:
    - "L01 EXECUTED 2026-05-04 (formal definition ρ ≡ -T^μ_μ/c_0²)"
    - "ax:c-ax:G axiomy (sek04_stale.tex)"
    - "M9.1'' canonical form (sek08a/sek08c)"
  impacts:
    - "PRIORITY_MATRIX EXT-1 P1 OTWARTE RYZYKO"
    - "TGP fenomenologia kosmologiczna (BBN, CMB)"
    - "L01 NEEDS N1, N2 promotion"
  source_of_status:
    - "audyt/EXTERNAL_REVIEW_2026-05-06.md §EXT-1 v2"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: 2026-05-06
---

# op-FRW-radiation-era-varying-c — H_TGP(z) z varying c, ℏ, G dla BBN/CMB

> **Cel:** Zamknąć **EXT-1 OTWARTE RYZYKO** (recenzja zewnętrzna
> 2026-05-06): pokazać explicit policzenie H_TGP(z) w erze radiacyjnej
> z varying c(Φ), ℏ(Φ), G(Φ), porównanie z BBN i CMB. **Decyzja
> dwustronna** — fenomenalny sukces (BBN+CMB w 1-2 wolnych params)
> lub silny sygnał konieczności pivotu L_mat.

## Kontekst i motywacja

### EXT-1 v2 diagnoza (z recenzji zewnętrznej)

[[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]:

L01 formal definition (cykl 2026-05-04) zamknął:
```
ρ(x) ≡ −T^μ_μ(x) / c_0²
```

→ **T^μ_μ_EM = 0** (Weyl-niezmienniczość 4D Lagrange'a EM) — strukturalnie
fotony NIE generują Φ przez L_mat.

ALE: aksjomaty `ax:c–ax:G` (`sek04_stale.tex` lin. 29, 250–254):
```
c(Φ) = c₀ √(Φ₀/Φ)     ℏ(Φ) = ℏ₀ √(Φ₀/Φ)     G(Φ) = G₀ Φ₀/Φ
```

→ c, ℏ, G NIE są stałymi fundamentalnymi, lecz **funkcjami pola Φ**.

**Dwustronność:** H_TGP(z) w erze radiacyjnej **NIE musi być równe
H_GR(z)** — pochodzi z Φ-EOM w FRW z varying constants, niezależnie
od standardowego mechanizmu ρ_rad ~ a⁻⁴.

### Open question

**Czy TGP odzyskuje fenomenologię BBN/CMB przez varying constants
zamiast przez ρ_rad?**

## Hipoteza nullowa H₀ vs alternatywna H₁

**H₀:** H_TGP(z) z varying c, ℏ, G **NIE** odzyskuje fenomenologii BBN+CMB
w 1-2 wolnych parametrach (γ, Φ₀ skalibrowane do innych obserwabli).

**H₁:** H_TGP(z) odzyskuje BBN+CMB w 1-2 wolnych parametrach z dokładnością:
- BBN ⁴He: |H_TGP(z=10⁹) - H_GR(z=10⁹)|/H_GR < **5%** (BBN tolerance)
- CMB l_1 = **220.0 ± 0.5** (pierwszy peak)
- r_s(z_drag) = **147 Mpc** (Planck consistent)
- N_eff = **3.046** (relativistic species)

**Subiektywna ocena recenzenta:** P(H₁) = **55-65%** (niepewność, NIE
przegrane a priori).

## 3-fazowy plan

### Phase 1 — Φ-EOM w FRW background z varying c, ℏ, G (5 sub-tests)

**Cel:** Wyprowadzić H_TGP(t, ψ) z Friedmann equation w TGP z explicit
ax:c-ax:G integration.

- **F1.1** — FRW background setup z ψ(t) cosmological evolution
- **F1.2** — Derivation H_TGP² = (8πG(Φ)/3c(Φ)²)·ρ_total + Φ-kinetic terms
- **F1.3** — ax:c(Φ), ℏ(Φ), G(Φ) substitution w równaniach
- **F1.4** — Limity asymptotyczne: radiation-dominated era (z >> z_eq)
- **F1.5** — Phase 1 GATE ≥4/5 PASS

**Falsyfikacja:** Jeśli Phase 1 ujawnia, że Φ-EOM w FRW nie ma
sensownego rozwiązania w erze radiacyjnej (np. singularność,
non-perturbative regime), Phase 2-3 abandon.

### Phase 2 — BBN + recombination kinematics (7 sub-tests)

**Cel:** Numerical computation H_TGP(z) dla z ∈ [10³, 10¹⁰] + porównanie
z BBN abundance + sound horizon.

- **F2.1** — H_TGP(z) numerical solver (z = 10³ do 10¹⁰)
- **F2.2** — H_TGP(z=10⁹) vs H_GR(z=10⁹) — drift < 5% gate
- **F2.3** — BBN ⁴He abundance Y_p (varying ℏ → Coulomb barrier modyfikacja)
- **F2.4** — D/H, ³He/H, ⁷Li/H abundance
- **F2.5** — σ_T·n_e·c kinematics (Saha equation, Thomson scattering, c⁻⁴)
- **F2.6** — sound horizon r_s(z_drag) — Planck 147 Mpc gate
- **F2.7** — Phase 2 GATE ≥6/7 PASS

**Falsyfikacja:** Jeśli żaden z (γ, Φ₀) ∈ realistic band nie daje BBN
zgodności w 5% AND sound horizon w 1%, ścieżka A FAILS → ścieżka D
(L_mat extension, narusza S04) lub ścieżka E (przyznanie zakresu).

### Phase 3 — CMB peaks + N_eff + falsification roadmap (6 sub-tests)

**Cel:** Pierwsza peak CMB l_1, N_eff, multi-channel falsification roadmap.

- **F3.1** — CMB pierwszy peak l_1 = 220.0 ± 0.5 gate
- **F3.2** — N_eff = 3.046 relativistic species
- **F3.3** — r_s(z_drag) cross-check Plancka 147 Mpc
- **F3.4** — alt-cosmologies cross-check (ΛCDM benchmark, z innymi varying-c models)
- **F3.5** — 4-channel falsification roadmap (BBN + CMB peaks + r_s + N_eff)
- **F3.6** — Phase 3 GATE ≥4/6 + decisive verdict (DERIVED / STRUCTURAL / FAIL)

**Falsyfikacja:** Phase 3 ≥4/6 + 4-channel convergence ≥3/4 → DERIVED;
3/6 + 2-3/4 → STRUCTURAL; <3/6 lub <2/4 → FAIL → ścieżka D/E.

## Decision tree post-Phase 3

```
Phase 1 + 2 + 3 → DECIZIJA TRZYSTANOWA:

(A) DERIVED: BBN 5% + CMB l_1 ±0.5 + r_s 147 Mpc + N_eff 3.046
    z 1-2 wolnymi parametrami
    → TGP odzyskuje fenomenologię standardową przez varying constants
    → EXT-1 ZAMKNIĘTE jako PHENOMENAL SUCCESS
    → L01 NEEDS N1, N2 demoted z P1 do P3

(B) STRUCTURAL: częściowa zgodność (np. BBN OK, CMB partial)
    → Wymagana dodatkowa analiza per-channel
    → EXT-1 STRUCTURAL CONDITIONAL — extension cyklu

(C) FAIL: brak zgodności w żadnej kombinacji (γ, Φ₀)
    → Pivot konieczny: ścieżka D (L_mat extension, narusza S04)
       lub ścieżka E (przyznanie zakresu post-recombination)
    → EXT-1 STRUCTURAL_NO_GO + decyzja autora pivot
```

## Phase 0 balance sheet (CALIBRATION_PROTOCOL §2 — MANDATORY)

Patrz [[Phase0_balance.md]].

Phase 6 ABSOLUTE BINDING gate (post-2026-05-06) wymaga **Phase0_balance.md
PRZED registry commit**. Plik gotowy z:
- External inputs (BBN PDG, CMB Planck, r_s, N_eff)
- Structural axioms TGP-internal LOCKED (L01, ax:c-ax:G, M9.1'')
- Derived outputs (H_TGP, Y_p, l_1, r_s, N_eff)
- Tautology test (czy varying constants nie kasują się?)
- Falsifiability test (explicit gates per-channel)
- Independent path (4 channels minimum)

## Co cykl NIE robi (z definicji scope)

- **Nie pivot L_mat** (ścieżka D rezerwowa, narusza S04 — wymaga osobnej decyzji)
- **Nie modyfikuje core/** (M03 framework: research-level, NIE core edits)
- **Nie rozstrzyga L01 NEEDS N1 (quantum EM trace anomaly)** — to osobny cykl
- **Nie rozstrzyga L01 NEEDS N2 (QCD vacuum)** — to osobny cykl
- **Nie pisze pełnego CAMB/CLASS port pod TGP** — to byłoby Phase 4+ (deferred)

## Open questions (trackable)

Patrz [[NEEDS.md]] — formalna lista N1...N5 do rozszerzenia po każdej Phase.

## Cross-references

- [[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2 — recenzja źródłowa
- [[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]] — folder audytowy EXT-1
- [[../../audyt/PRIORITY_MATRIX.md]] — EXT-1 P1 OTWARTE RYZYKO
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/]] — L01 EXECUTED 2026-05-04 + NEEDS N1, N2
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]] — N1 (quantum EM trace), N2 (QCD vacuum)
- [[../../core/sek04_stale/sek04_stale.tex]] — ax:c–ax:G axiomy
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]
  L01 formal definition lin. 148-183
- [[../../core/sek05_ciemna_energia/sek05_ciemna_energia.tex]] — sek05 cosmology
- [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 6 ABSOLUTE BINDING gate
- [[../op-M03-balance-sheet-retrofit-2026-05-06/template_Phase0_balance.md]] — template
- [[../op-M03-balance-sheet-retrofit-2026-05-06/Phase5_FULL_IMPLEMENTATION_2026-05-06.md]]

## Decyzja po Phase 1

- **Jeśli Phase 1 = 4/5 PASS** (Φ-EOM w FRW well-defined):
  → Phase 2 enabled (BBN + recombination)
- **Jeśli Phase 1 = <4/5** (singularność, non-perturbative regime):
  → Cykl ABANDONED, ścieżka D/E required, raport do user-a

## Status FINAL (post-2026-05-16 closure ceremony)

🟢 **CLOSED-SUPERSEDED** — patrz [[Phase_FINAL_close.md]] dla pełnej ceremonii zamknięcia.

**Werdykt cyklu:**
- **Path A** (varying-constants recovery w obecnym ax:c-ax:G framework): **STRUCTURAL_NO_GO**
  confirmed numerycznie w Phase 2 (H_TGP/H_GR ≈ 0.184%; Y_p_TGP ≈ 0.31% vs PDG 24.5%).
- **Path F** (pre-BBN inflation, NOVEL physics): **SUCCESS** via dedicated successor cycle
  [[../op-inflation-substrate-genesis-2026-05-11/]] (A−, 41/41 sympy PASS, F3 Starobinsky
  preferred Planck-compatible, Φ_eq chain across 6 cosmological epochs, 2026-05-13).
- **Path D** (L_mat extension + S04 re-open): **NIE PODJĘTA — niepotrzebna**; Path F
  resolved scope without naruszenia S04 (B9 MICROSCOPE preserved bezwarunkowo).
- **Path E** (TGP_FOUNDATIONS scope acknowledgment): **PENDING-COSMETIC** (deferred do
  dedicated housekeeping cycle); zreframed post-Path F success.

**EXT-1 P1 OTWARTE RYZYKO status:** 🟢 **CLOSED-SUPERSEDED** via Path F successor cycle.
