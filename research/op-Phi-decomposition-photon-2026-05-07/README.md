---
title: "op-Phi-decomposition-photon — Stage 2: foton jako mod δΦ"
date: 2026-05-07
type: research-cycle
status: SETUP_PHASE_0_PENDING_REVIEW
folder_status: active
classification: STRUCTURAL_EXTENSION_CANDIDATE
parent: "[[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]]"
related_audit:
  - "[[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]"
  - "[[../../audyt/L08_kink_fermion_closure/README.md]]"
  - "[[../../audyt/S07_M911_derivation/README.md]]"
related_research:
  - "[[../op-FRW-radiation-era-varying-c-2026-05-06/FINDINGS.md]]"
tgp_owner: research/op-Phi-decomposition-photon-2026-05-07
tags:
  - research-cycle
  - Stage2
  - photon-ontology
  - phi-decomposition
  - path-alpha
  - post-EXT1-pivot
---

# op-Phi-decomposition-photon — Stage 2 setup

## Geneza cyklu

EXT-1 ścieżka A (FRW radiation era z varying c bez dekompozycji
Φ = Φ̄ + δΦ) zakończona **STRUCTURAL_NO_GO**
([[../op-FRW-radiation-era-varying-c-2026-05-06/FINDINGS.md]]).

Diagnoza pivotowa: **Path α** — surgiczne rozdzielenie Φ na:

```
Φ(x,t) = Φ̄(t) + δΦ(x,t)
        ‾‾‾‾‾‾   ‾‾‾‾‾‾‾‾
        tło      perturbacja
        (kosmo)  (lokalna fizyka)
```

To rozwiązuje konflikt operacyjny w obecnym frameworku gdzie Φ jest
naprzemiennie traktowane jako:
- (i) globalne tło kosmologiczne
- (ii) lokalne źródło grawitacji (M9.1'')
- (iii) mod oscylacyjny (foton jako zaburzenie)

Bez explicit dekompozycji równania mieszają reżimy → Φ-EOM
w erze radiacyjnej staje się ill-defined (Phase 1 EXT-1: scenariusz (a)
ψ frozen, scenariusz (b) divergence).

## Hipoteza centralna Stage 2

**H1**: foton jest propagującym modem δΦ na tle Φ̄.

**Konsekwencje testowalne:**
- Z liniowego δΦ-EOM otrzymujemy `ω² = c²k²` (zwykła dyspersja
  bezdyspersyjna)
- Po kanonicznej kwantyzacji: `E = ℏω`, `p = ℏk` → `λ = hc/E` ✓
- c **wynika z tła Φ̄**, nie z amplitudy δΦ → fotony różnych energii
  lecą z tym samym c (rozwiązuje pytanie self-interaction)
- Polaryzacja: **OPEN PROBLEM** — skalarny δΦ ma 1 DOF, foton 2 transverse
  → wymaga modu wyższych pochodnych albo strukturalnego rozszerzenia

## Trzy reżimy δΦ (proponowana ontologia)

| Reżim | δΦ charakter | Manifestacja | Powiązanie audyt |
|-------|--------------|--------------|------------------|
| **Falowy** | ∂²δΦ/∂t² = c²∇²δΦ, oscylujące | Foton, promieniowanie | this cycle |
| **Solitonowy** | ∇δΦ ≠ 0, ∂t δΦ = 0, topologicznie chronione | Cząstka spoczynkowa (masa) | L08 kink-fermion |
| **Zerowy** | δΦ = 0 | Próżnia (samo Φ̄) | trivial |

**Wave-particle duality** wynika naturalnie: te same δΦ mogą być
delokalizowane (fala) albo zlokalizowane (cząstka), w zależności od
struktury topologicznej.

## Plan fazowy

### Phase 0 — Balance sheet (MANDATORY pre-derivation)
**Status:** [[Phase0_balance.md]] PENDING_REVIEW

8 ☑ gate per [[../../meta/CALIBRATION_PROTOCOL.md]] §2.

### Phase 1 — Formalna dekompozycja Φ = Φ̄ + δΦ

**Cele:**
- F1.1: zdefiniować separację skal (background vs perturbation)
- F1.2: wyprowadzić linearized δΦ-EOM z istniejącego pełnego Φ-EOM
  (`K_geo·□ψ + V'(ψ) = source`)
- F1.3: zidentyfikować że c² ~ K(Φ̄)/Φ̄ jest funkcją tła **only**
- F1.4: pokazać dyspersję ω² = c²k² + corrections (rzędu (δΦ/Φ̄)²)
- F1.5: weryfikacja sympy

**Output:** `Phase1_results.md` — formal δΦ-EOM + dyspersja

### Phase 2 — Foton jako mod δΦ

**Cele:**
- F2.1: kanoniczna kwantyzacja δΦ → operator pola
- F2.2: derivation E=ℏω, p=ℏk z ω=ck + zasad komutacyjnych
- F2.3: derivation λ=hc/E (już szkicowo udowodnione w Stage 2 intuicyjnym)
- F2.4: stress-energy T_μν fotonu jako mod δΦ
- F2.5: weryfikacja: T^μ_μ = 0 (Weyl-niezmienniczość 4D, spójność z L01)

**Output:** `Phase2_results.md` — photon ontology formal

### Phase 3 — Problem polaryzacji (KRYTYCZNY)

**Cele:**
- F3.1: zliczenie DOF skalarnego δΦ → 1 mod skalarny
- F3.2: alternatywy:
  - α: pojedyncza polaryzacja podłużna (problematyczne obserwacyjnie)
  - β: polaryzacja z gradientu ∇δΦ (3 DOF wektorowych, ale podłużny excluded)
  - γ: TT-mod ∂_i∂_j δΦ (6 → 2 DOF transversal-traceless)
- F3.3: empirical falsifier: dwie polaryzacje obserwowane **must** być
  reprodukowane
- F3.4: weryfikacja sympy + dyspersja

**Output:** `Phase3_results.md` — polarization mechanism (lub fail)

**Decision tree:**
- A) γ-path (TT-mod) działa → Stage 2 SUCCESS, kontynuujemy
- B) α-path (longitudinal only) → strukturalna porażka, Stage 2 FAIL
- C) β-path (gradient) marginal → STRUCTURAL CONDITIONAL

### Phase 4 — Powrót do BBN ρ_rad mechanism (CONDITIONAL on Phase 3 SUCCESS)

**Cele:**
- F4.1: jeśli foton = mod δΦ, to ρ_rad = energy density of δΦ-modes
- F4.2: re-derivation Friedmann eq w erze radiacyjnej z explicit
  rozdzieleniem Φ̄ (tło) + ρ_rad (δΦ-modes)
- F4.3: czy to **ratuje** EXT-1 ścieżkę A?
- F4.4: BBN drift ponownie — czy <5%?

**Output:** `Phase4_results.md` — EXT-1 retroactive review

**Decision tree:**
- A) BBN drift <5% → EXT-1 status revised, ścieżka A retroactively SAVED
- B) BBN drift 5-50% → STRUCTURAL CONDITIONAL na konkretne assumption
- C) BBN drift > 50% → EXT-1 STRUCTURAL_NO_GO **utrzymany**, ale
  Stage 2 photon ontology jako **niezależny rezultat strukturalny**

### Phase 5 — Cross-cycle integration (jeśli Phase 4 SUCCESS lub partial)

- L08 kink-fermion → matter ontology consistency
- S07 M9.1'' derivation → metric consistency post-decomposition
- Aktualizacja core LaTeX (sek04, sek08a, sek08c)

### Phase 6 — Honest reporting baseline (ABSOLUTE BINDING)

Per [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 6 gate enforced.

## Risk assessment

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| R1: dekompozycja Φ̄+δΦ łamie axiomy istniejące | 10-20% | HIGH | sympy weryfikacja Phase 1, możliwy re-open S04/S07 |
| R2: polaryzacja nie wychodzi z δΦ skalarnego | 40-60% | HIGH | path γ (TT-mod) jest non-trivial; ewentualny pivot do δΦ tensorowego |
| R3: BBN nadal FAILS post-decomposition | 30-40% | MEDIUM | EXT-1 STRUCTURAL_NO_GO już zaakceptowany — Stage 2 ma własną wartość niezależnie |
| R4: cykl wykryje deeper inconsistency M9.1'' | 5-10% | KRYTYCZNY | re-open S07 + cascade audit |

## Subiektywna probabilistyczna ocena (pre-Phase-1)

| Outcome | Probability |
|---------|-------------|
| Stage 2 → DERIVED FULL (foton + λ + polaryzacja) | 15-25% |
| Stage 2 → STRUCTURAL CONDITIONAL (foton + λ, polaryzacja partial) | 35-45% |
| Stage 2 → STRUCTURAL_NO_GO (polaryzacja blokuje) | 25-35% |
| Stage 2 → ratuje EXT-1 retroactively | 10-20% |
| Stage 2 → Stage 2 success ALE EXT-1 nadal STRUCTURAL_NO_GO | 35-45% |

## Cross-references

- [[Phase0_balance.md]] — MANDATORY pre-Phase-1
- [[NEEDS.md]] — open questions per phase
- [[../op-FRW-radiation-era-varying-c-2026-05-06/FINDINGS.md]] — geneza pivotu
- [[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]
- [[../../audyt/L08_kink_fermion_closure/README.md]] — soliton/kink ontology
- [[../../audyt/S07_M911_derivation/README.md]] — metric derivation ze substratu
- [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 6 ABSOLUTE BINDING gate
- [[../../core/sek04_stale/sek04_stale.tex]] — c-derivation source (linie 178-208)
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] — Φ-EOM
