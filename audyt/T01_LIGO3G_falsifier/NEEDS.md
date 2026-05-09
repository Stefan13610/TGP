---
title: "NEEDS — T01 LIGO 3G falsifier statement (|Δg_tt| = (5/6) U³)"
date: 2026-05-07
parent: "[[README.md]]"
type: needs
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - needs
  - audit
  - test-falsification
  - LIGO-3G
  - Einstein-Telescope
  - Cosmic-Explorer
  - M911
  - 3PN
  - ppE
  - T01
  - EXT-5
related:
  - "[[README.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[../EXTERNAL_REVIEW_2026-05-06.md]]"
  - "[[../PRIORITY_MATRIX.md]]"
  - "[[../S07_M911_derivation/README.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
---

# NEEDS — T01 LIGO 3G falsifier statement

> **Cel pliku.** Sprowadzić problem T01 ("|Δg_tt| = (5/6) U³ vs LIGO 3G")
> z opisu w `README.md` do **listy konkretnych, wykonalnych zadań** —
> z explicit lokalizacjami danych źródłowych, kandydatami dostawców,
> typem pracy i kryterium domknięcia. NEEDS jest **biletem do cykli
> badawczych**, nie samym cyklem. Audit-folder NIE wykonuje rachunków
> obserwacyjnych; przygotowuje *contract* dla cykli `research/op-…`.

## Stan pre-existing (co JEST w TGP_v1)

| Element | Lokalizacja | Status |
|---------|-------------|--------|
| Wyprowadzenie analityczne c_n (n=2..7) z α=2 vacuum Φ-EOM | [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] §2.1, §3.2 | DERIVED (sympy 5/5) |
| Rozwinięcie g_tt^TGP do U⁶ | M9_1_pp_P1_results.md §3.2 tabela | DERIVED |
| Różnica TGP − GR przy U³, U⁴, U⁵, U⁶ (explicit liczby) | M9_1_pp_P1_results.md §3.2 tabela | DERIVED |
| 1PN match β=γ=1 EXACT | M9_1_pp_P1_results.md §4.1 | LOCKED (P23 sympy 5/5) |
| Tabela detektywności Mercury / LLR / Cassini / GW170817 / EHT (skala U³) | M9_1_pp_P1_results.md §4.3 | OOM only |
| Zapis w `TGP_FOUNDATIONS.md` § 3 | [[../../TGP_FOUNDATIONS.md]] §3 ("M9.1'' przewiduje explicit \|Δg_tt\| = (5/6) U³ deviation od GR") | NARRATIVE |
| Zapis w sek08c | [[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]] M9.1'' P1 | NARRATIVE |

**Kluczowy fakt:** współczynniki PN są **już policzone** (M9_1_pp_P1).
T01 NIE musi wyprowadzać deviation od zera — musi (a) skonkretyzować
ją do *strain h(f)*, (b) zmapować na ppE, (c) zarejestrować jako
falsifier.

## Stan brakujący (co T01 wymaga zamknąć)

### Otwarte luki

| ID | Luka | Typ | Kandydat dostawcy | Path z README |
|----|------|-----|-------------------|---------------|
| N1 | Mapowanie c_n_TGP − c_n_GR → modyfikacja fazy waveformu Ψ(f) w stationary phase approximation (SPA) dla nieobrotowego BBH | derivation | nowy cykl `research/op-LIGO-3G-deviation/` Phase 1 | A |
| N2 | Mapowanie analityczne (5/6) U³ → (β_ppE, b_ppE) w Yunes–Pretorius parameterization (oraz (α_ppE, a_ppE) dla amplitudy) | derivation | nowy cykl `research/op-ppE-mapping/` Phase 1; preview: [[PPN_TO_PPE_MAPPING.md]] | B |
| N3 | Parameter degeneracy 1D (β_ppE_TGP) vs Fisher matrix dla ET-D / CE design sensitivities; chirp-mass × spin × β_ppE kowariancja | numerical | `op-LIGO-3G-deviation/` Phase 2; korzystać `pyCBC` lub `bilby` ppE plugin | A |
| N4 | Konkretne SNR-thresholds dla detekcji deviation > 5σ przy fixed M_chirp (10–60 M_⊙ BBH; 1.4 M_⊙ BNS) i fixed distance (40 Mpc / 400 Mpc / 4 Gpc) | numerical | `op-LIGO-3G-deviation/` Phase 3 | A |
| N5 | Wyższe rzędy: czy U⁴, U⁵, U⁶ deviations (też explicit w M9_1_pp_P1) trafiają w 4PN+ kontrolowane przez ET/CE? | analytical | `op-LIGO-3G-deviation/` Phase 4 — może zostać dodane jako sub-falsifiers `M911-P2`, `M911-P3` | A |
| N6 | Decyzja waveform-model: GR-base + (β_ppE) deviation (perturbative inspiral) **vs** pełny TGP numerical-relativity (NR-TGP — nie istnieje) | decision (autor) | autor + S07 status | A, decision |
| N7 | Reżim ważności: dla jakich (M_BBH, frequency band) U ≪ 1 trzyma się tak, że perturbative ppE jest zaufana (non-merger inspiral only)? Z explicit konwencją PN counting (zob. [[CONVENTION_DECISION.md]]: PHASE convention adoptowana, b_ppE = −1 dla M911-P1) | analytical | `op-LIGO-3G-deviation/` Phase 1 cutoff analysis | A |
| N8 | Konkretny wpis falsifier w `PREDICTIONS_REGISTRY.md` z ID `M911-P1` (lub `T01-FALSIFIER`), kategoria FALSIFIABLE, threshold X·10^(−Y) | registry edit | autor (decyzja); preview gotowy: [[FALSIFIER_STATEMENT_DRAFT.md]] | C |
| N9 | Decyzja: czy T01 publikować jako standalone short paper (Class. Quantum Grav. comm., ~6 stron) czy w ramach większego paperu o M9.1'' (PRD, full) | decision (autor) | autor | D |
| N10 | Powiązanie z S07 (derywacja M9.1'' z fundamentu) — gdy S07 zamknie się ścieżką A/B/D, T01 awansuje z "test ansatzu" do "test derywowanej predykcji" | structural | [[../S07_M911_derivation/README.md]] | meta |
| N11 | Spójność z Path B w sek08b ghost-free basenie ψ ∈ (0, 4/3): czy strong-field U → 1/2 jeszcze trzyma perturbative PPN, czy wymaga full-NR | analytical | `op-LIGO-3G-deviation/` Phase 1 cutoff (zob. też N7) | A |

### Pytania otwarte

- **Q1:** Jaki jest *najsłabszy* sygnał ET/CE, dla którego (5/6)U³
  daje detekcję 5σ? Konkretne: M_BBH, dłuższe inspirale (t_obs > 1h
  w ET), distance, network SNR.
  - Source: brak existing Fisher analysis dla TGP-specific β_ppE
  - Wpływa na: falsifier threshold X·10^(−Y) w PREDICTIONS_REGISTRY

- **Q2:** Czy (β_ppE_TGP = β_TGP, b_ppE_TGP = +1) jest zdegenerowane
  z którąś z 5 standardowych ppE modyfikacji (Bessel-like, scalar
  charge, dipole radiation, masive graviton, Lorentz-violating)?
  - Source: Yunes–Pretorius 2009; Yagi–Yunes 2016 review
  - Wpływa na: identyfikowalność (5/6)U³ jako *specifically TGP* signal

- **Q3 — CLOSED 2026-05-07** via [[../../research/op-GWTC3-reanalysis/Phase3_verdict.md]]:
  - **Verdict:** TGP M911-P1 jest **CONSISTENT z GWTC-3 within 1σ**;
    BF_TGP/GR ≈ 0.97 INCONCLUSIVE; 0.37σ deviation z observed
    centroid.
  - **Detection power (generic ppE-marginalized analysis):** ~230,000
    BBH events needed for 5σ — far beyond projected 3G era rates.
  - **CRITICAL caveat:** detection w 3G era wymaga **dedicated TGP-
    specific Bayes prior** (single-coefficient β = -5/64), nie
    generic LIGO ToGR multi-coefficient marginalization.
  - **Recommended next-session:** Workflow A (TGP-specific Bayes
    z public GWTC-3 events, ~1 sesja Python) z expected ~1-2σ
    tentative signal.
  - Source: Abbott et al. 2023 PRX 13:041039 + Phase3_verdict.md §3

- **Q4:** Czy U³ deviation w `g_tt` propaguje się 1:1 do fazy waveformu,
  czy istnieją kompensacje przez równolegle modyfikowane g_ij i
  energy-flux dE/dt? (ppE benchmark robi to w SPA dla orbitalnej
  prędkości; należy sprawdzić w pełnej M9.1'' kinematyce.)
  - Source: brak — to jest core analytic content N1
  - Wpływa na: czy "(5/6)U³" w metryce → "(5/6)·c_phase" w fazie
    z dokładnością do O(1) prefactor (oczekiwane)

- **Q5:** Czy zapis aksjomatyczny `c² = c₀²·(4-3ψ)²/ψ²` w M9.1''
  (zamiast standardowego c² = c₀²/ψ) modyfikuje propagację GW (c_T)
  inaczej niż przewidują GW1 (c_T = c_s exactly) i GW6 (m_s ≪ ω_LIGO)?
  - Source: M9_1_pp_setup.md §2.5; PREDICTIONS_REGISTRY GW1, GW6
  - Wpływa na: spójność T01 z istniejącymi predykcjami GW

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | T01 jest "test ansatzu" dopóki S07 nie zamknie M9.1'' z fundamentu (ścieżka A/B/D w [[../S07_M911_derivation/README.md]]) | S07 closure | EXT-3 |
| B2 | N6 (waveform-model decyzja) jest blokująca dla N1; bez decyzji nie wiadomo czy liczyć perturbative SPA czy próbować NR | autor decision | EXT-5 problem 5 |
| B3 | Cykl `research/op-LIGO-3G-deviation/` jeszcze NIE ISTNIEJE — T01 NEEDS zakłada że zostanie utworzony | autor + open_bridges declaration | nowy cykl |
| B4 | Cykl `research/op-ppE-mapping/` jeszcze NIE ISTNIEJE | autor + open_bridges declaration | nowy cykl |

## Closed needs

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| N1 | Mapowanie c_n_TGP − c_n_GR → modyfikacja fazy waveformu (SPA) | [[../../research/op-ppE-mapping/Phase1_results.md]] §6 (sympy LOCK 14/14, β_ppE^TGP^(b=-1) = -5/64 central, OOM window [5.5·10⁻², 1.2·10⁻¹]) | 2026-05-07 |
| N2 | Mapowanie analityczne (5/6) U³ → (β_ppE, b_ppE) | [[../../research/op-ppE-mapping/Phase1_results.md]] §6.2 + [[../../research/op-ppE-mapping/Phase3_paper_ready.md]] §1 | 2026-05-07 |
| N3 | Parameter degeneracy (Fisher matrix) | [[../../research/op-LIGO-3G-deviation/Phase2_results.md]] §4 (degeneracy_factor=5; Yagi-Yunes 2016) | 2026-05-07 |
| N4 | SNR thresholds dla detekcji 5σ ET-D/CE | [[../../research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] §2 (calibrated tabela: CE single ~26× margin, ET-D stack 10 BBH, LIGO-O5 stack ~250 BBH) | 2026-05-07 |
| N5 | Wyższe rzędy: U⁴, U⁵, U⁶ deviations | [[../../research/op-ppE-mapping/Phase1_results.md]] §7 (multi-coefficient pattern z 4 PN orders) + M911-P2 wpis w rejestrze | 2026-05-07 |
| N7 | Reżim ważności | [[../../research/op-LIGO-3G-deviation/Phase1_waveform_setup.md]] §1 (cutoff f_max = f_ISCO; perturbative inspiral OK) | 2026-05-07 |
| N8 | Konkretny wpis falsifier w `PREDICTIONS_REGISTRY.md` z ID `M911-P1` (LIVE status) | [[../../PREDICTIONS_REGISTRY.md]] Sector 2b, M911-P1 + M911-P2 + M911-P3 wpisy + Roadmap update 2027+/2035+ | 2026-05-07 |
| N11 | Spójność strong-field z ghost-free basenem ψ ∈ (0, 4/3) | [[../../research/op-LIGO-3G-deviation/Phase1_waveform_setup.md]] §1 (cutoff f_ISCO trzyma U ≪ 1) | 2026-05-07 |
| Q4 | Czy U³ deviation propaguje 1:1 do fazy z O(1) prefactor? | TAK — G_SPA = 1 (central, OOM 0.7-1.5); zob. [[../../research/op-ppE-mapping/Phase1_results.md]] §6.3 | 2026-05-07 |

## Mapowanie need → ścieżka domknięcia (z README)

| Path (README) | Co Path zamyka | Need IDs | Stan w audit-folderze (preview) |
|---------------|----------------|----------|----------------------------------|
| **A — explicit Δh(f) calc** (cykl `op-LIGO-3G-deviation`) | N1, N3, N4, N5, N7, N11 | preview: [[SENSITIVITY_BACK_OF_ENVELOPE.md]] (back-of-envelope, OOM) |
| **B — mapowanie ppE** (cykl `op-ppE-mapping`) | N2 | preview: [[PPN_TO_PPE_MAPPING.md]] (analytic dictionary draft) |
| **C — falsifier statement w `PREDICTIONS_REGISTRY.md`** | N8 | gotowy do adopcji: [[FALSIFIER_STATEMENT_DRAFT.md]] |
| **D — peer-review submit** | N9 | autorska decyzja, brak preview |
| **meta — S07 derywacja** | B1, N10 | [[../S07_M911_derivation/README.md]] |

## Kryterium domknięcia T01

T01 może zostać oznaczony jako **CLOSED-EXECUTED** kiedy spełnione
są jednocześnie:

1. **(C)** Wpis falsifier w `PREDICTIONS_REGISTRY.md` z konkretnym
   thresholdem (β_ppE^(3PN)_TGP_predicted) i mapowaniem na ET/CE
   sensitivity. **Lekkie** — jeden commit, ~30 linii, draft gotowy.
2. **(B)** Cykl `op-ppE-mapping/` Phase 1 zamknięty z:
   - explicit (β_ppE_TGP, b_ppE_TGP = +1) liczbowo,
   - argumentem że to jest *specyficzne dla M9.1''* (a nie generic 3PN).
3. **(A)** Cykl `op-LIGO-3G-deviation/` Phase 2 zamknięty z:
   - Fisher analysis dla ET-D + CE design sensitivities,
   - SNR_threshold dla detekcji 5σ przy ≥ 1 fizycznym scenariuszu BBH/BNS,
   - falsifier-clause: "Jeśli ET/CE przy SNR_network ≥ X i M_chirp ∈ [Y₁, Y₂]
     nie widzą β_ppE^(3PN) > β_TGP_predicted z 5σ, M9.1'' upada."

T01 może być oznaczony jako **CLOSED-PARTIAL (registered)** kiedy
samo (C) jest spełnione (najlżejszy ruch — nie wymaga cyklu badawczego).

## Polityka bezpieczeństwa

- T01 jest klasy **review-only** (`may_edit_core: false`).
- **Edycje w core/** (sek08c, sek08a, TGP_FOUNDATIONS) są **zabronione**
  z poziomu T01 — należą do S07 ścieżki domknięcia.
- **Edycja `PREDICTIONS_REGISTRY.md`** wymaga **autorskiej decyzji**;
  T01 dostarcza *gotowy draft tekstu* w [[FALSIFIER_STATEMENT_DRAFT.md]],
  który autor może wkleić bez modyfikacji.
- T01 **może wskazywać** open-bridges (`op-LIGO-3G-deviation`,
  `op-ppE-mapping`) i wymagać ich utworzenia, ale **nie tworzy ich
  sam** (tworzenie cykli `research/op-…` należy do autora).

## Cross-references

- [[README.md]] — diagnoza T01, ścieżki A/B/C/D
- [[FALSIFIER_STATEMENT_DRAFT.md]] — draft wpisu do `PREDICTIONS_REGISTRY.md` (Path C)
- [[PPN_TO_PPE_MAPPING.md]] — preview Path B (ppE dictionary)
- [[SENSITIVITY_BACK_OF_ENVELOPE.md]] — preview Path A (OOM detectability)
- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5 — recenzja źródłowa
- [[../S07_M911_derivation/README.md]] — meta-bloker B1
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] — pre-existing analytical c_n derivation
- [[../../research/op-newton-momentum/M9_1_pp_setup.md]] — M9.1'' setup
- [[../../PREDICTIONS_REGISTRY.md]] — target dla wpisu falsifier (Path C)
- [[../../TGP_FOUNDATIONS.md]] § 3 — narrative source predykcji
