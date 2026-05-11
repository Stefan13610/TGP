---
title: "Sensitivity OOM (Path A): detektywność (5/6) U³ — Phase 1.5 G_SPA=48 LOCK + GWTC-3 5.02σ FALSIFIED"
date: 2026-05-07 (original) / 2026-05-10 (v2 β=-15/4 LOCK update)
parent: "[[README.md]]"
type: analytical-preview-superseded-by-Phase1.5-and-GWTC3-RERUN
tgp_owner: audyt/T01_LIGO3G_falsifier
revision_history:
  - v1 (2026-05-07): preview Path A — Phase 1 OOM heuristic baseline (β=-5/64)
  - v2 (2026-05-07 sesja C-B-A-D): §4 patched z Phase 1 LOCK (TGP/bound 0.78–780 ratios)
  - v3 (2026-05-09): Phase 1.5 G_SPA=48 LOCK + GWTC-3 RE-RUN 5.02σ FALSIFIED-OBSERVATIONAL
  - v4 (2026-05-10): native-first reframe + thresholds recomputed z β=-15/4 prior;
                     aligns with [[ADDENDUM_2026-05-10_native_observables_first.md]]
tags:
  - sensitivity
  - back-of-envelope
  - ET-D
  - Cosmic-Explorer
  - LIGO-O5
  - LIGO-O3
  - SNR
  - Fisher
  - 2PN-phase
  - T01
  - EXT-5
  - G_SPA-locked
  - sympy-LOCK
  - Phase1.5
  - GWTC3-falsified
  - native-observables-first
related:
  - "[[README.md]]"
  - "[[ADDENDUM_2026-05-10_native_observables_first.md]] (parent native-first methodology)"
  - "[[NEEDS.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
  - "[[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (β=-15/4 sympy LOCK source)"
  - "[[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] (5.02σ FALSIFIED source)"
  - "[[../../research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] (Phase 3 Fisher LOCK)"
  - "[[../../meta/PPN_AS_PROJECTION.md]] (binding methodology)"
---

# Sensitivity back-of-envelope — Path A (post-Phase 1.5 LOCK + GWTC-3 RE-RUN)

## ⚠ STATUS UPDATE 2026-05-10 — Phase 1.5 + GWTC-3 RE-RUN incorporated

> **Two cascading findings invalidate the v1-v2 "preview" framing:**
>
> 1. **Phase 1.5 (2026-05-09):** G_SPA = 48 sympy-exact, NIE ≈ 1 jak Phase 1
>    heuristic zakładał. Source: [[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]].
>    → β_ppE^TGP^(b=-1) = -15/4 ≈ -3.75 (factor 48× larger than v2 assumed).
>
> 2. **GWTC-3 RE-RUN (2026-05-09):** TGP M9.1'' specific (4-3ψ)/ψ ansatz
>    **FALSIFIED at 5.02σ** by GWTC-3 combined posterior (BF = 3.5·10⁻⁶,
>    log10 BF = -5.45). The "borderline LIGO-O3 / decisive ET-CE 2035+"
>    framing is OBSOLETE — falsification has already occurred.
>    Source: [[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]].
>
> **What this file documents now (v4 2026-05-10):**
> - §0: post-Phase 1.5 + post-RE-RUN summary table (PRIMARY)
> - §1-§3: literature ppE Fisher bounds (still valid; depend only on detector
>   sensitivity, not on TGP β value)
> - §4: detection table RECOMPUTED with β=-15/4 prior — shows that in CURRENT
>   LIGO-O3 + GWTC-3 era, TGP/bound ratio ~5σ already (post-detection regime,
>   confirmed via TIGER framework Bayes RE-RUN)
> - §5-§6: caveats + falsification window (collapsed to "already past")
>
> **Native-first reframe (binding 2026-05-10+):** native L1 quantities
> (c_3 = -5/6, Δe_2 = -4/3, Δα_4 = -40) sympy-LOCKED; L2 chart projection
> (β_ppE, b_ppE, G_SPA) sympy-LOCKED; L3 falsifier (GWTC-3 BF, ET-D + CE
> thresholds) gives empirical verdict. M9.1'' specific f(ψ) FALSIFIED;
> recovery via emergent-metric Phase 4 zero-β region OR alternative S07
> f(ψ) reset remains viable.

## §0 — Post-Phase 1.5 + post-RE-RUN summary (PRIMARY 2026-05-10)

| Quantity | Phase 1 (v1, 2026-05-07) | Phase 1.5 (v3, 2026-05-09) | Comment |
|---|---|---|---|
| G_SPA | ≈ 1 (SYC 2013 heuristic, OBSOLETE) | **48** sympy-exact | factor 48× correction |
| β_ppE^TGP^(b=-1) at η=1/4 | -5/64 ≈ -0.078 | **-15/4 ≈ -3.75** | factor 48× larger |
| OOM window |β| | [0.055, 0.12] | **[2.81, 4.69]** | test-p ± 25% η-correction |
| LIGO-O3 status | borderline | **5.02σ FALSIFIED-OBSERVATIONAL** | GWTC-3 RE-RUN (TIGER) |
| LIGO-O5 prediction | "first decisive 2027+" | already falsified | obsolete framing |
| ET-D + CE prediction | "1000σ at single loud BBH" | already falsified | obsolete framing |
| Window of testability | "2027–2035" | already past (2019-2020 GWTC-3 closed it) | falsification occurred in past data |

**Native parameter audit (post-Phase 1.5):**

```
L1 native source (sympy-LOCKED, op-newton-momentum):
  c_3 (g_tt at U³)     = -5/6     [α=2 vacuum + M9.1'' f(ψ) FORCED]
  c_4 (g_tt at U⁴)     = -23/12   [structural]
  c_5 (g_tt at U⁵)     = +337/72  [structural]

L1 native dynamics (sympy-LOCKED, op-ppE-mapping Phase 1.5):
  Δe_2 (orbital)       = -4/3     [test-particle exact, η-correction ±25%]
  Δα_4 (TaylorF2)      = -40      [test-particle exact, sympy-exact rational]

L2 chart projection (sympy-LOCKED, op-ppE-mapping Phase 1.5):
  G_SPA                = 48       [test-particle exact, NIE 1 jak SYC 2013 heuristic]
  β_ppE^TGP^(b=-1)     = -15/4    [η=1/4, test-p ± 25%]
  b_ppE^TGP            = -1       [2PN-phase, U³ metric → phase]

L3 falsifier (empirical, op-GWTC3-reanalysis Phase 2 RE-RUN):
  GWTC-3 BF_TGP/GR     = 3.5·10⁻⁶ [combined posterior, TIGER framework]
  GWTC-3 σ-level       = 5.02σ    [FALSIFIED-OBSERVATIONAL]
  ET-D + CE windows    = N/A      [falsification already occurred in past data]

Falsifier scope: M9.1'' SPECIFIC ansatz f(ψ)=(4-3ψ)/ψ falsified.
Recovery: S07 alternative f(ψ) reset OR emergent-metric Phase 4
          zero-β region in (a_n, ξ_n, b_n, c_0, κ_σ) space.
```

## §0.1 — Cel pliku (HISTORICAL preview-era)

> **Original cel (2026-05-07).** Order-of-magnitude (OOM) szacowanie
> detektywności deviation (5/6) U³. Plik był previewem dla cyklu
> `research/op-LIGO-3G-deviation/`. Cykl został wykonany w sesji
> 2026-05-07 (Phase 0+1+2+3 Fisher matrix, calibrated thresholds).
> Phase 1.5 (2026-05-09) + GWTC-3 RE-RUN (2026-05-09) razem zamknęły
> wszystkie open assumptions A1-A6 i dały empirical verdict.
>
> **Status v4 (2026-05-10):** ANALYTICAL DOKUMENT (post-Phase 1.5 LOCK
> + post-RE-RUN). Liczby są teraz sympy-exact lub literature-locked.
> Sekcje §1-§3 zachowują strukturę OOM jako pedagogical reference;
> §4 recomputed z β=-15/4 prior; §0 (post-LOCK summary) jest primary
> post-2026-05-10.

## §1 — Reżim ważności U dla inspiralu BBH

### 1.1 Mapping orbital parameter U ↔ frequency f

Dla compact binary o całkowitej masie M_tot na okrągłej orbicie
o promieniu r:

```
U = G·M_tot / (r·c²)         (Newtonian potential at separation r)
v² ≈ U·c²                    (Keplerian: v² ≈ GM/r)
v/c ≈ √U
```

Dla GW częstotliwości f_GW = 2 · f_orb (quadrupole leading), Kepler:

```
f_GW = (1/π) · (G·M_tot)^(1/2) / r^(3/2)         × (multiplikator z factor 2)
       
       (πM·f_GW)^(2/3) ≈ U
                      ⇒  U(f) = (πM·f)^(2/3)         [M = G·M_tot/c³ in time units]
```

(używam tu `M` = chirp mass scale w jednostkach czasu, M = G·M_tot/c³).

### 1.2 Skala U dla typowych BBH/BNS w paśmie ET/CE

| Scenariusz | M_tot | f_GW (band) | M = GM_tot/c³ | U(f) leading | (5/6) U³ |
|------------|-------|-------------|----------------|--------------|----------|
| BBH equal-mass loud | 30 M_⊙ | 10–500 Hz | 1.5 · 10⁻⁴ s | 0.005 → 0.07 | 1·10⁻⁷ → 3·10⁻⁴ |
| BBH equal-mass heavy | 60 M_⊙ | 5–250 Hz | 3.0 · 10⁻⁴ s | 0.005 → 0.07 | 1·10⁻⁷ → 3·10⁻⁴ |
| BBH equal-mass GW150914 | 65 M_⊙ | 30–250 Hz, ISCO ~80 Hz | 3.2 · 10⁻⁴ s | 0.02 → 0.07 | 7·10⁻⁶ → 3·10⁻⁴ |
| BNS GW170817 | 2.74 M_⊙ | 30–1500 Hz | 1.4 · 10⁻⁵ s | 0.001 → 0.026 | 8·10⁻¹⁰ → 1.5·10⁻⁵ |
| EMRI (LISA, dla porównania) | 10⁶ M_⊙ + 10 M_⊙ | 10⁻⁴–10⁻¹ Hz | 5 · 10⁰ s | 0.1 → ~0.3 | 8·10⁻⁴ → 2·10⁻² |

**Punkty kluczowe:**
- Dla **loud BBH** w ET/CE band, U sięga ~0.07 tuż przed mergerem
  (przy ISCO; r ~ 6GM/c²).
- (5/6) U³ ≈ 3 · 10⁻⁴ przy peak (U=0.07), ale **kumuluje się
  w fazie** przez wiele cykli (N_cycles ≈ 10²–10³ w ET-D band).
- Faza waveformu przy 3PN deviation kumuluje:
  ```
  Δφ_total ≈ (5/6) · ⟨U³⟩ · N_cycles · O(1)
           ≈ 3·10⁻⁴ · 100 · 1 = 3·10⁻²  rad   [GW150914-like loud]
  ```
- 5σ detekcja wymaga Δφ_total > σ_φ (phase noise w detektorze).
  Dla SNR=100 typowo σ_φ ~ 1/SNR ~ 10⁻². Dla SNR=500 (ET-D loud
  events) σ_φ ~ 2·10⁻³.

**Wniosek pierwszego rzutu:** **(5/6) U³ jest detektywna w ET-D
dla loud BBH events przy SNR ≳ 100**, jeśli faza kumuluje O(N_cycles).
To zgodne z M9_1_pp_P1_results.md §4.3 OOM (waveform sensitivity range
10⁻³–10⁻⁵).

## §2 — SNR targets w detektorach 3G

### 2.1 Single-event SNR (literatura)

| Detektor | Era | Strain h_n(f) @ 100 Hz | SNR GW150914 (M=65 M_⊙, 410 Mpc) | SNR loud BBH (M=30 M_⊙, 200 Mpc) | Source |
|----------|-----|------------------------|-----------------------------------|-----------------------------------|--------|
| LIGO-O3 | 2019–2020 | ~5·10⁻²³ /√Hz | ~24 | ~30 | LIGO O3 noise budget |
| LIGO-O5 | 2027+ | ~3·10⁻²³ /√Hz | ~40 | ~80 | A+ design |
| LIGO-O5 ext | 2030+ | ~2·10⁻²³ /√Hz | ~60 | ~120 | A# upgrade |
| Einstein Telescope D | ~2035 | ~5·10⁻²⁴ /√Hz | ~250 | ~500 | Maggiore et al. 2020 |
| Cosmic Explorer | ~2035+ | ~2·10⁻²⁵ /√Hz @ low f | ~600 | ~1000 | Reitze et al. 2019 |
| ET + CE network | ~2035+ | combined | ~1000 | ~1500 | 1912.02622 §3 |

(strain numbers są round; precise depends on f-dependent
sensitivity curve.)

### 2.2 Stacked SNR przez N events

Dla **N independent events** z `np.sqrt(N)` scaling SNR-stack:

| Network | N events | SNR_stack | β_ppE_5σ_bound |
|---------|----------|-----------|-----------------|
| LIGO-O5 | 100 BBH/yr | ~800 | ~3·10⁻³ |
| ET + CE | ~10⁵ BBH/yr (!!) | ~3·10⁵ | ~10⁻⁵ |

**ET-D + CE detection rate** projektowany na ~10⁵ BBH/yr
(Maggiore 2020 §2; Reitze 2019). Po ~5 latach operacji stacked SNR
~5·10⁵ — pozwala constraint β_ppE^(3PN) na poziomie ~10⁻⁵–10⁻⁶.

## §3 — ppE Fisher bounds z literatury

### 3.1 Single-event β_ppE^(N PN) bound

Z Chamberlain & Yunes 2017 (arXiv:1704.08268, Tabela II–IV):

| PN order | β_ppE_5σ (LIGO-O5) | β_ppE_5σ (ET-D) | β_ppE_5σ (ET+CE) |
|----------|--------------------|------------------|--------------------|
| -1PN (b=-7) | ~10⁻⁵ | ~10⁻⁷ | ~3·10⁻⁸ |
| 0PN (b=-5) | ~10⁻⁴ | ~10⁻⁶ | ~3·10⁻⁷ |
| 1PN (b=-3) | ~10⁻² | ~10⁻⁴ | ~3·10⁻⁵ |
| 2PN (b=-1) | ~10⁻¹ | ~10⁻³ | ~3·10⁻⁴ |
| **3PN (b=+1)** | **~3·10⁻¹** | **~3·10⁻³** | **~10⁻³** |

(scenariusz: GW150914-like loud BBH, 5σ confidence; numbers ±factor 2).

### 3.2 Stacked bounds (~50–500 BBH events)

Z Maggiore et al. 2020 §3 i Yagi–Yunes 2016 (1602.04674):

| Detektor | N events | β_ppE_3PN bound |
|----------|----------|------------------|
| LIGO-O5 | 100 BBH | ~3·10⁻² |
| ET-D | 1000 BBH | ~3·10⁻⁴ |
| ET+CE | 5000 BBH | ~10⁻⁴ |
| ET+CE | 10⁵ BBH (5yr) | ~10⁻⁵ |

## §4 — Skonfrontowanie: TGP M9.1'' β_ppE^TGP vs detector bounds (RECOMPUTED v4 2026-05-10)

### 4.1 Skala β_ppE^TGP (LOCKED Phase 1.5, 2026-05-09)

**Phase 1.5 LOCK:** [[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] §4-§5
(sympy LOCK 5/5 PASS + 4-level verification: hand-calc + numerical sanity at
U=0.1 + alternative SPA orthogonal route → all converge to identical β=-15/4).

```
β_ppE^TGP^(b=-1) = -(3/(128·η)) · Δα_3_metric · G_SPA   [LOCKED form]

Equal-mass (η=1/4), G_SPA = 48 (sympy-exact test-particle):
   β_ppE^TGP^(b=-1) = -(3/(128·1/4)) · (-5/6) · 48 = -(3/32) · 40 = -15/4 ≈ -3.75

Test-particle ± 25% η-correction window:
   |β_ppE^TGP^(b=-1)| ∈ [2.81, 4.69]
```

**Phase 1 vs Phase 1.5:** Phase 1 szacował β ≈ -7.8·10⁻² (central z G_SPA ≈ 1
heuristic). Phase 1.5 sympy-LOCKED β = -3.75 (factor **48× larger**).
Phase 1 heuristic G_SPA ≈ 1 was applied OUTSIDE its regime of validity
(SYC 2013 small-perturbation framework; TGP M9.1'' jest structural O(1)
modification → SPA chain amplification G_SPA = 48).

### 4.2 Tabela detekcji (RECOMPUTED z β=-15/4 prior, v4 2026-05-10)

**Konwencja:** wszystkie wartości β_ppE_5σ_bound poniżej są dla
**2PN-phase, b_ppE = −1** (zob. [[CONVENTION_DECISION.md]]; bounds wzięte
z linii 2PN tabeli §3.1; Phase 3 Fisher LOCK
[[../../research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] §2
confirms detector thresholds independent of theory prediction).

| Detektor / scenariusz | β_ppE_5σ_bound^(b=−1) | β_ppE^TGP (LOCKED v3) | TGP/bound | Detekcja TGP? (v4 verdict) |
|-----------------------|------------------------|------------------------|-----------|------------------------------|
| **LIGO-O3 / GWTC-3 (TIGER reanalysis)** | β_5σ ≈ 0.78 (single-coef) | **3.75** | **~5.02σ** | **FALSIFIED-OBSERVATIONAL** ([[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]) |
| LIGO-O3 generic ToGR (multi-coef marginalized) | ~10⁻¹ | 3.75 | ~37 | would be >30σ but generic analysis ~50× weaker than dedicated TIGER |
| LIGO-O5 single | ~3·10⁻² | 3.75 | ~125 | (academic) >100σ if not already falsified |
| LIGO-O5 stack 100 BBH | ~3·10⁻³ | 3.75 | ~1250 | (academic) >1000σ |
| ET-D single (GW150914-like loud) | ~10⁻³ | 3.75 | ~3750 | (academic) >3000σ |
| ET-D stack 100 BBH | ~10⁻⁴ | 3.75 | ~37500 | (academic) >30000σ |
| CE single (loud BBH) | ~3·10⁻⁴ | 3.75 | ~12500 | (academic) >10000σ |
| ET+CE stack 5000 BBH | ~10⁻⁴ | 3.75 | ~37500 | (academic) >30000σ |

**Wszystkie wartości "(academic)" oznaczają: gdyby M9.1'' przetrwał GWTC-3
falsification, byłby trywialnie wykryty w późniejszych detektorach.
Ponieważ falsification już zaszedł w GWTC-3, te kolumny są historyczne.**

### 4.3 GWTC-3 RE-RUN verdict (post-Phase 1.5 corrected β, 2026-05-09)

Per [[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]:

```
GWTC-3 combined posterior (TIGER framework, ~90 BBH events):
  σ_β at b=-1 (2PN-phase)     ≈ 0.78        (single-coef Bayes prior)
  β_ppE^TGP (Phase 1.5 LOCK)  = -3.75       (η=1/4, ±25% test-p approx)
  TGP / σ_β                    = 4.81        (Z-score)
  Bayes factor BF_TGP/GR       = 3.5·10⁻⁶   (combined)
  log10 BF                     = -5.45      ("OVERWHELMING GR preference")
  σ-level vs combined          = 5.02σ       (FALSIFIED-OBSERVATIONAL)

Caveat: falsification applies to SPECIFIC M9.1'' ansatz f(ψ) = (4-3ψ)/ψ.
        Alternative f(ψ) structures (S07 reset) + emergent-metric Phase 4
        recovery framework (zero-β region) remain viable.
```

### 4.4 Główny insight (post-Phase 1.5 + post-RE-RUN)

1. **Falsification window collapsed.** Preview-era language "borderline LIGO-O3"
   był oparty na Phase 1 heuristic G_SPA ≈ 1. Phase 1.5 sympy-LOCK G_SPA = 48
   gives factor 48× larger β_TGP, czyniąc M9.1'' specific ansatz **już
   sfalsyfikowanym 5σ** w current LIGO-O3 / GWTC-3 data.

2. **"Window of testability 2027–2035" nieaktualne.** Falsifikacja już zaszła
   w 2019-2020 (LIGO O3 obs run, GWTC-3 catalog). RE-RUN 2026-05-09 z corrected
   β prior dokumentuje to retrospectively.

3. **Multi-coefficient pattern misleading (Phase 1 heuristic).** Phase 1 ratios
   {-23/10, -38/23, +337/228} are INCORRECT (Phase 1.5 alternative SPA gives
   β_3PN/β_2PN = -11161/504 ≈ -22.14, factor ~10× off). M911-P2 entry of
   PREDICTIONS_REGISTRY → WITHDRAWN-needs-rederivation.

4. **Recovery via emergent-metric Phase 4.** [[../../research/op-emergent-metric-from-interaction-2026-05-09]]
   Phase 4 STRUCTURAL DERIVED zero-β region in (a_n, ξ_n, b_n, c_0, κ_σ)
   parameter space — Path 1 (c_0 = 0 substrate) i Path 2 (κ_σ = 4/3 canonical
   coupling) viable post-falsification. Falsification dotyczy *specific
   f(ψ)*, NIE TGP framework całość.

5. **Methodological finding (independent of TGP fate).** SPA chain analysis
   for STRUCTURAL-modification theories (O(1) coupling, NOT small-perturbation)
   gives G_SPA ≠ 1 (factor 48 in this case). SYC 2013 "G_SPA ≈ 1" common
   assumption FAILS for such theories — notable methodological caveat for
   future ppE catalog work (Yagi-Yunes-Pretorius 2016 review etc.).

## §5 — Reżim ważności i caveat-e

### 5.1 Co OOM **NIE** uwzględnia (zostaw cyklowi `op-LIGO-3G-deviation/`)

| Effect | Magnitude | Treatment |
|--------|-----------|-----------|
| Spin precession | mała przy aligned spins, do x10× degradation przy precessing | Phase 2 Fisher z spin |
| Chirp-mass × spin × β kowariancja | typowo 30–80% degradacja β bound | Phase 2 Fisher z full param |
| Inclination averaging | factor √2 SNR fluctuation | sky-average |
| 3.5PN baseline waveform | TGP modyfikuje 3PN, nie 3.5PN | careful modeling |
| Higher-order TGP coefficients (4PN, 5PN) | własne Fisher detekcja | multi-coefficient ppE |
| Tidal effects (BNS) | Λ_tidal absorption może maskować | exclude BNS |

### 5.2 Co OOM **DOBRZE** szacuje

- Order of magnitude detektywności pojedynczego loud BBH ✓
- Hierarchy LIGO-O3 → O5 → ET → CE ✓
- General feasibility: M9.1'' P1 jest detektywna w 3G era ✓

### 5.3 Caveat: **konwencja PN counting**

Wszystkie liczby powyżej zakładają konwencję "phase 3PN = b_ppE = +1"
(zgodnie z Cutler–Flanagan / Yunes–Pretorius). Jeśli cykl
`op-ppE-mapping` zatwierdzi inną konwencję (energy 3PN ↔ phase 2PN
↔ b = -1), wtedy:

- bounds są **stronger** o ~10× (zob. tabela §3.1: b=-1 vs b=+1
  bounds różnią się o ~factor 10).
- TGP detection jest *jeszcze* clearer.

To jest **dobry argument dla TGP**: w ANY konwencji, deviation
jest detektywna.

## §6 — Konkluzja v4 (post-Phase 1.5 LOCK + post-GWTC-3 RE-RUN)

### 6.1 Recomputed answers do questions z [[NEEDS.md]]

| Pytanie | Odpowiedź v1-v2 (Phase 1) | Odpowiedź v4 (Phase 1.5 + RE-RUN) |
|---------|---------------------------|------------------------------------|
| Q1 — najsłabszy ET sygnał dla 5σ detekcji? | M_BBH ≥ 30 M_⊙, d ≤ 1 Gpc, SNR ≥ 100 | **N/A** — falsification already occurred in GWTC-3; ET-D not needed |
| Q3 — czy O3/O4 już *constraints*? | NIE — TGP w basenie ~10⁻¹ | **TAK** — TGP-specific TIGER reanalysis 5.02σ FALSIFIED ([[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]) |
| (komplement) — kiedy decyzyjne? | LIGO-O5 stack ~2027–2030 | **już zdecydowane (2019-2020 GWTC-3 catalog)** |

### 6.2 v4 verdict (post-Phase 1.5 + post-RE-RUN)

**(5/6) U³ deviation z M9.1'' specific f(ψ)=(4-3ψ)/ψ ansatz** jest:

1. **Empirycznie sfalsyfikowane 5.02σ** w GWTC-3 (combined ~90 BBH events,
   TIGER framework Bayes inference z corrected β=-15/4 prior).
2. **Methodologically robust** — sympy-LOCK 5/5 PASS at multiple verification
   levels (Phase 1 + Phase 1.5 + 4-level cross-check + alternative SPA route).
3. **Native parameter audit clean** — zero free parameters; falsification dotyczy
   *specific ansatz form*, NIE framework jako całość.

### 6.3 Recovery paths (post-falsification, viable)

| Path | Description | Status |
|------|-------------|--------|
| **S07 alternative f(ψ)** | Reset M9.1'' z innym hyperbolic / power ansatz; preserve native L1 c_n structure | OPEN — S07 audit cycle |
| **Emergent-metric Phase 4** | zero-β region in (a_n, ξ_n, b_n, c_0, κ_σ) substrate parameter space | STRUCTURAL DERIVED ([[../../research/op-emergent-metric-from-interaction-2026-05-09]]) |
| **Path 1: c_0 = 0** | substrate scaling vanishes → β-independent g_eff[Φ] | viable; in Phase 4 zero-β |
| **Path 2: κ_σ = 4/3 canonical** | substrate-coupling preserves PPN γ=β=1 + zero β at 2PN-phase | viable; in Phase 4 zero-β |

**Cykl `research/op-LIGO-3G-deviation/` Phase 3** (already executed,
2026-05-07) provides Fisher LOCK on detector thresholds — those are
**theory-independent** i remain valid as bounds on β_ppE^(b=-1) regardless
of which f(ψ) ansatz is tested.

### 6.4 Native-first scope clarification (binding 2026-05-10+)

Per [[../../meta/PPN_AS_PROJECTION.md]] §3.1, falsification result reads:

| Layer | Status |
|---|---|
| **L1 native** (Φ-EOM derivation, c_n coefficients) | **PRESERVED** — α=2 vacuum + (4-3ψ)/ψ produces sympy-LOCKED c_n |
| **L2 chart** (β_ppE projection) | **VALIDATED** — sympy-LOCK 5/5 + 4-level confirms Phase 1.5 derivation |
| **L3 falsifier** (GWTC-3 BF / σ-level) | **FALSIFIED** — 5.02σ on M9.1'' SPECIFIC ansatz only |

Falsification of L3 nie sfalsyfikowuje native L1 framework; sfalsyfikowuje
*specific functional form* mapping substrate Φ → effective metric. Recovery
via alternative L1 source structure (S07 reset) lub via emergent-metric Phase 4
zero-β substrate-parameter region.

## Cross-references

- [[README.md]] — diagnoza T01, ścieżki A/B/C/D
- [[NEEDS.md]] — N3, N4, N7 (Fisher, SNR thresholds, validity regime)
- [[FALSIFIER_STATEMENT_DRAFT.md]] §1 tabela — zaktualizować numbers
  z OOM na precise values po zamknięciu Phase 2
- [[PPN_TO_PPE_MAPPING.md]] — PN convention decision
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] §4.3 — concordant OOM (independent estimation)
- [[../../PREDICTIONS_REGISTRY.md]] BH5, ε.1 — orthogonal channel predictions

## Bibliografia

- LIGO Scientific Collaboration et al., Phys. Rev. D **103**, 122002 (2021), arXiv:2010.14529 (GWTC-2 ppE constraints)
- Maggiore et al., JCAP **03**, 050 (2020), arXiv:1912.02622 (ET science case)
- Reitze et al., Bull. Am. Astron. Soc. **51**, 035 (2019) (CE)
- Chamberlain & Yunes, Phys. Rev. D **96**, 084039 (2017), arXiv:1704.08268 (3G implications)
- Yagi & Yunes, arXiv:1602.04674 (review)
- Cutler & Flanagan, Phys. Rev. D **49**, 2658 (1994) (SPA)
- Will, *Living Rev. Relativ.* **17**, 4 (2014) (PPN review)
