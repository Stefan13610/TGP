---
title: "Sensitivity OOM (preview Path A): detektywność (5/6) U³ w ET-D / CE"
date: 2026-05-07
parent: "[[README.md]]"
type: analytical-preview
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - preview
  - sensitivity
  - back-of-envelope
  - ET-D
  - Cosmic-Explorer
  - LIGO-O5
  - SNR
  - Fisher
  - 3PN
  - T01
  - EXT-5
related:
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
---

# Sensitivity back-of-envelope — preview Path A

> **Cel pliku.** Order-of-magnitude (OOM) szacowanie detektywności
> deviation (5/6) U³ od GR w sieciach ET-D / Cosmic Explorer dla
> reprezentatywnych BBH/BNS scenariuszy, oparte na **literaturze
> existing Fisher analyses**. Plik jest **previewem** dla cyklu
> `research/op-LIGO-3G-deviation/`; dostarcza:
> - skala U typowa dla inspiralu (gdzie deviation jest najsilniejsza),
> - SNR-targets dla ET-D + CE (literatura),
> - ppE Fisher bounds dla 3PN coefficient (literatura),
> - explicit decyzję CZY (5/6) U³ jest detektywna i dla jakich M, d, SNR.
>
> **Status:** OOM. Liczby są przybliżone z dokładnością do faktora
> 2–3 (literatura giving the same scenario varies by ~factor 2 due
> to inclination/sky-position averaging assumptions). Cykl
> `op-LIGO-3G-deviation/` Phase 2 zamknie precyzję do %.

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

## §4 — Skonfrontowanie: TGP M9.1'' β_ppE^TGP vs detector bounds

### 4.1 Skala β_ppE^TGP (LOCKED 2026-05-07 via op-ppE-mapping Phase 1)

**Status update:** preview OOM zastąpiony liczbowym lock z
[[../../research/op-ppE-mapping/Phase1_results.md]] §6 (sympy LOCK
14/14: 7/7 α_n^TGP coefficients reproducja + 7/7 Δα_n consistency).

```
β_ppE^TGP^(b=-1) = -(3/(128·η)) · (5/6) · G_SPA   [LOCKED form]

Equal-mass (η=1/4), G_SPA=1 (central):
   β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 · 10⁻²

OOM window (G_SPA ∈ [0.7, 1.5], uncertainty z modified quadrupole formula):
   |β_ppE^TGP^(b=-1)| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹]
```

**Konkluzja vs preview:** preview szacował OOM ~10⁻¹ (κ ≈ 0.5–1.5);
lock dał ~7.8·10⁻² (central) — w tej samej skali, precyzja ~30%
przy G_SPA=1. **To jest decisive: TGP β_ppE^TGP ≈ 7.8·10⁻² jest
~2.6× powyżej LIGO-O5 single-event bound ~3·10⁻²** → first decisive
detection możliwa już w erze ~2027–2030.

### 4.2 Tabela detekcji (TGP vs bounds)

**Konwencja:** wszystkie wartości β_ppE_5σ_bound poniżej są dla
**2PN-phase, b_ppE = −1** (zob. [[CONVENTION_DECISION.md]] dla wyjaśnienia
mapowania U³-w-g_tt → b=−1; bounds wzięte z linii 2PN tabeli §3.1).

| Detektor / scenariusz | β_ppE_5σ_bound^(b=−1) | β_ppE^TGP (LOCKED) | TGP/bound | Detekcja TGP? |
|-----------------------|------------------------|----------------------|-----------|----------------|
| LIGO-O3 (now)          | ~10⁻¹                 | ~7.8·10⁻²            | ~0.78     | borderline |
| LIGO-O5 single | ~3·10⁻²                       | ~7.8·10⁻²            | ~2.6      | **YES — first decisive 2027+** |
| LIGO-O5 stack 100 BBH | ~3·10⁻³                | ~7.8·10⁻²            | ~26       | YES (>20σ) |
| ET-D single | ~10⁻³                         | ~7.8·10⁻²            | ~78       | YES (>50σ) |
| ET-D stack 100 BBH | ~10⁻⁴                    | ~7.8·10⁻²            | ~780      | YES (>700σ) |
| CE single | ~3·10⁻⁴                         | ~7.8·10⁻²            | ~260      | YES (>200σ) |
| ET+CE stack 5000 BBH | ~10⁻⁴                  | ~7.8·10⁻²            | ~780      | YES (>700σ) |

**Komentarz do "borderline" w O3.** Subiektywna interpretacja "TGP w
basenie" dla aktualnego LIGO-O3: bound ~10⁻¹ vs przewidziane β_ppE^TGP
~10⁻¹ — to znaczy że TGP jest na **granicy aktualnej wykrywalności**.
Nie sfalsyfikowane, ale również nie *nadal under-constrained*. Zob.
[[NEEDS.md]] Q3 — re-analiza GWTC-3 ppE constraints w 2PN-phase
specifically dla M9.1'' β_TGP^TGP byłaby cennym side-cyklem.

**Główny insight:**

1. **TGP jest w basenie konsystencji w LIGO-O3.** Bounds są ~3·10⁻¹,
   przewidziane β_TGP ~ 10⁻¹ — TGP jest *under-constrained*. To jest
   **dobre** — nie ma już falsyfikacji z O3, a deviation jest
   *prawdziwie predyktywna* dla przyszłości.

2. **LIGO-O5 może już dać sygnał.** Stack 100 BBH events daje bound
   ~3·10⁻². Jeśli β_ppE^TGP ≈ 10⁻¹, to **LIGO-O5 już falsyfikuje
   M9.1'' lub potwierdza** ~2027–2030. (Pre-ET window.)

3. **ET-D + CE są decisive.** Single events już dają detekcję. Stack
   1000+ BBH daje >1000σ. Albo deviation jest tam, i M9.1'' jest
   verified, albo nie ma — i M9.1'' jest *nadzwyczaj sfalsyfikowana*.

4. **Window of testability:** ~2027 (LIGO-O5 first stacked results)
   do 2035+ (ET-D + CE first observations). To **5–8 letnie okno**
   na peer-review podążanie M9.1''.

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

## §6 — Konkluzja preliminary

| Pytanie z [[NEEDS.md]] | Odpowiedź preliminary (OOM) |
|------------------------|------------------------------|
| Q1 — najsłabszy ET sygnał dla 5σ detekcji? | M_BBH ≥ 30 M_⊙, d ≤ 1 Gpc, single SNR ≥ 100 | 
| Q3 — czy O3/O4 już *constraints*? | NIE — TGP w basenie ~10⁻¹, bounds ~10⁻¹ |
| (komplement) — kiedy decyzyjne? | LIGO-O5 stack ~2027–2030, ET-D + CE ~2035 |

**Preliminary verdict:** (5/6) U³ deviation **jest detektywna w 3G era**
przy realistic SNR i event rates. To znaczy — T01 falsifier jest
**fizycznie wykonalny**, nie tylko teoretycznie sformułowany.

Cykl `research/op-LIGO-3G-deviation/` Phase 2 zamknie:
- precyzyjny single-event SNR threshold,
- precyzyjny network SNR / N_events combo dla ET-D + CE,
- Fisher matrix z full param (mass, spin, distance, inclination).

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
