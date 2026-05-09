---
title: "Draft falsifier statement dla PREDICTIONS_REGISTRY (M911-P1: |Δg_tt| = (5/6) U³)"
date: 2026-05-07
parent: "[[README.md]]"
type: registry-draft
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - draft
  - falsifier
  - registry
  - M911
  - 3PN
  - LIGO-3G
  - Einstein-Telescope
  - Cosmic-Explorer
  - ppE
  - T01
  - EXT-5
related:
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]]"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
---

# Draft falsifier statement — `M911-P1` (Path C, najlżejsze domknięcie T01)

> ## ⚠ CRITICAL UPDATE 2026-05-09 — STATUS CHANGE: FALSIFIED-OBSERVATIONAL
>
> **[[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] Phase 1.5
> derived G_SPA = 48 sympy-exact (test-particle, 4-level verified):**
> sympy LOCK 5/5 + independent hand-calculation E_TGP(U³)/m = 49/48 +
> numerical sanity at U=0.1 + alternative SPA derivation (orthogonal
> route) ALL CONFIRM β_ppE^TGP^(b=-1) = **-15/4 ≈ -3.75** at η=1/4
> (test-particle exact, ±25% η-correction estimate).
>
> This is **factor 48× LARGER** than Phase 1's heuristic value of
> -5/64 ≈ -0.078. Phase 1's "G_SPA ≈ 1" assumption from Sampson-Yunes-
> Cornish 2013 was applicable to SMALL-PERTURBATION regimes (BD with
> 1/ω_BD, dCS with ζ_dCS) but **DOES NOT APPLY** to TGP M9.1'' which
> has STRUCTURAL O(1) modifications via hyperbolic f(ψ) = (4-3ψ)/ψ.
>
> **GWTC-3 RE-RUN with corrected β** ([[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]):
> TGP M9.1'' **RULED OUT at 5.02σ** by GWTC-3 combined ~90 BBH posterior
> (BF_TGP/GR = 3.5·10⁻⁶, log10(BF) = -5.45 → "OVERWHELMING GR preference").
>
> **Falsifier verdict:** M9.1'' as currently formulated (specific
> hyperbolic ansatz f(ψ) = (4-3ψ)/ψ) is **observationally excluded**
> at 5σ by current GWTC-3 data in TIGER framework analysis.
>
> **Caveat:** falsification applies to specific (4-3ψ)/ψ form. Alternative
> f(ψ) structures via S07 audit (M9.1'' derivation) remain viable for
> exploration. Sections §1-§4 below describe the OUTDATED falsifier
> statement based on Phase 1 OOM heuristic (β = -5/64); they are
> historically informative but NO LONGER LIVE.
>
> Current operational status:
> - M911-P1: **FALSIFIED-OBSERVATIONAL** (PREDICTIONS_REGISTRY updated 2026-05-09).
> - M911-P2: **WITHDRAWN-NEEDS-REDERIVATION** (Phase 1 ratios incorrect; alternative SPA shows β_3PN/β_2PN = -11161/504 ≠ -23/10).
> - M911-P3: **PARTIAL-FALSIFIED** (2/4 channels invalid; remaining BH5/ε.1 still independent test of any future revised M9.1'').
> - Path D paper draft: requires major revision from "predictive forecast" to "negative result + factor-48 methodological finding".

---

> **Cel pliku.** Dostarczyć **gotowy do wklejenia** tekst wpisu do
> `PREDICTIONS_REGISTRY.md` z explicit warunkiem falsyfikacji predykcji
> M9.1'' P1 (`|Δg_tt| = (5/6) U³`). Bez tego wpisu predykcja jest
> *opisem konsekwencji* (T01 problem #2). Z nim — *kontraktem
> z obserwacjami*.
>
> **Autor TGP może wkleić §1 i §2 do `PREDICTIONS_REGISTRY.md`** bez
> modyfikacji; uzupełnić threshold liczbowy z cyklu `op-LIGO-3G-deviation/`
> Phase 2 kiedy zostanie zamknięty (placeholder `[β_th]` w tekście).
>
> **Polityka:** ten plik jest **draftem**. Jego *wklejenie* do
> `PREDICTIONS_REGISTRY.md` jest poza zakresem T01 (`may_edit_core: false`)
> i wymaga decyzji autora. T01 dostarcza tekst, autor go promuje.
>
> **Status as of 2026-05-09:** sections below describe Phase 1 baseline
> (now SUPERSEDED by Phase 1.5). For current locked falsifier verdict
> see CRITICAL UPDATE block above.

## §1 — Proponowana sekcja w `PREDICTIONS_REGISTRY.md` (do wklejenia)

```markdown
### M911-P1 — |Δg_tt| = (5/6) U³ deviation od GR (3PN)

**Predykcja.** W M9.1'' z α=2 i hiperbolicznej f(ψ)=(4-3ψ)/ψ
rozwinięcie statycznej, sferycznie symetrycznej metryki w
potencjale Newtonowskim U = GM/(rc²) daje:

```
g_tt^TGP / (-c²) = 1 − 2U + 2U²       − (7/3)U³ + (35/12)U⁴ − (91/24)U⁵ + …
g_tt^GR  / (-c²) = 1 − 2U + 2U²       − (3/2)U³ +    1·U⁴   − (5/8)U⁵   + …
                  └─ Newton ─┘└1PN┘ └─ 2PN ─┘└── 3PN ──┘└────── 4PN+ ──────┘
```

(zob. [[research/op-newton-momentum/M9_1_pp_P1_results.md]] §3.2,
sympy LOCK c_n=2..7 5/5 PASS).

Pierwsza strukturalna deviation od GR pojawia się w **3PN** (rząd U³):

```
Δg_tt / (-c²) = − (5/6) U³ + (23/12) U⁴ − (19/6) U⁵ + (337/72) U⁶ + …
```

**Status:** (W) wewnątrz M9.1'' (analitycznie wyprowadzone z α=2
vacuum Φ-EOM); (P) jako predykcja TGP dopóki S07 nie zamknie M9.1''
z fundamentu (zob. [[audyt/S07_M911_derivation/README.md]]).

**Konwencja PN counting.** "U³ deviation" w `g_tt` jest **3PN
w konwencji energy/PPN** (Will 1971; M9_1_pp_P1 cycle), co odpowiada
**2PN w konwencji phase** (Cutler–Flanagan 1994; standardowa GW
literature). Dla rejestru i ppE rekomendowana jest konwencja **phase**
(zob. [[audyt/T01_LIGO3G_falsifier/CONVENTION_DECISION.md]]).

**Mapowanie ppE (Yunes–Pretorius 2009).** Modyfikacja fazy
inspiralu w stationary phase approximation:

```
δΨ_TGP(f) = β_ppE^TGP · u^(b_ppE^TGP)        u = (πMf)^(1/3)
b_ppE^TGP = -1                                (2PN-phase, u⁻¹)

β_ppE^TGP^(b=-1) = -(3/(128·η)) · (5/6) · G_SPA   [LOCKED via op-ppE-mapping Phase 1]

Equal-mass (η=1/4), G_SPA=1 (central):
  β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 · 10⁻²

OOM window (G_SPA ∈ [0.7, 1.5], uncertainty z modified quadrupole formula):
  |β_ppE^TGP^(b=-1)| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹]
```

**Multi-coefficient TGP-distinguishing signature (M911-P2):**

```
β_(N+1)PN-phase / β_NPN-phase = {-23/10, -38/23, +337/228}
                              for N = 2, 3, 4 respectively

Wszystkie ratios FORCED przez α=2 + hyperbolic f(ψ) = (4-3ψ)/ψ;
0 free parameters w TGP. Inne modyfikacje GR (dCS, sGB, EÆ, BD)
dzielą b_ppE = −1 ale NIE reprodukują pattern bez fitting.
```

**Konkret: zob.** [[audyt/T01_LIGO3G_falsifier/PPN_TO_PPE_MAPPING.md]]
dla pełnego mapowania PN coefficients → ppE oraz
[[research/op-ppE-mapping/Phase1_results.md]] dla sympy LOCK
14/14 PASS (7/7 α_n^TGP + 7/7 Δα_n).

**Threshold detekcji (z β_ppE^TGP locked z op-ppE-mapping Phase 1).**
Bounds dla **2PN-phase, b_ppE=−1**; TGP prediction central
β_ppE^TGP ≈ -7.81·10⁻² (OOM 5.5–12·10⁻²):

| Detektor | Scenariusz | SNR | β_ppE_5σ_bound^(b=−1) | β_ppE^TGP / bound | Detekcja TGP |
|----------|-----------|-----|------------------------|---------------------|---------------|
| LIGO-O3 (now) | GW150914-like (M=65 M_⊙, 410 Mpc) | ~24 | ~10⁻¹ | ~0.78 | **borderline** |
| LIGO-O5 single (2027+) | Loud BBH (M=30 M_⊙, 200 Mpc) | ~80 | ~3·10⁻² | ~2.6 | **YES — first decisive 2027–2030** |
| LIGO-O5 stack 100 BBH | A+ design | ~800 | ~3·10⁻³ | ~26 | YES (>20σ) |
| ET-D single (~2035) | Loud BBH (M=30 M_⊙, 1 Gpc) | ~500 | ~10⁻³ | ~78 | YES (>50σ) |
| CE single (~2035+) | Loud BBH (M=30 M_⊙, 1 Gpc) | ~1000 | ~3·10⁻⁴ | ~260 | YES (>200σ) |
| ET+CE network stack 5000 BBH | (~5 yr) | network | ~10⁻⁴ | ~780 | YES (>700σ) |

**KEY UPDATE z B Phase 1 (2026-05-07):** TGP β_ppE^TGP ≈ 7.8·10⁻² jest
**>2× powyżej LIGO-O5 single-event bound** ~3·10⁻². To znaczy:
**LIGO-O5 *single event* w erze ~2027–2030 może już dać first
decisive detection M9.1''**, nie wymaga stack BBH events. Window of
testability accelerated od OOM "stack 2027–2030" do **"single event 2027+"**.

**Falsyfikacja (warunek explicit, z liczbowym threshold):**

> Jeśli ET-D + CE (lub LIGO-O5 single event z SNR ≥ 60), na próbie
> BBH events z **M_chirp ∈ [10, 50] M_⊙** w fazie inspiral
> (f ∈ [10, 100] Hz), **nie** detektuje ppE 2PN-phase coefficient
> β_ppE^(b=−1) zgodnego z M9.1'' predykcją:
>
>   |β_ppE^TGP^(b=-1)| ∈ [5.5·10⁻², 1.2·10⁻¹]   (OOM window, op-ppE-mapping Phase 1)
>   β_ppE^TGP^(b=-1) ≈ -5/64 ≈ -7.81·10⁻²       (central, η=1/4, G_SPA=1)
>
> z dokładnością ±5%, przy 5σ konfidencji, **M9.1'' jest sfalsyfikowana**.

**KRITYCZNY caveat (z GWTC-3 reanalysis 2026-05-07):**

> Detekcja M911-P1 wymaga **dedicated TGP-specific analysis pipeline**:
> single-coefficient Bayes prior z β_ppE^(b=-1) prior at -5/64 (lub
> alternativeness, multi-coefficient ratio test M911-P2 dla pattern
> {-23/10, -38/23, +337/228}). **Generic LIGO ToGR multi-coefficient
> marginalized analysis** (jak GWTC-2/GWTC-3 ToGR papers) jest factor
> ~50× weaker per coefficient i NIE detektuje TGP nawet w erze ET-D
> + CE. Wszystkie threshold values powyżej odnoszą się do
> *TGP-specific* analysis. Patrz [[../../research/op-GWTC3-reanalysis/Phase3_verdict.md]]
> §2 dla pełnego reconciliation z LIGO ToGR generic bounds (~factor 50× off).
> Falsyfikacja M9.1'' implikuje rewizję sek08c i sek08a (powrót do
> ścieżki S07 — derywacja metryki z fundamentu jako kolejna iteracja).

**Cykle odpowiedzialne:**
- `research/op-ppE-mapping/` — zamknięcie [β_th] (Phase 1 analytical;
  Phase 2 numerical Fisher cross-check)
- `research/op-LIGO-3G-deviation/` — Fisher analysis ET-D + CE,
  SNR thresholds, parameter degeneracy z spinem i M_chirp
- `audyt/T01_LIGO3G_falsifier/` — meta tracking, falsifier audit-bilet

**Cross-checks (orthogonal channels poza inspiralem 3PN):**
- BH5 (QNM ringdown δf/f ~ 8–16% w psi=1.20 ringdown) — zob.
  PREDICTIONS_REGISTRY: niezależny kanał strong-field deviation;
  spójny z M9.1'' jeśli M911-P1 trafiona
- ε.1 / E6 (5-channel ε photon-ring roadmap, ngEHT + LISA + LIGO O5)
- BH1 (EHT shadow size / brightness asymmetry)
- Wewnętrzny: cross-cycle R3 ≡ ψ=4/3 horizon (g₀_crit=1.874 w 4 cyfrach,
  zob. [[research/why_n3/PHASE1_psi_g0_identification.md]] Discovery 1)

**Dependency:**
- M911-P1 status (W) zależy od S07 closure (M9.1'' jako derywacja
  vs ansatz). Po S07 EXECUTED, M911-P1 awansuje do (W) bez zastrzeżeń.
  Przed S07 closure, M911-P1 jest **falsifier dla ansatzu**, nie dla
  fundamentu — ale to jest **wystarczające** dla testów ET/CE.

**Proponowana lokalizacja w `PREDICTIONS_REGISTRY.md`:**
- W sekcji "Strong-field gravitational tests" lub "Strong-field BH/GW
  predictions" (obok GW1, GW2, GW6, BH5).
- Sugerowany prefix: **M911-P1** (lub `T01-FALSIFIER`).
- Status w pierwszym wpisie: **LIVE** (target ET-D + CE), z notą
  "FALSIFIABLE-CONTRACT 2035+".

```

## §2 — Proponowana sekcja w `PREDICTIONS_REGISTRY.md` "Roadmap"

```markdown
| **2035+** | ET-D + CE inspirale BBH | **M911-P1** (3PN deviation \|Δg_tt\|=(5/6)U³ → β_ppE^(3PN); falsyfikacja jeśli inspiral 3PN coefficient niezgodny z β_ppE^TGP w 5%/5σ na N≥30 BBH events) |
```

(zaktualizować tabelę "Roadmap" na końcu rejestru — dodać wiersz po
istniejących E6 / UV6 / II6 wierszach).

## §3 — Komentarz do podejmowania decyzji autora

### Status epistemiczny tekstu §1, §2

- **Forma matematyczna deviation** ((5/6)U³ + wyższe rzędy):
  **DERIVED** — sympy LOCK 5/5, niezależny BVP residual test
  (M9_1_pp_P1_results.md §3.1).
- **Mapowanie na ppE** (β_ppE^TGP wartość liczbowa):
  **DRAFT** — dictionary w [[PPN_TO_PPE_MAPPING.md]] daje
  *strukturę* i *order of magnitude*, ale liczba [β_th] wymaga
  zamknięcia przez `op-ppE-mapping` Phase 1.
- **Threshold detekcji** (β_ppE Fisher bound dla ET-D / CE):
  **OOM** (order of magnitude) — back-of-envelope w
  [[SENSITIVITY_BACK_OF_ENVELOPE.md]] z literatury (Yagi–Yunes 2016,
  Chamberlain–Yunes 2017, Maggiore et al. ET science case 2020).
- **Falsifier clause** (warunek explicit):
  **NORMATIVE** — proponuje *jak* TGP weryfikuje M9.1'' przez ET/CE.
  Forma jest *kontraktem* na kształt sukcesu lub porażki teorii;
  liczbowy threshold zostanie skonkretyzowany z cykli A/B.

### Co autor musi zrobić, żeby wkleić §1, §2

1. **Zastąpić `[β_th]`** liczbą po zamknięciu `op-ppE-mapping/`
   Phase 1 (oczekiwane: rzędu 10⁻³–10⁻⁴ dla typowego loud BBH;
   exact value zależy od η, M_chirp normalization).
2. **Skopiować §1 jako nową sekcję** w `PREDICTIONS_REGISTRY.md`
   po wpisie BH5 lub w sekcji "Strong-field GW predictions".
3. **Skopiować §2** do tabeli "Roadmap 2030+ falsifiers".
4. **Bumpować counter** w PREDICTIONS_REGISTRY (obecnie ~856 → 857
   po dodaniu M911-P1).
5. (Opcjonalnie) **Dodać** w `TGP_FOUNDATIONS.md` § 3 link `M911-P1`
   po explicit "(5/6) U³" — żeby narrative było połączone z rejestrem.

### Argumenty za wpisaniem TERAZ vs po cyklu A/B

**Za (lekka adopcja teraz):**
- Path C w README.md jest **najlżejszą** akcją w T01.
- Predykcja matematyczna już jest DERIVED (M9_1_pp_P1) — wpis nie
  dodaje "fitted-anchor" do rejestru.
- T01 awansuje z `has_findings_file: false` na `has_findings_file:
  true` (audyt-folder reflectuje stan rzeczywisty).
- Falsifier-statement jest **publicznym kontraktem** TGP — sygnał
  do peer-review (Path D) i community.

**Za (wstrzymanie do zamknięcia A/B):**
- Bez liczbowego [β_th] wpis jest "pół-falsyfikowalny"; może być
  postrzegany jako preliminary.
- Pierwszy cytowalny wpis powinien być pełny.

**Rekomendacja:** wpis **częściowy** (z `[β_th]` jako placeholder
i etykietą `LIVE` + nota "[β_th] zostanie zamknięty z `op-ppE-mapping`
Phase 1; oczekiwana skala ~10⁻³–10⁻⁴") jest **akceptowalny** jako
honest reporting. Status `LIVE-PARTIAL` jest preferowany nad
`STRUCTURAL-only`, bo M9.1'' P1 jest *konkretną liczbą*, nie tylko
strukturą.

## §4 — Cytaty i pochodzenie

### Pochodzenie współczynników deviation

Wszystkie wartości w tabeli §1 (5/6, 23/12, 19/6, 337/72) są
**explicit policzone** w [[research/op-newton-momentum/M9_1_pp_P1_results.md]]
§3.2 z α=2 vacuum Φ-EOM, sympy LOCK 5/5 PASS, plus numerical
verification BVP residual (M9_1_pp_P1_results.md §3.1).

### Pochodzenie ppE framework

- Yunes & Pretorius, *Fundamental theoretical bias in gravitational-wave astronomy and the parameterized post-Einsteinian framework*, Phys. Rev. D **80**, 122003 (2009). arXiv:0909.3328.
- Yunes, Yagi, Pretorius, *Theoretical physics implications of the binary black-hole mergers GW150914 and GW151226*, Phys. Rev. D **94**, 084002 (2016). arXiv:1603.08955.
- Yagi, Yunes, *Approximate Universal Relations for Neutron Stars and Quark Stars*, Phys. Rep. **681**, 1 (2017).

### Pochodzenie sensitivity estimates

- LIGO Scientific Collaboration et al., *Tests of General Relativity with the Binary Black Hole Signals from the LIGO-Virgo Catalog GWTC-2*, Phys. Rev. D **103**, 122002 (2021). arXiv:2010.14529.
- Maggiore et al., *Science Case for the Einstein Telescope*, JCAP **03**, 050 (2020). arXiv:1912.02622.
- Reitze et al., *Cosmic Explorer: The U.S. Contribution to Gravitational-Wave Astronomy beyond LIGO*, Bull. Am. Astron. Soc. **51**, 035 (2019).
- Chamberlain & Yunes, *Theoretical Physics Implications of Gravitational Wave Observation with Future Detectors*, Phys. Rev. D **96**, 084039 (2017). arXiv:1704.08268.

## Cross-references

- [[README.md]] — diagnoza T01, ścieżki A/B/C/D
- [[NEEDS.md]] — N8 (registry edit) jako konkret tego draftu
- [[PPN_TO_PPE_MAPPING.md]] — Path B preview (detail dictionary)
- [[SENSITIVITY_BACK_OF_ENVELOPE.md]] — Path A preview (OOM detectability)
- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5 — recenzja źródłowa
- [[../../PREDICTIONS_REGISTRY.md]] — target dla wpisu (Path C)
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] — pochodzenie liczb 5/6, 23/12, …
- [[../../TGP_FOUNDATIONS.md]] § 3 — narrative source predykcji
- [[../S07_M911_derivation/README.md]] — meta-bloker B1
