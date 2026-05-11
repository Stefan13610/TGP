---
title: "Falsifier statement — M911-P1: post-Phase 1.5 + post-GWTC-3 RE-RUN native-first form"
date: 2026-05-10 (v2 native-first revision)
revision_history:
  - v1 (2026-05-07): Phase 1 OOM heuristic baseline (β=-5/64); SUPERSEDED
  - v1.5-update (2026-05-09): CRITICAL UPDATE block added; internally inconsistent
  - v2 (2026-05-10): native-first revision per `meta/PPN_AS_PROJECTION.md §3.1`;
                     OUTDATED Phase 1 baseline moved to APPENDIX A; current status
                     promoted to primary content
parent: "[[README.md]]"
type: registry-statement-native-first
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
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
  - native-observables-first
  - post-Phase1.5
  - post-GWTC3-RERUN
  - falsified-observational
related:
  - "[[README.md]]"
  - "[[ADDENDUM_2026-05-10_native_observables_first.md]] (native-first methodology binding)"
  - "[[NEEDS.md]]"
  - "[[PPN_TO_PPE_MAPPING.md]] (target of post-Phase 1.5 revision)"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]] (target of post-Phase 1.5 revision)"
  - "[[../../meta/PPN_AS_PROJECTION.md]] (parent methodology, binding 2026-05-10+)"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
  - "[[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (β=-15/4 sympy LOCK source)"
  - "[[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] (5.02σ FALSIFIED source)"
  - "[[../../research/op-emergent-metric-from-interaction-2026-05-09]] (recovery framework)"
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] (native c_n source)"
---

# Falsifier statement — M911-P1 (native-first v2, post-Phase 1.5 + post-GWTC-3 RE-RUN)

## §0 — Document status

**Current operational status (2026-05-10 v2):**

| Predykcja | Status | Source |
|---|---|---|
| **M911-P1** (g_tt = (5/6)U³ deviation, M9.1'' specific f(ψ)=(4-3ψ)/ψ form) | **FALSIFIED-OBSERVATIONAL 5.02σ** | GWTC-3 RE-RUN 2026-05-09 z β=-15/4 prior |
| **M911-P2** (multi-coefficient ratios {-23/10, -38/23, +337/228}) | **WITHDRAWN — needs re-derivation z corrected G_SPA** | Phase 1.5 alternative SPA shows β_3PN/β_2PN = -11161/504 ≠ -23/10 |
| **M911-P3** (4-channel orthogonal pattern) | **PARTIAL-FALSIFIED** (2/4 channels invalid; BH5/ε.1 remain) | per Phase 1.5 + GWTC-3 RE-RUN |
| **Recovery framework** (emergent-metric Phase 4) | **VIABLE** (zero-β region in (a_n, ξ_n) space) | op-emergent-metric Phase 4 STRUCTURAL DERIVED |

**Document organization:**
- §1: native-first falsifier statement (PRIMARY content, gotowy do PREDICTIONS_REGISTRY)
- §2: PREDICTIONS_REGISTRY recommended entries (corrected post-2026-05-09)
- §3: Roadmap update (post-falsification + recovery)
- §4: Cross-references
- **APPENDIX A**: OUTDATED Phase 1 baseline (β=-5/64) — historical reference only

## §1 — Native-first falsifier statement (CURRENT)

### §1.1 — Native predictions (L1)

TGP M9.1''-derived **native Taylor coefficients** of `g_tt[Φ]` w isotropic
weak-field expansion (sympy LOCK 5/5 PASS, op-newton-momentum/M9_1_pp_P1_results.md
§3.2):

```
g_tt^TGP / (-c²) = 1 − 2U·c_1 + 2U²·c_2 + 2U³·c_3 + 2U⁴·c_4 + 2U⁵·c_5 + ...

c_1 = +1                  c_2 = +1                  (forced z PPN γ=β=1)
c_3 = -5/6                c_4 = -23/12              c_5 = +337/72
                          (DEVIATIONS — derive z α=2 vacuum Φ-EOM, M9.1'' f(ψ)=(4-3ψ)/ψ)
```

Plus structural prefactor lock (op-ppE-mapping/Phase1.5_G_SPA_lock.md, sympy LOCK 5/5):

```
G_SPA = 48                (SPA chain prefactor in M9.1'' isotropic test-particle limit)
                          [NOT G_SPA ≈ 1 generic-ppE heuristic]
```

### §1.2 — ppE projection (L2)

Per Phase 1.5 sympy LOCK (and cross-validated z emergent-metric Phase 4
independent derivation):

```
β_ppE^TGP_(b=-1) = -(3/(128·η)) · Δα_3 · G_SPA
                 = -(3/(128·1/4)) · (-5/6) · 48
                 = -15/4 ≈ -3.75    [η=1/4 binary inspiral, M9.1'' specific point]

Equivalent ToGR fractional:
δφ̂_4_TGP = -0.865    (86% deviation w 2PN-phase coefficient)
```

Convention: PHASE-PN (Cutler-Flanagan 1994); M9.1'' convention dictionary w
[[CONVENTION_DECISION.md]]; b_ppE = -1 odpowiada 2PN-phase correction (= U³ in
energy convention, = "3PN energy" in M9.1'' historical naming).

### §1.3 — Falsification map (L3)

| Bound | Native coef constrained | Status |
|---|---|---|
| Cassini γ_PPN ≤ 2.3·10⁻⁵ | c_1 = 1 (M9.1'' classical) | PASS automatic |
| LLR β_PPN ≤ 8·10⁻⁵ | c_2 = 1 plus (a_2, ξ_2) | PASS automatic |
| **GWTC-3 RE-RUN 2026-05-09** (corrected β=-15/4 prior) | M9.1'' specific point {c_3, c_4, c_5} | **FALSIFIED 5.02σ** (BF=3.5·10⁻⁶, log10(BF)=-5.45) |
| Recovery Path 1: c_0=0 + ξ_3=(32-a_3)/32 | shifted (a_n, ξ_n) point | VIABLE (zero-β solution) |
| Recovery Path 2: M9.1'' a_n + c_0·κ_σ=4/3 | M9.1'' coefs + canonical σ-coupling | VIABLE (structurally preferred per emergent-metric §5) |

### §1.4 — Native-first falsifier statement (verbatim contract)

> **Predykcja**: TGP M9.1'' canonical f(ψ)=(4-3ψ)/ψ przewiduje native Taylor
> coefs `g_tt[Φ]`: c_3=-5/6, c_4=-23/12, c_5=+337/72 (sympy LOCK 5/5,
> op-newton-momentum/M9_1_pp_P1_results.md §3.2). Strukturalny prefactor
> SPA chain G_SPA=48 (sympy LOCK 5/5, op-ppE-mapping/Phase1.5).
>
> **Projekcja na ppE 2PN-phase chart:** β_ppE^TGP_(b=-1) = -15/4 ≈ -3.75
> (η=1/4 binary inspiral; cross-validated z emergent-metric Phase 4
> independent derivation: identical value, different route).
>
> **Falsyfikacja (current observational status):** GWTC-3 combined ~90 BBH
> posterior z corrected β=-15/4 prior (op-GWTC3-reanalysis Phase 2 RE-RUN
> 2026-05-09): BF_TGP/GR = 3.5·10⁻⁶, log10(BF) = -5.45, **TGP M9.1''
> RULED OUT at 5.02σ ("OVERWHELMING GR preference")**. M9.1'' specific
> form FALSIFIED-OBSERVATIONAL.
>
> **Recovery framework:** falsyfikacja punktu w przestrzeni native Taylor
> coefs (a_n, ξ_n, b_n, c_0, κ_σ) NIE wyklucza całej przestrzeni.
> op-emergent-metric-from-interaction-2026-05-09 (closed STRUCTURAL DERIVED,
> sympy 57/57 PASS) daje rodzinę Taylor expansions parametryzowanych przez
> te native coefs z **zero-β region** zgodnym z GWTC-3 1σ window
> |β_ppE_(b=-1)| ≤ 0.78. Path 1 (c_0=0 + 3PN tuning ξ_3=(32-a_3)/32) i
> Path 2 (M9.1'' a_n + canonical σ-coupling c_0·κ_σ=4/3) są strukturalnie
> viable.
>
> **Future falsifier (ET-D / CE 2027-2035):** confirm/refute *recovered
> specific point* (post-S07-derived alternative f(ψ) lub explicit σ-coupling
> via canonical c_0·κ_σ=4/3); NIE testuje już original M9.1'' canonical
> (które jest już wykluczone).
>
> **Falsifier kontrakt z observatoriumami:** jeśli przyszłe data (ET-D + CE
> 2027-2035) pokaże, że **żaden** punkt w native (a_n, ξ_n) space jest
> consistent z observed β_ppE_(b=-1) within 5σ, **emergent-metric framework
> byłby strukturalnie wykluczony** (nie tylko specific M9.1'' point).
> To byłaby strukturalna falsyfikacja TGP gravity sector.

## §2 — Proponowana sekcja w `PREDICTIONS_REGISTRY.md`

```markdown
### M911-P1 — TGP M9.1'' canonical metryka (FALSIFIED-OBSERVATIONAL)

**Native predykcja (preserved):** Taylor coefs `g_tt[Φ]` z M9.1'' f(ψ)=(4-3ψ)/ψ:
- c_3 = -5/6, c_4 = -23/12, c_5 = +337/72 (sympy LOCK 5/5)
- G_SPA = 48 (sympy LOCK 5/5, SPA chain in M9.1'' isotropic test-particle)

**ppE projection:** β_ppE^TGP_(b=-1) = -15/4 (η=1/4 binary inspiral)
**Cross-validation:** emergent-metric Phase 4 niezależna derywacja → identical -15/4

**Status:** **FALSIFIED-OBSERVATIONAL 5.02σ** (GWTC-3 RE-RUN 2026-05-09,
BF_TGP/GR = 3.5·10⁻⁶, "OVERWHELMING GR preference")

**Recovery framework:** op-emergent-metric-from-interaction-2026-05-09 daje
rodzinę Taylor expansions z zero-β region zgodnym z GWTC-3 1σ. Two viable paths:
- Path 1: c_0=0 + 3PN tuning ξ_3=(32-a_3)/32
- Path 2: M9.1'' a_n + canonical σ-coupling c_0·κ_σ=4/3 (structurally preferred)

**Pending closure:** S07 audyt (M9.1'' derivation z fundamentu) + emergent-metric
specific point identification (post-recovery).
```

## §3 — Roadmap update

```markdown
| **2026-05-09** | GWTC-3 RE-RUN | M9.1'' canonical FALSIFIED-OBSERVATIONAL 5.02σ; recovery via emergent-metric framework opened |
| **2027-2030** | LIGO-O5 stack | confirm/refute recovered point in (a_n, ξ_n) space (Path 1 lub Path 2) |
| **2035+** | ET-D + CE inspirale BBH | decisive 5σ test of recovered specific point post-S07 |
```

## §4 — Cross-references

### Native-first methodology binding
- [[../../meta/PPN_AS_PROJECTION.md]] — parent methodology (binding 2026-05-10+)
- [[ADDENDUM_2026-05-10_native_observables_first.md]] — T01 native-first overlay
  (chart correction Phase 1 → 1.5)

### Source documents (Phase 1.5 + GWTC-3 RE-RUN)
- [[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] — β=-15/4 sympy LOCK
- [[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]
  — 5.02σ FALSIFIED source
- [[../../research/op-emergent-metric-from-interaction-2026-05-09]] — recovery
  framework, sympy 57/57 PASS

### Native predictions sources
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] — c_n derivation,
  sympy 5/5 LOCK
- [[../../TGP_FOUNDATIONS.md]] § 3 — (5/6)U³ narrative

### Audit-level
- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5 — recenzja źródłowa T01
- [[../README.md]] — indeks audytu
- [[../S07_M911_derivation/README.md]] — meta-bloker (M9.1'' z fundamentu)

### Pre-revision (historical)
- See **APPENDIX A** below for OUTDATED Phase 1 baseline (β=-5/64).

---

## APPENDIX A — OUTDATED Phase 1 baseline (β=-5/64) — HISTORICAL ONLY

> **⚠ HISTORICAL REFERENCE ONLY:**
> Sections below describe the **OUTDATED Phase 1 OOM heuristic** baseline
> (β_ppE^TGP_(b=-1) ≈ -5/64 ≈ -7.81·10⁻², z G_SPA ≈ 1 założenia z
> Sampson-Yunes-Cornish 2013).
>
> **Diagnoza błędu:** Phase 1 cited SYC 2013 dla "G_SPA ≈ 1 for metric-only
> modifications". SYC 2013 framework był derived dla *small-perturbation*
> regimes (BD with 1/ω_BD coupling, dCS with ζ_dCS coupling), **NIE dla
> structural O(1) modifications**. TGP M9.1'' f(ψ)=(4-3ψ)/ψ z Δα_3=-5/6
> O(1) jest w drugiej kategorii. SPA chain prefactor amplifies (α_4=30·e_2
> + cross-terms) → G_SPA=48, factor 48× away from 1.
>
> Phase 1.5 sympy-LOCK derived G_SPA=48 explicit; β_ppE^TGP_(b=-1) corrected
> to -15/4. GWTC-3 RE-RUN z corrected prior dał 5.02σ falsification.
>
> Original §1-§4 content (preservowane dla audit trail; **NIE do PREDICTIONS_REGISTRY**)

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
