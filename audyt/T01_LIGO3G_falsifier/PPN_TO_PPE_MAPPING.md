---
title: "PPN → ppE mapping (Path B): (5/6) U³ deviation → β_ppE^(2PN-phase) — Phase 1.5 G_SPA=48 LOCK"
date: 2026-05-07 (original) / 2026-05-10 (v2 G_SPA=48 LOCK update)
parent: "[[README.md]]"
type: analytical-preview-superseded-by-Phase1.5
tgp_owner: audyt/T01_LIGO3G_falsifier
revision_history:
  - v1 (2026-05-07): preview Path B — Phase 1 OOM heuristic, β_ppE^TGP ≈ -5/64 (G_SPA ≈ 1)
  - v2 (2026-05-07 sesja C-B-A-D): preview SUPERSEDED przez Phase 1 LOCK (op-ppE-mapping/Phase1_results.md)
  - v3 (2026-05-09): Phase 1.5 G_SPA=48 sympy-LOCK derived → β_ppE^TGP = -15/4 ≈ -3.75 (factor 48× larger)
  - v4 (2026-05-10): native-first reframe + dictionary explicit G_SPA=48 update post-Phase 1.5 +
                     post-GWTC-3 RE-RUN; aligns with [[ADDENDUM_2026-05-10_native_observables_first.md]]
tags:
  - dictionary
  - analytical
  - ppE
  - PPN
  - 2PN-phase
  - M911
  - inspiral
  - phase-modification
  - T01
  - EXT-5
  - G_SPA-locked
  - sympy-LOCK
  - Phase1.5
  - native-observables-first
related:
  - "[[README.md]]"
  - "[[ADDENDUM_2026-05-10_native_observables_first.md]] (parent native-first methodology)"
  - "[[NEEDS.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
  - "[[../../research/op-ppE-mapping/Phase1_results.md]] (Phase 1 LOCK, β=-5/64 with G_SPA≈1 heuristic)"
  - "[[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (Phase 1.5 sympy-LOCK 5/5 PASS, G_SPA=48 sympy-exact)"
  - "[[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]] (GWTC-3 RE-RUN, 5.02σ FALSIFIED)"
  - "[[../../meta/PPN_AS_PROJECTION.md]] (binding methodology 2026-05-10+)"
---

# PPN → ppE mapping — Path B (post-Phase 1.5 LOCK)

## ⚠ STATUS UPDATE 2026-05-10 — Phase 1.5 G_SPA=48 LOCK incorporated

> **Co zmieniło się od v1 preview (2026-05-07):**
>
> 1. **G_SPA = 48 sympy-exact** (NIE ≈ 1 jak Phase 1 heuristic zakładał).
>    Phase 1.5 derivation 4-level verified (sympy LOCK 5/5 + independent
>    hand-calc + numerical sanity + alternative SPA orthogonal route).
>    Source: [[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] §4.5.
>
> 2. **β_ppE^TGP^(b=-1) = -15/4 ≈ -3.75** (test-particle exact, η=1/4 ±25%);
>    **factor 48× LARGER** niż Phase 1 OOM estimate -5/64 ≈ -0.078.
>
> 3. **Multi-coefficient ratios MISLEADING.** Phase 1 ratios {-23/10, -38/23,
>    +337/228} są INCORRECT; Phase 1.5 alternative SPA gives β_3PN/β_2PN =
>    -11161/504 ≈ -22.14 (factor ~10× different). M911-P2 multi-coefficient
>    pattern requires full re-derivation per Phase 1.5 §10 sign-off.
>
> 4. **GWTC-3 RE-RUN 2026-05-09 with corrected β**: TGP M9.1'' specific
>    (4-3ψ)/ψ ansatz **FALSIFIED at 5.02σ** (BF=3.5·10⁻⁶, log10 BF = -5.45).
>    Source: [[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]].
>
> 5. **Falsifier window collapsed:** "ET-D + CE 2035+" no longer applies;
>    falsification has occurred in CURRENT LIGO-O3 / GWTC-3 data. The
>    "preview Path B for op-ppE-mapping cycle" framing is HISTORICAL.
>
> **Dictionary entries below are UPDATED to Phase 1.5 LOCK values.**
> Sections preserved as analytical reference; the cycle `op-ppE-mapping/`
> Phase 1 + 1.5 closed all open assumptions. Native-first methodology
> per [[ADDENDUM_2026-05-10_native_observables_first.md]] (binding
> 2026-05-10+) treats `g_tt[Φ]` Taylor coefficients (c_3 = -5/6,
> c_4 = -23/12, c_5 = +337/72, …) as **L1 native** observables; ppE
> parameters (β_ppE, b_ppE, G_SPA) as **L2 projection chart**;
> GWTC-3 / ET-D / CE bounds as **L3 falsification map**.

## §0 — Native-first three-layer dictionary (PRIMARY post-2026-05-10)

| Layer | Quantity | TGP value | Source |
|---|---|---|---|
| **L1 native** | c_3 (g_tt Taylor at U³) | **-5/6** sympy-exact | M9_1_pp_P1 §3.2 LOCK 5/5 |
| **L1 native** | c_4 (g_tt Taylor at U⁴) | **-23/12** sympy-exact | M9_1_pp_P1 §3.2 LOCK 5/5 |
| **L1 native** | c_5 (g_tt Taylor at U⁵) | **+337/72** sympy-exact | M9_1_pp_P1 §3.2 LOCK 5/5 |
| **L1 native** | Δe_2 (orbital binding at v⁴) | **-4/3** sympy-exact | Phase 1.5 §2.4 LOCK L3 |
| **L1 native** | Δα_4 (TaylorF2 phase coef) | **-40** sympy-exact | Phase 1.5 §4.3 LOCK L5 |
| **L2 chart** | b_ppE (PN order encoding) | **-1** (2PN-phase) | Cutler-Flanagan / Yunes-Pretorius |
| **L2 chart** | G_SPA (SPA chain prefactor) | **48** sympy-exact (test-p) | Phase 1.5 §4.5 LOCK L5 |
| **L2 chart** | β_ppE^TGP^(b=-1) at η=1/4 | **-15/4 ≈ -3.75** ±25% | Phase 1.5 §4.4 LOCK L5 |
| **L3 falsifier** | GWTC-3 BF (TGP vs GR) | **3.5·10⁻⁶ → 5.02σ FALSIFIED** | RE-RUN 2026-05-09 |

**Native parameter audit:** native coefs {c_3, c_4, c_5, …} are forced
by α=2 vacuum Φ-EOM + M9.1'' f(ψ)=(4-3ψ)/ψ ansatz; **ZERO free
parameters**. ppE values are **derived projections**, not fitted.
This is what makes the falsification test sharp.

---

## §0.1 — Cel pliku (HISTORICAL preview-era language)

> **Original cel pliku (2026-05-07).** Dostarczyć **analityczny szkielet**
> mapowania (5/6) U³ deviation w `g_tt^TGP` na parametry ppE
> Yunes–Pretorius (β_ppE, b_ppE) dla fazy waveformu inspiralu BBH.
> Ten plik był **previewem** dla cyklu `research/op-ppE-mapping/`
> Phase 1; cykl został wykonany w sesji 2026-05-07 (Phase 0+1+2+3
> sympy LOCK 14/14) + Phase 1.5 (2026-05-09 sympy LOCK 5/5 +
> 4-level verification G_SPA=48).
>
> **Status:** ANALYTICAL DICTIONARY (post-Phase 1.5 LOCK). Liczby są
> teraz sympy-exact. Sekcje §1-§6 poniżej zachowują strukturę preview
> jako pedagogical reference; §0 (native-first dictionary) jest
> primary post-2026-05-10.

## §1 — ppE framework w skrócie

Yunes & Pretorius (2009, Phys. Rev. D 80:122003) parametryzują
modyfikacje GR w stationary phase approximation (SPA) inspiralu
BBH przez:

```
h̃_ppE(f) = h̃_GR(f) · (1 + α_ppE · u^a)  ·  exp[ i β_ppE · u^b ]
```

gdzie:
- **u = (πMf)^(1/3)** — bezwymiarowy parametr PN (M = chirp-like mass
  scale, dla równych mas ≈ M_total; precyzyjnie: M_chirp = (m₁m₂)^(3/5)/M^(1/5)·η^(3/5)·…),
- **α_ppE, a** — modyfikacja amplitudy,
- **β_ppE, b** — modyfikacja fazy,
- a, b liczby całkowite kodujące **PN order** modyfikacji.

### Mapping b_ppE ↔ PN order

GR phase w SPA (Cutler & Flanagan 1994):

```
Ψ_GR(f) = 2πf t_c − Φ_c − π/4
        + (3/(128 η)) · u⁻⁵ · [1 + α_2 u² + α_3 u³ + α_4 u⁴ + α_6 u⁶ + α_7 u⁷ + …]
                            ↑       ↑       ↑       ↑       ↑
                          0PN     1PN     1.5PN    2PN     3PN     3.5PN
```

(używam tu konwencji "absolute u-power"; współczynnik α_n PN = u^n
relative to leading u^(-5)).

Dla deviation w **N-th PN order** w fazie:
- **b_ppE = N − 5** (bo phase ~ u^(-5+N) i ppE wyjmuje to z normalization
  jako β · u^b).

**Lookup table:**

| PN order (in phase) | b_ppE | u^b power |
|---------------------|-------|-----------|
| 0PN (Newtonian)     | −5    | u⁻⁵       |
| 1PN                 | −3    | u⁻³       |
| 1.5PN               | −2    | u⁻²       |
| 2PN                 | −1    | u⁻¹       |
| 2.5PN               |  0    | u⁰        |
| **3PN**             | **+1** | **u¹**   |
| 3.5PN               | +2    | u²        |
| 4PN (target ET/CE)  | +3    | u³        |

**Konwencja audytu T01 README/M9_1_pp_P1.** Uwaga: README opisuje
"U³ to człon 3PN" w **konwencji metrycznej** (1PN = O(U), 2PN = O(U²),
3PN = O(U³) — power potencjału U). To **nie jest** standardowa
konwencja Cutler–Flanagan w fazie waveformu. Dla pełnej spójności:

| Power w U (metric coefficient) | PN w fazie waveformu | b_ppE |
|-------------------------------|----------------------|-------|
| U (0PN metric)                | 0PN phase            | −5    |
| U² (1PN metric)               | 1PN phase            | −3    |
| **U³ (2PN metric)**           | **2PN phase**        | **−1** |
| U⁴ (3PN metric)               | 3PN phase            | +1    |
| U⁵ (4PN metric)               | 4PN phase            | +3    |

**KRYTYCZNA UWAGA:** istnieją *dwie konwencje* PN counting (energy
vs phase). README T01 i M9_1_pp_P1 piszą "U³ → 3PN" w konwencji
energetycznej (pełny PN ≡ U^N relative to leading); standardowa
literatura inspiral piszą tę samą deviation jako "2PN phase
correction" w konwencji 0PN-phase = u⁻⁵. **Cykl `op-ppE-mapping/`
musi zatwierdzić jedną konwencję** dla rejestru i waveformów.

W tym dokumencie **przyjmuję konwencję PHASE** (b_ppE map standardowy
jak Yunes 2009): wtedy **U³ deviation w `g_tt` → 2PN phase correction
→ b_ppE = −1**. Cykl zwalidować.

## §2 — Skąd phase deviation z metric deviation?

Phase modification w SPA przepływa od metric → orbital energy E(v) →
luminosity dE/dt → df/dt → t(f) → Ψ(f). Dla "metric-only" deviation
(modyfikacja samego g_tt bez nowych pól radiacyjnych) struktura jest:

### 2.1 Standard SPA chain (Will–Wiseman 1996; Mishra et al. 2016)

```
E_orb(v) = -½ μ v² · [ 1 + e_2(η)·v² + e_3(η)·v³ + e_4(η)·v⁴ + … ]
                                ↑               ↑
                              1PN             2PN

dE/dt = -32/5 · η²·v¹⁰/G · [ 1 + f_2(η)·v² + f_3(η)·v³ + f_4(η)·v⁴ + … ]
                                          ↑                 ↑
                                       1PN              2PN

Ψ(f) = ∫ 2πf · (df/dt)⁻¹ df  (analytical w ν = (πMf)^(1/3) → u)
```

W M9.1'' modyfikacja jest tylko w `g_tt` (i konsekwentnie w
`g_ij` przez f·h=1). Energia orbitalna na okrągłej orbicie wokół
pojedynczej masy w hiperbolicznej metryce:

```
E_orb^TGP = -½ μ v² · [1 + (2/3) η v² + (... TGP-specific) v⁴ + …]
                                          ↑
                                  TGP odbiega od GR od v⁴ (2PN orbital)
```

Dla **dwóch** ciał potrzebny jest pełny rachunek M9.1'' two-body
problem — TO JEST ZAKRES `op-LIGO-3G-deviation/` Phase 1 (N1).
Ale dla preview szacujemy: deviation w E_orb przy v⁴ propaguje się
do fazy Ψ przy u⁴/u⁵ = u⁻¹, czyli **2PN phase deviation** (b_ppE = −1).

### 2.2 Heurystyka magnitude (HISTORICAL — superseded by Phase 1.5 LOCK)

> **Phase 1 heuristic (2026-05-07, OBSOLETE):** preview szacował
> κ_TGP ≈ 0.5–1.5 (G_SPA ≈ 1) na podstawie Sampson-Yunes-Cornish 2013
> "metric-only modifications" framework. **To było incorrectly applied
> outside SYC 2013 regime of validity** (small-perturbation; TGP
> M9.1'' jest structural O(1) modification). Phase 1.5 sympy-LOCK
> 4-level verified (2026-05-09): **G_SPA = 48 sympy-exact**
> (factor 48× larger than Phase 1 heuristic).
>
> Phase 1.5 derivation (sketch):
> ```
> Δe_2 (orbital binding at v⁴) = (8/5)·Δα_3_metric = (8/5)·(-5/6) = -4/3
> Δα_4 (TaylorF2 phase)         = 30·Δe_2 + cross-terms              = -40
> G_SPA = Δα_4 / Δα_3            = -40 / (-5/6)                       = 48
> β_ppE^TGP^(b=-1) at η=1/4      = -(3/(128·1/4)) · 40                = -15/4
> ```
> Source: [[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] §2-§4.

Dla M9.1'':

```
Δg_tt/c² = − (5/6) U³ + …      gdzie U ≈ v²/c² ≈ u²/c²
```

W phase modification (Phase 1.5 sympy-LOCKED SPA chain):

```
δΨ_TGP(u) = (3/(128 η)) · u⁻⁵ · κ_TGP^LOCKED · u⁴      [PHASE 1.5]

κ_TGP^LOCKED = (5/6) · G_SPA = (5/6) · 48 = 40           (test-p exact)
β_ppE^TGP^(b=-1) at η=1/4 = -(3/(128·η)) · κ_TGP = -15/4 ≈ -3.75
                                                    [test-p ± 25% η-correction]
```

**Phase 1.5 OOM window (test-particle ± η-correction):**
```
|β_ppE^TGP^(b=-1)| ∈ [2.81, 4.69]                         (NIE [0.055, 0.12])
```

**Phase 1.5 closed all preview-era assumptions** A1, A2, A3, A4, A5, A6
(see §5 below; status updated post-Phase 1.5).

## §3 — Dictionary: PN coefficients M9.1'' → ppE phase

Z M9_1_pp_P1_results.md §3.2 wszystkie deviation coefficients:

| Power U^N | Δ(α_n^TGP − α_n^GR) | Konwencja "metric PN" | Konwencja "phase PN" | b_ppE |
|-----------|---------------------|-----------------------|----------------------|-------|
| U³        | **− 5/6**           | 3PN (energetyczna)    | 2PN (phase)          | **−1** lub **+1** (decision) |
| U⁴        | + 23/12             | 4PN                   | 3PN                  | +1 lub +3 |
| U⁵        | − 19/6              | 5PN                   | 4PN                  | +3 lub +5 |
| U⁶        | + 337/72            | 6PN                   | 5PN                  | +5 lub +7 |

**Każdy z tych coefficients może być ppE entry** — multi-deviation
ansatz Yunes–Pretorius (2009) zezwala na sumę:

```
δΨ(u) = Σ_k β_k · u^(b_k)
```

To znaczy: TGP M9.1'' ma **całą "drabinę" ppE coefficients** o
specyficznym wzorze (5/6, 23/12, 19/6, 337/72, ...) — to jest
**charakterystyczny TGP signature** odróżniający od generic 3PN
modifications (które są pojedyncze).

**Implication dla falsifikacji:** ET-D + CE może testować *spójność*
M9.1'' nie tylko przez jeden b_ppE ale przez **wzajemny stosunek**
(β_3PN, β_4PN, β_5PN). Jeśli detektory widzą deviation na 3PN
zgodne z β_TGP^(3PN), ale *NIE* na 4PN zgodne z β_TGP^(4PN), to
falsyfikuje M9.1''. To jest **silniejszy** test niż single-coefficient.

## §4 — Konkurencja modeli (T01 problem #3)

### 4.1 Catalog of ppE-mapped 3PN modifications (Yunes–Yagi–Pretorius 2016)

| Theoretical bias | b_ppE | Physical origin |
|------------------|-------|-----------------|
| Brans–Dicke 1PN dipole | −7 | Scalar dipole radiation |
| Massive graviton | −3 | Yukawa-type modification |
| Lorentz-violation (vSME) | −5, −3, −1 | Cumulative ν-dependent |
| Quadratic gravity (sGB, dCS) | −1 | Scalar-tensor higher curvature |
| Einstein-Æther | −1 | Vector preferred-frame |
| Nonlocal gravity (1/□) | various | Memory effects |
| **TGP M9.1''** | **−1** (or +1) | **Hyperbolic g_tt potential structure** |

**Konkluzja:** TGP-specific signature **NIE** jest "świeży" b_ppE
sam w sobie (b = −1 dzielone z dCS, sGB, Einstein–Æther). Co czyni
TGP **distinguishable**:

1. **Wzorzec multi-coefficient** (5/6, 23/12, 19/6, …) — reprodukowalny
   *tylko* z α=2 hyperbolic ansatz, nie z innych modyfikacji.
2. **Cross-channel BH5** (QNM ringdown δf/f ~ 8–16% przy ψ=1.20)
   — TGP-specific, NIE w innych modelach.
3. **Brak parameter freedom** — w M9.1'' (5/6) jest *wymuszone*
   przez α=2 i hyperbolic f(ψ); nie ma "fitting freedom". W dCS i
   Einstein–Æther coefficient zależy od free coupling.

**Dictionary entry (LOCKED post-Phase 1.5, 2026-05-09):**

```
M9.1'' → ppE (Phase 1.5 sympy-LOCK 5/5 + 4-level verification):
  b_ppE^TGP            = -1                         (2PN-phase, U³ metric)
  G_SPA                = 48                         (sympy-exact, test-particle)
  β_ppE^TGP^(b=-1)     = -15/4 ≈ -3.75              (η=1/4, test-p ±25%)
  |β| OOM window       ∈ [2.81, 4.69]               (test-p ± η-correction)

Multi-coefficient (Phase 1.5 alternative SPA, M911-P2 needs full re-derivation):
  β_3PN/β_2PN          = -11161/504 ≈ -22.14        (NIE -23/10 jak Phase 1 heuristic
                                                      sugerował; Phase 1.5 §10 sign-off)
  β_4PN/β_3PN, β_5PN/β_4PN = TBD (M911-P2 WITHDRAWN-needs-rederivation)

Distinguishing signature: NIE jest zachowywany jako simple {-23/10, -38/23, +337/228}
  pattern (Phase 1 heuristic incorrect). Native multi-coefficient pattern requires
  full DJS 2-body Lagrangian (Phase 1.6+ work).
```

**Dictionary entry (HISTORICAL Phase 1 baseline, 2026-05-07 — preserved for revision history):**

```
M9.1'' → ppE [Phase 1 OOM heuristic, OBSOLETE]:
  b_ppE^TGP^(b=-1)  = -5/64 ≈ -0.078            (G_SPA ≈ 1 SYC 2013 heuristic, η=1/4)
  multi-pattern     = {-23/10, -38/23, +337/228}  (Phase 1 heuristic, INCORRECT)

Reason for retraction: SYC 2013 "G_SPA ≈ 1" applies only to small-perturbation
regimes (BD with 1/ω_BD coupling, dCS with ζ_dCS coupling). TGP M9.1'' f(ψ)=
(4-3ψ)/ψ is structural O(1) modification → SPA chain amplifies via α_4 = 30·e_2
+ cross-terms → G_SPA = 48 sympy-exact (Phase 1.5 §4.5).
```

## §5 — Założenia, które cykl `op-ppE-mapping/` musiał zwalidować — STATUS POST-PHASE 1.5

| ID | Założenie | Status (post-Phase 1.5 2026-05-09) |
|----|-----------|------------------------------------|
| A1 | M9.1'' two-body Lagrangian istnieje i daje finite Lagrangian (nie blow-up przy ψ → 4/3 horizon w external potencjale) | **PARTIAL** — test-particle limit OK (Phase 1.5 §2); equal-mass DJS 2-body deferred do Phase 1.6 (~25% uncertainty na η-correction) |
| A2 | Quadrupole formula w M9.1'' modyfikuje się w sposób reprezentatywny przez metric perturbation tylko (NIE nowymi dipole / scalar charges) | **CLOSED PASS** — Phase 1.5 §3.1 LOCK L4: F_TGP(v) = F_GR(v) at leading 2PN-orbital (cross-channel via GW1 c_T=c_s + GW2 3 DOF) |
| A3 | dE/dt luminosity w M9.1'' jest *modyfikowana* przy 2PN+ orbital w sposób spójny z metric coefficients (5/6) U³ | **CLOSED PASS** — Phase 1.5 §3.2: ΔF(v) = O(v⁴) ~ 3% at v_LSO, sub-leading dla 2PN-phase coef |
| A4 | Stationary phase approximation jest ważna w M9.1'' z tego samego powodu co w GR (adiabatic inspiral) | **CLOSED PASS** — Phase 1.5 §4 SPA chain inversion sympy-verified (cross-check Buonanno-Iyer 2009 GR test-p α_4) |
| A5 | TGP nie wprowadza nowych radiacyjnych DOF (no scalar mode) zgodnie z GW1 (c_T = c_s) i GW2 (3 DOF) | **DERIVED** w PREDICTIONS_REGISTRY (M911-P2/P3 partial; GW1/GW2 separately verified) |
| A6 | Konwencja PN counting "metric N-PN ↔ phase (N−1)-PN" jest poprawna dla M9.1'' (z hyperbolic ansatz) | **CLOSED CHOICE** — [[CONVENTION_DECISION.md]] adoptuje PHASE-PN; U³ metric = 2PN-phase = b_ppE = -1 |

**Phase 1.5 zamknął A1 (partial test-p), A2-A6 (full).** Pozostała ~25%
η-correction uncertainty na G_SPA(η=1/4) jest PRZEJRZYSTA — dla
falsification verdict (β_TGP/σ_β > 5σ) jest *nieistotna* (factor 48
correction już dominuje).

**Critical post-Phase 1.5 finding:** A1 partial-closure (test-p G_SPA=48
vs equal-mass G_SPA(η=1/4) ≤ ±25%) sufficient dla GWTC-3 RE-RUN, który
gives 5.02σ FALSIFIED-OBSERVATIONAL bez względu na η-correction precision.

## §6 — Status cyklu `op-ppE-mapping/` — POST-PHASE 1.5 LOCKED

**Phase 1 (2026-05-07):** β_ppE^TGP = -5/64 ≈ -0.078 (G_SPA ≈ 1 SYC 2013
heuristic) — **OBSOLETE** (regime-of-validity error).

**Phase 1.5 (2026-05-09):** β_ppE^TGP = -15/4 ≈ -3.75 (G_SPA = 48 sympy-exact
test-particle) — **LOCKED** (4-level verification: sympy LOCK 5/5 + hand-calc
+ numerical sanity + alternative SPA orthogonal route).

### Następne kroki post-Phase 1.5 (uporządkowane priority)

| Step | Description | Status |
|------|-------------|--------|
| **GWTC-3 RE-RUN** (highest priority) | Re-run TIGER-framework Bayes inference z β=-15/4 prior | **EXECUTED 2026-05-09** — BF=3.5·10⁻⁶, 5.02σ FALSIFIED-OBSERVATIONAL ([[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]) |
| **Falsifier statement revision** | Update FALSIFIER_STATEMENT_DRAFT z post-Phase 1.5 + post-RE-RUN values | **EXECUTED 2026-05-10 v2** — native-first format ([[FALSIFIER_STATEMENT_DRAFT.md]]) |
| **PREDICTIONS_REGISTRY update** | M911-P1 → FALSIFIED-OBSERVATIONAL; M911-P2 → WITHDRAWN; M911-P3 → PARTIAL | **EXECUTED 2026-05-09** |
| **Paper draft revision** | papers/M911_LIGO3G_paper "predictive forecast" → "negative result + recovery via emergent-metric" | **PENDING** (this propagation cycle 2026-05-10) |
| **SENSITIVITY_BACK_OF_ENVELOPE update** | Recompute thresholds z β=-15/4 prior; update detection table | **PENDING** (this propagation cycle 2026-05-10) |
| **PPN_TO_PPE_MAPPING dictionary** | Explicit G_SPA=48 dictionary entry | **EXECUTED** (this v4 update) |
| **Phase 1.6 — Equal-mass DJS** | Full DJS 2-body Lagrangian z hyperbolic f(ψ); η-correction precision 5% | **DEFERRED** (multi-session; not blocking — η-correction ±25% sufficient dla 5σ verdict) |
| **M911-P2 re-derivation** | Multi-coefficient pattern z corrected G_SPA chain | **DEFERRED** (Phase 1.6+; M911-P2 currently WITHDRAWN-needs-rederivation) |
| **Recovery framework** | emergent-metric Phase 4 zero-β region exploration | **STRUCTURAL DERIVED** ([[../../research/op-emergent-metric-from-interaction-2026-05-09]]) |

### Native-first reframe (post-2026-05-10)

Per [[../../meta/PPN_AS_PROJECTION.md]] binding methodology:

| Layer | Status post-Phase 1.5 + RE-RUN |
|---|---|
| **L1 native** (g_tt Taylor c_n) | sympy-LOCK 5/5 PASS, structurally forced by α=2 + M9.1'' f(ψ) |
| **L2 chart** (β_ppE, b_ppE, G_SPA) | sympy-LOCK 5/5 PASS post-Phase 1.5; G_SPA=48 sympy-exact |
| **L3 falsifier** (GWTC-3 + ET-D + CE bounds) | GWTC-3 5.02σ FALSIFIED on f(ψ)=(4-3ψ)/ψ specific ansatz |

Native-first methodology emphasizes że falsification dotyczy **specific
M9.1'' f(ψ)=(4-3ψ)/ψ form**, NIE TGP framework całość. Recovery via
S07 alternative f(ψ) ansatz exploration + emergent-metric Phase 4
zero-β region (a_n, ξ_n, b_n, c_0, κ_σ) parameter space remains viable.

## Cross-references

- [[README.md]] — diagnoza T01, ścieżki A/B/C/D
- [[NEEDS.md]] — N2 (ppE mapping) jako konkret tego previewu
- [[FALSIFIER_STATEMENT_DRAFT.md]] — używa β_ppE^TGP placeholder
- [[SENSITIVITY_BACK_OF_ENVELOPE.md]] — używa b_ppE = −1 lub +1
- [[../../research/op-newton-momentum/M9_1_pp_P1_results.md]] §3.2 — pochodzenie liczb 5/6, 23/12, …
- [[../../research/op-newton-momentum/M9_1_pp_setup.md]] §2.5 — modyfikowane c²(ψ)
- [[../../PREDICTIONS_REGISTRY.md]] GW1, GW2, GW6, BH5 — orthogonal channels
- [[../EXTERNAL_REVIEW_2026-05-06.md]] §EXT-5 — recenzja źródłowa

## Bibliografia

- Yunes & Pretorius, Phys. Rev. D **80**, 122003 (2009), arXiv:0909.3328
- Cutler & Flanagan, Phys. Rev. D **49**, 2658 (1994)
- Will & Wiseman, Phys. Rev. D **54**, 4813 (1996)
- Mishra et al., Phys. Rev. D **93**, 084054 (2016) (3PN waveform)
- Damour, Jaranowski, Schäfer, Phys. Rev. D **89**, 064058 (2014) (4PN ADM)
- Yunes, Yagi, Pretorius, Phys. Rev. D **94**, 084002 (2016) (catalog)
- Maggiore et al., JCAP **03**, 050 (2020) (ET science case)
- Chamberlain & Yunes, Phys. Rev. D **96**, 084039 (2017) (3G implications)
