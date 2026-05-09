---
title: "PPN → ppE mapping (preview Path B): (5/6) U³ deviation → β_ppE^(3PN)"
date: 2026-05-07
parent: "[[README.md]]"
type: analytical-preview
tgp_owner: audyt/T01_LIGO3G_falsifier
tags:
  - preview
  - analytical
  - ppE
  - PPN
  - 3PN
  - M911
  - inspiral
  - phase-modification
  - T01
  - EXT-5
related:
  - "[[README.md]]"
  - "[[NEEDS.md]]"
  - "[[FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[SENSITIVITY_BACK_OF_ENVELOPE.md]]"
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
---

# PPN → ppE mapping — preview Path B

> **Cel pliku.** Dostarczyć **analityczny szkielet** mapowania
> (5/6) U³ deviation w `g_tt^TGP` na parametry ppE Yunes–Pretorius
> (β_ppE, b_ppE) dla fazy waveformu inspiralu BBH. Ten plik jest
> **previewem** dla cyklu `research/op-ppE-mapping/` Phase 1; nie
> zastępuje go. Dostarcza:
> - dictionary "metric coefficient → SPA phase coefficient",
> - identyfikację b_ppE^TGP,
> - **structural** (nie precise) szacowanie β_ppE^TGP ze skali (5/6),
> - listę założeń, które cykl badawczy musi zwalidować.
>
> **Status:** SCHEMATYCZNY. Liczby są przybliżone z dokładnością
> O(1) prefactor. Cykl `op-ppE-mapping/` zamknie precyzję do %.

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

### 2.2 Heurystyka magnitude

Dla M9.1'':

```
Δg_tt/c² = − (5/6) U³ + …      gdzie U ≈ v²/c² ≈ u²/c²
```

W phase modification (heurystyczna SPA):

```
δΨ_TGP(u) ≈ (3/(128 η)) · u⁻⁵ · κ_TGP · u^(2N_phase)
```

gdzie N_phase = 2 (2PN phase) i κ_TGP zawiera (5/6) plus prefactory
z dE/dt i E_orb mappings. Konkretnie (heurystycznie, numerycznie):

```
κ_TGP ≈ (5/6) · (numerical prefactor O(1))
       ≈ 0.5 ÷ 1.5  (one-sigma; precise value to be locked in op-ppE-mapping)
```

Wtedy:

```
δΨ_TGP(u) ≈ (3/(128 η)) · κ_TGP · u⁻³        ← 1PN-like u⁻³ form
```

**Stop. Dla precise:** powyższy heurystyczny scaling NIE jest dokładny
— pełne mapowanie wymaga (a) M9.1'' two-body Lagrangian (DJS-like,
ale z hiperboliczną metryką), (b) dE/dt z multipole expansion w
M9.1'' (UWAGA: potencjalna modyfikacja luminosity przez modyfikację
quadrupole formula w M9.1''!). Cykl `op-ppE-mapping/` Phase 1 to
zamknie. Tutaj dostajemy *strukturę* i **expectation: β_ppE^TGP
to liczba O(1) razy (5/6), w b_ppE = −1 lub b_ppE = +1 zależnie
od konwencji PN counting**.

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

**Dictionary entry (proposed for `op-ppE-mapping/` Phase 1):**

```
M9.1'' → ppE:
  b_ppE^TGP = −1 (2PN phase, U³ metric)
  β_ppE^TGP = numerically_locked_from_E_orb_2body_M911
  cross-pattern: {β_(2PN), β_(3PN), β_(4PN)} = κ_TGP · {(5/6), (23/12)·c_4ratio, (19/6)·c_5ratio}
                 where c_Nratio = SPA chain prefactor at PN order N
                 (to be locked from 2-body Lagrangian)
```

## §5 — Założenia, które cykl `op-ppE-mapping/` musi zwalidować

| ID | Założenie | Status |
|----|-----------|--------|
| A1 | M9.1'' two-body Lagrangian istnieje i daje finite Lagrangian (nie blow-up przy ψ → 4/3 horizon w external potencjale) | OPEN — wymaga sek08c uogólnienia na dwa źródła |
| A2 | Quadrupole formula w M9.1'' modyfikuje się w sposób reprezentatywny przez metric perturbation tylko (NIE nowymi dipole / scalar charges) | OPEN — wymaga gauge-invariant analiza |
| A3 | dE/dt luminosity w M9.1'' jest *modyfikowana* przy 2PN+ orbital w sposób spójny z metric coefficients (5/6) U³ | OPEN |
| A4 | Stationary phase approximation jest ważna w M9.1'' z tego samego powodu co w GR (adiabatic inspiral) | LIKELY (no obvious mechanism breaks adiabaticity) |
| A5 | TGP nie wprowadza nowych radiacyjnych DOF (no scalar mode) zgodnie z GW1 (c_T = c_s) i GW2 (3 DOF) | LIVE w PREDICTIONS_REGISTRY |
| A6 | Konwencja PN counting "metric N-PN ↔ phase (N−1)-PN" jest poprawna dla M9.1'' (z hyperbolic ansatz) | OPEN — może być inna |

**Bez zwalidowania A1–A6 powyższe szacowania są PRELIMINARY.** Audit
folder T01 *nie zamyka* tych założeń — przekazuje je do
`op-ppE-mapping/` jako contract.

## §6 — Następny krok dla cyklu `op-ppE-mapping/` Phase 1

1. **Setup (Phase 1.1):** Two-body Lagrangian w M9.1'' z DJS approach
   uogólnionym na hyperbolic f(ψ). Reference: Damour–Jaranowski–Schäfer
   2014 dla GR; uogólnienie na M9.1'' wymaga substitution
   c² → c²·(4-3ψ)²/ψ² (z M9_1_pp_setup §2.5).

2. **Phase 1.2:** Dla equal-mass binary η = 1/4, wyprowadzić
   E_orb(v) do v⁸ (4PN orbital). Verify że E_orb^TGP odbiega od
   E_orb^GR przy v⁴ (2PN orbital), zgodnie z metric (5/6)U³.

3. **Phase 1.3:** Wyprowadzić dE/dt do v¹⁸ (4PN luminosity). Założenie
   A3 zwalidować lub modyfikować quadrupole formula.

4. **Phase 1.4:** SPA inversion → Ψ_TGP(f). Zidentyfikować
   β_ppE^(2PN_phase), β_ppE^(3PN_phase), β_ppE^(4PN_phase) liczbowo
   dla η=1/4 i równolegle dla η ∈ [0.1, 0.25].

5. **Phase 1.5:** Compare z catalog Yunes–Yagi–Pretorius 2016. Ustal
   uniqueness multi-coefficient pattern jako TGP-specific signature.

**Output Phase 1:** liczbowy [β_th] dla [[FALSIFIER_STATEMENT_DRAFT.md]]
§1 placeholder.

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
