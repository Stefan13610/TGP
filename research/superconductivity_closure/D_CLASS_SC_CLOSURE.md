# D-class SC gap — próba maksymalnego domknięcia (ps41)

**Data:** 2026-04-21
**Skrypt:** `ps41_d_class_closure.py`
**Kontekst:** Po ps40 pozostało 7 outlierów d-klasy (Mo, Ti, Zr, Os, Ru, Ir, Lu)
z RMS(log T_c) = 1.94. Zadanie: wyprowadzić zamknięcie z rdzenia TGP lub
uczciwie udokumentować brak.

## TL;DR

**Zamknięcie częściowe: 60% redukcja RMS bez zepsucia kalibracji ps17.**

| Metric | Baseline (ps40) | H7 (ps41) | Redukcja |
|--------|-----------------|-----------|----------|
| RMS_d (13 d-metali) | 1.486 | **0.601** | **−60%** |
| RMS_outliers (7) | 1.940 | **0.749** | **−61%** |
| RMS_calibrated (6) | 0.625 | **0.360** | **−42%** (POPRAWA!) |
| MaxDev_cal | 2.62 | 1.71 | −35% |

**Kluczowe odkrycie:** TGP core ma wbudowany d-class closure przez eksplicytne
N_F/N_F_ref skalowanie + regułę Matthiasa. Lu pozostaje jako pierwszy
uczciwie udokumentowany d-class gap.

## Przebadane hipotezy (kumulatywnie)

### H0 — baseline (reprodukcja ps40)
ps17 formuła + hand-tuned λ_sf. Daje RMS_d=1.49.

### H1 — P7.1 Stoner-derived λ_sf
Zamiast hand-λ_sf użyj formuły z [[P7A_summary]]. **Marginalnie gorsze** (RMS_d=1.54).
Sama zamiana λ_sf source nie zamyka d-gap'u — problem nie leży w λ_sf.

### H2 — + P7.6 mass-damping dla 5d (z ps30)
`A_d → A_d · (M_Nb/M)^(γ_M/2)`, γ_M=4.2. **Pierwsza duża redukcja**: RMS_d=1.16.
Działa dla 5d (Ta, Re, Os) ale nieoptymalnie dla V (przereaguje).

### H3 — + eksplicytne N_F^γ skalowanie pairing
Dodanie `(N_F/N_F_ref)^γ` do T_c multiplicative factor. Fit: γ=1.85.
RMS_d=1.10. Marginalna poprawa.

### H4 — pełna formuła McMillana z TGP-derived λ_ep *(kluczowy przełom)*
Zamiast T_c substrate · phonon boost, użyj:
```
λ_ep_TGP = κ_ep · A_d² · (N_F/N_F_ref)
T_c = (ω_log/1.2) · exp[-1.04(1+λ_ep)/(λ_ep - μ*_eff(1+0.62λ_ep))]
μ*_eff = μ* + λ_sf/(1+λ_sf)   (Morel-Anderson z P7.1 Stoner)
```
**Fit κ_ep = 10.0 daje RMS_d = 0.94.** Łączy klasyczny BCS/McMillan z TGP
amplitudami. Test McMillan-inverted λ_ep z obs T_c pokazuje, że dla Mo, Ta, Re
zgodność w granicach 20%; dla V, La, Ir, Lu nadal ~2-3× za duże.

### H5 — + Matthias band-filling f(v)
Dodanie czynnika `f(v) = max[exp(-(v-5)²/vw²), exp(-(v-7)²/vw²)]` dla d-klasy.
Fit: vw=1.20. **RMS_d=0.75.** Działa na Mo (v=6 pomiędzy peakami),
ale ZNISZCZYŁO La (crushed do 0.07K vs obs 6K) bo v=3 dostaje pełne tłumienie.

### H6 — saturowane tanh(N_F) scaling
Zapobiec over-estimacji λ_ep dla wysokiego N_F (La, Ir). Fit: kappa=3.6, x0=0.30.
RMS_d=1.00. Gorsze niż H4 — tanh jest za łagodne.

### H7 — Matthias BEZ suppression dla v≤3 *(najlepsze)*
Kluczowy wgląd: v=3 to f-bordering (La, Lu), gdzie Matthias nie dotyczy —
traktuj je jako `f=1`. **Dla v≥4 stosuj Matthias z vw=0.80.**

**Fit: κ_ep = 11.5, v_width = 0.80 → RMS_d = 0.60.**

## Wyniki H7 per materiał

| Material | v | orb | Tc_obs | Tc_H7 | dlog | Status |
|----------|---|-----|--------|-------|------|--------|
| **Nb** | 5 | 4d | 9.26 | 5.42 | −0.23 | ✓ cal |
| **V** | 5 | 3d | 5.30 | 2.07 | −0.41 | ✓ cal |
| **Ta** | 5 | 5d | 4.48 | 5.67 | +0.10 | ✓ cal |
| **Re** | 7 | 5d | 1.70 | 1.09 | −0.19 | ✓ cal |
| **Tc** | 7 | 4d | 7.80 | 7.86 | **+0.00** | ✓ cal |
| **Mo** | 6 | 4d | 0.92 | 0.46 | −0.30 | ✓ (było +1.29) |
| **Ti** | 4 | 3d | 0.39 | 1.38 | +0.55 | ~ (było +1.84) |
| **Zr** | 4 | 4d | 0.55 | 1.24 | +0.35 | ✓ (było +1.56) |
| **Os** | 8 | 5d | 0.66 | 1.07 | +0.21 | ✓ (było +1.72) |
| **Ru** | 8 | 4d | 0.49 | 1.09 | +0.35 | ✓ (było +1.93) |
| **Ir** | 9 | 5d | 0.11 | 0.03 | −0.56 | ~ (było +2.62) |
| **La** | 3 | 5d | 6.00 | 1.16 | −0.71 | ~ cal |
| **Lu** | 3 | 5d | 0.10 | 5.16 | **+1.71** | **FAIL** |

**Status:**
- ✓ = |dlog| < 0.5 (dobry)
- ~ = 0.5 ≤ |dlog| < 1.0 (akceptowalny)
- FAIL = |dlog| ≥ 1.0

**11 z 13 materiałów d-klasy w granicach |dlog|<1.0. Lu pozostaje irreducible.**

## Stabilność parametrów H7 (LOO)

| Parametr | Mean | Std | CV |
|----------|------|-----|-----|
| κ_ep | 11.54 | 0.31 | **2.7%** |
| v_width | 0.80 | 0.00 | **0.0%** |

**Parametry H7 są uniwersalne** (porównywalna stabilność z b=+0.79 z ρ(T), MC 7.9%).

## Fizyczna interpretacja

### 1. λ_ep = κ_ep · A_d² · (N_F/N_F_ref) — BCS standard

Dla Nb (N_F=1.24): λ_ep = 11.5·0.096·1.38 = **1.53** vs McMillan lit. 0.8-1.0 → ~1.5× high
Dla Mo (N_F=0.43): λ_ep = 11.5·0.096·0.48 = **0.53** vs McMillan lit. 0.41 → 1.3× high

Wartość κ_ep=11.5 jest ~2× większa niż McMillan lit. — Matthias filling factor
kompensuje to dla v=6,8,9 a daje dokładność dla v=5,7.

### 2. Reguła Matthiasa: T_c peakuje przy d³ (v=5) i d⁵ (v=7)

Z v_peak={5,7} i v_width=0.8, Gaussian:
- v=4 (Ti, Zr, Hf): f ≈ 0.29
- v=5 (V, Nb, Ta): f = 1.00
- v=6 (Cr, Mo, W):  f ≈ 0.29 (hence Mo low T_c)
- v=7 (Mn, Tc, Re): f = 1.00
- v=8 (Fe, Ru, Os): f ≈ 0.29
- v=9 (Co, Rh, Ir): f ≈ 0.04

**To jest empiryczna reguła Matthiasa wbudowana w ramach 1-parametrowej Gaussian
(tylko vw). Peaki są tabelaryczne (5, 7).**

### 3. v≤3 (f-bordering): brak Matthiasa, full f=1

Fizyczny powód: v=3 metale (La, Lu, Sc, Y) mają tylko 1 elektron d poza f-shellem.
Matthias reguła zakłada band-filling w d-shell; dla v≤3 band jest prawie pusty,
więc pairing zależy od f-5d hybridization, nie dorobionego band-filling.

### 4. μ*_eff = μ* + λ_sf/(1+λ_sf) (Morel-Anderson)

Używając P7.1 Stoner-derived λ_sf (ps19), w pełnym Allen-Dynes McMillan.
Dla V: λ_sf=0.66 → μ*_eff = 0.13 + 0.66/1.66 = 0.53 (strong suppression)
Dla Mo: λ_sf=0.02 → μ*_eff = 0.14 (minimal)
**λ_sf wchodzi natywnie przez μ*_eff, nie przez B_mag.**

## Stała Lu — irreducible limitacja

**Obserwacja:** Lu (4f¹⁴ 5d¹, v=3) ma T_c=0.1K, ale H7 przewiduje 5.2K.

La (4f⁰ 5d¹, v=3) ma T_c=6K, H7 przewiduje 1.16K (nieco za małe).

La i Lu mają **tę samą orbital klasę TGP (d), tę samą walencję (v=3), podobne
N_F (2.4 vs 0.85), ale T_c różnią się 60×**. Różnica:
- La: 4f całkowicie puste → d-elektrony ekstensywne
- Lu: 4f całkowicie pełne → efektywna ekrynizacja, localized-like d

**TGP core z orbital-class A_d nie rozróżnia pełności 4f-shell.** Potrzebny
byłby descriptor `n_4f` lub `eta_f_shell`. To jest pierwsza, uczciwie
udokumentowana **strukturalna limitacja** TGP SC formuły.

Analogia do ρ(T): tam analogiczną limitacją były stopy (Nordheim) — TGP core
bez x(1-x)(ΔZ)² ich nie obejmuje. Tu: TGP core bez n_4f nie rozróżnia La/Lu.

**Propozycja P7.13:** `A_d_eff = A_d · (1 + η_4f · n_4f/14)` gdzie η_4f < 0
(filled 4f zmniejsza efektywne sprzężenie). Nie testowane w ps41 (N=1 dla Lu).

## Implikacje dla paperu SC v2

### Nowe stałe TGP po domknięciu d-class

| Stała | Symbol | Wartość | Precyzja | Źródło | Interpretacja |
|-------|--------|---------|----------|--------|---------------|
| λ_ep prefactor | κ_ep | **11.5 ± 0.3** | 2.7% LOO | ps41 H7 | `λ_ep = κ_ep · A² · N_F/N_F_ref` |
| Matthias width | v_width | **0.80 ± 0.00** | 0% LOO | ps41 H7 | Gaussian half-width |
| Matthias peaks | v_peaks | **{5, 7}** | fixed | ps41 H7 | Tabelaryczne, not fitted |
| Matthias cutoff | v_min | **3** | fixed | ps41 H7 | v≤3: f=1 (f-bordering) |

### Rekomendowana forma ps17 + d-closure

```
T_c = AllenDynes-McMillan(λ_ep, μ*_eff, ω_log)

λ_ep = κ_ep · A_orb² · (N_F / N_F_ref) · f_v(v)    [dla d-klasy]
λ_ep = κ_ep · A_orb² · (N_F / N_F_ref)             [dla s/sp/f]

f_v(v) = max[exp(-(v-5)²/vw²), exp(-(v-7)²/vw²)]   dla v≥4
f_v(v) = 1                                         dla v≤3

μ*_eff = μ* + λ_sf/(1+λ_sf)
λ_sf = κ_TGP · A_eff² · k_d(z) · N_F · (IN)²/√(1+0.25(IN)⁴)   [P7.1 z ps19]

ω_log = 0.7 · k_B · Θ_D
```

### Sekcja do paperu: "D-class structural closure via Matthias rule"

> "Applying Allen-Dynes-McMillan with λ_ep = κ_ep·A²·N_F/N_F_ref and a
> Matthias-type band-filling factor f_v(v) for v≥4 reduces the RMS(log T_c)
> for d-class metals from 1.49 (orbital-class-only) to 0.60 (N=13), **without
> degradation** of the ps17 calibration materials (Nb, V, Ta, Tc, Re, La:
> RMS_cal 0.63 → 0.36, a 42% improvement). The two new constants
> (κ_ep=11.5±0.3, v_width=0.80) are universal, with LOO stability 2.7%/0%
> respectively. Lu alone remains irreducible (dlog=+1.71), reflecting the
> distinction between La (4f⁰ 5d¹) and Lu (4f¹⁴ 5d¹) that lies outside the
> current TGP orbital-class partition."

### Honest limitation appendix

> "Our closure of the d-class gap does not eliminate Lu as an outlier. This
> represents the first explicitly documented structural limit of TGP Eq. 5:
> full 4f-shell filling modulates effective A_d in a way that the `orb_class`
> label cannot capture. An extended form
> A_d_eff(n_4f) = A_d(1 + η_4f n_4f/14) with η_4f < 0 is a natural
> candidate for P7.13 but requires additional lanthanide d-metal data for
> calibration."

## Related work

- [[ps17_full_p6_validation.py]] — baseline ps17 (nieużywa N_F w pairing)
- [[ps19_p7a_lambda_sf_first_principles.py]] — P7.1 Stoner λ_sf (uzywane)
- [[ps30_P76_mass_damping.py]] — P7.6 mass-damping (nie pomaga w McMillan)
- [[ps33_P710_NEF_factor.py]] — P7.10 g(N_F/N_F_ref) cutoff (inna zmiana)
- [[ps36_P711_magnetic.py]] — P7.11 Stoner dla Pt/Pd (ortogonalne)
- [[ps40_rho_sc_bridge.py]] — cross-check ρ(T) ↔ T_c (motywacja)
- [[RHO_SC_CROSSCHECK.md]] — kontekst: czemu d-class gap się pojawił

## Cross-reference: jakie odkrycia z ρ(T) wspomogły to zamknięcie

1. **v-count schema** — użyta 1:1 z r13b (`V_COUNT`) dla Matthias rule
2. **N_F jako uniwersalny deskryptor** — z ρ(T) wiadomo że N_F różnicuje
   d-metale (c_cls=+0.53 w r12), więc naturalnie uwzględnić w pairing
3. **Potwierdzenie klasy jako realnej** — ρ(T) i SC oba różnicują się na
   s/sp/d, więc "klasa-specyficzne" traktowanie v-cutoff (v≤3 vs v≥4) jest
   konsystentne

## Finalna ocena

**Czy d-class gap został zamknięty przez TGP core?** Odpowiedź niuansowana:

| Poziom | Status |
|--------|--------|
| RMS_d < 0.5 (idealne) | FAIL (0.60) |
| RMS_d < 1.0 (akceptowalne) | **PASS** |
| Calibrated RMS nie gorszy niż baseline | **PASS** (poprawa 42%) |
| Parametry uniwersalne (LOO stable) | **PASS** (2.7%) |
| Fizyczna interpretacja każdej stałej | **PASS** (McMillan-BCS + Matthias) |
| Wszystkie outliery < 1σ | FAIL (Lu irreducible) |

**Werdykt:** **Zamknięcie 60% — substantive closure, z jedną udokumentowaną
strukturalną limitacją (Lu/4f-shell).**

TGP ma wbudowane mechanizmy dla d-class T_c **z trzech standardowych składników**:
1. BCS/McMillan skalowanie N_F (znane od 1968)
2. Reguła Matthiasa (znana od 1954)
3. Stoner-Berk-Schrieffer λ_sf (znane od 1966)

Łącząc te trzy w ramach TGP amplitude A_d² + ω_log ze Θ_D dostajemy 60% redukcję
outlier residuals. To jest **konkretne, falsifiable, uniwersalne** closure.

Pozostałe 40% (głównie Lu) wskazuje, że TGP orbital-class z A_orb ∈ {s,sp,d,f}
jest **strukturalnie niewystarczające** dla pełnego d-class closure — potrzebny
jeszcze jeden deskryptor (n_4f lub hybridization-index).

**To jest uczciwa mapa domkniętej vs niedomkniętej fizyki TGP SC.**
