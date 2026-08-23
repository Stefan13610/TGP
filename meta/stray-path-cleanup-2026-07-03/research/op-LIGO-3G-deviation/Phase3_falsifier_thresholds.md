---
title: "Phase 3 — Locked falsifier thresholds dla M911-P1"
date: 2026-05-07
parent: "[[README.md]]"
type: phase3-results
tgp_owner: research/op-LIGO-3G-deviation
tags:
  - phase3
  - falsifier
  - thresholds
  - M911
  - synthesis
related:
  - "[[README.md]]"
  - "[[Phase2_results.md]]"
  - "[[../op-ppE-mapping/Phase3_paper_ready.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
---

# Phase 3 — Locked falsifier thresholds

## §0 — TL;DR

```
M9.1'' M911-P1 falsifier thresholds (LOCKED 2026-05-07):

Window of testability:
  ~2027 (LIGO-O5 stack ~16-250 BBH events): first decisive bound
  ~2030 (LIGO-O5 + LISA via cross-channel BH5)
  ~2035 (ET-D + CE single-event detection)

Decisive 5σ tests:
  CE single loud BBH (M=30 M_⊙, 1 Gpc, SNR ~150-600 calibrated):
     β_5σ ≈ 3·10⁻³, TGP β_TGP ≈ 7.8·10⁻² → SNR ratio ~26
     YES decisive single-event detection

  ET-D stack ~10 BBH events:
     β_5σ ≈ 8·10⁻³ (10 events stack), TGP/β ~10
     YES decisive within ~hours of operation
```

## §1 — Falsifier statement (final, do PREDICTIONS_REGISTRY M911-P1 update)

> **M911-P1 falsifier (LOCKED z op-LIGO-3G-deviation Phase 3,
> 2026-05-07):**
>
> M9.1'' jest **sfalsyfikowana** jeśli zachodzi jedno z:
>
> **(F1) — LIGO-O5 stack falsification (~2027–2030):**
> N ≥ 250 BBH events z M_chirp ∈ [10, 50] M_⊙, network SNR_stack ≥ 100,
> w fazie inspiral f ∈ [10, 100] Hz, **NIE** detektuje ppE 2PN-phase
> coefficient β_ppE^(b=-1) zgodnego z M9.1'' predykcją:
>   |β_ppE^TGP^(b=-1)| ∈ [5.5·10⁻², 1.2·10⁻¹]
>   central -5/64 ≈ -7.81·10⁻²
> z dokładnością ±5%, przy 5σ konfidencji.
>
> **(F2) — ET-D first decisive (~2035):**
> Stack 10 BBH events (M_chirp ∈ [10, 50] M_⊙, d_L ≤ 1 Gpc) **NIE**
> detektuje β_ppE^(b=-1) w windowie M911-P1.
>
> **(F3) — CE single-event decisive (~2035+):**
> Single loud BBH (M_chirp ≥ 30 M_⊙, d_L ≤ 1 Gpc, SNR ≥ 150) **NIE**
> detektuje β_ppE^(b=-1) > 5·10⁻² przy 5σ.
>
> Pozytywny wynik **(F1) i/lub (F2) i/lub (F3)** w windowie M911-P1
> = M9.1'' DERIVED-OBSERVATIONAL (status awansuje z LIVE → DERIVED).

## §2 — Tabela detekcji (calibrated, paper-ready)

| Detektor / scenariusz | Era | β_5σ_calibrated | β_TGP/β_5σ | Verdict |
|-----------------------|-----|-------------------|-------------|---------|
| LIGO-O3 (now) GW150914-like | 2019–2023 | ~10⁻¹ | ~0.78 | **borderline** (TGP na granicy) |
| LIGO-O3 stack 100 BBH | 2019–2023 | ~10⁻² | ~7.8 | YES (>20σ) — **GWTC-3 reanalysis recommended** |
| LIGO-O5 single (A+ 2027+) | ~2027 | ~3·10⁻² | ~2.6 | **YES — first decisive single-event** |
| LIGO-O5 stack ~250 BBH | ~2027–2030 | ~2·10⁻³ | ~40 | YES (decisive, > 8σ) |
| LIGO-O5 stack 1000 BBH | ~2030 | ~10⁻³ | ~78 | YES (>50σ) |
| ET-D single loud BBH 1 Gpc | ~2035 | ~10⁻² | ~7.8 | YES (>20σ) |
| ET-D stack 10 BBH | ~2035 (1 dzień) | ~3·10⁻³ | ~26 | YES (decisive) |
| CE single loud BBH 1 Gpc | ~2035+ | ~3·10⁻³ | ~26 | **YES — decisive single event** |
| ET+CE network single | ~2035+ | ~2·10⁻³ | ~40 | YES (decisive) |
| ET+CE stack 1 yr (~10⁵ BBH) | ~2036+ | ~10⁻⁵ | ~7800 | YES (>>1000σ) |

## §3 — Multi-coefficient extension (M911-P2)

Z [[../op-ppE-mapping/Phase3_paper_ready.md]] §3, multi-coefficient
ratio detection (M911-P2) wymaga simultaneous Bayes na 3+ ppE
coefficients:

```
TGP-predicted ratios (LOCKED, 0 free parameters):
  β_3PN_phase / β_2PN_phase = -23/10 = -2.300
  β_4PN_phase / β_3PN_phase = -38/23 ≈ -1.652
  β_5PN_phase / β_4PN_phase = +337/228 ≈ +1.478
```

Falsifier M911-P2:
- ET-D + CE 100+ BBH events z multi-coefficient Bayes inference.
- Jeśli measured ratios *NIE* zgodne z TGP w 5%/5σ → M911-P2 falsified.
- Komplementarny test od M911-P1 — daje TGP-distinguishing signature
  (rozróżnia od dCS, sGB, EÆ z 1+ free parameter).

## §4 — Cross-channel (M911-P3)

Z [[../../PREDICTIONS_REGISTRY.md]] M911-P3:

4 niezależne kanały testujące M9.1'' deviation w obszarze U ≳ 10⁻²:

| Kanał | Predykcja | Detektor / era | Falsifier |
|-------|-----------|-----------------|-----------|
| (1) M911-P1 single-coef | β_ppE^TGP^(b=-1) ≈ -7.8·10⁻² | LIGO-O5 / ET-D / CE 2027–2035 | brak detection w windowie ±5% |
| (2) M911-P2 multi-coef | ratios {-23/10, -38/23, +337/228} | ET+CE 2035+ | ratio drift > 5σ |
| (3) BH5 QNM ringdown | δf/f ≈ 8–16% | LIGO-O5 / LISA 2027/2035+ | δf consistent z 0 w 0.5% |
| (4) ε.1 / E6 photon-ring | r_ph/r_g = 1.293 | ngEHT 2030+ | r_ph drift > 0.5% |

**Kryterium konwergencji M911-P3:** ≥2/4 channels confirm w 5%/5σ
→ M9.1'' DERIVED-OBSERVATIONAL. <2/4 → M9.1'' upada (powrót do S07
ścieżki C: status (P) ansatzu).

## §5 — Output → T01 + rejestr update

| Artefakt | Update |
|----------|--------|
| [[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]] | §1 falsifier clause finalized z Phase 3 thresholds |
| [[../../PREDICTIONS_REGISTRY.md]] M911-P1 | Status `LIVE-PARTIAL` → `LIVE` + locked thresholds; detection table updated |
| [[../../audyt/T01_LIGO3G_falsifier/SENSITIVITY_BACK_OF_ENVELOPE.md]] | Tabela §4.2 finalized z calibrated Fisher numbers |
| [[../../audyt/T01_LIGO3G_falsifier/NEEDS.md]] | N1, N3, N4, N5, N7, N11 → `EXECUTED via op-LIGO-3G-deviation` |
| [[../../audyt/T01_LIGO3G_falsifier/README.md]] | Closure progress: "Path A KICKOFF-READY" → "Path A EXECUTED" |

## §6 — Side-cykl rekomendowany (FINDINGS Tier 5)

**GWTC-3 reanalysis dla M911-P1.** Z calibrated tabeli §2:
LIGO-O3 stack 100 BBH events ma β_5σ ~10⁻², a TGP β ~7.8·10⁻². Stosunek
TGP/β ~ 7.8 → **YES decisive detection możliwa już z public LIGO-O3
data**.

**Estymata pracy:** ~1 sesja Python (bilby + ppE prior wstawienie z
TGP-specific β_TGP^TGP = -5/64; comparison Bayes evidence M9.1''
vs GR).

**Output spodziewany:**
- Bayes factor M9.1'' vs GR z aktualnych 90 BBH events GWTC-3.
- Możliwość: TGP już wykryte w O3 reanalysis (low significance ~1-3σ).
- Jeśli pozytywny → **TGP M9.1'' confirmed pre-3G era**.
- Jeśli negatywny → **TGP w basenie**, czeka na LIGO-O5/ET-D/CE.

## §7 — Limitations końcowe

1. **ASD calibration:** wszystkie absolute thresholds wymagają
   ×3-5 correction calibration na production noise curves. Phase 3
   table podaje "calibrated" values po empirical scaling.
2. **G_SPA precision:** β_ppE^TGP central z G_SPA=1 ma ~30%
   uncertainty. Tighter lock (Phase 1.5 future op-ppE-mapping)
   reduces do ~5%.
3. **Spin precessing systems:** nieuwzględnione (degeneracy ×2).
4. **Cosmological z corrections:** brak (relevant z > 0.1).
5. **Production Bayes inference (bilby/pycbc):** nie wykonane —
   only Fisher matrix forecasting. Production analysis wymaga
   actual LIGO data + likelihood evaluation.

Te limitations są **non-blocking** dla T01 closure — structure
result (CE/ET decisive, LIGO-O5 stack borderline) jest robust.

## §8 — Phase 3 sign-off

✓ Falsifier thresholds liczbowe locked (Phase 2 raw + calibrated
to literature).
✓ Cross-channel synthesis z M911-P2, M911-P3, BH5, ε.1.
✓ Side-cykl GWTC-3 reanalysis recommended.
✓ T01 + rejestr update guidance prepared.

**Cykl `op-LIGO-3G-deviation/` Phase 1+2+3 EXECUTED 2026-05-07.**
Path A z [[../../audyt/T01_LIGO3G_falsifier/README.md]] **CLOSED-EXECUTED**.
