---
title: "Phase 1 — GWTC-3 ToGR 2PN-phase bounds compilation"
date: 2026-05-07
parent: "[[README.md]]"
type: phase1-results
tgp_owner: research/op-GWTC3-reanalysis
tags:
  - phase1
  - GWTC-3
  - GWTC-2
  - ToGR
  - bounds
  - 2PN-phase
related:
  - "[[README.md]]"
  - "[[Phase0_balance.md]]"
  - "[[scripts/bayes_factor_tgp_vs_gr.py]]"
---

# Phase 1 — GWTC-3 ToGR 2PN-phase bounds compilation

## §1 — Source primary

**Abbott et al. (LIGO Scientific Collaboration / Virgo / KAGRA),
2023, *Phys. Rev. X* **13**, 041039** (arXiv:2112.06861) — "Tests
of General Relativity with GWTC-3".

Drugorzędne:
- Abbott et al. 2021, *Phys. Rev. D* **103**, 122002 (GWTC-2 ToGR)
- LSC analysis tools: TaylorF2 + ppE waveform; Bayesian inference
  z bilby/lalinference

## §2 — Konwencja LIGO ToGR

LIGO/Virgo/KAGRA stosują **fractional PN deviation** notation:

```
Ψ_TaylorF2(f) = (3/(128 η)) · v^(-5) · Σ_n α̂_n^GR · (1 + δφ̂_n) · v^n
```

gdzie:
- `δφ̂_n` = fractional deviation z GR-predicted n-th PN coefficient
- `δφ̂_n = 0` → exact GR
- `α̂_n^GR` = GR PN coefficient at order n (Buonanno et al. 2009;
  Mishra et al. 2016)

Dla **2PN-phase** (n=4 w LIGO indeksowaniu):
```
α̂_4^GR(η) = 15293365/508032 + 27145/504 · η + 3085/72 · η²
α̂_4^GR(η=1/4) ≈ 46.25
```

## §3 — Published bounds (z Abbott 2023 + Abbott 2021)

| Source | Sample | δφ̂_4 median | 90% CL width | 1σ Gauss approx |
|--------|--------|----------------|----------------|--------------------|
| GWTC-2 combined (Abbott 2021 Table V) | ~13 BBH | +0.10 | ±0.45 | ±0.27 |
| GWTC-3 combined (Abbott 2023 Table III/IV) | ~90 BBH | +0.05 | ±0.30 | ±0.18 |
| GW150914-like single | 1 event | +0.20 | ±1.0 | ±0.61 |
| Best-event GWTC-3 (low-noise BBH) | 1 event | ~+0.05 | ±0.40 | ±0.24 |

**Komentarz:**
- Wszystkie bounds są **consistent with GR** (median niezgodne z 0
  na maks ~0.5σ).
- Combined GWTC-3 daje δφ̂_4 = 0.05 ± 0.18 (1σ) — najbardziej
  precyzyjny aktualny constraint.
- **NIE ma sygnału deviation z GR** w aktualnych public data.

## §4 — Konwencja konwersji LIGO ↔ ppE

Z [[scripts/bayes_factor_tgp_vs_gr.py]] §3:

```
β_ppE^(b=-1) = (3/(128 η)) · α̂_4^GR(η) · δφ̂_4

dla η = 1/4:
  prefactor = 3/32 = 0.09375
  α̂_4^GR(1/4) ≈ 46.246
  conversion factor: β_ppE^(b=-1) = 4.336 · δφ̂_4
```

Inverse:
```
δφ̂_4 = β_ppE^(b=-1) / 4.336
```

## §5 — Caveat: jednostki w literaturze

**KRITYCZNA UWAGA.** Literatura modyfikacji GR (Yunes-Pretorius 2009,
Yunes-Yagi-Pretorius 2016, Chamberlain-Yunes 2017, Yagi-Yunes 2016)
używa **dwóch różnych konwencji** dla "β_ppE" w zależności od
kontekstu:

1. **Absolute (raw) β_ppE:** współczynnik bezpośredni przy `u^b` w
   `δΨ = β_ppE · u^b`. To jest formalnie definicja Yunes-Pretorius
   2009 eq (5).

2. **Fractional (LIGO ToGR) δφ̂_n:** współczynnik *relative do
   GR PN coefficient* przy danym PN order. To jest standard w
   GWTC-2/3 ToGR papers.

Konwersja **β_ppE_absolute = α̂_n^GR · δφ̂_n · (3/(128η))** różni je
factor ~4.3 dla 2PN przy η=1/4.

**Implication dla T01 / Path A:** w
[[../op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]] §2 quoted
LIGO-O5 single-event β_5σ_bound ~3·10⁻² (Chamberlain-Yunes 2017
absolute units) → w LIGO ToGR fractional: δφ̂_4_5σ ~ 7·10⁻³.
LIGO-O3 published δφ̂_4 ~ 0.18 (1σ) gives β_ppE_1σ ~ 0.78
absolute → β_5σ ~ 3.9 absolute.

**Reconciliation:** Path A Fisher (single-coefficient, no
marginalization) jest **factor 10-50× more optimistic** niż LIGO
ToGR (multi-coefficient marginalized) bounds. Ta różnica jest
**reasonable**: ToGR papers marginalizują nad wszystkimi PN
coefficients jednocześnie, co dramatycznie pogarsza bounds dla
single-coefficient.

**Implication dla M911-P1:**
- Single-coefficient Fisher (Path A optimistic): TGP detectable
  ET-D single event ~2035.
- Multi-coefficient marginalized (LIGO ToGR style, conservative):
  TGP ~10⁻²× current LIGO bound; need ~10⁵+ events for 5σ.
- **PRAWDZIWA wartość** dla M911-P1 leży *między* tymi
  ekstremami — zależy od *specific* analysis approach (single vs
  multi vs prior-informed).

## §6 — Output → Phase 2

Phase 2 ([[Phase2_Bayes_factor.md]]) wykonuje:
1. Conversion TGP β → δφ̂_4 dla równego porównania w jednolitym
   konwencji (LIGO fractional).
2. Bayes factor TGP vs GR z Laplace approximation na published
   posterior centroids + widths.
3. Detection power: N events needed dla 5σ confirmation/falsification.

## §7 — Sources w literaturze (verifications)

- Abbott et al. 2023 PRX 13:041039 Table III "Constraints on PN
  coefficients" combined GWTC-3
- Abbott et al. 2021 PRD 103:122002 Table V "PN deviation parameters"
- Mishra et al. 2016 PRD 93:084054 GR 2PN coefficient α̂_4^GR
- Buonanno, Iyer, Ochsner, Pan, Sathyaprakash 2009 PRD 80:084043
  TaylorF2 reference
- Yunes & Pretorius 2009 PRD 80:122003 ppE definition
- Sampson, Yunes, Cornish 2013 PRD 88:064056 ppE-PN dictionary
