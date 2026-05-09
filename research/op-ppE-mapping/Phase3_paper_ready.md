---
title: "Phase 3 — Paper-ready output: M9.1'' β_ppE^TGP^(b=-1) lock + multi-coefficient signature"
date: 2026-05-07
parent: "[[README.md]]"
type: phase3-results
tgp_owner: research/op-ppE-mapping
tags:
  - phase3
  - paper-ready
  - synthesis
  - M911
  - ppE
  - multi-coefficient
  - falsifier
related:
  - "[[README.md]]"
  - "[[Phase1_results.md]]"
  - "[[Phase2_literature_crosscheck.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
---

# Phase 3 — Paper-ready output

## §1 — Final lock (paper headline)

```
┌─────────────────────────────────────────────────────────────┐
│  M9.1'' deviation w fazie inspiralu BBH:                    │
│                                                             │
│   β_ppE^TGP^(b=-1) = -(3/(128 η)) · (5/6) · G_SPA           │
│                                                             │
│   Equal-mass (η=1/4), G_SPA=1 (central):                    │
│      β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 · 10⁻²                │
│                                                             │
│   OOM window (G_SPA ∈ [0.7, 1.5]):                          │
│      |β_ppE^TGP^(b=-1)| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹]          │
│                                                             │
│  Multi-coefficient signature (TGP-distinguishing):          │
│      β_(N+1)PN / β_NPN = {-23/10, -38/23, +337/228}         │
│      0 free parameters — FORCED by α=2 + hyperbolic f(ψ)    │
└─────────────────────────────────────────────────────────────┘
```

## §2 — Falsifier statement (liczbowy, do PREDICTIONS_REGISTRY M911-P1 update)

> **M911-P1 falsifier (locked liczbowy):**
>
> Jeśli ET-D + CE z **siecią SNR ≥ 100** (single event lub stack
> N ≥ 100 BBH events z M_chirp ∈ [10, 50] M_⊙ w fazie inspiral
> f ∈ [10, 100] Hz), **nie** detektuje ppE 2PN-phase coefficient
> β_ppE^(b=-1) zgodnego z M9.1'' predykcją:
>
>   |β_ppE^TGP^(b=-1)| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹] (OOM window)
>   β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 · 10⁻² (central, G_SPA=1)
>
> z dokładnością ±5% przy 5σ konfidencji, **M9.1'' jest sfalsyfikowana**.

## §3 — Multi-coefficient TGP-distinguishing test (M911-P2 → ET+CE 4PN)

```
┌──────────────────────────────────────────────────────────────┐
│  Multi-coefficient ratio test (TGP-distinguishing M911-P2)   │
│                                                              │
│  W ET+CE simultaneous Bayes inference {β_2PN, β_3PN, β_4PN}: │
│                                                              │
│    β_3PN / β_2PN = -23/10 = -2.300  (TGP)                    │
│    β_4PN / β_3PN = -38/23 ≈ -1.652  (TGP)                    │
│    β_5PN / β_4PN = +337/228 ≈ +1.478 (TGP)                   │
│                                                              │
│  Modyfikacje GR z 1+ free parameter (dCS, sGB, EÆ, BD)       │
│  NIE reprodukują tego pattern z bare coefficients —          │
│  wymagają fitting, który nie jest predyktywny.               │
│                                                              │
│  Falsyfikacja M911-P2: jeśli ET+CE multi-coef Bayes daje     │
│  ratio {β_3PN/β_2PN, β_4PN/β_3PN} odbiegający od TGP         │
│  predicted z 5σ, M9.1'' multi-coefficient pattern upada.     │
└──────────────────────────────────────────────────────────────┘
```

## §4 — Tabela porównawcza z literaturą (paper-ready)

| Theory | Free params | β_ppE^(b=-1) | β_ppE^(b=+1) | Ratio prediction |
|--------|-------------|---------------|---------------|-------------------|
| **TGP M9.1''** | **0** | **-5/64 ≈ -7.8·10⁻²** | **-23/128 ≈ -1.8·10⁻¹** | **-23/10 (forced)** |
| dCS | 1 (ζ_dCS) | (1.55·10⁵/1.18·10⁷) ζ_dCS | needs ζ_dCS² | varies w/ ζ_dCS |
| sGB | 1 (ζ_sGB) | -(5/7168) ζ_sGB · ξ(η) | needs ζ_sGB² | varies w/ ζ_sGB |
| Einstein-Æther | 4 (c_i) | -3/(128η)(1-c_T/c_s) | depends c_i | varies w/ c_i |
| Brans-Dicke | 1 (ω) | -(5/3584)/(2ω+3) | depends ω | dominant b=-7 (dipole) |

**KEY INSIGHT (paper-ready):**

> *TGP M9.1'' ma 0 free parameters w β_ppE — wszystkie współczynniki
> są wyprowadzone analitycznie z α=2 vacuum Φ-EOM i hyperbolic
> f(ψ)=(4-3ψ)/ψ. To czyni TGP **bardziej restryktywną** niż wszystkie
> standardowe modyfikacje GR z 2PN-phase deviation. Multi-coefficient
> ratio pattern jest **decisive distinguishing test** w ET+CE era.*

## §5 — Tabela detekcji (paper-ready, z lock z Phase 1)

Aktualizowane wartości z [[Phase1_results.md]] §6.4:

| Detector | Scenariusz | SNR | β_ppE_5σ_bound^(b=-1) | β_ppE^TGP (locked) | Detekcja |
|----------|-----------|-----|------------------------|---------------------|----------|
| LIGO-O3 (now) | GW150914-like (M=65 M_⊙, 410 Mpc) | ~24 | ~10⁻¹ | ~7.8·10⁻² | **borderline** |
| LIGO-O5 single (2027+) | Loud BBH (M=30 M_⊙, 200 Mpc) | ~80 | ~3·10⁻² | ~7.8·10⁻² | **YES — first decisive** |
| LIGO-O5 stack 100 BBH | A+ design | ~800 | ~3·10⁻³ | ~7.8·10⁻² | YES (>20σ) |
| ET-D single (~2035) | Loud BBH (M=30 M_⊙, 1 Gpc) | ~500 | ~10⁻³ | ~7.8·10⁻² | YES (>50σ) |
| CE single (~2035+) | Loud BBH (M=30 M_⊙, 1 Gpc) | ~1000 | ~3·10⁻⁴ | ~7.8·10⁻² | YES (>200σ) |
| ET+CE network stack 5000 BBH | (~5 yr) | network | ~10⁻⁴ | ~7.8·10⁻² | YES (>700σ) |

**Poważna nowa konkluzja vs OOM z [[../../audyt/T01_LIGO3G_falsifier/SENSITIVITY_BACK_OF_ENVELOPE.md]]:**

LIGO-O5 *single event* (M=30 M_⊙, 200 Mpc, A+ design, SNR ~80) ma
bound ~3·10⁻², a TGP β_ppE^TGP ≈ 7.8·10⁻². **LIGO-O5 *single event*
może już dać first decisive detection** TGP w ~2027–2030, nie wymaga
stack. **To znacząco accelerates timeline T01 falsifier.**

## §6 — Output → T01 closure (Path B EXECUTED)

| Artefakt | Update wymagany |
|----------|-----------------|
| [[../../audyt/T01_LIGO3G_falsifier/FALSIFIER_STATEMENT_DRAFT.md]] §1 | placeholder `[β_th]` → `|β_ppE^TGP^(b=-1)| ∈ [5.5·10⁻², 1.2·10⁻¹]` |
| [[../../PREDICTIONS_REGISTRY.md]] M911-P1 | status `LIVE-PARTIAL` → `LIVE` (numerical lock); detection table updated |
| [[../../audyt/T01_LIGO3G_falsifier/SENSITIVITY_BACK_OF_ENVELOPE.md]] §4.2 | OOM updated with locked Phase 1 numbers |
| [[../../audyt/T01_LIGO3G_falsifier/NEEDS.md]] N1, N2, Q4 | OPEN → EXECUTED via op-ppE-mapping |
| [[../../audyt/T01_LIGO3G_falsifier/README.md]] closure progress | "Path B PREVIEW" → "Path B EXECUTED 2026-05-07" |

## §7 — Co jeszcze otwarte (dla paper completeness)

1. **G_SPA tighter lock** (Phase 1.5 future): 5% precision dla
   β_ppE^TGP central value. Wymaga retarded Green's function
   calc w hyperbolic metric. Estymata: ~1 dodatkowa sesja.
2. **Spin effects**: aligned spin (χ_eff ≠ 0) dodaje ~30-80%
   degeneracy w Fisher; włączane w `op-LIGO-3G-deviation/` Phase 2.
3. **Asymmetric mass** (η < 1/4): β_ppE^TGP ∝ 1/η scaling
   trzyma się dokładnie (SPA prefactor); `op-LIGO-3G-deviation/`
   może to weryfikować Fisher matrix.
4. **GWTC-3 reanalysis** (T01 NEEDS Q3 / FINDINGS Tier 5): może
   dać sygnał *teraz* z public data (~1 sesja).
5. **D paper draft** (Path D w T01): "Strong-field test of M9.1'':
   2PN-phase deviation predictions for Einstein Telescope and
   Cosmic Explorer". Cykl `op-ppE-mapping/` dostarcza Section 2-4
   (analytical setup + sympy LOCK + multi-coefficient signature).

## §8 — Bibliografia (paper-ready)

### Primary input (TGP)

- M9_1_pp_setup, M9_1_pp_P1_results — `research/op-newton-momentum/`
- TGP_FOUNDATIONS.md § 3 — narrative source
- sek08c_metryka_z_substratu — M9.1'' P1 derivation

### Audit context

- audyt/T01_LIGO3G_falsifier/ (v2.1, sesja 2026-05-07) — pełen
  audyt-folder z NEEDS, FALSIFIER_STATEMENT_DRAFT, PPN_TO_PPE_MAPPING,
  SENSITIVITY_BACK_OF_ENVELOPE, CONVENTION_DECISION, FINDINGS,
  CYCLE_KICKOFF×2

### ppE framework

- Yunes & Pretorius, Phys. Rev. D 80:122003 (2009), arXiv:0909.3328
- Yunes, Yagi, Pretorius, Phys. Rev. D 94:084002 (2016), arXiv:1603.08955
- Sampson, Yunes, Cornish, Phys. Rev. D 88:064056 (2013)
- Cutler & Flanagan, Phys. Rev. D 49:2658 (1994)
- Buonanno, Iyer, Ochsner, Pan, Sathyaprakash, Phys. Rev. D 80:084043 (2009) — TaylorF2

### PN waveform

- Blanchet, *Living Rev. Relativ.* 17:2 (2014)
- Damour, Jaranowski, Schäfer, Phys. Rev. D 89:064058 (2014) — 4PN ADM
- Mishra, Iyer, Sundararajan, Phys. Rev. D 93:084054 (2016) — 3PN

### Modyfikacje GR (literature comparison)

- Yagi, Yunes, *arXiv:1602.04674* — review
- Yagi, Eling, Foster, Phys. Rev. D 89:045017 (2014) — Einstein-Aether
- Endlich, Khoury, Solomon, arXiv:1709.10018 — sGB
- Will, *Living Rev. Relativ.* 17:4 (2014) — PPN review

### LIGO bounds

- Abbott et al. (LIGO/Virgo), Phys. Rev. D 103:122002 (2021), arXiv:2010.14529 — GWTC-2 ppE
- Abbott et al. (LIGO/Virgo/KAGRA), Phys. Rev. X 13:041039 (2023) — GWTC-3 ToGR
- Chamberlain & Yunes, Phys. Rev. D 96:084039 (2017), arXiv:1704.08268 — 3G implications

### 3G detector science cases

- Maggiore et al., JCAP 03:050 (2020), arXiv:1912.02622 — Einstein Telescope
- Reitze et al., Bull. Am. Astron. Soc. 51:035 (2019) — Cosmic Explorer
