---
title: "NEEDS — L03 spectral stability (otwarte problemy)"
date: 2026-05-06
parent: "[[README.md]]"
type: needs
tgp_owner: research/op-L03-spectral-stability-2026-05-06
tags:
  - needs
  - L03
  - spectral-rigorous
  - functional-analysis
  - quantum-corrections
---

# NEEDS — L03 (otwarte problemy)

## Status

Po tym cyklu, audit L03 (`audyt/L03_K_phi_stability/`) jest substantialnie
**ZAMKNIĘTY** kombinacją pre-existing core (~70%) + 3 nowymi elementami
tego cyklu (~30%): mode_counting_Z2.md, tachyon_check_with_source.md,
spectral_synthesis.md.

**Pozostają otwarte problemy** w domenie *rigorous functional analysis*
(P3/P4), nie blokujące physical predictions.

## Otwarte luki

### N1: Rigorous lower bound spektrum gap

**Status:** OPEN — P3 (nie blokujący predykcji).

**Problem:** Asymptotyczny mass gap m_sp² = γ/K_geo > 0 jest dowiedziony
z pre-existing `prop:vacuum-stability` (linearization wokół jednorodnego ψ=1).
Brakuje **rigorous proof**, że dla *każdego* niejednorodnego tła Φ_eq[ρ]
(np. nontrivial Yukawa profile, soliton continuum):

```
inf σ(L̂_{Φ_eq}) = m_sp² = γ/K_geo
```

**Co zostało dowiedzione (po tym cyklu):** σ(L̂) ⊂ [0, ∞), continuum
zaczyna się przy m_sp² asymptotycznie.

**Co potrzebne:**

1. **Formal Rayleigh quotient analysis**: pokazać, że
   `inf_{δφ ≠ 0} ⟨δφ, L̂ δφ⟩/⟨δφ, δφ⟩ = m_sp²` (nie mniejszy).
2. **Min-max characterization**: λ_0(L̂_{Φ_eq}) jako infimum z funkcjonału
   energii.
3. **Bound state classification dla solitonów**: explicit count i value
   bound states λ_n^bound dla solitonów lekkich (e, μ, τ).

**Kandydat dostawcy:** dedicated cycle `op-L03-rigorous-spectral-bounds`
(estymata: ~3-4 tygodnie formal mathematical physics)

**Typ:** rigorous functional analysis

**Wpływ:** academic completeness; nie blokuje predykcji ani existential.

### N2: Functional analysis domain conditions (self-adjointness)

**Status:** OPEN — P4 (mathematical rigor).

**Problem:** Operator L̂ = -∂[K(φ_eq)·∂] + V''_eff(φ_eq) jest *formal*
operator różniczkowy. Dla rigorous spectral theory wymagana jest
explicit specyfikacja:

1. **Domain D(L̂) ⊂ Hilbert space H**: typowo `H = L²(ℝ³, w(x)dx)` z
   weight `w(x) = √(-g_eff(x))`
2. **Self-adjointness**: czy `L̂*` ma tę samą domain co L̂?
   - Sturm-Liouville singularny problem (singularność w r=0 jeśli soliton
     ma `g(0) > 1`) wymaga limit-point/limit-circle classification
     (Weyl).
3. **Sobolev space argumentation**: `D(L̂) ⊂ H¹(ℝ³)` lub W²,²?
4. **Boundary conditions w r=0 i r=∞**: regularność + decay.

**Co potrzebne:** referencja do podręcznika fizyki matematycznej
(Reed-Simon V2 §X.1 Sturm-Liouville, V4 spectral theory of partial
differential operators).

**Kandydat dostawcy:** krótki dodatek (1-2 dni) w `dodatekN_kwantyzacja.tex`
lub osobny dodatek `dodatekZ_spectral_theory.tex`.

**Typ:** mathematical rigor

**Wpływ:** academic, nie wpływa na fizyczne predykcje.

### N3: Quantum corrections to spectrum (1-loop)

**Status:** OPEN — P3 (relevant dla precyzji predykcji).

**Problem:** Klasyczny m_sp² = γ/K_geo. W kwantowej teorii pola na
curved background, 1-loop effective action ma korekcje:

```
m_sp²_eff = m_sp²_classical + (1/16π²)·[γ·R + ξ·R + α₁·m_sp⁴ + ...]
```

gdzie R jest curvature scalar, ξ jest non-minimal coupling.

**Konsekwencje:**

- Dla kosmologicznych perturbacji (FRW, R = 6(H'+2H²)): może dawać
  effective screening lub anti-screening m_sp²
- Dla strong-field regions (BH, neutron stars, magnetary): R ~ M/r³
  lokalnie duże, korekcje niemałe
- Dla GW propagation: dyspersja może być modyfikowana 1-loop

**Co potrzebne:**

1. Background field method dla L_TGP na curved background g_eff
2. Heat kernel expansion DeWitt-Schwinger dla operator L̂
3. Renormalization conditions (counter-terms γ_{0}, ξ_{0})
4. Phenomenological estimate korekcji w typowych physical regimes

**Source:** [[../../audyt/L03_K_phi_stability/README.md]] §"Co brakuje"
pkt. 1 (extension)

**Kandydat dostawcy:** dedicated cycle `op-quantum-corrections-Phi-spectrum`
(estymata: ~4-6 tygodni; requires QFT on curved spacetime expertise)

**Typ:** derivation (1-loop effective action)

**Wpływ:** może uściślić mass spectrum lepton predictions w 4-th znaczącej
cyfrze; obecnie PDG <0.01% już osiągnięte klasycznie.

### N4: Spectrum stability w ekstremalnych regimes

**Status:** OPEN — P3.

**Problem:** Spectral analysis tego cyklu zakłada physical sources ρ ≥ 0
i tła Φ_eq z `g_min ≥ 0.91 > g_ghost`. W *ekstremalnych* configurations:

- **Black hole horizon (M9.1'' Lorentzian)**: ψ_horizon = 4/3, ψ → ∞
  gdy r → r_g. Spectrum perturbacji w tej geometrii?
- **Big Bang (early universe)**: ψ → 0 limit (ale fizycznie unreachable
  wg N0-1)
- **Magnetar polar (B ~ 10¹¹ T, T-α regime)**: lokalnie strong gradient
  ∇Φ duże, czy wpadamy w ghost regime locally?

**Co potrzebne:**

1. Spectrum perturbacji na Schwarzschild-like background (TGP M9.1'')
2. Spectrum w pre-Big-Bang epoch (Φ → small)
3. Magnetar local field analysis

**Kandydat dostawcy:** część cyklu N3 (1-loop) lub osobne
`op-L03-extreme-regime-stability` (estymata: ~3 tygodnie).

**Typ:** numerical + asymptotic analysis

**Wpływ:** relevant dla τ.3 magnetar polar shift TT10 (nie krytyczny).

### N5: Spectrum gap dla soliton sektor (explicit values)

**Status:** OPEN — P3.

**Problem:** Dla solitonowych konfiguracji (lepton e, μ, τ), bound state
spektrum operatora L̂ wokół g(r) zawiera:

- 1 zero translation mode (gauge, niefizyczne)
- N bound states z λ_n ∈ (0, m_sp²) (analog "Higgs around kink")
- Continuum z λ ≥ m_sp²

**Co potrzebne:**

1. Numerical eigenvalue solver dla L̂_{soliton}
2. Tabela bound state eigenvalues per lepton (e, μ, τ)
3. Sprawdzenie że żaden bound state ma λ < 0 (no instability)

**Kandydat dostawcy:** krótki numerical cycle `op-L03-soliton-spectrum-numeric`
(estymata: 1 tydzień)

**Typ:** numerical

**Wpływ:** dodatkowa weryfikacja stabilności solitonów; obecnie
już domknięte przez `prop:nonlinear-stability` energetic argument.

## Pytania otwarte

- **Q1:** Czy ξ (non-minimal coupling) w L_field powinien być traktowany
  jako wolny parametr, czy fixed przez aksjomatykę N0?

- **Q2:** Dla M9.1'' Lorentzian black hole solution (TGP), czy istnieje
  Hawking-like radiation pole Φ?

- **Q3:** Spektrum gap m_sp² zależy od γ i K_geo. Czy γ/K_geo jest
  uniwersalną stałą fizyczną (jak Λ_QCD), czy może drift?

## Blokery

| ID | Bloker | Czeka na | Source |
|----|--------|----------|--------|
| B1 | N1 (rigorous spectral bounds) | dedicated math-physics cycle | functional analysis |
| B2 | N2 (functional analysis) | reference Reed-Simon V2/V4 | mathematical rigor |
| B3 | N3 (quantum corrections) | QFT on curved spacetime cycle | 1-loop effective action |
| B4 | N4 (extreme regimes) | curved background analysis | strong-field/early universe |
| B5 | N5 (soliton numerics) | numerical solver | implementation |

## Closed needs (po L03 cyklu 2026-05-06)

| ID | Co | Domknięte przez | Data |
|----|----|-----------------|------|
| Mode counting Z₂-broken | 1 massive scalar, 0 NGB, 0 ghost | mode_counting_Z2.md | 2026-05-06 |
| Tachyon check 4 profiles | brak tachyon mode na vacuum, FRW, Yukawa, soliton | tachyon_check_with_source.md | 2026-05-06 |
| S-L synthesis | thm:spectral-synthesis-L03 zsyntetyzowane | spectral_synthesis.md | 2026-05-06 |
| Ghost wall vs g_min | g_min ≥ 0.91 > g_ghost dla obu formulacji α=1, α=2 | sek08b cor:ghost-artifact + F7 | pre-existing + 2026-05-06 |

## Cross-references

- [[README.md]] — cykl L03 indeks
- [[FINDINGS.md]] — eksportowalne wyniki
- [[../../audyt/L03_K_phi_stability/]] — audit-source
- [[../../audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-05-06.md]] — będzie utworzony
- [[../op-quantum-corrections-Phi-spectrum]] (przyszły cykl, nie istniejący — N3)
- [[../op-L03-rigorous-spectral-bounds]] (przyszły cykl — N1)
- [[../op-L03-extreme-regime-stability]] (przyszły cykl — N4)
