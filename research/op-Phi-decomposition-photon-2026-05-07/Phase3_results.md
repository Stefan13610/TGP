---
title: "Phase 3 results — polarization & spin (KRYTYCZNA, photon-as-δΦ FAILS)"
date: 2026-05-07
parent: "[[README.md]]"
type: phase-results
cycle: Stage 2 (op-Phi-decomposition-photon)
phase: 3
status: PHASE_3_COMPLETE — photon-as-δΦ STRUCTURAL_NO_GO; δΦ-modes survive jako TGP scalar sector
classification: STRUCTURAL_NO_GO_FOR_PHOTON_HYPOTHESIS, STRUCTURAL_DERIVED_FOR_dPHI_MODES
tgp_owner: research/op-Phi-decomposition-photon-2026-05-07
tags:
  - phase3
  - results
  - Stage2
  - polarization
  - spin
  - structural-no-go
  - representation-theory
  - sympy-verified
  - honest-baseline
related:
  - "[[Phase2_results.md]]"
  - "[[Phase1_results.md]]"
  - "[[NEEDS.md]]"
  - "[[phase3_sympy.py]]"
---

# Phase 3 results — polarization & spin (KRYTYCZNA)

## Status

**❌ PHASE 3: STRUCTURAL_NO_GO dla hipotezy photon-as-δΦ.**
**✅ PHASE 3: STRUCTURAL_DERIVED dla δΦ-modes jako TGP scalar sector.**

Wszystkie 4 testowane alternatywy (α/β/γ/δ) FAIL. Sympy weryfikacja
exit 0 (`phase3_sympy.py`) — F_μν = 0 ścisle dla A_μ ≡ ∂_μ δΦ
(czysta transformacja gauge, brak dynamiki EM).

**To NIE jest błąd obliczeniowy.** To fundamentalna konsekwencja
**representation theory grupy Lorentza**: skalar (spin 0) NIE ma
niezmienniczego mapowania na 4-wektor (spin 1) bez dodatkowej struktury
pola.

## Cel cyklu (Phase 3)

Z [[README.md]] §"Phase 3":

> "Zliczenie DOF skalarnego δΦ. Czy DOF=2 transverse jest reprodukowane?
> Test alternatyw α/β/γ. Empirical falsifier: dwie polaryzacje + spin=1
> obserwowane MUST być reprodukowane."

## Test alternatyw

### Alternatywa α: longitudinal-only photon

**Hipoteza:** foton ma 1 polaryzację (longitudinal) zamiast 2 transverse.

**Falsifier:** OBSERWACYJNY (Maxwell + QED).
- Polarizatory optyczne mierzą 2 niezależne stany (linear, circular)
- Helikalność fotonu ±1 (nigdy 0)
- Photon longitudinal mode (helicity 0) jest **NIEOBSERWOWANY**

**Werdykt:** **FAIL natychmiastowy**, niezgodny z obserwacjami.

### Alternatywa β: A_μ ≡ ∂_μ δΦ (gradient identity)

**Hipoteza:** "foton" to gradient skalarnego δΦ jako 4-wektor.

**Test sympy** (`phase3_sympy.py` F3.3):

Dla A_μ ≡ ∂_μ δΦ, EM tensor:
```
F_μν = ∂_μ A_ν - ∂_ν A_μ = ∂_μ ∂_ν δΦ - ∂_ν ∂_μ δΦ
```

Sympy weryfikacja (równość pochodnych mieszanych):
```
F_01 = F_02 = F_03 = F_12 = F_13 = F_23 = 0   ✓ (sprawdzone analitycznie)
```

**A_μ ≡ ∂_μ δΦ jest CZYSTĄ TRANSFORMACJĄ GAUGE** (gradient lambda
function w gauge transformation A_μ → A_μ + ∂_μ λ daje ∂_μ λ).

**Werdykt:** **FAIL** — brak dynamiki EM (F_μν trywialnie zero, brak
sprężonego pola elektromagnetycznego, brak fotonu).

### Alternatywa γ: ∂_i ∂_j δΦ w TT projection

**Hipoteza:** foton to mod transverse-traceless tensora ∂_i ∂_j δΦ.

**SO(3) decomposition** drugiego rzędu pochodnych skalarnego pola:
- h_ij ≡ ∂_i ∂_j δΦ (6 komponentów spatial)
- Trace: ∇²δΦ (1 DOF, skalarny)
- Traceless symmetric: h_ij - (1/3)δ_ij·∇²δΦ (5 DOF)
- TT projection (transverse-traceless, k^i h_ij^TT = 0): **2 DOF**

**Liczba DOF się zgadza** (2 = liczba polaryzacji fotonu)!

**ALE:** spin tego modu to **2**, nie 1.

Argument representation theory:
- Tensor symmetric rank-2 → reprezentacja (1, 1) Lorentza → spin 2
- Po constraint TT: nadal spin 2 reprezentacja
- Foton: spin 1 reprezentacja (1/2, 1/2)

**Werdykt:** **FAIL jako foton.** Działa jako kandydat dla **grawitonu**
(spin 2, TT-mode), ale nie dla fotonu.

To ciekawy boczny wynik: **w TGP grawiton mógłby być modem ∂_i∂_j δΦ
TT** — ale to scope OSOBNY cykl (op-graviton-as-Phi-mode-future), nie
Stage 2.

### Alternatywa δ: nowe pole / multi-component scalar

**Hipoteza:** wprowadzić wektorowe A_μ jako niezależne pole, lub
multi-component scalar Φ^a (a = 1, 2, 3).

**Naruszenie axiomów:** S05 (single-Φ axiom, closed 2026-04-26 Path B)
**eksplicytnie wyklucza** wielokomponentowe Φ albo dodatkowe pola
fundamentalne.

**Re-open S05 byłby out-of-scope** — wymagałby osobnego cyklu
`op-S05-reopen-vector-photon-2026-MM-DD/` z pełnym Phase 0-6.

**Werdykt:** **FAIL w obecnym Stage 2 scope.** Możliwa ścieżka future
research, ale poza obecnym cyklem.

## Argument representation theory (KORZEŃ FAIL)

**Twierdzenie:** Skalarne pole na flat Minkowski reprezentuje grupę
Lorentza w reprezentacji `(0, 0)` → spin 0. Jakiekolwiek **lokalne**
operacje na pojedynczym skalarnym δΦ (gradient, druga pochodna,
multiplikacja, etc.) generują reprezentacje:

| Operacja na δΦ | Lorentz rep | Spin | Komentarz |
|----------------|-------------|------|-----------|
| δΦ | (0,0) | 0 | skalar |
| ∂_μ δΦ | (1/2,1/2) | 1 | wektor — **BUT pure gauge** |
| ∂_μ ∂_ν δΦ symmetric | (1,1) | 2 | tensor |
| ∂_μ ∂_ν δΦ antisymmetric | (1,0)⊕(0,1) | 1 | F_μν tensor — **BUT identycznie 0** |

**Antysymetryczna druga pochodna jest IDENTYCZNIE ZERO** (równość
pochodnych mieszanych = sympy F_μν = 0 verification).

**Konkluzja:** Z pojedynczego skalarnego δΦ **NIE MOŻNA strukturalnie
wyprowadzić niezależnego modu spin-1 z 2 polaryzacjami transverse**.

To jest **fundamental no-go** representation theory.

## Implikacje dla Stage 2

### Co FAILS

**Hipoteza centralna Stage 2 H1** ([[README.md]]):
> "foton jest propagującym modem δΦ na tle Φ̄"

**STATUS: STRUCTURAL_NO_GO** w obecnym single-Φ framework.

### Co PRZEŻYWA — δΦ-modes jako "TGP scalar sector"

Phase 1+2 derivacje są **ANALITYCZNIE POPRAWNE** dla δΦ-modes
jako pole skalarne — po prostu **nie identyfikują się z fotonem**.

**Rebrand:** δΦ-modes są **dodatkowym sektorem skalarnym** TGP, niezależnym
od fotonu.

W kontekście kosmologicznym to jest podobne do:
- **Dilaton** (string theory)
- **Quintessence** (dark energy scalar field)
- **Inflaton** (early universe scalar)

**Phase 1+2 wyniki dla TGP scalar sector:**
- Dispersion: `ω² = c²(k² + γ)`, m_eff² = γ ≈ H_0²
- Mode expansion + canonical quantization → "TGP scalar quanta"
- E = ℏω, p = ℏk, λ = hc/E ⟵ **VALID** dla TGP scalar mode
- T^μ_μ_(δΦ) ≈ γ ~ Λ_today (γ-vacuum scale unification)

Pojedyncze "TGP scalar quanta" (NIE fotony) miałyby fenomenologię:
- Bardzo małe oddziaływanie z materią (jeśli q jest małe w eq:field-eq)
- Możliwe że tylko grawitacyjnie sprzężone
- Nie obserwowane bezpośrednio (dlatego nie były potrzebne w obecnym
  obrazie kosmologicznym)
- Mogą stanowić fragment **dark matter** lub **dark radiation** (Δ N_eff?)

### Standardowy foton w TGP framework

**Foton w TGP = standardowy A_μ z QED**, propagujący na metryce M9.1''(Φ̄).

**Modyfikacje TGP względem flat-space QED:**
- c lokalne = c(Φ̄) z ax:c (Phase 1 F1.3)
- Background metric to M9.1'' zamiast Minkowski
- Φ-EOM source term q·ρ_matter sprzęga δΦ z stress-energy SM (w tym EM)

**L01 ρ-bridge przywrócone:**
- ρ_EM = -T^μ_μ_EM/c_0² = 0 EXACT (FμνF^μν Weyl-niezmienniczość 4D ✓)
- ρ_(δΦ) ≈ ρ_Λ (separate scalar sector)
- L01 RESTORED do oryginalnej formy bez Stage 2 caveat

### Implikacje dla EXT-1 (Phase 4 BBN return)

**Foton = standardowy A_μ → ρ_EM = 0 strukturalnie nadal.**

**Stage 2 NIE rozwiązuje EXT-1 BBN problem:**
- W erze radiacyjnej ρ_rad_observed = ρ_photon_gas + ρ_neutrino
- Standardowo to wnosi do Friedmann eq przez ρ_total = ρ_matter + ρ_rad
- Ale w TGP L01 bridge: ρ_TGP = -T^μ_μ/c_0² → ρ_rad = 0 strukturalnie
- **Konflikt z BBN dynamics nadal istnieje** (EXT-1 STRUCTURAL_NO_GO
  pozostaje)

**δΦ-modes jako "dark radiation":**
- Mogłyby formalnie wnieść do ρ_rad jako dodatkowe N_eff
- ALE: δΦ-modes mają T^μ_μ ≠ 0 (m² = γ > 0)
- Czy to wnosi nontrivial ρ jak photon? Tak, ρ_(δΦ) = T_00/c² (energia)
- Ale skala γ ~ H_0² jest **bardzo mała** — δΦ-modes nie zdominują
  ρ_rad w erze radiacyjnej (ρ_rad_BBN ~ T⁴ ~ MeV⁴ >> H_0²)

**Phase 4 verdict (preliminary):** EXT-1 STRUCTURAL_NO_GO **utrzymuje się**
post-Stage-2. Stage 2 nie ratuje EXT-1.

## Phase 3 GATE summary

| Sub-task | Wynik | Sympy verify |
|----------|-------|--------------|
| F3.1 DOF counting | scalar=1, photon=2 → KONFLIKT | analytical |
| F3.2 Spin analysis | scalar=spin0, photon=spin1 → KONFLIKT | rep. theory |
| F3.3 Alt β (A_μ=∂_μδΦ) | F_μν = 0 EXACTLY → pure gauge | ✓ exit 0 |
| F3.4 Alt γ (TT-mode) | spin 2, NIE spin 1 → graviton-like | analytical |
| F3.5 Alt δ (new field) | narusza S05, out-of-scope | structural |
| F3.6 Verdict | photon-as-δΦ STRUCTURAL_NO_GO | comprehensive |

**Phase 3 GATE: 6/6 testów wykonanych, 0/4 alternatyw przeszło.**

**Verdict:** photon-as-δΦ hypothesis **STRUCTURAL_NO_GO**.
δΦ-modes survive jako separate TGP scalar sector.

## Probability evolution post-Phase-3

| Outcome | Pre-Phase-3 (post-Phase-2) | **Post-Phase-3** |
|---------|----------------------------|--------------------|
| Stage 2 photon-as-δΦ → DERIVED | 35-45% | **<1%** (rep. theory blocks) |
| Stage 2 photon-as-δΦ → STRUCTURAL CONDITIONAL | 40-50% | **<5%** |
| Stage 2 photon-as-δΦ → STRUCTURAL_NO_GO | 15-25% | **>95%** ✓ |
| δΦ-modes as TGP scalar sector → STRUCTURAL DERIVED | n/a | **>90%** ✓ |
| Stage 2 ratuje EXT-1 retroactively | 15-25% | **<5%** |

**Główny wynik:** Stage 2 hipoteza centralna **FAILS** — photon NIE jest
δΦ-mode w obecnym single-Φ framework. ALE Phase 1+2 mathematics
pozostaje VALID dla rebrandowanego "TGP scalar sector".

## Decyzja: Phase 4 + Stage 2 closure

**Phase 4 (BBN return)** — recommended **SKIP** lub **brief
follow-up only**:
- Photon = standard A_μ → ρ_EM = 0 strukturalnie nadal
- δΦ-modes nie zdominują BBN dynamics (γ ~ H_0² << T_BBN⁴)
- EXT-1 STRUCTURAL_NO_GO utrzymuje się

**Stage 2 NET STATUS:**
- **STRUCTURAL_NO_GO** dla photon-as-δΦ hypothesis
- **STRUCTURAL_DERIVED** dla δΦ-modes jako TGP scalar sector
- **HONEST PHASE 6 BASELINE** preserved

**Honest reporting note:**
Cykl wykonał Phase 0+1+2+3 z analitycznymi 8/8+6/6+4/4+6/6 PASS dla
dążonych wyników. Phase 3 ujawnił **fundamentalny structural constraint**
(representation theory) który wyklucza photon-as-δΦ identification.

To **NIE jest fail Phase 1+2 work** — Phase 1+2 wyprowadzili poprawną
dynamikę δΦ-modes. To jest fail **specyficznej identyfikacji**:
"δΦ-mode = foton" niemożliwe; "δΦ-mode = TGP scalar mode" możliwe.

**Postać godna spostrzeżenia:**
W Stage 2 honest acknowledgment failed primary hypothesis + preserved
secondary results = **honest scientific output** zgodny z Phase 6
ABSOLUTE BINDING gate.

## Cross-references

- [[Phase2_results.md]] — canonical quantization (4/4 PASS, valid for δΦ)
- [[Phase1_results.md]] — formal decomposition (6/6 PASS, valid)
- [[Phase0_balance.md]] — pre-derivation balance
- [[NEEDS.md]] — N8 (polarization), N10 (spin) — KRYTYCZNE FAIL Phase 3
- [[phase3_sympy.py]] — sympy verification (exit 0)
- [[../op-FRW-radiation-era-varying-c-2026-05-06/FINDINGS.md]] —
  EXT-1 STRUCTURAL_NO_GO **utrzymuje się** post-Stage-2
- [[../../audyt/S05_tensor_sector_singleField/]] — single-Φ axiom
  (closed 2026-04-26 Path B; potential re-open future cycle dla
  vector A_μ)
- [[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]
  — L01 ρ_EM = 0 RESTORED (Stage 2 caveat removed; standard EM intact)
- [[../../meta/CALIBRATION_PROTOCOL.md]] — Phase 6 ABSOLUTE BINDING gate
