---
title: "L01 — mapping ρ na sektory SM (Dirac, scalar, EM, Yang-Mills)"
date: 2026-05-04
parent: "[[README.md]]"
type: sector-mapping
tgp_owner: research/op-L01-rho-stress-energy-bridge-2026-05-04
tags:
  - L01
  - SM-sectors
  - Dirac
  - EM
  - Yang-Mills
  - stress-energy
---

# Mapping `ρ = -T^μ_μ/c_0²` na sektory Standard Model

## 1. Sektor Dirac (fermiony) — m_e, m_μ, m_τ, kwarki

**Lagrangian:**

```
L_Dirac = ψ̄ · (i γ^μ D_μ - m) · ψ
```

z `D_μ` kowariantna pochodna na curved background `g_eff`.

**Stress-energy tensor:**

```
T^μ_ν[ψ_D] = (i/2) · [ψ̄ γ^μ D_ν ψ - (D_ν ψ̄) γ^μ ψ]
            + δ^μ_ν · L_Dirac
```

Klasycznie (on-shell, EOM `(i γ^μ D_μ - m) ψ = 0`):

```
L_Dirac |_on-shell = 0   (EOM trywialnie zeruje L)
T^μ_ν |_on-shell = (i/2) · [ψ̄ γ^μ D_ν ψ - h.c.]
```

**Trace** (po kontrakcji):

```
T = T^μ_μ |_on-shell = m · ψ̄ψ
```

(standardowy wynik: trace stress-energy fermion = m × scalar density)

**Definicja ρ:**

```
ρ_Dirac = -T^μ_μ / c_0² = -m · ψ̄ψ / c_0²
```

Dla fermion w stanie spoczynku: `ψ̄ψ = ρ_number · 2m/c_0` (normalizacja
fermion bilinear). Po podstawieniu:

```
ρ_Dirac = -m · (ρ_number · 2m/c_0) / c_0²
        = -2 · ρ_number · m² / c_0³
```

W jednostkach naturalnych `c_0=1`: `ρ_Dirac = -2 ρ_number · m²`.
Z konwencji minusa w `T = -ρ_rest c²`: ρ_Dirac > 0 (mass density positive).

**Konkretnie** dla pojedynczego elektronu:

```
ρ_Dirac(e) = m_e · |ψ_e|²    [M·L⁻³]
```

To jest **standardowa gęstość masy**, identyczna z relativistycznym
ρ_rest dla nierelatywistycznych fermionów. ✓ konsystentne z SPARC fits
(N-body uses `ρ = ρ_baryon`).

**Dla wszystkich generations:**

- `ρ_e ∝ m_e · |ψ_e|²`
- `ρ_μ ∝ m_μ · |ψ_μ|²` (massive lepton, większy ρ przy tym samym density)
- `ρ_τ ∝ m_τ · |ψ_τ|²`
- analogicznie kwarki (u, d, s, c, b, t)

## 2. Sektor skalarny (Higgs, inne skalary)

**Lagrangian** (general scalar with m i potential V):

```
L_scalar = ½ · g_eff^μν · ∂_μφ · ∂_νφ - V(φ)
        ≈ ½ · g_eff^μν · ∂_μφ · ∂_νφ - ½ m²φ² - λφ⁴/4 - ...
```

**Stress-energy tensor:**

```
T^μ_ν[φ_scalar] = ∂^μφ · ∂_νφ - δ^μ_ν · L_scalar
```

**Trace:**

```
T = T^μ_μ = (∂φ)² - 4·L_scalar
          = (∂φ)² - 2·(∂φ)² + 4·V(φ)
          = -(∂φ)² + 4·V(φ)
          = -(∂φ)² + 2m²φ² + λφ⁴ + ...
```

Dla Higgs scalar w stanie próżni `<φ> = v + h(x)`, expansion wokół v:

```
T |_vacuum = 4·V(v) = -4·m²v²/2 + 4·λv⁴/4 = -2m²v² + λv⁴
```

co dla Higgs SSB: `2m²v² = λv⁴` ⇒ `T = 0` w klasycznej próżni
(Higgs vacuum jest *traceless* klasycznie).

**Po kwantyzacji:** trace anomaly daje `T^μ_μ ~ m²·h² + ...` dla
fluktuacji Higgsa h(x). **Definicja ρ:**

```
ρ_Higgs = -T^μ_μ / c_0²
        ≈ ((∂h)² - m_H² · h²) / c_0²   (post-vacuum subtraction)
```

To jest **density energii fluktuacji Higgsa**, klasycznie cienka i
quantum-modyfikowana. **Open problem (NEEDS N1):** quantum trace anomaly
dla Higgs sektor.

## 3. Sektor EM (massless gauge field) — fotony

**Lagrangian:**

```
L_EM = -¼ · F_μν · F^μν   (with g_eff^μν metric)
```

**Stress-energy tensor:**

```
T^μ_ν[A] = -F^μλ · F_νλ + ¼ · δ^μ_ν · F_αβ·F^αβ
```

**Trace:**

```
T = T^μ_μ = -F^μλ · F_μλ + (1·F_αβ·F^αβ)
          = -F²+ F²
          = 0
```

**Klasycznie EM jest traceless tensor** — to jest standardowy wynik
(conformal invariance massless gauge field in 4D).

**Definicja ρ:**

```
ρ_EM = -T^μ_μ / c_0² = 0    [klasycznie]
```

### Implikacje dla TGP

1. **Fotony NIE generują pola Φ przez `L_mat`**: `L_mat ∝ ρ_EM = 0`.
2. **Fotony propagują się geodezyjnie po `g_eff`** (`prop:coupling-consequences`
   pkt 2): `c(Φ) = c_0 · e^(-U)` w weak-field approximation.
3. **GW170817 c_GW = c_EM exactly** (`prop:coupling-consequences` pkt 3):
   bo zarówno GW jak i EM używają tej samej metryki `g_eff`.
4. **Brak piątej siły między fotonami i innymi cząstkami**: foton nie
   ma mass-coupling do Φ przez ρ.

### Quantum trace anomaly EM

W QFT na curved background, `T^μ_μ_EM ≠ 0` w QFT 1-loop. Adler-Bardeen
formula:

```
T^μ_μ_EM,quantum = (β(α)/(2α)) · F_μν·F^μν + ... · R²
```

z `β(α)` running of fine structure. To jest *small correction* O(α²)
w typowych warunkach lab; w plasma dense regimes może być significant.

**Implikacje:** w wysokich gęstościach EM (gwiazdy neutronowe, wczesny
Wszechświat, lab Schwinger regime) `ρ_EM ~ α² · F²` może dawać małe
ale niezerowe sprzęganie EM z Φ. To jest **open problem (NEEDS N1)** —
relevantne dla τ.3 / ψ.1 / ω.1 cycles.

## 4. Sektor Yang-Mills (gluon, W, Z) — non-Abelian gauge fields

**Lagrangian:**

```
L_YM = -¼ · Tr(G_μν · G^μν)
```

z G_μν = ∂A - ∂A + g[A,A] (non-Abelian).

**Stress-energy tensor:** analogicznie do EM, ale z trace nad grupową
matrices.

**Trace klasyczny:**

```
T = T^μ_μ = 0   (klasycznie traceless, podobnie jak EM dla massless YM)
```

**Quantum trace anomaly (kluczowe dla QCD):**

```
T^μ_μ_YM = (β(g)/(2g)) · Tr(G_μν·G^μν) ≠ 0
```

z β(g) = -b₀g³ + ... (asymptotyczna swoboda dla SU(3): b₀ > 0).

W QCD: trace anomaly daje `T^μ_μ ~ Λ_QCD⁴` — to jest **non-perturbative
gluon condensate**.

**Definicja ρ:**

```
ρ_YM = -T^μ_μ_YM / c_0² ≈ Λ_QCD⁴ / c_0²    [quantum, gluon condensate]
```

To jest **gęstość energii QCD vacuum** O(Λ_QCD⁴) ~ (300 MeV)⁴ ≈
10⁻⁴ GeV⁴ ≈ standard QCD condensate.

### Implikacje dla TGP

1. **QCD vacuum sprzęga z TGP przez Λ_QCD⁴** — bardzo małe ale niezerowe
2. **W kontekście kosmologicznym** (wczesny Wszechświat): faza-przejście
   QCD wpływa na TGP `Φ` przez `ρ_YM`
3. **Open problem**: explicit *non-perturbative* derivation `ρ_YM` w
   running coupling regime → relevant dla cosmological closure

## 5. Płyn doskonały (perfect fluid) — kosmologia

**Stress-energy tensor:**

```
T^μ_ν[fluid] = (ρ_e + p) · u^μ · u_ν + p · δ^μ_ν
```

z `ρ_e` energy density, p ciśnienie, u^μ four-velocity.

**Trace:**

```
T = T^μ_μ = -ρ_e + 3p
```

(z signature konwencji).

**Definicja ρ_TGP fluidem:**

```
ρ_fluid_TGP = -T^μ_μ / c_0² = (ρ_e - 3p) / c_0²
```

Specjalne przypadki:
- **Dust** (p=0): `ρ_TGP = ρ_e / c_0² = ρ_rest` (standardowa gęstość masy)
- **Radiation** (p = ρ_e/3): `ρ_TGP = 0` (nie sprzęga z Φ przez ρ!)
- **Dark Energy** (p = -ρ_e): `ρ_TGP = 4·ρ_e / c_0²` (silne sprzęganie)

**Konsekwencja kosmologiczna:** w epoce radiation-dominated wczesnego
Wszechświata, *radiacja nie generuje pola Φ przez ρ* — Φ jest dynamicznie
fixed przez non-radiative components (matter + DE). To jest **predykcja
TGP**: w radiation era, Φ-evolution nie jest source-dominated.

## 6. Tabela zbiorcza

| Sektor | T^μ_μ klasycznie | ρ_TGP |
|--------|-------------------|--------|
| Dirac fermion (m≠0) | m·ψ̄ψ | m·|ψ|²/c_0² ≥ 0 |
| Massive scalar | -(∂φ)² + 4V | (∂φ)²-4V (signs depend on phase) |
| Higgs vacuum | 0 (klasycznie) | 0 (post-SSB cancellation) |
| Photon (massless EM) | **0** (conformal) | **0** (no coupling) |
| Yang-Mills classical | 0 | 0 |
| Yang-Mills quantum (QCD) | β·G²/(2g) ~ Λ_QCD⁴ | ~Λ_QCD⁴/c_0² (gluon condensate) |
| Dust | -ρ_e | ρ_e/c_0² |
| Radiation | 0 | **0** |
| Dark Energy (p=-ρ_e) | -4ρ_e | 4ρ_e/c_0² |

## 7. Kluczowe wnioski

### Co fotony robią a czego nie

- ✓ Propagują się geodezyjnie po `g_eff` (c(Φ) = c_0/√ψ w weak-field)
- ✓ GW170817: c_GW = c_EM exactly (bo wspólna metryka)
- ✗ NIE generują pola Φ przez `L_mat` (`ρ_EM = 0` klasycznie)
- ⚠ Quantum trace anomaly: małe ale niezerowe `ρ_EM ~ α²·F²` w wysokich
  gęstościach EM (cyclic ψ.1, τ.3, ω.1 territory)

### Co TGP przewiduje dla materii

- **Dirac fermiony**: standardowa `ρ = m·|ψ|²` — pełne sprzęganie z Φ
- **Higgs**: classical vacuum traceless, but fluctuations couple; tied
  to G.0 closure 2026-05-02 vacuum stability fix
- **QCD vacuum**: gluon condensate z trace anomaly daje ρ ~ Λ_QCD⁴
- **Radiation era**: Φ-evolution NIE jest source-dominated (predykcja
  kosmologiczna)
- **Dark Energy**: silne sprzęganie ρ ~ 4ρ_e — to jest mechanizm Λ
  z `closure_2026-04-26 T-Λ` (Φ_eq scale)

### Implikacje dla testów

- Eöt-Wash/MICROSCOPE: testują `ρ_Dirac` differences między atomami
  (Pt vs Ti) — TGP daje η_TGP = 1.32×10⁻²⁶ << bound 1.1×10⁻¹⁵ (B9 6/6 PASS)
- Tested EM dispersion w ψ.1.v2: post-A6 invalidation v1 → tylko v2
  (varying-α nie varying-c) — spójne z `ρ_EM = 0`
- GW170817: c_GW = c_EM ≤ 3×10⁻¹⁵ — automatycznie z ax:metric-coupling

## Cross-references

- [[README.md]]
- [[formal_definition.md]]
- [[photon_treatment.md]]
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`ax:metric-coupling` §`prop:coupling-consequences`
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]] (test ρ_Dirac kompozycji)
- [[../op-psi1-substrate-light-acceleration]] (ψ.1 v2 — varying-α)
- [[../closure_2026-04-26/Lambda_from_Phi0]] (T-Λ vacuum energy)
