---
title: "L04 — taxonomia klas ODE i selekcja α"
date: 2026-05-04
parent: "[[README.md]]"
type: classification
tgp_owner: research/op-L04-ODE-canonicalization-2026-05-04
tags:
  - L04
  - ODE-classification
  - alpha-selection
  - taxonomy
---

# Taxonomia klas operatorów kinetycznych i selekcja α

## 1. Klasa skalarnych operatorów II rzędu

Najogólniejszy lokalny, Lorentz-invariant operator kinetyczny:

```
L_kin = -½ K(φ)(∇φ)²
```

z gładką K: (0,∞) → ℝ_>0.

Po wariacji EL:

```
∇²φ + α(φ)·(∇φ)²/φ = 0,
α(φ) := K'(φ)·φ / (2K(φ))
```

## 2. Trzy warunki selekcji (`thm:D-uniqueness`)

[[../../core/sek08_formalizm/sek08_formalizm.tex]] lin. 956-1048.

### (C1) Stałość α

α(φ) = const ⟺ K'(φ)/(2K(φ)) = α/φ ⟺ K(φ) = C·φ^(2α)

**Implikacja:** każdy operator o stałym α odpowiada **K w postaci
potęgowej** φ^(2α).

### (C2) Warunek próżni N0-4: K(0) = 0

K(0) = C · 0^(2α) = 0 ⟺ 2α > 0 ⟺ α > 0

Wyklucza: α ≤ 0 (kinetic nieskończone w próżni Φ→0).

**Dopuszczalne po (C1)+(C2):** α ∈ (0, ∞).

### (C3) Geometric substrate coupling: K = K_geo · φ⁴

Z `prop:substrate-action` (dod. B `eq:Lkin-geo`):

```
K(φ) = K_geo · φ⁴
```

Porównując z (C1): C·φ^(2α) = K_geo·φ⁴ ⟹ **2α = 4** ⟹ **α = 2**.

## 3. Klasy α — co które dają

Dla każdej klasy K = C·φ^(2α) (po C1+C2):

| α | K(φ) | Kinetic Lagrangian | Status w TGP |
|---|------|---------------------|---------------|
| 0 | const | -½(∇φ)² | klasyczny Klein-Gordon (NIE TGP) |
| 1/2 | C·φ | -½φ(∇φ)² | exotic, brak fizycznej motywacji |
| 1 | C·φ² | -½φ²(∇φ)² | "substratowa" (R3 oryginalne) |
| 3/2 | C·φ³ | -½φ³(∇φ)² | inne (testowane w R3 α-scan) |
| **2** | **C·φ⁴** | **-½φ⁴(∇φ)²** | **TGP-canonical (Laplace-Beltrami konformalny)** |
| 5/2 | C·φ⁵ | -½φ⁵(∇φ)² | inne (R3 α-scan) |
| 3 | C·φ⁶ | -½φ⁶(∇φ)² | quartic kinetic (rare) |
| 4 | C·φ⁸ | -½φ⁸(∇φ)² | Hobart-Derrick balance (n(4)=0) |

**TGP-canonical α=2** odpowiada `K = K_geo·φ⁴`. Operator kinetyczny:

```
D_kin[φ] = ∇²φ + 2(∇φ)²/φ = (1/3φ²) ∇²(φ³)
```

To jest **laplasjan Laplace-Beltrami metryki konformalnej**
`h_ij = φ⁴·δ_ij` w R³ (`rem:LB-alpha2`).

## 4. Czemu właśnie φ⁴ z (C3)?

`prop:substrate-action` (dod. B) wynika z **post-pivotu v2 GL**
(2026-04-24):

- Pivot zmienił H_Γ z bilinearnego `−J Σ ŝᵢŝⱼ` na Ginzburg-Landau
  gradient bond `+J Σ A_ij ŝᵢ²ŝⱼ²(ŝⱼ²−ŝᵢ²)²` w Φ = ŝ²
- Coarse-graining tego H_Γ produkuje continuum K(φ) = K_geo · φ⁴

To jest *konkretna wewnętrzna struktura aksjomatu v2 GL*. Inne aksjomaty
(np. v1 bilinear) dawałyby inne K(φ).

## 5. „Substratowa α=1" jako historyczny artefakt

Pre-2026-04-24 (przed pivotem v2):

- v1 bilinear bond miał inne K(φ) — nie było jasne czy φ² czy φ⁴
- Numeryczne fits LP-4 (skrypty `lp4_*`, `r3_*`) używały K=g²
  (substratowa α=1) **bo była stabilna numerycznie** (kanoniczna K=g⁴
  niestabilna dla g₀>1.3)
- Fenomenologia (mass ratios, Koide) działała dla α=1 z formułą A⁴

Post-2026-04-24 (v2 GL):

- `prop:substrate-action` daje K(φ) = K_geo · φ⁴
- `thm:D-uniqueness` jednoznacznie wybiera α=2
- Phase 2 universal mass formula z α=2 daje <0.1% PDG
- α=1 staje się *historycznym fitującym przypadkiem*, nie kanoniczna
  formulacja

**Caveat (rem:alpha2-pivot-status-pl):**

> v1-style twierdzenie „α=2 wyprowadzone z bilinearnego wiązania
> -J Σ ŝᵢŝⱼ" **wycofane**. Nietrywialność `thm:D-uniqueness` polega na
> jednoznaczności α=2 wewnątrz klasy (C1-C3); rozluźnienie
> któregokolwiek z warunków (nielokalność, wyższe pochodne, inna
> Φ-kowariancja) otwierałoby inne wartości.

## 6. Klasyfikacja innych α (R3 PHASE6 alpha_em)

[[../why_n3/PHASE6_alpha_em_connection.md]]

**Specjalne wartości α z teorii:**

- α=2 (TGP-canonical) — kanoniczna z (C3)
- α=4 (Hobart-Derrick) — n(α)=0, m_obs=c·A² czyste; soliton w *balance*
  dimensional
- α=3 — n(3) = e²/4 ≈ 1.85; podejrzewana relacja z **α_em fine-structure**
  (4πα_em ≈ 0.0918, a (e²/4)·log(g₀^τ) ≈ ?) — open hypothesis
- α=1 — n(1) = 3e²/4 ≈ 5.54 (substratowa)

Linia α może być fundamentalnie kwantyzowana (α∈ℤ⁺/2 dla integral kinetic
exponents), choć nie jest to potwierdzone formalnie.

## 7. Bariera g₀_crit jako funkcja α

[[../why_n3/r3_alpha2_canonical_audit.txt]] Sekcja A.1:

| α | g₀_crit | Pochodzenie geometryczne |
|---|---------|---------------------------|
| 0.75 | 2.3703 | "geometryczne" |
| 1.00 | 2.2062 | substratowa (R3 default) |
| **2.00** | **1.8744** | **TGP-canonical, = 4/3·1.40554 (Lorentzian horizon M9.1'')** |

**Selection rule N=3** działa dla obu α=1 i α=2:
```
g₀^e (0.869) < g₀^μ (1.407) < g₀^τ (1.755) < g₀_crit < g₀^4
```

Dla α=2: g₀_crit = 1.874 *koincyduje numerycznie* z **horyzontem
Lorentzowskim metryki M9.1''**: `g_tt(ψ) = 0` dla ψ = 4/3, skąd via
identyfikację liniową ψ = 0.3814·g + 0.6186 otrzymujemy g₀ = 1.874
(why_n3 PHASE1, RESOLUTION 2026-05-01).

To jest **niezależne potwierdzenie geometrycznej kanoniczności α=2** —
bariera ODE pokrywa się z horyzontem metryki, co dla α=1 nie zachodzi.

## 8. Synteza: dlaczego TGP wybiera α=2

Trzy poziomy wyboru:

1. **Aksjomatyczny** (C3): post-pivot v2 GL `prop:substrate-action`
   wymaga K = K_geo·φ⁴, co z (C1)+(C2) daje α=2.
2. **Fenomenologiczny**: Phase 2 mass formula z α=2 jest spójna z PDG
   na poziomie 0.001% (m_μ/m_e), 0.085% (m_τ/m_e).
3. **Geometryczny**: bariera g₀_crit(α=2) = 1.874 koincyduje z
   horyzontem M9.1'' Lorentzowskim ψ=4/3 (4 cyfry znaczące).

α=1 substratowa jest historycznym fitującym przypadkiem (pre-v2 pivot),
zachowanym jako specjalny case Phase 2 universal mass formula. Nie
jest *równoważną* formulacją z α=2 — jest specjalizacją.

## Cross-references

- [[README.md]]
- [[m_obs_vs_M_full.md]]
- [[canonical_form_evidence.md]]
- [[mass_formula_unification.md]]
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`thm:D-uniqueness` lin. 956-1048
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`rem:alpha2-pivot-status-pl` lin. 1050
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`rem:LB-alpha2` lin. 1070+
- [[../why_n3/PHASE1_psi_g0_identification.md]] (ψ ↔ g₀ linear identification)
- [[../why_n3/r3_alpha2_canonical_audit.txt]] (bariera vs α scan)
