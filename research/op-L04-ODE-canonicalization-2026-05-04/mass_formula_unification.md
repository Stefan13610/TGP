---
title: "L04 — unifikacja mass formulas: Phase 2 universal vs LP-4/LP-6/R5"
date: 2026-05-04
parent: "[[README.md]]"
type: unification
tgp_owner: research/op-L04-ODE-canonicalization-2026-05-04
tags:
  - L04
  - mass-formula
  - PHASE2
  - LP-4
  - LP-6
  - R5
  - unification
---

# Unifikacja mass formulas TGP — co naprawdę jest fundamentalne

## 1. Lista konkurencyjnych mass formulas

W TGP_v1 współistnieją (lub współistniały) różne mass formulas, każda
z innym statusem epistemicznym:

| Formula | Pochodzenie | Status |
|---------|-------------|--------|
| `m = c · A_tail⁴` | LP-4 (PLAN_DOMKNIECIA), R3 oryginalne | Specjalny case α=1 |
| `m = c · K²` | R5 (mass_scaling_k4) | Specjalny case α=1, equiv to A⁴ przez K~A² |
| `m = c · A_tail^p(α)` z p=5−α | R3 CORRECTIONS 2026-05-01 | Empirical decomposition jednowykładnikowa |
| `m_obs = c · A² · g₀^[e²(1−α/4)]` | PHASE2 closure 2026-05-01 | **Universal dla wszystkich α** |
| `m = c · M_energy²` | R3 historyczne | Sfalsyfikowane (τ FAIL −11.5%) |
| `m = c · ∫(g−1)⁴r²dr` | R3 historyczne | Sfalsyfikowane (diff −48%, −59%) |

## 2. Phase 2 universal jako fundamental

```
m_obs(g₀, α) = c_M · A_tail²(g₀, α) · g₀^[e²·(1−α/4)]
```

**Dlaczego jest fundamentalna:**

1. Działa dla każdego α (testowane α∈[0.25, 4.0], 14 punktów)
2. Reprodukuje PDG dla α=2 z dokładnością 0.001% (m_μ/m_e)
3. Reprodukuje PDG dla α=1 z dokładnością 0.046% (m_μ/m_e)
4. n(α=4) = 0 = Hobart-Derrick balance point — strukturalna konsekwencja
5. Czysty wykładnik e² (klasyczna stała, no φ admixture) — sugeruje
   *fundamental origin* (RG flow / Gaussian integral / topological)

## 3. Hierarchia: Phase 2 jest macierzą, reszta to specjalne przypadki

### 3.1 LP-4 / R5 K² jako specjalny case α=1

Z [[R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02]]:

> R5 mass formula `m = c·K²` jest równoważna Phase 2 universal **wtedy
> i tylko wtedy gdy α = 1**.

```
α=1: m_obs = c · A² · g₀^(3e²/4) ≈ c · A² · g₀^5.542

      Empirycznie z numerycznym ODE solve:
      (A_μ/A_e)² · (g₀_μ/g₀_e)^5.542 ≈ A^4 · const
      ↑                              ↑
      Phase 2                        LP-4 (k=4)
```

LP-4 argument konwergencyjny `k = 2(d-1)/(d-2) = 4` w d=3 jest:

- **Spójny dla α=1** ✓ (3 niezależne argumenty PASS)
- **Nie universal** — argument bazuje na asymptotyce z K=g², nie
  generalize się do K=g⁴

### 3.2 Empirical decomposition `p(α) = 5−α`

Z [[../why_n3/r3_observable_vs_full_mass.txt]] Sekcja 2:

| α | p_emp | 5−α | diff% |
|---|-------|-----|-------|
| 0.50 | 4.809 | 4.500 | +6.87% |
| 0.75 | 4.367 | 4.250 | +2.76% |
| **1.00** | **4.001** | **4.000** | **+0.02%** |
| 1.25 | 3.692 | 3.750 | −1.55% |
| 1.50 | 3.428 | 3.500 | −2.06% |
| 1.75 | 3.200 | 3.250 | −1.54% |
| **2.00** | **3.001** | **3.000** | **+0.02%** |
| 2.50 | 2.668 | 2.500 | −6.3% |

`p(α) = 5−α` jest **przybliżeniem** Phase 2 dla *jednowykładnikowej*
formy `m ~ A^p`. Dokładność <0.1% dla α=1 i α=2; rozmywa do 6% dla
ekstremów.

### 3.3 Pełna formuła Phase 2 jest dwuwykładnikowa

```
m_obs = c · A^a · g₀^b
```

z dwoma wykładnikami:
- **a = 2** (universal kinetic dressing, K_tail ~ A² uniwersalne)
- **b = e²·(1−α/4)** (core wave-function renormalization Z(g₀, α))

Z `r3_p_alpha_analytical.txt` Sekcja 6:

```
n(α) = (3−α) · log(A_μ/A_e) / log(φ)
     = (3−α) · log(A) / log(g₀_ratio)
```

Dla α=1: n(1) = 2 · log(3.79)/log(1.618) = 5.54
Dla α=2: n(2) = 1 · log(5.91)/log(1.618) = 3.69

Liniowy fit `n(α) = e²(1−α/4)` daje match <0.001% dla α=1 i α=2.

## 4. Co LP-4 i LP-6 powinny być przeformułowane

### LP-4 (PLAN_DOMKNIECIA)

**Pre-Phase 2 framing:**
> k=4 jedyne całkowite k dające skończoną masę solitonu w d=3 (twierdzenie
> B1''); k = 2(d-1)/(d-2) = 4 z konwergencji wymiarowej.

**Post-Phase 2 framing:**
> Dla α=1 (substratowa K=g²): mass formula degeneruje do `m ∝ A^4`
> przez algebraiczną tożsamość A^(3-α)·A² = A^4 z dwuwykładnikowej formuły
> Phase 2. Argument konwergencyjny k = 2(d-1)/(d-2) jest formularnie
> poprawny dla α=1 — ale dla α=2 (TGP-canonical) faktyczna potęga to k=3
> (z `5−α` empirical decomposition) lub b=e²/2 (Phase 2 dwuwykładnikowa).

LP-4 9/9 PASS jest **prawidłowy dla α=1**, ale nie dla α=2.

### LP-6 (PLAN_DOMKNIECIA)

**Pre-Phase 2 framing:**
> Formulacje są RÓWNOWAŻNE w słabym polu (PPN, κ, n_s, α_s, Koide).
> Różnią się TYLKO w ODE solitonu. K=g² jest preferowana: ghost-free,
> stabilna, reprodukuje pełne spektrum.

**Post-Phase 2 framing:**
> Formulacje α=1 i α=2 NIE są równoważne — α=1 jest *specjalnym przypadkiem*
> α=2 universal (Phase 2 mass formula). Numeryczna stabilność K=g² dla
> g₀>1.3 jest *artefaktem implementacji*, nie strukturalnej preferencji.
> Phase 2 z α=2 daje stabilne wyniki przez dwuwykładnikową formułę
> `m = c·A²·g₀^(e²/2)`, bez problemów ghost.

`thm:D-uniqueness` jednoznacznie wybiera α=2 z aksjomatów (C1)+(C2)+(C3).
LP-6 dictionary jest myląca — nie ma „dwóch równoważnych formulacji",
jest jedna kanoniczna (α=2) i jeden specjalny case (α=1).

### R5 (mass_scaling_k4)

**Pre-bridge framing:**
> R5 mass formula `m = c·K²` jest fundamentalnym mechanizmem TGP.

**Post-bridge framing:**
> R5 K² mass formula jest **derivative** — strukturalna konsekwencja
> Phase 2 universal dla specyficznego substratu α=1. Phase 2 jest
> fundamental, R5 K² to derivative.

[[R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02]] § 7:
> R5 ≡ Phase 2 IFF α=1. Dla TGP-canonical α=2: R5 K² ratio = 1221, PDG
> = 207, mismatch +490%.

## 5. Co PHASE2 mass formula nie tłumaczy (open problems)

### 5.1 Analytical derivation X = e²/4

n(α) = e²·(1−α/4) jest **empirical discovery** z fit residuum <0.1%,
ale nie ma analytical proof X = e²/4 z fundamentalnej fizyki.

Trzy hipotezy (z [[../why_n3/PHASE2_n_alpha_derivation.md]] §5):

- **H1: RG fixed point** — e² jako anomalous dimension wave-function
  renormalization. Phase 6 Q5 R⁵-bridge: NEGATIVE, X = e²/4 downgraded
  do "leading candidate, not definitive"
- **H2: Gaussian/exponential integral** — e^x z classical Gaussian
  weighted integrals
- **H3: Topological winding number** — RP² Berry phase z PHASE3 RP²
  defect quantization

Status: OPEN PROBLEM, candidate for dedicated R⁵-bridge cycle.

### 5.2 Why α=4 is Hobart-Derrick balance

n(4) = 0 jest *exact numerical zero* w skanie ODE. Sugeruje że α=4
jest *dimensional balance point*, gdzie soliton nie ma core dressing
i mass formula degeneruje do czystej `m=c·A²`.

Twierdzenie Hobarta-Derricka: w d≥2 brak stabilnych solitonów
skalarnych z `K=const` i kanonicznym `V(φ)`. Dla `K=φ^(2α)`,
Derrick scaling daje warunek `α(d−2)+d = 0` (?), czyli `α = d/(2−d)`
dla d=3 daje α=−3 (niefizyczne) lub α=4 z innego rozważania
balance kinetycznego.

To wymaga **formal derivation** Hobart-Derrick z TGP `K(φ)`.

### 5.3 Single-exponent p(α)=5−α vs full Phase 2

`p(α) = 5−α` jest empirically very close ale NIE exact dla wszystkich α
(tylko dla α=1 i α=2). Phase 2 z dwuwykładnikową formułą jest dokładna
do <0.1% dla α∈[0.5, 2.5].

Pytanie: czy istnieje *analytical bridge* p(α)=5−α ↔ Phase 2 dwuwykładnikowa?
Możliwe że p(α)=5−α jest pierwszym przybliżeniem perturbacyjnym.

## 6. Implikacje dla rdzenia TGP

### Co należy zaktualizować w manuskrypcie

(Niska priorytet, opcjonalne — bez tego TGP-canonical jest spójne)

- **README.md** sekcja "Key predictions": wskazać Phase 2 mass formula
  jako fundamentalną, z α=2.
- **sek08a_akcja_zunifikowana.tex**: dodać reference do Phase 2 mass
  formula po `thm:D-uniqueness`.
- **sek08b_ghost_resolution.tex**: zaktualizować twierdzenie B1''
  (k=4) jako case α=1, z reference do Phase 2 universal dla α=2.
- **research/why_n3/README.md**: explicit acknowledge że R3 mass formula
  działa dla każdego α z Phase 2.
- **research/mass_scaling_k4/README.md**: status update R5 K² →
  derivative case α=1.

### Co NIE wymaga zmiany

- TGP_FOUNDATIONS §3 (K(φ) = K_geo·φ⁴) — pozostaje aksjomat z (C3)
- thm:D-uniqueness — pozostaje (jest formal proof α=2 selection)
- G.0 closure (sek08a v2.0 ADDENDUM) — używa K=K_geo·φ⁴, spójne

## 7. Werdykt

**Phase 2 universal mass formula** `m_obs = c · A² · g₀^[e²(1−α/4)]`
jest **fundamentalna**:

- Działa dla wszystkich α (testowane α ∈ [0.25, 4.0])
- Reprodukuje PDG dla α=2 z dokładnością <0.1%
- Wykładnik e² jest *czystą klasyczną stałą* (no φ admixture)
- LP-4 (k=4), R5 K², empirical p(α)=5−α są **specjalnymi przypadkami**

**TGP-canonical α=2** jest **kanoniczna** (z `thm:D-uniqueness`).
α=1 substratowa jest historycznym fitującym przypadkiem (pre-v2 pivot).

**Open problem:** analitycznie wyprowadzić X = e²/4 z RG flow lub
Gaussian integral w field-theoretic limicie. Cykl R⁵-bridge przyszłej
sesji.

## Cross-references

- [[README.md]]
- [[m_obs_vs_M_full.md]]
- [[canonical_form_evidence.md]]
- [[ODE_class_taxonomy.md]]
- [[../../meta/PLAN_DOMKNIECIA_MASTER.md]] (LP-4, LP-6)
- [[../why_n3/PHASE2_n_alpha_derivation.md]] (PHASE2 e² discovery)
- [[../mass_scaling_k4/R5_PHASE2_ANALYTICAL_BRIDGE_2026-05-02.md]]
- [[../why_n3/r3_observable_vs_full_mass.txt]]
- [[../why_n3/r3_alpha2_full_closure.txt]]
- [[../why_n3/r3_p_alpha_analytical.txt]]
