---
title: "Mode counting w fazie Z₂-broken — pojedynczy massive scalar, no Goldstone"
date: 2026-05-06
parent: "[[README.md]]"
type: derivation
tgp_owner: research/op-L03-spectral-stability-2026-05-06
tags:
  - L03
  - mode-counting
  - Z2-symmetry
  - Goldstone-theorem
  - discrete-symmetry-breaking
---

# Mode counting w fazie Z₂-broken

## Pytanie z audytu

[[../../audyt/L03_K_phi_stability/README.md]] §"Co brakuje" pkt. 2:

> **Mode counting w obecności złamania Z₂** — Z₂ symmetry breaking daje
> 1 Goldstone (dla ciągłego), ale Z₂ jest dyskretne ⇒ nie ma Goldstone
> bosona. Powinien być tylko 1 massive mode. Czy jest?

Krótka odpowiedź: **TAK, dokładnie 1 massive scalar mode, zero Goldstone.**

## Tło fizyczne

### Symetria Z₂ w TGP

[[../../TGP_FOUNDATIONS.md]] §1 (warstwa 0): pojedyncze fundamentalne
pole skalarne Φ z dyskretną symetrią Z₂:

```
Z₂: Φ(x) → -Φ(x)
```

Lub równoważnie w zmiennej `ψ = Φ/Φ_0`:

```
Z₂: ψ(x) → -ψ(x)
```

W zmiennej substratowej `g = Φ/Φ_0` (sek08b §ssec:substrate-resolution):
```
Z₂: g(x) → -g(x)  (per amplituda `φ_i = √g_i` z hat structure)
```

### Spontaniczne złamanie Z₂

W stanie próżni `ψ_vac = ±1`:

- Aksjomatyka N0-5 + warunek próżni `β = γ` daje `V'(1) = 0`
- Wybierając jedną gałąź (np. `ψ_vac = +1`), Z₂ jest *spontanicznie* złamane
- Konfiguracja "kink" interpoluje między `ψ = +1` a `ψ = -1` (sek01 ontologia,
  warstwa 3c domeny defektowe — ale to NIE są dynamiczne DOF perturbacyjne)

## Goldstone theorem dla ciągłych vs dyskretnych grup

### Klasyczny Goldstone theorem (Goldstone 1961, Nambu 1960, Goldstone-Salam-Weinberg 1962)

> **Twierdzenie:** Spontaniczne złamanie *ciągłej* symetrii globalnej
> generuje masową lukę: dla każdego złamanego generatora T^a istnieje
> bezmasowy boson Goldstone'a φ^a z relacjami dyspersji `ω² = c²k²`.

Liczba Goldstone bosonów = `dim(G) - dim(H)`, gdzie G jest grupą całkowitą
a H podgrupą stabilizującą próżnię.

### Z₂ jest dyskretne

**Z₂ = {1, -1}** ma 0 generatorów ciągłych (tylko inwolucję). W skutek tego:

```
dim(G) = 0  ⇒  liczba Goldstone bosonów = 0
```

Nie ma bezmasowego bosona z spontanicznego złamania Z₂. Domena/kink jest
*globalną* konfiguracją (defekt topologiczny w warstwie 3c), nie *lokalnym*
bezmasowym DOF.

Patrz: Peskin-Schroeder §11.1 (linear sigma model bez NGB dla discrete);
Weinberg V2 §19.2; Gell-Mann-Lévy 1960 (sigma model rozróżnienie linear
vs nonlinear).

## Liczenie DOF w TGP

### Pole klasyczne

`Φ(x)` jest pojedynczym rzeczywistym polem skalarnym (1 DOF lokalnie):

```
DOF_total = 1
```

### Po SSB Z₂ → identyfikacja

Wybranie `ψ_vac = +1` nie redukuje liczby DOF — Z₂ jest dyskretne, więc
nie ma orbit ciągłej symetrii. Cała 1 DOF pozostaje dynamiczna.

```
DOF_NGB = 0  (Z₂ dyskretne)
DOF_massive = 1  (pole δψ = ψ - 1)
DOF_total = 0 + 1 = 1  ✓
```

### Sprawdzenie z `prop:vacuum-stability`

Z [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`prop:vacuum-stability`
linearyzacja:

```
∇²(δΦ) + (2β - 3γ)·δΦ = -q·Φ_0·δρ
```

Dla `β = γ`: `m_sp² = -(2β - 3γ) = 3γ - 2β = γ > 0`.

Propagator: `G(k) = 1/(k² + m_sp²)`. **Jeden** propagator skalarny z masą
`m_sp = √γ`. Brak innych modes.

W zmiennej bezwymiarowej:
```
m_sp² = γ / K_geo  (z prop:vacuum-stability-G0 sek08a:889)
```

(Czynnik K_geo z `K(φ) = K_geo·φ⁴`, K(1) = K_geo.)

## Symetria dyskretna — analogia z innymi modelami

| Model | Symetria | Złamanie | NGB | Massive |
|-------|----------|----------|-----|---------|
| TGP | Z₂ dyskretne | spontaneous | 0 | 1 (m_sp = √γ) |
| Ising 3D (LPA) | Z₂ dyskretne | spontaneous | 0 | 1 (massive scalar near T_c) |
| φ⁴ tanh-kink (1D) | Z₂ dyskretne | spontaneous | 0 | 1 (small fluctuation around ±v) + bound state translation mode (zero mode na kinku, ale to mode kinku, nie bulk) |
| Linear σ-model O(N) | O(N) ciągłe | O(N) → O(N-1) | N-1 | 1 (σ) |
| QCD chiral | SU(2)_L × SU(2)_R | → SU(2)_V | 3 (π⁰, π±) | 1 (σ heavy) |

**Wzorzec:** TGP klastruje z 3D Ising (Z₂ discrete) — exact 1 massive scalar
spektrum, brak NGB.

### Notabene: kink jako *rozwiązanie* nie *DOF*

W TGP `Φ` w przejściu z `+1` do `-1` (np. wall ψ(z) = tanh(z·m_sp/√2))
jest *rozwiązaniem* solitonowym (sek01 warstwa 3c), nie nową masową cząstką.
Małe fluktuacje δψ wokół takiego kinku mają:

- 1 zero mode (translation kink) — nie jest fizycznym DOF, jest gauge translacji
- 1 massive mode — odpowiednik "Higgs around kink"
- Continuum modes z luką masową

W sektorze próżni (ψ_vac = +1 globalnie) wszystkie zero modes są nieobecne.

## Zero ghost modes

### Definicja

Ghost mode = mode z **odwrotnym znakiem normy/kinetic energy** w
funkcjonale akcji. Identyfikacja: kinetic prefactor < 0 dla niektórych
konfiguracji.

### W TGP — pre-existująca analiza sek08b

Z [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]]
`thm:ghost-free-soliton` + `rem:ghost-summary`:

| Sektor | Pozytywność |
|--------|-------------|
| Zmienne fundamentalne ψ | `K(ψ) = ψ⁴` ≥ 0 ∀ψ |
| Perturbacje MS-TGP | `Q_s = ψ_bg⁴ > 0` ∀ψ_bg > 0 |
| Solitony substrat | `K_sub(g) = K_geo·g²` > 0 ∀g > 0 |

**Wniosek:** dla *jakichkolwiek* perturbacji wokół fizycznego tła
ψ_bg ∈ (0, ∞), kinetic prefactor jest strict positive. **Zero ghost modes.**

### Komplementarność z formal Q_s sek08a

Sek08a v2.0 G.0 closure ma:

- `K(φ) = K_geo·φ⁴` — kinetic coupling
- `prop:vacuum-stability-G0`: `m_sp² = U_eff''(1)/K(1) = γ/K_geo > 0`

Both numerator (`U_eff''(1) > 0`) i denominator (`K(1) > 0`) **dodatnie**.
Iloraz dodatni. **No tachyon, no ghost.**

## Verification numerical

| Test | Skrypt | Wynik | Reference |
|------|--------|-------|-----------|
| Vacuum mass m_sp² > 0 | `phase2_P21_vacuum_uniqueness.py` | 4/5 PASS + bonus (tachion bug fix) | G.0 Phase 2 P21 |
| Mass spectrum lepton | `phase2_P22_mass_spectrum_verification.py` | 5/5 PASS (m_μ/m_e -0.0013%) | G.0 Phase 2 P22 |
| Soliton stability α=1 | `ex166_alpha1_vs_alpha2.py` | substratowa stabilna w [0.5, 2.0] | sek08b sssec:alpha-resolution |
| 12/12 PPN PASS | `lp6_formulation_dictionary.py` | 12/12 PASS | sek08b rem:formulation-dictionary |

## Werdykt mode counting

> **Z₂ dyskretne ⇒ liczba NGB = 0 (Goldstone theorem n/a).**
>
> **Pojedyncze pole skalarne Φ ⇒ DOF_total = 1.**
>
> **prop:vacuum-stability + prop:vacuum-stability-G0 ⇒ DOF_massive = 1
> z masą m_sp = √(γ/K_geo) > 0.**
>
> **sek08b 3-tier ghost-freedom ⇒ DOF_ghost = 0.**
>
> **Razem:** 1 DOF całkowite, 0 NGB, 1 massive scalar, 0 ghost.
>
> **Algebraicznie zamknięte. Audyt L03 pkt. (2) **ZAMKNIĘTY**.**

## Odpowiedzi na otwarte pytania

**Q (audyt L03 §"Co brakuje" pkt. 2):**
*"Z₂ jest dyskretne ⇒ nie ma Goldstone bosona. Powinien być tylko
1 massive mode. Czy jest?"*

**A:** TAK. Goldstone theorem nie aplikuje do dyskretnych grup
(`dim(G) = 0`). Pojedyncze pole skalarne Φ ma 1 DOF; po SSB Z₂ →
{+1} pozostaje 1 dynamiczne DOF z masą `m_sp = √(γ/K_geo)`. Pre-existing
`prop:vacuum-stability` (sek08_formalizm:1385) i `prop:vacuum-stability-G0`
(sek08a:889) explicit potwierdzają to formalnie.

## Cross-references

- [[README.md]] — werdykt + indeks
- [[tachyon_check_with_source.md]] — następny element (linearizacja na niejednorodnym tle)
- [[spectral_synthesis.md]] — Sturm-Liouville unifikacja
- [[../../core/sek08_formalizm/sek08_formalizm.tex]] §`prop:vacuum-stability` (lin. 1385)
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] §`prop:vacuum-stability-G0` (lin. 889)
- [[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]] §3-tier ghost-freedom (rem:ghost-summary)
- [[../../core/sek01_ontologia/sek01_ontologia.tex]] — Z₂ symmetry + warstwy 0/1/2/3
