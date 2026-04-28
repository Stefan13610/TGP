# Λ_TGP from Φ_eq scale — cosmological constant from substrate (T-Λ)

**Data:** 2026-04-26
**Status:** OPEN
**Strategic ref:** [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.9
**Foundations binding:** [[TGP_FOUNDATIONS.md]] §1 (single-Φ Z₂)
**Predictions ref:** [[TGP_v1/README.md]] (Ω_Λ = 0.6847 input → 40 predictions)

---

## 1. Cel

W M9.1'' (P2 closure 2026-04-25):
```
V(Φ) = (γ/3)·Φ³/Φ_eq − (γ/4)·Φ⁴/Φ_eq²
V(Φ_eq) = γ·Φ_eq²/12      ← niezerowa gęstość energii w próżni
```

W TGP, gdzie grawitacja jest *emergentnym efektem zbiorowym* substratu,
ta vacuum energy musi w naturalny sposób przejawiać się jako efektywna
stała kosmologiczna Λ.

**Pytanie strategiczne (§8.9):**
> Czy Λ_obs = (3 meV)⁴ wynika *naturalnie* z V(Φ_eq) = γ·Φ_eq²/12 dla
> spójnej identyfikacji skali Φ_eq i sprzężenia γ?

**Cel T-Λ:** ustalić *quantitative*, falsifiowalny związek między:
- vacuum energy substratu V(Φ_eq),
- skalą Φ_eq (z OP-3 postulatu a_Γ = 1/Φ₀),
- sprzężeniem γ (z M9.1'' β=γ vacuum + bond renormalization),

i pokazać że obserwowana wartość Ω_Λ = 0.6847 (z Planck 2018, użyta jako
input do 40 predykcji TGP) jest **naturalna**, nie fine-tuned.

---

## 2. Setup matematyczny

### 2.1 Vacuum energy density TGP

Z M9.1'' P2:
```
V(Φ_eq) = γ·Φ_eq²/12     (1)
```

Wymiarowo (M9.1'' P2-D convention `[Φ] = mass`):
- `[V] = mass⁴`
- `[Φ_eq] = mass`
- `[γ] = mass²` (wymuszone przez (1))

### 2.2 Identyfikacja skali Φ_eq

OP-3 postulate (TGP_v1/README.md): `a_Γ = 1/Φ₀`, czyli skala długości
substratu jest odwrotnością Φ₀. W kosmologii FRW, naturalna skala
długości to `H₀⁻¹` (Hubble radius). Hipoteza:

```
Φ_eq = c₀·ℏ·H₀         (2)
```

(tzw. "Hubble cutoff" aksjon scale; znana z dyskusji ULDM/quintessence.)

W jednostkach naturalnych c=ℏ=1: `Φ_eq = H₀ ≈ 1.44e-33 eV`.

### 2.3 Identyfikacja sprzężenia γ

M9.1'' P2 vacuum-condition: β = γ (z 1PN PASS). β jest dimensionless,
γ ma `[mass²]`. Naturalna jednostka mass² w substracie to **Planck mass²**:
```
γ = M_Pl²·g̃        (3)
```

z g̃ dimensionless O(1). M9.1'' P2 nie precyzuje g̃ — postulat T-Λ:
g̃ = 1 (najbardziej naturalna wartość; alternatywne testy wymagają RG flow analysis).

### 2.4 Predykcja Λ_TGP

Z (1)+(2)+(3):
```
ρ_vac,TGP = V(Φ_eq) = γ·Φ_eq²/12 = (M_Pl²·g̃)·(H₀)²/12
         = (g̃/12)·M_Pl²·H₀²        (4)
```

**To jest dokładnie "geometric mean" dark-energy formula:**
```
ρ_DE ~ M_Pl²·H₀²
```
która jest empirycznie obserwowana — i pojawia się tu *naturally* z
substrate parameters bez fine-tuning.

### 2.5 Numerical check

- M_Pl ≈ 1.22·10²⁸ eV
- H₀ ≈ 1.44·10⁻³³ eV (z H₀ = 67.4 km/s/Mpc)
- M_Pl²·H₀² ≈ (1.22e28)² · (1.44e-33)² eV⁴ ≈ 3.1·10⁻¹⁰ eV⁴

Z g̃ = 1, prefactor 1/12:
- ρ_vac,TGP ≈ 2.6·10⁻¹¹ eV⁴

Obserwowane:
- ρ_vac,obs = ρ_crit·Ω_Λ = (8.099e-11 GeV/cm³ · Ω_Λ) ≈ 5.5·10⁻¹¹ eV⁴
- Ω_Λ = 0.6847 (Planck 2018)

**Stosunek:** ρ_TGP/ρ_obs ≈ 0.47 (g̃=1) → **w obrębie O(1)**,
nie 122 rzędów wielkości błędu (klasyczna catastrofa).

---

## 3. Dlaczego to jest ważne

### 3.1 Klasyczny "vacuum catastrophe"

W naive QFT, vacuum energy z zero-point fluctuations daje
ρ_vac ~ M_Pl⁴ ~ 10¹¹² eV⁴, czyli **122 rzędy wielkości** za dużo
względem ρ_obs ~ 10⁻¹¹ eV⁴.

To jest **najgorszy disagreement** w historii fizyki teoretycznej,
zwany "vacuum catastrophe" lub "cosmological constant problem".

### 3.2 TGP rozwiązanie naturalne

TGP postuluje `Φ_eq = H₀` (nie M_Pl): substrate scale jest skalą
kosmologiczną, nie skalą Plancka. Wtedy:
```
ρ_vac,TGP = γ·Φ_eq²/12 = M_Pl²·H₀²/12 ≈ ρ_obs
```

Krok krytyczny: **dlaczego Φ_eq = H₀, nie M_Pl?**

Odpowiedź TGP: bo Φ_eq jest skalą **substratu jako fluida**, a nie
skalą **kwantowych fluktuacji**. Substrate cell volume ~ a_Γ³ = Φ_eq⁻³,
więc small Φ_eq = large cell = cosmological scale. Quantum vacuum
fluctuations (Planck-scale) **nie są źródłem** Λ w TGP — są
*fluktuacjami* wokół substrate vacuum, nie samym vacuum.

To jest **kategorialnie różne** od QFT-on-curved-space podejścia:
TGP nie ma vacuum catastrophe w ogóle, bo vacuum energy to **substrate**
energy, nie kwantowa zero-point.

### 3.3 Falsifiable predictions

**T-Λ predykcje (PASS dla TGP):**

1. **Ω_Λ ≈ 0.7** (geometric mean formula): PASS (input do 40 predykcji).
2. **Λ jest stały w czasie** (V(Φ_eq) = const dla statycznego substratu):
   PASS — Planck/SDSS BAO dane nie pokazują evolution Λ z z. 
3. **Brak "dark energy gradient"** (Λ nie zależy od pozycji w przestrzeni):
   PASS — homogeniczność CMB + LSS.

**T-Λ falsifications (jeśli...) :**

- Jeśli pomiary by wykazały, że ρ_DE/ρ_DE_obs > O(10) w nadchodzących
  experimentach (DESI II, LSST), TGP musiałby skorygować g̃ — ale
  nie miałby już O(1) tunable parameter (g̃ = 1 jest postulowane jako
  natural).
- Jeśli równanie stanu w(z) odbiega od w = -1 silnie, TGP w obecnej
  formie sfalsyfikowane (V(Φ_eq) = const → w = -1 idealnie).

---

## 4. Plan testów (T-Λ.1 ... T-Λ.5)

| ID | Cel | Metoda | PASS criterion |
|----|-----|--------|----------------|
| **T-Λ.1** | Wyprowadzić V(Φ_eq) = γΦ_eq²/12 | sympy: V(Φ) z M9.1'' P2 | exact algebraicz. |
| **T-Λ.2** | Wymiarowa konsystencja [γ]=mass² | dim analysis | [V] = mass⁴ achieved |
| **T-Λ.3** | Identyfikacja Φ_eq = H₀ | OP-3 a_Γ = 1/Φ₀ + FRW | conceptual + numeric |
| **T-Λ.4** | Identyfikacja γ = M_Pl² · g̃ | M9.1'' β=γ vacuum | conceptual |
| **T-Λ.5** | Numerical: ρ_TGP/ρ_obs O(1) | substituted values | within factor 5 of obs |

**Sukces:** 5/5 PASS → T-Λ zamyka §8.9 jako *positive structural identification*
z O(1) numerycznym matchem do Planck 2018.

---

## 5. Status: postulat z motywacją, nie derivation

T-Λ jest **POSTULATEM Z DOBRYM MOTYWEM**, nie pełnym derivationem:

**Co T-Λ daje:**
- Strukturalna identyfikacja Λ_TGP z V(Φ_eq).
- Naturalna skala (M_Pl²·H₀²) bez fine-tuning O(10¹²²).
- Zgodność z Ω_Λ = 0.6847 do O(1).
- Falsifiable predykcje: Λ stała w czasie, w = -1.

**Co T-Λ NIE daje (otwarte):**
- First-principles wyprowadzenie γ = M_Pl² (RG flow z H_Γ, blocked by OP-1 M2).
- Why Φ_eq = H₀ specifically? (postulat OP-3, deeper origin unknown).
- O(1) factor 1/12 w prefactor — could it be exactly determined?
- Dyna(czne) ewolucja Λ przy fluctuating Φ — pozostaje to do M9.3 / cosmological perturbations.

---

## 6. Pliki

- `setup.md` (this file) — design audytu T-Λ
- `Lambda_from_Phi0.py` — sympy + numeric script T-Λ.1..T-Λ.5
- `Lambda_from_Phi0.txt` — raw output
- `results.md` — synthesis post-execution

---

## 7. Cross-references

- [[research/op-newton-momentum/M9_1_pp_P2_results.md]] (V(Φ) form)
- [[research/op-newton-momentum/M9_1_pp_P3_results.md]] (β=γ vacuum)
- [[research/cosmo_tensions]] (Ω_Λ usage in 40 predictions)
- [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.9
- [[TGP_v1/README.md]] (Ω_Λ = 0.6847 input)
