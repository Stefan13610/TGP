# f(ψ) deeper principle — single unifying statement (T-FP)

**Data:** 2026-04-26
**Status:** OPEN
**Parent:** [[research/op-newton-momentum/M9_1_pp_P2_results.md]] (P2: triple convergence P2-C/D/E)
**Strategic ref:** [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.2
**Foundations binding:** [[TGP_FOUNDATIONS.md]] §1 (single-Φ Z₂)

---

## 1. Cel

M9.1'' P2 zamknął się z werdyktem "POZYTYWNY POSTULAT z potrójną motywacją
substratową" (P2-C/D/E wybierają tę samą formę). ALE eksplicit w P2 §6.3:

> **"Czy istnieje fundamentalna zasada substratowa, z której automatycznie
> wynikają P2-C, P2-D i P2-E? OTWARTY."**

T-FP zamyka tę lukę przez identyfikację **jednej zasady**, z której
P2-C, P2-D i P2-E są bezpośrednie konsekwencje:

> **Zasada T-FP (Substrate Polynomial-Degree Normalization):**
>
> W TGP potencjał substratowy V(Φ) jest wielomianem stopnia d w Φ. Dilatacja
> czasu f(ψ) jest jedynym **dobrze określonym dimensionless rationalem**
> typu V(Φ)/Φⁿ z n = d, znormalizowanym do próżni:
> ```
>     f(ψ) = [V(Φ)/Φ^d]  /  [V(Φ_eq)/Φ_eq^d]
> ```
> Wybór `n = d = deg V` jest *unique*, bo jest to **jedyna wartość**
> dla której:
> - (a) `f(ψ)` jest **bounded** dla ψ → ∞ (asymptotic boundedness),
> - (b) `f(ψ)` jest **singular** dla ψ → 0 (non-metric phase boundary),
> - (c) zera `V(Φ)/Φ^d` ⟺ zera `V(Φ)` poza początkiem (faza-boundary inheritance).

**Twierdzenie:** Z zasady T-FP, dla V(Φ) = (γ/3)·Φ³/Φ_eq − (γ/4)·Φ⁴/Φ_eq²
(M9.1'' P2 setup), automatycznie wynika:
```
f(ψ) = (4 − 3ψ)/ψ
```
**bez żadnego dodatkowego postulatu** poza single-Φ Z₂ axiom + struktura
V(Φ) jako polynomial stopnia 4.

P2-C, P2-D, P2-E stają się **trzema niezależnymi** *konsekwencjami*
T-FP, a nie trzema independent postulatami.

---

## 2. Plan testów (T-FP.1 ... T-FP.5)

| ID | Cel | Metoda | PASS |
|----|-----|--------|------|
| **T-FP.1** | Wykazać, że n = deg(V) = 4 jest jedynym wyborem dla bounded+singular+zero-inherit | sympy: skanowanie n ∈ {2,3,4,5,6} | tylko n=4 spełnia (a)+(b)+(c) |
| **T-FP.2** | Wynik f(ψ) = (4-3ψ)/ψ przy n=4 | sympy: V/Φ⁴ normalized | identyczne z M9.1'' hyperbolic |
| **T-FP.3** | T-FP ⇒ P2-C (boundary conditions) | analytic verification | f(1)=1, f(4/3)=0, f→∞ przy ψ→0 |
| **T-FP.4** | T-FP ⇒ P2-D (dimensional naturalness) | dim analysis [V]/[Φ]^d=dimensionless | n=d=4 unique dim choice |
| **T-FP.5** | T-FP ⇒ P2-E (T⁰⁰ correspondence) | static energy density | V(Φ) = T⁰⁰_static, f tracks ΔV |

**Sukces:** 5/5 PASS → T-FP zamyka P2 §6.3 jako *closed positive*.

---

## 3. Dlaczego n = deg(V)?

### 3.1 Dimensional argument

V(Φ) jest wielomianem stopnia d w Φ. Dla każdego n:
```
V(Φ)/Φ^n = Σ c_k Φ^(k-n)   dla k = 0, 1, ..., d
```

Asymptotyczne zachowanie:
- ψ → ∞: dominuje człon o najwyższym wykładniku, czyli Φ^(d-n).
  - n < d: f → ∞ (rośnie z ψ; nie-fizyczne dla "metric coefficient").
  - n = d: f → const (asymptotic boundedness ✓).
  - n > d: f → 0 (zbyt szybkie spadanie; traci informację o V).
- ψ → 0: dominuje człon o najniższym wykładniku, czyli Φ^(0-n) (= Φ^(-n)
  jeśli V ma człon stały, lub Φ^(k_min - n) gdzie k_min ≥ 1).
  Dla TGP V(Φ) bez członu stałego (V(0)=0), dominuje k_min = 3 (V ma najniższy
  człon kubiczny Φ³/Φ_eq), więc Φ^(3-n).
  - n ≤ 3: f → 0 lub const przy ψ→0; brak singularności na granicy fazowej.
  - n ≥ 4: f → ∞ przy ψ→0 (singularność, encoding non-metric phase ✓).

Jedyna wartość spełniająca **oba** warunki (bounded przy ∞ i singular przy 0):
**n = d = 4**.

### 3.2 Algebraic uniqueness

Z V(Φ) = (γ/3)·Φ³/Φ_eq − (γ/4)·Φ⁴/Φ_eq²:
```
V/Φ⁴ = (γ/(3·Φ_eq))·Φ^(-1) − (γ/(4·Φ_eq²))
     = (γ/(12·Φ_eq²))·(4·Φ_eq/Φ − 3)
     = (γ/(12·Φ_eq²))·(4/ψ − 3)
     = (γ/(12·Φ_eq²))·(4 − 3ψ)/ψ
```

Normalizacja: `f(1) = 1` ⇒ stała `(γ/(12·Φ_eq²))` znika po podzieleniu:
```
f(ψ) = [(γ/(12Φ_eq²))(4-3ψ)/ψ] / [(γ/(12Φ_eq²))·1] = (4 − 3ψ)/ψ ✓
```

To jest **literalnie M9.1'' hyperbolic form**, derived z **jednego postulatu**.

### 3.3 Single-Φ axiom preservation

T-FP nie wprowadza nowych pól (cała zawartość: Φ, V(Φ), Φ_eq).
Nie wprowadza nowych skali (V(Φ) i Φ_eq już są w aksjomatach M9.1'').
Single-Φ Z₂ axiom zachowany.

---

## 4. Dlaczego to jest "deeper" niż P2-C/D/E

| Aspekt P2 | Aspekt T-FP |
|-----------|-------------|
| P2-C postuluje *cztery* warunki brzegowe (E1-E4) | T-FP postuluje *jeden* warunek (n = deg V) |
| P2-D postuluje "naturalność wymiarową" jako odrębną zasadę | T-FP zawiera ją jako konsekwencję (n=d wymiarowo unique) |
| P2-E pokazuje *consistency* z T⁰⁰ ale nie jest unique | T-FP daje T⁰⁰ correspondence automatycznie (V = static energy density) |
| P2 trzy independent konwergencje | T-FP jedna zasada → trzy konsekwencje |

T-FP jest **silniejsze** niż P2: każdy z P2-C/D/E jest *wynikiem* T-FP, ale
nie każda z P2-konwergencji implikuje T-FP automatycznie. Hierarchia:
```
T-FP ⟹ {P2-C, P2-D, P2-E}    (T-FP fundamental)
{P2-C ∧ P2-D ∧ P2-E}  jako warunki minimalne, ale nie unique
```

---

## 5. Czemu V(Φ) ma stopień 4?

T-FP redukuje pytanie "czemu f(ψ) = (4-3ψ)/ψ?" do pytania "czemu V(Φ) ma
stopień 4?".

Odpowiedź **już jest w aksjomatach TGP**:
- Single-Φ Z₂ axiom (TGP_FOUNDATIONS §1) ⇒ V(ŝ) musi być Z₂-parzysty
  funkcjonał ŝ.
- Lokalność + renormalizowalność (jak w polach skalarnych d=4):
  V(ŝ) = m²ŝ²/2 + λŝ⁴/4 + ... w v0 substrate
- Pivot v2 GL-substrate (M9.1'' P2 setup): V(Φ) = (β/3)Φ³ − (γ/4)Φ⁴
  z Φ = ŝ² (Z₂-parzysty zmieniły zmienne; Φ³ ⇒ ŝ⁶, Φ⁴ ⇒ ŝ⁸).
- Triple-product term H_Γ ⇒ Φ-stopień max 4 (po fluktuacjach).

Więc **deg V = 4 jest determined by aksjomaty teorii**, a `n = 4` w T-FP
jest determined by deg V. **Łącznie: f(ψ) = (4-3ψ)/ψ jest derived
end-to-end z aksjomatów teorii**, bez free parameter w postaci metryki.

To jest właściwie "single-axiom derivation" *modulo* kwartyczny kształt V.

---

## 6. Pliki

- `setup.md` (this file) — design audytu
- `f_psi_principle.py` — sympy script T-FP.1..T-FP.5
- `f_psi_principle.txt` — raw output
- `results.md` — synthesis post-execution

---

## 7. Cross-references

- [[research/op-newton-momentum/M9_1_pp_P2_results.md]] (P2 triple convergence)
- [[research/op-newton-momentum/M9_1_pp_P1_results.md]] (PPN PASS at 1PN)
- [[research/op-newton-momentum/M9_1_pp_P3_results.md]] (PPN P3 weak-field)
- [[TGP_FOUNDATIONS.md]] §1 (single-Φ Z₂)
- [[../../../tgp-core-paper/paper/tgp_core.tex]] §B (V(Φ) potential)
- [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.2 (brainstorm)
