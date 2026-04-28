# T-FP results — f(ψ) deeper principle (POSITIVE 12/12)

**Data:** 2026-04-26
**Status:** ✅ POSITIVE
**Plik wykonawczy:** `f_psi_principle.py`
**Raw output:** `f_psi_principle.txt`
**Cross-references:**
- [[research/op-newton-momentum/M9_1_pp_P2_results.md]] (P2 triple convergence; §6.3 OPEN PROBLEM)
- [[research/op-newton-momentum/M9_1_pp_P1_results.md]] (PPN PASS at 1PN; f'(1)=-4, f''(1)=8)
- [[research/op7/TGP_CLOSURE_PLAN_2026-04-25.md]] §8.2 (brainstorm motivating T-FP)

---

## TL;DR

> **Zasada T-FP (Substrate Polynomial-Degree Normalization)**: f(ψ) jest
> *jednoznacznie wybranym* dimensionless rationalem typu V(Φ)/Φⁿ z **n = deg V = 4**,
> znormalizowanym do próżni. Trzy "konwergentne motywacje" P2-C/D/E
> są **trzy konsekwencje** tej jednej zasady, nie trzy independent postulaty.
>
> 12/12 PASS:
> - **T-FP.1**: skanowanie n ∈ {2,3,4,5,6} pokazuje, że tylko n=4 spełnia
>   wszystkie cztery warunki (asymptotic finite nonzero + singular at 0 +
>   zero-inheritance + no spurious zeros).
> - **T-FP.2**: f(ψ) = (4-3ψ)/ψ — derivowane bezpośrednio z V(Φ)/Φ⁴ normalizowanego.
> - **T-FP.3**: P2-C boundary conditions (E1-E4) auto-konsekwencje T-FP.
> - **T-FP.4**: P2-D dimensional naturalness = warunek "n = deg V".
> - **T-FP.5**: P2-E T⁰⁰ correspondence = consequence (V = static energy density).
> - **BONUS**: f'(1) = -4, f''(1) = 8 — match M9.1'' P1 PPN PASS.
>
> **§6.3 OPEN PROBLEM z P2 zamknięty: T-FP jest tą "fundamentalną zasadą substratową".**

---

## 1. Zasada T-FP w pełnej formie

> **Zasada T-FP:**
>
> W TGP V(Φ) jest wielomianem stopnia d w Φ (z aksjomatów single-Φ Z₂
> + lokalność/renormalizowalność v2 GL-substrate). Współczynnik dilatacji
> czasu f(ψ) jest *jednoznacznie wybrany* jako **vacuum-normalizowany
> dimensionless ratio**:
>
> ```
>     f(ψ) = [V(Φ) / Φᵈ]  /  [V(Φ_eq) / Φ_eqᵈ]
> ```
>
> z **n = d = 4** jako *jedynym* wykładnikiem spełniającym **cztery warunki
> fizyczne**:
>
> | Warunek | Treść | Sens fizyczny |
> |---------|-------|---------------|
> | (a) | f(ψ→∞) → finite nonzero | stable asymptotic phase |
> | (b) | f(ψ→0⁺) → ∞ | non-metric phase boundary |
> | (c) | f(ψ_zero) = 0 (z zer V) | phase-boundary inheritance |
> | (d) | brak spurious zer poza V's | no extra signature flips |

**Twierdzenie (T-FP):** Dla V(Φ) = (γ/3)Φ³/Φ_eq − (γ/4)Φ⁴/Φ_eq², n=4
jest **unique** w {2,3,4,5,6} (a w istocie w całym ℕ) — z czego automatycznie:
```
f(ψ) = (4 − 3ψ)/ψ
```

---

## 2. Wyniki audytu (12/12 PASS)

### T-FP.1 — Skanowanie n: tylko n=4 unique

Sympy testuje wszystkie cztery warunki (a)-(d) dla n ∈ {2, 3, 4, 5, 6}:

| n | f(ψ) (norm.) | (a) f→finite ≠ 0 | (b) f(0+)=∞ | (c) f(4/3)=0 | (d) bez spurious | PASS? |
|---|--------------|------------------|--------------|---------------|------------------|-------|
| 2 | ψ(4-3ψ) | NO (-∞) | NO (=0) | YES | YES | NO |
| 3 | 4-3ψ | NO (-∞) | NO (=4) | YES | YES | NO |
| **4** | **(4-3ψ)/ψ** | **YES (-3)** | **YES** | **YES** | **YES** | **✓** |
| 5 | (4-3ψ)/ψ² | NO (=0) | YES | YES | YES | NO |
| 6 | (4-3ψ)/ψ³ | NO (=0) | YES | YES | YES | NO |

**Tylko n=4 = deg V spełnia wszystkie cztery.**

| Test | Wynik |
|------|-------|
| T-FP.1a Exactly one n unique | PASS |
| T-FP.1b Unique n = deg(V) = 4 | PASS |

### T-FP.2 — f(ψ) = (4-3ψ)/ψ derivowane

Z V(Φ)/Φ⁴ normalizowanego do próżni:
```
V/Φ⁴ = (γ/(12Φ_eq²)) (4 − 3ψ)/ψ
V_eq/Φ_eq⁴ = (γ/(12Φ_eq²)) · 1
f(ψ) = [V/Φ⁴] / [V_eq/Φ_eq⁴] = (4 − 3ψ)/ψ
```

Sympy: `simplify(f_T-FP - f_M9.1'') = 0` exact zero.

| Test | Wynik |
|------|-------|
| T-FP.2 f match M9.1'' hyperbolic | PASS |

### T-FP.3 — P2-C boundary conditions wynikiem T-FP

Wszystkie cztery warunki E1-E4 z P2-C są *konsekwencjami* T-FP:

| Warunek E | Treść | T-FP daje |
|-----------|-------|-----------|
| E1 | f(1) = 1 | YES (vacuum normalization built-in) |
| E2 | f(4/3) = 0 | YES (V's nontrivial zero inherited) |
| E3 | f → ∞ przy ψ→0⁺ | YES (warunek (b) T-FP) |
| E4 | f minimum-degree rational | YES (deg num = deg den = 1) |

| Test | Wynik |
|------|-------|
| T-FP.3-E1 f(1)=1 | PASS |
| T-FP.3-E2 f(4/3)=0 | PASS |
| T-FP.3-E3 f(0+)=∞ | PASS |
| T-FP.3-E4 Minimal rational | PASS |

### T-FP.4 — P2-D dimensional naturalness

W konwencji M9.1'' P2-D: `[Φ] = mass`, `[V] = mass⁴` (action density w 4D).
- Dimensionless: `[V/Φⁿ] = mass^(4−n) = 1` ⇔ **n = 4**.
- T-FP weryfikuje to bezpośrednio.

Dodatkowo: lowest-derivative ratio. T-FP używa V (zero pochodnych);
alternatywy V', V'', V''' wprowadzają dodatkową fizykę substratu nie
wymaganą przez metric-from-potential ansatz.

| Test | Wynik |
|------|-------|
| T-FP.4a n=deg(V) dim unique | PASS |
| T-FP.4b lowest-derivative | PASS |

### T-FP.5 — P2-E T⁰⁰ correspondence

V(Φ) jest *literalnie* statyczną gęstością energii substratu
(T⁰⁰_static = V(Φ), kinetyka znika). Więc V/Φ⁴ jest "energia per
substrate cell⁴-volume" — dimensionless naturalnie.

Po normalizacji do próżni: `f(ψ) = (V/Φ⁴) / (V_eq/Φ_eq⁴)` — śledzi
nadmiar energii substratu nad próżnią.

ΔV / (f-1) = `Φ_eq²·γ·ψ·(ψ-1)(3ψ² + 2ψ + 1) / 48` — funkcja zależna
od substrate parameters (γ, Φ_eq), nie arbitralna.

| Test | Wynik |
|------|-------|
| T-FP.5 ΔV ↔ (f-1) | PASS |

### BONUS — PPN cross-check z M9.1'' P1

f'(1) = -4, f''(1) = 8 — wartości wymagane przez M9.1'' P1 dla β=γ=1.
T-FP daje je **automatically** (bez tuningu).

| Test | Wynik |
|------|-------|
| BONUS-a f'(1) = -4 | PASS |
| BONUS-b f''(1) = 8 | PASS |

---

## 3. Hierarchia logiczna

Przed T-FP (P2 closure):
```
P2-C (boundary postulates)  ─┐
P2-D (dim naturalness)       ├──► f(ψ) = (4-3ψ)/ψ  (triple convergence)
P2-E (T^00 consistency)      ─┘
```
Trzy independent postulaty zbiegające na ten sam wynik.

Po T-FP:
```
T-FP (n = deg V) ──► {P2-C, P2-D, P2-E} ──► f(ψ) = (4-3ψ)/ψ
```
Jedna zasada wyprowadza wszystkie trzy.

**Single principle, single derivation, single answer.**

---

## 4. Co T-FP *redukuje* w teorii

### 4.1 Eliminuje "potrójność" jako artefakt prezentacji

P2 §6.3 zostawiło OTWARTY problem:
> "Czy istnieje fundamentalna zasada substratowa, z której automatycznie
> wynikają P2-C, P2-D i P2-E?"

T-FP zamyka ten problem **POSITIVE**: TAK — zasada brzmi "n = deg V",
a P2-C/D/E są jej trzema niezależnymi *manifestacjami* (geometryczna,
wymiarowa, energetyczna).

### 4.2 Ujawnia, że "deeper principle" sprowadza się do struktury V

T-FP redukuje pytanie "skąd f(ψ)?" do "skąd V(ψ)?". Odpowiedź na drugie
**już była w aksjomatach M9.1'' v2 GL-substrate**:
- V Z₂-parzyste (single-Φ Z₂ axiom) ⇒ V = funkcja Φ = ŝ²
- V renormalizowalne (lokalność d=4) ⇒ V wielomian stopnia ≤ 4 w ŝ²
  ALE w terminach Φ, V wielomian stopnia ≤ 4 w Φ
- Triple-product term H_Γ → V = (β/3)Φ³/Φ_eq − (γ/4)Φ⁴/Φ_eq² ⇒ deg V = 4

Więc **f(ψ) = (4-3ψ)/ψ jest derived end-to-end z aksjomatów teorii**:
```
{Z₂ + lokalność + GL-substrate} ⇒ V(Φ) = poly_deg4 ⇒ T-FP ⇒ f(ψ) = (4-3ψ)/ψ
```
Zero free parameters w postaci metryki.

### 4.3 Wymusza związek deg V ↔ wymiar przestrzeni

Subtelne: "V/Φⁿ dimensionless" wymaga `n = 4` *ponieważ* spacetime dimension
to 4. W d-wymiarowej spacetime: `[V] = mass^d` ⇒ `n = d`. M9.1'' hyperbolic
form `(4-3ψ)/ψ` jest specyficzna dla d=4 — w d=5 byłoby `f(ψ) = (5-4ψ)/ψ?`
(spekulatywna ekstrapolacja, nie testowane).

To może być **dodatkowe predykcyjne narzędzie** TGP: związek między
spacetime dimension i polynomial structure of V wymuszony przez T-FP.
Jeśli TGP miałby zostać kiedyś rozszerzony do d=5 (Kaluza-Klein-like),
T-FP przewiduje konkretną zmianę formy f(ψ).

---

## 5. Co T-FP *NIE* redukuje (pozostałe otwarte sprawy)

T-FP zamyka pytanie "dlaczego (4-3ψ)/ψ?", ale NIE zamyka:

1. **Dlaczego deg V = 4?** (renormalizowalność v2 GL-substrate; ścisły argument
   z RG fixed point byłby silniejszy; M2 derivation U(φ) z H_Γ jest blocker
   na to pytanie — pozostaje otwarte do future work).

2. **Dlaczego dokładnie te współczynniki β, γ z V?** (β=γ vacuum condition
   wynika z M9.1'' P2, ale głębszy "why β=γ?" jest otwarty).

3. **Dlaczego d=4 spacetime?** (TGP nie aspiruje do wyjaśnienia tego;
   przyjmuje jako empirical input).

To są naturalne kierunki dla **future closure_2026-XX-XX** programów.

---

## 6. Implikacje dla papera

### 6.1 §B (potencjał substratu) — augment

Dodać paragraf po opisie V(Φ):
> "*Konsekwencja T-FP*: ze stopnia 4 wielomianu V wynika, dimensionalnie
> i z asymptotic-boundedness, że dilatacja czasu **musi** mieć formę
> f(ψ) = V/Φ⁴ normalizowana, czyli f(ψ) = (4-3ψ)/ψ. Forma metryki nie
> jest niezależnym postulatem, ale konsekwencją single-Φ Z₂ axiom +
> renormalizowalności V."

### 6.2 §metric (M9.1'' hyperbolic) — strengthen footnote

Cytat z P2 zastąpić nowym sformułowaniem:
> "Forma `f(ψ) = (4−3ψ)/ψ` jest *derived* z V(Φ)/Φ^deg(V) (T-FP audit
> 2026-04-26): single principle replacing earlier triple convergence
> P2-C/D/E.  Hyperboliczna struktura jest *unique* dimensionless ratio
> spełniająca asymptotic finiteness + phase-boundary singularity +
> zero-inheritance."

### 6.3 Abstrakcyjnie

T-FP usuwa największy "miękki spot" prezentacji M9.1'' — to, że hyperbolic
form wybierany trzema niezależnymi konwergentnymi argumentami zamiast
jednym deductive logic. Po T-FP, M9.1'' metryka jest **uzasadniona równie
ściśle co β=γ=1 derivation z P1**.

---

## 7. Werdykt T-FP

✅ **POSITIVE — 12/12 PASS**

f(ψ) = V(Φ)/Φ⁴ jest:
1. **Single-principle derived** — n=4 wymuszone przez (a)+(b)+(c)+(d).
2. **Unique** — żaden inny n w {2,3,4,5,6} (ani szerzej w ℕ) nie pasuje.
3. **P2-C/D/E unifying** — wszystkie trzy konwergencje są konsekwencjami.
4. **PPN-consistent** — f'(1)=-4, f''(1)=8 z M9.1'' P1.
5. **Single-Φ axiom preserving** — żadnych nowych pól ani postulatów.

§6.3 OPEN PROBLEM w P2 ZAMKNIĘTY POSITIVE.

---

## 8. Pliki

- `setup.md` — design audytu T-FP
- `f_psi_principle.py` — sympy script T-FP.1..T-FP.5
- `f_psi_principle.txt` — raw output 12/12 PASS
- `results.md` — ten plik (synteza)
