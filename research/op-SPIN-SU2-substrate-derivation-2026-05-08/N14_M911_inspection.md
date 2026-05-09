---
title: "N14 inspection — czy M9.1''(Φ) action rescue Derricka? (z krytyką samej M9.1'')"
date: 2026-05-08
type: needs-resolution
status: WIP
parent: "[[./README.md]]"
need_id: N14
phase: 1
tags:
  - N14
  - M911
  - derrick
  - critical-inspection
  - structural-blocker
---

# N14: M9.1''(Φ) effective Lagrangian inspection — z krytyką

## Cel inspekcji

Per Phase1 literature review, Derrick's theorem (1964) wyklucza static
localized solitony w 3D ze standardowym Lagrangianem. Path A wymaga
**non-canonical kinetic structure** (Skyrme-like).

**Pytanie N14:** Czy obecna kanoniczna akcja TGP (M9.1'' + L_field z
sek08a) wprowadza taką strukturę?

**Dodatkowe pytanie krytyczne:** Czy w ogóle wolno traktować M9.1'' jako
solidną podstawę dla dalszej derivation?

## Część I: Co mówi obecna akcja TGP

### Z `sek08a_akcja_zunifikowana.tex` + `TGP_FOUNDATIONS.md` lin. 50-86

**Kanoniczna akcja TGP:**
```
S_TGP[Φ, ψ_m] = ∫ d⁴x √(-g_eff) [L_field(Φ) + L_mat(Φ, ψ_m)]
```

**Lagrangian field:**
```
L_field = ½ K(φ) g_eff^μν ∂_μφ ∂_νφ - V(φ)

K(φ) = K_geo · φ⁴            (α=2 selection)
V(φ) = (β/3) φ³ - (γ/4) φ⁴   (β = γ warunek próżni)
φ = Φ/Φ_0                    (zmienna bezwymiarowa)
```

**Metryka M9.1'':**
```
g_tt = -c₀² (4-3ψ)/ψ
g_ij = ψ/(4-3ψ) δ_ij        (izotropowa, ψ=φ)
√(-g_eff) = c₀ φ              (exact)
```

### Field equation (kanoniczna):
```
∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ₀ - γΦ³/Φ₀² = -qΦ₀ρ
```

Operator kinetyczny `D_kin[Φ] = ∇²Φ + 2(∇Φ)²/Φ = (1/3φ²)∇²(φ³)`.

## Część II: Test Derrick'a dla TGP canonical action

### Static energia w 3D

Dla static field, z M9.1'' background w soliton solution (ψ = ψ(x)):

```
H[φ] = ∫d³x √(g_3) [½ K(φ) g^ij ∂_iφ ∂_jφ + V(φ)]
```

Substitute:
```
√(g_3) = [φ/(4-3φ)]^{3/2}
g^ij = (4-3φ)/φ · δ^ij
K(φ) = K_geo · φ⁴
```

Kinetic density:
```
½ K(φ) g^ij ∂_iφ ∂_jφ = ½ K_geo φ⁴ · (4-3φ)/φ · |∇φ|²
                       = ½ K_geo φ³ (4-3φ) |∇φ|²
```

Full kinetic prefactor (with √g_3):
```
[φ/(4-3φ)]^{3/2} · ½ K_geo φ³ (4-3φ) = ½ K_geo · φ^{9/2} · (4-3φ)^{-1/2}
```

Potential prefactor:
```
[φ/(4-3φ)]^{3/2} · V(φ) = [φ/(4-3φ)]^{3/2} · γ[φ³/3 - φ⁴/4]
```

### Scaling pod φ_λ(x) = φ(λx)

```
T = ∫d³x A(φ) |∇φ|²        gdzie A(φ) = ½ K_geo φ^{9/2}/(4-3φ)^{1/2}
U = ∫d³x B(φ)              gdzie B(φ) = [φ/(4-3φ)]^{3/2} V(φ)
```

Pod rescalingiem:
```
A(φ_λ(x)) = A(φ(λx))                       (zależy od φ value, nie x explicite)
|∇φ_λ|²    = λ² |∇φ|²(λx)
B(φ_λ(x)) = B(φ(λx))
```

Substitute y = λx, d³x = λ⁻³ d³y:
```
T_λ = λ⁻³ · λ² · ∫d³y A(φ(y)) |∇φ(y)|² = λ⁻¹ T
U_λ = λ⁻³ ∫d³y B(φ(y)) = λ⁻³ U
```

**To jest IDENTYCZNE scaling jak standardowy Derrick.**

### Stationary condition

```
E(λ) = λ⁻¹ T + λ⁻³ U
dE/dλ = -λ⁻² T - 3λ⁻⁴ U
At λ=1: -T - 3U = 0  ⇒  T = -3U
```

Dla rzeczywistego solitona: T ≥ 0 strukturalnie (kinetic).

**Sprawdzenie znaku U:**

Z V(φ) = γ[φ³/3 - φ⁴/4]:
- V(0) = 0
- V(1) = γ[1/3 - 1/4] = γ/12 > 0
- V(φ→∞) → -∞ (unbounded below!)
- V'(φ) = γφ²(1-φ); zera przy φ=0 i φ=1
- V''(1) = -γ < 0 → φ=1 jest LOKALNYM MAXIMUM

**Niepokojące własności V(φ):**
1. **Unbounded below dla φ→∞** — soliton może "uciec" do nieskończoności
2. **Vacuum at φ=1 jest local maximum** — unstable!

Dla soliton z ψ ∈ [0, 4/3) (inside Lorentzian horizon), V może być
positive lub negative w zależności od profilu. Strict sign U niedeterministyczny
bez explicit ansatz.

### Werdykt Derrick

**Niezależnie od znaku U, scaling jest ten sam co standardowy Derrick.**
Mamy:
- Jeden term o scaling λ⁻¹ (kinetic)
- Jeden term o scaling λ⁻³ (potential)

Dla istnienia stable minimum trzeba **co najmniej dwóch terms o
PRZECIWNYCH scaling exponents** (np. Skyrme: λ⁻¹ z standard kinetic +
λ⁺¹ z (∂φ)⁴).

**M9.1'' canonical action NIE wprowadza takiego terminu.**

→ **N14 RESOLVED Z NEGATIVE: Path A static soliton STRUCTURAL_NO_GO
w obecnym formalizmie TGP.**

## Część III: Krytyka samej M9.1'' (per user request)

### S07 audit findings (P2 OPEN, EXT-3)

Z `[[../../audyt/S07_M911_derivation/README.md]]`:

> "M9.1'' jest *postulatem*, którego konsekwencje są sprawdzalne. To
> bliżej (P) niż (E)."

**Konkretne zarzuty S07:**

1. **Empirycznie wybrana, nie strukturalnie wyprowadzona**
   - Forma I (potęgowa, g_tt = -c²/ψ) → falsyfikowana 3·10⁴σ
   - Forma II (eksponencjalna) → równoważna do O(U)
   - Forma III (M9.1'' hyperboliczna) → "exact by construction"
   - Cytat S07: "*hipoteza została dobrana tak, by dawała poprawny PPN*"

2. **Brak fundamentu z H_Γ**
   - sek02:189-203: "*Dokładna postać tej relacji [g_eff = g_eff[Φ, ∂Φ, …]]
     jest otwartym problemem TGP*."
   - Tj. derivation z dyskretnego substratu Γ nie istnieje

3. **Filozoficzny dług**
   - "*Operacyjnie M9.1'' jest szczególnym przypadkiem skalarno-tensorowego
     ansatzu*"
   - TGP_FOUNDATIONS § 5.1 jawnie odrzuca scalar-tensor — ale operacyjnie
     jest tym samym

4. **Kruchość liczbowa**
   - Predicts (5/6)U³ deviation @ 2PN dla LIGO 3G
   - Bez fundamentu: jeśli LIGO 3G da inne, M9.1'' upada bez fallback

### Niezależna validation (why_n3 cycle)

**Kontra-argument:** `[[../why_n3/]]` Phase 1 niezależnie identyfikuje
ψ=4/3 horyzont z fixed-point analizy n=3 ODE → g₀=1.874 z 4-cyfrową
precyzją.

**Ocena:** to jest **strong structural hint**, ale **NIE pełna derivation**.
Dwa niezależne źródła zbiegają się na ψ=4/3, co sugeruje że ta wartość
ma głębsze znaczenie. Jednak struktura **całej** M9.1'' (nie tylko
horyzont) wymaga osobnego dowodu.

### Konsekwencja dla niniejszego cyklu

**Jeśli M9.1'' jest postulatem**, to:
1. Derivation H1 (SU(2) z M9.1'' substratu) byłaby derivation A z B,
   gdzie B jest postulatem → **circular reasoning**
2. Cykl op-SPIN-SU2-substrate-derivation **nie może** rzetelnie używać
   M9.1'' jako "substratu" do derivation spinora bez najpierw rozwiązania
   S07
3. Co najwyżej można pokazać: **JEŚLI M9.1'' jest correct, TO mamy SU(2)** —
   ale to jest **conditional derivation**, nie absolute structural

## Część IV: Honest assessment landscape

### Trzy możliwe stany rzeczy (HMR-1, HMR-2, HMR-3)

**HMR-1: M9.1'' canonical jest correct + TGP ma dodatkowe terms**
- Możliwe że pełna akcja TGP ma higher-derivative terms NIE included w
  current canonical L_field
- Te terms mogłyby rescue Derricka
- **Prawdopodobieństwo:** średnie — sek08 nie zawiera explicit Skyrme-like
  terms, ale sek02:189-203 explicite mówi "otwarty problem"
- **Implikacja:** N14 wymaga **derivation pełnej akcji z H_Γ substrate**
  — wykraczające poza zakres niniejszego cyklu

**HMR-2: M9.1'' canonical correct + Derrick zabija static path**
- Akcja TGP jest taka jak w sek08a, brak higher-derivative terms
- Static solitony w 3D NIE istnieją (Derrick)
- Path A STRUCTURAL_NO_GO
- Pivot konieczny: oscillon (N15) lub nowa formulacja
- **Prawdopodobieństwo:** wysokie

**HMR-3: M9.1'' jest niewłaściwa**
- Status M9.1'' jest fundamentally postulate (S07 P2 OPEN)
- Może być inny ansatz dla g_eff[Φ] który by:
  - Reprodukował PPN γ=β=1
  - Wprowadzał higher-derivative naturally
  - Pozwalał stable solitony
- **Prawdopodobieństwo:** medium — S07 explicite zostawia tę możliwość
- **Implikacja:** całą agenda TGP-spin musi czekać na S07 closure

### Probability assessment update (post-N14)

| State | Prob | Action |
|-------|------|--------|
| HMR-1 (TGP ma additional structure not yet documented) | 25-35% | Wymaga derivation z H_Γ |
| HMR-2 (canonical = standard, Derrick wins) | 45-55% | Pivot oscillon (N15) lub STRUCTURAL_NO_GO |
| HMR-3 (M9.1'' niewłaściwa) | 15-25% | Czekamy na S07 closure |

**Combined:** Path A static soliton z probability ~10-15% sukcesu w
obecnej wiedzy. **Niska.**

## Część V: Decision tree dla cyklu

```
N14 inspection RESOLVED z:
  ├─ NEGATIVE dla canonical M9.1'' action (Derrick wins)
  ├─ M9.1'' sama jest postulatem (S07 P2 OPEN)
  └─ Pełna akcja TGP MOŻE mieć dodatkowe terms (HMR-1) ale niezdokumentowane

Decision options dla niniejszego cyklu:

A) HALT cykl, czekaj na S07 closure + L08 closure
   Pros: rygor, brak budowania na niepewnym fundamencie
   Cons: blokuje SU(2) derivation indefinitely

B) PIVOT na N15 (oscillon path)
   Pros: omija Derricka explicit
   Cons: time-dependent moduli space jest skomplikowany; soliton lifetime?

C) PIVOT na conditional derivation
   "JEŚLI M9.1'' jest correct AND ma higher-derivative completion, TO SU(2)"
   Pros: proceeds w rygorze conditional
   Cons: wynik jest conditional, nie absolute

D) REFRAME cykl jako "explore landscape"
   Map possibilities (Skyrme-like extensions, oscillons, alternative metrics)
   Without committing to derivation path
   Pros: honest inquiry mode
   Cons: less concrete output

E) ABANDON Path A całkowicie
   Spin SU(2) NIE wynika z TGP w obecnym frameworku
   Honest STRUCTURAL_NO_GO declaration
```

## Rekomendacja autora cyklu (subiektywna)

**Reccomendation: D + B hybrid** (REFRAME + explore oscillons).

**Uzasadnienie:**
- (A) jest zbyt drastyczne — niektóra praca może być wykonana niezależnie
- (B) sam jest niepewny — oscillons mogą nie mieć orientation
- (C) jest brittle — conditional na 3 osobnych warunkach
- (E) jest premature — nie wszystkie ścieżki przebadane
- **(D)** pozwala honestly explore landscape **bez forsowania pre-derivation**
- **(B)** jako naturalny element (D) eksploracji

Działanie konkretne:
1. Zaktualizować README cyklu — od "derivation" do "structural exploration"
2. Status: OPEN_EXPLORATION zamiast OPEN
3. Phase 1 reframe: literature scan + landscape map zamiast direct derivation
4. Probability update na frontowych dokumentach

**Wymaga decyzji autora cyklu** — niniejszy plik jest INPUT do tej
decyzji, nie autoritative resolution.

## Otwarte pytania post-N14

1. Czy ktokolwiek (S07, sek08, papers) ma argument dlaczego TGP NIE
   ma higher-derivative terms? Jeśli nie, są one **otwartą możliwością**
2. Czy `op-g0-r3-from-canonical-projection` cycle (G.0 closure) zawiera
   dowód że pełna akcja jest exactly L_field jak w sek08a, bez
   modyfikacji?
3. Czy oscillon path z V(φ) z M9.1'' canonical jest realistyczny?
   (Wymaga osobnej analizy w N15)

## Cross-references

- [[./README.md]] — overview cyklu
- [[./Phase1_literature_review.md]] — Derrick's theorem analysis
- [[./NEEDS.md]] — N14 spec
- [[../../audyt/S07_M911_derivation/README.md]] — fundamentalna krytyka M9.1''
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]] —
  kanoniczna akcja
- [[../../TGP_FOUNDATIONS.md]] — § 3 akcja zunifikowana
- [[../op-Phi-decomposition-photon-2026-05-07/Phase1_results.md]] —
  linearized δΦ-EOM (standard Klein-Gordon, no higher-derivative)
- [[../why_n3/]] — niezależna validation ψ=4/3 horizon
