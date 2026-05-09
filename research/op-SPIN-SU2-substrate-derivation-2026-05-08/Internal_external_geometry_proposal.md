---
title: "Internal/external geometry — rozbrojenie Derricka przez dualne ujęcie geometrii solitonu"
date: 2026-05-08
type: structural-proposal
status: WIP
parent: "[[./README.md]]"
phase: 1
tags:
  - structural-reframing
  - derrick-dissolution
  - internal-external-duality
  - relational-geometry
  - SU2-as-external
  - critical-analysis
---

# Internal / external geometry — proposal

## Source

Insight autora cyklu (sesja 2026-05-08, post-N14):

> "Geometria generowanej przestrzeni może być 'nietrywialna' wewnętrznie
> i dawać SU(2) dla obserwatora. Odległości i symetria w ramach przestrzeni
> wygenerowanej przez jeden soliton nie mają sensu bez odniesienia
> zewnętrznego obserwatora. Wewnętrzna konfiguracja przestrzeni w reżimie
> otaczającego pola powinna być widziana jako SU(2), ale nie musi a może
> nawet nie może być SU(2) we własnym układzie odniesienia."
>
> (autor cyklu, 2026-05-08)

## Centralne stwierdzenie

**SU(2) jest strukturą zewnętrznego widoku solitonu w tle Φ̄, NIE
strukturą wewnętrzną solitonu.**

Konsekwencje:
1. "Spin" nie jest własnością solitonu w izolacji
2. "Spin" jest **relacją** soliton ↔ otaczające Φ̄
3. Wewnętrzna geometria może być całkiem inna niż SU(2) — może nawet
   nie być meaningful manifold w external sense

## Część I: Krytyczna analiza idei

### I.1 Co dokładnie się twierdzi?

Można wyróżnić **trzy poziomy** twierdzenia, każdy o innej radykalności:

**Poziom 1 (słaby): SU(2) jest strukturą embedding map**
- Soliton ma jakąś wewnętrzną geometrię M_int
- Embedding M_int → M_ext (gdzie M_ext = otaczające R³ z metryką M9.1''(Φ̄))
- SU(2) charakteryzuje **klasę** embeddingów modulo equivalence
- M_int może być standardowym manifoldem; SU(2) powstaje z map structure

**Poziom 2 (średni): M_int jest nietrywialnym manifoldem**
- M_int sam ma topologię która wymusza SU(2) w external view
- Np. M_int jest non-orientable, lub ma π_1 nietrywialne
- External view "widzi" M_int przez Φ̄-mediated lens, dający SU(2)

**Poziom 3 (silny): M_int nie jest meaningful manifold w external sense**
- Wewnątrz solitonu pojęcia "odległości" i "kierunku" przestają być
  meaningful (degenerate metric? singular geometry?)
- External view SU(2) jest **emergent** z czegoś co wewnętrznie nie ma
  klasycznej struktury manifold
- To jest najbliższe analogii GR-czarnej dziury: wnętrze ma inny status
  ontologiczny niż exterior

User cytat sugeruje **Poziom 3** ("nie musi a może nawet nie może być
SU(2) we własnym układzie odniesienia").

### I.2 Czy to jest spójne z TGP?

**Argumenty ZA spójnością:**

1. **TGP_FOUNDATIONS § 5** explicite mówi że metryka jest **emergentna**
   z Φ-substratu. To oznacza że nawet "external" geometria nie jest
   fundamentalna — jest derived. Dla "internal" stosuje się to a fortiori.

2. **M9.1''(Φ)** ma horyzont przy ψ = 4/3 (Lorentzian horizon). Wewnątrz
   tego horyzontu metryka **zmienia signature** — z (-,+,+,+) na (+,-,-,-)
   lub staje się degenerate. To jest **strukturalna sygnalizacja** że
   wnętrze ma inny status.

3. **why_n3 cycle Phase 1**: ψ = 4/3 jest fixed-point R3 ODE i bariera
   strukturalna; wewnątrz wymaga osobnego traktowania.

4. **L08 kink-fermion** (audyt OPEN): emergent fermions powstają z topologii
   Φ-konfiguracji, nie z fundamentalnej SU(2) symmetry w action. To
   wspiera Poziom 1-2.

**Argumenty PRZECIW (lub: trudności):**

1. **Operacjonalna pustka.** Jeśli "wewnątrz nie ma meaningful odległości",
   jak w ogóle zdefiniować "wewnętrzną konfigurację"? Bez metryki nie
   ma framework dla equation of motion.

2. **Field theory wymaga manifoldu.** Φ jest skalarem **na manifoldzie**.
   Bez manifoldu wewnątrz, Φ-EOM się załamuje wewnątrz.

3. **Przejście wnętrze→exterior musi być smooth.** Field nie może mieć
   discontinuity. To wymaga że wewnętrzna geometria, jakkolwiek
   "nietrywialna", musi się **gładko łączyć** z external M9.1''.

4. **Co determinuje SU(2)?** Jeśli SU(2) emerguje "z embeddingu",
   konkretny mechanizm musi być wskazany. Bez tego twierdzenie jest puste.

### I.3 Czy to faktycznie dissolves Derricka?

**Naiwna odpowiedź:** TAK — Derrick mówi o static field na fixed
background; jeśli "static field on fixed background" jest nieadekwatnym
opisem (bo background jest emergent z field), Derrick scaling nie aplikuje
prosto.

**Krytyczna odpowiedź:** Czy aby na pewno?

Derrick jest **scaling argument**:
- Energia E[Φ] = T[Φ] + U[Φ] z konkretnymi scaling dimensions
- Pod Φ_λ(x) = Φ(λx) energia zmienia się
- Brak stable extremum → no stable static localized solution

Aby Derrick **nie aplikował**, potrzebne jest jedno z:
- (a) Energia E[Φ] nie istnieje w globalny sposób (no integral over R³)
- (b) Scaling Φ_λ(x) = Φ(λx) nie jest meaningful operacja
- (c) Dodatkowe terms w E [Φ] z innymi scaling dimensions

Internal/external duality **może** dać (a) lub (b):
- (a): jeśli energia "wewnątrz" nie ma sensu, total energy jest tylko
  external + boundary contribution. Może mieć inne scaling.
- (b): jeśli "x → λx" jest scaling **external coordinates**, ale **internal
  configuration jest invariantna**, scaling nie dotyka soliton.

**Ale to wymaga konkretnego matematycznego frameworku.** Sama idea bez
realizacji nie wyklucza Derricka — tylko wskazuje że **standardowy
Derrick może nie aplikować**.

### I.4 Najsilniejsza wersja proposal — formalizacja

**Hipoteza H2:** Soliton w TGP jest opisany przez **dwa coupled regimes**:

```
Region A (interior):  ψ ≥ 4/3  (lub inna threshold)
                     Geometria degeneruje / signature flips
                     Standardowy field-theoretic action niezdefiniowany

Region B (exterior):  ψ < 4/3
                     M9.1'' weak-field-like
                     Standard Φ-action stosuje się
                     Derrick aplikuje LOKALNIE w Region B
```

Boundary: ψ = 4/3 horizon — finite-area surface around soliton core.

**Energy decomposition:**
```
E_total = E_interior + E_boundary + E_exterior
```

- E_exterior: standardowy integral nad zewnętrzem, scaling λ⁻¹ T + λ⁻³ U
- E_boundary: surface term na horyzoncie, scaling λ⁰ (proportional to
  area · surface energy density)
- E_interior: niezdefiniowane przez standardowy formalizm

**Derrick-modified scaling:**
```
E(λ) = λ⁻¹ T_ext + λ⁻³ U_ext + λ⁰ E_boundary + E_interior(λ)
```

Jeśli E_interior(λ) jest nontrivial functional λ z odpowiednimi scaling
dimensions, można uzyskać balance → stable soliton.

**Specyficznie:** jeśli E_interior ~ λ^{+α} dla α > 0 (np. interior
"rosnie" przy compaktyfikacji external), mamy Skyrme-like balance.

**To jest WERSJA PROGRAMU**, nie dowód. Trzeba pokazać że internal
contribution faktycznie ma odpowiednie scaling.

### I.5 Operacjonalne wyzwanie: jak liczyć E_interior?

Jeśli wnętrze nie ma meaningful manifold structure, standard
∫ d³x [...] nie aplikuje. Możliwe podejścia:

**P1: Substrate counting (TGP-natywne)**
- Energia wewnątrz = "ile pulsacji substratu" jest used by interior
- Discrete count, finite
- Może mieć non-trivial λ dependence

**P2: Holographic / boundary-only**
- E_interior reducible do E_boundary (analog black hole entropy)
- Proportional to area horyzontu, nie volume
- Inny scaling dimension

**P3: Topological**
- E_interior = topologiczny invariant (winding, instanton number)
- λ-independent (Pontryagin-like)
- Stabilizes solution against rescaling collapse

**P3 jest najatrakcyjniejszy** dla TGP — daje stable soliton z natural
mechanism. Ale wymaga wskazania **co jest topologicznie konserwowane**
w TGP single-Φ framework.

## Część II: Co to robi z Derrickiem konkretnie?

### II.1 Standard Derrick re-examined

Derrick zakłada:
- (D1) Action S[Φ] = ∫d^Dx L(Φ, ∂Φ) — globalnie zdefiniowane
- (D2) Background spacetime fixed (np. R³)
- (D3) Field Φ ma jeden component (lub fixed multi-component structure)

Internal/external proposal kontestuje (D1) i (D2):
- (D1') Action może mieć contribution z interior, którego nie liczymy
  standardowymi środkami
- (D2') Background nie jest fixed — jest funkcja Φ przez M9.1''

### II.2 Czy (D2') sam wystarczy?

Nawet bez interior magic, **(D2')** sam może mieć efekt na Derrick.

Standard Derrick: pole na fixed background. M9.1''(Φ) means: gdy
rescaling Φ_λ, **background się zmienia**. Czyli Derrick scaling
**nie jest tylko rescaling field** — jest jednocześnie rescaling
background.

Sprawdzenie tego wymaga przepisania Derricka z **dynamic background**:
- T[Φ_λ] = ∫d³x √(-g_3(ψ_λ)) ½K(φ_λ) g^ij(ψ_λ) ∂_iφ_λ ∂_jφ_λ
- U[Φ_λ] = ∫d³x √(-g_3(ψ_λ)) V(φ_λ)

Z N14 inspection już wiemy: dla M9.1'' canonical form, wszystkie
zależności od ψ są **lokalne** w ψ (czyli zależą od ψ przy danym x).
Pod scaling ψ_λ(x) = ψ(λx) zmienia się tylko argument funkcji w prefactor;
po substitution y = λx scaling jest standardowy.

**Niestety: dynamic background w M9.1'' canonical NIE rescue Derricka.**

### II.3 Co to oznacza dla proposal?

Internal/external duality **wymaga prawdziwego internal contribution**
żeby rozbroić Derricka. Sam fakt "background zmienia się z field" nie
wystarczy.

To pcha nas w stronę:
- (P3) Topological invariant interior contribution
- (P2) Holographic boundary energy
- Lub: TGP ma fundamentalnie inną akcję niż canonical sek08a (HMR-1
  z N14 inspection)

## Część III: Re-interpretation SU(2) jako external structure

### III.1 Co znaczy "external observer" w TGP?

W TGP brak fundamentalnego observer (zgodnie z TGP-natywnym ujęciem
emergent everything). "External observer" operacjonalnie znaczy:
**otaczające Φ̄** (background field daleko od solitonu).

"Soliton widziany z zewnątrz" = Φ-konfiguracja taka jak ona wygląda z
perspektywy ψ → ψ̄ asymptotically.

### III.2 Jak SU(2) emerguje z external view?

**Mechanizm proponowany:**

Soliton ma:
- **Asymptotyczne profile** δΦ(r,θ,φ) → 0 dla r → ∞
- W skończonej "size" (≤ horyzont 4/3 bariery)

External observer widzi:
- Multipole expansion δΦ wokół centrum solitonu
- Każdy multipole transformuje się reprezentacją SO(3)
- Lift na quantum sektor → SU(2) reprezentacja

**Konkretna realizacja:**
- Najniższy multipole (l=0): scalar mode → spin 0 sektor (np. mass term)
- Drugi (l=1): vector mode → spin 1 sektor
- Halfsy (j=1/2): wymagają **dodatkowej struktury** boundary

Dla j=1/2:
- Standard QM: spinor jest "fundamental representation"
- Tu: spinor jest **emergent z embedding** soliton ↔ Φ̄
- Konkretnie: jak boundary configuration solitonu transformuje się pod
  rotation w external R³

### III.3 Klucz: rotacja jako embedding map

Niech soliton ma orientation **wewnętrzną** O ∈ SO(3) (lub SU(2)).
External rotation **zewnętrznej** koordynaty R ∈ SO(3).

**Match condition:** rotacja koordynat zewn R indukuje transformację
embedding parameters:
```
soliton(R · x) = soliton'(x)
gdzie soliton' jest tym samym solitonem, ale z O' = ρ(R) · O
```

Pytanie: czy ρ(R) jest reprezentacją SO(3) na "wewnętrznym" O?

- Jeśli O ∈ SO(3), ρ(R)O = R · O · R⁻¹ (adjoint rep) → spin 1
- Jeśli O ∈ SU(2) lift (double cover), ρ(R) działa przez SU(2) →
  spin 1/2 możliwe

**Krytyczne pytanie:** dlaczego O miałby być w SU(2), a nie SO(3)?

**Możliwa odpowiedź (Poziom 3):** wewnątrz solitonu nie ma meaningful
"3D rotation"; wewnętrzna geometria jest **inna** i pozwala na lift na
SU(2) bez sprzeczności wewnętrznej.

External observer widzi "spinor structure" jako **artefakt** patrzenia
na coś co wewnętrznie nie jest klasycznym 3D obiektem.

### III.4 Analogie referencyjne

- **GR czarna dziura:** zewnętrzny widok ma SO(3)-symmetric Schwarzschild
  metric z konkretnymi multipoles (mass, spin); wewnętrzny widok ma inną
  geometrię (singularity, signature change za horyzontem).
- **Holografia AdS/CFT:** boundary CFT ma jedną symmetrię, bulk ma
  inną; przejście jest non-trivial.
- **Spinor w geometrii Riemannowskiej:** spinor bundle nie istnieje na
  każdej rozmaitości — wymaga specjalnej struktury (spin structure).
  W TGP: czy soliton "ma spin structure" w sensie geometrycznym?

## Część IV: Przegląd ryzyk i niebezpieczeństw idei

### Risk 1: Idea jest pusta bez formalizacji
Bez konkretnego mechanizmu, "SU(2) emerguje z external view" jest
sloganem. Wymagana **konkretna mathematical structure** — np.:
- Boundary condition matching specifically generating SU(2)
- Topological invariant taking values in SU(2)/finite subgroup
- Holographic principle dla TGP-substrate

### Risk 2: "Internal nie ma manifoldu" jest zbyt ogólne
Twierdzenie że "inside soliton geometry doesn't make sense" może być:
- (a) Konkretne stwierdzenie matematyczne (degenerate metric, singularity)
- (b) Vague claim że "things are different inside"

(b) nie wystarcza do rozstrzygnięcia Derricka. (a) wymaga konkretnego
technicznego materiału.

### Risk 3: Boundary contribution może nie wystarczyć
Nawet jeśli E_boundary > 0 i scales jak λ⁰, to musi być **odpowiednia
wielkość** żeby zbalansować scaling λ⁻¹ T_ext i λ⁻³ U_ext do stationary
point. Numerical check niezbędny.

### Risk 4: Może dać niewłaściwą reprezentację
Embedding-based SU(2) może dawać reprezentację **adjoint** (spin 1),
nie **fundamental** (spin 1/2). Trzeba pokazać że konkretnie j=1/2
emerguje, nie j=1.

### Risk 5: Circular reasoning niebezpieczeństwo
Jeśli "SU(2) emerguje z embeddingu w Φ̄", a Φ̄ sama ma jakąś rotational
strukturę (wybiera oś?), to SU(2) jest tylko **przeniesienie** symetrii
otaczającej, nie własna struktura solitonu. To wraca do problemu Stage 2:
SU(2) wymaga additional structure.

## Część V: Ścieżki konkretne formalizacji

### Path-IE-A: Boundary-only energy
- E_total ≈ E_boundary
- E_boundary = surface integral nad horyzontem ψ=4/3
- Scaling: λ⁰ (intrinsic, niezależnie od external coords)
- Stability: jeśli E_boundary ma minimum przy konkretnej "size" boundary

**Plus:** mathematical structure jest standard (surface integrals)
**Minus:** wymaga wykazania że E_boundary ma odpowiednie minimum

### Path-IE-B: Topological invariant
- Soliton ma topologiczny "charge" Q ∈ Z (lub π_n)
- Energy ∝ |Q| (BPS-like)
- Stable jeśli Q ≠ 0

**Plus:** automatycznie stable, omija Derrick
**Minus:** wymaga **multi-component** field lub external structure dla
non-trivial topology — narusza S05?

### Path-IE-C: Effective field theory boundary
- Soliton jako quantum object zdefiniowany przez **boundary EFT**
- Wewnątrz nie liczymy explicite, używamy boundary action
- Standard quantum mechanics on boundary

**Plus:** dobrze rozwinięte machinery (spinor bundles, edge modes)
**Minus:** abstract; daleko od standardowego TGP language

### Path-IE-D: TGP-native substrate counting
- Wnętrze opisane discrete substratem H_Γ
- External description emergentna z coarse-graining
- "Internal degrees of freedom" = configurations H_Γ niewidoczne z external view

**Plus:** spójne z TGP foundations
**Minus:** wymaga substantial work na S07 / sek08c (substrate-to-metric)

## Część VI: Implikacje dla cyklu

### VI.1 Co się NIE zmienia
- Six requirements pozostają (1-6 spinor signatures)
- SU(2) jest właściwym matematycznym celem
- Phase 0 balance pozostaje valid

### VI.2 Co się zmienia FUNDAMENTALNIE
1. **Phase 1 reframe:** zamiast "find static localized soliton on M9.1'',"
   szukamy "describe soliton as interior+boundary+exterior coupled regime"
2. **Phase 2 reframe:** moduli space staje się **boundary configuration
   space**, nie interior orientation
3. **Phase 3 stays:** quantum lift na podwójne nakrycie (jeśli applicable)
4. **Probability assessment update:** zależy od ścieżki Path-IE-A..D

### VI.3 Probability re-update

**Pre-N14:**
- DERIVED: 15-20%, STRUCTURAL CONDITIONAL: 35-45%, NO_GO: 30-40%, EARLY_HALT: 10-15%

**Post-N14 (M9.1'' canonical NIE rescue Derrick):**
- DERIVED: 10-15%, STRUCTURAL: 30-40%, NO_GO: 40-50%, EARLY_HALT: 15-20%

**Post-IE-proposal (jeśli idea jest valid):**
- DERIVED: zależy od ścieżki konkretnej; bez konkretu trudno oszacować
- Jeśli Path-IE-B (topological): potencjalnie 25-35% (plus wymaga rozszerzenia
  S05 lub topological structure)
- Jeśli Path-IE-A (boundary-only): potencjalnie 15-25%
- **Dolna estymacja:** post-IE re-energizes cykl, ale bez konkretnego
  mathematical advance jest tylko **rebranding** problemu

### VI.4 Rekomendowany re-frame cyklu

```
ORIGINAL: op-SPIN-SU2-substrate-derivation
  → "wyprowadzenie SU(2) z Φ-EOM"

POST-N14: ABANDON / PIVOT decision pending

POST-IE-PROPOSAL:
  → "exploration internal/external duality jako framework dla SU(2)"
  → Status: STRUCTURAL_EXPLORATION (nie pure derivation)
  → Cel: konkretna realizacja Path-IE-A/B/C/D
```

## Część VII: Recommendations

### VII.1 Honest assessment
Idea autora jest **głęboka i ważna** — uderza w sedno problemu Derricka
nie atakując narzędzi (workarounds), ale założenia (czy field-theoretic
description w fixed background jest adekwatny dla solitonu w TGP).

ALE: **bez konkretnej formalizacji to jest framework, nie rozwiązanie**.
Ryzyko że idea jest pusta lub niewystarczająca jest realne.

### VII.2 Konkretne next steps

1. **Identyfikacja struktury matematycznej** (Path-IE-A..D):
   - Najmniej pretentious: Path-IE-A (boundary-only)
   - Najsilniej TGP-natywne: Path-IE-D (substrate)
   - Najbardziej standard QM: Path-IE-C (EFT)

2. **Sprawdzić czy któraś ścieżka jest już rozpoznawana w TGP literature:**
   - sek02:189-203 (otwarty problem g_eff)
   - L08 kink-fermion (topological kinks?)
   - why_n3 (ψ=4/3 horyzont jako boundary)

3. **Phase 1 re-scope:**
   - Zamiast "Φ-EOM solitony scan", zrobić "boundary structure analysis"
   - Konkretne claim: ψ=4/3 horyzont jest **boundary** wymagana dla
     internal/external duality
   - Sprawdzić czy boundary ma odpowiedniej topology

### VII.3 Co rekomendować autorowi cyklu
Idea wymaga:
- ☐ Wybór konkretnej Path-IE-A/B/C/D (lub kombinacji)
- ☐ Identyfikacja matematycznej struktury (boundary integrals, topological
  invariants, EFT, substrate counting)
- ☐ Test czy proposed structure faktycznie omija Derricka (compute scaling)
- ☐ Test czy daje SU(2) (compute symmetry of boundary/embedding)

Bez tych konkretów, idea pozostaje **przygotowaniem koncepcyjnym** dla
przyszłego cyklu, NIE rozstrzygnięciem N14.

## Cross-references

- [[./README.md]] — overview cyklu
- [[./N14_M911_inspection.md]] — werdykt że canonical M9.1'' NIE rescue Derricka
- [[./Phase1_literature_review.md]] — Derrick i workarounds
- [[../../audyt/S07_M911_derivation/README.md]] — postulate status M9.1''
- [[../../audyt/L08_kink_fermion_closure/]] — emergent fermions z topology
- [[../why_n3/]] — ψ=4/3 horizon (potential boundary)
- [[../../core/sek02/]] — open problem dla g_eff[Φ]

## Cytat autora preserwowany

> "Geometria generowanej przestrzeni może być 'nietrywialna' wewnętrznie
> i dawać SU(2) dla obserwatora. Odległości i symetria w ramach przestrzeni
> wygenerowanej przez jeden soliton nie mają sensu bez odniesienia
> zewnętrznego obserwatora. Mechanizm musi być bardziej złożony —
> wewnętrzna konfiguracja przestrzeni w reżimie otaczającego pola powinna
> być widziana jako SU(2) ale nie musi a może nawet nie może być SU(2)
> we własnym układzie odniesienia."
>
> — autor cyklu, 2026-05-08
