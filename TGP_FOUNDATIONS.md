# TGP — fundamenty interpretacyjne

**Cel dokumentu:** kanoniczny, jednostronicowy zapis ontologicznych
i interpretacyjnych zasad TGP, do którego sięga każda sesja agenta
przed rozpoczęciem pracy obliczeniowej. Ten dokument **nie zastępuje**
formalizmu w `core/sek08*.tex`, **definiuje natomiast jego znaczenie**.

**Data:** 2026-04-25.
**Autor:** ustalenia z dyskusji z autorem TGP (Mateusz, environmenstefan@gmail.com).
**Status:** binding dla wszystkich przyszłych prac w `research/`.

---

## 1. Ontologia TGP — jedno zdanie

> **TGP postuluje jedno fundamentalne pole skalarne `Φ` z symetrią Z₂.
> Wszystko inne — przestrzeń, czas, materia, grawitacja, interakcje,
> efektywne stopnie swobody — jest emergentne z dynamiki tego jednego
> pola.**

To jest **fundament nieruchomy**. Każda propozycja modyfikacji:
- dodania drugiego pola fundamentalnego (skalar, wektor, tensor),
- zmiany symetrii (poza Z₂),
- rezygnacji ze skalarnego charakteru,

**narusza fundament programu TGP** i jest a priori odrzucona, niezależnie
od tego, jak elegancko rozwiązuje konkretny problem techniczny. Pivot
substratu **w obrębie skalarnego Z₂** (np. zmiana potencjału, kinetyki,
struktury bondu) jest dozwolony jako narzędzie inżynieryjne.

## 2. Hierarchia formalizmu (cztery poziomy)

Patrz `rem:hierarchia-sektorow` w `core/sek08_formalizm.tex` (lin. 15–47).

| Poziom | Obiekt | Zawartość | Status |
|---|---|---|---|
| **0** | Substrat dyskretny `Γ = (V, E)` | Hamilton `H_Γ` (GL-bond, v2 2026-04-24), symetria Z₂, coarse-graining; `Φ = ⟨ŝ²⟩`, **`σ_ab = K_ab − (1/3)δ_ab Tr(K)`, `K_ab = ⟨(∂_a ŝ)(∂_b ŝ)⟩`** (gradient strain composite, OP-7 T2 2026-04-25) | (W) |
| **1** | Równanie pola Φ | Operator `D_kin`, samointerferencja `N`, akcja TGP, α=2, β=γ | (W) |
| **2** | Metryka efektywna `g_eff^μν[Φ, σ_ab]` | **Metryka hiperboliczna (M9.1'', 2026-04-25)** `g_tt = -c₀²(4-3ψ)/ψ`, emergencja w limicie GR, FRW, PPN (γ=β=1 exact) | (E) |
| **3** | Materia i pola cechowania | Sprzężenie metryczne, Dirac, U(1)×SU(2)×SU(3), generacje | (P) |

**Materia (poziom 3) sprzęga się z `Φ` wyłącznie przez `g_eff`**, nie
bezpośrednio przez `Φ` (aksjomat `ax:metric-coupling`,
`sek08_formalizm.tex` lin. 11115–11132).

## 3. Akcja zunifikowana

Z `sek08a_akcja_zunifikowana.tex` lin. 29–97:

```
S_TGP[Φ, ψ_m] = ∫ d⁴x √(-g_eff) [ L_field(Φ) + L_mat(Φ, ψ_m) ]
```

z:

- **L_field** = (1/2) K(φ) g_eff^μν ∂_μφ ∂_νφ - V(φ),
  - K(φ) = K_geo φ⁴ (z α=2 selection),
  - V(φ) = (β/3) φ³ - (γ/4) φ⁴, **β = γ** (warunek próżni),
  - φ = Φ/Φ_0 (zmienna bezwymiarowa).

- **L_mat** = -(q/Φ_0) φ ρ
  (sprzężenie minimalne źródła ρ z Φ przez czynnik φ).

- **Metryka efektywna** (forma hiperboliczna, M9.1'' 2026-04-25, izotropowa diagonalna):
  ```
  ds² = -c_0² (4-3ψ)/ψ · dt² + ψ/(4-3ψ) · δ_ij dx^i dx^j,    ψ = Φ/Φ_0
  ```
  Daje **γ_PPN = β_PPN = 1 exact** w 1PN; współczynniki PPN wyższych rzędów
  c₂=-1, c₃=+5/3, c₄=-10/3 (Schwarzschild reproduced, M9.1'' P1).

  **Forma potęgowa** `g_tt = -c²/ψ` (poprzednie M9.1) **została sfalsyfikowana**
  2026-04-25 (M9.1 T3: β_PPN=4 vs obs 1, 3·10⁴σ).

  **Forma eksponencjalna** `g_tt = -c² e^(-2U)` (pre-pivot 2026-04-24, weak-field)
  jest **równoważna M9.1'' do O(U)** ale różni się przy 2PN+:
  M9.1'' przewiduje explicit `|Δg_tt| = (5/6)U³` deviation od GR (testowalne LIGO 3G).

- **Element objętościowy:** √(-g_eff) = c_0 · φ (dokładnie).

- **Wariacja `δS_TGP/δΦ = 0`** daje równanie pola
  (`prop:field-eq-from-action`, `eq:field-eq-reproduced`):
  ```
  ∇²Φ + 2 (∇Φ)²/Φ + β Φ²/Φ_0 - γ Φ³/Φ_0² = -q Φ_0 ρ.
  ```
  Operator kinetyczny `D_kin[Φ] = ∇²Φ + 2(∇Φ)²/Φ = (1/3φ²) ∇²(φ³)`
  jest **kanoniczną postacią TGP**.

## 4. Co to jest "materia" i co to jest "źródło"

Trzy warstwy operacyjne, **każdą z innym statusem**:

| Warstwa | Co to jest | Sprzęga z Φ przez | Status |
|---|---|---|---|
| **3a — Pola materii** ψ_m (fermiony, cechowania, płyn) | Klasyczne i kwantowe pola standardowego modelu, żyjące na `g_eff` | `g_eff^μν[Φ]` w `L_mat`, **nie** `Φ` bezpośrednio | (P) |
| **3b — Gęstość ρ** | Skalarna gęstość w `L_mat = -(q/Φ_0) φ ρ` | Bezpośrednie sprzężenie minimalne (czynnik φ) | (W) — to jest źródło dla Φ-EOM |
| **3c — Kinki / defekty** (cząstka = radialny kink Φ + topologia chiralna) | **Hipoteza/roadmap** alternatywnego opisu fermionów jako struktur w samym Φ; otwarty problem propagatora Diraca | (otwarty problem) — `rem:materia-hierarchia` | (Hipoteza) |

**Operacyjne źródło dla Φ-EOM** to **ρ** (warstwa 3b). Punkt do
których wraca każdy rachunek typu Newton/PPN/GW.

Warstwa 3c (kinki) jest **długoterminową hipotezą** — TGP ma w
przyszłości pokazać, że ψ_m emerguje z topologii Φ — ale dziś **nie
zastępuje** ρ w bieżących rachunkach.

## 5. Grawitacja w TGP — co to jest, czym **nie** jest

### 5.1 Czym **nie** jest

- **Nie ma grawitonu** — ani fundamentalnego, ani kompozytowego.
- **Nie ma propagującej cząstki spin-2.**
- **Nie ma teorii skalarno-tensorowej** w stylu Branstr-Dicke / Horndeski
  (`rem:not-scalar-tensor`, `sek08a` lin. 129–151).
- **Metryka `g_eff` nie jest niezależną zmienną dynamiczną** — jest
  wynikowo zdefiniowana przez Φ.

### 5.2 Czym **jest**

> Grawitacja w TGP jest **kolektywnym efektem fluktuacji jednego pola
> substratu Φ**. Wszystkie zjawiska grawitacyjne (Newton, GR weak field,
> GW, soczewkowanie, dylatacja czasu, redshift) są **wzorcami statystycznymi
> długo-falowych zaburzeń pola Φ** wokół klasycznego rozwiązania
> równowagi `Φ_eq[ρ]`.

To jest tradycja **emergentnej grawitacji** (Sakharov 1968, Verlinde 2010,
Padmanabhan, Volovik). Brak grawitonu **nie jest** brakiem — jest
**cechą strukturalną** programu, eliminującą problem renormalizacji
grawitonu i automatycznie respektującą twierdzenie Weinberga-Wittena
(które dotyczy fundamentalnych cząstek, nie kompozytów ani efektów
kolektywnych).

### 5.3 Cel teoretyczny

> **GR w limicie**, jako analog numeryczny — nie izomorfizm analityczny.

TGP **nie musi** wyprowadzać równań Einsteina jako równanie Lagrange'a
dla `g_eff`. **Musi** dawać zgodne z GR predykcje liczbowe w odpowiednim
regimie (system słoneczny + LIGO/Virgo + EHT, w precyzji obserwacji).
**Odchylenia** od GR są **dozwolone** poza tym regimem (np. UV substrat,
kosmologiczne ekstremalne skale).

To jest **mocniejsza** pozycja teoretyczna niż "TGP = GR", nie słabsza:
otwiera drzwi falsyfikacji w skalach UV i daje konkretne predykcje
różnicujące.

## 6. Pęd i bezwładność — obraz Lenz-podobny

### 6.1 Fundamentalny obraz fizyczny (kanoniczny dla TGP)

> **Spoczynek = lokalna równowaga pola Φ wokół źródła ρ.** Przyspieszenie
> źródła łamie tę równowagę → pole reaguje zwrotnie → siła back-reakcji
> opiera się zmianie konfiguracji → manifestuje się jako bezwładność.
> Stały ruch = pole "ślizgające się" wraz ze źródłem, brak back-reakcji,
> zachowanie pędu (analog reguły Lenza).

### 6.2 Konsekwencje strukturalne

| Zjawisko | Standardowy obraz | Obraz TGP |
|---|---|---|
| **Bezwładność** | Aksjomat (Newton I) | Back-reakcja Φ-pola na zmianę konfiguracji równowagi |
| **Masa bezwładnościowa** | Niewyjaśniony parametr | Współczynnik back-reakcji, obliczalny z linearyzacji Φ-EOM wokół `Φ_eq` |
| **Masa grawitacyjna** | Sprzężenie z polem grawitacyjnym | Siła sprzężenia źródła z Φ (`q/Φ_0`) |
| **Zasada równoważności** (m_b = m_g) | Postulat / dyfeo. niezmienniczość | **Automatyczna** z konstrukcji jednopolowej (oba pochodzą z tej samej stałej q) |
| **Newton I (brak tarcia)** | Postulat | Translacyjna niezmienniczość minimum lokalnej równowagi → brak preferowanego układu odniesienia → brak back-reakcji w stanie ustalonym |
| **Newton II (F = ma)** | Postulat | Liniowa odpowiedź pola na zmianę konfiguracji równowagi |
| **Promieniowanie z przyspieszającego źródła** | Larmor (EM), kwadrupol (GR) | Energia uciekająca jako fluktuacje Φ na nieskończoność (analog) |

### 6.3 Skala mikro vs makro (interpretacja kwantowa — robocza)

W skali mikro fluktuacje pola Φ są porównywalne z lokalizacją źródła
(cząstki). "Pomiar pozycji cząstki" = pomiar centroidu fluktuującej
konfiguracji pola → emergentna nieoznaczoność `Δx · Δp ≥ ℏ/2`. To
jest stochastyczna interpretacja QM (Nelson; SED Marshalla/Boyera) —
**spekulatywna, ale spójna z architekturą TGP, do badania w późniejszych
cyklach**.

W skali makro wiele źródeł generuje uśrednione pole efektywne →
fluktuacje uśredniają się → klasyczny limit.

## 7. Status M3–M8 (archiwum)

W cyklach M3–M8 (kwiecień 2026) zbadano własności **wykładników
krytycznych** ratio `β/γ` w punkcie stałym Wilsona-Fishera 3D Isinga
metodami MK-RG (M3) i NPRG/Wetterich-LPA (M8), wraz z trzema kanałami
zamknięcia (M4 H-S Jacobian, M5 Z_Φ, M6/M7 wiązanie GL).

**Wnioski techniczne (zachowują wartość):**
- W jednoskładnikowym skalarnym Z₂ na poziomie LPA β/γ w punkcie stałym
  WF jest ujemne (≈ −0.3 do −0.5) w obu schematach RG.
- Wszystkie trzy single-channel kanały zamknięcia zostały sfalsyfikowane.

**Wnioski reinterpretacyjne (pod nową ramą "GR w limicie"):**
- M3–M8 odpowiadało na pytanie o **uniwersalność klasy krytycznej**
  (wykładniki w punkcie stałym WF), nie na pytanie o **fenomenologię
  TGP** (klasyczna dynamika pola średniego w fazie złamanej, na
  skończonej skali korelacji).
- TGP fizycznie nie żyje na T = T_c (faza krytyczna), tylko w **fazie
  złamanej** (`m₀² < 0`, `v² = |m₀²|/λ₀`), z **klasyczną dynamiką
  Φ-EOM** plus poprawki z fluktuacji.
- Wartość `β/γ` w punkcie stałym WF **nie jest** liczbą wiążącą dla
  fenomenologii TGP. **Wiążące jest:** liczbowa zgodność rozwiązań
  Φ-EOM z fenomenologią GR w odpowiednim regimie (Newton + γ/β_PPN +
  GW kwadrupol + c_GW=c).

**Status:** M3–M8 zachowane jako "dobrze wykonane na poprzednim torze",
ich wynik nie jest używany jako kryterium otwartego/zamkniętego dla
fenomenologii TGP. Sekwencja M9+ rozpoczyna nowy tor.

## 8. Aktualny program badawczy — M9+ ("klasyczna dynamika i pęd")

**Punkt wyjścia:** Φ-EOM (`eq:field-eq-reproduced`) jako **klasyczne
równanie pola w fazie złamanej**, ze źródłem `ρ` (warstwa 3b).

**Cel:** sprawdzenie liczbowe, że rozwiązania Φ-EOM odtwarzają
fenomenologię GR w precyzji obserwacji, w trzech sytuacjach:

- **M9.1 — Statyka i potencjał Newtona/PPN.** Sferyczne źródło,
  słabe pole, dopasowanie `q ↔ G`, weryfikacja `γ_PPN = 1` (dokładnie
  z konstrukcji metryki) i `β_PPN = 1` (z formy eksponencjalnej).

- **M9.2 — Pęd i bezwładność (Lenz-podobnie).** Statyczne i poruszające
  się źródło, back-reakcja pola na przyspieszenie, identyfikacja masy
  bezwładnościowej z back-reakcją, weryfikacja zasady równoważności.

- **M9.3 — Promieniowanie GW z dynamicznego źródła.** Linearyzacja
  Φ-EOM wokół `Φ_eq`, wzór kwadrupolowy, polaryzacje fluktuacji,
  porównanie z LIGO scalar-mode bound (< few %) i GW170817
  (`c_GW = c` do 10⁻¹⁵).

Pliki wykonawcze: `research/op-newton-momentum/M9_*.md` i `*.py`.

## 9. Reguły dla agentów

1. **Przed jakąkolwiek pracą obliczeniową w `research/`**: przeczytać
   ten dokument oraz odnośne `sek08*.tex`.
2. **Każda nowa hipoteza/test** musi się odnosić do hierarchii (poziom
   0/1/2/3) i deklarować, na którym poziomie operuje.
3. **Pivot substratu jest dozwolony** w obrębie aksjomatu §1 (jedno pole
   skalarne Z₂). Pivot poza ten aksjomat jest **odrzucony bez
   dyskusji** — wymaga osobnej rozmowy z autorem.
4. **GR jako numeryczny analog, nie izomorfizm.** Każdy test grawitacyjny
   formułować jako zgodność liczbową w precyzji obserwacji, nie jako
   próbę wyprowadzenia analitycznego równań Einsteina.
5. **Pęd jako Lenz-podobny back-reakcja Φ-pola.** Każda interpretacja
   inercji/masy/pędu musi być zgodna z §6 — w razie wątpliwości,
   poprosić autora o klaryfikację.
6. **Nie wprowadzać grawitonu** — ani fundamentalnego, ani
   kompozytowego (np. `σ_μν = ∂φ ⊗ ∂φ` jako "emergentny grawiton" jest
   **niezgodne** z §5: TGP nie używa cząstkowej ramy dla grawitacji).

## 10. Pliki referencyjne (canonical sources of truth)

- `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex` —
  akcja zunifikowana, wariacja → Φ-EOM, β=γ jako warunek próżni.
- `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex` —
  metryka efektywna z budżetu substratu, weryfikacja PPN.
- `core/sek08_formalizm/sek08_formalizm.tex` — pełny formalizm,
  hierarchia poziomów, sprzężenie metryczne, hipoteza materii-hierarchii.
- `KNOWN_ISSUES.md` — historia decyzji, statusy testów,
  sfalsyfikowane warianty.
- `research/op1-op2-op4/M{3..8}_*.md` — archiwum poprzedniego cyklu
  (krytyczne RG, jednoskładnikowy skalar, β/γ przy WF FP).
- `research/op-newton-momentum/M9_*.md` — bieżący cykl (klasyczna
  dynamika, pęd, GW jako analog numeryczny).

---

**Ostatnia aktualizacja:** 2026-04-25, po rozmowie ustalającej obraz
"grawitacja = efekt fluktuacji pola, brak grawitonu, pęd = Lenz-podobny
back-reakcja". Następne aktualizacje: po istotnych decyzjach
ontologicznych lub interpretacyjnych. **Nie** aktualizować po
pojedynczym wyniku liczbowym (te idą do `KNOWN_ISSUES.md` i
odnośnych `M*_results.md`).
