# M9 — Klasyczna dynamika Φ-pola, pęd, fenomenologia GR jako analog numeryczny

**Cykl badawczy:** następnik archiwalnego M3–M8.
**Data otwarcia:** 2026-04-25.
**Status post 2026-04-26:**
  - **M9.1''** ✅ CLOSED (3 PASS + 1 conditional + 1 open) — hiperboliczna metryka, β=γ=1
  - **M9.2** ✅ CLOSED 2026-04-26 — **5/5 PASS** (pęd, bezwładność, WEP) — [[M9_2_results.md]]
  - **M9.3** ✅ **CLOSED 2026-04-26 — 5/5 PASS** (linearyzacja, dyspersja, Peters-Mathews, polaryzacje, GW170817) — [[M9_3_results.md]]
  - **Cykl M9 (klasyczna grawitacja TGP) KOMPLETNY**
**Punkt wyjścia ontologiczny:** `TGP_FOUNDATIONS.md` (top-level), §5 (grawitacja
= efekt fluktuacji pola), §6 (pęd jako Lenz-podobny back-reakcja), §7 (M3–M8
archiwum).

**Ramy:** GR jako numeryczny analog, nie izomorfizm analityczny. Źródło dla
Φ-EOM = ρ (gęstość materii w sprzężeniu minimalnym, warstwa 3b). Klasyczna
dynamika pola w fazie złamanej `m₀² < 0`, wokół `Φ_eq[ρ]`.

---

## 1. Cele M9

Sprawdzenie liczbowe, że klasyczne rozwiązania równania pola TGP
(`prop:field-eq-from-action`, `eq:field-eq-reproduced`):

```
∇²Φ + 2 (∇Φ)²/Φ + β Φ²/Φ_0 - γ Φ³/Φ_0² = -q Φ_0 ρ          (Φ-EOM)
```

odtwarzają fenomenologię GR w precyzji obserwacji, w trzech reżimach:

- **M9.1** — statyka i potencjał Newtona/PPN.
- **M9.2** — dynamika ruchu i bezwładność (Lenz-podobnie).
- **M9.3** — promieniowanie GW z dynamicznych źródeł.

Każdy z trzech testów daje konkretny wynik liczbowy + porównanie z
binding-observation.

## 2. Tło teoretyczne

### 2.1 Akcja i Φ-EOM

Z `sek08a_akcja_zunifikowana.tex` lin. 29–171:

```
S_TGP = ∫ d⁴x √(-g_eff) [ (1/2) K(φ) g_eff^μν ∂_μφ ∂_νφ - V(φ)
                          - (q/Φ_0) φ ρ ]

K(φ) = K_geo φ⁴,    V(φ) = (β/3)φ³ - (γ/4)φ⁴,   β = γ (warunek próżni)

g_eff_μν = diag(-c_0²/φ, +φ, +φ, +φ),    √(-g_eff) = c_0 φ
```

Wariacja `δS/δΦ = 0` ⟹ Φ-EOM powyżej. Operator kinetyczny
`D_kin[Φ] = ∇²Φ + 2(∇Φ)²/Φ = (1/3φ²) ∇²(φ³)` jest **kanoniczny dla TGP**.

### 2.2 Linearyzacja wokół próżni

W zmiennej bezwymiarowej `φ = Φ/Φ_0`, `φ = 1 + ε`, `|ε| ≪ 1`:

- `O(1):` β - γ = 0 (warunek próżni, tożsamościowo spełniony)
- `O(ε):` `∇²ε - β ε = -q ρ` (równanie Helmholtza/Yukawy)
- `O(ε²):` poprawki nieliniowe od `(∇Φ)²/Φ` i potencjału

Dla **β → 0** (lub `m·r ≪ 1`, gdzie `m² = β`) — limita masless: pure Newton.
Dla skończonego β > 0: **promień Yukawy** `1/√β`, poza nim eksponencjalne
wytłumienie. Skala `√β` musi być w granicy obserwacyjnej (np. < H_0/c
dla zachowania Newtona na skalach kosmologicznych).

### 2.3 Identyfikacja z PPN (sek08c lin. 171–211)

Dla `ε(r) = 2U(r)/c_0²`, `U = GM/r`:

- `g_eff_tt = -c_0²/(1 + ε) = -c_0² + c_0² ε - c_0² ε² + ...
            = -c_0²(1 - 2U/c_0² + 4U²/c_0⁴ - ...)`
- `g_eff_rr = 1 + ε = 1 + 2U/c_0²`

Porównanie z PPN:
- `γ_PPN = 1` **dokładnie** (z formy potęgowej, sek08c lin. 189)
- `β_PPN = 2` w formie potęgowej, **β_PPN = 1** w równoważnej formie
  eksponencjalnej `g_tt = -c_0² e^(-2U)` (sek08c lin. 192–211)

Forma eksponencjalna jest preferowana fenomenologicznie. **Test M9.1
musi sprawdzić obie formy** i jakie warunki na Φ-EOM determinują
wybór.

### 2.4 Identyfikacja stałej `q`

Z `ε(r) = q M / (4π r)` (rozwiązanie Yukawy w limicie β → 0) i
`ε = 2U/c_0² = 2GM/(c_0² r)`:

```
q = 8π G / c_0²    (efektywna stała sprzężenia TGP ↔ Newton).
```

To jest **jeden parametr fenomenologiczny**, który determinuje skalę
sprzężenia źródła z polem. Nie jest aksjomatyczny — wynika z dopasowania
do Newtona w słabym polu.

## 3. M9.1 — Statyka i PPN

**Cel:** rozwiązać Φ-EOM dla statycznego sferycznie symetrycznego źródła
`ρ(r) = M δ³(x)` (lub gładkiego rozszerzonego), w słabym polu, i
zweryfikować:

- (a) `ε(r) = 2GM/(c_0² r)` w limicie β → 0 (Newton);
- (b) `γ_PPN = 1` z `g_eff_rr` na poziomie linearyzacji;
- (c) wybór między formą potęgową (β_PPN = 2) a eksponencjalną
  (β_PPN = 1) — który jest "naturalny" dla Φ-EOM przy uwzględnieniu
  poprawek O(ε²)?
- (d) skala Yukawy `1/√β` w stosunku do skal obserwacyjnych (system
  słoneczny < kpc < Mpc < H_0/c).

**Plik wykonawczy:** `nprg_m9_1_static.py`.

**Kroki:**

1. Wypisz Φ-EOM w formie sferycznie symetrycznej:
   `(1/r²)(r² ε')' - β ε + O(ε²) = -q ρ(r)`
   gdzie ε' = dε/dr. To jest ODE.
2. Rozwiązanie analityczne dla β → 0, ρ punktowe: `ε(r) = q M/(4π r)`
   z warunkiem brzegowym `ε(∞) = 0`.
3. Rozwiązanie numeryczne dla β > 0: Yukawa `ε(r) = qM/(4πr) · e^(-√β·r)`.
4. Włączenie poprawek O(ε²) z (∇Φ)²/Φ: rozwiązanie iteracyjne, sprawdzenie
   czy daje formę potęgową `1/φ` czy eksponencjalną `e^(-2U)`.
5. Wyciągnięcie γ_PPN, β_PPN z rozwiązania pełnego.
6. Porównanie z Cassini (γ_PPN = 1 do 10⁻⁵), Mercury perihelion
   (β_PPN do 10⁻⁴).

**Kryterium sukcesu:**
- (a) ✓ jeśli `ε(r)` ma poprawne 1/r asymptotyczne.
- (b) ✓ jeśli γ_PPN = 1 z linearyzacji.
- (c) Wynik **deklaratywny** — TGP daje formę potęgową lub eksponencjalną
  z dokładnego rachunku, nie wybór ad hoc.
- (d) Skala Yukawy w granicy obserwacyjnej.

## 4. M9.2 — Pęd, bezwładność, zasada równoważności

**Cel:** zweryfikować Lenz-podobny obraz pędu (TGP_FOUNDATIONS §6) na
poziomie obliczeniowym.

**Kroki:**

1. Statyczne źródło `ρ_0(x)` w spoczynku → rozwiązanie statyczne `Φ_eq(x)`
   z M9.1.
2. Małe stałe przyspieszenie `a` źródła: `ρ(x,t) = ρ_0(x - X(t))`,
   `X(t) = (1/2) a t²`. Linearyzacja Φ-EOM w `a`.
3. **Back-reakcja:** policzyć siłę pola na źródło
   `F_back = -∫ ρ ∇φ d³x` jako funkcję `a`. Wyciągnąć współczynnik
   liniowy: `F_back = -m_b a`, definicja `m_b` (masa bezwładnościowa).
4. **Masa grawitacyjna:** ze sprzężenia ρ↔Φ: `m_g = q · ∫ ρ d³x`
   (do współczynników c_0).
5. **Test zasady równoważności:** czy `m_b = m_g`?
6. **Energia radiowana:** różnica między dostarczoną energią
   przyspieszenia a kinetyczną źródła = energia uciekająca jako
   fluktuacje Φ na nieskończoność.

**Plik wykonawczy:** `m9_2_momentum.py`.

**Kryterium sukcesu:**
- (a) `m_b = m_g` z dokładnością obliczeniową.
- (b) Energia radiowana jest dodatnia i zgodna z analogiem wzoru
  Larmor/Einsteina (kwadrupol).
- (c) Stały ruch (a = 0) → F_back = 0 (Newton I).

## 5. M9.3 — Promieniowanie GW i polaryzacje

**Cel:** sprawdzić, czy linearyzowane Φ-EOM wokół `Φ_eq` z dynamicznym
kwadrupolem materii daje:
- wzór kwadrupolowy o postaci GR (Peters-Mathews) numerycznie;
- prędkość propagacji `c_GW = c_0` (zgodne z GW170817 do 10⁻¹⁵, z
  konstrukcji `g_eff` automatyczne);
- mix polaryzacji **falsyfikowalnie różny od GR** (TGP ma jedno pole,
  więc fluktuacje są "skalarne" w sensie pojedynczej funkcji
  amplitudowej, ale ich przejaw w `g_eff` poprzez nieliniowość
  `g_eff[φ]` może produkować efektywnie tensorowe wzorce w pomiarach).

**Plik wykonawczy:** `m9_3_gw.py`.

**Kroki:**

1. Linearyzacja Φ-EOM wokół `Φ_eq[ρ_static]`.
2. Źródło: `ρ(x,t) = ρ_0(x) + δρ(x,t)` z oscylującym kwadrupolem masy.
3. Daleko-polowe rozwiązanie `δφ(x,t)` w postaci fali bieżącej.
4. Odzwierciedlenie w `δg_eff_μν[δφ]` przez ekspansję `g_eff[φ_eq + δφ]`.
5. Rozkład `δg_eff` na nieprzywiedne reprezentacje SO(3) wokół kierunku
   propagacji (skalarne, wektorowe, tensorowe).
6. Wyciągnięcie wzoru promieniowania `dE/dt`: porównanie z Peters-Mathews
   `(32G/5c⁵) · ⟨Q̈_ij Q̈^ij⟩`.
7. Pomiar mixu polaryzacji: amplituda składowych skalarnych, wektorowych,
   tensorowych. Porównanie z LIGO scalar-mode bound (typowo < 5%).

**Kryterium sukcesu:**
- (a) Wzór kwadrupolowy odzyskany z dokładnością ≤ 10⁻¹.
- (b) `c_GW = c_0` automatyczne (z `g_eff` struktury).
- (c) Mix polaryzacji **w granicy LIGO** (skalarna < few %), albo —
  jeśli przewyższa — to jest **nowa falsyfikowalna predykcja TGP**, do
  testowania w przyszłych obserwacjach.

## 6. Pliki M9

- `M9_program.md` — ten plik.
- `M9_1_setup.md` — analityczny setup M9.1 (statyka, PPN).
- `m9_1_static.py` — rozwiązanie Φ-EOM dla statyki.
- `m9_1_results.txt` — raw wyniki.
- `M9_1_results.md` — werdykt M9.1.
- `M9_2_setup.md`, `m9_2_momentum.py`, `m9_2_results.txt`,
  `M9_2_results.md` — analogicznie.
- `M9_3_setup.md`, `m9_3_gw.py`, `m9_3_results.txt`,
  `M9_3_results.md` — analogicznie.

## 7. Cross-references

- `TGP_FOUNDATIONS.md` (top-level): ontologia, hierarchia, obraz pędu.
- `core/sek08a_akcja_zunifikowana.tex`: akcja, Φ-EOM, β=γ próżnia.
- `core/sek08c_metryka_z_substratu.tex`: metryka, weryfikacja PPN, wybór
  potęgowa vs eksponencjalna.
- `core/sek08_formalizm.tex`: hierarchia poziomów, sprzężenie metryczne.
- `KNOWN_ISSUES.md`: status decyzji, M3–M8 archiwum.
- `research/op1-op2-op4/M{3..8}_results.md`: poprzedni cykl, β/γ przy WF.

## 8. Bottom line po M9

Po wykonaniu M9.1, M9.2, M9.3:

- Jeśli wszystkie trzy są ✓ → **TGP fenomenologicznie zamyka grawitację
  w precyzji obserwacji**, jako numeryczny analog GR. OP-2b/OP-7 zamknięte
  w sensie reformułowanym (a nie sensie M3–M8 fixed-point).
- Jeśli M9.1 ✓ i M9.2 ✓ ale M9.3 — odchylenia w polaryzacjach > LIGO bound
  → **falsyfikacja TGP w obecnej formie**. Pivot substratu w obrębie
  jednoskładnikowego skalara Z₂ (dozwolony) lub strukturalna rewizja `g_eff`.
- Jeśli M9.1 ✗ (nie odzyskuje Newton) → **fundamentalny problem ze
  Φ-EOM** lub strukturą `g_eff`. To by była niespodzianka, biorąc pod
  uwagę że PPN γ=1 wynika z konstrukcji metryki.

Realistyczne oczekiwanie: M9.1 ✓ przejdzie łatwo (bo PPN γ=1 jest
strukturalnie wymuszone), M9.2 ✓ przejdzie po właściwej linearizacji
(Lenz-back-reakcja jest cechą każdego pola z lokalną równowagą), M9.3
będzie najtrudniejszy i potencjalnie najciekawszy.
