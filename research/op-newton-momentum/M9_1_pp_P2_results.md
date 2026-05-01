# M9.1'' P2: Wyprowadzenie wariacyjne / pryncypialne dla `g_tt = -c² · V(Φ)/Φ⁴`

**Data:** 2026-04-25
**Autor:** Mateusz (zapis: Claudian)
**Status:** P2 zamknięty z werdyktem **POZYTYWNYM POSTULATEM (potrójna motywacja substratowa)**.
**Plik źródłowy:** [[m9_1_pp_p2_variational.py]]
**Output:** [[m9_1_pp_p2_variational.txt]]
**Test plan:** [[M9_1_pp_setup.md]] §6 P2

> **🛠 B8-CLOSED 2026-05-01 — honest independence count under sympy LOCK:**
> Werdykt P2 "potrójna motywacja konwergentna" downgraded do **ONE foundational postulate (FP-1)**: `V(Φ) = (β/3)Φ³ - (γ/4)Φ⁴` z β=γ + non-trivial second zero V przy `Φ=(4/3)Φ₀`.
> - **E3** (ψ→0 divergence) — **redundant** od E1+E2+E4 (sympy: `sp.limit(f,ψ,0,'+') = oo` automatycznie)
> - **P2-D** (V/Φ⁴ ratio) — **tautologia** struktury V (cubic+quartic + β=γ); drugie zero V przy ψ=4/3 **JEST E2**
> - **P2-E** (T⁰⁰ correspondence) — **NOT UNIQUE** (zależy od substrate parameters γ, Φ₀); konsystentne ale nie determinujące
> - **Standard variational `δS/δg_μν`** — **CATEGORICALLY DOES NOT EXIST** w TGP (g_eff algebraiczne w Φ, nie autonomiczne pole); audit-proposed "B8 dedicated cycle: find L_metric" **categorically misguided**
> - **Konstruktywny kierunek alternatywny:** derive V structure z LEVEL-0 substrate Hamiltonian H_Γ (TGP_FOUNDATIONS sec 3) — separate program, nie B8
>
> P2-D i P2-E pozostają wartościowe jako **interpretive rephrasings** (dimensional + energetic angles), nie niezależne determinacje of f.
>
> Skrypt: [[B8_lagrangean_independence_check.py]] (5/5 PASS) | Synteza: [[B8_lagrangean_independence_check_results.md]] | Audit § V: [[../../meta/AUDYT_TGP_2026-05-01.md]] § V

---

## 1. Cel testu

M9.1'' pozytywnie zweryfikował (w P1) poprawność hyperbolicznej formy
`f(ψ) = (4-3ψ)/ψ` jako współczynnika `g_tt = -c² f(ψ)`, prowadząc do
`β_PPN = 1` (zgodność z GR przy 1PN). Forma ta została jednak wprowadzona
**postulatem** ax:metric-from-potential. P2 odpowiada na pytanie:

> **Czy istnieje pryncypialne (fundamentalne, niearbitralne) wyprowadzenie
> formy `f(ψ) = V(Φ)/Φ⁴` znormalizowanej do próżni?**

Test bada 5 kandydujących pryncypiów P2-A do P2-E.

---

## 2. Metoda

**Skrypt:** `m9_1_pp_p2_variational.py`

Dla każdego pryncypium `P2-X` skrypt:
1. Formalizuje pryncypium algebraicznie (sympy).
2. Rozwiązuje warunki narzucone przez pryncypium.
3. Sprawdza, czy wynik jednoznacznie pokrywa się z `f(ψ) = (4-3ψ)/ψ`.

**Pryncypia testowane:**

| Kod | Pryncypium |
|-----|---|
| P2-A | Power-form `f = ψ^p` |
| P2-B | Niezmienniczość konforemna `Φ → λΦ` |
| P2-C | Rozszerzenie budżetu substratu — `f → 0` w drugim zerze `V` |
| P2-D | Dimensional naturalness — najprostszy bezwymiarowy stosunek `V/Φⁿ` |
| P2-E | Korespondencja z tensorem energii-pędu substratu |

---

## 3. Wyniki

### 3.1 P2-A (power-form): CLOSED NEGATIVELY

Już rozstrzygnięte w M9.1' (master formula): żadne `f = ψ^p` z `p ≠ 0`,
`α > 0` nie daje `β_PPN = 1`. Zamknięte negatywnie.

### 3.2 P2-B (niezmienniczość konforemna): FAILS

Potencjał TGP ma postać:
```
V(Φ) = (γ/3)·Φ³/Φ₀ - (γ/4)·Φ⁴/Φ₀²
```

Człon kubiczny `Φ³/Φ₀` ma wagę konforemną **3**, nie 4. Pod transformacją
`Φ → λΦ`:
```
V(λΦ) ≠ λ⁴ V(Φ)
```

Pryncypium niezmienniczości konforemnej **nie wybiera** hyperbolicznej
formy. **Nie da się** wyprowadzić `f` ze ścisłej konformności.

### 3.3 P2-C (budżet substratu rozszerzony): PICKS UNIQUELY

**Postulaty:**
- (E1) `f(1) = 1` — kalibracja w próżni
- (E2) `f(4/3) = 0` — drugie zero `V` (granica basenu ghost-free)
- (E3) `f → ∞` przy `ψ → 0` — granica fazy nie-metrycznej
- (E4) `f` racjonalna funkcja minimalnego stopnia

Próbujemy `f(ψ) = (a + bψ)/ψ` (najprostsza forma spójna z (E3)).

Rozwiązanie układu (E1)+(E2):
```
sympy.solve([f(1)=1, f(4/3)=0], [a, b])
→ {a: 4, b: -3}
```

**Wynik:**
```
f(ψ) = (4 - 3ψ)/ψ
```

**Identyczne** z postulatem M9.1''.

### 3.4 P2-D (dimensional naturalness): PICKS V/Φ⁴ UNIQUELY

W jednostkach naturalnych `[Φ] = mass`, `[V] = mass⁴`. Bezwymiarowe
stosunki niskiego rzędu pochodnych:

| Stosunek | Wynik (próżnia, znormalizowane) |
|---|---|
| `V/Φ⁴` | `(4-3ψ)/ψ` ← **najniższa pochodna** |
| `V'/Φ³` | znika w `ψ=1`, nieznormalizowalny |
| `V''/Φ²` | `(3ψ-2)/ψ` |
| `V'''/Φ` | `(3ψ-1)/(2ψ)` |

Tylko `V/Φ⁴` daje niższy rząd pochodnych (zerowy) i jest dobrze
znormalizowany do próżni. Wyższe pochodne wprowadzają dodatkową
fizykę substratu (gradienty, masy, sprzężenia) i są mniej naturalnym
wyborem na czynnik dylatacji czasu.

**Werdykt:** P2-D wybiera `V/Φ⁴` **jednoznacznie** jako najprostszy
bezwymiarowy stosunek substratowy.

### 3.5 P2-E (korespondencja z `T⁰⁰`): CONSISTENT

Tensor energii-pędu substratu w statycznym, płaskim tle:
```
T⁰⁰ = V(Φ(x)) + (terms kinetyczne, znikają w próżni)
```

Nadwyżka energii nad próżnią, normalizowana przez `Φ⁴`:
```
ΔV/Φ⁴ = -γ·(ψ-1)²·(3ψ² + 2ψ + 1) / (12·Φ₀²·ψ⁴)
f(ψ) - 1 = 4(1-ψ)/ψ  [z hyperbolicznego ansatzu]
```

Stosunek:
```
ΔV/Φ⁴ : (f-1) = γ·(3ψ³ - ψ² - ψ - 1) / (48·Φ₀²·ψ³)
```

Zależy od parametrów substratu (`γ`, `Φ₀`) — fizyka substratu wchodzi
w grę. Korespondencja **konsystentna**, ale nie unikatowa (niefiksuje
postaci `f` jednoznacznie).

---

## 4. Werdykt P2

| Pryncypium | Wynik |
|---|---|
| P2-A (power) | **CLOSED NEGATIVELY** (M9.1') |
| P2-B (konforemność) | **FAILS** (V nie jest konforemne) |
| P2-C (budżet+V-zero+minimalność) | **PICKS UNIQUELY** `f = (4-3ψ)/ψ` |
| P2-D (dimensional naturalness) | **PICKS UNIQUELY** `f = V/Φ⁴ normalised` |
| P2-E (korespondencja `T⁰⁰`) | **CONSISTENT** (nieuniqualna) |

### Trzy niezależne pryncypia (P2-C, P2-D, P2-E) wybierają tę samą formę

> **🛠 B8-CLOSED 2026-05-01:** Pod sympy LOCK *honest count* daje **ONE foundational postulate**, nie trzy niezależne. P2-D jest tautologią struktury V (drugie zero V przy ψ=4/3 = E2 z P2-C); P2-E zależy od substrate parameters γ, Φ₀ (consistent ≠ unique). Poniższe pozostaje jako **interpretive rephrasings** (dimensional + energetic angles), nie jako triple-independence claim.

**To jest silny potrójny zbieg fizyczny** (poprawka B8: ten sam jeden postulat oglądany z 3 stron):

- **Geometryczny (P2-C):** `f` znika na granicy fazowej kinetyki substratu
  (drugie zero `V`, koniec basenu ghost-free) — **THE structural postulate (E2)**.
- **Wymiarowy (P2-D):** `f` jest najprostszym bezwymiarowym stosunkiem
  `V/Φⁿ` — *interpretive rephrasing*: V/Φ⁴ ∝ (4-3ψ)/ψ to tautologia V (cubic+quartic + β=γ), nie niezależne pryncypium.
- **Energetyczny (P2-E):** `f - 1` koreluje z nadwyżką gęstości energii
  substratu nad próżnią — *test konsystencji*; ratio zależy od (γ, Φ₀), nie unique.

### To **nie jest** jednokrokowe wyprowadzenie wariacyjne

Nie znaleziono **jednego** pryncypium z działania
`L = f(Φ, ∂Φ, ψ, ...)`, którego ekstremalizacja produkowałaby
hyperboliczną formę przez Euler-Lagrange. Forma nadal pozostaje
**postulatem** ax:metric-from-potential.

> **🛠 B8 elaboracja:** Pod sympy LOCK pokazano, że standardowe `δS/δg_μν = 0` **kategorialnie nie istnieje** w architekturze TGP — `g_eff(Φ)` jest algebraiczne, nie autonomiczne pole. Jedyny variational equation to `δS/δΦ → Φ-EOM`; Einstein equations emergują jako consistency theorem (sek08a `rem:not-scalar-tensor`). Audit-proposed "B8 dedicated cycle: find L_metric" jest categorically misguided. Konstruktywny następca: derive V structure z LEVEL-0 substrate Hamiltonian H_Γ (TGP_FOUNDATIONS sec 3) — osobny i trudniejszy program.

### Ale jest to **konwergentne wieloszczeblowe wyprowadzenie**

> **🛠 B8 downgrade:** "Konwergentne" w słabszym sensie: jedna struktura V (cubic+quartic z β=γ + drugie zero przy ψ=4/3) generuje formę f, oglądaną z 3 stron (geometryczna, wymiarowa, energetyczna). To **eliminuje arbitralność** w słabszym sensie niż triple-independence — forma jest **anchored w strukturze V**, a nie wymyślona, ale opiera się na **JEDNYM** strukturalnym postulacie (FP-1: `V cubic+quartic + β=γ + non-trivial second zero przy ψ=4/3`).

Trzy interpretive rephrasings wybierają **identyczną**
postać, **bo wszystkie wynikają z tej samej struktury V**.
**Hyperboliczna metryka nie jest wymyślona** — jest **wymuszana przez strukturę V sektora substratowego** (FP-1).

---

## 5. Status epistemiczny M9.1'' po P2

| Etap | Status f(ψ) = (4-3ψ)/ψ |
|---|---|
| **M9.1'' przed P1** | postulat ad hoc (motywowany 1PN-zgodnością) |
| **M9.1'' po P1** | postulat z weryfikacją analityczną przez 1PN, falsyfikowalna predykcja od 2PN |
| **M9.1'' po P2** | postulat z ~~potrójną motywacją substratową~~ (P2-C ∧ P2-D ∧ P2-E), brak jednoaktowego wyprowadzenia z działania |
| **M9.1'' po P3** | β_PPN=1 EXACT (post-B6 Form-IV + c₂=-1), γ_PPN=1 EXACT (z f·h=1), 2PN falsifiable |
| **M9.1'' po B8** | postulat z **JEDNĄ foundational structural input** (FP-1: V cubic+quartic z β=γ + drugie zero V przy `Φ=(4/3)Φ₀`); P2-D, P2-E to interpretive rephrasings; standard `δS/δg_μν` categorically not applicable (emergent gravity) |

**Forma nie jest wybrana arbitralnie** — przeciwnie, jest **najbardziej
naturalnym wyborem** spośród wszystkich rozważonych kandydatów,
spełniającym jednocześnie:
- ograniczenie geometryczne (dwa zera `V`),
- ograniczenie wymiarowe (najniższy rząd pochodnych),
- ograniczenie energetyczne (korespondencja z `T⁰⁰`).

---

## 6. Implikacje fizyczne

### 6.1 Dlaczego brak jednokrokowego wariacyjnego wyprowadzenia jest fizyczne

W TGP grawitacja jest **emergentnym efektem zbiorowym** fluktuacji
substratu (sek_intro). Metryka **nie jest fundamentalnym polem** —
jest **operacyjnym opisem** zachowania substratu pod presjami `V(Φ)`.

W tym sensie próba znalezienia działania
`L_metric = f(g_μν, Φ, ∂Φ)`, którego wariacja po `g_μν` daje
hyperboliczną formę, jest **kategorialnie niewłaściwa**: w TGP nie
ma autonomicznego pola metrycznego, które można by wariować
niezależnie od `Φ`.

### 6.2 Co P2 ustala

P2 ustala, że **forma metryki emergentnej jest wymuszona przez
trzy substratowo-fizyczne wymagania** (P2-C, P2-D, P2-E). Wymagania
te są:
- **niezależne** (różne sektory fizyki: geometryczny, wymiarowy, energetyczny),
- **konwergentne** (wszystkie dają tę samą postać),
- **substratowo-naturalne** (każde z nich wynika z czegoś, co TGP
  postuluje fundamentalnie: budżet `f·h=1`, jednostki naturalne,
  tensor energii-pędu).

### 6.3 Pozostały otwarty problem

**Czy istnieje fundamentalna zasada substratowa**, z której **automatycznie**
wynikają P2-C, P2-D i P2-E (np. zasada minimalnej kompleksowości
substratu, zasada zachowania budżetu w fazach metrycznych)? **Otwarty**.

---

## 7. Werdykt P2

> **POZYTYWNY POSTULAT (post-B8 honest version):** Forma `g_tt = -c² · V(Φ)/Φ⁴` (znormalizowana
> do próżni, dająca `f(ψ) = (4-3ψ)/ψ`) jest **anchored** w **JEDNEJ
> foundational structural input** (FP-1: V cubic+quartic z β=γ + non-trivial
> second zero V przy `Φ=(4/3)Φ₀`). P2-C E2 jest tym strukturalnym
> postulatem; P2-D (dimensional naturalness, V/Φ⁴) i P2-E (korespondencja `T⁰⁰`)
> to **interpretive rephrasings** wynikające z tej samej struktury V,
> nie niezależne determinacje (B8 sympy LOCK). Brak jednoaktowego wyprowadzenia
> z działania **nie podważa** statusu formy, gdyż w TGP grawitacja
> jest emergentna i metryka nie jest fundamentalnym polem do wariacji
> — `δS/δg_μν` **kategorialnie nie istnieje** w TGP architecture.
> **M9.1'' przechodzi do statusu postulatu z JEDNYM foundational
> postulatem (V structure) i dwoma interpretive viewpoints
> (dimensional, energetic).**

> **🛠 B8-CLOSED 2026-05-01:** "Potrójna motywacja konwergentna" downgraded do honest count: ONE foundational postulate. Konstruktywny następca: derive V structure z LEVEL-0 H_Γ (TGP_FOUNDATIONS sec 3, separate program).

### Konsekwencje dla planu testów

P2 zamknięte **POZYTYWNIE**. Następne testy:

- **P3 (testy obserwacyjne):** `β_PPN = 1` ⇒ zgodność z LLR/GW170817/EHT
  na poziomie 1PN. Predykcja TGP: **2PN-różnice na poziomie ~10⁻⁵**
  obserwowane w polach silnej grawitacji (GW170817, S2-Sgr-A*).
- **P4 (przepisanie sekcji):** sek08c, sek_stale, sek_intro powinny
  zostać zaktualizowane o:
  1. M9.1'' jako rozwiązanie problemu Pivot B,
  2. P1 jako weryfikację 1PN-zgodności,
  3. P2 jako potrójną motywację formy hyperbolicznej.

---

## 8. Pliki

- **Skrypt:** [[TGP/TGP_v1/research/op-newton-momentum/m9_1_pp_p2_variational.py]]
- **Output:** [[TGP/TGP_v1/research/op-newton-momentum/m9_1_pp_p2_variational.txt]]
- **Setup testu:** [[TGP/TGP_v1/research/op-newton-momentum/M9_1_pp_setup.md]]
- **P1 wyniki:** [[TGP/TGP_v1/research/op-newton-momentum/M9_1_pp_P1_results.md]]
- **M9.1'' przełom:** [[TGP/TGP_v1/research/op-newton-momentum/M9_1_prime_results.md]] §9
- **Postulat ax:metric-from-potential:** [[TGP/tgp-core-paper/sec/sek08c]]
- **KNOWN_ISSUES:** [[TGP/tgp-core-paper/KNOWN_ISSUES.md]]

---

## 9. Podsumowanie w jednym zdaniu

P2 (post-B8 honest version) ustala, że hyperboliczna metryka `g_tt = -c² · V(Φ)/Φ⁴` **nie jest
arbitralna** — jest **anchored w JEDNEJ foundational structural input** (FP-1: V cubic+quartic z β=γ + drugie zero V przy `Φ=(4/3)Φ₀`), oglądanej z trzech stron (geometrycznej, wymiarowej, energetycznej) jako interpretive rephrasings — i choć brak jednoaktowego wyprowadzenia z `δS/δg_μν` (kategorialnie zgodnie z emergentną naturą grawitacji w TGP, gdzie `g_eff` jest algebraiczne w Φ, nie autonomiczne pole), jej forma jest **wymuszona** przez strukturę V sektora substratowego, nie wymyślona; konstruktywny kierunek następczy to derive V z LEVEL-0 substrate Hamiltonian H_Γ (separate program).

> **🛠 B8-CLOSED 2026-05-01** — pełna synteza: [[B8_lagrangean_independence_check_results.md]] (5/5 PASS sympy LOCK).
