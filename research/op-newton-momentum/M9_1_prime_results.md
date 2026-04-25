# M9.1' — analiza pivotów po falsyfikacji w M9.1

**Cykl:** M9 (klasyczna dynamika); M9.1' = pivot-analysis po negatywnym M9.1.
**Data:** 2026-04-25.
**Solver:** `m9_1_prime_scan.py` (uogólniona Φ-EOM z parametrem α).
**Wyjście:** `m9_1_prime_scan.txt`.

Zob. `M9_1_results.md` dla wyniku M9.1 (β_PPN=4 sfalsyfikowane).
Zob. `M9_program.md` §3 dla nadrzędnego planu.

---

## 1. Kontekst i pytanie M9.1'

M9.1 zamknął test pełnej nieliniowej Φ-EOM TGP w formie kanonicznej:
- Metryka **boxed** w `sek08c` eq:metric-full-derived: `g_tt = -c²/ψ`,
  `g_rr = ψ` (forma potęgowa, p=-1).
- Operator kinetyczny **boxed** w `sek08_formalizm` thm:D-uniqueness:
  `α = 2`, K(φ) = K_geo·φ⁴ (geometryczne sprzężenie konforemne).

**Wynik M9.1: β_PPN = 4** (Mercury wymaga 1 ± 10⁻⁴ → 3·10⁴σ niezgodność).

**Pytanie M9.1':** czy minimalny pivot wewnątrz ontologii TGP
(jednopolowy Z₂ skalar Φ, brak grawitonu) odzyskuje β_PPN = 1?

## 2. Audyt status quo: które aksjomaty są WYMUSZONE?

### 2.1 sek08c (forma metryki) — status

`sek08c` zawiera **dwa niespójne wewnętrznie ansatze metryki**:

(I) **Boxed `eq:metric-full-derived` (linia 154 w sek08c):**
```
g_tt = -c²/ψ,   g_rr = ψ                  (p = -1, q = +1)
```
- Spełnia `f·h = 1` (substrate budget).
- W limicie liniowym `ε = 2U` daje γ_PPN = 1 ✓, β_PPN = 2 (potęgowe O(ε²)).
- M9.1 testował tę formę z pełną Φ-EOM → β_PPN = 4.

(II) **`thm:antipodal-uniqueness` (linia 345 w sek08c):**
```
f(ψ) = ψ^(-1/2),   h(ψ) = ψ^(+1/2)        (p = -1/2, q = +1/2)
```
- Twierdzenie wyprowadza s = 1 z `c² = c₀²·f/h = c₀²/ψ` (aksjomat ax:c)
  PLUS γ_PPN = 1.
- Forma RÓŻNA niż boxed (I).
- W remarku `rem:antipodal-implications` punkt 3 deklaruje
  „rozwinięcie do O(U²) daje β_PPN = 1 dokładnie" **bez dowodu**.
  Sprawdziliśmy bezpośrednio: forma (II) w O(ε²) daje **β_metric = 3**,
  a z dynamiką α=2 (c₂=-1) → **β_PPN = 7**, nie 1.

**Werdykt audytu §2.1:** sek08c jest wewnętrznie niespójny.
Boxed (I) i theorem (II) są dwiema RÓŻNYMI formami, oba zawodzą
β_PPN = 1 w pełnej teorii nieliniowej. Roszczenie „β_PPN=1 dokładnie"
z `rem:antipodal-implications` jest błędne (linearyzacja vs. pełne O(ε²)).

### 2.2 sek08_formalizm (kinetyka) — status

`sek08_formalizm` `thm:D-uniqueness` daje α = 2 + K = K_geo·φ⁴ z trzech
warunków:
- (C1) α = const (niezależne od φ),
- (C2) K(0) = 0 (warunek próżni N0-4),
- (C3) K(φ) = K_geo·φ⁴ (geometryczne sprzężenie z `prop:substrate-action`).

**Kluczowe: `rem:alpha2-pivot-status-pl` (datowany 2026-04-25):**
> Twierdzenie jest *wynikiem selekcji* w klasie (C1–C3) na ansatzu
> z v2-aksjomatem GL. Nie jest *wyprowadzeniem* α=2 z minimalnego
> bilinearnego substratu. Próba takiego wyprowadzenia (linia
> M2-a/b/c, M3-a/c) została **sfalsyfikowana** numerycznie i analitycznie.

Zatem **α = 2 jest WYBRANE** (warunek C3, postulat geometryczny),
**nie WYMUSZONE** ze struktury substratu. Inne α dają „inne wszechświaty"
w sensie `rem:coupling-map`, ale są dopuszczalnymi modyfikacjami TGP.

### 2.3 Stopnie swobody pivotu

Po audycie M9.1' identyfikujemy trzy parametry, które MOŻNA legalnie
modyfikować w obrębie ontologii Z₂-skalarnej:

1. **Wykładnik metryczny p** w `f(ψ) = ψ^p` (z `f·h = 1` → `h = ψ^(-p)`).
   - p = -1: boxed sek08c.
   - p = -1/2: thm:antipodal-uniqueness.
   - p = inne: nieprzebadane.
2. **Wykładnik kinetyczny α** w Φ-EOM `∇²ε + α(∇ε)²/(1+ε) - ...`.
   - α = 2: kanoniczne TGP (geometryczne).
   - α ≠ 2: zmienia C3.
3. **Aksjomat prędkości ax:c** (`c² = c₀²·f/h = c₀²/ψ`).
   - Modyfikacja ⇒ rewizja całego rozdziału `sek_stale`.

PARAMETRY ZAMROŻONE (ontologia):
- Pojedyncze pole skalarne Φ z symetrią Z₂.
- Brak grawitonu (fundamentalnego ani kompozytowego).
- Materia sprzęga się przez `g_eff` (axiom ax:metric-coupling).
- `f·h^s = 1` (zachowanie budżetu substratowego, prop:antipodal-from-budget).

## 3. Wzór master i weryfikacja numeryczna

### 3.1 Wyprowadzenie analityczne

Dla f(ψ) = ψ^p, h(ψ) = ψ^q z `f·h = 1` (s=1) wymuszającym q = -p:

**γ_PPN = -h'(1)/f'(1) = -(-p)/p = +1** (automatyczne dla każdego p).

**β_PPN = f''(1)/f'(1)² + 2c₂/f'(1) = (p-1)/p + 2c₂/p**.

Asymptotyczne rozwinięcie ε(r) = a₁/r + a₂/r² + ... w Φ-EOM próżniowej
`∇²ε + α(∇ε)²/(1+ε) = 0` (samokonsystentnie, r >> σ):

```
1/r⁴ :  2a₂ + α·a₁² = 0   ⟹   a₂ = -α·a₁²/2   ⟹   c₂ = a₂/a₁² = -α/2
```

Łącznie:

```
   β_PPN(p, α) = (p - 1 - α) / p           (dla f=ψ^p, h=ψ^-p, c₂=-α/2)

   γ_PPN = +1                              (automatyczne z f·h=1)
```

### 3.2 Werdykt analityczny

Wymóg β_PPN = 1:
```
(p - 1 - α)/p = 1   ⟺   p - 1 - α = p   ⟺   α = -1
```

**To jest niezależne od p (wykładnika metrycznego).** Każdy wybór
metryki potęgowej z f·h=1 wymaga **α = -1** dla zgodności z GR.

α = -1 oznacza K(φ) ∝ φ^(-2) (z eq:K-ode-formalizm: K = C·φ^(2α)),
co:
- **Łamie warunek N0-4** (K(0) = ∞ zamiast 0): brak fazy
  niemetrycznej, ontologia Z₂ rozpada się.
- Daje sygnaturowo niestabilne sprzężenie kinetyczne (operator nie
  ma dolnej granicy; duch).
- Sprzeczne z `prop:K0-from-substrate`: brak korelacji w fazie
  symetrycznej wymusza K(0) = 0.

**Konkluzja: w klasie metryk potęgowych z f·h=1, α=2 (lub jakiekolwiek
α>0) NIE MOŻE dać β_PPN = 1.** Pivot wymaga rezygnacji albo z (a)
formy potęgowej f, albo z f·h=1, albo z α>0.

### 3.3 Weryfikacja numeryczna (`m9_1_prime_scan.py`)

Solver z parametryzowanym α potwierdza c₂ = -α/2 z dokładnością ~0.8%:

| α | a₁ (A_eff) | a₂ | c₂ numeryczne | c₂ analityczne | błąd wzgl. |
|---|------------|-----|---------------|-----------------|------------|
| 0.0 | 0.0795 | -2.8·10⁻⁹ | -0.0000 | 0.0000 | szum (10⁻⁷) |
| 0.5 | 0.0813 | -1.64·10⁻³ | -0.2479 | -0.2500 | 0.82% |
| 1.0 | 0.0832 | -3.43·10⁻³ | -0.4959 | -0.5000 | 0.82% |
| **2.0** | **0.0872** | **-7.54·10⁻³** | **-0.9918** | **-1.0000** | **0.82%** |
| 3.0 | 0.0916 | -1.25·10⁻² | -1.4877 | -1.5000 | 0.82% |

Resztkowa różnica 0.82% to bias finite-grid (R_max=800), ten sam
mechanizm co w M9.1 §2.3 — zbiega do zera przy R_max → ∞.

**Master formula c₂ = -α/2 zweryfikowana numerycznie.** Linia α=0
(klasyczny laplasjan, c₂=0 → β_PPN = 2) reprezentuje granicę bez
korekty kinetycznej.

### 3.4 β_PPN dla wybranych (p, α=2) — numeryczne

| p | forma metryki | β_PPN analityczne | β_PPN z c₂ numerycznym |
|---|---------------|--------------------|--------------------------|
| -1 | f=1/ψ (boxed sek08c) | +4.000 | +3.984 |
| -1/2 | f=ψ⁻¹/² (antipodal) | +7.000 | +6.967 |
| -2 | f=1/ψ² | +2.500 | +2.492 |
| -1/3 | f=ψ⁻¹/³ | +10.00 | +9.951 |

**Minimum β_PPN po przestrzeni p (α=2) wynosi +2.5** (osiągalne dla
p → -∞, ale w sensie limitu: β_PPN(p,α)=1 - (1+α)/p → 1 wymaga
p → -∞ co nie jest dopuszczalną metryką bo divergent przy ψ→1).

Dokładnie: β_PPN(p, 2) = (p-3)/p = 1 - 3/p. Dla p < 0:
- p = -∞: β → 1 (nieosiągalne, divergent).
- p = -3: β = 2.
- p = -1: β = 4.
- p = -1/2: β = 7.

**Żadna skończona wartość p < 0 nie daje β_PPN = 1.**

## 4. Mapa pivotów: które łamią ontologię, które są legalne?

### 4.1 Pivot A — zmiana α (kinetyka)

α = -1: matematycznie spełnia β_PPN=1 dla każdego p. Łamie N0-4
(K(0)=0). **Wykluczone ontologicznie.**

Inne α > 0 z α ≠ 2: wymagałoby pivotu C3 (geometryczne sprzężenie
K=K_geo·φ⁴ → K=K_geo·φ^(2α)). To zmienia status `prop:substrate-action`
i `prop:K0-from-substrate`. Może być uzasadnione, **ale żadne α > 0
nie daje β_PPN=1 dla potęgowej f**.

### 4.2 Pivot B — zmiana formy f (metryki)

W obrębie f·h=1 + γ_PPN=1, dla f(ψ) ogólnego (nie potęgowego):
```
β_PPN = 1   ⟺   f''(1) = f'(1)² - 2 c₂ · f'(1)
            ⟺   f''(1) = f'(1)·(f'(1) - 2c₂)
```
Przy α=2 (c₂=-1): **f''(1) = f'(1)·(f'(1) + 2)**.

Skan kandydatów (`m9_1_prime_pivot_B.py`, wyjście
`m9_1_prime_pivot_B.txt`):

| Klasa | Forma | Wynik dla β_PPN=1 |
|-------|-------|-------------------|
| Potęgowa f=ψ^p | f'(1)=p, f''(1)=p(p-1) | tylko trywialne p=0 |
| Eksponencjalna f=exp(a·(ψ-1)) | f'(1)=a, f''(1)=a² | tylko trywialne a=0 |
| Wymierna f=1/(1+a·(ψ-1)) | f'(1)=-a, f''(1)=2a² | a=-2 (biegun ψ=3/2) |
| **Hiperboliczna f=1+a·(1-ψ)/ψ** | **f'(1)=-a, f''(1)=2a** | **a=+4** (zero ψ=4/3) |
| Logarytmiczna f=1+a·ln(ψ)+b·ln²(ψ) | f'(1)=a, f''(1)=-a+2b | b=(a²+3a)/2, 1-param |
| Wielomian f=1+a·(ψ-1)+b·(ψ-1)² | f'(1)=a, f''(1)=2b | b=a(a+2)/2, 1-param |

#### 4.2.1 Kandydat o specjalnym znaczeniu: forma (4) hiperboliczna

```
   f(ψ) = 1 + 4·(1-ψ)/ψ = (4 - 3·ψ)/ψ
   h(ψ) = 1/f(ψ) = ψ/(4 - 3·ψ)

   f(1) = 1,    f'(1) = -4,   f''(1) = +8        ✓ spełnia f''(1) = -4·(-4+2) = +8

   c² = c₀²/(...)?  -- zależnie od reinterpretacji ax:c
```

**Zero `f` przy ψ = 4/3** zbiega się dokładnie z **brzegiem basenu
ghost-free** w TGP (`prop:ghost-free-fundamental` w sek08_formalizm,
linia ~2562):

> Żadne z tych wyrażeń nie zmienia znaku w całym basenie ψ ∈ (0, 4/3).

Czyli forma (4) realizuje strukturę: **g_tt → 0 dokładnie tam, gdzie
sprzężenie kinetyczne substratu staje się patologiczne**. Geometrycznie
jest to horyzont współrzędnych (h(ψ=4/3) → ∞), interpretowalny jako
„brzeg sfery wpływu metrycznego" gdzie substrat opuszcza fazę
metryczną.

Dla porównania, forma kanoniczna f=1/ψ ma osobliwość przy ψ=0
(faza niemetryczna), ale nie ma żadnej cechy przy ψ=4/3 (kinetycznym
brzegu). Forma hiperboliczna jest „świadoma" obu granic basenu.

**Status kandydata (4):**
- Matematycznie spełnia β_PPN = γ_PPN = 1.
- Fizycznie wskazuje na dwa różne progi: kinetyczny (ψ=4/3) i
  konfiguracyjny (ψ=0).
- **Brak wyprowadzenia z fundamentu substratu** (`prop:substrate-action`,
  `prop:antipodal-from-budget`). Forma jest *post hoc* ekstrapolacją
  z dwóch znanych progów.
- Wymagałoby zrewidowania `sek08c` z nową propozycją: zamiast minimalizacji
  budżetu informacyjnego dającej f·h=1, wyprowadzić formę hiperboliczną
  z innego principium (np. „f musi mieć zero na brzegu basenu kinetycznego").

#### 4.2.2 Status Pivotu B

Istnieje co najmniej JEDNA forma hiperboliczna z naturalną interpretacją
fizyczną (forma (4)) spełniająca β_PPN = γ_PPN = 1 z α=2.

**Konsekwencje:**
- M9.1' nie zamyka jednoznacznie pivot B — istnieje matematycznie i
  fizycznie sensowny kandydat.
- ALE: ten kandydat NIE wynika z obecnego sek08c; wymaga przepisania
  `prop:antipodal-from-budget` z nową hipotezą.
- Obciążenie dowodu spoczywa na próbie wyprowadzenia formy (4) z
  niezależnych zasad substratu (np. budżet informacyjny + warunek
  brzegowy basenu kinetycznego).

To jest **otwarte zadanie badawcze** (M9.1'' — pivot B audit), poza
zakresem M9.1'.

### 4.3 Pivot C — zmiana c² axiom (ax:c)

Aksjomat `c² = c₀²·f/h = c₀²/ψ` w `thm:antipodal-uniqueness` (C1)
LUB `c_lok = c₀√f` w wcześniejszym wyprowadzeniu metryki — niespójność
między tymi dwoma definicjami w samym sek08c (zob. §2.1).

Bez aksjomatu c²: mamy tylko f·h = 1, czyli klasę 1-parametrową
(γ=1 automatyczne, niezależnie od p). Wtedy p jest wolnym parametrem
i można rozważać ogólne f. Sytuacja redukuje się do Pivotu B.

### 4.4 Pivot D — odejście od sferycznej symetrii (back-reaction
materii, momenty wyższe)

Możliwy mechanizm: stress-energy T_μν materii jako źródło Φ wyższego
rzędu. M9.1 sprzęga materię tylko przez ax:metric-coupling
(materia w geodezyjnej `g_eff`); nie ma sprzężenia stress-energy → Φ
poza skalarną gęstością ρ. Poszerzenie sprzężenia mogłoby
modyfikować c₂.

**Status:** niezbadane. Teoria nie ma obecnie mechanizmu, więc
trzeba go skonstruować ad hoc. Ryzyko: pivot poza ontologię
(stress-energy ma 10 składowych, nie 1 jak Φ — to wprowadza nowe
stopnie swobody).

### 4.5 Pivot E — zaakceptować falsyfikację

Sumaryczne:

| Pivot | Modyfikacja | Status ontologiczny | β_PPN=1? |
|-------|-------------|---------------------|----------|
| A: α | C3 thm:D-uniqueness | OK (warunkowo) | NIE (dla potęg) |
| A': α=-1 | C3 + N0-4 | **łamie próżnię** | TAK ale wykluczone |
| B: f spoza ψ^p | sek08c metryka | wymaga nowego wyprowadzenia | warunkowo TAK |
| C: ax:c | sek_stale | foundational | redukuje do B |
| D: back-react T_μν | ax:metric-coupling | ryzyko nowych d.o.f. | nieznane |
| **E: akceptacja** | OP-2b zamknięte | TGP fałszyfikowane na PN | – |

## 5. Werdykt M9.1'

**Brak minimalnego pivotu w obrębie ontologii TGP, który dawałby
β_PPN = 1 z naturalnego wyprowadzenia.**

Dokładniej:
1. W całej **klasie metryk potęgowych z f·h=1** (boxed sek08c +
   antipodal-uniqueness + dowolne p), **żadna kombinacja z α>0** nie
   daje β_PPN = 1.
2. Wartość α = -1 daje β_PPN = 1 dla każdego p, ale łamie N0-4
   (warunek próżni K(0)=0) i daje sprzężenie kinetyczne singularne
   przy φ=0 — **wyklucza Z₂-skalarną fazę niemetryczną**.
3. Pivot poza klasę potęgową f (Pivot B) jest matematycznie możliwy
   (1-parametrowa rodzina f spełnia β=γ=1 dla α=2), ale **wymaga
   uzasadnienia spoza obecnego fundamentu** — `sek08c` nie wskazuje
   takiej f bez nadbudowy postulatu.
4. Pivot stress-energy (D) jest spekulatywny i ryzykuje wprowadzenie
   stopni swobody niedostępnych w ontologii Z₂-skalarnej.

**Konkluzja silna:** TGP w obecnym sformułowaniu (sek08c + sek08_formalizm
+ ax:c + ax:metric-coupling) **jest obserwacyjnie sfalsyfikowane na
poziomie post-newtonowskim** w klasie metryk potęgowych. Mercury,
Cassini, LLR mierzą β_PPN = 1.000 ± 10⁻⁴; TGP w obrębie f=ψ^p z α=2
przewiduje β_PPN ∈ {2.5, 4, 7, ...}, nigdy 1.

**Konkluzja warunkowa (Pivot B):** Numeryczny skan rodzin kandydatów
(`m9_1_prime_pivot_B.py`) wykrył **formę hiperboliczną
f(ψ) = (4-3ψ)/ψ**, która:
- spełnia β_PPN = γ_PPN = 1 z α=2 dokładnie,
- ma zero przy ψ=4/3 — **brzegu basenu ghost-free TGP**
  (`prop:ghost-free-fundamental`),
- jest fizycznie interpretowalna jako „metryka zanika tam, gdzie
  substrat opuszcza fazę kinetycznie spójną".

**Ta forma NIE wynika z obecnego sek08c**, ale jej zgodność z dwoma
niezależnymi progami substratu (ψ=0, ψ=4/3) sugeruje że minimalna
rewizja `prop:antipodal-from-budget` z nowym warunkiem brzegowym
mogłaby ją wyprowadzić.

**M9.1'' (otwarte zadanie):** Czy forma hiperboliczna ma wyprowadzenie
z głębszego fundamentu substratu? Jeśli TAK — TGP zostaje uratowany
na poziomie post-newtonowskim. Jeśli NIE — falsyfikacja stoi.

## 6. Implikacje dla OP-2b i dalszego programu

### 6.1 OP-2b status

Po M9.1 + M9.1':
- M3–M8 (Wilson-Fisher, FP-uniwersalność): zamknięte negatywnie
  (brak grawitonu kompozytowego, β/γ < 0).
- M9.1: zamknięte negatywnie (β_PPN = 4).
- M9.1': zamknięte negatywnie w obrębie obecnego sformułowania.

**OP-2b (klasyczna grawitacja TGP) jest sfalsyfikowane na poziomie
strukturalnym.** Ratunek wymagałby pivotu zmieniającego co najmniej
jeden z:
- formę metryki (sek08c, postulat substratowy),
- sprzężenie kinetyczne (sek08_formalizm, postulat geometryczny),
- mechanizm sprzężenia materii (ax:metric-coupling),
- aksjomat prędkości światła (ax:c).

### 6.2 M9.2 i M9.3 — status

- **M9.2 (pęd, Lenz-like back-reaction):** wstrzymane.
  Bez zgodnej z GR statyki nie ma fundamentu do pęd-grawitacji.
  Sprzężenie pęd-Φ wykorzystywałoby tę samą strukturę kinetyczną
  (α=2), która zawiodła w M9.1 — niezależne testowanie dynamiki
  byłoby badaniem w sfalsyfikowanym reżimie.
- **M9.3 (fale grawitacyjne):** wstrzymane podobnie.
  Obserwacja GW170817 (c_T = c₀ z dokł. 10⁻¹⁵) i tak wyklucza
  K=φ^k dla k ≠ 4 (zob. `rem:coupling-map` w sek08_formalizm),
  więc nie jest niezależnym testem.

### 6.3 Następne kroki badawcze (poza M9)

Jeśli akceptujemy M9.1+M9.1' jako falsyfikację OP-2b:

1. **Eksternalna ocena** — kontynuować po `KNOWN_ISSUES.md`
   2026-04-25 (audyt α=2 retraction): podsumować całą sekwencję
   negatywnych wyników i ocenić, czy TGP jako program ma jeszcze
   sens.
2. **Pivot ontologii** — rozszerzenie o drugi stopień swobody
   (np. sprzężenie z polem elektromagnetycznym, wektorowy
   substrat) jest poza ontologią obecnego TGP. Stworzyłoby NOWĄ
   teorię, nie pivot TGP.
3. **Pivot postulatu substratowego** — rewizja
   `prop:substrate-action` (zamiast K=K_geo·φ⁴ inna geometria
   konformalna). To wymaga przepisania sek08_formalizm.
4. **Honest closure** — uznać OP-2b za negatywnie zamknięte i
   przepisać `tgp-core-paper` jako historię programu (TGP jako
   teoria, która zaproponowała konkretne struktury i została
   ich falsyfikowana).

## 7. Pliki

| Plik | Rola |
|------|------|
| `m9_1_prime_scan.py` | Solver z parametrem α, weryfikacja c₂=-α/2 |
| `m9_1_prime_scan.txt` | Wyjście numeryczne |
| `M9_1_prime_results.md` | Ten dokument — werdykt M9.1' |
| `M9_1_results.md` | Werdykt M9.1 (β_PPN=4) |
| `M9_1_setup.md` | Setup analityczny M9.1 |
| `m9_1_static.py` | Solver Φ-EOM kanoniczny (α=2, p=-1) |

## 7.1 Pliki Pivot B (skan kandydatów spoza f=ψ^p)

| Plik | Rola |
|------|------|
| `m9_1_prime_pivot_B.py` | Skan 6 klas funkcji f(ψ) na warunek β_PPN=1 |
| `m9_1_prime_pivot_B.txt` | Wyjście tabularyczne |

## 8. Podsumowanie jednoliniowe

> **β_PPN(p, α) = (p - 1 - α)/p** dla metryki `f(ψ)=ψ^p`, `h=ψ^-p`,
> kinetyki `α(∇ε)²/(1+ε)`. **Wymóg β_PPN=1 daje α=-1 niezależnie
> od p**, co łamie warunek próżni N0-4 i wyklucza ontologię
> Z₂-skalarną. Pivot poza klasę potęgową: forma hiperboliczna
> **f(ψ) = (4-3ψ)/ψ** matematycznie spełnia β=γ=1 z α=2 i ma zero
> przy ψ=4/3 (brzegu basenu ghost-free) — **kandydat warunkowy**
> wymagający wyprowadzenia z substratu. **W obecnym sformułowaniu
> sek08c TGP jest sfalsyfikowane; forma hiperboliczna jako pivot B
> jest otwartym problemem badawczym (M9.1'').**

---

## 9. M9.1'' — przełom: wyprowadzenie z V(Φ)/Φ⁴

**Data dodania:** 2026-04-25 (kontynuacja w tym samym dniu).
**Setup analityczny:** [[M9_1_pp_setup.md]].
**Skrypt weryfikacyjny:** `m9_1_pp_verify.py` (output: `m9_1_pp_verify.txt`).

### 9.1 Tożsamość algebraiczna

Potencjał TGP w warunku próżni (β=γ, `sek08a` `prop:vacuum-condition`):

```
   V(Φ) = (β/3)·Φ³/Φ₀ - (γ/4)·Φ⁴/Φ₀²
        = (γ/12)·Φ₀²·ψ³·(4 - 3ψ)
```

Definiując **znormalizowaną gęstość potencjału**:

```
   F(ψ) := V(Φ)/Φ⁴ = γ·(4 - 3ψ) / (12·Φ₀²·ψ)
   F(1)  = γ / (12·Φ₀²)

   f(ψ) := F(ψ)/F(1) = (4 - 3ψ)/ψ        ← TOŻSAMOŚĆ EXACT
```

To oznacza: **forma hiperboliczna z §4.2.1 NIE jest *post hoc* — to
dokładnie znormalizowane V(Φ)/Φ⁴.** Jeden krok algebraiczny od potencjału
TGP do GR-zgodnej metryki.

### 9.2 Numeryczna weryfikacja

`m9_1_pp_verify.py` wykorzystuje istniejące dane M9.1 ε(r):

```
Setup: M=q=σ=1, R_max=800, n_pts=5000
Fit asymptotyczny:  c_2 = -0.99180  (analityczny limit: -1.0)

Case (a) BOXED metryka f=1/ψ:
   β_PPN (numeric):  +3.9836
   β_PPN (c_2=-1):   +4.0000     [SFALSYFIKOWANE]

Case (b) HIPERBOLICZNA f=(4-3ψ)/ψ:
   β_PPN (numeric):  +0.9959
   β_PPN (c_2=-1):   +1.0000     [ZGODNE Z GR EXACT]
```

Resztkowe 0.4% to ten sam bias R_max=800 co w M9.1 §2.3 (zbiega do
zera w R_max → ∞). γ_PPN = 1 automatycznie (oba ansatze spełniają f·h=1).

### 9.3 Trzy progi fizyczne

Hiperboliczna forma jest „świadoma" trzech naturalnych progów substratu:

```
   ψ = 0    : V/Φ⁴ → ∞,  f → ∞     (faza niemetryczna, N0-4)
   ψ = 1    : V'(Φ) = 0, f = 1     (próżnia, kalibracja g_tt = -c²)
   ψ = 4/3  : V(Φ) = 0,  f = 0     (drugie zero V, brzeg basenu ghost-free)
```

Trzeci próg ψ=4/3 to dokładnie brzeg basenu kinetycznego z
`prop:ghost-free-fundamental` (sek08_formalizm linia ~2562). **Forma
metryki znika dokładnie tam, gdzie sprzężenie kinetyczne substratu
przestaje być dodatnio określone** — geometryczna interpretacja:
„brzeg sfery wpływu metrycznego".

### 9.4 Postulat M9.1''

Zaproponowany aksjomat (nadbudowa nad obecnym sek08c):

```
   ax:metric-from-potential
   ─────────────────────────
   g_tt(x) = -c² · (12·Φ₀²/γ) · V(Φ(x)) / Φ(x)⁴
   g_rr(x) = +δ_ij / [(12·Φ₀²/γ) · V(Φ(x)) / Φ(x)⁴]      (z f·h=1)
```

Implikacje:
- W formie znormalizowanej: f(ψ) = (4-3ψ)/ψ, h(ψ) = ψ/(4-3ψ).
- **β_PPN = γ_PPN = 1 dokładnie** (zgodne z Mercury, Cassini, LLR).
- Sprzężenie Newtona: q = 2πG/c² (4× mniejsze niż w boxed sek08c —
  szczegóły w `M9_1_pp_setup.md` §5.4).

### 9.5 Status logiczny

`ax:metric-from-potential` jest **nowym postulatem**, nie
wyprowadzeniem. Jego status badawczy:

| Test | Status | Cel |
|------|--------|-----|
| **P1**: wyższe rzędy PPN (c_3, c_4, ...) | **POZYTYWNE 2026-04-25** ([[M9_1_pp_P1_results.md]]) | spójność z GR poza β,γ |
| **P2**: wyprowadzenie wariacyjne | **POZYTYWNY POSTULAT 2026-04-25** ([[M9_1_pp_P2_results.md]]) | brak akcji jednokrokowej; potrójna motywacja substratowa (P2-C ∧ P2-D ∧ P2-E) |
| **P3**: testy obserwacyjne (LLR, GW170817) | otwarte | spójność dynamiczna |
| **P4**: rewrite sek08c, sek_stale | otwarte | dokumentacja nowej formy |

**P2 — wynik:** Trzy niezależne pryncypia substratowe (P2-C: rozszerzony budżet
z drugim zerem `V`; P2-D: dimensional naturalness `V/Φ⁴`; P2-E: korespondencja
z `T⁰⁰`) **wszystkie wybierają tę samą formę** `(4-3ψ)/ψ`. Brak jednokrokowego
wyprowadzenia z działania, ale w TGP grawitacja jest emergentna — nie ma
autonomicznego pola metrycznego do wariacji. **Postulat z potrójną motywacją
substratową jest naturalnym substytutem** dla `δS/δg_μν = 0` w teorii emergentnej.

### 9.6 Przeformułowanie werdyktu §5

**Aktualizacja §5 (po §9):**
- §5 stwierdzał: „brak minimalnego pivotu w obecnym sformułowaniu sek08c".
- §9 pokazuje: **istnieje pivot z naturalnym wyprowadzeniem algebraicznym
  z V(Φ)** (jeden krok od potencjału TGP do hiperbolicznej formy).
- Otwartym pytaniem nie jest „CZY forma istnieje" (istnieje, exact),
  lecz „CZY g_tt ∝ V/Φ⁴ jest wymuszone deeper principle (P2) lub
  ad-hoc selekcja".

**Konkluzja zaktualizowana:**
TGP w obecnym **boxed** sformułowaniu sek08c (g_tt = -c²/ψ) jest
**sfalsyfikowane** (β_PPN=4). Pivot do **g_tt = -c² · V/Φ⁴**
(po normalizacji: f=(4-3ψ)/ψ) **rescuje teorię na poziomie PPN**
(β_PPN=γ_PPN=1) i ma naturalne wyprowadzenie algebraiczne z
samego potencjału. Status M9.1'' = **otwarta hipoteza pivotu B**
do testowania (P1–P4).

### 9.7 Pliki M9.1''

| Plik | Rola |
|------|------|
| `M9_1_pp_setup.md` | Setup analityczny M9.1'' (postulat, derivation, test plan) |
| `m9_1_pp_verify.py` | Numeryczna weryfikacja β_PPN=0.996 z istniejącymi danymi M9.1 |
| `m9_1_pp_verify.txt` | Wyjście weryfikacji |
| `M9_1_prime_results.md` | Ten dokument (§9 = breakthrough M9.1'') |
