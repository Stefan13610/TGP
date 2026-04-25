# M9.1'' — wyprowadzenie metryki z potencjału substratu

**Cykl:** M9.1'' (rozszerzenie M9.1' Pivot B).
**Data:** 2026-04-25.
**Status:** propozycja teoretyczna z weryfikacją numeryczną.

Zob. `M9_1_prime_results.md` §4.2 dla kontekstu Pivotu B.
Zob. `m9_1_pp_verify.py` / `.txt` dla weryfikacji numerycznej.

---

## 1. Hipoteza centralna

**Postulat M9.1'' (`ax:metric-from-potential`, propozycja):**

Czynnik metryczny `f(ψ)` w `g_tt = -c²·f(ψ)` jest proporcjonalny do
**znormalizowanej gęstości potencjału substratu**:

```
                     V(Φ)              (12 Φ₀²)     (4 - 3ψ)
   f(ψ)  =  K · ─────────  =  ───────── · V(Φ)/Φ⁴  =  ─────────
                    Φ⁴                  γ                 ψ
```

gdzie:
- `V(Φ) = (β/3)Φ³/Φ₀ - (γ/4)Φ⁴/Φ₀²` (potencjał TGP, sek08a),
- `β = γ` (warunek próżni, sek08a `prop:vacuum-condition`),
- `K = 12Φ₀²/γ` to stała normalizacyjna ustawiona tak by `f(1) = 1`.

`h(ψ) = 1/f(ψ) = ψ/(4-3ψ)` (z budżetu substratowego f·h=1, sek08c
`prop:antipodal-from-budget`).

## 2. Pełne wyprowadzenie algebraiczne

### 2.1 Potencjał TGP

Z sek08a (warunek próżni β=γ) potencjał substratowy:

```
   V(Φ) = (β/3) Φ³/Φ₀ - (γ/4) Φ⁴/Φ₀²
        = (γ Φ₀²/12) · ψ³ · (4 - 3ψ)            [β=γ]
```

Potencjał ma:
- Zero potrójne przy `ψ = 0` (faza niemetryczna substratu).
- Minimum przy `ψ = 1` (próżnia: V(Φ₀) = -γΦ₀²/12 < 0,
  energia próżni).
- Punkt zerowy przy `ψ = 4/3` (drugi zero V, "brzeg basenu
  potencjałowego").
- V(Φ) → -∞ dla ψ → ∞ (artefakt obcięcia, zob. KNOWN_ISSUES C2).

### 2.2 Gęstość potencjału V(Φ)/Φ⁴

W kanonicznych zmiennych (`sek08_formalizm` `prop:substrate-action`)
naturalna miara to `Φ⁴` (sprzężenie konformalne K = K_geo·ψ⁴
generujące α=2). Stąd dimensionless density:

```
          V(Φ)         γ·Φ₀²·ψ³·(4-3ψ)/12             γ      4 - 3ψ
   F(ψ) := ─── =  ───────────────────────  =  ────────── · ──────────
          Φ⁴            Φ₀⁴ · ψ⁴               12·Φ₀²        ψ
```

przy ψ=1: `F(1) = γ/(12Φ₀²)` — stała próżni.

Normalizacja `f(ψ) = F(ψ)/F(1)`:

```
                         4 - 3ψ
   f(ψ)  =  ──────────  =  ────────  =  -3 + 4/ψ
            ψ                  ψ
```

### 2.3 Pochodne f w próżni

```
   f(1)    =  1                   ✓ normalizacja próżni
   f'(ψ)   =  -4/ψ²
   f'(1)   =  -4
   f''(ψ)  =  +8/ψ³
   f''(1)  =  +8
```

### 2.4 Sprawdzenie γ_PPN i β_PPN

**γ_PPN** (z `f·h=1`, automatyczne dla każdego f):
```
   γ_PPN = -h'(1)/f'(1) = -(-f'(1))/f'(1) = +1                ✓
```

**β_PPN** (master formula z M9.1' §3.1, c₂=-α/2=-1 dla α=2):
```
   β_PPN = f''(1)/f'(1)² + 2·c₂/f'(1)
         = 8/16 + 2·(-1)/(-4)
         = 0.5 + 0.5
         = 1.000     EXACTLY                                  ✓
```

### 2.5 Aksjomat prędkości ax:c — zmodyfikowany

W obecnym sek08c thm:antipodal-uniqueness:
```
   c² = c₀² · f/h = c₀²/ψ
```

W M9.1'' z f=(4-3ψ)/ψ, h=ψ/(4-3ψ):
```
   c² = c₀² · f/h = c₀² · (4-3ψ)²/ψ²
```

To **inna funkcjonalna postać** prędkości światła. Punkty kluczowe:
- ψ = 1: c = c₀ ✓ (kalibracja próżni).
- ψ = 4/3: c → 0 (prędkość światła zanika na brzegu basenu
  potencjałowego — fizyczna interpretacja: "sfera zatrzymanego
  światła" przy granicy substratu).
- ψ → 0: c → ∞·1/0 = niewykluczona patologia (zgodne z fazą
  niemetryczną).

Aksjomat ax:c wymaga przepisania: zamiast „pojedynczy power
form ψ^k" → „c² ∝ V(Φ)²/Φ⁸ × const", lub równoważnie c² ∝ f²(ψ).

## 3. Interpretacja fizyczna

### 3.1 Zasada „metryka z potencjału"

Hipoteza centralna mówi:

> **Lokalna miara czasu (g_tt) jest proporcjonalna do gęstości
> potencjału substratu V(Φ)/Φ⁴.**

Gdzie:
- V(Φ)/Φ⁴ to "potencjał na jednostkę kwadratu fazowej miary
  konformalnej" (skoro K(ψ)=ψ⁴, to ψ⁴ to naturalna miara objętości
  substratu w mocy 4).
- Zegar lokalny tyka tym wolniej, im niższa gęstość energii
  potencjalnej (nie odwrotnie — uwaga na znak).

### 3.2 Trzy granice fizyczne

| ψ | V(Φ) | V(Φ)/Φ⁴ | f(ψ) | g_tt | Interpretacja |
|---|------|---------|------|------|---------------|
| 0 | 0 (z V~Φ³) | ∞ (z 1/Φ) | ∞ | -∞ | faza niemetryczna substratu |
| 1 | min próżnia | γ/(12Φ₀²) | 1 | -c² | kalibracja próżni |
| 4/3 | 0 (drugi zero V) | 0 | 0 | 0 | brzeg basenu potencjału |

Trzy fizyczne progi:
- **ψ=0**: substrat całkowicie ekstynktowany; metryka się rozpada.
- **ψ=1**: stan próżni minimum potencjału; kalibracja c=c₀.
- **ψ=4/3**: punkt równowagi cubic-quartic V→0; horyzont
  współrzędnych (g_tt → 0, h → ∞).

### 3.3 Związek z basinem ghost-free

`prop:ghost-free-fundamental` (sek08_formalizm):
> Brak niestabilności gradientowych w basenie ψ ∈ (0, 4/3).

Forma hiperboliczna f(ψ) = (4-3ψ)/ψ jest niezerowo określona
DOKŁADNIE w tym samym basenie:
- f > 0 dla ψ ∈ (0, 4/3) ✓
- f = 0 przy ψ = 4/3 (brzeg basenu)
- f → ∞ przy ψ → 0 (drugi brzeg)

To **wewnętrzna spójność**: metryka jest dobrze określona dokładnie
tam, gdzie substrat jest ghost-free.

## 4. Status logiczny: czy to jest derivacja, czy curve fit?

### 4.1 Argumenty za derivacją

1. **Forma f(ψ) jest jednoznacznie wyznaczona** przez trzy zasady:
   - (a) `prop:antipodal-from-budget` (sek08c): f·h=1.
   - (b) γ_PPN=1 → automatyczne z (a).
   - (c) **„f znika w drugim zera potencjału V"** (NOWY postulat).
   - (d) Normalizacja f(1)=1 (kalibracja próżni).
   
   Zasady (a) i (d) są już w sek08c. Zasada (c) jest nowa, ale
   sharp i fizycznie motywowana (V→0 → metryka osobliwa).

2. **β_PPN=1 wynika** automatycznie z (a)-(d), nie jest osobnym
   postulatem.

3. **(4-3ψ) w f to TEN SAM czynnik** co w V(Φ); współczynniki β,γ
   wymuszają to z warunku próżni β=γ.

### 4.2 Argumenty za curve fit

1. Zasada (c) („f znika w drugim zera V") jest **nowa** — nie
   wynika z obecnego sek08c automatycznie. Trzeba ją postulować.

2. Inne dziwniejsze zasady mogłyby też dawać β_PPN=1 z różnymi
   formami f (np. forma wymierna z biegunem w ψ=3/2, kandydat (3)
   z `m9_1_prime_pivot_B.py`).

3. „Naturalność" wyboru V/Φ⁴ jako gęstości jest dyskusyjna —
   równie naturalna mogłaby być V/Φ², V/Φ³, czy inna miara.

### 4.3 Konkluzja statusu

Forma hiperboliczna jest **propozycją teoretyczną** z silnymi
wskazówkami ale bez ostatecznego wyprowadzenia. Trzy testy mogą
ją wzmocnić lub obalić:

1. **Test wariacyjny:** czy istnieje akcja z której wariacja po
   metryce daje `g_tt ∝ V(Φ)/Φ⁴`? (Hipoteza: tak, jeśli mnożymy
   action TGP przez `1/Φ⁴` i traktujemy f jako pomocniczy
   stopień swobody. Do sprawdzenia.)

2. **Test obserwacyjny dodatkowy:** czy forma hiperboliczna
   przewiduje DODATKOWE odchylenia od GR przy wyższych rzędach
   PN, które mogłyby być testowane (Cassini, LLR, binary pulsars)?
   Master formula daje O(1/r³), O(1/r⁴) inne niż GR — można
   policzyć.

3. **Test kosmologiczny:** zachowanie f w skalach kosmologicznych
   (G_eff(z), Λ_eff). Forma hiperboliczna może dawać różne
   przewidywania niż boxed forma.

## 5. Implikacje dla rdzenia TGP

### 5.1 Zmiany w sek08c

Wymagana **głęboka rewizja sek08c**:

1. **Boxed eq:metric-full-derived** (linia 154): zastąpić formą
   hiperboliczną:
```
   ds² = -c₀²·(4-3ψ)/ψ · dt² + ψ/(4-3ψ) · δ_ij · dx^i dx^j
```

2. **`thm:antipodal-uniqueness`** (linia 345): jest nieprawdziwe
   w obecnej postaci (s=1 nie wyznacza JEDNEJ formy bez
   dodatkowych warunków). Zmodyfikować na:
```
   Theorem (M9.1''): pod warunkami f·h=1, γ_PPN=1,
   f znika przy drugim zera V(Φ) (ψ=4/3), oraz f(1)=1,
   forma metryczna jest JEDNOZNACZNIE f(ψ) = (4-3ψ)/ψ.
```

3. **Element objętościowy** (sek08c §sssec:volume-element-check):
   przeliczyć dla nowej f. `√-g = c₀·f^(1/2)·h^(3/2) · dt·d³x =
   c₀·(...)`. Stała sprzężenia κ wymaga przeliczenia. Test LLR
   może się zmienić.

### 5.2 Zmiany w sek_stale (ax:c)

**`ax:c`** (aksjomat prędkości światła) wymaga uogólnienia:
- Zamiast: `c² = c₀²/ψ` (forma potęgowa).
- Na: `c² = c₀²·(4-3ψ)²/ψ²` (forma kwadratowa hiperboliczna).

Implikacje dla aksjomatów A6, A8 (skalowanie c, G, ℏ z Φ):
trzeba przeliczyć łańcuch z thm:metric-from-substrate-full §3.

### 5.3 Zmiany w sek08_formalizm

`sek08_formalizm` operuje na DZIAŁANIU TGP, nie na metryce.
Φ-EOM jest niezmienna (zależy od kinetyki K(ψ)=ψ⁴ i potencjału V,
nie od f). Zmiana to: re-interpretacja jak czytamy g_eff z Φ.

Pozostałe twierdzenia (D-uniqueness α=2, ghost-freeness) niezmienne.

### 5.4 Newton matching i stała sprzężenia

Z M9.1: `a₁ = q·M/(4π)` i `Newton: ε = 2GM/(c²r)`.
Z f'(1) = -4 (zamiast -1): teraz przy `g_tt = -c²·f(1+ε)`:
```
   g_tt ≈ -c²(1 + f'(1)·ε) = -c²(1 - 4ε)
```
Porównanie z PPN `g_tt = -c²(1 - 2U/c²)`:
```
   2U/c² = 4ε  ⟹  ε = U/(2c²)
```
Newton: U = GM/r, więc ε(r) = GM/(2c²r) = a₁/r → a₁ = GM/(2c²).
Z M9.1 mamy a₁ = qM/(4π), więc:
```
   q M/(4π) = GM/(2c²)  ⟹  q = 2πG/c²
```
Czyli q jest **4× mniejsze** niż w boxed metryce (gdzie
`q = 8πG/c²`). To jest tylko renormalizacja stałej sprzężenia
materii, nie wpływa na struktur teorii.

## 6. Plan testów M9.1''

### 6.1 P1 — analityczny

- [x] Wyprowadzić formę f(ψ) = V(Φ)/Φ⁴ symbolicznie. (✓ §2.2)
- [x] Sprawdzić β_PPN, γ_PPN dokładnie. (✓ §2.4)
- [ ] Sprawdzić wyższe rzędy PN: O(c₃, c₄, c₅) z formy hiperbolicznej
      vs. GR. Czy zgodne, czy są dodatkowe odchylenia?

### 6.2 P2 — wariacyjny

- [ ] Czy istnieje akcja S[Φ, g] z której g_tt ∝ V(Φ)/Φ⁴ wynika
      przez δS/δg^tt = 0?
- [ ] Czy ten S jest zgodny z istniejącym S_TGP (sek08_formalizm)?

### 6.3 P3 — obserwacyjny

- [ ] LLR (lunar laser ranging): zmieniona stała κ wymaga sprawdzenia
      ̇G/G < 2·10⁻² constraints. (zob. `rem:metric-action-consistency`
      w sek08c).
- [ ] Cassini γ_PPN-1 < 2.3·10⁻⁵: spełnione automatycznie (γ=1 z f·h=1).
- [ ] Mercury β_PPN-1 < 1·10⁻⁴: spełnione (β=1 dokładnie).
- [ ] Binary pulsars (PSR B1913+16): dynamiczne testy PPN. Wymaga
      uogólnienia na poruszające się masy.

### 6.4 P4 — pełne paper

- [ ] Przepisać sek08c z nową boxed formą.
- [ ] Przepisać sek_stale z nową ax:c.
- [ ] Sprawdzić spójność z sek09 (kosmologia).
- [ ] Spójność z rozdziałami eksperymentalnymi (LIGO, EHT, PSR).

## 7. Status

**M9.1'' jest OTWARTĄ propozycją teoretyczną.** β_PPN=1 jest
spełnione dokładnie (analitycznie i numerycznie z dokładnością
0.4% przy R_max=800). Status logiczny: "warunkowy pivot zgodny
z danymi" — ratuje OP-2b, ale wymaga (i) wariacyjnego wyprowadzenia
postulatu (c), (ii) rewizji sek08c, sek_stale, (iii) sprawdzenia
wyższych rzędów PN, (iv) testów obserwacyjnych dodatkowych.

Jeśli M9.1'' przejdzie test wariacyjny (P2), TGP zostaje
**uratowany na poziomie post-newtonowskim**. Jeśli nie — falsyfikacja
M9.1' zostaje.

## 8. Pliki

| Plik | Rola |
|------|------|
| `M9_1_pp_setup.md` | Ten dokument — analityczne wyprowadzenie |
| `m9_1_pp_verify.py` | Numeryczna weryfikacja β_PPN=1 z formy V/Φ⁴ |
| `m9_1_pp_verify.txt` | Wyjście numeryczne |
| `m9_1_prime_pivot_B.py` | Skan kandydatów z którego (4) hyperbolic wynikła |
| `M9_1_prime_results.md` | Werdykt M9.1' (kontekst) |
