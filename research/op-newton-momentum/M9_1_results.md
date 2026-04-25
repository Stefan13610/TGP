# M9.1 — Statyka i PPN: wyniki

**Cykl:** M9 ("klasyczna dynamika i pęd"), test 1 z 3.
**Data:** 2026-04-25.
**Solver:** `m9_1_static.py` (Φ-EOM BVP, scipy.integrate.solve_bvp).
**Surowe wyjście:** `m9_1_results.txt`.

Zob. `M9_1_setup.md` dla setupu analitycznego, `M9_program.md` §3 dla
nadrzędnego planu.

---

## 1. Streszczenie

Numerycznie rozwiązano statyczną sferycznie symetryczną Φ-EOM TGP:

```
∇²ε + 2(∇ε)²/(1+ε) - β(1+ε)²·ε = -q·ρ(r)        (φ = 1+ε, β=γ próżni)
```

w zmiennej `v(r) = r·ε(r)`, z BC `v(0)=0` (regularność) i
`v'(R_max)=0` (asymptotyczna płaskość). Wykonano testy T1–T5.

**Werdykt 5 testów:**

| Test | Wynik | Krótko |
|------|-------|--------|
| T1 | ✓ PASS | Newton odzyskany w limicie β→0; max błąd względny `3.4×10⁻³` |
| T2 | ✓ PASS | Linearyzacja: γ_PPN = 1.000 (dokładnie), β_PPN = 2.000 (forma potęgowa) |
| **T3** | **✗ FAIL** | **Pełne nieliniowe Φ-EOM: β_PPN = 3.984 ≈ 4 — wykluczone przez Mercury** |
| T4 | ✓ PASS | Yukawa potwierdzona dla β=0.01 i β=0.1 (rel. err `~3×10⁻³`) |
| T5 | ✓ PASS | Brak rozbieżności w r=0 dla rozszerzonego źródła Gaussa |

**Globalny werdykt M9.1:** ✗ — TGP w obecnej formie (jednopolowa Φ-EOM
z metryką `g_tt = -c²/φ`, `g_rr = φ`) **nie odzyskuje GR na poziomie
post-newtonowskim**. Próba pivotu w M9.1′ niżej.

## 2. T3 — analiza ilościowa (KEY RESULT)

### 2.1 Strategia ekstrakcji β_PPN

Dopasowano `ε(r)` w obszarze `r ∈ [8σ, 80σ]` do szeregu asymptotycznego:

```
ε(r) = a₁/r + a₂/r² + a₃/r³ + a₄/r⁴ + a₅/r⁵
```

`a₁ = A_eff` to renormalizowana stała sprzężenia (zawiera poprawkę
energii własnej `(1/2π)·∫(∇ε)²/(1+ε)·d³r` powyżej naiwnego `A₀=qM/(4π)`).

Definiujemy bezwymiarowy współczynnik: **c₂ := a₂/A_eff²**.

Mapowanie do PPN: `β_PPN = 2(1 - c₂)` (wyprowadzenie w `M9_1_setup.md` §4).

### 2.2 Analityczna predykcja TGP

Iteracyjne rozwiązanie EOM dla `r >> σ` daje (samokonsystentnie):

```
ε(r) = A/r - A²/r² + (5/3)A³/r³ - (10/3)A⁴/r⁴ + O(1/r⁵)
```

zatem:
- c₂ = -1 (dokładnie)
- c₃ = +5/3
- c₄ = -10/3
- → **β_PPN = 2(1 - (-1)) = 4** (analitycznie)

### 2.3 Studium zbieżności (R_max)

Warunek brzegowy `v'(R_max)=0` wprowadza obcięcie ogona `1/r²` na
brzegu siatki, biasując `a₂` o `O(A²/R_max)`. Wymagana zbieżność:

| R_max | n_pts | c₂ |
|-------|-------|-----|
| 100 | 1500 | -0.709 |
| 200 | 2000 | -0.869 |
| 400 | 3000 | -0.967 |
| **800** | **5000** | **-0.992** |

→ R_max = 800 wystarczy do dokładności ~1% na c₂. (Detale: `debug_rmax.py`.)

Sanity check: ten sam pipeline na **linearyzowanym solverze** (gdzie
analitycznie a₂ = 0) daje |a₂| ~ 10⁻¹¹, czyli ~12 rzędów wielkości
poniżej sygnału w pełnym rozwiązaniu — szum numeryczny nie tłumaczy
otrzymanego c₂ = -0.99.

### 2.4 Ostateczny pomiar (R_max=800, sweep po q)

| q | A₀ = qM/(4π) | A_eff (a₁) | a₂ | c₂ = a₂/A_eff² | β_PPN |
|---|--------------|-----------|-----|----------------|-------|
| 0.500 | 0.03979 | 0.041604 | -1.717·10⁻³ | -0.9918 | +3.9836 |
| 1.000 | 0.07958 | 0.087173 | -7.537·10⁻³ | -0.9918 | +3.9836 |
| 2.000 | 0.15915 | 0.192010 | -3.657·10⁻² | -0.9918 | +3.9837 |
| 3.000 | 0.23873 | 0.318587 | -1.007·10⁻¹ | -0.9919 | +3.9838 |

Ekstrapolacja liniowa do A₀ → 0 (granica słabego pola):

```
c₂_∞   = -0.9918
slope  = -0.0005   (c₂ praktycznie A-niezależne, jak powinno być)
β_PPN(A₀ → 0) = +3.984
```

Różnica do wartości analitycznej -1 (czyli ~1%) wynika ze skończonego
R_max=800; konwergencja do -1 obserwowana przy R_max → ∞ (§2.3).

### 2.5 Werdykt T3

| Forma | c₂ | β_PPN | Status |
|-------|-----|-------|--------|
| Potęgowa (linear ε=u, brak dyn. korekt) | 0 | 2 | wykluczona obs. |
| **TGP dynamiczna (ε = u - u² + ...)** | **-1** | **4** | **TGP, wykluczona** |
| Eksponencjalna (ε = e^u - 1) | +1/2 | 1 | GR-zgodna |

**TGP z pełną nieliniową Φ-EOM PRZEWIDUJE β_PPN = 4.**

Obserwacje (Will 2014, Cassini & Mercury & lunar laser ranging):
β_PPN = 1.000 ± 1×10⁻⁴.

Odchylenie: `(4 - 1)/10⁻⁴ = 3×10⁴`σ. **Mocne wykluczenie.**

## 3. Implikacje dla TGP

M9.1 ma trzy możliwe odczyty (zob. `M9_1_setup.md` §5):

- (A) Forma potęgowa surowa → β_PPN = 2 (wykluczona). Stało się.
- (B) Forma eksponencjalna emergentna → β_PPN = 1 (GR-zgodne). NIE wystąpiło.
- (C) Φ-EOM nie determinuje wyboru. NIE wystąpiło — Φ-EOM determinuje, ale w stronę c₂ = -1, β_PPN = 4.

**Faktyczny wynik to wariant (A′)**, jeszcze gorszy niż (A): nieliniowe
poprawki Φ-EOM dają jeszcze WIĘKSZĄ wartość β_PPN niż sama forma
potęgowa metryki bez poprawek dynamicznych.

### 3.1 Możliwe interpretacje

1. **TGP w obecnym sformułowaniu jest fałszyfikowane przez M9.1.**
   Mercury, Cassini i LLR zamykają OP-2b negatywnie.

2. **Brakuje mechanizmu po stronie materii.** TGP w obecnej formie
   uwzględnia materię tylko przez axiom `ax:metric-coupling` (materia
   sprzęga się przez `g_eff`). Może back-reaction materii (np.
   stress-energy modyfikujący Φ-EOM przez wyższe rzędy) zmieni
   nieliniową korektę. M9.2 (pęd/inercja) może to oświetlić.

3. **Substrat wymaga rozszerzenia.** Ale kluczowa zasada
   (`TGP_FOUNDATIONS.md` §1.2): substrat = jedno pole skalarne Z₂,
   nieprzenośna. Pivot poza tym = nowa teoria, nie TGP.

4. **Metryka wymaga rewizji ansatzu.** TGP postuluje konkretną
   strukturę `g_tt = -c²/φ`, `g_rr = φ`. Inny ansatz (np.
   eksponencjalny, lub wynikający z geometrycznego wyprowadzenia
   substratu) może dać c₂ = +1/2 (forma eksponencjalna). Jednak
   **`sek08c` wyprowadził obecny ansatz z minimalizacji budżetu
   substratowego** — zmiana ansatzu wymaga rewizji `sek08c`.

### 3.2 Plan natychmiastowy (M9.1 → M9.1′)

Przed kolejnymi kandydatami w M9.2/M9.3 należy rozstrzygnąć:

**M9.1′ — Test alternatywnych ansatzów metryki**:
- (a) Sprawdzić, czy ansatz eksponencjalny `g_tt = -c²·exp(-η)`,
  `g_rr = exp(+η)` (z `η = ε` lub innym mapowaniem) jest
  konsystentny z Φ-EOM po linearyzacji do O(ε²). Jeśli `η = ε`,
  to wymagałoby `1+ε = exp(η) ⟹ ε = exp(η) - 1` — czyli ε ma być
  funkcją wykładniczą jakiegoś bardziej fundamentalnego potencjału η,
  a sama Φ-EOM dla φ wciąż daje to co dała (c₂(ε) = -1). Pytanie:
  jak jest η powiązane z φ?
- (b) Alternatywnie: czy substrat dopuszcza nieliniową
  reparametryzację `Φ → Φ' = f(Φ)` która prowadziłaby do innej
  postaci metryki?

To jest aktywne otwarte pytanie do M9.1′. **Bez odpowiedzi na (a) lub
(b), M9.2 (pęd) i M9.3 (GW) tracą sens** — bez zgodnej z GR statyki,
nie ma fundamentu do pęd/GW.

## 4. Pliki

| Plik | Rola |
|------|------|
| `m9_1_static.py` | Solver Φ-EOM (BVP, scipy.solve_bvp) |
| `m9_1_results.txt` | Surowe numeryczne wyjście testów T1–T5 |
| `M9_1_setup.md` | Setup analityczny, wyprowadzenia |
| `M9_1_results.md` | Ten dokument — werdykt |
| `debug_rmax.py` | Studium zbieżności c₂ vs R_max (zob. §2.3) |
| `debug_t3.py` | Test pipeline'u dopasowania na linearyzowanym solverze |

## 5. Status M9 — propagacja do M9.2 i M9.3

- M9.1: ✗ (β_PPN = 4, wykluczone). Wymagany pivot M9.1′ przed dalszymi
  testami pęd/GW.
- M9.2 (pęd, Lenz-like back-reaction): wstrzymane do rozstrzygnięcia
  M9.1′.
- M9.3 (GW): wstrzymane do rozstrzygnięcia M9.1 + M9.2.

## 6. Konsekwencje dla OP-2b

Cykl M3–M8 zamknął warunek FP-uniwersalności (β/γ < 0 w Wilson-Fisher
3D). M9.1 zamyka klasyczną dynamikę post-newtonowską jako negatywną
**w obecnym sformułowaniu**. OP-2b wciąż otwarte:
- Drogą M9.1′ (alternatywne ansatze metryki).
- Drogą fizyki materii (back-reaction stress-energy modyfikująca
  Φ-EOM — niezbadane).
- Drogą rozszerzenia w obrębie Z₂ skalarnego, ale **bez gravitona**
  (`TGP_FOUNDATIONS.md` §5).

Wszystkie scenariusze poniżej mają wspólny wymóg: zachowują substrat
`Φ` jako jedyny stopień swobody, a grawitacja pozostaje kolektywnym
efektem fluktuacji (Sakharov/Verlinde/Volovik tradition, `TGP_FOUNDATIONS.md` §5).
