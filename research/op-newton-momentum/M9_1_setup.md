# M9.1 — Statyka i PPN: setup analityczny

**Cykl:** M9 ("klasyczna dynamika i pęd"), test 1 z 3.
**Data:** 2026-04-25.
**Cel:** rozwiązać Φ-EOM TGP dla statycznego sferycznie symetrycznego
źródła w słabym polu, zweryfikować odzyskanie Newtona, dopasować stałą
sprzężenia `q ↔ G`, przeczytać γ_PPN i β_PPN.

Zob. `M9_program.md` §3 dla nadrzędnego planu, `TGP_FOUNDATIONS.md`
§3 dla akcji i Φ-EOM, §5 dla interpretacji grawitacji.

---

## 1. Punkt wyjścia: Φ-EOM

Z `core/sek08a_akcja_zunifikowana.tex` `eq:field-eq-reproduced`:

```
∇²Φ + 2 (∇Φ)²/Φ + β Φ²/Φ_0 - γ Φ³/Φ_0² = -q Φ_0 ρ                    (1)
```

W zmiennej bezwymiarowej `φ = Φ/Φ_0`, dzieląc obie strony przez Φ_0:

```
∇²φ + 2 (∇φ)²/φ + β φ² - γ φ³ = -q ρ                                 (2)
```

z warunkiem próżni `β = γ` (`prop:vacuum-condition`).

## 2. Linearyzacja wokół próżni `φ = 1 + ε`

Zakładamy `|ε| ≪ 1` i rozwijamy.

### 2.1 Człony nieliniowe

- `(∇φ)² / φ = (∇ε)² / (1+ε) = (∇ε)² (1 - ε + ε² - ...)`
  → przy linearyzacji: `O(ε²)` (drops at leading order in ε).
- `φ² = (1+ε)² = 1 + 2ε + ε²`
- `φ³ = (1+ε)³ = 1 + 3ε + 3ε² + ε³`

### 2.2 Rozwinięcie po ε

**O(1):**
```
β · 1 - γ · 1 = β - γ = 0   ⟸   warunek próżni, automatycznie spełniony.
```

**O(ε):**
```
∇²ε + 2β · ε - 3γ · ε = -q · ρ                                       (3a)
```
Z `β = γ`:
```
∇²ε + (2β - 3β) ε = -q ρ
∇²ε - β · ε = -q · ρ                                                 (3)
```

Jest to **równanie Helmholtza/Yukawy** ze "skalą masową" `m² = β`
(zakładamy β > 0).

**O(ε²):**
```
∇²ε² + 2(∇ε)² · 1 + β(2ε² + ε²) - γ(3ε² + 3ε²·ε) + ...
```
Pełna nieliniowa korekta wymaga rozwiązania iteracyjnego (Krok 4
poniżej).

## 3. Rozwiązanie statyczne sferycznie symetryczne

W sferycznych współrzędnych, `ε = ε(r)`:

```
(1/r²) (r² ε')' - β ε = -q ρ(r)                                      (4)
```

### 3.1 Limit β → 0 (massless, czyste Newton)

```
(1/r²) (r² ε')' = -q ρ(r),    czyli   ∇²ε = -q ρ.
```

Dla punktowego źródła `ρ = M δ³(x)`, klasyczne rozwiązanie z warunkiem
brzegowym `ε(∞) = 0`:

```
ε(r) = q M / (4π r)                                                  (5)
```

To jest **dokładna postać Newtona** w zmiennej ε. Identyfikacja z
potencjałem grawitacyjnym U:

```
ε = 2 U / c_0²    (z metryki: g_tt = -c_0²/φ ≈ -c_0² (1 - ε))
U = G M / r       (Newton)
⟹   q M / (4π r) = 2 G M / (c_0² r)
⟹   q = 8π G / c_0²                                                  (6)
```

To jest **identyfikacja efektywnej stałej sprzężenia TGP ↔ Newton**.
Jeden parametr (q) determinuje skalę grawitacyjną w TGP, dopasowywany
do G obserwacyjnego.

### 3.2 Skończone β > 0 — promień Yukawy

Pełne rozwiązanie (4) z `m² = β`:

```
ε(r) = (q M / 4π r) · exp(-√β · r)                                  (7)
```

**Promień Yukawy** `R_Y = 1/√β`. Dla zachowania Newtona w obserwacjach
wymagamy `R_Y` >> skala obserwacyjna:

| Skala obserwacji | Wymóg na R_Y | Wymóg na √β |
|---|---|---|
| System słoneczny (~AU = 1.5×10¹¹ m) | R_Y ≫ AU | √β ≪ 6.7×10⁻¹² m⁻¹ |
| Galaktyki (~kpc = 3×10¹⁹ m) | R_Y ≫ kpc | √β ≪ 3×10⁻²⁰ m⁻¹ |
| Skala kosmologiczna (H_0/c ≈ 4×10⁻²⁷ m⁻¹) | R_Y ≫ c/H_0 | √β ≪ 4×10⁻²⁷ m⁻¹ |

**Wniosek:** TGP albo:
- (i) `β = 0` ścisle (skala masowa φ-fluktuacji = 0 w fazie złamanej),
  co daje czysty Newton w każdej skali; albo
- (ii) `β > 0` ale fenomenologicznie ograniczone przez horyzont:
  `√β < H_0/c ≈ 4×10⁻²⁷ m⁻¹`. To jest też potencjalna kosmologiczna
  signature TGP.

W M3–M8 (poprzedni cykl) liczby konkretne dla β nie były wyznaczone.
Tutaj pozostawiamy β jako parametr fenomenologiczny.

## 4. Linearyzacja metryki i PPN

Z `eq:g-eff-unified`:
```
g_eff_tt = -c_0² / φ = -c_0² / (1 + ε)
g_eff_rr = +φ = 1 + ε
```

Rozwijamy w ε:

```
g_eff_tt = -c_0² · (1 - ε + ε² - ε³ + ...)
g_eff_rr = +1 · (1 + ε)
```

Identyfikacja `ε = 2U/c_0²` (z M9.1 §3.1, eq. (6)):

```
g_eff_tt = -c_0² (1 - 2U/c_0² + 4U²/c_0⁴ - ...)                     (8a)
g_eff_rr = 1 + 2U/c_0²                                              (8b)
```

### 4.1 Standard PPN

Konwencja PPN (Will, "Theory and Experiment in Gravitational Physics"):

```
g_tt = -c² · (1 - 2U/c² + 2 β_PPN U²/c⁴ + ...)
g_rr = 1 + 2 γ_PPN U/c²
```

### 4.2 Odczytanie γ_PPN, β_PPN z (8a, 8b)

Z (8b): `2 γ_PPN U/c_0² = 2 U/c_0²` ⟹ **`γ_PPN = 1` dokładnie**.

Z (8a): `2 β_PPN U²/c_0⁴ = 4 U²/c_0⁴` ⟹ **`β_PPN = 2`** w formie potęgowej.

### 4.3 Forma eksponencjalna (sek08c lin. 192–211)

Alternatywna parametryzacja: `g_tt = -c_0² e^(-2U/c_0²)`,
`g_rr = e^(+2U/c_0²)`. Spełnia ten sam warunek `f · h = 1` (budżet
substratowy). Daje:

```
e^(-2U/c²) = 1 - 2U/c² + 2U²/c⁴ - ...   ⟹   β_PPN = 1.
```

Wniosek `sek08c`: substrat nie rozróżnia między formą potęgową a
eksponencjalną na poziomie fundamentalnym; warunek obserwacyjny
`β_PPN = 1` wybiera formę eksponencjalną.

## 5. Pytanie kluczowe dla M9.1

**Czy poprawki O(ε²) z Φ-EOM (3a) "dynamicznie" preferują formę
potęgową czy eksponencjalną?**

Możliwe odpowiedzi:
- **(A) Forma potęgowa wyłaniają się z liniowej linearyzacji.** Wtedy
  TGP daje `β_PPN = 2` na poziomie linearyzacji, co jest **w napięciu z
  obserwacjami** (Mercury perihelion mierzy efektywnie kombinację
  PPN parameters do precyzji ~10⁻⁴, β_PPN = 2 byłoby wykluczone).
- **(B) Forma eksponencjalna wyłaniają się przy uwzględnieniu poprawek
  O(ε²) w Φ-EOM.** Wtedy `β_PPN = 1` jest **wynikową predykcją**, nie
  arbitralnym wyborem reparametryzacji.
- **(C) Φ-EOM nie determinuje wyboru** — obie formy są równoważnie
  rozwiązaniami do dowolnego rzędu. Wtedy wybór jest fenomenologiczny i
  nie jest predykcją TGP.

**M9.1 musi rozróżnić te trzy przypadki.**

Strategia obliczeniowa:

1. Rozwiązać Φ-EOM (2) **iteracyjnie** do O(ε²):
   - ε⁽¹⁾(r) = rozwiązanie linearyzowane (eq. 3, 5).
   - ε⁽²⁾(r) = rozwiązanie z uwzględnieniem (∇Φ)²/Φ jako źródła i
     poprawek nieliniowych potencjału.
2. Sprawdzić, czy ε⁽²⁾ ma postać:
   - Potęgową: ε⁽²⁾(r) ≈ ε⁽¹⁾(r) - (1/2)[ε⁽¹⁾(r)]² + ... (rozwinięcie 1/(1+ε));
   - Eksponencjalną: ε⁽²⁾(r) ≈ ε⁽¹⁾(r) - [ε⁽¹⁾(r)]² + ... (rozwinięcie e^ε).
3. Z porównania współczynników wyciągnąć **predykcję β_PPN** z TGP.

## 6. Implementacja `m9_1_static.py`

**Strategia:** numerycznie rozwiązać ODE (4) z poprawkami nieliniowymi:

```python
def phi_eom_residual(phi, r, beta=0.0, q=1.0, rho=...):
    # ODE: (1/r^2)(r^2 phi')' + 2(phi')^2/phi + beta phi^2 - gamma phi^3 = -q rho
    # gamma = beta (warunek próżni)
    ...

def solve_static(beta, q, rho_profile, r_grid, ...):
    # ODE solver, condition phi(infty) = 1
    # Output: phi(r), eps(r) = phi - 1, identify beta_PPN, gamma_PPN
    return phi_array, gamma_PPN, beta_PPN
```

Zakres testów:
- (T1) Limit β → 0, ρ = punktowe: weryfikacja `ε = qM/(4πr)` (Newton).
- (T2) Linearyzacja (drop O(ε²)): odczytanie γ_PPN, β_PPN, porównanie z
  formą potęgową.
- (T3) Pełne nieliniowe rozwiązanie: forma `ε(r)` przy małym ale
  niezerowym `qM`. Dopasowanie do wzorów potęgowego i eksponencjalnego.
- (T4) Skończone β > 0: weryfikacja obcięcia Yukawy.
- (T5) Rozszerzone źródło `ρ(r)` (gładkie, np. gauss): brak rozbieżności
  w r → 0.

## 7. Output

Plik `m9_1_results.txt` powinien zawierać:

```
=== M9.1 — Statyka i PPN ===

T1: Limit β → 0, punktowe źródło
    ε_num(r) na siatce ⟶ porównanie z qM/(4πr)
    Błąd względny: ...

T2: Linearyzacja, odczyt PPN parameters
    γ_PPN = 1.0000   (oczekiwane: 1.0000)
    β_PPN = 2.0000   (oczekiwane: 2.0 w formie potęgowej)

T3: Pełne nieliniowe rozwiązanie do O(ε²)
    Dopasowanie do formy potęgowej: residualy = ...
    Dopasowanie do formy eksponencjalnej: residualy = ...
    Werdykt: forma potęgowa | eksponencjalna | nierozróżnialne

T4: Skończone β > 0
    R_Yukawa = 1/√β = ...
    Wytłumienie eksponencjalne potwierdzone: ...

T5: Rozszerzone źródło
    ρ(r) = ρ_0 exp(-r²/2σ²)
    ε(r) gładkie, brak osobliwości w r=0: ✓
```

## 8. Plik wynikowy `M9_1_results.md`

Werdykt na koniec:
- ✓ Newton odzyskany w limicie β → 0.
- ✓ γ_PPN = 1 dokładnie.
- Werdykt o β_PPN: 1 (eksponencjalna), 2 (potęgowa), albo
  nierozróżnialne — z odpowiedzią Φ-EOM.
- Skala Yukawy `√β` jako parametr fenomenologiczny do późniejszego
  wyznaczenia.
