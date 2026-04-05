# Problem 3 ciał w warstwie `nbody`

> **Ostatnia synchronizacja dokumentacji:** 2026-04-05  
> **Pełny stan teorii (repo):** `ANALIZA_SPOJNOSCI_v5.md` w korzeniu `TGP_v1/` — niniejszy plik opisuje wyłącznie warstwę wielociałową i kod w `nbody/`.

**Nawigacja:** [README.md](README.md) · [ZALOZENIA_NBODY.md](ZALOZENIA_NBODY.md) · [ANALIZA_NBODY_INTEGRACJA.md](ANALIZA_NBODY_INTEGRACJA.md) · [examples/README.md](examples/README.md) · [archiwum sesji](_archiwum_docs/WYNIKI_SESJI_2026_03_21.md)

---

## 1. Status materiałów w folderze `/nbody`

### 1.1 Pliki wartościowe (core)

| Plik | Zawartość | Status |
|------|-----------|--------|
| `pairwise.py` | Dokładny potencjał 2-ciałowy V₂(d) | **EXACT** |
| `three_body_terms.py` | Siły 3-ciałowe z Φ⁴ | **EXACT** (Coulomb), APPROX (Yukawa) |
| `dynamics_v2.py` | Integrator leapfrog + RK45 | **NUMERYCZNY** |
| `equilibria.py` | Finder równowag statycznych | **EXACT** (pairwise), APPROX (3b) |
| `stability.py` | Hessian + projekcja modów zerowych | **RYGORYSTYCZNY** |
| `configurations.py` | Konfiguracje startowe | Pomocniczy |
| `tgp_nbody_special_configurations.tex` | Redukcja algebraiczna ogólna | **KLUCZOWY WYNIK** |
| `../dodatekD_trojcialowe.tex` (korzeń repo) | Wyprowadzenie sił 3-ciałowych w głównym tomie; w nagłówku pliku odniesienie do historycznej nazwy `trojcialowe_sily.tex` | **KLUCZOWY** |

### 1.2 I_triple — status (bez martwych odnośników)

**Wniosek:** Całka potrójnego Yukawa \(I_{\text{triple}}\) **nie ma prostej zamkniętej formy** dla ogólnej geometrii (odpowiednik diagramu trójkątnego w 3D QFT euklidesowej).

**Obecny tor w kodzie (2026-04):**

- `three_body_terms.py`: `triple_overlap_numerical()` — siatka 3D (wolniejsza, uniwersalna).
- `three_body_force_exact.py` — całka Feynmana 2D (dokładniejszy tor numeryczny dla sił).

Wcześniejsze skrypty pomocnicze typu `tgp_yukawa_*_fit_helper.py` oraz notatki `tgp_shape_space_*.tex` / `tgp_yukawa_overlap_*.tex` **nie występują już w drzewie repozytorium** (historyczne próby fitu / patchwork).

### 1.3 Starsze skrypty `run_*.py`

Katalog `nbody/` nie zawiera już syntetycznych `run_*.py` w korzeniu. Eksperymenty regresyjne są w `examples/` oraz `examples/_archiwum/` (m.in. potwierdzenie bariery \(1/d^2\), dynamiki przy \(d_{\rm well}\), skali \(|V_3/V_2|\sim C\)).

---

## 2. Niespójności znalezione i korekty

### 2.1 Historyczny błąd w dokumentacji V₃ (skalowanie 1/d² zamiast 1/d)

**Stary komentarz w kodzie (błędny, naprawiony):** sugerował skalowanie \(1/d^2\) dla \(V_3\) przy trójkącie równobocznym.

**Poprawne:**
W granicy Coulomba (m_sp → 0), dla trójkąta równobocznego d₁₂=d₁₃=d₂₃=d:

$$V_3 = -6\gamma C^3 \cdot I_{\text{triple}} = -6\gamma C^3 \cdot \frac{8\pi^2}{3d} = -\frac{16\pi^2 \gamma C^3}{d}$$

Skalowanie: **1/d** (nie 1/d²). Człon V₃ w granicy Coulomba jest tego samego rzędu co V_grad (atrakcyjny człon 1/d), ale wzmocniony czynnikiem γ·C.

### 2.2 KLUCZOWE ODKRYCIE: Reżim Coulomba vs Yukawa dla V₃

Formuła V₃ zależy od stosunku m_sp·d_well. Użycie **złej granicy** daje błąd rzędu 5000×:

| Parametry | m_sp = √β | d_well | m_sp·d_well | Reżim | V₃/V₂ |
|---|---|---|---|---|---|
| β=1.5, C=0.3 | 1.225 | 3.95 | **4.84** | **Yukawa** | 0.086% ✓ |
| β=1.5, C=0.3 | 1.225 | 3.95 | 4.84 | Coulomba (błąd!) | 455% ✗ |
| β→0, γ~10⁻⁵² | ~0 | dowolne | ~0 | Coulomba | ~γ·C ≪ 1 ✓ |

**Reguła:**
- Jednostki kodowe (β=γ=O(1)): m_sp·d_well ≫ 1 → **używać Yukawa V₃**
- Fizyczny TGP (γ~Λ_eff~10⁻⁵²): m_sp·d ≪ 1 → **używać Coulomba V₃**

W obu fizycznych reżimach |V₃/V₂| ≪ 1 → sektor 3-ciałowy jest **perturbacyjny**.

### 2.3 Spójność znaku sił 3-ciałowych

W `dynamics_v2.py` (linia 161–165):
```python
forces[i] += coupling * dI_dP * dP_dxi
```
gdzie `coupling = 6 * gamma * Ci * Cj * Ck` i dla Coulomba `dI_dP = -8π²/P² < 0`.
Wynik: siła przyciągająca (ujemna w kierunku od mas) → **FIZYCZNIE POPRAWNE** (V₃ < 0 = atrakcja).

### 2.3 Siła 3-ciałowa — dokładność statusów

| Metoda | Status | Dokładność |
|--------|--------|------------|
| Coulomb (m_sp=0): I = 8π²/P | **DOKŁADNY** | Exactact |
| Yukawa saddle-point: I ~ (2π)^(3/2)/m² · e^(-ms)/s | **PRZYBLIŻONY** | ~5–30% dla m·d ~ 1 |
| Numeryczny 3D grid N=80 | **NUMERYCZNY** | <5% dla dostatecznym N |

Dla fizycznych parametrów TGP: m_sp = √(3γ-2β) ≈ √γ ≈ √(10⁻⁵²) m⁻¹ — masa Yukawy jest OGROMNIE mała. Dla separacji makroskopowych m_sp·d ≪ 1 → **granica Coulomba** jest realistyczna.

---

## 3. Zaktualizowane równania (zgodność z `pairwise.py` / `three_body_terms.py`)

### 3.1 Potencjał 2-ciałowy (DOKŁADNY)

$$V_2(d) = -\frac{4\pi C_1 C_2}{d} + \frac{8\pi \beta C_1 C_2}{d^2} - \frac{12\pi \gamma C_1 C_2 (C_1+C_2)}{d^3}$$

Dla równych mas C₁=C₂=C (i γ=β, warunek próżni):

$$V_2(d) = -\frac{4\pi C^2}{d} + \frac{8\pi \beta C^2}{d^2} - \frac{24\pi \beta C^3}{d^3}$$

Siła (F = -dV/dd, dodatnia = odpychanie):

$$F_2(d) = -\frac{4\pi C^2}{d^2} + \frac{16\pi \beta C^2}{d^3} - \frac{72\pi \beta C^3}{d^4}$$

### 3.2 Potencjał 3-ciałowy (źródło: człon Φ⁴, nielinearna TGP)

Z rozwinięcia $(1 + \delta_1 + \delta_2 + \delta_3)^4$ z współczynnikiem multinomialnym $4!/(1!1!1!1!) = 24$:

$$V_3 = -6\gamma C_1 C_2 C_3 \cdot I_{\text{triple}}$$

Gdzie całka potrójnego nakładania profili Yukawa:

$$I_{\text{triple}} = \int \frac{e^{-m r_1}}{r_1} \cdot \frac{e^{-m r_2}}{r_2} \cdot \frac{e^{-m r_3}}{r_3} \, d^3x$$

**Granica Coulomba** (m_sp → 0, praktyczna dla TGP fizycznego):

$$\boxed{I_{\text{triple}}^{\text{Coul}} = \frac{8\pi^2}{P}, \quad P = d_{12} + d_{13} + d_{23}}$$

Zatem:

$$V_3^{\text{Coul}} = -\frac{48\pi^2 \gamma C_1 C_2 C_3}{d_{12} + d_{13} + d_{23}}$$

Siła 3-ciałowa na cząstkę i (ANALITYCZNA w granicy Coulomba):

$$\vec{F}_3^{(i)} = \frac{48\pi^2 \gamma C_i C_j C_k}{P^2} \left(\hat{r}_{ij} + \hat{r}_{ik}\right)$$

gdzie $\hat{r}_{ij} = (\vec{x}_j - \vec{x}_i)/d_{ij}$ — wersor od i do j.

### 3.3 Warunek równowagi — Lemat ogólny (z tgp_nbody_special_configurations.tex)

Dla dowolnej **konfiguracji jednoskalowej** (geometria stała, parametr skali λ, sumy geometryczne S₁, S₂, S₃):

$$\boxed{S_1 \lambda^2 - 4\beta S_2 \lambda + 18\gamma C S_3 = 0}$$

Konkretne przypadki:

| Konfiguracja | Równanie równowagi |
|---|---|
| 2 ciała (d) | $d^2 - 4\beta d + 18\beta C = 0$ |
| Trójkąt równoboczny (d) | $d^2 - 4\beta d + 18\beta C = 0$ |
| Tetraedr foremny (a) | $a^2 - 4\beta a + 18\beta C = 0$ |
| Symetryczny kollinear (x) | $10x^2 - 36\beta x + 153\beta C = 0$ |
| Kwadrat (a) | $(4+\sqrt{2})a^2 - 20\beta a + 18\beta C(4+1/\sqrt{2}) = 0$ |

**Kluczowy wniosek:** 2 ciała, trójkąt i tetraedr dają IDENTYCZNE równanie (bo wszystkie pary są równoodległe, S₁=S₂=S₃=n_par). Kollinear i kwadrat dają INNE równania.

### 3.4 Rozwiązania i warunki istnienia

Dla 2 ciał / trójkąta / tetraedru (γ=β):

$$d_{1,2} = 2\beta \pm \sqrt{4\beta^2 - 18\beta C}$$

**Warunek istnienia równowag:**
$$\beta > \frac{9C}{2} \iff \frac{\beta}{C} > 4.5$$

Dla symetrycznego kollinear:

$$x_{1,2} = \frac{18\beta \pm \sqrt{324\beta^2 - 1530\beta C}}{10}$$

**Warunek:** dyskryminant > 0:
$$324\beta^2 > 1530\beta C \iff \frac{\beta}{C} > \frac{1530}{324} \approx 4.72$$

**Wniosek: trójkąt równoboczny jest STABILNIEJSZĄ konfigracją (niższy próg β/C = 4.5 vs 4.72)**

---

## 4. Co realnie można wyciągnąć z problemu 3 ciał TGP

### 4.1 WYNIKI DOKŁADNE (analityczne)

#### A. Diagram fazowy (β, C)

```
Regiony w przestrzeni (β, C):
  - β/C < 4.5:   brak równowag statycznych, tylko rozproszenie
  - β/C ∈ [4.5, 4.72): trójkąt równoboczny możliwy, kollinear niemożliwy
  - β/C > 4.72:  OBY geometrie mają równowagi
```

Dla każdego regionu: d_rep (bariera) i d_well (studnia) jako f(β, C).

#### B. Selekcja konfiguracji przez teorię

TGP **selekcjonuje** trójkąt równoboczny przed kollinearem. Jest to odwrotność Newtona gdzie obie Lagrange'a (L4/L5) i Eulera (L1/L2/L3) istnieją dla wszystkich stosunków mas.

W TGP z 3 RÓWNYMI masami:
- **Equilateral**: pierwsza konfiguracja, która pojawia się przy wzroście β/C
- **Collinear**: pojawia się przy wyższym β/C ≈ 4.72

#### C. Częstości oscylacji przy d_well (mody normalne)

Dla trójkąta równobocznego przy d_well, Hessian ma 3 fizyczne mody:

1. **Tryb oddechu** (breathing mode): wszystkie 3 masy zbliżają/oddalają się symetrycznie
   - $\omega_{\text{breath}}^2 = H_{\text{scalar}} / C_{\text{eff}}$
   - $H_{\text{scalar}} = 3 \cdot \left[-\frac{8\pi C^2}{d^3} + \frac{48\pi\beta C^2}{d^4} - \frac{288\pi\beta C^3}{d^5}\right]_{d=d_{\text{well}}}$

2. **Tryby ścinania** (shear modes, 2-krotnie zdegenerowany): deformacja trójkąta przy stałym obwodzie

Wzór analityczny dla trybu oddechu (używając warunku równowagi d_well):

$$\omega_{\text{breath}}^2 = \frac{12\pi C}{d_{\text{well}}^4} \cdot \left[4\beta d_{\text{well}} - 2d_{\text{well}}^2 - 18\beta C\right]$$

Po podstawieniu $d_{\text{well}} = 2\beta + \sqrt{4\beta^2-18\beta C}$: wyrażenie zależy tylko od β i C.

#### D. Kryterium Hill dla ograniczonych orbit

W granicy Coulomba (m_sp=0) dla trójkąta równobocznego:

$$V_{\text{total}}(d) = 3V_2(d) + V_3(d,d,d) = -\frac{12\pi C^2 + 16\pi^2\gamma C^3}{d} + \frac{24\pi\beta C^2}{d^2} - \frac{72\pi\beta C^3}{d^3}$$

Nota: V₃ modyfikuje człon atrakcyjny 1/d. Dla ograniczonej orbity:

$$E_{\text{total}} = T + V_{\text{total}} < \lim_{d\to\infty} V_{\text{total}}(d) = 0$$

Zatem: ruch ograniczony gdy $E < 0$ → $T < |V_{\text{total}}(d_0)|$.

### 4.2 WYNIKI PÓŁANALITYCZNE

#### E. Korekcja 3-ciałowa do równowagi

Perturbacyjna korekta d₀ → d* przy uwzględnieniu V₃:

$$\delta d = d^* - d_0 = -\frac{F_3}{H_2^{(2b)}}$$

Dla trójkąta równobocznego (Coulomb):

$$F_3 = \frac{d V_3}{d d}\bigg|_{d=d_0} = \frac{16\pi^2\gamma C^3}{d_0^2} \cdot 3 \quad (\text{siła odpychająca})$$

$$H_2^{(2b)} = 3 \left[-\frac{8\pi C^2}{d_0^3} + \frac{48\pi\beta C^2}{d_0^4} - \frac{288\pi\beta C^3}{d_0^5}\right]$$

Parametr małości: $\epsilon = |V_3/V_2| \approx 4\pi\gamma C / (d_{\text{well}})$. Dla $\gamma \ll 1$ (TGP fizyczny): korekta znikoma.

#### F. Porównanie z Yukawa

Dla m_sp·d > 1 (silne ekranowanie):
$$I_{\text{triple}}^{\text{Yuk}} \approx \frac{(2\pi)^{3/2}}{m_{\text{sp}}^2} \cdot \frac{e^{-m_{\text{sp}} s}}{s}, \quad s = \frac{P}{2}$$

Przybliżenie saddle-point dokładne do ~5% dla m_sp·d > 1. Dla TGP fizycznego (m_sp → 0) użyć granicy Coulomba.

### 4.3 WYNIKI NUMERYCZNE (z kodem)

#### G. Analiza stabilności przez Hessian

Kod `stability.py` (rigorous zero-mode projection) daje:
- Eigenvalues fizycznych modów po eliminacji translacji + rotacji
- `stable` gdy wszystkie eigenvalues > 0
- Częstości drgań ω = √(eigenvalue)

**Wynik potwierdzony numerycznie:**
Przy d_well (β/C > 4.5): trójkąt STABILNY → **naruszenie twierdzenia Earnshawa ✓**
Przy d_rep: trójkąt NIESTABILNY (siodło) → potwierdza fizykę bariery

#### H. Wyniki numeryczne (historyczny batch; oryginalny skrypt `run_3body_complete.py` już nie jest w repozytorium)

**Tabela równowag i częstości** (C=0.3, β=γ=vacuum):

| β/C | d_rep | d_well | V(d_well) | ω_breath | V₃/V₂(Yukawa) | d_coll |
|---|---|---|---|---|---|---|
| 4.6 | 2.353 | 3.167 | -0.405 | 0.2484 | 0.50% | — |
| 5.0 | 2.051 | 3.949 | -0.356 | 0.2894 | 0.09% | 3.336 |
| 6.0 | 1.800 | 5.400 | -0.279 | 0.2186 | ~0% | 4.735 |
| 7.0 | 1.690 | 6.710 | -0.232 | 0.1674 | ~0% | 5.936 |
| 10.0 | 1.550 | 10.450 | -0.154 | 0.0919 | ~0% | 9.323 |

**Korekcja 3-ciałowa (Yukawa):**
- β/C=5.0: δd/d_well = 0.61% ← perturbacyjna ✓
- β/C=7.0+: δd → 0 (eksponencjalnie stłumiona) ✓

**Zachowanie energii leapfrog:** max|ΔE/E₀| = 2×10⁻¹⁰ (doskonałe) ✓

#### I. Dynamika czasowa

Integrator leapfrog (symplektyczny):
- Doskonała długoterminowa konserwacja energii
- Brak kolapsów dla d > d_rep (bariera 1/d²)
- Trajektorie: oscylacje wokół d_well lub ucieczka przez d_rep

---

## 5. Porównanie TGP vs Newton — 3 ciała

| Właściwość | Newton | TGP |
|---|---|---|
| Statyczne równowagi | **NIEMOŻLIWE** (Earnshaw) | **MOŻLIWE** (β/C > 4.5) |
| Kolaps do r=0 | Nieunikniony | **Blokowany** przez V_beta ~ 1/d² |
| Siły 3-ciałowe | Brak (tylko pairwise 1/r²) | **V₃ ~ γC³/P** (niezerowe) |
| Zasada superpozycji | Dokładna | Przybliżona (V₃ ≠ 0) |
| Chaos | Silny (osobliwość r→0) | Zredukowany (regularna bariera) |
| Konfiguracja Lagrange'a | Metastabilna (wymaga mas) | **Stabilna** dla 3 równych mas |
| Konfiguracja Eulera | Zawsze niestabilna | Potencjalnie stabilna |

### Kluczowa predykcja TGP dla 3 ciał:

**Trzy równe masy w trójkącie równobocznym o boku d_well MOGĄ STAĆ W RÓWNOWADZE STATYCZNEJ.**
W Newtonie to NIEMOŻLIWE. Jest to bezpośrednia konsekwencja bariery repulsywnej i studni konfinującej TGP.

---

## 6. Hierarchia trudności: otwarte problemy

### Poziom 1: Rozwiązane (można wyciągnąć teraz)

- [x] Diagram fazowy (β, C) — algebraiczny
- [x] d_rep, d_well jako funkcje β, C — algebraiczne
- [x] Warunki istnienia równowag — algebraiczne
- [x] Stabilność przy równowagach — Hessian numeryczny
- [x] Częstości drgań — z Hessiana
- [x] Korekcja 3-ciałowa (perturbacyjna)

### Poziom 2: Dostępne numerycznie (wymaga symulacji)

- [ ] Lyapunov exponents — miara chaosu w TGP vs Newton
- [ ] Mapa Poincaré — struktura trajektorii
- [ ] Czas życia konfiguracji nierównowagowych
- [ ] Orbita ósemkowa Chenciner-Montgomery w TGP

### Poziom 3: Trudne (otwarte teoretycznie)

- [ ] Zamknięta forma I_triple dla Yukawa (ogólna geometria) — diagram trójkątny
- [ ] Kryterium Hill analog dla ogólnej konfiguracji 3-ciał
- [ ] Pełna klasyfikacja orbit zamkniętych w TGP
- [ ] Korekcja relatywistyczna (3-ciała z metryką TGP)

---

## 7. Konkretne Predykcje Obserwacyjne

### P1: Trojany TGP (analogia Trojany Jowisza)

W Newtonie trojany na L4/L5 są METASTABILNE (wymagają dużego stosunku mas). W TGP z 3 równymi masami:
- Równowaga istnieje gdy β/C > 4.5
- Jest **STATYCZNIE STABILNA** (nie wymaga obrotu!)
- Drgania wokół d_well: ω ~ √(H_scalar/C)

### P2: Skalowanie V₃/V₂

$$\frac{|V_3|}{|V_2|} = \frac{4\pi\gamma C}{1} \cdot \frac{\text{(geometryczny czynnik)}}{\text{(dla Coulomba)}}$$

Dla TGP fizycznego γ ~ 10⁻⁵² m⁻², C ~ masa/Φ₀² → V₃/V₂ ekstremalnie małe dla cząstek elementarnych, ale **niezerowe** i kumulatywne dla N → ∞ (galaktyki).

### P3: Brak chaosu dla d > d_rep

Dla konfiguracji z E < V(d_rep), trajektorie w TGP są **ograniczone i regularne**. Newton daje chaos. To jest FALSYFIKOWALNA PREDYKCJA dla symulacji N-body.

---

## 8. Kierunki pracy (zsynchronizowane z `PLAN_ROZWOJU_NBODY.md`)

Szczegóły priorytetów (Lyapunov, pełny Hessian, Poincaré, domknięcie \(I_{\text{triple}}\)): [PLAN_ROZWOJU_NBODY.md](PLAN_ROZWOJU_NBODY.md).

**Minimum pod dalsze numeryki w 3 ciałach:**

1. `dynamics_v2.py` + nowe skrypty w `examples/` — skany w przestrzeni \((\beta, C)\), mapy stabilności, porównanie TGP vs Newton przy tej samej energii.
2. Konfiguracje startowe z `equilibria.py` / `configurations.py` — zaburzenia wokół \(d_{\rm well}\).
3. Ewentualnie rozszerzenie `examples/_archiwum/ex45_lyapunov_tgp_vs_newton.py` lub nowy skrypt w `examples/` wg planu P1.

---

## 9. Równania do sprawdzenia spójności (rdzeń 3-ciałowy)

Wszystkie poniższe muszą być spełnione jednocześnie:

```
[1] V_gamma = -12π γ C₁C₂(C₁+C₂)/d³   ← z pairwise.py ✓
[2] F_gamma = -36π γ C₁C₂(C₁+C₂)/d⁴  ← F = -dV/dd ✓
[3] d² - 4βd + 18βC = 0               ← równowaga 2-ciał/trójkąt ✓
[4] V₃ = -6γC₁C₂C₃ · I_triple        ← z Φ⁴ nonlinearity ✓
[5] I_triple(Coulomb) = 8π²/P          ← dokładne ✓
[6] dI/dP = -8π²/P²                    ← pochodna ✓
[7] m_i = C_i                          ← aksjomat TGP (ekwiwalencja) ✓
[8] c(Φ) = c₀√(Φ₀/Φ)                 ← z metryki emergentnej ✓
```

Wszystkie [1]–[8] pozostają spójne z rdzeniem modułów Pythona w `nbody/` (stan 2026-04).

---

*Skrót utrzymaniowy: rdzeń `nbody/` = 8 modułów `.py` + 22 pliki `.tex` w tym katalogu; skrypty `examples/*.py` (setki eksperymentów, część w `_archiwum/`). Pełna mapa teorii: `ANALIZA_SPOJNOSCI_v5.md` w korzeniu `TGP_v1/`.*
