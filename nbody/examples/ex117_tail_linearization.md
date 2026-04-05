# ex117 — Linearyzacja ogona: \(h = g-1\) przy \(g \approx 1\)

Powiązanie z ansatzem ex114/ex106: \((g-1)r \approx B\cos r + C\sin r\).

---

## 1. Równanie w postaci ex114

Z `ex114_tail_phase_map.py` (zgodne z ex106):

\[
f(g)\,g'' + \frac{2 f(g)\,g'}{r}
= V'(g) - \frac{\alpha}{g}(g')^2,
\qquad
V'(g)=g^2(1-g),
\qquad
f(g)=1+2\alpha\ln g,
\qquad \alpha=2.
\]

Równoważnie:
\[
g'' = \frac{1}{f(g)}\left( V'(g) - \frac{\alpha}{g}(g')^2 - \frac{2 f(g)\,g'}{r} \right).
\]

---

## 2. Rozwinięcie przy \(g = 1 + h\), \(|h|\ll 1\)

**Potencjał:**
\[
V'(1+h) = (1+h)^2(-h) = -h - 2h^2 - h^3.
\]
W pierwszym rzędzie: \(V'(g) \approx -h\).

**Sprzężenie kinetyczne:**
\[
f(1+h) = 1 + 2\alpha\ln(1+h) = 1 + 2\alpha h + O(h^2).
\]

**Lewa strona** (do pierwszego rzędu w \(h\), przy traktowaniu \(h',h''\) jako tego samego rzędu co \(h\) w strefie asymptotycznej):
\[
f(g)g'' + \frac{2 f(g) g'}{r}
\approx g'' + \frac{2g'}{r}
= h'' + \frac{2h'}{r} + O(\alpha h\,h'',\, \alpha h\,h'/r).
\]

**Prawa strona:** człon \(\frac{\alpha}{g}(g')^2\) jest **kwadratowy** w gradientach; w strefie, gdzie \(h'\) jest małe w porównaniu ze skalą \(h\) (zewnętrzny ogon), pierwszym przybliżeniem jest pominięcie tego członu wobec \(V'\approx -h\).

**Prowadzące równanie liniowe:**
\[
\boxed{h'' + \frac{2}{r} h' + h \approx 0.}
\]

---

## 3. Rozwiązania i związek z fitowaniem ex114

Podstawienie \(u(r) = r\,h(r)\) daje
\[
u''(r) + u(r) = 0
\quad\Rightarrow\quad
u = A\cos r + B\sin r,
\]
czyli
\[
h(r) = \frac{A\cos r + B\sin r}{r},
\qquad
(g-1)\,r = A\cos r + B\sin r.
\]

To jest **dokładnie** forma używana w `fit_tail_detailed` (nazewnictwo \(B,C\) w kodzie odpowiada tym samym liniowym kombinacjom cos/sin).

---

## 4. Zakres ważności (epistemika)

- Linearyzacja \(V'\) i \(f\) wymaga \(|h|\ll 1\) — mierzone ex114 P9 / ex116 (średnie \(|g-1|\) w oknie).
- Pominięcie \(\frac{\alpha}{g}(g')^2\) wymaga, by w porównaniu do \(|h|\) było małe — to **hipoteza** na odcinku tail; weryfikacja: `ex117_linear_operator_residual.py` (ilorazy resztowe).
- Pełne \(f(g)\neq 1\) i człony \(O(\alpha h)\) dają poprawki do równania liniowego — odpowiadają części „RMSE” fitu cos/sin (ex114 P6).

**Korelacja empiryczna (triada, dwa pasma):** `ex118_rmse_vs_linear_residual.py` — na zbiorze 6 punktów współczynnik Pearsona między \(\zeta\) a RMSE/A bywa bardzo wysoki (\(\sim 1\)), co sugeruje, że na tym protokole reszta liniowa i jakość fitu ansatzu idą **w parze**.

**Człon krzyżowy \((\alpha/g)(g')^2\) vs \(V'(g)\):** `ex122_tail_cross_term_sketch.md` + `ex122_cross_term_ratio.py` — mediana i p99 ilorazu \(|cross|/|V'|\) na odcinkach (triada); uzasadnia ostrożne traktowanie „surowego” maksimum przy \(|V'|\to 0\).

### 4.1. Numeryczne potwierdzenie zakresu (ex122, triada, protokół ex114)

Na podstawie przebiegu `ex122_cross_term_ratio.py` (okna \(r\in[20,35]\) oraz \([25,32]\)):

- **Elektron:** mediana \(|cross|/|V'|\) jest rzędu \(10^{-2}\) — w tym sensie człon krzyżowy **nie dominuje** nad \(V'\) w typowym punkcie ogonu.
- **Muon:** mediana pozostaje \(\ll 1\) (rzędu kilku–kilkudziesięciu procent w skali 1).
- **Tau:** mediana bywa **większa** (np. \(\sim 0{,}15\) na \([20,35]\), nieco mniej na \([25,32]\)) — linearyzacja z pominięciem cross jest dla τ **najsłabszym ogniwem** triady; nadal **p99** (a nie surowe max) służy w ex122 do testów stabilnych wobec \(V'\to 0\).

Wniosek zachowawczy: algebra z §2–3 jest **najlepiej uzasadniona** dla e i μ w wybranym oknie; dla τ należy traktować ją jako przybliżenie z kontrolowaną resztą (\(\zeta\), RMSE, ex122).

---

*ex117 — 2026-04-02, spójne z ex106/ex114; §4.1 uzupełniono ex122.*
