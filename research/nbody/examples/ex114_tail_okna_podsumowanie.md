# ex114 — okna tail: wyprowadzenia algebraiczne i domknięcie operacyjne

Dokument towarzyszy `ex114_tail_phase_map.py` i wynikom w `_outputs/`.

---

## 1. Algebra fitu (bez ODE)

Przy ustalonym oknie \(r \in [r_L, r_R]\) definiujemy \(y(r) := (g(r)-1)\,r\). Ansatz regresji:

\[
y(r) = B\cos r + C\sin r + \varepsilon(r),
\]

gdzie \(\varepsilon\) jest resztą (nieliniowość, wyższe mod, skończony \(R_{\max}\), itd.).

**Amplituda i faza (w tej normalizacji):**

\[
A_{\rm tail} = \sqrt{B^2 + C^2}, \qquad
\phi = \operatorname{atan2}(C, B).
\]

**Odległość punktu \((B,C)\) od prostej \(B=C\)** (metryka euklidesowa w płaszczyźnie \((B,C)\)):

\[
\mathrm{dist}(B{=}C) = \frac{|B-C|}{\sqrt{2}}.
\]

Bezwymiarowy wskaźnik względem amplitudy: \(\mathrm{dist}/A_{\rm tail}\) (tak jak w CSV).

**Warunek \(C=0\)** to faza \(\phi \in \{0,\pi\}\) (oś \(B\)); **\(B=0\)** to \(\phi = \pm\pi/2\) (oś \(C\)).  
**\(B=C\)** to \(\phi = \pi/4 + k\pi\) (prosta pod \(45^\circ\)).

To są **czysto geometryczne** stwierdzenia — nie wymagają znajomości pełnego ODE.

---

## 2. Heurystyka asymptotyki (łącze z ODE)

Pełne równanie na \(g(r)\) w ujęciu ex106/ex114 jest **nieliniowe** (potencjał \(V'(g)\), człon \(\propto g'^2\), sprzężenie \(f(g)\), sferyczny \(2g'/r\)).  
W **zewnętrznej strefie**, gdzie \(g \approx 1\) i nieliniowość jest mała, oczekuje się zachowania zbliżonego do **oscylacji radialnej** modulowanej przez \(1/r\) — stąd naturalna próba dopasowania \(y = B\cos r + C\sin r\) na skończonym przedziale.

To **nie** jest twierdzenie „ogon = dokładnie cos/sin”; to **model pomiarowy** z jawną resztą \(\varepsilon\). Jakość fitu mierzymy przez **RMSE** w oknie i przez **stabilność wyniku przy przesunięciu okna** (patrz §3).

---

## 3. Domknięcie operacyjne wyboru okna

Ponieważ nie mamy zamkniętej formy \((B,C)(g_0)\) z równania, okno \([r_L,r_R]\) domykamy **procedurą**, nie jedną formułą:

| Kryterium | Implementacja w ex114 | Interpretacja |
|-----------|------------------------|----------------|
| Jakość ansatzu w oknie | P6: RMSE\(/A < 10\%\) | Reszta \(\varepsilon\) nie dominuje amplitudy |
| Stabilność amplitudy/fazy przy przesunięciu | P4–P5 na siatce \(g_0\) | Mały rozstrzał między \([20,35]\) a \([22,36]\) |
| Stabilność triady | P8: \(\lvert\Delta\phi\rvert_{\rm triada} < 0{,}35\) rad | e, μ, τ nie są artefaktem jednego przecięcia |
| Asymptotyka \(g\to 1\) w oknie | P9: średnie \(\lvert g-1\rvert\) na triadzie | Okno leży w strefie, gdzie tło jest bliskie próżni \(g=1\) |

Jeśli P4–P6, P8 i P9 przechodzą, **nie udowadniamy**, że \([20,35]\) jest optymalne globalnie — ale **uzasadniamy użycie** tego okna jako spójnego pomiaru wobec sąsiedniego okna i wobec małych odchyleń od \(g=1\).

### 3.1. \([20,35]\) vs „best” z ex115

| Pojęcie | Rola |
|--------|------|
| **Okno referencyjne \([20,35]\)** | Kotwica **protokołu** ex114 oraz ex116–ex120: jawne progi P4–P9, ten sam \(r_L,r_R\) przy skanach i overlayu. |
| **Okno „best” w ex115** | Wśród wierszy z `pass_all=True` wybierane jest minimum **max(RMSE/A)** po triadzie — empirycznie często **\([22{,}5,\,32{,}5]\)** (węższy pas, niższy RMSE/A), nadal PASS. |

**Decyzja dokumentowa:** łańcuch regresji (digest ex119, Pearson ex118, runner ex121) **nie jest** automatycznie przenoszony na „best” bez osobnej decyzji: \([20,35]\) utrzymuje **porównywalność** z wcześniejszymi eksperymentami i z kryteriami ε w ex116. Jeśli kiedyś ujednolicisz protokół na \([22{,}5,\,32{,}5]\), zrób to jako świadomy krok (stałe `R_TAIL_*` w ex114/ex106 + ponowny przebieg suite). Wizualnie: `ex120_ex115_scan_overlay.png`, CSV `ex115_window_scan.csv`.

---

## 4. Co dalej (jeśli chcesz „twardejsze” okno)

1. **Skan \((r_L, r_R)\)** — zaimplementowany w `ex115_tail_window_scan.py`: CSV `ex115_window_scan.csv`, mapa `ex115_window_scan_admissible.png` (maska `pass_all` + heatmapa max(RMSE/A)), triada vs fazy z okna referencyjnego \([20,35]\).  
2. **Warunek startu okna (średnia na odcinku):** `ex116_tail_rL_from_epsilon.py` — dla ustalonego \(r_R\) szuka **najmniejszego** \(r_L\), przy którym **jednocześnie** na \([r_L,r_R]\) średnia \(\lvert g-1\rvert < \varepsilon\) dla e, μ, τ; wynik w `ex116_rL_vs_epsilon.csv`. Uwaga: to kryterium **nie** jest tożsame z RMSE ansatzu cos/sin (ex114 P6).  
3. **Linearyzacja jawna** — `ex117_tail_linearization.md`: \(h'' + \frac{2}{r}h' + h \approx 0\) oraz związek z \((g-1)r \approx B\cos r + C\sin r\). Reszta \(L[h]\): `ex117_linear_operator_residual.py`. Korelacja \(\zeta\) ↔ RMSE/A: `ex118_rmse_vs_linear_residual.py` → `ex118_rmse_vs_zeta.csv`.

---

## 5. Digest łańcucha (jedno uruchomienie)

`ex119_tail_pipeline_digest.py` zapisuje `ex119_pipeline_digest.csv` i `ex119_pipeline_digest.png` — zbiór metryk z triady (ζ, RMSE/A, \(r_L^*(\varepsilon)\), Pearson) bez ręcznego sklejania wyników ex114–ex118.

`ex120_ex115_scan_overlay.py` wczytuje `ex115_window_scan.csv` (ew. odpala ex115), rysuje `ex120_ex115_overlay.png` oraz `ex120_ex115_overlay_meta.csv` — PASS/FAIL, punkt referencyjny \([20,35]\) i okno **najlepsze po min(RMSE/A)** wśród PASS.

**Pełny przebieg:** `ex121_tail_suite_runner.py` uruchamia ex114→ex120, **ex122**, **ex123**, **ex124** i zapisuje `ex121_tail_suite_log.txt` (czasy, kody wyjścia). Flaga `--dry-run` tylko wypisuje listę. **Plan fazowy:** `PLAN_TAIL_DOMKNIECIE_ROZWOJ.md`. **Gęsta siatka g₀ + dwa solvery:** `ex124_dense_g0_solver_compare.py` → `ex124_dense_solver_compare.csv`.

---

*Tekst roboczy, spójny z ex114 (2026-04-02).*
