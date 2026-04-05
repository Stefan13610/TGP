# Plan domknięcia i rozwoju — łańcuch tail (ex114–ex121) + analityka

**Status:** Faza C domknięta operacyjnie (2026-04-02); Faza D nadal poza `nbody` na żądanie.  
**Zakres:** `nbody/examples/` (ogon solitonowy, okna, linearyzacja); powiązania z `ex106` / Ścieżką 9.

---

## Faza A — Domknięcie teoretyczne (co „jest w LaTeX”, co w kodzie)

| ID | Cel | Stan | Następny krok |
|----|-----|------|----------------|
| A1 | Jednoznaczna postać ODE w notacji zgodnej z `main`/dodatekJ2 | Częściowo | Skopiować równoważność `f g'' + 2f g'/r = V' - (α/g)g'²` do jednego pliku `.tex` przy ogonie (poza `nbody` — gdy zdecydujesz) |
| A2 | Linearyzacja \(g=1+h\) → \(h''+2h'/r+h\approx 0\) | **Zamknięte** w `ex117_tail_linearization.md` + **§4.1** (ex122, triada) | — |
| A3 | Algebra \((B,C)\), fazy, proste \(B=C\) | **Zamknięte** w `ex114_tail_okna_podsumowanie.md` | — |
| A4 | Koide jako algebra przy danym \(r_{21}\) | Zewnętrzne względem tail | Utrzymać rozdział: Koide ≠ wyprowadzenie z ODE |

---

## Faza B — Domknięcie numeryczne (regresje, jeden runner)

| ID | Cel | Stan |
|----|-----|------|
| B1 | Mapa \((B,C,A)(g_0)\), triada, dwa okna | **ex114** |
| B2 | Skan \((r_L,r_R)\), PASS | **ex115** |
| B3 | \(r_L^*(\varepsilon)\) | **ex116** |
| B4 | \(\zeta = \|L[h]\|/\|h\|\) | **ex117** |
| B5 | Pearson \(\zeta\) ↔ RMSE/A | **ex118** |
| B6 | Digest CSV/PNG | **ex119** |
| B7 | Overlay ex115 | **ex120** |
| B8 | Runner całości | **ex121** |
| B9 | **Skala członu krzyżowego** \((\alpha/g)g'^2\) vs \(V'\) | **ex122** (`ex122_cross_term_ratio.py`, `ex122_tail_cross_term_sketch.md`) |

**Kryterium „domknięcia operacyjnego”:** `python ex121_tail_suite_runner.py` → exit 0, log w `_outputs/ex121_tail_suite_log.txt` (obejmuje ex122–ex124; ex124 ~1 min).

---

## Faza C — Rozwój badawczy (kolejność sugerowana)

1. ~~**Mapa fazowa na gęstszej siatce \(g_0\)** + drugi solver~~ **Zrobione:** `ex124_dense_g0_solver_compare.py` — CSV `ex124_dense_solver_compare.csv`, PNG `ex124_dense_A_compare.png` (reg vs \(K_{\mathrm{sub}}=g^2\)).  
2. ~~**Niezależność od Koida**~~ **Zrobione:** `ex123_koide_epistemics.py` + `ex123_koide_epistemics.md` — tabela \(r_{31}^{(\varphi^2)}\) vs \(r_{31}^{(\mathrm{Koide})}\) przy tym samym \(r_{21}\) ze Ścieżki 9.  
3. ~~**Powiązanie `best` z ex115**~~ **Zrobione:** §3.1 w `ex114_tail_okna_podsumowanie.md` — \([20,35]\) jako kotwica vs „best” (np. \([22.5,32.5]\)).  
4. ~~**Linearyzacja zachowawcza**~~ **Zrobione:** `ex117_tail_linearization.md` §4.1 (mediana/p99 ex122, e/μ vs τ).

---

## Faza D — Poza `nbody` (świadomie odłożone)

- Wpięcie skrótu planu lub odnośników do `dodatekJ2` / `ANALIZA_SPOJNOSCI` — tylko na Twoje polecenie.  
- Publikacja: jedna sekcja „Numerical protocol for tail extraction” z odnośnikiem do ex121.

---

*Ten plik można aktualizować po każdej zamkniętej fazie (data + krótka notka w tabeli).*
