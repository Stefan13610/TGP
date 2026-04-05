# ex122 — Człon krzyżowy \((\alpha/g)(g')^2\) względem \(V'(g)\)

Równanie (ex114 / ex106):

\[
g'' = \frac{1}{f(g)}\Bigl( V'(g) - \frac{\alpha}{g}(g')^2 - \frac{2 f(g)\,g'}{r} \Bigr).
\]

**Linearyzacja w \(h=g-1\)** w `ex117` pomija przy prowadzeniu równania na \(h\) nie tylko człony wyższego rzędu w \(V'\) i \(f\), ale traktuje **osobno** człon \(\mathcal{C} := (\alpha/g)(g')^2\): w „czystym” równaniu na \(h\) zakładamy, że w ogonie \(\mathcal{C}\) jest mały w porównaniu do skali \(V'(g)\approx -h\) (lub do \(\|L[h]\|\) po przeniesieniu na lewą stronę).

**Pierwszy sensowny porządek porównania (bez przenoszenia na \(h\)):**

\[
\mathcal{R}_{CV}(r) := \frac{\bigl|(\alpha/g)(g')^2\bigr|}{\bigl|V'(g)\bigr|}
= \frac{\alpha (g')^2}{g\,|g^2(1-g)|}.
\]

Dla \(g\) blisko 1 mianownik \(V'\) bywa bardzo mały — **surowe maksimum** \(\mathcal{R}_{CV}\) jest podatne na piki; w ex122 raportujemy też **percentyl 99** i medianę na \([r_L,r_R]\), z **\(g'\)** z solvera (nie z różnicowania szumu).

Jeśli mediana \(\mathcal{R}_{CV}\ll 1\) na \([25,32]\) dla triady, uzasadnia to traktowanie linearyzacji prowadzącej do \(h''+2h'/r+h\approx 0\) jako **pierwszego sensownego przybliżenia** w tym paśmie (obok już zmierzonego \(\zeta\) z ex117).

---

*Powiązane: `ex122_cross_term_ratio.py`, `ex117_tail_linearization.md`.*
