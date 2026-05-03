# Koide w TGP --- status po 23 sciezkach analizy (v47b)

**Status: Q_K = 3/2 jest KONSEKWENCJA STRUKTURALNA TGP, nie parametrem wejsciowym.**

---

## Lancuch dedukcyjny (zamkniety)

\[
d = 3 \;\to\; 2\text{ skladowe ogona} \;\to\; \chi^2(2) \;\to\; \mathrm{CV}(A^2) = 1 \;\to\; Q_K = 3/2
\]

### Szczegoly lancucha

1. **d=3 => oscylacyjne ogony** (prop:T-oscillatory-tail)
   Linearyzacja kanonicznego ODE wokol prozni daje sferyczne r. Bessela.
   W d=3: \(h(r) = (g(r)-1) \cdot r \to a\cos r + b\sin r\) — **dwie** niezalezne skladowe.

2. **Suma pitagorejska** (prop:T-pythagorean-selection)
   \(A^2 = a^2 + b^2\) — inwariancja obrotowa w plaszczyznie \((a,b)\).
   Masa \(m \sim A^4 = (a^2+b^2)^2\) zalezy od promienia, nie od fazy.

3. **N_gen = 3 z progu kolapsu** (prop:T-soliton-collapse)
   \(g_{0,\mathrm{crit}} = 8/5\) i \(\varphi^2 g_0^e = 2{,}27 > g_{0,\mathrm{crit}}\) =>
   tylko 3 generacje pod progiem.

4. **Zbieznosc chi2(2) i Brannena** (rem:T-cv1-origin)
   \[
   \underbrace{2 \text{ skladowe ogona } (a,b)}_{\chi^2(2)} = \underbrace{N-1 = 2 \text{ niezalezne mody}}_{\text{Brannen}}
   \]
   Obie "dwojki" to ta sama wielkosc. Rozklad \(\chi^2(2)\) ma CV = 1.

5. **CV = 1 <=> Q_K = 3/2** (thm:T-QK-CV)
   Algebraiczna rownowaznosc: \(Q_K = N/(1+\mathrm{CV}^2) = 3/2\).

### Zalozenie posredniczace

Hipoteza ekwipartycji (thm:T-QK-from-Ngen): \(N-1\) modow wnosi jednakowy wklad do wariancji.
Jest fizycznie motywowana i strukturalnie potwierdzona przez zbieznosc chi2(2)--Brannen.

---

## Pelny mechanizm selekcji

\[
g_0 \;\xrightarrow{\text{ODE}}\; (a(g_0), b(g_0)) \;\xrightarrow{\text{Pitagoras}}\; A^2 = a^2 + b^2 \;\xrightarrow{Q_K = 3/2}\; g_{0,\tau} = 1{,}5696
\]

- Mapa \(g_0 \to (a,b)\) kresli spirale w plaszczyznie \((a,b)\).
- Faza rotuje z przyspieszeniem ku kolapsowi (od 100 do 2000+ deg/dg0).
- \(A^2\) usrednia po fazie — dlatego **masa** (nie amplituda) spelnia Koide.
- Dodatkowy wynik: \(A^2 \propto I_{\mathrm{kin}}\) (calka kinetyczna) z wykladnikiem 1.000.

---

## Co zostalo wykluczone (zamkniete sciezki)

| Sciezka | Co testowano | Dlaczego zamknieta |
|---------|-------------|-------------------|
| 1-12 | 12 roznych podejsc analitycznych | Negatywne wyniki w manuskrypcie |
| 13 (observer) | Efekt obserwatora | Tautologiczna |
| 14-17 | Rownowaga polowa, cross-sector, stabilnosc | Negatywne |
| 18 (corridor) | Korytarz kolapsu | Q_K monotoniczne, brak ekstremum |
| 19 (Z3 phases) | Fazy Z3 | 150/256/106 deg, nie Z3 |
| 20 (substrate) | Momenty Isinga | Binder U4 = 0.62 != 1/3 |
| 21 (chi2 stat.) | Szum substratu | Nie stabilizuje Q_K |
| 21 (chi2 struct.) | Spirala ODE | Wyjasnia DLACZEGO masa — wlaczone do lancucha |
| 22-1..22-7 | Samopodoba, wariacyjne, entropia, RG, itp. | Wszystkie negatywne |
| 22-9 / 23 | Zbieznosc chi2/Brannen | **POZYTYWNA — zamyka lancuch** |

---

## Dwa wejscia TGP => pelne widmo mas

\[
\underbrace{K = g^4}_{\alpha = 2,\; g_{0,\mathrm{crit}} = 8/5,\; N_{\mathrm{gen}} = 3}
\;\;+\;\;
\underbrace{\varphi\text{-FP}}_{r_{21} = 206{,}77}
\;\;\xrightarrow{\;d=3,\;N=3\;\to\;Q_K = 3/2\;}
(m_e,\; m_\mu,\; m_\tau)
\]

---

## Pliki referencyjne

| Temat | Plik |
|-------|------|
| Manuskrypt formalny | `dodatekT_koide_atail_formal.tex` |
| Lancuch dedukcyjny | rem:T-cv1-origin |
| Mechanizm selekcji | prop:T-pythagorean-selection |
| Spirala ODE | rem:T-chi2-closure |
| Test ekwipartycji | rem:T-equipartition-test |
| Archiwum propozycji | `_archiwum/KOIDE_OPEN_PROPOSALS_v47b_archive.md` |

### Skrypty (sciezki 18-23)

| Sciezka | Skrypt |
|---------|--------|
| 18/18b | `koide_observer_v47b.py`, `koide_observer_critical_v47b.py` |
| 19 | `koide_collapse_corridor_v47b.py` |
| 20 | `koide_z3_correction_v47b.py` |
| 21 | `koide_substrate_cv_v47b.py` |
| 21b | `koide_chi2_closure_fast.py` |
| 22 | `koide_cv1_origin_v47b.py` |
| 23 | `koide_equipartition_v47b.py` |

---

*Q_K = 3/2 wynika z d=3 (dwie skladowe ogona) + N_gen=3 (dwa niezalezne mody) + ekwipartycja modow. Lancuch jest zamkniety.*
