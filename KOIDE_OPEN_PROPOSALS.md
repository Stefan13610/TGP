# Koide w TGP --- status po 22 sciezkach analizy (v47b)

**Jedyne otwarte pytanie:** Dlaczego \(\mathrm{CV}(A^2) = 1\)?

Rownowazne sformulowania:
- Dlaczego \(Q_K = 3/2\)?
- Dlaczego \(b = \sqrt{2}\) w parametryzacji Brannena?
- Dlaczego hipoteza ekwipartycji jest prawdziwa?

---

## Co zostalo ustalone (22 sciezki, pewne)

### Mechanizm selekcji (prop:T-pythagorean-selection)

Lancuch dedukcyjny:
\[
g_0 \;\xrightarrow{\text{ODE}}\; (a(g_0), b(g_0)) \;\xrightarrow{\text{Pitagoras}}\; A^2 = a^2 + b^2 \;\xrightarrow{Q_K = 3/2}\; g_{0,\tau} = 1{,}5696
\]

- **Krok 1:** Ogon solitonu w \(d=3\) oscyluje: \(h(r) \to a\cos r + b\sin r\), dwie skladowe.
- **Krok 2:** \(A^2 = a^2 + b^2\) -- inwariancja obrotowa w plaszczyznie \((a,b)\).
- **Krok 3:** Mapa \(g_0 \to (a,b)\) kresli spirale -- amplituda rosnie, faza rotuje.
- **Krok 4:** \(Q_K = 3/2\) jednoznacznie wyznacza \(g_{0,\tau}\) w korytarzu kolapsu \([1{,}404, 1{,}600)\).

### Dlaczego masa (nie amplituda) spelnia Koide?

Poniewaz \(A^2 = a^2 + b^2\) usrednia po fazie spiralnej. Inwariancja pitagorejska laczy Koide z geometria w plaszczyznie \((a,b)\), ktora istnieje tylko dla \(d \geq 3\) (oscylacyjne ogony).

### Zbieznosc dwoch argumentow (rem:T-cv1-origin)

\[
\underbrace{2 \text{ skladowe ogona } (a,b)}_{\chi^2(2)} = \underbrace{N-1 = 2 \text{ niezalezne mody}}_{\text{Brannen}}
\]

Oba "2" to ta sama wielkosc: wymiernosc faz w oscylacyjnym ogonie w \(d=3\).

---

## Co zostalo wykluczone (zamkniete sciezki)

| Sciezka | Co testowano | Dlaczego zamknieta |
|---------|-------------|-------------------|
| 1-12 (rem:T-QK-derivation-attempts) | 12 roznych podejsc analitycznych | Negatywne wyniki w manuskrypcie |
| 13 (observer critical) | Efekt obserwatora | Tautologiczna -- Q_K=3/2 bylo wymuszone jako warunek na \(g_{0,\tau}\) |
| 14-17 | Rownowaga polowa, cross-sector, stabilnosc laczna | Negatywne, patrz manuskrypt |
| 18 (corridor) | Q_K monotoniczne w korytarzu | Brak ekstremum, brak samoselekcji |
| 19 (Z3 phases) | Fazy ODE 150/256/106 deg | Nie Z3, overlap lamie Koide |
| 20 (substrate direct) | Momenty Isinga 3D | Binder \(U_4 = 0.62 \neq 1/3\), CV = 0.36 |
| 21 (chi2 statistical) | Szum substratu | Nie stabilizuje Q_K, duzy rozrzut |
| 21 (chi2 structural) | Spirala ODE | Wyjasnia DLACZEGO masa, ale nie derywuje CV=1 |
| 22-1 (self-similarity) | \(A^2(\varphi g_0)/A^2(g_0)\) | Zmienia sie od 0 do 165, nie samopodobne |
| 22-2 (variational) | Energia, entropia, AM/GM | Wszystkie monotoniczne, brak ekstremum |
| 22-3 (analytical) | \(A^2 \approx C(g_0-1)^{2.08}/(g_{0c}-g_0)^{0.37}\) | 0.4% dokladnosc, nie implikuje CV=1 |
| 22-4 (sum rules) | Szybkosci wzrostu | Brak regul sumowych |
| 22-5 (discrete RG) | \(\varphi^2 g_0^e > g_{0,\mathrm{crit}}\) | Krok RG lamie sie, A^2 nie geometryczne |
| 22-6 (collapse distance) | \(\delta_\mu/\delta_\tau = 6.45\) | Brak dopasowania do znanych stalych |

---

## Co nadal otwarte

### A. JADRO: Dowod hipotezy ekwipartycji

**Pytanie:** Dlaczego \(N-1\) niezaleznych modow generacyjnych wnosi jednakowy wklad do wariancji \(A^2\)?

Hipoteza ekwipartycji (thm:T-QK-from-Ngen) daje \(b = \sqrt{N-1}\), stad \(\mathrm{CV}^2 = (N-1)/2 = 1\) dla \(N=3\).

**Mozliwe podejscia (nie testowane):**
1. Symetria Liego ukladu trzech sprzezonych solitonow (dawna propozycja #12)
2. Warunek na mody normalne trojsolitonowego ukladu
3. Dowod z wlasnosci analitycznych A(g0) a nie numerycznych

### B. ORTOGONALNE (oddzielne pytania)

- **Biegnace masy:** Czy istnieje skala \(\mu^\star\) przy ktorej \(Q_K\) jest stacjonarne wzgledem RG?
- **Kwarki:** Czy istnieje rozszerzona regula laczaca \(Q_K = 3/2\) dla leptonow z odchyleniem kwarkow?

### C. MATEMATYCZNE (nie bezposrednio Koide)

- **Dowod analityczny \(g_{0,\mathrm{crit}}(\alpha, d)\) dla \(d > 1\)**

---

## Pliki referencyjne

| Temat | Plik |
|-------|------|
| Manuskrypt formalny | `dodatekT_koide_atail_formal.tex` |
| Mechanizm selekcji | prop:T-pythagorean-selection |
| Zbieznosc chi2/Brannen | rem:T-cv1-origin |
| Spirala ODE | rem:T-chi2-closure |
| Korytarz kolapsu | rem:T-collapse-corridor |
| Archiwum pelnych propozycji | `_archiwum/KOIDE_OPEN_PROPOSALS_v47b_archive.md` |

### Skrypty (sciezki 18-22)

| Sciezka | Skrypt | Output |
|---------|--------|--------|
| 18 | `koide_observer_v47b.py` | `koide_observer_output.txt` |
| 18b | `koide_observer_critical_v47b.py` | `koide_observer_critical_output.txt` |
| 19 | `koide_collapse_corridor_v47b.py` | `koide_collapse_corridor_output.txt` |
| 20 | `koide_z3_correction_v47b.py` | `koide_z3_correction_output.txt` |
| 21 | `koide_substrate_cv_v47b.py` | `koide_substrate_cv_output.txt` |
| 21b | `koide_chi2_closure_fast.py` | `koide_chi2_closure_output.txt` |
| 22 | `koide_cv1_origin_v47b.py` | `koide_cv1_origin_output.txt` |

---

*Po 22 sciezkach analizy, problem Koide w TGP jest zredukowany do jednego pytania: dlaczego ekwipartycja modow generacyjnych. Wszystkie inne podejscia sa zamkniete lub subsumowane.*
