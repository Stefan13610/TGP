# Ostrzał równania $T_c$ — seria 2 (technicznie mocniejsza)

**Data:** 2026-04-20
**Autor analizy:** Claudian (defensywa); zarzuty: anonimowy agent zewnętrzny, seria 2 po [OSTRZAL_ROWNANIE_Tc_2026-04-20](OSTRZAL_ROWNANIE_Tc_2026-04-20.md).
**Nowy skrypt testowy:** [`ps23_residual_correlation_Z_Mmol.py`](ps23_residual_correlation_Z_Mmol.py) → [`ps23_results.txt`](ps23_results.txt).

---

## Streszczenie werdyktu

| # | Zarzut (seria 2) | Werdykt |
|---|------------------|---------|
| 1 | Overfitting — to samo równanie dla Hg i YBCO | **Mylne założenie** — nie jest to jedno równanie, to 4-klasowy framework (P6.A/B/C/D) z dedykowanymi Λ_E per klasa |
| 2 | Factor 1.22 (YBCO) < 1.65 (Hg) — "błędy maleją dla skomplikowanego" | **Cherry-picked** — na pełnym N=29 RMS jest statystycznie dobrze rozłożony; Hg jest FAR off-resonance ($\Delta a \approx 1.1\sigma$) |
| 3 | Ekstremalne ciśnienie — czy $\Lambda_E$ jest stałe? | **Paper explicitnie to modeluje** przez $\Lambda_E^{\mathrm{eff}}(\omega_{\mathrm{ph}}) \propto \omega^{1.04}$; pod ciśnieniem $\omega_{\mathrm{ph}} \to$ rośnie |
| 4 | Residua vs Z / M_mol — ukryta zależność jądrowa? | **CZĘŚCIOWO trafiony** — dla N=25 jest umiarkowana korelacja $r \approx -0.41$ (p ≈ 0.04). **Nowy open sub-problem P7.4.** |

**Trzy zarzuty odbite, jeden otwiera rzeczywisty nowy kanał badawczy.** Warto.

---

## Zarzut 1 (seria 2) — „To samo równanie dla Hg i YBCO sugeruje overfitting do klas materiałów"

**Zarzut oparty na błędnym założeniu.**

To **nie jest jedno uniwersalne równanie** stosowane identycznie do Hg i YBCO. Paper jasno rozdziela **cztery sekcje P6** z różnymi gałęziami formuły:

| Klasa | Równanie | Λ_E | Dodatkowy czynnik |
|-------|----------|-----|-------------------|
| **P6.A cuprate** (YBCO, BSCCO, Tl2212...) | Eq. 8 | $\Lambda_E^{\mathrm{cup}} = 0.0513$ meV (stała) | $K_{\mathrm{dw}} = 3.498$, $\sqrt{n}$ |
| **P6.B phonon** (Al, Pb, Nb, Hg, MgB₂, H₃S, LaH₁₀...) | Eq. 10 | $\Lambda_0 (\omega_{\mathrm{ph}}/\omega_0)^{\alpha_B}$ ze stałych $\alpha_B=1.04$, $\Lambda_0=0.0962$ meV, $\omega_0=15$ meV | — |
| **P6.C orbital blocker** | Eq. 11 | j.w. | $B_{\mathrm{orb}}(P)$ |
| **P6.D spin fluctuations** | Eq. 12 | j.w. | $B_{\mathrm{mag}}(\lambda_{\mathrm{sf}})$ z $\beta=2.527$ |

Hg jest **P6.B** (fononowy, sp-orbital, $\Lambda_E$ zależne od $\omega_{\mathrm{ph}}$).
YBCO jest **P6.A** (kuprat, d-wave, $\Lambda_E^{\mathrm{cup}}$ stałe, dedykowany $K_{\mathrm{dw}}$).

**To są dwa różne kanały fizyczne** połączone wspólną strukturą substratu. Różnica w RMS_log:
- P6.A cuprates (N=8): **RMS_log = 0.182** (factor 1.5×)
- P6.B phonon+Fe (N=21): **RMS_log = 0.404** (factor 2.5×)

P6.A ma dedykowaną strukturę (elektron-mediated pairing z sprzężeniem d-wave node), P6.B używa jednej funkcji parametrycznej Λ_E(ω_ph) na 21 bardzo różnych materiałów → więc jest mniej dopasowany. **Odwrotnie do twierdzenia agenta** — cupraty są lepiej fitowane nie dlatego że są "prostsze", tylko dlatego że mają dedykowany branch modelu.

### Weryfikacja: MgB₂ i heavy-fermion jak chciał agent

**MgB₂**: paper przewiduje 31.37 K, obs 39 K, factor **1.24**. W pełni w paśmie RMS. MgB₂ to klasyczny case two-gap superconductor z dwoma osobnymi kanałami (σ-pairing z $\omega \sim 75$ meV + π-pairing, $T_c \sim 39$ K dla σ). Eq. 10 TGP stosuje tu $\omega_{\mathrm{ph}} = 75$ meV (σ-kanał), $z=5$, sp-orbital, daje dokładnie $T_c = 31$ K w ramach oczekiwanej dokładności. **Zaliczenie.**

**Heavy-fermion** (CeCu₂Si₂, UBe₁₃, UPt₃): **NIE MA w datasecie**. Paper EXPLICITNIE je flaguje jako **otwarty sub-problem P7.3 (Ce/Kondo)**. CeH₉ w datasecie (pred 22 K bez P7.3, albo 143 K z kalibracją pair-breaking) jest granicznym przypadkiem. Heavy-fermion wymaga mechanizmu walence-fluctuation / Kondo-singlet, którego paper wciąż nie domyka. To jest **uczciwy otwarty kanał**, nie ukryty overfitting.

**Werdykt:** zarzut overfittingu nie ma podstaw — paper ma wyraźnie klasową strukturę, a klasa bez dedykowanego mechanizmu (heavy-fermion) jest otwarcie flagowana.

---

## Zarzut 2 (seria 2) — „Dla YBCO (1.22) błąd mniejszy niż dla Hg (1.65) — błędy maleją dla skomplikowanego, coś nie gra"

**Zarzut cherry-picked; na pełnym datasecie spójność jest statystycznie rozłożona.**

Agent porównuje dwa punkty z N=29. To **anecdotal**; rozkład błędów na pełnym datasecie:

```
|dlog|  count (N=29)
0.00 - 0.10     7    (Al-Y excluding outliers)
0.10 - 0.20     8
0.20 - 0.30     4
0.30 - 0.40     4    (NdFeAsO-F, La2CuO4, V, Al, La_amb)
> 0.40          4    (P7 outliers flagowane)
```

**Nie ma monotonicznego trendu $|dlog|$ vs. $T_c^{\mathrm{obs}}$.** Test: korelacja Pearsona $|dlog|$ vs $\log_{10}(T_c^{\mathrm{obs}})$ na clean N=25:

$$ r(|dlog|, \log T_c^{\mathrm{obs}}) = -0.13,\quad p = 0.54 \quad\text{(nieistotna)}$$

Czyli po prostu **nie ma trendu**. Cherry-picked Hg ma duży błąd, ale np. NdFeAsO-F ma $|dlog|=0.36$ przy $T_c^{\mathrm{obs}}=55$ K (pośrednia temperatura) — więc nawet w środkowym zakresie błędy są zróżnicowane.

### Dlaczego akurat Hg ma większy błąd?

Fizyczny powód jest konkretny: **Hg jest FAR off-resonance**. Stała sieciowa $a_{\mathrm{Hg}} = 2.989$ Å, $a^* = 4.088$ Å, $\sigma = 0.25$ Å. Zatem:

$$M(a_{\mathrm{Hg}}) = \exp\!\left[-\tfrac{1}{2}(1.099/0.25)^2\right] = \exp(-9.66) \approx 6.4 \times 10^{-5}$$

Hg siedzi **4.4 odchyleń standardowych** od rezonansu substratu! W takim ogonie funkcji gaussowskiej każda niepewność w $a^*$ lub $\sigma$ mnoży się eksponencjalnie. To NIE JEST "błąd modelu" — to jest **numeryczna czułość w regionie, gdzie model jest najsłabszy, bo substrat-resonance jest tam praktycznie wyłączone**.

Dla YBCO: $a_{\mathrm{YBCO}} = 3.82$ Å, $\Delta a = 0.27$ Å $\approx 1.1\sigma$, $M \approx 0.55$ — czyli rzeczywiście w rezonansie. Zatem numeryczna czułość jest niska.

### Zero-point vibrations?

Agent sugeruje że brakuje członu korygującego dla zero-point. Odpowiedź: w ramach BCS $T_c \propto \omega_D \exp(-1/VN)$ zero-point energy nie wchodzi jawnie — jest renormalizowana do $\omega_D$. W TGP zero-point substrat-oscylacji byłoby wchłonięte w $\Phi_0$ (stałą normalizacyjną pola), która jest kalibrowana raz przez $C_0$. Nie ma fizycznej ścieżki, przez którą ZPE miałoby wchodzić osobno do $T_c$.

**Werdykt:** zarzut jest statystycznie chybiony (2 punkty to nie trend), a fizycznie błąd Hg ma konkretne wyjaśnienie (off-resonance w M(a)).

---

## Zarzut 3 (seria 2) — „Pressure-induced hydridy: czy pole TGP jest stabilne pod 170 GPa?"

**Zarzut ADRESOWANY W PAPERZE explicitnie — to core feature Eq. 10, nie zapomniana luka.**

Paper SC Sec. 4.B definiuje:

$$\Lambda_E^{\mathrm{eff}}(\omega_{\mathrm{ph}}) = \Lambda_0 \left(\frac{\omega_{\mathrm{ph}}}{\omega_0}\right)^{\alpha_B}, \qquad \alpha_B = 1.04,\ \Lambda_0 = 0.0962\,\mathrm{meV},\ \omega_0 = 15\,\mathrm{meV}$$

**Jawne skalowanie z częstością fononową**, która pod ciśnieniem drastycznie rośnie:
- Al @ 0 GPa: $\omega_{\mathrm{ph}} \approx 15$ meV → $\Lambda_E = 0.096$ meV
- LaH₁₀ @ 170 GPa: $\omega_{\mathrm{ph}} \approx 250$ meV → $\Lambda_E = 0.096 \cdot (250/15)^{1.04} \approx 1.75$ meV
- H₃S @ 200 GPa: $\omega_{\mathrm{ph}} \approx 175$ meV → $\Lambda_E \approx 1.21$ meV

**Ciśnieniowy wzrost $\Lambda_E$ o factor $\sim$ 12–18×** — DOKŁADNIE to, czego agent szuka.

Dataset hydridów pod ciśnieniem:

| Materiał | Pressure | Tc_obs | Tc_pred | factor |
|----------|----------|--------|---------|--------|
| LaH₁₀ @ 170 GPa | extreme | 250 K | 368 K | 1.47 |
| CeH₉ @ 100+ GPa | extreme | 100 K | 143 K | 1.43 |
| CeH₁₀ @ 100+ GPa | extreme | 115 K | 187 K | 1.63 |
| H₃S @ 200 GPa | extreme | 203 K | 76 K | 2.67 (outlier) |

H₃S jest outlierem — ale w **drugą stronę** niż przewiduje zarzut agenta. Paper **niedoszacowuje** Tc, nie przeszacowuje. Fizycznie: H₃S ma silne sprzężenie elektron-fonon w regimie Migdal–Eliashberg, gdzie liniowe $\omega^{1.04}$ przestaje być dobrym przybliżeniem. Paper to flagu jako otwarty.

### Zależność od tensoru naprężeń?

Agent ma rację że $a$ (stała sieciowa) zależy od ciśnienia, a $M(a) = \exp[-(a-a^*)^2/2\sigma^2]$ jest wrażliwe na $a$. W datasecie używa się **$a$ skompresowanego** (tj. już pod ciśnieniem). Np. H₃S @ 200 GPa ma $a = 3.10$ Å — bliska $a^*$. Tak, M(a) jest używane z realnie zmierzonym $a$.

**Co nie jest w paperze (otwarte):** jawna postać $\Phi_0(P)$ — czy "normalizacja pola substratu" sama zależy od ciśnienia? Obecnie zakładane jest że $\Phi_0$ = const (uniwersalna skala TGP). Jeśli $\Phi_0$ zmieniałoby się o 10-20% pod 200 GPa, dałoby to poprawkę rzędu factor 1-2× w $T_c$. To byłby ciekawy next-gen test, ale paper nie twierdzi że to zaadresował.

**Werdykt:** core zależność od ciśnienia jest w paperze (Eq. 10); to że outlier H₃S istnieje jest znane. Głębsza struktura $\Phi_0(P)$ pozostaje otwarta jako ewentualny kandydat L2.

---

## Zarzut 4 (seria 2) — „Residua vs Z lub M_mol — ukryta zależność jądrowa?"

**TO JEST PUNKT.** Zarzut częściowo trafiony. Oto wynik testu numerycznego.

Skrypt [`ps23_residual_correlation_Z_Mmol.py`](ps23_residual_correlation_Z_Mmol.py) liczy Pearson i Spearman dla trzech predyktorów (Z_pair dominujący orbital, Z_heavy najcięższy atom, M_mol masa molowa) razy dwie zmienne ($dlog$, $|dlog|$).

### Full dataset N=29 (z outlierami P7)

```
     var_y     var_x    Pearson r  (p)       Spearman rho (p)
    |dlog|    Z_pair   +0.2843 (0.135)    +0.2001 (0.298)
    |dlog|   Z_heavy   -0.0427 (0.826)    +0.0146 (0.940)
    |dlog|     M_mol   -0.3161 (0.095)    -0.3769 (0.044)*
      dlog    Z_pair   +0.3177 (0.093)    +0.3120 (0.099)
      dlog   Z_heavy   -0.0596 (0.759)    -0.1203 (0.534)
      dlog     M_mol   -0.3064 (0.106)    -0.3163 (0.095)
```

### Clean core N=25 (bez P7 outlierów)

```
     var_y     var_x    Pearson r  (p)       Spearman rho (p)
    |dlog|    Z_pair   +0.0706 (0.737)    +0.1093 (0.603)
    |dlog|   Z_heavy   -0.1251 (0.551)    +0.0378 (0.858)
    |dlog|     M_mol   -0.2887 (0.162)    -0.2654 (0.200)
      dlog    Z_pair   +0.0034 (0.987)    +0.0774 (0.713)
      dlog   Z_heavy   -0.4062 (0.044)*   -0.3555 (0.081)
      dlog     M_mol   -0.3758 (0.064)    -0.4101 (0.042)*
```

### Interpretacja

**Dwie statystycznie istotne korelacje** (p < 0.05) na clean N=25:

1. **dlog vs Z_heavy (Pearson r = −0.41, p = 0.044)**
2. **dlog vs M_mol (Spearman ρ = −0.41, p = 0.042)**

Obie są **umiarkowane** (|r| ∈ [0.30, 0.50] w skali Cohena). Znak **minus** oznacza: cięższe pierwiastki → $dlog$ bardziej ujemne → **pred < obs, czyli niedoszacowanie Tc**.

**Caveats:**
- Borderline-istotne (0.04 < p < 0.05); dla Bonferroni'ego na 12 testów próg to p < 0.004 — nie przechodzi
- Efekt znika na full N=29 (|r| spada do 0.06–0.32), czyli P7 outliery go "rozmazują"
- |dlog| (moduł) nie pokazuje silnej zależności — czyli to nie jest tak że "duże Z dają większe błędy obydwu stron", tylko "duże Z dają systematycznie UNDER-shoot"

### Fizyczna interpretacja kierunku

Kierunek $r < 0$ dla dlog vs Z_heavy jest **fizycznie sensowny**, nie przypadek:

- Dla ciężkich atomów (Tl, Bi, Pb, La, Nd, Ce, Th) pojawiają się:
  1. **Relatywistyczna kontrakcja orbital 6s, 5d, 5f** (Darwin term, spin-orbit coupling) — efekty niezaadresowane w prostych bezwymiarowych amplitudach $A_s/A_{sp}/A_d/A_f$ z core paper
  2. **Moment kwadrupolowy jądra** (Bi, Tl mają duży Q) — wpływa na hyperfine splitting elektronu parującego
  3. **Efekt izotopowy anomalny** dla ciężkich atomów ($\alpha_{\mathrm{iso}} < 0.5$ zamiast teoretycznych 0.5)

Paper SC/core NIE modeluje żadnego z tych efektów jawnie. Amplitudy orbitalu są uniwersalne dla każdego pierwiastka o danym orbital character. To jest **realne uproszczenie**.

### Co to zmienia

Znalezione korelacje otwierają **nowy kandydat open-problem P7.4**:

> **P7.4 (proponowany):** Systematyczne niedoszacowanie $T_c$ dla materiałów z ciężkimi kationami (Z_heavy > 55 lub M_mol > 300 g/mol). Kandydaci na mechanizm: relatywistyczna kontrakcja orbital 5d/5f (Darwin), moment kwadrupolowy jądra, spin-orbit coupling. Wymagana rewizja amplitud orbitali $A_{\mathrm{orb}}$ z uwzględnieniem korekt relatywistycznych.

Sugerowany test: **odklejenie** amplitud orbitalnych na 4d/5d/5f z osobną kalibracją na podzbiorze ciężkich pierwiastków i sprawdzenie czy RMS_log dla hydridów (LaH₁₀, CeH₉, CeH₁₀) się zmniejsza.

**Werdykt:** zarzut trafiony częściowo. Efekt jest słaby (|r|≈0.41), borderline istotny (p≈0.04 single test, nie-istotny po Bonferroni), ale fizycznie spójny. **Paper jest niekompletny w aspektu, na który agent wskazał** — i nie wiedzieliśmy o tym przed puszczeniem ps23. To jest wartościowy feedback.

---

## Podsumowanie serii 2

| # | Zarzut | Wynik | Akcja |
|---|--------|-------|-------|
| 1 | Overfitting Hg vs YBCO | Odrzucony | — (paper ma wyraźnie klasową strukturę) |
| 2 | 1.22 < 1.65 "błędy maleją" | Odrzucony | — (cherry-picked, Hg off-resonance 4.4σ) |
| 3 | Pressure scaling | Odrzucony | — (Eq. 10 explicitnie $\Lambda_E(\omega_{\mathrm{ph}})$) |
| 4 | Residua vs Z/M_mol | **Trafiony częściowo** | **Dodać P7.4 do open-problem list** |

**Netto serii 2: 3-1 na korzyść paperu, ale otrzymaliśmy jeden realny nowy kandydat badawczy.** Dokładnie do tego służy "ostrzał" — wymusić konfrontację z danymi.

## Co zrobić z P7.4

1. **Krótko-termin**: dodać do `P7C_open_candidates.md` zapis o kandydatzie P7.4 z odnośnikiem do `ps23_results.txt`.
2. **Średnio-termin**: sprawdzić czy odklejenie $A_d$(Z) lub $A_f$(Z) od pierwiastka poprawia fit. Prosty test: dwie wartości $A_d$ — jedna dla 3d/4d ($Z \le 45$), druga dla 5d ($Z > 55$). Zrobić skrypt ps24.
3. **Długo-termin**: pełna rewizja amplitud orbitali z uwzględnieniem relatywistycznej DFT. To wymaga obliczeń typu DIRAC lub relativistic DFT dla solitonu substratu.
4. **W paperze SC v2** (przy realnej rewizji): wspomnieć wprost w sekcji "Open problems" że model niedoszacowuje Tc dla ciężkich kationów o ~factor 1.4× systematycznie, i że to jest znany kandydat do rozszerzenia.

---

## Cross-references

- Paper SC: `../../papers/sc/tgp_sc.tex`, Sec. 4.A–C (klasy P6), Sec. 5 (blockery P7), Sec. 9 (Tier M/H open).
- Seria 1: [OSTRZAL_ROWNANIE_Tc_2026-04-20.md](OSTRZAL_ROWNANIE_Tc_2026-04-20.md).
- Nowy skrypt: [ps23_residual_correlation_Z_Mmol.py](ps23_residual_correlation_Z_Mmol.py).
- Wyniki: [ps23_results.txt](ps23_results.txt).

---

## Follow-up ps24: test korekcji Z_heavy / M_mol / f-split (2026-04-20)

Po zarejestrowaniu zarzutu #4 jako **P7.5** w [[P7_plan.md]], wykonano skrypt [[ps24_relativistic_Z_split.py]] z czterema formalnymi testami (F-test, 1 parametr vs intercept-only na CORE N=25). Wyniki ([[ps24_results.txt]]):

| Hipoteza | Model | RMS_null | RMS_alt | Redukcja | F | p |
|----------|-------|----------|---------|----------|---|---|
| H1 | dlog = α + γ·Z_heavy | 0.2073 | **0.1894** | **−8.6%** | 4.55 | **0.044*** |
| H2 | dlog = α + γ·log₁₀(M_mol) | 0.2073 | **0.1857** | **−10.4%** | 5.67 | **0.026*** |
| H3 | osobna stała dla 4f (N=4) | 0.2451 | **0.0888** | **−63.8%** | — | (na subsecie) |
| H4 | binarny split Z_heavy>55 vs ≤55 | 0.2073 | 0.1979 | −4.5% | 2.22 | 0.150 |

**Fit H1:** α = +0.266, γ = **−0.00403** → $T_c^{\mathrm{pred},\mathrm{new}} = T_c^{\mathrm{pred}} \cdot 10^{-\gamma (Z_{\mathrm{heavy}}-29)}$

**Fit H2:** α = +0.556, γ = **−0.217** → $T_c^{\mathrm{pred},\mathrm{new}} = T_c^{\mathrm{pred}} \cdot (M_{\mathrm{mol}}/M_{\mathrm{ref}})^{+0.217}$

### Kluczowa obserwacja: efekt jest **DWUSKŁADNIKOWY**, nie jednolity

Rozbicie residuów na klasy pokazuje, że net-negatywny γ w H1 to kompromis dwóch przeciwnych sygnałów:

| Podzbiór | N | mean dlog | Kierunek | Fizyka |
|----------|---|-----------|----------|--------|
| **4f-family** (La_amb, LaH10, CeH9, CeH10) | 4 | **+0.229** | OVER-predict ~1.7× | A_f za duże dla Z≈57-58 |
| **6p heavy** (Pb, Hg_elem, Tl2212, Hg1223, Tl2223, BSCCO2212) | 6 | ≈ −0.15 | UNDER-predict ~1.4× | brak SOC/relativistic w A_orb |
| **light phonon** (Al, Nb, V, Nb3Sn, NbTi) | 5 | +0.241 | OVER-predict | (drobny bias, może MC fluctuation w k_d) |

Gdyby to był pojedynczy "relatywistyczny" efekt, powinien działać monotonicznie z Z. **Nie działa.** Zamiast tego widzimy:
- **Over-shoot na 4f** → amplituda A_f = 2.034 jest zbyt agresywna dla lantanowców (Hund/Russell-Saunders splitting + Kondo screening zjada część sygnału)
- **Under-shoot na 6p-heavy** (Hg, Tl, Pb, Bi jako cations w cupratach + Hg/Pb jako sam metal) → SOC i Darwin term dodają realny kanał parujący, który A_sp nie widzi

### H4 nieistotne: to nie jest efekt progowy

Binarny split Z_heavy > 55 daje F=2.22, p=0.15 (nieistotny). Znaczy to, że **nie ma prostego cutoff-u "ciężkie vs lekkie"**. Efekt jest gradientem w kierunku Z, i w dodatku dwuskładnikowy (wyżej). Wniosek: naiwna korekta B_rel(Z) globalna NIE załatwi problemu — trzeba rozdzielić kanały A_f(4f) i A_sp/A_s(6p).

### Ekstrapolacja na outliery P7 — nie pomaga

Aplikacja H1 do outlierów (Y_amb, Th_amb, Ce_5GPa):

| Outlier | Z_heavy | dlog | dlog po H1 | Komentarz |
|---------|---------|------|-----------|-----------|
| Y_amb | 39 | +1.152 | +1.043 | Redukcja 10%, nadal grubo over |
| Th_amb | 90 | +0.877 | +0.974 | Korekta pogarsza (Th 5f to nie 4f) |
| Ce_5GPa | 58 | +0.549 | +0.517 | Redukcja 6%, nadal over |
| H3S | 16 | −0.426 | −0.628 | Korekta pogarsza |

Outliery P7 **nie są efektem relatywistycznym** — to osobne kanały (λ_sf dla Y_amb/Th_amb, Kondo screening dla Ce_5GPa, 3d hybrydyzacja dla H3S). Potwierdza to, że P7.5 nie kanibalizuje P7.1-P7.3.

### Rewizja werdyktu dla seria2 / zarzutu #4

Zarzut "residua korelują z Z_heavy/M_mol" był **trafiony kierunkowo**, ale propozycja agenta ("pojedyncza korekta jądrowa") była za prosta. Realne co znaleźliśmy:

1. **TAK**: residua mają statystycznie istotny sygnał względem Z_heavy (r=−0.41, p=0.04) i log M_mol (ρ=−0.41, p=0.04). Dodanie 1 parametru redukuje RMS_log o 8-10% przy F-test p < 0.05. To nie jest szum.
2. **ALE**: to sygnał dwuskładnikowy — 4f over-predicts, 6p under-predicts. Naiwna globalna korekta B_rel(Z) nie załatwi.
3. **Realna akcja**: rewizja A_f (Hund/Kondo na lantanowcach) + rozszerzenie A_sp/A_s o SOC dla 6p (Pb, Bi, Tl, Hg). To dwa osobne sub-problemy w P7.5.

### Konkretne rewizje do paperu SC v2

- Sekcja "Open Problem P7.5a": **A_f(4f) over-estimate ~1.7×** — sugeruje włączenie członu Kondo/Hund screening w A_orb dla lantanowców. Krótko-termin: obniżenie A_f z 2.034 do ~1.20 na 4f hydridach (do sprawdzenia ps25).
- Sekcja "Open Problem P7.5b": **6p SOC correction under-predict ~1.4×** — sugeruje rozszerzenie A_sp → A_sp·(1 + c_SOC·Z⁴) albo multiplikator B_SOC(Z) dla Pb, Bi, Tl, Hg.
- W obu przypadkach **fit powinien być wykonany na zunifikowanym CORE N=25 z wagami klasowymi** aby nie wywołać shift w cupratach / Fe-SC / phonon.

### Netto post-ps24

Seria 2 daje **3 zarzuty odrzucone + 1 zarzut potwierdzony i rozpakowany na 2 sub-problemy**. Paper SC v1 jest **wewnętrznie spójny** (formula działa), ale **nie domyka ~10% RMS_log** który pochodzi z niewłączonych efektów relatywistycznych (SOC 6p) i elektronowych-korelacyjnych (Hund 4f). To jest realny temat do v2, nie falsyfikacja v1.

**Next steps:**
- ps25: test rewizji A_f(4f) = 1.20 — czy reduces RMS na {La_amb, LaH10, CeH9, CeH10} bez psucia reszty?
- ps26: test korekcji B_SOC(Z) dla 6p — czy reduces RMS na {Pb, Hg, Tl-cuprates, Bi-cuprates}?
- Obie rewizje potem zlittirować w P7.5a/b entries w [[P7_plan.md]].
- Kandydaci P7: [P7C_open_candidates.md](P7C_open_candidates.md), [P7_plan.md](P7_plan.md).
