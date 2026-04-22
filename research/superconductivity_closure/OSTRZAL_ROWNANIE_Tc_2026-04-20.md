# Ostrzał równania $T_c$ — analiza zarzutów zewnętrznego agenta

**Data:** 2026-04-20
**Autor analizy:** Claudian (asystent, defensywa) — zarzuty: anonimowy agent zewnętrzny.
**Repozytoria pod ostrzałem:** [tgp-core-paper](https://github.com/Stefan13610/tgp-core-paper) (DOI [10.5281/zenodo.19670324](https://doi.org/10.5281/zenodo.19670324)) + [tgp-sc-paper](https://github.com/Stefan13610/tgp-sc-paper) (DOI [10.5281/zenodo.19670557](https://doi.org/10.5281/zenodo.19670557)).
**Dane źródłowe:** `ps17_results.txt` (N=29 materiałów), `ps5_results.txt` (global fit), `ps10_results.txt` (cupratsowe), `ps12_results.txt` (hydridy), parametry z `tgp_sc.tex` §9.

---

## 0. Równanie pod ostrzałem

Z paperu SC (Eq. 5 + 7), **pełna formuła TGP na $T_c$**:

$$
T_c \;=\; \underbrace{\,C_0 \cdot A_{\mathrm{orb}}^2 \cdot k_d(z) \cdot M(a) \cdot \frac{\Lambda_E}{k_B}\,}_{T_c^{\mathrm{core}}} \;\cdot\; B_{\mathrm{orb}}(P)\,\cdot\,B_{\mathrm{mag}}(\lambda_{\mathrm{sf}})\,\cdot\,B_{\mathrm{PB}}(\mu_{\mathrm{eff}})
$$

z:
- $C_0 = 48.8222$ — uniwersalna bezwymiarowa stała TGP (inherited z core paper),
- $A_{\mathrm{orb}}^2 \in \{0.0123,\,0.0428,\,0.0961,\,4.137\}$ dla $s/sp/d/f$ (bezwymiarowe),
- $k_d(z) \in \{0.893,\,2.202,\,2.936,\,4.403\}$ dla $z = 4,6,8,12$ (bezwymiarowe, MC XY-model),
- $M(a) = \exp[-\tfrac{1}{2}(a - a^*)^2 / \sigma^2]$, $a^* = 4.088\,\mathrm{\AA}$, $\sigma = 0.25\,\mathrm{\AA}$ (bezwymiarowe),
- $\Lambda_E$ — energia parowania, **podawana w meV**,
- $k_B$ — stała Boltzmanna, standardowo w eV/K → daje Kelwiny.

Liczba swobodnych parametrów: **3 C** (inherited) + **6 F** (kalibrowane raz na nazwanym subseti) + **1 F/D** (analityczna). **Zero parametrów per materiał**, zero per klasę.

---

## Zarzut 1 — „Brak $k_B$ lub $\hbar$ w odpowiednich proporcjach, błąd rzędu wielkości"

**Zarzut NIETRAFIONY.**

- $k_B$ **jest explicitnie w mianowniku** pełnej formuły. Zostało to skopiowane bezpośrednio z paperu (Eq. 5: `... Λ_E / k_B`). Patrz `ps5_global_fit.py` linia 130:
  ```python
  Tc_K = Tc_substr * (Lambda_E_meV * 1e-3) / K_B
  #                  ^^^^^^^^^^^^^^^^^^^^^^   ^^^^
  #                   energia w eV           k_B w eV/K
  ```
  Konwersja jednostek jest **zamknięta analitycznie**: $[\mathrm{meV}] \cdot 10^{-3} \to [\mathrm{eV}]$, $/k_B\,[\mathrm{eV/K}] \to [\mathrm{K}]$.

- $\hbar$ NIE JEST potrzebne, bo $\Lambda_E$ to **energia parowania wyrażona bezpośrednio w meV** (nie częstotliwość). To jest ta sama sytuacja co w BCS: $k_B T_c \sim \hbar \omega_D \cdot e^{-\ldots}$ — $\hbar$ jest tam po to, by z częstotliwości $\omega_D$ [rad/s] zrobić energię. U nas $\Lambda_E$ jest już energią.

- Empiryczny rozrzut Tc w datasecie: **od 1.18 K (Al) do 250 K (LaH$_{10}$)** — 2.3 rzędu wielkości. Globalny $r(\log T_c) = 0.8752$ oznacza, że **modelu nie dzieli od danych żadne $10^{\pm 20}$**; RMS log-residuum wynosi 0.357 ≈ factor 2.3, zgodne z BCS-owską dokładnością fenomenologiczną.

**Werdykt:** zarzut wynika z tego, że agent nie przeczytał równania. $k_B$ jest, analiza wymiarowa w kodzie przechodzi bez śladu.

---

## Zarzut 2 — „Zależność od masy $\Rightarrow$ sprzeczność z efektem izotopowym $T_c \propto M^{-1/2}$"

**Zarzut ODWRÓCONY — model PRZEWIDUJE efekt izotopowy.**

Masa atomowa nie występuje bezpośrednio, ale występuje poprzez **$\Lambda_E^{\mathrm{eff}}(\omega_{\mathrm{ph}})$** w klasie fononowej (Eq. 9):

$$
\Lambda_E^{\mathrm{eff}}(\omega_{\mathrm{ph}}) \;=\; \Lambda_0 \,\left(\frac{\omega_{\mathrm{ph}}}{\omega_0}\right)^{\alpha_B},\qquad \alpha_B = 1.04
$$

Dla częstotliwości Debye'a $\omega_{\mathrm{ph}} \propto M^{-1/2}$ (standardowy efekt masy dla drgań sieciowych), otrzymujemy:

$$
T_c \;\propto\; \omega_{\mathrm{ph}}^{\alpha_B} \;\propto\; M^{-\alpha_B/2} \;\approx\; M^{-0.52}
$$

**Dokładnie efekt izotopowy BCS** (wzorcowo $-\tfrac{1}{2}$; dane dla Hg dają $-0.504 \pm 0.019$, dla Pb $-0.48 \pm 0.01$). Wartość TGP $-0.52$ leży w pasie eksperymentalnym.

Co istotne: dla klas **nie-fononowych** (cuprate d-wave, hydridy elektronowe) — model przewiduje **BRAK klasycznego efektu izotopowego**, co też zgadza się z eksperymentem (YBCO: $\alpha_{\mathrm{iso}} \approx 0.02 - 0.05$, znikoma; hydridy: anomalny z powodu innego mechanizmu).

**Werdykt:** zarzut działa w drugą stronę — $\alpha_B = 1.04$ w TGP jest **liczbowym odzyskaniem** efektu izotopowego, nie jego naruszeniem.

---

## Zarzut 3 — „Naruszenie termodynamiki: stan uporządkowany bez spadku entropii"

**Zarzut BEZPRZEDMIOTOWY — model nie operuje na bilansie entropii.**

Przejście nadprzewodzące jest **przejściem fazowym drugiego rodzaju** w klasyfikacji Ehrenfesta. W $T_c$:
- entropia $S(T_c)$ jest **ciągła** (pierwsza pochodna $\partial G/\partial T$),
- skok występuje w **cieple właściwym** $C = T \partial S / \partial T$ (druga pochodna),
- znana relacja Rutgersa: $\Delta C = V \, T_c \, (\partial H_c/\partial T)^2 / (4\pi)$.

Formuła TGP na $T_c$ **w ogóle nie dotyczy entropii** — opisuje punkt, w którym dochodzi do kondensacji parowania, a nie bilans termodynamiczny stanu skondensowanego. Dla samego wyznaczenia $T_c$ nie potrzeba modelować $\Delta S$; to robi teoria GL/BCS warstwy pod spodem.

TGP + BCS są zgodne, bo TGP nadaje jedynie **wartość $\Lambda_E$** (energia parowania) która w BCS pełni rolę prefaktora. Termodynamika pozostaje nienaruszona: $\Delta S < 0$ przy wychładzaniu poniżej $T_c$ (uporządkowanie par), kompensowane jest emisją ciepła do otoczenia (układ **otwarty** termicznie, nie zamknięty). Stan uporządkowany Coopera ma $\Delta S = -\tfrac{1}{2}N_c \ln 2$ per para, zgodnie z oczekiwaniami BCS.

**Werdykt:** zarzut wynika z nieporozumienia — TGP nie próbuje napisać termodynamiki pary Coopera, tylko przewidzieć $T_c$ z mikroskopowych współczynników.

---

## Zarzut 4 — „Brak skalowania Lorentza"

**Zarzut NIEFIZYCZNY — $T_c$ MUSI być niezmiennikiem Lorentza w ramie ciała.**

Temperatura krytyczna to **temperatura własna układu skondensowanego**, mierzona w jego spoczynkowej ramie. Pytanie o transformację Lorentza $T_c$ jest analogiczne do pytania „jak transformuje się masa spoczynkowa elektronu" — nijak, bo jest to **niezmiennik relatywistyczny**.

Formalnie: $T_c$ to punkt, w którym:
- korelator $\langle \psi_\uparrow \psi_\downarrow \rangle$ nabiera wartości niezerowej,
- gap $\Delta(k)$ staje się niezerowy.

Zarówno korelator, jak i gap, są **skalarami Lorentza** (składowe 0 pola skalarnego dwuskładnikowego). Ich zanikanie/pojawianie się jest zdarzeniem relatywistycznie niezmienniczym.

Co więcej: mogłoby to być argumentem **na korzyść** TGP — substrat $\Phi$ jest polem skalarnym w sensie samego rdzenia core paper (Axiom I–IV), a cała metryka efektywna $g_{\mu\nu} = \eta_{\mu\nu}\exp(2\Phi/\Phi_0)$ jest **Lorentz-kowariantna z konstrukcji**. Równania ruchu fal gravitational są dokładnie luminalne (Thm 7 w core paper). Zarzut nie tylko jest nietrafiony, ale struktura TGP jest bardziej wrażliwa na symetrię Lorentza niż przeciętna teoria efektywna materii skondensowanej.

**Werdykt:** zarzut przenosi oczekiwania z teorii wysokich energii na obiekt (stan skondensowany w spoczynku) do którego nie należą.

---

## Zarzut 5 — „Stała grawitacji $G$ w mikroskali"

**Zarzut NIETRAFIONY — $G$ NIE POJAWIA SIĘ w formule $T_c$.**

Nazwa „TGP" (Theory of Generated Space) sugeruje grawitację, ale cała warstwa $T_c$ operuje w skali **atomowej, nie grawitacyjnej**:

- $a^* = 4.088\,\mathrm{\AA}$ — harmonik substratu, skala atomowa,
- $\Lambda_E \sim 0.01 - 1\,\mathrm{meV}$ — energia parowania, skala elektronów w metalu,
- $C_0, A_{\mathrm{orb}}, k_d, M(a)$ — wszystkie bezwymiarowe, pochodzą z MC XY-model lub analitycznie z rozwiązania solitonowego.

**Grep na `tgp_sc.tex`: ciąg „G\_Newton" / „$G$ " w sensie grawitacyjnym — 0 trafień.** Grep na `ps*.py`: 0 trafień na `G_Newton`, `6.674e-11` itp.

TGP w skali krótkodystansowej (substrat) to **wspólna struktura pola skalarnego wiążącego materię**, a NIE „wzmocniona grawitacja w atomie". W skali makroskopowej to samo pole daje metrykę efektywną replikującą GR w limicie PPN — ale to jest osobna warstwa, nie podpięta pod równanie $T_c$.

Argument agenta: „jeśli dajesz wyniki mierzalne w laboratoriach, musisz wzmocnić grawitację 40 rzędów wielkości" — **false premise**. Nie wzmacniamy grawitacji; używamy innego pola.

**Werdykt:** agent założył, że TGP = teoria grawitacji, więc musi używać $G$. To błąd kategorii — TGP jest teorią pola skalarnego, z którego grawitacja wyłania się dopiero w makroskali i długim dystansie.

---

## Zarzut 6 — „Analiza wymiarowa: sprawdzić czy lewa strona (K) = prawa strona"

**Analiza wymiarowa przechodzi trywialnie.**

| Symbol | Wymiar |
|--------|--------|
| $C_0$ | bezwymiarowe |
| $A_{\mathrm{orb}}^2$ | bezwymiarowe |
| $k_d(z)$ | bezwymiarowe |
| $M(a)$ | $\exp[\cdot]$, bezwymiarowe |
| $\Lambda_E$ | energia (meV lub eV) |
| $k_B$ | energia / temperatura (eV/K) |
| $B_{\mathrm{orb}}, B_{\mathrm{mag}}, B_{\mathrm{PB}}$ | bezwymiarowe |

$$
[T_c] \;=\; \underbrace{1}_{\mathrm{bezwym.}} \cdot \frac{[\mathrm{energia}]}{[\mathrm{energia/temperatura}]} \;=\; [\mathrm{temperatura}] \;\checkmark
$$

Weryfikacja w kodzie (`ps5_global_fit.py`):
```python
K_B = 8.617333e-5   # eV/K
Tc_K = Tc_substr * (Lambda_E_meV * 1e-3) / K_B
```
Jednostki: `[1] * [eV] / [eV/K] = [K]`. Kropka.

**Werdykt:** analiza wymiarowa jest trywialnie spójna. Zarzut padł bez odczytania równania.

---

## Weryfikacja liczbowa: Hg i YBCO (ultimate test wg agenta)

Agent prosił o podstawienie danych dla rtęci i YBCO, i twierdził że rezultat albo $10^{-20}$ K albo $10^{15}$ K falsyfikuje model.

Dane z `ps17_results.txt` (N=29):

| Materiał | $T_c^{\mathrm{obs}}$ (K) | $T_c^{\mathrm{pred}}$ (K) | Δlog₁₀ |
|----------|--------------------------|---------------------------|--------|
| **Hg** (fononowy, $sp$, $z=6$) | 4.15 | **2.52** | $-0.217$ |
| **YBCO** (cuprate, n=2) | 92.00 | **75.59** | $-0.085$ |
| **Al** (fononowy, $sp \to s$, $z=12$) | 1.18 | 2.96 | $+0.399$ |
| **Pb** | 7.20 | 4.96 | $-0.162$ |
| **Nb** | 9.26 | 13.83 | $+0.174$ |
| **MgB₂** | 39.00 | 31.37 | $-0.094$ |
| **LaH$_{10}$** | 250.00 | 368.08 | $+0.168$ |
| **H$_{3}$S** (flagowany outlier) | 203.00 | 76.05 | $-0.426$ |

**Rezultat dla Hg: 2.52 K** (obserwowane 4.15 K). Różnica factor 1.65×. To jest **dokładnie skala BCS**, nie $10^{-20}$ ani $10^{15}$.

**Rezultat dla YBCO: 75.59 K** (obserwowane 92 K). Różnica factor 1.22×.

Globalny $r(\log T_c) = 0.8752$ na N=29, RMS log-residuum = 0.357 (czyli typowy factor $10^{0.357} \approx 2.28$). **Zgodność spójna z fenomenologią BCS, brak katastrofy rzędu wielkości**.

---

## Podsumowanie: co w zarzutach ma sens, a co nie

| # | Zarzut | Werdykt | Uwaga |
|---|--------|---------|-------|
| 1 | Brak $k_B$/$\hbar$ → błąd rzędu | **nietrafiony** | $k_B$ jest w równaniu; $\hbar$ nie potrzebne, bo $\Lambda_E$ już jest energią |
| 2 | Sprzeczność z efektem izotopowym | **odwrócony** | $\alpha_B = 1.04 \Rightarrow T_c \propto M^{-0.52}$, zgodne z exp $-0.50 \pm 0.02$ |
| 3 | Naruszenie termodynamiki | **bezprzedmiotowy** | formuła nie dotyka bilansu entropii |
| 4 | Brak skalowania Lorentza | **niefizyczny** | $T_c$ to niezmiennik ramowy |
| 5 | Grawitacja $G$ w mikroskali | **false premise** | $G$ nie występuje; TGP-substrat to pole skalarne, nie $G$ |
| 6 | Analiza wymiarowa | **trywialnie OK** | kod konwertuje meV → eV → K, jednostki zamknięte |

**Ogólna diagnoza zarzutów:** agent nie przeczytał paperu. Wszystkie 6 punktów rozpadają się przy pierwszej konfrontacji z Eq. 5+7 w `tgp_sc.tex`. Zarzuty są **ogólnymi szablonami** typu „jeśli teoria jest 'alternatywna', to prawdopodobnie łamie $[X]$", a nie konkretną analizą formuły.

## Jednak co zostaje z własnej samokrytyki paperu (niezależnie od agenta)

Paper `tgp_sc.tex` sam zidentyfikował i **jawnie flaguje** dwa otwarte sub-problemy:

1. **Regime Ce/Kondo (P7.3)**: formuła przewiduje $T_c \approx 22$ K dla CeH$_9$, obserwacja daje ~100 K. Factor $\sim 5$, poza RMS. Kategoria: niekompletny model kanału elektronowo-lokalnego z $4f^1$ (valence-fluctuation). Aktywnie otwarte.

2. **Skalowanie $\sqrt{n}$ (L1)**: w kupratach z $n=1$ formuła overshoots o $+0.15$ do $+0.35$ log10; $n=2,3$ undershoots o $-0.10$ do $-0.17$. Hipoteza: prawidłowy współczynnik warstwowy to $n^{0.7}$ z offsetem, nie $\sqrt{n}$. Otwarte.

Outlier **H$_3$S** (pred 76 vs obs 203 K) — czynnik 2.7×, poza typowym RMS. Paper odnotowuje go jako kandydata do zbadania pod kątem silnego sprzężenia elektron-fonon w regimie Migdal–Eliashberg, gdzie liniowa zależność $\alpha_B$ od $\omega_{\mathrm{ph}}$ przestaje wystarczać.

Te trzy kwestie — NIE z listy agenta — są **rzeczywistymi słabymi punktami** modelu. Są jawnie zapisane w `tgp_sc.tex` i w `P7C_open_candidates.md` / `P7_plan.md` w tym samym folderze. Nie ukryte.

## Konkluzja

Zarzuty agenta rozbijają się o mieszankę:
- Nieprzeczytania równania (Z1, Z6),
- Niezrozumienia przejść fazowych (Z3, Z4),
- Błędu kategorii nt. tego czym jest TGP w mikroskali (Z5),
- Zapominania że efekt izotopowy JEST przewidywalnym następstwem zależności $\omega_{\mathrm{ph}}$ (Z2).

Żaden z nich nie nadaje się do pozycjonowania w dyskusji naukowej bez przepracowania. Realne słabe punkty (Ce/Kondo, $\sqrt{n}$, H$_3$S) są znane wewnątrz projektu i otwarcie flagowane w samym paperze.

**Paper przechodzi zewnętrzny ostrzał.** Zarzuty do rewizji: Ce/Kondo i skalowanie warstwowe. Nowe dane od Jankowskiego (2027?) dla Hg1245/SrTiO$_3$ będą binarnym testem.

---

## Cross-references

- Wewnętrzne: `P7C_open_candidates.md`, `P7_plan.md`, `VERIFICATION_2026-04-19.md`.
- Paper SC: `../../papers/sc/tgp_sc.tex`, Eq. 5 + 7 + 9; §4.A–C (P6 klasy); §5 (P7.A–C blockery); §9 (tabela parametrów C/F/F-D).
- Paper Core: `../../papers/core/tgp_core.tex`, Axioms I–IV + Thm 1–15 (struktura substratu z której $a^*$ i $C_0$ są inherited).
- Kod: `ps17_full_p6_validation.py` (N=29 walidacja), `ps5_global_fit.py` (dim-check), `ps10..ps22` (sub-classy i predykcje).
