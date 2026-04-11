# Koide w TGP — propozycje badań jeszcze nie domknięte ani nie obalone

**Kontekst (notacja).** W dodatku T (`dodatekT_koide_atail_formal.tex`) relacja Koide jest zapisana jako \(Q_K = 3/2\) w postaci

\[
Q_K = \frac{(\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau})^2}{m_e+m_\mu+m_\tau}.
\]

To jest równoważne „Brannen \(K = 2/3\)” w postaci \(K = (m_1+m_2+m_3)/(\sqrt{m_1}+\sqrt{m_2}+\sqrt{m_3})^2 = 2/Q_K\). W kodzie (`tgp_physical_consistency.py`) weryfikowany jest wariant bliski \(2/3\) w tej drugiej konwencji.

**Co jest już „zamknięte” logicznie (nie mylić z wyprowadzeniem \(Q_K\)).** Przy danym \(r_{21}\) warunek \(Q_K=3/2\) **jednoznacznie** wyznacza \(r_{31}\) (wzór kwadratowy / tw. T-r31). Problem O-L5 w manuskrypcie to **nie** „\(r_{31}\) z \(r_{21}\)”, lecz **skąd bierze się sama wartość \(Q_K=3/2\)** (albo równoważnie \(b=\sqrt{2}\) w okręgu Koide, albo \(\mathrm{CV}(\sqrt{m})=1\)).

Poniżej: **propozycje**, które **nie są** zawarte w wyczerpującej liście negatywnych ścieżek z `rem:T-QK-derivation-attempts` (12 pozycji) ani nie zostały jednoznacznie obalone w `rem:T-field-balance`, `rem:T-crosssector-joint`, `rem:T-observer-effect` ani w podsumowaniu `koide_derivation_v47b.py` / `scripts/koide_derivation_output.txt`.

> **Aktualizacja (v47b, ścieżki 13–18b):**  Ścieżka 18b (`koide_observer_critical_v47b.py`) wykazała, że hipoteza „efektu obserwatora" jest **tautologiczna** — \(Q_K(A^4)=3/2\) było wymuszone jako warunek na \(g_0^\tau\), nie odkryte.  Kluczowe odkrycie poboczne: czyste rozstawienie \(\varphi\)-FP daje \(g_0^\tau = \varphi^2 g_0^e = 2{,}27\), co **przekracza próg kolapsu** \(g_{0,\mathrm{crit}} = 8/5 = 1{,}6\).  Soliton \(\tau\) nie istnieje dla \(\varphi^2 g_0^e\) — trzecia generacja jest „wciśnięta" pod próg przez Koide (lub głębszy warunek).  To wzmacnia propozycje #3, #13 i ogranicza propozycję #7.

---

## Wyniki analizy numerycznej (v47b, ścieżki 19–20)

| # | Propozycja | Status | Skrypt |
|---|-----------|--------|--------|
| 13+1 | Korytarz kolapsu + ekwipartycja | Q_K monotoniczne w korytarzu; f=0.845 nie pasuje do stałej; d/d_max = 1/sqrt(2) na simpleksie | `koide_collapse_corridor_v47b.py` |
| 4 | Korekty Z₃ faz | Fazy ODE to 150/256/106 deg (nie Z₃); opt. fazy dla nierównych A to 0/180 deg; korekty overlap łamią Koide | `koide_z3_correction_v47b.py` |

---

## 1. Wyprowadzenie hipotezy ekwipartycji \(r=\sqrt{N-1}\) (Brannen dla \(N\) generacji)

W `thm:T-QK-from-Ngen` pojawia się **hipoteza ekwipartycji** amplitud w simpleksie \((N{-}1)\)-wymiarowym, dająca \(r=\sqrt{N-1}\) i stąd \(Q_K^{(N)}(\sqrt{N-1}) = 2N/(N+1)\), co dla \(N=3\) daje \(3/2\).

- **Dlaczego to nadal otwarte:** twierdzenie jest **warunkowe**; brak dowodu, że substrat TGP / liczenie modów solitonowych **wymusza** dokładnie tę wartość \(r\), a nie np. inną związaną z liczbą stopni swobody pola czy topologią \(\Gamma\).
- **Jak testować:** sformułować mikroskopowy model (np. statystyka fluktuacji na sieci przed coarse-grainingiem) i policzyć, czy rozkład pierwiastków mas (lub ich surrogate z \(A_{\mathrm{tail}}\)) ma **wymuszoną** wariancję zgodną z \(r=\sqrt{2}\) przy \(N=3\).
- **Wynik analizy (v47b):** Na simpleksie frakcji masowych stosunek odległości od centrum d/d_max = 1/\(\sqrt{2}\) = 0.7071 **dokładnie** przy \(Q_K = 3/2\). To jest geometryczna konsekwencja \(b = \sqrt{2}\), nie nowa informacja, ale potwierdza że ekwipartycja ma czyste geometryczne znaczenie. Brak dowodu **dlaczego** ta odległość jest preferowana.

---

## 2. „CV\((A^2)=1\) z symetrii” (T-OP1b) — bez importu Koide jako wejścia

W `cor:T-CV-Atail` pokazano, że \(Q_K=3/2\) jest równoważne warunkowi na drugie momenty \(A^2\) (odchylenie standardowe \(A^2\) równe średniej). Pytanie T-OP1b: czy **TGP z pierwszych zasad** narzuca ten warunek (np. przez \(\mathbb{Z}_3\) na fazach / modach ogona)?

- **Status:** algebraiczna równoważność jest; **implikacja substrat \(\Rightarrow\) CV\(=1\)** nie jest dowiedziona.
- **Jak testować:** szukać tożsamości typu sum rules na \(A_{\mathrm{tail}}(g_0^{(k)})\) wynikających z całek Noether / wiązów wielociałowych, **bez** wstępnego zakładania \(Q_K=3/2\).

---

## 3. Most między Ścieżką 9 (skalowanie \(A_{\mathrm{tail}}\)) a Krzywą Koide

`nbody/examples/ex123_koide_epistemics.md` rozdziela jawnie:

- domknięcie \(r_{31}\) przez **regułę geometryczną** na \(g_0\) (np. \(\varphi^2\) na \(g_0^*\)), od
- domknięcia przez **\(Q_K=3/2\)**.

Te dwie ścieżki **nie muszą** się pokrywać; pokrycie wymagałoby **nowego twierdzenia** łączącego oba języki.

- **Propozycja:** znaleźć transformację \(g_0^{(n)} \mapsto A_{\mathrm{tail}}\) taka, by **jednocześnie** (i) odtwarzała \(r_{21}\) z \(\varphi\)-FP oraz (ii) implikowała \(Q_K=3/2\) **bez** osobnego aksjomatu — np. przez ukrytą zmienną całkową ODE na całym łańcuchu generacji.
- **Ryzyko:** `prop:T-phi2-fail` pokazuje, że naiwne \(\varphi^2\) na \(g_0\) **nie** trafia w PDG; most musi obejmować **nie**-liniowe sprzężenie (Koide już to robi jako selekcja \(g_0^\tau\)).
- **Wzmocnienie (18b):** czyste \(\varphi^2 g_0^e = 2{,}27 > g_{0,\mathrm{crit}} = 1{,}6\) — soliton \(\tau\) **nie istnieje** dla naiwnego \(\varphi\)-FP. Most musi wyjaśniać **dlaczego** trzecia generacja jest ściśnięta pod próg kolapsu.

---

## 4. Drugi (i wyższy) rząd: „prawie \(\mathbb{Z}_3\)” fazy ogona i \(K=2/3\)

W `ex225_koide_geometric_origin.py` (§6) argument z **zerową interakcją** ogona przy sztywnym rozstawie faz \(2\pi/3\) prowadzi do **sprzeczności** z \(K=2/3\) (schematycznie: \(K\to 1\)). Autor skryptu zaznacza **korektę wyższego rzędu** w \(A_i/A_j\), która miałaby przywrócić \(2/3\).

- **Dlaczego nie obalone:** negatywny wynik dotyczy **modelu zerowego rzędu**; **nie ma** w repozytorium pełnego rozwiązania równania na poprawki faz/amplitud aż do uzyskania \(Q_K=3/2\).
- **Jak testować:** rozwinąć energię nakładania ogona \(E_{\mathrm{int}}[h_1,h_2,h_3]\) w TGP wokół konfiguracji \(\mathbb{Z}_3\) i sprawdzić, czy warunek stacjonarności **daje** \(b=\sqrt{2}\) (nie tylko \(Q_K\) bliskie \(3/2\)).
- **Wynik analizy (v47b, `koide_z3_correction_v47b.py`):**
  - Fazy ODE to 150°/256°/106° (nie 120° Z₃). Różnice zależą silnie od \(g_0\).
  - **Optymalne fazy** (min overlap) dla nierównych amplitud to ~0°/180° (nie Z₃!) — bo \(A_\tau\) dominuje i chce być antyrównoległy do \(A_\mu\).
  - **Zero overlap** (ortogonalność) definiuje krzywą w \((d_1, d_2)\), nie punkt. Potrzebny dodatkowy warunek.
  - Korekty overlap do mas **łamią** Koide — nawet przy \(\alpha = 0.01\) przesuwają \(Q_K\) od 3/2.
  - Nachylenie fazowe \(a \cdot \ln\varphi\): \(a\) zmienia się o rząd wielkości między \(g_0^e\) a \(g_0^\tau\). Warunek Z₃ nie jest generycznie spełniony.
  - **Status: osłabiony** — prosty model overlap nie prowadzi do Koide.

---

## 5. Wariacyjne zasady z **wiązami** (nie te same co „PATH 1” w `koide_derivation_v47b`)

Lista zawodów obejmuje m.in. ekstremum **całkowitej masy** \(M\) przy skanowaniu \(Q_K\) — \(M(Q_K)\) jest monotoniczne (**PATH 1 FAIL**). To **nie** wyklucza:

- ekstremum **innego funkcjonału** (np. energii pola z interakcją międzysolitonową, entropii **warunkowej** przy utrzymanym \(r_{21}\) z \(\varphi\)-FP, reszty wiązu Friedmanna/substratu),
- lub zasady typu **Lagrange** z mnożnikami na \(r_{21}\), \(N_{\mathrm{gen}}\), \(g_{0,\mathrm{crit}}\).

- **Jak testować:** sformułować jednoznaczny funkcjonał z fizycznie motywowanymi wiązami i sprawdzić analitycznie lub numerycznie, czy stacjonarny punkt ma \(Q_K=3/2\). (Świadomie inne niż „maksymalna entropia na prostym simpleksie” — ta wersja już pada w PATH 2.)

---

## 6. Warunek „samoenergia \(=\) \(4\times\) energia krzyżowa” w przestrzeni \(A_i\)

`rem:T-field-balance` podaje **algebraiczne** przeformułowanie \(Q_K=3/2\) jako relacji kwartetowej między \(A_i^4\) i iloczynami \(A_i^2A_j^2\). To jest tożsamość na masach, ale sugeruje **interpretację „równowagi pola”**.

- **Otwarte:** czy z **pełniejszej** akcji wielosolitonowej (nie tylko pojedynczy profil) wynika **dynamicznie** ten sam stosunek energii (np. z rozkładu na modowe „samo-” i „krzyżowe” przy zachowaniu wiązu \(\varphi\)-FP)?
- **Jak testować:** zdefiniować energię częściową dla superpozycji trzech solitonów (lub ich asymptotycznych ogona) i porównać iloraz z \(4\) **bez** podstawiania z góry mas PDG.

---

## 7. Fazy ogona: hipoteza \(a\cdot\ln\varphi = 2\pi/3\) (`ex225` / `ex226`)

Idea: jeśli \(\delta(g_0) \approx a\ln g_0 + b\) i \(g_0^{(k)}\) jest \(\varphi\)-drabinką, to \(\Delta\delta = a\ln\varphi\). Równość z \(2\pi/3\) łączyłaby złoty podział z kwantyzacją faz \(\mathbb{Z}_3\).

- **Status:** `ex226_z3_phase_verification.py` jest narzędziem do pomiaru \(\Delta\delta\) dla par \((g_0,\varphi g_0)\); **ostateczna zgodność z \(2\pi/3\)** zależy od kalibracji ODE (\(K=g^4\) vs \(K_{\mathrm{sub}}=g^2\), forma kanoniczna A/B, zakres \(g_0\)). To nie jest „zamknięte na papierze” jako twierdzenie — wymaga **powtarzalnej** numerologii faz w całym łańcuchu generacji, nie tylko dla jednej pary punktów.
- **Propozycja:** przepisać analizę na **te same** \(g_0^e,g_0^\mu,g_0^\tau\), które wychodzą z \(\varphi\)-FP + Koide (jak w `prop:T-g0tau-corrected`), i sprawdzić, czy **trzy** przyrosty faz są addytywne w sensie \(\mathbb{Z}_3\).

---

## 8. Koide na skali energetycznej / biegnące masy (nie tylko „pole masses”)

`scripts/ex171_running_masses_koide.py` i dyskusje w archiwum sugerują, że \(Q_K\) może się lekko przesuwać ze skalą. Jeśli „fizyczna” selekcja \(3/2\) dotyczy np. skali EWSB lub unifikacji, **nie** musi być widoczna w tym samym przybliżeniu co solitonowy \(A_{\mathrm{tail}}\) w \(d=3\).

- **Otwarte:** czy istnieje skala \(\mu^\star\), przy której \(Q_K\) jest **stacjonarne** względem RG w SM+TGP (jeśli TGP wprowadza wspólny mechanizm biegu)? To jest **inne pytanie** niż „czy \(Q_K\) jest niezmiennikiem uniwersalnego skalowania” (ostatnie tłumaczy tylko **stabilność** wartości po skalowaniu, nie jej pochodzenie — zauważone w `ex225` §8).

---

## 9. Geometria informacyjna / statystyka na simpleksie mas

Nie pojawia się na liście 12 zawodów: czy \(Q_K=3/2\) odpowiada **krzywiźnie**, **entropii Fishera**, lub **centrum entropii** rozkładu prawdopodobieństwa zdefiniowanego z \(\sqrt{m_i}/\sum\sqrt{m_j}\) przy naturalnej metryce Shahshahani / Fisher–Rao?

- **Dlaczego warto:** daje niezależną od „energii całkowitej” geometrię na przestrzeni tripletów dodatnich.
- **Ryzyko:** może okazać się czystą reparametryzacją istniejącej algebraicznej równoważności z \(b=\sqrt{2}\).

---

## 10. Substrat dyskretny \(\Gamma\) i momenty rozkładu pola (poza samym \(\Phi_0\))

`scripts/ex162_phi0_from_substrate.py` i powiązane materiały szukają związku \(\Phi_0\) z WF/Ising. Analogicznie: czy **momenty** rozkładu pola (lub jego kwadratu) na sieci w fazie krytycznej **wymuszają** na poziomie efektywnym relację typu \(\mathrm{CV}=1\) dla „efektywnych mas generacji”?

- **Otwarte:** nie widnieje jako zakończony negatywny wynik w `rem:T-QK-derivation-attempts`; wymaga **nowego** mapowania coarse-graining \(\Gamma \to (g_0^{(1)},g_0^{(2)},g_0^{(3)})\).

---

## 11. Kwarki i „Koide poza leptonami” (T-OP4 — nadal szkic)

`rem:T-crosssector-joint`: uniwersalne \(Q_K\) dla wszystkich sektorów **nie** przechodzi testów; to **nie zamyka** pytania, czy istnieje **rozszerzona** reguła (np. z koloru, podciągów \(m+m_0\), czy warstwy konstytuentów), która **redukuje się** do \(3/2\) tylko dla naładowanych leptonów.

- **Propozycja:** sformułować **jedną** zasadę na poziomie TGP (np. z \(N_c\) i śladami po \(\mathrm{SU}(3)_c\)), która przewiduje \(Q_{\mathrm{lept}}\approx 3/2\) oraz odległość kwarków od \(3/2\) **bez** osobnego dopasowania per sektor.

---

## 12. Symetrie Liego / redukcja dla **układu trzech** solitonów

Lista zawodów mówi o braku prostej symetrii skalowania pojedynczej ODE. Otwarte pozostaje przeszukanie **symetrii rozszerzonego** układu (np. grupy działającej na \((g^{(1)},g^{(2)},g^{(3)})\) z wiązami oddziaływania), niekoniecznie redukowalnej do skalowania \(g\to \varphi g\) na jednym profilu.

---

## 13. Analityczny dowód wzoru na \(g_{0,\mathrm{crit}}(\alpha,d)\) dla \(d>1\) i konsekwencje dla „korytarza” trzech generacji

`dodatekT` (dowód \(g_{0,\mathrm{crit}}\) dla \(d=1\), numer dla \(d>1\)): „Dowód analityczny dla \(d>1\) pozostaje otwarty.” To **nie jest** bezpośrednio Koide, ale **ściśle ogranicza** dopuszczalne \(g_0^\tau\) i łączy się z argumentem „trzy generacje pod progiem kolapsu”.

- **Propozycja:** domknąć dowód (Pohozaev–Derrick, perturbacja wokół \(d=1\), lub inna tożsamość) i sprawdzić, czy w **najwęższym** korytarzu dopuszczalnych \(g_0\) funkcja \(Q_K(g_0^\tau)\) ma **wyjątkową** cechę (np. ograniczenie z dołu/z góry zbliżone do \(3/2\)) — to byłaby **pośrednia** selekcja, nie pełne wyprowadzenie.
- **Wzmocnienie (18b):** \(g_0^\tau = 1{,}570\) leży zaledwie \(1{,}9\%\) poniżej \(g_{0,\mathrm{crit}} = 1{,}600\). Jednocześnie \(\varphi^2 g_0^e = 2{,}27\) jest 42 % ponad progiem. Koide (lub głębszy warunek) „ściska" \(g_0^\tau\) tuż pod próg — ścisły dowód \(g_{0,\mathrm{crit}}\) mógłby wyjaśnić, **jak wąski** jest dopuszczalny korytarz i jakie \(Q_K\) z niego wynikają.
- **Wynik analizy (v47b, `koide_collapse_corridor_v47b.py`):**
  - Korytarz: \(g_0^\tau \in [1{,}404, 1{,}600)\), szerokość 0.196.
  - \(Q_K\) jest **monotonicznie malejące** w korytarzu: \(Q_K \in [1{,}27, 2{,}13]\).
  - **Nie ma ekstremum** — \(dQ_K/dg_0^\tau = -7.1\) przy \(Q_K = 3/2\). Korytarz ogranicza zakres \(Q_K\) ale nie selekcjonuje 3/2.
  - Frakcja korytarza \(f = 0{,}845\) nie pasuje do żadnej znanej stałej (\(\varphi, \sqrt{2}, \pi\) itd.; min. odchylenie 16%).
  - \(G_\tau/G_e \approx \pi\) (1.3% błąd) — prawdopodobnie numerologia, nie struktura.
  - **Frakcja zmienia się z \(g_0^e\)** (0.73–0.90) — nie jest uniwersalna.
  - **Status: korytarz to konieczny warunek, nie wystarczający.** \(Q_K = 3/2\) wymaga dodatkowej selekcji.

---

## 14. Wysokie ryzyko spekulacji (świadomie na liście „nie obalone = nie sprawdzone”)

- **Resurgencja / struktura Borela:** czy stała \(\sqrt{2}\) pojawia się jako stosunek Stokesów w asymptotycznym rozwinięciu całek akcji (obecnie brak konstrukcji w TGP).
- **Kodowanie / \(\mathrm{GL}(n,\mathbb{F}_2)\):** `ex259_koide_from_action.py` zawiera eksplorację — warto rozstrzygnąć, czy istnieje **nie**-numerologiczny homomorphism łączący \(b^2=2\) z arytmetyką ciał skończonych i **sprawdziwalną** predykcją poza samym \(Q_K\).

---

## Mapa odniesień w repozytorium

| Temat | Gdzie szukać |
|-------|----------------|
| Algebra \(r_{31}(r_{21})\), status Koide | `dodatekT_koide_atail_formal.tex` |
| 12 zawodnych ścieżek + O-L5 | `rem:T-QK-derivation-attempts`, `rem:T-prediction-hierarchy` |
| Negatywne: równowaga polowa, cross-sector, observer | `rem:T-field-balance`, `rem:T-crosssector-joint`, `rem:T-observer-effect` |
| Negatywne: observer critical (tautologia) | `scripts/koide_observer_critical_v47b.py`, `scripts/koide_observer_critical_output.txt` |
| Skryptowe próby | `scripts/koide_derivation_v47b.py`, `scripts/koide_derivation_output.txt` |
| Epistemika Ścieżka 9 vs Koide | `nbody/examples/ex123_koide_epistemics.md` |
| Fazy \(\mathbb{Z}_3\) vs ODE | `nbody/examples/ex225_koide_geometric_origin.py`, `nbody/examples/ex226_z3_phase_verification.py` |
| Długa sesja analityczna (łańcuch \(\lambda\), OP-\(n\)) | `_archiwum/PLAN_ANALITYCZNY_KOIDE.md`, `_archiwum/ANALIZA_SPRZEZONY.md` |

---

*Plik pomocniczy: propozycje robocze dla dalszych prac nad pochodzeniem \(Q_K=3/2\) (lub \(K=2/3\)) w TGP; nie zastępuje formalnego statusu w manuskrypcie.*
