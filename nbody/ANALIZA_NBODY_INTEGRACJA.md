# TGP `nbody` — analiza integracyjna

*Ostatnia redakcja dokumentu: 2026-04-05. Zakres: folder `TGP_v1/nbody` + powiązanie z głównym tomem (LaTeX w korzeniu repo).*

**Nawigacja:** [README.md](README.md) · [ZALOZENIA_NBODY.md](ZALOZENIA_NBODY.md) · [ANALIZA_3CIALA.md](ANALIZA_3CIALA.md) · [examples/README.md](examples/README.md) · [PLAN_ROZWOJU_NBODY.md](PLAN_ROZWOJU_NBODY.md) · zrzuty: `examples/_outputs/` · snapshot sesji 2026-03: [_archiwum_docs/WYNIKI_SESJI_2026_03_21.md](_archiwum_docs/WYNIKI_SESJI_2026_03_21.md)

**Aktywny tor numeryczny (2026-04):** ścieżka 9 / ogon / Koide — **ex104–ex124** + kontynuacje **ex125–ex137** (szczegóły w `examples/README.md`). Starsze ex (Efimov, K13–K18, …) w `examples/_archiwum/`. Pełna mapa teorii: **`../ANALIZA_KRYTYCZNA_v6.md`**.

---

## 0. Cel i Zakres

Niniejszy dokument jest syntezą dorobku folderu `nbody/` w kontekście rdzenia TGP (numeryka + założenia N0). Odpowiada na trzy pytania:

1. **Czy materiały nbody zdryfowały od rdzenia teorii?** → Nie (patrz §1)
2. **Co analiza nbody potwierdziła, a co obaliła?** → §2–3
3. **Co analiza odkryła nieoczekiwanego?** → §4 (najważniejsza sekcja)

---

## 1. Analiza Dryfu od Aksjomatów N0

### 1.1 Aksjonaty N0 — status w kodzie nbody

| Aksjomat | Treść | Implementacja nbody | Status |
|---|---|---|---|
| **N0-1** | Przestrzeń = pole skalarne Φ | `tgp_field.py`: Φ = Φ₀(1+δ) | ✅ ZACHOWANY |
| **N0-2** | Próżnia: Φ = Φ₀ = const | `beta=gamma=1` domyślnie | ✅ ZACHOWANY |
| **N0-3** | Źródło: δ = C·e^{-m·r}/r (Yukawa) | `yukawa_profile()` w `tgp_field.py` | ✅ ZACHOWANY (Droga B) |
| **N0-4** | Ekwiwalencja masa-pole: C = m/(2√π) | `ex12_path_b_predictions.py` | ✅ ZACHOWANY |
| **N0-5** | V_eff = Φ³/3 − Φ⁴/4, β=γ | `energy_density()`, `screening_mass()` | ✅ ZACHOWANY |
| **N0-6** | m_sp = √(3γ−2β) | `screening_mass()` w `tgp_field.py` | ✅ ZACHOWANY |

**Wniosek §1:** Żaden plik w `/nbody` nie narusza aksjomatów N0. Kod poprawnie implementuje fizykę TGP bez wprowadzania obcych założeń (GR, standardowy model kwantowy, ciemna materia jako odrębny byt).

### 1.2 Jedyna pozorna rozbieżność — wyjaśniona

**Obserwacja:** Profil Yukawa (N0-3) jest w kodzie *narzucony* jako warunek graniczny dla źródła, nie *wyprowadzony* z równania pola TGP.

**Wyjaśnienie (kluczowy wynik ex11):** Linearyzacja ODE defektu przy Φ = Φ₀ daje:

$$\delta'' + \frac{2}{r}\delta' = -\delta \quad \Rightarrow \quad (\nabla^2 + 1)\delta = 0$$

Rozwiązanie: δ ~ sin(r)/r — **ogon OSCYLACYJNY, nie Yukawa.** Profil Yukawa musi być zadany zewnętrznie jako źródło δ-Diraca (Droga B):

$$(\nabla^2 - m_{sp}^2)\Phi = -4\pi C_i \, \delta^3(\mathbf{r} - \mathbf{r}_i)$$

Jest to wewnętrznie spójne z TGP — oznacza jedynie, że ciała są *defektami topologicznymi z nałożonym profile'm Yukawa*, nie swobodnymi rozwiązaniami nieklasyczną ODE. **Droga B = poprawne sformułowanie fizyki cząstkowej w TGP.**

---

## 2. Co nbody POTWIERDZIŁO z Rdzenia Teorii

### 2.1 Naruszenie twierdzenia Earnshawa ✅

> TGP pozwala na **statyczne równowagi** 3+ ciał w wolnej przestrzeni.

W Newtonie twierdzenie Earnshawa (1842) wyklucza statyczne równowagi ładunków w wolnej przestrzeni. W TGP bariera repulsywna V_β ~ 1/d² łamie ten zakaz.

**Wyniki dokładne (analytyczne):**
$$d^2 - 4\beta d + 18\beta C = 0 \implies d_{1,2} = 2\beta \pm \sqrt{4\beta^2 - 18\beta C}$$

Równanie identyczne dla: 2 ciał, trójkąta równobocznego i tetraedru foremnego.

| Konfiguracja | Próg β/C |
|---|---|
| 2 ciała / trójkąt / tetraedr | > **4.5** |
| Kollinear symetryczny | > 4.72 |

Trójkąt jest **geometrycznie preferowaną** pierwszą konfigurację — pojawia się przy niższym β/C.

### 2.2 Siły 3-ciałowe z Φ⁴ ✅

Z nielliniowości $(1+δ_1+δ_2+δ_3)^4$ wynika niezerowy człon 3-ciałowy:

$$V_3 = -6\gamma C_1 C_2 C_3 \cdot I_{\text{triple}}$$

**Granica Coulomba** (m_sp → 0, realistyczna dla TGP fizycznego):

$$I_{\text{triple}}^{\text{Coul}} = \frac{8\pi^2}{P}, \qquad P = d_{12} + d_{13} + d_{23}$$

$$\boxed{V_3^{\text{Coul}} = -\frac{48\pi^2 \gamma C_1 C_2 C_3}{d_{12}+d_{13}+d_{23}}}$$

Skalowanie: **1/P** (liniowe w odwrotności obwodu), nie 1/P² jak błędnie sądzono wcześniej.

### 2.3 Efekt konfinujący i hierarchia reżimów ✅

Trzy człony V₂ mają poprawną hierarchię w TGP:

| Człon | Skalowanie | Fizyczna rola |
|---|---|---|
| V_grad = −4πC²/d | 1/d | Atrakcja dalekozasięgowa |
| V_β = +8πβC²/d² | 1/d² | **Bariera repulsywna** (blokuje kolaps) |
| V_γ = −24πβC³/d³ | 1/d³ | Konfinowanie krótszozasięgowe |

Warunek β/C > 9/2 potwierdza tę hierarchię numerycznie (testy w `run_*.py` i `verify_all.py`).

### 2.4 Całka Feynmana: lepsza od siodłopunktowej ✅

Dokładna całka 2D przez parametryzację Feynmana:

$$I_Y = 2\int_{\Delta_2} \Delta^{-3/2} \, K_0\!\left(m\sqrt{Q/\Delta}\right) d\alpha$$

gdzie $Q = \alpha_2 d_{12}^2 + \alpha_1 d_{13}^2 + \alpha_3 d_{23}^2$, $\Delta = \alpha_1\alpha_2 + \alpha_1\alpha_3 + \alpha_2\alpha_3$.

Błąd aproksymacji siodłopunktowej dla małych t = m·d:

| t = m·d | Błąd siodłowy |
|---|---|
| 0.5 | +60% |
| 1.0 | +86% |
| 2.0 | +148% |
| 5.0 | +329% |

Aproksymacja siodłopunktowa wymaga t > √3 ≈ 1.73. Dla fizycznych parametrów TGP (m_sp mała) należy używać granicy Coulomba lub całki Feynmana.

---

## 3. Co nbody OBALIŁO

### 3.1 TGP jako Ciemna Materia — WYKLUCZONE ❌

> **Ex14 (krzywa rotacji NGC 3198):** TGP-Yukawa *zmniejsza* G_eff(r) = G·e^{-r/λ}.

| Model | v_rot(20 kpc) | χ²/N |
|---|---|---|
| Newton (baryony) | 71 km/s | 198 |
| CDM (ISO) | 197 km/s | 84 |
| TGP λ=100 kpc | **64 km/s** | 220 |

Obserwowane: ~147 km/s. TGP z Yukawa obniża krzywą rotacji, odwrotnie do potrzebnego efektu. Żadna wartość λ nie naprawia tego problemu — jest to konsekwencja struktury siły Yukawa (F ~ e^{-r/λ} · G/r² < Newton).

**Implikacja:** TGP nie zastępuje ciemnej materii w galaktykach. Efekty 3-ciałowe dla zwykłej materii: V₃/V₂ ~ C³/C² = C ~ 10⁻¹⁹ — całkowicie zaniedbywalny wkład.

**Okno zastosowania TGP:** m_sp ~ l_Pl⁻¹ → teoria planckowskiej fizyki; lub m_sp ~ kpc⁻¹ → modyfikacja grawitacji bez DM (wymaga innego mechanizmu płaskiej krzywej rotacji, np. solitonów Φ).

### 3.2 Nieskończona Wieża Stanów Efimova — WYKLUCZONA ❌

Analogia z Efimovem (1970) sugerowała nieskończoną wieżę stanów ze skalowaniem geometrycznym. W TGP:

- V₃ wygaszone Yukawa: I_Y ~ e^{-m·s}/s dla dużych separacji
- **Wynik (ex24, WKB):** max 1 stan WKB w oknie, brak wieży
- Kontrast z Efimowem 1970: tam V₃ ~ 1/R² (nie wygaszone) → nieskończona wieża

**Wniosek:** Stany Efimova TGP są jednopoziomowe — istnieje dokładnie jeden stan kwantowy trójcząstkowy w oknie C_crit(3B) < C < C_crit(2B).

### 3.3 Stacjonarny Profil Yukawa z ODE Pola — WYKLUCZONE ❌

Ex10 (zły znak) → ex11 (poprawione):

Ogon defektu topologicznego przy g → 1 jest **oscylacyjny**, nie Yukawa:

$$\delta'' + \frac{2}{r}\delta' = -\delta \implies \delta(r) \sim \frac{A\sin(r) + B\cos(r)}{r}$$

Profil Yukawa nie wynika ze swobodnych rozwiązań pola TGP — musi być zadany jako źródło punktowe (Droga B).

---

## 4. Kluczowe ODKRYCIA nbody — Integracja z Teorią TGP v24

> To najważniejsza sekcja. Poniższe wyniki nie były przewidywane a priori z rdzenia teorii.

---

### 4.1 ★★★ FAŁSZYWA PRÓŻNIA TGP (ex15)

**Odkrycie:** Potencjał $V(g) = g^3/3 - g^4/4$ ma przy g=1 (próżni TGP) **maksimum**, nie minimum:

$$V''(1) = 2 - 3 = -1 < 0$$

Struktura potencjału:
- g = 0: **prawdziwe minimum** V(0) = 0 ("nicość")
- g = 1: **maksimum** V(1) = 1/12 ("próżnia" TGP)
- g* = 3/4: infleksja (bariera kinetyczna)

**Diagram potencjału:**
```
V(g)
1/12 ┤                ◉ g=1 (próżnia TGP, MAKSIMUM!)
     │               /│\
     │              / │ \
     │             /  │  \
     │            /   │   \
   0 ┤◉──────────/    │    \──────────→ g
     g=0          g*≈0.78  g>1 (rozpad)
```

**Konsekwencje dla TGP:**

1. **Λ_eff > 0 (kosmologiczna):** Energia próżni TGP = V(1) = 1/12 > 0. Zgodne z testem M8 w v24: Λ_eff = V_eff(1) = 1/12. Wszechświat TGP **naturalnie ma dodatnią stałą kosmologiczną**.

2. **Ciała mają niższą energię niż próżnia:** V(g_defect) < V(1). Defekty (cząstki/ciała) leżą energetycznie PONIŻEJ próżni. "Materia = lokalne przejście do prawdziwej próżni g=0."

3. **Ogon oscylacyjny** wynika wprost z V''(1) < 0 (potwierdza ex11).

4. **Stabilizacja kinetyczna:** Bariera kinetyczna f(g) = 1 + 2α·ln(g) → ∞ dla g→0 uniemożliwia globalny kolaps próżni.

5. **Q-ball istnieje** dla ω < 0.5 (pole zespolone Φe^{iωt}) — potencjalny kandidat na bariony TGP.

**Integracja z v24:** Wynik ZGODNY z F-III-4 (Minimalna Zasada Generowania). Test M8 (Λ_eff = 1/12) jest dokładnie energią tej fałszywej próżni. Wynik nbody *potwierdza* i *fizycznie interpretuje* M8.

**Propozycja nowego aksjomatu N0-7 (opcjonalny):**
> *N0-7: Ciała materialne są defektami pola w stanie energetycznym niższym niż próżnia V(Φ₀). Ewolucja układu dąży do przejścia fazowego g: 1 → 0, blokowanego przez barierę kinetyczną f(g).*

---

### 4.2 ★★★ STANY EFIMOVA TGP — 3-Ciałowe Wiązanie bez Wiązania Parowego (ex23–ex33)

**Odkrycie:** Istnieje zakres parametrów C_crit(3B) < C < C_crit(2B), gdzie:
- para jest **NIEWIĄZANA**: E₀(2B+ZP) > 0
- trójka jest **ZWIĄZANA**: E₀(2B+3B+ZP) < 0

**To jest analog Efimova realizowany przez siły 3-ciałowe TGP (V₃ z Φ⁴).**

#### Wyniki — kwantowe okno (1D ekwilatcralny, metoda FD)

$$\boxed{C_{\rm Pl} = \frac{1}{2\sqrt{\pi}} \approx 0.282 \;\in\; \bigl(C_Q^{(3B)},\, C_Q^{(2B)}\bigr) \;\Longleftrightarrow\; m_{\rm sp}\in(0.076,\; 0.198)\,l_{\rm Pl}^{-1}}$$

Charakterystyka stanu kwantowego przy m_sp=0.1, C=C_Pl:
- E₀ = −0.0664 E_Pl (silnie związany)
- ⟨d⟩ = 4.52 l_Pl (3× dalej niż klasyczne minimum)
- σ_d = 2.10 l_Pl (szeroki — efekty kwantowe dominują)
- Frakcja tunelowania: 20.8%
- λ_dB/L_well = 3.03 → **WKB nieważne, wymaga FD**

#### Korekta 2D: Jacobi Schrödinger (ex32–ex33)

Ex30 wykazał, że konfiguracja równoboczna jest **punktem siodłowym** w pełnej przestrzeni kształtów (nie minimum energii dla dowolnych deformacji trójkąta). Właściwe minimum leży w podprzestrzeni izosceles (d12=d13=a, d23=b).

Ex32–ex33 (2D Jacobi) korygują wyniki 1D:

| Parametr | 1D (equilateral) | 2D (isosceles FD) | Korekta |
|---|---|---|---|
| E₀ | −0.00706 E_Pl | **−0.009 E_Pl** | −27% (głębszy) |
| Górna granica m_sp (C_Pl w oknie) | 0.198 l_Pl⁻¹ | **~0.12 l_Pl⁻¹** | zawężenie |

#### ⚠️ REWIZJA KRYTYCZNA: Pełny skan 2D (ex34)

Ex34 wykonał systematyczny skan 2D metodą Jacobi dla m_sp ∈ {0.076, 0.10, 0.12, 0.15, 0.198} i wykazał wynik sprzeczny z 1D:

| m_sp [l_Pl⁻¹] | E0_1D [E_Pl] | E0_2D [E_Pl] | 1D wiąże? | 2D wiąże? | C_Q(2B)_2D | C_Pl w oknie? |
|---|---|---|---|---|---|---|
| 0.076 | −0.106 | **−0.015** | TAK | TAK | 0.276 | ❌ (C_Pl=0.282 > 0.276) |
| 0.100 | −0.063 | **+0.023** | TAK | NIE | 0.307 | ❌ (3B też niewiązany) |
| 0.120 | −0.039 | **+0.043** | TAK | NIE | — | ❌ |
| 0.150 | −0.017 | **+0.062** | TAK | NIE | — | ❌ |
| 0.198 | −0.0004 | **+0.078** | TAK | NIE | — | ❌ |

**Korekta 2D przesuwa energię o +0.08 E_Pl** — wystarczająco dużo, żeby zamknąć okno dla wszystkich testowanych m_sp.

**Mechanizm:** Konfiguracja równoboczna (1D) jest siodłem — ma energię NIŻSZĄ niż prawdziwe minimum 2D izosceles. Dlatego 1D systematycznie zaniża E0 (=większe wiązanie). Korekta ΔE ≈ +0.08 E_Pl jest prawie stała w całym zakresie m_sp.

**Dwie sytuacje przy C=C_Pl=0.282:**

1. **m_sp=0.076:** 3B wiąże (E0_2D=−0.015), ale C_Q(2B)_2D=0.276 < C_Pl → **2B też wiąże** → brak okna Efimova (C_Pl powyżej okna)
2. **m_sp≥0.10:** 3B nie wiąże (E0_2D>0) → C_Pl poniżej progu 3B → **brak okna** (C_Pl poniżej okna)

**Możliwe okno dla C_Pl:** Gęstszy skan w zakresie m_sp ∈ (0.076, 0.090) zlecony ex34v2 (wyniki poniżej).

#### ✅ WYNIK PRECYZYJNY: Ex34v2 — Bisekcja granic okna (siatka Nb=75, Nh=55)

Ex34v2 wykonał gęstszy skan m_sp ∈ [0.076, 0.100], krok 0.003 l_Pl⁻¹, dla obu progów:

**KROK 1 — Skan E0_2D(C_Pl) (z V3=ON → próg 3B):**

| m_sp [l_Pl⁻¹] | E0_2D [E_Pl] | 3B wiąże? |
|---|---|---|
| 0.076 | −0.01298 | TAK |
| 0.079 | −0.00720 | TAK |
| 0.082 | −0.00179 | TAK |
| **0.085** | **+0.00327** | **NIE** |
| 0.088 | +0.00801 | NIE |
| 0.100 | +0.02428 | NIE |

**→ m_sp\* = 0.0831 l_Pl⁻¹** (interpolacja liniowa między 0.082 i 0.085)

**KROK 2 — Próg 2B (m_sp\*\*):**

Metoda V3=OFF w solverze izosceles jest **metodologicznie nieadekwatna** dla progu 2-ciałowego — geometria izosceles zmusza d₁₂=d₁₃=a, co przy b→∞ daje a→∞ (wszystkie odległości duże). Nie można wyizolować pary (1,2). Wynik E0_noV3 ≈ +0.08 dla wszystkich m_sp to artefakt geometrii, nie próg 2B.

**Poprawna wartość m_sp\*\*** pochodzi z interpolacji ex34 (C_Q(2B)_2D vs m_sp):
- C_Q(2B)_2D = 0.276 przy m_sp = 0.076
- C_Q(2B)_2D = 0.307 przy m_sp = 0.100
- Interpolacja dla C_Q(2B)_2D = C_Pl = 0.282:

$$m_{\rm sp}^{**} = 0.076 + \frac{0.282 - 0.276}{0.307 - 0.276} \times 0.024 = 0.0806\; l_{\rm Pl}^{-1}$$

**WERDYKT ex34v2:**

$$\boxed{m_{\rm sp}^{**} \approx 0.0806 < m_{\rm sp}^{*} \approx 0.0831 \implies \text{OKNO EFIMOVA ISTNIEJE!}}$$

$$\text{Okno: } m_{\rm sp} \in (0.0806,\; 0.0831)\; l_{\rm Pl}^{-1}, \quad \Delta m_{\rm sp} \approx 0.0025\; l_{\rm Pl}^{-1}$$

Wąskie, ale fizycznie realne. Wymaga potwierdzenia z:
1. Dedykowanym solverem 2-ciałowym dla precyzyjnego m_sp\*\*
2. Siatką Nb=120, Nh=80 (4× dokładniejszą) dla m_sp\*

**Zaktualizowane okno Efimova 2D (ex34 + ex34v2 + ex46):**

$$C_{\rm Pl} \in \text{okno Efimova 2D} \iff m_{\rm sp} \in (0,\; 0.0831)\; l_{\rm Pl}^{-1}$$

**Ex46 (solver 2-ciałowy, wynik definitywny):** C_Q(2B) ≈ 0.488–0.576 >> C_Pl = 0.282 dla wszystkich m_sp ∈ [0.065, 0.100]. Para TGP **NIGDY nie wiąże** przy C = C_Pl w tym zakresie. Efimov window jest szersze, niż zakładano z ex34.

| m_sp [l_Pl⁻¹] | C_Q(2B) [ex46] | C_Pl = 0.282 | Status pary |
|---|---|---|---|
| 0.065 | 0.488 | 0.282 | NIE WIĄŻE |
| 0.080 | 0.526 | 0.282 | NIE WIĄŻE |
| 0.083 | ~0.531 | 0.282 | NIE WIĄŻE |
| 0.100 | 0.576 | 0.282 | NIE WIĄŻE |

#### Status Kill-shot K13 (po ex34 + ex34v2 + ex46)

> **K13 [POTWIERDZONY]:** Dedykowany solver 2-ciałowy (ex46, 1D Schrödinger FD, μ=1/2) wykazał, że C_Q(2B) ≈ 0.49–0.58 >> C_Pl = 0.282 dla wszystkich m_sp ∈ [0.065, 0.100]. Para TGP jest zawsze niezwiązana przy C = C_Pl. Jednocześnie ex34v2 (V3=ON) pokazał, że trójka wiąże dla m_sp < 0.0831. Okno Efimova (para niewiązana, trójka wiązana) = **m_sp ∈ (0, 0.0831) l_Pl⁻¹** — dużo szersze niż zakładano.

**Interpretacja ex46:** Wcześniejszy szacunek m_sp** ≈ 0.0806 z ex34 interpolacji był oparty na błędnej metodzie (izosceles V3=OFF → problem geometryczny). Ex46 pokazuje, że m_sp** nie istnieje w zakresie [0.065, 0.100] — C_Q(2B) nigdy nie schodzi do C_Pl. Okno Efimova rozciąga się od m_sp=0 (gdzie 3B niezwiązana dla innych powodów) do m_sp* = 0.0831 (gdzie 3B też staje się niezwiązana).

---

### 4.3 ★★ Geometria wielokątów: Skalowanie N-ciałowe (ex28)

Dla regularnych wielokątów N=3,4,5 przy m_sp=0.1, C=0.155:

| N | E_min(3B) | d_eq | V3/V2 | C_crit(3B) |
|---|---|---|---|---|
| 3 (trójkąt) | −0.00706 E_Pl | 4.52 l_Pl | 1.5 | 0.128 |
| 4 (kwadrat) | −0.0685 E_Pl | 2.42 l_Pl | 2.7 | 0.102 |
| 5 (pentagon) | −0.226 E_Pl | 1.69 l_Pl | 3.8 | 0.091 |

**Pentagon wiąże 32× silniej niż trójkąt.** Okno 3B-only przesuwa się ku niższym C dla większych N — efekty N-ciałowe narastają super-ekstensywnie (~N³).

**Fizyczna implikacja:** W reżimie planckowskim, układy N ≫ 3 ciał mogą tworzyć głęboko związane kompaktowe klastery, nawet gdy N=3 jest słabo związany lub niewiązany.

---

### 4.4 ★★ Ograniczenie z Promieniowania (ex13)

Moc promieniowania TGP pola skalarnego:

$$P_{\rm TGP} = \frac{C^2 a^2}{6\pi}\sqrt{1 - \frac{m_{sp}^2}{\omega^2}}\;\theta(\omega - m_{sp})$$

Warunek supresji dla wszystkich obserwowanych systemów (LIGO, pulsary):

$$m_{sp} > \omega_{\rm LIGO} = 3.4\times10^{-41}\,[\text{Planck}]$$

Dla m_sp ~ 1 l_Pl⁻¹: P_TGP = 0 dla **wszystkich** f < f_Planck. TGP jest **automatycznie cichy** (brak promieniowania skalarnego) w całym obserwacyjnym paśmie.

---

### 4.5 ★★★ CIEMNA MATERIA FDM — Sektory TGP (ex35–ex42)

> Sesja ex35–ex42 odkryła, że TGP zawiera **dwa niezależne sektory N-body**, rozdzielone o ~50 rzędów wielkości w skali m_sp.

#### 4.5.1 Dwa Sektory N-body w TGP

| Sektor | Skala m_sp | Fizyka | Testy |
|---|---|---|---|
| **Efimov** | m_sp ~ 0.076–0.12 l_Pl⁻¹ (Planckowski) | 3-ciałowe klastry kwantowe, Earnshaw violation | ex23–ex33 |
| **FDM** | m_sp ~ 8.2×10⁻⁵¹ l_Pl⁻¹ (galaktyczny) | Solitony Φ jako ciemna materia, profil Schive | ex35–ex42 |

Sektory są **dynamicznie niezależne** — nie mieszają się. Paradoks γ (jak jeden m_sp może spełniać oba?) jest **ROZWIĄZANY** dla scenariusza F3: m_sp = m_boson_FDM jest jednym konkretnym parametrem galaktycznym, a sektor Efimov pozostaje otwartym kierunkiem fizyki planckowskiej.

#### 4.5.2 ε_th wywiedziony z N0-6 (ex39) — ZERO nowych parametrów

Kluczowy wynik ex39: próg stabilizacji solitonu ε_th można wyprowadzić **bezpośrednio z aksjomatu N0-6**:

$$m_{sp} = \sqrt{3\gamma - 2\beta} = \sqrt{\gamma} \quad (\beta = \gamma) \implies \varepsilon_{th} = \frac{m_{sp}^2}{2} = \frac{\gamma}{2}$$

Równoważnie: minimum V_mod przy g = φ (złoty podział):

$$V'_{mod}(g) = m_{sp}^2 \cdot g(1 + g - g^2) = 0 \implies g^2 - g - 1 = 0 \implies g = \varphi = \frac{1+\sqrt{5}}{2}$$

**Konsekwencje:**
- ε_th nie jest nowym parametrem — wynika wprost z N0-6
- Złoty podział φ pojawia się naturalnie jako minimim potencjału efektywnego przy progu
- m_boson = m_sp (kwant pola TGP jest bozonem FDM) — jeden parametr opisuje oba sektory

#### 4.5.3 Kill-Shot K18: Skalowanie r_c vs M_gal (ex41)

**Pytanie:** Jak zależy rdzeń solitonu r_c od masy galaktyki M_gal?

Dwa modele TGP:
- **F1 (środowiskowy):** m_boson ∝ M_gal → r_c ∝ M_gal⁻¹ (α = −1)
- **F3 (uniwersalny):** m_boson = m_sp = const → r_c ∝ M_sol^{−1/3} → α ≈ −1/9

**Test na 30 galaktykach SPARC (ex41) i 75 galaktykach SPARC+THINGS+LITTLE THINGS (ex43):**

| Próbka | n | α_obs | σ(α) | ΔAIC(F1−F3) | F3 OK? |
|---|---|---|---|---|---|
| SPARC (ex41) | 30 | −0.096 | 0.028 | +133 | ✅ (1.4σ) |
| SPARC (ex43 replikacja) | 30 | −0.096 | 0.024 | +211 | ✅ |
| THINGS (ex43) | 25 | −0.084 | 0.020 | +169 | ✅ |
| LITTLE THINGS (ex43) | 20 | −0.080 | 0.003 | +135 | ⚠️ |
| **SPARC+THINGS+LTTH. (ex43)** | **75** | **−0.086** | **0.021** | **+515** | **✅ (1.2σ)** |

**Analiza subsampli (ex43) — łamany profil mocy:**

| Subsample | n | α_obs | Dist. od F3 |
|---|---|---|---|
| Spirale (M>10¹⁰ Msun) | 19 | −0.141 | 2.7σ (strome) |
| Karłowe (M<5×10⁹ Msun) | 47 | −0.070 | 5.1σ (płytkie) |
| **Pełna próba** | **75** | **−0.086** | **1.2σ** ✅ |

Różnica spirale–karłowe = 5.2σ. Analiza ex44 (patrz poniżej) wyjaśnia to napięcie feedbackiem barionowym.

**Wynik:** F3 **BARDZO SILNIE PREFEROWANE** nad F1; ΔAIC(n=75) = **+515** (skalowanie ~7/gal, liniowe z n). Prognoza: ΔAIC(n=150) ≈ 1050.

#### ★★★ Wynik ex44 — Wyjaśnienie napięcia spirale/karłowe (K18 v3)

Ex44 przetestował 4 modele dla napięcia 5.2σ:

| Model | Parametry | AIC | ΔAIC vs F3 |
|---|---|---|---|
| M0: F3 (baseline) | 1 | −274.7 | 0.0 |
| M1: Wolny fit | 2 | −366.4 | −91.7 |
| M2: Łamany profil mocy | 4 | −453.2 | −178.5 |
| **M3: F3 + barionowy feedback** | **2** | **−568.5** | **−293.7** ← najlepszy |

**Zwycięzca: M3 — F3 z barionowym feedbackiem:**

$$r_c = \frac{1.61}{m_{22}} \left(\frac{M_{\rm sol}}{10^9 M_\odot}\right)^{-1/3} \cdot \left(\frac{f_{\rm bar}}{f_0}\right)^{0.127} \;\rm kpc$$

Gdzie β_bar = 0.127 i M_break = 5.13×10⁹ M_⊙ (granica spirale–karłowe).

**Wyjaśnienie fizyczne:**
- Spirale (duże f_bar): feedback barionowy (SN, AGN) modyfikuje rdzeń solitonu → r_c_obs ≠ r_c_Schive
- Karłowe (małe f_bar): słabe feedback → r_c_obs ≈ r_c_Schive (czyste F3)
- Efektem jest obserwowane: α_spiral = −0.141 < α_F3 = −1/9 < α_dwarf = −0.070
- F3 pozostaje **fundamentalnym predyktorem** — feedback to poprawka O(10%)

**Wniosek:** Napięcie 5.2σ NIE falsyfikuje F3. Jest to efekt barionowy przewidywany w literaturze (Di Cintio+2014). K18 status: **SILNIE POTWIERDZONY** z wyjaśnieniem subsample tension.

#### 4.5.4 Kill-Shot K14: Napięcie Lyman-α (ex40–ex48)

**Pierwotne napięcie (Irsic+2017):** m₂₂ > 20 vs TGP-FDM m₂₂ ≈ 1 → napięcie 20×

**Po analizie ex42 — napięcie ZREDUKOWANE do 2×:**

| Mechanizm | Status | Uzasadnienie |
|---|---|---|
| (A) Kameleon TGP | ❌ WYKLUCZONY | m_sp/H₀ = 7×10¹⁰ >> 1; pole zamrożone; m_eff(z) ≈ const |
| (B) Mieszana ciemna materia | ❌ WYKLUCZONY | CMB (Hlozek+2018): f_FDM < 4%; soliton za mały |
| (C) Rogers+Peiris 2021 | ✅ ZGODNY | ex48: 1D P_Lya = −0.4σ; m22_fit=1.10 |

**Ex48 — definitywna analiza 1D P_Lya:**

Właściwa obserwabla Lyman-α to **1D power spectrum**:
$$P_{\rm 1D}(k_\parallel) = \frac{1}{2\pi} \int_0^\infty P_{\rm 3D}\!\left(\sqrt{k_\parallel^2 + k_\perp^2}\right) \cdot T_{\rm FDM}^2 \cdot W_T \cdot W_J \cdot k_\perp\, dk_\perp$$

Kluczowy wynik (E-H P_CDM, b_T=20 km/s, k_J=15 h/Mpc):

| Model | S_1D(k=0.3, z=4.5) | Napięcie z S_obs=0.920±0.040 |
|---|---|---|
| m22=1.0 (std FDM) | 0.9308 | −0.3σ |
| m22=1.0 (+TGP δ=0.08) | 0.9348 | **−0.4σ** |
| m22=2.1 (Rogers limit) | 0.9588 | −1.0σ |

**Mechanizm rozwiązania:** Projekcja P_3D → P_1D integruje po k_⊥ — mody z dużym k_⊥ dominują CDM (P_FDM ≈ P_CDM), co rozcieńcza sygnał FDM. Pozorna supresja w P_1D << supresja w P_3D.

**Korekcja TGP:** δ_TGP = C²_Pl·β = 0.0796 → m22_eff = 1.09 (+9%). Rogers+2021 zmierzyłby m22_fit ≈ 1.10 dla TGP-FDM z m22=1.

> **K14 [POTWIERDZONY]:** Napięcie Rogers+2021 było artefaktem metodologicznym — ex47 używał 3D P(k) zamiast właściwej observabli 1D P_Lya. Poprawna analiza (ex48) daje −0.4σ napięcia. TGP-FDM (m22=1) jest ZGODNY z Lyman-α przy 1D projekcji.

#### 4.5.5 Widmo Mocy P(k) TGP-FDM (ex40)

Supresja małej skali FDM (funkcja transferu Irška+2017):

$$T_{FDM}(k) = \left[1 + \left(\alpha_0 k\right)^{2\mu}\right]^{-5/\mu}, \quad \alpha_0 = \frac{0.04}{m_{22}^{4/9}}, \quad \mu = 1.12$$

Skala supresji:
- k₁/₂(m₂₂=1) = 11.3 h/Mpc
- k₁/₂(m₂₂=2) = 15.4 h/Mpc
- k₁/₂(m₂₂=20) = 42.0 h/Mpc

Napięcie z Irsic (wymaga k ~ 42 h/Mpc przy m₂₂=1): napięcie ~3.7× w skali k, a nie 20× w skali masy. Rogers+2021 wymaga jedynie k ~ 15.4 h/Mpc → zgodne przy m₂₂ ≈ 2.

**Skalowanie:** k₁/₂ ∝ m₂₂^{4/9} (weryfikacja numeryczna: współczynnik 3.78 vs 3.79 analitycznie; błąd < 0.2%).

---

## 5. Identyfikacja Pozostałych Problemów Otwartych

### 5.1 Krytyczne (blokują predykcje)

| Problem | Status | Następny krok |
|---|---|---|
| **Q_new-1**: Pełna kwantowa 3D przestrzeń kształtów | ⚠️ OPEN | Ex33 dał 2D (isosceles). Brakuje pełnego solvera 3D w Jacobi |
| **Q_new-2**: Precyzyjny próg 2B dla okna Efimova | ✅ ZAMKNIĘTE | ex46: C_Q(2B)~0.49-0.58 >> C_Pl=0.282; para NIGDY niezwiązana; K13 POTWIERDZONY |
| **Q_new-3**: Mechanizm źródła (dlaczego ciało = Yukawa?) | ⚠️ OPEN | Droga C (zmodyfikowany Lagranżjan z masą skalara) |

### 5.2 Ważne (ciekawa fizyka)

| Problem | Status |
|---|---|
| **Q_new-4**: Lyapunov TGP vs Newton (chaos) | Wymaga długich symulacji |
| **Q_new-5**: Mapa Poincaré 3-ciał w TGP | Numeryczne |
| **Q_new-6**: Orbita ósemkowa Chenciner-Montgomery w TGP | Konfiguracja startowa jest w `configurations.py` |
| **Q_new-7**: Soliton TGP jako ciemna materia (FDM mod.) | ⚠️ CZĘŚCIOWO (ex35-ex42: m_boson=m_sp F3, K18 ✅, K14 ⚠️ marginalnie) |

### 5.3 Techniczne (zaległości)

| Problem | Status |
|---|---|
| Saddle-point warning w `dynamics_v2.py` | Do dodania |
| Integracja ex30–ex33 w `verify_all.py` | Brak testów dla nowych plików |
| Skan 2D kwantowego okna (skrypt ex34) | Do napisania |

---

## 6. Diagram Spójności: nbody ↔ TGP v24

```
TGP v24 Core                 nbody Wyniki               Status
─────────────────────────────────────────────────────────────────
S_Γ z V_eff = g³/3−g⁴/4 ──► V(g=1) jest MAKSIMUM      ✅ SPÓJNE
                               (fałszywa próżnia, Λ>0)    (M8 potwierdzone)

Φ⁴ człon w S_Γ          ──► V₃ = −6γC₁C₂C₃·I_triple   ✅ SPÓJNE
                               Efimov-like 3B states       (nowe predykcje)

m_sp = √(3γ−2β)          ──► Reżimy Coulomba/Yukawa      ✅ SPÓJNE
                               m_sp·d ≶ 1 (graniczne)     (korekcja V₃ skali)

C = m/(2√π) (N0-4)        ──► C_Pl = 0.282 ∈ okno        ✅ SPÓJNE
                               kwantowe Efimova            (K13 propozycja)

dV/dg(g=1) = 0 (M5)       ──► równowagi: d²−4βd+18βC=0  ✅ SPÓJNE
                               Earnshaw violation          (dokładnie)

c(Φ) = c₀√(Φ₀/Φ) (M4)    ──► v_W(χ) = c₀/√χ            ✅ SPÓJNE
                               (v23, O14 zamknięty)        (nie testowane w nbody)

Droga B (źródło Yukawa)   ──► Ogon oscylacyjny → Yukawa  ✅ SPÓJNE
                               musi być NARZUCONY          (wyjaśnia Drogę B)
```

---

## 7. Rekomendacje dla v25

### Priorytet 1 — Natychmiastowe (skrypty)

1. **ex34_corrected_efimov_window.py** — pełny skan 2D okna kwantowego Efimova dla m_sp ∈ [0.05, 0.30] z metodą FD isosceles. Wyznaczyć dokładne granice.

2. **Zaktualizować `verify_all.py`** — dodać testy dla ex30 (saddle-point equilateral), ex32 (isosceles minimum), ex33 (2D Jacobi energy).

### Priorytet 2 — Krótkoterminowe

3. **K13 jako Kill-shot v25:** Dodać do `ANALIZA_SPOJNOSCI_v25.md` predykcję Efimow-klastrów planckowskich jako nowy kill-shot K13.

4. **N0-7 opcjonalny aksjomat:** Fałszywa próżnia i bariera kinetyczna jako dodatkowy aksjomat (lub twierdzenie wynikające z S_Γ).

5. **Soliton TGP (Q_new-7):** Rozwiązanie ex16 (fałszywa próżnia) i ex15 sugeruje badanie stabilizacji próżni przez modyfikację V_eff.

### Priorytet 3 — Średnioterminowe

6. **Chaos (Lyapunov) TGP vs Newton** — 10⁴-krokowe symulacje z `dynamics_v2.py`. Predykcja: TGP ma mniejszy eksponent Lyapunova niż Newton (bariera 1/d² redukuje chaos).

7. **Orbita ósemkowa w TGP** — modyfikacja `configurations.py::figure_eight_initial()` + weryfikacja stabilności z V₃.

---

## 8. Tabela Zbiorcza: Wnioski nbody

| Wynik | Typ | Źródło | Integracja v25 |
|---|---|---|---|
| Earnshaw violation: β/C > 4.5 | Potwierdzony | `equilibria.py`, `stability.py` | ✅ Już w main.tex |
| V₃ ~ 1/P (nie 1/P²) | Korekta | `three_body_terms.py`, `ANALIZA_3CIALA.md` | ✅ Poprawione w equilibria.py |
| Fałszywa próżnia V(1) = max | Odkrycie | `ex15` | ⚠️ Brak w main.tex — **dodać do sek03/sek08** |
| Ogon oscylacyjny defektu | Odkrycie | `ex11` | ⚠️ Częściowe — **uzupełnić sek05** |
| Efimov 3B: C_Pl w oknie [0.076, ~0.12] | Odkrycie | `ex23–ex33` | ⚠️ Brak kill-shot — **K13 do sek09** |
| TGP ≠ ciemna materia | Obalenie | `ex14` | ⚠️ Brak w main.tex — **dodać do sek06** |
| Brak nieskończonej wieży Efimova | Obalenie | `ex24` | ⚠️ Dodać do tgp_efimov_analog.tex |
| Automatyczna supresja promieniowania | Odkrycie | `ex13` | ⚠️ Dodać do sek07 |
| Pentagon silniej związany 32× niż trójkąt | Odkrycie | `ex28` | ⚠️ Nowe predykcje do sek08 |
| 2D Jacobi: E₀ głębszy, okno węższe | Korekta | `ex32–ex33` | ⚠️ Skan ex34 potrzebny |
| ε_th = m_sp²/2 = γ/2 wywiedziony z N0-6 | Odkrycie | `ex39` | ✅ Już w v25 — ZERO nowych parametrów |
| Złoty podział φ jako minimum V_mod | Odkrycie | `ex39` | ✅ Zweryfikowane numerycznie |
| K18: F3 preferowane nad F1 (ΔAIC=133→515) | Weryfikacja | `ex41–ex43` | ✅ sek07 zaktualizowany; n=75 |
| K14: napięcie 2× nie 20× (Rogers+2021) | Rewizja | `ex42` | ⚠️ Otwarte — marginalnie napięte |
| m₂₂ = 2.0 ± 1.5 (predykcja TGP-FDM) | Predykcja | `ex40–ex42` | ⚠️ Falsyfikowalne przez DESI/Rubin |
| Paradoks γ ROZWIĄZANY dla F3 | Rozwiązanie | `ex39–ex42` | ✅ Jeden m_sp = m_boson, λ << 50 kpc |
| K13 POTWIERDZONY (solver 2B definitywny) | Potwierdzenie | `ex34+ex34v2+ex46` | ✅ C_Q(2B)~0.49>>C_Pl; okno Efimova = (0, 0.0831) l_Pl⁻¹ |
| ex47: DESI FDM forecast K14 | Prognoza | `ex47` | ✅ m22_crit=0.96; TGP-F3(m22=1): SNR=2.86 NIE wykluczone (3D aprox.) |
| ex48: K14 ROZWIĄZANY (1D P_Lya) | Rozwiązanie | `ex48` | ✅ 1D E-H: −0.4σ; m22_fit=1.10; K14 POTWIERDZONY; napięcie było artefaktem 3D→1D |
| THINGS+LITTLE THINGS K18 (n=75) | Weryfikacja | `ex43` | ✅ ΔAIC=515; α=−0.086≈−1/9; F1: 44σ |
| K18 napięcie spirale/karłowe wyjaśnione | Wyjaśnienie | `ex44` | ✅ M3(F3+feedback): ΔAIC=−294 najlepszy; β_bar=0.127; M_break=5.1×10⁹ M_⊙ |

---

## 9. Diagram spójności: `nbody` ↔ rdzeń TGP (skrót; szczegóły w `ANALIZA_SPOJNOSCI_v5`)

```
Rdzeń TGP (LaTeX + postulaty)  Wyniki w `nbody/`              Status
──────────────────────────────────────────────────────────────────────────
m_sp = sqrt(gamma) (N0-6) ──► ε_th = m_sp²/2 = γ/2 (ex39)   ✅ ZERO parametrów
                                Złoty podział φ = minimum      ✅ φ wynika z N0-6

Soliton TGP + Schive+2014 ──► r_c ∝ M_sol^{-1/3} (F3)        ✅ K18 PASS (α=-0.096)
                                ΔAIC(F1-F3) = 133              ✅ F3 silnie pref.

T_FDM(k, m22) transfer   ──► k_{1/2}(m22=1) = 11.3 h/Mpc     ⚠️ Rogers: 2× napięcie
                                k_{1/2}(m22=2) = 15.4 h/Mpc   ⚠️ Irsic: 20× (outlier)

V_eff maksimum (ex15)    ──► Λ_eff > 0 naturalnie             ✅ SPÓJNE
Efimov 3B window         ──► m_sp ∈ (0, 0.0831) l_Pl⁻¹        ✅ K13 POTWIERDZONY (ex46)
FDM sektor               ──► m_sp ~ 8.2e-51 l_Pl⁻¹            ✅ niezależny sektor
```

---

*Utrzymanie: dokument historycznie obejmował ex1–ex44; obecny pakiet ma 8 modułów `.py` w korzeniu `nbody/`, **22** pliki `.tex` w `nbody/` oraz **~92** aktywnych `examples/ex*.py` + **55** w `examples/_archiwum/`. Regresja Efimov/okien: `examples/verify_all.py` (**59/59**). Pełny stan mostów i korekt teorii: `ANALIZA_KRYTYCZNA_v6.md`.*
