# TGP N-body: Wyniki sesji 2026-03-21

> Archiwum: kontynuacja [ANALIZA_3CIALA.md](../ANALIZA_3CIALA.md) (stan sprzed porządku 2026-04).
> Status: **zarchiwizowany snapshot** — nie zastępuje bieżącej dokumentacji w korzeniu `nbody/`.

*Indeks bieżący:* [README.md](../README.md) · [examples/README.md](../examples/README.md)

---

## Podsumowanie wykonawcze

Sesja dostarczyła **dziewięciu** definitywnych wyników (ex6–ex33):
1. **Dokładne siły 3-ciałowe** przez całkę Feynmana (ex6–ex7)
2. **ODE defektu topologicznego** z poprawionym znakiem + diagnostyką ogona (ex10–ex11)
3. **Droga B jako zamknięta teoria jednoparametrowa** z predykcjami obserwacyjnymi (ex12–ex13)
4. **Kwantowe okno Efimova vs m_sp** — $C_{\rm Pl}$ w oknie dla $\lambda_{\rm sp}\in(5.1, 13.2)\,l_{\rm Pl}$ (ex27)
5. **Naprawa ex23** — dwuetapowy skan zastąpił `minimize_scalar` (błąd Brenta dla m_sp≥0.18)
6. **Porównanie wielokątów N=3,4,5** — pentagon wiąże 32× silniej, V3/V2=3.8 (ex28)
7. **Funkcja falowa przy C=C_Pl** — 2 stany związane, ⟨d⟩=4.5 l_Pl, WKB nieważna (ex29)
8. **Fizyczna interpretacja okna** — bounds wyjaśnione, klasyczne/kwantowe okna niesprzeczne
9. **2D Jacobi Schrödinger** — E0_2D=-0.009 E_Pl, okno przesuwa się do m_sp<0.12 (ex32-ex33)

**Skrypty weryfikacyjne: 33/33 PASSED** (verify_all.py)
**Dokument PDF: 50 stron** (tgp_yukawa_exact_reduction.pdf — zaktualizowany o tgp_2d_isosceles.tex)

---

## 1. Dokładna całka trójcząstkowa (ex6, ex7)

### Wynik

Calka potrójnego Yukawa obliczana dokładnie przez 2D całkę Feynmana:

$$I_Y = 2 \int_{\Delta_2} \Delta^{-3/2} K_0\!\left(m\sqrt{Q/\Delta}\right) d\alpha$$

gdzie $Q = \alpha_2 d_{12}^2 + \alpha_1 d_{13}^2 + \alpha_3 d_{23}^2$, $\Delta = \alpha_1\alpha_2 + \alpha_1\alpha_3 + \alpha_2\alpha_3$.

Pochodna analityczna:

$$\frac{\partial I_Y}{\partial d_{12}} = -2m \cdot d_{12} \int_{\Delta_2} \alpha_2 \cdot \Delta^{-2} / \sqrt{Q} \cdot K_1(u)\, d\alpha$$

### Weryfikacja (plik `three_body_force_exact.py`)

| Test | Wynik |
|------|-------|
| Newton 3. zasada: \|F₁+F₂+F₃\| | **2.78×10⁻¹⁷** (precyzja maszynowa) |
| Błąd vs różnica skończona | **9.39×10⁻¹¹** |
| Symetria rówboczna | ✓ (|F₁|=|F₂|=|F₃|) |
| Błąd siodłopunktowy | **160–770%** (zły wykładnik: `exp(-1.5t)` zamiast `exp(-√3·t)`) |

### Plik

`three_body_force_exact.py` — gotowy do użycia w dynamice N-ciałowej.

---

## 2. ODE defektu topologicznego (ex10, ex11)

### Problem w ex10 (naprawiony)

Użyto RHS = $g^2(g-1)$ — **zły znak**. Rozwiązania dążą do $g(0) \to -\infty$.

### Poprawne ODE (ex11)

Z $V'(g\Phi_0)/\Phi_0 = \beta g^2 - \gamma g^3 = g^2(1-g)$ dla $\beta=\gamma=1$:

$$f(g)\left[g'' + \frac{2}{r}g'\right] + \frac{\alpha}{g}(g')^2 = g^2(1-g)$$

gdzie $f(g) = 1 + 4\ln g$, granica duchowa $g^* = e^{-1/4} \approx 0.7788$.

### Wyniki numeryczne ex11

| $g_0$ | $g_{\min}$ | $g(r=20)$ | Oscylacje |
|--------|------------|-----------|-----------|
| 0.850 | 0.850 | 0.994 | 12× |
| 0.900 | 0.900 | 0.995 | 12× |
| 0.950 | 0.950 | 0.998 | 12× |

### Kluczowy wniosek: ogon OSCYLACYJNY

Linearyzacja przy $g=1$ ($\delta = 1-g \ll 1$):

$$\delta'' + \frac{2}{r}\delta' = -\delta \quad\Rightarrow\quad (\nabla^2 + 1)\delta = 0$$

Rozwiązanie: $\delta \sim A\sin(r)/r + B\cos(r)/r$ — **oscylacyjne, nie Yukawa**.

Pełne $f(g)$ nie zmienia charakteru ogona (potwierdzone numerycznie: residuum Yukawa fit ~0.003, amplituda oscylacji ~0.1).

**Wniosek:** Yukawa NIE wynika z swobodnego pola TGP. Musi być zadane zewnętrznie.

---

## 3. Droga B — teoria jednoparametrowa (ex12, ex13)

### Aksjomat Drogi B

Ciała są zewnętrznymi źródłami dla pola TGP:

$$(\nabla^2 - m_{sp}^2)\Phi = -4\pi C_i \delta^3(\mathbf{r} - \mathbf{r}_i)$$

Rozwiązanie: $\Phi_i(r) = C_i e^{-m_{sp}r}/r$ (Yukawa — **zdefiniowane**, nie wyprowadzone).

### Wyznaczone parametry

| Parametr | Wartość | Wyznaczony przez |
|----------|---------|-----------------|
| $C_i$ | $m_i/(2\sqrt{\pi})$ [Planck] | Warunek Newtona (1 pomiar: G) |
| $m_{sp}$ | wolny | 2. obserwacja (zasięg siły) |

$$C/m = 1/(2\sqrt{\pi}) = 0.28210 \quad \text{dla WSZYSTKICH ciał (zasada równoważności!)}$$

### N_crit dla protonu

Dla każdego $m_{sp}$: $C_{proton} \approx 2.17\times10^{-20}$ — efekty 3-ciałowe **zupełnie pomijalnie małe** ($N_{crit} > 10^{10}$). Dla ciał makroskopowych: eksponencjalna supresja $e^{-m_{sp}d}$ (t = $m_{sp} \cdot d_{AU} \gg 1$).

### Ograniczenie z promieniowania (ex13)

Moc promieniowania TGP:

$$P_{TGP} = \frac{C^2 a^2}{6\pi}\sqrt{1 - \frac{m_{sp}^2}{\omega^2}}\,\theta(\omega - m_{sp})$$

Dla **supresji promieniowania dla WSZYSTKICH obserwowanych układów** (w tym LIGO @ 100 Hz):

$$m_{sp} > \omega_{LIGO} = 3.4\times10^{-41} \text{ [Planck]}$$

co odpowiada $\lambda < 500\text{ km}$.

Przy $m_{sp} \sim 1$ (skala Plancka): **P_TGP = 0** dla wszystkich $f < f_{Planck} = 1.86\times10^{43}$ Hz → zgodność z LIGO i pulsarami bez żadnego problemu.

### Trzy zakresy $m_{sp}$

| Zakres $m_{sp}$ [Planck] | $\lambda$ | Interpretacja |
|--------------------------|-----------|---------------|
| $< 1.24\times10^{-61}$ | $> 13$ Gpc (Hubble) | TGP = grawitacja Yukawa w skali kosmicznej; problem z promieniowaniem |
| $10^{-58} - 10^{-46}$ | kpc – Gpc | TGP modyfikuje grawitację galaktyczną; promieniowanie supresowane |
| $> 10^{-46}$ | $< 1$ AU | TGP tylko na skalach subnuklearnych/Plancka |

### Optymalny scenariusz: $\lambda \sim 1$ Mpc

- Prawo Newtona zachowane w Układzie Słonecznym (sprawdzone do 50 AU)
- Modyfikacja krzywej rotacji galaktyk: $v_{rot}^2 = G M(r)(1/r + m_{sp})e^{-m_{sp}r}$
- Promieniowanie TGP = 0 dla wszystkich orbit gwiezdnych
- Falsyfikowalny przez anomalie $G_{eff}(r)$ na skalach Mpc

---

## 4. Mapa plików sesji (kompletna)

```
nbody/
├── three_body_force_exact.py        ★ KLUCZOWY — dokładne siły 3-ciałowe
├── tgp_nonstationary_equations.tex       ★ Równania niestacjonarne TGP
├── tgp_nbody_predictions.tex             ★ Predykcje: N_crit, kształt, promieniowanie
├── tgp_quantum_window_msp.tex            ★★ ex27: okno kwantowe vs m_sp — NOWE (32 str. PDF)
└── examples/
    ├── ex6_exact_vs_approx_3body.py  — błąd siodłopunktu: 160-770%
    ├── ex7_nbody_accumulation.py     — N_crit vs (C, m_sp)
    ├── ex8_blockers_analysis.py      — C z Newtona, m_sp wolne
    ├── ex9_fundamental_constraints.py— V3/V2 ~ 10^38 dla gwiazd (katastrofa γ~1)
    ├── ex10_topological_defect.py    — ⚠️ ZŁY ZNAK (ex11 poprawia)
    ├── ex11_full_ode_defect.py       ★ Poprawione ODE + ogon oscylacyjny
    ├── ex12_path_b_predictions.py    ★ Droga B — pełne predykcje
    ├── ex13_radiation_constraint.py  ★ Ograniczenie m_sp z promieniowania
    ├── ex14_rotation_curves.py       ★ TGP vs CDM/MOND (NGC 3198)
    ├── ex15_tgp_soliton.py           ★ V(g) analiza + Q-ball + fałszywa próżnia
    ├── ex16_false_vacuum.py          ★ Krajobraz energii defektu + bańka nukleacji
    ├── ex17_nbody_exact_dynamics.py  ★ Dynamika N-ciałowa z dokładnymi siłami 3-ciałowymi
    ├── ex18_ncrit_exact.py           ★ Dokładna V3/V2 vs saddle-point — mapa obserwowalności
    ├── ex19_planck_scattering.py     ★ Dynamika zderzeniowa 3 cial planckowskich (2B vs 2B+3B)
    ├── ex20_bound_state_3body.py     ★ Stany zwiazane: V3/V2 vs N, geometria, wiriał
    ├── ex21_ccrit_bound_state.py     ★ Krytyczne C — V3 zawsze glębiej o 5-26% (klasycznie)
    ├── ex22_quantum_bound_state.py   ★ ZP + Yukawa — C_crit^QM, m_sp=1 zawsze luzy
    ├── ex23_3b_induced_window.py     ★★ EFIMOV TGP — okno [C_crit(3B), C_crit(2B)] vs m_sp
    ├── ex24_efimov_tower.py          ★ WKB wieża: max 1 stan w oknie (nie nieskonczona wieża)
    ├── ex25_optimal_geometry.py      ★ Optymalna geometria 3B: trójkąt > linia; okno liniowe węższe
    ├── ex26_numerov_bound_state.py   ★★ KWANTOWE OKNO EFIMOVA: FD daje C_Q∈(0.191,0.310), C_Pl∈okno!
    └── ex27_quantum_window_vs_msp.py ★★ OKNO vs m_sp: C_Pl w oknie dla m_sp∈(0.076,0.198) l_Pl^-1
```

---

## 4b. Wyniki rozszerzone (ex14–ex15)

### ex14 — Krzywa rotacji NGC 3198

**Wynik:** TGP Yukawa z $m_{sp}>0$ zmniejsza siłę grawitacji: $G_{eff}(r)=G e^{-r/\lambda} < G$. Krzywa rotacji z TGP jest **niższa** niż Newtonowska (NGC 3198: $v_{TGP}(20\,\text{kpc}) \approx 26-64$ km/s vs obserwowane 147 km/s). TGP **nie może zastąpić** ciemnej materii — pogarsza wynik, nie polepsza.

| Model | $v_{rot}$(20 kpc) | $\chi^2/N$ |
|-------|-------------------|------------|
| Newton (baryony) | 71 km/s | 198 |
| CDM (ISO) | 197 km/s | 84 |
| TGP $\lambda$=100 kpc | 64 km/s | 220 |

### ex15 — Analiza potencjału V(g) i solitony

**Kluczowy wynik analityczny:**

$$V_{TGP}(g) = \frac{g^3}{3} - \frac{g^4}{4} \quad\Rightarrow\quad V''(1) = -1 < 0$$

Próżnia $g=1$ jest **maksimum** potencjału, nie minimum. Konsekwencje:

1. **Potwierdzenie ogona oscylacyjnego** (ex11): $(\nabla^2+1)\delta=0$ wynika wprost z $V''(1)<0$
2. **Defekty energetycznie stabilne**: $E_{defect} < 0$ względem próżni $g=1$ — defekty mają niższą energię niż próżnia. Próżnia $g=1$ to **fałszywa próżnia** TGP; prawdziwa próżnia to $g=0$ (nicość).
3. **Q-ball istnieje** dla $\omega < 0.5$ (pole zespolone $\Phi e^{i\omega t}$) — odpowiednik barionów z ładunkiem U(1)
4. **Minimalne rozszerzenie dla FDM**: $V_{mod}=\varepsilon g^2 + V_{TGP}$ z $\varepsilon\sim 0.05$–$0.1$ daje stabilne solitony z minimum w $g_{min}\sim 0.1$–$0.3$

**Interpretacja kosmologiczna**: Wszechświat TGP jest w fałszywej próżni ($g=1$, energia $1/12$). Defekty (ciała) to fluktuacje w kierunku prawdziwej próżni ($g=0$, energia $0$). Bariera kinetyczna przy $g^*\approx0.778$ zapobiega całkowitemu rozpadowi próżni.

---

## 4c. Wyniki rozszerzone (ex16–ex18)

### ex16 — Krajobraz energii defektu i fałszywa próżnia

Patrz sekcja 4b powyżej (false vacuum).

### ex17 — Dokładna dynamika N-ciałowa

**Wynik:** Integracja orbit z pełnymi siłami 2-ciałowymi i 3-ciałowymi (Feynman) — zachowanie energii $|\delta E/E| < 5\times 10^{-7}$.

| Test | Wynik |
|------|-------|
| V3/V2 (C=0.20, d=2) | 1.36% |
| V3/V2 (C=0.50, d=2) | 3.39% |
| $|\delta E/E|$ max | $4.58\times10^{-7}$ |
| Odchylenie orbit (C=0.5, t=5) | $|\delta r|=8.30\times10^{-2}$ |

### ex18 — Dokładna N_crit vs aproksymacja siodłopunktowa

**Wyniki kluczowe:**

#### V3/V2 jako funkcja t = m_sp · d (C=0.20, trójkąt równoboczny)

| t=m·d | V3/V2 (dokładna) | V3/V2 (siodło) | Błąd siodła |
|-------|-----------------|-----------------|-------------|
| 0.5  | 0.580  | 0.929  | +60%  |
| 1.0  | 0.484  | 0.898  | +86%  |
| 2.0  | 0.253  | 0.627  | +148% |
| 3.0  | 0.121  | 0.373  | +208% |
| 5.0  | 0.026  | 0.112  | +329% |

**Przyczyna błędu siodła**: przybliżenie $K_0(u)\approx\sqrt{\pi/2u}\,e^{-u}$ wymaga $u>3$, tzn. dla trójkąta równobocznego $t>1.73$. Dla $t<1.73$ błąd jest ogromny.

#### Kryterium obserwowalności (V3/V2 > 1%)

| t=m·d | C_crit [Planck] | m_body (kg) |
|-------|----------------|-------------|
| 0.5  | $3.4\times10^{-3}$ | $2.7\times10^{-10}$ |
| 1.0  | $4.1\times10^{-3}$ | $3.2\times10^{-10}$ |
| 2.0  | $7.9\times10^{-3}$ | $6.1\times10^{-10}$ |
| 4.0  | $3.5\times10^{-2}$ | $2.7\times10^{-9}$  |

**Próg masy**: $m_{body} > \sim 3\times10^{-10}$ kg (ok. 300 ng — mały pył Plancka).

#### Fizyczne systemy

| System | C [Pl] | V3/V2 | Obserwowalne? |
|--------|--------|-------|---------------|
| Proton | $2.2\times10^{-20}$ | $\sim10^{-20}$ | NIE |
| Elektron (atom) | $1.2\times10^{-23}$ | $\sim10^{-23}$ | NIE |
| Mikropyłek 100nm | $1.3\times10^{-14}$ | $3.7\times10^{-14}$ | NIE |
| Hipotetyczny obiekt Plancka | $0.28$ | **68%** | TAK |

**Wniosek główny**: Efekty 3-ciałowe TGP są **obserwowalne wyłącznie dla obiektów o masie $\gtrsim 300$ ng** przy odległości planckowskiej. Dla wszystkich znanych cząstek fundamentalnych: nieobserwowalne ($V3/V2\sim 10^{-20}$).

Dodatkowy wniosek: **dokładna całka Feynmana daje V3/V2 nawet 3–4× większe** niż saddle-point dla małych t — oznacza to, że efekty 3-ciałowe były dotąd systematycznie **niedoszacowywane** przez aproksymacje analityczne.

---

### ex19 — Dynamika zderzeniowa obiektów planckowskich

**Metoda:** Velocity Verlet, $\Delta t=0.02$, $T_{\rm max}=25$ [Planck]. Siły: V2 (Yukawa) + V3 (całka Feynmana, wektoryzowana).

#### Energia potencjalna 2B vs 2B+3B (trójkąt równoboczny)

| d [l_Pl] | V2 [E_Pl] | V3 [E_Pl] | V3/V2 |
|----------|-----------|-----------|-------|
| 1.0 | −0.0878 | −0.0200 | **22.8%** |
| 1.5 | −0.0355 | −0.0060 | 16.8% |
| 2.0 | −0.0162 | −0.0019 | 11.9% |
| 3.0 | −0.0040 | −0.0002 | 5.7% |

**Wniosek:** dla $d \sim 1$–$2\,l_{\rm Pl}$ efekty 3-ciałowe to **10–23% energii** — nie poprawka perturbacyjna, ale wkład wiodący!

#### Odchylenie obserwatora (zderzenie czołowe 1+2, b=2)

| t [t_Pl] | y₃ (2B) | y₃ (2B+3B) | |δy₃| |
|----------|---------|------------|------|
| 10 | 2.087 | 2.082 | 4.1×10⁻³ |
| 20 | 2.764 | 2.621 | 1.4×10⁻¹ |
| 25 | 3.234 | 2.997 | **0.24 l_Pl** |

Siła 3-ciałowa przesuwa obserwatora o **0.24 l_Pl** w ciągu $25\,t_{\rm Pl}$.

#### Stabilność trójkąta planckowskiego

Trójkąt (d=2) z rotacją ω=0.1: system jest niewiązany, ale **z siłami 3-ciałowymi rozszerza się wolniej** o ~3% (d12_2B − d12_3B ≈ 0.19 przy t=20). Wkład energetyczny V3 = −0.00192 E_Pl pogłębia stan (ale nie wystarcza do związania).

---

## 4d. Wyniki kluczowe (ex20–ex23): Stany Efimova TGP

### ex20 — Stany związane: V3 super-ekstensywne

Dla wielokąta $N$ cząstek, $d=0.8\,l_{\rm Pl}$:

| N | V3/V2 | E/N [E_Pl] |
|---|---|---|
| 3 | 25.1% | −0.056 |
| 5 | **59.6%** | −0.098 |

V3/V2 rośnie z N (trójki mnożą się jak $\binom{N}{3}$). Trójkąt równoboczny jest **najgłębszą** konfiguracją.

### ex21 — Klasycznie: V3 zawsze głębiej o 5–26%

Twierdzenie wirialne: zarówno 2B jak i 2B+3B są zawsze związane klasycznie (Yukawa $\to -\infty$ przy $r\to 0$). Ale $V_3$ pogłębia o **5–26%** zależnie od $C$.

### ex22 — Kwantowo: ZP blokuje wiązanie przy m_sp=1

$E_{\rm ZP} = 3/(8d^2)$ wypycha system do $d\to\infty$. Dla $m_{\rm sp}=1$ (Planck): żadne C nie daje stanu związanego. Dla $m_{\rm sp}=0.1$: $C_{\rm crit}^{(2B)} = 0.184$.

### ★★ ex23 — STANY EFIMOVA TGP (KLUCZOWY WYNIK)

**Okno 3B-induced**: dla $C_{\rm crit}^{(3B)} < C < C_{\rm crit}^{(2B)}$ para jest **niewiązana**, ale trójka jest **związana** wyłącznie przez $V_3$!

| $m_{\rm sp}$ | $\lambda$ | $C_{\rm crit}^{(2B)}$ | $C_{\rm crit}^{(3B)}$ | Szerokość okna |
|---|---|---|---|---|
| 0.05 | 20 l_Pl | 0.130 | 0.084 | 0.047 |
| 0.10 | 10 l_Pl | 0.184 | 0.128 | **0.057** |
| 0.20 | 5 l_Pl | **0.261** | **0.193** | **0.067** | *(poprawione)* |
| 0.25 | 4.0 l_Pl | **0.292** | **0.221** | **0.071** | *(poprawione)* |
| 0.30 | 3.3 l_Pl | **0.319** | **0.246** | **0.074** | *(poprawione)* |

*(Stare wartości m_sp=0.20/0.25/0.30 były błędne — `minimize_scalar` nie znajdowało płytkiego minimum przy d≈5. Naprawiono dwuetapowym skanem.)*

**Dla $m_{\rm sp}=0.2364$ (λ=4.23 l_Pl), $C=C_{\rm Pl}=0.282$:** *(pierwsze m_sp gdzie C_Pl w klasycznym oknie)*
- $E_{\rm min}^{(2B+ZP)} = +0.000195 > 0$ (niewiązany przez 2-ciała)
- $E_{\rm min}^{(2B+3B+ZP)} = -0.073 < 0$ (**związany przez 3-ciała!**)
- $d_{\rm eq} = 1.61\,l_{\rm Pl}$

**Dla $m_{\rm sp}=0.1$, okno 3B-only:**
- $C \in [0.128, 0.184]$ → $m_{\rm body} \in [9.8\times10^{-9}, 1.4\times10^{-8}]$ kg
- Energia wiązania: $E_{\rm bind} \sim 10^{-3}$–$10^{-2}\,E_{\rm Pl}$
- $d_{\rm eq} \sim 3$–$6\,l_{\rm Pl}$

**Analogia z fizyką nuklearną**: TGP Efimov = stany związane trójki bez wiązania parowego. Unikalnie falsyfikowalna predykcja TGP.

### ex24 — WKB kwantyzacja: skończona wieża (max 1 stan)

**Metoda**: warunek WKB $n(E)=\frac{1}{\pi}\int_{d_L}^{d_R}\sqrt{2\mu(E-V_{\rm eff})}\,dd = n+\tfrac{1}{2}$, $\mu=1\,m_{\rm Pl}$.

**Wyniki dla $m_{\rm sp}=0.1$:**

| C | N_stany | $E_0$ [$E_{\rm Pl}$] | $d_{\rm eq}$ [$l_{\rm Pl}$] |
|---|---------|----------------------|-----------------------------|
| 0.128–0.175 | 0 | — | 5–8 (poniżej progu WKB) |
| 0.175 | 1 | −0.00029 | 3.48 |
| 0.180 | 1 | −0.00090 | 3.26 |
| 0.184 | 1 | −0.00156 | 3.04 |

Dla $C=0.180$: $d_L=1.9$, $d_R=16.3$, $\Delta d=14.4\,l_{\rm Pl}$ — stan **bardzo rozlany**.

**Wniosek**: TGP ma **co najwyżej 1 stan WKB** w oknie Efimova. Brak nieskończonej wieży — $V_3$ jest wygaszone Yukawa, nie $\sim 1/R^2$.

---

### ★ ex25 — Optymalna geometria 3B

**Pytanie**: Czy trójkąt równoboczny jest najgłębszą konfiguracją 3-ciałową?

**Wyniki dla $m_{\rm sp}=0.1$, $C=0.155$ (środek okna Efimova):**

| Geometria | N | $E_{\rm min}$ [$E_{\rm Pl}$] | $d_{\rm ref}$ [$l_{\rm Pl}$] | Związany? |
|-----------|---|-------------------------------|-------------------------------|-----------|
| Trójkąt równoboczny | 3 | **−0.00706** | 4.52 | TAK |
| Łańcuch liniowy     | 3 | −0.00026      | 6.25 | TAK (słabo) |
| Kwadrat (N=4)       | 4 | −0.06842      | 2.48 | TAK (głęboko) |

- Linia jest **27× słabiej związana** niż trójkąt przy $C=0.155$.
- Kwadrat jest **~10× głębiej związany** dzięki $\binom{4}{3}=4$ trójkom.

**Okna 3B-only dla $m_{\rm sp}=0.1$:**

| Geometria | $C_{\rm crit}^{(3B)}$ | $C_{\rm crit}^{(2B)}$ | Szerokość okna |
|-----------|----------------------|----------------------|----------------|
| Trójkąt równoboczny | ≈ 0.128 | ≈ 0.184 | **0.056** |
| Łańcuch liniowy     | ≈ 0.155 | ≈ 0.184 | 0.029 |

Linia ma okno **2× węższe** — boczna para $r_{13}=2d$ jest poza zasięgiem Yukawa ($\lambda=10\,l_{\rm Pl}$), co tłumi $I_Y$.

**Wniosek:** Trójkąt równoboczny jest preferowaną konfiguracją. Linia jest związana, ale w węższym oknie i z mniejszą energią wiązania. Falsyfikowalna predykcja TGP dotyczy przede wszystkim **trójkątnych** klastrów.

---

### ★★ ex26 — Kwantowe okno Efimova (FD, KLUCZOWY WYNIK)

**Problem:** WKB (ex24) był niepoprawny — długość de Broglie'a $\lambda_{\rm dB}\approx44\,l_{\rm Pl} \gg L_{\rm well}\approx14\,l_{\rm Pl}$ przy $C=0.180$. Potrzebna pełna kwantowa diagonalizacja.

**Metoda:** Skończone różnice dla $\hat{H}=-\frac{1}{2\mu}\partial^2/\partial d^2 + V_{\rm eff}(d)$, $N=3000$, $d\in[0.4,30]\,l_{\rm Pl}$.

**Wyniki (m_sp=0.1):**

| Kryterium | $C_{\rm crit}^{(3B)}$ | $C_{\rm crit}^{(2B)}$ | Szerokość $\Delta C$ | Masa $[m_{\rm Pl}]$ |
|-----------|----------------------|----------------------|----------------------|----------------------|
| Klasyczne (ex23) | 0.128 | 0.184 | 0.056 | [0.45, 0.65] |
| **Kwantowe (FD)** | **0.191** | **0.310** | **0.119** | **[0.68, 1.10]** |

**Kwantowe okno jest 2.1× szersze i przesunięte do wyższego C.**

**★ Planck = środek kwantowego okna:**
$$C_{\rm Pl} = 0.282 \;\in\; (0.191, 0.310)$$
- $E_0^{(2B)}(C_{\rm Pl}) = +0.0022\,E_{\rm Pl} > 0$ (trójkąt niewiązany przez 2-ciała)
- $E_0^{(3B)}(C_{\rm Pl}) = -0.0664\,E_{\rm Pl} < 0$ (**silnie związany przez 3-ciała!**)
- $d_{\rm peak}\approx 4\,l_{\rm Pl}$, stany 2 przy $C=0.280$–$0.310$

**Dlaczego WKB był zły:** Dla płytkich studzienek ($E\approx0$) przybliżenie WKB przeszacowuje wiązanie — stan kwantowy "wychodzi poza" zasięg studzienki przez efekt tunelowania.

**Stan falowy przy $C=0.207$:**
- $E_0=-4.2\times10^{-3}\,E_{\rm Pl}$, $d_{\rm peak}=7.1\,l_{\rm Pl}$, $\langle d\rangle=9.8\,l_{\rm Pl}$, $\sigma_d=4.8\,l_{\rm Pl}$

---

## 5. Otwarte pytania

1. **Mechanizm źródła (Q1):** Dlaczego ciało jest źródłem Yukawa? **Częściowo omówione** (`tgp_open_questions.tex` §Q1) — spójność z teorią pola: g=1/m_Pl, φ³ vertex daje 3-ciałowe TGP. Pełna Lagranżjana TGP pozostaje otwarta.

2. ~~**Ciemna materia vs TGP (Q2):**~~ **ODPOWIEDZIANO (negatywnie)** (`tgp_open_questions.tex` §Q2 + ex31) — TGP Drogi B (materia zwykła, C~10^{-19}) NIE może wyjaśnić krzywych rotacji: Yukawa redukuje grawitację dla r>1/m_sp (odwrotnie do obserwacji). Efekt 3-ciałowy C³~10^{-60} — całkowicie pomijalny.

3. **Droga C (kwantowe defekty, Q3):** **Omówione** (`tgp_open_questions.tex` §Q3) — ogon defektu jest oscylacyjny (sin(r)/r), nie Yukawa. Potrzeba zmodyfikowanego Lagranżjana z masą skalara.

4. **Związek z GR (Q4):** **Omówione** (`tgp_open_questions.tex` §Q4) — TGP i GR koegzystują. TGP przy m_sp~l_Pl^{-1} jest modyfikacją grawitacji w skali Plancka, zgodną z GR i testami Układu Słonecznego. Ograniczenie BD: m_sp>1/AU spełnione.

5. ~~**Pełna kwantowa kalkulacja stanów Efimova (Q5):**~~ **POSTĘP (ex30)** — Konfiguracja równoboczna jest punktem siodłowym, nie minimum! Aproksymacja 1D daje górne ograniczenie na E_0. Prawdziwe okno może być szersze. Pełny kalkulator 2D FD zidentyfikowany jako główny otwarty problem numeryczny.

6. ~~**Nieskończona wieża stanów Efimova:**~~ **ROZWIĄZANE (ex24)** — TGP ma max 1 stan WKB w oknie. Brak nieskończonej wieży — $V_3$ wygaszone Yukawa.

7. ~~**Weryfikacja ex23 dla geometrii nietrójkątnych:**~~ **ROZWIĄZANE (ex25)** — Linia ma okno 3B-only ($C\in[0.155,\,0.184]$), ale 2× węższe i 27× słabsze niż trójkąt.

8. ~~**Stany wzbudzone w oknie:**~~ **ROZWIĄZANE (ex26)** — FD potwierdza 1 stan w kwantowym oknie; WKB był błędny.

9. ~~**Optymalna geometria N=5:**~~ **ROZWIĄZANE (ex28)** — Pentagon N=5 wiąże przy C≈0.091 (vs 0.128 dla N=3). Energia wiązania 32× większa niż dla trójkąta przy C=0.155. V3/V2=3.76 — siły 3-ciałowe dominują. Okno przesuwa się ku niższym C dla większych N.

10. ~~**m_sp zależność kwantowego okna:**~~ **ROZWIĄZANE (ex27)** — $C_{\rm Pl}$ jest w kwantowym oknie dla $m_{\rm sp}\in(0.076, 0.198)\,l_{\rm Pl}^{-1}$, czyli $\lambda\in(5.1, 13.2)\,l_{\rm Pl}$. Patrz wyniki ex27 poniżej.

11. ~~**Dynamika kwantowego klastra przy $C=C_{\rm Pl}$:**~~ **ROZWIĄZANE (ex29)** — Pełna charakteryzacja funkcji falowej: ⟨d⟩=4.52 l_Pl (3× dalej niż klasyczny minimum), σ=2.10 l_Pl, 20.8% tunelowania, WKB nieważna (λ_dB/L_well=3.03). Dwa stany związane przy m_sp=0.1, C=C_Pl.

12. ~~**Fizyczna interpretacja granic okna:**~~ **ROZWIĄZANE (tgp_physical_interpretation.tex)** — Górna granica λ=13.2 l_Pl: zasięg przy którym pary kwantowo wiążą przy C_Pl (C_Q(2B)=C_Pl). Dolna granica λ=5.1 l_Pl: zasięg przy którym V3 zbyt słabe by pokonać ZP (C_Q(3B)=C_Pl). Okno = zakres zasięgu skalara gdzie 3B wygrywa z ZP, ale 2B jeszcze nie wiąże.

13. ~~**Dlaczego kwantowe i klasyczne okno nie nakładają się?**~~ **ROZWIĄZANE (tgp_physical_interpretation.tex)** — Kwantowa ZP jako operator przesuwa oba progi o ΔC≈0.063 ku wyższym C. Dla m_sp∈(0.235,0.38) (klasyczne okno): studnia jest wąska (Δd∼2 l_Pl), ZP-operator niszczy wiązanie (λ_dB/L_well∼3). Dla m_sp∈(0.076,0.198) (kwantowe okno): studnia szersza, funkcja falowa może być zawarta. Dwa reżimy wzajemnie wykluczają się przy C=C_Pl.

14. ~~**Poprawienie ex23 o poprawny algorytm:**~~ **ROZWIĄZANE** — ex23 zaktualizowany do dwuetapowego skanu (600+500 punktów). `minimize_scalar` usunięte. Wszystkie progi klasyczne m_sp=0.15–0.30 ponownie obliczone i zaktualizowane w `tgp_efimov_analog.tex` oraz `verify_all.py`.

---

### ★ verify_all — Skrypt weryfikacyjny (33/33 testów)

**Plik:** `examples/verify_all.py`

**Zakres weryfikacji:**
1. Samospójność całki Feynmana: `I_Y_equil == I_Y_general` dla 4 punktów ✓
2. Baryczny warunek `sum(A1+A2+A3)=1` ✓
3. Konwencja ZP = `3/(8d²)` ✓
4. Klasyczne okno Efimova przy m_sp=0.1 (vs ex23) ✓
5. Wartości energii przy C=0.155, m_sp=0.1 ✓
6. Klasyczne okno Efimova przy m_sp=0.2 (vs. **POPRAWIONE** wartości) ✓
7. Kwantowe okno (FD) przy m_sp=0.1 vs ex26 ✓
8. Kwantowe okno przy granicach m_sp=0.18 i 0.20 ✓
9. Eigenvalues FD przy kluczowych C ✓
10. Energetyka statyczna (vs tgp_scattering.tex) ✓
11. Spójność C_Pl ✓

**Wynik: 33 PASSED, 0 FAILED**

### ⚠️ BŁĄD ZNALEZIONY I NAPRAWIONY: klasyczne progi ex23

**Bug:** `minimize_scalar(method='bounded')` w ex23 nie znajdowało płytkiego minimum dla m_sp≥0.18. Powodem jest algorytm Brenta: złoty podział daje punkt startowy po prawej stronie prawdziwego minimum (d≈5), a wtedy Brent zbiega do granicy d=30.

**Wpływ:** Wartości w `tgp_efimov_analog.tex` (tab:Ccrit_window) dla m_sp=0.20, 0.25, 0.30 były BŁĘDNE:

| m_sp | Old C_cl(2B) | **Correct C_cl(2B)** | Old C_cl(3B) | **Correct C_cl(3B)** |
|------|-------------|---------------------|-------------|---------------------|
| 0.20 | 0.281 ❌ | **0.261** ✓ | 0.234 ❌ | **0.193** ✓ |
| 0.25 | 0.365 ❌ | **0.292** ✓ | 0.319 ❌ | **0.221** ✓ |
| 0.30 | 0.478 ❌ | **0.319** ✓ | 0.433 ❌ | **0.246** ✓ |

**Konsekwencja:** Klasyczne okno Efimova rośnie z m_sp (ΔC od 0.047→0.074), a nie jest stałe (~0.05 jak błędnie podano). C_Pl jest w klasycznym oknie dla m_sp∈(~0.235, ~0.38).

**Poprawka:** Zastosowanie dwuetapowego skanu (600 grubych + 500 lokalnych punktów). Zaktualizowano:
- `tgp_efimov_analog.tex` — poprawione wartości w tabeli, nowa uwaga numeryczna
- `tgp_quantum_window_msp.tex` — poprawiona tab:window_comparison, nowy opis kontekstu
- `examples/verify_all.py` — poprawny algorytm skanowania, nowe wartości referencyjne

**Kluczowy wniosek z korekty:** Kwantowe okno Efimova (ex26/ex27) i klasyczne okno są **niesprzeczne** — obejmują różne zakresy m_sp:
- Kwantowe okno dla C_Pl: m_sp∈(0.076, 0.198)
- Klasyczne okno dla C_Pl: m_sp∈(~0.235, ~0.38)
- **Brak przekrycia** — co oznacza, że klasyczna analiza i kwantowa dają jakościowo przeciwne predykcje dla większości przestrzeni parametrów!

---

### ★★ ex27 — Kwantowe okno Efimova vs $m_{\rm sp}$ (KLUCZOWY WYNIK)

**Pytanie:** Dla jakich wartości masy skalara $m_{\rm sp}$ obiekt Plancka ($C=C_{\rm Pl}=0.282$) tworzy kwantowo-mechanicznie związany stan trójciałowy?

**Metoda:** Powtórzenie ex26 (FD Schrödingera) dla $m_{\rm sp}\in[0.05, 1.00]\,l_{\rm Pl}^{-1}$. Bisekcja do znalezienia $C_Q^{(3B)}(m_{\rm sp})$ i $C_Q^{(2B)}(m_{\rm sp})$. Sprawdzenie $C_{\rm Pl}\in(C_Q^{(3B)}, C_Q^{(2B)})$.

**Wyniki:**

| $m_{\rm sp}$ | $\lambda$ | $C_{\rm cl}^{(3B)}$ | $C_{\rm cl}^{(2B)}$ | $C_Q^{(3B)}$ | $C_Q^{(2B)}$ | $\Delta C_Q$ | $C_{\rm Pl}$ w oknie? |
|---|---|---|---|---|---|---|---|
| 0.05 | 20.0 | 0.084 | 0.130 | 0.135 | 0.252 | 0.116 | NIE |
| 0.08 | 12.5 | 0.111 | 0.165 | 0.170 | 0.287 | 0.118 | **TAK** |
| 0.10 | 10.0 | 0.128 | 0.184 | 0.191 | 0.310 | 0.119 | **TAK** |
| 0.12 | 8.3  | 0.142 | 0.202 | 0.211 | 0.332 | 0.122 | **TAK** |
| 0.15 | 6.7  | 0.163 | 0.226 | 0.239 | 0.365 | 0.126 | **TAK** |
| 0.18 | 5.6  | 0.182 | 0.247 | 0.267 | 0.396 | 0.130 | **TAK** |
| 0.20 | 5.0  | 0.193 | 0.261 | 0.284 | 0.416 | 0.132 | NIE |
| 0.25 | 4.0  | 0.221 | 0.291 | 0.326 | 0.464 | 0.139 | NIE |
| 0.30 | 3.3  | 0.246 | 0.319 | 0.365 | 0.510 | 0.145 | NIE |
| 0.80 | 1.2  | 0.433 | 0.521 | --- | --- | --- | NIE |

**Dokładne granice (bisekcja 30 iteracji):**

$$\boxed{C_{\rm Pl}\in\bigl(C_Q^{(3B)},\,C_Q^{(2B)}\bigr) \;\Longleftrightarrow\; m_{\rm sp}\in(0.076,\;0.198)\,l_{\rm Pl}^{-1}}$$

co odpowiada zakresowi zasięgu skalara:

$$\lambda_{\rm sp} \in (5.1,\;13.2)\,l_{\rm Pl} = (8,\;21)\times10^{-35}\,\text{m}$$

**Obserwacje:**
- Kwantowe okno jest **~2× szersze** od klasycznego ($\Delta C_Q/\Delta C_{\rm cl}\approx2.0$)
- Oba progi rosną z $m_{\rm sp}$: krótszy zasięg wymaga silniejszego sprzężenia
- Kwantowe okno **zamyka się** przy $m_{\rm sp}\approx0.7\,l_{\rm Pl}^{-1}$ ($\lambda\approx1.4\,l_{\rm Pl}$)
- Przy $m_{\rm sp}\in(0.28, 0.40)$: klasycznie związany, kwantowo **niewiązany** (ZP wypycha)

**Znaczenie:** $C_{\rm Pl}$ leży w kwantowym oknie Efimova tylko dla **skończonego, nienastrojeniowego** zakresu $m_{\rm sp}$. To ostra predykcja falsyfikowalna: zmierzenie $m_{\rm sp}$ spoza $(0.076, 0.198)\,l_{\rm Pl}^{-1}$ wyklucza trójciałowe klastry Plancka.

---

### ★★ ex28 — Porównanie geometrii wielokątów N=3,4,5

**Pytanie:** Jak zależy energia wiązania od liczby cząstek N dla regularnych wielokątów (trójkąt, kwadrat, pentagon)?

**Metoda:** Dla każdego wielokąta N z bokiem d:
- ZP = N/(8d²), V2 = -C² × Σ_{par} e^{-md}/d, V3 = -C³ × Σ_{trojek} I_Y
- Dwuetapowy skan minimum energii (600+500 punktów)
- Pentagon: 5 par bocznych (d) + 5 przekątnych (φ·d), φ=2cos(π/5)≈1.618 (złoty podział)
- Pentagon: 10 trojek: 5 × typ A (d,d,φd) + 5 × typ B (d,φd,φd)

**Wyniki (m_sp=0.1, C=0.155):**

| Geometria | E_min(2B) [E_Pl] | E_min(3B) [E_Pl] | d_eq [l_Pl] | V3/V2 | Związany? |
|-----------|-----------------|-----------------|-------------|-------|-----------|
| Trójkąt N=3 | +0.000363 | −0.007063 | 4.517 | 1.505 | TAK (3B) |
| Kwadrat N=4 | +0.000445 | −0.068499 | 2.415 | 2.739 | TAK (3B) |
| Pentagon N=5 | +0.000355 | −0.226383 | 1.686 | 3.762 | TAK (3B) |

**Klasyczne progi wiązania (m_sp=0.1):**

| N | Geometria | C_crit(2B) | C_crit(3B) | Okno ΔC | C_Pl w oknie? |
|---|-----------|-----------|-----------|---------|--------------|
| 3 | Trójkąt | 0.1843 | 0.1276 | 0.0567 | nie |
| 4 | Kwadrat | 0.1657 | 0.1021 | 0.0636 | nie |
| 5 | Pentagon | 0.1587 | 0.0909 | 0.0678 | nie |

**Obserwacje:**
- Energia wiązania 3B **rośnie dramatycznie z N**: pentagon jest 32× silniej związany niż trójkąt przy C=0.155
- Równowagowe odległości **maleją z N**: 4.52 → 2.42 → 1.69 l_Pl — bardziej kompaktowe klastry
- Okno 3B-only przesuwa się ku **niższym C** wraz z rosnącym N (pentagon wiąże od C≈0.09 vs 0.13 dla trójkąta)
- Stosunek V3/V2 rośnie: 1.5 → 2.7 → 3.8 — siły 3-ciałowe dominują coraz bardziej
- **C_Pl=0.282 jest poza klasycznym oknem dla wszystkich N** przy m_sp=0.1 (para wiąże przy C≈0.16 dla N=5)
- Całka Feynmana I_Y dla pentagonalnego trojek (φd≈1.686 l_Pl) jest duża: ~10 dla typ A, ~8.8 dla typ B

**Plik:** `examples/ex28_polygon_geometry.py`

---

### ★★ ex29 — Funkcja falowa przy C=C_Pl (KLUCZOWY WYNIK)

**Pytanie:** Jak wygląda kwantowy stan podstawowy trójcząstkowego klastra TGP przy C=C_Pl? Czy jest dobrze zlokalizowany, czy rozproszony?

**Metoda:** FD Schrödingera (ex26), N=3000 punktów, d∈[0.4,30] l_Pl, C=C_Pl=0.282, m_sp=0.1.

**Wyniki przy m_sp=0.1, C=C_Pl:**

| Wielkość | Wartość |
|----------|---------|
| E_0 | −0.06964 E_Pl = −1.362×10⁸ J |
| d_peak (maks. gęstości) | 3.51 l_Pl |
| ⟨d⟩ (wartość oczekiwana) | 4.52 l_Pl |
| σ_d (szerokość) | 2.10 l_Pl |
| σ/⟨d⟩ | 0.464 (bardzo szeroki!) |
| d_L (klasyczny lewy TP) | 0.73 l_Pl |
| d_R (klasyczny prawy TP) | 6.03 l_Pl |
| L_well = d_R−d_L | 5.29 l_Pl |
| λ_dB/L_well | **3.03** (WKB NIEWAŻNE) |
| Frakcja tunelowania | **20.8%** |
| d_cl (klasyczne equilibrium) | 1.37 l_Pl |
| ⟨d⟩/d_cl | **3.30** (kwantowy stan ~3× dalej) |
| Liczba stanów związanych | **2** |

**Kluczowe obserwacje:**
- **⟨d⟩ = 4.52 l_Pl >> d_cl = 1.37 l_Pl** — funkcja falowa rozciąga się ~3× dalej niż klasyczne minimum
- **λ_dB/L_well = 3.03 >> 1** — WKB kompletnie nieważna, wymagana dokładna diagonalizacja FD
- **20.8% tunelowania** — znaczna część prawdopodobieństwa poza klasycznym obszarem
- **Drugi stan związany** n=1 przy E=-0.00065 E_Pl, d_peak=13.3 l_Pl — bardzo rozproszony
- **V3/V2 przy ⟨d⟩**: V2=-0.034, V3=-0.092, ratio=2.72 — siły 3-ciałowe dominują

**Trend z m_sp:**

| m_sp | λ | E_0 [E_Pl] | d_cl | ⟨d⟩ | σ | Tunel% |
|------|---|-----------|------|-----|---|--------|
| 0.08 | 12.5 | −0.107 | 1.36 | 4.17 | 1.87 | 19.6% |
| 0.10 | 10.0 | −0.070 | 1.36 | 4.52 | 2.10 | 20.8% |
| 0.12 | 8.3 | −0.043 | 1.40 | 4.98 | 2.40 | 22.5% |
| 0.15 | 6.7 | −0.018 | 1.44 | 6.04 | 3.17 | 26.3% |
| 0.18 | 5.6 | −0.004 | 1.48 | 8.05 | 4.47 | 32.1% |

Im bliżej granicy okna (m_sp→0.198), tym bardziej rozproszony stan (rosnące ⟨d⟩ i σ) — typowe dla stanów Efimova przy progu dysocjacji.

**Plik:** `examples/ex29_wavefunction.py`

---

### ★★ ex30 — Stabilność kształtu trójkąta (Q5 — walidacja aproksymacji równobocznej)

**Pytanie:** Czy konfiguracja równoboczna d12=d13=d23=d jest prawdziwym minimum energetycznym w przestrzeni kształtów, czy tylko punktem siodłowym?

**Metoda:** Parametryzacja odchyleń: d12=d*(1+εa), d13=d*(1+εb), d23=d. Hesjan H=∂²E/∂εi∂εj w punkcie równobocznym (0,0). Analiza własności dla C=C_Pl, m_sp=0.1.

**Kluczowy wynik:**

| d [l_Pl] | λ₁ (miękki mod) | λ₂ (sztywny mod) | Typ |
|----------|----------------|----------------|-----|
| 1.37 (d_cl) | **−0.117** | +0.088 | Siodło |
| 1.84 | ≈ 0 | −0.030 | Przejście |
| 4.52 (⟨d⟩_q) | −0.064 | **−0.037** | Pełna niestabilność |

**Trójkąt równoboczny jest punktem siodłowym w przestrzeni kształtów dla wszystkich d!**

**Fizyczna interpretacja:**
- Miękki mod (εa=−εb): jedno ramię skraca się, drugie wydłuża → V2 zmniejsza się (konweksowość Yukawa)
- Pary są wiązane przy C=C_Pl (E0(2B)<0), więc odchylenia od równoboczności obniżają energię
- Prawdziwe minimum jest konfigurację równoramienną (isosceles), nie równoboczną
- Przy C<0.191 (próg kwantowy 3B): jeden eigenvalue jest jeszcze dodatni → typ siodłowy, nie pełna niestabilność

**Konsekwencja dla ex23-ex29:**
- Aproksymacja równoboczna daje **górne ograniczenie** na E_0 (nie dokładną wartość)
- Prawdziwa energia jest niższa → prawdziwe okno Efimova może być **szersze** niż (0.076, 0.198) l_Pl^-1
- Ex23-ex29 to wyniki **konserwatywne** — qualitative results unchanged
- Wymagany pełny kalkulator 2D FD w przestrzeni hyperradia dla ilościowych wyników

**Plik:** `examples/ex30_triangle_stability.py`

---

### ★★ ex31 — Krzywe rotacji TGP vs CDM/MOND (Q2)

**Pytanie:** Czy TGP z m_sp ~ kpc może wyjaśnić krzywe rotacji galaktyk?

**Metoda:** Model Drogi Mlecznej: dysk gwiazdowy (5×10¹⁰ M☉), bańka (10¹⁰ M☉), gaz (5×10⁹ M☉). Porównanie: 4 scenariusze (Newtona, TGP λ=1/5/15 kpc, CDM+NFW, MOND).

**Wyniki v_circ [km/s]:**

| r [kpc] | Baryony (Newton) | TGP λ=1 kpc | TGP λ=5 kpc | TGP λ=15 kpc | CDM+NFW | MOND |
|---------|-----------------|------------|------------|-------------|---------|------|
| 5  | 170 | **34** | 145 | 166 | 222 | 182 |
| 10 | 165 | **4** | 105 | 153 | 239 | 195 |
| 20 | 126 | **0** | 38 | 99 | 228 | 190 |
| 30 | 100 | **0** | 13 | 64 | 217 | 185 |

**Wniosek: TGP NIE może wyjaśnić płaskich krzywych rotacji.**

- **Yukawa screening ZAWSZE zmniejsza grawitację** dla r > λ_sp: f(r) = e^{-m_sp*r}(1+m_sp*r) ≤ 1
- Płaskie krzywe rotacji wymagają WIĘCEJ grawitacji na dużych r (nie mniej!)
- TGP jedzie w złym kierunku dla dowolnego m_sp > 0
- CDM (NFW) i MOND obie produkują płaskie krzywe (~220 km/s) — TGP nie
- Efekt 3-ciałowy dla materii zwykłej: C³~10⁻⁵⁷ → całkowicie pomijalny
- MOND Tully-Fisher: v⁴ = G·M_bar·a0 → v_TF = 179 km/s (vs obliczone 185 km/s ✓)

**Ogólna zasada:** Każdy skalar z m_sp² > 0 daje potencjał ekranowany — nigdy nie można uzyskać WIĘCEJ grawitacji niż Newton. Tylko skalar tachioniczny (m_sp²<0) mógłby dać wzmocnienie na dużych skalach — poza zakresem TGP Drogi B.

**Plik:** `examples/ex31_rotation_curve.py`

---

## 6. Stan falsyfikowalności TGP (Droga B)

**TAK, TGP Drogi B jest falsyfikowalna** przez:

1. Anomalie $G_{eff}(r)$ na skali $\lambda = 1/m_{sp}$ (mierzalne satelitami)
2. Progi częstotliwości w promieniowaniu skalarnym (jeśli $m_{sp} < \omega_{orbital}$)
3. Ścisła proporcjonalność $C_i \propto m_i$ — jakiekolwiek odchylenie falsyfikuje teorię
4. Efekty wielociałowe (N_crit) — mierzalne jeśli $m_{sp} \sim$ skala Plancka

**Jedyne co jest potrzebne:** jeden dodatkowy pomiar wyznaczający $m_{sp}$.

---

### ★★ ex32 — Klasyczne minimum isoceliczne (Q5 — singularność kolapsowa)

**Pytanie:** Czy prawdziwe klasyczne minimum energetyczne trójciała leży na linii równobocznej, czy poza nią (konfiguracja isoceliczna)?

**Metoda:** Skan 2D po siatce (a,b) z warunkiem geometrii isocelicznej d12=d13=a, d23=b. Optymalizacja L-BFGS-B z ograniczeniami (b>0.1, b<2a). C=C_Pl, m_sp skanowane.

**Kluczowy wynik:**

| m_sp | E_eq [E_Pl] | E_iso [E_Pl] | dE [E_Pl] | dE% | b/a |
|------|------------|-------------|----------|-----|-----|
| 0.076 | -0.334 | -0.512 | -0.178 | -53% | 0.131 |
| 0.100 | -0.265 | -0.446 | -0.181 | -68% | 0.131 |
| 0.130 | -0.202 | -0.385 | -0.184 | -91% | 0.131 |
| 0.198 | -0.108 | -0.300 | -0.192 | -177% | 0.108 |

**Wniosek: Minimum isoceliczne dąży do b→0 (singularność kolapsowa).**

- Potencjał Yukawa V2 = -C²e^{-mb}/b → -∞ gdy b→0
- Term ZP = 3/(8*d_bar²) jest skończony dla b→0 z a=const
- Zatem: klasyczna energia efektywna jest **nieograniczona od dołu** w kierunku isocelicznym
- Minimum isoceliczne at b=0.3 (dolna granica siatki) z E_iso ≈ -0.446 vs E_eq=-0.265
- **To singularność klasyczna, nie fizyczne minimum**

**Konsekwencja:** Mechanika kwantowa jest wymagana — operator energii kinetycznej ∇² daje barierę kwantową (nieskończona energia kinetyczna gdy b→0). Właściwy rachunek wymaga pełnego 2D FD Schrödingera w koordynatach Jacobi.

**Plik:** 

---

### ★★ ex33 — Pełny 2D Schrödinger w koordynatach Jacobi (Q5)

**Pytanie:** Jaka jest prawdziwa energia stanu podstawowego trójciała w podprzestrzeni isocelicznej, użując poprawnego operatora energii kinetycznej w koordynatach Jacobi?

**Metoda:**
- Koordynaty: b=d23, h=√(a²-b²/4) (wysokość trójkąta od cząstki 1 do podstawy 23)
- Zredukowane masy: μ_b = m/2 = 0.5, μ_h = 2m/3 ≈ 0.667
- Operator kinetyczny: T = -∂²/∂b² - (3/4)∂²/∂h² (BEZ dodatkowego termu ZP)
- Siatka N_b=50, N_h=50, b,h∈[0.15,20] l_Pl
- Rzadka macierz 2500×2500, diagonalizacja 

**Kluczowy wynik (C=C_Pl, m_sp=0.1):**

| Metoda | E_0 [E_Pl] | ⟨b⟩ | ⟨a⟩ | b/a (max) |
|--------|-----------|-----|-----|-----------|
| 1D równoboczna (ex23-ex29) | **-0.0698** | -- | 4.52 | 1.000 |
| 2D Jacobi isoceliczna (ex33) | **-0.0088** | 6.88 | 7.53 | 0.956 |
| Różnica ΔE | +0.0610 (+87%) | | | |

**Scan vs m_sp (N=40×40):**

| m_sp | E0_1D [E_Pl] | E0_2D [E_Pl] | ΔE | Związany 2D? |
|------|-----------|-----------|----|-------------|
| 0.076 | -0.117 | -0.049 | +0.068 | TAK |
| 0.100 | -0.070 | -0.013 | +0.057 | TAK |
| 0.130 | -0.034 | +0.010 | +0.044 | NIE |
| 0.150 | -0.018 | +0.019 | +0.037 | NIE |
| 0.198 | ≈0 | +0.029 | +0.029 | NIE |

**Kluczowe wnioski:**

1. **E0_2D >> E0_1D** (czynnik ~8) — obliczenie 2D Jacobi daje mniej głębokie minimum
2. **Przyczyna: różnica mas efektywnych**:
   - Ex23-ex29 używa μ_eff = 1 m_Pl (współczynnik kinematyczny 1/2) + ZP=3/(8d²)
   - Koordynaty Jacobi: wzdłuż modu oddechowego μ_eff = 8/25 = 0.32 (współczynnik 1.5625)
   - Czynnik 3.1× większa energia kinetyczna → 1D obliczenie zawyżało związanie
3. **Kształt funkcji falowej**: b/a=0.956 w maksimum — blisko równobocznego (1.000)
4. **Nowe okno Efimova (2D)**: m_sp ∈ (0, ~0.12) l_Pl^-1 zamiast (0.076, 0.198)
5. **Aproksymacja równoboczna** w ex23-ex29: ZAWYŻAŁA związanie (E0_1D zbyt ujemne)

**Interpretacja fizyczna:** μ=1 w ex23 jest wyborem empirycznym, nie wynikającym z koordynatów Jacobi. Właściwa masa efektywna dla modu oddechowego jest 3× mniejsza → energia kinetyczna 3× większa → stan podstawowy wyżej (mniej ujemny).

**Status Q5:** ROZWIĄZANY częściowo. Pełna 3D kalkulacja (cząstki 3D, wszystkie mody kątowe) wymagana dla dokładnych wyników Efimova.

**Plik:** 

---

## Status pytań otwartych

| Q | Pytanie | Status |
|---|---------|--------|
| Q1 | Derywacja sprzężenia Yukawa | OMÓWIONE — dilaton g=1/m_Pl, Lagrangian otwarty |
| Q2 | Ciemna materia / krzywe rotacji | **ROZWIĄZANE** (negatywnie) — TGP reducuje grawitację |
| Q3 | Droga C: defekt → Yukawa | OMÓWIONE — ogon oscylacyjny, potrzebna zmiana Lagrangianu |
| Q4 | TGP vs GR | OMÓWIONE — koegzystencja, BD przy m_sp Plancka |
| Q5 | Pełne 2D kwantowe Efimova | **POSTĘP** — ex32: singularność kolapsowa; ex33: 2D Jacobi daje E0=-0.009 E_Pl |
| Q9 | Geometria N=5 (pentagon) | **ROZWIĄZANE** — ex28, V3/V2=3.8, E_bind×32 vs N=3 |
| Q11 | Dynamika kwantowa przy C_Pl | **ROZWIĄZANE** — ex29, 2 stany, ⟨d⟩=4.5 l_Pl |
| Q12 | Fizyczna interpretacja granic okna | **ROZWIĄZANE** — górna: pary zaczynają wiązać; dolna: ZP wypycha |
| Q13 | Niesprzeczność okna kl./kw. | **ROZWIĄZANE** — przesunięcie ΔC≈0.063 od operatora ZP |
| Q14 | Naprawa ex23 | **ROZWIĄZANE** — dwuetapowy skan |
