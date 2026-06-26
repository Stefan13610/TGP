---
title: "Phase 3 — op-Csigma-coarse-graining: normalizacja σ_0 + rozstrzygnięcie C_σ vs σ_0. WYNIK KLUCZOWY: C_σ i σ_0 NIE są dwoma niezależnymi parametrami — istnieje DOKŁADNA redundancja przeskalowania pola σ→λσ (sympy 3/3 niezmienniki), więc fizyczny jest TYLKO iloczyn T=C_σσ_0² (tensor stiffness). To uzasadnia redukcję rem:param-counting 3→2. κ_E=8πG_0T/c³ przez JEDEN parametr; survival ⟺ T=5/6·(c³/8πG_0) miara zero; naturalna T=c³/8πG_0 (grawiton GR)→κ_E=1→FALSIFIED. WERDYKT: F-CG-C RESOLVED (param-counting)+PARTIAL (wartość); lean FALSIFIED."
type: phase_result
status: PHASE3_COMPLETE (PARTIAL; param-counting RESOLVED)
phase: 3
cycle: op-Csigma-coarse-graining-2026-06-14
created_date: 2026-06-14
authorization: "User 2026-06-14 (sesja #31): 'Działaj z fazą 3'"
sympy_script: "[[./Phase3_sympy.py]]"
sympy_output: "[[./Phase3_sympy.txt]]"
sympy: "PASS: redundancja σ→λσ, σ_0→λσ_0, C_σ→C_σ/λ² zostawia h^TT, kinetyk i flux niezmiennicze (3/3); T=C_σσ_0² jedyny fizyczny; κ_E=8πG_0T/c³; T_survive/T_natural=5/6"
flag_F_CG_C: "RESOLVED (param-counting: C_σ,σ_0 → 1 parametr T via field-rescaling redundancy) + PARTIAL (wartość T = GAP)"
flag_F_CG_D: "PARTIAL-lean-FALSIFIED (na 1 parametrze T; survival miara zero, naturalna κ_E=1→7/6 FALSIFIED)"
anti_lakatos_lock: PRESERVED
---

# Phase 3 — normalizacja $\sigma_0$ i rozstrzygnięcie $C_\sigma$ vs $\sigma_0$

## §0 — Werdykt fazy w skrócie
| Pytanie | Wynik |
|---|---|
| Czym jest $\sigma_0$? | Normalizacją korelatora kierunkowego (skala tła $\sigma_{ab}$, ~ $\Phi_0$ tego samego ŝ); **konwencja** |
| $C_\sigma$ vs $\sigma_0$: 1 czy 2 parametry? | **1** — istnieje **dokładna redundancja przeskalowania** $\sigma\to\lambda\sigma$ (sympy 3/3); fizyczny tylko $T\equiv C_\sigma\sigma_0^2$ |
| $\kappa_E$? | $\kappa_E=8\pi G_0\,T/c^3$ — przez **jeden** parametr (tensor stiffness $T$) |
| **F-CG-C** | **RESOLVED** (param-counting; uzasadnia rem:param-counting 3→2) **+ PARTIAL** (wartość $T$ = GAP) |
| **F-CG-D** | **PARTIAL-lean-FALSIFIED** — survival $T=\tfrac56\cdot\tfrac{c^3}{8\pi G_0}$ miara zero; naturalna $T=\tfrac{c^3}{8\pi G_0}\Rightarrow\kappa_E=1\Rightarrow7/6$ FALSIFIED |

## §1 — Identyfikacja $\sigma_0$
Z eq:h-TT-from-sigma: $h^{TT}_{ij}=2\sigma_{ij}/\sigma_0$. Zatem $\sigma_0$ jest **normalizacją**, która przeprowadza substratowy korelator kierunkowy $\sigma_{ab}\sim\langle\hat s\hat s\rangle$ (wymiar $\sim\Phi_0$) w bezwymiarową perturbację metryki $h^{TT}$ — dokładny analog $\Phi_0$ w module skalarnym ($U=\delta\Phi/\Phi_0$). Naturalna skala: $\sigma_0\sim\Phi_0$ (ten sam ŝ, ta sama amplituda VEV), z O(1) niepewnością. **Ale** — jak pokazuje §2 — wartość $\sigma_0$ osobno jest **konwencją**, nie fizyką.

## §2 — Wynik kluczowy: $C_\sigma$ i $\sigma_0$ to JEDEN parametr (redundancja przeskalowania)
Rozważmy przeskalowanie pola $\sigma_{ab}\to\lambda\sigma_{ab}$, z towarzyszącym $\sigma_0\to\lambda\sigma_0$ i $C_\sigma\to C_\sigma/\lambda^2$. Wszystkie trzy wielkości fizyczne są **niezmiennicze** (sympy, `Phase3_sympy.txt`, 3/3 PASS):
$$h^{TT}=\frac{2\sigma}{\sigma_0}\ \checkmark,\qquad \frac{C_\sigma}{2}(\partial\sigma)^2\ \checkmark,\qquad \text{flux}\sim C_\sigma\sigma_0^2\ \checkmark.$$
> **Twierdzenie (Phase 3):** $C_\sigma$ i $\sigma_0$ **nie są osobno fizyczne** — różnią się o gauge przeskalowania $\lambda$. Jedyny niezmienniczy (fizyczny) parametr radiacyjny sektora tensorowego to
> $$\boxed{\ T\equiv C_\sigma\,\sigma_0^2\quad(\text{„tensor stiffness"})\ }$$

To **rozstrzyga i uzasadnia** hedge rdzenia `rem:sigma-params` („$C_\sigma$ lub równoważnie $\sigma_0$") oraz redukcję `rem:param-counting` **3→2**: $C_\sigma$ i $\sigma_0$ liczone osobno były pozorną nadmiarowością; faktycznie to jeden DOF $T$. (Sektor skalarny: 2 parametry; tensorowy dodaje **1** = $T$; razem **2+1=3**, a redukcja do 2 wymaga policzenia *wartości* $T$ — nadal otwarte, §4.)

## §3 — Konsystencja konwencji normalizacyjnych (uczciwe bookkeeping)
Znalazłem **pozorną rozbieżność** w rdzeniu vs rodzicu, którą redundancja §2 wyjaśnia:
- **Core thm:amplitude-matching:** $\xi_{\rm eff}=4\pi G_0\sigma_0\Phi_0$ — amplituda $\propto\xi_{\rm eff}/\sigma_0$, **bez** $C_\sigma$.
- **Parent det J:** $O_{\rm amp}\sim\xi_{\rm eff}/(C_\sigma\sigma_0)$ — **z** $C_\sigma$.

**Rozwiązanie:** to ta sama fizyka w różnych normalizacjach pola. Core używa normalizacji **kanonicznej** ($\tilde\sigma=\sqrt{C_\sigma}\,\sigma$, $C_\sigma$ wchłonięte) — wtedy $C_\sigma$ znika z amplitudy, a h-normalizacja staje się $\sigma_0\sqrt{C_\sigma}$. Parent trzyma jawne $C_\sigma$. Pod redundancją §2 $O_{\rm amp}=\xi_{\rm eff}/(C_\sigma\sigma_0)$ jest niezmiennicza, gdy $\xi_{\rm eff}\to\xi_{\rm eff}/\lambda$ (sympy §2). Amplitude-matching **slave'uje** $\xi_{\rm eff}$ do $(C_\sigma\sigma_0)$; po matchingu jedyną swobodą radiacyjną zostaje $T$. **Brak realnej sprzeczności** — flaga do ujednolicenia notacji w rdzeniu (sugerowana edycja w §6).

## §4 — Złożenie $\kappa_E$ przez jeden parametr $T$ (F-CG-D)
$$\boxed{\ \kappa_E=\frac{8\pi G_0\,T}{c^3},\qquad T=C_\sigma\sigma_0^2\ }$$
Test value-blind (próg LOCKED Phase 0), sympy:
- **survival** ⟺ $\kappa_E=\tfrac56$ ⟺ $T=\tfrac{5}{6}\cdot\tfrac{c^3}{8\pi G_0}=\tfrac{5c^3}{48\pi G_0}$ (EXACT, miara zero);
- **naturalna** $\kappa_E=1$ ⟺ $T=\tfrac{c^3}{8\pi G_0}$ = **sztywność grawitonu GR** ⟹ total $7/6$ = gałąź B PR-025 (FALSIFIED);
- stosunek $T_{\rm survive}/T_{\rm natural}=5/6$.

## §5 — Czy $T$ jest zlockowane do grawitonu GR? (NIE — brak diff-inwariancji)
- **GR:** diffeomorfizm-inwariancja ⟹ jedno $G$ lockuje amplitudę i strumień ⟹ $\kappa_E=1$ **wymuszone** (graviton ma jedyną dozwoloną sztywność $c^3/8\pi G_0$).
- **TGP:** $\sigma_{ab}$ to **osobny, drugi DOF sieci** (aksjomat A1), **NIE** część diff-inwariantnego $G_{\mu\nu}$. $\det J=-\xi_{\rm eff}/C_\sigma\neq0$ (amplituda ⊥ strumień). **Brak tożsamości Warda** lockującej $T$ do $c^3/8\pi G_0$.

> **Wniosek §5:** „jeden parametr $T$" — TAK (redundancja §2), ale jego **wartość nie jest chroniona** żadną symetrią. $T=\tfrac{c^3}{8\pi G_0}$ ($\kappa_E=1$) jest „naturalnym" punktem (substratowe skale O(1)), nie wymuszonym; $T=\tfrac56\cdot\tfrac{c^3}{8\pi G_0}$ ($\kappa_E=5/6$) to **miara zero**. ⟹ lean **FALSIFIED**, niezmieniony, teraz na czystym jednym parametrze.

## §6 — Sugerowana spójność rdzenia (NIE wykonana — czeka na decyzję)
`rem:sigma-params` i `rem:param-counting` warto doprecyzować: „$C_\sigma$ i $\sigma_0$ wchodzą w obserwable
radiacyjne **wyłącznie** przez $T=C_\sigma\sigma_0^2$ (redundancja przeskalowania $\sigma\to\lambda\sigma$);
to jeden DOF, nie dwa". Oraz ujednolicić amplitude-matching (kanoniczna vs jawne $C_\sigma$, §3). **Nie
edytuję rdzenia w tym cyklu** (forbidden move #3: zakaz rewizji LOCKED bez osobnej autoryzacji) — zgłaszam jako rekomendację.

## §7 — Anti-Lakatos (checklist Phase 3)
✓ Redundancja przeskalowania **zweryfikowana sympy** (3/3 niezmienniki) — twierdzenie, nie zgadywanie.
✓ **Wartość $T$ NIE policzona** — pozostaje GAP (Phase 2 §9); zero strojenia do 5/6.
✓ Rozbieżność konwencji rdzeń/rodzic ujawniona i rozwiązana uczciwie (§3), nie zamieciona.
✓ Lean FALSIFIED oznaczony jako strukturalny (nie-twardy): $\det J\neq0$ ⟹ $T$ formalnie wolne.
✓ Rdzeń NIE edytowany (forbidden move #3); rekomendacja zgłoszona osobno (§6). Budżet stałych 0.

## §8 — Handoff do Phase FINAL
Po Phase 3 stan: sektor radiacyjny ma **dokładnie jeden** fizyczny parametr $T=C_\sigma\sigma_0^2$, z
$\kappa_E=8\pi G_0T/c^3$. **Phase FINAL** zbierze werdykt agregatowy (F-CG-E):
- F-CG-A = PASS-CONDITIONAL, F-CG-B = PARTIAL, F-CG-C = RESOLVED+PARTIAL, F-CG-D = PARTIAL-lean-FALSIFIED;
- agregat = **UNDERDETERMINED-fine-tuned (strukturalnie oczyszczony do 1 parametru $T$), lean FALSIFIED**;
- **residual GAP** do twardego werdyktu: liczbowa $T$ z **lattice-MC kierunkowego bąbla** $\langle O_{ab}O_{cd}\rangle$ na sieci 3D Ising (jedyna droga od PARTIAL do DERIVED/FALSIFIED). Phase FINAL wymaga „działaj".
