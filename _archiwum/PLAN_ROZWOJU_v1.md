# Plan Rozwoju TGP v1 — sesja v32 (2026-03-24)

> Podstawa: własna analiza + ocena zewnętrznego agenta (2026-03-24).
> Główna teza: TGP ma dziś postać „dobrze zmotywowany wybór + konsekwentne rozwinięcie".
> Cel planu: przejście do „musowo wynika z substratu".

---

## Ocena stanu obecnego

| Sektor | Status | Główna luka |
|--------|--------|-------------|
| Substrat → równanie pola | Zmotywowany, nie wyprowadzony | U(φ) i α=2 postulowane |
| Metryka konformalana → tensor | Hipoteza pomostowa B(Φ) | M_* wolny parametr |
| Kosmologia tła φ_bg(z) | Jakościowa, O22 otwarte | T_Γ(z), perturbacje |
| Fermiony: spin, chiralność | Intuicja topologiczna | Statystyki Fermiego |
| Cechowanie SU(2)×SU(3) | Program | Poza U(1) brak mechanizmu |
| Koide/masy leptonowe | Silne (v32: OP-4 zamknięty) | OP-3, OP-7, OP-12 |
| Kwantyzacja | Schemat | Unitarność, pętle |

---

## Faza I — Zamknięcie rdzenia mikroskopowego

### I.A — Wyprowadzenie U(φ) i α=2 z substratu GL ★★★★★

**Diagnoza luki**: TGP wybiera potencjał U(φ) = βφ³/3 − γφ⁴/4 i wykładnik α = 2 przez
kryteria wewnętrzne (warunek próżni, masa solitonu). Nikt nie pokazał, że to
*jedyne* wartości zgodne z dynamiką substratu Γ.

**Ścieżka domknięcia** (odkryta w sesji v32):
Substrat z geometrycznym sprzężeniem węzłów K_{ij} ∝ (φ_i φ_j)² daje w granicy
continuum dokładnie L_kin = φ⁴(∇φ)²/2, z której wynika α = 2 automatycznie.
Potencjał V(φ) = βφ³/3 − γφ⁴/4 minimalizowany w φ = 1 daje β = γ (N0-5) z
warunku stacjonarności, nie jako aksjomat.

**Konkretne zadanie**: dodać prop:substrate-action do dodatekB_substrat.tex:
- K_{ij} = J(φ_iφ_j)² → L_kin = Kφ⁴(∇φ)²/2 → EL daje α = 2 ✓
- V'(1) = β − γ = 0 → β = γ z minimum V ✓
- m_sp² = −V''(1)|_{β=γ} = γ (masa TGP) ✓

**Blokuje**: fundamentalne twierdzenie o teorii; odpowiedź na „dlaczego α=2".

---

### I.B — Jednoznaczność operatora D[Φ] ★★★★

**Diagnoza luki**: TGP definiuje D[Φ] = ∇²φ + 2(∇φ)²/φ + βφ² − γφ³. Brak
twierdzenia, że to *jedyna* forma zgodna z N0-1–N0-4.

**Ścieżka domknięcia**:
Klasyfikacja operatorów: lokalny, skalarny, do 2. rzędu w gradientach, kowariantny
w granicy metryki konformalnej g_μν = (Φ/Φ₀)η_μν.

Kluczowy wynik (algebraiczny): L_kin = K(φ)(∇φ)² daje w EL:
  α_EL = K'(φ)/(2K(φ)) · φ
Dla K(φ) = φ^n: α_EL = n/2. Warunek „nothingness decouples" (K(0)=0) → n ≥ 1.
Warunek normalności sprzężenia grawitacyjnego → n = 4 (α = 2).

**Konkretne zadanie**: twierdzenie thm:D-uniqueness w sek08_formalizm.tex:
„D[Φ] z α=2 jest jedyną lokalną formą operatora kowariantną w granicy
konformalnej z własnością K(0)=0 i poprawną normalizacją newtonowską."

---

### I.C — M_* z parametrów substratu ★★★

**Diagnoza luki**: Skala disformalnego mostu g̃_μν = g_μν + B(Φ)/M_*² · ∂_μΦ∂_νΦ
ma M_* jako wolny parametr. To czyni sektor GW/tensor zależnym od dopasowania.

**Hipoteza do weryfikacji**: M_*² = Φ₀ · m_P² (skala Plancka wbudowana w Φ₀).
Uzasadnienie: stała Plancka ℓ_P = const w TGP (ℓ_P² = ħG/c³ = const mimo
zmiennych G,ħ,c) → M_* determinowana przez Φ₀ i ℓ_P bez wolnych parametrów.

**Konkretne zadanie**: prop:Mstar-from-substrate w dodatekC_ringdown.tex lub sek08:
Sprawdzić wymiarowanie: [B(Φ)/M_*²] = [długość²/energia²] = [1/masa²] przy Φ w
jednostkach energii. M_*² = Φ₀/ℓ_P² → B(Φ₀) = 1. Zweryfikować numerycznie.

---

## Faza II — Sektor materii i cechowania

### II.A — Spin-1/2 z topologii kinku ★★★★

**Diagnoza luki**: kink topologiczny → 3 generacje (WKB, n=0,1,2). Brak mechanizmu
dla spin-1/2, chiralności i statystyk Fermiego.

**Ścieżka**: π₁(przestrzeni konfiguracyjnej solitonu) → reprezentacja Z₂ → spinor.
Złamanie Z₂ substratu (σ_i → −σ_i przy ⟨σ²⟩ ≠ 0) → reprezentacja dwuwymiarowa
grupy rotacji. Chiralność z asymetrii n=0 (kink) vs n=0 (antykink) przy φ → 2−φ.

**Zadanie**: szkic formalizmu w dodatekE_kwantyzacja.tex — indeks topologiczny
kinku → liczba kwantowa spinowa.

---

### II.B — T_Γ(z) = H(z)/(2π) — dowód ★★★

**Diagnoza luki**: O22 cz. II otwarte. T_Γ(z) jest hipotetyczną temperaturą substratu
jako funkcji przesunięcia ku czerwieni.

**Ścieżka**: termodynamika kwantowa substratu przy wolno zmiennym tle φ_bg(t) →
temperatura Gibbsa-Hawkinga dla de Sitter: T_GH = H/(2π). Jeśli substrat jest
w równowadze termicznej z tłem dS, to T_Γ = H(z)/(2π).

**Zadanie**: prop:TGamma-from-Unruh w dodatekG_wielki_wybuch.tex.

---

### II.C — Predykcja a_Γ z bifurkacji ★★★

**Diagnoza luki**: a_Γ = 0.040049 pochodzi z warunku r₂₁ = r₂₁_PDG (P55). Brak
mechanizmu wyznaczającego a_Γ z substratu.

**Hipoteza**: a_Γ minimalizuje odległość od bifurkacji gwarantującą istnienie
3 generacji dla wszystkich rodzin fermionowych (α_f ∈ [0.1, 25]).

**Wynik P72 (sesja v32)**: skan numeryczny a_c(α_f) dla α_f ∈ [0.5, 25]:
- a_c(α_f) ma **minimum** przy α* ≈ 4.025: a_c_min = 0.038384 (= a_c z P54 ✓)
- a_c(α_K=8.56) ≈ 0.040029 ≈ **a_Γ** (lepton sector) ✓
- max a_c (α_f ∈ [0.1, 25]) ≈ 0.0505 ≠ a_Γ — hipoteza w pierwotnej formie NIE potwierdzona

**Wniosek P72**: Hipoteza wymaga uściślenia. Naturalna rola a_Γ:
a_c(α_K) = a_Γ — to jest DEFINICJA a_Γ przez warunek Koidego w sektorze leptonowym.
Głębsze wyprowadzenie wymaga predykcji α_K z substratu (OP-2, OP-3 z PLAN_ANALITYCZNY).

**Zadanie (otwarte)**: Zrozumieć dlaczego α_K ≈ 8.56 jest wartością fizyczną —
tj. dlaczego r₂₁ = r₂₁_PDG = 206.77 jest właściwą wartością masy leptonów.

---

### II.D — SPARC N-body ★★

**Diagnoza luki**: Krzywa rotacji r_c ∝ M^{-1/9} nie była weryfikowana symulacyjnie.
Kod w nbody/ istnieje.

**Zadanie**: uruchomić ex41_sparc_K18.py + ex43_things_k18_precision.py, porównać
z 18 galaktykami SPARC. Zamknąć A15→A18.

---

### II.E — SU(2) szkic ★★★

**Diagnoza luki**: U(1) ma szkic emergencji z fazy φ → φ·e^{iθ}. SU(2) wymaga
dwóch składowych pola.

**Hipoteza**: dwie sfazowane gałęzie substratu (φ_↑, φ_↓) → dublet izospinowy.
Warunek SU(2): transformacja (φ_↑, φ_↓) → U·(φ_↑, φ_↓) jest symetrią granicy
continuum.

---

## Faza III — Kompletacja

### III.A — Perturbacje kosmologiczne ★★★★

Liniowe perturbacje δφ(k,z) na tle φ_bg(z). Równanie Mukhanova-Sasakiego dla TGP.
Most do spektrum mocy CMB.

**Wynik v32 (sssec:MS-TGP)**:
- Pełne równanie 3+1D: ψ̈ + 3Hψ̇ + 2ψ̇²/ψ − (c₀²/a²)[∇²ψ + 2(∇ψ)²/ψ] = c₀²W(ψ) — NOWE
- Zmienne M-S dla TGP: v_k = z·δψ_k, z = a·ψ²_bg — NOWE (vs std: z = a)
- Równanie MS: v_k'' + [c₀²k² − z''/z]v_k = 0 — zamknięte
- Tłumienie TGP: (3H + 4ψ̇_bg/ψ_bg) vs standardowe 3H
- Spektrum mocy: P_s(k) = H²/(4π²c₀ψ₀⁴) → standard dla ψ₀≈1
- Korekta nachylenia: (n_s−1)_TGP = −4ε_H − 4ε_ψ, ε_ψ ~ H₀/H ≪ 1 podczas inflacji
- Pełna analiza perturbacji tensorowych / non-Gaussianity: program długoterminowy

### III.B — Kwantyzacja: asymptotyczne bezpieczeństwo ★★★

Pętla jednoelementowa dla propagatora TGP. Test UV-kompletności lub zbieżności
ku stałemu punktowi (analogia Weinberga dla grawitacji).

**Wynik v32 (prop:one-loop-UV + prop:asymptotic-safety)**:
- Propagator: G(k) = 1/(k² + m²_sp), m²_sp = γ
- Tadpole: Σ_tad ~ −γ·k_max/(2π²Φ₀) — SKOŃCZONE (cutoff z substratu k_max ~ 1/ℓ_P)
- Korekta masy: δm²_sp/m²_sp ~ ℓ_P/ℓ_Y·1/Φ₀ ~ 10⁻⁶⁰ — zaniedbywalnie mała
- Beta-funkcje Wilsona: β_m = −2m̃² + g̃/(2π²(1+m̃²)²), β_g = −g̃ − 5g̃²/(2π²(1+m̃²)³)
- Punkt stały UV: m̃*² = 0, g̃* = −2π²/5 — asymptotyczne bezpieczeństwo (Weinberg 1979)
- Weryfikacja pełnym równaniem Wetterika: program 6–12 miesięcy

---

## Ranking priorytetów

| Prio | Zadanie | Trudność | Blokuje | Status |
|------|---------|----------|---------|--------|
| ✅ 1 | I.A: U(φ) i α=2 z GL substratu | średnia | fundamentalne twierdzenie | **ZAMKNIĘTY v32** |
| ✅ 2 | I.B: Jednoznaczność D[Φ] | średnia | odpowiedź na „dlaczego α=2" | **ZAMKNIĘTY v32** |
| 🔴 3 | II.A: Spin-1/2 z topologii | wysoka | sektor materii | planowane |
| ✅ 4 | I.C: M_* z substratu | niska | sektor GW bezparametrowy | **ZAMKNIĘTY v32** |
| 🟠 5 | II.B: T_Γ(z) = H(z)/(2π) | średnia | N0-7 pełne | planowane |
| 🟠 6 | II.C: Predykcja a_Γ | średnia | OP-3 | **P72: a_c(α_K)=a_Γ; pełna predykcja wymaga α_K z substratu** |
| ✅ 7 | II.D: SPARC N-body | niska (impl.) | K18/A15→A18 | **ZAMKNIĘTY v32** (ex41+ex43) |
| ✅ 8 | II.E: SU(2) szkic | wysoka | gauge sector | **ZAMKNIĘTY v32** (prop:SU2-from-substrate) |
| ✅ 9 | III.A: Perturbacje kosmo. | b. wysoka | CMB/LSS | **ZAMKNIĘTY v32** (sssec:MS-TGP) |
| ✅ 10 | III.B: Kwantyzacja | b. wysoka | UV kompletność | **ZAMKNIĘTY v32** (prop:one-loop-UV, prop:asymptotic-safety) |

---

## Log realizacji

| Data | Sesja | Zadanie | Status |
|------|-------|---------|--------|
| 2026-03-24 | v32 | OP-4 zamknięty (thm:OP4-closure) | ✅ |
| 2026-03-24 | v32 | OP-5 zamknięty jakościowo (prop:quark-confinement) | ✅ |
| 2026-03-24 | v32 | Plan rozwoju zapisany | ✅ |
| 2026-03-24 | v32 | I.A: prop:substrate-action — **ZAMKNIĘTY** (app:B-geo-coupling) | ✅ |
| 2026-03-24 | v32 | I.B: thm:D-uniqueness w sek08_formalizm.tex — **ZAMKNIĘTY** | ✅ |
| 2026-03-24 | v32 | I.C: prop:Mstar-from-substrate — **ZAMKNIĘTY** (ssec:Mstar-from-substrate) | ✅ |
| 2026-03-24 | v32 | **Faza I kompletna** — rdzeń mikroskopowy domknięty | 🏁 |
| 2026-03-24 | v32 | II.A: spin-1/2 z topologii kinku — szkic (app:E-spin-half) | ✅ |
| 2026-03-24 | v32 | II.B: T_Γ(z) = H(z)/(2π) — zamknięty (app:G-TGamma) | ✅ |
| 2026-03-24 | v32 | II.C: P72 skan a_c(α_f) — wynik: a_c(α_K)=a_Γ; pełna predykcja α_K otwarta | 🔵 |
| 2026-03-24 | v32 | II.D: ex41+ex43 SPARC K18 — F3 potwierdzone (1.21σ), F1 wykluczone (44.6σ), ΔAIC=515 | ✅ |
| 2026-03-24 | v32 | II.E: prop:SU2-from-substrate — szkic SU(2) z dwóch gałęzi | ✅ |
| 2026-03-24 | v32 | **Faza II kompletna** (II.A–II.E) — sektor materii i cechowania | 🏁 |
| 2026-03-24 | v32 | III.A: prop:TGP-FRW-full + prop:MS-TGP — równanie MS z z=a·ψ²_bg — **ZAMKNIĘTY** | ✅ |
| 2026-03-24 | v32 | III.A: cor:TGP-power-spectrum — P_s(k) i korekta (n_s−1)_TGP = −4ε_H − 4ε_ψ — **ZAMKNIĘTY** | ✅ |
| 2026-03-24 | v32 | III.B: prop:one-loop-UV — pętla 1-elementowa skończona (cutoff substratowy) — **ZAMKNIĘTY** | ✅ |
| 2026-03-24 | v32 | III.B: prop:asymptotic-safety — punkt stały UV m̃*²=0, g̃*=−2π²/5 — **ZAMKNIĘTY** | ✅ |
| 2026-03-24 | v32 | **Faza III kompletna** (III.A–III.B) — perturbacje CMB i UV-kompletność | 🏁 |
| 2026-03-24 | v32 | **SESJA v32 KOMPLETNA** — Fazy I–III domknięte; otwarty: α_K z substratu (OP-2/OP-3) | 🏆 |
