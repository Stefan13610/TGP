---
title: "Phase 0 — pre-registration: op-galactic-substrate-tail — nieekranowany dalekozasięgowy kanał substratowy (sektor fazowy U(1)) jako structural amendment dynamiki galaktycznej; falsyfikatory F-GST-A/B/C/D LOCKED"
type: phase0_balance
status: LOCKED
locked_date: 2026-06-13
cycle: op-galactic-substrate-tail-2026-06-13
authorization: "User 2026-06-13: 'jesteś ekspertem w dziedzinie fizyki teoretycznej; twoje zadanie rozpocząć cykl [HANDOFF_op-galactic-substrate-tail]' → fraza aktywacyjna; per precedens domu #15/#16/#17 (STATE.md sesja #17: 'fraza aktywacyjna pokrywa Phase 0 LOCK') aktywacja obejmuje Phase 0 LOCK. Lektury obowiązkowe HANDOFF §1 (×10) wykonane w zadanej kolejności PRZED niniejszym LOCK. Phase 1+ = osobne 'działaj' każdorazowo."
methodology_binding: "CALIBRATION_PROTOCOL §3.6 + §3.6.6–§3.6.13 BINDING; CYCLE_KICKOFF_TEMPLATE §1–§3 (L1 observable-first); PR-004 werdykt TRIGGERED-FALSIFIED LOCKED NIETYKALNY (ten cykl = structural-amendment path z if_recovery_exhausted, NIE recovery); pipeline porównawczy LOCKED z op-PR004-SPARC-fit-execution (reuse bez modyfikacji)"
anti_lakatos_lock: PRESERVED
PR_candidate: "PR-024 RESERVED (przydział numeru: niniejszy Phase 0; LOCK rejestru: wyłącznie user w Phase FINAL)"
queue_next_after_cycle: "op-nucleation-dimensionality (decyzja user 2026-06-13; osobny Phase 0; zakaz scope-creep)"
---

# Phase 0 — pre-registration LOCKED (przed jakimkolwiek rachunkiem)

## §0 — Tożsamość cyklu (nienegocjowalne)

1. **NIE-recovery:** PR-004 TRIGGERED-FALSIFIED (mechanism, 5.4σ) jest LOCKED i NIETYKALNY.
   Stare χ²_red (TGP Newton-baryon 578/85 vs MOND 50/10.5) pozostają w rejestrze na zawsze.
   Ten cykl realizuje dyspozycję kontraktową PR-004 `if_recovery_exhausted`: „framework needs
   structural amendment" — NOWY mechanizm ⇒ NOWY falsyfikator (kandydat **PR-024**), nigdy
   korekta starego.
2. **HONEST_NEGATIVE = pełnoprawny i wręcz OCZEKIWALNY wynik** (fast-kill §5). Negatywny
   werdykt Phase 1 ⇒ zamknięcie cyklu po Phase 1, przejście do kolejki
   (op-nucleation-dimensionality).
3. Każda faza = osobne „działaj"; po każdej fazie raport z decision menu; STATE.md wpis.

### §0.4 — Pre-flight methodology read confirmation (per KICKOFF §2.6)

- [x] Przeczytano STATE.md sesje #17 (BA) + #18 (PR-004 execution)
- [x] Przeczytano CAŁY cykl op-PR004-SPARC-fit-execution (Phase0 + Phase1_fit.py/.txt + FINAL; dane ./data/ potwierdzone lokalnie)
- [x] Przeczytano PRE_REGISTERED_FALSIFIERS §PR-004 z EXECUTION UPDATE 2026-06-13
- [x] Przeczytano op-CE-H-3D-native-interaction-retry (γ-1 retry CLEAN PASS: −2π log; G_θ = 1/(4πr); konwergencja exp(−m_σL)/L)
- [x] Przeczytano BA Phase4_derivation §1 (m_χ² = λ(w²−Φ₀²); Goldstone bezmasowy w bulku; parzystość w χ)
- [x] Przeczytano FM Phase2_derivation §2 K3 (F_substrat = 0 z ∇⟨Φ⟩ = 0 — MODUŁ; LOCKED)
- [x] Przeczytano CALIBRATION_PROTOCOL §3.6 + §3.6.6–§3.6.13 (BINDING)
- [x] Przeczytano CYCLE_KICKOFF_TEMPLATE §1–§3 (L1 observable-first)
- [x] Przeczytano REALITY_CONTACT_AUDIT_2026-06-12 (kontekst strategiczny)
- [x] Przeczytano SCOPING_op-nucleation-dimensionality (kolejka; NIE aktywować, NIE mieszać zakresów)

**Sign-off:** Claudian (agent wykonawczy cyklu GST) @ 2026-06-13

## §1 — Kontrakt kickoff L1 (observable-first)

```yaml
L1_native:
  output_observable: "v_rot(R) [km/s] — krzywa rotacji; χ²_red bezwymiarowy rezyduał"
  measurement_instrument: "SPARC 175 (Lelli+2016; 3391 pkt; dane lokalne LOCKED z op-PR004-execution)"
  native_coefs_constrained: ["sprzężenie soliton–χ (sektor fazowy U(1))", "klasa zasięgu V_int(r)", "a₀-analog (kombinacja stałych TGP, zero fitu)"]
  falsification_rule: "F-GST-C poniżej (progi liczbowe LOCKED PRZED re-runem)"
  pre_registration_date: "2026-06-13"
L2_framework_reduction:
  target_frameworks: ["Newton-limit (granica zerowego ogona — tożsamość z PR-004 pipeline)"]
  reduction_type: "analytical-exact (wymagana TOŻSAMOŚĆ: v_tail ≡ 0 ⇒ liczby PR-004 odtworzone EXACT)"
  failure_disposition: "L1-stands"
L3_falsification_map:
  - { bound: "SPARC 175 paired test", constrains: "mechanizm GST + a₀-analog", window: "progi F-GST-C", status: "pending" }
  - { bound: "PR-005 GW170817 |Δc/c| ≤ 7×10⁻¹⁶", constrains: "bezmasowy sektor fazowy (dyspersja/piąta siła)", window: "konsystencja deklaratywna w Phase 1", status: "pending" }
```

## §2 — Pytania (zbiory CLOSED; kolejność WYMUSZONA)

**Q1 (mechanizm — fast-kill first, Phase 1):** czy z akcji TGP (S05+Z₂+U(1)+RP²; LIVE
machinery, zero nowych pól) istnieje NIEEKRANOWANY dalekozasięgowy i PRZYCIĄGAJĄCY kanał
oddziaływania międzysolitonowego? Hipotezy (zbiór CLOSED):

- **H-GOLD** — wymiana bezmasowego modu fazowego χ (Goldstone U(1); m_χ²(Φ₀) = 0 EXACT,
  LOCKED BA P4 FP1): solitony z windingiem/teksturą U(1)/RP² sprzęgają się do χ; w 3D
  geometrii punktowej klasa 1/r, dla sprzężeń topologicznych liniowych — log (γ-1 LOCKED
  jako wsparcie klasy, z jawną deklaracją 2D-proxy ≠ 3D).
- **H-SCREEN** — brak użytecznego kanału ⇒ HONEST_NEGATIVE i zamknięcie cyklu po Phase 1.
  Pod-przypadki zadeklarowane TERAZ (przed rachunkiem; wszystkie ⇒ ta sama klasa werdyktu):
  (a) ekranowanie — jedyny kanał niesie efektywną masę ⇒ zanik exp na zasięgu 1/m_eff;
  (b) decoupling — sprzężenie statycznego solitonu do χ tożsamościowo zerowe
  (uwaga pre-derywacyjna §4.3: ładunek Noether statycznej konfiguracji = 0 — to jest
  NAJBARDZIEJ prawdopodobny nóż fast-killu);
  (c) zły znak — kanał istnieje, ale między jednoimiennymi solitonami materii jest
  ODPYCHAJĄCY (γ-1: same-sign repel) ⇒ nie może generować dodatkowego wiązania galaktycznego.

**Q2 (skala — TYLKO jeśli Q1 = H-GOLD_DERIVED):** wyprowadź a₀-analog z stałych TGP bez
jakiegokolwiek fitu. Kandydat strukturalny (JEDEN, zapisany tu przed rachunkiem,
INFORMATIONAL): klasa **a₀-analog ∈ O(1)·c/t = O(1)·cH** (γ-3: R = ct LOCKED; marginalność
FCR wiąże dynamikę lokalną z ct przez GM/(ct)). Współczynnik O(1) musi WYNIKAĆ z mechanizmu
Q1 — nie z menu {1/2π, 1/6, ...} dobieranego do a₀_obs. **a₀_obs = 1.2×10⁻¹⁰ m/s²
WYŁĄCZNIE comparison-only po LOCKu wartości pochodnej; NIGDY input** (circularity guard FP
w każdym sympy).

**Q3 (test — NOWY PR, Phase 3):** re-run SPARC 175 z mechanizmem Q1+Q2:
**identyczny zero-parametrowy pipeline LOCKED z op-PR004-execution** (te same Υ_d = 0.5,
Υ_b = 0.7; te same filtry e_V ≤ 0 / V_obs ≤ 0 / R ≤ 0; ta sama definicja χ²_g i agregatów
GLOBAL/MEDIAN; ten sam paired test + bootstrap seed 42, 10⁴; ta sama podpróba Q1+Q2;
te same pliki danych lokalne). Jedyna zmiana = model: v_model²(R) = v_bar²(R) + v_tail²(R)
z v_tail w pełni wyznaczonym przez Q1+Q2 (zero parametrów wolnych); **diff modelu jawny
w skrypcie** (reuse Phase1_fit.py). **Przegrana = koniec ścieżki bez recovery** (zapis
kontraktowy — patrz F-GST-C/FAIL).

## §3 — Falsyfikatory (zbiory CLOSED; progi LOCKED 2026-06-13)

### F-GST-A — mechanizm (Q1; werdykt Phase 1)

Klasy (CLOSED): **H-GOLD_DERIVED / H-SCREEN_NEGATIVE / GAP / INDETERMINATE.**

**H-GOLD_DERIVED wymaga WSZYSTKICH (mechanicznie):**
1. m_χ²(Φ₀) = 0 EXACT odtworzone symbolicznie (reuse LOCKED BA P4 FP1) ORAZ audyt
   U(1)-niezmienniczości akcji LIVE: zero monomianów łamiących U(1) (free-symbols audit;
   kompaktowość RP² rozliczona jawnie — §4.2);
2. sprzężenie soliton–χ wyprowadzone z akcji, symbolicznie NIEZEROWE dla konfiguracji
   statycznej lub topologicznej (nie postulowane);
3. klasa zasięgu V_int(r) w 3D = power-law (1/r-klasa) lub log — jeśli weryfikowana
   numerycznie: R²_form ≥ 0.95 przy porównaniu equal-param 2 vs 2 (§3.6.7) z przewagą
   ≥ 0.02 nad czystym exp; współczynnik wiodący zgodny z pre-derywacją §4 w ±5% (§3.6.9);
4. znak: kanał PRZYCIĄGAJĄCY między jednoimiennymi solitonami materii (konwencja §4.4;
   znak pre-derived, nie dobrany po wyniku);
5. zgodność z K3 LOCKED: wkład sektora fazowego do siły TŁA w jednorodnym nasyconym bulku
   = 0 EXACT (∇⟨Φ⟩ = 0 dotyczy MODUŁU; sektor fazowy może wnosić WYŁĄCZNIE człony
   międzysolitonowe znikające przy usunięciu defektów). Sprzeczność z K3 ⇒
   H-SCREEN_NEGATIVE (NIE reinterpretacja K3).

**H-SCREEN_NEGATIVE:** dowolny z pod-przypadków (a)/(b)/(c) §2-Q1 wykazany symbolicznie
⇒ HONEST_NEGATIVE; cykl ZAMYKA się po Phase 1 (Q2/Q3 nie są wykonywane).
**GAP:** wyprowadzenie wymaga elementu spoza LIVE (np. sektorów holonomii RP² / textured
configurations nie wyprowadzonych) — luka deklarowana, zero nowych pól.
**INDETERMINATE:** rachunek nierozstrzygalny w zadeklarowanych przybliżeniach (wyliczyć które).

### F-GST-B — skala (Q2; werdykt Phase 2; wykonywana TYLKO przy H-GOLD_DERIVED)

Klasy (CLOSED): **DERIVED / DERIVED_WITH_ANCHOR / GAP.**

- **DERIVED:** a₀-analog = zamknięta forma z {λ, Φ₀, c, t} i wielkości LOCKED (FCR/FM/γ-3)
  bez żadnego anchora obserwacyjnego i zero nowych stałych;
- **DERIVED_WITH_ANCHOR:** forma wymaga anchora klasy (γ) §3.6.13 (np. t_universe / H = 1/t
  z wieków gwiazd — jak γ-3); anchor zadeklarowany jawnie; NADAL zero fitu;
- **GAP:** kombinacja stałych TGP nie wycina skali przyspieszeniowej (brak = GAP, nie nowa stała).

**Protokół value-blind:** wartość liczbowa a₀-analog jest LOCKOWANA w Phase2_derivation.md
PRZED jakimkolwiek porównaniem z a₀_obs; porównanie (comparison-only, INFORMATIONAL) —
wyłącznie PO locku, w osobnej sekcji. Żaden wybór gałęzi/współczynnika nie może być
motywowany odległością od 1.2×10⁻¹⁰ (forbidden move #13).

### F-GST-C — SPARC re-run (Q3; werdykt Phase 3; progi liczbowe LOCKED)

Statystyka decyzyjna (identyczna konstrukcja jak PR-004-execution):
d_g = χ²_g(TGP-tail)/N_g − χ²_g(MOND simple)/N_g per galaktyka; t = mean(d)/SEM(d);
bootstrap 10⁴ (seed 42) sanity. Decyzyjny agregat = paired t na pełnej próbie;
Q1+Q2 = secondary (pre-deklarowany allowed direction, jak PR-004).

| Wynik | Warunek LOCKED | Dyspozycja |
|---|---|---|
| **PASS** | mean(d) ≤ 0 (TGP-tail ≤ MOND w teście sparowanym) | PR-024 LOCK eligible (decyzja wyłącznie user, Phase FINAL) |
| **PARTIAL** | 0 < t < 5 | pasmo pośrednie; klasyfikacja agregatem F-GST-D; żadnego łagodzenia ex post |
| **FAIL** | t ≥ 5 | **koniec ścieżki substrate-tail BEZ recovery** (zapis kontraktowy); werdykt propagowany do rejestru |

**Warunki ważności re-runu (FP-klasa, wszystkie wymagane):**
- **Tożsamość pipeline'u:** przy v_tail ≡ 0 skrypt MUSI odtworzyć liczby PR-004 EXACT
  (578.14 / 49.99 GLOBAL; 85.23 / 10.51 median; t = 5.4σ) — guard przeciw cichej zmianie
  pipeline'u;
- **Zero parametrów:** v_tail bez żadnej wielkości dopasowywanej (FP6-klasa: zero wywołań
  optymalizatora; Υ, a₀-analog FIXED z Phase 2 PRZED odczytem jakiegokolwiek χ²);
- **Circularity guard:** a₀_obs nieobecne w modelu TGP-tail (a₀ literaturowe występuje
  WYŁĄCZNIE w benchmarku MOND, jak w PR-004);
- **Raport pomocniczy (INFORMATIONAL, nie-werdyktowy):** paired t vs Newton-baryon baseline
  (czy ogon w ogóle poprawia względem PR-004) + rozkład wygranych galaktyk.

### F-GST-D — agregat (Phase FINAL)

Mapowanie (CLOSED):

| F-GST-A | F-GST-B | F-GST-C | F-GST-D |
|---|---|---|---|
| H-SCREEN_NEGATIVE | — (nie wykonywana) | — | **HONEST_NEGATIVE** (zamknięcie po Phase 1) |
| GAP / INDETERMINATE | — | — | **GAP_CLOSURE** (zamknięcie z luką deklarowaną) |
| H-GOLD_DERIVED | GAP | — | **MECHANISM_WITHOUT_SCALE** (zamknięcie; Q3 niewykonalne bez skali — zakaz fitu) |
| H-GOLD_DERIVED | DERIVED(±ANCHOR) | PASS | **TAIL_VIABLE — PR-024 LOCK eligible** (decyzja user) |
| H-GOLD_DERIVED | DERIVED(±ANCHOR) | PARTIAL | **TAIL_PARTIAL** (klasyfikacja + dyspozycja = user) |
| H-GOLD_DERIVED | DERIVED(±ANCHOR) | FAIL | **TAIL_FALSIFIED — koniec ścieżki bez recovery** |

Decyzja o LOCKu PR-024 w rejestrze = WYŁĄCZNIE user (Phase FINAL). Numer PR-024
zarezerwowany niniejszym; żaden append przed FINAL (forbidden move #14).

## §4 — Analytical pre-derivation (§3.6.1–§3.6.9; PRZED rachunkiem)

### §4.1 — Spektrum i propagatory (oczekiwane wartości; reuse LOCKED)

- Fluktuacje wokół nasyconego tła: Φ = (Φ₀ + h) + iχ (konwencja kanoniczna BA P4 §1.1).
  Oczekiwane (LOCKED BA P4 FP1): m_h² = 2λΦ₀² (≡ m_σ² = 2λv², CE-H), **m_χ²(Φ₀) = 0 EXACT**.
- Propagatory 3D (LOCKED γ-1 retry): **G_χ(r) = 1/(4πr)** (bezmasowy), G_h(r) = e^(−m_σ r)/(4πr)
  (ekranowany). Oczekiwanie klasowe: kanał modułowy ekranowany na 1/m_σ (mikroskopijny);
  jedyny kandydat dalekozasięgowy = sektor fazowy.

### §4.2 — Czy RP²/kompaktowość generuje masę χ? (do weryfikacji w Phase 1)

W LIVE akcji V = V(|Φ|²) — U(1) dokładne, zero monomianów łamiących fazę (BA P4 FP2-FP3:
akcja parzysta w χ; tu analogicznie: niezależna od przesunięcia fazy globalnej) ⇒
**perturbacyjnie masa χ = 0 w bulku, oczekiwane EXACT.** Kompaktowość (RP² = S²/Z₂;
identyfikacja antypodalna) nie generuje masy perturbacyjnej; może najwyżej zmieniać
spektrum defektów (sektory windingu Z zamiast pełnego U(1) — wpływa na ŁADUNKI, nie na
masę mediatora). Niespodzianka w tym punkcie (wyprowadzona masa > 0) ⇒ H-SCREEN(a).

### §4.3 — Sprzężenie soliton–χ (centralny nóż fast-killu; trzy kanały, zbiór CLOSED)

(i) **Kanał ładunku Noether:** Q = ∫ j⁰ = ∫ Φ₀² θ̇ d³x. Dla konfiguracji STATYCZNEJ θ̇ = 0
⇒ **oczekiwane Q = 0 EXACT ⇒ kanał Coulomb-przez-ładunek wymaga rotacji fazy w czasie
(klasa Q-ball) — w LIVE inwentarzu solitonów nie zadeklarowana** ⇒ kanał (i) prawdopodobnie
zamknięty. To jest uczciwie najkrótsza droga do H-SCREEN(b).
(ii) **Kanał topologiczny (winding/tekstura):** soliton z nietrywialnym uzwojeniem fazy
(π₁(U(1)) dla konfiguracji liniowych; π₂(RP²) dla punktowych) wymusza ∇θ ≠ 0 wokół defektu;
energia E ⊃ ½Φ₀²(∇θ)² konfiguracji dwudefektowej daje V_int(r). Oczekiwane klasy: geometria
liniowa ⇒ −2πΦ₀²n₁n₂ log(L/r₀) na jednostkę długości (LOCKED γ-1; **2D-proxy — deklaracja:
NIE jest to 3D claim**); geometria punktowa 3D ⇒ klasa 1/r z ładunkiem topologicznym.
Pre-derywacja zupełnia: czy bariony-solitony TGP (FFS: obiekty kwarkowe; defekty punktowe
π₂) NIOSĄ niezerowy ładunek tego kanału — to jest rachunek Phase 1, wynik OTWARTY.
(iii) **Kanał indukowany moduł–faza:** mieszanie h–χ wokół tła defektu — oczekiwane
ekranowane (każdy wkład przez h niesie e^(−m_σ r)) ⇒ nie-dalekozasięgowy.

**Zbiór kanałów CLOSED:** {(i), (ii), (iii)}. Żaden inny kanał nie może być dodany mid-cycle
(luka = GAP).

### §4.4 — Konwencje znaków (§3.6.6; LOCKED przed rachunkiem)

- **Konwencja:** V_int(r) = E(konfiguracja w separacji r) − E(∞); siła F(r) = −dV/dr;
  **przyciąganie ⟺ V rośnie z r ⟺ F < 0** (skierowana do zmniejszania r).
- **Zasada fizyczna (limit znany):** γ-1 LOCKED — jednoimienne ładunki 2D-Coulomb
  ODPYCHAJĄ (V = −2πv²n₁n₂log(L/r₀): dla n₁n₂ > 0 V maleje z L ⇒ F > 0 ⇒ repulsja ✓).
  Analogicznie dla wymiany Goldstone'a między statycznymi ładunkami topologicznymi:
  oczekiwanie generyczne = **jednoimienne ODPYCHAJĄ, różnoimienne przyciągają**
  (elektrostatyka-analog). Materia barionowa = solitony JEDNOIMIENNE (wspólny znak ładunku
  topologicznego per FFS) ⇒ **uczciwe oczekiwanie pre-derywacyjne: znak NIEKORZYSTNY
  (repulsja) — pod-przypadek H-SCREEN(c)**. Wyjątki możliwe (sprzężenia pochodne wyższego
  rzędu, kanały tensorowe tekstur) — rozstrzygnie rachunek; znak NIE może być przedefiniowany
  po wyniku.

### §4.5 — Geometria 2D-proxy vs 3D (deklaracja wiążąca)

Wynik γ-1 (log, −2π) uzyskano w geometrii wirowej z niezmienniczością z-translacyjną
(efektywnie 2D). Dla galaktyk relewantna jest geometria PUNKTOWYCH źródeł w 3D
(G_χ = 1/(4πr) ⇒ klasa 1/r, NIE log). Każde twierdzenie Phase 1+ musi jawnie deklarować
geometrię; log-form wolno cytować wyłącznie jako wsparcie KLASY (istnienie nieekranowanego
kanału w machinery), nie jako 3D-wynik (forbidden move #16).

### §4.6 — Enumeracja założeń (§3.6.8)

(a) tło: nasycony bulk ⟨|Φ|⟩ = Φ₀ = const, linearizacja wokół VEV; defekty = źródła
zlokalizowane; (b) normalizacja: L = ½|∂Φ|², jednostki naturalne w wyprowadzeniu (c = 1
przywracane wymiarowo na końcu); (c) granice: statyczność źródeł, separacja r ≫ rozmiar
rdzenia r₀, pole słabe (superpozycja liniowa w dalekim polu — F-γ-2 LOCKED wspiera);
(d) podstawienia efektywne: ładunki topologiczne n_i symboliczne (żadnego q = 1 implicit
bez deklaracji); (e) symetrie narzucone: sferyczność/osiowość konfiguracji testowych
deklarowana per-rachunek. W Phase 3: założenia pipeline'u = verbatim PR-004 Phase 0 §2.

### §4.7 — Precyzja (§3.6.9)

Standard ±5% dla każdej wielkości fitowanej numerycznie względem wartości analitycznej
(precedens γ-1: |B| vs 2π — 0.4%). Tożsamość pipeline (F-GST-C) = EXACT (odtworzenie liczb
PR-004 co do reprezentacji float).

### §4.8 — Klasyfikacja stałych (§3.6.13)

| Stała | Klasa | Nota |
|---|---|---|
| λ, Φ₀ | (α) TGP_FUNDAMENTAL | symboliczne; zero wartości liczbowych w wyprowadzeniu |
| m_σ = √(2λ)Φ₀ | (α) pochodna | LOCKED CE-H |
| c | (β) EMERGENT_FROM_PHI | reżim galaktyczny: c = c₀ uzasadnić deklaratywnie (saturated bulk, jak FM) |
| t (wiek) / H = 1/t | (γ) OBSERVATIONAL_ANCHOR | dozwolone TYLKO w F-GST-B klasie DERIVED_WITH_ANCHOR; jawnie |
| G | poza Lagrangianem TGP (AX-DECL declared limit) | wolno używać wyłącznie w księgowości porównawczej (jak FCR/FM), nie jako nowy element mechanizmu |
| a₀_obs = 1.2×10⁻¹⁰ | comparison-only | NIGDY input; benchmark MOND w pipeline = jedyne legalne wystąpienie |
| Υ_d = 0.5, Υ_b = 0.7 | LITERATURE_ANCHORED (LOCKED PR-004) | FIXED; nietykalne |
| **Nowe stałe fundamentalne** | **budżet: 0** | luka = GAP |

## §5 — Plan faz (fast-kill first)

| Faza | Zakres | Werdykt | Gate |
|---|---|---|---|
| **Phase 1 — FAST-KILL (Q1)** | najtańszy decydujący rachunek (precedens BA P1): U(1)-audyt + masa χ (FP1-FP2) → ładunek Noether statyczny (FP3) → sprzężenie topologiczne i klasa V_int(r) w 3D (FP4) → ZNAK (FP5) → zgodność K3 (FP6) → audyt ekranowania (FP7) → circularity guard (FP8) | F-GST-A | „działaj" |
| Phase 2 (TYLKO H-GOLD_DERIVED) | a₀-analog z stałych TGP; LOCK wartości PRZED porównaniem; porównanie comparison-only osobno | F-GST-B | „działaj" |
| Phase 3 (TYLKO B ∈ DERIVED*) | SPARC re-run: reuse Phase1_fit.py z podmienionym modelem (diff jawny); tożsamość pipeline przy v_tail ≡ 0; progi F-GST-C mechanicznie | F-GST-C | „działaj" |
| Phase FINAL | agregat F-GST-D; DOUBTS register; decyzja PR-024 = WYŁĄCZNIE user | F-GST-D | „działaj" |

Każda faza: PhaseN_derivation.md + PhaseN_sympy.py/.txt (PASS/FAIL per FP; 0 hardcoded
T_pass; circularity-guard FP: free-symbols audit na a₀_obs/V_obs/G_obs w każdym skrypcie).
HONEST_NEGATIVE w Phase 1 ⇒ od razu Phase FINAL (zamknięcie), bez Q2/Q3.

## §6 — Forbidden moves (LOCKED; 16)

1. Modyfikacja PR-004 TRIGGERED-FALSIFIED (oraz PR-022/FM/FCR/BA/γ-*/CE-H — wszystkie
   LOCKED read-only); nowy mechanizm = NOWY PR.
2. ρ_DM lub jakiekolwiek pole spoza S05+Z₂+U(1)+RP² (luka = GAP, nie nowe pole).
3. Użycie a₀_obs, V_obs, jakiejkolwiek danej SPARC w WYPROWADZENIU (dane wyłącznie
   w fazie porównawczej Q3 — analog zakazu G_obs; guard FP w każdym sympy).
4. Fitowanie czegokolwiek (Υ, a₀-analog, amplituda ogona, r₀...) — pipeline zero-parametrowy;
   zakaz wywołań optymalizatorów (FP6-klasa).
5. Zmiana pipeline'u porównawczego LOCKED z PR-004-execution (identyczność = wiarygodność
   porównania; guard tożsamości v_tail ≡ 0).
6. Reinterpretacja K3 (F_substrat = 0 z ∇⟨Φ⟩ = 0 — MODUŁ): sprzeczność sektora fazowego
   z tą granicą ⇒ HONEST_NEGATIVE, nie korekta K3.
7. Nowe stałe fundamentalne (budżet 0; λ, Φ₀ symboliczne).
8. Hardcoded T_pass.
9. Miękkie domknięcie / przeciąganie cyklu po negatywnym fast-killu (HONEST_NEGATIVE
   po Phase 1 = sukces metodologiczny — zamknąć i przejść do kolejki).
10. Scope-creep op-nucleation-dimensionality (następny cykl, osobna pre-rejestracja).
11. Zmiana klas/progów falsyfikatorów ex post (zbiory CLOSED; reguła IMMUTABLE).
12. Wybór konwencji znaku po obejrzeniu wyniku (§4.4 LOCKED).
13. Numerologia a₀: dobieranie kombinacji stałych/współczynników do 1.2×10⁻¹⁰ (value-blind;
    jeden kandydat strukturalny zapisany §2-Q2; współczynnik O(1) musi wynikać z mechanizmu).
14. Append/LOCK PR-024 przed Phase FINAL lub bez decyzji user.
15. Cytowanie H-SORT jako ustalonej bariogenezy (kalibracja epistemiczna user, SCOPING §2).
16. Cytowanie γ-1 log-form jako wyniku 3D bez jawnego wyprowadzenia geometrii punktowej
    (2D-proxy ≠ 3D claim).

## §7 — Risk register

| ID | Ryzyko | Severity | Mitygacja |
|---|---|---|---|
| R-GST-1 | Ekranowanie: jedyny kanał z masą efektywną ⇒ zanik exp | HIGH | fast-kill Phase 1 (FP7); HONEST_NEGATIVE tanio |
| R-GST-2 | Decoupling: statyczny soliton ma Q_Noether = 0 ⇒ brak sprzężenia do χ | HIGH | FP3 wprost; kanał (ii) topologiczny jako jedyna alternatywa — rachunek, nie życzenie |
| R-GST-3 | Zły znak: jednoimienne solitony się ODPYCHAJĄ (γ-1 precedens) ⇒ ogon pogarsza krzywe | HIGH | FP5; znak pre-derived §4.4; repulsja ⇒ H-SCREEN(c), koniec |
| R-GST-4 | Sprzeczność z K3 (F_substrat = 0): sektor fazowy wnosi siłę tła w bulku | HIGH | FP6; granica: wkład tła = 0 EXACT przy braku defektów; sprzeczność ⇒ NEGATIVE |
| R-GST-5 | 2D→3D: γ-1 log był w geometrii wirowej; w 3D punktowo 1/r — flat curves NIE z konstrukcji | MED-HIGH | deklaracja §4.5; geometria jawna per twierdzenie |
| R-GST-6 | Goldstone vs obserwacje: bezmasowy skalar ⇒ piąte siły / dyspersja GW (PR-005 Δc/c ≤ 7×10⁻¹⁶) | MED | konsystencja deklaratywna w Phase 1 (L3 map); konflikt raportowany, nie ukrywany |
| R-GST-7 | Numerologia a₀ (pokusa dobrania do 1.2×10⁻¹⁰) | MED-HIGH | value-blind protokół F-GST-B; LOCK przed porównaniem; forbidden #13 |
| R-GST-8 | Skalowanie BTFR: uniwersalny log/1/r-ogon daje v² ∝ N (masa) zamiast v⁴ ∝ M — napięcie z Tully-Fisherem wbudowane w klasę | MED-HIGH | zapisane PRZED rachunkiem; żadnego dorabiania skalowania ex post; rozstrzyga mechanicznie F-GST-C |
| R-GST-9 | Pokusa "rozwiązania ciemnej materii" — najsilniejsza w historii projektu | META/HIGH | progi LOCKED zanim policzony pierwszy χ²; anticipated outcomes §8; decyzje PR wyłącznie user |

## §8 — Anticipated outcomes (INFORMATIONAL; kierunki pokusy zapisane PRZED rachunkiem)

1. **Najbardziej prawdopodobny wynik: HONEST_NEGATIVE w Phase 1** — przez decoupling (R-GST-2)
   lub zły znak (R-GST-3). Historia projektu (rotacje galaktyk: wielokrotne porażki,
   PR-004 5.4σ) i generyka teorii pola (Goldstone globalny sprzęga się pochodnie; statyczne
   źródła bez ładunku) wskazują ten kierunek. Zapisujemy, by negatyw nie był łagodzony,
   a pozytyw wymagał wzmożonego audytu.
2. **Kierunek pokusy #1:** „mechanizm prawie działa" — naciąganie kanału topologicznego (ii)
   do przyciągania przez selektywny dobór konfiguracji. Antidotum: zbiór kanałów CLOSED,
   znak pre-derived, konfiguracje deklarowane przed rachunkiem.
3. **Kierunek pokusy #2:** numerologia a₀ (cH/2π vs cH/6 — oba „blisko"); antidotum:
   współczynnik musi wynikać z mechanizmu, porównanie po LOCKu.
4. **Kierunek pokusy #3:** przy PARTIAL w Q3 — reinterpretacja progu lub doszukiwanie
   podpróbek. Antidotum: podpróby zamknięte (pełna + Q1+Q2), progi IMMUTABLE.
5. **Sukces (TAIL_VIABLE) = „rozwiązanie ciemnej materii"** — najsilniejsza pokusa w historii
   projektu; właśnie dlatego wszystkie progi powyżej są LOCKED zanim policzymy pierwsze χ²,
   a decyzja PR-024 należy wyłącznie do user.

## §9 — Outcome sets per FP (§3.6.11)

Każdy FP Phase 1-3 deklaruje zbiór: {PASS / PARTIAL_compute / FAIL} (+ PARTIAL_concept_mismatch
tam, gdzie TGP-native nie ma odpowiednika — wymaga jawnego uzasadnienia per §3.6.11(b));
budżet PARTIAL_compute: max 1 na cykl. Werdykty F-GST-* wyłącznie w klasach CLOSED §3.

## §10 — Anti-Lakatos compliance (Phase 0)

Zbiory hipotez/klas/kanałów CLOSED przed rachunkiem ✓ · progi liczbowe LOCKED (F-GST-C) ✓ ·
analytical pre-derivation z konwencją znaku i enumeracją założeń ✓ · anticipated outcomes
z kierunkami pokusy ✓ · 0 nowych stałych ✓ · 0 modyfikacji werdyktów poprzedników (PR-004
NIETYKALNY) ✓ · a₀_obs/V_obs/G_obs comparison-only z guardem FP ✓ · PR-024 RESERVED,
append wyłącznie user w FINAL ✓ · kolejka po cyklu: op-nucleation-dimensionality (bez
scope-creep) ✓.

**Phase 0 LOCKED 2026-06-13. Następny krok: Phase 1 FAST-KILL — wymaga user „działaj".**
