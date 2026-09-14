---
title: "Phase0_balance (LOCK) — op-collapse-matter-source: status kolapsu do granicy dziedziny — czy korpusowe sprzężenie z materią (eq:L-mat-unified) deformuje próżnię Yukawa-podobnie i czy stabilizuje kolapsujące starty?"
date: 2026-09-14
type: phase0-balance
tgp_owner: research/op-collapse-matter-source-2026-09-14
status: LOCKED
anti_lakatos_lock: ACTIVE
related:
  - "[[HANDOFF_PROMPT.md]]"
  - "[[README.md]]"
  - "[[../op-r3-stationary-states-2026-09-14/NEEDS.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Phase 0 — LOCK (zero obliczeń przed zapisem tego pliku)

**Autoryzacja:** user-gate 2026-09-14 — wybór „N2: status kolapsu / sprzężenie
z materią" z [[../op-r3-stationary-states-2026-09-14/NEEDS.md]] (N2);
realizuje jednocześnie kandydata LOCKa N2 z
[[../op-metric-pair-M911-2026-09-02/NEEDS.md]] (sprzężenie z materią,
osobny LOCK) — w klasie dynamiki 2. rzędu (N3 tamtego cyklu).

## 0. Kontekst i pytania

Poprzednik (`op-r3-stationary-states`, Q-E-INCONCLUSIVE) ustalił nową twardą
obserwację: w czystej dynamice 2. rzędu gałęzi zdrowej **kolaps do granicy
dziedziny jest regułą** (6/10 startów, t≈2.3–19.8, zbieżnie), podczas gdy
w gradient flow był wyjątkiem. NEEDS N2 stawia pytanie: czy BREAKDOWN-BOUNDARY
to (a) artefakt braku członów stabilizujących (materia), (b) fizyka substratu,
(c) sygnał przeciw gałęzi zdrowej. Ten cykl testuje odczyt (a) JEDYNĄ formą
sprzężenia obecną w rdzeniu — bez wymyślania nowych członów (S05 nienaruszone;
materia = źródło zewnętrzne ρ, nie drugie pole dynamiczne).

- **Q-H1 (binarne pytanie M911-N2):** czy statyczne źródło ρ na tle próżni
  (i) tylko DEFORMUJE próżnię (Yukawa-podobnie, zgodnie z odpowiedzią
  zlinearyzowaną), czy (ii) INDUKUJE ucieczkę/kreację (wciąga pole pod próg
  dolny 5/6 / nad górny 7/6 lub do pasa granicznego)?
- **Q-H2 (centralne dla NEEDS-N2):** czy sprzężenie z materią STABILIZUJE
  starty, które w cyklu poprzednika kolabowały (w tym oba quasi-R3)?

**Zakres:** klasyfikacja odpowiedzi na źródło + stabilizacja kolapsu.
ZAKAZ claimów o masach leptonów i o oscylonach (to cykl-bliźniak
`op-oscillon-small-amplitude-2026-09-14`).

## 1. Model (formy FROZEN)

Sektor pola — CYTAT jak w poprzedniku (Q-D1-PASS, konwencja |g^tt|),
K_geo=γ=c₀=1: M=ψ⁶/(4−3ψ)², 𝒦=ψ⁴, 𝒰=ψ⁴/4−ψ³/3; π=Mψ̇; dziedzina (0,4/3).

**Człon materii — z RDZENIA, do WYPROWADZENIA w Phase 1 (nie założenia):**
`eq:L-mat-unified` (sek08a): L_mat=−(q/Φ₀)ψρ, ρ≡−T^μ_μ/c₀²≥0 (L01 formal
definition; decyzja post-audit 2026-05-01: czynnik ψ = konsekwencja elementu
objętości, nie dilaton B-D). Z √−g=c₀ψ/(4−3ψ) (eq:vol-element-M911):

S_mat = ∫√−g·L_mat ⟹ **𝒰_mat(ψ,r) = λ̃·ρ̂(r)·ψ²/(4−3ψ)** — forma oczekiwana
(P1-H1 wyprowadza ją sympy z zapisu literalnego; gate tożsamości),
gdzie λ̃ ≡ (q c₀/Φ₀)·ρ₀ (bezwymiarowa siła sprzężenia, λ̃>0 bo ρ≥0),
ρ̂(r)=exp(−r²/18) (gauss σ_ρ=3, FROZEN — skala rdzeni poprzednika).

EOM z materią: Mψ̈+½M′ψ̇² = (1/r²)(r²𝒦ψ′)′ − ½𝒦′ψ′² − 𝒰′ − ∂𝒰_mat/∂ψ.

## 2. Phase 1 — analityka pre-rejestrowana (sympy, simplify=0; PRZED numeryką)

- **P1-H1 (wyprowadzenie):** 𝒰_mat z literalnego √−g·(q/Φ₀)ψρ — gate: tożsame
  z λ̃ρ̂ψ²/(4−3ψ); ∂𝒰_mat/∂ψ = λ̃ρ̂·ψ(8−3ψ)/(4−3ψ)² — gate sympy 1e−12.
- **P1-H2 (predykcja pre-rejestrowana — odpowiedź zlinearyzowana):** wokół
  ψ=1: (−∇²+1)δψ = −5λ̃ρ̂ ⟹ δψ = −5λ̃·(G_Yuk∗ρ̂), G_Yuk=e^{−r}/(4πr).
  **Znak: δψ<0 (materia obniża ψ — pogłębia dylatację czasu)**; ogon Yukawy
  e^{−r}/r. Wyprowadzenie + ewaluacja numeryczna splotu (kwadratura, nie
  ewolucja) jako wzorzec dla gate'u P3-H1a.
- **P1-H3 (fakty brzegowe, pre-rejestrowane):** 𝒰_mat→+∞ przy ψ→4/3⁻
  (źródło ODPYCHA od sufitu tam, gdzie ρ̂>0); 𝒰_mat→0 przy ψ→0⁺ (podłoga
  bez bezpośredniej bariery od materii). Przewidywana ASYMETRIA: stabilizacja
  łatwiejsza dla kolapsów górnych niż dolnych — do konfrontacji w Q-H2.

## 3. Protokół numeryczny (FROZEN)

Silnik: adaptacja `../op-r3-stationary-states-2026-09-14/engine_core.py`
(kopia z cytatem; jedyna zmiana merytoryczna: człon −∂𝒰_mat/∂ψ w RHS +
𝒰_mat w ewaluatorze energii — z tożsamością bez kancelacji dla 𝒰(ψ)−𝒰(1)
jak u poprzednika). Siatka r_i=(i+½)h, h=0.05 (główna) / h=0.025
(potwierdzenia), R=200, sponge smootherstep γ₀=1 na [160,200], dt=0.005,
E_core r≤80, E_ref=E_core(50); zapis ψ(0,t), E_core(t) co dt_out=0.1.

**Phase 2 — bramka (FROZEN; FAIL ⟹ STOP):**
- P2a: próżnia BEZ źródła (λ̃=0), 100 T₀: ‖ψ−1‖∞ ≤ 1e−10 (obie siatki).
- P2b: regresja λ̃=0 — start qR3 a=−0.20 (jak u poprzednika): COLLAPSE
  z czasem zdarzenia 4.74 ±2% (h=0.05).
- P2c: zachowanie energii ZE źródłem statycznym (λ̃=0.05, start a=+0.05 σ=3,
  do t=700): dryf ≤1e−6/100 T₀ (energia z członem 𝒰_mat jest zachowana,
  bo ρ̂ nie zależy od t).

**Phase 3 — Q-H1 (źródło na próżni):** start ψ≡1, π₀=0, źródło włączone
od t=0; **λ̃∈{0.01, 0.05, 0.2, 0.5}**; t_max=300; h=0.05 (λ̃=0.01 i 0.5
dodatkowo h=0.025). Klasyfikacja per λ̃:
- **DEFORMATION:** pole osiada do statycznego profilu (‖ψ̇‖∞→ poziom szumu;
  średnia po oknie [250,300]); dla λ̃=0.01 gate ilościowy: profil vs
  −5λ̃(G_Yuk∗ρ̂) zgodny ≤5% względem max|δψ| (na r≤40).
- **THRESHOLD-PULL:** osiadły profil przekracza próg ψ<5/6 (detektory M911)
  w sposób trwały (okno [250,300]).
- **COLLAPSE:** pas graniczny (ψ<1e−6 / ψ>4/3−1e−6) lub niefinityczność
  (nadkategoria jak w cyklu-bliźniaku).
Deskryptywnie obowiązkowo: głębokość δψ(0) vs predykcja liniowa dla
wszystkich λ̃ (nieliniowe odchylenie = wynik, nie błąd).

**Phase 3 — Q-H2 (stabilizacja kolapsu):** 4 starty-reprezentanci kolapsu
poprzednika (FROZEN): qR3 a=−0.20, qR3 a=+0.20, gauss a=−0.30 σ=3,
gauss a=+0.15 σ=6 — identyczne profile początkowe jak u poprzednika,
π₀=0, źródło włączone od t=0; **λ̃∈{0.05, 0.2, 0.5}**; t_max=1000; h=0.05.
Kategorie: **COLLAPSE** (nadkategoria, podtyp raportowany) / **RADIATED**
(E_core<0.05·E_ref, bez zdarzenia brzegowego) / **STABILIZED** (przeżywa
t_max bez zdarzenia brzegowego i bez zaniku: E_core(t_max)≥0.05·E_ref;
podtypy deskryptywne: osiadły statycznie vs oscylujący) / INCONCLUSIVE-RUN.
**Potwierdzenia (FROZEN):** każdy bieg, którego kategoria RÓŻNI SIĘ od
kategorii λ̃=0 (baseline poprzednika = COLLAPSE) → h=0.025 ORAZ dt/2
(kategoria zgodna; czas zdarzenia/τ ±10%). Dodatkowo obowiązkowo jeden
bieg COLLAPSE (najniższe λ̃, qR3 a=+0.20) → h=0.025 (kontrola negatywu).

## 4. Werdykty (litera; INCONCLUSIVE ≠ pozytyw)

- **Q-H1-DEFORMATION:** wszystkie λ̃ z zamrożonej listy dają DEFORMATION
  zbieżnie, a gate liniowy (λ̃=0.01) przechodzi ≤5%. (= odpowiedź „Yukawa-
  podobna deformacja" na binarne pytanie M911-N2 w zbadanym zakresie.)
- **Q-H1-PULL:** ≥1 λ̃ daje THRESHOLD-PULL lub COLLAPSE zbieżnie
  (= źródło potrafi wciągnąć pole za progi; λ̃_crit raportowane opisowo).
- **Q-H1-INCONCLUSIVE:** inaczej.
- **Q-H2-PASS (stabilizacja istnieje):** ≥1 para (start, λ̃) z kategorią
  STABILIZED potwierdzoną (h/2 i dt/2), której baseline λ̃=0 = COLLAPSE.
- **Q-H2-FAIL:** wszystkie 12 par pozostają COLLAPSE zbieżnie (kategoria
  na h; potwierdzenie wg reguły §3).
- **Q-H2-INCONCLUSIVE:** inaczej.

## 5. Drzewo decyzyjne (pre-rejestrowane)

- **Q-H2-PASS** → kolaps NIE jest wewnętrzną własnością gałęzi zdrowej
  z materią: odczyt (a) z NEEDS-N2 wsparty; kandydat na cykl łączony
  (małe amplitudy + źródło) i na dopisek core o roli materii — user-gate.
- **Q-H2-FAIL ∧ Q-H1-DEFORMATION** → kolaps ODPORNY na korpusowe sprzężenie
  w zbadanym zakresie λ̃: wzmocnione odczyty (b) fizyka substratu
  (kolaps→poziom 0 / nasycenie) lub (c) sygnał przeciw gałęzi jako
  nośnikowi struktur — interpretacja user-gated (NEEDS).
- **Q-H1-PULL** → nowa fizyka programu: materia jako mechanizm indukcji
  przejść przez progi — konfrontacja z hipotezą kreacji M911 (drzewo
  tamtego LOCKa) — user-gate.
- INCONCLUSIVE → NEEDS metodologiczny (zakres λ̃, forma ρ̂, t_max).

## 6. Forbidden moves (egzekwowane)

Rdzeń `.tex`/STATE.md/git NIETYKANE; katalogi innych cykli tylko odczyt;
formy M,𝒦,𝒰 i 𝒰_mat (po wyprowadzeniu P1-H1) bez modyfikacji; lista λ̃,
ρ̂, starty, progi, detektory, sponge niezmienialne po pierwszym biegu
produkcyjnym; ZAKAZ podłóg/barier innych niż pas klasyfikacyjny; ZAKAZ
dodawania pól dynamicznych (S05 — ρ jest źródłem zewnętrznym, statycznym);
ZAKAZ claimów o masach leptonów/oscylonach; λ̃<0 poza zakresem (ρ≥0);
INCONCLUSIVE ≠ pozytyw; correction note tylko dla błędu implementacji,
PRZED użyciem wyniku, pierwotne outputy zachowane.

## 7. Deliverables

`Phase_method_decisions.md` (FROZEN) · `engine_core.py` (kopia+cytat+człon
materii) · `Phase1_matter_analytic.py`+output · `Phase2_gate.py`+output ·
`Phase3_qh1_response.py`+output · `Phase3_qh2_stabilize.py`+output +
`Phase3_results/` (json+npz, verdict.json) · `Phase_FINAL_close.md` ·
`NEEDS.md` · dopis logu `README.md`. Cykl bez FINAL+NEEDS+README NIE jest
zakończony.
