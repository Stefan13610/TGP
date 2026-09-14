---
title: "Phase_method_decisions — para metryczna (w, V_M9.1'', K=ψ⁴): rozjazd odczytów rozstrzygnięty CYTATEM (kinetyka BEZ w, bo w·g^ij(M9.1'')≡1 — eq:S-TGP-unified-M911-canonical + analogia eq:S-static-correct); 𝒰=w·V=γ(ψ⁴/4−ψ³/3) tożsamość z U_eff dowodu prop:V-M911-canonical; pas graniczny ψ>4/3−1e−6 = BREAKDOWN-BOUNDARY (klasyfikacja, zero barier); detektory w ψ (dolny <5/6, górny >7/6); σ bumpa=5.0 [INPUT-MD]; skalowanie sieci per-siatka do ψ_max=1.30; decyzje ZAMROŻONE przed startem obliczeń"
date: 2026-09-02
type: method-decisions
tgp_owner: research/op-metric-pair-M911-2026-09-02
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-metric-closure-relaxation-2026-09-02/Phase_method_decisions.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem cyklu)

Kryteria/progi/starty LOCKa (`Phase0_balance.md`) NIETKNIĘTE. Czynności
przed zamrożeniem: odczyt LOCKa, HANDOFF, dokumentów dziedziczonych
(MD + silnik + FINAL close poprzednika), cytaty z korpusu (sek08a),
weryfikacja środowiska importem (CPython 3.14.2, numpy 2.4.3,
scipy 1.17.1, sympy 1.14.0 — IDENTYCZNE z poprzednikiem), weryfikacja
npz tła 2π (mtime 2026-08-31 21:41:07 — zgodny z zapisem poprzednika,
READ-ONLY; klucze `2pi__A1.0__N{32,48}`; g∈[0.6064,1.4734] (N=32),
[0.6100,1.4688] (N=48) ⟹ ψ_max=g_max²≈2.171/2.157 — zgodne z uwagą
LOCKa „≈2.17"). Analiza wariacyjna §3 wykonana symbolicznie na papierze
PRZED kodem — część zamrożenia, nie obliczenie cyklu (weryfikacja
sympy w Phase 1).

**Rejestr WEJŚĆ (flagowany w każdym zależnym wyniku):**
- **K_geo = 1, γ = 1** [LOCK §1, bezwymiarowo];
- **seed = 20260903, amp szumu = 1e−3** [LOCK §1, start (i)];
- **ψ_max startów (ii)/(iii) = 1.30** [LOCK §1];
- **progi detektorów w ψ:** dolny ψ_thr_dn = 5/6 = 0.8333333, górny
  ψ_thr_up = 7/6 = 1.1666667 [LOCK §1, zalockowane];
- **pas graniczny:** ψ_band = 4/3 − 1e−6 [LOCK §1];
- **σ bumpa radialnego = 5.0** [INPUT-MD — LOCK specyfikuje amplitudę
  (ψ_max=1.30), geometrię (gauss, R=60, h∈{0.025,0.0125}), NIE
  szerokość; zamrażam σ=5.0: ≫h (≥400 komórek na σ), ≪R (12σ=60 —
  ogon w 1+3e−4·amp na brzegu Neumanna), skala rdzenia solitonowego
  linii ~O(10)];
- **skalowanie sieci (iii), procedura FROZEN:** ψ_raw = g² (mapowanie
  LOCKa), ψ₀ = 1 + s·(ψ_raw − 1) z s = 0.30/(max ψ_raw − 1)
  per siatka (warunek LOCKa ψ_max=1.30 DOKŁADNIE na każdej siatce;
  tła N=32/48 są niezależnymi stanami zrelaksowanymi, nie próbkami
  jednego pola — skala per-siatka; oryginalny kształt (min/max/s)
  raportowany obowiązkowo);
- siatki: geneza L=4π N∈{48,64}; bump radialny R=60 h∈{0.025,0.0125};
  sieć L=2π N∈{32,48} [LOCK §1];
- dt = 0.01 (kontrola dt/2 = 0.005), t_max = 200, stacjonarność
  ‖ψ̇‖∞ ≤ 1e−8 próbkowana co Δt=1 [LOCK §2 + dziedziczone];
- tło 2π: `../op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz`
  (READ-ONLY; mtime weryfikowany po KAŻDYM odczycie).

## 1. Formy korpusu — CYTATY (FROZEN; forbidden move a)

Z `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`:

1. **Potencjał** (prop:V-M911-canonical, eq:V-M911, lin. ~970–981):
   > „unikalnym (mod stała) potencjałem zapewniającym że static
   > spherically-symmetric EOM jest tożsame z R3 ODE jest:
   > V_{M911}(ψ) = −(γ/12) ψ² (4−3ψ)²."
2. **Element objętościowy / metryka** (rem:hyp-unified-action-M911-canonical,
   lin. ~349–360):
   > „Element objętościowy: √(−g_eff) = c₀ ψ/(4−3ψ) (M9.1'' canonical,
   > eq:vol-element-M911, sek08c body; A2 audit closure 2026-05-02)."
   > „Metryka efektywna: M9.1'' hiperboliczna ds² = −c₀²(4−3ψ)/ψ dt²
   > + ψ/(4−3ψ) δᵢⱼ dxⁱdxʲ (eq:metric-M911-canonical, sek08c body)."
   ⟹ **w(ψ) = ψ/(4−3ψ)**; przestrzenna inwersja metryki
   **g_eff^{ij} = (4−3ψ)/ψ · δ^{ij}**.
3. **Sprzężenie kinetyczne** (eq:K-coupling-unified, lin. ~169–172):
   > „K(φ) = K_geo φ⁴ (tw. thm:D-uniqueness, α = 2)."
4. **Akcja zunifikowana M9.1''** (eq:S-TGP-unified-M911-canonical,
   lin. ~363–373):
   > „S_TGP^{M9.1''} = ∫ d⁴x (c₀ψ/(4−3ψ)) [ ½ K_geo ψ⁴
   > g_eff^{μν(M9.1'')} ∂_μψ ∂_νψ − V_{M9.1''}(ψ) + L_mat ]."
5. **Reguła sektorowa** (przypis eq:V-selfinterference, lin. ~215–218):
   > „Reguła użycia (post-2026-05-09): cykle gravity-related (M9.1''
   > metryka, R3 ODE, Newton limit, PPN) MUSZĄ używać V_{M9.1''}."

Zakaz modyfikacji wykładników/współczynników (LOCK §3a). Sektor
materii (L_mat, V_orig) POZA modelem cyklu (LOCK §1: czysty sektor
grawitacyjny).

## 2. Rozjazd odczytów „czy w mnoży człon kinetyczny" — rozstrzygnięcie CYTATEM (PRZED startem)

**Odczyt A (literalny zapis LOCKa §1):** E[ψ] = ∫ w(ψ)·[½K_geo ψ⁴|∇ψ|²
+ V(ψ)] dx — czynnik w mnoży TAKŻE człon kinetyczny; brak czynnika
inwersji metryki. Efektywny współczynnik kinetyczny:
𝒦_A(ψ) = w·K_geo ψ⁴ = ψ⁵/(4−3ψ).

**Odczyt B (kanoniczny, z akcji):** w eq:S-TGP-unified-M911-canonical
(cytat §1.4) gradient jest kontraktowany z **g_eff^{μν(M9.1'')}**;
dla konfiguracji statycznej część przestrzenna: w(ψ)·g_eff^{ij}
= [ψ/(4−3ψ)]·[(4−3ψ)/ψ]·δ^{ij} = δ^{ij} — czynniki się DOKŁADNIE
kasują. Efektywny współczynnik kinetyczny: **𝒦(ψ) = K_geo ψ⁴** (bez w).

**Rozstrzygnięcie: PRIMARY = odczyt B.** Podstawa (cytaty):
- eq:S-TGP-unified-M911-canonical zawiera jawnie g_eff^{μν(M9.1'')}
  w członie kinetycznym (§1.4);
- dokładna analogia w korpusie dla formy I (lin. ~605–616, nota przy
  eq:S-unified-psi — HANDOFF wskazuje ten fragment jako rozstrzygający):
  > „efektywną miarę kinetyczną φ·φ⁴ = φ⁵. Gradient pola używa inwersji
  > metryki: g_eff^{ij} = φ^{−1}δ^{ij}, co daje dodatkowy czynnik φ^{−1},
  > stąd efektywny człon kinetyczny w całce wynosi
  > φ⁵·φ^{−1}·(∇ψ)²/2 = φ⁴(∇ψ)²/2."
  i skorygowana akcja statyczna eq:S-static-correct (lin. ~619–625):
  kinetyka ½ψ⁴(∇ψ)² BEZ czynnika objętościowego, potencjał Z czynnikiem
  (ψV). Identyczna kompozycja w M9.1'' daje 𝒦 = K_geo ψ⁴;
- dowód prop:V-M911-canonical (lin. ~984–1000) używa DOKŁADNIE tej
  struktury: EOM „ψ'' + (2/r)ψ' + (2/ψ)(ψ')² = −(1/K(ψ)) dU_eff/dψ,
  gdzie U_eff(ψ) = ψV(ψ)/(4−3ψ)" — tj. potencjał efektywny = w·V,
  współczynnik kinetyczny = K(ψ)=ψ⁴ (bez w).

**Odczyt A odnotowany, NIE realizowany** (byłby powtórzeniem błędu
klasy „hybryda" poprzednika: forma spoza akcji korpusu). Konsekwencja
dla kryteriów: ρ_eff = w·V jest WSPÓLNE dla obu odczytów (Q-A
identyczne); 𝒦_A i 𝒦 są oba dodatnie na (0,4/3) (znak współczynnika
kinetycznego — P1a — raportowany dla OBU odczytów).

## 3. Funkcjonał PRIMARY i gradient flow (FROZEN)

**E[ψ] = ∫ [ ½ 𝒦(ψ)|∇ψ|² + 𝒰(ψ) ] dx**, dziedzina ψ∈(0,4/3),
- 𝒦(ψ) = K_geo ψ⁴ = ψ⁴,
- **𝒰(ψ) = ρ_eff(ψ) = w(ψ)·V_{M911}(ψ) = −γψ³(4−3ψ)/12
  = γ(ψ⁴/4 − ψ³/3)** — tożsamość algebraiczna (biegun w kasuje się
  z podwójnym zerem V; dokładnie U_eff z dowodu prop:V-M911-canonical,
  lin. ~999: „U_eff(ψ) = γ(ψ⁴/4 − ψ³/3) + C"; weryfikacja sympy
  w Phase 1 jako część gate'u P1b). Implementacja numeryczna używa
  postaci WIELOMIANOWEJ (dokładnej) — zero dzielenia przez (4−3ψ)
  w silniku.
- Pochodne: 𝒦′ = 4ψ³; **𝒰′ = γψ²(ψ−1)** = w′V + wV′ (człony w′ i K′
  jawnie obecne w wyprowadzeniu; w′(ψ) = 4/(4−3ψ)²,
  V′ = −γψ(4−3ψ)(2−3ψ)/3 — zera V′ w {0, 2/3, 4/3}, stąd średnia
  progowa detektora dolnego 5/6 w LOCKu).
- **Wariacja:** δE/δψ = −∇·(𝒦∇ψ) + ½𝒦′|∇ψ|² + 𝒰′.
  **Gradient flow: ∂ψ/∂t = −δE/δψ.**
- **Stacjonarność próżni (analitycznie, pre-compute):** 𝒰′(1)
  = w′(1)V(1) + w(1)V′(1) = 4·(−1/12) + 1·(1/3) = 0 — ψ=1 JEST punktem
  krytycznym pełnego E (jednorodnie człony gradientowe znikają).
  Referencja próżniowa NIE jest potrzebna (w mnoży V w CAŁOŚCI;
  problem w′·U(1)≠0 poprzednika tu nie występuje). P2a używa ψ*=1
  warunkowo na potwierdzeniu w P1a (litera LOCKa).
- **Dyskretyzacja (dziedziczona dosłownie, silnik poprzednika):**
  struktura strumieniowa gradientu dyskretnej energii E_h —
  t_flux = 𝒦(ψ_mid)Δψ/h, t_quad = ¼𝒦′(ψ_mid)(Δψ/h)², człon lokalny
  h·𝒰′(ψ); radialnie miara r²dr + Neumann; 3D periodycznie.
  Podmiana (Keff,Keffp,Ueff,Ueffp) → (ψ⁴, 4ψ³, 𝒰, 𝒰′).
- **Krok czasowy (dziedziczony FROZEN):** semi-implicit Euler
  (I − dt·A_t·L)(ψ^{n+1}−ψ^n) = dt·rhs(ψ^n), A_t = 1.05·max 𝒦(ψ^n);
  L = laplasjan (3D rFFT / radialnie solve_banded). dt=0.01,
  kontrola dt/2; t_max=200; stacjonarność ‖rhs‖∞≤1e−8 co Δt=1.

## 4. Obsługa granic dziedziny (FROZEN; zero podłóg/barier — sedno Q-A)

- **ZERO członów dodanych do E** (forbidden move b): żadnych podłóg,
  sufitów, kar. Formy §3 są dokładne (wielomianowe) na całej osi —
  regularizacja bieguna w jest ZBĘDNA (biegun kasuje się analitycznie
  w 𝒰=wV; to własność pary korpusowej, nie decyzja).
- **Pas graniczny (LOCK §1):** wejście pola w pas ψ > ψ_band
  = 4/3 − 1e−6 (sprawdzane KAŻDY krok na max ψ) ⟹ status
  **BREAKDOWN-BOUNDARY** — bieg zatrzymany (poza pasem metryka M9.1''
  nie istnieje: g_tt zmienia znak; kontynuacja używałaby ekstrapolacji
  poza domenę ważności), ostatni stan zapisany, klasyfikacja
  deskryptywna („pole wybiera granicę"), NIE werdykt pozytywny.
- **Dolny kraniec (symetrycznie, deskryptywnie):** min ψ < 1e−6 ⟹
  **BREAKDOWN-BOUNDARY-LOWER** (dziedzina (0,4/3); klasyfikacja,
  zero sił dodanych). Niefinityczność ⟹ BREAKDOWN (numeryczne;
  obsługa błędów dziedziczona).

## 5. Detektory (FROZEN; reguły dziedziczone, progi w ψ z LOCKa)

- **Dolny:** maska m_dn = (ψ < 5/6); **górny:** maska m_up = (ψ > 7/6).
- Maszyneria dziedziczona DOSŁOWNIE: `scipy.ndimage.label`
  (3D 6-spójność + periodyczne sklejanie union-find przez pary ścian;
  radialnie 1D bez periodyczności); próbkowanie co Δt=1; N_seed_{dn,up}
  z t=0; NUKLEACJA_{dn,up}: pierwsze t₀ z N(t) > N_seed utrzymane
  dla wszystkich t∈{t₀,…,t₀+10} (≥10 j.cz. = 11 próbek); oba detektory
  równolegle, zdarzenie = pierwszy potwierdzony (kierunek raportowany).
- **Zbieżność (Q-B-PASS-NUCLEATION):** ten sam kierunek na OBU siatkach
  startu ORAZ w biegach dt/2 na obu siatkach; N_det zgodne ±1 między
  czterema biegami pary.
- Detektory NIEZMIENIALNE po pierwszym biegu (forbidden move c).

## 6. Starty i macierz biegów Phase 3 (FROZEN)

**Starty (LOCK §1, konstrukcje dziedziczone):**
- (i) **geneza:** ψ₀ = 1 + f; szum pasmowy dziedziczony verbatim
  (|n_i|≤8, w L=4π k=n/2≤4): rng(**20260903**) losuje
  standard_normal((17,17,17,2)), hermityzacja
  C_sym=(C+conj(C[::-1,::-1,::-1]))/2, TE SAME współczynniki wbudowane
  w N=48 i N=64, normalizacja max|f|=1e−3 ze zbudowanego N=64;
- (ii) **bump radialny:** ψ₀(r) = 1 + 0.30·exp(−r²/(2σ²)), σ=5.0
  [INPUT-MD], R=60, h∈{0.025,0.0125} (ψ_max=1.30 dokładnie; Neumann);
- (iii) **sieć 2π:** klucze `2pi__A1.0__N{32,48}`; ψ_raw=g²,
  ψ₀ = 1 + s·(ψ_raw−1), s=0.30/(max ψ_raw−1) per siatka (§ Rejestr);
  oryginalny kształt raportowany (min/max g, min/max ψ_raw, s).

**Macierz: 6 biegów głównych dt=0.01:**
gen×N∈{48,64}; bump×h∈{0.025,0.0125}; lat×N∈{32,48}.
**Biegi dt/2 „przy zdarzeniach" (LOCK §2):** warunkowo dla pary
(start), której DOWOLNY bieg główny kończy się NUCLEATION lub
BREAKDOWN-BOUNDARY (zbieżność werdyktu/kategorii); dla par czysto
STATIONARY/FAIL zbieżność ocenia porównanie siatkowe (≤5e−3,
wspólna podsiatka 16³ / interpolacja fine→coarse — dziedziczone;
geneza stride 3/4, sieć stride 2/3, radialnie interp).

**Werdykt Q-B (litera LOCKa §2, składanie dziedziczoną logiką):**
- Q-B-PASS-NUCLEATION: nukleacja zbieżna (obie siatki + dt/2, ±1)
  w ≥1 parze;
- inaczej Q-B-PASS-STATIC: para STATIONARY z min ‖ψ−const‖∞ ≥ 0.05
  (const = średnia pola; dev_const = ½(max−min) — dziedziczone)
  i zbieżnością podsiatkową ≤5e−3; sanity: ψ_max vs 4/3 raportowane;
- Q-B-FAIL: WSZYSTKIE pary STATIONARY jednorodne (dev<0.05) — raport
  wartości końcowej vs ψ*=1;
- Q-B-INCONCLUSIVE: pozostałe; BREAKDOWN-BOUNDARY zbieżny (obie
  siatki, ta sama klasyfikacja) = osobna kategoria deskryptywna.
- Deskryptywnie OBOWIĄZKOWO: los startu (iii) — czy struktura sieci
  przeżywa w sektorze metrycznym.

## 7. Phase 1 — specyfikacja (FROZEN)

- **P1a (sympy):** ρ_eff = w·V na (0,4/3): punkty krytyczne (solve
  ρ′=0), krzywizna ρ″ w minimach, wartość w próżni, granice 0⁺ i 4/3⁻,
  globalne minimum na dziedzinie; ograniczoność E z dołu (kinetyka
  ≥0 + inf ρ_eff > −∞ ⟹ E ≥ |Ω|·inf ρ_eff); znak współczynnika
  kinetycznego na (0,4/3) dla OBU odczytów (𝒦=ψ⁴, 𝒦_A=ψ⁵/(4−3ψ));
  stacjonarność ψ=1 w pełnym E (𝒰′(1)=0 sympy). Werdykt Q-A wg litery
  LOCKa: minimum ψ*∈(0,4/3) z ρ″>0 ORAZ ρ_eff(4/3)>ρ_eff(ψ*) ORAZ
  E ograniczone z dołu ⟹ Q-A-PASS; przeciwnie Q-A-FAIL (raport wprost,
  Phase 2–3 NADAL wykonywane z flagą — drzewo §5).
- **P1b (gate zgodności, osiągalny FAIL):** ρ_eff w {0.5, 1, 7/6, 1.3}
  sympy (exact→float) vs implementacja float silnika — zgodność 1e−12;
  dodatkowo tożsamość wielomianowa 𝒰 = wV i 𝒰′ = w′V+wV′ (sympy
  simplify == 0) jako warunek użycia postaci wielomianowej w silniku.
- Wykresy kontrolne: tabela ASCII ρ_eff/𝒰′/𝒦 na siatce dziedziny
  + `Phase1_landscape.png` (matplotlib, jeśli dostępny; brak = tylko
  tabela, odnotowane).

## 8. Phase 2 — specyfikacja (FROZEN)

- **P2a (próżnia zostaje):** start DOKŁADNIE ψ≡ψ* (=1 z P1a, litera
  LOCKa: jeśli P1a wykaże inaczej — ψ* z P1a, odnotować), bez zaburzeń,
  t=10, gate max_t ‖ψ−ψ*‖∞ ≤ 1e−10; trzy geometrie: radialna h=0.0125,
  3D L=2π N=32, 3D L=4π N=48.
- **P2b (detektory, osiągalny FAIL):** pola testowe analityczne
  (bez flow): (a) próżnia + zasiany obiekt DOLNY: dip gaussowski
  ψ = 1 − 0.4·exp(−r²/2σ_t²) (min 0.6 < 5/6), oczekiwane
  N_dn=1±0, N_up=0; (b) + zasiany obiekt GÓRNY: bump ψ = 1
  + 0.3·exp(−r²/2σ_t²) (max 1.3 > 7/6), oczekiwane N_up=1±0, N_dn=0;
  (c) czysta próżnia ψ≡1: zero alarmów obu detektorów. σ_t=1.0;
  geometrie: 3D L=4π N=48 (obiekt w centrum pudła) ORAZ radialna
  h=0.025 (obiekt w r=0) — 6 testów.
- **FAIL któregokolwiek ⟹ STOP (litera LOCKa §2).**

## 9. Phase 4 — specyfikacja warunkowa (FROZEN; tylko Q-B-PASS-STATIC)

Druga wariacja pełnego E wokół stanu ψ₀:
δ²E[φ] = ∫[𝒦|∇φ|² + 2𝒦′φ∇ψ₀·∇φ + (½𝒦″|∇ψ₀|² + 𝒰″)φ²]dx
(𝒦″=12ψ², 𝒰″=γψ(3ψ−2)); implementacja: DOKŁADNY Hessian dyskretnej
E_h (dziedziczona konstrukcja rzadka). Problem uogólniony
H_h φ = ω² M φ; **waga czasowa konsekwentnie z akcją i odczytem B:**
gęstość kinetyczna czasowa = ½·w·K·|g_eff^{tt}|·ψ̇² z |g^{tt}|
= ψ/(c₀²(4−3ψ)) ⟹ **M = diag(w(ψ₀)²·ψ₀⁴)** (c₀=1); w² dokładne
(bieg PASS-STATIC ma z definicji ψ<4/3; stan z pasem granicznym nie
jest PASS-STATIC). eigsh (najmniejsze), v0 deterministyczne
rng(20260903). Mody zerowe (translacyjne) identyfikowane PRZED
interpretacją (overlap z ∂ψ₀/∂x_i). Zbieżność siatkowa
≤0.05·max(|ω²_min|,0.1); **Q-C-PASS: ω²_min ≥ −1e−3.**
Przy Q-B-PASS-NUCLEATION zamiast widma: charakterystyka kaskady BEZ
progów (N_obj(t) obu detektorów, przyrosty/j.cz., rozkłady rozmiarów
w t₀, t₀+5, t₀+10 — forma dziedziczona).

## 10. Higiena wykonania

Pełne ścieżki bez `cd`; `ls` po każdym zapisie; runy >10 min w tle
z checkpointami npz co 10 j.cz. (proces ≤ ~55 min, wznawialne
`--resume`); wyniki per-bieg json/npz w `Phase3_results/`; `--verdict`
składa do `Phase3_output.txt` + `Phase3_relaxed_states.npz`.
INCOMPLETE raportowane z przyczyną. Rdzeń `.tex`, STATE.md, git,
katalogi innych cykli — NIETYKANE (npz READ-ONLY, mtime po każdym
odczycie). INCONCLUSIVE/BREAKDOWN-BOUNDARY ≠ pozytyw. Outputy do
`Phase{1,2,3,4}_output.txt`.

**FROZEN 2026-09-02, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
