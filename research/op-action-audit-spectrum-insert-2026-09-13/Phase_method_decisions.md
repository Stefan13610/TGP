---
title: "Phase_method_decisions — audyt jednej akcji: człon czasowy √−g·½K·|g^tt|ψ̇² ⟹ M(ψ)=ψ⁶/(4−3ψ)² (WYNIK Phase 1, tu tylko zapis formy i wariacji); odczyt B DZIEDZICZONY (𝒦=K_geo ψ⁴, w·g^ij≡1 — MD §2 poprzednika, nie rozstrzygany ponownie); 𝒰=w·V=γ(ψ⁴/4−ψ³/3); kaweat znakowanego g^tt raportowany RÓWNOLEGLE (deskryptywnie, tachion audytu ω²=k²−1); relaksacja radialna node-centered r_i=ih z więzem ψ(0)=A przez projekcję pinu po każdym kroku (residuum raportowane); pas ψ>4/3−1e−6 = BREAKDOWN-BOUNDARY; decyzje ZAMROŻONE przed startem obliczeń"
date: 2026-09-13
type: method-decisions
tgp_owner: research/op-action-audit-spectrum-insert-2026-09-13
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_method_decisions.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_FINAL_close.md]]"
  - "[[../op-metametric-boundary-2026-09-01/Phase_FINAL_close.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
  - "[[../../core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem cyklu)

Kryteria/progi/rodzina A/więz LOCKa (`Phase0_balance.md`) NIETKNIĘTE.
Czynności przed zamrożeniem: odczyt LOCKa W CAŁOŚCI, HANDOFF,
dokumentów dziedziczonych (MD §2 + FINAL close + silnik poprzednika
op-metric-pair-M911; FINAL close op-metametric-boundary), cytaty
z korpusu (sek08a, sek08c), weryfikacja środowiska importem
(CPython 3.14.2, numpy 2.4.3, scipy 1.17.1, sympy 1.14.0 — IDENTYCZNE
z oboma poprzednikami). Analiza wariacyjna §3–4 wykonana symbolicznie
na papierze PRZED kodem — część zamrożenia, nie obliczenie cyklu
(weryfikacja sympy w Phase 1 jest gate'em P1b/P1c).

**Rejestr WEJŚĆ (flagowany w każdym zależnym wyniku):**
- **K_geo = γ = c₀ = 1** [LOCK §1, bezwymiarowo];
- **rodzina amplitud A = {0.50, 0.70, 5/6, 7/6, 1.25, 1.30}** [LOCK §1,
  ZALOCKOWANA];
- **pudła R ∈ {60, 120}; siatki h ∈ {0.025, 0.0125}** [LOCK §1];
- **dt = 0.01, t_max = 200, stacjonarność ‖ψ̇‖∞ ≤ 1e−8** [LOCK §1];
- **starty deterministyczne:** ψ₀(r) = 1 + (A−1)·exp(−r²/(2σ²)),
  σ = 5.0 [INPUT-MD dziedziczony z poprzednika — LOCK: „gauss σ=5.0
  wokół 0 z ψ(∞)=1, amplituda dopasowana do A — jak bump poprzednika"];
  ψ₀(0) = A dokładnie; brak seeda (deterministycznie);
- **pas graniczny:** ψ_band = 4/3 − 1e−6; dolny kraniec deskryptywny
  ψ < 1e−6 (dziedziczone MD §4 poprzednika);
- **progi gate'ów:** P1c 1e−12; P3a dryf/ΔE ≤ 1e−10 (t=10);
  zbieżność h: |ΔE(h_c)−ΔE(h_f)|/max(|ΔE(h_f)|,1e−6) ≤ 5e−3
  (mianownik = |ΔE| siatki DROBNEJ, max z 1e−6 — doprecyzowanie
  implementacyjne litery LOCKa, zamrożone tu); R: |ΔE(120)−ΔE(60)|
  raportowane (brak członu objętościowego = warunek poprawności).

## 1. Formy korpusu — CYTATY (FROZEN; forbidden move a)

Z `core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex`:

1. **Akcja zunifikowana M9.1''** (eq:S-TGP-unified-M911-canonical,
   lin. ~363–373):
   > „S_TGP^{M9.1''}[Φ,ψ_m] = ∫ d⁴x (c₀ψ/(4−3ψ)) [ ½ K_geo ψ⁴
   > g_eff^{μν(M9.1'')} ∂_μψ ∂_νψ − V_{M9.1''}(ψ) + L_mat(ψ,ψ_m) ]"
   (sektor materii L_mat POZA modelem cyklu — czysty sektor
   grawitacyjny, LOCK §1). Forma bazowa: eq:S-TGP-unified (lin.
   ~148–156) z L_field = ½K(φ)g_eff^{μν}∂_μφ∂_νφ − V(φ)
   (eq:L-field, lin. ~160–166).
2. **Sprzężenie kinetyczne** (eq:K-coupling-unified, lin. ~169–172):
   > „K(φ) = K_geo φ⁴ (tw. thm:D-uniqueness, α = 2)."
3. **Potencjał** (prop:V-M911-canonical, eq:V-M911, lin. ~970–981):
   > „V_{M911}(ψ) = −(γ/12) ψ² (4−3ψ)²."
   Reguła sektorowa (przypis eq:V-selfinterference, lin. ~215–218):
   cykle gravity-related MUSZĄ używać V_{M9.1''}.
4. **Dowód prop:V-M911-canonical** (lin. ~984–1001): EOM statyczne
   „ψ'' + (2/r)ψ' + (2/ψ)(ψ')² = −(1/K(ψ)) dU_eff/dψ, gdzie
   U_eff(ψ) = ψV(ψ)/(4−3ψ)"; „U_eff(ψ) = γ(ψ⁴/4 − ψ³/3) + C".

Z `core/sek08c_metryka_z_substratu/sek08c_metryka_z_substratu.tex`:

5. **Pełna metryka** (eq:metric-M911-canonical, lin. ~420–424):
   > „ds² = −c₀²(4−3ψ)/ψ dt² + ψ/(4−3ψ) δᵢⱼ dxⁱ dxʲ, ψ ≡ Φ/Φ₀."
   ⟹ g_tt = −c₀²(4−3ψ)/ψ, **g^tt = −ψ/(c₀²(4−3ψ))**;
   g_ij = ψ/(4−3ψ)δᵢⱼ, **g^ij = (4−3ψ)/ψ δ^ij**.
6. **Element objętościowy** (rem:vol-element-M911,
   eq:vol-element-M911, lin. ~510–519):
   > „det(g_μν) = −c₀²ψ²/(4−3ψ)², √(−g_eff) = c₀ψ/(4−3ψ)."
   ⟹ **w(ψ) = ψ/(4−3ψ)**.

Zakaz modyfikacji wykładników/współczynników (LOCK §3a).

## 2. Odczyt kinetyki przestrzennej — DZIEDZICZONY (nie rozstrzygany)

**PRIMARY = odczyt B**, rozstrzygnięty CYTATEM w
`op-metric-pair-M911-2026-09-02/Phase_method_decisions.md` §2
(DZIEDZICZĘ per LOCK §1, nie rozstrzygam ponownie):
w(ψ)·g_eff^{ij} = [ψ/(4−3ψ)]·[(4−3ψ)/ψ]·δ^ij ≡ δ^ij ⟹
**𝒦(ψ) = K_geo ψ⁴** (bez w). Odczyt A (𝒦_A = w·K_geoψ⁴
= ψ⁵/(4−3ψ)) odnotowywany RÓWNOLEGLE tylko tam, gdzie różnicuje
(LOCK §2 Phase 2).

## 3. Człon czasowy akcji i wariacja (jawny zapis — LOCK §2 pkt 0)

**Człon czasowy (litera LOCKa §1, FROZEN):** gęstość kinetyczna
czasowa akcji = **√−g_eff · ½K(ψ) · |g_eff^{tt}| · ψ̇²**, z cytatami
form §1.5–1.6: |g^tt| = ψ/(c₀²(4−3ψ)), √−g = c₀ψ/(4−3ψ),
K = K_geoψ⁴. Stąd (do WYPROWADZENIA sympy w Phase 1 — M jest
WYNIKIEM, nie wejściem; tu tylko zapis formy):

  M(ψ) = √−g·K·|g^tt| = [c₀ψ/(4−3ψ)]·K_geoψ⁴·[ψ/(c₀²(4−3ψ))]
       = K_geo ψ⁶ / (c₀ (4−3ψ)²).

To ta sama waga czasowa, którą poprzednik zamroził w MD §9
(„M = diag(w(ψ₀)²·ψ₀⁴)"; w²ψ⁴ = ψ⁶/(4−3ψ)²) — spójność dziedziczenia.

**Gęstość Lagranżjanu (forma kanoniczna LOCKa §2 P1a):**

  L(ψ, ψ̇, ∇ψ) = ½M(ψ)ψ̇² − ½𝒦(ψ)|∇ψ|² − 𝒰(ψ),
  𝒦 = K_geoψ⁴ (odczyt B, §2), 𝒰 = w·V_{M911} = −γψ³(4−3ψ)/12
  = γ(ψ⁴/4 − ψ³/3) (tożsamość wielomianowa, weryfikacja sympy
  w Phase 1; dokładnie U_eff dowodu prop:V-M911-canonical, cytat §1.4
  — i DOKŁADNIE funkcjonał PRIMARY poprzednika: E[ψ] =
  ∫[½ψ⁴|∇ψ|² + 𝒰]dx, gate P1b).

**Wariacja (standardowa Lorentzowska, pole ψ; metryka podstawiona
jako funkcja ψ PRZED wariacją — forma złożona; g_eff NIE jest
zmienną niezależną, wariacja NIE jest euklidesowa ⟹ warunek STOP
z LOCKa §2 Phase 1 nie aktywuje się):**

  π = ∂L/∂ψ̇ = M(ψ)ψ̇;
  H[π,ψ] = ∫ [ π²/(2M(ψ)) + ½𝒦(ψ)|∇ψ|² + 𝒰(ψ) ] d³x;
  EOM (Euler–Lagrange, z członami M′, 𝒦′):
  M ψ̈ + ½M′ ψ̇² = ∇·(𝒦∇ψ) − ½𝒦′|∇ψ|² − 𝒰′  = −δE_PRIMARY/δψ.

  δH/δψ = −π²M′/(2M²) − ∇·(𝒦∇ψ) + ½𝒦′|∇ψ|² + 𝒰′;
  przy π=0: δH/δψ|_{π=0} = −∇·(𝒦∇ψ) + ½𝒦′|∇ψ|² + 𝒰′ — tożsamość
  z δE_PRIMARY/δψ do weryfikacji sympy simplify=0 (gate P1b,
  FAIL ⟹ STOP).

**Kaweat znaku (odnotowany PRZED obliczeniami; raport RÓWNOLEGŁY,
deskryptywny, NIE werdyktotwórczy):** metryka M9.1'' jest zapisana
w sygnaturze (−,+,+,+) (cytat §1.5), więc znakowany g^tt < 0;
kontrakcja literalnego zapisu +½K g^{μν}∂_μψ∂_νψ ze ZNAKOWANYM
g^{μν} dałaby człon czasowy ujemny (L ∝ −½Mψ̇² + ½𝒦|∇ψ|² − 𝒰),
czyli po linearyzacji ω² = k² − 1 — dokładnie tachion audytu
(rozdz. 2, stara hybryda). LOCK §1 zamraża odczyt członu czasowego
przez **|g_eff^{tt}|** (równoważnie: konwencja cząstkowa (+,−,−,−)
dla zapisu ½Kg^{μν}∂ψ∂ψ, standardowa dla L = T − V), a §2 P1a
zamraża formę kanoniczną L = ½Mψ̇² − ½𝒦|∇ψ|² − 𝒰. Oba odczyty
znaku będą policzone i zraportowane w Phase 1–2; **werdykt Q-D1
wg litery LOCKa używa M z zamrożonej formy |g^tt|**. Konsekwencja
uboczna do odnotowania deskryptywnie: statyczna granica EOM formy
kanonicznej daje δE_PRIMARY/δψ=0 (RHS (ψ−1)/ψ² po podzieleniu
przez 𝒦), podczas gdy dowód prop:V-M911-canonical wyprowadza R3 ODE
z RHS (1−ψ)/ψ² (odpowiada wariacji ze znakowanym g^{μν}) — napięcie
znaku na poziomie rdzenia, raportowane wprost w Phase 1/FINAL
(user-gate w NEEDS; ZERO samodzielnych napraw).

## 4. Phase 3 — schemat relaksacji z więzem ψ(0)=A (FROZEN)

**Funkcjonał (PRIMARY, identyczny z poprzednikiem):**
E[ψ] = 4π ∫₀^R r² [ ½𝒦(ψ)(∂_r ψ)² + 𝒰(ψ) ] dr, 𝒦=ψ⁴,
𝒰=γ(ψ⁴/4−ψ³/3) — postać WIELOMIANOWA dokładna (zero dzielenia
przez (4−3ψ) w silniku); 𝒰′=γψ²(ψ−1), 𝒦′=4ψ³.

**Dyskretyzacja (adaptacja silnika radialnego poprzednika pod pin
w r=0; HANDOFF: 4π, trapez, ghost points):**
- siatka WĘZŁOWA r_i = i·h, i = 0..N, N = R/h (węzeł r=0 istnieje ⟹
  więz ψ(0)=A literalnie na węźle 0); midpointy r_{i+1/2}=(i+½)h;
- energia dyskretna: E_h = 4π[ Σ_{i=0}^{N−1} h·r²_{i+1/2}·
  ½𝒦(ψ_{i+1/2})·((ψ_{i+1}−ψ_i)/h)² + Σ_{i=0}^{N} c_i·h·r_i²·𝒰(ψ_i) ],
  ψ_{i+1/2}=(ψ_i+ψ_{i+1})/2, wagi trapezu c_0=c_N=½, inaczej 1;
- gradient dyskretny dH_i = ∂E_h/∂ψ_i (struktura strumieniowa
  dziedziczona: t_flux = r²_{i+1/2}𝒦(mid)Δψ/h, t_quad =
  ¼h·r²_{i+1/2}𝒦′(mid)(Δψ/h)², człon lokalny c_i·h·r_i²·𝒰′);
- **flow:** ψ̇_i = −dH_i/(4π·m_i), masa węzła m_i = c_i·h·r_i²
  dla i≥1; **m_0 = h³/24** (dokładna ∫₀^{h/2} r²dr — waga trapezu
  ½h·r₀² = 0 czyni węzeł 0 bezmasowym; m_0 to parametryzacja
  PRĘDKOŚCI flow, nie zmienia E_h ani zbioru stanów stacjonarnych;
  w biegach z pinem węzeł 0 i tak jest nadpisywany);
- **brzegi:** jednorodny Neumann w 0 i R przez ghost points
  (ψ_{−1}=ψ_1, ψ_{N+1}=ψ_{N−1}) — realizowany naturalnie przez
  strukturę strumieniową (zerowy strumień przez r=0 i r=R);
- **krok czasowy (dziedziczony FROZEN):** semi-implicit Euler
  (I − dt·A_t·L)(ψ^{n+1}−ψ^n) = dt·rhs(ψ^n), A_t = 1.05·max 𝒦(ψ^n),
  L = dyskretny operator dyfuzyjny (Lu)_i = [a_i(u_{i+1}−u_i)
  − a_{i−1}(u_i−u_{i−1})]/m_i, a_i = r²_{i+1/2}/h (a_{−1}=a_N=0),
  solve_banded; dt=0.01;
- **więz pinu (LOCK §1, jedyny więz):** po KAŻDYM kroku
  **projekcja: ψ_0 ← A** (nadpisanie); residuum pinu
  |ψ_0^{pre-projekcja} − A| raportowane przy każdej próbce (co Δt=1)
  i na końcu biegu; w biegu P3a bez więzu projekcja WYŁĄCZONA;
- **stacjonarność:** ‖ψ̇‖∞ ≤ 1e−8 co Δt=1, liczona po projekcji
  (węzeł pinowany ma ψ̇=0 z projekcji; max po węzłach swobodnych);
  t_max=200;
- **ΔE_insert(A;R,h) = E_h[ψ_A^relax] − E_h[ψ≡1]** na IDENTYCZNEJ
  siatce/pudle/brzegu (ta sama E_h, ta sama tablica r).

**Obsługa granic dziedziny (FROZEN; zero podłóg/barier — forbidden
move c):** ZERO członów dodanych do E. Pas ψ > 4/3−1e−6 (sprawdzany
każdy krok na max ψ) ⟹ **BREAKDOWN-BOUNDARY** (stop, klasyfikacja
deskryptywna, ostatni stan zapisany, NIE werdykt pozytywny);
min ψ < 1e−6 ⟹ **BREAKDOWN-BOUNDARY-LOWER**; niefinityczność ⟹
**BREAKDOWN** (obsługa błędów dziedziczona z poprzednika,
z przechwyceniem wyjątków solvera).

**P3a (bramka, FROZEN):** (i) próżnia ψ≡1 BEZ więzu, t=10, wszystkie
4 pudła (R×h): dryf max_t‖ψ−1‖∞ ≤ 1e−10; (ii) relaksacja z pinem
A=1 (start gauss z A=1 ⟹ ψ₀≡1): ΔE_insert = 0 ± 1e−10, wszystkie
4 pudła. FAIL ⟹ STOP.

**P3b (macierz, FROZEN):** 6 A × 2 R × 2 h = 24 relaksacje;
werdykt per LOCK §2 (Q-D2-COST / Q-D2-CHANNEL / Q-D2-INCONCLUSIVE
wg litery); deskryptywnie obowiązkowo: tabela ΔE_insert(A)
+ monotonia + profile (rozciągłość r_half: największe r, gdzie
|ψ−1| ≥ ½|A−1|; ψ_max/ψ_min) + odniesienie definicyjne do liczb
Q1-POS poprzednika (−0.179 / +16156.6): wspólne tło vs różne tła,
NIE reprodukcja. Wyniki per bieg json (`Phase3_results/`), output
zbiorczy `Phase3_output.txt`. ZAKAZ wnioskowania o barierze kreacji
z ΔE_insert (forbidden move h).

## 5. Phase 1–2 — specyfikacja (FROZEN)

- **Phase 1 (`Phase1_canonical.py`, sympy, zero numeryki flow):**
  P1a: M, 𝒦, 𝒰 wyprowadzone symbolicznie z form §1 (M = WYNIK:
  simplify(√−g·K·|g^tt|)); π, H, EOM z członami M′, 𝒦′ (jawnie);
  P1b: tożsamość δH/δψ|_{π=0} ≡ δE_PRIMARY/δψ, simplify=0 —
  FAIL ⟹ STOP; dodatkowo tożsamość 𝒰 = wV (wielomianowo) simplify=0;
  P1c: M, 𝒦, 𝒰 w ψ∈{0.5, 1, 7/6, 1.3} sympy exact→float vs float
  implementacji silnika — zgodność 1e−12 (osiągalny FAIL).
  Deskryptywnie: wariant znakowany g^tt (kaweat §3) — L, EOM,
  statyka vs R3 ODE; NIE wpływa na gate'y.
- **Phase 2 (`Phase2_dispersion.py`):** linearyzacja EOM wokół ψ*=1
  (ψ = 1 + ε·e^{i(kx−ωt)}, sympy, wyprowadzić nie postulować) ⟹
  ω²(k) = [𝒦(1)k² + 𝒰″(1)]/M(1); m² = 𝒰″(1)/M(1),
  c_s² = 𝒦(1)/M(1); znaki M(ψ), 𝒦(ψ) na (0,4/3) (sympy solveset /
  analiza czynników + kontrola float na siatce dziedziny);
  odczyt A równolegle (𝒦_A = ψ⁵/(4−3ψ); w ψ*=1: 𝒦_A(1)=1 —
  odnotować gdzie różnicuje); wariant znakowany g^tt deskryptywnie
  (ω² = k² − 1). **Werdykt Q-D1 wg litery LOCKa §2** (M(1)>0 ∧
  𝒰″(1)>0 ∧ ω²(k)>0 ∀k≥0 ∧ M>0, 𝒦>0 na (0,4/3)).

## 6. Higiena wykonania

Pełne ścieżki bez `cd`; `ls` po każdym zapisie (artefakt zagnieżdżeń
`TGP/TGP_v1/TGP/...`); sandbox: bez /dev/null i heredoc — wszystko
przez pliki skryptów; batch Phase 3 w tle z logiem
(`Phase3_results_batch.log`), aktywne czekanie; outputy do
`Phase{1,2,3}_output.txt`; wyniki per bieg json w `Phase3_results/`.
Rdzeń `.tex`, STATE.md, git, katalogi innych cykli — NIETYKANE.
INCONCLUSIVE/BREAKDOWN-BOUNDARY ≠ pozytyw. Korekty wyłącznie dla
udokumentowanego błędu implementacji — correction_note PRZED użyciem
wyniku, pierwotne outputy zachowane.

**FROZEN 2026-09-13, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
