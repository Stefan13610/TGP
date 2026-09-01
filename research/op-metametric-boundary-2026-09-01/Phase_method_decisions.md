---
title: "Phase_method_decisions — mapowanie g_floor z dodatekB (z cytatem równań), gładka kara C² zachowująca strukturę wariacyjną (κ=100 FROZEN), stabilizowany semi-implicit gradient flow (dt=0.01), detektor nukleacji ndimage.label z periodycznym sklejaniem; decyzje ZAMROŻONE przed startem obliczeń"
date: 2026-09-01
type: method-decisions
tgp_owner: research/op-metametric-boundary-2026-09-01
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../../axioms/substrat/dodatekB_substrat.tex]]"
  - "[[../op-substrate-fluctuation-channel-2026-08-23/Phase_FINAL_close.md]]"
  - "[[../op-3d-canonical-lattice-2026-08-31/Phase_method_decisions.md]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem cyklu)

Kryteria/progi/punkty LOCKa (`Phase0_balance.md`) NIETKNIĘTE. Jedyne
obliczenia wykonane przed zamrożeniem: arytmetyka transkrypcyjna rejestru
wejść — pierwiastki trzech zalockowanych progów QB-2 (tabela §1 poniżej)
oraz odczyt kluczy `Phase2_backgrounds3d.npz` (weryfikacja dostępności
tła 2π; READ-ONLY). Środowisko sprawdzone importem: CPython 3.14.2,
numpy 2.4.3, scipy 1.17.1, sympy 1.14.0 — identyczne z op-3d-canonical.

**Rejestr WEJŚĆ (flagowany w każdym zależnym wyniku):**
- kotwica radialna λ_min = −1.646589 [INPUT, op-3d-canonical Phase 1,
  h=0.0125] wraz z tabelą siatek kotwicy h ∈ {0.05, 0.025, 0.0125};
- progi QB-2: Φ_c/Φ_vac ∈ {0.197, 0.298, 0.331} [INPUT,
  op-substrate-fluctuation-channel Phase 3];
- g₀_μ = φ·0.90548 = 1.4650974 (definicja INPUT: formuła g₀_μ := φ·g₀_e,
  jak MD op-3d-canonical);
- seed = 20260901, amplituda szumu 1e−3 [LOCK §2];
- β = γ = 1 [INPUT];
- tło 2π: `../op-3d-canonical-lattice-2026-08-31/Phase2_backgrounds3d.npz`
  (klucze `2pi__A1.0__N{32,48}`, `2pi__A0.7__N{32,48}`; READ-ONLY,
  mtime weryfikowany po użyciu).

## 1. Mapowanie g_floor (derywacja z dodatekB — FROZEN, z cytatem)

**Identyfikacja zmiennej.** Pole g akcji kanonicznej cyklu
(S = ∫dt d³x [½K(g)(g_t²−|∇g|²) − U(g)], K = g⁴, U = βg⁷/7 − γg⁸/8)
jest DOKŁADNIE zmienną φ z prop:substrate-action dodatku B:

> (eq:K-geometric) „K_ij = J(φ_i φ_j)², gdzie φ_i = (Φ_i/Φ₀)^{1/2}
> jest lokalną amplitudą pola w węźle i substratu."

> (eq:Lkin-geo, prop:substrate-action) „F_kin^geo[φ] = ∫d³x (K_geo/2)
> φ⁴(∇φ)², K_geo = 2dJ a_sub^{2−d}" — człon kinetyczny wagi φ⁴,
> identyczny z K(g) = g⁴ akcji kanonicznej.

> (eq:EL-alpha2) „∇²φ + 2(∇φ)²/φ = 0" — statyka EL akcji kanonicznej
> cyklu, g″+(2/r)g′+(2/g)g′² = g²(1−g), jest tym operatorem
> z potencjałem U (U′(g) = βg⁶−γg⁷ = K(g)·(βg²−γg³)).

**Normalizacja Φ₀ = Φ_vac.** Z thm:beta-eq-gamma-triple (Droga 1):

> „Definicja Φ₀ jako wartości próżniowej oznacza U′(1)=0, skąd β=γ.
> [...] jest konsekwencją normalizacji (φ_vac = 1)."

Próżnia kanoniczna g ≡ 1 ⇔ φ_vac = 1 ⇔ Φ = Φ₀ = Φ_vac. Zatem próg
gęstościowy QB-2 (spinodala rozrzedzenia w zmiennej Φ) przenosi się na
zmienną polową g przez pierwiastek:

**g_floor := (Φ_c/Φ_vac)^{1/2}** (jedyna wolność mapowania = wybór
Φ₀; zamknięty powyżej literą Drogi 1).

| Φ_c/Φ_vac (QB-2) | g_floor = √(Φ_c/Φ_vac) | g_thr = (1+g_floor)/2 |
|---|---|---|
| 0.197 | 0.4438468 | 0.7219234 |
| 0.298 | 0.5458938 | 0.7729469 |
| 0.331 | 0.5753260 | 0.7876630 |

(wartości = arytmetyka transkrypcyjna; skrypty liczą z formuły).

## 2. Implementacja podłogi (FROZEN: gładka kara C²)

Preferencja LOCKa (kara zachowująca hamiltonowość) przyjęta. Do
potencjału dodajemy karę aktywną WYŁĄCZNIE dla g < g_floor:

**V_fl(g) = (κ/3)·(g_floor − g)³ dla g < g_floor; 0 dla g ≥ g_floor;
κ = 100 (FROZEN, jedna wartość, bez strojenia).**

- V_fl ∈ C²: V_fl′ = −κ(g_floor−g)², V_fl″ = 2κ(g_floor−g) — obie
  znikają w g = g_floor ⟹ struktura wariacyjna/hamiltonowska zachowana
  (U_tot = U + V_fl to zwykły potencjał), gate energii osiągalny.
- Skala κ=100: głębokość penetracji δ z bilansu κδ² = |U′(g_floor)|
  ≈ 0.012 daje δ ≈ 0.011 ≪ g_floor (podłoga „twarda" na skali pola,
  ale gładka numerycznie). Jednorodny stan podłogowy:
  g* = g_floor − δ (LOCKowe „g ≡ g_floor" z dokładnością δ; przy
  klasyfikacji jednorodności używamy kryterium NIESTAŁOŚCI LOCKa
  ‖g−const‖∞ ≥ 0.05, więc δ ≈ 0.011 nie zmienia klasyfikacji).
- ZERO modyfikacji K(g); żadnej projekcji/clampu.

## 3. Schemat gradient flow (FROZEN)

**Równanie:** tłumiony (przetłumiony) gradient flow L²:
∂g/∂t = −δE/δg, E[g] = ∫ [½K(g)|∇g|² + U(g) + V_fl(g)],
tj. ∂g/∂t = ∇·(K∇g) − ½K′(g)|∇g|² − U′(g) − V_fl′(g).
Radialnie: miara r²dr (δE/δg z wagą r²), zerowy strumień na r=0
(automatycznie: r²=0) i r=R (naturalny warunek Neumanna). 3D:
periodycznie. Próżnia g≡1 jest dokładnie stacjonarna (U_tot′(1)=0).

**Dyskretyzacja przestrzenna:** dokładnie ta sama struktura strumieniowa
co RadialDynamicsCan / Dyn3D z op-3d-canonical (gradient dyskretnej
energii: t_flux = K(g_mid)Δg/h, t_quad = ¼K′(g_mid)(Δg/h)²; U1 → U1_tot)
— rhs_gf(g) := −δE_h/δg. K BEZ regularyzacji (w flow nie ma dzielenia
przez K; K ≥ K(g_floor−δ) > 0 w reżimie podłogi).

**Krok czasowy:** stabilizowany semi-implicit Euler:
(I − dt·A_t·L) (g^{n+1} − g^n) = dt·rhs_gf(g^n),
gdzie L = dyskretny laplasjan (3D: symbol FD w przestrzeni Fouriera,
rozwiązanie FFT; radialnie: trójdiagonalny FV z wagą r², solve_banded),
A_t = 1.05·K(max g^n) (parametr stabilizacyjny, adaptowany co krok —
wpływa tylko na stabilność; dokładność czasowa kontrolowana biegiem
dt/2). **dt = 0.01; kontrola nukleacji dt/2 = 0.005.** t_max = 200.

**Kryterium stacjonarności:** ‖rhs_gf(g)‖∞ ≤ 1e−8, próbkowane co Δt=1
(razem z detektorem). Niefinityczność ⟹ załamanie nie-nukleacyjne
(INCONCLUSIVE dla danego biegu). E(t) raportowane deskryptywnie
(monotonia = sanity flow).

**Siatki (FROZEN):** radialnie h ∈ {0.025, 0.0125} (dwie najdrobniejsze
z tabeli kotwicy), R = 60; 3D: pudło L = 2π, N ∈ {32, 48}. Porównania
siatkowe: radialnie — interpolacja liniowa fine→coarse na środkach
komórek coarse; 3D — wspólna podsiatka 16³ (stride 2 z N=32, 3 z N=48;
wzorzec Phase 2 poprzednika).

**Konfiguracje startowe (LOCK §2, geometria FROZEN):**
(i) soliton μ radialny: profil kanoniczny ODE (DOP853, rtol=1e−12,
    atol=1e−14, max_step=0.02 — verbatim op-3d-canonical) z g₀_μ,
    ewaluowany dense-output na obu siatkach radialnych;
(ii) sieć 2π z npz: klucze `2pi__A1.0__N32` i `2pi__A1.0__N48`
    (tag A1.0 = wariant PRIMARY poprzednika; A0.7 to duplikat tego
    samego tła — deskryptywnie odnotowany, nie relaksowany osobno);
(iii) próżnia + szum: 3D, L = 2π, N ∈ {32,48}; szum PASMOWO OGRANICZONY
    wspólny dla obu siatek: rng(20260901) losuje zespolone współczynniki
    C_n (standard_normal (17,17,17,2)) dla modów |n_i| ≤ 8,
    hermityzacja C_sym[n] = (C[n]+conj(C[−n]))/2, pole = Re iFFT
    (te same C_n wbudowane w siatkę N=32 i N=48 ⟹ identyczne pole
    ciągłe na obu siatkach), normalizacja max|f| = 1e−3; g₀ = 1 + f.

Start (i) biegnie w geometrii radialnej (obiekt sferyczny), starty
(ii)/(iii) w 3D (sieć/szum łamią symetrię radialną) — to jest
operacjonalizacja LOCKowego „(radialnie: h z tabeli kotwicy; 3D
N∈{32,48} dla sieci 2π)".

## 4. Detektor nukleacji (implementacja dokładna, FROZEN — LOCK §3)

- Maska: m = (g < g_thr), g_thr = (1+g_floor)/2 (tabela §1).
- Etykietowanie: `scipy.ndimage.label(m)` (3D: struktura domyślna
  6-spójność; sklejanie periodyczne: unifikacja etykiet przez pary
  ścian x/y/z algorytmem union-find; radialnie: label 1D maski w r,
  bez periodyczności).
- Próbkowanie co Δt = 1 (t = 0, 1, 2, …).
- **N_seed := liczba obiektów w t = 0** (zasiane).
- **NUKLEACJA:** pierwsze t₀ takie, że N_obj(t) > N_seed dla WSZYSTKICH
  t ∈ {t₀, t₀+1, …, t₀+10} (wzrost ponad zasiane utrzymany ≥10 jednostek
  czasu). Liczba obiektów w chwili detekcji := N_obj(t₀). Bieg może
  zakończyć się po potwierdzeniu okna.
- **Zbieżność nukleacji (werdykt Q2-PASS-NUCLEATION):** nukleacja obecna
  na OBU siatkach pary ORAZ w biegach dt/2 na obu siatkach; liczba
  obiektów w chwili detekcji zgodna ±1 między wszystkimi czterema
  biegami danego (start × podłoga).
- Detektor NIEZMIENIALNY po pierwszym biegu (forbidden move 2).

## 5. P2a — gate maszynerii podłogi (FROZEN)

- **P2a-i (stacjonarność próżni, osiągalny FAIL):** gradient flow
  z podłogą (każda z 3 wartości), start DOKŁADNIE g ≡ 1, bez zaburzeń,
  t = 10: wymóg ‖g − 1‖∞ ≤ 1e−12 przez cały bieg. Geometrie: radialna
  (h=0.0125) i 3D (N=32).
- **P2a-ii (gate energii tam, gdzie podłoga nieaktywna):** ewolucja
  HAMILTONOWSKA (RK4, K_ε=0.2 — wzorzec P1d poprzednika) z U_tot,
  start g = 1 + a·cos(x) (3D, L=2π, N=32, a=1e−3; radialnie: g = 1 +
  a·exp(−(r−10)²), a=1e−3), dt=0.004, t_end=4: wymóg |ΔE|/|E| ≤ 1e−6
  ORAZ min g > g_floor przez cały bieg (podłoga nieaktywna —
  weryfikowane). FAIL któregokolwiek ⟹ STOP (litera LOCKa).

## 6. P2c — kontrola sektora stabilnego (FROZEN)

Identyczna procedura (te same 3 starty × 3 podłogi × 2 siatki, ten sam
schemat, ten sam detektor, te same progi) z potencjałem sektora
stabilnego **U_stab(g) = −βg⁷/7 + γg⁸/8** (linearyzacja wokół g=1:
m² = +γ; forma zgodna z sektorem stabilnym P3c op-3d-canonical,
u-forma V = −(4βg−5γg²), V(1) = +1). Podłoga i detektor bez zmian.
Wymóg: wszystkie biegi → próżnia g≡1 (stacjonarnie), ZERO nukleacji.
Jakikolwiek alarm detektora ⟹ P2c FAIL ⟹ STOP (detektor nieważny).

## 7. Phase 1 — doprecyzowania (FROZEN)

- **P1a (sympy):** U(0)=0, U(1)=1/56, U″(1)=−1 exact; brak minimum
  lokalnego w (0,1): solveset U′=0 w (0,1) = ∅ (U′ = g⁶(β−γg) > 0
  na (0,1)) — formalnie, z wnioskiem o kierunku relaksacji.
- **P1b (kwadratura):** profil radialny μ jak §3(i); E_rel[g|ref] =
  ∫₀^R 4πr² [½K(g)g′² + U(g) − U(g_ref)] dr (midpoint na siatkach
  cell-centered, wszystkie trzy h z tabeli kotwicy — zbieżność
  raportowana); referencje: próżnia (U(1)=1/56) i stan pusty (U(0)=0).
- **P1c:** ε(2π) = (E_cell[g] − E_cell[1])/d³, E_cell = Σ[½K(g_mid)
  (Δg/h)² + U(g)]h³ (energia dyskretna Dyn3D, pi=0, K bez
  regularyzacji), dla wszystkich 4 tablic npz (A1.0/A0.7 × N32/N48),
  BEZ relaksowania.
- **Werdykt Q1 wg litery:** zbiór policzony = {ΔE_create(sol|vac),
  ΔE_create(sol|empty), ε(2π)}; granica istnieje ⟺ w zbiorze jest
  konfiguracja z μ=0 oddzielająca reżim μ<0 od μ>0 (zmiana znaku
  w policzonym zbiorze).

## 8. Phase 3 — operator z podłogą (FROZEN; tylko przy Q2-PASS-STATIC)

Druga wariacja E wokół stanu zrelaksowanego g₀ (spełniającego EL
z U_tot) daje operator wagi-K: −∇·(K∇φ) + Q_tot φ = ω² K φ,
**Q_tot = K(2βg−3γg² + 2|∇g|²/g²) + [V_fl″(g) − (4/g)·V_fl′(g)]**
(człony podłogowe z tej samej derywacji, która dała Q poprzednika:
Q = U_tot″ − (K′/K)U_tot′ + (K′²/2K − K″/2)|∇g|²). Maszyneria:
build_op3d + eigsh(tol=0, v0 det. rng(20260831)) + translation_analysis
verbatim z `Phase3_bloch3d.py` op-3d-canonical; |∇g|² centralne różnice.
- Tło 3D: Bloch Γ/X/M/R (fazy jak kpoints(d=2π)), N ∈ {32,48}.
- Tło radialne: forma FV wagi-K jak wK_lam_radial z Q→Q_tot,
  h ∈ {0.025, 0.0125}; mody translacyjne nieobecne w sektorze s
  (ℓ=1) — odnotowane; mody podłogowe identyfikowane przez lokalizację
  w regionie aktywnej kary (udział ∫_{g<g_floor}φ²/∫φ² ≥ 0.5).
- Zbieżność i tol wg LOCKa §3 Phase 3 (≤0.05·max(|ω²_min|,0.1);
  Q3-PASS: ω²_min ≥ −1e−3 zbieżnie, po odjęciu modów zerowych
  zidentyfikowanych PRZED interpretacją).
- Przy Q2-PASS-NUCLEATION zamiast spektrum: charakteryzacja kaskady
  BEZ progów — N_obj(t) (tempo kreacji = przyrosty na jednostkę czasu)
  i rozkład rozmiarów obiektów (liczności komórek per etykieta)
  w chwilach t₀, t₀+5, t₀+10.

## 9. Higiena wykonania

Pełne ścieżki bez `cd`; `ls` po każdym zapisie; runy >10 min w tle,
proces tła ≤ ~55 min (podział na etapy: `Phase2_floor_relax.py --stage k`
z wynikami per-run w json/npz; `--verdict` składa całość); zalockowane
siatki NIE zmniejszane; INCOMPLETE z przyczyną per bieg. Outputy:
`Phase1_output.txt`, `Phase2_output.txt` (+ `Phase2_relaxed_states.npz`,
per-stage `Phase2_stage*.txt` sklejane), warunkowo `Phase3_output.txt`.
Rdzeń `.tex`, STATE.md, git, katalogi innych cykli — NIETYKANE
(npz READ-ONLY, mtime weryfikowany po odczycie).

**FROZEN 2026-09-01, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
