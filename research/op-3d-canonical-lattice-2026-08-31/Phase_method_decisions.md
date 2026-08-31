---
title: "Phase_method_decisions — delta względem op-3d-lattice-bath-stability: cały cykl w akcji kanonicznej (profil radialny kanoniczny, kotwica w u-formie w1, most na operatorze wagi-K Phase 3, detektor t* kanoniczny, korekta eigsh tol=0 wszędzie); decyzje ZAMROŻONE przed startem obliczeń"
date: 2026-08-31
type: method-decisions
tgp_owner: research/op-3d-canonical-lattice-2026-08-31
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/Phase_method_decisions.md]]"
  - "[[../op-3d-lattice-bath-stability-2026-08-31/Phase_correction_note_eigsh.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_method_decisions.md]]"
---

# Decyzje metodyczne (delta; zamrożone PRZED jakimkolwiek obliczeniem cyklu)

Baza: `../op-3d-lattice-bath-stability-2026-08-31/Phase_method_decisions.md`
(dalej: MD-poprzednika) — **przejęta w całości tam, gdzie niżej nie ma
delty** (operator 3D §1, Bloch/eigsh §2, interpolacja §3, newton_krylov §4,
Phase 3 §6, Phase 4 §7, higiena §8). Kryteria/progi/punkty LOCKa
(Phase0_balance.md) NIETKNIĘTE. Jedyne obliczenie wykonane przed
zamrożeniem: arytmetyka rejestru wejść φ·0.90548 = 1.4650974 (weryfikacja
transkrypcji, nie rachunek cyklu).

**Rejestr WEJŚĆ (flagowany w każdym zależnym wyniku):** g₀_e=0.90548
[INPUT r₂₁/φ-FP, bath-two-sectors], φ=(1+√5)/2 (kalibracja μ),
g₀_μ = φ·g₀_e = 1.4650974, d*₁(μμ, M-P)=3.0790 [INPUT bath-two-sectors
Phase 1], ε=0.2 (kontrola 0.1) [INPUT], β=γ=1.

**Doprecyzowanie transkrypcyjne (PRZED obliczeniami):** LOCK §1/§2 podaje
„g₀=1.46507=φ·0.90548"; źródło INPUT (bath-two-sectors Phase1_output.txt:
„g0_mu=phi*g0_e=1.46510") i arytmetyka dają φ·0.90548=1.4650974. Wartość
1.46507 w LOCKu to błąd zaokrąglenia transkrypcji (odchył 3·10⁻⁵, bez
wpływu na jakąkolwiek bramkę ±5%/±15%). ZAMROŻONE: **definicją INPUT jest
formuła g₀_μ := φ·g₀_e** (litera źródła M-P), implementacja
`G0_MU = PHI*0.90548`.

## 0. Brak rozdziału modeli (delta zasadnicza vs MD-poprzednika §0)

CAŁY cykl w akcji kanonicznej: S = ∫dt d³x [½K(g)(g_t²−|∇g|²) − U(g)],
K=g⁴, U=βg⁷/7−γg⁸/8, β=γ=1. Statyka EL: ∇²g + (2/g)|∇g|² = g²(1−g)
(radialnie: g″+(2/r)g′+(2/g)g′² = g²(1−g)). Spektra BEZ regularyzacji;
ewolucje nieliniowe z K_ε (mapa x↦½(x+√(x²+ε²)) na K; ε=0.2, kontrola
0.1). Kotwice #63 (−1.3896/3.62) NIE są bramkami tego cyklu
(autoryzacja user-gate w LOCKu).

**Dziedziczenie P1a/P1d (LOCK §1):** legalne — model kanoniczny i
integratory VERBATIM z `Phase1_gate3d.py` poprzednika (te same funkcje
K, U, build_op3d, Dyn3D, RK4), środowisko identyczne (CPython 3.14.2,
numpy 2.4.3, scipy 1.17.1 — sprawdzone importem przed startem).
Odnotowywane w output z cytatem liczb, NIE liczone ponownie.

## 1. P1b-kan — kotwica radialna kanoniczna (doprecyzowania FROZEN)

- **Profil radialny μ:** DOP853 od r=10⁻⁶, y0=[g₀_μ, 0], rhs kanoniczny:
  g″ = g²(1−g) − (2/g)g′² − (2/r)g′; w r<10⁻¹⁰: g″ = [g²(1−g)−(2/g)g′²]/3;
  max_step=0.02, rtol=1e−12, atol=1e−14, dense_output. Siatka master:
  cell-centered r_j=(j+½)·0.0125, R=60 (N=4800). Siatki kotwicy
  h∈{0.05,0.025,0.0125}: ewaluacja dense output BEZPOŚREDNIO na każdej
  siatce (bez reinterpolacji).
- **Operacjonalizacja „residuum ≤1e−10" (LOCK):** ‖g(rtol=1e−12) −
  g(rtol=1e−13)‖∞ na siatce master ≤ 1e−10 (dokładność rozwiązania ODE;
  FAIL osiągalny). Deskryptywnie obok: residuum EL z g″ policzonego
  niezależnie splajnem z tablicy g′ (norma ∞ raportowana, informacyjnie —
  nosi błąd splajnu O(h³), nie jest bramką).
- **Kotwica λ_min(w1) — forma FROZEN: u-forma kanoniczna.** Podstawienie
  u=g³/3 (χ=g₀²φ) daje ściśle problem wagi 1:
  −∇²χ + V(g₀)χ = λχ, **V = 4βg₀ − 5γg₀²** (próżnia: 4−5=−1 ✓ —
  identycznie jak cross-check MD-poprzednika §1 i P3a tego cyklu).
  Mod s-falowy: u_r = r·χ, −u_r″ + V(r)u_r = λu_r, siatka cell-centered
  r_j=(j+½)h, Dirichlet (u_r=0 poza siatką), eigh_tridiagonal.
  **Kotwica := λ_min przy h=0.0125**; gate wewnętrzny LOCKa:
  |λ(0.025)−λ(0.0125)| ≤ 1e−3·|λ(0.0125)|. Deskryptywnie obok (nie
  bramka): λ_min radialnej formy wagi-K (FV jak build_tridiag_gen
  poprzednika z F→K, Q→K(2βg−3γg²+2g′²/g²), waga r²K) na tych samych h —
  kontrola równoważności form (rozjazd ≫O(h²) = czerwona flaga).
- **Diagnostyka kieszeni Q(r) (OBOWIĄZKOWA, LOCK):** profil pointwise
  V(r)=4βg−5γg² (u-forma; to jest potencjał linearyzacji wagi 1):
  min, r_min, głębokość D=V_min−V_vac (V_vac=−1), FWHM := szerokość
  przedziału {r: V(r) ≤ V_vac + D/2}; wniosek anty-pułapkowy: FWHM vs
  h=0.3 mostu. Deskryptywnie także Q/K formy ważonej:
  2βg−3γg²+2g′²/g².
- **t*_ref (ewolucja radialna kanoniczna):** siatka master (h=0.0125,
  R=60); dynamika = RadialDynamics poprzednika z (F_ε,W)→(K_ε,U)
  (struktura hamiltonowska verbatim, RK4); start g = g_eq + a·φ₁,
  φ₁ = (u₁/r)/g_eq² z modu kotwicy (h=0.0125), normalizacja
  ∫r²φ₁²dr·h=1, znak: wartość w argmax|φ₁| dodatnia; a=±0.01.
  **Detektor t* (FROZEN; delta vs #63 — runaway kanoniczny może iść
  w g→duże, U nieograniczone z dołu):** t* := pierwszy krok z
  (min g ≤ 0.01) LUB (max g ≥ 2·g₀_μ) LUB niefinityczność — co krok.
  Biegi: K_ε=0.2: (a=+0.01, dt=0.004), (a=+0.01, dt=0.002),
  (a=−0.01, dt=0.004); kontrola K_ε=0.1: (a=+0.01, dt=0.004).
  **dt-konsystencja (LOCK ≤1%):** |t*(0.004)−t*(0.002)|/t*(0.002) ≤ 0.01
  dla pary (K_ε=0.2, a=+). t*_ref := min po biegach K_ε=0.2 (konwencja
  min-po-biegach poprzednika); kontrola 0.1 raportowana deskryptywnie.
  t_end=20 (brak detekcji ⟹ t*=BRAK, raport wprost). Gate energii
  |ΔE|/E raportowany deskryptywnie (integrator certyfikowany w P1d
  dziedziczonym; brak progu w LOCKu dla tego biegu).

## 2. P1c-kan — most (doprecyzowania FROZEN)

- **Operator mostu = operator Phase 3 (waga K):** −∇·(K∇φ)+Qφ = λKφ,
  Q = K(2βg−3γg²+2|∇g|²/g²), build_op3d verbatim (7-punktowy,
  K na linkach ze średniej węzłowej, symetryzacja B^{−½}AB^{−½}),
  pudło L=30 periodyczne (fazy +1), soliton w centrum; |∇g|² = g′(r)²
  i ∇²g = g″+2g′/r ANALITYCZNIE ze splajnów radialnych (CubicSpline
  osobno dla g−1, g′, g″ na siatce master; r>r_end: próżnia; r→0:
  ∇²g=3g″(0⁺)) — MD-poprzednika §3 verbatim. Siatki: N=76 (h=0.3947,
  „h≈0.4" LOCKa jak u poprzednika) i N=100 (h=0.30).
  **Gate (litera LOCKa):** |λ_min(N=100) − kotwica| ≤ 0.05·|kotwica|
  ORAZ |λ_min(N=76) − kotwica| > |λ_min(N=100) − kotwica|.
  Deskryptywnie obok: λ_min 3D u-formy (waga 1: −∇²+(4βg−5γg²)) na obu
  siatkach — kontrola równoważności form w 3D.
- **t*_izo(3D):** siatka N=100 (h=0.3 — siatka bramki), dynamika
  Dyn3D(K_ε=0.2, U) verbatim; mod v: najniższy wektor eigsh operatora
  wagi-K przy N=100, φ=ψ/√K, normalizacja ‖φ‖_{L²(h³)}=1, znak:
  wartość w argmax|φ| dodatnia; start g₀3D + a·φ. Biegi:
  (a=+0.01, dt=0.02), (a=+0.01, dt=0.01), (a=−0.01, dt=0.02)
  (konwencja poprzednika; CFL: dt ≪ h=0.3). Detektor t* jak w §1.
  t_end = max(12, 2·t*_ref). **Gate (LOCK): KAŻDY bieg
  t* ∈ t*_ref·[0.85, 1.15]; t*_izo(3D) := min po biegach** (zasila
  Phase 4: t_end = 3·t*_izo(3D)).
- **eigsh wszędzie w cyklu: tol=0** (korekta 1 poprzednika,
  `Phase_correction_note_eigsh.md`, przyjęta od startu — dotyczy także
  Phase 1), k≥10, ncv=80, maxiter=200000, v0 deterministyczny
  pseudolosowy (rng(20260831), znormalizowany). Brak zbieżności ARPACK
  → punkt INCOMPLETE.

## 3. Phase 2 (delta względem skryptu poprzednika)

`Phase2_relax_lattice3d.py` SKOPIOWANY z
`../op-3d-lattice-bath-stability-2026-08-31/` (newton_krylov, precond
spilu, starty, klasyfikacja, dedup, fft_prolong, guard 0.05 — verbatim).
Delty (wyłącznie):
1. profil startowy = radialny profil KANONICZNY μ (ODE z §1; INPUT
   g₀_μ=φ·0.90548), nie profil #63 F_ε — obiekt PRIMARY LOCKa §2;
2. ścieżki BASE/OUT na katalog tego cyklu;
3. nagłówek/opisy (model kanoniczny od startu — bez zmiany równań:
   residuum EL było już kanoniczne u poprzednika).
Kryteria LOCKa bez zmian: ‖R‖∞≤1e−8, ‖g−1‖∞≥0.05, zbieżność 5e−3
(wspólna podsiatka 16³), d∈{π,2π,3.0790,3π,4π}×{1.0,0.7}×N∈{32,48}.

## 4. Phase 3 (delta względem skryptu poprzednika)

`Phase3_bloch3d.py` SKOPIOWANY verbatim (zawiera już korektę eigsh
tol=0 + kontrolę pokrycia translacji + P3b/P3c + cross-check u-formy).
Delty (wyłącznie): ścieżki BASE/NPZ/OUT_JSON na katalog tego cyklu;
nagłówek (INPUT g₀_μ). Kryteria/progi/punkty verbatim MD-poprzednika §6.

## 5. Phase 4 (warunkowa; doprecyzowania MD-poprzednika §7 przejęte)

Superkomórka 2×2×2 (N_sc=96), mod minimalny NIE-translacyjny z Phase 3
(N=48, argmin k), rozszerzenie φ_sc fazami ±1, ‖φ_sc‖∞=1;
a = 0.01·max|g₀−1|; ewolucja kanoniczna K_ε, RK4, dt=0.01, kontrole:
dt/2 przy (+, ε=0.2), ε=0.1 przy (+), oba znaki przy ε=0.2;
t_end = 3·t*_izo(3D). **Ucieczka (delta detektora — spójna z §1):**
(g ≤ 0.01 LUB max g ≥ 2·max(g₀) LUB niefinityczne — co krok) LUB
max|g−g₀| > 1.0·‖tło‖ (co próbkę 0.02). Gate energii ≤1e−6 do momentu
ucieczki. Ruling kwantyfikatora d przy istnieniu częściowym: zapisany w
`Phase3_verdict_ruling.md` PO Phase 3, PRZED Phase 4 (PRIMARY po tłach
istniejących, strict obok — wzorzec bloch-chain).

## 6. Higiena wykonania

Jak MD-poprzednika §8: pełne ścieżki, zero `cd`, `ls` po każdym zapisie;
runy >10 min w tle; zalockowane siatki NIE zmniejszane; INCOMPLETE
z przyczyną dozwolone per punkt. Outputy: `Phase1_output.txt`,
`Phase2_output.txt` + `Phase2_backgrounds3d.npz`, `Phase3_output.txt`
+ `Phase3_results3d.json`, warunkowo `Phase4_output.txt`. Rdzeń `.tex`,
STATE.md, git, katalogi innych cykli — NIETYKANE (odczyt wzorców
dozwolony).

**FROZEN 2026-08-31, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
