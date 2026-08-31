---
title: "Phase_method_decisions — decyzje implementacyjne spisane PRZED startem obliczeń (LOCK §2 wymaga tego dla wariantu wymiarowego Q2; reszta = higiena anti-Lakatos)"
date: 2026-08-23
type: method-note
tgp_owner: research/op-bath-two-sectors-2026-08-23
status: FROZEN-BEFORE-COMPUTATION
computations_performed_at_freeze: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[HANDOFF_PROMPT.md]]"
---

# Decyzje implementacyjne (zamrożone PRZED pierwszym uruchomieniem skryptu)

Wszystko poniżej to USZCZEGÓŁOWIENIE LOCKa (nie zmiana kryteriów).
Po starcie obliczeń zmiany tylko dla udokumentowanego błędu implementacji
(LOCK §4.1), opisane przed użyciem wyniku.

## Phase 1

1. **Solvery:** M-P = reuse `solve_m2` z `N1_provenance_check.py`
   (g″+(2/r)g′ = g²(1−g) − (2/g)g′², start r=0.01, rtol=1e−11, atol=1e−13,
   max_step=0.05); M-L = reuse `solve_p131` (F g″+(2/r)g′=g²(1−g),
   F=1+2α_eff ln g, α_eff=2/(1+η_K(g−1)²), η_K=181/15 **[INPUT]**).
2. **Okno PRIMARY** dla (A, δ) idących do Phase 2: **[120,260]**
   (bardziej asymptotyczne; konwencja N1). Okno [50,150] = kontrola
   rozjazdu (>3° ⇒ WINDOW-SENSITIVE, wg LOCK).
3. **R-kontrola:** „|Δδ| ≤ 1% zakresu" interpretuję jako 1% pełnego
   zakresu fazy 360°, tj. **|Δδ| ≤ 3.6°** (R=300 vs R=450).
4. **R² fitu:** R² = 1 − SS_res/SS_tot na (g−1)·r w oknie; gate ≥ 0.999.
5. **Tłumienie κ:** κ_eff = ln(A_[50,150]/A_[120,260]) / (r̄₂−r̄₁),
   r̄ = środki okien (Δr̄=90). Jeśli |κ_eff| > 1e−3 → używam zmierzonego
   κ w E_int i w kąpieli; inaczej κ:=0 (LOCK dopuszcza κ≈0).
6. **Konwencja E_int pary (i,j):**
   `E_int(d) = −A_i A_j · exp(−(κ_i+κ_j)d/2) · cos(d − δ_i − δ_j)/d`
   — suma faz z nakładki dwóch zespolonych amplitud ogona
   (A e^{iδ}); odstęp drabiny (2π) od konwencji NIE zależy, offset d* tak.
   Zakres poszukiwań minimów: d ∈ (2, 60), krok siatki 1e−4·(60−2).
7. **P1c (Yukawa):** `E(d) = −A_i A_j e^{−d}/d` (te same A, m=1);
   wymagane 0 minimów w (2,60).
8. **τ w M-P:** g₀_τ=1.5696 **[INPUT: Q_K=3/2]**; kontrola wrażliwości:
   g₀_τ·(1±0.02) = {1.5382, 1.6010}; 1.6010 > 8/5 — jeśli skolapsuje,
   raportuję kolaps jako wynik wrażliwości (τ leży 1.9% pod progiem).
9. Kalibracja g₀_e=0.90548, g₀_μ=φ·g₀_e **[INPUT: r₂₁/φ-FP]**.

## Phase 2 (Q1)

1. **Model dynamiczny:** M0-f_eps, ε=0.2 (jak #63 V3), kod reuse
   z `Phase3_nonlinear_evolution.py` (SoftWallDynamics, RK4, gate 1e−6).
2. **P2a baseline:** R=60, N=4000 (siatka #63), a=±0.01, dt=0.004
   i 0.002; oczekiwana reprodukcja BREAKDOWN t*≈3.62 (out #63);
   t*_izolacji := min t* z przebiegów baseline (konserwatywnie).
3. **Obiekt centralny per konfiguracja:** ee → soliton e (g₀=1.24915
   w f_eps, gałąź a3d **[INPUT kalibracji]**), μμ i eμ → soliton μ
   (g₀=φ·1.24915; jedyny audytowany reprezentant runaway #63).
   Jeśli profil e w f_eps nie dobiega R — konfiguracja ee raportowana
   jako NIEREPREZENTOWALNA w modelu dynamicznym (deskryptywnie),
   bez zamiany obiektu.
4. **Kąpiel (M-B):** średnia kątowa ogona sąsiada po sferze sąsiadów
   w odległości d (tożsamość Helmholtza dla ω=1, r<d):
   `δg_bath(r) = N_c · A_s · e^{−κ_s d} · cos(d−δ_s)/d · sin(r)/r`,
   gatunek sąsiada s: ee→e, μμ→μ, eμ→e (centralny μ);
   **N_c = 12** (koordynacja FCC/WS; decyzja implementacyjna).
   (A, δ, κ) WYŁĄCZNIE z Phase 1, model **M-P PRIMARY** (LOCK §2);
   przekaz przez `Phase1_results.json` (zero ręcznego przepisywania).
5. **Komórka:** promień komórki WS **R_cell ∈ {d/2, d}** (2 rozmiary);
   siatki **h ∈ {0.015, 0.0075}** (2 siatki); zero-flux BC (jak #63).
   Tło w komórce = obcięty profil izolowany + δg_bath (konfiguracja
   kąpieli; LOCK mówi „spektrum linearyzacji wokół KONFIGURACJI kąpieli"
   — to nie jest równowaga statyczna; raportowane wprost).
6. **ω²:** problem F-ważony (dokładna dynamika, #63 method note):
   L v = λ_gen F v, **ω²_min := λ_gen,min**; waga-1 λ_min raportowana
   obok. Zbieżność „do 2 cyfr" = max rozrzut ω²_min po {siatki×komórki}
   ≤ 0.01·max(|ω²_min|, 0.1).
7. **P2c:** spektrum próżni g≡1 w TEJ SAMEJ komórce (każdy rozmiar);
   mody komórki (przewidywanie ω²_n=((nπ/R)²−1)/F(1)) identyfikowane
   przed interpretacją; N_neg vs floor(R/π).
8. **Ewolucja:** perturbacja ± wzdłuż modu runaway = wektor własny
   λ_min spektrum KOMÓRKI (mod zlokalizowany; jeśli λ_min>0 — wektor
   własny najniższy), amp ±0.01; dodatkowo przebieg amp=0 (kontrola
   sloshingu niestacjonarnego tła — dodana, nieusuwana); t_end = 3·t*_izo;
   gate energii raportowany dla wszystkich przebiegów.
9. Reprodukcja izolacji w każdej komórce (amp bath = 0) jako komparator
   (dodatek dozwolony; punkty zalockowane nieusuwane).

## Phase 3 (Q2)

1. **WARIANT WYMIAROWY (wymóg LOCK §2, zapis przed startem):**
   **komórka 1D periodyczna** o okresie d (gęstość n ~ 1/d);
   uzasadnienie: dokładna sieć periodyczna źródeł, samouzgodnione tło
   przez Newtona bez artefaktów centrum sferycznego. NIE zmieniam po starcie.
2. **Źródła:** ρ̃(x) = suma obrazów gaussowskich, **σ_s = 0.5**,
   normalizacja ∫_komórka ρ̃ dx = 1; **q = 1.0 (PRIMARY)**, kontrola
   wrażliwości q = 0.3 przy d=4 (zapisana teraz).
3. **Równanie (M-S):** ψ″ + 2ψ′²/ψ + βψ² − γψ³ = −qρ̃, β=γ=1;
   Newton z tłumieniem + kontynuacja w q; start ψ≡1.
4. **d=∞ (kontrola P3a):** domena [−40,40], N=4000, Dirichlet ψ=1;
   gate kontroli: zanik monotoniczny + fit m z ln(ψ−1) w [8,20]:
   **m² = γ ± 5%** oraz ω²_min ≈ +γ (znak +). FAIL kontroli ⇒ STOP Q2
   (procedura nieważna) — osiągalny FAIL.
5. **Operator fluktuacji:** linearyzacja residuum:
   L̂φ = φ″ + (4ψ′/ψ)φ′ + [2ψ−3ψ² − 2ψ′²/ψ²]φ; forma samosprzężona
   z wagą **w = ψ⁴** (bo w′/w = 4ψ′/ψ):
   −(wφ′)′ − w·Q̃·φ = ω²·w·φ, Q̃ = 2ψ−3ψ²−2ψ′²/ψ².
   ω²_min = najniższa wartość własna uogólnionego problemu symetrycznego.
   Wokół ψ=1: ω²_min = γ + k² ≥ γ ✓ (kotwica konwencji).
6. **Siatki:** N/okres ∈ {400, 800} (2 siatki); d ∈ {∞, 8, 6, 4, 2}
   (LOCK). Zbieżność jak w Phase 2 pkt 6.
7. **Odpowiedź statyczna:** superkomórka 10 okresów, tło = tiling ψ_n;
   rozwiązanie L̂χ = −δ_h (dyskretna delta w centrum superkomórki);
   klasyfikacja: liczba zmian znaku χ na dystansie >2σ_s od źródła
   (0 → monotoniczna/Yukawa; ≥2 → oscylacyjna).

## Rejestr WEJŚĆ (LOCK §2)

Q_K=3/2 (g₀_τ), r₂₁/φ-FP (g₀_e, g₀_μ), η_K=181/15, g₀=1.24915 (gałąź
f_eps/a3d), ε=0.2, N_c=12, σ_s=0.5, q=1.0 — flagowane INPUT w outputach.

**Zamrożono 2026-08-23, przed pierwszym uruchomieniem jakiegokolwiek
skryptu tego cyklu.**
