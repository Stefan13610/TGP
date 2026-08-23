---
title: "Phase0_balance — LOCK cyklu op-wall-dynamics (interpretacja ściany + stabilność z więzem budżetu przestrzeni)"
date: 2026-07-03
type: phase0-lock
tgp_owner: research/op-wall-dynamics-2026-07-03
status: "LOCKED — realizacja: NASTĘPNA SESJA (osobny agent-implementator; Phase 1–3 user-authorized 2026-07-03: 'zajmiemy się tym w następnej sesji')"
anti_lakatos_lock: PRESERVED
tags: [wall-dynamics, constrained-stability, ghost-wall, Q-ball, L03, L04, korona]
related:
  - "[[../op-spectral-analysis-Phi-2026-07-03/README.md]]"
  - "[[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]]"
  - "[[../../audyt/L03_K_phi_stability/POST_ACTION_UPDATE_2026-07-03.md]]"
---

# Phase 0 — LOCK: op-wall-dynamics

## 0. Kontekst dla agenta (handoff — przeczytaj PRZED czymkolwiek)

**Skąd ten cykl.** CP-7 (`op-spectral-analysis-Phi-2026-07-03`, CLOSED)
wykonał diagonalizację operatora fluktuacji. Wyniki krytyczne:

1. Sektor grawitacyjny (K=ψ⁴+U_A): czysty, σ≥0 — NIE dotykać.
2. Sektor solitonowy (funkcjonał E_S = ∫½f(g)|∇g|²+W(g),
   f(g)=1+2α·ln g, α=2, W′(g)=g²(1−g); jego EL = ODE korony a3d/ls10):
   próżnia g=1 tachioniczna (kontinuum od −1), mody ZLOKALIZOWANE l=0:
   **e: 0, μ: 2 (−1,282; −1,057), τ: 3 (−4,216; −1,114; −1,010)**
   (jedn. γ=β=K_geo=1). μ/τ = siodła funkcjonału BEZ więzów.
3. Ghost wall g* = e^{−1/(2α)} ≈ 0,7788: profile μ/τ DOTYKAJĄ ściany
   (min g = 0,788/0,786; 1/3 odbicia; min f(g) ≈ 0,04); odbicie
   g′→−g′ przy g*+0,005 = regularyzacja ad-hoc (nie-EL). W substracie
   α=1 (bez ściany) profil τ KOLABUJE (min g=0,158, urywa się r≈3).
4. Narzędzia gotowe do reuse:
   `../op-spectral-analysis-Phi-2026-07-03/Phase2_bvp_spectrum.py`
   (funkcje: `soliton_profile`, `spectrum_on_background`,
   `assemble_and_solve`, `newton_background`; wszystkie zweryfikowane
   zbieżnością; stałe G0_E=1.24915, G0_MU=φ·G0_E, G0_TAU=3.18912).

**Hipoteza robocza cyklu (autor, 2026-07-03, zapisana PRZED analizą):**
> Ściana wynika z wewnętrznej energii solitonu — z ilości tworzonej
> przestrzeni, która w rdzeniu obiektu przekracza próg stabilności.

Formalizacja: „ilość tworzonej przestrzeni" B[g] = ∫(g−1)·w(r)·4πr²dr
jest wielkością budżetowaną ⇒ fizycznie dopuszczalne fluktuacje
zachowują budżet (δB=0) ⇒ stabilność = stabilność WARUNKOWA na
podprzestrzeni więzu (analogia Q-ball: bez więzu siodło, z więzem
minimum). Ściana dolna g* i górny ogranicznik rdzenia g_crit=8/5
(H7/H8 w `tgp_master_consistency_v47.py`) mogą być dwoma progami
jednego mechanizmu budżetowego.

## 1. Cel cyklu

Rozstrzygnąć, czy mody ujemne μ/τ z CP-7 znikają po nałożeniu
fizycznego więzu/dynamiki ściany — tj. czy stabilność korony jest
twierdzeniem (constrained), a nie założeniem.

## 2. Fazy i zadania (LOCKED)

### Phase 1 — stabilność z więzem budżetu przestrzeni (główna)

Operator L̂ z CP-7 na tłach e/μ/τ (F-S, α=2, identyczne profile —
reuse `soliton_profile`). Trzy kandydackie więzy (wszystkie policzyć,
żadnego nie wybierać post-hoc):
- **K1 (budżet prosty):** ⟨v⟩_w = ∫v·r²dr = 0 (w=1).
- **K2 (budżet metryczny):** ∫v·g_eq²(r)·r²dr = 0 (waga objętościowa
  √−g_eff ~ g²; uzasadnienie: przestrzeń tworzona ~ Φ, element
  objętości zależy od Φ).
- **K3 (budżet kinetyczny):** ∫v·f(g_eq(r))·r²dr = 0 (waga = miara
  S-L operatora).
Metoda: projekcja P L̂ P na podprzestrzeń ortogonalną do wektora więzu
(deflacja w dyskretyzacji tridiagonalnej z CP-7; uwaga na wagę B=r²
przy ortogonalności) LUB minimalizacja Rayleigha z mnożnikiem Lagrange'a.
Wynik: N_neg^constrained(tło, l=0, więz) + zbieżność N∈{2k,4k,8k}.

**Kryterium W1 (LOCKED):**
- N_neg^constrained(μ)=0 ORAZ N_neg^constrained(τ)≤1 przy dowolnym z K1–K3
  → „stabilizacja więzem POTWIERDZONA (dla tego więzu)"
  (τ: 1 rezydualny mod dozwolony TYLKO jeśli odpowiada kierunkowi
  wzdłuż rodziny profili g₀ — test: overlap z ∂g/∂g₀ > 0,9).
- Mody przetrwają wszystkie K1–K3 → **wynik NEGATYWNY zgłaszany wprost**
  (stabilizacja budżetem obalona w wersji liniowej; hipoteza autora
  wymaga wtedy więzu nieliniowego/innego ładunku — dokumentować).
- UWAGA metodologiczna: jeden więz skalarny usuwa co najwyżej 1 mod
  ujemny na mocy teorii przeplatania — μ ma 2, τ ma 3 mody. Jeśli po
  K1–K3 (1 więz) zostaje >0 modów, sprawdzić kombinacje (K_i ∧ K_j)
  oraz więz lokalnego budżetu w kuli rdzenia (K4: ∫_{r<r_core} v r²dr=0,
  r_core = pierwszy węzeł g_eq−1). To NIE jest zmiana kryterium —
  kombinacje zadeklarowane tu, w LOCK-u.

### Phase 2 — ściana jako warunek jednostronny + regularyzacje

- **W2a:** spektrum z warunkiem jednostronnym g_eq+v ≥ g* na zbiorze
  kontaktu (tam gdzie g_eq < g*+0,01): przybliżenie Dirichletem v=0 na
  zbiorze kontaktu + porównanie z pełnym problemem komplementarności
  (LCP, jeśli wykonalne; inaczej odnotować ograniczenie).
- **W2b (soft wall):** rodzina f_ε(g)=½[f(g)+√(f(g)²+ε²)] (gładka,
  >0); dla ε ∈ {0.2, 0.1, 0.05, 0.02}: re-solve profili μ/τ (ODE z f_ε)
  i spektrum. Pytania zalockowane: (i) czy λ_min(ε) zbiega przy ε→0,
  czy dywerguje (regularization-dependence z CP-7 — kwantyfikacja);
  (ii) czy r₂₁, r₃₁ (A_tail/Koide) są odporne na ε (dopuszczalny dryf
  <0,1% — inaczej mechanizm korony jest wrażliwy na model ściany:
  raportować wprost).

**Kryterium W2 (LOCKED):** rozstrzygnięcie „spektrum μ/τ jest / nie jest
zbieżne przy ε→0" + tabela dryfu r₂₁/r₃₁(ε). Dowolny kierunek = PASS
wykonania; brak zbieżności = wynik negatywny do Limitations.

### Phase 3 — budżet jako wspólne źródło obu progów (analityczna)

- **W3a:** zdefiniować E_core[g₀] i B_core[g₀] (energia i budżet
  przestrzeni w rdzeniu profilu) z profili ODE; sprawdzić, czy g*=e^{−1/4}
  (dolny) i g_crit=8/5 (górny, H7) odpowiadają ekstremom/utracie
  monotoniczności tej samej wielkości budżetowej (sympy + numeryka).
- **W3b (SPECULATIVE, zero claimów):** indeks siodłowy 0/2/3 vs
  generacja/odbicia — tylko dokumentacja korelacji, żadnych promocji.

**Kryterium W3 (LOCKED):** W3a rozstrzygnięte (tak/nie/nierozstrzygalne
z powodu X); W3b wyłącznie deskryptywne.

## 3. Forbidden moves (anti-Lakatos)

1. Zmiana kryteriów/tolerancji/listy więzów PO uruchomieniu obliczeń
   (kombinacje K_i∧K_j i K4 są PRE-deklarowane wyżej).
2. Tuning ε, g₀, siatek pod wynik; wszystkie skany raportowane w całości.
3. Odrzucanie modu bez testu zbieżności (3 siatki, jak CP-7).
4. Relitygacja: α=2/L04 (#49/#53), mechanizm Koide N=3 (#56),
   werdykty CP-7 (#60).
5. Edycje core/.tex TYLKO addytywne i TYLKO po wynikach (NEEDS →
   user-gate), z wyjątkiem literówek.
6. Zero nowych claimów pozytywnych bez zbieżności + kryterium z LOCK-a.

## 4. Definicje wspólne (jednostki, tolerancje)

γ=β=K_geo=1; TOL_NEG=−1e−6; zbieżność: N∈{2000,4000,8000}, R=60
(kontrola R∈{40,80} dla wyników granicznych); mod zlokalizowany:
λ < −1−10⁻³ (poniżej krawędzi kontinuum −1); solver: eigh_tridiagonal
po symetryzacji wagą r² (jak CP-7).

## 5. Deliverables

- `Phase1_constrained_spectrum.py` + output (+ tabela N_neg^constrained)
- `Phase2_wall_models.py` + output (λ_min(ε), dryf r₂₁/r₃₁(ε), kontakt)
- `Phase3_budget_thresholds.py`/`.md` + output
- `README.md` cyklu (werdykty W1–W3 względem LOCK-a), `NEEDS.md`
  (propozycje edycji core, user-gated), wpis STATE.md (#62+)
- Aktualizacja `audyt/L03_K_phi_stability/` (dyspozycja resztki F-S)
  i ew. T-OP4 w lepton paper (jeśli W1 pozytywne → upgrade z OPEN na
  CONDITIONAL-RESOLVED; jeśli negatywne → utrzymać OPEN + doprecyzować).

## 6. Czego ten cykl NIE robi

- Nie zmienia dopasowań mas korony (r₂₁/r₃₁) — tylko testuje ich
  odporność (W2b).
- Nie dotyka sektora grawitacyjnego (zamknięty w CP-7).
- Nie rozstrzyga L04 (dualizm formulacji) — wyniki mogą dostarczyć
  danych, dyspozycja L04 osobno.

## 7. Deklaracja LOCK

Kryteria, więzy, rodziny regularyzacji i kombinacje zapisano PRZED
jakimkolwiek obliczeniem Phase 1–3. Wyniki negatywne będą zgłoszone
wprost (STATE + README + Limitations).
