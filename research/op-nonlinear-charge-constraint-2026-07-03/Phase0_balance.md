---
title: "Phase0_balance — LOCK cyklu op-nonlinear-charge-constraint (stabilność μ/τ przy zachowanym ładunku: właściwy Q-ball zamiast liniowej deflacji)"
date: 2026-07-03
type: phase0-lock
tgp_owner: research/op-nonlinear-charge-constraint-2026-07-03
status: "LOCKED — realizacja: NASTĘPNA SESJA (osobny agent-implementator; user-authorized 2026-07-03: 'rozpisz cykl badawczy N4 dla nowego agenta')"
anti_lakatos_lock: PRESERVED
tags: [charge-constraint, Q-ball, Vakhitov-Kolokolov, orbital-stability, wall-dynamics, L03, L04, korona]
related:
  - "[[../op-wall-dynamics-2026-07-03/README.md]]"
  - "[[../op-wall-dynamics-2026-07-03/NEEDS.md]]"
  - "[[../op-spectral-analysis-Phi-2026-07-03/README.md]]"
  - "[[../../core/sek08b_ghost_resolution/sek08b_ghost_resolution.tex]]"
---

# Phase 0 — LOCK: op-nonlinear-charge-constraint

## 0. Kontekst dla agenta (handoff — przeczytaj PRZED czymkolwiek)

**Skąd ten cykl.** Łańcuch wyników:

1. **CP-7** (`op-spectral-analysis-Phi-2026-07-03`, CLOSED): solitony
   korony μ/τ (F-S: f(g)=1+2α·ln g, α=2; W′(g)=g²(1−g)) są punktami
   siodłowymi funkcjonału statycznego BEZ więzów — mody zlokalizowane
   l=0: e: 0, μ: 2 (−1,282; −1,057), τ: 3 (−4,216; −1,114; −1,010);
   próżnia g=1 tachioniczna (kontinuum od −1). Jedn. γ=β=K_geo=1.
2. **op-wall-dynamics** (#62, CLOSED): hipoteza budżetowa autora w
   wersji LINIOWEJ obalona — więzy K1–K4 i pary NIE zerują indeksu
   (min. osiągalne μ→1, τ→1; mod rezydualny ≠ kierunek rodziny,
   overlap ≈0,005). ALE dwa fakty strukturalne:
   (a) więz budżetu RDZENIA (K4, r<r_core) usuwa dokładnie mody
   GŁĘBOKIE (−1,282/−4,216) — kierunek najgłębszej niestabilności
   żyje w podprzestrzeni zmian budżetu rdzenia;
   (b) g_crit=8/5 (H7) = próg aktywacji ściany g* (0,71%) — dwa progi
   to jeden mechanizm ścienny.
   Ponadto: τ nie ma reprezentanta EL w gładkiej rodzinie f_ε (kolaps
   przy każdym ε); jedyny N-zbieżny punkt spektralny soft-wall: μ przy
   ε=0,2 (λ_min=−1,389).
3. **Dlaczego wersja nieliniowa może działać tam, gdzie liniowa
   zawiodła (fakt matematyczny, motywacja cyklu):** dla więzu
   LINIOWEGO B[g] człon mnożnika w drugiej wariacji znika (B″=0) —
   operator się nie zmienia, zmienia się tylko dopuszczalna
   podprzestrzeń (to policzył #62: za mało). Dla ŁADUNKU NIELINIOWEGO
   Q[g,ġ] (typu Noether): (i) tło spełnia E′=ωQ′ (profile zależą od
   ω), (ii) hessian dostaje człon −ωQ″, (iii) stabilność orbitalna
   rządzi się teorią Vakhitova–Kolokolova / Grillakisa–Shataha–
   Straussa (znak dQ/dω), a nie gołym indeksem siodłowym.

**Hipoteza robocza cyklu (autor, przeniesiona z #61/#62, zapisana
PRZED analizą):**
> Budżet tworzonej przestrzeni jest wielkością zachowywaną
> (ładunkiem); stabilność μ/τ jest stabilnością orbitalną przy
> ustalonym ładunku (właściwy Q-ball), a nie minimalizacją bez więzów.

**Narzędzia do reuse (wszystkie zweryfikowane w #60/#62):**
- `../op-wall-dynamics-2026-07-03/Phase1_constrained_spectrum.py`
  (funkcje: `soliton_profile`, `build_tridiag`, rachuba inercji
  Haynswortha + bisekcja, `constraint_matrix`, walidacja dense-vs-inertia)
- `../op-wall-dynamics-2026-07-03/Phase2_wall_models.py`
  (`make_feps`, `soliton_profile_soft`, `build_tridiag_gen`, `fit_tail`)
- `../../tooling/scripts/a3d_soliton_brannen_r.py` (konwencja
  A_tail → masy: m ∝ A_tail⁴, fit ogona r∈[20,35])
- Stałe: G0_E=1.24915, G0_MU=φ·G0_E, G0_TAU=3.18912; g*=e^{−1/4};
  guard odbicia g*+0,005 (konwencja a3d).

Kolejność czytania dla agenta:
1. Ten LOCK (całość).
2. `../op-wall-dynamics-2026-07-03/README.md` + Phase1/2 outputy.
3. `../op-spectral-analysis-Phi-2026-07-03/README.md` (CP-7).
4. `core/sek08b…tex`: `rem:wall-dynamics-2026-07-03`,
   `rem:spectral-CP7`, `rem:ghost-artifact-scope-CP7`.
5. STATE.md wpisy #60–#62.

## 1. Cel cyklu

Rozstrzygnąć, czy istnieje ZACHOWANY ładunek (z dynamiki substratu
lub jej minimalnego rozszerzenia), przy którego ustaleniu solitony
μ/τ są orbitalnie stabilne (kryterium VK/GSS) — tj. czy hipoteza
budżetowa autora jest prawdziwa na poziomie nieliniowym/ładunkowym.

## 2. Modele i kandydaci (ZAMKNIĘTA lista — żadnych dodatków post-hoc)

- **M0 (kanoniczne zanurzenie dynamiczne, pole rzeczywiste):**
  L_S = ½f(g)ġ² − ½f(g)|∇g|² − W(g). Spójne z E_S (część statyczna);
  wybór modelu FLAGOWANY (dualizm L04 — nie rozstrzygamy L04, tylko
  deklarujemy zanurzenie).
- **M1 (kompleksyfikacja, HIPOTEZA-ROZSZERZENIE teorii):**
  ψ = g·e^{iθ}, L = ½f(|ψ|)|∂_t ψ|² − ½f(|ψ|)|∇ψ|² − W(|ψ|);
  U(1) ⇒ ładunek Noether Q (postać wyprowadzić w Phase 1, sympy).
  Interpretacja kandydacka: faza = potencjał przepływu substratu
  (para Madelunga gęstość+przepływ); Q = liczba kwantów substratu =
  „budżet tworzonej przestrzeni". M1 NIE wchodzi do core bez
  user-gate (to rozszerzenie ontologii).
- **Kandydaci na ładunek w M0 (lista zamknięta):**
  C1: ∫(g−1)r²dr; C2: ∫f(g)ġ r²dr (pęd kanoniczny);
  C3: ∫(g−1)² r²dr; C4: ∫h(g)ġ r²dr dla h∈{1, g²};
  C5: ładunki skalowania/translacji (standardowy inwentarz Noether).
- **Model ściany (zamknięty):** baseline = hard-wall z odbiciem
  (konwencja CP-7/a3d); cross-check WYŁĄCZNIE f_ε przy ε=0,2 (jedyny
  N-zbieżny punkt z #62 W2b) + kontrola ε=0,1. Zakaz zmiany modelu
  ściany post-hoc.

## 3. Fazy i zadania (LOCKED)

### Phase 1 — inwentarz ładunków (analityczna, sympy; zero numeryki dopasowań)

- **P1a:** EOM dla M0 (sympy z L_S); weryfikacja exact zachowania
  energii; tożsamość statyczna EL(M0) = ODE korony (kontrola spójności
  z CP-7 Phase 1).
- **P1b:** dla KAŻDEGO kandydata C1–C5: sympy-werdykt
  ZACHOWANY/NIEZACHOWANY pod EOM M0 (d/dt + strumień brzegowy);
  pełna tabela, bez selekcji.
- **P1c:** M1: wyprowadzić ładunek Noether Q i EOM; sympy-check
  dQ/dt=0 exact; zredukować ansatz stacjonarny ψ=φ(r)e^{−iωt} do ODE
  profilu z parametrem ω; wyprowadzić Q(ω), E(ω) jako całki z φ_ω.
- **P1d (pytanie zalockowane):** czy człon ω² podnosi krawędź
  tachionicznego kontinuum wokół g=1 — wyznaczyć (sympy) krawędź
  σ_ess jako funkcję ω: c(ω)=W″(1)/f(1)+c₂ω²; czy istnieje
  ω_min z σ_ess ≥ 0?

**Kryterium V1 (LOCKED):** pełna tabela C1–C5 z werdyktami. Jeżeli
ŻADEN budżetopodobny ładunek nie jest zachowany w M0 → wynik zapisany
wprost: „hipoteza wymaga rozszerzenia M1" (dalej Phase 2 wyłącznie w
gałęzi M1, jawnie oznaczonej jako model-rozszerzenie).

### Phase 2 — rodzina Q-ball i kryterium VK (numeryczna, główna)

- **P2a:** rodziny profili φ_ω (shooting, DOP853 rtol 1e−10, konwencja
  ściany z §2) dla ω w siatce [0, ω_max] (ω_max = próg z P1d lub 1,0;
  krok ≤0,05; skan raportowany W CAŁOŚCI), starty z g₀^(e/μ/τ);
  ciągłość gałęzi ω→0 z profilami statycznymi CP-7.
- **P2b:** Q(ω), E(ω); test VK: znak dQ/dω na gałęziach; punkty zmiany
  znaku (jeśli są).
- **P2c:** spektrum operatorów fluktuacji wokół φ_ω (rozkład GSS na
  L₊/L₋ w M1; dyskretyzacja tridiagonalna jak CP-7; zbieżność
  N∈{2000,4000,8000}, R=60, kontrola R∈{40,80} dla wyników
  granicznych); N_loc z więzem ładunkowym (deflacja kierunku
  fazowego/rodziny — reuse rachuby inercji z #62).
- **P2d:** kontrola mas: w granicy ω→0 odtworzyć r₂₁, r₃₁ z dryfem
  <0,1% (baseline hard-wall #62: r₂₁=206,73, r₃₁=3479,6); przy
  ω=ω_k(stabilizującym, jeśli istnieje) — dryf r₂₁/r₃₁ raportowany
  w całości (bez progu akceptacji: to dane do interpretacji, nie
  kryterium).

**Kryterium V2 (LOCKED):** „stabilizacja ładunkowa POTWIERDZONA dla
generacji k" ⟺ koniunkcja:
(i) istnieje ω_k>0 z profilem na gałęzi ciągle połączonej z g₀^(k);
(ii) VK: dQ/dω < 0 w ω_k (konwencja: stabilna gałąź slope-negative;
     jeśli wyprowadzenie sympy w P1c da odwrotną konwencję znaku dla
     tej normalizacji — udokumentować PRZED obliczeniami P2b);
(iii) N_loc(L₊, deflacja fazy/rodziny) = 0 przy ω_k, zbieżne na
     3 siatkach.
Niespełnienie dla μ LUB τ przy wszystkich ω ze skanu → **wynik
NEGATYWNY zgłaszany wprost** (hipoteza budżetowa obalona także w
wersji ładunkowej M1; dokumentować, co dokładnie zawiodło).
UWAGA: kryterium ma sens TYLKO jeśli P1d/P2c wykażą usunięcie lub
podniesienie kontinuum; jeśli kontinuum pozostaje tachioniczne przy
każdym ω — raportować „nierozstrzygalne w M1 z powodu tła" (to też
wynik).

### Phase 3 — nieliniowy test dynamiczny (M0, niezależny od gałęzi M1)

- **P3a:** ewolucja czasowa μ w modelu f_{ε=0,2} (jedyny N-zbieżny;
  kontrola ε=0,1): g(r,0)=g_eq+a·v_deep, ġ(r,0)=0, a∈{±0,01;±0,03};
  schemat symplektyczny (leapfrog) lub RK4, zbieżność dt (2 kroki),
  gate: |ΔE|/E < 1e−6 przez cały run; R=60, r∈siatka jak CP-7.
  τ: BEZ reprezentanta EL w f_ε (kolaps, #62) — poza zakresem P3;
  odnotować jawnie.
- **P3b:** pomiar a(t)=⟨v_deep, δg⟩: (i) tempo wzrostu vs
  √|λ_min|=√1,389 z teorii liniowej (walidacja); (ii) saturacja:
  poziom nasycenia/rekurencje do t=200.

**Kryterium V3 (LOCKED):** (i) wzrost wykładniczy zgodny z λ_min
(±20%) i brak saturacji poniżej ‖δg‖~10% tła → „niestabilność
potwierdzona nieliniowo w M0-f_ε"; (ii) saturacja <10% + rekurencje →
„efektywna stabilizacja nieliniowa" (dokumentować mechanizm, zero
claimów bez analizy). Każdy kierunek = PASS wykonania.

## 4. Forbidden moves (anti-Lakatos)

1. Zmiana kryteriów/tolerancji/list (C1–C5, M0/M1, ω-siatka, model
   ściany) PO uruchomieniu obliczeń.
2. Tuning ω, g₀, ε, siatek pod wynik; wszystkie skany raportowane
   w całości; ZERO re-fitowania g₀ i stałych korony.
3. M1 (kompleksyfikacja) NIE wchodzi do core bez user-gate — to
   rozszerzenie ontologii; wyniki gałęzi M1 zawsze oznaczone
   „model-extension".
4. Odrzucanie/akceptowanie modu bez zbieżności 3 siatek (jak CP-7);
   niezbieżność raportowana JAKO niezbieżność.
5. Relitygacja: α=2/L04 (#49/#53), Koide N=3 (#56), werdykty CP-7
   (#60), werdykty W1–W3 (#62).
6. Edycje core/.tex TYLKO addytywne i TYLKO po wynikach
   (NEEDS → user-gate), z wyjątkiem literówek.
7. Zero nowych claimów pozytywnych bez spełnienia pełnej koniunkcji V2.

## 5. Definicje wspólne (jednostki, tolerancje)

γ=β=K_geo=1; TOL_NEG=−1e−6; mod zlokalizowany: λ < −1−10⁻³ względem
krawędzi kontinuum wyznaczonej dla danego ω (P1d!); zbieżność:
N∈{2000,4000,8000}, R=60 (kontrola R∈{40,80}); overlap kierunku
rodziny/fazy: próg 0,9; dryf mas ω→0: próg 0,1%; solver spektralny:
eigh_tridiagonal po symetryzacji wagą r² (jak CP-7); ODE: DOP853,
rtol 1e−10, atol 1e−13, max_step 0,02.

## 6. Deliverables

- `Phase1_charge_inventory.py` + output (tabela C1–C5, M1: Q, EOM,
  krawędź σ_ess(ω))
- `Phase2_qball_family.py` + output (gałęzie φ_ω, Q(ω), E(ω), VK,
  N_loc(ω), dryf mas)
- `Phase3_nonlinear_evolution.py` + output (a(t), gate energii,
  werdykt V3)
- `README.md` cyklu (werdykty V1–V3 względem LOCK-a), `NEEDS.md`
  (user-gated), wpis STATE.md (#63+), aktualizacja
  `audyt/L03_K_phi_stability/` (dyspozycja resztki F-S po tym cyklu)

## 7. Czego ten cykl NIE robi

- Nie zmienia dopasowań mas (r₂₁/r₃₁) — kontroluje tylko ich
  odtworzenie w granicy ω→0 i raportuje dryf przy ω>0.
- Nie rozstrzyga L04 (M0 to zadeklarowane zanurzenie, nie teza).
- Nie dotyka sektora grawitacyjnego F-A (zamknięty w CP-7).
- Nie wprowadza M1 do rdzenia (user-gate; najpierw wyniki).
- Nie rozstrzyga pierwszoprincypialnego modelu ściany (poza zakresem;
  konwencja zamknięta w §2).

## 8. Deklaracja LOCK

Modele (M0/M1), lista ładunków (C1–C5), siatka ω, konwencja ściany,
kryteria V1–V3 i progi zapisane PRZED jakimkolwiek obliczeniem.
Wyniki negatywne będą zgłoszone wprost (STATE + README + Limitations).
Jedyna dopuszczalna korekta przed startem obliczeń: udokumentowana
konwencja znaku VK z P1c (patrz V2(ii)) — po starcie P2b również
zamrożona.
