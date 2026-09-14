---
title: "Phase_method_decisions — klasa dynamiki (A: Cahn–Hilliard M≡1 spektralnie semi-implicit; B: RK4 na (ψ,ψ̇) z B(ψ)=w²K z lapse M9.1'', wyprowadzenie z cytatem); progi/okna/gate'y zamrożone; analiza stabilności dt_B=0.0025 (CFL do ψ≈0.07, sztywność B→0 pre-rejestrowana ψ_stiff=0.12); decyzje ZAMROŻONE przed startem obliczeń"
date: 2026-09-02
type: method-decisions
tgp_owner: research/op-dynamics-class-M911-2026-09-02
status: FROZEN-PRE-COMPUTE
computations_performed: ZERO
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/Phase_method_decisions.md]]"
  - "[[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]"
---

# Decyzje metodyczne (zamrożone PRZED jakimkolwiek obliczeniem cyklu)

Kryteria LOCKa NIETKNIĘTE. Przed zamrożeniem wyłącznie: odczyt LOCKa
i dokumentów poprzednika, analiza wariacyjna/stabilnościowa na papierze
(część zamrożenia; weryfikacje w bramkach Phase 1), arytmetyka
transkrypcyjna rejestru. Środowisko dziedziczone (CPython 3.14.2,
numpy 2.4.3, scipy 1.17.1, sympy 1.14.0).

**Rejestr WEJŚĆ (flagowany):** K_geo=γ=1; seed=20260904, amp=1e−3;
dip: amp 0.45, σ=1.5, ψ_min=0.55; sieć: skalowanie per-siatka do
ψ_max=1.30 (procedura FROZEN poprzednika); M≡1 (wariant A);
B(ψ)=w²K=ψ⁶/(4−3ψ)² (wariant B); dt_A=0.01, dt_B=0.0025 (dt/2:
0.005/0.00125); t_max: A=200, B=100; stacjonarność A ‖ψ̇‖∞≤1e−8;
okno B [80,100], obłożenie ≥80%; ε_H≤0.02 wzgl. max(|H_rel(0)|,0.01);
progi 5/6, 7/6; dev 0.05; podsiatka 5e−3 (tylko A); ψ_band=4/3−1e−6,
ψ_low=1e−6, ψ_stiff=0.12; siatki: gen/dip L=4π N∈{32,48}, lat L=2π
N∈{32,48}; próbkowanie serii Δt=1; checkpointy co 10 j.cz.

## 1. Formy i pochodzenie B(ψ) — CYTAT

Para (w,V,K), 𝒰=wV=γ(ψ⁴/4−ψ³/3), 𝒰′=γψ²(ψ−1), 𝒦=ψ⁴, odczyt B —
dziedziczone DOSŁOWNIE z MD poprzednika §1–3 (cytaty sek08a tamże;
tożsamości sympy: Phase1 poprzednika P1b PASS).

**Waga czasowa B (wyprowadzenie, pre-compute):** z metryki M9.1''
(cytat MD poprzednika §1.2): g_tt=−c₀²(4−3ψ)/ψ ⟹ |g^tt|=ψ/(c₀²(4−3ψ)).
Człon czasowy akcji: ½·w·K·|g^tt|·ψ̇² = ½·[ψ/(4−3ψ)]·ψ⁴·[ψ/(4−3ψ)]·ψ̇²
= ½B(ψ)ψ̇², **B=ψ⁶/(4−3ψ)²** (c₀=1) — identyczna waga jak M=diag(w²K)
zamrożona w MD poprzednika §9 dla Hessianu. B′ = 12ψ⁵(2−ψ)/(4−3ψ)³
(iloraz; weryfikacja sympy 1e−12 = bramka G2). EOM z EL dla
L=∫[½Bψ̇²−½𝒦|∇ψ|²−𝒰]: **Bψ̈ + ½B′ψ̇² = −δE/δψ**;
H=∫[½Bψ̇²+½𝒦|∇ψ|²+𝒰] zachowane (η=0).

## 2. Wariant A — dyskretyzacja (FROZEN)

μ = δE/δψ = −(rhs dziedziczonego silnika gradient flow) — dokładny
gradient dyskretnej E_h (struktura strumieniowa poprzednika, verbatim).
Krok spektralny semi-implicit (analog schematu linii, rząd biharmoniczny):
Δψ̂ = −dt·k²·μ̂ / (1 + dt·A_t·k⁴), A_t = 1.05·max𝒦(ψⁿ), k² = ksym
(dyskretny symbol laplasjanu — ten sam co u poprzednika).
**Zachowanie masy dokładne:** wiersz k=0 daje Δψ̂(0)=0 (gate G1
potwierdza ≤1e−13). ψ̇ raportowane jako ∇²μ (spektralnie);
stacjonarność ‖∇²μ‖∞≤1e−8. Sanity: E niemalejąco maleje
(dE/dt=−∫|∇μ|²; wzrost >1e−12 między próbkami = BREAKDOWN-NUMERIC,
klasyfikacja). Klasyfikacje granic dziedziczone (pas, dolny, niefinit.).

## 3. Wariant B — integrator i stabilność (FROZEN)

System (u,v)=(ψ,ψ̇): u̇=v, v̇=[−δE/δψ − ½B′v²]/B; **RK4**, dt=0.0025.
Start v≡0. Analiza stabilności (pre-compute, papier):
- prędkość falowa c²=𝒦/B=(4−3ψ)²/ψ² — rośnie przy słabym polu;
  granica RK4 (ω_max≈c·2√3/h, |dt·ω|≤2.8): dla h(N=48)=0.2618
  dt=0.0025 stabilne do c≈85, tj. ψ≳0.047;
- lokalna częstość potencjalna |𝒰″|/B ≈ 32/ψ⁵ (ψ→0) — rozdzielczość
  czasowa degraduje od ψ≈0.12 (dt·ω≈2.8) ⟹ **ψ_stiff=0.12**
  pre-rejestrowane w LOCKu (flaga STIFF; załamanie po fladze =
  INCONCLUSIVE-STIFF, rozstrzyga dt/2);
- w próżni ω=√(k²+1), dt·ω≪1 — dokładność RK4 O(dt⁴) (dryf H
  monitorowany co Δt=1; gate ε_H).
H liczone dyskretnie konsystentnie z E_h (człon ½Bv² per komórka).
Bez stopu stacjonarności; stop: nukleacja / pas / załamanie / t_max=100.
δE/δψ — TEN SAM dyskretny gradient co w A (jedna implementacja,
M911_common) ⟹ bramki testują wspólny rdzeń.

## 4. Detektory, starty, werdykty

Dziedziczone verbatim (label+sklejanie periodyczne, N_seed t=0, okno
11 próbek; progi 5/6 i 7/6). Starty: geneza (konstrukcja pasmowa
poprzednika, seed=20260904, wspólne współczynniki N=32/48, normalizacja
z N=48 — UWAGA: poprzednik normalizował z N=64; tu największa siatka
to 48, normalizacja z N=48 [decyzja FROZEN, transkrypcyjna]); dip
(centrum pudła, odległość periodyczna); lat (procedura FROZEN
poprzednika, mtime npz weryfikowany przy każdym odczycie).
Werdykty i kategorie — litera LOCKa §2 (logika składania jak u
poprzednika; dla B kryteria na skalarach okna per siatka — odstępstwo
od porównań polowych pre-rejestrowane w LOCKu). Zdarzenie ⟹ dt/2
na obu siatkach pary.

## 5. Higiena

Pełne ścieżki bez `cd`; `ls` po każdym zapisie; biegi B w tle
z checkpointami npz co 10 j.cz., wznawialne `--resume` (proces ≤~55 min
— batch per start); wyniki per-bieg json/npz w `Phase{2,3}_results/`;
`--verdict` składa do `Phase{2,3}_output.txt` + `*_relaxed_states.npz` /
`*_final_states.npz`. Rdzeń .tex/STATE/git nietykane; katalogi innych
cykli tylko odczyt. INCONCLUSIVE/STIFF/BOUNDARY ≠ pozytyw.

**FROZEN 2026-09-02, przed uruchomieniem jakiegokolwiek skryptu cyklu.**
