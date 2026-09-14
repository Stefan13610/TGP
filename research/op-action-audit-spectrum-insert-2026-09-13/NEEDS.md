---
title: "NEEDS — op-action-audit-spectrum-insert: wszystkie pozycje USER-GATED (drzewo §5 LOCKa); N1 poprawka definicji ΔE_insert (więz skończonej skali — punktowy pin ma zerową pojemność w 3D); N2 user-gate CORE: napięcie znaku wariacji (|g^tt| ⟹ ω²=k²+1 ale statyka (ψ−1)/ψ²; znakowana ⟹ R3 ODE ale tachion k²−1); N3 cykl dynamiki 2. rzędu z zalockowanym (M,𝒦,𝒰)"
date: 2026-09-13
type: needs
tgp_owner: research/op-action-audit-spectrum-insert-2026-09-13
status: USER-GATED
related:
  - "[[Phase0_balance.md]]"
  - "[[Phase_FINAL_close.md]]"
  - "[[../op-metric-pair-M911-2026-09-02/NEEDS.md]]"
---

# NEEDS (user-gated; żadna pozycja nie startuje bez decyzji usera)

Mapowanie na drzewo decyzyjne LOCKa §5: Q-D1-PASS → gałąź „operator
kinetyczny zalockowany dla dynamiki 2. rzędu"; Q-D2-INCONCLUSIVE →
gałąź „NEEDS metodologiczny: poprawka definicji (więz, brzeg, skala
pudła) — bez claimów o znaku".

## N1 (z Q-D2-INCONCLUSIVE): poprawka definicji ΔE_insert — więz skończonej skali

Punktowy więz ψ(0)=A ma zerową pojemność w 3D: minimizer warunkowy
zapada się do rdzenia skali siatki i ΔE∝h→0⁺ (potwierdzone dokładnymi
minimami Newtona, stosunek 0.4988 przy h→h/2 — Phase3_diag_output.txt).
Wymagany NOWY LOCK z definicją o dobrze określonej granicy kontinuum,
np. (do decyzji usera, nie przesądzam): (a) pin na kuli skończonego
promienia ψ|_{r≤ρ}=A z ρ zalockowanym (dodatnia pojemność), (b) więz
całkowy/projekcyjny na zamrożoną funkcję testową, (c) rodzina
zamrożonych profili bez relaksacji rdzenia. Aspekt R definicji
obecnej DZIAŁA (zero członu objętościowego dokładnie — zachować).
Uwaga dziedziczona: A=1.30 → BREAKDOWN-BOUNDARY-LOWER (mechanizm
½𝒦′|∇ψ|² na kolanie skali siatki) — więz skończonej skali usuwa
także ten artefakt ostrza.

## N2 (z kaweatu Phase 1–2; user-gate CORE — ZERO samowolnych napraw)

Napięcie znaku wariacji w rdzeniu, zlokalizowane precyzyjnie
(Phase1_output.txt, sympy obie redukcje):
- odczyt |g^tt| (zamrożony w LOCKu tego cyklu): L=½Mψ̇²−½𝒦|∇ψ|²−𝒰 ⟹
  ω²=k²+1 (zdrowa próżnia, Q-D1-PASS), ale statyka = (ψ−1)/ψ²
  (zaniki Yukawy) — NIE R3 ODE;
- literalna kontrakcja znakowana (sygnatura (−,+,+,+) zapisu
  eq:metric-M911-canonical z +½Kg^{μν}∂ψ∂ψ): odtwarza R3 ODE
  (1−ψ)/ψ² (fundament spektrum mas R3), ale ω²=k²−1 — dokładnie
  tachion audytu rozdz. 2, plus duch w członie czasowym.
Jedna rzeczywista wariacja Lorentzowska nie daje obu naraz. Decyzja
usera o procedurze wariacji/sygnaturze/konwencji zapisu w rdzeniu
(sek08a/sek08c) — dopiero potem ewentualny dopisek core (sesja
główna; rdzeń w tym cyklu NIETKNIĘTY).

## N3 (z Q-D1-PASS): cykl dynamiki 2. rzędu (NEEDS N3 poprzednika staje się dobrze postawiony)

Operator kinetyczny ZALOCKOWANY tym cyklem (warunkowo na konwencji
|g^tt| — patrz N2): M=K_geoψ⁶/(c₀(4−3ψ)²), 𝒦=K_geoψ⁴,
𝒰=γ(ψ⁴/4−ψ³/3); π=Mψ̇; EOM Mψ̈+½M′ψ̇²=−δE_PRIMARY/δψ; m²=c_s²=1
jako punkt odniesienia (nośnik propagacji P1.3, sondy P2). Osobny
LOCK (starty, dyssypacja vs zachowawczość, diagnostyki).

## N4 (metodologiczny, niższy priorytet): schemat relaksacji z więzem

Punkt stały zamrożonego schematu „krok semi-implicit → projekcja
pinu" różni się od dokładnego minimum warunkowego (tożsamość:
rhs₁=−d·A_t·a₀/m₁; przesunięcie ΔE +4.4% przy A=0.70, ~2× przy
A=1.25) i uniemożliwia formalne osiągnięcie progu stacjonarności
(statusy TMAX przy plateau E). Przy każdym przyszłym LOCKu z więzem:
egzekwowanie więzu WEWNĄTRZ kroku implicit (wiersz δ₀=0) albo
kryterium stacjonarności na dH węzłów swobodnych.
