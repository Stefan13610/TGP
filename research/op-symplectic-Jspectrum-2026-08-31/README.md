---
title: "op-symplectic-Jspectrum — mini-cykl: klasa dynamiki substratu a znak W (spektrum JL̂ na tłach łańcucha)"
date: 2026-08-31
type: research-cycle
status: CLOSED
verdict: "Q-FAIL (PRIMARY; wszystkie 4 tła: Re λ = +0.1396…+0.1434 zbieżnie ≫ tol; próżnia symplektyczna niestabilna analitycznie z max Re λ = γ/4; kontrole Gate-A/C1/C2/C3 czyste; odczyt literalny tol zdegenerowany — flagowany)"
tgp_owner: research/op-symplectic-Jspectrum-2026-08-31
authorization: "User 2026-08-31: 'ok, załuż nowy mini cykl'"
related:
  - "[[Phase0_balance.md]]"
  - "[[../op-bloch-chain-stability-2026-08-31/Phase_FINAL_close.md]]"
  - "[[../op-lattice-bath-runaway-2026-08-23/ANALIZA_N2_znak-W-z-akcji_2026-08-23.md]]"
---

# op-symplectic-Jspectrum (2026-08-31)

**Status: CLOSED — werdykt Q-FAIL** ([[Phase_FINAL_close.md]]).

Mini-cykl analityczno-decyzyjny (wzorzec L04) po weryfikacji konwencji
W_source/V_energy: relabeling jest widmowo inwariantny — statyka
maszynerii 2 (ogony oscylacyjne, A2) pinuje ω²=k²−γ w każdej dynamice
II rzędu i w gradient flow. Jedyna nieprzebadana klasa, w której ujemny
kierunek Hessianu NIE implikuje niestabilności, to dynamika
symplektyczna I rzędu (NLS/GP-podobna; mechanizm stabilności orbitalnej,
świat kryterium VK używanego w #63).

**Pytanie binarne:** czy tła łańcucha sektora tachionowego
(d ∈ {3π, 4π, 6π} z cyklu bloch, λ_min(L₊) ≈ −1.22) mają widmo
symplektyczne σ(JL̂) ⊂ iℝ?

Q-PASS ⟹ „dwa sektory znaku W" mogłyby być „dwiema DYNAMIKAMI jednej
akcji" — pierwszy pozytywny nośnik stabilności maszynerii 2, materiał
do decyzji aksjomatycznej autora. Q-FAIL ⟹ klasa symplektyczna też nie
ratuje 1D.

Kryteria, kontrole (NLS stabilny/niestabilny — osiągalne FAIL-e w obie
strony), progi i forbidden moves: [[Phase0_balance.md]].

## Log faz

- 2026-08-31: Phase 0 LOCK zapisany. Obliczenia: NIEROZPOCZĘTE.
- 2026-08-31: `Phase_method_decisions.md` FROZEN przed startem (operatory
  L₊/L₋ + czynnik ½, dyskretyzacja on-shell Q₋, ruling tol — odczyt
  PRIMARY pasmowy vs literalny, siatka 9 k, konstrukcja bramek).
- 2026-08-31: **Phase 1 PASS 16/16** (sympy-exact): P1a inwariantność
  ω²=k²−γ dla dowolnej wagi F; P1b gradient flow ⟺ ten sam Hessian;
  P1c L₋g_d=0 on-shell EXACT, Q₊ ≡ operator bloch, λ²_vac = −¼κ²(κ²−γ)
  ⟹ próżnia niestabilna też symplektycznie (max Re λ = γ/4).
  Korekta implementacyjna skryptu ([[Phase1_correction_note.md]];
  pierwotny output zachowany).
- 2026-08-31: **Phase 2 PASS 4/4**: kotwica λ_min(3π) = −1.222191
  (|δ|=1.19e−7); C1 stabilny: max Re λ = 1.5e−7 ≤ tol; C2 niestabilny
  wykryty: λ = 2.9152 zbieżnie; C3 próżnia vs analityka 6.8e−4 ≤ 1e−3
  (ratio 4.00). Korekta diagnostyki deskryptywnej
  ([[Phase2_correction_note.md]]; pierwotny output zachowany).
- 2026-08-31: **Phase 3 (centralna): Q-FAIL** — max Re λ (N=800, wszystkie
  ZBIEŻNE, płasko w k): 3π +0.143396; 4π +0.139810; 6π 2-garb +0.143413;
  6π 1-garb +0.139571 — wszystkie ≫ tol PRIMARY ≈ 1.2e−2; cross-check
  produktowy λ²=−ν/4 zgodny ≤ 3.4e−8. Odczyt literalny tol (zdegenerowany)
  dawałby formalny „Q-PASS" — oba raportowane, PRIMARY rządzi.
- 2026-08-31: Phase 4 NIE uruchomiona (LOCK: TYLKO przy Q-PASS).
- 2026-08-31: `Phase_FINAL_close.md`, `NEEDS.md` (user-gated). CLOSED.
