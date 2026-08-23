---
title: "op-CG-Kij-from-Hgamma: czy kwartyczne sprzezenie geometryczne K_ij=J(phi_i phi_j)^2 (zrodlo alpha=2) wyprowadza sie z mikro H_Gamma jakimkolwiek coarse-grainingiem, czy jest nieredukowalnym aksjomatem v2?"
date: 2026-06-27
type: research-cycle
folder_status: closed-resolved   # CLOSED 2026-06-27: F-CGK-D = NON-DERIVABLE (patrz Phase_FINAL_close.md)
parent: "[[../../TGP_FOUNDATIONS.md]]"

# ============== KICKOFF CONTRACT (structural cycle) ==============
contract:
  L1_native:
    output_observable: "BRAK (cykl strukturalny; output_type=structural) — werdykt: status derywowalnosci alpha=2 z H_Gamma (DERIVABLE | NON-DERIVABLE | UNDETERMINED)"
    measurement_instrument: "n/d (czysto strukturalny: ERG/wymiar operatora + bound unitarnosci)"
    falsification_rule: "patrz Phase0_LOCK §F (F-CGK-A..D, value-blind, pre-rejestrowane)"
    pre_registration_date: "2026-06-27"
  L2_framework_reduction:
    reduction_type: "not-attempted"
    failure_disposition: "L1-stands"
  L3_falsification_map: []

tgp_status:
  level: T0
  kind: audit            # audyt derywowalnosci (czy aksjomat redukowalny)
  output_type: structural
  core_compatibility: review-only
  may_edit_core: false
  has_needs_file: false
  has_findings_file: false
  exports_findings: false
  open_bridges: ["Gamma->Phi coarse-graining (CG-1..CG-4), most UV/CG"]
  depends_on:
    - "[[../op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] (#49, alpha_eff=-1/2, LOCKED)"
    - "[[../op-bond-order-RG-selection-2026-06-23/Phase_FINAL_close.md]] (#39, alpha=2 RG-irrelevant)"
    - "[[../op-amplitude-density-phase-bridge-2026-06-27/Phase1_Vsub_FINDINGS.md]] (T_A..T_D, V_sub, T_C)"
  impacts: ["status_map l.72 (selekcja na gestosci)", "HONEST_FRAMING_UV_CG_ROOTS (korzen alpha=2)", "#42 ledger (jesli Branch C wymaga nowego parametru)"]
  source_of_status: []

predecessors:
  - "[[../op-amplitude-density-phase-bridge-2026-06-27/]] (diagnostyka R1 + V_sub)"
  - "[[../op-CG-alpha-eff-convergence-2026-06-26/]] (#49)"
  - "[[../op-bond-order-RG-selection-2026-06-23/]] (#39)"
classification: AUDIT
goal: "Rozstrzygnac value-blind, czy alpha=2 (przez K_ij=J(phi_i phi_j)^2) jest wyprowadzalne z mikro H_Gamma (eq:B-H) jakimkolwiek legalnym coarse-grainingiem, czy nieredukowalnym aksjomatem v2. Domknac korzen alpha=2 mapa obstrukcji (kompletna), nie pojedynczym wynikiem."
target_window: "strukturalny; horizon: 1-2 sesje merytoryczne"

six_requirements_target:
  - "P1: werdykt WYLICZONY z flag F-CGK-A..D (0 hardcoded)"
  - "P2: #49 / #39 / status_map l.72 LOCKED — re-derywacja tylko jako baseline value-blind"
  - "P3: endpoint alpha=2 / +4 NIGDY nie wbudowany; wykladniki liczone"
  - "P4: nowy parametr (czteropolowy bond) — jesli potrzebny — JAWNIE liczony do ledgera #42"
  - "P5: bound unitarnosci na eta CYTOWANY (literatura), nie zalozony"
  - "P6: S05 single-Phi axiom preserved bezwarunkowo"

risk_flags:
  - "R1: re-litygacja #49 (mitygacja: F-CGK-A = baseline IMMUTABLE)"
  - "R2: pomylenie 'generowany pod RG' z 'relewantny pod RG' (forbidden #6)"
  - "R3: przemycenie +4 przez niemotywowana redefinicje kompozytu Phi=sigma^{1/3} (forbidden #2)"

phase_plan:
  Phase_0: "Balance sheet + pre-rejestracja F-CGK-A..D + forbidden + methodology read (Phase0_LOCK.md)"
  Phase_1: "F-CGK-A (baseline e=-1, sympy value-blind) + F-CGK-C1 (wymiar RG operatora O_kin=Phi^n(grad Phi)^2 przy 3D WF: relewantny/irrelewantny)"
  Phase_2: "F-CGK-B (wymagane eta vs bound unitarnosci 3D Z2) + F-CGK-C2 (wspolczynnik czteropolowego bondu: fiksowany czy wolny)"
  Phase_FINAL: "Agregat F-CGK-D + dyspozycja + anti-Lakatos sign-off"

tags:
  - cycle-scaffold-2026-06-27
  - CG-Kij-from-Hgamma
  - alpha2-derivability
  - coarse-graining
  - RG-relevance
  - unitarity-eta-bound
  - UV-CG-roots
  - anti-Lakatos-LOCKED
---

# op-CG-Kij-from-Hgamma

> **Cel:** rozstrzygnac, czy `alpha=2` (przez kwartyczne sprzezenie geometryczne
> `K_ij=J(phi_i phi_j)^2`) jest **wyprowadzalne** z mikroskopowego `H_Gamma` (eq:B-H)
> jakimkolwiek legalnym coarse-grainingiem — czy **nieredukowalnym aksjomatem v2**.

## §0 — Pytanie cyklu (jedno zdanie)

`H_Gamma = sum_i[(m0^2/2)s^2 + (lam0/4)s^4] - J sum_<ij> s_i s_j` (eq:B-H) jest 3D-Isingowy
(bilinearny bond → kinetyka kanoniczna w `s`). Manuskrypt v2 **postuluje** osobne sprzezenie
`K_ij=J(phi_i phi_j)^2` → `K(phi)=K_geo phi^4` → `alpha=2` (rem:B-v2-status: „bezposredni
skutek aksjomatu, bez wyprowadzenia MK→GL"). **Pytanie:** czy istnieje JAKIKOLWIEK
(non-Gaussian / strong-coupling / RG) coarse-graining, ktory generuje to kwartyczne sprzezenie
z `H_Gamma` — z fiksowanym wspolczynnikiem i bez nowego parametru?

## §1 — Trzy wyczerpujace galezie (pre-rejestrowane w Phase0_LOCK)

| Galaz | Pytanie | Oczekiwanie (NIE prog) |
|---|---|---|
| **A — baseline Gaussian** | bilinearny bond → jaki `e`? | `e=-1` (re-derywacja #49, value-blind anchor) |
| **B — ucieczka eta** | czy unitarny 3D FP daje `eta` na `De=5`? | REFUTED (eta_Ising≈0.036 << O(5); bound unitarnosci) |
| **C — relewancja / augmentacja** | czy `O_kin=Phi^n(grad Phi)^2` relewantny przy WF? czy czteropolowy bond ma fiksowany wspolczynnik? | AXIOM (Delta_O≈4.8>3 ⇒ irrelewantny; #39) |

**Agregat F-CGK-D:** `DERIVABLE` jesli (B-DERIVABLE ∨ C-DERIVABLE-CONDITIONAL);
`NON-DERIVABLE` (alpha=2 = nieredukowalny aksjomat) jesli (B-REFUTED ∧ C-AXIOM);
`UNDETERMINED` w pozostalych. **Werdykt wyliczany z flag, nie zakladany.**

## §2 — Win-win (dlaczego cykl ma wartosc niezaleznie od wyniku)

- **NON-DERIVABLE** (oczekiwane): domyka korzen `alpha=2` **kompletna mapa obstrukcji**
  (Gaussian → -1/2; eta-escape → unitarnosc; RG → irrelewancja; struktura → T_C). To pelnoprawny
  wynik metodologiczny (analog G_SPA=48), wzmacnia HONEST_FRAMING i `status_map` l.72.
- **DERIVABLE** (gdyby): przelom — `alpha=2` z pierwszych zasad, redukcja `N_axiom`, rewizja
  HONEST_FRAMING. Pre-rejestracja gwarantuje, ze taki wynik bylby przyjety (anti-Lakatos).
- **UNDETERMINED**: identyfikuje DOKLADNIE brakujacy rachunek (relewantny operator + nierozstrzygniety
  wspolczynnik) jako nastepny, wezszy cel.

## Status

🟢 **CLOSED-RESOLVED 2026-06-27 — F-CGK-D = NON-DERIVABLE.** Phase 0 (LOCK) + Phase 1 (A/C1) +
Phase 2 (C2/D) wykonane w 1 sesji. Werdykt: `alpha=2` nieredukowalny aksjomat v2 (kompletna mapa
obstrukcji). Patrz [[Phase_FINAL_close.md]]. #49/#39 LOCKED; rdzeń nietknięty; dyspozycja user-gated.

## Cross-references
- [[Phase0_LOCK.md]] (balance + F-CGK LOCK + forbidden)
- [[../op-amplitude-density-phase-bridge-2026-06-27/Phase1_Vsub_FINDINGS.md]] (V_sub, T_A..T_D)
- [[../op-CG-alpha-eff-convergence-2026-06-26/Phase_FINAL_close.md]] (#49, LOCKED)
- [[../op-bond-order-RG-selection-2026-06-23/Phase_FINAL_close.md]] (#39, RG-irrelevant)
- [[../../meta/HONEST_FRAMING_UV_CG_ROOTS.md]] · [[../../core/_meta_latex/status_map.tex]] l.72
