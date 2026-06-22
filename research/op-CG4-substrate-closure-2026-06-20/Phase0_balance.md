---
title: "Phase 0 balance sheet — op-CG4-substrate-closure (forward, mandatory gate)"
date: 2026-06-20
parent: "[[README.md]]"
type: balance-sheet
cycle: op-CG4-substrate-closure-2026-06-20
auditor: claudian
classification: PENDING (do wyliczenia w Phase FINAL z reguły LOCKED Phase0_lock §4)
anti_lakatos_lock: ACTIVATED
tags: [phase0, balance-sheet, gate, CG-4, kappa-E, anti-Lakatos]
---

# Phase 0 balance sheet — op-CG4-substrate-closure

> **Gate (CALIBRATION_GATE_ENFORCEMENT, ABSOLUTE BINDING):** ten plik istnieje **przed** jakimkolwiek
> zapisem do PREDICTIONS_REGISTRY. Klasyfikacja końcowa = wyliczona w Phase FINAL z reguły agregatu
> [[Phase0_lock.md]] §4. Tu: ramy 6 sekcji + testy zalockowane.

## 1. Co cykl twierdzi że robi

> Wyznaczyć **scheme-independent** $C_\sigma$ (sztywność kinetyczną tensora σ_ab) przez znalezienie
> niepatologicznego modelu substratu, i przez to **rozstrzygnąć** $\kappa_E=8\pi G_0 C_\sigma\sigma_0^2/c^3$
> względem progu przeżycia $5/6$ → twardy werdykt sektora radiacyjnego GW.

Główne claims (do weryfikacji, NIE z góry przyjęte):
- (C1) istnieje niepatologiczny substrat (stabilny + czysty punkt krytyczny Z₂);
- (C2) operator złożony $O_{ab}=\partial\hat s\,\partial\hat s$ renormalizuje się scheme-independent na tym substracie;
- (C3) wynikowe pasmo $\kappa_E$ rozróżnia $5/6$ od $1$ (faktor 1.2).

## 2. Phase 0 balance sheet (CALIBRATION_PROTOCOL §2)

### 2.1 External inputs (observational)
```
- GWTC-3 combined ~90 BBH posterior → constraint na modyfikację radiacyjną
  (parent: M9.1'' specific (4-3ψ)/ψ RULED OUT 5.02σ; β_ppE bound |β|≤0.78)   [obs, op-GWTC3-reanalysis]
- Próg przeżycia κ_E = 5/6: bilans Ṗ_b = κ_E P_GR + (1/6) P_GR = P_GR ⟺ κ_E=5/6
  (skalar konforemny 1/6 nieunikniony; #30)                                    [strukturalny, exact]
- G_0, c — stałe TGP w jednostkach LOCKED (D01 anchor manifest)               [CODATA-pochodne, locked]
```

### 2.2 Structural axioms (TGP-internal LOCKED — z source cyklem)
```
- Φ = ⟨ŝ²⟩ pole kanoniczne (gęstość)              [source: sek01 def:Phi]          status: LOCKED
- Φ_0 = 2ρ_0* = 0.0609, ν = 0.749                  [source: CG-2 tgp_erg_lpa_prime] status: DERIVED_NUM
- σ_ab = composite ⟨(∂_a ŝ)(∂_b ŝ)⟩, M²=2m_s²      [source: closure_2026-04-26 Path B] status: STRUCTURAL
- C_σ > 0, O(1) (znak + skaling)                   [source: op-Csigma-coarse-graining] status: STRUCTURAL
- κ_E = 8πG_0 C_σ σ_0²/c³                           [source: op-sigma-kinetic-Csigma]  status: STRUCTURAL
- α=2 = postulat na gęstości (NIE derywacja)        [source: sesja #32 rem:alpha2-pivot-status] status: AXIOM-SELECTION
- bąbel 3D Π(p)=1/(8πm) − p²/(96πm³)               [source: op-Csigma-lattice-MC]    status: DERIVED (sympy)
```

### 2.3 Derived outputs (the cycle claims)
```
- Output 1: C_σ scheme-independent (wartość + pasmo)     (claim: Phase 3)
- Output 2: κ_E band (central + [lo, hi])                (claim: Phase 3)
- Output 3: werdykt sektora radiacyjnego (SURVIVE / FALSIFIED-hard / PARTIAL / GAP)  (claim: Phase FINAL, wyliczony)
```

### 2.4 Tautology test (CRITICAL)
**Pytanie:** czy $\kappa_E$ jest wyrażalne tylko przez inputy/aksjomaty bez redukcji do tożsamości?
- $\kappa_E=8\pi G_0 C_\sigma\sigma_0^2/c^3$ zawiera $C_\sigma$ — **wielkość mierzoną niezależnie** (MC operatora na substracie), nie funkcję progu 5/6.
- **Krytyczne zabezpieczenie:** próg 5/6 NIE wchodzi do pomiaru $C_\sigma$ (forbidden #1 — zero strojenia). Pomiar i próg są rozłączne ⟹ porównanie ma treść.
- **Werdykt tautology test:** PASS (warunkowo — re-test w Phase 3 po pomiarze: czy $C_\sigma$ nie wciągnął ukrytego 5/6).

### 2.5 Falsifiability test (CRITICAL)
**Konkretny falsifier (dwustronny):**
```
Jeśli scheme-indep. κ_E (pasmo <1.2×) wyklucza 5/6 → sektor FALSIFIED-hard.
Jeśli pokrywa 5/6 i wyklucza 1            → sektor SURVIVE.
Aktualna pozycja: pasmo [0.04, 11.1] (#31) — NIErozróżnia (PARTIAL).
```
**Band check:** cel = $C_\sigma$ band < faktor 1.2 (odstęp 1↔5/6 = 1.2×). Falsyfikowalny **z definicji bramy C-D**.
- **Werdykt falsifiability test:** PASS (brama C-D wymusza falsyfikowalność; PARTIAL jeśli band za szeroki).

### 2.6 Independent-path cross-validation (CRITICAL dla DERIVED)
**Path 1 (numeryczny):** MC operatora złożonego na niepatologicznym substracie → $C_\sigma$ z structure factor + continuum.
**Path 2 (analityczny):** bąbel 3D $\Pi(p)$ (op-Csigma-lattice-MC) + renormalizacja perturbacyjna współczynnika $p^2$ → $C_\sigma$ z subtrakcją UV.
**Convergence:** wymagana zgodność Path 1 ↔ Path 2 w granicach pasma <1.2× (do sprawdzenia Phase 3).
- **Werdykt independent-path:** PENDING (≥2 ścieżki **zaplanowane**; PASS dopiero po konwergencji Phase 3).

## 3. Audit gate checklist (stan na Phase 0)
```
☑ Phase 0 balance sheet exists (this file)
☐ Tautology test PASS                      — wstępnie PASS; re-test Phase 3
☐ Falsifiability test PASS                 — brama C-D wymusza; weryfikacja Phase 3
☐ Independent-path (≥2 physical paths)      — zaplanowane (num + analit.); konwergencja Phase 3
☑ Alt-scan ≥4 kandydatów (M0–M3, Phase0_lock §7) z argumentem strukturalnym (forbidden #7)
☑ NIE post-hoc structural motivations       — kryteria zalockowane przed liczbami
☑ NIE constructed criterion                 — próg 5/6 odziedziczony (#30), nie skonstruowany tu
☑ NIE circular anchor                       — C_σ mierzone, nie funkcja 5/6 (forbidden #1)
☑ NIE inheriting drift > parent × 5         — budżet nowych stałych 0; reuse LOCKED
```
**Stan:** brama wymaga 3 testów Phase 3 (tautology re-test, falsifiability, independent-path). Do tego czasu
**max status = STRUCTURAL** (per gate: brak choćby jednego ☐ → max STRUCTURAL).

## 4. Klasyfikacja końcowa
**PENDING** — wyliczona w Phase FINAL z reguły [[Phase0_lock.md]] §4. Możliwe wyniki:
- werdykt sektora FALSIFIED-hard / SURVIVE → output klasy DERIVED_CONDITIONAL (po konwergencji 2 ścieżek);
- PARTIAL / GAP → output klasy STRUCTURAL (jak #31).

## 5. Recommended action (na teraz)
- [x] **PROCEED** — Phase 0 zalockowany; przejść do Phase 1 (silnik analityczny) po autoryzacji.
- [ ] registry write — **zablokowane** do Phase FINAL (gate).

## 6. Cross-references
- [[Phase0_lock.md]] — kryteria, reguła agregatu, forbidden moves
- [[README.md]] — opis cyklu
- [[../../meta/CALIBRATION_GATE_ENFORCEMENT.md]] — gate (ABSOLUTE BINDING)
- [[../op-Csigma-lattice-MC-2026-06-14/Phase_FINAL_close.md]] — parent (κ_E≈0.62, pasmo [0.04,11.1])
