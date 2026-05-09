---
title: "Phase 0 — balance sheet: emergent g_eff from many-body Φ-interaction"
date: 2026-05-09
type: phase-balance
status: 🟢 OPEN
parent: "[[./README.md]]"
phase: 0
gate: 8/8
tags:
  - phase0
  - balance-sheet
  - emergent-metric
---

# Phase 0 — balance sheet

## Cel cyklu (locked)

Wyprowadzić `g_eff^μν(x) = G[{Φ_i(x_i)}, σ_ab(x), Φ̄(x), x]` jako funkcjonał
konfiguracji wielu Φ-źródeł, *zamiast* postulować ją lokalnie jako f(ψ).
Mechanizm musi:
- (a) Zachować S05 (jedno fundamentalne pole skalarne Φ z Z₂).
- (b) Zachować §5.1 (NIE BD/Horndeski — g_eff jako funkcjonał, nie zmienna dynamiczna).
- (c) Realizować §5.2 (kolektywny efekt fluktuacji jednego pola).
- (d) Działać dla pędu (§6 Lenz) i grawitacji jednocześnie.
- (e) Dawać alternatywę do `β_ppE^TGP = -15/4` poza zakresem GWTC-3 falsyfikacji.

## 8/8 gate criteria

### Criterion 1: Single fundamental claim
✅ **Claim:** g_eff = G[{Φ_i}, σ_ab, Φ̄] jest **derivable** z wielociałowej
dynamiki Φ, nie postulowana per-poziom.

### Criterion 2: Anchor inputs identified
✅ **External anchors (PDG/CODATA/observational):**
- GWTC-3 ~90 BBH posterior, |δφ̂_4| ≲ 0.18 (1σ)
- Solar system PPN: γ = 1 + (2.1 ± 2.3)·10⁻⁵, β = 1 + (1.2 ± 1.1)·10⁻⁴
- Mercury perihelion, Cassini Shapiro delay
- LIGO scalar mode bound (< few %)
- GW170817 c_GW = c (do 10⁻¹⁵)

✅ **TGP structural anchors (LOCKED):**
- S05 single-Φ axiom (FOUNDATIONS §1)
- §5.1 no-BD/Horndeski demarcation
- §6 Lenz back-reaction picture
- σ_ab gradient strain composite (OP-7 T2 2026-04-25)
- Φ-EOM kanoniczna z FOUNDATIONS §3
- V(Φ) M9.1'' canonical multi-vacuum (Phase 1.6 op-Phi-vacuum-scale 15/15)
- Dynamic-equilibrium framework (SPIN-SU2 N16, sympy 3/3, cycle 47/47)

### Criterion 3: Output specification
✅ **Cycle outputs:**
- (O1) Formal definition g_eff = G[{Φ_i}] z wyprowadzeniem
- (O2) 1PN limit γ=β=1 z derivation (Phase 2)
- (O3) 2.5PN β_ppE^new alternative do -15/4 (Phase 3)
- (O4) GWTC-3 falsifier check (Phase 4)
- (O5) M9.2 Lenz sympy (Phase 5)
- (O6) Cross-consistency z SU(2) paths (Phase 6)

### Criterion 4: Falsifiability
✅ **Hard tests:**
- (F1) γ_PPN = β_PPN = 1 musi przetrwać 1PN limit (solar system)
- (F2) β_ppE^new musi mieścić się w |δφ̂_4| ≲ 0.18 1σ (GWTC-3)
- (F3) m_inertial = m_grav (zasada równoważności) musi być automatyczna z S05
- (F4) c_GW = c (GW170817) — emergentna metryka NIE może wprowadzić scalar mode niezgodnego z bound

Każdy z F1-F4 jest binarny gate. Failure → STRUCTURAL_NO_GO honestly.

### Criterion 5: Anti-pattern checks (M03)
✅ **Patterns under watch:**
- Multi-candidate fit: NIE — jeden mechanizm (interakcja gradientów)
- Constructed criterion: NIE — falsifiers F1-F4 ustalone *przed* derivation
- Drift hardening: NIE — jeśli β_ppE^new fail F2, akceptujemy STRUCTURAL_NO_GO
- Anchor borrowed from external: σ_ab i V(Φ) są TGP-natywne (poziom 0/1)
- Algebraic re-arrangement: weryfikujemy, że gradient cross-terms są
  *strukturalnie* nowe, NIE algebraicznym przepisaniem M9.1''
- Definitional tautology: g_eff = G[{Φ_i}] jest *funkcjonalna* zależność,
  formalna struktura G musi być wyprowadzona z akcji, nie ad hoc

### Criterion 6: Tautology check pre-emptive
✅ **Risk identified:** R2 z README — w 1-source limicie emergentna f(ψ)
może zwinąć się do (4-3ψ)/ψ trywialnie. **Counter:** Phase 1-2 musi
explicite zidentyfikować *które* terms w 2-source przypadku są nieobecne
w 1-source limicie. Te terms są *non-tautologicznym* deliverable.

### Criterion 7: Cross-cycle consistency
✅ **Required cross-checks:**
- Z SPIN-SU2 (closed, 47/47): mechanizm interakcji generujący SU(2) (3 paths)
  musi być zgodny z mechanizmem generującym g_eff. Jeśli nie — flagged R4.
- Z op-Phi-vacuum-scale (closed): multi-vacuum landscape (ψ ∈ {0, 2/3, 4/3})
  determinuje Φ̄, niniejszy cykl spójny z tym landscape.
- Z S07 audyt (OPEN): niniejszy cykl jest *kandydatem* na zamknięcie
  postulate-status M9.1''.

### Criterion 8: Cycle termination criteria
✅ **Termination triggers:**
- Phase 4 fail F2 → STRUCTURAL_NO_GO (honest)
- Phase 6 cross-consistency fail → STRUCTURAL CONDITIONAL (cykl zamknięty
  z explicit acknowledgment cross-poziom incompatibility)
- Phase 5 M9.2 sympy fail → cykl downgraded do "metric-only", pęd osobny cykl
- Wszystko PASS → STRUCTURAL DERIVED, foundations §3 update z M9.1''
  postulate replacement

## Anchor inventory (detailed)

### A. TGP poziom 0 (substrat)
- `H_Γ = (V, E)`: GL-bond Hamiltonian v2 (2026-04-24)
- `Φ = ⟨ŝ²⟩`: skalar field z coarse-graining
- `K_ab = ⟨(∂_aŝ)(∂_bŝ)⟩`: gradient strain tensor
- `σ_ab = K_ab − (1/3)δ_ab Tr(K)`: traceless part (OP-7 T2 2026-04-25)
  — *embrion* tensorowej struktury, **DOTĄD nieaktywowany jako źródło g_eff**

### B. TGP poziom 1 (Φ-EOM)
- `D_kin[Φ] = ∇²Φ + 2(∇Φ)²/Φ = (1/3φ²)∇²(φ³)`: kanoniczny operator kinetyczny
- `V(φ) = (β/3)φ³ − (γ/4)φ⁴, β=γ`: potencjał (warunek próżni)
- Φ-EOM: `∇²Φ + 2(∇Φ)²/Φ + βΦ²/Φ_0 − γΦ³/Φ_0² = −qΦ_0 ρ`
- α=2 selection, K(φ) = K_geo φ⁴

### C. TGP poziom 2 (g_eff) — **REDEFINIOWANY w niniejszym cyklu**
- *Postulat dotychczas:* `g_tt = -c²(4-3ψ)/ψ, g_ii = ψ/(4-3ψ)·δ_ij`
  (M9.1'' canonical, **5σ FALSIFIED-OBSERVATIONAL** GWTC-3 2026-05-09)
- *Cel cyklu:* `g_eff^μν = G[{Φ_i}, σ_ab, Φ̄, x]` jako funkcjonał

### D. TGP poziom 3 (materia)
- L_mat = -(q/Φ_0)φρ (sprzężenie minimalne)
- ax:metric-coupling: materia sprzęga z Φ tylko przez g_eff
- SU(2) z 3 ścieżek (path A bifurkacja, path B horizon multipole, path C
  external embedding) — *cross-consistency hook*

### E. Falsified specific forms (avoid)
- M9.1 g_tt = -c²/ψ (potęgowa) — sfalsyfikowana 2026-04-25 (β_PPN=4)
- M9.1'' g_tt = -c²(4-3ψ)/ψ (hyperbolic) — sfalsyfikowana 5σ GWTC-3
  2026-05-09 (β_TGP=-15/4 vs |δφ̂_4| ≲ 0.18)

## Claims do verification (Phase 1+)

| Claim | Phase | Sympy required |
|---|---|---|
| C1: g_eff redukuje do flat-Minkowski w {Φ_i}→∅ limit (demarcation R1) | 1 | yes |
| C2: 1PN limit γ=β=1 | 2 | yes |
| C3: 2-source case zawiera gradient cross-terms ∂_μΦ_i·∂_νΦ_j strukturalnie | 3 | yes |
| C4: β_ppE^new różny od -15/4 (R2 counter) | 3 | yes |
| C5: |β_ppE^new| ≲ 0.18 (GWTC-3 window, F2) | 4 | numerical |
| C6: m_inertial = współczynnik back-reakcji, sympy obliczalny | 5 | yes |
| C7: m_inertial = m_grav automatic (S05 consequence) | 5 | yes |
| C8: c_GW = c (no scalar mode pathology) | 4-5 | yes |
| C9: Cross-consistency z SU(2) path D mechanism | 6 | structural |

## Cross-references

- [[./README.md]] — overview
- [[../../TGP_FOUNDATIONS.md]] §1, §3, §5.1, §5.2, §6, §8
- [[../op-SPIN-SU2-substrate-derivation-2026-05-08/Dynamic_equilibrium_framework.md]]
- [[../op-SPIN-SU2-substrate-derivation-2026-05-08/Internal_external_geometry_proposal.md]]
- [[../op-SPIN-SU2-substrate-derivation-2026-05-08/Phase6_absolute_binding.md]]
- [[../op-ppE-mapping/Phase1.5_G_SPA_lock.md]] (G_SPA=48 sympy-exact)
- [[../op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]
- [[../../audyt/S07_M911_derivation/]]
- [[../../audyt/T01_LIGO3G_falsifier/]]

## Status

🟢 **8/8 gate PASS** — cycle authorized to proceed Phase 1.

Next: `NEEDS.md` enumeration N1..Nk, then Phase 1 (formal definition g_eff).
