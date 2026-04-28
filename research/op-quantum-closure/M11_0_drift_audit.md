---
title: "M11.0 — Drift audit: quantum/RG drafts inventory"
date: 2026-04-26
cycle: M11.0
phase: drift_audit
status: SCAFFOLDING
parent: "[[M11_program.md]]"
related:
  - "[[../continuum_limit/cg_results.txt]]"
  - "[[../op1-op2-op4/README.md]]"
  - "[[../closure_2026-04-26/KNOWN_ISSUES.md]]"
tags:
  - TGP
  - M11
  - drift-audit
  - quantization
  - RG
---

# M11.0 — Drift audit: quantum/RG drafts inventory

> **Cel sub-cyklu:** inwentaryzacja istniejących skryptów + notatek kwantowych/RG w vault, klasyfikacja status (closure-grade / draft / superseded / irrelevant), i decyzja audit-vs-rebuild dla każdego draftu w kontekście M11 (kwantyzacja Φ, 1-loop, RG flow γ(k)).

---

## 0. Naming conventions disambiguation

**WAŻNE:** w tym vault istnieją **dwie niezależne konwencje numerowania M**:

| Numbering | Folder | Scope | Status |
|-----------|--------|-------|:------:|
| **OP-1 M1...M8** | [[../op1-op2-op4/]] | "Continuum limit theorem" — bond-fluctuation derivations supporting CG-1/2/3/4 | various (M2a/M2b/M3-M8 mostly DRAFT) |
| **Cycle M9, M10, M11** | [[../op-newton-momentum/]], [[../op-cosmology-closure/]], [[../op-quantum-closure/]] | High-level research cycles successive after closure_2026-04-26 | M9 ✅, M10 ✅, M11 ⏳ |

**M11 (this cycle) uses cycle numbering.** OP-1 M1-M8 (op1-op2-op4 folder) jest **input** do M11 (jako audit targets), NIE część M11 numbering.

---

## 1. Inventory — quantum/RG drafts in vault

### 1.1 Closure-grade (already closed; M11 USES as input)

| Plik | Folder | Status | Wynik | Rola w M11 |
|------|--------|:------:|-------|------------|
| [[../continuum_limit/cg_strong_numerical.py]] | continuum_limit | **CLOSURE-GRADE** | CG-2 8/8 PASS, η = 0.044 | **Input**: η z LPA' fundament dla M11.2-3 |
| [[../continuum_limit/cg_results.txt]] | continuum_limit | **CLOSURE-GRADE** | Documented η = 0.044 | Source-of-truth |
| [[../../tooling/scripts/tgp_erg_eta_lpa_prime.py]] | tooling/scripts | **CLOSURE-GRADE** | LPA' Wetterich solver, η* = 0.044 | **Input**: M11.2 reproduces self-consistently |
| [[../op1-op2-op4/M3_results.md]] (NPRG LPA polynomial truncation) | op1-op2-op4 | **CLOSURE-GRADE** | ν_WF ≈ 0.649 stable across N=4-10 | Critical exponent reference (3D Ising) |
| [[../op1-op2-op4/nprg_lpa_3d.py]] | op1-op2-op4 | **CLOSURE-GRADE** | Polynomial truncation working | Cross-check tool |

### 1.2 Drafts (need audit in M11)

| Plik | Folder | Status | Treść | Decyzja M11 |
|------|--------|:------:|-------|-------------|
| [[../op1-op2-op4/M2b_loop_derivation.md]] | op1-op2-op4 | **DRAFT** | 1-loop V_eff(Φ) bond-fluctuation; analytical framework | **M11.1 audit**: numerical closure |
| [[../op1-op2-op4/m2b_loop.py]] | op1-op2-op4 | **DRAFT** | 1-loop momentum integrals I_n; Taylor coefficients δc_2, δc_3, δc_4 | **M11.1 audit**: verify against canonical sek08a |
| [[../op1-op2-op4/M2b_results.md]] | op1-op2-op4 | **DRAFT** | Numerical 1-loop results from m2b_loop.py | **M11.1 audit**: drift check vs sek08a action form |
| [[../op1-op2-op4/M2a_HS_derivation.md]] | op1-op2-op4 | **DRAFT** (predecessor to M2b) | Hubbard-Stratonovich derivation; Φ = ŝ² | M11.1 input — not separate audit |
| [[../op1-op2-op4/M2a_sanity_check_results.md]] | op1-op2-op4 | **DRAFT** | Sanity checks for M2a tree-level | M11.1 input — not separate audit |
| [[../op1-op2-op4/M8_results.md]] | op1-op2-op4 | **DRAFT** | Migdal-Kadanoff lattice RG | **PERIPHERAL**: lattice context, not continuum Φ — defer |

### 1.3 Peripheral (not Φ-quantization; outside M11 scope)

| Plik | Folder | Status | Treść | Decyzja M11 |
|------|--------|:------:|-------|-------------|
| [[../muon_g_minus_2/ps01_schwinger_TGP_1loop.py]] | muon_g_minus_2 | **DRAFT** | 1-loop QED Schwinger anomaly | **NOT M11 SCOPE**: SM particle 1-loop, not Φ scalar |
| [[../uv_completion/]] | uv_completion | **DEFERRED** | UV completion bonus track (α_i running) | **NOT M11 SCOPE**: SM gauge running, post-M11 |

### 1.4 OP-1 M1-M8 chain (continuum limit theorem program)

OP-1 M1-M8 jest osobnym programem (continuum limit theorem). M11 audit TYLKO M2a/M2b (1-loop V_eff); pozostałe M-pieces (M3-M8) są **input dla CG-1/2/3/4** i pozostają w op1-op2-op4 namespace:

| OP-1 piece | Treść | Relacja do M11 |
|------------|-------|----------------|
| M1 (potential inventory) | Catalog v2 substrate potentials | n/a |
| M2 (sketch) | Outline V_eff derivation | n/a |
| M2a (HS derivation) | Hubbard-Stratonovich; Φ=ŝ² | **M11.1 prerequisite** |
| M2b (1-loop) | Bond-fluctuation 1-loop V_eff | **M11.1 audit target** |
| M2c (irrelevance) | Higher operator irrelevance | M11.1 input (RG argument) |
| M3 (Cb,Cg numerical) | Numerical critical exponents | CG-2 input (closure-grade ν) |
| M4 (φ variable) | Φ=ŝ² composite field analysis | M11.1 input |
| M5 (Zφ) | Field renormalization | **M11.2 prerequisite** (η connection) |
| M6 (gl-bond) | g+L bond analysis | OP-1 M2 input (continuum limit theorem) |
| M7 (gl eigenvalue) | Eigenvalue of g+L bond | OP-1 M2 input |
| M8 (Migdal-Kadanoff) | Lattice block RG | **PERIPHERAL** — not Φ continuum |

---

## 2. KNOWN_ISSUES items adresowane przez M11

Z [[../closure_2026-04-26/KNOWN_ISSUES.md]]:

| ID | Item | Pre-M11 status | M11 sub-cycle | Hipoteza |
|---:|------|:---:|:---:|---|
| **C.3** | γ = M_Pl² z RG flow Migdala-Kadanoffa applied to H_Γ | BLOCKED by OP-1 M2 | M11.4 | Achievable z CG-2 closed + numerical FRG flow do IR |
| **B.5** | g̃ = 1 (T-Λ; alternative wartości wymagają RG) | POSTULATE (g̃=0.98 calibrated) | M11.4 | Numerical FRG flow → IR g̃ value with O(1) accuracy |
| **B.3** | α₀ ≈ 4 (T-α calibration) | CALIBRATED, not derived | M11.4 | Plausibly z RG running α(ψ) → α₀ at IR fixed-point |
| **B.2** | n = 2 (T-α smoothness + WEP) | STRUCTURAL POSTULATE | M11.4 | Possibly z RG; otherwise retains POSTULATE status |

---

## 3. Drift classification (per M2 audit target)

### 3.1 m2b_loop.py / M2b_loop_derivation.md

**Substrate model used:** v2 substrate Hamiltonian
```
H_Γ[Φ] = Σ_i V_onsite(Φ_i) + J Σ_⟨ij⟩ Φ_i Φ_j (Φ_j − Φ_i)²
V_onsite(Φ) = (m₀²/2) Φ + (λ₀/4) Φ² + (T/2) ln Φ   (last term = HS Jacobian)
```

**Sek08a canonical action:**
```
S = ∫ d⁴x √(-g_eff) [ ½K(φ)g_eff^μν ∂_μφ ∂_νφ - V(φ) - (q/Φ_0)φρ ]
K(φ) = K_geo·φ⁴
V(φ) = (β/3)φ³ - (γ/4)φ⁴
```

**Drift identification:**
- m2b uses **discrete lattice** Hamiltonian → continuum action via blocking/RG (CG-1/2/3/4)
- m2b `V_onsite` ≠ sek08a `V` directly: m2b has degree 2 on-site (m₀² Φ + λ₀ Φ²) plus log-Jacobian; sek08a has degree 3+4 (β/3 Φ³ - γ/4 Φ⁴)
- After RG flow (continuum limit), **expected** matching: c_2 → 0 at IR (β=γ vacuum); c_3 → β/3; c_4 → -γ/4
- m2b output `δc_2, δc_3, δc_4` (one-loop shifts) — should be measured against canonical sek08a values

**Drift verdict:** **YELLOW** — m2b is correct lattice 1-loop; need M11.1 to verify continuum limit matches sek08a.

### 3.2 NPRG LPA polynomial truncation

**Status:** CLOSURE-GRADE (M3 results.md, ν_WF ≈ 0.649)

**Drift verdict:** **GREEN** — uses standard FRG framework (Wetterich-Litim), polynomial truncation N=4-10 stable. Already closure-grade; M11.2 cross-checks.

### 3.3 LPA' (η = 0.044)

**Status:** CLOSURE-GRADE (CG-2 8/8 PASS)

**Drift verdict:** **GREEN** — self-consistent η solver, Wetterich-Litim, d=3, n=1. M11.2 reproduces; no drift.

### 3.4 Migdal-Kadanoff (M8)

**Drift verdict:** **PERIPHERAL** — lattice block RG ≠ continuum Φ-EOM. Different methodology. Defer to post-M11 (M12+ topology).

### 3.5 Schwinger 1-loop (ps01)

**Drift verdict:** **OUT OF SCOPE** — QED 1-loop for muon g-2, not Φ scalar. Defer to UV completion / SM extension cycle.

---

## 4. Audit/rebuild decisions

| Draft | Decyzja | Sub-cycle | Reason |
|-------|---------|:---------:|--------|
| m2b_loop.py | **AUDIT** | M11.1 | Verify continuum limit matches sek08a; quantify drift |
| M2b_loop_derivation.md | **AUDIT** | M11.1 | Verify analytical framework consistent z sek08a |
| nprg_lpa_3d.py | **CROSS-CHECK ONLY** | M11.2 | Already closure-grade; use as cross-check |
| tgp_erg_eta_lpa_prime.py | **CROSS-CHECK ONLY** | M11.2 | Already closure-grade; reproduce from FRG flow |
| M8_results.md (Migdal-Kadanoff) | **DEFER** | post-M11 | Different context; M12+ topology |
| ps01 Schwinger | **OUT OF SCOPE** | n/a | QED, not Φ |

---

## 5. M11 sub-cycle plan (refined post-audit)

| Sub-cycle | Scope | Audit targets | Pre-test hypothesis |
|---:|-------|---|---|
| **M11.1** | 1-loop V_eff(Φ) closure | m2b_loop.py, M2b_loop_derivation.md, M2b_results.md | YELLOW → GREEN: m2b lattice 1-loop matches sek08a continuum limit |
| **M11.2** | β-functions {β, γ, K_geo} z FRG | tgp_erg_eta_lpa_prime.py, nprg_lpa_3d.py (cross-check) | Reproduce η = 0.044; β-functions exhibit Wilson-Fisher FP |
| **M11.3** | γ(k) k-dependent running | M11.2 output | M10.5.4 prediction: γ(k_LSS)/γ(k_CMB) = 1.244 reproduced precisely |
| **M11.4** | RG-driven structural derivation | C.3 + B.3 + B.5 + B.2 | Achievable: γ=M_Pl²·g̃ from FRG IR (resolves C.3); other possibly |
| **M11.R** | Synthesis | aggregate M11.0-4 | 6/6 PASS R-checks parallel to M10.R |

---

## 6. Open questions (research-level)

1. **Czy M2b's lattice V_onsite uniquely → sek08a continuum V?** Wymaga CG-1/2/3/4 closure (CG-2 done; CG-1/3/4 open ≥6 months math).
2. **Czy LPA' η = 0.044 jest dokładnie 3D Ising η?** η_Ising = 0.0364 (literatura) vs CG-2 LPA' η = 0.044 — ~20% LPA' approximation error; consistent class but not identical.
3. **Czy γ = M_Pl² · g̃ z RG flow?** Otwarte; M11.4 spróbuje numerycznie.
4. **Czy α₀ ≈ 4 ma natural origin z RG?** Plausible (O(1) IR FP value), but not yet derived.

---

## 7. M11.0 verdict

| Item | Verdict |
|------|:-------:|
| Inventory complete | ✅ |
| Naming convention disambiguated (OP-1 M1-8 vs cycle M9-11) | ✅ |
| Closure-grade items identified | ✅ (3 items: CG-2, NPRG LPA, LPA') |
| Drafts identified for M11 audit | ✅ (m2b primarily) |
| Audit/rebuild decisions made | ✅ (1 audit, 2 cross-check, 2 defer) |
| KNOWN_ISSUES items mapped to M11.4 | ✅ (C.3, B.3, B.5, B.2) |
| Pre-test hypotheses formulated | ✅ |

**M11.0 verdict: ✅ CLOSED (scaffolding-grade).** Ready for M11.1 (1-loop V_eff audit).

---

*M11.0 closed 2026-04-26. Inventory complete; M11.1-4 plan refined. Awaiting user direction for M11.1 launch.*
