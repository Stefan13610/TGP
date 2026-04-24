# M1 — Inventory of U(φ) across the TGP corpus

**Status:** Draft seeded from grep over `TGP_v1/**/*.tex` on 2026-04-24.
Extend as analysis proceeds.

---

## 1. Canonical form (target of the derivation)

    U(φ) = (β/3) φ³ − (γ/4) φ⁴,        φ := Φ / Φ₀,  Φ := ŝ²        (♦)

with Z₂ symmetry `ŝ → −ŝ` making Φ Z₂-even and φ dimensionless.
Vacuum condition: `U'(1) = 0  ⇒  β = γ`.

---

## 2. Corpus occurrences

### Axiom layer (dodatek B)

- `axioms/substrat/dodatekB_substrat.tex`
  - **L157** — postulated minimal form of U(φ) (Route 1 of
    `thm:beta-eq-gamma-triple`).
  - **L235** — Taylor expansion identification
    `β_eff = 3 c_3, γ_eff = −4 c_4`. This fixes the
    parametrization convention.
  - **L289** — substrate → TGP map for β:
    `β = (λ₀ / a_sub²) · C_β`  (`eq:beta-from-substrate`).
  - **L292** — substrate → TGP map for γ:
    `γ = (λ₀² / (v² a_sub²)) · C_γ`  (`eq:gamma-from-substrate`).
  - **L135-147** — MK-RG ratio at WF fixed point:
    `β/γ |_* = 0.575 · C_β / C_γ`. Vacuum condition β=γ requires
    `C_β / C_γ ≈ 1.74`.
  - **L149-206** — Theorem `thm:beta-eq-gamma-triple` with three
    routes (variational, Z₂-structural, MK-RG + MC).
  - **L306-308** — symbolic map `Z₂(ŝ → −ŝ) ↔ (β = γ)` noted as
    "automatic vacuum condition".

### Polish core ontology (sek01)

- `core/sek01_ontologia/sek01_ontologia.tex`
  - **L868** — `V_eff(Φ) = (β/3) Φ³ / Φ₀ − (γ/4) Φ⁴ / Φ₀²`
    (dimensionful form, on Φ rather than φ).
  - **L873** — cubic term identified as "static potential of
    Formulation B".

### Polish core dark energy (sek05)

- `core/sek05_ciemna_energia/sek05_ciemna_energia.tex`
  - **L178-179** — same potential written with `ψ` as the field
    variable (Formulation B). Note: this is a **different
    naming convention** — `ψ`, `φ`, `Φ` are used in different
    sections. Inventory flag: convention to unify.

### Polish core gauge/formalism (sek08)

- `core/sek08_formalizm/sek08_formalizm.tex`
  - **L1481-1485** — kinetic-driven combination
    `W(1) = 7β/3 − 2γ = γ/3 > 0` under β = γ. Structural link
    between static (β, γ) and full non-stationary coefficients
    `(7β/3, 2γ)`.
  - **L4123-4126** — clarification that the non-linear dynamical
    coefficients differ from the static ones by factors `7/3` and
    `2` (memo for M2: any `U_eff(φ)` derivation must reproduce
    this relation).

### Polish core N₀ derivation (sek10)

- `core/sek10_N0_wyprowadzenie/sek10_N0_wyprowadzenie.tex`
  - **L243** — `U(ψ) = (β/3)ψ³ − (γ/4)ψ⁴` used as the base
    potential in N₀ extraction. Consumer only; not a derivation.

### English summary (papers_external)

- `papers_external/tgp_english_summary/sec08_predictions.tex`
  - **L154** — English restatement of (♦).

### Core paper (canonical)

- `tgp-core-paper/paper/tgp_core.tex`
  - Declared as postulate + OP-1/OP-2/OP-4 listed as open.
  - v2 pivot (2026-04-24) changed the microscopic bond but left
    U(φ) unchanged at the postulate layer.

---

## 3. Parameter status after M1

| Symbol | Role | Status |
|---|---|---|
| `Φ` | composite field `ŝ²`, Z₂-even | derived (from axiom); concrete |
| `Φ₀` | vacuum expectation of Φ | derived: `Φ₀ ∝ v² / Φ_scale`; fixes dim. normalization |
| `φ = Φ / Φ₀` | dimensionless order parameter | derived |
| `β` | cubic coefficient of U | **OP-1 form** + `β = (λ₀ / a_sub²) · C_β`, C_β unknown |
| `γ` | quartic coefficient of U | **OP-4** + `γ = (λ₀² / (v² a_sub²)) · C_γ`, C_γ unknown |
| `β = γ` at vacuum | definitional (Route 1) | exact (from φ_vac = 1) |
| `β = γ` at all scales | **OP-2** | proven structurally (Route 2) in 3D Ising universality; numerically open (Route 3) |
| `g₀ᵉ = 0.86941` | TGP coupling calibration | empirical; equivalent to fixed γ |

---

## 4. Dependency graph

```
          H_Γ (v2)
             │
      ┌──────┴──────┐
      │             │
  (kinetic)     (potential)
      │             │
  F_kin^geo[φ]   U_eff(φ)
  = prop:        = POSTULATE (♦)
  substrate-        │
  action            │
  (closed)     ┌────┴──────────┐
               │               │
           OP-1 form        OP-4 value of γ
               │               │
               └───┬───────────┘
                   │
                OP-2 protection
                   │
            ├───────────────┤
            │               │
       Route 1 (vac)   Route 3 (WF)
       EXACT           C_β/C_γ = ?
```

---

## 5. Notational inconsistencies (low priority; flag for later unification)

- `dodatek B` uses `φ` (Φ / Φ₀, dimensionless).
- `sek01` uses `Φ` (dimensional).
- `sek05, sek10` use `ψ` (in Formulation B).
- `sec08_predictions` (English) uses `φ`.

M2 will work in dimensionless `φ` throughout for consistency with
the core paper and dodatek B. Translation table into `Φ`, `ψ`:
`φ = Φ / Φ₀ = ψ` (Formulation B normalizes `ψ_vac = 1`).

---

## 6. What M1 buys us for M2 and M3

1. **M2 target operator basis.** The M2 derivation must produce
   *at least* a cubic `φ³` and a quartic `φ⁴` term, and *at most*
   irrelevant higher-order corrections. Existing dodatek B map
   tells us M2 also has to reproduce:
   - `β = (λ₀ / a_sub²) · C_β`,
   - `γ = (λ₀² / (v² a_sub²)) · C_γ`,
   - the kinetic / dynamical combination
     `7β/3 − 2γ = W(1)` when β = γ.
2. **M3 target.** `C_β / C_γ` is the computable output.
   Prediction: `≈ 1.74` for β = γ.
3. **Entanglement map.** OP-1 = operator basis completeness.
   OP-2 = Route 3 closure (`C_β / C_γ → 1.74`) together with
   Route 2 irrelevance of sextic+. OP-4 = value of γ via
   `λ₀² / (v² a_sub²) · C_γ` once `C_γ` is computed in M3.

Effectively: M2 closes OP-1 + Route 2 of OP-2; M3 closes Route 3
of OP-2 and delivers the input for M4 to close OP-4.

---

## 7. Open questions surfaced by the inventory

- **Q1.** Do sextic `φ⁶` and higher terms have computable
  contributions in M2? If nonzero, check irrelevance against 3D
  Ising bootstrap (`Δ_{ε²} ≈ 3.83`).
- **Q2.** Does the kinetic/dynamical relation
  `(7β/3, 2γ)` come out of M2's effective action automatically,
  or does it rely on an independent normalisation choice in
  `sek08`?
- **Q3.** Is `a_sub` (lattice spacing) the same object in the
  microscopic map of dodatek B and in the block-RG step of M3?
  A mismatch here would spoil M4 calibration.
- **Q4.** `v² = (Jz − m₀²)/λ₀` (dodatek B L300). Is this still
  valid for the v2 GL bond or only for the v1 bilinear bond? The
  GL bond has no single `J` prefactor acting on `ŝ_i ŝ_j` — the
  correct kinematic identification needs re-checking against v2.
  **Priority item for M2 kickoff.**

---

## 8. Status

M1 inventory complete enough to unblock M2 and M3. Pending
full-text cross-check and graph visualisation (optional).
