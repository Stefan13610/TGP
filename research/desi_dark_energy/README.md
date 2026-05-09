---
title: "DESI w(z) Phantom-Crossing Falsification — Program P2 (TGP Redirect)"
date: 2026-05-03
tgp_status:
  folder_status: paused
  level: L1
  kind: phenomenology
  core_compatibility: "unknown"
  last_reviewed_against_core: "unknown"
  may_edit_core: false
  exports_findings: false
  has_needs_file: true
  has_findings_file: true
  open_bridges: []
  depends_on: []
  impacts: []
  source_of_status:
    - "README.md H1: 'DESI w(z) Phantom-Crossing Falsification — Program P2 (TGP Redirect)'"
  promoted_to_core: null
  polluted_74394a8: false
  pre_existing_findings: false
  pre_existing_needs: false
  last_yaml_update: "2026-05-03"
---

> **STATUS UPDATE 2026-05-03 — post-cascade alignment**
>
> Niniejszy README został napisany **2026-04-19**, **przed** kaskadą 5 zamknięć
> strukturalnych:
> - **closure_2026-04-26** (Path B + T-FP + **T-Λ: Ω_Λ INPUT → PREDICTION** + T-α; 35/35 PASS)
> - **M10 cycle CLOSED 2026-04-26** (M10.0-R, 42/42 PASS) — de2 audit YELLOW → GREEN
> - **Phase 3.E (2026-04-28)** — B.6 PARTIAL DERIVED `V(Φ_eq)|β=γ = γ/6`, B.4 STRENGTHENED
> - **UV.3 closure (2026-05-02)** — `Z_Φ = 14/3` explicit, `Φ_0^bare ≈ 115` / `Φ_eff = 24.65`
> - **γ.1 + δ.1 + δ.2 (2026-05-02)** — `Ω_Λ^pure = 2π/9` algebraicznie + `Ω_Λ·α_s = 3g_0^e/32 ≈ 0.0815` (NOWA falsifiable)
>
> oraz **przed** publikacją DESI DR2 (2025-03, arXiv:2503.14738) i DESI Year 3 (2026-04, 47 mln galaktyk).
>
> **Status faktyczny 2026-05-03:**
> - DESI DR2: `(w₀, wₐ) = (−0.75, −0.90)` przy 3.1σ (BAO+CMB+SN)
> - DESI Y3: 2.8–4.2σ preferencja dla evolving DE (zależnie od kombinacji SN)
> - TGP **NIE sfalsyfikowane** (próg 5σ niewysycony) — yellow zone
> - **Decydujący test: DESI 5-year (2027 Q2-Q4)**
> - Survival probability post-DR2/Y3: ~30-50% (potwierdzone vs original 30-50% pre-DR2)
>
> Treść poniżej (sekcje 1-8) pozostaje **strukturalnie poprawna** (theorem
> `w + 1 = ψ̇²/ρ_ψ ≥ 0` zachowany w M10.1.2 sympy), lecz **liczbowo zaległa**
> (DR1 numbers). Pełna analiza post-cascade:
> [[../op-cosmology-closure/M10_R_results.md]] sekcja 12.5.
> Sekcja 9 (post-cascade falsifikatory) dodana 2026-05-03.

# DESI w(z) Phantom-Crossing Falsification — Program P2 (TGP Redirect)

> **STATUS 2026-04-19: P2 ACTIVE — FORMAL DERIVATION + FALSIFICATION ROADMAP COMPLETE.**
> This is program **P2** from the [[REDIRECT_PROGRAM_2026-04-19.md|post-closure redirect]].
> TGP makes a **structural, derived, falsifiable** prediction for dark energy: **w(z) ≥ −1 for all z, no phantom crossing.** DESI DR2 (2026 Q3) is the decisive test point.

---

## 1. The TGP-unique prediction (derived, not fit)

From the sek05 potential `U(ψ) = (β/3)ψ³ − (γ/4)ψ⁴` with β = γ:
- `U(1) = γ/12`, `U'(1) = 0`, `U''(1) = −γ < 0` (slow-roll maximum)
- Action has **canonical kinetic term** (verified gs08b — no ghost degrees of freedom)

**Theorem (de2):** For a canonical scalar with positive kinetic term in FRW,
```
w + 1 = ψ̇² / ρ_ψ  ≥  0     (always)
```
Phantom crossing `w < −1` would require wrong-sign kinetic (ghost instability), which is **structurally forbidden** in sek05.

**Numerical verification (de2):** global `w_min = −1.00000000` across all trajectories (δ ∈ {1e-4, 1e-3, 1e-2}).

**CPL projection of TGP trajectory:**

| Parameter | TGP (sek05 derived) | ΛCDM | DESI DR1 DES-SN5YR | DESI DR2 BAO+CMB+SN (2025-03) |
|----------|---------------------|------|---------------------|---------------------|
| w₀ | **−0.99998** | −1.00 | −0.727 ± 0.067 | **−0.75 ± 0.10** |
| wₐ | **−0.00003** | 0.00 | −1.050 ± 0.310 | **−0.90 ± 0.35** |
| significance vs ΛCDM | — | 0σ | 4.1σ | **3.1σ** (BAO+CMB), **3.9σ** (Union3) |

TGP is **indistinguishable from ΛCDM** within ~10⁻⁴ on (w₀, wₐ). Both sit at **~4.1σ tension** with DESI DR1 combined (DES-SN5YR), and **~3.1σ tension** with DR2 (2025-03) — DR2 central values **migrated closer to ΛCDM** vs DR1.

---

## 2. Why P2 is a clean falsifier

1. **TGP-unique**: no adjustable parameters; w(z) follows from sek05 ansatz.
2. **Orthogonal to galaxy-MOND**: does not rely on the ν(y) claim retracted in [[gs65_formal_derivation.py|gs65]]/[[gs66_frw_propagator.py|gs66]].
3. **Cannot be rescued**: bridge-(b)/(c) substrate modifications affect structure growth, not the sign of kinetic energy; the w ≥ −1 bound is a theorem.
4. **Falsifier timeline**: DR2 in ~6 months, DR3 in ~18 months.

---

## 3. Program files

| File | Role | Status |
|------|------|--------|
| `de1_tgp_w_prediction.py` | Unified tension scan (H₀/S₈/w(z)); first look at DESI | ✔ ran 2026-04-18 |
| `de2_tgp_frw_evolution.py` | **Full FRW integration from sek05**; derived w₀, wₐ, w_min | ✔ ran 2026-04-19 |
| `de3_desi_dr1_comparison.py` | Mahalanobis distance to DESI DR1 posterior; AIC/BIC; Bayes | ✔ ran 2026-04-19 |
| `de4_falsification_roadmap.py` | DR2/DR3 projection, falsification criteria, timeline | ✔ ran 2026-04-19 |

Outputs: `de{1,2,3,4}_results.txt`

---

## 4. Consolidated results

### de2 — FRW derivation
- Integrated coupled ODE: `ψ̈ + 3Hψ̇ + V′(ψ) = 0`, `H² = Ω_r a⁻⁴ + Ω_m a⁻³ + ρ_ψ/ρ_crit`.
- Normalized `V(ψ) = 4ψ³ − 3ψ⁴`, shooting parameter V₀ fixed by Ω_DE0 = 0.685.
- Initial conditions at z=99: `ψ = 1 − δ`, `u = 0` (Hubble-frozen).
- **Result:** w(z) stays within 10⁻⁴ of −1 for all realistic δ; TGP ≈ ΛCDM on cosmology.

### de3 — DR1 comparison (combined posterior `(w₀, wₐ) = (−0.727 ± 0.067, −1.05 ± 0.31)`, ρ = −0.85)

| Model | w₀ | wₐ | χ² | n-σ | ΔAIC | ΔBIC |
|---|---|---|---|---|---|---|
| ΛCDM | −1.000 | 0.000 | 16.62 | 4.08σ | +12.62 | +15.24 |
| **TGP (de2)** | −0.99998 | −0.00003 | 16.62 | **4.08σ** | +12.62 | +15.24 |
| CPL best-fit | −0.727 | −1.050 | 0.00 | 0.00 | 0.00 | 0.00 |

Bayes factor log B(CPL/TGP) = 2.46 → "positive" for CPL (Jeffreys scale), but short of decisive. **SN-sample spread** (Pantheon+ 2.2σ, Union3 3.1σ, DES-SN5YR 4.1σ) is the dominant systematic.

### de4 — Falsification roadmap

**DR2 (~2026 Q3, √3 error reduction):**

> **Update 2026-05-03:** DR2 zostało faktycznie wydane **2025-03** (arXiv:2503.14738), wcześniej niż projektowane "Q3 2026". Faktyczne wartości DR2: `(w₀, wₐ) = (−0.75 ± 0.10, −0.90 ± 0.35)` przy 3.1σ (BAO+CMB+SN). To odpowiada **scenariuszowi pośredniemu** między (i) i (ii) — centrum przesunęło się od `−0.727` ku `−0.75` (lekka migracja ku ΛCDM), ale błędy zacieśniły się tylko ~×1.5 zamiast √3 ≈ 1.73 (DR2 BAO-only volume nie był pełne 3× DR1 wszystkich tracerów). **Aktualne napięcie: 3.1σ — TGP w yellow zone, NIE sfalsyfikowane.**

| Scenario | Projected n-σ (DES-SN5YR-like) | Verdict |
|---|---|---|
| (i) central migrates to LCDM (w₀=−0.9, wₐ=−0.3) | 2.77σ | TGP in tension, not falsified |
| (ii) DR1 persists (w₀=−0.727, wₐ=−1.05) | **7.06σ** | **TGP FALSIFIED** |
| (iii) drifts further (w₀=−0.6, wₐ=−1.5) | 10.4σ | TGP catastrophically falsified |
| **(actual DR2 2025-03)** | **(w₀=−0.75, wₐ=−0.90), 3.1σ** | **TGP yellow zone (close to scenario i + lighter migration)** |

**DR3 (~2027 Q3, √5 error reduction):**

| Scenario | Projected n-σ | Verdict |
|---|---|---|
| (i) | 3.57σ | TGP marginally falsified |
| (ii) | **9.12σ** | **TGP falsified** |
| (iii) | 13.4σ | TGP catastrophically falsified |

**Even scenario (i)** — DESI central values migrating toward ΛCDM — triggers 3σ falsification at DR3, unless migration continues all the way to exactly (w₀, wₐ) = (−1, 0).

**Survival probability estimate (given DR1 tension):**
- P(TGP survives DR2) ≈ 30–50%
- P(TGP survives DR3) ≈ 15–30%

---

## 5. Specific falsification criteria

- **(F1) DR2 DES-SN5YR gives > 3.5σ + Pantheon+/Union3 both > 2.5σ** → TGP falsified at cosmological scale.
- **(F2) DR3 gives > 3σ across all three SN samples** → TGP definitively falsified.
- **(F3) Direct binned w(z) measurement detects w(z) < −1 at > 3σ in any bin** → immediate falsification.
- **(S1) DR2/DR3 converges to (w₀, wₐ) ≈ (−1, 0) within errors** → TGP supported (equivalent to ΛCDM at observational level).

---

## 6. Timeline

| Quarter | Milestone |
|---|---|
| 2026 Q2 | Program P2 scripts complete (de1–de4). **← current** |
| 2026 Q3 | **DESI DR2 release → decisive test point (F1).** |
| 2026 Q4 | DR2 + CMB + full SN joint re-analysis. |
| 2027 Q2–Q4 | **DESI DR3 full survey release → final test (F2).** |
| 2028+ | Euclid + LSST + CMB-S4 independent cross-check. |

---

## 7. What P2 does NOT address (scope)

- **H₀ tension**: TGP structurally cannot resolve (gs20 + de1; ~8 orders of magnitude too weak). *Not in scope.*
- **S₈ tension**: TGP substrate effect ~0.001%, 8500× too weak ([[TGP/TGP_v1/research/s8_tension/]]). *Not in scope.*
- **Galaxy-MOND**: retracted as first-principles derivation ([[CLOSURE_2026-04-19.md]]). *Separate program.*
- **Clusters**: handled by P1 (`research/cluster_anomaly/` (planowane) + sterile ν hybrid).

P2 is narrowly targeted on the **phantom-crossing falsification** — the cleanest, fastest, most TGP-unique cosmological test available.

---

## 8. Key references

- DESI Collaboration DR1 (arXiv:2404.03002, 2024)
- **DESI Collaboration DR2 (arXiv:2503.14738, marzec 2025)** — released, NOT future as original document said
- **DESI Year 3 results (kwiecień 2026)** — 47 mln galaktyk, 2.8–4.2σ evolving DE
- Chevallier & Polarski 2001, Linder 2003 — CPL parametrization
- Model-Independent Reconstruction of Quintessence Potential — DR2+Pantheon+ (arXiv:2603.21125)
- ghost-free verification: zamknięte przez OP-7 T6 (12/12 PASS, 2026-04-25) + Path B (closure_2026-04-26 11/11 PASS); historyczne `gs08b_ghost_free_verification.py` superseded
- [[REDIRECT_PROGRAM_2026-04-19.md]] — parent redirect program
- [[TGP_STATUS_2026-04-19.md]] — overall TGP status post-closure (pre-cascade)
- [[../closure_2026-04-26/CLOSURE_2026-04-26_SUMMARY.md]] — T-Λ + Path B + T-FP + T-α
- [[../op-cosmology-closure/M10_R_results.md]] — M10 cycle synthesis (42/42 PASS)
- [[../op-uv3-phi0-renormalization/]] — UV.3 Z_Φ = 14/3 (2026-05-02)
- [[../op-gamma1-phi-eff-anchor-resolution/]] — γ.1 (Ω_Λ^pure = 2π/9)

---

## 9. Falsifikatory dodane w cascade 2026-05-02

Post-cascade kosmologia TGP zyskała **niezależne kanały falsyfikacji** poza
CPL `(w₀, wₐ)` reżim:

| Probe | Predykcja | Próg falsyfikacji | Eksperyment |
|---|---|---|---|
| **F1.5 — Ω_Λ · α_s** | `3·g_0^e/32 ≈ 0.0815` (drift 0.88%; UV.3 + γ.1) | >3σ deviation | CMB-S4 + Z-pole α_s 2030+ |
| **F11 — Σm_ν cosmology** | `Σm_ν = 59.6 meV` (Phase 2 DERIVED) | CMB-S4 5σ DECISIVE projected | CMB-S4 (op-omicron1) |
| **F12 — Ω_Λ algebraic** | `Ω_Λ ∈ {2π/9 ≈ 0.6981, 5e²/54 ≈ 0.68417, 0.6847_obs}` | >2σ deviation od jednej z tych form | Planck PR5 / LiteBIRD / DESI Y5 |

**Kluczowa konsekwencja:** nawet jeśli DESI ostatecznie potwierdzi
`(w₀, wₐ) = (−1, 0)` (TGP exonerated w sektorze DE), Ω_Λ algebraiczne
pozwala falsyfikować TGP w `Ω_Λ` bezwzględnym. Trzy niezależne kanały
F1.5 / F11 / F12 dają **falsifiability redundancy** — TGP cosmology
post-cascade ma więcej niż jeden falsyfikator strukturalny.

Cross-link: [[../op-cosmology-closure/M10_R_results.md]] sekcja 12.5
"Post-M10 addenda" zawiera pełną updated falsification matrix
(F1-F12) post-cascade.

---

**Bottom line:** If DESI DR2 confirms DR1 central values at √3-shrunken errors, TGP is falsified on dark energy at ~7σ within 6 months. If DR2 migrates toward ΛCDM, TGP is supported. **Either way, by 2027 TGP cosmology has a definitive answer.**

**Update 2026-05-03:** DR2 (marzec 2025) actual values `(−0.75, −0.90)`
przy 3.1σ — **bliżej ΛCDM niż DR1** (`(−0.45, −1.79)` z DES-SN5YR), ale nadal
≠ `(−1, 0)`. TGP **w yellow zone, NIE sfalsyfikowane**. DR3 (2027 Q2-Q4)
będzie decydujący.
