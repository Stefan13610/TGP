# TGP Redirect Program — Post-Galaxy-Closure Research Agenda

**Date:** 2026-04-19
**Context:** After [[CLOSURE_2026-04-19.md|galaxy-scaling program closure]] (Path C), TGP predictive focus redirects to domains where the theory retains **unique, derived, falsifiable** predictions.

## Selection criteria

A research program is **HIGH PRIORITY** if:

1. **TGP-unique**: prediction differs from ΛCDM and MOND in a way attributable to TGP axioms (not shared phenomenology).
2. **Falsifiable**: there exists an observational test with timeline < 5 years.
3. **Derived, not fit**: prediction follows from existing sek02-sek10 formalism, not from post-hoc parameterization.
4. **Orthogonal to galaxy-MOND**: does NOT rely on the MOND-emergence claim retracted in gs65/gs66.

---

## Top-tier programs (begin immediately)

### P1. Cluster + Sterile Neutrino Hybrid

**Status:** Baseline established (gs55). Ready to extend.

**TGP-unique prediction:** Mass deficit in clusters is **32%** (gs54, gs55), vs MOND's **40%**. Therefore TGP requires **less** sterile neutrino mass to fill the gap than MOND does.

**Specific forecast (testable):**
- TGP+ν best-fit: `m_s ≈ 1.5-2.5 eV` (range from gs55)
- MOND+ν best-fit: `m_s ≈ 2.5-4 eV` (Sanders 2003)
- **Discriminant:** m_s < 2.5 eV favors TGP+ν over MOND+ν at cluster scale.
- **Galaxy-scale constraint:** v_thermal(m_s = 2 eV) >> v_esc(dwarf) → sterile ν must NOT cluster in galaxies. Retained.

**Observational test:**
- eROSITA DR1 X-ray + DES Y6 weak lensing on ~50 massive clusters (M₅₀₀ > 3×10¹⁴ M☉)
- KATRIN/Project 8 direct m_s bound
- Timeline: 2026-2027

**Research scripts needed:**
- `cluster_ng1_erosita_stacking.py` — stack 50-cluster profiles, fit (TGP ν(y)) + (NFW ν with m_s)
- `cluster_ng2_neutrino_constraint.py` — joint likelihood on m_s from stacked profiles + KATRIN prior
- `cluster_ng3_falsification.py` — quantify at what m_s value TGP+ν becomes worse than LCDM+CDM

**Falsifier:** If m_s > 3 eV required, or if no cluster deficit found after corrected analysis → TGP+ν hybrid rejected.

---

### P2. DESI w(z) Phantom-Crossing Falsification

**Status:** [[TGP/TGP_v1/research/desi_dark_energy/README.md|Placeholder exists]] with de1_tgp_w_prediction.py.

**TGP-unique prediction:** Dark-energy equation of state `w(z) ≥ −1` for all z, because sek05 potential `U(ψ) = (β/3)ψ³ − (γ/4)ψ⁴` has U'(1)=0 and bounded positive energy density. **No phantom crossing is possible in TGP.**

**Specific forecast (testable):**
- TGP: `w(z) = −1 + O(ψ̇²/ρ_DE)` ≥ −1, approaches −1 from above.
- CPL (LCDM+ext): w = w₀ + w_a·z/(1+z), can cross −1.
- DESI BAO-only result (2024-2025): `w₀ = −0.73 ± 0.11, w_a = −1.08 ± 0.33` → crosses −1 around z ~ 0.4 (3.9σ)
- **Discriminant:** Confirmed phantom crossing in DESI DR2/DR3 falsifies TGP. No crossing (DESI converges to w=−1) supports TGP.

**Observational test:**
- DESI DR2 (2026 Q2-Q3): BAO + RSD on full footprint
- DESI DR3 (2027): combined with Planck, DES, CMB-S4
- Timeline: 2026

**Research scripts needed:**
- `de2_tgp_frw_evolution.py` — solve sek05 Friedmann equations, compute w(z) from 0 to 2
- `de3_desi_comparison.py` — fit DESI DR1 data with TGP w(z) template, compare to CPL
- `de4_falsification_metric.py` — what σ-level phantom crossing in DESI DR2 falsifies TGP?

**Falsifier:** DESI DR2/DR3 confirms `w(z) < −1` at > 3σ in any redshift bin → TGP falsified on cosmology.

---

### P3. Black Hole Shadow / Ringdown

**Status:** Papers drafted ([[tgp_bh_shadow.tex]], [[dodatekC_ringdown.tex]]).

**TGP-unique prediction:** Shadow radius and ringdown QNM frequencies depend on substrate metric h(Φ)=Φ ansatz, which differs from pure Kerr at strong field. Specific predicted offsets at 0.1-1% level for M87* and Sgr A*.

**Observational test:**
- EHT high-cadence imaging (2025+): shadow radius at 1% precision
- LIGO O4/O5 (2025-2029): ringdown modes from BBH mergers, especially spectroscopy of (2,2,0), (2,2,1), (3,3,0) modes
- Timeline: ongoing, enhanced data arriving 2026-2029

**Research scripts needed:**
- Read + synthesize existing tex: extract quantitative predictions.
- `bh1_shadow_radius_numerical.py` — compute TGP shadow for Kerr+substrate, compare to GR
- `bh2_qnm_extraction.py` — solve perturbation equation in TGP metric, compute (ℓ,m,n) frequencies
- `bh3_eht_comparison.py` — compare to M87* / Sgr A* EHT data with error bars

**Falsifier:** EHT shadow radius within 0.5% of pure Kerr prediction → TGP deviation ruled out. Likewise for ringdown modes.

---

### P4. Particle Sector Closure

**Status:** Foundational results in [[TGP/TGP_v1/research/mass_scaling_k4/]], [[TGP/TGP_v1/research/brannen_sqrt2/]], [[TGP/TGP_v1/research/cabibbo_correction/]], [[TGP/TGP_v1/research/why_n3/]]. Three specific OPEN problems identified in STATUS §VIII.

**Open problems (purely analytical, no observational dependency):**

1. **α₃ closed form**: Current value 0.089722... (30 digits). NOT π²/110 (PSLQ fails). Needs new mathematical structure.
2. **g₀^τ determination**: Currently fit to Koide K=2/3. No first-principles derivation. Requires ODE analysis in Brannen setup.
3. **Quark Koide / QCD running**: K_up=0.85, K_down=0.73 currently fail K=2/3. Test: does 2-loop QCD running reconcile?
4. **Neutrino Koide**: K_ν ~ 0.58 < 2/3. Test: Majorana mass contribution?

**Observational test:** Not observational; purely mathematical. But if solved:
- Quark Koide reconciliation → testable at ~0.1% level via PDG masses + lattice QCD.
- Neutrino Koide → testable via KATRIN, DUNE, JUNO mass hierarchy.

**Research scripts needed:**
- `ps1_alpha3_hunting.py` — extend PSLQ to 60+ digits, search for Γ-function / elliptic integral identities
- `ps2_g0_tau_ODE_scan.py` — full Brannen ODE parameter sweep with bisection on τ mass
- `ps3_quark_koide_qcd.py` — 2-loop running from m_Z to μ_effective, Koide evaluation
- `ps4_neutrino_majorana.py` — include νR sector, check K ≥ 2/3 conditions

**Falsifier:** No explicit falsifier (these are structural closures). If α₃ cannot be closed → likely TGP-independent new constant. If g₀^τ cannot be derived → one free parameter admitted, weakens Koide achievement.

---

## Secondary programs (maintain, don't accelerate)

- **Hubble tension** (gs41 PASS at CMB, gs46 a₀(z)): maintain. TGP structurally cannot resolve H0 tension (gs20 verdict).
- **S8 tension** ([[TGP/TGP_v1/research/s8_tension/]]): maintain. TGP effect ~0.001%, 8500× too weak.
- **UV completion** (placeholder only): postpone unless needed for QG program.
- **QM continuum theorems** (CG-1, CG-3, CG-4): ongoing mathematical work, 6-12 month timeline.

---

## Explicitly abandoned

- **Galaxy-MOND first-principles derivation** (see gs65 + gs66 + CLOSURE). The ν(y) parameterization remains available for fitting but is NOT a TGP prediction.
- **Cluster deficit resolution via substrate mechanism** (gs52, gs54 ruled out). Only viable with sterile ν (gs55 path → P1).
- **Mech F substrate oscillation** (gs64 refuted).
- **Bridge (a) linear FRW** (gs66 no-go theorem).

---

## Timeline & milestones

| Quarter | Milestone |
|---|---|
| 2026 Q2 | CLUSTER: First 10 eROSITA-stacked cluster profiles fitted (P1). DESI: TGP w(z) template ready (P2). |
| 2026 Q3 | DESI DR2 released; phantom-crossing test performed (P2). BH: first draft of shadow-radius comparison (P3). |
| 2026 Q4 | α₃ literature/PSLQ deep search (P4). Cluster full-stack analysis (P1). |
| 2027 Q1 | g₀^τ ODE study completion (P4). Quark Koide + QCD running (P4). |
| 2027 Q2 | LIGO O4 ringdown spectroscopy data → compare to TGP QNM (P3). |
| 2027 Q3 | Euclid DR1 a₀(z) data → **optional**: if strong evidence for a₀∝H(z), reopen galaxy bridge (b)/(c). |
| 2027 Q4 | Integrated re-audit: does TGP remain falsifiable on at least 2 of 4 priority programs? |

---

## Publication strategy adjustment

The [[TGP_STATUS_2026-04-19.md]] §XI "Publication Strategy" must be updated:

- **Foundations paper**: unchanged. Core results (k=4, N=3, Koide, QM emergence) are intact.
- **Galaxy paper**: REWRITE. Present ν(y) phenomenology honestly — "TGP is *consistent with* MOND under parameterization" — not "TGP *derives* MOND". Cite gs65/gs66 explicitly as the formal result bounding the claim.
- **Cluster + ν paper**: NEW. Present P1 as primary cluster test with TGP's lower deficit than MOND.
- **Cosmology paper**: NEW. Present P2 as DESI falsification test.
- **BH paper**: PREPARE. Synthesize existing tex material + P3 numerics.
- **Particle sector paper**: CONSOLIDATE. Bundle lepton masses + Koide + Cabibbo + α₃ into one paper.

---

**Final remark:** Path C is a gain, not a loss. TGP keeps its strongest predictions (particle sector) and acquires cleaner falsifiers (P1, P2, P3) without the burden of defending an over-claim. The discipline of formal derivation applied to gs65/gs66 is now the standard for all future programs.
