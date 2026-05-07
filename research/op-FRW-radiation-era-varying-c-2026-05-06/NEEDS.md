---
title: "NEEDS — op-FRW-radiation-era-varying-c (open questions)"
date: 2026-05-06
parent: "[[README.md]]"
type: needs
tgp_owner: research/op-FRW-radiation-era-varying-c-2026-05-06
tags:
  - needs
  - open-questions
  - EXT-1
---

# NEEDS — open questions per Phase

## Phase 1 prerequisites (gotowe)

- ✅ L01 EXECUTED 2026-05-04 (formal definition ρ ≡ -T^μ_μ/c_0²)
- ✅ ax:c-ax:G axiomy w sek04_stale.tex
- ✅ M9.1'' canonical form (z S07 caveat)
- ✅ Phase 0 balance sheet 8/8 ☑ PASS

## Phase 1 NEEDS (do rozstrzygnięcia w Phase 1)

### N1 — Reżim ψ << 1 dla M9.1''

**Pytanie:** Czy M9.1'' (zaprojektowana dla weak field wokół ψ=1) jest
**valid w erze radiacyjnej** (z > 10¹⁰), gdzie ψ_then potencjalnie << 1?

**Status:** OPEN, **R1 critical risk** (Phase0_balance §7).

**Akcje:**
- F1.4 asymptotic analysis Φ-EOM dla z >> z_eq
- Identify czy ψ(z=10¹⁰) ∈ [0.1, 10] (M9.1'' valid range) albo ψ → 0/∞
  (M9.1'' łamie założenia perturbacyjne)

**Connect:** [[../../audyt/S07_M911_derivation/README.md]] — M9.1''
status (E) postulat vs (P) propozycja.

### N2 — Φ-EOM w FRW solver stability

**Pytanie:** Czy numerical solver Φ-EOM dla z ∈ [10³, 10¹⁰] jest
stable (no runaway, no singularity bez physical reason)?

**Status:** OPEN, do rozstrzygnięcia w F1.4.

## Phase 2 NEEDS (do rozstrzygnięcia w Phase 2)

### N3 — Quantum EM trace anomaly (z L01 NEEDS N1)

**Pytanie:** Czy `<T^μ_μ>_QED` quantum anomaly (NIE klasyczna T^μ_μ_EM = 0)
ma znaczący wkład do ρ w erze radiacyjnej?

**Status:** OPEN z L01 NEEDS N1 (op-L01-rho-stress-energy-bridge-2026-05-04).

**EXT-1 promotion:** Z P3 do **P1** (per Phase0_balance §7) — w erze
radiacyjnej kwantowy EM kondensat może modyfikować ρ_total i tym samym
H_TGP(z=10⁹).

**Akcje:**
- Sympy compute `<T^μ_μ>_QED` w 1-loop QED
- Check magnitude vs ρ_matter w obszarze BBN
- Phase 2 sub-test F2.5 (σ_T·n_e·c kinematics)

### N4 — QCD vacuum kondensat (z L01 NEEDS N2)

**Pytanie:** Czy QCD vacuum kondensat ~ Λ_QCD⁴ (anomalia śladu) ma
znaczący wkład w erze T ≥ 200 MeV (z ≈ 10¹²)?

**Status:** OPEN z L01 NEEDS N2.

**EXT-1 promotion:** Z P3 do **P1** (Phase0_balance §7).

**Akcje:**
- Compute Λ_QCD⁴ contribution to ρ_total dla T ≥ T_QCD
- Phase 1 limit asymptotyczny (F1.4) musi to uwzględnić

### N5 — Varying ℏ effects on nuclear cross-sections

**Pytanie:** Jak varying ℏ(Φ) wpływa na BBN nuclear cross-sections
(Coulomb barrier, kinetics)?

**Status:** OPEN, krytyczne dla F2.3 (BBN abundance ⁴He).

**Akcje:**
- Re-derive Coulomb barrier w funkcji ℏ(Φ) = ℏ_0·√(Φ_0/Φ)
- σ_nuclear(ℏ, energy) computation
- Phase 2 F2.3, F2.4 (BBN abundance)

## Phase 3 NEEDS (do rozstrzygnięcia w Phase 3)

### N6 — CMB sound horizon w varying-c cosmology

**Pytanie:** Jak r_s(z_drag) zmienia się z varying c(Φ) podczas
rekombinacji?

**Status:** OPEN, krytyczne dla F3.1 (l_1 = 220).

**Akcje:**
- r_s = ∫c_s(z) dz / H(z) z varying c(Φ_drag)
- Compare with Planck 147.05 ± 0.30 Mpc

### N7 — N_eff w TGP (relativistic species count)

**Pytanie:** Co liczy się jako "relativistic species" w TGP cosmology?
Czy Φ-kinetic terms wkład jak relativistic DOF?

**Status:** OPEN, krytyczne dla F3.2 (N_eff = 3.046).

**Akcje:**
- Identify Φ-kinetic terms behaviour w erze radiacyjnej
- Compute N_eff_TGP vs Planck 2.99 ± 0.17

## Decision tree NEEDS (post-Phase 3)

### N8 — Decision criteria (A/B/C)

**Pytanie:** Czy wynik Phase 3 jest jednoznaczny (A DERIVED / B
STRUCTURAL / C FAIL) czy wymaga dyskusji autora?

**Status:** OPEN, do rozstrzygnięcia post-Phase-3.

**Akcje:**
- Aggregate 5-channel falsification verdict
- Multi-path consistency check
- Determine final classification

### N9 — Pivot decision (jeśli FAIL)

**Pytanie:** Jeśli Phase 3 = FAIL, czy autor decyduje na ścieżkę D
(L_mat extension, **narusza S04 zamknięty 2026-05-04**) czy ścieżkę E
(przyznanie zakresu post-recombination)?

**Status:** OPEN, **wymaga decyzji autora** (NIE Claudian).

**Implications:**
- Ścieżka D: re-open S04, dodatkowe sprzężenie dilatonowe dla pól
  cechowania
- Ścieżka E: TGP jako "GR w limicie post-recombination", drastycznie
  obniża rangę

## Cross-references

- [[Phase0_balance.md]] §7 risk assessment (R1, R2, R3)
- [[../../audyt/EXTERNAL_REVIEW_2026-05-06.md]] §EXT-1 v2
- [[../../audyt/L01_rho_operational/EXT1_FRW_radiation_era_2026-05-06.md]]
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/NEEDS.md]] — N1, N2 source
- [[../../audyt/S07_M911_derivation/README.md]] — M9.1'' status
