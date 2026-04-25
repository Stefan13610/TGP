# OP-M92 Phase 0+ — Candidate D cosmology cross-check (results)

**Data:** 2026-04-25
**Status:** **Phase 0+ cosmology cross-check DONE — verdict POSITIVE**
**Sub-test:** Candidate D vs OP-DESI w(z) ≥ −1 prediction (de2 theorem)
**Następnik:** Phase 1 deferred item (d) "Cosmological consequences" — partial closure

---

## TL;DR

> Candidate D back-reaction term `α T^μν J_μ J_ν` w homogeneous FRW efektywnie
> przesuwa Φ kinetic coefficient o czynnik `(1 − 8πG α ρ_m)`. Dla α
> skalibrowanego do Sgr A* photon ring (`α_SI ≈ 180 s²`) ratio
> `8πG α ρ_m` wynosi:
>
> - **Dziś (z=0)**: 8.08e−34 (margines 1.24e+33× poniżej phantom threshold)
> - **Recombination (z=1090)**: 1.05e−24 (margines 9.53e+23×)
> - **BBN (z=4×10⁸)**: 5.17e−08 (margines 1.93e+07×)
> - **Phantom transition**: z ≈ 1.07×10¹¹ (pre-BBN, unobservable)
>
> **Verdict POSITIVE:** Candidate D **strukturalnie zachowuje** OP-DESI
> prediction `w(z) ≥ −1` przez całą observable cosmology. de2 theorem
> (`w_psi + 1 = ψ̇²/ρ_ψ ≥ 0`) pozostaje valid pod M9.2-D pivot.
>
> **Cross-program implication:** OP-DESI pozostaje **niezależnym** testem
> TGP. Phantom crossing observation (jeśli DR3 confirms) would falsify
> BOTH M9.1'' AND M9.2-D — nie jest redundant test.

---

## 1. Setup — Candidate D w FRW

W jednorodnej FRW z `Φ = Φ(t)`:
- `J_μ = ∂_μ Φ = (Φ̇, 0, 0, 0)`
- `T^μν = diag(ρ_m, p_m, p_m, p_m)` (perfect fluid)
- `T^μν J_μ J_ν = T^00 Φ̇² = ρ_m Φ̇²`

Modyfikacja Φ kinetic Lagrangian:

```
L_kin → −[1/(8πG)] Φ̇² + α ρ_m Φ̇² = −[1/(8πG) − α ρ_m] Φ̇²
```

⇒ effective kinetic coefficient `K_eff = −[1/(8πG)](1 − 8πG α ρ_m)`.

Phantom crossing wymaga `K_eff > 0`, czyli `8πG α ρ_m > 1`.

## 2. Calibration — α z Sgr A* photon ring

Z OP-M92 Phase 0+ sketch: α ≈ 0.1 (geometric units, M=1) tunuje scenario
(e) target przy r_ph^TGP=3.88M, ψ=1.168.

Dla Sgr A* (M ≈ 4.297×10⁶ M_⊙):
- R_S = 2GM/c² ≈ 1.27×10¹⁰ m
- α_SI = 0.1 × R_S² ≈ 1.61×10¹⁹ m²
- α (time²) = α_SI / c² ≈ 180 s²

## 3. Dimensionless ratio — dziś i przez epoki

`8πG α ρ_m = 3 α H₀² Ω_m × (1+z)³`

| Epoch | (1+z) | ratio | safety margin |
|-------|-------|-------|---------------|
| Today | 1.00 | 8.08e−34 | 1.24e+33× |
| DE-matter equality (z=0.3) | 1.30 | 1.78e−33 | 5.63e+32× |
| Recombination (z=1090) | 1.09e+3 | 1.05e−24 | 9.53e+23× |
| Matter-radiation equality (z=3400) | 3.40e+3 | 3.18e−23 | 3.15e+22× |
| BBN (z=4×10⁸) | 4.00e+8 | 5.17e−08 | 1.93e+07× |
| z=10¹⁵ | 1.00e+15 | 8.08e+11 | PHANTOM |
| z=10²⁵ (pre-Planck) | 1.00e+25 | 8.08e+41 | PHANTOM |

**Phantom transition:** `z_phantom = (1/ratio_today)^(1/3) − 1 ≈ 1.07×10¹¹`

## 4. Critical findings

### Finding 1: Observable cosmology fully protected

Najbardziej restrictive observable epoch to BBN (z ≈ 4×10⁸): ratio
5.17e−08 → margines 1.93×10⁷ poniżej phantom threshold. Dla wszystkich
obserwowalnych skal kosmologicznych (DE era, recombination, BBN) sign
`K_eff` pozostaje NEGATIVE (canonical kinetic structure niezmieniona).

⇒ de2 theorem `w_psi ≥ −1` valid w całej observable history.

### Finding 2: Phantom transition pre-BBN, ale post-Planck

Formalna phantom transition zachodzi przy z ≈ 10¹¹, co jest:
- **Po** BBN (z ≈ 4×10⁸) — beyond observable nucleosynthesis era
- **Przed** Planck epoch (z ≈ 10³²) — w régime gdzie classical FRW i tak
  się załamuje

W régime z > 10¹¹ Candidate D analiza klasyczna nie ma zastosowania
(quantum gravity, pre-equilibrium substrate dynamics). To jest rather
**limitation of classical cosmology**, nie prediction Candidate D.

### Finding 3: OP-DESI niezależność preserved

Pod Candidate D pivot:
- Φ field ewolucja wciąż canonical (w obserwowalnym régime)
- de2 theorem pozostaje active prediction
- DESI DR3 phantom crossing (jeśli zaobserwowane) **falsyfikuje**
  both M9.1'' AND M9.2-D simultaneously

⇒ OP-DESI nie staje się redundant pod M9.2 pivot. Phantom crossing
observation byłby ortogonalny falsification path do ngEHT shadow test.

### Finding 4: Calibration robustness

Result ratio_today ~ 10⁻³⁴ jest insensitive do calibration choice α:
- Even α ~ 10⁹ × current value (np. M_BH solar mass scale) keeps
  observable cosmology safe
- Phantom transition z scales as α^(1/3) — would need α ~ 10²⁰× larger
  to push z_phantom below recombination

⇒ Result is robust nawet jeśli future Phase 1 analysis revises α.

## 5. Phase 1 deferred items (cosmology-specific)

Phase 0+ cosmology cross-check zamyka **partial** Phase 1 item (d):

✅ **Closed (Phase 0+):**
- (d.1) Czy back-reaction wprowadza phantom crossing? **NIE** (margines >10⁷×)
- (d.2) OP-DESI niezależność: **PRESERVED**

⏳ **Still deferred (full Phase 1):**
- (d.3) DESI w(z) prediction shifts: czy α H² coupling wprowadza scale
  dependence w(z)? Wymaga full perturbation theory analysis.
- (d.4) Cross-impact OP-Hubble: czy α-coupling wpływa na H(z) kinematics?
- (d.5) Φ kinetic correction beyond-tree: 1-loop counterterms w FRW

## 6. Implications strategiczne

### Phase 0+ closure status

- ✅ **Phase 0+ kickoff** (rano 2026-04-25): structural sketch 5/5 POSITIVE
- ✅ **Phase 0+ cosmology** (popołudnie 2026-04-25): w(z) preservation POSITIVE
- ⏳ **Phase 1** (deferred 2026 Q3-Q4 lub post-ngEHT 2030+): full covariant
  derivation + photon ring + perturbations + 5th force

### F4 falsifiability hardening (post Phase 0+ cosmology)

> **TGP M9.1'' standalone falsifiable (ngEHT 2030+ + DESI DR3). M9.2-D
> pivot path:**
> - Strong-field: Candidate D structural sketch POSITIVE (α ≈ 0.1 tunes
>   scenario (e) target +1.46% deviation)
> - Weak-field: U⁴ auto-suppressed (margines 9e+12× nad Cassini)
> - Cosmology: w(z) ≥ −1 preserved (margines 1.24e+33× w obserwowalnym régime)
> - GW propagation: c_GW = c_0 vacuum exact (OP-7 unchanged)
> - Stability: no-ghost α > 0, Ostrogradsky-free at tree level
>
> **Cross-program independence preserved:** OP-DESI pozostaje niezależnym
> testem TGP (M9.1'' OR M9.2-D), bez redundancji.

## 7. Files (Phase 0+ cosmology)

- `op_m92_P0plus_candD_cosmology.py` — analytical cross-check script
- `op_m92_P0plus_candD_cosmology.txt` — raw output
- `OP_M92_P0plus_candD_cosmology_results.md` — this synthesis

## 8. Cross-references

- [[research/op-m92/OP_M92_P0plus_candD_results.md]] — Phase 0+ kickoff (parent)
- [[research/op-m92/OP_M92_readiness_summary.md]] — Phase 0 closure
- [[research/op-m92/OP_M92_setup.md]] — formal scope
- [[research/op-desi/OP_DESI_de2_theorem.md]] — w(z) ≥ −1 theorem (target)
- [[research/op-hubble/]] — H_0 cross-impact (Phase 1 item d.4)

---

## Bottom line

OP-M92 Phase 0+ cosmology cross-check **ZAMKNIĘTY POSITIVE**. Candidate D
back-reaction term `α T·J·J` w obserwowalnej kosmologii daje correction
~10⁻³⁴ (dziś) do ~10⁻⁸ (BBN), z phantom transition odepchniętą do z ≈ 10¹¹
(pre-BBN, unobservable régime).

**OP-DESI prediction `w(z) ≥ −1` STRUKTURALNIE PRESERVED** pod M9.2-D pivot.
DESI DR3 phantom crossing observation falsifies both M9.1'' i M9.2-D —
niezależny test, brak redundancji.

**M9.2 pivot path post-ngEHT 2030+ verdict:** Candidate D leads with
**both** strong-field tunability AND cosmological compatibility
strukturalnie pre-derived at sketch level. Response time 2-4 weeks
zamiast 2-3 lat.
