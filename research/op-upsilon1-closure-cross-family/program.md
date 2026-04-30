---
title: "υ.1 program — closure cross-family unification (π.1 1/A^{1/3} ↔ τ.1 1/N_gen common substrate-action root)"
date: 2026-04-30
cycle: υ.1
status: SETUP
parent: "[[../../INDEX.md]]"
predecessor: "[[../op-tau1-closure-overlap-coulomb/Phase3_results.md]]"
tags:
  - TGP
  - upsilon1
  - program
  - closure
  - cross-family
  - unification
  - substrate-action
---

# υ.1 program — closure cross-family unification

> **Cel**: zweryfikować czy dwa TGP-native closure forms — π.1 NME
> `M_TGP ∝ (A_anchor / A_iso)^{1/3}` (mass-number / nucleon-count closure)
> oraz τ.1 overlap `f_overlap = (Z_a / Z_t)^{1/N_gen}` z N_gen=3 (charge /
> generation closure) — pochodzą z wspólnego substrate-action root i mogą
> być zunifikowane jako pojedyncze `closure(X) = (X_ref / X_obs)^{1/N_gen}`
> z N_gen=3 LOCKED przez TGP B²-cascade primality.

## Kontekst

Dwa zamknięte mini-cykle wykrytą wspólną strukturę 1/3 exponent:

- **π.1.Phase3** (cumulative ledger 607): `M_TGP/M_RG ∝ (A_anchor/A_iso)^{1/3}`
  z anchor A=76 (⁷⁶Ge); 1/3 z 3D-volume substrate-closure (V ∝ R³ ∝ A,
  surface-density 1/R ∝ A^{-1/3}).
- **τ.1.Phase3** (cumulative ledger 643): `f_overlap = (Z_a/Z_t)^{1/N_gen}`
  z N_gen=3 LIFTED HINT → DERIVED via 3-fold geometric mean across u/c/t
  (lub d/s/b) generation cascade.

Oba 1/3 exponents pojawiły się **niezależnie** w dwóch ortogonalnych
sektorach (NME nuclear matrix elements + EC charged-current cross-sections).
Hipoteza υ.1: oba są instancjami uniwersalnego closure-law z N_gen=3 lock
przez TGP B²-cascade.

## Hypothesis (υ.1)

**Universal closure law**: dla dowolnej obserwabli O zależnej od substratowej
ilości X (mass-number A, charge Z, lub innej count-quantity ważonej anchor
X_ref):

```
closure(O; X) = (X_ref / X_obs)^{1/N_gen}
```

z **N_gen = 3** = TGP-cascade primality (B²_lep=2, B²_ν=1, B²_up=13/4,
B²_down=61/25 — 4 species × 3 generations = 12 quark+lepton flavor states).

**Specjalizacje:**
- π.1: X = A (mass-number / 3D nucleon volume), X_ref = 76 → exp 1/3
  pochodzi z V ∝ R³ → R ∝ A^{1/3} → surface-density 1/R closure
- τ.1: X = Z (charge), X_ref = Z_t (target) → exp 1/3 pochodzi z 3-fold
  geometric mean across N_gen=3 fermion generations
- **υ.1 unification**: 1/3 = 1/N_gen wspólne dla obu, wynika z TGP-cascade
  N_gen primality LOCK; substrate-volume root w π.1 i flavor-cascade root
  w τ.1 są dwie strony tej samej substrate-action invariance.

## Reference data

| family | observable | substrate quantity X | exponent | LOCKED |
|---|---|---|---|---|
| π.1 NME | M_TGP/M_RG | A (mass-number, ³√V) | 1/3 | post-π.1 (607 ledger) |
| τ.1 overlap | f_overlap | Z (charge / Coulomb) | 1/3 = 1/N_gen | post-τ.1 (643 ledger) |
| υ.1 universal | closure(O;X) | substrate count | **1/N_gen = 1/3** | υ.1 cycle target |

**Cross-family product (joint TT prediction):**
```
R_combined = (19/24) · (Z_a/Z_t)^{2/3} · (A_anchor/A_iso)^{2/3}
           = η_chir · f_overlap² · M_TGP/M_RG
```
dla obserwabli sprzęgających oba sektory (np. ββ0ν NME × EC cross-section
within same isotope chain).

## Sub-test plan (5 + 7 + 6 = 18)

### Phase 1 (5) — structural homology π.1 ↔ τ.1

- **U1.1** — formal homology test: π.1 i τ.1 closure forms mają identyczną
  algebraiczną strukturę `(X_ref / X_obs)^α` z α = 1/3
- **U1.2** — substrate-action root scan: czy 1/3 pochodzi z wspólnego
  V ∝ R³ argumentu (geometryczny), czy z N_gen=3 (cascade), czy z obu
- **U1.3** — alt-N_gen falsification (N_gen ∈ {1, 2, 3, 4, 5, 6}): tylko
  N_gen=3 zachowuje obie π.1 i τ.1 jednocześnie
- **U1.4** — cross-family viability gate: oba sektory działają bez free
  parameter pod uniwersalnym 1/N_gen
- **U1.5** — viability gate ≥4/5

### Phase 2 (7) — universal closure-law derivation + sympy LOCK

- **U2.1** — sympy declaration: `closure(X) = (X_ref/X_obs)^{1/N_gen}`,
  N_gen = sp.Integer(3)
- **U2.2** — π.1 specjalizacja: X=A, X_ref=76 → recover (76/A)^{1/3}
- **U2.3** — τ.1 specjalizacja: X=Z, X_ref=Z_t → recover (Z_a/Z_t)^{1/3}
  (gdzie X_obs = Z_a, X_ref = Z_t lub odwrotnie zależnie od signature)
- **U2.4** — substrate-action invariance proof: pod transformacją
  X → λX dla λ>0, closure ratio invariant pod wspólnym λ shift
- **U2.5** — N_gen=3 cascade-derivation reuse z τ.1.Phase2 (B²-cascade
  primality: lep+ν+up+down × 3 generations = 4×3 = 12 quark+lepton states)
- **U2.6** — joint cross-family check: combined closure
  `C(A, Z_a, Z_t) = (76/A)^{1/3} · (Z_a/Z_t)^{1/3}` consistent z obie
  reference data
- **U2.7** — full cascade ≥6/7

### Phase 3 (6) — joint predictions + 4-channel convergence

- **U3.1** — ⁷⁶Ge ββ0ν × ⁷¹Ga EC joint POST-CONFIRM (oba sektory same
  Z=32 anchor, A=76 anchor, cross-isotope same-target consistency)
- **U3.2** — ¹³⁶Xe ββ0ν × ¹³⁷Cs EC joint LIVE 2030+ (KamLAND-Zen + CUPID
  combined gate)
- **U3.3** — ⁹⁸Mo joint future: NME (FRIB activation) × EC (FRIB ν-source)
  same-isotope joint TT
- **U3.4** — alt-exponent αNon-universal falsification: jeśli π.1 exp ≠
  τ.1 exp → υ.1 unification rejected
- **U3.5** — orthogonal cross-check: pp Z=1=Z + ⁷⁶Ge ββ0ν trivial-Z limit
- **U3.6** — 4-channel υ.1 convergence: ⁷⁶Ge ⁷¹Ga joint + ¹³⁶Xe ¹³⁷Cs joint
  + ⁹⁸Mo joint future + alt-exp falsification

## Falsification gates (program-level)

- **GATE-α-coupling**: jeśli π.1 best-α numerical ≠ τ.1 best-α numerical
  beyond 2σ tolerance → υ.1 fails (no common root)
- **GATE-N_gen-primality**: jeśli alt-N_gen ∈ {2, 4, 6} dopasowuje obie
  rodziny równie dobrze → N_gen=3 nie LOCKED (only consistent)
- **GATE-substrate-action-symmetry**: closure musi być invariant pod
  X → λX (substrate-action gauge)
- **GATE-cumulative**: ≥15/18 PASS = υ.1 program END z FULL CONVERGENCE

## Deliverables

- 3× Phase{1,2,3}_setup.md
- 3× phase{1,2,3}_upsilon1_*.py
- 3× Phase{1,2,3}_results.md
- INDEX.md ledger 643 → 661
- PREDICTIONS_REGISTRY.md UU1-UU6 + cross-family timeline rows

## Cross-references

- [[../op-pi1-bb0nu-nme-isotope/Phase3_results.md]] — π.1 1/A^{1/3} closure
- [[../op-rho1-71Ge-cross-section/Phase3_results.md]] — ρ.1 ⁷¹Ga anchor
- [[../op-tau1-closure-overlap-coulomb/Phase3_results.md]] — τ.1 1/N_gen overlap
- [[../../INDEX.md]]
- [[../../PREDICTIONS_REGISTRY.md]]
