---
title: "M10.4 — CMB safety REBUILD (canonical scalar Φ, replaces gs41 f(R))"
date: 2026-04-26
cycle: M10.4
status: SETUP
predecessor: "[[M10_3_results.md]] (6/6 PASS, gs66 YELLOW → GREEN)"
audit_target: "[[../galaxy_scaling/gs41_cmb_compatibility.py]] (RED — uses f(R), to supersede)"
related:
  - "[[M10_program.md]]"
  - "[[M10_0_drift_audit.md]]"
  - "[[M10_3_results.md]] (FRW propagator + Yukawa screening)"
  - "[[../op-newton-momentum/M9_3_setup.md]] (M9.3.1 m_s² = β)"
  - "[[../closure_2026-04-26/Lambda_from_Phi0/results.md]] (T-Λ: V_eq=β/12, Φ_eq=H_0)"
tags:
  - TGP
  - M10
  - CMB
  - cosmology
  - rebuild
  - audit-cycle
---

# M10.4 — CMB safety REBUILD (canonical scalar Φ)

> **REBUILD, nie audit.** [[../galaxy_scaling/gs41_cmb_compatibility.py]] używa **f(R)** gravity — strukturalnie sprzeczne z TGP single-Φ axiom (sek08a). M10.4 zastępuje gs41 **canonical scalar Φ** formulation z fundamentalnymi zasadami M9.3.1 + T-Λ + M10.3.

---

## 1. Dlaczego REBUILD

Z [[M10_0_drift_audit.md]] gs41 status: **RED**.

| Aspect | gs41 (f(R)) | sek08a / canonical TGP |
|--------|-------------|------------------------|
| Framework | `f(R) = R + R₀^γ R^(1-γ) exp(-(R/R₀)^α)` | Single-Φ scalar field |
| Field | Spacetime curvature R modified | Scalar Φ on hyperbolic metric |
| Mass scale | `m_s²(R) = (1+f_R)/(3f_RR) - R/3` | `m_s² = β` (constant, from M9.3.1) |
| Screening | R-curvature exponential `exp(-(R/R₀)^α)` | Yukawa screening + Hubble friction |
| Vacuum | Implicit via f(R) profile | `V(1) = β/12` (T-Λ residual) |

**Diagnoza:** gs41 not derivable from sek08a action; uses different axioms. **NIE** poprawiamy gs41 — **odbudowujemy** w canonical scalar Φ.

---

## 2. Strukturalny punkt wyjścia (canonical TGP)

### 2.1 Field equations (M9.3.1 + M10.3 confirmed)

Linearized δΦ EOM around vacuum `Φ = Φ_0` (M10.3.1 sympy proof):
```
δΦ̈ + 3H·δΦ̇ − (1/a²)∇²δΦ + β·δΦ = source
```

Fourier mode `δΦ_k(t)`:
```
δΦ̈_k + 3H·δΦ̇_k + (k²/a² + β)·δΦ_k = source_k
```

### 2.2 Scale of β from T-Λ closure

Z T-Λ closure ([[../closure_2026-04-26/Lambda_from_Phi0/results.md]]):
- `V(Φ_0) = V_eq = β/12` (in normalized units)
- `ρ_vac,TGP = M_Pl²·H_0²/12` (matches Ω_Λ ≈ 0.685)
- ⇒ `β·M_Pl⁴/12 ~ M_Pl²·H_0²/12` (with V dimensionless multiplied by Φ_0⁴ ~ M_Pl⁴)
- ⇒ **`β ~ H_0²`** (in mass-squared dimension; `m_s ~ H_0`)

**Compton wavelength:** `λ_Compton ~ 1/√β ~ 1/H_0 ~ L_H` (Hubble radius today).

### 2.3 Reżim Hubble vs Yukawa

Two regimes on cosmological linear scale `k`:

| Regime | Condition | δΦ behavior |
|--------|-----------|-------------|
| **Hubble-frozen** | `k²/a² + β ≪ H²` | δΦ frozen by friction (slow-roll); follows ΛCDM background |
| **Sub-horizon Yukawa** | `k²/a² ≫ β, H²` | δΦ propagates as massive Yukawa; sub-horizon "fifth force" suppressed by `1/(k²/a² + β)` |

**Klucz CMB safety:** at recombination `H_rec ≈ H_0·√(Ω_m·z_rec³) ≈ H_0·10⁵`, we have `β ~ H_0² ≪ H_rec²` ⇒ **wszystkie super-horizon δΦ modes są Hubble-frozen** ⇒ no fifth-force on CMB scales.

---

## 3. CMB safety mechanisms (canonical TGP)

| Effect | gs41 (f(R)) | Canonical TGP |
|--------|-------------|---------------|
| **BBN** (T~MeV, z~10⁹) | R/R₀ huge → exp(−x^α) kills it | H_BBN ≫ √β, Φ frozen at Φ_0 (Hubble friction) |
| **Recombination** (z~1100) | R/R₀ ~ 10¹⁰ → exp suppression | H_rec ~ 10⁵ H_0 ≫ √β, δΦ frozen |
| **ISW** (z~0-2) | f_R correction sub-dominant | δΦ on Hubble scales O(δΦ/Φ_0) — testable! |
| **σ_8 / growth** | f_RR enhances scalaron force | Yukawa-screened sub-horizon, no enhancement |
| **CMB lensing** | Slip parameter η = (Φ−Ψ)/Φ ≈ 0 | Same: scalar Φ doesn't introduce slip at linear order |

---

## 4. Plan sub-testów

### M10.4.1 — Symbolic m_s² derivation (sympy)

Cel: pokazać że m_s² = β jest constant (NIE R-curvature dependent jak w f(R)) i scale ~ H_0² z T-Λ.

PASS: m_s² = β symbolicznie; numeryczne β ≈ H_0² consistent z T-Λ Ω_Λ=0.685.

### M10.4.2 — Background evolution: Φ at vacuum (numerical FRW)

Cel: pokazać że Φ pozostaje przy Φ_0 podczas cosmic history (radiation/matter eras).

Metoda:
- Solve Klein-Gordon w FRW: `Φ̈ + 3H·Φ̇ + V'(Φ)/K(Φ) = 0`
- IC: `Φ(z=10⁹) = Φ_0(1+ε)`, `Φ̇=0` (frozen at BBN)
- Track `δΦ/Φ_0` through z=10⁹ → 10³ → 0
- Background `H(z)` z standard ΛCDM (Ω_m=0.315, Ω_r=9.1e-5, Ω_Λ=0.6847)

PASS: `|δΦ(z)/Φ_0| < 10⁻⁵` at all z (Hubble-frozen).

### M10.4.3 — Linear perturbations at recombination

Cel: pokazać że super-horizon δΦ_k modes są Hubble-frozen przy z=1100.

Metoda: ODE solve dla δΦ_k(η) z `δΦ̈_k + 3H·δΦ̇_k + (k²/a² + β)δΦ_k = 0` dla scali k = {0.001, 0.01, 0.1, 1} h/Mpc.

PASS: δΦ_k frozen for super-horizon `k < a·H_rec`; subhorizon decay z Yukawa rate `√β + k/a`.

### M10.4.4 — ISW estimate at z<2

Cel: oszacować ΔT_ISW/T z δΦ-induced Φ_grav modification along photon path.

Formula (linearized): `ΔT_ISW/T = 2 ∫(Ψ̇ + Φ̇) dη` where `Φ` here = grav. potential, not TGP scalar.

W TGP: scalar δΦ_TGP couples to matter via `q·δρ` (sek08a coupling), gives modyfikację gravitational potential `δΨ ~ q·Φ_0·δΦ_TGP/M_Pl²`.

PASS: |ΔT_ISW^TGP / ΔT_ISW^ΛCDM − 1| < 0.05 (5% modification limit, currently within Planck error bars).

### M10.4.5 — Growth rate σ_8 and S_8

Cel: pokazać że TGP nie podnosi σ_8 vs ΛCDM (S_8 tension neutral, NIE worsened).

Metoda: linear growth `D(a)` solving `D̈ + 2H·Ḋ − (3/2)Ω_m H²·D = δF_TGP`; show `δF_TGP` Yukawa-suppressed sub-horizon ⇒ growth ≈ ΛCDM.

PASS: `|σ_8^TGP/σ_8^ΛCDM − 1| < 10⁻³` at z=0 for sub-horizon modes.

### M10.4.6 — Honest synthesis verdict

Cel: dokumentacja CMB-safety mechanism canonical TGP, falsifiable predictions.

PASS: 5/5 sub-tests OK + dokumentacja + statement że gs41 (f(R)) jest **superseded** przez M10.4 (canonical Φ).

---

## 5. Foundational consistency checks

| Constraint | Test sub | Expected |
|------------|----------|:--------:|
| Single-Φ axiom | M10.4.1-5 use scalar Φ only (no f(R)) | ✅ |
| β=γ vacuum cond. | V'(Φ_0) = 0 used in M10.4.2 background | ✅ |
| T-Λ residual `V_eq = β/12` | M10.4.1 derives `β ~ H_0²` from this | ✅ |
| K(φ) = K_geo·φ⁴ | M10.4.2 sub-leading near vacuum (z M10.1) | ✅ |
| M9.3.1 stable Yukawa m_s² = +β | M10.4.1 derives this; M10.4.3 uses for δΦ_k | ✅ |
| Hyperbolic metric M9.1'' | δΨ from δΦ uses linearization (M9.3 result) | ⏳ |

---

## 6. Hipoteza pre-test

**Oczekiwany verdict M10.4: 5/5 PASS** (lub 6/6 z synthesis):

- TGP CMB-safe **strukturalnie** przez:
  - Hubble friction freezing super-horizon Φ-modes (m_s ~ H_0 ≪ H przed today)
  - Yukawa screening sub-horizon δΦ-modes (Compton scale ~ L_H)
- ISW: O(δΦ/Φ_0) modification on Hubble scales — **falsifiable** by future surveys (LiteBIRD, CMB-S4)
- Growth rate: ΛCDM-like (no enhancement) — `S_8 tension nie pomaga` (chameleon-like screening)
- BBN: Φ frozen at Φ_0 — no modification to expansion rate

**Verdict gs41:** RED → **SUPERSEDED** (NIE upgraded; framework f(R) jest different theory). M10.4 staje się canonical reference for TGP CMB safety.

---

## 7. Files

| Plik | Status | Cel |
|------|--------|-----|
| [[M10_4_setup.md]] | NEW (this) | Plan + math foundation |
| [[m10_4_cmb.py]] | TO CREATE | REBUILD script (5+1 tests) |
| [[m10_4_cmb.txt]] | TO CREATE | Run log |
| [[M10_4_results.md]] | TO CREATE | Closure-grade synthesis |
| [[../galaxy_scaling/gs41_cmb_compatibility.py]] | UNCHANGED | Marked SUPERSEDED w results |

---

## 8. Następne (post-M10.4)

- **M10.5:** ct3+ct7 H₀/S₈ tensions audit (verify K=φ⁴ correction)
- **M10.R:** synthesis cyklu M10 + cross-check vs T-Λ + update PLAN_ROZWOJU_v4

---

*M10.4 setup ready 2026-04-26. Proceed with REBUILD script.*
