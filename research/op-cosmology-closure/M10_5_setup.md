---
title: "M10.5 — H_0/S_8 tensions audit (ct3 + ct7 closure)"
date: 2026-04-26
cycle: M10.5
status: SETUP
predecessor: "[[M10_4_results.md]] (6/6 PASS, gs41 RED → SUPERSEDED)"
audit_targets:
  - "[[../cosmo_tensions/ct3_dark_matter_backreaction.py]] (YELLOW — canonical K=1; misleading tachyonic interpretation)"
  - "[[../cosmo_tensions/ct7_soliton_cosmology.py]] (YELLOW — canonical K=1; HONEST verdict basis)"
related:
  - "[[M10_program.md]]"
  - "[[M10_0_drift_audit.md]]"
  - "[[M10_3_results.md]] (M_eff² = +β stable Yukawa)"
  - "[[M10_4_results.md]] (β ~ H_0² scale)"
  - "[[../op-newton-momentum/M9_3_setup.md]] (M9.3.1 stable linearization)"
tags:
  - TGP
  - M10
  - cosmology
  - tensions
  - audit-cycle
  - H0
  - S8
---

# M10.5 — H_0/S_8 tensions audit (ct3 + ct7)

> **Audyt łączony.** ct3 i ct7 oba mają drift YELLOW (canonical K=1 vs sek08a K=φ⁴). Najważniejsze: **ct3 jest misleadingly optimistic** ("amplification w right range dla H_0!"), podczas gdy **ct7 jest honestly negative** ("TGP scope = galaxy-scale, NIE cosmology"). M10.5 ustala że ct7 verdict jest **strukturalnie poprawny** po uwzględnieniu M9.3.1 stable linearization.

---

## 1. Stan wyjściowy: ct3 vs ct7

### 1.1 ct3 verdict (optimistic, oparty na błędnej interpretacji)

ct3 zakłada **tachyonic amplification + friction balance**:
- Linearizacja: `V(ψ) = -γ/2 (ψ-1)²` ⇒ "tachyonic mass" `m²_tach = -V''(1) = γ`
- Amplification: `exp(m_tach/H_0) ~ exp(3.5) ~ 33` per Hubble time
- Saturation z `3ψ̇²/ψ`: `δψ_sat ~ 1/3`, `B_ψ/H_0² ~ 0.03-0.3`
- **ct3 conclusion:** "B_ψ in right range for H_0 tension (0.17)!"

**Problem:** "Tachyonic mass" interpretation jest **misleading** — V''(1)=-γ jest **slow-roll maximum** w cosmological time (z K=K_geo·φ⁴ + hyperbolic g_eff), **NIE** spatial tachyonic. Spatial linearization foundations Φ-EOM (M10.3.1):
```
M_eff² = +β > 0  (stable Yukawa, M9.3.1)
```
Spatial fluctuations są **screened**, nie amplified.

### 1.2 ct7 verdict (honest, structural)

ct7 ostatecznie odrzuca tachyonic amplification (po ct5/ct6):
- `B_ψ/H_0² ~ 10⁻⁹` (8-9 orders below needed 0.17)
- Scale mismatch: tachyonic only at super-Hubble (5600 Mpc)
- Soliton dilution: `d/λ_C ~ 10¹⁸` — too dilute for collective effects
- RG running: `Δγ/γ ~ 0.5%` (vs needed 8% for H_0 shift)
- **Phantom crossing:** `B_ψ ≥ 0` ⇒ `w_eff > -1` always; TGP **strukturalnie nie może** produce DESI w₀ = -0.45

**ct7 conclusion:** "TGP scope = galaxy dynamics (y < 1), NOT cosmology (y >> 1). Honest position: TGP does not explain H_0/S_8/DESI tensions."

### 1.3 Reconciliation

| Aspekt | ct3 (optimistic) | ct7 (honest) | M10.5 verdict |
|--------|------------------|--------------|---------------|
| Linearization | tachyonic V''=-γ | tachyonic V''=-γ | **stable M_eff²=+β** (M9.3.1) |
| B_ψ/H_0² | ~0.03-0.3 (matches!) | ~10⁻⁹ (8-9 orders gap) | **~10⁻⁹** (ct7 correct) |
| Mechanism | tachyonic amplification | structural impossibility | **Yukawa screening** structural |
| Scope | partial cosmology success | galaxy-scale only | **galaxy-scale** (ct7 reaffirmed) |

---

## 2. Math foundation (canonical TGP)

### 2.1 Sek08a kinetic K=φ⁴ (verify drift correction)

Pełen kinetic term sek08a:
```
S_kin = ∫ d⁴x √(-g_eff) · (1/2) K(φ) g_eff^μν ∂_μφ ∂_νφ,    K(φ) = K_geo·φ⁴
```

Near vacuum φ = 1 + δ:
```
K(1+δ) = K_geo·(1 + 4δ + 6δ² + 4δ³ + δ⁴)
       = K_geo·(1 + 4δ + O(δ²))
```

**Linear order:** K renormalizes by const factor 1 + 4δ near vacuum; this enters EOM jako multiplikatywny faktor przed `□φ` term.

**Backreaction implication:** B_ψ pochodzi z `<K(φ)·(∂φ)²>` average. Z K = K_geo·(1 + 4δ + ...):
```
<K·(∂φ)²> = K_geo·<(∂δ)²> + 4·K_geo·<δ·(∂δ)²> + O(δ²·(∂δ)²)
```

Linear `<δ·(∂δ)²>` term jest cubic w fluctuations ⇒ vanishes for Gaussian random δ. Pierwsza nontrivial K=φ⁴ correction is `<δ²·(∂δ)²>` → O(σ⁴) = `O((δ_grav)⁴)` ~ 10⁻²⁰. **Sub-leading**, nie zmienia ct7 verdict.

### 2.2 Spatial M_eff² = +β (M9.3.1, M10.3 confirmed)

Foundations Φ-EOM linearization (M10.3.1):
```
∇²δ + (β/Φ_0)·(2 - 3) δ + O(δ²) = -q·ρ
∇²δ - β·δ = source  ⇒  M_eff² = +β  (stable Yukawa)
```

**Implications dla backreaction:**
- Spatial Yukawa: `δψ ~ exp(-r√β)/r` decays exponentially
- NO spatial fluctuation growth (NIE tachyonic instability)
- ct3 "amplification" mechanism falsified

### 2.3 Cosmological slow-roll (z K=φ⁴ + hyperbolic g_eff)

Cosmological perturbation equation (FRW background):
```
δ̈ + 3H·δ̇ + β·δ = source   (z M9.3.1 slow-roll max V''(1)=-γ but stable in M10.3)
```

To **cosmological slow-roll**, NIE spatial tachyonic. Hubble friction `3H·δ̇` damps fluctuations on `t ~ 1/(3H)`. Stable at vacuum.

### 2.4 B_ψ formula (canonical, post-M9.3.1)

Substrate backreaction (Buchert-style):
```
B_ψ = 3·<(δψ̇)²>/<ψ> - 3·<δψ̇>²/<ψ>
    ≈ 3·<(δψ̇)²>     (for <δψ̇>=0 at leading order)
```

Z linearized δ̈ + 3Hδ̇ + βδ = src, w stationary regime:
```
δψ̇_typical ~ H·δψ_typical    (from Hubble friction balance)
<(δψ̇)²> ~ H²·<(δψ)²> = H²·σ²_δ
```

Z `δψ ~ 2Φ_grav/c²`, σ_δ ~ 10⁻⁵ (from CMB potential):
```
B_ψ/H_0² ~ 3·σ²_δ ~ 3×10⁻¹⁰ ~ 10⁻⁹  (matches ct7 estimate!)
```

**Required dla H_0 = 73 (Δ ≈ 8.3%):** `B_ψ/H_0² ≈ 0.17`.
**Gap: 8 orders of magnitude.** Confirmed ct7 verdict.

### 2.5 w(z) phantom crossing

DESI DR1: `w_0 = -0.45, w_a = -1.79` ⇒ `w(z=0.5) ~ -1.5` (phantom).

TGP w_eff:
```
w_eff = -1 + ρ_kinetic/ρ_total
ρ_kinetic = (1/2)·K(φ)·φ̇² ≥ 0   (positive definite)
⇒ w_eff ≥ -1   strukturalnie
```

**Strukturalna impossibility:** TGP nie może produce w_eff < -1 dopóki `K(φ) ≥ 0`. Z K=K_geo·φ⁴ near φ=1: K = K_geo·(1+δ)⁴ > 0 wszystkie δ > -1.

⇒ **DESI phantom crossing strukturalnie wyklucza TGP DE** (jeśli sygnał potwierdzony przez DR2/DR3).

---

## 3. Plan sub-testów

### M10.5.1 — K=φ⁴ correction symbolic check (sympy)

Cel: pokazać że K=φ⁴ correction near vacuum jest sub-leading vs canonical K=1.

PASS: linear-in-δ correction = 4·K_geo·δ; `<δ·(∂δ)²> = 0` for symmetric Gaussian; pierwsze nontrivial K-correction to O(σ⁴).

### M10.5.2 — Spatial linearization correctness (sympy + numerical)

Cel: confirm M_eff² = +β > 0 stable (M9.3.1) supersedes ct3 "tachyonic amplification" interpretation.

PASS: foundations Φ-EOM linearization gives M_eff² = +β; Yukawa screening exp(-r√β); ZERO amplification z spatial fluctuations.

### M10.5.3 — B_ψ/H_0² numerical verification (corrected interpretation)

Cel: z poprawną M_eff²=+β linearization, oszacować B_ψ z grav potential variance.

Metoda:
- σ_δ = 2σ_Φ/c² ~ 6×10⁻⁵ (rms grav potential variance from CMB+LSS)
- δψ̇ ~ H·δψ (Hubble friction balance)
- B_ψ/H_0² ≈ 3·σ²_δ

PASS: B_ψ/H_0² ~ 10⁻⁹ (matches ct7 prediction; 8-9 orders below 0.17 needed dla H_0 tension).

### M10.5.4 — RG running γ(k) verification

Cel: verify ct7 RG calculation: `Δγ/γ ~ 0.5%` between CMB and local scales.

Metoda: γ(k) = γ_UV·(k/k_UV)^η, η = 0.044 (LPA' anomalous dim z CG-2). Range k_CMB → k_local: ratio ~140, ln(140)·η = 0.22, exp(0.22) = 1.25 → 25% wait, let me check again.

Actually η·ln(140) = 0.044·4.94 ≈ 0.22, exp(0.22)-1 ≈ 0.24 — 24% change. Hmm, ale ct7 daje 0.5%? Let me check ct7 numbers in audit script.

### M10.5.5 — w_eff ≥ -1 strukturalna ograniczenie (sympy)

Cel: pokazać że w_eff ≥ -1 wynika strukturalnie z K(φ) ≥ 0 i nie da się obejść w canonical TGP.

PASS: w_eff = -1 + ρ_kin/ρ_tot, ρ_kin ≥ 0 ⇒ w_eff ≥ -1. DESI phantom (w₀=-0.45, w_a=-1.79 ⇒ w(0.5)≈-1.5) jeśli się utrzyma w DR2/DR3, **strukturalnie falsyfikuje** TGP DE.

### M10.5.6 — Honest synthesis verdict

Cel: dokumentacja:
- ct3 optimistic estimate SUPERSEDED (oparte na misinterpretation "tachyonic")
- ct7 honest verdict CONFIRMED i WZMOCNIONE przez M9.3.1 stable linearization
- TGP scope = galaxy-scale, NIE cosmology tensions

PASS: 5/5 sub-tests OK + dokumentacja + final verdict YELLOW → GREEN dla ct7 (correct conclusion preserved); ct3 status YELLOW → GREEN dla honest reframing (B_ψ ~ 10⁻⁹, NIE ~0.1).

---

## 4. Foundational consistency checks

| Constraint | Test sub | Expected |
|------------|----------|:--------:|
| Single-Φ axiom | M10.5.1-5 use scalar Φ only | ✅ |
| β=γ vacuum cond. | V'(Φ_0) = 0 used in M10.5.2 lineariz | ✅ |
| K(φ) = K_geo·φ⁴ | M10.5.1 verifies sub-leading correction | ✅ |
| M9.3.1 stable Yukawa M_eff²=+β | M10.5.2-3 used (NIE ct3 tachyonic) | ✅ |
| T-Λ closure β ~ H_0² | M10.5.3 numerical scale | ✅ |
| Hyperbolic g_eff M9.1'' | implicit w cosmological slow-roll | ⏳ |

---

## 5. Hipoteza pre-test

**Oczekiwany verdict M10.5:** 5/5 PASS lub 6/6 z synthesis.

**Strukturalny finding:**
- ct3 "amplification mechanism" SUPERSEDED przez M9.3.1 (M_eff²=+β stable)
- ct7 honest verdict CONFIRMED: TGP scope = galaxy-scale
- B_ψ/H_0² ~ 10⁻⁹ << 0.17 (8-9 order gap, structural)
- w_eff ≥ -1 strukturalnie ⇒ DESI phantom crossing (jeśli potwierdzone) FALSYFIKUJE TGP DE

**Verdict ct3:** YELLOW → **GREEN-honest** (audit ujawnia że "amplification" was misleading; correct interpretation gives ~10⁻⁹ matching ct7).
**Verdict ct7:** YELLOW → **GREEN** (honest verdict reaffirmed and structurally grounded).

---

## 6. Files

| Plik | Status | Cel |
|------|--------|-----|
| [[M10_5_setup.md]] | NEW (this) | Plan + math foundation |
| [[m10_5_tensions.py]] | TO CREATE | Audit script (5+1 tests) |
| [[m10_5_tensions.txt]] | TO CREATE | Run log |
| [[M10_5_results.md]] | TO CREATE | Closure-grade synthesis |
| [[../cosmo_tensions/ct3_dark_matter_backreaction.py]] | UNCHANGED | Marked: optimistic estimate based on misinterpretation; correct value B_ψ ~ 10⁻⁹ |
| [[../cosmo_tensions/ct7_soliton_cosmology.py]] | UNCHANGED | Marked: HONEST verdict CONFIRMED |

---

## 7. Następne (post-M10.5)

- **M10.R:** synthesis cyklu M10 (M10.0/1/2/3/4/5 → final results) + cross-check vs T-Λ + update PLAN_ROZWOJU_v4
- **M11:** kwantyzacja Φ + RG flow γ(k) (deferred from M10.5; LPA' analysis już mamy w ct7)
- **N-body TGP:** dedicated full-sek08a code (out of M10 scope)

---

*M10.5 setup ready 2026-04-26. Proceed with audit script.*
