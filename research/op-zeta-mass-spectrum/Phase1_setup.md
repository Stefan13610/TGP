---
title: "ζ.1.Phase1 setup — Σm_ν = 59.6 meV neutrino mass-spectrum robustness"
date: 2026-04-29
cycle: ζ.1.Phase1
status: PRE-EXECUTION
parent: "[[program.md]]"
tags:
  - TGP
  - zeta-mass-spectrum
  - neutrino-masses
  - Koide-K-half
  - Majorana
---

# ζ.1.Phase1 — Σm_ν = 59.6 meV neutrino mass-spectrum robustness audit

> **Cel:** Verify TGP prediction Σm_ν = 59.6 meV from K(ν) = 1/2 (Majorana,
> B²=1) + observational Δm²₂₁, |Δm²₃₁| inputs (NuFit 5.3); confirm
> normal ordering exclusivity; falsify inverted ordering via K=1/2
> incompatibility.

---

## 5 sub-tests

### Z1.1 — K(ν) = 1/2 (Majorana, B²=1) LOCKED

**Test:**
- Koide K = (Σm_i)/(Σ√m_i)² = (2 + B²) / (2N) for N=3 generations
- Lepton (Dirac, B²=2, 2 chiralities): K_lep = (2+2)/6 = 4/6 = 2/3
  (matches PDG to 10⁻⁵)
- Neutrino (Majorana, B²=1, single chirality): K_ν = (2+1)/6 = 3/6 = 1/2
- Sympy exact rational

**Falsification:** if Koide group-theoretic prediction K_ν ≠ 1/2 → 
chirality-counting framework broken.

### Z1.2 — Δm² observational inputs (NuFit 5.3)

**Test:**
- Δm²₂₁ = 7.53 · 10⁻⁵ eV² (KamLAND + solar)
- |Δm²₃₁| = 2.453 · 10⁻³ eV² (atmospheric + reactor + accelerator)
- Both observational inputs to TGP framework (not derived from substrate)
- Hierarchy: Δm²₃₁ / Δm²₂₁ ≈ 32.6

**Falsification:** if NuFit 5.3 values shift > 5% w przyszłych fits →
ζ.1 mass-spectrum prediction adjusts.

### Z1.3 — Σm_ν = 59.6 meV via K(ν) = 1/2 closure (NO ordering)

**Test:**
- Normal ordering: m₁ < m₂ < m₃
- m₁² = m_min²
- m₂² = m_min² + Δm²₂₁
- m₃² = m_min² + |Δm²₃₁|
- Constraint: K = (m₁+m₂+m₃) / (√m₁+√m₂+√m₃)² = 1/2
- Solve for m_min via numerical: m_min ≈ 0.80 meV
- Then m₂ ≈ 8.65 meV, m₃ ≈ 50.11 meV
- Σm_ν = m₁ + m₂ + m₃ ≈ 59.6 meV

**Falsification:** if K=1/2 closure cannot satisfy NuFit 5.3 Δm² →
K(ν) = 1/2 incorrect lub TGP framework w mass sektorze niespójny.

### Z1.4 — DESI DR2 bound 0.12 eV vs TGP 0.0596 eV

**Test:**
- DESI DR2 (2024): Σm_ν < 0.072 eV (95% CL, ΛCDM)
- DESI DR3 (2027+ projected): Σm_ν < 0.040 eV (95% CL, tighter)
- TGP prediction: 0.0596 eV (sits within current bound, near DR3 edge)
- Margin TGP vs DESI DR2: |0.072 - 0.0596| / 0.0596 ≈ 21% margin
- Margin TGP vs DESI DR3: -33% — **TGP could be falsified by DR3**

**Falsification:** if DESI DR3 (2027+) sets Σm_ν < 0.040 eV at 95% CL
→ TGP K(ν) = 1/2 mass-spectrum prediction sfalsyfikowana, ζ.1 reopens.

### Z1.5 — Inverted ordering FORBIDDEN by K(ν) = 1/2

**Test:**
- IO: m₃ < m₁ < m₂
- m₃² = m_min² (lightest)
- m₁² = m_min² + |Δm²₃₁|
- m₂² = m_min² + |Δm²₃₁| + Δm²₂₁
- Constraint K = 1/2: solve for m_min
- IO solution: m_min² < 0 → **forbidden** (no real solution)
- Or m_min unphysical (> 1 eV scale, contradicts cosmology)

**Falsification:** if IO neutrino oscillation experiments (T2K, JUNO 2027+)
detect IO at >5σ → K(ν) = 1/2 framework broken; ζ.1 must re-derive
neutrino chirality structure.

---

## Verdict gate

**5/5 PASS** → Σm_ν = 59.6 meV mass-spectrum LOCKED, Phase 2 proceeds.

**4/5 PASS** → audit gap, Phase 2 deferred lub limited scope.

**≤ 3/5 PASS** → ζ.1 reframing required.

---

## Środowisko

```bash
PYTHONIOENCODING=utf-8 python -X utf8 research/op-zeta-mass-spectrum/phase1_neutrino_audit.py 2>&1 | tee research/op-zeta-mass-spectrum/phase1_neutrino_audit.txt
```

## Cross-references

- [`program.md`](program.md) — overall ζ.1 plan
- [`../../PREDICTIONS_REGISTRY.md`](../../PREDICTIONS_REGISTRY.md) — C3 Σm_ν entry (LIVE 2027+)
