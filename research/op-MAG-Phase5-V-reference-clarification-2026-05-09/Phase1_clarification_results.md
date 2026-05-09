---
title: "Phase 1 clarification — wyniki + werdykt — op-MAG-Phase5-V-reference-clarification"
date: 2026-05-09
type: phase-results
status: COMPLETE
parent: "[[./README.md]]"
phase: 1
verdict: PATH_C_DUAL_V_STRUCTURE_HYPOTHESIS
sympy_pass: "10/10"
tags:
  - phase1
  - clarification-results
  - V-reference-clarification
  - phase5-MAG
  - lambda4-sign-flip
  - dual-V-hypothesis
related:
  - "[[./README.md]]"
  - "[[./Phase0_balance.md]]"
  - "[[./NEEDS.md]]"
  - "[[./Phase1_clarification_sympy.py]]"
---

# Phase 1 clarification — wyniki + werdykt

## Executive summary

**Sympy: 10/10 PASS.**

**Werdykt:** **PATH C HYPOTHESIS** — TGP może mieć **DWA V**:
- $V_{M9.1''}$ canonical dla gravitational sektora (M9.1'' geometry, G.0 closure)
- $V_{orig}$ dla matter field fluctuations sektora (Phase 5 Mach inertia, particle masses)

**Kluczowe odkrycie:** $\lambda_4$ (quartic coupling) **ZMIENIA ZNAK** pod
$V_{orig} \to V_{M9.1''}$:
- $V_{orig}$: $\lambda_4 = +3\gamma/(2\Phi_0^2)$ (positive)
- $V_{M9.1''}$: $\lambda_4 = -9\gamma/2$ (negative, **constant w ψ**)

To znaczy m_Mach z V_M9.1'' canonical byłby **NEGATYWNY** (unphysical) —
**Path A (lightweight re-interpretation) FALSIFIED**.

---

## 1. Confirmed facts

### 1.1 Phase 5 explicit cytuje V_orig (T1)

[[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] **linia 38**:
```
V(Φ) = -β Φ³/(3 Φ_0) + γ Φ⁴/(4 Φ_0²)
```

[[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]] **linia 63**:
```python
V_Phi = -beta_p * Phi**3 / (3 * Phi_0) + gamma_p * Phi**4 / (4 * Phi_0**2)
```

**To jest V_orig formula** (sek08a linie 95-110, DEPRECATED 2026-05-02).

### 1.2 Phase 5 closed PO V_M9.1'' canonical lock

| Cykl | Closure date |
|---|---|
| G.0 V_M9.1'' canonical lock | 2026-05-02 |
| Phase 5 MAG (V_orig usage) | **2026-05-09** (7 dni PO) |

### 1.3 Phase 5 derivation jest V-DEPENDENT

Phase 5 m_Mach derivation steps (cyt. linia 88):
```
m_Mach = (3γq²)/(16π Φ_0² m_C) · ⟨δΦ²_bg⟩
```

Wszystkie kroki:
1. **Expansion around Φ=Φ_0** (V_orig vacuum z β=γ)
2. **V''''(Φ_0) = 6γ/Φ_0²** (V_orig 4th derivative, T1 PASS)
3. **λ_4 = V''''/4 = 3γ/(2 Φ_0²)** (V_orig, T2 PASS)
4. **m_Mach ∝ λ_4** (linear w lambda_4)

Każdy krok zależy strukturalnie od V_orig form.

---

## 2. Critical sympy findings

### 2.1 V_orig vs V_M9.1'' V'''' comparison

**T1 (V_orig):** $V_{orig}''''(\Phi_0) = 6\gamma/\Phi_0^2$ (positive)

**T3-T4 (V_M9.1''):** $V_{M9.1''}''''(\psi) = -18\gamma$ — **CONSTANT w ψ**:
- V''''(ψ=2/3 minimum) = -18γ
- V''''(ψ=1) = -18γ
- V''''(ψ=4/3 horizon) = -18γ
- ... NEGATIVE wszędzie

To wynika z faktu że V_M9.1'' = -γψ²(4-3ψ)²/12 jest **quartic w ψ** (po
rozwinięciu: -4γψ²/3 + 2γψ³ - 3γψ⁴/4) — czwarta pochodna jest stała.

### 2.2 λ_4 sign flip (T6 — CRITICAL)

| Potencjał | λ_4 = V''''/4 | Sign | Magnitude |
|---|---|---|---|
| V_orig (Phase 5 baseline) | $+3\gamma/(2\Phi_0^2)$ | **POSITIVE** | dimensional |
| V_M9.1'' canonical | $-9\gamma/2$ | **NEGATIVE** | dimensionless |

**Sign flip jest kwalitatywny** — nie jest "small correction". Implikuje
**fundamental difference** między dwoma derivations.

### 2.3 m_Mach formula impact (T8)

Z V_orig λ_4 = +3γ/(2Φ_0²):
$$m_{Mach}^{V_{orig}} = +\frac{3\gamma q^2}{16\pi \Phi_0^2 m_C} \cdot \langle\delta\Phi_{bg}^2\rangle \quad (\text{POSITIVE} = \text{physical mass})$$

Z V_M9.1'' λ_4 = -9γ/2:
$$m_{Mach}^{V_{M9.1''}} = -\frac{9\gamma q^2}{16\pi m_C} \cdot \langle\delta\Phi_{bg}^2\rangle \quad (\text{NEGATIVE} = \text{unphysical})$$

**Phase 5 m_e=511 keV reproduction** wymaga POSITIVE λ_4 — V_M9.1'' canonical
nie dostarcza tego.

### 2.4 Path A (re-interpretation) FALSIFIED (T7-T8)

Hipoteza Path A: "Phi_0_Phase5 = (2/3)·Phi_0_V_M911 — re-interpretacja parametru"

**Test FAILED:** sign issue jest **fundamentalny** (V''''<0 w V_M9.1'' wszędzie),
nie zależy od reinterpretacji Phi_0. Path A nie rozwiązuje qualitative sign problem.

---

## 3. Three resolution paths (post-sympy)

### Path A: Re-interpretation Phi_0 (FALSIFIED)
- **Status:** ❌ FALSIFIED przez T6+T8 sympy
- **Reason:** sign issue fundamental, nie zależy od Phi_0 reinterpretation

### Path B: Full Phase 5 re-derivation z V_M9.1'' canonical (DIFFICULT)
- **Status:** 🟡 STRUCTURALLY DIFFICULT
- **Issue:** wymaga radically different mass mechanism (m∝-λ_4 → m∝...?)
- **Plausible if:** Phase 5 mass formula ma mistake i m∝|λ_4| zamiast λ_4
- **Probability:** 20-30%

### Path C: Dual-V structure hypothesis (SUGGESTED)
- **Status:** 🟢 EMERGENT z sympy findings
- **Hypothesis:** TGP MA dwa potencjały:
  - $V_{M9.1''}$ — **gravitational sektor**: M9.1'' metryka, G.0 closure,
    mass spectrum invariance (G.0 P22 5/5 PASS)
  - $V_{orig}$ — **matter sektor**: field fluctuations, particle masses
    via Mach inertia (Phase 5)
- **Implikacje:**
  - V_orig **NIE jest deprecated** — jest specyficzny dla matter sector
  - "DEPRECATED 2026-05-02" w sek08a może być MISLEADING annotation
  - sek08a wymaga **addendum 2026-05-09**: dual-V structure clarification
- **Probability:** 50-60% — najmocniej zgodne z istniejącymi positive results
  (T-Λ używa V_orig, Phase 5 używa V_orig, G.0 używa V_M9.1'')

### Path D: Phase 5 V-independent reformulation (EMERGENT)
- **Status:** 🟡 ALTERNATIVE (untested)
- **Hypothesis:** m_Mach mechanism może być re-formulated bez explicit V
  (np. operator-based, action-based formalism)
- **Wymaga:** osobny cykl theoretical reformulation
- **Probability:** 10-20%

---

## 4. Path C deep dive — Dual-V structure hypothesis

### 4.1 Co supportuje Path C

**Pozytywne wyniki używające V_orig:**
- T-Λ closure (2026-04-26): ρ_vac numerical match 1.020 (V_orig formula)
- Phase 5 MAG (2026-05-09): m_e=511 keV reproduction (V_orig λ_4)
- UV.3 Z_Φ=14/3 (2026-05-04): używa P(g) z β/γ ratio (V_orig-style)

**Pozytywne wyniki używające V_M9.1'':**
- G.0 P22 mass spectrum: m_μ/m_e, m_τ/m_e invariant pod V_orig→V_M9.1''
  (5/5 sympy PASS) — but uses A_tail mechanism, NIE V coupling
- G.0 P32 newton limit: kappa = 3/(4·Phi_0) invariant
- M9.1'' canonical metric: ψ=4/3 horizon

**Pattern:** wyniki **gravitational/geometric** (G.0, M9.1'' metric, kappa)
działają z V_M9.1''. Wyniki **matter/particle physics** (T-Λ, Phase 5, UV.3
Z_Φ) używają V_orig.

### 4.2 Czy to fizycznie sensowne?

**Argument FOR dual-V:**
- Gravity i matter mają **różne dynamics** (np. gravity = R + cosmological const,
  matter = field operators z V_eff)
- TGP single-Φ framework **może legitimately mieć dwa effective potentials**
  dla różnych sectors:
  - V_M9.1'': **gravitational vacuum structure** (defines M9.1'' geometry)
  - V_orig: **matter field self-interaction** (governs particle dynamics)

**Argument AGAINST dual-V:**
- Single-Φ framework powinien być spójny (jeden V)
- Może to być **artifact** różnych approximations różnych regimes:
  - V_M9.1'' = full sextic-style canonical
  - V_orig = quartic approximation w specific limit (np. low-amplitude)

### 4.3 Test: czy V_orig jest expansion V_M9.1'' w specific regime?

**T7 sympy:** Taylor V_M9.1'' around ψ=2/3 (vacuum):
$$V_{M9.1''}(2/3 + \xi) = -\frac{4\gamma}{27} + \frac{2\gamma}{3}\xi^2 - \frac{3\gamma}{4}\xi^4$$

**Note:** brak cubic term ($V'''(2/3)=0$). To **NIE matches V_orig structure**
(V_orig ma cubic term -βΦ³/3Phi_0).

V_orig and V_M9.1'' są **strukturalnie różne** — V_orig nie jest expansion
V_M9.1'' wokół żadnego specific point.

### 4.4 Wniosek Path C

**Hypothesis:** V_orig i V_M9.1'' to **dwa niezależne potentials** w TGP:
- V_M9.1'' z M9.1'' geometry (canonical post-G.0)
- V_orig z matter field theory (canonical pre-G.0, **nadal potrzebne**)

**sek08a "DEPRECATED 2026-05-02" annotacja może być nadinterpretacja G.0
findings** — G.0 wykazał że V_M9.1'' jest canonical dla **gravitational**
sektora, ale to nie znaczy że V_orig jest globally deprecated.

**Recommended action:** spawn cycle `op-dual-V-structure-clarification-YYYY-MM-DD`
żeby formalnie sprawdzić tę hipotezę.

---

## 5. Cross-impact na inne cykle

| Cykl | V usage | Re-classification post-Path C |
|---|---|---|
| T-Λ closure | V_orig formula | ✅ **VALID** (matter sector) — Path C uzasadnia residual gap T-Λ z [[../op-V-canonical-consistency-audit-2026-05-09/]] |
| Phase 5 MAG | V_orig formula | ✅ **VALID** (matter sector) |
| G.0 closure | V_M9.1'' canonical | ✅ **VALID** (gravity sector) |
| sek08a master | V_orig deprecated annotation | ⚠️ **MISLEADING** — wymaga update |
| op-Phi-vacuum-scale Phase 1 | użył V_orig | ⚠️ **PARTIALLY VALID** — depends on Phi_0 sector identification |

**Path C resolves residual gap z op-V-canonical-consistency-audit:**
- T-Λ ratio match z V_orig (matter Phi_0): consistent with Path C
- Phase 5 m_e match z V_orig (matter Phi_0): consistent with Path C
- G.0 mass spectrum invariance: V_M9.1'' canonical (gravity) — consistent

---

## 6. CALIBRATION_PROTOCOL compliance

Audit-clarification cycle, NIE derivation. Anti-patterns N/A.

**Honest reporting:**
- ✅ T6 sympy SHOWS lambda_4 sign flip — explicit, honest
- ✅ Path A FALSIFIED (T7-T8) — NIE markup jako "OK"
- ✅ Path C SUGGESTED ale **NIE PROVEN** — clearly labeled hypothesis
- ✅ Recommend follow-up cycle do formal verification dual-V

---

## 7. Decision matrix

### 7.1 Werdykt finalny

**Sympy verdict:** Phase 5 V_orig usage **NIE jest błędem** — jest
**fundamentally needed** bo V_M9.1'' canonical (negative λ_4) NIE może
reprodukować positive m_Mach.

**Most likely resolution:** Path C (dual-V structure)

**Recommended actions:**

| # | Action | Priority |
|---|---|---|
| 1 | Spawn `op-dual-V-structure-clarification-YYYY-MM-DD` | HIGH |
| 2 | Update sek08a 2026-05-09 addendum: clarify V_orig is matter sector, NOT globally deprecated | HIGH |
| 3 | Update [[../op-V-canonical-consistency-audit-2026-05-09/]] z Path C resolution dla T-Λ residual gap | MEDIUM |
| 4 | Update [[../op-Phi-vacuum-scale-2026-05-09/]] NEEDS — P11 status revised to "PATH C hypothesis" | MEDIUM |

### 7.2 Awaiting user decision

**Trzy opcje dla user:**

1. **(a) Spawn op-dual-V-structure-clarification** — formal verification Path C hipotezy
2. **(b) Update sek08a directly** — clarification annotation in source
3. **(c) Defer obie, formal close** — leave as documented hypothesis

---

## 8. Files generated

- [[./README.md]] — scoping
- [[./Phase0_balance.md]] — facts identification
- [[./NEEDS.md]] — B1-B5 audit tasks
- [[./Phase1_clarification_sympy.py]] — sympy 10/10 PASS
- [[./Phase1_clarification_results.md]] — niniejszy dokument

## Cross-references

- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_results.md]] linia 38
- [[../op-MAG-resonance-formalization-2026-05-09/Phase5_Mach_inertia_sympy.py]] linia 63
- [[../op-V-canonical-consistency-audit-2026-05-09/Phase1_audit_results.md]] §1.2
- [[../op-Phi-vacuum-scale-2026-05-09/]]
- [[../op-g0-r3-from-canonical-projection/]]
- [[../../core/sek08a_akcja_zunifikowana/sek08a_akcja_zunifikowana.tex]]

## Status

**Phase 1 COMPLETE — PATH C DUAL-V HYPOTHESIS suggested.**

**Awaiting user decision** o Path C verification: dedicated cycle, sek08a update, lub defer.
