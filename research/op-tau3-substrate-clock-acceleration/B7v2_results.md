---
title: "τ.3 audit B7-v2 closure — geometric edge analysis + realistic detector volume + corrected unit conversions"
date: 2026-05-01
cycle: τ.3
status: B7-v2-CLOSED-NUMERICAL
parent: "[[program.md]]"
audit_link: "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
predecessor: "[[B7_greens_function_results.md]]"
tags:
  - TGP
  - tau3
  - audit-B7-v2-closure
  - geometric-edge-analysis
  - detector-volume-integration
  - TT7-TT12-revised
  - omega3-OPEN
---

> **B7-v2 closure NUMERICAL** (post-B7 STRUCTURAL): pełna re-derivation
> TT7-TT12 numerical thresholds z (i) corrected unit conversions, (ii)
> geometric edge analysis cylindrical field region, (iii) realistic detector
> volume integration. **KEY FINDING: B7 unit conversion bug** (E_to_GeV²
> off by ~10⁶) → wszystkie B7 numerical δω/ω były zawyżone o ~10¹². **Best-case
> lab feasibility (ELI-NP + edge + Λ=10 MeV)**: δω/ω ≈ 1.45·10⁻⁴¹ (22 OOM
> poniżej Sr/Yb precision 10⁻¹⁹). **τ.3 lab detection NIEMOŻLIWA z
> technologią 2030+** w default parameters; cykl wymaga **ω.3 super-light
> substrate pivot** lub astrofizycznej relocation.

# τ.3 audit B7-v2 closure — geometric edge + realistic volume + corrected units

## B7-v2 audit context

B7 STRUCTURAL CLOSURE (2026-05-01) ustanowiła:
- Sympy LOCK two-regime formulae: heavy J²/m_X² vs light J²L²/(16π²)
- KEY PHYSICS: default τ.3 parameters w heavy regime, bulk = 0, edge-only

B7 wymagała **dedicated B7-v2 cycle** dla:
- Geometric edge analysis (precise (∂lnX) profile across cylindrical boundary)
- Realistic detector volume integration (V_clock-averaging)
- Numerical TT7-TT12 thresholds re-derivation
- Comparison z Sr/Yb optical clock precision 2030+ projections

## Critical correction: B7 unit conversion bug

**Diagnoza**: B7 script (`B7_greens_function_derivation.py`) used
`E_to_GeV2 = 1.96·10⁻¹⁹` for converting electric field V/m → natural GeV²
units. **Correct value**: `6.5·10⁻²⁵` (off by factor ~3·10⁵).

**Derivation of correct E_to_GeV²**:
```
e·E_Schwinger = m_e²       (definition of critical Schwinger field)
m_e² = (0.511 MeV)² = 2.61·10⁻⁷ GeV²
e_natural = √(4π·α_em) ≈ 0.303
E_Schwinger_natural = m_e²/e_natural = 8.6·10⁻⁷ GeV²
E_Schwinger_SI = 1.32·10¹⁸ V/m
→ 1 V/m → 8.6·10⁻⁷ / 1.32·10¹⁸ = 6.5·10⁻²⁵ GeV²
```

**Sanity check (B7-v2 script)**: `E_Schwinger_SI × E_to_GeV² = 8.58e-7 GeV²`
matches `m_e²/e_natural = 8.62e-7 GeV²` within 0.5%. ✓

**Konsekwencja**: wszystkie B7 numerical δω/ω były **zawyżone o ~10¹²**
(skala (E_to_GeV²)² = 10⁶ × 10⁶ = 10¹² w J² ∝ E²). B7-v2 corrects this.

B-conversion (`B_to_GeV² = 1.95·10⁻¹⁶`) była już poprawna.

## Geometric edge analysis (heavy regime)

### Edge profile derivation

Dla heavy regime (m_X·L >> 1) z cylindrical field region (radius R, length L):
- Bulk interior (distance d od edge >> 1/m_X): lnX(d) → -J/m_X² (uniform)
- Near edge (0 ≤ d ≤ few/m_X): lnX(d) ≈ -(J/m_X²)·(1 - e^{-m_X·d})
- Outside (d < 0): lnX → 0

**Edge gradient**:
$$|\partial \ln X|(d) = \frac{J}{m_X} e^{-m_X d}$$

**Edge gradient squared**:
$$(\partial \ln X)^2(d) = \frac{J^2}{m_X^2} e^{-2 m_X d}$$

### Spatial integral over edge shell

$$\int_V (\partial \ln X)^2 dV = A_{\text{edge}} \cdot \int_0^\infty \frac{J^2}{m_X^2} e^{-2 m_X d}\,dd = \frac{J^2 \cdot A_{\text{edge}}}{2\,m_X^3}$$

Cylindrical surface: `A_edge = 2π R(L+R)` (side + 2 caps).

### Volume averaging over clock region V_clock

$$\langle (\partial \ln X)^2 \rangle_{\text{avg}} = \frac{1}{V_{\text{clock}}} \int_V (\partial \ln X)^2 dV = \frac{J^2 \cdot A_{\text{edge}}}{2\,m_X^3\,V_{\text{clock}}}$$

### Suppression factor vs naive bulk

**Naive bulk estimate** (uniform J²/m_X² assumed): `(∂lnX)²_naive = J²/m_X²`

**Real edge-integrated** (heavy regime):
$$\frac{\langle (\partial \ln X)^2 \rangle_{\text{real}}}{(\partial \ln X)^2_{\text{naive}}} = \frac{A_{\text{edge}}}{2\,m_X\,V_{\text{clock}}}$$

For cubic V_clock = V_field = πR²L z R = L:
$$\boxed{\text{Suppression}_{\text{cubic}} = \frac{2}{m_X \cdot L}}$$

For default τ.3 (m_X·L = 4.21·10⁹): **suppression ≈ 4.7·10⁻¹⁰**.

## Three detector geometries

| Geometry | V_clock | Description | Practical |
|---|---|---|---|
| **(a) V_clock = V_field** | πR²L | Atom cloud fills entire field region | Realistic dla cold atom traps |
| **(b) V_clock = 100·V_field** | 100·πR²L | Diluted clock (atoms span larger volume) | Less optimal |
| **(c) V_clock = V_edge** | A_edge/(2m_X) | Atoms positioned exactly at edge shell | **Niepraktyczne** (sub-fm precision wymagana) |

## Numerical TT7-TT12 re-derivation

### Field schedules (B12-aware)

| Schedule | E [V/m] | B [T] | L [m] | m_X·L | Regime |
|---|---|---|---|---|---|
| (i) Schwinger IDEAL [B12-flagged niefizyczne] | 10¹⁵ | 100 | 10⁻³ | 4.21·10⁹ | HEAVY |
| (ii) ELI-NP routine [B12-rec] | 10¹³ | 30 | 10⁻³ | 4.21·10⁹ | HEAVY |
| (iii) Magnetar polar (SGR 1806-20) | 10¹⁰ | 2·10¹¹ | 10⁴ | 4.21·10¹⁶ | HEAVY |
| (iv) Cosmological PMF | 0 | 10⁻¹³ | 3.1·10²³ | 1.30·10³⁶ | E·B=0 NULL |

### δω/ω results (B7-v2 corrected)

**Best-case w każdej geometrii** (Λ = 10 MeV):

| Schedule | (a) V_clock = V_field | (b) V_clock = 100×V_field | (c) V_clock = V_edge |
|---|---|---|---|
| Schwinger IDEAL | 7.6·10⁻⁴⁶ | 7.6·10⁻⁴⁸ | 1.6·10⁻³⁶ |
| ELI-NP routine | 6.9·10⁻⁵¹ | 6.9·10⁻⁵³ | 1.4·10⁻⁴¹ |
| Magnetar polar | 3.1·10⁻⁴⁴ | 3.1·10⁻⁴⁶ | 6.4·10⁻²⁸ |
| Cosmological PMF | 0 | 0 | 0 |

**Detection threshold (Sr/Yb optical clock 2030+ projection)**: ~10⁻¹⁹

### Falsifiability gap analysis

**Best-case lab-feasible** (ELI-NP routine + edge-positioning + Λ=10 MeV):
δω/ω ≈ **1.45·10⁻⁴¹** vs precision 10⁻¹⁹ → **22 OOM gap**.

**Best-case Schwinger-ideal** (NIEFIZYCZNE per B12) + edge + Λ=10 MeV:
δω/ω ≈ 1.6·10⁻³⁶ → **17 OOM gap** (mimo niefizycznego E·B).

**Best-case astrophysical** (Magnetar polar + edge + Λ=10 MeV):
δω/ω ≈ 6.4·10⁻²⁸ → **9 OOM gap** od atomic clock; **ale**: lab Sr/Yb nie ma
dostępu do magnetar fields. Wymaga astronomical signatures (X-ray spectroscopy
pulsar magnetospheres, atomic transitions w extreme fields).

**Werdykt**: τ.3 default parameter regime (g_ω.1 = 8.3·10⁻³, f_X = 100 MeV)
+ realistic lab fields **NIEMOŻLIWE do detekcji** z technologią 2030+.

## ω.3 light substrate sector — alternative pivot

### Light regime detection threshold

Light regime wymaga `m_X·L_lab << 1`:
- L_lab = 1 mm → m_X << 1/L_lab = 1.97·10⁻¹³ GeV = **0.2 µeV**
- WW8-anchored g = 8.3·10⁻³ → f_X << **24 µeV**

**Default τ.3 z f_X = 100 MeV jest 4·10⁹× za ciężki** dla light regime w lab.

### Light regime δω/ω (hypothesis: f_X = 24 µeV, m_X·L = 0.1)

Z bulk Coulomb formula `(∂lnX)² = J²L²/(16π²)`:

| Schedule | (∂lnX)² [GeV²] | Λ=10 MeV δω/ω |
|---|---|---|
| Schwinger IDEAL | 5.7·10¹⁹ | 5.7·10²³ (perturbative breakdown!) |
| ELI-NP routine | 5.1·10¹⁴ | 5.1·10¹⁸ (perturbative breakdown!) |

**Status**: Light regime sygnał gigantyczny → perturbative formula breaks
down → wymaga non-perturbative re-derivation w przypadku exotic super-light
substrate. **ω.3 OPEN dedicated cycle** required dla:
1. UV-derivation super-light f_X (mechanism dla 4·10⁹ hierarchy z TGP scale)
2. Non-perturbative δω/ω formula w bulk Coulomb regime
3. Constraint compatibility z 5th-force experiments (Eöt-Wash, MICROSCOPE)
   z light scalar mediator m << 0.2 µeV (= 1 mm⁻¹)

### Path forward — three options

| Option | Description | Pros | Cons |
|---|---|---|---|
| **(A) ω.3 super-light pivot** | f_X << 24 µeV exotic substrate | Light regime → bulk signal lab-detectable | Wymaga UV-mechanism dla 4·10⁹ hierarchy; 5th-force constraints critical |
| **(B) Edge-engineered geometry** | Sub-fm clock positioning na boundary | Heavy regime intact (default params) | **Niepraktyczne** (atomic localization 10 nm bez fundamental QM precision limits) |
| **(C) Astrophysical relocation** | Magnetar/pulsar/AGN signatures | Natural strong fields | NIE jest "lab clock" experiment; wymaga astronomical X-ray spectroscopy |

**Rekomendacja**: **Option (C) astrophysical** jako primary detection avenue
post-B7-v2. **Option (A) ω.3** jako structural extension TGP dla potencjalnej
lab-feasibility, ale wymaga UV-mechanism + 5th-force constraints addressed.

## TT7-TT12 status update post-B7-v2

| Prediction | B7-v2 status | Revised value | Notes |
|---|---|---|---|
| **TT7** (Schwinger ideal δω/ω) | REVISED | ~10⁻³⁶ | NIEFIZYCZNE per B12 |
| **TT8** (ELI-NP routine δω/ω) | REVISED | ~10⁻⁴¹ | UNDETECTABLE 2030+ |
| **TT9** (Λ-scan threshold @ 1 GeV) | REVISED | ~10⁻⁴⁵ | UNDETECTABLE |
| **TT10** (Λ-scan threshold @ 100 MeV) | REVISED | ~10⁻⁴³ | UNDETECTABLE |
| **TT11** (Λ-scan threshold @ 10 MeV) | REVISED | ~10⁻⁴¹ | UNDETECTABLE |
| **TT12** (Sagnac differential bound) | REVISED | NULL lab-default | Wymaga astrophysical |

**Status update**: Wszystkie TT7-TT12 lab predictions w default τ.3 parameters
**effectively NULL** (22 OOM poniżej precision). Cykl τ.3 wymaga structural
reconfiguration:
- **Numerical predictions reformulated** w kontekście astrofizycznym lub
  ω.3 super-light substrate.
- **Lab feasibility downgraded** z "challenging 2030+" do "infeasible without
  ω.3 cycle".

## Cross-channel implications

### ψ.1 (substrate-light acceleration)

ψ.1 inherits same B7 finding (heavy regime → bulk = 0, edge only). Implications:
- ψ.1 lab Sagnac predictions analogously revised downward
- ψ.1 c-variation observable wymaga edge-engineering lub astrophysical
- **Dedicated ψ.1-v3 cycle** post-B7-v2 pending

### ω.3 (axion decay constant)

B7-v2 motivuje **ω.3 OPEN dedicated cycle** z primary goal:
- Derivation super-light f_X consistent z TGP fundamentals
- 5th-force constraint compatibility (light scalar mediator m << µeV)
- Hierarchy mechanism (4·10⁹ ratio z f_X = 100 MeV → 24 µeV)

### τ.2 NULL prediction

τ.2 NULL pozostaje structurally consistent (no engineered E·B → no signal).
B7-v2 nie zmienia τ.2 status.

## Audit B7-v2 closure verdict

| element B7-v2 | status |
|---|---|
| Unit conversion bug fix (E_to_GeV²: 1.96e-19 → 6.5e-25) | **CORRECTED** |
| Geometric edge analysis sympy LOCK | **CLOSED** (J²·A_edge/(2 m_X³ V_clock) formula) |
| Realistic detector volume integration | **CLOSED** (3 geometries × 4 schedules × 3 Λ) |
| TT7-TT12 numerical re-derivation | **CLOSED** (revised values 10⁻³⁶...10⁻⁴⁵) |
| Falsifiability gap analysis | **CLOSED** (22 OOM gap dla ELI-NP + edge) |
| ω.3 light substrate path forward | **DOCUMENTED** (dedicated cycle pending) |
| Astrophysical alternative (magnetar) | **DOCUMENTED** (9 OOM gap, X-ray spectroscopy) |

→ **B7-v2 NUMERICAL CLOSED**. TT7-TT12 lab predictions effectively NULL
w default τ.3 parameters; cykl wymaga **ω.3 OPEN dedicated cycle** lub
**astrophysical relocation** dla maintenance falsifiability.

## Files

- `B7v2_edge_geometry_analysis.py` — sympy edge integrals + numerical 4×3×3 grid (uruchomiony ✓)
- `B7v2_results.md` — niniejszy plik (results documentation + KEY PHYSICS update)

## Cross-references

- [[B7_greens_function_results.md]] — B7 STRUCTURAL closure (predecessor)
- [[B7_greens_function_derivation.py]] — B7 sympy two-regime LOCK
- [[Phase3_results.md]] — TT7-TT12 numerical predictions (post-B7-v2 reconfig pending)
- [[../../meta/AUDYT_TGP_2026-05-01.md]] § O.6 (B7 STRUCTURAL) + § T (B7-v2 NUMERICAL closure)
- [[../../PREDICTIONS_REGISTRY.md]] TT7-TT12 entries (post-B7-v2 ⚠ markers pending)
