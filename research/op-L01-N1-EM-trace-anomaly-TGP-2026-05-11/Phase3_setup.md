---
title: "Phase 3 setup — Phenomenology: lab + magnetar regimes + MICROSCOPE/Eöt-Wash + QEP universality check (R6)"
date: 2026-05-11
parent: "[[./README.md]]"
type: phase-setup
phase: 3
status: 🟡 setup phase
sub_needs_addressed: [N0.7, N0.8, N0.9]
risks_addressed: [R5 full, R6]
predecessor: "[[./Phase2_results.md]] (8/8 PASS)"
tags:
  - phase3
  - phenomenology
  - lab-regime
  - magnetar-regime
  - MICROSCOPE
  - QEP-universality
---

# Phase 3 setup

## §0 — Cel Phase 3

Numerical estimates `ρ_EM_quantum[{Φ_i}]` w obserwowanych regimes:
- Lab (B~1 T, R_curv~10⁻⁵² m⁻²)
- Magnetar (B~10¹¹ T, B/B_QED~23 — granica perturbative)

i verification że bounds (MICROSCOPE η ≤ 1.1·10⁻¹⁵, Eöt-Wash, LLR Nordtvedt,
WEP universality) są **automatically satisfied** przez universal coupling
structure (S05 mechanism).

Sub-needs: N0.7 + N0.8 + N0.9.

Risks: R5 (perturbative breakdown w extreme magnetar) honestly documented;
R6 (QEP violations from trace anomaly z Asorey-2015) **verified to be cancelled**
przez universal coupling.

## §1 — Lab regime numerics

### §1.1 — Lab field strengths

Typical lab conditions:
- Magnetic field: B_lab ~ 1 T (electromagnet) ÷ 50 T (pulsed superconducting magnet)
- Electric field: E_lab ~ 10⁶ V/m (lab capacitor) ÷ 10⁸ V/m (FEL)
- F² = 2(B² - E²/c²); for B=1 T, F² ≈ 2 T² × c² ≈ 2 (1)² × (3·10⁸)² ≈ 2·10¹⁷ T²·m²/s²

Konkretnie w SI units:
- B² / (2μ_0) = energy density of magnetic field
- For B=1 T: u_mag = B²/(2μ_0) = 1²/(2·4π·10⁻⁷) ≈ 4·10⁵ J/m³

Background curvature (cosmological + local):
- R_cosmo ~ H₀² ≈ 10⁻⁵² m⁻²
- Local Earth gravity: R_Earth ~ G·M_Earth/r³ ≈ 6.67·10⁻¹¹ · 6·10²⁴ / (6.4·10⁶)³ ≈ 1.5·10⁻⁶ m⁻² ... wait, that's m⁻²?

Let me redo. Schwarzschild scalar Ricci R = 0 outside source (vacuum solution).
R_μν = 0 in vacuum. The relevant curvature is **Riemann tensor** components,
not R. For Earth surface:
- Tidal acceleration ~ GM/r³ ~ 1.5·10⁻⁶ s⁻²
- In geometric units, R_μνρσ ~ tidal/(c²) ~ 1.7·10⁻²³ m⁻²

For our analysis, R_local ~ R_cosmo (cosmological) — Earth's contribution to scalar
curvature R is zero (vacuum). Phase 2 §3.1 used R_cosmo correctly.

### §1.2 — ρ_EM_quantum lab estimate

Per Phase 2 §1.4, dominant term:
```
T^μ_μ_EM,1-loop_TGP ≈ (α/(3π)) · F²
```

Pomijając curvature × F² corrections (suppressed by R/m_e² ~ 10⁻⁷⁷ in lab):

```
T^μ_μ_EM,quantum_lab ≈ (α/(3π)) · F² ~ 7.74·10⁻⁴ · F²
```

Numerical dla B=1 T, E=0:
```
F² = -2 B²/c² (Lorentzian signature mostly)... wait
```

Konwencja: w Heaviside-Lorentz units w (-,+,+,+) signature:
```
F_{μν} F^{μν} = 2(B² - E²/c²)
```

z B² in Tesla² (units of A²·kg/m·s⁻²... actually let's just use energy density).

**Practical formula:**
```
F_{μν}F^{μν}_inSI = 2 (B² - E²/c²) = -(2/c²)·(E² - c²B²) = (2μ_0)·(u_mag - u_elec)
```

Dla B=1 T, E=0: F² = 2 B² = 2 T² (in SI: J·s/(C·m·m/s²·s) hmm... let's just do energy density).

Lepiej: trace anomaly **density** (energy per volume):
```
T^μ_μ ~ α/(3π) · 2·(u_mag - u_elec)
```

dla B=1 T, E=0: u_mag = B²/(2μ_0) ≈ 4·10⁵ J/m³.
```
T^μ_μ_quantum_lab ≈ 7.74·10⁻⁴ · 2 · 4·10⁵ ≈ 6.2·10² J/m³
```

```
ρ_EM_quantum_lab = -T^μ_μ / c_0² ≈ -6.2·10² / (9·10¹⁶) ≈ -7·10⁻¹⁵ kg/m³
```

(Note: sign — negative because trace anomaly contributes *negative* trace; ρ
density positive in *magnitude*, but contribution to source has sign-convention.)

**Numerical magnitude:** |ρ_EM_quantum_lab| ~ **7·10⁻¹⁵ kg/m³** for B=1 T.

For comparison: vacuum baseline ρ_baryon (interstellar) ~ 10⁻²¹ kg/m³ — ρ_EM_quantum
is **6 OOM larger** w polu B=1 T! But: this is *local* density inside the magnet,
NIE *cosmological* contribution.

### §1.3 — Implikacje dla MICROSCOPE / 5-th force / WEP

**MICROSCOPE η bound** dotyczy *différentielle* free-fall acceleration dla dwóch
materials z różnym EM-content. Pt vs Ti — różnica EM-content (Pt ma więcej
elektronów per atom, więcej F² magnetic energy bound w atomic structure).

**Z B9 closure baseline** (PASS 6/6 2026-05-01):
```
η_TGP_Dirac (Pt vs Ti) = 1.32·10⁻²⁶ (z mass-coupling structure)
```

**Quantum trace anomaly contribution** dla typical atomic Coulomb field:
- Atomic scale F²: dla hydrogen, E_atomic ~ e/(4πε_0·a_0²) ≈ 5.1·10¹¹ V/m
- F²_atomic ≈ ε_0 E²_atomic / c² = (8.85·10⁻¹²) · (5.1·10¹¹)² / (3·10⁸)²
                   ≈ (8.85·10⁻¹²) · 2.6·10²³ / 9·10¹⁶ ≈ 2.6·10⁻⁵ J/m³
  (wait, that's energy density in atomic vicinity)
- Per atom (volume ~ a_0³ ~ 10⁻³⁰ m³): u_atom ~ 2.6·10⁻⁵ · 10⁻³⁰ ≈ 2.6·10⁻³⁵ J ≈ 10⁻¹⁶ eV

(Hmm too small — atomic Coulomb energy should be ~ Rydberg ~ 13.6 eV per atom.)

Let me redo. Hydrogen Rydberg = 13.6 eV. Energy density inside Bohr radius:
```
u_atom = E_Coulomb / V_atom ≈ 13.6 eV / (4/3 π a_0³) 
       = 13.6 · 1.6·10⁻¹⁹ J / (1.4·10⁻³⁰ m³)
       ≈ 1.5·10⁹ J/m³
```

That's the EM energy density inside an atom. Trace anomaly contribution:
```
T^μ_μ_quantum_atom ~ α/(3π) · 2 · u_atom ≈ 7.74·10⁻⁴ · 3·10⁹ ≈ 2.3·10⁶ J/m³
```

Per atom (volume 10⁻³⁰ m³):
```
δm_atom ~ T^μ_μ_quantum_atom · V_atom / c² ≈ 2.3·10⁶ · 10⁻³⁰ / 9·10¹⁶
        ≈ 2.6·10⁻⁴⁰ kg ≈ 1.5·10⁻¹¹ eV
```

This is **Lamb-shift-like correction** to atomic mass per atom. For atom z mass
~ 1 GeV ~ 10⁻²⁷ kg, ratio:
```
δm/m_atom ≈ 2.6·10⁻⁴⁰ / 10⁻²⁷ = 2.6·10⁻¹³
```

For Pt (Z=78) vs Ti (Z=22), differential δ(δm/m) — to powinno być w wymiarze
α·Z×F²/(c²·m_atom) różnicy. Order of magnitude estimation:

```
δη_EM_quantum ~ (Z_Pt² - Z_Ti²) · α/(3π) · u_atom · V_atom / (m_atom·c²)
             ≈ (78² - 22²) · 7.74·10⁻⁴ · 1.5·10⁹ · 10⁻³⁰ / (10⁻²⁷ · 9·10¹⁶)
             = (5600) · 7.74·10⁻⁴ · 1.5·10⁹·10⁻³⁰ / (9·10⁻¹¹)
             = 5600 · 7.74·10⁻⁴ · 1.5·10⁻²¹ / 9·10⁻¹¹
             = 5600 · 7.74·10⁻⁴ · 1.7·10⁻¹¹
             ≈ 7.4·10⁻¹¹
```

Hmm — to wygląda **NA POZIOMIE bound 10⁻¹⁵** lub powyżej. Let me sprawdzić to
ostrożnie.

**Reality check:** Lamb shift for hydrogen ~ 1 GHz ~ 4·10⁻⁶ eV. Vs Rydberg ~13.6
eV → ratio 3·10⁻⁷. To znaczy że quantum EM correction do atomic energy z Lamb
shift jest ~10⁻⁷, NIE 10⁻¹³. So my estimate jest **zbyt szczegółowy**.

The proper comparison: trace anomaly contribution do ρ_atom jest *implicit* w
Lamb shift, które *jest już part* of measured atomic mass. Czyli MICROSCOPE
constrains *total* η (Pt vs Ti) ≤ 10⁻¹⁵, włączając już Lamb-shift-type
contributions.

The relevant question dla TGP: czy **dodatkowy** Φ-mediated coupling generates
różnicę? Odpowiedź:

**Universal coupling structure** (S05 mechanism): `L_mat = -(q/Φ_0)·φ·ρ` gdzie ρ
zawiera *wszystkie* contributions (rest mass + EM binding + Lamb shift +
trace anomaly). Coupling do Φ jest **uniwersalny** — taki sam dla wszystkich
ρ-componentów z atomu.

⇒ Pt i Ti mają *różne* ρ_atom (różne Z, różne EM binding), ale ich coupling do
Φ jest *proporcjonalny* do ρ_atom × q/Φ_0. **Nie generuje** differential
acceleration.

⇒ η_TGP_EM_quantum **automatically zero** structurally, NIE nontrivially
constrained przez MICROSCOPE.

To jest dokładnie **R6 closure** — universal coupling kasuje QEP violations
(Asorey-2015 type non-local Lagrangians dają QEP violations w **non-universal**
coupling theories; w TGP universal coupling z S05 jest immune).

### §1.4 — Sympy verification target

Phase3_sympy.py weryfikuje:

1. Lab regime ρ_EM_quantum estimate (B=1 T): magnitude check.
2. Atomic-scale ρ_EM_quantum (Pt vs Ti): per-atom estimate.
3. **R6 closure:** universal coupling structure → η_TGP_EM_quantum *exactly zero*
   strukturalnie (NIE numerically suppressed).
4. MICROSCOPE η bound check: η_TGP_EM_quantum ≤ 10⁻¹⁵ automatically.
5. Eöt-Wash bound check: similar.
6. LLR Nordtvedt η bound (4.4·10⁻⁴ for moon): automatically.
7. **R5 magnetar**: B ~ 10¹¹ T → B/B_QED ≈ 23, perturbative MARGINAL/BREAKS;
   estimate w B ≪ B_QED limit (linear extrapolation) z honest caveat.
8. Cross-cycle τ.3 mechanism decoupling: confirmed via independent estimate.

## §2 — Magnetar regime numerics

### §2.1 — Magnetar field strengths

Typical magnetar: B ~ 10⁹ to 10¹¹ T (surface). B_QED ≈ 4.4·10⁹ T (Schwinger).

Ratio B/B_QED:
- B = 10⁹ T → B/B_QED ≈ 0.23 (perturbative valid)
- B = 4.4·10⁹ T = B_QED (Schwinger boundary)
- B = 10¹⁰ T → B/B_QED ≈ 2.3 (perturbative starts to break)
- B = 10¹¹ T → B/B_QED ≈ 23 (deeply non-perturbative)

**R5 documentation:** dla B ≳ B_QED, perturbative QED expansion fails;
non-perturbative analysis (including pair production, Schwinger effect, etc.)
required. **Tego cyklu zasięg: B ≪ B_QED.**

### §2.2 — Estimate dla typical magnetar B ~ 10¹⁰ T

Linear extrapolation z perturbative formula:
```
T^μ_μ_quantum ≈ (α/(3π)) · 2 · u_mag = (α/(3π)) · B²/μ_0
```

For B = 10¹⁰ T:
```
u_mag = (10¹⁰)² / (2·4π·10⁻⁷) = 10²⁰ / (2.51·10⁻⁶) ≈ 4·10²⁵ J/m³
T^μ_μ_quantum ≈ 7.74·10⁻⁴ · 8·10²⁵ ≈ 6·10²² J/m³
ρ_EM_quantum = -T^μ_μ/c² ≈ -6·10²² / (9·10¹⁶) ≈ -7·10⁵ kg/m³
```

**Magnitude:** ρ_EM_quantum ~ **7·10⁵ kg/m³**.

Typowy magnetar surface ρ_NS ~ 4·10¹⁷ kg/m³ (neutron-star surface).

**Ratio:** ρ_EM_quantum / ρ_NS ~ 7·10⁵ / 4·10¹⁷ = **1.7·10⁻¹²**.

**Hmm.** That's compatible with L01 ADDENDUM Q3 estimate "10⁻¹²" — but
wait, L01 ADDENDUM had calculation:
```
B²/(2μ_0) ~ 4·10²⁸ J/m³ for B=10¹¹ T
α²/(3π) ≈ 7.7·10⁻⁷ (TYPO; actually α/(3π) ≈ 7.7·10⁻⁴)
```

L01 ADDENDUM used B=10¹¹ T (I'm using 10¹⁰ T) AND used `α²/(3π)≈10⁻⁷` (typo).
Correct estimate dla L01 ADDENDUM B=10¹¹ T z **correct** α/(3π)≈7.74·10⁻⁴:

```
u_mag (B=10¹¹) = (10¹¹)² / (2·4π·10⁻⁷) = 10²² / 2.51·10⁻⁶ ≈ 4·10²⁷ J/m³
T^μ_μ_quantum ≈ 7.74·10⁻⁴ · 2 · 4·10²⁷ ≈ 6·10²⁴ J/m³
ρ_EM_quantum ≈ -6·10²⁴ / 9·10¹⁶ ≈ -7·10⁷ kg/m³
ρ_EM_quantum/ρ_NS ~ 7·10⁷ / 4·10¹⁷ ≈ 1.7·10⁻¹⁰
```

**Corrected ratio for B=10¹¹ T (extreme magnetar):** ~10⁻¹⁰ (z corrected α/(3π)
prefactor).

To jest **2 OOM większe** niż L01 ADDENDUM cytuje (10⁻¹²) — z **correction
z typo**, ale dalej ≪ 1 ⇒ ρ_EM_quantum *nadal nieistotne* dla magnetar
gravitational dynamics.

### §2.3 — TT10 magnetar polar shift status

Per τ.3 ADDENDUM §2: TT10 magnetar X-ray timing testuje **L4 gradient-coupled
mass mechanism**, NIE quantum trace anomaly.

Z corrected estimate ρ_EM_quantum/ρ_NS ~ 10⁻¹⁰ (for B=10¹¹ T):
- Polar Φ-shift od ρ_EM_quantum: **2 OOM więcej niż L01 ADDENDUM**, ale wciąż
  ≪ surface ρ_NS contribution.
- TT10 *insensitive* na ρ_EM_quantum w typical magnetar — **τ.3 mechanism decoupling
  preserved** (10 OOM separation vs 8 OOM separation; wciąż mechanism decoupled).

**Action item Phase 4:** propagate corrected numerics do τ.3 ADDENDUM §2 (L01
ADDENDUM §3.2 typo correction).

### §2.4 — When ρ_EM_quantum becomes significant

Linear extrapolation breaks at B ≳ B_QED. Hypothetical *Schwinger-class lab*
(B ~ 10⁹ T w macroscopic volume — beyond 2030 zasięg):
- u_mag(B=10⁹) ~ 4·10²³ J/m³
- T^μ_μ_quantum ~ 7.74·10⁻⁴ · 8·10²³ ≈ 6·10²⁰ J/m³
- ρ_EM_quantum ~ 7·10³ kg/m³ — comparable do iron solid mass density.

**To jest scenario gdzie quantum trace anomaly Φ-shift staje się observable.**
Currently outside lab capability, ale future Schwinger-class lab byłby concrete
falsifier.

## §3 — R6 closure: QEP universality verification

### §3.1 — Asorey-2015 QEP violations argument

JHEP 05 (2015) 118 (Asorey, Bautista et al.): "QED trace anomaly, non-local
Lagrangians and quantum equivalence principle violations":
- W non-universal coupling theories, QED trace anomaly z non-local Lagrangians
  can induce *quantum equivalence principle* (QEP) violations.
- Constraint: MICROSCOPE η ≤ 10⁻¹⁵ is sensitive na takie violations.

### §3.2 — Why TGP is immune

**S05 + universal coupling structure:**

W TGP, action `L_mat = -(q/Φ_0)·φ·ρ` ma *uniwersalną strukturę*:
- ρ ≡ -T^μ_μ/c_0² zbiera **wszystkie** contributions z stress-energy:
  rest mass + EM binding + nuclear binding + Lamb shifts + quantum corrections.
- Coupling to Φ jest **proporcjonalne do total ρ** — nie ma sektor-specific
  prefactors że "EM ρ couples differently than rest-mass ρ".

⇒ Wszystkie atomy (independent of EM-content) widzą *ten sam* gradient Φ przez
*ten sam* universal coupling. **Nie ma differential acceleration.**

W Asorey-2015 framework: QEP violations require *non-universal* coupling z EM
sektorze. TGP ma *strict universal coupling* z S05 + ax:metric-coupling — nie
wprowadza takiego niewyraziktego sektor-specific coupling.

### §3.3 — Strukturalna konsekwencja

```
η_TGP_EM_quantum (Pt vs Ti) = 0   (structurally)
```

NIE 10⁻²⁶ jak η_TGP_Dirac z B9 — tutaj **exactly zero**, bo nie ma differential
coupling. Phase 3 sympy verification dostaje to explicit.

**R6 zamknięty strukturalnie.**

## §4 — Phase 3 deliverables

- [[Phase3_setup.md]] (this file)
- [[Phase3_results.md]] — numerical estimates lab + magnetar + MICROSCOPE/Eöt-Wash/LLR
  bound checks + R6 closure
- *Note: Phase 3 jest primarily numerical/dimensional analysis; sympy here is
  back-of-envelope verification, less formal niż Phase 1-2.*

## §5 — Cross-references

- [[./README.md]]
- [[./Phase0_balance.md]] §3
- [[./Phase1_results.md]] (β/(2α) prefactor LOCK)
- [[./Phase2_results.md]] (operator structure + GW170817)
- [[../op-L01-rho-stress-energy-bridge-2026-05-04/ADDENDUM_2026-05-10_native_observables_first.md]] §3.2 Q3 (typo source)
- [[../op-newton-momentum/B9_wep_microscope_composition_results.md]] (η_TGP_Dirac baseline)
- [[../op-tau3-substrate-clock-acceleration/ADDENDUM_2026-05-10_native_observables_first.md]] §2 (mechanism decoupling)
- JHEP 05 (2015) 118 (Asorey et al.) — QEP risk-input

---

**Phase 3 setup ready.** Proceeds with numerical analysis + R6 closure documentation.
