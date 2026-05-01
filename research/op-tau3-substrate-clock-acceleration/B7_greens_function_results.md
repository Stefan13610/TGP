---
title: "τ.3 audit B7 closure — explicit (∂lnX)² Greens function derivation"
date: 2026-05-01
cycle: τ.3
status: B7-CLOSED-STRUCTURAL
parent: "[[program.md]]"
audit_link: "[[../../meta/AUDYT_TGP_2026-05-01.md]]"
tags:
  - TGP
  - tau3
  - audit-B7-closure
  - Greens-function
  - omega1-EOM
  - Schwinger
---

> **B7 closure structural** (post-A5 multiplicative δω/ω): explicit Yukawa Greens
> function derivation of (∂lnX)² for 2 regimes (heavy/light substrate) +
> sympy-LOCKED two-regime formulae + numerical evaluation dla 4 field schedules.
> KEY PHYSICS finding: heavy-substrate regime (m_X·L >> 1, **default lab parameters
> g_ω.1 = 8.3e-3 + f_X = 100 MeV → m_X = 0.83 MeV → 1/m_X ≈ 240 fm** << L_lab ≈ 1 mm)
> daje BULK signal essentially ZERO — only screened edge shell of thickness 1/m_X
> contributes. **Lab clock immersed inside field volume otrzymuje wprost ZERO signal
> w heavy regime** unless geometrically positioned at field boundary.

# τ.3 audit B7 closure — explicit (∂lnX)² Greens function derivation

## B7 audit context

Audit B7 (HIGH severity, 2026-05-01):

> "(∂lnX)² nigdzie nie obliczone z ω.1 EOM × Schwinger E·B; 'stipulated'
>  ε~10⁻¹² dziedziczone między ψ.1 a τ.3"

This note documents the explicit derivation requested.

## Setup

ω.1 EOM (Lorenz gauge, post-W2.5):
$$\Box(\ln X) = -\frac{g}{f_X^2}\,\mathbf{E}\cdot\mathbf{B}$$

For massive substrate (m_X = g·f_X effective mass):
$$(\Box + m_X^2)(\ln X) = -\frac{g}{f_X^2}\,\mathbf{E}\cdot\mathbf{B}$$

Static field configuration:
$$(\nabla^2 - m_X^2)(\ln X) = -J(\mathbf{x}), \quad J = \frac{g\,\mathbf{E}\cdot\mathbf{B}}{f_X^2}$$

Yukawa Greens function:
$$G(r) = -\frac{e^{-m_X r}}{4\pi r}$$

Solution:
$$\ln X(\mathbf{x}) = -\int G(\mathbf{x} - \mathbf{x}')\,J(\mathbf{x}')\,d^3x'$$

## Two-regime sympy LOCK

### Heavy regime (m_X·L >> 1) — screened, edge-localized

For uniform source J inside cylindrical region (length L, m_X·L >> 1):
- **Inside bulk**: ln X(x) → -J/m_X² (uniform Yukawa-equilibrium value)
- **∂lnX inside bulk = 0** (uniform → no gradient)
- **Edge shell** (thickness ~ 1/m_X): drop from J/m_X² → 0
- **|∂lnX|_edge ~ J/m_X**

**Sympy LOCK heavy regime**:
$$\boxed{(\partial \ln X)^2_{\text{edge}} = \frac{J^2}{m_X^2} = \frac{g^2\,E^2 B^2}{f_X^4\,m_X^2}}$$

**KEY PHYSICS**: Lab clock immersed in bulk → ZERO signal (uniform lnX inside).
Edge-localized signal only — geometric detector positioning required.

### Light regime (m_X·L << 1) — Coulomb-like, bulk signal

For uniform source J in characteristic volume L³, m_X·L << 1:
- Greens reduces to Coulomb: G(r) → -1/(4π r)
- Center value: ln X(0) ~ -J·L²/(4π)
- Bulk gradient: |∂lnX| ~ J·L/(4π)

**Sympy LOCK light regime**:
$$\boxed{(\partial \ln X)^2_{\text{bulk}} = \frac{J^2\,L^2}{16\pi^2} = \frac{g^2\,E^2 B^2\,L^2}{16\pi^2\,f_X^4}}$$

**KEY PHYSICS**: Lab clock anywhere inside bulk → uniform gradient signal.
Detector positioning trivial.

## Post-A5 multiplicative δω/ω formulae

Post-A5 patched: δω/ω = (α_g/Λ²)(∂lnX)² (no 1/m_e factor).

### Heavy regime
$$\delta \omega/\omega_{\text{edge}} = \frac{\alpha_g\,g^2\,E^2 B^2}{f_X^4\,m_X^2\,\Lambda^2}$$

### Light regime
$$\delta \omega/\omega_{\text{bulk}} = \frac{\alpha_g\,g^2\,E^2 B^2\,L^2}{16\pi^2\,f_X^4\,\Lambda^2}$$

## Default τ.3 parameter regime

With WW8 winner anchor `g_ω.1 = α_em·E_TGP/(2π) ≈ 8.3·10⁻³` + `f_X = 100 MeV`:

$$m_X = g \cdot f_X = 8.3\cdot 10^{-3} \times 100\,\text{MeV} = 0.83\,\text{MeV}$$
$$1/m_X \approx 240\,\text{fm}$$

For lab L ~ 1 mm = 10¹² fm: **m_X·L ~ 4·10⁹** >> 1 → **HEAVY REGIME**.

For light regime in lab, would need either:
- f_X << 1 eV (super-light substrate, exotic), OR
- g << 10⁻⁹ (dramatic suppression of ω.1 coupling, conflicts WW8 anchor)

→ **Lab τ.3 default parameters are heavy-regime; bulk signal ZERO; only edge.**

## Field schedule evaluations (sympy two-regime)

| Schedule | E [V/m] | B [T] | L [m] | Regime | (∂lnX)²_relevant | Notes |
|---|---|---|---|---|---|---|
| (i) Schwinger IDEAL [B12-flagged] | 10¹⁵ | 100 | 10⁻³ | HEAVY | edge only | niefizyczne razem |
| (ii) ELI-NP routine [B12-rec] | 10¹³ | 30 | 10⁻³ | HEAVY | edge only | realistic 2030+ |
| (iii) Magnetar polar (SGR 1806-20) | 10¹⁰ | 2·10¹¹ | 10⁴ | HEAVY (m_X·L=4·10¹⁶) | edge only | astrophysical |
| (iv) Cosmological PMF (1 nG, 10 Mpc) | 0 | 10⁻¹³ | 3·10²³ | HEAVY (E·B=0) | trivial 0 | NULL ✓ |

**All 4 schedules HEAVY regime** with default TGP parameters. None give bulk signal.

## Physical interpretation

Audit B7 reveals critical structural finding:

1. **Default τ.3 numerical predictions TT7-TT12 (pre-A5) assumed bulk signal**
   without recognizing m_X·L >> 1 heavy regime → screening of bulk ∂lnX.

2. **Lab E·B engineering chain operates in heavy regime** with WW8-anchored g, f_X.
   Bulk signal ZERO; only EDGE shell of thickness 1/m_X ~ 240 fm contributes.

3. **Detection geometry constraint**: Lab clock must be positioned at field
   boundary (within ~240 fm of E·B region edge), NOT in bulk interior.
   This is **geometrically extreme** — requires sub-fm clock positioning in a
   field region with sharp edges.

4. **Realistic τ.3 lab observability**: heavy-regime edge-only signal,
   integrated over realistic clock volume (~mm³), has ~ (1/m_X / L_clock)
   geometric suppression ~ 10⁻¹² × bulk-regime estimate.

5. **Alternative pivot**: τ.3 mechanism may require **light substrate sector**
   (m_X << 1/L_lab → super-light X with f_X·g extremely small) — degenerate
   z m_field bound z ω.1 (light substrate degenerated z α_em).

## Audit B7 closure verdict

| element B7 | status |
|---|---|
| sympy LOCK two-regime formulae | **CLOSED** (heavy: J²/m_X²; light: J²L²/(16π²)) |
| ω.1 EOM × Schwinger E·B Greens function explicit | **CLOSED** (Yukawa kernel + cylindrical source) |
| Numerical (∂lnX)² dla 4 field schedules | **CLOSED** (all heavy regime, default τ.3 params) |
| TT7-TT12 numerical re-derivation | **B7-v2 dedicated cycle pending** (geometric edge analysis) |
| **KEY PHYSICS: heavy regime → bulk signal ZERO** | **DOCUMENTED** (geometric pivot required) |

→ **B7 STRUCTURAL CLOSED**. Numerical TT7-TT12 update require dedicated
B7-v2 cycle z proper geometric edge analysis (boundary clock positioning
+ realistic detector volume integration).

## Files

- `B7_greens_function_derivation.py` — sympy two-regime derivation script (uruchomiony ✓)
- `B7_greens_function_results.md` — niniejszy plik (results documentation)

## Cross-references

- [[Phase1_results.md]] — T1.3 m_e_eff multiplicative (post-A5)
- [[Phase2_results.md]] — T2.4 Λ-scan dual-table (audit-aware)
- [[Phase3_results.md]] — TT7-TT12 numerical re-derivation B7-v2 pending
- [[../../meta/AUDYT_TGP_2026-05-01.md]] § B.7 + § L (A5 closure) + § O (B7 closure pending)
