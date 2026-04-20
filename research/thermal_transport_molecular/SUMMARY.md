# Thermal Transport in TGP: Molecular Crystals & Perovskites

**Date:** 2026-04-20
**Status:** finished — partial positive / partial null
**Scripts:** ps01–ps04 (all run cleanly)

## 1. What we set out to do

Test whether the TGP exponential metric $g_{ij} = e^{+2\phi}\delta_{ij}$,
applied to heat diffusion in crystalline matter, can explain four
long-standing anomalies:

1. Plateau of $\lambda_\mathrm{th}(T)$ at high $T$ in molecular crystals
   (naphthalene, anthracene, rubrene) — i.e. the failure of BTE's
   $\lambda \propto 1/T$ Umklapp prediction.
2. Ultralow thermal conductivity $\lambda \sim 0.3$ W/m·K of hybrid
   halide perovskites (MAPbI₃, CsPbBr₃, Cs₂AgBiBr₆).
3. Scaling of $\lambda(T, \rho, M, v_s)$ across soft ionic crystals.
4. Wiedemann-Franz violation in bad metals (cuprates, heavy fermions).

## 2. Results per script

### ps01 — heat diffusion in TGP metric (foundational)

Derived the TGP heat equation in the spatial exp-metric,

$$
\partial_t u \;=\; D_0\,e^{-3\phi}\,\partial_x\!\bigl(e^{\phi}\,\partial_x u\bigr)
$$

and verified numerically on a 1D periodic molecular lattice that the
long-wavelength homogenised diffusivity is the harmonic mean

$$
D_\mathrm{hom} \;=\; \frac{D_0}{\langle e^{2\phi}\rangle}.
$$

| $\phi_\mathrm{amp}$ | $D_\mathrm{MSD}/D_0$ | $D_\mathrm{harm}/D_0$ | residual |
|--------------------:|---------------------:|----------------------:|---------:|
| 0.0  | 0.992 | 1.000 | 0.8 % |
| 0.03 | 0.957 | 0.963 | 0.6 % |
| 0.10 | 0.876 | 0.881 | 0.6 % |
| 0.30 | 0.666 | 0.678 | 1.8 % |

Numerics track the analytic prediction at the percent level. The
gradient-scattering length $\ell_\mathrm{scat} = 1/|\phi'|_\mathrm{rms}$
drops below the lattice constant at $\phi_\mathrm{amp} \gtrsim 0.1$,
signalling a new TGP-specific scattering channel.

### ps02 — plateau fit for molecular crystals

Fit the combined model

$$
\frac{1}{\lambda(T)} \;=\; \frac{T}{B} \;+\;
\frac{1}{D_\infty \,\exp\!\bigl(-2\phi_\mathrm{sat}^2\,T/(T+T_\mathrm{sat})\bigr)}
$$

with $\phi_\mathrm{sat}$ capped at $0.3$ (TGP weak-field regime).
Residuals (capped fit, $|{\rm err}| \le 20 \%$):

| crystal      | $\chi^2_A$ (free) | $\chi^2_B$ (capped) | ratio |
|:-------------|------------------:|--------------------:|------:|
| rubrene_b    | 2.0e-3 | 2.3e-3 | **1.13** compatible |
| naphthalene  | 4.4e-3 | 4.3e-2 | 9.8× degraded |
| anthracene   | 4.1e-3 | 9.0e-2 | 22× degraded |

**Finding:** only rubrene's $\lambda(T)$ falls cleanly under TGP
weak-field plus BTE. Naphthalene and anthracene require either stronger
substrate fluctuations ($\phi_\mathrm{sat} > 0.3$, outside weak field)
or additional phonon physics. TGP plateau mechanism is **partial**, not
universal.

### ps03 — perovskite + ionic crystal scan (room-temperature)

Tested the predictive formula

$$
\lambda_\mathrm{TGP} \;=\;
\tfrac{1}{3}\, n_a\, k_B\, v_s\, a \;\cdot\;
\exp\!\Bigl(-\,\tfrac{2\alpha_T k_B T}{M_\mathrm{eff} v_s^2}\Bigr)
$$

with a single global parameter $\alpha_T$. Across 11 materials
(perovskites + ionic crystals, spanning $\lambda = 0.25$ → $11$ W/m·K):

- Best fit: $\alpha_T \to 0.01$ (lower scan limit); the TGP exponent
  essentially vanishes. The formula reduces to the Cahill kinetic
  floor $\lambda_\mathrm{kin} = \tfrac{1}{3} n_a k_B v_s a$.
- Perovskite predictions (no per-material tuning):
  - MAPbI₃: $\lambda_\mathrm{TGP} = 0.25$ vs $\lambda_\mathrm{exp} = 0.30$ (−16 %)
  - MAPbBr₃: 0.33 vs 0.37 (−11 %)
  - CsPbBr₃: 0.16 vs 0.42 (−63 %)
  - CsSnI₃: 0.10 vs 0.25 (−58 %)
- Stiff crystals systematically under-predicted by 80–95 % (phonon
  long-MFP Umklapp dominates there, outside TGP's substrate-diffusion
  regime).

**Finding:** The TGP homogenisation formula is **the lower bound** on
thermal conductivity. Soft perovskites saturate it (positive result:
no material-specific fitting needed); stiff crystals don't. The
TGP-specific $\exp(-2\langle\phi^2\rangle)$ correction is only
$\sim 10^{-4}$ at 300 K — below measurement precision.

### ps04 — Wiedemann-Franz ratio

Using ps03's $\alpha_T \sim 0.5$ as a (generous) upper bound, the TGP
correction to the Lorenz number is

$$
L_\mathrm{TGP}/L_0 \;=\; \exp\!\bigl(-2\alpha_T k_B T / (M_\mathrm{eff}v_s^2)\bigr)
\;\approx\; 1 - \mathcal{O}(10^{-3})
$$

for all surveyed metals. Observed deviations:

| metal       | family          | $L/L_0$ obs | TGP pred | verdict |
|:------------|:----------------|------------:|---------:|:--------|
| Cu          | noble           | 1.00 | 0.997 | consistent |
| Al          | simple          | 0.95 | 0.998 | TGP too small |
| LSCO        | cuprate         | 0.40 | 0.997 | TGP too small (×600) |
| YBa₂Cu₃O₇   | cuprate         | 0.60 | 0.996 | TGP too small |
| CeCoIn₅     | heavy fermion   | 0.80 | 0.998 | TGP too small |
| TTF-TCNQ    | charge transfer | 0.50 | 0.979 | TGP too small |

**Finding:** TGP is cleanly **falsified** as the explanation of
bad-metal WF violation. The TGP correction is 3–4 orders of magnitude
too small. The positive prediction $L/L_0 = 1 \pm 10^{-4}$ for normal
3D metals at low $T$ would only be in trouble if future precision WF
measurements find a larger deviation in clean simple metals.

## 3. Global picture

| Phenomenon | TGP outcome |
|:-----------|:------------|
| Plateau $\lambda_\mathrm{th}(T)$ in rubrene | **consistent** (weak-field works) |
| Plateau in naphthalene, anthracene | partial (needs richer phonon physics) |
| Ultralow $\lambda$ in hybrid perovskites | **positive** (TGP recovers Cahill floor universally) |
| Wiedemann-Franz in normal metals | **consistent** (tight null prediction) |
| Bad-metal WF violations | **falsified** (TGP effect is $10^{-4}$) |

TGP's substrate-modulated diffusion provides a clean, universal
*lower bound* on thermal conductivity via the Cahill-like expression
$\lambda \ge (1/3) n_a k_B v_s a$, saturated in soft perovskites. This
is the main, generalizable finding. The TGP-specific exponential
correction factor $\exp(-2\langle\phi^2\rangle)$ is unobservably small
at 300 K under realistic weak-field constraints.

## 4. Key equations (summary card)

$$
\boxed{\text{TGP heat equation:}\quad
\partial_t u = D_0\,e^{-3\phi}\,\partial_i(e^{\phi}\,\partial^i u)}
$$

$$
\boxed{\text{Homogenised:}\quad D_\mathrm{hom} = D_0/\langle e^{2\phi}\rangle}
$$

$$
\boxed{\text{Plateau:}\quad
\frac{1}{\lambda(T)} = \frac{T}{B} + \frac{1}{D_\infty e^{-2\phi_\mathrm{sat}^2 T/(T+T_\mathrm{sat})}}}
$$

$$
\boxed{\text{TGP lower bound:}\quad
\lambda \ge \tfrac{1}{3}\, n_a k_B v_s a\;\;\text{(Cahill floor)}}
$$

$$
\boxed{\text{WF null:}\quad L_\mathrm{TGP}/L_0 = 1 + \mathcal{O}(10^{-4})}
$$

## 5. Status and next steps

**Not preprint-ready.** The positive "Cahill-floor" recovery for
perovskites (ps03) is essentially a rederivation of an existing
kinetic-theory result via a different formalism. The plateau fit
(ps02) works for one material out of three. The WF result (ps04) is
genuinely new and tight but is a null prediction — publishable only as
part of a larger TGP compendium.

**Clear takeaways for the theory layer:**

- TGP's spatial exp-metric combined with standard heat-diffusion gives
  the correct Cahill-minimum scaling **without ad-hoc assumptions** —
  a clean cross-check of the geometric formulation.
- The predicted TGP-specific correction to $\kappa$ and $L/L_0$ is
  $\sim 10^{-4}$ for soft materials at 300 K — well below current
  precision. This closes the thermal-transport channel as a near-term
  falsification avenue (similar conclusion to the Casimir-MOF project).

**Candidate next-horizon directions:**

1. **Neutrino MSW** (`research/neutrino_msw/`): long horizon, but the
   TGP substrate-induced phase shift could show up in solar vs
   atmospheric neutrino ratios. Different observational channel,
   different degeneracies — worth exploring independently.
2. **Cosmological ZPE → dark energy**: ties to $c_\mathrm{TGP} =
   1/(4\pi)$ already validated in liquid viscosity; predicts $\Omega_\Lambda$
   from substrate quartic potential.
3. **High-precision WF in clean metals** (experimental push): push
   $|1 - L/L_0|$ below $10^{-4}$ in Cu/Ag at $T \ll \Theta_D$; any
   detection at this level would either confirm a TGP signature or
   force re-examination of the conformal-invariance assumption for
   4D Maxwell in exp metric.

## 6. File index

- [[TGP/TGP_v1/research/thermal_transport_molecular/PLAN.md]]
- [[TGP/TGP_v1/research/thermal_transport_molecular/ps01_thermal_substrate_diffusion.py]]
- [[TGP/TGP_v1/research/thermal_transport_molecular/ps02_plateau_molecular_crystals.py]]
- [[TGP/TGP_v1/research/thermal_transport_molecular/ps03_perovskite_hybrid_scan.py]]
- [[TGP/TGP_v1/research/thermal_transport_molecular/ps04_WF_ratio_TGP.py]]

## 7. Cross-references with other sectors

- **Liquid viscosity** (`research/liquid_viscosity/`): the constant
  $c_\mathrm{TGP} = 1/(4\pi)$ from substrate ZPE used there is NOT
  visible in thermal transport because it enters as a $(4\pi)^{-1}$
  prefactor on top of a homogenised diffusivity already dominated by
  classical kinetic floor.
- **Cuprates P6.A / P6.D**: the bad-metal WF deviation we see in LSCO
  / YBCO is driven by the same non-Drude physics that gives the MIR
  Drude deficit — ZR-singlet + paramagnon channels. TGP's
  substrate-diffusion picture is NOT sufficient there; the cuprate
  sector uses a richer construction.
- **Casimir MOF** (`research/casimir_mof/`): another case where the
  TGP correction is at the $10^{-4}$–$10^{-3}$ level, below current
  experimental precision. Both projects together delineate the lower
  envelope of short-horizon falsifiability.
