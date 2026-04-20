# Casimir effect in MOF geometries — TGP analysis

**Status:** research complete; clean null prediction at current XRD precision.
**Date:** 2026-04-20
**Scripts:** [[ps01_casimir_parallel_plates.py]], [[ps02_casimir_cylindrical_pore.py]], [[ps03_casimir_MOF_cage.py]], [[ps04_osmotic_vs_casimir.py]], [[ps05_gradient_phi_correction.py]]

## Core finding

TGP predicts the Casimir effect in MOF cavities should DEVIATE from the classical
QED result by only the QED trace-anomaly correction:

$$
P_{\text{Cas}}^{\text{TGP}} \;=\; P_{\text{Cas}}^{\text{cl}} \; \Big(1 + \tfrac{\alpha}{\pi} \cdot \tfrac{a_\text{Bohr}}{R}\Big)
\;\approx\; P_{\text{Cas}}^{\text{cl}} \;(1 + 3 \times 10^{-4})
\quad \text{for MOF-5 cage } (R = 0.4\ \text{nm}).
$$

**Why so small:** the TGP exp-metric coupling $g_{\mu\nu} = \eta_{\mu\nu}\,e^{2\Phi/\Phi_0}$
is CLASSICALLY CONFORMAL for 4D Maxwell. Tree-level photon modes don't feel $\Phi$;
only the QED trace anomaly (1-loop) generates a small $\Phi$–$F^2$ coupling.

This gives a clean **null prediction**:
- $\Delta P / P \sim 3 \times 10^{-4}$ in MOF-5 small cage ($R = 0.395$ nm).
- Induced MOF-5 lattice strain: $\epsilon_{\text{TGP}} \sim 1.07$ ppm.
- Current synchrotron XRD precision: 1–2 ppm.
- → Borderline detectable at best-case precision; out of reach for routine XRD.
- Any measurement of a *larger* Casimir anomaly would **falsify** TGP's exp-metric form.

---

## Five-step numerical programme

### ps01 — Parallel plates (method verification)

Reproduced classical Casimir $E/A = -\pi^2\hbar c/(1440\,d^3)$ for Dirichlet scalar via three independent methods:
1. Analytic closed form.
2. Abel-Plana integral: $\int_0^\infty t^3/(e^{2\pi t}-1)\,dt = 1/240$ (verified to $2\times 10^{-16}$).
3. Zeta regularization: $\zeta(-3) = 1/120$ (verified to $6\times 10^{-16}$).

All three agree to machine precision. Method validated for more complex geometries.

**SI predictions** (classical EM, perfect conductors):
| $d$ (nm) | $P$ (Pa) | $P$ (GPa) |
|----------|----------|-----------|
| 0.5 | −2.08 × 10¹⁰ | −20.8 |
| 1.0 | −1.30 × 10⁹ | −1.30 |
| 1.1 (MOF-5) | −8.88 × 10⁸ | −0.888 |
| 1.5 (MOF-177) | −2.57 × 10⁸ | −0.257 |

### ps02 — Cylindrical pore (SBA-15, MCM-41)

Bessel-zero mode structure for Dirichlet scalar in a cylinder of radius $R$.
- Tabulated $\alpha_{m,n}$ for $m = 0..50$, $n = 1..80$.
- Weyl's law $K(\tau) \sim A/(4\pi\tau) - P/(8\sqrt{\pi\tau}) + \chi/6$ verified at small $\tau$.
- Heat-kernel / proper-time regularization implemented.
- Full numerical extraction of the Casimir coefficient is partial (higher Seeley terms needed); we use literature DeRaad-Milton $C_{\text{EM}} = -0.01356$ for SI predictions.

**Cylindrical pore Casimir pressures**:
| Material | $R$ (nm) | $P$ (Pa) | $P$ (MPa) |
|----------|----------|----------|-----------|
| SWNT (5,5) | 0.34 | −1.02 × 10¹⁰ | −10,200 |
| MOF channel | 0.50 | −2.18 × 10⁹ | −2,180 |
| MCM-41 | 1.25 | −5.59 × 10⁷ | −55.9 |
| SBA-15 | 3.50 | −9.09 × 10⁵ | −0.91 |

### ps03 — Spherical cavity (MOF cages)

Spherical-Bessel zero tabulation for the 3-ball Dirichlet Laplacian.
- Verified $j_{0,n} = n\pi$ to machine precision.
- $j_{1,1} = 4.4934$ (the famous "tan(x) = x" root), all correct.
- Weyl 3-ball law $K(\tau) \sim V/(4\pi\tau)^{3/2} - S/(16\pi\tau) + \dots$ verified.

**Used literature Boyer coefficient** $C_{\text{EM}} = +0.04618$ (repulsive, pushes walls OUT).

| MOF cage | $R$ (nm) | $P$ (Pa, naive Boyer) | $P$ (GPa) |
|----------|----------|------------------------|-----------|
| MOF-5 small | 0.395 | +4.77 × 10⁹ | +4.77 |
| MOF-5 large | 0.590 | +9.59 × 10⁸ | +0.96 |
| UiO-66 oct | 0.555 | +1.22 × 10⁹ | +1.22 |
| ZIF-8 | 0.575 | +1.06 × 10⁹ | +1.06 |
| MOF-210 | 1.120 | +7.38 × 10⁷ | +0.074 |

### ps04 — Osmotic vs Casimir separation (Ar@MOF-5, 87 K)

Critical finding: **the naive Boyer Casimir (4.8 GPa in MOF-5) is ~30% of
MOF-5 bulk modulus (15 GPa) — inconsistent with MOF stability**.

Therefore, the naive perfect-conductor Boyer coefficient overestimates by
$\sim 100\times$ in real MOFs. Suppressions:
- **Dielectric walls** (organic linkers, $\varepsilon_r \sim 2-4$): Lifshitz factor ~0.3-0.5.
- **Open-pore connectivity**: mode leakage → effective coefficient reduced ~5-10×.

**Realistic effective coefficient:** $C_{\text{eff}} \sim 5 \times 10^{-4}$ (100× Boyer suppression),
giving $P_{\text{Cas,MOF-5}} \sim 50$–$100$ MPa.

**Separation strategy:**
1. Measure lattice parameter $a(n)$ at varying Ar loadings $n$.
2. Osmotic contribution $P_{\text{osm}}(n) = n k_B T \cdot (1 + B_2 n)$ varies with loading;
   Casimir doesn't.
3. Extrapolate $a(n \to 0)$ to isolate pure Casimir strain.
4. Residual strain $\epsilon_{\text{Cas}} \sim 3-6$ ppm (realistic scenario).

**Alternative:** low-temperature thermal-expansion anomaly.
MOF-5 NTE ($\alpha = -13$ ppm/K, Zhou 2008) goes to zero at $T \to 0$; Casimir strain is T-independent.
Below ~20 K, lattice-parameter plateau should reveal Casimir.

### ps05 — TGP gradient-substrate correction

**Physics:** TGP exp-metric coupling is conformal at tree level for Maxwell in 4D:
$$
S_{\text{EM}} = -\tfrac{1}{4}\int d^4x\,\sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} = -\tfrac{1}{4}\int d^4x\,F^2
$$
No tree-level $\Phi$–photon coupling. The *only* non-zero TGP correction at leading order
comes from the **QED trace anomaly**:
$$
\mathcal{L}_{\text{anom}} = \tfrac{\beta(\alpha)}{4\alpha} \tfrac{\delta\Phi}{\Phi_0} F^2,
\qquad \tfrac{\beta}{4\alpha} = \tfrac{\alpha}{2\pi} \approx 1.16 \times 10^{-3}.
$$

**Prediction:** $P_{\text{Cas}}^{\text{TGP}} = P_{\text{Cas}}^{\text{cl}} \, (1 + (\alpha/\pi)(a_{\text{Bohr}}/R))$.

| MOF cage | $R$ (nm) | $\delta\Phi/\Phi_0$ | correction | strain (ppm) |
|----------|----------|---------------------|------------|--------------|
| MOF-5 small | 0.395 | 0.134 | +3.1 × 10⁻⁴ | +1.07 |
| MOF-5 large | 0.590 | 0.090 | +2.1 × 10⁻⁴ | +0.14 |
| MOF-177 | 0.735 | 0.072 | +1.7 × 10⁻⁴ | +0.05 |
| UiO-66 oct | 0.555 | 0.095 | +2.2 × 10⁻⁴ | +0.20 |
| ZIF-8 | 0.575 | 0.092 | +2.1 × 10⁻⁴ | +0.16 |

**Verdict:** strain ~1.07 ppm in MOF-5 small cage is at current XRD precision threshold (1-2 ppm).
TGP *marginally* detectable; most likely below noise in realistic setups.

---

## Publication readiness

**Status:** NOT preprint-ready. This is a NULL prediction, not a confirmed phenomenon.

**Summary statement for future reference:**
- TGP is *consistent* with classical Casimir in MOF geometries by virtue of 4D Maxwell conformal invariance.
- The only TGP-specific correction is $\sim \alpha/\pi \times (a_{\text{Bohr}}/R)$, giving ~1 ppm strain in MOF-5.
- Current XRD at 1-2 ppm precision cannot distinguish TGP from classical.
- FCC-era synchrotron (0.1 ppm) would enable a 10σ test.

**Sharp falsifier:** If any MOF Casimir measurement shows $\Delta P/P > 10^{-3}$ (0.1% level)
deviation from classical Lifshitz prediction (beyond the trace-anomaly bound), TGP's
exp-metric hypothesis is falsified.

## Suggested title for future work

*"Casimir effect in MOF geometries: null prediction from TGP's conformally-invariant exp-metric"*

## Identified open questions

1. **Boundary conditions on $\Phi$ at MOF walls.** We assumed Dirichlet for the
   substrate field at the material interface. Proper derivation requires extending
   TGP to include wall-substrate coupling. (Relates to Axiom A-IV: "space generated by matter".)

2. **$\delta\Phi/\Phi_0$ atomic profile.** Our dimensional estimate $a_\text{Bohr}/R$ is crude.
   A first-principles calculation from the TGP core action near a metal atom would pin down
   $\delta\Phi/\Phi_0$ more accurately. This is tied to the same substrate-matter coupling
   programme as the SC/muon physics.

3. **Casimir-Polder for guest molecules.** A single guest molecule inside a MOF cage
   sees a Casimir-Polder potential. The TGP correction to C-P in confined geometry
   hasn't been computed here — potentially another null check but might amplify if
   Phi gradients are larger near C-P test charges.

4. **Repulsive Casimir confirmation.** The Boyer repulsive sign is essential for our
   analysis. Experimental confirmation of Boyer's result in spherical microcavities
   (done only approximately so far, Chan et al. 2008 for a sphere-plate geometry).

## Relationship to other TGP programmes

| Programme | Connection |
|-----------|------------|
| [[../liquid_viscosity/SUMMARY.md]] | Uses substrate ZPE constant $c_{\text{TGP}} = 1/(4\pi)$; Casimir and viscosity both probe substrate zero-point physics. |
| [[../superconductivity_closure/P6_closure.md]] | Electron-phonon coupling via exp-metric provides analog to $\beta$-coefficient here. |
| [[../muon_g_minus_2/SUMMARY.md]] | Also uses QED trace anomaly (ps05 there); but for Casimir the anomaly is suppressed, whereas for $a_\mu$ it competes with direct Yukawa. |

## Next research direction

With Casimir-MOF now a documented null prediction, best-return TGP research directions:

1. **Neutrino MSW oscillations** (from NEW_DIRECTIONS_2026-04-20, status: pending).
   Substrate-induced matter-like phase shift for neutrinos might give observable
   deviation in solar/atmospheric ratio. More promising for positive TGP signal.

2. **Thermal transport in molecular solids** (NEW_DIRECTIONS, pending).
   Substrate-coupled phonon scattering from the TGP exp-metric has specific predictions
   for thermal conductivity minima. Connects to the $c_{\text{TGP}} = 1/(4\pi)$ constant.

3. **Cosmological substrate ZPE → dark energy.**
   If $c_{\text{TGP}} \times \hbar \omega_\Phi^4$ gives a vacuum energy density compatible
   with $\Lambda_{\text{obs}} \sim 10^{-122} M_{\text{Pl}}^4$, this is the biggest possible
   TGP closure: unify liquid viscosity, SC topology, Casimir bound, and cosmological constant
   under ONE substrate scale.

## Links

- [[ps01_casimir_parallel_plates.py]] — classical Casimir benchmark (method verified).
- [[ps02_casimir_cylindrical_pore.py]] — cylindrical pore mode structure.
- [[ps03_casimir_MOF_cage.py]] — spherical cavity (MOF cage).
- [[ps04_osmotic_vs_casimir.py]] — osmotic subtraction strategy.
- [[ps05_gradient_phi_correction.py]] — TGP trace-anomaly correction.
- [[PLAN.md]] — original research plan (2026-04-20).
- [[../NEW_DIRECTIONS_2026-04-20.md]] — index of active TGP programmes.
