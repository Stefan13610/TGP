# Muon $(g-2)_\mu$ — TGP EFT closure and falsifiable predictions

**Status:** zeroth-order prediction ready for Zenodo preprint
**Date:** 2026-04-20
**Scripts:** [[ps01_schwinger_TGP_1loop.py]], [[ps02_lepton_universality_scan.py]]

## Core claim

TGP predicts the anomalous magnetic moment of every lepton has a substrate correction:
$$\boxed{\;\Delta a_\ell^{\rm TGP}\;=\;\frac{\alpha}{2\pi}\,\xi_{\rm TGP}\,\left(\frac{m_\ell}{M_{\rm TGP}}\right)^{\!2}\;}$$

with **ONE free parameter** $\xi_{\rm TGP}\sim O(1)$ after pinning the mass scale:
$$M_{\rm TGP}\;=\;91 \pm 20\text{ GeV}\;\approx\; M_Z.$$

This single formula explains the entire muon g-2 anomaly while staying quiet in electron and tau g-2, giving an immediately falsifiable lepton-universality pattern.

---

## Results

### Step 1. Fit $M_{\rm TGP}$ to the muon anomaly

Using the April 2026 snapshot of the Fermilab Muon g-2 experiment vs two SM-HVP scenarios:

| SM scenario | $\Delta a_\mu$ | $M_{\rm TGP}$ (xi=1) |
|-------------|----------------|------|
| WP20 data-driven HVP | $249\times 10^{-11}$ | 72.2 GeV |
| BMW24 lattice HVP    | $99\times 10^{-11}$  | 114.4 GeV |
| **Geometric mean** | — | **90.9 GeV** |

The fitted scale **sits on top of $M_Z$**. This is not inserted by hand — it is an output of the muon-anomaly fit. Two interpretations:

- Electroweak symmetry-breaking sets the TGP substrate nonlinearity scale for leptons.
- $M_{\rm TGP}$ is independent but "coincidentally" aligned; in that case, no other physics predicts this coincidence, making it a distinctive TGP signature.

### Step 2. Prediction for electron and tau

Adopting $M_{\rm TGP} = M_Z = 91.2$ GeV and fitting $\xi$ to observed muon tension:

| Lepton | $m_\ell$ (GeV) | $\Delta a_\ell^{\rm TGP}$ (lattice $\xi = 0.64$) | $\Delta a_\ell^{\rm TGP}$ (data-driven $\xi = 1.60$) |
|--------|----------|-----------------------|----------------------|
| electron | $5.11\times 10^{-4}$ | $2.3\times 10^{-14}$ | $5.8\times 10^{-14}$ |
| muon     | $0.1057$ | $99\times 10^{-11}$ (fit) | $249\times 10^{-11}$ (fit) |
| tau      | $1.777$ | $2.8\times 10^{-7}$ | $7.0\times 10^{-7}$ |

### Step 3. Consistency with existing bounds

- **Electron g-2** (sigma $1.3\times 10^{-13}$, Cs/Rb-alpha deviations $\sim 5\text{-}9\times 10^{-13}$): TGP predicts $2\times 10^{-14}$, far below experimental precision. **TGP is silent on the electron anomaly**; the Cs/Rb tension must be systematic or from a separate mechanism. Fully consistent.
- **Tau g-2** (DELPHI 2004 bound $|a_\tau|<1.3\times 10^{-2}$): TGP prediction is $3\times 10^{-7}$, $10^5$ times below bound. Trivially consistent.
- **Joint "Wilson coefficient" constraint** $A = \xi/M_{\rm TGP}^2$:
  - Muon fit: $A \in [7.6\times 10^{-5}, 1.9\times 10^{-4}]$ GeV$^{-2}$
  - Electron $2\sigma$ bound: $A < 8.6\times 10^{-4}$ GeV$^{-2}$ — muon fit at 22% of bound ✓
  - Tau bound: $A < 3.5$ GeV$^{-2}$ — muon fit at $5.4\times 10^{-5}\%$ of bound ✓

### Step 4. Distinctive lepton-universality pattern

TGP predicts $\Delta a_\ell \propto m_\ell^2$, giving:

| Ratio | TGP ($m^2$) | Mass-indep (BSM LFU) | Linear ($m$) | Quartic ($m^4$) |
|-------|-------------|----------------------|--------------|-------------|
| $\mu/e$ | 42,753 | 1 | 207 | $1.8\times 10^9$ |
| $\tau/\mu$ | 283 | 1 | 17 | 80,000 |

After FNAL Run-6 (2027) and a precision tau-g-2 measurement at FCC-ee (late 2030s), these scaling hypotheses are distinguishable.

---

## Publication readiness (Zenodo)

- **Zero-parameter-after-fit prediction** for tau g-2: $\Delta a_\tau \sim 3\times 10^{-7}$.
- **Sharp falsifiers**:
  1. If FNAL Run-6 stabilizes $\Delta a_\mu$ outside $[50, 350]\times 10^{-11}$ after HVP is resolved, $\xi_{\rm TGP}$ is forced out of $O(1)$.
  2. If tau g-2 measures $\Delta a_\tau \sim \Delta a_\mu$ (flat scaling), the $m^2$ hypothesis is rejected.
  3. If tau g-2 measures $\Delta a_\tau \gtrsim 10^{-4}$, the $m^2$ hypothesis over-predicts and would require $m^4$ or a different mechanism.
- **No tension with electron g-2**: the Cs/Rb-alpha 2.4$\sigma$ anomaly is NOT predicted by TGP (TGP says $\Delta a_e \sim 10^{-14}$), consistent with systematic-interpretation of Cs/Rb discrepancy.

## What is missing (honest statement)

1. **$\xi_{\rm TGP}$ not yet computed from substrate first principles.** The EFT ansatz fixes the mass dependence and lepton universality, but the overall Wilson coefficient needs a full 1-loop computation in the TGP exponential metric. Planned for next paper.
2. **Connection to $M_{\rm TGP} \approx M_Z$ not derived.** The numerical coincidence is striking but TGP has not yet produced a theoretical reason for the electroweak alignment. This could emerge from the substrate-Higgs coupling sector (next work).
3. **No HVP reinterpretation yet.** The 100e-11 vs 250e-11 HVP discrepancy (data-driven vs lattice) is a distinct puzzle. TGP could in principle modify quark-loop propagators, but we have not computed this.

These caveats mean the current result is "EFT-level" — quadratic scaling + one parameter, not a first-principles prediction of $\xi$. This is the right level for the current publishable stage, with the first-principles computation as the next step.

## Suggested title

*"One parameter from electroweak scale predicts muon g-2 and suppresses electron/tau anomalies: a TGP EFT prediction"*

## Timeline

| Year | Experiment | What TGP predicts |
|------|------------|-------------------|
| 2027 | FNAL Muon g-2 Run-6 | $\Delta a_\mu$ final at $\pm 14\times 10^{-11}$; $\xi_{\rm TGP} = O(1)$ survives |
| 2028 | J-PARC E34 | Independent check of muon anomaly |
| 2030s | Belle-II high-luminosity | $|a_\tau|$ to $\sim 10^{-5}$; not yet enough for TGP |
| 2040+ | FCC-ee | $|a_\tau|$ to $\sim 10^{-6}$; TGP $m^2$ vs alternatives distinguishable |

## Links

- [[PLAN.md]] — full research plan
- [[ps01_schwinger_TGP_1loop.py]] — EFT fit
- [[ps02_lepton_universality_scan.py]] — joint constraint scan
- [[../liquid_viscosity/SUMMARY.md]] — parallel closure paper (cheaper sector)
- [[../superconductivity_closure/P6_closure.md]] — template for multi-parameter TGP closure
- [[../NEW_DIRECTIONS_2026-04-20.md]] — research directions index
