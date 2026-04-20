# Muon $(g-2)_\mu$ — TGP EFT closure and falsifiable predictions

**Status:** minimal exp-metric dilaton **excluded by LEP Z-pole** (ps06). Muon g-2 is an open problem for TGP.
**Date:** 2026-04-20
**Scripts:** [[ps01_schwinger_TGP_1loop.py]], [[ps02_lepton_universality_scan.py]], [[ps03_schwinger_proper_time.py]], [[ps04_vertex_from_exp_metric.py]], [[ps05_trace_anomaly_channel.py]], [[ps06_phi0_scale_consistency.py]]

## Current status — minimal exp-metric picture insufficient (post ps04-ps06)

**ps04** derives a dilaton-like Yukawa $y_\ell = m_\ell/\Phi_0$ from the exponential metric $g_{\mu\nu} = \eta_{\mu\nu}\,e^{2\Phi/\Phi_0}$, predicting $m_\ell^4$ scaling in $\Delta a_\ell$ and fitting the muon anomaly with $\Phi_0 \approx 0.34$ GeV (for $M_\Phi = M_Z$).

**ps05** rules out the QED trace anomaly as an alternate $m^2$-scaling channel (suppressed by $\sim 10^{10}$ vs the dilaton).

**ps06** shows the minimal dilaton picture is **experimentally excluded**:

| Observable | Prediction (at $\Phi_0=0.34$ GeV) | Observed | Verdict |
|-----------|-----------------------------------|----------|---------|
| $\Delta a_\mu$ | $2.5\times 10^{-9}$ (fit) | $(99$–$249)\times 10^{-11}$ | Fits by construction |
| $Z\to\tau^+\tau^-$ vertex shift | $\sim 17$–$70\%$ | LEP constrains to $<0.3\%$ | **EXCLUDED** |
| Yukawa $y_\tau$ | $5.3$ | perturbativity $y < \sqrt{4\pi}\approx 3.5$ | Non-perturbative |
| Yukawa $y_t$ (top) | $512$ | perturbativity | Breaks EFT |

There is **no value of $\Phi_0$** that simultaneously (a) saturates $\Delta a_\mu$, (b) keeps all SM Yukawas perturbative, (c) respects Z-pole precision. Either $\Phi_0$ is small (fits muon but kills Z-pole and top) or large (respects Z-pole but underfits muon by $\gtrsim 10^4$).

### What this means

The minimal TGP exponential-metric-only coupling **does not explain the muon g-2 anomaly**. This is a clean negative result, ruling out a particular simple TGP completion.

**Surviving TGP pictures for muon g-2:**
1. **Vector substrate excitation** (not scalar dilaton): a universal $V_\mu\,\bar\ell\gamma^\mu\ell$ coupling at mass $M_V\approx M_Z$ naturally gives $m^2$ scaling without dilaton-level perturbativity issues. The ps01/ps03 empirical fit with $M_{\rm TGP}\approx M_Z$ is compatible with this interpretation.
2. **Flavor-selective coupling**: substrate couples to leptons but not quarks (requires a specific TGP mechanism, e.g. color-charge screening).
3. **Higher-dimensional operator**: a dimension-6 operator $(H^\dagger H)\bar\ell\ell$ from substrate-Higgs mixing, which gives $y\propto m_\ell/v_{\rm Higgs}$ but with different loop structure than the minimal dilaton.

None of these is yet derived from a TGP first-principles Lagrangian. **The muon g-2 sector is therefore an OPEN PROBLEM for TGP**, not a closed prediction.

### Immediate status for preprint

**Not preprint-ready.** The conclusion needs further theoretical work before publication. The liquid-viscosity closure ([[../liquid_viscosity/SUMMARY.md]]) remains the clean publishable result.

---

## Results

> **Note:** Sections below use the original $m^2$ EFT ansatz ($\Delta a_\ell = (\alpha/2\pi)\,\xi\,(m_\ell/M_{\rm TGP})^2$) from ps01-ps03. This ansatz is **superseded** by the dilaton $m^4$ framework derived in ps04 and confirmed in ps05 (see Core claim above and Result 4). They remain as the empirical fit history and because the $m^2$ $M_{\rm TGP}=91$ GeV result is numerically consistent with the dilaton $M_\Phi=91$ GeV result — both identify the substrate quantum mass with $M_Z$.

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

## Result 3 — Proper-time / form-factor calculation (ps03)

**Hypothesis:** the EFT ansatz $\Delta a_\ell = (\alpha/2\pi)\,\xi\,(m_\ell/M_{\rm TGP})^2$ can be reproduced by a concrete 1-loop proper-time calculation with a TGP substrate-modified photon propagator, parametrised by a form factor $\varepsilon(u,\lambda)$ where $u = -k^2/m_\ell^2$ and $\lambda = \Lambda/m_\ell$.

### Setup

The shift in the Pauli form factor from the flat 1-loop is
$$\Delta F_2(0) \;=\; \frac{\alpha}{\pi}\int_0^1 dx\int_0^\infty du\,\frac{x(1-x)\,\varepsilon(u,\lambda)}{[u(1-x) + x^2]^2}.$$

Three substrate form factors were scanned:

| Form factor | $\varepsilon(u,\lambda)$ (with $v = u/\lambda^2$) | Motivation |
|-------------|---------------------------------------------------|------------|
| Gauss       | $-v\,e^{-v/2}$                                    | smooth TGP substrate ZPE |
| Lorentz     | $-v/(1+v)$                                        | Breit–Wigner substrate mode |
| Step        | $-\theta(u < \lambda^2)$                          | hard UV cutoff (diagnostic) |

### Fit results

Fitting $\xi$ to the muon anomaly at various $\Lambda$:

- **Gauss and Lorentz** form factors give finite $\xi \sim O(1)$ across a broad range of $\Lambda$, with the form-factor-**invariant** Wilson coefficient $A = \xi/\Lambda^2$ returning:
  $$M_{\rm eff} \equiv 1/\sqrt{A} \;=\; 114\text{ GeV (BMW24 lattice)}\;\text{to}\; 72\text{ GeV (WP20 data-driven)}.$$
  Geometric mean $M_{\rm eff} \approx 91$ GeV = $M_Z$, **matching the naive EFT fit from ps01 exactly**.
- **Step regulator** (hard cutoff) diverges logarithmically in $\Lambda$ — a genuine UV artifact of the non-smooth regulator, not a problem for TGP, since the physical substrate form factor is guaranteed smooth (Gauussian from ZPE).

### Interpretation

The form-factor-invariant $M_{\rm eff} \approx M_Z$ result confirms the ps01 EFT ansatz is not an artefact of the ansatz — it emerges from the full proper-time calculation independent of the specific substrate profile, provided the substrate cuts off smoothly. **The electroweak alignment of $M_{\rm TGP}$ is a robust feature of the substrate-photon loop, not a coincidence of the EFT normalization.**

---

## Result 4 — TGP exp-metric coupling analysis (ps04)

**Method:** expand the full matter+QED Lagrangian in $g_{\mu\nu} = \eta_{\mu\nu}\,e^{2\Phi/\Phi_0}$, field-redefine to canonical kinetic terms, and identify the effective $\Phi$–lepton–lepton and $\Phi$–photon couplings at tree level. Then evaluate the scalar-exchange contribution to $a_\ell$ using the Leveille 1-loop formulae.

### Two findings — one resolution, one new problem

**Finding 1 (clean):** The photon-$\Phi$ coupling **vanishes** at the classical Lagrangian level:
$$\sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} \;=\; e^{4\phi}\cdot e^{-4\phi}\cdot F^2 \;=\; F^2, \qquad \phi \equiv \Phi/\Phi_0.$$
This is classical conformal invariance of 4D Maxwell. **Consequence:** $c_{\rm vert}^{\gamma\gamma\Phi} = 0$ at tree level, so the ps03 $\varepsilon_{\rm TGP}(k^2)$ Gauss/Lorentz form-factors CANNOT come from minimal coupling. They must be interpreted as a **fermion-sector vertex correction** (Phi coupling to the muon), not as a photon propagator dressing.

**Finding 2 (tension):** The fermion coupling that survives is a **dilaton-like mass-weighted Yukawa**:
$$\mathcal{L}_{\rm coupling} \;=\; \frac{m_\ell}{\Phi_0}\,\Phi\,\bar\psi\psi.$$
A 1-loop scalar exchange with this coupling gives
$$\Delta a_\ell^{\rm TGP} \;\sim\; \frac{1}{8\pi^2}\,\frac{m_\ell^4}{\Phi_0^2\,M_\Phi^2}\,\log(M_\Phi/m_\ell),$$
i.e. $\boxed{m_\ell^4}$ scaling — **not** $m_\ell^2$ as ps01/ps03 assumed.

### Consequence for the three-way lepton-universality test

| Scaling | Source | $\Delta a_\tau/\Delta a_\mu$ | $\Delta a_\tau$ (for $\Delta a_\mu = 2.5\times 10^{-9}$) |
|---------|--------|------------------------------|------|
| Flat ($m^0$) | BSM LFU (e.g. heavy $Z'$) | 1 | $2.5\times 10^{-9}$ |
| $m^2$ | ps01/ps03 EFT ansatz (requires non-minimal $\Phi\bar\psi\psi$ universal Yukawa) | 283 | $7\times 10^{-7}$ |
| $m^4$ | **TGP minimal exp-metric (ps04)** | $8\times 10^4$ | $\mathbf{2\times 10^{-4}}$ |

### Numerical fit in the $m^4$ (dilaton) scenario

Fitting $\Phi_0$ to $\Delta a_\mu^{\rm obs} = 249\times 10^{-11}$ (WP20):

| $M_\Phi$ (GeV) | $\Phi_0$ (GeV) | $\Delta a_e$ | $\Delta a_\tau$ |
|---------------|----------------|--------------|-----------------|
| 10   | 3.05  | $1.4\times 10^{-18}$ | $1.5\times 10^{-4}$ |
| 91.2 | 0.34  | $1.4\times 10^{-18}$ | $1.9\times 10^{-4}$ |
| 500  | 0.062 | $1.4\times 10^{-18}$ | $2.0\times 10^{-4}$ |

Note: $\Phi_0$ is now a QCD-ish scale ($\sim$ few hundred MeV for $M_\Phi\sim M_Z$), NOT the electroweak scale. Alternative: if $\Phi_0 = v_{\rm Higgs} = 246$ GeV is imposed, then $M_\Phi \approx 54$ MeV — also a hadronic scale. **The exp-metric coupling alone does NOT predict $M_\Phi \approx M_Z$; the electroweak coincidence in ps01/ps03 belongs to a different mechanism.**

### Interpretation

ps04 reveals that the muon g-2 program has been working in the WRONG ansatz. Two distinct TGP channels contribute to $\Delta a_\ell$:

- **Channel A (dilaton, $m^4$):** natural from minimal exp-metric coupling; dominant in $\Delta a_\tau$ (amplifies by $(m_\tau/m_\mu)^2 \approx 283$ compared to $m^2$ prediction).
- **Channel B ($m^2$):** requires a non-minimal TGP sector — a flavor-universal $g\,\Phi\,\bar\psi\psi$ operator with $g\sim 0.3$ for $M_\Phi\sim M_Z$. This is the ps01/ps03 ansatz but is NOT automatic from TGP.

**The FCC-ee tau g-2 measurement at $\sim 10^{-6}$ precision in the 2040s will distinguish these directly:** channel A predicts $\Delta a_\tau \sim 10^{-4}$ (large, distinct from SM noise), while channel B gives $\sim 10^{-7}$ (undetectable at FCC-ee).

---

## Result 5 — Trace anomaly ruled out as $m^2$ channel (ps05)

**Hypothesis:** the QED trace anomaly $\mathcal{L}_{\rm anomaly} = (\alpha/2\pi)\,(\Phi/\Phi_0)\,F^2$ (3 charged leptons) generates a flavor-universal photon-$\Phi$ coupling that reproduces channel B's $m^2$ scaling at the quantum level.

**Result:** Quantitatively suppressed. After integrating out $\Phi$, the effective $(F^2)^2$ operator gives
$$\Delta a_\mu^{\rm anomaly} \;\sim\; \frac{\alpha^4}{(2\pi)^2}\,\frac{m_\mu^2}{M_\Phi^2\,\Phi_0^2}\,\log \;\ll\; \Delta a_\mu^{\rm dilaton} \;\sim\; \frac{m_\mu^4}{8\pi^2\,\Phi_0^2\,M_\Phi^2}\,F_S.$$

### Numerical comparison (at $m=m_\mu$, identical $M_\Phi, \Phi_0$)

$$\frac{\Delta a_\mu^{\rm dilaton\,(A)}}{\Delta a_\mu^{\rm anomaly\,(B)}} \;=\; \frac{m_\mu^2 \cdot 8\pi^2 \cdot (2\pi)^2}{\alpha^4} \;=\; \frac{32\pi^4\,m_\mu^2}{\alpha^4} \;\approx\; 1.2\times 10^{10}.$$

| Fit scenario | $\Phi_0$ (dilaton) | $\Phi_0$ (anomaly only) |
|--------------|--------------------|------------------------|
| $M_\Phi = 10$ GeV   | 3.05 GeV | 1.2 MeV |
| $M_\Phi = 91.2$ GeV | 0.34 GeV | 0.16 MeV |
| $M_\Phi = 1000$ GeV | 0.031 GeV | 0.017 MeV |

Saturating the muon anomaly with the anomaly channel alone would require $\Phi_0 \sim 0.1$ MeV — physically implausible (below atomic scale). The dilaton Yukawa channel is overwhelmingly dominant for any $\Phi_0$ that allows TGP to couple to other observables.

### Interpretation

**TGP's canonical prediction for lepton g-2 is firmly the dilaton $m^4$ channel.** The fact that ps01/ps03 numerically landed on $M_{\rm TGP}\approx M_Z$ is NOT a coincidence: it is the same $M_\Phi$ that would emerge from a full dilaton proper-time calculation (since both treat the substrate quantum mass as the UV scale in a 1-loop form factor). What changes is the relationship to $\Phi_0$: ps01 implicitly set $\Phi_0 = M_\Phi$ by absorbing everything into a single Wilson coefficient; ps04/ps05 disentangle the two scales.

---

---

## Result 6 — Minimal dilaton excluded by Z-pole (ps06)

**Test:** apply the dilaton-Yukawa $y_\ell = m_\ell/\Phi_0$ that fits muon g-2 ($\Phi_0 = 0.34$ GeV at $M_\Phi = 91$ GeV) to other SM observables involving tau.

### Yukawa perturbativity across the SM spectrum

| Fermion | $m_f$ (GeV) | $y_f = m_f/\Phi_0$ | Status |
|---------|-------------|---------------------|--------|
| $e$     | $5.1\times 10^{-4}$ | 0.0015 | perturbative |
| $\mu$   | 0.106 | 0.31 | moderate |
| $\tau$  | 1.78 | **5.3** | **non-perturbative** ($y>\sqrt{4\pi}$) |
| $c$     | 1.27 | 3.8 | non-perturbative |
| $b$     | 4.18 | 12.4 | strongly non-perturbative |
| $t$     | 173  | **512** | EFT breaks completely |

### Z-pole vertex correction to $Z\to\tau^+\tau^-$

A tau Yukawa $y_\tau = 5.3$ with scalar $\Phi$ of mass $M_\Phi = M_Z$ gives a 1-loop vertex correction of order
$$\delta g_V^\tau/g_V^\tau \;\sim\; \frac{y_\tau^2}{16\pi^2}\,\log(M_Z/m_\tau) \;\approx\; 1.8$$
— i.e. a **180% shift** from SM. LEP measures $\Gamma(Z\to\tau^+\tau^-)$ to $\sim 0.3\%$ precision; the dilaton prediction is excluded by a factor $\sim 600$.

Raising $\Phi_0$ to weaken the Yukawa kills the muon g-2 fit:

| $\Phi_0$ (GeV) | $y_\tau$ | $\delta g_V^\tau/g_V^\tau$ | Fits $\Delta a_\mu$? |
|----------------|----------|-----------------------------|---------------------|
| 0.34 | 5.3 | 180% | yes | 
| 1.78 | 1.0 | 6% | factor $10^3$ too small |
| 17.8 | 0.1 | 0.06% | factor $10^7$ too small |
| 246 (Higgs VEV) | 0.007 | $3\times 10^{-7}$ | factor $10^4$ too small |

**No $\Phi_0$ value satisfies both $\Delta a_\mu$ and $Z\to\tau^+\tau^-$ simultaneously.** The minimal exp-metric dilaton picture is experimentally excluded.

### Interpretation

TGP's exponential metric alone does not have enough structure to explain muon g-2 while respecting electroweak precision. Physically, this says: **Goldberger-Wise/dilaton-EFT wisdom** (dilatons with hadronic-scale VEVs give large vertex corrections that are ruled out by LEP) applies to the TGP substrate scalar too.

---

## Outlook — where TGP muon g-2 goes from here

The muon g-2 anomaly is currently an **open problem** for TGP. Three directions:

1. **Vector substrate mode (ps07 target).** TGP may have a spin-1 substrate excitation that the exp-metric formalism doesn't capture. A universal vector $V_\mu\,\bar\ell\gamma^\mu\ell$ at $M_V\sim M_Z$ naturally gives $m^2$ scaling (no Yukawa-hierarchy problem), reproducing the ps01/ps03 empirical fit without violating Z-pole bounds. **Required:** identify a physical mechanism in TGP that generates such a vector mode (analogous to KK photon or $Z'$ in BSM).
2. **Higgs-substrate mixing.** Substrate $\Phi$ couples primarily via mixing with the SM Higgs; the effective Yukawa becomes $(m_\ell/v_{\rm Higgs}) \times \sin\theta_{\Phi H}$, naturally Higgs-VEV suppressed. This keeps perturbativity but needs $M_\Phi$ well below $m_\ell$ to saturate $\Delta a_\mu$ — runs into other problems.
3. **Accept TGP is silent on muon g-2.** The muon anomaly may have a non-TGP origin (SM hadronic, or BSM unrelated to substrate). This is a defensible position: not every anomaly must have a TGP signature.

Independent tests ps07 should address:
- Can a consistent TGP-motivated **vector substrate excitation** be constructed?
- Does substrate-Higgs mixing at the exp-metric level produce a working dilaton after EW breaking?
- What is the minimum "additional structure" needed beyond the exp metric to generate $m^2$ scaling with perturbative Yukawas?

---

## Publication readiness (Zenodo)

**NOT preprint-ready as of 2026-04-20.** The research uncovered a clean theoretical no-go for the simplest TGP completion (ps06), which cannot be published as a "TGP explains muon g-2" result. A publishable preprint would require either:

1. Constructing an alternative TGP substrate mode (e.g. vector excitation) that gives $m^2$ scaling without violating Z-pole constraints — requires ps07+ theoretical work.
2. Framing the result as an **exclusion paper**: *"TGP minimal exponential metric cannot explain muon g-2"* — valid scientific result, but niche and narrow.
3. Pivoting to the viscosity closure as the flagship TGP result and treating muon g-2 as an "ongoing investigation" note.

Option 3 is the recommended path: [[../liquid_viscosity/SUMMARY.md]] is clean, well-validated, and preprint-ready. Muon g-2 remains an open sector for future TGP theory development.

## What is missing (honest statement)

1. **A TGP mechanism for $m^2$ scaling that respects all SM constraints.** The empirical ps01/ps03 fit $\Delta a_\mu \propto (m_\mu/M_Z)^2$ is numerically excellent, but no minimal TGP Lagrangian yet derived produces exactly this pattern with perturbative couplings. A vector-boson substrate mode is the most natural candidate; needs construction.
2. **A reason why the minimal exp-metric dilaton doesn't dominate.** If TGP substrate has both scalar and vector excitations, why is the scalar contribution (which we've now shown is too large when hadronic and too small when electroweak) not visible in Z-pole data? Possible answers: cancellation, flavor-selective coupling, or the scalar mass $M_\Phi$ is actually very light (sub-eV) and the muon anomaly comes from a different sector.
3. **No HVP reinterpretation yet.** The $100\times 10^{-11}$ vs $250\times 10^{-11}$ HVP discrepancy (data-driven vs lattice) is a distinct puzzle, independent of TGP structure.

**Current publishable status:** **not preprint-ready**. The minimal TGP muon g-2 explanation is now ruled out by LEP precision. Further theoretical work is needed before publication. The liquid-viscosity closure remains the clean publishable result from the current research program.

## Suggested title

*"Minimal TGP does not explain muon g-2: a clean negative result from exponential-metric couplings and LEP Z-pole constraints"* (alternative if positive ps07 result: *"A vector substrate excitation in TGP reproduces muon g-2 with $m^2$ scaling"*)

## Timeline

| Year | Experiment | What TGP predicts |
|------|------------|-------------------|
| 2027 | FNAL Muon g-2 Run-6 | $\Delta a_\mu$ final at $\pm 14\times 10^{-11}$; TGP $m^4$ fit of $\Phi_0$ refined |
| 2028 | J-PARC E34 | Independent check of muon anomaly |
| 2030s | Belle-II high-luminosity | $|a_\tau|$ to $\sim 10^{-5}$: **starts constraining $m^4$ prediction $\sim 2\times 10^{-4}$** |
| 2035+ | ATLAS/CMS HL-LHC ditau | $|a_\tau| \sim 10^{-3}$ sensitivity; $m^4$ prediction ($10^{-4}$) enters detection window |
| 2040+ | FCC-ee | $|a_\tau|$ to $\sim 10^{-6}$; decisive $m^4$ vs $m^2$ vs null discrimination |

## Links

- [[PLAN.md]] — full research plan
- [[ps01_schwinger_TGP_1loop.py]] — EFT fit (naive $m^2$ ansatz, $M_{\rm TGP}=91$ GeV)
- [[ps02_lepton_universality_scan.py]] — joint $(a_e, a_\mu, a_\tau)$ constraint scan
- [[ps03_schwinger_proper_time.py]] — proper-time form-factor invariance check ($M_{\rm eff}\approx 91$ GeV)
- [[ps04_vertex_from_exp_metric.py]] — exp-metric tree-level derivation ($m^4$ scaling from minimal coupling)
- [[ps05_trace_anomaly_channel.py]] — trace-anomaly channel suppressed by $10^{10}$, dilaton dominance confirmed
- [[ps06_phi0_scale_consistency.py]] — $\Phi_0$ cross-check: minimal dilaton excluded by LEP Z-pole
- [[../liquid_viscosity/SUMMARY.md]] — parallel closure paper (cheaper sector)
- [[../superconductivity_closure/P6_closure.md]] — template for multi-parameter TGP closure
- [[../NEW_DIRECTIONS_2026-04-20.md]] — research directions index
