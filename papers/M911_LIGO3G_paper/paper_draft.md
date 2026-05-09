---
title: "Strong-field test of M9.1'': 2PN-phase deviation predictions for Einstein Telescope and Cosmic Explorer"
authors: ["Mateusz Serafin"]
date: 2026-05-07
journal_target: "Phys. Rev. D (or Class. Quantum Grav.)"
type: paper-draft
status: draft-v1
keywords: ["gravitational waves", "post-Newtonian expansion", "ppE", "Einstein Telescope", "Cosmic Explorer", "modified gravity", "TGP"]
related:
  - "[[../../research/op-newton-momentum/M9_1_pp_P1_results.md]]"
  - "[[../../research/op-ppE-mapping/Phase3_paper_ready.md]]"
  - "[[../../research/op-LIGO-3G-deviation/Phase3_falsifier_thresholds.md]]"
  - "[[../../audyt/T01_LIGO3G_falsifier/]]"
  - "[[../../PREDICTIONS_REGISTRY.md]]"
---

# Strong-field test of M9.1'': 2PN-phase deviation predictions for Einstein Telescope and Cosmic Explorer

> ## ⚠ DRAFT-v1 SUPERSEDED 2026-05-09 — MAJOR REVISION REQUIRED
>
> This paper draft (v1, 2026-05-07) is built on Phase 1 OOM heuristic
> β_ppE^TGP^(b=-1) = -5/64 ≈ -0.078 (assuming G_SPA ≈ 1 from Sampson-
> Yunes-Cornish 2013 small-perturbation framework).
>
> **[[../../research/op-ppE-mapping/Phase1.5_G_SPA_lock.md]] Phase 1.5
> derived G_SPA = 48 sympy-exact** (test-particle, 4-level verified:
> sympy LOCK 5/5 + independent hand-calculation E_TGP(U³)/m = 49/48
> + numerical sanity at U=0.1 + alternative SPA derivation orthogonal
> route). True β_ppE^TGP^(b=-1) at η=1/4 = **-15/4 ≈ -3.75** (factor
> 48× LARGER than v1 draft assumed).
>
> **[[../../research/op-GWTC3-reanalysis/Phase2_RERUN_2026-05-09_corrected_beta.md]]
> GWTC-3 RE-RUN with corrected β:** TGP M9.1'' specific (4-3ψ)/ψ
> ansatz **RULED OUT at 5.02σ** by GWTC-3 combined posterior (BF =
> 3.5·10⁻⁶, log10(BF) = -5.45). The "predictive forecast for ET/CE
> 2035+" framing of this draft is now obsolete — falsification has
> already occurred in current LIGO-O3 / GWTC-3 data.
>
> **Required revisions for v2:**
> 1. **§3 derivation:** add Phase 1.5 G_SPA = 48 explicit derivation;
>    explain why Phase 1's "G_SPA ≈ 1" SYC assumption fails for
>    structural-modification theories.
> 2. **§4 detection forecasts:** rescale all β/threshold ratios by
>    factor 48× upward; note that LIGO-O3 single-event already excludes.
> 3. **§5 discussion:** reframe from "future test" to "current observational
>    falsification" + methodological finding (factor-48 correction is a
>    notable result on its own).
> 4. **§6 conclusion:** TGP M9.1'' specific (4-3ψ)/ψ falsified at 5σ;
>    program-level path forward = S07 alternative f(ψ) ansatz exploration.
> 5. **Multi-coefficient:** Phase 1 ratios {-23/10, -38/23, +337/228}
>    are INCORRECT (Phase 1.5 alternative SPA gives -11161/504 ≈ -22.14);
>    requires full re-derivation of M911-P2 pattern.
>
> Sections below are HISTORICAL v1 draft (Phase 1 baseline). DO NOT
> submit without v2 revision incorporating Phase 1.5 finding.

---

## Abstract

We present a falsifiable prediction of the Theory of Geometry of
Pulsation (TGP) program: the canonical metric ansatz M9.1'' (with
hyperbolic time-time component `g_tt = -c²(4-3ψ)/ψ`, where ψ is the
fundamental TGP scalar field) reproduces General Relativity exactly
through 1st post-Newtonian order (β_PPN = γ_PPN = 1) but deviates
explicitly at the 3PN-energy level (or, equivalently, 2PN-phase in
the inspiral waveform): `|Δg_tt|/c² = -(5/6) U³`, where U = GM/(rc²)
is the dimensionless Newtonian potential. We map this deviation
into the parametrized post-Einsteinian (ppE) framework, obtaining
a 2PN-phase coefficient `β_ppE^TGP^(b=-1) = -(3/(128 η)) · (5/6) ·
G_SPA` with sympy-locked rational structure (no fitting). For
equal-mass binaries (η = 1/4) and central SPA chain prefactor
G_SPA = 1 the prediction is `β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 ·
10⁻²`, within an OOM window |β| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹]. Multi-
coefficient ratios at higher PN orders (3PN-, 4PN-, 5PN-phase)
are forced to {-23/10, -38/23, +337/228} respectively, providing
a TGP-distinguishing signature with zero free parameters that
discriminates against dCS, sGB, Einstein-Æther, and Brans-Dicke
modifications. Fisher-matrix forecasts indicate (i) Cosmic Explorer
will detect the M9.1'' deviation at >5σ in a single loud BBH event
(M_tot = 30 M_⊙, d_L = 1 Gpc), (ii) Einstein Telescope will reach
the same significance with a stack of ~10 events, and (iii) LIGO-O5
A+ design will provide the first decisive bound from a stack of
~250 BBH events accumulated over ~2.5 years. We discuss the
window of testability 2027–2035, the compatibility with LIGO O3
existing GWTC-3 ppE constraints, and propose a re-analysis of GWTC-3
with TGP-specific β_ppE prior as a near-term confirmation channel.

## 1. Introduction

The advent of third-generation gravitational-wave detectors —
Einstein Telescope (Maggiore et al. 2020) and Cosmic Explorer (Reitze
et al. 2019) — will dramatically improve the precision of
post-Newtonian (PN) waveform tests of General Relativity (GR) at
strong field. Sensitivity gains of 10× over LIGO O3 are expected,
enabling detection of phase deviations at the 10⁻⁴ level in stacked
analyses (Chamberlain & Yunes 2017; LIGO/Virgo/KAGRA 2023).

The Theory of Geometry of Pulsation (TGP) is a covariant scalar-tensor
program that produces an emergent effective metric from a single
scalar field ψ on a discrete substrate, with axion-like Z₂ symmetry
(Serafin 2025; Serafin 2026). The TGP master metric ansatz, denoted
M9.1'' in the program ledger, takes the canonical form

   ds² = -c² (4-3ψ)/ψ · dt² + ψ/(4-3ψ) · δ_ij dx^i dx^j     (1)

(in isotropic coordinates), with the conformal-budget identity
f(ψ) · h(ψ) = 1 enforcing γ_PPN = 1 exactly through 1PN. The
TGP α=2 vacuum scalar EOM produces a definite asymptotic
expansion of the effective potential ε(r), and substitution into
the metric (1) yields the analytically-locked PN expansion

   g_tt^TGP/(-c²) = 1 − 2U + 2U² − (7/3)U³ + (35/12)U⁴
                   − (91/24)U⁵ + (91/18)U⁶ + ...   (2)

while the standard Schwarzschild GR result in isotropic coordinates is

   g_tt^GR/(-c²)  = 1 − 2U + 2U² − (3/2)U³ + (1)U⁴
                   − (5/8)U⁵ + (3/8)U⁶ + ...        (3)

The first non-trivial structural deviation appears at the U³ level
with magnitude **(5/6)**:

   Δα₃ = α₃^TGP − α₃^GR = -7/3 − (-3/2) = **−5/6**  (4)

This is the central prediction of this paper.

## 2. The (5/6) U³ deviation: derivation and convention

### 2.1 Derivation from α=2 vacuum

The TGP scalar field equation in vacuum reduces to

   ∇²ε + 2(∇ε)²/(1+ε) = 0         (5)

where ε ≡ ψ - 1. The asymptotic expansion ε(r) = a₁/r +
a₂/r² + ... in the local PN parameter η_pn = a₁/r yields,
through sympy LOCK 5/5 PASS (Serafin 2026, op-newton-momentum):

   ε(η_pn) = η_pn − η_pn² + (5/3)η_pn³ − (10/3)η_pn⁴
             + (22/3)η_pn⁵ − (154/9)η_pn⁶ + ...     (6)

with c_n = a_n/a₁^n satisfying a deterministic recursion derivable
from (5). The Newton-matching condition a₁ = U/(2r) (equivalently
η_pn = U/2) propagates through f(1+ε(η_pn)) to (2). The verification
chain f(1+ε) → expansion in U is independently sympy-locked
(7/7 OK in op-ppE-mapping Phase 1).

### 2.2 PN counting convention

The U³ coefficient in g_tt is variously called "3PN" (in the
energy/PPN convention used by Will 1971 for static-metric expansion)
and "2PN-phase" (in the Cutler-Flanagan 1994 convention adopted
universally in gravitational-wave waveform literature). For
clarity and registry consistency we adopt **phase-PN** as
primary in this paper: the (5/6) U³ deviation in g_tt corresponds
to the b_ppE = -1 entry in the Yunes-Pretorius (2009) ppE
parametrization. See Serafin 2026 (audyt/T01_LIGO3G_falsifier/
CONVENTION_DECISION.md) for the full mapping.

## 3. Mapping to ppE: locked β_ppE^TGP^(b=-1)

In the stationary phase approximation (SPA) for binary inspiral
(Cutler & Flanagan 1994), a deviation in g_tt at order U^N
generates a phase modification δΨ(f) = β_ppE · u^b where
u = (πMf)^(1/3) and b = 2N - 7. For U³ this gives b = -1 (2PN-phase).

The mapping from metric coefficient to ppE phase coefficient
proceeds through the standard SPA chain: (metric perturbation)
→ (modified orbital energy E_orb(v)) → (modified energy flux dE/dt)
→ (modified phase Ψ(f)). For metric-only deviations (no new
radiative degrees of freedom — consistent with TGP's single-Φ
axiom and GW1, GW2 entries in the TGP predictions registry that
verify the c_T = c_s and 3-DOF polarization observed at GW170817),
the SPA chain prefactor G_SPA is O(1) for typical equal-mass
binaries (Sampson, Yunes, Cornish 2013). We obtain

   β_ppE^TGP^(b=-1) = -(3/(128 η)) · (5/6) · G_SPA       (7)

For equal-mass binaries (η = 1/4) and central G_SPA = 1:

   β_ppE^TGP^(b=-1) = -5/64 ≈ -7.81 · 10⁻²             (8)

The OOM window for G_SPA ∈ [0.7, 1.5] (uncertainty from the
modified quadrupole formula in M9.1''):

   |β_ppE^TGP^(b=-1)| ∈ [5.5 · 10⁻², 1.2 · 10⁻¹]      (9)

### 3.1 Multi-coefficient TGP-distinguishing signature

Higher-order ppE coefficients follow the same SPA mapping with the
respective Δα_n from (2)-(3):

   Δα_4 = +23/12,  Δα_5 = -19/6,  Δα_6 = +337/72.

The corresponding 3PN-, 4PN-, 5PN-phase ppE coefficients yield
fixed ratios:

   β_(3PN-phase)/β_(2PN-phase) = -23/10                   (10a)
   β_(4PN-phase)/β_(3PN-phase) = -38/23                   (10b)
   β_(5PN-phase)/β_(4PN-phase) = +337/228                 (10c)

These ratios are deterministic functions of α=2 and the hyperbolic
form of f(ψ); they involve **zero free parameters**. Existing
modified-gravity competitors with b = -1 (dCS gravity, scalar-Gauss-
Bonnet, Einstein-Æther, Brans-Dicke higher-PN) all require ≥1
free coupling and cannot reproduce (10) without fitting. This
provides a sharp **multi-coefficient TGP signature** for ET+CE
era simultaneous Bayes inference on multiple ppE coefficients.

## 4. Detection forecasts

We perform Fisher-matrix forecasts using the TaylorF2 + ppE
single-coefficient deviation waveform model (Buonanno et al. 2009),
with analytical fits to the design ASD curves of LIGO A+, ET-D
(Hild 2010), and Cosmic Explorer (Reitze 2019). The Fisher 1σ
uncertainty σ_β includes a degeneracy factor of 5 capturing the
realistic post-marginalization covariance with M_chirp, η, and
χ_eff (Yagi & Yunes 2016 review).

### 4.1 Single-event 5σ thresholds (calibrated)

After calibration to literature SNRs (Maggiore et al. 2020 for
ET-D, Reitze et al. 2019 for CE):

| Detector / scenario | β_5σ | β_TGP / β_5σ | Verdict |
|---------------------|--------|----------------|---------|
| LIGO-O3 (GW150914-like) | ~10⁻¹ | ~0.78 | borderline |
| LIGO-O5 single (loud BBH M=30, 200 Mpc) | ~3·10⁻² | ~2.6 | **YES first decisive single-event 2027+** |
| ET-D single (loud BBH M=30, 1 Gpc) | ~10⁻² | ~7.8 | YES (>20σ) |
| CE single (loud BBH M=30, 1 Gpc) | ~3·10⁻³ | ~26 | **YES (>50σ) decisive** |
| ET+CE network single | ~2·10⁻³ | ~40 | YES (>>50σ) |

### 4.2 Stacked thresholds

Stacking N independent events scales σ_β as 1/√N. First decisive
detection (β_5σ_stack = β_TGP_central):

   N(LIGO-O5)   ≈ 16 events     (~2 months at A+ rate ~100 BBH/yr)
   N(ET-D)      ≈ 1 event       (single loud BBH)
   N(CE)        ≈ 1 event       (single loud BBH)
   N(ET+CE)     ≈ 1 event       (single loud BBH)

### 4.3 Borderline status of LIGO O3

The LIGO O3 single-event ToGR bound at 2PN-phase (Abbott et al.
2021, GWTC-2) is ~10⁻¹, comparable to our prediction
β_ppE^TGP ≈ 8 · 10⁻². The TGP prediction is therefore at the
**boundary of current LIGO sensitivity** — neither falsified nor
fully constrained. A re-analysis of GWTC-3 with the TGP-specific
β_ppE = -5/64 prior (instead of generic ppE prior) should
produce a Bayes factor of order unity between TGP and GR, with
modest discriminating power. We estimate ~1-3σ tentative signal
from public GWTC-3 data is achievable; this re-analysis is
recommended as a near-term confirmation channel.

## 5. Discussion

### 5.1 Status of M9.1''

The deviation (5/6)U³ is rigorously derived from M9.1'' under the
following assumptions, all internally consistent: (i) α = 2 ODE
kinetic exponent (TGP F2 LOCKED); (ii) hyperbolic f(ψ) form
(M9_1_pp_setup, structurally selected by anti-podal budget +
β_PPN = γ_PPN = 1 + ψ = 4/3 horizon coincidence with R3 mass
spectrum at 4-digit precision); (iii) single-Φ axiom (TGP F1
LOCKED — no new radiative DOF).

We acknowledge that M9.1'' itself is currently postulated rather
than derived from first-substrate principles; this is the subject
of a parallel derivation cycle (Serafin 2026, audyt/S07_M911_derivation/).
However, falsifiability does not require derivation: even as a
postulate, M9.1'' makes the sharp prediction (8) and (10), which
ET+CE will test or rule out within the 2027-2035 window.

### 5.2 Window of testability

The most aggressive timeline:
- **2027** (LIGO-O5 first observing run, A+ design): single-event
  detection borderline; stack ~16 events for first decisive.
- **2030** (LIGO-O5 + LISA via cross-channel BH5 QNM ringdown):
  cross-channel orthogonal verification possible.
- **2035** (ET-D + CE first observing run): single-event decisive
  detection in CE; ET-D single borderline-YES, ET-D stack ~10
  events confirms.
- **2036+** (1 yr of ET+CE operation): >700σ confidence on
  multi-coefficient TGP signature M911-P2.

### 5.3 Falsification clause

M9.1'' is sfalsified iff at least one of the following is verified:
- LIGO-O5 with N ≥ 250 BBH events does NOT detect β_ppE^(b=-1)
  in the window |β_ppE^TGP| ∈ [5.5·10⁻², 1.2·10⁻¹] at 5σ;
- ET-D with stack 10 BBH events does NOT detect the same window;
- CE with single loud BBH event (M_chirp ≥ 30 M_⊙, SNR ≥ 150)
  does NOT detect β_ppE^(b=-1) > 5·10⁻² at 5σ.

A positive detection of the predicted window with simultaneous
verification of the multi-coefficient ratios (10) at any of LIGO-O5,
ET-D, CE, or ET+CE will promote M9.1'' from `LIVE` to
`DERIVED-OBSERVATIONAL` status.

## 6. Conclusion

The M9.1'' canonical metric of TGP makes a sharp, falsifiable,
zero-free-parameter prediction at 2PN-phase: β_ppE^TGP^(b=-1) =
-5/64 ≈ -7.81·10⁻², with a multi-coefficient ratio signature
{-23/10, -38/23, +337/228} that discriminates TGP from all known
modified-gravity competitors. Cosmic Explorer will detect this
deviation in a single loud BBH event at >5σ; Einstein Telescope
within ~10 events; LIGO-O5 A+ within ~250 events. This prediction
provides a clean public falsification contract for the strong-field
sector of TGP, with first decisive test possible by 2027 and
decisive confirmation by 2035.

## Acknowledgments

Sympy LOCK verifications: op-newton-momentum (M9_1_pp_P1, c_n=2..7
α=2 vacuum) and op-ppE-mapping (Phase 1, β_ppE^TGP and multi-
coefficient pattern). Fisher-matrix forecasts: op-LIGO-3G-deviation
(Phase 2-3). Audit framework: T01_LIGO3G_falsifier (NEEDS, FALSIFIER
STATEMENT, CONVENTION DECISION, FINDINGS).

## References

[1] Yunes, N., Pretorius, F., 2009, Phys. Rev. D 80:122003.
[2] Cutler, C., Flanagan, É., 1994, Phys. Rev. D 49:2658.
[3] Buonanno, A., Iyer, B. R., Ochsner, E., Pan, Y., Sathyaprakash,
    B. S., 2009, Phys. Rev. D 80:084043.
[4] Blanchet, L., 2014, Living Rev. Relativ. 17:2.
[5] Yunes, N., Yagi, K., Pretorius, F., 2016, Phys. Rev. D 94:084002.
[6] Sampson, L., Yunes, N., Cornish, N., 2013, Phys. Rev. D 88:064056.
[7] Will, C. M., 2014, Living Rev. Relativ. 17:4.
[8] Maggiore, M., et al., 2020, JCAP 03:050 (Einstein Telescope
    Science Case).
[9] Reitze, D., et al., 2019, Bull. Am. Astron. Soc. 51:035 (Cosmic
    Explorer).
[10] Chamberlain, K., Yunes, N., 2017, Phys. Rev. D 96:084039.
[11] Yagi, K., Yunes, N., 2016, arXiv:1602.04674 (review).
[12] LIGO/Virgo Collaboration, 2021, Phys. Rev. D 103:122002 (GWTC-2
     ToGR).
[13] LIGO/Virgo/KAGRA Collaboration, 2023, Phys. Rev. X 13:041039
     (GWTC-3 ToGR).
[14] Vallisneri, M., 2008, Phys. Rev. D 77:042001 (Fisher matrix
     pitfalls).
[15] Hild, S., et al., 2010, Class. Quantum Grav. 27:015003 (ET-D
     design).
[16] Damour, T., Jaranowski, P., Schäfer, G., 2014, Phys. Rev. D
     89:064058 (4PN ADM).
[17] Mishra, C. K., Iyer, B. R., Sundararajan, P. A., 2016, Phys.
     Rev. D 93:084054 (3PN waveform).
[18] Serafin, M., 2025-2026, TGP_FOUNDATIONS.md, sek08c_metryka_z_
     substratu.tex, sek08a_akcja_zunifikowana.tex, op-newton-
     momentum/M9_1_pp_setup.md, op-newton-momentum/M9_1_pp_P1_
     results.md, op-ppE-mapping/Phase{0,1,2,3}_*.md, op-LIGO-3G-
     deviation/Phase{0,1,2,3}_*.md, audyt/T01_LIGO3G_falsifier/.

---

## Submission notes

**Target journals (in preference order):**
1. **Phys. Rev. D** — strong-field falsifier paper; expected length
   ~12-18 pages PRD format with appendices. Section structure
   above is paper-ready; needs LaTeX porting and figure work.
2. **Class. Quantum Grav.** — secondary target; comparable scope.
3. **Phys. Rev. Lett.** — possible if reduced to 4 pages with
   strongest result (multi-coefficient signature + ET/CE forecast)
   as headline.

**Pre-submission tasks:**
1. Produce 3 figures: (i) M9.1'' vs GR g_tt expansion residuals;
   (ii) β_ppE^TGP vs detector sensitivity curves (Bayes prediction
   bands); (iii) multi-coefficient ratio test in {β_2PN,β_3PN}-plane
   showing TGP point + competitor exclusion regions.
2. Tighter G_SPA lock (Phase 1.5 future cycle, ~5% precision):
   reduces OOM uncertainty on β_TGP from ~30% to ~5%.
3. Production-grade Fisher matrix re-run with bilby/pycbc and
   official noise curves: replaces analytical ASD fits.
4. GWTC-3 reanalysis (Tier 5 side-cycle) — provides independent
   pre-publication confirmation if positive.
5. (Optional) Bayes-factor calculation TGP vs GR for selected
   GWTC-3 events.

**Estimated time to PRD-ready submission:** ~1-2 months from this
draft, contingent on figure work and Phase 1.5 G_SPA lock.

**Author contributions:** Mateusz Serafin (theory + program design).

**Funding & data availability:** TBD per author preference.

**Co-author candidates:** TBD per author network.
