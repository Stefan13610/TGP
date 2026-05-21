---
title: "Phase 0 — Balance composite Higgs substrate attempt sesja-1-of-N"
date: 2026-05-18
parent: "[[./README.md]]"
type: phase-zero-balance
phase: 0
status: 🟢 ACTIVE
sympy_substance_ratio: "6 FP / 1 LIT / 1 DEC = 75% FP (cycle 1/2/7 STRICT pattern)"
hardcoded_T_pass: "1 (T8 DEC only — declarative substance budget)"
---

# Phase 0 — Balance composite Higgs substrate attempt

## §1 — Inputs

### §1.1 — TGP-native constants (LIVE)

| Constant | Value | Source | Status |
|---|---|---|---|
| m_X (substrate kink scale) | 60 MeV | L06 NUMERICAL ANCHOR | inherited |
| m_Pl (Planck mass) | 1.221·10¹⁹ GeV | PDG | external constant |
| H_0 (Hubble) | 67.4 km/s/Mpc → 1.4·10⁻³³ eV | PDG | external constant |
| m_ν (heaviest ν) | ~0.1 eV | PDG ν₃ | external constant |
| m_e (electron) | 0.511 MeV | PDG | external constant |
| **v_H (SM Higgs VEV)** | **246.22 GeV** | **PDG** | **TARGET** dla scale matching |
| α (fine structure) | 1/137.036 | PDG | external constant |
| G_F (Fermi constant) | 1.1664·10⁻⁵ GeV⁻² | PDG | dla cross-checks |

### §1.2 — TGP minimal axioms (LIVE LOCK)

- **S05:** single complex Φ field, 2 real DoF
- **Z₂:** Goldstone-Nambu symmetry (cycle L07 closure 2026-05-16)
- **U(1):** gauge from Φ phase rotation (standard QFT)
- **RP² topology:** for fermion sektor (cycle 1, 2 sesji 2026-05-17)

### §1.3 — Cycle 6 ruled-out paths (LIVE LOCK; NIE re-attempt)

| Path | Approach | Failure |
|---|---|---|
| α | Berry × spinor → SU(2) | RP² has 2 invariants; SU(2) needs 3 generators |
| β | π_n(RP²) higher homotopy | Invariants WITHIN gauge groups, not emergence |
| γ | Φ-Φ* doublet → SU(2) | TGP 2 real DoF vs SU(2) doublet 4 real DoF |
| δ | S05+Z₂ → emergent gauge | 1 continuous symmetry vs SM EW 4 generators |

**This cycle attempts 5th path: indirect emergence via composite (Higgs as bound state of strong dynamics).**

### §1.4 — Methodology inputs

- Sympy structural analysis (NOT numerical loop computation; that comes w sesjach 2+)
- Symmetry-counting (Goldstone theorem)
- Scale enumeration (log-space analysis of TGP-native scale combinations)
- S05 compatibility check (no new axioms allowed)

## §2 — Outputs

### §2.1 — Phase 1 deliverables (sesja-1 scope)

- **T1 LIT:** Composite Higgs literature anchors (3 sources + features inventory)
- **T2 FP:** TGP-native scale combination analysis (which combinations approach TeV?)
- **T3 FP:** Candidate "technicolor-like" TGP confining dynamics enumeration
- **T4 FP:** Goldstone counting (broken substrate symmetry → ≥4 Goldstones?)
- **T5 FP:** Hierarchy mechanism (m_H << Λ via composite analog)
- **T6 FP:** S05+Z₂+U(1) compatibility — does composite mechanism need NEW axiom?
- **T7 FP:** Sesja-1 verdict per pre-registered decision tree
- **T8 DEC:** S05 preservation declaration

### §2.2 — Phase FINAL deliverable

- Sesja-1 verdict (A- / B+ PARTIAL / HALT-B / HALT-A)
- Multi-session campaign positioning (sesja-1 of estimated 6-8)
- Updated path enumeration: 4 paths (α/β/γ/δ ruled out) + 1 composite (this sesja) → 5-path map
- Future sesji-2+ scope sketch if not HALT-B

## §3 — Risk register

| Risk | Severity | Mitigation |
|---|---|---|
| **R1** Multi-session scope discipline | HIGH | Sesja-1 narrow focus per §0.2 README; deferrals explicit |
| **R2** TeV scale absence w TGP | HIGH | T2 substantive scan; HALT-B akceptowalne wynik |
| **R3** Composite Higgs needs hidden gauge group | HIGH | T6 explicit check; if requires new axiom → HALT-B |
| **R4** Indirect emergence different from direct cycle 6 | medium | Honest disclosure: composite ≠ direct; different obstruction modes possible |
| **R5** Honest HALT-B akceptowalne | low | Pre-registered probabilities w README §0.3 (~30% HALT-B) |
| **R6** Sesja-1 B+ realistic | low | Pre-registered probabilities (~50% B+) |
| **R7** Post-hoc threshold adjustment | medium | Pre-registered thresholds strict enforcement |
| **R8** Hardcoded T_pass=True drift (sesja 2026-05-17 R1 lesson) | medium | Cycle 1/2/7 STRICT pattern: only T8 DEC hardcoded; all FP conditional |

## §4 — 8/8 gate

- [x] **G1:** Phase 0 balance sheet exists (this file)
- [x] **G2:** README z BINDING contract + pre-registered falsifier (§0.3)
- [x] **G3:** Sympy substance plan: 6 FP + 1 LIT + 1 DEC = 75% FP
- [x] **G4:** 0 hardcoded T_pass=True dla FP tests (only T8 DEC budget; cycle 1/2/7 strict)
- [x] **G5:** Pre-flight read confirmation §0.5 (8 sources)
- [x] **G6:** TGP-native check Q1-Q8 §0.4
- [x] **G7:** P-requirements 6/6 declared README
- [x] **G8:** Risk register §3 z 8 risks honestly identified

**Verdict:** 🟢 **8/8 PASS.** Cycle ready dla Phase 1 sympy.

## §5 — Scope discipline

**IN-SCOPE sesja-1:**
- Literature anchor extraction (T1)
- TGP-native scale enumeration (T2)
- Candidate confining dynamics identification (T3)
- Goldstone counting (T4)
- Hierarchy mechanism check (T5)
- S05 compatibility check (T6)
- Sesja-1 verdict per pre-registered tree (T7)
- S05 preservation declaration (T8)

**OUT-OF-SCOPE sesja-1 (deferred sesji 2+):**
- EWPO precision corrections (S, T, U parameters)
- Walking technicolor + fermion mass generation
- LHC Higgs coupling distinguishability tests
- Numerical Λ_compositeness prediction (full quantitative)
- Vacuum misalignment angle θ determination
- Strong-dynamics non-perturbative calculation
- Lattice methods for confining dynamics
- Specific UV-completion details

**Multi-session campaign plan (preliminary, sesja-1 may revise):**
- Sesja-1 (this): structural framework attempt + scale matching analysis
- Sesja-2: IF sesja-1 B+ — deepen mechanism w specific candidate (e.g., kink condensate)
- Sesja-3: Goldstone identification + Higgs doublet emergence (assuming sesja-2 progress)
- Sesja-4: W/Z mass generation via composite Higgs mechanism
- Sesja-5: EWPO compatibility check + walking technicolor candidate
- Sesja-6: Quantitative match to v_H = 246 GeV (if path viable)
- Sesja-7-8: Polish + alternative paths if obstructions

**Plan revisable per sesja-1 outcome:** HALT-B w sesji-1 → campaign closes 1-of-1 (composite ruled out); B+ → revise plan w sesji-2.

---

**Sign-off:** Claudian @ 2026-05-18 (sesja-1-of-N composite Higgs multi-session campaign).
