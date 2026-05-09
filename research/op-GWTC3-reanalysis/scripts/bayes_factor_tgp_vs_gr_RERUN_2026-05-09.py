"""
op-GWTC3-reanalysis Phase 2 RE-RUN with Phase 1.5 corrected β
=============================================================

Date: 2026-05-09
Purpose: Verify the Phase 1.5 critical finding (G_SPA = 48 not ≈ 1) by
re-running the original Phase 2 GWTC-3 Bayes factor calculation with the
corrected β_ppE^TGP^(b=-1) value.

Methodology IDENTICAL to original Phase 2 (bayes_factor_tgp_vs_gr.py):
  - Same published GWTC-3 ToGR δφ̂_4 posteriors (Abbott 2023 PRX 13:041039).
  - Same conversion β_ppE^(b=-1) = (3/(128 η)) · α̂_4^GR · δφ̂_4 = 4.336 · δφ̂_4 at η=1/4.
  - Same Laplace approximation Bayes factor.

Only change: β_TGP_central from -5/64 ≈ -7.81e-2 (Phase 1 OOM heuristic)
              to -15/4 ≈ -3.75 (Phase 1.5 sympy-exact, test-particle).

Compare verdicts side-by-side to see if Phase 1.5 finding implies TGP
RULED OUT or TGP TENTATIVELY DETECTED.

References: Phase1.5_G_SPA_lock.md, Phase2_Bayes_factor.md (original).
"""

import numpy as np

# ============================================================================
# §1 — Same published GWTC-3 ToGR 2PN-phase bounds (Abbott 2023)
# ============================================================================
GWTC3_combined_median = 0.05
GWTC3_combined_90CL = 0.30
GWTC3_combined_sigma = GWTC3_combined_90CL / 1.645  # ≈ 0.182

GW150914_median = 0.2
GW150914_90CL = 1.0
GW150914_sigma = GW150914_90CL / 1.645  # ≈ 0.608

GWTC2_combined_median = 0.10
GWTC2_combined_90CL = 0.45
GWTC2_combined_sigma = GWTC2_combined_90CL / 1.645  # ≈ 0.274

# ============================================================================
# §2 — TGP M911-P1 prediction: PHASE 1 vs PHASE 1.5
# ============================================================================
# Phase 1 OOM (heuristic, G_SPA ≈ 1 from SYC 2013 small-perturbation regime):
beta_TGP_Phase1 = -5.0 / 64.0  # ≈ -7.81e-2

# Phase 1.5 LOCKED (G_SPA = 48 sympy-exact, test-particle limit):
beta_TGP_Phase15 = -15.0 / 4.0  # = -3.75 exact

print("=" * 80)
print("op-GWTC3-reanalysis Phase 2 RE-RUN with Phase 1.5 corrected β")
print("Date: 2026-05-09")
print("=" * 80)
print()
print("§1 — Comparison Phase 1 vs Phase 1.5 β_ppE^TGP^(b=-1):")
print(f"  Phase 1 OOM (G_SPA ≈ 1):       β_TGP = {beta_TGP_Phase1:.6f} = -5/64")
print(f"  Phase 1.5 LOCKED (G_SPA = 48): β_TGP = {beta_TGP_Phase15:.6f} = -15/4")
print(f"  Ratio: Phase 1.5 / Phase 1 = {beta_TGP_Phase15 / beta_TGP_Phase1:.2f}× (factor 48)")
print()

# ============================================================================
# §3 — Conversion δφ̂_4 ↔ β_ppE^(b=-1) (same as Phase 2)
# ============================================================================
def alpha_hat_4_GR(eta):
    return 15293365/508032 + 27145/504 * eta + 3085/72 * eta**2

eta_eq = 0.25
alpha_4_GR_eq = alpha_hat_4_GR(eta_eq)
prefactor = 3 / (128 * eta_eq)
conversion_factor = prefactor * alpha_4_GR_eq

# Phase 1 prediction in LIGO units:
delta_phi4_TGP_Phase1 = beta_TGP_Phase1 / conversion_factor

# Phase 1.5 prediction in LIGO units:
delta_phi4_TGP_Phase15 = beta_TGP_Phase15 / conversion_factor

print("§2 — Conversion + TGP δφ̂_4 prediction in LIGO ToGR units:")
print(f"  α̂_4^GR(η=1/4) = {alpha_4_GR_eq:.4f}")
print(f"  Conversion: β_ppE = {conversion_factor:.4f} · δφ̂_4")
print(f"  Phase 1:    δφ̂_4_TGP = {delta_phi4_TGP_Phase1:.6f}")
print(f"  Phase 1.5:  δφ̂_4_TGP = {delta_phi4_TGP_Phase15:.6f}  (factor 48 larger)")
print()

# ============================================================================
# §4 — TGP overlay vs published GWTC-3 posteriors
# ============================================================================
print("=" * 80)
print("§3 — TGP overlay vs GWTC-3 (with Phase 1.5 corrected β)")
print("=" * 80)

scenarios = [
    ("GWTC-3 combined (~90 BBH)", GWTC3_combined_median, GWTC3_combined_sigma),
    ("GW150914-like single", GW150914_median, GW150914_sigma),
    ("GWTC-2 combined", GWTC2_combined_median, GWTC2_combined_sigma),
]

print()
print(f"PHASE 1 (β = {beta_TGP_Phase1:.4f}, δφ̂_TGP = {delta_phi4_TGP_Phase1:.4f}):")
print(f"{'Scenario':<32}{'observed':>12}{'σ_GW':>10}{'TGP':>12}{'(TGP-obs)/σ':>14}{'Verdict':>20}")
print("-" * 100)
for name, observed, sigma in scenarios:
    dev_Phase1 = (delta_phi4_TGP_Phase1 - observed) / sigma
    verdict_Phase1 = "CONSISTENT" if abs(dev_Phase1) < 1 else f"{abs(dev_Phase1):.2f}σ tension"
    print(f"{name:<32}{observed:>12.4f}{sigma:>10.4f}{delta_phi4_TGP_Phase1:>12.4f}"
          f"{dev_Phase1:>14.2f}{verdict_Phase1:>20}")
print()

print(f"PHASE 1.5 (β = {beta_TGP_Phase15:.4f}, δφ̂_TGP = {delta_phi4_TGP_Phase15:.4f}):")
print(f"{'Scenario':<32}{'observed':>12}{'σ_GW':>10}{'TGP':>12}{'(TGP-obs)/σ':>14}{'Verdict':>20}")
print("-" * 100)
for name, observed, sigma in scenarios:
    dev_Phase15 = (delta_phi4_TGP_Phase15 - observed) / sigma
    if abs(dev_Phase15) < 1:
        verdict = "CONSISTENT"
    elif abs(dev_Phase15) < 2:
        verdict = f"{abs(dev_Phase15):.2f}σ tension"
    elif abs(dev_Phase15) < 3:
        verdict = f"{abs(dev_Phase15):.2f}σ TENSION"
    elif abs(dev_Phase15) < 5:
        verdict = f"{abs(dev_Phase15):.2f}σ STRONG TENSION"
    else:
        verdict = f"{abs(dev_Phase15):.2f}σ RULED OUT"
    print(f"{name:<32}{observed:>12.4f}{sigma:>10.4f}{delta_phi4_TGP_Phase15:>12.4f}"
          f"{dev_Phase15:>14.2f}{verdict:>20}")
print()

# ============================================================================
# §5 — Bayes factor TGP vs GR (Laplace approximation), both phase versions
# ============================================================================
def compute_log_BF(observed, sigma, delta_phi_TGP):
    """Laplace approximation log Bayes factor TGP vs GR."""
    chi2_TGP = ((observed - delta_phi_TGP) / sigma) ** 2
    chi2_GR = ((observed - 0.0) / sigma) ** 2
    return (chi2_GR - chi2_TGP) / 2.0

print("=" * 80)
print("§4 — Bayes factor TGP vs GR (Laplace approximation)")
print("=" * 80)
print()

print(f"{'Scenario':<32}{'BF Phase 1':>13}{'log10(BF) Ph1':>15}{'BF Phase 1.5':>16}{'log10(BF) Ph1.5':>17}")
print("-" * 100)
for name, observed, sigma in scenarios:
    log_BF_P1 = compute_log_BF(observed, sigma, delta_phi4_TGP_Phase1)
    log_BF_P15 = compute_log_BF(observed, sigma, delta_phi4_TGP_Phase15)

    BF_P1 = np.exp(log_BF_P1)
    BF_P15 = np.exp(log_BF_P15)

    print(f"{name:<32}{BF_P1:>13.4f}{log_BF_P1/np.log(10):>15.3f}"
          f"{BF_P15:>16.6e}{log_BF_P15/np.log(10):>17.3f}")
print()

# Jeffreys interpretation
def jeffreys_verdict(log10_BF):
    if log10_BF > 1: return "STRONG TGP (>10× preferred)"
    elif log10_BF > 0.5: return "weak TGP preference"
    elif log10_BF > -0.5: return "INCONCLUSIVE"
    elif log10_BF > -1: return "weak GR preference"
    elif log10_BF > -2: return "STRONG GR (>10× preferred)"
    else: return f"OVERWHELMING GR (>10^{abs(int(log10_BF))} preferred)"

print()
print(f"{'Scenario':<32}{'Phase 1 verdict':<32}{'Phase 1.5 verdict':<35}")
print("-" * 100)
for name, observed, sigma in scenarios:
    log_BF_P1 = compute_log_BF(observed, sigma, delta_phi4_TGP_Phase1)
    log_BF_P15 = compute_log_BF(observed, sigma, delta_phi4_TGP_Phase15)
    v1 = jeffreys_verdict(log_BF_P1/np.log(10))
    v15 = jeffreys_verdict(log_BF_P15/np.log(10))
    print(f"{name:<32}{v1:<32}{v15:<35}")
print()

# ============================================================================
# §6 — Detection power: σ-significance of Phase 1.5 finding
# ============================================================================
print("=" * 80)
print("§5 — Significance of Phase 1.5 vs Phase 1 in current GWTC-3 data")
print("=" * 80)

print()
print("PHASE 1 prediction (β = -5/64):")
obs_p = GWTC3_combined_median
sig_p = GWTC3_combined_sigma
print(f"  GWTC-3 σ_combined ≈ {sig_p:.3f}")
print(f"  TGP δφ̂_4 = {delta_phi4_TGP_Phase1:.4f}")
print(f"  TGP/σ ratio = {abs(delta_phi4_TGP_Phase1)/sig_p:.3f} → ~{abs(delta_phi4_TGP_Phase1)/sig_p:.1f}σ marginal")
print(f"  σ-tension from observed: {abs((delta_phi4_TGP_Phase1 - obs_p)/sig_p):.2f}σ")
print(f"  → CONSISTENT (Phase 2 verdict 2026-05-07)")
print()

print("PHASE 1.5 prediction (β = -15/4):")
print(f"  GWTC-3 σ_combined ≈ {sig_p:.3f}")
print(f"  TGP δφ̂_4 = {delta_phi4_TGP_Phase15:.4f}")
print(f"  TGP/σ ratio = {abs(delta_phi4_TGP_Phase15)/sig_p:.3f} → would be ~{abs(delta_phi4_TGP_Phase15)/sig_p:.1f}σ signal")
print(f"  σ-tension from observed: {abs((delta_phi4_TGP_Phase15 - obs_p)/sig_p):.2f}σ")
log_BF_P15 = compute_log_BF(obs_p, sig_p, delta_phi4_TGP_Phase15)
print(f"  log10(BF_TGP/GR) = {log_BF_P15/np.log(10):.3f}")
print(f"  → STRONG GR preference (TGP RULED OUT at ~{abs((delta_phi4_TGP_Phase15 - obs_p)/sig_p):.1f}σ)")
print()

# ============================================================================
# §7 — Implication for TGP M9.1''
# ============================================================================
print("=" * 80)
print("§6 — IMPLICATION for TGP M9.1''")
print("=" * 80)
print()

dev_Phase15_combined = abs((delta_phi4_TGP_Phase15 - GWTC3_combined_median) / GWTC3_combined_sigma)

if dev_Phase15_combined > 5:
    print(f"  GWTC-3 combined RULES OUT Phase 1.5 corrected β at {dev_Phase15_combined:.1f}σ")
    print(f"  in GENERIC ToGR multi-coefficient marginalized analysis.")
    print()
    print("  Three possibilities:")
    print()
    print("  (A) TGP M9.1'' is FALSIFIED by current GWTC-3 data.")
    print("      → Phase 1.5 result correct + TGP M9.1'' ansatz observationally")
    print("        excluded. Need to either modify the ansatz (e.g., explore")
    print("        different f(ψ) structures) or abandon M9.1''.")
    print()
    print("  (B) Phase 1.5 has a derivation error (despite 5/5 LOCK PASS).")
    print("      → Should re-derive G_SPA via independent route (e.g., direct")
    print("        v² → ω² → x conversion in different gauge, full DJS 2-body,")
    print("        or independent literature cross-check).")
    print()
    print("  (C) The conversion β_ppE ↔ δφ̂_4 has unrecognized normalization.")
    print("      → If LIGO ToGR's δφ̂_4 normalization differs from Buonanno-Iyer")
    print("        2009 by some factor, the bound interpretation changes.")
    print("        Must verify Abbott 2023 §IV explicitly.")
    print()
    print("  Recommended: PAUSE downstream propagation. User decision needed on")
    print("  which possibility to investigate first (likely A, with B/C as cross-checks).")
elif dev_Phase15_combined > 2:
    print(f"  GWTC-3 combined STRONG TENSION with Phase 1.5 at {dev_Phase15_combined:.1f}σ")
elif dev_Phase15_combined > 1:
    print(f"  GWTC-3 combined mild tension with Phase 1.5 at {dev_Phase15_combined:.1f}σ")
else:
    print(f"  GWTC-3 combined CONSISTENT with Phase 1.5 at {dev_Phase15_combined:.1f}σ")
print()

print("=" * 80)
print("RE-RUN COMPLETE 2026-05-09")
print("=" * 80)
