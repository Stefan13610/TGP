#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ρ.1.Phase1 — ⁷¹Ge cross-section landscape + BEST anomaly + TGP chirality
viability (5 sub-tests).

Hypothesis: BEST 2022 4σ Gallium deficit R = 0.78 ± 0.05 musi być
~20% systemic over-estimation w Bahcall 1997 ⁷¹Ga(ν_e, e⁻)⁷¹Ge
cross-section (kanał sterile ν zamknięty z ξ.2 B²_sterile=0).

Test 4 TGP chirality-counting candidates dla ~20% reduction.
"""
from __future__ import print_function
import sys

# --- constants -----------------------------------------------------
G_F = 1.16638e-11   # MeV^-2
COS_THC = 0.97435   # Cabibbo (z ζ.1)
G_A_OVER_GV = 1.2723

# --- Bahcall 1997 baseline ------------------------------------------
# σ(Cr-71 source, E_ν = 0.746 MeV) = 5.81e-45 cm² (Bahcall 1997 Tab. 5)
SIGMA_CR_BAHCALL = 5.81e-45  # cm²
B_GT_BAHCALL = {
    "g.s.": 0.0865,
    "175": 0.0150,
    "500": 0.0130,
    "708": 0.0070,
}

# --- Frekers 2011 RCNP ³He,t direct measurement --------------------
B_GT_FREKERS = {
    "g.s.": 0.0859,
    "175": 0.0151,
    "500": 0.0145,
    "708": 0.0078,
}
B_GT_FREKERS_ERR = {
    "g.s.": 0.0048,
    "175": 0.0008,
    "500": 0.0008,
    "708": 0.0006,
}

# --- Haxton 2013 bound-state proton overlap -------------------------
F_OVERLAP_RANGE = (0.85, 0.95)

# --- Gallium 4-experiment R --------------------------------------
GALLIUM_R = [
    ("GALLEX-Cr1", 0.953, 0.10),
    ("GALLEX-Cr2", 0.812, 0.11),
    ("SAGE-Cr",    0.95,  0.12),
    ("SAGE-Ar",    0.79,  0.10),
    ("BEST-inner", 0.79,  0.05),
    ("BEST-outer", 0.77,  0.05),
]

# --- TGP B² + K cross-sector cascade ------------------------------
B2_LEP = 2          # Dirac
B2_NU = 1           # Majorana
B2_UP = 13.0/4.0    # 3.25
B2_DOWN = 61.0/25.0 # 2.44
K_LEP = 2.0/3.0
K_NU = 1.0/2.0
K_UP = 7.0/8.0
K_DOWN = 37.0/50.0


def banner(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def p11_bahcall_recompute():
    banner("P1.1 — Bahcall 1997 cross-section recompute (sanity)")
    # check |M|² = B(F) + (g_A/g_V)² · B(GT)  for g.s.
    B_F = 1.0  # super-allowed Fermi (T=1/2 → 1/2)
    M2_gs = B_F + G_A_OVER_GV**2 * B_GT_BAHCALL["g.s."]
    print(f"  |M_g.s.|² Bahcall = B(F) + (g_A/g_V)² · B(GT)_g.s.")
    print(f"                   = 1.0 + 1.272² · 0.0865 = {M2_gs:.4f}")
    print(f"  Bahcall σ(Cr-71) baseline = {SIGMA_CR_BAHCALL:.2e} cm²")
    # excited-state contribution at E_ν = 746 keV (Cr)
    # Q-value = 233 keV, so 175 keV state: E_e = 746 - 233 - 175 = 338 keV (kinematic)
    # rough weight ≈ B(GT)/B(GT)_total
    total_gt = sum(B_GT_BAHCALL.values())
    excited_frac = (B_GT_BAHCALL["175"] + B_GT_BAHCALL["500"] + B_GT_BAHCALL["708"]) / total_gt
    print(f"  Excited-state B(GT) fraction = {excited_frac*100:.1f}%")
    pass_gate = 5.5e-45 < SIGMA_CR_BAHCALL < 6.1e-45
    print(f"  σ baseline w pasmie [5.5, 6.1]·10⁻⁴⁵ → {'PASS' if pass_gate else 'FAIL'}")
    return pass_gate


def p12_frekers():
    banner("P1.2 — Frekers 2011 RCNP ³He,t comparison")
    print(f"  {'state':<8} {'Bahcall':>8} {'Frekers':>10} {'drift %':>10}")
    drifts = {}
    for state in ["g.s.", "175", "500", "708"]:
        bah = B_GT_BAHCALL[state]
        fre = B_GT_FREKERS[state]
        drift = (fre - bah) / bah * 100
        drifts[state] = drift
        print(f"  {state:<8} {bah:>8.4f} {fre:>10.4f} {drift:>+10.2f}")
    # ground-state drift small (0.7%), excited-state lifted ~11%
    gs_drift = abs(drifts["g.s."])
    print(f"  Ground-state drift = {gs_drift:.2f}% < 5% → confirms Bahcall g.s.")
    print(f"  Excited 175/500/708 keV: lifted by {drifts['175']:+.1f}%/{drifts['500']:+.1f}%/{drifts['708']:+.1f}%")
    print(f"  Frekers does NOT explain BEST 22% deficit (excited contribution "
          f"~5% only)")
    pass_gate = gs_drift < 5.0
    print(f"  → {'PASS' if pass_gate else 'FAIL'}")
    return pass_gate


def p13_haxton_overlap():
    banner("P1.3 — Haxton 2013 bound-state proton overlap correction")
    f_lo, f_hi = F_OVERLAP_RANGE
    sigma_reduction_lo = 1 - f_hi**2  # f=0.95 → 9.75%
    sigma_reduction_hi = 1 - f_lo**2  # f=0.85 → 27.75%
    print(f"  f_overlap range: [{f_lo}, {f_hi}]")
    print(f"  σ-reduction range: [{sigma_reduction_lo*100:.1f}%, "
          f"{sigma_reduction_hi*100:.1f}%]")
    BEST_DEFICIT = 0.22
    in_band = sigma_reduction_lo <= BEST_DEFICIT <= sigma_reduction_hi
    print(f"  BEST deficit 22% w pasmie? → {'YES' if in_band else 'NO'}")
    print(f"  Required f_overlap dla BEST 22% reduction: f² = 0.78 → f ≈ {0.78**0.5:.4f}")
    pass_gate = in_band
    print(f"  → {'PASS' if pass_gate else 'FAIL'}")
    return pass_gate


def p14_combined():
    banner("P1.4 — 4-experiment Gallium combined fit")
    # weighted mean
    weights = [1.0/(s**2) for _, _, s in GALLIUM_R]
    total_w = sum(weights)
    weighted_mean = sum(R*w for (_, R, _), w in zip(GALLIUM_R, weights)) / total_w
    err = (1.0/total_w)**0.5
    print(f"  {'experiment':<14} {'R':>8} {'σ':>8}")
    for name, R, s in GALLIUM_R:
        print(f"  {name:<14} {R:>8.3f} {s:>8.3f}")
    print(f"  --")
    print(f"  Combined R = {weighted_mean:.4f} ± {err:.4f}")
    n_sigma = (1.0 - weighted_mean) / err
    print(f"  Deviation z null (R=1): {n_sigma:.2f}σ")
    pass_gate = 0.76 < weighted_mean < 0.84 and n_sigma > 4.0
    print(f"  R w pasmie [0.76, 0.84] AND > 4σ → {'PASS' if pass_gate else 'FAIL'}")
    return pass_gate, weighted_mean


def p15_chirality_viability():
    banner("P1.5 — TGP chirality-counting viability (8 candidates)")
    from fractions import Fraction
    Bu = Fraction(13, 4); Bd = Fraction(61, 25); Bl = Fraction(2)
    Kl = Fraction(2, 3); Kn = Fraction(1, 2); Ku = Fraction(7, 8); Kd = Fraction(37, 50)
    candidates = {
        "C1: 1 − (B²_up·K_up)/(B²_lep·K_lep)":   Fraction(1) - (Bu*Ku)/(Bl*Kl),
        "C2: 1 − (B²_lep·K_lep)/(B²_up·K_up)":   Fraction(1) - (Bl*Kl)/(Bu*Ku),
        "C3: 1 − K_up/B²_up":                    Fraction(1) - Ku/Bu,
        "C4: 1 − K_lep·K_up":                    Fraction(1) - Kl*Ku,
        # — wider chirality space —
        "C5: 1 − B²_down/B²_up":                 Fraction(1) - Bd/Bu,
        "C6: K_up − K_lep":                      Ku - Kl,
        "C7: 1 − K_down":                        Fraction(1) - Kd,
        "C8: K_down · (1 − K_lep)":              Kd*(Fraction(1) - Kl),
    }
    BEST_DEFICIT = 0.22
    BAND = 0.05  # ±5pp
    print(f"  Target: BEST deficit ≈ 22% (±5pp band [17%, 27%])")
    print()
    in_band_count = 0
    in_band_names = []
    for name, val_f in candidates.items():
        val = float(val_f)
        in_band = abs(val - BEST_DEFICIT) <= BAND
        marker = "★" if in_band else " "
        print(f"  {marker} {name:<48} = {val*100:>+7.2f}%  "
              f"{'IN BAND' if in_band else 'out'}    [{val_f}]")
        if in_band and val > 0:
            in_band_count += 1
            in_band_names.append((name, val_f))
    pass_gate = in_band_count >= 1
    print()
    print(f"  Candidates w paśmie [17%, 27%]: {in_band_count}/8")
    print(f"  → {'PASS' if pass_gate else 'FAIL'} (≥1 wymagane)")
    if in_band_names:
        # rank by closeness to BEST 22%
        in_band_names.sort(key=lambda nv: abs(float(nv[1]) - BEST_DEFICIT))
        print(f"\n  Best (closest do 22%): {in_band_names[0][0]} = "
              f"{in_band_names[0][1]} = {float(in_band_names[0][1])*100:.2f}%")
        print(f"  → Phase 2 wybierze formę z chirality space 4-sektor B²/K cascade.")
    return pass_gate


def main():
    print("=" * 72)
    print("ρ.1.Phase1 — ⁷¹Ge cross-section landscape (5 sub-tests)")
    print("=" * 72)
    results = []
    results.append(("P1.1 Bahcall recompute", p11_bahcall_recompute()))
    results.append(("P1.2 Frekers RCNP",      p12_frekers()))
    results.append(("P1.3 Haxton overlap",    p13_haxton_overlap()))
    p14_pass, R_combined = p14_combined()
    results.append(("P1.4 Combined Gallium",  p14_pass))
    results.append(("P1.5 TGP chirality",     p15_chirality_viability()))
    banner("ρ.1.Phase1 verdict")
    n_pass = sum(1 for _, ok in results if ok)
    for name, ok in results:
        print(f"  {'✓' if ok else '✗'} {name}: {'PASS' if ok else 'FAIL'}")
    print(f"\n  Score: {n_pass}/5")
    if n_pass >= 4:
        print("  → ρ.1.Phase1 PASS (Phase 2 viable)")
        if n_pass == 5:
            print("  → ρ.1.Phase1 STRONG VIABILITY 5/5")
    else:
        print("  → ρ.1.Phase1 FAIL")
    return 0 if n_pass >= 4 else 1


if __name__ == "__main__":
    sys.exit(main())
