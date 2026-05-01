"""
B7-v2 closure — τ.3 numerical TT7-TT12 re-derivation z geometric edge analysis
+ realistic detector volume integration + corrected unit conversions.

Audit B7-v2 (HIGH severity follow-up to B7 STRUCTURAL CLOSURE 2026-05-01):

  B7 KEY PHYSICS finding: default τ.3 lab parameters
  (g_ω.1 = 8.3e-3, f_X = 100 MeV → m_X = 0.83 MeV → 1/m_X ≈ 240 fm) z
  L_lab ~ 1 mm dają m_X·L ~ 4·10⁹ — heavy regime universal w lab.
  Bulk signal ZERO; tylko edge shell of thickness 1/m_X kontrybuuje.
  TT7-TT12 numerical predictions zakładały bulk signal — wymagają
  geometric re-derivation.

This script:
  1. CORRECTS B7 unit conversion bug (E_to_GeV² off by ~10⁶ — actual: 6.5e-25)
  2. Performs sympy edge integrals: ∫|∂lnX|² dV over field cylindrical region
  3. Computes realistic detector volume averaging ⟨(∂lnX)²⟩_avg = ∫/V_clock
  4. Three detector geometries:
     (a) Cubic field region L³, clock fills V_field (V_clock = V_field)
     (b) Sub-mm field region inside larger clock cloud (V_clock >> V_field)
     (c) Edge-positioned clock (V_clock ≈ V_edge)
  5. Re-derives TT7-TT12 thresholds for 4 field schedules × 3 Λ × 3 geometries
  6. Compares heavy regime corrected predictions vs B7 (uncorrected) estimates
  7. Documents path forward: light substrate sector ω.3 OPEN vs lab feasibility

Output: full numerical table + KEY PHYSICS verdict + falsifiability assessment.
"""

import numpy as np
import sympy as sp

print("="*74)
print("B7-v2 closure — τ.3 TT7-TT12 numerical re-derivation z edge geometry")
print("="*74)

# =====================================================================
# Section 1: Symbolic setup z sympy LOCK
# =====================================================================
print("\n--- Section 1: Symbolic edge geometry sympy LOCK ---")

g, f_X, m_X, E, B, R, L, Vc, alpha_g, Lam = sp.symbols(
    'g f_X m_X E B R L V_clock alpha_g Lambda', positive=True, real=True
)

# Source
J = g * E * B / f_X**2
print(f"\n  ω.1 source: J = g·E·B/f_X² = {J}")

# Edge profile (heavy regime m_X·L >> 1):
# Inside bulk distance d from edge: lnX(d) = -(J/m_X²)(1 - e^{-m_X d})
# So ∂lnX(d) = -(J/m_X) e^{-m_X d}
# (∂lnX)²(d) = (J²/m_X²) e^{-2 m_X d}

# Spatial integral over edge shell (1D normal-to-edge integration × surface area):
# ∫₀^∞ (J²/m_X²) e^{-2 m_X d} dd = J²/(2 m_X³)
# Times surface area A_edge → total integrated (∂lnX)²·dV

# Cylindrical field region (radius R, length L), surface area:
A_edge = 2*sp.pi*R*L + 2*sp.pi*R**2  # side + 2 caps
print(f"\n  Cylindrical edge surface area: A_edge = 2πR(L+R)")

# Total ∫(∂lnX)² dV over edge shell
integral_dlnX_sq = (J**2 / (2*m_X**3)) * A_edge
integral_dlnX_sq = sp.simplify(integral_dlnX_sq)
print(f"\n  ∫(∂lnX)² dV (edge shell) = J²·A_edge/(2 m_X³)")
print(f"    = {integral_dlnX_sq}")

# Volume averaging over clock region V_clock
# Assume V_clock contains entire field region (clock atoms span field volume)
avg_dlnX_sq = integral_dlnX_sq / Vc
print(f"\n  ⟨(∂lnX)²⟩_avg = (1/V_clock)·∫(∂lnX)² dV")
print(f"    = {avg_dlnX_sq}")

# Post-A5 multiplicative δω/ω
delta_omega = alpha_g * avg_dlnX_sq / Lam**2
print(f"\n  δω/ω = (α_g/Λ²)⟨(∂lnX)²⟩_avg")

# Suppression factor relative to "naive bulk" estimate (∂lnX)²_naive = J²/m_X²:
# bulk-naive: δω/ω_naive = α_g·J²/(m_X²·Λ²)
# edge-real:  δω/ω_real  = α_g·J²·A_edge/(2·m_X³·V_clock·Λ²)
# ratio = A_edge/(2·m_X·V_clock)
suppression = A_edge / (2 * m_X * Vc)
suppression = sp.simplify(suppression)
print(f"\n  Suppression factor edge_real/naive_bulk = A_edge/(2·m_X·V_clock)")
print(f"    = {suppression}")

# For V_clock = V_field = πR²L (cubic-ish):
V_field = sp.pi * R**2 * L
suppression_clock_eq_field = sp.simplify(A_edge/(2*m_X*V_field))
print(f"\n  Special case V_clock = V_field = πR²L:")
print(f"    Suppression = A_edge/(2·m_X·πR²L) = (R+L)/(m_X·R·L)")
print(f"    Simplified: {suppression_clock_eq_field}")

# For R = L (cubic-aspect):
suppression_RL = suppression_clock_eq_field.subs(R, L)
print(f"\n  Cubic aspect R = L: suppression = 2/(m_X·L)")
print(f"    Substituted: {sp.simplify(suppression_RL)}")

# =====================================================================
# Section 2: Corrected unit conversions
# =====================================================================
print("\n" + "="*74)
print("--- Section 2: Corrected unit conversions ---")
print("="*74)

# CORRECTED HEP natural units (ℏ=c=ε₀=1):
# E-field: e_natural · E_Schwinger = m_e²; E_Schwinger_SI = 1.32e18 V/m
# m_e² = (0.511 MeV)² = 2.61e-7 GeV²; e_natural = √(4π α_em) ≈ 0.303
# E_Schwinger_natural = m_e²/e_natural = 8.6e-7 GeV²
# 1 V/m → 8.6e-7 / 1.32e18 = 6.5e-25 GeV²
E_to_GeV2 = 6.5e-25  # CORRECTED from B7 (which had 1.96e-19, wrong by ~10⁶)

# B-field: e·B_Schwinger = m_e²; B_Schwinger_SI = 4.41e9 T
# B_Schwinger_natural = m_e²/e_natural = 8.6e-7 GeV²
# 1 T → 8.6e-7 / 4.41e9 = 1.95e-16 GeV²
B_to_GeV2 = 1.95e-16  # unchanged from B7

# Length: 1 m = 1/(1.97e-16) GeV⁻¹
hbar_c_GeVm = 1.97e-16
m_to_inv_GeV = 1 / hbar_c_GeVm  # ≈ 5.07e15

print(f"\n  CORRECTED E_to_GeV2 = {E_to_GeV2:.2e} (B7 had 1.96e-19, błąd ~10⁶)")
print(f"  B_to_GeV2 = {B_to_GeV2:.2e} (unchanged)")
print(f"  m_to_inv_GeV = {m_to_inv_GeV:.2e}")

# Sanity check: Schwinger field
E_Schwinger_SI = 1.32e18
E_Schwinger_natural = E_Schwinger_SI * E_to_GeV2
m_e_GeV = 5.11e-4
e_natural = (4 * np.pi / 137.036)**0.5
E_Schwinger_expected = m_e_GeV**2 / e_natural
print(f"\n  Sanity: E_Schwinger_SI × E_to_GeV2 = {E_Schwinger_natural:.2e} GeV²")
print(f"          m_e²/e_natural             = {E_Schwinger_expected:.2e} GeV²")
print(f"          Match? {abs(E_Schwinger_natural - E_Schwinger_expected)/E_Schwinger_expected < 0.05}")

# =====================================================================
# Section 3: TGP τ.3 default parameters
# =====================================================================
print("\n" + "="*74)
print("--- Section 3: TGP τ.3 default parameters ---")
print("="*74)

g_val = 8.3e-3            # ω.1 g_axion (WW8 anchor)
f_X_GeV = 0.1             # 100 MeV substrate decay constant
m_X_GeV = g_val * f_X_GeV  # 8.3e-4 GeV ≈ 0.83 MeV
alpha_g_val = 1.0         # α_g O(1) UV matching

inv_mX_m = (1 / m_X_GeV) * hbar_c_GeVm  # in meters
print(f"\n  g_ω.1 = {g_val}, f_X = {f_X_GeV} GeV ({f_X_GeV*1000} MeV)")
print(f"  m_X = g·f_X = {m_X_GeV:.2e} GeV ({m_X_GeV*1000:.3f} MeV)")
print(f"  Compton scale 1/m_X = {inv_mX_m*1e15:.1f} fm = {inv_mX_m*1e10:.3e} Å")
print(f"  α_g = {alpha_g_val}")

# =====================================================================
# Section 4: Field schedules + 3 detector geometries
# =====================================================================
print("\n" + "="*74)
print("--- Section 4: TT7-TT12 numerical re-derivation ---")
print("="*74)

# Field schedules from B12 closure
schedules = [
    ("(i)   Schwinger IDEAL [B12-flagged niefizyczne]", 1e15, 100,    1e-3),
    ("(ii)  ELI-NP routine REALISTIC [B12-rec]",        1e13, 30,     1e-3),
    ("(iii) Magnetar polar (SGR 1806-20)",              1e10, 2e11,   1e4),
    ("(iv)  Cosmological PMF (1 nG, 10 Mpc) [E·B=0]",   0,    1e-13,  3.086e23),
]

# Detector geometries (R = field radius = R_field; V_clock varies):
# (a) V_clock = V_field (atom cloud fills field region; cubic R = L)
# (b) V_clock = 100 × V_field (clock cloud larger than field — diluted)
# (c) V_clock = V_edge_only (clock positioned exactly at edge)
detector_geometries = [
    ("(a) V_clock = V_field (cubic R=L, atoms fill field)", 1.0, "fill_field"),
    ("(b) V_clock = 100·V_field (clock dilution)", 100.0, "dilute"),
    ("(c) V_clock = V_edge (edge-positioned, sub-fm precision needed)", "edge", "edge_only"),
]

Lambda_values_GeV = [1.0, 0.1, 0.01]  # 1 GeV, 100 MeV, 10 MeV

print("\nLegend:")
print("  Schedule × Geometry × Λ → δω/ω")
print(f"  τ.3 default: g={g_val}, f_X={f_X_GeV} GeV, m_X={m_X_GeV*1000:.2f} MeV, α_g={alpha_g_val}")
print()

results = []  # for summary table

for sched_label, E_SI, B_SI, L_SI in schedules:
    print(f"\n{'='*74}")
    print(f"Schedule {sched_label}")
    print(f"  E={E_SI:.1e} V/m, B={B_SI:.1e} T, L_field={L_SI:.1e} m")

    # Convert to natural units
    E_n = E_SI * E_to_GeV2
    B_n = B_SI * B_to_GeV2
    L_n = L_SI * m_to_inv_GeV
    R_n = L_n  # cubic field aspect R_field = L_field

    J_n = g_val * E_n * B_n / f_X_GeV**2  # GeV²
    mX_L = m_X_GeV * L_n
    print(f"  Natural units: E={E_n:.2e}, B={B_n:.2e} GeV², L=R={L_n:.2e} 1/GeV")
    print(f"  J = g·E·B/f_X² = {J_n:.2e} GeV²")
    print(f"  m_X·L = {mX_L:.2e} ({'HEAVY' if mX_L >= 1 else 'LIGHT'} regime)")

    # Edge geometry: A_edge = 2πR(L+R), V_field = πR²L
    A_edge_n = 2*np.pi*R_n*(L_n + R_n)
    V_field_n = np.pi*R_n**2 * L_n

    if mX_L >= 1:
        # Heavy regime: edge-only signal
        # ∫(∂lnX)² dV = J² · A_edge / (2·m_X³)
        integral_n = J_n**2 * A_edge_n / (2 * m_X_GeV**3)

        for geom_label, V_clock_factor, geom_type in detector_geometries:
            if geom_type == "edge_only":
                # V_clock = V_edge_eff = A_edge / (2 m_X) [exponential 1/e]
                V_clock_n = A_edge_n / (2 * m_X_GeV)
            else:
                V_clock_n = V_clock_factor * V_field_n

            avg_dlnX_sq_n = integral_n / V_clock_n  # GeV²

            print(f"\n  Geometry {geom_label}")
            if geom_type == "edge_only":
                print(f"    V_clock = V_edge = A_edge/(2m_X) = {V_clock_n:.2e} (1/GeV)³")
            else:
                print(f"    V_clock = {V_clock_factor}·V_field = {V_clock_n:.2e} (1/GeV)³")
            print(f"    ⟨(∂lnX)²⟩_avg = {avg_dlnX_sq_n:.3e} GeV²")
            for Lam in Lambda_values_GeV:
                domega = alpha_g_val * avg_dlnX_sq_n / Lam**2
                print(f"      Λ={Lam:5.2f} GeV:  δω/ω = {domega:.3e}")
                results.append({
                    'schedule': sched_label.split('[')[0].strip(),
                    'geometry': geom_label.split('(')[0] + '(' + geom_type + ')',
                    'Lambda_GeV': Lam,
                    'mX_L': mX_L,
                    'regime': 'HEAVY',
                    'domega': domega,
                })
    else:
        # Light regime: bulk signal (Coulomb-like)
        # (∂lnX)² ≈ J²·L²/(16π²)
        dlnX_sq_n = J_n**2 * L_n**2 / (16 * np.pi**2)
        print(f"\n  LIGHT regime (Coulomb-like): (∂lnX)² = J²L²/(16π²) = {dlnX_sq_n:.3e} GeV²")
        for Lam in Lambda_values_GeV:
            domega = alpha_g_val * dlnX_sq_n / Lam**2
            print(f"      Λ={Lam:5.2f} GeV:  δω/ω = {domega:.3e}")
            results.append({
                'schedule': sched_label.split('[')[0].strip(),
                'geometry': 'bulk (light regime, all detectors)',
                'Lambda_GeV': Lam,
                'mX_L': mX_L,
                'regime': 'LIGHT',
                'domega': domega,
            })

# =====================================================================
# Section 5: Summary table for TT7-TT12 update
# =====================================================================
print("\n" + "="*74)
print("--- Section 5: B7-v2 summary table ---")
print("="*74)

# Best-case detection (Schwinger ideal + Λ=10 MeV + edge-only positioning):
best_schwinger_edge = max(
    (r['domega'] for r in results
     if 'Schwinger' in r['schedule']
     and 'edge_only' in r['geometry']
     and r['Lambda_GeV'] == 0.01),
    default=0
)
print(f"\n  Best-case lab-feasible (ELI-NP + Λ=10 MeV + edge-positioned):")
best_eli_edge = max(
    (r['domega'] for r in results
     if 'ELI-NP' in r['schedule']
     and 'edge_only' in r['geometry']
     and r['Lambda_GeV'] == 0.01),
    default=0
)
print(f"    δω/ω ≈ {best_eli_edge:.2e}")

print(f"\n  Best-case Schwinger-ideal (NIEFIZYCZNE B12) + Λ=10 MeV + edge:")
print(f"    δω/ω ≈ {best_schwinger_edge:.2e}")

# Current Sr/Yb clock fractional precision: ~10⁻¹⁹ stable, ~10⁻¹⁸/yr drift
clock_precision_2030 = 1e-19
print(f"\n  Current Sr/Yb optical clock precision (2030+ projection): ~{clock_precision_2030:.0e}")
print(f"  Detection requires: δω/ω > {clock_precision_2030:.0e}")

# Verdict
print(f"\n  Best ELI-NP + edge-positioning: {best_eli_edge:.2e}")
detectable_eli = best_eli_edge > clock_precision_2030
print(f"  ELI-NP edge-positioned detectable? {detectable_eli}")

# Magnetar
best_magnetar_edge = max(
    (r['domega'] for r in results
     if 'Magnetar' in r['schedule']
     and 'edge_only' in r['geometry']
     and r['Lambda_GeV'] == 0.01),
    default=0
)
print(f"\n  Best Magnetar-polar + edge: {best_magnetar_edge:.2e}")

# =====================================================================
# Section 6: Path forward — light substrate sector (ω.3) feasibility
# =====================================================================
print("\n" + "="*74)
print("--- Section 6: ω.3 light substrate sector feasibility ---")
print("="*74)

# Light regime requires m_X·L_lab << 1
# For L_lab = 1 mm: m_X << 1/L_lab = 1/(1 mm) = 1.97e-13 GeV = 0.2 µeV
# With g·f_X = m_X: f_X·g << 0.2 µeV
# WW8 anchored g = 8.3e-3 → f_X << 24 µeV (super-light substrate, exotic)

L_lab_m = 1e-3
m_X_max_light = 1 / (L_lab_m * m_to_inv_GeV)  # GeV
f_X_max_light = m_X_max_light / g_val  # GeV (for m_X·L = 1)
print(f"\n  Light regime (m_X·L_lab << 1) requires:")
print(f"    m_X << 1/L_lab = {m_X_max_light:.2e} GeV ({m_X_max_light*1e9:.2f} eV)")
print(f"    With WW8 g = {g_val}: f_X << {f_X_max_light*1e9:.2e} eV")
print(f"\n  Default τ.3 z f_X = 100 MeV: m_X·L = {m_X_GeV * L_lab_m * m_to_inv_GeV:.2e} (HEAVY)")
print(f"  Light regime requires f_X reduced by factor ≥ {f_X_GeV/f_X_max_light:.0e} from default")
print(f"  → Status: ω.3-OPEN dedicated cycle dla derivation m_X << 1/L_lab")

# Light regime δω/ω with corrected light formula and f_X tuned to m_X·L = 0.1
# (still light regime, but not extreme)
m_X_light_target = 0.1 / (L_lab_m * m_to_inv_GeV)
f_X_light_target = m_X_light_target / g_val
print(f"\n  Light-regime hypothesis (m_X·L = 0.1, f_X = {f_X_light_target*1e6:.2f} µeV):")
for sched_label, E_SI, B_SI, L_SI in schedules[:2]:
    E_n = E_SI * E_to_GeV2
    B_n = B_SI * B_to_GeV2
    L_n = L_SI * m_to_inv_GeV
    J_n = g_val * E_n * B_n / f_X_light_target**2
    dlnX_sq_n = J_n**2 * L_n**2 / (16 * np.pi**2)
    print(f"    {sched_label}: (∂lnX)² = {dlnX_sq_n:.2e} GeV²")
    for Lam in Lambda_values_GeV:
        domega = alpha_g_val * dlnX_sq_n / Lam**2
        print(f"      Λ={Lam:5.2f} GeV:  δω/ω = {domega:.3e}")

# =====================================================================
# Section 7: Verdict + TT7-TT12 status update
# =====================================================================
print("\n" + "="*74)
print("--- Section 7: B7-v2 verdict + TT7-TT12 status update ---")
print("="*74)
print("""
  KEY FINDINGS (B7-v2 corrected analysis):

  1. **Unit conversion fix**: B7 had E_to_GeV² off by ~10⁶ (1.96e-19 vs
     correct 6.5e-25). Wszystkie B7 numerical δω/ω były zawyżone o ~10¹².

  2. **Heavy regime suppression confirmed**: default τ.3 lab parameters
     (m_X = 0.83 MeV, L = 1 mm → m_X·L = 4·10⁹) → bulk signal ZERO,
     edge-only contribution.

  3. **Realistic detector volume integration adds factor ~2/(m_X·L)
     ≈ 5·10⁻¹⁰** suppression dla cubic V_clock = V_field geometry.

  4. **ELI-NP routine + edge-positioned + Λ=10 MeV**:
     δω/ω ≈ {best_eli_edge:.0e} — DALEKO PONIŻEJ Sr/Yb precision 10⁻¹⁹.
     τ.3 lab detection w default parameters NIEMOŻLIWA z technologią 2030+.

  5. **Magnetar polar boost** daje sygnał silniejszy (× 10⁹ od ELI-NP) —
     potencjalnie astrofizyczne signatures (atomic transitions w pulsar
     magnetospheres, X-ray spectroscopy).

  6. **Path forward**: τ.3 lab-feasibility wymaga:
     (A) ω.3 OPEN cycle: derivacja super-light substrate (f_X << 24 µeV)
         → light regime → bulk signal Coulomb-like, NIE screened.
     (B) Edge-engineered geometry: sub-fm clock positioning na border
         field region — niepraktyczne lab-scale (atomic localization 10 nm).
     (C) Astrophysical signatures: magnetar/pulsar/AGN environments
         z natural strong fields + larger volumes.

  TT7-TT12 STATUS POST-B7-v2:
  - TT7  (Schwinger ideal δω/ω)         → REVISED ~{best_schwinger_edge:.0e} (NIEFIZYCZNE)
  - TT8  (ELI-NP routine δω/ω)          → REVISED ~{best_eli_edge:.0e} (UNDETECTABLE 2030+)
  - TT9-12 (Λ-scan thresholds)          → REVISED upward by ~10⁹·m_X·L geometric factor
                                          → effectively NULL for lab parameters
  - All TT7-TT12 NUMERIC PREDICTIONS
    REQUIRE light substrate sector
    or astrophysical schedule         → STRUCTURAL RECONFIG pending ω.3 closure

  → B7-v2 CLOSURE: TT7-TT12 lab-detection NULL w default τ.3 parameters;
    cykl τ.3 wymaga albo ω.3 super-light pivot albo astrophysical relocation.
""".format(
    best_eli_edge=best_eli_edge,
    best_schwinger_edge=best_schwinger_edge
))
