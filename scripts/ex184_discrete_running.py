#!/usr/bin/env python3
"""
ex184_discrete_running.py
Sesja v45, 2026-04-05

HIPOTEZA: Jeśli Φ₀ = N_f², to formuła α_s = N_c³·g₀^e/(8·N_f²)
naturalnie koduje DYSKRETNY RUNNING — wartość α_s skacze na progach
flavorowych (m_c, m_b, m_t) przez zmianę N_f.

TEST:
  1. α_s^TGP(N_f) = 27·g₀^e/(8·N_f²) dla N_f = 3,4,5,6
  2. Porównanie z QCD running (1-loop, 2-loop) na progach flavorowych
  3. Ocena: czy TGP reprodukuje running bez ciągłego RGE?

PROGI FLAVOROWE:
  N_f=3: μ ∈ [Λ_QCD, m_c]  → próg przy m_c ≈ 1.27 GeV
  N_f=4: μ ∈ [m_c, m_b]    → próg przy m_b ≈ 4.18 GeV
  N_f=5: μ ∈ [m_b, m_t]    → próg przy M_Z = 91.19 GeV (reference)
  N_f=6: μ ∈ [m_t, ∞)      → próg przy m_t ≈ 172.7 GeV
"""
import numpy as np

PHI = (1 + np.sqrt(5)) / 2
N_c = 3

# === INPUTS ===
g0_e = 0.86941       # from phi-FP (ex180)
ALPHA_S_MZ = 0.1179  # PDG 2024, at M_Z
ALPHA_S_ERR = 0.0009
M_Z = 91.1876        # GeV

# Quark masses (PDG 2024, MSbar)
m_c = 1.27    # GeV (charm threshold)
m_b = 4.18    # GeV (bottom threshold)
m_t = 172.69  # GeV (top threshold)

print("=" * 72)
print("ex184: Discrete running hypothesis — α_s(N_f) = N_c³·g₀^e/(8·N_f²)")
print("=" * 72)

# ===== 1. TGP PREDICTIONS =====
print("\n--- 1. TGP predictions: α_s(N_f) = 27·g₀^e/(8·N_f²) ---\n")

nf_values = [3, 4, 5, 6]
alpha_tgp = {}

print(f"  {'N_f':>3s}  {'N_f²':>4s}  {'α_s^TGP':>10s}  {'threshold':>12s}")
print("  " + "-" * 40)

thresholds = {3: f"< m_c ({m_c} GeV)", 4: f"m_c–m_b",
              5: f"m_b–m_t (M_Z)", 6: f"> m_t ({m_t} GeV)"}

for nf in nf_values:
    a_s = N_c**3 * g0_e / (8 * nf**2)
    alpha_tgp[nf] = a_s
    print(f"  {nf:3d}  {nf**2:4d}  {a_s:10.5f}  {thresholds[nf]:>12s}")

# ===== 2. QCD RUNNING (1-loop) =====
print("\n--- 2. QCD running (1-loop) for comparison ---\n")

def beta0(nf):
    """1-loop QCD beta function coefficient."""
    return (11 * N_c - 2 * nf) / (12 * np.pi)

def alpha_s_1loop(mu, mu0, alpha0, nf):
    """1-loop running: 1/α(μ) = 1/α(μ₀) + β₀·ln(μ²/μ₀²)"""
    b0 = beta0(nf)
    return alpha0 / (1 + alpha0 * b0 * 2 * np.log(mu / mu0))

def beta1(nf):
    """2-loop QCD beta function coefficient."""
    return (102 - 38 * nf / 3) / (16 * np.pi**2)

# Run from M_Z downward through thresholds
print("  Running α_s(M_Z) = 0.1179 through flavor thresholds (1-loop):\n")

# M_Z → m_b (N_f=5)
alpha_at_mb_5 = alpha_s_1loop(m_b, M_Z, ALPHA_S_MZ, 5)
# Match at m_b: α_s continuous, N_f: 5→4
alpha_at_mb_4 = alpha_at_mb_5  # matching condition (1-loop)
# m_b → m_c (N_f=4)
alpha_at_mc_4 = alpha_s_1loop(m_c, m_b, alpha_at_mb_4, 4)
# Match at m_c: N_f: 4→3
alpha_at_mc_3 = alpha_at_mc_4
# m_c → 1 GeV (N_f=3) — representative low scale
alpha_at_1gev = alpha_s_1loop(1.0, m_c, alpha_at_mc_3, 3)

# M_Z → m_t (N_f=5→6)
alpha_at_mt_5 = alpha_s_1loop(m_t, M_Z, ALPHA_S_MZ, 5)
alpha_at_mt_6 = alpha_at_mt_5  # matching condition

# M_Z → 500 GeV (N_f=6)
alpha_at_500 = alpha_s_1loop(500, m_t, alpha_at_mt_6, 6)

print(f"  Scale       N_f  α_s(QCD 1-loop)  α_s^TGP(N_f)   ratio")
print("  " + "-" * 62)

# The question: at WHAT scale does α_s^QCD = α_s^TGP(N_f)?
# For N_f=5 (M_Z): α_s^TGP = 0.1174, α_s^QCD = 0.1179 ✓
# For N_f=3,4,6: compare at thresholds

comparisons = [
    ("m_c (1.27)", 3, alpha_at_mc_3),
    ("m_b (4.18)", 4, alpha_at_mb_4),
    ("M_Z (91.2)", 5, ALPHA_S_MZ),
    ("m_t (172.7)", 6, alpha_at_mt_6),
]

for label, nf, alpha_qcd in comparisons:
    a_tgp = alpha_tgp[nf]
    ratio = a_tgp / alpha_qcd
    dev = (ratio - 1) * 100
    print(f"  {label:<13s} {nf:3d}  {alpha_qcd:14.5f}  {a_tgp:12.5f}  {ratio:.4f} ({dev:+.1f}%)")

# ===== 3. SZCZEGÓŁOWA ANALIZA PROGOWA =====
print("\n--- 3. Detailed threshold analysis ---\n")

print("  Pytanie: Na jakiej skali μ* wartość α_s^QCD(μ*) = α_s^TGP(N_f)?")
print("  (czyli: gdzie TGP 'spoczywa' w danym reżimie flavorowym?)\n")

def find_matching_scale(alpha_target, mu_lo, mu_hi, mu_ref, alpha_ref, nf, tol=1e-6):
    """Find μ* where α_s(μ*) = alpha_target using bisection."""
    for _ in range(100):
        mu_mid = np.sqrt(mu_lo * mu_hi)  # geometric mean
        alpha_mid = alpha_s_1loop(mu_mid, mu_ref, alpha_ref, nf)
        if abs(alpha_mid - alpha_target) < tol:
            return mu_mid
        if alpha_mid > alpha_target:
            mu_lo = mu_mid  # need higher μ (lower α)
        else:
            mu_hi = mu_mid  # need lower μ (higher α)
    return mu_mid

# N_f = 5: we know it's near M_Z
mu_star_5 = find_matching_scale(alpha_tgp[5], 10, 500, M_Z, ALPHA_S_MZ, 5)
print(f"  N_f=5: α_s^TGP = {alpha_tgp[5]:.5f} → μ* = {mu_star_5:.1f} GeV (cf. M_Z = {M_Z:.1f})")

# N_f = 4: between m_c and m_b
mu_star_4 = find_matching_scale(alpha_tgp[4], m_c, m_b, m_b, alpha_at_mb_4, 4)
print(f"  N_f=4: α_s^TGP = {alpha_tgp[4]:.5f} → μ* = {mu_star_4:.2f} GeV (cf. m_c={m_c}, m_b={m_b})")

# N_f = 3: below m_c
# α_s^TGP(3) = 27*0.86941/72 = 0.3260
# This is large — need to check if 1-loop is reliable
if alpha_tgp[3] < alpha_at_mc_3:
    mu_star_3 = find_matching_scale(alpha_tgp[3], 0.5, m_c, m_c, alpha_at_mc_3, 3)
    print(f"  N_f=3: α_s^TGP = {alpha_tgp[3]:.5f} → μ* = {mu_star_3:.2f} GeV (cf. m_c={m_c})")
else:
    # α_s^TGP > α_s at m_c → TGP value is below m_c
    # Run further down
    mu_star_3 = find_matching_scale(alpha_tgp[3], 0.3, m_c, m_c, alpha_at_mc_3, 3)
    print(f"  N_f=3: α_s^TGP = {alpha_tgp[3]:.5f} → μ* = {mu_star_3:.2f} GeV (poniżej m_c!)")
    print(f"         UWAGA: 1-loop może nie być wiarygodny przy tak niskiej skali")

# N_f = 6: above m_t
if alpha_tgp[6] < alpha_at_mt_6:
    mu_star_6 = find_matching_scale(alpha_tgp[6], m_t, 10000, m_t, alpha_at_mt_6, 6)
    print(f"  N_f=6: α_s^TGP = {alpha_tgp[6]:.5f} → μ* = {mu_star_6:.0f} GeV (cf. m_t={m_t})")
else:
    print(f"  N_f=6: α_s^TGP = {alpha_tgp[6]:.5f} > α_s(m_t) = {alpha_at_mt_6:.5f}")
    print(f"         → TGP value corresponds to scale BELOW m_t (inconsistent)")

# ===== 4. GEOMETRIC MEAN HYPOTHESIS =====
print("\n--- 4. Geometric mean hypothesis ---\n")
print("  Jeśli TGP daje α_s na 'środkowej' skali danego reżimu N_f,")
print("  to μ* powinno być średnią geometryczną granic:\n")

geo_means = {
    3: ("Λ_QCD~0.3", "m_c", np.sqrt(0.3 * m_c)),
    4: ("m_c", "m_b", np.sqrt(m_c * m_b)),
    5: ("m_b", "m_t", np.sqrt(m_b * m_t)),
    6: ("m_t", "~1TeV", np.sqrt(m_t * 1000)),
}

print(f"  {'N_f':>3s}  {'granice':>15s}  {'μ_geo':>8s}  {'μ*_TGP':>8s}  {'ratio':>7s}")
print("  " + "-" * 50)

for nf in nf_values:
    lo_label, hi_label, mu_geo = geo_means[nf]
    if nf == 3:
        mu_star = mu_star_3
    elif nf == 4:
        mu_star = mu_star_4
    elif nf == 5:
        mu_star = mu_star_5
    elif nf == 6:
        mu_star = mu_star_6 if alpha_tgp[6] < alpha_at_mt_6 else float('nan')

    ratio = mu_star / mu_geo if not np.isnan(mu_star) else float('nan')
    print(f"  {nf:3d}  {lo_label}–{hi_label:>5s}  {mu_geo:8.2f}  {mu_star:8.2f}  {ratio:7.3f}")

# ===== 5. ALTERNATIVE: α_s AT THRESHOLD DIRECTLY =====
print("\n--- 5. α_s at threshold vs TGP (direct comparison) ---\n")

# PDG reference values from running (approximate, well-known)
# α_s(m_τ) ≈ 0.330 ± 0.014 (PDG, from τ decay)
# α_s(m_b) ≈ 0.225 (1-loop from M_Z)
# α_s(m_Z) = 0.1179 ± 0.0009
# α_s(m_t) ≈ 0.108

alpha_pdg_refs = {
    "m_τ (1.78)": (0.330, 0.014, 3),
    "m_c (1.27)": (alpha_at_mc_3, None, 3),   # 1-loop estimate
    "m_b (4.18)": (alpha_at_mb_4, None, 4),    # 1-loop estimate
    "M_Z (91.2)": (ALPHA_S_MZ, ALPHA_S_ERR, 5),
    "m_t (173)":  (alpha_at_mt_6, None, 6),
}

print(f"  {'scale':>12s}  {'N_f':>3s}  {'α_s(QCD)':>10s}  {'α_s^TGP':>10s}  {'dev':>8s}  {'comment':>10s}")
print("  " + "-" * 65)

for label, (a_qcd, a_err, nf) in alpha_pdg_refs.items():
    a_tgp = alpha_tgp[nf]
    dev = (a_tgp / a_qcd - 1) * 100
    if a_err:
        sig = abs(a_tgp - a_qcd) / a_err
        comment = f"{sig:.1f}σ"
    else:
        comment = "1-loop"
    print(f"  {label:>12s}  {nf:3d}  {a_qcd:10.5f}  {a_tgp:10.5f}  {dev:+7.1f}%  {comment:>10s}")

# ===== 6. RATIO TEST: α_s(N_f)/α_s(N_f') =====
print("\n--- 6. Ratio test: scaling with N_f ---\n")
print("  TGP predicts: α_s(N_f)/α_s(N_f') = (N_f'/N_f)²")
print("  This is a PARAMETER-FREE prediction (g₀^e cancels)!\n")

# QCD ratios (1-loop) at thresholds
ratio_pairs = [
    (5, 4, ALPHA_S_MZ, alpha_at_mb_4, "M_Z vs m_b"),
    (5, 3, ALPHA_S_MZ, alpha_at_mc_3, "M_Z vs m_c"),
    (5, 6, ALPHA_S_MZ, alpha_at_mt_6, "M_Z vs m_t"),
    (4, 3, alpha_at_mb_4, alpha_at_mc_3, "m_b vs m_c"),
]

print(f"  {'pair':>12s}  {'(N_f\'/N_f)²':>10s}  {'α_QCD ratio':>12s}  {'dev':>8s}")
print("  " + "-" * 50)

for nf1, nf2, a1, a2, label in ratio_pairs:
    tgp_ratio = (nf1 / nf2)**2  # TGP prediction: α(nf2)/α(nf1) = (nf1/nf2)²
    qcd_ratio = a2 / a1
    dev = (tgp_ratio / qcd_ratio - 1) * 100
    print(f"  {label:>12s}  {tgp_ratio:10.4f}  {qcd_ratio:12.4f}  {dev:+7.1f}%")

# ===== 7. PDG τ-DECAY VALUE =====
print("\n--- 7. Comparison with α_s(m_τ) from τ decay ---\n")

alpha_tau_pdg = 0.330
alpha_tau_err = 0.014
alpha_tgp_3 = alpha_tgp[3]

print(f"  α_s(m_τ, PDG) = {alpha_tau_pdg} ± {alpha_tau_err}")
print(f"  α_s^TGP(N_f=3) = {alpha_tgp_3:.5f}")
dev_tau = (alpha_tgp_3 / alpha_tau_pdg - 1) * 100
sig_tau = abs(alpha_tgp_3 - alpha_tau_pdg) / alpha_tau_err
print(f"  Odchylenie: {dev_tau:+.1f}%")
print(f"  sigma: {sig_tau:.1f}")
print()

if sig_tau < 2:
    print("  *** α_s^TGP(N_f=3) jest ZGODNY z α_s(m_τ) w granicach 2σ! ***")
    print("  To jest NIEZALEŻNY TEST — m_τ nie był używany w kalibracji!")
elif sig_tau < 3:
    print("  α_s^TGP(N_f=3) jest marginalnie zgodny z α_s(m_τ) (2-3σ)")
else:
    print("  α_s^TGP(N_f=3) NIE jest zgodny z α_s(m_τ) (>3σ)")

# ===== 8. COMPLETE FORMULA WITH N_f =====
print("\n--- 8. Complete TGP formula ---\n")

print("  Jeśli Φ₀ = N_f²:")
print()
print("  +---------------------------------------------------------+")
print("  |                                                         |")
print("  |  α_s(N_f) = N_c³ · g₀^e / (8 · N_f²)                  |")
print("  |           = 27 · g₀^e / (8 · N_f²)                     |")
print("  |                                                         |")
print("  |  N_f=3 (m_τ scale):  α_s = 3g₀^e/8  = {:.5f}        |".format(alpha_tgp[3]))
print("  |  N_f=4 (m_c–m_b):   α_s = 27g₀^e/128 = {:.5f}       |".format(alpha_tgp[4]))
print("  |  N_f=5 (M_Z scale):  α_s = 27g₀^e/200 = {:.5f}       |".format(alpha_tgp[5]))
print("  |  N_f=6 (>m_t):      α_s = 3g₀^e/32  = {:.5f}        |".format(alpha_tgp[6]))
print("  |                                                         |")
print("  |  Ratio (parameter-free): α(N_f)/α(N_f') = (N_f'/N_f)² |")
print("  |                                                         |")
print("  +---------------------------------------------------------+")

# ===== 9. INVERSE: g₀^e FROM α_s(m_τ) =====
print("\n--- 9. Inverse test: g₀^e from α_s(m_τ) ---\n")

g0_from_tau = alpha_tau_pdg * 8 * 9 / 27  # N_f=3, N_c=3
g0_from_MZ = ALPHA_S_MZ * 8 * 25 / 27      # N_f=5

print(f"  g₀^e z α_s(m_τ, N_f=3): {g0_from_tau:.5f}")
print(f"  g₀^e z α_s(M_Z, N_f=5): {g0_from_MZ:.5f}")
print(f"  g₀^e z φ-FP (r₂₁):      {g0_e:.5f}")
print()
dev_tau_g = (g0_from_tau / g0_e - 1) * 100
dev_MZ_g = (g0_from_MZ / g0_e - 1) * 100
print(f"  dev(m_τ): {dev_tau_g:+.2f}%")
print(f"  dev(M_Z): {dev_MZ_g:+.2f}%")
print()
print(f"  Średnia ważona g₀^e: {(g0_from_tau + g0_from_MZ)/2:.5f}")

# ===== WNIOSKI =====
print("\n" + "=" * 72)
print("WNIOSKI ex184")
print("=" * 72)
print(f"""
  1. DISCRETE RUNNING HYPOTHESIS:
     α_s(N_f) = N_c³·g₀^e/(8·N_f²) = 27·g₀^e/(8·N_f²)

     N_f=3: {alpha_tgp[3]:.5f} vs PDG(m_τ) {alpha_tau_pdg} ± {alpha_tau_err} → {sig_tau:.1f}σ
     N_f=5: {alpha_tgp[5]:.5f} vs PDG(M_Z) {ALPHA_S_MZ} ± {ALPHA_S_ERR} → {abs(alpha_tgp[5]-ALPHA_S_MZ)/ALPHA_S_ERR:.1f}σ

  2. PARAMETER-FREE RATIOS:
     α(N_f=3)/α(N_f=5) = (5/3)² = 25/9 = 2.778
     (testowalne BEZ znajomości g₀^e!)

  3. INTERPRETACJA:
     Φ₀ = N_f² oznacza, że substrat 'widzi' liczbę aktywnych flavorów.
     Formuła α_s = N_c³g₀^e/(8N_f²) łączy:
       - N_c = 3 (grupa cechowania)
       - N_f (aktywne flavory — zależy od skali!)
       - g₀^e (coupling leptonu z φ-FP)

  4. KLUCZOWE PYTANIE:
     Czy N_f-zależność Φ₀ jest fizyczna (discrete running)
     czy Φ₀=25 jest STAŁĄ (N_f=5 zamrożone na M_Z)?
     Test: porównanie z α_s(m_τ) rozstrzyga.
""")
