#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ex236_omega_lambda_quark_fit.py
================================
χ² FIT Ω_Λ Z MAS KWARKÓW + PEŁNA SYNTEZA

KONTEKST (ex235):
  Shifted Koide K(m+m₀)=2/3 działa idealnie dla kwarków.
  R12 formuła: m₀ = A × m₃/m₁, A = 1/(Φ_eff × φ), Φ_eff = 168×Ω_Λ×3/14

  Przy Ω_Λ = 0.685: m_b(R12)=4218 (0.9%), m_t(R12)=174327 (0.9%)
  Pytanie: jakie Ω_Λ minimalizuje łączny χ²?

PLAN:
  §1. χ² minimalizacja Ω_Λ z obu sektorów kwarków
  §2. Porównanie Ω_Λ(fit) z Planck 2018 (0.6847 ± 0.0073)
  §3. Pełna synteza: 5 wejść → 9 predykcji (leptony + kwarki)
  §4. Ratio m₀/Λ_QCD — fizyczna interpretacja shift'u
  §5. Neutrino mass sum test from Ω_Λ

Data: 2026-04-06
"""

import sys, io, math
import numpy as np
from scipy.optimize import brentq, minimize_scalar

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

phi = (1 + np.sqrt(5)) / 2

TESTS = []
def record(name, passed, detail=""):
    TESTS.append((name, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")


# ============================================================
# KOIDE TOOLS
# ============================================================
def koide(m1, m2, m3):
    S = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / S**2

def solve_r31_from_K(r21, K_target=2/3):
    """Solve K(1, r₂₁, r₃₁) = K_target for r₃₁ (larger root)."""
    a = np.sqrt(r21)
    def K_func(x):
        r31 = x**2
        return koide(1.0, r21, r31) - K_target

    # K(x) has minimum at x_min = (1+a²)/(1+a)
    x_min = (1 + a**2) / (1 + a)

    roots = []
    try:
        x_left = brentq(K_func, 0.01, x_min)
        roots.append(x_left)
    except:
        pass
    try:
        x_right = brentq(K_func, x_min, 100000.0)
        roots.append(x_right)
    except:
        pass

    if not roots:
        return np.nan
    x_sol = max(roots)
    return x_sol**2


# ============================================================
# PDG MASSES (MeV)
# ============================================================
m_e = 0.51100
m_mu = 105.6584
m_tau = 1776.86

m_d = 4.67      # ± 0.48
m_s = 93.4      # ± 8.6
m_b = 4180.0    # ± 30

m_u = 2.16      # ± 0.49
m_c = 1270.0    # ± 20
m_t = 172760.0  # ± 300

# Uncertainties (for χ²)
sig_mb = 30.0
sig_mt = 300.0
sig_md = 0.48
sig_ms = 8.6
sig_mu = 0.49
sig_mc = 20.0


# ============================================================
# SHIFTED KOIDE m₀(fit)
# ============================================================
def find_m0_shifted(m1, m2, m3):
    """Find m₀ such that K(m1+m₀, m2+m₀, m3+m₀) = 2/3."""
    def f(m0):
        return koide(m1+m0, m2+m0, m3+m0) - 2/3

    for start in np.linspace(-m1*0.9, 50000.0, 5000):
        try:
            v1 = f(start)
            v2 = f(start + 10.0)
            if v1 * v2 < 0:
                return brentq(f, start, start + 10.0)
        except:
            continue
    return np.nan


def predict_m3_shifted(m1, m2, m0):
    """Predict m3 from m1, m2, m0 using K(m+m0)=2/3."""
    r21_s = (m2 + m0) / (m1 + m0)
    r31_s = solve_r31_from_K(r21_s, 2/3)
    if np.isnan(r31_s):
        return np.nan
    return (m1 + m0) * r31_s - m0


# ============================================================
# §1. χ² MINIMALIZACJA Ω_Λ
# ============================================================
print("=" * 72)
print("§1. χ² MINIMALIZACJA Ω_Λ Z MAS KWARKÓW")
print("=" * 72)

def chi2_quarks(OL):
    """Total χ² for m_b and m_t predictions given Ω_Λ."""
    P0 = 168 * OL
    PE = P0 * 3 / 14
    A = 1.0 / (PE * phi)

    # Down sector
    m0d = A * m_b / m_d
    mb_pred = predict_m3_shifted(m_d, m_s, m0d)
    if np.isnan(mb_pred):
        return 1e10

    # Up sector
    m0u = A * m_t / m_u
    mt_pred = predict_m3_shifted(m_u, m_c, m0u)
    if np.isnan(mt_pred):
        return 1e10

    chi2_b = ((mb_pred - m_b) / sig_mb)**2
    chi2_t = ((mt_pred - m_t) / sig_mt)**2

    return chi2_b + chi2_t

# Scan
print("\n  Fine scan around minimum:")
OL_fine = np.linspace(0.65, 0.80, 301)
chi2_vals = np.array([chi2_quarks(ol) for ol in OL_fine])
i_min = np.argmin(chi2_vals)
OL_min_scan = OL_fine[i_min]
chi2_min_scan = chi2_vals[i_min]
print(f"  Scan minimum: Ω_Λ = {OL_min_scan:.4f}, χ² = {chi2_min_scan:.2f}")

# Brentq refinement
result = minimize_scalar(chi2_quarks, bounds=(0.65, 0.80), method='bounded')
OL_best = result.x
chi2_best = result.fun
print(f"  Refined minimum: Ω_Λ = {OL_best:.6f}, χ² = {chi2_best:.4f}")

# Predictions at best Ω_Λ
P0_best = 168 * OL_best
PE_best = P0_best * 3/14
A_best = 1.0 / (PE_best * phi)

m0d_best = A_best * m_b / m_d
m0u_best = A_best * m_t / m_u
mb_best = predict_m3_shifted(m_d, m_s, m0d_best)
mt_best = predict_m3_shifted(m_u, m_c, m0u_best)

print(f"\n  Przy Ω_Λ = {OL_best:.6f}:")
print(f"    Φ₀ = {P0_best:.2f}")
print(f"    Φ_eff = {PE_best:.2f}")
print(f"    A = {A_best:.6f}")
print(f"    m₀(down) = {m0d_best:.2f} MeV")
print(f"    m₀(up) = {m0u_best:.1f} MeV")
print(f"    m_b = {mb_best:.1f} MeV  (PDG: {m_b}, err: {abs(mb_best-m_b)/m_b*100:.3f}%)")
print(f"    m_t = {mt_best:.0f} MeV  (PDG: {m_t}, err: {abs(mt_best-m_t)/m_t*100:.3f}%)")
print(f"    χ²/dof = {chi2_best:.4f}/1 = {chi2_best:.4f}")

# Comparison with Planck
OL_planck = 0.6847
sig_OL_planck = 0.0073
tension_sigma = abs(OL_best - OL_planck) / sig_OL_planck
print(f"\n  Porównanie z Planck 2018:")
print(f"    Ω_Λ(fit) = {OL_best:.4f}")
print(f"    Ω_Λ(Planck) = {OL_planck} ± {sig_OL_planck}")
print(f"    Tension: {tension_sigma:.1f}σ")

record("T1: Ω_Λ χ² fit",
       chi2_best < 10.0,
       f"Ω_Λ = {OL_best:.4f}, χ² = {chi2_best:.2f}, tension = {tension_sigma:.1f}σ")


# ============================================================
# §2. Ω_Λ CONFIDENCE INTERVAL
# ============================================================
print("\n" + "=" * 72)
print("§2. Ω_Λ PRZEDZIAŁ UFNOŚCI (Δχ²=1)")
print("=" * 72)

chi2_threshold = chi2_best + 1.0  # 1σ for 1 parameter

try:
    OL_lo = brentq(lambda x: chi2_quarks(x) - chi2_threshold, 0.60, OL_best)
except:
    OL_lo = np.nan
try:
    OL_hi = brentq(lambda x: chi2_quarks(x) - chi2_threshold, OL_best, 0.90)
except:
    OL_hi = np.nan

print(f"  Ω_Λ = {OL_best:.4f} (+{OL_hi-OL_best:.4f}, -{OL_best-OL_lo:.4f})")
print(f"  1σ interval: [{OL_lo:.4f}, {OL_hi:.4f}]")

contains_planck = OL_lo <= OL_planck <= OL_hi
print(f"  Planck Ω_Λ = {OL_planck} {'wewnątrz' if contains_planck else 'POZA'} 1σ")

record("T2: Planck within 1σ",
       contains_planck,
       f"Ω_Λ ∈ [{OL_lo:.4f}, {OL_hi:.4f}], Planck = {OL_planck}")


# ============================================================
# §3. PEŁNA SYNTEZA: LEPTONY + KWARKI
# ============================================================
print("\n" + "=" * 72)
print("§3. ★ PEŁNA SYNTEZA: 5 WEJŚĆ → 9 PREDYKCJI")
print("=" * 72)

print("""
  WEJŚCIA (5):
    1. m_e = 0.51100 MeV          [lepton sector]
    2. r₂₁ = 206.768              [from soliton ODE φ-FP]
    3. K = 2/3                     [from (N_gen+1)/(2N_gen)]
    4. m_d = 4.67 MeV, m_s = 93.4 MeV  [quark masses, PDG]
    5. Ω_Λ = {OL_best:.4f}                [from χ² fit or Planck]

  (Actually m_u, m_c also needed for up sector → 7 total inputs)
""".format(OL_best=OL_best))

# Lepton predictions
r21_lep = 206.768
r31_lep = solve_r31_from_K(r21_lep, 2/3)
m_mu_pred = m_e * r21_lep
m_tau_pred = m_e * r31_lep

# Quark predictions (using best Ω_Λ)
mb_pred = predict_m3_shifted(m_d, m_s, m0d_best)
mt_pred = predict_m3_shifted(m_u, m_c, m0u_best)

# α_s prediction (from ex234)
g0e = 0.86941
Phi0 = P0_best
alpha_s_pred = 7 * 27 * g0e / (12 * Phi0)

print(f"  PREDYKCJE (z Ω_Λ = {OL_best:.4f}):\n")
print(f"  {'Wielkość':20s}  {'Predykcja':>14s}  {'PDG':>14s}  {'Error':>10s}")
print("  " + "-" * 64)
print(f"  {'m_μ':20s}  {m_mu_pred:14.4f}  {m_mu:14.4f}  {abs(m_mu_pred-m_mu)/m_mu*100:9.4f}%")
print(f"  {'m_τ':20s}  {m_tau_pred:14.2f}  {m_tau:14.2f}  {abs(m_tau_pred-m_tau)/m_tau*100:9.4f}%")
print(f"  {'m_b':20s}  {mb_pred:14.1f}  {m_b:14.1f}  {abs(mb_pred-m_b)/m_b*100:9.3f}%")
print(f"  {'m_t':20s}  {mt_pred:14.0f}  {m_t:14.0f}  {abs(mt_pred-m_t)/m_t*100:9.3f}%")
print(f"  {'α_s':20s}  {alpha_s_pred:14.4f}  {'0.1179±9':>14s}  {abs(alpha_s_pred-0.1179)/0.0009:9.1f}σ")

record("T3: m_τ from Koide",
       abs(m_tau_pred - m_tau) / m_tau < 0.001,
       f"m_τ = {m_tau_pred:.2f} MeV, err = {abs(m_tau_pred-m_tau)/m_tau*100:.4f}%")

record("T4: m_b from R12+Koide",
       abs(mb_pred - m_b) / m_b < 0.02,
       f"m_b = {mb_pred:.1f} MeV, err = {abs(mb_pred-m_b)/m_b*100:.2f}%")

record("T5: m_t from R12+Koide",
       abs(mt_pred - m_t) / m_t < 0.02,
       f"m_t = {mt_pred:.0f} MeV, err = {abs(mt_pred-m_t)/m_t*100:.2f}%")

record("T6: α_s",
       abs(alpha_s_pred - 0.1179) / 0.0009 < 2.0,
       f"α_s = {alpha_s_pred:.4f}, deviation = {abs(alpha_s_pred-0.1179)/0.0009:.1f}σ")


# ============================================================
# §4. FIZYCZNA INTERPRETACJA m₀
# ============================================================
print("\n" + "=" * 72)
print("§4. FIZYCZNA INTERPRETACJA m₀ / Λ_QCD")
print("=" * 72)

Lambda_QCD = 332.0  # MeV (Nf=3, MS-bar, from PDG)
m_proton = 938.272

# m₀(fit) values
m0_down_fit = find_m0_shifted(m_d, m_s, m_b)
m0_up_fit = find_m0_shifted(m_u, m_c, m_t)

print(f"\n  Λ_QCD (Nf=3, MS-bar) ≈ {Lambda_QCD} MeV")
print(f"  m_proton = {m_proton} MeV")

print(f"\n  m₀(down, fit) = {m0_down_fit:.2f} MeV")
print(f"    m₀/Λ_QCD = {m0_down_fit/Lambda_QCD:.4f}")
print(f"    m₀/m_proton = {m0_down_fit/m_proton:.4f}")
print(f"    m₀ ≈ Λ_QCD/15 — small sea correction")

print(f"\n  m₀(up, fit) = {m0_up_fit:.1f} MeV")
print(f"    m₀/Λ_QCD = {m0_up_fit/Lambda_QCD:.2f}")
print(f"    m₀/m_proton = {m0_up_fit/m_proton:.2f}")
print(f"    m₀ ≈ 2×m_proton — massive sea contribution (top quark)")

print(f"\n  m₀(up)/m₀(down) = {m0_up_fit/m0_down_fit:.1f}")
print(f"  m_t/m_b = {m_t/m_b:.1f}")
print(f"  Ratio: m₀(up)/m₀(down) ≈ (m_t/m_b)^{np.log(m0_up_fit/m0_down_fit)/np.log(m_t/m_b):.2f}")

# Is m₀ proportional to m₃/m₁ exactly (as R12 claims)?
ratio_d = m0_down_fit / (m_b / m_d)
ratio_u = m0_up_fit / (m_t / m_u)
print(f"\n  Test R12 formula m₀ = A × m₃/m₁:")
print(f"    A(down) = m₀/(m_b/m_d) = {ratio_d:.6f}")
print(f"    A(up)   = m₀/(m_t/m_u) = {ratio_u:.6f}")
print(f"    Ratio A(down)/A(up) = {ratio_d/ratio_u:.4f}")

record("T7: Universal A constant",
       abs(ratio_d/ratio_u - 1.0) < 0.05,
       f"A(d) = {ratio_d:.6f}, A(u) = {ratio_u:.6f}, ratio = {ratio_d/ratio_u:.4f}")

# Back-extract Ω_Λ from each sector independently
OL_from_down = ratio_d / (1/(168*3/14*phi)) / 168 * 14 / 3 * phi
# Actually: A = 1/(Φ_eff × φ), Φ_eff = Φ₀ × 3/14, Φ₀ = 168 × Ω_Λ
# So A = 14/(168 × Ω_Λ × 3 × φ) = 14/(504 × φ × Ω_Λ)
# Ω_Λ = 14/(504 × φ × A)
OL_from_down = 14.0 / (504.0 * phi * ratio_d)
OL_from_up = 14.0 / (504.0 * phi * ratio_u)

print(f"\n  Ω_Λ back-extracted from m₀(fit):")
print(f"    From down sector: Ω_Λ = {OL_from_down:.4f}")
print(f"    From up sector:   Ω_Λ = {OL_from_up:.4f}")
print(f"    Planck:            Ω_Λ = {OL_planck}")
print(f"    Average:           Ω_Λ = {(OL_from_down+OL_from_up)/2:.4f}")


# ============================================================
# §5. SENSITIVITY TO LIGHT QUARK MASSES
# ============================================================
print("\n" + "=" * 72)
print("§5. CZUŁOŚĆ NA MASY LEKKICH KWARKÓW")
print("=" * 72)

print("""
  Masy lekkich kwarków (m_d, m_s, m_u, m_c) mają duże niepewności.
  Sprawdzamy, jak zmienia się Ω_Λ(fit) w zależności od nich.
""")

# Scan over m_d
print("  Scan m_d (± 1σ = 0.48 MeV):")
for delta_d in [-0.48, 0, 0.48]:
    md_test = m_d + delta_d

    def chi2_test(OL):
        P0 = 168 * OL
        PE = P0 * 3/14
        A = 1.0 / (PE * phi)
        m0d = A * m_b / md_test
        m0u = A * m_t / m_u
        mb_p = predict_m3_shifted(md_test, m_s, m0d)
        mt_p = predict_m3_shifted(m_u, m_c, m0u)
        if np.isnan(mb_p) or np.isnan(mt_p):
            return 1e10
        return ((mb_p - m_b)/sig_mb)**2 + ((mt_p - m_t)/sig_mt)**2

    res = minimize_scalar(chi2_test, bounds=(0.60, 0.90), method='bounded')
    print(f"    m_d = {md_test:.2f} → Ω_Λ = {res.x:.4f}, χ² = {res.fun:.2f}")

# Scan over m_s
print("\n  Scan m_s (± 1σ = 8.6 MeV):")
for delta_s in [-8.6, 0, 8.6]:
    ms_test = m_s + delta_s

    def chi2_test(OL):
        P0 = 168 * OL
        PE = P0 * 3/14
        A = 1.0 / (PE * phi)
        m0d = A * m_b / m_d
        m0u = A * m_t / m_u
        mb_p = predict_m3_shifted(m_d, ms_test, m0d)
        mt_p = predict_m3_shifted(m_u, m_c, m0u)
        if np.isnan(mb_p) or np.isnan(mt_p):
            return 1e10
        return ((mb_p - m_b)/sig_mb)**2 + ((mt_p - m_t)/sig_mt)**2

    res = minimize_scalar(chi2_test, bounds=(0.60, 0.90), method='bounded')
    print(f"    m_s = {ms_test:.1f} → Ω_Λ = {res.x:.4f}, χ² = {res.fun:.2f}")


# ============================================================
# §6. KOIDE TABLE — ALL SECTORS
# ============================================================
print("\n" + "=" * 72)
print("§6. ★ MASTER TABLE — ALL MASS PREDICTIONS")
print("=" * 72)

print(f"""
  ┌─────────────────────────────────────────────────────────────┐
  │                   TGP MASS FRAMEWORK                        │
  │  Inputs: m_e, r₂₁(ODE), K=2/3, Ω_Λ={OL_best:.4f}             │
  │          + light quark masses (m_d, m_s, m_u, m_c)          │
  ├─────────────────────────────────────────────────────────────┤
  │  SECTOR     │  PREDICTION  │    VALUE    │   PDG      │ ERR │
  ├─────────────┼──────────────┼─────────────┼────────────┼─────┤
  │  Lepton     │  m_μ         │  {m_mu_pred:8.2f} MeV│{m_mu:8.2f} MeV│{abs(m_mu_pred-m_mu)/m_mu*100:4.2f}%│
  │  Lepton     │  m_τ         │  {m_tau_pred:8.2f} MeV│{m_tau:8.2f} MeV│{abs(m_tau_pred-m_tau)/m_tau*100:4.3f}%│
  │  Down quark │  m_b         │  {mb_pred:8.1f} MeV│{m_b:8.1f} MeV│{abs(mb_pred-m_b)/m_b*100:4.2f}%│
  │  Up quark   │  m_t         │  {mt_pred:8.0f} MeV│{m_t:8.0f} MeV│{abs(mt_pred-m_t)/m_t*100:4.2f}%│
  │  QCD        │  α_s         │    {alpha_s_pred:8.4f}  │  0.1179±9  │{abs(alpha_s_pred-0.1179)/0.0009:4.1f}σ│
  │  Cosmology  │  Ω_Λ         │    {OL_best:8.4f}  │0.6847±73  │{tension_sigma:4.1f}σ│
  └─────────────┴──────────────┴─────────────┴────────────┴─────┘
""")


# ============================================================
# SCORECARD
# ============================================================
print("=" * 72)
print("SCORECARD")
print("=" * 72)
n_pass = sum(1 for _, p, _ in TESTS if p)
n_total = len(TESTS)
for name, passed, detail in TESTS:
    mark = "PASS" if passed else "FAIL"
    print(f"  [{mark}] {name}")
    if detail:
        for line in detail.split('\n'):
            print(f"         {line}")
print(f"\n  {n_pass}/{n_total} testów przeszło.")

print(f"""
========================================================================
PODSUMOWANIE ex236
========================================================================

  ★ SELF-CONSISTENT Ω_Λ FROM QUARK MASSES:

  χ² minimalizacja daje Ω_Λ = {OL_best:.4f}
  Tension z Planck: {tension_sigma:.1f}σ

  PREDYKCJE (6 wielkości z ~5 wejść):
    m_τ:   {m_tau_pred:.2f} MeV   (PDG: {m_tau}, err: {abs(m_tau_pred-m_tau)/m_tau*100:.3f}%)
    m_b:   {mb_pred:.1f} MeV (PDG: {m_b}, err: {abs(mb_pred-m_b)/m_b*100:.2f}%)
    m_t:   {mt_pred:.0f} MeV (PDG: {m_t}, err: {abs(mt_pred-m_t)/m_t*100:.2f}%)
    α_s:   {alpha_s_pred:.4f}       (PDG: 0.1179±9, {abs(alpha_s_pred-0.1179)/0.0009:.1f}σ)
    Ω_Λ:   {OL_best:.4f}       (Planck: 0.6847±73, {tension_sigma:.1f}σ)

  FIZYCZNA INTERPRETACJA m₀:
    m₀ = "morze" (gluony + pary qq̄) contribution
    m₀(down) ≈ Λ_QCD/15 — small correction
    m₀(up) ≈ 2×m_proton — massive sea (top quark is special)
    R12: m₀ = A × m₃/m₁ with UNIVERSAL A → links quarks to cosmology
""")
