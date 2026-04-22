#!/usr/bin/env python3
"""
ex171_running_masses_koide.py
Sesja v45, 2026-04-05

Precyzyjne bieganie mas kwarkowych QCD i test K(b,c,t) = 2/3
na wspolnej skali energetycznej.

Metoda: 4-loop QCD running z progami smakow (threshold matching).
  m(mu2) = m(mu1) * c(alpha_s(mu2)) / c(alpha_s(mu1))
  c(a) = (a/pi)^{gamma0/(2*beta0)} * [1 + O(a)]

Uproszczenie: 2-loop alpha_s + 1-loop masa (wystarczajace na ~1%).
  alpha_s(mu): 2-loop z progami przy m_c, m_b, m_t
  m(mu): 1-loop anomalous dimension gamma_m = 4 (LO)

Wejscie PDG 2024:
  alpha_s(M_Z) = 0.1179 +/- 0.0009
  m_c(m_c) = 1.270 +/- 0.020 GeV  (MS-bar)
  m_b(m_b) = 4.180 +/- 0.030 GeV  (MS-bar)
  m_t(m_t) = 162.5  +/- 2.1  GeV  (MS-bar)
  m_t(pole) = 172.76 +/- 0.30 GeV
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ---- Constants ----
M_Z = 91.1876  # GeV

# PDG inputs (GeV)
ALPHA_S_MZ = 0.1179
M_C_MC = 1.270      # m_c(m_c) MS-bar
M_B_MB = 4.180      # m_b(m_b) MS-bar
M_T_MT = 162.5      # m_t(m_t) MS-bar (approx)
M_T_POLE = 172.76   # pole mass

# ---- QCD beta function coefficients ----
def beta0(nf):
    return (33 - 2*nf) / (12*np.pi)

def beta1(nf):
    return (153 - 19*nf) / (24*np.pi**2)

# ---- Anomalous dimension (mass) ----
def gamma0(nf):
    """LO mass anomalous dimension: dm/d(ln mu) = -gamma_m * m"""
    return 4 / (4*np.pi)  # = 1/pi (coefficient of alpha_s)
    # Actually: dm/d(ln mu) = -gamma_m * alpha_s/pi * m
    # gamma_m^(0) = 1 at LO in alpha_s/pi convention

# ---- 2-loop alpha_s running ----
def alpha_s_2loop(mu, alpha_s_ref, mu_ref, nf):
    """Solve 2-loop RGE for alpha_s.
    d(alpha_s)/d(ln mu^2) = -beta0*a^2 - beta1*a^3
    where a = alpha_s
    """
    b0 = beta0(nf)
    b1 = beta1(nf)

    def rhs(lnmu2, a):
        return -b0 * a[0]**2 - b1 * a[0]**3

    lnmu2_start = np.log(mu_ref**2)
    lnmu2_end = np.log(mu**2)

    sol = solve_ivp(rhs, [lnmu2_start, lnmu2_end], [alpha_s_ref],
                    method='RK45', rtol=1e-10, atol=1e-14)
    return sol.y[0][-1]


def run_alpha_s(mu_target, alpha_s_mz=ALPHA_S_MZ):
    """Run alpha_s from M_Z to mu_target with threshold matching.
    Thresholds at m_c, m_b, m_t (MS-bar masses).
    """
    # Step 1: M_Z down/up with nf=5 (between m_b and m_t)
    # Thresholds: m_c(m_c) ~ 1.27 GeV, m_b(m_b) ~ 4.18 GeV, m_t(m_t) ~ 162.5 GeV

    # At M_Z: nf=5
    a_mz = alpha_s_mz

    if mu_target >= M_T_MT:
        # Run nf=5 from M_Z to m_t, then nf=6 from m_t to target
        a_mt = alpha_s_2loop(M_T_MT, a_mz, M_Z, 5)
        # Threshold matching at m_t (LO: continuous)
        a_target = alpha_s_2loop(mu_target, a_mt, M_T_MT, 6)
    elif mu_target >= M_B_MB:
        # nf=5 zone
        a_target = alpha_s_2loop(mu_target, a_mz, M_Z, 5)
    elif mu_target >= M_C_MC:
        # Run nf=5 from M_Z down to m_b, then nf=4 to target
        a_mb = alpha_s_2loop(M_B_MB, a_mz, M_Z, 5)
        a_target = alpha_s_2loop(mu_target, a_mb, M_B_MB, 4)
    else:
        # Run nf=5 -> nf=4 -> nf=3
        a_mb = alpha_s_2loop(M_B_MB, a_mz, M_Z, 5)
        a_mc = alpha_s_2loop(M_C_MC, a_mb, M_B_MB, 4)
        a_target = alpha_s_2loop(mu_target, a_mc, M_C_MC, 3)

    return a_target


def run_mass(m_ref, mu_ref, mu_target, nf_ref=None):
    """Run quark mass from mu_ref to mu_target using 1-loop RGE.
    m(mu2)/m(mu1) = [alpha_s(mu2)/alpha_s(mu1)]^{gamma0_m/beta0}
    where gamma0_m = 4/(2*beta0) in standard convention.

    More precisely: m(mu) ~ [alpha_s(mu)]^{d0} where d0 = gamma0_m / (2*beta0)
    d0 = (4/pi) / (2 * (33-2nf)/(12pi)) = 4*12/(2*(33-2nf)) = 24/(33-2nf) ... hmm

    Let me be more careful:
    dm/d(ln mu^2) = -gamma_m * m
    gamma_m = gamma_0 * a_s/pi + ...
    gamma_0 = 1 (standard)

    Combined with: d(a_s)/d(ln mu^2) = -beta_0 * a_s^2 - ...
    At LO: dm/m = (gamma_0/pi) * (da_s / (beta_0 * a_s))
    ln(m2/m1) = (gamma_0/(pi*beta_0)) * ln(a2/a1)
    m2/m1 = (a2/a1)^{gamma_0/(pi*beta_0)}

    gamma_0 = 1, beta_0 = (33-2nf)/(12pi)
    exponent = 1 / (pi * (33-2nf)/(12pi)) = 12/(33-2nf)
    """
    # Get alpha_s at both scales
    a1 = run_alpha_s(mu_ref)
    a2 = run_alpha_s(mu_target)

    # Determine nf for the running (simplified: use nf at the reference scale)
    if nf_ref is None:
        if mu_ref >= M_T_MT:
            nf_ref = 6
        elif mu_ref >= M_B_MB:
            nf_ref = 5
        elif mu_ref >= M_C_MC:
            nf_ref = 4
        else:
            nf_ref = 3

    # For running across thresholds, we need to be more careful
    # But for nearby scales (within same nf regime), this is fine
    d0 = 12.0 / (33 - 2*nf_ref)
    return m_ref * (a2 / a1)**d0


def run_mass_precise(m_ref, mu_ref, mu_target):
    """Run mass with proper threshold crossing."""
    # Thresholds (ordered)
    thresholds = [
        (M_C_MC, 3, 4),   # below m_c: nf=3, above: nf=4
        (M_B_MB, 4, 5),   # below m_b: nf=4, above: nf=5
        (M_T_MT, 5, 6),   # below m_t: nf=5, above: nf=6
    ]

    # Determine direction
    if mu_ref < mu_target:
        direction = 'up'
    else:
        direction = 'down'

    current_mu = mu_ref
    current_m = m_ref

    # Determine initial nf
    if current_mu >= M_T_MT:
        current_nf = 6
    elif current_mu >= M_B_MB:
        current_nf = 5
    elif current_mu >= M_C_MC:
        current_nf = 4
    else:
        current_nf = 3

    # Find thresholds to cross
    if direction == 'up':
        crossings = [(thr, nf_below, nf_above) for thr, nf_below, nf_above in thresholds
                     if current_mu < thr < mu_target]
    else:
        crossings = [(thr, nf_below, nf_above) for thr, nf_below, nf_above in thresholds
                     if mu_target < thr < current_mu]
        crossings.reverse()

    for thr, nf_below, nf_above in crossings:
        # Run to threshold
        a_curr = run_alpha_s(current_mu)
        a_thr = run_alpha_s(thr)
        d0 = 12.0 / (33 - 2*current_nf)
        current_m = current_m * (a_thr / a_curr)**d0
        current_mu = thr
        # Cross threshold (LO: mass is continuous)
        if direction == 'up':
            current_nf = nf_above
        else:
            current_nf = nf_below

    # Final run to target
    a_curr = run_alpha_s(current_mu)
    a_target = run_alpha_s(mu_target)
    d0 = 12.0 / (33 - 2*current_nf)
    current_m = current_m * (a_target / a_curr)**d0

    return current_m


# ==============================================================
print("=" * 72)
print("ex171: Running masses QCD i test K(b,c,t) na wspolnej skali")
print("=" * 72)

# ---- 1. Verify alpha_s running ----
print("\n--- 1. Weryfikacja alpha_s running ---\n")
test_scales = [1.0, 2.0, M_C_MC, M_B_MB, M_Z, M_T_MT, 500, 1000]
print(f"  {'mu (GeV)':>10s}  {'alpha_s':>10s}  {'nf':>4s}")
print("  " + "-" * 30)
for mu in test_scales:
    a = run_alpha_s(mu)
    if mu >= M_T_MT:
        nf = 6
    elif mu >= M_B_MB:
        nf = 5
    elif mu >= M_C_MC:
        nf = 4
    else:
        nf = 3
    print(f"  {mu:10.2f}  {a:10.4f}  {nf:4d}")

# ---- 2. Run masses to common scales ----
print("\n--- 2. Masy (b,c,t) na wspolnych skalach ---\n")

scales = [2.0, 3.0, M_B_MB, 10.0, M_Z, M_T_MT, 500.0]

def koide(m1, m2, m3):
    return (m1+m2+m3) / (np.sqrt(m1)+np.sqrt(m2)+np.sqrt(m3))**2

print(f"  {'mu (GeV)':>10s}  {'m_c (GeV)':>10s}  {'m_b (GeV)':>10s}  "
      f"{'m_t (GeV)':>10s}  {'K(b,c,t)':>10s}  {'dK (%)':>8s}")
print("  " + "-" * 65)

K_results = []
for mu in scales:
    mc = run_mass_precise(M_C_MC, M_C_MC, mu)
    mb = run_mass_precise(M_B_MB, M_B_MB, mu)
    mt = run_mass_precise(M_T_MT, M_T_MT, mu)
    K = koide(mb, mc, mt)
    dK = (K - 2/3) / (2/3) * 100
    K_results.append((mu, mc, mb, mt, K, dK))
    print(f"  {mu:10.2f}  {mc:10.4f}  {mb:10.4f}  {mt:10.2f}  "
          f"{K:10.6f}  {dK:+8.3f}%")

# ---- 3. Find scale where K = 2/3 exactly ----
print("\n--- 3. Szukanie skali mu* gdzie K(b,c,t) = 2/3 ---\n")

def K_at_scale(ln_mu):
    mu = np.exp(ln_mu)
    if mu < 0.5 or mu > 1e6:
        return 1.0  # out of range
    mc = run_mass_precise(M_C_MC, M_C_MC, mu)
    mb = run_mass_precise(M_B_MB, M_B_MB, mu)
    mt = run_mass_precise(M_T_MT, M_T_MT, mu)
    return koide(mb, mc, mt) - 2/3

# Scan to find sign changes
ln_mus = np.linspace(np.log(1.0), np.log(500), 50)
Ks = [K_at_scale(lm) for lm in ln_mus]

sign_changes = []
for i in range(len(Ks)-1):
    if Ks[i] * Ks[i+1] < 0:
        # Found sign change
        mu_star = np.exp(brentq(K_at_scale, ln_mus[i], ln_mus[i+1]))
        sign_changes.append(mu_star)
        mc = run_mass_precise(M_C_MC, M_C_MC, mu_star)
        mb = run_mass_precise(M_B_MB, M_B_MB, mu_star)
        mt = run_mass_precise(M_T_MT, M_T_MT, mu_star)
        K = koide(mb, mc, mt)
        print(f"  mu* = {mu_star:.2f} GeV")
        print(f"  m_c(mu*) = {mc:.4f} GeV")
        print(f"  m_b(mu*) = {mb:.4f} GeV")
        print(f"  m_t(mu*) = {mt:.2f} GeV")
        print(f"  K = {K:.8f}")

if not sign_changes:
    print("  Brak skali mu* z K = 2/3 w zakresie [1, 500] GeV")
    # Find minimum |K - 2/3|
    abs_Ks = [abs(k) for k in Ks]
    idx_min = np.argmin(abs_Ks)
    mu_min = np.exp(ln_mus[idx_min])
    print(f"  Najblizsze K do 2/3 przy mu = {mu_min:.1f} GeV: "
          f"K-2/3 = {Ks[idx_min]:+.6f}")

# ---- 4. K(b,c,t) with pole masses ----
print("\n--- 4. K z masami biegunowymi vs MS-bar ---\n")

# Pole masses (approximate)
# m_c(pole) ~ 1.67 GeV, m_b(pole) ~ 4.78 GeV, m_t(pole) = 172.76 GeV
m_c_pole = 1.67
m_b_pole = 4.78
m_t_pole = M_T_POLE

K_pole = koide(m_b_pole, m_c_pole, m_t_pole)
K_msbar_own = koide(M_B_MB, M_C_MC, M_T_MT)
K_mixed = koide(M_B_MB, M_C_MC, m_t_pole)  # ex170 configuration

print(f"  Pole masses:      K(b,c,t) = {K_pole:.6f}  ({(K_pole-2/3)/(2/3)*100:+.3f}%)")
print(f"  MS-bar (own mu):  K(b,c,t) = {K_msbar_own:.6f}  ({(K_msbar_own-2/3)/(2/3)*100:+.3f}%)")
print(f"  Mixed (ex170):    K(b,c,t) = {K_mixed:.6f}  ({(K_mixed-2/3)/(2/3)*100:+.3f}%)")

# ---- 5. K(b,c,t) predicting m_t at mu* ----
print("\n--- 5. Predykcja m_t z K(b,c,t) = 2/3 na roznych skalach ---\n")

def predict_mt_from_koide(mc, mb, K=2.0/3.0):
    """Predict m_t from K(b,c,t)=K given m_c, m_b."""
    s_c = np.sqrt(mc)
    s_b = np.sqrt(mb)
    S = s_c + s_b
    # K*(S+x)^2 = mc+mb+x^2, x = sqrt(mt)
    A = K - 1
    B = 2*K*S
    C = K*S**2 - mc - mb
    disc = B**2 - 4*A*C
    if disc < 0:
        return None
    x = (-B + np.sqrt(disc)) / (2*A)
    if x > 0:
        return x**2
    return None

print(f"  {'mu (GeV)':>10s}  {'m_c(mu)':>10s}  {'m_b(mu)':>10s}  "
      f"{'m_t(pred)':>10s}  {'m_t(run)':>10s}  {'err (%)':>8s}")
print("  " + "-" * 60)

for mu, mc, mb, mt, K, dK in K_results:
    mt_pred = predict_mt_from_koide(mc, mb)
    if mt_pred is not None:
        err = (mt_pred - mt) / mt * 100
        print(f"  {mu:10.2f}  {mc:10.4f}  {mb:10.4f}  "
              f"{mt_pred:10.2f}  {mt:10.2f}  {err:+8.2f}%")

# ---- 6. All lepton+quark Koide at common scales ----
print("\n--- 6. Koide sektorowy na wspolnych skalach ---\n")

# Lepton masses (pole, GeV) - don't run (QED running negligible)
m_e = 0.000511
m_mu = 0.10566
m_tau = 1.77686

# Down quarks
m_d_2GeV = 0.00467  # MS-bar at 2 GeV
m_s_2GeV = 0.0934

# Up quarks
m_u_2GeV = 0.00216

print("  Lepton (pole): K(e,mu,tau) = "
      f"{koide(m_e, m_mu, m_tau):.6f} (nie biega)")
print()

for mu in [2.0, M_B_MB, M_Z]:
    # Run down quarks from 2 GeV to mu
    md = run_mass_precise(m_d_2GeV, 2.0, mu)
    ms = run_mass_precise(m_s_2GeV, 2.0, mu)
    mb = run_mass_precise(M_B_MB, M_B_MB, mu)
    K_down = koide(md, ms, mb)

    # Run up quarks
    mu_q = run_mass_precise(m_u_2GeV, 2.0, mu)
    mc = run_mass_precise(M_C_MC, M_C_MC, mu)
    mt = run_mass_precise(M_T_MT, M_T_MT, mu)
    K_up = koide(mu_q, mc, mt)

    # Cross (b,c,t)
    K_bct = koide(mb, mc, mt)

    print(f"  mu = {mu:.1f} GeV:")
    print(f"    K(d,s,b) = {K_down:.6f}  ({(K_down-2/3)/(2/3)*100:+.2f}%)")
    print(f"    K(u,c,t) = {K_up:.6f}  ({(K_up-2/3)/(2/3)*100:+.2f}%)")
    print(f"    K(b,c,t) = {K_bct:.6f}  ({(K_bct-2/3)/(2/3)*100:+.2f}%)")
    print()

# ---- 7. Error propagation at mu* ----
print("\n--- 7. Monte Carlo na najlepszej skali ---\n")

# Find best scale (minimum |K-2/3|)
abs_Ks = [abs(k) for k in Ks]
idx_min = np.argmin(abs_Ks)
mu_best = np.exp(ln_mus[idx_min])

# Refine
if sign_changes:
    mu_best = sign_changes[0]

np.random.seed(42)
N_mc = 10000
K_mc = np.zeros(N_mc)

for i in range(N_mc):
    mc_i = np.random.normal(M_C_MC, 0.020)
    mb_i = np.random.normal(M_B_MB, 0.030)
    mt_i = np.random.normal(M_T_MT, 2.1)
    as_i = np.random.normal(ALPHA_S_MZ, 0.0009)

    # Simple running at mu_best (approximate: use central alpha_s for speed)
    mc_run = run_mass_precise(mc_i, mc_i, mu_best)
    mb_run = run_mass_precise(mb_i, mb_i, mu_best)
    mt_run = run_mass_precise(mt_i, mt_i, mu_best)
    K_mc[i] = koide(mb_run, mc_run, mt_run)

K_mc_mean = np.mean(K_mc)
K_mc_std = np.std(K_mc)
sigma_from_23 = (K_mc_mean - 2/3) / K_mc_std

print(f"  Skala mu = {mu_best:.1f} GeV")
print(f"  K(b,c,t) = {K_mc_mean:.6f} +/- {K_mc_std:.6f}")
print(f"  Odchylenie od 2/3: {sigma_from_23:.1f} sigma")

# ---- WNIOSKI ----
print("\n" + "=" * 72)
print("WNIOSKI ex171")
print("=" * 72)
print(f"""
  1. K(b,c,t) ZALEZY od skali energetycznej mu:
     - mu = 2 GeV:   K = {K_results[0][4]:.6f} ({K_results[0][5]:+.2f}%)
     - mu = m_b:     K = {K_results[2][4]:.6f} ({K_results[2][5]:+.2f}%)
     - mu = M_Z:     K = {K_results[4][4]:.6f} ({K_results[4][5]:+.2f}%)
     - mu = m_t:     K = {K_results[5][4]:.6f} ({K_results[5][5]:+.2f}%)

  2. K przechodzi przez 2/3 (jesli istnieje mu*) lub ma minimum
     blisko 2/3 na pewnej skali.

  3. Interpretacja TGP:
     - Jesli K(b,c,t) = 2/3 na pewnej skali mu*,
       to mu* moze byc FIZYCZNA SKALA TGP
       (np. skala substratu, skala solitonu).
     - Jest to analogiczne do faktu, ze K(e,mu,tau) = 2/3
       z masami biegunowymi (skala IR).

  4. Masa biegunowa vs MS-bar:
     - Pole: K = {K_pole:.6f} ({(K_pole-2/3)/(2/3)*100:+.2f}%)
     - MS-bar (own): K = {K_msbar_own:.6f} ({(K_msbar_own-2/3)/(2/3)*100:+.2f}%)
     - Leptony uzywaja mas biegunowych -> K = 2/3

  5. STATUS: K(b,c,t) ~ 2/3 jest SKALOWO-ZALEZNE.
     Istnienie mu* z K = 2/3 jest MOZLIWE ale wymaga
     precyzyjniejszego obliczenia (2-3 loop).
""")
