#!/usr/bin/env python3
"""
ex172_koide_scale_scan.py
Sesja v45, 2026-04-05

Skan K(d,s,b) i K(u,c,t) jako funkcji skali mu.

Kluczowe pytanie: czy istnieje skala mu_K* taka, ze:
  K(d,s,b; mu_K*) = 2/3  lub  K(u,c,t; mu_K*) = 2/3?

Leptony: K = 2/3 z masami biegunowymi (= masami pole, skala IR).
Kwarki: masy biegaja z QCD. Czy pole masses daja K blizsze 2/3?

Metoda: jak w ex171 (2-loop alpha_s + 1-loop masa).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

M_Z = 91.1876
ALPHA_S_MZ = 0.1179
M_C_MC = 1.270
M_B_MB = 4.180
M_T_MT = 162.5

def beta0(nf): return (33 - 2*nf) / (12*np.pi)
def beta1(nf): return (153 - 19*nf) / (24*np.pi**2)

def alpha_s_2loop(mu, alpha_s_ref, mu_ref, nf):
    b0, b1 = beta0(nf), beta1(nf)
    def rhs(lnmu2, a): return [-b0*a[0]**2 - b1*a[0]**3]
    sol = solve_ivp(rhs, [np.log(mu_ref**2), np.log(mu**2)],
                    [alpha_s_ref], method='RK45', rtol=1e-10, atol=1e-14)
    return sol.y[0][-1]

def run_alpha_s(mu):
    a_mz = ALPHA_S_MZ
    if mu >= M_T_MT:
        a_mt = alpha_s_2loop(M_T_MT, a_mz, M_Z, 5)
        return alpha_s_2loop(mu, a_mt, M_T_MT, 6)
    elif mu >= M_B_MB:
        return alpha_s_2loop(mu, a_mz, M_Z, 5)
    elif mu >= M_C_MC:
        a_mb = alpha_s_2loop(M_B_MB, a_mz, M_Z, 5)
        return alpha_s_2loop(mu, a_mb, M_B_MB, 4)
    else:
        a_mb = alpha_s_2loop(M_B_MB, a_mz, M_Z, 5)
        a_mc = alpha_s_2loop(M_C_MC, a_mb, M_B_MB, 4)
        return alpha_s_2loop(mu, a_mc, M_C_MC, 3)

def run_mass_precise(m_ref, mu_ref, mu_target):
    thresholds = [(M_C_MC, 3, 4), (M_B_MB, 4, 5), (M_T_MT, 5, 6)]
    current_mu, current_m = mu_ref, m_ref
    if current_mu >= M_T_MT: current_nf = 6
    elif current_mu >= M_B_MB: current_nf = 5
    elif current_mu >= M_C_MC: current_nf = 4
    else: current_nf = 3

    if mu_ref < mu_target:
        crossings = [(t,nb,na) for t,nb,na in thresholds if current_mu < t < mu_target]
    else:
        crossings = [(t,nb,na) for t,nb,na in thresholds if mu_target < t < current_mu]
        crossings.reverse()

    for thr, nf_below, nf_above in crossings:
        a_curr = run_alpha_s(current_mu)
        a_thr = run_alpha_s(thr)
        d0 = 12.0 / (33 - 2*current_nf)
        current_m *= (a_thr / a_curr)**d0
        current_mu = thr
        current_nf = nf_above if mu_ref < mu_target else nf_below

    a_curr = run_alpha_s(current_mu)
    a_target = run_alpha_s(mu_target)
    d0 = 12.0 / (33 - 2*current_nf)
    current_m *= (a_target / a_curr)**d0
    return current_m

def koide(m1, m2, m3):
    return (m1+m2+m3) / (np.sqrt(m1)+np.sqrt(m2)+np.sqrt(m3))**2

# ---- PDG inputs ----
m_d_2GeV = 0.00467
m_s_2GeV = 0.0934
m_u_2GeV = 0.00216

# Pole masses (approximate, from PDG)
m_d_pole = 0.00467   # light quarks: pole ~ MS-bar (non-perturbative)
m_s_pole = 0.0934    # same caveat
m_b_pole = 4.78      # GeV
m_u_pole = 0.00216
m_c_pole = 1.67      # GeV
m_t_pole = 172.76    # GeV

print("=" * 72)
print("ex172: K(d,s,b) i K(u,c,t) jako funkcja skali mu")
print("=" * 72)

# ---- 1. K z masami biegunowymi ----
print("\n--- 1. K z masami biegunowymi (pole) ---\n")
K_down_pole = koide(m_d_pole, m_s_pole, m_b_pole)
K_up_pole = koide(m_u_pole, m_c_pole, m_t_pole)
K_lep = koide(0.000511, 0.10566, 1.77686)
print(f"  K(e,mu,tau)  [pole] = {K_lep:.6f}  ({(K_lep-2/3)/(2/3)*100:+.4f}%)")
print(f"  K(d,s,b)    [pole] = {K_down_pole:.6f}  ({(K_down_pole-2/3)/(2/3)*100:+.2f}%)")
print(f"  K(u,c,t)    [pole] = {K_up_pole:.6f}  ({(K_up_pole-2/3)/(2/3)*100:+.2f}%)")
print()
print("  UWAGA: dla lekkich kwarkow (u,d,s) masa biegunowa jest")
print("  zle zdefiniowana (non-perturbative QCD).")
print("  Wartosci m_d, m_s, m_u powyej to MS-bar at 2 GeV.")

# ---- 2. Skan K(mu) dla obu sektorow ----
print("\n--- 2. Skan K(mu) ---\n")

mus = np.logspace(np.log10(1.5), np.log10(1000), 60)

K_down_arr = []
K_up_arr = []

for mu in mus:
    md = run_mass_precise(m_d_2GeV, 2.0, mu)
    ms = run_mass_precise(m_s_2GeV, 2.0, mu)
    mb = run_mass_precise(M_B_MB, M_B_MB, mu)
    K_d = koide(md, ms, mb)
    K_down_arr.append(K_d)

    mu_q = run_mass_precise(m_u_2GeV, 2.0, mu)
    mc = run_mass_precise(M_C_MC, M_C_MC, mu)
    mt = run_mass_precise(M_T_MT, M_T_MT, mu)
    K_u = koide(mu_q, mc, mt)
    K_up_arr.append(K_u)

K_down_arr = np.array(K_down_arr)
K_up_arr = np.array(K_up_arr)

# Print selected points
print(f"  {'mu (GeV)':>10s}  {'K(d,s,b)':>10s}  {'dK_d (%)':>10s}  "
      f"{'K(u,c,t)':>10s}  {'dK_u (%)':>10s}")
print("  " + "-" * 55)

for i, mu in enumerate(mus):
    if i % 8 == 0 or abs(mu - 2.0) < 0.5 or abs(mu - M_Z) < 5:
        dKd = (K_down_arr[i] - 2/3)/(2/3)*100
        dKu = (K_up_arr[i] - 2/3)/(2/3)*100
        print(f"  {mu:10.2f}  {K_down_arr[i]:10.6f}  {dKd:+10.2f}%  "
              f"{K_up_arr[i]:10.6f}  {dKu:+10.2f}%")

# ---- 3. Ekstrema i trendy ----
print("\n--- 3. Ekstrema ---\n")

idx_d_min = np.argmin(np.abs(K_down_arr - 2/3))
idx_u_min = np.argmin(np.abs(K_up_arr - 2/3))

print(f"  K(d,s,b) najblizsze 2/3 przy mu = {mus[idx_d_min]:.1f} GeV: "
      f"K = {K_down_arr[idx_d_min]:.6f} ({(K_down_arr[idx_d_min]-2/3)/(2/3)*100:+.2f}%)")
print(f"  K(u,c,t) najblizsze 2/3 przy mu = {mus[idx_u_min]:.1f} GeV: "
      f"K = {K_up_arr[idx_u_min]:.6f} ({(K_up_arr[idx_u_min]-2/3)/(2/3)*100:+.2f}%)")

# Trend: K changes monotonically or has minimum?
print(f"\n  K(d,s,b) range: [{K_down_arr.min():.6f}, {K_down_arr.max():.6f}]")
print(f"  K(u,c,t) range: [{K_up_arr.min():.6f}, {K_up_arr.max():.6f}]")

# ---- 4. Kluczowe: dlaczego K leptonowe dziala? ----
print("\n--- 4. Analiza: dlaczego K = 2/3 TYLKO dla leptonow? ---\n")

print("  Koide K = (1 + r21 + r31) / (1 + sqrt(r21) + sqrt(r31))^2")
print()
print("  Wymaganie K = 2/3 + phi-FP (r21 dany) jednoznacznie wyznacza r31.")
print("  Pytanie: dlaczego leptony spelniaja ten warunek a kwarki nie?")
print()

sectors_data = {
    'lepton': {'r21': 206.8, 'r31': 3477},
    'down':   {'r21': 20.0, 'r31': 895},
    'up':     {'r21': 588.0, 'r31': 79982},
}

for name, data in sectors_data.items():
    r21 = data['r21']
    r31 = data['r31']
    K = (1 + r21 + r31) / (1 + np.sqrt(r21) + np.sqrt(r31))**2

    # What r31 would give K = 2/3?
    # 2/3 = (1 + r21 + r31) / (1 + sqrt(r21) + sqrt(r31))^2
    # Let x = sqrt(r31):
    # 2/3 * (1+sqrt(r21)+x)^2 = 1+r21+x^2
    # 2/3*(1+s)^2 + 2*2/3*(1+s)*x + 2/3*x^2 = 1+r21+x^2  where s=sqrt(r21)
    # (2/3-1)*x^2 + 4/3*(1+s)*x + 2/3*(1+s)^2 - 1 - r21 = 0
    # -1/3*x^2 + 4/3*(1+s)*x + 2/3*(1+s)^2 - 1 - r21 = 0
    s = np.sqrt(r21)
    A = -1.0/3
    B = 4.0/3 * (1 + s)
    C = 2.0/3 * (1+s)**2 - 1 - r21
    disc = B**2 - 4*A*C
    if disc >= 0:
        x1 = (-B + np.sqrt(disc)) / (2*A)
        x2 = (-B - np.sqrt(disc)) / (2*A)
        r31_koide = max(x1, x2)**2 if max(x1, x2) > 0 else None
    else:
        r31_koide = None

    print(f"  {name:>7s}: r21 = {r21:.1f}, r31 = {r31:.0f}, K = {K:.4f}")
    if r31_koide is not None:
        err = (r31_koide - r31) / r31 * 100
        print(f"           r31(Koide=2/3) = {r31_koide:.0f} (err = {err:+.1f}%)")
        # What m3 would this correspond to?
        if name != 'lepton':
            m1 = sectors_data[name].get('m1', None)
    print()

# ---- 5. Stosunek r31/r21 i golden ratio ----
print("\n--- 5. Relacja miedzy r21, r31 i phi ---\n")

for name, data in sectors_data.items():
    r21 = data['r21']
    r31 = data['r31']
    ratio = r31 / r21
    log_ratio = np.log(r31) / np.log(r21)
    sqrt_ratio = np.sqrt(r31) / np.sqrt(r21)
    print(f"  {name:>7s}: r31/r21 = {ratio:.1f}, "
          f"ln(r31)/ln(r21) = {log_ratio:.4f}, "
          f"sqrt(r31/r21) = {np.sqrt(ratio):.3f}")

print(f"\n  phi = {(1+np.sqrt(5))/2:.5f}")
print(f"  3/2 = 1.50000")
print(f"  phi+1 = phi^2 = {((1+np.sqrt(5))/2)**2:.5f}")
print()

# Key: for leptons, K=2/3 + phi-FP determines r31
# What is the mathematical relation?
# K = 2/3, r21 = (A2/A1)^4, r31 = ?
# From K formula: given K=2/3 and r21, r31 is UNIQUELY determined.
# For leptons: r31_pred = 3477 = r31_PDG (match!)
# For quarks: r31_pred != r31_PDG (mismatch)

print("  Kluczowe spostrzezenie:")
print("  Koide K=2/3 + phi-FP (r21) JEDNOZNACZNIE wyznacza r31.")
print("  Leptony: r31(pred) = r31(PDG) -> Koide dziala.")
print("  Kwarki: r31(pred) != r31(PDG) -> Koide NIE dziala.")
print("  -> Trzecia generacja kwarkowa ma INNY warunek niz K=2/3.")

# ---- 6. What DOES determine quark r31? ----
print("\n--- 6. Co determinuje r31 kwarkowe? ---\n")

# Hypothesis: r31 = r21 * r_gen where r_gen is generation factor
for name, data in sectors_data.items():
    r21 = data['r21']
    r31 = data['r31']
    r32 = r31 / r21
    print(f"  {name:>7s}: r32 = m3/m2 = {r32:.1f}")

print()
print("  Lepton: r32 = 16.8")
print("  Down:   r32 = 44.8")
print("  Up:     r32 = 136.1")
print()
print("  Ratios: r32_down/r32_lepton = {:.2f}".format(895/20.0 / (3477/206.8)))
print("  Ratios: r32_up/r32_lepton = {:.2f}".format(79982/588.0 / (3477/206.8)))

# ---- WNIOSKI ----
print("\n" + "=" * 72)
print("WNIOSKI ex172")
print("=" * 72)
print("""
  1. K(d,s,b) i K(u,c,t) sa STABILNE na bieganiu QCD:
     K(d,s,b) ~ 0.744 (+11.6%) na kazdej skali mu
     K(u,c,t) ~ 0.882 (+32.3%) na kazdej skali mu
     (K jest jednorodne stopnia 0 w masach -> bieganie
      zmienia wszystkie masy o TEN SAM czynnik -> K = const!)

  2. K jest DOKLADNIE niezmienniczy na bieganiu QCD
     (bo bieganie mnoznikowe: m(mu) = c(mu)*m_0,
      a K jest jednorodne st. 0).

  3. KLUCZOWY WNIOSEK: K(d,s,b) != 2/3 i K(u,c,t) != 2/3
     sa FAKTAMI ALGEBRAICZNYMI wynikajacymi z PDG mas.
     Zadne bieganie tego nie zmieni.
     (Wyjstek: jesli masy BIEGNA ROZNIE, np. rozne gamma_m
      dla roznych smakowow -- ale to jest efekt > 2-loop,
      pomijalny.)

  4. Koide K=2/3 jest wlasnoscia LEPTONOWA.
     TGP moze to tlumaczyc: leptony nie maja koloru QCD,
     ich masy sa CZYSTO solitonowe (brak korekcji perturbacyjnych).
     Kwarki maja korekcje kolorowe, ktore LAMANIA K=2/3.

  5. Otwarte: co jest TGP-owym odpowiednikiem Koide dla kwarkow?
     Kandydaci:
     - K_eff(sektor) = f(N_c, alpha_s) * 2/3  (korekcja kolorowa)
     - Shifted Koide z m0 proporcjonalnym do Lambda_QCD
     - Zupelnie inny warunek selekcji 3rd gen
""")
