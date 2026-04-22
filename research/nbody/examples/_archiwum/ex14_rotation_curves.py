"""
ex14_rotation_curves.py
========================
Krzywa rotacji galaktyk: TGP (Droga B) vs Newton vs MOND vs CDM.

PYTANIE:
  Czy TGP z Yukawa (m_sp ~ kpc^{-1}) jest odroznialne od CDM/MOND
  w obserwowanej krzywej rotacji?

MODELE:
  1. Newton only (baryony):    v^2 = G*M_bar(<r) / r
  2. CDM (NFW profile):        v^2 = v_Newton^2 + v_NFW^2
  3. MOND (Milgrom):           mu(a/a0)*a = a_Newton  [BM interpolation]
  4. TGP (Droga B):            G_eff(r) = G * f_TGP(m_sp, r)

GALAKTYKA REFERENCYJNA:
  NGC 3198 — klasyczny przyklad plasty krzywej rotacji.
  Dane: Begeman 1989 / de Blok et al 2008 (THINGS).
  Baryony: dysk eksponencjalny Rd ~ 3.2 kpc, v_bar_max ~ 100 km/s.
  Obserwowana: v_obs ~ 150 km/s dla r = 5-30 kpc (plaska).

DANE NGC 3198 (przyblizone, z Begeman 1989):
  r [kpc]:  0.5  1.0  2.0  3.0  5.0  7.0 10.0 15.0 20.0 25.0 30.0
  v [km/s]: 90  130  155  160  155  152  150  148  147  146  145
"""

import sys, os
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import k0, k1, i0, i1
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# ===========================================================================
# Dane obserwacyjne NGC 3198 (Begeman 1989 / THINGS)
# ===========================================================================
OBS_R_KPC = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0])
OBS_V_KMS = np.array([90., 130., 155., 160., 155., 152., 150., 148., 147., 146., 145.])
OBS_ERR   = np.array([ 8.,   6.,   5.,   5.,   4.,   4.,   4.,   4.,   5.,   6.,   7.])

# Parametry fizyczne
G_SI   = 6.674e-11      # m^3 kg^-1 s^-2
kpc_m  = 3.086e19       # m
M_sun  = 1.989e30       # kg
km_s   = 1e3            # m/s

# NGC 3198 parametry baryonowe (Begeman 1989 / de Blok 2008)
M_disk  = 2.0e10 * M_sun  # masa dysku [kg] (maximum disk fit)
R_d     = 3.2 * kpc_m     # skala dysku [m]

# ISO halo dla NGC 3198 (Begeman 1989 - klasyczny fit)
# Halo izotermiczne: rho(r) = rho0 / (1 + (r/rc)^2)
# Parametry: rho0 = 0.054 M_sun/pc^3, rc = 4.0 kpc (Begeman 1989)
pc_m       = 3.086e16       # m
rho0_ISO   = 0.054 * M_sun / pc_m**3   # [kg/m^3], Begeman 1989
r_c_ISO    = 4.0 * kpc_m               # [m] core radius
# NFW rowniez, dla porownan (de Blok 2008 THINGS)
rho0_NFW = 6.7e6 * M_sun / kpc_m**3  # [kg/m^3] = 6.7e6 M_sun/kpc^3
r_s_NFW  = 12.0 * kpc_m               # [m] skala NFW

print("=" * 70)
print("KRZYWA ROTACJI GALAKTYK: TGP vs CDM vs MOND")
print("=" * 70)
print()
print("Galaktyka: NGC 3198")
print(f"  Masa dysku: M_disk = {M_disk/M_sun:.2e} M_sun")
print(f"  Skala dysku: R_d   = {R_d/kpc_m:.1f} kpc")
print(f"  ISO: rho0 = {rho0_ISO*pc_m**3/M_sun:.4f} M_sun/pc^3, rc = {r_c_ISO/kpc_m:.1f} kpc")
print(f"  NFW: rho0 = {rho0_NFW*kpc_m**3/M_sun:.2e} M_sun/kpc^3, r_s = {r_s_NFW/kpc_m:.1f} kpc")
print()

# ===========================================================================
# Model 1: Dysk eksponencjalny (baryony only)
# ===========================================================================
def v_disk(r_m, M_disk, R_d):
    """
    Predkosc kryzowa dla dysku eksponencjalnego (Freeman 1970).
    v^2 = 4*pi*G*Sigma0*R_d * y^2 * [I0(y)*K0(y) - I1(y)*K1(y)]
    y = r / (2*R_d)
    """
    Sigma0 = M_disk / (2.0 * np.pi * R_d**2)
    y = r_m / (2.0 * R_d)
    y = np.where(y < 1e-10, 1e-10, y)
    # Bessele cylindryczne
    v2 = (4.0 * np.pi * G_SI * Sigma0 * R_d
          * y**2 * (i0(y)*k0(y) - i1(y)*k1(y)))
    return np.sqrt(np.maximum(v2, 0.0))

# ===========================================================================
# Model 2: ISO halo (izotermiczne) + NFW
# ===========================================================================
def v_ISO(r_m, rho0, r_c):
    """Predkosc kryzowa dla izotermicznego halo (Pseudo-ISO)."""
    v2 = 4.0*np.pi*G_SI*rho0*r_c**2 * (1.0 - r_c/r_m * np.arctan(r_m/r_c))
    return np.sqrt(np.maximum(v2, 0.0))

def M_NFW(r_m, rho0, r_s):
    """Masa NFW w kuli promienia r."""
    x = r_m / r_s
    return 4.0 * np.pi * rho0 * r_s**3 * (np.log(1.0 + x) - x/(1.0+x))

def v_NFW(r_m, rho0, r_s):
    """Predkosc kryzowa dla halo NFW."""
    M = M_NFW(r_m, rho0, r_s)
    return np.sqrt(G_SI * M / r_m)

# ===========================================================================
# Model 3: MOND (Milgrom 1983)
# ===========================================================================
a0_MOND = 1.2e-10  # m/s^2  (stala Milgroma)

def v_MOND(r_m, a_Newton_func):
    """
    Predkosc MOND z interpolacja BM (Begeman-Milgrom):
      mu(x) = x/sqrt(1+x^2), x = a/a0
    Rozwiazanie: a * mu(a/a0) = a_N => a^2 = a_N^2 + a_N*a0 (approx)
    Dokladniej: a^2(1 + a0^2/a^2)^{1/2} = a_N  -- iteracyjnie
    """
    # Prosta formula: v_MOND^2 = v_N^2 * (1 + sqrt(1 + 4*a0*r/v_N^2)) / 2
    # (z mu-interpolacji BM w granicy prostej)
    a_N = a_Newton_func(r_m)  # [m/s^2]
    # mu(x)*x = a_N, mu(x) = x/sqrt(1+x^2) -> x = a/a0
    # x^2/(1+x^2)^{1/2} = a_N/a0 => x^4 - (a_N/a0)^2 * (1+x^2) = 0
    # rozwiazanie: x^2 = 0.5*[(a_N/a0)^2 + sqrt((a_N/a0)^4 + 4*(a_N/a0)^2)]
    u = (a_N / a0_MOND)**2
    x2 = 0.5 * (u + np.sqrt(u**2 + 4.0*u))
    a_MOND = np.sqrt(x2) * a0_MOND
    v2 = a_MOND * r_m
    return np.sqrt(np.maximum(v2, 0.0))

# ===========================================================================
# Model 4: TGP z Yukawa (Droga B)
# ===========================================================================
def v_TGP(r_m, M_bar, R_d, m_sp_inv_m):
    """
    Predkosc TGP: potencjal parowy = Yukawa + dysk eksponencjalny.

    W TGP: G_eff(r) = G * exp(-r/lambda)
    Ale dysk nie jest punktem! Dla dysku eksponencjalnego z Yukawa:

    V_TGP(r) = -G_eff(r) * M_eff(r) / r  [uproszczenie]
    gdzie M_eff uwzglednia zanik Yukawa.

    Bardziej precyzyjnie: potencjal Yukawa dysku calkujemy numerycznie.
    v_TGP^2(r) = r * |dPhi/dr|

    Phi(r) = -G * Sigma0 * int_0^inf int_0^{2pi} Sigma(r')*exp(-m_sp*|r-r'|)/|r-r'| r' dr' dphi
    """
    if m_sp_inv_m == 0.0:
        # Granica Newtona
        return v_disk(r_m, M_bar, R_d)

    Sigma0 = M_bar / (2.0 * np.pi * R_d**2)
    lam = 1.0 / m_sp_inv_m  # dlugosc ekranowania [m]

    # Dla dysku eksponencjalnego z Yukawa:
    # v^2(r) = G * Sigma0 * r * int_0^inf Sigma(r')/r * exp(-m_sp*r) * f(r,r') dr'
    # Uzywamy dyskretnej quadratury po dysku

    # Siatka promieniowa
    r_max = 15.0 * R_d
    r_arr = np.linspace(0.01*R_d, r_max, 200)

    # Sigma(r) = Sigma0 * exp(-r/R_d)
    Sig = Sigma0 * np.exp(-r_arr / R_d)

    # Potencjal Yukawa z pierscienia o promieniu r' w punkcie r na plaszczyznie:
    # Phi_ring(r, r') = -2*G*Sigma*r' * int_0^pi exp(-m_sp*sqrt(r^2+r'^2-2rr'cos(phi))) /
    #                               sqrt(r^2+r'^2-2rr'cos(phi)) dphi
    # Ta calka to modyfikowana funkcja Bessela:
    # = -2*G*Sigma*r' * 2*K0(m_sp*sqrt|r-r'|) * exp(-m*sqrt{...}) -- brak prostego wzoru
    # Uzyjemy rozwinięcia calkowego numerycznego dla kilku r

    def phi_yukawa_at(r_eval):
        """Potencjal Yukawa dysku w punkcie r_eval."""
        n_phi = 32
        phi_arr = np.linspace(0, 2.0*np.pi, n_phi, endpoint=False)
        dphi = 2.0*np.pi / n_phi
        dr_r = r_arr[1] - r_arr[0]

        total = 0.0
        for rp in r_arr:
            sig_rp = Sigma0 * np.exp(-rp / R_d)
            for phi_v in phi_arr:
                dist = np.sqrt(r_eval**2 + rp**2 - 2.0*r_eval*rp*np.cos(phi_v))
                dist = max(dist, 0.01*R_d)
                total += -G_SI * sig_rp * rp * np.exp(-m_sp_inv_m*dist) / dist * dphi * dr_r
        return total

    # Sila i predkosc: v^2 = r * |dPhi/dr| przez rozniczkowanie
    # Uzywamy dwoch punktow
    h = 0.005 * R_d
    Phi_p = phi_yukawa_at(r_m + h)
    Phi_m = phi_yukawa_at(r_m - h) if r_m > h else phi_yukawa_at(h)
    dPhi = (Phi_p - Phi_m) / (2.0*h if r_m > h else h)
    v2 = max(r_m * (-dPhi), 0.0)
    return np.sqrt(v2)

def v_TGP_fast(r_m, M_bar, R_d, lam_m):
    """
    Szybka aproksymacja v_TGP dla dysku eksponencjalnego.

    Kluczowy wynik: dla lambda >> R_d (Yukawa nie zanika na skali dysku):
      v_TGP^2(r) ~ v_disk^2(r) * G_eff(r)/G * F(r, lambda)

    gdzie F(r,lambda) ~ (1 + r/lambda + (r/lambda)^2/2) * exp(-r/lambda) [geometryczne]

    Dla lambda << R_d: sila zanika na skali R_d => brak grawitacji!
    Dla lambda >> R_d: Newton recover.

    Precyzyjniejsza formula dla dysku z Yukawa (Milgrom & Sanders 2008 style):
    v_TGP^2(r) = G_eff(r)/G * v_disk^2(r)  [zeroth order]
    + korekta z calki powierzchniowej [order lambda/r]

    Tutaj uzywamy: v_TGP^2 = v_disk^2 * exp(-r/lambda)
    To jest DOLNE oszacowanie (faktyczny efekt nieco silniejszy dla r ~ lambda)
    """
    if lam_m == 0.0 or np.isinf(lam_m):
        return v_disk(r_m, M_bar, R_d)

    m_sp = 1.0 / lam_m
    # Dla dysku eksponencjalnego: calka analityczna przez Bessel
    # v_TGP^2(r) = G*Sigma0 * r * d/dr [-Phi_Y(r)]
    # Phi_Y(r) = int Sigma(r') * K0(m_sp * |r-r'|) * ... (nie ma prostej formy)
    #
    # Aproksymacja robocza (poprawna dla r >> R_d):
    # v_TGP^2(r) ~ v_disk^2(r) * exp(-m_sp * r) * (1 + m_sp*r)^{1/2}
    # (heurystyczne, ale odtwarza granice Newton i pelny zanik)
    v_d = v_disk(r_m, M_bar, R_d)
    t = m_sp * r_m
    # Efektywna modyfikacja: sredniowanie po dysku daje mniejsza supresje niz punkt
    # Uzyj: <exp(-m*|r-r'|)> ~ exp(-m*r) * cosh(m*R_d) dla r >> R_d
    corr = np.exp(-t) * np.cosh(min(m_sp * R_d, 50.0))
    corr = min(corr, 1.0)  # nie moze byc wieksze od Newtona
    return v_d * np.sqrt(corr)

# ===========================================================================
# Obliczenia dla siatki r
# ===========================================================================
r_kpc_arr = np.linspace(0.5, 32.0, 120)
r_m_arr   = r_kpc_arr * kpc_m

print("Obliczam modele predkosci...")

# Newton (baryony only)
v_bar_kms = v_disk(r_m_arr, M_disk, R_d) / km_s

# ISO halo (izotermiczne pseudo-sferyczne)
v_iso_kms = v_ISO(r_m_arr, rho0_ISO, r_c_ISO) / km_s

# CDM = dysk + ISO halo (Begeman 1989 fit)
v_cdm_kms = np.sqrt(v_bar_kms**2 + v_iso_kms**2)

# NFW osobno dla porownania
v_nfw_kms = v_NFW(r_m_arr, rho0_NFW, r_s_NFW) / km_s
v_cdm_nfw_kms = np.sqrt(v_bar_kms**2 + v_nfw_kms**2)

print(f"  v_bar(5kpc)={np.interp(5, r_kpc_arr, v_bar_kms):.1f},"
      f" v_iso(5kpc)={np.interp(5, r_kpc_arr, v_iso_kms):.1f},"
      f" v_cdm(5kpc)={np.interp(5, r_kpc_arr, v_cdm_kms):.1f} km/s")
print(f"  v_bar(20kpc)={np.interp(20, r_kpc_arr, v_bar_kms):.1f},"
      f" v_iso(20kpc)={np.interp(20, r_kpc_arr, v_iso_kms):.1f},"
      f" v_cdm(20kpc)={np.interp(20, r_kpc_arr, v_cdm_kms):.1f} km/s")

# MOND
def a_bar(r_m):
    vd = v_disk(r_m, M_disk, R_d)
    return vd**2 / r_m  # centrypetalne od baryonow

v_mond_kms = np.array([v_MOND(r, a_bar) / km_s for r in r_m_arr])

# TGP dla roznych lambda
lambda_vals_kpc = [1.0, 5.0, 10.0, 30.0, 100.0]  # kpc
v_tgp_models = {}
for lam_kpc in lambda_vals_kpc:
    lam_m = lam_kpc * kpc_m
    v_arr = np.array([v_TGP_fast(r, M_disk, R_d, lam_m) / km_s for r in r_m_arr])
    v_tgp_models[lam_kpc] = v_arr
    print(f"  TGP lambda={lam_kpc:5.0f} kpc: v(5kpc)={v_arr[np.argmin(np.abs(r_kpc_arr-5))]:.1f},"
          f" v(20kpc)={v_arr[np.argmin(np.abs(r_kpc_arr-20))]:.1f} km/s")

print()

# ===========================================================================
# Chi^2 dopasowania do danych NGC 3198
# ===========================================================================
print("=" * 70)
print("DOPASOWANIE DO DANYCH NGC 3198")
print("=" * 70)
print()

def chi2(v_model_kms, obs_r, obs_v, obs_err):
    """Chi^2 interpolujac model w punktach obserwacyjnych."""
    v_interp = np.interp(obs_r, r_kpc_arr, v_model_kms)
    return np.sum(((v_interp - obs_v) / obs_err)**2) / len(obs_v)

chi2_bar     = chi2(v_bar_kms,     OBS_R_KPC, OBS_V_KMS, OBS_ERR)
chi2_cdm     = chi2(v_cdm_kms,     OBS_R_KPC, OBS_V_KMS, OBS_ERR)
chi2_cdm_nfw = chi2(v_cdm_nfw_kms, OBS_R_KPC, OBS_V_KMS, OBS_ERR)
chi2_mond    = chi2(v_mond_kms,    OBS_R_KPC, OBS_V_KMS, OBS_ERR)

print(f"  {'Model':30s}  {'chi2/N':>10s}  {'Akceptowalny?':>15s}")
print("-" * 60)
def accept_str(c2):
    return 'TAK (chi2~1)' if c2 < 2 else ('OK (chi2<5)' if c2 < 5 else 'ZLE')

print(f"  {'Newton (baryony only)':30s}  {chi2_bar:>10.2f}  {'NIE (za wolna)':>15s}")
print(f"  {'CDM (ISO halo, Begeman)':30s}  {chi2_cdm:>10.2f}  {accept_str(chi2_cdm):>15s}")
print(f"  {'CDM (NFW halo)':30s}  {chi2_cdm_nfw:>10.2f}  {accept_str(chi2_cdm_nfw):>15s}")
print(f"  {'MOND':30s}  {chi2_mond:>10.2f}  {accept_str(chi2_mond):>15s}")

best_lam = None
best_chi2 = float('inf')
for lam_kpc, v_arr in v_tgp_models.items():
    c2 = chi2(v_arr, OBS_R_KPC, OBS_V_KMS, OBS_ERR)
    accept = 'TAK' if c2 < 2 else ('OK' if c2 < 5 else 'ZLE')
    print(f"  {'TGP lambda='+str(int(lam_kpc))+' kpc':30s}  {c2:>10.2f}  {accept:>15s}")
    if c2 < best_chi2:
        best_chi2 = c2
        best_lam = lam_kpc

print()
print(f"  Najlepsze TGP: lambda = {best_lam:.0f} kpc  (chi2/N = {best_chi2:.2f})")
print()

# ===========================================================================
# Analiza: co TGP musi dostarczyc
# ===========================================================================
print("=" * 70)
print("ANALIZA: Co TGP musi dostarczyc zamiast CDM?")
print("=" * 70)
print()
print("Problem: TGP Drogi B modyfikuje sile Newtona przez Yukawa G_eff(r).")
print("Ale G_eff(r) = G*exp(-r/lambda) MALEJE z r.")
print("=> Krzywa rotacji TGP jest NIZSZA od Newtona dla duzych r.")
print("=> Obserwowana plaska krzywa rotacji wymaga dodatkowej masy!")
print()
print("WNIOSEK 1: Sam Yukawa bez dodatkowej masy NIE wyjasmia plasty krzywej.")
print("           TGP ze skalarnym Yukawa daje MNIEJ grawitacji, nie wiecej.")
print()
print("WNIOSEK 2: TGP moze pomoc tylko jesli ODWROCI kierunek dzialania.")
print("           Mozliwosci:")
print("  (a) m_sp wyimaginowane (tachioniczne): Yukawa -> sinusoidalne")
print("      Potencjal: exp(+m|r|)/r (rosnacy!) - unphysical")
print()
print("  (b) Masa osłony TGP: Phi rownowazy sie z kondensatem =")
print("      'morze Yukawa' wokol galaktyki => efektywna masa ciemna")
print("      (ale to wymaga odrebnego pola kondensatu)")
print()
print("  (c) TGP nie jest teoria grawitacji [Scenariusz B/C].")
print("      Grawitacja jest GR (Newton). TGP to nowa sila ponad Newtonem.")
print("      Ale TGP nie daje plasty krzywej bezposrednio.")
print()
print("  (d) TGP z ujemnym m_sp^2 (tachion instability):")
print("      m_sp^2 < 0 => inna ODE => rosnacy potencjal? Wymaga analizy.")
print()

# ===========================================================================
# Model alternatywny: TGP jako skalarne halo kondensatu
# ===========================================================================
print("=" * 70)
print("MODEL ALTERNATYWNY: Pole TGP jako kondensatem ('fuzzy dark matter')")
print("=" * 70)
print()
print("Jesli pole Phi nie zeruje sie asymptotycznie, ale tworzy 'kondensate':")
print("  Phi_bg(r) = Phi_0 * J_0(m_sp*r) / r  [dla n=0]")
print("  (soliton skalarny, 'boson star' model)")
print()
print("Masa efektywna kondensatu TGP w kuli r:")
print("  M_cond(<r) = 4*pi * int_0^r rho_Phi(r') r'^2 dr'")
print("  rho_Phi = 1/2 * (nabla Phi)^2 + V(Phi)")
print()
print("To jest 'fuzzy dark matter' (FDM) z bosonu Phi_TGP.")
print("Predykcja FDM: v_rot wyplaszcza sie przy r ~ lambda/2.")
print()
print("UWAGA: FDM jest znana hipoteza (Hu & Sawicki 2000, Hui et al 2017).")
print("       TGP dostarcza NATURALNEGO kantydata na FDM!")
print("       m_boson_FDM ~ 1e-22 eV/c^2 => lambda ~ 1 kpc")
print()

m_FDM_eV = 1e-22   # eV
hbar_eV_s = 6.582e-16  # eV*s
c_SI_val = 3e8
m_FDM_kg  = m_FDM_eV * 1.602e-19 / c_SI_val**2
lam_deBroglie = hbar_eV_s * c_SI_val / (m_FDM_eV * 1e-15)  # Compton ~ de Broglie
print(f"  m_boson = {m_FDM_eV:.0e} eV/c^2")
print(f"  Lambda_Compton ~ hbar*c / (m*c^2) = {lam_deBroglie:.1e} m = {lam_deBroglie/kpc_m:.2e} kpc")
print()

# Masa TGP Plancka odpowiadajaca m_FDM
l_Planck_m = 1.616e-35
m_Planck_kg = 2.176e-8
m_sp_FDM_Planck = m_FDM_kg / m_Planck_kg * (c_SI_val/l_Planck_m)
# m_sp in Planck: m_sp = m_boson*c/hbar * l_Planck
m_sp_FDM_P = m_FDM_kg * c_SI_val * l_Planck_m / (1.055e-34)
print(f"  m_sp (TGP Planck) = {m_sp_FDM_P:.3e}")
print(f"  Czy zgodne z ograniczeniem ex13 (m_sp > {3.4e-41:.0e})?",
      "TAK" if m_sp_FDM_P > 3.4e-41 else "NIE")
print()

# ===========================================================================
# Porownanie graficzne
# ===========================================================================
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Lewy: krzywe rotacji
ax = axes[0]
ax.errorbar(OBS_R_KPC, OBS_V_KMS, yerr=OBS_ERR,
            fmt='ko', ms=5, capsize=3, label='NGC 3198 (obs.)', zorder=10)
ax.plot(r_kpc_arr, v_bar_kms,     'k--', lw=1.5, alpha=0.7, label='Newton (baryony)')
ax.plot(r_kpc_arr, v_cdm_kms,     'b-',  lw=2.5, label=f'CDM ISO ($\\chi^2$={chi2_cdm:.1f})')
ax.plot(r_kpc_arr, v_cdm_nfw_kms, 'b--', lw=1.5, alpha=0.7, label=f'CDM NFW ($\\chi^2$={chi2_cdm_nfw:.1f})')
ax.plot(r_kpc_arr, v_mond_kms,    'g-',  lw=2, label=f'MOND ($\\chi^2$={chi2_mond:.1f})')

colors = cm.Reds(np.linspace(0.4, 0.9, len(lambda_vals_kpc)))
for (lam_kpc, v_arr), col in zip(v_tgp_models.items(), colors):
    c2 = chi2(v_arr, OBS_R_KPC, OBS_V_KMS, OBS_ERR)
    ax.plot(r_kpc_arr, v_arr, '-', color=col, lw=1.5,
            label=rf'TGP $\lambda={lam_kpc:.0f}$ kpc ($\chi^2$={c2:.1f})')

ax.set_xlabel('r [kpc]', fontsize=12)
ax.set_ylabel('v [km/s]', fontsize=12)
ax.set_title('NGC 3198: Krzywa rotacji', fontsize=12)
ax.legend(fontsize=8, loc='lower right')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 33)
ax.set_ylim(0, 220)

# Prawy: chi2 vs lambda dla TGP
ax = axes[1]
lam_scan = np.logspace(-1, 3, 50)  # kpc
chi2_scan = []
for lam_kpc in lam_scan:
    lam_m = lam_kpc * kpc_m
    v_arr = np.array([v_TGP_fast(r, M_disk, R_d, lam_m) / km_s for r in r_m_arr])
    chi2_scan.append(chi2(v_arr, OBS_R_KPC, OBS_V_KMS, OBS_ERR))

ax.semilogx(lam_scan, chi2_scan, 'r-', lw=2, label='TGP')
ax.axhline(chi2_cdm,  color='b', ls='--', lw=1.5, label=f'CDM ($\\chi^2$={chi2_cdm:.2f})')
ax.axhline(chi2_mond, color='g', ls='--', lw=1.5, label=f'MOND ($\\chi^2$={chi2_mond:.2f})')
ax.axhline(2.0, color='gray', ls=':', lw=1, alpha=0.5, label='granica akceptacji')
ax.set_xlabel(r'$\lambda_{TGP}$ [kpc]', fontsize=12)
ax.set_ylabel(r'$\chi^2/N$', fontsize=12)
ax.set_title(r'Dopasowanie TGP do NGC 3198 vs $\lambda$', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 30)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), 'ex14_rotation_curves.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"Wykres: {out}")
print()

# ===========================================================================
# Podsumowanie diagnostyki
# ===========================================================================
print("=" * 70)
print("DIAGNOSTYKA: TGP jako teoria ciemnej materii?")
print("=" * 70)
print()
print("WYNIK: Proste Yukawa TGP NIE moze zastapic CDM.")
print("       Yukawa G_eff(r) < G zmniejsza sile, nie zwieksza.")
print()
print("DROGI WYJSCIA:")
print()
print("  1. TGP + CDM (dwa skladniki):")
print("     Baryony + CDM (NFW) jako oddzielne pola.")
print("     TGP dostarcza trojcialowe korekty do wewnetrznej struktury halo.")
print()
print("  2. TGP kondensatem (FDM):")
print("     Pole Phi tworzy soliton galaktyczny.")
print("     m_boson ~ 1e-22 eV odpowiada lambda ~ 1 kpc.")
print("     TGP boson = naturalny FDM kandydat.")
print("     Wymaga rozszerzonej teorii: V(Phi) musi miec minimum inne niz Phi=0.")
print()
print("  3. Zmodyfikowany V(Phi) (tachionowy):")
print("     Jesli m_sp^2 < 0: potencjal bez minimum w Phi=0.")
print("     Phi kondensuje na skali galaktycznej.")
print("     To jest scenariusz 'DE boson' lub 'chameleon field'.")
print()
print("  4. TGP tylko na skalach subnuklearnych (m_sp ~ 1):")
print("     Grawitacja = GR + baryony + standardowe CDM.")
print("     TGP = nowa sila dla procesow na skali Plancka.")
print("     Plaska krzywa rotacji wymaga innego mechanizmu.")
print()
print("REKOMENDACJA:")
print("  Najciekawszy scenariusz: TGP boson jako FDM z m ~ 1e-22 eV.")
print("  Wymaga analizy solitonowej (boson star profile).")
print("  Nastepny krok: ex15_tgp_soliton.py")
