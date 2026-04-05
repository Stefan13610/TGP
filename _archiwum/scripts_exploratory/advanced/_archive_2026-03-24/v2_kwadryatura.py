"""
v2_kwadryatura.py - Diagnoza: czy niezbieznosc energii M3 to blad kwadratury?

Problem z v2: calka energii dla M3 (K~47, psi_core~2442) nie zbiega z N.
Hipoteza: integrand ma ostre piecie przy r=a_gam (lambda*(psi-1)^6 ~ r^-6),
          wiec rownoomierna siatka wymaga bardzo duzego N.

Test: siatka logarytmiczna (gestsze punkty przy r=a_gam)
vs. rownoomierna => czy wyniki sie zgadzaja?
"""
import numpy as np
from scipy.integrate import quad
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0
LAM   = 5.514e-7
ALPHA = 7.0
A_GAM = 0.030

def V_mod(psi, lam=LAM):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

V1 = V_mod(1.0)

def integrand_Ek(r, K, alpha, msp=1.0):
    phi  = max(1.0 + K*np.exp(-msp*r)/r, 1e-10)
    dphi = K*np.exp(-msp*r)*(-msp*r - 1.0)/r**2
    return 4*np.pi * 0.5 * dphi**2 * (1 + alpha/phi) * r**2

def integrand_Ep(r, K, lam=LAM, msp=1.0):
    phi  = max(1.0 + K*np.exp(-msp*r)/r, 1e-10)
    psi  = phi  # Phi0 = 1
    return 4*np.pi * (V_mod(psi, lam) - V1) * r**2

def energy_uniform(K, alpha, a_gam, N, lam=LAM):
    r_max = max(60.0, 20.0)
    r = np.linspace(a_gam, r_max, N)
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    psi  = phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    Ep   = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
    return Ek, Ep, Ek+Ep

def energy_log(K, alpha, a_gam, N, lam=LAM):
    """Siatka logarytmiczna: gestsza przy r=a_gam."""
    r_max = max(60.0, 20.0)
    t     = np.linspace(0, 1, N)
    r     = a_gam * (r_max/a_gam)**t   # r = a_gam * exp(t*log(r_max/a_gam))
    phi   = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi  = K*np.exp(-r)*(-r - 1.0)/r**2
    psi   = phi
    Ek    = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
    Ep    = 4*np.pi*np.trapezoid((V_mod(psi,lam)-V1)*r**2, r)
    return Ek, Ep, Ek+Ep

def energy_quad(K, alpha, a_gam, lam=LAM, limit=200):
    """Adaptywna kwadryatura scipy.integrate.quad."""
    Ek, _ = quad(integrand_Ek, a_gam, 60.0, args=(K, alpha),
                 limit=limit, epsabs=1e-4, epsrel=1e-6)
    Ep, _ = quad(integrand_Ep, a_gam, 60.0, args=(K, lam),
                 limit=limit, epsabs=1e-4, epsrel=1e-6)
    return Ek, Ep, Ek+Ep

# ── Wartosci K do przetestowania ──────────────────────────────────────────
# K3 z v2_fine_scan (r31=3477): K3 ~ 47.32 (v2 best)
# Ale g(K) shows 3rd zero between K=60 and K=80 for alpha=7, a=0.03
# Wiec sprawdzimy oba
K_test_cases = [
    (47.32, "M3 z v2_fine_scan (r31=3477)"),
    (70.0,  "K ~ srodek okna 3. zera (K=60..80)"),
]

for K_test, label in K_test_cases:
    psi_c = 1.0 + K_test * np.exp(-A_GAM) / A_GAM
    print("=" * 70)
    print(f"K = {K_test}, psi_core ~ {psi_c:.1f}")
    print(f"Opis: {label}")
    print()

    # Rownoomierna siatka
    print(f"{'Siatka':12} {'N':>6} {'Ek':>14} {'Ep':>14} {'Etot':>14}")
    print("-" * 60)
    for N in [500, 1000, 2000, 4000, 8000, 16000]:
        Ek, Ep, Et = energy_uniform(K_test, ALPHA, A_GAM, N)
        print(f"{'rownoomierna':12} {N:>6} {Ek:>14.4e} {Ep:>14.4e} {Et:>14.4e}")

    print()
    # Siatka logarytmiczna
    for N in [500, 1000, 2000, 4000, 8000]:
        Ek, Ep, Et = energy_log(K_test, ALPHA, A_GAM, N)
        print(f"{'logarytmiczna':12} {N:>6} {Ek:>14.4e} {Ep:>14.4e} {Et:>14.4e}")

    print()
    # Adaptywna kwadryatura
    try:
        Ek, Ep, Et = energy_quad(K_test, ALPHA, A_GAM, limit=300)
        print(f"{'quad(adaptyw)':12} {'---':>6} {Ek:>14.4e} {Ep:>14.4e} {Et:>14.4e}")
    except Exception as e:
        print(f"quad error: {e}")
    print()

# ── Analiza integrandy przy r=a_gam ───────────────────────────────────────
print("=" * 70)
print("ANALIZA INTEGRANDY w poblizy r=a_gam:")
print()
print(f"{'r/a_gam':>10} {'psi':>10} {'Ep_integrand':>16} {'Ek_integrand':>16}")
print("-" * 55)
K_diag = 47.32
for frac in [1.0, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 2.0, 5.0, 10.0]:
    r = A_GAM * frac
    phi  = max(1.0 + K_diag * np.exp(-r) / r, 1e-10)
    dphi = K_diag * np.exp(-r) * (-r - 1.0) / r**2
    psi  = phi
    ep_int = 4*np.pi * (V_mod(psi) - V1) * r**2
    ek_int = 4*np.pi * 0.5 * dphi**2 * (1 + ALPHA/phi) * r**2
    print(f"{frac:>10.2f} {psi:>10.1f} {ep_int:>16.4e} {ek_int:>16.4e}")

print()
print("Jesli Ep_integrand zmienia znak gwaltownie => problem z kwadryatura")
print("Jesli Ep_integrand jest gladka => problem numeryczny jest gdzie indziej")

# ── Rozbicie calki na przedzial [a_gam, r_mid] i [r_mid, infty] ───────────
print()
print("=" * 70)
print("ROZBICIE CALKI Ep na przedzialy:")
K_diag = 47.32
lam    = LAM
r_max  = 60.0
N_pts  = 8000

r_full = np.linspace(A_GAM, r_max, N_pts)
phi_f  = np.maximum(1.0 + K_diag * np.exp(-r_full)/r_full, 1e-10)
psi_f  = phi_f
integrand_f = 4*np.pi * (V_mod(psi_f, lam) - V1) * r_full**2

# Znajdz gdzie integrand zmienia znak
sign_arr = np.sign(integrand_f)
zero_crossings = np.where(np.diff(sign_arr) != 0)[0]
print(f"Zmiana znaku integrandy Ep przy r = ", end="")
for idx in zero_crossings:
    print(f"{r_full[idx]:.4f} (psi={psi_f[idx]:.1f})", end="  ")
print()

# Calki na przedzialach
for i in range(len(zero_crossings)):
    i_start = 0 if i == 0 else zero_crossings[i-1]
    i_end   = zero_crossings[i]
    r_seg   = r_full[i_start:i_end+1]
    int_seg = integrand_f[i_start:i_end+1]
    seg_int = np.trapezoid(int_seg, r_seg)
    psi_lo  = psi_f[i_start]
    psi_hi  = psi_f[i_end]
    print(f"  Przedzial r=[{r_full[i_start]:.3f},{r_full[i_end]:.3f}], "
          f"psi=[{psi_hi:.1f},{psi_lo:.1f}]: Ep_seg = {seg_int:.4e}")

# Pozostaly przedzial
if len(zero_crossings) > 0:
    i_last = zero_crossings[-1]
    r_seg  = r_full[i_last:]
    int_seg = integrand_f[i_last:]
    seg_int = np.trapezoid(int_seg, r_seg)
    print(f"  Przedzial r=[{r_full[i_last]:.3f},{r_max:.1f}]: Ep_seg = {seg_int:.4e}")

total_ep = np.trapezoid(integrand_f, r_full)
print(f"  SUMA = {total_ep:.4e}")
print()
print(f"  psi_new ~ 1/sqrt(lam) = {1/np.sqrt(lam):.0f}")
print(f"  psi_core(M3) = {psi_f[0]:.1f}")
