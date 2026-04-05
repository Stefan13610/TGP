"""
p1_ode_shooting.py

Problem P1: Weryfikacja zasadnosci ansatzu Yukawa dla M3.

Pelne rownanie ODE dla profilu radialnego Phi(r) z potencjalem V_mod:
  (1 + alpha/phi) phi'' + (2/r)(1+alpha/phi)phi'
  - (alpha/(2*phi^2))(phi')^2 - V'_mod(phi) = 0

Granice:
  phi(r->inf) = Phi0 = 1  (proznia)
  phi'(0) = 0             (sferyczna symetria)

Metoda:
  1. Warunki poczatkowe z ansatzu Yukawa phi_Y(r)=1+K*exp(-r)/r
     na r = a_gam (unikamy osobliwosci r=0)
  2. Integracja ODE od r=a_gam do r=r_max = 60
  3. Porownanie profilu ODE z ansatzem Yukawa:
     - ksztalt phi(r) vs phi_Y(r)
     - energia E_ODE vs E_Yukawa
     - samospojnosc: E_ODE / K_eff - 4*pi

Parametry: alfa=8.5616, a_gam=0.040, lam*=5.501e-6, K3=34.1445 (najlepsze rozwiazanie).
Testujemy rowniez K1 i K2.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import warnings; warnings.filterwarnings('ignore')

# ---- Parametry najlepszego rozwiazania ----
ALPHA = 8.5616
A_GAM = 0.040
LAM   = 5.501357e-06
K1    = 0.009820
K2    = 2.032728
K3    = 34.14450
PHI0  = 1.0

GAMMA = 1.0

def V_mod(psi, lam=LAM):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def dV_mod(psi, lam=LAM):
    return GAMMA*psi**2 - GAMMA*psi**3 + lam*(psi-1)**5

V1 = V_mod(1.0)

# ---- Pelne ODE ----
def ode_rhs(r, y, alpha=ALPHA, lam=LAM):
    """
    y = [phi, phi']
    phi'' = V'_mod(phi)/(1+a/phi) + alpha*(phi')^2/(2*phi^2*(1+a/phi)) - (2/r)*phi'
    """
    phi, dphi = y
    phi = max(phi, 1e-10)
    a = 1.0 + alpha/phi
    ddphi = (dV_mod(phi, lam)/a
             + alpha*dphi**2/(2*phi**2*a)
             - (2.0/r)*dphi)
    return [dphi, ddphi]

# ---- Ansatz Yukawa ----
def yukawa(r, K):
    return 1.0 + K*np.exp(-r)/r

def yukawa_d(r, K):
    return K*np.exp(-r)*(-r-1.0)/r**2

# ---- Energia na siatce log ----
def energy_log(K, N=4000):
    t   = np.linspace(0, 1, N)
    r   = A_GAM * (60.0/A_GAM)**t
    phi = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K*np.exp(-r)*(-r - 1.0)/r**2
    psi  = phi
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep   = 4*np.pi*np.trapezoid((V_mod(psi)-V1)*r**2, r)
    return Ek, Ep, Ek+Ep

# ---- Integracja ODE od r=a_gam z warunkami poczatkowymi z ansatzu ----
def integrate_ode(K, r_max=60.0, n_eval=5000):
    r0 = A_GAM
    phi0 = yukawa(r0, K)
    dphi0 = yukawa_d(r0, K)
    y0 = [phi0, dphi0]

    # gesta siatka dla dokladnosci
    r_eval = A_GAM * (r_max/A_GAM)**np.linspace(0, 1, n_eval)

    sol = solve_ivp(
        ode_rhs, [r0, r_max], y0,
        t_eval=r_eval, method='DOP853',
        rtol=1e-10, atol=1e-12,
        dense_output=False
    )
    return sol

# ---- Energia profilu ODE ----
def energy_ode(sol):
    r   = sol.t
    phi = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return Ek, Ep, Ek+Ep

# ---- Dopasowanie K_eff z ODE: phi_ODE(r) ~ 1 + K_eff * exp(-r)/r ----
def fit_K_eff(sol):
    """Dopasuj K_eff przez minimalizacje bledu Yukawa dla r w [1, 10]."""
    r   = sol.t
    phi = sol.y[0]
    # Maska: r in [0.5, 15] dla dopasowania Yukawa
    mask = (r > 0.5) & (r < 15.0)
    r_fit  = r[mask]
    phi_fit = phi[mask]
    # Liniowe dopasowanie: phi-1 = K_eff * exp(-r)/r
    # => (phi-1)*r*exp(r) = K_eff
    K_vals = (phi_fit - 1.0) * r_fit * np.exp(r_fit)
    K_eff = np.mean(K_vals)
    K_std = np.std(K_vals)
    return K_eff, K_std

# ================================================================
# GLOWNA ANALIZA
# ================================================================
print("=" * 75)
print("P1: WERYFIKACJA ANSATZU YUKAWA — pelne ODE vs ansatz")
print(f"Parametry: alpha={ALPHA}, a_gam={A_GAM}, lam={LAM:.4e}")
print(f"Cel: 4*pi = {4*np.pi:.6f}")
print()

for K_label, K in [("K1", K1), ("K2", K2), ("K3", K3)]:
    print(f"{'='*75}")
    print(f"Kandydat: {K_label} = {K:.5f}")
    print()

    # 1. Ansatz Yukawa
    Ek_Y, Ep_Y, E_Y = energy_log(K, N=8000)
    g_Y = E_Y/K - 4*np.pi
    psi_core = yukawa(A_GAM, K)

    print(f"  Ansatz Yukawa (siatka log, N=8000):")
    print(f"    Ek={Ek_Y:.4e}, Ep={Ep_Y:.4e}, E={E_Y:.4e}")
    print(f"    g(K) = E/K - 4*pi = {g_Y:.4e}")
    print(f"    psi_core = phi({A_GAM:.3f}) = {psi_core:.2f}")
    print()

    # 2. Pelne ODE
    print(f"  Integracja pelnego ODE (DOP853, rtol=1e-10)...")
    sol = integrate_ode(K, r_max=60.0, n_eval=6000)

    phi_final = sol.y[0, -1]
    dphi_final = sol.y[1, -1]
    print(f"    phi(r_max=60) = {phi_final:.8f}  (cel: 1.0)")
    print(f"    phi'(r_max)   = {dphi_final:.4e}  (cel: 0.0)")
    print()

    # Energia ODE
    Ek_O, Ep_O, E_O = energy_ode(sol)
    g_O = E_O/K - 4*np.pi
    print(f"  Energia ODE:")
    print(f"    Ek={Ek_O:.4e}, Ep={Ep_O:.4e}, E={E_O:.4e}")
    print(f"    g(K) = E_ODE/K - 4*pi = {g_O:.4e}")
    print()

    # K_eff z dopasowania
    K_eff, K_std = fit_K_eff(sol)
    print(f"  Dopasowanie K_eff z profilu ODE (r in [0.5, 15]):")
    print(f"    K_eff = {K_eff:.5f} (wejsciowe K = {K:.5f})")
    print(f"    std[(phi-1)*r*e^r] = {K_std:.4e}")
    print(f"    blad K_eff: {100*(K_eff-K)/K:+.4f}%")
    print()

    # Maksymalna roznica profili
    r   = sol.t
    phi_ODE = sol.y[0]
    phi_Y_vec = yukawa(r, K)
    diff = phi_ODE - phi_Y_vec
    rel_diff = np.abs(diff / np.maximum(phi_Y_vec - 1.0, 1e-10))
    # Pomijamy r bliskie r_max gdzie phi~1 i roznica jest numeryczna
    mask = r < 30.0
    print(f"  Roznica |phi_ODE - phi_Yukawa| / |phi_Y - 1|:")
    print(f"    max (r < 30)   = {rel_diff[mask].max():.4e}")
    print(f"    mean (r < 30)  = {rel_diff[mask].mean():.4e}")
    print(f"    at r=a_gam     = {abs(diff[0])/(psi_core-1):.4e}")
    print()

print("=" * 75)

# ================================================================
# SZCZEGOLOWY TEST DLA K3: profil krok po kroku
# ================================================================
print()
print("=" * 75)
print(f"SZCZEGOLOWY PROFIL dla K3={K3:.4f}")
print(f"{'r':>8} {'phi_ODE':>12} {'phi_Yukawa':>12} {'roznica':>12} {'rel_diff%':>12}")
print("-" * 65)

sol3 = integrate_ode(K3, r_max=60.0, n_eval=6000)
r_pts = [A_GAM, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0]
for r_pt in r_pts:
    # Interpolacja do punktu r_pt
    idx = np.searchsorted(sol3.t, r_pt)
    idx = min(idx, len(sol3.t)-1)
    phi_o = sol3.y[0, idx]
    phi_y = yukawa(sol3.t[idx], K3)
    d = phi_o - phi_y
    rd = 100*abs(d)/max(abs(phi_y-1), 1e-8)
    print(f"{sol3.t[idx]:>8.4f} {phi_o:>12.6f} {phi_y:>12.6f} {d:>12.4e} {rd:>12.4f}%")

print()
print("=" * 75)
print()

# ================================================================
# TEST ZBIEZNOSCI: czy zmiana poczatkowych warunkow (K) zmienia profil?
# ================================================================
print("WRAZLIWOSC NA K: czy K3 jest izolowanym zerem ODE?")
print()
print(f"{'K':>10} {'phi(60)':>12} {'g_ODE':>14} {'K_eff':>10} {'blad_K%':>10}")
print("-" * 60)

for dK in [-1.0, -0.5, 0.0, +0.5, +1.0]:
    K_test = K3 + dK
    sol_t = integrate_ode(K_test, r_max=60.0, n_eval=3000)
    phi_fin = sol_t.y[0, -1]
    _, _, E_t = energy_ode(sol_t)
    g_t = E_t/K_test - 4*np.pi
    K_eff_t, _ = fit_K_eff(sol_t)
    print(f"{K_test:>10.4f} {phi_fin:>12.6f} {g_t:>14.4e} {K_eff_t:>10.4f} {100*(K_eff_t-K_test)/K_test:>+10.4f}%")

print()
print("=" * 75)
print("WNIOSKI:")
print("  1. Jesli phi_ODE(60) ~ 1.0 dla K3: ansatz Yukawa jest dobrym IC")
print("  2. Jesli g_ODE(K3) ~ 0: samospojnosc potwierdzona dla pelnego ODE")
print("  3. Jesli max(rel_diff) < 5%: ansatz Yukawa jest dobrym przyblizeniem profilu")
print("  4. Jesli K_eff ~ K przy malym std: profil jest faktycznie Yukawa-podobny")
