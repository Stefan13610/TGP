"""
p4_phi0_limit.py

Problem P4: Czy Phi0^(2) -> 2*pi gdy a_gam -> 0?

Mechanizm hierarchii tla (badanie_pi.py): dla kazdego (alpha, a_gam)
wyznaczamy Phi0^(2) takie ze M(Phi0^(2)) / M(Phi0=1) = 207.

Obserwacja: dla alpha~5.9-6.0 i malych a_gam: Phi0^(2) ~ 2*pi (+1-3%).

Pytanie:
  a) Czy Phi0^(2) -> 2*pi gdy a_gam -> 0 (dla stalego alpha)?
  b) Jaki alpha* daje DOKLADNE Phi0^(2) = 2*pi (dla a_gam -> 0)?
  c) Czy stosunek Phi0^(3)/Phi0^(2) -> const gdy a_gam -> 0?

Metoda analityczna:
  W limicie a_gam -> 0, soliton Yukawa staje sie punktowy:
    K*exp(-msp*r)/r -> K*delta(r)/(4*pi*r^2) * 4*pi*r^2 = K*delta(r)
  Energia: E(K; Phi0) ~ K^2 * f(alpha, msp)
  Warunek samospojnosci: g(K) = E/K - 4*pi = 0
  => K * f(alpha, msp) = 4*pi

Dla r >> a_gam, integrand dazy do gładkiej funkcji r.
Limit a_gam -> 0 odpowiada usunieciu regularyzacji UV.
"""
import numpy as np
from scipy.optimize import brentq
import warnings; warnings.filterwarnings('ignore')

GAMMA = 1.0
PHI0_BASE = 1.0
LAM = 5.501357e-6

def V_mod_psi(psi, lam=LAM):
    return GAMMA/3*psi**3 - GAMMA/4*psi**4 + lam/6*(psi-1)**6

def energy_soliton(K, Phi0, alpha, a_gam, N=2000):
    """Energia solitonu Yukawa w tle Phi0."""
    msp = np.sqrt(max(GAMMA / Phi0, 1e-9))
    r_max = max(60.0/msp, 20.0)
    t   = np.linspace(0, 1, N)
    r   = a_gam * (r_max/a_gam)**t
    phi  = np.maximum(Phi0 + K*np.exp(-msp*r)/r, 1e-10)
    dphi = K*np.exp(-msp*r)*(-msp*r - 1.0)/r**2
    psi  = phi/Phi0
    V1   = V_mod_psi(1.0)
    Ek   = 4*np.pi*np.trapezoid(0.5*dphi**2*(1.0 + alpha/(Phi0*psi))*r**2, r)
    Ep   = 4*np.pi*Phi0**2*np.trapezoid((V_mod_psi(psi)-V1)*r**2, r)
    return Ek + Ep

def g_func(M, Phi0, alpha, a_gam):
    if M <= 1e-12: return 0.0
    K = M/(4*np.pi*Phi0)
    return energy_soliton(K, Phi0, alpha, a_gam) - M

def find_M1(Phi0, alpha, a_gam, M_lo=1e-5, M_hi=1000.0, N=300):
    """Znajdz pierwsze zero g(M) (masa M1 = najmniejsza generacja)."""
    M_arr = np.concatenate([
        np.linspace(M_lo, 5.0, N//3),
        np.linspace(5.0, 100.0, N//3),
        np.linspace(100.0, M_hi, N//3)])
    g_arr = [g_func(M, Phi0, alpha, a_gam) for M in M_arr]
    for i in range(len(g_arr)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            try:
                return brentq(lambda M: g_func(M, Phi0, alpha, a_gam),
                              M_arr[i], M_arr[i+1], xtol=1e-8)
            except Exception:
                pass
    return None

def find_Phi0_for_ratio(target_ratio, M1, alpha, a_gam,
                         Phi0_lo=1.01, Phi0_hi=100.0):
    """Znajdz Phi0 takie ze M1(Phi0)/M1 = target_ratio."""
    def h(Phi0):
        M = find_M1(Phi0, alpha, a_gam)
        if M is None: return target_ratio
        return M/M1 - target_ratio
    try:
        return brentq(h, Phi0_lo, Phi0_hi, xtol=1e-5, maxiter=40)
    except Exception:
        return None

# ================================================================
print("=" * 70)
print("P4: LIMIT a_gam -> 0: Czy Phi0^(2) -> 2*pi?")
print(f"   2*pi = {2*np.pi:.6f}")
print()

# ================================================================
# CZESC 1: Limit a_gam -> 0 dla roznych alpha
# ================================================================
print("CZESC 1: Phi0^(2) vs a_gam (limit a_gam -> 0)")
print()

alpha_vals = [5.5, 5.9, 6.0, 6.5, 7.0]
a_gam_vals = [0.10, 0.08, 0.05, 0.03, 0.02, 0.015, 0.010]

print(f"{'alpha':>6} {'a_gam':>7} {'Phi0_2':>10} {'P2/2pi':>10} {'diff%':>8} {'Phi0_3':>10} {'P3/P2':>8}")
print("-" * 68)

limit_data = {}
for alpha in alpha_vals:
    for a_gam in a_gam_vals:
        M1 = find_M1(1.0, alpha, a_gam)
        if M1 is None: continue
        P2 = find_Phi0_for_ratio(207.0, M1, alpha, a_gam)
        if P2 is None: continue
        P3 = find_Phi0_for_ratio(3477.0, M1, alpha, a_gam,
                                   Phi0_lo=P2+0.01, Phi0_hi=500.0)
        r2pi = P2/(2*np.pi)
        diff = 100*(P2 - 2*np.pi)/(2*np.pi)
        s3 = f"{P3:.5f}" if P3 else "---"
        r32 = f"{P3/P2:.4f}" if P3 else "---"
        print(f"{alpha:>6.1f} {a_gam:>7.3f} {P2:>10.5f} {r2pi:>10.6f} {diff:>8.3f}% {s3:>10} {r32:>8}")
        if alpha not in limit_data:
            limit_data[alpha] = []
        limit_data[alpha].append({'a_gam': a_gam, 'P2': P2, 'P3': P3})
    print()

# ================================================================
# CZESC 2: Interpolacja do a_gam = 0
# ================================================================
print()
print("CZESC 2: Ekstrapolacja Phi0^(2) do a_gam=0 (fit liniowy w a_gam)")
print()
print(f"{'alpha':>6} {'P2(a_gam=0)':>13} {'P2(0)/2pi':>12} {'diff%':>8}")
print("-" * 45)

for alpha in alpha_vals:
    if alpha not in limit_data or len(limit_data[alpha]) < 3:
        continue
    data = limit_data[alpha]
    # Bierzemy 3 najmniejsze a_gam
    data_sorted = sorted(data, key=lambda x: x['a_gam'])[:4]
    ag = np.array([d['a_gam'] for d in data_sorted])
    P2 = np.array([d['P2'] for d in data_sorted])
    # Fit: P2 = P2_0 + c1*a_gam + c2*a_gam^2
    if len(ag) >= 3:
        coeffs = np.polyfit(ag, P2, 2)
        P2_limit = coeffs[2]  # wartosc przy a_gam=0
        r2pi = P2_limit / (2*np.pi)
        diff = 100*(P2_limit - 2*np.pi)/(2*np.pi)
        print(f"{alpha:>6.1f} {P2_limit:>13.6f} {r2pi:>12.7f} {diff:>8.4f}%")

# ================================================================
# CZESC 3: Szukamy alpha* takie ze P2(a_gam->0) = DOKLADNIE 2*pi
# ================================================================
print()
print("CZESC 3: Szukanie alpha* takie ze Phi0^(2)(a_gam->0) = 2*pi DOKLADNIE")
print()
# Uzyjemy malego a_gam = 0.010 jako proxy dla limitu

a_gam_proxy = 0.010
print(f"Uzywamy a_gam={a_gam_proxy} jako przyblizenia limitu a_gam->0")
print()

def P2_at_proxy(alpha, a_gam=a_gam_proxy):
    M1 = find_M1(1.0, alpha, a_gam)
    if M1 is None: return None
    P2 = find_Phi0_for_ratio(207.0, M1, alpha, a_gam)
    return P2

# Skan alpha wokol 5.0-7.0
print(f"{'alpha':>8} {'P2':>12} {'diff_od_2pi%':>14}")
for al in np.linspace(4.5, 7.0, 11):
    P2v = P2_at_proxy(al)
    if P2v:
        diff = 100*(P2v - 2*np.pi)/(2*np.pi)
        print(f"{al:>8.3f} {P2v:>12.5f} {diff:>14.4f}%")

# Bisekcja: znajdz alpha* gdzie P2 = 2*pi
def residual_2pi(alpha):
    P2v = P2_at_proxy(alpha, a_gam=a_gam_proxy)
    if P2v is None: return 1.0
    return P2v - 2*np.pi

# Sprawdz czy jest zmiana znaku
r_lo = residual_2pi(5.0)
r_hi = residual_2pi(7.0)
print()
print(f"P2(5.0) - 2pi = {r_lo:+.5f}")
print(f"P2(7.0) - 2pi = {r_hi:+.5f}")
if r_lo * r_hi < 0:
    try:
        alpha_star = brentq(residual_2pi, 5.0, 7.0, xtol=0.005, maxiter=20)
        P2_star = P2_at_proxy(alpha_star)
        print(f"\n==> alpha* = {alpha_star:.4f}")
        print(f"    Phi0^(2)(alpha*) = {P2_star:.6f} = {P2_star/(2*np.pi):.6f} * 2*pi")
        print(f"    Blad od 2*pi: {100*(P2_star-2*np.pi)/(2*np.pi):+.4f}%")
    except Exception as ex:
        print(f"Bisekcja nie zbiegla: {ex}")
else:
    print("Brak zmiany znaku — P2 nie przechodzi przez 2*pi w [5,7]")

print()
print("=" * 70)
print("CZESC 4: Stosunek Phi0^(3)/Phi0^(2) w limicie a_gam->0")
print()
for alpha in [5.9, 6.0, 7.0]:
    if alpha not in limit_data: continue
    data = sorted(limit_data[alpha], key=lambda x: x['a_gam'])
    for d in data[:3]:
        if d['P3']:
            print(f"  alpha={alpha:.1f}, a_gam={d['a_gam']:.3f}: "
                  f"P3/P2 = {d['P3']/d['P2']:.5f}  "
                  f"P3/P2 vs pi = {d['P3']/d['P2']/np.pi:.5f}  "
                  f"vs 4 = {d['P3']/d['P2']/4:.5f}")

print()
print("  Hipotezy dla P3/P2:")
print(f"  pi = {np.pi:.5f}")
print(f"  4  = 4.00000")
print(f"  e  = {np.e:.5f}")
print(f"  2pi-2 = {2*np.pi-2:.5f}")
