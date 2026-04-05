"""
p16_schwarzschild_selfconsist.py
=================================
Samospojnosc przy SKALOWANIU a_Gam = C*K (hipoteza Schwarzschilda).

SYNTEZA p10-p15:
  - p14c: przy stalym a_Gam=0.040 -> JEDEN soliton K*1=0.010414
  - Yukawa: 3 generacje wymagaja profili bliskich phi=1 (liniowa approx.)
  - Warunek phi~1: potrzebne K/a_Gam << 1, tj. a_Gam >> K

  KLUCZ: jesli a_Gam = C*K (proporcjonalne do masy):
    psi_core ~ 1 + (K/a_Gam)*exp(-m*a_Gam) = 1 + (1/C)*exp(-m*C*K)
    -> Dla duzego K: psi_core -> 1 (profil jest bliski 1!)
    -> Aproks. Yukawa POPRAWIA sie dla ciezszych generacji

  Pytanie: czy istnieje stala C taka ze wszystkie 3 K_Yukawa
           sa samospojne przy a_Gam = C*K?

CEL p16:
  1. Dla stalej C: wyznacz g(K) przy a_Gam(K) = C*K
  2. Znajdz C i K* takie ze g(K*)=0 dla 3 K wartosci (Yukawa K1,K2,K3)
  3. Sprawdz czy sa to te same K co Yukawa lub bliskie
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

ALPHA  = 8.5616
LAM    = 5.501357e-06
GAMMA  = 1.0
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX_BASE = 80.0
V1     = GAMMA/3 - GAMMA/4

print("P16: Samospojnosc przy a_Gam = C*K (skalowanie Schwarzschilda)")
print("="*65)
print(f"m_eff = {M_EFF:.4f},  1/m_eff = {1/M_EFF:.4f}")
print()

def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi)/kfac + ALPHA*dphi**2/(2.0*phi**2*kfac) - (2.0/r)*dphi)
    return [dphi, ddphi]

def phi_at_rmax(psi_core, K, a_gam, r_max, n_eval=2500):
    dphi0 = -K / a_gam**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    # Geometryczna siatka, clip zeby uniknac wartosci poza [a_gam, r_max]
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam * (1 + 1e-12), r_max * (1 - 1e-12))
    sol = solve_ivp(ode_rhs, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def find_psi_and_g(K, a_gam, psi_guess=None, n_eval=5000):
    """Znajdz psi_core i oblicz g dla danego (K, a_gam)."""
    if psi_guess is None:
        psi_guess = 1.0 + K / a_gam * np.exp(-M_EFF * a_gam)

    r_max = max(R_MAX_BASE, 30.0 / M_EFF, 20.0 * a_gam)

    # Skan psi_core wokol guess
    w = max(0.5, psi_guess * 0.3)
    psi_lo = max(1.001, psi_guess - w)
    psi_hi = psi_guess + w

    psi_z = np.nan
    for window in [w, 2*w, 5*w, 10*w]:
        psi_lo2 = max(1.001, psi_guess - window)
        psi_hi2 = psi_guess + window
        psi_scan = np.linspace(psi_lo2, psi_hi2, 50)
        F_scan = [phi_at_rmax(p, K, a_gam, r_max) - 1.0 for p in psi_scan]
        for i in range(len(psi_scan)-1):
            if np.isfinite(F_scan[i]) and np.isfinite(F_scan[i+1]):
                if F_scan[i] * F_scan[i+1] < 0:
                    try:
                        psi_z = brentq(
                            lambda p: phi_at_rmax(p, K, a_gam, r_max) - 1.0,
                            psi_scan[i], psi_scan[i+1],
                            xtol=1e-6, maxiter=30
                        )
                        break
                    except Exception:
                        pass
        if not np.isnan(psi_z):
            break

    if np.isnan(psi_z):
        return np.nan, np.nan

    # Oblicz energie
    dphi0 = -K / a_gam**2
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam*(1+1e-12), r_max*(1-1e-12))
    sol = solve_ivp(ode_rhs, [a_gam, r_max], [psi_z, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_eval)
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    g = (Ek+Ep)/(4*np.pi*K) - 1.0
    return psi_z, g

# ============================================================
# KROK 1: Wyznacz C z warunku K*1 = 0.010414 przy a_Gam = 0.040
# ============================================================
K1_star = 0.010414
A_GAM1  = 0.040
C_measured = A_GAM1 / K1_star
print(f"Stala C z K*1: C = a_Gam1/K*1 = {A_GAM1}/{K1_star} = {C_measured:.4f}")
print()

# ============================================================
# KROK 2: g(K) przy a_Gam = C*K dla C=C_measured
# ============================================================
print(f"KROK 2: g(K) przy a_Gam = {C_measured:.4f}*K:")
print()

K_vals = np.logspace(-3, 2, 40)

g_Sch   = []
psi_Sch = []

psi_prev = 1.24

for K in K_vals:
    a_gam = C_measured * K
    psi_Yukawa = 1.0 + (K/a_gam) * np.exp(-M_EFF * a_gam)
    psi_z, g = find_psi_and_g(K, a_gam, psi_guess=psi_Yukawa)
    g_Sch.append(g)
    psi_Sch.append(psi_z)

    tag = ''
    if not np.isnan(g):
        if abs(g) < 0.05: tag = ' <-- SAMOSPOJNY!'
        elif abs(g) < 0.5: tag = ' ~ bliski'
    print(f"  K={K:.4e}: a_Gam={a_gam:.4f}, psi_Yuk={psi_Yukawa:.4f}, "
          f"psi_ODE={psi_z:.4f}, g={g:.4e}{tag}")

g_Sch   = np.array(g_Sch)
psi_Sch = np.array(psi_Sch)

# Znajdz zera g
print()
print("Zera g(K):")
zeros_Sch = []
for i in range(len(K_vals)-1):
    gi, gj = g_Sch[i], g_Sch[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
        # Sprawdz ciaglosc
        dpsi = abs(psi_Sch[i+1] - psi_Sch[i]) if (not np.isnan(psi_Sch[i]) and not np.isnan(psi_Sch[i+1])) else 999
        if dpsi > 1.5:
            print(f"  [Skok psi: {psi_Sch[i]:.3f} -> {psi_Sch[i+1]:.3f} -- prawdopodobnie galaz zmiana]")
            continue
        K_approx = (K_vals[i]*K_vals[i+1])**0.5
        print(f"  K* ~ {K_approx:.4e}  (g: {gi:.4e} -> {gj:.4e})")
        try:
            def g_func(K_v):
                a_g = C_measured * K_v
                psi_Y = 1.0 + (K_v/a_g)*np.exp(-M_EFF*a_g)
                psi_z_v, g_v = find_psi_and_g(K_v, a_g, psi_guess=psi_Y)
                return g_v if not np.isnan(g_v) else 999.0
            K_zero = brentq(g_func, K_vals[i], K_vals[i+1],
                            xtol=1e-5, maxiter=20)
            a_gam_z = C_measured * K_zero
            psi_Y_z = 1.0 + (K_zero/a_gam_z)*np.exp(-M_EFF*a_gam_z)
            psi_z_z, g_z = find_psi_and_g(K_zero, a_gam_z, psi_guess=psi_Y_z)
            zeros_Sch.append({'K': K_zero, 'psi': psi_z_z, 'g': g_z, 'a_gam': a_gam_z})
            print(f"    -> K* = {K_zero:.6f}, a_Gam = {a_gam_z:.4f}, psi_core = {psi_z_z:.4f}, g = {g_z:.4e}")
        except Exception as ex:
            print(f"    [Blad: {ex}]")

print()
print(f"Lacznie {len(zeros_Sch)} samospojnych solitonow przy a_Gam = {C_measured:.4f}*K")
if zeros_Sch:
    print()
    print("Porownanie z Yukawa:")
    K_Yuk = [0.009820, 2.032728, 34.14450]
    for z in zeros_Sch:
        print(f"  K*(ODE) = {z['K']:.6f}  a_Gam = {z['a_gam']:.4f}")
    print(f"  K_Yukawa(1,2,3) = {K_Yuk}")

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(f'P16: Samospojnosc przy a_Gam = C*K (C={C_measured:.3f})\n'
             f'(alpha={ALPHA}, lam={LAM:.2e})', fontsize=11, fontweight='bold')

# g(K) przy a_Gam=C*K
ax = axes[0]
mask = np.isfinite(g_Sch)
ax.semilogx(K_vals[mask], np.clip(g_Sch[mask], -10, 20), 'b.-', lw=1.5, ms=5)
ax.axhline(0, color='red', lw=1.2, linestyle='--', label='g=0')
for z in zeros_Sch:
    ax.axvline(z['K'], color='green', lw=1.5, alpha=0.8, label=f"K*={z['K']:.3e}")
ax.axvline(0.009820, color='orange', lw=1.5, linestyle=':', alpha=0.6, label='K1_Yuk')
ax.axvline(2.032728, color='red', lw=1.5, linestyle=':', alpha=0.4, label='K2_Yuk')
ax.axvline(34.14450, color='purple', lw=1.5, linestyle=':', alpha=0.4, label='K3_Yuk')
ax.set_xlabel('K'); ax.set_ylabel('g = E/(4*pi*K) - 1')
ax.set_title(f'g(K) przy a_Gam = {C_measured:.3f}*K', fontsize=9)
ax.legend(fontsize=6); ax.grid(True, alpha=0.3); ax.set_ylim(-8, 10)

# psi_core(K)
ax = axes[1]
ax.semilogx(K_vals[mask], psi_Sch[mask], 'b.-', lw=1.5, ms=5, label='psi ODE')
K_line = np.logspace(-3, 2, 100)
psi_Yukawa_line = 1.0 + (1/C_measured) * np.exp(-M_EFF * C_measured * K_line)
ax.semilogx(K_line, psi_Yukawa_line, 'r--', lw=1.2, alpha=0.6, label='Yukawa approx')
ax.set_xlabel('K'); ax.set_ylabel('psi_core')
ax.set_title('psi_core(K) -- blisko 1?', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.axhline(1.5, color='gray', lw=0.7, linestyle=':', alpha=0.5)

# a_Gam skalowanie
ax = axes[2]
K_line2 = np.logspace(-3, 2, 100)
ax.loglog(K_line2, C_measured * K_line2, 'b-', lw=2, label=f'a_Gam={C_measured:.2f}*K')
ax.axvline(K1_star, color='green', lw=1.5, linestyle=':', label='K*1')
if zeros_Sch:
    for z in zeros_Sch:
        ax.scatter([z['K']], [z['a_gam']], s=100, c='green', zorder=5)
ax.set_xlabel('K'); ax.set_ylabel('a_Gam')
ax.set_title('Skalowanie Schwarzschilda', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
