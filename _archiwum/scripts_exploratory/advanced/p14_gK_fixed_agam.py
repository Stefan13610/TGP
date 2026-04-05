"""
p14_gK_fixed_agam.py
====================
Ile samospojnych solitonow istnieje przy stalym a_Gam=0.040?

PYTANIE FUNDAMENTALNE:
  Czy g(K) przy a_Gam=0.040 ma 1 zero (tylko K1*) czy 3 zera (jak w Yukawa)?

PROBLEM z poprzednimi skanami:
  - p6: skan K in [0.003, 0.1] na GLOWNEJ galezi psi_core~1.2 -> jedno zero K1*=0.010414
  - K2=2.033 i K3=34.14: ZERO skonczonych profili przy a_Gam=0.040 (nadmierny nachylenie)

PLAN p14:
  1. Skan K in [0.0001, 100] przy a_Gam=0.040
     - Dla kazdego K: szerokie okno psi_core (w tym duze wartosci!)
     - Szukaj WSZYSTKICH galezi i ich g
  2. Narysuj g_best(K) -- krzywa samospojnosci
  3. Policzyc zera g(K) = 0 -- tyle jest fizycznych generacji w TGP ODE

UWAGA na potencjal:
  V_mod(phi) = phi^3/3 - phi^4/4 + lam*(phi-1)^6/6
  Dla duzego phi: V_mod(phi) ~ -phi^4/4 -> -inf
  Wiec dla duzego psi_core: Ep ~ V_mod(psi_core)*V_core ~ -phi^4*r^3 -> glebokie ujemne!
  To moze sprawiac ze E < 0 dla K > K_threshold
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
A_GAM  = 0.040
LAM    = 5.501357e-06
GAMMA  = 1.0
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX  = 80.0
V1     = GAMMA/3 - GAMMA/4

print("P14: Mapa g(K) przy stalym a_Gam=0.040")
print("="*65)
print(f"Parametry: alpha={ALPHA}, a_Gam={A_GAM}, lam={LAM:.4e}")
print(f"V1 = {V1:.6f},  m_eff = {M_EFF:.4f}")
print()

def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi)/kfac + ALPHA*dphi**2/(2.0*phi**2*kfac) - (2.0/r)*dphi)
    return [dphi, ddphi]

def phi_at_rmax(psi_core, K, r_max=R_MAX, n_eval=2000):
    dphi0 = -K / A_GAM**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval = A_GAM * (r_max/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [A_GAM, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def compute_energy_full(psi_core, K, n_eval=5000):
    dphi0 = -K / A_GAM**2
    r_eval = A_GAM * (R_MAX/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [A_GAM, R_MAX], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_eval)
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return Ek+Ep, Ek, Ep

def best_g_for_K(K, psi_max=5.0, n_psi=120):
    """Najlepsze g dla danego K (wszystkie galezi)."""
    # Adaptacyjny psi_max: dla duzego K, Yukawa psi_core ~ 1 + K/A_GAM = duze
    psi_yukawa = 1.0 + K / A_GAM   # psi Yukawa przyblizone
    psi_max_eff = max(psi_max, min(psi_yukawa * 1.5, 3000.0))

    # Wieloskalowe skanowanie
    psi_vals = np.unique(np.concatenate([
        np.linspace(1.001, 3.0, 60),
        np.linspace(3.0, min(20.0, psi_max_eff), 40),
        np.linspace(min(20.0, psi_max_eff), psi_max_eff, 30)
    ]))

    F_vals = np.array([phi_at_rmax(p, K) - 1.0 for p in psi_vals])

    best_g   = np.nan
    best_psi = np.nan
    n_fin    = int(np.sum(np.isfinite(F_vals)))

    for i in range(len(psi_vals)-1):
        Fi, Fj = F_vals[i], F_vals[i+1]
        if np.isfinite(Fi) and np.isfinite(Fj) and Fi*Fj < 0:
            try:
                psi_z = brentq(
                    lambda p: phi_at_rmax(p, K) - 1.0,
                    psi_vals[i], psi_vals[i+1],
                    xtol=1e-4, maxiter=25
                )
                E, Ek, Ep = compute_energy_full(psi_z, K)
                g = E / (4*np.pi*K) - 1.0
                if np.isnan(best_g) or abs(g) < abs(best_g):
                    best_g   = g
                    best_psi = psi_z
            except Exception:
                pass

    return best_psi, best_g, n_fin

# ============================================================
# SKAN g(K) dla K in [1e-4, 50], stale a_Gam = 0.040
# ============================================================
# Rozne gestosci probkowania
K_scan = np.unique(np.concatenate([
    np.logspace(-4, -2, 20),    # K: 0.0001 do 0.01 (K1* region)
    np.logspace(-2, -1, 15),    # K: 0.01 do 0.1
    np.logspace(-1,  1, 15),    # K: 0.1 do 10
    np.logspace(1,   2, 10),    # K: 10 do 100 (K3 region)
]))

print(f"Skan {len(K_scan)} wartosci K w [{K_scan[0]:.4e}, {K_scan[-1]:.1f}]")
print(f"a_Gam = {A_GAM} (stalal)")
print()

g_best_arr  = np.full(len(K_scan), np.nan)
psi_best_arr = np.full(len(K_scan), np.nan)
n_fin_arr   = np.zeros(len(K_scan), dtype=int)

for i, K in enumerate(K_scan):
    psi_z, g, n_fin = best_g_for_K(K)
    g_best_arr[i]   = g
    psi_best_arr[i] = psi_z
    n_fin_arr[i]    = n_fin
    tag = ''
    if not np.isnan(g):
        if abs(g) < 0.1: tag = '  <-- SAMOSPOJNY!'
        elif abs(g) < 1.0: tag = '  ~ bliski'
    print(f"  K={K:.5e}  n_fin={n_fin:3d}  psi*={psi_z:.4f}  g={g:.4e}{tag}")

# Znajdz zera g(K)
print()
print("="*65)
print("ZERA g(K) = 0:")
zeros_K = []
for i in range(len(K_scan)-1):
    gi, gj = g_best_arr[i], g_best_arr[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
        try:
            def g_of_K(K_val):
                _, g_val, _ = best_g_for_K(K_val)
                return g_val if not np.isnan(g_val) else 999.0
            K_zero = brentq(g_of_K, K_scan[i], K_scan[i+1],
                            xtol=1e-4*K_scan[i], maxiter=15)
            psi_z, g_z, _ = best_g_for_K(K_zero)
            E, Ek, Ep = compute_energy_full(psi_z, K_zero)
            zeros_K.append({'K': K_zero, 'psi': psi_z, 'g': g_z, 'E': E})
            print(f"  K* = {K_zero:.6f}, psi_core = {psi_z:.4f}, g = {g_z:.4e}")
        except Exception as ex:
            print(f"  [Nie udalo sie wyznaczac zero miedzy {K_scan[i]:.4e} i {K_scan[i+1]:.4e}: {ex}]")

print()
print(f"Lacznie {len(zeros_K)} samospojnych solitonow przy a_Gam=0.040")

# ============================================================
# WYKRESY
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle(f'P14: g(K) przy stalym a_Gam={A_GAM}\n'
             f'(alpha={ALPHA}, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

# g(K) -- pelny skan
ax = axes[0]
mask_fin = np.isfinite(g_best_arr)
ax.semilogx(K_scan[mask_fin], np.clip(g_best_arr[mask_fin], -20, 20), 'b.-', lw=1.5, ms=5)
ax.axhline(0, color='red', lw=1.2, linestyle='--', label='g=0 (samospojne)')
for z in zeros_K:
    ax.axvline(z['K'], color='green', lw=1.5, alpha=0.8,
               label=f"K*={z['K']:.4f}")
ax.axvline(0.010414, color='green', lw=1.5, linestyle=':', alpha=0.5, label='K1*(ODE)')
ax.axvline(2.032728, color='orange', lw=1.5, linestyle=':', alpha=0.5, label='K2(Yukawa)')
ax.axvline(34.1445, color='red', lw=1.5, linestyle=':', alpha=0.5, label='K3(Yukawa)')
ax.set_xlabel('K'); ax.set_ylabel('g = E/(4*pi*K) - 1')
ax.set_title('g(K) -- wszystkie galezi\na_Gam=0.040', fontsize=9)
ax.legend(fontsize=6); ax.grid(True, alpha=0.3); ax.set_ylim(-20, 10)

# psi_core(K)
ax = axes[1]
ax.semilogx(K_scan[mask_fin], psi_best_arr[mask_fin], 'b.-', lw=1.5, ms=5)
ax.axvline(0.010414, color='green', lw=1.5, linestyle=':', label='K1*')
# Linia Yukawa: psi_core ~ 1 + K/a_Gam
K_yk = np.logspace(-4, 2, 100)
ax.semilogx(K_yk, 1.0 + K_yk/A_GAM, 'r--', lw=1.2, alpha=0.6, label='Yukawa: 1+K/a')
ax.set_xlabel('K'); ax.set_ylabel('psi_core (najlepsza galaz)')
ax.set_title('psi_core(K)', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# n_fin(K) -- ile profili dochodzi do R_max
ax = axes[2]
ax.semilogx(K_scan, n_fin_arr, 'g.-', lw=1.5, ms=5)
ax.axvline(0.010414, color='green', lw=1.5, linestyle=':', label='K1*')
ax.axvline(2.032728, color='orange', lw=1.5, linestyle=':', label='K2')
ax.axvline(34.1445, color='red', lw=1.5, linestyle=':', label='K3')
ax.set_xlabel('K'); ax.set_ylabel('Liczba skonczonych profili / 130')
ax.set_title('Ile profili phi(R_max)=1 istnieje?', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
print()

# ============================================================
# PODSUMOWANIE
# ============================================================
print("="*65)
print("PODSUMOWANIE p14: MAPA g(K) przy stalym a_Gam=0.040")
print()
print(f"Liczba samospojnych solitonow (g=0): {len(zeros_K)}")
if len(zeros_K) == 1:
    print("  -> Tylko JEDEN samospojny soliton: K1* (elektron)")
    print("  -> K2 i K3 (Yukawa) NIE odpowiadaja samospojnym rozwiazaniom")
    print("     ODE z a_Gam=0.040 i stalym lambda!")
    print()
    print("  INTERPRETACJA:")
    print("  Modele Yukawa przesadzaja -- masy generacji 2 i 3 wymagaja")
    print("  albo innego a_Gam, albo innego V_mod, albo sa artefaktem")
    print("  liniowego przyblizenia.")
elif len(zeros_K) == 3:
    print("  -> TRZY samospojne solitony (Gen 1, 2, 3)!")
    for z in zeros_K:
        print(f"     K* = {z['K']:.5e}, psi_core = {z['psi']:.4f}")
else:
    print(f"  -> {len(zeros_K)} samospojnych solitonow:")
    for z in zeros_K:
        print(f"     K* = {z['K']:.5e}, psi_core = {z['psi']:.4f}, g = {z['g']:.4e}")
