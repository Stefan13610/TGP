"""
p13b_flux_scaling.py
====================
Test hipotezy K/a_gam^2 = const (stala gestosc strumienia) dla K2 i K3.

OBSERWACJA z p13:
  K2=2.033 z a_Gam=0.040: phi'(a_Gam) = -K2/a_Gam^2 = -1271 (absurdalnie strome)
  -> ZERO skonczonych profili (100% crash do N0 lub eksplozja)

HIPOTEZA: K/a_Gam^2 = const = K1*/a_Gam1^2
  K1* = 0.010414, a_Gam1 = 0.040 -> K/a^2 = 0.010414/0.0016 = 6.509
  a_Gam^(n) = a_Gam1 * sqrt(K_n / K1*)

Oznacza to skalowanie a_Gam proportional sqrt(K) (N IE 1/sqrt(K) jak sugeruje p10!).
Sprawdzamy czy K2/K3 sa samospojne przy tej hipotezie.
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
R_MAX  = 80.0
V1     = GAMMA/3 - GAMMA/4

def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi)/kfac + ALPHA*dphi**2/(2.0*phi**2*kfac) - (2.0/r)*dphi)
    return [dphi, ddphi]

def phi_at_rmax(psi_core, K, a_gam, r_max=R_MAX, n_eval=2000):
    dphi0 = -K / a_gam**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e6
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def compute_energy(psi_core, K, a_gam, n_eval=5000):
    dphi0 = -K / a_gam**2
    r_eval = a_gam * (R_MAX/a_gam)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [a_gam, R_MAX], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_eval)
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return Ek+Ep, Ek, Ep

def best_branch(K, a_gam, psi_max=50.0, n_psi=200, label=''):
    psi_vals = np.concatenate([
        np.linspace(1.001, 5.0, 120),
        np.linspace(5.0, psi_max, 80)
    ])
    psi_vals = np.unique(psi_vals)
    F_vals = np.array([phi_at_rmax(p, K, a_gam) - 1.0 for p in psi_vals])
    n_fin = int(np.sum(np.isfinite(F_vals)))

    best_g    = np.nan
    best_psi  = np.nan
    all_zeros = []

    for i in range(len(psi_vals)-1):
        Fi, Fj = F_vals[i], F_vals[i+1]
        if np.isfinite(Fi) and np.isfinite(Fj) and Fi*Fj < 0:
            try:
                psi_z = brentq(
                    lambda p: phi_at_rmax(p, K, a_gam) - 1.0,
                    psi_vals[i], psi_vals[i+1],
                    xtol=1e-4, maxiter=30
                )
                E, Ek, Ep = compute_energy(psi_z, K, a_gam)
                g = E / (4*np.pi*K) - 1.0
                all_zeros.append({'psi': psi_z, 'g': g, 'E': E})
                if np.isnan(best_g) or abs(g) < abs(best_g):
                    best_g   = g
                    best_psi = psi_z
            except Exception:
                pass

    return best_psi, best_g, n_fin, all_zeros, psi_vals, F_vals

# ============================================================
# HIPOTEZA K/a_gam^2 = const
# ============================================================
K1_star = 0.010414
a_gam1  = 0.040
flux_density = K1_star / a_gam1**2

print("P13b: Test hipotezy K/a_Gam^2 = const")
print("="*65)
print(f"K1* = {K1_star}, a_Gam1 = {a_gam1}")
print(f"K1*/a_Gam1^2 = {flux_density:.4f}")
print()

K2 = 2.032728
K3 = 34.14450

a_gam2_flux = np.sqrt(K2 / flux_density)
a_gam3_flux = np.sqrt(K3 / flux_density)

print(f"Hipoteza K/a^2 = {flux_density:.4f}:")
print(f"  K2={K2:.4f}: a_Gam2 = sqrt({K2:.4f}/{flux_density:.4f}) = {a_gam2_flux:.4f}")
print(f"  K3={K3:.4f}: a_Gam3 = sqrt({K3:.4f}/{flux_density:.4f}) = {a_gam3_flux:.4f}")
print(f"  Skalowanie: a_Gam(n) = a_Gam1 * sqrt(K_n/K1*)")
print(f"    a_Gam2/a_Gam1 = {a_gam2_flux/a_gam1:.3f}")
print(f"    a_Gam3/a_Gam1 = {a_gam3_flux/a_gam1:.3f}")
print()
print(f"  phi'(a_Gam) = -K/a_Gam^2:")
print(f"    K1* at a1=0.040: phi' = {-K1_star/a_gam1**2:.4f}")
print(f"    K2  at a2={a_gam2_flux:.4f}: phi' = {-K2/a_gam2_flux**2:.4f}  [zamiast -1271]")
print(f"    K3  at a3={a_gam3_flux:.4f}: phi' = {-K3/a_gam3_flux**2:.4f}  [zamiast -21340]")
print()

# ============================================================
# TEST: K2 z a_gam2_flux
# ============================================================
print(f"Test K2={K2:.4f} z a_Gam2={a_gam2_flux:.4f}:")
psi2, g2, n2, zeros2, psi_arr2, F_arr2 = best_branch(K2, a_gam2_flux, psi_max=30.0, label='K2_flux')
print(f"  Profile skonczonych: {n2}/200")
print(f"  Galezi (zer F): {len(zeros2)}")
for z in zeros2:
    ok = 'SAMOSPOJNY' if abs(z['g'])<0.05 else ('~bliski' if abs(z['g'])<0.5 else 'niespojny')
    print(f"    psi_core={z['psi']:.4f}, g={z['g']:.4e}, E={z['E']:.4e}  {ok}")
print()

# ============================================================
# TEST: K3 z a_gam3_flux
# ============================================================
print(f"Test K3={K3:.4f} z a_Gam3={a_gam3_flux:.4f}:")
psi3, g3, n3, zeros3, psi_arr3, F_arr3 = best_branch(K3, a_gam3_flux, psi_max=100.0, label='K3_flux')
print(f"  Profile skonczonych: {n3}/200")
print(f"  Galezi (zer F): {len(zeros3)}")
for z in zeros3:
    ok = 'SAMOSPOJNY' if abs(z['g'])<0.05 else ('~bliski' if abs(z['g'])<0.5 else 'niespojny')
    print(f"    psi_core={z['psi']:.4f}, g={z['g']:.4e}, E={z['E']:.4e}  {ok}")
print()

# ============================================================
# WYKRESY
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
fig.suptitle(f'P13b: Hipoteza K/a_Gam^2=const={flux_density:.3f}\n'
             f'(alpha={ALPHA}, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

def plot_profile_ax(ax, psi_core, K, a_gam, label, color):
    dphi0 = -K/a_gam**2
    r_eval = a_gam*(R_MAX/a_gam)**np.linspace(0, 1, 3000)
    sol = solve_ivp(ode_rhs, [a_gam, R_MAX], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_eval)
    if len(sol.t) < 5: return
    r, phi = sol.t, sol.y[0]
    ax.plot(r, phi, '-', color=color, lw=1.5, label=label)
    ax.axhline(1.0, color='gray', lw=0.7, linestyle='--', alpha=0.7)
    ax.set_xlabel('r'); ax.set_ylabel('phi(r)')
    ax.grid(True, alpha=0.3)

# K2 flux -- F(psi)
ax = axes[0, 0]
mask2 = np.isfinite(F_arr2)
ax.plot(psi_arr2[mask2], np.clip(F_arr2[mask2], -5, 20), 'b-', lw=1.2)
ax.axhline(0, color='red', lw=1, linestyle='--')
for z in zeros2:
    ax.axvline(z['psi'], color='green', lw=1.2, alpha=0.7, label=f"g={z['g']:.2f}")
ax.set_title(f'K2={K2:.3f}, a_Gam2={a_gam2_flux:.4f}\nF(psi_core)', fontsize=9)
ax.set_xlabel('psi_core'); ax.set_ylabel('F'); ax.legend(fontsize=7)
ax.grid(True, alpha=0.3); ax.set_ylim(-5, 15)

# K2 flux -- profil
ax = axes[0, 1]
if not np.isnan(psi2):
    plot_profile_ax(ax, psi2, K2, a_gam2_flux, f'K2 psi={psi2:.3f}', 'orange')
ax.set_title(f'K2={K2:.3f}, a2={a_gam2_flux:.4f}, g={g2:.4e}', fontsize=9)
ax.set_xlim(0, 60)

# K2 flux vs K1 (porownanie)
ax = axes[0, 2]
plot_profile_ax(ax, 1.2419, K1_star, a_gam1, f'K1* a1=0.040', 'green')
if not np.isnan(psi2):
    plot_profile_ax(ax, psi2, K2, a_gam2_flux, f'K2 a2={a_gam2_flux:.3f}', 'orange')
ax.set_title('Porownanie K1* vs K2', fontsize=9)
ax.set_xlim(0, 60); ax.legend(fontsize=7)

# K3 flux -- F(psi)
ax = axes[1, 0]
mask3 = np.isfinite(F_arr3)
ax.plot(psi_arr3[mask3], np.clip(F_arr3[mask3], -5, 20), 'b-', lw=1.2)
ax.axhline(0, color='red', lw=1, linestyle='--')
for z in zeros3:
    ax.axvline(z['psi'], color='green', lw=1.2, alpha=0.7, label=f"g={z['g']:.2f}")
ax.set_title(f'K3={K3:.3f}, a_Gam3={a_gam3_flux:.4f}\nF(psi_core)', fontsize=9)
ax.set_xlabel('psi_core'); ax.set_ylabel('F'); ax.legend(fontsize=7)
ax.grid(True, alpha=0.3); ax.set_ylim(-5, 15)

# K3 flux -- profil
ax = axes[1, 1]
if not np.isnan(psi3):
    plot_profile_ax(ax, psi3, K3, a_gam3_flux, f'K3 psi={psi3:.3f}', 'red')
ax.set_title(f'K3={K3:.3f}, a3={a_gam3_flux:.4f}, g={g3:.4e}', fontsize=9)
ax.set_xlim(0, 70)

# Skalowanie a_Gam
ax = axes[1, 2]
Kvals = np.logspace(-3, 2, 100)
a_flux = np.sqrt(Kvals / flux_density)
ax.loglog(Kvals, a_flux, 'b-', lw=2, label='a_Gam prop sqrt(K) [K/a^2=const]')
ax.loglog([K1_star, K2, K3], [a_gam1, a_gam2_flux, a_gam3_flux],
          'ko', ms=10, zorder=5, label='Gen 1,2,3')
ax.axhline(a_gam1, color='gray', lw=0.8, linestyle=':', alpha=0.7)
ax.set_xlabel('K'); ax.set_ylabel('a_Gam')
ax.set_title('Skalowanie a_Gam(K) -- hipoteza K/a^2=const', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

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
print("PODSUMOWANIE p13b:")
print()
print(f"Hipoteza K/a^2 = {flux_density:.4f} (rowna poczatkowa gestosc strumienia):")
print(f"  K2={K2:.4f}: a_Gam2={a_gam2_flux:.4f}, n_fin={n2}, g_best={g2:.4e}")
print(f"  K3={K3:.4f}: a_Gam3={a_gam3_flux:.4f}, n_fin={n3}, g_best={g3:.4e}")
print()
if abs(g2) < 0.1 and abs(g3) < 0.1:
    print("SUKCES: Oba solitony samospojne przy K/a^2 = const!")
    print("  -> Fizyczna zasada: stala gestosc strumienia na granicy a_Gam")
elif abs(g2) < 1.0 or abs(g3) < 1.0:
    print("CZESCIOWY SUKCES: Przynajmniej jeden bliski samospojnosci")
    print("  -> Hipoteza czesciowo potwierdzona, wymaga fine-tuningu")
else:
    print("WYNIK NEGATYWNY: K/a^2 = const NIE daje samospojnosci K2/K3")
    print("  -> Hipoteza obalona, trzeba szukac innego skalowania")
