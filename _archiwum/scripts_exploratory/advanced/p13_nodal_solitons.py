"""
p13_nodal_solitons.py
=====================
Solitony wezlowe (nodal / excited-state solitons) w TGP.

MOTYWACJA:
  K1* = 0.010414 jest samospojny na glownej galezi ODE (psi_core~1.24, brak wezlow).
  K2 = 2.033 i K3 = 34.14 nie maja profilu phi(inf)=1 na tej galezi przy a_Gam=0.040.
  Hipoteza bańkowa (N0 w centrum) odrzucona przez p12.

NOWA HIPOTEZA: SOLITONY WEZLOWE (EXCITED STATES)
  W mechanice kwantowej stany wzbudzone maja wezly (miejsca, gdzie funkcja falowa
  przekracza wartosc rownowagowa i wraca). Analogicznie w TGP:

  - Gen 1 (K1*):  profil monotoniczny phi(r): z psi_core do 1 (brak wezlow)
  - Gen 2 (K2):   profil z 1 wezlem: phi nurkuje ponizej 1 przy r~r_node, wraca do 1
  - Gen 3 (K3):   profil z 2 wezlami

  Klucz: wezly sa w zakresie phi > 0 (nie dotykaja N0),
  wiec ODE pozostaje regularne (brak osobliwosci).

CEL p13:
  1. Skan psi_core dla K2 i K3 z szerszym oknem (psi_core moze byc ujemne lub < 1)
  2. Szukanie profili z phi(R_max) = 1 ALE z jednym lub wiecej "nurkowania" ponizej 1
  3. Policzenie wezlow i energie dla kazdej galezi
  4. Weryfikacja samospojnosci g(K) = E/(4*pi*K) - 1 = 0 na galezi wezlowej
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ============================================================
# PARAMETRY
# ============================================================
ALPHA  = 8.5616
A_GAM  = 0.040
LAM    = 5.501357e-06
GAMMA  = 1.0
M_EFF  = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX  = 80.0   # zwiekszony zakres dla wezlowych profili

print(f"Parametry: alpha={ALPHA}, a_Gam={A_GAM}, lam={LAM:.4e}")
print(f"m_eff = {M_EFF:.4f},  R_max = {R_MAX}")
print()

V1 = GAMMA/3 - GAMMA/4

def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi) / kfac
             + ALPHA * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]

def integrate_outward(psi_core, K, r_max=R_MAX, n_eval=3000):
    """Outward integration; returns full solution object."""
    phi0  = psi_core
    dphi0 = -K / A_GAM**2

    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e6
    ev_hi.terminal = True; ev_hi.direction = 1

    r_eval = A_GAM * (r_max / A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(
        ode_rhs, [A_GAM, r_max], [phi0, dphi0],
        method='DOP853', rtol=1e-9, atol=1e-11,
        t_eval=r_eval, events=[ev_lo, ev_hi]
    )
    return sol

def phi_at_rmax(psi_core, K, r_max=R_MAX, n_eval=2000):
    """phi(R_max) dla danego psi_core, K."""
    sol = integrate_outward(psi_core, K, r_max, n_eval)
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def count_nodes(sol, threshold=1.0):
    """Policz ile razy phi przekracza 'threshold' od gory (nurkowania ponizej)."""
    if len(sol.t) < 10:
        return 0
    phi = sol.y[0]
    crossings = 0
    for i in range(len(phi)-1):
        if phi[i] >= threshold and phi[i+1] < threshold:
            crossings += 1
    return crossings

def profile_energy(psi_core, K, n_eval=5000):
    """Energia solitonu."""
    dphi0 = -K / A_GAM**2
    r_eval = A_GAM * (R_MAX / A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(
        ode_rhs, [A_GAM, R_MAX], [psi_core, dphi0],
        method='DOP853', rtol=1e-9, atol=1e-11,
        t_eval=r_eval
    )
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return Ek+Ep, Ek, Ep

# ============================================================
# KROK 1: Szeroki skan psi_core dla K2 z rozszerzona skala
# ============================================================
print("="*65)
print("KROK 1: Skan F(psi_core) dla K2=2.033 -- profil + wezly")
print("="*65)

K2 = 2.032728

# Szerokie okno: moze byc male psi_core (np. 0.1) dla glebokich wezlow?
# Ale phi musi pozostac > 0, wiec minimalne psi_core jest ograniczone przez
# to by profil nie rozbil sie przez N0 zanim dotrze do R_max.
# Probujemy: psi_core od 1.001 (standard) az do 200 (duze)
# ORAZ psi_core < 1: profil startuje ponizej minimum -- profil wezlowy?

psi_ranges = [
    np.linspace(0.01, 1.0, 80),    # profil startuje ponizej 1
    np.linspace(1.001, 10.0, 150), # standard
    np.linspace(10.0, 100.0, 80),  # duze psi_core
]
psi_vals_K2 = np.unique(np.concatenate(psi_ranges))

print(f"  Skanowanie {len(psi_vals_K2)} wartosci psi_core dla K2")

F_K2 = []
nodes_K2 = []
for psi in psi_vals_K2:
    sol = integrate_outward(psi, K2)
    if len(sol.t) == 0 or sol.t[-1] < R_MAX * 0.99:
        F_K2.append(np.nan)
        nodes_K2.append(-1)
    else:
        F_K2.append(sol.y[0, -1] - 1.0)
        nodes_K2.append(count_nodes(sol, threshold=1.0))

F_K2 = np.array(F_K2)
nodes_K2 = np.array(nodes_K2)

# Znajdz zera dla kazdej liczby wezlow
print()
print(f"  Dystrybucja wezlow wsrod skonczonych profili:")
fin_mask = np.isfinite(F_K2)
for n_node in sorted(set(nodes_K2[fin_mask])):
    cnt = np.sum(nodes_K2[fin_mask] == n_node)
    print(f"    {n_node} wezlow: {cnt} profili")

print()
print(f"  Szukanie zer F(psi_core) per liczba wezlow:")

zeros_K2_by_nodes = {}
for n_node in sorted(set(nodes_K2[fin_mask])):
    mask_n = (nodes_K2 == n_node) & fin_mask
    zeros_n = []
    psi_n = psi_vals_K2[mask_n]
    F_n   = F_K2[mask_n]
    for i in range(len(psi_n)-1):
        if F_n[i]*F_n[i+1] < 0:
            try:
                psi_z = brentq(
                    lambda p: phi_at_rmax(p, K2) - 1.0,
                    psi_n[i], psi_n[i+1],
                    xtol=1e-4, maxiter=30
                )
                E, Ek, Ep = profile_energy(psi_z, K2)
                g = E / (4*np.pi*K2) - 1.0
                zeros_n.append({'psi_core': psi_z, 'E': E, 'g': g, 'n_nodes': n_node})
            except Exception:
                pass
    zeros_K2_by_nodes[n_node] = zeros_n
    if zeros_n:
        print(f"  [{n_node} wezlow]:")
        for z in zeros_n:
            tag = 'SAMOSPOJNY' if abs(z['g']) < 0.05 else ('~bliski' if abs(z['g']) < 0.5 else 'niespojny')
            print(f"    psi_core={z['psi_core']:.4f}, E={z['E']:.4e}, g={z['g']:.4e}  {tag}")

print()

# ============================================================
# KROK 2: To samo dla K3=34.14
# ============================================================
print("="*65)
print("KROK 2: Skan F(psi_core) dla K3=34.14 -- profil + wezly")
print("="*65)

K3 = 34.14450

psi_ranges_K3 = [
    np.linspace(0.01, 1.0, 60),
    np.linspace(1.001, 50.0, 150),
    np.linspace(50.0, 500.0, 80),
]
psi_vals_K3 = np.unique(np.concatenate(psi_ranges_K3))
print(f"  Skanowanie {len(psi_vals_K3)} wartosci psi_core dla K3")

F_K3 = []
nodes_K3 = []
for psi in psi_vals_K3:
    sol = integrate_outward(psi, K3)
    if len(sol.t) == 0 or sol.t[-1] < R_MAX * 0.99:
        F_K3.append(np.nan)
        nodes_K3.append(-1)
    else:
        F_K3.append(sol.y[0, -1] - 1.0)
        nodes_K3.append(count_nodes(sol, threshold=1.0))

F_K3 = np.array(F_K3)
nodes_K3 = np.array(nodes_K3)

fin_mask3 = np.isfinite(F_K3)
print()
print(f"  Dystrybucja wezlow wsrod skonczonych profili:")
for n_node in sorted(set(nodes_K3[fin_mask3])):
    cnt = np.sum(nodes_K3[fin_mask3] == n_node)
    print(f"    {n_node} wezlow: {cnt} profili")

print()
print(f"  Szukanie zer F(psi_core) per liczba wezlow:")

zeros_K3_by_nodes = {}
for n_node in sorted(set(nodes_K3[fin_mask3])):
    mask_n = (nodes_K3 == n_node) & fin_mask3
    zeros_n = []
    psi_n = psi_vals_K3[mask_n]
    F_n   = F_K3[mask_n]
    for i in range(len(psi_n)-1):
        if F_n[i]*F_n[i+1] < 0:
            try:
                psi_z = brentq(
                    lambda p: phi_at_rmax(p, K3) - 1.0,
                    psi_n[i], psi_n[i+1],
                    xtol=1e-4, maxiter=30
                )
                E, Ek, Ep = profile_energy(psi_z, K3)
                g = E / (4*np.pi*K3) - 1.0
                zeros_n.append({'psi_core': psi_z, 'E': E, 'g': g, 'n_nodes': n_node})
            except Exception:
                pass
    zeros_K3_by_nodes[n_node] = zeros_n
    if zeros_n:
        print(f"  [{n_node} wezlow]:")
        for z in zeros_n:
            tag = 'SAMOSPOJNY' if abs(z['g']) < 0.05 else ('~bliski' if abs(z['g']) < 0.5 else 'niespojny')
            print(f"    psi_core={z['psi_core']:.4f}, E={z['E']:.4e}, g={z['g']:.4e}  {tag}")

print()

# ============================================================
# KROK 3: Profile typowych wezlowych solitonow
# ============================================================
print("="*65)
print("KROK 3: Rysowanie profili wezlowych")
print("="*65)

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle('P13: Solitony wezlowe w TGP\n'
             f'(alpha={ALPHA}, a_Gam={A_GAM}, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

def plot_profile(ax, psi_core, K, label, color, style='-'):
    sol = integrate_outward(psi_core, K, r_max=R_MAX, n_eval=4000)
    if len(sol.t) == 0:
        ax.set_title(f'{label}\n[brak profilu]', fontsize=8)
        return
    r   = sol.t
    phi = sol.y[0]
    n_nodes = count_nodes(sol)
    ax.plot(r, phi, style, color=color, lw=1.5, label=f'psi={psi_core:.2f}, {n_nodes}w')
    ax.axhline(1.0, color='gray', lw=0.8, linestyle='--', alpha=0.7)
    ax.axhline(0.0, color='black', lw=0.5, linestyle=':', alpha=0.4)
    ax.set_xlabel('r')
    ax.set_ylabel('phi(r)')
    ax.set_title(f'{label}', fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 40)
    ax.set_ylim(-0.5, max(psi_core * 1.1, 2.0))

# K*_1 (referencja -- 0 wezlow)
ax = axes[0, 0]
plot_profile(ax, 1.2419, 0.010414, 'K*_1=0.0104 (ref: 0 wezlow)', 'green')

# K2 -- rozne galezi
ax = axes[0, 1]
ax.set_title('K2=2.033 -- F(psi_core) overview', fontsize=9)
mask_fin2 = np.isfinite(F_K2)
for n_node in [0, 1, 2]:
    mask_n = (nodes_K2 == n_node) & mask_fin2
    if mask_n.any():
        F_clip = np.clip(F_K2[mask_n], -5, 20)
        ax.plot(psi_vals_K2[mask_n], F_clip, '.', ms=2,
                label=f'{n_node} wezlow', alpha=0.6)
ax.axhline(0, color='red', lw=1, linestyle='--')
ax.set_xlabel('psi_core')
ax.set_ylabel('F = phi(R_max)-1')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)
ax.set_ylim(-5, 15)

# K2 -- najlepszy profil (min |g|) z kazdej galezi wezlowej
ax = axes[0, 2]
ax.set_title('K2=2.033 -- najlepsze galezi wezlowe', fontsize=9)
for n_node, zeros_n in zeros_K2_by_nodes.items():
    if zeros_n:
        best = min(zeros_n, key=lambda z: abs(z['g']))
        sol = integrate_outward(best['psi_core'], K2, r_max=R_MAX, n_eval=3000)
        if len(sol.t) > 10 and sol.t[-1] >= R_MAX*0.99:
            r, phi = sol.t, sol.y[0]
            ax.plot(r, phi, lw=1.2, label=f'{n_node}w, psi={best["psi_core"]:.1f}, g={best["g"]:.2f}')
ax.axhline(1.0, color='gray', lw=0.8, linestyle='--')
ax.set_xlabel('r'); ax.set_ylabel('phi(r)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
ax.set_xlim(0, 50); ax.set_ylim(-0.5, 3)

# K3 -- F(psi_core)
ax = axes[1, 0]
ax.set_title('K3=34.14 -- F(psi_core) overview', fontsize=9)
mask_fin3 = np.isfinite(F_K3)
for n_node in [0, 1, 2]:
    mask_n = (nodes_K3 == n_node) & mask_fin3
    if mask_n.any():
        F_clip = np.clip(F_K3[mask_n], -5, 20)
        ax.plot(psi_vals_K3[mask_n], F_clip, '.', ms=2,
                label=f'{n_node} wezlow', alpha=0.6)
ax.axhline(0, color='red', lw=1, linestyle='--')
ax.set_xlabel('psi_core'); ax.set_ylabel('F = phi(R_max)-1')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(-5, 15)

# K3 -- najlepszy profil z kazdej galezi
ax = axes[1, 1]
ax.set_title('K3=34.14 -- najlepsze galezi wezlowe', fontsize=9)
for n_node, zeros_n in zeros_K3_by_nodes.items():
    if zeros_n:
        best = min(zeros_n, key=lambda z: abs(z['g']))
        sol = integrate_outward(best['psi_core'], K3, r_max=R_MAX, n_eval=3000)
        if len(sol.t) > 10 and sol.t[-1] >= R_MAX*0.99:
            r, phi = sol.t, sol.y[0]
            ax.plot(r, phi, lw=1.2, label=f'{n_node}w, psi={best["psi_core"]:.1f}, g={best["g"]:.2f}')
ax.axhline(1.0, color='gray', lw=0.8, linestyle='--')
ax.set_xlabel('r'); ax.set_ylabel('phi(r)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
ax.set_xlim(0, 60); ax.set_ylim(-0.5, 5)

# Podsumowanie g dla wszystkich galezi
ax = axes[1, 2]
ax.set_title('g (samospojnosc) dla galezi wezlowych', fontsize=9)
K2_gs = [(z['n_nodes'], z['g']) for zeros in zeros_K2_by_nodes.values() for z in zeros]
K3_gs = [(z['n_nodes'], z['g']) for zeros in zeros_K3_by_nodes.values() for z in zeros]
if K2_gs:
    nn2 = [x[0] for x in K2_gs]
    g2  = [x[1] for x in K2_gs]
    ax.scatter(nn2, g2, c='orange', s=60, label='K2=2.033', zorder=5)
if K3_gs:
    nn3 = [x[0] for x in K3_gs]
    g3  = [x[1] for x in K3_gs]
    ax.scatter(nn3, g3, c='red', s=60, label='K3=34.14', zorder=5)
ax.axhline(0, color='black', lw=1, linestyle='--', alpha=0.7, label='g=0 (samospojne)')
ax.set_xlabel('Liczba wezlow'); ax.set_ylabel('g = E/(4*pi*K) - 1')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"  Zapisano: {out}")
print()

# ============================================================
# PODSUMOWANIE
# ============================================================
print("="*65)
print("PODSUMOWANIE p13: SOLITONY WEZLOWE")
print("="*65)
print()
print("K2 = 2.033:")
if any(zeros_K2_by_nodes.values()):
    for n_node, zeros_n in sorted(zeros_K2_by_nodes.items()):
        for z in zeros_n:
            ok = 'SAMOSPOJNY' if abs(z['g']) < 0.05 else ('~bliski' if abs(z['g']) < 1 else 'niespojny')
            print(f"  {n_node} wezlow: psi_core={z['psi_core']:.3f}, g={z['g']:.4e}  {ok}")
else:
    print("  Brak samospojnych galezi wezlowych")

print()
print("K3 = 34.14:")
if any(zeros_K3_by_nodes.values()):
    for n_node, zeros_n in sorted(zeros_K3_by_nodes.items()):
        for z in zeros_n:
            ok = 'SAMOSPOJNY' if abs(z['g']) < 0.05 else ('~bliski' if abs(z['g']) < 1 else 'niespojny')
            print(f"  {n_node} wezlow: psi_core={z['psi_core']:.3f}, g={z['g']:.4e}  {ok}")
else:
    print("  Brak samospojnych galezi wezlowych")

print()
print("INTERPRETACJA:")
print("  Jesli g=0 na galezi n-wezlowej dla K_n:")
print("    -> K_n to stan wzbudzony solitonu TGP (analogia do QM)")
print("    -> n=0: elektron (K1*), n=1: mion (K2?), n=2: taon (K3?)")
print("  Jesli brak g=0 dla wezlowych galezi:")
print("    -> Mechanizm wezlowy NIE jest odpowiedzialny za hierarchie")
print("    -> Wymaga nowego podejscia (lambda*(n), inna topologia)")
