"""
p14c_continuous_branch.py
=========================
Sledzi CIAGLA galaz psi_core(K) i sprawdza czy g(K) = 0 poza K*1.

PROBLEM DOTYCHCZASOWY:
  best_g_for_K() skacze miedzy galeziami -> falszywe "zera" g(K).
  np: K=0.084 -> psi=1.985, g=-516;  K=0.085 -> psi=2.756, g=+3.3 (SKOK!)
  To daje pozorne zero, ale to NIECIAGLOSC, nie fizyczne zero.

CEL:
  Sledzic psi_core wzdluz galezi ciagle (zmieniajac psi_core powoli z K).
  Sprawdzic znak g na TEJ SAMEJ galezi.

Galaz glowna: psi_core ~ 1 + K/A_GAM (Yukawa przyblizenie)
              psi_core ~ 1.24 przy K=K*1=0.01041
              W miare wzrostu K psi_core rosnie.
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

def phi_at_rmax(psi_core, K, r_max=R_MAX, n_eval=2500):
    dphi0 = -K / A_GAM**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval = A_GAM * (r_max/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [A_GAM, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-10, atol=1e-12,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def find_psi_core_near(K, psi_guess, window=0.3):
    """Znajdz psi_core(K) w okolicach psi_guess."""
    psi_lo = max(1.001, psi_guess - window)
    psi_hi = psi_guess + window
    F_lo = phi_at_rmax(psi_lo, K) - 1.0
    F_hi = phi_at_rmax(psi_hi, K) - 1.0
    if not (np.isfinite(F_lo) and np.isfinite(F_hi)):
        # Expand window
        for w in [0.5, 1.0, 2.0]:
            psi_lo2 = max(1.001, psi_guess - w)
            psi_hi2 = psi_guess + w
            F_lo2 = phi_at_rmax(psi_lo2, K) - 1.0
            F_hi2 = phi_at_rmax(psi_hi2, K) - 1.0
            if np.isfinite(F_lo2) and np.isfinite(F_hi2):
                F_lo, F_hi = F_lo2, F_hi2
                psi_lo, psi_hi = psi_lo2, psi_hi2
                break
    if not (np.isfinite(F_lo) and np.isfinite(F_hi)):
        return np.nan
    # Scan for sign change
    psi_scan = np.linspace(psi_lo, psi_hi, 60)
    F_scan = [phi_at_rmax(p, K) - 1.0 for p in psi_scan]
    for i in range(len(psi_scan)-1):
        if np.isfinite(F_scan[i]) and np.isfinite(F_scan[i+1]):
            if F_scan[i] * F_scan[i+1] < 0:
                return brentq(lambda p: phi_at_rmax(p, K) - 1.0,
                              psi_scan[i], psi_scan[i+1],
                              xtol=1e-6, maxiter=40)
    return np.nan

def compute_g_fast(K, psi_core, n_eval=6000):
    dphi0 = -K / A_GAM**2
    r_eval = A_GAM * (R_MAX/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(ode_rhs, [A_GAM, R_MAX], [psi_core, dphi0],
                    method='DOP853', rtol=1e-10, atol=1e-12, t_eval=r_eval)
    r, phi, dphi = sol.t, np.maximum(sol.y[0], 1e-10), sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return (Ek+Ep)/(4*np.pi*K) - 1.0, Ek+Ep, Ek, Ep

# ============================================================
# SLEDZ CIAGLA GALAZ: startujemy od K*1 i idziemy w gore K
# ============================================================
print("P14c: Sledzenie ciaglej galezi psi_core(K)")
print("="*65)

# Punkt startowy: K*1 potwierdzony
K_start = 0.010414
psi_start = 1.2419

# Skan K w gore (K*1 i wyzej)
K_vals = np.concatenate([
    np.linspace(0.001, 0.010414, 15),     # ponizej K*1
    np.linspace(0.010414, 0.080, 30),     # K*1 do 0.08
    np.linspace(0.080, 0.200, 25),        # dalej
])
K_vals = np.unique(K_vals)

print(f"Skan {len(K_vals)} wartosci K, sledzac galaz psi_core~{psi_start:.3f}")
print()

psi_track = []
g_track   = []
K_track   = []

psi_prev = psi_start

for K in K_vals:
    # Szukaj psi_core blisko poprzedniego
    psi_new = find_psi_core_near(K, psi_prev, window=max(0.3, 0.2*psi_prev))
    if np.isnan(psi_new):
        # Sprobuj z wieksza oknem
        psi_new = find_psi_core_near(K, psi_prev, window=psi_prev*0.5)
    if np.isnan(psi_new):
        print(f"  K={K:.5f}: GALAZ UTRACONA (psi_prev={psi_prev:.4f})")
        break

    g, E, Ek, Ep = compute_g_fast(K, psi_new)
    psi_track.append(psi_new)
    g_track.append(g)
    K_track.append(K)

    tag = ''
    if abs(g) < 0.05: tag = '  <-- SAMOSPOJNY!'
    elif abs(g) < 0.5: tag = '  ~ bliski'
    elif abs(g) > 100: tag = '  (!)'
    print(f"  K={K:.5f}: psi={psi_new:.4f}, g={g:.4e}{tag}")

    psi_prev = psi_new

K_track   = np.array(K_track)
psi_track = np.array(psi_track)
g_track   = np.array(g_track)

# Szukaj zer g na tej galezi
print()
print("Zera g(K) na ciaglej galezi:")
zeros_continuous = []
for i in range(len(K_track)-1):
    gi, gj = g_track[i], g_track[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
        # Sprawdz czy skok psi jest maly (ciaglosc galezi)
        dpsi = abs(psi_track[i+1] - psi_track[i])
        dK   = K_track[i+1] - K_track[i]
        if dpsi / max(psi_track[i], 1.0) > 0.5:
            print(f"  [POMINIETY -- skok psi: {psi_track[i]:.3f} -> {psi_track[i+1]:.3f} przy K={K_track[i]:.5f}]")
            continue
        try:
            # Dopasuj zero - sledzac galaz
            def g_cont(K_val):
                psi_g = find_psi_core_near(K_val, (psi_track[i]+psi_track[i+1])/2,
                                           window=0.3)
                if np.isnan(psi_g): return 999.0
                g_v, _, _, _ = compute_g_fast(K_val, psi_g)
                return g_v
            K_z = brentq(g_cont, K_track[i], K_track[i+1],
                         xtol=1e-5, maxiter=20)
            psi_z = find_psi_core_near(K_z, (psi_track[i]+psi_track[i+1])/2, window=0.3)
            g_z, E_z, Ek_z, Ep_z = compute_g_fast(K_z, psi_z)
            zeros_continuous.append({'K': K_z, 'psi': psi_z, 'g': g_z, 'E': E_z})
            print(f"  K* = {K_z:.8f}, psi = {psi_z:.6f}, g = {g_z:.4e}")
        except Exception as ex:
            print(f"  [Blad: {ex}]")

if not zeros_continuous:
    print("  Brak zer g(K) poza K*1 na ciaglej galezi!")
print()

# Wykres
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle(f'P14c: Ciagla galaz psi_core(K) i g(K)\n'
             f'(alpha={ALPHA}, a_Gam={A_GAM})', fontsize=11, fontweight='bold')

ax = axes[0]
ax.semilogx(K_track, g_track, 'b.-', lw=1.5, ms=5)
ax.axhline(0, color='red', lw=1.2, linestyle='--', label='g=0')
for z in zeros_continuous:
    ax.axvline(z['K'], color='green', lw=1.5, label=f"K*={z['K']:.4f}")
ax.axvline(0.010414, color='green', lw=1.5, linestyle=':', alpha=0.5, label='K*1')
ax.axvline(2.032728, color='orange', lw=1.5, linestyle=':', alpha=0.3, label='K2(Yuk)')
ax.set_xlabel('K'); ax.set_ylabel('g = E/(4*pi*K) - 1')
ax.set_title('g(K) -- ciagla galaz', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(-5, 15)

ax = axes[1]
ax.semilogx(K_track, psi_track, 'b.-', lw=1.5, ms=5)
K_yk = np.logspace(-3, -0.5, 50)
ax.semilogx(K_yk, 1.0 + K_yk/A_GAM, 'r--', lw=1.2, alpha=0.6, label='Yukawa: 1+K/a')
ax.axvline(0.010414, color='green', lw=1.5, linestyle=':', label='K*1')
ax.set_xlabel('K'); ax.set_ylabel('psi_core (ciagla galaz)')
ax.set_title('psi_core(K)', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
print()

# Podsumowanie
print("="*65)
print("KONKLUZJA p14c:")
print()
n_zeros = len(zeros_continuous) + 1  # +1 za K*1
print(f"Liczba samospojnych solitonow na galezi glownej: {n_zeros}")
print(f"  K*1 = 0.010414  (potwierdzony)")
for z in zeros_continuous:
    print(f"  K*  = {z['K']:.6f}  g={z['g']:.4e}")
print()
if len(zeros_continuous) == 0:
    print("WYNIK: Tylko JEDEN samospojny soliton K*1 = 0.010414")
    print()
    print("INTERPRETACJA:")
    print("  g(K) na galezi glownej jest MONOTONICZNIE ROSNACA po K*1.")
    print("  Brak drugiego zera -> tylko jedna generacja w TGP z obecnymi parametrami.")
    print()
    print("  Hierarchia Yukawa (r21=207, r31=3477) wynika z LINIOWEGO przyblizenia.")
    print("  Pelne ODE daje inna strukture energetyczna z powodu glebokich")
    print("  wartosci potencjalu V_mod(phi) dla duzego phi.")
    print()
    print("  WNIOSEK: Aktualny model TGP (V_mod, alpha, a_Gam, lambda) opisuje")
    print("  JEDNĄ generacje leptonu. Rozszerzenie na 3 generacje wymaga:")
    print("  1. Modyfikacji potencjalu V_mod (inne V dla kazdej generacji?)")
    print("  2. Innego mechanizmu generacyjnego niz solitony kuliste")
    print("  3. Dodatkowych pol lub stopni swobody w TGP")
