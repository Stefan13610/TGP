"""
p15_Vmod_multiGen.py
====================
Szukanie modyfikacji V_mod dajacych 3 samospojne solitony w pelnym ODE.

WYNIK p14c: przy V_mod_current = phi^3/3 - phi^4/4 + lam*(phi-1)^6/6
  g(K) ma JEDNO zero na galezi glownej -> tylko JEDEN soliton K*1.

PRZYCZYNA: V_mod(phi) ~ -phi^4/4 dla duzego phi -> E < 0 dla K > K*1.

CEL p15: Znalezc parametry dodatkowego czlonu V_extra(phi) takiego by:
  g(K) mialo 3 zera (generacje 1, 2, 3).

KRYTERIUM:
  1. V_mod(1) = 1/12 (zachowanie prozni)
  2. V'_mod(1) = 0   (1 jest ekstremum)
  3. V''_mod(1) = -1 (masa skalaru = 1)
  4. V_mod(phi) ma 3 samospojne solitony g=0

BADANE MODYFIKACJE:
  A) V_barrier: dodaje bariere phi^n / (1 + beta*phi^m) -- ogranicza zanurzenie w ujemne wartosci
  B) V_oscillate: dodaje oscylacje A*sin(omega*(phi-1))*exp(-mu*(phi-1)^2)
  C) V_soft: zastepuje -phi^4/4 czlonem ogranicz. -phi^4/(4*(1 + eps*phi^2))

Wszystkie modyfikacje musza zachowac:
  - V_mod(1), V'_mod(1), V''_mod(1) niezmienione (leptony sa wciaz tymi samymi cialami)
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
GAMMA  = 1.0
R_MAX  = 80.0

print("P15: Szukanie modyfikacji V_mod dla 3 generacji")
print("="*65)

# ============================================================
# ANALIZA: DLACZEGO AKTUALNY V_mod daje tylko 1 soliton
# ============================================================
print("ANALIZA biezacego potencjalu:")
print()

LAM_BASE = 5.501357e-06
def V_current(p):
    return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM_BASE/6*(p-1)**6
def dV_current(p):
    return GAMMA*p**2 - GAMMA*p**3 + LAM_BASE*(p-1)**5

V1 = V_current(1.0)
print(f"  V_mod(1) = {V1:.6f} (oczekiwane: 1/12={1/12:.6f})")
print(f"  V'_mod(1) = {dV_current(1.0):.2e}")
d2V = 2*GAMMA - 3*GAMMA + 5*LAM_BASE*0
print(f"  V''_mod(1) = {d2V:.4f} (oczekiwane: -1)")
print()

phi_test = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0])
print(f"  V_mod(phi) - V1 dla roznych phi:")
for phi in phi_test:
    print(f"    phi={phi:5.1f}: V_mod={V_current(phi):.4e}, V-V1={V_current(phi)-V1:.4e}")
print()
print("  -> V_mod(phi) ~ -phi^4/4 dla duzego phi: GLEBOKIE NEGATYWNE")
print("  -> Ep = 4*pi*int(V-V1)*r^2 dr << 0 dla duzego psi_core")
print("  -> g = E/(4*pi*K) - 1 ~ Ep/(4*pi*K) - 1 < 0 dla K > K*1")
print()

# ============================================================
# MODYFIKACJA C: V_soft -- ograniczenie ujemnej czesci
# ============================================================
# V_soft(phi) = phi^3/3 - phi^4/(4*(1 + eps*phi^2)) + lam*(phi-1)^6/6 + correction
# Chcemy zachowac V(1)=1/12, V'(1)=0, V''(1)=-1
# Correction term: musi skasowac zmiane w phi=1

print("MODYFIKACJA C: V_soft (saturacja -phi^4/4 czlonu):")
print()

def V_soft_base(p, eps):
    """Soft saturation: -phi^4/4*(1+eps*phi^2) -> -phi^4/4 dla eps->0."""
    return GAMMA/3*p**3 - GAMMA/4*p**4/(1.0 + eps*p**2) + LAM_BASE/6*(p-1)**6

def correction_for_V1(eps):
    """Wartosc korekty aby V_soft(1)=V_current(1)."""
    return V1 - V_soft_base(1.0, eps)

def V_soft(p, eps):
    """V_soft z korekcja V(1)=const."""
    corr = correction_for_V1(eps)
    return V_soft_base(p, eps) + corr

def dV_soft_base(p, eps):
    denom = (1.0 + eps*p**2)
    return GAMMA*p**2 - GAMMA*p**3/denom + GAMMA*p**5*eps*2.0/(4.0*denom**2) + LAM_BASE*(p-1)**5

def dV_soft(p, eps):
    return dV_soft_base(p, eps)  # korekcja jest stala, nie zmienia pochodnej

# Weryfikacja
print(f"  Warunki V_soft przy phi=1 dla roznych eps:")
for eps in [0.01, 0.1, 0.5, 1.0, 5.0]:
    V_1  = V_soft(1.0, eps)
    dV_1 = dV_soft(1.0, eps)
    d2V_1 = (V_soft(1.0+1e-5, eps) - 2*V_soft(1.0, eps) + V_soft(1.0-1e-5, eps)) / 1e-10
    print(f"    eps={eps:.2f}: V(1)={V_1:.6f}, V'(1)={dV_1:.4e}, V''(1)={d2V_1:.4f}")

print()
print(f"  V_soft(phi) - V1 dla roznych phi:")
for eps in [0.1, 1.0, 5.0]:
    print(f"  eps={eps}:")
    for phi in [1.0, 2.0, 5.0, 10.0, 20.0, 50.0]:
        diff = V_soft(phi, eps) - V_soft(1.0, eps)
        print(f"    phi={phi:5.1f}: V-V1={diff:.4e}", end="")
    print()
print()

# ============================================================
# g(K) DLA V_soft Z ROZNYM eps
# ============================================================
print("="*65)
print("Krzywa g(K) dla V_soft z roznym eps:")
print()

def ode_rhs_vmod(r, y, V1_eff, dV_func):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_func(phi) / kfac
             + ALPHA * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]

def phi_at_rmax_vmod(psi_core, K, dV_func, n_eval=2000):
    dphi0 = -K / A_GAM**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval = A_GAM * (R_MAX/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(lambda r, y: ode_rhs_vmod(r, y, None, dV_func),
                    [A_GAM, R_MAX], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < R_MAX * 0.99:
        return np.nan
    return sol.y[0, -1]

def compute_g_vmod(psi_core, K, V_func, dV_func, n_eval=5000):
    dphi0 = -K / A_GAM**2
    r_eval = A_GAM * (R_MAX/A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(lambda r, y: ode_rhs_vmod(r, y, None, dV_func),
                    [A_GAM, R_MAX], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_eval)
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    V1_eff = V_func(1.0)
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_func(phi)-V1_eff)*r**2, r)
    return (Ek+Ep)/(4*np.pi*K) - 1.0, Ek+Ep

def count_zeros_gK(V_func, dV_func, K_vals):
    """Zlicz zera g(K) na galezi glownej."""
    g_vals   = []
    psi_prev = 1.24  # start blisko K*1

    for K in K_vals:
        # Szukaj psi_core w okolicach psi_prev
        psi_lo = max(1.001, psi_prev - 0.4)
        psi_hi = psi_prev + 0.4
        F_lo = phi_at_rmax_vmod(psi_lo, K, dV_func) - 1.0
        F_hi = phi_at_rmax_vmod(psi_hi, K, dV_func) - 1.0

        psi_z = np.nan
        for w in [0.4, 0.8, 1.5, 3.0]:
            psi_lo2 = max(1.001, psi_prev - w)
            psi_hi2 = psi_prev + w
            F_lo2 = phi_at_rmax_vmod(psi_lo2, K, dV_func) - 1.0
            F_hi2 = phi_at_rmax_vmod(psi_hi2, K, dV_func) - 1.0
            if np.isfinite(F_lo2) and np.isfinite(F_hi2) and F_lo2 * F_hi2 < 0:
                try:
                    psi_z = brentq(
                        lambda p: phi_at_rmax_vmod(p, K, dV_func) - 1.0,
                        psi_lo2, psi_hi2, xtol=1e-5, maxiter=20
                    )
                    break
                except Exception:
                    pass
            elif np.isfinite(F_lo2) and np.isfinite(F_hi2):
                # Scan inside
                psi_scan = np.linspace(psi_lo2, psi_hi2, 30)
                F_scan = [phi_at_rmax_vmod(p, K, dV_func) - 1.0 for p in psi_scan]
                for i in range(len(psi_scan)-1):
                    if np.isfinite(F_scan[i]) and np.isfinite(F_scan[i+1]) and F_scan[i]*F_scan[i+1] < 0:
                        try:
                            psi_z = brentq(
                                lambda p: phi_at_rmax_vmod(p, K, dV_func) - 1.0,
                                psi_scan[i], psi_scan[i+1], xtol=1e-5, maxiter=20
                            )
                            break
                        except Exception:
                            pass
                if not np.isnan(psi_z):
                    break

        if not np.isnan(psi_z):
            dpsi = abs(psi_z - psi_prev)
            if dpsi > 0.8 * psi_prev:  # duzy skok = zmiana galezi
                g_vals.append(np.nan)
                continue
            g, E = compute_g_vmod(psi_z, K, V_func, dV_func)
            g_vals.append(g)
            psi_prev = psi_z
        else:
            g_vals.append(np.nan)

    return np.array(g_vals)

# K scan
K_vals = np.concatenate([
    np.logspace(-3, -2, 12),
    np.logspace(-2, -1, 20),
    np.logspace(-1, 0.5, 15),
])
K_vals = np.unique(K_vals)

eps_values = [0.0, 0.1, 0.5, 2.0, 10.0]
g_all = {}

for eps in eps_values:
    if eps == 0.0:
        V_f  = lambda p: V_current(p)
        dV_f = lambda p: dV_current(p)
        label = 'V_current (eps=0)'
    else:
        V_f  = lambda p, e=eps: V_soft(p, e)
        dV_f = lambda p, e=eps: dV_soft(p, e)
        label = f'V_soft eps={eps}'

    g_arr = count_zeros_gK(V_f, dV_f, K_vals)
    g_all[eps] = g_arr

    # Policz zera g
    zeros_count = 0
    for i in range(len(K_vals)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            dpsi_check = True  # zakladamy ciaglosc
            zeros_count += 1

    print(f"  eps={eps:.1f}: {zeros_count} zer g(K)  [{label}]")
    if zeros_count > 0:
        for i in range(len(K_vals)-1):
            gi, gj = g_arr[i], g_arr[i+1]
            if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
                K_approx = (K_vals[i]*K_vals[i+1])**0.5
                print(f"    K* ~ {K_approx:.4e}  (g: {gi:.3f} -> {gj:.3f})")

print()

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(f'P15: g(K) dla roznych modyfikacji V_mod\n'
             f'(alpha={ALPHA}, a_Gam={A_GAM})', fontsize=11, fontweight='bold')

colors = ['black', 'blue', 'green', 'orange', 'red']

ax = axes[0]
for eps, col in zip(eps_values, colors):
    g_arr = g_all[eps]
    mask  = np.isfinite(g_arr)
    lbl = 'oryg.' if eps == 0 else f'eps={eps}'
    ax.semilogx(K_vals[mask], np.clip(g_arr[mask], -10, 20), '.-',
                color=col, lw=1.5, ms=4, label=lbl)
ax.axhline(0, color='black', lw=1.2, linestyle='--', alpha=0.7)
ax.axvline(0.010414, color='gray', lw=1, linestyle=':', alpha=0.5, label='K*1')
ax.set_xlabel('K'); ax.set_ylabel('g = E/(4*pi*K) - 1')
ax.set_title('g(K) na galezi glownej -- rozne V_soft(eps)', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3); ax.set_ylim(-10, 15)

# Potencjaly
ax = axes[1]
phi_arr = np.linspace(0.5, 10.0, 200)
for eps, col in zip(eps_values, colors):
    if eps == 0.0:
        V_arr = np.array([V_current(p) - V1 for p in phi_arr])
        lbl = 'oryg.'
    else:
        V_arr = np.array([V_soft(p, eps) - V_soft(1.0, eps) for p in phi_arr])
        lbl = f'eps={eps}'
    ax.plot(phi_arr, np.clip(V_arr, -20, 10), '-', color=col, lw=1.5, label=lbl)
ax.axhline(0, color='black', lw=0.8, linestyle='--', alpha=0.5)
ax.axvline(1.0, color='gray', lw=0.8, linestyle=':', alpha=0.5)
ax.set_xlabel('phi'); ax.set_ylabel('V_mod(phi) - V_mod(1)')
ax.set_title('Ksztalt V_mod dla roznych eps', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_ylim(-15, 5); ax.set_xlim(0.5, 8)

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
print("PODSUMOWANIE p15:")
print()
print("V_soft modyfikacja: V_mod_new = phi^3/3 - phi^4/(4*(1+eps*phi^2)) + lam*(phi-1)^6/6")
print("Cel: eps duzy -> V_soft ograniczona -> E > 0 -> mozliwe wiecej zer g(K)")
print()
print("Wyniki (liczba zer g(K) na galezi glownej):")
for eps in eps_values:
    g_arr = g_all[eps]
    cnt = 0
    for i in range(len(K_vals)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            cnt += 1
    lbl = 'oryg.' if eps == 0 else f'eps={eps:.1f}'
    print(f"  {lbl}: {cnt} zer")
