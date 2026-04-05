"""
p20_large_C_limit.py
====================
HIPOTEZA: Dla C -> inf na linii Schwarzschilda a_Gam = C*K:
  - psi_core -> 1  (ODE zbiega do Yukawa)
  - r21 = K*2/K*1 -> 207  (hierarchia Yukawa)

UZASADNIENIE:
  psi_core(K; C) = 1 + (1/C)*exp(-M_EFF*C*K)
  Dla C=3.84:  delta = 25.7%  (nieliniowy, r21=28.7)
  Dla C=20:    delta =  4.7%  (blisko Yukawa)
  Dla C=100:   delta =  0.7%  (prawie Yukawa, r21->207)

CEL p20:
  Sprawdzic r21(C) dla C w [3.84, 500].
  Znalezc C* gdzie r21=207.
  Zbadac czy K*2(C) jest monotoniczne wzrastajace.

UWAGA TECHNICZNA:
  Dla duzego C, skan K musi siegnac do K*2(Yukawa) ~ 2.03.
  Zakres K: logspace(-3, 1.0, 50) = [0.001, 10]
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

ALPHA    = 8.5616
LAM      = 5.501357e-06
GAMMA    = 1.0
M_EFF    = 1.0 / np.sqrt(1.0 + ALPHA)
R_MAX_BASE = 80.0
V1       = GAMMA/3 - GAMMA/4

print("P20: Limit duzego C -- zbieznosc do Yukawa r21=207")
print("=" * 65)
print(f"  M_eff = {M_EFF:.4f},  K*1(ODE) = 0.010414,  K*2(Yukawa) = 2.0327")
print(f"  Hipoteza: r21(C) -> 207 dla C -> inf")
print()

# ============================================================
def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi)/kfac + ALPHA*dphi**2/(2.0*phi**2*kfac) - (2.0/r)*dphi)
    return [dphi, ddphi]

def phi_at_rmax(psi_core, K, a_gam, r_max, n_eval=2000):
    dphi0 = -K / a_gam**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam*(1+1e-12), r_max*(1-1e-12))
    sol = solve_ivp(ode_rhs, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def compute_g(psi_core, K, a_gam, n_eval=5000):
    dphi0 = -K / a_gam**2
    r_max = max(R_MAX_BASE, 20*a_gam, 30/M_EFF)
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam*(1+1e-12), r_max*(1-1e-12))
    sol = solve_ivp(ode_rhs, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-9, atol=1e-11, t_eval=r_eval)
    r   = sol.t
    phi = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return (Ek+Ep)/(4*np.pi*K) - 1.0

def find_psi_core(K, a_gam, n_scan=50):
    r_max = max(R_MAX_BASE, 20*a_gam, 30/M_EFF)
    psi_guess = 1.0 + K/a_gam * np.exp(-M_EFF*a_gam)
    psi_z = np.nan
    for window in [0.5, 1.5, 3.0, 5.0]:
        psi_lo = max(1.001, psi_guess - window)
        psi_hi = psi_guess + window
        psi_scan = np.linspace(psi_lo, psi_hi, n_scan)
        F_scan = [phi_at_rmax(p, K, a_gam, r_max) - 1.0 for p in psi_scan]
        for i in range(len(psi_scan)-1):
            if (np.isfinite(F_scan[i]) and np.isfinite(F_scan[i+1])
                    and F_scan[i]*F_scan[i+1] < 0):
                try:
                    psi_z = brentq(
                        lambda p: phi_at_rmax(p, K, a_gam, r_max) - 1.0,
                        psi_scan[i], psi_scan[i+1],
                        xtol=1e-5, maxiter=20)
                except Exception: pass
                if not np.isnan(psi_z): break
        if not np.isnan(psi_z): break
    return psi_z


# ============================================================
# GLOWNY SKAN: C od 3.84 do 500
# ============================================================
print("Skan C = [3.84, 5, 7, 10, 15, 20, 30, 50, 100, 200, 500]")
print("-" * 65)

# Dla kazdego C: liczymy g(K) na szerokim zakresie K
# K*1 ~ 0.010414 (tangentalne zero, maloznaczace)
# K*2 -- szukamy znaku g(K) dla K w [0.01, 10]

C_vals = [3.841, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]

# Szerszy skan K (do 10, bo K*2 Yukawa = 2.03)
K_scan_base = np.logspace(-2, 1.3, 40)  # 0.01 do ~20

results_C = {}

K1_ref = 0.010414  # z p14c

for C in C_vals:
    print(f"\n  C = {C:.1f}  (delta_core_K1 = {100/C*np.exp(-M_EFF*C*K1_ref):.2f}%)")
    g_vals  = []
    psi_vals = []

    for K in K_scan_base:
        a_gam = C * K
        psi_c = find_psi_core(K, a_gam)
        if np.isnan(psi_c):
            g_vals.append(np.nan); psi_vals.append(np.nan); continue
        g_v = compute_g(psi_c, K, a_gam)
        g_vals.append(g_v); psi_vals.append(psi_c)
        tag = ''
        if np.isfinite(g_v) and abs(g_v) < 0.05: tag = ' <z>'
        print(f"    K={K:.4f}: psi={psi_c:.4f}, g={g_v:+.4e}{tag}")

    g_arr   = np.array(g_vals)
    psi_arr = np.array(psi_vals)

    # Znajdz zera (sign changes)
    zeros = []
    for i in range(len(K_scan_base)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)): continue
        if gi*gj < 0:
            # Sprawdz ciaglosc psi
            if not (np.isnan(psi_arr[i]) or np.isnan(psi_arr[i+1])):
                if abs(psi_arr[i+1]-psi_arr[i]) > 1.5: continue
            K_z = (K_scan_base[i]*K_scan_base[i+1])**0.5
            zeros.append(K_z)

    # Tez sprawdz tangentalne zero (minimum g blisko K1_ref)
    # i czy g jest blisko zera w zakresie K~K1_ref
    g_near_K1 = g_arr[np.argmin(np.abs(K_scan_base - K1_ref))] if len(g_arr)>0 else np.nan
    print(f"    g w K~K*1: {g_near_K1:.4e}")
    print(f"    Zera sign-change: {[f'{z:.4e}' for z in zeros]}")

    # Oblicz ratio
    K2_cand = [z for z in zeros if z > 0.05]
    ratio = K2_cand[0]/K1_ref if K2_cand else np.nan

    results_C[C] = {
        'g': g_arr, 'psi': psi_arr, 'zeros': zeros,
        'K2': K2_cand[0] if K2_cand else np.nan,
        'ratio': ratio,
        'g_near_K1': g_near_K1
    }
    if not np.isnan(ratio):
        print(f"    K*2 ~ {K2_cand[0]:.4e}, r21 = K*2/K*1 = {ratio:.1f}")

# ============================================================
# WYKRES
# ============================================================
print("\nTworze wykres...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P20: Zbieznosc r21 do 207 dla C -> inf (pelne ODE)\n'
             f'(alpha={ALPHA}, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

ax = axes[0]
C_list   = [C for C in C_vals if not np.isnan(results_C[C]['ratio'])]
r21_list = [results_C[C]['ratio'] for C in C_list]
if C_list:
    ax.semilogx(C_list, r21_list, 'bo-', lw=2, ms=8, label='r21(C) ODE')
    ax.axhline(207.0, color='red',    lw=2,   linestyle='--', label='r21=207 (Yukawa)')
    ax.axhline(28.7,  color='orange', lw=1.5, linestyle=':', label='r21=28.7 (C=3.84)')
    ax.axvline(3.841, color='green',  lw=1.5, linestyle=':', alpha=0.7, label='C=3.841')
ax.set_xlabel('C (stala Schwarzschilda)'); ax.set_ylabel('K*2/K*1')
ax.set_title('r21 vs C -- czy zbieznosc do 207?', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

ax = axes[1]
# g(K) dla kilku C wartosci
for C in [3.841, 10.0, 50.0, 200.0]:
    if C not in results_C: continue
    g_arr = results_C[C]['g']
    mask = np.isfinite(g_arr)
    g_clip = np.clip(g_arr[mask], -5, 8)
    ax.semilogx(K_scan_base[mask], g_clip, '.-', lw=1.5, ms=3, label=f'C={C:.0f}')
ax.axhline(0, color='black', lw=1.2, linestyle='--')
ax.axvline(0.010414, color='green', lw=1, linestyle=':', alpha=0.5, label='K*1')
ax.axvline(2.0327,   color='red',   lw=1, linestyle=':', alpha=0.5, label='K*2(Yuk)')
ax.set_xlabel('K'); ax.set_ylabel('g')
ax.set_title('g(K) dla roznych C', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(-5, 8)

ax = axes[2]
# K*2(C) vs C
K2_list = [results_C[C]['K2'] for C in C_vals if not np.isnan(results_C[C].get('K2', np.nan))]
C_list2 = [C for C in C_vals if not np.isnan(results_C[C].get('K2', np.nan))]
if C_list2:
    ax.loglog(C_list2, K2_list, 'bo-', lw=2, ms=7, label='K*2(C) ODE')
ax.axhline(2.0327, color='red', lw=1.5, linestyle='--', label='K*2=2.033 (Yukawa)')
ax.axhline(0.2990, color='orange', lw=1.5, linestyle=':', label='K*2=0.299 (C=3.84)')
ax.set_xlabel('C'); ax.set_ylabel('K*2')
ax.set_title('K*2(C) -- zbieznosc do K*2(Yukawa)?', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")

# ============================================================
# PODSUMOWANIE
# ============================================================
print()
print("=" * 65)
print("PODSUMOWANIE p20: r21(C) i zbieznosc do Yukawa")
print("=" * 65)
print()
print(f"  {'C':>8}  {'delta_core%':>12}  {'K*2':>10}  {'r21':>8}")
print("  " + "-"*46)
for C in C_vals:
    d = 100/C*np.exp(-M_EFF*C*K1_ref)
    r = results_C[C]
    K2_str = f"{r['K2']:.4e}" if not np.isnan(r.get('K2', np.nan)) else "brak"
    r21_str = f"{r['ratio']:.1f}" if not np.isnan(r['ratio']) else "nan"
    print(f"  {C:>8.1f}  {d:>12.3f}%  {K2_str:>10}  {r21_str:>8}")

print()
print("  Yukawa limit (C->inf):")
print(f"    K*2 -> {2.0327:.4f},  r21 -> 207")
print()

# Sprawdz czy r21 jest monotoniczne
ratios_valid = [(C, results_C[C]['ratio']) for C in C_vals
                if not np.isnan(results_C[C]['ratio'])]
if len(ratios_valid) >= 2:
    r_vals = [r for _, r in ratios_valid]
    is_mono = all(r_vals[i] <= r_vals[i+1] for i in range(len(r_vals)-1))
    print(f"  r21 monotonicznie wzrastajace z C: {'TAK' if is_mono else 'NIE'}")
    max_r = max(r_vals)
    max_C = ratios_valid[r_vals.index(max_r)][0]
    print(f"  Maks. r21 znalezione: {max_r:.1f} przy C={max_C:.1f}")
    if max_r >= 180:
        print(f"  -> BLISKIE 207! C~{max_C:.0f} daje r21~{max_r:.1f}")
    elif max_r >= 100:
        print(f"  -> Wyrazna zbieznosc ku 207; wymaga wiekszego C")
    else:
        print(f"  -> Wciaz daleko od 207; zbieznosc powolna lub brak")
