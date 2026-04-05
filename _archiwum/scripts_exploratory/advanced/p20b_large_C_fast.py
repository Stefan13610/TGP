"""
p20b_large_C_fast.py
====================
HIPOTEZA: r21(C) -> 207 dla C -> inf (zbieznosc do Yukawa).

Zoptymalizowana wersja p20:
- C do 30 (r_max pozostaje rozumny)
- Twardy limit r_max = 400 (wydajnosc)
- Skupiony skan K wokol K*2(Yukawa) ~ 2.03

UZASADNIENIE MATEMATYCZNE:
  Na linii a_Gam = C*K:
  psi_core(K) = 1 + (1/C)*exp(-M_EFF*C*K)
  Dla C=3.84:  odchylenie 25.7%  -> r21=28.7
  Dla C=20:    odchylenie  4.7%  -> r21=???
  Dla C->inf:  odchylenie  0%    -> r21=207 (Yukawa limit)
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
V1       = GAMMA/3 - GAMMA/4

R_MAX_HARD = 400.0  # twardy limit (zamiast 8000 dla duzych C)

print("P20b: r21(C) -> 207? Zbieznosc do Yukawa (limit duzego C)")
print("=" * 65)
print(f"  M_eff={M_EFF:.4f}, K*1(ODE)=0.010414, K*2(Yukawa)=2.0327")
print(f"  r_max max: {R_MAX_HARD}")
print()

def V_mod(p): return GAMMA/3*p**3 - GAMMA/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA*p**2 - GAMMA*p**3 + LAM*(p-1)**5

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi)/kfac + ALPHA*dphi**2/(2.0*phi**2*kfac) - (2.0/r)*dphi)
    return [dphi, ddphi]

def phi_at_rmax(psi_core, K, a_gam, r_max, n_eval=1500):
    dphi0 = -K / a_gam**2
    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e7
    ev_hi.terminal = True; ev_hi.direction = 1
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam*(1+1e-12), r_max*(1-1e-12))
    sol = solve_ivp(ode_rhs, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-8, atol=1e-10,
                    t_eval=r_eval, events=[ev_lo, ev_hi])
    if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
        return np.nan
    return sol.y[0, -1]

def compute_g(psi_core, K, a_gam, r_max=None, n_eval=3000):
    if r_max is None:
        r_max = min(max(80.0, 20*a_gam, 30/M_EFF), R_MAX_HARD)
    dphi0 = -K / a_gam**2
    r_eval_raw = a_gam * (r_max/a_gam)**np.linspace(0, 1, n_eval)
    r_eval = np.clip(r_eval_raw, a_gam*(1+1e-12), r_max*(1-1e-12))
    sol = solve_ivp(ode_rhs, [a_gam, r_max], [psi_core, dphi0],
                    method='DOP853', rtol=1e-8, atol=1e-10, t_eval=r_eval)
    r   = sol.t
    phi = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return (Ek+Ep)/(4*np.pi*K) - 1.0

def find_psi_core(K, a_gam, r_max=None):
    if r_max is None:
        r_max = min(max(80.0, 20*a_gam, 30/M_EFF), R_MAX_HARD)
    psi_guess = 1.0 + K/a_gam * np.exp(-M_EFF*a_gam)
    psi_z = np.nan
    for window in [0.5, 1.5, 3.0]:
        psi_lo = max(1.001, psi_guess - window)
        psi_hi = psi_guess + window
        psi_scan = np.linspace(psi_lo, psi_hi, 40)
        F_scan = [phi_at_rmax(p, K, a_gam, r_max) - 1.0 for p in psi_scan]
        for i in range(len(psi_scan)-1):
            if (np.isfinite(F_scan[i]) and np.isfinite(F_scan[i+1])
                    and F_scan[i]*F_scan[i+1] < 0):
                try:
                    psi_z = brentq(
                        lambda p: phi_at_rmax(p, K, a_gam, r_max) - 1.0,
                        psi_scan[i], psi_scan[i+1], xtol=1e-5, maxiter=20)
                except Exception: pass
                if not np.isnan(psi_z): break
        if not np.isnan(psi_z): break
    return psi_z

# Yukawa model g(K) dla porownania
def g_yukawa(K, a_gam=0.040):
    """Aproksymacja Yukawa: E ~ 4*pi*K * K_Yukawa_factor."""
    msp = M_EFF  # masa ekranowania (pelne ODE)
    r_max = min(max(20.0, 20*a_gam, 30/msp), 300.0)
    r = np.linspace(a_gam, r_max, 3000)
    phi  = 1.0 + K * np.exp(-msp * r) / r
    phi  = np.maximum(phi, 1e-10)
    dphi = K * np.exp(-msp * r) * (-msp * r - 1.0) / r**2
    Ek = 4*np.pi * np.trapezoid(0.5 * dphi**2 * (1.0 + ALPHA/phi) * r**2, r)
    Ep = 4*np.pi * np.trapezoid((V_mod(phi)-V1) * r**2, r)
    return (Ek+Ep)/(4*np.pi*K) - 1.0

# ============================================================
# GLOWNY SKAN: C od 3.84 do 30
# ============================================================
# Skan K skupiony wokol [K*1, K*2_Yukawa] = [0.01, 2.03]
K_scan = np.unique(np.concatenate([
    np.logspace(-2.2, -0.5, 20),   # 0.006 .. 0.316
    np.logspace(-0.5, 0.5, 20),    # 0.316 .. 3.16
]))

C_vals = [3.841, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0]

K1_ref = 0.010414

results = {}

print(f"{'C':>8}  {'delta_core%':>12}  {'K*2':>10}  {'r21':>8}  {'psi@K*2':>10}")
print("-"*58)

for C in C_vals:
    d_core = 100/C * np.exp(-M_EFF*C*K1_ref)
    g_vals  = []
    psi_vals = []

    for K in K_scan:
        a_gam = C * K
        r_max = min(max(80.0, 20*a_gam, 30/M_EFF), R_MAX_HARD)
        psi_c = find_psi_core(K, a_gam, r_max)
        if np.isnan(psi_c):
            g_vals.append(np.nan); psi_vals.append(np.nan); continue
        g_v = compute_g(psi_c, K, a_gam, r_max)
        g_vals.append(g_v); psi_vals.append(psi_c)

    g_arr   = np.array(g_vals)
    psi_arr = np.array(psi_vals)

    # Znajdz zera sign-change
    zeros = []
    psi_at_zero = []
    for i in range(len(K_scan)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)): continue
        if gi*gj < 0:
            if not (np.isnan(psi_arr[i]) or np.isnan(psi_arr[i+1])):
                if abs(psi_arr[i+1]-psi_arr[i]) > 1.5: continue
            K_z = (K_scan[i]*K_scan[i+1])**0.5
            psi_z = 0.5*(psi_arr[i]+psi_arr[i+1]) if not (np.isnan(psi_arr[i]) or np.isnan(psi_arr[i+1])) else np.nan
            zeros.append(K_z)
            psi_at_zero.append(psi_z)

    K2 = zeros[0] if zeros else np.nan
    psi_K2 = psi_at_zero[0] if psi_at_zero else np.nan
    ratio = K2/K1_ref if not np.isnan(K2) else np.nan

    results[C] = {'K2': K2, 'ratio': ratio, 'psi_K2': psi_K2,
                  'g': g_arr, 'K': K_scan, 'psi': psi_arr, 'd_core': d_core}

    K2_str = f"{K2:.4e}" if not np.isnan(K2) else "brak"
    r21_str = f"{ratio:.2f}" if not np.isnan(ratio) else "nan"
    psi_str = f"{psi_K2:.4f}" if not np.isnan(psi_K2) else "nan"
    print(f"{C:>8.2f}  {d_core:>12.3f}%  {K2_str:>10}  {r21_str:>8}  {psi_str:>10}")

# ============================================================
# PRECYZYJNE K*2 PRZEZ BRENTQ DLA WARTOSCI GDZIE SIGN-CHANGE
# ============================================================
print()
print("Precyzyjne wyznaczenie K*2 przez brentq:")
print("-"*58)
for C, res in results.items():
    if np.isnan(res['K2']): continue
    g_arr = res['g']
    K_arr = res['K']
    psi_arr = res['psi']

    # Znajdz indeks sign-change
    for i in range(len(K_arr)-1):
        gi, gj = g_arr[i], g_arr[i+1]
        if not (np.isfinite(gi) and np.isfinite(gj)): continue
        if gi*gj < 0:
            if not (np.isnan(psi_arr[i]) or np.isnan(psi_arr[i+1])):
                if abs(psi_arr[i+1]-psi_arr[i]) > 1.5: continue
            K_lo, K_hi = K_arr[i], K_arr[i+1]
            try:
                def g_func_C(K_v, C_v=C):
                    a_g = C_v*K_v
                    r_mx = min(max(80.0, 20*a_g, 30/M_EFF), R_MAX_HARD)
                    psi_c = find_psi_core(K_v, a_g, r_mx)
                    if np.isnan(psi_c): return 999.0
                    return compute_g(psi_c, K_v, a_g, r_mx)

                K2_precise = brentq(g_func_C, K_lo, K_hi, xtol=1e-4, maxiter=15)
                a_g2 = C*K2_precise
                psi_c2 = find_psi_core(K2_precise, a_g2)
                ratio2 = K2_precise/K1_ref
                psi_dev = psi_c2 - 1.0 if not np.isnan(psi_c2) else np.nan
                results[C]['K2_precise'] = K2_precise
                results[C]['ratio_precise'] = ratio2
                results[C]['psi_K2_precise'] = psi_c2

                print(f"  C={C:5.1f}: K*2={K2_precise:.5f}, a_Gam={a_g2:.4f}, "
                      f"psi_core={psi_c2:.5f} (dev={100*psi_dev:.3f}%), r21={ratio2:.2f}")
            except Exception as ex:
                print(f"  C={C:5.1f}: brentq blad: {ex}")
            break

# ============================================================
# WYNIK: TABELA r21(C)
# ============================================================
print()
print("=" * 65)
print("TABELA: r21 = K*2/K*1 vs C  (zbieznosc do 207?)")
print("=" * 65)
print(f"  {'C':>8}  {'delta_core%':>12}  {'K*2':>10}  {'r21':>10}  {'psi_core@K*2':>14}")
print("  " + "-"*60)

for C in C_vals:
    res = results[C]
    d = res['d_core']
    K2 = res.get('K2_precise', res['K2'])
    r21 = res.get('ratio_precise', res['ratio'])
    psi = res.get('psi_K2_precise', res.get('psi_K2', np.nan))
    K2_str  = f"{K2:.5f}" if not np.isnan(K2) else "brak"
    r21_str = f"{r21:.2f}" if not np.isnan(r21) else "nan"
    psi_str = f"{psi:.5f}" if not np.isnan(psi) else "nan"
    print(f"  {C:>8.2f}  {d:>12.3f}%  {K2_str:>10}  {r21_str:>10}  {psi_str:>14}")

print()
print(f"  C=inf (Yukawa): K*2=2.03270, r21=207.00, psi_core=1.00000")
print()

# Zbieznosc
valid = [(C, results[C].get('ratio_precise', results[C]['ratio']))
         for C in C_vals if not np.isnan(results[C].get('ratio_precise', results[C]['ratio']))]
if len(valid) >= 2:
    C_arr = np.array([v[0] for v in valid])
    r_arr = np.array([v[1] for v in valid])
    is_mono = all(r_arr[i] <= r_arr[i+1] for i in range(len(r_arr)-1))
    print(f"  Monotoniczne wzrastanie r21 z C: {'TAK' if is_mono else 'NIE (anomalia)'}")
    max_r = max(r_arr)
    max_C = C_arr[np.argmax(r_arr)]
    print(f"  Maks. r21 przy C={max_C:.1f}: r21={max_r:.2f}")
    if max_r > 100:
        print(f"  -> Wyrazna zbieznosc do 207!")
    elif max_r > 50:
        print(f"  -> Umiarkowana zbieznosc; C>30 prawdopodobnie da r21->207")
    else:
        print(f"  -> Powolna zbieznosc lub brak.")

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P20b: Zbieznosc r21(C) -> 207 dla C -> inf (pelne ODE)\n'
             f'(alpha={ALPHA}, lam={LAM:.2e})',
             fontsize=11, fontweight='bold')

ax = axes[0]
C_list = [C for C in C_vals if not np.isnan(results[C].get('ratio_precise', results[C]['ratio']))]
r21_list = [results[C].get('ratio_precise', results[C]['ratio']) for C in C_list]
if C_list:
    ax.semilogx(C_list, r21_list, 'bo-', lw=2, ms=9, label='r21(C) ODE')
ax.axhline(207.0, color='red',    lw=2,   linestyle='--', label='r21=207 (Yukawa)')
ax.axhline(28.7,  color='orange', lw=1.5, linestyle=':', label='r21=28.7 (C=3.84)')
ax.axvline(3.841, color='green',  lw=1,   linestyle=':', alpha=0.7)
ax.set_xlabel('C'); ax.set_ylabel('r21 = K*2/K*1')
ax.set_title('r21(C) -- zbieznosc do 207', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[1]
for C in [3.841, 7.0, 15.0, 30.0]:
    if C not in results: continue
    g_arr = results[C]['g']
    K_arr = results[C]['K']
    mask = np.isfinite(g_arr)
    g_clip = np.clip(g_arr[mask], -5, 8)
    ax.semilogx(K_arr[mask], g_clip, '.-', lw=1.5, ms=3, label=f'C={C:.0f}')
ax.axhline(0, color='black', lw=1.2, linestyle='--')
ax.axvline(0.010414, color='green', lw=1, linestyle=':', alpha=0.5, label='K*1')
ax.axvline(2.0327,   color='red',   lw=1, linestyle=':', alpha=0.5, label='K*2(Yuk)')
ax.set_xlabel('K'); ax.set_ylabel('g')
ax.set_title('g(K) dla roznych C', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(-5, 8)

ax = axes[2]
K2_precise_list = [results[C].get('K2_precise', results[C]['K2']) for C in C_vals
                   if not np.isnan(results[C].get('K2_precise', results[C]['K2']))]
C_K2_list = [C for C in C_vals if not np.isnan(results[C].get('K2_precise', results[C]['K2']))]
if C_K2_list:
    ax.loglog(C_K2_list, K2_precise_list, 'bo-', lw=2, ms=7, label='K*2(C) ODE')
ax.axhline(2.0327, color='red',    lw=1.5, linestyle='--', label='K*2=2.033 (Yukawa)')
ax.axhline(0.2990, color='orange', lw=1.5, linestyle=':', label='K*2=0.299 (C=3.84)')
# Yukawa analityczne przewidywanie K*2~C^0:
ax.set_xlabel('C'); ax.set_ylabel('K*2 (samospojny)')
ax.set_title('K*2(C)', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
