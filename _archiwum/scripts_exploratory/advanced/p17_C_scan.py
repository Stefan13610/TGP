"""
p17_C_scan.py
=============
Jak stosunek K*2/K*1 zalezy od stalej C w a_Gam = C*K?

WYNIK p16: C=3.841 daje K*2/K*1 = 28.7 (nie 207!)
PYTANIE:   Czy dla innej C mozna uzyskac K*2/K*1 = 207?

BADAMY:
  Dla kazdego C: wyznacz g(K) na linii a_Gam=C*K
  Znajdz K*2 (drugie zero g=0) i oblicz K*2/K*1
  Sprawdz jaka C daje K*2/K*1 = 207

OCZEKIWANIE:
  Jezeli C = C* taki ze K*2/K*1 = 207, mozna to interpretowac jako
  "wymaganie fizyczne" TGP dla zachowania hierarchii Yukawa.
  Ale czy C* jest fizycznie uzasadnione?
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
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]
    Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA/phi)*r**2, r)
    Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return (Ek+Ep)/(4*np.pi*K) - 1.0

def find_K2_for_C(C, K_scan=None, verbose=False):
    """
    Dla danego C: znajdz drugie zero g(K) na linii a_Gam = C*K.
    Zwraca K*2 lub nan jesli nie ma.
    """
    if K_scan is None:
        K_scan = np.logspace(-2, 1.5, 35)

    g_vals = []
    psi_arr = []

    for K in K_scan:
        a_gam = C * K
        r_max_loc = max(R_MAX_BASE, 20*a_gam, 30/M_EFF)
        psi_guess = 1.0 + (K/a_gam)*np.exp(-M_EFF*a_gam)

        # Szukaj psi_core
        psi_z = np.nan
        for w in [0.4, 1.0, 3.0]:
            psi_lo2 = max(1.001, psi_guess - w)
            psi_hi2 = psi_guess + w
            psi_scan_loc = np.linspace(psi_lo2, psi_hi2, 30)
            F_scan_loc = [phi_at_rmax(p, K, a_gam, r_max_loc) - 1.0
                          for p in psi_scan_loc]
            for i in range(len(psi_scan_loc)-1):
                if np.isfinite(F_scan_loc[i]) and np.isfinite(F_scan_loc[i+1]):
                    if F_scan_loc[i] * F_scan_loc[i+1] < 0:
                        try:
                            psi_z = brentq(
                                lambda p: phi_at_rmax(p, K, a_gam, r_max_loc) - 1.0,
                                psi_scan_loc[i], psi_scan_loc[i+1],
                                xtol=1e-5, maxiter=20)
                            break
                        except Exception:
                            pass
            if not np.isnan(psi_z): break

        if not np.isnan(psi_z):
            g = compute_g(psi_z, K, a_gam)
            g_vals.append(g)
            psi_arr.append(psi_z)
        else:
            g_vals.append(np.nan)
            psi_arr.append(np.nan)

    g_vals = np.array(g_vals)
    psi_arr = np.array(psi_arr)

    # Znajdz zera g poza K*1 (szukamy w K > K*1/2)
    K1_approx = 0.010414
    K2_found = []

    for i in range(len(K_scan)-1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
            # Sprawdz ciaglosc
            dP = abs(psi_arr[i+1] - psi_arr[i]) if (not np.isnan(psi_arr[i]) and not np.isnan(psi_arr[i+1])) else 999
            if dP > 1.0:
                continue
            K_approx_z = (K_scan[i]*K_scan[i+1])**0.5
            K2_found.append(K_approx_z)
            if verbose:
                print(f"    K*2 ~ {K_approx_z:.4e}  g: {gi:.3f} -> {gj:.3f}")

    return K2_found, g_vals, K_scan

# ============================================================
# KROK 1: K*1(C) dla roznych C
# ============================================================
# K*1 wyznaczone wczesniej dla C=3.841 wynosi 0.010414
# Ale dla innych C, K*1 bedzie inne!
# Pytanie: jak K*2/K*1 zalezy od C?

print("P17: Skan C -> K*2/K*1")
print("="*65)
print()

C_vals = np.array([1.0, 2.0, 3.0, 3.84, 5.0, 7.0, 10.0, 15.0, 20.0])

results = {}

for C in C_vals:
    print(f"C = {C:.2f}:")
    K_scan_C = np.logspace(-4, 2, 50)

    # Najpierw znajdz K*1 (pierwsze zero g=0 przy malym K)
    # K*1 jest tam gdzie profil solitonu jest samospojny dla MALEGO K
    # Na linii a_Gam = C*K: dla malych K a_Gam jest male -> blisko K*1(a_Gam=0.040)

    K2_list, g_arr, K_s = find_K2_for_C(C, K_scan=K_scan_C, verbose=True)

    # K*1 to minimum g (czy zero dla malych K)
    mask_fin = np.isfinite(g_arr)
    K_fin = K_s[mask_fin]
    g_fin = g_arr[mask_fin]

    # Szukaj K*1: pierwsze zero od malych K
    K1_candidate = np.nan
    for i in range(len(K_fin)-1):
        if g_fin[i]*g_fin[i+1] < 0 and K_fin[i] < 0.1:
            K1_candidate = (K_fin[i]*K_fin[i+1])**0.5
            break

    # K*2: pierwsze zero powyzej K*1
    K2_candidates = [K for K in K2_list if K > 0.05]

    ratio = np.nan
    if not np.isnan(K1_candidate) and K2_candidates:
        ratio = K2_candidates[0] / K1_candidate
    elif K2_candidates and np.isnan(K1_candidate):
        # Jesli K*1 nie znalezione, uzyj K*1=0.010414
        ratio = K2_candidates[0] / 0.010414

    results[C] = {
        'K1': K1_candidate,
        'K2': K2_candidates[0] if K2_candidates else np.nan,
        'ratio': ratio,
        'g_arr': g_arr,
        'K_scan': K_s
    }
    K2_str = f"{K2_candidates[0]:.4e}" if K2_candidates else "brak"
    print(f"  K*1 ~ {K1_candidate:.4e}, K*2 ~ {K2_str}")
    print(f"  K*2/K*1 ~ {ratio:.2f}")
    print()

# ============================================================
# WYKRES
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('P17: Hierarchia mas K*2/K*1 vs. stala C (a_Gam = C*K)\n'
             f'(alpha={ALPHA}, lam={LAM:.2e})', fontsize=11, fontweight='bold')

ax = axes[0]
C_list = list(results.keys())
ratio_list = [results[C]['ratio'] for C in C_list]
K2_list_plot = [results[C]['K2'] for C in C_list]
K1_list_plot = [results[C]['K1'] for C in C_list]
ax.semilogy(C_list, ratio_list, 'bo-', lw=1.5, ms=7)
ax.axhline(207, color='red', lw=1.5, linestyle='--', label='r21=207 (Yukawa)')
ax.axhline(28.7, color='orange', lw=1.5, linestyle=':', label='r21=28.7 (C=3.84)')
ax.axvline(3.841, color='green', lw=1.5, linestyle=':', alpha=0.7, label='C=3.841 (K*1/a1)')
ax.set_xlabel('C (w a_Gam = C*K)'); ax.set_ylabel('K*2/K*1')
ax.set_title('Stosunek mas vs. C', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[1]
for C in C_vals[:6]:
    res = results[C]
    mask = np.isfinite(res['g_arr'])
    g_clip = np.clip(res['g_arr'][mask], -5, 10)
    ax.semilogx(res['K_scan'][mask], g_clip, '.-', lw=1.2, ms=3, label=f'C={C:.1f}')
ax.axhline(0, color='black', lw=1.2, linestyle='--')
ax.set_xlabel('K'); ax.set_ylabel('g')
ax.set_title('g(K) dla roznych C', fontsize=9)
ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(-5, 8)

ax = axes[2]
ax.semilogy(C_list, [r['K2'] if not np.isnan(r['K2']) else np.nan for r in results.values()],
            'bo-', lw=1.5, ms=7, label='K*2(ODE)')
ax.semilogy(C_list, [r['K1'] if not np.isnan(r['K1']) else np.nan for r in results.values()],
            'rs-', lw=1.5, ms=7, label='K*1(ODE)')
ax.axhline(0.010414, color='green', lw=1, linestyle=':', alpha=0.7, label='K*1=0.0104')
ax.set_xlabel('C'); ax.set_ylabel('K*')
ax.set_title('K*1 i K*2 vs. C', fontsize=9)
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
out = __file__.replace('.py', '.png')
plt.savefig(out, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out}")
print()

# Podsumowanie
print("="*65)
print("PODSUMOWANIE p17:")
print()
print(f"  C=3.841 (z K*1/a_Gam1): K*2/K*1 = 28.7")
print()
print(f"  Tabela C -> K*2/K*1:")
for C in C_vals:
    r = results[C]
    print(f"    C={C:5.2f}: K*1={r['K1']:.4e}, K*2={r['K2']:.4e}, r21={r['ratio']:.2f}")
print()
print(f"  Pytanie: czy istnieje C takie ze K*2/K*1 = 207?")
ratios = [(C, results[C]['ratio']) for C in C_vals if not np.isnan(results[C]['ratio'])]
ratios_arr = np.array([(c, r) for c, r in ratios])
if ratios_arr.size > 0:
    max_ratio = ratios_arr[np.argmax(ratios_arr[:, 1])]
    print(f"  Maks. znaleziony stosunek: {max_ratio[1]:.1f} przy C={max_ratio[0]:.2f}")
    if max_ratio[1] > 180:
        print(f"  -> MOZLIWE! C~{max_ratio[0]:.2f} daje r21~{max_ratio[1]:.1f}")
    else:
        print(f"  -> NIE: maks ratio {max_ratio[1]:.1f} < 207")
        print(f"     TGP z V_mod_current NIE daje r21=207 dla zadnego C")
