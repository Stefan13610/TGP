"""
p9_branch_diagnostics.py
========================
Diagnostyka gałęzi ψ_core(K) w pełnym ODE TGP.

PROBLEM z p6: find_psi_core() jest niedeterministyczne — przy K>0.034
F(ψ_core) = φ(R_max; ψ_core) − 1 ma WIELE zer, a skan 30-punktowy
losowo wybiera różne gałęzie, dając g(K) z przeskokami.

CEL p9:
  Dla wybranych K: zmapować CAŁĄ funkcję F(ψ_core) na szerokim zakresie,
  znaleźć WSZYSTKIE gałęzie, policzyć zera, ocenić g na każdej gałęzi.

WYNIKI p6 do weryfikacji:
  K₁=0.00982 → g=-9.2 (niespójny z ansatzem Yukawa)
  K₂=2.033   → brak ψ_core
  K₃=34.14   → g=-94
  K*≈0.033   → zmiana znaku g (prawdopodobnie jedyne zero)
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
R_MAX  = 50.0

print(f"Parametry: α={ALPHA}, a_Γ={A_GAM}, λ={LAM:.4e}")
print(f"m_eff = {M_EFF:.4f},  R_max = {R_MAX}")
print()

def V_mod(phi):
    return GAMMA/3*phi**3 - GAMMA/4*phi**4 + LAM/6*(phi-1)**6

def dV_mod(phi):
    return GAMMA*phi**2 - GAMMA*phi**3 + LAM*(phi-1)**5

V1 = V_mod(1.0)

def ode_rhs(r, y):
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac = 1.0 + ALPHA / phi
    ddphi = (dV_mod(phi) / kfac
             + ALPHA * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]

def phi_at_rmax(psi_core, K, r_max=R_MAX, n_eval=2000):
    """Outward integration, return φ(R_max)."""
    phi0  = psi_core
    dphi0 = -K / A_GAM**2

    def ev_lo(r, y): return y[0] - 1e-5
    ev_lo.terminal = True; ev_lo.direction = -1
    def ev_hi(r, y): return y[0] - 1e8
    ev_hi.terminal = True; ev_hi.direction = 1

    r_eval = A_GAM * (r_max / A_GAM)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(
        ode_rhs, [A_GAM, r_max], [phi0, dphi0],
        method='DOP853', rtol=1e-9, atol=1e-11,
        t_eval=r_eval, events=[ev_lo, ev_hi]
    )
    if len(sol.t) == 0:
        return np.nan
    if sol.t[-1] < r_max * 0.99:   # terminated early
        return np.nan
    return sol.y[0, -1]

def profile_energy(psi_core, K, n_eval=4000):
    """Compute soliton energy E[phi] on outward profile."""
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
    K_eff = r[0]**2 * abs(dphi[0])
    return Ek+Ep, Ek, Ep, K_eff


# ============================================================
# MAPA F(psi_core) DLA WYBRANYCH K
# ============================================================
test_Ks = {
    'K*_ODE ≈ 0.033':   0.033,
    'K1_Yukawa=0.00982': 0.009820,
    'K2_Yukawa=2.033':   2.032728,
    'K3_Yukawa=34.14':   34.14450,
}

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('p9: Mapa F(ψ_core) = φ(R_max;ψ_core)−1 dla różnych K\n'
             f'(α={ALPHA}, a_Γ={A_GAM}, λ={LAM:.2e})',
             fontsize=11, fontweight='bold')

all_branches = {}

for ax, (label, K) in zip(axes.flat, test_Ks.items()):
    print(f"{'='*60}")
    print(f"K = {K:.5e}  ({label})")
    print(f"{'='*60}")

    # Szeroki zakres ψ_core: od 1.001 do max(10, 3*K/A_GAM)
    psi_max = max(20.0, 5.0 * K / A_GAM)
    # Użyj skali logarytmicznej
    psi_vals = np.concatenate([
        np.linspace(1.001, 5.0, 200),
        np.linspace(5.0, min(50.0, psi_max), 100)
    ])
    psi_vals = np.unique(psi_vals)

    print(f"  Skanowanie {len(psi_vals)} wartości ψ_core ∈ [{psi_vals[0]:.3f}, {psi_vals[-1]:.1f}]")

    F_vals = np.array([phi_at_rmax(p, K) - 1.0 for p in psi_vals])

    # Znajdź wszystkie zmiany znaku
    zeros = []
    for i in range(len(F_vals)-1):
        Fi, Fj = F_vals[i], F_vals[i+1]
        if np.isfinite(Fi) and np.isfinite(Fj) and Fi * Fj < 0:
            try:
                psi_z = brentq(
                    lambda p: phi_at_rmax(p, K) - 1.0,
                    psi_vals[i], psi_vals[i+1],
                    xtol=1e-5, rtol=1e-5, maxiter=40
                )
                phi_end = phi_at_rmax(psi_z, K)
                E, Ek, Ep, K_eff = profile_energy(psi_z, K)
                g = E / (4*np.pi*K) - 1.0
                zeros.append({
                    'psi_core': psi_z,
                    'phi_end': phi_end,
                    'E': E, 'Ek': Ek, 'Ep': Ep,
                    'K_eff': K_eff, 'g': g
                })
            except Exception as e:
                pass

    all_branches[label] = zeros

    print(f"  Znaleziono {len(zeros)} gałęzi (zer F):")
    print(f"  {'ψ_core':>10}  {'φ(R)':>8}  {'E':>12}  {'g_ODE':>10}  {'samospójny?'}")
    print("  " + "-"*60)
    for z in zeros:
        ok = '✓' if abs(z['g']) < 0.1 else ('~' if abs(z['g']) < 1.0 else '✗')
        print(f"  {z['psi_core']:>10.4f}  {z['phi_end']:>8.5f}  {z['E']:>12.4e}  {z['g']:>10.4e}  {ok}")

    if zeros:
        g_vals_z = [z['g'] for z in zeros]
        g_min = min(g_vals_z, key=abs)
        print(f"\n  Najlepszy g (|g| min): {g_min:.4e}")
        closest = zeros[np.argmin([abs(z['g']) for z in zeros])]
        print(f"  → ψ_core = {closest['psi_core']:.4f},  K*={K:.4e},  g={closest['g']:.4e}")
    print()

    # Panel wykresu
    mask_fin = np.isfinite(F_vals)
    F_clip = np.clip(F_vals, -10, 50)
    ax.plot(psi_vals[mask_fin], F_clip[mask_fin], 'b-', lw=1.2, alpha=0.8)
    ax.axhline(0, color='red', lw=1, linestyle='--', label='F=0 (cel)')
    for z in zeros:
        g = z['g']
        color = 'green' if abs(g) < 0.1 else ('orange' if abs(g) < 2 else 'gray')
        ax.axvline(z['psi_core'], color=color, lw=1.2, alpha=0.7,
                   label=f"ψ*={z['psi_core']:.2f}, g={g:.2f}")
    ax.set_xlabel('ψ_core')
    ax.set_ylabel('F(ψ_core) = φ(R_max)−1')
    ax.set_title(f'{label}', fontsize=9)
    ax.legend(fontsize=7, loc='upper right')
    ax.set_ylim(-5, 20)
    ax.grid(True, alpha=0.3)


plt.tight_layout()
out_path = __file__.replace('.py', '.png')
plt.savefig(out_path, dpi=120, bbox_inches='tight')
plt.close()
print(f"Zapisano: {out_path}")
print()


# ============================================================
# ZBIORCZE PODSUMOWANIE
# ============================================================
print("=" * 65)
print("PODSUMOWANIE p9: WSZYSTKIE GAŁĘZIE ψ_core(K)")
print("=" * 65)
print()
for label, zeros in all_branches.items():
    print(f"K = {label}:")
    if not zeros:
        print("  — brak gałęzi F(ψ_core)=0")
    for z in zeros:
        g = z['g']
        tag = '✓ SAMOSPÓJNY' if abs(g) < 0.05 else (
              '~ BLISKI' if abs(g) < 0.5 else '✗ niespójny')
        print(f"  ψ_core={z['psi_core']:.3f}, E={z['E']:.3e}, g={g:.4e}  {tag}")
    print()

print()
print("INTERPRETACJA:")
print("  g ≈ 0 → soliton samospójny (K jest masą cząstki TGP)")
print("  g > 0 → E > 4πK: profil 'za ciężki'")
print("  g < 0 → E < 4πK: profil 'za lekki'")
print()
print("  Jeśli K_Yukawa mają gałęzie z g≈0 → hierarchia przeżywa w ODE")
print("  Jeśli brak → hierarchia Yukawa jest artefaktem aproksymacji")
