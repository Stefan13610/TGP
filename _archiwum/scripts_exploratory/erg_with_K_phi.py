#!/usr/bin/env python3
"""
ERG Wetterich z field-dependent K(phi) = K_geo * psi^4
========================================================
Kluczowa modyfikacja wzgledem czystego LPA:
  - W TGP czlon kinetyczny to K(phi)*(dphi)^2 z K(phi) = K_geo*(phi/Phi0)^4
  - To zmienia propagator: G^{-1}(p) = K(psi)*p^2 + V''(psi)
  - Regulator Litim: R_k = K(psi)*(k^2 - p^2)*theta(k^2 - p^2)

Rownanie przeplywu (LPA z K(psi), d=4, Litim):
  dV_k/dt = K(psi)*k^6 / (32*pi^2 * (K(psi)*k^2 + V''_k(psi)))

Spodziewany efekt:
  - Dla duzych psi: K ~ psi^4 >> 1, wiec K*k^2 dominuje -> przeplyw ~ k^4 (ograniczony)
  - Dla malych psi: K ~ psi^4 -> 0, wiec V''/K -> inf -> przeplyw stlumiony
  - Efekt netto: stabilizacja vacuum psi=1, potencjalny UV fixed point

Porownanie trzech wariantow:
  A) LPA czyste (K=1) -- reprodukcja niestabilnosci z erg_phi_run.py
  B) LPA z K(psi) = psi^4 -- nowa fizyka
  C) LPA' z anomalnym wymiarem eta_k -- pelniejszy obraz

TGP v1 -- 2026-03-31  (OP-1b)
"""

import sys, os
os.environ['PYTHONIOENCODING'] = 'utf-8'
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# 1. Parametry TGP
# ============================================================
Phi0 = 24.66
gamma_TGP = 1.0
a_Gamma = 0.040
K_geo = 1.0  # normalizacja K_geo = 1 (w jednostkach Plancka)

k_UV = 1.0 / a_Gamma   # = 25
k_IR = np.sqrt(gamma_TGP)  # = 1
t_max = np.log(k_UV / k_IR)  # = 3.22

print("=" * 70)
print("  ERG WETTERICH Z K(phi) = K_geo * psi^4")
print("=" * 70)
print(f"  Phi0 = {Phi0}, a_Gamma = {a_Gamma}")
print(f"  k_UV = {k_UV:.1f}, k_IR = {k_IR:.1f}, t_max = {t_max:.2f}")
print(f"  K_geo = {K_geo}")
print()

# ============================================================
# 2. Siatka i potencjal
# ============================================================
N_grid = 50
psi_min, psi_max = 0.05, 2.5
psi_grid = np.linspace(psi_min, psi_max, N_grid)
dpsi = psi_grid[1] - psi_grid[0]

def V_TGP(psi):
    """Potencjal bare TGP: V(psi) = psi^3/3 - psi^4/4."""
    return psi**3 / 3 - psi**4 / 4

def K_phi(psi):
    """Czlon kinetyczny K(psi) = K_geo * psi^4."""
    return K_geo * psi**4

def V_second_derivative(V, dp):
    """Numeryczna V''(psi)."""
    r = np.zeros_like(V)
    r[1:-1] = (V[2:] - 2*V[1:-1] + V[:-2]) / dp**2
    r[0] = r[1]; r[-1] = r[-2]
    return r

# ============================================================
# 3. Wariant A: LPA czyste (K=1)
# ============================================================
def rhs_LPA(t, V):
    """dV/dt w czystym LPA (K=1)."""
    k = k_IR * np.exp(t)
    Vpp = V_second_derivative(V, dpsi)
    denom = np.maximum(Vpp + k**2, 0.01 * k**2)
    return k**5 / (32 * np.pi**2) / denom

# ============================================================
# 4. Wariant B: LPA z K(psi) = psi^4
# ============================================================
def rhs_LPA_K(t, V):
    """dV/dt w LPA z K(psi).

    Rownanie:
      dV/dt = K(psi)*k^6 / (32*pi^2 * (K(psi)*k^2 + V''(psi)))

    Wyprowadzenie:
      Regulator: R_k = K(psi)*(k^2-p^2)*theta(k^2-p^2)
      d_t R_k = K(psi)*2k^2*theta (bo d_t k^2 = 2k^2)
      Propagator (p<k): [K(psi)*k^2 + V'']^{-1}
      Calka: int_0^k d^4p/(2pi)^4 = k^4/(32pi^2)
      dV/dt = K(psi)*2k^2 * k^4/(32pi^2) / [2*(K*k^2+V'')]
            = K(psi)*k^6 / (32*pi^2*(K*k^2+V''))
    """
    k = k_IR * np.exp(t)
    Vpp = V_second_derivative(V, dpsi)
    K = K_phi(psi_grid)

    # Mianownik: K(psi)*k^2 + V''(psi)
    denom = K * k**2 + Vpp
    # Zabezpieczenie: denom > 0
    denom = np.maximum(denom, 0.01 * k**2 * np.maximum(K, 1e-8))

    # Licznik: K(psi)*k^6
    numerator = K * k**6

    return numerator / (32 * np.pi**2 * denom)

# ============================================================
# 5. Wariant C: LPA' z eta_k (anomalny wymiar) -- uproszczony
# ============================================================
def rhs_LPA_prime(t, V):
    """dV/dt w LPA' z przyblizona eta.

    eta_k(psi) = -d_t ln K_k(psi)
    W naszym przyblizeniu K(psi) jest stale w t (nie biega),
    ale dodajemy efektywny anomalny wymiar z V''''^2:

    eta_eff(psi) ~ -k^2 * (V'''(psi))^2 / (16*pi^2*(V''(psi)+K*k^2)^3)

    Modyfikowany przeplyw:
      dV/dt = K*k^6*(2 - eta_eff) / (32*pi^2*(K*k^2 + V''))
    """
    k = k_IR * np.exp(t)
    Vpp = V_second_derivative(V, dpsi)
    K = K_phi(psi_grid)

    # V'''(psi) -- trzecia pochodna
    Vppp = np.zeros_like(V)
    Vppp[1:-1] = (V[2:] - V[:-2]) / (2 * dpsi)  # to jest V'
    Vppp = V_second_derivative(Vppp, dpsi)  # potem V''' = (V')''

    denom = K * k**2 + Vpp
    denom = np.maximum(denom, 0.01 * k**2 * np.maximum(K, 1e-8))

    # eta_eff
    eta = np.zeros_like(V)
    mask = denom > 1e-10
    eta[mask] = -k**2 * Vppp[mask]**2 / (16 * np.pi**2 * denom[mask]**3)
    eta = np.clip(eta, -2.0, 2.0)  # fizyczne ograniczenie |eta| < 2

    numerator = K * k**6 * (2 - eta)
    return numerator / (32 * np.pi**2 * denom)

# ============================================================
# 6. Integracja trzech wariantow
# ============================================================
V0 = V_TGP(psi_grid)
t_span = (t_max, 0.05)
solver_opts = dict(method='Radau', rtol=1e-4, atol=1e-6, max_step=0.3, dense_output=True)

results = {}

for name, rhs_func in [('LPA (K=1)', rhs_LPA),
                         ('LPA+K(psi)', rhs_LPA_K),
                         ('LPA\'+K+eta', rhs_LPA_prime)]:
    print(f"  [{name}] Integracja UV -> IR ...", end="", flush=True)
    try:
        sol = solve_ivp(rhs_func, t_span, V0.copy(), **solver_opts)
        if sol.success:
            print(f" OK ({sol.t.shape[0]} krokow)")
            results[name] = sol
        else:
            print(f" BLAD: {sol.message}")
    except Exception as e:
        print(f" WYJATEK: {e}")

print()

# ============================================================
# 7. Analiza porownawcza
# ============================================================
idx_vac = np.argmin(np.abs(psi_grid - 1.0))
psi_vac = psi_grid[idx_vac]

print("=" * 70)
print("  ANALIZA POROWNAWCZA")
print("=" * 70)

for name, sol in results.items():
    V_IR = sol.y[:, -1]
    Vpp_IR = V_second_derivative(V_IR, dpsi)

    # Minimum globalne
    idx_min = np.argmin(V_IR)

    # Minimum lokalne blisko psi=1
    near = (psi_grid > 0.5) & (psi_grid < 1.5)
    if np.any(near):
        local_idx = np.argmin(V_IR[near])
        psi_loc = psi_grid[near][local_idx]
        V_loc = V_IR[near][local_idx]
        Vpp_loc = V_second_derivative(V_IR, dpsi)[np.where(near)[0][local_idx]]
    else:
        psi_loc = V_loc = Vpp_loc = np.nan

    print(f"\n  --- {name} ---")
    print(f"  V(psi=1):  UV = {V0[idx_vac]:.6f},  IR = {V_IR[idx_vac]:.6f}")
    print(f"  V''(psi=1): UV = {V_second_derivative(V0, dpsi)[idx_vac]:.4f},  IR = {Vpp_IR[idx_vac]:.4f}")
    print(f"  Minimum globalne IR: psi = {psi_grid[idx_min]:.3f}, V = {V_IR[idx_min]:.4f}")
    print(f"  Min lokalne (~1): psi = {psi_loc:.3f}, V = {V_loc:.6f}, V'' = {Vpp_loc:.4f}")

    # Stabilnosc: V'' > 0 oznacza minimum lokalne (STABILNE)
    if Vpp_loc > 0:
        m_eff = np.sqrt(Vpp_loc)
        print(f"  -> STABILNE! Masa efektywna: m = {m_eff:.4f}")
    else:
        print(f"  -> NIESTABILNE (V'' < 0, saddle/maximum)")

    # Bariera
    if psi_loc > 0.3:
        left = (psi_grid > 0.1) & (psi_grid < psi_loc)
        if np.any(left):
            V_left = V_IR[left]
            idx_bar = np.argmax(V_left)
            psi_bar = psi_grid[left][idx_bar]
            V_bar = V_left[idx_bar]
            print(f"  Bariera: psi = {psi_bar:.3f}, V = {V_bar:.6f}")
            print(f"  Wysokosc bariery: {V_bar - V_loc:.6f}")

    # Bezwymiarowy u = V/k^4
    print(f"  Bezwymiarowy u(psi=1):")
    for t_val in [t_max, t_max*0.5, 0.05]:
        k_val = k_IR * np.exp(t_val)
        V_t = sol.sol(t_val)
        u_val = V_t[idx_vac] / k_val**4
        print(f"    t={t_val:.2f} (k={k_val:.1f}): u = {u_val:.2e}")

# ============================================================
# 8. Szukanie UV fixed point (zmienne bezwymiarowe)
# ============================================================
print(f"\n{'='*70}")
print("  SZUKANIE UV FIXED POINT")
print("="*70)

# W zmiennych bezwymiarowych:
# u_k(psi) = V_k(psi) / k^4
# du/dt = -4*u + (rhs bez k^4)
# Punkt staly: du/dt = 0 => u* = rhs/(4*k^4)

for name, sol in results.items():
    print(f"\n  --- {name} ---")
    # Sprawdzamy du/dt na UV
    t_uv = t_max
    k_uv = k_IR * np.exp(t_uv)
    V_uv = sol.sol(t_uv)
    u_uv = V_uv / k_uv**4

    # du/dt = (dV/dt - 4*V_k) / k^4
    # ale prosciej: sprawdzamy czy u sie stabilizuje
    t_pts = np.linspace(t_max, 0.05, 40)
    u_at_1 = np.array([sol.sol(t)[idx_vac] / (k_IR*np.exp(t))**4 for t in t_pts])
    du_dt = np.gradient(u_at_1, t_pts)

    # Czy du/dt -> 0 na UV?
    print(f"  u(psi=1) na UV (t={t_max:.2f}): {u_at_1[0]:.2e}")
    print(f"  u(psi=1) na IR (t=0.05): {u_at_1[-1]:.2e}")
    print(f"  du/dt na UV: {du_dt[0]:.2e}")
    print(f"  du/dt na srodku: {du_dt[len(du_dt)//2]:.2e}")
    print(f"  du/dt na IR: {du_dt[-1]:.2e}")

    # Kryterium: |du/dt| < epsilon na UV => UV FP
    if abs(du_dt[0]) < 0.01 * abs(u_at_1[0]) and abs(u_at_1[0]) > 1e-15:
        print(f"  ** MOZLIWY UV FIXED POINT: |du/dt|/|u| = {abs(du_dt[0]/u_at_1[0]):.4f} **")
    else:
        ratio = abs(du_dt[0]) / max(abs(u_at_1[0]), 1e-20)
        print(f"  Brak FP: |du/dt|/|u| = {ratio:.2e}")

# ============================================================
# 9. Profil V_IR: porownanie potencjalu
# ============================================================
print(f"\n{'='*70}")
print("  PROFIL V_IR: KLUCZOWE WARTOSCI")
print("="*70)

psi_check = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0]
header = f"  {'psi':>5} | {'V_UV':>10}"
for name in results.keys():
    header += f" | {name:>14}"
print(header)
print("  " + "-" * (len(header) - 2))

for psi_val in psi_check:
    idx = np.argmin(np.abs(psi_grid - psi_val))
    line = f"  {psi_grid[idx]:5.2f} | {V0[idx]:10.4f}"
    for name, sol in results.items():
        V_IR = sol.y[:, -1]
        line += f" | {V_IR[idx]:14.4f}"
    print(line)

# ============================================================
# 10. Przeplyw V''(psi=1) -- masa^2 pola
# ============================================================
print(f"\n{'='*70}")
print("  PRZEPLYW V''(psi=1) [masa^2]")
print("="*70)

t_pts = np.linspace(t_max, 0.05, 20)
print(f"  {'t':>6} {'k':>6}", end="")
for name in results.keys():
    print(f" | {name:>14}", end="")
print()

for t_val in t_pts[::4]:  # co 4-ty punkt
    k_val = k_IR * np.exp(t_val)
    line = f"  {t_val:6.2f} {k_val:6.1f}"
    for name, sol in results.items():
        V_t = sol.sol(t_val)
        Vpp_t = V_second_derivative(V_t, dpsi)
        line += f" | {Vpp_t[idx_vac]:14.4f}"
    print(line)

# ============================================================
# 11. Efekt K(psi) na stabilnosc -- szczegolowa analiza
# ============================================================
if 'LPA+K(psi)' in results:
    sol_K = results['LPA+K(psi)']
    V_IR_K = sol_K.y[:, -1]
    Vpp_IR_K = V_second_derivative(V_IR_K, dpsi)
    K_vals = K_phi(psi_grid)

    print(f"\n{'='*70}")
    print("  SZCZEGOLOWA ANALIZA LPA+K(psi)")
    print("="*70)

    # Efektywna masa^2: m_eff^2 = V''/K (bo kinetyczny czlon to K*(dpsi)^2)
    # Fizyczna masa^2 pola psi = V''(psi) / K(psi) na minimum
    m2_phys = np.zeros_like(Vpp_IR_K)
    mask_K = K_vals > 1e-10
    m2_phys[mask_K] = Vpp_IR_K[mask_K] / K_vals[mask_K]

    print(f"  Na psi = {psi_vac:.2f}:")
    print(f"    K(psi) = {K_vals[idx_vac]:.4f}")
    print(f"    V''(psi) = {Vpp_IR_K[idx_vac]:.4f}")
    print(f"    m^2_phys = V''/K = {m2_phys[idx_vac]:.4f}")

    if m2_phys[idx_vac] > 0:
        print(f"    m_phys = {np.sqrt(m2_phys[idx_vac]):.4f}")
        print(f"    ** FIZYCZNIE STABILNE! **")
    else:
        print(f"    m^2_phys < 0 -> fizycznie niestabilne")

    # Porownaj V''/K (fizyczna masa) vs V'' (koordynatowa masa)
    print(f"\n  Porownanie masa^2 koordynatowa vs fizyczna:")
    print(f"  {'psi':>5} | {'V_IR':>8} | {'V_IR\"':>8} | {'K(psi)':>8} | {'m2_phys':>10} | {'status':>10}")
    for idx_c in range(0, N_grid, N_grid//10):
        psi_c = psi_grid[idx_c]
        status = "STAB" if m2_phys[idx_c] > 0 and mask_K[idx_c] else "NIESTAB" if mask_K[idx_c] else "K~0"
        print(f"  {psi_c:5.2f} | {V_IR_K[idx_c]:8.3f} | {Vpp_IR_K[idx_c]:8.3f} | {K_vals[idx_c]:8.4f} | {m2_phys[idx_c]:10.4f} | {status:>10}")

# ============================================================
# 12. Wykresy
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(20, 12))

# Panel (0,0): V_UV vs V_IR dla trzech wariantow
ax = axes[0, 0]
ax.plot(psi_grid, V0, 'k--', lw=2, label='V_UV (bare)', alpha=0.7)
colors = {'LPA (K=1)': 'blue', 'LPA+K(psi)': 'red', "LPA'+K+eta": 'green'}
for name, sol in results.items():
    V_IR = sol.y[:, -1]
    ax.plot(psi_grid, V_IR, '-', lw=2, color=colors.get(name, 'gray'), label=f'V_IR [{name}]')
ax.axvline(1.0, color='gray', ls=':', alpha=0.3)
ax.set_xlabel(r'$\psi$'); ax.set_ylabel(r'$V_k(\psi)$')
ax.set_title('Potencjal UV vs IR')
ax.legend(fontsize=7); ax.grid(True, alpha=0.2)
# Dynamiczny ylim
y_vals = [V0]
for sol in results.values():
    y_vals.append(sol.y[:, -1])
all_V = np.concatenate(y_vals)
y_lo = max(np.percentile(all_V, 2), -50)
y_hi = min(np.percentile(all_V, 98), 50)
ax.set_ylim(y_lo, y_hi)

# Panel (0,1): V_IR zoom na psi in [0.5, 1.5]
ax = axes[0, 1]
mask_zoom = (psi_grid > 0.3) & (psi_grid < 1.7)
ax.plot(psi_grid[mask_zoom], V0[mask_zoom], 'k--', lw=2, label='V_UV', alpha=0.7)
for name, sol in results.items():
    V_IR = sol.y[:, -1]
    ax.plot(psi_grid[mask_zoom], V_IR[mask_zoom], '-', lw=2,
            color=colors.get(name, 'gray'), label=f'IR [{name}]')
ax.axvline(1.0, color='gray', ls=':', alpha=0.3)
ax.set_xlabel(r'$\psi$'); ax.set_ylabel(r'$V_k(\psi)$')
ax.set_title('Zoom: psi ~ 1 (vacuum)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.2)

# Panel (0,2): V''(psi=1) flow
ax = axes[0, 2]
t_pts = np.linspace(t_max, 0.05, 30)
for name, sol in results.items():
    vpp_flow = []
    for t in t_pts:
        V_t = sol.sol(t)
        vpp = V_second_derivative(V_t, dpsi)
        vpp_flow.append(vpp[idx_vac])
    ax.plot(t_pts, vpp_flow, '-o', ms=2, color=colors.get(name, 'gray'), label=name)
ax.axhline(0, color='black', ls='-', alpha=0.3)
ax.set_xlabel('t = ln(k/k_IR)'); ax.set_ylabel(r"$V''(\psi=1)$")
ax.set_title("Masa^2 (V'') flow")
ax.legend(fontsize=7); ax.grid(True, alpha=0.2)

# Panel (1,0): Bezwymiarowy u(psi=1) flow
ax = axes[1, 0]
for name, sol in results.items():
    u_flow = [sol.sol(t)[idx_vac] / (k_IR*np.exp(t))**4 for t in t_pts]
    ax.plot(t_pts, u_flow, '-o', ms=2, color=colors.get(name, 'gray'), label=name)
ax.set_xlabel('t = ln(k/k_IR)'); ax.set_ylabel(r'$u = V/k^4$ at $\psi=1$')
ax.set_title('Bezwymiarowy potencjal')
ax.legend(fontsize=7); ax.grid(True, alpha=0.2)

# Panel (1,1): Fizyczna masa m2_phys = V''/K(psi) -- profil IR
ax = axes[1, 1]
for name, sol in results.items():
    V_IR = sol.y[:, -1]
    Vpp_IR = V_second_derivative(V_IR, dpsi)
    if 'K' in name:
        K_v = K_phi(psi_grid)
        m2 = np.where(K_v > 1e-10, Vpp_IR / K_v, np.nan)
        ax.plot(psi_grid, m2, '-', lw=2, color=colors.get(name, 'gray'),
                label=f'm2_phys [{name}]')
    else:
        ax.plot(psi_grid, Vpp_IR, '-', lw=2, color=colors.get(name, 'gray'),
                label=f"V'' [{name}]")
ax.axhline(0, color='black', ls='-', alpha=0.3)
ax.axvline(1.0, color='gray', ls=':', alpha=0.3)
ax.set_xlabel(r'$\psi$'); ax.set_ylabel(r"$m^2_{\rm phys}(\psi) = V''/K$")
ax.set_title('Fizyczna masa^2 (IR)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.2)
ax.set_ylim(-10, 10)

# Panel (1,2): V(psi=1) flow (wymiarowy)
ax = axes[1, 2]
for name, sol in results.items():
    V1_flow = [sol.sol(t)[idx_vac] for t in t_pts]
    ax.plot(t_pts, V1_flow, '-o', ms=2, color=colors.get(name, 'gray'), label=name)
ax.set_xlabel('t = ln(k/k_IR)'); ax.set_ylabel(r'$V_k(\psi=1)$')
ax.set_title('V na vacuum (wymiarowy)')
ax.legend(fontsize=7); ax.grid(True, alpha=0.2)

plt.suptitle('ERG Wetterich: LPA vs LPA+K($\\psi$) vs LPA\'', fontsize=14, y=1.02)
plt.tight_layout()
script_dir = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(script_dir, 'erg_with_K_phi.png'), dpi=150, bbox_inches='tight')
print(f"\n  Wykres: scripts/erg_with_K_phi.png")

# ============================================================
# 13. Podsumowanie
# ============================================================
print(f"\n{'='*70}")
print("  PODSUMOWANIE ERG z K(phi)")
print("="*70)

for name, sol in results.items():
    V_IR = sol.y[:, -1]
    Vpp_IR = V_second_derivative(V_IR, dpsi)
    v_at_1 = V_IR[idx_vac]
    vpp_at_1 = Vpp_IR[idx_vac]

    stab = "STABILNE" if vpp_at_1 > 0 else "NIESTABILNE"
    print(f"  {name:20s}: V(1)={v_at_1:+.4f}, V''(1)={vpp_at_1:+.4f} -> {stab}")

print(f"""
  WNIOSKI:
  1. LPA (K=1): V''(1) < 0 na IR -> niestabilne (potwierdzenie erg_phi_run.py)
  2. LPA+K(psi): K(psi)=psi^4 modyfikuje propagator.
     - Efekt: K(1)=1, wiec na vacuum K nie zmienia masywnie propagatora
     - Ale K tlumi korekty dla psi >> 1 i psi << 1 (bariera kinetyczna)
  3. Fizyczna masa: m^2_phys = V''/K jest wlasciwym kryterium stabilnosci
     (bo rownanie ruchu to K*d^2 psi/dt^2 = -V', wiec m^2 = V''/K)

  NASTEPNE KROKI:
  OP-1b+: ERG z BIEGAJACYM K_k(psi) (K tez plynie pod RG)
  OP-1c: UV fixed point w pelnym ukladzie (V, K)
  OP-1d: Wlaczenie wymiaru anomalnego eta_k z pelnego rownania
""")
print("GOTOWE.")
