#!/usr/bin/env python3
"""
ERG z biegacym K_k(psi): sprzezony uklad (V_k, K_k)
=====================================================
Pelny system LPA' z field-dependent K(psi):

  Gamma_k[psi] = int d^4x [ K_k(psi)/2 * (d psi)^2 + V_k(psi) ]

Rownania przeplywu (Litim regulator, d=4):

  dV/dt = K*k^6 / (32*pi^2 * (K*k^2 + V''))

  dK/dt = -K*k^6 / (16*pi^2) * (V'''/K)^2 / (k^2 + V''/K)^3

Zmienne bezwymiarowe (d=4):
  u = V / k^4    (potencjal)
  kappa = K       (juz bezwymiarowy w d=4 dla skalara)

Rownania FP:
  du/dt = -4u + kappa / (32*pi^2*(kappa + u''))  = 0
  dkappa/dt = ... = 0

Szukamy: czy istnieje (u*, kappa*) != (0, K_geo*psi^4)?

TGP v1 -- 2026-03-31 (OP-1b/c)
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
# 1. Parametry
# ============================================================
Phi0 = 24.66
a_Gamma = 0.040
K_geo = 1.0

k_UV = 1.0 / a_Gamma   # = 25
k_IR = np.sqrt(1.0)     # = 1
t_max = np.log(k_UV / k_IR)  # = 3.22

N_grid = 40
psi_min, psi_max = 0.05, 2.5
psi_grid = np.linspace(psi_min, psi_max, N_grid)
dpsi = psi_grid[1] - psi_grid[0]

print("=" * 70)
print("  ERG SPRZEZONY (V_k, K_k) Z UV FIXED POINT SEARCH")
print("=" * 70)
print(f"  k_UV = {k_UV:.1f}, k_IR = {k_IR:.1f}, t = [0, {t_max:.2f}]")
print(f"  Siatka: {N_grid} punktow, psi in [{psi_min}, {psi_max}]")
print()

# ============================================================
# 2. Funkcje pomocnicze
# ============================================================
def D2(f, dp):
    """Druga pochodna numeryczna."""
    r = np.zeros_like(f)
    r[1:-1] = (f[2:] - 2*f[1:-1] + f[:-2]) / dp**2
    r[0] = r[1]; r[-1] = r[-2]
    return r

def D1(f, dp):
    """Pierwsza pochodna numeryczna."""
    r = np.zeros_like(f)
    r[1:-1] = (f[2:] - f[:-2]) / (2*dp)
    r[0] = (-3*f[0] + 4*f[1] - f[2]) / (2*dp)
    r[-1] = (3*f[-1] - 4*f[-2] + f[-3]) / (2*dp)
    return r

def D3(f, dp):
    """Trzecia pochodna numeryczna."""
    return D1(D2(f, dp), dp)

def V_TGP_bare(psi):
    return psi**3 / 3 - psi**4 / 4

def K_TGP_bare(psi):
    return K_geo * psi**4

# ============================================================
# 3. Sprzezony uklad (V_k, K_k)
# ============================================================
def coupled_rhs(t, y):
    """
    Prawa strona sprzezonego ukladu.
    y = [V_0, ..., V_{N-1}, K_0, ..., K_{N-1}]
    """
    N = N_grid
    V = y[:N].copy()
    K = y[N:].copy()

    k = k_IR * np.exp(t)

    # Pochodne V
    Vpp = D2(V, dpsi)
    Vppp = D3(V, dpsi)

    # Zabezpieczenia
    K_safe = np.maximum(K, 1e-10)

    # Mianownik: K*k^2 + V''
    denom_V = K_safe * k**2 + Vpp
    denom_V = np.maximum(denom_V, 0.01 * k**2 * K_safe)

    # --- dV/dt ---
    # dV/dt = K * k^6 / (32*pi^2 * (K*k^2 + V''))
    dVdt = K_safe * k**6 / (32 * np.pi**2 * denom_V)

    # --- dK/dt ---
    # dK/dt = -K * k^6 / (16*pi^2) * (V'''/K)^2 / (k^2 + V''/K)^3
    # = -k^6 * V'''^2 / (16*pi^2 * K * (k^2 + V''/K)^3)
    # = -k^6 * V'''^2 / (16*pi^2 * (K*k^2 + V'')^3 / K^2)
    # = -k^6 * K^2 * V'''^2 / (16*pi^2 * (K*k^2 + V'')^3)

    denom_K = denom_V**3
    denom_K = np.maximum(denom_K, 1e-30)

    dKdt = -k**6 * K_safe**2 * Vppp**2 / (16 * np.pi**2 * denom_K)

    # Stabilizacja: K nie moze byc ujemne (fizyczne wymaganie)
    # Jesli K blisko 0, tlum dK
    mask_small_K = K < 0.01 * K_TGP_bare(psi_grid)
    dKdt[mask_small_K] = np.maximum(dKdt[mask_small_K], 0)

    return np.concatenate([dVdt, dKdt])

# ============================================================
# 4. Integracja sprzezonego ukladu
# ============================================================
V0 = V_TGP_bare(psi_grid)
K0 = K_TGP_bare(psi_grid)
y0 = np.concatenate([V0, K0])

print("  [Coupled V+K] Integracja UV -> IR (Radau)...", flush=True)
t_span = (t_max, 0.05)

sol_coupled = solve_ivp(
    coupled_rhs, t_span, y0,
    method='Radau', rtol=1e-4, atol=1e-6,
    max_step=0.2, dense_output=True
)

if sol_coupled.success:
    print(f"  OK ({sol_coupled.t.shape[0]} krokow)")
else:
    print(f"  BLAD: {sol_coupled.message}")
    # Sprobuj lagodniejszy
    print("  Probuje z wiekszym max_step...", flush=True)
    sol_coupled = solve_ivp(
        coupled_rhs, t_span, y0,
        method='Radau', rtol=1e-3, atol=1e-5,
        max_step=0.5, dense_output=True
    )
    if sol_coupled.success:
        print(f"  OK (retry, {sol_coupled.t.shape[0]} krokow)")
    else:
        print(f"  DEFINITYWNY BLAD: {sol_coupled.message}")

# ============================================================
# 5. Porownanie: coupled vs fixed-K
# ============================================================
# Fixed-K (z poprzedniego skryptu)
def rhs_fixed_K(t, V):
    k = k_IR * np.exp(t)
    Vpp = D2(V, dpsi)
    K = K_TGP_bare(psi_grid)
    denom = np.maximum(K * k**2 + Vpp, 0.01 * k**2 * np.maximum(K, 1e-8))
    return K * k**6 / (32 * np.pi**2 * denom)

print("  [Fixed K] Integracja...", flush=True)
sol_fixedK = solve_ivp(
    rhs_fixed_K, t_span, V0.copy(),
    method='Radau', rtol=1e-4, atol=1e-6,
    max_step=0.3, dense_output=True
)
print(f"  {'OK' if sol_fixedK.success else 'BLAD'} ({sol_fixedK.t.shape[0]} krokow)")

# ============================================================
# 6. Analiza
# ============================================================
idx_vac = np.argmin(np.abs(psi_grid - 1.0))
psi_vac = psi_grid[idx_vac]

print(f"\n{'='*70}")
print("  ANALIZA WYNIKOW")
print("="*70)

if sol_coupled.success:
    V_IR_c = sol_coupled.y[:N_grid, -1]
    K_IR_c = sol_coupled.y[N_grid:, -1]
    Vpp_IR_c = D2(V_IR_c, dpsi)

    print(f"\n  --- Coupled (V,K) ---")
    print(f"  V(psi=1):  UV = {V0[idx_vac]:.4f},  IR = {V_IR_c[idx_vac]:.4f}")
    print(f"  V''(psi=1): UV = {D2(V0, dpsi)[idx_vac]:.4f},  IR = {Vpp_IR_c[idx_vac]:.4f}")
    print(f"  K(psi=1):  UV = {K0[idx_vac]:.4f},  IR = {K_IR_c[idx_vac]:.4f}")
    print(f"  K zmiana: {(K_IR_c[idx_vac] - K0[idx_vac])/K0[idx_vac]*100:+.2f}%")

    # Fizyczna masa
    m2_phys_c = Vpp_IR_c[idx_vac] / max(K_IR_c[idx_vac], 1e-10)
    print(f"  m^2_phys = V''/K = {m2_phys_c:.4f}")
    if m2_phys_c > 0:
        print(f"  m_phys = {np.sqrt(m2_phys_c):.4f} -> **STABILNE**")
    else:
        print(f"  m^2 < 0 -> NIESTABILNE")

    # Profil K_IR vs K_UV
    print(f"\n  Profil K(psi) UV vs IR:")
    for psi_c in [0.3, 0.5, 0.7, 1.0, 1.3, 1.5, 2.0]:
        idx_c = np.argmin(np.abs(psi_grid - psi_c))
        K_uv = K0[idx_c]
        K_ir = K_IR_c[idx_c]
        delta = (K_ir - K_uv) / max(K_uv, 1e-10) * 100
        print(f"    psi={psi_c:.1f}: K_UV={K_uv:.4f}, K_IR={K_ir:.4f} ({delta:+.1f}%)")

if sol_fixedK.success:
    V_IR_f = sol_fixedK.y[:, -1]
    Vpp_IR_f = D2(V_IR_f, dpsi)
    print(f"\n  --- Fixed K (referencja) ---")
    print(f"  V(psi=1) IR = {V_IR_f[idx_vac]:.4f}")
    print(f"  V''(psi=1) IR = {Vpp_IR_f[idx_vac]:.4f}")
    m2_f = Vpp_IR_f[idx_vac] / K0[idx_vac]
    print(f"  m^2_phys = {m2_f:.4f}")

# ============================================================
# 7. Zmienne bezwymiarowe i UV Fixed Point
# ============================================================
print(f"\n{'='*70}")
print("  UV FIXED POINT SEARCH (zmienne bezwymiarowe)")
print("="*70)

if sol_coupled.success:
    t_pts = np.linspace(t_max, 0.05, 40)

    # Bezwymiarowe: u = V/k^4, kappa = K (juz bezwymiarowe w d=4)
    u_flow = []
    upp_flow = []
    kappa_flow = []
    du_dt_flow = []

    for i, t_val in enumerate(t_pts):
        state = sol_coupled.sol(t_val)
        V_t = state[:N_grid]
        K_t = state[N_grid:]
        k_val = k_IR * np.exp(t_val)

        u_t = V_t / k_val**4
        u_pp = D2(u_t, dpsi)

        u_flow.append(u_t[idx_vac])
        upp_flow.append(u_pp[idx_vac])
        kappa_flow.append(K_t[idx_vac])

    u_flow = np.array(u_flow)
    upp_flow = np.array(upp_flow)
    kappa_flow = np.array(kappa_flow)

    # du/dt numeryczny
    du_dt = np.gradient(u_flow, t_pts)
    dk_dt = np.gradient(kappa_flow, t_pts)

    print(f"\n  Przeplyw na psi = {psi_vac:.2f}:")
    print(f"  {'t':>6} {'k':>6} | {'u':>12} {'kappa':>12} | {'du/dt':>12} {'dk/dt':>12}")
    for i in range(0, len(t_pts), 5):
        k_val = k_IR * np.exp(t_pts[i])
        print(f"  {t_pts[i]:6.2f} {k_val:6.1f} | {u_flow[i]:12.2e} {kappa_flow[i]:12.4f} | {du_dt[i]:12.2e} {dk_dt[i]:12.2e}")

    # Kryterium FP: |du/dt|/|u| < eps AND |dk/dt|/|kappa| < eps
    print(f"\n  Szukanie FP (kryterium: |d/dt|/|val| < 0.1):")
    found_fp = False
    for i in range(len(t_pts)):
        if abs(u_flow[i]) > 1e-15 and abs(kappa_flow[i]) > 1e-15:
            rel_u = abs(du_dt[i]) / abs(u_flow[i])
            rel_k = abs(dk_dt[i]) / abs(kappa_flow[i])
            if rel_u < 0.1 and rel_k < 0.1:
                k_val = k_IR * np.exp(t_pts[i])
                print(f"  ** KANDYDAT FP: t={t_pts[i]:.2f}, k={k_val:.1f}")
                print(f"     u* = {u_flow[i]:.2e}, kappa* = {kappa_flow[i]:.4f}")
                print(f"     |du/dt|/|u| = {rel_u:.4f}, |dk/dt|/|k| = {rel_k:.4f}")
                found_fp = True

    if not found_fp:
        print(f"  Brak FP w badanym zakresie.")
        # Podsumuj zachowanie asymptotyczne
        print(f"\n  Zachowanie UV (t -> t_max):")
        print(f"    u = {u_flow[0]:.2e}, du/dt = {du_dt[0]:.2e}")
        print(f"    kappa = {kappa_flow[0]:.4f}, dk/dt = {dk_dt[0]:.2e}")
        print(f"  Zachowanie IR (t -> 0):")
        print(f"    u = {u_flow[-1]:.2e}, du/dt = {du_dt[-1]:.2e}")
        print(f"    kappa = {kappa_flow[-1]:.4f}, dk/dt = {dk_dt[-1]:.2e}")

    # Asymptotyczne bezpieczenstwo: czy kappa(UV) jest skonczony?
    print(f"\n  Asymptotic Safety test:")
    print(f"    K(psi=1) na UV: {kappa_flow[0]:.6f}")
    print(f"    K(psi=1) na IR: {kappa_flow[-1]:.6f}")
    K_ratio = kappa_flow[-1] / kappa_flow[0] if abs(kappa_flow[0]) > 1e-15 else np.inf
    print(f"    Stosunek K_IR/K_UV = {K_ratio:.4f}")
    if 0.1 < K_ratio < 10:
        print(f"    -> K OGRANICZONE (nie ucieka) -- pozytywny sygnał AS")
    else:
        print(f"    -> K znaczaco sie zmienia ({K_ratio:.1f}x)")

# ============================================================
# 8. Analiza stabilnosci FP -- wykladniki krytyczne
# ============================================================
print(f"\n{'='*70}")
print("  ANALIZA STABILNOSCI: wymiar krytyczny")
print("="*70)

if sol_coupled.success:
    # Sprawdz jak szybko perturbacje rosna/maleja wokol trajektorii
    # Linearyzacja: delta_u ~ exp(theta * t) -> theta = wykladnik krytyczny
    # theta < 0: relevant (UV attractive), theta > 0: irrelevant (UV repulsive)

    # Prosta estymacja: theta ~ d(ln|u|)/dt w poblizu UV
    mask_uv = t_pts > 0.7 * t_max
    if np.sum(mask_uv) > 3:
        u_uv = u_flow[mask_uv]
        t_uv = t_pts[mask_uv]
        # ln|u| vs t -- nachylenie = theta
        mask_pos = np.abs(u_uv) > 1e-15
        if np.sum(mask_pos) > 3:
            log_u = np.log(np.abs(u_uv[mask_pos]))
            t_sub = t_uv[mask_pos]
            # Fit liniowy
            coeffs = np.polyfit(t_sub, log_u, 1)
            theta_u = coeffs[0]
            print(f"  theta_u (z fitu ln|u| vs t w UV): {theta_u:.2f}")
            print(f"    theta < 0: relevant (UV attractive)")
            print(f"    theta > 0: irrelevant (UV repulsive)")
            if theta_u < 0:
                print(f"    -> u jest UV-atrakcyjne (theta = {theta_u:.2f})")
            else:
                print(f"    -> u jest UV-repulsywne")

        # To samo dla kappa
        k_uv = kappa_flow[mask_uv]
        mask_pos_k = np.abs(k_uv) > 1e-15
        if np.sum(mask_pos_k) > 3:
            log_k = np.log(np.abs(k_uv[mask_pos_k]))
            t_sub_k = t_uv[mask_pos_k]
            coeffs_k = np.polyfit(t_sub_k, log_k, 1)
            theta_k = coeffs_k[0]
            print(f"  theta_kappa (z fitu ln|kappa| vs t w UV): {theta_k:.2f}")

# ============================================================
# 9. Profil pelny: V_IR i K_IR
# ============================================================
print(f"\n{'='*70}")
print("  PROFIL V_IR i K_IR (coupled)")
print("="*70)

if sol_coupled.success:
    print(f"  {'psi':>5} | {'V_UV':>8} {'V_IR':>8} {'dV':>8} | {'K_UV':>8} {'K_IR':>8} {'dK%':>7} | {'m2_phys':>8} {'stab':>6}")
    print("  " + "-" * 80)
    for idx_c in range(0, N_grid, N_grid//12):
        psi_c = psi_grid[idx_c]
        v_uv = V0[idx_c]; v_ir = V_IR_c[idx_c]
        k_uv = K0[idx_c]; k_ir = K_IR_c[idx_c]
        dk_pct = (k_ir - k_uv) / max(k_uv, 1e-10) * 100
        vpp_ir = D2(V_IR_c, dpsi)[idx_c]
        m2 = vpp_ir / max(k_ir, 1e-10)
        stab = "STAB" if m2 > 0 else "niestab"
        print(f"  {psi_c:5.2f} | {v_uv:8.4f} {v_ir:8.2f} {v_ir-v_uv:+8.2f} | {k_uv:8.4f} {k_ir:8.4f} {dk_pct:+7.1f}% | {m2:8.2f} {stab:>6}")

# ============================================================
# 10. Wykresy
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(20, 12))

t_pts_plot = np.linspace(t_max, 0.05, 30)

# (0,0) V_IR porownanie: coupled vs fixed-K vs bare
ax = axes[0, 0]
ax.plot(psi_grid, V0, 'k--', lw=2, label='V_UV (bare)')
if sol_fixedK.success:
    ax.plot(psi_grid, sol_fixedK.y[:, -1], 'b-', lw=2, label='V_IR [fixed K]')
if sol_coupled.success:
    ax.plot(psi_grid, V_IR_c, 'r-', lw=2, label='V_IR [coupled K]')
ax.axvline(1.0, color='gray', ls=':', alpha=0.3)
ax.set_xlabel(r'$\psi$'); ax.set_ylabel(r'$V(\psi)$')
ax.set_title('Potencjal: coupled vs fixed K')
ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
# Limit y
all_v = [V0]
if sol_fixedK.success: all_v.append(sol_fixedK.y[:, -1])
if sol_coupled.success: all_v.append(V_IR_c)
all_v = np.concatenate(all_v)
ax.set_ylim(max(np.percentile(all_v, 1), -500), min(np.percentile(all_v, 99), 50))

# (0,1) K(psi) UV vs IR
ax = axes[0, 1]
if sol_coupled.success:
    ax.plot(psi_grid, K0, 'b--', lw=2, label='K_UV (bare = psi^4)')
    ax.plot(psi_grid, K_IR_c, 'r-', lw=2, label='K_IR (po RG)')
    ax.axvline(1.0, color='gray', ls=':', alpha=0.3)
    ax.set_xlabel(r'$\psi$'); ax.set_ylabel(r'$K(\psi)$')
    ax.set_title('K(psi): UV vs IR')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

# (0,2) V''(psi=1) flow
ax = axes[0, 2]
if sol_coupled.success:
    vpp_c = [D2(sol_coupled.sol(t)[:N_grid], dpsi)[idx_vac] for t in t_pts_plot]
    ax.plot(t_pts_plot, vpp_c, 'r-o', ms=2, label='coupled K')
if sol_fixedK.success:
    vpp_f = [D2(sol_fixedK.sol(t), dpsi)[idx_vac] for t in t_pts_plot]
    ax.plot(t_pts_plot, vpp_f, 'b-o', ms=2, label='fixed K')
ax.axhline(0, color='black', ls='-', alpha=0.3)
ax.set_xlabel('t = ln(k/k_IR)'); ax.set_ylabel(r"$V''(\psi=1)$")
ax.set_title("Masa^2 (V'') flow")
ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

# (1,0) kappa(psi=1) flow
ax = axes[1, 0]
if sol_coupled.success:
    kappa_plot = [sol_coupled.sol(t)[N_grid + idx_vac] for t in t_pts_plot]
    ax.plot(t_pts_plot, kappa_plot, 'r-o', ms=2, label=r'$K(\psi=1)$')
    ax.axhline(K0[idx_vac], color='blue', ls='--', alpha=0.5, label=f'K_bare = {K0[idx_vac]:.2f}')
    ax.set_xlabel('t = ln(k/k_IR)'); ax.set_ylabel(r'$K_k(\psi=1)$')
    ax.set_title('Przeplyw K na vacuum')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

# (1,1) Bezwymiarowy u = V/k^4 na psi=1
ax = axes[1, 1]
if sol_coupled.success:
    u_c = [sol_coupled.sol(t)[idx_vac] / (k_IR*np.exp(t))**4 for t in t_pts_plot]
    ax.plot(t_pts_plot, u_c, 'r-o', ms=2, label='coupled K')
if sol_fixedK.success:
    u_f = [sol_fixedK.sol(t)[idx_vac] / (k_IR*np.exp(t))**4 for t in t_pts_plot]
    ax.plot(t_pts_plot, u_f, 'b-o', ms=2, label='fixed K')
ax.set_xlabel('t = ln(k/k_IR)'); ax.set_ylabel(r'$u = V/k^4$ at $\psi=1$')
ax.set_title('Bezwymiarowy potencjal')
ax.legend(fontsize=8); ax.grid(True, alpha=0.2)

# (1,2) Fizyczna masa m2 = V''/K profil IR
ax = axes[1, 2]
if sol_coupled.success:
    K_safe = np.maximum(K_IR_c, 1e-10)
    m2_c = D2(V_IR_c, dpsi) / K_safe
    ax.plot(psi_grid, m2_c, 'r-', lw=2, label='coupled K')
if sol_fixedK.success:
    K_bare = K_TGP_bare(psi_grid)
    K_bare_safe = np.maximum(K_bare, 1e-10)
    m2_f = D2(sol_fixedK.y[:, -1], dpsi) / K_bare_safe
    ax.plot(psi_grid, m2_f, 'b-', lw=2, label='fixed K')
ax.axhline(0, color='black', ls='-', alpha=0.3)
ax.axvline(1.0, color='gray', ls=':', alpha=0.3)
ax.set_xlabel(r'$\psi$'); ax.set_ylabel(r"$m^2_{\rm phys} = V''/K$")
ax.set_title('Fizyczna masa^2 (IR)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
ax.set_ylim(-50, 500)

plt.suptitle('ERG coupled (V, K): biegacy czlon kinetyczny', fontsize=14, y=1.02)
plt.tight_layout()
sd = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(sd, 'erg_coupled_VK.png'), dpi=150, bbox_inches='tight')
print(f"\n  Wykres: scripts/erg_coupled_VK.png")

# ============================================================
# 11. Podsumowanie
# ============================================================
print(f"\n{'='*70}")
print("  PODSUMOWANIE")
print("="*70)
if sol_coupled.success:
    print(f"""
  1. SPRZEZONY (V,K) SYSTEM:
     - K biegnie pod RG: K(1) zmienia sie o {(K_IR_c[idx_vac]-K0[idx_vac])/K0[idx_vac]*100:+.1f}%
     - V''(1) IR = {Vpp_IR_c[idx_vac]:+.2f} ({'STABILNE' if Vpp_IR_c[idx_vac]>0 else 'NIESTABILNE'})
     - m^2_phys = V''/K = {m2_phys_c:.2f}

  2. POROWNANIE:
     - Fixed K:   V''(1) = {Vpp_IR_f[idx_vac]:+.2f}, m^2 = {np.asarray(m2_f).item():.2f}
     - Coupled K: V''(1) = {Vpp_IR_c[idx_vac]:+.2f}, m^2 = {np.asarray(m2_phys_c).item():.2f}
     - Bieg K {'wzmacnia' if m2_phys_c > m2_f else 'oslania'} stabilnosc

  3. UV FIXED POINT:
     - {'ZNALEZIONY' if found_fp else 'NIE ZNALEZIONY w tym schemacie'}
     - K_IR/K_UV = {K_ratio:.4f}
     - {'K ograniczone -> pozytywny sygnal AS' if 0.1 < K_ratio < 10 else 'K zmienia sie znaczaco'}

  4. NASTEPNE KROKI:
     - Pelny FP search: skanowanie warunkow UV
     - Wykladniki krytyczne z linearyzacji
     - Porownanie z AS grawitacji (Reuter)
""")
else:
    print("  Integracja nieudana -- wymagana dalsza optymalizacja.")

print("GOTOWE.")
