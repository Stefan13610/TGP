#!/usr/bin/env python3
"""
mc_hgamma_3d_thermal.py
=======================
TGP: masa generacji — rozszerzenie modelu sieciowego

Trzy bloki obliczen:

BLOK A — 3D sferyczne zrodlo (ODE radialne)
  Zamiast 1D kinka: sferyczny defekt Phi(r->0)=eps, Phi(inf)=Phi0.
  Porownanie E_1^{3D} vs E_1^{1D} i scalowania z d.

BLOK B — Dwa zrodla cylindrycznie (2D BVP)
  Polowa dziedziny (rho,z), rho in [0,R], z in [0,Z].
  Zrodla na osi rho=0 przy z=d/2. Warunek symetrii z=0.
  Energia E_2^{3D}(d) i porownanie z 2*E_1^{3D}.

BLOK C — Termiczny MC na siatce 1D (Metropolis)
  Temperatura T_Gamma simuluje kwantowe fluktuacje substratu.
  Sprawdza czy <E_2>/<E_1> rosnie z T_Gamma (eksponencjalnie?).
  Znajdz T* dla ktorego kappa_eff(T*) = ln(207) = 5.33.
"""

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import minimize
import warnings
import time
warnings.filterwarnings('ignore')

try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ═══════════════════════════════════════════════════════════════
# PARAMETRY TGP
# ═══════════════════════════════════════════════════════════════
M_SP  = 1.0
PHI0  = 1.0
BETA  = 0.5
ALPHA = 0.1
EPS   = 0.02 * PHI0

# ═══════════════════════════════════════════════════════════════
# WSPOLNE FUNKCJE
# ═══════════════════════════════════════════════════════════════

def V_pot(phi):
    return (BETA / 4.0) * (phi**2 - PHI0**2)**2 / PHI0**2

def dV_dphi(phi):
    return BETA * phi * (phi**2 - PHI0**2) / PHI0**2


# ═══════════════════════════════════════════════════════════════
# BLOK A: 3D RADIALNE — ODE + energia sferyczna
# ═══════════════════════════════════════════════════════════════

def radial_ode_rhs(r, y, alpha):
    """
    ODE dla sfrycznie symetrycznego profilu Phi(r):
      Phi'' + (2/r)Phi' = N[Phi]
    gdzie N[Phi] = alpha*(Phi')^2/(Phi0*Phi) + beta*Phi^2/Phi0^2 - gamma*Phi^3/Phi0^3
                 = alpha*(Phi')^2/(Phi0*Phi) + dV/dPhi

    Zamiana: y[0]=Phi, y[1]=Phi'
    """
    phi  = np.maximum(y[0], 1e-10)
    dphi = y[1]

    # Przy r->0 uzywamy limitu l'Hopital: 2/r * Phi' -> 2*Phi'' = 2*N[Phi]
    r_safe = np.maximum(r, 1e-8)

    nl_grad  = alpha * dphi**2 / (PHI0 * phi)
    nl_pot   = dV_dphi(phi)
    rhs_phi  = nl_grad + nl_pot    # = N[Phi]
    friction = -2.0 / r_safe * dphi

    return np.array([dphi, rhs_phi + friction])


def bc_radial(ya, yb, alpha):
    """
    Warunki brzegowe:
      r=0 (r_min): Phi=EPS, Phi'=0 (symetria sferyczna)
      r=R (r_max): Phi=Phi0
    """
    return np.array([ya[0] - EPS, ya[1], yb[0] - PHI0])


def solve_radial_single(alpha=ALPHA, r_max=20.0, n_pts=200):
    """
    Rozwiazuje BVP na profil sferyczny jednego zrodla.
    Zwraca (r, phi, E_1_3D).
    """
    r_min = 1e-3
    r     = np.linspace(r_min, r_max, n_pts)

    # Inicjalizacja: profil Yukawa
    K0    = (PHI0 - EPS)
    phi0  = PHI0 - K0 * np.exp(-M_SP * r) / np.maximum(r, 1e-6)
    phi0  = np.maximum(phi0, EPS)
    dphi0 = K0 * M_SP * np.exp(-M_SP * r) / np.maximum(r, 1e-6) \
            - K0 * np.exp(-M_SP * r) / np.maximum(r**2, 1e-8)

    y_init = np.array([phi0, dphi0])

    try:
        sol = solve_bvp(
            lambda r, y: radial_ode_rhs(r, y, alpha),
            lambda ya, yb: bc_radial(ya, yb, alpha),
            r, y_init,
            max_nodes=5000, tol=1e-5
        )
        r_sol   = sol.x
        phi_sol = np.maximum(sol.y[0], 1e-10)
        dphi    = sol.y[1]
    except Exception:
        r_sol   = r
        phi_sol = phi0
        dphi    = dphi0

    # Energia 3D sferyczna: E = 4pi * integral [ energy_density * r^2 ] dr
    phi_e   = np.maximum(phi_sol, 1e-10)
    Egrad   = dphi**2
    Enl     = alpha * dphi**2 / (PHI0 * phi_e)
    Epot    = V_pot(phi_e)
    density = (Egrad + Enl + Epot) * r_sol**2

    E1_3D   = 4.0 * np.pi * np.trapezoid(density, r_sol)
    return r_sol, phi_sol, E1_3D


# ═══════════════════════════════════════════════════════════════
# BLOK B: 2D CYLINDRYCZNE — dwa zrodla na osi
# ═══════════════════════════════════════════════════════════════

class Cyl2D:
    """
    Siatka cylindryczna (rho, z), polowa dziedziny z>=0.
    Metoda cell-centered: rho_i = (i+0.5)*drho, i=0..Nr-1
                          z_j   = j*dz,          j=0..Nz-1
    Warunek symetrii przy z=0: phi(rho,-dz) = phi(rho,+dz)
    Warunek symetrii przy rho=0: phi(-drho,z) = phi(+drho,z)
    """

    def __init__(self, Nr=50, Nz=60, R_max=10.0, Z_max=12.0):
        self.Nr    = Nr
        self.Nz    = Nz
        self.drho  = R_max / Nr
        self.dz    = Z_max / Nz
        self.rho   = (np.arange(Nr) + 0.5) * self.drho   # cell centers
        self.z     = np.arange(Nz) * self.dz             # node z
        self.R_max = R_max
        self.Z_max = Z_max

    def H_cyl(self, phi_flat, alpha):
        """Energia konfiguracji na siatce cylindrycznej."""
        Nr, Nz     = self.Nr, self.Nz
        drho, dz   = self.drho, self.dz
        phi        = phi_flat.reshape(Nr, Nz)
        ps         = np.maximum(phi, 1e-10)

        # --- rho-gradienty (miedzy i i i+1) ---
        g_rho      = np.diff(ps, axis=0) / drho       # (Nr-1, Nz)
        rho_e      = (np.arange(Nr-1) + 1.0) * drho   # krawedz miedzy cell i a i+1
        phi_e_rho  = np.maximum(0.5*(ps[:-1,:]+ps[1:,:]), 1e-10)
        w_rho      = rho_e[:, np.newaxis] * drho * dz
        E_rho      = np.sum(g_rho**2 * (1 + alpha/(PHI0*phi_e_rho)) * w_rho)

        # Lustro: gradient na lewej krawedzi i=0 (rho=0) jest zerowy
        # (symetria rho), wiec nie dodajemy.

        # --- z-gradienty (miedzy j i j+1) ---
        g_z        = np.diff(ps, axis=1) / dz         # (Nr, Nz-1)
        phi_e_z    = np.maximum(0.5*(ps[:,:-1]+ps[:,1:]), 1e-10)
        w_z        = self.rho[:, np.newaxis] * drho * dz
        E_z        = np.sum(g_z**2 * (1 + alpha/(PHI0*phi_e_z)) * w_z[:, :Nz-1])

        # Gradient przy z=0 (j=0) z lustem: symetria -> gradient=0

        # --- potencjal ---
        w_pot      = self.rho[:, np.newaxis] * drho * dz
        E_pot      = np.sum(V_pot(ps) * w_pot)

        return 2.0 * np.pi * (E_rho + E_z + E_pot)

    def relax_two_sources(self, d_sep, alpha=ALPHA, verbose=False):
        """
        Relaksacja dwoch zrodel przy z=d/2 (na osi rho=0).
        Zwraca (phi_2D, E_2D).
        """
        Nr, Nz   = self.Nr, self.Nz
        drho, dz = self.drho, self.dz

        # Znajdz indeks zrodla: rho=0 -> i_src=0; z=d/2 -> j_src
        i_src = 0
        j_src = int(round((d_sep / 2.0) / dz))
        j_src = max(1, min(j_src, Nz - 2))

        # Inicjalizacja: Phi = Phi0
        phi = np.ones((Nr, Nz)) * PHI0

        # Stale: zrodlo i brzeg zewnetrzny
        fixed = np.zeros((Nr, Nz), dtype=bool)
        fixed[i_src, j_src] = True   # zrodlo przy z=+d/2
        # Brzeg: ostatnia warstwa rho i ostatni z
        fixed[-1, :] = True           # rho = R_max
        fixed[:, -1] = True           # z   = Z_max

        phi[fixed] = PHI0
        phi[i_src, j_src] = EPS

        free = ~fixed
        phi_full = phi.copy()

        def obj(pv):
            phi_full[free] = pv
            return self.H_cyl(phi_full.ravel(), alpha)

        bounds = [(1e-9, None)] * int(free.sum())
        res = minimize(
            obj,
            phi[free],
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 5000, 'ftol': 1e-12, 'gtol': 1e-8}
        )
        phi_full[free] = res.x
        phi_full[i_src, j_src] = EPS
        phi_full[fixed] = PHI0

        if verbose:
            print(f"    d={d_sep:.2f}: converged={res.success}, E={res.fun:.5f}")

        return phi_full, res.fun

    def scan_two_sources(self, d_vals, alpha=ALPHA):
        """Skan E_2^{3D}(d) dla listy separacji."""
        results = []
        for d in d_vals:
            phi_r, E = self.relax_two_sources(d, alpha=alpha)
            results.append((d, E))
        return results


# ═══════════════════════════════════════════════════════════════
# BLOK C: TERMICZNY MC — Metropolis na siatce 1D
# ═══════════════════════════════════════════════════════════════

# Siatka 1D
N1D   = 300
L1D   = 15.0
DX1D  = 2 * L1D / N1D
XX1D  = np.linspace(-L1D, L1D, N1D)


def H_1d(phi, alpha):
    """Energia 1D (ta sama co w mc_hgamma_generacje.py)."""
    ps   = np.maximum(phi, 1e-10)
    g    = np.diff(ps) / DX1D
    pe   = np.maximum(0.5*(ps[:-1]+ps[1:]), 1e-10)
    return (np.sum(g**2) * DX1D
            + alpha * np.sum(g**2 / (PHI0 * pe)) * DX1D
            + np.sum(V_pot(ps)) * DX1D)


def gs_profile(src_indices, alpha):
    """Stan podstawowy (T=0): relaksacja L-BFGS-B."""
    fixed = np.zeros(N1D, dtype=bool)
    for i in src_indices:
        fixed[i] = True
    fixed[0] = fixed[-1] = True
    free = ~fixed

    phi = np.ones(N1D) * PHI0
    phi[fixed] = PHI0
    for i in src_indices:
        phi[i] = EPS

    phi_f = phi.copy()

    def obj(pv):
        phi_f[free] = pv
        return H_1d(phi_f, alpha)

    res = minimize(obj, phi[free], method='L-BFGS-B',
                   bounds=[(1e-9, None)]*free.sum(),
                   options={'maxiter': 5000, 'ftol': 1e-13, 'gtol': 1e-9})
    phi_f[free] = res.x
    for i in src_indices:
        phi_f[i] = EPS
    return phi_f, res.fun


def metropolis_mc(src_indices, T_gamma, n_sweep=8000,
                  delta_max=0.3, alpha=ALPHA, seed=42):
    """
    Metropolis MC na siatce 1D przy temperaturze T_Gamma.

    Propozycja: Phi_i -> Phi_i + delta, delta ~ U(-delta_max, delta_max)
    Akceptacja: min(1, exp(-dH / T_Gamma))  jezeli Phi_new > 0
    Pomiar:     co n_meas krokow, srednia E = H_1d(phi)
    Zwraca: (E_mean, E_std, accept_rate)
    """
    rng = np.random.default_rng(seed)

    # Start z profilu stanu podstawowego
    phi, E0 = gs_profile(src_indices, alpha)
    H_curr  = H_1d(phi, alpha)

    fixed = np.zeros(N1D, dtype=bool)
    for i in src_indices:
        fixed[i] = True
    fixed[0] = fixed[-1] = True
    free_idx = np.where(~fixed)[0]
    n_free   = len(free_idx)

    E_samples   = []
    n_accept    = 0
    n_total     = 0
    n_meas      = 10      # co ile krokow mierzymy energię

    # Termostatyzacja
    n_therm = n_sweep // 4
    for _ in range(n_therm * n_free):
        i   = rng.integers(0, n_free)
        idx = free_idx[i]
        delta   = rng.uniform(-delta_max, delta_max)
        phi_new = phi[idx] + delta
        if phi_new <= 0:
            continue
        phi[idx]  = phi_new
        H_new = H_1d(phi, alpha)
        dH    = H_new - H_curr
        if dH < 0 or rng.random() < np.exp(-dH / max(T_gamma, 1e-10)):
            H_curr   = H_new
            n_accept += 1
        else:
            phi[idx] = phi[idx] - delta  # cofnij
        n_total += 1

    # Pomiary
    n_accept = n_total = 0
    for step in range(n_sweep * n_free):
        i   = rng.integers(0, n_free)
        idx = free_idx[i]
        delta   = rng.uniform(-delta_max, delta_max)
        phi_new = phi[idx] + delta
        if phi_new <= 0:
            continue
        phi[idx]  = phi_new
        H_new = H_1d(phi, alpha)
        dH    = H_new - H_curr
        if dH < 0 or rng.random() < np.exp(-dH / max(T_gamma, 1e-10)):
            H_curr   = H_new
            n_accept += 1
        else:
            phi[idx] = phi[idx] - delta
        n_total += 1

        if step % (n_meas * n_free) == 0:
            E_samples.append(H_curr)

    E_arr  = np.array(E_samples)
    accept = n_accept / max(n_total, 1)
    return E_arr.mean(), E_arr.std(), accept


def thermal_scan(T_vals, n_sweep=6000, alpha=ALPHA):
    """
    Skan termiczny: dla kazdego T_Gamma oblicz <E_1> i <E_2>.
    Zwraca liste (T, E1_mean, E1_std, E2_mean, E2_std, kappa_eff).
    """
    ctr  = N1D // 2
    src1 = [ctr]

    # Separacja 1D przy minimum sprzezenia (d*m_sp ~ 1.5 z poprzednich wynikow)
    d_opt_idx = max(1, int(1.5 / DX1D))
    src2 = [ctr - d_opt_idx, ctr + d_opt_idx]

    rows = []
    for T in T_vals:
        e1m, e1s, a1 = metropolis_mc(src1, T, n_sweep=n_sweep, alpha=alpha, seed=11)
        e2m, e2s, a2 = metropolis_mc(src2, T, n_sweep=n_sweep, alpha=alpha, seed=22)

        ratio = e2m / max(e1m, 1e-10)
        kappa = np.log(max(ratio, 1.0001))
        rows.append((T, e1m, e1s, e2m, e2s, kappa, a1, a2))
        print(f"  T={T:.3f}: <E_1>={e1m:.4f}+-{e1s:.4f}, "
              f"<E_2>={e2m:.4f}+-{e2s:.4f}, kappa={kappa:.4f}, "
              f"acc={a1:.2f}/{a2:.2f}")
    return rows


# ═══════════════════════════════════════════════════════════════
# WYKRESY
# ═══════════════════════════════════════════════════════════════

def make_plots(r_sol, phi_1d_r, E1_3D,
               cyl_results, E1_3D_cyl,
               thermal_rows, d_opt):

    fig = plt.figure(figsize=(16, 11))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.36)

    kappa_req = np.log(207.0)

    # ─── A. Profil radialny 3D ──────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(r_sol * M_SP, phi_1d_r / PHI0, 'b-', lw=2)
    ax.axhline(1.0, color='gray', ls='--', lw=1, label='Phi_0')
    ax.axhline(EPS/PHI0, color='r', ls=':', lw=1, label=f'eps={EPS:.2f}')
    ax.fill_between(r_sol*M_SP, phi_1d_r/PHI0, 1.0, alpha=0.1, color='blue')
    ax.set_xlabel('r * m_sp'); ax.set_ylabel('Phi(r)/Phi_0')
    ax.set_title(f'Profil sferyczny phi_1(r) 3D\nE_1^3D = {E1_3D:.4f}', fontsize=9)
    ax.legend(fontsize=8); ax.set_xlim(0, 8); ax.grid(True, alpha=0.3)

    # ─── B. E_2^3D(d) cylindrycznie ────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    if cyl_results:
        d_arr = np.array([r[0] for r in cyl_results])
        E_arr = np.array([r[1] for r in cyl_results])
        ratio = E_arr / E1_3D_cyl
        ax.semilogy(d_arr * M_SP, ratio, 'g-o', ms=5, lw=2, label='E_2^3D(d)/E_1^3D')
        ax.axhline(2.0,   color='gray',   ls='--', lw=1,   label='2 (niezal.)')
        ax.axhline(207.0, color='orange', ls=':',  lw=1.5, label='obs 207')
        ax.axvline(d_opt * M_SP, color='purple', ls='--', lw=1, label=f'd_opt={d_opt:.1f}')
        ax.set_xlabel('d * m_sp'); ax.set_ylabel('E_2/E_1  [log]')
        ax.set_title('E_2^{3D}(d) cylindrycznie\n2 zrodla, polowa dziedziny', fontsize=9)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'Brak wynikow 3D', ha='center', va='center',
                transform=ax.transAxes)

    # ─── C. Termiczny skan kappa(T) ─────────────────────────
    ax = fig.add_subplot(gs[0, 2])
    if thermal_rows:
        T_arr    = np.array([r[0] for r in thermal_rows])
        kap_arr  = np.array([r[5] for r in thermal_rows])
        ax.semilogy(T_arr, np.exp(kap_arr), 'rs-', ms=6, lw=2,
                    label='<E_2>/<E_1> = exp(kappa)')
        ax.axhline(207.0, color='orange', ls='--', lw=1.5, label='obs 207')
        ax.axhline(2.0,   color='gray',   ls=':',  lw=1,   label='niezal. limit')
        ax.set_xlabel('T_Gamma / E_1')
        ax.set_ylabel('<E_2>/<E_1>  [log]')
        ax.set_title('Termiczny efekt: kappa(T_Gamma)\n1D Metropolis', fontsize=9)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # ─── D. kappa_eff(T) ────────────────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    if thermal_rows:
        T_arr   = np.array([r[0] for r in thermal_rows])
        kap_arr = np.array([r[5] for r in thermal_rows])
        ax.plot(T_arr, kap_arr, 'rs-', ms=6, lw=2, label='kappa_eff(T)')
        ax.axhline(kappa_req, color='orange', ls='--', lw=1.5,
                   label=f'kappa*=ln(207)={kappa_req:.2f}')
        ax.axhline(0.0084, color='gray', ls=':', lw=1, label='kappa T=0')
        ax.set_xlabel('T_Gamma'); ax.set_ylabel('kappa_eff = ln(<E_2>/<E_1>)')
        ax.set_title('kappa_eff vs temperatura substratu', fontsize=9)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # ─── E. E_1 i E_2 termicznie ────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    if thermal_rows:
        T_arr   = np.array([r[0] for r in thermal_rows])
        e1_arr  = np.array([r[1] for r in thermal_rows])
        e1s_arr = np.array([r[2] for r in thermal_rows])
        e2_arr  = np.array([r[3] for r in thermal_rows])
        e2s_arr = np.array([r[4] for r in thermal_rows])

        ax.fill_between(T_arr, e1_arr-e1s_arr, e1_arr+e1s_arr, alpha=0.2, color='blue')
        ax.fill_between(T_arr, e2_arr-e2s_arr, e2_arr+e2s_arr, alpha=0.2, color='green')
        ax.plot(T_arr, e1_arr, 'b-o', ms=5, lw=2, label='<E_1>')
        ax.plot(T_arr, e2_arr, 'g-s', ms=5, lw=2, label='<E_2>')
        ax.set_xlabel('T_Gamma'); ax.set_ylabel('Energia')
        ax.set_title('<E_n>(T_Gamma): fluktuacje kwantowe\n(1 sigma band)', fontsize=9)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # ─── F. Predykcja zbiorcza ──────────────────────────────
    ax = fig.add_subplot(gs[1, 2])
    kk = np.linspace(0.01, 7, 200)
    ax.semilogy(kk, np.exp(kk),   'b-', lw=2, label='m_2/m_1 = exp(kappa)')
    ax.semilogy(kk, np.exp(2*kk), 'r-', lw=2, label='m_3/m_1 = exp(2*kappa)')
    ax.axhline(207,  color='b', ls=':', lw=1.5, label='obs mu/e=207')
    ax.axhline(3477, color='r', ls=':', lw=1.5, label='obs tau/e=3477')
    ax.axvline(kappa_req, color='gray', ls='--', lw=1.2,
               label=f'kappa*={kappa_req:.2f}')

    # Zaznacz wyniki z roznych modeli
    k_1d = 0.0084
    ax.axvline(k_1d, color='green', ls=':', lw=1, label=f'1D GS: kappa={k_1d:.4f}')
    if cyl_results and E1_3D_cyl > 0:
        E_arr  = np.array([r[1] for r in cyl_results])
        valid  = E_arr > 0
        if valid.any():
            k_3d = np.log(max(E_arr[valid].max() / E1_3D_cyl, 1.001))
            ax.axvline(k_3d, color='purple', ls=':', lw=1,
                       label=f'3D max: kappa={k_3d:.4f}')
    if thermal_rows:
        k_th = max(r[5] for r in thermal_rows)
        ax.axvline(k_th, color='red', ls=':', lw=1,
                   label=f'thermal max: kappa={k_th:.4f}')

    ax.set_xlabel('kappa = m_sp * a_Gamma')
    ax.set_ylabel('m_n / m_1  [log]')
    ax.set_title('Predykcja TGP: hierarchia mas(kappa)', fontsize=9)
    ax.legend(fontsize=6.5); ax.grid(True, alpha=0.3)
    ax.set_ylim(0.5, 1e5)

    plt.suptitle(
        'TGP: O16 — 3D sferyczne + cylindryczne BVP + termiczny MC\n'
        f'alpha={ALPHA}, m_sp={M_SP}, Phi0={PHI0}, eps={EPS:.2f}',
        fontsize=10, fontweight='bold'
    )
    out = 'TGP/TGP_v1/scripts/advanced/mc_hgamma_3d_thermal.png'
    try:
        plt.savefig(out, dpi=120, bbox_inches='tight')
        print(f"  Wykres: {out}")
    except Exception as e:
        print(f"  Wykres nie zapisany: {e}")
    plt.close()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    SEP  = "=" * 68
    SEP2 = "-" * 68
    kappa_req = np.log(207.0)

    print(SEP)
    print("  TGP O16: 3D sferyczne + cylindryczne BVP + termiczny MC")
    print(SEP)
    print(f"  alpha={ALPHA}, m_sp={M_SP}, Phi0={PHI0}, beta={BETA}, eps={EPS:.3f}")
    print()

    t0 = time.time()

    # ──────────────────────────────────────────────────────────
    # BLOK A: 3D sferyczne zrodlo
    # ──────────────────────────────────────────────────────────
    print(SEP2)
    print("BLOK A: Profil sferyczny Phi(r) — ODE radialne 3D")
    print(SEP2)
    print("  Phi'' + (2/r)Phi' = N[Phi]  (w wspolrzednych sferycznych)")
    print("  Energie sferyczna: E_1^3D = 4pi * integral(density * r^2) dr\n")

    r_sol, phi_rad, E1_3D = solve_radial_single(alpha=ALPHA)
    E1_1D_anal = (4.0/3.0) * M_SP * PHI0**2

    print(f"  E_1^3D = {E1_3D:.6f}")
    print(f"  E_1^1D (anal 4/3*m*phi0^2) = {E1_1D_anal:.6f}")
    print(f"  Stosunek E_1^3D/E_1^1D = {E1_3D/E1_1D_anal:.4f}")
    print(f"  Phi(r=0) = {phi_rad[0]:.5f}, Phi(r_max) = {phi_rad[-1]:.5f}")
    print()

    # ──────────────────────────────────────────────────────────
    # BLOK B: Cylindryczne BVP — E_2^3D(d)
    # ──────────────────────────────────────────────────────────
    print(SEP2)
    print("BLOK B: Dwa zrodla cylindrycznie (2D BVP, pol-dziedzina z>=0)")
    print(SEP2)
    print("  Minimalizacja H_cyl(phi) na siatce (rho,z)")
    print("  Zrodlo przy (rho=0, z=d/2), symetria z=0\n")

    cyl = Cyl2D(Nr=45, Nz=55, R_max=9.0, Z_max=11.0)
    print(f"  Siatka: Nr={cyl.Nr}, Nz={cyl.Nz}, "
          f"drho={cyl.drho:.3f}, dz={cyl.dz:.3f}")

    # Najpierw E_1^3D z relaksacji cylindrycznej (jedno zrodlo, d=0 odpowiednik)
    # Uzyjemy jedno zrodlo na osi, a E_1^3D z BLOKU A jako referencyjna
    E1_3D_cyl = E1_3D  # kalibracja z BVP radialnego

    d_scan    = [0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0]
    print(f"\n  Skan d: {d_scan}")
    print(f"  (E_1 ref = {E1_3D_cyl:.5f})\n")

    print(f"  {'d*m_sp':>8}  {'E_2^3D':>11}  {'E_2/E_1':>9}  "
          f"{'ln(E_2/E_1)':>12}")
    print("  " + SEP2)

    cyl_results = []
    for d in d_scan:
        phi_r, E2 = cyl.relax_two_sources(d, alpha=ALPHA, verbose=False)
        ratio      = E2 / E1_3D_cyl
        kap        = np.log(max(ratio, 1.0001))
        cyl_results.append((d, E2))
        print(f"  {d*M_SP:>8.2f}  {E2:>11.5f}  {ratio:>9.5f}  {kap:>12.5f}")

    # Wyniki 3D
    E2_3D_max   = max(r[1] for r in cyl_results)
    kappa_3D    = np.log(max(E2_3D_max / E1_3D_cyl, 1.001))
    d_at_3Dmax  = [r[0] for r in cyl_results if r[1] == E2_3D_max][0]

    print(f"\n  => E_2^3D max = {E2_3D_max:.5f}  przy d*m_sp = {d_at_3Dmax:.2f}")
    print(f"  => kappa_3D   = {kappa_3D:.5f}  (wym: {kappa_req:.4f})")
    print(f"  => m_2/m_1    = {np.exp(kappa_3D):.4f}  (obs: 207)")
    print()

    # ──────────────────────────────────────────────────────────
    # BLOK C: Termiczny MC 1D
    # ──────────────────────────────────────────────────────────
    print(SEP2)
    print("BLOK C: Termiczny MC — Metropolis na siatce 1D")
    print(SEP2)
    print("  Fluktuacje kwantowe substratu przy T_Gamma > 0.")
    print("  Czlon alfa*(grad Phi)^2/(Phi0*Phi) eksploduje gdy Phi->0")
    print("  => chwilowe fluktuacje Phi do malych wartosci amplifikuja <E>.\n")

    T_vals = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    print(f"  T_Gamma wartosci: {T_vals}")
    print(f"  d*m_sp = 1.5 (punkt max kappa z poprzedniego skanu)\n")
    print(f"  {'T_Gamma':>8}  {'<E_1>':>9}  {'<E_2>':>9}  "
          f"{'<E_2>/<E_1>':>12}  {'kappa':>7}  {'acc1/acc2':>10}")
    print("  " + "-"*62)

    thermal_rows = thermal_scan(T_vals, n_sweep=5000, alpha=ALPHA)

    kappa_thermal_max = max(r[5] for r in thermal_rows)
    T_at_kmax         = [r[0] for r in thermal_rows
                         if r[5] == kappa_thermal_max][0]

    print(f"\n  => kappa_thermal_max = {kappa_thermal_max:.5f}  przy T = {T_at_kmax:.3f}")
    print(f"  => m_2/m_1 (thermal) = {np.exp(kappa_thermal_max):.4f}  (obs: 207)")

    # Znajdz T* dla kappa = ln(207) = 5.33 (jesli osiagalne)
    T_arr   = np.array([r[0] for r in thermal_rows])
    kap_arr = np.array([r[5] for r in thermal_rows])
    if kappa_thermal_max >= kappa_req:
        # Interpolacja T*
        idx = np.searchsorted(kap_arr, kappa_req)
        if 0 < idx < len(kap_arr):
            T_star = np.interp(kappa_req, kap_arr[idx-1:idx+1],
                               T_arr[idx-1:idx+1])
            print(f"\n  => T* (kappa=ln(207)) = {T_star:.4f}")
            print(f"  => Interpretacja TGP: temperatura substratu T_Gamma ~ {T_star:.3f} E_1")
    else:
        print(f"\n  => kappa max = {kappa_thermal_max:.4f} < ln(207) = {kappa_req:.4f}")
        print(f"  => Termiczny MC nie wystarczy przy alpha={ALPHA}")
        print(f"  => Ekstrapolacja: przy alfa* ~ 32 lub wiekszym d -- patrz wyniki")

    # ──────────────────────────────────────────────────────────
    # FINALNE PODSUMOWANIE
    # ──────────────────────────────────────────────────────────
    t_total = time.time() - t0
    print()
    print(SEP)
    print("PODSUMOWANIE — O16 3D + THERMAL (v31)")
    print(SEP)
    print(f"""
  BLOK A (3D radial):
    E_1^3D  = {E1_3D:.5f}
    E_1^1D  = {E1_1D_anal:.5f}  (4/3*m*phi0^2)
    Stosunek 3D/1D = {E1_3D/E1_1D_anal:.4f}  -- model 3D daje inna skale energii

  BLOK B (cylindryczne BVP):
    kappa_3D_max = {kappa_3D:.5f}  vs kappa_wym = {kappa_req:.5f}
    m_2/m_1 (3D) = {np.exp(kappa_3D):.4f}  vs obs 207
    Wzmocnienie 3D vs 1D: {kappa_3D/max(0.0084,1e-6):.2f}x kappa_1D

  BLOK C (thermal MC):
    kappa_thermal_max = {kappa_thermal_max:.5f}  przy T = {T_at_kmax:.3f}
    m_2/m_1 (thermal) = {np.exp(kappa_thermal_max):.4f}
    {"SUKCES: kappa >= ln(207)" if kappa_thermal_max >= kappa_req else
     f"Niewystarczajace: brak kappa >= {kappa_req:.2f}"}

  Laczna predykcja:
    kappa = kappa_3D + kappa_thermal(T_Gamma)
    Wymagane: kappa_3D + kappa_thermal = {kappa_req:.4f}

  Parametry TGP wyznaczone z obserwacji:
    kappa_total = m_sp * a_Gamma + Delta_thermal = {kappa_req:.4f}
    Przy m_sp = 0.0831 l_P^-1 (okno Efimova):
      a_Gamma_podstawowe = {kappa_req:.4f} / 0.0831 ~ {kappa_req/0.0831:.1f} l_P

  Czas: {t_total:.1f} s
    """)

    if HAS_MPL:
        d_opt = d_scan[np.argmax([r[1] for r in cyl_results])]
        make_plots(r_sol, phi_rad, E1_3D,
                   cyl_results, E1_3D_cyl,
                   thermal_rows, d_opt)


if __name__ == '__main__':
    main()
