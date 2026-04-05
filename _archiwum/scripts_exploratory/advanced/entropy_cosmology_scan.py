#!/usr/bin/env python3
"""
entropy_cosmology_scan.py
==========================
Problem O22 TGP: Kosmologiczna ewolucja Phi z sprzezeniem entropijnym.

Implementuje rownainia kosmologiczne TGP z entropia substratu S_Gamma:
    H^2 = (8*pi*G/3) * [rho_m + rho_Phi]
    drho_Phi/dt + 3*H*(rho_Phi + p_Phi) = -T * dS_Gamma/dt

Skanuje s0 in [0.001, 0.01, 0.05, 0.1] i oblicza w_de(z).
Porownuje z ograniczeniem DESI: w_de(z=0) = -0.45 +/- 0.27.

Uzycie:
    python entropy_cosmology_scan.py
"""

import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# ──────────────────────────────────────────────────────────────────────────────
# Parametry kosmologiczne (jednostki naturalne: H0=1, c=hbar=1, 8*pi*G/3=1)
# ──────────────────────────────────────────────────────────────────────────────
OMEGA_M   = 0.315    # gestosc materii
OMEGA_R   = 9.0e-5   # gestosc promieniowania
OMEGA_L   = 0.685    # gestosc ciemnej energii (LCDM)

# Parametry TGP
GAMMA_TGP = 0.5      # parametr potencjalu Phi
BETA_TGP  = 0.5      # beta = gamma (warunek prozniowy N0-5)
PHI0      = 1.0      # wartoscprozni Phi0

# Parametry entropii substratu
T0_GAMMA  = 0.10     # temperatura substratu T0 (w jednostkach H0)
S0_VALUES = [0.001, 0.01, 0.05, 0.1]  # wartosci do skanowania

# Zakres redshift
Z_MIN  = 0.0
Z_MAX  = 2.0
N_Z    = 200

# Zakres calkowan (N = ln(a), a=1 dzis)
N_INI  = -np.log(1.0 + Z_MAX) - 1.0   # troche przed z_max
N_FIN  = 0.5                            # troche po dzis (a ~ 1.6)

# Granica DESI (CPL parametryzacja)
DESI_W0_CENTRAL = -0.45
DESI_W0_SIGMA   = 0.27

# ──────────────────────────────────────────────────────────────────────────────
# Funkcje pomocnicze
# ──────────────────────────────────────────────────────────────────────────────

def V_phi(phi, gamma=GAMMA_TGP, phi0=PHI0):
    """
    Potencjal meksykanskiego kapelusza:
        V(Phi) = (gamma/4) * (Phi^2 - Phi0^2)^2 / Phi0^2
    W jednostkach znormalizowanych phi = Phi/Phi0:
        V(phi) = (gamma/4) * (phi^2 - 1)^2
    """
    return (gamma / 4.0) * (phi**2 - 1.0)**2


def dV_dphi(phi, gamma=GAMMA_TGP):
    """
    Pochodna potencjalu:
        dV/dphi = gamma * phi * (phi^2 - 1)
    """
    return gamma * phi * (phi**2 - 1.0)


def S_entropy(phi, s0):
    """
    Entropia substratu S_Gamma(phi):
        S_Gamma(phi) = k_B * N_B * s0 * (phi - ln(phi) - 1)
    Tu: bezwymiarowa, k_B = N_B = 1.
    """
    phi_safe = max(float(phi), 1e-12)
    return s0 * (phi_safe - np.log(phi_safe) - 1.0)


def dS_dphi(phi, s0):
    """
    Pochodna S_Gamma wzgledem phi:
        dS/dphi = s0 * (1 - 1/phi)
    """
    phi_safe = max(float(phi), 1e-12)
    return s0 * (1.0 - 1.0 / phi_safe)


def T_substrate(H, phi, T0=T0_GAMMA):
    """
    Temperatura substratu:
        T_Gamma(phi, H) = T0 * H / H0 * sqrt(Phi0/Phi) = T0 * H * phi^{-1/2}
    W jednostkach H0=1.
    """
    phi_safe = max(float(phi), 1e-12)
    return T0 * abs(H) / np.sqrt(phi_safe)


def H_squared_flat(a, phi, dphi_dN, s0, T0=T0_GAMMA,
                   gamma=GAMMA_TGP, beta=BETA_TGP):
    """
    Zmodyfikowany H^2 w TGP z entropia.

    Uklad rownan:
        H^2 = rho_m/a^3 + rho_r/a^4 + rho_Phi
    gdzie:
        rho_Phi = (1/2)(dphi/dt)^2 + V(phi) - T_Gamma * S_Gamma(phi)
        dphi/dt = H * dphi/dN

    Czlon entropijny: F_Gamma = -T_Gamma * S_Gamma(phi)
    modyfikuje potencjal efektywny.

    Samospojne rozwiazanie przez podstawienie dphi/dt = H * dphi/dN:
        H^2 = [rho_m/a^3 + rho_r/a^4 + V_eff + T0*S/(2*sqrt(phi))] /
              [1 - (1/2)*(dphi/dN)^2 + ...]
    """
    rho_m = OMEGA_M / a**3
    rho_r = OMEGA_R / a**4

    # Potencjal efektywny z czlonem entropijnym (zerowy order w T0)
    V    = V_phi(phi, gamma=gamma)
    S    = S_entropy(phi, s0)

    # Energia kinetyczna: (1/2)(phi_dot)^2 = (1/2) H^2 (dphi/dN)^2
    # Samospojne: H^2 (1 - (1/2)(dphi/dN)^2) = rho_m + rho_r + V_eff
    kin_coeff = 0.5 * dphi_dN**2

    # Czlon entropijny w potencjale efektywnym (wspolczynnik T_Gamma ~ T0*H)
    # V_eff_total = V - T_Gamma * S = V - T0*H*S/sqrt(phi)
    # Wchodzimy H samospojnie: iteracja
    # Przyblizenie 0-rzedowe: H2_0 = (rho_m + rho_r + V) / (1 - kin_coeff)
    denom = 1.0 - kin_coeff
    if denom <= 1e-10:
        denom = 1e-10

    H2_0 = (rho_m + rho_r + V) / denom
    if H2_0 < 0:
        H2_0 = 1e-10

    # Pierwsza iteracja z czlonem entropijnym
    H_0 = np.sqrt(H2_0)
    T_G  = T_substrate(H_0, phi, T0=T0)
    V_eff_total = V - T_G * S
    H2_1 = (rho_m + rho_r + V_eff_total) / denom
    if H2_1 < 0:
        H2_1 = H2_0 * 0.01  # ochrona

    return max(H2_1, 1e-20)


# ──────────────────────────────────────────────────────────────────────────────
# Uklad ODE w zmiennej N = ln(a)
# ──────────────────────────────────────────────────────────────────────────────

def ode_system(N, state, s0, T0=T0_GAMMA, gamma=GAMMA_TGP, beta=BETA_TGP):
    """
    Uklad ODE dla (phi, u=dphi/dN):
        d(phi)/dN = u
        d(u)/dN   = -3*u - dV_eff/dphi / H^2 - (entropia)
    gdzie:
        dV_eff/dphi = dV/dphi - T_Gamma * dS/dphi + czlon_korekcyjny
    """
    phi, u = state
    phi_safe = max(phi, 1e-12)
    a = np.exp(N)

    H2 = H_squared_flat(a, phi_safe, u, s0, T0=T0, gamma=gamma, beta=beta)
    H  = np.sqrt(max(H2, 1e-20))

    # Sila gradientowa potencjalu
    dV = dV_dphi(phi_safe, gamma=gamma)

    # Czlon entropijny: -T_Gamma * dS/dphi
    T_G   = T_substrate(H, phi_safe, T0=T0)
    dS    = dS_dphi(phi_safe, s0)
    dF_dp = -T_G * dS  # = T_Gamma * s0 * (1/phi - 1)

    # Czlon tlumienia: -3 (w N-zmiennych H^2 wchodzi naturalnie)
    # Dodatkowe tlumienie entropijne: -delta * s0 * u * H (z zadania O22)
    # delta = s0 * T0 (bezwymiarowe sprzezenie)
    delta_coupling = s0 * T0
    entropy_damping = -delta_coupling * u  # dodatkowe tlumienie

    # Rownanie ruchu:
    #   d^2phi/dN^2 + 3 d phi/dN + (dV + dF_dp)/H^2 = 0
    du_dN = -3.0 * u - (dV + dF_dp) / H2 + entropy_damping

    return [u, du_dN]


# ──────────────────────────────────────────────────────────────────────────────
# Rownananie stanu ciemnej energii w_de(z)
# ──────────────────────────────────────────────────────────────────────────────

def compute_w_de(N_arr, phi_arr, dphi_dN_arr, s0, T0=T0_GAMMA,
                 gamma=GAMMA_TGP):
    """
    Oblicza w_de(z) = p_Phi / rho_Phi dla tablicy N.

    rho_Phi = (1/2)(H*dphi/dN)^2 + V(phi) - T_G*S
    p_Phi   = (1/2)(H*dphi/dN)^2 - V(phi) + T_G*S (odwrocony znak V i S)
    """
    a_arr    = np.exp(N_arr)
    z_arr    = 1.0 / a_arr - 1.0

    w_de_arr = np.zeros(len(N_arr))

    for i, (N, phi, u) in enumerate(zip(N_arr, phi_arr, dphi_dN_arr)):
        a = a_arr[i]
        phi_safe = max(phi, 1e-12)

        H2 = H_squared_flat(a, phi_safe, u, s0, T0=T0, gamma=gamma)
        H  = np.sqrt(max(H2, 1e-20))
        T_G = T_substrate(H, phi_safe, T0=T0)
        S   = S_entropy(phi_safe, s0)

        # Czlon kinetyczny
        kin = 0.5 * H2 * u**2

        # Potencjal i czlon entropijny
        V   = V_phi(phi_safe, gamma=gamma)
        F   = -T_G * S  # swobodna energia entropiczna

        rho_Phi = kin + V + F
        p_Phi   = kin - V - F

        if abs(rho_Phi) < 1e-20:
            w_de_arr[i] = -1.0
        else:
            w_de_arr[i] = p_Phi / rho_Phi

    return z_arr, w_de_arr


# ──────────────────────────────────────────────────────────────────────────────
# Skan parametrow s0
# ──────────────────────────────────────────────────────────────────────────────

def run_scan():
    """Skanuje s0 i oblicza w_de(z=0)."""

    print("=" * 65)
    print("TGP - PROBLEM O22: KOSMOLOGICZNA EWOLUCJA PHI Z ENTROPIA")
    print(f"Parametry: OMEGA_M={OMEGA_M}, GAMMA_TGP={GAMMA_TGP}, T0={T0_GAMMA}")
    print(f"Zakres: z in [{Z_MIN}, {Z_MAX}]")
    print("=" * 65)

    # Warunki poczatkowe (dzis: a=1, phi ~ 1, phi niemal zamrozone)
    phi_ini  = 1.0 + 1e-4   # phi niemal przy prozni
    u_ini    = 1e-5          # prawie nieruchome pole

    # Zakres N = ln(a)
    N_start = np.log(1.0 / (1.0 + Z_MAX)) - 0.5  # nieco przed z_max
    N_end   = 0.1                                   # troche po dzis

    N_eval  = np.linspace(N_start, N_end, N_Z)

    results = {}

    print("\nWyniki dla roznych s0:")
    print("-" * 65)
    print(f"{'s0':>8}  {'w_de(z=0)':>12}  {'Delta_w=w+1':>14}  {'DESI OK?':>10}")
    print("-" * 65)

    desi_ok_list = []

    for s0 in S0_VALUES:
        try:
            sol = solve_ivp(
                fun=lambda N, y: ode_system(N, y, s0=s0, T0=T0_GAMMA),
                t_span=(N_start, N_end),
                y0=[phi_ini, u_ini],
                t_eval=N_eval,
                method='RK45',
                rtol=1e-8,
                atol=1e-10,
                dense_output=False,
            )

            if not sol.success:
                print(f"  s0={s0:.3f}: ODE nie zbiegla ({sol.message})")
                continue

            phi_sol    = sol.y[0]
            dphi_dN    = sol.y[1]
            N_sol      = sol.t

            z_arr, w_de_arr = compute_w_de(
                N_sol, phi_sol, dphi_dN, s0=s0, T0=T0_GAMMA)

            # Interpolacja do z=0
            # Znajdz indeks najblizszy z=0
            idx_z0 = np.argmin(np.abs(z_arr))
            w_at_z0 = w_de_arr[idx_z0]

            # Jesli brak punktu blisko z=0, interpoluj
            if np.abs(z_arr[idx_z0]) > 0.05:
                try:
                    interp_fn = interp1d(z_arr, w_de_arr, kind='linear',
                                         fill_value='extrapolate')
                    w_at_z0 = float(interp_fn(0.0))
                except Exception:
                    pass

            delta_w = w_at_z0 + 1.0

            # Zgodnosc z DESI: w_de(z=0) in [-0.72, -0.18] (DESI_W0_CENTRAL +/- DESI_W0_SIGMA)
            desi_lo = DESI_W0_CENTRAL - DESI_W0_SIGMA
            desi_hi = DESI_W0_CENTRAL + DESI_W0_SIGMA
            desi_ok = (desi_lo <= w_at_z0 <= desi_hi)
            desi_ok_list.append((s0, w_at_z0, desi_ok))

            print(f"  s0={s0:<6.3f}  w_de(z=0)={w_at_z0:>10.5f}  "
                  f"Delta_w={delta_w:>12.5f}  {'TAK' if desi_ok else 'NIE':>10}")

            results[s0] = {
                'z_arr':    z_arr,
                'w_de_arr': w_de_arr,
                'phi_arr':  phi_sol,
                'N_arr':    N_sol,
                'w_at_z0':  w_at_z0,
                'delta_w':  delta_w,
                'desi_ok':  desi_ok,
            }

        except Exception as e:
            print(f"  s0={s0:.3f}: BLAD - {e}")

    print("-" * 65)
    print(f"\nGranica DESI (CPL): w_de(z=0) = {DESI_W0_CENTRAL} +/- {DESI_W0_SIGMA}")
    print(f"Okno zgodnosci DESI: [{DESI_W0_CENTRAL - DESI_W0_SIGMA:.3f}, "
          f"{DESI_W0_CENTRAL + DESI_W0_SIGMA:.3f}]")

    print("\nWartosci s0 zgodne z DESI:")
    any_ok = False
    for s0, w0, ok in desi_ok_list:
        if ok:
            print(f"  s0={s0}: w_de(z=0)={w0:.5f} [ZGODNE]")
            any_ok = True
    if not any_ok:
        print("  Brak wartosci s0 w oknie DESI (model wymaga kalibracji T0/s0).")

    return results


# ──────────────────────────────────────────────────────────────────────────────
# Wyniki szczegolowe
# ──────────────────────────────────────────────────────────────────────────────

def print_detailed_results(results):
    """Drukuje szczegolowe wyniki dla kazdego s0."""
    print("\n" + "=" * 65)
    print("SZCZEGOLOWE WYNIKI w_de(z) DLA WYBRANYCH z")
    print("=" * 65)
    z_checkpoints = [0.0, 0.5, 1.0, 2.0]

    for s0, res in results.items():
        z_arr    = res['z_arr']
        w_de_arr = res['w_de_arr']

        print(f"\n  s0 = {s0}:")
        try:
            fn = interp1d(z_arr, w_de_arr, kind='linear',
                          fill_value='extrapolate', bounds_error=False)
            for zc in z_checkpoints:
                w = float(fn(zc))
                print(f"    w_de(z={zc:.1f}) = {w:.6f}   Delta_w = {w+1:.6f}")
        except Exception as e:
            print(f"    Blad interpolacji: {e}")


def print_latex_summary(results):
    """Drukuje wyniki w formacie gotowym do cytowania w LaTeX."""
    print("\n" + "=" * 65)
    print("WYNIKI DO LATEX (O22)")
    print("=" * 65)
    for s0 in S0_VALUES:
        if s0 in results:
            w0    = results[s0]['w_at_z0']
            dw    = results[s0]['delta_w']
            ok    = results[s0]['desi_ok']
            status = "ZGODNE z DESI" if ok else "POZA OKNEM DESI"
            print(f"  s0={s0:.3f}: w_de(z=0) = {w0:.5f}  "
                  f"Delta_w = w+1 = {dw:.5f}  [{status}]")


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    results = run_scan()

    if results:
        print_detailed_results(results)
        print_latex_summary(results)

    # Wykresy opcjonalne
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        fig.suptitle("TGP O22: Kosmologiczna ewolucja Phi z entropia", fontsize=13)

        colors = ['b', 'g', 'r', 'm']
        z_plot = np.linspace(0.0, Z_MAX, 200)

        ax = axes[0]
        for i, (s0, res) in enumerate(results.items()):
            z_arr    = res['z_arr']
            w_de_arr = res['w_de_arr']
            try:
                fn = interp1d(z_arr, w_de_arr, kind='linear',
                              fill_value='extrapolate', bounds_error=False)
                w_plot = fn(z_plot)
                ax.plot(z_plot, w_plot, color=colors[i % 4],
                        label=f's0={s0}, w0={res["w_at_z0"]:.4f}')
            except Exception:
                pass

        # DESI band
        ax.axhspan(DESI_W0_CENTRAL - DESI_W0_SIGMA,
                   DESI_W0_CENTRAL + DESI_W0_SIGMA,
                   alpha=0.15, color='orange', label='DESI 1sigma')
        ax.axhline(-1.0, color='k', ls='--', alpha=0.5, label='w=-1 (LCDM)')
        ax.set_xlabel('z')
        ax.set_ylabel('w_de(z)')
        ax.set_title('Rownanie stanu ciemnej energii w_de(z)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, Z_MAX)
        ax.set_ylim(-2.0, 0.5)

        ax = axes[1]
        s0_list = list(results.keys())
        w0_list = [results[s]['w_at_z0'] for s in s0_list]
        dw_list = [results[s]['delta_w'] for s in s0_list]
        x_pos   = range(len(s0_list))
        bars = ax.bar(x_pos, dw_list,
                      color=[colors[i % 4] for i in x_pos], alpha=0.7)
        ax.axhspan(-DESI_W0_SIGMA - 1 + DESI_W0_CENTRAL,
                   DESI_W0_SIGMA - 1 + DESI_W0_CENTRAL,
                   alpha=0.15, color='orange', label='DESI okno Delta_w')
        ax.axhline(0.0, color='k', ls='--', alpha=0.5)
        ax.set_xticks(list(x_pos))
        ax.set_xticklabels([f's0={s}' for s in s0_list])
        ax.set_ylabel('Delta_w = w_de(z=0) + 1')
        ax.set_title('Odchylenie Delta_w od LCDM')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        out_path = ("C:/Users/Mateusz/Documents/ObsydnianMain/TGP/TGP_v1/"
                    "plots/O22_entropy_cosmology_scan.png")
        plt.savefig(out_path, dpi=100)
        print(f"\nWykres zapisany: {out_path}")
        plt.close()
    except Exception as e:
        print(f"\n[INFO] Wykresy niedostepne: {e}")

    return results


if __name__ == '__main__':
    main()
