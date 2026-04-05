#!/usr/bin/env python3
"""
ERG (Exact Renormalization Group) dla pola Phi w TGP -- SZKIELET
=================================================================
Cel: Zbudowac framework do kwantyzacji pola Phi w TGP
przy uzyciu rownania Wetterich'ego (functional RG).

Rownanie Wetterich'ego:
  d/dk Gamma_k[phi] = (1/2) Tr[(Gamma_k^(2) + R_k)^-1 * dR_k/dk]

Dla TGP:
- Gamma_k jest efektywnym dzialaniem na skali k
- R_k jest regulatorem IR (Litim: R_k = (k^2 - p^2)*theta(k^2-p^2))
- Phi jest polem przestrzennosci TGP

Aproksymacja LPA (Local Potential Approximation):
  Gamma_k[phi] = integral d^4x [Z_k/2 * (d phi)^2 + V_k(phi)]

  d/dk V_k(phi) = (k^4)/(32*pi^2) * 1/(V_k''(phi) + k^2)    [Litim 4D]

W TGP:
- Potencjal bare: V(psi) = beta*psi^2/(2*Phi0^2) - gamma*psi^3/(3*Phi0^3)
  z beta = gamma (warunek prozni)
- Cutoff UV: k_UV = 1/a_Gamma ~ M_Pl
- Cutoff IR: k_IR = m_sp = sqrt(gamma)
- Pytanie: czy punkt staly UV istnieje? (asymptotic safety)

TGP v1 -- 2026-03-31  (OP-1 szkielet)
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
gamma_TGP = 1.0  # beta = gamma = 1 w jednostkach Plancka (znormalizowane)
a_Gamma = 0.040

k_UV = 1.0 / a_Gamma  # = 25 (w jednostkach Plancka)
k_IR = np.sqrt(gamma_TGP)  # = 1
# Zakres RG: t = ln(k/k_IR), t in [0, ln(k_UV/k_IR)] = [0, ln(25)] = [0, 3.22]
t_max = np.log(k_UV / k_IR)

print("=" * 65)
print("  ERG (Wetterich) DLA POLA TGP -- SZKIELET OP-1")
print("=" * 65)
print(f"  Phi0 = {Phi0}")
print(f"  gamma = {gamma_TGP}")
print(f"  a_Gamma = {a_Gamma}")
print(f"  k_UV = 1/a_Gamma = {k_UV:.1f}")
print(f"  k_IR = sqrt(gamma) = {k_IR:.1f}")
print(f"  Zakres RG: t in [0, {t_max:.2f}]")
print()

# ============================================================
# 2. LPA: Dyskretyzacja potencjalu V_k(psi)
# ============================================================
# Potencjal V w zmiennej psi = Phi/Phi0:
# V(psi) = psi^3/3 - psi^4/4  (beta = gamma, znormalizowane)
# V'(psi) = psi^2 - psi^3 = psi^2*(1 - psi)
# V''(psi) = 2*psi - 3*psi^2

# Rownanie Wetterich w LPA (regulator Litim, d=4):
# dV_k/dk = k^4 / (32*pi^2) * 1 / (V_k''(psi) + k^2)
# lub w zmiennej t = ln(k/k_IR):
# dV_k/dt = k^5 / (32*pi^2) * 1 / (V_k''(psi) + k^2)
# Bo dk/dt = k, wiec dV/dt = (dV/dk) * k

# Dyskretyzacja: V_k(psi) na siatce psi_j, j=0,...,N
N_grid = 200
psi_min, psi_max = 0.01, 3.0
psi_grid = np.linspace(psi_min, psi_max, N_grid)
dpsi = psi_grid[1] - psi_grid[0]

def V_TGP(psi):
    """Potencjal bare TGP (beta=gamma)."""
    return psi**3 / 3 - psi**4 / 4

def V_TGP_pp(psi):
    """V''(psi) = 2*psi - 3*psi^2."""
    return 2*psi - 3*psi**2

def V_second_derivative(V_arr, dpsi):
    """Numeryczna druga pochodna."""
    Vpp = np.zeros_like(V_arr)
    Vpp[1:-1] = (V_arr[2:] - 2*V_arr[1:-1] + V_arr[:-2]) / dpsi**2
    Vpp[0] = Vpp[1]
    Vpp[-1] = Vpp[-2]
    return Vpp

# ============================================================
# 3. Przeplyw RG w LPA
# ============================================================
def rg_rhs(t, V_flat, psi_grid, dpsi, k_IR):
    """Prawa strona rownania Wetterich w LPA."""
    N = len(psi_grid)
    V = V_flat.copy()
    k = k_IR * np.exp(t)

    # V''(psi)
    Vpp = V_second_derivative(V, dpsi)

    # dV/dt = k^5 / (32*pi^2) * 1 / (Vpp + k^2)
    # (w 4D z regulatorem Litim)
    denominator = Vpp + k**2
    # Zabezpieczenie: denominator > 0 (stabilnosc)
    denominator = np.maximum(denominator, 1e-10)

    dVdt = k**5 / (32 * np.pi**2) * 1.0 / denominator
    return dVdt

# Warunek poczatkowy: V_k(t=t_max) = V_TGP(psi) [UV]
# Integrujemy od UV do IR: t biega od t_max do 0
V_UV = V_TGP(psi_grid)

print("--- Przeplyw RG: UV -> IR ---")
print(f"  Siatka: {N_grid} punktow, psi in [{psi_min}, {psi_max}]")
print(f"  Warunek UV: V_TGP(psi) = psi^3/3 - psi^4/4")
print()

# Integracja (od UV = t_max do IR = 0, wiec ujemny krok)
t_span = (t_max, 0.01)
t_eval = np.linspace(t_max, 0.01, 50)

try:
    sol = solve_ivp(
        lambda t, V: rg_rhs(t, V, psi_grid, dpsi, k_IR),
        t_span, V_UV, t_eval=t_eval,
        method='RK45', rtol=1e-6, atol=1e-8, max_step=0.1
    )

    if sol.success:
        print(f"  Integracja: SUKCES ({sol.t.shape[0]} krokow)")
        V_IR = sol.y[:, -1]  # potencjal na skali IR

        # ============================================================
        # 4. Analiza wyniku
        # ============================================================
        print(f"\n--- Analiza potencjalu IR ---")

        # Minimum V_IR
        idx_min = np.argmin(V_IR)
        psi_min_val = psi_grid[idx_min]
        V_min = V_IR[idx_min]
        print(f"  Minimum V_IR: psi = {psi_min_val:.4f}, V = {V_min:.6f}")

        # V''(psi) na minimum (masa^2 pola)
        Vpp_IR = V_second_derivative(V_IR, dpsi)
        m2_IR = Vpp_IR[idx_min]
        print(f"  V''(psi_min) = m^2 = {m2_IR:.6f}")
        if m2_IR > 0:
            print(f"  -> Masa pola: m = {np.sqrt(abs(m2_IR)):.4f} (stabilna)")
        else:
            print(f"  -> Masa^2 < 0: potencjal niestabilny na IR!")

        # Porownanie UV vs IR
        psi_vac = 1.0  # oczekiwane minimum TGP
        idx_vac = np.argmin(np.abs(psi_grid - psi_vac))
        print(f"\n  Na psi = 1 (proznia TGP):")
        print(f"    V_UV(1) = {V_UV[idx_vac]:.6f}")
        print(f"    V_IR(1) = {V_IR[idx_vac]:.6f}")
        print(f"    Zmiana  = {(V_IR[idx_vac] - V_UV[idx_vac]):.6f}")

        # Czy punkt staly UV istnieje?
        # Sprawdzamy: czy potencjal jest ograniczony na UV
        V_UV_max = np.max(np.abs(V_UV))
        V_IR_max = np.max(np.abs(V_IR[psi_grid < 2.0]))
        print(f"\n  Ograniczonosc:")
        print(f"    max|V_UV| = {V_UV_max:.4f}")
        print(f"    max|V_IR| (psi<2) = {V_IR_max:.4f}")

        # ============================================================
        # 5. Wykresy
        # ============================================================
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Potencjal UV vs IR
        ax1 = axes[0]
        ax1.plot(psi_grid, V_UV, 'b-', lw=2, label='V_UV (bare TGP)')
        ax1.plot(psi_grid, V_IR, 'r-', lw=2, label='V_IR (po RG flow)')
        ax1.axvline(1.0, color='gray', ls=':', alpha=0.5)
        ax1.set_xlabel(r'$\psi$')
        ax1.set_ylabel(r'$V_k(\psi)$')
        ax1.set_title('Potencjal: UV vs IR')
        ax1.legend()
        ax1.set_ylim(-2, 2)
        ax1.grid(True, alpha=0.3)

        # Przeplyw V(psi=1) od UV do IR
        ax2 = axes[1]
        V_at_1 = [sol.y[idx_vac, i] for i in range(sol.t.shape[0])]
        ax2.plot(sol.t, V_at_1, 'g-o', markersize=3)
        ax2.set_xlabel('t = ln(k/k_IR)')
        ax2.set_ylabel(r'$V_k(\psi=1)$')
        ax2.set_title('Przeplyw V na prozni')
        ax2.grid(True, alpha=0.3)

        # V'' (masa) od UV do IR
        ax3 = axes[2]
        Vpp_list = []
        for i in range(sol.t.shape[0]):
            Vpp_i = V_second_derivative(sol.y[:, i], dpsi)
            Vpp_list.append(Vpp_i[idx_vac])
        ax3.plot(sol.t, Vpp_list, 'm-o', markersize=3)
        ax3.set_xlabel('t = ln(k/k_IR)')
        ax3.set_ylabel(r"$V''_k(\psi=1)$")
        ax3.set_title('Masa^2 pola: przeplyw RG')
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plt.savefig(os.path.join(script_dir, 'erg_phi_skeleton.png'), dpi=150)
        print(f"\n  Wykres: scripts/erg_phi_skeleton.png")

    else:
        print(f"  Integracja: BLAD ({sol.message})")

except Exception as e:
    print(f"  WYJATEK: {e}")
    import traceback
    traceback.print_exc()

# ============================================================
# 6. Podsumowanie i dalsze kroki
# ============================================================
print(f"\n{'='*65}")
print(f"  PODSUMOWANIE ERG (SZKIELET)")
print(f"{'='*65}")
print(f"""
  Ten skrypt jest SZKIELETEM programu ERG dla TGP.

  Co robi:
  - LPA (Local Potential Approximation) z regulatorem Litim
  - Integracja Wetterich od UV (k=1/a_Gamma) do IR (k=m_sp)
  - Potencjal TGP V(psi) = psi^3/3 - psi^4/4 jako warunek UV

  Co NIE robi (otwarte problemy):
  1. Anomalny wymiar eta_k (potrzebna LPA')
  2. Pelny propagator z K(phi)*phi^4 kinetycznym
  3. Warunki brzegowe na UV (punkt staly AS?)
  4. Wielopolowe obciecia (beyond LPA)
  5. Sektor fermionowy i cechowania
  6. Porownanie z asymptotycznym bezpieczenstwem (AS) grawitacji

  Nastepne kroki:
  OP-1a: LPA' z eta_k (anomalny wymiar)
  OP-1b: Szukanie punktu stalego UV (asymptotic safety)
  OP-1c: Porownanie V_IR z obserwowanym V_eff
  OP-1d: 2-loop ERG (beyond LPA)
""")
print("GOTOWE.")
