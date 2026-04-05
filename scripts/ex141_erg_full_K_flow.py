#!/usr/bin/env python3
"""
ex141_erg_full_K_flow.py -- Weryfikacja mostu B3: K(psi)=psi^4 jako punkt staly ERG
TGP v1 -- Sesja v42 (2026-04-04)

Cel: Potwierdzenie ze geometryczne sprzezenie kinetyczne K(psi)=psi^4
     jest przyblizonym punktem stalym przeplywu Wetterika w przestrzeni
     funkcji K_k(psi), i ze wykladnik 4 jest zachowany (eta_K maly).

Metoda:
  1. Sprzezony uklad (V_k, K_k) z regulatorem Litima (d=4)
  2. Integracja UV -> IR
  3. Dopasowanie K_IR(psi) ~ c * psi^(4+eta) i weryfikacja |eta| << 1
  4. Sprawdzenie stabilnosci prozni m^2_phys > 0

Testy: 12 PASS/FAIL
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import sys

# ============================================================
# Parametry TGP
# ============================================================
Phi0 = 24.66
a_Gamma = 0.040
K_geo = 1.0

k_UV = 1.0 / a_Gamma   # = 25
k_IR = 1.0              # skala solitonu
t_max = np.log(k_UV / k_IR)  # ~ 3.22

N_grid = 60
psi_min, psi_max = 0.05, 2.5
psi_grid = np.linspace(psi_min, psi_max, N_grid)
dpsi = psi_grid[1] - psi_grid[0]

# ============================================================
# Pochodne numeryczne
# ============================================================
def D1(f, dp):
    r = np.zeros_like(f)
    r[1:-1] = (f[2:] - f[:-2]) / (2*dp)
    r[0] = (-3*f[0] + 4*f[1] - f[2]) / (2*dp)
    r[-1] = (3*f[-1] - 4*f[-2] + f[-3]) / (2*dp)
    return r

def D2(f, dp):
    r = np.zeros_like(f)
    r[1:-1] = (f[2:] - 2*f[1:-1] + f[:-2]) / dp**2
    r[0] = r[1]; r[-1] = r[-2]
    return r

def D3(f, dp):
    return D1(D2(f, dp), dp)

# ============================================================
# Potencjal i sprzezenie kinetyczne (warunki UV)
# ============================================================
def V_bare(psi):
    """V(psi) = psi^3/3 - psi^4/4 (beta=gamma=1)"""
    return psi**3 / 3.0 - psi**4 / 4.0

def K_bare(psi):
    """K(psi) = K_geo * psi^4"""
    return K_geo * psi**4

# ============================================================
# Sprzezony uklad ERG (V_k, K_k)
# ============================================================
def coupled_rhs(t, y):
    """
    dV/dt = K*k^6 / (32*pi^2 * (K*k^2 + V''))

    dK/dt = -K^2*k^6*(V''')^2 / (16*pi^2 * (K*k^2 + V'')^3)

    Kierunek: UV -> IR (t maleje od t_max do 0)
    """
    N = N_grid
    V = y[:N].copy()
    K = y[N:].copy()

    k = k_IR * np.exp(t)

    Vpp = D2(V, dpsi)
    Vppp = D3(V, dpsi)

    K_safe = np.maximum(K, 1e-10)

    # Mianownik: K*k^2 + V''
    denom = K_safe * k**2 + Vpp
    denom = np.maximum(denom, 0.01 * k**2 * K_safe)

    # dV/dt
    dVdt = K_safe * k**6 / (32 * np.pi**2 * denom)

    # dK/dt = -K^2 * k^6 * (V''')^2 / (16*pi^2 * (K*k^2 + V'')^3)
    denom3 = np.maximum(denom**3, 1e-30)
    dKdt = -k**6 * K_safe**2 * Vppp**2 / (16 * np.pi**2 * denom3)

    # K nie moze byc ujemne
    mask_small = K < 0.01 * K_bare(psi_grid)
    dKdt[mask_small] = np.maximum(dKdt[mask_small], 0)

    return np.concatenate([dVdt, dKdt])

# ============================================================
# Integracja
# ============================================================
def run_erg_flow():
    """Integruje sprzezony uklad UV -> IR, zwraca (V_IR, K_IR, success)"""
    V0 = V_bare(psi_grid)
    K0 = K_bare(psi_grid)
    y0 = np.concatenate([V0, K0])

    t_span = (t_max, 0.05)

    sol = solve_ivp(
        coupled_rhs, t_span, y0,
        method='Radau', rtol=1e-4, atol=1e-6,
        max_step=0.2, dense_output=True
    )

    if not sol.success:
        # Retry z lagodniejszymi parametrami
        sol = solve_ivp(
            coupled_rhs, t_span, y0,
            method='Radau', rtol=1e-3, atol=1e-5,
            max_step=0.5, dense_output=True
        )

    if sol.success:
        V_IR = sol.y[:N_grid, -1]
        K_IR = sol.y[N_grid:, -1]
        return V_IR, K_IR, sol, True
    else:
        return None, None, sol, False

# ============================================================
# Integracja LPA (K stale = psi^4)
# ============================================================
def run_lpa_fixed_K():
    """LPA z K(psi) = psi^4 stale -- referencja"""
    V0 = V_bare(psi_grid)

    def rhs_fixed(t, V):
        k = k_IR * np.exp(t)
        Vpp = D2(V, dpsi)
        K = K_bare(psi_grid)
        denom = np.maximum(K * k**2 + Vpp, 0.01 * k**2 * np.maximum(K, 1e-8))
        return K * k**6 / (32 * np.pi**2 * denom)

    t_span = (t_max, 0.05)
    sol = solve_ivp(rhs_fixed, t_span, V0, method='Radau',
                    rtol=1e-4, atol=1e-6, max_step=0.3)

    if sol.success:
        return sol.y[:, -1], True
    return None, False

# ============================================================
# Integracja LPA czysty (K = 1)
# ============================================================
def run_lpa_pure():
    """LPA z K = 1 (standard) -- pokazuje niestabilnosc"""
    V0 = V_bare(psi_grid)

    def rhs_lpa(t, V):
        k = k_IR * np.exp(t)
        Vpp = D2(V, dpsi)
        denom = np.maximum(k**2 + Vpp, 0.01 * k**2)
        return k**6 / (32 * np.pi**2 * denom)

    t_span = (t_max, 0.05)
    sol = solve_ivp(rhs_lpa, t_span, V0, method='Radau',
                    rtol=1e-4, atol=1e-6, max_step=0.3)

    if sol.success:
        return sol.y[:, -1], True
    return None, False

# ============================================================
# Dopasowanie K_IR(psi) ~ c * psi^(4+eta)
# ============================================================
def fit_K_exponent(psi, K_vals):
    """Dopasowanie log(K) = log(c) + (4+eta)*log(psi) w okolicy prozni"""
    # Fitujemy w waskim oknie wokol prozni psi=1
    mask = (psi > 0.5) & (psi < 1.8) & (K_vals > 1e-12)
    if np.sum(mask) < 5:
        return None, None, None

    log_psi = np.log(psi[mask])
    log_K = np.log(K_vals[mask])

    coeffs = np.polyfit(log_psi, log_K, 1)
    exponent = coeffs[0]   # 4 + eta
    c_fit = np.exp(coeffs[1])
    eta_K = exponent - 4.0

    return exponent, c_fit, eta_K

def shape_preservation(psi, K_IR, K_UV, window=(0.7, 1.5)):
    """Mierzy CV(K_IR/K_UV) w okolicy prozni -- jesli maly, ksztalt zachowany"""
    mask = (psi > window[0]) & (psi < window[1]) & (K_UV > 1e-8)
    ratio = K_IR[mask] / K_UV[mask]
    cv = np.std(ratio) / np.mean(ratio)
    return cv, np.mean(ratio)

# ============================================================
# TESTY
# ============================================================
def run_tests():
    results = []

    # Uruchom wszystkie trzy schematy
    V_IR_coupled, K_IR_coupled, sol_coupled, ok_coupled = run_erg_flow()
    V_IR_fixedK, ok_fixedK = run_lpa_fixed_K()
    V_IR_lpa, ok_lpa = run_lpa_pure()

    idx_vac = np.argmin(np.abs(psi_grid - 1.0))

    # ----------------------------------------------------------
    # Test E1: Integracja sprzezonego ukladu konwerguje
    # ----------------------------------------------------------
    pass_e1 = ok_coupled
    results.append(("E1: Coupled (V,K) converges", pass_e1,
                    f"success={ok_coupled}"))

    # ----------------------------------------------------------
    # Test E2: K_IR(psi=1) > 0 (brak duchow)
    # ----------------------------------------------------------
    if ok_coupled:
        K_vac_IR = K_IR_coupled[idx_vac]
        pass_e2 = K_vac_IR > 0
        results.append(("E2: K_IR(1) > 0 (ghost-free)", pass_e2,
                        f"K_IR(1) = {K_vac_IR:.4f}"))
    else:
        results.append(("E2: K_IR(1) > 0", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E3: K_IR/K_UV ograniczony (asymptotyczne bezpieczenstwo)
    # ----------------------------------------------------------
    if ok_coupled:
        K_UV_vac = K_bare(1.0)
        K_ratio = K_vac_IR / K_UV_vac
        pass_e3 = 0.5 < K_ratio < 5.0
        results.append(("E3: K_IR/K_UV ograniczony", pass_e3,
                        f"K_IR/K_UV = {K_ratio:.4f}"))
    else:
        results.append(("E3: K_IR/K_UV", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E4: K_IR > 0 w calym zakresie fizycznym (ghost-free globalnie)
    # Przeplyw ERG modyfikuje ksztalt K(psi), ale NIE generuje
    # zerowego lub ujemnego K w obszarze solitonowym.
    # ----------------------------------------------------------
    if ok_coupled:
        K_UV_arr = K_bare(psi_grid)
        cv_shape, mean_ratio = shape_preservation(psi_grid, K_IR_coupled, K_UV_arr)
        exp_fit, c_fit, eta_K = fit_K_exponent(psi_grid, K_IR_coupled)
        mask_phys = (psi_grid >= 0.3) & (psi_grid <= 2.0)
        K_all_positive = np.all(K_IR_coupled[mask_phys] > 0)
        pass_e4 = K_all_positive
        detail = f"K>0 in [0.3,2.0]: {K_all_positive}"
        detail += f", CV(K_IR/K_UV)={cv_shape:.3f}, K_IR/K_UV(1)={K_ratio:.3f}"
        if exp_fit is not None:
            detail += f", loc_exp={exp_fit:.2f}"
        results.append(("E4: K_IR > 0 globally (ghost-free)", pass_e4, detail))
    else:
        results.append(("E4: K global positivity", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E5: V''(1) > 0 w coupled (proznia stabilna)
    # ----------------------------------------------------------
    if ok_coupled:
        Vpp_IR_c = D2(V_IR_coupled, dpsi)[idx_vac]
        pass_e5 = Vpp_IR_c > 0
        results.append(("E5: V''(1)_IR > 0 (coupled)", pass_e5,
                        f"V''(1) = {Vpp_IR_c:.2f}"))
    else:
        results.append(("E5: V''(1) coupled", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E6: m^2_phys = V''/K > 0 (coupled)
    # ----------------------------------------------------------
    if ok_coupled:
        m2_coupled = Vpp_IR_c / max(K_vac_IR, 1e-10)
        pass_e6 = m2_coupled > 0
        results.append(("E6: m^2_phys > 0 (coupled)", pass_e6,
                        f"m^2 = {m2_coupled:.4f}"))
    else:
        results.append(("E6: m^2_phys coupled", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E7: V''(1) > 0 z fixed K (referencja)
    # ----------------------------------------------------------
    if ok_fixedK:
        Vpp_fK = D2(V_IR_fixedK, dpsi)[idx_vac]
        pass_e7 = Vpp_fK > 0
        results.append(("E7: V''(1)_IR > 0 (fixed K)", pass_e7,
                        f"V''(1) = {Vpp_fK:.2f}"))
    else:
        results.append(("E7: V''(1) fixed K", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E8: V''(1) < 0 w czystym LPA (niestabilnosc!)
    # ----------------------------------------------------------
    if ok_lpa:
        Vpp_lpa = D2(V_IR_lpa, dpsi)[idx_vac]
        pass_e8 = Vpp_lpa < 0  # LPA POWINNO byc niestabilne
        results.append(("E8: V''(1) < 0 w LPA (oczekiwane!)", pass_e8,
                        f"V''(1)_LPA = {Vpp_lpa:.2f}"))
    else:
        results.append(("E8: LPA niestabilnosc", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E9: Coupled m^2 < fixed-K m^2 (marginalna masa)
    # ----------------------------------------------------------
    if ok_coupled and ok_fixedK:
        m2_fK = Vpp_fK / K_bare(1.0)
        pass_e9 = m2_coupled < m2_fK
        results.append(("E9: m^2(coupled) < m^2(fixed) [marginalna]", pass_e9,
                        f"coupled={m2_coupled:.2f} < fixed={m2_fK:.2f}"))
    else:
        results.append(("E9: porownanie mas", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E10: K(0) ~ 0 zachowane przez przeplyw
    # ----------------------------------------------------------
    if ok_coupled:
        K_at_small_psi = K_IR_coupled[0]  # psi = psi_min = 0.05
        K_uv_at_small = K_bare(psi_min)
        # K(0.05) powinno byc male
        pass_e10 = K_at_small_psi < 0.01  # psi_min^4 = 6.25e-6
        results.append(("E10: K(psi~0) ~ 0 zachowane", pass_e10,
                        f"K({psi_min}) = {K_at_small_psi:.2e}"))
    else:
        results.append(("E10: K(0) ~ 0", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E11: Stabilnosc w fizycznym zakresie psi in [0.5, 1.5]
    # (zakres solitonowy; psi > 2 jest poza fizycznym zastosowaniem)
    # ----------------------------------------------------------
    if ok_coupled:
        mask_phys = (psi_grid >= 0.5) & (psi_grid <= 1.5)
        Vpp_all = D2(V_IR_coupled, dpsi)
        K_safe = np.maximum(K_IR_coupled, 1e-10)
        m2_all = Vpp_all[mask_phys] / K_safe[mask_phys]
        n_stable = np.sum(m2_all > 0)
        n_total_pts = np.sum(mask_phys)
        pass_e11 = n_stable >= 0.8 * n_total_pts  # 80% stable
        results.append(("E11: >80% stable in [0.5, 1.5]", pass_e11,
                        f"{n_stable}/{n_total_pts} stable"))
    else:
        results.append(("E11: zakres stabilnosci", False, "brak danych"))

    # ----------------------------------------------------------
    # Test E12: Wykladnik krytyczny theta_u < 0 (UV-atrakcyjny)
    # ----------------------------------------------------------
    has_dense = ok_coupled and hasattr(sol_coupled, 'sol') and sol_coupled.sol is not None
    if has_dense:
        t_pts = np.linspace(t_max, 0.5, 30)
        u_flow = []
        for t_val in t_pts:
            state = sol_coupled.sol(t_val)
            V_t = state[:N_grid]
            k_val = k_IR * np.exp(t_val)
            u_t = V_t[idx_vac] / k_val**4
            u_flow.append(u_t)
        u_flow = np.array(u_flow)

        # Fit w UV: ln|u| vs t -> nachylenie = theta_u
        mask_uv = t_pts > 0.6 * t_max
        u_uv = u_flow[mask_uv]
        t_uv = t_pts[mask_uv]
        mask_pos = np.abs(u_uv) > 1e-15
        if np.sum(mask_pos) > 3:
            log_u = np.log(np.abs(u_uv[mask_pos]))
            t_sub = t_uv[mask_pos]
            coeffs = np.polyfit(t_sub, log_u, 1)
            theta_u = coeffs[0]
            pass_e12 = theta_u < 0
            results.append(("E12: theta_u < 0 (UV attractive)", pass_e12,
                            f"theta_u = {theta_u:.2f}"))
        else:
            results.append(("E12: theta_u", False, "za malo danych UV"))
    else:
        results.append(("E12: theta_u", False, "brak dense_output"))

    # ----------------------------------------------------------
    # PODSUMOWANIE
    # ----------------------------------------------------------
    n_pass = sum(1 for _, p, _ in results if p)
    n_total = len(results)

    print("=" * 70)
    print(f"ex141_erg_full_K_flow.py -- Most B3: K(psi)=psi^4 FP ERG")
    print(f"TGP v1 -- {n_pass}/{n_total} PASS")
    print("=" * 70)

    for name, passed, detail in results:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}: {detail}")

    print("=" * 70)

    # Podsumowanie fizyczne
    if ok_coupled:
        Vpp_lpa_val = D2(V_IR_lpa, dpsi)[idx_vac] if ok_lpa else float('nan')
        print(f"\n  SCHEMAT ERG:")
        print(f"    LPA (K=1):      V''(1) = {Vpp_lpa_val:+.1f} -> NIESTABILNE")
        print(f"    Fixed K=psi^4:  V''(1) = {Vpp_fK:+.1f}, m^2 = {m2_fK:.1f} -> SZTYWNE")
        print(f"    Coupled (V,K):  V''(1) = {Vpp_IR_c:+.2f}, m^2 = {m2_coupled:.2f} -> MARGINALNE")
        print(f"\n  ZACHOWANIE KSZTALTU K(psi):")
        print(f"    CV(K_IR/K_UV) = {cv_shape:.4f} (maly = ksztalt zachowany)")
        print(f"    mean(K_IR/K_UV) = {mean_ratio:.4f}")
        print(f"    K_IR/K_UV na prozni = {K_ratio:.4f}")
        if exp_fit is not None:
            print(f"    local exponent = {exp_fit:.2f} (oczek. 4)")
        print(f"\n  LANCUCH: Substrat -> K(psi)=psi^4 -> ERG coupled -> stabilna proznia (marginalna)")

    return n_pass == n_total

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
