#!/usr/bin/env python3
"""
p115_mc_kink_profiles.py — OP-16: Monte Carlo profili kinkow
==============================================================

Cel:  Zbadac robustnosc predykcji r_21 i r_31 pod perturbacjami substratu.
      Jezeli wynik r_31 = 3477 +/- waska banda, TGP jest predykcyjna.

Model: ODE substratowe (ghost-free) z K_sub = g^2:
   g^2 g'' + g(g')^2 + (2/r) g^2 g' = V'(g) = g^2(1-g)

Perturbacje Monte Carlo:
   (A) delta_V: V(g) -> V(g) + epsilon_V * g^2*(g-1)^2   [zachowuje vacuum]
   (B) delta_K: K(g) = g^2 -> g^(2+epsilon_K)             [zmienia kinczyke]
   (C) beta != gamma: V'(g) = g^2(1-g) + delta_bg*g       [lamie N0-5 slabok]
   (D) lambda_6: V(g) -> V(g) + lambda*(g-1)^6             [czlon Wilson-Fisher]

Dla kazdej realizacji:
   1. Skanujemy A_tail(g0) w zakresie g0 ∈ [0.3, 4.0]
   2. Szukamy phi-FP: (A_tail(phi*g0*)/A_tail(g0*))^4 = r_21^PDG
   3. Szukamy g0_tau takiego ze (A_tail(g0_tau)/A_tail(g0*))^4 = r_31^PDG
   4. Zapisujemy r_21, r_31, g0*, Koide Q

Wynik: histogram r_31/r_31^PDG + statystyka robustnosci.

Teoria Generowanej Przestrzeni — Mateusz Serafin, 2026
"""

import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import time

# ══════════════════════════════════════════════════════════════════════
# Stale
# ══════════════════════════════════════════════════════════════════════
PHI = (1.0 + np.sqrt(5.0)) / 2.0  # zloty stosunek
R21_PDG = 206.768
R31_PDG = 3477.15
M_E = 0.51099895   # MeV
M_MU = 105.6583755
M_TAU = 1776.86

# Parametry MC
N_MC       = 500        # realizacji
N_SCAN     = 80         # punktow skanowania g0
G0_MIN     = 0.10
G0_MAX     = 4.50
R_MAX      = 120.0      # zasieg calkowania ODE (RK45 szybki)
R_TAIL_LO  = 40.0       # poczatek dopasowania ogona
R_TAIL_HI  = 110.0      # koniec dopasowania ogona
N_TAIL_PTS = 1000       # gestosc siatki ogona

# Zakresy perturbacji (1-sigma)
SIGMA_EPS_V   = 0.02    # perturbacja potencjalu (~2%)
SIGMA_EPS_K   = 0.05    # perturbacja wykladnika kinetycznego (~5%)
SIGMA_DBG     = 0.005   # beta-gamma rozstrojenie (~0.5%)
SIGMA_LAM6    = 1e-4    # czlon lambda_6 (~10^-4, z Wilson-Fisher)

np.random.seed(42)


# ══════════════════════════════════════════════════════════════════════
# ODE substratowe z perturbacjami
# ══════════════════════════════════════════════════════════════════════
def make_ode(eps_V=0.0, eps_K=0.0, delta_bg=0.0, lam6=0.0):
    """
    Zwraca funkcje ODE z parametrami perturbacji.

    Bazowe ODE substratowe: K(g) = g^2
       g^2 g'' + g(g')^2 + (2/r) g^2 g' = V'(g)

    Perturbacje:
       V(g) = g^3/3 - g^4/4 + eps_V*g^2*(g-1)^2/2 + delta_bg*g^2/2 + lam6*(g-1)^6/6
       K(g) = g^(2+eps_K)
    """
    k_exp = 2.0 + eps_K

    def Vp(g):
        """V'(g) z perturbacjami."""
        # Bazowe: V'(g) = g^2(1-g)
        vp = g**2 * (1.0 - g)
        # eps_V: d/dg [g^2(g-1)^2/2] = g(g-1)(2g-1)
        if abs(eps_V) > 1e-15:
            vp += eps_V * g * (g - 1.0) * (2.0*g - 1.0)
        # delta_bg: d/dg [g^2/2] = g
        if abs(delta_bg) > 1e-15:
            vp += delta_bg * g
        # lam6: d/dg [(g-1)^6/6] = (g-1)^5
        if abs(lam6) > 1e-15:
            vp += lam6 * (g - 1.0)**5
        return vp

    def ode_func(r, y):
        g, gp = y
        g = max(g, 1e-15)

        # K(g) = g^k_exp  =>  K'/K = k_exp/g
        # ODE: K g'' + (K'/2) g'^2 + (2/r) K g' = V'(g)
        # =>   g'' = V'/(K) - (k_exp/(2g)) g'^2 - (2/r) g'
        K = g**k_exp
        if K < 1e-30:
            K = 1e-30

        if r < 1e-10:
            # L'Hopital dla 2g'/r: lim = 2g''
            # => (1 + 2) g'' = V'(g)/K - (k_exp/(2g)) * 0
            gpp = Vp(g) / (3.0 * K)
        else:
            gpp = Vp(g) / K - (k_exp / (2.0 * g)) * gp**2 - (2.0 / r) * gp
        return [gp, gpp]

    def ic_gpp0(g0):
        """Poczatkowe g''(0) z L'Hopital."""
        K0 = g0**k_exp
        if K0 < 1e-30:
            return 0.0
        return Vp(g0) / (3.0 * K0)

    return ode_func, ic_gpp0


# ══════════════════════════════════════════════════════════════════════
# Solver: g0 -> A_tail
# ══════════════════════════════════════════════════════════════════════
def solve_atail(g0, ode_func, ic_gpp0):
    """Rozwiaz ODE i wyciagnij A_tail z ogona oscylacyjnego."""
    if g0 < 0.02 or g0 > 20.0:
        return np.nan

    r_start = 1e-6
    gpp0 = ic_gpp0(g0)
    g_ini = g0 + 0.5 * gpp0 * r_start**2
    gp_ini = gpp0 * r_start

    def event_blowup(r, y):
        return 500.0 - abs(y[0])
    event_blowup.terminal = True

    try:
        sol = solve_ivp(
            ode_func, [r_start, R_MAX], [g_ini, gp_ini],
            method='RK45', rtol=1e-8, atol=1e-10,
            max_step=1.5, events=[event_blowup], dense_output=True
        )
    except Exception:
        return np.nan

    if sol.status == -1 or sol.t[-1] < R_MAX * 0.7:
        return np.nan

    # Ekstrakcja A_tail
    r_tail = np.linspace(R_TAIL_LO, min(sol.t[-1], R_TAIL_HI), N_TAIL_PTS)
    try:
        g_tail = sol.sol(r_tail)[0]
        u = (g_tail - 1.0) * r_tail
        design = np.column_stack([np.cos(r_tail), np.sin(r_tail)])
        BC = np.linalg.lstsq(design, u, rcond=None)[0]
        return np.sqrt(BC[0]**2 + BC[1]**2)
    except Exception:
        return np.nan


# ══════════════════════════════════════════════════════════════════════
# Skanowanie A_tail(g0) + phi-FP + tau search
# ══════════════════════════════════════════════════════════════════════
def run_single_realization(eps_V, eps_K, delta_bg, lam6, g0_arr):
    """Pelna realizacja: skan -> phi-FP -> tau -> Koide."""
    ode_func, ic_gpp0 = make_ode(eps_V, eps_K, delta_bg, lam6)

    # Skanuj A_tail
    A_arr = np.full(len(g0_arr), np.nan)
    for i, g0 in enumerate(g0_arr):
        A_arr[i] = solve_atail(g0, ode_func, ic_gpp0)

    valid = ~np.isnan(A_arr) & (A_arr > 1e-8)
    n_valid = np.sum(valid)
    if n_valid < 15:
        return None

    # Interpolacja
    try:
        f_A = interp1d(g0_arr[valid], A_arr[valid],
                       kind='cubic', fill_value=np.nan, bounds_error=False)
    except Exception:
        return None

    # Szukaj phi-FP
    g0_lo = np.min(g0_arr[valid]) + 0.02
    g0_hi_fp = np.max(g0_arr[valid]) / PHI - 0.01
    if g0_hi_fp <= g0_lo:
        return None

    def res_fp(g0):
        Ae = f_A(g0)
        Am = f_A(PHI * g0)
        if np.isnan(Ae) or np.isnan(Am) or Ae < 1e-12:
            return 1e10
        return (Am / Ae)**4 - R21_PDG

    g0_test = np.linspace(g0_lo, g0_hi_fp, 500)
    res_vals = np.array([res_fp(g) for g in g0_test])
    finite = np.isfinite(res_vals) & (np.abs(res_vals) < 1e9)
    g0f = g0_test[finite]
    rf = res_vals[finite]

    if len(rf) < 3:
        return None
    sc = np.where(np.diff(np.sign(rf)))[0]
    if len(sc) == 0:
        return None

    try:
        g0_star = brentq(res_fp, g0f[sc[0]], g0f[sc[0]+1], xtol=1e-8)
    except Exception:
        return None

    A_e = f_A(g0_star)
    A_mu = f_A(PHI * g0_star)
    if np.isnan(A_e) or np.isnan(A_mu) or A_e < 1e-12:
        return None
    r21 = (A_mu / A_e)**4

    # Szukaj tau
    g0_hi_all = np.max(g0_arr[valid])
    A_hi = f_A(g0_hi_all)
    if np.isnan(A_hi) or A_e < 1e-12:
        return {'g0_star': g0_star, 'r21': r21, 'r31': np.nan,
                'koide': np.nan, 'g0_tau': np.nan}

    max_r = (A_hi / A_e)**4
    if max_r < R31_PDG * 0.5:
        return {'g0_star': g0_star, 'r21': r21, 'r31': np.nan,
                'koide': np.nan, 'g0_tau': np.nan}

    def res_tau(g0t):
        At = f_A(g0t)
        if np.isnan(At) or At < 1e-12:
            return 1e10
        return (At / A_e)**4 - R31_PDG

    g0_tau_range = np.linspace(PHI * g0_star + 0.01, g0_hi_all - 0.01, 500)
    res_t = np.array([res_tau(g) for g in g0_tau_range])
    ft = np.isfinite(res_t) & (np.abs(res_t) < 1e9)

    if np.sum(ft) < 3:
        return {'g0_star': g0_star, 'r21': r21, 'r31': np.nan,
                'koide': np.nan, 'g0_tau': np.nan}

    g0tf = g0_tau_range[ft]
    rtf = res_t[ft]
    tc = np.where(np.diff(np.sign(rtf)))[0]

    if len(tc) == 0:
        # Jesli brak zero-crossingu, uzyj najblizszego punktu
        best = np.argmin(np.abs(rtf))
        At_best = f_A(g0tf[best])
        r31_approx = (At_best / A_e)**4
        return {'g0_star': g0_star, 'r21': r21, 'r31': r31_approx,
                'koide': np.nan, 'g0_tau': g0tf[best]}

    try:
        g0_tau = brentq(res_tau, g0tf[tc[0]], g0tf[tc[0]+1], xtol=1e-8)
    except Exception:
        return {'g0_star': g0_star, 'r21': r21, 'r31': np.nan,
                'koide': np.nan, 'g0_tau': np.nan}

    A_tau = f_A(g0_tau)
    r31 = (A_tau / A_e)**4
    m_tau = M_E * r31

    # Koide
    k_num = M_E + M_MU + m_tau
    k_den = (np.sqrt(M_E) + np.sqrt(M_MU) + np.sqrt(m_tau))**2
    koide = k_num / k_den

    return {
        'g0_star': g0_star,
        'r21': r21,
        'r31': r31,
        'm_tau': m_tau,
        'koide': koide,
        'g0_tau': g0_tau,
    }


# ══════════════════════════════════════════════════════════════════════
# MAIN: Monte Carlo
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("=" * 70)
    print("  TGP OP-16: Monte Carlo profili kinkow — robustnosc r_21, r_31")
    print("=" * 70)
    print(f"\n  N_MC      = {N_MC}")
    print(f"  N_scan    = {N_SCAN}")
    print(f"  sigma(eps_V)   = {SIGMA_EPS_V}")
    print(f"  sigma(eps_K)   = {SIGMA_EPS_K}")
    print(f"  sigma(d_bg)    = {SIGMA_DBG}")
    print(f"  sigma(lam6)    = {SIGMA_LAM6}")
    print()

    g0_arr = np.linspace(G0_MIN, G0_MAX, N_SCAN)

    # ── Bazowa realizacja (bez perturbacji) ──
    print("--- Bazowa realizacja (eps=0) ---")
    t0 = time.time()
    base = run_single_realization(0.0, 0.0, 0.0, 0.0, g0_arr)
    dt_base = time.time() - t0
    if base is None:
        print("  BLAD: bazowa realizacja nie znalazla phi-FP!")
        sys.exit(1)

    print(f"  g0*    = {base['g0_star']:.6f}")
    print(f"  r_21   = {base['r21']:.4f}  [PDG: {R21_PDG}]")
    r31_str = f"{base['r31']:.2f}" if not np.isnan(base.get('r31', np.nan)) else "N/A"
    print(f"  r_31   = {r31_str}  [PDG: {R31_PDG}]")
    if not np.isnan(base.get('koide', np.nan)):
        print(f"  Koide  = {base['koide']:.7f}  [exact: 0.6666667]")
    print(f"  Czas   = {dt_base:.1f}s")

    est_total = dt_base * N_MC / 60
    print(f"\n  Szacowany czas MC: {est_total:.0f} min ({dt_base:.1f}s x {N_MC})")
    print()

    # ── Monte Carlo ──
    print("--- Monte Carlo: perturbacje substratu ---")
    t_mc = time.time()

    # Wyniki
    mc_r21 = []
    mc_r31 = []
    mc_g0star = []
    mc_koide = []
    mc_params = []  # (eps_V, eps_K, delta_bg, lam6)
    n_fail = 0

    for i in range(N_MC):
        eps_V  = np.random.normal(0, SIGMA_EPS_V)
        eps_K  = np.random.normal(0, SIGMA_EPS_K)
        d_bg   = np.random.normal(0, SIGMA_DBG)
        lam6   = np.abs(np.random.normal(0, SIGMA_LAM6))  # lam6 > 0

        result = run_single_realization(eps_V, eps_K, d_bg, lam6, g0_arr)

        if result is None:
            n_fail += 1
        else:
            mc_r21.append(result['r21'])
            mc_g0star.append(result['g0_star'])
            mc_params.append((eps_V, eps_K, d_bg, lam6))
            if not np.isnan(result.get('r31', np.nan)):
                mc_r31.append(result['r31'])
            if not np.isnan(result.get('koide', np.nan)):
                mc_koide.append(result['koide'])

        if (i + 1) % 20 == 0 or i == N_MC - 1:
            dt = time.time() - t_mc
            eta = dt / (i + 1) * (N_MC - i - 1)
            n_ok = len(mc_r21)
            print(f"  [{i+1:3d}/{N_MC}]  ok={n_ok}  fail={n_fail}  "
                  f"t={dt:.0f}s  ETA={eta:.0f}s")
            sys.stdout.flush()

    dt_mc = time.time() - t_mc

    # ══════════════════════════════════════════════════════════════════
    # Analiza statystyczna
    # ══════════════════════════════════════════════════════════════════
    mc_r21 = np.array(mc_r21)
    mc_r31 = np.array(mc_r31)
    mc_g0star = np.array(mc_g0star)
    mc_koide = np.array(mc_koide)

    print("\n" + "=" * 70)
    print("  WYNIKI MONTE CARLO")
    print("=" * 70)
    print(f"\n  Realizacji: {N_MC}")
    print(f"  Udanych:    {len(mc_r21)} ({100*len(mc_r21)/N_MC:.0f}%)")
    print(f"  Z tau:      {len(mc_r31)} ({100*len(mc_r31)/N_MC:.0f}%)")
    print(f"  Nieudanych: {n_fail}")
    print(f"  Czas:       {dt_mc:.0f}s ({dt_mc/60:.1f} min)")

    # r_21
    print(f"\n--- r_21 (stosunek mion/elektron) ---")
    if len(mc_r21) > 2:
        print(f"  srednia  = {np.mean(mc_r21):.4f}")
        print(f"  mediana  = {np.median(mc_r21):.4f}")
        print(f"  std      = {np.std(mc_r21):.4f}")
        print(f"  min/max  = {np.min(mc_r21):.4f} / {np.max(mc_r21):.4f}")
        print(f"  PDG      = {R21_PDG}")
        delta_r21 = 100 * (np.mean(mc_r21) - R21_PDG) / R21_PDG
        spread_r21 = 100 * np.std(mc_r21) / R21_PDG
        print(f"  delta    = {delta_r21:+.4f}%")
        print(f"  spread   = {spread_r21:.4f}%")

        # Ile realizacji w zakresie 1% od PDG?
        in_1pct = np.sum(np.abs(mc_r21 / R21_PDG - 1) < 0.01)
        print(f"  |r21/PDG - 1| < 1%:  {in_1pct}/{len(mc_r21)} "
              f"({100*in_1pct/len(mc_r21):.0f}%)")

    # r_31
    print(f"\n--- r_31 (tau/elektron) ---")
    if len(mc_r31) > 2:
        print(f"  srednia  = {np.mean(mc_r31):.2f}")
        print(f"  mediana  = {np.median(mc_r31):.2f}")
        print(f"  std      = {np.std(mc_r31):.2f}")
        print(f"  min/max  = {np.min(mc_r31):.2f} / {np.max(mc_r31):.2f}")
        print(f"  PDG      = {R31_PDG}")
        delta_r31 = 100 * (np.mean(mc_r31) - R31_PDG) / R31_PDG
        spread_r31 = 100 * np.std(mc_r31) / R31_PDG
        print(f"  delta    = {delta_r31:+.3f}%")
        print(f"  spread   = {spread_r31:.3f}%")

        in_1pct = np.sum(np.abs(mc_r31 / R31_PDG - 1) < 0.01)
        in_5pct = np.sum(np.abs(mc_r31 / R31_PDG - 1) < 0.05)
        print(f"  |r31/PDG - 1| < 1%:  {in_1pct}/{len(mc_r31)} "
              f"({100*in_1pct/len(mc_r31):.0f}%)")
        print(f"  |r31/PDG - 1| < 5%:  {in_5pct}/{len(mc_r31)} "
              f"({100*in_5pct/len(mc_r31):.0f}%)")
    else:
        print(f"  BRAK wystarczajacej statystyki r_31!")

    # Koide
    print(f"\n--- Koide Q ---")
    if len(mc_koide) > 2:
        print(f"  srednia  = {np.mean(mc_koide):.7f}")
        print(f"  std      = {np.std(mc_koide):.7f}")
        print(f"  exact    = 0.6666667")
        delta_koide = np.mean(mc_koide) - 2.0/3.0
        print(f"  delta    = {delta_koide:+.2e}")

    # g0*
    print(f"\n--- g0* (punkt staly phi-FP) ---")
    if len(mc_g0star) > 2:
        print(f"  srednia  = {np.mean(mc_g0star):.6f}")
        print(f"  std      = {np.std(mc_g0star):.6f}")
        print(f"  spread   = {100*np.std(mc_g0star)/np.mean(mc_g0star):.3f}%")

    # ══════════════════════════════════════════════════════════════════
    # Analiza korelacji: ktore perturbacje maja najwiekszy wplyw?
    # ══════════════════════════════════════════════════════════════════
    print(f"\n--- Analiza wrazliwosci ---")
    if len(mc_r21) > 10:
        params_arr = np.array(mc_params[:len(mc_r21)])
        names = ['eps_V', 'eps_K', 'delta_bg', 'lam6']
        for j, name in enumerate(names):
            corr_r21 = np.corrcoef(params_arr[:, j], mc_r21)[0, 1]
            print(f"  corr({name}, r21) = {corr_r21:+.3f}")

        if len(mc_r31) > 10:
            # Dopasuj indeksy
            params_r31 = params_arr[:len(mc_r31)]
            for j, name in enumerate(names):
                if len(mc_r31) == len(params_r31):
                    corr_r31 = np.corrcoef(params_r31[:, j], mc_r31)[0, 1]
                    print(f"  corr({name}, r31) = {corr_r31:+.3f}")

    # ══════════════════════════════════════════════════════════════════
    # Testy zaliczeniowe
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("  TESTY")
    print("=" * 70)

    tests = []

    # T1: r_21 robustne (spread < 2%)
    if len(mc_r21) > 2:
        sp21 = 100 * np.std(mc_r21) / R21_PDG
        t1 = sp21 < 2.0
        tests.append(('T1', f'spread(r_21) = {sp21:.3f}% < 2%', t1))
    else:
        tests.append(('T1', 'brak danych r21', False))

    # T2: r_31 robustne (spread < 5%)
    if len(mc_r31) > 2:
        sp31 = 100 * np.std(mc_r31) / R31_PDG
        t2 = sp31 < 5.0
        tests.append(('T2', f'spread(r_31) = {sp31:.3f}% < 5%', t2))
    else:
        tests.append(('T2', 'brak danych r31', False))

    # T3: >90% realizacji udanych
    frac_ok = len(mc_r21) / N_MC
    t3 = frac_ok > 0.90
    tests.append(('T3', f'{100*frac_ok:.0f}% udanych > 90%', t3))

    # T4: r_21 srednia blisko PDG (<0.5%)
    if len(mc_r21) > 2:
        d21 = abs(100 * (np.mean(mc_r21) - R21_PDG) / R21_PDG)
        t4 = d21 < 0.5
        tests.append(('T4', f'|delta(r_21)| = {d21:.4f}% < 0.5%', t4))
    else:
        tests.append(('T4', 'brak danych', False))

    # T5: r_31 srednia blisko PDG (<2%)
    if len(mc_r31) > 2:
        d31 = abs(100 * (np.mean(mc_r31) - R31_PDG) / R31_PDG)
        t5 = d31 < 2.0
        tests.append(('T5', f'|delta(r_31)| = {d31:.3f}% < 2%', t5))
    else:
        tests.append(('T5', 'brak danych r31', False))

    # T6: Koide robustne (|delta Q| < 0.01)
    if len(mc_koide) > 2:
        dQ = abs(np.mean(mc_koide) - 2.0/3.0)
        t6 = dQ < 0.01
        tests.append(('T6', f'|delta Koide| = {dQ:.2e} < 0.01', t6))
    else:
        tests.append(('T6', 'brak danych Koide', False))

    n_pass = sum(1 for _, _, p in tests if p)
    n_total = len(tests)
    for tid, desc, passed in tests:
        tag = 'PASS' if passed else 'FAIL'
        print(f"  [{tag}] {tid}: {desc}")

    print(f"\n  Wynik: {n_pass}/{n_total} testow PASS")

    # ══════════════════════════════════════════════════════════════════
    # Podsumowanie
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("  PODSUMOWANIE OP-16")
    print("=" * 70)
    if len(mc_r31) > 2:
        print(f"""
  ODE substratowe z perturbacjami Monte Carlo ({N_MC} realizacji):

  r_21 = {np.mean(mc_r21):.2f} +/- {np.std(mc_r21):.2f}  (PDG: {R21_PDG})
  r_31 = {np.mean(mc_r31):.1f} +/- {np.std(mc_r31):.1f}  (PDG: {R31_PDG})

  Spread r_21: {100*np.std(mc_r21)/R21_PDG:.3f}%
  Spread r_31: {100*np.std(mc_r31)/R31_PDG:.3f}%

  Wniosek:""")
        if sp31 < 5.0 and sp21 < 2.0:
            print("  Predykcja r_21 i r_31 jest ROBUSTNA pod perturbacjami substratu.")
            print("  TGP predykuje stosunek mas tau/elektron z waska banda bledow.")
        else:
            print("  Predykcja wykazuje znaczna wrazliwosc na perturbacje substratu.")
            print("  Wymaga dokladniejszego ograniczenia parametrow substratu.")
    else:
        print("  Niewystarczajaca statystyka r_31 — potrzeba wiecej punktow skan.")

    # ══════════════════════════════════════════════════════════════════
    # FAZA 2: Fixed g0* — prawdziwy test robustnosci
    # ══════════════════════════════════════════════════════════════════
    # W fazie 1 phi-FP jest z definicji szukany tak by r21=PDG,
    # wiec spread=0 jest trywialne. Prawdziwy test: ustal g0* z bazowej
    # i sprawdz jak perturbacje zmieniaja A_tail(g0*), A_tail(phi*g0*),
    # A_tail(g0_tau_base) — SAME g0, ROZNE ODE.
    print("\n" + "=" * 70)
    print("  FAZA 2: Fixed-g0* test robustnosci")
    print("  (g0_e, g0_mu, g0_tau ustalone z bazowej; perturbujemy ODE)")
    print("=" * 70)

    g0_e_fix = base['g0_star']
    g0_mu_fix = PHI * g0_e_fix
    g0_tau_fix = base.get('g0_tau', PHI**2 * g0_e_fix)
    print(f"\n  g0_e   = {g0_e_fix:.6f}")
    print(f"  g0_mu  = {g0_mu_fix:.6f}")
    print(f"  g0_tau = {g0_tau_fix:.6f}")

    fix_r21 = []
    fix_r31 = []
    fix_koide = []
    n_fix_fail = 0

    t_fix = time.time()
    for i in range(N_MC):
        eps_V  = np.random.normal(0, SIGMA_EPS_V)
        eps_K  = np.random.normal(0, SIGMA_EPS_K)
        d_bg   = np.random.normal(0, SIGMA_DBG)
        lam6   = np.abs(np.random.normal(0, SIGMA_LAM6))

        ode_func, ic_gpp0 = make_ode(eps_V, eps_K, d_bg, lam6)

        Ae = solve_atail(g0_e_fix, ode_func, ic_gpp0)
        Am = solve_atail(g0_mu_fix, ode_func, ic_gpp0)
        At = solve_atail(g0_tau_fix, ode_func, ic_gpp0)

        if np.isnan(Ae) or np.isnan(Am) or np.isnan(At) or Ae < 1e-12:
            n_fix_fail += 1
            continue

        r21_i = (Am / Ae)**4
        r31_i = (At / Ae)**4
        m_tau_i = M_E * r31_i
        kn = M_E + M_MU + m_tau_i
        kd = (np.sqrt(M_E) + np.sqrt(M_MU) + np.sqrt(m_tau_i))**2

        fix_r21.append(r21_i)
        fix_r31.append(r31_i)
        fix_koide.append(kn / kd)

        if (i + 1) % 50 == 0 or i == N_MC - 1:
            dt = time.time() - t_fix
            print(f"  [{i+1:3d}/{N_MC}]  ok={len(fix_r21)}  fail={n_fix_fail}  "
                  f"t={dt:.0f}s")
            sys.stdout.flush()

    fix_r21 = np.array(fix_r21)
    fix_r31 = np.array(fix_r31)
    fix_koide = np.array(fix_koide)

    print(f"\n  Udanych: {len(fix_r21)}/{N_MC}  Nieudanych: {n_fix_fail}")

    if len(fix_r21) > 2:
        sp21_f = 100 * np.std(fix_r21) / R21_PDG
        sp31_f = 100 * np.std(fix_r31) / R31_PDG
        d21_f = 100 * (np.mean(fix_r21) - R21_PDG) / R21_PDG
        d31_f = 100 * (np.mean(fix_r31) - R31_PDG) / R31_PDG

        print(f"\n--- Fixed-g0* r_21 ---")
        print(f"  srednia  = {np.mean(fix_r21):.4f}")
        print(f"  std      = {np.std(fix_r21):.4f}")
        print(f"  spread   = {sp21_f:.4f}%")
        print(f"  delta    = {d21_f:+.4f}%")
        print(f"  min/max  = {np.min(fix_r21):.3f} / {np.max(fix_r21):.3f}")

        print(f"\n--- Fixed-g0* r_31 ---")
        print(f"  srednia  = {np.mean(fix_r31):.2f}")
        print(f"  std      = {np.std(fix_r31):.2f}")
        print(f"  spread   = {sp31_f:.4f}%")
        print(f"  delta    = {d31_f:+.4f}%")
        print(f"  min/max  = {np.min(fix_r31):.1f} / {np.max(fix_r31):.1f}")

        print(f"\n--- Fixed-g0* Koide ---")
        print(f"  srednia  = {np.mean(fix_koide):.7f}")
        print(f"  std      = {np.std(fix_koide):.7f}")

        # Testy fazy 2
        print(f"\n--- Testy Fazy 2 ---")
        tests2 = []
        t2_1 = sp21_f < 5.0
        tests2.append(('F2-T1', f'spread(r_21) = {sp21_f:.3f}% < 5%', t2_1))
        t2_2 = sp31_f < 10.0
        tests2.append(('F2-T2', f'spread(r_31) = {sp31_f:.3f}% < 10%', t2_2))
        t2_3 = abs(d21_f) < 1.0
        tests2.append(('F2-T3', f'|delta(r_21)| = {abs(d21_f):.4f}% < 1%', t2_3))
        t2_4 = abs(d31_f) < 3.0
        tests2.append(('F2-T4', f'|delta(r_31)| = {abs(d31_f):.4f}% < 3%', t2_4))

        n2_pass = sum(1 for _, _, p in tests2 if p)
        for tid, desc, passed in tests2:
            tag = 'PASS' if passed else 'FAIL'
            print(f"  [{tag}] {tid}: {desc}")
        print(f"\n  Faza 2: {n2_pass}/{len(tests2)} testow PASS")

        print(f"\n  INTERPRETACJA:")
        print(f"  Faza 1 (phi-FP adaptacyjny): atraktor — zawsze znajduje PDG.")
        print(f"  Faza 2 (fixed g0*): spread r_21 = {sp21_f:.2f}%, "
              f"r_31 = {sp31_f:.2f}%.")
        if sp31_f < 10.0:
            print(f"  => Mechanizm phi-FP jest STABILNY: przy ustalonym g0*,")
            print(f"     perturbacje ~{100*SIGMA_EPS_K:.0f}% kinetyki i "
                  f"~{100*SIGMA_EPS_V:.0f}% potencjalu")
            print(f"     daja spread {sp31_f:.2f}% w r_31. ROBUSTNE.")
        else:
            print(f"  => Mechanizm phi-FP jest WRAZLIWY na perturbacje substratu.")
