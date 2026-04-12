"""
EX73: CZY alpha*_1 + alpha*_2 = 2pi - alpha_TGP/2 - beta/10 JEST PRAWEM?
STATUS: LEGACY-TRANSLATIONAL

This script is part of the older `alpha*` numerology / sum-rule exploration.
Treat it as historical exploratory context, not as synchronized canonical
`nbody`.

Idea:
  - W ex71/ex72 skanowalismy KINETYCZNY parametr alpha ODE i znalezlismy
    dwa zera delta(alpha) = g0*(alpha) - z0(alpha) przy:
      alpha*_1 = 2.44143, alpha*_2 = 2.74175
    ich suma = 2pi - 11/10 = 2pi - alpha_TGP/2 - beta_pot/10
    gdzie alpha_TGP=2 (bazowy) i beta_pot=1 (wspolczynnik potencjalu).

  - PYTANIE: Jesli zmieniamy beta_pot (wspolczynnik V'=beta*g^2*(1-g)),
    czy suma alpha*_1 + alpha*_2 przesuwa sie zgodnie z:
      suma = 2pi - 1 - beta_pot/10   ?
    (1 = alpha_TGP/2 = 2/2 = stale)

  - Mozliwy wynik:
    A) TAK: formula jest prawem ogolnym (nie numerologia!)
    B) NIE: formuła jest specyficzna dla beta=1

Schemat:
  Dla kazdej wartosci beta_pot: skanuj alpha (kinetyczny) w [2.0, 3.5]
  => znajdz dwa zera delta(alpha, beta_pot)
  => sprawdz czy ich suma = 2pi - 1 - beta_pot/10
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

PHI       = (1.0 + np.sqrt(5.0)) / 2.0
R21_EXP   = 206.768
G_OFF     = 0.005
R_MAX     = 100.0
WIN_LIST  = [16, 22, 28, 36, 46, 58, 72]
WIN_WIDTH = 14.0
N_CORES   = min(16, cpu_count())
ALPHA_TGP = 2.0   # bazowy kinetyczny parametr TGP

# ============================================================
# SOLVER: ODE z kinetycznym alpha i potencjalnym beta_pot
# ============================================================
def _integrate(g0, alpha_kin, beta_pot, maxb=8):
    """Calkuje ODE z kinetycznym alpha=alpha_kin i V'=beta_pot*g^2*(1-g)."""
    gb = np.exp(-1.0 / (2.0 * alpha_kin)) + G_OFF
    def rhs(r, y):
        g, gp = y
        g = max(g, gb + 1e-7)
        fg = 1.0 + 2.0 * alpha_kin * np.log(g)
        if abs(fg) < 1e-10: return [gp, 0.0]
        dr = beta_pot * g**2 * (1.0 - g)
        cr = (alpha_kin / g) * gp**2
        if r < 1e-10: return [gp, (dr - cr) / (3.0 * fg)]
        return [gp, (dr - cr - fg * 2.0 * gp / r) / fg]
    def ev(r, y): return y[0] - gb
    ev.terminal = True; ev.direction = -1
    y0 = [g0, 0.0]; r0 = 1e-10; ra, ga = [], []
    for _ in range(maxb + 1):
        sol = solve_ivp(rhs, (r0, R_MAX), y0, events=ev,
                        dense_output=True, rtol=1e-9, atol=1e-11, max_step=0.05)
        ra.append(sol.t); ga.append(sol.y[0])
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            rh = sol.t_events[0][-1]; st = sol.sol(rh)
            y0 = [st[0], -st[1]]; r0 = rh
        else:
            break
    r = np.concatenate(ra); g = np.concatenate(ga)
    idx = np.argsort(r)
    return r[idx], g[idx]

def _fit_win(r, g, rL, rR):
    mask = (r >= rL) & (r <= rR)
    if mask.sum() < 20: return 0., 0., 0.
    rf = r[mask]; df = (g[mask] - 1.0) * rf
    M = np.column_stack([np.cos(rf), np.sin(rf)])
    bc = np.linalg.lstsq(M, df, rcond=None)[0]
    return float(np.sqrt(bc[0]**2 + bc[1]**2)), float(bc[0]), float(bc[1])

def _A_inf(g0, alpha_kin, beta_pot):
    r, g = _integrate(g0, alpha_kin, beta_pot)
    Av, rv = [], []
    for rL in WIN_LIST:
        if rL + WIN_WIDTH > r[-1]: break
        A, _, _ = _fit_win(r, g, rL, rL + WIN_WIDTH)
        if A > 0.002: Av.append(A); rv.append(float(rL))
    if not Av: return 0.
    if len(Av) < 3: return Av[-1]
    try:
        p, _ = curve_fit(lambda x, ai, a: ai * (1 + a / x), rv, Av,
                         p0=[Av[-1], 0.], maxfev=2000)
        return float(p[0])
    except:
        return float(Av[-1])

def _B_coeff(g0, alpha_kin, beta_pot):
    r, g = _integrate(g0, alpha_kin, beta_pot)
    _, B, _ = _fit_win(r, g, WIN_LIST[2], WIN_LIST[2] + WIN_WIDTH)
    return B

def _f_ratio(g0, alpha_kin, beta_pot):
    Ae = _A_inf(g0, alpha_kin, beta_pot)
    Am = _A_inf(PHI * g0, alpha_kin, beta_pot)
    return (Am / Ae)**4 if Ae > 1e-6 else np.nan

def _find_z0(alpha_kin, beta_pot, lo=1.05, hi=2.00):
    gs = np.linspace(lo, hi, 40)
    Bs = [_B_coeff(g, alpha_kin, beta_pot) for g in gs]
    for i in range(len(Bs) - 1):
        if Bs[i] * Bs[i+1] < 0:
            try:
                return brentq(lambda x: _B_coeff(x, alpha_kin, beta_pot),
                              gs[i], gs[i+1], xtol=1e-7)
            except:
                pass
    return None

def _find_g0star(alpha_kin, beta_pot, target=R21_EXP, lo=1.05, hi=2.00):
    gs = np.linspace(lo, hi, 25)
    fs = [_f_ratio(g, alpha_kin, beta_pot) for g in gs]
    for i in range(len(fs) - 1):
        if not (np.isnan(fs[i]) or np.isnan(fs[i+1])):
            if (fs[i] - target) * (fs[i+1] - target) < 0:
                try:
                    return brentq(lambda x: _f_ratio(x, alpha_kin, beta_pot) - target,
                                  gs[i], gs[i+1], xtol=1e-7)
                except:
                    pass
    return None

def delta_scan(alpha_kin, beta_pot):
    """delta = g0*(alpha_kin, beta_pot) - z0(alpha_kin, beta_pot)."""
    z0  = _find_z0(alpha_kin, beta_pot)
    g0s = _find_g0star(alpha_kin, beta_pot)
    if z0 and g0s:
        return g0s - z0
    return np.nan

# ============================================================
# Dla ustalonego beta_pot: znajdz dwa zera delta(alpha_kin)
# ============================================================
def find_two_zeros(beta_pot, a_lo=2.0, a_hi=3.6, n_scan=28):
    """Skanuj alpha_kin szukajac zer delta(alpha_kin, beta_pot)."""
    a_arr = np.linspace(a_lo, a_hi, n_scan)
    deltas = [delta_scan(a, beta_pot) for a in a_arr]
    zeros = []
    for i in range(len(a_arr) - 1):
        d1, d2 = deltas[i], deltas[i+1]
        if not (np.isnan(d1) or np.isnan(d2)) and d1 * d2 < 0:
            try:
                z = brentq(lambda a: delta_scan(a, beta_pot),
                           a_arr[i], a_arr[i+1], xtol=1e-6)
                zeros.append(z)
            except:
                pass
    return zeros

# ============================================================
# Funkcja robocza dla jednego beta_pot
# ============================================================
def compute_case(beta_pot):
    zeros = find_two_zeros(beta_pot)
    formula = 2 * np.pi - ALPHA_TGP / 2.0 - beta_pot / 10.0
    if len(zeros) >= 2:
        a1, a2 = zeros[0], zeros[1]
        suma = a1 + a2
        diff = suma - formula
        return (beta_pot, a1, a2, suma, formula, diff, True)
    elif len(zeros) == 1:
        return (beta_pot, zeros[0], np.nan, np.nan, formula, np.nan, False)
    else:
        return (beta_pot, np.nan, np.nan, np.nan, formula, np.nan, False)

# ============================================================
# GLOWNY PROGRAM
# ============================================================
if __name__ == '__main__':
    print("=" * 68)
    print("EX73: SUMA alpha*_1 + alpha*_2 vs. beta_pot")
    print(f"  Testowana formula: suma = 2pi - {ALPHA_TGP}/2 - beta_pot/10")
    print(f"                           = 2pi - 1.0 - beta_pot/10")
    print("=" * 68)
    print(f"  Uzywam {N_CORES} rdzeni")
    print()

    beta_list = [0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]
    print(f"  Wartosci beta_pot: {beta_list}")
    print()

    with Pool(processes=N_CORES) as pool:
        results = pool.map(compute_case, beta_list)

    # Naglowek
    print(f"  {'beta_pot':>8}  {'a*_1':>9}  {'a*_2':>9}  "
          f"{'suma':>10}  {'2pi-1-b/10':>11}  {'diff':>10}  {'OK?':>4}")
    print(f"  {'-'*8}  {'-'*9}  {'-'*9}  "
          f"{'-'*10}  {'-'*11}  {'-'*10}  {'-'*4}")

    n_ok  = 0
    n_two = 0
    formula_diffs = []
    betas_two = []
    sums_two  = []
    formulas_two = []

    for r in results:
        beta_pot, a1, a2, suma, formula, diff, two = r
        if two:
            n_two += 1
            ok = not np.isnan(diff) and abs(diff) < 0.02
            if ok: n_ok += 1
            if not np.isnan(diff):
                formula_diffs.append(abs(diff))
                betas_two.append(beta_pot)
                sums_two.append(suma)
                formulas_two.append(formula)
            a1s   = f"{a1:.6f}"
            a2s   = f"{a2:.6f}"
            sums  = f"{suma:.7f}"
            fors  = f"{formula:.7f}"
            diffs = f"{diff:+.3e}" if not np.isnan(diff) else "N/A"
            oks   = "PASS" if ok else "FAIL"
        else:
            a1s = f"{a1:.6f}" if not np.isnan(a1) else "N/A"
            a2s = "N/A"; sums = "N/A"
            fors  = f"{formula:.7f}"
            diffs = "N/A"; oks = "---"
        print(f"  {beta_pot:8.2f}  {a1s:>9}  {a2s:>9}  "
              f"{sums:>10}  {fors:>11}  {diffs:>10}  {oks:>4}")

    print()
    print(f"  Przypadki z 2 zerami: {n_two}/{len(beta_list)}")
    print(f"  Zgodnosc (<0.02):     {n_ok}/{n_two}")

    # ============================================================
    # Sekcja 2: Trend sumy vs. beta_pot
    # ============================================================
    if len(betas_two) >= 3:
        print()
        print("-" * 68)
        print("--- Sekcja 2: Trend suma(beta_pot) ---")
        print()
        betas_arr  = np.array(betas_two)
        sums_arr   = np.array(sums_two)
        forms_arr  = np.array(formulas_two)

        # Dopasowanie liniowe suma = A + B*beta_pot
        from scipy.stats import linregress
        slope, intercept, r2, pval, se = linregress(betas_arr, sums_arr)
        print(f"  suma(beta) = {intercept:.6f} + {slope:.6f}*beta_pot")
        print(f"  R^2 = {r2**2:.6f},  oczekiwany slope = -1/10 = {-0.1:.6f}")
        print(f"  Rzeczywisty slope = {slope:.6f}   "
              f"diff od -0.1: {abs(slope+0.1):.4e}")
        print()
        print(f"  {'beta_pot':>8}  {'suma':>10}  {'formula':>10}  {'diff':>10}  {'suma-formula':>12}")
        for b, s, f in zip(betas_arr, sums_arr, forms_arr):
            print(f"  {b:8.2f}  {s:10.7f}  {f:10.7f}  {s-f:+10.5f}  "
                  f"{'CLOSE' if abs(s-f)<0.01 else 'FAR':>12}")

    # ============================================================
    # TESTY
    # ============================================================
    print()
    print("=" * 68)
    print("TESTY")
    print("=" * 68)

    T1 = n_two >= 5
    print(f"  T1: >= 5 przypadkow z 2 zerami: "
          f"{'PASS' if T1 else 'FAIL'}  ({n_two}/{len(beta_list)})")

    T2 = n_ok >= max(1, n_two - 1)
    print(f"  T2: >= n_two-1 PASS (<0.02):    "
          f"{'PASS' if T2 else 'FAIL'}  ({n_ok}/{n_two})")

    T3 = False
    if len(betas_two) >= 3:
        T3 = abs(slope + 0.1) < 0.02   # slope ≈ -0.1
        print(f"  T3: slope ~ -0.1 (diff<0.02):  "
              f"{'PASS' if T3 else 'FAIL'}  (slope={slope:.5f})")

    T4 = False
    # Bazowy beta=1.0 powinien miec diff < 1e-5
    base = [r for r in results if r[0] == 1.0]
    if base and base[0][6]:
        diff_base = base[0][5]
        T4 = not np.isnan(diff_base) and abs(diff_base) < 1e-4
        print(f"  T4: beta=1 diff < 1e-4:        "
              f"{'PASS' if T4 else 'FAIL'}  ({diff_base:.3e})")

    T5 = False
    if len(formula_diffs) >= 3:
        # Formula zawsze lepsza niz prosta stala (2pi - 1.1)
        stala_diffs = [abs(s - (2*np.pi - 1.1)) for s in sums_two]
        T5 = np.mean(formula_diffs) < np.mean(stala_diffs)
        print(f"  T5: formula lepsza niz stala 2pi-1.1: "
              f"{'PASS' if T5 else 'FAIL'}  "
              f"(mean_formula={np.mean(formula_diffs):.4e} "
              f"vs mean_stala={np.mean(stala_diffs):.4e})")

    n_pass = sum([T1, T2, T3, T4, T5])
    print(f"\nWYNIK: {n_pass}/5 testow przeszlo")

    # ============================================================
    # WNIOSEK
    # ============================================================
    print()
    print("=" * 68)
    print("WNIOSEK (EX73)")
    print("=" * 68)
    if T3 and T2:
        print("  >> PRAWO POTWIERDZONE!")
        print(f"     suma(alpha*_1 + alpha*_2) = 2pi - {ALPHA_TGP}/2 - beta_pot/10")
        print(f"     slope(suma, beta_pot) = {slope:.5f} ≈ -1/10 = -0.10000")
        print(f"     Formula jest ogolna — nie numerologia!")
    elif T3:
        print(f"  >> Slope = {slope:.5f} ≈ -0.1 (T3 PASS), ale precyzja nizsza.")
        print(f"     Formula jakosciowo poprawna.")
    else:
        print(f"  >> Slope = {slope:.5f if len(betas_two)>=3 else 'N/A'}, oczekiwano ≈ -0.1.")
        print(f"     Formula 2pi - 1 - beta/10 NIE jest ogolnym prawem.")
    print("=" * 68)
