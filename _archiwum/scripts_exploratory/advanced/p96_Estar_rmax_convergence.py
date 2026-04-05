#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p96_Estar_rmax_convergence.py  --  TGP v1
==========================================
Cel: Zbadac, czy minimum E*(alpha) = argmin[E*(alpha)]
     jest fizyczne (odporne na r_max) czy artefaktem r_max=40.

Kontekst:
  p82-p84: argmin E*(alpha, r_max=40) ~ 8.60-8.65, offset ~0.04-0.09 od alpha_K=8.5616
  p93:     K*(r_max) = K*(inf) + c/r_max  (konwergencja 1/r)
  p95:     K_nan_R(r_max=40) = artefakt => alpha_max nie fizyczny
  => Pytanie: czy argmin E* rowniez jest artefaktem?

Metoda:
  Skan alpha in [8.20, 8.70], 11 punktow
  Dla kazdego alpha: K*(alpha, r_max) przy r_max = 40, 60, 80
  E*(alpha, r_max) = 4*pi*K*(alpha, r_max)  [poniewaz g(K*)=0 => E=4piK*]
  Ekstrapolacja: K*(alpha, inf) z fitu K*(r) = a + b/r (dwa punkty r=60,80)
  argmin E*(alpha, r_max=40) vs argmin E*(alpha, r_max=inf)

  Dodatkowe sprawdzenie:
  - czy E*(alpha) jest monotoniczna (brak minimum) czy ma minimum?
  - polozenie minimum vs alpha_K=8.5616

Testy:
  P1: argmin E*(alpha, r_max=40) w [8.55, 8.70]  (zgodne z p82-p83)
  P2: argmin E*(alpha, r_max=60) blizej alpha_K niz argmin(r=40)
  P3: argmin E*(alpha, inf) ~ alpha_K (odchylenie < 0.05)
  P4: E*(alpha) monotonicznie maleje dla alpha < argmin => jest minimum
  P5: K*(alpha, inf) wykazuje minimum (nie monotoniczne)

Data: 2026-03-26
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import linregress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K    = 8.5616
A_GAM      = 0.040
LAM        = 5.501357e-06
GAMMA_V    = 1.0
V1         = GAMMA_V/3 - GAMMA_V/4

# Skan
ALPHAS = np.round(np.linspace(8.20, 8.70, 11), 4)  # [8.20, 8.25, ..., 8.70]
RMAXES = [40.0, 60.0, 80.0]

# Parametry numeryczne (balans szybkosc/dokladnosc)
N_PSI  = 80
N_EVAL = 500
N_K    = 25

# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):  return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, alpha, r_max):
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        return float(sol.y[0,-1]) if sol.t[-1] >= r_max*0.99 else np.nan
    except:
        return np.nan

def g_func(K, alpha, r_max):
    """g = E/(4piK) - 1 dla Galezy U (ostatnie zero psi w [1.001,3.0])."""
    psis = np.linspace(1.001, 3.0, N_PSI)
    Fv = [phi_at_rmax(p, K, alpha, r_max) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,alpha,r_max)-1.0,
                            psis[i], psis[i+1], xtol=1e-7, maxiter=40)
                roots.append(pz)
            except: pass
    if not roots: return np.nan, np.nan
    psi_z = roots[-1]
    # Energia
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0, 1, N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM, r_max],
                        [psi_z, dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
        Ek=4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep=4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return (Ek+Ep)/(4*np.pi*K)-1.0, psi_z
    except:
        return np.nan, psi_z

def find_K_star(alpha, r_max):
    """Szuka K* (zero g_func). Zwraca (K*, E*, psi*)."""
    Ks = np.linspace(0.006, 0.018, N_K)
    gv = []
    for K in Ks:
        g, _ = g_func(K, alpha, r_max)
        gv.append(g)
    K_star = np.nan; E_star = np.nan; psi_star = np.nan
    for i in range(len(gv)-1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i+1]) and gv[i]*gv[i+1] < 0:
            try:
                Kz = brentq(lambda K: g_func(K, alpha, r_max)[0],
                            Ks[i], Ks[i+1], xtol=5e-6, maxiter=20)
                K_star = Kz
                _, psi_star = g_func(Kz, alpha, r_max)
                # E* = 4*pi*K* (bo g=0 => E=4piK)
                # Ale liczymy E bezposrednio
                dphi0 = -K_star/A_GAM**2
                r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0,1,N_EVAL)
                sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM,r_max],
                                [psi_star, dphi0], method='DOP853',
                                rtol=1e-8, atol=1e-10, t_eval=r_ev)
                r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
                Ek=4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
                Ep=4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
                E_star = Ek + Ep
            except:
                K_star = np.nan
    return K_star, E_star, psi_star

# ─────────────────────────────────────────────────────────────────────────────
def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    print("=" * 68)
    print("TGP v1  *  p96_Estar_rmax_convergence.py")
    print(f"  alpha_K={ALPHA_K}  (OP-3 Sciezka 1: argmin E*(alpha) = alpha_K?)")
    print(f"  alphas = {list(ALPHAS)}")
    print(f"  r_max = {RMAXES},  n_psi={N_PSI}, n_K={N_K}, n_eval={N_EVAL}")
    print("=" * 68)

    PC=0; FC=0
    def record(label, ok, info=""):
        nonlocal PC,FC
        if ok: PC+=1
        else: FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    # ─────────────────────────────────────────────────────────────────────────
    # Tabela wynikow: results[alpha][r_max] = (K*, E*, psi*)
    results = {a: {} for a in ALPHAS}

    for alpha in ALPHAS:
        print(f"\n-- alpha={alpha:.4f} --")
        print(f"  {'r_max':>6}  {'K*':>12}  {'E*':>12}  {'psi*':>8}  t")
        for r_max in RMAXES:
            t0 = time.time()
            Ks, Es, ps = find_K_star(alpha, r_max)
            dt = time.time()-t0
            ks_s = f'{Ks:.7f}' if np.isfinite(Ks) else '   brak   '
            es_s = f'{Es:.6f}' if np.isfinite(Es) else '   brak  '
            ps_s = f'{ps:.4f}' if np.isfinite(ps) else ' brak '
            print(f"  {r_max:>6.0f}  {ks_s}  {es_s}  {ps_s}  {dt:.0f}s")
            results[alpha][r_max] = (Ks, Es, ps)

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== TABELA K*(alpha, r_max) ==")
    print(f"  {'alpha':>7}", end="")
    for r in RMAXES:
        print(f"  {'K*(r='+str(int(r))+')':>14}", end="")
    print(f"  {'K*(inf)':>14}  {'c':>10}")
    print("  " + "-"*80)

    K_inf = {}   # alpha -> K*(inf)
    for alpha in ALPHAS:
        print(f"  {alpha:.4f}", end="")
        Kvals = [(r, results[alpha][r][0]) for r in RMAXES if np.isfinite(results[alpha].get(r,(np.nan,)*3)[0])]
        for r in RMAXES:
            Ks = results[alpha].get(r, (np.nan,)*3)[0]
            ks_s = f'{Ks:.7f}' if np.isfinite(Ks) else '    brak   '
            print(f"  {ks_s}", end="")
        # Fit K*(r) = a + b/r z co najmniej 2 punktow
        if len(Kvals) >= 2:
            rs = np.array([x[0] for x in Kvals])
            Ks_arr = np.array([x[1] for x in Kvals])
            # Fit: K = a + b/r => K * r = a*r + b => linear in [r, 1]
            X = np.column_stack([rs, np.ones_like(rs)])
            # least squares: K = a + b/r
            inv_r = 1.0/rs
            X2 = np.column_stack([np.ones_like(rs), inv_r])
            from numpy.linalg import lstsq
            coef, _, _, _ = lstsq(X2, Ks_arr, rcond=None)
            a_fit, b_fit = coef
            K_inf[alpha] = a_fit
            print(f"  {a_fit:.7f}  {b_fit:.5f}")
        else:
            K_inf[alpha] = np.nan
            print(f"  {'brak':>14}  {'brak':>10}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== E*(alpha) = 4*pi*K*(alpha) dla roznych r_max ==")
    print(f"  {'alpha':>7}", end="")
    for r in RMAXES:
        print(f"  {'E*(r='+str(int(r))+')':>14}", end="")
    print(f"  {'E*(inf)':>14}")
    print("  " + "-"*74)

    E_by_rmax = {r: [] for r in RMAXES}   # r -> list of (alpha, E*)
    E_inf_list = []  # list of (alpha, E*(inf))

    for alpha in ALPHAS:
        print(f"  {alpha:.4f}", end="")
        for r in RMAXES:
            Es = results[alpha].get(r,(np.nan,)*3)[1]
            es_s = f'{Es:.6f}' if np.isfinite(Es) else '    brak   '
            print(f"  {es_s}", end="")
            if np.isfinite(Es):
                E_by_rmax[r].append((alpha, Es))
        # E*(inf) z K*(inf)
        Ki = K_inf.get(alpha, np.nan)
        if np.isfinite(Ki):
            Ei = 4*np.pi*Ki
            print(f"  {Ei:.6f}")
            E_inf_list.append((alpha, Ei))
        else:
            print(f"  {'brak':>14}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== MINIMUM E*(alpha) dla roznych r_max ==")

    argmin_vals = {}
    for r in RMAXES + ['inf']:
        if r == 'inf':
            pts = E_inf_list
        else:
            pts = E_by_rmax[r]
        if len(pts) < 3:
            print(f"  r_max={r}: za malo punktow")
            continue
        alphas_pts = np.array([x[0] for x in pts])
        E_pts = np.array([x[1] for x in pts])
        # Znajdz minimum
        idx_min = np.argmin(E_pts)
        alpha_min_raw = alphas_pts[idx_min]
        E_min = E_pts[idx_min]
        # Fit paraboliczny do 3 najblizszych punktow wokol minimum
        i0 = max(0, idx_min-1); i1 = min(len(E_pts), idx_min+2)
        seg_a = alphas_pts[i0:i1]; seg_E = E_pts[i0:i1]
        alpha_min_fit = alpha_min_raw
        if len(seg_a) >= 3:
            c = np.polyfit(seg_a, seg_E, 2)
            if c[0] > 0:  # concave up
                alpha_min_fit = -c[1]/(2*c[0])
        argmin_vals[r] = (alpha_min_raw, alpha_min_fit, E_min)
        delta = alpha_min_fit - ALPHA_K
        print(f"  r_max={r}: argmin_raw={alpha_min_raw:.4f}, "
              f"argmin_fit={alpha_min_fit:.4f}, E_min={E_min:.6f}, "
              f"delta_from_alphaK={delta:+.4f}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== TESTY ==")

    # P1: argmin(r_max=40) in [8.55, 8.70]
    am40 = argmin_vals.get(40.0, (np.nan,)*3)
    p1_ok = np.isfinite(am40[1]) and 8.55 <= am40[1] <= 8.70
    record("P1: argmin E*(alpha, r_max=40) in [8.55, 8.70] (zgod. z p82-p83)",
           p1_ok, f"argmin={am40[1]:.4f}" if np.isfinite(am40[1]) else "brak")

    # P2: argmin(r_max=60) blizej alpha_K niz argmin(r_max=40)?
    am60 = argmin_vals.get(60.0, (np.nan,)*3)
    if np.isfinite(am40[1]) and np.isfinite(am60[1]):
        d40 = abs(am40[1] - ALPHA_K)
        d60 = abs(am60[1] - ALPHA_K)
        p2_ok = d60 < d40
        record("P2: |argmin(r=60) - alpha_K| < |argmin(r=40) - alpha_K|",
               p2_ok, f"d40={d40:.4f}, d60={d60:.4f}")
    else:
        record("P2: zbieznosc argmin ku alpha_K", False, "brak danych")

    # P3: argmin(r_max=inf) ~ alpha_K (offset < 0.05)
    am_inf = argmin_vals.get('inf', (np.nan,)*3)
    if np.isfinite(am_inf[1]):
        delta_inf = abs(am_inf[1] - ALPHA_K)
        p3_ok = delta_inf < 0.05
        record("P3: |argmin E*(alpha, inf) - alpha_K| < 0.05",
               p3_ok, f"argmin_inf={am_inf[1]:.4f}, delta={delta_inf:.4f}")
    else:
        record("P3: argmin(inf) dostepne", False, "brak extrapolacji")

    # P4: E* ma minimum (nie monotoniczne) przy r_max=40
    if len(E_by_rmax[40.0]) >= 3:
        E_pts40 = [x[1] for x in E_by_rmax[40.0]]
        idx_min40 = np.argmin(E_pts40)
        # Jest minimum jesli nie na krawedzi
        p4_ok = (0 < idx_min40 < len(E_pts40)-1)
        record("P4: E*(alpha, r_max=40) ma minimum (nie na krawedzi)",
               p4_ok, f"idx_min={idx_min40}/{len(E_pts40)-1}")
    else:
        record("P4: brak danych r_max=40", False, "")

    # P5: K*(alpha, inf) ma minimum (nie monotoniczne)
    K_inf_list = [(a, K_inf[a]) for a in ALPHAS if np.isfinite(K_inf.get(a, np.nan))]
    if len(K_inf_list) >= 3:
        Kvals_inf = [x[1] for x in K_inf_list]
        idx_minK = np.argmin(Kvals_inf)
        p5_ok = (0 < idx_minK < len(Kvals_inf)-1)
        record("P5: K*(alpha, inf) ma minimum (nie monotoniczne)",
               p5_ok, f"idx_min={idx_minK}/{len(Kvals_inf)-1}, "
               f"a_minK={K_inf_list[idx_minK][0]:.4f}")
    else:
        record("P5: za malo punktow K*(inf)", False, "")

    # ─────────────────────────────────────────────────────────────────────────
    # Wykres
    print("\n== WYKRES E*(alpha) dla roznych r_max ==")
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        colors = {40.0: 'b', 60.0: 'g', 80.0: 'r', 'inf': 'k'}

        for r in RMAXES + ['inf']:
            if r == 'inf':
                pts = E_inf_list
                lbl = 'E*(inf) extrap'
                ls = '--'
            else:
                pts = E_by_rmax[r]
                lbl = f'E*(r_max={int(r)})'
                ls = '-'
            if pts:
                ax1.plot([x[0] for x in pts], [x[1] for x in pts],
                         f'o{ls}', color=colors[r], label=lbl)
                # Minimum marker
                am = argmin_vals.get(r)
                if am and np.isfinite(am[1]):
                    ax1.axvline(am[1], color=colors[r], alpha=0.3, ls=':')

        ax1.axvline(ALPHA_K, color='m', ls='--', lw=2, label=f'alpha_K={ALPHA_K}')
        ax1.set_xlabel('alpha'); ax1.set_ylabel('E*')
        ax1.set_title('E*(alpha) dla roznych r_max')
        ax1.legend(fontsize=9); ax1.grid(True, alpha=0.3)

        # K*(alpha, r_max)
        for r in RMAXES:
            pts = [(a, results[a][r][0]) for a in ALPHAS
                   if np.isfinite(results[a].get(r,(np.nan,)*3)[0])]
            if pts:
                ax2.plot([x[0] for x in pts], [x[1] for x in pts],
                         'o-', color=colors[r], label=f'K*(r={int(r)})')
        # K*(inf)
        Kinf_pts = [(a, K_inf[a]) for a in ALPHAS if np.isfinite(K_inf.get(a,np.nan))]
        if Kinf_pts:
            ax2.plot([x[0] for x in Kinf_pts], [x[1] for x in Kinf_pts],
                     's--', color='k', label='K*(inf)')
        ax2.axvline(ALPHA_K, color='m', ls='--', lw=2, label=f'alpha_K={ALPHA_K}')
        ax2.set_xlabel('alpha'); ax2.set_ylabel('K*')
        ax2.set_title('K*(alpha) dla roznych r_max')
        ax2.legend(fontsize=9); ax2.grid(True, alpha=0.3)

        plt.suptitle(f'p96: E*(alpha) zbieznosc z r_max  (alpha_K={ALPHA_K})')
        plt.tight_layout()
        out_path = 'p96_Estar_rmax_convergence.png'
        plt.savefig(out_path, dpi=120, bbox_inches='tight')
        plt.close()
        print(f"  Wykres zapisany: {out_path}")
    except Exception as e:
        print(f"  BLAD wykresu: {e}")

    # ─────────────────────────────────────────────────────────────────────────
    print()
    print("=" * 68)
    print("PODSUMOWANIE  p96_Estar_rmax_convergence.py")
    print("=" * 68)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")

    # Interpretacja OP-3 Sciezka 1
    print("\n== INTERPRETACJA OP-3 Sciezka 1 ==")
    if 'inf' in argmin_vals and np.isfinite(argmin_vals['inf'][1]):
        d = argmin_vals['inf'][1] - ALPHA_K
        if abs(d) < 0.03:
            print("  -> argmin E*(alpha, inf) ~ alpha_K: OP-3 Sciezka 1 POTWIERDZONA")
            print("     alpha_K jest wyznaczone przez minimum energii solitonu!")
        elif abs(d) < 0.10:
            print(f"  -> argmin E*(inf) = {argmin_vals['inf'][1]:.4f}, delta={d:+.4f}")
            print("     Minimum blisko alpha_K ale nie dokładnie — dalsze badania potrzebne")
        else:
            print(f"  -> argmin E*(inf) = {argmin_vals['inf'][1]:.4f}, delta={d:+.4f}")
            print("     Zbyt duze odchylenie — argmin E* nie jest mechanizmem alpha_K")

    # Czy argmin przesuwa sie z r_max?
    ams = [(r, argmin_vals[r][1]) for r in RMAXES + ['inf']
           if r in argmin_vals and np.isfinite(argmin_vals[r][1])]
    if len(ams) >= 2:
        rs_ams = [x[0] if x[0] != 'inf' else 1000 for x in ams]
        as_ams = [x[1] for x in ams]
        shift = as_ams[-1] - as_ams[0]
        print(f"\n  Przesuniecie argmin: r_max=40 ({as_ams[0]:.4f}) -> "
              f"r_max=inf ({as_ams[-1]:.4f}), shift={shift:+.4f}")
        if abs(shift) < 0.02:
            print("  -> argmin STABILNE wzgledem r_max (nie jest artefaktem!)")
        else:
            print("  -> argmin NIESTABILNE wzgledem r_max — podobne do K_nan_R")

    print()


if __name__ == '__main__':
    run_main()
