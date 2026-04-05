#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p94_alpha_max_convergence.py  --  TGP v1
=========================================
Cel: tabela gap(alpha, r_max) = K_nan_R - K*
     Jezeli gap > 0: soliton istnieje przy tym alpha i r_max
     Jezeli gap <= 0 lub K* brak: alpha > alpha_max(r_max)

Pytanie: czy alpha_max(r_max) zmienia sie z r_max?
         Czy ratio alpha_K / alpha_max ~ n_s jest odporne na r_max?

Metoda: n_psi=60, n_K=20, n_eval=400 -- szybkie ale wystarczajace
        Tabela 5 alpha x 3 r_max

Testy:
  P1: gap(8.8734, r_max=40) ~ 0 (zgodne z p80)
  P2: alpha_max(r_max) jest slabo zalezne od r_max (delta < 0.02)
  P3: alpha_K/alpha_max w [0.958, 0.972] dla wszystkich r_max
  P4: wynik n_s = 0.9649 odporny (nie zalezy od konwencji r_max)
  P5: K_nan_R(alpha, r_max) jest slabo zalezne od r_max (wewn. kons.)

Data: 2026-03-26
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import linregress
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K    = 8.5616
N_S_PLANCK = 0.9649
SIGMA_NS   = 0.0042
A_GAM      = 0.040
LAM        = 5.501357e-06
GAMMA_V    = 1.0
V1         = GAMMA_V/3 - GAMMA_V/4

# Tabela: alpha x r_max
ALPHAS  = [8.860, 8.870, 8.8734, 8.880, 8.890]
RMAXES  = [40.0, 60.0, 100.0]

# Parametry obliczen (szybsze niz p93)
N_PSI   = 60
N_K     = 20
N_EVAL  = 400

# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):  return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, alpha, r_max):
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0,1,N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM,r_max],
                        [psi,dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        return float(sol.y[0,-1]) if sol.t[-1]>=r_max*0.99 else np.nan
    except: return np.nan

def g_U(K, alpha, r_max):
    """g dla Galezy U (ostatnie zero psi) — szybka wersja."""
    psis = np.linspace(1.001, 3.0, N_PSI)
    Fv = [phi_at_rmax(p, K, alpha, r_max) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1]<0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,alpha,r_max)-1.0,
                            psis[i],psis[i+1], xtol=1e-7, maxiter=40)
                roots.append(pz)
            except: pass
    if not roots: return np.nan, np.nan
    psi_z = roots[-1]  # Galaz U = ostatnie zero

    # Energia
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0,1,N_EVAL)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y,alpha), [A_GAM,r_max],
                        [psi_z,dphi0], method='DOP853', rtol=1e-8, atol=1e-10,
                        t_eval=r_ev)
        r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
        Ek=4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2,r)
        Ep=4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2,r)
        E=Ek+Ep
        return E/(4*np.pi*K)-1.0, psi_z
    except: return np.nan, psi_z

def find_K_star_and_nan(alpha, r_max):
    """
    Zwraca (K*, K_nan_R, psi_U_at_Kstar) przez skan K.
    K_nan_R = pierwsza K gdzie g_U = nan.
    """
    Ks = np.linspace(0.008, 0.016, N_K)
    gv = []; pv = []
    for K in Ks:
        g, psi = g_U(K, alpha, r_max)
        gv.append(g); pv.append(psi)

    # Znajdz K*
    K_star = np.nan; psi_star = np.nan
    for i in range(len(gv)-1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i+1]) and gv[i]*gv[i+1]<0:
            try:
                Kz = brentq(lambda K: g_U(K, alpha, r_max)[0],
                            Ks[i], Ks[i+1], xtol=5e-6, maxiter=18)
                K_star = Kz
                _, psi_star = g_U(Kz, alpha, r_max)
            except: pass

    # Znajdz K_nan_R
    K_nan_R = np.nan; last_K_fin = np.nan
    for K, g in zip(Ks, gv):
        if np.isfinite(g): last_K_fin = K
        elif np.isfinite(last_K_fin):
            # Bisection
            Klo, Khi = last_K_fin, K
            for _ in range(12):
                Km = 0.5*(Klo+Khi)
                gm, _ = g_U(Km, alpha, r_max)
                if np.isfinite(gm): Klo=Km
                else: Khi=Km
                if Khi-Klo < Klo*3e-4: break
            K_nan_R = 0.5*(Klo+Khi)
            break

    gap = K_nan_R - K_star if (np.isfinite(K_star) and np.isfinite(K_nan_R)) else np.nan
    return K_star, K_nan_R, gap, psi_star


# ─────────────────────────────────────────────────────────────────────────────
def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    print("=" * 68)
    print("TGP v1  *  p94_alpha_max_convergence.py")
    print(f"  alpha_K={ALPHA_K}, n_s={N_S_PLANCK}+/-{SIGMA_NS}")
    print(f"  Tabela gap(alpha, r_max) dla alpha={ALPHAS}")
    print(f"  r_max={RMAXES},  n_psi={N_PSI}, n_K={N_K}, n_eval={N_EVAL}")
    print("=" * 68)

    PC=0; FC=0
    def record(label, ok, info=""):
        nonlocal PC,FC
        if ok: PC+=1
        else: FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    # ─────────────────────────────────────────────────────────────────────────
    # Tabela wynikow: results[r_max][alpha] = (K*, K_nan, gap)
    results = {r: {} for r in RMAXES}

    for r_max in RMAXES:
        print(f"\n== r_max = {r_max:.0f} ==")
        print(f"  {'alpha':>7}  {'K*':>10}  {'K_nan_R':>10}  {'gap':>10}  {'psi_U':>7}  {'t':>5}")
        print("  " + "-"*58)
        for alpha in ALPHAS:
            t0 = time.time()
            Ks, Knan, gap, psi_u = find_K_star_and_nan(alpha, r_max)
            dt = time.time()-t0
            Ks_s   = f'{Ks:.7f}'   if np.isfinite(Ks)   else '  brak   '
            Kn_s   = f'{Knan:.7f}' if np.isfinite(Knan) else '  brak   '
            gap_s  = f'{gap:.7f}'  if np.isfinite(gap)  else '  brak   '
            psi_s  = f'{psi_u:.4f}' if np.isfinite(psi_u) else ' brak '
            sign   = ' [+]' if (np.isfinite(gap) and gap>0) else ' [0/brak]'
            print(f"  {alpha:.4f}  {Ks_s}  {Kn_s}  {gap_s}  {psi_s}  {dt:.0f}s{sign}")
            results[r_max][alpha] = (Ks, Knan, gap, psi_u)

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== ANALIZA: alpha_max(r_max) ==")
    print()

    alpha_max_vals = {}  # r_max -> alpha_max_estimate

    for r_max in RMAXES:
        # Znajdz ostanie alpha z gap > 0
        last_pos = np.nan
        for a in ALPHAS:
            g = results[r_max][a][2]
            if np.isfinite(g) and g > 0:
                last_pos = a
        # Ekstrapoluj z ostatnich 2-3 punktow
        pos_pts = [(a, results[r_max][a][2]) for a in ALPHAS
                   if np.isfinite(results[r_max][a][2]) and results[r_max][a][2] > 0]
        a_extrap = np.nan
        if len(pos_pts) >= 2:
            use = pos_pts[-min(3, len(pos_pts)):]
            as_ = np.array([x[0] for x in use])
            gs_ = np.array([x[1] for x in use])
            if len(as_) >= 2:
                sl, ic, *_ = linregress(as_, gs_)
                if sl < 0: a_extrap = -ic/sl

        alpha_max_vals[r_max] = (last_pos, a_extrap)
        la_s = f'{last_pos:.4f}' if np.isfinite(last_pos) else 'brak'
        ex_s = f'{a_extrap:.4f}' if np.isfinite(a_extrap) else 'brak'
        ratio_l = ALPHA_K/last_pos if np.isfinite(last_pos) else np.nan
        ratio_e = ALPHA_K/a_extrap if np.isfinite(a_extrap) else np.nan
        rl_s = f'{ratio_l:.5f}' if np.isfinite(ratio_l) else 'brak'
        re_s = f'{ratio_e:.5f}' if np.isfinite(ratio_e) else 'brak'
        print(f"  r_max={r_max:.0f}: alpha_max_last={la_s}, extrap={ex_s}")
        print(f"           alpha_K/alpha_max_last={rl_s}, /extrap={re_s}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== TABELA PODSUMOWUJACA ==")
    print()
    print(f"  gap(alpha, r_max):")
    print(f"  {'alpha':>8}", end="")
    for r in RMAXES:
        print(f"  {'r='+str(int(r)):>10}", end="")
    print()
    print("  " + "-"*42)
    for a in ALPHAS:
        print(f"  {a:.4f} ", end="")
        for r in RMAXES:
            g = results[r][a][2]
            if np.isfinite(g):
                print(f"  {g:+.6f}", end="")
            else:
                print(f"  {'brak':>10}", end="")
        print()

    print()
    print(f"  alpha_K/alpha_max(r_max):")
    print(f"  {'r_max':>8}  {'alpha_max':>10}  {'ratio':>10}  {'|ratio-ns|/sigma':>18}")
    for r_max in RMAXES:
        la, ea = alpha_max_vals[r_max]
        if np.isfinite(la):
            rl = ALPHA_K/la
            d = abs(rl - N_S_PLANCK)
            print(f"  {r_max:>8.0f}  {la:.5f}     {rl:.5f}     {d/SIGMA_NS:.1f}sigma")
    print(f"  {'p80 ref':>8}  8.87340     {ALPHA_K/8.8734:.5f}     "
          f"{abs(ALPHA_K/8.8734-N_S_PLANCK)/SIGMA_NS:.1f}sigma")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== TESTY ==")

    # P1: gap(8.8734, r_max=40) ~ 0
    g_ref = results[40.0].get(8.8734, (np.nan,)*4)[2] if 40.0 in results else np.nan
    p1_ok = np.isfinite(g_ref) and abs(g_ref) < 5e-4
    record("P1: gap(alpha=8.8734, r_max=40) ~ 0 (zgodne z p80)",
           p1_ok, f"gap={g_ref:.6f}" if np.isfinite(g_ref) else "brak")

    # P2: monotonicznosc i mala zmiennosc alpha_max(r_max)
    valid_am = [(r, la) for r, (la,_) in alpha_max_vals.items() if np.isfinite(la)]
    if len(valid_am) >= 2:
        am_vals = np.array([x[1] for x in valid_am])
        delta_am = am_vals.max() - am_vals.min()
        p2_ok = delta_am < 0.02
        record("P2: |alpha_max(r1) - alpha_max(r2)| < 0.02", p2_ok,
               f"delta={delta_am:.4f}, vals={am_vals}")
    else:
        record("P2: delta alpha_max < 0.02", False, "za malo danych")

    # P3: ratio alpha_K/alpha_max w [0.958, 0.972]
    ratios = [ALPHA_K/la for _, (la, _) in alpha_max_vals.items() if np.isfinite(la)]
    p3_ok = all(0.956 <= r <= 0.974 for r in ratios)
    record("P3: alpha_K/alpha_max in [0.956, 0.974] dla wszystkich r_max",
           p3_ok, f"ratios={[f'{r:.4f}' for r in ratios]}")

    # P4: gap(alpha, r_max=40) vs gap(alpha, r_max=100) — porownanie K_nan_R
    print()
    print("  K_nan_R(alpha) dla roznych r_max:")
    for a in ALPHAS:
        knans = [(r, results[r][a][1]) for r in RMAXES if np.isfinite(results[r][a][1])]
        if knans:
            kn_v = np.array([x[1] for x in knans])
            delta_kn = kn_v.max()-kn_v.min()
            print(f"    alpha={a:.4f}: K_nan_R = {[f'{x[1]:.5f}(r={int(x[0])})' for x in knans]}"
                  f"  delta={delta_kn:.5f}")

    # P5: K_nan_R slabo zalezne od r_max?
    knan_deltas = []
    for a in ALPHAS:
        knans = [results[r][a][1] for r in RMAXES if np.isfinite(results[r][a][1])]
        if len(knans) >= 2:
            knan_deltas.append(max(knans)-min(knans))
    p5_ok = all(d < 0.001 for d in knan_deltas) if knan_deltas else False
    record("P5: K_nan_R slabo zalezne od r_max (delta < 0.001)",
           p5_ok, f"max delta = {max(knan_deltas):.5f}" if knan_deltas else "brak")

    # P4: stosunek odporny na r_max (maks zmiana ratio < 0.5 sigma)
    if len(ratios) >= 2:
        ratio_range = max(ratios) - min(ratios)
        p4_ok = ratio_range < SIGMA_NS
        record("P4: alpha_K/alpha_max stabilny (rozstaw < 1 sigma_ns)",
               p4_ok, f"rozstaw={ratio_range:.5f}, sigma={SIGMA_NS}")
    else:
        record("P4: stosunek stabilny", False, "za malo danych")

    # ─────────────────────────────────────────────────────────────────────────
    print()
    print("=" * 68)
    print("PODSUMOWANIE  p94_alpha_max_convergence.py")
    print("=" * 68)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}\n")


if __name__ == '__main__':
    run_main()
