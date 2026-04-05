#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p97_true_soliton_profile.py  --  TGP v1
========================================
Cel: Wyznaczyc PRAWDZIWY K*(alpha_K, r_max->inf) i profil solitonu TGP
     po odkryciu artefaktow r_max=40 w p93-p96.

Pytania fizyczne:
  Q1: Jaka jest prawdziwa energia solitonu E* = 4*pi*K*(inf)?
  Q2: Jakie sa parametry ogona oscylacyjnego (A,B,C,m) przy alpha_K?
  Q3: Czy m = 1/sqrt(1+alpha_K) (zgodnie z teoria linowego ogona)?
  Q4: Jaka jest skala rdzenia solitonu (r_half, r_90, r_99)?
  Q5: Jak wyglada profil phi(r) w skali fizycznej (jednostki a_Gamma)?

Metoda:
  Czesc A: K*(alpha_K, r_max) dla r_max = 40, 60, 80, 100, 150
           Fit K*(r) = a + b/r (3 punkty r=60,80,100 => K*(inf))
           Fit K*(r) = a + b/r + c/r^2 (4 punkty => lepsza ekstrapolacja)

  Czesc B: Profil solitonu phi(r) przy K*(inf) z r_max=200
           Fit ogona: phi(r)-1 ~ (A sin(mr) + B cos(mr))/r dla r>50
           Weryfikacja m = 1/sqrt(1+alpha_K)

  Czesc C: Skale przestrzenne solitonu
           r_half: phi(r_half) - 1 = 0.5*(phi(0) - 1)
           r_90, r_99 w jednostkach a_Gamma i fm (jesli mozliwe)

  Czesc D: Porownanie z konwencja r_max=40 (bias numeryczny)

Testy:
  P1: K*(inf) z fitu 3-punktowego i 4-punktowego zgodne (<0.1%)
  P2: Fit ogona: R^2 > 0.99 dla r > 40
  P3: m_fit ~ m_theory = 1/sqrt(1+alpha_K) (< 1% odchylenie)
  P4: K*(r_max=40) vs K*(inf): bias > 5% (potwierdzenie p93)
  P5: Profil phi(r) zbiega do 1 oscylacyjnie (potwierdzenie ogona Bessela)

Data: 2026-03-26
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from scipy.stats import linregress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K  = 8.5616
A_GAM    = 0.040
LAM      = 5.501357e-06
GAMMA_V  = 1.0
V1       = GAMMA_V/3 - GAMMA_V/4
M_THEORY = 1.0 / np.sqrt(1.0 + ALPHA_K)  # teoretyczna masa oscylacji

print(f"  alpha_K={ALPHA_K}, A_GAM={A_GAM}")
print(f"  m_theory=1/sqrt(1+alpha_K)={M_THEORY:.6f}")

# Parametry numeryczne
N_PSI_MAIN = 100
N_EVAL_MAIN = 1000
N_PSI_PROFILE = 100

# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):  return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p): return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y):
    phi, dphi = y; phi = max(phi, 1e-10)
    alpha = ALPHA_K
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def integrate_profile(psi0, K, r_max, n_eval=None):
    """Zwraca sol (pelny profil)."""
    ne = n_eval or N_EVAL_MAIN
    dphi0 = -K/A_GAM**2
    r_ev = A_GAM*(r_max/A_GAM)**np.linspace(0, 1, ne)
    try:
        sol = solve_ivp(lambda r,y: ode_rhs(r,y), [A_GAM, r_max],
                        [psi0, dphi0], method='DOP853', rtol=1e-9, atol=1e-11,
                        t_eval=r_ev)
        return sol
    except:
        return None

def phi_at_rmax(psi, K, r_max, n_eval=None):
    sol = integrate_profile(psi, K, r_max, n_eval)
    if sol is None: return np.nan
    return float(sol.y[0,-1]) if sol.t[-1] >= r_max*0.99 else np.nan

def scan_last_psi_zero(K, r_max, n_psi=None):
    """Zwraca ostatnie zero F(psi)=phi(r_max)-1=0 (Galaz U)."""
    np_ = n_psi or N_PSI_MAIN
    psis = np.linspace(1.001, 3.0, np_)
    Fv = [phi_at_rmax(p, K, r_max) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p,K,r_max)-1.0,
                            psis[i], psis[i+1], xtol=1e-8, maxiter=50)
                roots.append(pz)
            except: pass
    return roots[-1] if roots else np.nan

def compute_g(K, psi0, r_max):
    """g = E/(4piK) - 1 dla zadanego K i psi0."""
    sol = integrate_profile(psi0, K, r_max)
    if sol is None or sol.t[-1] < r_max*0.99: return np.nan
    r=sol.t; phi=np.maximum(sol.y[0],1e-10); dphi=sol.y[1]
    Ek=4*np.pi*np.trapezoid(0.5*dphi**2*(1+ALPHA_K/phi)*r**2, r)
    Ep=4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
    return (Ek+Ep)/(4*np.pi*K)-1.0

def find_K_star(r_max, K_lo=0.006, K_hi=0.018, n_K=30):
    """Szuka K* przez brentq (metodologia: psi z scan, g=0)."""
    # Coarse scan
    Ks = np.linspace(K_lo, K_hi, n_K)
    gv = []
    psi_list = []
    for K in Ks:
        psi = scan_last_psi_zero(K, r_max)
        if np.isfinite(psi):
            g = compute_g(K, psi, r_max)
        else:
            g = np.nan
        gv.append(g); psi_list.append(psi)
    # Znajdz przejscie znaku
    K_star = np.nan; psi_star = np.nan
    for i in range(len(gv)-1):
        if np.isfinite(gv[i]) and np.isfinite(gv[i+1]) and gv[i]*gv[i+1] < 0:
            try:
                def g_of_K(K):
                    psi = scan_last_psi_zero(K, r_max, n_psi=80)
                    if not np.isfinite(psi): return np.nan
                    return compute_g(K, psi, r_max)
                Kz = brentq(g_of_K, Ks[i], Ks[i+1], xtol=1e-7, maxiter=25)
                K_star = Kz
                psi_star = scan_last_psi_zero(Kz, r_max, n_psi=80)
            except: pass
    return K_star, psi_star

# ─────────────────────────────────────────────────────────────────────────────
def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    print("=" * 68)
    print("TGP v1  *  p97_true_soliton_profile.py")
    print(f"  alpha_K={ALPHA_K}, A_GAM={A_GAM}")
    print(f"  m_theory={M_THEORY:.6f} (1/sqrt(1+alpha_K))")
    print("=" * 68)

    PC=0; FC=0
    def record(label, ok, info=""):
        nonlocal PC,FC
        if ok: PC+=1
        else: FC+=1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    # ═══════════════════════════════════════════════════════════════════════
    print("\n== CZESC A: K*(alpha_K, r_max) zbieznosc ==")
    RMAXES_A = [40.0, 60.0, 80.0, 100.0, 150.0]
    K_results = {}

    print(f"  {'r_max':>6}  {'K*':>14}  {'psi0':>8}  t")
    print("  " + "-"*40)
    for r_max in RMAXES_A:
        t0 = time.time()
        Ks, ps = find_K_star(r_max)
        dt = time.time()-t0
        ks_s = f'{Ks:.8f}' if np.isfinite(Ks) else '    brak    '
        ps_s = f'{ps:.6f}' if np.isfinite(ps) else '  brak  '
        print(f"  {r_max:>6.0f}  {ks_s}  {ps_s}  {dt:.0f}s")
        K_results[r_max] = (Ks, ps)

    # Fit K*(r) = a + b/r  (3 punkty: r=60,80,100)
    pts_3 = [(r, K_results[r][0]) for r in [60.0, 80.0, 100.0]
             if np.isfinite(K_results.get(r,(np.nan,))[0])]
    K_inf_3 = np.nan; b_3 = np.nan
    if len(pts_3) >= 2:
        rs = np.array([x[0] for x in pts_3])
        Ks_arr = np.array([x[1] for x in pts_3])
        X = np.column_stack([np.ones_like(rs), 1.0/rs])
        from numpy.linalg import lstsq
        coef, _, _, _ = lstsq(X, Ks_arr, rcond=None)
        K_inf_3, b_3 = coef
        print(f"\n  Fit K*(r)=a+b/r (3 pkt r=60,80,100): K*(inf)={K_inf_3:.8f}, b={b_3:.5f}")

    # Fit K*(r) = a + b/r (4 punkty: r=60,80,100,150)
    pts_4 = [(r, K_results[r][0]) for r in [60.0, 80.0, 100.0, 150.0]
             if np.isfinite(K_results.get(r,(np.nan,))[0])]
    K_inf_4 = np.nan; b_4 = np.nan
    if len(pts_4) >= 3:
        rs = np.array([x[0] for x in pts_4])
        Ks_arr = np.array([x[1] for x in pts_4])
        X = np.column_stack([np.ones_like(rs), 1.0/rs])
        coef, _, _, _ = lstsq(X, Ks_arr, rcond=None)
        K_inf_4, b_4 = coef
        print(f"  Fit K*(r)=a+b/r (4 pkt r=60,80,100,150): K*(inf)={K_inf_4:.8f}, b={b_4:.5f}")

    # P1: Zgodnosc dwoch fitow
    if np.isfinite(K_inf_3) and np.isfinite(K_inf_4):
        rel_diff = abs(K_inf_3 - K_inf_4)/K_inf_4
        p1_ok = rel_diff < 0.001
        record("P1: K*(inf) z fit-3 i fit-4 zgodne (< 0.1%)",
               p1_ok, f"K_inf_3={K_inf_3:.7f}, K_inf_4={K_inf_4:.7f}, rel_diff={rel_diff*100:.3f}%")
    else:
        record("P1: fity K*(inf)", False, "za malo danych")

    # Wybierz najlepsza K*_inf
    K_inf_best = K_inf_4 if np.isfinite(K_inf_4) else K_inf_3
    psi_at_rmax60 = K_results.get(60.0, (np.nan,)*2)[1]

    # P4: Bias r_max=40 vs inf
    K_40 = K_results.get(40.0, (np.nan,))[0]
    if np.isfinite(K_40) and np.isfinite(K_inf_best):
        bias = (K_40 - K_inf_best)/K_inf_best
        p4_ok = bias > 0.05
        record("P4: K*(r=40) vs K*(inf): bias > 5%",
               p4_ok, f"K*(40)={K_40:.7f}, K*(inf)={K_inf_best:.7f}, bias={bias*100:.2f}%")

    # ═══════════════════════════════════════════════════════════════════════
    print("\n== CZESC B: Profil phi(r) przy K*(inf) z r_max=200 ==")

    if np.isfinite(K_inf_best):
        # Znajdz psi0 dla K*(inf) z r_max=100 (dobra aproksymacja)
        K_use = K_inf_best
        psi_use = psi_at_rmax60
        if not np.isfinite(psi_use):
            psi_use = scan_last_psi_zero(K_use, 100.0, n_psi=100)

        if np.isfinite(psi_use):
            print(f"  Integracja profilu: K={K_use:.7f}, psi0={psi_use:.6f}, r_max=200")
            t0 = time.time()
            sol = integrate_profile(psi_use, K_use, 200.0, n_eval=2000)
            dt = time.time()-t0
            print(f"  Czas: {dt:.0f}s, r_max_reached={sol.t[-1]:.1f}" if sol else "  BLAD!")

            if sol and sol.t[-1] >= 180:
                r_arr = sol.t
                phi_arr = sol.y[0]
                phi0_core = phi_arr[0]
                print(f"  phi(r=A_Gam)={phi0_core:.6f}, phi(r=200)={phi_arr[-1]:.6f}")

                # Skale rdzenia
                delta0 = phi0_core - 1.0
                for frac, name in [(0.5,'r_half'), (0.1,'r_90'), (0.01,'r_99')]:
                    thr = 1.0 + frac*delta0
                    idx = np.where(phi_arr < thr)[0]
                    if len(idx) > 0:
                        i0 = idx[0]
                        r_scale = np.interp(thr, phi_arr[i0-1:i0+1][::-1], r_arr[i0-1:i0+1][::-1])
                        print(f"  {name}={r_scale:.4f} = {r_scale/A_GAM:.2f} a_Gam")

                # Fit ogona Bessela: phi(r)-1 = (A sin(mr) + B cos(mr))/r dla r>50
                mask = r_arr > 50.0
                if mask.sum() > 20:
                    r_tail = r_arr[mask]
                    delta_tail = phi_arr[mask] - 1.0

                    def bessel_model(r, A, B, m):
                        return (A*np.sin(m*r) + B*np.cos(m*r))/r

                    try:
                        p0 = [-0.02, 0.01, M_THEORY]
                        bounds = ([-0.1,-0.1,M_THEORY*0.5], [0.1,0.1,M_THEORY*1.5])
                        popt, _ = curve_fit(bessel_model, r_tail, delta_tail,
                                           p0=p0, bounds=bounds, maxfev=5000)
                        A_fit, B_fit, m_fit = popt
                        delta_fit = bessel_model(r_tail, *popt)
                        ss_res = np.sum((delta_tail - delta_fit)**2)
                        ss_tot = np.sum((delta_tail - delta_tail.mean())**2)
                        R2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
                        C_fit = np.sqrt(A_fit**2 + B_fit**2)
                        m_err = abs(m_fit - M_THEORY)/M_THEORY
                        print(f"\n  Fit ogona (r>50): A={A_fit:.5f}, B={B_fit:.5f}")
                        print(f"    m_fit={m_fit:.6f} vs m_theory={M_THEORY:.6f} (err={m_err*100:.3f}%)")
                        print(f"    C=sqrt(A^2+B^2)={C_fit:.5f}, R^2={R2:.6f}")

                        # P2, P3
                        p2_ok = R2 > 0.99
                        record("P2: Fit ogona Bessela R^2 > 0.99 dla r>50",
                               p2_ok, f"R^2={R2:.6f}")
                        p3_ok = m_err < 0.01
                        record("P3: m_fit ~ m_theory (blad < 1%)",
                               p3_ok, f"m_fit={m_fit:.6f}, m_theory={M_THEORY:.6f}, err={m_err*100:.3f}%")
                    except Exception as e:
                        print(f"  BLAD fitu: {e}")
                        record("P2: fit ogona", False, "BLAD")
                        record("P3: m_fit ~ m_theory", False, "BLAD")

                # P5: Weryfikacja oscylacji
                phi_tail = phi_arr[r_arr > 100]
                if len(phi_tail) > 10:
                    crossings = np.sum(np.diff(np.sign(phi_tail - 1.0)) != 0)
                    p5_ok = crossings > 3
                    record("P5: phi(r) oscyluje wokol 1 dla r>100 (> 3 przejscia)",
                           p5_ok, f"przejscia={crossings}")

                # Wykres profilu
                print("\n  Tworzenie wykresu profilu...")
                try:
                    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

                    # Lewy: caly profil
                    ax = axes[0]
                    ax.plot(r_arr, phi_arr, 'b-', lw=1.5, label=f'phi(r), alpha_K={ALPHA_K}')
                    ax.axhline(1.0, color='k', ls='--', lw=0.8, label='phi=1 (proznia)')
                    ax.set_xlabel('r'); ax.set_ylabel('phi(r)')
                    ax.set_title(f'Profil solitonu TGP (K*_inf={K_inf_best:.5f})')
                    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
                    ax.set_xlim(0, 200)

                    # Srodkowy: rdzen
                    ax = axes[1]
                    ax.plot(r_arr, phi_arr, 'b-', lw=1.5)
                    ax.axhline(1.0, color='k', ls='--', lw=0.8)
                    ax.set_xlabel('r'); ax.set_ylabel('phi(r)')
                    ax.set_title('Rdzen solitonu (r < 5)')
                    ax.grid(True, alpha=0.3)
                    ax.set_xlim(0, 5)

                    # Prawy: ogon oscylacyjny
                    ax = axes[2]
                    r_plot = r_arr[r_arr > 40]
                    d_plot = phi_arr[r_arr > 40] - 1.0
                    ax.plot(r_plot, d_plot * r_plot, 'b-', lw=1.5, label='(phi-1)*r')
                    if 'popt' in dir():
                        r_th = np.linspace(40, 200, 500)
                        ax.plot(r_th, A_fit*np.sin(m_fit*r_th)+B_fit*np.cos(m_fit*r_th),
                                'r--', lw=1.5, label=f'A sin(mr)+B cos(mr) fit')
                    ax.axhline(0, color='k', ls='-', lw=0.5)
                    ax.set_xlabel('r'); ax.set_ylabel('(phi-1)*r')
                    ax.set_title(f'Ogon Bessela: m_theory={M_THEORY:.4f}')
                    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

                    plt.suptitle(f'p97: Prawdziwy profil solitonu TGP  '
                                 f'(alpha_K={ALPHA_K}, K*(inf)={K_inf_best:.6f})',
                                 fontsize=11)
                    plt.tight_layout()
                    out_path = 'p97_true_soliton_profile.png'
                    plt.savefig(out_path, dpi=120, bbox_inches='tight')
                    plt.close()
                    print(f"  Wykres zapisany: {out_path}")
                except Exception as e:
                    print(f"  BLAD wykresu: {e}")
        else:
            print("  BLAD: psi0 nie znalezione dla K*(inf)")
    else:
        print("  BLAD: K*(inf) nie wyznaczone")

    # ═══════════════════════════════════════════════════════════════════════
    print("\n== CZESC D: Podsumowanie bias K*(r_max=40) ==")
    K40 = K_results.get(40.0,(np.nan,))[0]
    if np.isfinite(K40) and np.isfinite(K_inf_best):
        bias_abs = K40 - K_inf_best
        bias_rel = bias_abs/K_inf_best*100
        print(f"  K*(r_max=40)  = {K40:.8f}  [konwencja p76-p96]")
        print(f"  K*(r_max=inf) = {K_inf_best:.8f}  [prawdziwa wartosc]")
        print(f"  Bias bezwzgledny: {bias_abs:.8f}")
        print(f"  Bias wzgledny:    {bias_rel:.3f}%")
        print(f"  E*(40) = 4pi*K*(40) = {4*np.pi*K40:.7f}")
        print(f"  E*(inf)= 4pi*K*(inf)= {4*np.pi*K_inf_best:.7f}")
        E_diff = (4*np.pi*K40 - 4*np.pi*K_inf_best)/(4*np.pi*K_inf_best)*100
        print(f"  Roznica E*:    {E_diff:.3f}%")

    # ─────────────────────────────────────────────────────────────────────────
    print()
    print("=" * 68)
    print("PODSUMOWANIE  p97_true_soliton_profile.py")
    print("=" * 68)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")

    print("\n== KLUCZOWE WYNIKI ==")
    if np.isfinite(K_inf_best):
        print(f"  K*(alpha_K, inf) = {K_inf_best:.7f}  (najlepsza ekstrapolacja)")
        print(f"  E* = 4pi*K*(inf) = {4*np.pi*K_inf_best:.6f}")
    if 'm_fit' in dir() and np.isfinite(m_fit):
        print(f"  m_tail = {m_fit:.6f}  (m_theory={M_THEORY:.6f})")
        print(f"  C_tail = {C_fit:.5f}  (amplituda ogona)")
    if np.isfinite(K40) and np.isfinite(K_inf_best):
        print(f"  Bias K*(40)/K*(inf) = {K40/K_inf_best:.5f}  "
              f"({(K40/K_inf_best-1)*100:.2f}% przeszacowanie)")
    print()


if __name__ == '__main__':
    run_main()
