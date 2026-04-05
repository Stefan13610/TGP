#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p92_two_soliton_profiles.py  --  TGP v1  *  Dwie klasy solitonow: profile i struktura
=======================================================================================
Motywacja (z p91):
  Istnieja DWA rodzaje fizycznych solitonow samozgodnych g(K)=0:
    - Galaz U (K*_B ~ 0.011): tradycyjny soliton, K niskie
    - Galaz L (K*_A ~ 0.014): nowa rodzina, K wyzsze
  Hipoteza uzytkownika: Galaz L to "obiekt rozlegly, grawitacyjnie aktywny,
  lecz bez mozliwosci formowania malej skali (ciemna materia?)".

Cel:
  Dla alpha=8.850 (obie galeze istnieja):
    A) Znajdz K*_B (Galaz U) i K*_A (Galaz L) precyzyjnie
    B) Oblicz profil phi(r) dla obu solitonow
    C) Oblicz gęstosc energii epsilon(r) = 0.5*(dphi)^2*(1+alpha/phi) + (V(phi)-V1)
    D) Diagnostyka "kompaktowosci":
       - phi_max = phi(a_gamma) = psi_0
       - r_half: promien gdzie phi(r) = 0.5*(psi_0 + 1.0)  [polowa amplitudy]
       - r_99: promien gdzie kumulatywna energia E(<r) = 0.99 * E_total
       - gradient wall: |dphi/dr| max = |dphi0| = K/a_gam^2
       - energia kinetyczna / potencjalna ratio
    E) Porownanie z: solitonem w alpha=0 (bez sprzezenia TGP)

Testy:
  P1: Oba solitony znalezione (K*_B i K*_A dla alpha=8.850)
  P2: Profile phi(r) monotonicznie malejace od psi_0 do 1
  P3: r_half(Galaz L) vs r_half(Galaz U): ktory bardziej "rozlegly"?
  P4: r_99 energia: porownanie stopnia kompaktowosci
  P5: Interpretacja fizyczna: czy Galaz L odpowiada "obiektu bez malej skali"

Workers: 1 (profil: sekwencyjnie).
Data: 2026-03-25
"""

import sys, io, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
ALPHA    = 8.840      # Wspolna wartosc alpha dla porownania (p91: obie galeze potwierdzone)
A_GAM    = 0.040
LAM      = 5.501357e-06
GAMMA_V  = 1.0
V1       = GAMMA_V/3 - GAMMA_V/4
R_MAX    = 60.0       # Dluzszy zasiag niz zwykle (chcemy zobaczyc "ogon")
N_EVAL   = 3000       # Gestszy skan r dla lepszego profilu

def V_mod(p):
    return GAMMA_V/3*p**3 - GAMMA_V/4*p**4 + LAM/6*(p-1)**6
def dV_mod(p):
    return GAMMA_V*p**2 - GAMMA_V*p**3 + LAM*(p-1)**5
def ode_rhs(r, y, alpha):
    phi, dphi = y; phi = max(phi, 1e-10)
    kfac = 1.0 + alpha/phi
    return [dphi, dV_mod(phi)/kfac + alpha*dphi**2/(2*phi**2*kfac) - 2/r*dphi]

def phi_at_rmax(psi, K, a_gam, alpha, r_max=R_MAX, n_eval=N_EVAL):
    dphi0 = -K/a_gam**2
    def ev(r, y): return y[0] - 1e-5
    ev.terminal = True; ev.direction = -1
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs(r, y, alpha), [a_gam, r_max],
                        [psi, dphi0], method='DOP853', rtol=1e-10, atol=1e-12,
                        t_eval=r_ev, events=[ev])
        return float(sol.y[0, -1]) if sol.t[-1] >= r_max*0.99 else np.nan
    except: return np.nan

def energy_at_psi(psi_z, K, a_gam, alpha, r_max=R_MAX, n_eval=N_EVAL):
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(lambda r, y: ode_rhs(r, y, alpha), [a_gam, r_max],
                        [psi_z, dphi0], method='DOP853', rtol=1e-10, atol=1e-12,
                        t_eval=r_ev)
        r = sol.t; phi = np.maximum(sol.y[0], 1e-10); dphi = sol.y[1]
        Ek = 4*np.pi*np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi*np.trapezoid((V_mod(phi)-V1)*r**2, r)
        return Ek+Ep, r, phi, dphi
    except: return np.nan, None, None, None

def full_profile(psi_z, K, a_gam, alpha, r_max=R_MAX, n_eval=N_EVAL):
    """Zwraca pelny profil (r, phi, dphi, eps_kin, eps_pot, eps_total, E_cumul)."""
    dphi0 = -K/a_gam**2
    r_ev = a_gam*(r_max/a_gam)**np.linspace(0, 1, n_eval)
    sol = solve_ivp(lambda r, y: ode_rhs(r, y, alpha), [a_gam, r_max],
                    [psi_z, dphi0], method='DOP853', rtol=1e-10, atol=1e-12,
                    t_eval=r_ev)
    r    = sol.t
    phi  = np.maximum(sol.y[0], 1e-10)
    dphi = sol.y[1]

    eps_kin  = 0.5*dphi**2*(1.0 + alpha/phi)
    eps_pot  = V_mod(phi) - V1
    eps_tot  = eps_kin + eps_pot

    # Skumulowana energia E(<r)
    dE = 4*np.pi*eps_tot*r**2
    E_cumul = np.zeros(len(r))
    for i in range(1, len(r)):
        E_cumul[i] = E_cumul[i-1] + 0.5*(dE[i]+dE[i-1])*(r[i]-r[i-1])

    return r, phi, dphi, eps_kin, eps_pot, eps_tot, E_cumul


def find_psi_branch(K, a_gam, alpha, branch='U', n_psi=150):
    """Znajdz psi_0 dla galezie U (druhie zero) lub L (pierwsze zero)."""
    psis = np.linspace(1.001, 3.0, n_psi)
    Fv = [phi_at_rmax(p, K, a_gam, alpha) - 1.0 for p in psis]
    roots = []
    for i in range(len(Fv)-1):
        if np.isfinite(Fv[i]) and np.isfinite(Fv[i+1]) and Fv[i]*Fv[i+1] < 0:
            try:
                pz = brentq(lambda p: phi_at_rmax(p, K, a_gam, alpha) - 1.0,
                            psis[i], psis[i+1], xtol=1e-9, maxiter=80)
                roots.append(pz)
            except: pass
    if branch == 'U' and len(roots) >= 2:
        return roots[1]
    elif branch == 'L' and len(roots) >= 1:
        return roots[0]
    return np.nan


def find_K_star(branch, K_lo, K_hi, a_gam, alpha, n_psi=150, tol=1e-5):
    """Bisection na K: szuka zera g_branch(K)."""
    def g_branch(K):
        psi_z = find_psi_branch(K, a_gam, alpha, branch, n_psi)
        if not np.isfinite(psi_z): return np.nan
        E, _, _, _ = energy_at_psi(psi_z, K, a_gam, alpha)
        if not np.isfinite(E): return np.nan
        return E / (4*np.pi*K) - 1.0

    glo = g_branch(K_lo); ghi = g_branch(K_hi)
    if not (np.isfinite(glo) and np.isfinite(ghi) and glo*ghi < 0):
        return np.nan, np.nan, np.nan
    try:
        Kz = brentq(g_branch, K_lo, K_hi, xtol=tol, rtol=tol, maxiter=40)
        psi_z = find_psi_branch(Kz, a_gam, alpha, branch, n_psi)
        g_val = g_branch(Kz)
        return Kz, psi_z, g_val
    except: return np.nan, np.nan, np.nan


def profile_diagnostics(r, phi, dphi, eps_kin, eps_tot, E_cumul, K, psi_0, label):
    """Oblicz i wyswietl diagnostyki kompaktowosci solitonu."""
    E_total = E_cumul[-1]

    # phi_max = psi_0 (warunek poczatkowy)
    phi_max  = psi_0
    phi_half = 0.5*(phi_max + 1.0)

    # r_half: gdzie phi = phi_half
    r_half = np.nan
    for i in range(len(phi)-1):
        if phi[i] >= phi_half >= phi[i+1]:
            # Interpolacja liniowa
            t = (phi[i] - phi_half) / (phi[i] - phi[i+1])
            r_half = r[i] + t*(r[i+1]-r[i])
            break

    # r_90, r_99: energia kumulatywna
    r_90 = r[np.searchsorted(E_cumul/E_total, 0.90)] if E_total > 0 else np.nan
    r_99 = r[np.searchsorted(E_cumul/E_total, 0.99)] if E_total > 0 else np.nan

    # gradient maksymalny (na scianie babelka)
    grad_max = abs(dphi[0])  # |-K/a_gam^2| przy r=a_gam

    # Ratio kinetic/potential
    Ek = np.trapezoid(eps_kin * r**2, r)
    Ep_pos = np.trapezoid(np.maximum(eps_tot - eps_kin, 0) * r**2, r)
    ratio_kp = Ek / (Ep_pos + 1e-30)

    # Kompaktosc: czy soliton "maly" wzgledem r_max?
    compactness = r_90 / R_MAX

    print(f"\n  --- {label} ---")
    print(f"  K*           = {K:.7f}")
    print(f"  psi_0        = {psi_0:.6f}")
    print(f"  E_total      = {E_total:.6f}")
    print(f"  g = E/(4piK)-1 = {E_total/(4*np.pi*K)-1:.2e}")
    print(f"  phi_max      = {phi_max:.6f}")
    print(f"  grad_wall    = |dphi(a_gam)| = {grad_max:.4f}  [= K/a_gam^2 = {K/A_GAM**2:.4f}]")
    print(f"  r_half (phi) = {r_half:.4f}  [w jedn. a_gam: {r_half/A_GAM:.1f}]")
    print(f"  r_90  (E)    = {r_90:.4f}  [w jedn. a_gam: {r_90/A_GAM:.1f}]")
    print(f"  r_99  (E)    = {r_99:.4f}  [w jedn. a_gam: {r_99/A_GAM:.1f}]")
    print(f"  E_kin/E_pot  = {ratio_kp:.4f}")
    print(f"  Kompaktosc   = r_90/R_max = {compactness:.4f}")

    # Profil phi(r) w tekstowym wykresie (20 punktow)
    r_plot = np.linspace(A_GAM, min(8.0, R_MAX), 40)
    phi_interp = np.interp(r_plot, r, phi)
    phi_norm = (phi_interp - 1.0) / (phi_max - 1.0)

    print(f"\n  Profil phi(r) [znorm. 0..1], r in [{A_GAM:.3f}, 8.0]:")
    print(f"  r/a_gam  |{'phi_norm':^30}| phi")
    for i in range(0, 40, 4):
        bar_len = int(phi_norm[i] * 28)
        bar_len = max(0, min(28, bar_len))
        bar = '#' * bar_len + '.' * (28-bar_len)
        print(f"  {r_plot[i]/A_GAM:6.1f}   |{bar}| {phi_interp[i]:.5f}")

    return {
        'K': K, 'psi_0': psi_0, 'E': E_total,
        'r_half': r_half, 'r_90': r_90, 'r_99': r_99,
        'grad_wall': grad_max, 'ratio_kp': ratio_kp
    }


def run_main():
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                      errors='replace')

    print("=" * 72)
    print("TGP v1  *  p92_two_soliton_profiles.py")
    print("  Dwie klasy solitonow: struktura i kompaktosc")
    print(f"  alpha = {ALPHA},  a_gam = {A_GAM},  R_max = {R_MAX}")
    print("=" * 72)

    PC = 0; FC = 0

    def record(label, ok, info=""):
        nonlocal PC, FC
        if ok: PC += 1
        else: FC += 1
        print(f"  [{'v' if ok else 'x'}] {label}: {info}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== KROK 1: Znajdz K*_B (Galaz U) ==")
    t0 = time.time()
    # K*_B ~ 0.01082 przy alpha=8.840 (z p91)
    K_B, psi_B, g_B = find_K_star('U', 0.0100, 0.0118, A_GAM, ALPHA)
    print(f"  K*_B  = {K_B:.7f}   psi_0 = {psi_B:.6f}   g = {g_B:.2e}   t={time.time()-t0:.1f}s")

    print("\n== KROK 2: Znajdz K*_A (Galaz L) ==")
    t0 = time.time()
    # K*_A ~ 0.01440 przy alpha=8.840 (z p91); K_nan_R=0.01539 -> K_hi < K_nan_R
    K_A, psi_A, g_A = find_K_star('L', 0.0125, 0.0152, A_GAM, ALPHA)
    print(f"  K*_A  = {K_A:.7f}   psi_0 = {psi_A:.6f}   g = {g_A:.2e}   t={time.time()-t0:.1f}s")

    found_B = np.isfinite(K_B) and abs(g_B) < 0.01
    found_A = np.isfinite(K_A) and abs(g_A) < 0.01
    record("P1: Oba solitony znalezione", found_B and found_A,
           f"K*_B={K_B:.6f} (g={g_B:.2e}), K*_A={K_A:.6f} (g={g_A:.2e})")

    if not found_B and not found_A:
        print("\n  BLAD: Brak obu solitonow. Sprawdz zakresy K.")
        return

    # ─────────────────────────────────────────────────────────────────────────
    print("\n== KROK 3: Pelne profile ==")
    diag = {}

    if found_B:
        t0 = time.time()
        r_B, phi_B, dphi_B, ek_B, ep_B, et_B, ec_B = full_profile(psi_B, K_B, A_GAM, ALPHA)
        print(f"  Profil Galaz U: t={time.time()-t0:.1f}s")
        ok_B = np.all(np.diff(phi_B) <= 0.001)  # monotonicznie (z tolerancja)
        record("P2a: Profil Galaz U monot. maleje", ok_B, f"phi: {phi_B[0]:.4f} -> {phi_B[-1]:.4f}")
        diag['B'] = profile_diagnostics(r_B, phi_B, dphi_B, ek_B, et_B, ec_B, K_B, psi_B,
                                         f"GALAZ U (K*_B={K_B:.6f}, 'tradycyjny soliton')")

    if found_A:
        t0 = time.time()
        r_A, phi_A, dphi_A, ek_A, ep_A, et_A, ec_A = full_profile(psi_A, K_A, A_GAM, ALPHA)
        print(f"\n  Profil Galaz L: t={time.time()-t0:.1f}s")
        ok_A = np.all(np.diff(phi_A) <= 0.001)
        record("P2b: Profil Galaz L monot. maleje", ok_A, f"phi: {phi_A[0]:.4f} -> {phi_A[-1]:.4f}")
        diag['A'] = profile_diagnostics(r_A, phi_A, dphi_A, ek_A, et_A, ec_A, K_A, psi_A,
                                         f"GALAZ L (K*_A={K_A:.6f}, 'nowy soliton')")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("POROWNANIE DWOCH KLAS SOLITONOW")
    print("=" * 72)

    if found_B and found_A:
        dB = diag['B']; dA = diag['A']

        print(f"\n  {'Wielkosc':24} {'Galaz U (K*_B)':>18} {'Galaz L (K*_A)':>18} {'Stosunek L/U':>14}")
        print(f"  {'-'*24}   {'-'*16}   {'-'*16}   {'-'*12}")

        def row(name, vB, vA, fmt='.5f'):
            ratio = vA/vB if vB != 0 else np.nan
            rs = f'{ratio:.3f}' if np.isfinite(ratio) else 'n/a'
            print(f"  {name:24}   {vB:{fmt}}       {vA:{fmt}}       {rs:>12}")

        row("K*",                  dB['K'],        dA['K'])
        row("psi_0 (amplituda)",   dB['psi_0'],    dA['psi_0'])
        row("E_total",             dB['E'],         dA['E'])
        row("|dphi| na scianie",   dB['grad_wall'], dA['grad_wall'])
        row("r_half / a_gam",      dB['r_half']/A_GAM, dA['r_half']/A_GAM, '.3f')
        row("r_90  / a_gam",       dB['r_90']/A_GAM,   dA['r_90']/A_GAM,   '.3f')
        row("r_99  / a_gam",       dB['r_99']/A_GAM,   dA['r_99']/A_GAM,   '.3f')
        row("E_kin/E_pot",         dB['ratio_kp'], dA['ratio_kp'])

        r_half_ratio = dA['r_half'] / dB['r_half'] if np.isfinite(dA['r_half']) and np.isfinite(dB['r_half']) else np.nan
        r_90_ratio   = dA['r_90']   / dB['r_90']   if np.isfinite(dA['r_90'])   and np.isfinite(dB['r_90'])   else np.nan
        r_99_ratio   = dA['r_99']   / dB['r_99']   if np.isfinite(dA['r_99'])   and np.isfinite(dB['r_99'])   else np.nan

        print(f"\n  Stosunek promieni (Galaz L / Galaz U):")
        rhs = f'{r_half_ratio:.3f}' if np.isfinite(r_half_ratio) else 'n/a'
        r90s = f'{r_90_ratio:.3f}' if np.isfinite(r_90_ratio) else 'n/a'
        r99s = f'{r_99_ratio:.3f}' if np.isfinite(r_99_ratio) else 'n/a'
        print(f"    r_half (L/U) = {rhs}")
        print(f"    r_90   (L/U) = {r90s}")
        print(f"    r_99   (L/U) = {r99s}")

        is_compact_B = np.isfinite(dB['r_90']) and dB['r_90'] < 5.0
        is_compact_A = np.isfinite(dA['r_90']) and dA['r_90'] < 5.0
        A_more_extended = np.isfinite(r_90_ratio) and r_90_ratio > 1.5

        record("P3: r_half(L) != r_half(U) — rozne rozmiary",
               np.isfinite(r_half_ratio) and abs(r_half_ratio - 1.0) > 0.1,
               f"r_half ratio = {rhs}")
        record("P4: Stopien kompaktowosci rozny (r_99)",
               A_more_extended or (is_compact_B != is_compact_A),
               f"r_90(U)={dB['r_90']:.3f}, r_90(L)={dA['r_90']:.3f}, ratio={r90s}")

        # Interpretacja
        print("\n" + "=" * 72)
        print("INTERPRETACJA FIZYCZNA")
        print("=" * 72)

        K_ratio = dA['K'] / dB['K']
        E_ratio = dA['E'] / dB['E']
        psi_ratio = (dA['psi_0']-1) / (dB['psi_0']-1)

        print(f"""
  Galaz U (K*_B = {dB['K']:.5f}):
    - Nizsze K -> lagodniejszy gradient na scianie babelka
    - Amplituda psi_0 - 1 = {dB['psi_0']-1:.5f}
    - Tradycyjny soliton (znany z p10, p76, p80...)
    - Kandydat: KOMPAKTOWY soliton (badzo skondensowany rdzen)

  Galaz L (K*_A = {dA['K']:.5f}):
    - Wyzsze K (o {(K_ratio-1)*100:.1f}%) -> silniejszy gradient
    - Amplituda psi_0 - 1 = {dA['psi_0']-1:.5f}  (roznica {(psi_ratio-1)*100:.1f}%)
    - Nowa rodzina, znana dopiero z p91
    - Energia wyzsze o {(E_ratio-1)*100:.1f}%
""")

        if A_more_extended:
            print("  WNIOSEK: Galaz L jest ROZLEGLEJSZA (r_90 wikszy) niz Galaz U.")
            print("  Wspiera hipoteze: Galaz L = 'obiekt bez malej skali',")
            print("  grawitacyjnie aktywny (E wyzsze), ale bez kompaktowego rdzenia.")
            print("  Potencjalny kandydat na TGP-analog ciemnej materii lub")
            print("  objetosci kosmologicznej o rozleglejszej strukturze.")
            interp_ok = True
        elif not A_more_extended and is_compact_B == is_compact_A:
            print("  WNIOSEK: Oba solitony maja podobny stopien kompaktowosci.")
            print("  Galaz L NIE jest wyraznie bardziej rozlegla.")
            print("  Roznica miedzy klasami jest subtelniejsza (K i energia).")
            interp_ok = False
        else:
            print("  WNIOSEK: Galaz U bardziej rozlegla niz L (niespodziewane!).")
            print("  Wymaga dalszej analizy.")
            interp_ok = False

        record("P5: Interpretacja fizyczna wsparta danymi", interp_ok,
               f"r_90 ratio = {r90s}")

    # ─────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("PODSUMOWANIE  p92_two_soliton_profiles.py")
    print("=" * 72)
    print(f"\n  PASS: {PC}/{PC+FC}")
    print(f"  FAIL: {FC}/{PC+FC}")


if __name__ == '__main__':
    run_main()
