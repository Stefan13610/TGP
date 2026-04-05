#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p76_alpha_K_bifurcation.py  --  TGP v1 * OP-3  (v3 - poprawione kryterium)
============================================================================
Cel: Numeryczne wyznaczenie alpha_K z warunkow substratu TGP.

Model: BUBBLE SOLITON z p10/p11 (kfac = 1 + alpha/phi).

Kryterium samospojnosci [POPRAWIONE v2]:
  g(K, a_Gam) = E[phi] / (4*pi*K) - 1 = 0
  gdzie E = 4*pi * integral [0.5*dphi^2*(1+alpha/phi) + V_mod(phi)-V1] r^2 dr

  --- NIE samo phi(r_max)=1 (to warunek brzegowy, zawsze trywialnie spelniony
      przy phi=1+epsilon i malym K).

Algorytm:
  1. compute_g(K, a_gam, alpha):
     - Skan psi_core -> znajdz phi(r_max)=1 (brentq)
     - Oblicz energie E -> zwroc g = E/(4*pi*K) - 1
  2. find_K_star(a_gam, alpha):
     - Skan K -> znajdz K gdzie g(K) = 0 (brentq)
  3. Bifurkacja a_c(alpha): min a_gam gdzie K*_1 istnieje

Testy PASS/FAIL:
  P1: K*_1 znaleziony przy alpha=8.5616, a_Gam=0.040
  P2: |K*_1_num - 0.010414| < 5%
  P3: K*_1(alpha) przechodzi przez 0.010414 przy alpha_c bliskim 8.5616
  P4: |alpha_c - 8.5616| < 25%  (interpolacja z Sekcji B)
  P5: K*_1(alpha) wykazuje monotoniczna tendencje

Data: 2026-03-24 (v3 - poprawione kryterium energetyczne)
"""

import sys
import io
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

# ─────────────────────────────────────────────────────────────────────────────
# Stale TGP
# ─────────────────────────────────────────────────────────────────────────────
ALPHA_K_REF = 8.5616
A_GAM_TGP   = 0.040049
K1_REF      = 0.010414
LAM         = 5.501357e-06
GAMMA_V     = 1.0

# ─────────────────────────────────────────────────────────────────────────────
# Potencjal V_mod
# ─────────────────────────────────────────────────────────────────────────────
def V_mod(p):
    return GAMMA_V/3 * p**3 - GAMMA_V/4 * p**4 + LAM/6 * (p - 1)**6

def dV_mod(p):
    return GAMMA_V * p**2 - GAMMA_V * p**3 + LAM * (p - 1)**5

V1 = GAMMA_V/3 - GAMMA_V/4   # = V_mod(1) = 1/12


def ode_rhs(r, y, alpha):
    """
    RHS dla ODE bubble solitonu TGP:
      phi'' = dV_mod(phi)/kfac + alpha*phi'^2/(2*phi^2*kfac) - (2/r)*phi'
      kfac = 1 + alpha/phi
    """
    phi, dphi = y
    phi = max(phi, 1e-10)
    kfac  = 1.0 + alpha / phi
    ddphi = (dV_mod(phi) / kfac
             + alpha * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]


def _phi_at_rmax(psi_core, K, a_gam, alpha, r_max, rtol=1e-7, atol=1e-9, n_eval=800):
    """
    Calkuje ODE od a_gam do r_max. Zwraca phi(r_max) lub nan.
    Uzywa natywnego ev_lo z terminal=True (nie przez lambda).
    """
    dphi0 = -K / a_gam**2

    def ev_lo(r, y):
        return y[0] - 1e-5
    ev_lo.terminal  = True
    ev_lo.direction = -1

    r_eval = a_gam * (r_max / a_gam) ** np.linspace(0, 1, n_eval)
    try:
        sol = solve_ivp(
            lambda r, y: ode_rhs(r, y, alpha),
            [a_gam, r_max], [psi_core, dphi0],
            method='DOP853', rtol=rtol, atol=atol,
            t_eval=r_eval, events=[ev_lo]
        )
        if len(sol.t) == 0 or sol.t[-1] < r_max * 0.99:
            return np.nan
        return float(sol.y[0, -1])
    except Exception:
        return np.nan


def _energy(psi_z, K, a_gam, alpha, r_max):
    """
    Oblicza energie E = Ek + Ep dla profilu phi(psi_z).
    Uzywa wysokiej dokladnosci (n_eval=4000).
    """
    dphi0 = -K / a_gam**2

    def ev_lo(r, y):
        return y[0] - 1e-5
    ev_lo.terminal  = True
    ev_lo.direction = -1

    r_ev = a_gam * (r_max / a_gam) ** np.linspace(0, 1, 4000)
    try:
        sol = solve_ivp(
            lambda r, y: ode_rhs(r, y, alpha),
            [a_gam, r_max], [psi_z, dphi0],
            method='DOP853', rtol=1e-9, atol=1e-11,
            t_eval=r_ev, events=[ev_lo]
        )
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        Ek = 4*np.pi * np.trapezoid(0.5 * dphi**2 * (1 + alpha/phi) * r**2, r)
        Ep = 4*np.pi * np.trapezoid((V_mod(phi) - V1) * r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan


def compute_g(K, a_gam, alpha, r_max=None, n_psi=120):
    """
    Oblicza g(K, a_gam, alpha) = E[phi] / (4*pi*K) - 1.

    Portowane 1:1 z p10.best_g() — znajdz WSZYSTKIE galezi phi(r_max)=1,
    oblicz g dla kazdej i zwroc galaz z NAJMNIEJSZYM |g|.
    (Zapobiega "skokowi" na inna galaz przy malym n_psi.)

    psi_max: jak w p10 -- max(3.0, 1.0 + 2.0*K/a_gam).
    """
    if r_max is None:
        M_eff = 1.0 / np.sqrt(1.0 + alpha)
        r_max = max(40.0, a_gam * 5, 6.0 / M_eff)

    psi_max  = max(3.0, 1.0 + 2.0 * K / a_gam)  # jak w p10!
    psi_vals = np.linspace(1.001, min(psi_max, 300.0), n_psi)

    F_vals = np.array([_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1.0
                       for p in psi_vals])

    # Zbierz WSZYSTKIE galezi (wszystkie zmiany znaku)
    branches = []
    for i in range(len(F_vals) - 1):
        Fi, Fj = F_vals[i], F_vals[i+1]
        if np.isfinite(Fi) and np.isfinite(Fj) and Fi * Fj < 0:
            try:
                psi_z = brentq(
                    lambda p: _phi_at_rmax(p, K, a_gam, alpha, r_max) - 1.0,
                    psi_vals[i], psi_vals[i+1],
                    xtol=1e-5, rtol=1e-5, maxiter=40
                )
                E = _energy(psi_z, K, a_gam, alpha, r_max)
                if np.isfinite(E):
                    g = E / (4*np.pi*K) - 1.0
                    branches.append((psi_z, g))
            except Exception:
                pass

    if not branches:
        return np.nan

    # Jak w p10: zwroc galaz z najmniejszym |g|
    best = min(branches, key=lambda x: abs(x[1]))
    return best[1]


def find_K_star(a_gam, alpha, K_min=3e-4, K_max=0.15, n_K=20, n_psi=40,
                verbose=False):
    """
    Znajduje K*_1 = pierwsze zero g(K) = 0 dla danych (a_gam, alpha).

    Algorytm:
      Skan K ∈ [K_min, K_max] (log-rownomierne).
      Dla kazdego K: compute_g(K, a_gam, alpha).
      Jesli g zmienia znak -> brentq -> K*_1.
      Zwraca K*_1 lub nan.
    """
    M_eff = 1.0 / np.sqrt(1.0 + alpha)
    r_max = max(80.0, a_gam * 5, 6.0 / M_eff)

    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)
    g_vals = []
    for K in K_scan:
        gv = compute_g(K, a_gam, alpha, r_max, n_psi=n_psi)
        g_vals.append(gv)
        if verbose:
            print(f"    K={K:.5f}  g={gv:.4f}" if np.isfinite(gv) else
                  f"    K={K:.5f}  g=nan")

    K_star = np.nan
    for i in range(len(g_vals) - 1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
            try:
                K_star = brentq(
                    lambda K: compute_g(K, a_gam, alpha, r_max, n_psi=n_psi),
                    K_scan[i], K_scan[i+1],
                    xtol=K_scan[i] * 0.02, rtol=0.01, maxiter=20
                )
                break
            except Exception:
                pass

    return K_star


# ─────────────────────────────────────────────────────────────────────────────
# PASS/FAIL
# ─────────────────────────────────────────────────────────────────────────────
PASS_COUNT = 0
FAIL_COUNT = 0


def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok:
        PASS_COUNT += 1
    else:
        FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")
    return ok


# ─────────────────────────────────────────────────────────────────────────────
# NAGLOWEK
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 68)
print("TGP v1  *  p76_alpha_K_bifurcation.py  --  OP-3  (v3)")
print("Model: bubble soliton ODE, kryterium: g = E/(4*pi*K) - 1 = 0")
print("=" * 68)
print(f"\n[PARAMETRY TGP]")
print(f"  alpha_K_ref = {ALPHA_K_REF}")
print(f"  A_Gam       = {A_GAM_TGP:.5f}  [H0/c0]")
print(f"  K*_1_ref    = {K1_REF:.6f}  (z p10/p11)")
print(f"  V_mod(1)    = {V1:.6f}  (1/12 = {1/12:.6f})")
print(f"  M_eff       = 1/sqrt(1+{ALPHA_K_REF}) = {1/np.sqrt(1+ALPHA_K_REF):.4f}")
print(f"\n[UWAGA] Kryterium v3: g(K) = E/(4pi*K) - 1 = 0 (jak w p10)")
print(f"  - Skan K ∈ [3e-4, 0.15] (n_K=20)")
print(f"  - compute_g: skan psi_core n_psi=40 + brentq + energia")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA A: Cross-check K*_1 przy alpha_K_ref i a_Gam = 0.040
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("SEKCJA A: Cross-check K*_1 przy alpha=8.5616, a_Gam=0.040 [z p10]")
print("=" * 68)

print(f"\n  Skan g(K) przy alpha={ALPHA_K_REF}, a_Gam=0.040 (n_K=20, n_psi=40)...")
print(f"  (Obliczanie energii dla kazdego K -- chwila...)\n")

K1_check = find_K_star(a_gam=0.040, alpha=ALPHA_K_REF,
                       K_min=3e-4, K_max=0.15, n_K=20, n_psi=40, verbose=True)

print(f"\n  Wynik:  K*_1 = {K1_check:.6f}  (ref: {K1_REF:.6f})")

ok_P1 = np.isfinite(K1_check)
record("P1: K*_1 znaleziony przy alpha=8.5616, a_Gam=0.040", ok_P1,
       f"K*_1 = {K1_check:.6f}" if ok_P1 else "brak rozwiazania g=0")

ok_P2 = ok_P1 and abs(K1_check - K1_REF) / K1_REF < 0.05
record("P2: |K*_1_num - 0.010414| < 5%", ok_P2,
       f"odch = {abs(K1_check-K1_REF)/K1_REF*100:.1f}%" if ok_P1
       else "K*_1 nie znaleziony")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA B: Skan alpha -- K*_1(alpha) przy a_Gam = 0.040
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("SEKCJA B: Skan alpha -- K*_1(alpha) przy a_Gam=0.040")
print("          (Szukamy alpha_c gdzie K*_1(alpha_c) = K1_REF = 0.010414)")
print("=" * 68)

alpha_scan = np.array([4.0, 5.0, 6.0, 7.0, 8.0, 8.5616, 9.0, 10.0, 12.0])
A_GAM_FIXED = 0.040

results_B = []
print(f"\n  a_Gam = {A_GAM_FIXED} (stale), n_K=15, n_psi=35\n")
print(f"  {'alpha':>8}  {'K*_1':>10}  {'K*_1/K1_ref':>12}  {'M_eff':>7}")
print("  " + "-" * 46)

for alpha in alpha_scan:
    print(f"  alpha={alpha:.4f} ...", flush=True)
    K_star = find_K_star(a_gam=A_GAM_FIXED, alpha=alpha,
                         K_min=3e-4, K_max=0.15, n_K=15, n_psi=35)
    M_eff = 1.0 / np.sqrt(1.0 + alpha)
    ratio = K_star / K1_REF if np.isfinite(K_star) else np.nan
    results_B.append((alpha, K_star, M_eff))
    if np.isfinite(K_star):
        print(f"  {alpha:>8.4f}  {K_star:>10.6f}  {ratio:>12.3f}  {M_eff:>7.4f}")
    else:
        print(f"  {alpha:>8.4f}  {'---':>10}  {'---':>12}  {M_eff:>7.4f}  [brak g=0]")

# Monotonicznosc K*_1(alpha)
valid_B = [(r[0], r[1]) for r in results_B if np.isfinite(r[1])]
ok_P5   = False
trend   = "brak danych"
if len(valid_B) >= 3:
    alphas_v = [r[0] for r in valid_B]
    K1s_v    = [r[1] for r in valid_B]
    diffs    = [K1s_v[i+1] - K1s_v[i] for i in range(len(K1s_v)-1)]
    n_pos    = sum(1 for d in diffs if d > 0)
    n_neg    = sum(1 for d in diffs if d < 0)
    ok_P5    = (n_pos >= len(diffs) * 0.7 or n_neg >= len(diffs) * 0.7)
    trend    = "rosnie" if n_pos > n_neg else "maleje"
    print(f"\n  K*_1 {trend} z alpha: {n_pos}/{len(diffs)} wzrostow, "
          f"{n_neg}/{len(diffs)} spadkow")

record("P5: K*_1(alpha) monotoniczna tendencja", ok_P5,
       trend if len(valid_B) >= 3 else "za malo punktow")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA B2: Szukanie alpha_c gdzie K*_1(alpha_c) = K1_REF = 0.010414
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "-" * 68)
print("SEKCJA B2: Interpolacja alpha_c (gdzie K*_1(alpha) = K1_REF)")
print("-" * 68)

alpha_c_interp = None
if len(valid_B) >= 2:
    alphas_v = np.array([r[0] for r in valid_B])
    K1s_v    = np.array([r[1] for r in valid_B])
    diff_K   = K1s_v - K1_REF

    crossings = np.where(np.diff(np.sign(diff_K)))[0]
    if len(crossings) > 0:
        i0, i1 = crossings[0], crossings[0] + 1
        a0, d0 = alphas_v[i0], diff_K[i0]
        a1, d1 = alphas_v[i1], diff_K[i1]
        alpha_c_interp = a0 + (-d0) * (a1 - a0) / (d1 - d0)
        print(f"\n  Interpolacja: alpha_c = {alpha_c_interp:.4f}")
        print(f"  (K*_1({a0:.2f})={K1s_v[i0]:.6f}, "
              f"K*_1({a1:.2f})={K1s_v[i1]:.6f}, ref={K1_REF:.6f})")
    else:
        print(f"\n  K*_1(alpha) nie przecina K1_REF w zakresie skanu.")
        if len(K1s_v) > 0:
            print(f"  K*_1 zakres: [{K1s_v.min():.6f}, {K1s_v.max():.6f}]  "
                  f"(K1_REF={K1_REF:.6f})")
            if K1s_v.min() > K1_REF:
                print("  -> K*_1 > K1_REF wszedzie: alpha_c ponizej zakresu skanu")
            elif K1s_v.max() < K1_REF:
                print("  -> K*_1 < K1_REF wszedzie: alpha_c powyzej zakresu skanu")


ok_P3 = alpha_c_interp is not None and np.isfinite(alpha_c_interp)
record("P3: alpha_c znaleziony (K*_1(alpha_c)=K1_REF)", ok_P3,
       f"alpha_c = {alpha_c_interp:.4f}" if ok_P3 else "nie znaleziony")

ok_P4 = ok_P3 and abs(alpha_c_interp - ALPHA_K_REF) / ALPHA_K_REF < 0.25
if ok_P3:
    record("P4: |alpha_c - 8.5616| < 25%", ok_P4,
           f"alpha_c={alpha_c_interp:.4f}, ref={ALPHA_K_REF}, "
           f"odch={abs(alpha_c_interp-ALPHA_K_REF)/ALPHA_K_REF*100:.1f}%")
else:
    record("P4: |alpha_c - 8.5616| < 25%", False, "alpha_c nie znaleziony")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA C: Bifurkacja a_c(alpha) — skan a_Gam przy alpha=8.5616
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("SEKCJA C: a_c(alpha) -- minimalne a_Gam dla solitonu (alpha=8.5616)")
print("          Czy K*_1 istnieje przy a_Gam < 0.040?")
print("=" * 68)

a_gam_test = np.array([0.010, 0.020, 0.030, 0.040, 0.060])
print(f"\n  Skan a_Gam ∈ [0.01, 0.06] przy alpha=8.5616:")
print(f"\n  {'a_Gam':>8}  {'K*_1':>10}  {'g=0?':>8}")
print("  " + "-" * 32)

results_C = []
for a_g in a_gam_test:
    print(f"  a_Gam={a_g:.3f} ...", flush=True)
    K_s = find_K_star(a_gam=a_g, alpha=ALPHA_K_REF,
                      K_min=1e-4, K_max=0.3, n_K=15, n_psi=35)
    results_C.append((a_g, K_s))
    if np.isfinite(K_s):
        print(f"  {a_g:>8.3f}  {K_s:>10.6f}  {'TAK':>8}")
    else:
        print(f"  {a_g:>8.3f}  {'---':>10}  {'NIE':>8}")

# Wnioski
a_c_arr = np.array([r[0] for r in results_C])
K_c_arr = np.array([r[1] for r in results_C])
exists_C = np.isfinite(K_c_arr)

if exists_C.any():
    idx_first = np.where(exists_C)[0][0]
    a_c_est   = a_c_arr[idx_first]
    print(f"\n  Minimalne a_Gam z K*_1: {a_c_est:.3f}")
    if a_c_est < A_GAM_TGP * 0.95:
        print(f"  -> Soliton istnieje dla a_Gam < {A_GAM_TGP:.3f}")
        print(f"     a_c({ALPHA_K_REF}) < A_Gam_TGP -- bifurkacja ponizej obserwacji")
    elif abs(a_c_est - A_GAM_TGP) < 0.005:
        print(f"  -> a_c ≈ A_Gam_TGP = {A_GAM_TGP:.4f} -- bifurkacja PRZY obserwacji!")
    else:
        print(f"  -> a_c > {A_GAM_TGP:.3f} -- soliton wymaga wiekszego a_Gam")
else:
    print(f"\n  Brak K*_1 dla zadnego a_Gam w skanowanym zakresie!")


# ─────────────────────────────────────────────────────────────────────────────
# SEKCJA D: Path 1 (Koide) — skrocony test K*_2, K*_3
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("SEKCJA D: Path 1 (Koide) — czy g=0 ma drugie zero ponad K*_1?")
print("=" * 68)

# Skan g(K) od K*_1 do K=3 przy alpha=8.5616, a_Gam=0.040
K_high_scan = np.logspace(
    np.log10(max(K1_check * 1.1 if np.isfinite(K1_check) else 0.015, 0.015)),
    np.log10(2.0), 12
)

print(f"\n  Skan g(K) ∈ [{K_high_scan[0]:.3f}, 2.0] przy alpha=8.5616, a_Gam=0.040")
print(f"  (szukamy drugiego zera g ponad K*_1)\n")
print(f"  {'K':>10}  {'g(K)':>10}  {'status':>10}")
print("  " + "-" * 36)

M_eff_K = 1.0 / np.sqrt(1.0 + ALPHA_K_REF)
r_max_D  = max(80.0, 0.040 * 5, 6.0 / M_eff_K)

g_D_vals = []
for K_d in K_high_scan:
    gd = compute_g(K_d, 0.040, ALPHA_K_REF, r_max_D, n_psi=35)
    g_D_vals.append(gd)
    tag = "OK" if np.isfinite(gd) else "nan"
    gstr = f"{gd:.4f}" if np.isfinite(gd) else "---"
    print(f"  {K_d:>10.5f}  {gstr:>10}  {tag:>10}")

# Czy g zmienia znak (sugerujac K*_2)?
K2_found = False
for i in range(len(g_D_vals)-1):
    gi, gj = g_D_vals[i], g_D_vals[i+1]
    if np.isfinite(gi) and np.isfinite(gj) and gi*gj < 0:
        K2_found = True
        print(f"\n  [!] Zmiana znaku g miedzy K={K_high_scan[i]:.4f} i "
              f"K={K_high_scan[i+1]:.4f} -- mozliwe K*_2!")
        break

if not K2_found:
    print(f"\n  Brak zmian znaku g -- K*_2 nie istnieje w biezacym V_mod.")
    print(f"  Path 1 (Koide) wymaga modyfikacji V_mod (p15).")


# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("PODSUMOWANIE  OP-3  --  p76_alpha_K_bifurcation.py  (v3)")
print("=" * 68)

print(f"\n  Parametr            | Obliczone   | Referencja  | Status")
print(f"  " + "-" * 58)
Kp1_s = f"{K1_check:.6f}" if np.isfinite(K1_check) else "---"
print(f"  K*_1 (cross-check)  | {Kp1_s:>11} | {K1_REF:>11.6f} | "
      f"{'OK (<5%)' if ok_P2 else 'ROZNICA'}")
if alpha_c_interp is not None:
    print(f"  alpha_c (K*_1=K1ref)| {alpha_c_interp:>11.4f} | {ALPHA_K_REF:>11.4f} | "
          f"{'OK (<25%)' if ok_P4 else 'ROZNICA >25%'}")
else:
    print(f"  alpha_c (K*_1=K1ref)| {'---':>11} | {ALPHA_K_REF:>11.4f} | nie znaleziony")
if results_C and exists_C.any():
    print(f"  a_c(8.5616)         | {a_c_est:>11.3f} | {A_GAM_TGP:>11.4f} | "
          f"{'~OK' if abs(a_c_est - A_GAM_TGP) < 0.005 else 'patrz wniosek'}")
print(f"  K*_2, K*_3          | {'N/A (p15)':>11} | {'wymagane':>11} | Path 1 otwarty")

print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

print("\n  [INTERPRETACJA OP-3]")
print("  Model: ODE bubble soliton, kryterium g=E/(4piK)-1=0 (poprawione)")
if ok_P1 and ok_P2:
    print(f"\n  P1/P2 PASS: K*_1 = {K1_check:.6f} odtworzone z dokladnoscia < 5%.")
    print(f"  Kryterium energetyczne dziala poprawnie.")
if ok_P3 and ok_P4:
    print(f"\n  P3/P4 PASS: alpha_c = {alpha_c_interp:.4f} (K*_1=K1_REF przy tym alpha)")
    print(f"  => Wartosc alpha_K_ref = {ALPHA_K_REF} jest WYZNACZALNA z substratu!")
    print(f"  => OP-3 Path 2 CZEŚCIOWO POTWIERDZONA (K*_1 = K1_REF jest kryterium).")
elif ok_P3:
    print(f"\n  P3 PASS, P4 FAIL: alpha_c = {alpha_c_interp:.4f} (odch > 25% od 8.5616).")
    print(f"  Model odtwarza alpha_c ale niezgodny z alpha_K_ref.")
else:
    print(f"\n  P3/P4 FAIL: K*_1(alpha) stale lub brak przeciecia z K1_REF.")
    print(f"  Mozliwe: K*_1 slabo zalezne od alpha w badanym zakresie.")
    print(f"  Wymagany szerszy skan lub inna definicja bifurkacji.")
print(f"\n  Path 1 (Koide): OTWARTA -- V_mod z 1 solitonem, K*_2 wymaga p15.")
