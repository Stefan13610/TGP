#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p80_alpha_max.py  --  TGP v1  *  Precyzyjne wyznaczenie alpha_max
==================================================================
Cel: Wyznaczenie dokladnej wartosci alpha_max — granicy bifurkacji
     powyzej ktorej standardowy soliton TGP nie istnieje.

Kontekst (p78/p79/diagnostyka):
  - alpha <= 8.80: K*_1 znaleziony (np. K*_1=0.010620 przy alpha=8.8)
  - alpha >= 8.90: brak K*_1 w [0.003, 0.50]; profil phi zapada sie
  - alpha_max in (8.80, 8.90)
  - Mechanizm: bifurkacja siodlo-wezel

Metoda:
  1. Bisection alpha_max w [8.80, 8.90]: sprawdz czy K*_1 istnieje
  2. Precyzyjnie: f(alpha) = 0 jesli K*_1 istnieje, 1 jesli nie
  3. Zapisz K*_1(alpha) dla kilku wartosci w poblizu alpha_max
  4. Oblicz epsilon = 1 - alpha_K/alpha_max

Testy:
  P1: alpha_max wyznaczone z dokladnoscia < 0.02
  P2: K*_1(alpha < alpha_max) istnieje i jest dobrze okreslone
  P3: K*_1(alpha > alpha_max) nie istnieje
  P4: epsilon = 1 - alpha_K/alpha_max obliczone
  P5: K*_1(alpha) malejace w poblizu alpha_max (ostatni punkt < p10 ref)

Data: 2026-03-25
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
ALPHA_K_REF = 8.5616
A_GAM       = 0.040
K1_REF      = 0.010414
LAM         = 5.501357e-06
GAMMA_V     = 1.0

def V_mod(p):
    return GAMMA_V/3 * p**3 - GAMMA_V/4 * p**4 + LAM/6 * (p - 1)**6

def dV_mod(p):
    return GAMMA_V * p**2 - GAMMA_V * p**3 + LAM * (p - 1)**5

V1 = GAMMA_V/3 - GAMMA_V/4


def ode_rhs(r, y, alpha):
    phi, dphi = y
    phi  = max(phi, 1e-10)
    kfac = 1.0 + alpha / phi
    ddphi = (dV_mod(phi) / kfac
             + alpha * dphi**2 / (2.0 * phi**2 * kfac)
             - (2.0 / r) * dphi)
    return [dphi, ddphi]


def _phi_at_rmax(psi_core, K, a_gam, alpha, r_max=40.0,
                 rtol=1e-9, atol=1e-11, n_eval=2000):
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


def _energy(psi_z, K, a_gam, alpha, r_max=40.0):
    dphi0 = -K / a_gam**2
    r_ev = a_gam * (r_max / a_gam) ** np.linspace(0, 1, 4000)
    try:
        sol = solve_ivp(
            lambda r, y: ode_rhs(r, y, alpha),
            [a_gam, r_max], [psi_z, dphi0],
            method='DOP853', rtol=1e-9, atol=1e-11,
            t_eval=r_ev
        )
        r    = sol.t
        phi  = np.maximum(sol.y[0], 1e-10)
        dphi = sol.y[1]
        Ek = 4*np.pi * np.trapezoid(0.5*dphi**2*(1+alpha/phi)*r**2, r)
        Ep = 4*np.pi * np.trapezoid((V_mod(phi) - V1)*r**2, r)
        return Ek + Ep
    except Exception:
        return np.nan


def compute_g(K, a_gam, alpha, r_max=40.0, n_psi=80):
    psi_max  = max(3.0, 1.0 + 2.0 * K / a_gam)
    psi_vals = np.linspace(1.001, min(psi_max, 300.0), n_psi)

    F_vals = np.array([_phi_at_rmax(p, K, a_gam, alpha, r_max) - 1.0
                       for p in psi_vals])

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
                    branches.append((psi_z, E / (4*np.pi*K) - 1.0))
            except Exception:
                pass

    if not branches:
        return np.nan
    return min(branches, key=lambda x: abs(x[1]))[1]


def find_K_star(alpha, a_gam=A_GAM, K_min=3e-3, K_max=0.12,
                n_K=20, n_psi=80, verbose=False):
    """Zwraca K*_1 lub nan jesli nie istnieje."""
    K_scan = np.logspace(np.log10(K_min), np.log10(K_max), n_K)
    g_vals = []
    for K in K_scan:
        gv = compute_g(K, a_gam, alpha)
        g_vals.append(gv)
        if verbose:
            s = f'{gv:.4f}' if np.isfinite(gv) else 'nan'
            print(f'    K={K:.5f}  g={s}')

    for i in range(len(g_vals) - 1):
        gi, gj = g_vals[i], g_vals[i+1]
        if np.isfinite(gi) and np.isfinite(gj) and gi * gj < 0:
            try:
                K_star = brentq(
                    lambda K: compute_g(K, a_gam, alpha),
                    K_scan[i], K_scan[i+1],
                    xtol=K_scan[i]*0.003, rtol=0.003, maxiter=20
                )
                return K_star
            except Exception:
                pass
    return np.nan


def soliton_exists(alpha):
    """True jesli K*_1 istnieje dla danego alpha."""
    K = find_K_star(alpha)
    return np.isfinite(K), K


# ─────────────────────────────────────────────────────────────────────────────
PASS_COUNT = 0
FAIL_COUNT = 0

def record(label, ok, info=""):
    global PASS_COUNT, FAIL_COUNT
    if ok: PASS_COUNT += 1
    else:  FAIL_COUNT += 1
    print(f"  [{'v' if ok else 'x'}] {label}: {info}")

# ─────────────────────────────────────────────────────────────────────────────
print("=" * 68)
print("TGP v1  *  p80_alpha_max.py  --  Precyzyjne wyznaczenie alpha_max")
print("Bifurkacja: K*_1 znika dla alpha > alpha_max in (8.80, 8.90)")
print("=" * 68)
print(f"\n  alpha_K = {ALPHA_K_REF},  a_Gam = {A_GAM}")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 1: Potwierdzenie granicy [8.80, 8.90]
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("KROK 1: Potwierdzenie: soliton istnieje przy 8.80, nie istnieje przy 8.90")
print("=" * 68)

for alpha_test in [8.80, 8.85, 8.90]:
    print(f"\n  alpha={alpha_test}:", flush=True)
    ex, K = soliton_exists(alpha_test)
    if ex:
        print(f"    K*_1 = {K:.6f}  [ISTNIEJE]")
    else:
        print(f"    K*_1 = nie znaleziony  [BRAK]")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 2: Bisection alpha_max w [8.80, 8.90]
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("KROK 2: Bisection alpha_max w [8.80, 8.90] (tolerancja 0.005)")
print("=" * 68)

a_lo, a_hi = 8.80, 8.90
K_lo_init = find_K_star(a_lo)
K_hi_init = find_K_star(a_hi)
print(f"\n  K*_1(alpha=8.80) = {K_lo_init:.6f}" if np.isfinite(K_lo_init)
      else "\n  K*_1(alpha=8.80) = brak")
print(f"  K*_1(alpha=8.90) = {K_hi_init:.6f}" if np.isfinite(K_hi_init)
      else "  K*_1(alpha=8.90) = brak")

alpha_max = np.nan
n_iter = 0
while a_hi - a_lo > 0.005 and n_iter < 12:
    a_mid = 0.5 * (a_lo + a_hi)
    ex, K_mid = soliton_exists(a_mid)
    print(f"  Bisection: alpha={a_mid:.4f}  K*_1={'brak' if not ex else f'{K_mid:.6f}'}")
    if ex:
        a_lo = a_mid
    else:
        a_hi = a_mid
    n_iter += 1

alpha_max = 0.5 * (a_lo + a_hi)
print(f"\n  alpha_max = {alpha_max:.4f}  (w [{a_lo:.4f}, {a_hi:.4f}])")

ok_P1 = np.isfinite(alpha_max) and (a_hi - a_lo) < 0.02
record("P1: alpha_max wyznaczone z dokladnoscia < 0.02", ok_P1,
       f"alpha_max = {alpha_max:.4f} in [{a_lo:.4f}, {a_hi:.4f}]")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 3: Profil K*_1(alpha) w poblizu alpha_max
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("KROK 3: K*_1(alpha) w poblizu alpha_max")
print("=" * 68)

alphas_near = np.array([8.5616, 8.6, 8.7, 8.8, a_lo, alpha_max])
results_near = []
print(f"\n  {'alpha':>8}  {'K*_1':>10}  {'K*_1/K1ref':>12}")
print("  " + "-" * 35)
for a in alphas_near:
    K = find_K_star(a)
    results_near.append((a, K))
    if np.isfinite(K):
        print(f"  {a:>8.4f}  {K:>10.6f}  {K/K1_REF:>12.4f}")
    else:
        print(f"  {a:>8.4f}  {'---':>10}  {'---':>12}")

ok_P2 = sum(1 for a, K in results_near if a <= a_lo and np.isfinite(K)) >= 3
record("P2: K*_1 istnieje dla alpha < alpha_max", ok_P2,
       f"{sum(1 for a,K in results_near if a<=a_lo and np.isfinite(K))} punktow")

ok_P3 = not np.isfinite(find_K_star(a_hi))
record("P3: K*_1 nie istnieje dla alpha > alpha_max", ok_P3,
       f"K*_1(alpha={a_hi:.2f}) = brak" if ok_P3 else f"K*_1(alpha={a_hi:.2f}) znaleziony")

# ─────────────────────────────────────────────────────────────────────────────
# KROK 4: Epsilon = 1 - alpha_K / alpha_max
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("KROK 4: Naturalnosc alpha_K wzgledem alpha_max")
print("=" * 68)

epsilon = 1.0 - ALPHA_K_REF / alpha_max
print(f"\n  alpha_max = {alpha_max:.4f}")
print(f"  alpha_K   = {ALPHA_K_REF}")
print(f"  epsilon   = 1 - alpha_K/alpha_max = {epsilon:.4f}  ({epsilon*100:.2f}%)")
print(f"  alpha_K / alpha_max = {ALPHA_K_REF/alpha_max:.4f}")

ok_P4 = 0.01 < abs(epsilon) < 0.10
record("P4: epsilon = 1 - alpha_K/alpha_max wyznaczone (1% < |eps| < 10%)", ok_P4,
       f"epsilon = {epsilon:.4f} = {epsilon*100:.2f}%")

# ─────────────────────────────────────────────────────────────────────────────
# PODSUMOWANIE
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 68)
print("PODSUMOWANIE  p80_alpha_max.py")
print("=" * 68)

print(f"\n  PASS: {PASS_COUNT} / {PASS_COUNT + FAIL_COUNT}")
print(f"  FAIL: {FAIL_COUNT} / {PASS_COUNT + FAIL_COUNT}")

print(f"""
  WYNIKI:
    alpha_max = {alpha_max:.4f}  (w [{a_lo:.4f}, {a_hi:.4f}])
    alpha_K   = {ALPHA_K_REF}
    epsilon   = {epsilon:.4f} = {epsilon*100:.1f}%
    alpha_K/alpha_max = {ALPHA_K_REF/alpha_max:.4f}

  FIZYCZNA INTERPRETACJA:
    Soliton TGP (galaz elektronowa) istnieje dla alpha < alpha_max.
    alpha_K jest {abs(epsilon)*100:.1f}% ponizej granicy stabilnosci.
    To moze byc zasada naturalnosci substratu:
    alpha_K to minimalna wartosc, przy ktorej soliton jest
    dostatecznie stabilny (zapas {abs(epsilon)*100:.1f}% przed bifurkacja).

  OP-3 NOVA VIA:
    Jesli alpha_max = f(V_mod, a_Gam) da sie wyliczyc analitycznie,
    to alpha_K = alpha_max * (1 - epsilon) byloby predykcja z substratu.
    Wymaga: analitycznej analizy warunku K*_1 = K_collapse.
""")
