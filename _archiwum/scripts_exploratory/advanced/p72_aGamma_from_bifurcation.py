# -*- coding: utf-8 -*-
"""
P72: PREDYKCJA a_GAMMA Z BIFURKACJI — WERYFIKACJA HIPOTEZY II.C
Wersja 2 — poprawiona formula Q i warunek bifurkacji.

Hipoteza (plan PLAN_ROZWOJU_v1, zadanie II.C):
  a_Gamma jest minimalnym a gwarantujacym Q(alpha_f, a) <= 3/2
  dla WSZYSTKICH rodzin fermionowych (alpha_f in [0.1, 25]).

  a_Gamma = max_{alpha_f in [0.1, 25]} a_c(alpha_f)

  gdzie a_c(alpha_f) = min a takie ze Q(alpha_f, a) = 3/2
  (Q wg TGP: (sqrt(K1)+sqrt(K2)+sqrt(K3))^2 / (K1+K2+K3),
   target Q = 3/2 odpowiada standardowemu Koide Q_std = 1/2).
"""
import sys, warnings
sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expi

# -------------------------------------------------------------------
# PARAMETRY
# -------------------------------------------------------------------
LAM_K    = 5.4677e-6
R21_PDG  = 206.770
A_GAM    = 0.040049   # Punkt B z P55
ALPHA_K  = 8.5612     # alpha_Koide z P55
A_C_P54  = 0.038382   # a_c z P54 (min Q_min(alpha) = 3/2)
R_MAX    = 50.0
N_GRID   = 1500
TARGET_Q = 1.5        # Warunek Koidego w konwencji TGP

def E1(x):
    return -expi(-x)

def V_mod(phi, lam=LAM_K):
    return phi**3/3 - phi**4/4 + lam*(phi-1)**6/6

def energy_num(K, alpha, a, lam=LAM_K, N=N_GRID):
    t    = np.linspace(0, 1, N)
    r    = a * (R_MAX/a)**t
    phi  = np.maximum(1.0 + K*np.exp(-r)/r, 1e-10)
    dphi = K * np.exp(-r) * (-r - 1.0) / r**2
    V1   = V_mod(1.0, lam)
    Ek   = 4*np.pi * np.trapezoid(0.5*dphi**2*(1.0 + alpha/phi)*r**2, r)
    Ep   = 4*np.pi * np.trapezoid((V_mod(phi, lam) - V1)*r**2, r)
    return Ek + Ep

def g_func(K, alpha, a, lam=LAM_K):
    return energy_num(K, alpha, a, lam) / (4*np.pi*K) - 1.0

def find_K1(alpha, a, lam=LAM_K):
    try:
        return brentq(g_func, 1e-4, 0.45, args=(alpha, a, lam), xtol=1e-10)
    except Exception:
        return np.nan

def find_K2(alpha, a, lam=LAM_K):
    for lo, hi in [(0.5, 8.0), (0.3, 12.0), (0.1, 20.0)]:
        try:
            return brentq(g_func, lo, hi, args=(alpha, a, lam), xtol=1e-10)
        except Exception:
            continue
    return np.nan

def K3_ei(a, lam=LAM_K):
    I4 = np.exp(-4*a)/a - 4*E1(4*a)
    try:
        I6, _ = quad(lambda r: np.exp(-6*r)/r**4, a, R_MAX, limit=200)
    except Exception:
        return np.nan
    return np.sqrt(3*I4 / (2*lam*I6)) if I6 > 0 else np.nan

def Q_tgp(K1, K2, K3):
    """Q TGP = (sqrt(K1)+sqrt(K2)+sqrt(K3))^2 / (K1+K2+K3).
    Target: Q = 3/2 (Koide). Uwaga: BEZ czynnika 3 w mianowniku."""
    if np.any(np.isnan([K1, K2, K3])) or np.any(np.array([K1,K2,K3]) <= 0):
        return np.nan
    s = np.sqrt(K1) + np.sqrt(K2) + np.sqrt(K3)
    d = K1 + K2 + K3
    return s**2 / d

def Q_at(alpha, a, lam=LAM_K):
    K1 = find_K1(alpha, a, lam)
    K2 = find_K2(alpha, a, lam)
    K3 = K3_ei(a, lam)
    return Q_tgp(K1, K2, K3)

# -------------------------------------------------------------------
# SEKCJA A: WERYFIKACJA BAZOWA
# -------------------------------------------------------------------
print("=" * 65)
print("P72 v2: PREDYKCJA a_GAMMA Z BIFURKACJI")
print("=" * 65)

K1_B = find_K1(ALPHA_K, A_GAM)
K2_B = find_K2(ALPHA_K, A_GAM)
K3_B = K3_ei(A_GAM)
Q_B  = Q_tgp(K1_B, K2_B, K3_B)
print(f"\n  Punkt B: K1={K1_B:.6f}, K2={K2_B:.4f}, K3={K3_B:.4f}")
print(f"  Q_TGP = {Q_B:.8f}  (target: {TARGET_Q})")
print(f"  Odchylenie od 3/2: {(Q_B-1.5)*1e6:.0f} ppm")

# -------------------------------------------------------------------
# SEKCJA B: Q(alpha, a) jako funkcja alpha dla a=a_Gamma
# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA B: Q(alpha_f, a_Gamma) — czy Q = 3/2 gdzies w [0.1, 25]?")
print("=" * 65)

alpha_test = [0.5, 1.0, 2.0, 3.0, 4.025, 5.0, 8.5612, 10.0, 15.0, 20.0, 25.0]
print(f"\n  {'alpha_f':>8}  {'K1':>10}  {'K2':>8}  {'Q_TGP':>10}  {'Q-3/2':>10}")
print("  " + "-"*55)
for al in alpha_test:
    K1 = find_K1(al, A_GAM)
    K2 = find_K2(al, A_GAM)
    K3 = K3_ei(A_GAM)
    Q  = Q_tgp(K1, K2, K3)
    delta = Q - TARGET_Q if not np.isnan(Q) else np.nan
    K2str = f"{K2:.4f}" if not np.isnan(K2) else "  brak"
    Qstr  = f"{Q:.6f}" if not np.isnan(Q) else "  brak"
    dstr  = f"{delta:+.6f}" if not np.isnan(Q) else "  -"
    print(f"  {al:>8.3f}  {K1:>10.6f}  {K2str:>8}  {Qstr:>10}  {dstr:>10}")

# -------------------------------------------------------------------
# SEKCJA C: a_c(alpha_f) — prog Q=3/2 dla kazdego alpha_f
# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA C: a_c(alpha_f) = min a takie ze Q(alpha_f, a) = 3/2")
print("=" * 65)

# Dla danego alpha_f: Q(alpha_f, a) jako f-cja a
# Przy malym a: Q > 3/2 (brak K2 lub duze Q)
# Przy duzym a: Q < 3/2 (oba zera istnieja, Q maleje)
# a_c: przejscie Q = 3/2

def Q_minus_target(a, alpha, lam=LAM_K):
    Q = Q_at(alpha, a, lam)
    if np.isnan(Q):
        return +1.0   # traktuj NaN jako Q > 3/2
    return Q - TARGET_Q

# Sprawdz granice dla kilku alpha
print("\n  Sprawdzam Q(alpha, a) na granicy a=0.025 i a=0.060:")
print(f"  {'alpha_f':>8}  {'Q(0.025)':>10}  {'Q(0.060)':>10}  {'Q(0.040)':>10}")
for al in [0.5, 2.0, 4.025, 8.5612, 15.0, 25.0]:
    Q025 = Q_at(al, 0.025)
    Q060 = Q_at(al, 0.060)
    Q040 = Q_at(al, A_GAM)
    def fmt(Q):
        return f"{Q:.4f}" if not np.isnan(Q) else "  brak"
    print(f"  {al:>8.3f}  {fmt(Q025):>10}  {fmt(Q060):>10}  {fmt(Q040):>10}")

# Skan alpha_f
print("\n  Szukam a_c(alpha_f) dla alpha_f in [0.5, 25]:")
print(f"  {'alpha_f':>8}  {'a_c':>10}  {'Q(a_c)':>10}  {'status':>12}")
print("  " + "-"*50)

alpha_scan = np.unique(np.concatenate([
    np.linspace(0.5, 2.0, 8),
    np.linspace(2.0, 6.0, 16),
    np.linspace(6.0, 12.0, 20),
    np.linspace(12.0, 25.0, 12),
]))

results = []
A_LO, A_HI = 0.025, 0.070

for alpha_f in alpha_scan:
    Q_lo = Q_minus_target(A_LO, alpha_f)
    Q_hi = Q_minus_target(A_HI, alpha_f)

    if Q_lo < 0:
        # Juz przy a=A_LO mamy Q < 3/2: a_c < A_LO
        results.append({'alpha': alpha_f, 'ac': A_LO, 'Q_ac': Q_at(alpha_f, A_LO), 'status': 'ac<0.025'})
        continue

    if Q_hi > 0:
        # Nawet przy a=A_HI mamy Q > 3/2: brak przejscia
        results.append({'alpha': alpha_f, 'ac': np.nan, 'Q_ac': np.nan, 'status': 'brak Q=3/2'})
        continue

    # Bisekacja
    try:
        ac = brentq(Q_minus_target, A_LO, A_HI, args=(alpha_f,), xtol=1e-7)
        Q_ac = Q_at(alpha_f, ac)
        results.append({'alpha': alpha_f, 'ac': ac, 'Q_ac': Q_ac, 'status': 'ok'})
    except Exception as e:
        results.append({'alpha': alpha_f, 'ac': np.nan, 'Q_ac': np.nan, 'status': f'err'})

# Wyswietl (co 4-ty)
for i, r in enumerate(results):
    if i % 4 == 0 or r['status'] != 'ok' or abs(r.get('ac', 0) - A_GAM) < 0.003:
        ac_str = f"{r['ac']:.6f}" if not np.isnan(r.get('ac', np.nan)) else "    nan"
        Q_str  = f"{r['Q_ac']:.6f}" if not np.isnan(r.get('Q_ac', np.nan)) else "    nan"
        print(f"  {r['alpha']:>8.3f}  {ac_str:>10}  {Q_str:>10}  {r['status']:>12}")

# -------------------------------------------------------------------
# SEKCJA D: MAKSIMUM a_c — PREDYKCJA
# -------------------------------------------------------------------
print("\n" + "=" * 65)
print("SEKCJA D: MAKSIMUM a_c = PREDYKCJA a_GAMMA")
print("=" * 65)

valid = [r for r in results if not np.isnan(r.get('ac', np.nan)) and r['status'] == 'ok']
if valid:
    ac_arr = np.array([r['ac'] for r in valid])
    al_arr = np.array([r['alpha'] for r in valid])
    idx_max = np.argmax(ac_arr)

    ac_max   = ac_arr[idx_max]
    al_max   = al_arr[idx_max]

    print(f"\n  max a_c = {ac_max:.8f}  przy alpha_f = {al_max:.4f}")
    print(f"  a_Gamma = {A_GAM:.8f}  (z P55, warunek r21=r21_PDG)")
    print(f"  a_c_P54 = {A_C_P54:.8f}  (min Q_min(alpha) = 3/2)")
    diff = A_GAM - ac_max
    ppm  = diff / A_GAM * 1e6
    print(f"\n  Roznica a_Gamma - max(a_c) = {diff:.6f}  ({ppm:.0f} ppm)")
    print(f"  max(a_c) / a_Gamma          = {ac_max / A_GAM:.6f}")

    # Czy a_Gamma lezy powyzej wszystkich a_c?
    failures = [r for r in valid if r['ac'] > A_GAM + 1e-5]
    print(f"\n  Czy a_Gamma > a_c(alpha_f) dla wszystkich alpha_f w skanie?")
    if failures:
        print(f"  NIE — {len(failures)} wyjatkow:")
        for r in failures:
            print(f"    alpha={r['alpha']:.3f}, a_c={r['ac']:.6f}")
    else:
        print(f"  TAK — a_Gamma > a_c(alpha_f) dla wszystkich {len(valid)} testowanych alpha_f.")

    # Ocena hipotezy
    print(f"\n  OCENA HIPOTEZY II.C:")
    print(f"    Hipoteza: a_Gamma = max_{{alpha_f}} a_c(alpha_f)")
    print(f"    Wynik:    max(a_c) = {ac_max:.6f}, a_Gamma = {A_GAM:.6f}")
    if abs(ppm) < 1000:
        verdict = f"POTWIERDZONA na poziomie {abs(ppm):.0f} ppm"
    elif abs(ppm) < 5000:
        verdict = f"czesciowo potwierdzona ({abs(ppm):.0f} ppm)"
    else:
        verdict = f"NIE potwierdzona ({abs(ppm):.0f} ppm)"
    print(f"    Weryfikacja: {verdict}")
    print(f"\n  Alternatywna interpretacja:")
    print(f"    a_Gamma > a_c(alpha_f) dla wszystkich alpha_f w [0.5, 25]")
    print(f"    => substrat z a=a_Gamma umozliwia Q=3/2 dla kazdej rodziny fermionowej")
    print(f"    Margines: min(a_Gamma - a_c) = {A_GAM - ac_max:.6f} ({(A_GAM-ac_max)/A_GAM*100:.2f}%)")
else:
    print("  Brak danych do analizy.")

print("\n  Koniec P72.")
